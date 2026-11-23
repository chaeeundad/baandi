//! G-FILTER-REPLAY: fixed-potential, resident-block Chebyshev experiment.
//!
//! No SCF caller or environment switch. The caller validates snapshot identity,
//! supplies the captured stopping profile and executes Davidson if `fallback`
//! is set. Work includes failed sweeps; reference solves belong to the harness.
//! All energies and residuals are Hartree. ORC rows are unchanged; applications
//! use the same distributed operator interface (including paired Gamma FFTs).
use crate::{
    SolverError,
    distributed_davidson::{DistributedRealOperator, PanelGs, panel_column},
};
use dft_mpi::Collective;
use dft_numerics::{
    RealGemmWorkspace, RealMatrix, linalg::DEPENDENCE_TOLERANCE, real::RealEighWorkspace,
};

#[derive(Clone, Copy, Debug)]
pub struct Profile {
    pub tolerance: f64,
    pub strict_bands: usize,
    pub buffer_tolerance: f64,
    pub ritz: bool,
    pub occupied: usize,
}
impl Profile {
    fn threshold(self, j: usize) -> f64 {
        if j < self.strict_bands - usize::from(self.ritz) {
            self.tolerance
        } else {
            self.buffer_tolerance
        }
    }
}

/// Optional replay instrumentation. No clocks, I/O or logging live in the solver.
#[derive(Clone, Copy, Debug)]
#[repr(usize)]
pub enum Phase {
    Setup,
    Bounds,
    Filter,
    Qr,
    TerminalApply,
    RayleighRitz,
}
struct PhaseSpan {
    observer: Option<fn(Phase, bool)>,
    phase: Phase,
}
impl PhaseSpan {
    fn new(observer: Option<fn(Phase, bool)>, phase: Phase) -> Self {
        if let Some(f) = observer {
            f(phase, true);
        }
        Self { observer, phase }
    }
}
impl Drop for PhaseSpan {
    fn drop(&mut self) {
        if let Some(f) = self.observer {
            f(self.phase, false);
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Options {
    /// None derives the degree from the current spectrum and requested error.
    /// Explicit degrees are experiment arms, never material-specific defaults.
    pub degree: Option<usize>,
    pub max_sweeps: usize,
    pub bound_steps: usize,
    pub observer: Option<fn(Phase, bool)>,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            degree: None,
            max_sweeps: 24,
            bound_steps: 12,
            observer: None,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct Work {
    /// Indexed by exact width: number of calls, not number of columns.
    pub h_calls_by_width: Vec<u64>,
    pub bound_h_columns: u64,
    pub terminal_h_columns: u64,
    pub qr_calls: u64,
    pub qr_flops: u64,
    pub gram_flops: u64,
    pub rotate_flops: u64,
    pub recurrence_flops: u64,
    /// Actual reduced eigenproblem orders, including the initial RR.
    pub eigen_orders: Vec<usize>,
    pub reduction_calls: u64,
    pub reduction_doubles: u64,
    /// Explicit block copies (read + write bytes); excludes kernel internals.
    pub copy_bytes: u64,
}
impl Work {
    pub fn h_columns(&self) -> u64 {
        self.h_calls_by_width
            .iter()
            .enumerate()
            .map(|(w, &n)| w as u64 * n)
            .sum()
    }
    fn reduction(&mut self, n: usize) {
        self.reduction_calls += 1;
        self.reduction_doubles += n as u64;
    }
}

#[derive(Clone, Debug)]
pub struct Sweep {
    pub degree: usize,
    pub max_residual_ha: f64,
    pub max_occupied_change_ha: f64,
    pub gap_ha: f64,
    pub converged: bool,
    pub h_columns: u64,
    pub failure: Option<&'static str>,
}
#[derive(Debug)]
pub struct Attempt {
    pub vectors: RealMatrix,
    pub eigenvalues: Vec<f64>,
    pub residuals: Vec<f64>,
    pub sweeps: Vec<Sweep>,
    pub work: Work,
    pub lower_ha: f64,
    pub upper_estimate_ha: f64,
    pub cutoff_ha: f64,
    pub fallback: Option<&'static str>,
}

fn mpi(e: dft_mpi::MpiError) -> SolverError {
    SolverError::Application(e.to_string())
}
fn sum<C: Collective>(c: &C, a: &[f64], b: &mut [f64], w: &mut Work) -> Result<(), SolverError> {
    w.reduction(a.len());
    c.tree_sum_f64(a, b).map_err(mpi)
}
fn apply<C: Collective, O: DistributedRealOperator>(
    c: &C,
    op: &mut O,
    x: &RealMatrix,
    y: &mut RealMatrix,
    w: &mut Work,
) -> Result<(), SolverError> {
    w.h_calls_by_width[x.cols()] += 1;
    op.apply_block(c, x, y)
}

struct Scratch {
    gemm: RealGemmWorkspace,
    eigen: RealEighWorkspace,
    partial: RealMatrix,
    gram: RealMatrix,
    coeff: RealMatrix,
    rotated: RealMatrix,
    local: Vec<f64>,
    sums: Vec<f64>,
}
impl Scratch {
    fn new<O: DistributedRealOperator>(op: &O, b: usize) -> Result<Self, SolverError> {
        let n = op.local_rows().len();
        Ok(Self {
            gemm: RealGemmWorkspace::with_pool(
                n.max(b),
                b,
                op.worker_pool().unwrap_or_else(|| {
                    dft_numerics::pool::WorkerPool::new(op.inner_threads()).expect("worker pool")
                }),
            )?,
            eigen: RealEighWorkspace::new(b)?,
            partial: RealMatrix::zeros(b, b)?,
            gram: RealMatrix::zeros(b, b)?,
            coeff: RealMatrix::zeros(b, b)?,
            rotated: RealMatrix::zeros(n, b)?,
            local: Vec::with_capacity(b + 1),
            sums: Vec::with_capacity(b + 1),
        })
    }
    fn gram<C: Collective>(
        &mut self,
        c: &C,
        x: &RealMatrix,
        y: &RealMatrix,
        w: &mut Work,
    ) -> Result<(), SolverError> {
        crate::distributed::gram_slices(
            c,
            &mut self.gemm,
            x.as_slice(),
            x.cols(),
            y.as_slice(),
            y.cols(),
            x.rows(),
            &mut self.partial,
            &mut self.gram,
        )?;
        w.reduction(x.cols() * y.cols());
        w.gram_flops += (2 * x.rows() * x.cols() * y.cols()) as u64;
        Ok(())
    }
    /// Existing two-pass panel QR, with shared rank threshold. Each input
    /// column is normalized first so the rank decision is scale-independent.
    fn qr<C: Collective>(
        &mut self,
        c: &C,
        x: &mut RealMatrix,
        w: &mut Work,
    ) -> Result<bool, SolverError> {
        w.qr_calls += 1;
        let (n, b) = (x.rows(), x.cols());
        for j in 0..b {
            let (q, tail) = x.as_mut_slice().split_at_mut(j * n);
            let v = &mut tail[..n];
            let a = [v.iter().map(|v| v * v).sum()];
            let mut z = [0.0];
            sum(c, &a, &mut z, w)?;
            let before = z[0].sqrt();
            if !before.is_finite() || before == 0.0 {
                return Ok(false);
            }
            for v in v.iter_mut() {
                *v /= before;
            }
            let len = panel_column(
                c,
                q,
                n,
                j,
                v,
                &mut self.local,
                &mut self.sums,
                PanelGs::Classic,
                None,
            )?;
            if j > 0 {
                w.reduction(j);
                w.reduction(j);
            }
            w.reduction(1);
            w.qr_flops += (8 * n * j + 6 * n) as u64;
            if !len.is_finite() || len <= DEPENDENCE_TOLERANCE {
                return Ok(false);
            }
            for v in v {
                *v /= len;
            }
        }
        Ok(true)
    }
    fn rr<C: Collective>(
        &mut self,
        c: &C,
        x: &mut RealMatrix,
        hx: &mut RealMatrix,
        values: &mut [f64],
        w: &mut Work,
    ) -> Result<(), SolverError> {
        let (n, b) = (x.rows(), x.cols());
        self.gram(c, x, hx, w)?;
        for j in 0..b {
            for i in 0..j {
                let a = 0.5 * (self.gram[(i, j)] + self.gram[(j, i)]);
                self.gram[(i, j)] = a;
                self.gram[(j, i)] = a;
            }
        }
        w.eigen_orders.push(b);
        self.eigen.solve_into(&self.gram, values, &mut self.coeff)?;
        for matrix in [x, hx] {
            crate::distributed::rotate_local_with(
                &mut self.gemm,
                matrix,
                &self.coeff,
                &mut self.rotated,
                true,
            )?;
            std::mem::swap(matrix, &mut self.rotated);
            w.rotate_flops += (2 * n * b * b) as u64;
        }
        Ok(())
    }
}

fn residual<C: Collective>(
    c: &C,
    x: &RealMatrix,
    hx: &RealMatrix,
    values: &[f64],
    out: &mut [f64],
    local: &mut [f64],
    w: &mut Work,
) -> Result<(), SolverError> {
    let n = x.rows();
    for j in 0..x.cols() {
        local[j] = (0..n)
            .map(|i| (hx[(i, j)] - values[j] * x[(i, j)]).powi(2))
            .sum();
    }
    sum(c, local, out, w)?;
    for r in out {
        *r = r.sqrt();
    }
    Ok(())
}

/// Heuristic upper estimate from shifted block power iteration (two columns
/// preserve Gamma pairing). The 20% padding is dimensionless, not a guarantee
/// of spectral completeness. The replay compares to an independent strict
/// solve; failure or an unresolved edge must never become an SCF result.
fn upper<C: Collective, O: DistributedRealOperator>(
    c: &C,
    op: &mut O,
    lower: f64,
    steps: usize,
    w: &mut Work,
) -> Result<f64, SolverError> {
    let n = op.local_rows().len();
    let mut x = RealMatrix::zeros(n, 2)?;
    let mut y = RealMatrix::zeros(n, 2)?;
    for j in 0..2 {
        for (i, &g) in op.local_rows().iter().enumerate() {
            let mut z = (g as u64 + 1).wrapping_add((j as u64 + 1) * 0x9e3779b9);
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
            x[(i, j)] = ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64 - 0.5;
        }
    }
    let mut bound = lower;
    for k in 0..=steps {
        let a = [
            (0..n).map(|i| x[(i, 0)].powi(2)).sum(),
            (0..n).map(|i| x[(i, 1)].powi(2)).sum(),
        ];
        let mut norms = [0.0; 2];
        sum(c, &a, &mut norms, w)?;
        if norms.iter().any(|a| !a.is_finite() || *a <= 0.0) {
            return Ok(f64::NAN);
        }
        if k > 0 {
            bound = lower + 1.2 * norms.iter().map(|a| a.sqrt()).fold(0.0, f64::max);
        }
        if k == steps {
            break;
        }
        for j in 0..2 {
            for i in 0..n {
                x[(i, j)] /= norms[j].sqrt();
            }
        }
        apply(c, op, &x, &mut y, w)?;
        w.bound_h_columns += 2;
        for (a, b) in y.as_mut_slice().iter_mut().zip(x.as_slice()) {
            *a -= lower * b;
        }
        std::mem::swap(&mut x, &mut y);
    }
    Ok(bound)
}

/// Degree from attenuation at the occupied edge; no material-dependent table.
/// A degree above 64 asks for fallback rather than silently weakening accuracy.
pub fn degree_from_spectrum(
    lower: f64,
    upper: f64,
    cutoff: f64,
    edge: f64,
    error: f64,
    tol: f64,
) -> Option<usize> {
    let radius = (upper - cutoff) * 0.5;
    let z = ((upper + cutoff) * 0.5 - edge) / radius;
    if !(lower < edge && edge < cutoff && cutoff < upper && tol > 0.0 && error.is_finite()) {
        return None;
    }
    let d = ((2.0 * (error / tol).max(1.0)).acosh() / z.acosh())
        .ceil()
        .max(2.0) as usize;
    (d <= 64).then_some(d.next_multiple_of(2))
}

/// Attempt the filter. A failed attempt returns its last finite diagnostics and
/// all charged work. The caller MUST run the original Davidson on `fallback`.
/// Arrays and journal capacity are reserved before the sweep loop.
pub fn solve<C: Collective, O: DistributedRealOperator>(
    c: &C,
    op: &mut O,
    start: &RealMatrix,
    kept_hx: Option<&RealMatrix>,
    profile: Profile,
    options: Options,
) -> Result<Attempt, SolverError> {
    let setup_span = PhaseSpan::new(options.observer, Phase::Setup);
    let (n, b) = (start.rows(), start.cols());
    if n != op.local_rows().len()
        || b < 2
        || profile.occupied == 0
        || profile.occupied >= b
        || profile.strict_bands == 0
        || profile.strict_bands > b
        || profile.tolerance <= 0.0
        || !profile.tolerance.is_finite()
        || !profile.buffer_tolerance.is_finite()
        || profile.buffer_tolerance < profile.tolerance
        || options.max_sweeps == 0
        || options.bound_steps == 0
        || options.degree.is_some_and(|d| !(2..=64).contains(&d))
    {
        return Err(SolverError::InvalidParameter(
            "filter dimensions/profile/options",
        ));
    }
    let mut a = Attempt {
        vectors: start.clone(),
        eigenvalues: vec![0.0; b],
        residuals: vec![0.0; b],
        sweeps: Vec::with_capacity(options.max_sweeps),
        work: Work {
            h_calls_by_width: vec![0; b.max(2) + 1],
            eigen_orders: Vec::with_capacity(options.max_sweeps + 1),
            copy_bytes: (16 * n * b) as u64,
            ..Work::default()
        },
        lower_ha: f64::NAN,
        upper_estimate_ha: f64::NAN,
        cutoff_ha: f64::NAN,
        fallback: None,
    };
    let mut s = Scratch::new(op, b)?;
    let mut hx = RealMatrix::zeros(n, b)?;
    s.gram(c, &a.vectors, &a.vectors, &mut a.work)?;
    let defect = (0..b)
        .flat_map(|j| (0..b).map(move |i| (i, j)))
        .map(|(i, j)| (s.gram[(i, j)] - f64::from(u8::from(i == j))).abs())
        .fold(0.0, f64::max);
    if !defect.is_finite() || defect > 1e-12 {
        if !s.qr(c, &mut a.vectors, &mut a.work)? {
            a.fallback = Some("initial_rank_loss");
            return Ok(a);
        }
        if apply(c, op, &a.vectors, &mut hx, &mut a.work).is_err() {
            a.fallback = Some("initial_operator_failure");
            return Ok(a);
        }
    } else if let Some(h) = kept_hx {
        if h.rows() != n || h.cols() != b {
            return Err(SolverError::InvalidParameter("kept HX shape"));
        }
        hx.as_mut_slice().copy_from_slice(h.as_slice());
        a.work.copy_bytes += (16 * n * b) as u64;
    } else {
        if apply(c, op, &a.vectors, &mut hx, &mut a.work).is_err() {
            a.fallback = Some("initial_operator_failure");
            return Ok(a);
        }
    }
    let mut local = vec![0.0; b];
    if s.rr(c, &mut a.vectors, &mut hx, &mut a.eigenvalues, &mut a.work)
        .is_err()
    {
        a.fallback = Some("initial_nonfinite_or_eigen_failure");
        return Ok(a);
    }
    residual(
        c,
        &a.vectors,
        &hx,
        &a.eigenvalues,
        &mut a.residuals,
        &mut local,
        &mut a.work,
    )?;
    a.lower_ha = a.eigenvalues[0] - a.residuals.iter().copied().fold(0.0, f64::max);
    drop(setup_span);
    let bound_span = PhaseSpan::new(options.observer, Phase::Bounds);
    a.upper_estimate_ha = match upper(c, op, a.lower_ha, options.bound_steps, &mut a.work) {
        Ok(bound) => bound,
        Err(_) => {
            a.fallback = Some("bound_operator_failure");
            return Ok(a);
        }
    };
    drop(bound_span);
    a.cutoff_ha = a.eigenvalues[b - 1]
        + a.residuals[b - 1].max((a.eigenvalues[b - 1] - a.eigenvalues[0]) / (b - 1) as f64);
    if !a.upper_estimate_ha.is_finite()
        || a.upper_estimate_ha <= a.cutoff_ha
        || a.lower_ha >= a.cutoff_ha
    {
        a.fallback = Some("invalid_spectral_interval");
        return Ok(a);
    }
    let mut prev = RealMatrix::zeros(n, b)?;
    let mut curr = RealMatrix::zeros(n, b)?;
    let mut previous = vec![0.0; b];
    let mut best = a.residuals[..].iter().copied().fold(0.0, f64::max);
    let mut stalled = 0;
    for _ in 0..options.max_sweeps {
        let degree = options.degree.or_else(|| {
            degree_from_spectrum(
                a.lower_ha,
                a.upper_estimate_ha,
                a.cutoff_ha,
                a.eigenvalues[profile.occupied - 1],
                best,
                profile.tolerance,
            )
        });
        let Some(degree) = degree else {
            a.fallback = Some("degree_budget_or_unresolved_edge");
            break;
        };
        let mut row = Sweep {
            degree,
            max_residual_ha: f64::NAN,
            max_occupied_change_ha: f64::NAN,
            gap_ha: f64::NAN,
            converged: false,
            h_columns: 0,
            failure: None,
        };
        let before = a.work.h_columns();
        let filter_span = PhaseSpan::new(options.observer, Phase::Filter);
        previous.copy_from_slice(&a.eigenvalues);
        prev.as_mut_slice().copy_from_slice(a.vectors.as_slice());
        a.work.copy_bytes += (16 * n * b + 16 * b) as u64;
        let e = (a.upper_estimate_ha - a.cutoff_ha) * 0.5;
        let center = (a.upper_estimate_ha + a.cutoff_ha) * 0.5;
        let alpha = (a.lower_ha - center) / e;
        let mut sigma = 1.0 / alpha;
        for i in 0..n * b {
            curr.as_mut_slice()[i] = sigma / e * (hx.as_slice()[i] - center * prev.as_slice()[i]);
        }
        for _ in 1..degree {
            if apply(c, op, &curr, &mut hx, &mut a.work).is_err() {
                row.failure = Some("operator_failure");
                break;
            }
            let next = 1.0 / (2.0 * alpha - sigma);
            for i in 0..n * b {
                prev.as_mut_slice()[i] = 2.0 * next / e
                    * (hx.as_slice()[i] - center * curr.as_slice()[i])
                    - next * sigma * prev.as_slice()[i];
            }
            sigma = next;
            std::mem::swap(&mut prev, &mut curr);
        }
        a.work.recurrence_flops += ((3 + 6 * (degree - 1)) * n * b) as u64;
        drop(filter_span);
        if row.failure.is_none() {
            std::mem::swap(&mut a.vectors, &mut curr);
            let qr_span = PhaseSpan::new(options.observer, Phase::Qr);
            let rank_ok = s.qr(c, &mut a.vectors, &mut a.work)?;
            drop(qr_span);
            if !rank_ok {
                row.failure = Some("rank_loss_or_nonfinite");
            } else {
                a.work.terminal_h_columns += b as u64;
                let terminal_span = PhaseSpan::new(options.observer, Phase::TerminalApply);
                let applied = apply(c, op, &a.vectors, &mut hx, &mut a.work);
                drop(terminal_span);
                let rr_span = PhaseSpan::new(options.observer, Phase::RayleighRitz);
                let rr_failed = applied.is_err()
                    || s.rr(c, &mut a.vectors, &mut hx, &mut a.eigenvalues, &mut a.work)
                        .is_err();
                drop(rr_span);
                if rr_failed {
                    row.failure = Some("terminal_operator_or_eigen_failure");
                } else {
                    residual(
                        c,
                        &a.vectors,
                        &hx,
                        &a.eigenvalues,
                        &mut a.residuals,
                        &mut local,
                        &mut a.work,
                    )?;
                    row.max_residual_ha = a.residuals[..profile.occupied]
                        .iter()
                        .copied()
                        .fold(0.0, f64::max);
                    row.max_occupied_change_ha = (0..profile.occupied)
                        .map(|j| (a.eigenvalues[j] - previous[j]).abs())
                        .fold(0.0, f64::max);
                    row.gap_ha =
                        a.eigenvalues[profile.occupied] - a.eigenvalues[profile.occupied - 1];
                    row.converged = (0..b).all(|j| {
                        if profile.ritz {
                            (a.eigenvalues[j] - previous[j]).abs() < profile.threshold(j)
                        } else {
                            a.residuals[j] <= profile.threshold(j)
                        }
                    });
                    if !row.max_residual_ha.is_finite() {
                        row.failure = Some("nonfinite_residual");
                    } else if row.converged
                        && row.gap_ha
                            <= a.residuals[profile.occupied - 1] + a.residuals[profile.occupied]
                    {
                        row.failure = Some("unresolved_occupied_edge");
                        row.converged = false;
                    } else if a.residuals.iter().copied().fold(0.0, f64::max) >= best {
                        stalled += 1;
                    } else {
                        best = a.residuals.iter().copied().fold(0.0, f64::max);
                        stalled = 0;
                    }
                    if stalled >= 3 && !row.converged {
                        row.failure = Some("stalled_progress");
                    }
                }
            }
        }
        if row.converged && !profile.ritz && row.failure.is_none() {
            let _certificate_span = PhaseSpan::new(options.observer, Phase::TerminalApply);
            a.work.terminal_h_columns += b as u64;
            if apply(c, op, &a.vectors, &mut hx, &mut a.work).is_err() {
                row.failure = Some("certificate_operator_failure");
                row.converged = false;
            } else {
                residual(
                    c,
                    &a.vectors,
                    &hx,
                    &a.eigenvalues,
                    &mut a.residuals,
                    &mut local,
                    &mut a.work,
                )?;
                row.max_residual_ha = a.residuals[..profile.occupied]
                    .iter()
                    .copied()
                    .fold(0.0, f64::max);
                row.converged = (0..b).all(|j| a.residuals[j] <= profile.threshold(j));
            }
        }
        row.h_columns = a.work.h_columns() - before;
        let (done, failure) = (row.converged, row.failure);
        a.sweeps.push(row);
        if let Some(f) = failure {
            a.fallback = Some(f);
            break;
        }
        if done {
            return Ok(a);
        }
    }
    if a.fallback.is_none() {
        a.fallback = Some("sweep_budget");
    }
    Ok(a)
}

#[cfg(test)]
mod tests {
    use super::*;
    struct Diagonal {
        rows: Vec<usize>,
        values: Vec<f64>,
    }
    impl DistributedRealOperator for Diagonal {
        fn global_len(&self) -> usize {
            self.values.len()
        }
        fn local_rows(&self) -> &[usize] {
            &self.rows
        }
        fn apply_block<C: Collective>(
            &mut self,
            _: &C,
            x: &RealMatrix,
            y: &mut RealMatrix,
        ) -> Result<(), SolverError> {
            for j in 0..x.cols() {
                for (i, &g) in self.rows.iter().enumerate() {
                    y[(i, j)] = self.values[g] * x[(i, j)];
                }
            }
            Ok(())
        }
    }
    fn profile(ritz: bool) -> Profile {
        Profile {
            tolerance: 1e-9,
            strict_bands: 4,
            buffer_tolerance: 1e-8,
            ritz,
            occupied: 3,
        }
    }
    fn run(ranks: usize, ritz: bool, degree: usize) {
        let outputs = dft_mpi::simulate(ranks, |c| {
            let rows: Vec<_> = (0..24).filter(|i| i % ranks == c.rank()).collect();
            let values: Vec<_> = (0..24)
                .map(|i| if i < 3 { -1.0 } else { i as f64 * 0.5 })
                .collect();
            let mut op = Diagonal { rows, values };
            let mut x = RealMatrix::zeros(op.rows.len(), 6).unwrap();
            for j in 0..6 {
                for (i, &g) in op.rows.iter().enumerate() {
                    x[(i, j)] = if g == j {
                        1.0
                    } else {
                        0.002 * ((g * 7 + j * 3 + 1) as f64).sin()
                    };
                }
            }
            let a = solve(
                c,
                &mut op,
                &x,
                None,
                profile(ritz),
                Options {
                    degree: Some(degree),
                    max_sweeps: 40,
                    ..Options::default()
                },
            )
            .unwrap();
            assert_eq!(a.fallback, None, "{:?}", a.sweeps);
            assert_eq!(a.work.bound_h_columns, 24);
            assert_eq!(a.work.eigen_orders.len(), a.sweeps.len() + 1);
            assert_eq!(a.work.qr_calls, a.sweeps.len() as u64 + 1);
            assert_eq!(
                a.work.h_columns(),
                30 + a.sweeps.len() as u64 * degree as u64 * 6 + if ritz { 0 } else { 6 }
            );
            for j in 0..3 {
                assert!((a.eigenvalues[j] + 1.0).abs() < 1e-10);
            }
            if !ritz {
                assert!(a.residuals[..4].iter().all(|r| *r <= 1e-9));
            }
            (a.eigenvalues, a.work.h_columns())
        });
        assert!(outputs.windows(2).all(|v| v[0] == v[1]));
    }
    #[test]
    fn filter_matches_degenerate_spectrum_on_uneven_and_empty_ranks() {
        for p in [1, 2, 5, 25] {
            for d in [4, 8, 12, 16] {
                run(p, false, d);
            }
        }
    }
    #[test]
    fn captured_ritz_profile_is_distinct_from_residual_profile() {
        run(3, true, 8);
    }
    #[test]
    fn rank_loss_and_missing_edge_request_fallback() {
        dft_mpi::simulate(2, |c| {
            let mut op = Diagonal {
                rows: (0..8).filter(|g| g % 2 == c.rank()).collect(),
                values: vec![1.0; 8],
            };
            let x = RealMatrix::zeros(op.rows.len(), 6).unwrap();
            let a = solve(c, &mut op, &x, None, profile(false), Options::default()).unwrap();
            assert_eq!(a.fallback, Some("initial_rank_loss"));
            assert_eq!(a.work.qr_calls, 1);
        });
        assert_eq!(degree_from_spectrum(-1., 10., 0., 0., 1., 1e-8), None);
    }
    #[test]
    fn degree_is_shift_and_scale_invariant() {
        let a = degree_from_spectrum(-2., 10., 2., -1., 0.1, 1e-6);
        let b = degree_from_spectrum(1., 25., 9., 3., 0.2, 2e-6);
        assert_eq!(a, b);
        assert!(a.is_some());
    }
}
