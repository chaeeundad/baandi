//! One bounded RMM-DIIS sweep for the reference Gamma-real SCF driver.
//!
//! Derived solely from `docs/research/qe-study/qe-rmm-davidson.md` (own-words
//! note). Hartree energies and orthonormal real coordinates throughout. This
//! is deliberately NOT a converged eigensolve or a spectral certificate.
use crate::ref_davidson::{self as dav, ReducedPlacement, Solution};
use crate::{SolverError, distributed_davidson::DistributedRealOperator, real_davidson::lean};
use dav::accounting::{ApplyHScope, Counted, clock, elapsed};
use dav::workspace::{Workspace, agree_apply, rotate};
use dft_mpi::Collective;
use dft_numerics::{GemmSite, RealGemmWorkspace, RealMatrix, RealOp, real::RealEighWorkspace};

/// Same iteration index as the outer SCF loop, including a first-step redo.
pub fn scheduled(iteration: usize) -> bool {
    iteration >= 3 && !iteration.is_multiple_of(5)
}

#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub bands: usize,
    pub threshold_ha: f64,
    pub occupied: usize,
    pub steps: usize,
    pub block: usize,
    pub profile: bool,
    /// Reporting context only; never used to choose the solver here.
    pub scf_iteration: usize,
}
impl Options {
    pub fn new(bands: usize, threshold_ha: f64, occupied: usize) -> Self {
        Self {
            bands,
            threshold_ha,
            occupied,
            steps: 4,
            block: 16,
            profile: false,
            scf_iteration: 0,
        }
    }
}

pub const STAGES: [&str; 7] = [
    "rotation",
    "dense_solve",
    "diis",
    "precondition",
    "apply_h",
    "line_search",
    "gram_schmidt",
];
/// Full-diagnostic rank-local record. GEMM FLOPs exclude scalar band work and
/// backend EVD internals; H columns include exactly the rotation and directions.
#[derive(Clone, Debug, Default)]
pub struct Record {
    pub scf_iteration: usize,
    pub seconds: [f64; 7],
    pub total_seconds: f64,
    pub gemm_flops: [u64; 7],
    pub rotation_columns: usize,
    pub step_columns: Vec<usize>,
    pub newest_fallbacks: usize,
    pub clamped_low: usize,
    pub clamped_high: usize,
    pub gs_fallbacks: usize,
    pub min_gs_norm: f64,
    pub unconverged: usize,
    pub communication: dav::Communication,
}
thread_local! {
    static RECORDS: std::cell::RefCell<Vec<Record>> = const { std::cell::RefCell::new(Vec::new()) };
}
/// Consume only at a reporting boundary. Disabled solves retain no records.
pub fn take_records() -> Vec<Record> {
    RECORDS.with(|s| std::mem::take(&mut *s.borrow_mut()))
}
fn mpi(e: dft_mpi::MpiError) -> SolverError {
    SolverError::Application(e.to_string())
}
fn invalid() -> SolverError {
    SolverError::InvalidParameter("reference RMM dimensions/threshold/kinetic/start")
}

/// Teter-Payne-Allan / Kresse-Furthmueller diagonal in Ha^-1. Reciprocal
/// evaluation avoids overflowing y^4. The T=0 limit is -1/t for t>0; at
/// G=0 the exact normalized constant state's residual component is zero.
fn precondition(t: f64, kinetic: f64) -> f64 {
    if t > 1.5 * kinetic {
        let q = 1.5 * kinetic / t;
        let p = 8. + q * (12. + q * (18. + 27. * q));
        -2. / t * p / (16. + q * p)
    } else if kinetic > 0. {
        let y = t / (1.5 * kinetic);
        let p = 27. + y * (18. + y * (12. + 8. * y));
        -(4. / 3.) / kinetic * p / (p + 16. * y.powi(4))
    } else {
        0.
    }
}
fn step_length(a: &[f64], record: &mut Record) -> (f64, f64, f64) {
    let [a11, a12, a22, b11, b12, b22] = <[f64; 6]>::try_from(a).unwrap();
    let e0 = a11 / b11;
    let slope = 2. * (a12 * b11 - a11 * b12) / b11.powi(2);
    let e1 = (a11 + a12 + 0.25 * a22) / (b11 + b12 + 0.25 * b22);
    let curvature = 4. * (e1 - e0 - 0.5 * slope);
    // The note's 1e-16 Ry curvature zero test, converted to Ha.
    let raw = if curvature.abs() > 5e-17 {
        -slope / (2. * curvature)
    } else if slope < 0. {
        1.
    } else {
        0.05
    };
    record.clamped_low += usize::from(raw < 0.05);
    record.clamped_high += usize::from(raw > 1.);
    let t = raw.clamp(0.05, 1.);
    let norm2 = b11 + 2. * t * b12 + t * t * b22;
    (t, norm2.sqrt(), (a11 + 2. * t * a12 + t * t * a22) / norm2)
}

struct Diis {
    a: Vec<RealMatrix>,
    b: Vec<RealMatrix>,
    beta: Vec<f64>,
    lambda: Vec<f64>,
    coeff: Vec<f64>,
    metric: RealMatrix,
    u: RealMatrix,
    w: RealMatrix,
    aw: RealMatrix,
    small: RealMatrix,
    y: RealMatrix,
    evd: RealEighWorkspace,
}
impl Diis {
    fn new(n: usize, m: usize) -> Result<Self, SolverError> {
        let matrices = || {
            (0..n)
                .map(|_| RealMatrix::zeros(m, m))
                .collect::<Result<Vec<_>, _>>()
        };
        Ok(Self {
            a: matrices()?,
            b: matrices()?,
            beta: vec![0.; m],
            lambda: vec![0.; m],
            coeff: vec![0.; m],
            metric: RealMatrix::zeros(m, m)?,
            u: RealMatrix::zeros(m, m)?,
            w: RealMatrix::zeros(m, m)?,
            aw: RealMatrix::zeros(m, m)?,
            small: RealMatrix::zeros(m, m)?,
            y: RealMatrix::zeros(m, m)?,
            evd: RealEighWorkspace::new(m)?,
        })
    }
    fn coefficients(
        &mut self,
        band: usize,
        k: usize,
        g: &mut RealGemmWorkspace,
        row: &mut dav::Iteration,
    ) -> Result<bool, SolverError> {
        self.coeff.fill(0.);
        if k == 1 {
            self.coeff[0] = 1.;
            return Ok(false);
        }
        self.metric.reshape_within_storage(k, k)?;
        self.u.reshape_within_storage(k, k)?;
        for j in 0..k {
            for i in 0..k {
                self.metric[(i, j)] = self.b[band][(i, j)];
            }
        }
        self.evd
            .solve_into(&self.metric, &mut self.beta[..k], &mut self.u)?;
        let kept = self.beta[..k].iter().filter(|&&v| v > 1e-14).count();
        if kept <= 1 {
            self.coeff[k - 1] = 1.;
            return Ok(true);
        }
        self.w.reshape_within_storage(k, kept)?;
        let mut col = 0;
        for j in 0..k {
            if self.beta[j] > 1e-14 {
                for i in 0..k {
                    self.w[(i, col)] = self.u[(i, j)] / self.beta[j].sqrt();
                }
                col += 1;
            }
        }
        self.aw.reshape_within_storage(k, kept)?;
        self.small.reshape_within_storage(kept, kept)?;
        self.y.reshape_within_storage(kept, kept)?;
        // Pack the active leading square; history matrices have fixed stride m.
        for j in 0..k {
            for i in 0..k {
                self.metric[(i, j)] = self.a[band][(i, j)];
            }
        }
        product(g, &self.metric, RealOp::Identity, &self.w, &mut self.aw)?;
        row.gemm("rmm_diis_aw", k, kept, k);
        product(g, &self.w, RealOp::Transpose, &self.aw, &mut self.small)?;
        row.gemm("rmm_diis_waw", kept, kept, k);
        symmetrize(&mut self.small);
        self.evd
            .solve_into(&self.small, &mut self.lambda[..kept], &mut self.y)?;
        let root = (0..kept)
            .min_by(|&a, &b| self.lambda[a].abs().total_cmp(&self.lambda[b].abs()))
            .unwrap();
        for i in 0..k {
            self.coeff[i] = (0..kept).map(|j| self.w[(i, j)] * self.y[(j, root)]).sum();
        }
        let mut norm2 = 0.;
        for j in 0..k {
            for i in 0..k {
                norm2 += self.coeff[i] * self.b[band][(i, j)] * self.coeff[j];
            }
        }
        if !norm2.is_finite() || norm2 <= 0. {
            return Err(invalid());
        }
        for v in &mut self.coeff[..k] {
            *v /= norm2.sqrt();
        }
        Ok(false)
    }
}
fn product(
    g: &mut RealGemmWorkspace,
    a: &RealMatrix,
    op: RealOp,
    b: &RealMatrix,
    out: &mut RealMatrix,
) -> Result<(), SolverError> {
    g.multiply_slices_into(
        GemmSite::Other,
        a.as_slice(),
        a.rows(),
        a.cols(),
        a.rows(),
        op,
        b.as_slice(),
        b.cols(),
        b.rows(),
        out.as_mut(),
    )?;
    Ok(())
}
fn symmetrize(a: &mut RealMatrix) {
    for j in 0..a.cols() {
        for i in 0..j {
            let v = 0.5 * (a[(i, j)] + a[(j, i)]);
            a[(i, j)] = v;
            a[(j, i)] = v;
        }
    }
}

/// Apply H afresh, rotate in the measured overlap, then accept one sweep.
/// `kinetic_ha` contains t_G for local ORC rows, NOT the Hamiltonian diagonal.
/// All options and row partition contracts are identical on every rank.
pub fn solve<C: Collective, O: DistributedRealOperator>(
    comm: &C,
    op: &mut O,
    start: &RealMatrix,
    kinetic_ha: &[f64],
    o: Options,
) -> Result<Solution, SolverError> {
    let total = clock(o.profile);
    let c = Counted::new(comm, o.profile);
    let n = o.bands;
    let rows = start.rows();
    let m = o.steps;
    let bad = n == 0
        || n > op.global_len()
        || n != start.cols()
        || rows != op.local_rows().len()
        || o.occupied > n
        || !(1..=16).contains(&m)
        || !(1..=64).contains(&o.block)
        || !o.threshold_ha.is_finite()
        || o.threshold_ha <= 0.
        || !(5. * o.threshold_ha).is_finite()
        || kinetic_ha.len() != rows
        || kinetic_ha.iter().any(|v| !v.is_finite() || *v < 0.)
        || start.as_slice().iter().any(|v| !v.is_finite());
    let mut any = [0.];
    c.tree_sum_f64(&[f64::from(bad)], &mut any).map_err(mpi)?;
    if any[0] != 0. {
        return Err(invalid());
    }
    let apply_scope = ApplyHScope::new(o.profile, m + 1);
    let mut common = dav::Record::new(m, n, o.profile, ReducedPlacement::Replicated);
    for row in &mut common.iterations {
        row.gemms.reserve(2 * n);
    }
    common.finish.gemms.reserve(6 * n.div_ceil(o.block));
    let mut record = Record {
        scf_iteration: o.scf_iteration,
        rotation_columns: n,
        step_columns: Vec::with_capacity(m),
        min_gs_norm: f64::INFINITY,
        ..Record::default()
    };
    let mut opts = dav::Options::new(n, o.threshold_ha, o.occupied);
    opts.profile = o.profile;
    let (mut x, mut hx, mut values, mut g) = {
        let mut w = Workspace::new(op, n, n.max(m))?;
        w.v.reshape_within_storage(rows, n)?;
        w.hv.reshape_within_storage(rows, n)?;
        w.v.as_mut_slice().copy_from_slice(start.as_slice());
        let t = clock(o.profile);
        w.apply(&c, op, start, true, &mut common.setup)?;
        record.seconds[4] += elapsed(t);
        let t = clock(o.profile);
        w.border(&c, 0, n, false, &mut common.setup, o.profile)?;
        if !w.factor(0, n, &mut common.setup, o.profile)? {
            return Err(SolverError::InvalidParameter(
                "reference RMM dependent start",
            ));
        }
        w.border(&c, 0, n, true, &mut common.setup, o.profile)?;
        record.seconds[0] += elapsed(t);
        let t = clock(o.profile);
        w.eigen(&c, 0, n, &mut common.setup, &opts)?;
        record.seconds[1] += elapsed(t);
        let t = clock(o.profile);
        w.finish(n, true, &mut common.setup, o.profile)?;
        record.seconds[0] += elapsed(t);
        for p in &common.setup.gemms {
            record.gemm_flops[if matches!(p.stage, "reduced_transform" | "ritz_coefficients") {
                1
            } else {
                0
            }] += p.flops;
        }
        common.setup.communication = c.snapshot();
        // Drop unused Davidson row scratch before allocating the histories.
        (w.x, w.hx, w.values, w.gemm)
    };
    values.truncate(n);
    let mut phi = RealMatrix::zeros(rows, n * m)?;
    let mut hphi = RealMatrix::zeros(rows, n * m)?;
    let mut residual = RealMatrix::zeros(rows, n * m)?;
    let mut p = RealMatrix::zeros(rows, n)?;
    let mut hp = RealMatrix::zeros(rows, n)?;
    let mut diis = Diis::new(n, m)?;
    let mut active = Vec::with_capacity(n);
    let mut frozen = vec![false; n];
    let mut local = vec![0.; (2 * m * n).max(6 * n)];
    let mut reduced = local.clone();
    let mut changes = if o.profile {
        vec![vec![0.; n]; m]
    } else {
        Vec::new()
    };
    let mut ritz = if o.profile {
        vec![vec![[0.; 2]; n]; m]
    } else {
        Vec::new()
    };
    let mut used = 0;
    // Allocate QR scratch before entering the sweep; no history allocation in it.
    let mut gs = GramSchmidt::new(rows, n, o.block.min(n))?;
    for step in 0..m {
        active.clear();
        active.extend((0..n).filter(|&i| !frozen[i]));
        if active.is_empty() {
            break;
        }
        used += 1;
        let k = step + 1;
        let u = active.len();
        let before = c.snapshot();
        let row = &mut common.iterations[step];
        row.active_bands.extend_from_slice(&active);
        row.basis_before = n;
        record.step_columns.push(u);
        let t = clock(o.profile);
        for (a, &band) in active.iter().enumerate() {
            let new = band * m + step;
            for i in 0..rows {
                phi[(i, new)] = x[(i, band)];
                hphi[(i, new)] = hx[(i, band)];
                residual[(i, new)] = hx[(i, band)] - values[band] * x[(i, band)];
            }
            for j in 0..k {
                let old = band * m + j;
                local[a * 2 * k + j] = (0..rows)
                    .map(|i| residual[(i, new)] * residual[(i, old)])
                    .sum();
                local[a * 2 * k + k + j] = (0..rows).map(|i| phi[(i, new)] * phi[(i, old)]).sum();
            }
        }
        c.tree_sum_f64(&local[..2 * k * u], &mut reduced[..2 * k * u])
            .map_err(mpi)?;
        for (a, &band) in active.iter().enumerate() {
            for j in 0..k {
                let av = reduced[a * 2 * k + j];
                let bv = reduced[a * 2 * k + k + j];
                diis.a[band][(j, step)] = av;
                diis.a[band][(step, j)] = av;
                diis.b[band][(j, step)] = bv;
                diis.b[band][(step, j)] = bv;
            }
            record.newest_fallbacks += usize::from(diis.coefficients(band, k, &mut g, row)?);
            for i in 0..rows {
                let mut xv = 0.;
                let mut hv = 0.;
                let mut rv = 0.;
                for j in 0..k {
                    let v = diis.coeff[j];
                    let col = band * m + j;
                    xv += v * phi[(i, col)];
                    hv += v * hphi[(i, col)];
                    rv += v * residual[(i, col)];
                }
                x[(i, band)] = xv;
                hx[(i, band)] = hv;
                p[(i, a)] = rv;
            }
        }
        record.seconds[2] += elapsed(t);
        let t = clock(o.profile);
        for (a, &band) in active.iter().enumerate() {
            local[a] = (0..rows)
                .map(|i| kinetic_ha[i] * x[(i, band)].powi(2))
                .sum();
        }
        c.tree_sum_f64(&local[..u], &mut reduced[..u])
            .map_err(mpi)?;
        for a in 0..u {
            for i in 0..rows {
                p[(i, a)] *= precondition(kinetic_ha[i], reduced[a]);
            }
        }
        record.seconds[3] += elapsed(t);
        let t = clock(o.profile);
        p.reshape_within_storage(rows, u)?;
        hp.reshape_within_storage(rows, u)?;
        let applied = op.apply_block(&c, &p, &mut hp);
        agree_apply(&c, applied, &hp, o.profile)?;
        row.h_calls = 1;
        row.h_columns = u;
        record.seconds[4] += elapsed(t);
        let t = clock(o.profile);
        local[..6 * u].fill(0.);
        for (a, &band) in active.iter().enumerate() {
            for i in 0..rows {
                let xv = x[(i, band)];
                let hv = hx[(i, band)];
                let pv = p[(i, a)];
                local[6 * a] += xv * hv;
                local[6 * a + 1] += pv * hv;
                local[6 * a + 2] += pv * hp[(i, a)];
                local[6 * a + 3] += xv * xv;
                local[6 * a + 4] += pv * xv;
                local[6 * a + 5] += pv * pv;
            }
        }
        c.tree_sum_f64(&local[..6 * u], &mut reduced[..6 * u])
            .map_err(mpi)?;
        for (a, &band) in active.iter().enumerate() {
            let (alpha, norm, value) = step_length(&reduced[6 * a..6 * a + 6], &mut record);
            if !norm.is_finite() || norm <= 0. || !value.is_finite() {
                return Err(invalid());
            }
            let old = values[band];
            frozen[band] = (value - old).abs() < opts.band_threshold(band);
            if o.profile {
                changes[step][band] = (value - old).abs();
                ritz[step][band] = [old, value];
            }
            values[band] = value;
            for i in 0..rows {
                x[(i, band)] = (x[(i, band)] + alpha * p[(i, a)]) / norm;
                hx[(i, band)] = (hx[(i, band)] + alpha * hp[(i, a)]) / norm;
            }
        }
        record.seconds[5] += elapsed(t);
        record.gemm_flops[2] += row.gemms.iter().map(|p| p.flops).sum::<u64>();
        row.unconverged = frozen.iter().filter(|&&v| !v).count();
        row.communication = c.snapshot().difference(before);
        // Restore storage shape before packing the next active block.
        p.reshape_within_storage(rows, n)?;
    }
    common.iterations.truncate(used);
    changes.truncate(used);
    ritz.truncate(used);
    record.unconverged = frozen.iter().filter(|&&v| !v).count();
    let t = clock(o.profile);
    let before = c.snapshot();
    gs.apply(&c, &mut g, &mut x, &mut hx, &mut common.finish, &mut record)?;
    local[..n].fill(0.);
    for j in 0..n {
        local[j] = (0..rows).map(|i| x[(i, j)] * hx[(i, j)]).sum();
    }
    c.tree_sum_f64(&local[..n], &mut values[..n]).map_err(mpi)?;
    // Deterministic in-place sort; carry HX and flags with every swap.
    for j in 0..n {
        let k = (j..n)
            .min_by(|&a, &b| values[a].total_cmp(&values[b]))
            .unwrap();
        values.swap(j, k);
        frozen.swap(j, k);
        for i in 0..rows {
            x.as_mut_slice().swap(i + j * rows, i + k * rows);
            hx.as_mut_slice().swap(i + j * rows, i + k * rows);
        }
    }
    record.seconds[6] += elapsed(t);
    record.gemm_flops[6] = common.finish.gemms.iter().map(|p| p.flops).sum();
    common.finish.communication = c.snapshot().difference(before);
    common.communication = c.snapshot();
    record.communication = common.communication;
    common.summarize();
    common.apply_h = apply_scope.finish();
    common.total_seconds = elapsed(total);
    record.total_seconds = common.total_seconds;
    if o.profile {
        RECORDS.with(|s| s.borrow_mut().push(record));
    }
    Ok(Solution {
        eigenvalues: values,
        eigenvectors: x,
        hx: Some(hx),
        converged: frozen,
        ritz_changes_ha: changes,
        ritz_values_ha: ritz,
        band_thresholds_ha: (0..n).map(|i| opts.band_threshold(i)).collect(),
        record: common,
    })
}

struct GramSchmidt {
    block: usize,
    q: RealMatrix,
    hq: RealMatrix,
    gram: RealMatrix,
    sum: Vec<f64>,
    coefficients: RealMatrix,
    update: RealMatrix,
    factor: lean::LeanScratch,
}
impl GramSchmidt {
    fn new(rows: usize, n: usize, b: usize) -> Result<Self, SolverError> {
        Ok(Self {
            block: b,
            q: RealMatrix::zeros(rows, b)?,
            hq: RealMatrix::zeros(rows, b)?,
            gram: RealMatrix::zeros(b, n)?,
            sum: vec![0.; b * n],
            coefficients: RealMatrix::zeros(b, n)?,
            update: RealMatrix::zeros(rows, n)?,
            factor: lean::LeanScratch::new(b, b)?,
        })
    }
    #[allow(clippy::too_many_arguments)]
    fn apply<C: Collective>(
        &mut self,
        c: &C,
        g: &mut RealGemmWorkspace,
        x: &mut RealMatrix,
        hx: &mut RealMatrix,
        row: &mut dav::Iteration,
        record: &mut Record,
    ) -> Result<(), SolverError> {
        let r = x.rows();
        let n = x.cols();
        for start in (0..n).step_by(self.block) {
            let b = self.block.min(n - start);
            let end = start + b;
            self.q.reshape_within_storage(r, b)?;
            self.hq.reshape_within_storage(r, b)?;
            self.q
                .as_mut_slice()
                .copy_from_slice(&x.as_slice()[r * start..r * end]);
            self.hq
                .as_mut_slice()
                .copy_from_slice(&hx.as_slice()[r * start..r * end]);
            self.gram.reshape_within_storage(b, b)?;
            product(g, &self.q, RealOp::Transpose, &self.q, &mut self.gram)?;
            row.gemm("rmm_gs_panel_gram", b, b, r);
            c.tree_sum_f64(self.gram.as_slice(), &mut self.sum[..b * b])
                .map_err(mpi)?;
            for j in 0..b {
                for i in 0..b {
                    self.factor.s[(i, j)] = self.sum[i + j * b];
                }
            }
            // Cholesky QR is the positive-diagonal, band-order classical GS
            // factorization. One panel reduction replaces per-band reductions.
            // Reuse the guarded factor adapter; a rejected panel uses literal GS.
            let (ok, _) = lean::extend_factor(g, &mut self.factor, 0, b)?;
            if ok {
                self.coefficients.reshape_within_storage(b, b)?;
                for j in 0..b {
                    record.min_gs_norm = record.min_gs_norm.min(1. / self.factor.t[(j, j)]);
                    for i in 0..b {
                        self.coefficients[(i, j)] = self.factor.t[(j, i)];
                    }
                }
                rotate(
                    g,
                    &self.q,
                    &self.coefficients,
                    &mut self.update,
                    b,
                    b,
                    row,
                    "rmm_gs_panel_x",
                )?;
                self.q
                    .as_mut_slice()
                    .copy_from_slice(self.update.as_slice());
                rotate(
                    g,
                    &self.hq,
                    &self.coefficients,
                    &mut self.update,
                    b,
                    b,
                    row,
                    "rmm_gs_panel_hx",
                )?;
                self.hq
                    .as_mut_slice()
                    .copy_from_slice(self.update.as_slice());
            } else {
                record.gs_fallbacks += 1;
                for j in 0..b {
                    for k in 0..j {
                        self.sum[k] = (0..r).map(|i| self.q[(i, k)] * self.q[(i, j)]).sum();
                    }
                    c.tree_sum_f64(&self.sum[..j], &mut self.gram.as_mut_slice()[..j])
                        .map_err(mpi)?;
                    for k in 0..j {
                        let v = self.gram.as_slice()[k];
                        for i in 0..r {
                            self.q[(i, j)] -= v * self.q[(i, k)];
                            self.hq[(i, j)] -= v * self.hq[(i, k)];
                        }
                    }
                    let local = (0..r).map(|i| self.q[(i, j)].powi(2)).sum::<f64>();
                    c.tree_sum_f64(&[local], &mut self.sum[..1]).map_err(mpi)?;
                    let norm = self.sum[0].sqrt();
                    record.min_gs_norm = record.min_gs_norm.min(norm);
                    if !norm.is_finite() || norm < 1e-16 {
                        return Err(SolverError::InvalidParameter(
                            "reference RMM linearly dependent bands",
                        ));
                    }
                    for i in 0..r {
                        self.q[(i, j)] /= norm;
                        self.hq[(i, j)] /= norm;
                    }
                }
            }
            x.as_mut_slice()[r * start..r * end].copy_from_slice(self.q.as_slice());
            hx.as_mut_slice()[r * start..r * end].copy_from_slice(self.hq.as_slice());
            let later = n - end;
            if later == 0 {
                continue;
            }
            self.gram.reshape_within_storage(b, later)?;
            g.multiply_slices_into(
                GemmSite::Other,
                self.q.as_slice(),
                r,
                b,
                r,
                RealOp::Transpose,
                &x.as_slice()[r * end..],
                later,
                r,
                self.gram.as_mut(),
            )?;
            row.gemm("rmm_gs_cross", b, later, r);
            c.tree_sum_f64(self.gram.as_slice(), &mut self.sum[..b * later])
                .map_err(mpi)?;
            self.coefficients.reshape_within_storage(b, later)?;
            self.coefficients
                .as_mut_slice()
                .copy_from_slice(&self.sum[..b * later]);
            rotate(
                g,
                &self.q,
                &self.coefficients,
                &mut self.update,
                b,
                later,
                row,
                "rmm_gs_later_x",
            )?;
            for (dst, v) in x.as_mut_slice()[r * end..]
                .iter_mut()
                .zip(self.update.as_slice())
            {
                *dst -= v;
            }
            rotate(
                g,
                &self.hq,
                &self.coefficients,
                &mut self.update,
                b,
                later,
                row,
                "rmm_gs_later_hx",
            )?;
            for (dst, v) in hx.as_mut_slice()[r * end..]
                .iter_mut()
                .zip(self.update.as_slice())
            {
                *dst -= v;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    struct Dense {
        h: RealMatrix,
        rows: Vec<usize>,
        counts: Vec<usize>,
        full: Vec<f64>,
    }
    impl DistributedRealOperator for Dense {
        fn global_len(&self) -> usize {
            self.h.rows()
        }
        fn local_rows(&self) -> &[usize] {
            &self.rows
        }
        fn apply_block<C: Collective>(
            &mut self,
            c: &C,
            x: &RealMatrix,
            y: &mut RealMatrix,
        ) -> Result<(), SolverError> {
            for j in 0..x.cols() {
                c.allgatherv_f64(
                    &x.as_slice()[j * x.rows()..(j + 1) * x.rows()],
                    &self.counts,
                    &mut self.full,
                )
                .map_err(mpi)?;
                for (i, &row) in self.rows.iter().enumerate() {
                    y[(i, j)] = (0..self.h.rows())
                        .map(|k| self.h[(row, k)] * self.full[k])
                        .sum();
                }
            }
            Ok(())
        }
    }
    fn fixed(ranks: usize, cluster: bool, profile: bool) -> Vec<(Vec<f64>, usize, u64)> {
        dft_mpi::simulate(ranks, |c| {
            let dim = 24;
            let bands = 6;
            let counts: Vec<_> = (0..ranks)
                .map(|rank| dim / ranks + usize::from(rank < dim % ranks))
                .collect();
            let first: usize = counts[..c.rank()].iter().sum();
            let rows: Vec<_> = (first..first + counts[c.rank()]).collect();
            let mut h = RealMatrix::zeros(dim, dim).unwrap();
            // Orthogonal DCT eigenbasis with a complete degenerate low cluster.
            let q = |i: usize, j: usize| {
                if j == 0 {
                    1. / (dim as f64).sqrt()
                } else {
                    (2. / dim as f64).sqrt()
                        * (std::f64::consts::PI * (i as f64 + 0.5) * j as f64 / dim as f64).cos()
                }
            };
            let e = |j: usize| {
                if cluster && (2..=3).contains(&j) {
                    0.6
                } else {
                    0.2 * j as f64
                }
            };
            for j in 0..dim {
                for i in 0..dim {
                    h[(i, j)] = (0..dim).map(|k| q(i, k) * e(k) * q(j, k)).sum();
                }
            }
            let (oracle, _) = dft_numerics::symmetric_eigh(&h).unwrap();
            let mut x = RealMatrix::zeros(rows.len(), bands).unwrap();
            for j in 0..bands {
                for (i, &row) in rows.iter().enumerate() {
                    x[(i, j)] = q(row, j) + 0.03 * q(row, j + bands);
                }
            }
            let kinetic: Vec<_> = rows.iter().map(|&i| 0.2 + 0.05 * i as f64).collect();
            let mut op = Dense {
                h,
                rows,
                counts,
                full: vec![0.; dim],
            };
            let mut options = Options::new(bands, 1e-16, bands);
            options.block = 2;
            options.profile = profile;
            take_records();
            let mut hcols = 0;
            let mut flops = 0;
            for _ in 0..30 {
                let s = solve(c, &mut op, &x, &kinetic, options).unwrap();
                assert_eq!(
                    s.record.h_columns,
                    bands
                        + s.record
                            .iterations
                            .iter()
                            .map(|r| r.h_columns)
                            .sum::<usize>()
                );
                assert!(s.record.h_columns <= 5 * bands);
                let mut fresh = RealMatrix::zeros(x.rows(), bands).unwrap();
                op.apply_block(c, &s.eigenvectors, &mut fresh).unwrap();
                for (a, b) in fresh
                    .as_slice()
                    .iter()
                    .zip(s.hx.as_ref().unwrap().as_slice())
                {
                    assert!((a - b).abs() < 2e-13);
                }
                hcols += s.record.h_columns;
                flops += s.record.gemm_flops;
                x = s.eigenvectors;
            }
            let mut hx = RealMatrix::zeros(x.rows(), bands).unwrap();
            op.apply_block(c, &x, &mut hx).unwrap();
            let mut local = vec![0.; bands * bands + bands];
            let mut sums = local.clone();
            for j in 0..bands {
                for k in 0..bands {
                    local[j * bands + k] = (0..x.rows()).map(|i| x[(i, j)] * x[(i, k)]).sum();
                }
                local[bands * bands + j] = (0..x.rows()).map(|i| x[(i, j)] * hx[(i, j)]).sum();
            }
            c.tree_sum_f64(&local, &mut sums).unwrap();
            for j in 0..bands {
                assert!(
                    (sums[bands * bands + j] - oracle[j]).abs() < 1e-10,
                    "j={j} E={} oracle={}",
                    sums[bands * bands + j],
                    oracle[j]
                );
                for k in 0..bands {
                    assert!((sums[j * bands + k] - f64::from(j == k)).abs() < 2e-12);
                }
            }
            let records = take_records();
            assert_eq!(records.len(), if profile { 30 } else { 0 });
            if profile {
                assert_eq!(
                    records
                        .iter()
                        .map(|r| r.gemm_flops.iter().sum::<u64>())
                        .sum::<u64>(),
                    flops
                );
                assert!(records.iter().all(|r| r.min_gs_norm > 0.9));
            }
            (sums[bands * bands..].to_vec(), hcols, flops)
        })
    }
    #[test]
    fn fixed_h_dense_oracle_repeatability_and_hx() {
        for cluster in [false, true] {
            for ranks in [1, 3, 8] {
                let first = fixed(ranks, cluster, true);
                let second = fixed(ranks, cluster, false);
                assert_eq!(first, second, "profiling/rank repeatability");
                for row in &first[1..] {
                    assert_eq!(row, &first[0]);
                }
            }
        }
    }
    #[test]
    fn exact_constant_band_freezes_and_empty_ranks_agree_errors() {
        dft_mpi::simulate(8, |c| {
            let counts: Vec<_> = (0..8).map(|rank| usize::from(rank < 3)).collect();
            let rows = if c.rank() < 3 {
                vec![c.rank()]
            } else {
                Vec::new()
            };
            let mut h = RealMatrix::zeros(3, 3).unwrap();
            for i in 0..3 {
                h[(i, i)] = 1. + i as f64;
            }
            let mut x = RealMatrix::zeros(rows.len(), 1).unwrap();
            if c.rank() == 0 {
                x[(0, 0)] = 1.;
            }
            let kinetic: Vec<_> = rows.iter().map(|&i| i as f64).collect();
            let mut op = Dense {
                h,
                rows,
                counts,
                full: vec![0.; 3],
            };
            let o = Options::new(1, 1e-12, 1);
            let s = solve(c, &mut op, &x, &kinetic, o).unwrap();
            assert_eq!(s.eigenvalues, vec![1.]);
            assert_eq!(s.converged, vec![true]);
            assert_eq!(s.record.h_columns, 2);
            assert_eq!(s.record.iterations.len(), 1);
            let mut bad = kinetic;
            if c.rank() == 2 {
                bad[0] = f64::NAN;
            }
            assert!(solve(c, &mut op, &x, &bad, o).is_err());
            // All ranks reject a dependent start, including empty row owners.
            x.as_mut_slice().fill(0.);
            let kinetic = vec![0.; x.rows()];
            assert!(solve(c, &mut op, &x, &kinetic, o).is_err());
        });
    }

    #[test]
    fn guarded_panel_falls_back_to_classical_gs_carrying_hx() {
        dft_mpi::simulate(1, |c| {
            let pool = dft_numerics::pool::WorkerPool::new(1).unwrap();
            let mut g = RealGemmWorkspace::with_pool(3, 2, pool).unwrap();
            let mut x = RealMatrix::new(3, 2, vec![1., 0., 0., 1., 1e-5, 0.]).unwrap();
            let mut hx = RealMatrix::new(3, 2, vec![1., 0., 0., 1., 2e-5, 0.]).unwrap();
            let mut gs = GramSchmidt::new(3, 2, 2).unwrap();
            let mut common = dav::Record::new(0, 2, false, ReducedPlacement::Replicated);
            let mut record = Record {
                min_gs_norm: f64::INFINITY,
                ..Default::default()
            };
            gs.apply(c, &mut g, &mut x, &mut hx, &mut common.finish, &mut record)
                .unwrap();
            assert_eq!(record.gs_fallbacks, 1);
            assert_eq!(record.min_gs_norm, 1e-5);
            assert_eq!(x.as_slice(), &[1., 0., 0., 0., 1., 0.]);
            assert_eq!(hx.as_slice(), &[1., 0., 0., 0., 2., 0.]);
        });
    }

    #[test]
    fn schedule_threshold_units_and_polynomial() {
        assert_eq!(
            (1..=15).filter(|&i| !scheduled(i)).collect::<Vec<_>>(),
            vec![1, 2, 5, 10, 15]
        );
        let o = dav::Options::new(2, 1e-8, 1);
        assert_eq!(o.band_threshold(0), 1e-8);
        assert_eq!(o.band_threshold(1), 5e-6);
        for kinetic in [0.01_f64, 1., 10.] {
            for t in [0., 0.5, 1., 10., 1e5] {
                let y = t / (1.5 * kinetic);
                let p = 27. + 18. * y + 12. * y * y + 8. * y * y * y;
                let expected = -(4. / 3.) / kinetic * p / (p + 16. * y.powi(4));
                assert!((precondition(t, kinetic) - expected).abs() < 1e-12 * expected.abs());
                assert!(
                    (2. * precondition(2. * t, 2. * kinetic) - precondition(t, kinetic)).abs()
                        < 1e-12
                );
            }
        }
        assert_eq!(precondition(0., 0.), 0.);
        assert_eq!(precondition(2., 0.), -0.5);
    }
    #[test]
    fn canonical_diis_newest_fallback_and_constraint() {
        let pool = dft_numerics::pool::WorkerPool::new(1).unwrap();
        let mut g = RealGemmWorkspace::with_pool(4, 4, pool).unwrap();
        let mut d = Diis::new(1, 4).unwrap();
        let mut r = dav::Record::new(0, 4, false, ReducedPlacement::Replicated);
        for j in 0..4 {
            for i in 0..4 {
                d.b[0][(i, j)] = 1.;
                d.a[0][(i, j)] = 1.;
            }
        }
        assert!(d.coefficients(0, 4, &mut g, &mut r.setup).unwrap());
        assert_eq!(d.coeff, vec![0., 0., 0., 1.]);
        for j in 0..4 {
            for i in 0..4 {
                d.b[0][(i, j)] = f64::from(i == j);
                d.a[0][(i, j)] = if i == j { (j + 1) as f64 } else { 0. };
            }
        }
        assert!(!d.coefficients(0, 4, &mut g, &mut r.setup).unwrap());
        assert!((d.coeff[0].abs() - 1.).abs() < 1e-14);
        assert!(d.coeff[1..].iter().all(|v| v.abs() < 1e-14));
        let mut rec = Record::default();
        let (t, n, e) = step_length(&[1., -0.1, 2., 1., 0., 0.1], &mut rec);
        assert!((0.05..=1.).contains(&t));
        assert!((n * n - (1. + 0.1 * t * t)).abs() < 1e-14);
        assert!((e - (1. - 0.2 * t + 2. * t * t) / (1. + 0.1 * t * t)).abs() < 1e-14);
    }
}
