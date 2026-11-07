//! Frozen P33A baseline solver (73f59e3), benchmark wrapper only.
//! Run: cargo test --release -p dft-solvers --test real_workspace -- --ignored --nocapture
#![allow(dead_code)]
use dft_numerics::{
    linalg::{DEPENDENCE_TOLERANCE, LinalgError},
    pool::WorkerPool,
    real::{RealGemmWorkspace, RealMatrix, RealMatrixRef, RealOp, symmetric_eigh},
};
use dft_solvers::{BandTolerances, RealLinearOperator, SolverError};

/// Same accuracy flags and work units as the complex EigenSolution.
#[derive(Debug, Clone, PartialEq)]
pub struct RealEigenSolution {
    pub eigenvalues: Vec<f64>,
    pub eigenvectors: RealMatrix,
    pub residuals: Vec<f64>,
    pub converged: Vec<bool>,
    pub buffer: Vec<bool>,
    pub residuals_certified: bool,
    pub iterations: usize,
    pub history: Vec<f64>,
    /// Vector columns applied, including true-residual certificates.
    pub operator_applications: usize,
    pub restarts: usize,
    pub contracted_restarts: usize,
    pub certification_attempts: usize,
    /// Largest live basis width, for capacity and dispatch validation.
    pub max_basis_dimension: usize,
}

/// Block size is bands; default capacity is min(n, max(max_subspace, 6*bands)).
/// The research entry point replaces 6 with its validated capacity multiplier.
/// `diagonal` arguments are the kinetic preconditioner in ORC order, in Ha.
#[derive(Debug, Clone)]
pub struct RealDavidson {
    pub bands: usize,
    /// Absolute residual norm (strict) or Ritz-value change (standard), in Ha.
    pub tolerance: f64,
    pub max_iterations: usize,
    pub max_subspace: usize,
    pub seed: u64,
}
impl RealDavidson {
    pub fn solve(
        &self,
        operator: &dyn RealLinearOperator,
        initial: Option<&RealMatrix>,
        diagonal: &[f64],
    ) -> Result<RealEigenSolution, SolverError> {
        self.solve_with_band_tolerances(
            operator,
            initial,
            diagonal,
            BandTolerances {
                strict_bands: self.bands,
                buffer_tolerance: self.tolerance,
            },
        )
    }

    /// A successful strict solve always has a fresh true-residual certificate.
    /// Only at the iteration limit may buffers fail their requested threshold.
    pub fn solve_with_band_tolerances(
        &self,
        operator: &dyn RealLinearOperator,
        initial: Option<&RealMatrix>,
        diagonal: &[f64],
        policy: BandTolerances,
    ) -> Result<RealEigenSolution, SolverError> {
        self.solve_with_band_tolerances_and_certification(operator, initial, diagonal, policy, true)
    }

    /// With `certify=false`, ordinary returns accept recurrence residuals.
    /// As in complex Davidson, iteration-limit returns always certify.
    pub fn solve_with_band_tolerances_and_certification(
        &self,
        operator: &dyn RealLinearOperator,
        initial: Option<&RealMatrix>,
        diagonal: &[f64],
        policy: BandTolerances,
        certify: bool,
    ) -> Result<RealEigenSolution, SolverError> {
        self.solve_policy(operator, initial, diagonal, policy, false, certify)
    }

    /// Standard-profile semantics: the last mandatory root is the first empty
    /// root and uses the buffer threshold; never a residual certificate.
    pub fn solve_ritz_changes(
        &self,
        operator: &dyn RealLinearOperator,
        initial: Option<&RealMatrix>,
        diagonal: &[f64],
        policy: BandTolerances,
    ) -> Result<RealEigenSolution, SolverError> {
        self.solve_policy(operator, initial, diagonal, policy, true, false)
    }

    fn solve_policy(
        &self,
        op: &dyn RealLinearOperator,
        initial: Option<&RealMatrix>,
        diagonal: &[f64],
        policy: BandTolerances,
        ritz: bool,
        certify: bool,
    ) -> Result<RealEigenSolution, SolverError> {
        self.solve_with_capacity_per_band(op, initial, diagonal, policy, ritz, certify, 6)
    }

    /// Research capacity control. The existing entry points retain multiplier 6.
    /// Minimum 2 allows one full correction block; values above 12 are rejected.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_with_capacity_per_band(
        &self,
        op: &dyn RealLinearOperator,
        initial: Option<&RealMatrix>,
        diagonal: &[f64],
        policy: BandTolerances,
        ritz: bool,
        certify: bool,
        capacity_per_band: usize,
    ) -> Result<RealEigenSolution, SolverError> {
        if !(2..=12).contains(&capacity_per_band) {
            return Err(SolverError::InvalidParameter(
                "davidson_capacity_per_band must be an integer from 2 to 12 (minimum 2 for a full correction block)",
            ));
        }
        let n = op.dimension();
        validate_bands(n, self.bands)?;
        if diagonal.len() != n {
            return Err(LinalgError::Dimension.into());
        }
        finite(diagonal)?;
        if !self.tolerance.is_finite()
            || self.tolerance <= 0.
            || self.max_iterations == 0
            || policy.strict_bands == 0
            || policy.strict_bands > self.bands
            || !policy.buffer_tolerance.is_finite()
            || policy.buffer_tolerance < self.tolerance
        {
            return Err(SolverError::InvalidParameter(
                "invalid tolerance, iteration limit or strict prefix",
            ));
        }
        let capacity = self
            .max_subspace
            .max(self.bands.saturating_mul(capacity_per_band))
            .min(n);
        let pool = match op.worker_pool() {
            Some(pool) => pool,
            None => WorkerPool::new(op.inner_threads())
                .map_err(|e| SolverError::Application(e.to_string()))?,
        };
        // Existing R1 adapter: 256-row panels, ordered sums, persistent pool.
        // 262144 complex MAC * 8 flops / 2 flops per real MAC = 1048576
        // real MAC per dispatch. This is a work estimate, not a speed claim.
        let mut gemm = RealGemmWorkspace::with_pool(n, capacity, pool)?;
        let mut projection = Projection::new(n, capacity)?;
        let mut rng = Random(self.seed);
        let mut basis = if let Some(x) = initial {
            if x.rows() != n || x.cols() > capacity {
                return Err(LinalgError::Dimension.into());
            }
            projection.run(x, None, &mut gemm)?
        } else {
            free_electron_initial_guess(diagonal, self.bands)?
        };
        while basis.cols() < self.bands {
            let extra = projection.run(
                &rng.block(n, self.bands - basis.cols())?,
                Some(&basis),
                &mut gemm,
            )?;
            if extra.cols() == 0 {
                return Err(SolverError::InvalidParameter(
                    "unable to complete initial basis",
                ));
            }
            basis = append(&basis, &extra)?;
        }
        let threshold = |j| {
            if j < policy.strict_bands - usize::from(ritz) {
                self.tolerance
            } else {
                policy.buffer_tolerance
            }
        };
        let meets_policy = |r: &[f64]| r.iter().enumerate().all(|(j, &x)| x <= threshold(j));
        let meets_strict = |r: &[f64]| {
            r[..policy.strict_bands]
                .iter()
                .all(|&x| x <= self.tolerance)
        };
        let mut history = Vec::with_capacity(self.max_iterations);
        let mut applications = 0;
        let mut work = Trace::new(
            self.bands,
            policy.strict_bands,
            self.tolerance,
            policy.buffer_tolerance,
        );
        let mut restarts = 0;
        let mut contractions = 0;
        let mut certificates = 0;
        let mut products = apply(op, &basis, &mut applications, &mut work)?;
        let mut reduced = None;
        let mut previous = vec![f64::NAN; self.bands];
        let mut max_basis_dimension = 0;
        for iteration in 1..=self.max_iterations {
            max_basis_dimension = max_basis_dimension.max(basis.cols());
            let small = reduced_matrix(&basis, &products, reduced.as_ref(), &mut gemm)?;
            let (mut values, coefficients) = {
                let _timer = Timer::new("eigh");
                symmetric_eigh(&small)?
            };
            reduced = Some(small);
            values.truncate(self.bands);
            let coefficients = coefficients.as_ref().subcols(0, self.bands);
            let vectors = multiply(&basis, RealOp::Identity, coefficients, &mut gemm)?;
            let rotated_products = multiply(&products, RealOp::Identity, coefficients, &mut gemm)?;
            let (mut residual_block, mut residuals) =
                residual(&vectors, &rotated_products, &values)?;
            let mut certified = false;
            let mut trusted = true;
            if !ritz
                && (certify || iteration == self.max_iterations)
                && (meets_policy(&residuals)
                    || (iteration == self.max_iterations && meets_strict(&residuals)))
            {
                certificates += 1;
                let actual = apply(op, &vectors, &mut applications, &mut work)?;
                (residual_block, residuals) = residual(&vectors, &actual, &values)?;
                certified = true;
                trusted = meets_policy(&residuals);
            }
            let maximum = residuals.iter().copied().fold(0., f64::max);
            history.push(maximum);
            work.residuals(iteration, &residuals);
            let converged: Vec<bool> = (0..self.bands)
                .map(|j| {
                    if ritz {
                        iteration > 1 && (values[j] - previous[j]).abs() < threshold(j)
                    } else {
                        residuals[j] <= threshold(j)
                    }
                })
                .collect();
            previous.copy_from_slice(&values);
            let finished = converged.iter().all(|&b| b)
                || (iteration == self.max_iterations
                    && converged[..policy.strict_bands].iter().all(|&b| b));
            if finished {
                work.returned();
                return Ok(RealEigenSolution {
                    eigenvalues: values,
                    eigenvectors: vectors,
                    residuals,
                    converged,
                    buffer: (0..self.bands).map(|j| j >= policy.strict_bands).collect(),
                    residuals_certified: certified,
                    iterations: iteration,
                    history,
                    operator_applications: applications,
                    restarts,
                    contracted_restarts: contractions,
                    certification_attempts: certificates,
                    max_basis_dimension,
                });
            }
            if iteration == self.max_iterations {
                return Err(SolverError::NotConverged {
                    iterations: iteration,
                    max_residual: maximum,
                    residuals,
                });
            }
            if ritz && iteration == 1 && basis.cols() == n {
                continue;
            }
            let selected: Vec<_> = (0..self.bands)
                .filter(|&j| !converged[j])
                .take(self.bands.min(capacity - self.bands))
                .collect();
            let mut corrections = RealMatrix::zeros(n, selected.len())?;
            for (k, &j) in selected.iter().enumerate() {
                for ((dst, &r), &d) in column_mut(&mut corrections, k)
                    .iter_mut()
                    .zip(column(&residual_block, j))
                    .zip(diagonal)
                {
                    *dst = precondition(r, d, values[j])?;
                }
            }
            if !trusted || basis.cols() + corrections.cols() > capacity {
                restarts += 1;
                work.restart();
                reduced = None;
                if trusted && orthogonality_defect(&vectors) <= 1e-12 {
                    contractions += 1;
                    basis = vectors;
                    products = rotated_products;
                } else {
                    basis = projection.run(&vectors, None, &mut gemm)?;
                    products = apply(op, &basis, &mut applications, &mut work)?;
                }
            }
            let mut extra = projection.run(&corrections, Some(&basis), &mut gemm)?;
            if extra.cols() == 0 {
                let random = rng.block(n, self.bands.min(capacity - basis.cols()))?;
                extra = projection.run(&random, Some(&basis), &mut gemm)?;
            }
            if extra.cols() == 0 {
                return Err(SolverError::NotConverged {
                    iterations: iteration,
                    max_residual: maximum,
                    residuals,
                });
            }
            products = append(&products, &apply(op, &extra, &mut applications, &mut work)?)?;
            basis = append(&basis, &extra)?;
        }
        unreachable!("validated positive iteration limit")
    }
}

/// Lowest kinetic coordinate units; ties retain ORC half-list order (0, cos, sin).
/// The diagonal is kinetic, not the full Hamiltonian diagonal.
pub fn free_electron_initial_guess(
    kinetic: &[f64],
    bands: usize,
) -> Result<RealMatrix, SolverError> {
    validate_bands(kinetic.len(), bands)?;
    finite(kinetic)?;
    let mut indices: Vec<_> = (0..kinetic.len()).collect();
    indices.sort_by(|&i, &j| kinetic[i].total_cmp(&kinetic[j]).then(i.cmp(&j)));
    let mut x = RealMatrix::zeros(kinetic.len(), bands)?;
    for (j, &i) in indices[..bands].iter().enumerate() {
        x[(i, j)] = 1.;
    }
    Ok(x)
}
fn finite(x: &[f64]) -> Result<(), LinalgError> {
    if x.iter().all(|x| x.is_finite()) {
        Ok(())
    } else {
        Err(LinalgError::NonFinite)
    }
}
fn column(x: &RealMatrix, j: usize) -> &[f64] {
    &x.as_slice()[j * x.rows()..(j + 1) * x.rows()]
}
fn column_mut(x: &mut RealMatrix, j: usize) -> &mut [f64] {
    let rows = x.rows();
    &mut x.as_mut_slice()[j * rows..(j + 1) * rows]
}
fn norm(x: &[f64]) -> f64 {
    let mut sums = [0.; 4];
    for chunk in x.chunks(4) {
        for (k, &v) in chunk.iter().enumerate() {
            sums[k] += v * v;
        }
    }
    let sum = (sums[0] + sums[1]) + (sums[2] + sums[3]);
    if sum.is_finite() && sum > f64::MIN_POSITIVE * (x.len() as f64).max(1.) {
        sum.sqrt()
    } else {
        x.iter().fold(0_f64, |a, &b| a.hypot(b))
    }
}
fn dot(x: &[f64], y: &[f64]) -> f64 {
    let mut sums = [0.; 4];
    for (a, b) in x.chunks(4).zip(y.chunks(4)) {
        for k in 0..a.len() {
            sums[k] += a[k] * b[k];
        }
    }
    (sums[0] + sums[1]) + (sums[2] + sums[3])
}
fn precondition(r: f64, d: f64, lambda: f64) -> Result<f64, LinalgError> {
    let delta = d - lambda;
    if !delta.is_finite() {
        return Err(LinalgError::NonFinite);
    }
    // Existing kinetic floor, in Ha; zero has positive sign.
    let denominator = if delta < 0. {
        -delta.abs().max(0.1)
    } else {
        delta.max(0.1)
    };
    let result = -r / denominator;
    finite(&[result])?;
    Ok(result)
}
fn prefix(x: &RealMatrix, count: usize) -> Result<RealMatrix, LinalgError> {
    let _timer = Timer::new("prefix");
    RealMatrix::new(x.rows(), count, x.as_slice()[..x.rows() * count].to_vec())
}
fn append(a: &RealMatrix, b: &RealMatrix) -> Result<RealMatrix, LinalgError> {
    let _timer = Timer::new("append");
    let mut data = Vec::with_capacity(a.as_slice().len() + b.as_slice().len());
    data.extend_from_slice(a.as_slice());
    data.extend_from_slice(b.as_slice());
    RealMatrix::new(a.rows(), a.cols() + b.cols(), data)
}
fn apply(
    op: &dyn RealLinearOperator,
    x: &RealMatrix,
    count: &mut usize,
    work: &mut Trace,
) -> Result<RealMatrix, SolverError> {
    let _timer = Timer::new("apply");
    work.apply(x.cols(), None);
    for pair in (0..x.cols()).step_by(2) {
        let a = norm(column(x, pair));
        let b = if pair + 1 < x.cols() {
            norm(column(x, pair + 1))
        } else {
            a
        };
        debug_assert!(
            a > 0. && b > 0. && a.max(b) <= 1e3 * a.min(b),
            "paired H columns must have comparable norms"
        );
    }
    let mut y = RealMatrix::zeros(x.rows(), x.cols())?;
    op.apply_block(x, &mut y)?;
    if y.rows() != x.rows() || y.cols() != x.cols() {
        return Err(LinalgError::Dimension.into());
    }
    finite(y.as_slice())?;
    *count += x.cols();
    Ok(y)
}
fn multiply(
    a: &RealMatrix,
    op: RealOp,
    b: RealMatrixRef<'_>,
    ws: &mut RealGemmWorkspace,
) -> Result<RealMatrix, LinalgError> {
    let _timer = Timer::new("multiply");
    let rows = if op == RealOp::Identity {
        a.rows()
    } else {
        a.cols()
    };
    let mut out = RealMatrix::zeros(rows, b.ncols())?;
    timed_gemm(ws, a.as_ref(), op, b, out.as_mut())?;
    Ok(out)
}
fn reduced_matrix(
    v: &RealMatrix,
    hv: &RealMatrix,
    old: Option<&RealMatrix>,
    ws: &mut RealGemmWorkspace,
) -> Result<RealMatrix, LinalgError> {
    let _timer = Timer::new("reduced_matrix");
    let m = v.cols();
    let prefix = old.map_or(0, RealMatrix::cols);
    let mut out = RealMatrix::zeros(m, m)?;
    if let Some(old) = old {
        for j in 0..prefix {
            column_mut(&mut out, j)[..prefix].copy_from_slice(column(old, j));
        }
    }
    // One old/new rectangle and the new square, in one GEMM. Preserve the
    // old prefix exactly; the newly computed upper triangle defines symmetry.
    timed_gemm(
        ws,
        v.as_ref(),
        RealOp::Transpose,
        hv.as_ref().subcols(prefix, m - prefix),
        out.as_mut().subcols_mut(prefix, m - prefix),
    )?;
    for j in prefix..m {
        let (earlier, current) = out.as_mut_slice().split_at_mut(j * m);
        for (col, &entry) in earlier.chunks_exact_mut(m).zip(&current[..j]) {
            col[j] = entry;
        }
    }
    Ok(out)
}
fn residual(
    v: &RealMatrix,
    hv: &RealMatrix,
    values: &[f64],
) -> Result<(RealMatrix, Vec<f64>), LinalgError> {
    let _timer = Timer::new("residual");
    let mut r = hv.clone();
    let mut norms = Vec::with_capacity(values.len());
    for (j, &lambda) in values.iter().enumerate() {
        for (r, &v) in column_mut(&mut r, j).iter_mut().zip(column(v, j)) {
            *r -= lambda * v;
        }
        norms.push(norm(column(&r, j)));
    }
    finite(&norms)?;
    Ok((r, norms))
}
fn orthogonality_defect(v: &RealMatrix) -> f64 {
    let _timer = Timer::new("orthogonality_defect");
    let mut defect = 0_f64;
    for j in 0..v.cols() {
        for i in 0..=j {
            let error = (dot(column(v, i), column(v, j)) - if i == j { 1. } else { 0. }).abs();
            if !error.is_finite() {
                return f64::INFINITY;
            }
            defect = defect.max(error);
        }
    }
    defect
}

// Preallocated projection buffers; each pass is Q^T X, Q times overlaps,
// subtraction, then two-pass MGS and normalization within the incoming block.
struct Projection {
    normalized: RealMatrix,
    overlap: RealMatrix,
    projected: RealMatrix,
    output: RealMatrix,
    column: Vec<f64>,
}
impl Projection {
    fn new(n: usize, cap: usize) -> Result<Self, LinalgError> {
        Ok(Self {
            normalized: RealMatrix::zeros(n, cap)?,
            overlap: RealMatrix::zeros(cap, cap)?,
            projected: RealMatrix::zeros(n, cap)?,
            output: RealMatrix::zeros(n, cap)?,
            column: vec![0.; n],
        })
    }
    fn run(
        &mut self,
        x: &RealMatrix,
        against: Option<&RealMatrix>,
        ws: &mut RealGemmWorkspace,
    ) -> Result<RealMatrix, LinalgError> {
        let _timer = Timer::new("run");
        finite(x.as_slice())?;
        let n = x.rows();
        let mut count = x.cols();
        for j in 0..count {
            let len = norm(column(x, j));
            finite(&[len])?;
            for (dst, &src) in column_mut(&mut self.normalized, j)
                .iter_mut()
                .zip(column(x, j))
            {
                *dst = if len == 0. { 0. } else { src / len };
            }
        }
        let against = against.filter(|q| q.cols() > 0);
        for pass in 0..if against.is_some() { 2 } else { 1 } {
            if pass > 0 {
                // Both buffers have the same capacity; the previous MGS result
                // becomes input, and the old input storage becomes output.
                std::mem::swap(&mut self.normalized, &mut self.output);
            }
            if let Some(q) = against {
                timed_gemm(
                    ws,
                    q.as_ref(),
                    RealOp::Transpose,
                    self.normalized.as_ref().subcols(0, count),
                    self.overlap.as_mut().submatrix_mut(0, 0, q.cols(), count),
                )?;
                timed_gemm(
                    ws,
                    q.as_ref(),
                    RealOp::Identity,
                    self.overlap.as_ref().submatrix(0, 0, q.cols(), count),
                    self.projected.as_mut().subcols_mut(0, count),
                )?;
                for (dst, src) in self.normalized.as_mut_slice()[..n * count]
                    .iter_mut()
                    .zip(&self.projected.as_slice()[..n * count])
                {
                    *dst -= src;
                }
            }
            {
                let _timer = Timer::new("mgs");
                let mut kept = 0;
                for j in 0..count {
                    self.column.copy_from_slice(column(&self.normalized, j));
                    for _ in 0..2 {
                        for k in 0..kept {
                            let p = dot(&self.column, column(&self.output, k));
                            for (dst, &src) in self.column.iter_mut().zip(column(&self.output, k)) {
                                *dst -= p * src;
                            }
                        }
                    }
                    let len = norm(&self.column);
                    finite(&[len])?;
                    // Same relative rank floor as complex MGS; no CholeskyQR2.
                    if len > DEPENDENCE_TOLERANCE {
                        for (dst, &src) in column_mut(&mut self.output, kept)
                            .iter_mut()
                            .zip(&self.column)
                        {
                            *dst = src / len;
                        }
                        kept += 1;
                    }
                }
                count = kept;
            }
        }
        prefix(&self.output, count)
    }
}
struct Random(u64);
impl Random {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        (((z ^ (z >> 31)) >> 11) as f64) * (1. / ((1_u64 << 53) as f64)) - 0.5
    }
    fn block(&mut self, n: usize, m: usize) -> Result<RealMatrix, LinalgError> {
        RealMatrix::new(n, m, (0..n * m).map(|_| self.next()).collect())
    }
}

fn validate_bands(n: usize, bands: usize) -> Result<(), SolverError> {
    if bands == 0 || bands > n {
        return Err(SolverError::InvalidParameter("invalid bands"));
    }
    Ok(())
}
struct Trace;
impl Trace {
    fn new(_: usize, _: usize, _: f64, _: f64) -> Self {
        Self
    }
    fn apply(&mut self, _: usize, _: Option<usize>) {}
    fn residuals(&mut self, _: usize, _: &[f64]) {}
    fn returned(&mut self) {}
    fn restart(&mut self) {}
}
use std::{
    cell::RefCell,
    collections::BTreeMap,
    sync::Arc,
    time::{Duration, Instant},
};
#[derive(Default)]
struct Timings {
    stack: Vec<Duration>,
    totals: BTreeMap<&'static str, (usize, Duration)>,
}
thread_local! { static TIMES: RefCell<Timings> = RefCell::new(Timings::default()); }
struct Timer {
    name: &'static str,
    start: Instant,
}
impl Timer {
    fn new(name: &'static str) -> Self {
        TIMES.with_borrow_mut(|t| t.stack.push(Duration::ZERO));
        Self {
            name,
            start: Instant::now(),
        }
    }
}
impl Drop for Timer {
    fn drop(&mut self) {
        let elapsed = self.start.elapsed();
        TIMES.with_borrow_mut(|t| {
            let child = t.stack.pop().unwrap();
            if let Some(parent) = t.stack.last_mut() {
                *parent += elapsed;
            }
            let total = t.totals.entry(self.name).or_default();
            total.0 += 1;
            total.1 += elapsed.saturating_sub(child);
        });
    }
}
struct Synthetic {
    n: usize,
    pool: Arc<WorkerPool>,
}
impl RealLinearOperator for Synthetic {
    fn dimension(&self) -> usize {
        self.n
    }
    fn worker_pool(&self) -> Option<Arc<WorkerPool>> {
        Some(self.pool.clone())
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), SolverError> {
        for (i, dst) in y.iter_mut().enumerate() {
            *dst = (1. + i as f64 * 0.01) * x[i];
            if i > 0 {
                *dst += 0.03 * x[i - 1];
            }
            if i + 1 < self.n {
                *dst += 0.03 * x[i + 1];
            }
        }
        Ok(())
    }
}
#[test]
#[ignore = "large controlled operator; release benchmark; timers only in this wrapper"]
fn routine_split() {
    for (n, bands) in [(2969, 20), (11877, 68)] {
        for threads in [1, 2, 4, 8, 12] {
            let op = Synthetic {
                n,
                pool: WorkerPool::new(threads).unwrap(),
            };
            let diagonal: Vec<_> = (0..n).map(|i| 1. + i as f64 * 0.01).collect();
            let solver = RealDavidson {
                bands,
                tolerance: 1e-10,
                max_iterations: 400,
                max_subspace: 6 * bands,
                seed: 42,
            };
            TIMES.with_borrow_mut(|t| *t = Timings::default());
            let start = Instant::now();
            let mut last = None;
            for _ in 0..5 {
                last = Some(solver.solve(&op, None, &diagonal).unwrap());
            }
            let elapsed = start.elapsed();
            let result = last.unwrap();
            println!(
                "SHAPE n={n} bands={bands} threads={threads} solves=5 seconds={} steps={} H={}",
                elapsed.as_secs_f64(),
                result.iterations,
                result.operator_applications
            );
            TIMES.with_borrow(|t| {
                for (name, (calls, duration)) in &t.totals {
                    println!(
                        "ROUTINE {name} calls={calls} seconds={:.9} share={:.6}",
                        duration.as_secs_f64(),
                        duration.as_secs_f64() / elapsed.as_secs_f64()
                    );
                }
            });
            let current = dft_solvers::RealDavidson {
                bands,
                tolerance: 1e-10,
                max_iterations: 400,
                max_subspace: 6 * bands,
                seed: 42,
            }
            .solve(&op, None, &diagonal)
            .unwrap();
            assert_eq!(result.eigenvalues, current.eigenvalues);
            assert_eq!(result.eigenvectors, current.eigenvectors);
            assert_eq!(result.residuals, current.residuals);
            assert_eq!(result.history, current.history);
            assert_eq!(result.iterations, current.iterations);
            assert_eq!(result.operator_applications, current.operator_applications);
        }
    }
}

fn timed_gemm(
    ws: &mut RealGemmWorkspace,
    a: RealMatrixRef<'_>,
    op: RealOp,
    b: RealMatrixRef<'_>,
    out: dft_numerics::real::RealMatrixMut<'_>,
) -> Result<(), LinalgError> {
    let _timer = Timer::new(if op == RealOp::Transpose {
        "gram"
    } else {
        "rotate"
    });
    ws.multiply_into(a, op, b, out)
}

// Test-only allocator adapter explicitly requested by the P33A allocation gate.
// Production code contains no unsafe and no allocation instrumentation.
#[allow(unsafe_code)]
mod allocation_probe {
    use std::{
        alloc::{GlobalAlloc, Layout, System},
        sync::atomic::{AtomicBool, AtomicUsize, Ordering},
    };
    pub static TRACK: AtomicBool = AtomicBool::new(false);
    pub static COUNT: AtomicUsize = AtomicUsize::new(0);
    pub struct Allocator;
    fn record() {
        if TRACK.load(Ordering::Relaxed) {
            COUNT.fetch_add(1, Ordering::Relaxed);
        }
    }
    // SAFETY: all pointers/layouts are forwarded unchanged to the System allocator;
    // counters use atomics and never allocate or alter allocation semantics.
    unsafe impl GlobalAlloc for Allocator {
        unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
            record();
            // SAFETY: caller supplies a valid allocation layout.
            unsafe { System.alloc(layout) }
        }
        unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
            record();
            // SAFETY: caller supplies a valid allocation layout.
            unsafe { System.alloc_zeroed(layout) }
        }
        unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
            // SAFETY: pointer and original layout come from this forwarding allocator.
            unsafe { System.dealloc(ptr, layout) }
        }
        unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, size: usize) -> *mut u8 {
            record();
            // SAFETY: pointer/layout/new size retain the GlobalAlloc caller contract.
            unsafe { System.realloc(ptr, layout, size) }
        }
    }
}
#[global_allocator]
static ALLOCATOR: allocation_probe::Allocator = allocation_probe::Allocator;
struct CountingOperator(Synthetic);
impl RealLinearOperator for CountingOperator {
    fn dimension(&self) -> usize {
        self.0.n
    }
    fn worker_pool(&self) -> Option<Arc<WorkerPool>> {
        Some(self.0.pool.clone())
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), SolverError> {
        self.0.apply(x, y)
    }
    fn apply_block(&self, x: &RealMatrix, y: &mut RealMatrix) -> Result<(), SolverError> {
        for (x, y) in x
            .as_slice()
            .chunks_exact(self.0.n)
            .zip(y.as_mut_slice().chunks_exact_mut(self.0.n))
        {
            self.apply(x, y)?;
        }
        // Initial H application is after all solver workspace construction.
        // Track from its return through every iteration, certification and return.
        allocation_probe::TRACK.store(true, std::sync::atomic::Ordering::SeqCst);
        Ok(())
    }
}
#[test]
#[ignore = "global allocator probe: release, --test-threads=1; includes large dimensions"]
fn no_allocations_after_initial_h() {
    use allocation_probe::{COUNT, TRACK};
    use std::sync::atomic::Ordering;
    let mut violations = Vec::new();
    for (n, bands, warm, ritz) in [
        (2969, 20, false, false),
        (11877, 68, false, false),
        (97, 8, true, false),
        (97, 8, true, true),
    ] {
        for threads in [1, 2, 4, 8, 12] {
            let op = CountingOperator(Synthetic {
                n,
                pool: WorkerPool::new(threads).unwrap(),
            });
            let diagonal: Vec<_> = (0..n).map(|i| 1. + i as f64 * 0.01).collect();
            let initial = warm.then(|| Random(17).block(n, bands).unwrap());
            let solver = dft_solvers::RealDavidson {
                bands,
                tolerance: 1e-10,
                max_iterations: 400,
                max_subspace: 6 * bands,
                seed: 42,
            };
            COUNT.store(0, Ordering::SeqCst);
            let result = if ritz {
                solver.solve_ritz_changes(
                    &op,
                    initial.as_ref(),
                    &diagonal,
                    BandTolerances {
                        strict_bands: bands,
                        buffer_tolerance: 1e-10,
                    },
                )
            } else {
                solver.solve(&op, initial.as_ref(), &diagonal)
            };
            TRACK.store(false, Ordering::SeqCst);
            let count = COUNT.load(Ordering::SeqCst);
            let result = result.unwrap();
            println!(
                "ALLOC n={n} bands={bands} threads={threads} warm={warm} ritz={ritz} allocations={count} steps={} restarts={} certificates={}",
                result.iterations, result.restarts, result.certification_attempts
            );
            if count != 0 {
                violations.push((n, bands, threads, count));
            }
        }
    }
    assert!(
        violations.is_empty(),
        "backend allocations remain: {violations:?}"
    );
}

#[test]
#[ignore = "allocation localization diagnostic, release and one test thread"]
fn localize_allocations() {
    use allocation_probe::{COUNT, TRACK};
    use std::sync::atomic::Ordering;
    for n in [20, 40, 88, 160] {
        let mut a = Random(17).block(n, n).unwrap();
        for j in 0..n {
            for i in 0..j {
                a[(i, j)] = a[(j, i)];
            }
        }
        let mut ws = dft_numerics::real::RealEighWorkspace::new(n).unwrap();
        let mut v = RealMatrix::zeros(n, n).unwrap();
        let mut values = vec![0.; n];
        COUNT.store(0, Ordering::SeqCst);
        TRACK.store(true, Ordering::SeqCst);
        let result = ws.solve_into(&a, &mut values, &mut v);
        TRACK.store(false, Ordering::SeqCst);
        result.unwrap();
        println!("EVD n={n} allocations={}", COUNT.load(Ordering::SeqCst));
    }
    for op in [RealOp::Identity, RealOp::Transpose] {
        let a = Random(17).block(2969, 40).unwrap();
        let b = Random(18)
            .block(if op == RealOp::Identity { 40 } else { 2969 }, 20)
            .unwrap();
        let mut out =
            RealMatrix::zeros(if op == RealOp::Identity { 2969 } else { 40 }, 20).unwrap();
        let mut ws = RealGemmWorkspace::with_pool(2969, 40, WorkerPool::new(1).unwrap()).unwrap();
        COUNT.store(0, Ordering::SeqCst);
        TRACK.store(true, Ordering::SeqCst);
        let result = ws.multiply_into(a.as_ref(), op, b.as_ref(), out.as_mut());
        TRACK.store(false, Ordering::SeqCst);
        result.unwrap();
        println!(
            "GEMM op={op:?} allocations={}",
            COUNT.load(Ordering::SeqCst)
        );
    }
}
