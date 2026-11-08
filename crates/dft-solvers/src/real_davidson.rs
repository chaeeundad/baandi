//! Isolated ORC Davidson; no complex-path dispatch or arithmetic is changed.
use crate::{BandTolerances, RealLinearOperator, SolverError};
use dft_numerics::{
    linalg::{DEPENDENCE_TOLERANCE, LinalgError},
    pool::WorkerPool,
    real::{RealEighWorkspace, RealGemmWorkspace, RealMatrix, RealOp},
};

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
/// Persistent scratch only: each solve resets its mathematical state. Backend
/// EVD may still allocate internally (see P33A); this does not promise otherwise.
pub struct RealDavidsonWorkspace {
    n: usize,
    bands: usize,
    capacity: usize,
    pool: std::sync::Arc<WorkerPool>,
    gemm: RealGemmWorkspace,
    projection: Projection,
    basis: RealMatrix,
    products: RealMatrix,
    reduced: RealMatrix,
    eigensolver: RealEighWorkspace,
    values: Vec<f64>,
    coefficients: RealMatrix,
    vectors: RealMatrix,
    rotated_products: RealMatrix,
    residual_block: RealMatrix,
    residuals: Vec<f64>,
    actual: RealMatrix,
    corrections: RealMatrix,
    extra: RealMatrix,
    extra_products: RealMatrix,
    random: RealMatrix,
    selected: Vec<usize>,
    converged: Vec<bool>,
    previous: Vec<f64>,
}
impl RealDavidsonWorkspace {
    fn new(
        n: usize,
        bands: usize,
        capacity: usize,
        pool: std::sync::Arc<WorkerPool>,
    ) -> Result<Self, SolverError> {
        Ok(Self {
            gemm: RealGemmWorkspace::with_pool(n, capacity, pool.clone())?,
            projection: Projection::with_pool(n, capacity, pool.clone())?,
            basis: RealMatrix::zeros(n, capacity)?,
            products: RealMatrix::zeros(n, capacity)?,
            reduced: RealMatrix::zeros(capacity, capacity)?,
            eigensolver: RealEighWorkspace::new(capacity)?,
            values: vec![0.; capacity],
            coefficients: RealMatrix::zeros(capacity, capacity)?,
            vectors: RealMatrix::zeros(n, capacity)?,
            rotated_products: RealMatrix::zeros(n, capacity)?,
            residual_block: RealMatrix::zeros(n, bands)?,
            residuals: vec![0.; bands],
            actual: RealMatrix::zeros(n, bands)?,
            corrections: RealMatrix::zeros(n, bands)?,
            extra: RealMatrix::zeros(n, capacity)?,
            extra_products: RealMatrix::zeros(n, capacity)?,
            random: RealMatrix::zeros(n, bands)?,
            selected: Vec::with_capacity(bands),
            converged: vec![false; bands],
            previous: vec![f64::NAN; bands],
            n,
            bands,
            capacity,
            pool,
        })
    }
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
        self.solve_reusing(
            op,
            initial,
            diagonal,
            policy,
            ritz,
            certify,
            capacity_per_band,
            &mut None,
        )
    }

    /// Reuse storage across solves; shape or executor changes replace it. Options
    /// and starting vectors are read afresh on every call, never cached here.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_reusing(
        &self,
        op: &dyn RealLinearOperator,
        initial: Option<&RealMatrix>,
        diagonal: &[f64],
        policy: BandTolerances,
        ritz: bool,
        certify: bool,
        capacity_per_band: usize,
        workspace: &mut Option<RealDavidsonWorkspace>,
    ) -> Result<RealEigenSolution, SolverError> {
        self.solve_reusing_with_preconditioner(
            op,
            initial,
            diagonal,
            policy,
            ritz,
            certify,
            capacity_per_band,
            workspace,
            crate::DavidsonPreconditioner::Default,
        )
    }

    /// Research damping selection; the old entry point is bit-preserving.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_reusing_with_preconditioner(
        &self,
        op: &dyn RealLinearOperator,
        initial: Option<&RealMatrix>,
        diagonal: &[f64],
        policy: BandTolerances,
        ritz: bool,
        certify: bool,
        capacity_per_band: usize,
        workspace: &mut Option<RealDavidsonWorkspace>,
        preconditioner: crate::DavidsonPreconditioner,
    ) -> Result<RealEigenSolution, SolverError> {
        if !(2..=12).contains(&capacity_per_band) {
            return Err(SolverError::InvalidParameter(
                "davidson_capacity_per_band must be an integer from 2 to 12 (minimum 2 for a full correction block)",
            ));
        }
        let n = op.dimension();
        crate::validate_bands(n, self.bands)?;
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
        if workspace.as_ref().is_none_or(|w| {
            w.n != n
                || w.bands != self.bands
                || w.capacity != capacity
                || !std::sync::Arc::ptr_eq(&w.pool, &pool)
        }) {
            *workspace = Some(RealDavidsonWorkspace::new(
                n,
                self.bands,
                capacity,
                pool.clone(),
            )?);
        }
        let RealDavidsonWorkspace {
            gemm,
            projection,
            basis,
            products,
            reduced,
            eigensolver,
            values,
            coefficients,
            vectors,
            rotated_products,
            residual_block,
            residuals,
            actual,
            corrections,
            extra,
            extra_products,
            random,
            selected,
            converged,
            previous,
            ..
        } = workspace.as_mut().unwrap();
        previous.fill(f64::NAN);
        selected.clear();
        let mut rng = Random(self.seed);
        if let Some(x) = initial {
            if x.rows() != n || x.cols() > capacity {
                return Err(LinalgError::Dimension.into());
            }
            projection.run_into(x, None, gemm, basis)?;
        } else {
            let starting_basis = free_electron_initial_guess(diagonal, self.bands)?;
            basis.reshape_within_storage(n, starting_basis.cols())?;
            basis
                .as_mut_slice()
                .copy_from_slice(starting_basis.as_slice());
        }
        while basis.cols() < self.bands {
            let extra =
                projection.run(&rng.block(n, self.bands - basis.cols())?, Some(basis), gemm)?;
            if extra.cols() == 0 {
                return Err(SolverError::InvalidParameter(
                    "unable to complete initial basis",
                ));
            }
            append_into(basis, &extra)?;
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
        let mut work = crate::work_accounting::Trace::new(
            self.bands,
            policy.strict_bands,
            self.tolerance,
            policy.buffer_tolerance,
        );
        let mut restarts = 0;
        let mut contractions = 0;
        let mut certificates = 0;
        let mut reduced_width = 0;
        let buffer = (0..self.bands).map(|j| j >= policy.strict_bands).collect();
        // Owned result storage is prepared before the step loop; scratch stays
        // with the workspace even when a caller retains the previous solution.
        let mut result_values = vec![0.; self.bands];
        let mut result_vectors = RealMatrix::zeros(n, self.bands)?;
        let mut result_residuals = vec![0.; self.bands];
        let mut result_converged = vec![false; self.bands];
        apply_into(op, basis, &mut applications, &mut work, products)?;
        let mut max_basis_dimension = 0;
        for iteration in 1..=self.max_iterations {
            max_basis_dimension = max_basis_dimension.max(basis.cols());
            let m = basis.cols();
            reduced_matrix_into(basis, products, reduced_width, gemm, reduced)?;
            reduced_width = m;
            coefficients.reshape_within_storage(m, m)?;
            eigensolver.solve_into(reduced, &mut values[..m], coefficients)?;
            let coefficients = coefficients.as_ref().subcols(0, self.bands);
            vectors.reshape_within_storage(n, self.bands)?;
            rotated_products.reshape_within_storage(n, self.bands)?;
            gemm.rotate_pair(
                basis.as_ref(),
                products.as_ref(),
                coefficients,
                vectors.as_mut(),
                rotated_products.as_mut(),
            )?;
            residual_into(
                vectors,
                rotated_products,
                &values[..self.bands],
                residual_block,
                residuals,
                &pool,
            )?;
            let mut certified = false;
            let mut trusted = true;
            if !ritz
                && (certify || iteration == self.max_iterations)
                && (meets_policy(residuals)
                    || (iteration == self.max_iterations && meets_strict(residuals)))
            {
                certificates += 1;
                apply_into(op, vectors, &mut applications, &mut work, actual)?;
                residual_into(
                    vectors,
                    actual,
                    &values[..self.bands],
                    residual_block,
                    residuals,
                    &pool,
                )?;
                certified = true;
                trusted = meets_policy(residuals);
            }
            let maximum = residuals.iter().copied().fold(0., f64::max);
            history.push(maximum);
            work.residuals(iteration, residuals);
            for (j, flag) in converged.iter_mut().enumerate() {
                *flag = if ritz {
                    iteration > 1 && (values[j] - previous[j]).abs() < threshold(j)
                } else {
                    residuals[j] <= threshold(j)
                };
            }
            previous.copy_from_slice(&values[..self.bands]);
            let finished = converged.iter().all(|&b| b)
                || (iteration == self.max_iterations
                    && converged[..policy.strict_bands].iter().all(|&b| b));
            if finished {
                work.returned();
                result_values.copy_from_slice(&values[..self.bands]);
                result_vectors
                    .as_mut_slice()
                    .copy_from_slice(vectors.as_slice());
                result_residuals.copy_from_slice(residuals);
                result_converged.copy_from_slice(converged);
                return Ok(RealEigenSolution {
                    eigenvalues: result_values,
                    eigenvectors: result_vectors,
                    residuals: result_residuals,
                    converged: result_converged,
                    buffer,
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
                result_residuals.copy_from_slice(residuals);
                return Err(SolverError::NotConverged {
                    iterations: iteration,
                    max_residual: maximum,
                    residuals: result_residuals,
                });
            }
            if ritz && iteration == 1 && basis.cols() == n {
                continue;
            }
            selected.clear();
            selected.extend(
                (0..self.bands)
                    .filter(|&j| !converged[j])
                    .take(self.bands.min(capacity - self.bands)),
            );
            corrections.reshape_within_storage(n, selected.len())?;
            pool.install(|| {
                column_work(corrections.as_mut_slice(), n, 0, &|k, dst| {
                    let j = selected[k];
                    // The default damping keeps the original tight loop; research
                    // variants take their own loops (no per-element dispatch).
                    match preconditioner {
                        crate::DavidsonPreconditioner::Default => {
                            for ((dst, &r), &d) in
                                dst.iter_mut().zip(column(residual_block, j)).zip(diagonal)
                            {
                                *dst = precondition(r, d, values[j])?;
                            }
                        }
                        crate::DavidsonPreconditioner::Tpa => {
                            let scale = crate::kinetic_scale(
                                column(vectors, j)
                                    .iter()
                                    .zip(diagonal)
                                    .map(|(&c, &t)| (c * c, t)),
                            )?;
                            for ((dst, &r), &d) in
                                dst.iter_mut().zip(column(residual_block, j)).zip(diagonal)
                            {
                                *dst = -r * crate::tpa_weight(d, scale);
                                finite(&[*dst])?;
                            }
                        }
                        crate::DavidsonPreconditioner::RitzShift => {
                            for ((dst, &r), &d) in
                                dst.iter_mut().zip(column(residual_block, j)).zip(diagonal)
                            {
                                let delta = d - values[j];
                                finite(&[delta])?;
                                *dst = -r / delta.max(0.1);
                                finite(&[*dst])?;
                            }
                        }
                    }
                    Ok(())
                })
            })?;
            if !trusted || basis.cols() + corrections.cols() > capacity {
                restarts += 1;
                work.restart();
                reduced_width = 0;
                if trusted && orthogonality_defect(vectors) <= 1e-12 {
                    contractions += 1;
                    std::mem::swap(basis, vectors);
                    std::mem::swap(products, rotated_products);
                } else {
                    projection.run_into(vectors, None, gemm, basis)?;
                    apply_into(op, basis, &mut applications, &mut work, products)?;
                }
            }
            projection.run_into(corrections, Some(basis), gemm, extra)?;
            if extra.cols() == 0 {
                random.reshape_within_storage(n, self.bands.min(capacity - basis.cols()))?;
                for value in random.as_mut_slice() {
                    *value = rng.next();
                }
                projection.run_into(random, Some(basis), gemm, extra)?;
            }
            if extra.cols() == 0 {
                result_residuals.copy_from_slice(residuals);
                return Err(SolverError::NotConverged {
                    iterations: iteration,
                    max_residual: maximum,
                    residuals: result_residuals,
                });
            }
            apply_into(op, extra, &mut applications, &mut work, extra_products)?;
            append_into(products, extra_products)?;
            append_into(basis, extra)?;
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
    crate::validate_bands(kinetic.len(), bands)?;
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
#[cfg(test)]
fn prefix(x: &RealMatrix, count: usize) -> Result<RealMatrix, LinalgError> {
    RealMatrix::new(x.rows(), count, x.as_slice()[..x.rows() * count].to_vec())
}
#[cfg(test)]
fn append(a: &RealMatrix, b: &RealMatrix) -> Result<RealMatrix, LinalgError> {
    let mut data = Vec::with_capacity(a.as_slice().len() + b.as_slice().len());
    data.extend_from_slice(a.as_slice());
    data.extend_from_slice(b.as_slice());
    RealMatrix::new(a.rows(), a.cols() + b.cols(), data)
}
fn apply_into(
    op: &dyn RealLinearOperator,
    x: &RealMatrix,
    count: &mut usize,
    work: &mut crate::work_accounting::Trace,
    y: &mut RealMatrix,
) -> Result<(), SolverError> {
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
    y.reshape_within_storage(x.rows(), x.cols())?;
    op.apply_block(x, y)?;
    if y.rows() != x.rows() || y.cols() != x.cols() {
        return Err(LinalgError::Dimension.into());
    }
    finite(y.as_slice())?;
    *count += x.cols();
    Ok(())
}
fn append_into(a: &mut RealMatrix, b: &RealMatrix) -> Result<(), LinalgError> {
    let old_len = a.as_slice().len();
    a.reshape_within_storage(a.rows(), a.cols() + b.cols())?;
    a.as_mut_slice()[old_len..].copy_from_slice(b.as_slice());
    Ok(())
}
#[cfg(test)]
fn reduced_matrix(
    v: &RealMatrix,
    hv: &RealMatrix,
    old: Option<&RealMatrix>,
    ws: &mut RealGemmWorkspace,
) -> Result<RealMatrix, LinalgError> {
    let mut out = RealMatrix::zeros(v.cols(), v.cols())?;
    let prefix = old.map_or(0, RealMatrix::cols);
    if let Some(old) = old {
        out.as_mut_slice()[..prefix * prefix].copy_from_slice(old.as_slice());
    }
    reduced_matrix_into(v, hv, prefix, ws, &mut out)?;
    Ok(out)
}
fn reduced_matrix_into(
    v: &RealMatrix,
    hv: &RealMatrix,
    prefix: usize,
    ws: &mut RealGemmWorkspace,
    out: &mut RealMatrix,
) -> Result<(), LinalgError> {
    let m = v.cols();
    out.reshape_within_storage(m, m)?;
    // Expand the contiguous leading dimension in place, last column first.
    // Keep the exact compact layout passed to the old sequential EVD.
    if prefix != m {
        for j in (0..prefix).rev() {
            out.as_mut_slice()
                .copy_within(j * prefix..(j + 1) * prefix, j * m);
        }
    }
    // One old/new rectangle and the new square, in one GEMM. Preserve the
    // old prefix exactly; the newly computed upper triangle defines symmetry.
    ws.multiply_into(
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
    Ok(())
}
// Fixed pieces of four whole columns, independent of the pool budget.
// Every scalar reduction remains inside its original column, in its old order.
fn column_work<F>(data: &mut [f64], n: usize, first: usize, f: &F) -> Result<(), LinalgError>
where
    F: Fn(usize, &mut [f64]) -> Result<(), LinalgError> + Sync,
{
    let count = data.len() / n;
    if count > 4 {
        let split = (count.div_ceil(4) / 2) * 4;
        let (left, right) = data.split_at_mut(split * n);
        let (a, b) = WorkerPool::join(
            || column_work(left, n, first, f),
            || column_work(right, n, first + split, f),
        );
        a?;
        b
    } else {
        for (j, dst) in data.chunks_exact_mut(n).enumerate() {
            f(first + j, dst)?;
        }
        Ok(())
    }
}
fn residual_columns(
    r: &mut [f64],
    v: &[f64],
    hv: &[f64],
    values: &[f64],
    norms: &mut [f64],
    n: usize,
) {
    if values.len() > 4 {
        let split = (values.len().div_ceil(4) / 2) * 4;
        let (rl, rr) = r.split_at_mut(split * n);
        let (nl, nr) = norms.split_at_mut(split);
        WorkerPool::join(
            || {
                residual_columns(
                    rl,
                    &v[..split * n],
                    &hv[..split * n],
                    &values[..split],
                    nl,
                    n,
                )
            },
            || {
                residual_columns(
                    rr,
                    &v[split * n..],
                    &hv[split * n..],
                    &values[split..],
                    nr,
                    n,
                )
            },
        );
    } else {
        for (j, &lambda) in values.iter().enumerate() {
            let r = &mut r[j * n..(j + 1) * n];
            for ((r, &h), &v) in r
                .iter_mut()
                .zip(&hv[j * n..(j + 1) * n])
                .zip(&v[j * n..(j + 1) * n])
            {
                *r = h - lambda * v;
            }
            norms[j] = norm(r);
        }
    }
}
fn residual_into(
    v: &RealMatrix,
    hv: &RealMatrix,
    values: &[f64],
    r: &mut RealMatrix,
    norms: &mut [f64],
    pool: &WorkerPool,
) -> Result<(), LinalgError> {
    pool.install(|| {
        residual_columns(
            r.as_mut_slice(),
            v.as_slice(),
            hv.as_slice(),
            values,
            norms,
            v.rows(),
        )
    });
    finite(norms)
}
fn orthogonality_defect(v: &RealMatrix) -> f64 {
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
    pool: std::sync::Arc<WorkerPool>,
    normalized: RealMatrix,
    overlap: RealMatrix,
    projected: RealMatrix,
    output: RealMatrix,
    column: Vec<f64>,
}
impl Projection {
    #[cfg(test)]
    fn new(n: usize, cap: usize) -> Result<Self, LinalgError> {
        Self::with_pool(n, cap, WorkerPool::new(1).unwrap())
    }
    fn with_pool(
        n: usize,
        cap: usize,
        pool: std::sync::Arc<WorkerPool>,
    ) -> Result<Self, LinalgError> {
        Ok(Self {
            pool,
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
        let mut out = RealMatrix::zeros(self.output.rows(), self.output.cols())?;
        self.run_into(x, against, ws, &mut out)?;
        Ok(out)
    }
    fn run_into(
        &mut self,
        x: &RealMatrix,
        against: Option<&RealMatrix>,
        ws: &mut RealGemmWorkspace,
        out: &mut RealMatrix,
    ) -> Result<(), LinalgError> {
        finite(x.as_slice())?;
        let n = x.rows();
        let mut count = x.cols();
        self.pool.install(|| {
            column_work(
                &mut self.normalized.as_mut_slice()[..n * count],
                n,
                0,
                &|j, dst| {
                    let len = norm(column(x, j));
                    finite(&[len])?;
                    for (dst, &src) in dst.iter_mut().zip(column(x, j)) {
                        *dst = if len == 0. { 0. } else { src / len };
                    }
                    Ok(())
                },
            )
        })?;
        let against = against.filter(|q| q.cols() > 0);
        for pass in 0..if against.is_some() { 2 } else { 1 } {
            if pass > 0 {
                // Both buffers have the same capacity; the previous MGS result
                // becomes input, and the old input storage becomes output.
                std::mem::swap(&mut self.normalized, &mut self.output);
            }
            if let Some(q) = against {
                ws.multiply_into(
                    q.as_ref(),
                    RealOp::Transpose,
                    self.normalized.as_ref().subcols(0, count),
                    self.overlap.as_mut().submatrix_mut(0, 0, q.cols(), count),
                )?;
                ws.multiply_into(
                    q.as_ref(),
                    RealOp::Identity,
                    self.overlap.as_ref().submatrix(0, 0, q.cols(), count),
                    self.projected.as_mut().subcols_mut(0, count),
                )?;
                self.pool.install(|| {
                    column_work(
                        &mut self.normalized.as_mut_slice()[..n * count],
                        n,
                        0,
                        &|j, dst| {
                            for (dst, &src) in dst.iter_mut().zip(column(&self.projected, j)) {
                                *dst -= src;
                            }
                            Ok(())
                        },
                    )
                })?;
            }
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
        let cap = self.output.cols();
        std::mem::swap(&mut self.output, out);
        out.reshape_within_storage(n, count)?;
        self.output.reshape_within_storage(n, cap)?;
        Ok(())
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

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn incremental_prefix_and_projection_rank() {
        let n = 257;
        let mut ws = RealGemmWorkspace::with_pool(n, 8, WorkerPool::new(1).unwrap()).unwrap();
        let mut projection = Projection::new(n, 8).unwrap();
        let raw = Random(17).block(n, 8).unwrap();
        let q = projection.run(&raw, None, &mut ws).unwrap();
        let mut hq = q.clone();
        for j in 0..8 {
            for i in 0..n {
                hq[(i, j)] *= i as f64 / n as f64;
            }
        }
        let old = reduced_matrix(
            &prefix(&q, 3).unwrap(),
            &prefix(&hq, 3).unwrap(),
            None,
            &mut ws,
        )
        .unwrap();
        let full = reduced_matrix(&q, &hq, None, &mut ws).unwrap();
        let incremental = reduced_matrix(&q, &hq, Some(&old), &mut ws).unwrap();
        for j in 0..8 {
            for i in 0..8 {
                assert!((full[(i, j)] - incremental[(i, j)]).abs() < 1e-14);
                assert_eq!(incremental[(i, j)].to_bits(), incremental[(j, i)].to_bits());
                if i < 3 && j < 3 {
                    assert_eq!(old[(i, j)].to_bits(), incremental[(i, j)].to_bits());
                }
            }
        }
        // Repeated and zero columns are dropped; normalization protects tiny
        // independent columns from an absolute-scale rank decision.
        let mut incoming = RealMatrix::zeros(n, 3).unwrap();
        for i in 0..n {
            incoming[(i, 0)] = q[(i, 0)];
            incoming[(i, 1)] = q[(i, 4)] * 1e-200;
        }
        let extra = projection
            .run(&incoming, Some(&prefix(&q, 3).unwrap()), &mut ws)
            .unwrap();
        assert_eq!(extra.cols(), 1);
        assert!(orthogonality_defect(&append(&prefix(&q, 3).unwrap(), &extra).unwrap()) < 1e-12);
    }
    #[test]
    fn sign_preserving_floor() {
        assert_eq!(precondition(1., 0., 0.).unwrap(), -10.);
        assert_eq!(precondition(1., 0., 0.01).unwrap(), 10.);
        assert_eq!(precondition(1., 0.01, 0.).unwrap(), -10.);
        assert_eq!(precondition(1., 0., 2.).unwrap(), 0.5);
        assert!(precondition(1., f64::MAX, -f64::MAX).is_err());
    }
}
