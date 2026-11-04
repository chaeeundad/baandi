//! Backend-independent dense and matrix-free Hermitian eigensolvers.
#[path = "mixing/kerker.rs"]
pub mod kerker;
pub mod mixing;
pub mod rmm;
#[doc(hidden)]
pub use dft_numerics::pool::analysis;
pub use rmm::{
    HybridPhase, HybridRmm, HybridSolution, RmmDiis, RmmError, RmmFailureReason, SolverChoice,
};

use dft_numerics::linalg::{
    ComplexMatrix, GemmWorkspace, LinalgError,
    MatrixOp::{Adjoint, Identity},
    ProjectionWorkspace, hermitian_eigh, norm, orthonormalize,
};
use num_complex::Complex64;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum SolverError {
    #[error(transparent)]
    Linalg(#[from] LinalgError),
    #[error("invalid solver configuration: {0}")]
    InvalidParameter(&'static str),
    #[error("operator/preconditioner failure: {0}")]
    Application(String),
    #[error("dense dimension {dimension} exceeds limit {limit}")]
    DenseLimit { dimension: usize, limit: usize },
    #[error("not converged after {iterations} iterations; max residual {max_residual:e}")]
    NotConverged {
        iterations: usize,
        max_residual: f64,
        residuals: Vec<f64>,
    },
}

/// Hermitian by contract. Implementations must overwrite all of y, without changing x.
/// Both slices have length dimension(); errors must be returned, not hidden.
pub trait LinearOperator {
    fn dimension(&self) -> usize;
    /// Available intra-k workers; implementations must not exceed the total budget.
    fn worker_pool(&self) -> Option<std::sync::Arc<dft_numerics::pool::WorkerPool>> {
        None
    }
    fn inner_threads(&self) -> usize {
        1
    }
    fn apply(&self, x: &[Complex64], y: &mut [Complex64]) -> Result<(), SolverError>;
    fn apply_block(&self, x: &ComplexMatrix, y: &mut ComplexMatrix) -> Result<(), SolverError> {
        if x.rows() != self.dimension() || y.rows() != x.rows() || y.cols() != x.cols() {
            return Err(LinalgError::Dimension.into());
        }
        for j in 0..x.cols() {
            self.apply(x.column(j), y.column_mut(j))?;
        }
        y.validate_finite()?;
        Ok(())
    }
}
/// Transform residual columns, using the matching Ritz values; overwrite output.
pub trait Preconditioner {
    fn apply(
        &self,
        eigenvalues: &[f64],
        residuals: &ComplexMatrix,
        output: &mut ComplexMatrix,
    ) -> Result<(), SolverError>;
}

/// Correction -r/(d-lambda), with sign-preserving denominator floor (zero takes +).
#[derive(Debug, Clone)]
pub struct DiagonalPreconditioner {
    diagonal: Vec<f64>,
    denominator_floor: f64,
}
impl DiagonalPreconditioner {
    /// The positive floor has the same units as the diagonal/eigenvalues.
    pub fn new(diagonal: Vec<f64>, denominator_floor: f64) -> Result<Self, SolverError> {
        if diagonal.iter().any(|d| !d.is_finite())
            || !denominator_floor.is_finite()
            || denominator_floor <= 0.0
        {
            return Err(SolverError::InvalidParameter(
                "finite diagonal and positive finite denominator floor required",
            ));
        }
        Ok(Self {
            diagonal,
            denominator_floor,
        })
    }
}
impl Preconditioner for DiagonalPreconditioner {
    fn apply(
        &self,
        eigenvalues: &[f64],
        residuals: &ComplexMatrix,
        output: &mut ComplexMatrix,
    ) -> Result<(), SolverError> {
        if residuals.rows() != self.diagonal.len()
            || eigenvalues.len() != residuals.cols()
            || output.rows() != residuals.rows()
            || output.cols() != residuals.cols()
        {
            return Err(LinalgError::Dimension.into());
        }
        residuals.validate_finite()?;
        if eigenvalues.iter().any(|x| !x.is_finite()) {
            return Err(LinalgError::NonFinite.into());
        }
        for (j, &lambda) in eigenvalues.iter().enumerate() {
            for (i, &d) in self.diagonal.iter().enumerate() {
                let delta = d - lambda;
                if !delta.is_finite() {
                    return Err(LinalgError::NonFinite.into());
                }
                let denom = if delta < 0.0 {
                    -delta.abs().max(self.denominator_floor)
                } else {
                    delta.max(self.denominator_floor)
                };
                output[(i, j)] = -residuals[(i, j)] / denom;
            }
        }
        output.validate_finite()?;
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct EigenSolution {
    pub eigenvalues: Vec<f64>,
    pub eigenvectors: ComplexMatrix,
    pub residuals: Vec<f64>,
    /// Whether each residual meets its band's requested tolerance.
    pub converged: Vec<bool>,
    /// Bands outside the strict prefix; excluded from strict certification.
    pub buffer: Vec<bool>,
    /// True when residuals were recomputed by applying the operator to Ritz vectors.
    pub residuals_certified: bool,
    pub iterations: usize,
    /// Maximum residual 2-norm at each Rayleigh-Ritz iteration. In recurrence
    /// mode this is the recurrence norm, replaced by the true norm whenever
    /// certification is attempted (including failed certification).
    pub history: Vec<f64>,
    /// Number of vector columns passed to the operator (including certification).
    pub operator_applications: usize,
}

/// Oracle only: rejects dimensions above 1024 before allocating or applying H.
#[derive(Debug, Default, Clone, Copy)]
pub struct DenseEigensolver;
impl DenseEigensolver {
    pub const MAX_DIMENSION: usize = 1024;
    pub fn solve(
        &self,
        operator: &dyn LinearOperator,
        bands: usize,
    ) -> Result<EigenSolution, SolverError> {
        let n = operator.dimension();
        if n > Self::MAX_DIMENSION {
            return Err(SolverError::DenseLimit {
                dimension: n,
                limit: Self::MAX_DIMENSION,
            });
        }
        validate_bands(n, bands)?;
        let mut unit = ComplexMatrix::zeros(n, n)?;
        for i in 0..n {
            unit[(i, i)] = Complex64::new(1.0, 0.0);
        }
        let matrix = apply(operator, &unit)?;
        let (mut values, vectors) = hermitian_eigh(&matrix)?;
        values.truncate(bands);
        let vectors = first_columns(&vectors, bands)?;
        let (_, residuals) = residual(operator, &vectors, &values)?;
        let maximum = max_residual(&residuals);
        Ok(EigenSolution {
            eigenvalues: values,
            eigenvectors: vectors,
            residuals,
            converged: vec![true; bands],
            buffer: vec![false; bands],
            residuals_certified: true,
            iterations: 1,
            history: vec![maximum],
            operator_applications: n + bands,
        })
    }
}

/// Additive accuracy policy. Residual norms use the operator's eigenvalue units.
#[derive(Debug, Clone, Copy)]
pub struct BandTolerances {
    pub strict_bands: usize,
    pub buffer_tolerance: f64,
}

#[derive(Clone, Copy)]
struct SolveMode {
    recurrence: bool,
    contract_restart: bool,
    certify: bool,
}

#[derive(Debug, Clone)]
pub struct Davidson {
    pub bands: usize,
    /// Absolute residual norm, in the operator's eigenvalue units.
    pub tolerance: f64,
    pub max_iterations: usize,
    pub max_subspace: usize,
    /// Maximum number of residual corrections appended per iteration.
    pub block_size: usize,
    pub seed: u64,
}
impl Davidson {
    pub fn solve(
        &self,
        operator: &dyn LinearOperator,
        initial: Option<&ComplexMatrix>,
        preconditioner: Option<&dyn Preconditioner>,
    ) -> Result<EigenSolution, SolverError> {
        self.solve_with_recurrence(operator, initial, preconditioner, true)
    }

    /// Select stored-product residuals, always certified with a true application.
    /// `false` retains the bit-identical cached-basis algorithm for operators
    /// with column-independent block rounding; `solve` uses `true`.
    pub fn solve_with_recurrence(
        &self,
        operator: &dyn LinearOperator,
        initial: Option<&ComplexMatrix>,
        preconditioner: Option<&dyn Preconditioner>,
        recurrence: bool,
    ) -> Result<EigenSolution, SolverError> {
        self.solve_with_restart(operator, initial, preconditioner, recurrence, recurrence)
    }

    /// Retain the pre-contraction path as an arithmetic oracle.
    pub fn solve_with_restart(
        &self,
        operator: &dyn LinearOperator,
        initial: Option<&ComplexMatrix>,
        preconditioner: Option<&dyn Preconditioner>,
        recurrence: bool,
        contract_restart: bool,
    ) -> Result<EigenSolution, SolverError> {
        self.solve_policy(
            operator,
            initial,
            preconditioner,
            SolveMode {
                recurrence,
                contract_restart,
                certify: true,
            },
            BandTolerances {
                strict_bands: self.bands,
                buffer_tolerance: self.tolerance,
            },
        )
    }

    /// Strict prefix plus optional buffer. At the iteration limit, success still
    /// requires a fresh true-residual certificate for every strict band.
    pub fn solve_with_band_tolerances(
        &self,
        operator: &dyn LinearOperator,
        initial: Option<&ComplexMatrix>,
        preconditioner: Option<&dyn Preconditioner>,
        policy: BandTolerances,
    ) -> Result<EigenSolution, SolverError> {
        self.solve_with_band_tolerances_and_certification(
            operator,
            initial,
            preconditioner,
            policy,
            true,
        )
    }

    /// With `certify=false`, ordinary successful returns may report recurrence
    /// residuals. Intended only for intermediate SCF maps; final solves must
    /// certify. An iteration-limit return always checks true residuals.
    pub fn solve_with_band_tolerances_and_certification(
        &self,
        operator: &dyn LinearOperator,
        initial: Option<&ComplexMatrix>,
        preconditioner: Option<&dyn Preconditioner>,
        policy: BandTolerances,
        certify: bool,
    ) -> Result<EigenSolution, SolverError> {
        self.solve_policy(
            operator,
            initial,
            preconditioner,
            SolveMode {
                recurrence: true,
                contract_restart: true,
                certify,
            },
            policy,
        )
    }

    fn solve_policy(
        &self,
        operator: &dyn LinearOperator,
        initial: Option<&ComplexMatrix>,
        preconditioner: Option<&dyn Preconditioner>,
        mode: SolveMode,
        policy: BandTolerances,
    ) -> Result<EigenSolution, SolverError> {
        let _phase = dft_numerics::pool::analysis::Span::new("davidson");
        let SolveMode {
            recurrence,
            contract_restart,
            certify,
        } = mode;
        let n = operator.dimension();
        validate_bands(n, self.bands)?;
        if policy.strict_bands == 0
            || policy.strict_bands > self.bands
            || !policy.buffer_tolerance.is_finite()
            || policy.buffer_tolerance < self.tolerance
        {
            return Err(SolverError::InvalidParameter(
                "invalid strict prefix or buffer tolerance",
            ));
        }
        let threshold = |j: usize| {
            if j < policy.strict_bands {
                self.tolerance
            } else {
                policy.buffer_tolerance
            }
        };
        let meets_policy = |rs: &[f64]| rs.iter().enumerate().all(|(j, &r)| r <= threshold(j));
        let meets_strict = |rs: &[f64]| {
            rs[..policy.strict_bands]
                .iter()
                .all(|&r| r <= self.tolerance)
        };
        let capacity = self.max_subspace.min(n);
        if !self.tolerance.is_finite()
            || self.tolerance <= 0.0
            || self.max_iterations == 0
            || self.block_size == 0
            || capacity < self.bands
            || (capacity == self.bands && self.bands < n)
        {
            return Err(SolverError::InvalidParameter(
                "positive tolerance/iterations/block size and subspace larger than bands required (unless bands=n)",
            ));
        }
        let mut gemm_work = match operator.worker_pool() {
            Some(pool) => GemmWorkspace::with_pool(n, capacity, operator.inner_threads(), pool)?,
            None => GemmWorkspace::new(n, capacity, operator.inner_threads())?,
        };
        let mut projection = ProjectionWorkspace::new(n, capacity)?;
        let mut rng = SplitMix64(self.seed);
        let mut basis = if let Some(initial) = initial {
            if initial.rows() != n || initial.cols() > capacity {
                return Err(LinalgError::Dimension.into());
            }
            projection
                .orthonormalize(initial, None, &mut gemm_work)?
                .vectors
        } else {
            ComplexMatrix::zeros(n, 0)?
        };
        basis.reserve_columns(capacity)?;
        let initial_width = self.bands.max(self.block_size).min(capacity);
        // Fill missing rank with independent seeded columns, including deficient user guesses.
        while basis.cols() < initial_width {
            let random = rng.block(n, initial_width - basis.cols())?;
            let extra = projection
                .orthonormalize(&random, Some(&basis), &mut gemm_work)?
                .vectors;
            if extra.cols() == 0 {
                return Err(SolverError::InvalidParameter(
                    "unable to complete initial basis",
                ));
            }
            basis.append_columns(&extra)?;
        }
        let mut history = Vec::new();
        let mut operator_applications = 0;
        // Fixed-capacity cache; only the first basis.cols() columns are live.
        let mut cached = ComplexMatrix::zeros(n, 0)?;
        cached.reserve_columns(capacity)?;
        let mut cache_valid = false;
        let mut restart_cache = ComplexMatrix::zeros(n, 0)?;
        restart_cache.reserve_columns(capacity)?;
        // The operator API still owns a contiguous output block. Reuse one
        // narrow buffer for new columns; no full-subspace packing is needed.
        let mut appended_products = ComplexMatrix::zeros(n, 0)?;
        appended_products.reserve_columns(self.block_size.min(capacity - self.bands))?;
        for iteration in 1..=self.max_iterations {
            let _iteration = dft_numerics::pool::analysis::Span::new("iteration");
            cached.resize_columns(basis.cols())?;
            if !cache_valid {
                operator_applications += basis.cols();
                operator.apply_block(&basis, &mut cached)?;
                if cached.rows() != n || cached.cols() != basis.cols() {
                    return Err(LinalgError::Dimension.into());
                }
                cached.validate_finite()?;
                // Third-party operators may replace the output allocation.
                // The Hamiltonian preserves it, so this is normally a no-op.
                cached.reserve_columns(capacity)?;
            }
            // Borrow the live prefix directly; the backend sees its original
            // column stride and no first_columns allocation or copy.
            let h_basis = &cached;
            let projected = gemm_work.multiply(&basis, Adjoint, h_basis, Identity)?;
            let (mut values, mut coefficients) = {
                let _phase = dft_numerics::pool::analysis::Span::serial("evd");
                hermitian_eigh(&projected)?
            };
            values.truncate(self.bands);
            coefficients.resize_columns(self.bands)?;
            let (vectors, products) = if recurrence {
                let (vectors, products) = gemm_work.rotate_pair(&basis, h_basis, &coefficients)?;
                (vectors, Some(products))
            } else {
                (
                    gemm_work.multiply(&basis, Identity, &coefficients, Identity)?,
                    None,
                )
            };
            // A contracted restart repeats exactly this H-basis rotation.
            // Preserve its result before residual formation only when the
            // current basis can overflow; no extra allocation or arithmetic.
            let mut restart_products_saved = false;
            let (mut residual_block, mut residuals) = if let Some(products) = products {
                if contract_restart
                    && basis.cols() + self.block_size.min(capacity - self.bands) > capacity
                {
                    restart_cache.replace_columns(&products)?;
                    restart_products_saved = true;
                }
                residual_from_products(products, &vectors, &values)?
            } else {
                operator_applications += self.bands;
                residual(operator, &vectors, &values)?
            };
            let mut products_trusted = true;
            let mut residuals_certified = !recurrence;
            if recurrence
                && (certify || iteration == self.max_iterations)
                && (meets_policy(&residuals)
                    || (iteration == self.max_iterations && meets_strict(&residuals)))
            {
                operator_applications += self.bands;
                (residual_block, residuals) = residual(operator, &vectors, &values)?;
                residuals_certified = true;
                products_trusted = meets_policy(&residuals);
                // A failed certificate supplies actual correction directions and
                // its actual maximum is recorded in this iteration's history.
            }
            let maximum = max_residual(&residuals);
            history.push(maximum);
            if meets_policy(&residuals)
                || (iteration == self.max_iterations && meets_strict(&residuals))
            {
                return Ok(EigenSolution {
                    eigenvalues: values,
                    eigenvectors: vectors,
                    converged: residuals
                        .iter()
                        .enumerate()
                        .map(|(j, &r)| r <= threshold(j))
                        .collect(),
                    residuals_certified,
                    buffer: (0..self.bands).map(|j| j >= policy.strict_bands).collect(),
                    residuals,
                    iterations: iteration,
                    history,
                    operator_applications,
                });
            }
            if iteration == self.max_iterations {
                return Err(SolverError::NotConverged {
                    iterations: iteration,
                    max_residual: maximum,
                    residuals,
                });
            }
            let selected: Vec<_> = residuals
                .iter()
                .enumerate()
                .filter(|(j, r)| **r > threshold(*j))
                .map(|(j, _)| j)
                .take(self.block_size.min(capacity - self.bands))
                .collect();
            let mut corrections = ComplexMatrix::zeros(n, selected.len())?;
            let mut shifts = Vec::new();
            dft_numerics::pool::analysis::copied(n * selected.len() * size_of::<Complex64>());
            for (k, &j) in selected.iter().enumerate() {
                corrections
                    .column_mut(k)
                    .copy_from_slice(residual_block.column(j));
                shifts.push(values[j]);
            }
            if let Some(p) = preconditioner {
                let mut output = ComplexMatrix::zeros(n, corrections.cols())?;
                {
                    let _phase = dft_numerics::pool::analysis::Span::serial("preconditioner");
                    p.apply(&shifts, &corrections, &mut output)?;
                }
                if output.rows() != n || output.cols() != corrections.cols() {
                    return Err(LinalgError::Dimension.into());
                }
                output.validate_finite()?;
                corrections = output;
            }
            cache_valid = products_trusted && basis.cols() + corrections.cols() <= capacity;
            if !cache_valid {
                let (restarted, valid) = if restart_products_saved
                    && products_trusted
                    && orthogonality_defect(&vectors) <= 1e-12
                {
                    (vectors, true)
                } else {
                    restart(
                        vectors,
                        h_basis,
                        &coefficients,
                        &mut restart_cache,
                        contract_restart && products_trusted,
                        &mut gemm_work,
                    )?
                };
                basis.replace_columns(&restarted)?;
                cache_valid = valid;
                std::mem::swap(&mut cached, &mut restart_cache);
            }
            let mut extra = projection
                .orthonormalize(&corrections, Some(&basis), &mut gemm_work)?
                .vectors;
            // A stalled direction can occur with an exact preconditioner or an invariant
            // initial space. Try fresh independent directions, never infer convergence.
            if extra.cols() == 0 {
                let width = self.block_size.min(capacity - basis.cols());
                let random = rng.block(n, width)?;
                extra = projection
                    .orthonormalize(&random, Some(&basis), &mut gemm_work)?
                    .vectors;
            }
            if extra.cols() == 0 {
                return Err(SolverError::NotConverged {
                    iterations: iteration,
                    max_residual: maximum,
                    residuals,
                });
            }
            if cache_valid {
                appended_products.resize_columns(extra.cols())?;
                let products = &mut appended_products;
                operator.apply_block(&extra, products)?;
                if products.rows() != n || products.cols() != extra.cols() {
                    return Err(LinalgError::Dimension.into());
                }
                products.validate_finite()?;
                products.reserve_columns(self.block_size.min(capacity - self.bands))?;
                operator_applications += extra.cols();
                dft_numerics::pool::analysis::copied(
                    extra.rows() * extra.cols() * size_of::<Complex64>(),
                );
                cached.resize_columns(basis.cols() + extra.cols())?;
                for j in 0..extra.cols() {
                    cached
                        .column_mut(basis.cols() + j)
                        .copy_from_slice(products.column(j));
                }
            }
            basis.append_columns(&extra)?;
        }
        unreachable!("positive iteration limit validated above")
    }
}
fn validate_bands(n: usize, bands: usize) -> Result<(), SolverError> {
    if bands == 0 || bands > n {
        Err(SolverError::InvalidParameter(
            "require 1 <= bands <= dimension",
        ))
    } else {
        Ok(())
    }
}
fn apply(op: &dyn LinearOperator, x: &ComplexMatrix) -> Result<ComplexMatrix, SolverError> {
    let mut y = ComplexMatrix::zeros(op.dimension(), x.cols())?;
    op.apply_block(x, &mut y)?;
    if y.rows() != op.dimension() || y.cols() != x.cols() {
        return Err(LinalgError::Dimension.into());
    }
    y.validate_finite()?;
    Ok(y)
}
fn residual(
    op: &dyn LinearOperator,
    vectors: &ComplexMatrix,
    values: &[f64],
) -> Result<(ComplexMatrix, Vec<f64>), SolverError> {
    residual_from_products(apply(op, vectors)?, vectors, values)
}
fn residual_from_products(
    mut r: ComplexMatrix,
    vectors: &ComplexMatrix,
    values: &[f64],
) -> Result<(ComplexMatrix, Vec<f64>), SolverError> {
    let _phase = dft_numerics::pool::analysis::Span::serial("residual_norms");
    let mut norms = Vec::new();
    for (j, &lambda) in values.iter().enumerate() {
        for i in 0..r.rows() {
            r[(i, j)] -= vectors[(i, j)] * lambda;
        }
        let magnitude = norm(r.column(j));
        if !magnitude.is_finite() {
            return Err(LinalgError::NonFinite.into());
        }
        norms.push(magnitude);
    }
    Ok((r, norms))
}
fn max_residual(r: &[f64]) -> f64 {
    r.iter().copied().fold(0.0, f64::max)
}
fn first_columns(a: &ComplexMatrix, count: usize) -> Result<ComplexMatrix, LinalgError> {
    let _phase = dft_numerics::pool::analysis::Span::serial("first_columns");
    dft_numerics::pool::analysis::copied(a.rows() * count * size_of::<Complex64>());
    ComplexMatrix::new(a.rows(), count, a.as_slice()[..a.rows() * count].to_vec())
}
#[cfg(test)]
fn concatenate(a: &ComplexMatrix, b: &ComplexMatrix) -> Result<ComplexMatrix, LinalgError> {
    let _phase = dft_numerics::pool::analysis::Span::serial("concatenate");
    if a.rows() != b.rows() {
        return Err(LinalgError::Dimension);
    }
    dft_numerics::pool::analysis::copied(
        (a.as_slice().len() + b.as_slice().len()) * size_of::<Complex64>(),
    );
    let mut data = Vec::with_capacity(a.as_slice().len() + b.as_slice().len());
    data.extend_from_slice(a.as_slice());
    data.extend_from_slice(b.as_slice());
    ComplexMatrix::new(a.rows(), a.cols() + b.cols(), data)
}
/// SplitMix64, deterministic and local: no process-global random state.
struct SplitMix64(u64);
impl SplitMix64 {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        ((z ^ (z >> 31)) >> 11) as f64 * (1.0 / ((1_u64 << 53) as f64)) - 0.5
    }
    fn block(&mut self, rows: usize, cols: usize) -> Result<ComplexMatrix, LinalgError> {
        let mut block = ComplexMatrix::zeros(rows, cols)?;
        for j in 0..cols {
            for i in 0..rows {
                block[(i, j)] = Complex64::new(self.next(), self.next());
            }
        }
        Ok(block)
    }
}
#[cfg(test)]
mod tests;

#[cfg(test)]
mod legacy;

// Maximum entrywise defect of V^H V from identity, without a Gram allocation.
fn orthogonality_defect(v: &ComplexMatrix) -> f64 {
    let mut maximum: f64 = 0.0;
    for j in 0..v.cols() {
        for k in 0..=j {
            let dot: Complex64 = v
                .column(k)
                .iter()
                .zip(v.column(j))
                .map(|(a, b)| a.conj() * b)
                .sum();
            let defect = (dot - if j == k { 1.0 } else { 0.0 }).norm();
            if !defect.is_finite() {
                return f64::INFINITY;
            }
            maximum = maximum.max(defect);
        }
    }
    maximum
}

// A valid cache means products may be reused and only appended columns need H.
fn restart(
    vectors: ComplexMatrix,
    h_basis: &ComplexMatrix,
    coefficients: &ComplexMatrix,
    cached: &mut ComplexMatrix,
    contract: bool,
    gemm_work: &mut GemmWorkspace,
) -> Result<(ComplexMatrix, bool), SolverError> {
    if contract && orthogonality_defect(&vectors) <= 1e-12 {
        // Derived from study note section 1: H(V C) = (H V) C.
        cached.resize_columns(vectors.cols())?;
        gemm_work.multiply_into(h_basis, Identity, coefficients, Identity, cached)?;
        Ok((vectors, true))
    } else {
        Ok((orthonormalize(&vectors, None)?.vectors, false))
    }
}

#[cfg(test)]
use dft_numerics::linalg::gemm;
