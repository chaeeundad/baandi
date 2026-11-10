//! Backend-independent dense and matrix-free Hermitian eigensolvers.
mod bordered;
pub mod davidson_stages;
pub mod distributed;
pub mod distributed_davidson;
pub mod distributed_fft;
#[path = "mixing/kerker.rs"]
pub mod kerker;
pub mod locking;
pub mod mixing;
pub mod real_davidson;
pub use real_davidson::{RealDavidson, RealEigenSolution};
#[doc(hidden)]
pub mod replay;
pub mod rmm;
mod spectral;
pub use spectral::SpectralCompleteness;
#[doc(hidden)]
pub mod work_accounting;
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
    /// Davidson supplies current Ritz vectors and selected column indices.
    /// Other preconditioners retain their original residual-only contract.
    fn apply_with_vectors(
        &self,
        eigenvalues: &[f64],
        residuals: &ComplexMatrix,
        _vectors: &ComplexMatrix,
        _selected: &[usize],
        output: &mut ComplexMatrix,
    ) -> Result<(), SolverError> {
        self.apply(eigenvalues, residuals, output)
    }
    fn apply(
        &self,
        eigenvalues: &[f64],
        residuals: &ComplexMatrix,
        output: &mut ComplexMatrix,
    ) -> Result<(), SolverError>;
}

/// Real symmetric operator in orthonormal real coordinates (ORC).
/// Overwrite every output entry; blocks supplied by RealDavidson have unit-norm
/// columns. A physics adapter may pair adjacent columns for FFT application.
pub trait RealLinearOperator {
    fn dimension(&self) -> usize;
    fn worker_pool(&self) -> Option<std::sync::Arc<dft_numerics::pool::WorkerPool>> {
        None
    }
    fn inner_threads(&self) -> usize {
        1
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), SolverError>;
    fn apply_block(
        &self,
        x: &dft_numerics::real::RealMatrix,
        y: &mut dft_numerics::real::RealMatrix,
    ) -> Result<(), SolverError> {
        let n = self.dimension();
        if x.rows() != n || y.rows() != n || x.cols() != y.cols() {
            return Err(LinalgError::Dimension.into());
        }
        for j in 0..x.cols() {
            self.apply(
                &x.as_slice()[j * n..(j + 1) * n],
                &mut y.as_mut_slice()[j * n..(j + 1) * n],
            )?;
        }
        if y.as_slice().iter().any(|x| !x.is_finite()) {
            return Err(LinalgError::NonFinite.into());
        }
        Ok(())
    }
}

/// Switchable kinetic Davidson damping; legacy entry points retain Default.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum DavidsonPreconditioner {
    #[default]
    Default,
    Tpa,
    RitzShift,
}
impl DavidsonPreconditioner {
    pub fn name(self) -> &'static str {
        match self {
            Self::Default => "default",
            Self::Tpa => "tpa",
            Self::RitzShift => "ritz_shift",
        }
    }
}

/// Ordered expectation of kinetic energy, with the existing 0.1 Ha floor.
fn kinetic_scale(weights: impl Iterator<Item = (f64, f64)>) -> Result<f64, LinalgError> {
    let (mut norm, mut mean) = (0.0, 0.0);
    for (weight, kinetic) in weights {
        if !weight.is_finite() || weight < 0.0 || !kinetic.is_finite() || kinetic < 0.0 {
            return Err(LinalgError::NonFinite);
        }
        norm += weight;
        mean += weight * kinetic;
    }
    let scale = mean / norm;
    if !scale.is_finite() {
        return Err(LinalgError::NonFinite);
    }
    Ok(scale.max(0.1))
}

/// Positive TPA damping. Reciprocal evaluation avoids polynomial overflow.
fn tpa_weight(kinetic: f64, scale: f64) -> f64 {
    let x = kinetic / scale;
    if x <= 1.0 {
        let p = 27.0 + x * (18.0 + x * (12.0 + x * 8.0));
        p / (p + 16.0 * x * x * x * x)
    } else {
        let y = scale / kinetic;
        let p = y * (8.0 + y * (12.0 + y * (18.0 + y * 27.0)));
        p / (p + 16.0)
    }
}

/// Correction -r/(d-lambda), with sign-preserving denominator floor (zero takes +).
#[derive(Debug, Clone)]
pub struct DiagonalPreconditioner {
    diagonal: Vec<f64>,
    denominator_floor: f64,
    variant: DavidsonPreconditioner,
}
impl DiagonalPreconditioner {
    /// Select a research variant; its constants are fixed independently of the legacy floor.
    pub fn with_variant(mut self, variant: DavidsonPreconditioner) -> Self {
        self.variant = variant;
        self
    }

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
            variant: DavidsonPreconditioner::Default,
        })
    }
}
impl Preconditioner for DiagonalPreconditioner {
    fn apply_with_vectors(
        &self,
        eigenvalues: &[f64],
        residuals: &ComplexMatrix,
        vectors: &ComplexMatrix,
        selected: &[usize],
        output: &mut ComplexMatrix,
    ) -> Result<(), SolverError> {
        if self.variant != DavidsonPreconditioner::Tpa {
            return self.apply(eigenvalues, residuals, output);
        }
        if residuals.rows() != self.diagonal.len()
            || vectors.rows() != residuals.rows()
            || selected.len() != residuals.cols()
            || eigenvalues.len() != residuals.cols()
            || output.rows() != residuals.rows()
            || output.cols() != residuals.cols()
            || selected.iter().any(|&j| j >= vectors.cols())
        {
            return Err(LinalgError::Dimension.into());
        }
        residuals.validate_finite()?;
        if eigenvalues.iter().any(|x| !x.is_finite()) {
            return Err(LinalgError::NonFinite.into());
        }
        for (k, &j) in selected.iter().enumerate() {
            let scale = kinetic_scale(
                vectors
                    .column(j)
                    .iter()
                    .zip(&self.diagonal)
                    .map(|(c, &t)| (c.norm_sqr(), t)),
            )?;
            for (i, &t) in self.diagonal.iter().enumerate() {
                output[(i, k)] = -residuals[(i, k)] * tpa_weight(t, scale);
            }
        }
        output.validate_finite()?;
        Ok(())
    }

    fn apply(
        &self,
        eigenvalues: &[f64],
        residuals: &ComplexMatrix,
        output: &mut ComplexMatrix,
    ) -> Result<(), SolverError> {
        if self.variant == DavidsonPreconditioner::Tpa {
            return Err(SolverError::InvalidParameter(
                "TPA requires current Ritz vectors",
            ));
        }
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
                let denom = if self.variant == DavidsonPreconditioner::RitzShift {
                    delta.max(0.1)
                } else if delta < 0.0 {
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
    pub locking: Option<locking::LockingReport>,
    pub spectral_completeness: Option<SpectralCompleteness>,
    pub eigenvalues: Vec<f64>,
    pub eigenvectors: ComplexMatrix,
    pub residuals: Vec<f64>,
    /// Whether each root meets its requested stopping criterion: residual norm
    /// in strict solves, consecutive Ritz-value change in standard solves.
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
            locking: None,
            spectral_completeness: None,
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
    ritz: bool,
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
                ritz: false,
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
                ritz: false,
                recurrence: true,
                contract_restart: true,
                certify,
            },
            policy,
        )
    }

    /// Standard-profile Ritz-value-change stopping in eigenvalue units.
    /// Empty roots use the buffer threshold, but the required prefix (occupied
    /// plus first empty) must converge even at the iteration limit. No locking:
    /// every rotated root is tested again. This is not a residual certificate.
    pub fn solve_ritz_changes(
        &self,
        operator: &dyn LinearOperator,
        initial: Option<&ComplexMatrix>,
        preconditioner: Option<&dyn Preconditioner>,
        policy: BandTolerances,
    ) -> Result<EigenSolution, SolverError> {
        self.solve_policy(
            operator,
            initial,
            preconditioner,
            SolveMode {
                ritz: true,
                recurrence: true,
                contract_restart: true,
                certify: false,
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
        // Diagnostic SCF oracle, deliberately not a production solver option.
        if option_env!("P15C_DENSE_SCF").is_some() {
            let solution = DenseEigensolver.solve(operator, self.bands)?;
            if replay::ENABLED {
                let mut trace = replay::Trace::new(self.bands, 1);
                trace.step(
                    operator.dimension(),
                    solution.operator_applications,
                    &solution.eigenvalues,
                    &solution.residuals,
                );
                replay::store(replay::Record {
                    dense_values: Some(solution.eigenvalues.clone()),
                    traces: vec![(
                        self.tolerance,
                        policy.strict_bands,
                        policy.buffer_tolerance,
                        true,
                        trace,
                    )],
                });
            }
            return Ok(solution);
        }
        let enriched;
        let initial = if option_env!("P15C_RICH_START").is_some()
            && initial.is_some_and(|v| {
                (0..v.cols())
                    .all(|j| v.column(j).iter().filter(|z| z.norm_sqr() > 0.0).count() == 1)
            }) {
            let mut start = initial.unwrap().clone();
            let width = (2 * self.bands)
                .min(self.max_subspace)
                .min(operator.dimension());
            start.reserve_columns(width)?;
            let random = SplitMix64(self.seed ^ 0x9e3779b97f4a7c15)
                .block(operator.dimension(), width.saturating_sub(start.cols()))?;
            start.append_columns(&random)?;
            enriched = start;
            Some(&enriched)
        } else {
            initial
        };
        if !replay::ENABLED {
            return self.solve_policy_traced::<false>(
                operator,
                initial,
                preconditioner,
                mode,
                policy,
                &mut replay::Trace::default(),
            );
        }
        let mut trace = replay::Trace::new(self.bands, self.max_iterations);
        let result = self.solve_policy_traced::<true>(
            operator,
            initial,
            preconditioner,
            mode,
            policy,
            &mut trace,
        );
        let mut record = replay::Record {
            dense_values: None,
            traces: vec![(
                self.tolerance,
                policy.strict_bands,
                policy.buffer_tolerance,
                result.is_ok(),
                trace,
            )],
        };
        if replay::SPECTRAL {
            record.dense_values = replay::dense_oracle(operator, self.bands)?;
            if option_env!("P15C_CANDIDATES").is_some() && record.dense_values.is_some() {
                for rich in [true, false] {
                    let solver = Self {
                        tolerance: 1e-10,
                        seed: self.seed ^ 0x9e3779b97f4a7c15,
                        ..self.clone()
                    };
                    let mut rng = SplitMix64(solver.seed);
                    let width = (2 * self.bands)
                        .min(self.max_subspace)
                        .min(operator.dimension());
                    let mut start = if rich {
                        initial
                            .cloned()
                            .unwrap_or(ComplexMatrix::zeros(operator.dimension(), 0)?)
                    } else {
                        ComplexMatrix::zeros(operator.dimension(), 0)?
                    };
                    start.reserve_columns(width)?;
                    start.append_columns(
                        &rng.block(operator.dimension(), width.saturating_sub(start.cols()))?,
                    )?;
                    let mut trace = replay::Trace::new(self.bands, self.max_iterations);
                    let checked = solver.solve_policy_traced::<true>(
                        operator,
                        Some(&start),
                        preconditioner,
                        SolveMode {
                            ritz: false,
                            certify: true,
                            ..mode
                        },
                        BandTolerances {
                            strict_bands: self.bands,
                            buffer_tolerance: 1e-10,
                        },
                        &mut trace,
                    );
                    record
                        .traces
                        .push((1e-10, self.bands, 1e-10, checked.is_ok(), trace));
                }
            }
            replay::store(record);
            return result;
        }
        // Independent solves from identical starts, with the same expansion,
        // restart and preconditioner. The last solve is an all-root reference.
        for tolerance in [1e-2_f64, 1e-3, 1e-4, 1e-5, 1e-6, 1e-10] {
            let solver = Self {
                tolerance,
                ..self.clone()
            };
            let strict = if tolerance == 1e-10 {
                self.bands
            } else {
                policy.strict_bands
            };
            let buffer = if strict == self.bands {
                tolerance
            } else {
                tolerance.max(1e-5)
            };
            let mut trace = replay::Trace::new(self.bands, self.max_iterations);
            let replayed = solver.solve_policy_traced::<true>(
                operator,
                initial,
                preconditioner,
                SolveMode {
                    ritz: false,
                    certify: tolerance == 1e-10 || mode.certify,
                    ..mode
                },
                BandTolerances {
                    strict_bands: strict,
                    buffer_tolerance: buffer,
                },
                &mut trace,
            );
            record
                .traces
                .push((tolerance, strict, buffer, replayed.is_ok(), trace));
        }
        // An independent full-space oracle for the small synthetic fixtures.
        // This diagnostic work is not charged to any replay trajectory.
        record.dense_values = replay::dense_oracle(operator, self.bands)?;
        replay::store(record);
        result
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_policy_traced<const TRACE: bool>(
        &self,
        operator: &dyn LinearOperator,
        initial: Option<&ComplexMatrix>,
        preconditioner: Option<&dyn Preconditioner>,
        mode: SolveMode,
        policy: BandTolerances,
        trace: &mut replay::Trace,
    ) -> Result<EigenSolution, SolverError> {
        let _phase = dft_numerics::pool::analysis::Span::new("davidson");
        let SolveMode {
            ritz,
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
            if j < policy.strict_bands - usize::from(ritz) {
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
        let mut work = work_accounting::Trace::new(
            self.bands,
            policy.strict_bands,
            self.tolerance,
            policy.buffer_tolerance,
        );
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
        let mut previous_values = if ritz {
            vec![f64::NAN; self.bands]
        } else {
            Vec::new()
        };
        let mut changes = if ritz {
            vec![f64::INFINITY; self.bands]
        } else {
            Vec::new()
        };
        let mut reduced_cache = None;
        for iteration in 1..=self.max_iterations {
            gemm_work.accounting_step(iteration);
            let _iteration = dft_numerics::pool::analysis::Span::new("iteration");
            cached.resize_columns(basis.cols())?;
            if !cache_valid {
                work.apply(
                    basis.cols(),
                    (iteration == 1 && basis.cols() == self.bands)
                        .then_some(self.bands - policy.strict_bands),
                );
                operator_applications += basis.cols();
                if TRACE {
                    trace.blocks += 1;
                }
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
            gemm_work.accounting_site("reduced_hamiltonian");
            let projected = if option_env!("P24A_RECTANGLE_ONLY").is_some() {
                gemm_work.reduced_hamiltonian_hermitian(&basis, h_basis, reduced_cache.as_ref())?
            } else {
                gemm_work.reduced_hamiltonian(&basis, h_basis, reduced_cache.as_ref())?
            };
            let (mut values, mut coefficients) = {
                let _phase = dft_numerics::pool::analysis::Span::serial("evd");
                hermitian_eigh(&projected)?
            };
            reduced_cache = Some(projected);
            values.truncate(self.bands);
            coefficients.resize_columns(self.bands)?;
            gemm_work.accounting_site("ritz_and_residual");
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
                work.apply(self.bands, Some(self.bands - policy.strict_bands));
                operator_applications += self.bands;
                {
                    if TRACE {
                        trace.blocks += 1;
                    }
                    residual(operator, &vectors, &values)?
                }
            };
            let mut products_trusted = true;
            let mut residuals_certified = !recurrence;
            if !ritz
                && recurrence
                && (certify || iteration == self.max_iterations)
                && (meets_policy(&residuals)
                    || (iteration == self.max_iterations && meets_strict(&residuals)))
            {
                work.apply(self.bands, Some(self.bands - policy.strict_bands));
                operator_applications += self.bands;
                (residual_block, residuals) = {
                    if TRACE {
                        trace.blocks += 1;
                    }
                    residual(operator, &vectors, &values)?
                };
                residuals_certified = true;
                products_trusted = meets_policy(&residuals);
                // A failed certificate supplies actual correction directions and
                // its actual maximum is recorded in this iteration's history.
            }
            if TRACE {
                trace.step(basis.cols(), operator_applications, &values, &residuals);
            }
            let maximum = max_residual(&residuals);
            work.residuals(iteration, &residuals);
            history.push(maximum);
            if ritz {
                for j in 0..self.bands {
                    changes[j] = if iteration == 1 {
                        f64::INFINITY
                    } else {
                        (values[j] - previous_values[j]).abs()
                    };
                }
                previous_values.copy_from_slice(&values);
            }
            let stopping = if ritz { &changes } else { &residuals };
            let root_converged = |j: usize| {
                if ritz {
                    stopping[j] < threshold(j)
                } else {
                    stopping[j] <= threshold(j)
                }
            };
            let finished = if ritz {
                (0..self.bands).all(root_converged)
                    || (iteration == self.max_iterations
                        && (0..policy.strict_bands).all(root_converged))
            } else {
                meets_policy(&residuals)
                    || (iteration == self.max_iterations && meets_strict(&residuals))
            };
            if finished {
                work.returned();
                return Ok(EigenSolution {
                    locking: None,
                    spectral_completeness: None,
                    eigenvalues: values,
                    eigenvectors: vectors,
                    converged: (0..self.bands).map(root_converged).collect(),
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
            // A full-space first solve has no expansion direction. Repeat its
            // small solve once to establish an actual consecutive Ritz change.
            if ritz && iteration == 1 && basis.cols() == n {
                continue;
            }
            let selected: Vec<_> = (0..self.bands)
                .filter(|&j| !root_converged(j))
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
                    p.apply_with_vectors(&shifts, &corrections, &vectors, &selected, &mut output)?;
                }
                if output.rows() != n || output.cols() != corrections.cols() {
                    return Err(LinalgError::Dimension.into());
                }
                output.validate_finite()?;
                corrections = output;
            }
            cache_valid = products_trusted && basis.cols() + corrections.cols() <= capacity;
            if !cache_valid {
                reduced_cache = None;
                work.restart();
                if TRACE {
                    trace.restarts += 1;
                }
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
            let mut correction_origins_known = true;
            // A stalled direction can occur with an exact preconditioner or an invariant
            // initial space. Try fresh independent directions, never infer convergence.
            if extra.cols() == 0 {
                correction_origins_known = false;
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
                if work_accounting::ENABLED {
                    work.apply(
                        extra.cols(),
                        (correction_origins_known && extra.cols() == selected.len()).then(|| {
                            selected
                                .iter()
                                .filter(|&&j| j >= policy.strict_bands)
                                .count()
                        }),
                    );
                }
                appended_products.resize_columns(extra.cols())?;
                let products = &mut appended_products;
                if TRACE {
                    trace.blocks += 1;
                }
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
        gemm_work.accounting_site("restart_contraction");
        gemm_work.multiply_into(h_basis, Identity, coefficients, Identity, cached)?;
        Ok((vectors, true))
    } else {
        Ok((orthonormalize(&vectors, None)?.vectors, false))
    }
}

#[cfg(test)]
use dft_numerics::linalg::gemm;
