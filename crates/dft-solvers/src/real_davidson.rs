//! Isolated ORC Davidson; no complex-path dispatch or arithmetic is changed.
use crate::bordered::MAX_STREAK;
use crate::davidson_stages::{B1Count, B5Count, Stage, Timers, b1_add, b5_add};
pub use crate::davidson_stages::{
    B1Counters, B1Mode, B5Counters, ResearchOptions, SketchKind, StageSeconds, WIDTH_LABELS,
    b1_counters, b5_counters, stage_seconds,
};
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
    /// Column of each root in the Ritz and residual blocks of the iteration.
    column_of: Vec<usize>,
    /// B1 scratch, allocated at the first B1-active solve only.
    b1: Option<B1Scratch>,
    /// B5 base decomposition and scratch, allocated at the first B5 solve.
    b5: Option<crate::bordered::Bordered>,
    /// B6 projector-product cache, allocated at the first B6-active solve.
    b6: Option<B6Cache>,
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
            column_of: (0..bands).collect(),
            b1: None,
            b5: None,
            b6: None,
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
        self.solve_research(
            op,
            initial,
            diagonal,
            policy,
            ritz,
            certify,
            capacity_per_band,
            workspace,
            preconditioner,
            ResearchOptions::from_env(),
        )
    }

    /// Research entry point with explicit switches (the other entry points read
    /// them from the environment, see `davidson_stages`). With
    /// `ResearchOptions::default()` this is the default path, bit for bit.
    ///
    /// B1 (selected-root products): in a narrow iteration only the roots that
    /// can receive a correction get Ritz vectors, H images and residual vectors;
    /// restart, return, the iteration limit and wide iterations (at least half
    /// the bands) take the all-band path. Standard (Ritz-change) solves use the
    /// exact Ritz-value flags. Strict (residual) solves flag a root converged for
    /// one iteration when its sketch estimate ||S H V c - theta S V c|| is below
    /// half its threshold; every other root gets an explicit residual. A return
    /// is always decided on the all-band path with explicit residuals and the
    /// unchanged certification. In narrow iterations `history` holds, for
    /// standard solves, the largest residual of the corrected roots and, for
    /// strict solves, the largest of the explicit residuals and the estimates.
    /// B1 is inactive while `BAANDI_SECOND_VECTOR` is set.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_research(
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
        options: ResearchOptions,
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
            column_of,
            b1: b1_scratch,
            b5: b5_scratch,
            b6: b6_scratch,
            ..
        } = workspace.as_mut().unwrap();
        // B6 (`BAANDI_B6=on`, default off): cache P^T V per basis column so the
        // return certificate skips its first projector product.
        let mut b6 = B6Cache::activate(b6_scratch, op, capacity, self.bands, options.b6);
        let timers = Timers::new(options.stage_timers);
        let second_vector = SecondVector::from_env();
        let b1_active = options.b1.acts_on(ritz) && second_vector == SecondVector::Off;
        let sketching = b1_active && !ritz;
        if b1_active {
            b1_add(B1Count::Solves, 1);
            if b1_scratch
                .as_ref()
                .is_none_or(|s| !s.matches(self.seed, sketching, options.sketch))
            {
                *b1_scratch = Some(B1Scratch::new(
                    n,
                    self.bands,
                    capacity,
                    self.seed,
                    sketching.then_some(options.sketch),
                )?);
            }
        }
        let mut b1 = if b1_active { b1_scratch.as_mut() } else { None };
        if options.b5 && b5_scratch.is_none() {
            *b5_scratch = Some(crate::bordered::Bordered::new(capacity));
        }
        let mut b5 = if options.b5 {
            b5_scratch.as_mut()
        } else {
            None
        };
        if let Some(b) = b5.as_mut() {
            b.invalidate();
        }
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
        apply_into(op, basis, &mut applications, &mut work, products, timers)?;
        if let Some(c) = b6.as_mut() {
            c.load(op, 0, basis.cols());
        }
        if let Some(sk) = b1.as_mut().and_then(|b| b.sketch.as_mut()) {
            let t = timers.start();
            sk.refresh(gemm, basis, products)?;
            timers.stop(t, Stage::Sketch);
        }
        let mut previous_corrections = if second_vector == SecondVector::PreviousCorrection {
            vec![0.; n * self.bands]
        } else {
            Vec::new()
        };
        let mut has_previous = vec![false; self.bands];
        let mut max_basis_dimension = 0;
        for iteration in 1..=self.max_iterations {
            let iteration_clock = timers.start();
            max_basis_dimension = max_basis_dimension.max(basis.cols());
            let m = basis.cols();
            let t = timers.start();
            reduced_matrix_into(basis, products, reduced_width, gemm, reduced)?;
            timers.stop(t, Stage::Reduced);
            let prefix = reduced_width;
            reduced_width = m;
            coefficients.reshape_within_storage(m, m)?;
            let t = timers.start();
            // B5: one appended column on an unbroken base -> bordered update.
            let bordered = match b5.as_mut() {
                Some(b) if prefix + 1 == m && b.order == prefix && b.streak < MAX_STREAK => {
                    let done = b.update(gemm, reduced, &mut values[..m], coefficients)?;
                    if !done {
                        b5_add(B5Count::Declined);
                    }
                    done
                }
                _ => false,
            };
            if !bordered {
                eigensolver.solve_into(reduced, &mut values[..m], coefficients)?;
            }
            if let Some(b) = b5.as_mut() {
                if bordered {
                    b5_add(B5Count::Bordered);
                    b.streak += 1;
                } else {
                    b5_add(B5Count::Full);
                    b.streak = 0;
                }
                b.store(&values[..m], coefficients);
            }
            timers.stop(t, Stage::Eigen);
            // B1: plan a narrow iteration (never at the iteration limit).
            let mut narrow = false;
            if let Some(b1) = b1.as_mut() {
                if let Some(sk) = b1.sketch.as_mut() {
                    let t = timers.start();
                    sk.estimate(gemm, coefficients, m, &values[..self.bands], &mut b1.est)?;
                    timers.stop(t, Stage::Sketch);
                    b1_add(B1Count::Estimates, self.bands);
                }
                b1.set.clear();
                if iteration < self.max_iterations {
                    if ritz {
                        if !(iteration == 1 && basis.cols() == n) {
                            b1.set.extend(
                                (0..self.bands)
                                    .filter(|&j| {
                                        !(iteration > 1
                                            && (values[j] - previous[j]).abs() < threshold(j))
                                    })
                                    .take(self.bands.min(capacity - self.bands)),
                            );
                        }
                    } else {
                        let est = &b1.est;
                        b1.set
                            .extend((0..self.bands).filter(|&j| est[j] >= 0.5 * threshold(j)));
                    }
                    narrow = !b1.set.is_empty()
                        && 2 * b1.set.len() < self.bands
                        && basis.cols() + b1.set.len() <= capacity;
                }
                if narrow {
                    let k = b1.set.len();
                    for (slot, &j) in b1.set.iter().enumerate() {
                        b1.coeff[slot * m..(slot + 1) * m]
                            .copy_from_slice(&coefficients.as_slice()[j * m..(j + 1) * m]);
                        b1.set_values[slot] = values[j];
                        column_of[j] = slot;
                    }
                    vectors.reshape_within_storage(n, k)?;
                    rotated_products.reshape_within_storage(n, k)?;
                    let t = timers.start();
                    gemm.rotate_pair_slices(
                        basis.as_slice(),
                        products.as_slice(),
                        m,
                        &b1.coeff[..m * k],
                        k,
                        vectors.as_mut_slice(),
                        rotated_products.as_mut_slice(),
                    )?;
                    timers.stop(t, Stage::Rotate);
                    residual_block.reshape_within_storage(n, k)?;
                    let t = timers.start();
                    residual_into(
                        vectors,
                        rotated_products,
                        &b1.set_values[..k],
                        residual_block,
                        &mut b1.set_norms[..k],
                        &pool,
                    )?;
                    timers.stop(t, Stage::Residual);
                    if !ritz {
                        residuals.copy_from_slice(&b1.est);
                    }
                    for (slot, &j) in b1.set.iter().enumerate() {
                        residuals[j] = b1.set_norms[slot];
                    }
                    // Strict: if every explicit residual passes, all roots are
                    // flagged converged; a return is decided on the full path.
                    if !ritz && b1.set.iter().all(|&j| residuals[j] <= threshold(j)) {
                        narrow = false;
                        b1_add(B1Count::Fallbacks, 1);
                    }
                }
                if narrow {
                    b1_add(B1Count::Narrow, 1);
                    b1_add(B1Count::NarrowColumns, b1.set.len());
                    b1_add(B1Count::Avoided, self.bands - b1.set.len());
                } else {
                    b1_add(B1Count::Full, 1);
                }
            }
            let mut certified = false;
            let mut trusted = true;
            if !narrow {
                for (j, c) in column_of.iter_mut().enumerate() {
                    *c = j;
                }
                vectors.reshape_within_storage(n, self.bands)?;
                rotated_products.reshape_within_storage(n, self.bands)?;
                let t = timers.start();
                gemm.rotate_pair_slices(
                    basis.as_slice(),
                    products.as_slice(),
                    m,
                    &coefficients.as_slice()[..m * self.bands],
                    self.bands,
                    vectors.as_mut_slice(),
                    rotated_products.as_mut_slice(),
                )?;
                timers.stop(t, Stage::Rotate);
                residual_block.reshape_within_storage(n, self.bands)?;
                let t = timers.start();
                residual_into(
                    vectors,
                    rotated_products,
                    &values[..self.bands],
                    residual_block,
                    residuals,
                    &pool,
                )?;
                timers.stop(t, Stage::Residual);
                if sketching && let Some(b1) = b1.as_ref() {
                    let missed = (0..self.bands)
                        .filter(|&j| b1.est[j] < 0.5 * threshold(j) && residuals[j] > threshold(j))
                        .count();
                    b1_add(B1Count::Disagreements, missed);
                }
                if !ritz
                    && (certify || iteration == self.max_iterations)
                    && (meets_policy(residuals)
                        || (iteration == self.max_iterations && meets_strict(residuals)))
                {
                    certificates += 1;
                    match b6.as_mut().filter(|c| c.covers(m)) {
                        Some(c) => {
                            let pc = c.rotate(coefficients.as_slice(), m, self.bands);
                            apply_with_projections_into(
                                op,
                                vectors,
                                pc,
                                &mut applications,
                                &mut work,
                                actual,
                                timers,
                            )?
                        }
                        None => {
                            apply_into(op, vectors, &mut applications, &mut work, actual, timers)?
                        }
                    }
                    let t = timers.start();
                    residual_into(
                        vectors,
                        actual,
                        &values[..self.bands],
                        residual_block,
                        residuals,
                        &pool,
                    )?;
                    timers.stop(t, Stage::Residual);
                    certified = true;
                    trusted = meets_policy(residuals);
                }
            }
            let maximum = match b1.as_ref() {
                Some(b1) if narrow && ritz => {
                    b1.set.iter().map(|&j| residuals[j]).fold(0., f64::max)
                }
                _ => residuals.iter().copied().fold(0., f64::max),
            };
            history.push(maximum);
            work.residuals(iteration, residuals);
            for (j, flag) in converged.iter_mut().enumerate() {
                *flag = if ritz {
                    iteration > 1 && (values[j] - previous[j]).abs() < threshold(j)
                } else if narrow {
                    // Not in the explicit set: the sketch estimate is below
                    // half the threshold.
                    b1.as_ref().is_some_and(|b| b.est[j] < 0.5 * threshold(j))
                        || residuals[j] <= threshold(j)
                } else {
                    residuals[j] <= threshold(j)
                };
            }
            previous.copy_from_slice(&values[..self.bands]);
            let finished = converged.iter().all(|&b| b)
                || (iteration == self.max_iterations
                    && converged[..policy.strict_bands].iter().all(|&b| b));
            if narrow && finished {
                return Err(SolverError::Application(
                    "B1 internal error: narrow iteration reached a return".into(),
                ));
            }
            if finished {
                timers.iteration(iteration_clock, 0);
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
                timers.iteration(iteration_clock, 0);
                result_residuals.copy_from_slice(residuals);
                return Err(SolverError::NotConverged {
                    iterations: iteration,
                    max_residual: maximum,
                    residuals: result_residuals,
                });
            }
            if ritz && iteration == 1 && basis.cols() == n {
                timers.iteration(iteration_clock, 0);
                continue;
            }
            selected.clear();
            selected.extend(
                (0..self.bands)
                    .filter(|&j| !converged[j])
                    .take(self.bands.min(capacity - self.bands)),
            );
            corrections.reshape_within_storage(n, selected.len())?;
            let t = timers.start();
            let column_slots: &[usize] = &column_of[..];
            pool.install(|| {
                column_work(corrections.as_mut_slice(), n, 0, &|k, dst| {
                    let j = selected[k];
                    // Column of root j in the Ritz and residual blocks (j on the
                    // all-band path, its slot in a B1 narrow iteration).
                    let c = column_slots[j];
                    // The default damping keeps the original tight loop; research
                    // variants take their own loops (no per-element dispatch).
                    match preconditioner {
                        crate::DavidsonPreconditioner::Default => {
                            for ((dst, &r), &d) in
                                dst.iter_mut().zip(column(residual_block, c)).zip(diagonal)
                            {
                                *dst = precondition(r, d, values[j])?;
                            }
                        }
                        crate::DavidsonPreconditioner::Tpa => {
                            let scale = crate::kinetic_scale(
                                column(vectors, c)
                                    .iter()
                                    .zip(diagonal)
                                    .map(|(&c, &t)| (c * c, t)),
                            )?;
                            for ((dst, &r), &d) in
                                dst.iter_mut().zip(column(residual_block, c)).zip(diagonal)
                            {
                                *dst = -r * crate::tpa_weight(d, scale);
                                finite(&[*dst])?;
                            }
                        }
                        crate::DavidsonPreconditioner::RitzShift => {
                            for ((dst, &r), &d) in
                                dst.iter_mut().zip(column(residual_block, c)).zip(diagonal)
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
            timers.stop(t, Stage::Precondition);
            if second_vector != SecondVector::Off {
                add_second_vector(
                    second_vector,
                    n,
                    SecondVectorContext {
                        selected,
                        converged,
                        residuals,
                        strict_bands: policy.strict_bands,
                        bands: self.bands,
                        values: &values[..self.bands],
                        diagonal,
                        residual_block,
                    },
                    corrections,
                    previous_corrections.as_mut_slice(),
                    has_previous.as_mut_slice(),
                )?;
            }
            if !trusted || basis.cols() + corrections.cols() > capacity {
                if narrow {
                    return Err(SolverError::Application(
                        "B1 internal error: narrow iteration reached a restart".into(),
                    ));
                }
                restarts += 1;
                work.restart();
                reduced_width = 0;
                if trusted && orthogonality_defect(vectors) <= 1e-12 {
                    contractions += 1;
                    std::mem::swap(basis, vectors);
                    std::mem::swap(products, rotated_products);
                    if let Some(c) = b6.as_mut() {
                        // Same rotation as the basis: P^T V <- (P^T V) C.
                        c.contract(coefficients.as_slice(), m, self.bands);
                    }
                    if let Some(sk) = b1.as_mut().and_then(|b| b.sketch.as_mut()) {
                        // Same rotation as the basis: A <- A C, B <- B C.
                        sk.contract(self.bands);
                    }
                } else {
                    let t = timers.start();
                    projection.run_into(vectors, None, gemm, basis)?;
                    timers.stop(t, Stage::Orthonormalize);
                    apply_into(op, basis, &mut applications, &mut work, products, timers)?;
                    if let Some(c) = b6.as_mut() {
                        c.load(op, 0, basis.cols());
                    }
                    if let Some(sk) = b1.as_mut().and_then(|b| b.sketch.as_mut()) {
                        let t = timers.start();
                        sk.refresh(gemm, basis, products)?;
                        timers.stop(t, Stage::Sketch);
                    }
                }
            }
            let t = timers.start();
            projection.run_into(corrections, Some(basis), gemm, extra)?;
            if extra.cols() == 0 {
                random.reshape_within_storage(n, self.bands.min(capacity - basis.cols()))?;
                for value in random.as_mut_slice() {
                    *value = rng.next();
                }
                projection.run_into(random, Some(basis), gemm, extra)?;
            }
            timers.stop(t, Stage::Orthonormalize);
            if extra.cols() == 0 {
                timers.iteration(iteration_clock, 0);
                result_residuals.copy_from_slice(residuals);
                return Err(SolverError::NotConverged {
                    iterations: iteration,
                    max_residual: maximum,
                    residuals: result_residuals,
                });
            }
            apply_into(
                op,
                extra,
                &mut applications,
                &mut work,
                extra_products,
                timers,
            )?;
            if let Some(c) = b6.as_mut() {
                c.load(op, basis.cols(), extra.cols());
            }
            if let Some(sk) = b1.as_mut().and_then(|b| b.sketch.as_mut()) {
                let t = timers.start();
                sk.append(gemm, basis.cols(), extra, extra_products)?;
                timers.stop(t, Stage::Sketch);
            }
            append_into(products, extra_products)?;
            append_into(basis, extra)?;
            timers.iteration(iteration_clock, corrections.cols());
        }
        unreachable!("validated positive iteration limit")
    }
}

/// Lowest kinetic coordinate units; ties retain ORC half-list order (0, cos, sin).
/// The diagonal is kinetic, not the full Hamiltonian diagonal.
/// Experimental (docs/roadmap/second-vector-plan.md): fill the free slot of an
/// odd-sized correction block. Off unless `BAANDI_SECOND_VECTOR=v1|v2|v3`.
#[derive(Clone, Copy, PartialEq, Eq)]
enum SecondVector {
    Off,
    /// Previous correction of the selected band with the largest residual.
    PreviousCorrection,
    /// Preconditioned residual of the converged band with the largest residual.
    ConvergedBand,
    /// Unpreconditioned residual of the worst strict band.
    WorstResidual,
    /// As `WorstResidual`, but only when the correction block has exactly one
    /// column (the single-column tail); never in wider blocks.
    WorstResidualTail,
}
static SECOND_VECTOR_COLUMNS: std::sync::atomic::AtomicUsize =
    std::sync::atomic::AtomicUsize::new(0);
/// Number of second-vector columns inserted so far in this process (for the run
/// report: lets a measurement check that the gate itself fired).
pub fn second_vector_columns() -> usize {
    SECOND_VECTOR_COLUMNS.load(std::sync::atomic::Ordering::Relaxed)
}
impl SecondVector {
    fn from_env() -> Self {
        match std::env::var("BAANDI_SECOND_VECTOR").as_deref() {
            Ok("v1") => Self::PreviousCorrection,
            Ok("v2") => Self::ConvergedBand,
            Ok("v3") => Self::WorstResidual,
            Ok("v4") => Self::WorstResidualTail,
            _ => Self::Off,
        }
    }
}
struct SecondVectorContext<'a> {
    selected: &'a [usize],
    converged: &'a [bool],
    residuals: &'a [f64],
    strict_bands: usize,
    bands: usize,
    values: &'a [f64],
    diagonal: &'a [f64],
    residual_block: &'a RealMatrix,
}
fn add_second_vector(
    variant: SecondVector,
    n: usize,
    c: SecondVectorContext<'_>,
    corrections: &mut RealMatrix,
    previous: &mut [f64],
    has_previous: &mut [bool],
) -> Result<(), SolverError> {
    let k = c.selected.len();
    // V1 remembers this iteration's corrections for the next one.
    if variant == SecondVector::PreviousCorrection {
        let stored: Vec<(usize, usize)> = c.selected.iter().copied().enumerate().collect();
        let mut extra: Option<usize> = None;
        if k % 2 == 1 && k < c.bands {
            extra = c
                .selected
                .iter()
                .copied()
                .filter(|&j| has_previous[j])
                .max_by(|&a, &b| c.residuals[a].total_cmp(&c.residuals[b]));
        }
        let mut column_extra = None;
        if let Some(j) = extra {
            column_extra = Some(previous[j * n..(j + 1) * n].to_vec());
        }
        for (slot, j) in stored {
            previous[j * n..(j + 1) * n]
                .copy_from_slice(&corrections.as_slice()[slot * n..(slot + 1) * n]);
            has_previous[j] = true;
        }
        if let Some(col) = column_extra {
            SECOND_VECTOR_COLUMNS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            corrections.reshape_within_storage(n, k + 1)?;
            corrections.as_mut_slice()[k * n..(k + 1) * n].copy_from_slice(&col);
        }
        return Ok(());
    }
    if k.is_multiple_of(2) || k >= c.bands {
        return Ok(());
    }
    if variant == SecondVector::WorstResidualTail && k != 1 {
        return Ok(());
    }
    let chosen = match variant {
        SecondVector::ConvergedBand => (0..c.bands)
            .filter(|&j| c.converged[j])
            .max_by(|&a, &b| c.residuals[a].total_cmp(&c.residuals[b])),
        SecondVector::WorstResidual | SecondVector::WorstResidualTail => {
            (0..c.strict_bands.min(c.bands))
                .max_by(|&a, &b| c.residuals[a].total_cmp(&c.residuals[b]))
        }
        _ => None,
    };
    let Some(j) = chosen else {
        return Ok(());
    };
    corrections.reshape_within_storage(n, k + 1)?;
    SECOND_VECTOR_COLUMNS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let dst = &mut corrections.as_mut_slice()[k * n..(k + 1) * n];
    let r = column(c.residual_block, j);
    match variant {
        SecondVector::ConvergedBand => {
            for ((dst, &r), &d) in dst.iter_mut().zip(r).zip(c.diagonal) {
                *dst = precondition(r, d, c.values[j])?;
            }
        }
        _ => dst.copy_from_slice(r),
    }
    Ok(())
}

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
fn mgs_dot(wide: bool, x: &[f64], y: &[f64]) -> f64 {
    if wide && x.len() == y.len() {
        dot_wide(x, y)
    } else {
        dot(x, y)
    }
}
/// WIDE-DOT: eight independent lanes (element i adds into lane i % 8, in
/// ascending order) combined pairwise. The four-lane `dot` is bound by the
/// latency of its four ordered add chains; eight chains let the core overlap
/// them (Mac observation 5 -> 10 elements/ns at n = 11877). Differs from `dot`
/// by rounding only; independent of the thread count (serial per column).
fn dot_wide(x: &[f64], y: &[f64]) -> f64 {
    const L: usize = 8;
    let mut s = [0.; L];
    let ((cx, tx), (cy, ty)) = (x.as_chunks::<L>(), y.as_chunks::<L>());
    for (a, b) in cx.iter().zip(cy) {
        for k in 0..L {
            s[k] += a[k] * b[k];
        }
    }
    for ((s, &a), &b) in s.iter_mut().zip(tx).zip(ty) {
        *s += a * b;
    }
    ((s[0] + s[4]) + (s[2] + s[6])) + ((s[1] + s[5]) + (s[3] + s[7]))
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
    timers: Timers,
) -> Result<(), SolverError> {
    let clock = timers.start();
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
    timers.stop(clock, Stage::HApply);
    Ok(())
}
/// `apply_into` for the B6 certificate: the operator receives P^T X.
fn apply_with_projections_into(
    op: &dyn RealLinearOperator,
    x: &RealMatrix,
    projections: &[f64],
    count: &mut usize,
    work: &mut crate::work_accounting::Trace,
    y: &mut RealMatrix,
    timers: Timers,
) -> Result<(), SolverError> {
    let clock = timers.start();
    work.apply(x.cols(), None);
    y.reshape_within_storage(x.rows(), x.cols())?;
    op.apply_block_with_projections(x, projections, y)?;
    if y.rows() != x.rows() || y.cols() != x.cols() {
        return Err(LinalgError::Dimension.into());
    }
    finite(y.as_slice())?;
    *count += x.cols();
    timers.stop(clock, Stage::HApply);
    Ok(())
}
/// B6 cache: P^T V for the live basis columns (column-major, `count` rows).
/// Every basis change updates it (append: copied from the operator; contraction:
/// the same rotation; rebuild: reloaded); any missing copy invalidates it, and
/// an invalid cache falls back to the plain apply. Products use fixed-order
/// loops, so results do not depend on the thread count.
struct B6Cache {
    count: usize,
    pv: Vec<f64>,
    pc: Vec<f64>,
    cols: usize,
    valid: bool,
}
impl B6Cache {
    fn activate<'a>(
        slot: &'a mut Option<B6Cache>,
        op: &dyn RealLinearOperator,
        capacity: usize,
        bands: usize,
        enabled: bool,
    ) -> Option<&'a mut B6Cache> {
        let count = op.projector_count();
        if !enabled || count == 0 {
            return None;
        }
        if slot
            .as_ref()
            .is_none_or(|c| c.count != count || c.pv.len() != count * capacity)
        {
            *slot = Some(B6Cache {
                count,
                pv: vec![0.; count * capacity],
                pc: vec![0.; count * bands],
                cols: 0,
                valid: false,
            });
        }
        let cache = slot.as_mut().unwrap();
        cache.cols = 0;
        cache.valid = false;
        Some(cache)
    }
    /// Columns `first..first + k` of the basis were just applied by `op`.
    fn load(&mut self, op: &dyn RealLinearOperator, first: usize, k: usize) {
        if first == 0 {
            self.valid = true;
        }
        let range = first * self.count..(first + k) * self.count;
        self.valid = self.valid
            && first == self.cols
            && range.end <= self.pv.len()
            && op.copy_last_projections(k, &mut self.pv[range]);
        self.cols = first + k;
    }
    fn covers(&self, m: usize) -> bool {
        self.valid && self.cols == m
    }
    /// pc = pv[:, :m] * C[:, :bands] (C column-major, leading dimension m).
    fn rotate(&mut self, c: &[f64], m: usize, bands: usize) -> &[f64] {
        let count = self.count;
        for b in 0..bands {
            let out = &mut self.pc[b * count..(b + 1) * count];
            out.fill(0.);
            for k in 0..m {
                let factor = c[b * m + k];
                for (o, &p) in out.iter_mut().zip(&self.pv[k * count..(k + 1) * count]) {
                    *o += p * factor;
                }
            }
        }
        &self.pc[..bands * count]
    }
    fn contract(&mut self, c: &[f64], m: usize, bands: usize) {
        if !self.covers(m) {
            self.valid = false;
            return;
        }
        self.rotate(c, m, bands);
        let len = bands * self.count;
        self.pv[..len].copy_from_slice(&self.pc[..len]);
        self.cols = bands;
    }
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
    let rows = v.rows();
    ws.multiply_slices_into(
        dft_numerics::GemmSite::Reduced,
        v.as_slice(),
        rows,
        v.cols(),
        rows,
        RealOp::Transpose,
        &hv.as_slice()[prefix * rows..],
        m - prefix,
        rows,
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
/// Panel width of the blocked Gram-Schmidt switch.
const BLOCK_COLUMNS: usize = 16;
struct Projection {
    pool: std::sync::Arc<WorkerPool>,
    normalized: RealMatrix,
    overlap: RealMatrix,
    projected: RealMatrix,
    output: RealMatrix,
    column: Vec<f64>,
    /// Panel-wise block Gram-Schmidt (default on MKL builds, `BAANDI_MGS` overrides).
    blocked: bool,
    /// WIDE-DOT (opt-in `BAANDI_WIDE_DOT=on`, default off): the in-block MGS
    /// dot products use eight lane accumulators (`dot_wide`), rounding class.
    wide: bool,
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
            // Blocked Gram-Schmidt is the default on MKL builds only (dev08
            // evidence MGS-B R2); `BAANDI_MGS=classic` or `blocked` overrides.
            blocked: match std::env::var("BAANDI_MGS").as_deref() {
                Ok("blocked") => true,
                Ok("classic") => false,
                _ => dft_numerics::MKL_BACKEND,
            },
            wide: std::env::var("BAANDI_WIDE_DOT").as_deref() == Ok("on"),
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
                ws.multiply_slices_into(
                    dft_numerics::GemmSite::ProjectionOverlap,
                    q.as_slice(),
                    q.rows(),
                    q.cols(),
                    q.rows(),
                    RealOp::Transpose,
                    self.normalized.as_slice(),
                    count,
                    n,
                    self.overlap.as_mut().submatrix_mut(0, 0, q.cols(), count),
                )?;
                let overlap_rows = self.overlap.rows();
                ws.multiply_slices_into(
                    dft_numerics::GemmSite::ProjectionApply,
                    q.as_slice(),
                    q.rows(),
                    q.cols(),
                    q.rows(),
                    RealOp::Identity,
                    self.overlap.as_slice(),
                    count,
                    overlap_rows,
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
            if self.blocked {
                // Panels of BLOCK_COLUMNS columns: two block projections against
                // the kept columns (GEMM), then the same two-pass MGS and rank
                // floor inside the panel.
                let mut start = 0;
                while start < count {
                    let w = BLOCK_COLUMNS.min(count - start);
                    for _ in 0..2 {
                        if kept == 0 {
                            break;
                        }
                        ws.multiply_slices_into(
                            dft_numerics::GemmSite::Other,
                            self.output.as_slice(),
                            n,
                            kept,
                            n,
                            RealOp::Transpose,
                            &self.normalized.as_slice()[start * n..(start + w) * n],
                            w,
                            n,
                            self.overlap.as_mut().submatrix_mut(0, 0, kept, w),
                        )?;
                        let overlap_rows = self.overlap.rows();
                        ws.multiply_slices_into(
                            dft_numerics::GemmSite::Other,
                            self.output.as_slice(),
                            n,
                            kept,
                            n,
                            RealOp::Identity,
                            self.overlap.as_slice(),
                            w,
                            overlap_rows,
                            self.projected.as_mut().subcols_mut(0, w),
                        )?;
                        let projected = &self.projected;
                        self.pool.install(|| {
                            column_work(
                                &mut self.normalized.as_mut_slice()[start * n..(start + w) * n],
                                n,
                                0,
                                &|j, dst| {
                                    for (dst, &src) in dst.iter_mut().zip(column(projected, j)) {
                                        *dst -= src;
                                    }
                                    Ok(())
                                },
                            )
                        })?;
                    }
                    let panel_start = kept;
                    for jj in start..start + w {
                        self.column.copy_from_slice(column(&self.normalized, jj));
                        for _ in 0..2 {
                            for k in panel_start..kept {
                                let p = mgs_dot(self.wide, &self.column, column(&self.output, k));
                                for (dst, &src) in
                                    self.column.iter_mut().zip(column(&self.output, k))
                                {
                                    *dst -= p * src;
                                }
                            }
                        }
                        let len = norm(&self.column);
                        finite(&[len])?;
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
                    start += w;
                }
            } else {
                for j in 0..count {
                    self.column.copy_from_slice(column(&self.normalized, j));
                    for _ in 0..2 {
                        for k in 0..kept {
                            let p = mgs_dot(self.wide, &self.column, column(&self.output, k));
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
/// Rows of the B1 strict sketch (capped by the subspace capacity, because the
/// shared GEMM workspace bounds the short dimension). Relative estimate error is
/// about 1/sqrt(2 s): 6 percent at s = 128.
const SKETCH_ROWS: usize = 128;
/// Stream offset of the sketch generator, so the solver's own random fallback
/// sequence is unchanged by B1.
const SKETCH_STREAM: u64 = 0x5b1c_e7c4_2f6a_9d01;

/// Scratch of B1 (selected-root products); allocated once per workspace.
struct B1Scratch {
    seed: u64,
    /// Strict-profile sketch; absent for Ritz-change solves.
    sketch: Option<Sketch>,
    /// Sketch estimates per root (strict), in Ha.
    est: Vec<f64>,
    /// Gathered Ritz coefficients of the narrow set, m x k.
    coeff: Vec<f64>,
    /// Roots formed explicitly in a narrow iteration, ascending.
    set: Vec<usize>,
    set_values: Vec<f64>,
    set_norms: Vec<f64>,
}
impl B1Scratch {
    fn new(
        n: usize,
        bands: usize,
        capacity: usize,
        seed: u64,
        sketch: Option<SketchKind>,
    ) -> Result<Self, SolverError> {
        Ok(Self {
            seed,
            sketch: match sketch {
                Some(kind) => Some(Sketch::new(n, bands, capacity, seed, kind)?),
                None => None,
            },
            est: vec![0.; bands],
            coeff: vec![0.; capacity * bands],
            set: Vec::with_capacity(bands),
            set_values: vec![0.; bands],
            set_norms: vec![0.; bands],
        })
    }
    fn matches(&self, seed: u64, sketching: bool, kind: SketchKind) -> bool {
        self.seed == seed
            && match &self.sketch {
                Some(s) => sketching && s.kind() == kind,
                None => !sketching,
            }
    }
}

/// Fixed-seed sign sketch S (s x N) with the running products A = S V and
/// B = S (H V), one column per basis column. Dense form (proposal): entries
/// +-1/sqrt(s), stored transposed as N x s, s N multiply-adds per column. Count
/// form (CountSketch): one +-1 entry per coordinate in a seeded bucket, N adds
/// per column; same expectation and the same variance bound 2 ||x||^4 / s of
/// ||S x||^2. All products use fixed-order loops or the fixed-order kernels of
/// the shared GEMM workspace, so estimates do not depend on the thread budget.
struct Sketch {
    rows: usize,
    /// Requested form (the identity is used for small spaces of either kind).
    kind: SketchKind,
    form: SketchForm,
    /// S V and S H V, s x capacity (live width = basis width).
    a: RealMatrix,
    b: RealMatrix,
    /// S V C and S H V C for the first `bands` Ritz coefficients, s x bands.
    ac: RealMatrix,
    bc: RealMatrix,
}
impl Sketch {
    fn new(
        n: usize,
        bands: usize,
        capacity: usize,
        seed: u64,
        kind: SketchKind,
    ) -> Result<Self, SolverError> {
        // The dense form is bounded by the GEMM workspace's short dimension;
        // the count form is not. When the space is no larger than the sketch,
        // S is the identity (exact residual norms; no sketch error at all).
        let rows = match kind {
            SketchKind::Dense => SKETCH_ROWS.min(capacity),
            SketchKind::Count => SKETCH_ROWS.min(n),
        };
        let mut rng = Random(seed ^ SKETCH_STREAM);
        let form = match kind {
            _ if n <= rows => SketchForm::Count {
                bucket: (0..n).collect(),
                sign: vec![1.; n],
            },
            SketchKind::Dense => {
                let scale = 1. / (rows as f64).sqrt();
                let mut transposed = RealMatrix::zeros(n, rows)?;
                for v in transposed.as_mut_slice() {
                    *v = if rng.next() < 0. { -scale } else { scale };
                }
                SketchForm::Dense(transposed)
            }
            SketchKind::Count => {
                let mut bucket = Vec::with_capacity(n);
                let mut sign = Vec::with_capacity(n);
                for _ in 0..n {
                    // next() is uniform in [-0.5, 0.5).
                    let b = ((rng.next() + 0.5) * rows as f64) as usize;
                    bucket.push(b.min(rows - 1));
                    sign.push(if rng.next() < 0. { -1. } else { 1. });
                }
                SketchForm::Count { bucket, sign }
            }
        };
        Ok(Self {
            rows,
            kind,
            form,
            a: RealMatrix::zeros(rows, capacity)?,
            b: RealMatrix::zeros(rows, capacity)?,
            ac: RealMatrix::zeros(rows, bands)?,
            bc: RealMatrix::zeros(rows, bands)?,
        })
    }
    fn kind(&self) -> SketchKind {
        self.kind
    }
    /// out[:, first..first + k] = S X for the k columns of `x`.
    fn project(
        form: &SketchForm,
        gemm: &mut RealGemmWorkspace,
        x: &RealMatrix,
        out: &mut RealMatrix,
        first: usize,
    ) -> Result<(), LinalgError> {
        if x.cols() == 0 {
            return Ok(());
        }
        let transposed = match form {
            SketchForm::Dense(t) => t,
            SketchForm::Count { bucket, sign } => {
                let (n, rows) = (x.rows(), out.rows());
                for c in 0..x.cols() {
                    let dst = &mut out.as_mut_slice()[(first + c) * rows..(first + c + 1) * rows];
                    dst.fill(0.);
                    for ((&v, &b), &s) in column(x, c)[..n].iter().zip(bucket).zip(sign) {
                        dst[b] += s * v;
                    }
                }
                return finite(&out.as_slice()[first * rows..(first + x.cols()) * rows]);
            }
        };
        let (n, rows) = (transposed.rows(), transposed.cols());
        gemm.multiply_slices_into(
            dft_numerics::GemmSite::Other,
            transposed.as_slice(),
            n,
            rows,
            n,
            RealOp::Transpose,
            x.as_slice(),
            x.cols(),
            n,
            out.as_mut().submatrix_mut(0, first, rows, x.cols()),
        )
    }
    /// Recompute A and B for the whole basis (start and non-contracting restart).
    fn refresh(
        &mut self,
        gemm: &mut RealGemmWorkspace,
        basis: &RealMatrix,
        products: &RealMatrix,
    ) -> Result<(), LinalgError> {
        Self::project(&self.form, gemm, basis, &mut self.a, 0)?;
        Self::project(&self.form, gemm, products, &mut self.b, 0)
    }
    /// Append the sketches of new basis columns starting at column `first`.
    fn append(
        &mut self,
        gemm: &mut RealGemmWorkspace,
        first: usize,
        extra: &RealMatrix,
        extra_products: &RealMatrix,
    ) -> Result<(), LinalgError> {
        Self::project(&self.form, gemm, extra, &mut self.a, first)?;
        Self::project(&self.form, gemm, extra_products, &mut self.b, first)
    }
    /// Contraction restart: the basis becomes V C[:, :bands], so A and B become
    /// the products formed by the last `estimate`.
    fn contract(&mut self, bands: usize) {
        let len = self.rows * bands;
        self.a.as_mut_slice()[..len].copy_from_slice(&self.ac.as_slice()[..len]);
        self.b.as_mut_slice()[..len].copy_from_slice(&self.bc.as_slice()[..len]);
    }
    /// est_j = || B c_j - theta_j A c_j || for the first `bands` Ritz pairs.
    fn estimate(
        &mut self,
        gemm: &mut RealGemmWorkspace,
        coefficients: &RealMatrix,
        m: usize,
        values: &[f64],
        est: &mut [f64],
    ) -> Result<(), LinalgError> {
        let rows = self.rows;
        let bands = values.len();
        for (side, out) in [(&self.a, &mut self.ac), (&self.b, &mut self.bc)] {
            gemm.multiply_slices_into(
                dft_numerics::GemmSite::Other,
                side.as_slice(),
                rows,
                m,
                rows,
                RealOp::Identity,
                coefficients.as_slice(),
                bands,
                m,
                out.as_mut(),
            )?;
        }
        for (j, (e, &theta)) in est.iter_mut().zip(values).enumerate() {
            let a = &self.ac.as_slice()[j * rows..(j + 1) * rows];
            let b = &self.bc.as_slice()[j * rows..(j + 1) * rows];
            let mut sums = [0.; 4];
            for (i, (&b, &a)) in b.iter().zip(a).enumerate() {
                let r = b - theta * a;
                sums[i % 4] += r * r;
            }
            *e = ((sums[0] + sums[1]) + (sums[2] + sums[3])).sqrt();
        }
        finite(est)
    }
}
enum SketchForm {
    /// S^T, N x s.
    Dense(RealMatrix),
    /// Bucket and sign of every coordinate.
    Count { bucket: Vec<usize>, sign: Vec<f64> },
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
    /// WIDE-DOT agrees with the four-lane dot to rounding for every tail length.
    #[test]
    fn wide_dot_matches_dot_to_rounding() {
        for len in [0, 1, 3, 7, 8, 9, 15, 16, 17, 1001, 11877] {
            let x: Vec<f64> = (0..len).map(|i| ((i + 1) as f64 * 0.37).sin()).collect();
            let y: Vec<f64> = (0..len).map(|i| ((i + 2) as f64 * 1.91).cos()).collect();
            let bound: f64 = x.iter().zip(&y).map(|(a, b)| (a * b).abs()).sum();
            let (a, b) = (dot(&x, &y), dot_wide(&x, &y));
            assert!(
                (a - b).abs() <= 4. * f64::EPSILON * bound,
                "{len}: {a:e} {b:e}"
            );
            assert_eq!(mgs_dot(false, &x, &y).to_bits(), a.to_bits());
            assert_eq!(mgs_dot(true, &x, &y).to_bits(), b.to_bits());
        }
    }
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
