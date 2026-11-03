//! Blocked residual minimization with independent, short per-band DIIS histories.
//! See docs/research/rmm-diis-design.md for choices and limitations.
use super::*;
use dft_numerics::linalg::{HermitianWorkspace, gemm_into, orthonormalize_into};
use std::cell::Cell;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RmmFailureReason {
    RankLoss,
    Stagnation,
    IterationLimit,
}

/// Davidson phase associated with a terminal hybrid nonconvergence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HybridPhase {
    Warmup,
    Fallback,
}

#[derive(Debug, Error)]
pub enum RmmError {
    #[error(
        "Davidson {phase:?} failed after {operator_applications} total operator columns: {source}"
    )]
    Davidson {
        phase: HybridPhase,
        operator_applications: usize,
        /// Includes all per-band residuals for numerical nonconvergence.
        source: SolverError,
    },
    #[error(transparent)]
    Solver(#[from] SolverError),
    #[error("RMM refinement failed: {reason:?}")]
    Refinement {
        reason: RmmFailureReason,
        /// Last full-rank, rotated and truly evaluated state; absent on initial rank loss.
        last: Option<Box<EigenSolution>>,
    },
}
impl From<LinalgError> for RmmError {
    fn from(e: LinalgError) -> Self {
        Self::Solver(e.into())
    }
}

/// Experimental settings; none changes the Davidson default. Tolerances are
/// supplied separately for every band in ascending Ritz order (eigenvalue units).
#[derive(Debug, Clone)]
pub struct RmmDiis {
    pub max_sweeps: usize,
    pub history_size: usize,
    pub stagnation_sweeps: usize,
}
impl Default for RmmDiis {
    fn default() -> Self {
        Self {
            max_sweeps: 80,
            history_size: 4,
            stagnation_sweeps: 10,
        }
    }
}

fn checked_apply(
    op: &dyn LinearOperator,
    x: &ComplexMatrix,
    y: &mut ComplexMatrix,
) -> Result<(), SolverError> {
    op.apply_block(x, y)?;
    if y.rows() != x.rows() || y.cols() != x.cols() {
        return Err(LinalgError::Dimension.into());
    }
    y.validate_finite()?;
    Ok(())
}
fn dot(a: &[Complex64], b: &[Complex64]) -> Complex64 {
    a.iter().zip(b).map(|(a, b)| a.conj() * b).sum()
}
fn copy(a: &ComplexMatrix, b: &mut ComplexMatrix) {
    for j in 0..a.cols() {
        b.column_mut(j).copy_from_slice(a.column(j));
    }
}
fn validate_tolerances(n: usize, tolerances: &[f64]) -> Result<(), SolverError> {
    validate_bands(n, tolerances.len())?;
    if tolerances.iter().any(|v| !v.is_finite() || *v <= 0.0) {
        return Err(SolverError::InvalidParameter(
            "positive finite per-band tolerances required",
        ));
    }
    Ok(())
}

impl RmmDiis {
    pub fn solve(
        &self,
        operator: &dyn LinearOperator,
        initial: &ComplexMatrix,
        preconditioner: Option<&dyn Preconditioner>,
        tolerances: &[f64],
    ) -> Result<EigenSolution, RmmError> {
        let n = operator.dimension();
        let m = tolerances.len();
        validate_tolerances(n, tolerances)?;
        if initial.rows() != n || initial.cols() != m {
            return Err(LinalgError::Dimension.into());
        }
        if self.max_sweeps == 0
            || !(2..=16).contains(&self.history_size)
            || self.stagnation_sweeps == 0
        {
            return Err(SolverError::InvalidParameter(
                "positive sweep budgets and history size 2..=16 required",
            )
            .into());
        }
        initial.validate_finite()?;
        let h = self.history_size;
        let mh = m.checked_mul(h).ok_or(LinalgError::Dimension)?;
        // All refinement storage is allocated before the sweep/inner/band loops.
        let mut x = initial.clone();
        let mut hx = ComplexMatrix::zeros(n, m)?;
        let mut work = ComplexMatrix::zeros(n, m)?;
        let mut r = ComplexMatrix::zeros(n, m)?;
        let mut mixed = ComplexMatrix::zeros(n, m)?;
        let mut mixed_r = ComplexMatrix::zeros(n, m)?;
        let mut p = ComplexMatrix::zeros(n, m)?;
        let mut hp = ComplexMatrix::zeros(n, m)?;
        let mut history_x = ComplexMatrix::zeros(n, mh)?;
        let mut history_r = ComplexMatrix::zeros(n, mh)?;
        let mut small = ComplexMatrix::zeros(m, m)?;
        let mut rotation = ComplexMatrix::zeros(m, m)?;
        let mut rr = HermitianWorkspace::new(m)?;
        let mut metric = ComplexMatrix::zeros(h, h)?;
        let mut metric_vectors = ComplexMatrix::zeros(h, h)?;
        let mut metric_values = vec![0.0; h];
        let mut diis = HermitianWorkspace::new(h)?;
        let mut weights = vec![0.0; h];
        let mut scratch = vec![Complex64::default(); n];
        let mut values = vec![0.0; m];
        let mut residuals = vec![f64::INFINITY; m];
        let mut history = Vec::with_capacity(self.max_sweeps);
        let mut applications = 0;
        let mut best = f64::INFINITY;
        let mut stale = 0;
        // Preserve last certified state for duplicate/rank-loss recovery.
        let mut safe_x = ComplexMatrix::zeros(n, m)?;
        let mut safe_values = vec![0.0; m];
        let mut safe_residuals = vec![f64::INFINITY; m];
        for sweep in 1..=self.max_sweeps {
            let rank = orthonormalize_into(&x, &mut work, &mut scratch)?;
            if rank != m {
                let last = if sweep == 1 {
                    None
                } else {
                    Some(Box::new(EigenSolution {
                        eigenvalues: safe_values,
                        eigenvectors: safe_x,
                        residuals: safe_residuals,
                        iterations: sweep - 1,
                        history,
                        operator_applications: applications,
                    }))
                };
                return Err(RmmError::Refinement {
                    reason: RmmFailureReason::RankLoss,
                    last,
                });
            }
            std::mem::swap(&mut x, &mut work);
            checked_apply(operator, &x, &mut hx)?;
            applications += m;
            gemm_into(&x, Adjoint, &hx, Identity, &mut small)?;
            rr.solve(&small, &mut values, &mut rotation)?;
            gemm_into(&x, Identity, &rotation, Identity, &mut work)?;
            std::mem::swap(&mut x, &mut work);
            // True certification, not a recurrence residual.
            checked_apply(operator, &x, &mut hx)?;
            applications += m;
            for j in 0..m {
                for i in 0..n {
                    r[(i, j)] = hx[(i, j)] - x[(i, j)] * values[j];
                }
                residuals[j] = norm(r.column(j));
            }
            r.validate_finite()?;
            let maximum = max_residual(&residuals);
            history.push(maximum);
            let defect = orthogonality_defect(&x);
            if defect > 1e-10 {
                history.pop();
                let last = if sweep == 1 {
                    None
                } else {
                    Some(Box::new(EigenSolution {
                        eigenvalues: safe_values,
                        eigenvectors: safe_x,
                        residuals: safe_residuals,
                        iterations: sweep - 1,
                        history,
                        operator_applications: applications,
                    }))
                };
                return Err(RmmError::Refinement {
                    reason: RmmFailureReason::RankLoss,
                    last,
                });
            }
            if residuals.iter().zip(tolerances).all(|(r, t)| r <= t) {
                return Ok(EigenSolution {
                    eigenvalues: values,
                    eigenvectors: x,
                    residuals,
                    iterations: sweep,
                    history,
                    operator_applications: applications,
                });
            }
            // Relative improvement is used only to trigger fallback, never success.
            if maximum < best * 0.99 {
                best = maximum;
                stale = 0;
            } else {
                stale += 1;
            }
            if sweep == self.max_sweeps || stale >= self.stagnation_sweeps {
                return Err(RmmError::Refinement {
                    reason: if sweep == self.max_sweeps {
                        RmmFailureReason::IterationLimit
                    } else {
                        RmmFailureReason::Stagnation
                    },
                    last: Some(Box::new(EigenSolution {
                        eigenvalues: values,
                        eigenvectors: x,
                        residuals,
                        iterations: sweep,
                        history,
                        operator_applications: applications,
                    })),
                });
            }
            copy(&x, &mut safe_x);
            safe_values.copy_from_slice(&values);
            safe_residuals.copy_from_slice(&residuals);
            // A rotation changes band identities: histories are restarted every sweep.
            for inner in 0..h {
                if inner != 0 {
                    checked_apply(operator, &x, &mut hx)?;
                    applications += m;
                    for j in 0..m {
                        for i in 0..n {
                            r[(i, j)] = hx[(i, j)] - values[j] * x[(i, j)];
                        }
                    }
                }
                r.validate_finite()?;
                for j in 0..m {
                    history_x
                        .column_mut(j * h + inner)
                        .copy_from_slice(x.column(j));
                    history_r
                        .column_mut(j * h + inner)
                        .copy_from_slice(r.column(j));
                    weights.fill(0.0);
                    weights[inner] = 1.0;
                    if inner > 0 {
                        let active = inner + 1;
                        let mut scale: f64 = 0.0;
                        for a in 0..h {
                            for b in 0..h {
                                let value = if a < active && b < active {
                                    dot(history_r.column(j * h + a), history_r.column(j * h + b)).re
                                } else {
                                    0.0
                                };
                                metric[(a, b)] = value.into();
                                if a == b {
                                    scale = scale.max(value);
                                }
                            }
                        }
                        if scale.is_finite() && scale > f64::MIN_POSITIVE {
                            for a in 0..h {
                                for b in 0..h {
                                    metric[(a, b)] /= scale;
                                }
                            }
                            // Padded inactive modes are excluded from the constraint.
                            for a in active..h {
                                metric[(a, a)] = 1.0.into();
                            }
                            diis.solve(&metric, &mut metric_values, &mut metric_vectors)?;
                            weights.fill(0.0);
                            for (k, &e) in metric_values.iter().enumerate() {
                                let overlap: f64 =
                                    (0..active).map(|a| metric_vectors[(a, k)].re).sum();
                                for a in 0..active {
                                    weights[a] +=
                                        metric_vectors[(a, k)].re * overlap / e.max(1e-12);
                                }
                            }
                            let sum: f64 = weights[..active].iter().sum();
                            if sum.is_finite() && sum.abs() > f64::MIN_POSITIVE {
                                for w in &mut weights[..active] {
                                    *w /= sum;
                                }
                            }
                            if !sum.is_finite()
                                || sum.abs() <= f64::MIN_POSITIVE
                                || weights.iter().any(|w| !w.is_finite() || w.abs() > 10.0)
                            {
                                weights.fill(0.0);
                                weights[inner] = 1.0;
                            }
                        }
                    }
                    for i in 0..n {
                        mixed[(i, j)] = (0..=inner)
                            .map(|a| history_x[(i, j * h + a)] * weights[a])
                            .sum();
                        mixed_r[(i, j)] = (0..=inner)
                            .map(|a| history_r[(i, j * h + a)] * weights[a])
                            .sum();
                    }
                }
                if let Some(pre) = preconditioner {
                    pre.apply(&values, &mixed_r, &mut p)?;
                    if p.rows() != n || p.cols() != m {
                        return Err(LinalgError::Dimension.into());
                    }
                    p.validate_finite()?;
                } else {
                    for j in 0..m {
                        for i in 0..n {
                            p[(i, j)] = -mixed_r[(i, j)];
                        }
                    }
                }
                checked_apply(operator, &p, &mut hp)?;
                applications += m;
                for j in 0..m {
                    for i in 0..n {
                        hp[(i, j)] -= values[j] * p[(i, j)];
                    }
                    let qnorm = norm(hp.column(j));
                    let alpha = if qnorm > f64::MIN_POSITIVE {
                        -dot(hp.column(j), mixed_r.column(j)).re / qnorm / qnorm
                    } else {
                        0.0
                    };
                    if !alpha.is_finite() {
                        return Err(LinalgError::NonFinite.into());
                    }
                    for i in 0..n {
                        work[(i, j)] = mixed[(i, j)] + alpha * p[(i, j)];
                    }
                    let length = norm(work.column(j));
                    if length > 1e-12 && length.is_finite() {
                        for i in 0..n {
                            scratch[i] = (mixed_r[(i, j)] + alpha * hp[(i, j)]) / length;
                        }
                        // Reject a deteriorating DIIS extrapolation; next sweep resets it.
                        if norm(&scratch) <= norm(r.column(j)) * (1.0 + 1e-12) {
                            for i in 0..n {
                                x[(i, j)] = work[(i, j)] / length;
                            }
                        }
                    }
                }
                x.validate_finite()?;
            }
        }
        unreachable!("positive sweep count validated")
    }
}

/// Tolerance-based Davidson warm-up, RMM refinement, then at most one Davidson
/// fallback. Warm-up/fallback retain the existing allocation-owning Davidson.
#[derive(Debug, Clone)]
pub struct HybridRmm {
    pub refinement: RmmDiis,
    pub warmup_tolerance: f64,
}
impl Default for HybridRmm {
    fn default() -> Self {
        Self {
            refinement: RmmDiis::default(),
            warmup_tolerance: 1e-3,
        }
    }
}
#[derive(Debug)]
pub struct HybridSolution {
    pub solution: EigenSolution,
    pub warmed_up: bool,
    pub fallback: Option<RmmFailureReason>,
    pub rmm_sweeps: usize,
}
struct Counted<'a> {
    op: &'a dyn LinearOperator,
    count: Cell<usize>,
}
impl LinearOperator for Counted<'_> {
    fn dimension(&self) -> usize {
        self.op.dimension()
    }
    fn apply(&self, x: &[Complex64], y: &mut [Complex64]) -> Result<(), SolverError> {
        self.count.set(self.count.get() + 1);
        self.op.apply(x, y)
    }
    fn apply_block(&self, x: &ComplexMatrix, y: &mut ComplexMatrix) -> Result<(), SolverError> {
        self.count.set(self.count.get() + x.cols());
        self.op.apply_block(x, y)
    }
}
impl HybridRmm {
    pub fn solve(
        &self,
        davidson: &Davidson,
        op: &dyn LinearOperator,
        initial: Option<&ComplexMatrix>,
        pre: Option<&dyn Preconditioner>,
        tolerances: &[f64],
    ) -> Result<HybridSolution, RmmError> {
        validate_tolerances(op.dimension(), tolerances)?;
        if davidson.bands != tolerances.len()
            || !self.warmup_tolerance.is_finite()
            || self.warmup_tolerance <= 0.0
        {
            return Err(SolverError::InvalidParameter(
                "hybrid bands and positive warm-up tolerance required",
            )
            .into());
        }
        let counted = Counted {
            op,
            count: Cell::new(0),
        };
        let mut strict = davidson.clone();
        strict.tolerance = tolerances.iter().copied().fold(f64::INFINITY, f64::min);
        let mut warm = strict.clone();
        warm.tolerance = strict.tolerance.max(self.warmup_tolerance);
        // Orthogonality alone does not imply a useful starting eigenspace.
        // Screen true Rayleigh residuals before bypassing Davidson warm-up.
        let mut warmed_up = true;
        if let Some(v) = initial {
            if v.rows() != op.dimension() {
                return Err(LinalgError::Dimension.into());
            }
            v.validate_finite()?;
            if v.cols() == davidson.bands && orthogonality_defect(v) <= 1e-10 {
                let mut products = ComplexMatrix::zeros(v.rows(), v.cols())?;
                checked_apply(&counted, v, &mut products)?;
                let mut maximum: f64 = 0.0;
                for j in 0..v.cols() {
                    let value = dot(v.column(j), products.column(j)).re;
                    for i in 0..v.rows() {
                        products[(i, j)] -= value * v[(i, j)];
                    }
                    maximum = maximum.max(norm(products.column(j)));
                }
                products.validate_finite()?;
                warmed_up = maximum > warm.tolerance;
            }
        }
        let warm_solution = if warmed_up {
            Some(
                warm.solve(&counted, initial, pre)
                    .map_err(|e| hybrid_failure(e, HybridPhase::Warmup, counted.count.get()))?,
            )
        } else {
            None
        };
        let guess = warm_solution
            .as_ref()
            .map(|s| &s.eigenvectors)
            .or(initial)
            .expect("warm-up supplies guess");
        let mut prior_history = Vec::new();
        if let Some(s) = &warm_solution {
            prior_history.extend_from_slice(&s.history);
        }
        let attempt = self.refinement.solve(&counted, guess, pre, tolerances);
        let (mut solution, fallback, rmm_sweeps) = match attempt {
            Ok(s) => {
                let sweeps = s.iterations;
                (s, None, sweeps)
            }
            Err(RmmError::Refinement { reason, last }) => {
                let sweeps = last.as_ref().map_or(0, |s| s.iterations);
                if let Some(s) = &last {
                    prior_history.extend_from_slice(&s.history);
                }
                let recovery = last.as_ref().map(|s| &s.eigenvectors).unwrap_or(guess);
                (
                    strict.solve(&counted, Some(recovery), pre).map_err(|e| {
                        hybrid_failure(e, HybridPhase::Fallback, counted.count.get())
                    })?,
                    Some(reason),
                    sweeps,
                )
            }
            Err(e) => return Err(e),
        };
        prior_history.extend_from_slice(&solution.history);
        solution.history = prior_history;
        solution.iterations = solution.history.len();
        solution.operator_applications = counted.count.get();
        Ok(HybridSolution {
            solution,
            warmed_up,
            fallback,
            rmm_sweeps,
        })
    }
}

/// Additive caller-controlled selection. The default preserves Davidson exactly.
#[derive(Debug, Clone, Default)]
pub enum SolverChoice {
    #[default]
    Davidson,
    Hybrid(HybridRmm),
}
impl SolverChoice {
    pub fn solve(
        &self,
        davidson: &Davidson,
        operator: &dyn LinearOperator,
        initial: Option<&ComplexMatrix>,
        preconditioner: Option<&dyn Preconditioner>,
    ) -> Result<EigenSolution, SolverError> {
        match self {
            Self::Davidson => davidson.solve(operator, initial, preconditioner),
            Self::Hybrid(policy) => policy
                .solve(
                    davidson,
                    operator,
                    initial,
                    preconditioner,
                    &vec![davidson.tolerance; davidson.bands],
                )
                .map(|s| s.solution)
                .map_err(|e| match e {
                    RmmError::Solver(e) | RmmError::Davidson { source: e, .. } => e,
                    RmmError::Refinement { last: Some(s), .. } => SolverError::NotConverged {
                        iterations: s.iterations,
                        max_residual: max_residual(&s.residuals),
                        residuals: s.residuals,
                    },
                    RmmError::Refinement { last: None, .. } => {
                        SolverError::InvalidParameter("no independent RMM initial state")
                    }
                }),
        }
    }
}

fn hybrid_failure(
    source: SolverError,
    phase: HybridPhase,
    operator_applications: usize,
) -> RmmError {
    match source {
        source @ SolverError::NotConverged { .. } => RmmError::Davidson {
            phase,
            operator_applications,
            source,
        },
        e => RmmError::Solver(e),
    }
}
