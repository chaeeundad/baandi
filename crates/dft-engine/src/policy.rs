//! Pure preparation of boundary policy transactions. No numerical state is mutated.
//!
//! Schedule and divergence interventions are explicitly rejected by the P11-COMMON
//! allow-list. Model fields (XC, species, charge, geometry, grid, bands, threads,
//! seed and pseudopotentials) cannot be expressed as a policy change.
use crate::{Mixing, ScfOptions};
use dft_solvers::SolverChoice;

/// Select the existing implementations with their existing defaults. Arbitrary
/// hybrid tuning parameters are not a runtime intervention.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RuntimeEigensolver {
    Davidson,
    Hybrid,
}

#[derive(Clone, Debug)]
pub enum PolicyChange {
    Mixing(Mixing),
    MixingBeta(f64),
    MixingHistory(usize),
    ResetHistory,
    Eigensolver(RuntimeEigensolver),
    DensityTolerance(f64),
    EnergyToleranceHartree(f64),
    EigenToleranceHartree(f64),
    MaxIterations(usize),
    ScheduleFlags { residual_eigen_schedule: bool },
    DivergenceParameters { factor: f64, patience: usize },
}

impl PartialEq for PolicyChange {
    fn eq(&self, other: &Self) -> bool {
        use PolicyChange::*;
        match (self, other) {
            (Mixing(a), Mixing(b)) => mixing_key(a) == mixing_key(b),
            (MixingBeta(a), MixingBeta(b))
            | (DensityTolerance(a), DensityTolerance(b))
            | (EnergyToleranceHartree(a), EnergyToleranceHartree(b))
            | (EigenToleranceHartree(a), EigenToleranceHartree(b)) => a.to_bits() == b.to_bits(),
            (MixingHistory(a), MixingHistory(b)) | (MaxIterations(a), MaxIterations(b)) => a == b,
            (ResetHistory, ResetHistory) => true,
            (Eigensolver(a), Eigensolver(b)) => a == b,
            (
                ScheduleFlags {
                    residual_eigen_schedule: a,
                },
                ScheduleFlags {
                    residual_eigen_schedule: b,
                },
            ) => a == b,
            (
                DivergenceParameters {
                    factor: a,
                    patience: p,
                },
                DivergenceParameters {
                    factor: b,
                    patience: q,
                },
            ) => a.to_bits() == b.to_bits() && p == q,
            _ => false,
        }
    }
}

/// Immutable deployment metadata, not new scientific defaults. Lower bounds
/// prevent unbounded tightening; the launch and certificate ceilings separately
/// prevent loosening. Without evidence-backed floors tolerance changes reject.
#[derive(Clone, Copy, Debug)]
pub struct ToleranceFloors {
    pub density_residual: f64,
    pub energy_hartree: f64,
    pub eigen_hartree: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PolicyError {
    Empty,
    NoChange,
    TooManyChanges,
    DuplicateField,
    NonFiniteOrNonPositive,
    InvalidBeta,
    HistoryQuota,
    HistoryRequiresPulay,
    IterationQuota,
    BudgetMustIncrease,
    MissingToleranceFloors,
    InvalidToleranceFloors,
    BelowToleranceFloor,
    ToleranceRelaxation,
    CertificateThreshold,
    RestrictedSchedule,
    RestrictedDivergence,
}

/// Prepared options are deliberately read-only. A session integration must
/// prepare a replacement mixer before committing these options and the audit.
#[derive(Clone, Debug)]
pub struct PreparedPolicy {
    options: ScfOptions,
    reset_history: bool,
}
impl PreparedPolicy {
    pub fn options(&self) -> &ScfOptions {
        &self.options
    }
    pub fn resets_history(&self) -> bool {
        self.reset_history
    }
}

fn positive(v: f64) -> bool {
    v.is_finite() && v > 0.0
}
fn validate_mixing(mixing: &Mixing) -> Result<(), PolicyError> {
    let (beta, history, q0) = match *mixing {
        Mixing::Linear { beta } => (beta, None, None),
        Mixing::Pulay { beta, history } => (beta, Some(history), None),
        Mixing::KerkerPulay {
            beta,
            history,
            q0_inverse_bohr,
        } => (beta, Some(history), Some(q0_inverse_bohr)),
    };
    if !positive(beta) || beta > 1.0 {
        return Err(PolicyError::InvalidBeta);
    }
    if history.is_some_and(|n| n == 0 || n > 64) {
        return Err(PolicyError::HistoryQuota);
    }
    if q0.is_some_and(|q| !positive(q)) {
        return Err(PolicyError::NonFiniteOrNonPositive);
    }
    Ok(())
}
fn mixing_key(mixing: &Mixing) -> (u8, u64, usize, u64) {
    match *mixing {
        Mixing::Linear { beta } => (0, beta.to_bits(), 0, 0),
        Mixing::Pulay { beta, history } => (1, beta.to_bits(), history, 0),
        Mixing::KerkerPulay {
            beta,
            history,
            q0_inverse_bohr,
        } => (2, beta.to_bits(), history, q0_inverse_bohr.to_bits()),
    }
}
fn solver_equal(a: &SolverChoice, b: &SolverChoice) -> bool {
    match (a, b) {
        (SolverChoice::Davidson, SolverChoice::Davidson) => true,
        (SolverChoice::Hybrid(a), SolverChoice::Hybrid(b)) => {
            a.warmup_tolerance.to_bits() == b.warmup_tolerance.to_bits()
                && a.refinement.max_sweeps == b.refinement.max_sweeps
                && a.refinement.history_size == b.refinement.history_size
                && a.refinement.stagnation_sweeps == b.refinement.stagnation_sweeps
        }
        _ => false,
    }
}
fn tighten(new: f64, old: f64, floor: f64, ceiling: f64) -> Result<(), PolicyError> {
    if !positive(new) {
        return Err(PolicyError::NonFiniteOrNonPositive);
    }
    if !positive(floor) || floor > old {
        return Err(PolicyError::InvalidToleranceFloors);
    }
    if new > old {
        return Err(PolicyError::ToleranceRelaxation);
    }
    if new > ceiling {
        return Err(PolicyError::CertificateThreshold);
    }
    if new < floor {
        return Err(PolicyError::BelowToleranceFloor);
    }
    Ok(())
}

/// Validate a whole patch against the current committed boundary. Duplicate or
/// overlapping mixer assignments reject, avoiding order-dependent patches.
/// There is no unsafe-relaxed escape hatch. Tolerance ceilings use the existing
/// certificate defaults, and all live tolerance changes must also tighten.
pub fn prepare_policy(
    current: &ScfOptions,
    iteration_done: usize,
    floors: Option<ToleranceFloors>,
    changes: &[PolicyChange],
) -> Result<PreparedPolicy, PolicyError> {
    if changes.is_empty() {
        return Err(PolicyError::Empty);
    }
    if changes.len() > 10 {
        return Err(PolicyError::TooManyChanges);
    }
    let mut next = current.clone();
    let mut fields = 0u16;
    let mut reset = false;
    let defaults = ScfOptions::default();
    for change in changes {
        let field = match change {
            PolicyChange::Mixing(_)
            | PolicyChange::MixingBeta(_)
            | PolicyChange::MixingHistory(_) => 0,
            PolicyChange::ResetHistory => 1,
            PolicyChange::Eigensolver(_) => 2,
            PolicyChange::DensityTolerance(_) => 3,
            PolicyChange::EnergyToleranceHartree(_) => 4,
            PolicyChange::EigenToleranceHartree(_) => 5,
            PolicyChange::MaxIterations(_) => 6,
            PolicyChange::ScheduleFlags { .. } => return Err(PolicyError::RestrictedSchedule),
            PolicyChange::DivergenceParameters { .. } => {
                return Err(PolicyError::RestrictedDivergence);
            }
        };
        if fields & (1 << field) != 0 {
            return Err(PolicyError::DuplicateField);
        }
        fields |= 1 << field;
        match *change {
            PolicyChange::Mixing(ref m) => next.mixing = m.clone(),
            PolicyChange::MixingBeta(beta) => match &mut next.mixing {
                Mixing::Linear { beta: b }
                | Mixing::Pulay { beta: b, .. }
                | Mixing::KerkerPulay { beta: b, .. } => *b = beta,
            },
            PolicyChange::MixingHistory(n) => match &mut next.mixing {
                Mixing::Linear { .. } => return Err(PolicyError::HistoryRequiresPulay),
                Mixing::Pulay { history, .. } | Mixing::KerkerPulay { history, .. } => *history = n,
            },
            PolicyChange::ResetHistory => reset = true,
            PolicyChange::Eigensolver(s) => {
                next.eigensolver = match s {
                    RuntimeEigensolver::Davidson => SolverChoice::Davidson,
                    RuntimeEigensolver::Hybrid => SolverChoice::Hybrid(Default::default()),
                }
            }
            PolicyChange::DensityTolerance(v) => {
                tighten(
                    v,
                    current.density_tolerance,
                    floors
                        .ok_or(PolicyError::MissingToleranceFloors)?
                        .density_residual,
                    defaults.density_tolerance,
                )?;
                next.density_tolerance = v;
            }
            PolicyChange::EnergyToleranceHartree(v) => {
                tighten(
                    v,
                    current.energy_tolerance,
                    floors
                        .ok_or(PolicyError::MissingToleranceFloors)?
                        .energy_hartree,
                    defaults.energy_tolerance,
                )?;
                next.energy_tolerance = v;
            }
            PolicyChange::EigenToleranceHartree(v) => {
                tighten(
                    v,
                    current.eigen_tolerance,
                    floors
                        .ok_or(PolicyError::MissingToleranceFloors)?
                        .eigen_hartree,
                    defaults.eigen_tolerance,
                )?;
                next.eigen_tolerance = v;
            }
            PolicyChange::MaxIterations(n) => {
                if n > 10000 {
                    return Err(PolicyError::IterationQuota);
                }
                if n <= current.max_iterations || n <= iteration_done {
                    return Err(PolicyError::BudgetMustIncrease);
                }
                next.max_iterations = n;
            }
            PolicyChange::ScheduleFlags { .. } | PolicyChange::DivergenceParameters { .. } => {
                unreachable!()
            }
        }
    }
    validate_mixing(&next.mixing)?;
    let mixing_changed = mixing_key(&next.mixing) != mixing_key(&current.mixing);
    if !reset
        && !mixing_changed
        && solver_equal(&next.eigensolver, &current.eigensolver)
        && next.density_tolerance == current.density_tolerance
        && next.energy_tolerance == current.energy_tolerance
        && next.eigen_tolerance == current.eigen_tolerance
        && next.max_iterations == current.max_iterations
    {
        return Err(PolicyError::NoChange);
    }
    Ok(PreparedPolicy {
        options: next,
        reset_history: reset || mixing_changed,
    })
}
