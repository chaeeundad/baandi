//! Shared runtime lifecycle and policy payload contracts.
use super::Mixing;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RuntimeState {
    Created,
    Initialized,
    Running,
    PauseRequested,
    Paused,
    Restarting,
    Converged,
    Exhausted,
    Failed,
    Cancelled,
}
impl RuntimeState {
    pub fn terminal(self) -> bool {
        matches!(
            self,
            Self::Converged | Self::Exhausted | Self::Failed | Self::Cancelled
        )
    }
}
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

pub(crate) fn mixing_key(mixing: &Mixing) -> (u8, u64, usize, u64) {
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
