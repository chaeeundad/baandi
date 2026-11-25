//! Shared SCF failures.
use super::IterationRecord;
use dft_physics::{
    data::{FftError, LinalgError},
    hamiltonian,
};
use std::fmt;

#[derive(Debug)]
pub enum ScfError {
    InvalidInput(&'static str),
    UnsupportedFunctional(String),
    OddElectronCount(f64),
    /// The sampled indirect gap does not certify fixed integer occupations.
    NonInsulating {
        gap_hartree: f64,
        required_hartree: f64,
    },
    InsufficientBands {
        bands: usize,
        occupied: usize,
    },
    Numerical(String),
    NonFinite,
    /// Bounded observation storage cannot reserve another complete step.
    TelemetryExhausted,
    NotConverged {
        history: Vec<IterationRecord>,
    },
    Diverged {
        reason: String,
        history: Vec<IterationRecord>,
    },
}
impl fmt::Display for ScfError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Self::UnsupportedFunctional(label) = self {
            write!(
                f,
                "unsupported UPF functional {label:?}: requires SLA PZ or SLA PZ NOGX NOGC"
            )
        } else {
            write!(f, "{self:?}")
        }
    }
}
impl std::error::Error for ScfError {}
macro_rules! numerical_error { ($($t:ty),*) => { $(impl From<$t> for ScfError { fn from(e:$t)->Self { Self::Numerical(e.to_string()) } })* }; }
numerical_error!(dft_basis::BasisError, FftError);
impl From<dft_physics::PhysicsError> for ScfError {
    fn from(e: dft_physics::PhysicsError) -> Self {
        match e {
            dft_physics::PhysicsError::Unrepresentable => Self::NonFinite,
            dft_physics::PhysicsError::Density(n) if !n.is_finite() => Self::NonFinite,
            _ => Self::Numerical(e.to_string()),
        }
    }
}
impl From<hamiltonian::HamiltonianError> for ScfError {
    fn from(e: hamiltonian::HamiltonianError) -> Self {
        match e {
            hamiltonian::HamiltonianError::NonFinite => Self::NonFinite,
            _ => Self::Numerical(e.to_string()),
        }
    }
}
impl From<dft_solvers::SolverError> for ScfError {
    fn from(e: dft_solvers::SolverError) -> Self {
        match e {
            dft_solvers::SolverError::Linalg(LinalgError::NonFinite) => Self::NonFinite,
            _ => Self::Numerical(e.to_string()),
        }
    }
}
