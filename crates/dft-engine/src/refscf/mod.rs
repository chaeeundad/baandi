//! Independent reciprocal-density SCF (issue #25, patch 1).
//!
//! Conventions: docs/research/qe-study/qe-density-mixing-spec.md, an own-words
//! specification, never QE source. Method: D. D. Johnson, Phys. Rev. B 38,
//! 12807 (1988), doi:10.1103/PhysRevB.38.12807; difference-space Anderson/Pulay
//! mixing, P. Pulay, Chem. Phys. Lett. 73, 393 (1980); G. Kresse and
//! J. Furthmueller, Phys. Rev. B 54, 11169 (1996).
//! No legacy session or ladder control flow is used here. The map adapter still
//! reuses the Hamiltonian/density kernels and the new reference Davidson.
pub mod api;
pub mod broyden;
pub mod checkpoint;
pub mod delivered;
pub mod density;
pub mod driver;
pub mod electronic;
mod electronic_map;
pub mod k_driver;
mod k_kept;
mod k_map;
pub mod kept_hx;
pub mod ladder;
pub mod lazy;
mod map;
pub mod memory;
pub mod metric;
pub mod nscf;
pub mod profile;
pub mod prolong;
pub mod qe_reference;
pub mod record;
pub mod runtime;
pub mod stepper;
pub mod symmetry;

use std::fmt;
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    Invalid(String),
    Numerical(String),
    Singular {
        iteration: usize,
        min: f64,
        max: f64,
        condition: f64,
    },
}
impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(s) | Self::Numerical(s) => write!(f, "reference SCF: {s}"),
            Self::Singular {
                iteration,
                min,
                max,
                condition,
            } => write!(
                f,
                "reference SCF: singular history A at mixing iteration {iteration}: eigenvalues [{min:e}, {max:e}] Ry, condition {condition:e}; no regularization applied"
            ),
        }
    }
}
impl std::error::Error for Error {}
impl From<Error> for crate::ScfError {
    fn from(e: Error) -> Self {
        Self::Numerical(e.to_string())
    }
}
pub(crate) fn invalid(s: impl Into<String>) -> Error {
    Error::Invalid(s.into())
}
pub(crate) fn numerical(e: impl fmt::Display) -> Error {
    Error::Numerical(e.to_string())
}
