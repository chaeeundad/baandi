//! Weighted Bloch, spin-unpolarized insulating LDA-PZ SCF and energy audit.
pub mod energy;
pub mod model;
#[cfg(feature = "strict-oracle")]
pub mod oracle;
#[cfg(feature = "strict-oracle")]
pub mod scf;
#[cfg(not(feature = "strict-oracle"))]
pub use model::*;
#[cfg(feature = "strict-oracle")]
pub use scf::*;

/// Stable input facade for clients; numerical implementations remain private to
/// their owning crates. All faer calls use explicit `Par::Seq`, and FFTs are serial.
pub mod input {
    pub use dft_geometry::kpoints;
    pub use dft_geometry::{Atom, Cell};
    pub use dft_upf::{Pseudopotential, parse_upf};
}
/// The backend explicitly selects sequential execution (no global thread pool).
pub const EXECUTION_THREADS: usize = 1;

pub mod functional;
pub use functional::{Functional, Pw92Prefactor};

#[cfg(feature = "strict-oracle")]
use oracle::band_density;

#[cfg(feature = "strict-oracle")]
pub use nscf::*;
#[cfg(feature = "strict-oracle")]
use oracle::nscf;

#[cfg(feature = "strict-oracle")]
pub use events::{ExecutionConfig, ExecutionMode};
#[cfg(feature = "strict-oracle")]
pub use oracle::events;
#[cfg(feature = "strict-oracle")]
pub use oracle::policy;
#[cfg(feature = "strict-oracle")]
pub use oracle::runtime;
#[cfg(feature = "strict-oracle")]
pub use oracle::state_view;
#[cfg(feature = "strict-oracle")]
pub use session::{CompletionMode, ScfSession, StepReport, StepStatus, assess_fixed_point};

#[cfg(feature = "strict-oracle")]
pub use oracle::diagnostics;
pub mod phase_clock;
// Optional exclusive phase profile and SCF snapshots.
pub mod profile;
#[cfg(feature = "strict-oracle")]
pub use oracle::snapshot;

#[cfg(feature = "strict-oracle")]
pub use oracle::checkpoint;

mod spectral;
pub use spectral::{SpectralEvent, SpectralGuard, SpectralReport};

pub mod forces;
pub use forces::{ElectronicState, ForceResult, ForceState, forces};

mod force_identity;

pub mod stress;
pub use stress::{PulayAssessmentMode, StressResult, StressState, stress};

pub use dft_solvers::locking::Locking;

#[cfg(feature = "strict-oracle")]
pub use oracle::gamma_real;
pub mod grid_dist;
#[cfg(feature = "strict-oracle")]
pub use oracle::m2;
#[cfg(feature = "strict-oracle")]
pub use oracle::mix_dist;
pub mod parallel;
pub use dft_solvers::real_davidson::{
    StageSeconds as StageSecondsReport, WIDTH_LABELS as DAVIDSON_WIDTH_LABELS, algebra_counters,
    b1_counters, b5_counters, flop_counters, lean_counters, ortho_counters, second_vector_columns,
    stage_seconds, stage_seconds_by_run_stage,
};
#[cfg(feature = "strict-oracle")]
pub use gamma_real::{GammaReal, GammaRealDiagnostics, nl_audit_counters};

/// Research radial conventions; default variants retain the legacy path.
pub use dft_upf::radial::{LocalMesh, RadialEvenRule};

mod atomic_start;
pub mod refscf;
/// Final-cutoff, explicitly selected distributed reference initialization.
pub mod refstart;
/// R2 dielectric response (chi0, head/wings, epsilon_M) on NSCF bands.
pub mod response;

#[cfg(test)]
#[allow(dead_code)]
#[path = "../tests/common/mod.rs"]
mod pool_fixture;
