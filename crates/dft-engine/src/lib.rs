//! Weighted Bloch, spin-unpolarized insulating LDA-PZ SCF and energy audit.
pub mod energy;
pub mod scf;
pub use scf::*;

/// Stable input facade for clients; numerical implementations remain private to
/// their owning crates. All faer calls use explicit `Par::Seq`, and FFTs are serial.
pub mod input {
    pub use dft_geometry::{Atom, Cell};
    pub use dft_upf::{Pseudopotential, parse_upf};
}
/// The backend explicitly selects sequential execution (no global thread pool).
pub const EXECUTION_THREADS: usize = 1;

pub mod functional;
pub use functional::{Functional, Pw92Prefactor};

mod band_density;

mod nscf;
pub use nscf::*;

pub mod events;
pub mod policy;
pub mod runtime;
pub mod state_view;
pub use events::{ExecutionConfig, ExecutionMode};
pub use session::{CompletionMode, ScfSession, StepReport, StepStatus};

pub mod diagnostics;
pub mod phase_clock;

pub mod checkpoint;

mod spectral;
pub use spectral::{SpectralEvent, SpectralGuard, SpectralReport};

pub mod forces;
pub use forces::{ElectronicState, ForceResult, ForceState, forces};

mod force_identity;

pub mod stress;
pub use stress::{StressResult, StressState, stress};

pub use dft_solvers::locking::Locking;

pub mod gamma_real;
pub mod m2;
pub mod parallel;
pub use dft_solvers::real_davidson::{
    StageSeconds as StageSecondsReport, WIDTH_LABELS as DAVIDSON_WIDTH_LABELS, algebra_counters,
    b1_counters, b5_counters, ortho_counters, second_vector_columns, stage_seconds,
    stage_seconds_by_run_stage,
};
pub use gamma_real::{GammaReal, GammaRealDiagnostics};

/// Research radial conventions; default variants retain the legacy path.
pub use dft_upf::radial::{LocalMesh, RadialEvenRule};
