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
pub use functional::Functional;

mod band_density;

mod nscf;
pub use nscf::*;

pub mod events;
pub mod state_view;
pub use events::{ExecutionConfig, ExecutionMode};
pub use session::{CompletionMode, ScfSession, StepReport, StepStatus};

pub mod diagnostics;
