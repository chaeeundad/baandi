//! Constant-size immutable projection; large arrays require explicit accessors.
use crate::{ScfOptions, energy::EnergyAudit, events::ExecutionMode, session::StepStatus};
#[derive(Debug)]
pub struct StateSnapshot<'a> {
    pub schema_version: &'static str,
    pub calc_id: &'a str,
    pub mode: ExecutionMode,
    pub sequence: u64,
    pub state_version: u64,
    pub iteration: usize,
    pub status: StepStatus,
    /// May exceed `iteration` after a failed partial step; never implies completion.
    pub record_iteration: Option<usize>,
    pub energy_change_hartree: Option<f64>,
    pub requested_eigen_tolerance_hartree: Option<f64>,
    pub input_charge_error_relative: Option<f64>,
    /// None before the first evaluated record. All energy terms are in Ha.
    pub energy: Option<&'a EnergyAudit>,
    /// Normalized density L2 residual, Bohr^-3/2.
    pub density_residual: Option<f64>,
    pub max_eigen_residual_hartree: Option<f64>,
    pub charge_error_relative: Option<f64>,
    pub options: &'a ScfOptions,
}
