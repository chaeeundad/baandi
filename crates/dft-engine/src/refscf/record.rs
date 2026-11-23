//! Plain in-memory diagnostics; formatting and file I/O belong to the CLI.
use super::broyden::Diagnostics;
use crate::energy::EnergyAudit;
/// Discrete work before/after each Ritz stopping decision; no timers or payloads.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CorrectionRecord {
    pub active_bands: Vec<usize>,
    pub basis_before: usize,
    pub h_columns: usize,
    pub reduced_order: usize,
    pub repairs: usize,
    pub dropped_corrections: usize,
    /// May include contraction after the stopping decision.
    pub restarts: usize,
}
#[derive(Clone, Debug)]
pub struct SolveRecord {
    pub ethr_ry: f64,
    pub first_step: bool,
    pub iterations: usize,
    pub active_bands: Vec<usize>,
    pub h_columns: usize,
    pub reduced_orders: Vec<usize>,
    pub restarts: usize,
    pub unconverged: usize,
    pub ritz_changes_ha: Vec<Vec<f64>>,
    pub ritz_values_ha: Vec<Vec<[f64; 2]>>,
    pub band_thresholds_ha: Vec<f64>,
    pub corrections: Vec<CorrectionRecord>,
    pub repairs: usize,
    pub dropped_corrections: usize,
}
#[derive(Clone, Debug)]
pub struct Record {
    pub iteration: usize,
    pub dr2_ry: f64,
    pub ethr_requested_ry: f64,
    pub ethr_ry: f64,
    pub first_iteration_redo: bool,
    pub beta: f64,
    pub mixing: Diagnostics,
    pub calls: Vec<SolveRecord>,
    /// Calls in Run::verification_calls belonging to this map, including retries.
    pub verification_calls: std::ops::Range<usize>,
    pub davidson_iterations: usize,
    pub h_columns: usize,
    pub map_s: f64,
    pub solve_s: f64,
    pub density_s: f64,
    pub potential_s: f64,
    pub energy_s: f64,
    pub density_transform_s: f64,
    pub metric_s: f64,
    pub mixing_s: f64,
    pub correction_s: f64,
    pub verification_s: f64,
    /// - integral (x_next-y) (V_H[x_next]+V_xc[x_next]) dr, in Ha.
    pub descf_ha: f64,
    pub variational_energy_ha: f64,
    pub energy: EnergyAudit,
    /// E_H and E_xc at the density used in the QE-style variational energy.
    pub variational_hartree_ha: f64,
    pub variational_xc_ha: f64,
}
