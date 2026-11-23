//! Plain in-memory diagnostics; formatting and file I/O belong to the CLI.
use super::broyden::Diagnostics;
#[derive(Clone, Debug)]
pub struct Record {
    pub iteration: usize,
    pub dr2_ry: f64,
    pub ethr_ry: f64,
    pub first_iteration_redo: bool,
    pub beta: f64,
    pub mixing: Diagnostics,
    pub davidson_iterations: usize,
    pub h_columns: usize,
    pub map_s: f64,
    pub density_transform_s: f64,
    pub metric_s: f64,
    pub mixing_s: f64,
    pub correction_s: f64,
    /// - integral (x_next-y) (V_H[x_next]+V_xc[x_next]) dr, in Ha.
    pub descf_ha: f64,
    pub variational_energy_ha: f64,
}
