//! Shared species and geometry-verified symmetry inputs.
use super::ScfError;
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use dft_upf::Pseudopotential;

/// Species labels resolve explicitly against Atom::species; the PP must declare SLA PZ or SLA PZ NOGX NOGC.
pub struct SpeciesPotential<'a> {
    pub label: &'a str,
    pub potential: &'a Pseudopotential,
}

/// Immutable geometry-verified symmetry preparation shared with input adapters.
/// Private fields prevent callers from supplying an unverified orbit map.
#[derive(Debug)]
pub struct PreparedSymmetry {
    pub(crate) lattice: [[f64; 3]; 3],
    pub(crate) atoms: Vec<Atom>,
    full: Vec<KPoint>,
    tolerance_bohr: f64,
    pub(crate) time_reversal: bool,
    map: dft_geometry::symmetry::KOrbitMap,
    crystal_operations: Vec<dft_geometry::symmetry::SymmetryOperation>,
    /// Explicit FFT grid the group was restricted to (GS1C), if any.
    grid: Option<[usize; 3]>,
    /// Detected operations dropped because they do not permute `grid`.
    dropped_for_grid: usize,
}
impl PreparedSymmetry {
    pub fn new(
        cell: &Cell,
        atoms: &[Atom],
        full: &[KPoint],
        tolerance_bohr: f64,
        time_reversal: bool,
    ) -> Result<Self, ScfError> {
        Self::new_for_grid(cell, atoms, full, tolerance_bohr, time_reversal, None)
    }
    /// GS1C: with an explicit FFT grid keep only the largest subgroup whose
    /// operations permute that grid (QE-like), instead of rejecting the grid.
    pub fn new_for_grid(
        cell: &Cell,
        atoms: &[Atom],
        full: &[KPoint],
        tolerance_bohr: f64,
        time_reversal: bool,
        grid: Option<[usize; 3]>,
    ) -> Result<Self, ScfError> {
        let detected = dft_geometry::symmetry::detect(cell, atoms, tolerance_bohr)
            .map_err(|e| ScfError::Numerical(e.to_string()))?;
        let crystal_operations = match grid {
            None => detected.clone(),
            Some(dims) => {
                dft_geometry::symmetry::grid_subgroup(cell, &detected, dims, tolerance_bohr)
            }
        };
        let dropped_for_grid = detected.len() - crystal_operations.len();
        let map =
            dft_geometry::symmetry::reduce_mesh(full, &crystal_operations, time_reversal, false)
                .map_err(|e| ScfError::Numerical(e.to_string()))?;
        Ok(Self {
            lattice: cell.lattice_columns_bohr(),
            atoms: atoms.to_vec(),
            full: full.to_vec(),
            tolerance_bohr,
            time_reversal,
            map,
            crystal_operations,
            grid,
            dropped_for_grid,
        })
    }
    /// Detected operations dropped for the explicit FFT grid (GS1C).
    pub fn dropped_for_grid(&self) -> usize {
        self.dropped_for_grid
    }
    /// Same as `matches`, for a preparation restricted to an explicit grid.
    pub fn matches_for_grid(
        &self,
        cell: &Cell,
        atoms: &[Atom],
        full: &[KPoint],
        tolerance_bohr: f64,
        time_reversal: bool,
        grid: Option<[usize; 3]>,
    ) -> bool {
        self.grid == grid && self.matches_geometry(cell, atoms, full, tolerance_bohr, time_reversal)
    }
    /// Public configurations may change during relaxation. Reuse only for the
    /// exact geometry, mesh and detection settings that produced this record.
    pub fn matches(
        &self,
        cell: &Cell,
        atoms: &[Atom],
        full: &[KPoint],
        tolerance_bohr: f64,
        time_reversal: bool,
    ) -> bool {
        self.grid.is_none()
            && self.matches_geometry(cell, atoms, full, tolerance_bohr, time_reversal)
    }
    /// Geometry, mesh and detection settings only (any grid restriction).
    pub fn matches_geometry(
        &self,
        cell: &Cell,
        atoms: &[Atom],
        full: &[KPoint],
        tolerance_bohr: f64,
        time_reversal: bool,
    ) -> bool {
        self.lattice == cell.lattice_columns_bohr()
            && self.atoms == atoms
            && self.full == full
            && self.tolerance_bohr.to_bits() == tolerance_bohr.to_bits()
            && self.time_reversal == time_reversal
    }
    pub fn map(&self) -> &dft_geometry::symmetry::KOrbitMap {
        &self.map
    }
    pub fn density_operations(
        &self,
        crystal: bool,
    ) -> Vec<dft_geometry::symmetry::SymmetryOperation> {
        if crystal {
            self.crystal_operations.clone()
        } else {
            self.map
                .operations
                .iter()
                .filter(|o| !o.antiunitary)
                .cloned()
                .collect()
        }
    }
}
