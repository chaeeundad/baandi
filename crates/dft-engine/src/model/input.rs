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
}
impl PreparedSymmetry {
    pub fn new(
        cell: &Cell,
        atoms: &[Atom],
        full: &[KPoint],
        tolerance_bohr: f64,
        time_reversal: bool,
    ) -> Result<Self, ScfError> {
        let crystal_operations = dft_geometry::symmetry::detect(cell, atoms, tolerance_bohr)
            .map_err(|e| ScfError::Numerical(e.to_string()))?;
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
        })
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
