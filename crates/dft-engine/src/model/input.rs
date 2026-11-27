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
    /// Requested k reduction policy (GT2, `K_REDUCTION`); `map().reduction`
    /// records the one applied.
    requested: dft_geometry::symmetry::KReduction,
    /// GS1D-B (`GRID_INCOMPATIBLE_OPS=reciprocal`): operations whose rotation
    /// does not map `grid` were kept.
    rotations_free: bool,
}
/// GS1D-B: `GRID_INCOMPATIBLE_OPS=reciprocal` in the calling thread's
/// configuration (ADR 0049).
pub fn grid_rotations_free() -> bool {
    dft_numerics::switches::choice(dft_numerics::switches::Id::GRID_INCOMPATIBLE_OPS)
        == "reciprocal"
}
/// Whether a symmetrized run on `dims` must project its local potential onto
/// `operations` in G space (GS1D-B): some operation does not permute the
/// real-space grid, so the XC quadrature is not invariant.
pub fn needs_potential_projection(
    dims: [usize; 3],
    operations: &[dft_geometry::symmetry::SymmetryOperation],
) -> bool {
    !operations.is_empty()
        && grid_rotations_free()
        && !dft_geometry::symmetry::preserves_fft_grid(dims, operations)
}
/// The grid check of a symmetrized run: every operation permutes `dims`, or
/// (GS1D-B) every translation is commensurate and the potential is projected.
pub fn grid_supports(
    dims: [usize; 3],
    operations: &[dft_geometry::symmetry::SymmetryOperation],
) -> bool {
    let free = grid_rotations_free();
    !dims.contains(&0)
        && operations
            .iter()
            .all(|o| dft_geometry::symmetry::grid_keeps(dims, o, free))
}
/// The `K_REDUCTION` switch of the calling thread's configuration (GT2,
/// ADR 0051).
pub fn requested_k_reduction() -> dft_geometry::symmetry::KReduction {
    use dft_geometry::symmetry::KReduction;
    match dft_numerics::switches::choice(dft_numerics::switches::Id::K_REDUCTION) {
        "mesh" => KReduction::MeshSubgroup,
        "unreduced" => KReduction::Unreduced,
        _ => KReduction::Pointwise,
    }
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
        let rotations_free = grid_rotations_free();
        let crystal_operations = match grid {
            None => detected.clone(),
            Some(dims) => dft_geometry::symmetry::grid_subgroup_with(
                cell,
                &detected,
                dims,
                tolerance_bohr,
                rotations_free,
            ),
        };
        let dropped_for_grid = detected.len() - crystal_operations.len();
        let requested = requested_k_reduction();
        let map = dft_geometry::symmetry::reduce_mesh_with(
            full,
            &crystal_operations,
            time_reversal,
            requested,
        )
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
            requested,
            rotations_free,
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
        self.grid == grid
            && self.requested == requested_k_reduction()
            && self.rotations_free == grid_rotations_free()
            && self.matches_geometry(cell, atoms, full, tolerance_bohr, time_reversal)
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
            && self.requested == requested_k_reduction()
            && self.rotations_free == grid_rotations_free()
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
