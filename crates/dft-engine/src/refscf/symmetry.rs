//! Geometry-checked symmetry input and owned-stick density projection.
use crate::{
    PreparedSymmetry, ScfError,
    refstart::{StartFields, StartInput},
};
use dft_geometry::{Atom, Cell, kpoints::KPoint, symmetry::SymmetryOperation};
use dft_mpi::Collective;
use dft_physics::{data::Complex64, symmetry::DistributedDensitySymmetrizer};

/// A checked view of the same preparation used by legacy SCF. The geometry
/// identity is retained so a caller cannot reuse a map after moving atoms.
pub struct Prepared {
    pub(super) time_reversal: bool,
    lattice: [[f64; 3]; 3],
    atoms: Vec<Atom>,
    points: Vec<KPoint>,
    mesh_operations: Vec<SymmetryOperation>,
    crystal_operations: Vec<SymmetryOperation>,
}
impl Prepared {
    pub fn new(
        cell: &Cell,
        atoms: &[Atom],
        full: &[KPoint],
        tolerance_bohr: f64,
        time_reversal: bool,
    ) -> Result<Self, ScfError> {
        let prepared = PreparedSymmetry::new(cell, atoms, full, tolerance_bohr, time_reversal)?;
        Self::from_prepared(cell, atoms, full, tolerance_bohr, time_reversal, &prepared)
    }
    pub fn from_prepared(
        cell: &Cell,
        atoms: &[Atom],
        full: &[KPoint],
        tolerance_bohr: f64,
        time_reversal: bool,
        prepared: &PreparedSymmetry,
    ) -> Result<Self, ScfError> {
        if !prepared.matches(cell, atoms, full, tolerance_bohr, time_reversal) {
            return Err(ScfError::InvalidInput(
                "reference symmetry: stale prepared geometry or mesh",
            ));
        }
        Ok(Self {
            time_reversal,
            lattice: cell.lattice_columns_bohr(),
            atoms: atoms.to_vec(),
            points: prepared.map().points.clone(),
            mesh_operations: prepared.density_operations(false),
            crystal_operations: prepared.density_operations(true),
        })
    }
    pub fn points(&self) -> &[KPoint] {
        &self.points
    }
    pub fn operations(&self, crystal: bool) -> &[SymmetryOperation] {
        if crystal {
            &self.crystal_operations
        } else {
            &self.mesh_operations
        }
    }
    pub(super) fn validate(
        &self,
        input: &StartInput<'_>,
        points: &[KPoint],
    ) -> Result<(), ScfError> {
        if self.lattice != input.cell.lattice_columns_bohr()
            || self.atoms != input.atoms
            || self.points != points
        {
            return Err(ScfError::InvalidInput(
                "reference symmetry: geometry or irreducible k list differs",
            ));
        }
        Ok(())
    }
    pub(super) fn grid(&self, input: &StartInput<'_>) -> Result<[usize; 3], ScfError> {
        let db = input.options.basis(input.cell, [0.; 3], 4.)?;
        let mut dims = [0; 3];
        for p in &self.points {
            let b = input.options.basis(input.cell, p.fractional, 1.)?;
            let d = dft_basis::FftGridSpec::new(&b, &db, input.options.grid_override)?.dims();
            for i in 0..3 {
                dims[i] = dims[i].max(d[i]);
            }
        }
        let ops = self.operations(input.options.crystal_density_group);
        if input.options.grid_override.is_some() {
            if !dft_geometry::symmetry::preserves_fft_grid(dims, ops) {
                return Err(ScfError::InvalidInput(
                    "reference symmetry: incompatible explicit FFT grid",
                ));
            }
            Ok(dims)
        } else {
            dft_geometry::symmetry::compatible_fft_grid(dims, ops)
                .map_err(|e| super::numerical(e).into())
        }
    }
}
#[derive(Debug, Clone)]
pub struct Residual {
    pub iteration: usize,
    pub verification: bool,
    /// sqrt(integral sum_spin |rho - P rho|^2), electrons/Bohr^(3/2).
    pub norm: f64,
    /// ||rho - P rho|| / ||rho|| before projection, dimensionless. Zero for
    /// zero density; an irreducible-k sum need not be invariant at convergence.
    pub pre_projection_relative_norm: f64,
    /// ||P(P rho) - P rho||, electrons/Bohr^(3/2). Diagnostic only: the
    /// second projection never replaces the physical density.
    pub post_projection_invariance_norm: f64,
}
#[derive(Debug, Clone)]
pub struct Report {
    pub group_order: usize,
    pub irreducible_k_count: usize,
    pub channels: usize,
    pub maps: Vec<Residual>,
}
pub(super) struct Workspace {
    projector: DistributedDensitySymmetrizer,
    pub report: Report,
    pub iteration: usize,
    pub verification: bool,
    volume: f64,
    invariance_scratch: Vec<Vec<Complex64>>,
}
impl Workspace {
    pub fn new<K: Collective>(
        input: &StartInput<'_>,
        prepared: &Prepared,
        dims: [usize; 3],
        fields: &StartFields,
        channels: usize,
        comm: &K,
    ) -> Result<Self, ScfError> {
        prepared.validate(input, prepared.points())?;
        let ops = prepared.operations(input.options.crystal_density_group);
        let db = input.options.basis(input.cell, [0.; 3], 4.)?;
        let projector =
            DistributedDensitySymmetrizer::new(&db, dims, ops, comm.rank(), comm.size(), channels)?;
        if projector.local_indices() != fields.reciprocal_indices() {
            return Err(ScfError::InvalidInput(
                "reference symmetry: stick layout mismatch",
            ));
        }
        let invariance_scratch =
            vec![vec![Complex64::default(); projector.local_indices().len()]; channels];
        Ok(Self {
            projector,
            invariance_scratch,
            report: Report {
                group_order: ops.len(),
                irreducible_k_count: prepared.points.len(),
                channels,
                maps: Vec::with_capacity(input.options.max_iterations * 2),
            },
            iteration: 0,
            verification: false,
            volume: input.cell.volume_bohr3(),
        })
    }
    pub fn apply<K: Collective>(
        &mut self,
        fields: &mut StartFields,
        comm: &K,
        density: &mut [Vec<f64>],
    ) -> Result<(), ScfError> {
        let mut coefficients = density
            .iter()
            .map(|n| fields.forward(comm, n))
            .collect::<Result<Vec<_>, _>>()?;
        let density_norm2: f64 = coefficients.iter().flatten().map(|z| z.norm_sqr()).sum();
        let local = self.projector.apply(comm, &mut coefficients)?;
        for (scratch, projected) in self.invariance_scratch.iter_mut().zip(&coefficients) {
            scratch.copy_from_slice(projected);
        }
        let post = self.projector.apply(comm, &mut self.invariance_scratch)?;
        let mut total = [0.; 3];
        comm.tree_sum_f64(&[local, density_norm2, post], &mut total)
            .map_err(super::numerical)?;
        for (n, g) in density.iter_mut().zip(&coefficients) {
            *n = fields.inverse(comm, g)?;
        }
        self.report.maps.push(Residual {
            iteration: self.iteration,
            verification: self.verification,
            norm: (self.volume * total[0]).sqrt(),
            pre_projection_relative_norm: if total[1] == 0. {
                0.
            } else {
                (total[0] / total[1]).sqrt()
            },
            post_projection_invariance_norm: (self.volume * total[2]).sqrt(),
        });
        Ok(())
    }
    /// Legacy checks the initial density without replacing it. Retain the same
    /// 1e-12 pointwise allowance and leave the original values untouched.
    pub fn check_start<K: Collective>(
        &mut self,
        fields: &mut StartFields,
        comm: &K,
        density: &[Vec<f64>],
    ) -> Result<(), ScfError> {
        let mut projected = density.to_vec();
        self.apply(fields, comm, &mut projected)?;
        self.report.maps.clear();
        let bad = density
            .iter()
            .flatten()
            .zip(projected.iter().flatten())
            .any(|(a, b)| !a.is_finite() || (a - b).abs() > 1e-12);
        let mut agreed = [0.];
        comm.tree_sum_f64(&[f64::from(bad)], &mut agreed)
            .map_err(super::numerical)?;
        if agreed[0] != 0. {
            return Err(ScfError::InvalidInput(
                "symmetry.density: input breaks the selected symmetry; use full mesh",
            ));
        }
        Ok(())
    }
}
