//! Post-SCF derivative boundary. The SCF stays distributed; current derivative
//! kernels consume replicated complex orbitals, materialized once on delivery.
use crate::{
    ElectronicState, ForceState, model::ScfError, model::ScfOutcome, refstart::StartInput,
};

fn validate_model(
    input: &StartInput<'_>,
    result: &crate::model::ScfResult,
    operations: &[dft_geometry::symmetry::SymmetryOperation],
) -> Result<(), ScfError> {
    use crate::force_identity::{ModelDescription, OccupationIdentity};
    let dims = result
        .density
        .shape()
        .try_into()
        .map_err(|_| ScfError::InvalidInput("reference delivery: density shape"))?;
    let mut model = ModelDescription::new(
        input.cell,
        input.atoms,
        input.potentials,
        input.options,
        dims,
        result.k_points.iter().map(|k| k.point).collect(),
    );
    model.symmetry = operations.to_vec();
    let blocks = if result.gamma_real_used {
        1
    } else {
        result.k_points.len()
    };
    let policy = OccupationIdentity::new(&vec![result.occupations.as_slice(); blocks], None, None);
    if !result.force_identity.matches(&model, &policy)
        || input.options.grid_override.is_some_and(|grid| grid != dims)
    {
        return Err(ScfError::InvalidInput(
            "reference delivery: model or occupations differ from delivered state",
        ));
    }
    Ok(())
}

/// Validate delivery and materialize all k blocks collectively, once. The
/// returned view borrows the original model and delivered (not check-map) state.
/// G4 can use the existing collinear/finite-temperature ElectronicState variants
/// with their per-spin/per-k occupations; no derivative kernel assumes f=2.
pub fn force_state<'a>(
    input: StartInput<'a>,
    outcome: &'a mut ScfOutcome,
) -> Result<ForceState<'a>, ScfError> {
    force_state_with_operations(input, outcome, &[])
}
/// Same model identity and derivative projector as the legacy symmetry path.
pub fn force_state_prepared<'a>(
    input: StartInput<'a>,
    outcome: &'a mut ScfOutcome,
    prepared: &'a super::symmetry::Prepared,
) -> Result<ForceState<'a>, ScfError> {
    prepared.validate(
        &input,
        &outcome
            .result
            .k_points
            .iter()
            .map(|k| k.point)
            .collect::<Vec<_>>(),
    )?;
    let operations = prepared.operations(input.options.crystal_density_group);
    force_state_with_operations(input, outcome, operations)
}
fn force_state_with_operations<'a>(
    input: StartInput<'a>,
    outcome: &'a mut ScfOutcome,
    operations: &'a [dft_geometry::symmetry::SymmetryOperation],
) -> Result<ForceState<'a>, ScfError> {
    if !outcome.converged {
        return Err(ScfError::InvalidInput(
            "reference derivatives require a converged delivery",
        ));
    }
    let r = &mut outcome.result;
    validate_model(&input, r, operations)?;
    let verified = r.verification.as_ref().is_some_and(|v| v.passed());
    let opted_out =
        input.options.verification == crate::model::Verification::Off && r.verification.is_none();
    if r.standard.is_none() || (!verified && !opted_out) {
        return Err(ScfError::InvalidInput(
            "reference derivatives require verification or explicit verify=off",
        ));
    }
    r.materialize_orbitals()?;
    r.force_identity.unverified_delivery = opted_out && !verified;
    r.force_identity.reference_delivery = true;
    Ok(ForceState {
        cell: input.cell,
        atoms: input.atoms,
        potentials: input.potentials,
        options: input.options,
        electronic: ElectronicState::Fixed(r),
        density_floor: None,
        symmetry: operations,
    })
}

/// Stress shares the force model, occupations and verification contract.
/// The retained options include the frozen reference cell: derivative kernels
/// rebuild that ordered Miller support at input.cell, not a new cutoff sphere.
/// The separate Pulay assessment keeps its nominal-cutoff scan protocol.
pub use force_state as stress_state;
pub use force_state_prepared as stress_state_prepared;

/// Export a real-Gamma ionic warm start in the existing owned ladder format.
/// Replicated delivery arrays are read locally; no orbital collective is added.
/// The next start rebuilds ionic fields and reorthogonalizes these orbitals.
pub fn gamma_seed<K: dft_mpi::Collective>(
    input: StartInput<'_>,
    result: &crate::model::ScfResult,
    density: &[f64],
    comm: &K,
) -> Result<super::ladder::Seed, ScfError> {
    use dft_basis::{FftGridSpec, HalfSphere};
    use dft_numerics::{Fft3d, RealMatrix};
    use dft_physics::data::{Array3, Complex64, FftPlan3d};
    use dft_solvers::distributed_fft::position_owners;
    validate_model(&input, result, &[])?;
    if !result.orbitals_delivered()
        || result.k_points.len() != 1
        || result.k_points[0].point.fractional != [0.; 3]
        || !result.gamma_real_used
    {
        return Err(ScfError::InvalidInput(
            "reference ionic seed requires materialized real Gamma orbitals",
        ));
    }
    let dims: [usize; 3] = result
        .density
        .shape()
        .try_into()
        .map_err(|_| ScfError::InvalidInput("reference ionic seed: density shape"))?;
    if density.len() != result.density.len() || density.iter().any(|x| !x.is_finite()) {
        return Err(ScfError::InvalidInput(
            "reference ionic seed: invalid density",
        ));
    }
    let charge = density.iter().sum::<f64>() * input.cell.volume_bohr3() / density.len() as f64;
    if !charge.is_finite() || (charge - result.charge).abs() > 1e-8 * result.charge.abs().max(1.) {
        return Err(ScfError::InvalidInput(
            "reference ionic seed: density charge mismatch",
        ));
    }
    let basis = input.options.basis(input.cell, [0.; 3], 1.)?;
    let db = input.options.basis(input.cell, [0.; 3], 4.)?;
    let grid = FftGridSpec::new(&basis, &db, Some(dims))?;
    let half = HalfSphere::new(&basis, &grid).map_err(super::numerical)?;
    let owners = position_owners(dims, grid.g_to_fft(), comm.size()).map_err(super::numerical)?;
    let rows: Vec<_> = (0..basis.vectors().len())
        .filter(|&r| owners[half.plus()[r.div_ceil(2)] as usize] == comm.rank())
        .collect();
    let full = &result.k_points[0].solution.eigenvectors;
    if full.rows() != basis.vectors().len() || full.cols() != result.occupations.len() {
        return Err(ScfError::InvalidInput(
            "reference ionic seed: orbital shape",
        ));
    }
    let mut orbitals = RealMatrix::zeros(rows.len(), full.cols()).map_err(super::numerical)?;
    let mut packed = vec![0.; full.rows()];
    for band in 0..full.cols() {
        half.pack(
            &full.as_slice()[band * full.rows()..(band + 1) * full.rows()],
            &mut packed,
        )
        .map_err(super::numerical)?;
        for (i, &r) in rows.iter().enumerate() {
            orbitals.as_mut_slice()[band * rows.len() + i] = packed[r];
        }
    }
    let positions: Vec<_> = (0..dims[0])
        .flat_map(|x| (0..dims[1]).flat_map(move |y| (0..dims[2]).map(move |z| [x, y, z])))
        .collect();
    let owners = position_owners(dims, &positions, comm.size()).map_err(super::numerical)?;
    let indices: Vec<_> = (0..density.len())
        .filter(|&i| owners[i] == comm.rank())
        .collect();
    let real = Array3::from_shape_vec(
        dims,
        density.iter().map(|&v| Complex64::new(v, 0.)).collect(),
    )
    .map_err(super::numerical)?;
    let mut reciprocal = Array3::zeros(dims);
    FftPlan3d::new(dims)?.to_reciprocal(&real, &mut reciprocal)?;
    let flat = reciprocal
        .as_slice()
        .ok_or(ScfError::InvalidInput("ionic seed: noncontiguous density"))?;
    let density_g = indices.iter().map(|&i| flat[i]).collect();
    Ok(super::ladder::Seed {
        basis,
        grid,
        half,
        orbitals,
        rows,
        density_g,
        indices,
    })
}

/// Fixed-cell ionic continuation. K seeds retain the delivered full complex
/// blocks; injection selects the current communicator's owned rows.
#[derive(Clone)]
pub enum IonicSeed {
    Gamma(Box<super::ladder::Seed>),
    K(KSeed),
}
#[derive(Clone)]
pub struct KSeed {
    lattice: [[f64; 3]; 3],
    strain_reference: Option<Box<[[f64; 3]; 3]>>,
    cutoff: f64,
    dims: [usize; 3],
    points: Vec<dft_geometry::kpoints::KPoint>,
    orbitals: Vec<dft_numerics::linalg::ComplexMatrix>,
    density: Vec<f64>,
    charge: f64,
}

pub fn ionic_seed<K: dft_mpi::Collective>(
    input: StartInput<'_>,
    result: &crate::model::ScfResult,
    density: &[f64],
    comm: &K,
) -> Result<IonicSeed, ScfError> {
    if result.gamma_real_used {
        return gamma_seed(input, result, density, comm)
            .map(|seed| IonicSeed::Gamma(Box::new(seed)));
    }
    k_seed(input, result, density, &[])
}

/// Export a seed with the delivered symmetry identity intact. Prepared solves
/// use complex k blocks, including Gamma; no representation conversion occurs.
pub fn ionic_seed_prepared(
    input: StartInput<'_>,
    result: &crate::model::ScfResult,
    density: &[f64],
    prepared: &super::symmetry::Prepared,
) -> Result<IonicSeed, ScfError> {
    prepared.validate(
        &input,
        &result.k_points.iter().map(|k| k.point).collect::<Vec<_>>(),
    )?;
    if result.gamma_real_used {
        return Err(ScfError::InvalidInput(
            "prepared ionic seed requires complex k blocks",
        ));
    }
    let operations = prepared.operations(input.options.crystal_density_group);
    k_seed(input, result, density, operations)
}

fn k_seed(
    input: StartInput<'_>,
    result: &crate::model::ScfResult,
    density: &[f64],
    operations: &[dft_geometry::symmetry::SymmetryOperation],
) -> Result<IonicSeed, ScfError> {
    validate_model(&input, result, operations)?;
    let charge = density.iter().sum::<f64>() * input.cell.volume_bohr3() / density.len() as f64;
    if !result.orbitals_delivered()
        || density.len() != result.density.len()
        || density.iter().any(|x| !x.is_finite())
        || !charge.is_finite()
        || (charge - result.charge).abs() > 1e-8 * result.charge.abs().max(1.)
    {
        return Err(ScfError::InvalidInput(
            "reference k ionic seed: invalid delivery or density",
        ));
    }
    Ok(IonicSeed::K(KSeed {
        lattice: input.cell.lattice_columns_bohr(),
        strain_reference: input
            .options
            .strain_reference
            .as_ref()
            .map(|cell| Box::new(cell.lattice_columns_bohr())),
        cutoff: input.options.cutoff_hartree,
        dims: result
            .density
            .shape()
            .try_into()
            .map_err(|_| ScfError::InvalidInput("ionic seed: grid"))?,
        points: result.k_points.iter().map(|k| k.point).collect(),
        orbitals: result
            .k_points
            .iter()
            .map(|k| k.solution.eigenvectors.clone())
            .collect(),
        density: density.to_vec(),
        charge: result.charge,
    }))
}
impl IonicSeed {
    /// Affine fractional-grid transfer within one frozen-Miller segment.
    /// Coefficients are normalized in the cell-volume-normalized plane-wave
    /// basis, so only density scales by Omega_old/Omega_new. A rebuild starts
    /// fresh instead of transferring across different reference supports.
    pub fn in_cell(
        mut self,
        source: &dft_geometry::Cell,
        target: &dft_geometry::Cell,
        options: &crate::model::ScfOptions,
    ) -> Result<Self, ScfError> {
        if options.strain_reference.is_none() || options.grid_override.is_none() {
            return Err(ScfError::InvalidInput(
                "cell seed requires frozen basis and grid",
            ));
        }
        let ratio = source.volume_bohr3() / target.volume_bohr3();
        match &mut self {
            Self::Gamma(seed) => {
                let basis = options.basis(target, [0.; 3], 1.)?;
                if seed.basis.cutoff().0 != options.cutoff_hartree
                    || seed.grid.dims() != options.grid_override.unwrap()
                    || seed
                        .basis
                        .vectors()
                        .iter()
                        .map(|g| g.miller)
                        .ne(basis.vectors().iter().map(|g| g.miller))
                {
                    return Err(ScfError::InvalidInput(
                        "cell seed: frozen Miller map changed",
                    ));
                }
                seed.basis = basis;
                for value in &mut seed.density_g {
                    *value *= ratio;
                }
            }
            Self::K(seed) => {
                if seed.lattice != source.lattice_columns_bohr()
                    || seed.strain_reference.as_deref().copied()
                        != options
                            .strain_reference
                            .as_ref()
                            .map(|c| c.lattice_columns_bohr())
                    || seed.cutoff != options.cutoff_hartree
                    || Some(seed.dims) != options.grid_override
                {
                    return Err(ScfError::InvalidInput("cell seed: source segment changed"));
                }
                seed.lattice = target.lattice_columns_bohr();
                for value in &mut seed.density {
                    *value *= ratio;
                }
            }
        }
        Ok(self)
    }

    /// Replace only the density of an ionic seed (C3 predictor). Orbitals and
    /// their Miller map are retained; the input is on the source cell's grid.
    pub fn with_density(
        mut self,
        density: &[f64],
        volume: f64,
        charge: f64,
    ) -> Result<Self, ScfError> {
        use dft_numerics::Fft3d;
        use dft_physics::data::{Array3, Complex64, FftPlan3d};
        let dims = match &self {
            Self::Gamma(s) => s.grid.dims(),
            Self::K(s) => s.dims,
        };
        let integral = density.iter().sum::<f64>() * volume / density.len() as f64;
        if !volume.is_finite()
            || volume <= 0.
            || !charge.is_finite()
            || charge <= 0.
            || density.len() != dims.iter().product::<usize>()
            || density.iter().any(|v| !v.is_finite() || *v < 0.)
            || !integral.is_finite()
            || (integral - charge).abs() > 1e-8 * charge.abs().max(1.)
        {
            return Err(ScfError::InvalidInput(
                "ionic predictor: invalid density or charge",
            ));
        }
        match &mut self {
            Self::Gamma(seed) => {
                let real = Array3::from_shape_vec(
                    dims,
                    density.iter().map(|&v| Complex64::new(v, 0.)).collect(),
                )
                .map_err(super::numerical)?;
                let mut reciprocal = Array3::zeros(dims);
                FftPlan3d::new(dims)?.to_reciprocal(&real, &mut reciprocal)?;
                let flat = reciprocal
                    .as_slice()
                    .ok_or(ScfError::InvalidInput("ionic predictor: density layout"))?;
                seed.density_g = seed.indices.iter().map(|&i| flat[i]).collect();
            }
            Self::K(seed) => seed.density = density.to_vec(),
        }
        Ok(self)
    }
}

impl KSeed {
    pub(super) fn inject(
        self,
        input: &StartInput<'_>,
        state: &mut crate::refstart::KStartState,
    ) -> Result<(), ScfError> {
        if self.lattice != input.cell.lattice_columns_bohr()
            || self.strain_reference.as_deref().copied()
                != input
                    .options
                    .strain_reference
                    .as_ref()
                    .map(|cell| cell.lattice_columns_bohr())
            || self.cutoff != input.options.cutoff_hartree
            || self.dims != state.grid.dims()
            || self.points != state.blocks.iter().map(|b| b.point).collect::<Vec<_>>()
            || (self.charge - state.electrons).abs() > 1e-8 * state.electrons.max(1.)
        {
            return Err(ScfError::InvalidInput(
                "reference k ionic seed: incompatible cell, basis, grid, mesh or charge",
            ));
        }
        for ((full, block), owned) in self
            .orbitals
            .iter()
            .zip(&state.blocks)
            .zip(&mut state.orbitals.0)
        {
            if full.rows() != block.basis.vectors().len()
                || full.cols() != owned.cols()
                || full
                    .as_slice()
                    .iter()
                    .any(|z| !z.re.is_finite() || !z.im.is_finite())
            {
                return Err(ScfError::InvalidInput(
                    "reference k ionic seed: orbital shape or values",
                ));
            }
            for band in 0..owned.cols() {
                for (i, &g) in block.operator.local_rows().iter().enumerate() {
                    owned[(i, band)] = full[(g, band)];
                }
            }
        }
        state.density = state
            .fields
            .plane_offsets()
            .iter()
            .map(|&i| self.density[i])
            .collect();
        Ok(())
    }
}

/// Dispatch an ionic warm start with the active communicator. Fresh reference
/// runs retain the regular Gamma/k dispatch and default switch configuration.
pub fn run(
    input: StartInput<'_>,
    profile: super::profile::Profile,
    points: &[dft_geometry::kpoints::KPoint],
    seed: Option<IonicSeed>,
) -> Result<super::driver::Run, ScfError> {
    match seed {
        Some(IonicSeed::Gamma(seed)) => {
            if points.len() != 1
                || points[0].fractional != [0.; 3]
                || input.options.gamma_real != crate::model::GammaReal::On
            {
                return Err(ScfError::InvalidInput(
                    "reference ionic seed: Gamma/k mismatch",
                ));
            }
            if let Some(comm) = crate::parallel::SharedComm::active() {
                super::driver::run_seeded_with(input, profile, &comm, *seed)
            } else {
                super::driver::run_seeded_with(input, profile, &super::metric::Serial, *seed)
            }
        }
        Some(IonicSeed::K(seed)) => {
            if let Some(comm) = crate::parallel::SharedComm::active() {
                super::k_driver::run_seeded_with(input, profile, points, &comm, seed)
            } else {
                super::k_driver::run_seeded_with(
                    input,
                    profile,
                    points,
                    &super::metric::Serial,
                    seed,
                )
            }
        }
        None => super::k_driver::run(input, profile, points),
    }
}

#[cfg(test)]
mod cell_seed_tests {
    use super::*;
    use dft_basis::{FftGridSpec, HalfSphere};
    use dft_geometry::Cell;
    use dft_numerics::RealMatrix;
    use dft_physics::data::Complex64;

    #[test]
    fn cell_seed_affine_charge_and_frozen_orbitals() {
        let source =
            Cell::from_user_rows_bohr([[5., 0., 0.], [0.4, 4., 0.], [0., 0.3, 6.]]).unwrap();
        let target =
            Cell::from_user_rows_bohr([[5.2, 0., 0.], [0.5, 3.9, 0.], [0., 0.2, 6.1]]).unwrap();
        let mut options = crate::model::ScfOptions {
            cutoff_hartree: 1.,
            strain_reference: Some(Box::new(source.clone())),
            ..Default::default()
        };
        let basis = options.basis(&source, [0.; 3], 1.).unwrap();
        let db = options.basis(&source, [0.; 3], 4.).unwrap();
        let grid = FftGridSpec::new(&basis, &db, None).unwrap();
        options.grid_override = Some(grid.dims());
        let half = HalfSphere::new(&basis, &grid).unwrap();
        let rows: Vec<_> = (0..basis.vectors().len()).collect();
        let orbitals = RealMatrix::zeros(rows.len(), 1).unwrap();
        let dims = grid.dims();
        let count = dims.iter().product();
        let seed = IonicSeed::Gamma(Box::new(crate::refscf::ladder::Seed {
            basis,
            grid,
            half,
            orbitals,
            rows,
            density_g: vec![Complex64::default(); count],
            indices: (0..count).collect(),
        }));
        let charge = 2.;
        let density = vec![charge / source.volume_bohr3(); count];
        let seed = seed
            .with_density(&density, source.volume_bohr3(), charge)
            .unwrap();
        let moved = seed.clone().in_cell(&source, &target, &options).unwrap();
        let (IonicSeed::Gamma(old), IonicSeed::Gamma(new)) = (&seed, &moved) else {
            panic!()
        };
        assert_eq!(old.orbitals.as_slice(), new.orbitals.as_slice());
        assert_eq!(old.rows, new.rows);
        assert_eq!(old.half.plus(), new.half.plus());
        assert!((new.density_g[0].re * target.volume_bohr3() - charge).abs() < 1e-13);
        assert!(new.density_g[1..].iter().all(|z| z.norm() < 1e-14));
        assert!(
            old.basis.vectors().iter().map(|g| g.miller).eq(new
                .basis
                .vectors()
                .iter()
                .map(|g| g.miller))
        );
        assert_ne!(
            old.basis.vectors()[1].norm_squared,
            new.basis.vectors()[1].norm_squared
        );
        let k = IonicSeed::K(KSeed {
            lattice: source.lattice_columns_bohr(),
            strain_reference: Some(Box::new(source.lattice_columns_bohr())),
            cutoff: options.cutoff_hartree,
            dims,
            points: vec![],
            orbitals: vec![],
            density: density.clone(),
            charge,
        });
        let IonicSeed::K(k_moved) = k.clone().in_cell(&source, &target, &options).unwrap() else {
            panic!()
        };
        assert_eq!(k_moved.lattice, target.lattice_columns_bohr());
        assert!(
            (k_moved.density.iter().sum::<f64>() * target.volume_bohr3() / count as f64 - charge)
                .abs()
                < 1e-12
        );
        assert!(k.clone().in_cell(&target, &source, &options).is_err());
        let mut changed_reference = options.clone();
        changed_reference.strain_reference = Some(Box::new(target.clone()));
        assert!(k.in_cell(&source, &target, &changed_reference).is_err());
        let mut invalid = options.clone();
        invalid.strain_reference = None;
        assert!(seed.clone().in_cell(&source, &target, &invalid).is_err());
        invalid = options.clone();
        invalid.cutoff_hartree *= 2.;
        assert!(seed.clone().in_cell(&source, &target, &invalid).is_err());
        assert!(
            seed.clone()
                .with_density(&density, source.volume_bohr3(), 3.)
                .is_err()
        );
        let mut invalid_density = density;
        invalid_density[0] = f64::NAN;
        assert!(
            seed.with_density(&invalid_density, source.volume_bohr3(), charge)
                .is_err()
        );
    }
}
