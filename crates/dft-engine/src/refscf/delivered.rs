//! Post-SCF derivative boundary. The SCF stays distributed; current derivative
//! kernels consume replicated complex orbitals, materialized once on delivery.
use crate::{ElectronicState, ForceState, ScfError, ScfOutcome, refstart::StartInput};

fn validate_model(
    input: &StartInput<'_>,
    result: &crate::ScfResult,
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
        input.options.verification == crate::Verification::Off && r.verification.is_none();
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
    result: &crate::ScfResult,
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
    result: &crate::ScfResult,
    density: &[f64],
    comm: &K,
) -> Result<IonicSeed, ScfError> {
    if result.gamma_real_used {
        return gamma_seed(input, result, density, comm)
            .map(|seed| IonicSeed::Gamma(Box::new(seed)));
    }
    validate_model(&input, result, &[])?;
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
                || input.options.gamma_real != crate::GammaReal::On
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
