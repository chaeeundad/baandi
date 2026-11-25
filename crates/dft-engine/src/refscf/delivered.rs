//! Post-SCF derivative boundary. The SCF stays distributed; current derivative
//! kernels consume replicated complex orbitals, materialized once on delivery.
use crate::{ElectronicState, ForceState, ScfError, ScfOutcome, refstart::StartInput};

fn validate_model(input: &StartInput<'_>, result: &crate::ScfResult) -> Result<(), ScfError> {
    use crate::force_identity::{ModelDescription, OccupationIdentity};
    let dims = result
        .density
        .shape()
        .try_into()
        .map_err(|_| ScfError::InvalidInput("reference delivery: density shape"))?;
    let model = ModelDescription::new(
        input.cell,
        input.atoms,
        input.potentials,
        input.options,
        dims,
        result.k_points.iter().map(|k| k.point).collect(),
    );
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
    if !outcome.converged {
        return Err(ScfError::InvalidInput(
            "reference derivatives require a converged delivery",
        ));
    }
    let r = &mut outcome.result;
    validate_model(&input, r)?;
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
        symmetry: &[],
    })
}

/// Stress shares the force model, occupations and verification contract.
pub use force_state as stress_state;

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
    validate_model(&input, result)?;
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

/// Dispatch an ionic warm start with the active communicator. Fresh reference
/// runs retain the regular Gamma/k dispatch and default switch configuration.
pub fn run(
    input: StartInput<'_>,
    profile: super::profile::Profile,
    points: &[dft_geometry::kpoints::KPoint],
    seed: Option<super::ladder::Seed>,
) -> Result<super::driver::Run, ScfError> {
    if let Some(seed) = seed {
        if points.len() != 1
            || points[0].fractional != [0.; 3]
            || input.options.gamma_real != crate::GammaReal::On
        {
            return Err(ScfError::InvalidInput(
                "reference k ionic warm seed requires a k-driver seed interface",
            ));
        }
        if let Some(comm) = crate::parallel::SharedComm::active() {
            super::driver::run_seeded_with(input, profile, &comm, seed)
        } else {
            super::driver::run_seeded_with(input, profile, &super::metric::Serial, seed)
        }
    } else {
        super::k_driver::run(input, profile, points)
    }
}
