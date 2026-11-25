//! Residual-certified spectra of a fixed supplied density. No SCF map, density
//! update, mixing, occupations solve or self-consistency certificate is run.
use crate::{
    Functional,
    model::ScfError,
    refstart::{self, KStartState, StartInput},
};
use dft_geometry::{Cell, kpoints::KPoint};
use dft_mpi::Collective;
use dft_numerics::linalg::ComplexMatrix;
use dft_physics::data::Array3;
use dft_solvers::{ref_davidson, ref_davidson_complex};

#[derive(Debug)]
pub struct Bands {
    /// Hartree, ascending band order.
    pub eigenvalues: Vec<f64>,
    /// Fresh ||H psi - epsilon psi||, Hartree, for every requested band.
    pub residuals: Vec<f64>,
    /// This rank's owned rows, column-major, normalized over all ranks.
    pub eigenvectors: ComplexMatrix,
    /// Global plane-wave indices corresponding to the local eigenvector rows.
    pub rows: Vec<usize>,
}
#[derive(Debug)]
pub struct NscfResult {
    /// Original coordinates/weights, without wrapping path coordinates.
    pub points: Vec<KPoint>,
    /// Channel-major then k-major. One total-density or two up/down channels.
    pub channels: Vec<Vec<Bands>>,
    pub electrons: f64,
    pub tolerance_hartree: f64,
}

/// `density` contains one total-density array or two collinear arrays (up/down)
/// in electrons/Bohr^3. The supplied grid is retained without interpolation,
/// clipping or charge normalization. `bands` is the total count per channel.
/// Every rank processes every k point sequentially, including complex Gamma.
/// All arguments must be identical on every rank of the active communicator.
pub fn nscf(
    input: StartInput<'_>,
    density: &[Array3<f64>],
    points: &[KPoint],
    bands: usize,
    eigen_tolerance: f64,
) -> Result<NscfResult, ScfError> {
    if let Some(comm) = crate::parallel::SharedComm::active() {
        nscf_with(input, density, points, bands, eigen_tolerance, &comm)
    } else {
        nscf_with(
            input,
            density,
            points,
            bands,
            eigen_tolerance,
            &super::metric::Serial,
        )
    }
}

/// Explicit communicator variant, also used by deterministic collective tests.
pub fn nscf_with<K: Collective>(
    input: StartInput<'_>,
    density: &[Array3<f64>],
    points: &[KPoint],
    bands: usize,
    eigen_tolerance: f64,
    comm: &K,
) -> Result<NscfResult, ScfError> {
    if !eigen_tolerance.is_finite() || eigen_tolerance <= 0. || bands == 0 {
        return Err(ScfError::InvalidInput(
            "reference NSCF: invalid bands/tolerance",
        ));
    }
    let mut options = input.options.clone();
    options.bands = Some(bands);
    let (mut state, potentials) = refstart::start_kpoints_from_density(
        StartInput {
            options: &options,
            ..input
        },
        density,
        points,
        comm,
    )?;
    let dims = state.grid.dims();
    let mut channels = Vec::with_capacity(potentials.len());
    for (channel, potential) in potentials.iter().enumerate() {
        let mut average = [0.];
        comm.tree_sum_f64(&[potential.iter().sum::<f64>()], &mut average)
            .map_err(super::numerical)?;
        let v0 = average[0] / dims.iter().product::<usize>() as f64;
        let mut spectra = Vec::with_capacity(points.len());
        for k in 0..points.len() {
            let op = &mut state.blocks[k].operator;
            let mut planes = crate::grid_dist::PlaneMap::new(
                dims,
                comm.rank(),
                comm.size(),
                state.common.fields.first_axis(),
                op.first_axis(),
            );
            let mut local = Vec::new();
            planes.a_to_b(comm, potential, &mut local)?;
            op.set_owned_potential(&local, channel as u64 + 2)
                .map_err(super::numerical)?;
            let diagonal: Vec<_> = op.diagonal_static().iter().map(|v| v + v0).collect();
            let result = solve(
                comm,
                op,
                &state.orbitals.0[k],
                &diagonal,
                bands,
                eigen_tolerance,
            )
            .map_err(|e| {
                ScfError::Numerical(format!(
                    "reference NSCF channel {channel}, k index {k}: {e}"
                ))
            })?;
            spectra.push(result);
        }
        channels.push(spectra);
    }
    Ok(NscfResult {
        points: points.to_vec(),
        channels,
        electrons: state.electrons,
        tolerance_hartree: eigen_tolerance,
    })
}

pub(crate) fn potentials<K: Collective>(
    state: &mut KStartState,
    comm: &K,
    cell: &Cell,
    functional: Functional,
    density: &[Vec<f64>],
) -> Result<Vec<Vec<f64>>, ScfError> {
    let mut workspace =
        super::electronic_map::Workspace::new(state, comm, cell, density.len() == 2, functional)?;
    let fields =
        super::electronic_map::fields(&mut workspace, state, comm, cell, functional, density)?;
    Ok(fields
        .hxc
        .iter()
        .map(|channel| {
            channel
                .iter()
                .zip(&state.ionic)
                .map(|(x, v)| x + v)
                .collect()
        })
        .collect())
}

fn solve<K: Collective>(
    comm: &K,
    op: &mut dft_physics::distributed_complex::DistributedComplexHamiltonian,
    start: &ComplexMatrix,
    diagonal: &[f64],
    bands: usize,
    tolerance: f64,
) -> Result<Bands, ScfError> {
    // Ritz changes are only an inner stopping rule. All bands are marked
    // occupied to disable the SCF empty-band relaxation. Restarting the solve
    // reactivates all bands even when a previous Ritz-change test passed.
    // A minimum positive normal threshold makes Ritz stagnation (zero change),
    // rather than an SCF energy tolerance, the inner stopping condition. Only
    // the fresh residual below controls NSCF acceptance.
    let mut options = ref_davidson::Options::new(bands, f64::MIN_POSITIVE, bands);
    options.return_hx = false;
    let mut vectors = start.clone();
    let mut hx = ComplexMatrix::zeros(start.rows(), bands).map_err(super::numerical)?;
    let mut local = vec![0.; bands];
    let mut residuals = vec![0.; bands];
    for _ in 0..64 {
        let solution = ref_davidson_complex::solve(comm, op, &vectors, diagonal, options)
            .map_err(super::numerical)?;
        vectors = solution.eigenvectors;
        // A fresh application certifies the returned vectors, independent of
        // cached/rotated HX and the solver's Ritz-change convergence flags.
        op.apply_block(comm, &vectors, &mut hx)
            .map_err(super::numerical)?;
        for (j, value) in local.iter_mut().enumerate() {
            *value = (0..vectors.rows())
                .map(|i| (hx[(i, j)] - solution.eigenvalues[j] * vectors[(i, j)]).norm_sqr())
                .sum();
        }
        comm.tree_sum_f64(&local, &mut residuals)
            .map_err(super::numerical)?;
        for r in &mut residuals {
            *r = r.sqrt();
        }
        if solution.eigenvalues.iter().all(|e| e.is_finite())
            && residuals.iter().all(|r| r.is_finite() && *r <= tolerance)
        {
            return Ok(Bands {
                eigenvalues: solution.eigenvalues,
                residuals,
                eigenvectors: vectors,
                rows: op.local_rows().to_vec(),
            });
        }
    }
    Err(ScfError::Numerical(format!(
        "unconverged band residual after 64 calls: max {:e} Ha, requested {tolerance:e} Ha",
        residuals.iter().copied().fold(0., f64::max)
    )))
}
