//! Residual-certified spectra of a fixed supplied density. No SCF map, density
//! update, mixing, occupations solve or self-consistency certificate is run.
//!
//! R1 (`docs/roadmap/rpa-crpa-plan.md` 3.2, ticket R1): [`nscf_bands_with`]
//! adds two many-band solvers next to the reference Davidson. `Dense` builds
//! the explicit H_k (kinetic diagonal, local V(G - G') from one DFT of the
//! gathered potential, nonlocal P D P^H) and diagonalizes it completely with
//! the faer Hermitian EVD: the small-cell production path and the oracle of
//! the iterative path. `Chebyshev` is the filtered subspace iteration of
//! `dft_solvers::chebyshev_complex` with buffer columns (never certified,
//! discarded) and locking. Every returned band carries a fresh
//! `||H psi - epsilon psi||` from an independent application, whatever the
//! solver; a per-k [`BandStore`] lets a host checkpoint completed k points and
//! resume them (resumed bands are re-certified, never trusted).
use crate::{
    Functional,
    model::ScfError,
    refstart::{self, KStartState, StartInput},
};
use dft_geometry::{Cell, kpoints::KPoint};
use dft_mpi::Collective;
use dft_numerics::linalg::ComplexMatrix;
use dft_numerics::linalg::hermitian_eigh;
use dft_physics::data::Array3;
use dft_physics::data::Complex64;
use dft_physics::distributed_complex::DistributedComplexHamiltonian;
use dft_solvers::{chebyshev_complex, ref_davidson, ref_davidson_complex};

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
    pub solver: BandSolver,
    /// Channel-major then k-major, like `channels`.
    pub stats: Vec<Vec<BandStats>>,
}

/// Eigensolver of the fixed-density band problem.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BandSolver {
    /// Reference block Davidson, every band strict (the default band path).
    Davidson,
    /// Explicit H_k and a complete Hermitian EVD (all eigenpairs).
    Dense,
    /// Chebyshev-filtered subspace iteration with buffer bands and locking.
    Chebyshev,
}
impl BandSolver {
    pub fn name(self) -> &'static str {
        match self {
            Self::Davidson => "davidson",
            Self::Dense => "dense",
            Self::Chebyshev => "chebyshev",
        }
    }
    pub fn parse(name: &str) -> Option<Self> {
        [Self::Davidson, Self::Dense, Self::Chebyshev]
            .into_iter()
            .find(|s| s.name() == name)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct BandOptions {
    /// Certified bands per channel and k.
    pub bands: usize,
    /// Fresh-residual certificate (Ha) of every returned band.
    pub tolerance: f64,
    pub solver: BandSolver,
    /// Chebyshev: uncertified buffer columns; `None` = max(8, bands / 10).
    /// Clamped to the k basis size.
    pub buffer: Option<usize>,
    /// Chebyshev polynomial degree per sweep.
    pub degree: usize,
    /// Chebyshev sweeps before failure.
    pub max_iterations: usize,
}
impl BandOptions {
    /// Chebyshev defaults (degree 10, 300 sweeps) are initial values from the
    /// synthetic and Si fixtures of R1, not yet calibrated on the generality set.
    pub fn new(bands: usize, tolerance: f64, solver: BandSolver) -> Self {
        Self {
            bands,
            tolerance,
            solver,
            buffer: None,
            degree: 10,
            max_iterations: 300,
        }
    }
}

/// Per channel and k: basis size and solver work (zero when resumed).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct BandStats {
    pub basis_size: usize,
    pub resumed: bool,
    /// Chebyshev sweeps or Davidson calls; 1 for the dense EVD.
    pub iterations: usize,
    /// H columns applied by the solver (dense: 0), certificate excluded.
    pub h_columns: u64,
}

/// Host persistence of completed k points (checkpoint and resume). Every
/// method is collective: all ranks call it in the same order with the same
/// channel/k. `load` returns this rank's owned rows (row order of `rows`)
/// and the eigenvalues of a stored unit, `None` when the unit is absent, and
/// an error when a stored unit does not match the run. The engine agrees the
/// presence over ranks and re-certifies loaded bands with a fresh application.
pub trait BandStore<K: Collective> {
    fn load(
        &mut self,
        comm: &K,
        channel: usize,
        k: usize,
        rows: &[usize],
        bands: usize,
    ) -> Result<Option<(Vec<f64>, ComplexMatrix)>, ScfError>;
    fn save(
        &mut self,
        comm: &K,
        channel: usize,
        k: usize,
        global_len: usize,
        bands: &Bands,
    ) -> Result<(), ScfError>;
}
/// No persistence.
pub struct NoStore;
impl<K: Collective> BandStore<K> for NoStore {
    fn load(
        &mut self,
        _: &K,
        _: usize,
        _: usize,
        _: &[usize],
        _: usize,
    ) -> Result<Option<(Vec<f64>, ComplexMatrix)>, ScfError> {
        Ok(None)
    }
    fn save(&mut self, _: &K, _: usize, _: usize, _: usize, _: &Bands) -> Result<(), ScfError> {
        Ok(())
    }
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
    nscf_bands_with(
        input,
        density,
        points,
        BandOptions::new(bands, eigen_tolerance, BandSolver::Davidson),
        comm,
        &mut NoStore,
    )
}

/// [`nscf_bands_with`] on the registered communicator (or serially).
pub fn nscf_bands<S: BandStore<crate::parallel::SharedComm> + BandStore<super::metric::Serial>>(
    input: StartInput<'_>,
    density: &[Array3<f64>],
    points: &[KPoint],
    options: BandOptions,
    store: &mut S,
) -> Result<NscfResult, ScfError> {
    if let Some(comm) = crate::parallel::SharedComm::active() {
        nscf_bands_with(input, density, points, options, &comm, store)
    } else {
        nscf_bands_with(
            input,
            density,
            points,
            options,
            &super::metric::Serial,
            store,
        )
    }
}

/// Fixed-density bands with a selectable solver. `options.bands` lowest
/// eigenpairs per channel and k, each with a fresh residual at most
/// `options.tolerance`. Every rank processes every k point in order. With the
/// Davidson solver this is exactly [`nscf_with`]. Dense and Chebyshev require
/// `bands` at most the basis size of every k (the complete basis is allowed).
pub fn nscf_bands_with<K: Collective, S: BandStore<K>>(
    input: StartInput<'_>,
    density: &[Array3<f64>],
    points: &[KPoint],
    band_options: BandOptions,
    comm: &K,
    store: &mut S,
) -> Result<NscfResult, ScfError> {
    let BandOptions {
        bands,
        tolerance: eigen_tolerance,
        solver,
        ..
    } = band_options;
    if !eigen_tolerance.is_finite()
        || eigen_tolerance <= 0.
        || bands == 0
        || band_options.degree == 0
        || band_options.max_iterations == 0
    {
        return Err(ScfError::InvalidInput(
            "reference NSCF: invalid bands/tolerance",
        ));
    }
    // The many-band solvers need the start only for the operator (and, for
    // Chebyshev, the occupied-like leading columns); the start block is
    // capped near the occupied count so that a large band request does not
    // run the atomic-start subspace solve over thousands of columns. A
    // request below occupied + 4 still reaches the start's own rejection.
    let start_bands = match solver {
        BandSolver::Davidson => bands,
        BandSolver::Dense | BandSolver::Chebyshev => {
            let electrons: f64 = input
                .atoms
                .iter()
                .filter_map(|a| {
                    input
                        .potentials
                        .iter()
                        .find(|p| p.label == a.species)
                        .map(|p| p.potential.z_valence())
                })
                .sum();
            if electrons.is_finite() && electrons >= 0. {
                bands.min(electrons.ceil() as usize + 4)
            } else {
                bands
            }
        }
    };
    let mut options = input.options.clone();
    options.bands = Some(start_bands);
    let (mut state, potentials) = refstart::start_kpoints_from_density(
        StartInput {
            options: &options,
            ..input
        },
        density,
        points,
        comm,
    )?;
    if state.blocks.iter().any(|b| b.operator.global_len() < bands) {
        return Err(ScfError::InvalidInput(
            "reference NSCF: too many bands for a k basis",
        ));
    }
    let dims = state.grid.dims();
    let mut channels = Vec::with_capacity(potentials.len());
    let mut stats = Vec::with_capacity(potentials.len());
    for (channel, potential) in potentials.iter().enumerate() {
        let mut average = [0.];
        comm.tree_sum_f64(&[potential.iter().sum::<f64>()], &mut average)
            .map_err(super::numerical)?;
        let v0 = average[0] / dims.iter().product::<usize>() as f64;
        let mut spectra = Vec::with_capacity(points.len());
        let mut channel_stats = Vec::with_capacity(points.len());
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
            let context = |e: ScfError| {
                ScfError::Numerical(format!(
                    "reference NSCF channel {channel}, k index {k}: {e}"
                ))
            };
            let basis_size = op.global_len();
            let loaded = store
                .load(comm, channel, k, op.local_rows(), bands)
                .and_then(|l| {
                    let ok = l.as_ref().is_none_or(|(values, vectors)| {
                        values.len() == bands
                            && values.iter().all(|v| v.is_finite())
                            && vectors.rows() == op.local_rows().len()
                            && vectors.cols() == bands
                    });
                    let mut flags = [0.; 2];
                    comm.tree_sum_f64(&[f64::from(l.is_some()), f64::from(!ok)], &mut flags)
                        .map_err(super::numerical)?;
                    if flags[1] != 0. || (flags[0] != 0. && flags[0] != comm.size() as f64) {
                        return Err(ScfError::InvalidInput(
                            "reference NSCF: stored bands disagree with this run",
                        ));
                    }
                    Ok(l)
                })
                .map_err(context)?;
            let (result, stat) = if let Some((eigenvalues, eigenvectors)) = loaded {
                let residuals = certify(comm, op, &eigenvectors, &eigenvalues).map_err(context)?;
                if let Some(r) = residuals
                    .iter()
                    .find(|r| r.is_nan() || **r > eigen_tolerance)
                {
                    return Err(context(ScfError::Numerical(format!(
                        "stored bands fail the residual certificate: {r:e} Ha > {eigen_tolerance:e} Ha"
                    ))));
                }
                (
                    Bands {
                        eigenvalues,
                        residuals,
                        eigenvectors,
                        rows: op.local_rows().to_vec(),
                    },
                    BandStats {
                        basis_size,
                        resumed: true,
                        ..Default::default()
                    },
                )
            } else {
                let diagonal: Vec<_> = op.diagonal_static().iter().map(|v| v + v0).collect();
                let solved = match solver {
                    BandSolver::Davidson => solve(
                        comm,
                        op,
                        &state.orbitals.0[k],
                        &diagonal,
                        bands,
                        eigen_tolerance,
                    ),
                    BandSolver::Dense => solve_dense(comm, op, bands, eigen_tolerance),
                    BandSolver::Chebyshev => {
                        solve_chebyshev(comm, op, &state.orbitals.0[k], k, &band_options)
                    }
                }
                .map_err(context)?;
                let (mut bands_k, mut stat) = solved;
                stat.basis_size = basis_size;
                if solver != BandSolver::Davidson {
                    bands_k.residuals =
                        certify(comm, op, &bands_k.eigenvectors, &bands_k.eigenvalues)
                            .map_err(context)?;
                    if let Some(r) = bands_k
                        .residuals
                        .iter()
                        .find(|r| r.is_nan() || **r > eigen_tolerance)
                    {
                        return Err(context(ScfError::Numerical(format!(
                            "{} bands fail the residual certificate: {r:e} Ha > {eigen_tolerance:e} Ha",
                            solver.name()
                        ))));
                    }
                }
                store
                    .save(comm, channel, k, basis_size, &bands_k)
                    .map_err(context)?;
                (bands_k, stat)
            };
            spectra.push(result);
            channel_stats.push(stat);
        }
        channels.push(spectra);
        stats.push(channel_stats);
    }
    Ok(NscfResult {
        points: points.to_vec(),
        channels,
        electrons: state.electrons,
        tolerance_hartree: eigen_tolerance,
        solver,
        stats,
    })
}

/// Operator columns per application in the many-band paths.
const CHUNK: usize = 64;

/// Fresh `||H x_j - epsilon_j x_j||` (Ha) from an independent application in
/// column chunks.
fn certify<K: Collective>(
    comm: &K,
    op: &mut DistributedComplexHamiltonian,
    vectors: &ComplexMatrix,
    eigenvalues: &[f64],
) -> Result<Vec<f64>, ScfError> {
    op.ensure_capacity(CHUNK).map_err(super::numerical)?;
    let chunk = op.capacity();
    let rows = vectors.rows();
    let mut local = vec![0.; vectors.cols()];
    let mut start = 0;
    while start < vectors.cols() {
        let end = (start + chunk).min(vectors.cols());
        let mut x = ComplexMatrix::zeros(rows, end - start).map_err(super::numerical)?;
        for j in start..end {
            x.column_mut(j - start).copy_from_slice(vectors.column(j));
        }
        let mut hx = ComplexMatrix::zeros(rows, end - start).map_err(super::numerical)?;
        op.apply_block(comm, &x, &mut hx)
            .map_err(super::numerical)?;
        for j in start..end {
            local[j] = x
                .column(j - start)
                .iter()
                .zip(hx.column(j - start))
                .map(|(a, b)| (b - eigenvalues[j] * a).norm_sqr())
                .sum();
        }
        start = end;
    }
    let mut residuals = vec![0.; local.len()];
    comm.tree_sum_f64(&local, &mut residuals)
        .map_err(super::numerical)?;
    Ok(residuals.into_iter().map(f64::sqrt).collect())
}

/// Complete EVD of the explicit H_k (identical bits on every rank, so each
/// rank keeps its own rows of the same eigenvectors).
fn solve_dense<K: Collective>(
    comm: &K,
    op: &mut DistributedComplexHamiltonian,
    bands: usize,
    _tolerance: f64,
) -> Result<(Bands, BandStats), ScfError> {
    let h = op.dense_matrix(comm).map_err(super::numerical)?;
    let (values, vectors) = hermitian_eigh(&h).map_err(super::numerical)?;
    drop(h);
    let rows = op.local_rows();
    let mut local = ComplexMatrix::zeros(rows.len(), bands).map_err(super::numerical)?;
    for j in 0..bands {
        for (i, &g) in rows.iter().enumerate() {
            local[(i, j)] = vectors[(g, j)];
        }
    }
    Ok((
        Bands {
            eigenvalues: values[..bands].to_vec(),
            residuals: Vec::new(),
            eigenvectors: local,
            rows: rows.to_vec(),
        },
        BandStats {
            iterations: 1,
            ..Default::default()
        },
    ))
}

/// SplitMix64 finalizer: a rank-independent pseudo-random value per
/// (global row, column, k).
fn mix(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}
fn uniform(seed: u64) -> f64 {
    (mix(seed) >> 11) as f64 / (1u64 << 53) as f64 * 2. - 1.
}

/// Start block: the start orbitals, then deterministic random columns keyed
/// by global row (rank-count independent) and damped by 1 / (1 + T_G) toward
/// low kinetic energy.
fn solve_chebyshev<K: Collective>(
    comm: &K,
    op: &mut DistributedComplexHamiltonian,
    start: &ComplexMatrix,
    k: usize,
    options: &BandOptions,
) -> Result<(Bands, BandStats), ScfError> {
    let bands = options.bands;
    let buffer = options.buffer.unwrap_or((bands / 10).max(8));
    let block = (bands + buffer).min(op.global_len());
    let rows = op.local_rows().to_vec();
    let kinetic = op.kinetic_diagonal().to_vec();
    let mut x = ComplexMatrix::zeros(rows.len(), block).map_err(super::numerical)?;
    let kept = start.cols().min(block);
    for j in 0..kept {
        x.column_mut(j).copy_from_slice(start.column(j));
    }
    for j in kept..block {
        for (i, &g) in rows.iter().enumerate() {
            let seed = ((g as u64) << 24 ^ j as u64) ^ ((k as u64) << 52);
            x[(i, j)] =
                Complex64::new(uniform(2 * seed), uniform(2 * seed + 1)) / (1. + kinetic[i]);
        }
    }
    let upper = op.spectral_upper_bound(comm).map_err(super::numerical)?;
    op.ensure_capacity(CHUNK).map_err(super::numerical)?;
    let solution = chebyshev_complex::solve(
        comm,
        op,
        &x,
        chebyshev_complex::Options {
            bands,
            tolerance: options.tolerance,
            degree: options.degree,
            max_iterations: options.max_iterations,
            upper_bound: upper,
            chunk: op.capacity(),
        },
    )
    .map_err(super::numerical)?;
    Ok((
        Bands {
            eigenvalues: solution.eigenvalues,
            residuals: solution.residuals,
            eigenvectors: solution.eigenvectors,
            rows,
        },
        BandStats {
            iterations: solution.iterations,
            h_columns: solution.h_columns,
            ..Default::default()
        },
    ))
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
) -> Result<(Bands, BandStats), ScfError> {
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
    for call in 1..=64 {
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
            return Ok((
                Bands {
                    eigenvalues: solution.eigenvalues,
                    residuals,
                    eigenvectors: vectors,
                    rows: op.local_rows().to_vec(),
                },
                BandStats {
                    iterations: call,
                    ..Default::default()
                },
            ));
        }
    }
    Err(ScfError::Numerical(format!(
        "unconverged band residual after 64 calls: max {:e} Ha, requested {tolerance:e} Ha",
        residuals.iter().copied().fold(0., f64::max)
    )))
}
