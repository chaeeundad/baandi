//! R2 (`docs/roadmap/rpa-crpa-plan.md`, ticket R2): dielectric response of an
//! insulator in the optical limit q -> 0. A fixed density (normally the
//! verified reference SCF density) is solved for `bands` bands on the
//! supplied time-reversal-symmetric k set by `refscf::nscf` (dense EVD by
//! default, every band residual-certified); chi0(q -> 0, i w) with head and
//! wings (velocity including the nonlocal commutator) is accumulated by the
//! sum over states of `dft-response`; the macroscopic dielectric tensor is
//! formed with and without local fields at w = 0 and on a mapped
//! Gauss-Legendre imaginary-frequency grid (decisions D1, D5).
//!
//! MPI level 1 (plan 3.4, decision D7): the NSCF keeps its owned-row layout;
//! for each k the eigenvectors are gathered (a permutation, identical bits on
//! every rank) and the transition units (k, block of empty bands) are dealt
//! round robin to ranks in a fixed order; every rank accumulates a partial
//! replicated chi0 for all frequencies, combined by `tree_sum_f64` (fixed
//! association). One rank is bit-identical to the threads-only run; other rank
//! counts agree to roundoff. Metals reject (decision D6); spin-polarized
//! densities reject (R3).
use crate::{
    model::ScfError,
    refscf::nscf::{BandOptions, BandSolver, Bands, NoStore, nscf_bands_with},
    refstart::StartInput,
};
use dft_geometry::kpoints::KPoint;
use dft_mpi::Collective;
use dft_numerics::linalg::ComplexMatrix;
use dft_physics::data::{Array3, Complex64};
use dft_physics::hamiltonian::Species;
use dft_response::{
    chi0::Accumulator,
    dielectric::{OpticalEpsilon, optical},
    frequency::{FrequencyGrid, auto_omega0},
    pair::PairFft,
    transitions::{KOrbitals, Transitions, Velocity},
};
use dft_units::Hartree;

pub mod rpa;
pub mod subspace;
pub mod wannier;

/// Settings of [`dielectric_with`] (Hartree atomic units).
#[derive(Clone, Copy, Debug)]
pub struct DielectricOptions {
    /// Response sphere |G|^2 / 2 <= ecut_chi (Ha), at most 4 x the
    /// wavefunction cutoff (the density sphere).
    pub ecut_chi_hartree: f64,
    /// Bands per k (occupied + empty); `None` = N_G(ecut_chi), clamped to the
    /// smallest k basis (plan 3.6 "auto").
    pub bands: Option<usize>,
    pub solver: BandSolver,
    /// Residual certificate of every NSCF band (Ha).
    pub eigen_tolerance: f64,
    /// Gauss-Legendre points of the imaginary-frequency grid.
    pub frequency_points: usize,
    /// Mapping scale (Ha); `None` = mean of the lowest direct transition and
    /// the largest included transition (plan 2.4).
    pub omega0_hartree: Option<f64>,
    /// Empty bands per transition unit (the MPI work unit).
    pub band_block: usize,
    /// R2-VEL-CHECK diagnostic: which k points get a gauge-invariant dipole
    /// record ([`DielectricResult::dipoles`]); default none.
    pub dipole_dump: DipoleDump,
}

/// k points of the diagnostic dipole record.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DipoleDump {
    Off,
    /// Index 0 (Gamma of a Gamma-centred mesh) and index 1 (the first other k).
    GammaAndFirst,
    All,
}

/// Bands whose consecutive eigenvalues differ by less than this (Ha) form one
/// degenerate group of the dipole record.
pub const DIPOLE_GROUP_TOLERANCE: f64 = 1e-5;

/// One (occupied group V, empty group C) pair of the dipole record at one k.
/// Band indices 0-based, inclusive. `velocity[a] = sum_{v in V, c in C}
/// |<u_v| dH_k/dk_a |u_c>|^2` (Ha^2 Bohr^2, the full velocity with the
/// nonlocal commutator); `position[a] = sum |<u_v|dH/dk_a|u_c>|^2 /
/// (e_c - e_v)^2` = `sum |<v| r_a |c>|^2` (Bohr^2); `*_local` with the
/// kinetic velocity k + G only (no [V_NL, r]).
#[derive(Clone, Debug)]
pub struct DipoleGroupPair {
    pub v: (usize, usize),
    pub c: (usize, usize),
    /// Group mean eigenvalues (Ha).
    pub e_v: f64,
    pub e_c: f64,
    /// The empty group ends at the last computed band (it may be incomplete).
    pub c_at_band_limit: bool,
    pub velocity: [f64; 3],
    pub position: [f64; 3],
    pub velocity_local: [f64; 3],
    pub position_local: [f64; 3],
}

/// Dipole record of one k point.
#[derive(Clone, Debug)]
pub struct KDipoles {
    pub k_index: usize,
    pub k_fractional: [f64; 3],
    pub k_cartesian: [f64; 3],
    pub pairs: Vec<DipoleGroupPair>,
}

/// Consecutive runs of `e[range]` with spacing below the group tolerance.
fn groups(e: &[f64], range: std::ops::Range<usize>) -> Vec<(usize, usize)> {
    let mut out: Vec<(usize, usize)> = Vec::new();
    for i in range {
        match out.last_mut() {
            Some(g) if e[i] - e[g.1] < DIPOLE_GROUP_TOLERANCE => g.1 = i,
            _ => out.push((i, i)),
        }
    }
    out
}

/// R2-GAMMA-VEL diagnostics of the velocity at k + G = 0 (switch
/// `RESPONSE_VELOCITY`): `no-g0` zeroes the projector k-gradient row there
/// (only the l = 1 columns are nonzero, `beta_1'(0) Y_1m(e_b)`); `no-g0-all`
/// also zeroes the projector values of that row (l = 0 only), i.e. the
/// commutator without the k + G = 0 plane wave at all. Both are physically
/// wrong (the FD checks fail); they model other codes' treatments. Other
/// values: no change.
fn drop_k_plus_g_zero(
    basis: &dft_basis::PlaneWaveBasis,
    projectors: &mut dft_physics::stress::ProjectorKGradient,
    terms: &str,
) {
    if terms != "no-g0" && terms != "no-g0-all" {
        return;
    }
    for (ig, g) in basis.momenta().iter().enumerate() {
        if g.norm_squared == 0. {
            if terms == "no-g0-all" {
                for c in 0..projectors.values.cols() {
                    projectors.values[(ig, c)] = Complex64::default();
                }
            }
            for m in projectors.gradients.iter_mut() {
                for c in 0..m.cols() {
                    m[(ig, c)] = Complex64::default();
                }
            }
        }
    }
}

/// The dipole record of one k (diagnostic; outside the chi0 accumulation).
#[allow(clippy::too_many_arguments)]
fn dipole_record(
    k_index: usize,
    basis: &dft_basis::PlaneWaveBasis,
    volume: f64,
    species: &[Species<'_>],
    full: &ComplexMatrix,
    energies: &[f64],
    occupied: usize,
    bands: usize,
    velocity_terms: &str,
) -> Result<KDipoles, ScfError> {
    let n = full.rows();
    let pick = |r: std::ops::Range<usize>| -> Result<ComplexMatrix, ScfError> {
        let mut m = ComplexMatrix::zeros(n, r.len()).map_err(numerical)?;
        for (j, b) in r.enumerate() {
            m.column_mut(j).copy_from_slice(full.column(b));
        }
        Ok(m)
    };
    let (occ, emp) = (pick(0..occupied)?, pick(occupied..bands)?);
    let momenta: Vec<[f64; 3]> = basis.momenta().iter().map(|g| g.cartesian).collect();
    let mut projectors =
        dft_physics::stress::projector_k_gradient(basis, volume, species).map_err(numerical)?;
    drop_k_plus_g_zero(basis, &mut projectors, velocity_terms);
    let full_p = dft_response::velocity::velocity(&momenta, Some(&projectors), &occ, &emp)
        .map_err(numerical)?;
    let local_p =
        dft_response::velocity::velocity(&momenta, None, &occ, &emp).map_err(numerical)?;
    let mean = |g: (usize, usize)| energies[g.0..=g.1].iter().sum::<f64>() / (g.1 - g.0 + 1) as f64;
    let mut pairs = Vec::new();
    for v in groups(energies, 0..occupied) {
        for c in groups(energies, occupied..bands) {
            let mut pair = DipoleGroupPair {
                v,
                c,
                e_v: mean(v),
                e_c: mean(c),
                c_at_band_limit: c.1 + 1 == bands,
                velocity: [0.; 3],
                position: [0.; 3],
                velocity_local: [0.; 3],
                position_local: [0.; 3],
            };
            for a in 0..3 {
                for i in v.0..=v.1 {
                    for j in c.0..=c.1 {
                        let d2 = (energies[j] - energies[i]).powi(2);
                        let f = full_p[a][(i, j - occupied)].norm_sqr();
                        let l = local_p[a][(i, j - occupied)].norm_sqr();
                        pair.velocity[a] += f;
                        pair.position[a] += f / d2;
                        pair.velocity_local[a] += l;
                        pair.position_local[a] += l / d2;
                    }
                }
            }
            pairs.push(pair);
        }
    }
    Ok(KDipoles {
        k_index,
        k_fractional: basis.k_fractional(),
        k_cartesian: basis.k_cartesian(),
        pairs,
    })
}

impl DielectricOptions {
    pub fn new(ecut_chi_hartree: f64) -> Self {
        Self {
            ecut_chi_hartree,
            bands: None,
            solver: BandSolver::Dense,
            eigen_tolerance: 1e-6,
            frequency_points: 16,
            omega0_hartree: None,
            band_block: 32,
            dipole_dump: DipoleDump::Off,
        }
    }
}

#[derive(Clone, Debug)]
pub struct DielectricResult {
    /// Response sphere size N_G (G = 0 included).
    pub response_size: usize,
    pub bands: usize,
    pub occupied: usize,
    pub electrons: f64,
    pub k_points: usize,
    /// Smallest and largest plane-wave basis over k.
    pub basis_sizes: (usize, usize),
    /// Evaluation frequencies (Ha): 0, then the Gauss-Legendre nodes.
    pub frequencies: Vec<f64>,
    /// Quadrature weights of `int_0^inf dw` for each frequency (0 for w = 0).
    pub weights: Vec<f64>,
    pub omega0: f64,
    /// Head tensor L (Bohr^-1 Ha^-1, `chi0_00 = q^2 q^.L.q^`) per frequency.
    pub head: Vec<[[Complex64; 3]; 3]>,
    /// Macroscopic tensors per frequency.
    pub epsilon: Vec<OpticalEpsilon>,
    /// Band edges (Ha) and transitions of the NSCF spectrum.
    pub valence_maximum: f64,
    pub conduction_minimum: f64,
    pub lowest_direct_gap: f64,
    pub largest_transition: f64,
    pub max_residual: f64,
    pub transitions: u64,
    pub units: usize,
    /// NSCF eigenvalues (Ha) per k.
    pub eigenvalues: Vec<Vec<f64>>,
    /// Diagnostic dipole records (`options.dipole_dump`), k order; every rank.
    pub dipoles: Vec<KDipoles>,
}

/// [`dielectric_with`] on the registered communicator (or serially).
pub fn dielectric(
    input: StartInput<'_>,
    density: &[Array3<f64>],
    points: &[KPoint],
    options: DielectricOptions,
) -> Result<DielectricResult, ScfError> {
    if let Some(comm) = crate::parallel::SharedComm::active() {
        dielectric_with(input, density, points, options, &comm)
    } else {
        dielectric_with(
            input,
            density,
            points,
            options,
            &crate::refscf::metric::Serial,
        )
    }
}

fn invalid(message: &'static str) -> ScfError {
    ScfError::InvalidInput(message)
}
fn numerical(e: impl std::fmt::Display) -> ScfError {
    ScfError::Numerical(format!("response: {e}"))
}

/// Every rank's owned rows of `bands` in global plane-wave order (`n_pw`
/// rows): a permutation of the owned values, identical bits on every rank.
fn gather<K: Collective>(comm: &K, bands: &Bands, n_pw: usize) -> Result<ComplexMatrix, ScfError> {
    let nb = bands.eigenvectors.cols();
    let width = 1 + 2 * nb;
    let mut lens = vec![0.; comm.size()];
    comm.allgatherv_f64(&[bands.rows.len() as f64], &vec![1; comm.size()], &mut lens)
        .map_err(numerical)?;
    let counts: Vec<usize> = lens.iter().map(|&l| l as usize * width).collect();
    if counts.iter().sum::<usize>() != n_pw * width {
        return Err(numerical("eigenvector row layout"));
    }
    let mut local = Vec::with_capacity(bands.rows.len() * width);
    for (i, &g) in bands.rows.iter().enumerate() {
        local.push(g as f64);
        for b in 0..nb {
            let z = bands.eigenvectors[(i, b)];
            local.extend([z.re, z.im]);
        }
    }
    let mut all = vec![0.; n_pw * width];
    comm.allgatherv_f64(&local, &counts, &mut all)
        .map_err(numerical)?;
    let mut full = ComplexMatrix::zeros(n_pw, nb).map_err(numerical)?;
    let mut seen = vec![false; n_pw];
    for row in all.chunks_exact(width) {
        let g = row[0] as usize;
        if g >= n_pw || seen[g] {
            return Err(numerical("eigenvector row ownership"));
        }
        seen[g] = true;
        for b in 0..nb {
            full[(g, b)] = Complex64::new(row[1 + 2 * b], row[2 + 2 * b]);
        }
    }
    Ok(full)
}

/// Agree a rank-local outcome before the next collective.
fn agree<K: Collective>(comm: &K, local: Result<(), ScfError>) -> Result<(), ScfError> {
    let mut flag = [0.];
    comm.tree_sum_f64(&[f64::from(local.is_err())], &mut flag)
        .map_err(numerical)?;
    match local {
        Err(e) => Err(e),
        Ok(()) if flag[0] != 0. => Err(numerical("a transition unit failed on another rank")),
        Ok(()) => Ok(()),
    }
}

/// Optical (q -> 0) dielectric response; collective (all arguments
/// identical on every rank). `density` must be one spin-unpolarized field;
/// `points` a time-reversal-symmetric k set (normally the full
/// Gamma-centred mesh) with positive weights.
pub fn dielectric_with<K: Collective>(
    input: StartInput<'_>,
    density: &[Array3<f64>],
    points: &[KPoint],
    options: DielectricOptions,
    comm: &K,
) -> Result<DielectricResult, ScfError> {
    if density.len() != 1 {
        return Err(invalid(
            "dielectric: spin-polarized response is not implemented (R3)",
        ));
    }
    let cutoff = input.options.cutoff_hartree;
    let ecut_chi = options.ecut_chi_hartree;
    if !ecut_chi.is_finite() || ecut_chi <= 0. || ecut_chi > 4. * cutoff {
        return Err(invalid(
            "dielectric: ecut_chi must be positive and at most 4 x the wavefunction cutoff",
        ));
    }
    if options.frequency_points == 0 || options.band_block == 0 || points.is_empty() {
        return Err(invalid("dielectric: frequency points, band block or k set"));
    }
    if points
        .iter()
        .any(|p| !p.weight.is_finite() || p.weight <= 0.)
    {
        return Err(invalid("dielectric: k weights must be positive"));
    }
    // Time reversal: -k must be in the set (mod G) with the same weight.
    for p in points {
        let partner = points.iter().any(|o| {
            (0..3).all(|i| {
                let s = p.fractional[i] + o.fractional[i];
                (s - s.round()).abs() < 1e-9
            }) && (o.weight - p.weight).abs() <= 1e-12 * p.weight
        });
        if !partner {
            return Err(invalid(
                "dielectric: the k set must be time-reversal symmetric (use the full Gamma-centred mesh)",
            ));
        }
    }
    let cell = input.cell;
    let volume = cell.volume_bohr3();
    let response = dft_basis::PlaneWaveBasis::new(cell, Hartree(ecut_chi)).map_err(numerical)?;
    let n_g = response.vectors().len();
    if response.vectors()[0].norm_squared != 0. {
        return Err(numerical("response sphere must start with G = 0"));
    }
    let mut basis_sizes = (usize::MAX, 0);
    for p in points {
        let b = dft_basis::PlaneWaveBasis::at_k(cell, Hartree(cutoff), p.fractional)
            .map_err(numerical)?;
        basis_sizes.0 = basis_sizes.0.min(b.vectors().len());
        basis_sizes.1 = basis_sizes.1.max(b.vectors().len());
    }
    let bands = options.bands.unwrap_or(n_g.min(basis_sizes.0));
    let band_options = BandOptions::new(bands, options.eigen_tolerance, options.solver);
    let threads = input.options.threads.max(1);
    let (atoms, potentials) = (input.atoms, input.potentials);
    let nscf = nscf_bands_with(input, density, points, band_options, comm, &mut NoStore)?;
    let half = nscf.electrons / 2.;
    if (half - half.round()).abs() > 1e-8 || half.round() < 1. {
        return Err(ScfError::OddElectronCount(nscf.electrons));
    }
    let occupied = half.round() as usize;
    if bands <= occupied {
        return Err(ScfError::InsufficientBands { bands, occupied });
    }
    let spectra = &nscf.channels[0];
    let (mut vbm, mut cbm) = (f64::NEG_INFINITY, f64::INFINITY);
    let (mut direct, mut largest, mut max_residual) = (f64::INFINITY, 0_f64, 0_f64);
    for b in spectra {
        let e = &b.eigenvalues;
        vbm = vbm.max(e[occupied - 1]);
        cbm = cbm.min(e[occupied]);
        direct = direct.min(e[occupied] - e[occupied - 1]);
        largest = largest.max(e[bands - 1] - e[0]);
        max_residual = b.residuals.iter().copied().fold(max_residual, f64::max);
    }
    // A numerical zero for touching bands, not a material threshold.
    if cbm - vbm <= 1e-6 {
        return Err(ScfError::NonInsulating {
            gap_hartree: cbm - vbm,
            required_hartree: 1e-6,
        });
    }
    let omega0 = match options.omega0_hartree {
        Some(w) => w,
        None => auto_omega0(direct, largest).map_err(numerical)?,
    };
    let grid = FrequencyGrid::mapped_gauss_legendre(options.frequency_points, omega0)
        .map_err(numerical)?;
    let mut frequencies = vec![0.];
    frequencies.extend(&grid.points);
    let mut weights = vec![0.];
    weights.extend(&grid.weights);
    let dims = nscf.grid_dims;
    for g in response.vectors() {
        dft_basis::FftGridSpec::map_miller(g.miller, dims).map_err(|_| {
            invalid("dielectric: the response sphere does not fit the density grid")
        })?;
    }
    // Species for the projector k gradients (Cartesian Bohr positions).
    let positions: Vec<Vec<[f64; 3]>> = potentials
        .iter()
        .map(|p| {
            atoms
                .iter()
                .filter(|a| a.species == p.label)
                .map(|a| cell.fractional_to_cartesian_bohr(a.fractional_position))
                .collect()
        })
        .collect();
    let species: Vec<Species<'_>> = potentials
        .iter()
        .zip(&positions)
        .filter(|(_, r)| !r.is_empty())
        .map(|(p, r)| Species {
            potential: p.potential,
            positions_bohr: r,
        })
        .collect();
    // R2-VEL-CHECK diagnostic (switch RESPONSE_VELOCITY, default full): which
    // part of d H_k / d k enters the dipoles. `column_l` is the angular
    // momentum of every projector column of `projector_k_gradient` (species,
    // atom, projector, m order).
    let velocity_terms =
        dft_numerics::switches::choice(dft_numerics::switches::Id::RESPONSE_VELOCITY);
    let column_l: Vec<usize> = species
        .iter()
        .flat_map(|s| {
            let per_atom: Vec<usize> = s
                .potential
                .projectors()
                .iter()
                .flat_map(|p| std::iter::repeat_n(p.l(), 2 * p.l() + 1))
                .collect();
            std::iter::repeat_n(per_atom, s.positions_bohr.len()).flatten()
        })
        .collect();
    let response_miller: Vec<[i32; 3]> = response.vectors().iter().map(|g| g.miller).collect();
    let mut fft = PairFft::new(dims, &response_miller).map_err(numerical)?;
    let mut accumulator = Accumulator::new(n_g, &frequencies, true, threads).map_err(numerical)?;
    let total_weight: f64 = points.iter().map(|p| p.weight).sum();
    let blocks: Vec<std::ops::Range<usize>> = (occupied..bands)
        .step_by(options.band_block)
        .map(|m| m..(m + options.band_block).min(bands))
        .collect();
    let mut unit = 0;
    let mut dipoles = Vec::new();
    for (k, (spectrum, basis)) in spectra.iter().zip(&nscf.bases).enumerate() {
        let n_pw = basis.vectors().len();
        let full = gather(comm, spectrum, n_pw)?;
        let mine: Vec<_> = blocks
            .iter()
            .enumerate()
            .filter(|(i, _)| (unit + i) % comm.size() == comm.rank())
            .map(|(_, r)| r.clone())
            .collect();
        unit += blocks.len();
        let local = (|| -> Result<(), ScfError> {
            if mine.is_empty() {
                return Ok(());
            }
            let miller: Vec<[i32; 3]> = basis.vectors().iter().map(|g| g.miller).collect();
            let momenta: Vec<[f64; 3]> = basis.momenta().iter().map(|g| g.cartesian).collect();
            let mut projectors = dft_physics::stress::projector_k_gradient(basis, volume, &species)
                .map_err(numerical)?;
            if velocity_terms == "sp" {
                // Zero the l >= 2 columns: D couples only equal (l, m), so this
                // removes exactly their terms from d V_NL / d k.
                for (c, &l) in column_l.iter().enumerate() {
                    if l >= 2 {
                        for m in std::iter::once(&mut projectors.values)
                            .chain(projectors.gradients.iter_mut())
                        {
                            m.column_mut(c).fill(Complex64::default());
                        }
                    }
                }
            }
            drop_k_plus_g_zero(basis, &mut projectors, velocity_terms);
            let orbitals = KOrbitals {
                miller: &miller,
                coefficients: &full,
                energies: &spectrum.eigenvalues,
            };
            let mut transitions =
                Transitions::new(dims, orbitals, occupied, &miller, [0; 3]).map_err(numerical)?;
            let scale = 2. * points[k].weight / (total_weight * volume);
            for block in mine {
                transitions
                    .add(
                        &mut accumulator,
                        &mut fft,
                        orbitals,
                        block,
                        scale,
                        Some(Velocity {
                            momenta: &momenta,
                            projectors: (velocity_terms != "local").then_some(&projectors),
                        }),
                    )
                    .map_err(numerical)?;
            }
            Ok(())
        })();
        agree(comm, local)?;
        let dump = match options.dipole_dump {
            DipoleDump::Off => false,
            DipoleDump::GammaAndFirst => k < 2,
            DipoleDump::All => true,
        };
        if dump {
            // Replicated data on every rank: same result (or error) everywhere.
            dipoles.push(dipole_record(
                k,
                basis,
                volume,
                &species,
                &full,
                &spectrum.eigenvalues,
                occupied,
                bands,
                velocity_terms,
            )?);
        }
    }
    let chi0 = accumulator.reduce(comm).map_err(numerical)?;
    let norms: Vec<f64> = response
        .vectors()
        .iter()
        .map(|g| g.norm_squared.sqrt())
        .collect();
    let epsilon = (0..frequencies.len())
        .map(|w| optical(&chi0, w, &norms).map_err(numerical))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(DielectricResult {
        response_size: n_g,
        bands,
        occupied,
        electrons: nscf.electrons,
        k_points: points.len(),
        basis_sizes,
        frequencies,
        weights,
        omega0,
        head: chi0.head.clone(),
        epsilon,
        valence_maximum: vbm,
        conduction_minimum: cbm,
        lowest_direct_gap: direct,
        largest_transition: largest,
        max_residual,
        transitions: chi0.transitions,
        units: unit,
        eigenvalues: spectra.iter().map(|b| b.eigenvalues.clone()).collect(),
        dipoles,
    })
}
