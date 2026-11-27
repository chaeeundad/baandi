//! R3 (`docs/roadmap/rpa-crpa-plan.md`, ticket R3; ADR 0055): RPA
//! correlation energy of an insulator from the ACFDT on the imaginary axis,
//!
//! `E_c = sum_q w_q (1 / 2 pi) int_0^inf dw Tr[ln(1 - v chi0(q, i w)) + v chi0(q, i w)]`
//!
//! (Hartree per cell). Pipeline: fixed-density NSCF on the FULL Gamma-centred
//! k mesh (`refscf::nscf`, every band residual-certified), q on the same mesh
//! reduced by time reversal or by the mesh-preserving crystal group
//! (`reduce_mesh`, exact orbits; E_c(q) is invariant; with one |G| sphere
//! for every q, q is the first-zone image, `g_sphere_q_units`), chi0(q, i w) for every
//! irreducible q by the R2 sum over states with `k + q = k' + G0` umklapps,
//! `ln det` by Cholesky of `1 - v^1/2 chi0 v^1/2`, the q = 0 term as the
//! Lebedev direction average of the k.p head/wings limit
//! (`dft_response::acfdt`), and a cutoff series `E_cut^chi` evaluated from
//! one pass over the pair densities (each cutoff keeps its own accumulator on
//! its sub-sphere and its own band count, by default `N_b = N_G(E_cut^chi)`),
//! with the `(E_cut^chi)^{-3/2}` extrapolation.
//!
//! Spin: one density (factor 2 in chi0) or two collinear densities of a
//! fixed-moment insulator (`chi0 = chi0_up + chi0_down`, integer electron
//! count per channel, each channel gapped).
//!
//! MPI level 1 (D7): transition units (q, channel, k, band block) are dealt
//! round robin in a fixed order; per q every cutoff accumulator is reduced
//! by `tree_sum_f64`; the log det runs replicated (identical bits). One rank
//! is bit-identical to the threads-only run.
//!
//! Checkpoint: one [`QContribution`] per irreducible q through a host
//! [`QStore`] (the CLI writes JSON records, ADR 0055); a loaded record is
//! checked against the run (q, weight, frequencies, response sizes) and its
//! chi0 is not recomputed.
use super::{agree, gather, invalid, numerical};
use crate::{
    model::ScfError,
    refscf::nscf::{BandOptions, BandSolver, NoStore, nscf_bands_with},
    refstart::StartInput,
};
use dft_basis::PlaneWaveBasis;
use dft_geometry::{
    kpoints::KPoint,
    symmetry::{SpinAction, SymmetryOperation, reduce_mesh},
};
use dft_mpi::Collective;
use dft_numerics::linalg::ComplexMatrix;
use dft_physics::data::Array3;
use dft_physics::hamiltonian::Species;
use dft_response::{
    acfdt::{self, Extrapolation},
    chi0::Accumulator,
    frequency::{FrequencyGrid, auto_omega0},
    pair::PairFft,
    transitions::{KOrbitals, Transitions, Velocity},
};
use dft_units::Hartree;
use std::collections::HashMap;

/// Unit directions with weights summing to one.
type Directions = Vec<([f64; 3], f64)>;
/// (index, representative fractional q, weight, members) of one q orbit.
type QEntry = (usize, [f64; 3], f64, usize);

/// q-point reduction of the mesh sum.
#[derive(Clone, Debug)]
pub enum QSymmetry {
    /// Every mesh point is its own q (the oracle of the reductions).
    Full,
    /// q and -q merged (identity plus time reversal).
    TimeReversal,
    /// Orbits of the mesh-preserving subgroup of these crystal operations
    /// plus time reversal (`reduce_mesh`, exact for the mesh sum).
    Crystal(Vec<SymmetryOperation>),
}
impl QSymmetry {
    pub fn label(&self) -> &'static str {
        match self {
            Self::Full => "none",
            Self::TimeReversal => "time-reversal",
            Self::Crystal(_) => "crystal",
        }
    }
}

/// Settings of [`rpa_with`] (Hartree atomic units).
#[derive(Clone, Debug)]
pub struct RpaOptions {
    /// Response cutoffs E_cut^chi (Ha), ascending, each at most 4 x the
    /// wavefunction cutoff.
    pub ecut_chi_series: Vec<f64>,
    /// Bands per k for every cutoff; `None` = N_G(E_cut^chi) at q = 0 per
    /// cutoff, clamped to the smallest k basis (plan 2.5: bands follow the basis).
    pub bands: Option<usize>,
    pub solver: BandSolver,
    pub eigen_tolerance: f64,
    pub frequency_points: usize,
    pub omega0_hartree: Option<f64>,
    /// Empty bands per transition unit (the MPI work unit).
    pub band_block: usize,
    /// Lebedev order of the q -> 0 direction average (6, 14, 26 or 50).
    pub q0_directions: usize,
    /// Gamma-centred k mesh (= q mesh).
    pub k_mesh: [usize; 3],
    pub q_symmetry: QSymmetry,
    /// Response sphere at finite q: `false` (default) = {G : |q+G|^2/2 <=
    /// E_cut^chi} (ADR 0054); `true` = the q = 0 Miller set {|G|^2/2 <=
    /// E_cut^chi} for every q (oracle matching for codes that use one
    /// G sphere for all q; v still uses |q+G|). In this mode E_c(q)
    /// depends on the q vector, not only on q modulo the reciprocal
    /// lattice: each mesh q is evaluated at its first-Brillouin-zone
    /// image(s) (mean over ties on zone faces, `first_zone_images`), the
    /// stored `QContribution::fractional` is that evaluated vector, and the
    /// reductions group (q, image) units related by S or -S without G0
    /// (R3-QIMAGE, ADR 0055 item 2).
    pub sphere_at_gamma: bool,
}
impl RpaOptions {
    pub fn new(ecut_chi_series: Vec<f64>, k_mesh: [usize; 3]) -> Self {
        Self {
            ecut_chi_series,
            bands: None,
            solver: BandSolver::Dense,
            eigen_tolerance: 1e-6,
            frequency_points: 16,
            omega0_hartree: None,
            band_block: 32,
            q0_directions: 26,
            k_mesh,
            q_symmetry: QSymmetry::TimeReversal,
            sphere_at_gamma: false,
        }
    }
}

/// Per cutoff at one q.
#[derive(Clone, Debug, PartialEq)]
pub struct QCutoff {
    /// Response sphere size at this q.
    pub response_size: usize,
    /// `Tr[ln(1 - X) + X]` per frequency (q = 0: the headline direction average).
    pub terms: Vec<f64>,
    /// `(1 / 2 pi) sum_w weight_w terms_w` (Ha per cell, not q-weighted).
    pub energy: f64,
    /// Smallest squared Cholesky pivot of `1 - X` over the frequencies (>= 1).
    pub min_pivot: f64,
    /// Largest relative Hermiticity defect of X over the frequencies.
    pub hermiticity: f64,
    /// q = 0 only: energy with each direction set (label, Ha per cell):
    /// Lebedev 6, 14, 26, 50 and a 24 x 48 product rule.
    pub direction_energies: Vec<(String, f64)>,
    /// q = 0 only: energy of the G != 0 body alone (diagnostic).
    pub body_only_energy: Option<f64>,
}

/// One irreducible q (the checkpoint unit).
#[derive(Clone, Debug, PartialEq)]
pub struct QContribution {
    pub index: usize,
    /// Representative q (fractional reciprocal coordinates of the mesh):
    /// the mesh value in [0, 1) (q+G sphere; evaluated at its [-1/2, 1/2)
    /// wrap, equivalent there), or the evaluated first-zone image (G sphere).
    pub fractional: [f64; 3],
    /// Weight in the mesh sum (orbit size / N_q; G sphere: sum over the
    /// orbit's units of 1 / (N_q x images of the unit's mesh q)).
    pub weight: f64,
    /// Members of the orbit (mesh points, or (q, image) units in G mode).
    pub multiplicity: usize,
    pub optical: bool,
    pub transitions: u64,
    /// Mapping scale of the frequency grid used (Ha).
    pub omega0: f64,
    pub cutoffs: Vec<QCutoff>,
}

/// Host persistence of completed q (collective: all ranks call in order).
pub trait QStore<K: Collective> {
    fn load(&mut self, comm: &K, index: usize) -> Result<Option<QContribution>, ScfError>;
    fn save(&mut self, comm: &K, record: &QContribution) -> Result<(), ScfError>;
}
/// No persistence.
pub struct NoQStore;
impl<K: Collective> QStore<K> for NoQStore {
    fn load(&mut self, _: &K, _: usize) -> Result<Option<QContribution>, ScfError> {
        Ok(None)
    }
    fn save(&mut self, _: &K, _: &QContribution) -> Result<(), ScfError> {
        Ok(())
    }
}

#[derive(Clone, Debug)]
pub struct RpaResult {
    pub cutoffs: Vec<f64>,
    /// N_G at q = 0 per cutoff.
    pub response_sizes: Vec<usize>,
    /// Bands entering chi0 per cutoff.
    pub bands: Vec<usize>,
    /// Bands of the NSCF (the largest of `bands`).
    pub nscf_bands: usize,
    /// Occupied bands per spin channel.
    pub occupied: Vec<usize>,
    pub electrons: f64,
    pub atoms: usize,
    pub k_points: usize,
    pub basis_sizes: (usize, usize),
    /// FFT grid of the operators and pair densities (the density grid).
    pub grid_dims: [usize; 3],
    pub frequencies: Vec<f64>,
    pub weights: Vec<f64>,
    pub omega0: f64,
    pub q_reduction: &'static str,
    pub q_points: Vec<QContribution>,
    /// q points taken from the store.
    pub resumed: usize,
    /// E_c per cutoff (Ha per cell).
    pub energies: Vec<f64>,
    /// `sum_q w_q terms_w` per cutoff and frequency (dimensionless).
    pub frequency_terms: Vec<Vec<f64>>,
    /// Fit over the cutoff series in `(E_cut^chi)^{-3/2}` (at least two cutoffs).
    pub extrapolation: Option<Extrapolation>,
    /// Diagnostic fit of the same energies in `1 / N_G(q = 0)` (the discrete
    /// shell form of the same law; slope in Ha x N_G).
    pub extrapolation_inverse_ng: Option<Extrapolation>,
    pub valence_maximum: f64,
    pub conduction_minimum: f64,
    pub lowest_direct_gap: f64,
    pub largest_transition: f64,
    pub max_residual: f64,
    pub min_pivot: f64,
    pub max_hermiticity: f64,
    pub transitions: u64,
}

/// [`rpa_with`] on the registered communicator (or serially).
pub fn rpa<S: QStore<crate::parallel::SharedComm> + QStore<crate::refscf::metric::Serial>>(
    input: StartInput<'_>,
    density: &[Array3<f64>],
    options: &RpaOptions,
    store: &mut S,
) -> Result<RpaResult, ScfError> {
    if let Some(comm) = crate::parallel::SharedComm::active() {
        rpa_with(input, density, options, &comm, store)
    } else {
        rpa_with(
            input,
            density,
            options,
            &crate::refscf::metric::Serial,
            store,
        )
    }
}

/// Gamma-centred full mesh, i-major (index `(i n1 + j) n2 + l`).
pub fn gamma_mesh(n: [usize; 3]) -> Vec<KPoint> {
    let total = (n[0] * n[1] * n[2]) as f64;
    let mut points = Vec::with_capacity(n.iter().product());
    for i in 0..n[0] {
        for j in 0..n[1] {
            for l in 0..n[2] {
                points.push(KPoint {
                    fractional: [
                        i as f64 / n[0] as f64,
                        j as f64 / n[1] as f64,
                        l as f64 / n[2] as f64,
                    ],
                    weight: 1. / total,
                });
            }
        }
    }
    points
}

fn mesh_index(n: [usize; 3], fractional: [f64; 3]) -> Result<[usize; 3], ScfError> {
    let mut out = [0; 3];
    for i in 0..3 {
        let x = fractional[i] * n[i] as f64;
        if (x - x.round()).abs() > 1e-9 {
            return Err(numerical("q is not a point of the k mesh"));
        }
        out[i] = (x.round() as i64).rem_euclid(n[i] as i64) as usize;
    }
    Ok(out)
}
fn flat(n: [usize; 3], i: [usize; 3]) -> usize {
    (i[0] * n[1] + i[1]) * n[2] + i[2]
}

/// Per k plane-wave data shared by all q.
struct KData {
    miller: Vec<[i32; 3]>,
    momenta: Vec<[f64; 3]>,
    wrapped: [f64; 3],
    radius: f64,
    lo: [i32; 3],
    hi: [i32; 3],
}

/// Alias guard of the pair products at (k, k + q): every product component
/// P = (M - G0) - G1 and response G satisfy P - G = n N (grid multiples) only
/// for n = 0. A nonzero n is excluded when its box range is empty or when the
/// grid superlattice vector is longer than the largest |q+P| + |q+G|.
fn alias_free(
    reciprocal: &[[f64; 3]; 3],
    dims: [usize; 3],
    left: &KData,
    right: &KData,
    shift: [i32; 3],
    response: (&[i32; 3], &[i32; 3], f64),
) -> bool {
    let (rlo, rhi, rr) = response;
    let mut range = [(0_i64, 0_i64); 3];
    for i in 0..3 {
        let plo = i64::from(right.lo[i] - shift[i] - left.hi[i]);
        let phi = i64::from(right.hi[i] - shift[i] - left.lo[i]);
        let n = dims[i] as i64;
        let lo = plo - i64::from(rhi[i]);
        let hi = phi - i64::from(rlo[i]);
        range[i] = (
            lo.div_euclid(n) + i64::from(lo.rem_euclid(n) != 0),
            hi.div_euclid(n),
        );
    }
    let reach = left.radius + right.radius + rr;
    for a in range[0].0..=range[0].1 {
        for b in range[1].0..=range[1].1 {
            for c in range[2].0..=range[2].1 {
                if a == 0 && b == 0 && c == 0 {
                    continue;
                }
                let m = [
                    a as f64 * dims[0] as f64,
                    b as f64 * dims[1] as f64,
                    c as f64 * dims[2] as f64,
                ];
                let l: f64 = (0..3)
                    .map(|i| (0..3).map(|j| reciprocal[i][j] * m[j]).sum::<f64>().powi(2))
                    .sum::<f64>()
                    .sqrt();
                if l <= reach * (1. + 1e-12) {
                    return false;
                }
            }
        }
    }
    true
}

fn identity(atoms: usize) -> SymmetryOperation {
    SymmetryOperation {
        rotation: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        translation: [0.; 3],
        atom_permutation: (0..atoms).collect(),
        spin_action: SpinAction::Identity,
        antiunitary: false,
        model_invariance_proof: "identity",
    }
}

/// Relative tolerance of equal |q| among the images of one mesh q (a
/// floating-point zero: images of a mesh point on a zone face have exactly
/// equal length; other image lengths differ by O(|B|^2 / N_k)).
const IMAGE_TIE: f64 = 1e-10;

/// First-Brillouin-zone images of the mesh point `index` (integer mesh
/// coordinates, `fractional = index / mesh`): every q = index / mesh - G0
/// of minimal |q| (Cartesian metric of the cell), as integer numerators
/// over `mesh`, ascending lexicographically. One image inside the zone;
/// two or more on a zone face, edge or corner (ties).
pub fn first_zone_images(
    cell: &dft_geometry::Cell,
    mesh: [usize; 3],
    index: [usize; 3],
) -> Vec<[i64; 3]> {
    let reciprocal = cell.reciprocal_columns_inverse_bohr();
    let lattice = cell.lattice_columns_bohr();
    let n: [i64; 3] = std::array::from_fn(|i| mesh[i] as i64);
    let norm2 = |k: [i64; 3]| -> f64 {
        (0..3)
            .map(|i| {
                (0..3)
                    .map(|j| reciprocal[i][j] * (k[j] as f64 / n[j] as f64))
                    .sum::<f64>()
                    .powi(2)
            })
            .sum()
    };
    // Start from the [-1/2, 1/2) wrap; a shorter image q' = w - c has
    // |c_i| = |a_i . (B w - B q')| / 2 pi <= |a_i| 2 |B w| / 2 pi.
    let wrapped: [i64; 3] = std::array::from_fn(|i| {
        let m = (index[i] as i64).rem_euclid(n[i]);
        if 2 * m < n[i] { m } else { m - n[i] }
    });
    let length = norm2(wrapped).sqrt();
    let reach: [i64; 3] = std::array::from_fn(|i| {
        let a = (0..3).map(|c| lattice[c][i].powi(2)).sum::<f64>().sqrt();
        (a * 2. * length / std::f64::consts::TAU * (1. + 1e-12)).floor() as i64
    });
    let mut candidates = Vec::new();
    for c0 in -reach[0]..=reach[0] {
        for c1 in -reach[1]..=reach[1] {
            for c2 in -reach[2]..=reach[2] {
                let c = [c0, c1, c2];
                let k: [i64; 3] = std::array::from_fn(|i| wrapped[i] - c[i] * n[i]);
                candidates.push((norm2(k), k));
            }
        }
    }
    let shortest = candidates.iter().map(|c| c.0).fold(f64::INFINITY, f64::min);
    let mut images: Vec<[i64; 3]> = candidates
        .into_iter()
        .filter(|c| c.0 <= shortest * (1. + IMAGE_TIE))
        .map(|c| c.1)
        .collect();
    images.sort_unstable();
    images
}

/// q units of the G-sphere mode (`sphere_at_gamma`). With one |G| sphere for
/// every q, E_c(q) is a function of the vector q, not of q modulo the
/// reciprocal lattice: q and q - G0 see different sets {q + G}. Convention:
/// the mesh point q contributes the mean over its first-Brillouin-zone
/// images (`first_zone_images`; one image except on zone faces), i.e. one
/// unit (q, image) of weight 1 / (N_q x images). This definition is
/// invariant under every lattice rotation and time reversal (they map the
/// |G| sphere onto itself and the image set of q onto that of S q), so the
/// units are reduced exactly by `group` (mesh-preserving operations with
/// their antiunitary partners, or the identity): an orbit is a set of
/// units whose image vectors are related by S without any G0. Returns
/// (orbit index, representative image, orbit weight, units in the orbit).
fn g_sphere_q_units(
    cell: &dft_geometry::Cell,
    mesh: [usize; 3],
    group: &[SymmetryOperation],
) -> Result<Vec<QEntry>, ScfError> {
    let total = mesh.iter().product::<usize>() as f64;
    let mut units: Vec<([i64; 3], f64)> = Vec::new();
    for i in 0..mesh[0] {
        for j in 0..mesh[1] {
            for l in 0..mesh[2] {
                let images = first_zone_images(cell, mesh, [i, j, l]);
                let w = 1. / (total * images.len() as f64);
                units.extend(images.into_iter().map(|k| (k, w)));
            }
        }
    }
    let lookup: HashMap<[i64; 3], usize> = units
        .iter()
        .enumerate()
        .map(|(u, (k, _))| (*k, u))
        .collect();
    let fractional =
        |k: [i64; 3]| -> [f64; 3] { std::array::from_fn(|i| k[i] as f64 / mesh[i] as f64) };
    let mut assigned = vec![false; units.len()];
    let mut out = Vec::new();
    for u in 0..units.len() {
        if assigned[u] {
            continue;
        }
        let q = fractional(units[u].0);
        let mut members: Vec<usize> = Vec::with_capacity(group.len());
        for op in group {
            let image = op.reciprocal_action(q);
            let mut key = [0_i64; 3];
            for i in 0..3 {
                let x = image[i] * mesh[i] as f64;
                if (x - x.round()).abs() > 1e-6 {
                    return Err(numerical("rpa q units: operation leaves the q mesh"));
                }
                key[i] = x.round() as i64;
            }
            let m = *lookup.get(&key).ok_or_else(|| {
                numerical("rpa q units: rotated first-zone image is not a first-zone image")
            })?;
            if !members.contains(&m) {
                members.push(m);
            }
        }
        if !members.contains(&u) {
            return Err(numerical("rpa q units: the operations lack the identity"));
        }
        let mut weight = 0.;
        members.sort_unstable();
        for &m in &members {
            if assigned[m] {
                return Err(numerical("rpa q units: overlapping orbits"));
            }
            assigned[m] = true;
            weight += units[m].1;
        }
        out.push((out.len(), q, weight, members.len()));
    }
    Ok(out)
}

/// Response sphere at q: the largest cutoff's sphere and, per cutoff, its
/// rows in that sphere and |q+G|.
struct Sphere {
    miller: Vec<[i32; 3]>,
    rows: Vec<Vec<usize>>,
    norms: Vec<Vec<f64>>,
    lo: [i32; 3],
    hi: [i32; 3],
    radius: f64,
}
fn sphere(
    cell: &dft_geometry::Cell,
    cutoffs: &[f64],
    q: [f64; 3],
    at_gamma: bool,
) -> Result<Sphere, ScfError> {
    let max = cutoffs.last().copied().unwrap_or(0.);
    // The Miller set is chosen at q (or at Gamma); |q+G| always uses q.
    let center = if at_gamma { [0.; 3] } else { q };
    let norms_of = |b: &PlaneWaveBasis| -> Vec<f64> {
        let reciprocal = cell.reciprocal_columns_inverse_bohr();
        b.vectors()
            .iter()
            .map(|g| {
                let f: [f64; 3] = std::array::from_fn(|i| f64::from(g.miller[i]) + q[i]);
                (0..3)
                    .map(|i| (0..3).map(|j| reciprocal[i][j] * f[j]).sum::<f64>().powi(2))
                    .sum::<f64>()
                    .sqrt()
            })
            .collect()
    };
    let basis = PlaneWaveBasis::at_k(cell, Hartree(max), center).map_err(numerical)?;
    let miller: Vec<[i32; 3]> = basis.vectors().iter().map(|g| g.miller).collect();
    let index: HashMap<[i32; 3], usize> = miller.iter().enumerate().map(|(i, m)| (*m, i)).collect();
    let mut rows = Vec::with_capacity(cutoffs.len());
    let mut norms = Vec::with_capacity(cutoffs.len());
    for &c in cutoffs {
        let b = PlaneWaveBasis::at_k(cell, Hartree(c), center).map_err(numerical)?;
        let mut r = Vec::with_capacity(b.vectors().len());
        for g in b.vectors() {
            r.push(
                *index
                    .get(&g.miller)
                    .ok_or_else(|| numerical("response spheres are not nested"))?,
            );
        }
        norms.push(norms_of(&b));
        rows.push(r);
    }
    let mut lo = [i32::MAX; 3];
    let mut hi = [i32::MIN; 3];
    for m in &miller {
        for i in 0..3 {
            lo[i] = lo[i].min(m[i]);
            hi[i] = hi[i].max(m[i]);
        }
    }
    let radius = norms_of(&basis).into_iter().fold(0., f64::max);
    Ok(Sphere {
        miller,
        rows,
        norms,
        lo,
        hi,
        radius,
    })
}

/// Copy `rows` of `source` into `target` (reused capacity).
fn gather_rows(
    source: &ComplexMatrix,
    rows: &[usize],
    target: &mut ComplexMatrix,
) -> Result<(), ScfError> {
    let cols = source.cols();
    if target.rows() != rows.len() {
        *target = ComplexMatrix::zeros(rows.len(), cols).map_err(numerical)?;
    } else {
        if cols > target.cols() {
            target.reserve_columns(cols).map_err(numerical)?;
        }
        target.resize_columns(cols).map_err(numerical)?;
    }
    for j in 0..cols {
        let from = source.column(j);
        for (t, &r) in target.column_mut(j).iter_mut().zip(rows) {
            *t = from[r];
        }
    }
    Ok(())
}

/// RPA correlation energy; collective (all arguments identical on every rank).
pub fn rpa_with<K: Collective, S: QStore<K>>(
    input: StartInput<'_>,
    density: &[Array3<f64>],
    options: &RpaOptions,
    comm: &K,
    store: &mut S,
) -> Result<RpaResult, ScfError> {
    let cutoff = input.options.cutoff_hartree;
    let series = &options.ecut_chi_series;
    if series.is_empty()
        || series
            .iter()
            .any(|c| !c.is_finite() || *c <= 0. || *c > 4. * cutoff)
        || series.windows(2).any(|w| w[1] <= w[0])
    {
        return Err(invalid(
            "rpa: ecut_chi series must be ascending, positive and at most 4 x the wavefunction cutoff",
        ));
    }
    if options.frequency_points == 0
        || options.band_block == 0
        || options.k_mesh.contains(&0)
        || !matches!(density.len(), 1 | 2)
    {
        return Err(invalid(
            "rpa: frequency points, band block, k mesh or density channels",
        ));
    }
    let directions = acfdt::lebedev(options.q0_directions).map_err(numerical)?;
    let mut direction_sets: Vec<(String, Directions)> = Vec::new();
    for n in [6, 14, 26, 50] {
        direction_sets.push((
            format!("lebedev-{n}"),
            acfdt::lebedev(n).map_err(numerical)?,
        ));
    }
    direction_sets.push((
        "product-24x48".to_owned(),
        acfdt::product_directions(24).map_err(numerical)?,
    ));
    let cell = input.cell;
    let volume = cell.volume_bohr3();
    let mesh = options.k_mesh;
    let points = gamma_mesh(mesh);
    let nk = points.len();
    let reciprocal = cell.reciprocal_columns_inverse_bohr();
    // q = 0 sphere sizes and the band rule.
    let gamma = sphere(cell, series, [0.; 3], false)?;
    let response_sizes: Vec<usize> = gamma.rows.iter().map(Vec::len).collect();
    let mut basis_sizes = (usize::MAX, 0);
    for p in &points {
        let b = PlaneWaveBasis::at_k(cell, Hartree(cutoff), p.fractional).map_err(numerical)?;
        basis_sizes.0 = basis_sizes.0.min(b.vectors().len());
        basis_sizes.1 = basis_sizes.1.max(b.vectors().len());
    }
    let bands: Vec<usize> = response_sizes
        .iter()
        .map(|&n| options.bands.unwrap_or(n.min(basis_sizes.0)))
        .collect();
    let nscf_bands = *bands.iter().max().expect("non-empty series");
    // Electrons per channel from the density (fixed-moment insulators).
    let grid_points = density[0].len() as f64;
    let channel_electrons: Vec<f64> = density
        .iter()
        .map(|d| d.sum() * volume / grid_points)
        .collect();
    let spin_factor = if density.len() == 1 { 2. } else { 1. };
    let mut occupied = Vec::with_capacity(density.len());
    for n in &channel_electrons {
        let per_band = n / spin_factor;
        if (per_band - per_band.round()).abs() > 1e-6 || per_band.round() < 0. {
            return Err(ScfError::OddElectronCount(*n));
        }
        occupied.push(per_band.round() as usize);
    }
    let max_occupied = *occupied.iter().max().expect("channels");
    if max_occupied == 0 {
        return Err(ScfError::OddElectronCount(0.));
    }
    if bands.iter().any(|&b| b <= max_occupied) {
        return Err(ScfError::InsufficientBands {
            bands: *bands.iter().min().expect("series"),
            occupied: max_occupied,
        });
    }
    let threads = input.options.threads.max(1);
    let (atoms, potentials) = (input.atoms, input.potentials);
    let natoms = atoms.len();
    // q reduction (exact orbits of the mesh).
    let reduction = match &options.q_symmetry {
        QSymmetry::Full => None,
        symmetry => {
            let ops: Vec<SymmetryOperation> = match symmetry {
                QSymmetry::Crystal(ops) => ops
                    .iter()
                    .filter(|o| o.spin_action == SpinAction::Identity && !o.antiunitary)
                    .cloned()
                    .collect(),
                _ => vec![identity(natoms)],
            };
            Some(
                reduce_mesh(&points, &ops, true, false)
                    .map_err(|e| ScfError::Numerical(format!("rpa q reduction: {e}")))?,
            )
        }
    };
    let q_list: Vec<QEntry> = if options.sphere_at_gamma {
        // One |G| sphere for every q: E_c depends on the q vector itself, so
        // q is the first-Brillouin-zone image (convention in `g_sphere_q_units`).
        let group = reduction
            .as_ref()
            .map_or_else(|| vec![identity(natoms)], |m| m.operations.clone());
        g_sphere_q_units(cell, mesh, &group)?
    } else {
        match reduction {
            None => points
                .iter()
                .enumerate()
                .map(|(i, p)| (i, p.fractional, p.weight, 1))
                .collect(),
            Some(map) => {
                let mut multiplicity = vec![0; map.points.len()];
                for e in &map.entries {
                    multiplicity[e.representative_id] += 1;
                }
                map.points
                    .iter()
                    .zip(multiplicity)
                    .enumerate()
                    .map(|(i, (p, m))| (i, p.fractional, p.weight, m))
                    .collect()
            }
        }
    };
    let weight_sum: f64 = q_list.iter().map(|q| q.2).sum();
    if (weight_sum - 1.).abs() > 1e-12 {
        return Err(numerical("q weights do not sum to one"));
    }

    let band_options = BandOptions::new(nscf_bands, options.eigen_tolerance, options.solver);
    let nscf = nscf_bands_with(input, density, &points, band_options, comm, &mut NoStore)?;
    if nscf.channels.len() != density.len() {
        return Err(numerical("NSCF channel count"));
    }
    let (mut vbm, mut cbm) = (f64::NEG_INFINITY, f64::INFINITY);
    let (mut direct, mut largest, mut max_residual) = (f64::INFINITY, 0_f64, 0_f64);
    for (spectra, &occ) in nscf.channels.iter().zip(&occupied) {
        // An empty channel (fully polarized moment) has no transitions.
        if occ == 0 {
            continue;
        }
        let (mut v, mut c) = (f64::NEG_INFINITY, f64::INFINITY);
        for b in spectra {
            let e = &b.eigenvalues;
            v = v.max(e[occ - 1]);
            c = c.min(e[occ]);
            direct = direct.min(e[occ] - e[occ - 1]);
            largest = largest.max(e[nscf_bands - 1] - e[0]);
            max_residual = b.residuals.iter().copied().fold(max_residual, f64::max);
        }
        // A numerical zero for touching bands, not a material threshold.
        if c - v <= 1e-6 {
            return Err(ScfError::NonInsulating {
                gap_hartree: c - v,
                required_hartree: 1e-6,
            });
        }
        vbm = vbm.max(v);
        cbm = cbm.min(c);
    }
    let omega0 = match options.omega0_hartree {
        Some(w) => w,
        None => auto_omega0(direct, largest).map_err(numerical)?,
    };
    let grid = FrequencyGrid::mapped_gauss_legendre(options.frequency_points, omega0)
        .map_err(numerical)?;
    let nw = grid.points.len();
    let dims = nscf.grid_dims;
    let kdata: Vec<KData> = nscf
        .bases
        .iter()
        .map(|b| {
            let miller: Vec<[i32; 3]> = b.vectors().iter().map(|g| g.miller).collect();
            let mut lo = [i32::MAX; 3];
            let mut hi = [i32::MIN; 3];
            for m in &miller {
                for i in 0..3 {
                    lo[i] = lo[i].min(m[i]);
                    hi[i] = hi[i].max(m[i]);
                }
            }
            KData {
                momenta: b.momenta().iter().map(|g| g.cartesian).collect(),
                radius: b
                    .momenta()
                    .iter()
                    .map(|g| g.norm_squared.sqrt())
                    .fold(0., f64::max),
                wrapped: b.k_fractional(),
                miller,
                lo,
                hi,
            }
        })
        .collect();
    // Species for the projector k gradients (q = 0 velocity).
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
    // Band blocks per channel: split at every cutoff's band count.
    let blocks: Vec<Vec<std::ops::Range<usize>>> = occupied
        .iter()
        .map(|&occ| {
            let mut cuts: Vec<usize> = bands.clone();
            cuts.extend((occ..nscf_bands).step_by(options.band_block));
            cuts.push(nscf_bands);
            cuts.retain(|&c| c >= occ);
            cuts.sort_unstable();
            cuts.dedup();
            cuts.windows(2).map(|w| w[0]..w[1]).collect()
        })
        .collect();
    let scale_factor = spin_factor / volume;

    let mut q_points = Vec::with_capacity(q_list.len());
    let mut resumed = 0;
    let mut unit = 0_usize;
    let mut subset = ComplexMatrix::zeros(0, 0).map_err(numerical)?;
    for &(index, fractional, weight, multiplicity) in &q_list {
        let q_mesh = mesh_index(mesh, fractional)?;
        let optical = q_mesh == [0; 3];
        if let Some(record) = store.load(comm, index)? {
            let ok = record.index == index
                && (0..3).all(|i| (record.fractional[i] - fractional[i]).abs() <= 1e-12)
                && (record.weight - weight).abs() <= 1e-14
                && record.optical == optical
                && record.omega0.to_bits() == omega0.to_bits()
                && record.cutoffs.len() == series.len()
                && record.cutoffs.iter().all(|c| c.terms.len() == nw);
            if !ok {
                return Err(numerical(format!(
                    "rpa checkpoint of q index {index} does not match this run"
                )));
            }
            resumed += 1;
            q_points.push(record);
            continue;
        }
        // The evaluated q vector: the stored first-zone image in G-sphere
        // mode (E_c depends on it), the [-1/2, 1/2) wrap in q+G mode (where
        // the sphere follows q and every image gives the same chi0 up to a
        // row permutation).
        let q_wrapped = if options.sphere_at_gamma {
            fractional
        } else {
            dft_geometry::kpoints::wrap_fractional(fractional).map_err(|_| numerical("q wrap"))?
        };
        let sph = sphere(cell, series, q_wrapped, options.sphere_at_gamma)?;
        for m in &sph.miller {
            dft_basis::FftGridSpec::map_miller(*m, dims)
                .map_err(|_| invalid("rpa: the response sphere does not fit the density grid"))?;
        }
        if optical && (sph.rows.iter().any(|r| r[0] != 0) || sph.norms.iter().any(|n| n[0] != 0.)) {
            return Err(numerical("q = 0 sphere must start with G = 0"));
        }
        let mut fft = PairFft::new(dims, &sph.miller).map_err(numerical)?;
        let mut accumulators = sph
            .rows
            .iter()
            .map(|r| Accumulator::new(r.len(), &grid.points, optical, threads))
            .collect::<Result<Vec<_>, _>>()
            .map_err(numerical)?;
        for (channel, spectra) in nscf.channels.iter().enumerate() {
            let occ = occupied[channel];
            if occ == 0 {
                continue;
            }
            let channel_blocks = &blocks[channel];
            for k in 0..nk {
                let km = mesh_index(mesh, points[k].fractional)?;
                let kq = flat(mesh, std::array::from_fn(|i| (km[i] + q_mesh[i]) % mesh[i]));
                let (left, right) = (&kdata[k], &kdata[kq]);
                let mut shift = [0; 3];
                for i in 0..3 {
                    let s = left.wrapped[i] + q_wrapped[i] - right.wrapped[i];
                    if (s - s.round()).abs() > 1e-9 {
                        return Err(numerical("k + q is not on the mesh"));
                    }
                    shift[i] = s.round() as i32;
                }
                if !alias_free(
                    &reciprocal,
                    dims,
                    left,
                    right,
                    shift,
                    (&sph.lo, &sph.hi, sph.radius),
                ) {
                    return Err(invalid(
                        "rpa: the FFT grid cannot hold the pair products at finite q without aliasing (use a larger grid)",
                    ));
                }
                let full_left = gather(comm, &spectra[k], left.miller.len())?;
                let full_right = if kq == k {
                    None
                } else {
                    Some(gather(comm, &spectra[kq], right.miller.len())?)
                };
                let mine: Vec<_> = channel_blocks
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| (unit + i) % comm.size() == comm.rank())
                    .map(|(_, r)| r.clone())
                    .collect();
                unit += channel_blocks.len();
                let local = (|| -> Result<(), ScfError> {
                    if mine.is_empty() {
                        return Ok(());
                    }
                    let lo = KOrbitals {
                        miller: &left.miller,
                        coefficients: &full_left,
                        energies: &spectra[k].eigenvalues,
                    };
                    let ro = KOrbitals {
                        miller: &right.miller,
                        coefficients: full_right.as_ref().unwrap_or(&full_left),
                        energies: &spectra[kq].eigenvalues,
                    };
                    let projectors = if optical {
                        Some(
                            dft_physics::stress::projector_k_gradient(
                                &nscf.bases[k],
                                volume,
                                &species,
                            )
                            .map_err(numerical)?,
                        )
                    } else {
                        None
                    };
                    let mut transitions =
                        Transitions::new(dims, lo, occ, &right.miller, shift).map_err(numerical)?;
                    let scale = scale_factor * points[k].weight;
                    for block in mine {
                        let velocity = projectors.as_ref().map(|p| Velocity {
                            momenta: &left.momenta,
                            projectors: Some(p),
                        });
                        let b = transitions
                            .block(&mut fft, ro, block.clone(), velocity)
                            .map_err(numerical)?;
                        for ((acc, rows), &nb) in accumulators.iter_mut().zip(&sph.rows).zip(&bands)
                        {
                            if block.end > nb {
                                continue;
                            }
                            if rows.len() == sph.miller.len() {
                                acc.add(b.pairs, b.delta, scale, b.slopes)
                                    .map_err(numerical)?;
                            } else {
                                gather_rows(b.pairs, rows, &mut subset)?;
                                acc.add(&subset, b.delta, scale, b.slopes)
                                    .map_err(numerical)?;
                            }
                        }
                    }
                    Ok(())
                })();
                agree(comm, local)?;
            }
        }
        let mut cutoffs = Vec::with_capacity(series.len());
        let mut transitions = 0;
        for (acc, norms) in accumulators.into_iter().zip(&sph.norms) {
            let chi0 = acc.reduce(comm).map_err(numerical)?;
            transitions += chi0.transitions;
            let mut terms = Vec::with_capacity(nw);
            let (mut min_pivot, mut hermiticity) = (f64::INFINITY, 0_f64);
            let mut set_terms = vec![Vec::with_capacity(nw); direction_sets.len()];
            let mut body_terms = Vec::with_capacity(nw);
            for w in 0..nw {
                if optical {
                    let t = acfdt::optical_term(&chi0, w, norms).map_err(numerical)?;
                    terms.push(t.average(&directions).map_err(numerical)?);
                    for ((_, set), out) in direction_sets.iter().zip(&mut set_terms) {
                        out.push(t.average(set).map_err(numerical)?);
                    }
                    body_terms.push(t.body_only());
                    min_pivot = min_pivot.min(t.min_pivot);
                    hermiticity = hermiticity.max(t.hermiticity);
                } else {
                    let t = acfdt::finite_q_term(&chi0.body[w], norms).map_err(numerical)?;
                    terms.push(t.value);
                    min_pivot = min_pivot.min(t.min_pivot);
                    hermiticity = hermiticity.max(t.hermiticity);
                }
            }
            let integral = |t: &[f64]| acfdt::frequency_integral(&grid.weights, t);
            let energy = integral(&terms).map_err(numerical)?;
            let (direction_energies, body_only_energy) = if optical {
                (
                    direction_sets
                        .iter()
                        .zip(&set_terms)
                        .map(|((label, _), t)| Ok((label.clone(), integral(t)?)))
                        .collect::<Result<Vec<_>, dft_response::ResponseError>>()
                        .map_err(numerical)?,
                    Some(integral(&body_terms).map_err(numerical)?),
                )
            } else {
                (Vec::new(), None)
            };
            cutoffs.push(QCutoff {
                response_size: norms.len(),
                terms,
                energy,
                min_pivot,
                hermiticity,
                direction_energies,
                body_only_energy,
            });
        }
        let record = QContribution {
            index,
            fractional,
            weight,
            multiplicity,
            optical,
            transitions,
            omega0,
            cutoffs,
        };
        store.save(comm, &record)?;
        q_points.push(record);
    }

    let mut energies = vec![0.; series.len()];
    let mut frequency_terms = vec![vec![0.; nw]; series.len()];
    let (mut min_pivot, mut max_hermiticity, mut transitions) = (f64::INFINITY, 0_f64, 0);
    for q in &q_points {
        transitions += q.transitions;
        for (c, qc) in q.cutoffs.iter().enumerate() {
            energies[c] += q.weight * qc.energy;
            for (f, t) in frequency_terms[c].iter_mut().zip(&qc.terms) {
                *f += q.weight * t;
            }
            min_pivot = min_pivot.min(qc.min_pivot);
            max_hermiticity = max_hermiticity.max(qc.hermiticity);
        }
    }
    let (extrapolation, extrapolation_inverse_ng) = if series.len() >= 2 {
        // c = N_G^{2/3} makes the fitted abscissa c^{-3/2} = 1 / N_G.
        let shells: Vec<f64> = response_sizes
            .iter()
            .map(|&n| (n as f64).powf(2. / 3.))
            .collect();
        (
            Some(acfdt::extrapolate(series, &energies).map_err(numerical)?),
            acfdt::extrapolate(&shells, &energies).ok(),
        )
    } else {
        (None, None)
    };
    Ok(RpaResult {
        cutoffs: series.clone(),
        response_sizes,
        bands,
        nscf_bands,
        occupied,
        electrons: nscf.electrons,
        atoms: natoms,
        k_points: nk,
        basis_sizes,
        grid_dims: dims,
        frequencies: grid.points.clone(),
        weights: grid.weights.clone(),
        omega0,
        q_reduction: options.q_symmetry.label(),
        q_points,
        resumed,
        energies,
        frequency_terms,
        extrapolation,
        extrapolation_inverse_ng,
        valence_maximum: vbm,
        conduction_minimum: cbm,
        lowest_direct_gap: direct,
        largest_transition: largest,
        max_residual,
        min_pivot,
        max_hermiticity,
        transitions,
    })
}
