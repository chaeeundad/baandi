//! R4a (`docs/roadmap/rpa-crpa-plan.md` 2.8, 3.1, 3.6; ADR 0059): atomic-orbital
//! projection onto a band window, Loewdin orthonormalization, subspace weights
//! `p_mk`, projected populations and Wigner-Seitz real-space hoppings `H(R)`.
//!
//! Pipeline: fixed-density NSCF on a full Gamma-centered `N_1 x N_2 x N_3` k
//! mesh (`refscf::nscf`, every band residual-certified), atomic trial orbitals
//! built from UPF `PSWFC` radial wavefunctions with the same Fourier/phase
//! conventions as `refstart`, projection and Loewdin orthonormalization via
//! `dft_response::subspace::project_k`, and Fourier transform to `H(R)`.
use super::{gather, invalid, numerical};
use crate::{
    model::ScfError,
    refscf::nscf::{BandOptions, BandSolver, NoStore, nscf_bands_with},
    refstart::StartInput,
};
use dft_basis::PlaneWaveBasis;
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use dft_mpi::Collective;
use dft_numerics::linalg::ComplexMatrix;
use dft_physics::data::{Array3, Complex64};
pub use dft_response::subspace::{DEFAULT_SINGULAR_TOLERANCE, RealSpaceHamiltonian};
use dft_response::subspace::{
    KProjection, build_real_space_hamiltonian, project_k, projected_populations,
    raw_overlap_all_bands,
};
use dft_upf::{
    AtomicWavefunction, Pseudopotential,
    radial::{RadialEvenRule, projector_transform_with_rule, real_spherical_harmonic},
};
use std::collections::BTreeMap;
use std::f64::consts::PI;

/// Band window selection for subspace projection.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BandWindow {
    /// Energy window `[min_ha, max_ha]` relative to the reference energy
    /// (insulator VBM across the mesh, Ha).
    EnergyRelative { min_ha: f64, max_ha: f64 },
    /// Explicit 0-based inclusive band index range `first..=last`.
    BandRange { first: usize, last: usize },
    /// All solved bands.
    All,
}

/// One atomic-orbital projection selector (expanded into concrete [`TrialOrbital`]s).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ProjectionSelector {
    /// Species label as declared in `[system.atoms]` / `[pseudopotentials]`.
    pub species: String,
    /// Optional 0-based indices into `atoms` (must match `species`); `None` selects
    /// every atom of `species` in input order.
    pub atoms: Option<Vec<usize>>,
    /// Optional angular momentum filter `l in 0..=3`. When `orbitals` is empty,
    /// expands all `2l + 1` real harmonics in Wannier90 standard order.
    pub l: Option<usize>,
    /// Optional UPF `PP_CHI` channel label (e.g. `"3D"`, `"2S"`) to disambiguate
    /// multiple radial shells of the same `l`.
    pub channel: Option<String>,
    /// Orbital names or shorthands (e.g. `["dxy", "dyz", "dxz"]`, `["t2g"]`,
    /// `["sp3"]`, `["s", "p"]`). Empty requires `l` to be set.
    pub orbitals: Vec<String>,
}

/// One term `(channel_index, l, m, coefficient)` of a normalized trial orbital
/// on a single atom.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OrbitalTerm {
    /// Index into `Pseudopotential::atomic_wavefunctions()`.
    pub channel_index: usize,
    pub l: usize,
    pub m: isize,
    pub coeff: f64,
}

/// One resolved trial orbital in the subspace basis.
#[derive(Clone, Debug)]
pub struct TrialOrbital {
    pub index: usize,
    pub atom_index: usize,
    pub species: String,
    pub position_fractional: [f64; 3],
    pub position_cartesian_bohr: [f64; 3],
    pub orbital_name: String,
    pub channel_labels: Vec<String>,
    pub l_values: Vec<usize>,
    pub terms: Vec<OrbitalTerm>,
}

/// Settings for [`subspace_with`] (Hartree atomic units).
#[derive(Clone, Debug)]
pub struct SubspaceOptions {
    pub k_mesh: [usize; 3],
    /// NSCF bands per k; `None` defaults to `max(occupied + 8, occupied + n_orb)`
    /// clamped to the smallest k basis size.
    pub bands: Option<usize>,
    pub solver: BandSolver,
    pub eigen_tolerance: f64,
    pub window: BandWindow,
    pub projections: Vec<ProjectionSelector>,
    pub singular_tolerance: f64,
}

impl SubspaceOptions {
    pub fn new(k_mesh: [usize; 3], projections: Vec<ProjectionSelector>) -> Self {
        Self {
            k_mesh,
            bands: None,
            solver: BandSolver::Dense,
            eigen_tolerance: 1e-6,
            window: BandWindow::All,
            projections,
            singular_tolerance: DEFAULT_SINGULAR_TOLERANCE,
        }
    }
}

/// Complete R4a subspace projection result.
#[derive(Clone, Debug)]
pub struct SubspaceResult {
    pub k_mesh: [usize; 3],
    pub k_points: Vec<KPoint>,
    pub bands: usize,
    pub occupied: usize,
    pub electrons: f64,
    pub reference_energy_ha: f64,
    pub valence_maximum_ha: f64,
    pub conduction_minimum_ha: f64,
    pub max_residual_ha: f64,
    pub basis_sizes: (usize, usize),
    pub orbitals: Vec<TrialOrbital>,
    /// NSCF eigenvalues per k (Ha).
    pub eigenvalues: Vec<Vec<f64>>,
    /// Per-k projection and Loewdin orthonormalization records.
    pub projections: Vec<KProjection>,
    /// Per-orbital raw population summed over all solved occupied bands.
    pub raw_populations_all: Vec<f64>,
    /// Per-orbital raw population summed over occupied bands in the window.
    pub raw_populations_window: Vec<f64>,
    /// Per-orbital Loewdin-orthonormalized population in the window.
    pub lowdin_populations: Vec<f64>,
    /// Worst `orthonormality_defect` across all k-points.
    pub max_orthonormality_defect: f64,
    /// Minimum overlap eigenvalue across all k-points.
    pub min_overlap_eigenvalue: f64,
    /// Mean spilling parameter across all k-points.
    pub mean_spilling: f64,
    /// Real-space tight-binding Hamiltonian `H(R)` on the Wigner-Seitz supercell (Ha).
    pub real_space: RealSpaceHamiltonian,
    /// Maximum `|H_ij(R) - H_ji^*(-R)|` (Ha).
    pub hr_hermiticity_defect: f64,
}

/// Full Gamma-centered uniform `n[0] x n[1] x n[2]` k-mesh with equal weights
/// summing to 1.
pub fn gamma_mesh(n: [usize; 3]) -> Result<Vec<KPoint>, ScfError> {
    if n.contains(&0) {
        return Err(invalid("subspace: k_mesh entries must be positive"));
    }
    let total = (n[0] * n[1] * n[2]) as f64;
    let mut points = Vec::with_capacity(n[0] * n[1] * n[2]);
    for i in 0..n[0] {
        for j in 0..n[1] {
            for k in 0..n[2] {
                points.push(KPoint {
                    fractional: [
                        i as f64 / n[0] as f64,
                        j as f64 / n[1] as f64,
                        k as f64 / n[2] as f64,
                    ],
                    weight: 1.0 / total,
                });
            }
        }
    }
    Ok(points)
}

type ExpandedOrbital = (String, Vec<(usize, isize, f64)>);

/// Parse one orbital name or shorthand into `(canonical_name, Vec<(l, m, coeff)>)`.
/// Follows `dft_upf::radial::real_spherical_harmonic` conventions (ADR 0059).
fn expand_orbital_token(token: &str) -> Result<Vec<ExpandedOrbital>, ScfError> {
    let norm = token.trim().to_ascii_lowercase();
    let single = |name: &str, l: usize, m: isize| vec![(name.to_owned(), vec![(l, m, 1.0)])];
    let out = match norm.as_str() {
        "s" => single("s", 0, 0),
        "pz" => single("pz", 1, 0),
        "px" => single("px", 1, 1),
        "py" => single("py", 1, -1),
        "p" => vec![
            ("pz".into(), vec![(1, 0, 1.0)]),
            ("px".into(), vec![(1, 1, 1.0)]),
            ("py".into(), vec![(1, -1, 1.0)]),
        ],
        "dz2" => single("dz2", 2, 0),
        "dxz" => single("dxz", 2, 1),
        "dyz" => single("dyz", 2, -1),
        "dx2-y2" => single("dx2-y2", 2, 2),
        "dxy" => single("dxy", 2, -2),
        "d" => vec![
            ("dz2".into(), vec![(2, 0, 1.0)]),
            ("dxz".into(), vec![(2, 1, 1.0)]),
            ("dyz".into(), vec![(2, -1, 1.0)]),
            ("dx2-y2".into(), vec![(2, 2, 1.0)]),
            ("dxy".into(), vec![(2, -2, 1.0)]),
        ],
        "t2g" => vec![
            ("dxy".into(), vec![(2, -2, 1.0)]),
            ("dyz".into(), vec![(2, -1, 1.0)]),
            ("dxz".into(), vec![(2, 1, 1.0)]),
        ],
        "eg" => vec![
            ("dz2".into(), vec![(2, 0, 1.0)]),
            ("dx2-y2".into(), vec![(2, 2, 1.0)]),
        ],
        "f0" | "fz3" => single("f0", 3, 0),
        "f+1" | "f1" | "fxz2" => single("f+1", 3, 1),
        "f-1" | "fyz2" => single("f-1", 3, -1),
        "f+2" | "f2" | "fz(x2-y2)" => single("f+2", 3, 2),
        "f-2" | "fxyz" => single("f-2", 3, -2),
        "f+3" | "f3" | "fx(x2-3y2)" => single("f+3", 3, 3),
        "f-3" | "fy(3x2-y2)" => single("f-3", 3, -3),
        "f" => vec![
            ("f0".into(), vec![(3, 0, 1.0)]),
            ("f+1".into(), vec![(3, 1, 1.0)]),
            ("f-1".into(), vec![(3, -1, 1.0)]),
            ("f+2".into(), vec![(3, 2, 1.0)]),
            ("f-2".into(), vec![(3, -2, 1.0)]),
            ("f+3".into(), vec![(3, 3, 1.0)]),
            ("f-3".into(), vec![(3, -3, 1.0)]),
        ],
        // Standard Wannier90 tetrahedral sp3 hybrids (normalized to 1).
        "sp3" => vec![
            (
                "sp3-1".into(),
                vec![(0, 0, 0.5), (1, 1, 0.5), (1, -1, 0.5), (1, 0, 0.5)],
            ),
            (
                "sp3-2".into(),
                vec![(0, 0, 0.5), (1, 1, 0.5), (1, -1, -0.5), (1, 0, -0.5)],
            ),
            (
                "sp3-3".into(),
                vec![(0, 0, 0.5), (1, 1, -0.5), (1, -1, 0.5), (1, 0, -0.5)],
            ),
            (
                "sp3-4".into(),
                vec![(0, 0, 0.5), (1, 1, -0.5), (1, -1, -0.5), (1, 0, 0.5)],
            ),
        ],
        _ => {
            return Err(ScfError::Numerical(format!(
                "subspace: unknown orbital specification '{token}'"
            )));
        }
    };
    Ok(out)
}

fn default_shorthand_for_l(l: usize) -> Result<&'static str, ScfError> {
    match l {
        0 => Ok("s"),
        1 => Ok("p"),
        2 => Ok("d"),
        3 => Ok("f"),
        _ => Err(invalid("subspace: angular momentum l must be in 0..=3")),
    }
}

/// Select the UPF `PP_CHI` channel for angular momentum `l` (and optional label).
fn select_channel<'a>(
    species: &str,
    wfcs: &'a [AtomicWavefunction],
    l: usize,
    channel: Option<&str>,
) -> Result<(usize, &'a AtomicWavefunction), ScfError> {
    let matches: Vec<(usize, &AtomicWavefunction)> =
        wfcs.iter().enumerate().filter(|(_, w)| w.l == l).collect();
    if matches.is_empty() {
        return Err(ScfError::Numerical(format!(
            "subspace: pseudopotential for species '{species}' has no PSWFC channel with l = {l}"
        )));
    }
    if let Some(wanted) = channel {
        let wanted_norm = wanted.trim().to_ascii_lowercase();
        let exact: Vec<_> = matches
            .iter()
            .copied()
            .filter(|(_, w)| {
                w.label
                    .as_deref()
                    .is_some_and(|lbl| lbl.trim().to_ascii_lowercase() == wanted_norm)
            })
            .collect();
        if exact.len() == 1 {
            return Ok(exact[0]);
        }
        return Err(ScfError::Numerical(format!(
            "subspace: species '{species}' l = {l} channel '{wanted}' matched {} PSWFC entries",
            exact.len()
        )));
    }
    // Without an explicit channel label, choose the valence channel: highest
    // occupation (breaking ties by taking the last matching entry).
    let best = matches
        .into_iter()
        .max_by(|(idx_a, a), (idx_b, b)| {
            let oa = a.occupation.unwrap_or(0.0);
            let ob = b.occupation.unwrap_or(0.0);
            oa.total_cmp(&ob).then_with(|| idx_a.cmp(idx_b))
        })
        .expect("non-empty matches");
    Ok(best)
}

/// Resolve all [`ProjectionSelector`]s against `atoms` and `potentials` into a
/// deterministic ordered list of [`TrialOrbital`]s.
pub fn resolve_trial_orbitals(
    cell: &Cell,
    atoms: &[Atom],
    potentials: &[crate::SpeciesPotential<'_>],
    selectors: &[ProjectionSelector],
) -> Result<Vec<TrialOrbital>, ScfError> {
    if selectors.is_empty() {
        return Err(invalid(
            "subspace: at least one projection selector is required",
        ));
    }
    let pp_by_species: BTreeMap<&str, &Pseudopotential> =
        potentials.iter().map(|p| (p.label, p.potential)).collect();

    let mut orbitals = Vec::new();
    for sel in selectors {
        let pp = pp_by_species
            .get(sel.species.as_str())
            .copied()
            .ok_or_else(|| {
                ScfError::Numerical(format!(
                    "subspace: unknown species '{}' in projection selector",
                    sel.species
                ))
            })?;
        let wfcs = pp.atomic_wavefunctions().ok_or_else(|| {
            ScfError::Numerical(format!(
                "subspace: pseudopotential for species '{}' has no PP_PSWFC atomic orbitals",
                sel.species
            ))
        })?;

        let atom_indices: Vec<usize> = match &sel.atoms {
            Some(list) => {
                if list.is_empty() {
                    return Err(invalid("subspace: selector atoms list cannot be empty"));
                }
                for &idx in list {
                    let atom = atoms.get(idx).ok_or_else(|| {
                        ScfError::Numerical(format!(
                            "subspace: atom index {idx} out of range (atoms = {})",
                            atoms.len()
                        ))
                    })?;
                    if atom.species != sel.species {
                        return Err(ScfError::Numerical(format!(
                            "subspace: atom {idx} has species '{}', expected '{}'",
                            atom.species, sel.species
                        )));
                    }
                }
                list.clone()
            }
            None => atoms
                .iter()
                .enumerate()
                .filter(|(_, a)| a.species == sel.species)
                .map(|(i, _)| i)
                .collect(),
        };
        if atom_indices.is_empty() {
            return Err(ScfError::Numerical(format!(
                "subspace: no atoms found for species '{}'",
                sel.species
            )));
        }

        let tokens: Vec<String> = if sel.orbitals.is_empty() {
            let l = sel.l.ok_or_else(|| {
                invalid("subspace: projection selector needs either 'l' or 'orbitals'")
            })?;
            vec![default_shorthand_for_l(l)?.to_owned()]
        } else {
            sel.orbitals.clone()
        };

        let mut expanded = Vec::new();
        for tok in &tokens {
            expanded.extend(expand_orbital_token(tok)?);
        }

        // Validate optional `l` filter against expanded orbital terms.
        if let Some(expected_l) = sel.l {
            for (name, terms) in &expanded {
                if terms.iter().any(|(l, _, _)| *l != expected_l) {
                    return Err(ScfError::Numerical(format!(
                        "subspace: orbital '{name}' does not match selector l = {expected_l}"
                    )));
                }
            }
        }

        for atom_index in atom_indices {
            let pos_frac = atoms[atom_index].fractional_position;
            let pos_cart = cell.fractional_to_cartesian_bohr(pos_frac);
            for (name, raw_terms) in &expanded {
                let mut terms = Vec::with_capacity(raw_terms.len());
                let mut channel_labels = Vec::new();
                let mut l_values = Vec::new();
                for &(l, m, coeff) in raw_terms {
                    let (ch_idx, wfc) =
                        select_channel(&sel.species, wfcs, l, sel.channel.as_deref())?;
                    terms.push(OrbitalTerm {
                        channel_index: ch_idx,
                        l,
                        m,
                        coeff,
                    });
                    let lbl = wfc
                        .label
                        .clone()
                        .unwrap_or_else(|| format!("chi{}_l{l}", ch_idx + 1));
                    if !channel_labels.contains(&lbl) {
                        channel_labels.push(lbl);
                    }
                    if !l_values.contains(&l) {
                        l_values.push(l);
                    }
                }
                orbitals.push(TrialOrbital {
                    index: orbitals.len(),
                    atom_index,
                    species: sel.species.clone(),
                    position_fractional: pos_frac,
                    position_cartesian_bohr: pos_cart,
                    orbital_name: name.clone(),
                    channel_labels,
                    l_values,
                    terms,
                });
            }
        }
    }
    Ok(orbitals)
}

/// Normalize a UPF radial wavefunction `chi(r) = r R(r)` on the radial mesh so
/// that `int_0^inf |chi(r)|^2 dr = 1`.
fn normalized_radial(pp: &Pseudopotential, channel_index: usize) -> Result<Vec<f64>, ScfError> {
    let wfc = &pp
        .atomic_wavefunctions()
        .ok_or_else(|| numerical("missing PSWFC"))?[channel_index];
    let sq: Vec<f64> = wfc.radial_times_r.iter().map(|&v| v * v).collect();
    let norm_sq = pp.mesh().integrate(&sq).map_err(numerical)?;
    if !norm_sq.is_finite() || norm_sq <= 1e-16 {
        return Err(numerical(format!(
            "atomic wavefunction channel {} has zero or non-finite radial norm",
            channel_index
        )));
    }
    let inv = 1.0 / norm_sq.sqrt();
    Ok(wfc.radial_times_r.iter().map(|&v| v * inv).collect())
}

/// Evaluate the plane-wave matrix `G(k)` (`N_pw x N_orb`) of all resolved trial
/// atomic orbitals at one k-point basis (ADR 0059 Section 1).
/// Shared between R4a (`subspace_with`) and R4b (Wannier90 `.amn` export).
pub fn build_trial_matrix_at_k(
    basis: &PlaneWaveBasis,
    volume_bohr3: f64,
    potentials: &[crate::SpeciesPotential<'_>],
    orbitals: &[TrialOrbital],
    radial_rule: RadialEvenRule,
) -> Result<ComplexMatrix, ScfError> {
    let n_pw = basis.vectors().len();
    let n_orb = orbitals.len();
    let pp_by_species: BTreeMap<&str, &Pseudopotential> =
        potentials.iter().map(|p| (p.label, p.potential)).collect();

    // Pre-normalize each referenced (species, channel_index) radial wavefunction once.
    let mut radials: BTreeMap<(&str, usize), Vec<f64>> = BTreeMap::new();
    for orb in orbitals {
        let pp = pp_by_species[orb.species.as_str()];
        for t in &orb.terms {
            if let std::collections::btree_map::Entry::Vacant(e) =
                radials.entry((orb.species.as_str(), t.channel_index))
            {
                let chi = normalized_radial(pp, t.channel_index)?;
                e.insert(chi);
            }
        }
    }

    // Cache radial transforms by (species, channel_index, |k+G|^2 bits).
    let inv_sqrt_vol = 1.0 / volume_bohr3.sqrt();
    let mut radial_cache: BTreeMap<(&str, usize, u64), f64> = BTreeMap::new();
    let mut out = ComplexMatrix::zeros(n_pw, n_orb).map_err(numerical)?;

    for (col, orb) in orbitals.iter().enumerate() {
        let pp = pp_by_species[orb.species.as_str()];
        let pos = orb.position_cartesian_bohr;
        for (ig, g) in basis.momenta().iter().enumerate() {
            let angle = g.cartesian[0] * pos[0] + g.cartesian[1] * pos[1] + g.cartesian[2] * pos[2];
            let phase = Complex64::from_polar(1.0, -angle);
            let q2_bits = g.norm_squared.to_bits();
            let q = g.norm_squared.sqrt();
            let mut sum = Complex64::default();
            for term in &orb.terms {
                let r_val = if let Some(&v) =
                    radial_cache.get(&(orb.species.as_str(), term.channel_index, q2_bits))
                {
                    v
                } else {
                    let chi = &radials[&(orb.species.as_str(), term.channel_index)];
                    let v = projector_transform_with_rule(pp.mesh(), term.l, chi, q, radial_rule)
                        .map_err(numerical)?
                        * inv_sqrt_vol;
                    radial_cache.insert((orb.species.as_str(), term.channel_index, q2_bits), v);
                    v
                };
                let y = if g.norm_squared == 0.0 {
                    if term.l == 0 {
                        1.0 / (4.0 * PI).sqrt()
                    } else {
                        0.0
                    }
                } else {
                    real_spherical_harmonic(term.l, term.m, g.cartesian).map_err(numerical)?
                };
                let il = Complex64::new(0.0, -1.0).powu(term.l as u32);
                sum += il * (term.coeff * r_val * y);
            }
            out[(ig, col)] = phase * sum;
        }
    }
    Ok(out)
}

/// Select 0-based ascending window band indices at one k-point.
pub fn select_window_bands(
    eigenvalues: &[f64],
    window: BandWindow,
    reference_ha: f64,
) -> Result<Vec<usize>, ScfError> {
    match window {
        BandWindow::All => Ok((0..eigenvalues.len()).collect()),
        BandWindow::BandRange { first, last } => {
            if first > last || last >= eigenvalues.len() {
                return Err(ScfError::Numerical(format!(
                    "subspace: band_range ({}..={}) out of bounds for {} solved bands",
                    first + 1,
                    last + 1,
                    eigenvalues.len()
                )));
            }
            Ok((first..=last).collect())
        }
        BandWindow::EnergyRelative { min_ha, max_ha } => {
            if !min_ha.is_finite() || !max_ha.is_finite() || min_ha >= max_ha {
                return Err(invalid("subspace: window_ha must satisfy min < max"));
            }
            let lo = reference_ha + min_ha;
            let hi = reference_ha + max_ha;
            Ok(eigenvalues
                .iter()
                .enumerate()
                .filter(|(_, e)| **e >= lo && **e <= hi)
                .map(|(i, _)| i)
                .collect())
        }
    }
}

/// [`subspace_with`] on the registered communicator (or serially).
pub fn subspace(
    input: StartInput<'_>,
    density: &[Array3<f64>],
    options: &SubspaceOptions,
) -> Result<SubspaceResult, ScfError> {
    if let Some(comm) = crate::parallel::SharedComm::active() {
        subspace_with(input, density, options, &comm)
    } else {
        subspace_with(input, density, options, &crate::refscf::metric::Serial)
    }
}

/// Run R4a atomic-orbital subspace projection on a fixed density (collective;
/// arguments identical on every rank).
pub fn subspace_with<K: Collective>(
    input: StartInput<'_>,
    density: &[Array3<f64>],
    options: &SubspaceOptions,
    comm: &K,
) -> Result<SubspaceResult, ScfError> {
    if density.len() != 1 {
        return Err(invalid(
            "subspace: spin-polarized projection is deferred to R5b",
        ));
    }
    let points = gamma_mesh(options.k_mesh)?;
    let orbitals = resolve_trial_orbitals(
        input.cell,
        input.atoms,
        input.potentials,
        &options.projections,
    )?;
    let n_orb = orbitals.len();
    let cutoff = input.options.cutoff_hartree;
    let mut basis_sizes = (usize::MAX, 0usize);
    for p in &points {
        let b = PlaneWaveBasis::at_k(input.cell, dft_units::Hartree(cutoff), p.fractional)
            .map_err(numerical)?;
        basis_sizes.0 = basis_sizes.0.min(b.vectors().len());
        basis_sizes.1 = basis_sizes.1.max(b.vectors().len());
    }

    // Estimate occupied bands from valence electrons before solving NSCF so
    // default `bands` always covers occupied states plus trial subspace.
    let z_total: f64 = input
        .atoms
        .iter()
        .map(|a| {
            input
                .potentials
                .iter()
                .find(|p| p.label == a.species)
                .map(|p| p.potential.z_valence())
                .ok_or_else(|| invalid("subspace: missing species potential"))
        })
        .collect::<Result<Vec<_>, _>>()?
        .iter()
        .sum();
    let half_est = (z_total / 2.0).ceil() as usize;
    let default_bands = (half_est + n_orb.max(8)).min(basis_sizes.0);
    let bands = options.bands.unwrap_or(default_bands);
    if bands < n_orb {
        return Err(ScfError::Numerical(format!(
            "subspace: requested bands ({bands}) is smaller than trial orbital count ({n_orb})"
        )));
    }
    if bands > basis_sizes.0 {
        return Err(ScfError::Numerical(format!(
            "subspace: requested bands ({bands}) exceeds smallest k-point basis ({})",
            basis_sizes.0
        )));
    }

    let cell = input.cell;
    let potentials = input.potentials;
    let radial_rule = input.options.radial_even_rule;
    let band_options = BandOptions::new(bands, options.eigen_tolerance, options.solver);
    let nscf = nscf_bands_with(input, density, &points, band_options, comm, &mut NoStore)?;
    let half = nscf.electrons / 2.0;
    if (half - half.round()).abs() > 1e-8 || half.round() < 1.0 {
        return Err(ScfError::OddElectronCount(nscf.electrons));
    }
    let occupied = half.round() as usize;
    if bands < occupied {
        return Err(ScfError::InsufficientBands { bands, occupied });
    }

    let spectra = &nscf.channels[0];
    let mut vbm = f64::NEG_INFINITY;
    let mut cbm = f64::INFINITY;
    let mut max_residual = 0.0_f64;
    for b in spectra {
        vbm = vbm.max(b.eigenvalues[occupied - 1]);
        if bands > occupied {
            cbm = cbm.min(b.eigenvalues[occupied]);
        }
        max_residual = b.residuals.iter().copied().fold(max_residual, f64::max);
    }
    let reference_energy_ha = vbm;

    let volume = cell.volume_bohr3();
    let mut projections = Vec::with_capacity(points.len());
    let mut raw_all_pop = vec![0.0_f64; n_orb];
    let mut hamiltonians_k = Vec::with_capacity(points.len());
    let mut k_fracs = Vec::with_capacity(points.len());
    let mut k_weights = Vec::with_capacity(points.len());
    let mut occupations_owned: Vec<Vec<f64>> = Vec::with_capacity(points.len());
    let mut eigenvalues_out = Vec::with_capacity(points.len());

    for (k_idx, (spectrum, basis)) in spectra.iter().zip(&nscf.bases).enumerate() {
        let n_pw = basis.vectors().len();
        let full = gather(comm, spectrum, n_pw)?;
        let trial = build_trial_matrix_at_k(basis, volume, potentials, &orbitals, radial_rule)?;
        let window_bands =
            select_window_bands(&spectrum.eigenvalues, options.window, reference_energy_ha)?;
        let proj = project_k(
            &full,
            &spectrum.eigenvalues,
            &trial,
            &window_bands,
            options.singular_tolerance,
        )
        .map_err(|e| {
            ScfError::Numerical(format!(
                "subspace at k-point {k_idx} ({:?}): {e}",
                points[k_idx].fractional
            ))
        })?;

        let raw_all = raw_overlap_all_bands(&full, &trial).map_err(numerical)?;
        let wk = points[k_idx].weight;
        for m in 0..occupied {
            for i in 0..n_orb {
                raw_all_pop[i] += 2.0 * wk * raw_all[(m, i)].norm_sqr();
            }
        }

        let mut occ = vec![0.0_f64; bands];
        for f in occ.iter_mut().take(occupied) {
            *f = 1.0;
        }
        occupations_owned.push(occ);
        k_fracs.push(points[k_idx].fractional);
        k_weights.push(wk);
        hamiltonians_k.push(proj.hamiltonian_k.clone());
        eigenvalues_out.push(spectrum.eigenvalues.clone());
        projections.push(proj);
    }

    let occ_slices: Vec<&[f64]> = occupations_owned.iter().map(Vec::as_slice).collect();
    let (raw_populations_window, lowdin_populations) =
        projected_populations(&k_weights, &occ_slices, &projections, 2.0).map_err(numerical)?;

    let max_orthonormality_defect = projections
        .iter()
        .map(|p| p.orthonormality_defect)
        .fold(0.0_f64, f64::max);
    let min_overlap_eigenvalue = projections
        .iter()
        .map(|p| p.overlap_eigenvalues[0])
        .fold(f64::INFINITY, f64::min);
    let mean_spilling =
        projections.iter().map(|p| p.spilling).sum::<f64>() / (projections.len() as f64);

    let cols = cell.lattice_columns_bohr();
    let lattice_vectors = std::array::from_fn(|j| [cols[0][j], cols[1][j], cols[2][j]]);
    let real_space =
        build_real_space_hamiltonian(lattice_vectors, options.k_mesh, &k_fracs, &hamiltonians_k)
            .map_err(numerical)?;
    let hr_hermiticity_defect = real_space.hermiticity_defect();

    Ok(SubspaceResult {
        k_mesh: options.k_mesh,
        k_points: points,
        bands,
        occupied,
        electrons: nscf.electrons,
        reference_energy_ha,
        valence_maximum_ha: vbm,
        conduction_minimum_ha: cbm,
        max_residual_ha: max_residual,
        basis_sizes,
        orbitals,
        eigenvalues: eigenvalues_out,
        projections,
        raw_populations_all: raw_all_pop,
        raw_populations_window,
        lowdin_populations,
        max_orthonormality_defect,
        min_overlap_eigenvalue,
        mean_spilling,
        real_space,
        hr_hermiticity_defect,
    })
}
