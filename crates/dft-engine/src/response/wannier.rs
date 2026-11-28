//! R4b (`docs/roadmap/rpa-crpa-plan.md` 2.8, 3.6; ADR 0060): Wannier90 export
//! (`.win`/`.eig`/`.amn`/`.mmn` data) and gauge import (`_u.mat` / `_u_dis.mat`)
//! engine drivers.
//!
//! Reuses R4a (`super::subspace`) atomic trial orbitals and Wigner-Seitz
//! real-space Hamiltonian construction so projected and Wannier90 subspaces
//! share one numerical representation for R5 (cRPA).
use super::{
    gather, invalid, numerical,
    subspace::{
        ProjectionSelector, TrialOrbital, build_trial_matrix_at_k, gamma_mesh,
        resolve_trial_orbitals,
    },
};
use crate::{
    model::ScfError,
    refscf::nscf::{BandOptions, BandSolver, NoStore, nscf_bands_with},
    refstart::StartInput,
};
use dft_basis::PlaneWaveBasis;
use dft_geometry::{Cell, kpoints::KPoint};
use dft_mpi::Collective;
use dft_numerics::linalg::{ComplexMatrix, MatrixOp, gemm};
use dft_physics::data::Array3;
use dft_response::subspace::{
    DEFAULT_SINGULAR_TOLERANCE, RealSpaceHamiltonian, build_real_space_hamiltonian, project_k,
    raw_overlap_all_bands,
};
pub use dft_response::wannier::{
    BVector, BVectorSet, DEFAULT_GAUGE_TOLERANCE, ImportedKProjection, KNeighbor, SpreadReport,
};
use dft_response::wannier::{
    bloch_overlap, evaluate_spreads, find_bvector_shells, mesh_neighbors, project_with_gauge,
};
use std::collections::HashMap;

fn cell_lattice_vectors_bohr(cell: &Cell) -> [[f64; 3]; 3] {
    let cols = cell.lattice_columns_bohr();
    std::array::from_fn(|j| [cols[0][j], cols[1][j], cols[2][j]])
}

fn cell_reciprocal_vectors_inv_bohr(cell: &Cell) -> [[f64; 3]; 3] {
    let cols = cell.reciprocal_columns_inverse_bohr();
    std::array::from_fn(|j| [cols[0][j], cols[1][j], cols[2][j]])
}

/// Settings for [`wannier_export_with`].
#[derive(Clone, Debug)]
pub struct WannierExportOptions {
    pub k_mesh: [usize; 3],
    pub bands: Option<usize>,
    pub solver: BandSolver,
    pub eigen_tolerance: f64,
    pub projections: Vec<ProjectionSelector>,
    /// Optional explicit neighbor list parsed from `seedname.nnkp`
    /// (`explicit_neighbors[k]` for each mesh k-point `k`). When `None`,
    /// automatic Marzari-Vanderbilt B1 shells are constructed.
    pub explicit_neighbors: Option<Vec<Vec<KNeighbor>>>,
}

impl WannierExportOptions {
    pub fn new(k_mesh: [usize; 3], projections: Vec<ProjectionSelector>) -> Self {
        Self {
            k_mesh,
            bands: None,
            solver: BandSolver::Dense,
            eigen_tolerance: 1e-6,
            projections,
            explicit_neighbors: None,
        }
    }
}

/// Complete Wannier90 export dataset (Hartree atomic units).
#[derive(Clone, Debug)]
pub struct WannierExportResult {
    pub k_mesh: [usize; 3],
    pub k_points: Vec<KPoint>,
    pub bands: usize,
    pub occupied: usize,
    pub electrons: f64,
    pub valence_maximum_ha: f64,
    pub conduction_minimum_ha: f64,
    pub max_residual_ha: f64,
    pub orbitals: Vec<TrialOrbital>,
    /// NSCF eigenvalues per k (Ha).
    pub eigenvalues: Vec<Vec<f64>>,
    /// Raw projection overlap `A_mn(k) = <psi_mk | g_nk>` (`bands x num_wann`).
    pub amn: Vec<ComplexMatrix>,
    /// Standalone B1 shell specification (`None` when `--nnkp` was supplied).
    pub bset: Option<BVectorSet>,
    /// Neighbor table `neighbors[k][ib]` matching `mmn[k][ib]`.
    pub neighbors: Vec<Vec<KNeighbor>>,
    /// Bloch overlap matrices `M_mn^(k,b) = <u_mk | u_{n,k+b}>` (`bands x bands`).
    pub mmn: Vec<Vec<ComplexMatrix>>,
    /// `max |M^(k', -b) - (M^(k, b))^H|` across all reversible neighbor pairs.
    pub max_adjoint_defect: f64,
    /// Initial Marzari-Vanderbilt spreads of the R4a Loewdin projection (when
    /// standalone B1 shells are active and projection over all bands succeeds).
    pub initial_spreads: Option<SpreadReport>,
}

/// Compute Bloch overlaps `M^(k,b)` and adjointness defect for a gathered NSCF
/// state on `points`.
fn compute_mmn(
    points: &[KPoint],
    bases: &[PlaneWaveBasis],
    full_evecs: &[ComplexMatrix],
    neighbors: &[Vec<KNeighbor>],
) -> Result<(Vec<Vec<ComplexMatrix>>, f64), ScfError> {
    let nk = points.len();
    if neighbors.len() != nk {
        return Err(ScfError::Numerical(format!(
            "wannier: neighbor table has {} k-points, expected {nk}",
            neighbors.len()
        )));
    }
    let millers: Vec<Vec<[i32; 3]>> = bases
        .iter()
        .map(|b| b.vectors().iter().map(|g| g.miller).collect())
        .collect();
    let k_bases: Vec<[f64; 3]> = bases.iter().map(|b| b.k_fractional()).collect();

    let mut mmn = Vec::with_capacity(nk);
    for (k, neighs) in neighbors.iter().enumerate() {
        let mut per_b = Vec::with_capacity(neighs.len());
        for nb in neighs {
            if nb.k_index != k || nb.neighbor_k_index >= nk {
                return Err(ScfError::Numerical(format!(
                    "wannier: invalid neighbor entry {:?} at k-point {k}",
                    nb
                )));
            }
            let kn = nb.neighbor_k_index;
            let mut dg = [0_i32; 3];
            for a in 0..3 {
                // b = k_mesh[kn] + g_shift - k_mesh[k]
                // dG_basis = k_basis[k] + b - k_basis[kn]
                let b_frac =
                    points[kn].fractional[a] + (nb.g_shift[a] as f64) - points[k].fractional[a];
                let shift_f = k_bases[k][a] + b_frac - k_bases[kn][a];
                let shift_i = shift_f.round();
                if (shift_f - shift_i).abs() > 1e-6 {
                    return Err(ScfError::Numerical(format!(
                        "wannier: non-integer basis shift {shift_f} between k={k} and kn={kn}"
                    )));
                }
                dg[a] = shift_i as i32;
            }
            let m = bloch_overlap(
                &millers[k],
                &full_evecs[k],
                &millers[kn],
                &full_evecs[kn],
                dg,
            )
            .map_err(numerical)?;
            per_b.push(m);
        }
        mmn.push(per_b);
    }

    // Check Hermitian adjoint symmetry M^(kn, -b) == (M^(k, b))^H for every pair
    // whose reverse entry appears in the neighbor list.
    let mut lookup: HashMap<(usize, usize, [i32; 3]), usize> = HashMap::new();
    for (k, neighs) in neighbors.iter().enumerate() {
        for (ib, nb) in neighs.iter().enumerate() {
            lookup.insert((k, nb.neighbor_k_index, nb.g_shift), ib);
        }
    }
    let mut max_adjoint_defect = 0.0_f64;
    for (k, neighs) in neighbors.iter().enumerate() {
        for (ib, nb) in neighs.iter().enumerate() {
            let rev_g = [-nb.g_shift[0], -nb.g_shift[1], -nb.g_shift[2]];
            if let Some(&ib_rev) = lookup.get(&(nb.neighbor_k_index, k, rev_g)) {
                let fwd = &mmn[k][ib];
                let rev = &mmn[nb.neighbor_k_index][ib_rev];
                for j in 0..fwd.cols() {
                    for i in 0..fwd.rows() {
                        max_adjoint_defect =
                            max_adjoint_defect.max((fwd[(i, j)] - rev[(j, i)].conj()).norm());
                    }
                }
            }
        }
    }
    Ok((mmn, max_adjoint_defect))
}

/// [`wannier_export_with`] on the registered communicator (or serially).
pub fn wannier_export(
    input: StartInput<'_>,
    density: &[Array3<f64>],
    options: &WannierExportOptions,
) -> Result<WannierExportResult, ScfError> {
    if let Some(comm) = crate::parallel::SharedComm::active() {
        wannier_export_with(input, density, options, &comm)
    } else {
        wannier_export_with(input, density, options, &crate::refscf::metric::Serial)
    }
}

/// Compute complete Wannier90 `.win`/`.eig`/`.amn`/`.mmn` data on a fixed
/// density (collective; identical arguments on every rank).
pub fn wannier_export_with<K: Collective>(
    input: StartInput<'_>,
    density: &[Array3<f64>],
    options: &WannierExportOptions,
    comm: &K,
) -> Result<WannierExportResult, ScfError> {
    if density.len() != 1 {
        return Err(invalid(
            "wannier-export: spin-polarized export is deferred to R5b",
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
    let mut min_basis = usize::MAX;
    for p in &points {
        let b = PlaneWaveBasis::at_k(input.cell, dft_units::Hartree(cutoff), p.fractional)
            .map_err(numerical)?;
        min_basis = min_basis.min(b.vectors().len());
    }

    let z_total: f64 = input
        .atoms
        .iter()
        .map(|a| {
            input
                .potentials
                .iter()
                .find(|p| p.label == a.species)
                .map(|p| p.potential.z_valence())
                .ok_or_else(|| invalid("wannier-export: missing species potential"))
        })
        .collect::<Result<Vec<_>, _>>()?
        .iter()
        .sum();
    let half_est = (z_total / 2.0).ceil() as usize;
    let default_bands = (half_est + n_orb.max(8)).min(min_basis);
    let bands = options.bands.unwrap_or(default_bands);
    if bands < n_orb {
        return Err(ScfError::Numerical(format!(
            "wannier-export: bands ({bands}) is smaller than trial orbital count ({n_orb})"
        )));
    }
    if bands > min_basis {
        return Err(ScfError::Numerical(format!(
            "wannier-export: bands ({bands}) exceeds smallest k-point basis ({min_basis})"
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

    let volume = cell.volume_bohr3();
    let mut full_evecs = Vec::with_capacity(points.len());
    let mut amn = Vec::with_capacity(points.len());
    let mut loewdin_gauges: Option<Vec<ComplexMatrix>> = Some(Vec::with_capacity(points.len()));
    let all_window: Vec<usize> = (0..bands).collect();

    for (spectrum, basis) in spectra.iter().zip(&nscf.bases) {
        let full = gather(comm, spectrum, basis.vectors().len())?;
        let trial = build_trial_matrix_at_k(basis, volume, potentials, &orbitals, radial_rule)?;
        let a = raw_overlap_all_bands(&full, &trial).map_err(numerical)?;
        if let Some(ref mut gs) = loewdin_gauges {
            match project_k(
                &full,
                &spectrum.eigenvalues,
                &trial,
                &all_window,
                DEFAULT_SINGULAR_TOLERANCE,
            ) {
                Ok(p) => gs.push(p.gauge),
                Err(_) => loewdin_gauges = None,
            }
        }
        amn.push(a);
        full_evecs.push(full);
    }

    let (bset, neighbors) = match &options.explicit_neighbors {
        Some(list) => (None, list.clone()),
        None => {
            let b = find_bvector_shells(cell_reciprocal_vectors_inv_bohr(cell), options.k_mesh)
                .map_err(numerical)?;
            let n = mesh_neighbors(options.k_mesh, &b);
            (Some(b), n)
        }
    };

    let (mmn, max_adjoint_defect) = compute_mmn(&points, &nscf.bases, &full_evecs, &neighbors)?;
    let initial_spreads = match (&bset, &loewdin_gauges) {
        (Some(b), Some(gs)) => evaluate_spreads(b, &neighbors, &mmn, gs).ok(),
        _ => None,
    };

    Ok(WannierExportResult {
        k_mesh: options.k_mesh,
        k_points: points,
        bands,
        occupied,
        electrons: nscf.electrons,
        valence_maximum_ha: vbm,
        conduction_minimum_ha: cbm,
        max_residual_ha: max_residual,
        orbitals,
        eigenvalues: spectra.iter().map(|b| b.eigenvalues.clone()).collect(),
        amn,
        bset,
        neighbors,
        mmn,
        max_adjoint_defect,
        initial_spreads,
    })
}

/// Settings for [`subspace_from_wannier_with`] (importing Wannier90 `_u.mat` /
/// optional `_u_dis.mat` gauges).
#[derive(Clone, Debug)]
pub struct WannierImportOptions {
    pub k_mesh: [usize; 3],
    pub solver: BandSolver,
    pub eigen_tolerance: f64,
    /// Fractional k-points read from `_u.mat`, checked against `gamma_mesh(k_mesh)`.
    pub k_fractional: Vec<[f64; 3]>,
    /// `u_matrices[k]` (`num_wann x num_wann`) from `<seedname>_u.mat`.
    pub u_matrices: Vec<ComplexMatrix>,
    /// Optional `u_dis_matrices[k]` (`num_bands x num_wann`) from `<seedname>_u_dis.mat`.
    pub u_dis_matrices: Option<Vec<ComplexMatrix>>,
    /// Maximum acceptable semi-unitarity defect `max |U^H U - I|`.
    pub gauge_tolerance: f64,
}

/// Complete imported-gauge subspace result (R4b -> R5 cRPA compatible).
#[derive(Clone, Debug)]
pub struct ImportedSubspaceResult {
    pub k_mesh: [usize; 3],
    pub k_points: Vec<KPoint>,
    pub bands: usize,
    pub num_wann: usize,
    pub occupied: usize,
    pub electrons: f64,
    pub reference_energy_ha: f64,
    pub valence_maximum_ha: f64,
    pub conduction_minimum_ha: f64,
    pub max_residual_ha: f64,
    pub eigenvalues: Vec<Vec<f64>>,
    pub projections: Vec<ImportedKProjection>,
    /// Per-Wannier-orbital occupied populations (`2 sum_k w_k sum_m f_mk |U_mj|^2`).
    pub populations: Vec<f64>,
    pub max_orthonormality_defect: f64,
    pub real_space: RealSpaceHamiltonian,
    pub hr_hermiticity_defect: f64,
    pub spreads: SpreadReport,
    pub max_mmn_adjoint_defect: f64,
}

/// [`subspace_from_wannier_with`] on the registered communicator (or serially).
pub fn subspace_from_wannier(
    input: StartInput<'_>,
    density: &[Array3<f64>],
    options: &WannierImportOptions,
) -> Result<ImportedSubspaceResult, ScfError> {
    if let Some(comm) = crate::parallel::SharedComm::active() {
        subspace_from_wannier_with(input, density, options, &comm)
    } else {
        subspace_from_wannier_with(input, density, options, &crate::refscf::metric::Serial)
    }
}

/// Build subspace weights `p_mk`, populations, Wigner-Seitz `H(R)` and exact
/// Marzari-Vanderbilt spreads from imported Wannier90 `_u.mat` / `_u_dis.mat`
/// gauge matrices (collective; identical arguments on every rank).
pub fn subspace_from_wannier_with<K: Collective>(
    input: StartInput<'_>,
    density: &[Array3<f64>],
    options: &WannierImportOptions,
    comm: &K,
) -> Result<ImportedSubspaceResult, ScfError> {
    if density.len() != 1 {
        return Err(invalid(
            "subspace(wannier90): spin-polarized import is deferred to R5b",
        ));
    }
    let points = gamma_mesh(options.k_mesh)?;
    let nk = points.len();
    if options.k_fractional.len() != nk || options.u_matrices.len() != nk {
        return Err(ScfError::Numerical(format!(
            "subspace(wannier90): expected {nk} k-points for mesh {:?}, got {} coords and {} U matrices",
            options.k_mesh,
            options.k_fractional.len(),
            options.u_matrices.len()
        )));
    }
    // Verify k-point ordering matches the uniform Gamma-centered mesh mod Z^3.
    for (k, (got, exp)) in options.k_fractional.iter().zip(&points).enumerate() {
        for a in 0..3 {
            let d = got[a] - exp.fractional[a];
            if (d - d.round()).abs() > 1e-6 {
                return Err(ScfError::Numerical(format!(
                    "subspace(wannier90): k-point {k} mismatch: imported {:?} vs mesh {:?}",
                    got, exp.fractional
                )));
            }
        }
    }

    let num_wann = options.u_matrices[0].cols();
    if num_wann == 0 {
        return Err(invalid("subspace(wannier90): empty U matrix"));
    }
    // Combine U^(tot)(k) = U_dis(k) * U(k) if disentangled, else U(k).
    let mut u_total = Vec::with_capacity(nk);
    let bands = match &options.u_dis_matrices {
        Some(u_dis) => {
            if u_dis.len() != nk {
                return Err(ScfError::Numerical(format!(
                    "subspace(wannier90): _u_dis.mat has {} k-points, expected {nk}",
                    u_dis.len()
                )));
            }
            let nb = u_dis[0].rows();
            for (k, (ud, u)) in u_dis.iter().zip(&options.u_matrices).enumerate() {
                if ud.rows() != nb
                    || ud.cols() != num_wann
                    || u.rows() != num_wann
                    || u.cols() != num_wann
                {
                    return Err(ScfError::Numerical(format!(
                        "subspace(wannier90): shape mismatch at k={k}: U_dis {}x{}, U {}x{}",
                        ud.rows(),
                        ud.cols(),
                        u.rows(),
                        u.cols()
                    )));
                }
                let prod =
                    gemm(ud, MatrixOp::Identity, u, MatrixOp::Identity).map_err(numerical)?;
                u_total.push(prod);
            }
            nb
        }
        None => {
            for (k, u) in options.u_matrices.iter().enumerate() {
                if u.rows() != num_wann || u.cols() != num_wann {
                    return Err(ScfError::Numerical(format!(
                        "subspace(wannier90): U at k={k} has shape {}x{}, expected {num_wann}x{num_wann}",
                        u.rows(),
                        u.cols()
                    )));
                }
                u_total.push(u.clone());
            }
            num_wann
        }
    };

    let cell = input.cell;
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

    let mut full_evecs = Vec::with_capacity(nk);
    let mut projections = Vec::with_capacity(nk);
    let mut hamiltonians_k = Vec::with_capacity(nk);
    let mut k_fracs = Vec::with_capacity(nk);
    let mut populations = vec![0.0_f64; num_wann];
    let mut gauges_for_spread = Vec::with_capacity(nk);

    for (k, ((spectrum, basis), u_k)) in spectra.iter().zip(&nscf.bases).zip(u_total).enumerate() {
        let full = gather(comm, spectrum, basis.vectors().len())?;
        let proj = project_with_gauge(&spectrum.eigenvalues, u_k, options.gauge_tolerance)
            .map_err(|e| ScfError::Numerical(format!("subspace(wannier90) k={k}: {e}")))?;
        let wk = points[k].weight;
        for m in 0..occupied {
            for (j, pop) in populations.iter_mut().enumerate() {
                *pop += 2.0 * wk * proj.gauge[(m, j)].norm_sqr();
            }
        }
        hamiltonians_k.push(proj.hamiltonian_k.clone());
        k_fracs.push(points[k].fractional);
        gauges_for_spread.push(proj.gauge.clone());
        projections.push(proj);
        full_evecs.push(full);
    }

    let max_orthonormality_defect = projections
        .iter()
        .map(|p| p.orthonormality_defect)
        .fold(0.0_f64, f64::max);

    let real_space = build_real_space_hamiltonian(
        cell_lattice_vectors_bohr(cell),
        options.k_mesh,
        &k_fracs,
        &hamiltonians_k,
    )
    .map_err(numerical)?;
    let hr_hermiticity_defect = real_space.hermiticity_defect();

    let bset = find_bvector_shells(cell_reciprocal_vectors_inv_bohr(cell), options.k_mesh)
        .map_err(numerical)?;
    let neighbors = mesh_neighbors(options.k_mesh, &bset);
    let (mmn, max_mmn_adjoint_defect) = compute_mmn(&points, &nscf.bases, &full_evecs, &neighbors)?;
    let spreads =
        evaluate_spreads(&bset, &neighbors, &mmn, &gauges_for_spread).map_err(numerical)?;

    Ok(ImportedSubspaceResult {
        k_mesh: options.k_mesh,
        k_points: points,
        bands,
        num_wann,
        occupied,
        electrons: nscf.electrons,
        reference_energy_ha: vbm,
        valence_maximum_ha: vbm,
        conduction_minimum_ha: cbm,
        max_residual_ha: max_residual,
        eigenvalues: spectra.iter().map(|b| b.eigenvalues.clone()).collect(),
        projections,
        populations,
        max_orthonormality_defect,
        real_space,
        hr_hermiticity_defect,
        spreads,
        max_mmn_adjoint_defect,
    })
}
