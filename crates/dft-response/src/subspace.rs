//! R4a (`docs/roadmap/rpa-crpa-plan.md` 2.8, 2.11; ADR 0059): atomic-orbital
//! projection onto a band window, Loewdin orthonormalization, subspace weights
//! `p_mk` (for R5 cRPA projector exclusion), Wigner-Seitz supercell hoppings
//! `H(R)` and Wannier band interpolation.
//!
//! Pure numerical core: no I/O, serde or logging. Hartree atomic units
//! throughout (energies in Ha, lengths in Bohr, fractional k in `[0, 1)`).
use crate::{ResponseError, input, numerical};
use dft_numerics::linalg::{ComplexMatrix, MatrixOp, gemm, hermitian_eigh};
use num_complex::Complex64;
use std::collections::BTreeMap;
use std::f64::consts::TAU;

/// Default lower bound on the smallest eigenvalue of the projection overlap
/// `S = A^H A`. Overlaps at or below this value indicate a rank-deficient trial
/// orbital set in the chosen band window and are rejected.
pub const DEFAULT_SINGULAR_TOLERANCE: f64 = 1e-10;

/// Projection and Loewdin orthonormalization result at a single k-point.
#[derive(Clone, Debug)]
pub struct KProjection {
    /// 0-based band indices included in the projection window (strictly ascending).
    pub window_bands: Vec<usize>,
    /// Raw overlap matrix `A_mn = <psi_m | g_n>` (`N_win x N_orb`, dimensionless).
    pub overlap_raw: ComplexMatrix,
    /// Loewdin-orthonormalized gauge `U = A S^{-1/2}` (`N_win x N_orb`, semi-unitary).
    pub gauge: ComplexMatrix,
    /// Eigenvalues of `S = A^H A` in ascending order (dimensionless).
    pub overlap_eigenvalues: Vec<f64>,
    /// `max_{i,j} |(U^H U - I)_ij|`.
    pub orthonormality_defect: f64,
    /// `1 - (1 / N_orb) Tr(S)` (fraction of trial orbitals outside the window).
    pub spilling: f64,
    /// Per-band subspace weights `p_m = sum_i |U_mi|^2` for every band at this k
    /// (`0.0` outside `window_bands`). Exact sum rule: `sum_m p_m = N_orb`.
    pub band_weights: Vec<f64>,
    /// Projected Hamiltonian `H^(W)(k) = U^H diag(eps_win) U` (`N_orb x N_orb`, Ha,
    /// explicitly Hermitized).
    pub hamiltonian_k: ComplexMatrix,
    /// Maximum antisymmetric entry dropped when Hermitizing `H^(W)(k)` (Ha).
    pub hermiticity_defect: f64,
}

/// Project `trial` orbitals (`N_pw x N_orb`) onto the `window_bands` columns of
/// `bloch` (`N_pw x N_bands`) and Loewdin-orthonormalize via `S = A^H A`.
pub fn project_k(
    bloch: &ComplexMatrix,
    eigenvalues: &[f64],
    trial: &ComplexMatrix,
    window_bands: &[usize],
    singular_tol: f64,
) -> Result<KProjection, ResponseError> {
    let n_pw = bloch.rows();
    let n_bands = bloch.cols();
    let n_orb = trial.cols();
    let n_win = window_bands.len();
    if n_pw == 0 || n_bands == 0 || n_orb == 0 {
        return Err(input("subspace: empty basis, band set or trial orbitals"));
    }
    if trial.rows() != n_pw {
        return Err(input(format!(
            "subspace: trial row count {} does not match Bloch basis {}",
            trial.rows(),
            n_pw
        )));
    }
    if eigenvalues.len() != n_bands {
        return Err(input(format!(
            "subspace: eigenvalue count {} does not match band columns {}",
            eigenvalues.len(),
            n_bands
        )));
    }
    if !singular_tol.is_finite() || singular_tol <= 0.0 {
        return Err(input(
            "subspace: singular tolerance must be finite and positive",
        ));
    }
    if n_win < n_orb {
        return Err(input(format!(
            "subspace: window has {n_win} bands, fewer than {n_orb} trial orbitals"
        )));
    }
    if eigenvalues.iter().any(|e| !e.is_finite()) {
        return Err(input("subspace: non-finite eigenvalue"));
    }
    for (idx, &b) in window_bands.iter().enumerate() {
        if b >= n_bands {
            return Err(input(format!(
                "subspace: window band index {b} out of range (n_bands = {n_bands})"
            )));
        }
        if idx > 0 && b <= window_bands[idx - 1] {
            return Err(input(
                "subspace: window band indices must be strictly ascending",
            ));
        }
    }

    // Gather window Bloch columns into a contiguous matrix C_win (N_pw x N_win).
    let mut c_win = ComplexMatrix::zeros(n_pw, n_win).map_err(numerical)?;
    for (col, &b) in window_bands.iter().enumerate() {
        c_win.column_mut(col).copy_from_slice(bloch.column(b));
    }

    // A = C_win^H * G  (N_win x N_orb)
    let overlap_raw =
        gemm(&c_win, MatrixOp::Adjoint, trial, MatrixOp::Identity).map_err(numerical)?;

    // S = A^H * A  (N_orb x N_orb), explicitly Hermitized.
    let mut s = gemm(
        &overlap_raw,
        MatrixOp::Adjoint,
        &overlap_raw,
        MatrixOp::Identity,
    )
    .map_err(numerical)?;
    hermitize_in_place(&mut s);

    let trace_s: f64 = (0..n_orb).map(|i| s[(i, i)].re).sum();
    let spilling = 1.0 - trace_s / (n_orb as f64);

    let (overlap_eigenvalues, evecs) = hermitian_eigh(&s).map_err(numerical)?;
    let min_ev = overlap_eigenvalues[0];
    if !min_ev.is_finite() || min_ev <= singular_tol {
        return Err(numerical(format!(
            "projection overlap S = A^H A is rank-deficient (min eigenvalue {min_ev:.3e} <= {singular_tol:.3e})"
        )));
    }

    // S^{-1/2} = V diag(lambda^{-1/2}) V^H
    let mut v_scaled = ComplexMatrix::zeros(n_orb, n_orb).map_err(numerical)?;
    for j in 0..n_orb {
        let factor = overlap_eigenvalues[j].powf(-0.5);
        for i in 0..n_orb {
            v_scaled[(i, j)] = evecs[(i, j)] * factor;
        }
    }
    let inv_sqrt =
        gemm(&v_scaled, MatrixOp::Identity, &evecs, MatrixOp::Adjoint).map_err(numerical)?;

    // U = A * S^{-1/2}  (N_win x N_orb)
    let gauge = gemm(
        &overlap_raw,
        MatrixOp::Identity,
        &inv_sqrt,
        MatrixOp::Identity,
    )
    .map_err(numerical)?;

    // Verify U^H U = I
    let u_gram = gemm(&gauge, MatrixOp::Adjoint, &gauge, MatrixOp::Identity).map_err(numerical)?;
    let mut orthonormality_defect = 0.0_f64;
    for j in 0..n_orb {
        for i in 0..n_orb {
            let target = if i == j {
                Complex64::new(1.0, 0.0)
            } else {
                Complex64::default()
            };
            orthonormality_defect = orthonormality_defect.max((u_gram[(i, j)] - target).norm());
        }
    }

    // Per-band subspace weights p_m = sum_i |U_mi|^2
    let mut band_weights = vec![0.0_f64; n_bands];
    for (row, &b) in window_bands.iter().enumerate() {
        let mut w = 0.0;
        for i in 0..n_orb {
            w += gauge[(row, i)].norm_sqr();
        }
        band_weights[b] = w;
    }

    // H^(W)(k) = U^H diag(eps_win) U
    let mut eu = ComplexMatrix::zeros(n_win, n_orb).map_err(numerical)?;
    for (row, &b) in window_bands.iter().enumerate() {
        let e = eigenvalues[b];
        for j in 0..n_orb {
            eu[(row, j)] = gauge[(row, j)] * e;
        }
    }
    let mut hamiltonian_k =
        gemm(&gauge, MatrixOp::Adjoint, &eu, MatrixOp::Identity).map_err(numerical)?;
    let hermiticity_defect = hermitize_in_place(&mut hamiltonian_k);

    Ok(KProjection {
        window_bands: window_bands.to_vec(),
        overlap_raw,
        gauge,
        overlap_eigenvalues,
        orthonormality_defect,
        spilling,
        band_weights,
        hamiltonian_k,
        hermiticity_defect,
    })
}

/// Raw projection matrix over ALL bands (`N_bands x N_orb`), used for total
/// atomic-orbital populations (`projwfc.x`-comparable) independent of window.
pub fn raw_overlap_all_bands(
    bloch: &ComplexMatrix,
    trial: &ComplexMatrix,
) -> Result<ComplexMatrix, ResponseError> {
    if bloch.rows() != trial.rows() || bloch.rows() == 0 {
        return Err(input("subspace: incompatible Bloch and trial dimensions"));
    }
    gemm(bloch, MatrixOp::Adjoint, trial, MatrixOp::Identity).map_err(numerical)
}

/// Symmetrize a square matrix in place to `(M + M^H)/2`, returning the largest
/// dropped antisymmetric entry norm.
fn hermitize_in_place(m: &mut ComplexMatrix) -> f64 {
    let n = m.rows();
    debug_assert_eq!(n, m.cols());
    let mut defect = 0.0_f64;
    for j in 0..n {
        let d = m[(j, j)];
        defect = defect.max(d.im.abs());
        m[(j, j)] = Complex64::new(d.re, 0.0);
        for i in (j + 1)..n {
            let a = m[(i, j)];
            let b = m[(j, i)];
            defect = defect.max((a - b.conj()).norm() * 0.5);
            let avg = (a + b.conj()) * 0.5;
            m[(i, j)] = avg;
            m[(j, i)] = avg.conj();
        }
    }
    defect
}

/// Projected populations per orbital across a k-point mesh.
/// Returns `(raw_in_window, lowdin_in_window)` for each orbital.
/// `k_weights` should sum to 1 (or any positive normalized weights).
/// `occupations[k]` has length `N_bands` with values in `[0, 1]` (per spin
/// channel); `spin_factor` is `2.0` for spin-unpolarized, `1.0` for collinear.
pub fn projected_populations(
    k_weights: &[f64],
    occupations: &[&[f64]],
    projections: &[KProjection],
    spin_factor: f64,
) -> Result<(Vec<f64>, Vec<f64>), ResponseError> {
    if k_weights.is_empty()
        || k_weights.len() != occupations.len()
        || k_weights.len() != projections.len()
    {
        return Err(input("subspace: mismatched k-point counts for populations"));
    }
    if !spin_factor.is_finite() || spin_factor <= 0.0 {
        return Err(input("subspace: spin factor must be positive"));
    }
    let n_orb = projections[0].gauge.cols();
    let mut raw = vec![0.0_f64; n_orb];
    let mut lowdin = vec![0.0_f64; n_orb];
    for ((&wk, occ), proj) in k_weights.iter().zip(occupations).zip(projections) {
        if proj.gauge.cols() != n_orb {
            return Err(input(
                "subspace: inconsistent orbital count across k-points",
            ));
        }
        for (row, &b) in proj.window_bands.iter().enumerate() {
            let f = *occ
                .get(b)
                .ok_or_else(|| input("subspace: band index exceeds occupation slice"))?;
            let w = spin_factor * wk * f;
            for i in 0..n_orb {
                raw[i] += w * proj.overlap_raw[(row, i)].norm_sqr();
                lowdin[i] += w * proj.gauge[(row, i)].norm_sqr();
            }
        }
    }
    Ok((raw, lowdin))
}

/// One Wigner-Seitz lattice vector `R` of the `N_1 x N_2 x N_3` supercell and
/// its boundary degeneracy `g(R) >= 1`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WsPoint {
    pub r: [i32; 3],
    pub degeneracy: usize,
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn vec_from_int(a: [[f64; 3]; 3], r: [i32; 3]) -> [f64; 3] {
    let (r0, r1, r2) = (r[0] as f64, r[1] as f64, r[2] as f64);
    [
        r0 * a[0][0] + r1 * a[1][0] + r2 * a[2][0],
        r0 * a[0][1] + r1 * a[1][1] + r2 * a[2][1],
        r0 * a[0][2] + r1 * a[1][2] + r2 * a[2][2],
    ]
}

/// Exact Wigner-Seitz lattice vectors `R in Z^3` and degeneracies `g(R)` for
/// the supercell `(N_0 a_0, N_1 a_1, N_2 a_2)`, where `lattice_vectors_bohr[i]`
/// is the `i`-th primitive lattice vector `a_i` (Bohr).
///
/// Search bounds are derived from the longest supercell half-diagonal and the
/// dual (reciprocal) vector norms, so arbitrarily skewed cells are handled
/// without truncation. Verifies `sum_R 1/g(R) = N_0 N_1 N_2`.
pub fn wigner_seitz_points(
    lattice_vectors_bohr: [[f64; 3]; 3],
    mesh: [usize; 3],
) -> Result<Vec<WsPoint>, ResponseError> {
    if mesh.contains(&0) {
        return Err(input("subspace: k-mesh dimensions must be positive"));
    }
    for v in &lattice_vectors_bohr {
        if v.iter().any(|x| !x.is_finite()) {
            return Err(input("subspace: non-finite lattice vector"));
        }
    }
    let a = lattice_vectors_bohr;
    let vol = dot(a[0], cross(a[1], a[2])).abs();
    if !vol.is_finite() || vol <= 1e-12 {
        return Err(input("subspace: degenerate lattice vectors"));
    }
    // Dual vectors b_i / (2 pi) satisfy (b_i / 2pi) . a_j = delta_ij.
    let signed_vol = dot(a[0], cross(a[1], a[2]));
    let dual = [
        cross(a[1], a[2]).map(|x| x / signed_vol),
        cross(a[2], a[0]).map(|x| x / signed_vol),
        cross(a[0], a[1]).map(|x| x / signed_vol),
    ];
    let sc = [
        a[0].map(|x| x * mesh[0] as f64),
        a[1].map(|x| x * mesh[1] as f64),
        a[2].map(|x| x * mesh[2] as f64),
    ];
    // Longest half-diagonal of the supercell parallelepiped bounds |r| for any
    // point in the Wigner-Seitz cell.
    let mut r_max = 0.0_f64;
    for s0 in [-1.0_f64, 1.0] {
        for s1 in [-1.0_f64, 1.0] {
            for s2 in [-1.0_f64, 1.0] {
                let diag = [
                    0.5 * (s0 * sc[0][0] + s1 * sc[1][0] + s2 * sc[2][0]),
                    0.5 * (s0 * sc[0][1] + s1 * sc[1][1] + s2 * sc[2][1]),
                    0.5 * (s0 * sc[0][2] + s1 * sc[1][2] + s2 * sc[2][2]),
                ];
                r_max = r_max.max(norm(diag));
            }
        }
    }
    // Integer search bounds for R (|R_i| <= r_max |dual_i| + 1) and supercell
    // shifts T = (t_0 N_0, t_1 N_1, t_2 N_2) (|t_i| <= 2 r_max |dual_i| / N_i + 1).
    let mut r_bound = [0_i32; 3];
    let mut t_bound = [0_i32; 3];
    for i in 0..3 {
        let d_norm = norm(dual[i]);
        r_bound[i] = (r_max * d_norm).ceil() as i32 + 1;
        t_bound[i] = ((2.0 * r_max * d_norm) / (mesh[i] as f64)).ceil() as i32 + 1;
    }

    let n = [mesh[0] as i32, mesh[1] as i32, mesh[2] as i32];
    let mut shifts: Vec<[i32; 3]> = Vec::new();
    for t0 in -t_bound[0]..=t_bound[0] {
        for t1 in -t_bound[1]..=t_bound[1] {
            for t2 in -t_bound[2]..=t_bound[2] {
                shifts.push([t0 * n[0], t1 * n[1], t2 * n[2]]);
            }
        }
    }

    let mut points = Vec::new();
    for r0 in -r_bound[0]..=r_bound[0] {
        for r1 in -r_bound[1]..=r_bound[1] {
            for r2 in -r_bound[2]..=r_bound[2] {
                let r = [r0, r1, r2];
                let d0 = {
                    let v = vec_from_int(a, r);
                    dot(v, v)
                };
                if d0 > (r_max + 1e-8).powi(2) {
                    continue;
                }
                let tol = 1e-10 * d0.max(1.0);
                let mut min_d = f64::INFINITY;
                let mut count = 0_usize;
                for &t in &shifts {
                    let diff = [r[0] - t[0], r[1] - t[1], r[2] - t[2]];
                    let v = vec_from_int(a, diff);
                    let d = dot(v, v);
                    if d + tol < min_d {
                        min_d = d;
                        count = 1;
                    } else if (d - min_d).abs() <= tol {
                        count += 1;
                    }
                }
                if (d0 - min_d).abs() <= tol {
                    points.push(WsPoint {
                        r,
                        degeneracy: count,
                    });
                }
            }
        }
    }

    let sum_inv_g: f64 = points.iter().map(|p| 1.0 / (p.degeneracy as f64)).sum();
    let expected = (mesh[0] * mesh[1] * mesh[2]) as f64;
    if (sum_inv_g - expected).abs() > 1e-10 * expected {
        return Err(numerical(format!(
            "Wigner-Seitz partition of unity failed: sum(1/g) = {sum_inv_g:.12}, expected {expected}"
        )));
    }
    Ok(points)
}

/// Real-space tight-binding Hamiltonian `H(R)` in the projected Wannier basis
/// on the Wigner-Seitz supercell of a Gamma-centered uniform k-mesh (Hartree).
#[derive(Clone, Debug)]
pub struct RealSpaceHamiltonian {
    pub num_wann: usize,
    pub mesh: [usize; 3],
    pub points: Vec<WsPoint>,
    /// `hoppings[p]` corresponds to `points[p]`, each `num_wann x num_wann` (Ha).
    pub hoppings: Vec<ComplexMatrix>,
}

/// Fourier transform `H^(W)(k)` from a full Gamma-centered uniform k-mesh into
/// the Wigner-Seitz real-space representation `H(R)`.
pub fn build_real_space_hamiltonian(
    lattice_vectors_bohr: [[f64; 3]; 3],
    mesh: [usize; 3],
    k_fractional: &[[f64; 3]],
    hamiltonians_k: &[ComplexMatrix],
) -> Result<RealSpaceHamiltonian, ResponseError> {
    let nk = mesh[0]
        .checked_mul(mesh[1])
        .and_then(|v| v.checked_mul(mesh[2]))
        .ok_or_else(|| input("subspace: k-mesh overflow"))?;
    if nk == 0 || k_fractional.len() != nk || hamiltonians_k.len() != nk {
        return Err(input(format!(
            "subspace: expected {nk} k-points for mesh {:?}, got {} coords and {} Hamiltonians",
            mesh,
            k_fractional.len(),
            hamiltonians_k.len()
        )));
    }
    let num_wann = hamiltonians_k[0].rows();
    if num_wann == 0 {
        return Err(input("subspace: empty Wannier Hamiltonian"));
    }
    for h in hamiltonians_k {
        if h.rows() != num_wann || h.cols() != num_wann {
            return Err(input(
                "subspace: inconsistent Wannier Hamiltonian dimensions",
            ));
        }
    }

    let points = wigner_seitz_points(lattice_vectors_bohr, mesh)?;
    let inv_nk = 1.0 / (nk as f64);
    let mut hoppings = Vec::with_capacity(points.len());
    for pt in &points {
        let mut h_r = ComplexMatrix::zeros(num_wann, num_wann).map_err(numerical)?;
        for (kf, hk) in k_fractional.iter().zip(hamiltonians_k) {
            let phase_arg =
                -TAU * (kf[0] * pt.r[0] as f64 + kf[1] * pt.r[1] as f64 + kf[2] * pt.r[2] as f64);
            let factor = Complex64::from_polar(inv_nk, phase_arg);
            for j in 0..num_wann {
                for i in 0..num_wann {
                    h_r[(i, j)] += hk[(i, j)] * factor;
                }
            }
        }
        hoppings.push(h_r);
    }

    Ok(RealSpaceHamiltonian {
        num_wann,
        mesh,
        points,
        hoppings,
    })
}

impl RealSpaceHamiltonian {
    /// Wannier-interpolated Hamiltonian `H^(W)(k)` (Ha, Hermitized) at an
    /// arbitrary fractional wavevector `k_fractional`.
    pub fn interpolate_matrix(
        &self,
        k_fractional: [f64; 3],
    ) -> Result<ComplexMatrix, ResponseError> {
        if k_fractional.iter().any(|x| !x.is_finite()) {
            return Err(input("subspace: non-finite k-point for interpolation"));
        }
        let n = self.num_wann;
        let mut h_k = ComplexMatrix::zeros(n, n).map_err(numerical)?;
        for (pt, h_r) in self.points.iter().zip(&self.hoppings) {
            let phase_arg = TAU
                * (k_fractional[0] * pt.r[0] as f64
                    + k_fractional[1] * pt.r[1] as f64
                    + k_fractional[2] * pt.r[2] as f64);
            let weight = Complex64::from_polar(1.0 / (pt.degeneracy as f64), phase_arg);
            for j in 0..n {
                for i in 0..n {
                    h_k[(i, j)] += h_r[(i, j)] * weight;
                }
            }
        }
        hermitize_in_place(&mut h_k);
        Ok(h_k)
    }

    /// Wannier-interpolated eigenvalues (Ha, ascending) at `k_fractional`.
    pub fn interpolate_eigenvalues(
        &self,
        k_fractional: [f64; 3],
    ) -> Result<Vec<f64>, ResponseError> {
        let h_k = self.interpolate_matrix(k_fractional)?;
        let (evals, _) = hermitian_eigh(&h_k).map_err(numerical)?;
        Ok(evals)
    }

    /// Maximum Hermiticity symmetry defect `max_{R,i,j} |H_ij(R) - H_ji^*(-R)|` (Ha).
    pub fn hermiticity_defect(&self) -> f64 {
        let index: BTreeMap<[i32; 3], usize> = self
            .points
            .iter()
            .enumerate()
            .map(|(idx, p)| (p.r, idx))
            .collect();
        let n = self.num_wann;
        let mut defect = 0.0_f64;
        for (idx, pt) in self.points.iter().enumerate() {
            let neg = [-pt.r[0], -pt.r[1], -pt.r[2]];
            let Some(&neg_idx) = index.get(&neg) else {
                return f64::INFINITY;
            };
            let h_pos = &self.hoppings[idx];
            let h_neg = &self.hoppings[neg_idx];
            for j in 0..n {
                for i in 0..n {
                    defect = defect.max((h_pos[(i, j)] - h_neg[(j, i)].conj()).norm());
                }
            }
        }
        defect
    }
}
