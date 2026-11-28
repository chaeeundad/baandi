//! R4b (`docs/roadmap/rpa-crpa-plan.md` 2.8, 3.6; ADR 0060): Wannier90
//! numerical core — Marzari-Vanderbilt B1 neighbor shells, Bloch-orbital
//! overlaps `M_mn^(k,b) = <u_mk | u_{n,k+b}>`, Marzari-Vanderbilt spread
//! diagnostics (`Omega_I`, `Omega`, centers `r_bar_n`), and imported-gauge
//! subspace projection.
//!
//! Pure numerical core: no I/O, serde or logging. Hartree atomic units
//! throughout (lengths in Bohr, reciprocal vectors in Bohr^-1, energies in Ha,
//! spreads in Bohr^2).
use crate::{ResponseError, input, numerical};
use dft_numerics::linalg::{ComplexMatrix, MatrixOp, gemm};
use num_complex::Complex64;
use std::collections::HashMap;

/// Maximum acceptable residual on the Marzari-Vanderbilt B1 completeness
/// identity `max_{a,c} |sum_b w_b b_a b_c - delta_ac|`.
pub const B1_COMPLETENESS_TOLERANCE: f64 = 1e-10;

/// Maximum semi-unitarity defect `max |U^H U - I|` accepted for an imported
/// Wannier90 gauge matrix.
pub const DEFAULT_GAUGE_TOLERANCE: f64 = 1e-8;

/// One finite-difference neighbor vector `b` on the `N_0 x N_1 x N_2` k-mesh.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BVector {
    /// Integer step in units of the sub-mesh vectors `(b_0/N_0, b_1/N_1, b_2/N_2)`.
    pub mesh_step: [i32; 3],
    /// Fractional displacement in primitive reciprocal coordinates (`mesh_step[i] / N_i`).
    pub fractional: [f64; 3],
    /// Cartesian displacement (Bohr^-1).
    pub cartesian: [f64; 3],
    /// B1 finite-difference weight `w_b` (Bohr^2).
    pub weight: f64,
    /// 0-based shell index.
    pub shell: usize,
}

/// Complete set of neighbor `b`-vectors solving the Marzari-Vanderbilt B1
/// completeness condition `sum_b w_b b_a b_c = delta_ac`.
#[derive(Clone, Debug)]
pub struct BVectorSet {
    pub mesh: [usize; 3],
    pub vectors: Vec<BVector>,
    pub shells: usize,
    /// `max_{a,c} |sum_b w_b b_a b_c - delta_ac|`.
    pub completeness_residual: f64,
}

/// One `(ik, ik_nn, g0)` neighbor entry for `.mmn` / spread evaluation:
/// `k_mesh[ik] + b = k_mesh[ik_nn] + g0`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct KNeighbor {
    pub k_index: usize,
    pub neighbor_k_index: usize,
    pub g_shift: [i32; 3],
}

/// Marzari-Vanderbilt spread decomposition (Bohr^2) and Wannier centers (Bohr).
#[derive(Clone, Debug)]
pub struct SpreadReport {
    /// Gauge-invariant spread `Omega_I` (Bohr^2).
    pub omega_i_bohr2: f64,
    /// Gauge-dependent spread `Omega_tilde = Omega - Omega_I` (Bohr^2).
    pub omega_tilde_bohr2: f64,
    /// Total spread `Omega = Omega_I + Omega_tilde` (Bohr^2).
    pub omega_total_bohr2: f64,
    /// Per-orbital total spread `Omega_n` (Bohr^2).
    pub orbital_spreads_bohr2: Vec<f64>,
    /// Per-orbital Wannier centers `r_bar_n` (Cartesian, Bohr).
    pub centers_bohr: Vec<[f64; 3]>,
}

/// Projection summary built from an imported Wannier90 gauge `U^(tot)(k)`
/// (`N_bands x N_wann`).
#[derive(Clone, Debug)]
pub struct ImportedKProjection {
    pub gauge: ComplexMatrix,
    pub orthonormality_defect: f64,
    pub band_weights: Vec<f64>,
    pub hamiltonian_k: ComplexMatrix,
    pub hermiticity_defect: f64,
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

/// Solve a small real symmetric positive-definite linear system `A x = rhs`
/// (up to 6x6) via Gaussian elimination with partial pivoting.
fn solve_small_linear(a: &[Vec<f64>], rhs: &[f64]) -> Option<Vec<f64>> {
    let n = rhs.len();
    let mut m: Vec<Vec<f64>> = a
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let mut r = row.clone();
            r.push(rhs[i]);
            r
        })
        .collect();
    for col in 0..n {
        let mut pivot = col;
        let mut best = m[col][col].abs();
        for (r, row) in m.iter().enumerate().skip(col + 1) {
            let v = row[col].abs();
            if v > best {
                best = v;
                pivot = r;
            }
        }
        if !best.is_finite() || best <= 1e-13 {
            return None;
        }
        m.swap(col, pivot);
        let diag = m[col][col];
        for r in (col + 1)..n {
            let factor = m[r][col] / diag;
            for c in col..=n {
                m[r][c] -= factor * m[col][c];
            }
        }
    }
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut s = m[i][n];
        for j in (i + 1)..n {
            s -= m[i][j] * x[j];
        }
        x[i] = s / m[i][i];
        if !x[i].is_finite() {
            return None;
        }
    }
    Some(x)
}

/// Solve least squares `min_w ||C w - target||_2` where `C` has 6 rows and
/// `ns <= 6` columns (each column is a shell's `[xx, yy, zz, xy, xz, yz]`),
/// returning `(weights, max_residual)` if full column rank.
fn solve_b1_least_squares(cols: &[[f64; 6]]) -> Option<(Vec<f64>, f64)> {
    let target = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0];
    let ns = cols.len();
    if ns == 0 || ns > 6 {
        return None;
    }
    // Normal equations: (C^T C) w = C^T target.
    let mut ata = vec![vec![0.0_f64; ns]; ns];
    let mut rhs = vec![0.0_f64; ns];
    for j in 0..ns {
        for r in 0..6 {
            rhs[j] += cols[j][r] * target[r];
        }
        for k in 0..ns {
            let mut s = 0.0;
            for r in 0..6 {
                s += cols[j][r] * cols[k][r];
            }
            ata[j][k] = s;
        }
    }
    let w = solve_small_linear(&ata, &rhs)?;
    let mut res = 0.0_f64;
    for (r, &tgt) in target.iter().enumerate() {
        let mut val = 0.0;
        for j in 0..ns {
            val += cols[j][r] * w[j];
        }
        res = res.max((val - tgt).abs());
    }
    Some((w, res))
}

/// Check if `cand` is linearly independent of `basis` (columns in R^6) via
/// modified Gram-Schmidt.
fn is_independent_r6(basis: &[[f64; 6]], cand: [f64; 6]) -> bool {
    let mut v = cand;
    let orig_norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    if !orig_norm.is_finite() || orig_norm <= 1e-14 {
        return false;
    }
    let mut orth: Vec<[f64; 6]> = Vec::with_capacity(basis.len());
    for &b in basis {
        let mut q = b;
        for o in &orth {
            let d: f64 = q.iter().zip(o).map(|(x, y)| x * y).sum();
            for (x, y) in q.iter_mut().zip(o) {
                *x -= d * y;
            }
        }
        let nrm = q.iter().map(|x| x * x).sum::<f64>().sqrt();
        if nrm > 1e-14 {
            for x in &mut q {
                *x /= nrm;
            }
            orth.push(q);
        }
    }
    for o in &orth {
        let d: f64 = v.iter().zip(o).map(|(x, y)| x * y).sum();
        for (x, y) in v.iter_mut().zip(o) {
            *x -= d * y;
        }
    }
    let rem = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    rem > 1e-8 * orig_norm
}

/// Find the minimal set of finite-difference neighbor shells solving the
/// Marzari-Vanderbilt B1 completeness condition `sum_b w_b b_a b_c = delta_ac`.
/// `reciprocal_vectors_inv_bohr[i]` is the `i`-th primitive reciprocal vector
/// `b_i` (Bohr^-1, columns of `B = 2 pi A^{-T}`).
pub fn find_bvector_shells(
    reciprocal_vectors_inv_bohr: [[f64; 3]; 3],
    mesh: [usize; 3],
) -> Result<BVectorSet, ResponseError> {
    if mesh.contains(&0) {
        return Err(input("wannier: k-mesh entries must be positive"));
    }
    for v in &reciprocal_vectors_inv_bohr {
        if v.iter().any(|x| !x.is_finite()) {
            return Err(input("wannier: non-finite reciprocal lattice vector"));
        }
    }
    // Sub-mesh step vectors d_i = b_i / N_i.
    let d = [
        reciprocal_vectors_inv_bohr[0].map(|x| x / mesh[0] as f64),
        reciprocal_vectors_inv_bohr[1].map(|x| x / mesh[1] as f64),
        reciprocal_vectors_inv_bohr[2].map(|x| x / mesh[2] as f64),
    ];
    let vol = dot3(d[0], cross3(d[1], d[2])).abs();
    if !vol.is_finite() || vol <= 1e-14 {
        return Err(input("wannier: degenerate reciprocal lattice"));
    }

    // Search sub-mesh integer lattice vectors in [-search, search]^3 (expanded
    // dynamically for skewed cells up to radius 4).
    let max_d = norm3(d[0]).max(norm3(d[1])).max(norm3(d[2]));
    let dual = [
        cross3(d[1], d[2]).map(|x| x / dot3(d[0], cross3(d[1], d[2]))),
        cross3(d[2], d[0]).map(|x| x / dot3(d[0], cross3(d[1], d[2]))),
        cross3(d[0], d[1]).map(|x| x / dot3(d[0], cross3(d[1], d[2]))),
    ];
    let radius = 2.5 * max_d;
    let bounds: [i32; 3] =
        std::array::from_fn(|i| ((radius * norm3(dual[i])).ceil() as i32).clamp(2, 6));

    let mut cands: Vec<([i32; 3], [f64; 3], f64)> = Vec::new();
    for n0 in -bounds[0]..=bounds[0] {
        for n1 in -bounds[1]..=bounds[1] {
            for n2 in -bounds[2]..=bounds[2] {
                if n0 == 0 && n1 == 0 && n2 == 0 {
                    continue;
                }
                let cart = [
                    n0 as f64 * d[0][0] + n1 as f64 * d[1][0] + n2 as f64 * d[2][0],
                    n0 as f64 * d[0][1] + n1 as f64 * d[1][1] + n2 as f64 * d[2][1],
                    n0 as f64 * d[0][2] + n1 as f64 * d[1][2] + n2 as f64 * d[2][2],
                ];
                let r2 = dot3(cart, cart);
                cands.push(([n0, n1, n2], cart, r2));
            }
        }
    }
    cands.sort_by(|a, b| a.2.total_cmp(&b.2).then_with(|| a.0.cmp(&b.0)));

    // Group into shells of equal |b|^2.
    let mut shells: Vec<Vec<([i32; 3], [f64; 3])>> = Vec::new();
    let mut shell_r2: Vec<f64> = Vec::new();
    for (step, cart, r2) in cands {
        if let Some(&last_r2) = shell_r2.last()
            && (r2 - last_r2).abs() <= 1e-7 * last_r2.max(1e-12)
        {
            shells.last_mut().unwrap().push((step, cart));
        } else {
            shell_r2.push(r2);
            shells.push(vec![(step, cart)]);
        }
    }

    // Greedily pick linearly independent shells until B1 is solved.
    let mut selected_shells: Vec<usize> = Vec::new();
    let mut shell_cols: Vec<[f64; 6]> = Vec::new();
    let mut solution: Option<(Vec<f64>, f64)> = None;

    for (s_idx, members) in shells.iter().enumerate() {
        let mut col = [0.0_f64; 6];
        for &(_, b) in members {
            col[0] += b[0] * b[0];
            col[1] += b[1] * b[1];
            col[2] += b[2] * b[2];
            col[3] += b[0] * b[1];
            col[4] += b[0] * b[2];
            col[5] += b[1] * b[2];
        }
        if !is_independent_r6(&shell_cols, col) {
            continue;
        }
        selected_shells.push(s_idx);
        shell_cols.push(col);
        if let Some((weights, res)) = solve_b1_least_squares(&shell_cols)
            && res <= B1_COMPLETENESS_TOLERANCE
        {
            solution = Some((weights, res));
            break;
        }
        if shell_cols.len() == 6 {
            break;
        }
    }

    let (weights, completeness_residual) = solution.ok_or_else(|| {
        numerical("failed to solve Marzari-Vanderbilt B1 completeness shells for this cell/mesh")
    })?;

    let mut vectors = Vec::new();
    for (local_shell, (&s_idx, &w)) in selected_shells.iter().zip(&weights).enumerate() {
        for &(step, cart) in &shells[s_idx] {
            vectors.push(BVector {
                mesh_step: step,
                fractional: [
                    step[0] as f64 / mesh[0] as f64,
                    step[1] as f64 / mesh[1] as f64,
                    step[2] as f64 / mesh[2] as f64,
                ],
                cartesian: cart,
                weight: w,
                shell: local_shell,
            });
        }
    }

    Ok(BVectorSet {
        mesh,
        vectors,
        shells: selected_shells.len(),
        completeness_residual,
    })
}

/// Build the full neighbor list `neighbors[ik][ib]` for a uniform Gamma-centered
/// mesh ordered as `(i0 * N1 + i1) * N2 + i2`.
pub fn mesh_neighbors(mesh: [usize; 3], bset: &BVectorSet) -> Vec<Vec<KNeighbor>> {
    let (n0, n1, n2) = (mesh[0] as i32, mesh[1] as i32, mesh[2] as i32);
    let nk = mesh[0] * mesh[1] * mesh[2];
    let mut out = Vec::with_capacity(nk);
    for i0 in 0..n0 {
        for i1 in 0..n1 {
            for i2 in 0..n2 {
                let k_index = ((i0 as usize) * mesh[1] + (i1 as usize)) * mesh[2] + (i2 as usize);
                let mut list = Vec::with_capacity(bset.vectors.len());
                for bv in &bset.vectors {
                    let s0 = i0 + bv.mesh_step[0];
                    let s1 = i1 + bv.mesh_step[1];
                    let s2 = i2 + bv.mesh_step[2];
                    let g0 = s0.div_euclid(n0);
                    let g1 = s1.div_euclid(n1);
                    let g2 = s2.div_euclid(n2);
                    let j0 = s0.rem_euclid(n0) as usize;
                    let j1 = s1.rem_euclid(n1) as usize;
                    let j2 = s2.rem_euclid(n2) as usize;
                    let neighbor_k_index = (j0 * mesh[1] + j1) * mesh[2] + j2;
                    list.push(KNeighbor {
                        k_index,
                        neighbor_k_index,
                        g_shift: [g0, g1, g2],
                    });
                }
                out.push(list);
            }
        }
    }
    out
}

/// Exact Bloch overlap matrix `M_mn^(k,b) = <u_mk | u_{n,k+b}>` (`N_bands x N_bands`)
/// computed in reciprocal space as `sum_G c_mk^*(G) c_{n,k'}(G + dG_basis)`,
/// where `dG_basis = round(k_basis + b - k'_basis)` accounts for arbitrary
/// wrapping of `k_basis` / `k'_basis`.
pub fn bloch_overlap(
    left_miller: &[[i32; 3]],
    left_coeffs: &ComplexMatrix,
    right_miller: &[[i32; 3]],
    right_coeffs: &ComplexMatrix,
    dg_basis: [i32; 3],
) -> Result<ComplexMatrix, ResponseError> {
    if left_miller.len() != left_coeffs.rows() || right_miller.len() != right_coeffs.rows() {
        return Err(input(
            "wannier: Miller list length does not match eigenvector rows",
        ));
    }
    let nb_l = left_coeffs.cols();
    let nb_r = right_coeffs.cols();
    if nb_l == 0 || nb_r == 0 {
        return Err(input("wannier: empty band set for Bloch overlap"));
    }
    let right_map: HashMap<[i32; 3], usize> = right_miller
        .iter()
        .enumerate()
        .map(|(idx, &m)| (m, idx))
        .collect();
    let mut pairs: Vec<(usize, usize)> = Vec::with_capacity(left_miller.len());
    for (il, &g) in left_miller.iter().enumerate() {
        let shifted = [g[0] + dg_basis[0], g[1] + dg_basis[1], g[2] + dg_basis[2]];
        if let Some(&ir) = right_map.get(&shifted) {
            pairs.push((il, ir));
        }
    }
    if pairs.is_empty() {
        return ComplexMatrix::zeros(nb_l, nb_r).map_err(numerical);
    }
    let mut l_sub = ComplexMatrix::zeros(pairs.len(), nb_l).map_err(numerical)?;
    let mut r_sub = ComplexMatrix::zeros(pairs.len(), nb_r).map_err(numerical)?;
    for (row, &(il, ir)) in pairs.iter().enumerate() {
        for b in 0..nb_l {
            l_sub[(row, b)] = left_coeffs[(il, b)];
        }
        for b in 0..nb_r {
            r_sub[(row, b)] = right_coeffs[(ir, b)];
        }
    }
    gemm(&l_sub, MatrixOp::Adjoint, &r_sub, MatrixOp::Identity).map_err(numerical)
}

/// Evaluate Marzari-Vanderbilt spreads (`Omega_I`, `Omega_tilde`, `Omega`) and
/// Wannier centers `r_bar_n` from the Bloch overlaps `overlaps[k][ib]`,
/// neighbor table `neighbors[k][ib]`, b-vectors `bset`, and semi-unitary gauge
/// matrices `gauges[k]` (`N_bands x N_wann`).
pub fn evaluate_spreads(
    bset: &BVectorSet,
    neighbors: &[Vec<KNeighbor>],
    overlaps: &[Vec<ComplexMatrix>],
    gauges: &[ComplexMatrix],
) -> Result<SpreadReport, ResponseError> {
    let nk = gauges.len();
    if nk == 0 || neighbors.len() != nk || overlaps.len() != nk {
        return Err(input(
            "wannier: mismatched k-point counts for spread evaluation",
        ));
    }
    let n_wann = gauges[0].cols();
    if n_wann == 0 {
        return Err(input("wannier: empty Wannier gauge"));
    }
    let nb = bset.vectors.len();
    let inv_nk = 1.0 / (nk as f64);

    // Transform M^(W)(k, b) = U(k)^H M(k, b) U(k') once per (k, b).
    let mut mw: Vec<Vec<ComplexMatrix>> = Vec::with_capacity(nk);
    for k in 0..nk {
        if neighbors[k].len() != nb || overlaps[k].len() != nb {
            return Err(input("wannier: inconsistent neighbor count per k-point"));
        }
        let mut row = Vec::with_capacity(nb);
        for ib in 0..nb {
            let kn = neighbors[k][ib].neighbor_k_index;
            let mu = gemm(
                &overlaps[k][ib],
                MatrixOp::Identity,
                &gauges[kn],
                MatrixOp::Identity,
            )
            .map_err(numerical)?;
            let m_w =
                gemm(&gauges[k], MatrixOp::Adjoint, &mu, MatrixOp::Identity).map_err(numerical)?;
            row.push(m_w);
        }
        mw.push(row);
    }

    // Pass 1: Wannier centers r_n = -(1/N_k) sum_{k,b} w_b b Im ln M^(W)_nn(k,b)
    // and gauge-invariant spread Omega_I.
    let mut centers = vec![[0.0_f64; 3]; n_wann];
    let mut omega_i = 0.0_f64;
    for mw_k in &mw {
        for (ib, bv) in bset.vectors.iter().enumerate() {
            let m = &mw_k[ib];
            let mut frob2 = 0.0_f64;
            for j in 0..n_wann {
                for i in 0..n_wann {
                    frob2 += m[(i, j)].norm_sqr();
                }
            }
            omega_i += inv_nk * bv.weight * ((n_wann as f64) - frob2);
            for (n, center) in centers.iter_mut().enumerate() {
                let diag = m[(n, n)];
                let phase = diag.im.atan2(diag.re);
                let factor = -inv_nk * bv.weight * phase;
                for a in 0..3 {
                    center[a] += factor * bv.cartesian[a];
                }
            }
        }
    }

    // Pass 2: total spread per orbital.
    let mut orbital_spreads = vec![0.0_f64; n_wann];
    for mw_k in &mw {
        for (ib, bv) in bset.vectors.iter().enumerate() {
            let m = &mw_k[ib];
            for n in 0..n_wann {
                let diag = m[(n, n)];
                let phase = diag.im.atan2(diag.re);
                let b_dot_r = dot3(bv.cartesian, centers[n]);
                let q = phase + b_dot_r;
                orbital_spreads[n] += inv_nk * bv.weight * ((1.0 - diag.norm_sqr()) + q * q);
            }
        }
    }
    let omega_total: f64 = orbital_spreads.iter().sum();
    let omega_tilde = omega_total - omega_i;

    Ok(SpreadReport {
        omega_i_bohr2: omega_i,
        omega_tilde_bohr2: omega_tilde,
        omega_total_bohr2: omega_total,
        orbital_spreads_bohr2: orbital_spreads,
        centers_bohr: centers,
    })
}

/// Build an [`ImportedKProjection`] from an imported Wannier90 gauge matrix
/// `gauge` (`N_bands x N_wann`, already multiplied `U_dis * U` if disentangled)
/// and NSCF `eigenvalues` (`N_bands`, Ha).
pub fn project_with_gauge(
    eigenvalues: &[f64],
    gauge: ComplexMatrix,
    orth_tol: f64,
) -> Result<ImportedKProjection, ResponseError> {
    let n_bands = eigenvalues.len();
    let n_wann = gauge.cols();
    if n_bands == 0 || n_wann == 0 || gauge.rows() != n_bands || n_wann > n_bands {
        return Err(input(format!(
            "wannier: imported gauge shape ({}x{}) incompatible with {} bands",
            gauge.rows(),
            n_wann,
            n_bands
        )));
    }
    let gram = gemm(&gauge, MatrixOp::Adjoint, &gauge, MatrixOp::Identity).map_err(numerical)?;
    let mut orthonormality_defect = 0.0_f64;
    for j in 0..n_wann {
        for i in 0..n_wann {
            let target = if i == j {
                Complex64::new(1.0, 0.0)
            } else {
                Complex64::default()
            };
            orthonormality_defect = orthonormality_defect.max((gram[(i, j)] - target).norm());
        }
    }
    if !orthonormality_defect.is_finite() || orthonormality_defect > orth_tol {
        return Err(numerical(format!(
            "imported Wannier90 gauge is not semi-unitary: max |U^H U - I| = {orthonormality_defect:.3e} > {orth_tol:.3e}"
        )));
    }

    let mut band_weights = vec![0.0_f64; n_bands];
    for (m, w) in band_weights.iter_mut().enumerate() {
        *w = (0..n_wann).map(|j| gauge[(m, j)].norm_sqr()).sum();
    }

    let mut eu = ComplexMatrix::zeros(n_bands, n_wann).map_err(numerical)?;
    for m in 0..n_bands {
        let e = eigenvalues[m];
        for j in 0..n_wann {
            eu[(m, j)] = gauge[(m, j)] * e;
        }
    }
    let mut hamiltonian_k =
        gemm(&gauge, MatrixOp::Adjoint, &eu, MatrixOp::Identity).map_err(numerical)?;
    let mut hermiticity_defect = 0.0_f64;
    for j in 0..n_wann {
        let d = hamiltonian_k[(j, j)];
        hermiticity_defect = hermiticity_defect.max(d.im.abs());
        hamiltonian_k[(j, j)] = Complex64::new(d.re, 0.0);
        for i in (j + 1)..n_wann {
            let a = hamiltonian_k[(i, j)];
            let b = hamiltonian_k[(j, i)];
            hermiticity_defect = hermiticity_defect.max((a - b.conj()).norm() * 0.5);
            let avg = (a + b.conj()) * 0.5;
            hamiltonian_k[(i, j)] = avg;
            hamiltonian_k[(j, i)] = avg.conj();
        }
    }

    Ok(ImportedKProjection {
        gauge,
        orthonormality_defect,
        band_weights,
        hamiltonian_k,
        hermiticity_defect,
    })
}
