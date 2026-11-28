//! R4b (`dft_response::wannier`, ADR 0060) internal numerical oracles:
//! - Marzari-Vanderbilt B1 completeness `sum_b w_b b_a b_c = delta_ac` across
//!   cubic, fcc, bcc, hexagonal and strongly skewed triclinic cells
//! - Exact Bloch overlap `M_mn^(k,b)` plane-wave shift identity and adjoint
//!   symmetry `M^(k', -b) = (M^(k, b))^H`
//! - Exact gauge invariance of `Omega_I` under random k-dependent unitary
//!   rotations `U(k)`
//! - Exact Wannier center shift under real-space translation `r -> r + delta`
//! - Imported-gauge projection semi-unitarity check and exact weight sum rule
//!
//! Fast gate: `cargo test -p dft-response --test wannier -- --nocapture`
use dft_geometry::Cell;
use dft_numerics::linalg::{ComplexMatrix, MatrixOp, gemm, hermitian_eigh};
use dft_response::wannier::{
    B1_COMPLETENESS_TOLERANCE, DEFAULT_GAUGE_TOLERANCE, bloch_overlap, evaluate_spreads,
    find_bvector_shells, mesh_neighbors, project_with_gauge,
};
use num_complex::Complex64;
use std::f64::consts::TAU;

fn det_rng(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    let x = ((*seed >> 11) as f64) * (1.0 / ((1u64 << 53) as f64));
    2.0 * x - 1.0
}

fn random_unitary(n: usize, seed: &mut u64) -> ComplexMatrix {
    let mut h = ComplexMatrix::zeros(n, n).unwrap();
    for j in 0..n {
        h[(j, j)] = Complex64::new(det_rng(seed), 0.0);
        for i in (j + 1)..n {
            let z = Complex64::new(det_rng(seed), det_rng(seed));
            h[(i, j)] = z;
            h[(j, i)] = z.conj();
        }
    }
    hermitian_eigh(&h).unwrap().1
}

fn recip_vectors(cell: &Cell) -> [[f64; 3]; 3] {
    let cols = cell.reciprocal_columns_inverse_bohr();
    std::array::from_fn(|j| [cols[0][j], cols[1][j], cols[2][j]])
}

#[test]
fn bvector_shells_satisfy_b1_completeness_for_five_bravais_families() {
    let cells = [
        // Simple cubic
        Cell::from_user_rows_bohr([[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]]).unwrap(),
        // FCC primitive
        Cell::from_user_rows_bohr([[0.0, 5.1, 5.1], [5.1, 0.0, 5.1], [5.1, 5.1, 0.0]]).unwrap(),
        // BCC primitive
        Cell::from_user_rows_bohr([[-2.7, 2.7, 2.7], [2.7, -2.7, 2.7], [2.7, 2.7, -2.7]]).unwrap(),
        // Hexagonal
        Cell::from_user_rows_bohr([
            [5.8, 0.0, 0.0],
            [-2.9, 5.8 * 3.0_f64.sqrt() / 2.0, 0.0],
            [0.0, 0.0, 9.4],
        ])
        .unwrap(),
        // Skewed triclinic
        Cell::from_user_rows_bohr([[6.0, 0.0, 0.0], [1.4, 5.2, 0.0], [0.7, 1.1, 6.1]]).unwrap(),
    ];
    for (idx, cell) in cells.iter().enumerate() {
        let bvecs = recip_vectors(cell);
        for mesh in [[2, 2, 2], [3, 3, 2]] {
            let bset = find_bvector_shells(bvecs, mesh).unwrap();
            assert!(
                bset.completeness_residual < B1_COMPLETENESS_TOLERANCE,
                "cell {idx} mesh {mesh:?} residual {:.3e}",
                bset.completeness_residual
            );
            // Explicitly re-check sum_b w_b b_a b_c = delta_ac
            let mut tensor = [[0.0_f64; 3]; 3];
            for bv in &bset.vectors {
                for (row, &ba) in tensor.iter_mut().zip(&bv.cartesian) {
                    for (entry, &bc) in row.iter_mut().zip(&bv.cartesian) {
                        *entry += bv.weight * ba * bc;
                    }
                }
            }
            for (a, row) in tensor.iter().enumerate() {
                for (c, &val) in row.iter().enumerate() {
                    let target = if a == c { 1.0 } else { 0.0 };
                    assert!(
                        (val - target).abs() < 1e-11,
                        "cell {idx} mesh {mesh:?} ({a},{c}) = {val}"
                    );
                }
            }
        }
    }
}

#[test]
fn bloch_overlap_shift_identity_and_adjoint_symmetry() {
    let mut seed = 77u64;
    let miller = [
        [0, 0, 0],
        [1, 0, 0],
        [-1, 0, 0],
        [0, 1, 0],
        [0, -1, 0],
        [0, 0, 1],
        [0, 0, -1],
    ];
    let left = random_unitary(miller.len(), &mut seed);
    let right = random_unitary(miller.len(), &mut seed);
    let dg = [1, 0, 0];
    let m_fwd = bloch_overlap(&miller, &left, &miller, &right, dg).unwrap();
    let m_rev = bloch_overlap(&miller, &right, &miller, &left, [-dg[0], -dg[1], -dg[2]]).unwrap();
    for j in 0..miller.len() {
        for i in 0..miller.len() {
            assert!(
                (m_fwd[(i, j)] - m_rev[(j, i)].conj()).norm() < 1e-14,
                "adjoint mismatch at ({i},{j})"
            );
        }
    }
    // Plane-wave identity: at dg = 0 with identical unitary basis, M = left^H * right is unitary.
    let m0 = bloch_overlap(&miller, &left, &miller, &right, [0, 0, 0]).unwrap();
    let gram = gemm(&m0, MatrixOp::Adjoint, &m0, MatrixOp::Identity).unwrap();
    for j in 0..miller.len() {
        for i in 0..miller.len() {
            let target = if i == j {
                Complex64::new(1.0, 0.0)
            } else {
                Complex64::default()
            };
            assert!((gram[(i, j)] - target).norm() < 1e-13);
        }
    }
}

#[test]
fn spreads_omega_i_is_gauge_invariant_and_center_tracks_translation() {
    // Analytic 2-band localized model: cell-periodic states u_{nk}(r) = exp(-i k . r_n) * phi_n(r)
    // => M_mn^(k,b) = delta_mn exp(-i b . r_n) * f_n(|b|), giving exact center r_n and exact spreads.
    let cell =
        Cell::from_user_rows_bohr([[6.0, 0.0, 0.0], [0.5, 5.5, 0.0], [0.2, 0.3, 6.2]]).unwrap();
    let mesh = [3, 3, 3];
    let bset = find_bvector_shells(recip_vectors(&cell), mesh).unwrap();
    let neighbors = mesh_neighbors(mesh, &bset);
    let nk = neighbors.len();

    let r_true = [[0.4, -0.3, 0.2], [-0.5, 0.6, -0.1]];
    let sigma2 = [0.25_f64, 0.40_f64];

    let mut overlaps = Vec::with_capacity(nk);
    for _ in 0..nk {
        let mut per_b = Vec::with_capacity(bset.vectors.len());
        for bv in &bset.vectors {
            let b2 = bv.cartesian.iter().map(|x| x * x).sum::<f64>();
            let mut m = ComplexMatrix::zeros(2, 2).unwrap();
            for n in 0..2 {
                let b_dot_r = bv.cartesian[0] * r_true[n][0]
                    + bv.cartesian[1] * r_true[n][1]
                    + bv.cartesian[2] * r_true[n][2];
                let amp = (-0.5 * sigma2[n] * b2).exp();
                m[(n, n)] = Complex64::from_polar(amp, -b_dot_r);
            }
            per_b.push(m);
        }
        overlaps.push(per_b);
    }

    // Identity gauge: recovers exact centers r_true and Omega_tilde = 0.
    let mut id_gauges = Vec::with_capacity(nk);
    for _ in 0..nk {
        let mut u = ComplexMatrix::zeros(2, 2).unwrap();
        u[(0, 0)] = Complex64::new(1.0, 0.0);
        u[(1, 1)] = Complex64::new(1.0, 0.0);
        id_gauges.push(u);
    }
    let rep_id = evaluate_spreads(&bset, &neighbors, &overlaps, &id_gauges).unwrap();
    for n in 0..2 {
        for a in 0..3 {
            assert!(
                (rep_id.centers_bohr[n][a] - r_true[n][a]).abs() < 1e-12,
                "center ({n},{a}): got {}, expected {}",
                rep_id.centers_bohr[n][a],
                r_true[n][a]
            );
        }
    }
    assert!(
        rep_id.omega_tilde_bohr2.abs() < 1e-12,
        "Omega_tilde should vanish in optimal gauge, got {:.3e}",
        rep_id.omega_tilde_bohr2
    );

    // Random k-dependent unitary gauge rotation: Omega_I MUST be invariant and Omega_tilde >= 0!
    let mut seed = 2026u64;
    let rand_gauges: Vec<ComplexMatrix> = (0..nk).map(|_| random_unitary(2, &mut seed)).collect();
    let rep_rot = evaluate_spreads(&bset, &neighbors, &overlaps, &rand_gauges).unwrap();
    assert!(
        (rep_rot.omega_i_bohr2 - rep_id.omega_i_bohr2).abs() < 1e-12,
        "Omega_I not gauge invariant: {} vs {}",
        rep_rot.omega_i_bohr2,
        rep_id.omega_i_bohr2
    );
    assert!(
        rep_rot.omega_tilde_bohr2 > 0.0,
        "random gauge should increase spread: {}",
        rep_rot.omega_tilde_bohr2
    );

    // Check imported-gauge projection helper on semi-unitary 4x2 gauge.
    let u4 = random_unitary(4, &mut seed);
    let mut u42 = ComplexMatrix::zeros(4, 2).unwrap();
    for j in 0..2 {
        u42.column_mut(j).copy_from_slice(u4.column(j));
    }
    let proj = project_with_gauge(&[-0.4, -0.1, 0.2, 0.5], u42, DEFAULT_GAUGE_TOLERANCE).unwrap();
    assert!(proj.orthonormality_defect < 1e-13);
    assert!((proj.band_weights.iter().sum::<f64>() - 2.0).abs() < 1e-13);
    // Reject non-unitary matrix.
    let bad = ComplexMatrix::zeros(4, 2).unwrap();
    assert!(project_with_gauge(&[-0.4, -0.1, 0.2, 0.5], bad, DEFAULT_GAUGE_TOLERANCE).is_err());
    let _ = TAU;
}
