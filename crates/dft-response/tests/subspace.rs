//! R4a (`dft_response::subspace`, ADR 0059) internal numerical oracles:
//! - Loewdin orthonormality defect <= 1e-12 for isolated and entangled windows
//! - Exact KS eigenvalue recovery for isolated manifolds (`|W| = N_orb`)
//! - Exact subspace weight sum rule (`sum_m p_mk = N_orb`) and bounds `[0, 1]`
//! - Exact Wigner-Seitz partition of unity (`sum_R 1/g(R) = N_k`) on cubic,
//!   fcc, hexagonal and strongly skewed triclinic cells
//! - Real-space hopping Hermiticity `H_ij(R) = H_ji^*(-R)` and exact
//!   interpolation back onto the uniform k-mesh
//! - Analytic tight-binding model exact recovery at arbitrary off-mesh k-points
//! - Singular overlap rejection
//!
//! Fast gate: `cargo test -p dft-response --test subspace -- --nocapture`
use dft_numerics::linalg::{ComplexMatrix, hermitian_eigh};
use dft_response::subspace::{
    DEFAULT_SINGULAR_TOLERANCE, build_real_space_hamiltonian, project_k, projected_populations,
    wigner_seitz_points,
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

/// Deterministic orthonormal Bloch matrix via Hermitian EVD of a random matrix.
fn random_orthonormal(n_pw: usize, n_bands: usize, seed: &mut u64) -> ComplexMatrix {
    let mut h = ComplexMatrix::zeros(n_pw, n_pw).unwrap();
    for j in 0..n_pw {
        h[(j, j)] = Complex64::new(det_rng(seed) * 2.0, 0.0);
        for i in (j + 1)..n_pw {
            let z = Complex64::new(det_rng(seed), det_rng(seed));
            h[(i, j)] = z;
            h[(j, i)] = z.conj();
        }
    }
    let (_, evecs) = hermitian_eigh(&h).unwrap();
    let mut out = ComplexMatrix::zeros(n_pw, n_bands).unwrap();
    for j in 0..n_bands {
        out.column_mut(j).copy_from_slice(evecs.column(j));
    }
    out
}

#[test]
fn isolated_manifold_is_unitary_and_recovers_eigenvalues() {
    let mut seed = 42u64;
    let n_pw = 24;
    let n_bands = 10;
    let bloch = random_orthonormal(n_pw, n_bands, &mut seed);
    let eigenvalues: Vec<f64> = (0..n_bands).map(|i| -0.5 + 0.13 * i as f64).collect();
    let window = [2, 3, 4, 5];
    let n_orb = window.len();

    // Build trial orbitals as a non-orthogonal, full-rank linear combination
    // of the window Bloch states plus an orthogonal component outside the window.
    let mut trial = ComplexMatrix::zeros(n_pw, n_orb).unwrap();
    for col in 0..n_orb {
        for (row_idx, &b) in window.iter().enumerate() {
            let c = if row_idx == col {
                Complex64::new(0.85, 0.1 * col as f64)
            } else {
                Complex64::new(0.15 * det_rng(&mut seed), 0.1 * det_rng(&mut seed))
            };
            for g in 0..n_pw {
                trial[(g, col)] += bloch[(g, b)] * c;
            }
        }
        // Add small leakage into band 8 (outside window) so spilling > 0.
        for g in 0..n_pw {
            trial[(g, col)] += bloch[(g, 8)] * Complex64::new(0.12, -0.05);
        }
    }

    let proj = project_k(
        &bloch,
        &eigenvalues,
        &trial,
        &window,
        DEFAULT_SINGULAR_TOLERANCE,
    )
    .unwrap();
    assert!(
        proj.orthonormality_defect < 1e-13,
        "orthonormality defect = {:.3e}",
        proj.orthonormality_defect
    );
    assert!(
        proj.hermiticity_defect < 1e-14,
        "hermiticity defect = {:.3e}",
        proj.hermiticity_defect
    );
    // Exact eigenvalue preservation for an isolated manifold.
    let (h_evals, _) = hermitian_eigh(&proj.hamiltonian_k).unwrap();
    for (&got, &b) in h_evals.iter().zip(&window) {
        assert!(
            (got - eigenvalues[b]).abs() < 1e-13,
            "band {b}: got {got}, expected {}",
            eigenvalues[b]
        );
    }
    // Every band in an isolated manifold has weight 1.0, outside 0.0.
    for (b, &w) in proj.band_weights.iter().enumerate() {
        let expected = if window.contains(&b) { 1.0 } else { 0.0 };
        assert!(
            (w - expected).abs() < 1e-13,
            "band {b} weight = {w}, expected {expected}"
        );
    }
}

#[test]
fn entangled_window_weights_and_populations() {
    let mut seed = 1337u64;
    let n_pw = 20;
    let n_bands = 10;
    let bloch = random_orthonormal(n_pw, n_bands, &mut seed);
    let eigenvalues: Vec<f64> = (0..n_bands).map(|i| -0.4 + 0.11 * i as f64).collect();
    let window = [1, 2, 3, 4, 5, 6];
    let n_orb = 3;

    let mut trial = ComplexMatrix::zeros(n_pw, n_orb).unwrap();
    for col in 0..n_orb {
        for g in 0..n_pw {
            trial[(g, col)] = Complex64::new(det_rng(&mut seed), det_rng(&mut seed));
        }
        let nrm: f64 = trial
            .column(col)
            .iter()
            .map(|z| z.norm_sqr())
            .sum::<f64>()
            .sqrt();
        for z in trial.column_mut(col) {
            *z /= nrm;
        }
    }

    let proj = project_k(
        &bloch,
        &eigenvalues,
        &trial,
        &window,
        DEFAULT_SINGULAR_TOLERANCE,
    )
    .unwrap();
    assert!(proj.orthonormality_defect < 1e-13);
    assert!(proj.spilling >= 0.0 && proj.spilling < 1.0);
    let weight_sum: f64 = proj.band_weights.iter().sum();
    assert!(
        (weight_sum - n_orb as f64).abs() < 1e-13,
        "sum p_m = {weight_sum}, expected {n_orb}"
    );
    for (b, &w) in proj.band_weights.iter().enumerate() {
        if window.contains(&b) {
            assert!((-1e-14..=1.0 + 1e-14).contains(&w), "band {b} weight {w}");
        } else {
            assert_eq!(w, 0.0);
        }
    }
    // Populations when all window bands are fully occupied: sum_i lowdin_i = spin * N_orb.
    let occ = vec![1.0; n_bands];
    let (_, lowdin) = projected_populations(&[1.0], &[occ.as_slice()], &[proj], 2.0).unwrap();
    let total_lowdin: f64 = lowdin.iter().sum();
    assert!((total_lowdin - 2.0 * n_orb as f64).abs() < 1e-13);
}

#[test]
fn singular_overlap_rejects() {
    let mut seed = 9u64;
    let bloch = random_orthonormal(12, 6, &mut seed);
    let eigenvalues = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
    let mut trial = ComplexMatrix::zeros(12, 2).unwrap();
    // Duplicate columns -> rank 1 overlap matrix.
    trial.column_mut(0).copy_from_slice(bloch.column(0));
    trial.column_mut(1).copy_from_slice(bloch.column(0));
    let err = project_k(
        &bloch,
        &eigenvalues,
        &trial,
        &[0, 1, 2],
        DEFAULT_SINGULAR_TOLERANCE,
    )
    .unwrap_err();
    assert!(
        err.to_string().contains("rank-deficient"),
        "unexpected error: {err}"
    );
}

#[test]
fn wigner_seitz_partition_of_unity_across_cell_shapes() {
    let cells = [
        // Simple cubic
        [[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]],
        // FCC primitive
        [[0.0, 5.1, 5.1], [5.1, 0.0, 5.1], [5.1, 5.1, 0.0]],
        // Hexagonal
        [
            [5.8, 0.0, 0.0],
            [-2.9, 5.8 * 3.0_f64.sqrt() / 2.0, 0.0],
            [0.0, 0.0, 9.4],
        ],
        // Strongly skewed triclinic
        [[6.0, 0.0, 0.0], [3.5, 4.2, 0.0], [2.1, 2.8, 5.1]],
    ];
    for cell in cells {
        for mesh in [[2, 2, 2], [3, 3, 3], [4, 3, 2]] {
            let pts = wigner_seitz_points(cell, mesh).unwrap();
            let sum: f64 = pts.iter().map(|p| 1.0 / p.degeneracy as f64).sum();
            let expected = (mesh[0] * mesh[1] * mesh[2]) as f64;
            assert!(
                (sum - expected).abs() < 1e-12,
                "cell={cell:?} mesh={mesh:?} sum={sum} expected={expected}"
            );
            // Origin R = [0,0,0] must be present with degeneracy 1.
            let origin = pts.iter().find(|p| p.r == [0, 0, 0]).unwrap();
            assert_eq!(origin.degeneracy, 1);
        }
    }
}

#[test]
fn real_space_hamiltonian_exact_on_mesh_and_analytic_tb_model() {
    // Analytic 2-band nearest-neighbor tight-binding model on a skewed cell:
    // H(k) = H0 + H_x e^{i 2pi k_x} + H_x^H e^{-i 2pi k_x} + H_y e^{i 2pi k_y} + H_y^H e^{-i 2pi k_y}
    // A 3x3x3 mesh resolves nearest-neighbor hoppings R = +/-e_x, +/-e_y with
    // zero aliasing, so Wannier interpolation must match H(k') at ANY off-mesh k'.
    let cell = [[6.0, 0.0, 0.0], [1.2, 5.5, 0.0], [0.4, 0.6, 6.2]];
    let mesh = [3, 3, 3];
    let h0 = ComplexMatrix::new(
        2,
        2,
        vec![
            Complex64::new(-0.30, 0.0),
            Complex64::new(0.05, 0.02),
            Complex64::new(0.05, -0.02),
            Complex64::new(0.20, 0.0),
        ],
    )
    .unwrap();
    let hx = ComplexMatrix::new(
        2,
        2,
        vec![
            Complex64::new(-0.08, 0.01),
            Complex64::new(0.03, -0.04),
            Complex64::new(-0.02, 0.01),
            Complex64::new(0.06, -0.01),
        ],
    )
    .unwrap();
    let hy = ComplexMatrix::new(
        2,
        2,
        vec![
            Complex64::new(-0.05, -0.02),
            Complex64::new(0.01, 0.02),
            Complex64::new(0.04, -0.01),
            Complex64::new(0.07, 0.0),
        ],
    )
    .unwrap();

    let eval_tb = |k: [f64; 3]| -> ComplexMatrix {
        let px = Complex64::from_polar(1.0, TAU * k[0]);
        let py = Complex64::from_polar(1.0, TAU * k[1]);
        let mut h = h0.clone();
        for j in 0..2 {
            for i in 0..2 {
                h[(i, j)] += hx[(i, j)] * px
                    + hx[(j, i)].conj() * px.conj()
                    + hy[(i, j)] * py
                    + hy[(j, i)].conj() * py.conj();
            }
        }
        h
    };

    let mut ks = Vec::new();
    let mut hks = Vec::new();
    for i in 0..mesh[0] {
        for j in 0..mesh[1] {
            for l in 0..mesh[2] {
                let k = [
                    i as f64 / mesh[0] as f64,
                    j as f64 / mesh[1] as f64,
                    l as f64 / mesh[2] as f64,
                ];
                ks.push(k);
                hks.push(eval_tb(k));
            }
        }
    }

    let hr = build_real_space_hamiltonian(cell, mesh, &ks, &hks).unwrap();
    assert!(
        hr.hermiticity_defect() < 1e-14,
        "H(R) Hermiticity defect = {:.3e}",
        hr.hermiticity_defect()
    );

    // Exact on every mesh k-point AND at arbitrary off-mesh k-points!
    for &k in ks
        .iter()
        .chain(&[[0.17, 0.43, 0.29], [0.5, 0.25, 0.125], [0.333, 0.667, 0.0]])
    {
        let expected_h = eval_tb(k);
        let got_h = hr.interpolate_matrix(k).unwrap();
        for j in 0..2 {
            for i in 0..2 {
                assert!(
                    (got_h[(i, j)] - expected_h[(i, j)]).norm() < 1e-13,
                    "k={k:?} ({i},{j}): got {}, expected {}",
                    got_h[(i, j)],
                    expected_h[(i, j)]
                );
            }
        }
        let (exp_e, _) = hermitian_eigh(&expected_h).unwrap();
        let got_e = hr.interpolate_eigenvalues(k).unwrap();
        for (a, b) in got_e.iter().zip(&exp_e) {
            assert!((a - b).abs() < 1e-13);
        }
    }
}
