//! R3 analytic oracle of the ACFDT pipeline (accumulator -> Cholesky log det
//! -> mapped Gauss-Legendre frequency integral) without any UPF.
//!
//! For T transitions with pair vectors rho_j (N_G rows), gaps Delta_j > 0 and
//! weight s, X(i w) = -R F(w) R^H with R = v^1/2 [rho_j],
//! F = diag(2 s Delta_j / (Delta_j^2 + w^2)). The determinant lemma gives
//! det(1 - X) = prod_i (Omega_i^2 + w^2) / prod_j (Delta_j^2 + w^2) with
//! Omega_i^2 the eigenvalues of D^2 + 2 D^1/2 K D^1/2, K = s R^H R, and
//! int_0^inf ln((w^2 + a^2)/(w^2 + b^2)) dw = pi (a - b), so the exact
//! result is the plasmon formula
//!   E_c = (1/2) (sum_i Omega_i - sum_j Delta_j - Tr K).
//! Fast gate (debug): `cargo test -p dft-response --test acfdt -- --nocapture`
use dft_numerics::linalg::{ComplexMatrix, hermitian_eigh};
use dft_response::{acfdt, chi0::Accumulator, frequency::FrequencyGrid};
use num_complex::Complex64 as C;
use std::f64::consts::PI;

struct Model {
    pairs: ComplexMatrix,
    delta: Vec<f64>,
    norms: Vec<f64>,
    scale: f64,
}

fn model(n_g: usize, t: usize, seed: u64) -> Model {
    let mut state = seed;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut pairs = ComplexMatrix::zeros(n_g, t).unwrap();
    for z in pairs.as_mut_slice() {
        *z = C::new(next() - 0.5, next() - 0.5) * 0.3;
    }
    let delta = (0..t).map(|_| -(0.15 + 1.5 * next())).collect();
    let norms = (0..n_g).map(|_| 0.6 + 2. * next()).collect();
    Model {
        pairs,
        delta,
        norms,
        scale: 0.2,
    }
}

/// Closed form (plasmon formula), Ha per cell.
fn exact(m: &Model) -> f64 {
    let t = m.delta.len();
    let v: Vec<f64> = m.norms.iter().map(|g| 4. * PI / (g * g)).collect();
    let gaps: Vec<f64> = m.delta.iter().map(|d| -d).collect();
    let mut k = ComplexMatrix::zeros(t, t).unwrap();
    for i in 0..t {
        for j in 0..t {
            k[(i, j)] = m.scale
                * (0..v.len())
                    .map(|g| m.pairs[(g, i)].conj() * v[g] * m.pairs[(g, j)])
                    .sum::<C>();
        }
    }
    let mut a = ComplexMatrix::zeros(t, t).unwrap();
    for i in 0..t {
        for j in 0..t {
            a[(i, j)] = 2. * gaps[i].sqrt() * k[(i, j)] * gaps[j].sqrt();
        }
        a[(i, i)] += gaps[i] * gaps[i];
    }
    // Hermitize the roundoff of the product.
    let mut h = a.clone();
    for i in 0..t {
        for j in 0..t {
            h[(i, j)] = 0.5 * (a[(i, j)] + a[(j, i)].conj());
        }
    }
    let (omega2, _) = hermitian_eigh(&h).unwrap();
    let trace_k: f64 = (0..t).map(|i| k[(i, i)].re).sum();
    0.5 * (omega2.iter().map(|x| x.sqrt()).sum::<f64>() - gaps.iter().sum::<f64>() - trace_k)
}

fn acfdt_energy(m: &Model, points: usize, omega0: f64, split: usize) -> f64 {
    let grid = FrequencyGrid::mapped_gauss_legendre(points, omega0).unwrap();
    let n_g = m.norms.len();
    let mut acc = Accumulator::new(n_g, &grid.points, false, 3).unwrap();
    // Several add calls (as the band blocks of the driver).
    let t = m.delta.len();
    let mut start = 0;
    while start < t {
        let end = (start + split).min(t);
        let mut block = ComplexMatrix::zeros(n_g, end - start).unwrap();
        for j in start..end {
            block
                .column_mut(j - start)
                .copy_from_slice(m.pairs.column(j));
        }
        acc.add(&block, &m.delta[start..end], m.scale, None)
            .unwrap();
        start = end;
    }
    let chi0 = acc.reduce(&dft_mpi::Comm::Single).unwrap();
    let terms: Vec<f64> = chi0
        .body
        .iter()
        .map(|b| {
            let t = acfdt::finite_q_term(b, &m.norms).unwrap();
            assert!(t.min_pivot >= 1. - 1e-12);
            t.value
        })
        .collect();
    acfdt::frequency_integral(&grid.weights, &terms).unwrap()
}

#[test]
fn acfdt_matches_the_plasmon_formula() {
    for (n_g, t, seed) in [(1, 1, 1), (3, 7, 2), (12, 30, 3), (25, 10, 4)] {
        let m = model(n_g, t, seed);
        let reference = exact(&m);
        let gaps: Vec<f64> = m.delta.iter().map(|d| -d).collect();
        let lo = gaps.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = gaps.iter().copied().fold(0., f64::max);
        let omega0 = 0.5 * (lo + hi);
        let mut previous = f64::INFINITY;
        for points in [4, 8, 16, 32, 64] {
            let e = acfdt_energy(&m, points, omega0, 4);
            let err = (e - reference).abs();
            println!(
                "N_G={n_g} T={t}: N_w={points:2} E_c={e:.15e} exact={reference:.15e} |err|={err:.3e} Ha"
            );
            // Absolute roundoff floor: ln det(1 - X) + Tr X cancels to
            // ~1e-16 of Tr X per node, and the tail nodes carry large weights.
            let floor = 1e-12 * gaps.iter().sum::<f64>();
            if points >= 16 {
                assert!(err <= previous * 1.0001 + floor, "not converging");
            }
            previous = err;
        }
        assert!(reference < 0.);
        assert!(
            previous <= (1e-10 * reference.abs()).max(1e-12 * gaps.iter().sum::<f64>()),
            "64-point error {previous:e} vs {reference:e}"
        );
        // Block split does not change the result beyond roundoff.
        let a = acfdt_energy(&m, 32, omega0, 1);
        let b = acfdt_energy(&m, 32, omega0, 64);
        assert!((a - b).abs() <= 1e-14 * a.abs().max(1e-12) + 1e-17);
    }
}

#[test]
fn eigenvalue_path_agrees_on_accumulated_chi0() {
    let m = model(20, 60, 7);
    let mut acc = Accumulator::new(20, &[0.01, 0.4, 3.0], false, 1).unwrap();
    acc.add(&m.pairs, &m.delta, m.scale, None).unwrap();
    let chi0 = acc.reduce(&dft_mpi::Comm::Single).unwrap();
    for b in &chi0.body {
        let c = acfdt::finite_q_term(b, &m.norms).unwrap();
        let e = acfdt::finite_q_term_eigen(b, &m.norms).unwrap();
        println!(
            "log det {:.6e} trace {:.6e} value {:.15e} eigen {e:.15e} herm {:.1e}",
            c.log_det, c.trace, c.value, c.hermiticity
        );
        assert!((c.value - e).abs() <= 1e-12 * c.log_det.max(1.));
        assert!(c.hermiticity <= 1e-14);
    }
}
