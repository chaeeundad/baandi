//! B5 gate microbenchmark: the existing reduced eigensolve
//! (`RealEighWorkspace::solve_into`, the backend the real Davidson uses) at the
//! orders the solver meets, and one m x m product for the bordered-update cost
//! model. Observation only (Mac or dev08), no assertion on time.
//!
//! `cargo test --release -p dft-solvers --test b5_reduced_eigen -- --ignored --nocapture`
//!
//! The non-ignored tests check B5 (`ResearchOptions::b5`) inside the real
//! Davidson: same converged roots (1e-10 Ha), certified strict returns, work
//! within 10 percent, thread-count independence, and that the bordered path
//! was taken.
use dft_numerics::pool::WorkerPool;
use dft_numerics::real::{RealEighWorkspace, RealMatrix};
use dft_solvers::real_davidson::{RealEigenSolution, ResearchOptions, b5_counters};
use dft_solvers::{
    BandTolerances, DavidsonPreconditioner, RealDavidson, RealLinearOperator, SolverError,
};
use std::sync::Arc;
use std::time::Instant;

struct Operator {
    a: RealMatrix,
    pool: Arc<WorkerPool>,
}
impl RealLinearOperator for Operator {
    fn dimension(&self) -> usize {
        self.a.rows()
    }
    fn worker_pool(&self) -> Option<Arc<WorkerPool>> {
        Some(self.pool.clone())
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), SolverError> {
        let n = self.a.rows();
        y.fill(0.);
        for (j, &xj) in x.iter().enumerate() {
            for (y, &a) in y.iter_mut().zip(&self.a.as_slice()[j * n..(j + 1) * n]) {
                *y += a * xj;
            }
        }
        Ok(())
    }
}

/// Kinetic-like diagonal (degenerate pairs) with decaying coupling, Ha.
fn kinetic(n: usize, coupling: f64, threads: usize) -> (Operator, Vec<f64>) {
    let diag: Vec<f64> = (0..n).map(|i| 0.2 * (i / 2) as f64).collect();
    let mut a = RealMatrix::zeros(n, n).unwrap();
    for j in 0..n {
        for i in 0..=j {
            let v = if i == j {
                diag[i] + 0.3
            } else {
                coupling * ((i * 31 + j * 17 + 1) as f64).sin() / ((j - i) as f64).sqrt()
            };
            a[(i, j)] = v;
            a[(j, i)] = v;
        }
    }
    (
        Operator {
            a,
            pool: WorkerPool::new(threads).unwrap(),
        },
        diag,
    )
}

fn solve(op: &Operator, diag: &[f64], ritz: bool, b5: bool) -> RealEigenSolution {
    let s = RealDavidson {
        bands: 16,
        tolerance: if ritz { 1e-10 } else { 1e-9 },
        max_iterations: 400,
        max_subspace: 96,
        seed: 42,
    };
    let policy = BandTolerances {
        strict_bands: 12,
        buffer_tolerance: if ritz { 1e-8 } else { 1e-6 },
    };
    s.solve_research(
        op,
        None,
        diag,
        policy,
        ritz,
        !ritz,
        6,
        &mut None,
        DavidsonPreconditioner::Default,
        ResearchOptions {
            b5,
            ..ResearchOptions::default()
        },
    )
    .unwrap()
}

#[test]
fn b5_keeps_roots_certificates_and_thread_identity() {
    let before = b5_counters();
    for coupling in [0.025, 0.08] {
        let (op1, diag) = kinetic(220, coupling, 1);
        let (op4, _) = kinetic(220, coupling, 4);
        for ritz in [false, true] {
            let base = solve(&op1, &diag, ritz, false);
            let fast = solve(&op1, &diag, ritz, true);
            assert_eq!(fast, solve(&op4, &diag, ritz, true), "thread identity");
            let near = |a: usize, b: usize| a.abs_diff(b) * 10 <= a;
            assert!(near(base.iterations, fast.iterations));
            assert!(near(base.operator_applications, fast.operator_applications));
            for (a, b) in base.eigenvalues.iter().zip(&fast.eigenvalues) {
                assert!((a - b).abs() <= 1e-10, "{a} vs {b}");
            }
            assert!(fast.converged.iter().all(|&c| c));
            if !ritz {
                assert!(fast.residuals_certified);
                assert!(fast.residuals[..12].iter().all(|&r| r <= 1e-9));
            }
        }
    }
    let after = b5_counters();
    assert!(after.bordered > before.bordered, "bordered path used");
}

fn symmetric(m: usize, seed: u64) -> RealMatrix {
    let mut state = seed;
    let mut a = RealMatrix::zeros(m, m).unwrap();
    for j in 0..m {
        for i in 0..=j {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            let r = ((state >> 11) as f64) / (1_u64 << 53) as f64 - 0.5;
            // Ritz-like: growing diagonal, small coupling.
            let v = if i == j { 0.05 * i as f64 } else { 1e-3 * r };
            a[(i, j)] = v;
            a[(j, i)] = v;
        }
    }
    a
}

#[test]
#[ignore = "B5 gate timing: release, --nocapture"]
fn reduced_eigensolve_timing() {
    println!("order | eigensolve ms (median of 9) | m x m product ms (median of 9)");
    for m in [48, 96, 136, 204, 272, 340, 408] {
        let a = symmetric(m, 11 + m as u64);
        let mut ws = RealEighWorkspace::new(m).unwrap();
        let mut values = vec![0.; m];
        let mut vectors = RealMatrix::zeros(m, m).unwrap();
        let mut eig = Vec::new();
        for _ in 0..9 {
            let t = Instant::now();
            ws.solve_into(&a, &mut values, &mut vectors).unwrap();
            eig.push(t.elapsed().as_secs_f64() * 1e3);
        }
        let mut product = Vec::new();
        let mut out = vec![0.; m * m];
        for _ in 0..9 {
            let t = Instant::now();
            // Plain column-major m x m x m product (the bordered update's
            // C [Y] step), cache-friendly loop order, no BLAS.
            out.fill(0.);
            let c = vectors.as_slice();
            let y = a.as_slice();
            for j in 0..m {
                for k in 0..m {
                    let s = y[j * m + k];
                    for (o, &x) in out[j * m..(j + 1) * m]
                        .iter_mut()
                        .zip(&c[k * m..(k + 1) * m])
                    {
                        *o += x * s;
                    }
                }
            }
            product.push(t.elapsed().as_secs_f64() * 1e3);
            assert!(out.iter().all(|x| x.is_finite()));
        }
        eig.sort_by(f64::total_cmp);
        product.sort_by(f64::total_cmp);
        println!("{m} | {:.3} | {:.3}", eig[4], product[4]);
    }
}
