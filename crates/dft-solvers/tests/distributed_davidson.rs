//! M1c: distributed real Davidson on simulated ranks against the serial solver.
use dft_mpi::simulate;
use dft_numerics::{pool::WorkerPool, real::RealMatrix};
use dft_solvers::{BandTolerances, RealDavidson, RealLinearOperator, SolverError};
use std::sync::Arc;

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
        for (i, y) in y.iter_mut().enumerate() {
            *y = x.iter().enumerate().map(|(j, x)| self.a[(i, j)] * x).sum();
        }
        Ok(())
    }
}

fn fixture(n: usize) -> (Operator, Vec<f64>) {
    let kinetic: Vec<f64> = (0..n).map(|i| 0.2 * (i / 2) as f64).collect();
    let mut a = RealMatrix::zeros(n, n).unwrap();
    for j in 0..n {
        for i in 0..=j {
            let v = if i == j {
                kinetic[i] + 0.3
            } else {
                0.025 * ((i * 31 + j * 17 + 1) as f64).sin() / ((j - i) as f64).sqrt()
            };
            a[(i, j)] = v;
            a[(j, i)] = v;
        }
    }
    (
        Operator {
            a,
            pool: WorkerPool::new(1).unwrap(),
        },
        kinetic,
    )
}
fn solver(bands: usize) -> RealDavidson {
    RealDavidson {
        bands,
        tolerance: 1e-11,
        max_iterations: 400,
        max_subspace: 6 * bands,
        seed: 42,
    }
}
fn policy(bands: usize) -> BandTolerances {
    BandTolerances {
        strict_bands: bands,
        buffer_tolerance: 1e-11,
    }
}

#[test]
fn distributed_matches_serial_eigenvalues_and_residuals() {
    let (n, bands) = (260usize, 9usize);
    let (op, diag) = fixture(n);
    let s = solver(bands);
    let serial = s.solve(&op, None, &diag).unwrap();
    for ranks in [2usize, 3, 4, 5] {
        let results = simulate(ranks, |comm| {
            s.solve_distributed(comm, &op, None, &diag, policy(bands), false, true, 6)
                .unwrap()
        });
        for r in &results {
            for (a, b) in r.eigenvalues.iter().zip(&serial.eigenvalues) {
                assert!((a - b).abs() < 1e-9, "P={ranks}: {a} vs {b}");
            }
            // Independent residual check with the full operator.
            let x = &r.eigenvectors;
            for j in 0..bands {
                let col = &x.as_slice()[j * n..(j + 1) * n];
                let mut y = vec![0.0; n];
                op.apply(col, &mut y).unwrap();
                let res: f64 = y
                    .iter()
                    .zip(col)
                    .map(|(y, c)| (y - r.eigenvalues[j] * c).powi(2))
                    .sum::<f64>()
                    .sqrt();
                assert!(res < 1e-8, "P={ranks} band {j} residual {res}");
            }
            assert_eq!(
                r.eigenvalues, results[0].eigenvalues,
                "ranks disagree at P={ranks}"
            );
        }
    }
}

#[test]
fn distributed_is_reproducible_at_fixed_rank_count() {
    let (n, bands) = (200usize, 6usize);
    let (op, diag) = fixture(n);
    let s = solver(bands);
    let run = || {
        simulate(3, |comm| {
            s.solve_distributed(comm, &op, None, &diag, policy(bands), false, true, 6)
                .unwrap()
        })
    };
    let (a, b) = (run(), run());
    let bits = |m: &RealMatrix| m.as_slice().iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    assert_eq!(bits(&a[0].eigenvectors), bits(&b[0].eigenvectors));
    assert_eq!(a[0].iterations, b[0].iterations);
}
