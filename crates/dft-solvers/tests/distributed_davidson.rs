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

#[test]
#[ignore = "diagnostic: prints iteration counts"]
fn iteration_counts() {
    for (n, bands) in [(260usize, 9usize), (500, 20), (700, 32)] {
        let (op, diag) = fixture(n);
        let s = solver(bands);
        let serial = s.solve(&op, None, &diag).unwrap();
        let d = simulate(2, |c| {
            s.solve_distributed(c, &op, None, &diag, policy(bands), false, true, 6)
                .unwrap()
        });
        println!(
            "serial hist {:?}\n   dist hist {:?}",
            serial.history, d[0].history
        );
        println!(
            "n={n} bands={bands} serial it={} apps={} restarts={} | P2 it={} apps={} restarts={}",
            serial.iterations,
            serial.operator_applications,
            serial.restarts,
            d[0].iterations,
            d[0].operator_applications,
            d[0].restarts
        );
    }
}

fn spectrum_fixture(n: usize, group: usize) -> (Operator, Vec<f64>) {
    // Orthogonal Q from deterministic Householder reflections; known clustered spectrum.
    let mut q = RealMatrix::zeros(n, n).unwrap();
    for i in 0..n {
        q[(i, i)] = 1.;
    }
    let mut state = 753_u64;
    for _ in 0..24 {
        let v: Vec<_> = (0..n)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                ((state >> 11) as f64) / (1_u64 << 53) as f64 - 0.5
            })
            .collect();
        let scale = 2. / v.iter().map(|x| x * x).sum::<f64>();
        for j in 0..n {
            let d = (0..n).map(|i| v[i] * q[(i, j)]).sum::<f64>() * scale;
            for i in 0..n {
                q[(i, j)] -= v[i] * d;
            }
        }
    }
    let values: Vec<f64> = (0..n)
        .map(|i| 0.05 * ((i / group) as f64).powf(1.3) - 0.3)
        .collect();
    let mut a = RealMatrix::zeros(n, n).unwrap();
    for j in 0..n {
        for i in 0..=j {
            let v: f64 = (0..n).map(|k| q[(i, k)] * values[k] * q[(j, k)]).sum();
            a[(i, j)] = v;
            a[(j, i)] = v;
        }
    }
    let diag: Vec<f64> = (0..n).map(|i| a[(i, i)]).collect();
    (
        Operator {
            a,
            pool: WorkerPool::new(1).unwrap(),
        },
        diag,
    )
}

#[test]
#[ignore = "diagnostic: clustered spectrum, strict prefix like the si32 case"]
fn clustered_strict_prefix_iterations() {
    let (n, bands, strict) = (800usize, 68usize, 65usize);
    for group in [1usize, 4] {
        let (op, diag) = spectrum_fixture(n, group);
        let s = RealDavidson {
            bands,
            tolerance: 1e-9,
            max_iterations: 400,
            max_subspace: 6 * bands,
            seed: 42,
        };
        let pol = BandTolerances {
            strict_bands: strict,
            buffer_tolerance: 1e-5,
        };
        let serial = s.solve_with_band_tolerances(&op, None, &diag, pol).unwrap();
        println!("serial hist {:?}", serial.history);
        for p in [2usize, 4, 8] {
            let d = simulate(p, |c| {
                s.solve_distributed(c, &op, None, &diag, pol, false, true, 6)
                    .unwrap()
            });
            if p == 2 {
                println!("   P2 hist {:?}", d[0].history);
            }
            println!(
                "group={group} serial it={} apps={} | P={p} it={} apps={} restarts={}/{}",
                serial.iterations,
                serial.operator_applications,
                d[0].iterations,
                d[0].operator_applications,
                d[0].restarts,
                serial.restarts
            );
        }
    }
}

#[test]
#[ignore = "diagnostic: iteration-count distributions over fixture shifts"]
fn clustered_distribution() {
    let (n, bands, strict) = (800usize, 68usize, 65usize);
    for group in [3usize, 4, 6] {
        for shift in [0usize, 1, 2] {
            let (mut op, _) = spectrum_fixture(n, group);
            // Shift the potential landscape slightly so each case is a different problem.
            for i in 0..n {
                op.a[(i, i)] += 0.001 * (shift * (i % 7)) as f64;
            }
            let diag: Vec<f64> = (0..n).map(|i| op.a[(i, i)]).collect();
            let s = RealDavidson {
                bands,
                tolerance: 1e-9,
                max_iterations: 400,
                max_subspace: 6 * bands,
                seed: 42,
            };
            let pol = BandTolerances {
                strict_bands: strict,
                buffer_tolerance: 1e-5,
            };
            let serial = s.solve_with_band_tolerances(&op, None, &diag, pol).unwrap();
            let d4 = simulate(4, |c| {
                s.solve_distributed(c, &op, None, &diag, pol, false, true, 6)
                    .unwrap()
            });
            println!(
                "DIST group={group} shift={shift} serial={} P4={}",
                serial.iterations, d4[0].iterations
            );
        }
    }
}
