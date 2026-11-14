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

struct Counting {
    inner: Operator,
    calls: std::sync::atomic::AtomicUsize,
}
impl RealLinearOperator for Counting {
    fn dimension(&self) -> usize {
        self.inner.a.rows()
    }
    fn worker_pool(&self) -> Option<Arc<WorkerPool>> {
        Some(self.inner.pool.clone())
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), SolverError> {
        self.calls
            .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        self.inner.apply(x, y)
    }
}

#[test]
fn operator_columns_are_divided_among_ranks() {
    let (n, bands) = (400usize, 16usize);
    let (inner, diag) = fixture(n);
    let op = Counting {
        inner,
        calls: 0.into(),
    };
    let s = solver(bands);
    let serial = s.solve(&op, None, &diag).unwrap();
    let serial_calls = op.calls.swap(0, std::sync::atomic::Ordering::Relaxed);
    let d = simulate(4, |c| {
        s.solve_distributed(c, &op, None, &diag, policy(bands), false, true, 6)
            .unwrap()
    });
    let total = op.calls.load(std::sync::atomic::Ordering::Relaxed);
    println!(
        "serial apply calls {serial_calls} (reported {}) ; 4 ranks total {total} (reported {})",
        serial.operator_applications, d[0].operator_applications
    );
    assert!(
        total <= serial_calls + serial_calls / 5 + 8,
        "columns are applied redundantly"
    );
}

#[test]
#[ignore = "diagnostic: serial iteration counts of clustered cases (set BAANDI_SECOND_VECTOR)"]
fn serial_clustered_counts() {
    let (n, bands, strict) = (800usize, 68usize, 65usize);
    let mut line = String::new();
    for group in [1usize, 3, 4, 6] {
        for shift in [0usize, 1, 2] {
            let (mut op, _) = spectrum_fixture(n, group);
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
            let r = s.solve_with_band_tolerances(&op, None, &diag, pol).unwrap();
            line += &format!(
                " g{group}s{shift}:{}/{}",
                r.iterations, r.operator_applications
            );
        }
    }
    println!("SERIALCOUNTS{line}");
}

// ---- M2c: Davidson directly on an arbitrary (non-contiguous) row partition ----
use dft_mpi::Collective;
use dft_solvers::distributed_davidson::DistributedRealOperator;

fn owner(r: usize, ranks: usize) -> usize {
    (r / 3 + r % 2 + r / 17) % ranks
}

struct RowsOp<'a> {
    a: &'a RealMatrix,
    ranks: usize,
    rows: Vec<usize>,
    pool: Arc<WorkerPool>,
}
impl<'a> RowsOp<'a> {
    fn new(a: &'a RealMatrix, rank: usize, ranks: usize) -> Self {
        let rows = (0..a.rows()).filter(|&r| owner(r, ranks) == rank).collect();
        Self {
            a,
            ranks,
            rows,
            pool: WorkerPool::new(1).unwrap(),
        }
    }
}
impl DistributedRealOperator for RowsOp<'_> {
    fn global_len(&self) -> usize {
        self.a.rows()
    }
    fn local_rows(&self) -> &[usize] {
        &self.rows
    }
    fn worker_pool(&self) -> Option<Arc<WorkerPool>> {
        Some(self.pool.clone())
    }
    fn apply_block<C: Collective>(
        &mut self,
        comm: &C,
        x: &RealMatrix,
        y: &mut RealMatrix,
    ) -> Result<(), SolverError> {
        let n = self.a.rows();
        let m = x.cols();
        let rows_of =
            |r: usize| -> Vec<usize> { (0..n).filter(|&i| owner(i, self.ranks) == r).collect() };
        let counts: Vec<usize> = (0..self.ranks).map(|r| rows_of(r).len() * m).collect();
        let mut flat = vec![0.0; n * m];
        comm.allgatherv_f64(x.as_slice(), &counts, &mut flat)
            .unwrap();
        let mut full = vec![0.0; n * m];
        let mut at = 0;
        for r in 0..self.ranks {
            let rows_r = rows_of(r);
            let c = rows_r.len();
            for j in 0..m {
                for (k, &g) in rows_r.iter().enumerate() {
                    full[j * n + g] = flat[at + j * c + k];
                }
            }
            at += c * m;
        }
        for j in 0..m {
            for (l, &gr) in self.rows.iter().enumerate() {
                y[(l, j)] = (0..n).map(|k| self.a[(gr, k)] * full[j * n + k]).sum();
            }
        }
        Ok(())
    }
}

fn bits_of(m: &RealMatrix) -> Vec<u64> {
    m.as_slice().iter().map(|x| x.to_bits()).collect()
}

fn check_rows_case(op: &Operator, diag: &[f64], bands: usize, pol: BandTolerances, tol: f64) {
    check_rows_case_tuned(op, diag, bands, pol, tol, DistributedTuning::default());
}

/// Returns the rank-0 distributed solution of every rank count for comparisons.
fn check_rows_case_tuned(
    op: &Operator,
    diag: &[f64],
    bands: usize,
    pol: BandTolerances,
    tol: f64,
    tuning: DistributedTuning,
) -> Vec<dft_solvers::RealEigenSolution> {
    let n = op.dimension();
    let mut out = Vec::new();
    let s = RealDavidson {
        bands,
        tolerance: tol,
        max_iterations: 400,
        max_subspace: 6 * bands,
        seed: 42,
    };
    let serial = s.solve_with_band_tolerances(op, None, diag, pol).unwrap();
    for ranks in [2usize, 3, 4, 5] {
        let run = || {
            simulate(ranks, |comm| {
                let mut rop = RowsOp::new(&op.a, comm.rank(), ranks);
                let dl: Vec<f64> = rop.rows.iter().map(|&r| diag[r]).collect();
                s.solve_distributed_rows_tuned(
                    comm, &mut rop, None, &dl, pol, false, true, 6, tuning,
                )
                .unwrap()
            })
        };
        let results = run();
        for r in &results {
            for (a, b) in r.eigenvalues.iter().zip(&serial.eigenvalues) {
                assert!((a - b).abs() < 1e-9, "P={ranks}: {a} vs {b}");
            }
            for j in 0..pol.strict_bands {
                let col = &r.eigenvectors.as_slice()[j * n..(j + 1) * n];
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
            assert_eq!(bits_of(&r.eigenvectors), bits_of(&results[0].eigenvectors));
            assert_eq!(r.eigenvalues, results[0].eigenvalues);
        }
        let again = run();
        assert_eq!(
            bits_of(&results[0].eigenvectors),
            bits_of(&again[0].eigenvectors),
            "not reproducible at P={ranks}"
        );
        out.push(results.into_iter().next().unwrap());
    }
    out
}

#[test]
fn rows_solver_matches_serial_on_a_scattered_partition() {
    let (op, diag) = fixture(260);
    check_rows_case(&op, &diag, 9, policy(9), 1e-11);
}

#[test]
fn rows_solver_matches_serial_on_clustered_spectra() {
    for group in [1usize, 4] {
        let (op, diag) = spectrum_fixture(64, group);
        let pol = BandTolerances {
            strict_bands: 4,
            buffer_tolerance: 1e-5,
        };
        check_rows_case(&op, &diag, 6, pol, 1e-9);
    }
}

// ---- MPI entry cost (docs/performance/mpi-entry-cost-analysis.md) ----
use dft_solvers::distributed_davidson::DistributedTuning;

fn work_of(r: &dft_solvers::RealEigenSolution) -> (usize, usize) {
    (r.iterations, r.operator_applications)
}

/// M2_LEAN is bit-identical; M2_DGEMM is bit-identical on pure builds (the wide
/// route is the same product there; on MKL builds it is a rounding change and only
/// the accuracy checks of `check_rows_case_tuned` apply).
#[test]
fn lean_and_wide_tuning_keep_the_classic_bits() {
    let (op, diag) = fixture(260);
    let base = check_rows_case_tuned(
        &op,
        &diag,
        9,
        policy(9),
        1e-11,
        DistributedTuning::default(),
    );
    let lean = DistributedTuning {
        lean: true,
        ..Default::default()
    };
    let wide = DistributedTuning {
        wide_gemm: true,
        ..Default::default()
    };
    for (name, tuning, exact) in [
        ("lean", lean, true),
        ("wide", wide, !dft_numerics::MKL_BACKEND),
    ] {
        let got = check_rows_case_tuned(&op, &diag, 9, policy(9), 1e-11, tuning);
        for (a, b) in got.iter().zip(&base) {
            if exact {
                assert_eq!(bits_of(&a.eigenvectors), bits_of(&b.eigenvectors), "{name}");
                assert_eq!(work_of(a), work_of(b), "{name}");
            }
        }
    }
}

/// M2_ORTHO=serial: the serial solver's orthonormalization order (one projection,
/// DGKS, panel Gram-Schmidt) keeps eigenvalues within 1e-9 of the serial solver,
/// residuals below 1e-8, all ranks bit-equal and every rank count reproducible
/// (checked inside `check_rows_case_tuned`), and the Davidson iterations within
/// 10 percent plus two of the classic order, for the panel rules whole block,
/// fixed 4 (several panels) and the shape rule, with and without DGKS.
#[test]
fn serial_order_orthonormalization_matches_serial_eigenpairs() {
    let (a, d) = fixture(260);
    let (b, e) = spectrum_fixture(64, 4);
    let clustered = BandTolerances {
        strict_bands: 4,
        buffer_tolerance: 1e-5,
    };
    let cases = [
        ("scattered", &a, &d, 9, policy(9), 1e-11),
        ("clustered", &b, &e, 6, clustered, 1e-9),
    ];
    for (name, op, diag, bands, pol, tol) in cases {
        let classic =
            check_rows_case_tuned(op, diag, bands, pol, tol, DistributedTuning::default());
        for panel in [None, Some(4), Some(0)] {
            for dgks in [None, Some(0.1)] {
                let tuning = DistributedTuning {
                    serial_ortho: true,
                    lean: true,
                    wide_gemm: true,
                    panel,
                    dgks,
                };
                let got = check_rows_case_tuned(op, diag, bands, pol, tol, tuning);
                for (g, c) in got.iter().zip(&classic) {
                    let (gi, ga) = work_of(g);
                    let (ci, ca) = work_of(c);
                    println!(
                        "{name} panel {panel:?} dgks {dgks:?}: iterations {gi} vs {ci}, applications {ga} vs {ca}"
                    );
                    assert!(
                        gi <= ci + ci / 10 + 2,
                        "{name} {panel:?} {dgks:?}: {gi} vs {ci}"
                    );
                }
            }
        }
    }
}
