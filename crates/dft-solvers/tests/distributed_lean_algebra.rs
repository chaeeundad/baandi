//! MPI-OVERLAP-DAVIDSON (`BAANDI_M2_LEAN_ALGEBRA=on`): the serial lean algebra
//! on the distributed stick-row solver (docs/design/mpi-overlap-davidson.md).
//! Synthetic operators, 1 to 4 simulated ranks with an interleaved row
//! partition:
//! - eigenvalues agree with the orthonormal (panel Gram-Schmidt) path within
//!   the strict class (1e-10 Ha, Ritz-change tolerance 1e-11) on a cubic-like,
//!   an exactly degenerate (two identical blocks) and a skewed spectrum, cold
//!   and warm, capacities 3 and 6 per band (contractions included); returned
//!   vectors orthonormal to 1e-11 and reported residuals equal to the true
//!   residuals;
//! - every rank returns the same bits and a fixed rank count is bit-reproducible;
//! - at most two border tree sums per iteration (S border with the residual
//!   norms before H is applied, H border after), fewer tree sums per solve
//!   than the orthonormal path, no H column applied to a rejected block;
//! - near-dependent corrections (the ill-conditioned fixture of the serial lean
//!   test) take the same rank decisions as the serial lean solver: the same
//!   appends, stage-1 re-projections, stage-2 fallbacks, iterations and H
//!   columns (no wasted application), with both fallback stages exercised (stage 2 down to the random
//!   block on an exactly diagonal operator).
use dft_mpi::{Collective, SimRank, simulate};
use dft_numerics::{pool::WorkerPool, real::RealMatrix};
use dft_solvers::distributed_davidson::{
    DistributedKeptStart, DistributedRealOperator, DistributedStart, DistributedTuning,
    lean_algebra_stats, reduction_stats,
};
use dft_solvers::real_davidson::{RealEigenSolution, ResearchOptions, lean_counters};
use dft_solvers::{
    BandTolerances, DavidsonPreconditioner, RealDavidson, RealLinearOperator, SolverError,
};
use std::sync::{Arc, Mutex};

/// Process-wide counters (reductions, lean statistics): one test at a time.
static SERIAL: Mutex<()> = Mutex::new(());
fn serial() -> std::sync::MutexGuard<'static, ()> {
    SERIAL.lock().unwrap_or_else(|e| e.into_inner())
}

fn symmetric(n: usize, entry: impl Fn(usize, usize) -> f64) -> RealMatrix {
    let mut a = RealMatrix::zeros(n, n).unwrap();
    for j in 0..n {
        for i in 0..=j {
            let v = entry(i, j);
            a[(i, j)] = v;
            a[(j, i)] = v;
        }
    }
    a
}

/// Kinetic shells with a decaying coupling (Ha).
fn cubic(n: usize, shift: f64) -> RealMatrix {
    symmetric(n, |i, j| {
        if i == j {
            0.2 * (i / 2) as f64 + 0.3 + shift * ((i * 7 + 3) as f64).sin()
        } else {
            0.025 * ((i * 31 + j * 17 + 1) as f64).sin() / ((j - i) as f64).sqrt()
        }
    })
}

/// Two identical blocks (rows interleaved): every eigenvalue exactly doubly
/// degenerate.
fn degenerate(n: usize, shift: f64) -> RealMatrix {
    symmetric(n, |i, j| {
        if i % 2 != j % 2 {
            return 0.0;
        }
        let (p, q) = (i / 2, j / 2);
        if p == q {
            0.15 * p as f64 + 0.1 + shift * ((p * 5 + 1) as f64).cos()
        } else {
            0.03 * ((p * 11 + q * 3 + 2) as f64).cos() / (1.0 + (q - p) as f64)
        }
    })
}

/// Skewed-cell plane waves with a smooth local potential.
fn skew(n: usize, shift: f64) -> RealMatrix {
    let b = [[1.0, 0.0, 0.0], [0.34, 0.94, 0.0], [-0.25, 0.3, 1.35]];
    let mut waves = Vec::new();
    for h in -5_i32..=5 {
        for k in -5_i32..=5 {
            for l in -4_i32..=4 {
                let g: [f64; 3] = std::array::from_fn(|c| {
                    f64::from(h) * b[0][c] + f64::from(k) * b[1][c] + f64::from(l) * b[2][c]
                });
                waves.push((0.5 * (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]), g));
            }
        }
    }
    waves.sort_by(|a, b| a.0.total_cmp(&b.0));
    waves.truncate(n);
    symmetric(n, |i, j| {
        let d: [f64; 3] = std::array::from_fn(|c| waves[i].1[c] - waves[j].1[c]);
        let q2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        let local = -(0.35 + shift) * (-0.6 * q2).exp() * (0.7 * d[0] - 0.4 * d[2]).cos();
        if i == j { waves[i].0 + local } else { local }
    })
}

/// The serial lean test's ill-conditioned fixture: an almost diagonal operator
/// whose diagonal is the exact preconditioner (corrections nearly dependent).
fn ill_conditioned(n: usize, coupling: f64) -> RealMatrix {
    symmetric(n, |i, j| {
        if i == j {
            0.3 * i as f64
        } else {
            coupling * ((i * 13 + j * 7) as f64).cos()
        }
    })
}

fn owner(r: usize, ranks: usize) -> usize {
    (r / 3 + r % 2 + r / 17) % ranks
}

struct RowsOp<'a> {
    a: &'a RealMatrix,
    rows: Vec<usize>,
    all_rows: Vec<Vec<usize>>,
    pool: Arc<WorkerPool>,
}
impl<'a> RowsOp<'a> {
    fn new(a: &'a RealMatrix, rank: usize, ranks: usize) -> Self {
        let n = a.rows();
        let all_rows: Vec<Vec<usize>> = (0..ranks)
            .map(|r| (0..n).filter(|&i| owner(i, ranks) == r).collect())
            .collect();
        Self {
            a,
            rows: all_rows[rank].clone(),
            all_rows,
            pool: WorkerPool::new(1).unwrap(),
        }
    }
    fn diagonal(&self) -> Vec<f64> {
        self.rows.iter().map(|&r| self.a[(r, r)]).collect()
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
        let counts: Vec<usize> = self.all_rows.iter().map(|r| r.len() * m).collect();
        let mut flat = vec![0.0; n * m];
        comm.allgatherv_f64(x.as_slice(), &counts, &mut flat)
            .unwrap();
        let mut full = vec![0.0; n * m];
        let mut at = 0;
        for rows_r in &self.all_rows {
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

fn tuning(lean: bool) -> DistributedTuning {
    // The adopted MPI defaults (serial order, dgemm route, slice loops,
    // workspace, selected roots, fused panel Gram-Schmidt) plus the switch.
    DistributedTuning {
        lean_algebra: lean,
        ..DistributedTuning::from_env()
    }
}

/// True residual norms and orthonormality defect of a gathered solution.
fn check(a: &RealMatrix, sol: &RealEigenSolution) -> (Vec<f64>, f64) {
    let n = a.rows();
    let x = &sol.eigenvectors;
    let mut res = Vec::new();
    for (j, &theta) in sol.eigenvalues.iter().enumerate() {
        let r: f64 = (0..n)
            .map(|i| {
                let ax: f64 = (0..n).map(|k| a[(i, k)] * x[(k, j)]).sum();
                (ax - theta * x[(i, j)]).powi(2)
            })
            .sum();
        res.push(r.sqrt());
    }
    let mut defect = 0_f64;
    for i in 0..x.cols() {
        for j in 0..=i {
            let d: f64 = (0..n).map(|r| x[(r, i)] * x[(r, j)]).sum();
            defect = defect.max((d - if i == j { 1. } else { 0. }).abs());
        }
    }
    (res, defect)
}

struct Run {
    cold: RealEigenSolution,
    warm: RealEigenSolution,
}

/// Cold solve of `a`, then a warm solve of `b` from the cold eigenvectors.
fn solve(
    ranks: usize,
    a: &RealMatrix,
    b: &RealMatrix,
    bands: usize,
    per_band: usize,
    lean: bool,
) -> Vec<Run> {
    let solver = RealDavidson {
        bands,
        tolerance: 1e-11,
        max_iterations: 400,
        max_subspace: 2 * bands,
        seed: 42,
    };
    let policy = BandTolerances {
        strict_bands: bands - 1,
        buffer_tolerance: 1e-9,
    };
    let t = tuning(lean);
    simulate(ranks, |sim: &SimRank| {
        let mut op = RowsOp::new(a, sim.rank(), ranks);
        let diag = op.diagonal();
        let cold = solver
            .solve_distributed_rows_tuned(
                sim, &mut op, None, &diag, policy, true, false, per_band, t,
            )
            .unwrap();
        let mut op_b = RowsOp::new(b, sim.rank(), ranks);
        let diag_b = op_b.diagonal();
        let warm = solver
            .solve_distributed_rows_tuned(
                sim,
                &mut op_b,
                Some(&cold.eigenvectors),
                &diag_b,
                policy,
                true,
                false,
                per_band,
                t,
            )
            .unwrap();
        Run { cold, warm }
    })
}

fn bits(s: &RealEigenSolution) -> Vec<u64> {
    s.eigenvalues
        .iter()
        .chain(s.eigenvectors.as_slice())
        .chain(&s.residuals)
        .map(|v| v.to_bits())
        .collect()
}

fn tree_sums() -> u64 {
    reduction_stats().1.iter().map(|s| s.calls).sum()
}

#[test]
fn lean_algebra_agrees_with_orthonormal_path() {
    let _guard = serial();
    type Fixture = fn(usize, f64) -> RealMatrix;
    let fixtures: [(&str, Fixture, usize, usize); 3] = [
        ("cubic", cubic, 96, 6),
        ("degenerate", degenerate, 96, 6),
        ("skew", skew, 120, 8),
    ];
    for (name, f, n, bands) in fixtures {
        let (a, b) = (f(n, 0.0), f(n, 1e-3));
        let oracle_a = dft_numerics::real::symmetric_eigh(&a).unwrap().0;
        let oracle_b = dft_numerics::real::symmetric_eigh(&b).unwrap().0;
        for per_band in [3, 6] {
            for ranks in 1..=4 {
                let sums0 = tree_sums();
                let plain = solve(ranks, &a, &b, bands, per_band, false);
                let sums1 = tree_sums();
                let stats0 = lean_algebra_stats();
                let lean = solve(ranks, &a, &b, bands, per_band, true);
                let stats1 = lean_algebra_stats();
                let sums2 = tree_sums();
                let tag = format!("{name} cap {per_band} ranks {ranks}");
                let hash = |s: &RealEigenSolution| {
                    bits(s).into_iter().fold(0xcbf29ce484222325u64, |h, b| {
                        (h ^ b).wrapping_mul(0x100000001b3)
                    })
                };
                eprintln!(
                    "BASELINE {tag}: {:016x} {:016x} {:?} sums {}",
                    hash(&lean[0].cold),
                    hash(&lean[0].warm),
                    (
                        lean[0].cold.iterations,
                        lean[0].warm.iterations,
                        lean[0].cold.operator_applications,
                        lean[0].warm.operator_applications,
                        stats1.appends - stats0.appends,
                        stats1.reprojections - stats0.reprojections,
                        stats1.fallbacks - stats0.fallbacks
                    ),
                    sums2 - sums1
                );
                for (p, l) in [
                    (&plain[0].cold, &lean[0].cold),
                    (&plain[0].warm, &lean[0].warm),
                ] {
                    for (x, y) in p.eigenvalues.iter().zip(&l.eigenvalues) {
                        assert!((x - y).abs() <= 1e-10, "{tag}: {x} vs {y}");
                    }
                    let matrix = if std::ptr::eq(p, &plain[0].cold) {
                        &a
                    } else {
                        &b
                    };
                    let oracle = if std::ptr::eq(matrix, &a) {
                        &oracle_a
                    } else {
                        &oracle_b
                    };
                    for (j, &value) in l.eigenvalues.iter().enumerate() {
                        assert!(
                            (value - oracle[j]).abs() <= if l.buffer[j] { 1e-9 } else { 1e-10 },
                            "{tag}: dense oracle band {j}: value {value:e}, exact {:e}, delta {:e}",
                            oracle[j],
                            value - oracle[j]
                        );
                    }
                    let (res, defect) = check(matrix, l);
                    // Ritz-change stopping is not residual certification. This
                    // independent 1e-4 Ha bound covers ALL requested roots;
                    // the original 1e-9 Ha buffer Ritz-change policy is kept.
                    // Tighter 1e-6 Ha all-band bounds live in lean_harden and
                    // the forced contraction regression.
                    assert!(
                        res.iter().all(|&r| r <= 1e-4),
                        "{tag}: true residuals {res:?}"
                    );
                    assert!(defect <= 1e-11, "{tag}: defect {defect:e}");
                    for (r, reported) in res.iter().zip(&l.residuals) {
                        assert!(
                            (r - reported).abs() <= 1e-11,
                            "{tag}: {r:e} vs {reported:e}"
                        );
                    }
                    assert!(l.converged.iter().all(|&c| c), "{tag}");
                }
                // Every rank returns the same bits.
                for r in &lean[1..] {
                    assert_eq!(bits(&r.cold), bits(&lean[0].cold), "{tag}");
                    assert_eq!(bits(&r.warm), bits(&lean[0].warm), "{tag}");
                }
                // At most two border reductions per iteration (S border with
                // the residual norms, H border), one more per fallback stage;
                // the rank decisions are checked below.
                let d = |f: fn(&dft_solvers::distributed_davidson::LeanAlgebraStats) -> u64| {
                    f(&stats1) - f(&stats0)
                };
                assert_eq!(d(|s| s.rebuilds), 0, "{tag}");
                assert_eq!(d(|s| s.wasted_columns), 0, "{tag}");
                assert_eq!(d(|s| s.solves), 2 * ranks as u64, "{tag}");
                let iterations = (lean[0].cold.iterations + lean[0].warm.iterations) as u64;
                let extra = d(|s| s.reprojections) + d(|s| s.fallbacks);
                assert!(
                    d(|s| s.border_reductions) <= 2 * iterations * ranks as u64 + extra,
                    "{tag}"
                );
                // Fewer tree sums than the orthonormal path.
                let (plain_sums, lean_sums) = (sums1 - sums0, sums2 - sums1);
                assert!(lean_sums < plain_sums, "{tag}: {lean_sums} vs {plain_sums}");
                eprintln!(
                    "LEAN-ALGEBRA {tag}: iterations {}+{} / {}+{}, H columns {}+{} / {}+{}, \
                     tree sums {lean_sums} / {plain_sums} (lean / orthonormal, all ranks)",
                    lean[0].cold.iterations,
                    lean[0].warm.iterations,
                    plain[0].cold.iterations,
                    plain[0].warm.iterations,
                    lean[0].cold.operator_applications,
                    lean[0].warm.operator_applications,
                    plain[0].cold.operator_applications,
                    plain[0].warm.operator_applications,
                );
            }
            // Cold solve: the rank decisions and counts of the serial lean
            // solver at every rank count.
            let serial = serial_lean(&a, bands, per_band);
            let hash = bits(&serial.0)
                .into_iter()
                .fold(0xcbf29ce484222325u64, |h, b| {
                    (h ^ b).wrapping_mul(0x100000001b3)
                });
            eprintln!(
                "SERIAL-BASELINE {name} cap {per_band}: {hash:016x} {:?}",
                (
                    serial.0.iterations,
                    serial.0.operator_applications,
                    serial.1,
                    serial.2,
                    serial.3
                )
            );
            for ranks in 1..=4 {
                let s0 = lean_algebra_stats();
                let out = cold_lean(ranks, &a, bands, per_band);
                let s1 = lean_algebra_stats();
                let r = ranks as u64;
                let tag = format!("{name} cap {per_band} ranks {ranks} (cold vs serial lean)");
                assert_eq!((s1.appends - s0.appends) / r, serial.1, "{tag}");
                assert_eq!((s1.reprojections - s0.reprojections) / r, serial.2, "{tag}");
                assert_eq!((s1.fallbacks - s0.fallbacks) / r, serial.3, "{tag}");
                assert_eq!(out.iterations, serial.0.iterations, "{tag}");
                assert_eq!(
                    out.operator_applications, serial.0.operator_applications,
                    "{tag}"
                );
                for (x, y) in out.eigenvalues.iter().zip(&serial.0.eigenvalues) {
                    assert!((x - y).abs() <= 1e-10, "{tag}: {x} vs {y}");
                }
                eprintln!(
                    "LEAN-ALGEBRA {tag}: appends {}, reprojections {}, fallbacks {}, wasted H {}",
                    serial.1,
                    serial.2,
                    serial.3,
                    (s1.wasted_columns - s0.wasted_columns) / r
                );
            }
            // A fixed rank count is bit-reproducible.
            let first = solve(3, &a, &b, bands, per_band, true);
            let again = solve(3, &a, &b, bands, per_band, true);
            assert_eq!(bits(&first[0].warm), bits(&again[0].warm), "{name}");
            assert_eq!(bits(&first[0].cold), bits(&again[0].cold), "{name}");
        }
    }
}

fn lean_solver(bands: usize) -> (RealDavidson, BandTolerances) {
    (
        RealDavidson {
            bands,
            tolerance: 1e-11,
            max_iterations: 400,
            max_subspace: 2 * bands,
            seed: 42,
        },
        BandTolerances {
            strict_bands: bands - 1,
            buffer_tolerance: 1e-9,
        },
    )
}

/// Serial lean cold solve: (solution, appends, re-projections, fallbacks).
fn serial_lean(
    a: &RealMatrix,
    bands: usize,
    per_band: usize,
) -> (RealEigenSolution, u64, u64, u64) {
    let (solver, policy) = lean_solver(bands);
    let op = Serial {
        a: a.clone(),
        pool: WorkerPool::new(1).unwrap(),
    };
    let diag: Vec<f64> = (0..a.rows()).map(|i| a[(i, i)]).collect();
    let before = lean_counters();
    let sol = solver
        .solve_research(
            &op,
            None,
            &diag,
            policy,
            true,
            false,
            per_band,
            &mut None,
            DavidsonPreconditioner::Default,
            ResearchOptions {
                lean: true,
                ..ResearchOptions::from_env()
            },
        )
        .unwrap();
    let after = lean_counters();
    (
        sol,
        after.appends - before.appends,
        after.reprojections - before.reprojections,
        after.fallbacks - before.fallbacks,
    )
}

/// Distributed lean cold solve (rank 0's solution).
fn cold_lean(ranks: usize, a: &RealMatrix, bands: usize, per_band: usize) -> RealEigenSolution {
    let (solver, policy) = lean_solver(bands);
    let mut out = simulate(ranks, |sim: &SimRank| {
        let mut op = RowsOp::new(a, sim.rank(), ranks);
        let diag = op.diagonal();
        solver
            .solve_distributed_rows_tuned(
                sim,
                &mut op,
                None,
                &diag,
                policy,
                true,
                false,
                per_band,
                tuning(true),
            )
            .unwrap()
    });
    out.swap_remove(0)
}

/// Near-dependent corrections: the rank decisions (appends, stage-1
/// re-projections, stage-2 fallbacks), iterations and H columns of the serial
/// lean solver, at every rank count; both fallback stages fire.
#[test]
fn rank_decisions_match_the_serial_lean_solver() {
    let _guard = serial();
    let bands = 6;
    let solver = RealDavidson {
        bands,
        tolerance: 1e-12,
        max_iterations: 400,
        max_subspace: 2 * bands,
        seed: 42,
    };
    let policy = BandTolerances {
        strict_bands: 5,
        buffer_tolerance: 1e-10,
    };
    for (coupling, stage2) in [(2e-6, false), (0.0, true)] {
        let a = ill_conditioned(150, coupling);
        let op = Serial {
            a: a.clone(),
            pool: WorkerPool::new(1).unwrap(),
        };
        let diag: Vec<f64> = (0..150).map(|i| a[(i, i)]).collect();
        let before = lean_counters();
        let reference = solver
            .solve_research(
                &op,
                None,
                &diag,
                policy,
                true,
                false,
                6,
                &mut None,
                DavidsonPreconditioner::Default,
                ResearchOptions {
                    lean: true,
                    ..ResearchOptions::from_env()
                },
            )
            .unwrap();
        let after = lean_counters();
        let serial_reproj = after.reprojections - before.reprojections;
        let serial_fallbacks = after.fallbacks - before.fallbacks;
        let serial_appends = after.appends - before.appends;
        assert!(serial_reproj > 0, "{coupling}: serial stage 1");
        if stage2 {
            assert!(serial_fallbacks > 0, "{coupling}: serial stage 2");
        }
        for ranks in 1..=4 {
            let s0 = lean_algebra_stats();
            let out = simulate(ranks, |sim: &SimRank| {
                let mut rows = RowsOp::new(&a, sim.rank(), ranks);
                let d = rows.diagonal();
                let mut x = RealMatrix::zeros(0, 0).unwrap();
                let mut kept = DistributedKeptStart::new();
                solver
                    .solve_distributed_rows_local(
                        sim,
                        &mut rows,
                        None::<DistributedStart<'_>>,
                        &d,
                        policy,
                        true,
                        false,
                        6,
                        tuning(true),
                        &mut x,
                        &mut kept,
                        true,
                    )
                    .unwrap()
            });
            let s1 = lean_algebra_stats();
            let r = ranks as u64;
            let tag = format!("coupling {coupling} ranks {ranks}");
            assert_eq!((s1.appends - s0.appends) / r, serial_appends, "{tag}");
            assert_eq!(
                (s1.reprojections - s0.reprojections) / r,
                serial_reproj,
                "{tag}"
            );
            assert_eq!((s1.fallbacks - s0.fallbacks) / r, serial_fallbacks, "{tag}");
            assert_eq!(out[0].iterations, reference.iterations, "{tag}");
            assert_eq!(
                out[0].operator_applications, reference.operator_applications,
                "{tag}: H columns of the serial lean solver"
            );
            for (x, y) in out[0].eigenvalues.iter().zip(&reference.eigenvalues) {
                assert!((x - y).abs() <= 1e-9, "{tag}: {x} vs {y}");
            }
            let oracle = dft_numerics::real::symmetric_eigh(&a).unwrap().0;
            for (j, &value) in out[0].eigenvalues.iter().enumerate() {
                assert!(
                    (value - oracle[j]).abs() <= 1e-10,
                    "{tag}: dense oracle band {j}: value {value:e}, exact {:e}, delta {:e}",
                    oracle[j],
                    value - oracle[j]
                );
            }
            let (res, defect) = check(&a, &out[0]);
            assert!(
                res.iter().all(|&r| r <= 1e-4),
                "{tag}: true residuals {res:?}"
            );
            assert!(defect <= 1e-11, "{tag}: defect {defect:e}");
            for (t, reported) in res.iter().zip(&out[0].residuals) {
                assert!((t - reported).abs() <= 1e-11, "{tag}");
            }
            for o in &out[1..] {
                assert_eq!(bits(o), bits(&out[0]), "{tag}");
            }
            eprintln!(
                "LEAN-ALGEBRA fallback {tag}: appends {}, reprojections {}, fallbacks {}, \
                 wasted H columns {}, H columns {} (serial {}), min pivot {:e}",
                (s1.appends - s0.appends) / r,
                (s1.reprojections - s0.reprojections) / r,
                (s1.fallbacks - s0.fallbacks) / r,
                (s1.wasted_columns - s0.wasted_columns) / r,
                out[0].operator_applications,
                reference.operator_applications,
                s1.min_pivot
            );
        }
    }
}

/// Strict (residual) solves never take the lean algebra.
#[test]
fn strict_solves_are_unchanged() {
    let _guard = serial();
    let a = cubic(96, 0.0);
    let solver = RealDavidson {
        bands: 6,
        tolerance: 1e-9,
        max_iterations: 400,
        max_subspace: 12,
        seed: 42,
    };
    let policy = BandTolerances {
        strict_bands: 6,
        buffer_tolerance: 1e-9,
    };
    let run = |lean: bool| {
        simulate(2, |sim: &SimRank| {
            let mut op = RowsOp::new(&a, sim.rank(), 2);
            let d = op.diagonal();
            bits(
                &solver
                    .solve_distributed_rows_tuned(
                        sim,
                        &mut op,
                        None,
                        &d,
                        policy,
                        false,
                        true,
                        4,
                        tuning(lean),
                    )
                    .unwrap(),
            )
        })
    };
    let s0 = lean_algebra_stats();
    assert_eq!(run(true), run(false));
    assert_eq!(lean_algebra_stats().solves, s0.solves);
}

struct Serial {
    a: RealMatrix,
    pool: Arc<WorkerPool>,
}
impl RealLinearOperator for Serial {
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
