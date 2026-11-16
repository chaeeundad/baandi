//! F05 (`BAANDI_M2_SELECTED_ROOTS=on`): selected-root rotation in the
//! distributed Davidson for Ritz-change (standard) solves. Against the
//! full-rotation path on synthetic problems (scattered, clustered and
//! degenerate spectra; default and contracting capacities; 1 to 4 simulated
//! ranks): identical iteration, H-column, restart and contraction counts,
//! identical convergence flags and root ordering, eigenvalues and every band's
//! returned residual within rounding, every returned residual equal to an
//! independent residual with the full matrix; fewer rotation columns, GEMM
//! flops and reduction payload. Strict (residual) solves are unchanged bit for bit.
use dft_mpi::{Collective, MpiError, SimRank, simulate};
use dft_numerics::{pool::WorkerPool, real::RealMatrix};
use dft_solvers::distributed_davidson::{
    DistributedRealOperator, DistributedTuning, rotation_stats,
};
use dft_solvers::{BandTolerances, RealDavidson, RealEigenSolution, SolverError};
use std::cell::Cell;
use std::sync::Arc;

/// Counts the tree sums (calls and payload doubles) of one rank.
struct Counting<'a> {
    inner: &'a SimRank,
    sums: Cell<u64>,
    payload: Cell<u64>,
}
impl Collective for Counting<'_> {
    fn rank(&self) -> usize {
        self.inner.rank()
    }
    fn size(&self) -> usize {
        self.inner.size()
    }
    fn barrier(&self) -> Result<(), MpiError> {
        self.inner.barrier()
    }
    fn allgatherv_f64(
        &self,
        local: &[f64],
        counts: &[usize],
        out: &mut [f64],
    ) -> Result<(), MpiError> {
        self.inner.allgatherv_f64(local, counts, out)
    }
    fn tree_sum_f64(&self, local: &[f64], out: &mut [f64]) -> Result<(), MpiError> {
        self.sums.set(self.sums.get() + 1);
        self.payload.set(self.payload.get() + local.len() as u64);
        self.inner.tree_sum_f64(local, out)
    }
    fn bcast_f64(&self, data: &mut [f64], root: usize) -> Result<(), MpiError> {
        self.inner.bcast_f64(data, root)
    }
    fn alltoall_f64(&self, send: &[f64], block: usize, recv: &mut [f64]) -> Result<(), MpiError> {
        self.inner.alltoall_f64(send, block, recv)
    }
}

fn fixture(n: usize) -> (RealMatrix, Vec<f64>) {
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
    (a, kinetic)
}

/// Spectrum in groups of `group` equal eigenvalues (exact degeneracy) in a
/// rotated basis.
fn spectrum_fixture(n: usize, group: usize) -> (RealMatrix, Vec<f64>) {
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
    (a, diag)
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

fn entry_stack() -> DistributedTuning {
    DistributedTuning {
        serial_ortho: true,
        wide_gemm: true,
        lean: true,
        panel: Some(4),
        dgks: Some(0.1),
        workspace: true,
        ..Default::default()
    }
}

struct Run {
    solution: RealEigenSolution,
    tree_sums: u64,
    payload: u64,
}

#[allow(clippy::too_many_arguments)]
fn run(
    a: &RealMatrix,
    diag: &[f64],
    bands: usize,
    policy: BandTolerances,
    tol: f64,
    ritz: bool,
    capacity: usize,
    ranks: usize,
    tuning: DistributedTuning,
) -> Vec<Run> {
    let s = RealDavidson {
        bands,
        tolerance: tol,
        max_iterations: 400,
        max_subspace: capacity * bands,
        seed: 42,
    };
    simulate(ranks, |sim| {
        let comm = Counting {
            inner: sim,
            sums: Cell::new(0),
            payload: Cell::new(0),
        };
        let mut op = RowsOp::new(a, sim.rank(), ranks);
        let dl: Vec<f64> = op.rows.iter().map(|&r| diag[r]).collect();
        let solution = s
            .solve_distributed_rows_tuned(
                &comm, &mut op, None, &dl, policy, ritz, !ritz, capacity, tuning,
            )
            .unwrap();
        Run {
            solution,
            tree_sums: comm.sums.get(),
            payload: comm.payload.get(),
        }
    })
}

fn true_residual(a: &RealMatrix, x: &[f64], lambda: f64) -> f64 {
    let n = a.rows();
    (0..n)
        .map(|i| {
            let ax: f64 = (0..n).map(|k| a[(i, k)] * x[k]).sum();
            (ax - lambda * x[i]).powi(2)
        })
        .sum::<f64>()
        .sqrt()
}

struct Problem {
    name: &'static str,
    a: RealMatrix,
    diag: Vec<f64>,
    bands: usize,
    policy: BandTolerances,
    tol: f64,
}

fn problems() -> Vec<Problem> {
    let (a, d) = fixture(200);
    let (b, e) = spectrum_fixture(72, 1);
    let (c, f) = spectrum_fixture(72, 4);
    let (g, h) = spectrum_fixture(80, 3);
    vec![
        Problem {
            name: "scattered",
            a,
            diag: d,
            bands: 9,
            policy: BandTolerances {
                strict_bands: 8,
                buffer_tolerance: 1e-6,
            },
            tol: 1e-10,
        },
        Problem {
            name: "clustered",
            a: b,
            diag: e,
            bands: 7,
            policy: BandTolerances {
                strict_bands: 6,
                buffer_tolerance: 1e-6,
            },
            tol: 1e-9,
        },
        Problem {
            name: "degenerate4",
            a: c,
            diag: f,
            bands: 8,
            policy: BandTolerances {
                strict_bands: 8,
                buffer_tolerance: 1e-9,
            },
            tol: 1e-9,
        },
        Problem {
            name: "degenerate3-split",
            a: g,
            diag: h,
            bands: 7,
            policy: BandTolerances {
                strict_bands: 7,
                buffer_tolerance: 1e-9,
            },
            tol: 1e-9,
        },
    ]
}

/// `rotation_stats` is process-wide: tests that read it run one at a time.
static SERIAL: std::sync::Mutex<()> = std::sync::Mutex::new(());

fn compare(ranks_list: &[usize]) {
    let _guard = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    let full = entry_stack();
    let selected = DistributedTuning {
        selected_roots: true,
        ..entry_stack()
    };
    for p in problems() {
        for capacity in [6usize, 3] {
            for &ranks in ranks_list {
                let r0 = rotation_stats();
                let off = run(
                    &p.a, &p.diag, p.bands, p.policy, p.tol, true, capacity, ranks, full,
                );
                let r1 = rotation_stats();
                let on = run(
                    &p.a, &p.diag, p.bands, p.policy, p.tol, true, capacity, ranks, selected,
                );
                let r2 = rotation_stats();
                let (a, b) = (&off[0].solution, &on[0].solution);
                let tag = format!("{} capacity {capacity} P={ranks}", p.name);
                // Work and control flow.
                assert_eq!(a.iterations, b.iterations, "{tag}: iterations");
                assert_eq!(
                    a.operator_applications, b.operator_applications,
                    "{tag}: H columns"
                );
                assert_eq!(a.restarts, b.restarts, "{tag}: restarts");
                assert_eq!(
                    a.contracted_restarts, b.contracted_restarts,
                    "{tag}: contractions"
                );
                assert_eq!(a.max_basis_dimension, b.max_basis_dimension, "{tag}");
                assert_eq!(a.converged, b.converged, "{tag}: convergence flags");
                assert!(!b.residuals_certified);
                // Root ordering and degeneracy: ascending values, same order,
                // equal within rounding.
                for w in b.eigenvalues.windows(2) {
                    assert!(w[0] <= w[1], "{tag}: root ordering");
                }
                let mut exact = true;
                for (j, (x, y)) in a.eigenvalues.iter().zip(&b.eigenvalues).enumerate() {
                    exact &= x.to_bits() == y.to_bits();
                    assert!(
                        (x - y).abs() <= 1e-12 * x.abs().max(1.0),
                        "{tag}: band {j} {x} vs {y}"
                    );
                }
                // Every band's returned residual: rounding agreement with the
                // full path and with an independent residual from the matrix.
                let n = p.a.rows();
                for j in 0..p.bands {
                    let (ra, rb) = (a.residuals[j], b.residuals[j]);
                    assert!(
                        (ra - rb).abs() <= 1e-10 * ra.max(1e-6),
                        "{tag}: band {j} residual {ra:e} vs {rb:e}"
                    );
                    let x = &b.eigenvectors.as_slice()[j * n..(j + 1) * n];
                    let t = true_residual(&p.a, x, b.eigenvalues[j]);
                    assert!(
                        (t - rb).abs() <= 1e-9 * t.max(1e-6) + 1e-12,
                        "{tag}: band {j} returned residual {rb:e}, true {t:e}"
                    );
                    exact &= ra.to_bits() == rb.to_bits();
                }
                // Ranks agree bit for bit.
                for r in &on {
                    assert_eq!(r.solution.eigenvalues, b.eigenvalues, "{tag}: ranks");
                }
                let (cols_off, cols_on) = (r1.columns - r0.columns, r2.columns - r1.columns);
                let (fl_off, fl_on) = (r1.flops - r0.flops, r2.flops - r1.flops);
                println!(
                    "{tag}: iterations {} H {} contractions {}; rotation columns {} -> {} (summed over ranks), flops {:.3e} -> {:.3e}; tree sums per rank {} -> {}, payload doubles {} -> {}; bit-identical {exact}",
                    b.iterations,
                    b.operator_applications,
                    b.contracted_restarts,
                    cols_off,
                    cols_on,
                    fl_off as f64,
                    fl_on as f64,
                    off[0].tree_sums,
                    on[0].tree_sums,
                    off[0].payload,
                    on[0].payload
                );
                assert!(cols_on <= cols_off && fl_on <= fl_off, "{tag}: more work");
                assert_eq!(off[0].tree_sums, on[0].tree_sums, "{tag}: reduction count");
                assert!(on[0].payload <= off[0].payload, "{tag}: payload");
            }
        }
    }
}

/// Ranks 2 and 3 (fast subset for `scripts/check.sh`).
#[test]
fn selected_roots_match_the_full_rotation_path() {
    compare(&[2, 3]);
}

/// Ranks 1 to 4.
#[test]
#[ignore = "about 20 s in debug; release: cargo test --release -p dft-solvers --test distributed_selected_roots -- --include-ignored --nocapture"]
fn selected_roots_match_the_full_rotation_path_at_every_rank_count() {
    compare(&[1, 2, 3, 4]);
}

/// Strict (residual) solves ignore the switch: bit-identical results.
#[test]
fn strict_solves_are_unchanged() {
    let _guard = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    let selected = DistributedTuning {
        selected_roots: true,
        ..entry_stack()
    };
    for p in problems().into_iter().take(2) {
        for ranks in [2usize, 3] {
            let a = run(
                &p.a,
                &p.diag,
                p.bands,
                p.policy,
                p.tol,
                false,
                6,
                ranks,
                entry_stack(),
            );
            let b = run(
                &p.a, &p.diag, p.bands, p.policy, p.tol, false, 6, ranks, selected,
            );
            let bits =
                |m: &RealMatrix| m.as_slice().iter().map(|x| x.to_bits()).collect::<Vec<_>>();
            assert_eq!(
                bits(&a[0].solution.eigenvectors),
                bits(&b[0].solution.eigenvectors),
                "{}",
                p.name
            );
            assert_eq!(a[0].solution.history, b[0].solution.history);
            assert_eq!(a[0].tree_sums, b[0].tree_sums);
        }
    }
}

/// F08: the per-site reduction counters see exactly the tree sums the
/// communicator saw (count-only; order and results unchanged).
#[test]
fn reduction_counters_match_the_communicator() {
    use dft_solvers::distributed_davidson::reduction_stats;
    let _guard = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    let p = &problems()[1];
    for ranks in [2usize, 3] {
        let (s0, before) = reduction_stats();
        let runs = run(
            &p.a,
            &p.diag,
            p.bands,
            p.policy,
            p.tol,
            true,
            6,
            ranks,
            entry_stack(),
        );
        let (s1, after) = reduction_stats();
        assert_eq!(s1 - s0, ranks as u64);
        let calls: u64 = after
            .iter()
            .zip(&before)
            .map(|(a, b)| a.calls - b.calls)
            .sum();
        let payload: u64 = after
            .iter()
            .zip(&before)
            .map(|(a, b)| a.payload - b.payload)
            .sum();
        let seen: u64 = runs.iter().map(|r| r.tree_sums).sum();
        let seen_payload: u64 = runs.iter().map(|r| r.payload).sum();
        assert_eq!(calls, seen);
        assert_eq!(payload, seen_payload);
        for (a, b) in after.iter().zip(&before) {
            println!(
                "P={ranks} {}: calls {} payload {}",
                a.site,
                (a.calls - b.calls) / ranks as u64,
                (a.payload - b.payload) / ranks as u64
            );
        }
        let site = |n: &str| {
            let i = after.iter().position(|s| s.site == n).unwrap();
            after[i].calls - before[i].calls
        };
        assert!(site("panel_dot") > 0 && site("panel_norm") > 0 && site("reduced_gram") > 0);
    }
}
