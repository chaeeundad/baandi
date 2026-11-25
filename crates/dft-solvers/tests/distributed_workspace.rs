//! Frozen legacy workspace: repeatability across rank counts, cold and warm
//! starts, and consecutive solves of different shapes.
use dft_mpi::{Collective, simulate};
use dft_numerics::{pool::WorkerPool, real::RealMatrix};
use dft_solvers::distributed_davidson::{
    DistributedRealOperator, DistributedTuning, release_workspace, workspace_stats,
};
use dft_solvers::{BandTolerances, RealDavidson, RealEigenSolution, SolverError};
use std::sync::Arc;

// ---- fixtures ----
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

/// Clustered spectrum (groups of `group` equal eigenvalues) in a rotated basis.
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

/// Dense operator on a scattered row partition; the gather scratch is kept, so
/// a steady-state apply allocates nothing.
struct RowsOp<'a> {
    a: &'a RealMatrix,
    ranks: usize,
    rows: Vec<usize>,
    all_rows: Vec<Vec<usize>>,
    pool: Arc<WorkerPool>,
    counts: Vec<usize>,
    flat: Vec<f64>,
    full: Vec<f64>,
}
impl<'a> RowsOp<'a> {
    fn new(a: &'a RealMatrix, rank: usize, ranks: usize) -> Self {
        let n = a.rows();
        let all_rows: Vec<Vec<usize>> = (0..ranks)
            .map(|r| (0..n).filter(|&i| owner(i, ranks) == r).collect())
            .collect();
        Self {
            a,
            ranks,
            rows: all_rows[rank].clone(),
            all_rows,
            pool: WorkerPool::new(1).unwrap(),
            counts: Vec::with_capacity(ranks),
            flat: vec![0.0; n * n],
            full: vec![0.0; n * n],
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
        self.counts.clear();
        self.counts
            .extend(self.all_rows.iter().map(|rows| rows.len() * m));
        comm.allgatherv_f64(x.as_slice(), &self.counts, &mut self.flat[..n * m])
            .unwrap();
        let mut at = 0;
        for r in 0..self.ranks {
            let rows_r = &self.all_rows[r];
            let c = rows_r.len();
            for j in 0..m {
                for (k, &g) in rows_r.iter().enumerate() {
                    self.full[j * n + g] = self.flat[at + j * c + k];
                }
            }
            at += c * m;
        }
        for j in 0..m {
            for (l, &gr) in self.rows.iter().enumerate() {
                y[(l, j)] = (0..n).map(|k| self.a[(gr, k)] * self.full[j * n + k]).sum();
            }
        }
        Ok(())
    }
}

fn entry_stack() -> DistributedTuning {
    DistributedTuning {
        wide_gemm: true,

        panel: Some(4),
        dgks: Some(0.1),
    }
}

struct Case {
    a: RealMatrix,
    diag: Vec<f64>,
    bands: usize,
    policy: BandTolerances,
    tol: f64,
    ritz: bool,
}

fn cases() -> Vec<Case> {
    let (a, d) = fixture(160);
    let (b, e) = spectrum_fixture(64, 4);
    let clustered = BandTolerances {
        strict_bands: 4,
        buffer_tolerance: 1e-5,
    };
    let strict = BandTolerances {
        strict_bands: 7,
        buffer_tolerance: 1e-11,
    };
    vec![
        Case {
            a: a.clone(),
            diag: d.clone(),
            bands: 7,
            policy: strict,
            tol: 1e-11,
            ritz: false,
        },
        Case {
            a,
            diag: d,
            bands: 7,
            policy: BandTolerances {
                strict_bands: 6,
                buffer_tolerance: 1e-6,
            },
            tol: 1e-10,
            ritz: true,
        },
        Case {
            a: b,
            diag: e,
            bands: 6,
            policy: clustered,
            tol: 1e-9,
            ritz: false,
        },
    ]
}

fn solver(c: &Case, capacity: usize) -> RealDavidson {
    RealDavidson {
        bands: c.bands,
        tolerance: c.tol,
        max_iterations: 400,
        max_subspace: capacity * c.bands,
        seed: 42,
    }
}

/// Everything a solve reports, as bits.
fn fingerprint(r: &RealEigenSolution) -> Vec<u64> {
    let mut out: Vec<u64> = r.eigenvalues.iter().map(|x| x.to_bits()).collect();
    out.extend(r.eigenvectors.as_slice().iter().map(|x| x.to_bits()));
    out.extend(r.residuals.iter().map(|x| x.to_bits()));
    out.extend(r.history.iter().map(|x| x.to_bits()));
    out.extend(r.converged.iter().map(|&b| b as u64));
    out.extend([
        r.iterations as u64,
        r.operator_applications as u64,
        r.restarts as u64,
        r.contracted_restarts as u64,
        r.certification_attempts as u64,
        r.max_basis_dimension as u64,
        r.residuals_certified as u64,
    ]);
    out
}

/// Every case, cold and then warm-started from a perturbed result, one after
/// another on one rank thread (one workspace for all of them).
fn sequence(comm: &dft_mpi::SimRank, ranks: usize, tuning: DistributedTuning) -> Vec<Vec<u64>> {
    let mut out = Vec::new();
    for case in cases() {
        // Capacity 3 forces contraction restarts on these sizes; 6 is the default.
        for capacity in [6usize, 3] {
            let s = solver(&case, capacity);
            let mut op = RowsOp::new(&case.a, comm.rank(), ranks);
            let dl: Vec<f64> = op.rows.iter().map(|&r| case.diag[r]).collect();
            let cold = s
                .solve_distributed_rows_tuned(
                    comm,
                    &mut op,
                    None,
                    &dl,
                    case.policy,
                    case.ritz,
                    !case.ritz,
                    capacity,
                    tuning,
                )
                .unwrap();
            let mut warm = cold.eigenvectors.clone();
            for (i, v) in warm.as_mut_slice().iter_mut().enumerate() {
                *v += 1e-3 * ((i * 7 + 3) as f64).sin();
            }
            let again = s
                .solve_distributed_rows_tuned(
                    comm,
                    &mut op,
                    Some(&warm),
                    &dl,
                    case.policy,
                    case.ritz,
                    !case.ritz,
                    capacity,
                    tuning,
                )
                .unwrap();
            out.push(fingerprint(&cold));
            out.push(fingerprint(&again));
        }
    }
    release_workspace();
    out
}

fn check_identity(stacks: &[(&str, DistributedTuning)], rank_counts: &[usize]) {
    for &(name, base) in stacks {
        for &ranks in rank_counts {
            let first = simulate(ranks, |c| sequence(c, ranks, DistributedTuning { ..base }));
            let again = simulate(ranks, |c| sequence(c, ranks, DistributedTuning { ..base }));
            for (rank, (a, b)) in first.iter().zip(&again).enumerate() {
                assert_eq!(a.len(), b.len());
                for (k, (x, y)) in a.iter().zip(b).enumerate() {
                    assert!(x == y, "{name}: P={ranks} rank {rank} solve {k} differs");
                }
            }
            // Ranks agree with rank 0.
            for b in &again {
                assert_eq!(b, &again[0], "{name}: P={ranks} ranks disagree");
            }
        }
    }
}

/// The adopted entry-cost stack (default on every build), ranks 1 to 4.
#[test]
fn workspace_is_bit_identical_on_the_entry_stack() {
    let before = workspace_stats();
    check_identity(&[("entry", entry_stack())], &[1, 2, 3, 4]);
    let after = workspace_stats();
    assert!(after.solves > before.solves);
    assert!(
        after.warm_solves > before.warm_solves,
        "repeated solves of one shape must reuse the workspace"
    );
}

/// Remaining shared Gram-Schmidt policies, ranks 1 to 4.
#[test]
#[ignore = "about 30 s in debug; release: cargo test --release -p dft-solvers --test distributed_workspace -- --include-ignored"]
fn workspace_is_reproducible_with_shared_orthogonalization_policies() {
    let stacks = [
        ("default", DistributedTuning::default()),
        (
            "serial whole-block two-pass",
            DistributedTuning {
                panel: None,
                dgks: None,
                ..entry_stack()
            },
        ),
    ];
    check_identity(&stacks, &[1, 2, 3, 4]);
}
