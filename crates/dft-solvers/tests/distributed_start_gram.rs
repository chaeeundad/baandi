//! Frozen Gram start: strict solves decline nonorthonormal starts; lean
//! Ritz-change solves factor nonsingular starts. Rank repeatability is preserved.
use dft_mpi::{Collective, SimRank, simulate};
use dft_numerics::{pool::WorkerPool, real::RealMatrix};
use dft_solvers::distributed_davidson::{
    DistributedRealOperator, DistributedTuning, start_gram_stats,
};
use dft_solvers::{BandTolerances, RealDavidson, SolverError};
use std::sync::Arc;

fn matrix(n: usize, shift: &[f64]) -> RealMatrix {
    let mut a = RealMatrix::zeros(n, n).unwrap();
    for j in 0..n {
        for i in 0..=j {
            let v = if i == j {
                0.2 * (i / 2) as f64 + 0.3 + shift[i]
            } else {
                0.025 * ((i * 31 + j * 17 + 1) as f64).sin() / ((j - i) as f64).sqrt()
            };
            a[(i, j)] = v;
            a[(j, i)] = v;
        }
    }
    a
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

fn solve(ranks: usize, scale: f64, ritz: bool) -> Vec<(Vec<u64>, usize)> {
    let n = 96;
    let bands = 6;
    let a = matrix(n, &vec![0.0; n]);
    let shift: Vec<f64> = (0..n).map(|i| 1e-3 * ((i * 7 + 3) as f64).sin()).collect();
    let b = matrix(n, &shift);
    let solver = RealDavidson {
        bands,
        tolerance: 1e-11,
        max_iterations: 400,
        max_subspace: 4 * bands,
        seed: 42,
    };
    let tuning = DistributedTuning {
        wide_gemm: true,

        ..Default::default()
    };
    let policy = BandTolerances {
        strict_bands: bands,
        buffer_tolerance: 1e-11,
    };
    simulate(ranks, |sim: &SimRank| {
        let mut op = RowsOp::new(&a, sim.rank(), ranks);
        let diag: Vec<f64> = op.rows.iter().map(|&r| a[(r, r)]).collect();
        let cold = solver
            .solve_distributed_rows_tuned(sim, &mut op, None, &diag, policy, ritz, !ritz, 4, tuning)
            .unwrap();
        let mut warm = cold.eigenvectors.clone();
        for v in warm.as_mut_slice() {
            *v *= scale;
        }
        let mut op_b = RowsOp::new(&b, sim.rank(), ranks);
        let diag_b: Vec<f64> = op_b.rows.iter().map(|&r| b[(r, r)]).collect();
        let next = solver
            .solve_distributed_rows_tuned(
                sim,
                &mut op_b,
                Some(&warm),
                &diag_b,
                policy,
                ritz,
                !ritz,
                4,
                tuning,
            )
            .unwrap();
        let bits = cold
            .eigenvalues
            .iter()
            .chain(&next.eigenvalues)
            .map(|v| v.to_bits())
            .collect();
        (bits, next.operator_applications)
    })
}

#[test]
fn gram_start_takes_orthonormal_starts_and_keeps_the_eigenvalues() {
    for ranks in 1..=4 {
        for ritz in [false, true] {
            let plain = solve(ranks, 1.0, ritz);
            let before = start_gram_stats();
            let gram = solve(ranks, 1.0, ritz);
            let after = start_gram_stats();
            // Cold (free-electron unit vectors) and warm (returned Ritz vectors)
            // starts are orthonormal: both taken on every rank. Other tests of this
            // binary do not run concurrently, so the counters move only here.
            assert_eq!(after.0 - before.0, 2 * ranks as u64, "ranks {ranks}");
            assert_eq!(after.1, before.1, "ranks {ranks}");
            for (p, g) in plain.iter().zip(&gram) {
                for (x, y) in p.0.iter().zip(&g.0) {
                    let d = (f64::from_bits(*x) - f64::from_bits(*y)).abs();
                    assert!(d <= 1e-10, "ranks {ranks}: {d:e}");
                }
            }
            // Every rank returns the same bits; a fixed rank count is reproducible.
            assert!(gram.iter().all(|r| r.0 == gram[0].0));
            assert_eq!(solve(ranks, 1.0, ritz)[0].0, gram[0].0);
            // Strict declines the scaled start; lean factors its nonsingular Gram.
            let before = start_gram_stats();
            let scaled = solve(ranks, 1.0 + 1e-6, ritz);
            let after = start_gram_stats();
            assert_eq!(
                after.1 - before.1,
                if ritz { 0 } else { ranks as u64 },
                "ranks {ranks}"
            );
            assert_eq!(
                after.0 - before.0,
                if ritz { 2 * ranks as u64 } else { ranks as u64 },
                "ranks {ranks}"
            );
            for (x, y) in scaled[0].0.iter().zip(&plain[0].0) {
                assert!((f64::from_bits(*x) - f64::from_bits(*y)).abs() <= 1e-10);
            }
        }
    }
}
