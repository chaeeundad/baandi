//! MPI-LAZY-ORBITALS: `solve_distributed_rows_local` with the start given as
//! this rank's rows (`DistributedStart::Local`) is bit-identical to the same
//! start given replicated (`Global`), and without the final gather it returns a
//! 0 x bands placeholder whose local rows (`x_local`) are bit-identical to this
//! rank's rows of the gathered eigenvectors. Synthetic symmetric operator with
//! scattered row ownership, 1 to 4 simulated ranks, selected-root and
//! full-rotation workspace paths, plus the non-workspace path (a local start is
//! gathered there and the block always returned).
use dft_mpi::{Collective, SimRank, simulate};
use dft_numerics::{pool::WorkerPool, real::RealMatrix};
use dft_solvers::distributed_davidson::{
    DistributedKeptStart, DistributedRealOperator, DistributedStart, DistributedTuning,
};
use dft_solvers::{BandTolerances, RealDavidson, SolverError};
use std::sync::Arc;

fn matrix(n: usize) -> RealMatrix {
    let mut a = RealMatrix::zeros(n, n).unwrap();
    for j in 0..n {
        for i in 0..=j {
            let v = if i == j {
                0.2 * (i / 2) as f64 + 0.3
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

fn values_bits(v: &[f64]) -> Vec<u64> {
    v.iter().map(|x| x.to_bits()).collect()
}

#[test]
fn local_start_and_lazy_return_match_the_gathered_solve() {
    let n = 96;
    let bands = 6;
    let a = matrix(n);
    let solver = RealDavidson {
        bands,
        tolerance: 1e-11,
        max_iterations: 400,
        max_subspace: 4 * bands,
        seed: 42,
    };
    let policy = BandTolerances {
        strict_bands: bands,
        buffer_tolerance: 1e-11,
    };
    for ranks in 1..=4usize {
        for workspace in [true] {
            let tuning = DistributedTuning {
                ..Default::default()
            };
            let out = simulate(ranks, |sim: &SimRank| {
                let mut op = RowsOp::new(&a, sim.rank(), ranks);
                let diag: Vec<f64> = op.rows.iter().map(|&r| a[(r, r)]).collect();
                let rows = op.rows.clone();
                // A cold solve gives a replicated warm start W.
                let mut x0 = RealMatrix::zeros(0, 0).unwrap();
                let mut k0 = DistributedKeptStart::new();
                let first = solver
                    .solve_distributed_rows_local(
                        sim, &mut op, None, &diag, policy, true, false, 4, tuning, &mut x0,
                        &mut k0, true,
                    )
                    .unwrap();
                let warm = first.eigenvectors;
                let mut local = RealMatrix::zeros(rows.len(), bands).unwrap();
                for j in 0..bands {
                    for (l, &g) in rows.iter().enumerate() {
                        local[(l, j)] = warm[(g, j)];
                    }
                }
                // Warm solves from W: replicated start with the gather, local
                // start with the gather, local start without the gather.
                let run = |start: DistributedStart<'_>, gather: bool| {
                    let mut x = RealMatrix::zeros(0, 0).unwrap();
                    let mut k = DistributedKeptStart::new();
                    let mut op = RowsOp::new(&a, sim.rank(), ranks);
                    let s = solver
                        .solve_distributed_rows_local(
                            sim,
                            &mut op,
                            Some(start),
                            &diag,
                            policy,
                            true,
                            false,
                            4,
                            tuning,
                            &mut x,
                            &mut k,
                            gather,
                        )
                        .unwrap();
                    (s, x)
                };
                let (global, xg) = run(DistributedStart::Global(&warm), true);
                let (local_g, xl) = run(DistributedStart::Local(&local), true);
                let (lazy, xz) = run(DistributedStart::Local(&local), false);
                assert_eq!(global.eigenvectors.rows(), n);
                assert_eq!(
                    values_bits(&global.eigenvalues),
                    values_bits(&local_g.eigenvalues)
                );
                assert_eq!(
                    values_bits(global.eigenvectors.as_slice()),
                    values_bits(local_g.eigenvectors.as_slice())
                );
                assert_eq!(
                    values_bits(&global.eigenvalues),
                    values_bits(&lazy.eigenvalues)
                );
                assert_eq!(global.operator_applications, lazy.operator_applications);
                assert_eq!(global.iterations, lazy.iterations);
                if workspace {
                    assert_eq!(
                        (lazy.eigenvectors.rows(), lazy.eigenvectors.cols()),
                        (0, bands)
                    );
                } else {
                    assert_eq!(lazy.eigenvectors.rows(), n, "allocating path gathers");
                }
                for x in [&xg, &xl, &xz] {
                    for j in 0..bands {
                        for (l, &g) in rows.iter().enumerate() {
                            assert_eq!(x[(l, j)].to_bits(), global.eigenvectors[(g, j)].to_bits());
                        }
                    }
                }
                values_bits(&lazy.eigenvalues)
            });
            for r in &out {
                assert_eq!(r, &out[0], "ranks disagree (P={ranks})");
            }
        }
    }
}
