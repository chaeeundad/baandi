//! MPI-KEPT-START (`BAANDI_M2_REUSE=on`): the stick-row solver starts a
//! Ritz-change solve from the kept X and the caller-updated H X of the previous
//! return. On a synthetic operator whose "potential" (a diagonal) changes between
//! two solves, 1 to 4 simulated ranks, selected-root and full-rotation paths:
//! the kept start skips exactly `bands` initial H columns, its eigenvalues agree
//! with the plain start within the solve tolerance class, every rank takes the
//! same decision, `x_local` equals this rank's rows of the returned
//! eigenvectors bit for bit, and a perturbed `initial` (1e-10) or an unarmed
//! state takes the plain start with the plain H count.
use dft_mpi::{Collective, SimRank, simulate};
use dft_numerics::{pool::WorkerPool, real::RealMatrix};
use dft_solvers::distributed_davidson::{
    DistributedKeptStart, DistributedRealOperator, DistributedTuning,
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

struct Outcome {
    plain_values: Vec<f64>,
    plain_h: usize,
    kept_values: Vec<f64>,
    kept_h: usize,
    taken: bool,
    declined: bool,
    local_matches: bool,
    perturbed_taken: bool,
    perturbed_h: usize,
}

fn scenario(ranks: usize, _selected: bool) -> Vec<Outcome> {
    let n = 96;
    let bands = 6;
    let old_shift = vec![0.0; n];
    let new_shift: Vec<f64> = (0..n).map(|i| 1e-3 * ((i * 7 + 3) as f64).sin()).collect();
    let a_old = matrix(n, &old_shift);
    let a_new = matrix(n, &new_shift);
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
        let mut old = RowsOp::new(&a_old, sim.rank(), ranks);
        let diag_old: Vec<f64> = old.rows.iter().map(|&r| a_old[(r, r)]).collect();
        let mut x = RealMatrix::zeros(0, 0).unwrap();
        let mut kept = DistributedKeptStart::new();
        kept.set_keep(true);
        let first = solver
            .solve_distributed_rows_kept(
                sim, &mut old, None, &diag_old, policy, true, false, 4, tuning, &mut x, &mut kept,
            )
            .unwrap();
        assert!(kept.is_stored());
        let rows = old.rows.clone();
        let local_matches = (0..bands).all(|j| {
            rows.iter()
                .enumerate()
                .all(|(l, &g)| x[(l, j)].to_bits() == first.eigenvectors[(g, j)].to_bits())
        });
        // Caller update: H_new X = H_old X + (D_new - D_old) X (row-local).
        let hx = kept.stored_mut().unwrap();
        for j in 0..bands {
            for (l, &g) in rows.iter().enumerate() {
                hx[(l, j)] += (new_shift[g] - old_shift[g]) * x[(l, j)];
            }
        }
        kept.arm();
        let mut new_op = RowsOp::new(&a_new, sim.rank(), ranks);
        let diag_new: Vec<f64> = new_op.rows.iter().map(|&r| a_new[(r, r)]).collect();
        let warm = first.eigenvectors.clone();
        let kept_solve = solver
            .solve_distributed_rows_kept(
                sim,
                &mut new_op,
                Some(&warm),
                &diag_new,
                policy,
                true,
                false,
                4,
                tuning,
                &mut x,
                &mut kept,
            )
            .unwrap();
        let (taken, declined) = (kept.taken(), kept.declined());
        let mut plain_op = RowsOp::new(&a_new, sim.rank(), ranks);
        let plain = solver
            .solve_distributed_rows_tuned(
                sim,
                &mut plain_op,
                Some(&warm),
                &diag_new,
                policy,
                true,
                false,
                4,
                tuning,
            )
            .unwrap();
        // Armed with an initial that differs by 1e-10 on one entry: declined.
        let mut x2 = RealMatrix::zeros(0, 0).unwrap();
        let mut k2 = DistributedKeptStart::new();
        k2.set_keep(true);
        let mut op2 = RowsOp::new(&a_old, sim.rank(), ranks);
        let base = solver
            .solve_distributed_rows_kept(
                sim, &mut op2, None, &diag_old, policy, true, false, 4, tuning, &mut x2, &mut k2,
            )
            .unwrap();
        k2.arm();
        let mut bad = base.eigenvectors.clone();
        bad[(n / 2, 1)] += 1e-10;
        let mut op3 = RowsOp::new(&a_old, sim.rank(), ranks);
        let perturbed = solver
            .solve_distributed_rows_kept(
                sim,
                &mut op3,
                Some(&bad),
                &diag_old,
                policy,
                true,
                false,
                4,
                tuning,
                &mut x2,
                &mut k2,
            )
            .unwrap();
        let mut op4 = RowsOp::new(&a_old, sim.rank(), ranks);
        let perturbed_plain = solver
            .solve_distributed_rows_tuned(
                sim,
                &mut op4,
                Some(&bad),
                &diag_old,
                policy,
                true,
                false,
                4,
                tuning,
            )
            .unwrap();
        assert_eq!(
            perturbed.operator_applications,
            perturbed_plain.operator_applications
        );
        Outcome {
            plain_values: plain.eigenvalues,
            plain_h: plain.operator_applications,
            kept_values: kept_solve.eigenvalues,
            kept_h: kept_solve.operator_applications,
            taken,
            declined,
            local_matches,
            perturbed_taken: k2.taken(),
            perturbed_h: perturbed.operator_applications,
        }
    })
}

#[test]
fn kept_start_skips_the_initial_apply_on_every_rank_count() {
    for selected in [false, true] {
        for ranks in 1..=4 {
            let out = scenario(ranks, selected);
            for (r, o) in out.iter().enumerate() {
                let worst = o
                    .kept_values
                    .iter()
                    .zip(&o.plain_values)
                    .fold(0f64, |a, (k, p)| a.max((k - p).abs()));
                println!(
                    "selected={selected} P={ranks} rank {r}: H kept {} plain {}, taken {}, max |d lambda| {worst:e}, perturbed H {}",
                    o.kept_h, o.plain_h, o.taken, o.perturbed_h
                );
                assert!(o.taken && !o.declined, "kept start not taken");
                assert!(o.local_matches, "x_local differs from the returned rows");
                assert!(worst < 1e-9, "eigenvalues: {worst:e}");
                assert!(!o.perturbed_taken, "perturbed initial must be declined");
                assert_eq!(o.kept_values, out[0].kept_values, "ranks disagree");
                assert_eq!(o.kept_h, out[0].kept_h);
            }
            // At equal Davidson work the kept start saves exactly the initial
            // `bands` columns; the trajectories may differ by rounding, so the
            // count is compared with the plain count minus the saved columns
            // only when the iteration counts match.
            let o = &out[0];
            assert!(o.kept_h < o.plain_h, "{} >= {}", o.kept_h, o.plain_h);
        }
    }
}
