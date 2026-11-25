//! Agreement storage lifetime, growth and failure recovery (eigenvalues in Ha).
use dft_mpi::{Collective, simulate};
use dft_numerics::real::RealMatrix;
use dft_solvers::distributed_davidson::{
    DistributedRealOperator, DistributedTuning, release_workspace, scratch_allocation_stats,
    set_scratch_allocation_profile,
};
use dft_solvers::{BandTolerances, RealDavidson, SolverError};
use std::time::Duration;

struct Diagonal {
    n: usize,
    rows: Vec<usize>,
    fault: u8,
}
impl DistributedRealOperator for Diagonal {
    fn global_len(&self) -> usize {
        self.n
    }
    fn local_rows(&self) -> &[usize] {
        &self.rows
    }
    fn apply_block<C: Collective>(
        &mut self,
        comm: &C,
        x: &RealMatrix,
        y: &mut RealMatrix,
    ) -> Result<(), SolverError> {
        let mut sum = [0.0];
        comm.tree_sum_f64(&[1.0], &mut sum).unwrap();
        for j in 0..x.cols() {
            for (i, &r) in self.rows.iter().enumerate() {
                y[(i, j)] = (r + 1) as f64 * x[(i, j)];
            }
        }
        if comm.rank() == 1 {
            if self.fault == 1 {
                return Err(SolverError::Application(
                    "rank-local injected failure".into(),
                ));
            }
            if self.fault == 2 && !y.as_slice().is_empty() {
                y[(0, 0)] = f64::NAN;
            }
        }
        Ok(())
    }
}

#[test]
fn warm_growth_shrink_and_failed_then_healthy_solves() {
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        // Rank 3 is empty for the small shape. Exercise the frozen workspace.
        {
            simulate(4, |comm| {
                release_workspace();
                set_scratch_allocation_profile(true);
                let mut healthy = std::collections::BTreeMap::new();
                for (n, bands, fault, growth) in [
                    (3, 2, 0, true),
                    (3, 2, 0, false),
                    (9, 3, 1, true),
                    (9, 3, 0, false),
                    (9, 3, 2, false),
                    (9, 3, 0, false),
                    (3, 2, 1, false),
                    (3, 2, 0, false),
                ] {
                    let mut op = Diagonal {
                        n,
                        rows: (0..n).filter(|r| r % 4 == comm.rank()).collect(),
                        fault,
                    };
                    let diag: Vec<_> = op.rows.iter().map(|r| (r + 1) as f64).collect();
                    let mut start = RealMatrix::zeros(n, bands).unwrap();
                    for j in 0..bands {
                        start[(j, j)] = 1.0;
                    }
                    let result = RealDavidson {
                        bands,
                        tolerance: 1e-12,
                        max_iterations: 20,
                        max_subspace: n,
                        seed: 42,
                    }
                    .solve_distributed_rows_tuned(
                        comm,
                        &mut op,
                        Some(&start),
                        &diag,
                        BandTolerances {
                            strict_bands: bands,
                            buffer_tolerance: 1e-12,
                        },
                        true,
                        false,
                        6,
                        DistributedTuning {
                            ..DistributedTuning::from_env()
                        },
                    );
                    let stats = scratch_allocation_stats().unwrap();
                    assert_eq!(stats.allocations, usize::from(growth));
                    if !growth {
                        assert_eq!(stats.previous_allocations, 0);
                        assert_eq!(stats.previous_bytes, 0);
                    }
                    assert_eq!(
                        stats.bytes,
                        if growth {
                            (2 * n * n + 2 * n + 1) * 8
                        } else {
                            0
                        }
                    );
                    if fault != 0 {
                        assert!(result.is_err(), "every rank must reject the failed apply");
                    } else {
                        let s = result.unwrap();
                        for (j, &value) in s.eigenvalues.iter().enumerate() {
                            assert!((value - (j + 1) as f64).abs() < 1e-12);
                        }
                        let bits: Vec<_> = s
                            .eigenvalues
                            .iter()
                            .chain(s.eigenvectors.as_slice())
                            .chain(&s.residuals)
                            .map(|x| x.to_bits())
                            .collect();
                        if let Some(before) = healthy.insert(n, bits.clone()) {
                            assert_eq!(before, bits, "reuse must preserve every output bit");
                        }
                    }
                }
                set_scratch_allocation_profile(false);
                assert_eq!(scratch_allocation_stats(), None);
                release_workspace();
            });
        }
        tx.send(()).unwrap();
    });
    rx.recv_timeout(Duration::from_secs(20))
        .expect("agreement must not hang");
}
