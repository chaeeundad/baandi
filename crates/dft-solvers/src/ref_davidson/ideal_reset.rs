//! Guarded contraction experiment, derived from the own-words Davidson note.
//! No guard is an energy certificate; independent SCF verification stays required.
use dft_numerics::switches::{Id, RunConfig};
use std::cell::Cell;

// General safeguards, not material-tuned parameters. Match the existing
// restart orthogonality contract; limit unmeasured off-diagonal drift to three
// consecutive ideal contractions, even when all column norms remain one.
pub(super) const ORTHOGONALITY_LIMIT: f64 = 1e-12;
pub(super) const REFRESH_INTERVAL: usize = 4;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ResetStats {
    pub attempts: usize,
    pub ideal: usize,
    pub periodic_refreshes: usize,
    /// Guard trips, including off-diagonal defects found at periodic refresh.
    pub guard_fallbacks: usize,
    /// Subsequent measured restarts after a guard trip in the same solve.
    pub disabled_refreshes: usize,
    /// Maximum entrywise |S-I| at periodic full measurements, dimensionless.
    pub max_measured_defect: f64,
}
impl ResetStats {
    pub(super) fn add(&mut self, other: Self) {
        self.attempts += other.attempts;
        self.ideal += other.ideal;
        self.periodic_refreshes += other.periodic_refreshes;
        self.guard_fallbacks += other.guard_fallbacks;
        self.disabled_refreshes += other.disabled_refreshes;
        self.max_measured_defect = self.max_measured_defect.max(other.max_measured_defect);
    }
}
thread_local! {
    static TOTAL: Cell<(u64, ResetStats)> = const { Cell::new((0, ResetStats {
        attempts: 0, ideal: 0, periodic_refreshes: 0, guard_fallbacks: 0,
        disabled_refreshes: 0, max_measured_defect: 0.0,
    })) };
}
/// Consecutive successful solves on this thread for one run ID; changing IDs
/// resets this bounded summary. Per-solve Record counters remain complete.
/// Publish outside the numerical loop through the existing run manifest, so
/// CLI consumers can read fallback counts without a new serialization contract.
pub(super) fn publish(stats: ResetStats) {
    let config = RunConfig::current();
    if !config.on(Id::REF_DAV_IDEAL_RESET) {
        return;
    }
    let total = TOTAL.with(|slot| {
        let (id, mut total) = slot.get();
        if id != config.run_id() {
            total = ResetStats::default();
        }
        total.add(stats);
        slot.set((config.run_id(), total));
        total
    });
    config.note_path(
        Id::REF_DAV_IDEAL_RESET,
        &format!("successful solves, rank-local: attempts={} ideal={} periodic={} guard={} disabled={} max_measured_defect={:.17e}",
            total.attempts, total.ideal, total.periodic_refreshes,
            total.guard_fallbacks, total.disabled_refreshes, total.max_measured_defect),
        (total.guard_fallbacks > 0).then_some("orthogonality guard: measured refresh until solve end"),
    );
}

#[cfg(test)]
thread_local! {
    static AUDIT: std::cell::RefCell<Option<Vec<f64>>> = const { std::cell::RefCell::new(None) };
}
#[cfg(test)]
pub(super) fn start_audit() {
    AUDIT.with(|s| *s.borrow_mut() = Some(Vec::new()));
}
#[cfg(test)]
pub(super) fn take_audit() -> Vec<f64> {
    AUDIT.with(|s| s.borrow_mut().take().unwrap())
}
#[cfg(test)]
pub(super) fn audit_restart<C: dft_mpi::Collective>(c: &C, v: &dft_numerics::RealMatrix) {
    if !AUDIT.with(|s| s.borrow().is_some()) {
        return;
    }
    let n = v.cols();
    let mut local = vec![0.0; n * n];
    let mut global = vec![0.0; n * n];
    for j in 0..n {
        for i in 0..n {
            local[i + j * n] = (0..v.rows()).map(|k| v[(k, i)] * v[(k, j)]).sum();
        }
    }
    c.tree_sum_f64(&local, &mut global).unwrap();
    let defect = global
        .iter()
        .enumerate()
        .map(|(k, x)| {
            assert!(x.is_finite());
            (x - f64::from(k % n == k / n)).abs()
        })
        .fold(0.0, f64::max);
    AUDIT.with(|s| s.borrow_mut().as_mut().unwrap().push(defect));
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ref_davidson::{Options, ReducedPlacement, solve};
    use crate::{SolverError, distributed_davidson::DistributedRealOperator};
    use dft_mpi::{Collective, block_counts, simulate};
    use dft_numerics::{RealMatrix, symmetric_eigh};

    struct Dense {
        h: RealMatrix,
        rows: Vec<usize>,
        counts: Vec<usize>,
        full: Vec<f64>,
    }
    impl DistributedRealOperator for Dense {
        fn global_len(&self) -> usize {
            self.h.rows()
        }
        fn local_rows(&self) -> &[usize] {
            &self.rows
        }
        fn apply_block<C: Collective>(
            &mut self,
            c: &C,
            x: &RealMatrix,
            y: &mut RealMatrix,
        ) -> Result<(), SolverError> {
            let dim = self.h.rows();
            for j in 0..x.cols() {
                c.allgatherv_f64(
                    &x.as_slice()[j * x.rows()..(j + 1) * x.rows()],
                    &self.counts,
                    &mut self.full,
                )
                .unwrap();
                for (i, &row) in self.rows.iter().enumerate() {
                    y[(i, j)] = (0..dim).map(|k| self.h[(row, k)] * self.full[k]).sum();
                }
            }
            Ok(())
        }
    }
    // Same DCT spectrum and start as tests/ref_davidson.rs, including the
    // historically rejected degenerate occupied-edge fixture (32 rows, 8 bands).
    fn fixture(dim: usize, cluster: bool) -> RealMatrix {
        let q = |i: usize, k: usize| {
            if k == 0 {
                1.0 / (dim as f64).sqrt()
            } else {
                (2.0 / dim as f64).sqrt()
                    * (std::f64::consts::PI * (i as f64 + 0.5) * k as f64 / dim as f64).cos()
            }
        };
        let values: Vec<_> = (0..dim)
            .map(|i| {
                if cluster {
                    match i {
                        2 | 3 => 0.3,
                        4 | 5 => 0.5,
                        _ => 0.1 * i as f64,
                    }
                } else {
                    0.2 * i as f64
                }
            })
            .collect();
        let mut h = RealMatrix::zeros(dim, dim).unwrap();
        for j in 0..dim {
            for i in 0..dim {
                h[(i, j)] = (0..dim).map(|k| q(i, k) * values[k] * q(j, k)).sum();
            }
        }
        h
    }

    #[test]
    fn ideal_reset_accuracy_decisions_and_restart_orthogonality() {
        assert!(
            RunConfig::from_pairs(&[])
                .unwrap()
                .on(Id::REF_DAV_IDEAL_RESET) // default on since 2026-11-25
        );
        for (dim, bands, occupied, capacity, cluster) in [
            (32, 8, 6, 16, false),
            (32, 8, 6, 16, true),
            (32, 8, 6, 24, true),
            (7, 4, 3, 7, true),
        ] {
            let h = fixture(dim, cluster);
            let (exact, _) = symmetric_eigh(&h).unwrap();
            for ranks in [1, 3, 8] {
                for placement in [ReducedPlacement::Root, ReducedPlacement::Replicated] {
                    simulate(ranks, |c| {
                        let counts = block_counts(dim, ranks);
                        let offset: usize = counts[..c.rank()].iter().sum();
                        let rows: Vec<_> = (offset..offset + counts[c.rank()]).collect();
                        let mut x = RealMatrix::zeros(rows.len(), bands).unwrap();
                        for j in 0..bands {
                            for (i, &row) in rows.iter().enumerate() {
                                x[(i, j)] = if row == j {
                                    1.0
                                } else {
                                    0.03 * ((row * 31 + j * 17 + 3) as f64).sin()
                                };
                            }
                        }
                        let diagonal: Vec<_> = rows.iter().map(|&i| h[(i, i)]).collect();
                        let mut op = Dense {
                            h: h.clone(),
                            rows,
                            counts,
                            full: vec![0.0; dim],
                        };
                        let mut options = Options::new(bands, 1e-13, occupied);
                        options.capacity = capacity;
                        options.first_step = true;
                        options.max_iterations = 400;
                        options.placement = placement;
                        options.return_hx = true;
                        let run = |op: &mut Dense, on| {
                            let config =
                                RunConfig::from_pairs(&[("ref_dav_ideal_reset", on)]).unwrap();
                            let _scope = config.enter();
                            start_audit();
                            let s = solve(c, op, &x, &diagonal, options).unwrap();
                            let defects = take_audit();
                            assert_eq!(defects.len(), s.record.restarts);
                            assert!(s.converged.iter().all(|&v| v));
                            // Fresh H independently checks the stored rotated products.
                            let mut fresh = RealMatrix::zeros(x.rows(), bands).unwrap();
                            op.apply_block(c, &s.eigenvectors, &mut fresh).unwrap();
                            let hx = s.hx.as_ref().unwrap();
                            assert!(
                                fresh
                                    .as_slice()
                                    .iter()
                                    .zip(hx.as_slice())
                                    .all(|(a, b)| (a - b).abs() < 2e-11)
                            );
                            if on == "on" {
                                assert!(
                                    config
                                        .path_note(Id::REF_DAV_IDEAL_RESET)
                                        .unwrap()
                                        .path
                                        .contains("attempts=")
                                );
                                assert_eq!(s.record.resets.attempts, s.record.restarts);
                            } else {
                                assert_eq!(s.record.resets, ResetStats::default());
                            }
                            (s, defects)
                        };
                        let (a, ad) = run(&mut op, "off");
                        let (b, bd) = run(&mut op, "on");
                        let delta = a
                            .eigenvalues
                            .iter()
                            .zip(&b.eigenvalues)
                            .map(|(a, b)| (a - b).abs())
                            .fold(0.0, f64::max);
                        let dense = exact
                            .iter()
                            .zip(&b.eigenvalues)
                            .map(|(a, b)| (a - b).abs())
                            .fold(0.0, f64::max);
                        let energy = |s: &crate::ref_davidson::Solution| {
                            2.0 * s.eigenvalues[..occupied].iter().sum::<f64>()
                        };
                        let de = (energy(&a) - energy(&b)).abs();
                        assert!(
                            delta < 1e-10 && dense < 1e-10 && de < 1e-10,
                            "delta={delta:e} dense={dense:e} band_energy={de:e} Ha"
                        );
                        let defect = bd.iter().copied().fold(0.0, f64::max);
                        assert!(defect < 1e-10, "restart defect={defect:e}");
                        let decisions = |s: &crate::ref_davidson::Solution| {
                            s.record
                                .iterations
                                .iter()
                                .map(|r| {
                                    (
                                        r.active_bands.clone(),
                                        r.reduced_order,
                                        r.restarts,
                                        r.repairs,
                                        r.dropped_corrections,
                                        r.unconverged,
                                    )
                                })
                                .collect::<Vec<_>>()
                        };
                        let same = decisions(&a) == decisions(&b);
                        assert!(same, "fixture decisions changed: dim={dim} ranks={ranks}");
                        assert_eq!(a.record.h_columns, b.record.h_columns);
                        assert_eq!(a.record.h_calls, b.record.h_calls);
                        if c.rank() == 0 {
                            eprintln!(
                                "RESET dim={dim} bands={bands} cap={capacity} cluster={cluster} ranks={ranks} placement={placement:?} H={}/{} iterations={}/{} restarts={}/{} decisions_equal={same} root_delta_Ha={delta:.17e} dense_delta_Ha={dense:.17e} band_energy_delta_Ha={de:.17e} resets={:?} orth_off={ad:?} orth_on={bd:?}",
                                a.record.h_columns,
                                b.record.h_columns,
                                a.record.iterations.len(),
                                b.record.iterations.len(),
                                a.record.restarts,
                                b.record.restarts,
                                b.record.resets
                            );
                        }
                    });
                }
            }
        }
    }
}
