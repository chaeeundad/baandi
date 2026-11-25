use super::*;
use crate::ref_davidson::{Record, Solution, solve, solve_with_initial_hx};
use dft_numerics::switches::RunConfig;

struct Operator {
    dim: usize,
    rows: Vec<usize>,
    clustered: bool,
    threads: usize,
}
impl Operator {
    fn diagonal(&self, row: usize) -> f64 {
        if self.clustered && row < 2 {
            1.
        } else {
            1. + 0.4 * row as f64
        }
    }
}
impl DistributedRealOperator for Operator {
    fn global_len(&self) -> usize {
        self.dim
    }
    fn local_rows(&self) -> &[usize] {
        &self.rows
    }
    fn inner_threads(&self) -> usize {
        self.threads
    }
    fn apply_block<C: Collective>(
        &mut self,
        c: &C,
        x: &RealMatrix,
        y: &mut RealMatrix,
    ) -> Result<(), SolverError> {
        for j in 0..x.cols() {
            let dot: f64 = self
                .rows
                .iter()
                .enumerate()
                .map(|(i, &row)| ((row + 1) as f64).sin() * x[(i, j)])
                .sum();
            let mut sum = [0.];
            c.tree_sum_f64(&[dot], &mut sum).map_err(error)?;
            for (i, &row) in self.rows.iter().enumerate() {
                y[(i, j)] =
                    self.diagonal(row) * x[(i, j)] + 0.03 * ((row + 1) as f64).sin() * sum[0];
            }
        }
        Ok(())
    }
}
fn config(mask: usize) -> RunConfig {
    RunConfig::from_pairs(&[
        ("ref_dav_triangle", if mask & 1 != 0 { "on" } else { "off" }),
        (
            "ref_dav_stack_rotate",
            if mask & 2 != 0 { "on" } else { "off" },
        ),
        (
            "ref_dav_setup_fused",
            if mask & 4 != 0 { "on" } else { "off" },
        ),
        // Count cold allocation traffic consistently in each arm.
        ("ref_dav_reuse_rows", "off"),
    ])
    .unwrap()
}
fn decisions(s: &Solution) -> String {
    let r = &s.record;
    let trace: Vec<_> = std::iter::once(&r.setup)
        .chain(&r.iterations)
        .chain(std::iter::once(&r.finish))
        .map(|r| {
            (
                &r.active_bands,
                r.basis_before,
                r.h_calls,
                r.h_columns,
                r.reduced_order,
                r.restarts,
                r.repairs,
                r.dropped_corrections,
                r.unconverged,
                (
                    r.resets.ideal,
                    r.resets.guard_fallbacks,
                    r.resets.periodic_refreshes,
                    r.resets.disabled_refreshes,
                ),
            )
        })
        .collect();
    format!(
        "{:?} {:?} {trace:?}",
        s.converged,
        (
            r.h_calls,
            r.h_columns,
            r.reduced_solves,
            r.backend_solves,
            r.communication
        )
    )
}
fn assert_repeat(a: &Solution, b: &Solution) {
    assert_eq!(decisions(a), decisions(b));
    for (a, b) in [
        (&a.eigenvalues[..], &b.eigenvalues[..]),
        (a.eigenvectors.as_slice(), b.eigenvectors.as_slice()),
        (
            a.hx.as_ref().unwrap().as_slice(),
            b.hx.as_ref().unwrap().as_slice(),
        ),
    ] {
        assert_eq!(
            a.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            b.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
        );
    }
    assert_eq!(a.ritz_values_ha, b.ritz_values_ha);
}

#[test]
fn p4_arms_preserve_decisions_h_columns_repeats_and_threads() {
    let mut max_delta = 0.0_f64;
    for (dim, n, cap, clustered) in [
        (12, 3, 6, false),
        (17, 5, 15, true),
        (7, 4, 7, true),
        (513, 3, 6, false),
    ] {
        for ranks in [1, 3, 8] {
            let deltas = dft_mpi::simulate(ranks, |c| {
                let mut op = Operator {
                    dim,
                    rows: (c.rank()..dim).step_by(ranks).collect(),
                    clustered,
                    threads: 1,
                };
                let rows = op.rows.len();
                let mut x = RealMatrix::zeros(rows, n).unwrap();
                for j in 0..n {
                    for (i, &row) in op.rows.iter().enumerate() {
                        x[(i, j)] = if row == j {
                            1.
                        } else {
                            0.03 * ((row * 7 + j * 13 + 1) as f64).sin()
                        };
                    }
                }
                let diag: Vec<_> = op
                    .rows
                    .iter()
                    .map(|&row| op.diagonal(row) + 0.03 * ((row + 1) as f64).sin().powi(2))
                    .collect();
                let mut local_delta = 0.0_f64;
                for placement in [ReducedPlacement::Root, ReducedPlacement::Replicated] {
                    let mut o = Options::new(n, 1e-10, n - 1);
                    o.capacity = cap;
                    o.profile = true;
                    o.return_hx = true;
                    o.placement = placement;
                    o.max_iterations = if dim > 500 { 2 } else { 40 };
                    work::take_solves();
                    let a = {
                        let _g = config(0).enter();
                        solve(c, &mut op, &x, &diag, o).unwrap()
                    };
                    let aw = work::take_solves().pop().unwrap();
                    for mask in [1, 2, 4, 7] {
                        let _g = config(mask).enter();
                        let b = solve(c, &mut op, &x, &diag, o).unwrap();
                        let bw = work::take_solves().pop().unwrap();
                        assert_eq!(
                            decisions(&a),
                            decisions(&b),
                            "dim={dim} ranks={ranks} mask={mask}; Ritz operands: {:?} vs {:?}",
                            a.ritz_values_ha,
                            b.ritz_values_ha
                        );
                        for (a, b) in a.eigenvalues.iter().zip(&b.eigenvalues) {
                            local_delta = local_delta.max((a - b).abs());
                            assert!((a - b).abs() < 1e-10, "eigenvalue delta in Ha");
                        }
                        if mask == 4 {
                            assert_repeat(&a, &b);
                        }
                        if mask & 2 != 0 {
                            for stage in [Stage::Restart, Stage::Return] {
                                let a = aw.stages[stage as usize];
                                let b = bw.stages[stage as usize];
                                assert_eq!(a.gemm_flops, b.gemm_flops);
                                assert_eq!(a.gemm_calls, 2 * b.gemm_calls);
                            }
                        }
                        let triangle_saved: u64 = b
                            .record
                            .setup
                            .gemms
                            .iter()
                            .chain(b.record.iterations.iter().flat_map(|r| &r.gemms))
                            .filter(|g| g.stage.ends_with("_gemmt"))
                            .map(|g| g.k as u64 * g.m as u64 * g.m.saturating_sub(1) as u64)
                            .sum();
                        assert_eq!(a.record.gemm_flops - b.record.gemm_flops, triangle_saved);
                        assert_eq!(
                            bw.stages.iter().map(|s| s.gemm_flops).sum::<u64>(),
                            b.record.gemm_flops
                        );
                        if mask == 4 {
                            assert_eq!(
                                aw.stages[Stage::Setup as usize].streamed_bytes
                                    - bw.stages[Stage::Setup as usize].streamed_bytes,
                                8 * (rows * n + rows) as u64
                            );
                        }
                        let again = solve(c, &mut op, &x, &diag, o).unwrap();
                        assert_repeat(&b, &again);
                        work::take_solves();
                        op.threads = 4;
                        assert_repeat(&b, &solve(c, &mut op, &x, &diag, o).unwrap());
                        op.threads = 1;
                        work::take_solves();
                    }
                }
                local_delta
            });
            for delta in deltas {
                max_delta = max_delta.max(delta);
            }
        }
    }
    eprintln!(
        "P4-DAVREM 12 shape/rank fixtures, 2 placements, 4 arms: max eigenvalue delta={max_delta:e} Ha"
    );
}

#[test]
fn p4_setup_kept_hx_and_collective_invalid_inputs() {
    for ranks in [1, 8] {
        dft_mpi::simulate(ranks, |c| {
            let mut op = Operator {
                dim: 7,
                rows: (c.rank()..7).step_by(ranks).collect(),
                clustered: false,
                threads: 1,
            };
            let rows = op.rows.len();
            let mut x = RealMatrix::zeros(rows, 3).unwrap();
            for (i, &row) in op.rows.iter().enumerate() {
                if row < 3 {
                    x[(i, row)] = 1.;
                }
            }
            let mut hx = RealMatrix::zeros(rows, 3).unwrap();
            op.apply_block(c, &x, &mut hx).unwrap();
            let diag: Vec<_> = op.rows.iter().map(|&row| op.diagonal(row)).collect();
            let mut o = Options::new(3, 1e-10, 3);
            o.profile = true;
            o.return_hx = true;
            work::take_solves();
            let a = {
                let _g = config(0).enter();
                solve_with_initial_hx(c, &mut op, &x, &diag, o, Some(&hx)).unwrap()
            };
            let aw = work::take_solves().pop().unwrap();
            let _g = config(4).enter();
            let b = solve_with_initial_hx(c, &mut op, &x, &diag, o, Some(&hx)).unwrap();
            let bw = work::take_solves().pop().unwrap();
            assert_repeat(&a, &b);
            assert_eq!(
                aw.stages[Stage::Setup as usize].streamed_bytes
                    - bw.stages[Stage::Setup as usize].streamed_bytes,
                8 * (2 * rows * 3 + rows) as u64
            );
            for field in 0..4 {
                let mut badx = x.clone();
                let mut badhx = hx.clone();
                let mut badd = diag.clone();
                if c.rank() == 0 {
                    match field {
                        0 => badx[(0, 0)] = f64::NAN,
                        1 => badhx[(0, 0)] = f64::INFINITY,
                        2 => badd[0] = f64::NAN,
                        _ => {
                            badhx.reshape_within_storage(rows, 2).unwrap();
                        }
                    }
                }
                assert!(solve_with_initial_hx(c, &mut op, &badx, &badd, o, Some(&badhx)).is_err());
                assert!(work::take_solves().is_empty());
            }
        });
    }
}

#[test]
fn p4_upper_border_matches_independent_dot_and_ignores_stale_lower() {
    for rows in [0, 1, 257] {
        let op = Operator {
            dim: rows.max(1),
            rows: (0..rows).collect(),
            clustered: false,
            threads: 1,
        };
        dft_mpi::simulate(1, |c| {
            let mut w = Workspace::new(&op, 3, 7, true).unwrap();
            w.triangle = true;
            for j in 0..7 {
                for i in 0..rows {
                    w.v[(i, j)] = ((i + 3 * j) % 7) as f64 - 3.;
                    w.hv[(i, j)] = ((2 * i + j) % 5) as f64 - 2.;
                }
            }
            for from in [0, 3, 6] {
                for is_h in [false, true] {
                    w.border.as_mut_slice().fill(f64::NAN);
                    let mut r = Record::new(1, 7, true, ReducedPlacement::Root);
                    w.border(c, from, 7, is_h, &mut r.setup, true).unwrap();
                    let target = if is_h { &w.hr } else { &w.lean.s };
                    let b = if is_h { &w.hv } else { &w.v };
                    for j in from..7 {
                        for i in 0..=j {
                            let expected: f64 = (0..rows).map(|q| w.v[(q, i)] * b[(q, j)]).sum();
                            assert_eq!(target[(i, j)], expected);
                            assert_eq!(target[(j, i)], expected);
                        }
                    }
                    let flops: u64 = r.setup.gemms.iter().map(|g| g.flops).sum();
                    let k = 7 - from;
                    let expected = if dft_numerics::MKL_BACKEND {
                        2 * rows * (7 * k - k * (k - 1) / 2)
                    } else {
                        2 * rows * 7 * k
                    };
                    assert_eq!(flops, expected as u64);
                }
            }
        });
    }
}

#[test]
fn p4_stacked_rotation_scalar_oracle_traffic_and_x_only_return() {
    let _g = config(2).enter();
    for rows in [0, 1, 257] {
        let op = Operator {
            dim: rows.max(1),
            rows: (0..rows).collect(),
            clustered: false,
            threads: 1,
        };
        let (m, n) = (5, 3);
        let mut w = Workspace::new(&op, n, 7, true).unwrap();
        w.v.reshape_within_storage(rows, m).unwrap();
        w.hv.reshape_within_storage(rows, m).unwrap();
        w.lean.c.reshape_within_storage(m, n).unwrap();
        for j in 0..m {
            for i in 0..rows {
                w.v[(i, j)] = ((i + 3 * j) % 7) as f64 - 3.;
                w.hv[(i, j)] = ((2 * i + j) % 5) as f64 - 2.;
            }
        }
        for j in 0..n {
            for i in 0..m {
                w.lean.c[(i, j)] = ((i + j) % 3) as f64 - 1.;
            }
        }
        work::take_solves();
        let scope = work::Scope::new(true, 1, 7);
        let mut r = Record::new(1, 7, true, ReducedPlacement::Root);
        w.finish(m, true, &mut r.finish, true).unwrap();
        scope.finish();
        let counted = work::take_solves().pop().unwrap().stages[Stage::Return as usize];
        assert_eq!(counted.gemm_calls, 1);
        assert_eq!(counted.gemm_flops, 4 * (rows * m * n) as u64);
        assert_eq!(counted.streamed_bytes, 32 * (rows * (m + n)) as u64);
        for j in 0..n {
            for i in 0..rows {
                let x: f64 = (0..m).map(|q| w.v[(i, q)] * w.lean.c[(q, j)]).sum();
                let hx: f64 = (0..m).map(|q| w.hv[(i, q)] * w.lean.c[(q, j)]).sum();
                assert_eq!(w.x[(i, j)], x);
                assert_eq!(w.hx[(i, j)], hx);
            }
        }
        w.hx.as_mut_slice().fill(123.);
        let mut r = Record::new(1, 7, true, ReducedPlacement::Root);
        w.finish(m, false, &mut r.finish, true).unwrap();
        assert_eq!(r.finish.gemms.len(), 1);
        assert_eq!(r.finish.gemms[0].stage, "return_x");
        assert!(w.hx.as_slice().iter().all(|v| *v == 123.));
    }
}
