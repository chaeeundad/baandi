//! Fixed-H oracles: sorted roots AND projectors, never residuals alone.
use dft_mpi::{Collective, simulate};
use dft_numerics::{RealMatrix, symmetric_eigh};
use dft_solvers::distributed_davidson::DistributedRealOperator;
use dft_solvers::ref_davidson::{self, Options, ReducedPlacement};
use dft_solvers::{RealDavidson, RealLinearOperator, SolverError};

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
        let n = self.h.rows();
        for j in 0..x.cols() {
            c.allgatherv_f64(
                &x.as_slice()[j * x.rows()..(j + 1) * x.rows()],
                &self.counts,
                &mut self.full,
            )
            .unwrap();
            for (i, &row) in self.rows.iter().enumerate() {
                y[(i, j)] = (0..n).map(|k| self.h[(row, k)] * self.full[k]).sum();
            }
        }
        Ok(())
    }
}
struct Serial(RealMatrix);
impl RealLinearOperator for Serial {
    fn dimension(&self) -> usize {
        self.0.rows()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), SolverError> {
        for (i, v) in y.iter_mut().enumerate() {
            *v = x.iter().enumerate().map(|(j, x)| self.0[(i, j)] * x).sum();
        }
        Ok(())
    }
}
fn fixture(n: usize, cluster: bool) -> RealMatrix {
    // Orthogonal discrete cosine basis; two degenerate pairs, including one
    // ending at the occupied edge (bands 0..6). An unrelated close empty pair.
    let vals: Vec<_> = (0..n)
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
    let q = |i: usize, k: usize| {
        if k == 0 {
            1.0 / (n as f64).sqrt()
        } else {
            (2.0 / n as f64).sqrt()
                * (std::f64::consts::PI * (i as f64 + 0.5) * k as f64 / n as f64).cos()
        }
    };
    let mut h = RealMatrix::zeros(n, n).unwrap();
    for j in 0..n {
        for i in 0..n {
            h[(i, j)] = (0..n).map(|k| q(i, k) * vals[k] * q(j, k)).sum();
        }
    }
    h
}
fn start(n: usize, b: usize) -> RealMatrix {
    let mut x = RealMatrix::zeros(n, b).unwrap();
    // Full generic overlap, deterministic and independent of the partition.
    for j in 0..b {
        for i in 0..n {
            x[(i, j)] = if i == j {
                1.0
            } else {
                0.03 * ((i * 31 + j * 17 + 3) as f64).sin()
            };
        }
    }
    x
}
fn projector_error(a: &RealMatrix, b: &RealMatrix, k: usize) -> f64 {
    let mut error: f64 = 0.0;
    for j in 0..a.rows() {
        for i in 0..a.rows() {
            let d: f64 = (0..k)
                .map(|l| a[(i, l)] * a[(j, l)] - b[(i, l)] * b[(j, l)])
                .sum();
            error = error.max(d.abs());
        }
    }
    error
}
#[test]
fn dense_legacy_rank_and_degenerate_edge_oracles() {
    let dim = 32;
    let bands = 8;
    let occupied = 6;
    for cluster in [false, true] {
        let h = fixture(dim, cluster);
        let initial = start(dim, bands);
        let (exact, q) = symmetric_eigh(&h).unwrap();
        let diagonal: Vec<_> = (0..dim).map(|i| h[(i, i)]).collect();
        let old = RealDavidson {
            bands,
            tolerance: 1e-11,
            max_iterations: 400,
            max_subspace: dim,
            seed: 42,
        }
        .solve(&Serial(h.clone()), Some(&initial), &diagonal)
        .unwrap();
        for ranks in [1, 2, 3, 8] {
            for placement in [ReducedPlacement::Replicated, ReducedPlacement::Root] {
                let results = simulate(ranks, |c| {
                    let counts = dft_mpi::block_counts(dim, ranks);
                    let offset = counts[..c.rank()].iter().sum::<usize>();
                    let rows: Vec<_> = (offset..offset + counts[c.rank()]).collect();
                    let mut x = RealMatrix::zeros(rows.len(), bands).unwrap();
                    for j in 0..bands {
                        for (i, &global) in rows.iter().enumerate() {
                            x[(i, j)] = initial[(global, j)];
                        }
                    }
                    let kinetic: Vec<_> = rows.iter().map(|&i| diagonal[i]).collect();
                    let mut op = Dense {
                        h: h.clone(),
                        rows,
                        counts: counts.clone(),
                        full: vec![0.0; dim],
                    };
                    let mut options = Options::new(bands, 1e-13, occupied);
                    options.first_step = true;
                    options.max_iterations = 400;
                    options.placement = placement;
                    options.return_hx = true;
                    let sol = ref_davidson::solve(c, &mut op, &x, &kinetic, options).unwrap();
                    assert!(sol.converged.iter().all(|&v| v), "{:?}", sol.eigenvalues);
                    let mut full = RealMatrix::zeros(dim, bands).unwrap();
                    for j in 0..bands {
                        c.allgatherv_f64(
                            &sol.eigenvectors.as_slice()[j * x.rows()..(j + 1) * x.rows()],
                            &counts,
                            &mut full.as_mut_slice()[j * dim..(j + 1) * dim],
                        )
                        .unwrap();
                    }
                    assert!(sol.record.restarts > 0);
                    assert_eq!(
                        sol.record.h_columns,
                        bands
                            + sol
                                .record
                                .iterations
                                .iter()
                                .map(|r| r.h_columns)
                                .sum::<usize>()
                    );
                    assert_eq!(sol.record.setup.h_columns, bands);
                    assert_eq!(
                        sol.record.iterations[0].active_bands,
                        (0..bands).collect::<Vec<_>>()
                    );
                    assert_eq!(
                        sol.record.reduced_solves,
                        1 + sol
                            .record
                            .iterations
                            .iter()
                            .filter(|r| r.h_columns > 0)
                            .count()
                    );
                    assert_eq!(sol.record.total_seconds, 0.0);
                    let mut fresh = RealMatrix::zeros(x.rows(), bands).unwrap();
                    op.apply_block(c, &sol.eigenvectors, &mut fresh).unwrap();
                    let hx = sol.hx.as_ref().unwrap();
                    assert!(
                        fresh
                            .as_slice()
                            .iter()
                            .zip(hx.as_slice())
                            .all(|(a, b)| (a - b).abs() < 2e-11)
                    );
                    (sol, full)
                });
                let (s, x) = &results[0];
                let de = s
                    .eigenvalues
                    .iter()
                    .zip(&exact)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0, f64::max);
                let dp = projector_error(x, &q, occupied);
                let dl = projector_error(x, &old.eigenvectors, occupied);
                eprintln!(
                    "cluster={cluster} ranks={ranks} placement={placement:?} iterations={} H={} dE={de:.3e} Ha projector={dp:.3e} legacy={dl:.3e} active={:?}",
                    s.record.iterations.len(),
                    s.record.h_columns,
                    s.record
                        .iterations
                        .iter()
                        .map(|r| r.active_bands.len())
                        .collect::<Vec<_>>()
                );
                assert!(de < 1e-11, "{de}");
                assert!(dp < 2e-6, "{dp}");
                assert!(dl < 2e-6, "{dl}");
                for (other, _) in &results {
                    assert_eq!(other.eigenvalues, s.eigenvalues);
                }
            }
        }
    }
}
#[test]
fn thresholds_preconditioner_and_invalid_input() {
    let mut o = Options::new(4, 2e-8, 2);
    assert_eq!(o.band_threshold(0), 2e-8);
    assert_eq!(o.band_threshold(3), 5e-6);
    o.first_step = true;
    assert_eq!(o.band_threshold(3), 2e-8);
    for d in [-1e8, -10.0, -1.0, 0.0, 0.5, 1.0, 100.0] {
        let p = ref_davidson::preconditioner_denominator(d, 0.0);
        assert!(p.is_finite() && p > 0.0);
        let direct = 0.25 * (1.0 + 2.0 * d + (1.0 + (2.0 * d - 1.0).powi(2)).sqrt());
        assert!((p - direct).abs() < 1e-8);
    }
    // Independent Ry expression, including shifts from local/nonlocal terms.
    for (h, e) in [(-0.7, -0.2), (0.35, -0.15), (2.1, 0.6)] {
        let x: f64 = 2.0 * (h - e);
        let d_ry = (1.0 + x + (1.0 + (x - 1.0).powi(2)).sqrt()) / 2.0;
        let actual = ref_davidson::preconditioner_denominator(h, e);
        assert!((2.0 * actual - d_ry).abs() < 2e-15);
        assert!(
            (actual - ref_davidson::preconditioner_denominator(h + 0.25, e + 0.25)).abs() < 2e-15
        );
    }
    let results = simulate(3, |c| {
        let mut op = Dense {
            h: fixture(9, false),
            rows: (3 * c.rank()..3 * c.rank() + 3).collect(),
            counts: vec![3; 3],
            full: vec![0.0; 9],
        };
        let mut x = start(3, 2);
        if c.rank() == 1 {
            x[(0, 0)] = f64::NAN;
        }
        ref_davidson::solve(c, &mut op, &x, &[1.0; 3], Options::new(2, 1e-8, 2)).is_err()
    });
    assert!(results.iter().all(|&v| v));
}
#[test]
fn exact_start_and_iteration_limit_are_explicit() {
    simulate(1, |c| {
        let dim = 12;
        let bands = 3;
        let mut h = RealMatrix::zeros(dim, dim).unwrap();
        for i in 0..dim {
            h[(i, i)] = i as f64;
        }
        let mut op = Dense {
            h,
            rows: (0..dim).collect(),
            counts: vec![dim],
            full: vec![0.0; dim],
        };
        let mut x = RealMatrix::zeros(dim, bands).unwrap();
        for j in 0..bands {
            x[(j, j)] = 1.0;
        }
        let mut o = Options::new(bands, 1e-12, bands);
        o.return_hx = true;
        let s = ref_davidson::solve(c, &mut op, &x, &vec![1.0; dim], o).unwrap();
        assert!(s.converged.iter().all(|&v| v));
        assert_eq!(s.record.h_columns, bands);
        assert_eq!(s.record.iterations[0].dropped_corrections, bands);
        o.max_iterations = 1;
        o.profile = true;
        let s = ref_davidson::solve(c, &mut op, &start(dim, bands), &vec![1.0; dim], o).unwrap();
        assert!(s.converged.iter().any(|&v| !v));
        assert_eq!(s.record.iterations.len(), 1);
        assert!(s.record.total_seconds > 0.0);
    });
}

#[test]
fn dependent_start_and_operator_failure_are_collective_errors() {
    struct Failing {
        rows: Vec<usize>,
    }
    impl DistributedRealOperator for Failing {
        fn global_len(&self) -> usize {
            9
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
            if c.rank() == 1 {
                return Err(SolverError::Application("injected rank failure".into()));
            }
            y.as_mut_slice().copy_from_slice(x.as_slice());
            Ok(())
        }
    }
    let failures = simulate(3, |c| {
        let mut op = Failing {
            rows: (3 * c.rank()..3 * c.rank() + 3).collect(),
        };
        ref_davidson::solve(
            c,
            &mut op,
            &start(3, 2),
            &[1.0; 3],
            Options::new(2, 1e-8, 2),
        )
        .is_err()
    });
    assert!(failures.into_iter().all(|x| x));
    let dependent = simulate(3, |c| {
        let mut op = Dense {
            h: fixture(9, false),
            rows: (3 * c.rank()..3 * c.rank() + 3).collect(),
            counts: vec![3; 3],
            full: vec![0.0; 9],
        };
        let mut x = RealMatrix::zeros(3, 2).unwrap();
        x.as_mut_slice().fill(1.0);
        ref_davidson::solve(c, &mut op, &x, &[1.0; 3], Options::new(2, 1e-8, 2)).is_err()
    });
    assert!(dependent.into_iter().all(|x| x));
}

#[test]
fn zero_row_ranks_and_optional_products() {
    let results = simulate(8, |c| {
        let counts = dft_mpi::block_counts(3, 8);
        let offset = counts[..c.rank()].iter().sum::<usize>();
        let rows: Vec<_> = (offset..offset + counts[c.rank()]).collect();
        let mut h = RealMatrix::zeros(3, 3).unwrap();
        h[(1, 1)] = 1.0;
        h[(2, 2)] = 2.0;
        let mut x = RealMatrix::zeros(rows.len(), 1).unwrap();
        for (i, &g) in rows.iter().enumerate() {
            x[(i, 0)] = if g == 0 { 1.0 } else { 0.0 };
        }
        let kinetic = vec![1.0; rows.len()];
        let mut op = Dense {
            h,
            rows,
            counts,
            full: vec![0.0; 3],
        };
        let s = ref_davidson::solve(c, &mut op, &x, &kinetic, Options::new(1, 1e-10, 1)).unwrap();
        assert!(s.hx.is_none());
        assert_eq!(s.eigenvalues, vec![0.0]);
        let reductions = s.record.setup.communication.reductions
            + s.record
                .iterations
                .iter()
                .map(|r| r.communication.reductions)
                .sum::<u64>();
        assert_eq!(reductions, s.record.communication.reductions);
    });
    assert_eq!(results.len(), 8);
}
