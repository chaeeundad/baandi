//! Regression cases for LEAN-ALGEBRA-HARDEN. Energies and residuals in Ha.
use dft_mpi::{Collective, simulate};
use dft_numerics::real::RealMatrix;
use dft_solvers::distributed_davidson::{
    DistributedKeptStart, DistributedRealOperator, DistributedStart, DistributedTuning,
};
use dft_solvers::{BandTolerances, RealDavidson, SolverError};
use std::time::Duration;

struct Diagonal {
    n: usize,
    rows: Vec<usize>,
    fault: u8,
    calls: usize,
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
        self.calls += 1;
        // Complete the operator's own collective before the injected failure.
        let mut sum = [0.0];
        comm.tree_sum_f64(&[1.0], &mut sum).unwrap();
        for j in 0..x.cols() {
            for (i, &r) in self.rows.iter().enumerate() {
                y[(i, j)] = (r + 1) as f64 * x[(i, j)];
            }
        }
        if comm.rank() == 1 {
            if self.fault == 1 || (self.fault == 4 && self.calls == 2) {
                return Err(SolverError::Application("injected failure".into()));
            }
            if (self.fault == 2 || (self.fault == 5 && self.calls == 2)) && !y.as_slice().is_empty()
            {
                y[(0, 0)] = f64::NAN;
            }
        }
        Ok(())
    }
}
fn run(
    ranks: usize,
    n: usize,
    bands: usize,
    start: RealMatrix,
    fault: u8,
) -> Vec<Result<dft_solvers::real_davidson::RealEigenSolution, SolverError>> {
    simulate(ranks, |comm| {
        let mut op = Diagonal {
            n,
            rows: (0..n).filter(|r| r % ranks == comm.rank()).collect(),
            fault,
            calls: 0,
        };
        let diag: Vec<_> = op.rows.iter().map(|r| (r + 1) as f64).collect();
        let mut local = RealMatrix::zeros(op.rows.len(), start.cols()).unwrap();
        for j in 0..start.cols() {
            for (i, &r) in op.rows.iter().enumerate() {
                local[(i, j)] = start[(r, j)];
            }
        }
        if fault == 3 && comm.rank() == 1 {
            local = RealMatrix::zeros(0, start.cols()).unwrap();
        }

        RealDavidson {
            bands,
            tolerance: 1e-13,
            max_iterations: 100,
            max_subspace: n,
            seed: 42,
        }
        .solve_distributed_rows_local(
            comm,
            &mut op,
            Some(DistributedStart::Local(&local)),
            &diag,
            BandTolerances {
                strict_bands: bands,
                buffer_tolerance: 1e-13,
            },
            true,
            false,
            6,
            DistributedTuning {
                ..DistributedTuning::from_env()
            },
            &mut RealMatrix::zeros(0, 0).unwrap(),
            &mut DistributedKeptStart::new(),
            true,
        )
    })
}
fn check(sol: &dft_solvers::real_davidson::RealEigenSolution) {
    let x = &sol.eigenvectors;
    for (j, &lambda) in sol.eigenvalues.iter().enumerate() {
        // Independent exact dense diagonal spectrum, including every requested band.
        assert!(
            (lambda - (j + 1) as f64).abs() < 1e-10,
            "eigenvalue {j}: {lambda}"
        );
        let residual = (0..x.rows())
            .map(|i| (((i + 1) as f64 - lambda) * x[(i, j)]).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!(residual < 1e-6, "true residual {j}: {residual:e}");
        for k in 0..x.cols() {
            let gram = (0..x.rows()).map(|i| x[(i, j)] * x[(i, k)]).sum::<f64>();
            assert!(
                (gram - if j == k { 1.0 } else { 0.0 }).abs() <= 1e-12,
                "Gram {j},{k}: {gram:e}"
            );
        }
    }
}
#[test]
fn conditioning_chain_and_empty_partition() {
    let n = 8;
    let mut start = RealMatrix::zeros(n, n).unwrap();
    start[(0, 0)] = 1.0;
    for j in 1..n {
        start[(j - 1, j)] = 0.99_f64.sqrt();
        start[(j, j)] = 0.1;
    }
    for ranks in 1..=4 {
        for sol in run(ranks, n, n, start.clone(), 0) {
            check(&sol.unwrap());
        }
    }
    let mut unit = RealMatrix::zeros(2, 2).unwrap();
    unit[(0, 0)] = 1.0;
    unit[(1, 1)] = 1.0;
    for sol in run(4, 2, 2, unit, 0) {
        check(&sol.unwrap());
    }
}
#[test]
fn scaled_partial_start_is_normalized_before_completion() {
    for ranks in 1..=4 {
        let mut normal = RealMatrix::zeros(4, 1).unwrap();
        normal[(0, 0)] = 1.0;
        let mut huge = normal.clone();
        huge[(0, 0)] = 1e40;
        let reference = run(ranks, 4, 2, normal, 0);
        let out = run(ranks, 4, 2, huge, 0);
        for (a, b) in out.into_iter().zip(reference) {
            let (a, b) = (a.unwrap(), b.unwrap());
            check(&a);
            for (a, b) in a.eigenvalues.iter().zip(b.eigenvalues) {
                assert!((a - b).abs() < 1e-12);
            }
        }
    }
}
fn rank_local_failure(fault: u8) {
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let mut start = RealMatrix::zeros(4, 2).unwrap();
        start[(0, 0)] = 1.0;
        start[(1, 1)] = 1.0;
        if fault >= 4 {
            start[(2, 0)] = 0.2;
            start[(3, 1)] = 0.3;
        }
        let out = run(2, 4, 2, start, fault);
        tx.send(
            out.into_iter()
                .map(|r| r.err().map(|e| e.to_string()))
                .collect::<Vec<_>>(),
        )
        .unwrap();
    });
    let errors = rx
        .recv_timeout(Duration::from_secs(3))
        .expect("ranks hung after a local error");
    assert!(errors[0].is_some());
    assert_eq!(errors[0], errors[1]);
}
#[test]
fn rank_local_error_agrees() {
    rank_local_failure(1);
}
#[test]
fn rank_local_nan_agrees() {
    rank_local_failure(2);
}
#[test]
fn rank_local_bad_start_shape_agrees() {
    rank_local_failure(3);
}

#[test]
fn empty_partition_can_repair_rank_deficient_start() {
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let mut start = RealMatrix::zeros(2, 2).unwrap();
        start[(0, 0)] = 1.0;
        start[(0, 1)] = 1.0;
        tx.send(run(4, 2, 2, start, 0)).unwrap();
    });
    for result in rx
        .recv_timeout(Duration::from_secs(3))
        .expect("empty partition hung during Gram-Schmidt")
    {
        check(&result.unwrap());
    }
}

#[test]
fn rank_local_error_after_append_agrees() {
    rank_local_failure(4);
}
#[test]
fn rank_local_nan_after_append_agrees() {
    rank_local_failure(5);
}

impl dft_solvers::RealLinearOperator for Diagonal {
    fn dimension(&self) -> usize {
        self.n
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), SolverError> {
        for (i, (x, y)) in x.iter().zip(y).enumerate() {
            *y = (i + 1) as f64 * x;
        }
        Ok(())
    }
}

#[test]
fn serial_conditioned_and_scaled_partial_starts_have_exact_spectrum() {
    for (n, bands, partial) in [(8, 8, false), (4, 2, true)] {
        let mut start = RealMatrix::zeros(n, if partial { 1 } else { n }).unwrap();
        start[(0, 0)] = if partial { 1e40 } else { 1.0 };
        if !partial {
            for j in 1..n {
                start[(j - 1, j)] = 0.99_f64.sqrt();
                start[(j, j)] = 0.1;
            }
        }
        let op = Diagonal {
            n,
            rows: (0..n).collect(),
            fault: 0,
            calls: 0,
        };
        let diagonal: Vec<_> = (1..=n).map(|i| i as f64).collect();
        let sol = RealDavidson {
            bands,
            tolerance: 1e-13,
            max_iterations: 100,
            max_subspace: n,
            seed: 42,
        }
        .solve_research(
            &op,
            Some(&start),
            &diagonal,
            BandTolerances {
                strict_bands: bands,
                buffer_tolerance: 1e-13,
            },
            true,
            false,
            6,
            &mut None,
            dft_solvers::DavidsonPreconditioner::Default,
            dft_solvers::real_davidson::ResearchOptions {
                ..dft_solvers::real_davidson::ResearchOptions::from_env()
            },
        )
        .unwrap();
        check(&sol);
    }
}
