use dft_numerics::linalg::{
    ComplexMatrix,
    MatrixOp::{Adjoint, Identity},
    gemm, norm,
};
use dft_solvers::*;
use num_complex::Complex64;
use std::cell::Cell;
use std::time::Instant;

struct Matrix {
    a: ComplexMatrix,
    calls: Cell<usize>,
}
impl LinearOperator for Matrix {
    fn dimension(&self) -> usize {
        self.a.rows()
    }
    fn apply(&self, x: &[Complex64], y: &mut [Complex64]) -> Result<(), SolverError> {
        self.calls.set(self.calls.get() + 1);
        for (i, v) in y.iter_mut().enumerate() {
            *v = x.iter().enumerate().map(|(j, x)| self.a[(i, j)] * x).sum();
        }
        Ok(())
    }
}
fn fixture(n: usize, kind: usize) -> Matrix {
    let mut u = ComplexMatrix::zeros(n, n).unwrap();
    for j in 0..n {
        for i in 0..n {
            u[(i, j)] = Complex64::from_polar(
                1.0 / (n as f64).sqrt(),
                std::f64::consts::TAU * (i * j) as f64 / n as f64,
            );
        }
    }
    let mut ud = u.clone();
    for j in 0..n {
        let e = match kind {
            0 => j as f64,
            1 => (j / 3) as f64,
            2 => (j / 3) as f64 + (j % 3) as f64 * 1e-7,
            _ => (j as f64).powi(2),
        };
        for z in ud.column_mut(j) {
            *z *= e;
        }
    }
    Matrix {
        a: gemm(&ud, Identity, &u, Adjoint).unwrap(),
        calls: Cell::new(0),
    }
}
fn config(n: usize, m: usize) -> Davidson {
    Davidson {
        bands: m,
        tolerance: 1e-9,
        max_iterations: 800,
        max_subspace: n.min(6 * m),
        block_size: m,
        seed: 2026,
    }
}
fn gram(v: &ComplexMatrix) -> f64 {
    let g = gemm(v, Adjoint, v, Identity).unwrap();
    (0..g.cols())
        .flat_map(|j| (0..g.rows()).map(move |i| (i, j)))
        .map(|(i, j)| (g[(i, j)] - if i == j { 1.0 } else { 0.0 }).norm())
        .fold(0.0, f64::max)
}
fn check(op: &Matrix, solution: &EigenSolution, exact: &EigenSolution, tol: &[f64]) -> f64 {
    let mut y = vec![Complex64::default(); op.dimension()];
    for (j, t) in tol.iter().enumerate() {
        op.apply(solution.eigenvectors.column(j), &mut y).unwrap();
        for (i, z) in y.iter_mut().enumerate() {
            *z -= solution.eigenvectors[(i, j)] * solution.eigenvalues[j];
        }
        assert!(norm(&y) <= *t, "band {j}: {} > {t}", norm(&y));
        assert_eq!(norm(&y), solution.residuals[j]);
    }
    let error = solution
        .eigenvalues
        .iter()
        .zip(&exact.eigenvalues)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    assert!(error < 1e-8, "eigenvalue error {error}");
    assert!(gram(&solution.eigenvectors) < 1e-10);
    error
}
#[test]
fn spectra_and_partial_degenerate_boundary() {
    for kind in 0..4 {
        let op = fixture(24, kind);
        let cfg = config(24, 8); // boundary lies inside the threefold cluster
        let dense = DenseEigensolver.solve(&op, 8).unwrap();
        let mut warm = cfg.clone();
        warm.tolerance = 1e-2;
        warm.max_subspace = 16;
        let guess = warm.solve(&op, None, None).unwrap();
        let tol = vec![1e-9; 8];
        let solver = RmmDiis {
            max_sweeps: 200,
            ..Default::default()
        };
        let s = match solver.solve(&op, &guess.eigenvectors, None, &tol) {
            Ok(s) => s,
            Err(RmmError::Refinement {
                reason: RmmFailureReason::Stagnation,
                last: Some(last),
            }) if kind == 2 => {
                assert!(last.residuals.iter().any(|r| *r > 1e-9));
                assert!(gram(&last.eigenvectors) < 1e-10);
                println!(
                    "near-degenerate boundary: explicit stagnation at {:e}; testing Davidson recovery",
                    last.residuals.iter().copied().fold(0.0, f64::max)
                );
                let recovered = HybridRmm {
                    warmup_tolerance: 1.0,
                    ..Default::default()
                }
                .solve(&cfg, &op, Some(&guess.eigenvectors), None, &tol)
                .unwrap();
                check(&op, &recovered.solution, &dense, &tol);
                assert_eq!(recovered.fallback, Some(RmmFailureReason::Stagnation));
                continue;
            }
            Err(e) => panic!("{e}"),
        };
        let err = check(&op, &s, &dense, &tol);
        println!(
            "RMM spectrum={kind} error={err:e} residual={:e} orth={:e} sweeps={} calls={}",
            s.residuals.iter().copied().fold(0.0, f64::max),
            gram(&s.eigenvectors),
            s.iterations,
            s.operator_applications
        );
        let repeated = solver.solve(&op, &guess.eigenvectors, None, &tol).unwrap();
        assert_eq!(s, repeated);
    }
}
#[test]
fn per_band_tolerances_and_explicit_failures() {
    let op = fixture(20, 0);
    let dense = DenseEigensolver.solve(&op, 6).unwrap();
    let mut cfg = config(20, 6);
    cfg.tolerance = 1e-3;
    cfg.max_subspace = 12;
    let guess = cfg.solve(&op, None, None).unwrap();
    let tol = [1e-9, 1e-9, 1e-9, 1e-5, 1e-5, 1e-5];
    let s = RmmDiis::default()
        .solve(&op, &guess.eigenvectors, None, &tol)
        .unwrap();
    check(&op, &s, &dense, &tol);
    let failure = RmmDiis {
        max_sweeps: 1,
        ..Default::default()
    }
    .solve(&op, &guess.eigenvectors, None, &[1e-14; 6])
    .unwrap_err();
    assert!(
        matches!(failure,RmmError::Refinement { reason:RmmFailureReason::IterationLimit,last:Some(ref s) } if s.residuals.len()==6)
    );
    let duplicate = ComplexMatrix::new(
        20,
        6,
        (0..6)
            .flat_map(|_| dense.eigenvectors.column(0).iter().copied())
            .collect(),
    )
    .unwrap();
    assert!(matches!(
        RmmDiis::default().solve(&op, &duplicate, None, &tol),
        Err(RmmError::Refinement {
            reason: RmmFailureReason::RankLoss,
            ..
        })
    ));
    assert!(
        RmmDiis::default()
            .solve(&op, &guess.eigenvectors, None, &[f64::NAN; 6])
            .is_err()
    );
    let wrong = ComplexMatrix::zeros(19, 6).unwrap();
    assert!(RmmDiis::default().solve(&op, &wrong, None, &tol).is_err());
}
#[test]
fn hybrid_fallback_counts_all_applications() {
    let op = fixture(24, 2);
    let cfg = config(24, 8);
    let dense = DenseEigensolver.solve(&op, 8).unwrap();
    let mut warm = cfg.clone();
    warm.tolerance = 1e-2;
    warm.max_subspace = 16;
    let guess = warm.solve(&op, None, None).unwrap();
    let policy = HybridRmm {
        warmup_tolerance: 1.0,
        refinement: RmmDiis {
            max_sweeps: 1,
            ..Default::default()
        },
    };
    op.calls.set(0);
    let s = policy
        .solve(&cfg, &op, Some(&guess.eigenvectors), None, &[1e-9; 8])
        .unwrap();
    assert_eq!(s.fallback, Some(RmmFailureReason::IterationLimit));
    assert_eq!(s.solution.operator_applications, op.calls.get());
    assert_eq!(s.solution.history.len(), s.solution.iterations);
    check(&op, &s.solution, &dense, &[1e-9; 8]);
    op.calls.set(0);
    let cold = policy.solve(&cfg, &op, None, None, &[1e-9; 8]).unwrap();
    assert!(cold.warmed_up);
    assert_eq!(cold.solution.operator_applications, op.calls.get());
}
#[test]
fn operator_errors_are_not_hidden_by_fallback() {
    struct Bad;
    impl LinearOperator for Bad {
        fn dimension(&self) -> usize {
            4
        }
        fn apply(&self, _: &[Complex64], _: &mut [Complex64]) -> Result<(), SolverError> {
            Err(SolverError::Application("injected".into()))
        }
    }
    assert!(matches!(
        HybridRmm::default().solve(&config(4, 2), &Bad, None, None, &[1e-9; 2]),
        Err(RmmError::Solver(SolverError::Application(_)))
    ));
}

/// Reproduce: cargo test --release -p dft-solvers --test rmm benchmark -- --ignored --nocapture
#[test]
#[ignore = "release synthetic band-count/guess-quality benchmark"]
fn benchmark() {
    let bands: Vec<usize> = std::env::var("RMM_BANDS")
        .unwrap_or("8,16,32,64,128".into())
        .split(',')
        .map(|s| s.parse().unwrap())
        .collect();
    let repeats: usize = std::env::var("RMM_REPEATS")
        .unwrap_or("3".into())
        .parse()
        .unwrap();
    println!(
        "kind,bands,quality,repeat,solver,seconds,applications,sweeps,fallback,eigen_error,residual,orth"
    );
    for m in bands {
        for kind in 0..4 {
            let n = 2 * m + 8;
            let op = fixture(n, kind);
            let cfg = config(n, m);
            let exact = DenseEigensolver.solve(&op, m).unwrap();
            // Controlled oracle-contaminated guess: no narrow-subspace Davidson
            // pre-run, whose own failures would censor the benchmark matrix.
            let mut guess = exact.eigenvectors.clone();
            for j in 0..m {
                for i in 0..n {
                    guess[(i, j)] += Complex64::from_polar(
                        1e-3 / (n * n) as f64 / (n as f64).sqrt(),
                        std::f64::consts::TAU * (i * (j + m)) as f64 / n as f64,
                    );
                }
            }
            let guess = dft_numerics::linalg::orthonormalize(&guess, None)
                .unwrap()
                .vectors;
            for quality in ["cold", "warm"] {
                for rep in 0..repeats {
                    for hybrid in if rep % 2 == 0 {
                        [false, true]
                    } else {
                        [true, false]
                    } {
                        let selected = std::env::var("RMM_SOLVER").unwrap_or_default();
                        if (selected == "hybrid" && !hybrid) || (selected == "davidson" && hybrid) {
                            continue;
                        }
                        if selected == "rmm" && (!hybrid || quality == "cold") {
                            continue;
                        }
                        let initial = if quality == "warm" {
                            Some(&guess)
                        } else {
                            None
                        };
                        let start = Instant::now();
                        let (s, fallback, sweeps) = if selected == "rmm" {
                            match RmmDiis::default().solve(
                                &op,
                                &guess,
                                None,
                                &vec![cfg.tolerance; m],
                            ) {
                                Ok(s) => {
                                    let sweeps = s.iterations;
                                    (s, None, sweeps)
                                }
                                Err(RmmError::Refinement { reason, last }) => {
                                    println!(
                                        "RMM_FAILURE,{kind},{m},{quality},{rep},{reason:?},{}",
                                        last.as_ref().map_or(f64::INFINITY, |s| s
                                            .residuals
                                            .iter()
                                            .copied()
                                            .fold(0.0, f64::max))
                                    );
                                    continue;
                                }
                                Err(e) => panic!("unexpected RMM failure: {e}"),
                            }
                        } else if hybrid {
                            let r = HybridRmm::default()
                                .solve(&cfg, &op, initial, None, &vec![cfg.tolerance; m])
                                .unwrap();
                            (r.solution, r.fallback, r.rmm_sweeps)
                        } else {
                            (cfg.solve(&op, initial, None).unwrap(), None, 0)
                        };
                        let elapsed = start.elapsed().as_secs_f64();
                        let error = check(&op, &s, &exact, &vec![cfg.tolerance; m]);
                        if std::env::var_os("RMM_HISTORY").is_some() {
                            println!("HISTORY,{kind},{m},{quality},{hybrid},{:?}", s.history);
                        }
                        println!(
                            "{kind},{m},{quality},{rep},{},{elapsed:.9},{},{sweeps},{fallback:?},{error:e},{:e},{:e}",
                            if selected == "rmm" {
                                "rmm"
                            } else if hybrid {
                                "hybrid"
                            } else {
                                "davidson"
                            },
                            s.operator_applications,
                            s.residuals.iter().copied().fold(0.0, f64::max),
                            gram(&s.eigenvectors)
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn malformed_overrides_and_stagnation_are_explicit() {
    struct Override {
        shape: bool,
    }
    impl LinearOperator for Override {
        fn dimension(&self) -> usize {
            4
        }
        fn apply(&self, _: &[Complex64], _: &mut [Complex64]) -> Result<(), SolverError> {
            unreachable!()
        }
        fn apply_block(&self, _: &ComplexMatrix, y: &mut ComplexMatrix) -> Result<(), SolverError> {
            if self.shape {
                *y = ComplexMatrix::zeros(3, 2)?;
            } else {
                y[(0, 0)] = f64::NAN.into();
            }
            Ok(())
        }
    }
    let mut x = ComplexMatrix::zeros(4, 2).unwrap();
    x[(0, 0)] = 1.0.into();
    x[(1, 1)] = 1.0.into();
    for shape in [false, true] {
        assert!(matches!(
            RmmDiis::default().solve(&Override { shape }, &x, None, &[1e-9; 2]),
            Err(RmmError::Solver(SolverError::Linalg(_)))
        ));
    }
    struct Zero;
    impl Preconditioner for Zero {
        fn apply(
            &self,
            _: &[f64],
            _: &ComplexMatrix,
            out: &mut ComplexMatrix,
        ) -> Result<(), SolverError> {
            for j in 0..out.cols() {
                out.column_mut(j).fill(Complex64::default());
            }
            Ok(())
        }
    }
    let op = fixture(12, 0);
    let mut cfg = config(12, 4);
    cfg.max_subspace = 8;
    cfg.tolerance = 0.01;
    let guess = cfg.solve(&op, None, None).unwrap();
    let result = RmmDiis {
        stagnation_sweeps: 2,
        ..Default::default()
    }
    .solve(&op, &guess.eigenvectors, Some(&Zero), &[1e-12; 4]);
    assert!(
        matches!(result,Err(RmmError::Refinement{reason:RmmFailureReason::Stagnation,last:Some(s)}) if s.residuals.iter().any(|r|*r>1e-12))
    );
}

#[test]
fn loose_buffer_can_stop_without_relaxing_strict_bands() {
    let mut a = ComplexMatrix::zeros(6, 6).unwrap();
    for i in 0..6 {
        a[(i, i)] = (i as f64).into();
    }
    let op = Matrix {
        a,
        calls: Cell::new(0),
    };
    let mut x = ComplexMatrix::zeros(6, 3).unwrap();
    x[(0, 0)] = 1.0.into();
    x[(1, 1)] = 1.0.into();
    x[(2, 2)] = 1.0.into();
    x[(3, 2)] = 1e-4.into();
    let solver = RmmDiis {
        max_sweeps: 1,
        ..Default::default()
    };
    let s = solver.solve(&op, &x, None, &[1e-12, 1e-12, 1e-3]).unwrap();
    assert_eq!(s.iterations, 1);
    assert!(s.residuals[0] <= 1e-12 && s.residuals[1] <= 1e-12);
    assert!(s.residuals[2] > 1e-12 && s.residuals[2] < 1e-3);
    assert!(matches!(
        solver.solve(&op, &x, None, &[1e-12; 3]),
        Err(RmmError::Refinement {
            reason: RmmFailureReason::IterationLimit,
            ..
        })
    ));
}

#[test]
fn duplicate_trial_vectors_trigger_safe_davidson_fallback() {
    struct Collapse {
        initial: ComplexMatrix,
        calls: Cell<usize>,
    }
    impl Preconditioner for Collapse {
        fn apply(
            &self,
            _: &[f64],
            r: &ComplexMatrix,
            out: &mut ComplexMatrix,
        ) -> Result<(), SolverError> {
            let first = self.calls.get() == 0;
            self.calls.set(self.calls.get() + 1);
            for j in 0..out.cols() {
                for i in 0..out.rows() {
                    out[(i, j)] = if first {
                        Complex64::new(if i == 1 { 1.0 } else { 0.0 }, 0.0) - self.initial[(i, j)]
                    } else {
                        -r[(i, j)]
                    };
                }
            }
            Ok(())
        }
    }
    let mut a = ComplexMatrix::zeros(4, 4).unwrap();
    for i in 0..4 {
        a[(i, i)] = (i as f64).into();
    }
    let op = Matrix {
        a,
        calls: Cell::new(0),
    };
    let mut x = ComplexMatrix::zeros(4, 2).unwrap();
    x[(0, 0)] = (0.5_f64.sqrt()).into();
    x[(2, 0)] = (0.5_f64.sqrt()).into();
    x[(1, 1)] = 1.0.into();
    let pre = Collapse {
        initial: x.clone(),
        calls: Cell::new(0),
    };
    let r = RmmDiis::default().solve(&op, &x, Some(&pre), &[1e-9; 2]);
    assert!(matches!(
        r,
        Err(RmmError::Refinement {
            reason: RmmFailureReason::RankLoss,
            last: Some(_)
        })
    ));
    pre.calls.set(0);
    let policy = HybridRmm {
        warmup_tolerance: 2.0,
        ..Default::default()
    };
    let s = policy
        .solve(&config(4, 2), &op, Some(&x), Some(&pre), &[1e-9; 2])
        .unwrap();
    assert_eq!(s.fallback, Some(RmmFailureReason::RankLoss));
    let exact = DenseEigensolver.solve(&op, 2).unwrap();
    check(&op, &s.solution, &exact, &[1e-9; 2]);
}

#[test]
fn terminal_hybrid_failures_identify_phase_and_keep_residuals() {
    let op = fixture(24, 0);
    let mut cfg = config(24, 4);
    cfg.max_iterations = 1;
    let error = HybridRmm::default()
        .solve(&cfg, &op, None, None, &[1e-12; 4])
        .unwrap_err();
    assert!(
        matches!(error,RmmError::Davidson{phase:HybridPhase::Warmup,operator_applications,source:SolverError::NotConverged{residuals,..}} if operator_applications>0 && residuals.len()==4)
    );
    cfg.max_iterations = 800;
    cfg.max_subspace = 12;
    cfg.tolerance = 1e-3;
    let initial = cfg.solve(&op, None, None).unwrap();
    cfg.max_iterations = 1;
    let policy = HybridRmm {
        warmup_tolerance: 1.0,
        refinement: RmmDiis {
            max_sweeps: 1,
            ..Default::default()
        },
    };
    let error = policy
        .solve(&cfg, &op, Some(&initial.eigenvectors), None, &[1e-12; 4])
        .unwrap_err();
    assert!(
        matches!(error,RmmError::Davidson{phase:HybridPhase::Fallback,operator_applications,source:SolverError::NotConverged{residuals,..}} if operator_applications>8 && residuals.len()==4)
    );
}
