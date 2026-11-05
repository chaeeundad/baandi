use super::*;
struct MatrixOperator(ComplexMatrix);
impl LinearOperator for MatrixOperator {
    fn dimension(&self) -> usize {
        self.0.rows()
    }
    fn apply(&self, x: &[Complex64], y: &mut [Complex64]) -> Result<(), SolverError> {
        if x.len() != self.dimension() || y.len() != self.dimension() {
            return Err(LinalgError::Dimension.into());
        }
        for (i, y) in y.iter_mut().enumerate() {
            *y = x.iter().enumerate().map(|(j, x)| self.0[(i, j)] * x).sum();
        }
        Ok(())
    }
    fn apply_block(&self, x: &ComplexMatrix, y: &mut ComplexMatrix) -> Result<(), SolverError> {
        *y = gemm(&self.0, Identity, x, Identity)?;
        Ok(())
    }
}
struct Diagonal(Vec<f64>);
impl LinearOperator for Diagonal {
    fn dimension(&self) -> usize {
        self.0.len()
    }
    fn apply(&self, x: &[Complex64], y: &mut [Complex64]) -> Result<(), SolverError> {
        for ((y, x), d) in y.iter_mut().zip(x).zip(&self.0) {
            *y = x * d;
        }
        Ok(())
    }
}
fn config(bands: usize, n: usize) -> Davidson {
    Davidson {
        bands,
        tolerance: 1e-11,
        max_iterations: 500,
        max_subspace: n,
        block_size: bands,
        seed: 42,
    }
}
fn random_matrix(n: usize, seed: u64) -> ComplexMatrix {
    let mut rng = SplitMix64(seed);
    let mut a = ComplexMatrix::zeros(n, n).unwrap();
    for j in 0..n {
        for i in 0..=j {
            let z = Complex64::new(rng.next(), if i == j { 0.0 } else { rng.next() });
            a[(i, j)] = z;
            a[(j, i)] = z.conj();
        }
    }
    a
}
fn fourier(d: &[f64]) -> ComplexMatrix {
    let n = d.len();
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
    for (j, d) in d.iter().enumerate() {
        for z in ud.column_mut(j) {
            *z *= d;
        }
    }
    gemm(&ud, Identity, &u, Adjoint).unwrap()
}
fn assess(label: &str, op: &dyn LinearOperator, s: &EigenSolution, expected: &[f64], tol: f64) {
    let error = s
        .eigenvalues
        .iter()
        .zip(expected)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    let (_, r) = residual(op, &s.eigenvectors, &s.eigenvalues).unwrap();
    let maximum = max_residual(&r);
    let gram = gemm(&s.eigenvectors, Adjoint, &s.eigenvectors, Identity).unwrap();
    let mut orth: f64 = 0.0;
    for j in 0..gram.cols() {
        for i in 0..gram.rows() {
            orth = orth.max((gram[(i, j)] - if i == j { 1.0 } else { 0.0 }).norm());
        }
    }
    println!(
        "{label}: eigenvalue_error={error:e} residual={maximum:e} orth={orth:e} iterations={} (dimensionless)",
        s.iterations
    );
    assert!(error < 1e-9);
    assert!(maximum < tol);
    assert!(orth < 1e-10);
    assert_eq!(r, s.residuals);
    assert_eq!(s.history.len(), s.iterations);
    assert!(s.eigenvalues.windows(2).all(|v| v[0] <= v[1]));
}
#[test]
fn small_analytic_and_diagonal() {
    let two = ComplexMatrix::new(
        2,
        2,
        vec![
            2.0.into(),
            Complex64::new(1.0, 1.0),
            Complex64::new(1.0, -1.0),
            3.0.into(),
        ],
    )
    .unwrap();
    for (a, expected) in [
        (two, vec![1.0, 4.0]),
        (fourier(&[-2.0, 0.5, 3.0, 7.0]), vec![-2.0, 0.5, 3.0, 7.0]),
    ] {
        let op = MatrixOperator(a);
        let n = op.dimension();
        for bands in [1, n] {
            let dense = DenseEigensolver.solve(&op, bands).unwrap();
            let iterative = config(bands, n).solve(&op, None, None).unwrap();
            for s in [&dense, &iterative] {
                assess("analytic", &op, s, &expected[..bands], 1e-12);
                assert!(
                    s.eigenvalues
                        .iter()
                        .zip(&expected)
                        .all(|(a, b)| (a - b).abs() < 1e-12)
                );
            }
        }
    }
    let op = Diagonal((0..20).rev().map(|i| i as f64).collect());
    let s = config(4, 20).solve(&op, None, None).unwrap();
    assess("diagonal", &op, &s, &[0.0, 1.0, 2.0, 3.0], 1e-11);
}
#[test]
fn degeneracy_and_near_degeneracy() {
    for spectrum in [
        vec![-1.0, -1.0, -1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        vec![-1.0, -1.0 + 1e-6, 0.0, 2.0, 3.0, 4.0, 5.0, 6.0],
    ] {
        let op = MatrixOperator(fourier(&spectrum));
        for bands in [2, 3] {
            let s = config(bands, 8).solve(&op, None, None).unwrap();
            assess("degenerate", &op, &s, &spectrum[..bands], 1e-11);
        }
    }
}
#[test]
fn restart_initial_guess_and_determinism() {
    let op = MatrixOperator(random_matrix(18, 83));
    let dense = DenseEigensolver.solve(&op, 3).unwrap();
    let mut cfg = config(3, 9);
    cfg.tolerance = 1e-10;
    let mut initial = ComplexMatrix::zeros(18, 3).unwrap();
    initial[(0, 0)] = 1.0.into();
    initial[(0, 1)] = 1.0.into();
    let a = cfg.solve(&op, Some(&initial), None).unwrap();
    assess("restart", &op, &a, &dense.eigenvalues, cfg.tolerance);
    assert!(a.iterations > 3); // Capacity forces a restart after at most three expansions.
    let b = cfg.solve(&op, Some(&initial), None).unwrap();
    assert_eq!(a, b); // Includes vectors, residuals and history, not just eigenvalues.
    for (a, b) in a.eigenvalues.iter().zip(&b.eigenvalues) {
        assert_eq!(a.to_bits(), b.to_bits());
    }
    for (x, y) in a
        .eigenvectors
        .as_slice()
        .iter()
        .zip(b.eigenvectors.as_slice())
    {
        assert_eq!(x.re.to_bits(), y.re.to_bits());
        assert_eq!(x.im.to_bits(), y.im.to_bits());
    }
    for (x, y) in a
        .residuals
        .iter()
        .chain(&a.history)
        .zip(b.residuals.iter().chain(&b.history))
    {
        assert_eq!(x.to_bits(), y.to_bits());
    }
    cfg.seed = 91;
    let c = cfg.solve(&op, Some(&initial), None).unwrap();
    assess("different seed", &op, &c, &a.eigenvalues, cfg.tolerance);
}
#[test]
fn explicit_failure_and_invalid_inputs() {
    let op = MatrixOperator(random_matrix(15, 12));
    let mut cfg = config(2, 8);
    cfg.max_iterations = 1;
    match cfg.solve(&op, None, None) {
        Err(SolverError::NotConverged {
            iterations,
            max_residual,
            residuals,
        }) => {
            println!(
                "expected failure iterations={iterations} max_residual={max_residual:e} residuals={residuals:?} (dimensionless)"
            );
            assert_eq!(iterations, 1);
            assert_eq!(residuals.len(), 2);
            assert_eq!(max_residual, super::max_residual(&residuals));
            assert!(max_residual > cfg.tolerance);
        }
        other => panic!("expected NotConverged, got {other:?}"),
    }
    cfg.max_iterations = 0;
    assert!(matches!(
        cfg.solve(&op, None, None),
        Err(SolverError::InvalidParameter(_))
    ));
    cfg.max_iterations = 5;
    cfg.max_subspace = 2;
    assert!(cfg.solve(&op, None, None).is_err());
    assert!(matches!(
        DenseEigensolver.solve(&Diagonal(vec![0.0; 1025]), 1),
        Err(SolverError::DenseLimit { .. })
    ));
    let bad = MatrixOperator(
        ComplexMatrix::new(2, 2, vec![0.0.into(), 1.0.into(), 0.0.into(), 0.0.into()]).unwrap(),
    );
    assert!(matches!(
        DenseEigensolver.solve(&bad, 1),
        Err(SolverError::Linalg(LinalgError::NotHermitian { .. }))
    ));
    assert!(DiagonalPreconditioner::new(vec![0.0], 0.0).is_err());
    let p = DiagonalPreconditioner::new(vec![1.0, 1.0 - 1e-8, 1.0 + 1e-8], 0.1).unwrap();
    let r = ComplexMatrix::new(3, 1, vec![1.0.into(); 3]).unwrap();
    let mut out = r.clone();
    p.apply(&[1.0], &r, &mut out).unwrap();
    assert_eq!(
        out.as_slice(),
        &[(-10.0).into(), 10.0.into(), (-10.0).into()]
    );
    struct Broken;
    impl LinearOperator for Broken {
        fn dimension(&self) -> usize {
            2
        }
        fn apply(&self, _: &[Complex64], _: &mut [Complex64]) -> Result<(), SolverError> {
            Err(SolverError::Application("test failure".into()))
        }
    }
    assert!(matches!(
        config(1, 2).solve(&Broken, None, None),
        Err(SolverError::Application(_))
    ));
    let nan = Diagonal(vec![f64::NAN; 2]);
    assert!(matches!(
        config(1, 2).solve(&nan, None, None),
        Err(SolverError::Linalg(LinalgError::NonFinite))
    ));
}
/// Run: cargo test --offline --release -p dft-solvers -- --ignored --nocapture
#[test]
#[ignore = "N=200 reference comparison; release numerical acceptance test"]
fn random_two_hundred() {
    let op = MatrixOperator(random_matrix(200, 942));
    let mut rng = SplitMix64(9);
    let x = rng.block(200, 1).unwrap();
    let y = rng.block(200, 1).unwrap();
    let hx = apply(&op, &x).unwrap();
    let hy = apply(&op, &y).unwrap();
    let lhs: Complex64 = x
        .as_slice()
        .iter()
        .zip(hy.as_slice())
        .map(|(x, hy)| x.conj() * hy)
        .sum();
    let rhs: Complex64 = hx
        .as_slice()
        .iter()
        .zip(y.as_slice())
        .map(|(hx, y)| hx.conj() * y)
        .sum();
    let defect = (lhs - rhs).norm() / lhs.norm().max(rhs.norm());
    println!("N=200 operator Hermiticity relative={defect:e} (dimensionless)");
    assert!(defect < 1e-12);
    let dense = DenseEigensolver.solve(&op, 8).unwrap();
    let mut cfg = config(8, 80);
    cfg.tolerance = 1e-10;
    let s = cfg.solve(&op, None, None).unwrap();
    assess(
        "N=200 Davidson vs dense",
        &op,
        &s,
        &dense.eigenvalues,
        cfg.tolerance,
    );
}
/// Run: cargo test --offline --release -p dft-solvers -- --ignored --nocapture
#[test]
#[ignore = "N=300 preconditioner comparison; release numerical acceptance test"]
fn kinetic_three_hundred() {
    let n = 300;
    let diagonal: Vec<_> = (0..n).map(|k| 0.5 * (k * k) as f64).collect();
    let mut a = random_matrix(n, 932);
    for j in 0..n {
        for i in 0..n {
            a[(i, j)] *= 0.02;
        }
        a[(j, j)] += diagonal[j];
    }
    let op = MatrixOperator(a);
    let dense = DenseEigensolver.solve(&op, 4).unwrap();
    let p = DiagonalPreconditioner::new(diagonal, 0.1).unwrap();
    let mut cfg = config(4, 160);
    cfg.tolerance = 1e-9;
    cfg.max_iterations = 3000;
    let pre = cfg.solve(&op, None, Some(&p)).unwrap();
    let plain = cfg.solve(&op, None, None).unwrap();
    assess(
        "N=300 preconditioned",
        &op,
        &pre,
        &dense.eigenvalues,
        cfg.tolerance,
    );
    assess(
        "N=300 unpreconditioned",
        &op,
        &plain,
        &dense.eigenvalues,
        cfg.tolerance,
    );
    println!(
        "N=300 iteration comparison: preconditioned={} unpreconditioned={}",
        pre.iterations, plain.iterations
    );
    assert!(pre.iterations < plain.iterations);
}

#[test]
fn dependent_corrections_and_shape_errors() {
    struct ZeroCorrections;
    impl Preconditioner for ZeroCorrections {
        fn apply(
            &self,
            _: &[f64],
            _: &ComplexMatrix,
            _: &mut ComplexMatrix,
        ) -> Result<(), SolverError> {
            Ok(())
        }
    }
    let op = MatrixOperator(random_matrix(8, 11));
    let expected = DenseEigensolver.solve(&op, 2).unwrap();
    let result = config(2, 8)
        .solve(&op, None, Some(&ZeroCorrections))
        .unwrap();
    assess(
        "dependent correction fallback",
        &op,
        &result,
        &expected.eigenvalues,
        1e-11,
    );
    struct WrongShape;
    impl LinearOperator for WrongShape {
        fn dimension(&self) -> usize {
            2
        }
        fn apply(&self, _: &[Complex64], _: &mut [Complex64]) -> Result<(), SolverError> {
            unreachable!()
        }
        fn apply_block(&self, _: &ComplexMatrix, y: &mut ComplexMatrix) -> Result<(), SolverError> {
            *y = ComplexMatrix::zeros(0, 0)?;
            Ok(())
        }
    }
    assert!(matches!(
        config(1, 2).solve(&WrongShape, None, None),
        Err(SolverError::Linalg(LinalgError::Dimension))
    ));
    let op = Diagonal(vec![3.0]);
    let result = config(1, 1).solve(&op, None, None).unwrap();
    assess("scalar", &op, &result, &[3.0], 1e-12);
}

#[path = "product_tests.rs"]
mod product_tests;

#[test]
fn contracted_restarts_match_oracle_and_degeneracies() {
    for (name, matrix) in [
        ("random", random_matrix(18, 83)),
        (
            "degenerate",
            fourier(&[-1., -1., -1., 2., 3., 4., 5., 6., 7., 8., 9., 10.]),
        ),
    ] {
        let op = MatrixOperator(matrix);
        let v = orthonormalize(&SplitMix64(71).block(op.dimension(), 6).unwrap(), None)
            .unwrap()
            .vectors;
        let hv = apply(&op, &v).unwrap();
        let projected = gemm(&v, Adjoint, &hv, Identity).unwrap();
        let (_, coefficients) = hermitian_eigh(&projected).unwrap();
        let coefficients = first_columns(&coefficients, 3).unwrap();
        let ritz = gemm(&v, Identity, &coefficients, Identity).unwrap();
        let true_products = apply(&op, &ritz).unwrap();
        let mut cached = ComplexMatrix::zeros(op.dimension(), 3).unwrap();
        let gram = orthogonality_defect(&ritz);
        let mut work = GemmWorkspace::new(op.dimension(), 6, 1).unwrap();
        let (_, valid) = restart(ritz, &hv, &coefficients, &mut cached, true, &mut work).unwrap();
        assert!(valid);
        let error = cached
            .as_slice()
            .iter()
            .zip(true_products.as_slice())
            .map(|(a, b)| (a - b).norm())
            .fold(0.0, f64::max);
        assert!(error < 1e-14);
        println!("contraction {name}: product error {error:e}, restart Gram defect {gram:e}");
        let mut cfg = config(3, 9);
        cfg.tolerance = 1e-10;
        let old = cfg
            .solve_with_restart(&op, None, None, true, false)
            .unwrap();
        let new = cfg.solve(&op, None, None).unwrap();
        let dense = DenseEigensolver.solve(&op, 3).unwrap();
        assess(name, &op, &new, &dense.eigenvalues, cfg.tolerance);
        for (a, b) in old.eigenvalues.iter().zip(&new.eigenvalues) {
            assert!((a - b).abs() < 1e-12);
        }
        let defect = orthogonality_defect(&new.eigenvectors);
        assert!(defect < 1e-12);
        println!(
            "contraction {name}: applications {} -> {}, iterations {} -> {}, Gram defect {defect:e}",
            old.operator_applications, new.operator_applications, old.iterations, new.iterations
        );
        assert!(new.operator_applications < old.operator_applications);
    }
}

#[test]
fn contraction_guard_detects_loss_of_orthogonality() {
    let mut v = ComplexMatrix::zeros(4, 2).unwrap();
    v[(0, 0)] = 1.0.into();
    v[(1, 1)] = 1.0.into();
    assert_eq!(orthogonality_defect(&v), 0.0);
    v[(0, 1)] = 2e-12.into();
    assert!(orthogonality_defect(&v) > 1e-12);
    let mut coefficients = ComplexMatrix::zeros(2, 2).unwrap();
    coefficients[(0, 0)] = 1.0.into();
    coefficients[(1, 1)] = 1.0.into();
    let products = v.clone();
    let mut work = GemmWorkspace::new(4, 2, 1).unwrap();
    let mut cache = ComplexMatrix::zeros(4, 2).unwrap();
    let (repaired, valid) = restart(
        v.clone(),
        &products,
        &coefficients,
        &mut cache,
        true,
        &mut work,
    )
    .unwrap();
    assert!(!valid, "fallback must force a fresh application");
    assert!(orthogonality_defect(&repaired) < 1e-15);
    let (_, valid) = restart(
        repaired.clone(),
        &products,
        &coefficients,
        &mut cache,
        false,
        &mut work,
    )
    .unwrap();
    assert!(!valid, "untrusted products must not be contracted");
    let (_, valid) = restart(
        repaired,
        &products,
        &coefficients,
        &mut cache,
        true,
        &mut work,
    )
    .unwrap();
    assert!(valid);
    assert_eq!(cache, products);
    v[(0, 1)] = f64::NAN.into();
    assert_eq!(orthogonality_defect(&v), f64::INFINITY);
}

#[test]
fn subspace_multipliers_preserve_lowest_eigenpairs() {
    let op = MatrixOperator(random_matrix(18, 73));
    let dense = DenseEigensolver.solve(&op, 3).unwrap();
    for multiplier in [2, 3, 4, 6] {
        let mut cfg = config(3, 3 * multiplier);
        cfg.tolerance = 1e-9;
        let s = cfg.solve(&op, None, None).unwrap();
        assess("capacity", &op, &s, &dense.eigenvalues, cfg.tolerance);
    }
}

#[test]
fn degenerate_buffer_policy_and_iteration_limit() {
    let n = 40;
    let mut spectrum = vec![-3.0, -2.0, -1.0];
    spectrum.extend((3..n).map(|i| {
        if i < 9 {
            1.0
        } else {
            1.0 + 10f64.powf((i - 9) as f64 / 7.0 - 6.0)
        }
    }));
    let op = Diagonal(spectrum);
    let mut initial = ComplexMatrix::zeros(n, 4).unwrap();
    for j in 0..3 {
        initial[(j, j)] = Complex64::new(1.0, 0.0);
    }
    initial[(3, 3)] = Complex64::new(1.0, 0.0);
    for i in 4..n {
        initial[(i, 3)] = Complex64::new(1e-6, 0.0);
    }
    let cfg = Davidson {
        bands: 4,
        tolerance: 1e-12,
        max_iterations: 400,
        max_subspace: 6,
        block_size: 4,
        seed: 1,
    };
    let old = cfg.solve(&op, Some(&initial), None);
    println!(
        "cluster old: {:?}",
        old.as_ref().map(|s| (s.iterations, &s.residuals))
    );
    let new = cfg
        .solve_with_band_tolerances(
            &op,
            Some(&initial),
            None,
            BandTolerances {
                strict_bands: 3,
                buffer_tolerance: 1e-5,
            },
        )
        .unwrap();
    println!("cluster new: {} {:?}", new.iterations, new.residuals);
    assert!(matches!(
        old,
        Err(SolverError::NotConverged {
            iterations: 400,
            ..
        })
    ));
    assert_eq!(new.buffer, [false, false, false, true]);
    assert!(new.converged.iter().all(|&c| c));
    assert!(new.residuals[..3].iter().all(|&r| r <= cfg.tolerance));
    assert!(new.residuals[3] > cfg.tolerance && new.residuals[3] <= 1e-5);
    let (_, actual) = residual(&op, &new.eigenvectors, &new.eigenvalues).unwrap();
    assert_eq!(actual, new.residuals);
    assert!(new.residuals_certified);
    let intermediate = cfg
        .solve_with_band_tolerances_and_certification(
            &op,
            Some(&initial),
            None,
            BandTolerances {
                strict_bands: 3,
                buffer_tolerance: 1e-5,
            },
            false,
        )
        .unwrap();
    assert!(!intermediate.residuals_certified);
    assert_eq!(intermediate.eigenvectors, new.eigenvectors);
    assert_eq!(intermediate.eigenvalues, new.eigenvalues);
    assert_eq!(
        intermediate.operator_applications + cfg.bands,
        new.operator_applications
    );

    let limited = Davidson {
        max_iterations: 1,
        ..cfg.clone()
    };
    let result = limited
        .solve_with_band_tolerances(
            &op,
            Some(&initial),
            None,
            BandTolerances {
                strict_bands: 3,
                buffer_tolerance: cfg.tolerance,
            },
        )
        .unwrap();
    assert_eq!(result.converged, [true, true, true, false]);
    assert_eq!(result.iterations, 1);
    let limit_uncertified = limited
        .solve_with_band_tolerances_and_certification(
            &op,
            Some(&initial),
            None,
            BandTolerances {
                strict_bands: 3,
                buffer_tolerance: cfg.tolerance,
            },
            false,
        )
        .unwrap();
    assert!(limit_uncertified.residuals_certified);
    assert_eq!(limit_uncertified, result);

    assert!(matches!(
        limited.solve_with_band_tolerances(
            &op,
            Some(&initial),
            None,
            BandTolerances {
                strict_bands: 4,
                buffer_tolerance: cfg.tolerance
            }
        ),
        Err(SolverError::NotConverged { .. })
    ));
    for policy in [
        BandTolerances {
            strict_bands: 0,
            buffer_tolerance: 1e-5,
        },
        BandTolerances {
            strict_bands: 5,
            buffer_tolerance: 1e-5,
        },
        BandTolerances {
            strict_bands: 3,
            buffer_tolerance: f64::NAN,
        },
        BandTolerances {
            strict_bands: 3,
            buffer_tolerance: 1e-14,
        },
    ] {
        assert!(matches!(
            cfg.solve_with_band_tolerances(&op, None, None, policy),
            Err(SolverError::InvalidParameter(_))
        ));
    }
}

/// P15C mechanism fixture (Ha): a zero residual in an invariant start does not
/// count eigenvalues below its Ritz boundary. The diagonal preconditioner
/// preserves this subspace, so no correction can supply the missing direction.
#[test]
fn spectral_completeness_is_not_a_residual_certificate() {
    let op = Diagonal(vec![-3.0, -3.0, -1.0, 2.0, 4.0, 6.0, 8.0, 10.0]);
    let mut start = ComplexMatrix::zeros(8, 2).unwrap();
    start[(1, 0)] = 1.0.into();
    start[(2, 1)] = 1.0.into();
    let cfg = config(2, 8);
    let preconditioner = DiagonalPreconditioner::new(op.0.clone(), 0.1).unwrap();
    for pre in [None, Some(&preconditioner as &dyn Preconditioner)] {
        let missed = cfg.solve(&op, Some(&start), pre).unwrap();
        assert!(missed.residuals_certified);
        assert_eq!(missed.residuals, vec![0.0, 0.0]);
        assert_eq!(missed.eigenvalues, vec![-3.0, -1.0]);
        let dense = DenseEigensolver.solve(&op, 2).unwrap();
        assert!(dense.eigenvalues.iter().all(|v| (v + 3.0).abs() < 1e-9));
        assert!((missed.eigenvalues[1] - dense.eigenvalues[1] - 2.0).abs() < 1e-9);
        // Independent seeded start restores both members of the doublet.
        let recovered = cfg.solve(&op, None, pre).unwrap();
        assess(
            "independent spectral start",
            &op,
            &recovered,
            &dense.eigenvalues,
            1e-10,
        );
    }
}

#[test]
fn spectral_guard_recovers_missing_degeneracy_and_preserves_complete_state() {
    let op = Diagonal(vec![-3.0, -3.0, -1.0, 2.0, 4.0, 6.0, 8.0, 10.0]);
    let mut start = ComplexMatrix::zeros(8, 2).unwrap();
    start[(1, 0)] = 1.0.into();
    start[(2, 1)] = 1.0.into();
    let cfg = config(2, 8);
    let missed = cfg.solve(&op, Some(&start), None).unwrap();
    let mut invalid = missed.clone();
    invalid.eigenvalues.clear();
    assert!(matches!(
        cfg.guard_spectrum(&op, None, invalid),
        Err(SolverError::InvalidParameter(_))
    ));
    let fixed = cfg.guard_spectrum(&op, None, missed).unwrap();
    assert!(fixed.eigenvalues.iter().all(|v| (v + 3.0).abs() < 1e-9));
    let evidence = fixed.spectral_completeness.as_ref().unwrap();
    assert!(evidence.passed);
    assert_eq!(evidence.retries, 1);
    assert!(evidence.extra_operator_applications > 10); // Dense detection plus enlarged retry.
    let checked = cfg.guard_spectrum(&op, None, fixed.clone()).unwrap();
    assert_eq!(checked.eigenvalues, fixed.eigenvalues);
    assert_eq!(checked.eigenvectors, fixed.eigenvectors);
    assert_eq!(checked.residuals, fixed.residuals);
    assert_eq!(checked.spectral_completeness.unwrap().retries, 0);
}

/// Release-only matrix-free guard regression; exercises the complementary
/// operator rather than the small full-space fallback.
/// cargo test --release -p dft-solvers spectral_complement_guard -- --ignored --nocapture
#[test]
#[ignore = "independent spectral searches; run in release"]
fn spectral_complement_guard() {
    let n = 144;
    let op = Diagonal((0..n).map(|i| i as f64 * 0.05).collect());
    let mut start = ComplexMatrix::zeros(n, 2).unwrap();
    start[(10, 0)] = 1.0.into();
    start[(11, 1)] = 1.0.into();
    let cfg = config(2, n);
    let missed = cfg.solve(&op, Some(&start), None).unwrap();
    let checked = cfg.guard_spectrum(&op, None, missed).unwrap();
    assert!(checked.eigenvalues[0].abs() < 1e-9);
    assert!((checked.eigenvalues[1] - 0.05).abs() < 1e-9);
    assert_eq!(checked.spectral_completeness.as_ref().unwrap().retries, 1);
    let complete = cfg.guard_spectrum(&op, None, checked.clone()).unwrap();
    assert_eq!(complete.eigenvalues, checked.eigenvalues);
    assert_eq!(complete.eigenvectors, checked.eigenvectors);
    assert_eq!(complete.spectral_completeness.as_ref().unwrap().retries, 0);
    println!(
        "complement guard: eigenvalues={:?} Ha, extra H={}",
        complete.eigenvalues,
        complete
            .spectral_completeness
            .unwrap()
            .extra_operator_applications
    );
}

#[test]
fn ritz_rule_requires_two_values_and_required_empty_root() {
    let op = MatrixOperator(fourier(&[0.1, 0.8, 1.4, 2.0, 4.0, 6.0, 9.0, 12.0]));
    let mut solver = config(3, 8);
    solver.tolerance = 1e-9;
    let policy = BandTolerances {
        strict_bands: 3,
        buffer_tolerance: 5e-6,
    };
    solver.max_iterations = 1;
    assert!(matches!(
        solver.solve_ritz_changes(&op, None, None, policy),
        Err(SolverError::NotConverged { .. })
    ));
    solver.max_iterations = 100;
    let result = solver.solve_ritz_changes(&op, None, None, policy).unwrap();
    assert!(result.iterations >= 2);
    assert!(result.converged.iter().all(|x| *x));
    assert!(!result.residuals_certified);
    for (a, b) in result.eigenvalues.iter().zip([0.1, 0.8, 1.4]) {
        assert!((a - b).abs() < 1e-5);
    }
    // Full-space first solve cannot expand; it must still establish two Ritz values.
    let full = Davidson {
        bands: 8,
        max_subspace: 8,
        block_size: 8,
        ..solver
    };
    let r = full
        .solve_ritz_changes(
            &op,
            None,
            None,
            BandTolerances {
                strict_bands: 8,
                buffer_tolerance: 5e-6,
            },
        )
        .unwrap();
    assert_eq!(r.iterations, 2);
}

#[test]
fn hard_locking_certifies_all_roots_against_dense() {
    let n = 40;
    let mut diagonal: Vec<_> = (0..n).map(|i| i as f64 * 0.2).collect();
    diagonal[0] = -2.0;
    diagonal[2] = diagonal[1];
    let op = Diagonal(diagonal.clone());
    let mut guess = SplitMix64(31).block(n, 4).unwrap();
    guess.column_mut(0).fill(Complex64::default());
    guess[(0, 0)] = Complex64::new(1.0, 0.0);
    let mut solver = config(4, n);
    solver.max_subspace = 20;
    solver.tolerance = 1e-9;
    let solution = solver
        .solve_locked(
            &op,
            Some(&guess),
            None,
            BandTolerances {
                strict_bands: 4,
                buffer_tolerance: solver.tolerance,
            },
        )
        .unwrap();
    assess("locked", &op, &solution, &diagonal[..4], 1e-8);
    assert!(solution.residuals_certified);
    assert!(solution.locking.as_ref().unwrap().locked_roots > 0);
}

#[test]
fn hard_locking_unlocks_after_a_failed_true_residual_audit() {
    use std::sync::atomic::{AtomicBool, Ordering};
    // Fault injection: only the initial cached H product has a shifted first
    // root. A true application is always the unshifted operator. This tests
    // cache/certificate recovery, not a time-dependent physical Hamiltonian.
    struct InitialCacheFault {
        diagonal: Vec<f64>,
        first: AtomicBool,
    }
    impl LinearOperator for InitialCacheFault {
        fn dimension(&self) -> usize {
            self.diagonal.len()
        }
        fn apply(&self, x: &[Complex64], y: &mut [Complex64]) -> Result<(), SolverError> {
            for ((y, x), d) in y.iter_mut().zip(x).zip(&self.diagonal) {
                *y = x * d;
            }
            Ok(())
        }
        fn apply_block(&self, x: &ComplexMatrix, y: &mut ComplexMatrix) -> Result<(), SolverError> {
            let fault = self.first.swap(false, Ordering::Relaxed);
            for j in 0..x.cols() {
                self.apply(x.column(j), y.column_mut(j))?;
                if fault {
                    y[(0, j)] += 0.01 * x[(0, j)];
                }
            }
            Ok(())
        }
    }
    let n = 40;
    let diagonal: Vec<_> = (0..n).map(|j| j as f64 * 0.2 - 2.0).collect();
    let op = InitialCacheFault {
        diagonal: diagonal.clone(),
        first: AtomicBool::new(true),
    };
    let mut guess = SplitMix64(31).block(n, 4).unwrap();
    guess.column_mut(0).fill(Complex64::default());
    guess[(0, 0)] = Complex64::new(1.0, 0.0);
    let mut solver = config(4, n);
    solver.max_subspace = 20;
    solver.tolerance = 1e-9;
    let result = solver
        .solve_locked(
            &op,
            Some(&guess),
            None,
            BandTolerances {
                strict_bands: 4,
                buffer_tolerance: solver.tolerance,
            },
        )
        .unwrap();
    assess("unlock", &op, &result, &diagonal[..4], 1e-8);
    let report = result.locking.unwrap();
    assert!(report.unlocks > 0);
    assert!(
        report
            .failures
            .iter()
            .any(|&(_, root, residual)| root == 0 && residual > 0.009)
    );
}
