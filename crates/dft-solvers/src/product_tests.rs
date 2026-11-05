use super::*;
use std::cell::Cell;

struct Counted<'a> {
    inner: &'a dyn LinearOperator,
    columns: Cell<usize>,
    blocks: Cell<usize>,
}
impl LinearOperator for Counted<'_> {
    fn dimension(&self) -> usize {
        self.inner.dimension()
    }
    fn apply(&self, x: &[Complex64], y: &mut [Complex64]) -> Result<(), SolverError> {
        self.columns.set(self.columns.get() + 1);
        self.inner.apply(x, y)
    }
    fn apply_block(&self, x: &ComplexMatrix, y: &mut ComplexMatrix) -> Result<(), SolverError> {
        self.blocks.set(self.blocks.get() + 1);
        self.columns.set(self.columns.get() + x.cols());
        self.inner.apply_block(x, y)
    }
}
// Keep operator arithmetic independent of block width when comparing the
// frozen legacy solver with cached products and the new projection arithmetic.
struct ColumnMatrix(ComplexMatrix);
impl LinearOperator for ColumnMatrix {
    fn dimension(&self) -> usize {
        self.0.rows()
    }
    fn apply(&self, x: &[Complex64], y: &mut [Complex64]) -> Result<(), SolverError> {
        for (i, y) in y.iter_mut().enumerate() {
            *y = x.iter().enumerate().map(|(j, x)| self.0[(i, j)] * x).sum();
        }
        Ok(())
    }
}
fn compare(
    label: &str,
    cfg: Davidson,
    op: &dyn LinearOperator,
    initial: Option<&ComplexMatrix>,
    p: Option<&dyn Preconditioner>,
) {
    let counter = Counted {
        inner: op,
        columns: Cell::new(0),
        blocks: Cell::new(0),
    };
    let old = cfg.solve_legacy(&counter, initial, p).unwrap();
    assert_eq!(counter.columns.get(), old.operator_applications);
    counter.columns.set(0);
    let new = cfg
        .solve_with_recurrence(&counter, initial, p, false)
        .unwrap();
    assert_eq!(counter.columns.get(), new.operator_applications);
    // P10I changes projection/GEMM rounding. The unchanged legacy solver is
    // still an independent numerical oracle, but no longer a bitwise oracle.
    let dense = DenseEigensolver.solve(op, cfg.bands).unwrap();
    assess(label, op, &old, &dense.eigenvalues, cfg.tolerance);
    assess(label, op, &new, &dense.eigenvalues, cfg.tolerance);
    for (a, b) in old.eigenvalues.iter().zip(&new.eigenvalues) {
        assert!((a - b).abs() < 1e-12, "legacy eigenvalue agreement");
    }
    assert!(new.operator_applications < old.operator_applications);
    println!(
        "{label}: applications {} -> {}, ratio {:.6}, iterations {}, numerical agreement",
        old.operator_applications,
        new.operator_applications,
        old.operator_applications as f64 / new.operator_applications as f64,
        new.iterations
    );
    let recurrent = cfg.solve_with_recurrence(op, initial, p, true).unwrap();
    let dense = DenseEigensolver.solve(op, cfg.bands).unwrap();
    assess(label, op, &recurrent, &dense.eigenvalues, cfg.tolerance);
    let delta = |a: &[f64], b: &[f64]| {
        a.iter()
            .zip(b)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max)
    };
    println!(
        "{label}: stage2 applications {}, iterations {} -> {}, eigenvalue delta {:e}, residual delta {:e} (dimensionless)",
        recurrent.operator_applications,
        new.iterations,
        recurrent.iterations,
        delta(&new.eigenvalues, &recurrent.eigenvalues),
        delta(&new.residuals, &recurrent.residuals)
    );
}

/// Run in a source-only measurement build, isolated from other trace producers:
/// P14A_WORK_ACCOUNTING=1 cargo test --release -p dft-solvers work_accounting_counts -- --ignored
#[test]
#[ignore = "requires P14A_WORK_ACCOUNTING=1 measurement build"]
#[allow(clippy::assertions_on_constants)] // Fail explicitly if invoked in the wrong build.
fn work_accounting_counts() {
    assert!(work_accounting::ENABLED);
    let op = ColumnMatrix(random_matrix(18, 83));
    let counted = Counted {
        inner: &op,
        columns: Cell::new(0),
        blocks: Cell::new(0),
    };
    let solution = config(3, 9).solve(&counted, None, None).unwrap();
    let records = work_accounting::take();
    assert_eq!(records.len(), 1);
    let record = &records[0];
    assert_eq!(record.block_calls, counted.blocks.get());
    assert_eq!(record.vector_applications, counted.columns.get());
    assert_eq!(record.vector_applications, solution.operator_applications);
    assert_eq!(record.iterations, solution.iterations);
    assert!(record.restarts > 0);
    assert!(record.returned);
    assert!(record.first_below.iter().all(|&i| i > 0));
    assert!(record.last_above.iter().all(|&i| i < record.iterations));
}

#[test]
fn cached_products_oracle_and_counts() {
    let mut a = random_matrix(36, 149);
    for j in 0..36 {
        a[(j, j)] += 20.0 / (j + 1) as f64;
    }
    compare(
        "decaying diagonal",
        config(4, 36),
        &ColumnMatrix(a),
        None,
        None,
    );
    let op = ColumnMatrix(random_matrix(18, 83));
    let mut cfg = config(3, 9);
    cfg.tolerance = 1e-10;
    let mut initial = ComplexMatrix::zeros(18, 3).unwrap();
    initial[(0, 0)] = 1.0.into();
    initial[(0, 1)] = 1.0.into();
    compare("restart deficient guess", cfg, &op, Some(&initial), None);
    let mut a = ComplexMatrix::zeros(48, 48).unwrap();
    for i in 0..48 {
        a[(i, i)] = (0.5 * (i * i) as f64).into();
        if i > 0 {
            a[(i, i - 1)] = Complex64::new(0.02, 0.01);
            a[(i - 1, i)] = a[(i, i - 1)].conj();
        }
    }
    let p =
        DiagonalPreconditioner::new((0..48).map(|i| 0.5 * (i * i) as f64).collect(), 0.1).unwrap();
    compare(
        "sparse kinetic preconditioned",
        config(4, 24),
        &ColumnMatrix(a),
        None,
        Some(&p),
    );
    compare(
        "degenerate",
        config(3, 8),
        &ColumnMatrix(fourier(&[-1., -1., -1., 2., 3., 4., 5., 6.])),
        None,
        None,
    );
}

#[test]
fn failed_true_residual_certificate_cannot_return() {
    // First product is a stale identity product. Subsequent applications use
    // diag(2,3), so the recurrence initially reports zero but truth disagrees.
    struct Drift {
        calls: Cell<usize>,
    }
    impl LinearOperator for Drift {
        fn dimension(&self) -> usize {
            2
        }
        fn apply(&self, x: &[Complex64], y: &mut [Complex64]) -> Result<(), SolverError> {
            let first = self.calls.get() == 0;
            self.calls.set(self.calls.get() + 1);
            for i in 0..2 {
                y[i] = x[i] * if first { 1.0 } else { 2.0 + i as f64 };
            }
            Ok(())
        }
    }
    let mut initial = ComplexMatrix::zeros(2, 1).unwrap();
    initial[(0, 0)] = 1.0.into();
    let mut cfg = config(1, 2);
    cfg.max_iterations = 1;
    let op = Drift {
        calls: Cell::new(0),
    };
    match cfg.solve_with_recurrence(&op, Some(&initial), None, true) {
        Err(SolverError::NotConverged {
            iterations,
            max_residual,
            ..
        }) => {
            assert_eq!(iterations, 1);
            assert_eq!(max_residual, 1.0);
        }
        other => panic!("premature certificate: {other:?}"),
    }
    cfg.max_iterations = 20;
    let op = Drift {
        calls: Cell::new(0),
    };
    let s = cfg
        .solve_with_recurrence(&op, Some(&initial), None, true)
        .unwrap();
    assert!(s.iterations > 1);
    assert_eq!(s.history[0], 1.0);
    assert!(s.residuals[0] <= cfg.tolerance);
    assert!((s.eigenvalues[0] - 2.0).abs() < cfg.tolerance);
    println!(
        "drift certificate: history {:?}, applications {}",
        s.history, s.operator_applications
    );
}

#[test]
fn replay_trace_matches_independent_counts_and_dense_spectrum() {
    let n = 12;
    let mut matrix = ComplexMatrix::zeros(n, n).unwrap();
    for i in 0..n {
        matrix[(i, i)] = Complex64::new(i as f64 * 0.4, 0.0);
        if i > 0 {
            matrix[(i, i - 1)] = Complex64::new(0.07, 0.02);
            matrix[(i - 1, i)] = matrix[(i, i - 1)].conj();
        }
    }
    let op = ColumnMatrix(matrix);
    let counter = Counted {
        inner: &op,
        columns: Cell::new(0),
        blocks: Cell::new(0),
    };
    let solver = Davidson {
        bands: 3,
        tolerance: 1e-10,
        max_iterations: 80,
        max_subspace: 10,
        block_size: 3,
        seed: 42,
    };
    let mode = SolveMode {
        ritz: false,
        recurrence: true,
        contract_restart: true,
        certify: true,
    };
    let policy = BandTolerances {
        strict_bands: 3,
        buffer_tolerance: 1e-10,
    };
    let mut trace = replay::Trace::new(3, 80);
    let measured = solver
        .solve_policy_traced::<true>(&counter, None, None, mode, policy, &mut trace)
        .unwrap();
    assert_eq!(trace.blocks, counter.blocks.get());
    assert_eq!(trace.rows.last().unwrap().1, counter.columns.get());
    assert_eq!(trace.rows.len(), measured.iterations);
    assert_eq!(trace.values.len(), 3 * measured.iterations);
    let ordinary = solver
        .solve_policy_traced::<false>(&op, None, None, mode, policy, &mut replay::Trace::default())
        .unwrap();
    assert_eq!(ordinary.eigenvalues, measured.eigenvalues);
    assert_eq!(ordinary.residuals, measured.residuals);
    assert_eq!(ordinary.history, measured.history);
    assert_eq!(
        ordinary.operator_applications,
        measured.operator_applications
    );
    let dense = DenseEigensolver.solve(&op, 3).unwrap();
    for (a, b) in dense.eigenvalues.iter().zip(&measured.eigenvalues) {
        assert!((a - b).abs() < 1e-12);
    }
    assert!(
        trace
            .rows
            .windows(2)
            .all(|p| p[1].1 >= p[0].1 && p[1].2 >= p[0].2)
    );
}
