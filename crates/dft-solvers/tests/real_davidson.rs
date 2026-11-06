use dft_numerics::{linalg::ComplexMatrix, pool::WorkerPool, real::RealMatrix};
use dft_solvers::{
    BandTolerances, Davidson, DenseEigensolver, DiagonalPreconditioner, LinearOperator,
    RealDavidson, RealLinearOperator, SolverError, real_davidson::free_electron_initial_guess,
};
use num_complex::Complex64;
use std::{sync::Arc, time::Instant};

struct Operator {
    a: RealMatrix,
    pool: Arc<WorkerPool>,
}
impl RealLinearOperator for Operator {
    fn dimension(&self) -> usize {
        self.a.rows()
    }
    fn worker_pool(&self) -> Option<Arc<WorkerPool>> {
        Some(self.pool.clone())
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), SolverError> {
        for (i, y) in y.iter_mut().enumerate() {
            *y = x.iter().enumerate().map(|(j, x)| self.a[(i, j)] * x).sum();
        }
        Ok(())
    }
}
impl LinearOperator for Operator {
    fn dimension(&self) -> usize {
        self.a.rows()
    }
    fn worker_pool(&self) -> Option<Arc<WorkerPool>> {
        Some(self.pool.clone())
    }
    fn apply(&self, x: &[Complex64], y: &mut [Complex64]) -> Result<(), SolverError> {
        for (i, y) in y.iter_mut().enumerate() {
            *y = x.iter().enumerate().map(|(j, x)| self.a[(i, j)] * x).sum();
        }
        Ok(())
    }
}
fn complex(x: &RealMatrix) -> ComplexMatrix {
    ComplexMatrix::new(
        x.rows(),
        x.cols(),
        x.as_slice()
            .iter()
            .map(|&x| Complex64::new(x, 0.))
            .collect(),
    )
    .unwrap()
}
fn solver(bands: usize) -> RealDavidson {
    RealDavidson {
        bands,
        tolerance: 1e-11,
        max_iterations: 400,
        max_subspace: 6 * bands,
        seed: 42,
    }
}
fn complex_solver(s: &RealDavidson) -> Davidson {
    Davidson {
        bands: s.bands,
        tolerance: s.tolerance,
        max_iterations: s.max_iterations,
        max_subspace: s.max_subspace.max(6 * s.bands),
        block_size: s.bands,
        seed: s.seed,
    }
}
fn fixture(kind: usize, n: usize) -> (RealMatrix, Vec<f64>, Option<Vec<f64>>) {
    let kinetic: Vec<_> = (0..n).map(|i| 0.2 * (i / 2) as f64).collect();
    let mut a = RealMatrix::zeros(n, n).unwrap();
    if kind == 1 {
        for j in 0..n {
            for i in 0..=j {
                let v = if i == j {
                    kinetic[i] + 0.3
                } else {
                    0.025 * ((i * 31 + j * 17 + 1) as f64).sin() / ((j - i) as f64).sqrt()
                };
                a[(i, j)] = v;
                a[(j, i)] = v;
            }
        }
        return (a, kinetic, None);
    }
    // Deterministic pseudo-random Householder reflections; Q is orthogonal and
    // Q diag(lambda) Q^T has a known spectrum independently of any EVD backend.
    let mut q = RealMatrix::zeros(n, n).unwrap();
    for i in 0..n {
        q[(i, i)] = 1.;
    }
    let mut state = 753_u64;
    for _ in 0..n {
        let v: Vec<_> = (0..n)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                ((state >> 11) as f64) / (1_u64 << 53) as f64 - 0.5
            })
            .collect();
        let scale = 2. / v.iter().map(|x| x * x).sum::<f64>();
        for j in 0..n {
            let d = (0..n).map(|i| v[i] * q[(i, j)]).sum::<f64>() * scale;
            for i in 0..n {
                q[(i, j)] -= v[i] * d;
            }
        }
    }
    let values: Vec<_> = (0..n)
        .map(|i| {
            if kind == 2 {
                (i / 3) as f64 * 0.2
            } else {
                i as f64 * 0.13 - 0.7
            }
        })
        .collect();
    for j in 0..n {
        for i in 0..=j {
            let v = (0..n).map(|k| q[(i, k)] * values[k] * q[(j, k)]).sum();
            a[(i, j)] = v;
            a[(j, i)] = v;
        }
    }
    (a, kinetic, Some(values))
}
fn check(
    op: &Operator,
    result: &dft_solvers::RealEigenSolution,
    oracle: &[f64],
    tolerance: f64,
) -> (f64, f64, f64) {
    let n = op.a.rows();
    let x = &result.eigenvectors;
    let mut error = 0_f64;
    let mut residual = 0_f64;
    let mut orth = 0_f64;
    for j in 0..x.cols() {
        error = error.max((result.eigenvalues[j] - oracle[j]).abs());
        let mut r2 = 0.;
        for i in 0..n {
            let r = (0..n).map(|k| op.a[(i, k)] * x[(k, j)]).sum::<f64>()
                - result.eigenvalues[j] * x[(i, j)];
            r2 += r * r;
        }
        residual = residual.max(r2.sqrt());
        for k in 0..x.cols() {
            orth = orth.max(
                ((0..n).map(|i| x[(i, j)] * x[(i, k)]).sum::<f64>() - if j == k { 1. } else { 0. })
                    .abs(),
            );
        }
    }
    assert!(error < 1e-12, "eigenvalue error {error:e} Ha");
    assert!(residual <= tolerance, "residual {residual:e} Ha");
    // Same dimensionless restart-orthogonality threshold as complex Davidson.
    assert!(orth < 1e-12, "orthogonality {orth:e}");
    (error, residual, orth)
}
#[test]
fn fixed_operators_dense_and_complex() {
    for kind in 0..3 {
        let (a, diagonal, known) = fixture(kind, 49);
        let op = Operator {
            a,
            pool: WorkerPool::new(1).unwrap(),
        };
        let s = solver(if kind == 2 { 4 } else { 3 });
        let initial = free_electron_initial_guess(&diagonal, s.bands).unwrap();
        let dense = DenseEigensolver.solve(&op, s.bands).unwrap();
        if let Some(known) = known {
            for (a, b) in dense.eigenvalues.iter().zip(known) {
                assert!((a - b).abs() < 1e-12);
            }
        }
        let start = Instant::now();
        let real = s.solve(&op, Some(&initial), &diagonal).unwrap();
        let real_time = start.elapsed().as_secs_f64();
        let start = Instant::now();
        let other = complex_solver(&s)
            .solve(
                &op,
                Some(&complex(&initial)),
                Some(&DiagonalPreconditioner::new(diagonal, 0.1).unwrap()),
            )
            .unwrap();
        let complex_time = start.elapsed().as_secs_f64();
        let (error, residual, orth) = check(&op, &real, &dense.eigenvalues, s.tolerance);
        for (a, b) in real.eigenvalues.iter().zip(&other.eigenvalues) {
            assert!((a - b).abs() < 1e-12);
        }
        assert!(real.residuals_certified && other.residuals_certified);
        assert!(real.certification_attempts > 0);
        if kind == 0 {
            assert!(real.contracted_restarts > 0);
        }
        println!(
            "kind={kind} real/complex steps={}/{} H={}/{} restarts={} contractions={} eigen_error={error:e} Ha residual={residual:e} Ha orth={orth:e} Mac_seconds={real_time:.6}/{complex_time:.6}",
            real.iterations,
            other.iterations,
            real.operator_applications,
            other.operator_applications,
            real.restarts,
            real.contracted_restarts
        );
    }
}
#[test]
fn buffer_and_ritz_semantics() {
    let (a, diagonal, _) = fixture(1, 49);
    let op = Operator {
        a,
        pool: WorkerPool::new(1).unwrap(),
    };
    let mut s = solver(4);
    s.tolerance = 1e-9;
    let policy = BandTolerances {
        strict_bands: 2,
        buffer_tolerance: 1e-6,
    };
    let initial = free_electron_initial_guess(&diagonal, s.bands).unwrap();
    let p = DiagonalPreconditioner::new(diagonal.clone(), 0.1).unwrap();
    let real = s
        .solve_with_band_tolerances(&op, Some(&initial), &diagonal, policy)
        .unwrap();
    let other = complex_solver(&s)
        .solve_with_band_tolerances(&op, Some(&complex(&initial)), Some(&p), policy)
        .unwrap();
    assert_eq!(real.buffer, other.buffer);
    assert_eq!(real.converged, other.converged);
    assert!(real.residuals_certified && other.residuals_certified);
    for (j, &r) in real.residuals.iter().enumerate() {
        assert!(
            r <= if j < 2 {
                s.tolerance
            } else {
                policy.buffer_tolerance
            }
        );
    }
    for (a, b) in real.eigenvalues.iter().zip(&other.eigenvalues) {
        assert!((a - b).abs() < 1e-12);
    }
    println!(
        "buffer real/complex steps={}/{} H={}/{}",
        real.iterations, other.iterations, real.operator_applications, other.operator_applications
    );
    let real = s
        .solve_ritz_changes(&op, Some(&initial), &diagonal, policy)
        .unwrap();
    let other = complex_solver(&s)
        .solve_ritz_changes(&op, Some(&complex(&initial)), Some(&p), policy)
        .unwrap();
    assert!(!real.residuals_certified && !other.residuals_certified);
    assert_eq!(real.converged, other.converged);
    for (a, b) in real.eigenvalues.iter().zip(&other.eigenvalues) {
        assert!((a - b).abs() < 1e-12);
    }
    println!(
        "ritz real/complex steps={}/{} H={}/{}",
        real.iterations, other.iterations, real.operator_applications, other.operator_applications
    );
}
#[test]
fn cold_start_errors_and_full_space() {
    let d = [0., 2., 2., 1., 1.];
    let x = free_electron_initial_guess(&d, 4).unwrap();
    for (j, i) in [0, 3, 4, 1].into_iter().enumerate() {
        assert_eq!(x[(i, j)], 1.);
    }
    assert!(free_electron_initial_guess(&[f64::NAN], 1).is_err());
    let mut a = RealMatrix::zeros(5, 5).unwrap();
    for i in 0..5 {
        a[(i, i)] = d[i];
    }
    let op = Operator {
        a,
        pool: WorkerPool::new(1).unwrap(),
    };
    let s = solver(5);
    let strict = s.solve(&op, None, &d).unwrap();
    assert_eq!(strict.iterations, 1);
    assert_eq!(strict.operator_applications, 10);
    let standard = s
        .solve_ritz_changes(
            &op,
            None,
            &d,
            BandTolerances {
                strict_bands: 5,
                buffer_tolerance: s.tolerance,
            },
        )
        .unwrap();
    assert_eq!(standard.iterations, 2);
    assert!(!standard.residuals_certified);
    assert!(s.solve(&op, None, &d[..4]).is_err());
    let mut invalid = s.clone();
    invalid.tolerance = f64::NAN;
    assert!(invalid.solve(&op, None, &d).is_err());
}
#[test]
fn small_thread_identity() {
    let (a, d, _) = fixture(1, 65);
    let s = solver(5);
    let mut first = None;
    for threads in [1, 2, 4, 8] {
        let op = Operator {
            a: a.clone(),
            pool: WorkerPool::new(threads).unwrap(),
        };
        let result = s.solve(&op, None, &d).unwrap();
        let bits: Vec<_> = result
            .eigenvectors
            .as_slice()
            .iter()
            .chain(&result.eigenvalues)
            .chain(&result.residuals)
            .chain(&result.history)
            .map(|x| x.to_bits())
            .collect();
        if let Some((ref old, steps, apps)) = first {
            assert_eq!(bits, *old);
            assert_eq!(result.iterations, steps);
            assert_eq!(result.operator_applications, apps);
        } else {
            first = Some((bits, result.iterations, result.operator_applications));
        }
    }
}

// Sparse kinetic+potential operator exercises real GEMM dispatch without a
// quadratic dense H cost. Every supplied block must have normalized columns.
struct Large {
    pool: Arc<WorkerPool>,
    n: usize,
}
impl RealLinearOperator for Large {
    fn dimension(&self) -> usize {
        self.n
    }
    fn worker_pool(&self) -> Option<Arc<WorkerPool>> {
        Some(self.pool.clone())
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), SolverError> {
        let norm = x.iter().map(|x| x * x).sum::<f64>();
        assert!((norm - 1.).abs() < 1e-12);
        for i in 0..self.n {
            y[i] = (i as f64 * 0.02 + 0.5) * x[i];
            if i > 0 {
                y[i] += 0.015 * x[i - 1];
            }
            if i + 1 < self.n {
                y[i] += 0.015 * x[i + 1];
            }
            y[i] += 0.001 * x[(i + 513) % self.n];
            y[i] += 0.001 * x[(i + self.n - 513) % self.n];
        }
        Ok(())
    }
}
#[test]
#[ignore = "release-mode pool-dispatch test: cargo test --release -p dft-solvers --test real_davidson parallel_dispatch -- --ignored --nocapture"]
fn parallel_dispatch() {
    let n = 1025;
    let d: Vec<_> = (0..n).map(|i| i as f64 * 0.02).collect();
    let s = solver(17);
    let mut first = None;
    for threads in [1, 2, 4, 8] {
        let pool = WorkerPool::new(threads).unwrap();
        let op = Large {
            pool: pool.clone(),
            n,
        };
        let created = pool.workers_created();
        let t = Instant::now();
        let result = s.solve(&op, None, &d).unwrap();
        let elapsed = t.elapsed().as_secs_f64();
        assert_eq!(pool.workers_created(), created);
        assert!(n * result.max_basis_dimension * s.bands >= 1_048_576);
        let bits: Vec<_> = result
            .eigenvectors
            .as_slice()
            .iter()
            .chain(&result.eigenvalues)
            .chain(&result.residuals)
            .chain(&result.history)
            .map(|x| x.to_bits())
            .collect();
        if let Some((ref old, steps, apps)) = first {
            assert_eq!(bits, *old);
            assert_eq!(result.iterations, steps);
            assert_eq!(result.operator_applications, apps);
        } else {
            first = Some((bits, result.iterations, result.operator_applications));
        }
        println!(
            "parallel T={threads} steps={} H={} Mac_seconds={elapsed:.6}",
            result.iterations, result.operator_applications
        );
    }
}

#[test]
fn final_limit_requires_strict_certificate_but_allows_failed_buffers() {
    let mut a = RealMatrix::zeros(7, 7).unwrap();
    for i in 0..7 {
        a[(i, i)] = i as f64;
    }
    a[(2, 3)] = 0.1;
    a[(3, 2)] = 0.1;
    let d: Vec<_> = (0..7).map(|i| i as f64).collect();
    let op = Operator {
        a,
        pool: WorkerPool::new(1).unwrap(),
    };
    let mut s = solver(3);
    s.max_iterations = 1;
    let policy = BandTolerances {
        strict_bands: 2,
        buffer_tolerance: 1e-6,
    };
    let initial = free_electron_initial_guess(&d, 3).unwrap();
    let r = s
        .solve_with_band_tolerances(&op, Some(&initial), &d, policy)
        .unwrap();
    let c = complex_solver(&s)
        .solve_with_band_tolerances(
            &op,
            Some(&complex(&initial)),
            Some(&DiagonalPreconditioner::new(d.clone(), 0.1).unwrap()),
            policy,
        )
        .unwrap();
    assert_eq!(r.converged, vec![true, true, false]);
    assert_eq!(r.converged, c.converged);
    assert_eq!(r.eigenvalues, c.eigenvalues);
    assert!(r.residuals_certified && c.residuals_certified);
    assert_eq!(r.operator_applications, c.operator_applications);
    // Making the unconverged root mandatory must fail, not silently certify.
    assert!(s.solve(&op, Some(&initial), &d).is_err());
}

// Deliberately different first-block rounding: the recurrence appears exact,
// while all subsequent applications use the real symmetric matrix. This
// adversarial oracle verifies that certification rejects a stale product.
struct CertificateProbe {
    op: Operator,
    calls: std::cell::Cell<usize>,
}
impl RealLinearOperator for CertificateProbe {
    fn dimension(&self) -> usize {
        self.op.a.rows()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), SolverError> {
        RealLinearOperator::apply(&self.op, x, y)
    }
    fn apply_block(&self, x: &RealMatrix, y: &mut RealMatrix) -> Result<(), SolverError> {
        let calls = self.calls.get();
        self.calls.set(calls + 1);
        if calls == 0 {
            for j in 0..x.cols() {
                for i in 0..x.rows() {
                    y[(i, j)] = self.op.a[(i, i)] * x[(i, j)];
                }
            }
            Ok(())
        } else {
            RealLinearOperator::apply_block(&self.op, x, y)
        }
    }
}
impl LinearOperator for CertificateProbe {
    fn dimension(&self) -> usize {
        self.op.a.rows()
    }
    fn apply(&self, x: &[Complex64], y: &mut [Complex64]) -> Result<(), SolverError> {
        LinearOperator::apply(&self.op, x, y)
    }
    fn apply_block(&self, x: &ComplexMatrix, y: &mut ComplexMatrix) -> Result<(), SolverError> {
        let calls = self.calls.get();
        self.calls.set(calls + 1);
        if calls == 0 {
            for j in 0..x.cols() {
                for i in 0..x.rows() {
                    y[(i, j)] = self.op.a[(i, i)] * x[(i, j)];
                }
            }
            Ok(())
        } else {
            LinearOperator::apply_block(&self.op, x, y)
        }
    }
}
#[test]
fn failed_certificate_restarts_with_true_residuals() {
    let mut a = RealMatrix::zeros(13, 13).unwrap();
    for i in 0..13 {
        a[(i, i)] = i as f64;
    }
    a[(0, 3)] = 0.03;
    a[(3, 0)] = 0.03;
    a[(1, 4)] = 0.02;
    a[(4, 1)] = 0.02;
    let op = CertificateProbe {
        op: Operator {
            a,
            pool: WorkerPool::new(1).unwrap(),
        },
        calls: std::cell::Cell::new(0),
    };
    let d: Vec<_> = (0..13).map(|i| i as f64).collect();
    let s = solver(2);
    let x = free_electron_initial_guess(&d, 2).unwrap();
    let r = s.solve(&op, Some(&x), &d).unwrap();
    assert!(r.certification_attempts >= 2);
    assert!(r.restarts > r.contracted_restarts);
    assert!(r.residuals_certified);
    op.calls.set(0);
    let c = complex_solver(&s)
        .solve(
            &op,
            Some(&complex(&x)),
            Some(&DiagonalPreconditioner::new(d, 0.1).unwrap()),
        )
        .unwrap();
    for (a, b) in r.eigenvalues.iter().zip(&c.eigenvalues) {
        assert!((a - b).abs() < 1e-12);
    }
    assert!(c.residuals_certified);
    println!(
        "failed certificate real/complex steps={}/{} H={}/{} certificates={}",
        r.iterations,
        c.iterations,
        r.operator_applications,
        c.operator_applications,
        r.certification_attempts
    );
}

#[test]
fn recurrence_return_skips_final_application_but_limit_certifies() {
    let mut a = RealMatrix::zeros(7, 7).unwrap();
    for i in 0..7 {
        a[(i, i)] = i as f64;
    }
    let op = CertificateProbe {
        op: Operator {
            a,
            pool: WorkerPool::new(1).unwrap(),
        },
        calls: std::cell::Cell::new(0),
    };
    let d: Vec<_> = (0..7).map(|i| i as f64).collect();
    let policy = BandTolerances {
        strict_bands: 2,
        buffer_tolerance: 1e-6,
    };
    for (certify, limit, blocks, vectors) in
        [(false, 400, 1, 3), (true, 400, 2, 6), (false, 1, 2, 6)]
    {
        op.calls.set(0);
        let mut s = solver(3);
        s.max_iterations = limit;
        let r = s
            .solve_with_band_tolerances_and_certification(&op, None, &d, policy, certify)
            .unwrap();
        assert_eq!(op.calls.get(), blocks);
        assert_eq!(r.operator_applications, vectors);
        assert_eq!(r.certification_attempts, blocks - 1);
        assert_eq!(r.residuals_certified, blocks == 2);
        assert_eq!(r.iterations, 1);
        assert_eq!(r.eigenvalues, vec![0., 1., 2.]);
        assert!(r.converged.iter().all(|&x| x));
    }
}
