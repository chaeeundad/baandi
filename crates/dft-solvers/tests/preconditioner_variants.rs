use dft_numerics::{linalg::ComplexMatrix, pool::WorkerPool, real::RealMatrix};
use dft_solvers::{
    BandTolerances, Davidson, DavidsonPreconditioner as Kind, DiagonalPreconditioner,
    LinearOperator, Preconditioner, RealDavidson, RealLinearOperator, SolverError,
};
use num_complex::Complex64 as C;
use std::sync::Arc;

#[test]
fn tpa_uses_selected_band_expectation_and_positive_shift_floor() {
    let d = vec![0., 0.1, 0.2, 1e200];
    let p = DiagonalPreconditioner::new(d.clone(), 0.1)
        .unwrap()
        .with_variant(Kind::Tpa);
    let r = ComplexMatrix::new(4, 1, vec![C::new(1., 0.); 4]).unwrap();
    let mut v = ComplexMatrix::zeros(4, 2).unwrap();
    v[(2, 0)] = C::new(2., 0.); // Normalization must not affect the expectation.
    v[(0, 1)] = C::new(1., 0.); // Pure G=0 uses the 0.1 Ha floor.
    let mut y = ComplexMatrix::zeros(4, 1).unwrap();
    p.apply_with_vectors(&[0.], &r, &v, &[0], &mut y).unwrap();
    assert_eq!(y[(0, 0)].re, -1.);
    assert!((y[(2, 0)].re + 65. / 81.).abs() < 1e-15);
    assert!(y[(3, 0)].re < 0. && y[(3, 0)].re.is_finite());
    p.apply_with_vectors(&[0.], &r, &v, &[1], &mut y).unwrap();
    assert!((y[(1, 0)].re + 65. / 81.).abs() < 1e-15);
    assert!(p.apply(&[0.], &r, &mut y).is_err());
    assert!(p.apply_with_vectors(&[0.], &r, &v, &[2], &mut y).is_err());
    let p = DiagonalPreconditioner::new(d, 0.1)
        .unwrap()
        .with_variant(Kind::RitzShift);
    p.apply(&[0.2], &r, &mut y).unwrap();
    for i in 0..3 {
        assert_eq!(y[(i, 0)].re, -10.);
    }
    assert!(y[(3, 0)].re < 0.);
}

// A tridiagonal chain has independently known eigenvalues d+2t cos(k*pi/(n+1)).
struct Chain {
    n: usize,
    pool: Arc<WorkerPool>,
}
impl LinearOperator for Chain {
    fn dimension(&self) -> usize {
        self.n
    }
    fn worker_pool(&self) -> Option<Arc<WorkerPool>> {
        Some(self.pool.clone())
    }
    fn apply(&self, x: &[C], y: &mut [C]) -> Result<(), SolverError> {
        for i in 0..self.n {
            y[i] = 2. * x[i];
            if i > 0 {
                y[i] += 0.2 * x[i - 1];
            }
            if i + 1 < self.n {
                y[i] += 0.2 * x[i + 1];
            }
        }
        Ok(())
    }
}
impl RealLinearOperator for Chain {
    fn dimension(&self) -> usize {
        self.n
    }
    fn worker_pool(&self) -> Option<Arc<WorkerPool>> {
        Some(self.pool.clone())
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), SolverError> {
        for i in 0..self.n {
            y[i] = 2. * x[i];
            if i > 0 {
                y[i] += 0.2 * x[i - 1];
            }
            if i + 1 < self.n {
                y[i] += 0.2 * x[i + 1];
            }
        }
        Ok(())
    }
}
#[test]
fn analytic_spectrum_and_thread_identity_both_representations() {
    let n = 24;
    let bands = 3;
    let d: Vec<_> = (0..n).map(|i| i as f64 * 0.2).collect();
    let mut initial = RealMatrix::zeros(n, bands).unwrap();
    for j in 0..bands {
        initial[(j, j)] = 1.;
    }
    let c_initial = ComplexMatrix::new(
        n,
        bands,
        initial.as_slice().iter().map(|&v| C::new(v, 0.)).collect(),
    )
    .unwrap();
    for kind in [Kind::Default, Kind::Tpa, Kind::RitzShift] {
        let mut previous = None;
        for threads in [1, 2, 4, 8, 12] {
            let op = Chain {
                n,
                pool: WorkerPool::new(threads).unwrap(),
            };
            let real = RealDavidson {
                bands,
                tolerance: 1e-10,
                max_iterations: 400,
                max_subspace: n,
                seed: 42,
            };
            let complex = Davidson {
                bands,
                tolerance: 1e-10,
                max_iterations: 400,
                max_subspace: n,
                block_size: bands,
                seed: 42,
            };
            let policy = BandTolerances {
                strict_bands: bands,
                buffer_tolerance: 1e-10,
            };
            let r = real
                .solve_reusing_with_preconditioner(
                    &op,
                    Some(&initial),
                    &d,
                    policy,
                    false,
                    true,
                    6,
                    &mut None,
                    kind,
                )
                .unwrap();
            let p = DiagonalPreconditioner::new(d.clone(), 0.1)
                .unwrap()
                .with_variant(kind);
            let c = complex.solve(&op, Some(&c_initial), Some(&p)).unwrap();
            assert!(r.residuals_certified && c.residuals_certified);
            for j in 0..bands {
                let expected =
                    2. - 0.4 * (((j + 1) as f64 * std::f64::consts::PI) / (n + 1) as f64).cos();
                assert!((r.eigenvalues[j] - expected).abs() < 1e-10);
                assert!((c.eigenvalues[j] - expected).abs() < 1e-10);
            }
            let bits = (
                r.eigenvalues
                    .iter()
                    .map(|v| v.to_bits())
                    .collect::<Vec<_>>(),
                r.eigenvectors
                    .as_slice()
                    .iter()
                    .map(|v| v.to_bits())
                    .collect::<Vec<_>>(),
                c.eigenvalues
                    .iter()
                    .map(|v| v.to_bits())
                    .collect::<Vec<_>>(),
                c.eigenvectors
                    .as_slice()
                    .iter()
                    .map(|v| (v.re.to_bits(), v.im.to_bits()))
                    .collect::<Vec<_>>(),
                r.iterations,
                r.operator_applications,
                c.iterations,
                c.operator_applications,
            );
            if let Some(ref before) = previous {
                assert_eq!(&bits, before, "{} T={threads}", kind.name());
            }
            previous = Some(bits);
        }
    }
}
