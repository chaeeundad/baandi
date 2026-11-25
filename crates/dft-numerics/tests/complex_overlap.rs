//! Dimensionless generalized-Hermitian contracts, independent explicit products.
use dft_numerics::linalg::{ComplexMatrix as Matrix, ComplexOverlapFactor, LinalgError};
use num_complex::Complex64 as C;

#[test]
fn complex_generalized_residual_metric_and_rank_rule() {
    let mut max_residual = 0_f64;
    let mut max_metric = 0_f64;
    let mut factor = ComplexOverlapFactor::with_capacity(9).unwrap();
    for n in [1, 2, 5, 9, 2, 1] {
        for scale in [1e-3, 1., 1e3] {
            let mut b = Matrix::zeros(n, n).unwrap();
            for j in 0..n {
                for i in 0..=j {
                    b[(i, j)] = if i == j {
                        C::new(scale * (1. + i as f64 / n as f64), 0.)
                    } else {
                        scale * C::new(0.04 * (i + j + 1) as f64, 0.03 * (j - i) as f64)
                    };
                }
            }
            let mut s = Matrix::zeros(n, n).unwrap();
            let mut h = s.clone();
            for j in 0..n {
                for i in 0..n {
                    for k in 0..n {
                        s[(i, j)] += b[(k, i)].conj() * b[(k, j)];
                        h[(i, j)] += b[(k, i)].conj() * b[(k, j)] * (k as f64 - 2.);
                    }
                }
            }
            assert!(factor.factor(&s, f64::EPSILON / 1e-12).unwrap());
            let (e, x) = factor.solve_lowest(&h, n).unwrap();
            for j in 0..n {
                assert!((e[j] - (j as f64 - 2.)).abs() < 1e-12);
                for i in 0..n {
                    let r: C = (0..n)
                        .map(|k| (h[(i, k)] - e[j] * s[(i, k)]) * x[(k, j)])
                        .sum();
                    max_residual = max_residual.max(r.norm() / scale);
                    let m: C = (0..n)
                        .flat_map(|k| (0..n).map(move |l| (k, l)))
                        .map(|(k, l)| x[(k, i)].conj() * s[(k, l)] * x[(l, j)])
                        .sum();
                    max_metric = max_metric.max((m - C::new(f64::from(i == j), 0.)).norm());
                }
            }
        }
    }
    assert!(max_residual < 1e-12, "{max_residual:e}");
    assert!(max_metric < 1e-12, "{max_metric:e}");
    println!(
        "complex generalized relative residual={max_residual:e}, metric={max_metric:e} (dimensionless)"
    );
    let guard = f64::EPSILON / 1e-12;
    for (ratio, accepted) in [(guard * 0.5, false), (guard * 2., true)] {
        let off = C::new(0., (1. - ratio).sqrt());
        let s = Matrix::new(2, 2, vec![C::new(1., 0.), off.conj(), off, C::new(1., 0.)]).unwrap();
        assert_eq!(
            ComplexOverlapFactor::new(&s, guard).unwrap().is_some(),
            accepted
        );
    }
}
#[test]
fn complex_overlap_rejects_invalid_inputs() {
    assert!(matches!(
        ComplexOverlapFactor::new(&Matrix::zeros(2, 1).unwrap(), 1e-4),
        Err(LinalgError::NotSquare)
    ));
    let mut s = Matrix::zeros(2, 2).unwrap();
    assert!(ComplexOverlapFactor::new(&s, 1e-4).unwrap().is_none());
    s[(0, 0)] = C::new(1., 0.1);
    assert!(matches!(
        ComplexOverlapFactor::new(&s, 1e-4),
        Err(LinalgError::NotHermitian { .. })
    ));
    s[(0, 0)] = C::new(f64::NAN, 0.);
    assert!(matches!(
        ComplexOverlapFactor::new(&s, 1e-4),
        Err(LinalgError::NonFinite)
    ));
}
