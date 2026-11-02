use super::*;
use MatrixOp::{Adjoint, Identity};
fn matrix(n: usize, m: usize, f: impl Fn(usize, usize) -> Complex64) -> ComplexMatrix {
    let mut a = ComplexMatrix::zeros(n, m).unwrap();
    for j in 0..m {
        for i in 0..n {
            a[(i, j)] = f(i, j);
        }
    }
    a
}
fn quality(a: &ComplexMatrix, w: &[f64], v: &ComplexMatrix) -> (f64, f64) {
    let av = gemm(a, Identity, v, Identity).unwrap();
    let gram = gemm(v, Adjoint, v, Identity).unwrap();
    let r = (0..v.cols())
        .map(|j| {
            norm(
                &(0..v.rows())
                    .map(|i| av[(i, j)] - w[j] * v[(i, j)])
                    .collect::<Vec<_>>(),
            )
        })
        .fold(0.0, f64::max);
    let o = (0..v.cols())
        .flat_map(|j| (0..v.cols()).map(move |i| (i, j)))
        .map(|(i, j)| (gram[(i, j)] - Complex64::new(if i == j { 1.0 } else { 0.0 }, 0.0)).norm())
        .fold(0.0, f64::max);
    (r, o)
}
#[test]
fn analytic_and_fourier() {
    let a = ComplexMatrix::new(
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
    let d = [-2.0, 0.5, 3.0, 7.0];
    let u = matrix(4, 4, |i, j| {
        Complex64::from_polar(0.5, std::f64::consts::TAU * (i * j) as f64 / 4.0)
    });
    let ud = matrix(4, 4, |i, j| u[(i, j)] * d[j]);
    let b = gemm(&ud, Identity, &u, Adjoint).unwrap();
    for (a, expected) in [(a, vec![1.0, 4.0]), (b, d.to_vec())] {
        let (w, v) = hermitian_eigh(&a).unwrap();
        let error = w
            .iter()
            .zip(expected)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        let (r, o) = quality(&a, &w, &v);
        println!(
            "numerics analytic N={} values={w:?} error={error:e} residual={r:e} orth={o:e} (dimensionless)",
            a.rows()
        );
        assert!(error < 1e-12 && r < 1e-12 && o < 1e-12);
    }
}
#[test]
fn seeded_sixty() {
    let mut seed = 123_u64;
    let mut next = || {
        seed = seed.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = seed;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        ((z ^ (z >> 31)) >> 11) as f64 / (1_u64 << 53) as f64 - 0.5
    };
    let mut a = ComplexMatrix::zeros(60, 60).unwrap();
    for j in 0..60 {
        for i in 0..=j {
            let z = Complex64::new(next(), if i == j { 0.0 } else { next() });
            a[(i, j)] = z;
            a[(j, i)] = z.conj();
        }
    }
    let (w, v) = hermitian_eigh(&a).unwrap();
    let (r, o) = quality(&a, &w, &v);
    println!("numerics N=60 residual={r:e} orth={o:e} (dimensionless)");
    assert!(r < 1e-10 && o < 1e-10);
}
#[test]
fn gemm_rectangular_all_operations() {
    for oa in [Identity, Adjoint] {
        for ob in [Identity, Adjoint] {
            let f = |i, j| Complex64::new((i + 2 * j) as f64 / 7.0, (2 * i + j + 1) as f64 / 11.0);
            let a = match oa {
                Identity => matrix(3, 4, f),
                Adjoint => matrix(4, 3, f),
            };
            let b = match ob {
                Identity => matrix(4, 2, f),
                Adjoint => matrix(2, 4, f),
            };
            let c = gemm(&a, oa, &b, ob).unwrap();
            let mut error: f64 = 0.0;
            for j in 0..2 {
                for i in 0..3 {
                    let expected: Complex64 = (0..4)
                        .map(|k| {
                            let x = match oa {
                                Identity => a[(i, k)],
                                Adjoint => a[(k, i)].conj(),
                            };
                            let y = match ob {
                                Identity => b[(k, j)],
                                Adjoint => b[(j, k)].conj(),
                            };
                            x * y
                        })
                        .sum();
                    error = error.max((c[(i, j)] - expected).norm());
                }
            }
            println!("gemm {oa:?}/{ob:?} error={error:e} (dimensionless)");
            assert!(error < 1e-12);
        }
    }
}
#[test]
fn rank_filter_and_validation() {
    let q = matrix(4, 1, |i, _| if i == 0 { 1.0.into() } else { 0.0.into() });
    let a = matrix(4, 4, |i, j| {
        if i == 0 || i == 1 + j.min(1) {
            1.0.into()
        } else {
            0.0.into()
        }
    });
    let result = orthonormalize(&a, Some(&q)).unwrap();
    assert_eq!(result.dropped, vec![2, 3]);
    let gram = gemm(&result.vectors, Adjoint, &result.vectors, Identity).unwrap();
    let cross = gemm(&q, Adjoint, &result.vectors, Identity).unwrap();
    let mut error = cross
        .as_slice()
        .iter()
        .map(|z| z.norm())
        .fold(0.0, f64::max);
    for j in 0..2 {
        for i in 0..2 {
            error = error.max((gram[(i, j)] - if i == j { 1.0 } else { 0.0 }).norm());
        }
    }
    println!(
        "orthonormalize dropped={:?} error={error:e} (dimensionless)",
        result.dropped
    );
    assert!(error < 1e-12);
    let bad = matrix(2, 2, |i, j| Complex64::new((i + 2 * j) as f64, 0.0));
    assert!(matches!(
        hermitian_eigh(&bad),
        Err(LinalgError::NotHermitian { .. })
    ));
    assert!(matches!(hermitian_eigh(&q), Err(LinalgError::NotSquare)));
    assert!(ComplexMatrix::new(2, 2, vec![]).is_err());
    assert!(ComplexMatrix::zeros(usize::MAX, 2).is_err());
    assert!(gemm(&q, Identity, &q, Identity).is_err());
    let nan = matrix(1, 1, |_, _| Complex64::new(f64::NAN, 0.0));
    assert!(matches!(hermitian_eigh(&nan), Err(LinalgError::NonFinite)));
    let empty = ComplexMatrix::zeros(0, 0).unwrap();
    assert!(hermitian_eigh(&empty).unwrap().0.is_empty());
}
