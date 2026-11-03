#![cfg(feature = "mkl")]
use dft_mkl::{Error, Matrix, Plan3d, gemm};
use dft_numerics::{Fft3d, RustFftPlan3d};
use ndarray::Array3;
use num_complex::Complex64 as C;

fn data(n: usize) -> Vec<C> {
    let mut seed = 19u64;
    let mut next = || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5
    };
    (0..n).map(|_| C::new(next(), next())).collect()
}
#[test]
fn fft_rust_agreement_roundtrip_parseval_and_moved_plan() {
    for shape in [
        [1, 1, 1],
        [2, 3, 5],
        [4, 5, 6],
        [7, 4, 3],
        [24, 24, 24],
        [25, 25, 25],
    ] {
        let input = data(shape.iter().product());
        let a = Array3::from_shape_vec(shape, input.clone()).unwrap();
        let mut b = a.clone();
        let mut rust = RustFftPlan3d::new(shape).unwrap();
        let mut mkl = Plan3d::new(shape).unwrap();
        let mut out = input.clone();
        for inverse in [false, true] {
            mkl.transform(&input, &mut out, inverse).unwrap();
            if inverse {
                rust.to_real(&a, &mut b).unwrap();
            } else {
                rust.to_reciprocal(&a, &mut b).unwrap();
            }
            let mut max: f64 = 0.0;
            for (actual, expected) in out.iter().zip(&b) {
                assert!(actual.re.is_finite() && actual.im.is_finite());
                assert!(expected.re.is_finite() && expected.im.is_finite());
                let error = (*actual - *expected).norm();
                let rel = if error == 0.0 {
                    0.0
                } else {
                    error / expected.norm()
                };
                max = max.max(rel);
            }
            println!("shape={shape:?} inverse={inverse} max_element_relative={max:e}");
            assert!(max <= 1e-13);
        }
        mkl.transform(&input, &mut out, false).unwrap();
        let p: f64 = input.iter().map(|z| z.norm_sqr()).sum();
        let q: f64 = out.iter().map(|z| z.norm_sqr()).sum();
        let parseval = (q * input.len() as f64 / p - 1.0).abs();
        let mut back = input.clone();
        mkl.transform(&out, &mut back, true).unwrap();
        assert!(back.iter().all(|z| z.re.is_finite() && z.im.is_finite()));
        let roundtrip = input
            .iter()
            .zip(back)
            .map(|(x, y)| (*x - y).norm())
            .fold(0.0, f64::max);
        println!("shape={shape:?} roundtrip_max_abs={roundtrip:e} parseval_relative={parseval:e}");
        assert!(roundtrip <= 1e-13 && parseval <= 1e-13);
        // Test Send by transferring an already committed descriptor to a worker.
        let actual = std::thread::spawn(move || {
            let mut result = input.clone();
            mkl.transform(&input, &mut result, false).unwrap();
            result
        })
        .join()
        .unwrap();
        assert!(
            actual
                .iter()
                .zip(out)
                .all(|(a, b)| a.re.to_bits() == b.re.to_bits() && a.im.to_bits() == b.im.to_bits())
        );
    }
}
#[test]
fn dimensions_and_lengths_fail_before_ffi() {
    assert!(matches!(Plan3d::new([0, 2, 3]), Err(Error::Dimension)));
    assert!(matches!(
        Plan3d::new([usize::MAX, 2, 3]),
        Err(Error::Dimension)
    ));
    let mut plan = Plan3d::new([2, 3, 4]).unwrap();
    assert_eq!(plan.transform(&[], &mut [], true), Err(Error::Dimension));
    assert!(matches!(Matrix::new(&[], 1, 1), Err(Error::Dimension)));
    assert!(matches!(
        Matrix::new(&[], i32::MAX as usize + 1, 0),
        Err(Error::Dimension)
    ));
    assert_eq!(
        gemm(
            Matrix::new(&[], 0, 2).unwrap(),
            false,
            Matrix::new(&[], 3, 0).unwrap(),
            false,
            &mut []
        ),
        Err(Error::Dimension)
    );
    let mut c = vec![C::new(9.0, 3.0); 6];
    gemm(
        Matrix::new(&[], 2, 0).unwrap(),
        false,
        Matrix::new(&[], 0, 3).unwrap(),
        false,
        &mut c,
    )
    .unwrap();
    assert_eq!(c, vec![C::default(); 6]);
}
#[test]
fn all_gemm_adjoint_variants_match_scalar_oracle() {
    for (m, k, n) in [(1, 1, 1), (3, 5, 2), (17, 9, 11)] {
        for aa in [false, true] {
            for ab in [false, true] {
                let (ar, ac) = if aa { (k, m) } else { (m, k) };
                let (br, bc) = if ab { (n, k) } else { (k, n) };
                let a = data(ar * ac);
                let b = data(br * bc);
                let mut c = vec![C::default(); m * n];
                gemm(
                    Matrix::new(&a, ar, ac).unwrap(),
                    aa,
                    Matrix::new(&b, br, bc).unwrap(),
                    ab,
                    &mut c,
                )
                .unwrap();
                for j in 0..n {
                    for i in 0..m {
                        let mut expected = C::default();
                        for q in 0..k {
                            let x = if aa {
                                a[q + i * ar].conj()
                            } else {
                                a[i + q * ar]
                            };
                            let y = if ab {
                                b[j + q * br].conj()
                            } else {
                                b[q + j * br]
                            };
                            expected += x * y;
                        }
                        assert!(
                            (c[i + j * m] - expected).norm() <= 1e-13 * expected.norm().max(1.0)
                        );
                    }
                }
            }
        }
    }
}
