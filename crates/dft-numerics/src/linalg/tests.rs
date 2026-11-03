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

fn legacy_norm(x: &[Complex64]) -> f64 {
    x.iter().fold(0.0_f64, |s, z| s.hypot(z.re).hypot(z.im))
}
fn legacy_orthonormalize(
    block: &ComplexMatrix,
    against: Option<&ComplexMatrix>,
) -> Result<Orthonormalized, LinalgError> {
    block.validate_finite()?;
    if let Some(q) = against {
        q.validate_finite()?;
        if q.rows != block.rows {
            return Err(LinalgError::Dimension);
        }
    }
    let mut data = Vec::new();
    let mut dropped = Vec::new();
    let mut count = 0;
    for j in 0..block.cols {
        let mut v = block.column(j).to_vec();
        let original = legacy_norm(&v);
        if !original.is_finite() {
            return Err(LinalgError::NonFinite);
        }
        // Normalize first so rank decisions also work for very small input columns.
        if original == 0.0 {
            dropped.push(j);
            continue;
        }
        for z in &mut v {
            *z /= original;
        }
        for _ in 0..2 {
            if let Some(q) = against {
                for k in 0..q.cols {
                    legacy_subtract_projection(&mut v, q.column(k));
                }
            }
            for k in 0..count {
                legacy_subtract_projection(&mut v, &data[k * block.rows..(k + 1) * block.rows]);
            }
        }
        let remaining = legacy_norm(&v);
        if !remaining.is_finite() {
            return Err(LinalgError::NonFinite);
        }
        if remaining <= DEPENDENCE_TOLERANCE {
            dropped.push(j);
        } else {
            data.extend(v.iter().map(|z| z / remaining));
            count += 1;
        }
    }
    Ok(Orthonormalized {
        vectors: ComplexMatrix::new(block.rows, count, data)?,
        dropped,
    })
}
fn legacy_subtract_projection(v: &mut [Complex64], q: &[Complex64]) {
    let dot: Complex64 = q.iter().zip(v.iter()).map(|(q, v)| q.conj() * v).sum();
    for (v, q) in v.iter_mut().zip(q) {
        *v -= dot * q;
    }
}

#[test]
fn norm_extremes_and_orthogonalization_oracle() {
    for scale in [0.0, 1e-300, 1e-100, 1.0, 1e300] {
        let x = [
            Complex64::new(3.0 * scale, 4.0 * scale),
            Complex64::new(1e-300, 0.0),
        ];
        let expected = legacy_norm(&x);
        assert!((norm(&x) / expected - 1.0).abs() <= 5e-16);
    }
    assert_eq!(norm(&[]), 0.0);
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let a = matrix(2, 1, |_, _| Complex64::new(bad, 0.0));
        assert_eq!(
            orthonormalize(&a, None).unwrap_err(),
            LinalgError::NonFinite
        );
    }
    // Orthogonal Fourier modes give singular values sqrt(8+eps^2), eps.
    for epsilon in [
        1.0,
        (8.0_f64 / (1e12 - 1.0)).sqrt(),
        (8.0_f64 / (1e20 - 1.0)).sqrt(),
        0.0,
    ] {
        for scale in [1e-300, 1.0, 1e300] {
            let a = matrix(32, 8, |i, j| {
                let phase = std::f64::consts::TAU * i as f64 / 32.0;
                scale
                    * (Complex64::from_polar(1.0, phase)
                        + epsilon * Complex64::from_polar(1.0, phase * (j + 2) as f64))
            });
            let old = legacy_orthonormalize(&a, None).unwrap();
            let new = orthonormalize(&a, None).unwrap();
            assert_eq!(old.dropped, new.dropped);
            let gram = gemm(&new.vectors, Adjoint, &new.vectors, Identity).unwrap();
            let mut defect: f64 = 0.0;
            for j in 0..gram.cols() {
                for i in 0..gram.rows() {
                    defect = defect.max((gram[(i, j)] - if i == j { 1.0 } else { 0.0 }).norm());
                }
            }
            println!("epsilon={epsilon:e} scale={scale:e} orthogonality={defect:e}");
            assert!(defect <= 2e-15);
        }
    }
}

#[test]
fn random_blocks_against_oracle() {
    let mut state = 17_u64;
    for _ in 0..12 {
        let mut a = ComplexMatrix::zeros(64, 12).unwrap();
        for z in &mut a.data {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            z.re = (state >> 11) as f64 / (1_u64 << 53) as f64 - 0.5;
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            z.im = (state >> 11) as f64 / (1_u64 << 53) as f64 - 0.5;
        }
        let q = orthonormalize(&a, None).unwrap().vectors;
        let old = legacy_orthonormalize(&a, Some(&q)).unwrap();
        let new = orthonormalize(&a, Some(&q)).unwrap();
        assert_eq!(old.dropped, new.dropped);
        assert_eq!(new.vectors.cols(), 0);
        let gram = gemm(&q, Adjoint, &q, Identity).unwrap();
        for j in 0..12 {
            for i in 0..12 {
                assert!((gram[(i, j)] - if i == j { 1.0 } else { 0.0 }).norm() <= 2e-15);
            }
        }
    }
}

/// Run: cargo test --release -p dft-numerics kernel_timing -- --ignored --nocapture
#[test]
#[ignore]
fn kernel_timing() {
    use std::{hint::black_box, time::Instant};
    let x: Vec<_> = (0..2048)
        .map(|i| Complex64::new((i as f64).sin(), (i as f64).cos()))
        .collect();
    for (name, f) in [
        ("old norm", legacy_norm as fn(&[Complex64]) -> f64),
        ("new norm", norm),
    ] {
        let t = Instant::now();
        for _ in 0..20000 {
            black_box(f(black_box(&x)));
        }
        println!("{name}: {} s", t.elapsed().as_secs_f64());
    }
    for (name, f) in [
        (
            "old projection",
            legacy_subtract_projection as fn(&mut [Complex64], &[Complex64]),
        ),
        ("new projection", subtract_projection),
    ] {
        let mut v = x.clone();
        let q: Vec<_> = x.iter().map(|z| z / norm(&x)).collect();
        let t = Instant::now();
        for _ in 0..20000 {
            f(black_box(&mut v), black_box(&q));
        }
        println!("{name}: {} s", t.elapsed().as_secs_f64());
    }
}

#[test]
fn reusable_evd_gemm_and_rank_filter() {
    let a = matrix(7, 7, |i, j| {
        Complex64::new(
            if i == j {
                i as f64 + 1.0
            } else {
                0.03 * (i + j) as f64
            },
            0.02 * (i as f64 - j as f64),
        )
    });
    let (expected, _) = hermitian_eigh(&a).unwrap();
    let mut workspace = HermitianWorkspace::new(7).unwrap();
    let mut values = vec![0.0; 7];
    let mut vectors = ComplexMatrix::zeros(7, 7).unwrap();
    let vector_pointer = vectors.data.as_ptr();
    let diagonal_pointer = workspace.diagonal.as_ref().column_vector().as_ptr();
    let mut product = ComplexMatrix::zeros(7, 7).unwrap();
    for _ in 0..3 {
        workspace.solve(&a, &mut values, &mut vectors).unwrap();
        assert_eq!(values, expected);
        assert_eq!(vectors.data.as_ptr(), vector_pointer);
        assert_eq!(
            workspace.diagonal.as_ref().column_vector().as_ptr(),
            diagonal_pointer
        );
        for aa in [MatrixOp::Identity, MatrixOp::Adjoint] {
            for bb in [MatrixOp::Identity, MatrixOp::Adjoint] {
                gemm_into(&a, aa, &vectors, bb, &mut product).unwrap();
                let expected = gemm(&a, aa, &vectors, bb).unwrap();
                for (x, y) in product.as_slice().iter().zip(expected.as_slice()) {
                    assert!((*x - *y).norm() < 1e-12);
                }
            }
        }
    }
    let mut block = matrix(7, 4, |i, j| if i == j { 1.0.into() } else { 0.0.into() });
    for i in 0..7 {
        block[(i, 3)] = block[(i, 1)];
    }
    let oracle = orthonormalize(&block, None).unwrap();
    let mut output = ComplexMatrix::zeros(7, 4).unwrap();
    let mut scratch = vec![Complex64::default(); 7];
    assert_eq!(
        orthonormalize_into(&block, &mut output, &mut scratch).unwrap(),
        3
    );
    for j in 0..3 {
        assert_eq!(output.column(j), oracle.vectors.column(j));
    }
    assert!(output.column(3).iter().all(|z| *z == Complex64::default()));
    let wrong = ComplexMatrix::zeros(6, 6).unwrap();
    assert!(workspace.solve(&wrong, &mut values, &mut vectors).is_err());
    let mut bad = a.clone();
    bad[(0, 1)] += 1.0;
    assert!(matches!(
        workspace.solve(&bad, &mut values, &mut vectors),
        Err(LinalgError::NotHermitian { .. })
    ));
    bad[(0, 0)] = f64::NAN.into();
    assert!(matches!(
        workspace.solve(&bad, &mut values, &mut vectors),
        Err(LinalgError::NonFinite)
    ));
}

#[test]
fn reusable_adapters_rectangular_and_empty_shapes() {
    for aa in [MatrixOp::Identity, MatrixOp::Adjoint] {
        for bb in [MatrixOp::Identity, MatrixOp::Adjoint] {
            let (ar, ac) = match aa {
                MatrixOp::Identity => (2, 3),
                MatrixOp::Adjoint => (3, 2),
            };
            let (br, bc) = match bb {
                MatrixOp::Identity => (3, 4),
                MatrixOp::Adjoint => (4, 3),
            };
            let a = matrix(ar, ac, |i, j| {
                Complex64::new((i + 2 * j) as f64, 0.3 * (i as f64 - j as f64))
            });
            let b = matrix(br, bc, |i, j| {
                Complex64::new((2 * i + j) as f64, 0.2 * (i as f64 + j as f64))
            });
            let mut c = ComplexMatrix::zeros(2, 4).unwrap();
            gemm_into(&a, aa, &b, bb, &mut c).unwrap();
            let oracle = gemm(&a, aa, &b, bb).unwrap();
            for (v, w) in c.as_slice().iter().zip(oracle.as_slice()) {
                assert!((*v - *w).norm() < 1e-12);
            }
            let mut bad = ComplexMatrix::zeros(4, 2).unwrap();
            assert!(matches!(
                gemm_into(&a, aa, &b, bb, &mut bad),
                Err(LinalgError::Dimension)
            ));
        }
    }
    let mut workspace = HermitianWorkspace::new(0).unwrap();
    let empty = ComplexMatrix::zeros(0, 0).unwrap();
    let mut output = ComplexMatrix::zeros(0, 0).unwrap();
    workspace.solve(&empty, &mut [], &mut output).unwrap();
    gemm_into(
        &empty,
        MatrixOp::Identity,
        &empty,
        MatrixOp::Identity,
        &mut output,
    )
    .unwrap();
    assert_eq!(
        orthonormalize_into(&empty, &mut output, &mut []).unwrap(),
        0
    );
}
