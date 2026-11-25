//! GEMMT contract tests: execute on dev08 with `cargo test -p dft-mkl --features mkl --test dgemmt`.
#![cfg(feature = "mkl")]
use dft_mkl::{Error, dgemmt_upper};

#[test]
fn upper_only_padded_offset_and_distinct_operands() {
    for (n, k) in [(1, 1), (3, 7), (17, 259)] {
        let (lda, ldb, ldc) = (k + 2, k + 3, n + 4);
        let mut a = vec![f64::NAN; lda * n];
        let mut b = vec![f64::NAN; ldb * n];
        for j in 0..n {
            for i in 0..k {
                a[i + j * lda] = ((i + 3 * j) % 7) as f64 - 3.;
                b[i + j * ldb] = ((2 * i + j) % 5) as f64 - 2.;
            }
        }
        let sentinel = f64::from_bits(0x7ff8_0000_0000_1234);
        let mut out = vec![sentinel; 2 + ldc * n];
        dgemmt_upper(n, k, &a, lda, &b, ldb, &mut out[2..], ldc).unwrap();
        assert!(out[..2].iter().all(|v| v.to_bits() == sentinel.to_bits()));
        for j in 0..n {
            for i in 0..ldc {
                if i <= j {
                    let expected: f64 = (0..k).map(|q| a[q + i * lda] * b[q + j * ldb]).sum();
                    assert_eq!(out[2 + i + j * ldc], expected);
                } else {
                    assert_eq!(out[2 + i + j * ldc].to_bits(), sentinel.to_bits());
                }
            }
        }
    }
}

#[test]
fn empty_products_and_bad_storage_fail_before_ffi() {
    let mut c = [-0., 7., 8., -0., 11., 12.];
    dgemmt_upper(2, 0, &[], 0, &[], 0, &mut c, 3).unwrap();
    assert_eq!(c, [0., 7., 8., 0., 0., 12.]);
    dgemmt_upper(0, 3, &[], 0, &[], 0, &mut [], 0).unwrap();
    for (n, k, lda, ldb, ldc, alen, blen, clen) in [
        (2, 3, 2, 3, 2, 6, 6, 4),
        (2, 3, 3, 2, 2, 6, 6, 4),
        (2, 3, 3, 3, 1, 6, 6, 4),
        (2, 3, 3, 3, 2, 5, 6, 4),
        (2, 3, 3, 3, 2, 6, 5, 4),
        (2, 3, 3, 3, 2, 6, 6, 3),
        (i32::MAX as usize + 1, 0, 0, 0, 0, 0, 0, 0),
        (0, usize::MAX, 0, 0, 0, 0, 0, 0),
    ] {
        let mut out = vec![13.; clen];
        assert_eq!(
            dgemmt_upper(
                n,
                k,
                &vec![1.; alen],
                lda,
                &vec![1.; blen],
                ldb,
                &mut out,
                ldc
            ),
            Err(Error::Dimension)
        );
        assert!(out.iter().all(|v| *v == 13.));
    }
}
