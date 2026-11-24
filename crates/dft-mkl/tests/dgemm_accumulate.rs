//! LP64, storage, and beta contracts. Run on dev08 with --features mkl.
#![cfg(feature = "mkl")]
use dft_mkl::{DgemmSite, Error, dgemm_add_at, dgemm_at};

#[test]
fn accumulate_padded_operands_and_transpose_match_integer_oracle() {
    for transpose in [false, true] {
        for (m, n, k) in [(1, 1, 1), (3, 5, 2), (257, 3, 7)] {
            let (ar, ac) = if transpose { (k, m) } else { (m, k) };
            let (lda, ldb) = (ar + 2, k + 3);
            let mut a = vec![f64::NAN; lda * ac];
            let mut b = vec![f64::NAN; ldb * n];
            for j in 0..ac {
                for i in 0..ar {
                    a[i + j * lda] = ((i + j * 3) % 7) as f64 - 3.;
                }
            }
            for j in 0..n {
                for i in 0..k {
                    b[i + j * ldb] = ((i * 2 + j) % 5) as f64 - 2.;
                }
            }
            let mut product = vec![0.; m * n];
            for j in 0..n {
                for i in 0..m {
                    for q in 0..k {
                        let av = if transpose {
                            a[q + i * lda]
                        } else {
                            a[i + q * lda]
                        };
                        product[i + j * m] += av * b[q + j * ldb];
                    }
                }
            }
            let mut c = vec![11.; m * n];
            for repeats in [1, 2] {
                dgemm_add_at(
                    DgemmSite::Other,
                    transpose,
                    m,
                    n,
                    k,
                    &a,
                    lda,
                    &b,
                    ldb,
                    &mut c,
                )
                .unwrap();
                for (&got, &ab) in c.iter().zip(&product) {
                    assert_eq!(got, 11. + repeats as f64 * ab);
                }
            }
            // Replacement still ignores the old C, including NaN, and beta=0
            // remains the API used by all existing callers.
            c.fill(f64::NAN);
            dgemm_at(
                DgemmSite::Other,
                transpose,
                m,
                n,
                k,
                &a,
                lda,
                &b,
                ldb,
                &mut c,
            )
            .unwrap();
            assert_eq!(c, product);
        }
    }
}

#[test]
fn zero_reduction_preserves_accumulator_bits_but_replace_clears() {
    let mut c = [-0., f64::from_bits(0x7ff8_0000_0000_1234), 7., -8.];
    let before = c.map(f64::to_bits);
    dgemm_add_at(DgemmSite::Other, false, 2, 2, 0, &[], 2, &[], 0, &mut c).unwrap();
    assert_eq!(before, c.map(f64::to_bits));
    dgemm_at(DgemmSite::Other, false, 2, 2, 0, &[], 2, &[], 0, &mut c).unwrap();
    assert_eq!(c.map(f64::to_bits), [0; 4]);
    for (m, n) in [(0, 3), (3, 0), (0, 0)] {
        dgemm_add_at(DgemmSite::Other, false, m, n, 2, &[], m, &[], 2, &mut []).unwrap();
    }
}

#[test]
fn invalid_accumulate_dimensions_fail_without_touching_c() {
    // m,n,k,lda,ldb,A length,B length,C length. Include LP64 and
    // multiplication overflow without constructing huge slices.
    let huge = i32::MAX as usize + 1;
    for (m, n, k, lda, ldb, alen, blen, clen) in [
        (2, 3, 4, 1, 4, 8, 12, 6),
        (2, 3, 4, 2, 3, 8, 12, 6),
        (2, 3, 4, 2, 4, 7, 12, 6),
        (2, 3, 4, 2, 4, 8, 11, 6),
        (2, 3, 4, 2, 4, 8, 12, 5),
        (huge, 0, 0, huge, 0, 0, 0, 0),
        (0, huge, 0, 0, 0, 0, 0, 0),
        (0, 0, huge, 0, huge, 0, 0, 0),
        (0, 0, 0, huge, 0, 0, 0, 0),
        (0, 0, 0, 0, huge, 0, 0, 0),
        (usize::MAX, usize::MAX, 0, 0, 0, 0, 0, 0),
    ] {
        let mut c = vec![17.; clen];
        assert_eq!(
            dgemm_add_at(
                DgemmSite::Other,
                false,
                m,
                n,
                k,
                &vec![1.; alen],
                lda,
                &vec![1.; blen],
                ldb,
                &mut c
            ),
            Err(Error::Dimension)
        );
        assert!(c.iter().all(|&v| v == 17.));
    }
    // Transposed A stores k rows, not m rows.
    let mut c = [17.; 6];
    assert_eq!(
        dgemm_add_at(
            DgemmSite::Other,
            true,
            2,
            3,
            4,
            &[1.; 8],
            2,
            &[1.; 12],
            4,
            &mut c
        ),
        Err(Error::Dimension)
    );
    assert_eq!(c, [17.; 6]);
}
