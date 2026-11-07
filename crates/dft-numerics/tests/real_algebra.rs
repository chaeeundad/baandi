//! R8 old-kernel oracle and panel alternatives. No SCF tolerance is changed.
use dft_numerics::{
    pool::WorkerPool,
    real::{RealGemmWorkspace, RealMatrix, RealMatrixMut, RealMatrixRef, RealOp},
};
use faer::{Accum, Par, linalg::matmul::matmul};

fn fixture(r: usize, c: usize) -> RealMatrix {
    RealMatrix::new(
        r,
        c,
        (0..r * c)
            .map(|i| ((i * 13 + 7) as f64).sin() / 32.)
            .collect(),
    )
    .unwrap()
}
// Previous implementation, retained verbatim as the rounding oracle.
fn previous(mut out: RealMatrixMut<'_>, a: RealMatrixRef<'_>, b: RealMatrixRef<'_>) {
    const MICRO: usize = 16;
    for j in (0..out.ncols()).step_by(MICRO) {
        let n = MICRO.min(out.ncols() - j);
        for i in (0..out.nrows()).step_by(MICRO) {
            let m = MICRO.min(out.nrows() - i);
            for k in (0..a.ncols()).step_by(MICRO) {
                let len = MICRO.min(a.ncols() - k);
                matmul(
                    out.as_mut().submatrix_mut(i, j, m, n),
                    if k == 0 { Accum::Replace } else { Accum::Add },
                    a.submatrix(i, k, m, len),
                    b.submatrix(k, j, len, n),
                    1.,
                    Par::Seq,
                );
            }
        }
    }
}
fn columns(mut out: RealMatrixMut<'_>, a: RealMatrixRef<'_>, b: RealMatrixRef<'_>) {
    for j in 0..b.ncols() {
        matmul(
            out.as_mut().subcols_mut(j, 1),
            Accum::Replace,
            a,
            b.subcols(j, 1),
            1.,
            Par::Seq,
        );
    }
}
fn slices(mut out: RealMatrixMut<'_>, a: RealMatrixRef<'_>, b: RealMatrixRef<'_>) {
    if a.row_stride() == 1 {
        let a = a.try_as_col_major().unwrap();
        for j in 0..b.ncols() {
            let dst = out
                .as_mut()
                .col_mut(j)
                .try_as_col_major_mut()
                .unwrap()
                .as_slice_mut();
            dst.fill(0.);
            for k in 0..a.ncols() {
                let weight = b[(k, j)];
                for (dst, &src) in dst.iter_mut().zip(a.col(k).as_slice()) {
                    *dst += src * weight;
                }
            }
        }
    } else {
        let at = a.transpose().try_as_col_major().unwrap();
        let b = b.try_as_col_major().unwrap();
        for j in 0..b.ncols() {
            let rhs = b.col(j).as_slice();
            let dst = out
                .as_mut()
                .col_mut(j)
                .try_as_col_major_mut()
                .unwrap()
                .as_slice_mut();
            for (i, dst) in dst.iter_mut().enumerate() {
                let lhs = at.col(i).as_slice();
                let mut sums = [0.; 4];
                let (lhs_blocks, lhs_tail) = lhs.as_chunks::<4>();
                let (rhs_blocks, rhs_tail) = rhs.as_chunks::<4>();
                for (x, y) in lhs_blocks.iter().zip(rhs_blocks) {
                    for lane in 0..4 {
                        sums[lane] += x[lane] * y[lane];
                    }
                }
                let mut sum = (sums[0] + sums[1]) + (sums[2] + sums[3]);
                for (&x, &y) in lhs_tail.iter().zip(rhs_tail) {
                    sum += x * y;
                }
                *dst = sum;
            }
        }
    }
}
fn shaped<const M: usize, const N: usize, const K: usize>(
    mut out: RealMatrixMut<'_>,
    a: RealMatrixRef<'_>,
    b: RealMatrixRef<'_>,
) {
    for j in (0..out.ncols()).step_by(N) {
        let n = N.min(out.ncols() - j);
        for i in (0..out.nrows()).step_by(M) {
            let m = M.min(out.nrows() - i);
            for k in (0..a.ncols()).step_by(K) {
                let len = K.min(a.ncols() - k);
                matmul(
                    out.as_mut().submatrix_mut(i, j, m, n),
                    if k == 0 { Accum::Replace } else { Accum::Add },
                    a.submatrix(i, k, m, len),
                    b.submatrix(k, j, len, n),
                    1.,
                    Par::Seq,
                );
            }
        }
    }
}
fn large(out: RealMatrixMut<'_>, a: RealMatrixRef<'_>, b: RealMatrixRef<'_>) {
    matmul(out, Accum::Replace, a, b, 1., Par::Seq);
}
fn error(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt()
        / b.iter()
            .map(|b| b * b)
            .sum::<f64>()
            .sqrt()
            .max(f64::MIN_POSITIVE)
}
#[test]
fn previous_naive_and_threads() {
    for (rows, width, cols) in [
        (1, 1, 1),
        (255, 7, 3),
        (256, 17, 9),
        (257, 33, 17),
        (1025, 33, 37),
        (17, 257, 3),
        (513, 9, 5),
    ] {
        for op in [RealOp::Identity, RealOp::Transpose] {
            let a = fixture(rows, width);
            let b = fixture(if op == RealOp::Identity { width } else { rows }, cols);
            let m = if op == RealOp::Identity { rows } else { width };
            let mut old = RealMatrix::zeros(m, cols).unwrap();
            let mut naive = old.clone();
            for row in (0..rows).step_by(256) {
                let len = 256.min(rows - row);
                if op == RealOp::Identity {
                    previous(
                        old.as_mut().subrows_mut(row, len),
                        a.as_ref().subrows(row, len),
                        b.as_ref(),
                    );
                } else {
                    let mut partial = RealMatrix::zeros(m, cols).unwrap();
                    previous(
                        partial.as_mut(),
                        a.as_ref().subrows(row, len).transpose(),
                        b.as_ref().subrows(row, len),
                    );
                    for (dst, src) in old.as_mut_slice().iter_mut().zip(partial.as_slice()) {
                        *dst += src;
                    }
                }
            }
            for j in 0..cols {
                for i in 0..m {
                    for k in 0..b.rows() {
                        naive[(i, j)] += if op == RealOp::Identity {
                            a[(i, k)]
                        } else {
                            a[(k, i)]
                        } * b[(k, j)];
                    }
                }
            }
            let mut first = None;
            for threads in [1, 2, 4, 8] {
                let mut ws = RealGemmWorkspace::with_pool(
                    rows,
                    width.max(cols),
                    WorkerPool::new(threads).unwrap(),
                )
                .unwrap();
                let mut out = RealMatrix::zeros(m, cols).unwrap();
                ws.multiply_into(a.as_ref(), op, b.as_ref(), out.as_mut())
                    .unwrap();
                // Only association inside fixed tiles may change. Relative Frobenius
                // error avoids ill-defined entrywise relative error at cancellation zeros.
                let e_old = error(out.as_slice(), old.as_slice());
                let e_naive = error(out.as_slice(), naive.as_slice());
                assert!(e_old < 1e-13, "old {e_old:e}");
                assert!(e_naive < 1e-14, "naive {e_naive:e}");
                let bits: Vec<_> = out.as_slice().iter().map(|v| v.to_bits()).collect();
                if let Some(ref first) = first {
                    assert_eq!(&bits, first)
                } else {
                    first = Some(bits)
                }
                if threads == 1 {
                    println!("{op:?} {rows}/{width}/{cols}: old={e_old:e} naive={e_naive:e}");
                }
            }
        }
    }
}
#[test]
#[ignore = "release-only panel timing: cargo test --release -p dft-numerics --test real_algebra panel_alternatives -- --ignored --nocapture"]
fn panel_alternatives() {
    for (rows, width, cols) in [(256, 20, 20), (256, 70, 20), (256, 120, 20), (256, 216, 36)] {
        let a = fixture(rows, width);
        for transpose in [false, true] {
            let b = fixture(if transpose { rows } else { width }, cols);
            let a = if transpose {
                a.as_ref().transpose()
            } else {
                a.as_ref()
            };
            let mut out = RealMatrix::zeros(a.nrows(), cols).unwrap();
            for (name, kernel) in [
                (
                    "previous",
                    previous
                        as for<'a, 'b, 'c> fn(
                            RealMatrixMut<'a>,
                            RealMatrixRef<'b>,
                            RealMatrixRef<'c>,
                        ),
                ),
                ("columns", columns),
                ("slices", slices),
                ("256x4x4", shaped::<256, 4, 4>),
                ("64x4x16", shaped::<64, 4, 16>),
                ("32x8x16", shaped::<32, 8, 16>),
                ("4x4x256", shaped::<4, 4, 256>),
                ("8x4x128", shaped::<8, 4, 128>),
                ("large_allocating", large),
            ] {
                let mut samples = Vec::new();
                for _ in 0..7 {
                    let start = std::time::Instant::now();
                    for _ in 0..100 {
                        kernel(out.as_mut(), a, b.as_ref());
                        std::hint::black_box(out.as_slice());
                    }
                    samples.push(start.elapsed().as_secs_f64() / 100.);
                }
                samples.sort_by(f64::total_cmp);
                println!(
                    "panel {rows}/{width}/{cols} transpose={transpose} {name} median_s={:.9}",
                    samples[3]
                );
            }
        }
    }
}

#[test]
fn strided_views_and_zero_reduction() {
    // Transposed storage exercises non-unit row strides and padded output.
    let a = fixture(7, 257);
    let b = fixture(5, 7);
    let mut backing = RealMatrix::zeros(5, 257).unwrap();
    let mut ws = RealGemmWorkspace::with_pool(257, 7, WorkerPool::new(4).unwrap()).unwrap();
    ws.multiply_into(
        a.as_ref().transpose(),
        RealOp::Identity,
        b.as_ref().transpose(),
        backing.as_mut().transpose_mut(),
    )
    .unwrap();
    for j in 0..5 {
        for i in 0..257 {
            let expected: f64 = (0..7).map(|k| a[(k, i)] * b[(j, k)]).sum();
            assert!((backing[(j, i)] - expected).abs() < 1e-14);
        }
    }
    let a = RealMatrix::zeros(257, 0).unwrap();
    let b = RealMatrix::zeros(0, 5).unwrap();
    ws.multiply_into(
        a.as_ref(),
        RealOp::Identity,
        b.as_ref(),
        backing.as_mut().transpose_mut(),
    )
    .unwrap();
    assert!(backing.as_slice().iter().all(|&v| v == 0.));
}
