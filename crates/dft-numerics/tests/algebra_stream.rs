//! Algebra-stream kernels of the real GEMM workspace (docs/performance/requests/
//! algebra-stream.md): the fused Ritz-and-residual pass is bit-identical to the
//! pair rotation plus the elementwise residual, the late single-side rotation is
//! bit-identical to the matching side of the pair, the narrow projection kernels
//! agree with a naive product to rounding, and every kernel gives the same bits
//! at every pool size. Dimensionless synthetic matrices.
use dft_numerics::{
    pool::WorkerPool,
    real::{RealGemmWorkspace, RealMatrix, RealOp, RotateKernel},
};

fn fixture(r: usize, c: usize, seed: u64) -> RealMatrix {
    let mut state = seed;
    RealMatrix::new(
        r,
        c,
        (0..r * c)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                ((state >> 11) as f64) / (1_u64 << 53) as f64 - 0.5
            })
            .collect(),
    )
    .unwrap()
}
fn bits(x: &[f64]) -> Vec<u64> {
    x.iter().map(|v| v.to_bits()).collect()
}

#[test]
fn fused_ritz_residual_is_bit_identical_to_pair_rotation() {
    for (m, k, n) in [
        (1, 3, 2),
        (255, 5, 3),
        (256, 7, 4),
        (257, 33, 17),
        (1025, 40, 9),
        (4100, 64, 17),
    ] {
        let v = fixture(m, k, 1 + m as u64);
        let hv = fixture(m, k, 2 + k as u64);
        let c = fixture(k, n, 3 + n as u64);
        let values: Vec<f64> = (0..n).map(|j| 0.37 * j as f64 - 1.1).collect();
        let mut reference: Option<(Vec<u64>, Vec<u64>, Vec<u64>)> = None;
        for threads in [1, 2, 4, 8] {
            let pool = WorkerPool::new(threads).unwrap();
            // A panel of H X must fit one tile of the workspace (short^2).
            let short = k
                .max(n)
                .max(((256.min(m) * n) as f64).sqrt().ceil() as usize);
            let mut ws = RealGemmWorkspace::with_pool(m, short, pool).unwrap();
            let (mut pv, mut ph) = (vec![0.; m * n], vec![0.; m * n]);
            ws.rotate_pair_slices(
                v.as_slice(),
                hv.as_slice(),
                k,
                c.as_slice(),
                n,
                &mut pv,
                &mut ph,
            )
            .unwrap();
            let expected_r: Vec<f64> = (0..m * n).map(|e| ph[e] - values[e / m] * pv[e]).collect();
            let (mut fv, mut fr) = (vec![0.; m * n], vec![0.; m * n]);
            assert!(
                ws.rotate_residual_slices(
                    v.as_slice(),
                    hv.as_slice(),
                    k,
                    c.as_slice(),
                    n,
                    &values,
                    &mut fv,
                    &mut fr
                )
                .unwrap()
            );
            assert_eq!(bits(&fv), bits(&pv), "X {m}x{k}x{n} T={threads}");
            assert_eq!(bits(&fr), bits(&expected_r), "R {m}x{k}x{n} T={threads}");
            let mut late = vec![0.; m * n];
            ws.rotate_slices(hv.as_slice(), k, c.as_slice(), n, &mut late)
                .unwrap();
            assert_eq!(bits(&late), bits(&ph), "late HX {m}x{k}x{n} T={threads}");
            let now = (bits(&fv), bits(&fr), bits(&late));
            match &reference {
                Some(r) => assert_eq!(&now, r, "thread identity {m}x{k}x{n} T={threads}"),
                None => reference = Some(now),
            }
            // The research kernel declines the fused form; its late rotation is
            // the Blocked side of the pair.
            ws.set_rotate_kernel(RotateKernel::Blocked);
            if !dft_numerics::MKL_BACKEND {
                assert!(
                    !ws.rotate_residual_slices(
                        v.as_slice(),
                        hv.as_slice(),
                        k,
                        c.as_slice(),
                        n,
                        &values,
                        &mut fv,
                        &mut fr
                    )
                    .unwrap()
                );
                ws.rotate_pair_slices(
                    v.as_slice(),
                    hv.as_slice(),
                    k,
                    c.as_slice(),
                    n,
                    &mut pv,
                    &mut ph,
                )
                .unwrap();
                ws.rotate_slices(hv.as_slice(), k, c.as_slice(), n, &mut late)
                    .unwrap();
                assert_eq!(bits(&late), bits(&ph));
            }
        }
    }
}

#[test]
fn fused_rejects_bad_shapes_and_non_finite() {
    let pool = WorkerPool::new(1).unwrap();
    // Tiles of 4 x 4 cannot hold a 16 x 2 panel: the fused form declines.
    let mut small = RealGemmWorkspace::with_pool(16, 4, pool.clone()).unwrap();
    let v = fixture(16, 4, 5);
    let hv = fixture(16, 4, 6);
    let c = fixture(4, 2, 7);
    let (mut x, mut r) = (vec![0.; 32], vec![0.; 32]);
    assert!(
        !small
            .rotate_residual_slices(
                v.as_slice(),
                hv.as_slice(),
                4,
                c.as_slice(),
                2,
                &[0., 1.],
                &mut x,
                &mut r
            )
            .unwrap()
    );
    let mut ws = RealGemmWorkspace::with_pool(16, 8, pool).unwrap();
    let v = fixture(16, 4, 5);
    let mut hv = fixture(16, 4, 6);
    let c = fixture(4, 2, 7);
    let (mut x, mut r) = (vec![0.; 32], vec![0.; 32]);
    assert!(
        ws.rotate_residual_slices(
            v.as_slice(),
            hv.as_slice(),
            4,
            c.as_slice(),
            2,
            &[0.],
            &mut x,
            &mut r
        )
        .is_err()
    );
    hv[(3, 1)] = f64::NAN;
    assert!(
        ws.rotate_residual_slices(
            v.as_slice(),
            hv.as_slice(),
            4,
            c.as_slice(),
            2,
            &[0., 1.],
            &mut x,
            &mut r
        )
        .is_err()
    );
}

#[test]
fn narrow_kernels_match_naive_products_and_thread_count() {
    // Shapes above and below the parallel dispatch, partial row blocks and
    // pieces, one to eight block columns.
    for (n, mq, k) in [(37, 5, 1), (1025, 17, 3), (5000, 300, 1), (4099, 260, 8)] {
        let q = fixture(n, mq, 11 + n as u64);
        let x = fixture(n, k, 12 + mq as u64);
        let mut naive_s = vec![0.; mq * k];
        for c in 0..k {
            for j in 0..mq {
                naive_s[c * mq + j] = (0..n).map(|i| q[(i, j)] * x[(i, c)]).sum();
            }
        }
        let mut reference: Option<(Vec<u64>, Vec<u64>)> = None;
        for threads in [1, 2, 3, 8] {
            let pool = WorkerPool::new(threads).unwrap();
            let ws = RealGemmWorkspace::with_pool(n, mq.max(k), pool).unwrap();
            let mut s = RealMatrix::zeros(mq, k).unwrap();
            ws.narrow_transpose_into(q.as_slice(), n, mq, x.as_slice(), k, s.as_mut())
                .unwrap();
            let bound = n as f64 * 1e-16;
            for (a, b) in s.as_slice().iter().zip(&naive_s) {
                assert!((a - b).abs() <= bound, "{n}x{mq}x{k}: {a:e} vs {b:e}");
            }
            let mut y = x.clone();
            ws.narrow_subtract(q.as_slice(), n, mq, s.as_ref(), y.as_mut())
                .unwrap();
            // Reference update with the same s: x - Q s (summation order differs).
            let smax = s.as_slice().iter().fold(1_f64, |a, v| a.max(v.abs()));
            for c in 0..k {
                for i in 0..n {
                    let e: f64 = x[(i, c)] - (0..mq).map(|j| q[(i, j)] * s[(j, c)]).sum::<f64>();
                    assert!((y[(i, c)] - e).abs() <= mq as f64 * 1e-16 * smax);
                }
            }
            let now = (bits(s.as_slice()), bits(y.as_slice()));
            match &reference {
                Some(r) => assert_eq!(&now, r, "thread identity {n}x{mq}x{k} T={threads}"),
                None => reference = Some(now),
            }
        }
    }
}

#[test]
fn column_groups_match_the_serial_product() {
    for (rows, k, n) in [(1, 1, 1), (95, 95, 96), (407, 407, 408), (300, 120, 33)] {
        let a = fixture(rows, k, 21 + rows as u64);
        let b = fixture(k, n, 22 + n as u64);
        let mut serial = RealMatrix::zeros(rows, n).unwrap();
        let mut ws1 =
            RealGemmWorkspace::with_pool(rows.max(k), k.max(n), WorkerPool::new(1).unwrap())
                .unwrap();
        ws1.multiply_into(a.as_ref(), RealOp::Identity, b.as_ref(), serial.as_mut())
            .unwrap();
        let mut reference: Option<Vec<u64>> = None;
        for threads in [1, 2, 4, 8] {
            let pool = WorkerPool::new(threads).unwrap();
            let mut ws = RealGemmWorkspace::with_pool(rows.max(k), k.max(n), pool).unwrap();
            let mut out = vec![0.; rows * n];
            ws.multiply_column_groups(a.as_slice(), rows, k, rows, b.as_slice(), n, k, &mut out)
                .unwrap();
            if dft_numerics::MKL_BACKEND {
                for (x, y) in out.iter().zip(serial.as_slice()) {
                    assert!((x - y).abs() <= k as f64 * 1e-16);
                }
            } else {
                // Groups are multiples of the four-column micro tile: the faer
                // product is the serial one bit for bit.
                assert_eq!(bits(&out), bits(serial.as_slice()), "{rows}x{k}x{n}");
            }
            match &reference {
                Some(r) => assert_eq!(&bits(&out), r, "thread identity T={threads}"),
                None => reference = Some(bits(&out)),
            }
        }
    }
}

fn median_ms(mut f: impl FnMut(), reps: usize) -> f64 {
    let mut t: Vec<f64> = (0..reps)
        .map(|_| {
            let s = std::time::Instant::now();
            f();
            s.elapsed().as_secs_f64() * 1e3
        })
        .collect();
    t.sort_by(f64::total_cmp);
    t[reps / 2]
}

/// Kernel timing at the si32 shape (n = 11877 rows, basis m = 408, 68 bands):
/// panel GEMM against the narrow kernels for the projection (k = 1, 2, 4, 8),
/// the reduced border, and the pair rotation plus residual against the fused
/// pass. Observation only.
/// `cargo test --release -p dft-numerics --test algebra_stream stream_timing -- --ignored --nocapture`
#[test]
#[ignore = "timing observation, release"]
fn stream_timing() {
    let (n, m, bands) = (11877, 408, 68);
    let q = fixture(n, m, 1);
    let hq = fixture(n, m, 2);
    let c = fixture(m, bands, 3);
    let values: Vec<f64> = (0..bands).map(|j| 0.01 * j as f64).collect();
    for threads in [1, 8] {
        let pool = WorkerPool::new(threads).unwrap();
        let mut ws = RealGemmWorkspace::with_pool(n, m, pool).unwrap();
        println!("T={threads}");
        for k in [1, 2, 4, 8, 16] {
            let x = fixture(n, k, 4);
            let mut s = RealMatrix::zeros(m, k).unwrap();
            let mut p = RealMatrix::zeros(n, k).unwrap();
            let mut y = x.clone();
            let gemm = median_ms(
                || {
                    ws.multiply_slices_into(
                        dft_numerics::GemmSite::ProjectionOverlap,
                        q.as_slice(),
                        n,
                        m,
                        n,
                        RealOp::Transpose,
                        x.as_slice(),
                        k,
                        n,
                        s.as_mut(),
                    )
                    .unwrap();
                    ws.multiply_slices_into(
                        dft_numerics::GemmSite::ProjectionApply,
                        q.as_slice(),
                        n,
                        m,
                        n,
                        RealOp::Identity,
                        s.as_slice(),
                        k,
                        m,
                        p.as_mut(),
                    )
                    .unwrap();
                    for (a, b) in y.as_mut_slice().iter_mut().zip(p.as_slice()) {
                        *a -= b;
                    }
                },
                15,
            );
            let narrow = median_ms(
                || {
                    ws.narrow_transpose_into(q.as_slice(), n, m, x.as_slice(), k, s.as_mut())
                        .unwrap();
                    ws.narrow_subtract(q.as_slice(), n, m, s.as_ref(), y.as_mut())
                        .unwrap();
                },
                15,
            );
            println!("  projection k={k}: panel GEMM {gemm:.3} ms, narrow {narrow:.3} ms");
        }
        let hx = fixture(n, 1, 5);
        let mut border = RealMatrix::zeros(m, 1).unwrap();
        let gemm = median_ms(
            || {
                ws.multiply_slices_into(
                    dft_numerics::GemmSite::Reduced,
                    q.as_slice(),
                    n,
                    m,
                    n,
                    RealOp::Transpose,
                    hx.as_slice(),
                    1,
                    n,
                    border.as_mut(),
                )
                .unwrap()
            },
            15,
        );
        let narrow = median_ms(
            || {
                ws.narrow_transpose_into(q.as_slice(), n, m, hx.as_slice(), 1, border.as_mut())
                    .unwrap()
            },
            15,
        );
        println!("  reduced border k=1: panel GEMM {gemm:.3} ms, narrow {narrow:.3} ms");
        let (mut xv, mut xh, mut r) = (
            vec![0.; n * bands],
            vec![0.; n * bands],
            vec![0.; n * bands],
        );
        let pair = median_ms(
            || {
                ws.rotate_pair_slices(
                    q.as_slice(),
                    hq.as_slice(),
                    m,
                    c.as_slice(),
                    bands,
                    &mut xv,
                    &mut xh,
                )
                .unwrap();
                for (j, &l) in values.iter().enumerate() {
                    for i in j * n..(j + 1) * n {
                        r[i] = xh[i] - l * xv[i];
                    }
                }
            },
            7,
        );
        let fused = median_ms(
            || {
                ws.rotate_residual_slices(
                    q.as_slice(),
                    hq.as_slice(),
                    m,
                    c.as_slice(),
                    bands,
                    &values,
                    &mut xv,
                    &mut r,
                )
                .unwrap();
            },
            7,
        );
        println!("  rotate+residual 68 bands: pair {pair:.3} ms, fused {fused:.3} ms");
    }
}
