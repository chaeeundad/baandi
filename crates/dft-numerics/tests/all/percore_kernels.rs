//! PERCORE-KERNELS (docs/performance/percore-kernels-2026-11-17.md): kernel
//! microbenchmark of the lean Davidson algebra at the per-iteration shapes of
//! the stage-0 si32 / si64 solves, and the bit-identity checks of the two fused
//! forms it measures.
//!
//! * Border: the overlap border V^T E and the reduced-matrix border V^T (H E) as
//!   two transposed products of width k, against one transposed product of
//!   width 2k on the contiguous block [E | H E] (copies included). On the
//!   aarch64 pure kernel every output entry is one ascending fused multiply-add
//!   chain inside its fixed 256-row panel, so the fused product gives the same
//!   bits (asserted below); MKL and x86 faer are rounding class.
//! * Ritz: V C_s, (H V) C_s and the residual H X_s - X_s diag(theta) by the pair
//!   rotation plus the elementwise residual pass, against the existing fused
//!   one-pass kernel `rotate_residual_slices` (bit-identical, see
//!   tests/all/algebra_stream.rs).
//!
//! Dimensionless synthetic matrices; timings are Mac observations, never
//! acceptance evidence.
use dft_numerics::{
    GemmSite,
    pool::WorkerPool,
    real::{RealGemmWorkspace, RealMatrix, RealOp},
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

/// Two border products (separate) or one on the packed block (fused).
#[allow(clippy::too_many_arguments)]
fn border(
    ws: &mut RealGemmWorkspace,
    v: &RealMatrix,
    e: &[f64],
    he: &[f64],
    k: usize,
    pair: &mut [f64],
    fused: bool,
    s: &mut RealMatrix,
    h: &mut RealMatrix,
    sh: &mut RealMatrix,
) {
    let (n, m) = (v.rows(), v.cols());
    if fused {
        pair[..n * k].copy_from_slice(e);
        pair[n * k..2 * n * k].copy_from_slice(he);
        sh.reshape_within_storage(m, 2 * k).unwrap();
        ws.multiply_slices_into(
            GemmSite::Gram,
            v.as_slice(),
            n,
            m,
            n,
            RealOp::Transpose,
            &pair[..2 * n * k],
            2 * k,
            n,
            sh.as_mut(),
        )
        .unwrap();
    } else {
        for (b, out) in [(e, &mut *s), (he, &mut *h)] {
            out.reshape_within_storage(m, k).unwrap();
            ws.multiply_slices_into(
                GemmSite::Gram,
                v.as_slice(),
                n,
                m,
                n,
                RealOp::Transpose,
                b,
                k,
                n,
                out.as_mut(),
            )
            .unwrap();
        }
    }
}

#[test]
fn fused_border_matches_two_products() {
    // Bit identity holds on the aarch64 pure kernel (one FMA chain per entry,
    // independent of the output width); elsewhere the fused product is checked
    // to the rounding bound of a reordered sum.
    let exact = cfg!(target_arch = "aarch64") && !dft_numerics::MKL_BACKEND;
    for (n, m, k) in [(300, 7, 1), (1025, 40, 3), (2049, 64, 17), (777, 33, 33)] {
        let v = fixture(n, m, 1 + n as u64);
        let e = fixture(n, k, 2 + m as u64);
        let he = fixture(n, k, 3 + k as u64);
        let mut reference: Option<Vec<u64>> = None;
        for threads in [1, 2, 4] {
            let pool = WorkerPool::new(threads).unwrap();
            let mut ws = RealGemmWorkspace::with_pool(n, m.max(2 * k), pool).unwrap();
            let mut pair = vec![0.; 2 * n * k];
            let mut s = RealMatrix::zeros(m, k).unwrap();
            let mut h = RealMatrix::zeros(m, k).unwrap();
            let mut sh = RealMatrix::zeros(m, 2 * k).unwrap();
            border(
                &mut ws,
                &v,
                e.as_slice(),
                he.as_slice(),
                k,
                &mut pair,
                false,
                &mut s,
                &mut h,
                &mut sh,
            );
            border(
                &mut ws,
                &v,
                e.as_slice(),
                he.as_slice(),
                k,
                &mut pair,
                true,
                &mut s,
                &mut h,
                &mut sh,
            );
            let (fs, fh) = sh.as_slice().split_at(m * k);
            if exact {
                assert_eq!(bits(fs), bits(s.as_slice()), "S border n={n} m={m} k={k}");
                assert_eq!(bits(fh), bits(h.as_slice()), "H border n={n} m={m} k={k}");
            } else {
                for (x, y) in fs
                    .iter()
                    .chain(fh)
                    .zip(s.as_slice().iter().chain(h.as_slice()))
                {
                    assert!((x - y).abs() <= 1e-13 * n as f64, "n={n} m={m} k={k}");
                }
            }
            let got = bits(sh.as_slice());
            match &reference {
                Some(r) => assert_eq!(r, &got, "thread budget {threads}"),
                None => reference = Some(got),
            }
        }
    }
}

fn median_us(mut f: impl FnMut(), reps: usize) -> f64 {
    f();
    let mut t: Vec<f64> = (0..7)
        .map(|_| {
            let start = std::time::Instant::now();
            for _ in 0..reps {
                f();
            }
            start.elapsed().as_secs_f64() * 1e6 / reps as f64
        })
        .collect();
    t.sort_by(f64::total_cmp);
    t[3]
}

#[test]
#[ignore = "release microbenchmark: cargo test --release -p dft-numerics --test all -- --ignored --nocapture percore_kernels::percore_kernel_microbench"]
fn percore_kernel_microbench() {
    let threads: usize = std::env::var("PERCORE_BENCH_THREADS")
        .ok()
        .and_then(|t| t.parse().ok())
        .unwrap_or(1);
    let pool = WorkerPool::new(threads).unwrap();
    println!(
        "threads {threads}, MKL backend {}",
        dft_numerics::MKL_BACKEND
    );
    // (case, n, basis columns m (after the append, borders), block k)
    // from the stage-0 traces: si32 warm m 68 -> 244, si64 warm m 132 -> ...
    let shapes: [(&str, usize, usize, usize); 12] = [
        ("si32 wide", 11877, 136, 68),
        ("si32 wide2", 11877, 200, 64),
        ("si32 narrow", 11877, 210, 10),
        ("si32 narrow", 11877, 224, 7),
        ("si32 narrow", 11877, 237, 5),
        ("si32 tail", 11877, 244, 2),
        ("si32 tail", 11877, 219, 1),
        ("si64 wide", 23871, 264, 132),
        ("si64 wide2", 23871, 396, 124),
        ("si64 narrow", 23871, 420, 12),
        ("si64 narrow", 23871, 440, 6),
        ("si64 tail", 23871, 460, 2),
    ];
    println!("BORDER  two TN products of width k vs one of width 2k (copies included)");
    println!(
        "case          n      m    k   sep_us  sep_GF/s  fused_us  fused_GF/s  fused/sep  AI_sep  AI_fused"
    );
    for (name, n, m, k) in shapes {
        let v = fixture(n, m, 1);
        let e = fixture(n, k, 2);
        let he = fixture(n, k, 3);
        let mut ws = RealGemmWorkspace::with_pool(n, m.max(2 * k), pool.clone()).unwrap();
        let mut pair = vec![0.; 2 * n * k];
        let mut s = RealMatrix::zeros(m, k).unwrap();
        let mut h = RealMatrix::zeros(m, k).unwrap();
        let mut sh = RealMatrix::zeros(m, 2 * k).unwrap();
        let flops = 2. * 2. * (n * m * k) as f64;
        let reps = ((1e9 / flops) as usize).clamp(2, 200);
        let sep = median_us(
            || {
                border(
                    &mut ws,
                    &v,
                    e.as_slice(),
                    he.as_slice(),
                    k,
                    &mut pair,
                    false,
                    &mut s,
                    &mut h,
                    &mut sh,
                )
            },
            reps,
        );
        let fused = median_us(
            || {
                border(
                    &mut ws,
                    &v,
                    e.as_slice(),
                    he.as_slice(),
                    k,
                    &mut pair,
                    true,
                    &mut s,
                    &mut h,
                    &mut sh,
                )
            },
            reps,
        );
        // Arithmetic intensity (flop per byte of compulsory traffic): separate
        // reads V twice plus E and HE; fused reads V once plus E, HE, and the
        // copy (read + write of 2 n k).
        let bytes_sep = 8. * (2 * n * m + 2 * n * k + 2 * m * k) as f64;
        let bytes_fused = 8. * (n * m + 2 * n * k + 4 * n * k + 2 * m * k) as f64;
        println!(
            "{name:<11} {n:>6} {m:>5} {k:>4} {sep:>8.0} {:>9.1} {fused:>9.0} {:>11.1} {:>10.3} {:>7.2} {:>9.2}",
            flops / sep / 1e3,
            flops / fused / 1e3,
            fused / sep,
            flops / bytes_sep,
            flops / bytes_fused
        );
    }
    println!("RITZ  pair rotation + residual pass vs fused one-pass kernel (loop site)");
    println!("case          n      m    k   pair_us  pair_GF/s  fused_us  fused_GF/s  fused/pair");
    for (name, n, m, k) in shapes {
        let v = fixture(n, m, 4);
        let hv = fixture(n, m, 5);
        let c = fixture(m, k, 6);
        let values: Vec<f64> = (0..k).map(|j| 0.01 * j as f64).collect();
        let mut ws = RealGemmWorkspace::with_pool(n, m.max(2 * k), pool.clone()).unwrap();
        let (mut x, mut hx, mut r) = (vec![0.; n * k], vec![0.; n * k], vec![0.; n * k]);
        let mut norms = vec![0.; k];
        let flops = 2. * 2. * (n * m * k) as f64;
        let reps = ((1e9 / flops) as usize).clamp(2, 200);
        let pair = median_us(
            || {
                ws.rotate_pair_slices(
                    v.as_slice(),
                    hv.as_slice(),
                    m,
                    c.as_slice(),
                    k,
                    &mut x,
                    &mut hx,
                )
                .unwrap();
                for j in 0..k {
                    let mut sum = 0.;
                    for ((r, &h), &x) in r[j * n..(j + 1) * n]
                        .iter_mut()
                        .zip(&hx[j * n..(j + 1) * n])
                        .zip(&x[j * n..(j + 1) * n])
                    {
                        *r = h - values[j] * x;
                        sum += *r * *r;
                    }
                    norms[j] = sum.sqrt();
                }
            },
            reps,
        );
        let fused = median_us(
            || {
                assert!(
                    ws.rotate_residual_slices(
                        v.as_slice(),
                        hv.as_slice(),
                        m,
                        c.as_slice(),
                        k,
                        &values,
                        &mut x,
                        &mut r,
                    )
                    .unwrap()
                );
                for j in 0..k {
                    norms[j] = r[j * n..(j + 1) * n]
                        .iter()
                        .map(|x| x * x)
                        .sum::<f64>()
                        .sqrt();
                }
            },
            reps,
        );
        println!(
            "{name:<11} {n:>6} {m:>5} {k:>4} {pair:>9.0} {:>10.1} {fused:>9.0} {:>11.1} {:>11.3}",
            flops / pair / 1e3,
            flops / fused / 1e3,
            fused / pair
        );
    }
}
