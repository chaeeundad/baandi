//! P48B packed pure-Rust GEMM (`BAANDI_PURE_GEMM`): agreement with the nano
//! kernel within the rounding class, bit identity across thread budgets, the
//! exact ascending fused multiply-add chain of the SIMD kernel, and an ignored
//! per-site microbenchmark (release command in its attribute).
use dft_numerics::{
    pool::WorkerPool,
    real::{PureGemm, RealGemmWorkspace, RealMatrix, RealOp},
};
use std::sync::Arc;

const KERNELS: [PureGemm; 3] = [PureGemm::Packed, PureGemm::Simd, PureGemm::Library];

fn random(rows: usize, cols: usize, seed: u64) -> RealMatrix {
    let mut state = seed.wrapping_mul(0x9e37_79b9_7f4a_7c15) | 1;
    RealMatrix::new(
        rows,
        cols,
        (0..rows * cols)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                (state >> 11) as f64 / (1u64 << 53) as f64 - 0.5
            })
            .collect(),
    )
    .unwrap()
}
/// Ill-conditioned: magnitudes over 16 decades and columns of the right operand
/// built so that sums cancel (each entry repeated with opposite sign).
fn ill(rows: usize, cols: usize, seed: u64) -> RealMatrix {
    let base = random(rows, cols, seed);
    let data = base
        .as_slice()
        .iter()
        .enumerate()
        .map(|(i, &x)| {
            let scale = 10f64.powi((i % 17) as i32 - 8);
            let sign = if (i / rows).is_multiple_of(2) {
                1.
            } else {
                -1.
            };
            sign * scale * (1. + x * 1e-9)
        })
        .collect();
    RealMatrix::new(rows, cols, data).unwrap()
}

/// Per entry: |x - y| <= 2 k eps sum_p |a_ip b_pj|, the bound for two
/// summation orders of the same k products. Returns the largest ratio
/// |x - y| / (eps sum |a||b|) (the difference in units of the error scale).
fn check(
    got: &RealMatrix,
    want: &RealMatrix,
    a: &RealMatrix,
    transpose: bool,
    b: &RealMatrix,
) -> f64 {
    let (m, n) = (want.rows(), want.cols());
    let k = if transpose { a.rows() } else { a.cols() };
    let mut worst = 0_f64;
    for j in 0..n {
        for i in 0..m {
            let scale: f64 = (0..k)
                .map(|p| {
                    let x = if transpose { a[(p, i)] } else { a[(i, p)] };
                    (x * b[(p, j)]).abs()
                })
                .sum();
            let diff = (got[(i, j)] - want[(i, j)]).abs();
            let bound = 2. * k as f64 * f64::EPSILON * scale;
            assert!(
                diff <= bound,
                "entry ({i},{j}) differs by {diff:e} > bound {bound:e} (k {k})"
            );
            if scale > 0. {
                worst = worst.max(diff / (f64::EPSILON * scale));
            }
        }
    }
    worst
}

fn workspace(threads: usize, long: usize, short: usize, kernel: PureGemm) -> RealGemmWorkspace {
    let mut ws =
        RealGemmWorkspace::with_pool(long, short, WorkerPool::new(threads).unwrap()).unwrap();
    ws.set_pure_gemm(kernel);
    ws
}

fn product(
    threads: usize,
    kernel: PureGemm,
    a: &RealMatrix,
    op: RealOp,
    b: &RealMatrix,
) -> RealMatrix {
    let m = if op == RealOp::Transpose {
        a.cols()
    } else {
        a.rows()
    };
    let short = a.cols().max(b.cols());
    let mut ws = workspace(threads, a.rows(), short, kernel);
    let mut out = RealMatrix::zeros(m, b.cols()).unwrap();
    ws.multiply_into(a.as_ref(), op, b.as_ref(), out.as_mut())
        .unwrap();
    out
}

const SHAPES: [(usize, usize, usize); 9] = [
    (1, 1, 1),
    (7, 5, 3),
    (9, 6, 5),
    (255, 7, 5),
    (257, 40, 27),
    (300, 270, 9),
    (1025, 40, 27),
    (2969, 126, 36),
    (3001, 33, 1),
];

#[test]
fn packed_agrees_with_nano_within_rounding_class() {
    let mut worst = [0_f64; 3];
    for (s, &(rows, k, n)) in SHAPES.iter().enumerate() {
        for ill_conditioned in [false, true] {
            let make = if ill_conditioned { ill } else { random };
            // A B: rows x k times k x n.
            let a = make(rows, k, 11 + s as u64);
            let b = make(k, n, 23 + s as u64);
            let want = product(1, PureGemm::Nano, &a, RealOp::Identity, &b);
            // A^T B: rows x k (transposed) times rows x n.
            let bt = make(rows, n, 31 + s as u64);
            let want_t = product(1, PureGemm::Nano, &a, RealOp::Transpose, &bt);
            for (w, kernel) in KERNELS.into_iter().enumerate() {
                let got = product(1, kernel, &a, RealOp::Identity, &b);
                worst[w] = worst[w].max(check(&got, &want, &a, false, &b));
                let got = product(1, kernel, &a, RealOp::Transpose, &bt);
                worst[w] = worst[w].max(check(&got, &want_t, &a, true, &bt));
            }
        }
    }
    println!(
        "largest |packed - nano| / (eps sum|a||b|): packed {:.3}, simd {:.3}, library {:.3}",
        worst[0], worst[1], worst[2]
    );
}

#[test]
fn simd_entries_are_one_ascending_fma_chain() {
    for &(rows, k, n) in &[(9, 600, 6), (257, 300, 5), (40, 513, 13)] {
        let a = random(rows, k, 5);
        let b = random(k, n, 6);
        let got = product(1, PureGemm::Simd, &a, RealOp::Identity, &b);
        let bt = random(rows, n, 7);
        let got_t = product(1, PureGemm::Simd, &a, RealOp::Transpose, &bt);
        for j in 0..n {
            for i in 0..rows {
                let chain = (0..k).fold(0_f64, |acc, p| a[(i, p)].mul_add(b[(p, j)], acc));
                assert_eq!(got[(i, j)].to_bits(), chain.to_bits(), "A B ({i},{j})");
            }
            for i in 0..k {
                // Gram: per 256-term panel one chain, panels summed in ascending order.
                let mut total = 0_f64;
                for (t, start) in (0..rows).step_by(256).enumerate() {
                    let end = (start + 256).min(rows);
                    let part =
                        (start..end).fold(0_f64, |acc, p| a[(p, i)].mul_add(bt[(p, j)], acc));
                    total = if t == 0 { 0. + part } else { total + part };
                }
                assert_eq!(got_t[(i, j)].to_bits(), total.to_bits(), "A^T B ({i},{j})");
            }
        }
    }
}

fn bits(m: &RealMatrix) -> Vec<u64> {
    m.as_slice().iter().map(|x| x.to_bits()).collect()
}

/// Every entry point of the workspace, run with one thread budget.
fn all_sites(threads: usize, kernel: PureGemm) -> Vec<Vec<u64>> {
    let (rows, k, n) = (2969, 126, 36);
    let v = random(rows, k, 1);
    let hv = random(rows, k, 2);
    let c = random(k, n, 3);
    let mut ws = workspace(threads, rows, k, kernel);
    let mut out = Vec::new();
    // Rotation pair (slice form).
    let mut x = RealMatrix::zeros(rows, n).unwrap();
    let mut y = RealMatrix::zeros(rows, n).unwrap();
    ws.rotate_pair_slices(
        v.as_slice(),
        hv.as_slice(),
        k,
        c.as_slice(),
        n,
        x.as_mut_slice(),
        y.as_mut_slice(),
    )
    .unwrap();
    out.push(bits(&x));
    out.push(bits(&y));
    // Fused Ritz and residual.
    let values: Vec<f64> = (0..n).map(|j| j as f64 * 0.1 - 1.).collect();
    let mut r = RealMatrix::zeros(rows, n).unwrap();
    assert!(
        ws.rotate_residual_slices(
            v.as_slice(),
            hv.as_slice(),
            k,
            c.as_slice(),
            n,
            &values,
            x.as_mut_slice(),
            r.as_mut_slice(),
        )
        .unwrap()
    );
    out.push(bits(&x));
    out.push(bits(&r));
    // Gram and projection overlap (A^T B with leading dimension).
    let mut g = RealMatrix::zeros(k, n).unwrap();
    ws.multiply_slices_into(
        dft_numerics::real::GemmSite::Gram,
        v.as_slice(),
        rows,
        k,
        rows,
        RealOp::Transpose,
        hv.as_slice(),
        n,
        rows,
        g.as_mut(),
    )
    .unwrap();
    out.push(bits(&g));
    // Projection apply (A B with leading dimension).
    let mut p = RealMatrix::zeros(rows, n).unwrap();
    ws.multiply_slices_into(
        dft_numerics::real::GemmSite::ProjectionApply,
        v.as_slice(),
        rows,
        k,
        rows,
        RealOp::Identity,
        c.as_slice(),
        n,
        k,
        p.as_mut(),
    )
    .unwrap();
    out.push(bits(&p));
    // Column groups.
    let mut q = RealMatrix::zeros(rows, n).unwrap();
    ws.multiply_column_groups(
        v.as_slice(),
        rows,
        k,
        rows,
        c.as_slice(),
        n,
        k,
        q.as_mut_slice(),
    )
    .unwrap();
    out.push(bits(&q));
    out
}

#[test]
fn packed_is_bit_identical_across_thread_budgets() {
    for kernel in KERNELS {
        let reference = all_sites(1, kernel);
        for threads in [2, 4, 8] {
            assert!(
                all_sites(threads, kernel) == reference,
                "{kernel:?} differs at {threads} threads"
            );
        }
    }
}

#[test]
fn default_follows_switch_and_is_nano() {
    let ws = RealGemmWorkspace::with_pool(8, 4, WorkerPool::new(1).unwrap()).unwrap();
    // The test process does not set BAANDI_PURE_GEMM.
    assert_eq!(ws.pure_gemm(), PureGemm::Nano);
    let expected = if cfg!(target_arch = "x86_64") {
        PureGemm::Library
    } else {
        PureGemm::Simd
    };
    assert_eq!(PureGemm::Packed.resolve(), expected);
}

/// Per-site shapes (rows = half-sphere plane waves, k = basis columns, n =
/// output columns): si16, si32 (dev08 profile shape 5939 x 408 x 68) and si64
/// estimates; nonlocal P^T X and P W with 8 projectors per atom.
#[test]
#[ignore = "release microbenchmark: cargo test --release -p dft-numerics --test pure_gemm -- --ignored --nocapture pure_gemm_microbench"]
fn pure_gemm_microbench() {
    let threads: usize = std::env::var("PURE_GEMM_BENCH_THREADS")
        .ok()
        .and_then(|t| t.parse().ok())
        .unwrap_or(1);
    let cases: [(&str, RealOp, usize, usize, usize); 12] = [
        ("si16 rotate", RealOp::Identity, 2969, 216, 36),
        ("si32 rotate", RealOp::Identity, 5939, 408, 68),
        ("si64 rotate", RealOp::Identity, 11878, 816, 136),
        ("si16 gram", RealOp::Transpose, 2969, 216, 36),
        ("si32 gram", RealOp::Transpose, 5939, 408, 68),
        ("si64 gram", RealOp::Transpose, 11878, 816, 136),
        ("si32 narrow gram", RealOp::Transpose, 5939, 408, 4),
        ("si32 reduced", RealOp::Transpose, 5939, 68, 68),
        ("si64 reduced", RealOp::Transpose, 11878, 136, 136),
        ("si32 nonlocal PtX", RealOp::Transpose, 5939, 256, 68),
        ("si32 nonlocal PW", RealOp::Identity, 5939, 256, 68),
        ("si64 nonlocal PtX", RealOp::Transpose, 11878, 512, 136),
    ];
    let pool: Arc<WorkerPool> = WorkerPool::new(threads).unwrap();
    println!("threads {threads}; GFLOP/s, median of 7 (2 m n k flops)");
    println!("case                     rows    k    n   nano   packed   simd   library");
    for (name, op, rows, k, n) in cases {
        let a = random(rows, k, 1);
        let b = if op == RealOp::Transpose {
            random(rows, n, 2)
        } else {
            random(k, n, 2)
        };
        let m = if op == RealOp::Transpose { k } else { rows };
        let mut out = RealMatrix::zeros(m, n).unwrap();
        let mut rates = Vec::new();
        for kernel in [
            PureGemm::Nano,
            PureGemm::Packed,
            PureGemm::Simd,
            PureGemm::Library,
        ] {
            let mut ws = RealGemmWorkspace::with_pool(rows, k.max(n), pool.clone()).unwrap();
            ws.set_pure_gemm(kernel);
            let flops = 2. * (rows * k * n) as f64;
            let reps = ((2e9 / flops) as usize).clamp(1, 200);
            ws.multiply_into(a.as_ref(), op, b.as_ref(), out.as_mut())
                .unwrap();
            let mut times: Vec<f64> = (0..7)
                .map(|_| {
                    let start = std::time::Instant::now();
                    for _ in 0..reps {
                        ws.multiply_into(a.as_ref(), op, b.as_ref(), out.as_mut())
                            .unwrap();
                    }
                    start.elapsed().as_secs_f64() / reps as f64
                })
                .collect();
            times.sort_by(f64::total_cmp);
            rates.push(flops / times[3] / 1e9);
        }
        println!(
            "{name:<22} {rows:>6} {k:>4} {n:>4} {:>6.1} {:>8.1} {:>6.1} {:>9.1}",
            rates[0], rates[1], rates[2], rates[3]
        );
    }
}
