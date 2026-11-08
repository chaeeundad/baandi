//! Blocked rotation kernel: agreement with the legacy kernel, independence of
//! the thread budget and row partition, and an ignored timing comparison.
use dft_numerics::{
    pool::WorkerPool,
    real::{RealGemmWorkspace, RealMatrix, RotateKernel},
};
fn matrix(n: usize, m: usize, phase: usize) -> RealMatrix {
    RealMatrix::new(
        n,
        m,
        (0..n * m)
            .map(|i| ((i * 7 + phase) as f64).sin() * (1.0 + (i % 13) as f64 * 0.01))
            .collect(),
    )
    .unwrap()
}
fn rotate(
    pool: &std::sync::Arc<WorkerPool>,
    kernel: RotateKernel,
    shape: (usize, usize, usize),
) -> (RealMatrix, RealMatrix) {
    let (n, k, m) = shape;
    let v = matrix(n, k, 1);
    let hv = matrix(n, k, 3);
    let c = matrix(k, m, 7);
    let mut ws = RealGemmWorkspace::with_pool(n, k.max(m), pool.clone()).unwrap();
    ws.set_rotate_kernel(kernel);
    let mut x = matrix(n, m, 0);
    let mut y = x.clone();
    pool.install(|| ws.rotate_pair(v.as_ref(), hv.as_ref(), c.as_ref(), x.as_mut(), y.as_mut()))
        .unwrap();
    (x, y)
}
#[test]
fn blocked_matches_legacy_and_is_thread_independent() {
    let shapes = [
        (1, 3, 2),
        (7, 5, 3),
        (8, 5, 4),
        (9, 6, 5),
        (255, 7, 5),
        (257, 40, 27),
        (1025, 40, 27),
        (2969, 126, 36),
    ];
    for shape in shapes {
        let reference = rotate(&WorkerPool::new(1).unwrap(), RotateKernel::Legacy, shape);
        let blocked1 = rotate(&WorkerPool::new(1).unwrap(), RotateKernel::Blocked, shape);
        let scale = reference
            .0
            .as_slice()
            .iter()
            .fold(0_f64, |a, &b| a.max(b.abs()))
            .max(1.0);
        let k = shape.1 as f64;
        for (got, want) in blocked1
            .0
            .as_slice()
            .iter()
            .chain(blocked1.1.as_slice())
            .zip(reference.0.as_slice().iter().chain(reference.1.as_slice()))
        {
            assert!(
                (got - want).abs() <= 4.0 * k * f64::EPSILON * scale,
                "{shape:?}"
            );
        }
        for threads in [2, 4, 8, 12] {
            let other = rotate(
                &WorkerPool::new(threads).unwrap(),
                RotateKernel::Blocked,
                shape,
            );
            for (a, b) in other
                .0
                .as_slice()
                .iter()
                .chain(other.1.as_slice())
                .zip(blocked1.0.as_slice().iter().chain(blocked1.1.as_slice()))
            {
                assert_eq!(a.to_bits(), b.to_bits(), "{shape:?} threads {threads}");
            }
        }
    }
}
/// `cargo test --release -p dft-numerics --test real_rotate_kernel timing -- --ignored --nocapture`
/// Observation only: legacy versus blocked rotation at Davidson shapes, 1 thread.
#[test]
#[ignore]
fn timing() {
    let pool = WorkerPool::new(1).unwrap();
    for (label, shape) in [
        ("si8 ", (1480, 96, 20)),
        ("si16", (2958, 126, 36)),
        ("si16 s=216", (2958, 216, 36)),
        ("si32", (5939, 238, 68)),
        ("si32 s=408", (5939, 408, 68)),
    ] {
        let (n, k, m) = shape;
        let v = matrix(n, k, 1);
        let hv = matrix(n, k, 3);
        let c = matrix(k, m, 7);
        let mut out = Vec::new();
        for kernel in [RotateKernel::Legacy, RotateKernel::Blocked] {
            let mut ws = RealGemmWorkspace::with_pool(n, k.max(m), pool.clone()).unwrap();
            ws.set_rotate_kernel(kernel);
            let mut x = matrix(n, m, 0);
            let mut y = x.clone();
            let mut best = f64::INFINITY;
            for _ in 0..9 {
                let t = std::time::Instant::now();
                for _ in 0..4 {
                    pool.install(|| {
                        ws.rotate_pair(v.as_ref(), hv.as_ref(), c.as_ref(), x.as_mut(), y.as_mut())
                    })
                    .unwrap();
                }
                best = best.min(t.elapsed().as_secs_f64() / 4.0);
            }
            let gflops = 2.0 * 2.0 * n as f64 * k as f64 * m as f64 / best / 1e9;
            out.push(format!(
                "{kernel:?} {:.3} ms {:.2} GFLOP/s",
                best * 1e3,
                gflops
            ));
        }
        println!("{label} (n={n}, k={k}, b={m}): {}", out.join(" | "));
    }
}
