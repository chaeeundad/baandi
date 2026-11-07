use dft_numerics::{
    pool::WorkerPool,
    real::{RealGemmWorkspace, RealMatrix, RealOp},
};
fn matrix(n: usize, m: usize, phase: usize) -> RealMatrix {
    RealMatrix::new(
        n,
        m,
        (0..n * m).map(|i| ((i * 7 + phase) as f64).sin()).collect(),
    )
    .unwrap()
}
#[test]
fn paired_rotation_preserves_every_panel_bit() {
    for threads in [1, 2, 4, 8, 12] {
        let pool = WorkerPool::new(threads).unwrap();
        for (n, k, m) in [
            (0, 3, 2),
            (1, 3, 2),
            (255, 7, 5),
            (256, 7, 5),
            (257, 7, 5),
            (1025, 40, 27),
            (2969, 40, 20),
            (1025, 0, 3),
        ] {
            let v = matrix(n, k, 1);
            let hv = matrix(n, k, 3);
            let c = matrix(k, m, 7);
            let mut ws = RealGemmWorkspace::with_pool(n, k.max(m), pool.clone()).unwrap();
            let mut a = matrix(n, m, 0);
            let mut b = a.clone();
            ws.multiply_into(v.as_ref(), RealOp::Identity, c.as_ref(), a.as_mut())
                .unwrap();
            ws.multiply_into(hv.as_ref(), RealOp::Identity, c.as_ref(), b.as_mut())
                .unwrap();
            let mut x = a.clone();
            let mut y = b.clone();
            x.as_mut_slice().fill(f64::NAN);
            y.as_mut_slice().fill(f64::NAN);
            // Exercise the resident call path as in SCF.
            pool.install(|| {
                ws.rotate_pair(v.as_ref(), hv.as_ref(), c.as_ref(), x.as_mut(), y.as_mut())
            })
            .unwrap();
            for (got, expected) in x
                .as_slice()
                .iter()
                .chain(y.as_slice())
                .zip(a.as_slice().iter().chain(b.as_slice()))
            {
                assert_eq!(got.to_bits(), expected.to_bits());
            }
        }
    }
}
