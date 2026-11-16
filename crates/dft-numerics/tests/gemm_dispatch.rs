//! F11 dispatch counters (`dft_numerics::gemm_dispatch`): every real GEMM
//! wrapper call is counted once under its effective route with its shape and
//! 2 m n k flops; a pure fallback inside a wrapper is not counted again; the
//! counted calls do not change the product (same bits as the uncounted kernel
//! path). One test function: the counters are process-wide.
use dft_numerics::gemm_dispatch::{DispatchRow, rows, rows_as_if_mkl, wide_requests};
use dft_numerics::pool::WorkerPool;
use dft_numerics::real::{GemmSite, RealGemmWorkspace, RealMatrix, RealOp};

fn fixture(r: usize, c: usize, seed: usize) -> RealMatrix {
    RealMatrix::new(
        r,
        c,
        (0..r * c)
            .map(|i| ((i * 13 + seed) as f64 * 0.31).sin())
            .collect(),
    )
    .unwrap()
}
fn row(site: &str, route: &str) -> Option<DispatchRow> {
    rows()
        .into_iter()
        .find(|r| r.site == site && r.route == route)
}

#[test]
fn dispatch_counters_follow_the_route() {
    // short^2 = 4096: a 256-row panel of the 11-column output fits a tile,
    // the whole 600 x 11 output does not.
    let (n, short) = (600, 64);
    let pool = WorkerPool::new(2).unwrap();
    let mut ws = RealGemmWorkspace::with_pool(n, short, pool).unwrap();
    let a = fixture(n, 10, 1);
    let b = fixture(n, 7, 2);
    // A^T B through the slice entry (gram site): output 10 x 7, reduction n.
    let mut out = RealMatrix::zeros(10, 7).unwrap();
    ws.multiply_slices_into(
        GemmSite::Gram,
        a.as_slice(),
        n,
        10,
        n,
        RealOp::Transpose,
        b.as_slice(),
        7,
        n,
        out.as_mut(),
    )
    .unwrap();
    // The same product through the direct pure entry: same bits.
    let mut direct = RealMatrix::zeros(10, 7).unwrap();
    ws.multiply_into(a.as_ref(), RealOp::Transpose, b.as_ref(), direct.as_mut())
        .unwrap();
    if !dft_numerics::MKL_BACKEND {
        assert_eq!(out.as_slice(), direct.as_slice());
    }
    let gram = if dft_numerics::MKL_BACKEND {
        row("gram", "mkl")
    } else {
        row("gram", "pure_build")
    }
    .expect("gram row");
    assert_eq!(gram.calls, 1);
    assert_eq!(gram.flops, 2 * 10 * 7 * n as u64);
    assert_eq!((gram.max_m, gram.max_n, gram.max_k), (10, 7, n as u64));
    // The direct call is counted once, as pure_direct (no second count from
    // the slice entry's fallback).
    let d = row("unlabelled", "pure_direct").expect("direct row");
    assert_eq!(d.calls, 1);
    // A wide apply Q S (output n x 11 = 6600 > short^2 entries).
    let s = fixture(10, 11, 3);
    let mut long = RealMatrix::zeros(n, 11).unwrap();
    ws.multiply_slices_wide_into(
        GemmSite::ProjectionApply,
        a.as_slice(),
        n,
        10,
        n,
        RealOp::Identity,
        s.as_slice(),
        11,
        10,
        long.as_mut(),
    )
    .unwrap();
    let wide = rows()
        .into_iter()
        .filter(|r| r.site == "projection_apply")
        .collect::<Vec<_>>();
    assert_eq!(wide.len(), 1);
    assert_eq!(
        wide[0].route,
        if dft_numerics::MKL_BACKEND {
            "mkl_wide"
        } else {
            "pure_build"
        }
    );
    assert_eq!(wide_requests(), vec![("projection_apply", 1)]);
    // The same shape through the square-bound entry falls back on MKL builds.
    ws.multiply_slices_into(
        GemmSite::ProjectionApply,
        a.as_slice(),
        n,
        10,
        n,
        RealOp::Identity,
        s.as_slice(),
        11,
        10,
        long.as_mut(),
    )
    .unwrap();
    let route = if dft_numerics::MKL_BACKEND {
        "pure_output_bound"
    } else {
        "pure_build"
    };
    let r = row("projection_apply", route).unwrap();
    assert!(r.calls >= 1);
    assert!(!r.is_mkl());
    // Pure builds also count the route an MKL build would take.
    if !dft_numerics::MKL_BACKEND {
        let as_if = rows_as_if_mkl();
        let find = |route: &str| {
            as_if
                .iter()
                .find(|r| r.site == "projection_apply" && r.route == route)
                .map(|r| r.calls)
        };
        assert_eq!(find("mkl_wide"), Some(1));
        assert_eq!(find("pure_output_bound"), Some(1));
        assert_eq!(
            as_if.iter().find(|r| r.site == "gram").map(|r| r.route),
            Some("mkl")
        );
    } else {
        assert!(rows_as_if_mkl().is_empty());
    }
    // Rotation pair: two products counted under rotate.
    let c = fixture(10, 4, 4);
    let mut xv = vec![0.; n * 4];
    let mut xh = vec![0.; n * 4];
    ws.rotate_pair_slices(
        a.as_slice(),
        a.as_slice(),
        10,
        c.as_slice(),
        4,
        &mut xv,
        &mut xh,
    )
    .unwrap();
    let rot = rows()
        .into_iter()
        .filter(|r| r.site == "rotate")
        .map(|r| r.calls)
        .sum::<u64>();
    assert_eq!(rot, 2);
}
