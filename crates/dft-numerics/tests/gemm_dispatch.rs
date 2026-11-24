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
    single_products();
}

fn single_products() {
    use dft_numerics::linalg::LinalgError;
    use dft_numerics::switches::{Id, RunConfig};
    let default = RunConfig::from_pairs(&[]).unwrap();
    assert!(!default.on(Id::REF_DAV_SINGLE_GEMM));
    let on = RunConfig::from_pairs(&[("ref_dav_single_gemm", "on")]).unwrap();
    let _scope = on.enter();
    let mut calls = 0;
    let mut flops = 0;
    let mut max_error = 0.0_f64;
    // Includes empty owners, panel boundaries, padded inputs, a strided output,
    // and enough work for the original panel entry to use the worker pool.
    for (k, m, n) in [
        (0, 3, 2),
        (1, 3, 2),
        (255, 7, 5),
        (256, 7, 5),
        (257, 7, 5),
        (1025, 40, 27),
    ] {
        let lda = k + 3;
        let ldb = k + 5;
        let a = fixture(lda, m, 7);
        let b = fixture(ldb, n, 13);
        let mut bits = None;
        for threads in [1, 2, 4, 8] {
            let mut ws =
                RealGemmWorkspace::with_pool(k, 64, WorkerPool::new(threads).unwrap()).unwrap();
            let mut legacy = RealMatrix::zeros(m, n).unwrap();
            ws.multiply_slices_into(
                GemmSite::Other,
                a.as_slice(),
                k,
                m,
                lda,
                RealOp::Transpose,
                b.as_slice(),
                n,
                ldb,
                legacy.as_mut(),
            )
            .unwrap();
            // The new switch alone cannot reroute any old entry.
            let old = legacy.clone();
            {
                let _scope = default.enter();
                ws.multiply_slices_into(
                    GemmSite::Other,
                    a.as_slice(),
                    k,
                    m,
                    lda,
                    RealOp::Transpose,
                    b.as_slice(),
                    n,
                    ldb,
                    legacy.as_mut(),
                )
                .unwrap();
            }
            assert_eq!(
                old.as_slice()
                    .iter()
                    .map(|x| x.to_bits())
                    .collect::<Vec<_>>(),
                legacy
                    .as_slice()
                    .iter()
                    .map(|x| x.to_bits())
                    .collect::<Vec<_>>()
            );
            for _ in 0..2 {
                let mut out = RealMatrix::new(m + 2, n + 2, vec![91.0; (m + 2) * (n + 2)]).unwrap();
                ws.multiply_slices_single_into(
                    GemmSite::ProjectionOverlap,
                    a.as_slice(),
                    k,
                    m,
                    lda,
                    RealOp::Transpose,
                    b.as_slice(),
                    n,
                    ldb,
                    out.as_mut().submatrix_mut(1, 1, m, n),
                )
                .unwrap();
                if k > 0 {
                    calls += 1;
                    flops += 2 * m * n * k;
                }
                let mut actual = Vec::new();
                for j in 0..n + 2 {
                    for i in 0..m + 2 {
                        if i == 0 || i == m + 1 || j == 0 || j == n + 1 {
                            assert_eq!(out[(i, j)], 91.0);
                        } else {
                            let expected: f64 = (0..k).map(|p| a[(p, i - 1)] * b[(p, j - 1)]).sum();
                            max_error = max_error.max((out[(i, j)] - expected).abs());
                            assert!((out[(i, j)] - expected).abs() < 1e-11);
                            assert!((out[(i, j)] - legacy[(i - 1, j - 1)]).abs() < 1e-11);
                            actual.push(out[(i, j)].to_bits());
                        }
                    }
                }
                if let Some(ref prior) = bits {
                    assert_eq!(&actual, prior);
                }
                bits = Some(actual);
            }
        }
    }
    let table = if dft_numerics::MKL_BACKEND {
        rows()
    } else {
        rows_as_if_mkl()
    };
    let single = table
        .iter()
        .find(|r| r.site == "projection_overlap" && r.route == "mkl_single")
        .unwrap();
    assert!(single.is_mkl());
    assert_eq!(single.calls, calls);
    assert_eq!(single.flops, flops as u64);
    assert_eq!((single.max_m, single.max_n, single.max_k), (40, 27, 1025));
    let empty = table
        .iter()
        .find(|r| r.site == "projection_overlap" && r.route == "pure_layout")
        .unwrap();
    assert_eq!((empty.calls, empty.flops), (8, 0));
    println!(
        "single GEMM: {calls} products, {flops} FLOPs, max scalar error {max_error:e} (dimensionless); T=1/2/4/8 repeats bit-identical; MKL executed={}",
        dft_numerics::MKL_BACKEND
    );

    let mut ws = RealGemmWorkspace::with_pool(257, 8, WorkerPool::new(1).unwrap()).unwrap();
    let a = fixture(257, 3, 1);
    let b = fixture(257, 2, 2);
    let mut out = RealMatrix::zeros(3, 2).unwrap();
    assert_eq!(
        ws.multiply_slices_single_into(
            GemmSite::Other,
            &a.as_slice()[..10],
            257,
            3,
            257,
            RealOp::Transpose,
            b.as_slice(),
            2,
            257,
            out.as_mut()
        ),
        Err(LinalgError::Dimension)
    );
    assert_eq!(
        ws.multiply_slices_single_into(
            GemmSite::Other,
            a.as_slice(),
            257,
            3,
            256,
            RealOp::Transpose,
            b.as_slice(),
            2,
            257,
            out.as_mut()
        ),
        Err(LinalgError::Dimension)
    );
    let mut bad = a.clone();
    bad[(0, 0)] = f64::NAN;
    assert_eq!(
        ws.multiply_slices_single_into(
            GemmSite::Other,
            bad.as_slice(),
            257,
            3,
            257,
            RealOp::Transpose,
            b.as_slice(),
            2,
            257,
            out.as_mut()
        ),
        Err(LinalgError::NonFinite)
    );
    // Explicit vendor opt-out remains the same pure result even on MKL builds.
    let pure = RunConfig::from_pairs(&[("mkl_dgemm", "off")]).unwrap();
    let _scope = pure.enter();
    let mut ws = RealGemmWorkspace::with_pool(257, 8, WorkerPool::new(1).unwrap()).unwrap();
    ws.multiply_slices_single_into(
        GemmSite::Other,
        a.as_slice(),
        257,
        3,
        257,
        RealOp::Transpose,
        b.as_slice(),
        2,
        257,
        out.as_mut(),
    )
    .unwrap();
    let first = out.clone();
    ws.multiply_slices_into(
        GemmSite::Other,
        a.as_slice(),
        257,
        3,
        257,
        RealOp::Transpose,
        b.as_slice(),
        2,
        257,
        out.as_mut(),
    )
    .unwrap();
    assert_eq!(first, out);
}
