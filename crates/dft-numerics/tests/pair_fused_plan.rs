use dft_numerics::{
    PairDensityOperator, PairDensityWorkspace, PairFusedPlan, PairMap, PairPlan,
    fused::{FusedOperator, OrderedFusedPlan},
    fused_pair::PairLocalWorkspace,
    pool::WorkerPool,
};
use num_complex::Complex64 as C;
// Synthetic reciprocal metrics for cubic, skewed, elongated and slab cells.
// Enumerate a physical ellipsoid, sorted by |G|^2 then Miller index, just as
// the basis layer does. No reverse dependency on dft-basis is introduced.
fn fixture(shape: [usize; 3], metric: [[f64; 3]; 3]) -> (Vec<[usize; 3]>, PairMap) {
    let mut gs = Vec::new();
    for x in -((shape[0] as i32 - 1) / 2)..=((shape[0] as i32 - 1) / 2) {
        for y in -((shape[1] as i32 - 1) / 2)..=((shape[1] as i32 - 1) / 2) {
            for z in -((shape[2] as i32 - 1) / 2)..=((shape[2] as i32 - 1) / 2) {
                let m = [x, y, z];
                let q: f64 = metric
                    .iter()
                    .map(|r| {
                        r.iter()
                            .zip(m)
                            .map(|(a, b)| a * b as f64)
                            .sum::<f64>()
                            .powi(2)
                    })
                    .sum();
                if q <= 9. {
                    gs.push((q, m));
                }
            }
        }
    }
    gs.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
    let indices: Vec<_> = gs
        .iter()
        .map(|(_, m)| std::array::from_fn(|d| m[d].rem_euclid(shape[d] as i32) as usize))
        .collect();
    let mut plus = vec![0];
    let mut minus = vec![0];
    for (i, (_, m)) in gs.iter().enumerate().skip(1) {
        if *m > [0; 3] {
            plus.push(i as u32);
            minus.push(
                gs.iter()
                    .position(|(_, other)| *other == m.map(|v| -v))
                    .unwrap() as u32,
            );
        }
    }
    let map = PairMap::new(shape, &indices, &plus, &minus).unwrap();
    (indices, map)
}
fn states(n: usize, bands: usize) -> Vec<f64> {
    let mut out: Vec<_> = (0..n * bands)
        .map(|i| ((i * 17 + 11) as f64).sin())
        .collect();
    for col in out.chunks_mut(n) {
        let norm = col.iter().map(|x| x * x).sum::<f64>().sqrt();
        for x in col {
            *x /= norm;
        }
    }
    out
}

fn relative(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt()
        / b.iter()
            .map(|v| v * v)
            .sum::<f64>()
            .sqrt()
            .max(f64::MIN_POSITIVE)
}
fn bits(v: &[f64]) -> Vec<u64> {
    v.iter().map(|v| v.to_bits()).collect()
}

#[test]
fn boundaries_and_validation() {
    assert!(PairFusedPlan::new([0, 1, 1], &[]).is_err());
    assert!(PairFusedPlan::new([usize::MAX, 2, 1], &[]).is_err());
    assert!(PairFusedPlan::new([2; 3], &[[2, 0, 0]]).is_err());
    for indices in [
        vec![],
        vec![[0, 0, 0]],
        vec![[1, 0, 0], [0, 0, 0], [1, 0, 0]],
    ] {
        let mut p = PairFusedPlan::new([3, 1, 1], &indices).unwrap();
        let mut q = OrderedFusedPlan::new([3, 1, 1], &indices).unwrap();
        let c = vec![C::new(0.3, 0.7); indices.len()];
        let mut a = c.clone();
        let mut b = c.clone();
        p.apply_local(&c, &[0.4; 3], &mut a).unwrap();
        q.apply_local(&c, &[0.4; 3], &mut b).unwrap();
        for (a, b) in a.iter().zip(b) {
            assert!((*a - b).norm() < 1e-14);
        }
        let saved = a.clone();
        assert!(p.apply_local(&c, &[], &mut a).is_err());
        assert_eq!(a, saved);
        let mut rho = [3.; 3];
        p.accumulate_density_pair(&c, [0.; 2], 10., &mut rho)
            .unwrap();
        assert_eq!(rho, [3.; 3]);
        assert!(
            p.accumulate_density_pair(&c, [f64::NAN, 1.], 10., &mut rho)
                .is_err()
        );
        assert!(p.accumulate_density(&c, 1., 0., &mut rho).is_err());
        assert_eq!(rho, [3.; 3]);
    }
}

#[test]
#[ignore = "release matrix: cargo test --release -p dft-numerics --test pair_fused_plan -- --include-ignored --nocapture --test-threads=1"]
fn nine_shapes_and_worker_determinism() {
    for (case, shape) in [
        [20, 20, 20],
        [25, 25, 25],
        [36, 36, 36],
        [48, 48, 48],
        [72, 36, 36],
        [36, 36, 18],
        [36, 36, 72],
        [15, 20, 24],
        [36, 36, 36],
        [20, 40, 24],
    ]
    .into_iter()
    .enumerate()
    {
        let sc = shape.map(|n| 12. / n as f64);
        let mut metric = [[sc[0], 0., 0.], [0., sc[1], 0.], [0., 0., sc[2]]];
        if case == 6 {
            metric[2][2] = 0.5;
        }
        if case == 8 {
            metric[0][1] = 0.35 * sc[1];
        }
        let (indices, map) = fixture(shape, metric);
        let n = indices.len();
        let size = shape.iter().product();
        let v: Vec<_> = (0..size)
            .map(|i| 0.7 + 0.2 * (i as f64 * 0.19).cos())
            .collect();
        let x = states(n, 9);
        let mut phi = vec![C::default(); n];
        map.pack(&x[..n], Some(&x[n..2 * n]), &mut phi).unwrap();
        let mut p = PairFusedPlan::new(shape, &indices).unwrap();
        let mut q = OrderedFusedPlan::new(shape, &indices).unwrap();
        let mut a = phi.clone();
        let mut b = phi.clone();
        p.apply_local(&phi, &v, &mut a).unwrap();
        q.apply_local(&phi, &v, &mut b).unwrap();
        let err = a
            .iter()
            .zip(&b)
            .map(|(a, b)| (*a - *b).norm_sqr())
            .sum::<f64>()
            .sqrt()
            / b.iter().map(|v| v.norm_sqr()).sum::<f64>().sqrt();
        assert!(err < 1e-13, "{shape:?} {err}");
        let mut da = vec![0.; size];
        let mut db = da.clone();
        p.accumulate_density_pair(&phi, [0.3, 0.7], 123., &mut da)
            .unwrap();
        q.accumulate_density_pair(&phi, [0.3, 0.7], 123., &mut db)
            .unwrap();
        let de = relative(&da, &db);
        assert!(de < 1e-13);
        da.fill(0.);
        db.fill(0.);
        p.accumulate_density(&phi, 0.6, 123., &mut da).unwrap();
        q.accumulate_density(&phi, 0.6, 123., &mut db).unwrap();
        assert!(relative(&da, &db) < 1e-13);
        let mut pair =
            PairPlan::new(PairFusedPlan::new(shape, &indices).unwrap(), map.clone()).unwrap();
        let mut reference = vec![0.; x.len()];
        pair.apply_block(&x, &v, &mut reference).unwrap();
        let mut rho_reference = None;
        for threads in [1, 2, 4, 8, 12] {
            let pool = WorkerPool::new(threads).unwrap();
            let mut local = PairLocalWorkspace::new(map.clone(), 9, pool.clone(), || {
                PairFusedPlan::new(shape, &indices)
            })
            .unwrap();
            let mut density = PairDensityWorkspace::new(map.clone(), 9, pool, || {
                PairFusedPlan::new(shape, &indices)
            })
            .unwrap();
            let mut out = vec![0.; x.len()];
            let mut rho = vec![0.; size];
            for _ in 0..2 {
                local.apply_block(&x, &v, &mut out).unwrap();
                assert_eq!(bits(&out), bits(&reference));
                density.density(&x, &[0.4; 9], 123., &mut rho).unwrap();
                if let Some(ref r) = rho_reference {
                    assert_eq!(&bits(&rho), r);
                } else {
                    rho_reference = Some(bits(&rho));
                }
            }
        }
        println!(
            "fixture={case} shape={shape:?} npw={n} local_relative={err:.17e} density_relative={de:.17e}"
        );
    }
}

// Test-only System allocator adapter for the requested no-allocation gate.
// Run alone with one test thread so unrelated tests cannot enter the count.
#[allow(unsafe_code)]
mod allocation_probe {
    use std::{
        alloc::{GlobalAlloc, Layout, System},
        sync::atomic::{AtomicBool, AtomicUsize, Ordering},
    };
    pub static TRACK: AtomicBool = AtomicBool::new(false);
    pub static COUNT: AtomicUsize = AtomicUsize::new(0);
    pub struct Allocator;
    fn record() {
        if TRACK.load(Ordering::Relaxed) {
            COUNT.fetch_add(1, Ordering::Relaxed);
        }
    }
    // SAFETY: all pointers/layouts are forwarded unchanged to the System allocator;
    // counters use atomics and never allocate or alter allocation semantics.
    unsafe impl GlobalAlloc for Allocator {
        unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
            record();
            // SAFETY: caller supplies a valid allocation layout.
            unsafe { System.alloc(layout) }
        }
        unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
            record();
            // SAFETY: caller supplies a valid allocation layout.
            unsafe { System.alloc_zeroed(layout) }
        }
        unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
            // SAFETY: pointer and original layout come from this forwarding allocator.
            unsafe { System.dealloc(ptr, layout) }
        }
        unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, size: usize) -> *mut u8 {
            record();
            // SAFETY: pointer/layout/new size retain the GlobalAlloc caller contract.
            unsafe { System.realloc(ptr, layout, size) }
        }
    }
}
#[global_allocator]
static ALLOCATOR: allocation_probe::Allocator = allocation_probe::Allocator;

#[test]
#[ignore = "release allocation probe; run this integration test with --test-threads=1"]
fn no_allocation_after_construction() {
    use std::sync::atomic::Ordering;
    for shape in [
        [25; 3],
        [36; 3],
        [48; 3],
        [72, 36, 36],
        [15, 20, 24],
        [20, 40, 24],
    ] {
        let sc = shape.map(|d| 12. / d as f64);
        let (indices, map) = fixture(shape, [[sc[0], 0., 0.], [0., sc[1], 0.], [0., 0., sc[2]]]);
        let n = indices.len();
        let size = shape.iter().product();
        let x = states(n, 2);
        let mut p = PairPlan::new(PairFusedPlan::new(shape, &indices).unwrap(), map).unwrap();
        let mut out = vec![0.; x.len()];
        let v = vec![0.7; size];
        let mut rho = vec![0.; size];
        allocation_probe::COUNT.store(0, Ordering::SeqCst);
        allocation_probe::TRACK.store(true, Ordering::SeqCst);
        for _ in 0..3 {
            p.apply_block(&x, &v, &mut out).unwrap();
            p.accumulate_density_pair(&x[..n], Some(&x[n..]), [0.3, 0.7], 123., &mut rho)
                .unwrap();
        }
        allocation_probe::TRACK.store(false, Ordering::SeqCst);
        let count = allocation_probe::COUNT.load(Ordering::SeqCst);
        println!("allocation shape={shape:?} count={count}");
        assert_eq!(count, 0);
    }
}
