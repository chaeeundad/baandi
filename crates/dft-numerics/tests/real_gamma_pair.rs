use dft_numerics::{
    PairDensityOperator, PairDensityWorkspace, PairMap, PairPlan, RustFftPlan3d,
    fused::{FusedOperator, GridFusedPlan, OrderedFusedPlan},
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
fn completed(map: &PairMap, x: &[f64]) -> Vec<C> {
    let mut out = vec![C::default(); map.len()];
    map.pack(x, None, &mut out).unwrap();
    out
}
fn norm_error(a: &[f64], b: &[f64]) -> f64 {
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
fn check<P: PairDensityOperator + Send>(
    shape: [usize; 3],
    indices: &[[usize; 3]],
    map: PairMap,
    make: impl Fn() -> P,
) {
    let n = map.len();
    let grid: usize = shape.iter().product();
    let bands = 19;
    let x = states(n, bands);
    let v: Vec<_> = (0..grid)
        .map(|i| 0.7 + 0.2 * (i as f64 * 0.19).cos())
        .collect();
    let mut pair = PairPlan::new(make(), map.clone()).unwrap();
    let mut got = vec![0.; x.len()];
    pair.apply_block(&x, &v, &mut got).unwrap();
    // R7: old serial method is the byte-for-byte oracle. Reuse the same
    // workspace at changing active widths, including an empty block and tails.
    for threads in [1, 2, 4, 8, 12] {
        let pool = WorkerPool::new(threads).unwrap();
        let mut parallel = dft_numerics::fused_pair::PairLocalWorkspace::new(
            map.clone(),
            bands,
            pool.clone(),
            || Ok(make()),
        )
        .unwrap();
        for width in [19, 0, 1, 2, 3, 8, 18, 19] {
            let input = &x[..width * n];
            let mut serial = vec![0.; input.len()];
            let mut output = vec![0.; input.len()];
            pair.apply_block(input, &v, &mut serial).unwrap();
            pool.install(|| parallel.apply_block(input, &v, &mut output))
                .unwrap();
            assert!(
                serial
                    .iter()
                    .zip(&output)
                    .all(|(a, b)| a.to_bits() == b.to_bits())
            );
        }
        assert!(
            parallel
                .apply_block(&x, &v[..v.len() - 1], &mut got)
                .is_err()
        );
        assert!(
            parallel
                .apply_block(&x[..x.len() - 1], &v, &mut got)
                .is_err()
        );
        assert_eq!(
            pool.workers_created(),
            if threads == 1 { 0 } else { threads }
        );
    }
    pair.apply_block(&x, &v, &mut got).unwrap();
    let mut single = make();
    let mut expected = vec![0.; x.len()];
    let mut ref_density = vec![0.; grid];
    let weights: Vec<_> = (0..bands)
        .map(|i| {
            if (4..8).contains(&i) {
                0.
            } else {
                (i % 3) as f64 * 0.35
            }
        })
        .collect();
    for (band, col) in x.chunks(n).enumerate() {
        let full = completed(&map, col);
        let mut h = full.clone();
        single.apply_local(&full, &v, &mut h).unwrap();
        map.unpack(&h, &mut expected[band * n..(band + 1) * n], None)
            .unwrap();
        // Compare completed pair output to the full complex result as well,
        // so the parity projection cannot hide imaginary/G=0 defects.
        let recovered = completed(&map, &got[band * n..(band + 1) * n]);
        let error = recovered
            .iter()
            .zip(&h)
            .map(|(a, b)| (a - b).norm_sqr())
            .sum::<f64>()
            .sqrt()
            / h.iter().map(C::norm_sqr).sum::<f64>().sqrt();
        assert!(error < 1e-14, "full-complex application {error:e}");
        single
            .accumulate_density(&full, weights[band], 120., &mut ref_density)
            .unwrap();
    }
    // 1e-14 relative admits FFT rounding of the combined unit-norm functions,
    // ORC sqrt(2) factors and the different fixed density grouping, not physics error.
    let apply_error = norm_error(&got, &expected);
    assert!(apply_error < 1e-14);
    let mut bits = None;
    let mut density_error = 0.0_f64;
    for threads in [1, 2, 4, 8] {
        let pool = WorkerPool::new(threads).unwrap();
        let mut ws =
            PairDensityWorkspace::new(map.clone(), bands, pool.clone(), || Ok(make())).unwrap();
        let mut density = vec![0.; grid];
        pool.install(|| ws.density(&x, &weights, 120., &mut density))
            .unwrap();
        density_error = density_error.max(norm_error(&density, &ref_density));
        assert!(density_error < 1e-14, "density {density_error:e}");
        let current: Vec<_> = density.iter().map(|x| x.to_bits()).collect();
        if let Some(ref bits) = bits {
            assert_eq!(bits, &current);
        } else {
            bits = Some(current);
        }
        ws.density(&x, &weights, 120., &mut density).unwrap();
        assert_eq!(
            bits.as_ref().unwrap(),
            &density.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
        );
    }
    // Exactly two zero weights skip even nonfinite coefficient contents.
    let nan = vec![f64::NAN; n];
    let mut untouched = vec![-0.; grid];
    pair.accumulate_density_pair(&nan, Some(&nan), [0., -0.], 120., &mut untouched)
        .unwrap();
    assert!(
        untouched
            .iter()
            .all(|x| x.to_bits() == (-0.0_f64).to_bits())
    );
    assert!(pair.apply_block(&x[..x.len() - 1], &v, &mut got).is_err());
    assert!(
        pair.accumulate_density_pair(&x[..n], None, [1., 0.], 0., &mut untouched)
            .is_err()
    );
    println!(
        "shape={shape:?} n={} apply_relative={apply_error:e} density_relative={density_error:e} budgets=1,2,4,8 identical",
        indices.len()
    );
}
#[test]
fn pairs_four_shapes_both_backends() {
    for (shape, metric) in [
        ([11, 11, 11], [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]),
        (
            [11, 12, 13],
            [[1., 0.2, 0.], [0., 0.9, 0.2], [0.1, 0., 1.1]],
        ),
        ([11, 12, 19], [[1.4, 0., 0.], [0., 1.2, 0.], [0., 0., 0.6]]),
        ([13, 14, 21], [[1., 0., 0.], [0., 1., 0.], [0., 0., 0.4]]),
    ] {
        let (indices, map) = fixture(shape, metric);
        check(shape, &indices, map.clone(), || {
            OrderedFusedPlan::new(shape, &indices).unwrap()
        });
        check(shape, &indices, map, || {
            GridFusedPlan::from_fft(RustFftPlan3d::new(shape).unwrap(), &indices).unwrap()
        });
    }
}
#[test]
fn parity_zero_and_validation() {
    let indices = [[0, 0, 0], [1, 0, 0], [4, 0, 0]];
    let map = PairMap::new([5, 1, 1], &indices, &[0, 1], &[0, 2]).unwrap();
    let a = [0.5, 0.3, -0.2];
    let b = [-0.4, 0.1, 0.7];
    let mut phi = vec![C::default(); 3];
    map.pack(&a, Some(&b), &mut phi).unwrap();
    let mut ra = [0.; 3];
    let mut rb = [0.; 3];
    map.unpack(&phi, &mut ra, Some(&mut rb)).unwrap();
    assert!(norm_error(&ra, &a) < 1e-14);
    assert!(norm_error(&rb, &b) < 1e-14);
    phi[0] = C::new(2., 1e-12);
    map.unpack(&phi, &mut ra, None).unwrap();
    assert_eq!(ra[0], 2.);
    assert_eq!(completed(&map, &ra)[0].im.to_bits(), 0.0_f64.to_bits());
    map.unpack(&phi, &mut ra, Some(&mut rb)).unwrap();
    assert_eq!(rb[0], 1e-12);
    assert_eq!(completed(&map, &rb)[0].im.to_bits(), 0.0_f64.to_bits());
    assert!(PairMap::new([5, 1, 1], &indices, &[0, 1], &[0, 1]).is_err());
    assert!(PairMap::new([6, 1, 1], &indices, &[0, 1], &[0, 2]).is_err());
    assert!(PairMap::new([5, 1, 1], &[], &[], &[]).is_err());
    let zero_map = PairMap::new([1; 3], &[[0; 3]], &[0], &[0]).unwrap();
    let mut plan =
        PairPlan::new(OrderedFusedPlan::new([1; 3], &[[0; 3]]).unwrap(), zero_map).unwrap();
    let mut out = [0.; 3];
    plan.apply_block(&[1., 2., 3.], &[4.], &mut out).unwrap();
    assert_eq!(out, [4., 8., 12.]);
}
#[test]
#[cfg(debug_assertions)]
#[should_panic(expected = "within a factor 1e3")]
fn rejects_unbalanced_nonzero_pair() {
    let map = PairMap::new([1; 3], &[[0; 3]], &[0], &[0]).unwrap();
    map.pack(&[1.], Some(&[1e-4]), &mut [C::default()]).unwrap();
}

/// Observation only (the machine that runs it; the printed label is `observation`); release command documented in the R2 report.
/// Grids include a non-cubic 72x36x36 (the si16 grid) next to the cubic ones.
#[test]
#[ignore]
fn pair_microbenchmark() {
    use std::{hint::black_box, time::Instant};
    for shape in [[25, 25, 25], [36, 36, 36], [48, 48, 48], [72, 36, 36]] {
        let points = shape[0] * shape[1] * shape[2];
        let sc = shape.map(|d| 12. / d as f64);
        let (indices, map) = fixture(shape, [[sc[0], 0., 0.], [0., sc[1], 0.], [0., 0., sc[2]]]);
        let n = map.len();
        let x = states(n, 2);
        let a = completed(&map, &x[..n]);
        let b = completed(&map, &x[n..]);
        let potential: Vec<_> = (0..points)
            .map(|i| 0.7 + 0.2 * (i as f64 * 0.19).cos())
            .collect();
        let mut pair =
            PairPlan::new(OrderedFusedPlan::new(shape, &indices).unwrap(), map.clone()).unwrap();
        let mut mirror = PairPlan::new(
            dft_numerics::PairFusedPlan::new(shape, &indices).unwrap(),
            map,
        )
        .unwrap();
        let mut single = OrderedFusedPlan::new(shape, &indices).unwrap();
        let mut out = vec![0.; 2 * n];
        let mut oa = vec![C::default(); n];
        let mut ob = oa.clone();
        let mut samples = [Vec::new(), Vec::new(), Vec::new()];
        for round in 0..32 {
            for mode in if round % 2 == 0 { [0, 1, 2] } else { [2, 1, 0] } {
                let start = Instant::now();
                for _ in 0..20 {
                    if mode == 0 {
                        pair.apply_block(black_box(&x), black_box(&potential), &mut out)
                            .unwrap();
                        black_box(&out);
                    } else if mode == 2 {
                        mirror
                            .apply_block(black_box(&x), black_box(&potential), &mut out)
                            .unwrap();
                        black_box(&out);
                    } else {
                        single
                            .apply_local(black_box(&a), black_box(&potential), &mut oa)
                            .unwrap();
                        single
                            .apply_local(black_box(&b), black_box(&potential), &mut ob)
                            .unwrap();
                        black_box((&oa, &ob));
                    }
                }
                if round >= 2 {
                    samples[mode].push(start.elapsed().as_secs_f64() * 1e6 / 40.);
                }
            }
        }
        println!("raw grid={shape:?} library_single_mirror_us_per_band={samples:?}");
        for s in &mut samples {
            s.sort_by(f64::total_cmp);
        }
        let med = samples.each_ref().map(|s| (s[14] + s[15]) * 0.5);
        println!(
            "observation grid={}x{}x{} n={n} samples=30 pair_us_per_band={:.6} single_us_per_band={:.6} ratio={:.6}",
            shape[0],
            shape[1],
            shape[2],
            med[0],
            med[1],
            med[0] / med[1]
        );
        println!(
            "mirror_observation grid={shape:?} samples=30 library_us_per_band={:.6} mirror_us_per_band={:.6} mirror_over_library={:.6}",
            med[0],
            med[2],
            med[2] / med[0]
        );
    }
}

#[test]
fn fallback_transform_counts_and_zero_weights() {
    use dft_numerics::{Fft3d, FftError};
    use ndarray::Array3;
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };
    struct Counted {
        fft: RustFftPlan3d,
        calls: Arc<[AtomicUsize; 2]>,
    }
    impl Fft3d for Counted {
        fn shape(&self) -> [usize; 3] {
            self.fft.shape()
        }
        fn to_real(&mut self, a: &Array3<C>, b: &mut Array3<C>) -> Result<(), FftError> {
            self.calls[0].fetch_add(1, Ordering::Relaxed);
            self.fft.to_real(a, b)
        }
        fn to_reciprocal(&mut self, a: &Array3<C>, b: &mut Array3<C>) -> Result<(), FftError> {
            self.calls[1].fetch_add(1, Ordering::Relaxed);
            self.fft.to_reciprocal(a, b)
        }
    }
    let calls = Arc::new([AtomicUsize::new(0), AtomicUsize::new(0)]);
    let indices = [[0; 3]];
    let map = PairMap::new([1; 3], &indices, &[0], &[0]).unwrap();
    let backend = GridFusedPlan::from_fft(
        Counted {
            fft: RustFftPlan3d::new([1; 3]).unwrap(),
            calls: calls.clone(),
        },
        &indices,
    )
    .unwrap();
    let mut pair = PairPlan::new(backend, map).unwrap();
    let mut out = [0.; 3];
    pair.apply_block(&[1., 2., 3.], &[2.], &mut out).unwrap();
    assert_eq!(out, [2., 4., 6.]);
    assert_eq!(calls.each_ref().map(|n| n.load(Ordering::Relaxed)), [2, 2]);
    let mut density = [0.];
    for weights in [[0., 0.], [0., 0.5], [0.25, 0.], [0., 1e-300]] {
        pair.accumulate_density_pair(&[1.], Some(&[2.]), weights, 2., &mut density)
            .unwrap();
    }
    assert_eq!(calls.each_ref().map(|n| n.load(Ordering::Relaxed)), [5, 2]);
    assert_eq!(density, [1.125]);
    let before = density;
    assert!(
        pair.accumulate_density_pair(&[1.], Some(&[2.]), [f64::NAN, 0.], 2., &mut density)
            .is_err()
    );
    assert_eq!(density, before);
    assert_eq!(calls.each_ref().map(|n| n.load(Ordering::Relaxed)), [5, 2]);
}

#[test]
fn empty_density_block() {
    let map = PairMap::new([1; 3], &[[0; 3]], &[0], &[0]).unwrap();
    let mut ws = PairDensityWorkspace::new(map, 0, WorkerPool::new(1).unwrap(), || {
        OrderedFusedPlan::new([1; 3], &[[0; 3]])
    })
    .unwrap();
    let mut out = [9.];
    ws.density(&[], &[], 1., &mut out).unwrap();
    assert_eq!(out, [0.]);
}
