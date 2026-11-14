//! P0.2 work timers (`BAANDI_DAVIDSON_STAGES=on`): timed and untimed local
//! applications are bit-identical, counts are exact and independent of the
//! pool budget, and the run-stage index attributes work. One test function:
//! the counters are process-global, so this binary runs nothing else.
use dft_numerics::{
    PairDensityWorkspace, PairFusedPlan, PairMap, RustFftPlan3d,
    fused::{FusedOperator, GridFusedPlan, OrderedFusedPlan},
    fused_pair::PairLocalWorkspace,
    pool::WorkerPool,
    work_timers::{self, LocalLaps, RunStage},
};
use num_complex::Complex64 as C;

// Sorted ellipsoid of a skewed metric (the construction of real_gamma_pair.rs).
fn fixture(shape: [usize; 3]) -> (Vec<[usize; 3]>, PairMap) {
    let metric = [[1.0, 0.2, 0.0], [0.0, 1.1, 0.3], [0.1, 0.0, 0.9]];
    let mut gs = Vec::new();
    let h = |n: usize| (n as i32 - 1) / 2;
    for x in -h(shape[0])..=h(shape[0]) {
        for y in -h(shape[1])..=h(shape[1]) {
            for z in -h(shape[2])..=h(shape[2]) {
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
fn bits(x: &[f64]) -> Vec<u64> {
    x.iter().map(|v| v.to_bits()).collect()
}
fn cbits(x: &[C]) -> Vec<(u64, u64)> {
    x.iter().map(|v| (v.re.to_bits(), v.im.to_bits())).collect()
}
fn laps_identical(plan: &mut impl FusedOperator, c: &[C], potential: &[f64]) {
    let mut plain = vec![C::default(); c.len()];
    let mut timed = vec![C::default(); c.len()];
    plan.apply_local(c, potential, &mut plain).unwrap();
    let mut laps = LocalLaps::start();
    plan.apply_local_laps(c, potential, &mut timed, &mut laps)
        .unwrap();
    assert_eq!(cbits(&plain), cbits(&timed));
}

#[test]
fn timers_are_exact_counts_and_change_no_bit() {
    let shape = [9, 8, 10];
    let (indices, map) = fixture(shape);
    let n = indices.len();
    let grid: usize = shape.iter().product();
    let potential: Vec<f64> = (0..grid).map(|i| 0.2 * ((i * 7) as f64).sin()).collect();
    let c: Vec<C> = (0..n)
        .map(|i| C::new(((i * 3) as f64).sin(), ((i * 5) as f64).cos()))
        .collect();

    // Plans: the timed body is the untimed body (same operations, same order).
    laps_identical(
        &mut PairFusedPlan::new(shape, &indices).unwrap(),
        &c,
        &potential,
    );
    laps_identical(
        &mut GridFusedPlan::from_fft(RustFftPlan3d::new(shape).unwrap(), &indices).unwrap(),
        &c,
        &potential,
    );
    laps_identical(
        &mut OrderedFusedPlan::new(shape, &indices).unwrap(),
        &c,
        &potential,
    );

    // Block applications: off records nothing; on is bit-identical and counts
    // one pair per two columns (odd tail with a zero partner), at every budget.
    let bands = 5;
    let x: Vec<f64> = (0..n * bands)
        .map(|i| ((i * 17 + 11) as f64).sin() / (n as f64).sqrt())
        .collect();
    work_timers::set_enabled(false);
    let before = work_timers::snapshot();
    let mut reference = vec![0.; x.len()];
    let pool = WorkerPool::new(1).unwrap();
    let mut ws = PairLocalWorkspace::new(map.clone(), bands, pool, || {
        PairFusedPlan::new(shape, &indices)
    })
    .unwrap();
    ws.apply_block(&x, &potential, &mut reference).unwrap();
    assert_eq!(work_timers::snapshot(), before, "off must record nothing");
    work_timers::set_enabled(true);
    for threads in [1, 2, 4] {
        let pool = WorkerPool::new(threads).unwrap();
        let mut ws = PairLocalWorkspace::new(map.clone(), bands, pool, || {
            PairFusedPlan::new(shape, &indices)
        })
        .unwrap();
        let start = work_timers::snapshot()[RunStage::Fine as usize];
        let mut out = vec![0.; x.len()];
        ws.apply_block(&x, &potential, &mut out).unwrap();
        assert_eq!(bits(&out), bits(&reference), "T={threads}");
        let w = work_timers::snapshot()[RunStage::Fine as usize];
        assert_eq!(w.pair_transforms - start.pair_transforms, 3);
        assert_eq!(w.pair_columns - start.pair_columns, bands as u64);
        assert!(w.pair_inverse > start.pair_inverse && w.pair_forward > start.pair_forward);
        assert_eq!(
            w.pair_unsplit, start.pair_unsplit,
            "mirror-slab plan splits"
        );
    }

    // Density pieces count one inverse transform per pair, in the run stage
    // entered by the caller; the guard restores the previous stage.
    let pool = WorkerPool::new(2).unwrap();
    let mut density_ws =
        PairDensityWorkspace::new(map, bands, pool, || PairFusedPlan::new(shape, &indices))
            .unwrap();
    let mut density = vec![0.; grid];
    let weights = [2., 2., 2., 0., 1.];
    let coarse_before = work_timers::snapshot()[RunStage::Coarse as usize];
    {
        let _stage = work_timers::enter(RunStage::Coarse);
        assert!(!work_timers::in_fine_stage());
        density_ws
            .density(&x, &weights, 10.0, &mut density)
            .unwrap();
    }
    assert!(work_timers::in_fine_stage());
    let coarse = work_timers::snapshot()[RunStage::Coarse as usize];
    // Pairs (0,1), (2,3) and the tail (4): three inverse transforms, five bands.
    assert_eq!(
        coarse.density_transforms - coarse_before.density_transforms,
        3
    );
    assert_eq!(coarse.density_columns - coarse_before.density_columns, 5);
    assert_eq!(coarse.pair_transforms, coarse_before.pair_transforms);
    work_timers::set_enabled(false);
}
