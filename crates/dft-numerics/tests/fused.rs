use dft_numerics::{
    Fft3d, RustFftPlan3d as FftPlan3d, SpherePattern,
    fused::{FusedOperator, GridFusedPlan, OrderedFusedPlan, PlaneFusedPlan},
};
use ndarray::Array3;
use num_complex::Complex64 as C;
use std::{hint::black_box, time::Instant};

fn random(n: usize, seed: u64) -> Vec<C> {
    let mut state = seed;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 11) as f64 / (1_u64 << 53) as f64 - 0.5
    };
    (0..n).map(|_| C::new(next(), next())).collect()
}
fn sphere(shape: [usize; 3], radius: f64, shift: [f64; 3]) -> Vec<[usize; 3]> {
    let signed = |i: usize, n: usize| {
        if i <= n / 2 {
            i as f64
        } else {
            i as f64 - n as f64
        }
    };
    let mut indices = Vec::new();
    for x in 0..shape[0] {
        for y in 0..shape[1] {
            for z in 0..shape[2] {
                let p = [x, y, z];
                if (0..3)
                    .map(|a| (signed(p[a], shape[a]) + shift[a]).powi(2))
                    .sum::<f64>()
                    <= radius * radius
                {
                    indices.push(p);
                }
            }
        }
    }
    indices.reverse(); // The API must not assume sorted sphere order.
    indices
}
struct Oracle {
    fft: FftPlan3d,
    pattern: SpherePattern,
    indices: Vec<[usize; 3]>,
    reciprocal: Array3<C>,
    real: Array3<C>,
}
impl Oracle {
    fn new(shape: [usize; 3], indices: &[[usize; 3]]) -> Self {
        Self {
            fft: FftPlan3d::new(shape).unwrap(),
            pattern: SpherePattern::from_indices(shape, indices).unwrap(),
            indices: indices.to_vec(),
            reciprocal: Array3::zeros(shape),
            real: Array3::zeros(shape),
        }
    }
    fn scatter(&mut self, c: &[C]) {
        self.reciprocal.fill(C::default());
        for (&i, &v) in self.indices.iter().zip(c) {
            self.reciprocal[i] = v;
        }
    }
    fn local(&mut self, c: &[C], v: &[f64], out: &mut [C]) {
        self.scatter(c);
        self.fft
            .to_real_pruned(&self.reciprocal, &self.pattern, &mut self.real)
            .unwrap();
        for (z, &v) in self.real.iter_mut().zip(v) {
            *z *= v;
        }
        self.fft
            .to_reciprocal_pruned(&self.real, &self.pattern, &mut self.reciprocal)
            .unwrap();
        for (dst, &i) in out.iter_mut().zip(&self.indices) {
            *dst = self.reciprocal[i];
        }
    }
    fn density(&mut self, c: &[C], weight: f64, volume: f64, out: &mut [f64]) {
        self.scatter(c);
        // The pre-change density path uses the full inverse, independently of
        // the Hamiltonian's pruned oracle.
        self.fft.to_real(&self.reciprocal, &mut self.real).unwrap();
        for (dst, z) in out.iter_mut().zip(&self.real) {
            *dst += weight * z.norm_sqr() / volume;
        }
    }
}
fn bits(values: &[C]) -> Vec<(u64, u64)> {
    values
        .iter()
        .map(|z| (z.re.to_bits(), z.im.to_bits()))
        .collect()
}
fn real_bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|z| z.to_bits()).collect()
}
fn compare(shape: [usize; 3], indices: &[[usize; 3]]) -> (f64, f64) {
    let size: usize = shape.iter().product();
    let mut oracle = Oracle::new(shape, indices);
    let mut ordered = OrderedFusedPlan::new(shape, indices).unwrap();
    let mut grid = GridFusedPlan::from_fft(FftPlan3d::new(shape).unwrap(), indices).unwrap();
    let mut plane = PlaneFusedPlan::new(shape, indices).unwrap();
    let mut dense = PlaneFusedPlan::with_pruning(shape, indices, false).unwrap();
    let v: Vec<_> = random(size, 198).iter().map(|z| z.re).collect();
    let mut max_local: f64 = 0.0;
    let mut max_density: f64 = 0.0;
    for seed in [37, 519] {
        let c = random(indices.len(), seed);
        let mut expected = vec![C::default(); c.len()];
        oracle.local(&c, &v, &mut expected);
        let mut actual = expected.clone();
        ordered.apply_local(&c, &v, &mut actual).unwrap();
        assert_eq!(bits(&actual), bits(&expected), "shape={shape:?}");
        grid.apply_local(&c, &v, &mut actual).unwrap();
        assert_eq!(bits(&actual), bits(&expected));
        plane.apply_local(&c, &v, &mut actual).unwrap();
        let scale = expected
            .iter()
            .map(|z| z.norm())
            .fold(0., f64::max)
            .max(1e-300);
        let error = actual
            .iter()
            .zip(&expected)
            .map(|(a, b)| (*a - *b).norm())
            .fold(0., f64::max)
            / scale;
        max_local = max_local.max(error);
        assert!(
            error <= 1e-13,
            "plane relative local error={error:e}, shape={shape:?}"
        );
        let mut control = actual.clone();
        dense.apply_local(&c, &v, &mut control).unwrap();
        assert_eq!(bits(&actual), bits(&control));
        let mut rho = vec![0.; size];
        let mut got = rho.clone();
        let mut reordered = rho.clone();
        for weight in [1.3, 0.25, 0.0] {
            oracle.density(&c, weight, 137.9, &mut rho);
            ordered
                .accumulate_density(&c, weight, 137.9, &mut got)
                .unwrap();
            plane
                .accumulate_density(&c, weight, 137.9, &mut reordered)
                .unwrap();
        }
        assert_eq!(real_bits(&got), real_bits(&rho));
        let scale = rho.iter().copied().fold(0., f64::max).max(1e-300);
        let error = reordered
            .iter()
            .zip(&rho)
            .map(|(a, b)| (a - b).abs())
            .fold(0., f64::max)
            / scale;
        max_density = max_density.max(error);
        assert!(error <= 1e-13, "plane relative density error={error:e}");
    }
    (max_local, max_density)
}
#[test]
fn fused_edge_cases_and_arbitrary_order() {
    for shape in [[1, 1, 1], [1, 3, 5], [5, 1, 3], [3, 5, 1], [7, 9, 11]] {
        for radius in [0., 1., 3., 100.] {
            compare(shape, &sphere(shape, radius, [0.; 3]));
        }
        compare(shape, &[]);
        compare(shape, &[[0, 0, 0], [0, 0, 0]]);
    }
}
#[test]
fn fused_errors_do_not_write_and_zero_weight_skips() {
    for invalid in [[0, 2, 3], [usize::MAX, 2, 3]] {
        assert!(OrderedFusedPlan::new(invalid, &[]).is_err());
        assert!(PlaneFusedPlan::new(invalid, &[]).is_err());
    }
    assert!(OrderedFusedPlan::new([3; 3], &[[3, 0, 0]]).is_err());
    assert!(PlaneFusedPlan::new([3; 3], &[[3, 0, 0]]).is_err());
    for mut plan in [
        Box::new(OrderedFusedPlan::new([3; 3], &[[0; 3]]).unwrap()) as Box<dyn FusedOperator>,
        Box::new(PlaneFusedPlan::new([3; 3], &[[0; 3]]).unwrap()),
    ] {
        let mut out = [C::new(7., -3.)];
        assert!(plan.apply_local(&[], &[1.; 27], &mut out).is_err());
        assert_eq!(out, [C::new(7., -3.)]);
        let mut rho = [-0.; 27];
        assert!(
            plan.accumulate_density(&[C::default()], 1., 0., &mut rho)
                .is_err()
        );
        plan.accumulate_density(&[C::new(f64::NAN, 0.)], 0., 1., &mut rho)
            .unwrap();
        assert_eq!(real_bits(&rho), vec![(-0.0f64).to_bits(); 27]);
    }
}
#[test]
#[ignore = "full shape/support sweep: cargo test --release -p dft-numerics --test fused sweep -- --ignored --nocapture"]
fn sweep() {
    let mut worst = (0.0_f64, 0.0_f64);
    for shape in [
        [24, 24, 24],
        [25, 25, 25],
        [27, 27, 27],
        [36, 36, 36],
        [15, 20, 24],
        [7, 9, 11],
    ] {
        for radius in [0., 1., 3., 6., 10., shape[0] as f64 * 0.7, 100.] {
            for shift in [[0.; 3], [0.25, -0.25, 0.125]] {
                let indices = sphere(shape, radius, shift);
                let e = compare(shape, &indices);
                worst.0 = worst.0.max(e.0);
                worst.1 = worst.1.max(e.1);
            }
        }
    }
    println!(
        "ordered: bit-identical; plane max relative local={:.17e}, density={:.17e}",
        worst.0, worst.1
    );
}
fn samples(mut f: impl FnMut(), calls: usize) -> Vec<f64> {
    f();
    (0..5)
        .map(|_| {
            let t = Instant::now();
            for _ in 0..calls {
                f();
            }
            t.elapsed().as_secs_f64() * 1e6 / calls as f64
        })
        .collect()
}
#[test]
#[ignore = "timing: cargo test --release -p dft-numerics --test fused timing -- --ignored --nocapture"]
fn timing() {
    for shape in [
        [24, 24, 24],
        [25, 25, 25],
        [27, 27, 27],
        [36, 36, 36],
        [15, 20, 24],
    ] {
        let indices = sphere(shape, shape[0] as f64 / 4., [0.25, -0.25, 0.125]);
        let n = shape.iter().product();
        let c = random(indices.len(), 32);
        let v: Vec<_> = random(n, 8).iter().map(|z| z.re).collect();
        let mut out = vec![C::default(); c.len()];
        let mut rho = vec![0.; n];
        let mut oracle = Oracle::new(shape, &indices);
        let mut ordered = OrderedFusedPlan::new(shape, &indices).unwrap();
        let mut plane = PlaneFusedPlan::new(shape, &indices).unwrap();
        println!("shape={shape:?} coefficients={}", c.len());
        println!(
            "old_local_us={:?}",
            samples(
                || {
                    oracle.local(black_box(&c), &v, &mut out);
                    black_box(&out);
                },
                64
            )
        );
        println!(
            "ordered_local_us={:?}",
            samples(
                || {
                    ordered.apply_local(black_box(&c), &v, &mut out).unwrap();
                    black_box(&out);
                },
                64
            )
        );
        println!(
            "plane_local_us={:?}",
            samples(
                || {
                    plane.apply_local(black_box(&c), &v, &mut out).unwrap();
                    black_box(&out);
                },
                64
            )
        );
        println!(
            "old_density_us={:?}",
            samples(
                || {
                    oracle.density(black_box(&c), 1., 100., &mut rho);
                    black_box(&rho);
                },
                64
            )
        );
        println!(
            "ordered_density_us={:?}",
            samples(
                || {
                    ordered
                        .accumulate_density(black_box(&c), 1., 100., &mut rho)
                        .unwrap();
                    black_box(&rho);
                },
                64
            )
        );
        println!(
            "plane_density_us={:?}",
            samples(
                || {
                    plane
                        .accumulate_density(black_box(&c), 1., 100., &mut rho)
                        .unwrap();
                    black_box(&rho);
                },
                64
            )
        );
        // Scheduling groups, with distinct band data; no claimed potential reuse.
        let bands: Vec<_> = (0..32).map(|i| random(c.len(), i + 1)).collect();
        for group in [1, 2, 4, 8, 16, 32] {
            let times = samples(
                || {
                    for chunk in bands.chunks(group) {
                        for c in chunk {
                            plane.apply_local(black_box(c), &v, &mut out).unwrap();
                            black_box(&out);
                        }
                    }
                },
                8,
            )
            .iter()
            .map(|t| t / 32.)
            .collect::<Vec<_>>();
            println!("plane_group={group} us_per_band={times:?}");
        }
        for fraction in [0.05, 0.3, 0.6] {
            let nx = (shape[0] as f64 * fraction).floor().max(1.) as usize;
            let indices: Vec<_> = (0..nx)
                .flat_map(|x| (0..shape[1] / 2).map(move |y| [x, y, 0]))
                .collect();
            let c = random(indices.len(), 32);
            let mut out = vec![C::default(); c.len()];
            for prune in [false, true] {
                let mut p = PlaneFusedPlan::with_pruning(shape, &indices, prune).unwrap();
                println!(
                    "columns={nx}/{} prune={prune} us={:?}",
                    shape[0],
                    samples(
                        || {
                            p.apply_local(&c, &v, &mut out).unwrap();
                            black_box(&out);
                        },
                        64
                    )
                );
            }
        }
    }
}

#[test]
fn fused_charge_and_hermiticity() {
    let shape = [15, 20, 24];
    let indices = sphere(shape, 3.7, [0.25, -0.25, 0.125]);
    let mut a = random(indices.len(), 8);
    let norm = a.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
    for z in &mut a {
        *z /= norm;
    }
    let mut b = random(indices.len(), 14);
    let overlap: C = a.iter().zip(&b).map(|(a, b)| a.conj() * b).sum();
    for (b, a) in b.iter_mut().zip(&a) {
        *b -= a * overlap;
    }
    let norm = b.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
    for z in &mut b {
        *z /= norm;
    }
    let orthogonality: C = a.iter().zip(&b).map(|(a, b)| a.conj() * b).sum();
    let mut plan = OrderedFusedPlan::new(shape, &indices).unwrap();
    let mut rho = vec![0.; shape.iter().product()];
    plan.accumulate_density(&a, 2., 137.9, &mut rho).unwrap();
    let charge_error = (rho.iter().sum::<f64>() * 137.9 / rho.len() as f64 - 2.).abs();
    let v: Vec<_> = random(rho.len(), 19).iter().map(|z| z.re).collect();
    let mut va = a.clone();
    let mut vb = b.clone();
    plan.apply_local(&a, &v, &mut va).unwrap();
    plan.apply_local(&b, &v, &mut vb).unwrap();
    let ab: C = a.iter().zip(&vb).map(|(a, b)| a.conj() * b).sum();
    let ba: C = va.iter().zip(&b).map(|(a, b)| a.conj() * b).sum();
    assert!(charge_error < 1e-12);
    assert!(orthogonality.norm() < 1e-14);
    assert!((ab - ba).norm() < 1e-13);
    println!(
        "fused charge error={charge_error:.17e} electrons; orthogonality={:.17e}; Hermiticity error={:.17e} Ha",
        orthogonality.norm(),
        (ab - ba).norm()
    );
}
