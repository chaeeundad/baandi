//! Run on dev08; timings include the public adapter copies, never just DFTI.
use dft_mkl::Plan3d;
use dft_numerics::{Fft3d, RustFftPlan3d, SpherePattern};
use ndarray::Array3;
use num_complex::Complex64 as C;
use std::time::Instant;

fn random(n: usize) -> Vec<C> {
    let mut state = 12345u64;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 11) as f64 / (1u64 << 53) as f64 - 0.5
    };
    (0..n).map(|_| C::new(next(), next())).collect()
}
fn relative(a: &[C], b: &[C]) -> f64 {
    assert!(
        a.iter()
            .chain(b)
            .all(|z| z.re.is_finite() && z.im.is_finite())
    );
    a.iter()
        .zip(b)
        .map(|(a, b)| {
            let err = (*a - *b).norm();
            if err == 0.0 { 0.0 } else { err / b.norm() }
        })
        .fold(0.0, f64::max)
}
fn timed(mut f: impl FnMut(), reps: usize) -> f64 {
    for _ in 0..3 {
        f();
    }
    let start = Instant::now();
    for _ in 0..reps {
        f();
    }
    start.elapsed().as_secs_f64() * 1e6 / reps as f64
}
// Independent direct DFT with compensated complex summation. No FFT routines.
// f64 sin/cos phase evaluation remains a limitation; this is not an arbitrary-
// precision oracle. Small cases below test ALL coefficients; large ones sample.
fn reference(input: &[C], shape: [usize; 3], q: [usize; 3]) -> C {
    let mut sum = C::default();
    let mut compensation = C::default();
    for (p, z) in input.iter().enumerate() {
        let r = [
            p / (shape[1] * shape[2]),
            p / shape[2] % shape[1],
            p % shape[2],
        ];
        let phase = -std::f64::consts::TAU
            * (0..3)
                .map(|a| (q[a] * r[a]) as f64 / shape[a] as f64)
                .sum::<f64>();
        let term = *z * C::from_polar(1.0, phase) - compensation;
        let next = sum + term;
        compensation = (next - sum) - term;
        sum = next;
    }
    sum / input.len() as f64
}
fn small_reference() {
    for shape in [[1, 2, 9], [2, 3, 5], [4, 5, 6], [7, 4, 3]] {
        let input = random(shape.iter().product());
        let mut output = input.clone();
        Plan3d::new(shape)
            .unwrap()
            .transform(&input, &mut output, false)
            .unwrap();
        assert!(output.iter().all(|z| z.re.is_finite() && z.im.is_finite()));
        let mut max: f64 = 0.0;
        for (p, z) in output.iter().enumerate() {
            let q = [
                p / (shape[1] * shape[2]),
                p / shape[2] % shape[1],
                p % shape[2],
            ];
            max = max.max((*z - reference(&input, shape, q)).norm());
        }
        println!("reference_small,{shape:?},max_abs={max:e}");
        assert!(max < 1e-13);
    }
}
fn main() {
    small_reference();
    let max_mib: usize = std::env::var("MKL_BENCH_MAX_MIB")
        .unwrap_or("1024".into())
        .parse()
        .unwrap();
    let reps: usize = std::env::var("MKL_BENCH_REPS")
        .unwrap_or("10".into())
        .parse()
        .unwrap();
    assert!(reps > 0);
    println!("backend,edge,plan_us,forward_us,inverse_us");
    let mut passed = true;
    for n in [24usize, 25, 27, 32, 36, 48, 64, 96, 128, 192] {
        let len = n * n * n;
        // Conservative user budget estimate, including patterns and snapshots;
        // MKL internal workspace is opaque. Set the limit below available RAM.
        if len.saturating_mul(256) > max_mib.saturating_mul(1024 * 1024) {
            println!("skip,{n},estimated_bytes={}", len * 256);
            continue;
        }
        let data = random(len);
        let a = Array3::from_shape_vec((n, n, n), data.clone()).unwrap();
        let mut b = a.clone();
        let t = Instant::now();
        let mut rust = RustFftPlan3d::new([n; 3]).unwrap();
        let rust_plan = t.elapsed().as_secs_f64() * 1e6;
        let rf = timed(|| rust.to_reciprocal(&a, &mut b).unwrap(), reps);
        let expected = b.as_slice().unwrap().to_vec();
        let ri = timed(|| rust.to_real(&a, &mut b).unwrap(), reps);
        println!("rustfft,{n},{rust_plan},{rf},{ri}");
        let t = Instant::now();
        let mut mkl = Plan3d::new([n; 3]).unwrap();
        let mkl_plan = t.elapsed().as_secs_f64() * 1e6;
        let mut out = vec![C::default(); len];
        let mf = timed(|| mkl.transform(&data, &mut out, false).unwrap(), reps);
        let cross = relative(&out, &expected);
        let norm_in: f64 = data.iter().map(|z| z.norm_sqr()).sum();
        let norm_out: f64 = out.iter().map(|z| z.norm_sqr()).sum();
        let parseval = (norm_out * len as f64 / norm_in - 1.0).abs();
        let mut back = vec![C::default(); len];
        mkl.transform(&out, &mut back, true).unwrap();
        assert!(back.iter().all(|z| z.re.is_finite() && z.im.is_finite()));
        let roundtrip = data
            .iter()
            .zip(&back)
            .map(|(a, b)| (*a - *b).norm())
            .fold(0.0, f64::max);
        let mut reference_error: f64 = 0.0;
        for q in [[0, 0, 0], [1, 2, 3], [n / 2, n / 3, n / 4], [n - 1, 0, 1]] {
            reference_error = reference_error
                .max((out[(q[0] * n + q[1]) * n + q[2]] - reference(&data, [n; 3], q)).norm());
        }
        // Each worker owns and creates its plan. Compare the same inputs bitwise,
        // across 1/2/4 outer workers; no thread-global MKL settings are changed.
        for threads in [1, 2, 4] {
            std::thread::scope(|scope| {
                let jobs: Vec<_> = (0..threads)
                    .map(|_| {
                        let data = &data;
                        scope.spawn(move || {
                            let mut p = Plan3d::new([n; 3]).unwrap();
                            let mut result = vec![C::default(); len];
                            p.transform(data, &mut result, false).unwrap();
                            result
                        })
                    })
                    .collect();
                for job in jobs {
                    let actual = job.join().unwrap();
                    assert!(
                        actual
                            .iter()
                            .zip(&out)
                            .all(|(a, b)| a.re.to_bits() == b.re.to_bits()
                                && a.im.to_bits() == b.im.to_bits())
                    );
                }
            });
        }
        let mi = timed(|| mkl.transform(&data, &mut out, true).unwrap(), reps);
        println!("mkl,{n},{mkl_plan},{mf},{mi}");
        println!(
            "accuracy,{n},cross_max_relative={cross:e},roundtrip_max_abs={roundtrip:e},parseval_relative={parseval:e},sample_kahan_max_abs={reference_error:e},threads_1_2_4=identical"
        );
        passed &=
            cross <= 1e-13 && roundtrip <= 1e-13 && parseval <= 1e-13 && reference_error <= 1e-13;
        // Sphere occupies ~6.5% of coefficients; same supported input for both.
        let indices: Vec<_> = (0..len)
            .filter_map(|p| {
                let q = [p / (n * n), p / n % n, p % n];
                let radius2: usize = q.iter().map(|&x| x.min(n - x).pow(2)).sum();
                (radius2 <= (n / 4).pow(2)).then_some(q)
            })
            .collect();
        let t = Instant::now();
        let pattern = SpherePattern::from_indices([n; 3], &indices).unwrap();
        println!(
            "pattern,{n},creation_us={}",
            t.elapsed().as_secs_f64() * 1e6
        );
        let mut supported = Array3::<C>::zeros((n, n, n));
        for q in &indices {
            supported[*q] = a[*q];
        }
        let pf = timed(
            || rust.to_reciprocal_pruned(&a, &pattern, &mut b).unwrap(),
            reps,
        );
        let pi = timed(
            || rust.to_real_pruned(&supported, &pattern, &mut b).unwrap(),
            reps,
        );
        println!("rustfft_pruned,{n},{rust_plan},{pf},{pi}");
        let mut present = vec![false; len];
        for q in &indices {
            present[(q[0] * n + q[1]) * n + q[2]] = true;
        }
        let pmf = timed(
            || {
                mkl.transform(&data, &mut out, false).unwrap();
                for (z, p) in out.iter_mut().zip(&present) {
                    if !p {
                        *z = C::default();
                    }
                }
            },
            reps,
        );
        let pmi = timed(
            || {
                mkl.transform(supported.as_slice().unwrap(), &mut out, true)
                    .unwrap()
            },
            reps,
        );
        println!("mkl_full_mask,{n},{mkl_plan},{pmf},{pmi}");
    }
    assert!(
        passed,
        "strict FFT accuracy gate failed; retain measured errors and investigate, do not relax tolerances"
    );
    println!("PASS FFT accuracy and thread determinism");
}
