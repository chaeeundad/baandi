use super::*;
use ndarray::ShapeBuilder;
use std::f64::consts::TAU;

// Dimensionless tolerance for all relative max-norm and energy comparisons.
// Unit-amplitude mode, symmetry, and zero-leakage checks use absolute tolerance.
const TOL: f64 = 1e-12;
fn error(a: &Array3<Complex64>, b: &Array3<Complex64>) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).norm())
        .fold(0., f64::max)
}
fn norm(a: &Array3<Complex64>) -> f64 {
    a.iter().map(|z| z.norm()).fold(0., f64::max)
}
fn data(d: (usize, usize, usize)) -> Array3<Complex64> {
    Array3::from_shape_fn(d, |(i, j, k)| {
        let x = ((i * d.1 + j) * d.2 + k) as f64;
        Complex64::new((x * 0.371).sin(), (x * 0.217).cos())
    })
}
fn wave(d: (usize, usize, usize), q: [i32; 3]) -> Array3<Complex64> {
    Array3::from_shape_fn(d, |(i, j, k)| {
        Complex64::from_polar(
            1.,
            TAU * (q[0] as f64 * i as f64 / d.0 as f64
                + q[1] as f64 * j as f64 / d.1 as f64
                + q[2] as f64 * k as f64 / d.2 as f64),
        )
    })
}
#[test]
fn round_trips_parseval_and_reuse() {
    for d in [(4, 5, 6), (15, 20, 24), (16, 16, 16), (27, 27, 27)] {
        let mut plan = FftPlan3d::new([d.0, d.1, d.2]).unwrap();
        let input = data(d);
        let mut out = Array3::zeros(d);
        let mut back = out.clone();
        let mut repeat = out.clone();
        plan.to_real(&input, &mut out).unwrap();
        plan.to_reciprocal(&out, &mut back).unwrap();
        let g_round = error(&input, &back) / norm(&input);
        let energy_g: f64 = input.iter().map(|z| z.norm_sqr()).sum();
        let energy_r: f64 = out.iter().map(|z| z.norm_sqr()).sum();
        let parseval = (energy_r / input.len() as f64 - energy_g).abs() / energy_g;
        plan.to_real(&input, &mut repeat).unwrap();
        assert_eq!(out, repeat);
        plan.to_reciprocal(&input, &mut out).unwrap();
        plan.to_real(&out, &mut back).unwrap();
        let r_round = error(&input, &back) / norm(&input);
        plan.to_reciprocal(&input, &mut repeat).unwrap();
        assert_eq!(out, repeat);
        println!(
            "FFT {d:?}: G round={g_round:.16e} r round={r_round:.16e} Parseval={parseval:.16e} reuse=0"
        );
        assert!(g_round < TOL && r_round < TOL && parseval < TOL);
    }
}
#[test]
fn constant_functions() {
    // Dimensionless relative tolerance for G=0 and real-space values;
    // absolute tolerance for every non-origin coefficient: 1e-14.
    const CONSTANT_TOL: f64 = 1e-14;
    for d in [(4, 5, 6), (15, 20, 24)] {
        let mut plan = FftPlan3d::new([d.0, d.1, d.2]).unwrap();
        for value in [Complex64::new(1., 0.), Complex64::new(0.3, -0.7)] {
            let constant = Array3::from_elem(d, value);
            let mut coefficients = Array3::zeros(d);
            plan.to_reciprocal(&constant, &mut coefficients).unwrap();
            let origin_error = (coefficients[[0, 0, 0]] - value).norm() / value.norm();
            let leakage = coefficients
                .indexed_iter()
                .filter(|(index, _)| *index != (0, 0, 0))
                .map(|(_, z)| z.norm())
                .fold(0., f64::max);
            let mut origin_only = Array3::zeros(d);
            origin_only[[0, 0, 0]] = value;
            let mut real = Array3::zeros(d);
            plan.to_real(&origin_only, &mut real).unwrap();
            let real_error = error(&real, &constant) / value.norm();
            println!(
                "constant {value:?} grid={d:?}: origin_rel={origin_error:.16e} leakage_abs={leakage:.16e} real_rel={real_error:.16e}"
            );
            assert!(origin_error < CONSTANT_TOL);
            assert!(leakage < CONSTANT_TOL);
            assert!(real_error < CONSTANT_TOL);
        }
    }
}
#[test]
fn signed_modes() {
    let d = (4, 5, 6);
    let mut plan = FftPlan3d::new([4, 5, 6]).unwrap();
    let mut coefficients = Array3::zeros(d);
    let mut expected = Array3::zeros(d);
    let mut real = Array3::zeros(d);
    // Even-axis index 2 on length 4 represents -2, not +2.
    for (q, index) in [([1, -1, 2], [1, 4, 2]), ([-2, 1, -3], [2, 1, 3])] {
        let analytic = wave(d, q);
        expected.fill(Complex64::new(0., 0.));
        expected[index] = Complex64::new(1., 0.);
        plan.to_reciprocal(&analytic, &mut coefficients).unwrap();
        plan.to_real(&expected, &mut real).unwrap();
        let forward = error(&coefficients, &expected);
        let inverse = error(&real, &analytic);
        println!(
            "mode {q:?} index={index:?}: coefficient_abs={forward:.16e} real_abs={inverse:.16e}"
        );
        assert!(forward < TOL && inverse < TOL);
    }
}
#[test]
fn symmetry_and_linearity() {
    let d = (4, 5, 6);
    let mut plan = FftPlan3d::new([4, 5, 6]).unwrap();
    let a = data(d);
    let b = wave(d, [1, -1, 2]);
    let real = a.mapv(|z| Complex64::new(z.re, 0.));
    let mut out = Array3::zeros(d);
    plan.to_reciprocal(&real, &mut out).unwrap();
    let mut symmetry = 0.0_f64;
    for ((i, j, k), z) in out.indexed_iter() {
        symmetry = symmetry
            .max((*z - out[[(d.0 - i) % d.0, (d.1 - j) % d.1, (d.2 - k) % d.2]].conj()).norm());
    }
    assert!(symmetry < TOL);
    let alpha = Complex64::new(0.3, -0.7);
    let beta = Complex64::new(-0.2, 0.9);
    let sum = a.mapv(|z| alpha * z) + b.mapv(|z| beta * z);
    let mut fa = out.clone();
    let mut fb = out.clone();
    for inverse in [false, true] {
        plan.transform(&a, &mut fa, inverse).unwrap();
        plan.transform(&b, &mut fb, inverse).unwrap();
        plan.transform(&sum, &mut out, inverse).unwrap();
        let expected = fa.mapv(|z| alpha * z) + fb.mapv(|z| beta * z);
        let linearity = error(&out, &expected) / norm(&expected);
        println!("inverse={inverse}: linearity_rel={linearity:.16e} symmetry_abs={symmetry:.16e}");
        assert!(linearity < TOL);
    }
}
#[test]
fn band_limited_coefficients_across_grids() {
    let mut peaks = Vec::new();
    for d in [(4, 5, 6), (15, 20, 24)] {
        let mut plan = FftPlan3d::new([d.0, d.1, d.2]).unwrap();
        let f = wave(d, [1, -1, 2]).mapv(|z| z * Complex64::new(0.3, 0.2))
            + wave(d, [-1, 1, -1]).mapv(|z| z * Complex64::new(-0.7, 0.1));
        let mut out = Array3::zeros(d);
        plan.to_reciprocal(&f, &mut out).unwrap();
        let mut expected = Array3::zeros(d);
        expected[[1, d.1 - 1, 2]] = Complex64::new(0.3, 0.2);
        expected[[d.0 - 1, 1, d.2 - 1]] = Complex64::new(-0.7, 0.1);
        let err = error(&out, &expected) / norm(&expected);
        println!("band-limited {d:?}: coefficient_rel={err:.16e}");
        assert!(err < TOL);
        peaks.push([out[[1, d.1 - 1, 2]], out[[d.0 - 1, 1, d.2 - 1]]]);
    }
    for (a, b) in peaks[0].iter().zip(&peaks[1]) {
        assert!((a - b).norm() / a.norm() < TOL);
    }
}
#[test]
fn grid_layout_and_errors() {
    assert!(matches!(
        FftPlan3d::new([0, 1, 1]),
        Err(FftError::InvalidGrid)
    ));
    assert!(matches!(
        FftPlan3d::new([usize::MAX, 2, 2]),
        Err(FftError::InvalidGrid)
    ));
    let mut plan = FftPlan3d::new([4, 5, 6]).unwrap();
    assert_eq!(plan.shape(), [4, 5, 6]);
    let input = data((4, 5, 6));
    for ((i, j, k), z) in input.indexed_iter() {
        assert_eq!(*z, input.as_slice().unwrap()[(i * 5 + j) * 6 + k]);
    }
    let mut bad = Array3::zeros((4, 6, 5));
    assert_eq!(
        plan.to_real(&input, &mut bad),
        Err(FftError::IncompatibleArray)
    );
    let fortran = Array3::zeros((4, 5, 6).f());
    let mut out = Array3::zeros((4, 5, 6));
    assert_eq!(
        plan.to_reciprocal(&fortran, &mut out),
        Err(FftError::IncompatibleArray)
    );
    let mut bad_layout = fortran;
    assert_eq!(
        plan.to_real(&input, &mut bad_layout),
        Err(FftError::IncompatibleArray)
    );
    let mut single = FftPlan3d::new([1, 1, 1]).unwrap();
    let input = data((1, 1, 1));
    let mut out = input.clone();
    single.to_real(&input, &mut out).unwrap();
    assert_eq!(out, input);
    single.to_reciprocal(&input, &mut out).unwrap();
    assert_eq!(out, input);
}
