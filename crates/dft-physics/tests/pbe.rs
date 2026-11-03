use dft_geometry::Cell;
use dft_numerics::FftPlan3d;
use dft_physics::gga::{DENSITY_THRESHOLD, evaluate};
use ndarray::Array3;
use std::f64::consts::PI;
use xcx::{Functional, FunctionalId, Spin, XcInput};

// Independent published PBE96/PW92 formulas, no xcx helper calls.
fn oracle(n: f64, sigma: f64) -> [f64; 2] {
    let rs = (3.0 / (4.0 * PI * n)).cbrt();
    let kf = (3.0 * PI * PI * n).cbrt();
    let s2 = sigma / (4.0 * kf * kf * n * n);
    let ex = -0.75
        * (3.0 * n / PI).cbrt()
        * (1.0 + 0.804 - 0.804 / (1.0 + 0.219_514_972_764_517_1 * s2 / 0.804));
    let a = 0.0310907;
    let ec = -2.0
        * a
        * (1.0 + 0.21370 * rs)
        * (1.0
            + 1.0
                / (2.0
                    * a
                    * (7.5957 * rs.sqrt()
                        + 3.5876 * rs
                        + 1.6382 * rs.powf(1.5)
                        + 0.49294 * rs * rs)))
            .ln();
    let gamma = (1.0 - 2.0_f64.ln()) / (PI * PI);
    let beta = 0.066_724_550_603_149_22;
    let t2 = sigma / (4.0 * (4.0 * kf / PI) * n * n);
    let aa = beta / (gamma * ((-ec / gamma).exp_m1()));
    let h = gamma
        * (1.0 + beta / gamma * t2 * (1.0 + aa * t2) / (1.0 + aa * t2 + aa * aa * t2 * t2)).ln();
    [ex, ec + h]
}
#[test]
fn closed_form_oracle() {
    // Original seven spike points, supplied by Claude and recorded in ADR 0018.
    let points = [
        (0.01, 1e-6),
        (0.03, 4e-4),
        (0.03, 0.0),
        (0.1, 0.05),
        (0.2, 2.0),
        (0.3, 0.2),
        (1.0, 3.0),
    ];
    for (part, id) in [FunctionalId::GgaXPbe, FunctionalId::GgaCPbe]
        .into_iter()
        .enumerate()
    {
        let f = Functional::new(id, Spin::Unpolarized).unwrap();
        let mut errors = [0.0_f64; 3];
        for (n, s) in points {
            let r = f.eval(1, &XcInput::gga(&[n], &[s])).unwrap();
            let hn = n * 1e-5;
            // Scale the step by sigma or the natural s^2=1 gradient scale.
            // A fixed 1e-10 step at sigma=0 subtracts nearly identical energies
            // and loses the requested relative accuracy. Central differences
            // use the analytic extension to negative sigma near zero only.
            let sigma_scale = 4.0 * (3.0 * PI * PI * n).powf(2.0 / 3.0) * n * n;
            let hs = 1e-5 * s.max(sigma_scale);
            let dn = ((n + hn) * oracle(n + hn, s)[part] - (n - hn) * oracle(n - hn, s)[part])
                / (2.0 * hn);
            let ds = n * (oracle(n, s + hs)[part] - oracle(n, s - hs)[part]) / (2.0 * hs);
            let e = [
                (r.exc[0] - oracle(n, s)[part]).abs(),
                (r.vrho[0] - dn).abs() / r.vrho[0].abs(),
                (r.vsigma[0] - ds).abs() / r.vsigma[0].abs(),
            ];
            println!(
                "{id:?} n={n} sigma={s}: exc={:.16e} Ha, vrho={:.16e} Ha, vsigma={:.16e} Ha Bohr^5; errors={e:?}",
                r.exc[0], r.vrho[0], r.vsigma[0]
            );
            for j in 0..3 {
                errors[j] = errors[j].max(e[j]);
            }
            assert!(e[0] < 1e-14 && e[1] < 1e-8 && e[2] < 1e-8);
        }
        println!("{id:?} max errors={errors:?}");
    }
}
#[test]
fn variational_grid_translation_and_uniform() {
    for (dims, rows) in [
        (
            [12, 14, 16],
            [[6.0, 0.0, 0.0], [0.0, 7.0, 0.0], [0.0, 0.0, 8.0]],
        ),
        (
            [15, 12, 10],
            [[5.0, 0.0, 0.0], [1.2, 6.0, 0.0], [0.7, 0.9, 7.0]],
        ),
    ] {
        let cell = Cell::from_user_rows_bohr(rows).unwrap();
        let mut fft = FftPlan3d::new(dims).unwrap();
        let n = Array3::from_shape_fn(dims, |(i, j, k)| {
            let p = [
                i as f64 / dims[0] as f64,
                j as f64 / dims[1] as f64,
                k as f64 / dims[2] as f64,
            ];
            let mut n = 0.05;
            for center in [[0.2, 0.3, 0.4], [0.65, 0.6, 0.7]] {
                for a in -1..=1 {
                    for b in -1..=1 {
                        for c in -1..=1 {
                            let r = cell.fractional_to_cartesian_bohr(std::array::from_fn(|x| {
                                p[x] - center[x] - [a, b, c][x] as f64
                            }));
                            n += 0.1 * (-r.iter().map(|x| x * x).sum::<f64>()).exp();
                        }
                    }
                }
            }
            n
        });
        // Include Nyquist content in the perturbation as well as smooth modes.
        let dn = Array3::from_shape_fn(dims, |(i, j, k)| {
            0.01 * ((2.0 * PI * i as f64 / dims[0] as f64).cos()
                + (2.0 * PI * j as f64 / dims[1] as f64).sin()
                + if k % 2 == 0 { 0.3 } else { -0.3 })
        });
        let x = evaluate(&n, &cell, &mut fft).unwrap();
        let eps = 1e-4;
        let plus = &n + &dn.mapv(|x| eps * x);
        let minus = &n - &dn.mapv(|x| eps * x);
        let fd = (evaluate(&plus, &cell, &mut fft).unwrap().energy.0
            - evaluate(&minus, &cell, &mut fft).unwrap().energy.0)
            / (2.0 * eps);
        let analytic = cell.volume_bohr3() / n.len() as f64
            * dn.iter().zip(&x.potential).map(|(d, v)| d * v).sum::<f64>();
        let shifted = Array3::from_shape_fn(dims, |(i, j, k)| {
            n[[(i + 2) % dims[0], (j + 3) % dims[1], (k + 1) % dims[2]]]
        });
        let translation =
            (evaluate(&shifted, &cell, &mut fft).unwrap().energy.0 - x.energy.0).abs();
        println!(
            "PBE grid={dims:?}: E={:.16e} Ha, derivative error={:.9e} Ha, translation={translation:.9e} Ha",
            x.energy.0,
            (fd - analytic).abs()
        );
        assert!((fd - analytic).abs() < 1e-8);
        assert!(translation < 1e-12);
        let uniform = evaluate(&Array3::from_elem(dims, 0.2), &cell, &mut fft).unwrap();
        let expected = 0.2 * cell.volume_bohr3() * oracle(0.2, 0.0).iter().sum::<f64>();
        println!(
            "uniform energy error={:.9e} Ha",
            (uniform.energy.0 - expected).abs()
        );
        assert!((uniform.energy.0 - expected).abs() < 1e-12);
    }
}
#[test]
fn density_policy() {
    let cell =
        Cell::from_user_rows_bohr([[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]]).unwrap();
    let mut fft = FftPlan3d::new([4; 3]).unwrap();
    for n in [-1e-16, -1.0, f64::NAN, f64::INFINITY] {
        assert!(evaluate(&Array3::from_elem([4; 3], n), &cell, &mut fft).is_err());
    }
    for n in [0.0, DENSITY_THRESHOLD / 2.0, DENSITY_THRESHOLD] {
        let x = evaluate(&Array3::from_elem([4; 3], n), &cell, &mut fft).unwrap();
        assert_eq!(x.energy.0, 0.0);
        assert!(x.potential.iter().all(|v| *v == 0.0));
    }
    assert!(
        evaluate(
            &Array3::from_elem([4; 3], 2.0 * DENSITY_THRESHOLD),
            &cell,
            &mut fft
        )
        .unwrap()
        .energy
        .0 < 0.0
    );
}
