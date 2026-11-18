use dft_physics::pbe_pw92::{Pw92Prefactor as A, correlation};
use std::f64::consts::PI;
use xcx::{Functional, FunctionalId, Spin, XcInput};

#[test]
fn dense_random_library_and_derivative_oracles() {
    let mut state = 270027_u64;
    let mut random = || {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        (state >> 11) as f64 / ((1_u64 << 53) as f64)
    };
    let mut n = Vec::new();
    let mut sigma = Vec::new();
    for i in 0..20000 {
        let rs = 0.1 * 200_f64.powf(random());
        let rho = 3.0 / (4.0 * PI * rs.powi(3));
        let kf = (3.0 * PI * PI * rho).cbrt();
        let reduced = if i % 100 == 0 { 0.0 } else { 3.0 * random() };
        n.push(rho);
        sigma.push(4.0 * kf * kf * rho * rho * reduced * reduced);
    }
    let c = Functional::new(FunctionalId::GgaCPbe, Spin::Unpolarized)
        .unwrap()
        .eval(n.len(), &XcInput::gga(&n, &sigma))
        .unwrap();
    let mut maxima = [0.0_f64; 3];
    let mut fdmax = [0.0_f64; 2];
    for (i, (&n, &s)) in n.iter().zip(&sigma).enumerate() {
        let point = correlation(n, s, A::Modified);
        for (j, (actual, expected)) in [
            (point.exc, c.exc[i]),
            (point.vrho, c.vrho[i]),
            (point.vsigma, c.vsigma[i]),
        ]
        .into_iter()
        .enumerate()
        {
            let rel = (actual - expected).abs() / expected.abs().max(1e-30);
            maxima[j] = maxima[j].max(rel);
            assert!(
                rel <= 1e-12,
                "i={i} component={j} n={n} sigma={s} rel={rel:e}"
            );
        }
        let p = correlation(n, s, A::Original);
        // Central differences: h ~ eps^(1/3), with natural sigma scale to
        // resolve cancellation at sigma=0. Relative cap 1e-7 balances O(h^2)
        // truncation and f64 subtraction; fixed before observing results.
        let hn = n * 1e-5;
        let kf = (3.0 * PI * PI * n).cbrt();
        let hs = (s + 16.0 * kf * n * n / PI) * 1e-5;
        let energy = |n, s| n * correlation(n, s, A::Original).exc;
        let fd = [
            (energy(n + hn, s) - energy(n - hn, s)) / (2.0 * hn),
            (energy(n, s + hs) - energy(n, s - hs)) / (2.0 * hs),
        ];
        for (j, actual) in [p.vrho, p.vsigma].into_iter().enumerate() {
            let err = (fd[j] - actual).abs() / actual.abs().max(1e-30);
            fdmax[j] = fdmax[j].max(err);
            assert!(err < 1e-7, "FD component={j} n={n} sigma={s} err={err:e}");
        }
    }
    println!(
        "20000 points rs=0.1..20 s=0..3: relative maxima exc/vrho/vsigma={maxima:?}; original FD={fdmax:?}"
    );
}

#[test]
fn grid_potential_and_variational_derivative() {
    use dft_geometry::Cell;
    use dft_numerics::FftPlan3d;
    use dft_physics::gga;
    use ndarray::Array3;
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [1., 6., 0.], [0.5, 1., 7.]]).unwrap();
    for dims in [[8, 10, 12], [9, 7, 5]] {
        let n = Array3::from_shape_fn(dims, |(i, j, k)| {
            0.1 + 0.02 * (2. * PI * i as f64 / dims[0] as f64).cos() + 0.01 * ((j + k) as f64).sin()
        });
        let d = Array3::from_shape_fn(dims, |(i, j, k)| ((i + 2 * j + 3 * k) as f64).cos());
        let mut fft = FftPlan3d::new(dims).unwrap();
        let old = gga::evaluate(&n, &cell, &mut fft).unwrap();
        let new = gga::evaluate_pw92(A::Modified, &n, &cell, &mut fft).unwrap();
        let verr = old
            .potential
            .iter()
            .zip(&new.potential)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0_f64, f64::max);
        assert!((old.energy.0 - new.energy.0).abs() < 1e-12);
        assert!(verr < 1e-12);
        let orig = gga::evaluate_pw92(A::Original, &n, &cell, &mut fft).unwrap();
        let analytic = cell.volume_bohr3() / n.len() as f64 * (&orig.potential * &d).sum();
        let h = 1e-6;
        let plus = gga::evaluate_pw92(A::Original, &(&n + &(&d * h)), &cell, &mut fft)
            .unwrap()
            .energy
            .0;
        let minus = gga::evaluate_pw92(A::Original, &(&n - &(&d * h)), &cell, &mut fft)
            .unwrap()
            .energy
            .0;
        let err = ((plus - minus) / (2. * h) - analytic).abs();
        assert!(err < 1e-7 * analytic.abs().max(1.));
        println!(
            "grid {dims:?}: energy error={} Ha, potential max={verr:e} Ha, FD={err:e} Ha",
            old.energy.0 - new.energy.0
        );
    }
}
