mod common;
use dft_upf::radial::{
    LocalTransform, RadialMesh, projector_transform, real_spherical_harmonic as y,
    spherical_bessel as j,
};
use std::f64::consts::PI;
fn logmesh(n: usize, end: f64) -> RadialMesh {
    let (r, rab) = common::mesh(n, end);
    RadialMesh::new(r, rab).unwrap()
}
#[test]
fn simpson_convergence_and_parity() {
    for n in [101, 201, 401, 801, 802] {
        let mesh = logmesh(n, 20.0);
        let exp: Vec<_> = mesh.r().iter().map(|r| (-r).exp()).collect();
        let gauss: Vec<_> = mesh.r().iter().map(|r| (-r * r).exp()).collect();
        let e = (mesh.integrate(&exp).unwrap() - (1.0 - (-20.0_f64).exp())).abs();
        let g = (mesh.integrate(&gauss).unwrap() - PI.sqrt() / 2.0).abs();
        println!("Simpson n={n}: exponential absolute error={e:.16e} Bohr, Gaussian={g:.16e} Bohr");
        assert!(e < 2e-7 && g < 2e-7);
        if n >= 801 {
            assert!(e < 1e-10 && g < 1e-10);
        }
    }
    for n in 2..10 {
        let mesh = RadialMesh::new((0..n).map(|i| i as f64).collect(), vec![1.0; n]).unwrap();
        let end = (n - 1) as f64;
        for power in 0..=if n == 2 { 1 } else { 3 } {
            let v: Vec<_> = mesh.r().iter().map(|r| r.powi(power)).collect();
            assert!(
                (mesh.integrate(&v).unwrap() - end.powi(power + 1) / (power + 1) as f64).abs()
                    < 1e-12
            );
        }
    }
}
#[test]
fn bessel_identities() {
    let mut worst = 0.0_f64;
    for x in [0.2_f64, 0.8, 1.5, 3.9, 4.0, 4.1, 7.1, 20.2, 100.3, 10000.1] {
        let exact = [
            x.sin() / x,
            x.sin() / x.powi(2) - x.cos() / x,
            (3.0 / x.powi(3) - 1.0 / x) * x.sin() - 3.0 / x.powi(2) * x.cos(),
            (15.0 / x.powi(4) - 6.0 / x.powi(2)) * x.sin() + (1.0 / x - 15.0 / x.powi(3)) * x.cos(),
        ];
        // Elementary closed forms cancel at small x for l>=2; use them only at x>=4.
        for (l, expected) in exact.into_iter().enumerate() {
            if l < 2 || x >= 4.0 {
                let err = ((j(l, x).unwrap() - expected) / expected).abs();
                worst = worst.max(err);
                assert!(err < 1e-13, "l={l} x={x} err={err}");
            }
        }
        for l in 1..=2 {
            let lhs = j(l - 1, x).unwrap() + j(l + 1, x).unwrap();
            let rhs = (2 * l + 1) as f64 / x * j(l, x).unwrap();
            let scale = j(l - 1, x).unwrap().abs() + j(l + 1, x).unwrap().abs() + rhs.abs();
            assert!((lhs - rhs).abs() / scale < 1e-13);
        }
    }
    for l in 0..=3 {
        let df = (1..=l).fold(1.0, |p, k| p * (2 * k + 1) as f64);
        assert_eq!(j(l, 0.0).unwrap(), if l == 0 { 1.0 } else { 0.0 });
        for x in [1e-8_f64, 1e-20, 1e-90] {
            let limit = x.powi(l as i32) / df;
            assert!((j(l, x).unwrap() / limit - 1.0).abs() < 1e-13);
        }
    }
    println!("Bessel closed forms max relative error={worst:.16e}");
}
// Independent test quadrature: Newton roots of Legendre P_n, n=12.
fn gauss_legendre(n: usize) -> Vec<(f64, f64)> {
    (0..n)
        .map(|i| {
            let mut x = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
            let mut derivative = 0.0;
            for _ in 0..30 {
                let (mut p0, mut p1) = (1.0, x);
                for k in 2..=n {
                    let p = ((2 * k - 1) as f64 * x * p1 - (k - 1) as f64 * p0) / k as f64;
                    p0 = p1;
                    p1 = p;
                }
                derivative = n as f64 * (x * p1 - p0) / (x * x - 1.0);
                let dx = p1 / derivative;
                x -= dx;
                if dx.abs() < 1e-16 {
                    break;
                }
            }
            (x, 2.0 / ((1.0 - x * x) * derivative * derivative))
        })
        .collect()
}
#[test]
fn harmonic_orthonormality_and_addition() {
    let indices: Vec<_> = (0..=3)
        .flat_map(|l| (-(l as isize)..=l as isize).map(move |m| (l, m)))
        .collect();
    let mut matrix = [[0.0; 16]; 16];
    for (z, w) in gauss_legendre(12) {
        for k in 0..32 {
            let phi = 2.0 * PI * k as f64 / 32.0;
            let s = (1.0 - z * z).sqrt();
            let a = [s * phi.cos(), s * phi.sin(), z];
            let values: Vec<_> = indices.iter().map(|(l, m)| y(*l, *m, a).unwrap()).collect();
            for i in 0..16 {
                for j in 0..16 {
                    matrix[i][j] += w * 2.0 * PI / 32.0 * values[i] * values[j];
                }
            }
        }
    }
    let mut worst = 0.0_f64;
    for (i, row) in matrix.iter().enumerate() {
        for (j, v) in row.iter().enumerate() {
            worst = worst.max((v - if i == j { 1.0 } else { 0.0 }).abs());
        }
    }
    assert!(worst < 1e-12);
    let mut seed = 91_u64;
    let mut direction = || {
        let mut a = [0.0; 3];
        for v in &mut a {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            *v = 2.0 * ((seed >> 11) as f64 / (1_u64 << 53) as f64) - 1.0;
        }
        let n = a.iter().map(|v| v * v).sum::<f64>().sqrt();
        a.map(|v| v / n)
    };
    let mut addition = 0.0_f64;
    for _ in 0..100 {
        let a = direction();
        let b = direction();
        let dot = a.iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
        let legendre = [
            1.0,
            dot,
            (3.0 * dot * dot - 1.0) / 2.0,
            (5.0 * dot.powi(3) - 3.0 * dot) / 2.0,
        ];
        for (l, p) in legendre.into_iter().enumerate() {
            let sum = (-(l as isize)..=l as isize)
                .map(|m| y(l, m, a).unwrap() * y(l, m, b).unwrap())
                .sum::<f64>();
            addition = addition.max((sum - (2 * l + 1) as f64 / (4.0 * PI) * p).abs());
        }
    }
    assert!(addition < 1e-13);
    // Signs/order: p_x,p_y,p_z are positive along their Cartesian axes.
    let norm = (3.0 / (4.0 * PI)).sqrt();
    assert!((y(1, 1, [1.0, 0.0, 0.0]).unwrap() - norm).abs() < 1e-15);
    assert!((y(1, -1, [0.0, 1.0, 0.0]).unwrap() - norm).abs() < 1e-15);
    assert!((y(1, 0, [0.0, 0.0, 1.0]).unwrap() - norm).abs() < 1e-15);
    println!(
        "Y_lm 12x32 quadrature max absolute error={worst:.16e}; addition theorem={addition:.16e}"
    );
}
fn potential(mesh: &RadialMesh, z: f64, rc: f64, a: f64, s: f64) -> Vec<f64> {
    mesh.r()
        .iter()
        .map(|r| {
            let coulomb = if *r == 0.0 {
                -2.0 * z / (PI.sqrt() * rc)
            } else {
                -z * libm::erf(r / rc) / r
            };
            coulomb + a * (-r * r / (s * s)).exp()
        })
        .collect()
}
#[test]
fn local_closed_forms_and_zero_bookkeeping() {
    let mesh = logmesh(12001, 14.0);
    let z = 4.0;
    let volume = 120.0;
    for rc in [0.3_f64, 0.5] {
        for a in [0.0, -2.0] {
            let s = 0.4;
            let v = potential(&mesh, z, rc, a, s);
            let t = LocalTransform::new(&mesh, &v, z, volume, None).unwrap();
            let mut worst = 0.0_f64;
            for k in 1..=120 {
                let g = k as f64 * 0.1;
                // Independent 3D Gaussian Fourier transform and Poisson solution.
                let exact = (-4.0 * PI * z * (-g * g * rc * rc / 4.0).exp() / (g * g)
                    + a * PI.powf(1.5) * s.powi(3) * (-g * g * s * s / 4.0).exp())
                    / volume;
                let err = ((t.at_g(g).unwrap() - exact) / exact).abs();
                worst = worst.max(err);
                assert!(err < 1e-8, "rc={rc} A={a} g={g} error={err}");
            }
            let alpha = (PI * z * rc * rc + a * PI.powf(1.5) * s.powi(3)) / volume;
            let err = (t.alpha_z().unwrap() - alpha).abs();
            assert!(err < 1e-11);
            let direct: Vec<_> = mesh
                .r()
                .iter()
                .zip(&v)
                .map(|(r, v)| 4.0 * PI / volume * (r * r * v + z * r))
                .collect();
            assert!((mesh.integrate(&direct).unwrap() - t.alpha_z().unwrap()).abs() < 1e-11);
            println!(
                "Local rc={rc} Bohr A={a} Ha, G=0.1..12 Bohr^-1: max relative={worst:.16e}, alpha-Z={:.16e} Ha, alpha error={err:.16e} Ha",
                t.alpha_z().unwrap()
            );
        }
    }
    let v = potential(&mesh, z, 1.2, -0.3, 0.7);
    let exact = (-4.0 * PI * z * (-1.2_f64.powi(2) / 4.0).exp()
        - 0.3 * PI.powf(1.5) * 0.7_f64.powi(3) * (-0.7_f64.powi(2) / 4.0).exp())
        / volume;
    let mut previous = f64::INFINITY;
    for cutoff in [3.0, 5.0, 8.0, 10.0] {
        let t = LocalTransform::new(&mesh, &v, z, volume, Some(cutoff)).unwrap();
        let error = (t.at_g(1.0).unwrap() - exact).abs();
        println!(
            "Local r_max={cutoff} Bohr effective={:.8} Bohr: G=1 absolute error={error:.16e} Ha",
            t.effective_r_max()
        );
        assert!(error <= previous + 1e-14);
        previous = error;
    }
    assert!(previous < 1e-11);
}
#[test]
fn gaussian_projectors_and_parseval() {
    let mesh = logmesh(4001, 12.0);
    let mut worst = 0.0_f64;
    for s in [0.4_f64, 0.8] {
        for l in 0..=2 {
            let u: Vec<_> = mesh
                .r()
                .iter()
                .map(|r| r.powi(l as i32 + 1) * (-r * r / (2.0 * s * s)).exp())
                .collect();
            for q in [0.0_f64, 0.1, 0.7, 1.0, 3.0, 6.0] {
                // Gaussian solid-harmonic FT from differentiation of Cartesian Gaussian.
                let exact = 4.0
                    * PI
                    * (PI / 2.0).sqrt()
                    * s.powi(2 * l as i32 + 3)
                    * q.powi(l as i32)
                    * (-s * s * q * q / 2.0).exp();
                let got = projector_transform(&mesh, l, &u, q).unwrap();
                if exact == 0.0 {
                    assert_eq!(got, 0.0);
                } else {
                    let err = ((got - exact) / exact).abs();
                    worst = worst.max(err);
                    assert!(err < 1e-9, "l={l} s={s} q={q} err={err}");
                }
            }
            // Independent linear q mesh, integrate measured transform, not its closed form.
            let nq = 401;
            let qmax = 12.0 / s;
            let dq = qmax / (nq - 1) as f64;
            let qmesh =
                RadialMesh::new((0..nq).map(|i| i as f64 * dq).collect(), vec![dq; nq]).unwrap();
            let integrand: Vec<_> = qmesh
                .r()
                .iter()
                .map(|q| q * q * projector_transform(&mesh, l, &u, *q).unwrap().powi(2))
                .collect();
            let left = 2.0 / PI * qmesh.integrate(&integrand).unwrap();
            // Integral r^(2l+2) exp(-r^2/s^2) dr = s^(2l+3) Gamma(l+3/2)/2.
            let gamma = (0..=l).fold(PI.sqrt(), |g, k| g * (k as f64 + 0.5));
            let right = 16.0 * PI * PI * 0.5 * s.powi(2 * l as i32 + 3) * gamma;
            let error = (left / right - 1.0).abs();
            assert!(error < 1e-9);
            println!("Projector Parseval l={l}, s={s} Bohr relative error={error:.16e}");
        }
    }
    println!("Gaussian projector max relative error={worst:.16e}");
}
#[test]
fn invalid_radial_arguments() {
    assert!(RadialMesh::new(vec![0.0, 1.0], vec![0.1, -1.0]).is_err());
    assert!(RadialMesh::new(vec![0.0, f64::NAN], vec![0.1, 1.0]).is_err());
    assert!(j(4, 1.0).is_err());
    assert!(j(0, -1.0).is_err());
    assert!(y(1, 2, [1.0, 0.0, 0.0]).is_err());
    assert!(y(0, 0, [0.0; 3]).is_err());
    let mesh = logmesh(101, 12.0);
    assert!(mesh.integrate(&[1.0]).is_err());
    let v = vec![0.0; 101];
    let t = LocalTransform::new(&mesh, &v, 1.0, 1.0, None).unwrap();
    assert!(t.at_g(0.0).is_err());
    assert!(t.at_g(f64::NAN).is_err());
    assert!(LocalTransform::new(&mesh, &v, 1.0, 0.0, None).is_err());
    assert!(projector_transform(&mesh, 4, &v, 1.0).is_err());
    // Public parsed transforms use exactly the documented radial storage convention.
    let p = dft_upf::parse_upf(&common::fixture()).unwrap();
    assert!(
        p.local_transform(100.0, None)
            .unwrap()
            .at_g(1.0)
            .unwrap()
            .is_finite()
    );
    assert_eq!(
        p.projector_transform(0, 1.0).unwrap(),
        projector_transform(p.mesh(), 0, p.projectors()[0].radial_times_r(), 1.0).unwrap()
    );
    assert!(p.projector_transform(3, 1.0).is_err());
}
