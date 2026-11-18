use dft_geometry::Cell;
use dft_physics::{data::*, gga, hartree, spin::*};
use std::f64::consts::PI;
#[test]
fn discrete_spin_functional_derivative_and_workspace_reuse() {
    for dims in [[12, 14, 16], [15, 12, 10]] {
        let cell =
            Cell::from_user_rows_bohr([[6.0, 0.0, 0.0], [1.0, 7.0, 0.0], [0.5, 1.0, 9.0]]).unwrap();
        let rho = std::array::from_fn(|s| {
            Array3::from_shape_fn(dims, |(i, j, k)| {
                0.03 + 0.01 * s as f64
                    + 0.004 * (2.0 * PI * i as f64 / dims[0] as f64).sin()
                    + 0.003 * ((s + 1) as f64 * 2.0 * PI * j as f64 / dims[1] as f64).cos()
                    + 0.002 * (2.0 * PI * k as f64 / dims[2] as f64).sin()
            })
        });
        let delta: [Array3<f64>; 2] = std::array::from_fn(|s| {
            Array3::from_shape_fn(dims, |(i, j, _)| {
                0.01 * ((s + 1) as f64 * 2.0 * PI * i as f64 / dims[0] as f64).sin()
                    + 0.005 * (2.0 * PI * j as f64 / dims[1] as f64).cos()
            })
        });
        for kind in [SpinFunctional::LdaPw, SpinFunctional::Pbe] {
            let mut xc = SpinXc::new(kind, &cell, dims).unwrap();
            let mut response = SpinXcResponse::new(dims);
            let pointers = response.potentials.each_ref().map(|n| n.as_ptr());
            xc.evaluate(&rho, &mut response).unwrap();
            let weight = cell.volume_bohr3() / rho[0].len() as f64;
            let expected = weight
                * (0..2)
                    .map(|s| {
                        delta[s]
                            .iter()
                            .zip(&response.potentials[s])
                            .map(|(d, v)| d * v)
                            .sum::<f64>()
                    })
                    .sum::<f64>();
            let base_energy = response.energy_ha;
            let mut slopes = [0.0; 3];
            for (i, h) in [1e-2, 5e-3, 2.5e-3].into_iter().enumerate() {
                let plus = std::array::from_fn(|s| &rho[s] + &delta[s].mapv(|v| h * v));
                let minus = std::array::from_fn(|s| &rho[s] - &delta[s].mapv(|v| h * v));
                xc.evaluate(&plus, &mut response).unwrap();
                let ep = response.energy_ha;
                xc.evaluate(&minus, &mut response).unwrap();
                slopes[i] = (ep - response.energy_ha) / (2.0 * h);
            }
            let cap = 1e-10_f64.max(1e-8 * expected.abs());
            let error = ((4.0 * slopes[2] - slopes[1]) / 3.0 - expected).abs();
            assert!(error < cap, "{kind:?} {dims:?}: {error:e} vs {cap:e}");
            let d0 = (slopes[0] - slopes[1]).abs();
            let d1 = (slopes[1] - slopes[2]).abs();
            assert!(
                d1 <= 0.4 * d0
                    || (slopes.iter().all(|v| (v - expected).abs() < cap)
                        && d0.max(d1) < cap / 3.0)
            );
            xc.evaluate(&rho, &mut response).unwrap();
            assert_eq!(base_energy, response.energy_ha);
            assert_eq!(pointers, response.potentials.each_ref().map(|n| n.as_ptr()));
            println!("{kind:?} {dims:?}: grid directional derivative error={error:e} Ha");
        }
        let total = &rho[0] + &rho[1];
        let mut work = SpinHartree::new(&cell, dims).unwrap();
        work.evaluate(&total).unwrap();
        let reference =
            hartree::evaluate(&total, &cell, &mut FftPlan3d::new(dims).unwrap()).unwrap();
        assert_eq!(work.energy_ha.to_bits(), reference.energy.0.to_bits());
        for (a, b) in work.potential.iter().zip(&reference.potential) {
            assert_eq!(a.to_bits(), b.re.to_bits());
        }
    }
}
#[test]
fn grid_zero_spin_reduction_and_density_tag() {
    let dims = [8, 10, 12];
    let cell =
        Cell::from_user_rows_bohr([[5.0, 0.0, 0.0], [0.5, 6.0, 0.0], [0.0, 0.2, 8.0]]).unwrap();
    let n = Array3::from_shape_fn(dims, |(i, j, _)| {
        0.1 + 0.01 * (2.0 * PI * i as f64 / 8.0).sin() + 0.02 * (2.0 * PI * j as f64 / 10.0).cos()
    });
    let channels = [n.mapv(|v| v / 2.0), n.mapv(|v| v / 2.0)];
    let tagged = SpinDensity::Collinear {
        grid_id: "density".into(),
        channels: channels.clone(),
    };
    let mut total = Array3::zeros(dims);
    tagged.total_into(&mut total).unwrap();
    assert_eq!(total, n);
    let mut magnetization = Array3::zeros(dims);
    tagged.magnetization_into(&mut magnetization).unwrap();
    assert!(magnetization.iter().all(|&v| v == 0.0));
    for kind in [SpinFunctional::LdaPw, SpinFunctional::Pbe] {
        let mut xc = SpinXc::new(kind, &cell, dims).unwrap();
        let mut response = SpinXcResponse::new(dims);
        xc.evaluate(&channels, &mut response).unwrap();
        let u = if kind == SpinFunctional::LdaPw {
            unpolarized_pw(&n, cell.volume_bohr3()).unwrap()
        } else {
            gga::evaluate(&n, &cell, &mut FftPlan3d::new(dims).unwrap()).unwrap()
        };
        assert!((u.energy.0 - response.energy_ha).abs() < 1e-10);
        for s in 0..2 {
            for (a, b) in response.potentials[s].iter().zip(&u.potential) {
                assert!((a - b).abs() < 1e-10);
            }
        }
    }
    assert!(
        SpinDensity::Collinear {
            grid_id: "density".into(),
            channels: [Array3::zeros(dims), Array3::zeros([1, 2, 3])]
        }
        .validate()
        .is_err()
    );
    assert!(
        SpinDensity::Unpolarized {
            grid_id: "".into(),
            total: n
        }
        .validate()
        .is_err()
    );
}
