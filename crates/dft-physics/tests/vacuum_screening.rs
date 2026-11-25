use dft_geometry::Cell;
use dft_numerics::{FftPlan3d, switches::RunConfig};
use dft_physics::{
    data::Array3,
    gga,
    spin::{SpinFunctional, SpinXc, SpinXcResponse},
    xc,
};

fn cell() -> Cell {
    Cell::from_user_rows_bohr([[8., 0., 0.], [0.3, 7., 0.], [0.2, 0.1, 6.]]).unwrap()
}
fn config(on: bool) -> RunConfig {
    RunConfig::from_pairs(&[("xc_screen_negative", if on { "on" } else { "off" })]).unwrap()
}

#[test]
fn signed_vacuum_is_xc_only_and_budgeted() {
    let dims = [4, 3, 5];
    let mut n = Array3::from_elem(dims, 0.1);
    n[[1, 2, 3]] = -2e-7;
    let original = n.clone();
    let c = cell();
    let mut fft = FftPlan3d::new(dims).unwrap();
    {
        let _scope = config(false).enter();
        let e = gga::evaluate(&n, &c, &mut fft).err().unwrap().to_string();
        assert!(
            e.contains("channel total") && e.contains("grid point 28"),
            "{e}"
        );
    }
    let _scope = config(true).enter();
    xc::clear_screening_records();
    let out = gga::evaluate(&n, &c, &mut fft).unwrap();
    assert_eq!(n, original);
    assert_eq!(out.potential[[1, 2, 3]], 0.);
    let r = xc::screening_records();
    assert_eq!(r[0].points, 1);
    assert_eq!(r[0].minimum_grid_index, 28);
    assert!((r[0].charge_electrons - 2e-7 * c.volume_bohr3() / 60.).abs() < 1e-20);
    n[[1, 2, 3]] = -1.;
    assert!(gga::evaluate(&n, &c, &mut fft).is_err());
    n[[1, 2, 3]] = f64::NAN;
    assert!(gga::evaluate(&n, &c, &mut fft).is_err());
}

#[test]
fn positive_fields_keep_all_bits_and_screened_pbe_is_variational() {
    let dims = [4, 3, 5];
    let c = cell();
    let mut fft = FftPlan3d::new(dims).unwrap();
    let mut n = Array3::from_shape_fn(dims, |(i, j, k)| {
        0.1 + 0.01 * ((i + 2 * j + 3 * k) as f64).sin()
    });
    let off = {
        let _s = config(false).enter();
        gga::evaluate(&n, &c, &mut fft).unwrap()
    };
    let _s = config(true).enter();
    let on = gga::evaluate(&n, &c, &mut fft).unwrap();
    assert_eq!(off.energy.0.to_bits(), on.energy.0.to_bits());
    assert_eq!(off.potential, on.potential);
    n[[1, 2, 3]] = -2e-7;
    let v = gga::evaluate(&n, &c, &mut fft).unwrap();
    let direction = Array3::from_shape_fn(dims, |(i, j, k)| {
        if (i, j, k) == (1, 2, 3) {
            0.
        } else {
            ((2 * i + j + 4 * k) as f64).cos()
        }
    });
    let h = 1e-5;
    let ep = gga::evaluate(&(&n + &(&direction * h)), &c, &mut fft)
        .unwrap()
        .energy
        .0;
    let em = gga::evaluate(&(&n - &(&direction * h)), &c, &mut fft)
        .unwrap()
        .energy
        .0;
    let analytic = c.volume_bohr3() / 60. * (&v.potential * &direction).sum();
    let error = ((ep - em) / (2. * h) - analytic).abs();
    eprintln!("screened PBE derivative error={error:e} Ha");
    assert!(error < 1e-8, "{error:e}");
}

#[test]
fn spin_core_split_and_lda_screens() {
    let dims = [4, 3, 5];
    let c = cell();
    let mut channels = [Array3::from_elem(dims, 0.06), Array3::from_elem(dims, 0.04)];
    let mut core = Array3::from_elem(dims, 0.01);
    channels[0][[1, 2, 3]] = 1e-8;
    channels[1][[1, 2, 3]] = 2e-8;
    core[[1, 2, 3]] = -4.2e-7;
    let _s = config(true).enter();
    for kind in [SpinFunctional::LdaPw, SpinFunctional::Pbe] {
        let mut xc = SpinXc::new(kind, &c, dims).unwrap();
        let mut out = SpinXcResponse::new(dims);
        xc.evaluate_with_core(&channels, Some(&core), &mut out)
            .unwrap();
        assert_eq!(out.potentials[0][[1, 2, 3]], 0.);
        assert_eq!(out.potentials[1][[1, 2, 3]], 0.);
        assert!(out.energy_ha.is_finite());
    }
    let mut total = Array3::from_elem(dims, 0.1);
    total[[1, 2, 3]] = -2e-7;
    assert_eq!(
        xc::evaluate(&total, c.volume_bohr3()).unwrap().potential[[1, 2, 3]],
        0.
    );
    assert_eq!(
        dft_physics::spin::unpolarized_pw(&total, c.volume_bohr3())
            .unwrap()
            .potential[[1, 2, 3]],
        0.
    );
}

#[test]
fn fixed_row_statistics_are_identical_at_one_three_eight_ranks() {
    use dft_mpi::Collective;
    let dims = [4, 3, 5];
    let n: Vec<_> = (0..60)
        .map(|i| {
            if i % 7 == 0 {
                -2e-7
            } else {
                0.1 + 1e-5 * i as f64
            }
        })
        .collect();
    let serial = xc::screening_rows(&n, &(0..60).collect::<Vec<_>>(), dims);
    for ranks in [1, 3, 8] {
        dft_mpi::simulate(ranks, |comm| {
            let indices: Vec<_> = (0..60).filter(|i| (i / 5) % ranks == comm.rank()).collect();
            let local: Vec<_> = indices.iter().map(|&i| n[i]).collect();
            let row = xc::screening_rows(&local, &indices, dims);
            let mut all = vec![0.; row.len()];
            comm.tree_sum_f64(&row, &mut all).unwrap();
            assert_eq!(all, serial);
        });
    }
}

#[test]
fn positive_spin_core_fields_keep_all_bits() {
    let dims = [4, 3, 5];
    let c = cell();
    let density = [Array3::from_elem(dims, 0.06), Array3::from_elem(dims, 0.04)];
    let core = Array3::from_shape_fn(dims, |(i, j, k)| 0.01 + 0.001 * ((i + j + k) as f64).sin());
    for kind in [SpinFunctional::LdaPw, SpinFunctional::Pbe] {
        let mut old = None;
        for on in [false, true] {
            let _s = config(on).enter();
            let mut adapter = SpinXc::new(kind, &c, dims).unwrap();
            let mut out = SpinXcResponse::new(dims);
            adapter
                .evaluate_with_core(&density, Some(&core), &mut out)
                .unwrap();
            let bits = (
                out.energy_ha.to_bits(),
                out.valence_potential_integral_ha.to_bits(),
                out.core_potential_integral_ha.to_bits(),
                out.potentials,
            );
            if let Some(old) = &old {
                assert_eq!(&bits, old);
            } else {
                old = Some(bits);
            }
        }
    }
}
