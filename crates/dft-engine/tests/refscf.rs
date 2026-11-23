//! Independent analytic oracles for issue #25; synthetic data are never QE evidence.
use dft_engine::refscf::{
    Error,
    broyden::Broyden,
    density::{DensityPlan, half},
    metric::{Metric, Serial},
    profile::{Preconditioning, Profile, next_ethr_ry},
    qe_reference::{QeDensity, printed_accuracy},
};
use dft_mpi::Collective;
use dft_physics::data::{Array3, Complex64 as C};
use std::f64::consts::{PI, TAU};
fn c(v: f64) -> C {
    C::new(v, 0.)
}
fn metric(g2: &[f64], smooth: &[bool]) -> Metric {
    Metric::new(1. / (8. * PI), g2, smooth, 1, 350).unwrap()
}
#[test]
fn two_dimensional_linear_fixed_point_closed_form() {
    // F(x)=(2,3)-diag(2,3)x, unit weights. Exact fixed point is (1,1).
    let p = Profile {
        beta: 0.5,
        ..Default::default()
    };
    let mut mix = Broyden::new(p, &[1., 1.], &[true; 2], 1., 2.).unwrap();
    let mut metric = metric(&[1., 1.], &[true; 2]);
    let mut x = vec![c(0.), c(0.)];
    let mut next = x.clone();
    let d = mix
        .update(&Serial, &mut metric, &x, &[c(2.), c(3.)], &mut next)
        .unwrap();
    assert_eq!(next, vec![c(1.), c(1.5)]);
    assert_eq!(d.occupancy, 0);
    x = next.clone();
    let d = mix
        .update(&Serial, &mut metric, &x, &[c(0.), c(-1.5)], &mut next)
        .unwrap();
    assert_eq!(d.occupancy, 1);
    assert!((mix.gamma()[0] + 27. / 97.).abs() < 1e-14);
    assert!((next[0].re - 1.).abs() < 1e-14);
    assert!((next[1].re - 93. / 97.).abs() < 1e-14);
    x = next.clone();
    let f = [c(2. - 2. * x[0].re), c(3. - 3. * x[1].re)];
    mix.update(&Serial, &mut metric, &x, &f, &mut next).unwrap();
    assert!((next[0].re - 1.).abs() < 1e-14);
    assert!((next[1].re - 1.).abs() < 1e-14);
}
#[test]
fn ring_eviction_keeps_only_last_two_secants() {
    let mut mix = Broyden::new(
        Profile {
            history: 2,
            beta: 0.5,
            ..Default::default()
        },
        &[1.; 3],
        &[true; 3],
        1.,
        2.,
    )
    .unwrap();
    let mut metric = metric(&[1.; 3], &[true; 3]);
    let mut f = vec![c(3.), c(4.), c(5.)];
    let mut out = f.clone();
    let mut fs = Vec::new();
    let mut xs = Vec::new();
    for k in 0..7 {
        if k > 0 {
            f[(k - 1) % 3] -= c(1.);
        }
        let x = vec![c(k as f64), c((k * k) as f64), c((k * k * k) as f64)];
        let d = mix.update(&Serial, &mut metric, &x, &f, &mut out).unwrap();
        assert_eq!(d.occupancy, k.min(2));
        assert_eq!(d.inserted_slot, (k > 0).then(|| (k - 1) % 2));
        // Differences are distinct unit vectors, so A=I and gamma is exactly
        // the corresponding component of the current residual. No matrix code.
        let mut expected: Vec<_> = x.iter().zip(&f).map(|(x, f)| x + 0.5 * f).collect();
        fs.push(f.clone());
        xs.push(x.clone());
        for t in k.saturating_sub(1).max(1)..=k {
            let axis = (t - 1) % 3;
            let gamma = f[axis].re;
            for j in 0..3 {
                expected[j] -= gamma * (xs[t - 1][j] - xs[t][j] + 0.5 * (fs[t - 1][j] - fs[t][j]));
            }
        }
        for (a, b) in out.iter().zip(&expected) {
            assert!((*a - *b).norm() < 1e-12, "k={k} {out:?} {expected:?}");
        }
    }
}
#[test]
fn metric_half_sphere_zero_mode_and_high_g() {
    let mut m = Metric::new(2., &[0., 1., 4., 9.], &[true, true, true, false], 1, 4).unwrap();
    let a = [c(900.), C::new(1., 2.), C::new(2., -1.), c(40.)];
    let b = [c(-500.), C::new(3., -1.), C::new(-1., 2.), c(70.)];
    let mut out = [0.; 2];
    m.batch(&Serial, &[(&a, &b), (&a, &a)], &mut out).unwrap();
    // Re(conj(a)b): 1 and -4; 1 + (-4)/4 = 0.
    assert!(out[0].abs() < 1e-13);
    assert!((out[1] - 100. * PI).abs() < 1e-12);
}
#[test]
fn high_g_simple_mixing_and_tf_is_opt_in() {
    let p = Profile {
        preconditioning: Preconditioning::ThomasFermi,
        beta: 0.7,
        ..Default::default()
    };
    let mut mixer = Broyden::new(p, &[0., 1., 9.], &[true, true, false], 10., 2.).unwrap();
    let mut m = metric(&[0., 1., 9.], &[true, true, false]);
    let mut out = [c(0.); 3];
    mixer
        .update(
            &Serial,
            &mut m,
            &[c(2.), c(3.), c(4.)],
            &[c(0.), c(1.), c(2.)],
            &mut out,
        )
        .unwrap();
    let rs = (30. / (8. * PI)).cbrt();
    let a = (12. / PI).powf(2. / 3.) / rs;
    assert_eq!(out[0], c(2.));
    assert!((out[1].re - (3. + 0.7 / (1. + a))).abs() < 1e-14);
    assert_eq!(out[2], c(5.4));
}
#[test]
fn singular_history_is_an_error_not_a_reset() {
    let mut mix = Broyden::new(Profile::default(), &[1.], &[true], 1., 2.).unwrap();
    let mut m = metric(&[1.], &[true]);
    let mut out = [c(0.)];
    mix.update(&Serial, &mut m, &[c(0.)], &[c(1.)], &mut out)
        .unwrap();
    assert!(matches!(
        mix.update(&Serial, &mut m, &[c(1.)], &[c(1.)], &mut out),
        Err(Error::Singular { iteration: 2, .. })
    ));
}
#[test]
fn threshold_units_reset_and_floor() {
    assert_eq!(Profile::from_standard_ha(5e-9).unwrap().threshold_ry, 1e-8);
    assert_eq!(next_ethr_ry(1e-2, 2., 10.), 1e-2);
    assert_eq!(next_ethr_ry(1e-3, 0.01, 10.), 1e-4);
    assert_eq!(next_ethr_ry(1e-6, 0., 10.), 1e-13);
}
#[test]
fn distributed_density_round_trip_and_broyden() {
    for ranks in [1, 2, 3, 8] {
        let answers = dft_mpi::simulate(ranks, |comm| {
            let cell = dft_geometry::Cell::from_user_rows_bohr([
                [TAU, 0., 0.],
                [0., TAU, 0.],
                [0., 0., TAU],
            ])
            .unwrap();
            let mut plan = DensityPlan::new(comm, &cell, [5; 3], 1., 1.).unwrap();
            let mut coefficients: Vec<_> = plan
                .miller
                .iter()
                .map(|&m| match m {
                    [0, 0, 0] => c(2.),
                    [1, 0, 0] => c(0.05),
                    [0, 1, 0] => C::new(0., -0.1),
                    _ => c(0.),
                })
                .collect();
            let mut planes = vec![0.; plan.plane_len()];
            plan.inverse(comm, &coefficients, &mut planes).unwrap();
            let real = Array3::from_shape_fn([5; 3], |(i, j, _)| {
                2. + 0.1 * (TAU * i as f64 / 5.).cos() + 0.2 * (TAU * j as f64 / 5.).sin()
            });
            let mut analytic = vec![0.; planes.len()];
            plan.planes_from_full(&real, &mut analytic).unwrap();
            for (a, b) in planes.iter().zip(&analytic) {
                assert!((a - b).abs() < 2e-15);
            }
            let expected = coefficients.clone();
            plan.forward(comm, &planes, &mut coefficients).unwrap();
            for (a, b) in coefficients.iter().zip(&expected) {
                assert!((*a - *b).norm() < 1e-15);
            }
            let mut metric = Metric::new(
                cell.volume_bohr3(),
                &plan.g2,
                &plan.smooth,
                comm.size(),
                350,
            )
            .unwrap();
            let mut mix = Broyden::new(
                Profile {
                    history: 2,
                    ..Default::default()
                },
                &plan.g2,
                &plan.smooth,
                cell.volume_bohr3(),
                2.,
            )
            .unwrap();
            let mut next = coefficients.clone();
            let f: Vec<_> = plan
                .miller
                .iter()
                .map(|m| if *m == [1, 0, 0] { c(0.1) } else { c(0.) })
                .collect();
            mix.update(comm, &mut metric, &coefficients, &f, &mut next)
                .unwrap();
            let dr2 = metric.norm(comm, &f).unwrap();
            assert!((dr2 - 8. * PI * cell.volume_bohr3() * 0.01).abs() < 1e-12);
            (plan.miller.len(), dr2)
        });
        assert_eq!(answers.iter().map(|r| r.0).sum::<usize>(), 10);
        assert!(answers.iter().all(|r| r.1 == answers[0].1));
    }
}
fn tiny_file(big: bool) -> Vec<u8> {
    let mut out = Vec::new();
    let i = |v: i32| {
        if big {
            v.to_be_bytes()
        } else {
            v.to_le_bytes()
        }
    };
    let d = |v: f64| {
        if big {
            v.to_be_bytes()
        } else {
            v.to_le_bytes()
        }
    };
    let mut record = |p: Vec<u8>| {
        out.extend(i(p.len() as i32));
        out.extend(&p);
        out.extend(i(p.len() as i32));
    };
    record([1, 2, 1].into_iter().flat_map(i).collect());
    record(
        [1., 0., 0., 0., 2., 0., 0., 0., 3.]
            .into_iter()
            .flat_map(d)
            .collect(),
    );
    record([1, 0, 0, 0, 0, 0].into_iter().flat_map(i).collect());
    record([0.1, -0.2, 0.5, 0.].into_iter().flat_map(d).collect());
    out
}
#[test]
fn four_record_reader_both_endians_and_miller_mapping() {
    for big in [false, true] {
        let bytes = tiny_file(big);
        let q = QeDensity::decode(&bytes).unwrap();
        assert_eq!(q.reciprocal[1][1], 2.);
        assert_eq!(
            q.aligned(&[[0; 3], [1, 0, 0]]).unwrap(),
            [c(0.5), C::new(0.1, -0.2)]
        );
        assert_eq!(q.order(), [[0; 3], [1, 0, 0]]);
        for n in 0..bytes.len() {
            assert!(QeDensity::decode(&bytes[..n]).is_err(), "truncation {n}");
        }
        let mut bad = bytes.clone();
        *bad.last_mut().unwrap() ^= 1;
        assert!(QeDensity::decode(&bad).is_err());
    }
    assert!(half([0, 0, 0]));
    assert!(half([1, -2, -3]));
    assert!(!half([0, -1, 3]));
    let p = printed_accuracy("estimated scf accuracy < 0.00000231 Ry\n").unwrap();
    assert_eq!(p.ry, 2.31e-6);
    assert!(p.upper_bound);
    assert_eq!(p.rounding_ry, 5e-9);
}
#[test]
fn qe_si8_committed_trajectory() {
    use dft_engine::refscf::{density::g_squared, qe_reference::difference};
    let root =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/reference/qe/i25-si8");
    if !root.exists() {
        eprintln!("SKIP QE si8 reference: fixture directory absent; no synthetic substitute");
        return;
    }
    // Once a delivery exists, missing/corrupt individual files must fail.
    let dir = root.join("si8-np1");
    let read = |name: &str| {
        QeDensity::decode(&std::fs::read(dir.join(name).join("charge-density.dat")).unwrap())
            .unwrap()
    };
    let first = read("k1-out");
    let order = first.order();
    assert_eq!(order.len(), 11936);
    let cell = dft_geometry::Cell::from_user_rows_bohr([
        [10.261212856717933, 0., 0.],
        [0., 10.261212856717933, 0.],
        [0., 0., 10.261212856717933],
    ])
    .unwrap();
    first
        .validate_cell(cell.reciprocal_columns_inverse_bohr())
        .unwrap();
    let align = |q: QeDensity| {
        q.validate_cell(first.reciprocal).unwrap();
        q.aligned(&order).unwrap()
    };
    let y1 = align(first.clone());
    let mix1 = align(read("k1-mix"));
    let recover = |m: &[C], beta: f64| {
        m.iter()
            .zip(&y1)
            .map(|(m, y)| (*m - beta * y) / (1. - beta))
            .collect::<Vec<_>>()
    };
    let mut x = recover(&mix1, 0.7);
    let recovered03 = recover(&align(read("k1-mix-beta0.3")), 0.3);
    let all: Vec<_> = (0..order.len()).collect();
    let check = |label: &str, a: &[C], b: &[C], idx: &[usize]| {
        let d = difference(a, b, idx);
        println!(
            "QE_DIFF label={label} count={} max={:.17e} rms={:.17e} relative_max={:.17e} relative_rms={:.17e}",
            idx.len(),
            d.max,
            d.rms,
            d.relative_max,
            d.relative_rms
        );
        // Measured QE 7.6 Si8: max 3.47e-18; shell relative RMS < 9.70e-16.
        assert!(d.max <= 1e-17 + 1e-12 * d.scale, "{label}: {d:?}");
        assert!(d.rms <= 1e-17 + 1e-12 * d.scale, "{label} RMS: {d:?}");
        d
    };
    check("startup_beta03", &x, &recovered03, &all);
    check("startup_beta1", &y1, &align(read("k1-mix-beta1.0")), &all);
    let zero = order.iter().position(|m| *m == [0; 3]).unwrap();
    let volume = cell.volume_bohr3();
    assert!((x[zero].re * volume - 32.).abs() < 1e-12);
    let g2: Vec<_> = order
        .iter()
        .map(|&m| g_squared(first.reciprocal, m))
        .collect();
    let mut shells = Vec::<f64>::new();
    for &g in &g2 {
        if g > 0.
            && shells
                .last()
                .is_none_or(|p| (g - p).abs() > 1e-10 * g.max(1.))
        {
            shells.push(g);
            if shells.len() == 10 {
                break;
            }
        }
    }
    assert_eq!(shells.len(), 10);
    let mut metric = Metric::new(volume, &g2, &vec![true; g2.len()], 1, 350).unwrap();
    let mut mixer =
        Broyden::new(Profile::default(), &g2, &vec![true; g2.len()], volume, 32.).unwrap();
    let mut next = vec![C::default(); g2.len()];
    for k in 1..=3 {
        let y = align(read(&format!("k{k}-out")));
        let expected = align(read(&format!("k{k}-mix")));
        let f: Vec<_> = y.iter().zip(&x).map(|(y, x)| y - x).collect();
        let dr2 = metric.norm(&Serial, &f).unwrap();
        let stdout = std::fs::read_to_string(dir.join(format!("k{k}-mix/stdout.txt"))).unwrap();
        let printed =
            printed_accuracy(&stdout).expect("fixture must have a printed estimated scf accuracy");
        // QE prints '<', but this fixture specifically establishes agreement
        // with the rounded number, not merely an upper bound.
        assert!(
            (dr2 - printed.ry).abs() <= printed.rounding_ry + 1e-14,
            "k={k} dr2={dr2:.17e} printed={printed:?}"
        );
        let d = mixer
            .update(&Serial, &mut metric, &x, &f, &mut next)
            .unwrap();
        assert_eq!(d.occupancy, k - 1);
        check(&format!("k{k}_all"), &next, &expected, &all);
        check(&format!("k{k}_g0"), &next, &expected, &[zero]);
        for (i, &g) in shells.iter().enumerate() {
            let idx: Vec<_> = g2
                .iter()
                .enumerate()
                .filter_map(|(j, &v)| ((v - g).abs() <= 1e-10 * g.max(1.)).then_some(j))
                .collect();
            check(&format!("k{k}_shell{}", i + 1), &next, &expected, &idx);
        }
        let charge = next[zero].re * volume;
        let delta_charge = (next[zero] - expected[zero]).re * volume;
        assert!((charge - 32.).abs() < 1e-12);
        assert!(delta_charge.abs() < 1e-12);
        println!(
            "QE_STEP k={k} dr2_ry={dr2:.17e} printed_ry={:.17e} dr2_delta_ry={:.17e} charge_e={charge:.17e} charge_delta_e={delta_charge:.17e} condition={:?}",
            printed.ry,
            dr2 - printed.ry,
            d.condition
        );
        x = expected; // QE's own input, never a propagated BAANDI replay error.
    }
}

#[allow(dead_code)]
mod common;
#[test]
fn driver_owned_maps_and_real_verification() {
    use dft_engine::{
        GammaReal, ScfOptions, SpeciesPotential, StartDensity, Verification, refscf::driver,
    };
    use std::sync::Arc;
    // A broad Gaussian is resolved at this deliberately small test cutoff.
    // A narrow rc=0.7 fixture produces negative Fourier ringing after mixing
    // and is correctly rejected; no clipping or tolerance relaxation is used.
    let pp = common::potential_with_atomic_start(4., 1.8, true, 161, "SLA PZ");
    let cell =
        dft_geometry::Cell::from_user_rows_bohr([[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]])
            .unwrap();
    let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
    let run = |verify, max_iterations, threshold_ry| {
        let species = [SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let options = ScfOptions {
            start_density: StartDensity::Atomic,
            gamma_real: GammaReal::On,
            cutoff_hartree: 2.,
            bands: Some(3),
            max_iterations,
            verification: verify,
            ..Default::default()
        };
        // Large stop threshold is a control-flow fixture, never a production
        // convergence tolerance. It forces the verification boundary in one map.
        driver::run(
            dft_engine::refstart::StartInput {
                cell: &cell,
                atoms: &atoms,
                potentials: &species,
                options: &options,
            },
            Profile {
                threshold_ry,
                ..Default::default()
            },
        )
        .unwrap()
    };
    let off = run(Verification::Off, 1, 1e6);
    assert!(off.reference_stopped && off.outcome.converged);
    let on = run(Verification::FixedPoint, 1, 1e6);
    assert!(on.reference_stopped && on.outcome.converged);
    assert!(on.outcome.result.verification.as_ref().unwrap().passed());
    assert!(
        on.outcome.result.history[0]
            .fixed_point_certificate
            .is_some()
    );
    let values = dft_mpi::simulate(2, |comm| {
        dft_engine::parallel::set_communicator(Some(Arc::new(comm.clone())));
        let run = run(Verification::Off, 1, 1e6);
        dft_engine::parallel::set_communicator(None);
        assert_eq!(run.records.len(), 1);
        (run.outcome.result.total_hartree(), run.records[0].dr2_ry)
    });
    for (e, dr2) in values {
        assert!((e - off.outcome.result.total_hartree()).abs() < 1e-10);
        assert!((dr2 - off.records[0].dr2_ry).abs() < 1e-10);
    }
    let serial = run(Verification::Off, 2, 1e-20);
    let values = dft_mpi::simulate(2, |comm| {
        dft_engine::parallel::set_communicator(Some(Arc::new(comm.clone())));
        let run = run(Verification::Off, 2, 1e-20);
        dft_engine::parallel::set_communicator(None);
        assert_eq!(run.records.len(), 2);
        assert_eq!(run.records[0].mixing.occupancy, 0);
        assert_eq!(run.records[1].mixing.occupancy, 1);
        assert!(run.records.iter().all(|r| r.descf_ha.is_finite()));
        run.outcome.result.total_hartree()
    });
    for e in values {
        assert!((e - serial.outcome.result.total_hartree()).abs() < 1e-10);
    }
}
