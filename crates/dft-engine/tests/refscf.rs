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
    let config =
        dft_numerics::switches::RunConfig::from_pairs(&[("refscf_verify_energy", "off")]).unwrap();
    let _scope = config.enter();
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

fn stop_to_verified_clocks(run: &dft_engine::refscf::driver::Run) -> (f64, f64, f64, f64) {
    let attempts = &run.outcome.result.verification.as_ref().unwrap().attempts;
    let first = run.first_reference_stop.unwrap();
    let last = attempts.last().unwrap();
    assert!(last.passed);
    let start = run.records[first - 1].estimator_s;
    let end = run.records[last.iteration - 1].verification_end_s.unwrap();
    let (mut verify, mut recovery, mut previous_end) = (0., 0., 0.);
    let mut failed = false;
    for r in &run.records {
        assert!(previous_end <= r.map_start_s && r.map_start_s <= r.estimator_s);
        assert!(r.estimator_s <= r.map_end_s);
        previous_end = r.map_end_s;
        match (r.verification_start_s, r.verification_end_s) {
            (Some(a), Some(b)) => {
                assert!(r.estimator_s <= a && a <= b && b <= r.map_end_s);
                if failed {
                    recovery += a - r.map_start_s;
                }
                verify += b - a;
                let attempt = attempts
                    .iter()
                    .find(|a| a.iteration == r.iteration)
                    .unwrap();
                if !attempt.passed {
                    failed = true;
                    recovery += r.map_end_s - b;
                }
            }
            (None, None) => {
                if failed {
                    recovery += r.map_end_s - r.map_start_s;
                }
            }
            _ => panic!("partial attempt clock"),
        }
    }
    let interval = end - start;
    let other = interval - verify - recovery;
    assert!(other >= -64. * f64::EPSILON * end.max(1.));
    assert!((verify + recovery + other - interval).abs() <= 64. * f64::EPSILON * end.max(1.));
    (interval, verify, recovery, other)
}

fn stop_to_verified_science(run: &dft_engine::refscf::driver::Run) -> String {
    let r = &run.outcome.result;
    let mut records = run.records.clone();
    for map in &mut records {
        map.map_start_s = 0.;
        map.estimator_s = 0.;
        map.verification_start_s = map.verification_start_s.map(|_| 0.);
        map.verification_end_s = map.verification_end_s.map(|_| 0.);
        map.map_end_s = 0.;
        map.map_s = 0.;
        map.solve_s = 0.;
        map.density_s = 0.;
        map.potential_s = 0.;
        map.energy_s = 0.;
        map.density_transform_s = 0.;
        map.metric_s = 0.;
        map.mixing_s = 0.;
        map.correction_s = 0.;
        map.verification_s = 0.;
    }
    let mut calls = run.verification_calls.clone();
    for c in records
        .iter_mut()
        .flat_map(|r| &mut r.calls)
        .chain(&mut calls)
    {
        c.inner_clocks = None;
        c.ritz_changes_ha.clear();
        c.ritz_values_ha.clear();
    }
    let mut history = r.history.clone();
    for h in &mut history {
        h.timings = Default::default();
        if let Some(c) = &mut h.fixed_point_certificate {
            c.timings = Default::default();
        }
    }
    let mut verdict = r.verification.clone();
    for a in &mut verdict.as_mut().unwrap().attempts {
        a.wall = Default::default();
    }
    format!(
        "{:?}",
        (
            run.outcome.converged,
            run.first_reference_stop,
            &r.energy,
            r.density.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            r.orbitals
                .as_slice()
                .iter()
                .map(|v| (v.re.to_bits(), v.im.to_bits()))
                .collect::<Vec<_>>(),
            &r.eigenvalues,
            &r.occupations,
            history,
            verdict,
            records,
            calls
        )
    )
}

#[test]
fn stop_to_verified_first_pass_all_diagnostics() {
    use dft_engine::{GammaReal, ScfOptions, SpeciesPotential, StartDensity, refscf::driver};
    use dft_numerics::switches::RunConfig;
    let pp = common::potential_with_atomic_start(4., 1.8, true, 161, "SLA PZ");
    let cell =
        dft_geometry::Cell::from_user_rows_bohr([[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]])
            .unwrap();
    let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        start_density: StartDensity::Atomic,
        gamma_real: GammaReal::On,
        cutoff_hartree: 2.,
        bands: Some(3),
        max_iterations: 1,
        ..Default::default()
    };
    let mut previous = None;
    for level in ["off", "on"] {
        let scope =
            RunConfig::from_pairs(&[("davidson_stages", level), ("refscf_verify_energy", "off")])
                .unwrap();
        let _guard = scope.enter();
        // Control-flow fixture only: the existing first-pass case, not a new default.
        let run = driver::run(
            dft_engine::refstart::StartInput {
                cell: &cell,
                atoms: &atoms,
                potentials: &species,
                options: &options,
            },
            Profile {
                threshold_ry: 1e6,
                ..Default::default()
            },
        )
        .unwrap();
        let (interval, verification, recovery, other) = stop_to_verified_clocks(&run);
        assert_eq!(run.first_reference_stop, Some(1));
        assert_eq!(
            run.outcome
                .result
                .verification
                .as_ref()
                .unwrap()
                .attempts
                .len(),
            1
        );
        assert_eq!(recovery, 0.);
        assert!(verification > 0.);
        let science = stop_to_verified_science(&run);
        if let Some(p) = previous {
            assert_eq!(science, p);
        }
        previous = Some(science);
        println!(
            "STOP_TO_VERIFIED first_pass diagnostics={level} first=1 passing=1 attempts=1 interval_s={interval:.9} verification_s={verification:.9} recovery_s={recovery:.9} unassigned_s={other:.9}"
        );
    }
}

/// dev08 only: DFT_SI_UPF=/absolute/Si.pz-vbc.UPF cargo test --release
/// -p dft-engine --test refscf stop_to_verified_si8_displaced -- --ignored --nocapture
#[test]
#[ignore = "external Si UPF and full SCFs; release on dev08, never a Mac matrix"]
fn stop_to_verified_si8_displaced() {
    use dft_engine::{GammaReal, ScfOptions, SpeciesPotential, StartDensity, refscf::driver};
    use dft_numerics::switches::RunConfig;
    let pp = dft_upf::parse_upf(
        &std::fs::read_to_string(std::env::var("DFT_SI_UPF").expect("DFT_SI_UPF required"))
            .unwrap(),
    )
    .unwrap();
    let a = 10.261_212_856_717_933;
    let cell =
        dft_geometry::Cell::from_user_rows_bohr([[a, 0., 0.], [0., a, 0.], [0., 0., a]]).unwrap();
    let base = [
        [0., 0., 0.],
        [0., 0.5, 0.5],
        [0.5, 0., 0.5],
        [0.5, 0.5, 0.],
        [0.25, 0.25, 0.25],
        [0.25, 0.75, 0.75],
        [0.75, 0.25, 0.75],
        [0.75, 0.75, 0.25],
    ];
    let atoms: Vec<_> = base
        .into_iter()
        .enumerate()
        .map(|(i, p)| {
            dft_geometry::Atom::new(
                "Si",
                [
                    p[0] + 1e-4 * ((i * 7 + 1) as f64).sin(),
                    p[1] + 1e-4 * ((i * 11 + 3) as f64).sin(),
                    p[2] + 1e-4 * ((i * 13 + 5) as f64).sin(),
                ],
            )
        })
        .collect();
    let species = [SpeciesPotential {
        label: "Si",
        potential: &pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: 15.,
        bands: Some(20),
        start_density: StartDensity::Atomic,
        gamma_real: GammaReal::On,
        max_iterations: 200,
        ..Default::default()
    };
    let mut science = None;
    let mut baseline = None;
    for (variant, level, pairs) in [
        (
            "off",
            "off",
            vec![
                ("refscf_verify_reuse", "off"),
                ("refscf_verify_gate", "off"),
            ],
        ),
        (
            "off",
            "on",
            vec![
                ("refscf_verify_reuse", "off"),
                ("refscf_verify_gate", "off"),
            ],
        ),
        ("default", "off", vec![]),
        ("default", "on", vec![]),
        ("first", "off", vec![("refscf_verify_gate_first", "on")]),
        (
            "reuse",
            "off",
            vec![("refscf_verify_reuse", "on"), ("refscf_verify_gate", "off")],
        ),
        ("gate", "off", vec![("refscf_verify_gate", "on")]),
        (
            "both",
            "off",
            vec![("refscf_verify_reuse", "on"), ("refscf_verify_gate", "on")],
        ),
    ] {
        let mut pairs = pairs;
        pairs.push(("davidson_stages", level));
        pairs.push(("refscf_verify_energy", "off"));
        let scope = RunConfig::from_pairs(&pairs).unwrap();
        let reuse = if scope.on(dft_numerics::switches::Id::REFSCF_VERIFY_REUSE) {
            "on"
        } else {
            "off"
        };
        let gate = if scope.on(dft_numerics::switches::Id::REFSCF_VERIFY_GATE) {
            "on"
        } else {
            "off"
        };
        let _guard = scope.enter();
        let run = driver::run(
            dft_engine::refstart::StartInput {
                cell: &cell,
                atoms: &atoms,
                potentials: &species,
                options: &options,
            },
            Profile::default(),
        )
        .unwrap();
        let attempts = &run.outcome.result.verification.as_ref().unwrap().attempts;
        verify_cost_contract(&run);
        let energy = run.outcome.result.total_hartree();
        let base = *baseline.get_or_insert(energy);
        assert!((energy - base).abs() < 5e-9);
        if reuse == "off" && gate == "off" {
            assert_eq!(run.first_reference_stop, Some(10));
            assert_eq!(
                attempts
                    .iter()
                    .map(|a| (a.iteration, a.passed))
                    .collect::<Vec<_>>(),
                [(10, false), (11, false), (12, true)]
            );
            let value = stop_to_verified_science(&run);
            if let Some(s) = &science {
                assert_eq!(&value, s);
            }
            science = Some(value);
        }
        let reused = run
            .records
            .iter()
            .filter(|r| {
                r.verification_cost
                    .as_ref()
                    .is_some_and(|c| c.reused_from.is_some())
            })
            .count();
        let skipped = run
            .records
            .iter()
            .filter(|r| r.verification_cost.as_ref().is_some_and(|c| c.skipped))
            .count();
        let solves =
            run.records.iter().map(|r| r.calls.len()).sum::<usize>() + run.verification_calls.len();
        let ordinary_requests = run.records.len() - reused
            + run
                .records
                .iter()
                .filter(|r| r.first_iteration_redo)
                .count();
        let map_requests = ordinary_requests + attempts.len();
        if reuse == "on" {
            assert!(reused > 0);
        }
        let (interval, verification, recovery, other) = stop_to_verified_clocks(&run);
        assert!(recovery > 0.);
        println!(
            "STOP_TO_VERIFIED si8 variant={variant} diagnostics={level} reuse={reuse} gate={gate} first={:?} passing={} attempts={} evaluated_maps={solves} map_requests={map_requests} solves={solves} reused={reused} skipped={skipped} energy_ha={energy:.17e} delta_ha={:.9e} interval_s={interval:.9} verification_s={verification:.9} recovery_s={recovery:.9} unassigned_s={other:.9}",
            run.first_reference_stop,
            run.records.len(),
            attempts.len(),
            energy - base
        );
    }
}

#[test]
fn verification_cost_predictor_fallback_and_refresh() {
    use dft_engine::refscf::driver::VerificationPredictor;
    let mut p = VerificationPredictor::default();
    assert_eq!(p.predict(1e-9), None);
    assert!(!p.skip(None, 5e-9));
    p.observe(1e-8, 1e-7);
    assert_eq!(p.predict(5e-9), Some(5e-8));
    assert!(p.skip(p.predict(5e-9), 5e-9));
    assert!(!p.skip(p.predict(4e-9), 5e-9)); // refresh after one skip
    p.observe(4e-9, 2e-8);
    assert!(!p.skip(p.predict(1e-10), 5e-9));
    for invalid in [0., -1., f64::NAN, f64::INFINITY, 4e-9, 5e-9] {
        assert_eq!(p.predict(invalid), None);
        assert!(!p.skip(p.predict(invalid), 5e-9));
    }
    for invalid in [0., -1., f64::NAN, f64::INFINITY] {
        p.observe(invalid, 1e-7);
        assert_eq!(p.predict(1e-9), None);
    }
    p.observe(1e-8, f64::NAN);
    assert_eq!(p.predict(1e-9), None);
    assert!(!p.skip(Some(f64::NAN), 5e-9));
    assert!(p.skip(Some(5e-9), 5e-9)); // equality cannot pass the strict '<'
}

#[test]
fn first_prediction_linear_oracle_and_read_only_history() {
    use dft_engine::refscf::broyden::PredictionScratch;
    use dft_engine::refscf::driver::contraction_prediction;
    let profile = Profile {
        history: 2,
        beta: 0.7,
        ..Default::default()
    };
    let new = || Broyden::new(profile, &[1.; 2], &[true; 2], 1., 2.).unwrap();
    let mut mix = new();
    let mut control = new();
    let mut metric = Metric::new(1. / (8. * PI), &[1.; 2], &[true; 2], 1, 350).unwrap();
    let mut scratch = PredictionScratch::new(2);
    // F(x)=b-2x amplifies every direction by -2. Its residual f=b-3x
    // gives an independent exact oracle, including an initially failing
    // candidate followed by the mandatory next-candidate attempt.
    let residual = |x: &[C]| [c(2.) - 3. * x[0], c(3.) - 3. * x[1]];
    let mut out = [c(0.); 2];
    let mut expected = out;
    assert_eq!(
        mix.predict_output_error(&Serial, &mut metric, &out, &residual(&out), &mut scratch)
            .unwrap(),
        None
    );
    for x in [[c(0.), c(0.)], [c(1.), c(0.)]] {
        mix.update(&Serial, &mut metric, &x, &residual(&x), &mut out)
            .unwrap();
        control
            .update(&Serial, &mut metric, &x, &residual(&x), &mut expected)
            .unwrap();
    }
    let x = [c(1.), c(1.)];
    let f = residual(&x);
    let actual = residual(&[x[0] + f[0], x[1] + f[1]]);
    let actual_ha = metric.norm(&Serial, &actual).unwrap() / 2.;
    let gamma = mix.gamma().to_vec();
    let predicted = mix
        .predict_output_error(&Serial, &mut metric, &x, &f, &mut scratch)
        .unwrap()
        .unwrap();
    assert!((predicted - actual_ha).abs() < 1e-14);
    let threshold = actual_ha / 2.;
    assert!(metric.norm(&Serial, &f).unwrap() / 2. < threshold);
    let mut gate = dft_engine::refscf::driver::VerificationPredictor::default();
    assert!(gate.skip(Some(predicted), threshold));
    assert!(!gate.skip(Some(predicted), threshold));
    assert_eq!(mix.gamma(), gamma);
    // A zero density secant makes the observed gain undefined.
    assert_eq!(
        mix.predict_output_error(&Serial, &mut metric, &[c(1.), c(0.)], &f, &mut scratch)
            .unwrap(),
        None
    );
    assert_eq!(
        mix.predict_output_error(&Serial, &mut metric, &x, &[c(f64::NAN); 2], &mut scratch)
            .unwrap(),
        None
    );
    mix.update(&Serial, &mut metric, &x, &f, &mut out).unwrap();
    control
        .update(&Serial, &mut metric, &x, &f, &mut expected)
        .unwrap();
    assert_eq!(out, expected);
    assert_eq!(mix.gamma(), control.gamma());
    assert_eq!(mix.occupancy(), control.occupancy());
    println!("FIRST_PREDICTION linear predicted_ha={predicted:.17e} actual_ha={actual_ha:.17e}");
    assert_eq!(contraction_prediction(1e-7, 1e-8), Some(5e-10));
    for bad in [0., -1., f64::NAN, f64::INFINITY, 1e-8] {
        assert_eq!(contraction_prediction(bad, 1e-8), None);
        assert_eq!(contraction_prediction(1e-8, bad), None);
    }
}

#[test]
fn verification_cost_small_fixture() {
    verification_cost_fixture(false);
}

#[test]
fn first_prediction_eviction_and_preconditioner_independence() {
    use dft_engine::refscf::broyden::PredictionScratch;
    for preconditioning in [Preconditioning::Plain, Preconditioning::ThomasFermi] {
        let profile = Profile {
            history: 1,
            preconditioning,
            ..Default::default()
        };
        let mut mixer = Broyden::new(profile, &[1.], &[true], 1., 2.).unwrap();
        let mut metric = metric(&[1.], &[true]);
        let mut scratch = PredictionScratch::new(1);
        let mut out = [c(0.)];
        mixer
            .update(&Serial, &mut metric, &[c(0.)], &[c(4.)], &mut out)
            .unwrap();
        // Delta x=-1, Delta f=3, hence squared gain 4 and energy 4/2.
        assert_eq!(
            mixer
                .predict_output_error(&Serial, &mut metric, &[c(1.)], &[c(1.)], &mut scratch)
                .unwrap(),
            Some(2.)
        );
        mixer
            .update(&Serial, &mut metric, &[c(1.)], &[c(1.)], &mut out)
            .unwrap();
        // The new secant replaces that gain: Delta x=-1, Delta f=0.5,
        // squared gain 0.25, current norm 0.25, predicted Ha=0.03125.
        assert_eq!(
            mixer
                .predict_output_error(&Serial, &mut metric, &[c(2.)], &[c(0.5)], &mut scratch)
                .unwrap(),
            Some(0.03125)
        );
        assert_eq!(mixer.occupancy(), 1);
    }
}

#[test]
fn verification_cost_failed_attempt_fixture() {
    verification_cost_fixture(true);
}

fn verification_cost_fixture(active: bool) {
    use dft_engine::{GammaReal, ScfOptions, SpeciesPotential, StartDensity, refscf::driver};
    use dft_numerics::switches::{Id, RunConfig};
    assert!(
        !RunConfig::from_pairs(&[])
            .unwrap()
            .on(Id::REFSCF_VERIFY_REUSE)
    );
    // Default on since 2026-11-24 (dev08 A/B); each arm below sets it explicitly.
    assert!(
        RunConfig::from_pairs(&[])
            .unwrap()
            .on(Id::REFSCF_VERIFY_GATE)
    );
    let (z, rc, length, bands) = if active {
        (4., 1.2, 11., 3)
    } else {
        (2., 1.8, 7., 2)
    };
    let profile = Profile {
        history: if active { 4 } else { 8 },
        ..Default::default()
    };
    let pp = common::potential_with_atomic_start(z, rc, active, 161, "SLA PZ");
    let cell = dft_geometry::Cell::from_user_rows_bohr([
        [length, 0., 0.],
        [0.2, length, 0.],
        [0., 0.1, length],
    ])
    .unwrap();
    let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        start_density: StartDensity::Atomic,
        gamma_real: GammaReal::On,
        cutoff_hartree: 2.,
        bands: Some(bands),
        max_iterations: 100,
        ..Default::default()
    };
    let mut baseline = None;
    let mut fingerprint = None;
    let mut off_fingerprint = None;
    for (name, pairs) in [
        ("default", vec![]),
        ("first_off", vec![("refscf_verify_gate_first", "off")]),
        ("first", vec![("refscf_verify_gate_first", "on")]),
        (
            "off",
            vec![
                ("refscf_verify_reuse", "off"),
                ("refscf_verify_gate", "off"),
            ],
        ),
        (
            "first_inert",
            vec![
                ("refscf_verify_gate", "off"),
                ("refscf_verify_gate_first", "on"),
            ],
        ),
        (
            "reuse",
            vec![("refscf_verify_reuse", "on"), ("refscf_verify_gate", "off")],
        ),
        ("gate", vec![("refscf_verify_gate", "on")]),
        (
            "both",
            vec![("refscf_verify_reuse", "on"), ("refscf_verify_gate", "on")],
        ),
    ] {
        let mut pairs = pairs;
        pairs.push(("refscf_verify_energy", "off"));
        let config = RunConfig::from_pairs(&pairs).unwrap();
        let _guard = config.enter();
        let run = driver::run(
            dft_engine::refstart::StartInput {
                cell: &cell,
                atoms: &atoms,
                potentials: &species,
                options: &options,
            },
            profile,
        )
        .unwrap();
        assert!(run.outcome.converged, "{name}: {}", run.exit_reason);
        let energy = run.outcome.result.total_hartree();
        let base = *baseline.get_or_insert(energy);
        assert!(
            (energy - base).abs() < 5e-9,
            "{name}: delta={} Ha",
            energy - base
        );
        if name == "default" {
            fingerprint = Some(stop_to_verified_science(&run));
        }
        // The default equals the gate-on arm since the 2026-11-24 adoption.
        if name == "gate" || name == "first_off" {
            assert_eq!(
                fingerprint.as_ref().unwrap(),
                &stop_to_verified_science(&run)
            );
        }
        if name == "off" {
            off_fingerprint = Some(stop_to_verified_science(&run));
        }
        if name == "first_inert" {
            assert_eq!(
                off_fingerprint.as_ref().unwrap(),
                &stop_to_verified_science(&run)
            );
        }
        assert!(
            !RunConfig::from_pairs(&[])
                .unwrap()
                .on(Id::REFSCF_VERIFY_GATE_FIRST)
        );
        verify_cost_contract(&run);
        let attempts = &run.outcome.result.verification.as_ref().unwrap().attempts;
        let reused = run
            .records
            .iter()
            .filter(|r| {
                r.verification_cost
                    .as_ref()
                    .is_some_and(|c| c.reused_from.is_some())
            })
            .count();
        let skipped = run
            .records
            .iter()
            .filter(|r| r.verification_cost.as_ref().is_some_and(|c| c.skipped))
            .count();
        let solves: usize =
            run.records.iter().map(|r| r.calls.len()).sum::<usize>() + run.verification_calls.len();
        let expected = if !active {
            (6, 7, 1, 0, 0)
        } else {
            match name {
                "reuse" | "both" => (29, 30, 2, 1, 0),
                "gate" | "default" | "first_off" => (28, 30, 2, 0, 1),
                "first" => (27, 28, 1, 0, 1),
                _ => (27, 29, 2, 0, 0),
            }
        };
        assert_eq!(
            (run.records.len(), solves, attempts.len(), reused, skipped),
            expected
        );
        if active && name == "first" {
            assert_eq!(run.first_reference_stop, Some(23));
            assert_eq!(attempts[0].iteration, 27);
            assert!(run.records[22].verification_cost.as_ref().unwrap().skipped);
            assert!(
                run.records[26]
                    .verification_cost
                    .as_ref()
                    .unwrap()
                    .predicted_error_ha
                    .unwrap()
                    >= 5e-9
            );
        }
        for r in run
            .records
            .iter()
            .filter(|r| r.dr2_ry < profile.threshold_ry)
        {
            if let Some(cost) = &r.verification_cost {
                println!(
                    "FIRST_PREDICTION active={active} variant={name} iteration={} dr2_ry={:.9e} prediction_ha={:?} secant_ha={:?} contraction_ha={:?} actual_ha={:?} skipped={}",
                    r.iteration,
                    r.dr2_ry,
                    cost.predicted_error_ha,
                    cost.first_prediction_ha,
                    cost.contraction_prediction_ha,
                    cost.actual_error_ha,
                    cost.skipped
                );
            }
        }
        if active && (name == "gate" || name == "reuse" || (name == "first" && skipped > 0)) {
            let max_iterations = if name == "gate" || name == "first" {
                run.records
                    .iter()
                    .find(|r| r.verification_cost.as_ref().unwrap().skipped)
                    .unwrap()
                    .iteration
            } else {
                attempts[0].iteration
            };
            let limited = ScfOptions {
                max_iterations,
                ..options.clone()
            };
            let exhausted = driver::run(
                dft_engine::refstart::StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &species,
                    options: &limited,
                },
                profile,
            )
            .unwrap();
            assert!(!exhausted.outcome.converged);
            assert!(
                !exhausted
                    .outcome
                    .result
                    .verification
                    .as_ref()
                    .unwrap()
                    .passed()
            );
            assert_eq!(exhausted.records.len(), max_iterations);
            if name == "gate" || name == "first" {
                assert!(
                    exhausted
                        .records
                        .last()
                        .unwrap()
                        .verification_cost
                        .as_ref()
                        .unwrap()
                        .skipped
                );
            }
        }
        if active && name == "reuse" {
            assert!(reused > 0);
        }
        if active && name == "gate" {
            assert!(skipped > 0);
        }
        println!(
            "VERIFY_COST active={active} variant={name} energy_ha={energy:.17e} delta_ha={:.9e} iterations={} solves={solves} attempts={} reused={reused} skipped={skipped} error_ha={:.9e}",
            energy - base,
            run.records.len(),
            attempts.len(),
            attempts.last().unwrap().coulomb_error_ha
        );
    }
}

fn verify_cost_contract(run: &dft_engine::refscf::driver::Run) {
    let attempts = &run.outcome.result.verification.as_ref().unwrap().attempts;
    let pass = attempts.last().unwrap();
    assert!(pass.passed);
    assert_eq!(pass.iteration, run.records.len());
    assert_eq!(
        run.outcome.result.total_hartree().to_bits(),
        run.records.last().unwrap().energy.total_hartree().to_bits()
    );
    for attempt in attempts {
        assert_eq!(attempt.threshold_ha, 5e-9);
        assert_eq!(
            attempt.passed,
            attempt.coulomb_error_ha.is_finite()
                && attempt.coulomb_error_ha < 5e-9
                && attempt
                    .residual_bound_ha
                    .is_none_or(|b| b.is_finite() && b < 5e-9)
        );
        let map = &run.records[attempt.iteration - 1];
        assert_eq!(
            attempt.eigen_tolerance_ha,
            next_ethr_ry(
                map.ethr_ry,
                map.dr2_ry,
                run.outcome.result.occupations.iter().sum()
            ) / 2.
        );
        if let Some(c) = &map.verification_cost {
            assert_eq!(c.actual_error_ha, Some(attempt.coulomb_error_ha));
        }
    }
    let mut skipped_previous_candidate = false;
    for r in &run.records {
        if let Some(c) = &r.verification_cost {
            if r.dr2_ry < 1e-8 {
                assert!(!(skipped_previous_candidate && c.skipped));
                skipped_previous_candidate = c.skipped;
            }
            if c.skipped {
                assert!(c.gate_enabled && c.predicted_error_ha.unwrap() >= 5e-9);
                assert!(
                    c.actual_error_ha.is_none()
                        && r.verification_start_s.is_none()
                        && r.verification_end_s.is_none()
                );
                assert!(r.verification_calls.is_empty());
            }
            if let Some(source) = c.reused_from {
                assert_eq!(source + 1, r.iteration);
                let a = attempts.iter().find(|a| a.iteration == source).unwrap();
                assert!(!a.passed);
                assert_eq!(r.ethr_ry, 2. * a.eigen_tolerance_ha);
                assert!(r.calls.is_empty());
                assert_eq!(r.h_columns, 0);
                assert_eq!(r.map_s, 0.);
                assert_eq!(run.records[source - 1].descf_ha, 0.);
                assert!((r.dr2_ry / 2. - a.coulomb_error_ha).abs() < 1e-20);
            }
        }
    }
    stop_to_verified_clocks(run);
}

#[test]
fn setup_phase_clocks_cover_small_start() {
    let pp = common::potential_with_atomic_start(4., 1.8, true, 161, "SLA PZ");
    for (label, rows) in [
        ("cubic", [[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]]),
        ("skew", [[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]]),
        ("elongated", [[5., 0., 0.], [0., 5., 0.], [0., 0., 10.]]),
    ] {
        let cell = dft_geometry::Cell::from_user_rows_bohr(rows).unwrap();
        let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
        let potentials = [dft_engine::SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let options = dft_engine::ScfOptions {
            gamma_real: dft_engine::GammaReal::On,
            start_density: dft_engine::StartDensity::Atomic,
            cutoff_hartree: 2.,
            bands: Some(5),
            max_iterations: 1,
            verification: dft_engine::Verification::Off,
            ..Default::default()
        };
        let run = dft_engine::refscf::driver::run_with(
            dft_engine::refstart::StartInput {
                cell: &cell,
                atoms: &atoms,
                potentials: &potentials,
                options: &options,
            },
            Profile::default(),
            &Serial,
        )
        .unwrap();
        assert!(run.setup.total_s > 0.);
        assert!(run.setup.collective_s >= 0. && run.setup.collective_s < run.setup.total_s);
        for key in [
            "basis_gvectors",
            "field_fft_plan",
            "distribution_owned_rows",
            "nonlocal_tables",
            "atomic_density",
            "atomic_orbital_fill",
            "atomic_orthonormalize_complete",
            "atomic_diagonalize_rotate",
            "density_fft_plan",
        ] {
            assert!(run.setup.phases[key] > 0., "{key}");
        }
        assert!(run.records[0].map_start_s >= run.setup.total_s);
        assert!(run.output.phases["density_gather"] > 0.);
        assert!(run.output.phases["orbital_gather_unpack"] > 0.);
        eprintln!(
            "P3-SETUP observation {label}: {:?}; output {:?}",
            run.setup, run.output
        );
    }
}

#[test]
fn setup_phase_reuse_preserves_owned_start_bits() {
    use dft_engine::refstart::{StartInput, start};
    use dft_numerics::switches::RunConfig;
    let shapes = [
        [[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]],
        [[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]],
        [[5., 0., 0.], [0., 5., 0.], [0., 0., 10.]],
        [[8., 0., 0.], [0., 8., 0.], [0., 0., 4.]],
    ];
    for (shape, rows) in shapes.into_iter().enumerate() {
        let pp = common::potential_with_atomic_start(4., 1.8, shape % 2 == 0, 161, "SLA PZ");
        let cell = dft_geometry::Cell::from_user_rows_bohr(rows).unwrap();
        let atoms = [
            dft_geometry::Atom::new("Xx", [0., 0., 0.]),
            dft_geometry::Atom::new("Xx", [0.43, 0.51, 0.32]),
        ];
        let potentials = [dft_engine::SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        for ranks in [1, 3] {
            let build = |reuse, threads| {
                dft_mpi::simulate(ranks, |comm| {
                    let config = RunConfig::from_pairs(&[
                        ("refscf_setup_phase_reuse", reuse),
                        ("stage_timers", "off"),
                        ("phase_profile", "off"),
                        ("davidson_stages", "off"),
                    ])
                    .unwrap();
                    let _scope = config.enter();
                    let options = dft_engine::ScfOptions {
                        cutoff_hartree: 2.,
                        bands: Some(9),
                        threads,
                        gamma_real: dft_engine::GammaReal::On,
                        ..Default::default()
                    };
                    let s = start(
                        StartInput {
                            cell: &cell,
                            atoms: &atoms,
                            potentials: &potentials,
                            options: &options,
                        },
                        comm,
                    )
                    .unwrap();
                    assert!(s.setup.total_s > 0. && s.setup.phases["atomic_orbital_fill"] > 0.);
                    assert!(s.setup.collective_s >= 0. && s.setup.collective_s < s.setup.total_s);
                    let bits: Vec<_> = s
                        .orbitals
                        .as_slice()
                        .iter()
                        .chain(&s.ritz_values_hartree)
                        .chain(&s.density)
                        .chain(&s.ionic)
                        .chain(&s.potential)
                        .chain(&s.diagonal_static)
                        .map(|x| x.to_bits())
                        .chain(
                            s.density_g
                                .iter()
                                .flat_map(|z| [z.re.to_bits(), z.im.to_bits()]),
                        )
                        .collect();
                    (
                        bits,
                        s.setup.phases["atomic_orbital_fill"],
                        s.setup.phases["atomic_local_tables_structure_factors"],
                    )
                })
            };
            let off = build("off", 1);
            for threads in [1, 2] {
                let on = build("on", threads);
                for (rank, (a, b)) in off.iter().zip(&on).enumerate() {
                    assert_eq!(
                        a.0, b.0,
                        "shape={shape} ranks={ranks} rank={rank} threads={threads}"
                    );
                }
            }
        }
    }
}

#[test]
fn setup_phase_reuse_preserves_verified_science_at_both_diagnostic_levels() {
    use dft_numerics::switches::RunConfig;
    let pp = common::potential_with_atomic_start(2., 1.8, false, 161, "SLA PZ");
    let cell =
        dft_geometry::Cell::from_user_rows_bohr([[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]])
            .unwrap();
    let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
    let potentials = [dft_engine::SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = dft_engine::ScfOptions {
        cutoff_hartree: 2.,
        bands: Some(2),
        max_iterations: 100,
        gamma_real: dft_engine::GammaReal::On,
        ..Default::default()
    };
    let mut baseline = None;
    for reuse in ["off", "on"] {
        for diagnostics in ["off", "on"] {
            let config = RunConfig::from_pairs(&[
                ("refscf_setup_phase_reuse", reuse),
                ("davidson_stages", diagnostics),
            ])
            .unwrap();
            let _scope = config.enter();
            let run = dft_engine::refscf::driver::run_with(
                dft_engine::refstart::StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &options,
                },
                Profile::default(),
                &Serial,
            )
            .unwrap();
            assert!(
                run.outcome.converged && run.outcome.result.verification.as_ref().unwrap().passed()
            );
            assert!(run.setup.total_s > 0. && run.output.total_s > 0.);
            let science = stop_to_verified_science(&run);
            if let Some(previous) = &baseline {
                assert_eq!(&science, previous);
            } else {
                baseline = Some(science);
            }
            eprintln!(
                "P3-SETUP verified reuse={reuse} diagnostics={diagnostics}: E={:.17e} Ha maps={}",
                run.outcome.result.energy.total_hartree(),
                run.records.len()
            );
        }
    }
}

#[test]
fn verification_energy_attempts_and_failure_recovery() {
    use dft_engine::{
        GammaReal, ScfOptions, SpeciesPotential, StartDensity, Verification, refscf::driver,
    };
    use dft_numerics::switches::{Id, RunConfig};
    // Default off since 2026-11-25 (dev08 screen: no earlier pass, slower on si128); tested explicitly.
    assert!(
        !RunConfig::from_pairs(&[])
            .unwrap()
            .on(Id::REFSCF_VERIFY_ENERGY)
    );
    for active in [false, true] {
        let (z, rc, length, bands) = if active {
            (4., 1.2, 11., 3)
        } else {
            (2., 1.8, 7., 2)
        };
        let pp = common::potential_with_atomic_start(z, rc, active, 161, "SLA PZ");
        let cell = dft_geometry::Cell::from_user_rows_bohr([
            [length, 0., 0.],
            [0.2, length, 0.],
            [0., 0.1, length],
        ])
        .unwrap();
        let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
        let species = [SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let options = ScfOptions {
            start_density: StartDensity::Atomic,
            gamma_real: GammaReal::On,
            cutoff_hartree: 2.,
            bands: Some(bands),
            max_iterations: 100,
            ..Default::default()
        };
        let profile = Profile {
            history: if active { 4 } else { 8 },
            ..Default::default()
        };
        let run = |pairs: &[(&str, &str)], options: &ScfOptions, profile| {
            let config = RunConfig::from_pairs(pairs).unwrap();
            let _scope = config.enter();
            driver::run(
                dft_engine::refstart::StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &species,
                    options,
                },
                profile,
            )
            .unwrap()
        };
        let old = run(&[("refscf_verify_energy", "off")], &options, profile);
        let energy = run(&[("refscf_verify_energy", "on")], &options, profile);
        let gates = run(
            &[
                ("refscf_verify_energy", "on"),
                ("refscf_verify_gate", "on"),
                ("refscf_verify_gate_first", "on"),
            ],
            &options,
            profile,
        );
        assert_eq!(
            stop_to_verified_science(&energy),
            stop_to_verified_science(&gates)
        );
        assert!(energy.outcome.converged);
        let attempts = &energy
            .outcome
            .result
            .verification
            .as_ref()
            .unwrap()
            .attempts;
        assert_eq!(attempts.len(), if active { 2 } else { 1 });
        assert_eq!(Some(attempts[0].iteration), energy.first_reference_stop);
        assert_eq!(energy.records.len(), attempts.last().unwrap().iteration);
        assert!(!energy.verification_calls.is_empty());
        assert!(energy.records.iter().all(|r| r.verification_cost.is_none()));
        let last = energy.records.last().unwrap();
        let e = last.verification_energy.as_ref().unwrap();
        assert_eq!(
            e.delivered_ha.to_bits(),
            last.variational_energy_ha.to_bits()
        );
        assert_eq!(last.descf_ha, 0.);
        assert!(e.passed());
        assert_eq!(attempts[0].threshold_ha, 5e-9);
        assert!(attempts[0].residual_bound_ha.unwrap().is_finite());
        assert!(attempts[0].density_residual.is_finite());
        let delta =
            (energy.outcome.result.total_hartree() - old.outcome.result.total_hartree()).abs();
        assert!(delta < 5e-9, "delivered difference {delta} Ha");
        if active {
            assert!(
                old.outcome
                    .result
                    .verification
                    .as_ref()
                    .unwrap()
                    .attempts
                    .len()
                    > 1
            );
            assert!(attempts[0].coulomb_error_ha > 5e-9);
            assert!(!attempts[0].passed);
            assert!(
                energy.records[attempts[0].iteration - 1]
                    .verification_energy
                    .as_ref()
                    .unwrap()
                    .change_ha()
                    >= 5e-9
            );
            assert!(energy.records.len() < old.records.len());
        }
        // Verification must not replace any part of the candidate state,
        // whether the first attempt passes or fails.
        let candidate = run(
            &[],
            &ScfOptions {
                max_iterations: attempts[0].iteration,
                ..options.clone()
            },
            profile,
        );
        let unverified = run(
            &[],
            &ScfOptions {
                verification: Verification::Off,
                ..options.clone()
            },
            profile,
        );
        assert_eq!(
            format!("{:?}", candidate.outcome.result.energy),
            format!("{:?}", unverified.outcome.result.energy)
        );
        assert_eq!(
            candidate.outcome.result.density,
            unverified.outcome.result.density
        );
        assert_eq!(
            candidate.outcome.result.orbitals,
            unverified.outcome.result.orbitals
        );
        assert_eq!(
            candidate.outcome.result.eigenvalues,
            unverified.outcome.result.eigenvalues
        );
        stop_to_verified_clocks(&energy);
        println!(
            "VERIFY_ENERGY active={active} first_stop={} old_maps={} new_maps={} attempts={} delta_old_ha={delta:.17e} energy_change_ha={:.17e} coulomb_ha={:.17e}",
            attempts[0].iteration,
            old.records.len(),
            energy.records.len(),
            attempts.len(),
            e.change_ha(),
            attempts[0].coulomb_error_ha
        );
        // Force two early stops: the fixed energy threshold must still reject,
        // charge both attempts, perform recovery and leave the run unverified.
        let failed = run(
            &[("refscf_verify_energy", "on")],
            &ScfOptions {
                max_iterations: 2,
                ..options.clone()
            },
            Profile {
                threshold_ry: 1e6,
                ..profile
            },
        );
        assert!(!failed.outcome.converged);
        let attempts = &failed
            .outcome
            .result
            .verification
            .as_ref()
            .unwrap()
            .attempts;
        assert_eq!(attempts.len(), 2);
        assert!(attempts.iter().all(|a| !a.passed && a.threshold_ha == 5e-9));
        assert!(
            failed
                .records
                .iter()
                .all(|r| r.verification_energy.as_ref().is_some_and(|e| !e.passed()))
        );
        assert!(failed.records[0].mixing_s > 0.);
        assert!(
            failed
                .records
                .iter()
                .all(|r| r.verification_start_s.unwrap() >= r.estimator_s
                    && r.verification_end_s.unwrap() >= r.verification_start_s.unwrap()
                    && r.map_end_s >= r.verification_end_s.unwrap())
        );
    }
}

#[test]
fn verification_energy_threshold_is_strict_and_finite() {
    use dft_engine::refscf::record::VerificationEnergy;
    let assess = |delivered_ha, checked_ha| {
        VerificationEnergy {
            delivered_ha,
            checked_ha,
        }
        .passed()
    };
    assert!(assess(0., 4.999e-9));
    assert!(assess(0., -4.999e-9));
    assert!(!assess(0., 5e-9));
    assert!(!assess(0., -5e-9));
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(!assess(bad, 0.));
        assert!(!assess(0., bad));
        assert!(!assess(bad, bad));
    }
}

/// P3-SETUP2: compare every owned start value, including off-diagonal projector
/// couplings, mixed species, missing-band completion and J > B Ritz truncation.
#[test]
fn setup2_owned_start_identity_and_work() {
    use dft_engine::refstart::{StartInput, start};
    use dft_numerics::switches::RunConfig;
    let shapes = [
        ("cubic", [[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]]),
        ("skew", [[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]]),
        ("elongated", [[5., 0., 0.], [0., 5., 0.], [0., 0., 10.]]),
        ("slab", [[8., 0., 0.], [0., 8., 0.], [0., 0., 4.]]),
    ];
    for (shape, (label, rows)) in shapes.into_iter().enumerate() {
        let pp = common::potential_with_atomic_start(2., 1.8, shape != 3, 161, "SLA PZ");
        let pp2 = common::potential_with_atomic_start(2., 2., false, 161, "SLA PZ");
        let cell = dft_geometry::Cell::from_user_rows_bohr(rows).unwrap();
        let mut atoms = vec![
            dft_geometry::Atom::new("Xx", [0., 0., 0.]),
            dft_geometry::Atom::new("Xx", [0.43, 0.51, 0.32]),
        ];
        let mut potentials = vec![dft_engine::SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        if shape == 1 {
            atoms.push(dft_geometry::Atom::new("Yy", [0.1, 0.7, 0.8]));
            potentials.push(dft_engine::SpeciesPotential {
                label: "Yy",
                potential: &pp2,
            });
        }
        for ranks in [1, 3, 8] {
            let build = |diagonal, orbital, threads| {
                dft_mpi::simulate(ranks, |comm| {
                    let config = RunConfig::from_pairs(&[
                        // This test pins the QR start path (P3-SETUP2 phases).
                        ("refstart_generalized", "off"),
                        ("refstart_compact_g", "off"),
                        ("refstart_diagonal_reuse", diagonal),
                        ("refstart_orbital_reuse", orbital),
                        (
                            "refscf_setup_phase_reuse",
                            if shape % 2 == 0 { "on" } else { "off" },
                        ),
                    ])
                    .unwrap();
                    let _scope = config.enter();
                    let options = dft_engine::ScfOptions {
                        cutoff_hartree: 2.,
                        bands: Some([6, 13, 8, 9][shape]),
                        threads,
                        gamma_real: dft_engine::GammaReal::On,
                        ..Default::default()
                    };
                    let s = start(
                        StartInput {
                            cell: &cell,
                            atoms: &atoms,
                            potentials: &potentials,
                            options: &options,
                        },
                        comm,
                    )
                    .unwrap();
                    let bits: Vec<_> = s
                        .orbitals
                        .as_slice()
                        .iter()
                        .chain(&s.ritz_values_hartree)
                        .chain(&s.density)
                        .chain(&s.ionic)
                        .chain(&s.potential)
                        .chain(&s.diagonal_static)
                        .map(|v| v.to_bits())
                        .chain(
                            s.density_g
                                .iter()
                                .flat_map(|z| [z.re.to_bits(), z.im.to_bits()]),
                        )
                        .collect();
                    assert!(!s.setup.phases.contains_key("atomic_qr_mgs_fallback"));
                    let p = &s.setup.phases;
                    assert!(
                        p["atomic_qr_gram"]
                            + p["atomic_qr_factor_inverse"]
                            + p["atomic_qr_rotate"]
                            + p["atomic_completion"]
                            <= p["atomic_orthonormalize_complete"]
                    );
                    assert!(p["atomic_h_apply"] + p["atomic_h_gram"] <= p["atomic_h_apply_gram"]);
                    if comm.rank() == 0 && threads == 1 {
                        let m = &s.memory;
                        let (l, j, b) = (m.owned_rows, m.atomic_columns, m.bands);
                        let k = j.max(b);
                        eprintln!(
                            "P3-SETUP2 {label} ranks={ranks} D={diagonal} O={orbital}: G={} L={l} J={j} B={b} F={} qr_flops={} hgram_flops={} rotate_flops={} qr_reduce_bytes={} hgram_reduce_bytes={} H_columns={k} {:?}",
                            m.wave_rows,
                            m.grid_points,
                            8 * l * j * j,
                            2 * l * k * k,
                            2 * l * k * b,
                            16 * j * j,
                            8 * k * k,
                            s.setup
                        );
                    }
                    bits
                })
            };
            let baseline = build("off", "off", 1);
            for (diagonal, orbital) in [("on", "off"), ("off", "on"), ("on", "on")] {
                for threads in [1, 2] {
                    assert_eq!(
                        baseline,
                        build(diagonal, orbital, threads),
                        "{label} ranks={ranks} D={diagonal} O={orbital} threads={threads}"
                    );
                }
            }
        }
    }
}

/// P3-SETUP2: the complete verified trajectory, not only the starting energy.
#[test]
fn setup2_verified_trajectory_identity() {
    use dft_numerics::switches::RunConfig;
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    for (label, rows) in [
        ("cubic", [[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]]),
        ("skew", [[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]]),
        ("elongated", [[5., 0., 0.], [0., 5., 0.], [0., 0., 10.]]),
    ] {
        let cell = dft_geometry::Cell::from_user_rows_bohr(rows).unwrap();
        let atoms = [
            dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32]),
            dft_geometry::Atom::new("Xx", [0.57, 0.68, 0.71]),
        ];
        let potentials = [dft_engine::SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let options = dft_engine::ScfOptions {
            cutoff_hartree: 2.,
            bands: Some(5),
            max_iterations: 100,
            gamma_real: dft_engine::GammaReal::On,
            ..Default::default()
        };
        let mut baseline = None;
        for (diagonal, orbital) in [("off", "off"), ("on", "off"), ("off", "on"), ("on", "on")] {
            let config = RunConfig::from_pairs(&[
                ("refstart_diagonal_reuse", diagonal),
                ("refstart_orbital_reuse", orbital),
                ("refscf_verify_energy", "off"),
            ])
            .unwrap();
            let _scope = config.enter();
            let run = dft_engine::refscf::driver::run_with(
                dft_engine::refstart::StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &options,
                },
                Profile::default(),
                &Serial,
            )
            .unwrap();
            assert!(
                run.outcome.converged && run.outcome.result.verification.as_ref().unwrap().passed()
            );
            let science = stop_to_verified_science(&run);
            if let Some(previous) = &baseline {
                assert_eq!(&science, previous);
            } else {
                baseline = Some(science);
            }
            let calls: Vec<_> = run
                .records
                .iter()
                .flat_map(|r| &r.calls)
                .chain(&run.verification_calls)
                .collect();
            eprintln!(
                "P3-SETUP2 verified {label} D={diagonal} O={orbital}: E={:.17e} Ha maps={} Dav={} H={} setup={:?}",
                run.outcome.result.energy.total_hartree(),
                run.records.len(),
                calls.iter().map(|c| c.iterations).sum::<usize>(),
                calls.iter().map(|c| c.h_columns).sum::<usize>(),
                run.setup
            );
        }
    }
}

/// Cross the 256-row GEMM panel boundary with more atoms/bands than the small
/// trajectory fixtures. Timings are observations; this test has no speed gate.
#[test]
fn setup2_larger_start_profile() {
    use dft_engine::refstart::{StartInput, start};
    use dft_numerics::switches::RunConfig;
    let pp = common::potential_with_atomic_start(4., 1.8, true, 161, "SLA PZ");
    let atoms: Vec<_> = (0..4)
        .flat_map(|i| {
            (0..2).flat_map(move |j| {
                (0..2).map(move |k| {
                    dft_geometry::Atom::new("Xx", [i as f64 / 4., j as f64 / 2., k as f64 / 2.])
                })
            })
        })
        .collect();
    let potentials = [dft_engine::SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = dft_engine::ScfOptions {
        cutoff_hartree: 4.,
        bands: Some(34),
        gamma_real: dft_engine::GammaReal::On,
        ..Default::default()
    };
    for (label, rows) in [
        ("large-cubic", [[14., 0., 0.], [0., 14., 0.], [0., 0., 14.]]),
        (
            "large-skew",
            [[14., 0., 0.], [0.7, 14., 0.], [0.2, 0.3, 14.]],
        ),
        (
            "large-elongated",
            [[10., 0., 0.], [0., 10., 0.], [0., 0., 24.]],
        ),
    ] {
        let cell = dft_geometry::Cell::from_user_rows_bohr(rows).unwrap();
        let mut baseline = None;
        for (diagonal, orbital) in [("off", "off"), ("on", "off"), ("off", "on"), ("on", "on")] {
            let config = RunConfig::from_pairs(&[
                ("refstart_diagonal_reuse", diagonal),
                ("refstart_orbital_reuse", orbital),
            ])
            .unwrap();
            let _scope = config.enter();
            let s = start(
                StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &options,
                },
                &Serial,
            )
            .unwrap();
            let bits: Vec<_> = s
                .orbitals
                .as_slice()
                .iter()
                .chain(&s.ritz_values_hartree)
                .chain(&s.diagonal_static)
                .chain(&s.density)
                .chain(&s.ionic)
                .chain(&s.potential)
                .map(|v| v.to_bits())
                .chain(
                    s.density_g
                        .iter()
                        .flat_map(|z| [z.re.to_bits(), z.im.to_bits()]),
                )
                .collect();
            if let Some(previous) = &baseline {
                assert_eq!(&bits, previous);
            } else {
                baseline = Some(bits);
            }
            let k = s.memory.atomic_columns.max(s.memory.bands);
            let b = s.memory.bands;
            let h_collectives = (0..k)
                .step_by(b)
                .map(|first| s.operator.collectives_per_apply(b.min(k - first)))
                .fold((0, 0), |a, b| (a.0 + b.0, a.1 + b.1));
            eprintln!(
                "P3-SETUP2 {label} D={diagonal} O={orbital}: memory={:?} H_collectives={h_collectives:?} setup={:?}",
                s.memory, s.setup
            );
        }
    }
}

/// P3-SETUP3: compare invariant Ritz projectors (signs and degenerate
/// eigenvectors are arbitrary), full spectra and all non-orbital start bits.
#[test]
fn setup3_generalized_start_and_compact_g() {
    use dft_engine::refstart::{StartInput, start};
    use dft_numerics::switches::{Id, RunConfig};
    for id in [Id::REFSTART_GENERALIZED, Id::REFSTART_COMPACT_G] {
        assert!(RunConfig::from_pairs(&[]).unwrap().on(id)); // default on since 2026-11-25
    }
    for (shape, rows) in [
        [[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]],
        [[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]],
        [[5., 0., 0.], [0., 5., 0.], [0., 0., 10.]],
        [[8., 0., 0.], [0., 8., 0.], [0., 0., 4.]],
    ]
    .into_iter()
    .enumerate()
    {
        let cell = dft_geometry::Cell::from_user_rows_bohr(rows).unwrap();
        let pp = common::potential_with_atomic_start(2., 1.8, shape != 3, 161, "SLA PZ");
        let pp2 = common::potential_with_atomic_start(2., 2., false, 161, "SLA PZ");
        let mut atoms = vec![
            dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32]),
            dft_geometry::Atom::new("Xx", [0.57, 0.68, 0.71]),
        ];
        let mut potentials = vec![dft_engine::SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        if shape == 1 {
            atoms.push(dft_geometry::Atom::new("Yy", [0.1, 0.7, 0.8]));
            potentials.push(dft_engine::SpeciesPotential {
                label: "Yy",
                potential: &pp2,
            });
        }
        for bands in [5, 4 * atoms.len(), 4 * atoms.len() + 1] {
            for ranks in [1, 3, 8] {
                dft_mpi::simulate(ranks, |comm| {
                    let build = |generalized, compact, threads| {
                        let config = RunConfig::from_pairs(&[
                            ("refstart_generalized", generalized),
                            ("refstart_compact_g", compact),
                        ])
                        .unwrap();
                        let _scope = config.enter();
                        let options = dft_engine::ScfOptions {
                            cutoff_hartree: 2.,
                            bands: Some(bands),
                            threads,
                            gamma_real: dft_engine::GammaReal::On,
                            ..Default::default()
                        };
                        start(
                            StartInput {
                                cell: &cell,
                                atoms: &atoms,
                                potentials: &potentials,
                                options: &options,
                            },
                            comm,
                        )
                        .unwrap()
                    };
                    let base = build("off", "off", 1);
                    let bits = |s: &dft_engine::refstart::StartState| -> Vec<u64> {
                        s.density
                            .iter()
                            .chain(&s.ionic)
                            .chain(&s.potential)
                            .chain(&s.diagonal_static)
                            .map(|x| x.to_bits())
                            .chain(
                                s.density_g
                                    .iter()
                                    .flat_map(|z| [z.re.to_bits(), z.im.to_bits()]),
                            )
                            .collect()
                    };
                    for (generalized, compact) in [("off", "on"), ("on", "off"), ("on", "on")] {
                        let mut worker_bits = None;
                        for threads in [1, 2] {
                            let candidate = build(generalized, compact, threads);
                            assert_eq!(bits(&base), bits(&candidate));
                            assert_eq!(base.grid.g_to_fft(), candidate.grid.g_to_fft());
                            assert_eq!(base.basis.vectors(), candidate.basis.vectors());
                            let active = generalized == "on" && bands <= base.memory.atomic_columns;
                            assert_eq!(
                                candidate
                                    .setup
                                    .phases
                                    .contains_key("atomic_generalized_diagonalize_rotate"),
                                active
                            );
                            assert_eq!(
                                candidate.setup.phases.contains_key("atomic_qr_gram"),
                                !active
                            );
                            let ev_error = base
                                .ritz_values_hartree
                                .iter()
                                .zip(&candidate.ritz_values_hartree)
                                .map(|(a, b)| (a - b).abs())
                                .fold(0., f64::max);
                            assert!(
                                ev_error < 1e-12,
                                "shape={shape} B={bands} ranks={ranks} eigenvalues {ev_error:e} Ha"
                            );
                            let mut local = vec![0.; 2 * bands * bands];
                            for j in 0..bands {
                                for i in 0..bands {
                                    for row in 0..base.orbitals.rows() {
                                        local[i + j * bands] +=
                                            base.orbitals[(row, i)] * candidate.orbitals[(row, j)];
                                        local[bands * bands + i + j * bands] += candidate.orbitals
                                            [(row, i)]
                                            * candidate.orbitals[(row, j)];
                                    }
                                }
                            }
                            let mut overlaps = vec![0.; local.len()];
                            comm.tree_sum_f64(&local, &mut overlaps).unwrap();
                            let mut orth_error: f64 = 0.;
                            for j in 0..bands {
                                for i in 0..bands {
                                    orth_error = orth_error.max(
                                        (overlaps[bands * bands + i + j * bands]
                                            - f64::from(i == j))
                                        .abs(),
                                    );
                                }
                            }
                            let span_error = setup3_projection_residual(
                                comm,
                                &base.orbitals,
                                &candidate.orbitals,
                                &overlaps[..bands * bands],
                            );
                            assert!(orth_error < 1e-12, "orthogonality {orth_error:e}");
                            assert!(
                                span_error < 1e-12,
                                "shape={shape} B={bands} ranks={ranks} subspace {span_error:e}"
                            );
                            let orbital_bits: Vec<_> = candidate
                                .orbitals
                                .as_slice()
                                .iter()
                                .chain(&candidate.ritz_values_hartree)
                                .map(|v| v.to_bits())
                                .collect();
                            if !active {
                                assert_eq!(base.orbitals.as_slice(), candidate.orbitals.as_slice());
                                assert_eq!(base.ritz_values_hartree, candidate.ritz_values_hartree);
                            }
                            if let Some(previous) = &worker_bits {
                                assert_eq!(&orbital_bits, previous);
                            }
                            worker_bits = Some(orbital_bits);
                            if active && threads == 1 && compact == "on" && comm.rank() == 0 {
                                eprintln!(
                                    "P3-SETUP3 start shape={shape} P={ranks} G={} L={} J={} B={bands} eigen_error={ev_error:.9e} Ha orth_error={orth_error:.9e} span_error={span_error:.9e}",
                                    base.memory.wave_rows,
                                    base.memory.owned_rows,
                                    base.memory.atomic_columns
                                );
                            }
                        }
                    }
                });
            }
        }
    }
}

/// Default convergence thresholds and independent verification remain enabled.
/// Count comparisons cover every solve, redo, restart, repair and verification.
#[test]
fn setup3_generalized_verified_trajectories() {
    use dft_numerics::switches::RunConfig;
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    for (label, rows) in [
        ("cubic", [[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]]),
        ("skew", [[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]]),
        ("elongated", [[5., 0., 0.], [0., 5., 0.], [0., 0., 10.]]),
    ] {
        let cell = dft_geometry::Cell::from_user_rows_bohr(rows).unwrap();
        let atoms = [
            dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32]),
            dft_geometry::Atom::new("Xx", [0.57, 0.68, 0.71]),
        ];
        let potentials = [dft_engine::SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let options = dft_engine::ScfOptions {
            cutoff_hartree: 2.,
            bands: Some(5),
            max_iterations: 100,
            gamma_real: dft_engine::GammaReal::On,
            ..Default::default()
        };
        for ranks in [1, 3, 8] {
            dft_mpi::simulate(ranks, |comm| {
                let run = |generalized, compact| {
                    let config = RunConfig::from_pairs(&[
                        ("refstart_generalized", generalized),
                        ("refstart_compact_g", compact),
                        ("refscf_verify_energy", "off"),
                    ])
                    .unwrap();
                    let _scope = config.enter();
                    dft_engine::refscf::driver::run_with(
                        dft_engine::refstart::StartInput {
                            cell: &cell,
                            atoms: &atoms,
                            potentials: &potentials,
                            options: &options,
                        },
                        Profile::default(),
                        comm,
                    )
                    .unwrap()
                };
                let base = run("off", "off");
                for (generalized, compact) in [("off", "on"), ("on", "off"), ("on", "on")] {
                    let candidate = run(generalized, compact);
                    assert!(
                        candidate.outcome.converged
                            && candidate
                                .outcome
                                .result
                                .verification
                                .as_ref()
                                .unwrap()
                                .passed()
                    );
                    let counts = |r: &dft_engine::refscf::driver::Run| {
                        let maps: Vec<_> = r
                            .records
                            .iter()
                            .map(|s| (s.first_iteration_redo, s.verification_calls.clone()))
                            .collect();
                        let solves: Vec<_> = r
                            .records
                            .iter()
                            .flat_map(|s| &s.calls)
                            .chain(&r.verification_calls)
                            .map(|s| {
                                (
                                    s.iterations,
                                    s.h_columns,
                                    s.restarts,
                                    s.repairs,
                                    s.dropped_corrections,
                                    s.reduced_orders.clone(),
                                    s.active_bands.clone(),
                                )
                            })
                            .collect();
                        (maps, solves)
                    };
                    let delta = candidate.outcome.result.energy.total_hartree()
                        - base.outcome.result.energy.total_hartree();
                    let dr2 = candidate.records[0].dr2_ry - base.records[0].dr2_ry;
                    let calls: Vec<_> = candidate
                        .records
                        .iter()
                        .flat_map(|s| &s.calls)
                        .chain(&candidate.verification_calls)
                        .collect();
                    if comm.rank() == 0 {
                        eprintln!(
                            "P3-SETUP3 verified {label} P={ranks} G={generalized} C={compact}: E={:.17e} Ha delta={delta:.9e} Ha first_dr2={:.17e} Ry delta_dr2={dr2:.9e} Ry maps={} Dav={} H={}",
                            candidate.outcome.result.energy.total_hartree(),
                            candidate.records[0].dr2_ry,
                            candidate.records.len(),
                            calls.iter().map(|s| s.iterations).sum::<usize>(),
                            calls.iter().map(|s| s.h_columns).sum::<usize>()
                        );
                    }
                    assert!(delta.abs() < 1e-12);
                    // First-map dr2 is a diagnostic, not a converged energy:
                    // the ticket requires its signed delta, not a 1e-12 Ry
                    // equality gate on the loose first Davidson solve.
                    assert!(dr2.is_finite());
                    assert_eq!(
                        counts(&base),
                        counts(&candidate),
                        "{label} G={generalized} C={compact}"
                    );
                    if generalized == "off" {
                        assert_eq!(
                            stop_to_verified_science(&base),
                            stop_to_verified_science(&candidate)
                        );
                    }
                }
            });
        }
    }
}

#[test]
fn setup3_guarded_fallback_rank_failure_and_empty_owners() {
    use dft_engine::refstart::{StartInput, start};
    use dft_numerics::switches::RunConfig;
    let cell = dft_geometry::Cell::from_user_rows_bohr([[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]])
        .unwrap();
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    let potentials = [dft_engine::SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    // A seven-row wave sphere on eight ranks exercises empty owners. Small
    // replicated algebra must still have J-row capacity on every rank.
    let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
    let empty = dft_mpi::simulate(8, |comm| {
        let config =
            RunConfig::from_pairs(&[("refstart_generalized", "on"), ("refstart_compact_g", "on")])
                .unwrap();
        let _scope = config.enter();
        let options = dft_engine::ScfOptions {
            cutoff_hartree: 0.5,
            bands: Some(2),
            ..Default::default()
        };
        let s = start(
            StartInput {
                cell: &cell,
                atoms: &atoms,
                potentials: &potentials,
                options: &options,
            },
            comm,
        )
        .unwrap();
        assert!(
            s.setup
                .phases
                .contains_key("atomic_generalized_diagonalize_rotate")
        );
        s.orbitals.rows() == 0
    });
    assert!(empty.iter().any(|&e| e));
    // Near-coincident atomic centers stress the overlap metric; coincident
    // centers are rank deficient and must not be silently truncated or repaired.
    for displacement in [1e-5, 0.] {
        let atoms = [
            dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32]),
            dft_geometry::Atom::new("Xx", [0.13 + displacement, 0.21, 0.32]),
        ];
        dft_mpi::simulate(3, |comm| {
            let build = |value| {
                let config = RunConfig::from_pairs(&[("refstart_generalized", value)]).unwrap();
                let _scope = config.enter();
                let options = dft_engine::ScfOptions {
                    cutoff_hartree: 2.,
                    bands: Some(5),
                    ..Default::default()
                };
                start(
                    StartInput {
                        cell: &cell,
                        atoms: &atoms,
                        potentials: &potentials,
                        options: &options,
                    },
                    comm,
                )
            };
            let base = build("off");
            let candidate = build("on");
            if displacement == 0. {
                assert!(base.is_err() && candidate.is_err());
            } else {
                let base = base.unwrap();
                let candidate = candidate.unwrap();
                assert!(candidate.setup.phases.contains_key("atomic_overlap_gram"));
                assert!(candidate.setup.phases.contains_key("atomic_qr_gram"));
                assert!(
                    !candidate
                        .setup
                        .phases
                        .contains_key("atomic_generalized_diagonalize_rotate")
                );
                assert_eq!(base.orbitals.as_slice(), candidate.orbitals.as_slice());
                assert_eq!(base.ritz_values_hartree, candidate.ritz_values_hartree);
            }
        });
    }
}

/// More than one GEMM row panel and 64 raw atomic columns; no timing gate.
#[test]
fn setup3_generalized_large_start() {
    use dft_engine::refstart::{StartInput, start};
    use dft_numerics::switches::RunConfig;
    let pp = common::potential_with_atomic_start(4., 1.8, true, 161, "SLA PZ");
    let atoms: Vec<_> = (0..16)
        .map(|n| {
            dft_geometry::Atom::new(
                "Xx",
                [
                    (n / 4) as f64 / 4. + 0.001 * (n % 3) as f64,
                    (n / 2 % 2) as f64 / 2. + 0.002 * (n % 5) as f64,
                    (n % 2) as f64 / 2. + 0.003 * (n % 7) as f64,
                ],
            )
        })
        .collect();
    let potentials = [dft_engine::SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    for (label, rows) in [
        ("cubic", [[14., 0., 0.], [0., 14., 0.], [0., 0., 14.]]),
        ("skew", [[14., 0., 0.], [0.7, 14., 0.], [0.2, 0.3, 14.]]),
        ("elongated", [[10., 0., 0.], [0., 10., 0.], [0., 0., 24.]]),
    ] {
        let cell = dft_geometry::Cell::from_user_rows_bohr(rows).unwrap();
        let build = |generalized, threads| {
            let config = RunConfig::from_pairs(&[
                ("refstart_generalized", generalized),
                ("refstart_compact_g", "on"),
            ])
            .unwrap();
            let _scope = config.enter();
            let options = dft_engine::ScfOptions {
                cutoff_hartree: 4.,
                bands: Some(34),
                threads,
                ..Default::default()
            };
            start(
                StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &options,
                },
                &Serial,
            )
            .unwrap()
        };
        let base = build("off", 1);
        let candidate = build("on", 1);
        let threaded = build("on", 2);
        assert_eq!(candidate.orbitals.as_slice(), threaded.orbitals.as_slice());
        assert_eq!(candidate.ritz_values_hartree, threaded.ritz_values_hartree);
        assert!(
            candidate
                .setup
                .phases
                .contains_key("atomic_generalized_diagonalize_rotate")
        );
        let eigen_error = base
            .ritz_values_hartree
            .iter()
            .zip(&candidate.ritz_values_hartree)
            .map(|(a, b)| (a - b).abs())
            .fold(0., f64::max);
        let q = &candidate.orbitals;
        let mut orth_error: f64 = 0.;
        let mut projection = vec![0.; q.cols() * q.cols()];
        for j in 0..q.cols() {
            for i in 0..q.cols() {
                let norm: f64 = (0..q.rows()).map(|r| q[(r, i)] * q[(r, j)]).sum();
                orth_error = orth_error.max((norm - f64::from(i == j)).abs());
                projection[i + j * q.cols()] = (0..q.rows())
                    .map(|r| base.orbitals[(r, i)] * q[(r, j)])
                    .sum();
            }
        }
        let span_error = setup3_projection_residual(&Serial, &base.orbitals, q, &projection);
        eprintln!(
            "P3-SETUP3 large {label} G={} J=64 B=34 eigen_error={eigen_error:.9e} Ha orth_error={orth_error:.9e} span_error={span_error:.9e}",
            q.rows()
        );
        assert!(eigen_error < 1e-12 && orth_error < 1e-12 && span_error < 1e-12);
    }
}

// Linear (not squared-angle) subspace error: largest column norm of B-A(A^t B).
// Using only I-(A^t B)^t(A^t B) would hide first-order errors below sqrt(eps).
fn setup3_projection_residual<K: Collective>(
    comm: &K,
    a: &dft_numerics::RealMatrix,
    b: &dft_numerics::RealMatrix,
    overlap: &[f64],
) -> f64 {
    let mut local = vec![0.; b.cols()];
    for j in 0..b.cols() {
        for r in 0..b.rows() {
            let projected: f64 = (0..a.cols())
                .map(|k| a[(r, k)] * overlap[k + j * a.cols()])
                .sum();
            local[j] += (b[(r, j)] - projected).powi(2);
        }
    }
    let mut global = vec![0.; b.cols()];
    comm.tree_sum_f64(&local, &mut global).unwrap();
    global.into_iter().map(f64::sqrt).fold(0., f64::max)
}

// P3-RMM: deterministic synthetic PP oracles, never QE reference evidence.
// Both arms use the same 1e-12 Ry estimator; the delivered-energy gate is
// 1e-10 Ha. Default-estimator deviations are a documented acceptance limitation.
fn rmm_verified_fixture(label: &str, rows: [[f64; 3]; 3], nonlocal: bool, threshold_ry: f64) {
    use dft_numerics::switches::{Id, RunConfig};
    use dft_solvers::{ref_davidson::work, ref_rmm};
    assert!(!RunConfig::from_pairs(&[]).unwrap().on(Id::REFSCF_RMM));
    let pp = common::potential_with_atomic_start(2., 1.8, nonlocal, 161, "SLA PZ");
    let cell = dft_geometry::Cell::from_user_rows_bohr(rows).unwrap();
    let atoms = [
        dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32]),
        dft_geometry::Atom::new("Xx", [0.57, 0.68, 0.71]),
    ];
    let potentials = [dft_engine::SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = dft_engine::ScfOptions {
        cutoff_hartree: 2.,
        bands: Some(5),
        max_iterations: 120,
        gamma_real: dft_engine::GammaReal::On,
        ..Default::default()
    };
    let mut rank_energies = Vec::new();
    for ranks in [1, 3, 8] {
        let outputs = dft_mpi::simulate(ranks, |comm| {
            let run = |rmm, level, threads| {
                let options = dft_engine::ScfOptions {
                    threads,
                    ..options.clone()
                };
                let config = RunConfig::from_pairs(&[
                    ("refscf_rmm", rmm),
                    ("davidson_stages", level),
                    ("refscf_verify_energy", "on"),
                ])
                .unwrap();
                let _scope = config.enter();
                ref_rmm::take_records();
                work::take_solves();
                let run = dft_engine::refscf::driver::run_with(
                    dft_engine::refstart::StartInput {
                        cell: &cell,
                        atoms: &atoms,
                        potentials: &potentials,
                        options: &options,
                    },
                    Profile {
                        threshold_ry,
                        ..Default::default()
                    },
                    comm,
                )
                .unwrap();
                assert!(run.outcome.converged, "{label} ranks={ranks} arm={rmm}");
                assert!(run.outcome.result.verification.as_ref().unwrap().passed());
                let records = ref_rmm::take_records();
                assert!(records.iter().all(|r| ref_rmm::scheduled(r.scf_iteration)));
                if rmm == "on" && level == "on" {
                    assert!(!records.is_empty());
                    let expected: Vec<_> = run
                        .records
                        .iter()
                        .filter(|r| ref_rmm::scheduled(r.iteration) && !r.calls.is_empty())
                        .map(|r| r.iteration)
                        .collect();
                    assert_eq!(
                        records.iter().map(|r| r.scf_iteration).collect::<Vec<_>>(),
                        expected
                    );
                    for r in &records {
                        assert!(r.min_gs_norm > 1e-6);
                        let map = &run.records[r.scf_iteration - 1];
                        assert_eq!(map.calls.len(), 1, "one sweep; no inner retry");
                        assert_eq!(
                            map.calls[0].h_columns,
                            r.rotation_columns + r.step_columns.iter().sum::<usize>()
                        );
                    }
                } else {
                    assert!(records.is_empty());
                }
                let calls: Vec<_> = run
                    .records
                    .iter()
                    .flat_map(|r| &r.calls)
                    .chain(&run.verification_calls)
                    .collect();
                let hcols: usize = calls.iter().map(|c| c.h_columns).sum();
                let dense: u64 = work::take_solves()
                    .iter()
                    .flat_map(|w| &w.stages)
                    .map(|s| s.gemm_flops)
                    .sum::<u64>()
                    + records
                        .iter()
                        .map(|r| r.gemm_flops.iter().sum::<u64>())
                        .sum::<u64>();
                let energy = run.outcome.result.energy.total_hartree();
                if level == "on" && comm.rank() == 0 {
                    eprintln!(
                        "P3-RMM {label} threshold={threshold_ry:e} Ry P={ranks} arm={rmm} E={energy:.17e} Ha iterations={} maps={} verification={} H={hcols} GEMM_FLOPs_rank0={dense} RMM={} min_GS={:.9e}",
                        run.records.len(),
                        calls.len(),
                        run.verification_calls.len(),
                        records.len(),
                        records.iter().map(|r| r.min_gs_norm).fold(1., f64::min)
                    );
                }
                let science = stop_to_verified_science(&run);
                (energy, hcols, dense, science)
            };
            let d = run("off", "on", 1);
            let m = run("on", "on", 1);
            let repeat = run("on", "off", 2);
            assert_eq!(m.0.to_bits(), repeat.0.to_bits());
            assert_eq!(m.1, repeat.1);
            assert_eq!(m.3, repeat.3, "diagnostics must not change science");
            let delta = m.0 - d.0;
            if comm.rank() == 0 {
                eprintln!(
                    "P3-RMM {label} threshold={threshold_ry:e} Ry P={ranks} delta_M_D={delta:.9e} Ha"
                );
            }
            assert!(delta.abs() < 1e-10, "{label} P={ranks}: {delta:e} Ha");
            (d.0, m.0)
        });
        for out in &outputs[1..] {
            assert_eq!(out, &outputs[0]);
        }
        rank_energies.push(outputs[0]);
    }
    for &(d, m) in &rank_energies[1..] {
        assert!((d - rank_energies[0].0).abs() < 1e-10);
        assert!((m - rank_energies[0].1).abs() < 1e-10);
    }
}
#[test]
fn rmm_verified_cubic() {
    rmm_verified_fixture(
        "cubic",
        [[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]],
        false,
        1e-12,
    );
}
#[test]
fn rmm_verified_skew() {
    rmm_verified_fixture(
        "skew",
        [[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]],
        false,
        1e-12,
    );
}
#[test]
fn rmm_verified_elongated() {
    rmm_verified_fixture(
        "elongated",
        [[5., 0., 0.], [0., 5., 0.], [0., 0., 10.]],
        false,
        1e-12,
    );
}
#[test]
fn rmm_verified_nonlocal() {
    rmm_verified_fixture(
        "nonlocal",
        [[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]],
        true,
        1e-12,
    );
}

// Retain the failed 1e-10 Ha acceptance experiment as an explicit reproducer.
// These are known failures, not runtime exclusions or passing acceptance gates.
// cargo test -p dft-engine --test refscf rmm_default_ -- --ignored --nocapture
#[test]
#[ignore = "known default-estimator strict energy-gate failure; see rmm-dev08-request.md"]
fn rmm_default_cubic_strict_energy_gate() {
    rmm_verified_fixture(
        "cubic-default",
        [[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]],
        false,
        1e-8,
    );
}
#[test]
#[ignore = "known default-estimator strict energy-gate failure; see rmm-dev08-request.md"]
fn rmm_default_skew_strict_energy_gate() {
    rmm_verified_fixture(
        "skew-default",
        [[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]],
        false,
        1e-8,
    );
}
#[test]
#[ignore = "known default-estimator strict energy-gate failure; see rmm-dev08-request.md"]
fn rmm_default_elongated_strict_energy_gate() {
    rmm_verified_fixture(
        "elongated-default",
        [[5., 0., 0.], [0., 5., 0.], [0., 0., 10.]],
        false,
        1e-8,
    );
}
#[test]
#[ignore = "known default-estimator strict energy-gate failure; see rmm-dev08-request.md"]
fn rmm_default_nonlocal_strict_energy_gate() {
    rmm_verified_fixture(
        "nonlocal-default",
        [[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]],
        true,
        1e-8,
    );
}
