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
    let config = dft_numerics::switches::RunConfig::from_pairs(&[]).unwrap();
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
        let scope = RunConfig::from_pairs(&[("davidson_stages", level)]).unwrap();
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
        ("off", "off", vec![]),
        ("off", "on", vec![]),
        ("default", "off", vec![]),
        ("default", "on", vec![]),
        ("first", "off", vec![]),
        ("reuse", "off", vec![]),
        ("gate", "off", vec![]),
        ("both", "off", vec![]),
    ] {
        let mut pairs = pairs;
        pairs.push(("davidson_stages", level));
        let scope = RunConfig::from_pairs(&pairs).unwrap();
        let reuse = if false { "on" } else { "off" };
        let gate = if true { "on" } else { "off" };
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
            let build = |_reuse, threads| {
                dft_mpi::simulate(ranks, |comm| {
                    let config = RunConfig::from_pairs(&[
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
            let config = RunConfig::from_pairs(&[("davidson_stages", diagnostics)]).unwrap();
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
                    let config = RunConfig::from_pairs(&[]).unwrap();
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
                    if p.contains_key("atomic_qr_gram") {
                        assert!(
                            p["atomic_qr_gram"]
                                + p["atomic_qr_factor_inverse"]
                                + p["atomic_qr_rotate"]
                                + p["atomic_completion"]
                                <= p["atomic_orthonormalize_complete"]
                        );
                    }
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
            let config = RunConfig::from_pairs(&[]).unwrap();
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
            let config = RunConfig::from_pairs(&[]).unwrap();
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
    use dft_numerics::switches::RunConfig;

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
                    let build = |_generalized, _compact, threads| {
                        let config = RunConfig::from_pairs(&[]).unwrap();
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
                            let active = bands <= base.memory.atomic_columns;
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
                let run = |_generalized, _compact| {
                    let config = RunConfig::from_pairs(&[]).unwrap();
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
        let config = RunConfig::from_pairs(&[]).unwrap();
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
            let build = |_value| {
                let config = RunConfig::from_pairs(&[]).unwrap();
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
        let build = |_generalized, threads| {
            let config = RunConfig::from_pairs(&[]).unwrap();
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

#[test]
fn lean_energy_preserves_delivery_and_counts() {
    use dft_engine::refscf::driver::run_with;
    use dft_engine::{ScfOptions, SpeciesPotential, Verification, refstart::StartInput};
    use dft_numerics::switches::RunConfig;
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    for rows in [
        [[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]],
        [[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]],
        [[5., 0., 0.], [0., 5., 0.], [0., 0., 12.]],
    ] {
        let cell = dft_geometry::Cell::from_user_rows_bohr(rows).unwrap();
        for (verification, max_iterations, _energy_rule) in [
            (Verification::Off, 100, "off"),
            (Verification::FixedPoint, 100, "off"),
            (Verification::FixedPoint, 100, "on"),
            (Verification::Off, 1, "off"),
        ] {
            dft_mpi::simulate(2, |comm| {
                let options = ScfOptions {
                    cutoff_hartree: 1.,
                    bands: Some(3),
                    gamma_real: dft_engine::GammaReal::On,
                    max_iterations,
                    verification,
                    ..Default::default()
                };
                let run = |_lean| {
                    let config = RunConfig::from_pairs(&[]).unwrap();
                    let _scope = config.enter();
                    run_with(
                        StartInput {
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
                let a = run("off");
                let b = run("on");
                assert_eq!(a.first_reference_stop, b.first_reference_stop);
                assert_eq!(a.verification_calls.len(), b.verification_calls.len());
                for (a, b) in a.records.iter().zip(&b.records) {
                    assert_eq!(
                        a.variational_energy_ha.to_bits(),
                        b.variational_energy_ha.to_bits()
                    );
                    assert_eq!(a.dr2_ry.to_bits(), b.dr2_ry.to_bits());
                    if !b.energy.kinetic.is_finite() {
                        assert!(b.energy.band_sum.is_finite());
                        assert!(b.energy.double_counting.is_finite());
                        assert!(b.energy.functional_difference_bound.is_nan());
                    }
                }
                let ar = &a.outcome.result;
                let br = &b.outcome.result;
                assert_eq!(a.outcome.converged, b.outcome.converged);
                assert_eq!(
                    ar.energy.total_hartree().to_bits(),
                    br.energy.total_hartree().to_bits()
                );
                assert_eq!(ar.density, br.density);
                assert_eq!(ar.history.len(), br.history.len());
                for (a, b) in ar.history.iter().zip(&br.history) {
                    assert_eq!(a.davidson_iterations, b.davidson_iterations);
                    assert_eq!(a.operator_applications, b.operator_applications);
                    assert_eq!(a.standard_error_ha, b.standard_error_ha);
                }
                assert!(
                    br.history
                        .last()
                        .unwrap()
                        .energy
                        .total_hartree()
                        .is_finite()
                );
            });
        }
    }
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
                let config =
                    RunConfig::from_pairs(&[("refscf_rmm", rmm), ("davidson_stages", level)])
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

/// Plans are cached per run, not globally or by grid shape alone. Exercise
/// different sparse supports, axis choices, empty ranks and private workspaces.
#[test]
fn setup4_fft_plan_cache_and_workspace_identity() {
    use dft_numerics::switches::RunConfig;
    use dft_solvers::distributed_fft::{DistributedPairFft, PlanCache};
    for shape in [[4, 4, 4], [3, 5, 7]] {
        let indices: Vec<_> = (0..shape[0])
            .flat_map(|i| (0..shape[1]).flat_map(move |j| (0..shape[2]).map(move |k| [i, j, k])))
            .collect();
        for ranks in [1, 8] {
            dft_mpi::simulate(ranks, |comm| {
                for vendor in [false, true] {
                    let config = RunConfig::from_pairs(&[
                        ("fft_fused_rows", "on"),
                        ("fft_sparse_clear", "on"),
                    ])
                    .unwrap();
                    let _config = config.enter();
                    let cache = PlanCache::default();
                    for sparse in [false, true] {
                        let used = if sparse {
                            &indices[..shape[2]]
                        } else {
                            &indices[..]
                        };
                        let mut private =
                            DistributedPairFft::with_lanes(shape, used, comm.rank(), ranks, vendor)
                                .unwrap();
                        assert!(private.shared_workspace().is_none());
                        let scope = cache.enter();
                        let mut shared =
                            DistributedPairFft::with_lanes(shape, used, comm.rank(), ranks, vendor)
                                .unwrap();
                        let counts = cache.counts();
                        let repeated =
                            DistributedPairFft::with_lanes(shape, used, comm.rank(), ranks, vendor)
                                .unwrap();
                        assert_eq!(cache.counts(), counts, "no repeated kernel plans");
                        let mut fork = shared.shared_workspace().unwrap();
                        assert_eq!(cache.counts(), counts, "no planning in workspace copy");
                        drop(repeated);
                        drop(scope);
                        let local: Vec<_> = private
                            .local_positions()
                            .iter()
                            .map(|&i| C::new((i as f64 * 0.37).sin(), (i as f64 * 0.21).cos()))
                            .collect();
                        let potential: Vec<_> = (0..shape.iter().product())
                            .map(|i| (i as f64 * 0.03).cos())
                            .collect();
                        for scale in [1., 0.2, -0.7] {
                            let input: Vec<_> = local.iter().map(|z| z * scale).collect();
                            let mut expected = input.clone();
                            private
                                .apply_local_batch(comm, &mut expected, &potential, 1)
                                .unwrap();
                            for fft in [&mut shared, &mut fork] {
                                let mut actual = input.clone();
                                fft.apply_local_batch(comm, &mut actual, &potential, 1)
                                    .unwrap();
                                let bits = |v: &[C]| {
                                    v.iter()
                                        .flat_map(|z| [z.re.to_bits(), z.im.to_bits()])
                                        .collect::<Vec<_>>()
                                };
                                assert_eq!(bits(&actual), bits(&expected));
                            }
                        }
                    }
                    if !vendor || !dft_numerics::MKL_BACKEND {
                        assert_eq!(
                            cache.counts(),
                            (if shape[0] == shape[1] { 2 } else { 6 }, 0)
                        );
                    } else {
                        assert!(cache.counts().1 > 0);
                    }
                    assert_eq!(PlanCache::default().counts(), (0, 0));
                }
            });
        }
    }
}

/// Cubic, skew multi-species, elongated and slab starts, including J>B and
/// kinetic completion, odd pairs, nonlocal projectors and empty row owners.
#[test]
fn setup4_start_switches_preserve_fields_and_ritz() {
    use dft_engine::refstart::{StartInput, start};
    use dft_numerics::switches::RunConfig;
    // FFT_SHARE and SINGLE_GEMM default on since 2026-11-26; LOCAL_SHELLS stays off.

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
        let atoms = [
            dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32]),
            dft_geometry::Atom::new("Yy", [0.57, 0.68, 0.71]),
        ];
        let potentials = [
            dft_engine::SpeciesPotential {
                label: "Xx",
                potential: &pp,
            },
            dft_engine::SpeciesPotential {
                label: "Yy",
                potential: if shape == 1 { &pp2 } else { &pp },
            },
        ];
        for ranks in [1, 3, 8] {
            dft_mpi::simulate(ranks, |comm| {
                let build = |_fft, _gram, _shells, threads| {
                    let config = RunConfig::from_pairs(&[]).unwrap();
                    let _scope = config.enter();
                    let options = dft_engine::ScfOptions {
                        cutoff_hartree: if shape == 0 { 8. } else { 2. },
                        bands: Some(if shape == 3 { 9 } else { 5 }),
                        gamma_real: dft_engine::GammaReal::On,
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
                        comm,
                    )
                    .unwrap()
                };
                let base = build("off", "off", "off", 1);
                if shape == 0 && ranks == 1 {
                    assert!(
                        base.orbitals.rows() > 256,
                        "exercise a multi-panel Gram reduction"
                    );
                }
                let fields = |s: &dft_engine::refstart::StartState| {
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
                        .collect::<Vec<_>>()
                };
                for (fft, gram, shells) in [
                    ("on", "off", "off"),
                    ("off", "off", "on"),
                    ("off", "on", "off"),
                    ("on", "on", "on"),
                ] {
                    let mut previous = None;
                    for threads in [1, 2] {
                        let mut candidate = build(fft, gram, shells, threads);
                        assert_eq!(fields(&base), fields(&candidate));
                        let delta = base
                            .ritz_values_hartree
                            .iter()
                            .zip(&candidate.ritz_values_hartree)
                            .map(|(a, b)| (a - b).abs())
                            .fold(0., f64::max);
                        assert!(delta < 1e-12, "Ritz error {delta:e} Ha");
                        let bands = candidate.orbitals.cols();
                        let mut dots = vec![0.; 2 * bands * bands];
                        for j in 0..bands {
                            for i in 0..bands {
                                for r in 0..candidate.orbitals.rows() {
                                    dots[i + j * bands] +=
                                        base.orbitals[(r, i)] * candidate.orbitals[(r, j)];
                                    dots[bands * bands + i + j * bands] +=
                                        candidate.orbitals[(r, i)] * candidate.orbitals[(r, j)];
                                }
                            }
                        }
                        let mut global = vec![0.; dots.len()];
                        comm.tree_sum_f64(&dots, &mut global).unwrap();
                        for j in 0..bands {
                            for i in 0..bands {
                                assert!(
                                    (global[bands * bands + i + j * bands] - f64::from(i == j))
                                        .abs()
                                        < 1e-12
                                );
                            }
                        }
                        assert!(
                            setup3_projection_residual(
                                comm,
                                &base.orbitals,
                                &candidate.orbitals,
                                &global[..bands * bands]
                            ) < 1e-12
                        );

                        let bits = candidate
                            .orbitals
                            .as_slice()
                            .iter()
                            .chain(&candidate.ritz_values_hartree)
                            .map(|v| v.to_bits())
                            .collect::<Vec<_>>();
                        if gram == "off" || !dft_numerics::MKL_BACKEND {
                            assert_eq!(
                                base.orbitals
                                    .as_slice()
                                    .iter()
                                    .chain(&base.ritz_values_hartree)
                                    .map(|v| v.to_bits())
                                    .collect::<Vec<_>>(),
                                bits
                            );
                        }
                        if let Some(prev) = &previous {
                            assert_eq!(prev, &bits);
                        }
                        previous = Some(bits);
                        // The same tuned operator remains usable after the
                        // constructor cache scope has ended.
                        let mut hx = dft_numerics::RealMatrix::zeros(
                            candidate.orbitals.rows(),
                            candidate.orbitals.cols(),
                        )
                        .unwrap();
                        candidate
                            .operator
                            .apply_block(comm, &candidate.orbitals, &mut hx)
                            .unwrap();
                        assert!(hx.as_slice().iter().all(|v| v.is_finite()));
                    }
                }
            });
        }
    }
}

#[test]
fn setup4_shared_density_plan_verified_trajectory() {
    use dft_numerics::switches::RunConfig;
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    for rows in [
        [[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]],
        [[5., 0., 0.], [0., 5., 0.], [0., 0., 10.]],
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
        for ranks in [1, 3] {
            dft_mpi::simulate(ranks, |comm| {
                let run = |_value| {
                    let config = RunConfig::from_pairs(&[]).unwrap();
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
                let base = run("off");
                let candidate = run("on");
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
                assert_eq!(
                    base.outcome.result.energy.total_hartree().to_bits(),
                    candidate.outcome.result.energy.total_hartree().to_bits()
                );
                assert_eq!(
                    base.outcome
                        .result
                        .density
                        .as_slice()
                        .unwrap()
                        .iter()
                        .map(|v| v.to_bits())
                        .collect::<Vec<_>>(),
                    candidate
                        .outcome
                        .result
                        .density
                        .as_slice()
                        .unwrap()
                        .iter()
                        .map(|v| v.to_bits())
                        .collect::<Vec<_>>()
                );
                assert_eq!(base.records.len(), candidate.records.len());
                for (a, b) in base.records.iter().zip(&candidate.records) {
                    assert_eq!(a.dr2_ry.to_bits(), b.dr2_ry.to_bits());
                    assert_eq!(a.h_columns, b.h_columns);
                    assert_eq!(a.calls.len(), b.calls.len());
                }
                if comm.rank() == 0 {
                    eprintln!(
                        "P3-SETUP4 verified P={ranks}: delta_E=0 Ha density_bits_equal=true maps={}",
                        candidate.records.len()
                    );
                }
            });
        }
    }
}

#[test]
fn permap_preserves_delivery_and_counts() {
    use dft_engine::refscf::driver::run_with;
    use dft_engine::{ScfOptions, SpeciesPotential, Verification, refstart::StartInput};
    use dft_numerics::switches::RunConfig;
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    for rows in [
        [[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]],
        [[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]],
        [[5., 0., 0.], [0., 5., 0.], [0., 0., 12.]],
    ] {
        let cell = dft_geometry::Cell::from_user_rows_bohr(rows).unwrap();
        for (verification, max_iterations, _energy_rule) in [
            (Verification::Off, 100, "off"),
            (Verification::FixedPoint, 100, "off"),
            (Verification::FixedPoint, 100, "on"),
            (Verification::Off, 1, "off"),
        ] {
            dft_mpi::simulate(2, |comm| {
                let options = ScfOptions {
                    cutoff_hartree: 1.,
                    bands: Some(3),
                    gamma_real: dft_engine::GammaReal::On,
                    max_iterations,
                    verification,
                    ..Default::default()
                };
                let run = |permap| {
                    let config = RunConfig::from_pairs(&[("refscf_permap", permap)]).unwrap();
                    let _scope = config.enter();
                    run_with(
                        StartInput {
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
                let a = run("off");
                let b = run("on");
                assert_eq!(a.first_reference_stop, b.first_reference_stop);
                assert_eq!(a.verification_calls.len(), b.verification_calls.len());
                assert_eq!(a.records.len(), b.records.len());
                for (a, b) in a.records.iter().zip(&b.records) {
                    assert_eq!(a.descf_ha.to_bits(), b.descf_ha.to_bits());
                    assert_eq!(
                        a.variational_energy_ha.to_bits(),
                        b.variational_energy_ha.to_bits()
                    );
                    assert_eq!(a.dr2_ry.to_bits(), b.dr2_ry.to_bits());
                    if !b.energy.kinetic.is_finite() {
                        assert!(b.energy.band_sum.is_finite());
                        assert!(b.energy.double_counting.is_finite());
                        assert!(b.energy.functional_difference_bound.is_nan());
                    }
                }
                let ar = &a.outcome.result;
                let br = &b.outcome.result;
                assert_eq!(a.outcome.converged, b.outcome.converged);
                assert_eq!(
                    ar.energy.total_hartree().to_bits(),
                    br.energy.total_hartree().to_bits()
                );
                assert_eq!(
                    ar.density.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                    br.density.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
                );
                assert_eq!(ar.history.len(), br.history.len());
                for (a, b) in ar.history.iter().zip(&br.history) {
                    assert_eq!(a.davidson_iterations, b.davidson_iterations);
                    assert_eq!(a.operator_applications, b.operator_applications);
                    assert_eq!(a.standard_error_ha, b.standard_error_ha);
                }
                assert!(
                    br.history
                        .last()
                        .unwrap()
                        .energy
                        .total_hartree()
                        .is_finite()
                );
            });
        }
    }
}

#[test]
fn permap_gram_cache_ring_and_metric_invalidation() {
    use dft_numerics::switches::RunConfig;
    use std::cell::RefCell;
    struct Counted(RefCell<Vec<usize>>);
    impl Collective for Counted {
        fn rank(&self) -> usize {
            0
        }
        fn size(&self) -> usize {
            1
        }
        fn barrier(&self) -> Result<(), dft_mpi::MpiError> {
            Ok(())
        }
        fn allgatherv_f64(
            &self,
            a: &[f64],
            n: &[usize],
            b: &mut [f64],
        ) -> Result<(), dft_mpi::MpiError> {
            self.0.borrow_mut().push(a.len());
            Serial.allgatherv_f64(a, n, b)
        }
        fn tree_sum_f64(&self, a: &[f64], b: &mut [f64]) -> Result<(), dft_mpi::MpiError> {
            Serial.tree_sum_f64(a, b)
        }
        fn bcast_f64(&self, a: &mut [f64], root: usize) -> Result<(), dft_mpi::MpiError> {
            Serial.bcast_f64(a, root)
        }
        fn alltoall_f64(
            &self,
            a: &[f64],
            n: usize,
            b: &mut [f64],
        ) -> Result<(), dft_mpi::MpiError> {
            Serial.alltoall_f64(a, n, b)
        }
    }
    let run = |arm| {
        let cfg = RunConfig::from_pairs(&[("refscf_permap", arm)]).unwrap();
        let _scope = cfg.enter();
        let comm = Counted(RefCell::new(Vec::new()));
        let g2: Vec<_> = (0..24).map(|i| i as f64).collect();
        let smooth: Vec<_> = (0..24).map(|i| i < 22).collect();
        let mut mix = Broyden::new(
            Profile {
                history: 4,
                ..Default::default()
            },
            &g2,
            &smooth,
            10.,
            2.,
        )
        .unwrap();
        let mut m = Metric::new(10., &g2, &smooth, 1, 350).unwrap();
        let mut outputs = Vec::new();
        for k in 0..12 {
            if k == 7 {
                // New geometry, same vector lengths: cached Gram is invalid.
                m = Metric::new(11., &g2, &smooth, 1, 350).unwrap();
            }
            let x: Vec<_> = (0..24)
                .map(|i| C::new(((i + 1) * (k + 1)) as f64 * 0.03, (i as f64 * 0.2).cos()))
                .collect();
            let f: Vec<_> = (0..24)
                .map(|i| C::new(((i + 1) * (k + 1)) as f64, 0.2).sin())
                .collect();
            let mut next = x.clone();
            let d = mix.update(&comm, &mut m, &x, &f, &mut next).unwrap();
            outputs.push((next, mix.gamma().to_vec(), d.condition));
        }
        (outputs, comm.0.into_inner())
    };
    let (a, old) = run("off");
    let (b, new) = run("on");
    assert_eq!(
        a, b,
        "identical arithmetic including ring overwrite and changed metric"
    );
    assert_eq!(old[3], 14); // p=4: ten Gram forms plus four RHS forms.
    assert_eq!(new[3], 8); // Only the newest Gram row and the RHS.
    assert_eq!(new[6], 14); // Replaced metric forces a full refresh.
    assert!(new.iter().sum::<usize>() < old.iter().sum::<usize>());
    eprintln!(
        "P3-PERMAP Gram forms: {} -> {} (12 updates, history 4, one metric replacement)",
        old.iter().sum::<usize>(),
        new.iter().sum::<usize>()
    );
}

/// P3-PERMAP regression for dev08 0afbe868: the default shared rustfft
/// constructor must not bypass calibration. Inspect the computation's frozen
/// manifest BEFORE any JSON serializer could mark a switch as consulted.
#[test]
fn reference_scf_records_calibrated_fft_choices() {
    use dft_numerics::{
        fused_pair_plan::kernel_choices,
        switches::{Id, RunConfig},
    };
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
    let potentials = [dft_engine::SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    for rows in [
        [[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]],
        [[5., 0., 0.], [0.2, 5., 0.], [0., 0.1, 12.]],
    ] {
        let cell = dft_geometry::Cell::from_user_rows_bohr(rows).unwrap();
        for ranks in [1, 3, 8] {
            dft_mpi::simulate(ranks, |comm| {
                let options = dft_engine::ScfOptions {
                    cutoff_hartree: 1.,
                    bands: Some(3),
                    max_iterations: 100,
                    gamma_real: dft_engine::GammaReal::On,
                    verification: dft_engine::Verification::Off,
                    ..Default::default()
                };
                let run = |adaptive, _share| {
                    let config = RunConfig::from_pairs(&[
                        ("fft_kernel_by_length", adaptive),
                        ("m2_fft", "rust"),
                        ("fft_fused_rows", "off"),
                    ])
                    .unwrap();
                    let _scope = config.enter();
                    let output = dft_engine::refscf::driver::run_with(
                        dft_engine::refstart::StartInput {
                            cell: &cell,
                            atoms: &atoms,
                            potentials: &potentials,
                            options: &options,
                        },
                        Profile::default(),
                        comm,
                    )
                    .unwrap();
                    config.freeze_consulted();
                    if adaptive == "on" {
                        assert!(
                            config.consulted(Id::FFT_KERNEL_BY_LENGTH),
                            "reference computation bypassed the selector"
                        );
                        let choices = kernel_choices();
                        assert!(
                            !choices.is_empty(),
                            "rank {} has no kernel choices",
                            comm.rank()
                        );
                        for inverse in [false, true] {
                            assert!(choices.iter().any(|c| c.inverse == inverse && c.batch > 0));
                        }
                        assert!(choices.iter().all(|c| c.length > 0 && c.batch > 0));
                        if dft_numerics::MKL_BACKEND {
                            assert!(
                                choices
                                    .iter()
                                    .all(|c| c.rust_seconds.is_some() && c.mkl_seconds.is_some())
                            );
                        } else {
                            assert!(
                                choices
                                    .iter()
                                    .all(|c| !c.vendor && c.calibration_seconds == 0.)
                            );
                        }
                    }
                    assert!(output.outcome.converged);
                    output
                };
                let base = run("off", "on");
                for share in ["on", "off"] {
                    let candidate = run("on", share);
                    let a = &base.outcome.result;
                    let b = &candidate.outcome.result;
                    assert!((a.energy.total_hartree() - b.energy.total_hartree()).abs() < 5e-9);
                    if !dft_numerics::MKL_BACKEND {
                        assert_eq!(
                            a.energy.total_hartree().to_bits(),
                            b.energy.total_hartree().to_bits()
                        );
                        assert_eq!(
                            a.density.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                            b.density.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
                        );
                        assert_eq!(base.records.len(), candidate.records.len());
                    }
                }
            });
        }
    }
}

/// Start fields and H use the shared constructor; forks used by the density
/// transform must retain the selected kernels without replanning/calibrating.
#[test]
fn calibrated_shared_fft_keeps_choices_and_workspace_reuse() {
    use dft_numerics::{
        fused_pair_plan::kernel_choices,
        switches::{Id, RunConfig},
    };
    use dft_solvers::distributed_fft::{DistributedPairFft, PlanCache};
    for shape in [[4, 4, 4], [5, 7, 11], [3, 5, 24]] {
        for ranks in [1, 3, 8] {
            dft_mpi::simulate(ranks, |comm| {
                let config = RunConfig::from_pairs(&[
                    ("fft_kernel_by_length", "on"),
                    ("m2_fft", "rust"),
                    ("fft_fused_rows", "on"),
                ])
                .unwrap();
                let _config = config.enter();
                let cache = PlanCache::default();
                let _plans = cache.enter();
                let indices = vec![[0, 0, 0], [1, 1, 0], [2, 0, 1]];
                let mut original =
                    DistributedPairFft::new(shape, &indices, comm.rank(), ranks).unwrap();
                assert!(config.consulted(Id::FFT_KERNEL_BY_LENGTH));
                let choices = kernel_choices();
                assert!(!choices.is_empty());
                let counts = cache.counts();
                let mut fork = original.shared_workspace().expect("shared cache retained");
                let _repeated =
                    DistributedPairFft::new(shape, &indices, comm.rank(), ranks).unwrap();
                assert_eq!(cache.counts(), counts, "no new shared batch plans");
                let after = kernel_choices();
                assert_eq!(choices.len(), after.len(), "no repeated calibration");
                for (a, b) in choices.iter().zip(&after) {
                    assert_eq!(
                        (
                            a.vendor,
                            a.rust_seconds,
                            a.mkl_seconds,
                            a.calibration_seconds
                        ),
                        (
                            b.vendor,
                            b.rust_seconds,
                            b.mkl_seconds,
                            b.calibration_seconds
                        )
                    );
                }
                let input: Vec<_> = original
                    .local_positions()
                    .iter()
                    .map(|&i| C::new((i as f64 + 0.3).sin(), (i as f64 + 0.2).cos()))
                    .collect();
                let mut a = input.clone();
                let mut b = input.clone();
                let potential = vec![1.; shape.iter().product()];
                original
                    .apply_local_batch(comm, &mut a, &potential, 1)
                    .unwrap();
                fork.apply_local_batch(comm, &mut b, &potential, 1).unwrap();
                let bits = |v: &[C]| {
                    v.iter()
                        .map(|z| (z.re.to_bits(), z.im.to_bits()))
                        .collect::<Vec<_>>()
                };
                assert_eq!(
                    bits(&a),
                    bits(&b),
                    "same selected kernels in cloned density workspace"
                );
                for (actual, expected) in a.iter().zip(&input) {
                    assert!((*actual - expected).norm() < 2e-14);
                }
            });
        }
    }
}

fn ladder_fixture(label: &str, rows: [[f64; 3]; 3], nonlocal: bool) {
    use dft_engine::{GammaReal, ScfOptions, ScfProblem, SpeciesPotential, Verification};
    use dft_numerics::switches::{Id, RunConfig};
    assert!(!RunConfig::from_pairs(&[]).unwrap().on(Id::REFSCF_LADDER));
    let pp = common::potential_with_atomic_start(2., 1.8, nonlocal, 161, "SLA PZ");
    let cell = dft_geometry::Cell::from_user_rows_bohr(rows).unwrap();
    let atoms = [
        dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32]),
        dft_geometry::Atom::new("Xx", [0.57, 0.68, 0.71]),
    ];
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: 2.,
        bands: Some(5),
        max_iterations: 120,
        gamma_real: GammaReal::On,
        ..Default::default()
    };
    let strict_options = ScfOptions {
        // Independent deterministic random start avoids the legacy atomic
        // seed's incomplete occupied subspace on synthetic nonlocal projectors.
        free_electron_initial_guess: false,
        start_wavefunctions: dft_engine::StartWavefunctions::FreeElectron,
        cutoff_hartree: 2.,
        bands: Some(5),
        max_iterations: 120,
        gamma_real: GammaReal::On,
        ..ScfOptions::strict_reference()
    };
    let strict = ScfProblem::new(&cell, &atoms, &potentials, strict_options)
        .unwrap()
        .run()
        .unwrap();
    assert!(strict.certified);
    let reference = strict.energy.total_hartree();
    let mut rank_results = Vec::new();
    for ranks in [1, 3, 8] {
        let results = dft_mpi::simulate(ranks, |comm| {
            let run = |arm, threshold_ry, verify| {
                let config = RunConfig::from_pairs(&[("refscf_ladder", arm)]).unwrap();
                let _scope = config.enter();
                let opts = ScfOptions {
                    verification: verify,
                    ..options.clone()
                };
                let run = dft_engine::refscf::driver::run_with(
                    dft_engine::refstart::StartInput {
                        cell: &cell,
                        atoms: &atoms,
                        potentials: &potentials,
                        options: &opts,
                    },
                    Profile {
                        threshold_ry,
                        ..Default::default()
                    },
                    comm,
                )
                .unwrap();
                assert!(
                    run.outcome.converged,
                    "{label} P={ranks} {arm}: {}",
                    run.exit_reason
                );
                if verify == Verification::FixedPoint {
                    assert!(run.outcome.result.verification.as_ref().unwrap().passed());
                } else {
                    assert!(run.outcome.result.verification.is_none());
                }
                assert!(run.records.iter().all(|r| r.stage == "fine"));
                assert_eq!(run.records[0].mixing.occupancy, 0);
                assert_eq!(run.records[0].ethr_requested_ry, 1e-2);
                let coarse = run.ladder.as_ref().map_or(0, |l| {
                    assert!(
                        l.coarse_records
                            .iter()
                            .all(|r| r.stage == "coarse" && r.verification_calls.is_empty())
                    );
                    assert_eq!(l.coarse_threshold_ry, 100. * threshold_ry);
                    assert!(run.records[0].map_start_s >= l.coarse_seconds);
                    l.coarse_records.len()
                });
                assert_eq!(run.ladder.is_some(), arm == "on");
                let e = run.outcome.result.energy.total_hartree();
                assert!(
                    (e - reference).abs() < 5e-9,
                    "{label} {arm} P={ranks} threshold={threshold_ry:e} E-Estrict={:e}",
                    e - reference
                );
                if comm.rank() == 0 {
                    eprintln!(
                        "P3-LADDER {label} P={ranks} arm={arm} threshold={threshold_ry:e} Ry verify={verify:?} coarse={coarse} fine={} E={e:.17e} Ha E-Estrict={:.9e} Ha",
                        run.records.len(),
                        e - reference
                    );
                }
                (e, run.records.len())
            };
            let d = run("off", 1e-12, Verification::FixedPoint);
            let l = run("on", 1e-12, Verification::FixedPoint);
            assert!(
                (d.0 - l.0).abs() < 1e-10,
                "{label}: tight L-D={:e} Ha",
                l.0 - d.0
            );
            run("on", 1e-8, Verification::FixedPoint);
            run("on", 1e-8, Verification::Off);
            (d, l)
        });
        for r in &results[1..] {
            assert_eq!(r, &results[0]);
        }
        rank_results.push(results[0]);
    }
    for r in &rank_results[1..] {
        assert!((r.0.0 - rank_results[0].0.0).abs() < 1e-10);
        assert!((r.1.0 - rank_results[0].1.0).abs() < 1e-10);
    }
}
#[test]
fn ladder_verified_cubic() {
    ladder_fixture("cubic", [[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]], false);
}
#[test]
fn ladder_verified_skew() {
    ladder_fixture("skew", [[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]], false);
}
#[test]
fn ladder_verified_elongated() {
    ladder_fixture(
        "elongated",
        [[5., 0., 0.], [0., 5., 0.], [0., 0., 10.]],
        false,
    );
}
#[test]
fn ladder_verified_nonlocal() {
    ladder_fixture("nonlocal", [[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]], true);
}

#[test]
fn fast_order_kept_hx_delta_and_invalidation() {
    use dft_engine::{GammaReal, ScfOptions, SpeciesPotential, refscf::kept_hx::KeptHx, refstart};
    use dft_numerics::{RealMatrix, switches::RunConfig};
    use dft_solvers::ref_davidson::{self, Options};
    // default on since 2026-11-26
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    let cell =
        dft_geometry::Cell::from_user_rows_bohr([[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]])
            .unwrap();
    let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: 1.,
        bands: Some(3),
        gamma_real: GammaReal::On,
        ..Default::default()
    };
    dft_mpi::simulate(2, |comm| {
        let _scope = RunConfig::from_pairs(&[]).unwrap().enter();
        let mut state = refstart::start(
            refstart::StartInput {
                cell: &cell,
                atoms: &atoms,
                potentials: &potentials,
                options: &options,
            },
            comm,
        )
        .unwrap();
        let x = state.orbitals.clone();
        let b = x.cols();
        let offsets = state.operator.plane_grid_offsets();
        let old: Vec<_> = offsets.iter().map(|&g| 0.03 * (g as f64).sin()).collect();
        let new: Vec<_> = offsets
            .iter()
            .map(|&g| 0.04 * (0.7 * g as f64).cos())
            .collect();
        state.operator.set_owned_potential(&old, 1).unwrap();
        let mut kept = KeptHx::new(true);
        assert!(
            kept.prepare(comm, &mut state.operator, &x, &old)
                .unwrap()
                .is_none()
        );
        let mut hx = RealMatrix::zeros(x.rows(), b).unwrap();
        state.operator.apply_block(comm, &x, &mut hx).unwrap();
        let mut density = vec![0.; old.len()];
        kept.density(
            comm,
            &mut state.operator,
            &x,
            &hx,
            &old,
            &[2., 0., 0.],
            cell.volume_bohr3(),
            &mut density,
        )
        .unwrap();
        state.operator.set_owned_potential(&new, 2).unwrap();
        let reused = kept
            .prepare(comm, &mut state.operator, &x, &new)
            .unwrap()
            .unwrap();
        let mut direct = RealMatrix::zeros(x.rows(), b).unwrap();
        state.operator.apply_block(comm, &x, &mut direct).unwrap();
        let local = [
            reused
                .as_slice()
                .iter()
                .zip(direct.as_slice())
                .map(|(a, b)| (a - b).powi(2))
                .sum(),
            direct.as_slice().iter().map(|v| v * v).sum(),
        ];
        let mut total = [0.; 2];
        comm.tree_sum_f64(&local, &mut total).unwrap();
        assert!(total[1] > 0.);
        let relative = (total[0] / total[1]).sqrt();
        assert!(relative <= 1e-13, "relative Frobenius error={relative:e}");
        let mut opts = Options::new(b, 1e-12, b);
        opts.max_iterations = 1;
        let fresh =
            ref_davidson::solve(comm, &mut state.operator, &x, &state.diagonal_static, opts)
                .unwrap();
        let hit = ref_davidson::solve_with_initial_hx(
            comm,
            &mut state.operator,
            &x,
            &state.diagonal_static,
            opts,
            Some(&reused),
        )
        .unwrap();
        assert_eq!(fresh.record.setup.h_columns, b);
        assert_eq!(hit.record.setup.h_columns, 0);
        assert_eq!(fresh.record.h_columns - hit.record.h_columns, b);
        let report = kept.report(comm).unwrap();
        assert_eq!(report.cache_hits, 1);
        assert_eq!(report.full_h_starting_columns_replaced, b);
        assert_eq!(report.delta_v_columns, b);
        assert!(
            report
                .retained_slab_bytes_per_rank
                .iter()
                .all(|&n| n > 0 && n <= dft_engine::refscf::kept_hx::SLAB_LIMIT_BYTES)
        );
        // A change on only one rank invalidates the entire distributed block.
        kept.density(
            comm,
            &mut state.operator,
            &x,
            &direct,
            &new,
            &[2., 0., 0.],
            cell.volume_bohr3(),
            &mut density,
        )
        .unwrap();
        let mut changed = x.clone();
        if comm.rank() == 0 {
            changed[(0, 0)] += 1e-8;
        }
        state.operator.set_owned_potential(&old, 3).unwrap();
        let invalid = kept
            .prepare(comm, &mut state.operator, &changed, &old)
            .unwrap();
        assert!(invalid.is_none());
        let miss = ref_davidson::solve_with_initial_hx(
            comm,
            &mut state.operator,
            &changed,
            &state.diagonal_static,
            opts,
            invalid.as_ref(),
        )
        .unwrap();
        assert_eq!(miss.record.setup.h_columns, b);
        assert_eq!(kept.report(comm).unwrap().cache_hits, 1);
        // An equal-shaped rebuilt operator must not accept the old token.
        state.operator.set_owned_potential(&new, 4).unwrap();
        kept.density(
            comm,
            &mut state.operator,
            &x,
            &direct,
            &new,
            &[2., 0., 0.],
            cell.volume_bohr3(),
            &mut density,
        )
        .unwrap();
        let mut rebuilt = refstart::start(
            refstart::StartInput {
                cell: &cell,
                atoms: &atoms,
                potentials: &potentials,
                options: &options,
            },
            comm,
        )
        .unwrap();
        rebuilt.operator.set_owned_potential(&old, 5).unwrap();
        assert!(
            kept.prepare(comm, &mut rebuilt.operator, &x, &old)
                .unwrap()
                .is_none()
        );
        // Rank-local invalidation also forces a collective miss everywhere.
        kept.density(
            comm,
            &mut state.operator,
            &x,
            &direct,
            &new,
            &[2., 0., 0.],
            cell.volume_bohr3(),
            &mut density,
        )
        .unwrap();
        if comm.rank() == 1 {
            kept.invalidate();
        }
        state.operator.set_owned_potential(&old, 5).unwrap();
        assert!(
            kept.prepare(comm, &mut state.operator, &x, &old)
                .unwrap()
                .is_none()
        );
        // Supplied products cannot bypass rank or finite-input checks.
        let zero = RealMatrix::zeros(x.rows(), b).unwrap();
        assert!(
            ref_davidson::solve_with_initial_hx(
                comm,
                &mut state.operator,
                &zero,
                &state.diagonal_static,
                opts,
                Some(&zero)
            )
            .is_err()
        );
        let mut invalid_hx = reused.clone();
        if comm.rank() == 0 {
            invalid_hx[(0, 0)] = f64::NAN;
        }
        assert!(
            ref_davidson::solve_with_initial_hx(
                comm,
                &mut state.operator,
                &x,
                &state.diagonal_static,
                opts,
                Some(&invalid_hx)
            )
            .is_err()
        );
        assert!(
            ref_davidson::solve_with_initial_hx(
                comm,
                &mut state.operator,
                &x,
                &state.diagonal_static,
                opts,
                (comm.rank() == 0).then_some(&reused)
            )
            .is_err()
        );
        if comm.rank() == 0 {
            eprintln!(
                "KEPT_HX relative_Frobenius={relative:.17e} hit_start_H=0 fresh_start_H={b} invalid_start_H={b} replaced={b} delta_V={b} slab_bytes={:?}",
                report.retained_slab_bytes_per_rank
            );
        }
    });
}

#[test]
fn fast_order_ladder_exact_injection_and_verified_delivery() {
    use dft_engine::{
        GammaReal, ScfOptions, SpeciesPotential, Verification, refscf::ladder::Seed, refstart,
    };
    use dft_numerics::{RealMatrix, switches::RunConfig};
    let pp = common::potential_with_atomic_start(2., 1.8, false, 161, "SLA PZ");
    let cell =
        dft_geometry::Cell::from_user_rows_bohr([[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]])
            .unwrap();
    let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let fine_options = ScfOptions {
        cutoff_hartree: 2.,
        bands: Some(3),
        max_iterations: 120,
        gamma_real: GammaReal::On,
        verification: Verification::FixedPoint,
        ..Default::default()
    };
    let coarse_options = ScfOptions {
        cutoff_hartree: 1.,
        ..fine_options.clone()
    };
    dft_mpi::simulate(2, |comm| {
        let _scope = RunConfig::from_pairs(&[]).unwrap().enter();
        let input = |options| refstart::StartInput {
            cell: &cell,
            atoms: &atoms,
            potentials: &potentials,
            options,
        };
        let mut coarse = refstart::start(input(&coarse_options), comm).unwrap();
        let mut fine = refstart::start(input(&fine_options), comm).unwrap();
        let gather = |n: usize, rows: &[usize], x: &RealMatrix| {
            // Transport without addition so even signed-zero bits survive.
            let width = n * (x.cols() + 1);
            let mut local = vec![0.; width];
            for (i, &r) in rows.iter().enumerate() {
                local[r] = 1.;
                for j in 0..x.cols() {
                    local[n + j * n + r] = x[(i, j)];
                }
            }
            let mut packets = vec![0.; width * comm.size()];
            comm.allgatherv_f64(&local, &vec![width; comm.size()], &mut packets)
                .unwrap();
            let mut all = vec![0.; n * x.cols()];
            for r in 0..n {
                let owners: Vec<_> = packets.chunks_exact(width).filter(|p| p[r] == 1.).collect();
                assert_eq!(owners.len(), 1);
                for j in 0..x.cols() {
                    all[j * n + r] = owners[0][n + j * n + r];
                }
            }
            all
        };
        let nc = coarse.basis.vectors().len();
        let nf = fine.basis.vectors().len();
        let before = gather(nc, coarse.operator.local_rows(), &coarse.orbitals);
        let mapping: Vec<_> = coarse
            .basis
            .vectors()
            .iter()
            .map(|g| {
                fine.basis
                    .vectors()
                    .iter()
                    .position(|f| f.miller == g.miller)
                    .unwrap()
            })
            .collect();
        let mut orc_map = vec![0; nc];
        for p in 1..coarse.half.plus().len() {
            let g = mapping[coarse.half.plus()[p] as usize];
            let q = fine
                .half
                .plus()
                .iter()
                .position(|&v| v as usize == g)
                .unwrap();
            assert_eq!(
                mapping[coarse.half.minus()[p] as usize],
                fine.half.minus()[q] as usize
            );
            orc_map[2 * p - 1] = 2 * q - 1;
            orc_map[2 * p] = 2 * q;
        }
        let dims = coarse.grid.dims();
        let density: Vec<_> = coarse
            .fields
            .plane_offsets()
            .iter()
            .map(|&g| {
                (2. + 0.1 * (TAU * (g / (dims[1] * dims[2])) as f64 / dims[0] as f64).cos())
                    / cell.volume_bohr3()
            })
            .collect();
        let density_g = coarse.fields.forward(comm, &density).unwrap();
        let seed = Seed {
            basis: coarse.basis,
            grid: coarse.grid,
            half: coarse.half,
            rows: coarse.operator.local_rows().to_vec(),
            orbitals: coarse.orbitals,
            indices: coarse.fields.reciprocal_indices().to_vec(),
            density_g,
        };
        let delivery_seed = seed.clone();
        let (g, x) = seed
            .inject(
                &fine.basis,
                &fine.grid,
                &fine.half,
                fine.operator.local_rows(),
                &fine.fields,
                comm,
            )
            .unwrap();
        // Inspect raw injected ORC values before any fine-stage QR or eigensolve.
        let after = gather(nf, fine.operator.local_rows(), &x);
        let mut expected = vec![0f64; nf * 3];
        for b in 0..3 {
            for (i, &r) in orc_map.iter().enumerate() {
                expected[b * nf + r] = before[b * nc + i];
            }
        }
        for (a, b) in after.iter().zip(&expected) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
        let mut gram_change = 0f64;
        for i in 0..3 {
            for j in 0..3 {
                let a: f64 = (0..nc)
                    .map(|r| before[i * nc + r] * before[j * nc + r])
                    .sum();
                let b: f64 = (0..nf).map(|r| after[i * nf + r] * after[j * nf + r]).sum();
                gram_change = gram_change.max((a - b).abs());
            }
        }
        assert!(gram_change <= 1e-12);
        for b in 0..3 {
            let mut complex = vec![C::default(); nf];
            fine.half
                .unpack(&after[b * nf..(b + 1) * nf], &mut complex)
                .unwrap();
            for (&p, &m) in fine.half.plus().iter().zip(fine.half.minus()).skip(1) {
                assert_eq!(
                    complex[p as usize].re.to_bits(),
                    complex[m as usize].re.to_bits()
                );
                assert_eq!(
                    (-complex[p as usize].im).to_bits(),
                    complex[m as usize].im.to_bits()
                );
            }
        }
        let planes = fine.fields.inverse(comm, &g).unwrap();
        let local = planes.iter().sum::<f64>() * cell.volume_bohr3()
            / fine.grid.dims().iter().product::<usize>() as f64;
        let mut charge = [0.];
        comm.tree_sum_f64(&[local], &mut charge).unwrap();
        let charge_error = (charge[0] - 2.).abs();
        assert!(charge_error <= 1e-12);
        let run = |ladder| {
            let _scope = RunConfig::from_pairs(&[("refscf_ladder", ladder)])
                .unwrap()
                .enter();
            let run = dft_engine::refscf::driver::run_with(
                input(&fine_options),
                Profile {
                    threshold_ry: 1e-12,
                    ..Default::default()
                },
                comm,
            )
            .unwrap();
            assert!(run.outcome.converged, "{}", run.exit_reason);
            assert!(run.outcome.result.verification.as_ref().unwrap().passed());
            assert!(run.kept_hx.cache_hits > 0);
            assert_eq!(
                run.kept_hx.full_h_starting_columns_replaced,
                3 * run.kept_hx.cache_hits
            );
            run.outcome.result.energy.total_hartree()
        };
        let direct = run("off");
        let ladder = run("on");
        let seeded = {
            let _scope = RunConfig::from_pairs(&[]).unwrap().enter();
            let run = dft_engine::refscf::driver::run_seeded_with(
                input(&fine_options),
                Profile {
                    threshold_ry: 1e-12,
                    ..Default::default()
                },
                comm,
                delivery_seed,
            )
            .unwrap();
            assert!(run.outcome.converged, "{}", run.exit_reason);
            assert!(run.outcome.result.verification.as_ref().unwrap().passed());
            assert!(run.kept_hx.cache_hits > 0);
            run.outcome.result.energy.total_hartree()
        };
        let seeded_delta = (direct - seeded).abs();
        assert!(
            seeded_delta <= 1e-10,
            "inspected seed-direct={seeded_delta:e} Ha"
        );
        let delta = (direct - ladder).abs();
        assert!(delta <= 1e-10, "ladder-direct={delta:e} Ha");
        if comm.rank() == 0 {
            eprintln!(
                "LADDER copied_and_zero_bits=identical charge_error={charge_error:.17e} e max_Gram_change={gram_change:.17e} delta_E={delta:.17e} Ha inspected_seed_delta_E={seeded_delta:.17e} Ha direct_E={direct:.17e} Ha ladder_E={ladder:.17e} Ha"
            );
        }
    });
}

/// Isolate the process-wide counter enable flag from parallel SCF tests.
#[test]
fn fast_order_keeping_density_fft_counts_and_bits() {
    const CHILD: &str = "BAANDI_TEST_KEPTHX2_COUNTER_CHILD";
    if std::env::var_os(CHILD).is_none() {
        let status = std::process::Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "fast_order_keeping_density_fft_counts_and_bits",
                "--nocapture",
            ])
            .env(CHILD, "1")
            .status()
            .unwrap();
        assert!(status.success());
        return;
    }
    use dft_engine::{GammaReal, ScfOptions, SpeciesPotential, refstart};
    use dft_numerics::{RealMatrix, stage_clock::local};
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    let cell =
        dft_geometry::Cell::from_user_rows_bohr([[7., 0., 0.], [0.2, 8., 0.], [0., 0.1, 6.]])
            .unwrap();
    let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let bits = |v: &[f64]| v.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
    for ranks in [1, 2, 5] {
        for bands in [3, 4] {
            let options = ScfOptions {
                cutoff_hartree: 1.,
                bands: Some(bands),
                gamma_real: GammaReal::On,
                ..Default::default()
            };
            dft_mpi::simulate(ranks, |comm| {
                let mut state = refstart::start(
                    refstart::StartInput {
                        cell: &cell,
                        atoms: &atoms,
                        potentials: &potentials,
                        options: &options,
                    },
                    comm,
                )
                .unwrap();
                local::enable(true);
                let x = &state.orbitals;
                let op = &mut state.operator;
                let points = op.plane_points();
                let weights: Vec<_> = (0..bands)
                    .map(|b| if b == 1 { 0. } else { 2. / (b + 1) as f64 })
                    .collect();
                let mut plain = vec![0.; points];
                let before = local::snapshot();
                op.band_density_sum(comm, x, &weights, cell.volume_bohr3(), &mut plain)
                    .unwrap();
                let plain_counts = local::snapshot().since(&before);
                assert_eq!(
                    plain_counts.calls[local::OP_COUNT + 3],
                    bands.div_ceil(2) as u64
                );
                assert_eq!(plain_counts.calls[local::OP_COUNT + 2], 0);
                // The unchanged non-keeping value path supplies an independent
                // slab/layout oracle, including zero-occupation buffer bands.
                let mut values = vec![0.; bands * points];
                op.band_values(comm, x, &mut values).unwrap();
                for pairs in [0, 1, bands / 2] {
                    let oracle: Vec<_> = (0..pairs)
                        .flat_map(|t| {
                            let values = &values;
                            (0..points).map(move |p| {
                                C::new(values[2 * t * points + p], values[(2 * t + 1) * points + p])
                            })
                        })
                        .collect();
                    let mut keep = vec![C::new(f64::NAN, f64::NAN); pairs * points];
                    let mut density = vec![f64::NAN; points];
                    let before = local::snapshot();
                    op.band_density_sum_keeping(
                        comm,
                        x,
                        &weights,
                        cell.volume_bohr3(),
                        &mut density,
                        &mut keep,
                    )
                    .unwrap();
                    let counts = local::snapshot().since(&before);
                    for slot in [local::OP_COUNT + 2, local::OP_COUNT + 3] {
                        assert_eq!(counts.calls[slot], plain_counts.calls[slot]);
                    }
                    assert_eq!(bits(&density), bits(&plain));
                    for (a, b) in keep.iter().zip(&oracle) {
                        assert_eq!(
                            (a.re.to_bits(), a.im.to_bits()),
                            (b.re.to_bits(), b.im.to_bits())
                        );
                    }
                    let delta: Vec<_> = (0..op.potential().len())
                        .map(|g| 0.02 * (g as f64).sin())
                        .collect();
                    let mut kept_change = RealMatrix::zeros(x.rows(), bands).unwrap();
                    let mut old_change = kept_change.clone();
                    op.apply_local_change(comm, x, &keep, pairs, &delta, &mut kept_change)
                        .unwrap();
                    op.apply_local_change(comm, x, &oracle, pairs, &delta, &mut old_change)
                        .unwrap();
                    assert_eq!(bits(kept_change.as_slice()), bits(old_change.as_slice()));
                }
            });
            eprintln!(
                "KEPTHX2 ranks={ranks} bands={bands}: plain/keeping inverse={} forward=0; density/slab/delta bits identical",
                bands.div_ceil(2)
            );
        }
    }
}

#[test]
fn fast_order_kept_hx_delivery_fingerprint() {
    use dft_engine::{GammaReal, ScfOptions, SpeciesPotential, refstart};
    use dft_numerics::switches::RunConfig;
    use std::hash::{Hash, Hasher};
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    let cell =
        dft_geometry::Cell::from_user_rows_bohr([[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]])
            .unwrap();
    let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: 1.,
        bands: Some(3),
        gamma_real: GammaReal::On,
        ..Default::default()
    };
    dft_mpi::simulate(2, |comm| {
        let _scope = RunConfig::from_pairs(&[]).unwrap().enter();
        let mut previous = None;
        for _ in 0..2 {
            let run = dft_engine::refscf::driver::run_with(
                refstart::StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &options,
                },
                Profile {
                    threshold_ry: 1e-12,
                    ..Default::default()
                },
                comm,
            )
            .unwrap();
            assert!(run.outcome.converged, "{}", run.exit_reason);
            let result = &run.outcome.result;
            assert!(result.verification.as_ref().unwrap().passed());
            assert!(run.kept_hx.cache_hits > 0);
            let mut hash = std::collections::hash_map::DefaultHasher::new();
            for v in result
                .density
                .iter()
                .chain(&result.eigenvalues)
                .chain(&result.occupations)
            {
                v.to_bits().hash(&mut hash);
            }
            for v in result.orbitals.as_slice() {
                v.re.to_bits().hash(&mut hash);
                v.im.to_bits().hash(&mut hash);
            }
            result.total_hartree().to_bits().hash(&mut hash);
            result.charge.to_bits().hash(&mut hash);
            run.records.len().hash(&mut hash);
            run.kept_hx.cache_hits.hash(&mut hash);
            run.kept_hx.full_h_starting_columns_replaced.hash(&mut hash);
            let fingerprint = hash.finish();
            if let Some(previous) = previous {
                assert_eq!(fingerprint, previous);
            } else if comm.rank() == 0 {
                eprintln!(
                    "KEPTHX2 delivery fingerprint={fingerprint:016x} energy={:.17e} Ha maps={} hits={} replaced={}",
                    result.total_hartree(),
                    run.records.len(),
                    run.kept_hx.cache_hits,
                    run.kept_hx.full_h_starting_columns_replaced
                );
            }
            previous = Some(fingerprint);
        }
    });
}

// P3-VERWARM: compare the complete delivered trajectory, including failed checks.
fn verify_warm_fixture(label: &str, rows: [[f64; 3]; 3], nonlocal: bool, recovery: bool) {
    use dft_engine::{GammaReal, ScfOptions, SpeciesPotential, refscf::driver, refstart};
    use dft_numerics::switches::RunConfig;

    let pp = common::potential_with_atomic_start(
        if recovery { 4. } else { 2. },
        if recovery { 1.2 } else { 1.8 },
        nonlocal,
        161,
        "SLA PZ",
    );
    let cell = dft_geometry::Cell::from_user_rows_bohr(rows).unwrap();
    let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: if recovery { 2. } else { 1. },
        bands: Some(3),
        max_iterations: 120,
        gamma_real: GammaReal::On,
        ..Default::default()
    };
    let profile = Profile {
        history: if recovery { 4 } else { 8 },
        ..Default::default()
    };
    let mut rank_energy: Option<f64> = None;
    for ranks in [1, 3, 8] {
        let outputs = dft_mpi::simulate(ranks, |comm| {
            let run = || {
                let _scope = RunConfig::from_pairs(&[]).unwrap().enter();
                driver::run_with(
                    refstart::StartInput {
                        cell: &cell,
                        atoms: &atoms,
                        potentials: &potentials,
                        options: &options,
                    },
                    profile,
                    comm,
                )
                .unwrap()
            };
            let first = run();
            let repeat = run();
            for r in [&first, &repeat] {
                assert!(r.outcome.converged, "{label}: {}", r.exit_reason);
                assert!(r.outcome.result.verification.as_ref().unwrap().passed());
                stop_to_verified_clocks(r);
                verify_cost_contract(r);
            }
            let attempts = |r: &driver::Run| {
                r.outcome
                    .result
                    .verification
                    .as_ref()
                    .unwrap()
                    .attempts
                    .iter()
                    .map(|a| (a.iteration, a.passed))
                    .collect::<Vec<_>>()
            };
            assert_eq!(attempts(&first), attempts(&repeat), "{label} ranks={ranks}");
            assert_eq!(first.first_reference_stop, repeat.first_reference_stop);
            assert_eq!(first.records.len(), repeat.records.len());
            let delta =
                repeat.outcome.result.total_hartree() - first.outcome.result.total_hartree();
            assert!(delta.abs() < 1e-12, "{label} ranks={ranks}: {delta:e} Ha");
            let density_delta = repeat
                .outcome
                .result
                .density
                .iter()
                .zip(&first.outcome.result.density)
                .map(|(a, b)| (a - b).abs())
                .fold(0., f64::max);
            assert!(
                density_delta < 1e-10,
                "density delta {density_delta:e} Bohr^-3"
            );
            if !recovery {
                assert_eq!(
                    repeat.outcome.result.total_hartree().to_bits(),
                    first.outcome.result.total_hartree().to_bits()
                );
                assert!(
                    repeat
                        .outcome
                        .result
                        .density
                        .iter()
                        .zip(&first.outcome.result.density)
                        .all(|(a, b)| a.to_bits() == b.to_bits())
                );
            }
            let cols = |r: &driver::Run| {
                r.verification_calls
                    .iter()
                    .map(|c| c.h_columns)
                    .sum::<usize>()
            };
            assert_eq!(cols(&repeat), cols(&first));
            assert!(
                first
                    .records
                    .iter()
                    .filter_map(|r| r.verification_warm)
                    .all(|v| v)
            );
            assert!(
                repeat
                    .records
                    .iter()
                    .filter_map(|r| r.verification_warm)
                    .all(|v| v)
            );
            for (d, w) in first
                .outcome
                .result
                .verification
                .as_ref()
                .unwrap()
                .attempts
                .iter()
                .zip(
                    &repeat
                        .outcome
                        .result
                        .verification
                        .as_ref()
                        .unwrap()
                        .attempts,
                )
            {
                assert_eq!(d.threshold_ha, w.threshold_ha);
                // The input estimator can differ by rounding after a recovery;
                // every attempt still uses exactly the same next-ratchet rule.
                for r in [&first, &repeat] {
                    let a = &r.outcome.result.verification.as_ref().unwrap().attempts;
                    let a = a.iter().find(|a| a.iteration == d.iteration).unwrap();
                    let m = &r.records[a.iteration - 1];
                    assert_eq!(
                        a.eigen_tolerance_ha,
                        next_ethr_ry(
                            m.ethr_ry,
                            m.dr2_ry,
                            r.outcome.result.occupations.iter().sum()
                        ) / 2.
                    );
                }
            }
            if recovery {
                let failed = attempts(&repeat)
                    .into_iter()
                    .find(|(_, passed)| !passed)
                    .unwrap()
                    .0;
                assert!(first.records[failed].warm_start);
                assert!(repeat.records[failed].warm_start);
                assert_eq!(
                    repeat.records[failed].h_columns, first.records[failed].h_columns,
                    "repeated recovery must retain identical H work"
                );
            }
            if comm.rank() == 0 {
                eprintln!(
                    "VERWARM {label} ranks={ranks} maps={} attempts={:?} verification_H={}->{} dE={delta:.17e} Ha max_drho={density_delta:.3e} Bohr^-3",
                    repeat.records.len(),
                    attempts(&repeat),
                    cols(&first),
                    cols(&repeat)
                );
            }
            repeat.outcome.result.total_hartree()
        });
        if let Some(e) = rank_energy {
            assert!((outputs[0] - e).abs() < 1e-12);
        } else {
            rank_energy = Some(outputs[0]);
        }
    }
}
#[test]
fn verify_warm_cubic() {
    verify_warm_fixture(
        "cubic",
        [[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]],
        false,
        false,
    );
}
#[test]
fn verify_warm_skew() {
    verify_warm_fixture(
        "skew",
        [[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]],
        false,
        false,
    );
}
#[test]
fn verify_warm_elongated() {
    verify_warm_fixture(
        "elongated",
        [[5., 0., 0.], [0., 5., 0.], [0., 0., 10.]],
        false,
        false,
    );
}
#[test]
fn verify_warm_nonlocal() {
    verify_warm_fixture(
        "nonlocal",
        [[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]],
        true,
        false,
    );
}
#[test]
fn verify_warm_failed_recovery() {
    verify_warm_fixture(
        "recovery",
        [[11., 0., 0.], [0.2, 11., 0.], [0., 0.1, 11.]],
        true,
        true,
    );
}

/// P4: generalized triangular solves plus exact construction reuse; all four
/// cell shapes, multi-species/projectors, odd tails and kinetic completion.
#[test]
fn p4_start_switches_preserve_fields_and_ritz() {
    use dft_engine::refstart::{StartInput, start};
    use dft_numerics::switches::RunConfig;
    // Default on since 2026-11-26 (dev08 P4 screen).

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
        let atoms = [
            dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32]),
            dft_geometry::Atom::new(if shape == 1 { "Yy" } else { "Xx" }, [0.57, 0.68, 0.71]),
        ];
        let potentials = [
            dft_engine::SpeciesPotential {
                label: "Xx",
                potential: &pp,
            },
            dft_engine::SpeciesPotential {
                label: "Yy",
                potential: if shape == 1 { &pp2 } else { &pp },
            },
        ];
        for ranks in [1, 3, 8] {
            dft_mpi::simulate(ranks, |comm| {
                let build = |_fft, _gram, _shells, threads| {
                    let config = RunConfig::from_pairs(&[]).unwrap();
                    let _scope = config.enter();
                    let options = dft_engine::ScfOptions {
                        cutoff_hartree: if shape == 0 { 8. } else { 2. },
                        bands: Some(if shape == 3 { 9 } else { 5 }),
                        gamma_real: dft_engine::GammaReal::On,
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
                        comm,
                    )
                    .unwrap()
                };
                let base = build("off", "off", "off", 1);
                if shape == 0 && ranks == 1 {
                    assert!(
                        base.orbitals.rows() > 256,
                        "exercise a multi-panel Gram reduction"
                    );
                }
                let fields = |s: &dft_engine::refstart::StartState| {
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
                        .collect::<Vec<_>>()
                };
                for (fft, gram, shells) in [
                    ("on", "off", "off"),
                    ("off", "off", "on"),
                    ("off", "on", "off"),
                    ("on", "on", "on"),
                ] {
                    let mut previous = None;
                    for threads in [1, 2] {
                        let mut candidate = build(fft, gram, shells, threads);
                        assert!(candidate.operator.batched());
                        if gram == "on" && shape != 3 {
                            assert!(
                                candidate
                                    .setup
                                    .phases
                                    .contains_key("atomic_overlap_factor_solve")
                            );
                            assert!(
                                !candidate
                                    .setup
                                    .phases
                                    .contains_key("atomic_overlap_factor_inverse")
                            );
                            assert!(
                                !candidate.setup.phases.contains_key("atomic_qr_gram"),
                                "candidate must exercise the generalized solve, not fallback"
                            );
                        }
                        assert_eq!(fields(&base), fields(&candidate));
                        let delta = base
                            .ritz_values_hartree
                            .iter()
                            .zip(&candidate.ritz_values_hartree)
                            .map(|(a, b)| (a - b).abs())
                            .fold(0., f64::max);
                        assert!(delta < 1e-12, "Ritz error {delta:e} Ha");
                        if comm.rank() == 0 && threads == 1 {
                            println!(
                                "P4 shape={shape} ranks={ranks} partition={fft} cholesky={gram} fill={shells}: Ritz_delta={delta:e} Ha"
                            );
                        }
                        let bands = candidate.orbitals.cols();
                        let mut dots = vec![0.; 2 * bands * bands];
                        for j in 0..bands {
                            for i in 0..bands {
                                for r in 0..candidate.orbitals.rows() {
                                    dots[i + j * bands] +=
                                        base.orbitals[(r, i)] * candidate.orbitals[(r, j)];
                                    dots[bands * bands + i + j * bands] +=
                                        candidate.orbitals[(r, i)] * candidate.orbitals[(r, j)];
                                }
                            }
                        }
                        let mut global = vec![0.; dots.len()];
                        comm.tree_sum_f64(&dots, &mut global).unwrap();
                        for j in 0..bands {
                            for i in 0..bands {
                                assert!(
                                    (global[bands * bands + i + j * bands] - f64::from(i == j))
                                        .abs()
                                        < 1e-12
                                );
                            }
                        }
                        assert!(
                            setup3_projection_residual(
                                comm,
                                &base.orbitals,
                                &candidate.orbitals,
                                &global[..bands * bands]
                            ) < 1e-12
                        );

                        let bits = candidate
                            .orbitals
                            .as_slice()
                            .iter()
                            .chain(&candidate.ritz_values_hartree)
                            .map(|v| v.to_bits())
                            .collect::<Vec<_>>();
                        if gram == "off" {
                            assert_eq!(
                                base.orbitals
                                    .as_slice()
                                    .iter()
                                    .chain(&base.ritz_values_hartree)
                                    .map(|v| v.to_bits())
                                    .collect::<Vec<_>>(),
                                bits
                            );
                        }
                        if let Some(prev) = &previous {
                            assert_eq!(prev, &bits);
                        }
                        previous = Some(bits);
                        // The same tuned operator remains usable after the
                        // constructor cache scope has ended.
                        let mut hx = dft_numerics::RealMatrix::zeros(
                            candidate.orbitals.rows(),
                            candidate.orbitals.cols(),
                        )
                        .unwrap();
                        candidate
                            .operator
                            .apply_block(comm, &candidate.orbitals, &mut hx)
                            .unwrap();
                        assert!(hx.as_slice().iter().all(|v| v.is_finite()));
                    }
                }
            });
        }
    }
}

/// Independent verification still succeeds after the rounding-class start.
#[test]
fn p4_start_verified_trajectories_and_dependent_rejection() {
    use dft_engine::refstart::{StartInput, start};
    use dft_numerics::switches::RunConfig;
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    let potentials = [dft_engine::SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let atoms = [
        dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32]),
        dft_geometry::Atom::new("Xx", [0.57, 0.68, 0.71]),
    ];
    let options = dft_engine::ScfOptions {
        cutoff_hartree: 2.,
        bands: Some(5),
        max_iterations: 100,
        gamma_real: dft_engine::GammaReal::On,
        ..Default::default()
    };
    for skew in [0., 0.2] {
        let cell =
            dft_geometry::Cell::from_user_rows_bohr([[7., 0., 0.], [skew, 7., 0.], [0., skew, 7.]])
                .unwrap();
        for ranks in [1, 3, 8] {
            dft_mpi::simulate(ranks, |comm| {
                let mut base = None;
                for _arm in ["off", "on"] {
                    let config = RunConfig::from_pairs(&[]).unwrap();
                    let _scope = config.enter();
                    let run = dft_engine::refscf::driver::run_with(
                        StartInput {
                            cell: &cell,
                            atoms: &atoms,
                            potentials: &potentials,
                            options: &options,
                        },
                        Profile::default(),
                        comm,
                    )
                    .unwrap();
                    assert!(run.outcome.converged);
                    assert!(run.outcome.result.verification.as_ref().unwrap().passed());
                    let e = run.outcome.result.energy.total_hartree();
                    if let Some(old) = base {
                        let delta: f64 = e - old;
                        assert!(delta.abs() < 1e-12, "{delta:e} Ha");
                        if comm.rank() == 0 {
                            println!(
                                "P4 verified skew={skew} ranks={ranks} E={e:.17e} Ha delta={delta:e} Ha maps={}",
                                run.records.len()
                            );
                        }
                    } else {
                        base = Some(e);
                    }
                    // Same-position duplicate atomic functions must fail on
                    // every rank, including the candidate's guarded fallback.
                    let dependent = [atoms[0].clone(), atoms[0].clone()];
                    assert!(
                        start(
                            StartInput {
                                cell: &cell,
                                atoms: &dependent,
                                potentials: &potentials,
                                options: &options
                            },
                            comm
                        )
                        .is_err()
                    );
                }
            });
        }
    }
}

/// Memory-only replay of the existing density/store traversal. This separates
/// slab allocation/first touch from repeated stores; no FFT or SCF speed claim.
/// Run release, one arm/process with /usr/bin/time -l for process page faults:
/// DENSKEEP_PROBE_ARM=plain|fused|split DENSKEEP_PROBE_MIB=185|500
/// cargo test --release -p dft-engine --test refscf denskeep_store_probe -- --ignored --nocapture
#[test]
#[ignore = "large memory diagnostic; use the documented release command"]
fn denskeep_store_probe() {
    use std::{hint::black_box, time::Instant};
    let arm = std::env::var("DENSKEEP_PROBE_ARM").unwrap_or_else(|_| "fused".into());
    assert!(["plain", "fused", "split"].contains(&arm.as_str()));
    let mib: usize = std::env::var("DENSKEEP_PROBE_MIB")
        .unwrap_or_else(|_| "185".into())
        .parse()
        .unwrap();
    let maps: usize = std::env::var("DENSKEEP_PROBE_MAPS")
        .unwrap_or_else(|_| "16".into())
        .parse()
        .unwrap();
    assert!(mib > 0 && mib <= 512 && maps > 0);
    let pairs = 128;
    let line = 4096;
    let planes = (mib * 1024 * 1024 / (16 * pairs)).div_ceil(line);
    let points = planes * line;
    let elements = pairs * points;
    let slab: Vec<_> = (0..line)
        .map(|j| C::new((j as f64 * 0.3).sin(), (j as f64 * 0.7).cos()))
        .collect();
    let mut density = vec![0.; points];
    let mut keep = Vec::new();
    let start = Instant::now();
    if arm != "plain" {
        keep.resize(elements, C::default());
    }
    // Force the initialized allocation to be observable separately from writes.
    black_box(&mut keep);
    let allocation = start.elapsed().as_secs_f64();
    let mut samples = Vec::new();
    for _ in 0..maps {
        density.fill(0.);
        let start = Instant::now();
        for t in 0..pairs {
            for i in 0..planes {
                let base = i * line;
                let source = black_box(&slab);
                if arm == "fused" {
                    for (j, (dst, v)) in keep[t * points + base..t * points + base + line]
                        .iter_mut()
                        .zip(source)
                        .enumerate()
                    {
                        *dst = *v;
                        density[base + j] += 2. * v.re * v.re / 343.;
                        density[base + j] += 2. * v.im * v.im / 343.;
                    }
                } else {
                    if arm == "split" {
                        keep[t * points + base..t * points + base + line].copy_from_slice(source);
                    }
                    for (j, v) in source.iter().enumerate() {
                        density[base + j] += 2. * v.re * v.re / 343.;
                    }
                    for (j, v) in source.iter().enumerate() {
                        density[base + j] += 2. * v.im * v.im / 343.;
                    }
                }
            }
        }
        black_box((&mut keep, &mut density));
        samples.push(start.elapsed().as_secs_f64());
    }
    let mut expected = vec![0.; line];
    for _ in 0..pairs {
        for (d, v) in expected.iter_mut().zip(&slab) {
            *d += 2. * v.re * v.re / 343.;
            *d += 2. * v.im * v.im / 343.;
        }
    }
    assert!(density.chunks(line).all(|d| {
        d.iter()
            .zip(&expected)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    }));
    let first = samples[0];
    samples.sort_by(f64::total_cmp);
    let median = samples[samples.len() / 2];
    eprintln!(
        "DENSKEEP_STORE arm={arm} maps={maps} pairs={pairs} points={points} line={line} payload_bytes={} store_bytes_per_map={} allocation_s={allocation:.9} first_s={first:.9} median_s={median:.9} logical_store_GB_s={:.6} density_bits=identical",
        keep.len() * 16,
        keep.len() * 16,
        keep.len() as f64 * 16. / median / 1e9
    );
}

// P4-DENSKEEP: recomputation changes storage/FFT work, not precision. Check the
// actual cache path (including zero-plane ranks and an odd unoccupied tail).
#[test]
fn fast_order_denskeep_recompute_products_and_counts() {
    use dft_engine::{GammaReal, ScfOptions, SpeciesPotential, refscf::kept_hx::KeptHx, refstart};
    use dft_numerics::{
        RealMatrix,
        stage_clock::local,
        switches::{Id, RunConfig},
    };
    assert!(
        !RunConfig::from_pairs(&[])
            .unwrap()
            .on(Id::REFSCF_KEPT_RECOMPUTE)
    );
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    let cell =
        dft_geometry::Cell::from_user_rows_bohr([[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]])
            .unwrap();
    let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    for ranks in [1, 2, 9] {
        for bands in [3, 4] {
            let options = ScfOptions {
                cutoff_hartree: 1.,
                bands: Some(bands),
                gamma_real: GammaReal::On,
                ..Default::default()
            };
            dft_mpi::simulate(ranks, |comm| {
                let mut state = refstart::start(
                    refstart::StartInput {
                        cell: &cell,
                        atoms: &atoms,
                        potentials: &potentials,
                        options: &options,
                    },
                    comm,
                )
                .unwrap();
                local::enable(true);
                let x = &state.orbitals;
                let op = &mut state.operator;
                let old: Vec<_> = op
                    .plane_grid_offsets()
                    .iter()
                    .map(|&g| 0.03 * (g as f64).sin())
                    .collect();
                let new: Vec<_> = op
                    .plane_grid_offsets()
                    .iter()
                    .map(|&g| 0.04 * (0.7 * g as f64).cos())
                    .collect();
                op.set_owned_potential(&old, 1).unwrap();
                let mut hx = RealMatrix::zeros(x.rows(), bands).unwrap();
                op.apply_block(comm, x, &mut hx).unwrap();
                let mut kept = Vec::new();
                let mut densities = Vec::new();
                for arm in ["off", "on"] {
                    let _scope = RunConfig::from_pairs(&[("refscf_kept_recompute", arm)])
                        .unwrap()
                        .enter();
                    let mut cache = KeptHx::new(true);
                    let mut density = vec![0.; old.len()];
                    let before = local::snapshot();
                    let mut weights = vec![0.; bands];
                    weights[0] = 2.;
                    cache
                        .density(
                            comm,
                            op,
                            x,
                            &hx,
                            &old,
                            &weights,
                            cell.volume_bohr3(),
                            &mut density,
                        )
                        .unwrap();
                    let counts = local::snapshot().since(&before);
                    assert_eq!(counts.calls[local::OP_COUNT + 3], bands.div_ceil(2) as u64);
                    assert_eq!(counts.calls[local::OP_COUNT + 2], 0);
                    kept.push(cache);
                    densities.push(density);
                }
                assert_eq!(densities[0], densities[1]);
                op.set_owned_potential(&new, 2).unwrap();
                let mut products = Vec::new();
                for (i, cache) in kept.iter_mut().enumerate() {
                    let before = local::snapshot();
                    products.push(cache.prepare(comm, op, x, &new).unwrap().unwrap());
                    let counts = local::snapshot().since(&before);
                    assert_eq!(counts.calls[local::OP_COUNT], 0, "no full H");
                    assert_eq!(counts.calls[local::OP_COUNT + 1], bands as u64);
                    assert_eq!(counts.calls[local::OP_COUNT + 2], bands.div_ceil(2) as u64);
                    assert_eq!(
                        counts.calls[local::OP_COUNT + 3],
                        if i == 0 { bands % 2 } else { bands.div_ceil(2) } as u64
                    );
                    let report = cache.report(comm).unwrap();
                    assert_eq!(report.full_h_starting_columns_replaced, bands);
                    assert_eq!(report.cache_hits, 1);
                    if i == 1 {
                        assert!(report.retained_slab_bytes_per_rank.iter().all(|&b| b == 0));
                    }
                }
                assert!(
                    products[0]
                        .as_slice()
                        .iter()
                        .zip(products[1].as_slice())
                        .all(|(a, b)| a.to_bits() == b.to_bits()),
                    "recomputed HX bits"
                );
                let mut direct = RealMatrix::zeros(x.rows(), bands).unwrap();
                op.apply_block(comm, x, &mut direct).unwrap();
                let local_error = [
                    products[1]
                        .as_slice()
                        .iter()
                        .zip(direct.as_slice())
                        .map(|(a, b)| (a - b).powi(2))
                        .sum(),
                    direct.as_slice().iter().map(|v| v * v).sum(),
                ];
                let mut total = [0.; 2];
                comm.tree_sum_f64(&local_error, &mut total).unwrap();
                let relative = (total[0] / total[1]).sqrt();
                assert!(relative < 1e-13);
                if comm.rank() == 0 {
                    eprintln!(
                        "DENSKEEP_PRODUCT ranks={ranks} bands={bands} relative_HX={relative:.17e} recompute_delta=0 inverse_per_hit={}->{} replaced={bands} slab_bytes=0",
                        bands % 2,
                        bands.div_ceil(2)
                    );
                }
                // A miss cannot turn an invalid block into a hit in either arm.
                if comm.rank() == 0 {
                    kept[1].invalidate();
                }
                assert!(kept[1].prepare(comm, op, x, &new).unwrap().is_none());
            });
        }
    }
}

fn denskeep_delivery(label: &str, rows: [[f64; 3]; 3], nonlocal: bool) {
    use dft_engine::{GammaReal, ScfOptions, SpeciesPotential, refscf::driver, refstart};
    use dft_numerics::switches::RunConfig;
    let pp = common::potential_with_atomic_start(2., 1.8, nonlocal, 161, "SLA PZ");
    let cell = dft_geometry::Cell::from_user_rows_bohr(rows).unwrap();
    let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: 1.,
        bands: Some(3),
        max_iterations: 120,
        gamma_real: GammaReal::On,
        ..Default::default()
    };
    dft_mpi::simulate(2, |comm| {
        let run = |arm| {
            let _scope = RunConfig::from_pairs(&[("refscf_kept_recompute", arm)])
                .unwrap()
                .enter();
            driver::run_with(
                refstart::StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &options,
                },
                Profile {
                    threshold_ry: 1e-12,
                    ..Default::default()
                },
                comm,
            )
            .unwrap()
        };
        let baseline = run("off");
        let candidate = run("on");
        for r in [&baseline, &candidate] {
            assert!(r.outcome.converged, "{label}: {}", r.exit_reason);
            assert!(r.outcome.result.verification.as_ref().unwrap().passed());
            assert!(r.kept_hx.cache_hits > 0);
        }
        let a = &baseline.outcome.result;
        let b = &candidate.outcome.result;
        let delta = b.total_hartree() - a.total_hartree();
        assert!(delta.abs() <= 1e-10, "{label}: {delta:e} Ha");
        assert_eq!(a.total_hartree().to_bits(), b.total_hartree().to_bits());
        assert!(
            a.density
                .iter()
                .zip(&b.density)
                .all(|(a, b)| a.to_bits() == b.to_bits())
        );
        assert_eq!(baseline.records.len(), candidate.records.len());
        for (a, b) in baseline.records.iter().zip(&candidate.records) {
            assert_eq!(
                (a.davidson_iterations, a.h_columns, a.warm_start),
                (b.davidson_iterations, b.h_columns, b.warm_start)
            );
            assert_eq!(a.calls.len(), b.calls.len());
            for (a, b) in a.calls.iter().zip(&b.calls) {
                assert_eq!(a.corrections, b.corrections);
            }
        }
        assert_eq!(
            baseline.verification_calls.len(),
            candidate.verification_calls.len()
        );
        for (a, b) in baseline
            .verification_calls
            .iter()
            .zip(&candidate.verification_calls)
        {
            assert_eq!((a.iterations, a.h_columns), (b.iterations, b.h_columns));
        }
        assert_eq!(baseline.kept_hx.cache_hits, candidate.kept_hx.cache_hits);
        assert_eq!(
            baseline.kept_hx.full_h_starting_columns_replaced,
            candidate.kept_hx.full_h_starting_columns_replaced
        );
        assert!(
            candidate
                .kept_hx
                .retained_slab_bytes_per_rank
                .iter()
                .all(|&n| n == 0)
        );
        if comm.rank() == 0 {
            eprintln!(
                "DENSKEEP_DELIVERY {label} E={:.17e} Ha dE={delta:.1e} Ha maps={} Davidson={} full_H={} hits={} replaced={} slabs={:?}->0 verification=passed density_bits=identical",
                b.total_hartree(),
                candidate.records.len(),
                candidate
                    .records
                    .iter()
                    .map(|r| r.davidson_iterations)
                    .sum::<usize>(),
                candidate.records.iter().map(|r| r.h_columns).sum::<usize>(),
                candidate.kept_hx.cache_hits,
                candidate.kept_hx.full_h_starting_columns_replaced,
                baseline.kept_hx.retained_slab_bytes_per_rank
            );
        }
    });
}
#[test]
fn fast_order_denskeep_delivery_cubic() {
    denskeep_delivery("cubic", [[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]], false);
}
#[test]
fn fast_order_denskeep_delivery_skew_nonlocal() {
    denskeep_delivery(
        "skew-nonlocal",
        [[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]],
        true,
    );
}
#[test]
fn fast_order_denskeep_delivery_elongated() {
    denskeep_delivery(
        "elongated",
        [[5., 0., 0.], [0., 5., 0.], [0., 0., 10.]],
        false,
    );
}
#[test]
fn fast_order_denskeep_delivery_slab() {
    // The thin direction still resolves the p_z atomic starting orbital at
    // this cutoff; z=4 Bohr would make the pre-existing start rank deficient.
    denskeep_delivery("slab", [[10., 0., 0.], [0., 10., 0.], [0., 0., 5.]], false);
}

/// Common-input density replay on a real synthetic operator (FFT included).
/// cargo test --release -p dft-engine --test refscf denskeep_density_probe -- --ignored --nocapture
#[test]
#[ignore = "density timing diagnostic; run in release separately from builds"]
fn denskeep_density_probe() {
    use dft_engine::{GammaReal, ScfOptions, SpeciesPotential, refscf::kept_hx::KeptHx, refstart};
    use dft_numerics::{RealMatrix, stage_clock::local, switches::RunConfig};
    use std::time::Instant;
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    let cell =
        dft_geometry::Cell::from_user_rows_bohr([[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]])
            .unwrap();
    let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: 1.,
        bands: Some(16),
        grid_override: Some([64, 64, 64]),
        gamma_real: GammaReal::On,
        ..Default::default()
    };
    dft_mpi::simulate(2, |comm| {
        let mut state = refstart::start(
            refstart::StartInput {
                cell: &cell,
                atoms: &atoms,
                potentials: &potentials,
                options: &options,
            },
            comm,
        )
        .unwrap();
        let x = &state.orbitals;
        let op = &mut state.operator;
        let potential = vec![0.03; op.plane_points()];
        op.set_owned_potential(&potential, 1).unwrap();
        let mut hx = RealMatrix::zeros(x.rows(), x.cols()).unwrap();
        op.apply_block(comm, x, &mut hx).unwrap();
        let mut out = vec![0.; potential.len()];
        let mut weights = vec![0.; x.cols()];
        weights[0] = 2.;
        let mut oracle = Vec::new();
        local::enable(true);
        for arm in ["off", "on"] {
            let _scope = RunConfig::from_pairs(&[("refscf_kept_recompute", arm)])
                .unwrap()
                .enter();
            let mut kept = KeptHx::new(true);
            let before = local::snapshot();
            let mut times = Vec::new();
            for _ in 0..16 {
                let start = Instant::now();
                kept.density(
                    comm,
                    op,
                    x,
                    &hx,
                    &potential,
                    &weights,
                    cell.volume_bohr3(),
                    &mut out,
                )
                .unwrap();
                times.push(start.elapsed().as_secs_f64());
            }
            let counters = local::snapshot().since(&before);
            assert_eq!(counters.calls[local::OP_COUNT + 3], 128);
            assert_eq!(counters.calls[local::OP_COUNT + 2], 0);
            let bits: Vec<_> = out.iter().map(|v| v.to_bits()).collect();
            if arm == "off" {
                oracle = bits;
            } else {
                assert_eq!(oracle, bits);
            }
            let report = kept.report(comm).unwrap();
            let first = times.remove(0);
            times.sort_by(f64::total_cmp);
            eprintln!(
                "DENSKEEP_DENSITY rank={} recompute={arm} grid=64x64x64 bands=16 maps=16 inverse=128 forward=0 bytes_per_map={} first_s={first:.9} warm_median_s={:.9} density_bits=identical",
                comm.rank(),
                report.retained_slab_bytes_per_rank[comm.rank()],
                times[times.len() / 2]
            );
        }
    });
}

#[allow(dead_code)]
#[path = "common/rank_rule.rs"]
mod nlcc_rank_rule;

#[test]
fn nlcc_owned_start_matches_legacy_and_gaussian_fourier() {
    use dft_engine::{Functional, GammaReal, ScfOptions, ScfProblem, SpeciesPotential};
    use dft_numerics::switches::RunConfig;
    let cell =
        dft_geometry::Cell::from_user_rows_bohr([[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 8.]])
            .unwrap();
    let atoms = [
        dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32]),
        dft_geometry::Atom::new("Yy", [0.57, 0.68, 0.71]),
    ];
    for functional in [Functional::LdaPz, Functional::LdaPw] {
        let name = if functional == Functional::LdaPz {
            "SLA PZ"
        } else {
            "SLA PW"
        };
        let pp = common::potential_with_nlcc(true, name);
        let no_core = common::potential_with_atomic_start(2., 1.8, false, 401, name);
        let species = [
            SpeciesPotential {
                label: "Xx",
                potential: &pp,
            },
            SpeciesPotential {
                label: "Yy",
                potential: &no_core,
            },
        ];
        let options = ScfOptions {
            functional,
            start_density: dft_engine::StartDensity::Atomic,
            cutoff_hartree: 1.,
            bands: Some(3),
            gamma_real: GammaReal::On,
            ..Default::default()
        };
        let old = ScfProblem::new(&cell, &atoms, &species, options.clone()).unwrap();
        let density = old.starting_density();
        let core = old.core_density().unwrap();
        let mut fft = dft_physics::data::FftPlan3d::new(old.grid_dimensions()).unwrap();
        let h = dft_physics::hartree::evaluate(&density, &cell, &mut fft).unwrap();
        let xc = functional
            .evaluate_with_core(
                &density,
                Some(core),
                &mut Array3::zeros(density.raw_dim()),
                &cell,
                &mut fft,
            )
            .unwrap();
        let ionic = old.ionic_potential().clone();
        for _shells in ["off", "on"] {
            for ranks in [1, 3, 8] {
                dft_mpi::simulate(ranks, |comm| {
                    let config = RunConfig::from_pairs(&[]).unwrap();
                    let _scope = config.enter();
                    let mut state = dft_engine::refstart::start(
                        dft_engine::refstart::StartInput {
                            cell: &cell,
                            atoms: &atoms,
                            potentials: &species,
                            options: &options,
                        },
                        comm,
                    )
                    .unwrap();
                    let actual = state
                        .fields
                        .gather(comm, state.core_density.as_ref().unwrap())
                        .unwrap();
                    for (a, b) in actual.iter().zip(core) {
                        assert!((a - b).abs() < 2e-14);
                    }
                    let n = state.fields.gather(comm, &state.density).unwrap();
                    for (a, b) in n.iter().zip(&density) {
                        assert!((a - b).abs() < 2e-14);
                    }
                    let v = state.fields.gather(comm, &state.potential).unwrap();
                    for (((a, ion), h), x) in
                        v.iter().zip(&ionic).zip(&h.potential).zip(&xc.potential)
                    {
                        assert!((a - (ion + h.re + x)).abs() < 2e-13);
                    }
                    // Independent analytic Fourier coefficients, including G=0.
                    let cg = state
                        .fields
                        .forward(comm, state.core_density.as_ref().unwrap())
                        .unwrap();
                    let dims = state.grid.dims();
                    let b = cell.reciprocal_columns_inverse_bohr();
                    for (&z, &g) in cg.iter().zip(state.fields.reciprocal_indices()) {
                        let ijk = [g / (dims[1] * dims[2]), g / dims[2] % dims[1], g % dims[2]];
                        let m: [f64; 3] = std::array::from_fn(|a| {
                            if ijk[a] <= (dims[a] - 1) / 2 {
                                ijk[a] as f64
                            } else {
                                ijk[a] as f64 - dims[a] as f64
                            }
                        });
                        let cart: [f64; 3] =
                            std::array::from_fn(|a| (0..3).map(|k| b[a][k] * m[k]).sum());
                        let g2: f64 = cart.iter().map(|v| v * v).sum();
                        // Only the density sphere is populated (cutoff=4 Ha).
                        if g2 < 7.9 {
                            let angle = -TAU
                                * m.iter()
                                    .zip(atoms[0].fractional_position)
                                    .map(|(m, r)| m * r)
                                    .sum::<f64>();
                            let expected =
                                C::from_polar((-g2 / 2.).exp() / cell.volume_bohr3(), angle);
                            assert!(
                                (z - expected).norm() < 2e-14,
                                "G={m:?} error={:e}",
                                (z - expected).norm()
                            );
                        }
                    }
                });
            }
        }
    }
}

fn nlcc_delivered_fixture(nonlocal: bool, functional: dft_engine::Functional) {
    use dft_engine::{GammaReal, ScfOptions, ScfProblem, SpeciesPotential, StartWavefunctions};
    use dft_numerics::switches::RunConfig;
    let name = if functional == dft_engine::Functional::LdaPz {
        "SLA PZ"
    } else {
        "SLA PW"
    };
    let pp = common::potential_with_nlcc(nonlocal, name);
    let cell = dft_geometry::Cell::from_user_rows_bohr(if nonlocal {
        [[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 8.]]
    } else {
        [[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]]
    })
    .unwrap();
    let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        functional,
        cutoff_hartree: 2.,
        bands: Some(3),
        max_iterations: 120,
        gamma_real: GammaReal::On,
        ..Default::default()
    };
    let old = ScfProblem::new(
        &cell,
        &atoms,
        &species,
        ScfOptions {
            standard: None,
            verification: dft_engine::Verification::Off,
            start_wavefunctions: StartWavefunctions::FreeElectron,
            free_electron_initial_guess: false,
            eigen_tolerance: 1e-11,
            energy_tolerance: 1e-12,
            density_tolerance: 1e-9,
            ..options.clone()
        },
    )
    .unwrap()
    .run()
    .unwrap();
    assert!(old.certified);
    for threshold_ry in [1e-8, 1e-12] {
        let mut baseline = None;
        for ranks in [1, 3, 8] {
            let run_layout = || {
                dft_mpi::simulate(ranks, |comm| {
                    let config = RunConfig::from_pairs(&[("davidson_stages", "on")]).unwrap();
                    let _scope = config.enter();
                    dft_engine::refscf::driver::run_with(
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
                        comm,
                    )
                    .unwrap()
                })
            };
            let runs = run_layout();
            let repeats = run_layout();
            for (run, repeat) in runs.iter().zip(&repeats) {
                assert_eq!(
                    nlcc_rank_rule::fingerprint(run),
                    nlcc_rank_rule::fingerprint(repeat)
                );
                let branches = nlcc_rank_rule::accept(
                    baseline.as_ref().unwrap_or(&runs[0]),
                    run,
                    old.total_hartree(),
                    2.,
                )
                .unwrap();
                for branch in branches {
                    eprintln!("G1_NLCC_BRANCH {branch:?}");
                }
                let r = &run.outcome.result;
                let de = (r.total_hartree() - old.total_hartree()).abs();
                assert!(de < if threshold_ry < 1e-8 { 1e-10 } else { 5e-9 });
                let dn = r
                    .density
                    .iter()
                    .zip(&old.density)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0., f64::max);
                // Empty-band Ritz threshold is 1e-5 Ry; retain a separate cap.
                let ev = (r.eigenvalues[0] - old.eigenvalues[0]).abs();
                let all_ev = r
                    .eigenvalues
                    .iter()
                    .zip(&old.eigenvalues)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0., f64::max);
                assert!(all_ev < 1e-5, "all eigenvalues {all_ev:e} Ha");
                assert!(
                    dn < if threshold_ry < 1e-8 { 1e-7 } else { 1e-6 },
                    "density {dn:e} e/Bohr^3"
                );
                assert!(ev < 1e-5, "occupied eigenvalue {ev:e} Ha");
                let ledger = r.energy.nlcc.as_ref().unwrap();
                assert!((ledger.core_electrons - 1.).abs() < 1e-12);
                assert!(
                    (ledger.core_potential_integral_ha
                        - old.energy.nlcc.as_ref().unwrap().core_potential_integral_ha)
                        .abs()
                        < 1e-5
                );
                assert!((r.energy.band_form - r.total_hartree()).abs() < 1e-10);
                if std::ptr::eq(run, &runs[0]) {
                    eprintln!(
                        "G1_NLCC {functional:?} nonlocal={nonlocal} ranks={ranks} threshold={threshold_ry:e} Ry E={:.17e} Ha dE={de:e} Ha dn={dn:e} e/Bohr^3 deig={ev:e} Ha all_deig={all_ev:e} Ha maps={}",
                        r.total_hartree(),
                        run.records.len()
                    );
                }
            }
            if baseline.is_none() {
                baseline = runs.into_iter().next();
            }
        }
    }
}

// Slow full convergence/repeatability matrices: run each in release with
// cargo test --release -p dft-engine --test refscf nlcc_delivered -- --ignored --nocapture
#[test]
#[ignore = "full convergence and rank-count matrix; run in release"]
fn nlcc_delivered_pz_local() {
    nlcc_delivered_fixture(false, dft_engine::Functional::LdaPz);
}
#[test]
#[ignore = "full convergence and rank-count matrix; run in release"]
fn nlcc_delivered_pw_nonlocal() {
    nlcc_delivered_fixture(true, dft_engine::Functional::LdaPw);
}

/// Before/after build fingerprint for the external no-core Si fixture. The
/// pseudopotential is never written into this repository.
#[test]
#[ignore = "external UPF and full SCF; run in release with DFT_SI_UPF"]
fn nlcc_no_core_silicon_fingerprint() {
    use dft_engine::{GammaReal, ScfOptions, SpeciesPotential};
    let pp = dft_upf::parse_upf(
        &std::fs::read_to_string(std::env::var("DFT_SI_UPF").expect("DFT_SI_UPF required"))
            .unwrap(),
    )
    .unwrap();
    assert!(pp.core_density().is_none());
    let cell = common::cell();
    let atoms = common::atoms();
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: 4.,
        bands: Some(5),
        gamma_real: GammaReal::On,
        ..Default::default()
    };
    for ranks in [1, 3, 8] {
        let runs = dft_mpi::simulate(ranks, |comm| {
            let config =
                dft_numerics::switches::RunConfig::from_pairs(&[("davidson_stages", "on")])
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
        });
        let run = &runs[0];
        assert!(run.outcome.converged);
        assert!(run.outcome.result.verification.as_ref().unwrap().passed());
        let fingerprint = nlcc_rank_rule::fingerprint(run);
        let hash = fingerprint.bytes().fold(0xcbf29ce484222325_u64, |h, b| {
            (h ^ u64::from(b)).wrapping_mul(0x100000001b3)
        });
        eprintln!(
            "G1_SI_NO_CORE ranks={ranks} hash={hash:016x} E={:.17e} Ha",
            run.outcome.result.total_hartree()
        );
    }
}

/// Run this exact test alone under mpirun, with dft-mpi/mpi enabled. Without
/// that feature it also checks the harness through the serial World adapter.
#[test]
#[ignore = "dev08 real MPI harness; one exact test per process invocation"]
fn nlcc_real_mpi_legacy_comparison() {
    use dft_engine::{GammaReal, ScfOptions, ScfProblem, SpeciesPotential};
    let universe = dft_mpi::Universe::init().unwrap();
    let world = universe.world();
    let config =
        dft_numerics::switches::RunConfig::from_pairs(&[("davidson_stages", "on")]).unwrap();
    let _scope = config.enter();
    for nonlocal in [false, true] {
        let functional = if nonlocal {
            dft_engine::Functional::LdaPw
        } else {
            dft_engine::Functional::LdaPz
        };
        let pp = common::potential_with_nlcc(nonlocal, if nonlocal { "SLA PW" } else { "SLA PZ" });
        let cell = dft_geometry::Cell::from_user_rows_bohr(if nonlocal {
            [[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 8.]]
        } else {
            [[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]]
        })
        .unwrap();
        let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
        let species = [SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let options = ScfOptions {
            functional,
            cutoff_hartree: 2.,
            bands: Some(3),
            max_iterations: 120,
            gamma_real: GammaReal::On,
            ..Default::default()
        };
        // World is Copy in the serial adapter, Clone in the MPI adapter.
        #[allow(clippy::clone_on_copy)]
        dft_engine::parallel::set_communicator(Some(std::sync::Arc::new(world.clone())));
        let old = ScfProblem::new(
            &cell,
            &atoms,
            &species,
            ScfOptions {
                standard: None,
                verification: dft_engine::Verification::Off,
                start_wavefunctions: dft_engine::StartWavefunctions::FreeElectron,
                free_electron_initial_guess: false,
                eigen_tolerance: 1e-11,
                energy_tolerance: 1e-12,
                density_tolerance: 1e-9,
                ..options.clone()
            },
        )
        .unwrap()
        .run()
        .unwrap();
        dft_engine::parallel::set_communicator(None);
        for threshold_ry in [1e-8, 1e-12] {
            let run = dft_engine::refscf::driver::run_with(
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
                &world,
            )
            .unwrap();
            assert!(nlcc_rank_rule::delivered(&run, old.total_hartree(), 2.));
            let r = &run.outcome.result;
            let de = (r.total_hartree() - old.total_hartree()).abs();
            assert!(de < if threshold_ry < 1e-8 { 1e-10 } else { 5e-9 });
            let dn = r
                .density
                .iter()
                .zip(&old.density)
                .map(|(a, b)| (a - b).abs())
                .fold(0., f64::max);
            let ev = (r.eigenvalues[0] - old.eigenvalues[0]).abs();
            assert!(
                dn < 1e-6 && ev < 1e-5,
                "dE={de:e} Ha dn={dn:e} e/Bohr^3 deig={ev:e} Ha"
            );
            assert!(r.energy.nlcc.is_some());
            if world.rank() == 0 {
                eprintln!(
                    "G1_NLCC_MPI ranks={} nonlocal={nonlocal} threshold={threshold_ry:e} Ry E={:.17e} Ha dE={de:e} Ha dn={dn:e} e/Bohr^3 deig={ev:e} Ha",
                    world.size(),
                    r.total_hartree()
                );
            }
        }
    }
}

// G2: full-grid objects below are independent test oracles only. Production
// gradients and XC retain exclusively owned slabs and reciprocal sticks.
#[test]
fn pbe_distributed_gradients_start_and_variational_potential() {
    // Historical fixture contract; screening has separate charge/mask tests.
    let screening =
        dft_numerics::switches::RunConfig::from_pairs(&[("xc_screen_negative", "off")]).unwrap();
    let _screen_scope = screening.enter();
    use dft_engine::{
        Functional, ScfOptions, SpeciesPotential,
        refstart::{self, StartInput},
    };
    use dft_physics::{data::FftPlan3d, gga};
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "PBE");
    let cell =
        dft_geometry::Cell::from_user_rows_bohr([[7., 0., 0.], [0.5, 6., 0.], [0.3, 0.2, 5.]])
            .unwrap();
    let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    for dims in [[8, 10, 12], [9, 11, 7], [7, 8, 9]] {
        let options = ScfOptions {
            functional: Functional::Pbe,
            cutoff_hartree: 1.,
            bands: Some(3),
            grid_override: Some(dims),
            ..Default::default()
        };
        let n = Array3::from_shape_fn(dims, |(i, j, k)| {
            let q = [i, j, k];
            0.03 + (0..3)
                .map(|a| 0.003 * (TAU * q[a] as f64 / dims[a] as f64).sin())
                .sum::<f64>()
                + if dims[0] % 2 == 0 {
                    0.002 * (PI * i as f64).cos() * (TAU * j as f64 / dims[1] as f64).sin()
                } else {
                    0.
                }
        });
        let mut fft = FftPlan3d::new(dims).unwrap();
        let serial = gga::evaluate(&n, &cell, &mut fft).unwrap();
        let delta = Array3::from_shape_fn(dims, |(i, j, _)| {
            (TAU * i as f64 / dims[0] as f64).sin()
                + 0.3 * (TAU * (i as f64 / dims[0] as f64 + j as f64 / dims[1] as f64)).cos()
        });
        let h = 1e-6;
        let derivative = (gga::evaluate(&(&n + &delta.mapv(|v| h * v)), &cell, &mut fft)
            .unwrap()
            .energy
            .0
            - gga::evaluate(&(&n - &delta.mapv(|v| h * v)), &cell, &mut fft)
                .unwrap()
                .energy
                .0)
            / (2. * h);
        assert!(
            derivative.abs() > 1e-3,
            "nonzero directional derivative oracle"
        );
        for ranks in [1, 3, 8] {
            let errors = dft_mpi::simulate(ranks, |comm| {
                let _screen_scope = screening.enter();
                let mut max_gradient = 0.0_f64;
                let mut max_potential = 0.0_f64;
                let mut state = refstart::start(
                    StartInput {
                        cell: &cell,
                        atoms: &atoms,
                        potentials: &potentials,
                        options: &options,
                    },
                    comm,
                )
                .unwrap();
                // Start potential must also use the GGA divergence, not vrho alone.
                let start_n = state.fields.gather(comm, &state.density).unwrap();
                let mut fft = FftPlan3d::new(dims).unwrap();
                let sx = gga::evaluate(&start_n, &cell, &mut fft).unwrap();
                let sh = dft_physics::hartree::evaluate(&start_n, &cell, &mut fft).unwrap();
                let offsets = state.fields.plane_offsets();
                for (j, &i) in offsets.iter().enumerate() {
                    assert!(
                        (state.potential[j]
                            - state.ionic[j]
                            - sh.potential.as_slice().unwrap()[i].re
                            - sx.potential.as_slice().unwrap()[i])
                            .abs()
                            < 2e-12
                    );
                }
                let local: Vec<_> = offsets.iter().map(|&i| n.as_slice().unwrap()[i]).collect();
                let coeff = state.fields.forward(comm, &local).unwrap();
                let owned_comm = PbeNoGather(comm);
                let grad =
                    gga::distributed::gradient(&mut state.fields, &owned_comm, &cell, &coeff)
                        .unwrap();
                let b = cell.reciprocal_columns_inverse_bohr();
                for axis in 0..3 {
                    for (j, &index) in offsets.iter().enumerate() {
                        let q = [
                            index / (dims[1] * dims[2]),
                            index / dims[2] % dims[1],
                            index % dims[2],
                        ];
                        let mut exact: f64 = (0..3)
                            .map(|a| {
                                0.003 * b[axis][a] * (TAU * q[a] as f64 / dims[a] as f64).cos()
                            })
                            .sum();
                        if dims[0] % 2 == 0 {
                            exact += 0.002
                                * b[axis][1]
                                * (PI * q[0] as f64).cos()
                                * (TAU * q[1] as f64 / dims[1] as f64).cos();
                        }
                        max_gradient = max_gradient.max((grad[axis][j] - exact).abs());
                        assert!(
                            (grad[axis][j] - exact).abs() < 3e-15,
                            "P={ranks} dims={dims:?} axis={axis} gradient {} vs {exact}",
                            grad[axis][j]
                        );
                    }
                }
                let x = gga::distributed::evaluate(
                    &mut state.fields,
                    &owned_comm,
                    &cell,
                    &local,
                    &coeff,
                )
                .unwrap();
                for (&v, &i) in x.potential.iter().zip(&offsets) {
                    max_potential =
                        max_potential.max((v - serial.potential.as_slice().unwrap()[i]).abs());
                    assert!((v - serial.potential.as_slice().unwrap()[i]).abs() < 2e-12);
                }
                let dv = cell.volume_bohr3() / n.len() as f64;
                let dot = dv
                    * x.potential
                        .iter()
                        .zip(&offsets)
                        .map(|(v, &i)| v * delta.as_slice().unwrap()[i])
                        .sum::<f64>();
                let mut sums = [0.; 2];
                comm.tree_sum_f64(&[x.energy, dot], &mut sums).unwrap();
                assert!((sums[0] - serial.energy.0).abs() < 2e-12);
                assert!(
                    (sums[1] - derivative).abs() < 1e-8,
                    "variational derivative Ha"
                );
                // One-rank invalid density must fail collectively before any FFT.
                let mut bad = local.clone();
                if comm.rank() == 0 && !bad.is_empty() {
                    bad[0] = -1e-15;
                }
                assert!(
                    gga::distributed::evaluate(&mut state.fields, &owned_comm, &cell, &bad, &coeff)
                        .is_err()
                );
                [
                    max_gradient,
                    max_potential,
                    (sums[0] - serial.energy.0).abs(),
                    (sums[1] - derivative).abs(),
                ]
            });
            let maxima: [f64; 4] =
                std::array::from_fn(|i| errors.iter().map(|v| v[i]).fold(0., f64::max));
            eprintln!(
                "G2-PBE grid={dims:?} P={ranks} max gradient={} Bohr^-4 potential={} Ha energy={} Ha derivative={} Ha",
                maxima[0], maxima[1], maxima[2], maxima[3]
            );
        }
    }
}

fn pbe_legacy_fixture(label: &str, cell: dft_geometry::Cell, one_atom: bool) {
    use dft_engine::{
        Functional, GammaReal, ScfOptions, ScfProblem, SpeciesPotential, StartWavefunctions,
        refstart::StartInput,
    };
    use dft_numerics::switches::RunConfig;
    let pp = common::potential_with_atomic_start(2., 1.8, false, 161, "PBE");
    let all_atoms = [
        dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32]),
        dft_geometry::Atom::new("Xx", [0.57, 0.68, 0.71]),
    ];
    let atoms = if one_atom {
        &all_atoms[..1]
    } else {
        &all_atoms[..]
    };
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let base = ScfOptions {
        functional: Functional::Pbe,
        cutoff_hartree: 2.,
        bands: Some(if one_atom { 3 } else { 5 }),
        max_iterations: 160,
        gamma_real: GammaReal::On,
        ..Default::default()
    };
    let strict = ScfProblem::new(
        &cell,
        atoms,
        &potentials,
        ScfOptions {
            functional: base.functional,
            cutoff_hartree: base.cutoff_hartree,
            bands: base.bands,
            max_iterations: base.max_iterations,
            gamma_real: base.gamma_real,
            free_electron_initial_guess: false,
            start_wavefunctions: StartWavefunctions::FreeElectron,
            ..ScfOptions::strict_reference()
        },
    )
    .unwrap()
    .run()
    .unwrap();
    assert!(strict.certified);
    let reference = strict.energy.total_hartree();
    let mut energies = Vec::new();
    for ranks in [1, 3, 8] {
        let answers = dft_mpi::simulate(ranks, |comm| {
            let mut values = Vec::new();
            for threshold_ry in [1e-8, 1e-12] {
                let mut full = None;
                for lean in ["off", "on"] {
                    let config = RunConfig::from_pairs(&[]).unwrap();
                    let _scope = config.enter();
                    let run = dft_engine::refscf::driver::run_with(
                        StartInput {
                            cell: &cell,
                            atoms,
                            potentials: &potentials,
                            options: &base,
                        },
                        Profile {
                            threshold_ry,
                            ..Default::default()
                        },
                        comm,
                    )
                    .unwrap();
                    assert!(
                        run.outcome.converged,
                        "{label} P={ranks}: {}",
                        run.exit_reason
                    );
                    let r = &run.outcome.result;
                    assert!(r.verification.as_ref().unwrap().passed());
                    let e = r.energy.total_hartree();
                    let cap = if threshold_ry == 1e-8 { 5e-9 } else { 1e-10 };
                    assert!(
                        (e - reference).abs() < cap,
                        "{label} P={ranks}: E-Estrict={} Ha",
                        e - reference
                    );
                    assert!((e - r.energy.band_form).abs() < 1e-10);
                    assert!(
                        (e - r.energy.harris_foulkes).abs()
                            <= r.energy.functional_difference_bound + 1e-12
                    );
                    if let Some((energy, density)) = &full {
                        assert_eq!(*energy, e.to_bits(), "lean delivery energy");
                        assert_eq!(density, &r.density, "lean delivery density");
                    } else {
                        full = Some((e.to_bits(), r.density.clone()));
                    }
                    if comm.rank() == 0 {
                        eprintln!(
                            "G2-PBE {label} P={ranks} lean={lean} thr={threshold_ry:e} Ry E={e:.17e} Ha delta={:.9e} Ha maps={}",
                            e - reference,
                            run.records.len()
                        );
                    }
                    values.push(e);
                }
            }
            values
        });
        for a in &answers[1..] {
            assert_eq!(a, &answers[0]);
        }
        energies.push(answers[0].clone());
    }
    for e in &energies[1..] {
        for (a, b) in e.iter().zip(&energies[0]) {
            assert!((a - b).abs() < 1e-10);
        }
    }
}

#[test]
#[ignore = "SCF matrix; cargo test --release -p dft-engine --test refscf pbe_legacy -- --ignored --nocapture"]
fn pbe_legacy_silicon_geometry_synthetic() {
    pbe_legacy_fixture("synthetic-diamond", common::cell(), false);
}
#[test]
#[ignore = "SCF matrix; cargo test --release -p dft-engine --test refscf pbe_legacy -- --ignored --nocapture"]
fn pbe_legacy_skew_local() {
    pbe_legacy_fixture(
        "synthetic-skew-local",
        dft_geometry::Cell::from_user_rows_bohr([[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]])
            .unwrap(),
        true,
    );
}

// Field evaluation can transpose owned values and reduce scalars, not gather.
struct PbeNoGather<'a, K>(&'a K);
impl<K: Collective> Collective for PbeNoGather<'_, K> {
    fn rank(&self) -> usize {
        self.0.rank()
    }
    fn size(&self) -> usize {
        self.0.size()
    }
    fn barrier(&self) -> Result<(), dft_mpi::MpiError> {
        self.0.barrier()
    }
    fn tree_sum_f64(&self, a: &[f64], b: &mut [f64]) -> Result<(), dft_mpi::MpiError> {
        assert!(a.len() <= 2, "no field reductions");
        self.0.tree_sum_f64(a, b)
    }
    fn bcast_f64(&self, _: &mut [f64], _: usize) -> Result<(), dft_mpi::MpiError> {
        panic!("no field broadcasts");
    }
    fn allgatherv_f64(
        &self,
        _: &[f64],
        _: &[usize],
        _: &mut [f64],
    ) -> Result<(), dft_mpi::MpiError> {
        panic!("no field gathers");
    }
    fn alltoall_f64(&self, a: &[f64], n: usize, b: &mut [f64]) -> Result<(), dft_mpi::MpiError> {
        self.0.alltoall_f64(a, n, b)
    }
}

#[test]
fn pbe_negative_mixed_density_fails_collectively() {
    let screening =
        dft_numerics::switches::RunConfig::from_pairs(&[("xc_screen_negative", "off")]).unwrap();
    let _scope = screening.enter();
    use dft_engine::{Functional, GammaReal, ScfOptions, SpeciesPotential, refstart::StartInput};
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "PBE");
    let cell =
        dft_geometry::Cell::from_user_rows_bohr([[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]])
            .unwrap();
    let atoms = [
        dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32]),
        dft_geometry::Atom::new("Xx", [0.57, 0.68, 0.71]),
    ];
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        functional: Functional::Pbe,
        cutoff_hartree: 2.,
        bands: Some(5),
        gamma_real: GammaReal::On,
        ..Default::default()
    };
    dft_mpi::simulate(3, |comm| {
        let _scope = screening.enter();
        let result = dft_engine::refscf::driver::run_with(
            StartInput {
                cell: &cell,
                atoms: &atoms,
                potentials: &potentials,
                options: &options,
            },
            Profile::default(),
            comm,
        );
        let Err(error) = result else {
            panic!("negative PBE density must not be silently repaired")
        };
        assert!(
            error.to_string().contains("finite nonnegative density"),
            "{error}"
        );
    });
}

/// Allocation ownership only: every scientific bit and solver decision must
/// survive trimming. Simulated ranks cover uneven/empty owners, not MPI speed.
#[test]
fn memory_trim_verified_bits_and_rank_census() {
    memory_trim_case(
        "skew",
        [[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]],
        true,
        3,
    );
}

/// Broader shape/thread matrix (under one second in the Mac debug gate).
#[test]
fn memory_trim_generality() {
    for (name, rows, nl, ranks) in [
        (
            "cubic",
            [[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]],
            false,
            1,
        ),
        (
            "elongated",
            [[6., 0., 0.], [0., 6., 0.], [0., 0., 11.]],
            true,
            8,
        ),
        (
            "slab",
            [[10., 0., 0.], [0.2, 10., 0.], [0., 0., 5.]],
            false,
            3,
        ),
    ] {
        memory_trim_case(name, rows, nl, ranks);
    }
}
fn memory_trim_case(name: &str, rows: [[f64; 3]; 3], nonlocal: bool, ranks: usize) {
    use dft_engine::{
        GammaReal, ScfOptions, SpeciesPotential,
        refscf::{driver, memory},
        refstart,
    };
    use dft_numerics::switches::RunConfig;
    let pp = common::potential_with_atomic_start(2., 1.8, nonlocal, 161, "SLA PZ");
    let cell = dft_geometry::Cell::from_user_rows_bohr(rows).unwrap();
    let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    dft_mpi::simulate(ranks, |comm| {
        let run = |arm, threads| {
            let _scope = RunConfig::from_pairs(&[("refscf_memory_trim", arm)])
                .unwrap()
                .enter();
            let options = ScfOptions {
                cutoff_hartree: 1.,
                bands: Some(3),
                max_iterations: 120,
                gamma_real: GammaReal::On,
                threads,
                ..Default::default()
            };
            driver::run_with(
                refstart::StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &options,
                },
                Profile {
                    threshold_ry: 1e-12,
                    ..Default::default()
                },
                comm,
            )
            .unwrap()
        };
        let baseline = run("off", 1);
        for threads in [1, 2] {
            let candidate = run("on", threads);
            for r in [&baseline, &candidate] {
                assert!(r.outcome.converged, "{name}: {}", r.exit_reason);
                assert!(r.outcome.result.verification.as_ref().unwrap().passed());
                assert_eq!(r.memory_breakdown.components_per_rank.len(), ranks);
                assert_eq!(r.memory_breakdown.davidson_per_rank.len(), ranks);
                assert!(
                    r.memory_breakdown
                        .components_per_rank
                        .iter()
                        .all(|r| r.len() == memory::NAMES.len())
                );
            }
            let a = &baseline.outcome.result;
            let b = &candidate.outcome.result;
            assert_eq!(a.total_hartree().to_bits(), b.total_hartree().to_bits());
            let bits = |s: &[f64]| s.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
            assert_eq!(
                bits(a.density.as_slice().unwrap()),
                bits(b.density.as_slice().unwrap())
            );
            assert_eq!(a.orbitals, b.orbitals);
            assert_eq!(bits(&a.eigenvalues), bits(&b.eigenvalues));
            assert_eq!(baseline.records.len(), candidate.records.len());
            for (a, b) in baseline.records.iter().zip(&candidate.records) {
                assert_eq!(
                    (a.dr2_ry.to_bits(), a.davidson_iterations, a.h_columns),
                    (b.dr2_ry.to_bits(), b.davidson_iterations, b.h_columns)
                );
                for (a, b) in a.calls.iter().zip(&b.calls) {
                    assert_eq!(a.corrections, b.corrections);
                }
            }
            assert_eq!(
                baseline.verification_calls.len(),
                candidate.verification_calls.len()
            );
            for (a, b) in baseline
                .verification_calls
                .iter()
                .zip(&candidate.verification_calls)
            {
                assert_eq!((a.iterations, a.h_columns), (b.iterations, b.h_columns));
                assert_eq!(a.corrections, b.corrections);
            }
            let idx = |name| memory::NAMES.iter().position(|&n| n == name).unwrap();
            for (a, b) in baseline
                .memory_breakdown
                .components_per_rank
                .iter()
                .zip(&candidate.memory_breakdown.components_per_rank)
            {
                assert_eq!(a[idx("start_vectors_released")], 0);
                assert!(b[idx("start_vectors_released")] > 0);
                assert_eq!(a[idx("kept_slabs")], b[idx("kept_slabs")]);
            }
            for (a, b) in baseline
                .memory_breakdown
                .davidson_per_rank
                .iter()
                .zip(&candidate.memory_breakdown.davidson_per_rank)
            {
                assert_eq!(a[0..6], b[0..6]); // Same live solver workspace.
                assert_eq!(a[8], 2 * b[8]); // Returned X/HX use n rather than 2n.
            }
            if comm.rank() == 0 {
                eprintln!(
                    "P4_MEM {name} ranks={ranks} threads={threads} E={:.17e} Ha dE=0 Ha density/orbitals/decisions=bit-identical returned_bytes={:?}->{:?} released_start_bytes={:?}",
                    b.total_hartree(),
                    baseline
                        .memory_breakdown
                        .davidson_per_rank
                        .iter()
                        .map(|r| r[8])
                        .collect::<Vec<_>>(),
                    candidate
                        .memory_breakdown
                        .davidson_per_rank
                        .iter()
                        .map(|r| r[8])
                        .collect::<Vec<_>>(),
                    candidate
                        .memory_breakdown
                        .components_per_rank
                        .iter()
                        .map(|r| r[idx("start_vectors_released")])
                        .collect::<Vec<_>>()
                );
            }
        }
        // Rollback changes the cache shape/key; repeating off must restore bits.
        let rollback = run("off", 1);
        assert_eq!(
            baseline.outcome.result.orbitals,
            rollback.outcome.result.orbitals
        );
    });
}

// G2B: core Fourier coefficients and gradients are compared with a full serial
// FFT of the independently built legacy total density, including skew/Nyquist.
#[test]
fn pbe_nlcc_total_gradient_and_start_match_serial() {
    use dft_engine::{Functional, ScfOptions, ScfProblem, SpeciesPotential, refstart::StartInput};
    use dft_numerics::Fft3d;
    use dft_physics::{data::FftPlan3d, gga};
    let pp = common::potential_with_nlcc(true, "PBE");
    let cell =
        dft_geometry::Cell::from_user_rows_bohr([[5., 0., 0.], [0.2, 5., 0.], [0.1, 0.2, 5.]])
            .unwrap();
    let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    for dims in [[8, 10, 12], [9, 11, 7], [7, 8, 9]] {
        let options = ScfOptions {
            functional: Functional::Pbe,
            cutoff_hartree: 3.,
            start_density: dft_engine::StartDensity::Atomic,
            bands: Some(3),
            grid_override: Some(dims),
            ..Default::default()
        };
        let old = ScfProblem::new(&cell, &atoms, &species, options.clone()).unwrap();
        let n = old.starting_density();
        let core = old.core_density().unwrap();
        let total = &n + core;
        let mut fft = FftPlan3d::new(dims).unwrap();
        let x = gga::evaluate(&total, &cell, &mut fft).unwrap();
        let h = dft_physics::hartree::evaluate(&n, &cell, &mut fft).unwrap();
        let mut coeff = Array3::zeros(dims);
        fft.to_reciprocal(&total.mapv(c), &mut coeff).unwrap();
        let b = cell.reciprocal_columns_inverse_bohr();
        let gradients: Vec<_> = (0..3)
            .map(|axis| {
                let mut derivative = coeff.clone();
                for ((i, j, k), z) in derivative.indexed_iter_mut() {
                    let indices = [i, j, k];
                    let wavevector: f64 = (0..3)
                        .map(|a| {
                            let m = if dims[a] % 2 == 0 && indices[a] == dims[a] / 2 {
                                0.
                            } else if indices[a] <= (dims[a] - 1) / 2 {
                                indices[a] as f64
                            } else {
                                indices[a] as f64 - dims[a] as f64
                            };
                            b[axis][a] * m
                        })
                        .sum();
                    *z *= C::new(0., wavevector);
                }
                let mut out = Array3::zeros(dims);
                fft.to_real(&derivative, &mut out).unwrap();
                out.mapv(|z| z.re)
            })
            .collect();
        for ranks in [1, 3, 8] {
            let errors = dft_mpi::simulate(ranks, |comm| {
                let mut state = dft_engine::refstart::start(
                    StartInput {
                        cell: &cell,
                        atoms: &atoms,
                        potentials: &species,
                        options: &options,
                    },
                    comm,
                )
                .unwrap();
                let offsets = state.fields.plane_offsets();
                let cached = state.pbe_core.as_ref().unwrap().density_g();
                let total_g: Vec<_> = state
                    .density_g
                    .iter()
                    .zip(cached)
                    .map(|(n, c)| n + c)
                    .collect();
                let owned = PbeNoGather(comm);
                let actual =
                    gga::distributed::gradient(&mut state.fields, &owned, &cell, &total_g).unwrap();
                let mut max_gradient = 0_f64;
                let mut max_potential = 0_f64;
                for (local, &global) in offsets.iter().enumerate() {
                    assert!(
                        (state.core_density.as_ref().unwrap()[local]
                            - core.as_slice().unwrap()[global])
                            .abs()
                            < 2e-14
                    );
                    for axis in 0..3 {
                        max_gradient = max_gradient.max(
                            (actual[axis][local] - gradients[axis].as_slice().unwrap()[global])
                                .abs(),
                        );
                    }
                    let expected = state.ionic[local]
                        + h.potential.as_slice().unwrap()[global].re
                        + x.potential.as_slice().unwrap()[global];
                    max_potential = max_potential.max((state.potential[local] - expected).abs());
                }
                assert!(max_gradient < 2e-13, "gradient {max_gradient:e} Bohr^-4");
                assert!(
                    max_potential < 2e-12,
                    "start potential {max_potential:e} Ha"
                );
                [max_gradient, max_potential]
            });
            let maxima: [f64; 2] =
                std::array::from_fn(|i| errors.iter().map(|v| v[i]).fold(0., f64::max));
            eprintln!(
                "G2B_NLCC grid={dims:?} ranks={ranks} gradient={} Bohr^-4 start_potential={} Ha",
                maxima[0], maxima[1]
            );
        }
    }
}

fn pbe_nlcc_delivered_fixture(nonlocal: bool) {
    use dft_engine::{
        Functional, GammaReal, ScfOptions, ScfProblem, SpeciesPotential, StartWavefunctions,
        refstart::StartInput,
    };
    use dft_numerics::switches::RunConfig;
    let pp = common::potential_with_nlcc(nonlocal, "PBE");
    let cell = dft_geometry::Cell::from_user_rows_bohr(if nonlocal {
        [[5., 0., 0.], [0.2, 5., 0.], [0., 0.1, 6.]]
    } else {
        [[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]]
    })
    .unwrap();
    let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        functional: Functional::Pbe,
        cutoff_hartree: if nonlocal { 3. } else { 2. },
        start_density: dft_engine::StartDensity::Atomic,
        bands: Some(3),
        max_iterations: 160,
        gamma_real: GammaReal::On,
        ..Default::default()
    };
    let old = ScfProblem::new(
        &cell,
        &atoms,
        &species,
        ScfOptions {
            standard: None,
            verification: dft_engine::Verification::Off,
            start_wavefunctions: StartWavefunctions::FreeElectron,
            free_electron_initial_guess: false,
            eigen_tolerance: 1e-11,
            energy_tolerance: 1e-12,
            density_tolerance: 1e-9,
            ..options.clone()
        },
    )
    .unwrap();
    let core = old.core_density().unwrap().clone();
    let dims = old.grid_dimensions();
    let old = old.run().unwrap();
    assert!(old.certified);
    let mut rank_energies = Vec::new();
    for ranks in [1, 3, 8] {
        let answers = dft_mpi::simulate(ranks, |comm| {
            let mut values = Vec::new();
            for threshold_ry in [1e-8, 1e-12] {
                let mut full = None;
                for lean in ["off", "on"] {
                    let config = RunConfig::from_pairs(&[]).unwrap();
                    let _scope = config.enter();
                    let run = dft_engine::refscf::driver::run_with(
                        StartInput {
                            cell: &cell,
                            atoms: &atoms,
                            potentials: &species,
                            options: &options,
                        },
                        Profile {
                            threshold_ry,
                            ..Default::default()
                        },
                        comm,
                    )
                    .unwrap();
                    assert!(run.outcome.converged, "{}", run.exit_reason);
                    let r = &run.outcome.result;
                    assert!(r.verification.as_ref().unwrap().passed());
                    let energy = r.total_hartree();
                    let de = (energy - old.total_hartree()).abs();
                    assert!(
                        de < if threshold_ry == 1e-8 { 5e-9 } else { 1e-10 },
                        "dE={de:e} Ha"
                    );
                    assert!((energy - r.energy.band_form).abs() < 1e-10);
                    assert!(
                        (energy - r.energy.harris_foulkes).abs()
                            <= r.energy.functional_difference_bound + 1e-12
                    );
                    let ledger = r.energy.nlcc.as_ref().unwrap();
                    assert!((ledger.core_electrons - 1.).abs() < 1e-12);
                    assert!(
                        (ledger.core_potential_integral_ha
                            - old.energy.nlcc.as_ref().unwrap().core_potential_integral_ha)
                            .abs()
                            < 1e-5
                    );
                    // Independently recompute XC on the delivered valence + fixed
                    // core, including the valence-only double-counting integral.
                    let mut fft = dft_physics::data::FftPlan3d::new(dims).unwrap();
                    let xc = Functional::Pbe
                        .evaluate_with_core(
                            &r.density,
                            Some(&core),
                            &mut Array3::zeros(r.density.raw_dim()),
                            &cell,
                            &mut fft,
                        )
                        .unwrap();
                    assert!((r.energy.xc - xc.energy.0).abs() < 2e-12);
                    assert!(
                        (ledger.valence_potential_integral_ha - xc.potential_integral.0).abs()
                            < 2e-12
                    );
                    let dv = cell.volume_bohr3() / r.density.len() as f64;
                    let core_integral = dv
                        * core
                            .iter()
                            .zip(&xc.potential)
                            .map(|(c, v)| c * v)
                            .sum::<f64>();
                    assert!((ledger.core_potential_integral_ha - core_integral).abs() < 2e-12);
                    assert!((dv * r.density.sum() - 2.).abs() < 1e-12);
                    if let Some((e, n)) = &full {
                        assert_eq!(*e, energy.to_bits(), "lean energy");
                        assert_eq!(n, &r.density, "lean density");
                    } else {
                        full = Some((energy.to_bits(), r.density.clone()));
                    }
                    if comm.rank() == 0 {
                        eprintln!(
                            "G2B_NLCC nonlocal={nonlocal} ranks={ranks} lean={lean} threshold={threshold_ry:e} Ry E={energy:.17e} Ha dE={de:e} Ha maps={}",
                            run.records.len()
                        );
                    }
                    values.push(energy);
                }
            }
            values
        });
        for answer in &answers[1..] {
            assert_eq!(answer, &answers[0]);
        }
        rank_energies.push(answers[0].clone());
    }
    for energies in &rank_energies[1..] {
        for (a, b) in energies.iter().zip(&rank_energies[0]) {
            assert!((a - b).abs() < 1e-10);
        }
    }
}

#[test]
#[ignore = "SCF matrix; cargo test --release -p dft-engine --test refscf pbe_nlcc_delivered -- --ignored --nocapture"]
fn pbe_nlcc_delivered_local() {
    pbe_nlcc_delivered_fixture(false);
}

#[test]
#[ignore = "SCF matrix; cargo test --release -p dft-engine --test refscf pbe_nlcc_delivered -- --ignored --nocapture"]
fn pbe_nlcc_delivered_nonlocal() {
    pbe_nlcc_delivered_fixture(true);
}

#[test]
fn pbe_nlcc_negative_total_start_fails_collectively() {
    let screening =
        dft_numerics::switches::RunConfig::from_pairs(&[("xc_screen_negative", "off")]).unwrap();
    let _scope = screening.enter();
    use dft_engine::{
        Functional, ScfOptions, ScfProblem, SpeciesPotential, StartDensity, refstart::StartInput,
    };
    let pp = common::potential_with_nlcc(true, "PBE");
    let cell =
        dft_geometry::Cell::from_user_rows_bohr([[7., 0., 0.], [0.5, 6., 0.], [0.3, 0.2, 5.]])
            .unwrap();
    let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        functional: Functional::Pbe,
        cutoff_hartree: 1.,
        bands: Some(3),
        start_density: StartDensity::Atomic,
        grid_override: Some([8, 10, 12]),
        ..Default::default()
    };
    // This deliberately under-resolved core rings below zero where the atomic
    // valence was clipped to zero. PBE must reject the total, never repair it.
    let old = ScfProblem::new(&cell, &atoms, &species, options.clone()).unwrap();
    let density = old.starting_density();
    let total = &density + old.core_density().unwrap();
    let minimum = total.iter().copied().fold(f64::INFINITY, f64::min);
    assert!(minimum < 0.);
    let mut fft = dft_physics::data::FftPlan3d::new(old.grid_dimensions()).unwrap();
    assert!(Functional::Pbe.evaluate(&total, &cell, &mut fft).is_err());
    eprintln!("G2B_NLCC rejected total minimum={minimum:e} e/Bohr^3");
    for ranks in [1, 3, 8] {
        dft_mpi::simulate(ranks, |comm| {
            let _scope = screening.enter();
            let result = dft_engine::refstart::start(
                StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &species,
                    options: &options,
                },
                comm,
            );
            let Err(error) = result else {
                panic!("negative total density accepted");
            };
            assert!(error.to_string().contains("XC failed"), "{error}");
        });
    }
}

#[test]
fn g3_complex_k_smoke_and_kept_state() {
    use dft_engine::{ScfOptions, SpeciesPotential, refstart::StartInput};
    use dft_geometry::{Atom, Cell, kpoints::KPoint};
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    let cell = Cell::from_user_rows_bohr([[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: 1.,
        bands: Some(3),
        max_iterations: 100,
        ..Default::default()
    };
    let points = [
        KPoint {
            fractional: [0.13, -0.21, 0.07],
            weight: 0.4,
        },
        KPoint {
            fractional: [0.25; 3],
            weight: 0.6,
        },
    ];
    dft_mpi::simulate(3, |comm| {
        let mut previous: Option<f64> = None;
        for kept in ["off", "on"] {
            let config = dft_numerics::switches::RunConfig::from_pairs(&[]).unwrap();
            let _scope = config.enter();
            let run = dft_engine::refscf::k_driver::run_with(
                StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &options,
                },
                Profile::default(),
                &points,
                comm,
            )
            .unwrap();
            assert!(run.outcome.converged, "{}", run.exit_reason);
            assert!(run.outcome.result.verification.as_ref().unwrap().passed());
            assert_eq!(run.outcome.result.k_points.len(), 2);
            assert!((run.outcome.result.charge - 2.).abs() < 1e-12);
            assert!(
                run.records
                    .iter()
                    .all(|r| r.calls.len() == r.k_call_indices.len())
            );
            if kept == "on" {
                assert!(run.kept_hx.cache_hits > 0);
            }
            let energy = run.outcome.result.total_hartree();
            if let Some(old) = previous {
                assert!((energy - old).abs() < 1e-10);
            }
            previous = Some(energy);
            if comm.rank() == 0 {
                println!(
                    "G3 smoke kept={kept}: E={energy:.17e} Ha maps={} hits={}",
                    run.records.len(),
                    run.kept_hx.cache_hits
                );
            }
        }
    });
}

/// External Si UPF, never copied into the repository. Run in release:
/// DFT_SI_UPF=/absolute/Si.pz-vbc.UPF cargo test --release -p dft-engine
/// --test refscf g3_k_mesh_legacy_matrix -- --ignored --nocapture
#[test]
#[ignore = "external Si UPF and full SCF matrix; run the documented release command"]
fn g3_k_mesh_legacy_matrix() {
    use dft_engine::{
        GammaReal, OccupationPolicy, ScfOptions, ScfProblem, SpeciesPotential, refstart::StartInput,
    };
    use dft_geometry::{
        Cell,
        kpoints::{KPoint, monkhorst_pack, reduce_time_reversal},
    };
    let pp = dft_upf::parse_upf(
        &std::fs::read_to_string(std::env::var("DFT_SI_UPF").expect("DFT_SI_UPF required"))
            .unwrap(),
    )
    .unwrap();
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let atoms = common::atoms();
    let cubic = common::cell();
    let mut lattice = cubic.lattice_columns_bohr();
    lattice[0][1] += 0.17;
    lattice[1][2] -= 0.11;
    let skew = Cell::from_user_rows_bohr(std::array::from_fn(|i| {
        std::array::from_fn(|j| lattice[j][i])
    }))
    .unwrap();
    for (name, cell) in [("Si2", cubic), ("skew-Si2", skew)] {
        let meshes = [
            (
                "222",
                reduce_time_reversal(&monkhorst_pack([2; 3], [0.; 3]).unwrap()).unwrap(),
            ),
            (
                "222s",
                reduce_time_reversal(&monkhorst_pack([2; 3], [0.5; 3]).unwrap()).unwrap(),
            ),
            (
                "list",
                reduce_time_reversal(&[
                    KPoint {
                        fractional: [0.13, -0.21, 0.07],
                        weight: 0.2,
                    },
                    KPoint {
                        fractional: [-0.13, 0.21, -0.07],
                        weight: 0.2,
                    },
                    KPoint {
                        fractional: [0.25; 3],
                        weight: 0.6,
                    },
                ])
                .unwrap(),
            ),
        ];
        for (mesh, points) in meshes {
            let base = ScfOptions {
                cutoff_hartree: 4.,
                bands: Some(5),
                gamma_real: GammaReal::Off,
                max_iterations: 180,
                ..Default::default()
            };
            let strict = ScfProblem::with_kpoints(
                &cell,
                &atoms,
                &potentials,
                ScfOptions {
                    cutoff_hartree: base.cutoff_hartree,
                    bands: base.bands,
                    gamma_real: GammaReal::Off,
                    max_iterations: base.max_iterations,
                    free_electron_initial_guess: false,
                    ..ScfOptions::strict_reference()
                },
                &points,
                OccupationPolicy::Insulator {
                    minimum_gap_hartree: 1e-7,
                },
            )
            .unwrap()
            .run()
            .unwrap();
            assert!(strict.certified, "legacy {name}/{mesh}");
            let reference = strict.total_hartree();
            for ranks in [1, 3, 8] {
                let answers = dft_mpi::simulate(ranks, |comm| {
                    let mut answers = Vec::new();
                    for threshold_ry in [1e-8, 1e-12] {
                        let run = dft_engine::refscf::k_driver::run_with(
                            StartInput {
                                cell: &cell,
                                atoms: &atoms,
                                potentials: &potentials,
                                options: &base,
                            },
                            Profile {
                                threshold_ry,
                                ..Default::default()
                            },
                            &points,
                            comm,
                        )
                        .unwrap();
                        assert!(
                            run.outcome.converged,
                            "{name}/{mesh} P={ranks}: {}",
                            run.exit_reason
                        );
                        let r = &run.outcome.result;
                        assert!(r.verification.as_ref().unwrap().passed());
                        let energy = r.total_hartree();
                        let delta = (energy - reference).abs();
                        let cap = if threshold_ry == 1e-8 { 5e-9 } else { 1e-10 };
                        assert!(
                            delta < cap,
                            "{name}/{mesh} P={ranks} threshold={threshold_ry:e}: delta={delta:e} Ha"
                        );
                        assert!((r.charge - 8.).abs() < 1e-11);
                        assert!((energy - r.energy.band_form).abs() < 1e-9);
                        let mut eig = 0_f64;
                        let mut occupied = 0_f64;
                        for (k, want) in r.k_points.iter().zip(&strict.k_points) {
                            assert_eq!(k.point, want.point);
                            for (b, (a, z)) in k
                                .solution
                                .eigenvalues
                                .iter()
                                .zip(&want.solution.eigenvalues)
                                .enumerate()
                            {
                                eig = eig.max((a - z).abs());
                                if b < 4 {
                                    occupied = occupied.max((a - z).abs());
                                }
                            }
                        }
                        assert!(
                            (r.gap_hartree.unwrap() - strict.gap_hartree.unwrap()).abs()
                                <= 2. * eig + 1e-12,
                            "sampled indirect gap"
                        );
                        // Empty roots retain the reference 5e-6 Ha Ritz-change floor.
                        assert!(eig < 5e-5, "per-k eigenvalue error={eig:e} Ha");
                        assert!(
                            occupied < if threshold_ry == 1e-8 { 1e-5 } else { 1e-7 },
                            "occupied eigenvalue error={occupied:e} Ha"
                        );
                        if comm.rank() == 0 {
                            println!(
                                "G3 {name}/{mesh} P={ranks} threshold={threshold_ry:e} Ry E={energy:.17e} Ha legacy_delta={delta:e} Ha eig={eig:e} Ha occupied={occupied:e} Ha maps={}",
                                run.records.len()
                            );
                        }
                        answers.push(energy);
                    }
                    answers
                });
                for answer in &answers[1..] {
                    assert_eq!(answer, &answers[0]);
                }
            }
        }
    }
}

#[test]
fn g3_gamma_dispatch_is_bit_identical() {
    use dft_engine::{GammaReal, ScfOptions, SpeciesPotential, refstart::StartInput};
    use dft_geometry::{Atom, Cell, kpoints::KPoint};
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    let cell = Cell::from_user_rows_bohr([[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: 1.,
        bands: Some(3),
        gamma_real: GammaReal::On,
        ..Default::default()
    };
    for ranks in [1, 3, 8] {
        dft_mpi::simulate(ranks, |comm| {
            let input = || StartInput {
                cell: &cell,
                atoms: &atoms,
                potentials: &potentials,
                options: &options,
            };
            let a =
                dft_engine::refscf::driver::run_with(input(), Profile::default(), comm).unwrap();
            let b = dft_engine::refscf::k_driver::run_with(
                input(),
                Profile::default(),
                &[KPoint {
                    fractional: [0.; 3],
                    weight: 1.,
                }],
                comm,
            )
            .unwrap();
            assert_eq!(
                nlcc_rank_rule::fingerprint(&a),
                nlcc_rank_rule::fingerprint(&b)
            );
            assert_eq!(a.outcome.result.orbitals, b.outcome.result.orbitals);
        });
    }
}

/// The same test is injected into an archived base checkout for the G3 audit.
#[test]
fn g3_gamma_base_fingerprint() {
    use dft_engine::{GammaReal, ScfOptions, SpeciesPotential, refstart::StartInput};
    use dft_geometry::{Atom, Cell};
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    let cell = Cell::from_user_rows_bohr([[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: 1.,
        bands: Some(3),
        gamma_real: GammaReal::On,
        ..Default::default()
    };
    for ranks in [1, 3, 8] {
        dft_mpi::simulate(ranks, |comm| {
            let run = dft_engine::refscf::driver::run_with(
                StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &options,
                },
                Profile::default(),
                comm,
            )
            .unwrap();
            assert!(run.outcome.converged);
            // Normalize only the additive empty Gamma k-index metadata. Every
            // pre-existing scientific record and all numeric bits remain hashed.
            let mut fingerprint = nlcc_rank_rule::fingerprint(&run)
                .replace("k_call_indices: [], ", "")
                .into_bytes();
            for z in run.outcome.result.orbitals.as_slice() {
                fingerprint.extend_from_slice(&z.re.to_le_bytes());
                fingerprint.extend_from_slice(&z.im.to_le_bytes());
            }
            if comm.rank() == 0
                && let Some(dir) = std::env::var_os("G3_GAMMA_DUMP_DIR")
            {
                std::fs::write(
                    std::path::Path::new(&dir).join(format!("rank-{ranks}.bin")),
                    &fingerprint,
                )
                .unwrap();
            }
            let hash = fingerprint.iter().fold(0xcbf29ce484222325_u64, |h, &b| {
                (h ^ u64::from(b)).wrapping_mul(0x100000001b3)
            });
            if comm.rank() == 0 {
                println!(
                    "G3_GAMMA P={ranks} hash={hash:016x} E={:.17e} Ha",
                    run.outcome.result.total_hartree()
                );
            }
        });
    }
}

#[test]
fn g3_k_start_does_not_require_a_gamma_band_space() {
    use dft_engine::{ScfOptions, SpeciesPotential, refstart::StartInput};
    use dft_geometry::{Atom, Cell, kpoints::KPoint};
    let pp = common::potential_with_atomic_start(2., 1.8, false, 161, "SLA PZ");
    let cell = Cell::from_user_rows_bohr([[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: 0.35,
        bands: Some(3),
        ..Default::default()
    };
    assert_eq!(
        dft_basis::PlaneWaveBasis::new(&cell, dft_units::Hartree(0.35))
            .unwrap()
            .vectors()
            .len(),
        1
    );
    dft_mpi::simulate(8, |comm| {
        let state = dft_engine::refstart::start_kpoints(
            StartInput {
                cell: &cell,
                atoms: &atoms,
                potentials: &potentials,
                options: &options,
            },
            &[KPoint {
                fractional: [0.25; 3],
                weight: 1.,
            }],
            comm,
        )
        .unwrap();
        assert_eq!(state.orbitals.cols(), 3);
        assert!(state.blocks[0].basis.vectors().len() >= 3);
    });
}

#[test]
fn g4_spin_metric_uniform_and_unpolarized_limits() {
    for ranks in [1, 3, 8] {
        dft_mpi::simulate(ranks, |comm| {
            let global = [0., 2., 5., 9.];
            let indices: Vec<_> = (0..4).filter(|i| i % ranks == comm.rank()).collect();
            let g2: Vec<_> = indices.iter().map(|&i| global[i]).collect();
            let mask = vec![true; g2.len()];
            let mut metric = Metric::spin(comm, 17., &g2, &mask, 8).unwrap();
            let mut scalar = Metric::new(17., &g2, &mask, ranks, 8).unwrap();
            let charge: Vec<_> = indices
                .iter()
                .map(|&i| c(if i == 0 { 0. } else { 0.01 * i as f64 }))
                .collect();
            let mut stacked = charge.clone();
            stacked.extend(vec![c(0.); g2.len()]);
            assert!(
                (metric.norm(comm, &stacked).unwrap() - scalar.norm(comm, &charge).unwrap()).abs()
                    < 1e-14
            );
            stacked.fill(c(0.));
            for (j, &i) in indices.iter().enumerate() {
                if i == 0 {
                    stacked[g2.len() + j] = c(0.03);
                }
            }
            let expected = 4. * PI * 17. / 2. * 0.03_f64.powi(2);
            assert!((metric.norm(comm, &stacked).unwrap() - expected).abs() < 1e-14);
        });
    }
}

#[test]
fn g4_scalar_fd_gamma_legacy() {
    assert!(g4_fixture("scalar", false, &[1], &[1e-12]).is_empty());
}

#[test]
#[ignore = "release matrix: cargo test --release -p dft-engine --test refscf g4_synthetic_matrix -- --ignored --nocapture"]
fn g4_synthetic_matrix() {
    let mut failures = Vec::new();
    for mode in ["scalar", "fixed"] {
        for mesh in [false, true] {
            failures.extend(g4_fixture(mode, mesh, &[1, 3, 8], &[1e-8, 1e-12]));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
#[test]
// R6a qualifies moment policies with registry defaults; keep the long matrix in release.
#[ignore = "release matrix: cargo test --release -p dft-engine --test refscf g4b_moment_acceptance -- --ignored --nocapture"]
fn g4b_moment_acceptance() {
    let mut failures = Vec::new();
    for mode in ["free", "moment"] {
        for mesh in [false, true] {
            failures.extend(g4_fixture(mode, mesh, &[1, 3, 8], &[1e-8, 1e-12]));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

fn g4_fixture(mode: &str, mesh: bool, ranks: &[usize], thresholds: &[f64]) -> Vec<String> {
    g4_fixture_profile(mode, mesh, ranks, thresholds, Profile::default())
}

#[test]
#[ignore = "diagnostic only: cargo test --release -p dft-engine --test refscf g4_free_moment_tight_short_history -- --ignored --nocapture"]
fn g4_free_moment_tight_short_history() {
    // This diagnoses the default-history failures, not a replacement for the
    // default-profile acceptance matrix or a change of the accuracy bounds.
    let mut failures = Vec::new();
    for mesh in [false, true] {
        failures.extend(g4_fixture_profile(
            "free",
            mesh,
            &[1, 3, 8],
            &[1e-12],
            Profile {
                history: 4,
                beta: 0.5,
                ..Default::default()
            },
        ));
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

fn g4_fixture_profile(
    mode: &str,
    mesh: bool,
    ranks: &[usize],
    thresholds: &[f64],
    profile: Profile,
) -> Vec<String> {
    // Historical fixture contract; screening has separate charge/mask tests.
    let screening =
        dft_numerics::switches::RunConfig::from_pairs(&[("xc_screen_negative", "off")]).unwrap();
    let _screen_scope = screening.enter();
    use dft_engine::{
        Functional, GammaReal, ScfOptions, SpeciesPotential, StartDensity,
        collinear::{CollinearProblem, SpinSector},
        finite_temperature::FiniteTProblem,
        refscf::electronic::{self, Policy},
        refstart::StartInput,
    };
    use dft_geometry::{
        Atom, Cell,
        kpoints::{KPoint, monkhorst_pack},
    };
    use dft_physics::occupations::Constraint;
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0.2, 5., 0.], [0., 0.1, 5.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    let functional = if mode == "scalar" {
        Functional::LdaPz
    } else {
        Functional::LdaPw
    };
    let pp = common::potential_with_atomic_start(
        1.,
        1.,
        true,
        161,
        if mode == "scalar" { "SLA PZ" } else { "SLA PW" },
    );
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let points = if mesh {
        monkhorst_pack([2; 3], [0.; 3]).unwrap()
    } else {
        vec![KPoint {
            fractional: [0.; 3],
            weight: 1.,
        }]
    };
    let sector = SpinSector {
        n_up: 1,
        n_down: 0,
        axis_cartesian: [0., 0., 1.],
        initial_moment_seed: Some(0.3),
    };
    let policy = match mode {
        "scalar" => Policy::FermiDirac {
            theta_ha: 0.01,
            band_cap: 12,
        },
        "fixed" => Policy::FixedSpin(sector.clone()),
        _ => Policy::SpinFermiDirac {
            sector: sector.clone(),
            theta_ha: 0.01,
            constraint: if mode == "moment" {
                Constraint::FixedMoment(0.3)
            } else {
                Constraint::Total
            },
            band_cap: 12,
        },
    };
    let options = ScfOptions {
        functional,
        cutoff_hartree: 2.5,
        bands: Some(if mode == "fixed" { 4 } else { 12 }),
        gamma_real: GammaReal::Off,
        start_density: StartDensity::Atomic,
        max_iterations: 400,
        // The legacy mixer imposes the output channel means even on a damped
        // proposal. At free moment that can make a minority channel negative.
        // A full linear step preserves the positive FD output; the reference
        // driver's Broyden policy and all certificate tolerances stay fixed.
        mixing: if mode == "free" {
            dft_engine::Mixing::Linear { beta: 1.0 }
        } else {
            dft_engine::Mixing::default()
        },
        residual_eigen_schedule: false,
        initial_eigen_tolerance: 1e-10,
        eigen_tolerance: 1e-10,
        density_tolerance: 1e-10,
        energy_tolerance: 1e-11,
        ..ScfOptions::strict_reference()
    };
    let (oracle_energy, oracle_moment, oracle_mu, oracle_occupations) = if mode == "scalar" {
        let r = FiniteTProblem::new(
            &cell,
            &atoms,
            &potentials,
            options.clone(),
            &points,
            0.01,
            12,
        )
        .unwrap()
        .run()
        .unwrap();
        assert!(r.converged);
        (
            r.thermodynamics.free_energy_ha,
            0.,
            r.occupations
                .evaluation
                .unwrap()
                .roots
                .map(|r| r.map(|r| r.mu_ha)),
            Some(r.occupations.occupations().to_vec()),
        )
    } else {
        let p = if mode == "fixed" {
            CollinearProblem::new(&cell, &atoms, &potentials, options.clone(), &points, sector)
        } else {
            CollinearProblem::new_fermi_dirac(
                &cell,
                &atoms,
                &potentials,
                options.clone(),
                &points,
                sector,
                0.01,
                if mode == "moment" {
                    Constraint::FixedMoment(0.3)
                } else {
                    Constraint::Total
                },
                12,
            )
        };
        let r = p.unwrap().run().unwrap();
        assert!(
            r.converged,
            "legacy {mode} mesh={mesh} moment={} last={:?} cert={:?}",
            r.moment_mu_b,
            r.history.last(),
            r.finite_certificate
        );
        (
            r.thermodynamics
                .map_or(r.energy.total_hartree(), |e| e.free_energy_ha),
            r.moment_mu_b,
            r.finite_occupations.as_ref().map_or([None; 2], |f| {
                f.evaluation.unwrap().roots.map(|r| r.map(|r| r.mu_ha))
            }),
            r.finite_occupations
                .as_ref()
                .map(|f| f.occupations().to_vec()),
        )
    };
    let options = ScfOptions {
        standard: Some(dft_engine::StandardProfile { threshold_ha: 5e-9 }),
        verification: dft_engine::Verification::FixedPoint,
        ..options
    };
    let mut failures = Vec::new();
    for &threshold in thresholds {
        for &rank in ranks {
            let results = dft_mpi::simulate(rank, |comm| {
                let _screen_scope = screening.enter();
                let trace_dir = std::env::var_os("G4B_TRACE_DIR");
                let mut trace = electronic::Trace::default();
                let result = electronic::run_with_trace(
                    StartInput {
                        cell: &cell,
                        atoms: &atoms,
                        potentials: &potentials,
                        options: &options,
                    },
                    Profile {
                        threshold_ry: threshold,
                        ..profile
                    },
                    &points,
                    policy.clone(),
                    comm,
                    trace_dir.as_ref().map(|_| &mut trace),
                );
                if comm.rank() == 0
                    && let Some(dir) = trace_dir
                {
                    let mut text = String::from(
                        "iteration,verification,redo,charge_Ry,spin_Ry,uniform_Ry,dF_Ha,dM_muB,dmu_Ha,docc_e,number_residual_e,moment_residual_e,condition,max_pair_cosine,pair_norm2_Ry,stiffness_Ha,dropped_pairs,dropped_conditions\n",
                    );
                    use std::fmt::Write;
                    for t in trace.records() {
                        let mu = t
                            .chemical_potentials_ha
                            .iter()
                            .zip(oracle_mu)
                            .filter_map(|(a, b)| a.zip(b))
                            .map(|(a, b)| (a - b).abs())
                            .fold(0., f64::max);
                        let docc = oracle_occupations.as_ref().map_or(0., |old| {
                            t.occupations
                                .iter()
                                .zip(old)
                                .enumerate()
                                .map(|(i, (a, b))| {
                                    points[(i / 12) % points.len()].weight * (a - b).abs()
                                })
                                .sum::<f64>()
                        });
                        writeln!(
                            text,
                            "{},{},{},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{:e},\"{:?}\",{:e},{},\"{:?}\"",
                            t.iteration,
                            t.verification,
                            t.redo,
                            t.estimator_terms_ry[0],
                            t.estimator_terms_ry[1],
                            t.estimator_terms_ry[2],
                            (t.free_energy_ha - oracle_energy).abs(),
                            (t.moment_mu_b - oracle_moment).abs(),
                            mu,
                            docc,
                            t.electron_residual,
                            t.moment_residual.unwrap_or(0.),
                            t.condition.unwrap_or(0.),
                            t.history_max_cosine,
                            t.history_pair_norm2_ry,
                            t.stiffness_ha,
                            t.dropped_pairs,
                            t.dropped_conditions
                        )
                        .unwrap();
                    }
                    std::fs::write(
                        std::path::Path::new(&dir)
                            .join(format!("{mode}-{mesh}-{rank}-{threshold:e}.csv")),
                        text,
                    )
                    .unwrap();
                }
                let r = match result {
                    Ok(r) => r,
                    Err(e) => {
                        return Some(format!(
                            "{mode} mesh={mesh} ranks={rank} threshold={threshold:e} Ry: {e}"
                        ));
                    }
                };
                if !r.converged {
                    return Some(format!(
                        "{mode} mesh={mesh} ranks={rank} threshold={threshold:e} Ry: exhausted"
                    ));
                }
                if comm.rank() == 0 {
                    g4b_dump(&r, format!("{mode}-{mesh}-{rank}-{threshold:e}.txt"));
                }
                assert!(r.verification.as_ref().unwrap().passed());
                let delta = (r.thermodynamics.free_energy_ha - oracle_energy).abs();
                let moment = (r.moment_mu_b - oracle_moment).abs();
                let mu = r.occupations.map_or(0., |e| {
                    e.roots
                        .iter()
                        .zip(oracle_mu)
                        .filter_map(|(a, b)| a.zip(b))
                        .map(|(a, b)| (a.mu_ha - b).abs())
                        .fold(0., f64::max)
                });
                let occupation_change: f64 = oracle_occupations.as_ref().map_or(0., |old| {
                    let nb = r.blocks[0][0].occupations.len();
                    r.blocks
                        .iter()
                        .flatten()
                        .flat_map(|b| &b.occupations)
                        .zip(old)
                        .enumerate()
                        .map(|(i, (a, b))| points[(i / nb) % points.len()].weight * (a - b).abs())
                        .sum()
                });
                if comm.rank() == 0 {
                    println!(
                        "G4 mode={mode} mesh={mesh} ranks={rank} threshold={threshold:e} Ry F={:.17e} Ha dF={delta:e} Ha dM={moment:e} mu_B dmu={mu:e} Ha docc={occupation_change:e} e iterations={} M={:.12e} mu_B verification_attempts={}",
                        r.thermodynamics.free_energy_ha,
                        r.history.len(),
                        r.moment_mu_b,
                        r.verification.as_ref().unwrap().attempts.len()
                    );
                }
                assert!((r.channel_electrons.iter().sum::<f64>() - 1.).abs() < 1e-10);
                // Collect accuracy failures after all collectives so one bad
                // cell does not hide the remaining rank/threshold evidence.
                (delta > if threshold < 1e-8 { 1e-10 } else { 5e-9 }
                    || moment >= 1e-4 || mu >= 1e-4 || occupation_change > 1e-4).then(|| format!(
                    "{mode} mesh={mesh} ranks={rank} threshold={threshold:e} Ry: dF={delta:e} Ha dM={moment:e} mu_B dmu={mu:e} Ha docc={occupation_change:e} e"
                ))
            });
            failures.extend(results.into_iter().take(1).flatten());
        }
    }
    failures
}

#[test]
#[ignore = "external UPF/release: DFT_SI_UPF=/path/Si.pz-vbc.UPF cargo test --release -p dft-engine --test refscf g4_si_gap -- --ignored --nocapture"]
fn g4_si_gap() {
    use dft_engine::{
        GammaReal, OccupationPolicy, ScfOptions, ScfProblem, SpeciesPotential,
        finite_temperature::FiniteTProblem,
        refscf::electronic::{self, Policy},
        refstart::StartInput,
    };
    use dft_geometry::kpoints::{KPoint, monkhorst_pack};
    let pp = dft_upf::parse_upf(
        &std::fs::read_to_string(std::env::var("DFT_SI_UPF").expect("DFT_SI_UPF required"))
            .unwrap(),
    )
    .unwrap();
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let cell = common::cell();
    let atoms = common::atoms();
    for mesh in [false, true] {
        let points = if mesh {
            monkhorst_pack([2; 3], [0.; 3]).unwrap()
        } else {
            vec![KPoint {
                fractional: [0.; 3],
                weight: 1.,
            }]
        };
        let tight = ScfOptions {
            cutoff_hartree: 4.,
            bands: Some(8),
            gamma_real: GammaReal::Off,
            max_iterations: 200,
            eigen_tolerance: 1e-10,
            density_tolerance: 1e-10,
            energy_tolerance: 1e-11,
            ..ScfOptions::strict_reference()
        };
        let integer = ScfProblem::with_kpoints(
            &cell,
            &atoms,
            &potentials,
            tight.clone(),
            &points,
            OccupationPolicy::Insulator {
                minimum_gap_hartree: 1e-7,
            },
        )
        .unwrap()
        .run()
        .unwrap();
        let finite = FiniteTProblem::new(&cell, &atoms, &potentials, tight, &points, 1e-4, 8)
            .unwrap()
            .run()
            .unwrap();
        assert!(integer.certified && finite.converged);
        assert!((integer.total_hartree() - finite.thermodynamics.free_energy_ha).abs() < 1e-10);
        let options = ScfOptions {
            cutoff_hartree: 4.,
            bands: Some(8),
            gamma_real: GammaReal::Off,
            max_iterations: 200,
            ..Default::default()
        };
        for threshold in [1e-8, 1e-12] {
            for ranks in [1, 3, 8] {
                dft_mpi::simulate(ranks, |comm| {
                    let run = electronic::run_with(
                        StartInput {
                            cell: &cell,
                            atoms: &atoms,
                            potentials: &potentials,
                            options: &options,
                        },
                        Profile {
                            threshold_ry: threshold,
                            ..Default::default()
                        },
                        &points,
                        Policy::FermiDirac {
                            theta_ha: 1e-4,
                            band_cap: 8,
                        },
                        comm,
                    )
                    .unwrap();
                    if comm.rank() == 0 {
                        g4b_dump(&run, format!("si-{mesh}-{ranks}-{threshold:e}.txt"));
                    }
                    assert!(run.converged && run.verification.as_ref().unwrap().passed());
                    let delta_integer =
                        (run.thermodynamics.free_energy_ha - integer.total_hartree()).abs();
                    let delta_fd = (run.thermodynamics.free_energy_ha
                        - finite.thermodynamics.free_energy_ha)
                        .abs();
                    let limit = if threshold == 1e-8 { 5e-9 } else { 1e-10 };
                    assert!(
                        delta_integer <= limit && delta_fd <= limit,
                        "Si mesh={mesh} ranks={ranks} d_integer={delta_integer:e} d_fd={delta_fd:e}"
                    );
                    if comm.rank() == 0 {
                        println!(
                            "G4 Si mesh={mesh} ranks={ranks} threshold={threshold:e} Ry F={:.17e} Ha d_integer={delta_integer:e} Ha d_fd={delta_fd:e} Ha TS={:e} Ha",
                            run.thermodynamics.free_energy_ha, run.thermodynamics.entropy_term_ha
                        );
                    }
                });
            }
        }
    }
}

#[test]
fn g4_rejections_and_unverified_delivery() {
    use dft_engine::{
        Functional, ScfOptions, SpeciesPotential, Verification,
        collinear::SpinSector,
        refscf::electronic::{self, Policy},
        refstart::StartInput,
    };
    use dft_geometry::{Atom, Cell, kpoints::KPoint};
    use dft_physics::occupations::Constraint;
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0., 5., 0.], [0., 0., 5.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.; 3])];
    let pp = common::potential_with_atomic_start(1., 1., false, 161, "SLA PZ");
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let points = [KPoint {
        fractional: [0.; 3],
        weight: 1.,
    }];
    let options = ScfOptions {
        cutoff_hartree: 2.5,
        bands: Some(12),
        verification: Verification::Off,
        ..Default::default()
    };
    let run = |options: &ScfOptions, policy| {
        electronic::run_with(
            StartInput {
                cell: &cell,
                atoms: &atoms,
                potentials: &potentials,
                options,
            },
            Profile::default(),
            &points,
            policy,
            &Serial,
        )
    };
    let good = run(
        &options,
        Policy::FermiDirac {
            theta_ha: 0.01,
            band_cap: 12,
        },
    )
    .unwrap();
    assert!(good.converged && good.verification.is_none());
    assert!(good.blocks[0][0].occupations[0] > 0.99 && good.blocks[0][0].occupations[0] < 1.01);
    let verified = run(
        &ScfOptions {
            verification: Verification::FixedPoint,
            ..options.clone()
        },
        Policy::FermiDirac {
            theta_ha: 0.01,
            band_cap: 12,
        },
    )
    .unwrap();
    assert!(verified.verification.as_ref().unwrap().passed());
    assert_eq!(verified.history.len(), good.history.len());
    assert_eq!(verified.density, good.density);
    assert_eq!(
        verified.thermodynamics.free_energy_ha.to_bits(),
        good.thermodynamics.free_energy_ha.to_bits()
    );
    assert_eq!(
        verified.blocks[0][0].eigenvalues_ha,
        good.blocks[0][0].eigenvalues_ha
    );
    assert_eq!(
        verified.blocks[0][0].occupations,
        good.blocks[0][0].occupations
    );
    assert_eq!(
        verified.blocks[0][0].orbitals.as_slice(),
        good.blocks[0][0].orbitals.as_slice()
    );
    assert_eq!(
        verified.occupations.unwrap().roots[0]
            .unwrap()
            .mu_ha
            .to_bits(),
        good.occupations.unwrap().roots[0].unwrap().mu_ha.to_bits()
    );
    assert!(
        run(
            &options,
            Policy::FermiDirac {
                theta_ha: 0.,
                band_cap: 12
            }
        )
        .is_err()
    );
    assert!(
        run(
            &options,
            Policy::FermiDirac {
                theta_ha: 0.5,
                band_cap: 12
            }
        )
        .is_err()
    );
    let sector = SpinSector {
        n_up: 1,
        n_down: 0,
        axis_cartesian: [0., 0., 1.],
        initial_moment_seed: None,
    };
    assert!(run(&options, Policy::FixedSpin(sector.clone())).is_err());
    let pw = ScfOptions {
        functional: Functional::LdaPw,
        allow_functional_override: true,
        ..options.clone()
    };
    assert!(
        run(
            &pw,
            Policy::SpinFermiDirac {
                sector: sector.clone(),
                theta_ha: 0.01,
                constraint: Constraint::FixedMoment(1.),
                band_cap: 12
            }
        )
        .is_err()
    );
    let pbe = ScfOptions {
        functional: Functional::Pbe,
        ..pw
    };
    assert!(run(&pbe, Policy::FixedSpin(sector)).unwrap().converged);
}

#[test]
fn g4b_conditioning_preserves_newest_pair_and_off_rejects() {
    let profile = Profile {
        history: 3,
        ..Default::default()
    };
    let mut metric = Metric::new(1., &[1., 1.], &[true; 2], 1, 350).unwrap();
    let mut mix = Broyden::new(profile, &[1., 1.], &[true; 2], 1., 1.).unwrap();
    mix.condition_history();
    let x0 = [C::new(0., 0.); 2];
    let x1 = [C::new(1., 0.), C::new(2., 0.)];
    let x2 = [C::new(2., 0.), C::new(3., 0.)];
    let f0 = [C::new(1., 0.), C::new(0., 0.)];
    let f1 = [C::new(0., 0.); 2];
    let f2 = [C::new(-1., 0.), C::new(1e-6, 0.)];
    let mut next = [C::default(); 2];
    mix.update(&Serial, &mut metric, &x0, &f0, &mut next)
        .unwrap();
    mix.update(&Serial, &mut metric, &x1, &f1, &mut next)
        .unwrap();
    let d = mix
        .update(&Serial, &mut metric, &x2, &f2, &mut next)
        .unwrap();
    assert_eq!(
        (mix.conditioning_diagnostics().dropped_pairs, d.occupancy),
        (1, 1)
    );
    assert!(mix.conditioning_diagnostics().dropped_conditions[0] > 1e12);
    assert_eq!(d.condition, Some(1.));
    // The remaining newest df=-f2, so gamma=-1 and x_tilde=x1.
    for i in 0..2 {
        assert!((next[i] - x1[i]).norm() < 1e-13);
    }
    // Exact duplicate/zero pairs are also removed down to the linear step.
    let d = mix
        .update(&Serial, &mut metric, &x2, &f2, &mut next)
        .unwrap();
    assert_eq!(
        (mix.conditioning_diagnostics().dropped_pairs, d.occupancy),
        (2, 0)
    );
    for i in 0..2 {
        assert_eq!(next[i], x2[i] + profile.beta * f2[i]);
    }
    let mut legacy = Broyden::new(profile, &[1., 1.], &[true; 2], 1., 1.).unwrap();
    legacy
        .update(&Serial, &mut metric, &x2, &f2, &mut next)
        .unwrap();
    assert!(matches!(
        legacy.update(&Serial, &mut metric, &x2, &f2, &mut next),
        Err(Error::Singular { .. })
    ));
}

// Stable complete scientific result stream, excluding wall clocks only.
fn g4b_fingerprint(run: &dft_engine::refscf::electronic::Run) -> String {
    use std::fmt::Write;
    let mut out = format!(
        "{:?} {:?} {:?} {:?} {:?} {:?} {:?}",
        run.converged,
        run.energy,
        run.thermodynamics,
        run.occupations,
        run.channel_electrons,
        run.moment_mu_b,
        run.absolute_moment_mu_b
    );
    for density in &run.density {
        for x in density {
            write!(out, " {:016x}", x.to_bits()).unwrap();
        }
    }
    for b in run.blocks.iter().flatten() {
        for x in b
            .eigenvalues_ha
            .iter()
            .chain(&b.occupations)
            .chain(&b.residuals_ha)
        {
            write!(out, " {:016x}", x.to_bits()).unwrap();
        }
        for x in b.orbitals.as_slice() {
            write!(out, " {:016x}{:016x}", x.re.to_bits(), x.im.to_bits()).unwrap();
        }
    }
    for h in &run.history {
        let mut calls = h.calls.clone();
        for (_, _, c) in &mut calls {
            c.inner_clocks = None;
        }
        write!(
            out,
            " {:?} {:?}",
            (
                h.iteration,
                h.dr2_ry,
                h.ethr_ry,
                h.first_iteration_redo,
                h.energy_ha,
                h.entropy_term_ha,
                h.free_energy_ha,
                h.moment_mu_b,
                h.residual_guard_ha,
                h.mixing,
                h.history_conditioning,
                h.mixing_fraction
            ),
            calls
        )
        .unwrap();
    }
    let mut verification = run.verification.clone();
    if let Some(v) = &mut verification {
        for a in &mut v.attempts {
            a.wall = std::time::Duration::ZERO;
        }
    }
    let mut calls = run.verification_calls.clone();
    for (_, _, c) in &mut calls {
        c.inner_clocks = None;
    }
    write!(out, " {verification:?} {calls:?}").unwrap();
    out
}

fn g4b_dump(run: &dft_engine::refscf::electronic::Run, name: String) {
    if let Some(dir) = std::env::var_os("G4B_FINGERPRINT_DIR") {
        std::fs::write(std::path::Path::new(&dir).join(name), g4b_fingerprint(run)).unwrap();
    }
}

// G6: synthetic checks are internal oracles, never QE evidence.
#[test]
fn g6_symmetry_geometry_guard_and_small_map() {
    use dft_engine::{
        ScfOptions, SpeciesPotential,
        refscf::{k_driver, symmetry::Prepared},
        refstart::StartInput,
    };
    use dft_geometry::{Atom, Cell, kpoints::monkhorst_pack};
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0.2, 5., 0.], [0., 0.1, 5.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.25; 3])];
    let full = monkhorst_pack([2; 3], [0.5; 3]).unwrap();
    let prepared = Prepared::new(&cell, &atoms, &full, 1e-7, true).unwrap();
    let pp = common::potential_with_atomic_start(2., 1.8, true, 81, "SLA PZ");
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: 1.,
        bands: Some(3),
        max_iterations: 1,
        verification: dft_engine::Verification::Off,
        ..Default::default()
    };
    let spin_tr = dft_engine::refscf::electronic::run_prepared(
        StartInput {
            cell: &cell,
            atoms: &atoms,
            potentials: &potentials,
            options: &options,
        },
        Profile::default(),
        &prepared,
        dft_engine::refscf::electronic::Policy::FixedSpin(dft_engine::collinear::SpinSector {
            n_up: 1,
            n_down: 1,
            axis_cartesian: [0., 0., 1.],
            initial_moment_seed: None,
        }),
    );
    assert!(
        matches!(spin_tr, Err(dft_engine::ScfError::InvalidInput(message)) if message.contains("disable time reversal"))
    );
    let moved = [Atom::new("Xx", [0.125; 3])];
    assert!(
        k_driver::run_prepared(
            StartInput {
                cell: &cell,
                atoms: &moved,
                potentials: &potentials,
                options: &options
            },
            Profile::default(),
            &prepared
        )
        .is_err()
    );
    // Clipping an under-resolved atomic start introduces out-of-sphere
    // components. Both drivers must reject it at the same invariant check.
    let narrow = common::potential_with_atomic_start(2., 1., true, 81, "SLA PZ");
    let narrow_potentials = [SpeciesPotential {
        label: "Xx",
        potential: &narrow,
    }];
    let narrow_options = ScfOptions {
        grid_override: Some([8; 3]),
        start_density: dft_engine::StartDensity::Atomic,
        ..options.clone()
    };
    let bad = k_driver::run_prepared(
        StartInput {
            cell: &cell,
            atoms: &atoms,
            potentials: &narrow_potentials,
            options: &narrow_options,
        },
        Profile::default(),
        &prepared,
    );
    assert!(
        matches!(bad,Err(dft_engine::ScfError::InvalidInput(message)) if message.contains("symmetry.density"))
    );
    let mut old = dft_engine::ScfProblem::with_kpoints(
        &cell,
        &atoms,
        &narrow_potentials,
        narrow_options,
        prepared.points(),
        dft_engine::OccupationPolicy::Insulator {
            minimum_gap_hartree: 1e-7,
        },
    )
    .unwrap();
    old.enable_symmetry(&full, 1e-7, true).unwrap();
    assert!(
        matches!(old.run(),Err(dft_engine::ScfError::InvalidInput(message)) if message.contains("symmetry.density"))
    );
    for ranks in [1, 3, 8] {
        dft_mpi::simulate(ranks, |comm| {
            let run = k_driver::run_prepared_with(
                StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &options,
                },
                Profile::default(),
                &prepared,
                comm,
            )
            .unwrap();
            let report = run.symmetry.unwrap();
            assert_eq!(report.group_order, 2);
            assert_eq!(report.irreducible_k_count, 4);
            assert!(!report.maps.is_empty());
            assert!(
                report
                    .maps
                    .iter()
                    .all(|r| r.norm.is_finite() && r.iteration == 1 && !r.verification)
            );
            assert!((run.outcome.result.charge - 2.).abs() < 1e-12);
        });
    }
}

/// External UPF stays outside git. Full mesh, reduced reference and strict
/// reduced legacy share cutoff/grid/start. Release-only (ADR 0004):
/// DFT_SI_UPF=/absolute/Si.pz-vbc.UPF cargo test --release -p dft-engine --test refscf g6_si444 -- --ignored --nocapture
#[test]
#[ignore = "Si 4x4x4 SCF matrix; release command above"]
fn g6_si444_symmetry_full_legacy_and_derivatives() {
    use dft_engine::{
        GammaReal, OccupationPolicy, PreparedSymmetry, ScfOptions, ScfProblem, SpeciesPotential,
        StartDensity,
        refscf::{delivered, k_driver, symmetry::Prepared},
        refstart::StartInput,
    };
    use dft_geometry::{Cell, kpoints::monkhorst_pack};
    let pp = dft_upf::parse_upf(
        &std::fs::read_to_string(std::env::var("DFT_SI_UPF").expect("DFT_SI_UPF required"))
            .unwrap(),
    )
    .unwrap();
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let atoms = common::atoms();
    let cubic = common::cell();
    let mut a = cubic.lattice_columns_bohr();
    a[0][1] += 0.17;
    a[1][2] -= 0.11;
    let skew = Cell::from_user_rows_bohr(std::array::from_fn(|i| std::array::from_fn(|j| a[j][i])))
        .unwrap();
    for (name, cell) in [("diamond", cubic), ("skew", skew)] {
        let full = monkhorst_pack([4; 3], [0.; 3]).unwrap();
        let legacy_prepared = PreparedSymmetry::new(&cell, &atoms, &full, 1e-7, true).unwrap();
        let prepared =
            Prepared::from_prepared(&cell, &atoms, &full, 1e-7, true, &legacy_prepared).unwrap();
        let operations = prepared.operations(false);
        assert!(
            operations
                .iter()
                .any(|o| o.translation.iter().any(|t| *t != 0.))
        );
        if name == "diamond" {
            assert_eq!(operations.len(), 48);
        } else {
            assert!(operations.len() <= 4);
        }
        let b = dft_basis::PlaneWaveBasis::new(&cell, dft_units::Hartree(4.)).unwrap();
        let db = dft_basis::PlaneWaveBasis::new(&cell, dft_units::Hartree(16.)).unwrap();
        let grid = dft_geometry::symmetry::compatible_fft_grid(
            dft_basis::FftGridSpec::new(&b, &db, None).unwrap().dims(),
            operations,
        )
        .unwrap();
        let options = ScfOptions {
            cutoff_hartree: 4.,
            bands: Some(5),
            gamma_real: GammaReal::Off,
            start_density: StartDensity::Atomic,
            max_iterations: 200,
            grid_override: Some(grid),
            ..Default::default()
        };
        let strict = ScfOptions {
            standard: None,
            verification: dft_engine::Verification::Off,
            eigen_tolerance: 1e-11,
            density_tolerance: 1e-11,
            energy_tolerance: 1e-12,
            ..options.clone()
        };
        let mut legacy = ScfProblem::with_kpoints(
            &cell,
            &atoms,
            &potentials,
            strict.clone(),
            prepared.points(),
            OccupationPolicy::Insulator {
                minimum_gap_hartree: 1e-7,
            },
        )
        .unwrap();
        legacy.enable_prepared_symmetry(&legacy_prepared).unwrap();
        let oracle = legacy.run().unwrap();
        assert!(oracle.certified);
        let view = dft_engine::ForceState {
            cell: &cell,
            atoms: &atoms,
            potentials: &potentials,
            options: &strict,
            electronic: dft_engine::ElectronicState::Fixed(&oracle),
            density_floor: None,
            symmetry: operations,
        };
        let old_forces = dft_engine::forces(&view).unwrap().values;
        let old_stress = dft_engine::stress::frozen_stress(&view).unwrap().values;
        for threshold_ry in [1e-8, 1e-12] {
            let profile = Profile {
                threshold_ry,
                ..Default::default()
            };
            let full_run = k_driver::run_with(
                StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &options,
                },
                profile,
                &full,
                &Serial,
            )
            .unwrap();
            assert!(
                full_run.outcome.converged,
                "full {name}: {}",
                full_run.exit_reason
            );
            let full_result = &full_run.outcome.result;
            for ranks in [1, 3, 8] {
                dft_mpi::simulate(ranks, |comm| {
                    let mut run = k_driver::run_prepared_with(
                        StartInput {
                            cell: &cell,
                            atoms: &atoms,
                            potentials: &potentials,
                            options: &options,
                        },
                        profile,
                        &prepared,
                        comm,
                    )
                    .unwrap();
                    assert!(run.outcome.converged, "{name}/{ranks}: {}", run.exit_reason);
                    let r = &run.outcome.result;
                    assert!(r.verification.as_ref().unwrap().passed());
                    let de = (r.total_hartree() - full_result.total_hartree()).abs();
                    let dl = (r.total_hartree() - oracle.total_hartree()).abs();
                    let dn = r
                        .density
                        .iter()
                        .zip(&full_result.density)
                        .map(|(a, b)| (a - b).abs())
                        .fold(0., f64::max);
                    let cap = if threshold_ry == 1e-8 { 5e-9 } else { 1e-10 };
                    assert!(
                        de < cap && dl < cap,
                        "{name}/{ranks}: full={de:e}, legacy={dl:e} Ha"
                    );
                    assert!(
                        dn < if threshold_ry == 1e-8 { 1e-5 } else { 1e-7 },
                        "density={dn:e} e/Bohr^3"
                    );
                    let report = run.symmetry.as_ref().unwrap();
                    assert!(report.maps.iter().any(|r| r.verification));
                    if comm.rank() == 0 {
                        println!(
                            "G6 {name} P={ranks} threshold={threshold_ry:e} Ry group={} nk={} E={:.17e} Ha full_delta={de:e} Ha legacy_delta={dl:e} Ha density_max={dn:e} e/Bohr^3",
                            report.group_order,
                            report.irreducible_k_count,
                            r.total_hartree()
                        );
                    }
                    if threshold_ry == 1e-12 {
                        let view = delivered::force_state_prepared(
                            StartInput {
                                cell: &cell,
                                atoms: &atoms,
                                potentials: &potentials,
                                options: &options,
                            },
                            &mut run.outcome,
                            &prepared,
                        )
                        .unwrap();
                        let forces = dft_engine::forces::forces(&view).unwrap();
                        let stress = dft_engine::stress::frozen_stress(&view).unwrap();
                        let df = forces
                            .values
                            .iter()
                            .flatten()
                            .zip(old_forces.iter().flatten())
                            .map(|(a, b)| (a - b).abs())
                            .fold(0., f64::max);
                        let ds = stress
                            .values
                            .iter()
                            .flatten()
                            .zip(old_stress.iter().flatten())
                            .map(|(a, b)| (a - b).abs())
                            .fold(0., f64::max);
                        assert!(
                            df < 1e-6 && ds < 1e-6,
                            "symmetrized derivatives dF={df:e} Ha/Bohr dS={ds:e} Ha/Bohr^3"
                        );
                        if comm.rank() == 0 {
                            println!(
                                "G6 derivatives {name} P={ranks} dF={df:e} Ha/Bohr dS={ds:e} Ha/Bohr^3"
                            );
                        }
                    }
                });
            }
        }
    }
}

/// cargo test --release -p dft-engine --test refscf g6_spin_fd -- --ignored --nocapture
#[test]
#[ignore = "spin/FD full SCF matrix; release command above"]
fn g6_spin_fd_symmetry_matches_full_and_legacy() {
    use dft_engine::{
        Functional, GammaReal, PreparedSymmetry, ScfOptions, SpeciesPotential, StartDensity,
        collinear::{CollinearProblem, SpinSector},
        finite_temperature::FiniteTProblem,
        refscf::{
            electronic::{self, Policy},
            symmetry::Prepared,
        },
        refstart::StartInput,
    };
    use dft_geometry::{Atom, Cell, kpoints::monkhorst_pack};
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0.2, 5., 0.], [0., 0.1, 5.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.25; 3])];
    let full = monkhorst_pack([2; 3], [0.5; 3]).unwrap();
    for spin in [false, true] {
        let old = PreparedSymmetry::new(&cell, &atoms, &full, 1e-7, !spin).unwrap();
        let prepared = Prepared::from_prepared(&cell, &atoms, &full, 1e-7, !spin, &old).unwrap();
        let pp = common::potential_with_atomic_start(
            1.,
            1.8,
            true,
            161,
            if spin { "SLA PW" } else { "SLA PZ" },
        );
        let potentials = [SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let options = ScfOptions {
            functional: if spin {
                Functional::LdaPw
            } else {
                Functional::LdaPz
            },
            cutoff_hartree: 2.5,
            bands: Some(if spin { 4 } else { 12 }),
            gamma_real: GammaReal::Off,
            start_density: StartDensity::Atomic,
            max_iterations: 300,
            grid_override: Some([16; 3]),
            ..Default::default()
        };
        let strict = ScfOptions {
            standard: None,
            verification: dft_engine::Verification::Off,
            eigen_tolerance: 1e-10,
            density_tolerance: 1e-10,
            energy_tolerance: 1e-11,
            ..options.clone()
        };
        let sector = SpinSector {
            n_up: 1,
            n_down: 0,
            axis_cartesian: [0., 0., 1.],
            initial_moment_seed: None,
        };
        let policy = if spin {
            Policy::FixedSpin(sector.clone())
        } else {
            Policy::FermiDirac {
                theta_ha: 0.01,
                band_cap: 12,
            }
        };
        let legacy_energy = if spin {
            let mut p = CollinearProblem::new(
                &cell,
                &atoms,
                &potentials,
                strict,
                prepared.points(),
                sector,
            )
            .unwrap();
            p.enable_prepared_symmetry(&old).unwrap();
            let r = p.run().unwrap();
            assert!(r.converged);
            r.energy.total_hartree()
        } else {
            let mut p = FiniteTProblem::new(
                &cell,
                &atoms,
                &potentials,
                strict,
                prepared.points(),
                0.01,
                12,
            )
            .unwrap();
            p.enable_prepared_symmetry(&old).unwrap();
            let r = p.run().unwrap();
            assert!(r.converged);
            r.thermodynamics.free_energy_ha
        };
        for threshold_ry in [1e-8, 1e-12] {
            let profile = Profile {
                threshold_ry,
                ..Default::default()
            };
            let full_run = electronic::run_with(
                StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &options,
                },
                profile,
                &full,
                policy.clone(),
                &Serial,
            )
            .unwrap();
            assert!(full_run.converged);
            for ranks in [1, 3, 8] {
                dft_mpi::simulate(ranks, |comm| {
                    let r = electronic::run_prepared_with(
                        StartInput {
                            cell: &cell,
                            atoms: &atoms,
                            potentials: &potentials,
                            options: &options,
                        },
                        profile,
                        &prepared,
                        policy.clone(),
                        comm,
                    )
                    .unwrap();
                    assert!(r.converged && r.verification.as_ref().unwrap().passed());
                    let delta = (r.thermodynamics.free_energy_ha - legacy_energy).abs();
                    let full_delta = (r.thermodynamics.free_energy_ha
                        - full_run.thermodynamics.free_energy_ha)
                        .abs();
                    let cap = if threshold_ry == 1e-8 { 5e-9 } else { 1e-10 };
                    assert!(
                        delta < cap && full_delta < cap,
                        "spin={spin} legacy={delta:e} full={full_delta:e} Ha"
                    );
                    let dn = r
                        .density
                        .iter()
                        .flatten()
                        .zip(full_run.density.iter().flatten())
                        .map(|(a, b)| (a - b).abs())
                        .fold(0., f64::max);
                    assert!(dn < if threshold_ry == 1e-8 { 1e-5 } else { 1e-7 });
                    let report = r.symmetry.unwrap();
                    assert_eq!(report.channels, if spin { 2 } else { 1 });
                    assert!(report.maps.iter().any(|r| r.verification));
                    if comm.rank() == 0 {
                        println!(
                            "G6 spin={spin} P={ranks} threshold={threshold_ry:e} Ry F={:.17e} Ha legacy_delta={delta:e} Ha full_delta={full_delta:e} Ha density_max={dn:e} e/Bohr^3",
                            r.thermodynamics.free_energy_ha
                        );
                    }
                });
            }
        }
    }
}

// G4C: one-electron synthetic open shell. This is not a material/QE reference.
fn g4c_potential(with_core: bool) -> dft_upf::Pseudopotential {
    let base = common::potential_with_atomic_start(1., 1., true, 401, "PBE");
    let r = base.mesh().r();
    let array = |v: &[f64]| {
        v.iter()
            .map(|v| format!("{v:.17e}"))
            .collect::<Vec<_>>()
            .join(" ")
    };
    let sample = |f: &dyn Fn(f64) -> f64| array(&r.iter().map(|&r| f(r)).collect::<Vec<_>>());
    let mut beta = String::new();
    for (j, l) in [0, 0, 1].into_iter().enumerate() {
        let i = j + 1;
        beta.push_str(&format!("<PP_BETA.{i} index=\"{i}\" angular_momentum=\"{l}\" cutoff_radius=\"12\" cutoff_radius_index=\"401\">{}</PP_BETA.{i}>", sample(&|r| r.powi(l+1)*(-r*r/i as f64).exp())));
    }
    let core = if with_core {
        format!(
            "<PP_NLCC>{}</PP_NLCC>",
            sample(&|r| (0.5 / PI).powf(1.5) * (-0.5 * r * r).exp())
        )
    } else {
        String::new()
    };
    dft_upf::parse_upf(&format!(r#"<UPF version="2.0.1"><PP_INFO>SYNTHETIC G4C Gaussian ion</PP_INFO>
<PP_HEADER generated="synthetic" author="tests" date="2026" comment="synthetic" element="Xx" pseudo_type="NC" relativistic="scalar" is_ultrasoft="F" is_paw="F" is_coulomb="F" has_so="F" has_wfc="F" has_gipaw="F" core_correction="{}" functional="PBE" z_valence="1" total_psenergy="0" wfc_cutoff="0" rho_cutoff="0" l_max="1" l_max_rho="0" l_local="-1" mesh_size="401" number_of_wfc="2" number_of_proj="3"/>
<PP_MESH mesh="401"><PP_R>{}</PP_R><PP_RAB>{}</PP_RAB></PP_MESH>
<PP_LOCAL>{}</PP_LOCAL><PP_NONLOCAL>{beta}<PP_DIJ>2 0.4 0 0.4 4 0 0 0 6</PP_DIJ></PP_NONLOCAL>
<PP_RHOATOM>{}</PP_RHOATOM>{core}<PP_PSWFC>
<PP_CHI.1 l="0" label="s" occupation="1">{}</PP_CHI.1><PP_CHI.2 l="1" label="p" occupation="0">{}</PP_CHI.2>
</PP_PSWFC></UPF>"#, if with_core {"T"} else {"F"}, array(r), array(base.mesh().rab()),
        array(&base.local_hartree().iter().map(|v|2.*v).collect::<Vec<_>>()),
        sample(&|r|4.*PI*r*r*(1./PI).powf(1.5)*(-r*r).exp()),
        sample(&|r|r*(-0.6*r*r).exp()), sample(&|r|r*r*(-0.5*r*r).exp()))).unwrap()
}

#[test]
#[ignore = "SCF matrix; cargo test --release -p dft-engine --test refscf g4c_spin_pbe_delivered -- --ignored --nocapture"]
fn g4c_spin_pbe_delivered() {
    g4c_pbe_matrix(&["fixed", "moment"]);
}

#[test]
#[ignore = "G4D unqualified: expected standard no-core free-moment accuracy failure; cargo test --release -p dft-engine --test refscf g4d_free_moment_pbe_acceptance -- --ignored --nocapture"]
fn g4d_free_moment_pbe_acceptance() {
    g4c_pbe_matrix(&["free"]);
}

fn g4c_pbe_matrix(modes: &[&str]) {
    use dft_engine::{
        Functional, GammaReal, ScfOptions, SpeciesPotential, StartDensity,
        collinear::{CollinearProblem, SpinSector},
        refscf::electronic::{self, Policy},
        refstart::StartInput,
    };
    use dft_geometry::{Atom, Cell, kpoints::KPoint};
    use dft_numerics::switches::RunConfig;
    use dft_physics::occupations::Constraint;
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0.2, 5., 0.], [0., 0.1, 5.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    let points = [KPoint {
        fractional: [0.; 3],
        weight: 1.,
    }];
    // Pin the qualified G4B combination to isolate PBE/NLCC coverage from routing defaults.
    let config = RunConfig::from_pairs(&[("refscf_spin_fd_moment", "on")]).unwrap();
    let _scope = config.enter();
    let mut failures = Vec::new();
    for core in [false, true] {
        let pp = g4c_potential(core);
        let potentials = [SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        for &mode in modes {
            let sector = SpinSector {
                n_up: 1,
                n_down: 0,
                axis_cartesian: [0., 0., 1.],
                initial_moment_seed: Some(0.3),
            };
            let options = ScfOptions {
                functional: Functional::Pbe,
                cutoff_hartree: 2.5,
                bands: Some(if mode == "fixed" { 4 } else { 12 }),
                gamma_real: GammaReal::Off,
                start_density: StartDensity::Atomic,
                max_iterations: 400,
                mixing: if mode == "free" {
                    dft_engine::Mixing::Linear { beta: 1. }
                } else {
                    dft_engine::Mixing::default()
                },
                residual_eigen_schedule: false,
                initial_eigen_tolerance: 1e-10,
                eigen_tolerance: 1e-10,
                density_tolerance: 1e-10,
                energy_tolerance: 1e-11,
                ..ScfOptions::strict_reference()
            };
            let constraint = if mode == "moment" {
                Constraint::FixedMoment(0.3)
            } else {
                Constraint::Total
            };
            let legacy = if mode == "fixed" {
                CollinearProblem::new(
                    &cell,
                    &atoms,
                    &potentials,
                    options.clone(),
                    &points,
                    sector.clone(),
                )
            } else {
                CollinearProblem::new_fermi_dirac(
                    &cell,
                    &atoms,
                    &potentials,
                    options.clone(),
                    &points,
                    sector.clone(),
                    0.01,
                    constraint,
                    12,
                )
            }
            .unwrap()
            .run()
            .unwrap();
            assert!(legacy.converged, "legacy {mode} core={core}");
            eprintln!(
                "G4C ORACLE mode={mode} core={core} M={} mu_B F={:?}",
                legacy.moment_mu_b, legacy.thermodynamics
            );
            let oracle = legacy
                .thermodynamics
                .map_or(legacy.energy.total_hartree(), |t| t.free_energy_ha);
            let policy = if mode == "fixed" {
                Policy::FixedSpin(sector)
            } else {
                Policy::SpinFermiDirac {
                    sector,
                    theta_ha: 0.01,
                    constraint,
                    band_cap: 12,
                }
            };
            let options = ScfOptions {
                standard: Some(dft_engine::StandardProfile { threshold_ha: 5e-9 }),
                verification: dft_engine::Verification::FixedPoint,
                ..options
            };
            for threshold in [1e-8, 1e-12] {
                for ranks in [1, 3, 8] {
                    let deltas = dft_mpi::simulate(ranks, |comm| {
                        let _scope = config.enter();
                        let diagnose = mode == "free" && !core && threshold == 1e-8 && ranks == 1;
                        let mut trace = electronic::Trace::default();
                        let r = electronic::run_with_trace(
                            StartInput {
                                cell: &cell,
                                atoms: &atoms,
                                potentials: &potentials,
                                options: &options,
                            },
                            Profile {
                                threshold_ry: threshold,
                                ..Profile::default()
                            },
                            &points,
                            policy.clone(),
                            comm,
                            diagnose.then_some(&mut trace),
                        )
                        .unwrap();
                        if diagnose {
                            let old_down: f64 =
                                legacy.finite_occupations.as_ref().unwrap().occupations()[12..]
                                    .iter()
                                    .sum();
                            eprintln!(
                                "G4D_ORACLE F={oracle:.17e} Ha M={:.17e} mu_B Ndown={old_down:.17e} e",
                                legacy.moment_mu_b
                            );
                            for t in trace.records() {
                                let down: f64 = t.occupations[12..].iter().sum();
                                eprintln!(
                                    "G4D_TRACE {},{},{},{:.10e},{:.10e},{:.10e},{:.10e},{:.10e},{:.10e},{:.10e}",
                                    t.iteration,
                                    t.verification,
                                    t.redo,
                                    t.estimator_terms_ry[0],
                                    t.estimator_terms_ry[1],
                                    t.estimator_terms_ry[2],
                                    t.stiffness_ha,
                                    t.free_energy_ha - oracle,
                                    t.moment_mu_b - legacy.moment_mu_b,
                                    down - old_down
                                );
                            }
                        }
                        assert!(
                            r.converged && r.verification.as_ref().unwrap().passed(),
                            "{mode} core={core}"
                        );
                        let delta = (r.thermodynamics.free_energy_ha - oracle).abs();
                        if comm.rank() == 0 {
                            eprintln!("G4C delivered M={} mu_B", r.moment_mu_b);
                            eprintln!(
                                "G4C {mode} core={core} P={ranks} threshold={threshold:e} Ry F={:.17e} Ha delta={delta:e} Ha",
                                r.thermodynamics.free_energy_ha
                            );
                        }
                        assert!((r.channel_electrons.iter().sum::<f64>() - 1.).abs() < 1e-10);
                        assert_eq!(r.energy.nlcc.is_some(), core);
                        delta
                    });
                    let delta = deltas.into_iter().fold(0., f64::max);
                    if delta >= if threshold == 1e-8 { 5e-9 } else { 1e-10 } {
                        failures.push(format!("{mode} core={core} P={ranks} threshold={threshold:e} Ry delta={delta:e} Ha"));
                    }
                }
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// Synthetic diamond: automatic-grid trajectory fingerprint for G6B.
#[test]
fn g6b_symmetry_diamond_diagnostics() {
    use dft_engine::{
        ScfOptions, SpeciesPotential,
        refscf::{k_driver, symmetry::Prepared},
        refstart::StartInput,
    };
    let cell = common::cell();
    let atoms = common::atoms();
    let full = dft_geometry::kpoints::monkhorst_pack([2; 3], [0.5; 3]).unwrap();
    let prepared = Prepared::new(&cell, &atoms, &full, 1e-7, true).unwrap();
    let pp = common::potential_with_atomic_start(2., 1.8, true, 81, "SLA PZ");
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: 1.,
        bands: Some(3),
        max_iterations: 3,
        verification: dft_engine::Verification::Off,
        ..Default::default()
    };
    let run = k_driver::run_prepared(
        StartInput {
            cell: &cell,
            atoms: &atoms,
            potentials: &potentials,
            options: &options,
        },
        Profile::default(),
        &prepared,
    )
    .unwrap();
    let density_hash = run
        .outcome
        .result
        .density
        .iter()
        .fold(0xcbf29ce484222325u64, |hash, value| {
            (hash ^ value.to_bits()).wrapping_mul(0x100000001b3)
        });
    println!(
        "G6B automatic diamond energy={:.17e} Ha density_fnv={density_hash:016x}",
        run.outcome.result.energy.total_hartree()
    );
    for r in &run.symmetry.unwrap().maps {
        assert!(r.pre_projection_relative_norm > 1e-6);
        assert!(r.post_projection_invariance_norm < 1e-13);
        println!(
            "G6B map={} pre_norm={:.17e} e/Bohr^(3/2) relative={:.17e} post={:.17e} e/Bohr^(3/2)",
            r.iteration, r.norm, r.pre_projection_relative_norm, r.post_projection_invariance_norm
        );
    }
}

/// Direct reference API must agree with the CLI's explicit-grid rule.
#[test]
fn g6b_symmetry_direct_explicit_grid() {
    use dft_engine::{
        ScfOptions, SpeciesPotential,
        refscf::{k_driver, symmetry::Prepared},
        refstart::StartInput,
    };
    let cell = common::cell();
    let atoms = common::atoms();
    let full = dft_geometry::kpoints::monkhorst_pack([2; 3], [0.5; 3]).unwrap();
    let prepared = Prepared::new(&cell, &atoms, &full, 1e-7, true).unwrap();
    let pp = common::potential_with_atomic_start(2., 1.8, true, 81, "SLA PZ");
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    for dims in [[28; 3], [25; 3]] {
        let options = ScfOptions {
            cutoff_hartree: 1.,
            bands: Some(3),
            max_iterations: 1,
            grid_override: Some(dims),
            verification: dft_engine::Verification::Off,
            ..Default::default()
        };
        dft_mpi::simulate(3, |comm| {
            let result = k_driver::run_prepared_with(
                StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &options,
                },
                Profile::default(),
                &prepared,
                comm,
            );
            if dims == [25; 3] {
                assert!(
                    matches!(result, Err(dft_engine::ScfError::InvalidInput(message)) if message.contains("incompatible explicit FFT grid"))
                );
            } else {
                let run = result.unwrap();
                assert_eq!(run.outcome.result.density.shape(), &dims);
                for r in &run.symmetry.unwrap().maps {
                    assert!(r.pre_projection_relative_norm > 1e-6);
                    assert!(r.post_projection_invariance_norm < 1e-13);
                }
            }
        });
    }
}

/// Synthetic isolated one-electron Gaussian with a Fourier-truncated NLCC core.
/// No external pseudopotential or QE data are used.
#[test]
fn g4e_vacuum() {
    use dft_engine::{
        Functional, GammaReal, ScfOptions, SpeciesPotential, StartDensity,
        collinear::{CollinearProblem, SpinSector},
        refscf::electronic::{self, Policy},
        refstart::StartInput,
    };
    use dft_geometry::{Atom, Cell, kpoints::KPoint};
    use dft_numerics::switches::RunConfig;
    let cell = Cell::from_user_rows_bohr([[10., 0., 0.], [0., 10., 0.], [0., 0., 10.]]).unwrap();
    let atoms = [Atom::new("Xx", [0., 0., 0.])];
    let points = [KPoint {
        fractional: [0.; 3],
        weight: 1.,
    }];
    let pp = g4c_potential(true);
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let sector = SpinSector {
        n_up: 1,
        n_down: 0,
        axis_cartesian: [0., 0., 1.],
        initial_moment_seed: Some(1.),
    };
    let options = ScfOptions {
        functional: Functional::Pbe,
        cutoff_hartree: 2.5,
        bands: Some(4),
        gamma_real: GammaReal::Off,
        start_density: StartDensity::Uniform,
        max_iterations: 300,
        residual_eigen_schedule: false,
        initial_eigen_tolerance: 1e-11,
        eigen_tolerance: 1e-11,
        density_tolerance: 1e-10,
        energy_tolerance: 1e-11,
        ..ScfOptions::strict_reference()
    };
    let cfg = |on| {
        RunConfig::from_pairs(&[("xc_screen_negative", if on { "on" } else { "off" })]).unwrap()
    };
    {
        let config = cfg(false);
        let _s = config.enter();
        let err = CollinearProblem::new(
            &cell,
            &atoms,
            &potentials,
            options.clone(),
            &points,
            sector.clone(),
        )
        .unwrap()
        .run()
        .err()
        .expect("negative NLCC must fail")
        .to_string();
        assert!(
            err.contains("channel") && err.contains("minimum") && err.contains("grid point"),
            "{err}"
        );
        let reference_options = ScfOptions {
            standard: Some(dft_engine::StandardProfile { threshold_ha: 5e-9 }),
            verification: dft_engine::Verification::FixedPoint,
            ..options.clone()
        };
        dft_mpi::simulate(3, |comm| {
            let _s = config.enter();
            let err = electronic::run_with(
                StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &reference_options,
                },
                Profile::default(),
                &points,
                Policy::FixedSpin(sector.clone()),
                comm,
            )
            .err()
            .expect("negative NLCC must fail")
            .to_string();
            assert!(
                err.contains("channel") && err.contains("minimum") && err.contains("grid point"),
                "{err}"
            );
        });
    }
    let config = cfg(true);
    let _s = config.enter();
    let legacy = CollinearProblem::new(
        &cell,
        &atoms,
        &potentials,
        options.clone(),
        &points,
        sector.clone(),
    )
    .unwrap()
    .run()
    .unwrap();
    assert!(legacy.converged);
    let oracle = legacy.energy.total_hartree();
    assert!(
        dft_physics::xc::screening_records()
            .iter()
            .any(|r| r.points > 0 && r.minimum_bohr3 < 0.)
    );
    let options = ScfOptions {
        standard: Some(dft_engine::StandardProfile { threshold_ha: 5e-9 }),
        verification: dft_engine::Verification::FixedPoint,
        ..options
    };
    for threshold in [1e-8, 1e-12] {
        for ranks in [1, 3, 8] {
            dft_mpi::simulate(ranks, |comm| {
                let _s = config.enter();
                dft_physics::xc::clear_screening_records();
                let r = electronic::run_with(
                    StartInput {
                        cell: &cell,
                        atoms: &atoms,
                        potentials: &potentials,
                        options: &options,
                    },
                    Profile {
                        threshold_ry: threshold,
                        ..Profile::default()
                    },
                    &points,
                    Policy::FixedSpin(sector.clone()),
                    comm,
                )
                .unwrap();
                assert!(r.converged && r.verification.as_ref().unwrap().passed());
                let delta = (r.energy.total_hartree() - oracle).abs();
                if comm.rank() == 0 {
                    eprintln!(
                        "G4E P={ranks} threshold={threshold:e} Ry E={:.17e} Ha legacy_delta={delta:e} Ha",
                        r.energy.total_hartree()
                    );
                }
                assert!(
                    delta < if threshold == 1e-8 { 5e-9 } else { 1e-10 },
                    "{delta:e}"
                );
                let records = dft_physics::xc::screening_records();
                assert!(records.iter().any(|r| r.points > 0 && r.minimum_bohr3 < 0.));
            });
        }
    }
}

#[test]
fn g4e_positive_density_driver_bits() {
    use dft_engine::{
        Functional, GammaReal, ScfOptions, ScfProblem, SpeciesPotential, StartDensity,
        refstart::StartInput,
    };
    use dft_geometry::{Atom, Cell};
    use dft_numerics::switches::RunConfig;
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0.1, 5., 0.], [0., 0.2, 5.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    for (functional, label) in [
        (Functional::LdaPz, "SLA PZ"),
        (Functional::LdaPw, "SLA PW"),
        (Functional::Pbe, "PBE"),
    ] {
        let pp = common::potential_with_atomic_start(2., 1.8, false, 161, label);
        let potentials = [SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let options = ScfOptions {
            functional,
            cutoff_hartree: 1.5,
            bands: Some(3),
            gamma_real: GammaReal::On,
            start_density: StartDensity::Uniform,
            ..Default::default()
        };
        let mut old = None;
        for setting in ["off", "on"] {
            let config = RunConfig::from_pairs(&[("xc_screen_negative", setting)]).unwrap();
            let _s = config.enter();
            let legacy = ScfProblem::new(&cell, &atoms, &potentials, options.clone())
                .unwrap()
                .run()
                .unwrap();
            let reference = dft_engine::refscf::driver::run_with(
                StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &options,
                },
                Profile::default(),
                &Serial,
            )
            .unwrap();
            assert!(reference.outcome.converged);
            let bits = (
                legacy.energy.total_hartree().to_bits(),
                legacy
                    .density
                    .iter()
                    .map(|n| n.to_bits())
                    .collect::<Vec<_>>(),
                reference.outcome.result.energy.total_hartree().to_bits(),
                reference
                    .outcome
                    .result
                    .density
                    .iter()
                    .map(|n| n.to_bits())
                    .collect::<Vec<_>>(),
            );
            if let Some(old) = &old {
                assert_eq!(&bits, old, "{functional:?}");
            } else {
                old = Some(bits);
            }
            if setting == "on" {
                assert!(
                    dft_physics::xc::screening_records()
                        .iter()
                        .all(|r| r.points == 0)
                );
            }
        }
    }
}

// Synthetic isolated open shell with Fourier-truncated NLCC; no external UPF.
// Release command: cargo test --release -p dft-engine --test refscf g4f_vacuum_fd -- --include-ignored --nocapture
#[test]
#[ignore = "isolated FD atom matrix in release"]
fn g4f_vacuum_fd_fixed() {
    g4f_vacuum_fd_case(false);
}

#[test]
#[ignore = "isolated FD atom matrix in release"]
fn g4f_vacuum_fd_free() {
    g4f_vacuum_fd_case(true);
}

fn g4f_vacuum_fd_case(free: bool) {
    use dft_engine::{
        Functional, GammaReal, ScfOptions, SpeciesPotential, StartDensity,
        collinear::{CollinearProblem, SpinSector},
        refscf::electronic::{self, Policy},
        refstart::StartInput,
    };
    use dft_geometry::{Atom, Cell, kpoints::KPoint};
    use dft_numerics::switches::RunConfig;
    use dft_physics::occupations::Constraint;
    let size = 12.;
    let cell = Cell::from_user_rows_bohr([[size, 0., 0.], [0., size, 0.], [0., 0., size]]).unwrap();
    let atoms = [Atom::new("Xx", [0., 0., 0.])];
    let points = [KPoint {
        fractional: [0.; 3],
        weight: 1.,
    }];
    let pp = g4f_potential();
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let sector = SpinSector {
        n_up: 1,
        n_down: 0,
        axis_cartesian: [0., 0., 1.],
        initial_moment_seed: Some(1.0),
    };
    let options = ScfOptions {
        functional: if free {
            Functional::LdaPw
        } else {
            Functional::Pbe
        },
        allow_functional_override: true,
        cutoff_hartree: 2.6,
        bands: Some(12),
        gamma_real: GammaReal::Off,
        start_density: StartDensity::Uniform,
        max_iterations: 400,
        mixing: dft_engine::Mixing::Linear { beta: 1.0 },
        residual_eigen_schedule: false,
        initial_eigen_tolerance: 1e-11,
        eigen_tolerance: 1e-11,
        density_tolerance: 1e-10,
        energy_tolerance: 1e-11,
        ..ScfOptions::strict_reference()
    };
    let config = RunConfig::from_pairs(&[("xc_screen_negative", "on")]).unwrap();
    let _s = config.enter();
    let constraint = if free {
        Constraint::Total
    } else {
        Constraint::FixedMoment(0.9999999)
    };
    let legacy = CollinearProblem::new_fermi_dirac(
        &cell,
        &atoms,
        &potentials,
        options.clone(),
        &points,
        sector.clone(),
        0.005,
        constraint,
        12,
    )
    .unwrap()
    .run()
    .unwrap();
    assert!(legacy.converged);
    let oracle = legacy.thermodynamics.unwrap().free_energy_ha;
    let oracle_moment = legacy.moment_mu_b;
    let options = ScfOptions {
        standard: Some(dft_engine::StandardProfile { threshold_ha: 5e-9 }),
        verification: dft_engine::Verification::FixedPoint,
        ..options
    };
    for ranks in [1, 3, 8] {
        dft_mpi::simulate(ranks, |comm| {
            let _s = config.enter();
            dft_physics::xc::clear_screening_records();
            let r = electronic::run_with(
                StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &options,
                },
                Profile::default(),
                &points,
                Policy::SpinFermiDirac {
                    sector: sector.clone(),
                    theta_ha: 0.005,
                    constraint,
                    band_cap: 12,
                },
                comm,
            )
            .unwrap();
            assert!(r.converged && r.verification.as_ref().unwrap().passed());
            let delta = (r.thermodynamics.free_energy_ha - oracle).abs();
            assert!(delta < 5e-9, "{delta:e} Ha");
            assert!((r.moment_mu_b - oracle_moment).abs() < 1e-4);
            assert!((r.channel_electrons.iter().sum::<f64>() - 1.).abs() < 1e-10);
            assert!(
                dft_physics::xc::screening_records()
                    .iter()
                    .any(|r| r.points > 0 && r.minimum_bohr3 < 0.)
            );
            if comm.rank() == 0 {
                eprintln!(
                    "G4F free={free} P={ranks} F={:.17e} Ha legacy_delta={delta:e} Ha moment={} mu_B",
                    r.thermodynamics.free_energy_ha, r.moment_mu_b
                );
            }
        });
    }
}

// A short-range 10 Ha Gaussian well localizes the open shell while leaving
// the -1/r tail unchanged. Its vacuum slabs expose the erroneous local budget.
fn g4f_potential() -> dft_upf::Pseudopotential {
    let with_core = true;
    let base = common::potential_with_atomic_start(1., 1., true, 401, "PBE");
    let r = base.mesh().r();
    let array = |v: &[f64]| {
        v.iter()
            .map(|v| format!("{v:.17e}"))
            .collect::<Vec<_>>()
            .join(" ")
    };
    let sample = |f: &dyn Fn(f64) -> f64| array(&r.iter().map(|&r| f(r)).collect::<Vec<_>>());
    let mut beta = String::new();
    for (j, l) in [0, 0, 1].into_iter().enumerate() {
        let i = j + 1;
        beta.push_str(&format!("<PP_BETA.{i} index=\"{i}\" angular_momentum=\"{l}\" cutoff_radius=\"12\" cutoff_radius_index=\"401\">{}</PP_BETA.{i}>", sample(&|r| r.powi(l+1)*(-r*r/i as f64).exp())));
    }
    let core = if with_core {
        format!(
            "<PP_NLCC>{}</PP_NLCC>",
            sample(&|r| (0.5 / PI).powf(1.5) * (-0.5 * r * r).exp())
        )
    } else {
        String::new()
    };
    dft_upf::parse_upf(&format!(r#"<UPF version="2.0.1"><PP_INFO>SYNTHETIC G4F localized Gaussian ion</PP_INFO>
<PP_HEADER generated="synthetic" author="tests" date="2026" comment="synthetic" element="Xx" pseudo_type="NC" relativistic="scalar" is_ultrasoft="F" is_paw="F" is_coulomb="F" has_so="F" has_wfc="F" has_gipaw="F" core_correction="{}" functional="PBE" z_valence="1" total_psenergy="0" wfc_cutoff="0" rho_cutoff="0" l_max="1" l_max_rho="0" l_local="-1" mesh_size="401" number_of_wfc="2" number_of_proj="3"/>
<PP_MESH mesh="401"><PP_R>{}</PP_R><PP_RAB>{}</PP_RAB></PP_MESH>
<PP_LOCAL>{}</PP_LOCAL><PP_NONLOCAL>{beta}<PP_DIJ>2 0.4 0 0.4 4 0 0 0 6</PP_DIJ></PP_NONLOCAL>
<PP_RHOATOM>{}</PP_RHOATOM>{core}<PP_PSWFC>
<PP_CHI.1 l="0" label="s" occupation="1">{}</PP_CHI.1><PP_CHI.2 l="1" label="p" occupation="0">{}</PP_CHI.2>
</PP_PSWFC></UPF>"#, if with_core {"T"} else {"F"}, array(r), array(base.mesh().rab()),
        array(&base.local_hartree().iter().zip(r).map(|(v,r)|2.*(v - 10. * (-r*r).exp())).collect::<Vec<_>>()),
        sample(&|r|4.*PI*r*r*(1./PI).powf(1.5)*(-r*r).exp()),
        sample(&|r|r*(-0.6*r*r).exp()), sample(&|r|r*r*(-0.5*r*r).exp()))).unwrap()
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
        }
    }
    stop_to_verified_clocks(run);
}
