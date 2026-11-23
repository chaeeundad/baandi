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
