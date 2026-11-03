mod common;
use common::*;
use dft_engine::*;
use dft_geometry::Atom;
use dft_upf::Pseudopotential;
fn problem(pp: &Pseudopotential, options: ScfOptions) -> ScfProblem {
    ScfProblem::new(
        &cell(),
        &atoms(),
        &[SpeciesPotential {
            label: "Xx",
            potential: pp,
        }],
        options,
    )
    .unwrap()
}
fn check(result: &ScfResult, options: &ScfOptions) {
    let last = result.history.last().unwrap();
    assert!(last.density_residual <= options.density_tolerance);
    assert!(last.energy_change.unwrap().abs() <= options.energy_tolerance);
    assert!(last.max_eigen_residual <= options.eigen_tolerance);
    assert!(last.requested_eigen_tolerance <= options.eigen_tolerance);
    for r in &result.history {
        assert!(r.charge_error < 1e-10 && r.input_charge_error < 1e-10);
    }
    assert!(last.min_density >= -1e-14); // Bohr^-3, existing P3A roundoff rule.
    let max_charge = result
        .history
        .iter()
        .map(|r| r.charge_error.max(r.input_charge_error))
        .fold(0.0, f64::max);
    println!(
        "maximum iteration charge error={max_charge:.9e} relative; Pulay feasibility resets={}",
        result.history.iter().filter(|r| r.mixer_reset).count()
    );
    let ew = &result.ewald_breakdown;
    assert_eq!(
        ew.energy.0,
        ew.real.0 + ew.reciprocal.0 + ew.self_energy.0 + ew.background.0
    );
    println!(
        "Ewald audit: real={:.12e}, reciprocal={:.12e}, self={:.12e}, background={:.12e}, tail bound={:.3e} Ha",
        ew.real.0, ew.reciprocal.0, ew.self_energy.0, ew.background.0, ew.tail_bound.0
    );
    let e = &result.energy;
    assert_eq!(
        e.total_hartree(),
        e.kinetic + e.local + e.nonlocal + e.hartree + e.xc + e.ewald
    );
    assert_eq!(e.total_rydberg(), 2.0 * e.total_hartree());
    // Roundoff allowance covers independently accumulated scalar sums only.
    assert!((e.total_hartree() - e.harris_foulkes).abs() <= e.functional_difference_bound + 1e-12);
    assert!((e.total_hartree() - e.band_form).abs() <= 8.0 * last.max_eigen_residual + 1e-12);
    println!(
        "SCF: iterations={}, E={:.16e} Ha, E={:.16e} Ry, residual={:.9e}, eigen residual={:.9e} Ha, charge error={:.3e}, min n={:.9e} Bohr^-3",
        result.history.len(),
        result.total_hartree(),
        result.total_rydberg(),
        last.density_residual,
        last.max_eigen_residual,
        last.charge_error,
        last.min_density
    );
    println!(
        "audit: T={:.16e}, loc={:.16e}, loc0={:.16e}, NL={:.16e}, H={:.16e}, xc={:.16e}, Ewald={:.16e}, band={:.16e}, DC={:.16e}, band form={:.16e}, HF={:.16e}, bound={:.9e} Ha",
        e.kinetic,
        e.local,
        e.local_zero_mode,
        e.nonlocal,
        e.hartree,
        e.xc,
        e.ewald,
        e.band_sum,
        e.double_counting,
        e.band_form,
        e.harris_foulkes,
        e.functional_difference_bound
    );
}
#[test]
fn invalid_electrons_bands_and_options() {
    let pp = potential(3.0, 1.0, false, 401);
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    assert!(matches!(
        ScfProblem::new(
            &cell(),
            &[Atom::new("Xx", [0.0; 3])],
            &species,
            ScfOptions::default()
        ),
        Err(ScfError::OddElectronCount(3.0))
    ));
    let options = ScfOptions {
        bands: Some(2),
        ..Default::default()
    };
    assert!(matches!(
        ScfProblem::new(&cell(), &atoms(), &species, options),
        Err(ScfError::InsufficientBands {
            bands: 2,
            occupied: 3
        })
    ));
    assert!(
        ScfProblem::new(
            &cell(),
            &atoms(),
            &species,
            ScfOptions {
                density_tolerance: f64::NAN,
                ..Default::default()
            }
        )
        .is_err()
    );
}
#[test]
#[ignore = "multi-run SCF audit; run release acceptance command in ADR 0009"]
fn synthetic_scf_audit() {
    let pp = potential(4.0, 1.0, true, 1201);
    let options = ScfOptions {
        max_iterations: 600,
        // Preserve the seeded trajectory that exercises the positivity safeguard.
        residual_eigen_schedule: false,
        free_electron_initial_guess: false,
        ..Default::default()
    };
    let linear_options = ScfOptions {
        mixing: Mixing::Linear { beta: 0.05 },
        ..options.clone()
    };
    let linear = problem(&pp, linear_options.clone()).run().unwrap();
    check(&linear, &linear_options);
    let pulay = problem(&pp, options.clone()).run().unwrap();
    check(&pulay, &options);
    println!(
        "mixers: linear={}, Pulay={}, energy difference={:.9e} Ha",
        linear.history.len(),
        pulay.history.len(),
        (linear.total_hartree() - pulay.total_hartree()).abs()
    );
    assert!((linear.total_hartree() - pulay.total_hartree()).abs() < 1e-8);
    assert!(pulay.history.len() < linear.history.len());
    let repeated = problem(&pp, options.clone()).run().unwrap();
    assert_eq!(
        repeated.total_hartree().to_bits(),
        pulay.total_hartree().to_bits()
    );
    assert_eq!(repeated.density, pulay.density);
    assert_eq!(repeated.orbitals, pulay.orbitals);
    assert_eq!(repeated.eigenvalues, pulay.eigenvalues);
    assert_eq!(repeated.history.len(), pulay.history.len());
    for (a, b) in repeated.history.iter().zip(&pulay.history) {
        assert_eq!(a.density_residual.to_bits(), b.density_residual.to_bits());
        assert_eq!(
            a.energy.total_hartree().to_bits(),
            b.energy.total_hartree().to_bits()
        );
    }
    let other = problem(
        &pp,
        ScfOptions {
            seed: 997,
            ..options.clone()
        },
    )
    .run()
    .unwrap();
    check(&other, &options);
    assert!(other.history.iter().any(|r| r.mixer_reset));
    println!(
        "seed energy difference={:.9e} Ha",
        (other.total_hartree() - pulay.total_hartree()).abs()
    );
    assert!((other.total_hartree() - pulay.total_hartree()).abs() < 1e-8);
    let mut external = problem(&pp, options.clone());
    let step = external
        .step(&pulay.density, Some(&pulay.orbitals), 1e-10)
        .unwrap();
    let ev = step.solution.eigenvalues[..4]
        .iter()
        .zip(&pulay.eigenvalues)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    let dn = density_distance(&step.density, &pulay.density, cell().volume_bohr3(), 8.0);
    println!(
        "external fixed point: occupied eigenvalue change={ev:.9e} Ha, density change={dn:.9e}"
    );
    assert!(ev < 1e-7 && dn <= options.density_tolerance);
    // Controlled displacements about converged density: fixed orbitals are not
    // used; each perturbed potential is independently solved to 1e-11 Ha.
    let uniform = external.uniform_density();
    let mut last: Option<(f64, f64)> = None;
    for amplitude in [0.04, 0.02, 0.01, 0.005] {
        let n = &pulay.density + &((&uniform - &pulay.density) * amplitude);
        let s = external.step(&n, Some(&pulay.orbitals), 1e-11).unwrap();
        let residual = density_distance(&s.density, &n, cell().volume_bohr3(), 8.0);
        let delta = (s.energy.total_hartree() - s.energy.harris_foulkes).abs();
        println!(
            "functional trend: amplitude={amplitude}, residual={residual:.9e}, |direct-HF|={delta:.9e} Ha"
        );
        if let Some((r, e)) = last {
            let power = (delta / e).ln() / (residual / r).ln();
            println!("functional trend exponent={power:.6}");
            assert!((1.7..2.3).contains(&power));
        }
        last = Some((residual, delta));
    }
}
#[test]
#[ignore = "SCF invariances and cutoff sweep; run release acceptance command in ADR 0009"]
fn invariances_and_cutoff_trend() {
    let pp = potential(4.0, 1.0, true, 1201);
    let options = ScfOptions {
        grid_override: Some([32; 3]),
        ..Default::default()
    };
    let base = problem(&pp, options.clone()).run().unwrap();
    check(&base, &options);
    let coarse_problem = problem(&pp, ScfOptions::default());
    let coarse_dims = coarse_problem.grid_dimensions();
    let coarse = coarse_problem.run().unwrap();
    println!(
        "XC grid quadrature: default {coarse_dims:?} -> [32,32,32] energy difference={:.12e} Ha",
        (coarse.total_hartree() - base.total_hartree()).abs()
    );
    let mut shifted = atoms();
    for a in &mut shifted {
        for (r, d) in a.fractional_position.iter_mut().zip([0.137, -0.213, 0.319]) {
            *r += d;
        }
    }
    let translated = ScfProblem::new(
        &cell(),
        &shifted,
        &[SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }],
        options.clone(),
    )
    .unwrap()
    .run()
    .unwrap();
    let larger = problem(
        &pp,
        ScfOptions {
            grid_override: Some([36; 3]),
            ..options.clone()
        },
    )
    .run()
    .unwrap();
    println!(
        "invariances: translation={:.9e} Ha, grid 32->36={:.9e} Ha",
        (translated.total_hartree() - base.total_hartree()).abs(),
        (larger.total_hartree() - base.total_hartree()).abs()
    );
    assert!((translated.total_hartree() - base.total_hartree()).abs() < 1e-8);
    assert!((larger.total_hartree() - base.total_hartree()).abs() < 1e-8);
    for cutoff in [4.0, 5.0, 6.0] {
        let result = problem(
            &pp,
            ScfOptions {
                cutoff_hartree: cutoff,
                ..options.clone()
            },
        )
        .run()
        .unwrap();
        println!(
            "cutoff={} Ry: E={:.16e} Ha, gap={:.12e} Ha",
            2.0 * cutoff,
            result.total_hartree(),
            result.eigenvalues[4] - result.eigenvalues[3]
        );
    }
}
#[test]
#[ignore = "SCF error paths; run release acceptance command in ADR 0009"]
fn divergence_and_iteration_budget() {
    let pp = potential(4.0, 1.0, true, 801);
    let err = problem(
        &pp,
        ScfOptions {
            mixing: Mixing::Linear { beta: 100.0 },
            ..Default::default()
        },
    )
    .run()
    .unwrap_err();
    match err {
        ScfError::Diverged { reason, history } => println!(
            "oversized mixing: divergence after {} records: {reason}",
            history.len()
        ),
        other => panic!("wrong error: {other}"),
    }
    let err = problem(
        &pp,
        ScfOptions {
            max_iterations: 1,
            ..Default::default()
        },
    )
    .run()
    .unwrap_err();
    assert!(matches!(err,ScfError::NotConverged {history} if history.len()==1));
}

#[test]
fn small_synthetic_fixed_point() {
    let pp = potential(2.0, 1.0, false, 401);
    let c = dft_geometry::Cell::from_user_rows_bohr([
        [4.0, 0.0, 0.0],
        [0.0, 4.0, 0.0],
        [0.0, 0.0, 4.0],
    ])
    .unwrap();
    let options = ScfOptions {
        cutoff_hartree: 1.5,
        ..Default::default()
    };
    let result = ScfProblem::new(
        &c,
        &[Atom::new("Xx", [0.13, 0.21, 0.32])],
        &[SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }],
        options.clone(),
    )
    .unwrap()
    .run()
    .unwrap();
    check(&result, &options);
    let cert = result
        .history
        .last()
        .unwrap()
        .fixed_point_certificate
        .as_ref()
        .unwrap();
    assert!(cert.density_residual <= options.density_tolerance);
    assert!(cert.occupied_eigenvalue_change <= 1e-7);
}

#[test]
fn functional_labels_preserve_energy_and_reject_other_functionals() {
    let c = dft_geometry::Cell::from_user_rows_bohr([
        [4.0, 0.0, 0.0],
        [0.0, 4.0, 0.0],
        [0.0, 0.0, 4.0],
    ])
    .unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    let options = ScfOptions {
        cutoff_hartree: 1.5,
        ..Default::default()
    };
    let mut energies = Vec::new();
    for label in ["SLA PZ", " SLA  PZ   NOGX NOGC "] {
        let pp = potential_with_functional(2.0, 1.0, false, 401, label);
        let result = ScfProblem::new(
            &c,
            &atoms,
            &[SpeciesPotential {
                label: "Xx",
                potential: &pp,
            }],
            options.clone(),
        )
        .unwrap()
        .run()
        .unwrap();
        energies.push(result.total_hartree());
    }
    assert_eq!(energies[0].to_bits(), energies[1].to_bits()); // Exact equality, Ha.
    println!(
        "functional labels: energies={energies:?} Ha; difference={} Ha",
        energies[1] - energies[0]
    );
    let pp = potential_with_functional(2.0, 1.0, false, 401, "SLA PW");
    let err = ScfProblem::new(
        &c,
        &atoms,
        &[SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }],
        options,
    )
    .err()
    .expect("unsupported functional must fail");
    assert!(matches!(err, ScfError::UnsupportedFunctional(_)));
    assert_eq!(
        err.to_string(),
        "unsupported UPF functional \"SLA PW\": requires SLA PZ or SLA PZ NOGX NOGC"
    );
}

#[test]
fn residual_schedule_matches_legacy_energy_and_finishes_tightly() {
    let pp = potential(2.0, 1.0, false, 401);
    let c = dft_geometry::Cell::from_user_rows_bohr([
        [4.0, 0.0, 0.0],
        [0.0, 4.0, 0.0],
        [0.0, 0.0, 4.0],
    ])
    .unwrap();
    let mut energies = Vec::new();
    for residual_eigen_schedule in [false, true] {
        let options = ScfOptions {
            cutoff_hartree: 1.5,
            residual_eigen_schedule,
            ..Default::default()
        };
        let result = ScfProblem::new(
            &c,
            &[Atom::new("Xx", [0.13, 0.21, 0.32])],
            &[SpeciesPotential {
                label: "Xx",
                potential: &pp,
            }],
            options.clone(),
        )
        .unwrap()
        .run()
        .unwrap();
        check(&result, &options);
        assert!(
            result.history.last().unwrap().requested_eigen_tolerance <= options.eigen_tolerance
        );
        assert!(
            result
                .history
                .last()
                .unwrap()
                .fixed_point_certificate
                .is_some()
        );
        assert_eq!(
            result.history[0].requested_eigen_tolerance,
            if residual_eigen_schedule {
                1e-2
            } else {
                options.initial_eigen_tolerance
            }
        );
        energies.push(result.total_hartree());
    }
    println!(
        "schedule energy difference = {:.16e} Ha",
        (energies[0] - energies[1]).abs()
    );
    assert!((energies[0] - energies[1]).abs() < 1e-12);
}

#[test]
#[ignore = "full synthetic schedule comparison; cargo test --release -p dft-engine --test scf full_schedule_equivalence -- --ignored --nocapture"]
fn full_schedule_equivalence() {
    let pp = potential(4.0, 1.0, true, 1201);
    let legacy_options = ScfOptions {
        residual_eigen_schedule: false,
        free_electron_initial_guess: false,
        ..Default::default()
    };
    let legacy = problem(&pp, legacy_options.clone()).run().unwrap();
    check(&legacy, &legacy_options);
    let options = ScfOptions::default();
    let result = problem(&pp, options.clone()).run().unwrap();
    check(&result, &options);
    let difference = (result.total_hartree() - legacy.total_hartree()).abs();
    println!("full schedule energy difference = {difference:.16e} Ha");
    assert!(difference <= 1e-12);
}
