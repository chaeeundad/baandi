#[allow(dead_code)]
mod common;
use dft_engine::{GammaReal, ScfOptions, ScfProblem, ScfResult, SpeciesPotential};
use dft_geometry::{Atom, Cell};

fn model(kind: usize) -> (Cell, Vec<Atom>, dft_upf::Pseudopotential) {
    let rows = match kind {
        0 => [[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]],
        1 => [[5., 0., 0.], [0.5, 4., 0.], [0.2, 0.4, 6.]],
        2 => [[4., 0., 0.], [0., 4., 0.], [0., 0., 8.]],
        _ => [[5., 0., 0.], [0., 5., 0.], [0., 0., 10.]],
    };
    let atoms = if kind < 2 {
        vec![Atom::new("Xx", [0.13, 0.21, 0.32])]
    } else {
        vec![
            Atom::new("Xx", [0.13, 0.21, 0.17]),
            Atom::new("Xx", [0.6, 0.7, 0.68]),
        ]
    };
    (
        Cell::from_user_rows_bohr(rows).unwrap(),
        atoms,
        common::potential(2., 1., kind != 2, 401),
    )
}
fn options(mode: GammaReal, threads: usize, seed: u64) -> ScfOptions {
    ScfOptions {
        gamma_real: mode,
        threads,
        seed,
        cutoff_hartree: 2.5,
        bands: Some(5),
        ..Default::default()
    }
}
fn problem(kind: usize, opts: ScfOptions) -> ScfProblem {
    let (cell, atoms, pp) = model(kind);
    ScfProblem::new(
        &cell,
        &atoms,
        &[SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }],
        opts,
    )
    .unwrap()
}
fn run(kind: usize, mode: GammaReal, threads: usize, seed: u64) -> ScfResult {
    problem(kind, options(mode, threads, seed)).run().unwrap()
}
fn differences(a: &ScfResult, b: &ScfResult, volume: f64) -> (f64, f64, f64) {
    let occupied = a.occupations.iter().filter(|&&x| x > 0.).count();
    (
        (a.total_hartree() - b.total_hartree()).abs(),
        a.eigenvalues[..occupied]
            .iter()
            .zip(&b.eigenvalues)
            .map(|(x, y)| (x - y).abs())
            .fold(0., f64::max),
        (a.density
            .iter()
            .zip(&b.density)
            .map(|(x, y)| (x - y).powi(2))
            .sum::<f64>()
            * volume
            / a.density.len() as f64)
            .sqrt(),
    )
}
#[test]
#[ignore = "release calibration before real comparison; see real-gamma-r5-results.md"]
fn calibrate_complex_noise() {
    for kind in 0..4 {
        let mut opts = options(GammaReal::Off, 1, 2026);
        opts.free_electron_initial_guess = false;
        let a = problem(kind, opts).run().unwrap();
        for (threads, seed) in [(1, 2027), (8, 2026), (8, 2027)] {
            let mut opts = options(GammaReal::Off, threads, seed);
            opts.free_electron_initial_guess = false;
            let b = problem(kind, opts).run().unwrap();
            let d = differences(&a, &b, model(kind).0.volume_bohr3());
            println!(
                "CALIBRATION kind={kind} threads={threads} seed={seed} E={:.17e} eigen={:.17e} density={:.17e} SCF={}/{}",
                d.0,
                d.1,
                d.2,
                a.history.len(),
                b.history.len()
            );
        }
    }
}
#[test]
fn off_is_default() {
    assert_eq!(ScfOptions::default().gamma_real, GammaReal::Off);
}

#[test]
#[ignore = "release SCF matrix; cargo test --release -p dft-engine --test gamma_real_scf -- --include-ignored --nocapture"]
fn real_complex_and_thread_matrix() {
    for kind in 0..4 {
        let complex = run(kind, GammaReal::Off, 1, 2026);
        let real = run(kind, GammaReal::On, 1, 2026);
        assert!(real.gamma_real_used);
        assert!(real.gamma_real_fallback_reason.is_none());
        let d = differences(&complex, &real, model(kind).0.volume_bohr3());
        let occupied = real.occupations.iter().filter(|&&x| x > 0.).count();
        println!(
            "REAL kind={kind} E={:.17e} deltaE={:.17e} eigen={:.17e} density={:.17e} SCF={}/{} gap={:.9e} steps={}/{} H={}/{}",
            real.total_hartree(),
            d.0,
            d.1,
            d.2,
            complex.history.len(),
            real.history.len(),
            real.eigenvalues[occupied] - real.eigenvalues[occupied - 1],
            complex
                .history
                .iter()
                .map(|h| h.davidson_iterations)
                .sum::<usize>(),
            real.history
                .iter()
                .map(|h| h.davidson_iterations)
                .sum::<usize>(),
            complex
                .history
                .iter()
                .map(|h| h.operator_applications)
                .sum::<usize>(),
            real.history
                .iter()
                .map(|h| h.operator_applications)
                .sum::<usize>()
        );
        // Fixed design energy gate and complex-noise calibration, frozen in the results note.
        assert!(d.0 < 1e-10 && d.1 < 3e-8 && d.2 < 7e-7, "{d:?}");
        for threads in [2, 4, 8, 12] {
            let other = run(kind, GammaReal::On, threads, 2026);
            // Scheduling must change no arithmetic, including completed orbital coefficients.
            assert_eq!(
                real.total_hartree().to_bits(),
                other.total_hartree().to_bits()
            );
            assert!(
                real.density
                    .iter()
                    .zip(&other.density)
                    .all(|(a, b)| a.to_bits() == b.to_bits())
            );
            assert_eq!(real.k_points[0].solution, other.k_points[0].solution);
            assert_eq!(real.history.len(), other.history.len());
        }
    }
}

#[test]
fn fallback_modes_and_warm_vectors() {
    use dft_engine::{Locking, OccupationPolicy, SpectralGuard, StandardProfile};
    use dft_geometry::kpoints::KPoint;
    for case in 0..6 {
        let (cell, atoms, pp) = model(0);
        let mut opts = options(GammaReal::Auto, 1, 2026);
        let mut points = vec![KPoint {
            fractional: [0.; 3],
            weight: 1.,
        }];
        match case {
            0 => points[0].fractional[0] = 0.1,
            1 => {
                points = vec![
                    KPoint {
                        fractional: [0.25; 3],
                        weight: 0.5,
                    },
                    KPoint {
                        fractional: [-0.25; 3],
                        weight: 0.5,
                    },
                ]
            }
            2 => opts.standard = Some(StandardProfile::default()),
            3 => opts.locking = Locking::Auto,
            4 => opts.spectral_guard = SpectralGuard::Final,
            _ => opts.eigensolver = dft_solvers::SolverChoice::Hybrid(Default::default()),
        }
        let create = |o| {
            ScfProblem::with_kpoints(
                &cell,
                &atoms,
                &[SpeciesPotential {
                    label: "Xx",
                    potential: &pp,
                }],
                o,
                &points,
                OccupationPolicy::Insulator {
                    minimum_gap_hartree: 1e-7,
                },
            )
        };
        let mut auto = create(opts.clone()).unwrap();
        let n = auto.uniform_density();
        let a = auto.step_kpoints(&n, None, 1e-8).unwrap();
        assert!(
            auto.gamma_real_diagnostics()
                .all(|d| !d.gamma_real_used && d.fallback_reason.is_some())
        );
        opts.gamma_real = GammaReal::Off;
        let b = create(opts.clone())
            .unwrap()
            .step_kpoints(&n, None, 1e-8)
            .unwrap();
        // Fallback must preserve the entire complex solver result, not just energy tolerance.
        assert_eq!(a.density, b.density);
        assert_eq!(
            a.energy.total_hartree().to_bits(),
            b.energy.total_hartree().to_bits()
        );
        for (a, b) in a.k_points.iter().zip(&b.k_points) {
            assert_eq!(a.solution, b.solution);
        }
        opts.gamma_real = GammaReal::On;
        let error = match create(opts) {
            Err(e) => e.to_string(),
            Ok(mut p) => p.step_kpoints(&n, None, 1e-8).err().unwrap().to_string(),
        };
        assert!(error.contains("gamma_real=on"), "{error}");
    }
    let mut p = problem(0, options(GammaReal::Off, 1, 2026));
    let n = p.uniform_density();
    let cold = p.step_kpoints(&n, None, 1e-8).unwrap();
    let mut bad = cold.solution.eigenvectors;
    bad[(0, 0)].im = 0.1;
    let guesses = [bad];
    let mut auto = problem(0, options(GammaReal::Auto, 1, 2026));
    let a = auto.step_kpoints(&n, Some(&guesses), 1e-8).unwrap();
    assert!(
        auto.gamma_real_diagnostics()
            .next()
            .unwrap()
            .fallback_reason
            .as_ref()
            .unwrap()
            .contains("warm start")
    );
    let b = p.step_kpoints(&n, Some(&guesses), 1e-8).unwrap();
    assert_eq!(a.solution, b.solution);
    assert!(
        problem(0, options(GammaReal::On, 1, 2026))
            .step_kpoints(&n, Some(&guesses), 1e-8)
            .is_err()
    );
}

#[test]
fn spin_fallback_is_explicit_and_identical() {
    use dft_engine::collinear::{CollinearProblem, SpinSector};
    use dft_geometry::kpoints::KPoint;
    let (cell, atoms, _) = model(0);
    let pp = common::potential_with_functional(2., 1., true, 401, "SLA PW");
    let create = |mode| {
        let mut opts = options(mode, 1, 2026);
        opts.functional = dft_engine::Functional::LdaPw;
        CollinearProblem::new(
            &cell,
            &atoms,
            &[SpeciesPotential {
                label: "Xx",
                potential: &pp,
            }],
            opts,
            &[KPoint {
                fractional: [0.; 3],
                weight: 1.,
            }],
            SpinSector {
                n_up: 2,
                n_down: 0,
                axis_cartesian: [0., 0., 1.],
                initial_moment_seed: None,
            },
        )
    };
    assert!(
        create(GammaReal::On)
            .err()
            .unwrap()
            .to_string()
            .contains("spin")
    );
    let auto = create(GammaReal::Auto).unwrap().run().unwrap();
    let off = create(GammaReal::Off).unwrap().run().unwrap();
    assert!(!auto.gamma_real_used);
    assert!(auto.gamma_real_fallback_reason.unwrap().contains("spin"));
    // No spin numerical code is routed through the real bridge in R5.
    assert_eq!(
        auto.energy.total_hartree().to_bits(),
        off.energy.total_hartree().to_bits()
    );
    assert_eq!(format!("{:?}", auto.density), format!("{:?}", off.density));
    for (a, b) in auto
        .k_points
        .iter()
        .flatten()
        .zip(off.k_points.iter().flatten())
    {
        assert_eq!(a.solution, b.solution);
    }
}

#[test]
#[ignore = "release converged derivative comparison"]
fn forces_and_stress_through_completion() {
    use dft_engine::forces::{ElectronicState, ForceState, forces};
    use dft_engine::stress::frozen_stress;
    for kind in [0, 1, 2, 3] {
        let (cell, atoms, pp) = model(kind);
        let species = [SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let mut derivatives = Vec::new();
        for mode in [GammaReal::Off, GammaReal::On] {
            let mut opts = options(mode, 1, 2026);
            opts.density_tolerance = 1e-10;
            opts.energy_tolerance = 1e-12;
            opts.eigen_tolerance = 1e-11;
            let r = ScfProblem::new(&cell, &atoms, &species, opts.clone())
                .unwrap()
                .run()
                .unwrap();
            let state = ForceState {
                cell: &cell,
                atoms: &atoms,
                potentials: &species,
                options: &opts,
                electronic: ElectronicState::Fixed(&r),
                density_floor: None,
                symmetry: &[],
            };
            derivatives.push((
                forces(&state).unwrap().values,
                frozen_stress(&state).unwrap().values,
            ));
        }
        let a = &derivatives[0];
        let b = &derivatives[1];
        let force =
            a.0.iter()
                .flatten()
                .zip(b.0.iter().flatten())
                .map(|(x, y)| (x - y).abs())
                .fold(0., f64::max);
        let stress =
            a.1.iter()
                .flatten()
                .zip(b.1.iter().flatten())
                .map(|(x, y)| (x - y).abs())
                .fold(0., f64::max);
        println!(
            "DERIVATIVES kind={kind} force={force:.17e} Ha/Bohr stress={stress:.17e} Ha/Bohr^3"
        );
        // Stricter than existing synthetic FD caps (1e-6 force, 1e-7 stress),
        // matching the 1e-9 state-identity comparisons in the derivative suite.
        assert!(force < 1e-9 && stress < 1e-9);
    }
}

#[allow(dead_code)]
#[path = "common/checkpoint.rs"]
mod checkpoint_fixture;
#[test]
fn restart_path_identity_and_warm_interoperability() {
    use dft_engine::{CompletionMode, ExecutionConfig, ScfSession, checkpoint::*};
    let (model, exec) = checkpoint_fixture::identity(0, 1);
    for source in [GammaReal::Off, GammaReal::On] {
        let target = if source == GammaReal::On {
            GammaReal::Off
        } else {
            GammaReal::Auto
        };
        let mut s = ScfSession::new(
            problem(0, options(source, 1, 2026)),
            CompletionMode::Strict,
            ExecutionConfig::default(),
        )
        .unwrap();
        s.initialize().unwrap();
        for _ in 0..3 {
            s.step().unwrap();
        }
        let exact = s
            .checkpoint(
                CheckpointKind::ExactContinuation,
                model.clone(),
                exec.clone(),
            )
            .unwrap();
        let error =
            ScfSession::restore_exact(problem(0, options(target, 1, 2026)), &exact, &model, &exec)
                .err()
                .unwrap();
        assert_eq!(error, Rejection::Options);
        let same =
            ScfSession::restore_exact(problem(0, options(source, 1, 2026)), &exact, &model, &exec)
                .unwrap();
        let warm = s
            .checkpoint(CheckpointKind::WarmStart, model.clone(), exec.clone())
            .unwrap();
        let resumed = ScfSession::warm_start(
            problem(0, options(target, 1, 2026)),
            &warm,
            &model,
            &exec,
            CompletionMode::Strict,
            ExecutionConfig {
                calc_id: "new-warm-run".into(),
                ..Default::default()
            },
        )
        .unwrap();
        let direct = checkpoint_fixture::finish(s);
        let restored = checkpoint_fixture::finish(same);
        // Exact restart must retain all numerical bits on the same path.
        checkpoint_fixture::assert_identical(&direct, &restored);
        let resumed = checkpoint_fixture::finish(resumed);
        let d = differences(&direct, &resumed, 64.);
        // Warm restarts use the frozen SCF comparison gates, not exact identity.
        assert!(d.0 < 1e-10 && d.1 < 3e-8 && d.2 < 7e-7, "{d:?}");
        println!(
            "RESTART source={source:?} target={target:?} delta={d:?} real_used={}",
            resumed.gamma_real_used
        );
    }
}

#[test]
fn finite_temperature_is_excluded_without_mislabeling_spin() {
    use dft_engine::finite_temperature::FiniteTProblem;
    use dft_geometry::kpoints::KPoint;
    let (cell, atoms, pp) = model(0);
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let create = |mode| {
        FiniteTProblem::new(
            &cell,
            &atoms,
            &species,
            options(mode, 1, 2026),
            &[KPoint {
                fractional: [0.; 3],
                weight: 1.,
            }],
            0.001,
            5,
        )
    };
    let error = create(GammaReal::On).err().unwrap().to_string();
    assert!(error.contains("Fermi-Dirac"));
    let auto = create(GammaReal::Auto).unwrap().run().unwrap();
    let off = create(GammaReal::Off).unwrap().run().unwrap();
    assert!(!auto.gamma_real_used);
    assert_eq!(
        auto.gamma_real_fallback_reason.as_deref(),
        Some("Fermi-Dirac occupations")
    );
    assert_eq!(auto.state.density, off.state.density);
    assert_eq!(auto.state.solution, off.state.solution);
}
