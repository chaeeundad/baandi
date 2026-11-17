#[allow(dead_code)]
mod common;
use common::potential_with_functional;
use dft_engine::{collinear::*, *};
use dft_geometry::{Atom, Cell, kpoints::KPoint};

fn fixture_with(threads: usize, counts: [usize; 2], seed: f64, finite: bool) -> CollinearProblem {
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0.5, 5., 0.], [0.3, 0.2, 5.]]).unwrap();
    let pp = potential_with_functional((counts[0] + counts[1]) as f64, 1., true, 401, "SLA PW");
    let options = ScfOptions {
        functional: Functional::LdaPw,
        cutoff_hartree: 2.5,
        bands: Some(4),
        threads,
        max_iterations: 200,
        gamma_real: GammaReal::Off,
        ..ScfOptions::default()
    };
    let sector = SpinSector {
        n_up: counts[0],
        n_down: counts[1],
        axis_cartesian: [0., 0., 1.],
        initial_moment_seed: Some(seed),
    };
    let atoms = [Atom::new("Xx", [0.; 3])];
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let points = [KPoint {
        fractional: [0.; 3],
        weight: 1.,
    }];
    if finite {
        return CollinearProblem::new_fermi_dirac(
            &cell,
            &atoms,
            &species,
            options,
            &points,
            sector,
            0.01,
            dft_physics::occupations::Constraint::Total,
            8,
        )
        .unwrap();
    }
    CollinearProblem::new(
        &cell,
        &[Atom::new("Xx", [0.; 3])],
        &[SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }],
        options,
        &points,
        sector,
    )
    .unwrap()
}

fn fixture(threads: usize) -> CollinearProblem {
    fixture_with(threads, [1, 1], 0., false)
}

// Recorded BEFORE engine edits on b1c95a61282e66e741e98d1ef0199545faf1cae2,
// pure Rust, default switches, synthetic skewed cell, standard+verification.
const BASELINE_ENERGY_BITS: u64 = 13829793526668695510;

#[test]
fn off_preserves_main_bits() {
    let mut p = fixture(1);
    p.force_spin_share(false);
    let r = p.run().unwrap();
    assert!(r.converged);
    assert_eq!(r.energy.total_hartree().to_bits(), BASELINE_ENERGY_BITS);
    assert_eq!(r.spin_share.shared_maps, 0);
    assert!(r.spin_share.fallback_counts[0] > 0);
}

fn equal_states(a: &CollinearResult, b: &CollinearResult) {
    assert!(a.converged && b.converged);
    assert_eq!(a.history.len(), b.history.len());
    assert!((a.energy.total_hartree() - b.energy.total_hartree()).abs() <= 1e-12);
    assert_eq!(
        a.energy.total_hartree().to_bits(),
        b.energy.total_hartree().to_bits()
    );
    assert_eq!(a.moment_mu_b.to_bits(), 0.0_f64.to_bits());
    assert_eq!(b.moment_mu_b.to_bits(), 0.0_f64.to_bits());
    assert_eq!(
        a.spin_share.davidson_iterations,
        b.spin_share.davidson_iterations
    );
    for (a, b) in a.k_points.iter().flatten().zip(b.k_points.iter().flatten()) {
        assert_eq!(a.solution.eigenvalues, b.solution.eigenvalues);
        assert_eq!(a.solution.eigenvectors, b.solution.eigenvectors);
        assert_eq!(a.solution.residuals, b.solution.residuals);
        assert_eq!(a.solution.iterations, b.solution.iterations);
    }
    let (
        dft_physics::spin::SpinDensity::Collinear { channels: a, .. },
        dft_physics::spin::SpinDensity::Collinear { channels: b, .. },
    ) = (&a.density, &b.density)
    else {
        panic!()
    };
    for (a, b) in a
        .iter()
        .flat_map(|x| x.iter())
        .zip(b.iter().flat_map(|x| x.iter()))
    {
        assert_eq!(a.to_bits(), b.to_bits());
    }
}

#[test]
fn shares_every_map_including_verification_and_halves_h_work() {
    let mut off = fixture(1);
    off.force_spin_share(false);
    let off = off.run().unwrap();
    let mut previous = None;
    for threads in [1, 4] {
        let mut on = fixture(threads);
        on.force_spin_share(true);
        let on = on.run().unwrap();
        equal_states(&off, &on);
        assert_eq!(on.spin_share.fallback_counts.iter().sum::<usize>(), 0);
        assert_eq!(on.spin_share.shared_maps, off.spin_share.fallback_counts[0]);
        assert_eq!(
            2 * on.spin_share.operator_applications,
            off.spin_share.operator_applications
        );
        let a = off.verification.as_ref().unwrap();
        let b = on.verification.as_ref().unwrap();
        assert_eq!(a.attempts.len(), b.attempts.len());
        assert!(!b.attempts.is_empty());
        for (a, b) in a.attempts.iter().zip(&b.attempts) {
            assert!(a.passed && b.passed);
            assert_eq!(a.davidson_iterations, b.davidson_iterations);
            assert_eq!(a.operator_applications, 2 * b.operator_applications);
            assert_eq!(a.coulomb_error_ha.to_bits(), b.coulomb_error_ha.to_bits());
        }
        assert_eq!(
            on.k_points[0][0].solution.eigenvalues,
            on.k_points[1][0].solution.eigenvalues
        );
        assert_eq!(
            on.k_points[0][0].solution.eigenvectors,
            on.k_points[1][0].solution.eigenvectors
        );
        println!(
            "T={threads}: E={:.17} Ha, SCF={}, logical Davidson={}, H columns {} -> {}, shared maps={}",
            on.energy.total_hartree(),
            on.history.len(),
            on.spin_share.davidson_iterations,
            off.spin_share.operator_applications,
            on.spin_share.operator_applications,
            on.spin_share.shared_maps
        );
        if let Some(previous) = &previous {
            equal_states(previous, &on);
        }
        previous = Some(on);
    }
}

#[test]
fn symmetry_breaking_and_finite_temperature_fall_back() {
    for (counts, seed, finite, reason) in [
        ([1, 1], 0.5, false, 4),
        ([1, 0], 0., false, 5),
        ([1, 1], 0., true, 3),
    ] {
        let mut p = fixture_with(1, counts, seed, finite);
        p.force_spin_share(true);
        let r = p.run().unwrap();
        assert!(r.converged);
        assert_eq!(r.spin_share.shared_maps, 0);
        assert!(r.spin_share.fallback_counts[reason] > 0);
        assert_eq!(
            r.spin_share.fallback_counts.iter().sum::<usize>(),
            r.spin_share.fallback_counts[reason]
        );
        assert!(
            r.k_points
                .iter()
                .flatten()
                .all(|k| k.solution.operator_applications > 0)
        );
        if finite {
            assert!(
                r.finite_occupations
                    .as_ref()
                    .unwrap()
                    .occupations()
                    .iter()
                    .any(|&f| f > 0. && f < 1.)
            );
        }
    }
}

#[test]
fn environment_switch_is_used_without_override() {
    // Process isolation exercises the real cached switch without unsafe set_var
    // or interference between concurrent tests.
    const CHILD: &str = "SPIN_SHARE_TEST_CHILD";
    if let Ok(mode) = std::env::var(CHILD) {
        let r = fixture(1).run().unwrap();
        assert_eq!(r.energy.total_hartree().to_bits(), BASELINE_ENERGY_BITS);
        assert_eq!(r.spin_share.shared_maps > 0, mode == "on");
        return;
    }
    for mode in ["on", "off"] {
        let result = std::process::Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "environment_switch_is_used_without_override",
                "--nocapture",
            ])
            .env(CHILD, mode)
            .env("BAANDI_SPIN_SHARE", mode)
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stdout)
        );
    }
}

#[test]
fn atomic_start_sharing_preserves_lda_and_pbe_states() {
    let cell = Cell::from_user_rows_bohr([[7., 0., 0.], [0.6, 6.8, 0.], [0.4, 0.3, 7.2]]).unwrap();
    let atoms = [Atom::new("Xx", [0.; 3]), Atom::new("Xx", [0.2, 0., 0.])];
    for (functional, label) in [(Functional::LdaPw, "SLA PW"), (Functional::Pbe, "PBE")] {
        let pp = common::potential_with_orbitals(1., 1., true, 401, label);
        let run = |share| {
            let mut p = CollinearProblem::new(
                &cell,
                &atoms,
                &[SpeciesPotential {
                    label: "Xx",
                    potential: &pp,
                }],
                ScfOptions {
                    functional,
                    cutoff_hartree: 3.,
                    gamma_real: GammaReal::Off,
                    start_wavefunctions: StartWavefunctions::Atomic,
                    ..ScfOptions::default()
                },
                &[KPoint {
                    fractional: [0.; 3],
                    weight: 1.,
                }],
                SpinSector {
                    n_up: 1,
                    n_down: 1,
                    axis_cartesian: [0., 0., 1.],
                    initial_moment_seed: Some(0.),
                },
            )
            .unwrap();
            p.force_spin_share(share);
            p.run().unwrap()
        };
        let off = run(false);
        let on = run(true);
        equal_states(&off, &on);
        assert!(on.atomic_start.initial_h_applications > 0);
        assert_eq!(
            off.atomic_start.initial_h_applications,
            2 * on.atomic_start.initial_h_applications
        );
        assert_eq!(
            off.spin_share.operator_applications,
            2 * on.spin_share.operator_applications
        );
        assert_eq!(on.spin_share.fallback_counts.iter().sum::<usize>(), 0);
        println!(
            "atomic {functional:?}: dE=0 Ha, H columns {} -> {}, atomic H {} -> {}",
            off.spin_share.operator_applications,
            on.spin_share.operator_applications,
            off.atomic_start.initial_h_applications,
            on.atomic_start.initial_h_applications
        );
    }
}
