use super::*;
use crate::scf::pool_fixture as common;
fn run(theta: f64, bands: usize, electrons: f64) -> FiniteTResult {
    let cell =
        Cell::from_user_rows_bohr([[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]]).unwrap();
    let pp = common::potential(electrons, 1.0, false, 401);
    FiniteTProblem::new(
        &cell,
        &[Atom::new("Xx", [0.0; 3])],
        &[SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }],
        ScfOptions {
            cutoff_hartree: if electrons == 2.0 { 4.5 } else { 2.5 },
            bands: Some(bands),
            density_tolerance: 1e-10,
            energy_tolerance: 1e-11,
            eigen_tolerance: 1e-10,
            max_iterations: 200,
            ..Default::default()
        },
        &[KPoint {
            fractional: [0.0; 3],
            weight: 1.0,
        }],
        theta,
        bands,
    )
    .unwrap()
    .run()
    .unwrap()
}
#[test]
fn synthetic_metal() {
    let r = run(0.01, 12, 1.0);
    assert!(r.converged);
    let e = r.occupations.evaluation.unwrap();
    assert!(e.electron_residual.abs() < 1e-12);
    assert!((r.state.energy.total_hartree() - r.state.energy.band_form).abs() < 1e-10);
    println!(
        "synthetic N=1 E={:.16e} Ha F={:.16e} Ha mu={:.16e} Ha entropy={:.16e} iterations={}",
        r.thermodynamics.internal_energy_ha,
        r.thermodynamics.free_energy_ha,
        e.roots[0].unwrap().mu_ha,
        e.entropy_over_kb,
        r.iterations
    );
}

#[test]
fn empty_band_refinement() {
    let a = run(0.01, 12, 1.0);
    for bands in [16, 20] {
        let b = run(0.01, bands, 1.0);
        assert!(b.converged);
        let df = (a.thermodynamics.free_energy_ha - b.thermodynamics.free_energy_ha).abs();
        let dn = density_distance(&a.state.density, &b.state.density, 125.0, 1.0);
        println!(
            "band refinement 12 -> {bands}: delta F={df:.16e} Ha, density={dn:.16e} Bohr^-3/2"
        );
        assert!(df <= 1e-8 && dn <= 1e-8);
    }
}
#[test]
fn stationary_entropy_derivative() {
    let theta = 0.05;
    let a = run(theta, 24, 1.0);
    assert!(a.converged);
    // Independent binary entropy, from the emitted occupations (Gamma, g=2).
    let entropy = a
        .occupations
        .occupations()
        .iter()
        .map(|f| {
            let p = f / 2.0;
            if p == 0.0 || p == 1.0 {
                0.0
            } else {
                -2.0 * (p * p.ln() + (1.0 - p) * (1.0 - p).ln())
            }
        })
        .sum::<f64>();
    assert!((entropy - a.occupations.evaluation.unwrap().entropy_over_kb).abs() < 1e-10);
    let mut errors = Vec::new();
    for h in [0.004, 0.002, 0.001] {
        let plus = run(theta + h, 24, 1.0);
        let minus = run(theta - h, 24, 1.0);
        assert!(plus.converged && minus.converged);
        let err = ((plus.thermodynamics.free_energy_ha - minus.thermodynamics.free_energy_ha)
            / (2.0 * h)
            + entropy)
            .abs();
        println!("entropy derivative h={h} Ha: error={err:.16e} dimensionless");
        errors.push(err);
    }
    assert!(
        errors[1] < 0.4 * errors[0] && errors[2] < 0.4 * errors[1],
        "{errors:?}"
    );
}
#[test]
fn spin_finite_roots_and_scf() {
    let cell =
        Cell::from_user_rows_bohr([[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]]).unwrap();
    let pp = common::potential_with_functional(1.0, 1.0, false, 401, "SLA PW");
    for constraint in [Constraint::Total, Constraint::FixedMoment(0.3)] {
        let r = super::super::collinear::CollinearProblem::new_fermi_dirac(
            &cell,
            &[Atom::new("Xx", [0.0; 3])],
            &[SpeciesPotential {
                label: "Xx",
                potential: &pp,
            }],
            ScfOptions {
                functional: Functional::LdaPw,
                cutoff_hartree: 2.5,
                eigen_tolerance: 1e-10,
                density_tolerance: 1e-9,
                max_iterations: 200,
                ..Default::default()
            },
            &[KPoint {
                fractional: [0.0; 3],
                weight: 1.0,
            }],
            super::super::collinear::SpinSector {
                n_up: 1,
                n_down: 0,
                axis_cartesian: [0.0, 0.0, 1.0],
                initial_moment_seed: None,
            },
            0.01,
            constraint,
            12,
        )
        .unwrap()
        .run()
        .unwrap();
        assert!(r.converged);
        assert!(r.finite_certificate.as_ref().unwrap().passed);
        assert!((r.energy.total_hartree() - r.energy.band_form).abs() < 1e-9);
        let e = r.finite_occupations.as_ref().unwrap().evaluation.unwrap();
        assert!(
            e.electron_residual.abs() < 1e-12 && e.moment_residual.unwrap_or(0.0).abs() < 1e-12
        );
        if let Constraint::FixedMoment(m) = constraint {
            assert!((r.moment_mu_b - m).abs() < 1e-10);
        }
        println!(
            "spin {constraint:?}: F={:.16e} Ha M={:.16e} iterations={}",
            r.thermodynamics.unwrap().free_energy_ha,
            r.moment_mu_b,
            r.history.len()
        );
    }
}
#[test]
#[ignore = "external Al NC data; cargo test --release -p dft-engine real_al_finite -- --ignored --nocapture"]
fn real_al_finite() {
    let dir = std::env::var("DFT_NC_DIR").expect("set DFT_NC_DIR to qualified PseudoDojo NC files");
    let path = std::path::Path::new(&dir).join("Al.upf");
    let pp = dft_upf::parse_upf(&std::fs::read_to_string(path).unwrap()).unwrap();
    let cell =
        Cell::from_user_rows_bohr([[0.0, 3.8, 3.8], [3.8, 0.0, 3.8], [3.8, 3.8, 0.0]]).unwrap();
    let points = dft_geometry::kpoints::monkhorst_pack([2; 3], [0.5; 3]).unwrap();
    let r = FiniteTProblem::new(
        &cell,
        &[Atom::new("Al", [0.0; 3])],
        &[SpeciesPotential {
            label: "Al",
            potential: &pp,
        }],
        ScfOptions {
            functional: Functional::Pbe,
            cutoff_hartree: 20.0,
            max_iterations: 200,
            ..Default::default()
        },
        &points,
        0.01,
        20,
    )
    .unwrap()
    .run()
    .unwrap();
    assert!(r.converged);
    println!(
        "Al smoke theta=0.01 Ha E={:.16e} Ha mu={:.16e} Ha iterations={}",
        r.thermodynamics.internal_energy_ha,
        r.occupations.evaluation.unwrap().roots[0].unwrap().mu_ha,
        r.iterations
    );
}

#[test]
fn gapped_zero_temperature_limit() {
    let cell =
        Cell::from_user_rows_bohr([[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]]).unwrap();
    let pp = common::potential(2.0, 1.0, false, 401);
    let fixed = ScfProblem::new(
        &cell,
        &[Atom::new("Xx", [0.0; 3])],
        &[SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }],
        ScfOptions {
            cutoff_hartree: 4.5,
            bands: Some(40),
            density_tolerance: 1e-10,
            eigen_tolerance: 1e-10,
            energy_tolerance: 1e-11,
            max_iterations: 200,
            ..Default::default()
        },
    )
    .unwrap()
    .run()
    .unwrap();
    let mut last = f64::INFINITY;
    for theta in [0.1, 0.01, 0.001] {
        let r = run(theta, 40, 2.0);
        assert!(r.converged);
        let delta = (r.thermodynamics.internal_energy_ha - fixed.energy.total_hartree()).abs();
        let dn = density_distance(&r.state.density, &fixed.density, 125.0, 2.0);
        println!("T->0 theta={theta} Ha delta E={delta:.16e} Ha density={dn:.16e} Bohr^-3/2");
        assert!(delta <= last + 1e-12);
        if theta <= 0.01 {
            assert!(delta <= 1e-8 && dn <= 1e-8);
        }
        last = delta;
    }
}

#[test]
fn insufficient_empty_band_cap_fails_closed() {
    let cell =
        Cell::from_user_rows_bohr([[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]]).unwrap();
    let pp = common::potential(1.0, 1.0, false, 401);
    let result = FiniteTProblem::new(
        &cell,
        &[Atom::new("Xx", [0.0; 3])],
        &[SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }],
        ScfOptions {
            cutoff_hartree: 2.5,
            bands: Some(1),
            ..Default::default()
        },
        &[KPoint {
            fractional: [0.0; 3],
            weight: 1.0,
        }],
        0.01,
        1,
    )
    .unwrap()
    .run();
    assert!(
        matches!(result,Err(ScfError::InvalidInput(message)) if message.contains("insufficient empty-state coverage"))
    );
}
