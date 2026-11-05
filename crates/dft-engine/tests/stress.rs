#![allow(clippy::needless_range_loop)]
#[allow(dead_code)]
mod common;
#[path = "../../dft-physics/tests/stress_support/mod.rs"]
mod stress_support;
use dft_engine::stress::{energy_at, frozen_stress};
use dft_engine::{ElectronicState, ForceState, ScfOptions, ScfProblem, SpeciesPotential};
use dft_geometry::{Atom, Cell};
fn strained(cell: &Cell, a: usize, b: usize, h: f64) -> Cell {
    let c = cell.lattice_columns_bohr();
    let mut d = c;
    for j in 0..3 {
        d[a][j] += h * c[b][j];
        if a != b {
            d[b][j] += h * c[a][j];
        }
    }
    Cell::from_user_rows_bohr(std::array::from_fn(|i| std::array::from_fn(|j| d[j][i]))).unwrap()
}
fn oracle(state: &ForceState<'_>) {
    let stress = frozen_stress(state).unwrap();
    println!(
        "analytic {:?} pressure={} antisymmetric={}",
        stress.values, stress.pressure, stress.antisymmetric_norm
    );
    let dims = match &state.electronic {
        ElectronicState::Fixed(r) => r.density.shape(),
        ElectronicState::FiniteTemperature(r) => r.state.density.shape(),
        ElectronicState::Collinear(r) => match &r.density {
            dft_physics::spin::SpinDensity::Collinear { channels, .. } => channels[0].shape(),
            _ => unreachable!(),
        },
    };
    let mut options = state.options.clone();
    options.strain_reference = Some(Box::new(state.cell.clone()));
    options.grid_override = Some(dims.try_into().unwrap());
    for a in 0..3 {
        for b in a..3 {
            let mut errors = [0.; 3];
            let mut values = [0.; 3];
            for (i, h) in [0.002, 0.001, 0.0005].into_iter().enumerate() {
                values[i] = (energy_at(state, &strained(state.cell, a, b, h), options.clone())
                    .unwrap()
                    - energy_at(state, &strained(state.cell, a, b, -h), options.clone()).unwrap())
                    / (2. * h * state.cell.volume_bohr3() * if a == b { 1. } else { 2. });
                errors[i] = (values[i] - stress.values[a][b]).abs();
            }
            println!("{a}{b} FD={values:?} errors={errors:?} Ha/Bohr^3");
            assert!(errors[2] < 1e-7, "{errors:?}");
            assert!(
                errors[2] < 0.15 * errors[0] || errors.iter().all(|e| *e < 1e-9),
                "no O(h^2) trend: {errors:?}"
            );
        }
    }
}
#[test]
#[ignore = "release strain SCF sweep: cargo test --release -p dft-engine --test stress -- --ignored --nocapture"]
fn synthetic_local_nonlocal() {
    for nonlocal in [false, true] {
        let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [0., 8., 0.], [0., 0., 8.]]).unwrap();
        let pp = stress_support::potential(1., 1.5, nonlocal, 401);
        let species = [SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let atoms = vec![
            Atom::new("Xx", [0.; 3]),
            Atom::new("Xx", [0.27, 0.23, 0.19]),
        ];
        let options = ScfOptions {
            cutoff_hartree: 1.,
            density_tolerance: 1e-11,
            energy_tolerance: 1e-13,
            eigen_tolerance: 1e-12,
            max_iterations: 500,
            grid_override: Some([24; 3]),
            ..Default::default()
        };
        let r = ScfProblem::new(&cell, &atoms, &species, options.clone())
            .unwrap()
            .run()
            .unwrap();
        println!("nonlocal={nonlocal}");
        oracle(&ForceState {
            cell: &cell,
            atoms: &atoms,
            potentials: &species,
            options: &options,
            electronic: ElectronicState::Fixed(&r),
            density_floor: None,
            symmetry: &[],
        });
    }
}

use dft_physics::data::erf;
use std::f64::consts::PI;
fn synthetic_text(core: bool, functional: &str) -> String {
    let n = 1001;
    let h = 0.012;
    let a = |f: &dyn Fn(f64) -> f64| {
        (0..n)
            .map(|i| format!("{:.17e}", f(i as f64 * h)))
            .collect::<Vec<_>>()
            .join(" ")
    };
    let nlcc = if core {
        format!(
            "<PP_NLCC>{}</PP_NLCC>",
            a(&|r| (0.5 / PI).powf(1.5) * (-0.5 * r * r).exp())
        )
    } else {
        String::new()
    };
    format!(
        r#"<UPF version="2.0.1"><PP_INFO>Synthetic Gaussian core</PP_INFO>
<PP_HEADER generated="synthetic" author="test" date="2026" comment="synthetic" element="Xx" pseudo_type="NC" relativistic="scalar" is_ultrasoft="F" is_paw="F" is_coulomb="F" has_so="F" has_wfc="F" has_gipaw="F" core_correction="{}" functional="{functional}" z_valence="2" total_psenergy="0" wfc_cutoff="0" rho_cutoff="0" l_max="0" l_max_rho="0" l_local="-1" mesh_size="{n}" number_of_wfc="0" number_of_proj="0"/>
<PP_MESH dx="{h}" mesh="{n}" xmin="0" rmax="12" zmesh="1"><PP_R>{}</PP_R><PP_RAB>{}</PP_RAB></PP_MESH>{nlcc}<PP_LOCAL>{}</PP_LOCAL><PP_NONLOCAL/><PP_RHOATOM>{}</PP_RHOATOM></UPF>"#,
        if core { "T" } else { "F" },
        a(&|r| r),
        a(&|_| h),
        a(&|r| if r == 0.0 {
            -4.0 / PI.sqrt()
        } else {
            -4.0 * erf(r / 2.0) / r
        }),
        a(&|r| 8.0 * PI * r * r * (0.5 / PI).powf(1.5) * (-0.5 * r * r).exp())
    )
}
fn synthetic(core: bool, functional: &str) -> dft_upf::Pseudopotential {
    dft_upf::parse_upf(&synthetic_text(core, functional)).unwrap()
}

#[test]
#[ignore = "release strain SCF sweep"]
fn synthetic_models_and_shapes() {
    use dft_engine::{
        Functional,
        collinear::{CollinearProblem, SpinSector},
        finite_temperature::FiniteTProblem,
    };
    use dft_geometry::kpoints::KPoint;
    for mode in ["skew", "slab", "nlcc", "pbe", "spin", "spin-pbe", "fd"] {
        if std::env::var("DFT_STRESS_CASE").is_ok_and(|selected| selected != mode) {
            continue;
        }
        let cell = Cell::from_user_rows_bohr(if mode == "slab" {
            [[7., 0., 0.], [0.5, 8., 0.], [0.2, 0.3, 16.]]
        } else {
            [[8., 0., 0.], [1., 7., 0.], [0.5, 1., 9.]]
        })
        .unwrap();
        let core = matches!(mode, "nlcc" | "pbe" | "spin-pbe");
        let functional = if mode.contains("pbe") {
            Functional::Pbe
        } else {
            Functional::LdaPw
        };
        let pp = if core {
            synthetic(
                true,
                if functional == Functional::Pbe {
                    "PBE"
                } else {
                    "SLA PW"
                },
            )
        } else {
            common::potential_with_functional(1., 1.5, true, 401, "SLA PW")
        };
        let species = [SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let atoms = vec![
            Atom::new("Xx", [0.; 3]),
            Atom::new("Xx", [0.27, 0.23, 0.19]),
        ];
        let options = ScfOptions {
            functional,
            cutoff_hartree: if core || mode == "fd" { 2. } else { 1. },
            davidson_subspace: if mode == "fd" { 96 } else { 48 },
            density_tolerance: 1e-11,
            energy_tolerance: 1e-13,
            eigen_tolerance: 1e-12,
            max_iterations: 500,
            grid_override: Some(if mode == "slab" {
                [24, 24, 48]
            } else if core {
                [48; 3]
            } else {
                [24; 3]
            }),
            ..Default::default()
        };
        let points = [KPoint {
            fractional: [0.; 3],
            weight: 1.,
        }];
        let evaluate = |electronic| {
            println!("case={mode}");
            oracle(&ForceState {
                cell: &cell,
                atoms: &atoms,
                potentials: &species,
                options: &options,
                electronic,
                density_floor: None,
                symmetry: &[],
            });
        };
        if mode.starts_with("spin") {
            let sector = SpinSector {
                n_up: if core { 3 } else { 1 },
                n_down: 1,
                axis_cartesian: [0., 0., 1.],
                initial_moment_seed: None,
            };
            let r =
                CollinearProblem::new(&cell, &atoms, &species, options.clone(), &points, sector)
                    .unwrap()
                    .run()
                    .unwrap();
            evaluate(ElectronicState::Collinear(&r));
            if !core {
                let scalar = ScfProblem::new(&cell, &atoms, &species, options.clone())
                    .unwrap()
                    .run()
                    .unwrap();
                let a = frozen_stress(&ForceState {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &species,
                    options: &options,
                    electronic: ElectronicState::Collinear(&r),
                    density_floor: None,
                    symmetry: &[],
                })
                .unwrap();
                let b = frozen_stress(&ForceState {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &species,
                    options: &options,
                    electronic: ElectronicState::Fixed(&scalar),
                    density_floor: None,
                    symmetry: &[],
                })
                .unwrap();
                for (a, b) in a.values.iter().flatten().zip(b.values.iter().flatten()) {
                    assert!((a - b).abs() < 1e-9);
                }
            }
        } else if mode == "fd" {
            let r =
                FiniteTProblem::new(&cell, &atoms, &species, options.clone(), &points, 0.05, 48)
                    .unwrap()
                    .run()
                    .unwrap();
            assert!(r.thermodynamics.entropy_term_ha > 1e-6);
            println!(
                "FD theta=0.05 Ha TS={} Ha F={} Ha E={} Ha",
                r.thermodynamics.entropy_term_ha,
                r.thermodynamics.free_energy_ha,
                r.thermodynamics.internal_energy_ha
            );
            evaluate(ElectronicState::FiniteTemperature(&r));
        } else {
            let r = ScfProblem::new(&cell, &atoms, &species, options.clone())
                .unwrap()
                .run()
                .unwrap();
            evaluate(ElectronicState::Fixed(&r));
        }
    }
}

#[test]
#[ignore = "external DFT_NC_DIR Si/MgO strain and convergence scans; release only"]
fn real_files() {
    let dir = std::path::PathBuf::from(std::env::var_os("DFT_NC_DIR").expect("set DFT_NC_DIR"));
    for (material, a, cutoff) in [("Si", 10.2612, 20.), ("Si", 10.6, 20.), ("MgO", 7.96, 40.)] {
        if std::env::var("DFT_STRESS_MATERIAL").is_ok_and(|selected| selected != material) {
            continue;
        }
        let labels = if material == "Si" {
            ["Si", "Si"]
        } else {
            ["Mg", "O"]
        };
        let p0 = dft_upf::from_file(dir.join(format!("{}.upf", labels[0]))).unwrap();
        let p1 = dft_upf::from_file(dir.join(format!("{}.upf", labels[1]))).unwrap();
        let mut species = vec![
            SpeciesPotential {
                label: labels[0],
                potential: &p0,
            },
            SpeciesPotential {
                label: labels[1],
                potential: &p1,
            },
        ];
        if material == "Si" {
            species.truncate(1);
        }
        let cell = Cell::from_user_rows_bohr([
            [0., a / 2., a / 2.],
            [a / 2., 0., a / 2.],
            [a / 2., a / 2., 0.],
        ])
        .unwrap();
        let atoms = vec![
            Atom::new(labels[0], [0.; 3]),
            Atom::new(labels[1], [if material == "Si" { 0.25 } else { 0.5 }; 3]),
        ];
        let options = ScfOptions {
            functional: dft_engine::Functional::Pbe,
            cutoff_hartree: cutoff,
            density_tolerance: 1e-9,
            energy_tolerance: 1e-11,
            eigen_tolerance: 1e-10,
            max_iterations: 300,
            ..Default::default()
        };
        let run = |options: ScfOptions| {
            let mut p = ScfProblem::new(&cell, &atoms, &species, options).unwrap();
            p.set_density_floor(1e-10).unwrap();
            p.run()
        };
        let r = match run(options.clone()) {
            Ok(r) => r,
            Err(e) => {
                println!("INCONCLUSIVE {material} a={a}: {e}");
                continue;
            }
        };
        let state = ForceState {
            cell: &cell,
            atoms: &atoms,
            potentials: &species,
            options: &options,
            electronic: ElectronicState::Fixed(&r),
            density_floor: Some(1e-10),
            symmetry: &[],
        };
        let tensor = frozen_stress(&state).unwrap();
        println!(
            "REAL {material} a={a} cutoff={cutoff} tensor={:?} pressure={} Ha/Bohr^3",
            tensor.values, tensor.pressure
        );
        let mut uncertainty = 0_f64;
        for refinement in [0, 1] {
            let mut refined = options.clone();
            if refinement == 0 {
                refined.cutoff_hartree *= 1.5;
            } else {
                refined.grid_override = Some(std::array::from_fn(|i| r.density.shape()[i] * 2));
            }
            match run(refined.clone()) {
                Ok(rr) => {
                    let t = frozen_stress(&ForceState {
                        options: &refined,
                        electronic: ElectronicState::Fixed(&rr),
                        ..state_view(&state)
                    })
                    .unwrap();
                    for (a, b) in t
                        .values
                        .iter()
                        .flatten()
                        .zip(tensor.values.iter().flatten())
                    {
                        uncertainty = uncertainty.max((a - b).abs());
                    }
                }
                Err(e) => {
                    println!("INCONCLUSIVE refinement: {e}");
                    uncertainty = f64::INFINITY;
                }
            }
        }
        // Gamma vs a full 2x2x2 quadrature is a separate uncertainty source.
        let points = dft_geometry::kpoints::monkhorst_pack([2; 3], [0.; 3]).unwrap();
        let mut kp = ScfProblem::with_kpoints(
            &cell,
            &atoms,
            &species,
            options.clone(),
            &points,
            dft_engine::OccupationPolicy::Insulator {
                minimum_gap_hartree: 1e-6,
            },
        )
        .unwrap();
        kp.set_density_floor(1e-10).unwrap();
        match kp.run() {
            Ok(rr) => {
                let t = frozen_stress(&ForceState {
                    electronic: ElectronicState::Fixed(&rr),
                    ..state_view(&state)
                })
                .unwrap();
                for (a, b) in t
                    .values
                    .iter()
                    .flatten()
                    .zip(tensor.values.iter().flatten())
                {
                    uncertainty = uncertainty.max((a - b).abs());
                }
            }
            Err(e) => {
                println!("INCONCLUSIVE k refinement: {e}");
                uncertainty = f64::INFINITY;
            }
        }
        println!(
            "REAL uncertainty={uncertainty:e} Ha/Bohr^3; {}",
            if uncertainty <= 1e-5 / 3. {
                "qualified for real cap"
            } else {
                "INCONCLUSIVE cutoff/k/grid convergence"
            }
        );
        let mut opts = options.clone();
        opts.strain_reference = Some(Box::new(cell.clone()));
        opts.grid_override = Some(r.density.shape().try_into().unwrap());
        for a in 0..3 {
            for b in a..3 {
                let mut errors = [0.; 3];
                for (i, h) in [0.002, 0.001, 0.0005].into_iter().enumerate() {
                    let plus = energy_at(&state, &strained(&cell, a, b, h), opts.clone());
                    let minus = energy_at(&state, &strained(&cell, a, b, -h), opts.clone());
                    if let (Ok(p), Ok(m)) = (plus, minus) {
                        errors[i] = ((p - m)
                            / (2. * h * cell.volume_bohr3() * if a == b { 1. } else { 2. })
                            - tensor.values[a][b])
                            .abs();
                    } else {
                        errors[i] = f64::INFINITY;
                    }
                }
                println!("REAL {material} {a}{b} errors={errors:?} Ha/Bohr^3");
                if uncertainty <= 1e-5 / 3. {
                    assert!(errors[2] < 1e-5);
                    assert!(errors[2] < 0.15 * errors[0] || errors.iter().all(|e| *e < 1e-5 / 3.));
                }
            }
        }
        let mut pressures = [0.; 3];
        for (i, h) in [0.002, 0.001, 0.0005].into_iter().enumerate() {
            let make = |h: f64| {
                Cell::from_user_rows_bohr(
                    cell.lattice_columns_bohr().map(|r| r.map(|v| v * (1. + h))),
                )
                .unwrap()
            };
            pressures[i] = -(energy_at(&state, &make(h), opts.clone()).unwrap()
                - energy_at(&state, &make(-h), opts.clone()).unwrap())
                / (6. * h * cell.volume_bohr3());
        }
        println!(
            "REAL pressure slopes={pressures:?} analytic={}",
            tensor.pressure
        );
    }
}
fn state_view<'a>(s: &ForceState<'a>) -> ForceState<'a> {
    ForceState {
        cell: s.cell,
        atoms: s.atoms,
        potentials: s.potentials,
        options: s.options,
        electronic: match s.electronic {
            ElectronicState::Fixed(r) => ElectronicState::Fixed(r),
            ElectronicState::Collinear(r) => ElectronicState::Collinear(r),
            ElectronicState::FiniteTemperature(r) => ElectronicState::FiniteTemperature(r),
        },
        density_floor: s.density_floor,
        symmetry: s.symmetry,
    }
}

#[test]
#[ignore = "release certified symmetry/preconditions and cutoff estimate"]
fn symmetry_preconditions_and_pulay() {
    let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [0., 8., 0.], [0., 0., 8.]]).unwrap();
    let pp = common::potential(2., 1.5, true, 401);
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let atoms = [Atom::new("Xx", [0.; 3])];
    let options = ScfOptions {
        cutoff_hartree: 1.,
        grid_override: Some([24; 3]),
        density_tolerance: 1e-10,
        energy_tolerance: 1e-12,
        eigen_tolerance: 1e-11,
        ..Default::default()
    };
    let points = [dft_geometry::kpoints::KPoint {
        fractional: [0.; 3],
        weight: 1.,
    }];
    let mut p = ScfProblem::new(&cell, &atoms, &species, options.clone()).unwrap();
    p.enable_symmetry(&points, 1e-7, false).unwrap();
    let mut r = p.run().unwrap();
    let operations = dft_geometry::symmetry::detect(&cell, &atoms, 1e-7).unwrap();
    let make = |r| ForceState {
        cell: &cell,
        atoms: &atoms,
        potentials: &species,
        options: &options,
        electronic: ElectronicState::Fixed(r),
        density_floor: None,
        symmetry: &operations,
    };
    let bits = r.total_hartree().to_bits();
    let density = r.density.clone();
    let stress = dft_engine::stress(&make(&r)).unwrap();
    assert_eq!(bits, r.total_hartree().to_bits());
    assert_eq!(density, r.density);
    for a in 0..3 {
        for b in 0..3 {
            assert!(
                (stress.values[a][b] - if a == b { stress.values[0][0] } else { 0. }).abs() < 1e-9
            );
        }
    }
    assert!(stress.antisymmetric_norm < 1e-10);
    assert_eq!(stress.ledger.len(), 8);
    assert!(
        stress
            .ledger
            .iter()
            .all(|e| !e.active || e.stress_available)
    );
    let estimate = stress.pulay.unwrap();
    println!(
        "Pulay correction={} uncertainty={} warning={}",
        estimate.correction, estimate.uncertainty, estimate.warning
    );
    assert!(
        estimate.warning,
        "deliberately low-cutoff fixture must warn"
    );
    let mut standard = options.clone();
    standard.standard = Some(dft_engine::StandardProfile { threshold_ha: 1e-6 });
    assert!(
        frozen_stress(&ForceState {
            options: &standard,
            ..make(&r)
        })
        .unwrap_err()
        .to_string()
        .contains("strict")
    );
    let altered = [Atom::new("Xx", [0.01, 0., 0.])];
    assert!(
        frozen_stress(&ForceState {
            atoms: &altered,
            ..make(&r)
        })
        .unwrap_err()
        .to_string()
        .contains("model")
    );
    let mut changed = options.clone();
    changed.strain_reference = Some(Box::new(cell.clone()));
    assert!(
        frozen_stress(&ForceState {
            options: &changed,
            ..make(&r)
        })
        .unwrap_err()
        .to_string()
        .contains("model")
    );
    r.certified = false;
    assert!(
        frozen_stress(&ForceState {
            cell: &cell,
            atoms: &atoms,
            potentials: &species,
            options: &options,
            electronic: ElectronicState::Fixed(&r),
            density_floor: None,
            symmetry: &operations
        })
        .unwrap_err()
        .to_string()
        .contains("certificate")
    );
}

#[test]
#[ignore = "release fixed-moment finite-T stress oracle"]
fn spin_fd_free_energy() {
    use dft_engine::{
        Functional,
        collinear::{CollinearProblem, SpinSector},
    };
    let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [1., 7., 0.], [0.5, 1., 9.]]).unwrap();
    let pp = common::potential_with_functional(1., 1.5, true, 401, "SLA PW");
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let atoms = [
        Atom::new("Xx", [0.; 3]),
        Atom::new("Xx", [0.27, 0.23, 0.19]),
    ];
    let options = ScfOptions {
        functional: Functional::LdaPw,
        cutoff_hartree: 2.,
        grid_override: Some([24; 3]),
        density_tolerance: 1e-11,
        energy_tolerance: 1e-13,
        eigen_tolerance: 1e-12,
        davidson_subspace: 96,
        max_iterations: 500,
        ..Default::default()
    };
    let points = [dft_geometry::kpoints::KPoint {
        fractional: [0.; 3],
        weight: 1.,
    }];
    let sector = SpinSector {
        n_up: 1,
        n_down: 1,
        axis_cartesian: [0., 0., 1.],
        initial_moment_seed: None,
    };
    let r = CollinearProblem::new_fermi_dirac(
        &cell,
        &atoms,
        &species,
        options.clone(),
        &points,
        sector,
        0.05,
        dft_physics::occupations::Constraint::FixedMoment(0.4),
        48,
    )
    .unwrap()
    .run()
    .unwrap();
    assert!((r.moment_mu_b - 0.4).abs() < 1e-10);
    assert!(r.thermodynamics.unwrap().entropy_term_ha > 1e-6);
    println!("case=spin-fd");
    println!(
        "TS={} Ha M={} electrons",
        r.thermodynamics.unwrap().entropy_term_ha,
        r.moment_mu_b
    );
    oracle(&ForceState {
        cell: &cell,
        atoms: &atoms,
        potentials: &species,
        options: &options,
        electronic: ElectronicState::Collinear(&r),
        density_floor: None,
        symmetry: &[],
    });
}

#[test]
#[ignore = "release polarized PBE channel-swap pair"]
fn spin_pbe_swap() {
    use dft_engine::{
        Functional,
        collinear::{CollinearProblem, SpinSector},
    };
    let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [1., 7., 0.], [0.5, 1., 9.]]).unwrap();
    let pp = synthetic(true, "PBE");
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let atoms = [
        Atom::new("Xx", [0.; 3]),
        Atom::new("Xx", [0.27, 0.23, 0.19]),
    ];
    let options = ScfOptions {
        functional: Functional::Pbe,
        cutoff_hartree: 2.,
        grid_override: Some([48; 3]),
        density_tolerance: 1e-11,
        energy_tolerance: 1e-13,
        eigen_tolerance: 1e-12,
        max_iterations: 500,
        ..Default::default()
    };
    let points = [dft_geometry::kpoints::KPoint {
        fractional: [0.; 3],
        weight: 1.,
    }];
    let mut tensors = Vec::new();
    for (n_up, n_down) in [(3, 1), (1, 3)] {
        let sector = SpinSector {
            n_up,
            n_down,
            axis_cartesian: [0., 0., 1.],
            initial_moment_seed: None,
        };
        let r = CollinearProblem::new(&cell, &atoms, &species, options.clone(), &points, sector)
            .unwrap()
            .run()
            .unwrap();
        tensors.push(
            frozen_stress(&ForceState {
                cell: &cell,
                atoms: &atoms,
                potentials: &species,
                options: &options,
                electronic: ElectronicState::Collinear(&r),
                density_floor: None,
                symmetry: &[],
            })
            .unwrap()
            .values,
        );
    }
    let error = tensors[0]
        .iter()
        .flatten()
        .zip(tensors[1].iter().flatten())
        .map(|(a, b)| (a - b).abs())
        .fold(0_f64, f64::max);
    println!("polarized PBE+NLCC swap max error={error:e} Ha/Bohr^3");
    assert!(error < 1e-9);
}
