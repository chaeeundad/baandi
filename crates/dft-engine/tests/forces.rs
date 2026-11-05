#![allow(clippy::needless_range_loop)] // Explicit atom/Cartesian axes in independent oracles.
#[allow(dead_code)]
mod common;
use common::potential;
use dft_engine::{
    ScfOptions, ScfProblem, SpeciesPotential,
    forces::{ForceState, forces},
};
use dft_geometry::{Atom, Cell};
use dft_physics::ewald::{self, EwaldOptions, Ion};

fn trend(values: [f64; 3], force: f64) {
    let errors = values.map(|v| (v - force).abs());
    assert!(
        errors[2] < 1e-6,
        "forces={force:e}, differences={values:?}, errors={errors:?} Ha/Bohr"
    );
    assert!(
        errors[2] < errors[0] * 0.15 || errors.iter().all(|e| *e < 1e-9),
        "missing second-order trend: {errors:?}"
    );
    if force.abs() > 1e-6 {
        assert_eq!(values[2].is_sign_positive(), force.is_sign_positive());
    }
}
#[test]
fn ewald_all_atoms_axes_and_eta() {
    for rows in [
        [[8., 0., 0.], [0., 8., 0.], [0., 0., 8.]],
        [[8., 0., 0.], [1., 7., 0.], [0.5, 1., 9.]],
    ] {
        let cell = Cell::from_user_rows_bohr(rows).unwrap();
        let ions = [
            Ion {
                charge: 2.,
                fractional_position: [0.03, 0.02, 0.01],
            },
            Ion {
                charge: 1.,
                fractional_position: [0.27, 0.23, 0.19],
            },
        ];
        for eta in [0.25, 0.5] {
            let opts = EwaldOptions {
                eta: Some(eta),
                ..Default::default()
            };
            let force = ewald::forces(&cell, &ions, opts).unwrap();
            for a in 0..3 {
                assert!((force[0][a] + force[1][a]).abs() < 1e-12);
            }
            for atom in 0..2 {
                for axis in 0..3 {
                    let mut fd = [0.; 3];
                    for (i, h) in [0.01, 0.005, 0.0025].into_iter().enumerate() {
                        let shift = cell.cartesian_bohr_to_fractional(std::array::from_fn(|a| {
                            if a == axis { h } else { 0. }
                        }));
                        let mut plus = ions;
                        let mut minus = ions;
                        for a in 0..3 {
                            plus[atom].fractional_position[a] += shift[a];
                            minus[atom].fractional_position[a] -= shift[a];
                        }
                        fd[i] = -(ewald::energy(&cell, &plus, opts).unwrap().energy.0
                            - ewald::energy(&cell, &minus, opts).unwrap().energy.0)
                            / (2. * h);
                    }
                    trend(fd, force[atom][axis]);
                    println!(
                        "Ewald eta={eta} atom={atom} axis={axis} error={:e} Ha/Bohr",
                        (fd[2] - force[atom][axis]).abs()
                    );
                }
            }
        }
    }
}
#[test]
#[ignore = "release-mode converged displacement scan; cargo test --release -p dft-engine --test forces converged -- --ignored --nocapture"]
fn converged_local_and_nonlocal() {
    for cutoff in [2., 3.] {
        for nonlocal in [false, true] {
            let cell =
                Cell::from_user_rows_bohr([[8., 0., 0.], [0., 8., 0.], [0., 0., 8.]]).unwrap();
            let pp = potential(1., 1.5, nonlocal, 801);
            let species = [SpeciesPotential {
                label: "Xx",
                potential: &pp,
            }];
            let atoms = vec![
                Atom::new("Xx", [0., 0., 0.]),
                Atom::new("Xx", [0.27, 0.23, 0.19]),
            ];
            let opts = ScfOptions {
                cutoff_hartree: cutoff,
                grid_override: Some([if cutoff == 2. { 48 } else { 96 }; 3]),
                density_tolerance: 1e-10,
                energy_tolerance: 1e-12,
                eigen_tolerance: 1e-11,
                max_iterations: 300,
                ..Default::default()
            };
            let run = |atoms: &[Atom]| {
                ScfProblem::new(&cell, atoms, &species, opts.clone())
                    .unwrap()
                    .run()
                    .unwrap()
            };
            let mut state = run(&atoms);
            let bits = state.total_hartree().to_bits();
            let force = forces(&ForceState {
                cell: &cell,
                atoms: &atoms,
                potentials: &species,
                options: &opts,
                electronic: dft_engine::forces::ElectronicState::Fixed(&state),
                density_floor: None,
                symmetry: &[],
            })
            .unwrap();
            assert_eq!(bits, state.total_hartree().to_bits());
            state.certified = false;
            assert!(
                forces(&ForceState {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &species,
                    options: &opts,
                    electronic: dft_engine::forces::ElectronicState::Fixed(&state),
                    density_floor: None,
                    symmetry: &[]
                })
                .unwrap_err()
                .to_string()
                .contains("certificate")
            );
            state.certified = true;
            let mut standard = opts.clone();
            standard.standard = Some(dft_engine::StandardProfile { threshold_ha: 1e-6 });
            assert!(
                forces(&ForceState {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &species,
                    options: &standard,
                    electronic: dft_engine::forces::ElectronicState::Fixed(&state),
                    density_floor: None,
                    symmetry: &[]
                })
                .unwrap_err()
                .to_string()
                .contains("strict convergence")
            );

            assert!(
                force.net_force.iter().all(|f| f.abs() < 1e-8),
                "net force {:?}",
                force.net_force
            );
            for atom in 0..2 {
                for axis in 0..3 {
                    let mut fd = [0.; 3];
                    for (i, h) in [0.01, 0.005, 0.0025].into_iter().enumerate() {
                        let shift = cell.cartesian_bohr_to_fractional(std::array::from_fn(|a| {
                            if a == axis { h } else { 0. }
                        }));
                        let mut plus = atoms.clone();
                        let mut minus = atoms.clone();
                        for a in 0..3 {
                            plus[atom].fractional_position[a] += shift[a];
                            minus[atom].fractional_position[a] -= shift[a];
                        }
                        fd[i] =
                            -(run(&plus).total_hartree() - run(&minus).total_hartree()) / (2. * h);
                    }
                    trend(fd, force.values[atom][axis]);
                    println!(
                        "converged cutoff={cutoff} nonlocal={nonlocal} atom={atom} axis={axis} force={:e} error={:e} Ha/Bohr",
                        force.values[atom][axis],
                        (fd[2] - force.values[atom][axis]).abs()
                    );
                }
            }
        }
    }
}

#[test]
fn frozen_local_and_coupled_projectors() {
    use dft_basis::PlaneWaveBasis;
    use dft_physics::{
        data::*,
        forces as partial,
        hamiltonian::{self, NonlocalOperator, Species},
    };
    let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [1., 7., 0.], [0.5, 1., 9.]]).unwrap();
    let pp = potential(1., 1.5, true, 401);
    let basis = PlaneWaveBasis::at_k(&cell, Hartree(1.), [0.12, 0.07, 0.03]).unwrap();
    let density_basis = PlaneWaveBasis::new(&cell, Hartree(4.)).unwrap();
    let positions = [[0.2, 0.1, 0.3], [2.3, 1.7, 1.5]];
    let species = [Species {
        potential: &pp,
        positions_bohr: &positions,
    }];
    let field: Vec<_> = density_basis
        .vectors()
        .iter()
        .map(|g| {
            Complex64::from_polar(
                0.001 * (-g.norm_squared).exp(),
                0.2 * g.cartesian[0] + 0.1 * g.cartesian[1],
            )
        })
        .collect();
    let mut orbitals = ComplexMatrix::zeros(basis.vectors().len(), 2).unwrap();
    for j in 0..2 {
        for i in 0..orbitals.rows() {
            orbitals[(i, j)] =
                Complex64::new(((i + 3 * j) as f64).sin(), ((2 * i + j) as f64).cos())
                    / (basis.vectors().len() as f64).sqrt();
        }
    }

    let occupations = [0.7, 0.3];
    let operator = NonlocalOperator::new(&basis, cell.volume_bohr3(), &species).unwrap();
    let nl = partial::nonlocal(&operator, &basis, &orbitals, &occupations).unwrap();
    let loc = partial::local_or_core(&density_basis, cell.volume_bohr3(), &species, &field, false)
        .unwrap();
    let energy = |positions: &[[f64; 3]; 2]| {
        let species = [Species {
            potential: &pp,
            positions_bohr: positions,
        }];
        let local =
            hamiltonian::ionic_local(&density_basis, cell.volume_bohr3(), &species).unwrap();
        let el = cell.volume_bohr3()
            * local
                .iter()
                .zip(&field)
                .map(|(v, n)| (n.conj() * v).re)
                .sum::<f64>();
        let op = NonlocalOperator::new(&basis, cell.volume_bohr3(), &species).unwrap();
        let mut en = 0.;
        for b in 0..2 {
            let x: Vec<_> = (0..orbitals.rows()).map(|i| orbitals[(i, b)]).collect();
            let mut y = vec![Complex64::default(); x.len()];
            let mut scratch = vec![Complex64::default(); op.scratch_len()];
            op.apply_into(&x, &mut y, &mut scratch).unwrap();
            en += occupations[b] * x.iter().zip(y).map(|(x, y)| (x.conj() * y).re).sum::<f64>();
        }
        [el, en]
    };
    for atom in 0..2 {
        for axis in 0..3 {
            let mut fd = [[0.; 3]; 2];
            for (i, h) in [0.01, 0.005, 0.0025].into_iter().enumerate() {
                let mut plus = positions;
                let mut minus = positions;
                plus[atom][axis] += h;
                minus[atom][axis] -= h;
                let ep = energy(&plus);
                let em = energy(&minus);
                for term in 0..2 {
                    fd[term][i] = -(ep[term] - em[term]) / (2. * h);
                }
            }
            trend(fd[0], loc[atom][axis]);
            trend(fd[1], nl[atom][axis]);
            println!(
                "frozen atom={atom} axis={axis} local_error={:e} nonlocal_error={:e} Ha/Bohr",
                (fd[0][2] - loc[atom][axis]).abs(),
                (fd[1][2] - nl[atom][axis]).abs()
            );
        }
    }
}

#[test]
fn ledger_lists_every_missing_derivative() {
    use dft_engine::forces::*;
    let mut ledger = derivative_ledger();
    for e in &mut ledger {
        if ["local", "core"].contains(&e.term_id) {
            e.force_available = false;
        }
    }
    let error = require_complete(&ledger).unwrap_err();
    assert!(error.contains("local, core"));
    assert!(ledger.iter().all(|e| !e.stress_available));
}

use dft_physics::data::erf;
use std::f64::consts::PI;
fn synthetic_text(core: bool, functional: &str) -> String {
    let n = 2001;
    let h = 0.006;
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
#[ignore = "release-mode spin/NLCC/free-energy scans; cargo test --release -p dft-engine --test forces spin_core_and_finite -- --ignored --nocapture"]
fn spin_core_and_finite() {
    use dft_engine::{
        collinear::{CollinearProblem, SpinSector},
        finite_temperature::FiniteTProblem,
        forces::ElectronicState,
    };
    use dft_geometry::kpoints::KPoint;
    for mode in ["nlcc", "pbe", "spin", "fd"] {
        if std::env::var("DFT_FORCE_CASE").is_ok_and(|selected| selected != mode) {
            continue;
        }
        let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [1., 7., 0.], [0.5, 1., 9.]]).unwrap();
        let pp = if mode == "pbe" {
            synthetic(true, "PBE")
        } else if mode == "nlcc" {
            synthetic(true, "SLA PW")
        } else {
            common::potential_with_functional(1., 1.5, true, 801, "SLA PW")
        };
        let species = [SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let atoms = vec![
            Atom::new("Xx", [0., 0., 0.]),
            Atom::new("Xx", [0.27, 0.23, 0.19]),
        ];
        let opts = ScfOptions {
            cutoff_hartree: 2.,
            grid_override: Some(if mode == "pbe" {
                [96, 96, 192]
            } else {
                [48; 3]
            }),
            functional: if mode == "pbe" {
                dft_engine::Functional::Pbe
            } else {
                dft_engine::Functional::LdaPw
            },
            density_tolerance: 1e-10,
            energy_tolerance: 1e-12,
            eigen_tolerance: 1e-11,
            davidson_subspace: if mode == "fd" { 96 } else { 48 },
            max_iterations: 400,
            ..Default::default()
        };
        let points = [KPoint {
            fractional: [0.; 3],
            weight: 1.,
        }];
        let run = |atoms: &[Atom],
                   get_force: bool|
         -> (f64, Option<dft_engine::forces::ForceResult>) {
            let evaluate = |electronic| {
                if get_force {
                    Some(
                        forces(&ForceState {
                            cell: &cell,
                            atoms,
                            potentials: &species,
                            options: &opts,
                            electronic,
                            density_floor: None,
                            symmetry: &[],
                        })
                        .unwrap(),
                    )
                } else {
                    None
                }
            };
            if mode == "spin" {
                let r = CollinearProblem::new(
                    &cell,
                    atoms,
                    &species,
                    opts.clone(),
                    &points,
                    SpinSector {
                        n_up: 1,
                        n_down: 1,
                        axis_cartesian: [0., 0., 1.],
                        initial_moment_seed: None,
                    },
                )
                .unwrap()
                .run()
                .unwrap();
                let spin_force = evaluate(ElectronicState::Collinear(&r));
                if get_force {
                    let scalar = ScfProblem::new(&cell, atoms, &species, opts.clone())
                        .unwrap()
                        .run()
                        .unwrap();
                    let scalar_force = evaluate(ElectronicState::Fixed(&scalar)).unwrap();
                    let difference = scalar_force
                        .values
                        .iter()
                        .flatten()
                        .zip(spin_force.as_ref().unwrap().values.iter().flatten())
                        .map(|(a, b)| (a - b).abs())
                        .fold(0., f64::max);
                    assert!(
                        difference < 1e-8,
                        "spin/unpolarized force mismatch {difference:e}"
                    );
                    println!("spin reduction force_difference={difference:e} Ha/Bohr");
                }
                (r.energy.total_hartree(), spin_force)
            } else if mode == "fd" {
                let r =
                    FiniteTProblem::new(&cell, atoms, &species, opts.clone(), &points, 0.05, 48)
                        .unwrap()
                        .run()
                        .unwrap();
                if get_force {
                    let fractions = r.occupations.occupations();
                    assert!(
                        fractions.iter().any(|f| *f > 1e-3 && *f < 2. - 1e-3),
                        "FD fixture must have resolved fractional occupations"
                    );
                    assert!(
                        r.thermodynamics.entropy_term_ha > 1e-5,
                        "FD fixture must distinguish F from E"
                    );
                    println!(
                        "FD theta=0.05 Ha TS={} Ha occupations={fractions:?}",
                        r.thermodynamics.entropy_term_ha
                    );
                }
                (
                    r.thermodynamics.free_energy_ha,
                    evaluate(ElectronicState::FiniteTemperature(&r)),
                )
            } else {
                let r = ScfProblem::new(&cell, atoms, &species, opts.clone())
                    .unwrap()
                    .run()
                    .unwrap();
                (r.total_hartree(), evaluate(ElectronicState::Fixed(&r)))
            }
        };
        let (_, force) = run(&atoms, true);
        let force = force.unwrap();
        println!("{mode} net_force={:?} Ha/Bohr", force.net_force);
        assert!(
            force.net_force.iter().all(|f| f.abs() < 1e-8),
            "{mode} net {:?}",
            force.net_force
        );
        for atom in 0..2 {
            for axis in 0..3 {
                let mut fd = [0.; 3];
                for (i, h) in [0.01, 0.005, 0.0025].into_iter().enumerate() {
                    let shift = cell.cartesian_bohr_to_fractional(std::array::from_fn(|a| {
                        if a == axis { h } else { 0. }
                    }));
                    let mut plus = atoms.clone();
                    let mut minus = atoms.clone();
                    for a in 0..3 {
                        plus[atom].fractional_position[a] += shift[a];
                        minus[atom].fractional_position[a] -= shift[a];
                    }
                    fd[i] = -(run(&plus, false).0 - run(&minus, false).0) / (2. * h);
                }
                trend(fd, force.values[atom][axis]);
                println!(
                    "{mode} atom={atom} axis={axis} error={:e} Ha/Bohr",
                    (fd[2] - force.values[atom][axis]).abs()
                );
            }
        }
    }
}

#[test]
fn independent_molecule_ewald_and_local_energy() {
    // Test-side direct cubic sums: fixed integer boxes, independently chosen
    // screening eta, explicit trigonometric pair sum (no production Ewald calls).
    use dft_basis::PlaneWaveBasis;
    use dft_physics::{data::*, forces as partial, hamiltonian::Species};
    let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [0., 8., 0.], [0., 0., 8.]]).unwrap();
    let positions = [[0.2, 0.3, 0.4], [2.3, 1.7, 1.6]];
    let pp = potential(1., 1.5, false, 801);
    let basis = PlaneWaveBasis::new(&cell, Hartree(4.)).unwrap();
    let field: Vec<_> = basis
        .vectors()
        .iter()
        .map(|g| Complex64::from_polar(0.001 * (-g.norm_squared).exp(), 0.2 * g.cartesian[0]))
        .collect();
    let species = [Species {
        potential: &pp,
        positions_bohr: &positions,
    }];
    let local = partial::local_or_core(&basis, 512., &species, &field, false).unwrap();
    let ions: Vec<_> = positions
        .iter()
        .map(|p| Ion {
            charge: 1.,
            fractional_position: p.map(|x| x / 8.),
        })
        .collect();
    let ew = ewald::forces(&cell, &ions, EwaldOptions::default()).unwrap();
    let energy = |positions: [[f64; 3]; 2]| {
        let eta = 0.4;
        let mut e = -2. * eta / PI.sqrt() - PI * 4. / (2. * 512. * eta * eta);
        for i in -6..=6 {
            for j in -6..=6 {
                for k in -6..=6 {
                    let m = [i, j, k];
                    for a in 0..2 {
                        for b in 0..2 {
                            if a == b && m == [0; 3] {
                                continue;
                            }
                            let r = (0..3)
                                .map(|d| {
                                    (positions[a][d] - positions[b][d] + 8. * m[d] as f64).powi(2)
                                })
                                .sum::<f64>()
                                .sqrt();
                            e += 0.5 * erfc(eta * r) / r;
                        }
                    }
                    if m == [0; 3] {
                        continue;
                    }
                    let g = m.map(|n| 2. * PI * n as f64 / 8.);
                    let g2 = g.iter().map(|x| x * x).sum::<f64>();
                    let phase = (0..3)
                        .map(|d| g[d] * (positions[0][d] - positions[1][d]))
                        .sum::<f64>();
                    e += 2. * PI / 512. * (-g2 / (4. * eta * eta)).exp() / g2
                        * (2. + 2. * phase.cos());
                }
            }
        }
        for (g, n) in basis.vectors().iter().zip(&field) {
            if g.miller == [0; 3] {
                continue;
            }
            // Analytic Gaussian-screened Coulomb ion, independent radial transform.
            let v = -4. * PI * (-1.5_f64.powi(2) * g.norm_squared / 4.).exp() / g.norm_squared;
            for p in positions {
                let angle = -(0..3).map(|d| g.cartesian[d] * p[d]).sum::<f64>();
                e += v * (n.re * angle.cos() + n.im * angle.sin());
            }
        }
        e
    };
    for atom in 0..2 {
        for axis in 0..3 {
            let mut fd = [0.; 3];
            for (i, h) in [0.01, 0.005, 0.0025].into_iter().enumerate() {
                let mut plus = positions;
                let mut minus = positions;
                plus[atom][axis] += h;
                minus[atom][axis] -= h;
                fd[i] = -(energy(plus) - energy(minus)) / (2. * h);
            }
            trend(fd, ew[atom][axis] + local[atom][axis]);
            println!(
                "independent atom={atom} axis={axis} error={:e} Ha/Bohr",
                (fd[2] - ew[atom][axis] - local[atom][axis]).abs()
            );
        }
    }
}

#[test]
#[ignore = "external real-file convergence/force study; DFT_NC_DIR=... cargo test --release -p dft-engine --test forces real_files -- --ignored --nocapture"]
fn real_files() {
    let dir = std::path::PathBuf::from(std::env::var_os("DFT_NC_DIR").expect("set DFT_NC_DIR"));
    for material in ["Si", "MgO"] {
        let labels = if material == "Si" {
            ["Si", "Si"]
        } else {
            ["Mg", "O"]
        };
        let pp0 = dft_upf::from_file(dir.join(format!("{}.upf", labels[0]))).unwrap();
        let pp1 = dft_upf::from_file(dir.join(format!("{}.upf", labels[1]))).unwrap();
        let mut species = vec![
            SpeciesPotential {
                label: labels[0],
                potential: &pp0,
            },
            SpeciesPotential {
                label: labels[1],
                potential: &pp1,
            },
        ];
        if material == "Si" {
            species.truncate(1);
        }
        let a = if material == "Si" { 10.2612 } else { 7.96 };
        let cell = Cell::from_user_rows_bohr([
            [0., a / 2., a / 2.],
            [a / 2., 0., a / 2.],
            [a / 2., a / 2., 0.],
        ])
        .unwrap();
        let mut atoms = vec![
            Atom::new(labels[0], [0.; 3]),
            Atom::new(labels[1], [if material == "Si" { 0.25 } else { 0.5 }; 3]),
        ];
        let shift = cell.cartesian_bohr_to_fractional([0.03, 0.02, 0.01]);
        for d in 0..3 {
            atoms[1].fractional_position[d] += shift[d];
        }
        let opts = ScfOptions {
            cutoff_hartree: if material == "Si" { 20. } else { 40. },
            functional: dft_engine::Functional::Pbe,
            density_tolerance: 1e-9,
            energy_tolerance: 1e-11,
            eigen_tolerance: 1e-10,
            max_iterations: 200,
            ..Default::default()
        };
        let run = |atoms: &[Atom], opts: &ScfOptions| {
            let mut p = ScfProblem::new(&cell, atoms, &species, opts.clone()).unwrap();
            p.set_density_floor(1e-10).unwrap();
            p.run()
                .unwrap_or_else(|e| panic!("INCONCLUSIVE {material}: SCF failed: {e}"))
        };
        let state = run(&atoms, &opts);
        let eval = |r: &dft_engine::ScfResult, opts: &ScfOptions| {
            forces(&ForceState {
                cell: &cell,
                atoms: &atoms,
                potentials: &species,
                options: opts,
                electronic: dft_engine::forces::ElectronicState::Fixed(r),
                density_floor: Some(1e-10),
                symmetry: &[],
            })
            .unwrap()
        };
        let force = eval(&state, &opts);
        let refined = ScfOptions {
            cutoff_hartree: if material == "Si" { 30. } else { 60. },
            ..opts.clone()
        };
        let reference = eval(&run(&atoms, &refined), &refined);
        let mut grid_options = opts.clone();
        grid_options.grid_override = Some(std::array::from_fn(|i| {
            state.density.shape()[i] * if material == "MgO" { 4 } else { 2 }
        }));
        let grid_force = eval(&run(&atoms, &grid_options), &grid_options);
        println!(
            "real {material} refined_grid_net_force={:?} Ha/Bohr",
            grid_force.net_force
        );
        assert!(
            grid_force.net_force.iter().all(|f| f.abs() < 1e-8),
            "INCONCLUSIVE {material}: refined grid translation residual {:?}",
            grid_force.net_force
        );
        let tight_options = ScfOptions {
            density_tolerance: 1e-10,
            energy_tolerance: 1e-12,
            eigen_tolerance: 1e-11,
            ..opts.clone()
        };
        let tight_force = eval(&run(&atoms, &tight_options), &tight_options);
        let difference = |other: &dft_engine::ForceResult| {
            force
                .values
                .iter()
                .flatten()
                .zip(other.values.iter().flatten())
                .map(|(a, b)| (a - b).abs())
                .fold(0., f64::max)
        };
        let grid_uncertainty = difference(&grid_force);
        let scf_uncertainty = difference(&tight_force);
        println!(
            "real {material} grid_uncertainty={grid_uncertainty:e} scf_uncertainty={scf_uncertainty:e} net_force={:?} Ha/Bohr",
            force.net_force
        );
        let uncertainty = force
            .values
            .iter()
            .flatten()
            .zip(reference.values.iter().flatten())
            .map(|(a, b)| (a - b).abs())
            .fold(0., f64::max);
        println!("real {material} cutoff_uncertainty={uncertainty:e} Ha/Bohr");
        assert!(
            uncertainty + grid_uncertainty + scf_uncertainty < 1e-4 / 3.,
            "INCONCLUSIVE {material}: cutoff uncertainty {uncertainty:e} exceeds one third of cap; refine cutoff/grid before acceptance"
        );
        for atom in 0..2 {
            for axis in 0..3 {
                let mut fd = [0.; 3];
                for (i, h) in [0.01, 0.005, 0.0025].into_iter().enumerate() {
                    let shift = cell.cartesian_bohr_to_fractional(std::array::from_fn(|d| {
                        if d == axis { h } else { 0. }
                    }));
                    let mut plus = atoms.clone();
                    let mut minus = atoms.clone();
                    for d in 0..3 {
                        plus[atom].fractional_position[d] += shift[d];
                        minus[atom].fractional_position[d] -= shift[d];
                    }
                    fd[i] = -(run(&plus, &opts).total_hartree()
                        - run(&minus, &opts).total_hartree())
                        / (2. * h);
                }
                let error = (fd[2] - force.values[atom][axis]).abs();
                println!("real {material} atom={atom} axis={axis} error={error:e} Ha/Bohr");
                assert!(error < 1e-4);
                assert!((fd[2] - fd[1]).abs() < (fd[1] - fd[0]).abs() * 0.4 || error < 1e-6);
                if force.values[atom][axis].abs() > 1e-6 {
                    assert_eq!(
                        fd[2].is_sign_positive(),
                        force.values[atom][axis].is_sign_positive()
                    );
                }
            }
        }
    }
}
