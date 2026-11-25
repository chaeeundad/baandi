//! G5 derivative oracles. Synthetic fixtures are not QE reference evidence.
//! External Si2 matrix: cargo test --release -p dft-engine --test refscf_forces -- --include-ignored --nocapture
#[allow(dead_code)]
mod common;
use dft_engine::{
    ElectronicState, ForceState, GammaReal, ScfOptions, ScfProblem, SpeciesPotential, StartDensity,
    StartWavefunctions, Verification,
    refscf::{delivered, k_driver, profile::Profile},
    refstart::StartInput,
};
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use dft_numerics::switches::RunConfig;
use dft_physics::data::{Complex64, ComplexMatrix};

#[test]
fn band_first_identity_complex_k_fractional_occupations() {
    use dft_basis::PlaneWaveBasis;
    use dft_physics::hamiltonian::{NonlocalOperator, Species};
    let pp = common::potential(2., 1.8, true, 161);
    let positions = [[0.3, 0.7, 1.1], [2.1, 3.2, 1.7]];
    let species = [Species {
        potential: &pp,
        positions_bohr: &positions,
    }];
    for k in [[0.; 3], [0.2, -0.1, 0.3]] {
        let cell =
            Cell::from_user_rows_bohr([[7., 0., 0.], [0.3, 6., 0.], [0.2, 0.4, 8.]]).unwrap();
        let basis = PlaneWaveBasis::at_k(&cell, dft_units::Hartree(2.), k).unwrap();
        let operator = NonlocalOperator::new(&basis, cell.volume_bohr3(), &species).unwrap();
        let mut orbitals = ComplexMatrix::zeros(basis.vectors().len(), 5).unwrap();
        for (i, z) in orbitals.as_mut_slice().iter_mut().enumerate() {
            *z = Complex64::new((i as f64 * 0.13).sin(), (i as f64 * 0.37).cos()) / 20.;
        }
        let occupations = [2., 1.3, 0.7, 0., 0.01];
        let old =
            dft_physics::forces::nonlocal(&operator, &basis, &orbitals, &occupations).unwrap();
        let new =
            dft_engine::forces::nonlocal_band_first(&operator, &basis, &orbitals, &occupations)
                .unwrap();
        let error = delta(old.iter().flatten(), new.iter().flatten());
        println!("band-first k={k:?} max error={error:e} Ha/Bohr");
        assert!(error < 1e-12);
    }
}
fn delta<'a>(a: impl Iterator<Item = &'a f64>, b: impl Iterator<Item = &'a f64>) -> f64 {
    a.zip(b).map(|(a, b)| (a - b).abs()).fold(0., f64::max)
}
fn options() -> ScfOptions {
    ScfOptions {
        standard: Some(dft_engine::StandardProfile {
            threshold_ha: 5e-13,
        }),
        cutoff_hartree: 2.,
        bands: Some(3),
        gamma_real: GammaReal::On,
        start_density: StartDensity::Atomic,
        start_wavefunctions: StartWavefunctions::Atomic,
        max_iterations: 200,
        ..ScfOptions::default()
    }
}
fn compare(
    cell: &Cell,
    atoms: &[Atom],
    pp: &dft_upf::Pseudopotential,
    opts: ScfOptions,
    points: &[KPoint],
) {
    let species = [SpeciesPotential {
        label: "Xx",
        potential: pp,
    }];
    let mut strict = opts.clone();
    strict.standard = None;
    strict.verification = Verification::Off;
    strict.density_tolerance = 1e-9;
    strict.energy_tolerance = 1e-12;
    strict.eigen_tolerance = 1e-11;
    let legacy = ScfProblem::with_kpoints(
        cell,
        atoms,
        &species,
        strict.clone(),
        points,
        dft_engine::OccupationPolicy::Insulator {
            minimum_gap_hartree: 1e-8,
        },
    )
    .unwrap()
    .run()
    .unwrap();
    let old = ForceState {
        cell,
        atoms,
        potentials: &species,
        options: &strict,
        electronic: ElectronicState::Fixed(&legacy),
        density_floor: None,
        symmetry: &[],
    };
    let old_f = dft_engine::forces(&old).unwrap();
    let old_s = dft_engine::stress::frozen_stress(&old).unwrap();
    for ranks in [1, 3, 8] {
        let results = dft_mpi::simulate(ranks, |comm| {
            let mut run = k_driver::run_with(
                StartInput {
                    cell,
                    atoms,
                    potentials: &species,
                    options: &opts,
                },
                Profile::from_standard_ha(opts.standard.unwrap().threshold_ha).unwrap(),
                points,
                comm,
            )
            .unwrap();
            let energy = run.outcome.result.total_hartree();
            let view = delivered::force_state(
                StartInput {
                    cell,
                    atoms,
                    potentials: &species,
                    options: &opts,
                },
                &mut run.outcome,
            )
            .unwrap();
            let force = dft_engine::forces(&view).unwrap();
            let stress = dft_engine::stress::frozen_stress(&view).unwrap();
            let config = RunConfig::from_pairs(&[("force_band_first", "on")]).unwrap();
            let _scope = config.enter();
            let fast = dft_engine::forces(&view).unwrap();
            assert!(delta(force.values.iter().flatten(), fast.values.iter().flatten()) < 1e-12);
            (energy, force.values, stress.values)
        });
        if ranks == 1 {
            check_energy_derivatives(
                StartInput {
                    cell,
                    atoms,
                    potentials: &species,
                    options: &opts,
                },
                points,
                &results[0].1,
                &results[0].2,
            );
        }
        for (energy, f, s) in results {
            let df = delta(f.iter().flatten(), old_f.values.iter().flatten());
            let ds = delta(s.iter().flatten(), old_s.values.iter().flatten());
            println!(
                "ranks={ranks} k={} dE={:e} Ha dF={df:e} Ha/Bohr dS={ds:e} Ha/Bohr^3",
                points.len(),
                energy - legacy.total_hartree()
            );
            assert!((energy - legacy.total_hartree()).abs() < 5e-9);
            assert!(df < 1e-6 && ds < 1e-6);
        }
    }
}
fn check_energy_derivatives(
    input: StartInput<'_>,
    points: &[KPoint],
    force: &[[f64; 3]],
    stress: &[[f64; 3]; 3],
) {
    // Keep the quadrature fixed for both displacement and frozen-Miller strain.
    let db = dft_basis::PlaneWaveBasis::new(
        input.cell,
        dft_units::Hartree(4. * input.options.cutoff_hartree),
    )
    .unwrap();
    let mut dims = [0; 3];
    for point in points {
        let basis = dft_basis::PlaneWaveBasis::at_k(
            input.cell,
            dft_units::Hartree(input.options.cutoff_hartree),
            point.fractional,
        )
        .unwrap();
        let grid = dft_basis::FftGridSpec::new(&basis, &db, input.options.grid_override).unwrap();
        for (d, n) in dims.iter_mut().zip(grid.dims()) {
            *d = (*d).max(n);
        }
    }
    let mut options = input.options.clone();
    options.grid_override = Some(dims);
    let input = StartInput {
        options: &options,
        ..input
    };
    let cell = input.cell;
    let atoms = input.atoms;
    let profile = Profile::from_standard_ha(input.options.standard.unwrap().threshold_ha).unwrap();
    let h = 0.002;
    let mut energies = Vec::new();
    for sign in [-1., 1.] {
        let mut displaced = atoms.to_vec();
        let shift = cell.cartesian_bohr_to_fractional([sign * h, 0., 0.]);
        for (v, dx) in displaced[1].fractional_position.iter_mut().zip(shift) {
            *v += dx;
        }
        let r = k_driver::run(
            StartInput {
                atoms: &displaced,
                ..input
            },
            profile,
            points,
        )
        .unwrap();
        assert!(r.outcome.converged);
        energies.push(r.outcome.result.total_hartree());
    }
    let fd = -(energies[1] - energies[0]) / (2. * h);
    println!(
        "reference FD force error={:e} Ha/Bohr",
        (fd - force[1][0]).abs()
    );
    assert!((fd - force[1][0]).abs() < 1e-6);
    // A tiny strain away from cutoff shells preserves the Miller set; assert
    // that fact independently before using the reference energy difference.
    let h = 0.0001;
    let mut energies = Vec::new();
    for sign in [-1., 1.] {
        let mut columns = cell.lattice_columns_bohr();
        for v in &mut columns[0] {
            *v *= 1. + sign * h;
        }
        let rows = std::array::from_fn(|i| std::array::from_fn(|j| columns[j][i]));
        let strained = Cell::from_user_rows_bohr(rows).unwrap();
        let strained_db = dft_basis::PlaneWaveBasis::new(
            &strained,
            dft_units::Hartree(4. * input.options.cutoff_hartree),
        )
        .unwrap();
        let mut density_millers: Vec<_> = db.vectors().iter().map(|g| g.miller).collect();
        let mut strained_millers: Vec<_> = strained_db.vectors().iter().map(|g| g.miller).collect();
        density_millers.sort();
        strained_millers.sort();
        assert_eq!(
            density_millers, strained_millers,
            "strain crossed a density cutoff shell"
        );
        for point in points {
            let millers = |cell: &Cell| {
                let basis = dft_basis::PlaneWaveBasis::at_k(
                    cell,
                    dft_units::Hartree(input.options.cutoff_hartree),
                    point.fractional,
                )
                .unwrap();
                let mut list: Vec<_> = basis.vectors().iter().map(|g| g.miller).collect();
                list.sort();
                list
            };
            assert_eq!(
                millers(cell),
                millers(&strained),
                "strain crossed a cutoff shell"
            );
        }
        let r = k_driver::run(
            StartInput {
                cell: &strained,
                ..input
            },
            profile,
            points,
        )
        .unwrap();
        assert!(r.outcome.converged);
        energies.push(r.outcome.result.total_hartree());
    }
    let fd = (energies[1] - energies[0]) / (2. * h * cell.volume_bohr3());
    println!(
        "reference FD stress error={:e} Ha/Bohr^3",
        (fd - stress[0][0]).abs()
    );
    assert!((fd - stress[0][0]).abs() < 1e-6);
}

#[test]
fn synthetic_reference_legacy_ranks_and_kpoints() {
    let pp = common::potential_with_atomic_start(1., 1.5, true, 161, "SLA PZ");
    for skew in [0., 0.3] {
        let cell =
            Cell::from_user_rows_bohr([[7., 0., 0.], [skew, 7., 0.], [0., skew, 7.]]).unwrap();
        let atoms = [
            Atom::new("Xx", [0.13, 0.21, 0.32]),
            Atom::new("Xx", [0.40, 0.44, 0.51]),
        ];
        compare(
            &cell,
            &atoms,
            &pp,
            options(),
            &[KPoint {
                fractional: [0.; 3],
                weight: 1.,
            }],
        );
        let mut opts = options();
        opts.gamma_real = GammaReal::Off;
        compare(
            &cell,
            &atoms,
            &pp,
            opts,
            &[
                KPoint {
                    fractional: [0.25, 0., 0.],
                    weight: 0.5,
                },
                KPoint {
                    fractional: [-0.25, 0., 0.],
                    weight: 0.5,
                },
            ],
        );
    }
}
#[test]
#[ignore = "real Si2 SCF matrix; set DFT_SI_UPF and use the release command"]
fn silicon_displaced_and_skew() {
    let path = std::env::var("DFT_SI_UPF").expect("DFT_SI_UPF must name an external Si UPF");
    let pp = dft_upf::parse_upf(&std::fs::read_to_string(path).unwrap()).unwrap();
    let mut atoms = common::atoms();
    atoms[1].fractional_position[0] += 0.003;
    let mut opts = options();
    opts.cutoff_hartree = 4.;
    opts.bands = Some(5);
    compare(
        &common::cell(),
        &atoms,
        &pp,
        opts.clone(),
        &[KPoint {
            fractional: [0.; 3],
            weight: 1.,
        }],
    );
    let mut k_options = opts.clone();
    k_options.gamma_real = GammaReal::Off;
    compare(
        &common::cell(),
        &atoms,
        &pp,
        k_options,
        &[
            KPoint {
                fractional: [0.25, 0., 0.],
                weight: 0.5,
            },
            KPoint {
                fractional: [-0.25, 0., 0.],
                weight: 0.5,
            },
        ],
    );
    let cell = Cell::from_user_rows_bohr([[0.1, 5.13, 5.13], [5.13, 0., 5.13], [5.13, 5.13, 0.2]])
        .unwrap();
    compare(
        &cell,
        &atoms,
        &pp,
        opts,
        &[KPoint {
            fractional: [0.; 3],
            weight: 1.,
        }],
    );
}

#[test]
fn certificate_opt_out_rejection_and_reference_finite_difference() {
    let pp = common::potential_with_atomic_start(1., 1.5, true, 161, "SLA PZ");
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let cell = Cell::from_user_rows_bohr([[7., 0., 0.], [0.3, 7., 0.], [0., 0.2, 7.]]).unwrap();
    let atoms = [
        Atom::new("Xx", [0.1, 0.2, 0.3]),
        Atom::new("Xx", [0.37, 0.43, 0.49]),
    ];
    let mut opts = options();
    let points = [KPoint {
        fractional: [0.; 3],
        weight: 1.,
    }];
    let profile = Profile::from_standard_ha(5e-12).unwrap();
    opts.standard = Some(dft_engine::StandardProfile {
        threshold_ha: profile.threshold_ry / 2.,
    });
    let input = || StartInput {
        cell: &cell,
        atoms: &atoms,
        potentials: &species,
        options: &opts,
    };
    let mut run = k_driver::run(input(), profile, &points).unwrap();
    let delivered_density = run.outcome.result.density.clone();
    let delivered_energy = run.outcome.result.total_hartree();
    let mut changed = atoms.clone();
    changed[0].fractional_position[0] += 0.01;
    assert!(
        delivered::force_state(
            StartInput {
                atoms: &changed,
                ..input()
            },
            &mut run.outcome
        )
        .is_err()
    );
    {
        let view = delivered::force_state(input(), &mut run.outcome).unwrap();
        let f = dft_engine::forces(&view).unwrap();
        let s = dft_engine::stress::frozen_stress(&view).unwrap();
        assert_eq!(f.verification, "verified");
        check_energy_derivatives(input(), &points, &f.values, &s.values);
    }
    // Adapter and both derivatives preserve the stop map, including repeated delivery.
    delivered::stress_state(input(), &mut run.outcome).unwrap();
    assert_eq!(run.outcome.result.density, delivered_density);
    assert_eq!(
        run.outcome.result.total_hartree().to_bits(),
        delivered_energy.to_bits()
    );
    run.outcome.converged = false;
    assert!(delivered::force_state(input(), &mut run.outcome).is_err());
    run.outcome.converged = true;
    run.outcome
        .result
        .verification
        .as_mut()
        .unwrap()
        .attempts
        .last_mut()
        .unwrap()
        .passed = false;
    assert!(delivered::force_state(input(), &mut run.outcome).is_err());
    let mut off = opts.clone();
    off.verification = Verification::Off;
    assert!(
        delivered::force_state(
            StartInput {
                options: &off,
                ..input()
            },
            &mut run.outcome
        )
        .is_err()
    );
    {
        let config = RunConfig::from_pairs(&[("refscf_verify_energy", "on")]).unwrap();
        let _scope = config.enter();
        let mut energy_verified = k_driver::run(input(), profile, &points).unwrap();
        let view = delivered::force_state(input(), &mut energy_verified.outcome).unwrap();
        assert_eq!(dft_engine::forces(&view).unwrap().verification, "verified");
    }

    let mut run = k_driver::run(
        StartInput {
            options: &off,
            ..input()
        },
        profile,
        &points,
    )
    .unwrap();
    let view = delivered::force_state(
        StartInput {
            options: &off,
            ..input()
        },
        &mut run.outcome,
    )
    .unwrap();
    assert_eq!(
        dft_engine::forces(&view).unwrap().verification,
        "unverified"
    );
    assert_eq!(
        dft_engine::stress::frozen_stress(&view)
            .unwrap()
            .verification,
        "unverified"
    );
}

#[test]
fn distributed_ionic_seed_matches_cold_delivery() {
    let pp = common::potential_with_atomic_start(1., 1.5, true, 161, "SLA PZ");
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let cell = Cell::from_user_rows_bohr([[7., 0., 0.], [0.3, 7., 0.], [0., 0.2, 7.]]).unwrap();
    let old_atoms = [
        Atom::new("Xx", [0.1, 0.2, 0.3]),
        Atom::new("Xx", [0.37, 0.43, 0.49]),
    ];
    let mut atoms = old_atoms.clone();
    atoms[1].fractional_position[0] += 0.001;
    let opts = options();
    let profile = Profile::from_standard_ha(opts.standard.unwrap().threshold_ha).unwrap();
    let points = [KPoint {
        fractional: [0.; 3],
        weight: 1.,
    }];
    for ranks in [1, 3, 8] {
        dft_mpi::simulate(ranks, |comm| {
            let input = || StartInput {
                cell: &cell,
                atoms: &old_atoms,
                potentials: &species,
                options: &opts,
            };
            let mut initial = k_driver::run_with(input(), profile, &points, comm).unwrap();
            delivered::force_state(input(), &mut initial.outcome).unwrap();
            let density: Vec<_> = initial.outcome.result.density.iter().copied().collect();
            let seed =
                delivered::gamma_seed(input(), &initial.outcome.result, &density, comm).unwrap();
            let next = || StartInput {
                atoms: &atoms,
                ..input()
            };
            let mut warm =
                dft_engine::refscf::driver::run_seeded_with(next(), profile, comm, seed).unwrap();
            let mut cold = k_driver::run_with(next(), profile, &points, comm).unwrap();
            let de =
                (warm.outcome.result.total_hartree() - cold.outcome.result.total_hartree()).abs();
            let warm_view = delivered::force_state(next(), &mut warm.outcome).unwrap();
            let cold_view = delivered::force_state(next(), &mut cold.outcome).unwrap();
            let wf = dft_engine::forces(&warm_view).unwrap();
            let cf = dft_engine::forces(&cold_view).unwrap();
            let df = delta(wf.values.iter().flatten(), cf.values.iter().flatten());
            println!("ionic seed ranks={ranks} dE={de:e} Ha dF={df:e} Ha/Bohr");
            assert!(de < 5e-9 && df < 1e-6);
        });
    }
}

#[test]
fn lazy_delivery_materializes_once_and_preserves_derivatives() {
    use std::sync::Arc;
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let cell = Cell::from_user_rows_bohr([[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    let mut opts = options();
    opts.cutoff_hartree = 1.;
    let points = [KPoint {
        fractional: [0.; 3],
        weight: 1.,
    }];
    for ranks in [3, 8] {
        dft_mpi::simulate(ranks, |comm| {
            dft_engine::parallel::set_communicator(Some(Arc::new(comm.clone())));
            let input = || StartInput {
                cell: &cell,
                atoms: &atoms,
                potentials: &species,
                options: &opts,
            };
            let mut baseline = None;
            for lazy in ["off", "on"] {
                let config = RunConfig::from_pairs(&[
                    ("refscf_lazy_output", lazy),
                    ("m2_lazy_result", "on"),
                ])
                .unwrap();
                let _scope = config.enter();
                dft_engine::gamma_real::set_delivered_orbitals(false);
                let mut run = k_driver::run(
                    input(),
                    Profile::from_standard_ha(opts.standard.unwrap().threshold_ha).unwrap(),
                    &points,
                )
                .unwrap();
                assert_eq!(run.outcome.result.orbitals_delivered(), lazy == "off");
                let view = delivered::force_state(input(), &mut run.outcome).unwrap();
                let force = dft_engine::forces(&view).unwrap();
                let stress = dft_engine::stress::frozen_stress(&view).unwrap();
                let pointer = run.outcome.result.k_points[0]
                    .solution
                    .eigenvectors
                    .as_slice()
                    .as_ptr();
                delivered::stress_state(input(), &mut run.outcome).unwrap();
                assert_eq!(
                    pointer,
                    run.outcome.result.k_points[0]
                        .solution
                        .eigenvectors
                        .as_slice()
                        .as_ptr()
                );
                let science = (
                    run.outcome.result.total_hartree().to_bits(),
                    force.values,
                    stress.values,
                );
                if let Some(old) = &baseline {
                    assert_eq!(&science, old);
                } else {
                    baseline = Some(science);
                }
            }
            dft_engine::gamma_real::set_delivered_orbitals(true);
            dft_engine::parallel::set_communicator(None);
        });
    }
}
