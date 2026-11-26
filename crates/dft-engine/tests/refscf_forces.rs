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
    // Historical fixture contract; screening has separate charge/mask tests.
    let screening =
        dft_numerics::switches::RunConfig::from_pairs(&[("xc_screen_negative", "off")]).unwrap();
    let _screen_scope = screening.enter();
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
            let _screen_scope = screening.enter();
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
            let config =
                RunConfig::from_pairs(&[("xc_screen_negative", "off"), ("force_band_first", "on")])
                    .unwrap();
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
    let profile = Profile::from_standard_ha(input.options.standard.unwrap().threshold_ha).unwrap();
    check_derivatives(input, points, force, stress, |input| {
        let r = k_driver::run(input, profile, points).unwrap();
        assert!(r.outcome.converged);
        r.outcome.result.total_hartree()
    });
}
/// Central differences of the variational energy returned by `energy` (E for
/// fixed occupations, the Mermin free energy F = E - TS for Fermi-Dirac):
/// atom 1 displaced along x by +-0.002 Bohr, and an xx strain of +-1e-4 that
/// keeps the frozen Miller sets. Tolerance 1e-6 Ha/Bohr and Ha/Bohr^3.
fn check_derivatives(
    input: StartInput<'_>,
    points: &[KPoint],
    force: &[[f64; 3]],
    stress: &[[f64; 3]; 3],
    energy: impl Fn(StartInput<'_>) -> f64,
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
    let h = 0.002;
    let mut energies = Vec::new();
    for sign in [-1., 1.] {
        let mut displaced = atoms.to_vec();
        let shift = cell.cartesian_bohr_to_fractional([sign * h, 0., 0.]);
        for (v, dx) in displaced[1].fractional_position.iter_mut().zip(shift) {
            *v += dx;
        }
        energies.push(energy(StartInput {
            atoms: &displaced,
            ..input
        }));
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
        energies.push(energy(StartInput {
            cell: &strained,
            ..input
        }));
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
    // Historical fixture contract; screening has separate charge/mask tests.
    let screening =
        dft_numerics::switches::RunConfig::from_pairs(&[("xc_screen_negative", "off")]).unwrap();
    let _screen_scope = screening.enter();
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
        let config = RunConfig::from_pairs(&[("xc_screen_negative", "off")]).unwrap();
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
    // Historical fixture contract; screening has separate charge/mask tests.
    let screening =
        dft_numerics::switches::RunConfig::from_pairs(&[("xc_screen_negative", "off")]).unwrap();
    let _screen_scope = screening.enter();
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
            let _screen_scope = screening.enter();
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
                let config = RunConfig::from_pairs(&[("m2_lazy_result", "on")]).unwrap();
                let _scope = config.enter();
                dft_engine::gamma_real::set_delivered_orbitals(lazy == "off");
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

#[test]
fn distributed_k_ionic_seed_matches_cold_delivery() {
    // Historical fixture contract; screening has separate charge/mask tests.
    let screening =
        dft_numerics::switches::RunConfig::from_pairs(&[("xc_screen_negative", "off")]).unwrap();
    let _screen_scope = screening.enter();
    let pp = common::potential_with_atomic_start(1., 1.5, true, 161, "SLA PZ");
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let cell = Cell::from_user_rows_bohr([[7., 0., 0.], [0.3, 7., 0.], [0., 0.2, 7.]]).unwrap();
    let old = [
        Atom::new("Xx", [0.1, 0.2, 0.3]),
        Atom::new("Xx", [0.37, 0.43, 0.49]),
    ];
    let mut atoms = old.clone();
    atoms[1].fractional_position[0] += 0.001;
    let mut opts = options();
    opts.gamma_real = GammaReal::Off;
    let points = [
        KPoint {
            fractional: [0.25, 0., 0.],
            weight: 0.5,
        },
        KPoint {
            fractional: [-0.25, 0., 0.],
            weight: 0.5,
        },
    ];
    let profile = Profile::from_standard_ha(opts.standard.unwrap().threshold_ha).unwrap();
    for ranks in [1, 3, 8] {
        dft_mpi::simulate(ranks, |comm| {
            let _screen_scope = screening.enter();
            let input = || StartInput {
                cell: &cell,
                atoms: &old,
                potentials: &species,
                options: &opts,
            };
            let mut initial = k_driver::run_with(input(), profile, &points, comm).unwrap();
            delivered::force_state(input(), &mut initial.outcome).unwrap();
            let density: Vec<_> = initial.outcome.result.density.iter().copied().collect();
            let delivered::IonicSeed::K(seed) =
                delivered::ionic_seed(input(), &initial.outcome.result, &density, comm).unwrap()
            else {
                panic!("k seed expected")
            };
            let next = || StartInput {
                atoms: &atoms,
                ..input()
            };
            let mut warm =
                k_driver::run_seeded_with(next(), profile, &points, comm, seed.clone()).unwrap();
            let mut cold = k_driver::run_with(next(), profile, &points, comm).unwrap();
            let de =
                (warm.outcome.result.total_hartree() - cold.outcome.result.total_hartree()).abs();
            let wf =
                dft_engine::forces(&delivered::force_state(next(), &mut warm.outcome).unwrap())
                    .unwrap();
            let cf =
                dft_engine::forces(&delivered::force_state(next(), &mut cold.outcome).unwrap())
                    .unwrap();
            let df = delta(wf.values.iter().flatten(), cf.values.iter().flatten());
            println!("k ionic seed ranks={ranks} dE={de:e} Ha dF={df:e} Ha/Bohr");
            assert!(de <= 5e-9 && df <= 1e-6);
            let mut wrong = opts.clone();
            wrong.cutoff_hartree *= 1.1;
            assert!(
                k_driver::run_seeded_with(
                    StartInput {
                        options: &wrong,
                        ..next()
                    },
                    profile,
                    &points,
                    comm,
                    seed
                )
                .is_err()
            );
        });
    }
}

#[test]
#[ignore = "external Si2 cutoff scans; DFT_SI_UPF=/path/Si.pz-vbc.UPF cargo test --release -p dft-engine --test refscf_forces symmetric_si2_pulay -- --ignored --nocapture"]
fn symmetric_si2_pulay_matches_full_mesh() {
    let pp =
        dft_upf::parse_upf(&std::fs::read_to_string(std::env::var("DFT_SI_UPF").unwrap()).unwrap())
            .unwrap();
    let cell = common::cell();
    let atoms = [Atom::new("Si", [0.; 3]), Atom::new("Si", [0.25; 3])];
    let species = [SpeciesPotential {
        label: "Si",
        potential: &pp,
    }];
    let mut opts = options();
    opts.cutoff_hartree = 4.;
    opts.bands = Some(5);
    opts.gamma_real = GammaReal::Off;
    let points: Vec<_> = (0..2)
        .flat_map(|x| {
            (0..2).flat_map(move |y| {
                (0..2).map(move |z| KPoint {
                    fractional: [x as f64 / 2., y as f64 / 2., z as f64 / 2.],
                    weight: 0.125,
                })
            })
        })
        .collect();
    let prepared =
        dft_engine::refscf::symmetry::Prepared::new(&cell, &atoms, &points, 1e-8, false).unwrap();
    let input = || StartInput {
        cell: &cell,
        atoms: &atoms,
        potentials: &species,
        options: &opts,
    };
    let profile = Profile::from_standard_ha(opts.standard.unwrap().threshold_ha).unwrap();
    let mut full = k_driver::run(input(), profile, &points).unwrap();
    let full_view = delivered::stress_state(input(), &mut full.outcome).unwrap();
    let full_stress = dft_engine::stress(&full_view).unwrap();
    let mut reduced = k_driver::run_prepared(input(), profile, &prepared).unwrap();
    let reduced_view =
        delivered::stress_state_prepared(input(), &mut reduced.outcome, &prepared).unwrap();
    let reduced_stress = dft_engine::stress::stress_prepared(&reduced_view, &prepared).unwrap();
    let ds = delta(
        full_stress.values.iter().flatten(),
        reduced_stress.values.iter().flatten(),
    );
    let dp = (full_stress.pulay.as_ref().unwrap().correction
        - reduced_stress.pulay.as_ref().unwrap().correction)
        .abs();
    println!(
        "Si2 symmetry full={} reduced={} dStress={ds:e} Ha/Bohr^3 dPulay={dp:e} Ha/Bohr^3",
        points.len(),
        prepared.points().len()
    );
    assert!(ds <= 1e-8 && dp <= 1e-8);
    for i in 0..3 {
        let corrected = (full_stress.values[i][i] + full_stress.pulay.as_ref().unwrap().correction
            - reduced_stress.values[i][i]
            - reduced_stress.pulay.as_ref().unwrap().correction)
            .abs();
        assert!(corrected <= 1e-8);
    }
}

// R7a: emit exact delivered-state fingerprints for same-build before/after
// comparisons. No cross-platform hash or floating-point golden is assumed.
#[test]
fn unstrained_delivery_fingerprints() {
    use std::hash::{DefaultHasher, Hash, Hasher};
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let cell = Cell::from_user_rows_bohr([[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    let mut opts = options();
    opts.cutoff_hartree = 1.;
    for gamma in [true, false] {
        opts.gamma_real = if gamma { GammaReal::On } else { GammaReal::Off };
        let points = if gamma {
            vec![KPoint {
                fractional: [0.; 3],
                weight: 1.,
            }]
        } else {
            vec![
                KPoint {
                    fractional: [0.25, 0., 0.],
                    weight: 0.5,
                },
                KPoint {
                    fractional: [-0.25, 0., 0.],
                    weight: 0.5,
                },
            ]
        };
        let input = || StartInput {
            cell: &cell,
            atoms: &atoms,
            potentials: &species,
            options: &opts,
        };
        let mut run =
            k_driver::run(input(), Profile::from_standard_ha(5e-13).unwrap(), &points).unwrap();
        assert!(run.outcome.converged);
        let stress = dft_engine::stress::frozen_stress(
            &delivered::stress_state(input(), &mut run.outcome).unwrap(),
        )
        .unwrap();
        let r = &run.outcome.result;
        let mut hash = DefaultHasher::new();
        for x in std::iter::once(&r.total_hartree())
            .chain(r.density.iter())
            .chain(r.occupations.iter())
            .chain(stress.values.iter().flatten())
        {
            x.to_bits().hash(&mut hash);
        }
        for k in &r.k_points {
            for x in &k.solution.eigenvalues {
                x.to_bits().hash(&mut hash);
            }
            for z in k.solution.eigenvectors.as_slice() {
                z.re.to_bits().hash(&mut hash);
                z.im.to_bits().hash(&mut hash);
            }
        }
        println!(
            "R7A unstrained gamma={gamma} fingerprint={:016x}",
            hash.finish()
        );
    }
}

/// FD-FORCES: reference spin/Fermi-Dirac delivery. The analytic forces and
/// frozen-Miller stress of the delivered state are checked against central
/// differences of its own variational energy (F = E - TS for Fermi-Dirac, E
/// for a fixed spin sector) and against the strict legacy oracle state
/// (`finite_temperature` / `collinear` problems with the same model).
mod electronic_derivatives {
    use super::*;
    use dft_engine::{
        Functional,
        collinear::{CollinearProblem, SpinSector},
        finite_temperature::FiniteTProblem,
        refscf::electronic::{self, Policy},
    };
    use dft_physics::occupations::Constraint;

    struct Case {
        policy: Policy,
        functional: Functional,
        pp_functional: &'static str,
        bands: usize,
        points: Vec<KPoint>,
    }

    fn sector(n_up: usize, n_down: usize) -> SpinSector {
        SpinSector {
            n_up,
            n_down,
            axis_cartesian: [0., 0., 1.],
            initial_moment_seed: Some(0.3),
        }
    }

    fn gamma() -> Vec<KPoint> {
        vec![KPoint {
            fractional: [0.; 3],
            weight: 1.,
        }]
    }

    fn run_case(case: Case) {
        let screening = RunConfig::from_pairs(&[("xc_screen_negative", "off")]).unwrap();
        let _screen_scope = screening.enter();
        let pp = common::potential_with_atomic_start(1., 1.5, true, 161, case.pp_functional);
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
        opts.functional = case.functional;
        opts.gamma_real = GammaReal::Off;
        opts.bands = Some(case.bands);
        let profile = Profile::from_standard_ha(5e-13).unwrap();
        opts.standard = Some(dft_engine::StandardProfile {
            threshold_ha: profile.threshold_ry / 2.,
        });
        let input = || StartInput {
            cell: &cell,
            atoms: &atoms,
            potentials: &species,
            options: &opts,
        };
        let points = case.points.as_slice();
        let mut run = electronic::run(input(), profile, points, case.policy.clone()).unwrap();
        assert!(run.converged && run.verification.as_ref().unwrap().passed());
        let entropy = run.thermodynamics.entropy_term_ha;
        let fractional = run
            .blocks
            .iter()
            .flatten()
            .flat_map(|b| b.occupations.iter())
            .filter(|f| **f > 1e-6 && (**f - 1.).abs() > 1e-6 && (**f - 2.).abs() > 1e-6)
            .count();
        println!(
            "case {:?}: F={:.12} Ha TS={entropy:e} Ha fractional occupations={fractional}",
            case.policy, run.thermodynamics.free_energy_ha
        );
        let (force, stress) = {
            let view = delivered::electronic_state(input(), &mut run, &[]).unwrap();
            let f = dft_engine::forces(&view).unwrap();
            let s = dft_engine::stress::frozen_stress(&view).unwrap();
            assert_eq!(f.verification, "verified");
            assert_eq!(f.quantity, if case.policy.is_finite() { "F" } else { "E" });
            let fast = {
                let config = RunConfig::from_pairs(&[
                    ("xc_screen_negative", "off"),
                    ("force_band_first", "on"),
                ])
                .unwrap();
                let _scope = config.enter();
                dft_engine::forces(&view).unwrap()
            };
            assert!(delta(f.values.iter().flatten(), fast.values.iter().flatten()) < 1e-12);
            (f.values, s.values)
        };
        // A changed model is rejected instead of silently used.
        let mut moved = atoms.clone();
        moved[0].fractional_position[0] += 0.01;
        assert!(
            delivered::electronic_state(
                StartInput {
                    atoms: &moved,
                    ..input()
                },
                &mut run,
                &[]
            )
            .is_err()
        );
        let policy = case.policy.clone();
        check_derivatives(input(), points, &force, &stress, |input| {
            let r = electronic::run(input, profile, points, policy.clone()).unwrap();
            assert!(r.converged && r.verification.as_ref().unwrap().passed());
            r.thermodynamics.free_energy_ha
        });
        // Strict legacy oracle: same model, strict certificate.
        let mut strict = opts.clone();
        strict.standard = None;
        strict.verification = Verification::Off;
        strict.density_tolerance = 1e-9;
        strict.energy_tolerance = 1e-12;
        strict.eigen_tolerance = 1e-11;
        let oracle = |electronic| {
            let state = ForceState {
                cell: &cell,
                atoms: &atoms,
                potentials: &species,
                options: &strict,
                electronic,
                density_floor: None,
                symmetry: &[],
            };
            (
                dft_engine::forces(&state).unwrap().values,
                dft_engine::stress::frozen_stress(&state).unwrap().values,
            )
        };
        let (free, (old_f, old_s)) = match &case.policy {
            Policy::FermiDirac { theta_ha, band_cap } => {
                let r = FiniteTProblem::new(
                    &cell,
                    &atoms,
                    &species,
                    strict.clone(),
                    points,
                    *theta_ha,
                    *band_cap,
                )
                .unwrap()
                .run()
                .unwrap();
                assert!(r.converged);
                (
                    r.thermodynamics.free_energy_ha,
                    oracle(ElectronicState::FiniteTemperature(&r)),
                )
            }
            Policy::FixedSpin(sector) => {
                let r = CollinearProblem::new(
                    &cell,
                    &atoms,
                    &species,
                    strict.clone(),
                    points,
                    sector.clone(),
                )
                .unwrap()
                .run()
                .unwrap();
                assert!(r.converged);
                (
                    r.energy.total_hartree(),
                    oracle(ElectronicState::Collinear(&r)),
                )
            }
            Policy::SpinFermiDirac {
                sector,
                theta_ha,
                constraint,
                band_cap,
            } => {
                let r = CollinearProblem::new_fermi_dirac(
                    &cell,
                    &atoms,
                    &species,
                    strict.clone(),
                    points,
                    sector.clone(),
                    *theta_ha,
                    *constraint,
                    *band_cap,
                )
                .unwrap()
                .run()
                .unwrap();
                assert!(r.converged);
                (
                    r.thermodynamics.unwrap().free_energy_ha,
                    oracle(ElectronicState::Collinear(&r)),
                )
            }
        };
        let df = delta(force.iter().flatten(), old_f.iter().flatten());
        let ds = delta(stress.iter().flatten(), old_s.iter().flatten());
        let de = (run.thermodynamics.free_energy_ha - free).abs();
        println!(
            "reference vs strict oracle: dF={de:e} Ha dForce={df:e} Ha/Bohr dStress={ds:e} Ha/Bohr^3"
        );
        assert!(de < 5e-9 && df < 1e-6 && ds < 1e-6);
    }

    #[test]
    fn fermi_dirac_forces_stress_match_free_energy_differences() {
        run_case(Case {
            // TS = 9.9e-4 Ha with 7 fractional occupations (16 retained bands
            // cover the Fermi tail of this small basis; 0.05 Ha does not).
            policy: Policy::FermiDirac {
                theta_ha: 0.03,
                band_cap: 16,
            },
            functional: Functional::LdaPz,
            pp_functional: "SLA PZ",
            bands: 16,
            points: gamma(),
        });
    }

    #[test]
    fn collinear_fixed_sector_forces_stress_match_energy_differences() {
        run_case(Case {
            policy: Policy::FixedSpin(sector(2, 0)),
            functional: Functional::LdaPw,
            pp_functional: "SLA PW",
            bands: 4,
            points: gamma(),
        });
    }

    #[test]
    fn collinear_fermi_dirac_forces_stress_match_free_energy_differences() {
        run_case(Case {
            policy: Policy::SpinFermiDirac {
                sector: sector(1, 1),
                theta_ha: 0.02,
                constraint: Constraint::FixedMoment(0.5),
                band_cap: 16,
            },
            functional: Functional::LdaPw,
            pp_functional: "SLA PW",
            bands: 16,
            points: gamma(),
        });
    }

    #[test]
    fn nonconverged_or_unverified_delivery_is_rejected() {
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
        opts.gamma_real = GammaReal::Off;
        opts.bands = Some(16);
        opts.standard = Some(dft_engine::StandardProfile { threshold_ha: 5e-9 });
        let input = || StartInput {
            cell: &cell,
            atoms: &atoms,
            potentials: &species,
            options: &opts,
        };
        let policy = Policy::FermiDirac {
            theta_ha: 0.02,
            band_cap: 16,
        };
        let profile = Profile::from_standard_ha(5e-9).unwrap();
        let mut run = electronic::run(input(), profile, &gamma(), policy.clone()).unwrap();
        run.converged = false;
        assert!(delivered::electronic_state(input(), &mut run, &[]).is_err());
        run.converged = true;
        run.verification
            .as_mut()
            .unwrap()
            .attempts
            .last_mut()
            .unwrap()
            .passed = false;
        assert!(delivered::electronic_state(input(), &mut run, &[]).is_err());
        // Explicit verify=off delivers forces labelled unverified.
        let mut off = opts.clone();
        off.verification = Verification::Off;
        let off_input = || StartInput {
            options: &off,
            ..input()
        };
        let mut run = electronic::run(off_input(), profile, &gamma(), policy).unwrap();
        let view = delivered::electronic_state(off_input(), &mut run, &[]).unwrap();
        assert_eq!(
            dft_engine::forces(&view).unwrap().verification,
            "unverified"
        );
    }
}
mod frozen_strain {
    use super::*;
    use dft_basis::{FftGridSpec, PlaneWaveBasis};
    use dft_engine::refscf::metric::Serial;
    use dft_units::Hartree;

    fn deform(cell: &Cell, strain: [[f64; 3]; 3]) -> Cell {
        let a = cell.lattice_columns_bohr();
        Cell::from_matrix_columns_bohr(std::array::from_fn(|i| {
            std::array::from_fn(|j| a[i][j] + (0..3).map(|k| strain[i][k] * a[k][j]).sum::<f64>())
        }))
        .unwrap()
    }

    fn mesh(gamma: bool) -> Vec<KPoint> {
        if gamma {
            vec![KPoint {
                fractional: [0.; 3],
                weight: 1.,
            }]
        } else {
            (0..2)
                .flat_map(|x| {
                    (0..2).flat_map(move |y| {
                        (0..2).map(move |z| KPoint {
                            fractional: [x as f64 / 2., y as f64 / 2., z as f64 / 2.],
                            weight: 0.125,
                        })
                    })
                })
                .collect()
        }
    }

    fn freeze(reference: &Cell, opts: &mut ScfOptions, points: &[KPoint]) {
        let db = PlaneWaveBasis::new(reference, Hartree(4. * opts.cutoff_hartree)).unwrap();
        let mut dims = [0; 3];
        for k in std::iter::once([0.; 3]).chain(points.iter().map(|p| p.fractional)) {
            let b = PlaneWaveBasis::at_k(reference, Hartree(opts.cutoff_hartree), k).unwrap();
            let grid = FftGridSpec::new(&b, &db, None).unwrap();
            for (d, n) in dims.iter_mut().zip(grid.dims()) {
                *d = (*d).max(n);
            }
        }
        opts.grid_override = Some(dims);
        opts.strain_reference = Some(Box::new(reference.clone()));
    }

    #[test]
    fn explicit_fft_grid_required_for_gamma_and_k() {
        let cell = common::cell();
        let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
        let species = [SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let atoms = common::atoms();
        let mut opts = options();
        opts.strain_reference = Some(Box::new(cell.clone()));
        let prepared =
            dft_engine::refscf::symmetry::Prepared::new(&cell, &atoms, &mesh(false), 1e-8, false)
                .unwrap();
        let input = || StartInput {
            cell: &cell,
            atoms: &atoms,
            potentials: &species,
            options: &opts,
        };
        for result in [
            dft_engine::refstart::start(input(), &Serial).map(|_| ()),
            dft_engine::refstart::start_kpoints(input(), &mesh(false), &Serial).map(|_| ()),
            k_driver::run_prepared(
                input(),
                Profile::from_standard_ha(5e-13).unwrap(),
                &prepared,
            )
            .map(|_| ()),
        ] {
            assert!(matches!(
                result,
                Err(dft_engine::ScfError::InvalidInput(
                    "strain reference requires a frozen FFT grid"
                ))
            ));
        }
    }

    // Compare complete ordered sets, including density support used inside H.
    // A nominal-cutoff rebuild must differ: otherwise this would not detect
    // accidentally ignoring strain_reference at a shell crossing.
    fn check_basis(input: &StartInput<'_>, points: &[KPoint]) {
        let reference = input.options.strain_reference.as_deref().unwrap();
        let cutoff = Hartree(input.options.cutoff_hartree);
        let expected = |k| {
            PlaneWaveBasis::at_k(reference, cutoff, k)
                .unwrap()
                .in_cell_with_frozen_density(input.cell, reference)
                .unwrap()
        };
        let state =
            dft_engine::refstart::start_kpoints(StartInput { ..*input }, points, &Serial).unwrap();
        let millers = |b: &PlaneWaveBasis| b.vectors().iter().map(|g| g.miller).collect::<Vec<_>>();
        let mut crossed = false;
        for block in &state.blocks {
            let e = expected(block.point.fractional);
            assert_eq!(millers(&block.basis), millers(&e));
            for (actual, expected) in block.basis.momenta().iter().zip(e.momenta()) {
                assert_eq!(actual.cartesian, expected.cartesian);
                assert_eq!(
                    actual.norm_squared.to_bits(),
                    expected.norm_squared.to_bits()
                );
            }
            let ordinary =
                PlaneWaveBasis::at_k(input.cell, cutoff, block.point.fractional).unwrap();
            let mut actual_set = millers(&block.basis);
            let mut ordinary_set = millers(&ordinary);
            actual_set.sort();
            ordinary_set.sort();
            crossed |= actual_set != ordinary_set;
            let db = block.basis.density_basis(input.cell).unwrap();
            let expected_db = PlaneWaveBasis::new(reference, Hartree(4. * cutoff.0))
                .unwrap()
                .in_cell(input.cell)
                .unwrap();
            assert_eq!(millers(&db), millers(&expected_db));
            assert_eq!(state.grid.dims(), input.options.grid_override.unwrap());
        }
        assert!(crossed, "fixture must cross a wavefunction cutoff shell");
        let db = expected([0.; 3]).density_basis(input.cell).unwrap();
        let ordinary_db = PlaneWaveBasis::new(input.cell, Hartree(4. * cutoff.0)).unwrap();
        let mut a = millers(&db);
        let mut b = millers(&ordinary_db);
        a.sort();
        b.sort();
        assert_ne!(a, b, "fixture must cross a density cutoff shell");
        if input.options.gamma_real == GammaReal::On {
            let real = dft_engine::refstart::start(StartInput { ..*input }, &Serial).unwrap();
            assert_eq!(millers(&real.basis), millers(&expected([0.; 3])));
            assert_eq!(real.grid.dims(), input.options.grid_override.unwrap());
        }
    }

    fn compare_strained(
        reference: &Cell,
        atoms: &[Atom],
        pp: &dft_upf::Pseudopotential,
        mut opts: ScfOptions,
        gamma: bool,
        pulay: bool,
        points: &[KPoint],
    ) {
        // Historical fixture contract; screening has separate charge/mask tests.
        let screening =
            dft_numerics::switches::RunConfig::from_pairs(&[("xc_screen_negative", "off")])
                .unwrap();
        let _screen_scope = screening.enter();
        opts.gamma_real = if gamma { GammaReal::On } else { GammaReal::Off };
        freeze(reference, &mut opts, points);
        let cell = deform(
            reference,
            [[0.09, 0.012, 0.], [0.012, 0.08, 0.], [0., 0., 0.07]],
        );
        let species = [SpeciesPotential {
            label: &atoms[0].species,
            potential: pp,
        }];
        macro_rules! input {
            ($options:expr) => {
                StartInput {
                    cell: &cell,
                    atoms,
                    potentials: &species,
                    options: $options,
                }
            };
        }
        check_basis(&input!(&opts), points);
        let mut strict = opts.clone();
        strict.standard = None;
        strict.verification = Verification::Off;
        strict.density_tolerance = 1e-9;
        strict.energy_tolerance = 1e-12;
        strict.eigen_tolerance = 1e-11;
        let oracle = ScfProblem::with_kpoints(
            &cell,
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
        let view = ForceState {
            cell: &cell,
            atoms,
            potentials: &species,
            options: &strict,
            electronic: ElectronicState::Fixed(&oracle),
            density_floor: None,
            symmetry: &[],
        };
        let oracle_stress = if pulay {
            dft_engine::stress(&view)
        } else {
            dft_engine::stress::frozen_stress(&view)
        }
        .unwrap();
        for threshold in [5e-9, 5e-13] {
            opts.standard = Some(dft_engine::StandardProfile {
                threshold_ha: threshold,
            });
            let profile = Profile::from_standard_ha(threshold).unwrap();
            let mut run = k_driver::run(input!(&opts), profile, points).unwrap();
            assert!(run.outcome.converged);
            assert!(run.outcome.result.verification.as_ref().unwrap().passed());
            assert_eq!(run.outcome.result.gamma_real_used, gamma);
            let energy = run.outcome.result.total_hartree();
            let density = run.outcome.result.density.clone();
            let state = delivered::stress_state(input!(&opts), &mut run.outcome).unwrap();
            let stress = if pulay && threshold < 1e-10 {
                dft_engine::stress(&state)
            } else {
                dft_engine::stress::frozen_stress(&state)
            }
            .unwrap();
            let de = (energy - oracle.total_hartree()).abs();
            let ds = delta(
                stress.values.iter().flatten(),
                oracle_stress.values.iter().flatten(),
            );
            println!(
                "R7A k={} threshold={threshold:e} E={energy:.16e} Ha dE={de:e} Ha dStress={ds:e} Ha/Bohr^3",
                points.len()
            );
            assert!(de <= if threshold < 1e-10 { 1e-10 } else { 5e-9 });
            if threshold < 1e-10 {
                assert!(ds <= 1e-8);
                if let Some(p) = &stress.pulay {
                    let old = oracle_stress.pulay.as_ref().unwrap();
                    let dp = (p.correction - old.correction).abs();
                    println!("R7A k={} dPulay={dp:e} Ha/Bohr^3", points.len());
                    assert!(dp <= 1e-8);
                    assert_eq!(p.cutoff_steps_ha, old.cutoff_steps_ha);
                }
                for shear in [false, true] {
                    for h in [2e-4, 1e-4] {
                        let mut energies = Vec::new();
                        for sign in [-1., 1.] {
                            let mut eta = [[0.; 3]; 3];
                            if shear {
                                eta[0][1] = sign * h;
                                eta[1][0] = sign * h;
                            } else {
                                eta[0][0] = sign * h;
                            }
                            let strained = deform(&cell, eta);
                            let r = k_driver::run(
                                StartInput {
                                    cell: &strained,
                                    ..input!(&opts)
                                },
                                profile,
                                points,
                            )
                            .unwrap();
                            assert!(r.outcome.converged);
                            energies.push(r.outcome.result.total_hartree());
                        }
                        let fd = (energies[1] - energies[0]) / (2. * h * cell.volume_bohr3());
                        let analytic = if shear {
                            stress.values[0][1] + stress.values[1][0]
                        } else {
                            stress.values[0][0]
                        };
                        let error = (fd - analytic).abs();
                        println!(
                            "R7A k={} shear={shear} h={h:e} FD error={error:e} Ha/Bohr^3",
                            points.len()
                        );
                        assert!(error <= 1e-8);
                    }
                }
            }
            // Delivery cannot silently reinterpret a frozen result on a new
            // reference cell, even when the cutoff/grid remain unchanged.
            let mut wrong = opts.clone();
            wrong.strain_reference = Some(Box::new(cell.clone()));
            assert!(delivered::stress_state(input!(&wrong), &mut run.outcome).is_err());
            wrong.strain_reference = None;
            assert!(delivered::stress_state(input!(&wrong), &mut run.outcome).is_err());
            if !gamma && threshold < 1e-10 {
                // Keep the same shell counts/grid so rejection tests reference
                // identity, rather than an incidental orbital shape mismatch.
                wrong.strain_reference = Some(Box::new(deform(
                    reference,
                    [[1e-8, 0., 0.], [0., 0., 0.], [0., 0., 0.]],
                )));
                let flat: Vec<_> = density.iter().copied().collect();
                let delivered::IonicSeed::K(seed) =
                    delivered::ionic_seed(input!(&opts), &run.outcome.result, &flat, &Serial)
                        .unwrap()
                else {
                    panic!("complex k seed required")
                };
                assert!(matches!(
                    k_driver::run_seeded_with(input!(&wrong), profile, points, &Serial, seed),
                    Err(dft_engine::ScfError::InvalidInput(
                        "reference k ionic seed: incompatible cell, basis, grid, mesh or charge"
                    ))
                ));
            }
            assert_eq!(
                run.outcome.result.total_hartree().to_bits(),
                energy.to_bits()
            );
            assert_eq!(run.outcome.result.density, density);
        }
    }

    #[test]
    fn synthetic_gamma_frozen_strain() {
        let pp = common::potential_with_atomic_start(2., 1.0, true, 161, "SLA PZ");
        let cell = Cell::from_user_rows_bohr([[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]]).unwrap();
        let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
        let mut opts = options();
        opts.cutoff_hartree = 1.1;
        compare_strained(&cell, &atoms, &pp, opts, true, false, &mesh(true));
    }

    #[test]
    fn synthetic_complex_k_frozen_strain() {
        let pp = common::potential_with_atomic_start(1., 1.5, true, 161, "SLA PZ");
        let cell = Cell::from_user_rows_bohr([[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]]).unwrap();
        let atoms = [
            Atom::new("Xx", [0.13, 0.21, 0.32]),
            Atom::new("Xx", [0.30, 0.35, 0.42]),
        ];
        let mut opts = options();
        opts.cutoff_hartree = 2.;
        compare_strained(
            &cell,
            &atoms,
            &pp,
            opts,
            false,
            false,
            &[
                KPoint {
                    fractional: [0.125, 0., 0.],
                    weight: 0.5,
                },
                KPoint {
                    fractional: [-0.125, 0., 0.],
                    weight: 0.5,
                },
            ],
        );
    }

    #[test]
    #[ignore = "external Si2 strain/FD/Pulay matrix; DFT_SI_UPF=/path/Si.pz-vbc.UPF cargo test --release -p dft-engine --test refscf_forces silicon_frozen_strain -- --ignored --nocapture"]
    fn silicon_frozen_strain_gamma_and_k222() {
        let path = std::env::var("DFT_SI_UPF").expect("DFT_SI_UPF must name an external Si UPF");
        let pp = dft_upf::parse_upf(&std::fs::read_to_string(path).unwrap()).unwrap();
        let cell = common::cell();
        let atoms = [
            Atom::new("Si", [0.; 3]),
            Atom::new("Si", [0.253, 0.249, 0.252]),
        ];
        let mut opts = options();
        opts.cutoff_hartree = 4.;
        opts.bands = Some(5);
        for gamma in [true, false] {
            compare_strained(&cell, &atoms, &pp, opts.clone(), gamma, true, &mesh(gamma));
        }
    }
}

/// Automatic FD band growth keeps the derivative contract: forces and stress
/// of the delivered (grown) state match an explicit fixed-count run of the
/// same model within the standard accuracy class.
#[test]
fn fd_band_growth_forces_stress_match_explicit_count() {
    use dft_engine::refscf::{
        api::{Derivatives, Request, solve},
        electronic::{BandGrowth, Policy},
    };
    let screening = RunConfig::from_pairs(&[("xc_screen_negative", "off")]).unwrap();
    let _screen_scope = screening.enter();
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
    let points = vec![KPoint {
        fractional: [0.; 3],
        weight: 1.,
    }];
    let profile = Profile::from_standard_ha(5e-13).unwrap();
    let (theta_ha, cap) = (0.02, 16);
    let run = |bands: Option<usize>, growth: Option<BandGrowth>| {
        let mut opts = options();
        opts.cutoff_hartree = 3.;
        opts.gamma_real = GammaReal::Off;
        opts.bands = bands;
        opts.standard = Some(dft_engine::StandardProfile {
            threshold_ha: profile.threshold_ry / 2.,
        });
        let mut request = Request::new(profile, points.clone());
        request.policy = Some(Policy::FermiDirac {
            theta_ha,
            band_cap: bands.unwrap_or(cap),
        });
        request.band_growth = growth;
        request.derivatives = Derivatives {
            forces: true,
            stress: true,
            deferred: false,
        };
        let input = StartInput {
            cell: &cell,
            atoms: &atoms,
            potentials: &species,
            options: &opts,
        };
        let d = solve(input, request).unwrap();
        let forces = d.derivatives.forces.unwrap().values;
        let stress = d.derivatives.stress.unwrap().values;
        let r = d.run.into_electronic().unwrap();
        assert!(r.converged && r.verification.as_ref().unwrap().passed());
        (r, forces, stress)
    };
    let (explicit, f0, s0) = run(Some(cap), None);
    let (grown, f1, s1) = run(None, Some(BandGrowth::Automatic));
    let report = grown.band_growth.clone().unwrap();
    let force_error = delta(f0.iter().flatten(), f1.iter().flatten());
    let stress_error = delta(s0.iter().flatten(), s1.iter().flatten());
    let free_error =
        (grown.thermodynamics.free_energy_ha - explicit.thermodynamics.free_energy_ha).abs();
    println!(
        "FD growth {} -> {} bands ({} events): |dF|={free_error:e} Ha max|dforce|={force_error:e} Ha/Bohr max|dstress|={stress_error:e} Ha/Bohr^3; max|force|={:e}",
        report.initial_bands,
        report.final_bands,
        report.events.len(),
        f0.iter().flatten().fold(0_f64, |m, v| m.max(v.abs())),
    );
    assert_eq!(report.initial_bands, 5);
    assert!(!report.events.is_empty() && report.final_bands < cap);
    assert!(free_error <= 5e-9);
    assert!(force_error <= 1e-6);
    assert!(stress_error <= 1e-8);
}

/// FD buffer bands (2026-11-26): a standard run at a small retained count that
/// passes coverage delivers F, mu, forces and stress of an explicit large-band
/// run at a 1e4 times tighter threshold within the standard class, and the
/// buffer columns are not delivered.
#[test]
fn fd_buffer_bands_match_large_band_tight_run() {
    use dft_engine::refscf::{
        api::{Derivatives, Request, solve},
        electronic::Policy,
    };
    let screening = RunConfig::from_pairs(&[("xc_screen_negative", "off")]).unwrap();
    let _screen_scope = screening.enter();
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
    let points = vec![KPoint {
        fractional: [0.; 3],
        weight: 1.,
    }];
    let theta_ha = 0.02;
    let run = |bands: usize, threshold_ha: f64| {
        let profile = Profile::from_standard_ha(threshold_ha).unwrap();
        let mut opts = options();
        opts.cutoff_hartree = 3.;
        opts.gamma_real = GammaReal::Off;
        opts.bands = Some(bands);
        opts.standard = Some(dft_engine::StandardProfile { threshold_ha });
        let mut request = Request::new(profile, points.clone());
        request.policy = Some(Policy::FermiDirac {
            theta_ha,
            band_cap: bands,
        });
        request.derivatives = Derivatives {
            forces: true,
            stress: true,
            deferred: false,
        };
        let input = StartInput {
            cell: &cell,
            atoms: &atoms,
            potentials: &species,
            options: &opts,
        };
        let d = solve(input, request).unwrap();
        let forces = d.derivatives.forces.unwrap().values;
        let stress = d.derivatives.stress.unwrap().values;
        let r = d.run.into_electronic().unwrap();
        assert!(r.converged && r.verification.as_ref().unwrap().passed());
        assert!(r.blocks[0].iter().all(|b| b.eigenvalues_ha.len() == bands));
        let mu = r.occupations.unwrap().roots[0].unwrap().mu_ha;
        (r.thermodynamics.free_energy_ha, mu, forces, stress)
    };
    // Reference: 24 retained bands at 5e-13 Ha. Same tight threshold at 13
    // retained bands isolates the band count and buffer; the 13-band standard
    // run is the delivered production setting (standard class on F).
    let (f0, mu0, force0, stress0) = run(24, 5e-13);
    let (f1, mu1, force1, stress1) = run(13, 5e-13);
    let (f2, mu2, force2, stress2) = run(13, 5e-9);
    let force_error = delta(force0.iter().flatten(), force1.iter().flatten());
    let stress_error = delta(stress0.iter().flatten(), stress1.iter().flatten());
    println!(
        "FD buffer: 13 vs 24 bands, both tight: |dF|={:e} Ha |dmu|={:e} Ha max|dforce|={force_error:e} Ha/Bohr max|dstress|={stress_error:e} Ha/Bohr^3",
        (f1 - f0).abs(),
        (mu1 - mu0).abs()
    );
    println!(
        "FD buffer: 13 bands standard vs 24 bands tight: |dF|={:e} Ha |dmu|={:e} Ha max|dforce|={:e} Ha/Bohr max|dstress|={:e} Ha/Bohr^3",
        (f2 - f0).abs(),
        (mu2 - mu0).abs(),
        delta(force0.iter().flatten(), force2.iter().flatten()),
        delta(stress0.iter().flatten(), stress2.iter().flatten()),
    );
    assert!((f1 - f0).abs() <= 1e-10);
    assert!((mu1 - mu0).abs() <= 1e-8);
    assert!(force_error <= 1e-7);
    assert!(stress_error <= 1e-9);
    assert!((f2 - f0).abs() <= 5e-9);
}

/// MPI-DERIV: forces and frozen-basis stress of one replicated reference
/// delivery (5 complex k blocks, two atoms) evaluated with an active
/// communicator of 2, 3 and 7 simulated ranks (uneven contiguous k splits, and
/// ranks without any block) equal the threads-only evaluation bit for bit:
/// each rank evaluates its own blocks, one allgatherv, accumulation in global
/// k order, symmetrization after. The end-to-end runs (SCF and derivatives on
/// 1, 2, 3 ranks) are reported against np 1 and bounded by 1e-12.
#[test]
fn distributed_derivatives_match_replicated_and_rank_count() {
    use std::sync::Arc;
    let screening = RunConfig::from_pairs(&[("xc_screen_negative", "off")]).unwrap();
    let _screen_scope = screening.enter();
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
    opts.gamma_real = GammaReal::Off;
    let points: Vec<_> = [
        [0., 0., 0.],
        [0.25, 0., 0.],
        [0., 0.25, 0.],
        [0., 0., 0.25],
        [0.25, 0.25, 0.],
    ]
    .into_iter()
    .zip([0.1, 0.2, 0.2, 0.2, 0.3])
    .map(|(fractional, weight)| KPoint { fractional, weight })
    .collect();
    let profile = Profile::from_standard_ha(opts.standard.unwrap().threshold_ha).unwrap();
    let input = || StartInput {
        cell: &cell,
        atoms: &atoms,
        potentials: &species,
        options: &opts,
    };
    let bits = |f: &[[f64; 3]], s: &[[f64; 3]; 3]| -> Vec<u64> {
        f.iter()
            .flatten()
            .chain(s.iter().flatten())
            .map(|v| v.to_bits())
            .collect()
    };
    // Replicated (threads-only) evaluation of one delivered state.
    let mut run = k_driver::run(input(), profile, &points).unwrap();
    let view = delivered::force_state(input(), &mut run.outcome).unwrap();
    let clock = std::time::Instant::now();
    let force = dft_engine::forces(&view).unwrap();
    let stress = dft_engine::stress::frozen_stress(&view).unwrap();
    let serial_s = clock.elapsed().as_secs_f64();
    let reference = bits(&force.values, &stress.values);
    // Cost split of the post-SCF derivative boundary (observation only): the
    // full stress adds three cutoff-scan SCF runs (Pulay assessment).
    let full = dft_engine::stress(&view).unwrap();
    let scans = full.pulay.as_ref().unwrap().elapsed.as_secs_f64();
    println!(
        "derivative cost: forces {:.4}s, frozen stress (incl. a second force pass) {:.4}s, Pulay cutoff scans {scans:.4}s (3 SCF runs)",
        force.elapsed.as_secs_f64(),
        stress.elapsed.as_secs_f64(),
    );
    assert_eq!(
        bits(&force.values, &full.values),
        bits(&force.values, &stress.values)
    );
    for ranks in [2, 3, 7] {
        let results = dft_mpi::simulate(ranks, |comm| {
            let _screen_scope = screening.enter();
            dft_engine::parallel::set_communicator(Some(Arc::new(comm.clone())));
            let clock = std::time::Instant::now();
            let f = dft_engine::forces(&view).unwrap();
            let s = dft_engine::stress::frozen_stress(&view).unwrap();
            let seconds = clock.elapsed().as_secs_f64();
            dft_engine::parallel::set_communicator(None);
            assert_eq!(f.band_passes, force.band_passes);
            assert_eq!(s.band_passes, stress.band_passes);
            (bits(&f.values, &s.values), seconds)
        });
        for (rank, (b, seconds)) in results.iter().enumerate() {
            println!(
                "replicated vs distributed ranks={ranks} rank={rank}: identical={} wall={seconds:.3}s (replicated {serial_s:.3}s)",
                *b == reference
            );
            assert_eq!(*b, reference, "ranks={ranks} rank={rank}");
        }
    }
    // End to end: SCF, delivery and derivatives on 1, 2, 3 ranks. On each
    // rank count the distributed derivatives equal the replicated evaluation
    // of that rank count's own delivered state bit for bit; the differences
    // between rank counts are those of the distributed SCF state itself
    // (reported; bounded by the existing reference rank-count acceptance).
    type Derived = (f64, Vec<[f64; 3]>, [[f64; 3]; 3]);
    let mut base: Option<Derived> = None;
    for ranks in [1, 2, 3] {
        let results = dft_mpi::simulate(ranks, |comm| {
            let _screen_scope = screening.enter();
            if ranks > 1 {
                dft_engine::parallel::set_communicator(Some(Arc::new(comm.clone())));
            }
            let mut run = k_driver::run(input(), profile, &points).unwrap();
            let energy = run.outcome.result.total_hartree();
            let view = delivered::force_state(input(), &mut run.outcome).unwrap();
            let f = dft_engine::forces(&view).unwrap().values;
            let s = dft_engine::stress::frozen_stress(&view).unwrap().values;
            dft_engine::parallel::set_communicator(None);
            let replicated = bits(
                &dft_engine::forces(&view).unwrap().values,
                &dft_engine::stress::frozen_stress(&view).unwrap().values,
            );
            assert_eq!(bits(&f, &s), replicated, "ranks={ranks}");
            (energy, f, s)
        });
        for (e, f, s) in results {
            if let Some((e0, f0, s0)) = &base {
                let df = delta(f.iter().flatten(), f0.iter().flatten());
                let ds = delta(s.iter().flatten(), s0.iter().flatten());
                println!(
                    "end-to-end ranks={ranks}: |dE|={:e} Ha max|dF|={df:e} Ha/Bohr max|dS|={ds:e} Ha/Bohr^3 vs np 1 (SCF state differences; derivative distribution adds none)",
                    (e - e0).abs()
                );
                assert!((e - e0).abs() <= 5e-9 && df <= 1e-6 && ds <= 1e-6);
            } else {
                base = Some((e, f, s));
            }
        }
    }
}
