#[allow(dead_code)]
mod common;
use common::potential_with_functional;
use dft_engine::{collinear::*, *};
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use dft_physics::spin::SpinDensity;
fn run_case(
    functional: Functional,
    rows: [[f64; 3]; 3],
    counts: [usize; 2],
    nonlocal: bool,
    points: &[KPoint],
) -> CollinearResult {
    let cell = Cell::from_user_rows_bohr(rows).unwrap();
    let label = if functional == Functional::Pbe {
        "PBE"
    } else {
        "SLA PW"
    };
    let pp = potential_with_functional((counts[0] + counts[1]) as f64, 1.0, nonlocal, 401, label);
    let opts = options(functional);
    CollinearProblem::new(
        &cell,
        &[Atom::new("Xx", [0.0; 3])],
        &[SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }],
        opts,
        points,
        SpinSector {
            n_up: counts[0],
            n_down: counts[1],
            axis_cartesian: [0.0, 0.0, 1.0],
            initial_moment_seed: None,
        },
    )
    .unwrap()
    .run()
    .unwrap()
}
fn options(functional: Functional) -> ScfOptions {
    ScfOptions {
        functional,
        cutoff_hartree: 2.5,
        bands: Some(4),
        eigen_tolerance: 1e-10,
        density_tolerance: 1e-10,
        energy_tolerance: 1e-11,
        max_iterations: 200,
        // Spin runs use the complex path; the unpolarized reference must too.
        gamma_real: dft_engine::GammaReal::Off,
        ..Default::default()
    }
}
fn gamma() -> Vec<KPoint> {
    vec![KPoint {
        fractional: [0.0; 3],
        weight: 1.0,
    }]
}
fn check(r: &CollinearResult, volume: f64) {
    assert!(
        r.converged,
        "iterations={}, last={:?}",
        r.history.len(),
        r.history.last()
    );
    let c = r.certificate.as_ref().unwrap();
    assert!(c.passed);
    let SpinDensity::Collinear { channels, .. } = &r.density else {
        panic!()
    };
    for (s, n) in channels.iter().enumerate() {
        let count = n.sum() * volume / n.len() as f64;
        let target = if s == 0 {
            r.sector.n_up
        } else {
            r.sector.n_down
        } as f64;
        assert!((count - target).abs() < 1e-10);
        let from_bands: f64 = r.k_points[s]
            .iter()
            .map(|k| k.point.weight * r.occupations[s].iter().sum::<f64>())
            .sum();
        assert!((count - from_bands).abs() < 1e-10);
    }
    assert!((r.moment_mu_b - r.sector.moment()).abs() < 1e-10);
    assert!((r.energy.total_hartree() - r.energy.band_form).abs() < 1e-9);
    assert!(
        (r.energy.total_hartree() - r.energy.harris_foulkes).abs()
            <= r.energy.functional_difference_bound + 1e-11
    );
    println!(
        "E={:.16e} Ha, M={:.12}, iterations={}, certificate charge={:.3e}, spin={:.3e}, gaps={:?}",
        r.energy.total_hartree(),
        r.moment_mu_b,
        r.history.len(),
        c.charge_density_residual,
        c.spin_density_residual,
        c.gap_ha
    );
}
#[test]
fn odd_electron_certificate_and_swap() {
    let rows = [[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]];
    for functional in [Functional::LdaPw, Functional::Pbe] {
        let a = run_case(functional, rows, [1, 0], false, &gamma());
        let b = run_case(functional, rows, [0, 1], false, &gamma());
        check(&a, 125.0);
        check(&b, 125.0);
        assert!((a.energy.total_hartree() - b.energy.total_hartree()).abs() < 1e-12);
        let (
            SpinDensity::Collinear { channels: aa, .. },
            SpinDensity::Collinear { channels: bb, .. },
        ) = (&a.density, &b.density)
        else {
            panic!()
        };
        assert_eq!(aa[0], bb[1]);
        assert_eq!(aa[1], bb[0]);
        assert_eq!(a.certificate.as_ref().unwrap().gap_ha[1], None);
        assert!(a.legacy_scalar_result().is_err());
    }
}
#[test]
fn zero_moment_reduces_to_same_unpolarized_functional() {
    let rows = [[5.0, 0.0, 0.0], [0.5, 5.0, 0.0], [0.3, 0.2, 5.0]];
    let cell = Cell::from_user_rows_bohr(rows).unwrap();
    for functional in [Functional::LdaPw, Functional::Pbe] {
        let label = if functional == Functional::Pbe {
            "PBE"
        } else {
            "SLA PW"
        };
        let pp = potential_with_functional(2.0, 1.0, true, 401, label);
        let u = ScfProblem::new(
            &cell,
            &[Atom::new("Xx", [0.0; 3])],
            &[SpeciesPotential {
                label: "Xx",
                potential: &pp,
            }],
            options(functional),
        )
        .unwrap()
        .run()
        .unwrap();
        let r = run_case(functional, rows, [1, 1], true, &gamma());
        check(&r, cell.volume_bohr3());
        let SpinDensity::Collinear { channels, .. } = &r.density else {
            panic!()
        };
        let total = &channels[0] + &channels[1];
        let de = (r.energy.total_hartree() - u.total_hartree()).abs();
        let dn = density_distance(&total, &u.density, cell.volume_bohr3(), 2.0);
        println!(
            "zero moment {functional:?}: dE={de:e} Ha, dn={dn:e}, iterations spin={} unpolarized={}",
            r.history.len(),
            u.history.len()
        );
        assert!(de <= 1e-10);
        assert!(dn <= 1e-9);
        // Frozen reduction rule: at most one iteration for independent rounding
        // and the additional spin/strict-empty-band convergence gate.
        assert!(r.history.len().abs_diff(u.history.len()) <= 1);
    }
}
#[test]
fn sector_validation() {
    let rows = [[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]];
    let cell = Cell::from_user_rows_bohr(rows).unwrap();
    let pp = potential_with_functional(1.0, 1.0, false, 101, "SLA PZ");
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let atoms = [Atom::new("Xx", [0.0; 3])];
    let sector = SpinSector {
        n_up: 1,
        n_down: 0,
        axis_cartesian: [0.0, 0.0, 1.0],
        initial_moment_seed: None,
    };
    let error = CollinearProblem::new(
        &cell,
        &atoms,
        &species,
        ScfOptions::default(),
        &gamma(),
        sector.clone(),
    )
    .err()
    .unwrap()
    .to_string();
    assert!(error.contains("required polarized capability unavailable for lda-pz"));
    assert!(
        SpinSector {
            n_down: 1,
            ..sector.clone()
        }
        .validate(1.0)
        .is_err()
    );
    assert!(
        SpinSector {
            axis_cartesian: [0.0; 3],
            ..sector
        }
        .validate(1.0)
        .is_err()
    );
}

/// Complete geometry/mesh coverage exceeds the debug gate budget.
/// cargo test --release -p dft-engine --test collinear mf03_synthetic_matrix -- --ignored --nocapture
#[test]
#[ignore = "release-only MF03 geometry, mesh, projector and band matrix"]
fn mf03_synthetic_matrix() {
    let geometries = [
        [[8.0, 0.0, 0.0], [0.0, 8.0, 0.0], [0.0, 0.0, 8.0]],
        [[8.0, 0.0, 0.0], [1.0, 7.0, 0.0], [0.5, 1.0, 9.0]],
        [[6.0, 0.0, 0.0], [0.0, 6.0, 0.0], [0.0, 0.0, 18.0]],
        [[8.0, 0.0, 0.0], [0.0, 8.0, 0.0], [0.0, 0.0, 24.0]],
    ];
    let quadratures = [
        gamma(),
        dft_geometry::kpoints::monkhorst_pack([2, 2, 2], [0.5; 3]).unwrap(),
        vec![
            KPoint {
                fractional: [0.0; 3],
                weight: 0.2,
            },
            KPoint {
                fractional: [0.13, 0.07, 0.0],
                weight: 0.3,
            },
            KPoint {
                fractional: [-0.17, 0.03, 0.11],
                weight: 0.5,
            },
        ],
    ];
    for (g, rows) in geometries.into_iter().enumerate() {
        let cell = Cell::from_user_rows_bohr(rows).unwrap();
        let atoms = if g < 2 {
            vec![
                Atom::new("Xx", [0.0; 3]),
                Atom::new("Xx", [0.27, 0.23, 0.19]),
            ]
        } else {
            vec![
                Atom::new("Xx", [0.0; 3]),
                Atom::new("Xx", [0.27, 0.23, 0.19]),
                Atom::new("Xx", [0.5, 0.5, 0.5]),
                Atom::new("Xx", [0.77, 0.73, 0.69]),
            ]
        };
        for (mesh, points) in quadratures.iter().enumerate() {
            for functional in [Functional::LdaPw, Functional::Pbe] {
                let charge = if g % 2 == 0 { 1.0 } else { 2.0 };
                let occupied = atoms.len() * (charge as usize) / 2;
                let label = if functional == Functional::Pbe {
                    "PBE"
                } else {
                    "SLA PW"
                };
                let pp = potential_with_functional(charge, 1.0, g % 2 != 0, 401, label);
                let species = [SpeciesPotential {
                    label: "Xx",
                    potential: &pp,
                }];
                let mut opts = options(functional);
                opts.bands = Some(if g < 2 { 8 } else { 16 });
                let r = CollinearProblem::new(
                    &cell,
                    &atoms,
                    &species,
                    opts.clone(),
                    points,
                    SpinSector {
                        n_up: occupied,
                        n_down: occupied,
                        axis_cartesian: [0.0, 1.0, 0.0],
                        initial_moment_seed: Some(0.0),
                    },
                )
                .unwrap()
                .run();
                println!("matrix geometry={g} mesh={mesh} xc={functional:?}");
                let u = ScfProblem::with_kpoints(
                    &cell,
                    &atoms,
                    &species,
                    opts,
                    points,
                    OccupationPolicy::Insulator {
                        minimum_gap_hartree: 1e-7,
                    },
                )
                .unwrap()
                .run();
                let (r, u) = match (r, u) {
                    (Ok(r), Ok(u)) => {
                        assert!(
                            !(g == 1 && mesh == 2 && functional == Functional::LdaPw),
                            "known noninsulating fixture was certified"
                        );
                        (r, u)
                    }
                    (
                        Err(ScfError::NonInsulating { gap_hartree: a, .. }),
                        Err(ScfError::NonInsulating { gap_hartree: b, .. }),
                    ) => {
                        assert!(
                            g == 1 && mesh == 2 && functional == Functional::LdaPw,
                            "unexpected loss of insulating gap"
                        );
                        println!(
                            "matched noninsulating rejection spin_gap={a:e} Ha legacy_gap={b:e} Ha"
                        );
                        assert!((a - b).abs() < 1e-10);
                        continue;
                    }
                    (a, b) => panic!(
                        "inconsistent outcomes: spin error={:?}, legacy error={:?}",
                        a.err(),
                        b.err()
                    ),
                };
                check(&r, cell.volume_bohr3());
                let SpinDensity::Collinear { channels, .. } = &r.density else {
                    panic!()
                };
                let de = (r.energy.total_hartree() - u.total_hartree()).abs();
                let dn = density_distance(
                    &(&channels[0] + &channels[1]),
                    &u.density,
                    cell.volume_bohr3(),
                    (2 * occupied) as f64,
                );
                println!(
                    "matrix reduction dE={de:e} Ha dn={dn:e} iterations={}/{}",
                    r.history.len(),
                    u.history.len()
                );
                assert!(de <= 1e-10 && dn <= 1e-9);
                assert!(r.history.len().abs_diff(u.history.len()) <= 1);
            }
        }
    }
}

#[test]
fn occupied_channels_and_existing_mixers() {
    let cell =
        Cell::from_user_rows_bohr([[5.0, 0.0, 0.0], [0.7, 5.4, 0.0], [0.3, 0.5, 6.0]]).unwrap();
    let pp = potential_with_functional(3.0, 1.0, true, 401, "SLA PW");
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let atoms = [Atom::new("Xx", [0.0; 3])];
    let mut reference: Option<f64> = None;
    for mixing in [
        Mixing::Pulay {
            beta: 0.3,
            history: 6,
        },
        Mixing::Linear { beta: 0.3 },
        Mixing::KerkerPulay {
            beta: 0.3,
            history: 6,
            q0_inverse_bohr: 1.0,
        },
    ] {
        let mut opts = options(Functional::LdaPw);
        opts.mixing = mixing;
        let r = CollinearProblem::new(
            &cell,
            &atoms,
            &species,
            opts,
            &gamma(),
            SpinSector {
                n_up: 2,
                n_down: 1,
                axis_cartesian: [1.0, 0.0, 0.0],
                initial_moment_seed: Some(0.0),
            },
        )
        .unwrap()
        .run()
        .unwrap();
        check(&r, cell.volume_bohr3());
        let e = r.energy.total_hartree();
        if let Some(old) = reference {
            assert!((e - old).abs() < 1e-10);
        } else {
            reference = Some(e);
        }
        assert!(
            r.certificate
                .as_ref()
                .unwrap()
                .gap_ha
                .iter()
                .all(Option::is_some)
        );
    }
}
