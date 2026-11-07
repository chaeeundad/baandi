// Included as a unit module so the forced-rebuild oracle has no production flag.
use super::*;

#[test]
fn real_state_reuse_matches_forced_rebuild_each_scf_map() {
    let potential = pool_fixture::potential(2.0, 1.0, true, 101);
    let cell = Cell::from_user_rows_bohr([[4., 0., 0.], [0.2, 4., 0.], [0., 0.1, 5.]]).unwrap();
    for threads in [1, 4] {
        for standard in [false, true] {
            let options = ScfOptions {
                gamma_real: crate::GammaReal::On,
                threads,
                standard: standard.then(StandardProfile::default),
                cutoff_hartree: 2.5,
                bands: Some(3),
                ..Default::default()
            };
            let make = || {
                ScfProblem::new(
                    &cell,
                    &[Atom::new("Xx", [0.13, 0.21, 0.32])],
                    &[SpeciesPotential {
                        label: "Xx",
                        potential: &potential,
                    }],
                    options.clone(),
                )
                .unwrap()
            };
            let mut reused = make();
            let mut fresh = make();
            let mut density = reused.uniform_density();
            let mut guess = None;
            for iteration in 0..6 {
                fresh.states[0].real_solve = None;
                // A supported solver option change at an iteration boundary.
                if iteration == 3 {
                    reused.options.seed += 1;
                    fresh.options.seed += 1;
                    reused.options.davidson_subspace += 2;
                    fresh.options.davidson_subspace += 2;
                }
                let a = reused.step(&density, guess.as_ref(), 1e-9).unwrap();
                let b = fresh.step(&density, guess.as_ref(), 1e-9).unwrap();
                assert_eq!(format!("{:?}", a.solution), format!("{:?}", b.solution));
                assert_eq!(format!("{:?}", a.energy), format!("{:?}", b.energy));
                assert_eq!(
                    a.density.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                    b.density.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
                );
                assert!(reused.states[0].real_solve.is_some());
                density = 0.7 * density + 0.3 * a.density;
                guess = Some(a.solution.eigenvectors);
            }
        }
    }
}

#[test]
#[ignore = "full SCF reuse/rebuild matrix; run release --lib real_state_full_scf -- --ignored"]
fn real_state_full_scf_matches_forced_rebuild() {
    let potential = pool_fixture::potential(2.0, 1.0, true, 101);
    for rows in [
        [[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]],
        [[5., 0., 0.], [0.5, 4., 0.], [0.2, 0.4, 6.]],
    ] {
        let cell = Cell::from_user_rows_bohr(rows).unwrap();
        for threads in [1, 4] {
            for standard in [false, true] {
                let make = || {
                    ScfProblem::new(
                        &cell,
                        &[Atom::new("Xx", [0.13, 0.21, 0.32])],
                        &[SpeciesPotential {
                            label: "Xx",
                            potential: &potential,
                        }],
                        ScfOptions {
                            gamma_real: crate::GammaReal::On,
                            threads,
                            standard: standard.then(StandardProfile::default),
                            cutoff_hartree: 2.5,
                            bands: Some(5),
                            ..Default::default()
                        },
                    )
                    .unwrap()
                };
                let a = make().run().unwrap();
                let mut rebuilt = make();
                rebuilt.force_real_rebuild = true;
                let b = rebuilt.run().unwrap();
                assert_eq!(a.total_hartree().to_bits(), b.total_hartree().to_bits());
                assert_eq!(
                    a.density.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                    b.density.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
                );
                assert_eq!(format!("{:?}", a.k_points), format!("{:?}", b.k_points));
                assert_eq!(a.history.len(), b.history.len());
                for (a, b) in a.history.iter().zip(&b.history) {
                    assert_eq!(a.davidson_iterations, b.davidson_iterations);
                    assert_eq!(a.operator_applications, b.operator_applications);
                }
            }
        }
    }
}
