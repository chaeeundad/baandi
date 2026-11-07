//! Release-only reproduction of the existing odd-electron spin and synthetic
//! Fermi-Dirac metal fixtures across all P36B budgets.
#[allow(dead_code)]
mod common;
use dft_engine::{collinear::*, finite_temperature::*, *};
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use dft_physics::spin::SpinDensity;

fn snapshot(
    energy: f64,
    density: impl Iterator<Item = f64>,
    points: &[KPointSolution],
) -> Vec<u64> {
    let mut bits = vec![energy.to_bits()];
    bits.extend(density.map(f64::to_bits));
    for p in points {
        bits.extend([
            p.solution.iterations as u64,
            p.solution.operator_applications as u64,
        ]);
        bits.extend(p.solution.residuals.iter().map(|v| v.to_bits()));
        bits.extend(p.solution.eigenvalues.iter().map(|v| v.to_bits()));
        for j in 0..p.solution.eigenvectors.cols() {
            for v in p.solution.eigenvectors.column(j) {
                bits.extend([v.re.to_bits(), v.im.to_bits()]);
            }
        }
    }
    bits
}

#[test]
#[ignore = "release matrix: cargo test --release -p dft-engine --test hybrid_budget -- --ignored --nocapture"]
fn spin_and_fermi_dirac_budget_identity() {
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0., 5., 0.], [0., 0., 5.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.; 3])];
    let points = [KPoint {
        fractional: [0.; 3],
        weight: 1.,
    }];
    for spin in [true, false] {
        let pp = common::potential_with_functional(
            1.,
            1.,
            false,
            401,
            if spin { "SLA PW" } else { "SLA PZ" },
        );
        let refs = [SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let mut reference = None;
        for threads in [1, 2, 3, 4, 6, 8, 10, 12, 16] {
            let options = ScfOptions {
                threads,
                functional: if spin {
                    Functional::LdaPw
                } else {
                    Functional::LdaPz
                },
                cutoff_hartree: 2.5,
                bands: Some(if spin { 4 } else { 12 }),
                eigen_tolerance: 1e-10,
                density_tolerance: 1e-10,
                energy_tolerance: 1e-11,
                max_iterations: 200,
                ..Default::default()
            };
            let (mut bits, work) = if spin {
                let r = CollinearProblem::new(
                    &cell,
                    &atoms,
                    &refs,
                    options,
                    &points,
                    SpinSector {
                        n_up: 1,
                        n_down: 0,
                        axis_cartesian: [0., 0., 1.],
                        initial_moment_seed: None,
                    },
                )
                .unwrap()
                .run()
                .unwrap();
                assert!(r.converged);
                let SpinDensity::Collinear { channels, .. } = &r.density else {
                    panic!()
                };
                let mut bits = snapshot(
                    r.energy.total_hartree(),
                    channels.iter().flat_map(|n| n.iter().copied()),
                    &r.k_points[0],
                );
                bits.extend(snapshot(0., std::iter::empty(), &r.k_points[1]));
                (bits, format!("{:?} {:?}", r.history, r.certificate))
            } else {
                let r = FiniteTProblem::new(&cell, &atoms, &refs, options, &points, 0.01, 12)
                    .unwrap()
                    .run()
                    .unwrap();
                assert!(r.converged);
                (
                    snapshot(
                        r.thermodynamics.free_energy_ha,
                        r.state.density.iter().copied(),
                        &r.state.k_points,
                    ),
                    format!(
                        "{} {:?} {:?}",
                        r.iterations, r.certificate, r.thermodynamics
                    ),
                )
            };
            let records = work_accounting::take();
            if work_accounting::ENABLED {
                let counts = [
                    records.len() as u64,
                    records.iter().map(|r| r.iterations as u64).sum(),
                    records.iter().map(|r| r.vector_applications as u64).sum(),
                ];
                bits.extend(counts);
                println!("accounting spin={spin} T={threads} solves/Davidson/H={counts:?}");
            }
            println!(
                "spin={spin} T={threads} energy_bits={} values={}",
                bits[0],
                bits.len()
            );
            if let Some(dir) = std::env::var_os("P36B_SNAPSHOT_DIR") {
                let dir = std::path::PathBuf::from(dir);
                std::fs::create_dir_all(&dir).unwrap();
                let mut bytes: Vec<u8> = bits.iter().flat_map(|v| v.to_le_bytes()).collect();
                bytes.extend(work.as_bytes());
                std::fs::write(dir.join(format!("spin-{spin}-T{threads}.bin")), bytes).unwrap();
            }
            if let Some((ref b, ref w)) = reference {
                assert_eq!(&bits, b, "spin={spin} T={threads}");
                assert_eq!(&work, w, "spin={spin} T={threads}");
            } else {
                reference = Some((bits, work));
            }
        }
    }
}
