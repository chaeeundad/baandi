//! R2 dielectric response driver (`docs/roadmap/rpa-crpa-plan.md`, ticket
//! R2) on synthetic in-memory pseudopotentials with nonlocal projectors and
//! a fixed nonuniform density (not an SCF fixed point and not QE reference
//! evidence): NSCF bands -> chi0(q -> 0, i w) with head and wings -> eps_M
//! tensors; MPI level 1 (np = 1 bit-identical to the serial adapter, other
//! rank counts within roundoff); insulator, spin and k-set guards.
//!
//! Fast gate (debug): `cargo test -p dft-engine --test response_r2`
#[allow(dead_code)]
mod common;
use dft_engine::refscf::nscf::BandSolver;
use dft_engine::response::{DielectricOptions, DielectricResult, dielectric_with};
use dft_engine::{refstart::StartInput, *};
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use dft_physics::data::Array3;

struct System<'a> {
    cell: Cell,
    atoms: Vec<Atom>,
    potentials: Vec<SpeciesPotential<'a>>,
    options: ScfOptions,
    density: Array3<f64>,
    points: Vec<KPoint>,
}
impl System<'_> {
    fn input(&self) -> StartInput<'_> {
        StartInput {
            cell: &self.cell,
            atoms: &self.atoms,
            potentials: &self.potentials,
            options: &self.options,
        }
    }
}

fn gamma_mesh(n: [usize; 3]) -> Vec<KPoint> {
    let total = (n[0] * n[1] * n[2]) as f64;
    let mut points = Vec::new();
    for i in 0..n[0] {
        for j in 0..n[1] {
            for k in 0..n[2] {
                points.push(KPoint {
                    fractional: [
                        i as f64 / n[0] as f64,
                        j as f64 / n[1] as f64,
                        k as f64 / n[2] as f64,
                    ],
                    weight: 1. / total,
                });
            }
        }
    }
    points
}

/// Charge-normalized smooth density on the problem's grid (atomic
/// superposition from the synthetic PP_RHOATOM).
fn density(system: &System<'_>) -> Array3<f64> {
    let problem = ScfProblem::with_kpoints(
        &system.cell,
        &system.atoms,
        &system.potentials,
        system.options.clone(),
        &system.points,
        OccupationPolicy::Insulator {
            minimum_gap_hartree: 1e-7,
        },
    )
    .unwrap();
    problem.starting_density()
}

fn bits(r: &DielectricResult) -> Vec<u64> {
    r.epsilon
        .iter()
        .flat_map(|e| {
            e.with_local_fields
                .iter()
                .chain(&e.without_local_fields)
                .flatten()
                .map(|v| v.to_bits())
                .collect::<Vec<_>>()
        })
        .chain(
            r.head
                .iter()
                .flatten()
                .flatten()
                .flat_map(|z| [z.re.to_bits(), z.im.to_bits()]),
        )
        .collect()
}

fn with_system(test: impl FnOnce(&System<'_>)) {
    let p = common::potential_with_atomic_start_dij(
        2.,
        0.8,
        true,
        401,
        "SLA PZ",
        "-1 0.1 0 0.1 -0.5 0 0 0 0.8",
    );
    let atoms = vec![Atom::new("Xx", [0.13, 0.21, 0.37])];
    let potentials = vec![SpeciesPotential {
        label: "Xx",
        potential: &p,
    }];
    let options = ScfOptions {
        cutoff_hartree: 4.0,
        bands: Some(8),
        eigen_tolerance: 1e-8,
        functional: Functional::LdaPz,
        ..ScfOptions::strict_reference()
    };
    let mut system = System {
        cell: Cell::from_user_rows_bohr([[6.0, 0., 0.], [0.7, 5.5, 0.], [0.3, 0.4, 6.5]]).unwrap(),
        atoms,
        potentials,
        options,
        density: Array3::zeros((1, 1, 1)),
        points: gamma_mesh([2, 2, 2]),
    };
    system.density = density(&system);
    test(&system);
}

fn options() -> DielectricOptions {
    DielectricOptions {
        bands: Some(24),
        frequency_points: 4,
        band_block: 7,
        ..DielectricOptions::new(3.0)
    }
}

#[test]
fn dielectric_tensor_and_rank_contract() {
    with_system(|system| {
        let serial = dielectric_with(
            system.input(),
            std::slice::from_ref(&system.density),
            &system.points,
            options(),
            &refscf::metric::Serial,
        )
        .unwrap();
        let e0 = &serial.epsilon[0];
        println!(
            "synthetic skewed one-atom Xx, 2x2x2: N_G {}, bands {}, occupied {}, gap {:.4} Ha, transitions {}, units {}, omega0 {:.4} Ha",
            serial.response_size,
            serial.bands,
            serial.occupied,
            serial.conduction_minimum - serial.valence_maximum,
            serial.transitions,
            serial.units,
            serial.omega0
        );
        println!(
            "eps_M(0) without LF {:.6?}\neps_M(0) with LF {:.6?}\nmin body eigenvalue {:.6}, asymmetry {:.2e}",
            e0.without_local_fields, e0.with_local_fields, e0.min_body_eigenvalue, e0.max_asymmetry
        );
        assert!(serial.max_residual <= 1e-6);
        assert!(e0.max_asymmetry <= 1e-10);
        assert!(e0.min_body_eigenvalue >= 1. - 1e-12);
        for a in 0..3 {
            assert!(e0.without_local_fields[a][a] > 1.);
            assert!(e0.with_local_fields[a][a] > 1.);
        }
        // eps_M(i w) decreases monotonically toward 1 on the imaginary axis.
        for pair in serial.epsilon.windows(2) {
            assert!(pair[1].with_local_fields[0][0] <= pair[0].with_local_fields[0][0] + 1e-12);
            assert!(pair[1].with_local_fields[0][0] >= 1.);
        }
        let reference = bits(&serial);
        for ranks in [1, 2, 3] {
            let runs = dft_mpi::simulate(ranks, |comm| {
                dielectric_with(
                    system.input(),
                    std::slice::from_ref(&system.density),
                    &system.points,
                    options(),
                    comm,
                )
                .unwrap()
            });
            let first = bits(&runs[0]);
            let mut spread = 0_f64;
            for run in &runs {
                assert_eq!(bits(run), first, "ranks disagree at np={ranks}");
                for (a, b) in run.epsilon.iter().zip(&serial.epsilon) {
                    for (x, y) in a
                        .with_local_fields
                        .iter()
                        .flatten()
                        .chain(a.without_local_fields.iter().flatten())
                        .zip(
                            b.with_local_fields
                                .iter()
                                .flatten()
                                .chain(b.without_local_fields.iter().flatten()),
                        )
                    {
                        spread = spread.max((x - y).abs());
                    }
                }
            }
            if ranks == 1 {
                assert_eq!(
                    first, reference,
                    "np=1 not bit-identical to the serial adapter"
                );
            }
            println!("np={ranks}: max |eps - serial| {spread:.3e}");
            assert!(spread <= 1e-12, "rank spread {spread:e}");
        }
    });
}

/// R2-VEL-CHECK diagnostic dipole record: results unchanged, every (v, c)
/// band pair in exactly one group pair, nonnegative sums, the nonlocal
/// commutator visible, and the r sums consistent with the velocity sums.
#[test]
fn dipole_record_is_diagnostic_only_and_complete() {
    use dft_engine::response::{DIPOLE_GROUP_TOLERANCE, DipoleDump};
    with_system(|system| {
        let run = |dump| {
            dielectric_with(
                system.input(),
                std::slice::from_ref(&system.density),
                &system.points,
                DielectricOptions {
                    dipole_dump: dump,
                    ..options()
                },
                &refscf::metric::Serial,
            )
            .unwrap()
        };
        let off = run(DipoleDump::Off);
        let all = run(DipoleDump::All);
        let first = run(DipoleDump::GammaAndFirst);
        assert!(off.dipoles.is_empty());
        assert_eq!(bits(&off), bits(&all), "the record must not change results");
        assert_eq!(all.dipoles.len(), system.points.len());
        assert_eq!(first.dipoles.len(), 2);
        let mut nonlocal = 0_f64;
        for (k, record) in all.dipoles.iter().enumerate() {
            assert_eq!(record.k_index, k);
            let e = &all.eigenvalues[k];
            let pairs: usize = record
                .pairs
                .iter()
                .map(|p| (p.v.1 - p.v.0 + 1) * (p.c.1 - p.c.0 + 1))
                .sum();
            assert_eq!(pairs, all.occupied * (all.bands - all.occupied));
            for p in &record.pairs {
                assert!(p.v.1 < all.occupied && p.c.0 >= all.occupied && p.c.1 < all.bands);
                for g in [p.v, p.c] {
                    assert!(e[g.1] - e[g.0] < DIPOLE_GROUP_TOLERANCE * (g.1 - g.0 + 1) as f64);
                }
                for a in 0..3 {
                    assert!(p.velocity[a] >= 0. && p.velocity_local[a] >= 0.);
                    // r sums between the bounds set by the extreme gaps.
                    let (lo, hi) = (e[p.c.0] - e[p.v.1], e[p.c.1] - e[p.v.0]);
                    assert!(p.position[a] <= p.velocity[a] / (lo * lo) * (1. + 1e-12));
                    assert!(p.position[a] >= p.velocity[a] / (hi * hi) * (1. - 1e-12));
                    nonlocal = nonlocal.max((p.velocity[a] - p.velocity_local[a]).abs());
                }
            }
        }
        assert!(nonlocal > 1e-6, "the fixture must show the nonlocal term");
    });
}

#[test]
fn guards_reject_unsupported_inputs() {
    with_system(|system| {
        let run = |points: &[KPoint], density: &[Array3<f64>], o: DielectricOptions| {
            dielectric_with(system.input(), density, points, o, &refscf::metric::Serial)
        };
        let d = std::slice::from_ref(&system.density);
        // Not time-reversal symmetric.
        let half = [KPoint {
            fractional: [0.25, 0., 0.],
            weight: 1.,
        }];
        assert!(run(&half, d, options()).is_err());
        // Spin channels.
        let two = [system.density.clone(), system.density.clone()];
        assert!(run(&system.points, &two, options()).is_err());
        // Response sphere beyond the density sphere.
        assert!(run(&system.points, d, DielectricOptions::new(17.0)).is_err());
        // Chebyshev bands agree with the dense default.
        let cheb = run(
            &system.points,
            d,
            DielectricOptions {
                solver: BandSolver::Chebyshev,
                ..options()
            },
        )
        .unwrap();
        let dense = run(&system.points, d, options()).unwrap();
        let diff = (cheb.epsilon[0].with_local_fields[0][0]
            - dense.epsilon[0].with_local_fields[0][0])
            .abs();
        println!("Chebyshev vs dense eps_xx(0) with LF: {diff:.3e}");
        assert!(diff <= 1e-8, "{diff:e}");
    });
}
