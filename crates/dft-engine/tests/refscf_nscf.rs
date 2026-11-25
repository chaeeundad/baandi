//! Fixed-density spectra; synthetic data are not QE reference evidence.
//! Full matrix: cargo test --release -p dft-engine --test refscf_nscf -- --include-ignored --nocapture
//! Real Si2 additionally requires DFT_SI_UPF=/absolute/path/Si.pz-vbc.UPF.
#[allow(dead_code)]
mod common;
use dft_engine::{refscf::nscf::nscf_with, refstart::StartInput, *};
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use dft_physics::data::Array3;

fn path() -> Vec<KPoint> {
    [[0.; 3], [0.5, 0., 0.5], [0.5; 3]]
        .into_iter()
        .map(|fractional| KPoint {
            fractional,
            weight: 1. / 3.,
        })
        .collect()
}
fn compare(
    name: &str,
    cell: &Cell,
    atoms: &[Atom],
    pp: &[SpeciesPotential<'_>],
    options: ScfOptions,
    ranks: &[usize],
) {
    let points = path();
    let problem = ScfProblem::with_kpoints(
        cell,
        atoms,
        pp,
        options.clone(),
        &points,
        OccupationPolicy::Insulator {
            minimum_gap_hartree: 1e-7,
        },
    )
    .unwrap();
    let mut density = problem.starting_density();
    // Use a nonuniform, exactly charge-normalized supplied density, independent
    // of either solver's atomic start. This is intentionally not a fixed point.
    let shape = density.shape().to_vec();
    let scale = problem.electron_count() / cell.volume_bohr3();
    for ((i, j, k), n) in density.indexed_iter_mut() {
        *n = scale
            * (1.
                + 0.1 * (std::f64::consts::TAU * i as f64 / shape[0] as f64).cos()
                + 0.04
                    * (std::f64::consts::TAU
                        * (j as f64 / shape[1] as f64 + k as f64 / shape[2] as f64))
                        .sin());
    }
    let before = density.clone();
    // Legacy NSCF retains SCF buffer-band stopping. Request its full basis
    // at each k so the wrapper's own all-band residual guard certifies the
    // oracle. K weights do not enter a fixed-density Hamiltonian.
    let legacy: Vec<_> = points
        .iter()
        .map(|point| {
            let mut oracle_options = options.clone();
            oracle_options.bands = Some(
                dft_basis::PlaneWaveBasis::at_k(
                    cell,
                    dft_physics::data::Hartree(options.cutoff_hartree),
                    point.fractional,
                )
                .unwrap()
                .vectors()
                .len(),
            );
            dft_engine::nscf(
                cell,
                atoms,
                pp,
                &density,
                &[KPoint {
                    weight: 1.,
                    ..*point
                }],
                oracle_options,
            )
            .unwrap()
            .solutions
            .remove(0)
        })
        .collect();
    assert!(
        legacy
            .iter()
            .all(|s| s.solution.residuals.iter().all(|r| *r <= 1e-8))
    );
    for &ranks in ranks {
        let answers = dft_mpi::simulate(ranks, |comm| {
            nscf_with(
                StartInput {
                    cell,
                    atoms,
                    potentials: pp,
                    options: &options,
                },
                std::slice::from_ref(&density),
                &points,
                options.bands.unwrap(),
                1e-8,
                comm,
            )
            .unwrap()
        });
        for answer in answers {
            assert_eq!(answer.points, points);
            assert_eq!(answer.channels.len(), 1);
            let mut error = 0_f64;
            let mut residual = 0_f64;
            for (new, old) in answer.channels[0].iter().zip(&legacy) {
                assert_eq!(new.eigenvalues.len(), options.bands.unwrap());
                for (a, b) in new.eigenvalues.iter().zip(&old.solution.eigenvalues) {
                    error = error.max((a - b).abs());
                }
                residual = new.residuals.iter().copied().fold(residual, f64::max);
            }
            assert!(
                error <= 1e-8,
                "{name} ranks={ranks}: eigenvalue error {error:e} Ha"
            );
            assert!(residual <= 1e-8);
            println!(
                "{name} ranks={ranks}: max eigenvalue error {error:e} Ha; max residual {residual:e} Ha"
            );
        }
    }
    assert_eq!(density, before);
}
fn options(functional: Functional, bands: usize) -> ScfOptions {
    ScfOptions {
        cutoff_hartree: 1.5,
        bands: Some(bands),
        eigen_tolerance: 1e-8,
        functional,
        ..ScfOptions::strict_reference()
    }
}
#[test]
fn scalar_fixed_density_matches_legacy() {
    let cell = Cell::from_user_rows_bohr([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]]).unwrap();
    let pp = common::potential_with_orbitals(2., 1.5, true, 101, "SLA PZ");
    compare(
        "cubic",
        &cell,
        &[Atom::new("Xx", [0.13, 0.21, 0.37])],
        &[SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }],
        options(Functional::LdaPz, 5),
        &[1],
    );
}
#[test]
fn synthetic_rank_matrix() {
    let cell = Cell::from_user_rows_bohr([[6., 0., 0.], [0.7, 5.5, 0.], [0.3, 0.4, 6.5]]).unwrap();
    let p = common::potential_with_atomic_start(2., 1.5, true, 101, "SLA PZ");
    let q = common::potential_with_atomic_start(2., 1.8, false, 101, "SLA PZ");
    compare(
        "skew-two-species",
        &cell,
        &[
            Atom::new("Xx", [0.13, 0.21, 0.37]),
            Atom::new("Yy", [0.61, 0.52, 0.74]),
        ],
        &[
            SpeciesPotential {
                label: "Xx",
                potential: &p,
            },
            SpeciesPotential {
                label: "Yy",
                potential: &q,
            },
        ],
        options(Functional::LdaPz, 6),
        &[1, 3, 8],
    );
    let pbe = common::potential_with_nlcc(true, "PBE");
    compare(
        "PBE-NLCC",
        &cell,
        &[Atom::new("Xx", [0.13, 0.21, 0.37])],
        &[SpeciesPotential {
            label: "Xx",
            potential: &pbe,
        }],
        options(Functional::Pbe, 5),
        &[1, 3, 8],
    );
}
#[test]
#[ignore = "external Si2 UPF; DFT_SI_UPF=... cargo test --release -p dft-engine --test refscf_nscf silicon_path -- --ignored --nocapture"]
fn silicon_path() {
    let pp = dft_upf::parse_upf(
        &std::fs::read_to_string(std::env::var("DFT_SI_UPF").expect("external Si UPF required"))
            .unwrap(),
    )
    .unwrap();
    let mut o = options(Functional::LdaPz, 8);
    o.cutoff_hartree = 3.;
    compare(
        "Si2-Gamma-X-L",
        &common::cell(),
        &common::atoms(),
        &[SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }],
        o,
        &[1, 3, 8],
    );
}
#[test]
fn rejects_invalid_input_before_solving() {
    let cell = Cell::from_user_rows_bohr([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]]).unwrap();
    let pp = common::potential_with_atomic_start(2., 1.5, true, 101, "SLA PZ");
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let atoms = [Atom::new("Xx", [0.; 3])];
    let options = options(Functional::LdaPz, 5);
    let input = || StartInput {
        cell: &cell,
        atoms: &atoms,
        potentials: &potentials,
        options: &options,
    };
    let good = Array3::from_elem([9; 3], 2. / cell.volume_bohr3());
    for bad in [
        good.mapv(|x| 2. * x),
        good.mapv(|_| f64::NAN),
        good.mapv(|_| -1.),
    ] {
        assert!(nscf_with(input(), &[bad], &path(), 5, 1e-8, &refscf::metric::Serial).is_err());
    }
    for tolerance in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(
            nscf_with(
                input(),
                std::slice::from_ref(&good),
                &path(),
                5,
                tolerance,
                &refscf::metric::Serial
            )
            .is_err()
        );
    }
    assert!(
        nscf_with(
            input(),
            std::slice::from_ref(&good),
            &path(),
            4,
            1e-8,
            &refscf::metric::Serial
        )
        .is_err()
    );
    assert!(
        nscf_with(
            input(),
            std::slice::from_ref(&good),
            &[],
            5,
            1e-8,
            &refscf::metric::Serial
        )
        .is_err()
    );
    assert!(nscf_with(input(), &[], &path(), 5, 1e-8, &refscf::metric::Serial).is_err());
    assert!(
        nscf_with(
            input(),
            &[Array3::zeros([0, 9, 9])],
            &path(),
            5,
            1e-8,
            &refscf::metric::Serial
        )
        .is_err()
    );
}

#[test]
#[ignore = "spin dense-oracle rank matrix; run in release with --include-ignored --nocapture"]
fn spin_channels_match_serial_dense_oracle() {
    use dft_numerics::linalg::hermitian_eigh;
    use dft_physics::data::{Complex64, ComplexMatrix, Hartree};
    use dft_physics::{
        hamiltonian::{Hamiltonian, NonlocalOperator, Species},
        spin::{SpinFunctional, SpinHartree, SpinXc, SpinXcResponse},
    };
    use dft_solvers::LinearOperator;
    let cell = Cell::from_user_rows_bohr([[6., 0., 0.], [0.4, 5.5, 0.], [0.2, 0.3, 6.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.17, 0.23, 0.31])];
    let positions = [cell.fractional_to_cartesian_bohr(atoms[0].fractional_position)];
    let points = path();
    for (functional, spin_functional, name) in [
        (Functional::LdaPw, SpinFunctional::LdaPw, "SLA PW"),
        (Functional::Pbe, SpinFunctional::Pbe, "PBE"),
    ] {
        let pp = common::potential_with_nlcc(true, name);
        let potentials = [SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let options = options(functional, 6);
        let problem = ScfProblem::with_kpoints(
            &cell,
            &atoms,
            &potentials,
            options.clone(),
            &points,
            OccupationPolicy::Insulator {
                minimum_gap_hartree: 1e-7,
            },
        )
        .unwrap();
        let dims = problem.grid_dimensions();
        let channels = [1.25, 0.75].map(|charge| {
            Array3::from_shape_fn(dims, |(i, j, k)| {
                charge / cell.volume_bohr3()
                    * (1.
                        + 0.1
                            * (std::f64::consts::TAU
                                * (i as f64 / dims[0] as f64
                                    + j as f64 / dims[1] as f64
                                    + k as f64 / dims[2] as f64))
                                .cos())
            })
        });
        let before = channels.clone();
        // Legacy collinear's public serial field adapters, followed by a full
        // dense Hermitian solve (legacy nscf itself has no spin API).
        let mut hartree = SpinHartree::new(&cell, dims).unwrap();
        hartree.evaluate(&(&channels[0] + &channels[1])).unwrap();
        let mut xc = SpinXc::new(spin_functional, &cell, dims).unwrap();
        let mut response = SpinXcResponse::new(dims);
        xc.evaluate_with_core(&channels, problem.core_density(), &mut response)
            .unwrap();
        let mut expected = Vec::new();
        for xc_potential in &response.potentials {
            let v = problem.ionic_potential() + &hartree.potential + xc_potential;
            let mut spectra = Vec::new();
            for point in &points {
                let basis = dft_basis::PlaneWaveBasis::at_k(
                    &cell,
                    Hartree(options.cutoff_hartree),
                    point.fractional,
                )
                .unwrap();
                let nl = NonlocalOperator::new(
                    &basis,
                    cell.volume_bohr3(),
                    &[Species {
                        potential: &pp,
                        positions_bohr: &positions,
                    }],
                )
                .unwrap();
                let h = Hamiltonian::new(&cell, &basis, v.clone(), nl).unwrap();
                let n = basis.vectors().len();
                let mut identity = ComplexMatrix::zeros(n, n).unwrap();
                for i in 0..n {
                    identity[(i, i)] = Complex64::new(1., 0.);
                }
                let mut dense = ComplexMatrix::zeros(n, n).unwrap();
                h.apply_block(&identity, &mut dense).unwrap();
                let (values, vectors) = hermitian_eigh(&dense).unwrap();
                let mut hx = ComplexMatrix::zeros(n, n).unwrap();
                h.apply_block(&vectors, &mut hx).unwrap();
                for j in 0..6 {
                    let r = (0..n)
                        .map(|i| (hx[(i, j)] - values[j] * vectors[(i, j)]).norm_sqr())
                        .sum::<f64>()
                        .sqrt();
                    assert!(r < 1e-8, "serial dense residual {r:e} Ha");
                }
                spectra.push(values[..6].to_vec());
            }
            expected.push(spectra);
        }
        assert!((expected[0][0][0] - expected[1][0][0]).abs() > 1e-4);
        for ranks in [1, 3, 8] {
            let answers = dft_mpi::simulate(ranks, |comm| {
                nscf_with(
                    StartInput {
                        cell: &cell,
                        atoms: &atoms,
                        potentials: &potentials,
                        options: &options,
                    },
                    &channels,
                    &points,
                    6,
                    1e-8,
                    comm,
                )
                .unwrap()
            });
            for answer in answers {
                assert_eq!(answer.channels.len(), 2);
                assert!(
                    answer
                        .channels
                        .iter()
                        .all(|channel| channel.len() == points.len())
                );
                let mut error = 0_f64;
                let mut residual = 0_f64;
                for (actual, expected) in answer.channels.iter().zip(&expected) {
                    for (actual, expected) in actual.iter().zip(expected) {
                        for (a, b) in actual.eigenvalues.iter().zip(expected) {
                            error = error.max((a - b).abs());
                        }
                        residual = actual.residuals.iter().copied().fold(residual, f64::max);
                    }
                }
                assert!(error < 1e-8, "spin {name} ranks={ranks}: {error:e} Ha");
                assert!(residual <= 1e-8);
                println!(
                    "spin {name} ranks={ranks}: max eigenvalue error {error:e} Ha; max residual {residual:e} Ha"
                );
            }
        }
        assert_eq!(channels, before);
    }
}
