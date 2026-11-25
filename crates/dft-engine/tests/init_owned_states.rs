//! F-INIT-STATES: sparse construction, lazy fallbacks and the atomic H Q bridge.
//! Wide SCF sweep: cargo test --release -p dft-engine --test init_owned_states
//! -- --include-ignored --nocapture (uses simulated ranks, no MPI installation).
#[allow(dead_code)]
mod common;
use dft_basis::{FftGridSpec, HalfSphere, PlaneWaveBasis};
use dft_engine::{GammaReal, ScfOptions, ScfProblem, SpeciesPotential, parallel};
use dft_geometry::{Atom, Cell};
use dft_mpi::{Collective, SimRank, simulate};
use dft_numerics::{
    linalg::ComplexMatrix,
    pool::WorkerPool,
    switches::{RunConfig, RunCounter},
};
use dft_physics::{
    data::{Array3, Complex64 as C},
    hamiltonian::{Hamiltonian, NonlocalOperator, Species},
};
use dft_solvers::LinearOperator;
use dft_units::Hartree;
use std::sync::Arc;

fn cell(skew: bool) -> Cell {
    Cell::from_user_rows_bohr(if skew {
        [[6.2, 0., 0.], [1.1, 5.7, 0.], [0.4, 0.7, 6.6]]
    } else {
        [[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]]
    })
    .unwrap()
}

#[test]
fn owned_rows_atomic_image_and_lazy_full_table_on_one_to_seven_ranks() {
    let pp = common::potential(4., 1., true, 101);
    let mut max_error = 0_f64;
    for skew in [false, true] {
        // More than one 1024-G projection tile, with appreciable projectors
        // beyond the first tile (large cell rather than only a high cutoff).
        let cell = Cell::from_user_rows_bohr(if skew {
            [[9.3, 0., 0.], [1.65, 8.55, 0.], [0.6, 1.05, 9.9]]
        } else {
            [[9., 0., 0.], [0., 9., 0.], [0., 0., 9.]]
        })
        .unwrap();
        let basis = PlaneWaveBasis::new(&cell, Hartree(12.)).unwrap();
        let density = basis.density_basis(&cell).unwrap();
        let grid = FftGridSpec::new(&basis, &density, None).unwrap();
        let half = HalfSphere::new(&basis, &grid).unwrap();
        let positions = [[0.3, 0.7, 1.2], [2.1, 2.7, 3.4], [4., 1., 2.]];
        let species = [Species {
            potential: &pp,
            positions_bohr: &positions,
        }];
        let n = basis.momenta().len();
        assert!(n > 1024, "exercise the fixed tile boundary");
        let mut q = ComplexMatrix::zeros(n, 3).unwrap();
        for b in 0..3 {
            let mut real: Vec<_> = (0..n)
                .map(|i| ((i * 17 + b * 7 + 1) as f64).sin())
                .collect();
            let norm = real.iter().map(|x| x * x).sum::<f64>().sqrt();
            real.iter_mut().for_each(|x| *x /= norm);
            half.unpack(&real, q.column_mut(b)).unwrap();
        }
        let potential = Array3::from_shape_fn(grid.dims(), |(i, j, k)| {
            0.1 * ((i + j * 2 + k * 3) as f64).cos()
        });
        let eager = Hamiltonian::new_with_grid(
            &basis,
            grid.clone(),
            potential.clone(),
            NonlocalOperator::new(&basis, cell.volume_bohr3(), &species).unwrap(),
        )
        .unwrap();
        let mut expected = ComplexMatrix::zeros(n, 3).unwrap();
        eager.apply_block(&q, &mut expected).unwrap();
        let mut first_image: Option<Vec<C>> = None;
        let mut first_weights: Option<Vec<u64>> = None;
        for ranks in 1..=7 {
            let images = simulate(ranks, |comm: &SimRank| {
                let run = RunConfig::from_pairs(&[("init_owned_states", "on")]).unwrap();
                let _scope = run.enter();
                let indices = dft_physics::distributed_real_gamma::owned_projector_indices(
                    &basis,
                    &grid,
                    comm.rank(),
                    comm.size(),
                )
                .unwrap();
                let nl = NonlocalOperator::new_owned_with_rule(
                    &basis,
                    cell.volume_bohr3(),
                    &species,
                    dft_upf::radial::RadialEvenRule::Simpson3_8,
                    &indices,
                )
                .unwrap();
                assert!(!nl.full_table_built());
                assert_eq!(nl.table_bytes(), indices.len() * nl.projector_count() * 16);
                let h = Hamiltonian::new_with_grid(&basis, grid.clone(), potential.clone(), nl)
                    .unwrap();
                // Distributed construction must not call columns() even for count.
                let _dist = dft_physics::distributed_real_gamma::DistributedRealHamiltonian::new(
                    &h,
                    &basis,
                    4,
                    comm.rank(),
                    ranks,
                    WorkerPool::new(1).unwrap(),
                )
                .unwrap();
                assert!(!h.nonlocal().full_table_built());
                let weights = dft_engine::m2::owned_atomic_weights(comm, &h, &q).unwrap();
                let weight_bits: Vec<_> = weights
                    .as_slice()
                    .iter()
                    .flat_map(|z| [z.re.to_bits(), z.im.to_bits()])
                    .collect();
                let mut actual = ComplexMatrix::zeros(n, 3).unwrap();
                h.with_atomic_weights(weights, || h.apply_block(&q, &mut actual))
                    .unwrap()
                    .unwrap();
                let error = actual
                    .as_slice()
                    .iter()
                    .zip(expected.as_slice())
                    .map(|(a, b)| (*a - *b).norm())
                    .fold(0., f64::max);
                assert!(
                    error <= 1e-12,
                    "H Q error {error:e} Ha, ranks {ranks}, skew {skew}"
                );
                for b in 0..3 {
                    let energy: C = q
                        .column(b)
                        .iter()
                        .zip(actual.column(b).iter().zip(expected.column(b)))
                        .map(|(q, (a, b))| q.conj() * (a - b))
                        .sum();
                    assert!(energy.norm() <= 1e-12, "Rayleigh energy {energy} Ha");
                }
                assert_eq!(run.counters()[RunCounter::StatesFullBuilds as usize], 0);
                let before_fallback = run.counters();
                // Error exits must also clear the scoped weights.
                let injected = h.with_atomic_weights(
                    ComplexMatrix::zeros(h.nonlocal().projector_count(), 3).unwrap(),
                    || Err::<(), _>("injected"),
                );
                assert_eq!(injected.unwrap(), Err("injected"));
                // A complex apply is the lazy fallback. Its columns must be the
                // original values bit for bit, and repeated access builds once.
                let mut fallback = ComplexMatrix::zeros(n, 3).unwrap();
                h.apply_block(&q, &mut fallback).unwrap();
                assert_eq!(fallback, expected);
                h.apply_block(&q, &mut fallback).unwrap();
                assert_eq!(run.counters()[RunCounter::StatesFullBuilds as usize], 1);
                assert!(h.nonlocal().full_table_built());
                (
                    actual.as_slice().to_vec(),
                    error,
                    before_fallback,
                    weight_bits,
                )
            });
            for (image, error, _, weights) in &images {
                if let Some(first) = &first_weights {
                    assert_eq!(
                        first, weights,
                        "fixed chunks: rank-count-independent weights"
                    );
                } else {
                    first_weights = Some(weights.clone());
                }
                max_error = max_error.max(*error);
                if let Some(ref first) = first_image {
                    let d = first
                        .iter()
                        .zip(image)
                        .map(|(a, b)| (*a - *b).norm())
                        .fold(0., f64::max);
                    assert!(d <= 1e-12, "cross-rank H Q {d:e} Ha");
                } else {
                    first_image = Some(image.clone());
                }
            }
            let bytes: Vec<_> = images
                .iter()
                .map(|r| r.2[RunCounter::StatesTableBytes as usize])
                .collect();
            println!(
                "skew={skew} ranks={ranks} owned_bytes={bytes:?} replicated_bytes={}",
                images[0].2[RunCounter::StatesReplicatedBytes as usize]
            );
        }
    }
    println!("max H Q error {max_error:.16e} Ha");
}

type ScfOutput = (
    dft_engine::ScfResult,
    [u64; dft_numerics::switches::RUN_COUNTERS],
    Vec<u64>,
);
fn scf(ranks: usize, skew: bool, on: bool, fallback: bool) -> Vec<ScfOutput> {
    scf_threads(ranks, skew, on, fallback, 1)
}
fn scf_threads(
    ranks: usize,
    skew: bool,
    on: bool,
    fallback: bool,
    threads: usize,
) -> Vec<ScfOutput> {
    let pp = common::potential_with_orbitals(4., 1., true, 101, "SLA PZ");
    simulate(ranks, |comm: &SimRank| {
        let run =
            RunConfig::from_pairs(&[("init_owned_states", if on { "on" } else { "off" })]).unwrap();
        let _scope = run.enter();
        parallel::set_m2_override(Some(true));
        parallel::set_communicator(Some(Arc::new(comm.clone())));
        dft_engine::gamma_real::set_m2_dist_disabled_for_tests(fallback);
        // Exercise M1 fallback with a full potential, as in m2_lazy_serial.rs.
        dft_engine::grid_dist::set_owned_override(if fallback {
            Some((false, false))
        } else {
            None
        });
        let problem = ScfProblem::new(
            &cell(skew),
            &[
                Atom::new("Xx", [0.13, 0.21, 0.32]),
                Atom::new("Xx", [0.61, 0.57, 0.74]),
            ],
            &[SpeciesPotential {
                label: "Xx",
                potential: &pp,
            }],
            ScfOptions {
                gamma_real: GammaReal::Auto,
                cutoff_hartree: 3.,
                bands: Some(6),
                threads,
                ..ScfOptions::strict_reference()
            },
        )
        .unwrap();
        let ionic_bits = problem
            .ionic_potential()
            .iter()
            .map(|v| v.to_bits())
            .collect();
        let result = problem.run().unwrap();
        parallel::set_communicator(None);
        parallel::set_m2_override(None);
        dft_engine::gamma_real::set_m2_dist_disabled_for_tests(false);
        dft_engine::grid_dist::set_owned_override(None);
        (result, run.counters(), ionic_bits)
    })
}
fn compare_scf(ranks: usize, skew: bool, fallback: bool) {
    let a = scf(ranks, skew, false, fallback);
    let b = scf(ranks, skew, true, fallback);
    for ((a, _, ionic_a), (b, counts, ionic_b)) in a.iter().zip(&b) {
        assert_eq!(
            ionic_a, ionic_b,
            "root-only ionic transform is bit-identical"
        );
        let delta = (a.total_hartree() - b.total_hartree()).abs();
        assert!(
            delta <= 1e-12,
            "SCF delta {delta:e} Ha, ranks={ranks} skew={skew} fallback={fallback}"
        );
        assert_eq!(a.history.len(), b.history.len(), "SCF counts");
        assert_eq!(
            a.atomic_start.initial_h_applications,
            b.atomic_start.initial_h_applications
        );
        assert!(a.atomic_start.atomic_orbitals > 0);
        assert!(b.certified);
        if ranks > 1 {
            assert_eq!(
                counts[RunCounter::StatesFullBuilds as usize],
                u64::from(fallback)
            );
        }
        println!(
            "SCF ranks={ranks} skew={skew} fallback={fallback} dE={delta:.16e} Ha steps={} full_builds={}",
            b.history.len(),
            counts[RunCounter::StatesFullBuilds as usize]
        );
    }
}
#[test]
fn production_and_m1_fallback() {
    compare_scf(2, false, false);
    compare_scf(2, false, true);
}
#[test]
fn scf_one_to_seven_ranks_including_skew() {
    for ranks in 1..=7 {
        for skew in [false, true] {
            compare_scf(ranks, skew, false);
        }
    }
}

#[test]
fn recipe_preserves_both_radial_rules_empty_channels_and_multiple_species() {
    let cell = cell(true);
    let basis = PlaneWaveBasis::new(&cell, Hartree(3.)).unwrap();
    let a = common::potential(4., 1., true, 100);
    let b = common::potential(2., 0.7, false, 100);
    let positions = [[0.2, 0.3, 0.4], [1., 2., 3.]];
    let species = [
        Species {
            potential: &a,
            positions_bohr: &positions,
        },
        Species {
            potential: &b,
            positions_bohr: &positions,
        },
    ];
    for rule in [
        dft_upf::radial::RadialEvenRule::Simpson3_8,
        dft_upf::radial::RadialEvenRule::SimpsonQuadraticEnd,
    ] {
        let eager =
            NonlocalOperator::new_with_rule(&basis, cell.volume_bohr3(), &species, rule).unwrap();
        for selected in [
            vec![],
            vec![0],
            (0..basis.momenta().len()).step_by(3).collect(),
        ] {
            let owned = NonlocalOperator::new_owned_with_rule(
                &basis,
                cell.volume_bohr3(),
                &species,
                rule,
                &selected,
            )
            .unwrap();
            let mut column = 0;
            for atom in eager.atoms() {
                for col in atom.columns() {
                    for &i in &selected {
                        let z = owned.owned_value(column, i).unwrap();
                        assert_eq!(
                            (z.re.to_bits(), z.im.to_bits()),
                            (col[i].re.to_bits(), col[i].im.to_bits())
                        );
                    }
                    column += 1;
                }
            }
            assert!(!owned.full_table_built());
            for (a, b) in eager.atoms().iter().zip(owned.atoms()) {
                assert_eq!(a.coupling_hartree(), b.coupling_hartree());
                for (a, b) in a.columns().iter().zip(b.columns()) {
                    let bits = |v: &[C]| {
                        v.iter()
                            .flat_map(|z| [z.re.to_bits(), z.im.to_bits()])
                            .collect::<Vec<_>>()
                    };
                    assert_eq!(bits(a), bits(b));
                }
            }
        }
    }
}

#[test]
fn owned_construction_survives_pool_dispatch() {
    let one = scf_threads(2, true, true, false, 1);
    let two = scf_threads(2, true, true, false, 2);
    for ((a, ca, ia), (b, cb, ib)) in one.iter().zip(&two) {
        assert_eq!(ia, ib);
        assert_eq!(a.total_hartree().to_bits(), b.total_hartree().to_bits());
        assert_eq!(a.density, b.density);
        assert_eq!(a.history.len(), b.history.len());
        for c in [ca, cb] {
            assert_eq!(c[RunCounter::StatesFullBuilds as usize], 0);
            assert!(
                c[RunCounter::StatesTableBytes as usize]
                    < c[RunCounter::StatesReplicatedBytes as usize]
            );
        }
    }
}

#[test]
fn initial_atomic_ritz_energies_match() {
    let pp = common::potential_with_orbitals(4., 1., true, 101, "SLA PZ");
    let mut worst = 0_f64;
    for skew in [false, true] {
        let cell = cell(skew);
        let atoms = [
            Atom::new("Xx", [0.13, 0.21, 0.32]),
            Atom::new("Xx", [0.61, 0.57, 0.74]),
        ];
        let positions: Vec<_> = atoms
            .iter()
            .map(|a| cell.fractional_to_cartesian_bohr(a.fractional_position))
            .collect();
        for ranks in 1..=7 {
            let probe = |on| {
                simulate(ranks, |comm: &SimRank| {
                    let config = RunConfig::from_pairs(&[(
                        "init_owned_states",
                        if on { "on" } else { "off" },
                    )])
                    .unwrap();
                    let _scope = config.enter();
                    parallel::set_m2_override(Some(true));
                    parallel::set_communicator(Some(Arc::new(comm.clone())));
                    let mut problem = ScfProblem::new(
                        &cell,
                        &atoms,
                        &[SpeciesPotential {
                            label: "Xx",
                            potential: &pp,
                        }],
                        ScfOptions {
                            gamma_real: GammaReal::Auto,
                            cutoff_hartree: 3.,
                            bands: Some(6),
                            threads: 1,
                            ..ScfOptions::strict_reference()
                        },
                    )
                    .unwrap();
                    let rho = problem.uniform_density();
                    let guesses = problem.atomic_start_vectors(&rho).unwrap().unwrap();
                    assert_eq!(config.counters()[RunCounter::StatesFullBuilds as usize], 0);
                    // Independent Rayleigh evaluation of the returned atomic Ritz
                    // vectors: uniform density has zero Hartree potential (G=0
                    // omitted), no NLCC, and a constant PZ XC potential.
                    let xc = dft_physics::xc::lda_pz(rho[[0, 0, 0]]).unwrap().v_xc();
                    let potential = problem.ionic_potential().mapv(|v| v + xc);
                    let basis = PlaneWaveBasis::new(&cell, Hartree(3.)).unwrap();
                    let grid = FftGridSpec::new(
                        &basis,
                        &basis.density_basis(&cell).unwrap(),
                        Some(problem.grid_dimensions()),
                    )
                    .unwrap();
                    let reference = Hamiltonian::new_with_grid(
                        &basis,
                        grid,
                        potential,
                        NonlocalOperator::new(
                            &basis,
                            cell.volume_bohr3(),
                            &[Species {
                                potential: &pp,
                                positions_bohr: &positions,
                            }],
                        )
                        .unwrap(),
                    )
                    .unwrap();
                    let q = &guesses[0];
                    let mut out = vec![C::default(); q.rows()];
                    let energies: Vec<_> = (0..q.cols())
                        .map(|b| {
                            reference.apply(q.column(b), &mut out).unwrap();
                            q.column(b)
                                .iter()
                                .zip(&out)
                                .map(|(x, y)| (x.conj() * y).re)
                                .sum::<f64>()
                        })
                        .collect();
                    parallel::set_communicator(None);
                    parallel::set_m2_override(None);
                    energies
                })
            };
            let a = probe(false);
            let b = probe(true);
            for (a, b) in a.iter().flatten().zip(b.iter().flatten()) {
                let error = (a - b).abs();
                assert!(
                    error <= 1e-12,
                    "initial atomic Ritz energy {error:e} Ha, ranks={ranks}, skew={skew}"
                );
                worst = worst.max(error);
            }
        }
    }
    println!("max initial atomic Ritz energy difference {worst:.16e} Ha");
}
