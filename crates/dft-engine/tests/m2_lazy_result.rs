//! LAZY-RESULT (`BAANDI_M2_LAZY_RESULT`, `gamma_real::set_delivered_orbitals`):
//! the production standard SCF on simulated ranks delivers its orbitals lazily
//! when the caller declared that it reads none. Own test binary: it reads
//! process-wide M2 counters.
#[allow(dead_code)]
mod common;
use dft_engine::{GammaReal, ScfOptions, ScfProblem, ScfResult, SpeciesPotential, parallel};
use dft_geometry::{Atom, Cell};
use dft_mpi::{Collective, SimRank, simulate};
use std::sync::Arc;

fn m2_count(name: &str) -> usize {
    parallel::m2_stats()
        .0
        .into_iter()
        .find(|(n, _)| n == name)
        .map_or(0, |(_, v)| v)
}

fn cell(kind: usize) -> Cell {
    Cell::from_user_rows_bohr(if kind == 0 {
        [[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]]
    } else {
        [[6.5, 0., 0.], [0.6, 5.8, 0.], [0.3, 0.5, 7.]]
    })
    .unwrap()
}

fn problem(kind: usize, pp: &dft_upf::Pseudopotential) -> ScfProblem {
    let cell = cell(kind);
    let atoms = vec![
        Atom::new("Xx", [0.13, 0.21, 0.32]),
        Atom::new("Xx", [0.61, 0.57, 0.74]),
    ];
    ScfProblem::new(
        &cell,
        &atoms,
        &[SpeciesPotential {
            label: "Xx",
            potential: pp,
        }],
        ScfOptions {
            gamma_real: GammaReal::Auto,
            cutoff_hartree: 3.,
            bands: Some(6),
            threads: 1,
            ..ScfOptions::strict_reference()
        },
    )
    .unwrap()
}

fn orbital_bits(r: &ScfResult) -> (Vec<u64>, Vec<u64>) {
    let flat = |m: &dft_numerics::linalg::ComplexMatrix| {
        m.as_slice()
            .iter()
            .flat_map(|z| [z.re.to_bits(), z.im.to_bits()])
            .collect::<Vec<_>>()
    };
    (
        flat(&r.orbitals),
        flat(&r.k_points[0].solution.eigenvectors),
    )
}

/// With `set_delivered_orbitals(false)` the M2 result keeps its lazy blocks
/// (counter `lazy_results`, `orbitals_delivered() == false`, forces refuse
/// it), and its energy, density, eigenvalues and history are bit-identical to
/// the default run; `materialize_orbitals` (collective) then yields the
/// default run's orbitals bit for bit. Cubic and skewed cells, 2 and 3 ranks.
#[test]
fn lazy_result_matches_the_gathered_result() {
    let pp = common::potential(4., 1., true, 401);
    for kind in 0..2 {
        for ranks in [2usize, 3] {
            let arm = |needed: bool| {
                simulate(ranks, |comm: &SimRank| {
                    let comm: Arc<dyn Collective + Send + Sync> = Arc::new(comm.clone());
                    parallel::set_m2_override(Some(true));
                    dft_engine::gamma_real::set_m2_lazy_override(Some(true));
                    parallel::set_communicator(Some(comm));
                    dft_engine::gamma_real::set_delivered_orbitals(needed);
                    let mut r = problem(kind, &pp).run().unwrap();
                    let delivered = r.orbitals_delivered();
                    let refused = !delivered
                        && dft_engine::forces::forces(&dft_engine::forces::ForceState {
                            cell: &cell(kind),
                            atoms: &[],
                            potentials: &[],
                            options: &ScfOptions::strict_reference(),
                            electronic: dft_engine::forces::ElectronicState::Fixed(&r),
                            density_floor: None,
                            symmetry: &[],
                        })
                        .is_err();
                    let before = orbital_bits(&r);
                    r.materialize_orbitals().unwrap();
                    dft_engine::gamma_real::set_delivered_orbitals(true);
                    parallel::set_communicator(None);
                    dft_engine::gamma_real::set_m2_lazy_override(None);
                    parallel::set_m2_override(None);
                    (r, delivered, refused, before)
                })
            };
            let full = arm(true);
            let results0 = m2_count("lazy_results");
            let lazy = arm(false);
            assert_eq!(m2_count("lazy_results") - results0, ranks, "one per rank");
            for ((a, da, _, _), (b, db, refused, before)) in full.iter().zip(&lazy) {
                assert!(*da, "default run delivers orbitals");
                assert!(!*db, "lazy run keeps them distributed");
                assert!(*refused, "forces refuse a lazy result");
                assert_eq!(before.0.len(), 0, "placeholder until materialized");
                assert_eq!(
                    a.energy.total_hartree().to_bits(),
                    b.energy.total_hartree().to_bits()
                );
                assert!(
                    a.density
                        .iter()
                        .zip(&b.density)
                        .all(|(x, y)| x.to_bits() == y.to_bits())
                );
                assert_eq!(
                    a.eigenvalues
                        .iter()
                        .map(|x| x.to_bits())
                        .collect::<Vec<_>>(),
                    b.eigenvalues
                        .iter()
                        .map(|x| x.to_bits())
                        .collect::<Vec<_>>()
                );
                assert_eq!(a.history.len(), b.history.len());
                for (x, y) in a.history.iter().zip(&b.history) {
                    assert_eq!(x.density_residual.to_bits(), y.density_residual.to_bits());
                    assert_eq!(x.operator_applications, y.operator_applications);
                }
                assert_eq!(orbital_bits(a), orbital_bits(b));
                assert!(b.orbitals_delivered());
            }
        }
    }
}
