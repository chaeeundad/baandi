//! M1d: SCF with the distributed real Davidson on simulated ranks agrees with
//! the threads-only run within the accuracy class, and every rank agrees exactly.
#[allow(dead_code)]
mod common;
use dft_engine::{GammaReal, ScfOptions, ScfProblem, ScfResult, SpeciesPotential, parallel};
use dft_geometry::{Atom, Cell};
use dft_mpi::{Collective, SimRank, simulate};
use std::sync::Arc;

fn run(kind: usize, comm: Option<Arc<dyn Collective + Send + Sync>>) -> ScfResult {
    parallel::set_communicator(comm);
    let rows = if kind == 0 {
        [[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]
    } else {
        [[5., 0., 0.], [0.5, 4., 0.], [0.2, 0.4, 6.]]
    };
    let cell = Cell::from_user_rows_bohr(rows).unwrap();
    let atoms = vec![Atom::new("Xx", [0.13, 0.21, 0.32])];
    let pp = common::potential(2., 1., true, 401);
    let opts = ScfOptions {
        gamma_real: GammaReal::Auto,
        cutoff_hartree: 2.5,
        bands: Some(5),
        eigen_tolerance: 1e-10,
        ..Default::default()
    };
    let result = ScfProblem::new(
        &cell,
        &atoms,
        &[SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }],
        opts,
    )
    .unwrap()
    .run()
    .unwrap();
    parallel::set_communicator(None);
    result
}

#[test]
fn distributed_scf_matches_threads_only() {
    for kind in 0..2 {
        let serial = run(kind, None);
        for ranks in [2usize, 3] {
            let before = parallel::distributed_solves();
            let results = simulate(ranks, |comm: &SimRank| {
                run(kind, Some(Arc::new(comm.clone())))
            });
            assert!(
                parallel::distributed_solves() > before,
                "distributed path not taken"
            );
            for r in &results {
                let de = (r.total_hartree() - serial.total_hartree()).abs();
                assert!(de < 5e-9, "kind {kind} P={ranks}: dE {de:e}");
                assert_eq!(
                    r.total_hartree().to_bits(),
                    results[0].total_hartree().to_bits(),
                    "ranks disagree"
                );
            }
        }
    }
}
