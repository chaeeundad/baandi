//! M1d: SCF with the distributed real Davidson on simulated ranks agrees with
//! the threads-only run within the accuracy class, and every rank agrees exactly.
#[allow(dead_code)]
mod common;
use dft_engine::{GammaReal, ScfOptions, ScfProblem, ScfResult, SpeciesPotential, parallel};
use dft_geometry::{Atom, Cell};
use dft_mpi::{Collective, SimRank, simulate};
use std::sync::Arc;

fn run(kind: usize, comm: Option<Arc<dyn Collective + Send + Sync>>) -> ScfResult {
    run_with(kind, comm, 1, 5, 2.5)
}
fn run_with(
    kind: usize,
    comm: Option<Arc<dyn Collective + Send + Sync>>,
    threads: usize,
    bands: usize,
    cutoff: f64,
) -> ScfResult {
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
        cutoff_hartree: cutoff,
        bands: Some(bands),
        threads,
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

#[test]
fn distributed_with_threads_takes_the_distributed_path_and_keeps_iterations() {
    // threads >= 2 used to move the solve to a pool thread where the communicator
    // is not visible, silently running the serial solver on every rank.
    let serial = run_with(1, None, 2, 24, 4.0);
    let before = parallel::distributed_solves();
    let results = simulate(3, |comm: &SimRank| {
        run_with(1, Some(Arc::new(comm.clone())), 2, 24, 4.0)
    });
    assert!(
        parallel::distributed_solves() > before,
        "distributed path not taken with threads=2"
    );
    for r in &results {
        assert!((r.total_hartree() - serial.total_hartree()).abs() < 5e-9);
        println!(
            "davidson iterations: serial {} distributed {}",
            total_davidson(&serial),
            total_davidson(r)
        );
        assert!(
            total_davidson(r) <= total_davidson(&serial) + total_davidson(&serial) / 10 + 2,
            "distributed Davidson needs many more iterations: {} vs {}",
            total_davidson(r),
            total_davidson(&serial)
        );
    }
}

fn total_davidson(r: &ScfResult) -> usize {
    r.history.iter().map(|h| h.davidson_iterations).sum()
}
