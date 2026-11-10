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

fn run_m2(kind: usize, comm: Option<Arc<dyn Collective + Send + Sync>>) -> ScfResult {
    parallel::set_m2_override(Some(true));
    let r = run(kind, comm);
    parallel::set_m2_override(None);
    r
}

#[test]
fn m2c_stick_row_scf_matches_threads_only() {
    // Cubic and skewed cells; the M2c path is selected per thread (no env mutation).
    assert!(!parallel::m2_enabled(), "M2 must be off by default");
    for kind in 0..2 {
        let serial = run(kind, None);
        for ranks in [2usize, 3] {
            let before = parallel::distributed_rows_solves();
            let results = simulate(ranks, |comm: &SimRank| {
                run_m2(kind, Some(Arc::new(comm.clone())))
            });
            assert!(
                parallel::distributed_rows_solves() > before,
                "M2c path not taken"
            );
            for r in &results {
                let de = (r.total_hartree() - serial.total_hartree()).abs();
                println!("M2c kind {kind} P={ranks}: dE {de:e} Ha");
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

// ---- M2d: distributed density, XC, nonlocal energy on the M2 path ----

fn m2_count(name: &str) -> usize {
    parallel::m2_stats()
        .0
        .into_iter()
        .find(|(n, _)| n == name)
        .map_or(0, |(_, v)| v)
}

fn run_m2_profile(
    kind: usize,
    comm: Option<Arc<dyn Collective + Send + Sync>>,
    standard: bool,
) -> ScfResult {
    parallel::set_m2_override(Some(true));
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
        threads: 1,
        eigen_tolerance: 1e-10,
        standard: standard.then(dft_engine::StandardProfile::default),
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
    parallel::set_m2_override(None);
    result
}

#[test]
fn m2d_distributed_scf_stages_match_threads_only() {
    // Two cell shapes, strict (1e-10 Ha) and standard (5e-9 Ha) classes, 2..4 ranks.
    for (standard, limit) in [(false, 1e-10), (true, 5e-9)] {
        for kind in 0..2 {
            let serial = run_m2_profile(kind, None, standard);
            for ranks in [2usize, 3, 4] {
                let names = [
                    "rows_solves",
                    "batched_applies",
                    "band_density_builds",
                    "density_builds",
                    "xc_evaluations",
                    "nonlocal_energy_audits",
                ];
                let before: Vec<usize> = names.iter().map(|n| m2_count(n)).collect();
                let go = || {
                    simulate(ranks, |comm: &SimRank| {
                        run_m2_profile(kind, Some(Arc::new(comm.clone())), standard)
                    })
                };
                let results = go();
                for (name, b) in names.iter().zip(&before) {
                    assert!(m2_count(name) > *b, "M2d counter {name} did not move");
                }
                let again = go();
                for (r, res) in results.iter().enumerate() {
                    let de = (res.total_hartree() - serial.total_hartree()).abs();
                    println!(
                        "M2d standard={standard} kind {kind} P={ranks} rank {r}: dE {de:e} Ha, iterations {} vs {}",
                        res.history.len(),
                        serial.history.len()
                    );
                    assert!(de < limit, "kind {kind} P={ranks}: dE {de:e}");
                    assert_eq!(
                        res.total_hartree().to_bits(),
                        results[0].total_hartree().to_bits(),
                        "ranks disagree"
                    );
                    assert_eq!(
                        res.total_hartree().to_bits(),
                        again[r].total_hartree().to_bits(),
                        "not bit-reproducible at fixed rank count"
                    );
                    assert!(
                        res.density
                            .iter()
                            .zip(&again[r].density)
                            .all(|(a, b)| a.to_bits() == b.to_bits()),
                        "density not bit-reproducible at fixed rank count"
                    );
                }
            }
        }
    }
}

#[test]
fn m2_off_and_one_rank_keep_the_threads_only_bits() {
    // M2 requested on one rank (no communicator) must be the threads-only path.
    let serial = run(0, None);
    let one = run_m2_profile(0, None, false);
    assert_eq!(
        serial.total_hartree().to_bits(),
        one.total_hartree().to_bits()
    );
    assert!(
        serial
            .density
            .iter()
            .zip(&one.density)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
}
