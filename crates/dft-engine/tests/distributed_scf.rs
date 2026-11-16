//! M1d: SCF with the distributed real Davidson on simulated ranks agrees with
//! the threads-only run within the accuracy class, and every rank agrees exactly.
#[allow(dead_code)]
mod common;
use dft_engine::{GammaReal, ScfOptions, ScfProblem, ScfResult, SpeciesPotential, parallel};
use dft_geometry::{Atom, Cell};
use dft_mpi::{Collective, SimRank, simulate};
use std::sync::Arc;

/// The distributed solver keeps process-wide reduction counters and a process-wide
/// tuning override; concurrent tests in this binary would mix them, so every test
/// takes this lock (the tests are small).
fn serial() -> std::sync::MutexGuard<'static, ()> {
    static LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());
    LOCK.lock().unwrap_or_else(|e| e.into_inner())
}

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
        ..ScfOptions::strict_reference()
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
    let _serial = serial();
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
    let _serial = serial();
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

/// Review F02 (2026-11-16): with an active communicator and more than one k
/// point there is no k-point pool and no distributed solve: every rank solves
/// every k point on the complex path (replicated execution, reported by the CLI
/// as `parallel.model = "mpi_replicated"`). The replicated run is safe: the
/// real-Gamma (distributed) path is not taken, every rank agrees bit for bit and
/// the energy equals the threads-only run within the strict class (1e-10 Ha).
#[test]
fn multiple_k_points_run_replicated_on_every_rank() {
    let _serial = serial();
    use dft_engine::OccupationPolicy;
    use dft_geometry::kpoints::KPoint;
    let run_k = |comm: Option<Arc<dyn Collective + Send + Sync>>| {
        parallel::set_communicator(comm);
        let cell = Cell::from_user_rows_bohr([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]).unwrap();
        let pp = common::potential(2., 1., true, 401);
        let r = ScfProblem::with_kpoints(
            &cell,
            &[Atom::new("Xx", [0.13, 0.21, 0.32])],
            &[SpeciesPotential {
                label: "Xx",
                potential: &pp,
            }],
            ScfOptions {
                gamma_real: GammaReal::Auto,
                cutoff_hartree: 2.5,
                bands: Some(5),
                threads: 1,
                eigen_tolerance: 1e-10,
                ..ScfOptions::strict_reference()
            },
            &[
                KPoint {
                    fractional: [0., 0., 0.],
                    weight: 0.5,
                },
                KPoint {
                    fractional: [0.25, 0., 0.],
                    weight: 0.5,
                },
            ],
            OccupationPolicy::Insulator {
                minimum_gap_hartree: 1e-6,
            },
        )
        .unwrap()
        .run()
        .unwrap();
        parallel::set_communicator(None);
        r
    };
    let serial = run_k(None);
    let results = simulate(2, |comm: &SimRank| run_k(Some(Arc::new(comm.clone()))));
    for r in &results {
        assert!(
            !r.gamma_real_used,
            "multi-k must not take the real-Gamma path"
        );
        assert_eq!(
            r.total_hartree().to_bits(),
            results[0].total_hartree().to_bits(),
            "ranks disagree"
        );
        let de = (r.total_hartree() - serial.total_hartree()).abs();
        assert!(de < 1e-10, "replicated dE {de:e} Ha");
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
    let _serial = serial();
    // Cubic and skewed cells; the M2c path is selected per thread (no env mutation).
    assert!(
        parallel::m2_enabled(),
        "M2 is on by default since the adoption rule"
    );
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
        ..ScfOptions::strict_reference()
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
    let _serial = serial();
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
    let _serial = serial();
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

/// MPI entry cost: the distributed solver with every tuning switch on
/// (`BAANDI_M2_ORTHO=serial` with the serial panel and DGKS rules, `M2_DGEMM`,
/// `M2_LEAN`; forced per simulated rank) stays within the accuracy classes of the
/// threads-only run on two cell shapes, strict and standard, 2 and 3 ranks; every
/// rank agrees bit for bit and a fixed rank count is bit-reproducible.
#[test]
fn m2_tuned_solver_matches_threads_only() {
    let _serial = serial();
    use dft_solvers::distributed_davidson::{DistributedTuning, set_tuning_override};
    let tuned = DistributedTuning {
        serial_ortho: true,
        wide_gemm: true,
        lean: true,
        ..DistributedTuning::from_env()
    };
    for (standard, limit) in [(false, 1e-10), (true, 5e-9)] {
        for kind in 0..2 {
            let serial = run_m2_profile(kind, None, standard);
            for ranks in [2usize, 3] {
                let go = || {
                    simulate(ranks, |comm: &SimRank| {
                        set_tuning_override(Some(tuned));
                        let r = run_m2_profile(kind, Some(Arc::new(comm.clone())), standard);
                        set_tuning_override(None);
                        r
                    })
                };
                let results = go();
                let again = go();
                for (r, res) in results.iter().enumerate() {
                    let de = (res.total_hartree() - serial.total_hartree()).abs();
                    println!(
                        "tuned standard={standard} kind {kind} P={ranks} rank {r}: dE {de:e} Ha, Davidson {} vs {}",
                        total_davidson(res),
                        total_davidson(&serial)
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
                }
            }
        }
    }
}

/// `BAANDI_M2_PANEL_GS=pythag|fused` (docs/performance/profile-2026-11-16-mpi.md):
/// on the default distributed stack the SCF energy stays within the accuracy
/// classes of the threads-only run (strict 1e-10 Ha, standard 5e-9 Ha), ranks
/// agree bit for bit, a fixed rank count is bit-reproducible, and the in-panel
/// reductions fall (fused: fewer tree sums than classic). Cubic and skewed
/// cells, 2 and 3 ranks, strict and standard.
#[test]
fn m2_panel_gs_variants_keep_classes() {
    let _serial = serial();
    use dft_solvers::distributed_davidson::{
        DistributedTuning, PanelGs, reduction_stats, set_tuning_override,
    };
    let base = DistributedTuning::from_env();
    let tree_sums = || -> u64 { reduction_stats().1.iter().map(|s| s.calls).sum() };
    for (standard, limit) in [(false, 1e-10), (true, 5e-9)] {
        for kind in 0..2 {
            let serial = run_m2_profile(kind, None, standard);
            for ranks in [2usize, 3] {
                let mut calls = Vec::new();
                for mode in [PanelGs::Classic, PanelGs::Pythag, PanelGs::Fused] {
                    let tuning = DistributedTuning {
                        panel_gs: mode,
                        ..base
                    };
                    let go = || {
                        simulate(ranks, |comm: &SimRank| {
                            set_tuning_override(Some(tuning));
                            let r = run_m2_profile(kind, Some(Arc::new(comm.clone())), standard);
                            set_tuning_override(None);
                            r
                        })
                    };
                    let before = tree_sums();
                    let results = go();
                    calls.push(tree_sums() - before);
                    let again = go();
                    for (r, res) in results.iter().enumerate() {
                        let de = (res.total_hartree() - serial.total_hartree()).abs();
                        println!(
                            "panel_gs {mode:?} standard={standard} kind {kind} P={ranks} rank {r}: dE {de:e} Ha, Davidson {} vs {}",
                            total_davidson(res),
                            total_davidson(&serial)
                        );
                        assert!(de < limit, "{mode:?} kind {kind} P={ranks}: dE {de:e}");
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
                    }
                }
                println!("tree sums classic/pythag/fused: {calls:?}");
                assert!(
                    calls[2] < calls[0],
                    "fused must reduce the tree sums: {calls:?}"
                );
            }
        }
    }
}

/// MPI-SCALING-GAP (docs/performance/mpi-scaling-gap.md): `M2_STICKS=mirror`
/// (mirror-paired stick partition, no mirror exchanges) and `M2_START=gram`
/// (orthonormal starts taken after one Gram reduction), alone and together, stay
/// within the accuracy classes of the threads-only run (strict 1e-10 Ha,
/// standard 5e-9 Ha); ranks agree bit for bit and a fixed rank count is
/// bit-reproducible. Cubic and skewed cells, 2 and 3 ranks, strict and standard.
#[test]
fn m2_scaling_gap_switches_keep_classes() {
    let _serial = serial();
    use dft_solvers::distributed_davidson::set_start_gram_override;
    use dft_solvers::distributed_fft::set_mirror_sticks_override;
    for (standard, limit) in [(false, 1e-10), (true, 5e-9)] {
        for kind in 0..2 {
            let serial = run_m2_profile(kind, None, standard);
            for ranks in [2usize, 3] {
                for (mirror, gram) in [(true, false), (false, true), (true, true)] {
                    let go = || {
                        simulate(ranks, |comm: &SimRank| {
                            set_mirror_sticks_override(Some(mirror));
                            set_start_gram_override(Some(gram));
                            let r = run_m2_profile(kind, Some(Arc::new(comm.clone())), standard);
                            set_mirror_sticks_override(None);
                            set_start_gram_override(None);
                            r
                        })
                    };
                    let results = go();
                    let again = go();
                    for (r, res) in results.iter().enumerate() {
                        let de = (res.total_hartree() - serial.total_hartree()).abs();
                        println!(
                            "sticks mirror={mirror} start gram={gram} standard={standard} kind {kind} P={ranks} rank {r}: dE {de:e} Ha"
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
                    }
                }
            }
        }
    }
}

/// MPI reduced eigenproblem (docs/performance/mpi-reduced-eigen-analysis.md):
/// on top of the entry-cost stack, `M2_EIGEN=root` (rank 0 solves, broadcast)
/// gives bit-identical energies and densities to the replicated solve, and
/// `M2_B5=on` (bordered update) stays within the accuracy classes of the
/// threads-only run (strict 1e-10 Ha, standard 5e-9 Ha); ranks agree bit for
/// bit and a fixed rank count is bit-reproducible. Cubic and skewed cells, 2 and
/// 3 ranks, strict and standard.
#[test]
fn m2_reduced_eigen_switches_keep_classes() {
    let _serial = serial();
    use dft_solvers::distributed_davidson::{DistributedTuning, set_tuning_override};
    let stack = DistributedTuning {
        serial_ortho: true,
        wide_gemm: true,
        lean: true,
        ..DistributedTuning::from_env()
    };
    let arm = |tuning: DistributedTuning, kind: usize, ranks: usize, standard: bool| {
        simulate(ranks, |comm: &SimRank| {
            set_tuning_override(Some(tuning));
            let r = run_m2_profile(kind, Some(Arc::new(comm.clone())), standard);
            set_tuning_override(None);
            r
        })
    };
    let density_bits = |r: &ScfResult| r.density.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    for (standard, limit) in [(false, 1e-10), (true, 5e-9)] {
        for kind in 0..2 {
            let serial = run_m2_profile(kind, None, standard);
            for ranks in [2usize, 3] {
                let replicated = arm(stack, kind, ranks, standard);
                let root = DistributedTuning {
                    eigen_root: true,
                    ..stack
                };
                let rooted = arm(root, kind, ranks, standard);
                for (a, b) in rooted.iter().zip(&replicated) {
                    assert_eq!(
                        a.total_hartree().to_bits(),
                        b.total_hartree().to_bits(),
                        "root broadcast changed the energy"
                    );
                    assert_eq!(density_bits(a), density_bits(b));
                }
                let b5 = DistributedTuning {
                    bordered: true,
                    eigen_root: true,
                    ..stack
                };
                let bordered = arm(b5, kind, ranks, standard);
                let again = arm(b5, kind, ranks, standard);
                for (r, res) in bordered.iter().enumerate() {
                    let de = (res.total_hartree() - serial.total_hartree()).abs();
                    println!(
                        "B5+root standard={standard} kind {kind} P={ranks} rank {r}: dE {de:e} Ha, Davidson {} vs {}",
                        total_davidson(res),
                        total_davidson(&serial)
                    );
                    assert!(de < limit, "kind {kind} P={ranks}: dE {de:e}");
                    assert_eq!(
                        res.total_hartree().to_bits(),
                        bordered[0].total_hartree().to_bits(),
                        "ranks disagree"
                    );
                    assert_eq!(
                        res.total_hartree().to_bits(),
                        again[r].total_hartree().to_bits(),
                        "not bit-reproducible at fixed rank count"
                    );
                }
            }
        }
    }
}

/// F06 (`BAANDI_M2_WORKSPACE=on`): the distributed solver on its persistent
/// per-rank workspace (reused across the solves of the whole SCF) gives energies
/// and densities bit-identical to the allocating path; cubic and skewed cells,
/// strict and standard, 2 and 3 ranks, on the adopted stack from the switches.
#[test]
fn m2_workspace_is_bit_identical_in_scf() {
    let _serial = serial();
    use dft_solvers::distributed_davidson::{
        DistributedTuning, release_workspace, set_tuning_override, workspace_stats,
    };
    let stack = DistributedTuning::from_env();
    let arm = |workspace: bool, kind: usize, ranks: usize, standard: bool| {
        simulate(ranks, |comm: &SimRank| {
            set_tuning_override(Some(DistributedTuning { workspace, ..stack }));
            let r = run_m2_profile(kind, Some(Arc::new(comm.clone())), standard);
            set_tuning_override(None);
            release_workspace();
            r
        })
    };
    let density_bits = |r: &ScfResult| r.density.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    let before = workspace_stats();
    for standard in [false, true] {
        for kind in 0..2 {
            for ranks in [2usize, 3] {
                let off = arm(false, kind, ranks, standard);
                let on = arm(true, kind, ranks, standard);
                for (a, b) in on.iter().zip(&off) {
                    assert_eq!(
                        a.total_hartree().to_bits(),
                        b.total_hartree().to_bits(),
                        "workspace changed the energy (standard={standard} kind {kind} P={ranks})"
                    );
                    assert_eq!(density_bits(a), density_bits(b));
                    assert_eq!(total_davidson(a), total_davidson(b));
                }
            }
        }
    }
    let after = workspace_stats();
    println!("workspace stats before {before:?} after {after:?}");
    assert!(after.warm_solves > before.warm_solves);
}
