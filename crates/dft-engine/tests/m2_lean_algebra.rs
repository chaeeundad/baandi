//! MPI-OVERLAP-DAVIDSON (`BAANDI_M2_LEAN_ALGEBRA=on`,
//! docs/design/mpi-overlap-davidson.md) in the production standard SCF
//! (fixed-point verification) on the M2 path: on the adopted defaults and on
//! the pending stack (mirror sticks, Gram start, kept start), cubic and skewed
//! cells, 2 and 3 simulated ranks: standard class against the threads-only run
//! (5e-9 Ha), converged and verified, ranks agree bit for bit and a fixed rank
//! count is bit-reproducible. Own test binary with a lock (process-wide
//! counters and per-thread overrides).
#[allow(dead_code)]
mod common;
use dft_engine::{GammaReal, ScfOptions, ScfProblem, ScfResult, SpeciesPotential, parallel};
use dft_geometry::{Atom, Cell};
use dft_mpi::{Collective, SimRank, simulate};
use dft_solvers::distributed_davidson::{
    DistributedTuning, lean_algebra_stats, set_start_gram_override, set_tuning_override,
};
use std::sync::{Arc, Mutex};

static LOCK: Mutex<()> = Mutex::new(());

fn cell(kind: usize) -> Cell {
    Cell::from_user_rows_bohr(if kind == 0 {
        [[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]]
    } else {
        [[6.5, 0., 0.], [0.6, 5.8, 0.], [0.3, 0.5, 7.]]
    })
    .unwrap()
}

fn run(kind: usize, pp: &dft_upf::Pseudopotential) -> ScfResult {
    let atoms = vec![
        Atom::new("Xx", [0.13, 0.21, 0.32]),
        Atom::new("Xx", [0.61, 0.57, 0.74]),
    ];
    ScfProblem::new(
        &cell(kind),
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
            ..ScfOptions::default()
        },
    )
    .unwrap()
    .run()
    .unwrap()
}

#[derive(Clone, Copy)]
struct Arm {
    lean: bool,
    stack: bool,
}

fn results(kind: usize, ranks: usize, a: Arm, pp: &dft_upf::Pseudopotential) -> Vec<ScfResult> {
    if ranks == 1 {
        return vec![run(kind, pp)];
    }
    simulate(ranks, |comm: &SimRank| {
        use dft_solvers::distributed_fft::set_mirror_sticks_override;
        parallel::set_m2_override(Some(true));
        set_tuning_override(Some(DistributedTuning {
            lean_algebra: a.lean,
            ..DistributedTuning::from_env()
        }));
        set_mirror_sticks_override(Some(a.stack));
        set_start_gram_override(Some(a.stack));
        dft_engine::gamma_real::set_m2_reuse_override(Some(a.stack));
        let comm: Arc<dyn Collective + Send + Sync> = Arc::new(comm.clone());
        parallel::set_communicator(Some(comm));
        let r = run(kind, pp);
        parallel::set_communicator(None);
        dft_engine::gamma_real::set_m2_reuse_override(None);
        set_start_gram_override(None);
        set_mirror_sticks_override(None);
        set_tuning_override(None);
        parallel::set_m2_override(None);
        r
    })
}

fn density_bits(r: &ScfResult) -> Vec<u64> {
    r.density.iter().map(|x| x.to_bits()).collect()
}

fn work(r: &ScfResult) -> (usize, usize, usize, usize) {
    (
        r.history.len(),
        r.history.iter().map(|h| h.davidson_iterations).sum(),
        r.history.iter().map(|h| h.operator_applications).sum(),
        r.verification.as_ref().map_or(0, |v| v.attempts.len()),
    )
}

fn verified(r: &ScfResult) -> bool {
    r.verification.as_ref().is_some_and(|v| v.passed())
}

#[test]
fn lean_algebra_keeps_the_standard_class_in_scf() {
    let _lock = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let pp = common::potential(4., 1., true, 401);
    for kind in 0..2 {
        let reference = results(
            kind,
            1,
            Arm {
                lean: false,
                stack: false,
            },
            &pp,
        )
        .remove(0);
        for ranks in [2usize, 3] {
            for stack in [false, true] {
                let plain = results(kind, ranks, Arm { lean: false, stack }, &pp);
                let before = lean_algebra_stats();
                let lean = results(kind, ranks, Arm { lean: true, stack }, &pp);
                let after = lean_algebra_stats();
                let again = results(kind, ranks, Arm { lean: true, stack }, &pp);
                assert!(after.solves > before.solves, "the lean algebra did not run");
                for r in 0..ranks {
                    let de = (lean[r].total_hartree() - reference.total_hartree()).abs();
                    let dp = (lean[r].total_hartree() - plain[r].total_hartree()).abs();
                    println!(
                        "kind {kind} P={ranks} stack {stack} rank {r}: lean dE vs threads-only {de:e} Ha, vs orthonormal M2 {dp:e} Ha, work {:?} (orthonormal {:?}, threads-only {:?}); lean solves {}, border sums {}, re-projections {}, fallbacks {}, rebuilds {}, wasted H {}, factored starts {}",
                        work(&lean[r]),
                        work(&plain[r]),
                        work(&reference),
                        (after.solves - before.solves) / ranks as u64,
                        (after.border_reductions - before.border_reductions) / ranks as u64,
                        (after.reprojections - before.reprojections) / ranks as u64,
                        (after.fallbacks - before.fallbacks) / ranks as u64,
                        (after.rebuilds - before.rebuilds) / ranks as u64,
                        (after.wasted_columns - before.wasted_columns) / ranks as u64,
                        (after.factored_starts - before.factored_starts) / ranks as u64,
                    );
                    assert!(verified(&lean[r]));
                    assert!(de < 5e-9, "standard class: dE {de:e} Ha");
                    assert_eq!(
                        lean[r].total_hartree().to_bits(),
                        lean[0].total_hartree().to_bits(),
                        "ranks disagree"
                    );
                    assert_eq!(density_bits(&lean[r]), density_bits(&again[r]));
                    assert_eq!(
                        lean[r].total_hartree().to_bits(),
                        again[r].total_hartree().to_bits()
                    );
                }
            }
        }
    }
}
