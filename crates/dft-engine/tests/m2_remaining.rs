//! MPI-REMAINING (docs/performance/mpi-remaining.md): the replicated grid-stage
//! switches `BAANDI_MIX_INCREMENTAL` and `BAANDI_STOP_METRIC`, and the combined
//! stack of the pending MPI candidates (mirror sticks, Gram start, fused panel
//! Gram-Schmidt, kept start, lazy orbitals) on the production standard SCF
//! (fixed-point verification). Own test binary with a lock: it reads
//! process-wide counters.
#[allow(dead_code)]
mod common;
use dft_engine::{GammaReal, ScfOptions, ScfProblem, ScfResult, SpeciesPotential, parallel};
use dft_geometry::{Atom, Cell};
use dft_mpi::{Collective, SimRank, simulate};
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

/// Two atoms of four valence electrons, six bands, production standard path.
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

#[derive(Clone, Copy, Default)]
struct Arm {
    incremental: bool,
    kept_metric: bool,
    mirror: bool,
    gram: bool,
    reuse: bool,
    lazy: bool,
}

fn with_arm<R>(
    comm: Option<Arc<dyn Collective + Send + Sync>>,
    a: Arm,
    f: impl FnOnce() -> R,
) -> R {
    use dft_solvers::distributed_davidson::set_start_gram_override;
    use dft_solvers::distributed_fft::set_mirror_sticks_override;
    let m2 = comm.is_some();
    if m2 {
        parallel::set_m2_override(Some(true));
    }
    dft_solvers::mixing::set_incremental_override(Some(a.incremental));
    dft_engine::set_stop_metric_override(Some(a.kept_metric));
    set_mirror_sticks_override(Some(a.mirror));
    set_start_gram_override(Some(a.gram));
    dft_engine::gamma_real::set_m2_reuse_override(Some(a.reuse));
    dft_engine::gamma_real::set_m2_lazy_override(Some(a.lazy));
    parallel::set_communicator(comm);
    let r = f();
    parallel::set_communicator(None);
    dft_engine::gamma_real::set_m2_lazy_override(None);
    dft_engine::gamma_real::set_m2_reuse_override(None);
    set_start_gram_override(None);
    set_mirror_sticks_override(None);
    dft_engine::set_stop_metric_override(None);
    dft_solvers::mixing::set_incremental_override(None);
    parallel::set_m2_override(None);
    r
}

/// Results of every rank (one entry on the threads-only path).
fn results(kind: usize, ranks: usize, a: Arm, pp: &dft_upf::Pseudopotential) -> Vec<ScfResult> {
    if ranks == 1 {
        return vec![with_arm(None, a, || run(kind, pp))];
    }
    simulate(ranks, |comm: &SimRank| {
        with_arm(Some(Arc::new(comm.clone())), a, || run(kind, pp))
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

/// `BAANDI_MIX_INCREMENTAL=on` is bit-identical to off (energy, density, SCF /
/// Davidson / H / verification counts) on the threads-only path and on 2
/// simulated ranks; `BAANDI_STOP_METRIC=kept` takes the estimate from kept
/// coefficients (counter), stays in the standard class against fft (here the
/// trajectory is unchanged: |dE| <= 1e-12 Ha), is verified, ranks agree and
/// a fixed rank count is bit-reproducible. Cubic and skewed cells.
#[test]
fn grid_stage_switches_keep_the_state() {
    let _lock = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let pp = common::potential(4., 1., true, 401);
    for kind in 0..2 {
        for ranks in [1usize, 2] {
            let base = results(kind, ranks, Arm::default(), &pp);
            let inc = results(
                kind,
                ranks,
                Arm {
                    incremental: true,
                    ..Arm::default()
                },
                &pp,
            );
            let before = dft_engine::stop_metric_counters();
            let kept = results(
                kind,
                ranks,
                Arm {
                    kept_metric: true,
                    ..Arm::default()
                },
                &pp,
            );
            let after = dft_engine::stop_metric_counters();
            let again = results(
                kind,
                ranks,
                Arm {
                    kept_metric: true,
                    ..Arm::default()
                },
                &pp,
            );
            let hits = (after[0] - before[0]) as usize / ranks;
            let fallbacks = (after[1] - before[1]) as usize / ranks;
            for r in 0..ranks {
                assert!(verified(&base[r]) && verified(&inc[r]) && verified(&kept[r]));
                assert_eq!(
                    inc[r].total_hartree().to_bits(),
                    base[r].total_hartree().to_bits(),
                    "incremental mixer changed the energy: kind {kind} P={ranks}"
                );
                assert_eq!(density_bits(&inc[r]), density_bits(&base[r]));
                assert_eq!(work(&inc[r]), work(&base[r]));
                let de = (kept[r].total_hartree() - base[r].total_hartree()).abs();
                println!(
                    "kind {kind} P={ranks} rank {r}: kept-metric dE {de:e} Ha, work {:?} -> {:?}, estimates kept {hits} fft {fallbacks}",
                    work(&base[r]),
                    work(&kept[r])
                );
                assert!(de <= 1e-12, "kept metric dE {de:e} Ha");
                assert_eq!(work(&kept[r]), work(&base[r]));
                assert_eq!(
                    kept[r].total_hartree().to_bits(),
                    kept[0].total_hartree().to_bits(),
                    "ranks disagree"
                );
                assert_eq!(density_bits(&kept[r]), density_bits(&again[r]));
            }
            assert!(hits > 0, "no stop-rule estimate used kept coefficients");
        }
    }
}

/// MPI-REMAINING item (d): the pending MPI candidates together (mirror sticks,
/// Gram start, fused panel Gram-Schmidt from the default, kept start, lazy
/// orbitals) plus the two grid switches stay in the standard class against the
/// threads-only run (5e-9 Ha; observed far below), converge and verify on 2 and
/// 3 ranks, ranks agree bit for bit and a fixed rank count is bit-reproducible;
/// the stack without the grid switches matches it to 1e-12 Ha.
#[test]
fn combined_mpi_stack_keeps_the_standard_class() {
    let _lock = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let pp = common::potential(4., 1., true, 401);
    let stack = Arm {
        incremental: true,
        kept_metric: true,
        mirror: true,
        gram: true,
        reuse: true,
        lazy: true,
    };
    let solver_only = Arm {
        incremental: false,
        kept_metric: false,
        ..stack
    };
    for kind in 0..2 {
        let reference = results(kind, 1, Arm::default(), &pp).remove(0);
        for ranks in [2usize, 3] {
            let all = results(kind, ranks, stack, &pp);
            let again = results(kind, ranks, stack, &pp);
            let solver = results(kind, ranks, solver_only, &pp);
            for r in 0..ranks {
                let de = (all[r].total_hartree() - reference.total_hartree()).abs();
                let ds = (all[r].total_hartree() - solver[r].total_hartree()).abs();
                println!(
                    "kind {kind} P={ranks} rank {r}: stack dE vs threads-only {de:e} Ha, vs solver-only stack {ds:e} Ha, work {:?} (threads-only {:?})",
                    work(&all[r]),
                    work(&reference)
                );
                assert!(verified(&all[r]) && verified(&solver[r]));
                assert!(de < 5e-9, "standard class: dE {de:e} Ha");
                assert!(ds <= 1e-12, "grid switches moved the stack: {ds:e} Ha");
                assert_eq!(
                    all[r].total_hartree().to_bits(),
                    all[0].total_hartree().to_bits(),
                    "ranks disagree"
                );
                assert_eq!(density_bits(&all[r]), density_bits(&again[r]));
                assert_eq!(
                    all[r].total_hartree().to_bits(),
                    again[r].total_hartree().to_bits()
                );
            }
        }
    }
}
