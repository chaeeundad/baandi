//! Reproducible PACKED_REDUCE scientific identity oracle. Copy temporarily to
//! crates/dft-engine/tests/lane_packed_oracle.rs to resolve existing fixtures.
//! Mac output times are discarded; the exported trace tests counter transport only.
#[allow(dead_code)]
mod common;
use dft_engine::{GammaReal, ScfOptions, ScfProblem, ScfResult, SpeciesPotential, parallel};
use dft_geometry::{Atom, Cell};
use dft_mpi::{Collective, SimRank, simulate};
use dft_solvers::distributed_davidson::{
    DistributedTuning, set_start_gram_override, set_tuning_override,
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

fn results(kind: usize, ranks: usize, a: Arm, packed: bool, pp: &dft_upf::Pseudopotential) -> Vec<ScfResult> {
    if ranks == 1 {
        return vec![run(kind, pp)];
    }
    simulate(ranks, |comm: &SimRank| {
        let config = dft_numerics::switches::RunConfig::from_pairs(&[("davidson_packed_reduce", if packed { "on" } else { "off" }), ("davidson_profile", "on"), ("davidson_upper_border", "off"), ("davidson_rotate_columns", "off")]).unwrap();
        let _scope = config.enter();
        dft_engine::profile::davidson_lane::reset();
        use dft_solvers::distributed_fft::set_mirror_sticks_override;
        parallel::set_m2_override(Some(true));
        set_tuning_override(Some(DistributedTuning {
            lean_algebra: a.lean,
            ..DistributedTuning::from_env()
        }));
        set_mirror_sticks_override(Some(a.stack));
        set_start_gram_override(Some(a.stack));
        dft_engine::gamma_real::set_m2_reuse_override(Some(a.stack));
        let shared: Arc<dyn Collective + Send + Sync> = Arc::new(comm.clone());
        parallel::set_communicator(Some(shared));
        let r = run(kind, pp);
        use dft_engine::profile::davidson_lane as lane;
        let trace = lane::take(comm.rank());
        let traces = lane::gather(comm, trace).unwrap();
        if comm.rank() == 0 {
            if let Ok(dir) = std::env::var("PACKED_ORACLE_OUT") {
                let rows: Vec<String> = traces.iter().map(|t| format!(
                    "{{\"rank\":{},\"dropped\":{},\"records\":{:?}}}",
                    t.rank, t.dropped, t.records)).collect();
                let json = format!("{{\"run\":{{\"completed\":true}},\"davidson_profile\":{{\"schema\":\"davidson-profile-1\",\"fields\":{:?},\"operations\":{:?},\"stages\":{:?},\"kinds\":{:?},\"reduction_sites\":{:?},\"per_rank\":[{}]}}}}",
                    lane::FIELDS, lane::OPERATIONS, lane::STAGES, lane::KINDS,
                    lane::REDUCTION_SITES, rows.join(","));
                std::fs::write(format!("{dir}/packed-{kind}-np{ranks}-{packed}.json"), json).unwrap();
            }
        }
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
fn packed_reduction_full_scf_identity() {
    let _lock = LOCK.lock().unwrap();
    let pp = common::potential(4., 1., true, 401);
    for kind in 0..2 {
        for ranks in [2, 3] {
            let arm = Arm { lean: true, stack: false };
            let off = results(kind, ranks, arm, false, &pp);
            let on = results(kind, ranks, arm, true, &pp);
            for (a, b) in off.iter().zip(&on) {
                assert!(verified(a) && verified(b));
                assert_eq!(a.total_hartree().to_bits(), b.total_hartree().to_bits());
                assert_eq!(density_bits(a), density_bits(b));
                assert_eq!(work(a), work(b));
                assert_eq!(science_debug(&a.history), science_debug(&b.history));
                assert_eq!(science_debug(&a.verification), science_debug(&b.verification));
                assert_eq!(a.total_hartree().to_bits(), off[0].total_hartree().to_bits());
            }
            println!("packed kind={kind} ranks={ranks}: energy/density/history/verification bit-identical; work={:?}", work(&on[0]));
        }
    }
}

// Scientific debug records include wall-clock metadata. Exclude only those
// PhaseTimings blocks; every scientific field remains compared verbatim.
fn science_debug(value: &impl std::fmt::Debug) -> String {
    let mut text = format!("{value:?}");
    while let Some(start) = text.find("timings: PhaseTimings {") {
        let end = start + text[start..].find('}').unwrap() + 1;
        text.replace_range(start..end, "");
    }
    while let Some(start) = text.find("wall: ") {
        let end = start + text[start..].find(',').unwrap() + 1;
        text.replace_range(start..end, "");
    }
    text
}
