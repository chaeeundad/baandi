//! BREAKTHROUGH-4-OTHER-AND-INIT: the switches that remove replicated work
//! found by the finer phase-profile labels are bit-identical, and the new
//! labels classify the solve wall.
//!
//! - `BAANDI_M2_LIGHT_STATE` (per-thread override): the serial real Gamma
//!   Hamiltonian of an M2 solve state is built with the minimal capacity.
//! - `BAANDI_M2_LAYOUT_CACHE` (per-thread override): the distributed solver
//!   reuses the row-index layout of its previous solve.
//! - `BAANDI_EWALD_CACHE` (per-thread override): a construction with bit-equal
//!   lattice and ions reuses the Ewald energy.
//!
//! Production standard SCF (fixed-point verification) on two simulated ranks
//! (M2 stick rows, kept start, lazy orbitals as by default). Own test binary
//! with a lock (process-wide forced switches and counters).
#[allow(dead_code)]
mod common;
use dft_engine::{GammaReal, ScfOptions, ScfProblem, ScfResult, SpeciesPotential, parallel};
use dft_geometry::{Atom, Cell};
use dft_mpi::{Collective, SimRank, simulate};
use std::sync::{Arc, Mutex};

static LOCK: Mutex<()> = Mutex::new(());

fn problem(pp: &dft_upf::Pseudopotential, skew: bool) -> ScfProblem {
    let rows = if skew {
        [[6.2, 0., 0.], [1.1, 5.7, 0.], [0.4, 0.7, 6.6]]
    } else {
        [[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]]
    };
    let cell = Cell::from_user_rows_bohr(rows).unwrap();
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

#[derive(Clone, Copy)]
struct Arm {
    light: bool,
    layout: bool,
}

fn results(ranks: usize, pp: &dft_upf::Pseudopotential, skew: bool, arm: Arm) -> Vec<ScfResult> {
    simulate(ranks, |comm: &SimRank| {
        parallel::set_m2_override(Some(true));

        dft_engine::gamma_real::set_m2_light_state_override(Some(arm.light));
        dft_solvers::distributed_davidson::set_layout_cache_override(Some(arm.layout));
        let shared: Arc<dyn Collective + Send + Sync> = Arc::new(comm.clone());
        parallel::set_communicator(Some(shared));
        let r = problem(pp, skew).run().unwrap();
        dft_engine::profile::finish_rank(comm.rank());
        parallel::set_communicator(None);
        dft_solvers::distributed_davidson::set_layout_cache_override(None);
        dft_engine::gamma_real::set_m2_light_state_override(None);

        parallel::set_m2_override(None);
        r
    })
}

fn bits(r: &ScfResult) -> (u64, Vec<u64>, Vec<u64>) {
    (
        r.total_hartree().to_bits(),
        r.density.iter().map(|x| x.to_bits()).collect(),
        r.k_points[0]
            .solution
            .eigenvalues
            .iter()
            .map(|x| x.to_bits())
            .collect(),
    )
}

#[test]
fn light_state_and_layout_cache_are_bit_identical() {
    let _lock = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let pp = common::potential(4., 1., true, 401);
    for skew in [false, true] {
        let plain = results(
            2,
            &pp,
            skew,
            Arm {
                light: false,
                layout: false,
            },
        );
        let (reused_before, _) = dft_solvers::distributed_davidson::layout_cache_stats();
        let switched = results(
            2,
            &pp,
            skew,
            Arm {
                light: true,
                layout: true,
            },
        );
        let (reused_after, stored) = dft_solvers::distributed_davidson::layout_cache_stats();
        let _ = dft_engine::profile::take_registry();
        for r in 0..2 {
            assert_eq!(
                bits(&plain[r]),
                bits(&switched[r]),
                "skew {skew} rank {r}: results changed"
            );
            assert!(switched[r].certified);
            assert_eq!(plain[r].history.len(), switched[r].history.len());
        }
        assert_eq!(bits(&switched[0]), bits(&switched[1]), "ranks agree");
        assert!(stored > 0, "layout stored");
        assert!(
            reused_after > reused_before,
            "layout reused by later solves"
        );
    }
}

#[test]
fn ewald_cache_reuses_bit_equal_constructions_only() {
    let _lock = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let pp = common::potential(4., 1., true, 401);
    dft_engine::set_ewald_cache_override(Some(false));
    let reference = problem(&pp, true).run().unwrap();
    let skew_reference = problem(&pp, false).run().unwrap();
    dft_engine::set_ewald_cache_override(Some(true));
    let [hits0, _] = dft_engine::ewald_cache_counters();
    let first = problem(&pp, true).run().unwrap();
    let second = problem(&pp, true).run().unwrap();
    let [hits1, _] = dft_engine::ewald_cache_counters();
    // Another cell: computed, not taken from the cache.
    let other = problem(&pp, false).run().unwrap();
    let [hits2, _] = dft_engine::ewald_cache_counters();
    dft_engine::set_ewald_cache_override(None);
    assert_eq!(hits1 - hits0, 1, "the second identical construction hits");
    assert_eq!(hits2, hits1, "a different cell misses");
    for r in [&first, &second] {
        assert_eq!(
            r.total_hartree().to_bits(),
            reference.total_hartree().to_bits()
        );
        assert_eq!(
            r.ewald_breakdown.energy.0.to_bits(),
            reference.ewald_breakdown.energy.0.to_bits()
        );
    }
    assert_eq!(
        other.total_hartree().to_bits(),
        skew_reference.total_hartree().to_bits()
    );
}

#[test]
fn finer_labels_classify_the_m2_solve() {
    let _lock = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let pp = common::potential(4., 1., true, 401);
    dft_engine::profile::force(Some(true));
    let _ = dft_engine::profile::take_registry();
    let _ = results(
        2,
        &pp,
        false,
        Arm {
            light: false,
            layout: false,
        },
    );
    dft_engine::profile::force(Some(false));
    let profiles = dft_engine::profile::take_registry();
    assert_eq!(profiles.len(), 2);
    let report = dft_engine::profile::Report::from_ranks(&profiles);
    print!("{}", report.table(None));
    let entries = report.entries();
    let fine = &entries
        .iter()
        .find(|(b, _)| *b == "fine")
        .expect("fine bucket")
        .1;
    let get = |label: &str| {
        fine.iter()
            .find(|(l, _, _)| *l == label)
            .unwrap_or_else(|| panic!("missing {label}"))
            .1
    };
    let solve = get("solve");
    assert!(solve > 0.0);
    assert!(get("solve.solver_total") > 0.0, "solver container timed");
    assert!(get("solve.real_state") > 0.0, "real state labelled");
    assert!(get("solve.setup") > 0.0, "solver setup labelled");
    assert_eq!(get("solve.kept_update"), 0.0, "strict does not keep H X");
    // Exclusive split: per rank other = other_solver + other_engine; the
    // report takes each maximum over ranks separately, so only bounds hold.
    let other = get("solve.other");
    let tiny = 1e-12 * solve.max(1.0) + 1e-15;
    assert!(get("solve.other_solver") <= other + tiny);
    assert!(get("solve.other_engine") <= other + tiny);
    assert!(other <= get("solve.other_solver") + get("solve.other_engine") + tiny);
    assert!(other <= solve);
    // Init parts of both constructions are recorded (none on a missing clock).
    let init = &entries.iter().find(|(b, _)| *b == "init").unwrap().1;
    let part = |label: &str| init.iter().find(|(l, _, _)| *l == label).unwrap().1;
    assert!(part("init.ewald") > 0.0 && part("init.states") > 0.0);
}
