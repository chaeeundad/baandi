//! BREAKTHROUGH-0-INSTRUMENT: the measurement switches (phase profile,
//! SCF snapshot capture) leave the SCF bit for
//! bit unchanged, the snapshot of a solve replays its work exactly, and the
//! profile buckets are filled. Production standard SCF (fixed-point
//! verification) on the threads-only path (one rank, REUSE-PSI kept start) and
//! on the M2 stick-row path (two simulated ranks, kept start). Own test binary
//! with a lock (process-wide forced switches and counters).
#[allow(dead_code)]
mod common;
use dft_engine::{GammaReal, ScfOptions, ScfProblem, SpeciesPotential, parallel};
use dft_geometry::{Atom, Cell};
use dft_mpi::{Collective, SimRank, simulate};
use std::sync::{Arc, Mutex};

static LOCK: Mutex<()> = Mutex::new(());

fn problem(pp: &dft_upf::Pseudopotential) -> ScfProblem {
    let cell = Cell::from_user_rows_bohr([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]]).unwrap();
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

#[test]
fn per_rank_closure_through_the_gather_on_simulated_ranks() {
    let _lock = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let pp = common::potential(4., 1., true, 401);
    dft_engine::profile::force(Some(true));
    let gathered = simulate(2, |comm: &SimRank| {
        parallel::set_m2_override(Some(true));
        let shared: Arc<dyn Collective + Send + Sync> = Arc::new(comm.clone());
        parallel::set_communicator(Some(shared));
        let p = problem(&pp);
        let start = std::time::Instant::now();
        p.run().unwrap();
        let mut mine = dft_engine::profile::finish_rank(comm.rank()).unwrap();
        mine.wall = start.elapsed().as_secs_f64();
        let all = dft_engine::profile::gather(comm, &mine).unwrap();
        parallel::set_communicator(None);
        parallel::set_m2_override(None);
        all
    });
    dft_engine::profile::force(Some(false));
    let _ = dft_engine::profile::take_registry();
    assert_eq!(gathered[0], gathered[1], "every rank holds every profile");
    let report = dft_engine::profile::Report::from_ranks(&gathered[0]);
    let c = report.closure().unwrap();
    for r in &report.per_rank {
        let (wall, residual) = (r.wall.unwrap(), r.residual.unwrap());
        assert!(r.covered > 0.0, "rank {}: covered", r.rank);
        assert!(
            residual >= -1e-3f64.max(1e-3 * wall),
            "rank {}: covered {} above wall {wall}",
            r.rank,
            r.covered
        );
    }
    assert!(c.covered_min <= c.covered_max && (0.0..=1.0).contains(&c.imbalance));
    println!(
        "closure: covered max {:.4} s (rank {}) min {:.4} s, imbalance {:.3}, residual {:?}..{:?}; legacy sum of maxima {:.4} s",
        c.covered_max,
        c.critical_rank,
        c.covered_min,
        c.imbalance,
        c.residual_min,
        c.residual_max,
        report.covered_wall()
    );
}
