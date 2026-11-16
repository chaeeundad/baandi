//! F21 (`BAANDI_STAGE_TIMERS=on|off`, default on): the MPI stage and collective
//! counters in count-only mode count exactly the calls and payload bytes of the
//! time mode, read no clock (all seconds stay 0) and leave the SCF result bit for
//! bit unchanged. Simulated ranks (two threads) on a one-atom Gaussian fixture
//! with nonlocal projectors. One test function: the mode and counters are
//! process-wide.
#[allow(dead_code)]
#[path = "../../dft-cli/tests/common/mod.rs"]
mod common;
use dft_engine::{GammaReal, ScfOptions, ScfProblem, ScfResult, SpeciesPotential, parallel};
use dft_geometry::{Atom, Cell};
use dft_mpi::{Collective, SimRank, simulate};
use std::sync::Arc;

fn run(comm: Arc<dyn Collective + Send + Sync>) -> ScfResult {
    parallel::set_m2_override(Some(true));
    parallel::set_communicator(Some(comm));
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0., 5.2, 0.], [0.3, 0., 5.4]]).unwrap();
    let pp = dft_upf::parse_upf(&common::gaussian_with(8., 1., true, 401)).unwrap();
    let result = ScfProblem::new(
        &cell,
        &[Atom::new("Xx", [0.13, 0.21, 0.32])],
        &[SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }],
        ScfOptions {
            gamma_real: GammaReal::Auto,
            cutoff_hartree: 3.,
            bands: Some(6),
            threads: 1,
            eigen_tolerance: 1e-10,
            ..ScfOptions::strict_reference()
        },
    )
    .unwrap()
    .run()
    .unwrap();
    parallel::set_communicator(None);
    parallel::set_m2_override(None);
    result
}

/// Process counters: (stage calls, stage seconds, collective calls, seconds,
/// bytes, M2 counters incl. `calls.*`, M2 stage seconds).
type Snap = (
    Vec<u64>,
    Vec<f64>,
    Vec<usize>,
    Vec<f64>,
    Vec<u64>,
    Vec<usize>,
    Vec<f64>,
);
fn snapshot() -> Snap {
    let (counts, stages) = parallel::m2_stats();
    (
        parallel::stage_calls().into_iter().map(|x| x.1).collect(),
        parallel::stage_seconds().into_iter().map(|x| x.1).collect(),
        parallel::collective_stats()
            .into_iter()
            .map(|x| x.1)
            .collect(),
        parallel::collective_stats()
            .into_iter()
            .map(|x| x.2)
            .collect(),
        parallel::collective_bytes()
            .into_iter()
            .map(|x| x.1)
            .collect(),
        counts.into_iter().map(|x| x.1).collect(),
        stages.into_iter().map(|x| x.1).collect(),
    )
}
fn du(a: &[u64], b: &[u64]) -> Vec<u64> {
    a.iter().zip(b).map(|(x, y)| y - x).collect()
}
fn dz(a: &[usize], b: &[usize]) -> Vec<usize> {
    a.iter().zip(b).map(|(x, y)| y - x).collect()
}
fn df(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b).map(|(x, y)| y - x).collect()
}

#[test]
fn count_only_mode_counts_the_same_and_reads_no_clock() {
    let arm = |timed: bool| {
        dft_numerics::stage_clock::force(Some(timed));
        let before = snapshot();
        let results = simulate(2, |comm: &SimRank| run(Arc::new(comm.clone())));
        let after = snapshot();
        dft_numerics::stage_clock::force(None);
        (results, before, after)
    };
    let (off, b0, a0) = arm(false);
    let (on, b1, a1) = arm(true);
    // Same work in both modes: stage entries, collective calls and bytes, M2 counters.
    assert_eq!(du(&b0.0, &a0.0), du(&b1.0, &a1.0), "solver stage calls");
    assert_eq!(dz(&b0.2, &a0.2), dz(&b1.2, &a1.2), "collective calls");
    assert_eq!(du(&b0.4, &a0.4), du(&b1.4, &a1.4), "collective bytes");
    assert_eq!(dz(&b0.5, &a0.5), dz(&b1.5, &a1.5), "M2 counters");
    assert!(du(&b0.0, &a0.0).iter().sum::<u64>() > 0, "no stage counted");
    assert!(dz(&b0.2, &a0.2).iter().sum::<usize>() > 0, "no collective");
    assert!(du(&b0.4, &a0.4).iter().sum::<u64>() > 0, "no payload bytes");
    // Count-only: no time anywhere; time mode: some time.
    for d in [df(&b0.1, &a0.1), df(&b0.3, &a0.3), df(&b0.6, &a0.6)] {
        assert!(
            d.iter().all(|&s| s == 0.),
            "count-only mode read a clock: {d:?}"
        );
    }
    assert!(df(&b1.1, &a1.1).iter().sum::<f64>() > 0.);
    assert!(df(&b1.3, &a1.3).iter().sum::<f64>() > 0.);
    // The result does not depend on the mode.
    for (x, y) in off.iter().zip(&on) {
        assert_eq!(x.total_hartree().to_bits(), y.total_hartree().to_bits());
        assert_eq!(x.history.len(), y.history.len());
        assert!(
            x.density
                .iter()
                .zip(&y.density)
                .all(|(a, b)| a.to_bits() == b.to_bits())
        );
    }
    println!(
        "collective calls {:?}, payload bytes {:?}",
        dz(&b0.2, &a0.2),
        du(&b0.4, &a0.4)
    );
}
