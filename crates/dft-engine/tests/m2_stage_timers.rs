//! M2-STAGE-TIMERS (2026-11-18): `BAANDI_DAVIDSON_STAGES=on` fills the serial
//! Davidson stage buckets and the h_apply parts on the MPI M2 path (distributed
//! solver and operator), with collective counts inside H applications, and
//! never changes the arithmetic. Two simulated ranks (threads) on a one-atom
//! Gaussian fixture with nonlocal projectors. One test function: the timer
//! switches and counters are process-wide.
#[allow(dead_code)]
#[path = "../../dft-cli/tests/common/mod.rs"]
mod common;
use dft_engine::{GammaReal, ScfOptions, ScfProblem, ScfResult, SpeciesPotential, parallel};
use dft_geometry::{Atom, Cell};
use dft_mpi::{Collective, SimRank, simulate};
use dft_numerics::work_timers::{self, WorkTimes};
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

fn pooled() -> WorkTimes {
    work_timers::snapshot()
        .iter()
        .fold(WorkTimes::default(), |a, b| a.plus(b))
}
/// (h_apply, reduced, rotate, orthonormalize seconds, timed solves) summed
/// over the run stages.
fn stages() -> (f64, f64, f64, f64, u64) {
    let s = dft_solvers::davidson_stages::stage_seconds();
    (
        s.h_apply,
        s.reduced,
        s.rotate,
        s.orthonormalize,
        s.timed_solves,
    )
}

#[test]
fn m2_path_fills_stage_and_h_apply_timers_without_changing_results() {
    dft_numerics::stage_clock::force(Some(true));
    // Timers off: nothing recorded, no clock.
    work_timers::set_enabled(false);
    let (w0, s0) = (pooled(), stages());
    let off = simulate(2, |comm: &SimRank| run(Arc::new(comm.clone())));
    assert_eq!(pooled(), w0, "work timers recorded with the switch off");
    assert_eq!(stages(), s0, "stage timers recorded with the switch off");

    // Timers on.
    work_timers::set_enabled(true);
    let on = simulate(2, |comm: &SimRank| run(Arc::new(comm.clone())));
    let (w, s1) = (pooled(), stages());
    work_timers::set_enabled(false);
    dft_numerics::stage_clock::force(None);

    // Serial stage buckets of the distributed solver.
    assert!(s1.0 > 0. && s1.1 > 0. && s1.2 > 0. && s1.3 > 0., "{s1:?}");
    assert!(s1.4 >= 2, "timed solves {}", s1.4);
    // h_apply parts: all positive, summing to the whole call.
    assert!(w.apply_calls > 0 && w.apply_columns >= w.apply_calls);
    assert!(w.apply > 0. && w.local > 0., "{w:?}");
    assert!(w.nonlocal_projection > 0. && w.nonlocal_back > 0., "{w:?}");
    let parts = w.local + w.nonlocal_projection + w.nonlocal_back;
    let remainder = w.apply - parts;
    assert!(
        remainder >= -1e-9 && remainder <= 0.05 * w.apply + 1e-4,
        "apply {} parts {parts} remainder {remainder}",
        w.apply
    );
    // Distributed FFT pairs and density transforms, collectives inside H.
    assert!(w.pair_transforms > 0 && w.pair_columns > 0 && w.pair_unsplit > 0.);
    assert!(w.density_transforms > 0 && w.density > 0.);
    let (alltoall, tree_sum) = (3, 1);
    assert!(w.apply_coll_calls[alltoall] > 0 && w.apply_coll[alltoall] > 0.);
    assert!(w.apply_coll_calls[tree_sum] > 0, "projector tree sum");
    for i in 0..work_timers::COLLECTIVES {
        assert!(w.apply_coll_calls[i] <= w.coll_calls[i]);
    }
    println!(
        "h_apply: calls {} cols {} apply {:.3e} s local {:.3e} proj {:.3e} back {:.3e} rem {:.3e}; \
         pairs {} fft {:.3e} s; in-apply alltoall {} tree_sum {}; stages {s1:?}",
        w.apply_calls,
        w.apply_columns,
        w.apply,
        w.local,
        w.nonlocal_projection,
        w.nonlocal_back,
        remainder,
        w.pair_transforms,
        w.pair_unsplit,
        w.apply_coll_calls[alltoall],
        w.apply_coll_calls[tree_sum],
    );

    // Timers never change the arithmetic.
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
}
