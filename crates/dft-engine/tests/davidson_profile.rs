//! Lane profile transport and scientific identity on simulated MPI ranks.
#[allow(dead_code)]
#[path = "../../dft-cli/tests/common/mod.rs"]
mod common;
use dft_engine::profile::davidson_lane as lane;
use dft_engine::{GammaReal, ScfOptions, ScfProblem, SpeciesPotential, parallel};
use dft_geometry::{Atom, Cell};
use dft_mpi::{Collective, simulate};
use dft_numerics::switches::RunConfig;
use std::sync::Arc;

#[test]
fn lane_trace_is_bit_identical_and_rank_local() {
    let run = |on| {
        simulate(3, |comm| {
            let config =
                RunConfig::from_pairs(&[("davidson_profile", if on { "on" } else { "off" })])
                    .unwrap();
            let _scope = config.enter();
            lane::reset();
            parallel::set_m2_override(Some(true));
            parallel::set_communicator(Some(Arc::new(comm.clone())));
            let cell =
                Cell::from_user_rows_bohr([[5., 0., 0.], [0., 5.2, 0.], [0.3, 0., 5.4]]).unwrap();
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
            let trace = lane::take(comm.rank());
            let trace = lane::gather(comm, trace).unwrap();
            parallel::set_communicator(None);
            parallel::set_m2_override(None);
            (result, trace)
        })
    };
    let off = run(false);
    let on = run(true);
    for ((a, ta), (b, tb)) in off.iter().zip(&on) {
        assert_eq!(a.total_hartree().to_bits(), b.total_hartree().to_bits());
        assert_eq!(a.history.len(), b.history.len());
        assert!(
            a.density
                .iter()
                .zip(&b.density)
                .all(|(a, b)| a.to_bits() == b.to_bits())
        );
        assert!(ta.iter().all(|r| r.records.is_empty() && r.dropped == 0));
        assert_eq!(tb.len(), 3);
        for (rank, trace) in tb.iter().enumerate() {
            assert_eq!(trace.rank, rank);
            assert_eq!(trace.dropped, 0);
            assert!(!trace.records.is_empty());
            assert!(
                trace
                    .records
                    .iter()
                    .any(|r| r[6] == lane::Op::TreeSum as usize as f64)
            );
            assert!(trace.records.iter().all(|r| r[1] > 0.0 && r[2] > 0.0));
        }
    }
}
