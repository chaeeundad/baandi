//! G-HANDOFF-OWNED pricing gate: profiling cannot change the delivered state.
//! The default-off full route is the oracle for the owned transfer and profiling.
//! Broader simulated-rank matrix (debug budget):
//! cargo test --release -p dft-engine --test owned_handoff -- --include-ignored
#[allow(dead_code)]
mod common;
use dft_engine::{ScfOptions, ScfProblem, SpeciesPotential, cutoff_ladder, parallel, profile};
use dft_geometry::{Atom, Cell};
use dft_mpi::{Collective, SimRank, simulate};
use std::sync::{Arc, Mutex};

static LOCK: Mutex<()> = Mutex::new(());

fn compare(ranks: usize, skew: bool, terminal: bool) {
    compare_case(ranks, skew, terminal, 0);
}

fn compare_case(ranks: usize, skew: bool, terminal: bool, fallback: usize) {
    let pp = common::potential(2., 1., true, 201);
    let cell = Cell::from_user_rows_bohr(if skew {
        [[6., 0., 0.], [0.4, 7., 0.], [0.2, 0.3, 6.]]
    } else {
        [[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]]
    })
    .unwrap();
    let atoms = [Atom::new("Xx", [0.17, 0.23, 0.31])];
    let run = |enabled, owned| {
        profile::force(Some(enabled));
        simulate(ranks, |comm: &SimRank| {
            let config = dft_numerics::switches::RunConfig::from_pairs(&[(
                "m2_owned_handoff",
                if owned { "on" } else { "off" },
            )])
            .unwrap();
            let _config = config.enter();
            parallel::set_m2_override(Some(true));
            dft_engine::gamma_real::set_m2_lazy_override((fallback == 2).then_some(false));
            let shared: Arc<dyn Collective + Send + Sync> = Arc::new(comm.clone());
            parallel::set_communicator(Some(shared));
            let potentials = [SpeciesPotential {
                label: "Xx",
                potential: &pp,
            }];
            let options = ScfOptions {
                cutoff_hartree: 2.,
                bands: Some(2),
                max_iterations: if terminal { 1 } else { 80 },
                grid_override: Some(if skew { [15, 17, 15] } else { [16, 16, 16] }),
                ..if fallback == 3 {
                    ScfOptions::strict_reference()
                } else {
                    ScfOptions::default()
                }
            };
            let mut problem = ScfProblem::new(&cell, &atoms, &potentials, options).unwrap();
            if fallback == 1 {
                problem.decline_owned(dft_engine::grid_dist::OwnedDecline::Checkpoint);
            }
            let ladder = cutoff_ladder::CutoffLadderOptions {
                standard: true,
                coarse_max_iterations: 1,
                ..cutoff_ladder::CutoffLadderOptions::default()
            };
            let outcome = cutoff_ladder::run_ladder(problem, &potentials, &ladder, true).unwrap();
            assert!(outcome.report.used, "{:?}", outcome.report.fallback_reason);
            if !terminal {
                assert!(outcome.outcome.converged, "fine fixture must converge");
            }
            let mut result = outcome.outcome.result;
            if !terminal && fallback != 3 {
                assert!(result.verification.as_ref().is_some_and(|v| v.passed()));
            }
            let measured = cutoff_ladder::handoff_profile();
            result.materialize_orbitals().unwrap();
            let mut bits = vec![result.total_hartree().to_bits()];
            bits.extend(result.density.iter().map(|x| x.to_bits()));
            for k in &result.k_points {
                bits.extend(k.solution.eigenvalues.iter().map(|x| x.to_bits()));
                bits.extend(
                    k.solution
                        .eigenvectors
                        .as_slice()
                        .iter()
                        .flat_map(|z| [z.re.to_bits(), z.im.to_bits()]),
                );
            }
            bits.extend(
                result
                    .history
                    .iter()
                    .map(|r| r.operator_applications as u64),
            );
            parallel::set_communicator(None);
            parallel::set_m2_override(None);
            dft_engine::gamma_real::set_m2_lazy_override(None);
            (bits, measured)
        })
    };
    let off = run(false, false);
    let on = run(true, false);
    let direct = run(true, true);
    let unprofiled_direct = run(false, true);
    for (a, b) in direct.iter().zip(&unprofiled_direct) {
        assert_eq!(a.0, b.0, "owned setup depends on profiling");
        assert_eq!(b.1.collectives, [0; 5]);
    }
    profile::force(None);
    for (a, b) in on.iter().zip(&direct) {
        assert_eq!(
            a.0, b.0,
            "owned handoff changes fine result P={ranks} skew={skew} terminal={terminal}"
        );
        if ranks > 1 && !terminal && fallback == 0 {
            assert_eq!(b.1.owned[0], 1);
            assert!(
                b.1.collectives.iter().sum::<u64>() <= 2,
                "whole handoff P={ranks}: {:?}",
                b.1.collectives
            );
            assert_eq!(b.1.collectives[1], 0, "no handoff tree reductions");
            assert_eq!(b.1.materializations, 0);
            assert_eq!(b.1.gathered_bytes, 0);
            assert!(b.1.owned[1] > 0);
        } else {
            assert_eq!(b.1.owned[0], 0, "unsupported/cached consumer fallback");
        }
    }
    for (a, b) in off.iter().zip(&on) {
        assert_eq!(a.0, b.0, "profile changed fine result, P={ranks}");
        assert_eq!(b.0, on[0].0, "rank disagreement");
        assert_eq!(a.1.seconds(), 0.);
        assert_eq!(a.1.collectives, [0; 5]);
        assert_eq!(a.1.full_matrix_bytes, 0);
        assert!(b.1.seconds() > 0.);
        assert!(b.1.full_matrix_bytes > 0);
        assert!(b.1.coefficient_peak_bound_bytes >= b.1.full_matrix_bytes);
        if ranks > 1 {
            // A terminal result may already have been materialized by the
            // retained-result audit inside its final map, outside the handoff.
            assert_eq!(
                b.1.materializations,
                u64::from(!terminal && fallback != 2 && fallback != 3)
            );
            if b.1.materializations == 1 {
                assert!(b.1.gathered_bytes > 0);
                assert!(b.1.collectives[0] >= 3);
            }
        }
    }
}

#[test]
fn profile_is_observational_at_active_and_terminal_handoffs() {
    let _lock = LOCK.lock().unwrap();
    compare(2, false, false);
    compare(2, true, true);
}

#[test]
#[ignore = "release-only simulated rank/cell/grid matrix"]
fn profile_matrix_one_through_seven_ranks() {
    let _lock = LOCK.lock().unwrap();
    for ranks in 1..=7 {
        for skew in [false, true] {
            compare(ranks, skew, false);
        }
    }
}

#[test]
fn whole_handoff_uses_at_most_two_collectives_on_two_through_seven_ranks() {
    let _lock = LOCK.lock().unwrap();
    for ranks in 2..=7 {
        compare(ranks, false, false);
    }
}

#[test]
fn unsupported_consumers_keep_the_full_route() {
    let _lock = LOCK.lock().unwrap();
    for fallback in 1..=3 {
        compare_case(3, true, false, fallback);
    }
}

#[test]
fn rank_local_consumer_mismatch_is_an_agreed_error() {
    let _lock = LOCK.lock().unwrap();
    let pp = common::potential(2., 1., true, 201);
    for ranks in 2..=7 {
        let errors = simulate(ranks, |comm: &SimRank| {
            let config =
                dft_numerics::switches::RunConfig::from_pairs(&[("m2_owned_handoff", "on")])
                    .unwrap();
            let _config = config.enter();
            parallel::set_m2_override(Some(true));
            parallel::set_communicator(Some(Arc::new(comm.clone())));
            let cell =
                Cell::from_user_rows_bohr([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]]).unwrap();
            let potentials = [SpeciesPotential {
                label: "Xx",
                potential: &pp,
            }];
            let mut problem = ScfProblem::new(
                &cell,
                &[Atom::new("Xx", [0.1, 0.2, 0.3])],
                &potentials,
                ScfOptions {
                    cutoff_hartree: 2.,
                    bands: Some(2),
                    ..ScfOptions::default()
                },
            )
            .unwrap();
            if comm.rank() == 0 {
                problem.decline_owned(dft_engine::grid_dist::OwnedDecline::Checkpoint);
            }
            let result = cutoff_ladder::run_ladder(
                problem,
                &potentials,
                &cutoff_ladder::CutoffLadderOptions {
                    standard: true,
                    ..Default::default()
                },
                true,
            );
            let error = match result {
                Ok(_) => panic!("rank-local eligibility must reject"),
                Err(e) => e.to_string(),
            };
            parallel::set_communicator(None);
            parallel::set_m2_override(None);
            error
        });
        assert!(
            errors
                .iter()
                .all(|e| e == &errors[0] && e.contains("eligibility differs"))
        );
    }
}
