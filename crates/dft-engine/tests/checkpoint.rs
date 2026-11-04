#[allow(dead_code)]
mod common;
#[path = "common/checkpoint.rs"]
mod fixture;
use dft_engine::{checkpoint::*, *};
use fixture::*;
fn start(kind: usize, mixing: Mixing, threads: usize, mode: ExecutionMode) -> ScfSession {
    let mut s = ScfSession::new(
        fixture(kind, mixing, threads, 200),
        CompletionMode::Strict,
        ExecutionConfig {
            mode,
            ..Default::default()
        },
    )
    .unwrap();
    s.initialize().unwrap();
    drain(&mut s);
    s
}
#[test]
fn named_compatibility_rejections() {
    let (m, e) = identity(0, 1);
    let s = start(0, Mixing::default(), 1, ExecutionMode::Classic);
    let cp = s
        .checkpoint(CheckpointKind::ExactContinuation, m.clone(), e.clone())
        .unwrap();
    macro_rules! rejects {
        ($field:ident, $value:expr, $why:ident) => {{
            let mut target = m.clone();
            target.$field = $value;
            assert_eq!(
                cp.compatible(CheckpointKind::ExactContinuation, &target, &e),
                Err(Rejection::$why)
            );
        }};
    }
    rejects!(cell_bohr, [[5., 0., 0.], [0., 4., 0.], [0., 0., 4.]], Cell);
    rejects!(atoms_fractional, vec![[0.1, 0.2, 0.3]], Atoms);
    rejects!(species, vec!["Yy".into()], Species);
    rejects!(
        pseudopotential_sha256,
        vec!["a".repeat(64)],
        Pseudopotential
    );
    rejects!(cutoff_hartree, 3.0, Cutoff);
    rejects!(kpoints, vec![[0.1, 0., 0., 1.]], KMesh);
    rejects!(electrons, 4.0, ElectronCount);
    rejects!(xc, "pbe".into(), Xc);
    rejects!(spin, "polarized".into(), Spin);
    rejects!(environment, "other".into(), Environment);
    let mut bad = cp.clone();
    bad.schema_version += 1;
    assert_eq!(
        bad.compatible(CheckpointKind::ExactContinuation, &m, &e),
        Err(Rejection::SchemaVersion)
    );
    for (target, reason) in [
        (
            ExecutionIdentity {
                threads: 4,
                ..e.clone()
            },
            Rejection::Threads,
        ),
        (
            ExecutionIdentity {
                build_sha256: "a".repeat(64),
                ..e.clone()
            },
            Rejection::Build,
        ),
        (
            ExecutionIdentity {
                platform: "other".into(),
                ..e.clone()
            },
            Rejection::Platform,
        ),
    ] {
        assert_eq!(
            cp.compatible(CheckpointKind::ExactContinuation, &m, &target),
            Err(reason)
        );
    }
    assert_eq!(
        cp.compatible(CheckpointKind::WarmStart, &m, &e),
        Err(Rejection::Kind)
    );
}
#[test]
fn invalid_state_rejected_without_touching_source() {
    let (m, e) = identity(0, 1);
    let s = start(0, Mixing::default(), 1, ExecutionMode::Classic);
    let cp = s
        .checkpoint(CheckpointKind::ExactContinuation, m.clone(), e.clone())
        .unwrap();
    let restore = |c: &Checkpoint| {
        ScfSession::restore_exact(fixture(0, Mixing::default(), 1, 200), c, &m, &e)
            .err()
            .unwrap()
    };
    let mut bad = cp.clone();
    bad.electronic_state.density[0] = f64::NAN;
    assert_eq!(restore(&bad), Rejection::Density);
    let mut bad = cp.clone();
    bad.electronic_state.density[0] += 1.;
    assert_eq!(restore(&bad), Rejection::Density);
    let mut bad = cp.clone();
    bad.electronic_state.grid[0] += 1;
    assert_eq!(restore(&bad), Rejection::Grid);
    let mut bad = cp.clone();
    bad.continuation.as_mut().unwrap().iteration = 1;
    bad.continuation.as_mut().unwrap().state_version = 2;
    assert_eq!(restore(&bad), Rejection::History);
    let mut bad = cp.clone();
    bad.continuation.as_mut().unwrap().options.mixing = Mixing::Linear { beta: 0.3 };
    assert_eq!(restore(&bad), Rejection::Options);
    let mut observed = start(0, Mixing::default(), 1, ExecutionMode::Observed);
    observed.step().unwrap();
    assert_eq!(
        observed
            .checkpoint(CheckpointKind::ExactContinuation, m.clone(), e.clone())
            .unwrap_err(),
        Rejection::UndrainedEvents
    );
    drain(&mut observed);
    assert!(
        observed
            .checkpoint(CheckpointKind::ExactContinuation, m, e)
            .is_ok()
    );
}
#[test]
#[ignore = "release matrix: cargo test --release -p dft-engine --test checkpoint -- --include-ignored --nocapture"]
fn exact_continuation_matrix() {
    for kind in 0..3 {
        for threads in [1, 4] {
            for mixing in [
                Mixing::Linear { beta: 0.3 },
                Mixing::default(),
                Mixing::KerkerPulay {
                    beta: 0.3,
                    history: 6,
                    q0_inverse_bohr: 1.0,
                },
            ] {
                let (m, e) = identity(kind, threads);
                let baseline = finish(start(kind, mixing.clone(), threads, ExecutionMode::Classic));
                for mode in [
                    ExecutionMode::Classic,
                    ExecutionMode::Observed,
                    ExecutionMode::AgentControlled,
                ] {
                    for cut in [0, 3, 8] {
                        let mut s = start(kind, mixing.clone(), threads, mode);
                        for _ in 0..cut {
                            s.step().unwrap();
                            drain(&mut s);
                        }
                        let cp = s
                            .checkpoint(CheckpointKind::ExactContinuation, m.clone(), e.clone())
                            .unwrap();
                        drop(s);
                        let restored = ScfSession::restore_exact(
                            fixture(kind, mixing.clone(), threads, 200),
                            &cp,
                            &m,
                            &e,
                        )
                        .unwrap();
                        assert_eq!(restored.snapshot().sequence, cp.sequence);
                        assert_identical(&baseline, &finish(restored));
                    }
                }
                println!(
                    "exact kind={kind} threads={threads} mixing={mixing:?} energy={:.17} Ha iterations={}",
                    baseline.total_hartree(),
                    baseline.history.len()
                );
            }
        }
    }
}
#[test]
#[ignore = "release convergence test"]
fn warm_start_different_mixing() {
    for kind in 0..3 {
        for threads in [1, 4] {
            let (m, e) = identity(kind, threads);
            let baseline = finish(start(
                kind,
                Mixing::default(),
                threads,
                ExecutionMode::Classic,
            ));
            let mut s = start(
                kind,
                Mixing::Linear { beta: 0.3 },
                threads,
                ExecutionMode::Classic,
            );
            for _ in 0..3 {
                s.step().unwrap();
            }
            let cp = s
                .checkpoint(CheckpointKind::WarmStart, m.clone(), e.clone())
                .unwrap();
            let restored = ScfSession::warm_start(
                fixture(kind, Mixing::default(), threads, 200),
                &cp,
                &m,
                &e,
                CompletionMode::Strict,
                ExecutionConfig {
                    calc_id: "new-calculation".into(),
                    ..Default::default()
                },
            )
            .unwrap();
            assert!(restored.history().is_empty());
            assert_eq!(restored.snapshot().iteration, 0);
            let result = finish(restored);
            let delta = (result.total_hartree() - baseline.total_hartree()).abs();
            println!("warm kind={kind} threads={threads} delta={delta:.17e} Ha");
            assert!(delta < 1e-10);
        }
    }
}

#[test]
fn runtime_boundaries_and_consumed_autonomy_survive() {
    use dft_engine::policy::PolicyChange;
    use dft_engine::runtime::*;
    let (m, e) = identity(0, 1);
    let session = ScfSession::new(
        fixture(0, Mixing::default(), 1, 200),
        CompletionMode::Strict,
        ExecutionConfig {
            mode: ExecutionMode::AgentControlled,
            ..Default::default()
        },
    )
    .unwrap();
    let config = RuntimeConfig {
        autonomy: Some(AutonomyPolicy {
            max_applied_total: 1,
            rules: vec![AutonomyRule {
                action: AllowedAction::Beta {
                    minimum: 0.1,
                    maximum: 0.3,
                    only_decrease: true,
                },
                max_count: 1,
                cooldown_iterations: 0,
                max_version_lag: 0,
                allow_revalidation: false,
            }],
        }),
        ..Default::default()
    };
    let mut runtime = Runtime::new(session, config).unwrap();
    runtime.initialize().unwrap();
    while runtime.events().unwrap().pop().is_some() {}
    let request = |r: &Runtime, id: &str, beta: f64| ControlRequest {
        schema_version: 1,
        calc_id: "scf".into(),
        id: id.into(),
        observed_state_version: r.snapshot().state_version,
        observed_policy_revision: r.policy_revision(),
        action: ControlAction::Policy(vec![PolicyChange::MixingBeta(beta)]),
        preconditions: vec![],
        expiry: Expiry {
            max_iteration: Some(100),
            deadline: None,
        },
        stale_policy: StalePolicy::Reject,
        reason: "test".into(),
    };
    runtime
        .submit(
            request(&runtime, "first", 0.2),
            "host".into(),
            Authority::PendingApproval,
        )
        .unwrap();
    assert_eq!(
        runtime.checkpoint(m.clone(), e.clone()).unwrap_err(),
        Rejection::PendingControls
    );
    runtime.deny("first").unwrap();
    runtime
        .submit(
            request(&runtime, "autonomous", 0.2),
            "host".into(),
            Authority::Autonomous,
        )
        .unwrap();
    runtime.safe_point().unwrap();
    assert_eq!(runtime.policy_revision(), 1);
    let cp = runtime.checkpoint(m.clone(), e.clone()).unwrap();
    assert!(
        ScfSession::restore_exact(
            fixture(
                0,
                Mixing::Pulay {
                    beta: 0.2,
                    history: 6
                },
                1,
                200
            ),
            &cp,
            &m,
            &e
        )
        .is_err()
    );
    let mut restored = Runtime::restore_exact(
        fixture(
            0,
            Mixing::Pulay {
                beta: 0.2,
                history: 6,
            },
            1,
            200,
        ),
        &cp,
        &m,
        &e,
    )
    .unwrap();
    assert_eq!(
        restored
            .checkpoint(m.clone(), e.clone())
            .unwrap()
            .runtime
            .unwrap()
            .autonomous_total,
        1
    );
    restored
        .submit(
            request(&restored, "over-quota", 0.1),
            "host".into(),
            Authority::Autonomous,
        )
        .unwrap();
    restored.safe_point().unwrap();
    assert_eq!(restored.policy_revision(), 1);
    assert!(matches!(
        restored.audit().last().unwrap().outcome,
        AuditOutcome::Rejected(ControlError::Quota)
    ));
    restored.begin_step().unwrap();
    assert_eq!(restored.checkpoint(m, e).unwrap_err(), Rejection::Boundary);
}

#[test]
fn mixer_memory_is_validated_before_restore() {
    let (m, e) = identity(0, 1);
    let mixing = Mixing::KerkerPulay {
        beta: 0.3,
        history: 6,
        q0_inverse_bohr: 1.0,
    };
    let mut s = start(0, mixing.clone(), 1, ExecutionMode::Classic);
    s.step().unwrap();
    let original = s
        .checkpoint(CheckpointKind::ExactContinuation, m.clone(), e.clone())
        .unwrap();
    let mut bad = original.clone();
    if let MixerState::Kerker { diagnostics, .. } = &mut bad.continuation.as_mut().unwrap().mixer {
        diagnostics.resets = usize::MAX - 1;
    } else {
        panic!("wrong mixer");
    }
    assert_eq!(
        ScfSession::restore_exact(fixture(0, mixing.clone(), 1, 200), &bad, &m, &e)
            .err()
            .unwrap(),
        Rejection::Mixer
    );
    let mut bad = original.clone();
    if let MixerState::Kerker { pulay, .. } = &mut bad.continuation.as_mut().unwrap().mixer {
        pulay.residuals[0].pop();
    }
    assert_eq!(
        ScfSession::restore_exact(fixture(0, mixing.clone(), 1, 200), &bad, &m, &e)
            .err()
            .unwrap(),
        Rejection::Mixer
    );
    let mut bad = original.clone();
    bad.continuation.as_mut().unwrap().iteration = usize::MAX;
    assert_eq!(
        ScfSession::restore_exact(fixture(0, mixing.clone(), 1, 200), &bad, &m, &e)
            .err()
            .unwrap(),
        Rejection::Identity
    );
    assert!(ScfSession::restore_exact(fixture(0, mixing, 1, 200), &original, &m, &e).is_ok());
}
