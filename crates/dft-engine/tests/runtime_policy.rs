use dft_engine::{Mixing, ScfOptions, policy::*};

fn floors() -> Option<ToleranceFloors> {
    // Validator fixtures only, not adopted numerical defaults.
    Some(ToleranceFloors {
        density_residual: 1e-12,
        energy_hartree: 1e-13,
        eigen_hartree: 1e-13,
    })
}

#[test]
fn rejected_batch_does_not_modify_input() {
    let original = ScfOptions::default();
    let before = format!("{original:?}").into_bytes();
    let cases = [
        (PolicyChange::MixingBeta(f64::NAN), PolicyError::InvalidBeta),
        (
            PolicyChange::MixingBeta(f64::INFINITY),
            PolicyError::InvalidBeta,
        ),
        (PolicyChange::MixingBeta(0.0), PolicyError::InvalidBeta),
        (PolicyChange::MixingBeta(1.1), PolicyError::InvalidBeta),
        (PolicyChange::MixingHistory(65), PolicyError::HistoryQuota),
        (PolicyChange::MixingHistory(0), PolicyError::HistoryQuota),
        (
            PolicyChange::DensityTolerance(1e-6),
            PolicyError::ToleranceRelaxation,
        ),
        (
            PolicyChange::EnergyToleranceHartree(1e-7),
            PolicyError::ToleranceRelaxation,
        ),
        (
            PolicyChange::EigenToleranceHartree(1e-7),
            PolicyError::ToleranceRelaxation,
        ),
        (
            PolicyChange::DensityTolerance(1e-14),
            PolicyError::BelowToleranceFloor,
        ),
        (
            PolicyChange::MaxIterations(10001),
            PolicyError::IterationQuota,
        ),
        (
            PolicyChange::ScheduleFlags {
                residual_eigen_schedule: false,
            },
            PolicyError::RestrictedSchedule,
        ),
        (
            PolicyChange::DivergenceParameters {
                factor: 30.,
                patience: 4,
            },
            PolicyError::RestrictedDivergence,
        ),
    ];
    for (change, expected) in cases {
        let result = prepare_policy(
            &original,
            3,
            floors(),
            &[PolicyChange::ResetHistory, change],
        );
        assert_eq!(result.unwrap_err(), expected);
        assert_eq!(format!("{original:?}").as_bytes(), before);
    }
}

#[test]
fn reset_semantics_and_merged_mixers() {
    let original = ScfOptions::default();
    for mixing in [
        Mixing::Linear { beta: 0.2 },
        Mixing::Pulay {
            beta: 0.2,
            history: 8,
        },
        Mixing::KerkerPulay {
            beta: 0.2,
            history: 8,
            q0_inverse_bohr: 1.0,
        },
    ] {
        let prepared =
            prepare_policy(&original, 5, floors(), &[PolicyChange::Mixing(mixing)]).unwrap();
        assert!(prepared.resets_history());
        assert_eq!(prepared.options().max_iterations, original.max_iterations);
    }
    for change in [
        PolicyChange::ResetHistory,
        PolicyChange::MixingBeta(0.2),
        PolicyChange::MixingHistory(8),
    ] {
        assert!(
            prepare_policy(&original, 5, floors(), &[change])
                .unwrap()
                .resets_history()
        );
    }
    let prepared = prepare_policy(
        &original,
        5,
        floors(),
        &[
            PolicyChange::MaxIterations(300),
            PolicyChange::DensityTolerance(1e-8),
            PolicyChange::EnergyToleranceHartree(1e-9),
            PolicyChange::EigenToleranceHartree(1e-9),
        ],
    )
    .unwrap();
    assert!(!prepared.resets_history());
    assert_eq!(prepared.options().density_tolerance, 1e-8);
    assert_eq!(prepared.options().max_iterations, 300);
    assert_eq!(original.max_iterations, 200);
}

#[test]
fn quotas_noops_and_ambiguous_batches_reject() {
    let original = ScfOptions::default();
    assert_eq!(
        prepare_policy(&original, 0, floors(), &[]).unwrap_err(),
        PolicyError::Empty
    );
    assert_eq!(
        prepare_policy(
            &original,
            0,
            floors(),
            &vec![PolicyChange::ResetHistory; 11]
        )
        .unwrap_err(),
        PolicyError::TooManyChanges
    );
    for change in [
        PolicyChange::MixingBeta(0.3),
        PolicyChange::Eigensolver(RuntimeEigensolver::Davidson),
    ] {
        assert_eq!(
            prepare_policy(&original, 0, floors(), &[change]).unwrap_err(),
            PolicyError::NoChange
        );
    }
    assert_eq!(
        prepare_policy(
            &original,
            0,
            floors(),
            &[
                PolicyChange::MixingBeta(0.2),
                PolicyChange::MixingHistory(8)
            ]
        )
        .unwrap_err(),
        PolicyError::DuplicateField
    );
    for (done, budget) in [(0, 200), (250, 250), (400, 300)] {
        assert_eq!(
            prepare_policy(
                &original,
                done,
                floors(),
                &[PolicyChange::MaxIterations(budget)]
            )
            .unwrap_err(),
            PolicyError::BudgetMustIncrease
        );
    }
    let linear = ScfOptions {
        mixing: Mixing::Linear { beta: 0.3 },
        ..original
    };
    assert_eq!(
        prepare_policy(&linear, 0, floors(), &[PolicyChange::MixingHistory(8)]).unwrap_err(),
        PolicyError::HistoryRequiresPulay
    );
}

#[test]
fn tolerance_metadata_and_certificate_ceiling_are_required() {
    let original = ScfOptions::default();
    assert_eq!(
        prepare_policy(&original, 0, None, &[PolicyChange::DensityTolerance(1e-8)]).unwrap_err(),
        PolicyError::MissingToleranceFloors
    );
    for value in [f64::NAN, f64::INFINITY, -1., 0.] {
        assert_eq!(
            prepare_policy(
                &original,
                0,
                floors(),
                &[PolicyChange::EigenToleranceHartree(value)]
            )
            .unwrap_err(),
            PolicyError::NonFiniteOrNonPositive
        );
    }
    let relaxed_launch = ScfOptions {
        density_tolerance: 1e-5,
        ..original
    };
    assert_eq!(
        prepare_policy(
            &relaxed_launch,
            0,
            floors(),
            &[PolicyChange::DensityTolerance(1e-6)]
        )
        .unwrap_err(),
        PolicyError::CertificateThreshold
    );
    let invalid = Some(ToleranceFloors {
        density_residual: 1e-3,
        ..floors().unwrap()
    });
    assert_eq!(
        prepare_policy(
            &relaxed_launch,
            0,
            invalid,
            &[PolicyChange::DensityTolerance(1e-7)]
        )
        .unwrap_err(),
        PolicyError::InvalidToleranceFloors
    );
}

#[allow(dead_code)]
mod common;
use dft_engine::{runtime::*, *};
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use std::time::{Duration, SystemTime};

fn problem(kind: usize, budget: usize) -> ScfProblem {
    let p = common::potential(2.0, 1.0, true, 401);
    let q = common::potential(2.0, 1.3, false, 401);
    let cell = Cell::from_user_rows_bohr(match kind {
        2 => [[5., 0., 0.], [0.5, 4., 0.], [0.2, 0.4, 6.]],
        _ => [[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]],
    })
    .unwrap();
    let mut atoms = vec![Atom::new("Xx", [0.13, 0.21, 0.32])];
    if kind == 2 {
        atoms.push(Atom::new("Yy", [0.6, 0.7, 0.8]));
    }
    let species = [
        SpeciesPotential {
            label: "Xx",
            potential: &p,
        },
        SpeciesPotential {
            label: "Yy",
            potential: &q,
        },
    ];
    let points: Vec<_> = (0..if kind == 1 { 4 } else { 1 })
        .map(|i| KPoint {
            fractional: [0.025 * i as f64, 0.015 * i as f64, 0.],
            weight: if kind == 1 { 0.25 } else { 1. },
        })
        .collect();
    ScfProblem::with_kpoints(
        &cell,
        &atoms,
        &species[..if kind == 2 { 2 } else { 1 }],
        ScfOptions {
            cutoff_hartree: 2.5,
            bands: Some(if kind == 2 { 4 } else { 3 }),
            max_iterations: budget,
            ..Default::default()
        },
        &points,
        OccupationPolicy::Insulator {
            minimum_gap_hartree: 1e-7,
        },
    )
    .unwrap()
}
fn runtime_config(
    mode: ExecutionMode,
    budget: usize,
    kind: usize,
    config: RuntimeConfig,
) -> Runtime {
    Runtime::new(
        ScfSession::new(
            problem(kind, budget),
            CompletionMode::RetainLast,
            ExecutionConfig {
                mode,
                calc_id: "runtime-test".into(),
                ..Default::default()
            },
        )
        .unwrap(),
        config,
    )
    .unwrap()
}
fn runtime(mode: ExecutionMode, budget: usize, kind: usize) -> Runtime {
    runtime_config(
        mode,
        budget,
        kind,
        RuntimeConfig {
            tolerance_floors: floors(),
            ..Default::default()
        },
    )
}
fn request(r: &Runtime, id: &str, action: ControlAction) -> ControlRequest {
    ControlRequest {
        schema_version: 1,
        calc_id: "runtime-test".into(),
        id: id.into(),
        observed_state_version: r.snapshot().state_version,
        observed_policy_revision: r.policy_revision(),
        action,
        preconditions: vec![],
        expiry: Expiry {
            max_iteration: Some(10000),
            deadline: None,
        },
        stale_policy: StalePolicy::Reject,
        reason: "test boundary transaction".into(),
    }
}
fn approval() -> Approval {
    Approval {
        id: "review-1".into(),
        approver: "operator".into(),
        expires_at: SystemTime::now() + Duration::from_secs(3600),
        allow_revalidation: true,
    }
}
fn control(r: &mut Runtime, id: &str, action: ControlAction) {
    let req = request(r, id, action);
    r.submit(req, "operator".into(), Authority::Approved(approval()))
        .unwrap();
    r.safe_point().unwrap();
    assert!(
        r.audit()
            .iter()
            .any(|a| a.request.id == id && a.outcome == AuditOutcome::Applied)
    );
}
fn rejection(r: &mut Runtime, req: ControlRequest, expected: ControlError) {
    let before = format!("{:?}", r.snapshot());
    let density: Vec<_> = r
        .session()
        .density()
        .unwrap()
        .iter()
        .map(|v| v.to_bits())
        .collect();
    r.submit(req, "agent".into(), Authority::Approved(approval()))
        .unwrap();
    r.safe_point().unwrap();
    assert_eq!(
        r.audit().last().unwrap().outcome,
        AuditOutcome::Rejected(expected)
    );
    assert_eq!(format!("{:?}", r.snapshot()).as_bytes(), before.as_bytes());
    assert_eq!(
        r.session()
            .density()
            .unwrap()
            .iter()
            .map(|v| v.to_bits())
            .collect::<Vec<_>>(),
        density
    );
}
fn finish(mut r: Runtime) -> RuntimeResult {
    while !r.state().terminal() {
        r.step().unwrap();
        if let Some(q) = r.events() {
            while q.pop().is_some() {}
        }
    }
    r.into_result()
}
fn normalize(mut o: ScfOutcome) -> ScfOutcome {
    o.result.setup_time = Duration::ZERO;
    o.result.total_time = Duration::ZERO;
    for h in &mut o.result.history {
        h.timings = PhaseTimings::default();
        if let Some(c) = &mut h.fixed_point_certificate {
            c.timings = PhaseTimings::default();
        }
    }
    o
}
fn identical(a: ScfOutcome, b: ScfOutcome) {
    let a = normalize(a);
    let b = normalize(b);
    assert_eq!(a.converged, b.converged);
    assert_eq!(format!("{:?}", a.result), format!("{:?}", b.result));
    assert_eq!(
        a.result.total_hartree().to_bits(),
        b.result.total_hartree().to_bits()
    );
    assert_eq!(a.result.density.shape(), b.result.density.shape());
    for (x, y) in a.result.density.iter().zip(&b.result.density) {
        assert_eq!(x.to_bits(), y.to_bits());
    }
    for (x, y) in a.result.k_points.iter().zip(&b.result.k_points) {
        assert_eq!(x.solution, y.solution);
        for j in 0..x.solution.eigenvectors.cols() {
            for (v, w) in x
                .solution
                .eigenvectors
                .column(j)
                .iter()
                .zip(y.solution.eigenvectors.column(j))
            {
                assert_eq!(v.re.to_bits(), w.re.to_bits());
                assert_eq!(v.im.to_bits(), w.im.to_bits());
            }
        }
    }
}

#[test]
fn complete_transition_table() {
    use RuntimeState::*;
    use Trigger::*;
    let states = [
        Created,
        Initialized,
        Running,
        PauseRequested,
        Paused,
        Restarting,
        Converged,
        Exhausted,
        Failed,
        Cancelled,
    ];
    let triggers = [
        Initialize,
        Start,
        Complete,
        CertificatePassed,
        BudgetExhausted,
        Fail,
        Pause,
        Resume,
        ResumeComplete,
        Cancel,
    ];
    let legal = [
        (Created, Initialize, Initialized),
        (Created, Fail, Failed),
        (Initialized, Start, Running),
        (Initialized, Pause, Paused),
        (Initialized, Cancel, Cancelled),
        (Running, Complete, Initialized),
        (Running, CertificatePassed, Converged),
        (Running, BudgetExhausted, Exhausted),
        (Running, Fail, Failed),
        (Running, Pause, PauseRequested),
        (PauseRequested, Complete, Initialized),
        (PauseRequested, CertificatePassed, Converged),
        (PauseRequested, BudgetExhausted, Exhausted),
        (PauseRequested, Fail, Failed),
        (Paused, Resume, Restarting),
        (Paused, Cancel, Cancelled),
        (Restarting, ResumeComplete, Initialized),
        (Restarting, Fail, Failed),
    ];
    for state in states {
        for trigger in triggers {
            let expected = legal.iter().find(|(s, t, _)| *s == state && *t == trigger);
            match expected {
                Some((_, _, next)) => assert_eq!(transition(state, trigger), Ok(*next)),
                None => assert_eq!(
                    transition(state, trigger),
                    Err(InvalidTransition {
                        from: state,
                        trigger
                    })
                ),
            }
        }
    }
}

#[test]
fn transactional_rejection_preserves_hidden_continuation() {
    let mut changed = runtime(ExecutionMode::AgentControlled, 5, 0);
    let mut baseline = runtime(ExecutionMode::AgentControlled, 5, 0);
    changed.initialize().unwrap();
    baseline.initialize().unwrap();
    for _ in 0..2 {
        changed.step().unwrap();
        baseline.step().unwrap();
    }
    let cases = [
        (
            ControlAction::Policy(vec![
                PolicyChange::ResetHistory,
                PolicyChange::MixingBeta(f64::NAN),
            ]),
            ControlError::Policy(PolicyError::InvalidBeta),
        ),
        (
            ControlAction::Policy(vec![
                PolicyChange::MixingBeta(0.2),
                PolicyChange::EnergyToleranceHartree(1e-6),
            ]),
            ControlError::Policy(PolicyError::ToleranceRelaxation),
        ),
        (
            ControlAction::Policy(vec![PolicyChange::ScheduleFlags {
                residual_eigen_schedule: false,
            }]),
            ControlError::Policy(PolicyError::RestrictedSchedule),
        ),
        (
            ControlAction::Checkpoint,
            ControlError::UnsupportedCapability,
        ),
        (
            ControlAction::Policy(vec![PolicyChange::Mixing(Mixing::KerkerPulay {
                beta: 0.2,
                history: 6,
                q0_inverse_bohr: 1e300,
            })]),
            ControlError::PrepareFailed,
        ),
    ];
    for (i, (action, expected)) in cases.into_iter().enumerate() {
        let req = request(&changed, &format!("reject-{i}"), action);
        rejection(&mut changed, req, expected);
    }
    let changed = finish(changed);
    let baseline = finish(baseline);
    assert_eq!(changed.policy_changes, 0);
    identical(changed.outcome.unwrap(), baseline.outcome.unwrap());
}

#[test]
fn stale_expiry_preconditions_and_approval_are_enforced() {
    let mut r = runtime(ExecutionMode::AgentControlled, 8, 0);
    r.initialize().unwrap();
    r.step().unwrap();
    let action = ControlAction::Policy(vec![PolicyChange::MixingBeta(0.2)]);
    let mut q = request(&r, "stale", action.clone());
    q.observed_state_version -= 1;
    rejection(&mut r, q, ControlError::StaleState);
    let mut q = request(&r, "future", action.clone());
    q.observed_state_version += 1;
    rejection(&mut r, q, ControlError::FutureState);
    let mut q = request(&r, "expired", action.clone());
    q.expiry.max_iteration = Some(0);
    rejection(&mut r, q, ControlError::Expired);
    let mut q = request(&r, "deadline", action.clone());
    q.expiry.deadline = Some(SystemTime::UNIX_EPOCH);
    rejection(&mut r, q, ControlError::Expired);
    let mut q = request(&r, "predicate", action.clone());
    q.preconditions.push(Precondition::DiagnosticActive);
    rejection(&mut r, q, ControlError::PreconditionFailed);
    let mut q = request(&r, "wrong-calc", action.clone());
    q.calc_id = "other".into();
    rejection(&mut r, q, ControlError::WrongCalculation);
    let mut q = request(&r, "pending", action);
    q.stale_policy = StalePolicy::Revalidate;
    r.submit(q, "agent".into(), Authority::PendingApproval)
        .unwrap();
    r.step().unwrap(); // No wait and no numerical intervention without approval.
    assert_eq!(r.snapshot().iteration, 2);
    assert_eq!(r.policy_revision(), 0);
    assert!(r.audit().iter().all(|a| a.outcome != AuditOutcome::Applied));
    r.approve("pending", approval()).unwrap();
    r.safe_point().unwrap();
    assert_eq!(r.policy_revision(), 1);
    let a = &r.audit()[r.audit().len() - 2..];
    assert_eq!(a[0].outcome, AuditOutcome::Accepted);
    assert_eq!(a[1].outcome, AuditOutcome::Applied);
    assert_eq!(a[1].request.observed_state_version, 2);
    assert_eq!(a[1].before.state_version, 3);
    assert_eq!(a[1].after.state_version, 4);
    assert!(a[1].history_reset);
    let mut q = request(&r, "stale-policy", ControlAction::Pause);
    q.observed_policy_revision = 0;
    q.stale_policy = StalePolicy::Revalidate;
    rejection(&mut r, q, ControlError::StalePolicy);
}

#[test]
fn pause_cancel_and_terminal_precedence() {
    let mut r = runtime(ExecutionMode::AgentControlled, 4, 0);
    r.initialize().unwrap();
    r.begin_step().unwrap();
    let mut q = request(&r, "pause", ControlAction::Pause);
    q.stale_policy = StalePolicy::Revalidate;
    r.submit(q, "operator".into(), Authority::Approved(approval()))
        .unwrap();
    assert_eq!(r.state(), RuntimeState::PauseRequested);
    assert!(r.safe_point().is_err());
    assert_eq!(r.snapshot().iteration, 0);
    r.finish_step().unwrap();
    assert_eq!(r.state(), RuntimeState::Paused);
    assert_eq!(r.snapshot().iteration, 1);
    assert!(r.step().is_err());
    control(&mut r, "resume", ControlAction::Resume);
    assert_eq!(r.state(), RuntimeState::Initialized);
    r.begin_step().unwrap();
    for (id, action) in [
        ("pause-2", ControlAction::Pause),
        ("cancel", ControlAction::Cancel),
    ] {
        let mut q = request(&r, id, action);
        q.stale_policy = StalePolicy::Revalidate;
        r.submit(q, "operator".into(), Authority::Approved(approval()))
            .unwrap();
    }
    r.finish_step().unwrap();
    assert_eq!(r.state(), RuntimeState::Cancelled);
    assert_eq!(r.snapshot().iteration, 2);
    assert!(r.session().result().is_none());
    assert_eq!(
        r.audit().last().unwrap().outcome,
        AuditOutcome::Rejected(ControlError::Terminal)
    );
    assert!(r.step().is_err());
    let mut terminal = runtime(ExecutionMode::AgentControlled, 1, 0);
    terminal.initialize().unwrap();
    terminal.begin_step().unwrap();
    let mut q = request(&terminal, "late-cancel", ControlAction::Cancel);
    q.stale_policy = StalePolicy::Revalidate;
    terminal
        .submit(q, "operator".into(), Authority::Approved(approval()))
        .unwrap();
    terminal.finish_step().unwrap();
    assert_eq!(terminal.state(), RuntimeState::Exhausted);
    assert_eq!(
        terminal.audit()[0].outcome,
        AuditOutcome::Rejected(ControlError::Terminal)
    );
    assert!(!terminal.into_result().outcome.unwrap().converged);
}

#[test]
fn mode_switches_preserve_numerics_and_revoke_pending_controls() {
    let mut r = runtime(ExecutionMode::Classic, 3, 0);
    control(
        &mut r,
        "observe-before-start",
        ControlAction::SwitchMode(ExecutionMode::Observed),
    );
    r.initialize().unwrap();
    let req = request(
        &r,
        "readonly",
        ControlAction::Policy(vec![PolicyChange::MixingBeta(0.2)]),
    );
    rejection(&mut r, req, ControlError::ModeDisallowsControl);
    control(
        &mut r,
        "enable",
        ControlAction::SwitchMode(ExecutionMode::AgentControlled),
    );
    let q = request(
        &r,
        "waiting",
        ControlAction::Policy(vec![PolicyChange::MixingBeta(0.2)]),
    );
    r.submit(q, "agent".into(), Authority::PendingApproval)
        .unwrap();
    control(
        &mut r,
        "disable",
        ControlAction::SwitchMode(ExecutionMode::Observed),
    );
    // Last record can be pending proposal rejection after the applied mode change.
    assert_eq!(r.policy_revision(), 0);
    let req = request(
        &r,
        "classic-forbidden",
        ControlAction::SwitchMode(ExecutionMode::Classic),
    );
    rejection(&mut r, req, ControlError::ModeUnavailableAfterStart);
    let mut baseline = runtime(ExecutionMode::Classic, 3, 0);
    baseline.initialize().unwrap();
    identical(
        finish(r).outcome.unwrap(),
        finish(baseline).outcome.unwrap(),
    );
    let mut classic = runtime(ExecutionMode::Classic, 1, 0);
    classic.initialize().unwrap();
    let req = request(
        &classic,
        "too-late",
        ControlAction::SwitchMode(ExecutionMode::Observed),
    );
    rejection(&mut classic, req, ControlError::ModeUnavailableAfterStart);
}

#[test]
fn autonomous_allowlist_counts_ranges_cooldown_and_revocation() {
    let config = RuntimeConfig {
        tolerance_floors: floors(),
        autonomy: Some(AutonomyPolicy {
            max_applied_total: 2,
            rules: vec![AutonomyRule {
                action: AllowedAction::Beta {
                    minimum: 0.1,
                    maximum: 0.3,
                    only_decrease: true,
                },
                max_count: 2,
                cooldown_iterations: 2,
                max_version_lag: 2,
                allow_revalidation: true,
            }],
        }),
        ..Default::default()
    };
    let mut r = runtime_config(ExecutionMode::AgentControlled, 8, 0, config);
    r.initialize().unwrap();
    for (id, beta, expected) in [
        (
            "range",
            0.8,
            AuditOutcome::Rejected(ControlError::AutonomousDenied),
        ),
        ("first", 0.25, AuditOutcome::Applied),
        (
            "cooldown",
            0.2,
            AuditOutcome::Rejected(ControlError::Cooldown),
        ),
    ] {
        let q = request(
            &r,
            id,
            ControlAction::Policy(vec![PolicyChange::MixingBeta(beta)]),
        );
        r.submit(q, "agent".into(), Authority::Autonomous).unwrap();
        r.safe_point().unwrap();
        assert_eq!(r.audit().last().unwrap().outcome, expected);
    }
    r.step().unwrap();
    r.step().unwrap();
    let q = request(
        &r,
        "second",
        ControlAction::Policy(vec![PolicyChange::MixingBeta(0.2)]),
    );
    r.submit(q, "agent".into(), Authority::Autonomous).unwrap();
    r.safe_point().unwrap();
    assert_eq!(r.policy_revision(), 2);
    let q = request(&r, "third", ControlAction::Pause);
    r.submit(q, "agent".into(), Authority::Autonomous).unwrap();
    r.safe_point().unwrap();
    assert_eq!(
        r.audit().last().unwrap().outcome,
        AuditOutcome::Rejected(ControlError::Quota)
    );
    control(
        &mut r,
        "observe",
        ControlAction::SwitchMode(ExecutionMode::Observed),
    );
    control(
        &mut r,
        "control",
        ControlAction::SwitchMode(ExecutionMode::AgentControlled),
    );
    let q = request(&r, "revoked", ControlAction::Pause);
    r.submit(q, "agent".into(), Authority::Autonomous).unwrap();
    r.safe_point().unwrap();
    assert_eq!(
        r.audit().last().unwrap().outcome,
        AuditOutcome::Rejected(ControlError::AutonomousDenied)
    );
}

#[test]
fn idempotency_quota_and_failed_journal_do_not_commit() {
    let mut r = runtime(ExecutionMode::AgentControlled, 3, 0);
    r.initialize().unwrap();
    let q = request(
        &r,
        "once",
        ControlAction::Policy(vec![PolicyChange::MixingBeta(0.2)]),
    );
    r.submit(
        q.clone(),
        "operator".into(),
        Authority::Approved(approval()),
    )
    .unwrap();
    r.safe_point().unwrap();
    let version = r.snapshot().state_version;
    r.submit(
        q.clone(),
        "operator".into(),
        Authority::Approved(approval()),
    )
    .unwrap();
    r.safe_point().unwrap();
    assert_eq!(r.audit().len(), 2);
    assert_eq!(r.snapshot().state_version, version);
    let mut conflict = q;
    conflict.action = ControlAction::Cancel;
    assert_eq!(
        r.submit(conflict, "operator".into(), Authority::Approved(approval())),
        Err(ControlError::Conflict)
    );
    struct FailedSink;
    impl AuditSink for FailedSink {
        fn append(&mut self, _: &[AuditRecord]) -> Result<(), JournalError> {
            Err(JournalError)
        }
    }
    let mut failed = runtime(ExecutionMode::AgentControlled, 3, 0)
        .with_audit_sink(Box::new(FailedSink))
        .unwrap();
    failed.initialize().unwrap();
    let before = format!("{:?}", failed.snapshot());
    let q = request(
        &failed,
        "write-failure",
        ControlAction::Policy(vec![PolicyChange::MixingBeta(0.2)]),
    );
    failed
        .submit(q, "operator".into(), Authority::Approved(approval()))
        .unwrap();
    assert_eq!(failed.safe_point(), Err(ControlError::JournalFailed));
    assert_eq!(format!("{:?}", failed.snapshot()), before);
    assert_eq!(failed.policy_revision(), 0);
    assert!(matches!(
        failed.step(),
        Err(RuntimeError::Control(ControlError::JournalFailed))
    ));
    let mut bounded = runtime_config(
        ExecutionMode::AgentControlled,
        3,
        0,
        RuntimeConfig {
            max_audit_records: 2,
            ..Default::default()
        },
    );
    bounded.initialize().unwrap();
    control(&mut bounded, "one", ControlAction::Pause);
    let q = request(&bounded, "two", ControlAction::Resume);
    assert_eq!(
        bounded.submit(q, "operator".into(), Authority::Approved(approval())),
        Err(ControlError::Capacity)
    );
    assert_eq!(bounded.audit().len(), 2);
}

#[test]
#[ignore = "SCF replay matrix; cargo test --release -p dft-engine --test runtime_policy -- --include-ignored --nocapture"]
fn recorded_interventions_replay_bit_identically_and_certify() {
    for kind in 0..3 {
        let mut r = runtime(ExecutionMode::AgentControlled, 200, kind);
        r.initialize().unwrap();
        r.step().unwrap();
        control(
            &mut r,
            "at-1",
            ControlAction::Policy(vec![PolicyChange::Mixing(Mixing::Linear { beta: 0.2 })]),
        );
        r.step().unwrap();
        control(
            &mut r,
            "at-2",
            ControlAction::Policy(vec![PolicyChange::Mixing(Mixing::KerkerPulay {
                beta: 0.3,
                history: 6,
                q0_inverse_bohr: 1.0,
            })]),
        );
        control(&mut r, "pause", ControlAction::Pause);
        control(&mut r, "resume", ControlAction::Resume);
        r.step().unwrap();
        control(
            &mut r,
            "at-3",
            ControlAction::Policy(vec![
                PolicyChange::ResetHistory,
                PolicyChange::MaxIterations(250),
            ]),
        );
        let log = r.audit().to_vec();
        let first = finish(r);
        let mut replay = runtime(ExecutionMode::AgentControlled, 200, kind);
        replay.replay(&log).unwrap();
        let second = finish(replay);
        assert_eq!(first.policy_changes, 3);
        assert_eq!(second.policy_changes, 3);
        assert!(!first.unsafe_relaxed);
        let a = first.outcome.unwrap();
        let b = second.outcome.unwrap();
        assert!(a.converged && b.converged);
        let certificate = a
            .result
            .history
            .last()
            .unwrap()
            .fixed_point_certificate
            .as_ref()
            .unwrap();
        assert!(certificate.density_residual <= 1e-7);
        assert!(certificate.occupied_eigenvalue_change <= 1e-7);
        println!(
            "runtime replay fixture={kind}: energy={:.17} Ha; iterations={}; certificate density={:.9e} Bohr^-3/2; eigenvalue change={:.9e} Ha; policy_changes=3; full bits identical",
            a.result.total_hartree(),
            a.result.history.len(),
            certificate.density_residual,
            certificate.occupied_eigenvalue_change
        );
        identical(a, b);
        let mut corrupted = log.clone();
        corrupted[1].after.options.density_tolerance *= 10.;
        let mut rejected = runtime(ExecutionMode::AgentControlled, 200, kind);
        assert!(matches!(
            rejected.replay(&corrupted),
            Err(RuntimeError::Control(ControlError::ReplayMismatch))
        ));
        assert_eq!(rejected.policy_revision(), 0);
    }
}

#[test]
fn boundary_schema_approval_and_audit_contracts() {
    let mut r = runtime(ExecutionMode::AgentControlled, 6, 0);
    r.initialize().unwrap();
    let mut q = request(&r, "oversized", ControlAction::Pause);
    q.reason = "x".repeat(513);
    assert_eq!(
        r.submit(q, "agent".into(), Authority::PendingApproval),
        Err(ControlError::InvalidSchema)
    );
    assert!(r.audit().is_empty());
    let q = request(&r, "approval-expired", ControlAction::Pause);
    let mut a = approval();
    a.expires_at = SystemTime::UNIX_EPOCH;
    r.submit(q, "agent".into(), Authority::Approved(a)).unwrap();
    r.safe_point().unwrap();
    assert_eq!(
        r.audit().last().unwrap().outcome,
        AuditOutcome::Rejected(ControlError::ApprovalExpired)
    );
    let q = request(&r, "denied", ControlAction::Pause);
    r.submit(q, "agent".into(), Authority::PendingApproval)
        .unwrap();
    r.deny("denied").unwrap();
    assert_eq!(
        r.audit().last().unwrap().outcome,
        AuditOutcome::Rejected(ControlError::ApprovalDenied)
    );
    r.begin_step().unwrap();
    let mut q = request(&r, "unapproved-pause", ControlAction::Pause);
    q.stale_policy = StalePolicy::Revalidate;
    r.submit(q, "agent".into(), Authority::PendingApproval)
        .unwrap();
    r.finish_step().unwrap();
    assert_eq!(r.state(), RuntimeState::Initialized);
    assert_eq!(r.snapshot().iteration, 1);
    r.deny("unapproved-pause").unwrap();
    let mut q = request(&r, "no-stale-approval", ControlAction::Pause);
    q.observed_state_version -= 1;
    q.stale_policy = StalePolicy::Revalidate;
    let mut a = approval();
    a.allow_revalidation = false;
    r.submit(q, "agent".into(), Authority::Approved(a)).unwrap();
    r.safe_point().unwrap();
    assert_eq!(
        r.audit().last().unwrap().outcome,
        AuditOutcome::Rejected(ControlError::StaleState)
    );
    control(
        &mut r,
        "reset",
        ControlAction::Policy(vec![PolicyChange::ResetHistory]),
    );
    let q = request(
        &r,
        "second-same-boundary",
        ControlAction::Policy(vec![PolicyChange::MixingBeta(0.2)]),
    );
    rejection(&mut r, q, ControlError::Superseded);
    for (i, a) in r.audit().iter().enumerate() {
        assert_eq!(a.sequence, i as u64 + 1);
        assert_eq!(a.schema_version, 1);
        assert_eq!(a.calc_id, "runtime-test");
        assert!(!a.requester.is_empty());
        assert!(a.request.reason.len() <= 512);
        match a.outcome {
            AuditOutcome::Accepted | AuditOutcome::Rejected(_) => {
                assert_eq!(format!("{:?}", a.before), format!("{:?}", a.after))
            }
            AuditOutcome::Applied => assert_eq!(a.after.state_version, a.before.state_version + 1),
        }
    }
}

#[test]
fn solver_tolerance_and_mode_log_replay_without_clock_dependency() {
    let mut r = runtime(ExecutionMode::Observed, 5, 0);
    r.initialize().unwrap();
    control(
        &mut r,
        "enable",
        ControlAction::SwitchMode(ExecutionMode::AgentControlled),
    );
    r.step().unwrap();
    control(
        &mut r,
        "solver",
        ControlAction::Policy(vec![PolicyChange::Eigensolver(RuntimeEigensolver::Hybrid)]),
    );
    r.step().unwrap();
    control(
        &mut r,
        "tighten",
        ControlAction::Policy(vec![
            PolicyChange::DensityTolerance(5e-8),
            PolicyChange::EnergyToleranceHartree(5e-9),
            PolicyChange::EigenToleranceHartree(5e-9),
        ]),
    );
    let changed = r.audit().last().unwrap();
    assert!(!changed.history_reset);
    control(&mut r, "pause", ControlAction::Pause);
    control(&mut r, "resume", ControlAction::Resume);
    control(
        &mut r,
        "observe",
        ControlAction::SwitchMode(ExecutionMode::Observed),
    );
    let mut log = r.audit().to_vec();
    for a in &mut log {
        a.timestamp = SystemTime::UNIX_EPOCH;
        a.request.expiry.deadline = Some(SystemTime::UNIX_EPOCH);
    }
    let expected = finish(r);
    let mut replay = runtime(ExecutionMode::Observed, 5, 0);
    replay.replay(&log).unwrap();
    let actual = finish(replay);
    assert_eq!(actual.policy_changes, 2);
    identical(expected.outcome.unwrap(), actual.outcome.unwrap());
}

#[test]
fn silent_supervised_agent_and_classic_wrapper_match_the_existing_session() {
    let reference = problem(0, 4).run_outcome().unwrap();
    for mode in [
        ExecutionMode::Classic,
        ExecutionMode::Observed,
        ExecutionMode::AgentControlled,
    ] {
        let mut r = runtime(mode, 4, 0);
        r.initialize().unwrap();
        if mode == ExecutionMode::AgentControlled {
            let mut q = request(
                &r,
                "silent",
                ControlAction::Policy(vec![PolicyChange::MixingBeta(0.1)]),
            );
            q.stale_policy = StalePolicy::Revalidate;
            r.submit(q, "agent".into(), Authority::PendingApproval)
                .unwrap();
        }
        let output = finish(r);
        assert_eq!(output.policy_changes, 0);
        if mode == ExecutionMode::Classic {
            assert!(output.audit.is_empty());
        }
        if mode == ExecutionMode::AgentControlled {
            assert_eq!(
                output.audit.last().unwrap().outcome,
                AuditOutcome::Rejected(ControlError::Terminal)
            );
        }
        identical(
            problem(0, 4).run_outcome().unwrap(),
            output.outcome.unwrap(),
        );
    }
    assert!(!reference.converged);
}

#[test]
fn certificate_completion_wins_over_pending_cancel() {
    let baseline = problem(0, 200).run_outcome().unwrap();
    assert!(baseline.converged);
    let iterations = baseline.result.history.len();
    let mut r = runtime(ExecutionMode::AgentControlled, 200, 0);
    r.initialize().unwrap();
    for _ in 1..iterations {
        r.step().unwrap();
    }
    r.begin_step().unwrap();
    let mut q = request(&r, "cancel-at-certificate", ControlAction::Cancel);
    q.stale_policy = StalePolicy::Revalidate;
    r.submit(q, "operator".into(), Authority::Approved(approval()))
        .unwrap();
    r.finish_step().unwrap();
    assert_eq!(r.state(), RuntimeState::Converged);
    assert_eq!(
        r.audit().last().unwrap().outcome,
        AuditOutcome::Rejected(ControlError::Terminal)
    );
    identical(baseline, r.into_result().outcome.unwrap());
}

#[test]
fn same_mode_is_an_audited_idempotent_response() {
    let mut r = runtime(ExecutionMode::Observed, 3, 0);
    r.initialize().unwrap();
    let before = format!("{:?}", r.snapshot());
    control(
        &mut r,
        "same-mode",
        ControlAction::SwitchMode(ExecutionMode::Observed),
    );
    assert_eq!(format!("{:?}", r.snapshot()), before);
    assert_eq!(r.audit().len(), 2);
}

#[test]
fn replay_preserves_autonomous_quotas_and_empty_log_runs() {
    let config = RuntimeConfig {
        autonomy: Some(AutonomyPolicy {
            max_applied_total: 1,
            rules: vec![AutonomyRule {
                action: AllowedAction::ResetHistory,
                max_count: 1,
                cooldown_iterations: 0,
                max_version_lag: 0,
                allow_revalidation: false,
            }],
        }),
        ..Default::default()
    };
    let mut r = runtime_config(ExecutionMode::AgentControlled, 3, 0, config.clone());
    r.initialize().unwrap();
    let q = request(
        &r,
        "reset-autonomously",
        ControlAction::Policy(vec![PolicyChange::ResetHistory]),
    );
    r.submit(q, "agent".into(), Authority::Autonomous).unwrap();
    r.safe_point().unwrap();
    let mut replay = runtime_config(ExecutionMode::AgentControlled, 3, 0, config);
    replay.replay(r.audit()).unwrap();
    replay.step().unwrap();
    let q = request(
        &replay,
        "cannot-reset-again",
        ControlAction::Policy(vec![PolicyChange::ResetHistory]),
    );
    replay
        .submit(q, "agent".into(), Authority::Autonomous)
        .unwrap();
    replay.safe_point().unwrap();
    assert_eq!(
        replay.audit().last().unwrap().outcome,
        AuditOutcome::Rejected(ControlError::Quota)
    );
    assert!(replay.audit()[1].replayed);
    let mut empty = runtime(ExecutionMode::Classic, 3, 0);
    empty.replay(&[]).unwrap();
    identical(
        problem(0, 3).run_outcome().unwrap(),
        finish(empty).outcome.unwrap(),
    );
}
