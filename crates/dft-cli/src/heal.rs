//! Offline demo adapter, compiled by the self_healing integration harness.
//! No numerical algorithms: all decisions pass through the engine runtime.
use dft_engine::{
    diagnostics::{Sample, Thresholds, diagnose},
    events::EventKind,
    policy::PolicyChange,
    runtime::*,
    *,
};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{
    fs,
    path::Path,
    time::{Duration, Instant, SystemTime},
};

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Decision {
    pub schema_version: u32,
    pub calc_id: String,
    pub request_id: String,
    pub observed_state_version: u64,
    pub observed_policy_revision: u64,
    pub apply_after_iteration: usize,
    pub max_iterations: usize,
}
impl Decision {
    pub fn request(&self) -> ControlRequest {
        ControlRequest {
            schema_version: self.schema_version,
            calc_id: self.calc_id.clone(),
            id: self.request_id.clone(),
            observed_state_version: self.observed_state_version,
            observed_policy_revision: self.observed_policy_revision,
            action: ControlAction::Policy(vec![PolicyChange::MaxIterations(self.max_iterations)]),
            preconditions: vec![
                Precondition::State(RuntimeState::Initialized),
                Precondition::PolicyRevision(0),
                Precondition::Mode(ExecutionMode::AgentControlled),
            ],
            expiry: Expiry {
                max_iteration: Some(self.apply_after_iteration),
                deadline: None,
            },
            stale_policy: StalePolicy::Reject,
            reason: "budget_preflight_v1: initial budget cannot certify".into(),
        }
    }
}
pub fn read_decision(path: &Path) -> Decision {
    assert!(
        fs::metadata(path).expect("decision metadata").len() <= 4096,
        "decision size limit"
    );
    let bytes = fs::read(path).expect("read recorded decision");
    assert!(bytes.len() <= 4096, "decision size limit");
    serde_json::from_slice(&bytes).expect("strict recorded decision schema")
}
pub fn hash(bytes: impl AsRef<[u8]>) -> String {
    format!("{:x}", Sha256::digest(bytes.as_ref()))
}
pub fn write_json(path: &Path, value: &Value) {
    fs::write(
        path,
        format!("{}\n", serde_json::to_string_pretty(value).unwrap()),
    )
    .unwrap();
}
fn stamp(t: SystemTime) -> String {
    t.duration_since(SystemTime::UNIX_EPOCH)
        .unwrap()
        .as_nanos()
        .to_string()
}
fn snapshot(s: &BoundarySnapshot) -> Value {
    json!({"state":format!("{:?}",s.state), "mode":format!("{:?}",s.mode), "iteration":s.iteration,
        "state_version":s.state_version,"policy_revision":s.policy_revision,"max_iterations":s.options.max_iterations,
        "expanded_options_debug":format!("{:?}",s.options)})
}
fn audit(a: &AuditRecord) -> Value {
    json!({"schema_version":a.schema_version,"calc_id":a.calc_id,"sequence":a.sequence,
        "timestamp_unix_ns":stamp(a.timestamp), "requester":a.requester,"request_debug":format!("{:?}",a.request),
        "request_id":a.request.id,"observed_state_version":a.request.observed_state_version,
        "observed_policy_revision":a.request.observed_policy_revision,"reason":a.request.reason,
        "expiry_max_iteration":a.request.expiry.max_iteration,
        "authority_debug":format!("{:?}",a.authority),"approval_debug":format!("{:?}",a.approval),
        "outcome":format!("{:?}",a.outcome),"before":snapshot(&a.before),"after":snapshot(&a.after),
        "history_reset":a.history_reset,"replayed":a.replayed,"autonomy_rule_index":a.autonomy_rule_index})
}
// Atomic batch publication before runtime mutation. This is a demo journal, not
// an exact-restart loader. A failed write propagates through JournalError.
struct Journal {
    path: std::path::PathBuf,
    bytes: Vec<u8>,
    previous: Option<String>,
}
impl AuditSink for Journal {
    fn append(&mut self, records: &[AuditRecord]) -> Result<(), JournalError> {
        use std::io::Write;
        let mut next = self.bytes.clone();
        let mut previous = self.previous.clone();
        for record in records {
            let mut value = audit(record);
            value["previous_record_sha256"] = json!(previous);
            let digest = hash(serde_json::to_vec(&value).map_err(|_| JournalError)?);
            value["record_sha256"] = json!(digest);
            previous = Some(digest);
            serde_json::to_writer(&mut next, &value).map_err(|_| JournalError)?;
            next.push(b'\n');
        }
        let temp = self.path.with_extension("tmp");
        let mut f = fs::File::create(&temp).map_err(|_| JournalError)?;
        f.write_all(&next).map_err(|_| JournalError)?;
        f.sync_all().map_err(|_| JournalError)?;
        fs::rename(temp, &self.path).map_err(|_| JournalError)?;
        self.bytes = next;
        self.previous = previous;
        Ok(())
    }
}
pub struct Run {
    pub outcome: ScfOutcome,
    pub summary: Value,
}
fn density_hash(r: &Runtime) -> String {
    hash(
        r.session()
            .density()
            .unwrap()
            .iter()
            .flat_map(|v| v.to_bits().to_le_bytes())
            .collect::<Vec<_>>(),
    )
}
fn drain(
    r: &mut Runtime,
    events: &mut Vec<Value>,
    samples: &mut Vec<Sample>,
    findings: &mut Vec<Value>,
) {
    let revision = r.policy_revision();
    let thresholds = Thresholds::from(&r.snapshot().options);
    while let Some(e) = r.events().unwrap().pop() {
        let mut scalars = json!({});
        if let EventKind::IterationFinished {
            energy_hartree,
            density_residual,
            charge_error_relative,
            mixing_applied,
        } = e.kind
        {
            scalars = json!({"energy_hartree":energy_hartree,"density_residual_bohr_neg_3_over_2":density_residual,
                "charge_error_relative":charge_error_relative,"mixing_applied":mixing_applied});
        }
        assert!(
            !matches!(e.kind, EventKind::EventsDropped { .. }),
            "demo must retain all events"
        );
        if let Some(s) = Sample::from_event(&e, revision) {
            samples.push(s);
            let f = diagnose(samples, &[], thresholds).unwrap();
            findings.push(json!({"iteration":e.iteration,"policy_revision":revision,"findings_debug":format!("{f:?}")}));
        }
        events.push(json!({"schema_version":e.schema_version,"calc_id":e.calc_id.as_ref(),"sequence":e.sequence,
            "timestamp_unix_ns":stamp(e.timestamp),"iteration":e.iteration,"state_version":e.state_version,
            "kind_debug":format!("{:?}",e.kind),"scalars":scalars}));
    }
}
pub fn run(
    problem: ScfProblem,
    dir: &Path,
    decision: Option<&Decision>,
    actor: &str,
    unsafe_probe: bool,
) -> Run {
    fs::create_dir_all(dir).unwrap();
    fs::write(dir.join("audit.jsonl"), "").unwrap();
    let start = Instant::now();
    let mut r = Runtime::new(
        ScfSession::new(
            problem,
            CompletionMode::RetainLast,
            ExecutionConfig {
                mode: ExecutionMode::AgentControlled,
                calc_id: "self-healing-demo".into(),
                ..Default::default()
            },
        )
        .unwrap(),
        RuntimeConfig {
            autonomy: Some(AutonomyPolicy {
                max_applied_total: 1,
                rules: vec![AutonomyRule {
                    action: AllowedAction::ExtendBudget { maximum: 200 },
                    max_count: 1,
                    cooldown_iterations: 0,
                    max_version_lag: 0,
                    allow_revalidation: false,
                }],
            }),
            ..Default::default()
        },
    )
    .unwrap()
    .with_audit_sink(Box::new(Journal {
        path: dir.join("audit.jsonl"),
        bytes: vec![],
        previous: None,
    }))
    .unwrap();
    r.initialize().unwrap();
    let initial = snapshot(&r.snapshot());
    let initial_density = density_hash(&r);
    let mut events = vec![];
    let mut samples = vec![];
    let mut findings = vec![];
    drain(&mut r, &mut events, &mut samples, &mut findings);
    let budget_inadequate = r.snapshot().options.max_iterations < 2;
    findings.push(json!({"rule":"budget_preflight_v1","iteration":0,"max_iterations":r.snapshot().options.max_iterations,
        "triggered":budget_inadequate,"evidence":"first iteration has no energy change; cannot certify","hypothesis":null}));
    if let Some(d) = decision {
        assert!(budget_inadequate);
        if unsafe_probe {
            let before = format!("{:?}", r.snapshot());
            let mut q = d.request();
            q.id = "unsafe-beta".into();
            q.action = ControlAction::Policy(vec![PolicyChange::MixingBeta(1.1)]);
            q.reason = "negative control: invalid beta even with operator approval".into();
            r.submit(
                q,
                "negative-control".into(),
                Authority::Approved(Approval {
                    id: "negative-control-approval".into(),
                    approver: "demo-host".into(),
                    expires_at: SystemTime::now() + Duration::from_secs(3600),
                    allow_revalidation: false,
                }),
            )
            .unwrap();
            r.safe_point().unwrap();
            assert_eq!(
                r.audit().last().unwrap().outcome,
                AuditOutcome::Rejected(ControlError::Policy(
                    dft_engine::policy::PolicyError::InvalidBeta
                ))
            );
            assert_eq!(before, format!("{:?}", r.snapshot()));
            assert_eq!(initial_density, density_hash(&r));
        }
        r.submit(d.request(), actor.into(), Authority::Autonomous)
            .unwrap();
        r.safe_point().unwrap();
        assert_eq!(r.audit().last().unwrap().outcome, AuditOutcome::Applied);
        assert_eq!(r.policy_revision(), 1);
    }
    while !r.state().terminal() {
        r.step().unwrap();
        drain(&mut r, &mut events, &mut samples, &mut findings);
    }
    let result = r.into_result();
    assert!(!result.unsafe_relaxed && !result.journal_failed);
    let outcome = result.outcome.unwrap();
    let last = outcome.result.history.last().unwrap();
    let cert = last.fixed_point_certificate.as_ref();
    let summary = json!({"state":format!("{:?}",result.state),"converged":outcome.converged,
        "energy_hartree":outcome.result.total_hartree(),"certified_energy_hartree":outcome.converged.then(||outcome.result.total_hartree()),
        "iterations":outcome.result.history.len(),"certificate_solves":outcome.result.history.iter().filter(|h|h.fixed_point_certificate.is_some()).count(),
        "operator_applications_including_certificates":outcome.result.history.iter().map(|h|h.operator_applications+h.fixed_point_certificate.as_ref().map_or(0,|c|c.operator_applications)).sum::<usize>(),
        "density_residual_bohr_neg_3_over_2":last.density_residual,"charge_error_relative":last.charge_error,
        "input_charge_error_relative":last.input_charge_error,"min_density_bohr_neg_3":last.min_density,
        "certificate_density_residual_bohr_neg_3_over_2":cert.map(|c|c.density_residual),
        "certificate_occupied_eigenvalue_change_hartree":cert.map(|c|c.occupied_eigenvalue_change),
        "policy_changes":result.policy_changes,"rejected":result.audit.iter().filter(|a|matches!(a.outcome,AuditOutcome::Rejected(_))).count(),
        "wall_seconds":start.elapsed().as_secs_f64(),"initial":initial,"initial_density_sha256":initial_density,
        "events_dropped":0,"observer_overhead":null,
        "grant_debug":format!("{:?}",result.control_config.autonomy),
        "grant_sha256":hash(format!("{:?}",result.control_config.autonomy))});
    let jsonl = events.iter().map(|v| format!("{v}\n")).collect::<String>();
    fs::write(dir.join("events.jsonl"), jsonl).unwrap();
    write_json(&dir.join("diagnostics.json"), &json!(findings));
    write_json(&dir.join("summary.json"), &summary);
    Run { outcome, summary }
}
pub fn normalized_result(outcome: &mut ScfOutcome) -> String {
    outcome.result.setup_time = Duration::ZERO;
    outcome.result.total_time = Duration::ZERO;
    for h in &mut outcome.result.history {
        h.timings = PhaseTimings::default();
        if let Some(c) = &mut h.fixed_point_certificate {
            c.timings = PhaseTimings::default();
        }
    }
    format!("{:?}", outcome.result)
}
