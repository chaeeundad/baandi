//! Single-owner control around the unchanged SCF step. Authentication, transport,
//! policy-file loading and durable journal storage belong to trusted adapters.
//! No controller is installed in the ordinary `ScfSession`/Classic driver.
use crate::policy::{PolicyChange, PolicyError, RuntimeEigensolver, ToleranceFloors};
use crate::session::{BoundaryMode, BoundaryPolicy};
use crate::{ExecutionMode, ScfError, ScfOptions, ScfOutcome, ScfSession, StepReport, StepStatus};
use std::time::{Instant, SystemTime};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RuntimeState {
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
}
impl RuntimeState {
    pub fn terminal(self) -> bool {
        matches!(
            self,
            Self::Converged | Self::Exhausted | Self::Failed | Self::Cancelled
        )
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Trigger {
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
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct InvalidTransition {
    pub from: RuntimeState,
    pub trigger: Trigger,
}

/// Pure transition table. It cannot mutate a runtime or confer certification.
pub fn transition(from: RuntimeState, trigger: Trigger) -> Result<RuntimeState, InvalidTransition> {
    use RuntimeState::*;
    use Trigger::*;
    let next = match (from, trigger) {
        (Created, Initialize) => Initialized,
        (Created | Running | PauseRequested | Restarting, Fail) => Failed,
        (Initialized, Start) => Running,
        (Running, Complete) => Initialized,
        (PauseRequested, Complete) => Initialized,
        (Running | PauseRequested, CertificatePassed) => Converged,
        (Running | PauseRequested, BudgetExhausted) => Exhausted,
        (Initialized, Pause) => Paused,
        (Running, Pause) => PauseRequested,
        (Paused, Resume) => Restarting,
        (Restarting, ResumeComplete) => Initialized,
        (Initialized | Paused, Cancel) => Cancelled,
        _ => return Err(InvalidTransition { from, trigger }),
    };
    Ok(next)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ControlError {
    InvalidSchema,
    WrongCalculation,
    StaleState,
    FutureState,
    StalePolicy,
    Expired,
    PreconditionFailed,
    ApprovalRequired,
    ApprovalExpired,
    ApprovalDenied,
    AutonomousDenied,
    Cooldown,
    Quota,
    Capacity,
    Conflict,
    ModeDisallowsControl,
    ModeUnavailableAfterStart,
    Terminal,
    NoChange,
    UnsupportedCapability,
    PrepareFailed,
    JournalFailed,
    ReplayMismatch,
    Superseded,
    Policy(PolicyError),
    InvalidTransition(InvalidTransition),
}
impl From<InvalidTransition> for ControlError {
    fn from(e: InvalidTransition) -> Self {
        Self::InvalidTransition(e)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum ControlAction {
    Policy(Vec<PolicyChange>),
    Pause,
    Resume,
    Cancel,
    /// Administrative only, never on the autonomous allow-list.
    SwitchMode(ExecutionMode),
    /// Reserved for P11C. Explicitly rejected, never advertised as implemented.
    Checkpoint,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StalePolicy {
    Reject,
    Revalidate,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Scalar {
    DensityResidual,
    ChargeErrorRelative,
    EnergyChangeAbsHartree,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Comparison {
    Less,
    LessOrEqual,
    Greater,
    GreaterOrEqual,
}
#[derive(Clone, Debug, PartialEq)]
pub enum Precondition {
    State(RuntimeState),
    PolicyRevision(u64),
    Mode(ExecutionMode),
    Scalar {
        field: Scalar,
        comparison: Comparison,
        value: f64,
    },
    /// Diagnostics are P11D; unavailable evidence always fails closed.
    DiagnosticActive,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Expiry {
    pub max_iteration: Option<usize>,
    pub deadline: Option<SystemTime>,
}
#[derive(Clone, Debug, PartialEq)]
pub struct ControlRequest {
    pub schema_version: u32,
    pub calc_id: String,
    pub id: String,
    pub observed_state_version: u64,
    pub observed_policy_revision: u64,
    pub action: ControlAction,
    pub preconditions: Vec<Precondition>,
    pub expiry: Expiry,
    pub stale_policy: StalePolicy,
    /// Display-only text, at most 512 UTF-8 bytes. Never interpreted as authority.
    pub reason: String,
}

/// Supplied by the authenticated driver, never deserialized from agent parameters.
/// Approval is bound to the runtime's owned immutable request, not a reusable token.
#[derive(Clone, Debug)]
pub struct Approval {
    pub id: String,
    pub approver: String,
    pub expires_at: SystemTime,
    pub allow_revalidation: bool,
}
#[derive(Clone, Debug)]
pub enum Authority {
    PendingApproval,
    Approved(Approval),
    Autonomous,
}

#[derive(Clone, Debug)]
pub enum AllowedAction {
    Beta {
        minimum: f64,
        maximum: f64,
        only_decrease: bool,
    },
    ResetHistory,
    Eigensolver(RuntimeEigensolver),
    TightenTolerances,
    ExtendBudget {
        maximum: usize,
    },
    Pause,
    Resume,
    Cancel,
}
#[derive(Clone, Debug)]
pub struct AutonomyRule {
    pub action: AllowedAction,
    pub max_count: usize,
    pub cooldown_iterations: usize,
    pub max_version_lag: u64,
    pub allow_revalidation: bool,
}
#[derive(Clone, Debug)]
pub struct AutonomyPolicy {
    pub max_applied_total: usize,
    pub rules: Vec<AutonomyRule>,
}
#[derive(Clone, Debug)]
pub struct RuntimeConfig {
    pub tolerance_floors: Option<ToleranceFloors>,
    /// None is Supervised. Installation is a trusted launch-time operation.
    pub autonomy: Option<AutonomyPolicy>,
    pub max_audit_records: usize,
    pub max_pending: usize,
}
impl Default for RuntimeConfig {
    fn default() -> Self {
        Self {
            tolerance_floors: None,
            autonomy: None,
            max_audit_records: 4096,
            max_pending: 16,
        }
    }
}

#[derive(Clone, Debug)]
pub struct BoundarySnapshot {
    pub state: RuntimeState,
    pub mode: ExecutionMode,
    pub iteration: usize,
    pub state_version: u64,
    pub policy_revision: u64,
    pub options: ScfOptions,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AuditOutcome {
    Accepted,
    Applied,
    Rejected(ControlError),
}
#[derive(Clone, Debug)]
pub struct AuditRecord {
    pub schema_version: u32,
    pub calc_id: String,
    pub sequence: u64,
    pub timestamp: SystemTime,
    pub request: ControlRequest,
    pub requester: String,
    pub approval: Option<Approval>,
    pub authority: Authority,
    pub replayed: bool,
    pub autonomy_rule_index: Option<usize>,
    pub before: BoundarySnapshot,
    pub after: BoundarySnapshot,
    pub outcome: AuditOutcome,
    pub history_reset: bool,
}

/// Optional durable adapter. `append` must atomically append the entire batch,
/// or return Err without publishing any of it. Called only at a control boundary.
/// A failed sink stops further driving; no numerical/policy mutation is committed.
/// The in-memory default makes no crash-recovery or durable-storage promise.
pub trait AuditSink {
    fn append(&mut self, records: &[AuditRecord]) -> Result<(), JournalError>;
}
#[derive(Clone, Copy, Debug)]
pub struct JournalError;

struct Pending {
    request: ControlRequest,
    requester: String,
    authority: Authority,
    deadline: Option<Instant>,
    approval_deadline: Option<Instant>,
    replayed: bool,
}
enum Mutation {
    Noop,
    Policy(BoundaryPolicy),
    Mode(BoundaryMode),
    State(RuntimeState),
}
/// Complete controlled result, including policy-change provenance even when no
/// certified numerical outcome exists (cancelled, failed or strict exhaustion).
pub struct RuntimeResult {
    pub outcome: Option<ScfOutcome>,
    pub state: RuntimeState,
    pub policy_changes: u64,
    pub unsafe_relaxed: bool,
    pub journal_failed: bool,
    pub control_config: RuntimeConfig,
    pub audit: Vec<AuditRecord>,
}
pub struct Runtime {
    session: ScfSession,
    state: RuntimeState,
    revision: u64,
    config: RuntimeConfig,
    pending: Vec<Pending>,
    audit: Vec<AuditRecord>,
    sink: Option<Box<dyn AuditSink>>,
    journal_failed: bool,
    autonomous_total: usize,
    autonomy_enabled: bool,
    last_policy_iteration: Option<usize>,
    rule_counts: Vec<usize>,
    rule_last_iteration: Vec<Option<usize>>,
}
impl Runtime {
    pub fn new(session: ScfSession, config: RuntimeConfig) -> Result<Self, ControlError> {
        if session.snapshot().status != StepStatus::Created
            || config.max_pending == 0
            || config.max_pending > 16
            || config.max_audit_records < 2
            || config.max_audit_records > 65536
        {
            return Err(ControlError::InvalidSchema);
        }
        if let Some(f) = config.tolerance_floors {
            let o = session.snapshot().options;
            for (floor, initial) in [
                (f.density_residual, o.density_tolerance),
                (f.energy_hartree, o.energy_tolerance),
                (f.eigen_hartree, o.eigen_tolerance),
            ] {
                if !floor.is_finite() || floor <= 0.0 || floor > initial {
                    return Err(ControlError::Policy(PolicyError::InvalidToleranceFloors));
                }
            }
        }
        let rules = config.autonomy.as_ref().map_or(0, |p| p.rules.len());
        if let Some(p) = &config.autonomy {
            if p.max_applied_total == 0 || p.max_applied_total > 32 || rules > 16 {
                return Err(ControlError::InvalidSchema);
            }
            for r in &p.rules {
                if r.max_count == 0 || r.max_count > 32 {
                    return Err(ControlError::InvalidSchema);
                }
                match r.action {
                    AllowedAction::Beta {
                        minimum, maximum, ..
                    } if !minimum.is_finite()
                        || !maximum.is_finite()
                        || minimum <= 0.0
                        || minimum > maximum
                        || maximum > 1.0 =>
                    {
                        return Err(ControlError::InvalidSchema);
                    }
                    AllowedAction::ExtendBudget { maximum } if maximum == 0 || maximum > 10000 => {
                        return Err(ControlError::InvalidSchema);
                    }
                    _ => {}
                }
            }
        }
        Ok(Self {
            session,
            state: RuntimeState::Created,
            revision: 0,
            pending: Vec::new(),
            audit: Vec::new(),
            sink: None,
            journal_failed: false,
            autonomous_total: 0,
            autonomy_enabled: config.autonomy.is_some(),
            last_policy_iteration: None,
            rule_counts: vec![0; rules],
            rule_last_iteration: vec![None; rules],
            config,
        })
    }
    /// Install the adapter before any initialization or control. No agent may call this.
    pub fn with_audit_sink(mut self, sink: Box<dyn AuditSink>) -> Result<Self, ControlError> {
        if self.state != RuntimeState::Created || !self.audit.is_empty() || !self.pending.is_empty()
        {
            return Err(ControlError::InvalidSchema);
        }
        self.sink = Some(sink);
        Ok(self)
    }
    pub fn state(&self) -> RuntimeState {
        self.state
    }
    pub fn policy_revision(&self) -> u64 {
        self.revision
    }
    pub fn session(&self) -> &ScfSession {
        &self.session
    }
    pub fn events(&mut self) -> Option<&mut crate::events::EventQueue> {
        self.session.events()
    }
    pub fn audit(&self) -> &[AuditRecord] {
        &self.audit
    }
    pub fn snapshot(&self) -> BoundarySnapshot {
        let s = self.session.snapshot();
        BoundarySnapshot {
            state: self.state,
            mode: s.mode,
            iteration: s.iteration,
            state_version: s.state_version,
            policy_revision: self.revision,
            options: s.options.clone(),
        }
    }
    pub fn into_result(self) -> RuntimeResult {
        RuntimeResult {
            outcome: self.session.into_result(),
            state: self.state,
            policy_changes: self.revision,
            unsafe_relaxed: false,
            journal_failed: self.journal_failed,
            control_config: self.config,
            audit: self.audit,
        }
    }
    pub fn initialize(&mut self) -> Result<(), RuntimeError> {
        let next = transition(self.state, Trigger::Initialize).map_err(ControlError::from)?;
        self.healthy()?;
        match self.session.initialize() {
            Ok(()) => {
                self.state = next;
                Ok(())
            }
            Err(e) => {
                self.state = RuntimeState::Failed;
                Err(RuntimeError::Scf(e))
            }
        }
    }
    /// Splitting start/finish lets a driver queue controls while an iteration is
    /// scheduled. The session is exclusively borrowed during the actual solve;
    /// asynchronous transports deliver their messages at the next boundary.
    pub fn begin_step(&mut self) -> Result<(), ControlError> {
        self.safe_point()?;
        self.state = transition(self.state, Trigger::Start)?;
        Ok(())
    }
    pub fn finish_step(&mut self) -> Result<StepReport, RuntimeError> {
        self.healthy()?;
        if !matches!(
            self.state,
            RuntimeState::Running | RuntimeState::PauseRequested
        ) {
            return Err(ControlError::InvalidTransition(InvalidTransition {
                from: self.state,
                trigger: Trigger::Complete,
            })
            .into());
        }
        let result = self.session.step();
        let trigger = match &result {
            Ok(r) => match r.status {
                StepStatus::Ready => Trigger::Complete,
                StepStatus::Converged => Trigger::CertificatePassed,
                StepStatus::NotConverged => Trigger::BudgetExhausted,
                _ => unreachable!("step result is terminal or ready"),
            },
            Err(_) => Trigger::Fail,
        };
        // PauseRequested describes a pending proposal, not approval. Only a
        // successful boundary transaction can actually enter Paused.
        self.state = transition(self.state, trigger).map_err(ControlError::from)?;
        self.safe_point()?;
        result.map_err(RuntimeError::Scf)
    }
    /// With no incoming requests this is deterministic continue_current; it
    /// never waits for an agent, approval or subscriber. Classic allocates/logs
    /// no control data on this path.
    pub fn step(&mut self) -> Result<StepReport, RuntimeError> {
        self.begin_step()?;
        self.finish_step()
    }
    fn healthy(&self) -> Result<(), ControlError> {
        if self.journal_failed {
            Err(ControlError::JournalFailed)
        } else {
            Ok(())
        }
    }
    fn schema(request: &ControlRequest, requester: &str) -> Result<(), ControlError> {
        if request.schema_version != 1
            || !valid_label(&request.id)
            || !valid_label(&request.calc_id)
            || !valid_label(requester)
            || request.reason.len() > 512
            || request.preconditions.len() > 8
            || (request.expiry.max_iteration.is_none() && request.expiry.deadline.is_none())
            || matches!(&request.action, ControlAction::Policy(c) if c.len() > 10)
            || request
                .preconditions
                .iter()
                .any(|p| matches!(p, Precondition::Scalar { value, .. } if !value.is_finite()))
        {
            return Err(ControlError::InvalidSchema);
        }
        Ok(())
    }
    /// Admission is not semantic acceptance. Pending proposals do not stop SCF.
    /// Malformed/oversized input rejects before admission and consumes no journal.
    pub fn submit(
        &mut self,
        request: ControlRequest,
        requester: String,
        authority: Authority,
    ) -> Result<(), ControlError> {
        self.healthy()?;
        Self::schema(&request, &requester)?;
        for p in &self.pending {
            if p.request.id == request.id {
                return if p.request == request && p.requester == requester {
                    Ok(())
                } else {
                    Err(ControlError::Conflict)
                };
            }
        }
        if let Some(record) = self.audit.iter().rev().find(|r| r.request.id == request.id) {
            if record.request != request || record.requester != requester {
                return Err(ControlError::Conflict);
            }
            return match record.outcome {
                AuditOutcome::Rejected(e) => Err(e),
                _ => Ok(()),
            };
        }
        if self.pending.len() >= self.config.max_pending
            || self.audit.len() + 2 * (self.pending.len() + 1) > self.config.max_audit_records
        {
            return Err(ControlError::Capacity);
        }
        let deadline = request
            .expiry
            .deadline
            .map(|t| {
                Instant::now()
                    .checked_add(t.duration_since(SystemTime::now()).unwrap_or_default())
                    .ok_or(ControlError::InvalidSchema)
            })
            .transpose()?;
        let approval_deadline = match &authority {
            Authority::Approved(a) => Some(monotonic_deadline(a.expires_at)?),
            _ => None,
        };
        if self.state == RuntimeState::Running && matches!(request.action, ControlAction::Pause) {
            self.state = transition(self.state, Trigger::Pause)?;
        }
        self.pending.push(Pending {
            request,
            requester,
            authority,
            deadline,
            approval_deadline,
            replayed: false,
        });
        Ok(())
    }
    /// Host-authenticated approval binds the already-owned payload. Reusing an
    /// id with different parameters is a conflict even after rejection/application.
    pub fn approve(&mut self, id: &str, approval: Approval) -> Result<(), ControlError> {
        self.healthy()?;
        let p = self
            .pending
            .iter_mut()
            .find(|p| p.request.id == id)
            .ok_or(ControlError::Conflict)?;
        p.approval_deadline = Some(monotonic_deadline(approval.expires_at)?);
        p.authority = Authority::Approved(approval);
        Ok(())
    }
    pub fn deny(&mut self, id: &str) -> Result<(), ControlError> {
        self.healthy()?;
        if matches!(
            self.state,
            RuntimeState::Running | RuntimeState::PauseRequested | RuntimeState::Restarting
        ) {
            return Err(ControlError::InvalidTransition(InvalidTransition {
                from: self.state,
                trigger: Trigger::Complete,
            }));
        }
        let i = self
            .pending
            .iter()
            .position(|p| p.request.id == id)
            .ok_or(ControlError::Conflict)?;
        let p = self.pending.remove(i);
        self.reject(&p, ControlError::ApprovalDenied)
    }
    /// Complete-boundary arbitration: approved cancellation, approved controls,
    /// then autonomous proposals; stable receipt order within each class.
    /// One numerical transaction per boundary; terminal completion wins.
    pub fn safe_point(&mut self) -> Result<(), ControlError> {
        self.healthy()?;
        if matches!(
            self.state,
            RuntimeState::Running | RuntimeState::PauseRequested | RuntimeState::Restarting
        ) {
            return Err(ControlError::InvalidTransition(InvalidTransition {
                from: self.state,
                trigger: Trigger::Complete,
            }));
        }
        if self.pending.is_empty() {
            return Ok(());
        }
        self.pending
            .sort_by_key(|p| match (&p.authority, &p.request.action) {
                (Authority::Approved(_), ControlAction::Cancel) => 0,
                (Authority::Approved(_), _) => 1,
                (Authority::Autonomous, _) => 2,
                _ => 3,
            });
        let mut numerical_applied =
            self.last_policy_iteration == Some(self.session.snapshot().iteration);
        let mut i = 0;
        while i < self.pending.len() {
            let validation = self.validate(&self.pending[i], false);
            if validation == Err(ControlError::ApprovalRequired) {
                i += 1;
                continue;
            }
            let p = self.pending.remove(i);
            let numerical = matches!(p.request.action, ControlAction::Policy(_));
            let validation = validation.and_then(|rule| {
                if numerical && numerical_applied {
                    Err(ControlError::Superseded)
                } else {
                    Ok(rule)
                }
            });
            match validation {
                Err(e) => self.reject(&p, e)?,
                Ok(rule) => match self.apply(&p, rule) {
                    Ok(()) => numerical_applied |= numerical,
                    Err(ControlError::JournalFailed) => return Err(ControlError::JournalFailed),
                    Err(e) => self.reject(&p, e)?,
                },
            }
        }
        Ok(())
    }
    fn validate(&self, p: &Pending, replay: bool) -> Result<Option<usize>, ControlError> {
        let r = &p.request;
        let s = self.session.snapshot();
        if r.calc_id != s.calc_id {
            return Err(ControlError::WrongCalculation);
        }
        if self.state.terminal() {
            return Err(ControlError::Terminal);
        }
        if !matches!(r.action, ControlAction::SwitchMode(_))
            && s.mode != ExecutionMode::AgentControlled
        {
            return Err(ControlError::ModeDisallowsControl);
        }
        if replay {
            return Ok(None);
        }
        if r.expiry.max_iteration.is_some_and(|n| s.iteration > n)
            || p.deadline.is_some_and(|t| Instant::now() >= t)
            || r.expiry.deadline.is_some_and(|t| SystemTime::now() >= t)
        {
            return Err(ControlError::Expired);
        }
        if r.observed_policy_revision != self.revision {
            return Err(ControlError::StalePolicy);
        }
        if r.observed_state_version > s.state_version {
            return Err(ControlError::FutureState);
        }
        let stale = r.observed_state_version != s.state_version;
        if stale && r.stale_policy == StalePolicy::Reject {
            return Err(ControlError::StaleState);
        }
        if !r.preconditions.iter().all(|p| self.predicate(p)) {
            return Err(ControlError::PreconditionFailed);
        }
        match &p.authority {
            Authority::PendingApproval => Err(ControlError::ApprovalRequired),
            Authority::Approved(a) => {
                if !valid_label(&a.id) || !valid_label(&a.approver) {
                    return Err(ControlError::InvalidSchema);
                }
                if SystemTime::now() >= a.expires_at
                    || p.approval_deadline.is_some_and(|t| Instant::now() >= t)
                {
                    return Err(ControlError::ApprovalExpired);
                }
                if stale && !a.allow_revalidation {
                    return Err(ControlError::StaleState);
                }
                Ok(None)
            }
            Authority::Autonomous => {
                if !self.autonomy_enabled {
                    return Err(ControlError::AutonomousDenied);
                }
                let policy = self
                    .config
                    .autonomy
                    .as_ref()
                    .ok_or(ControlError::AutonomousDenied)?;
                if self.autonomous_total >= policy.max_applied_total {
                    return Err(ControlError::Quota);
                }
                for (i, rule) in policy.rules.iter().enumerate() {
                    if !rule_matches(&rule.action, &r.action, s.options) {
                        continue;
                    }
                    if self.rule_counts[i] >= rule.max_count {
                        return Err(ControlError::Quota);
                    }
                    if self.rule_last_iteration[i]
                        .is_some_and(|n| s.iteration.saturating_sub(n) < rule.cooldown_iterations)
                    {
                        return Err(ControlError::Cooldown);
                    }
                    if stale
                        && (!rule.allow_revalidation
                            || s.state_version - r.observed_state_version > rule.max_version_lag)
                    {
                        return Err(ControlError::StaleState);
                    }
                    return Ok(Some(i));
                }
                Err(ControlError::AutonomousDenied)
            }
        }
    }
    fn predicate(&self, p: &Precondition) -> bool {
        let s = self.session.snapshot();
        match *p {
            Precondition::State(state) => state == self.state,
            Precondition::PolicyRevision(revision) => revision == self.revision,
            Precondition::Mode(mode) => mode == s.mode,
            Precondition::DiagnosticActive => false,
            Precondition::Scalar {
                field,
                comparison,
                value,
            } => {
                let measured = match field {
                    Scalar::DensityResidual => s.density_residual,
                    Scalar::ChargeErrorRelative => s.charge_error_relative,
                    Scalar::EnergyChangeAbsHartree => s.energy_change_hartree.map(f64::abs),
                };
                measured.is_some_and(|v| {
                    v.is_finite()
                        && match comparison {
                            Comparison::Less => v < value,
                            Comparison::LessOrEqual => v <= value,
                            Comparison::Greater => v > value,
                            Comparison::GreaterOrEqual => v >= value,
                        }
                })
            }
        }
    }
    fn prepare(
        &self,
        action: &ControlAction,
    ) -> Result<(Mutation, BoundarySnapshot, bool), ControlError> {
        self.session.check_control_version()?;
        let mut after = self.snapshot();
        let mut reset = false;
        let mutation = match action {
            ControlAction::Policy(changes) => {
                if self.last_policy_iteration == Some(after.iteration) {
                    return Err(ControlError::Superseded);
                }
                if !matches!(self.state, RuntimeState::Initialized | RuntimeState::Paused) {
                    return Err(ControlError::InvalidTransition(InvalidTransition {
                        from: self.state,
                        trigger: Trigger::Complete,
                    }));
                }
                let prepared = crate::policy::prepare_policy(
                    &after.options,
                    after.iteration,
                    self.config.tolerance_floors,
                    changes,
                )
                .map_err(ControlError::Policy)?;
                after.options = prepared.options().clone();
                reset = prepared.resets_history();
                after.policy_revision =
                    self.revision.checked_add(1).ok_or(ControlError::Capacity)?;
                Mutation::Policy(
                    self.session
                        .prepare_boundary_policy(changes, self.config.tolerance_floors)?,
                )
            }
            ControlAction::SwitchMode(mode) => {
                after.mode = *mode;
                if self.session.snapshot().mode == *mode {
                    Mutation::Noop
                } else {
                    Mutation::Mode(self.session.prepare_boundary_mode(*mode)?)
                }
            }
            ControlAction::Pause => {
                after.state = transition(self.state, Trigger::Pause)?;
                Mutation::State(after.state)
            }
            ControlAction::Cancel => {
                after.state = transition(self.state, Trigger::Cancel)?;
                Mutation::State(after.state)
            }
            ControlAction::Resume => {
                let restarting = transition(self.state, Trigger::Resume)?;
                after.state = transition(restarting, Trigger::ResumeComplete)?;
                Mutation::State(after.state)
            }
            ControlAction::Checkpoint => return Err(ControlError::UnsupportedCapability),
        };
        if !matches!(mutation, Mutation::Noop) {
            after.state_version += 1;
        }
        Ok((mutation, after, reset))
    }
    fn record(
        &self,
        p: &Pending,
        outcome: AuditOutcome,
        after: BoundarySnapshot,
        history_reset: bool,
        offset: usize,
    ) -> AuditRecord {
        AuditRecord {
            schema_version: 1,
            calc_id: self.session.snapshot().calc_id.into(),
            sequence: (self.audit.len() + offset + 1) as u64,
            timestamp: SystemTime::now(),
            request: p.request.clone(),
            requester: p.requester.clone(),
            approval: match &p.authority {
                Authority::Approved(a) => Some(a.clone()),
                _ => None,
            },
            authority: p.authority.clone(),
            replayed: p.replayed,
            autonomy_rule_index: None,
            before: self.snapshot(),
            after,
            outcome,
            history_reset,
        }
    }
    fn publish(&mut self, batch: Vec<AuditRecord>) -> Result<(), ControlError> {
        if self.audit.len() + batch.len() > self.config.max_audit_records {
            return Err(ControlError::Capacity);
        }
        self.audit
            .try_reserve(batch.len())
            .map_err(|_| ControlError::Capacity)?;
        if let Some(sink) = &mut self.sink
            && sink.append(&batch).is_err()
        {
            self.journal_failed = true;
            // The failed durable write cannot be reported to that sink.
            // Retain a bounded local failure disposition for the result.
            let mut failed = batch.last().expect("nonempty control batch").clone();
            failed.sequence = self.audit.len() as u64 + 1;
            failed.outcome = AuditOutcome::Rejected(ControlError::JournalFailed);
            failed.after = failed.before.clone();
            failed.history_reset = false;
            self.audit.push(failed);
            return Err(ControlError::JournalFailed);
        }
        self.audit.extend(batch);
        Ok(())
    }
    fn reject(&mut self, p: &Pending, error: ControlError) -> Result<(), ControlError> {
        let record = self.record(p, AuditOutcome::Rejected(error), self.snapshot(), false, 0);
        self.publish(vec![record])
    }
    fn apply(&mut self, p: &Pending, rule: Option<usize>) -> Result<(), ControlError> {
        let (mutation, after, reset) = self.prepare(&p.request.action)?;
        let mut accepted = self.record(p, AuditOutcome::Accepted, self.snapshot(), false, 0);
        let mut applied = self.record(p, AuditOutcome::Applied, after.clone(), reset, 1);
        accepted.autonomy_rule_index = rule;
        applied.autonomy_rule_index = rule;
        // Both records and replacement state are prepared before publication.
        // No fallible operation remains after a successful atomic sink append.
        self.publish(vec![accepted, applied])?;
        match mutation {
            Mutation::Noop => {}
            Mutation::Policy(prepared) => {
                self.session.commit_boundary_policy(prepared);
                self.last_policy_iteration = Some(after.iteration);
            }
            Mutation::Mode(prepared) => {
                self.session.commit_boundary_mode(prepared);
                if after.mode != ExecutionMode::AgentControlled {
                    self.autonomy_enabled = false;
                }
            }
            Mutation::State(state) => {
                self.state = state;
                self.session.commit_control_version();
            }
        }
        self.revision = after.policy_revision;
        if let Some(i) = rule {
            self.autonomous_total += 1;
            self.rule_counts[i] += 1;
            self.rule_last_iteration[i] = Some(after.iteration);
        }
        Ok(())
    }
    /// Replay committed decisions in journal order at recorded boundaries. This
    /// ignores wall time/approval timing but checks the complete before/after
    /// projection and revalidates numerical invariants. Same model/build/threads
    /// and initial state are required; this is not checkpoint restoration.
    pub fn replay(&mut self, records: &[AuditRecord]) -> Result<(), RuntimeError> {
        if self.state != RuntimeState::Created || !self.audit.is_empty() || !self.pending.is_empty()
        {
            return Err(ControlError::ReplayMismatch.into());
        }
        let mut sequence = 0;
        for r in records {
            if r.schema_version != 1
                || r.sequence <= sequence
                || r.calc_id != self.session.snapshot().calc_id
            {
                return Err(ControlError::ReplayMismatch.into());
            }
            sequence = r.sequence;
            if r.outcome != AuditOutcome::Applied {
                continue;
            }
            if self.state == RuntimeState::Created && r.before.state != RuntimeState::Created {
                self.initialize()?;
            }
            while self.session.snapshot().iteration < r.before.iteration {
                self.step()?;
            }
            if !same_boundary(&self.snapshot(), &r.before) {
                return Err(ControlError::ReplayMismatch.into());
            }
            Self::schema(&r.request, &r.requester)?;
            let p = Pending {
                request: r.request.clone(),
                requester: r.requester.clone(),
                authority: r.authority.clone(),
                deadline: None,
                approval_deadline: None,
                replayed: true,
            };
            self.validate(&p, true)?;
            let rule = r.autonomy_rule_index;
            if let Some(i) = rule {
                let policy = self
                    .config
                    .autonomy
                    .as_ref()
                    .ok_or(ControlError::ReplayMismatch)?;
                let grant = policy.rules.get(i).ok_or(ControlError::ReplayMismatch)?;
                if !self.autonomy_enabled
                    || !matches!(r.authority, Authority::Autonomous)
                    || self.autonomous_total >= policy.max_applied_total
                    || self.rule_counts[i] >= grant.max_count
                    || self.rule_last_iteration[i].is_some_and(|n| {
                        r.before.iteration.saturating_sub(n) < grant.cooldown_iterations
                    })
                    || !rule_matches(
                        &grant.action,
                        &r.request.action,
                        self.session.snapshot().options,
                    )
                {
                    return Err(ControlError::ReplayMismatch.into());
                }
            } else if matches!(r.authority, Authority::Autonomous) {
                return Err(ControlError::ReplayMismatch.into());
            }
            let (_, after, _) = self.prepare(&p.request.action)?;
            if !same_boundary(&after, &r.after) {
                return Err(ControlError::ReplayMismatch.into());
            }
            self.apply(&p, rule)?;
        }
        if self.state == RuntimeState::Created {
            self.initialize()?;
        }
        Ok(())
    }
}

#[derive(Debug)]
pub enum RuntimeError {
    Control(ControlError),
    Scf(ScfError),
}
impl From<ControlError> for RuntimeError {
    fn from(e: ControlError) -> Self {
        Self::Control(e)
    }
}
fn valid_label(s: &str) -> bool {
    !s.is_empty() && s.len() <= 128 && s.bytes().all(|b| b.is_ascii_graphic())
}
fn same_boundary(a: &BoundarySnapshot, b: &BoundarySnapshot) -> bool {
    a.state == b.state
        && a.mode == b.mode
        && a.iteration == b.iteration
        && a.state_version == b.state_version
        && a.policy_revision == b.policy_revision
        && format!("{:?}", a.options) == format!("{:?}", b.options)
}
fn rule_matches(rule: &AllowedAction, action: &ControlAction, options: &ScfOptions) -> bool {
    match (rule, action) {
        (AllowedAction::Pause, ControlAction::Pause)
        | (AllowedAction::Resume, ControlAction::Resume)
        | (AllowedAction::Cancel, ControlAction::Cancel) => true,
        (_, ControlAction::Policy(changes)) if changes.len() == 1 => match (rule, &changes[0]) {
            (
                AllowedAction::Beta {
                    minimum,
                    maximum,
                    only_decrease,
                },
                PolicyChange::MixingBeta(beta),
            ) => {
                let old = match options.mixing {
                    crate::Mixing::Linear { beta }
                    | crate::Mixing::Pulay { beta, .. }
                    | crate::Mixing::KerkerPulay { beta, .. } => beta,
                };
                beta >= minimum && beta <= maximum && (!only_decrease || *beta <= old)
            }
            (AllowedAction::ResetHistory, PolicyChange::ResetHistory) => true,
            (AllowedAction::Eigensolver(a), PolicyChange::Eigensolver(b)) => a == b,
            (
                AllowedAction::TightenTolerances,
                PolicyChange::DensityTolerance(_)
                | PolicyChange::EnergyToleranceHartree(_)
                | PolicyChange::EigenToleranceHartree(_),
            ) => true,
            (AllowedAction::ExtendBudget { maximum }, PolicyChange::MaxIterations(n)) => {
                n <= maximum
            }
            _ => false,
        },
        _ => false,
    }
}

fn monotonic_deadline(wall: SystemTime) -> Result<Instant, ControlError> {
    Instant::now()
        .checked_add(wall.duration_since(SystemTime::now()).unwrap_or_default())
        .ok_or(ControlError::InvalidSchema)
}
