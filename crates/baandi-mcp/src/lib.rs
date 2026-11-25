//! Minimal MCP adapter over the existing safe-point runtime. No filesystem tools.
mod reference;
pub use reference::ReferenceInput;
use reference::{Action, Worker};
pub mod schema;
mod transport;
use dft_engine::{
    refscf::{
        runtime::{Control, Snapshot},
        stepper::StepReport,
    },
    runtime::{Approval, RuntimeState},
};
use serde_json::{Value, json};
use std::{
    collections::BTreeMap,
    path::{Path, PathBuf},
    time::SystemTime,
};

#[derive(Clone, Debug)]
pub struct Quotas {
    pub calculations: usize,
    pub requests: usize,
    pub total_steps: usize,
    pub frame_bytes: usize,
}
impl Default for Quotas {
    fn default() -> Self {
        Self {
            calculations: 8,
            requests: 4096,
            total_steps: 1000,
            frame_bytes: 65536,
        }
    }
}
/// Host owns preparation costs and must preflight basis/grid/band memory before
/// registering. Agent inputs can only select these bounded, single-use objects.
pub struct Server {
    workspace: PathBuf,
    quotas: Quotas,
    inputs: BTreeMap<String, ReferenceInput>,
    runtimes: BTreeMap<String, Calculation>,
    requests: usize,
    steps: usize,
    initialized: bool,
    ready: bool,
}
fn error(e: impl std::fmt::Debug) -> String {
    format!("{e:?}")
}
impl Server {
    pub fn new(workspace: impl AsRef<Path>, quotas: Quotas) -> Result<Self, String> {
        let workspace = workspace.as_ref().canonicalize().map_err(error)?;
        if !workspace.is_dir()
            || quotas.calculations == 0
            || quotas.calculations > 64
            || quotas.requests == 0
            || quotas.frame_bytes < 1024
            || quotas.frame_bytes > 1_048_576
        {
            return Err("invalid workspace or quotas".into());
        }
        Ok(Self {
            workspace,
            quotas,
            inputs: BTreeMap::new(),
            runtimes: BTreeMap::new(),
            requests: 0,
            steps: 0,
            initialized: false,
            ready: false,
        })
    }
    pub fn workspace(&self) -> &Path {
        &self.workspace
    }
    pub fn register_input(&mut self, id: &str, problem: ReferenceInput) -> Result<(), String> {
        problem.request.profile.validate().map_err(error)?;
        if id.is_empty()
            || id.len() > 128
            || !id.bytes().all(|b| b.is_ascii_graphic())
            || self.inputs.contains_key(id)
            || self.inputs.len() + self.runtimes.len() >= self.quotas.calculations
        {
            return Err("input identity or capacity".into());
        }
        self.inputs.insert(id.into(), problem);
        Ok(())
    }
    pub fn runtime(&self, id: &str) -> Option<&Calculation> {
        self.runtimes.get(id)
    }
    /// Only the authenticated embedding host calls this, outside MCP dispatch.
    pub fn approve(&mut self, calc: &str, request: &str, approval: Approval) -> Result<(), String> {
        let r = self.runtimes.get_mut(calc).ok_or("unknown calculation")?;
        let pending = r.pending.remove(request).ok_or("unknown pending request")?;
        let before = r.snapshot();
        let result = if approval.id.is_empty()
            || approval.approver.is_empty()
            || approval.expires_at <= SystemTime::now()
        {
            Err("invalid or expired host approval".into())
        } else {
            r.validate(&pending.control)
                .and_then(|_| r.worker.apply(pending.action))
        };
        r.record(
            &pending.control,
            pending.action,
            before,
            &result,
            Some(&approval),
        );
        result
    }
    pub fn deny(&mut self, calc: &str, request: &str) -> Result<(), String> {
        let r = self.runtimes.get_mut(calc).ok_or("unknown calculation")?;
        let pending = r.pending.remove(request).ok_or("unknown pending request")?;
        r.record(
            &pending.control,
            pending.action,
            r.snapshot(),
            &Err("denied by host".into()),
            None,
        );
        Ok(())
    }
    fn snapshot(r: &Calculation) -> Value {
        let s = r.snapshot();
        json!({"schema_version":2,"driver":"reference","calc_id":r.id,"state":format!("{:?}",s.state),"mode":format!("{:?}",r.mode),"iteration":s.last_step.map_or(0, |r|r.iteration),"completed_maps":s.completed_maps,"state_version":s.state_version,"policy_revision":s.policy_revision,"safe_point":s.safe_point,"energy_hartree":s.last_step.and_then(|r|r.energy_ha),"last_step":s.last_step.map(step_json)})
    }
    fn propose(&mut self, p: schema::Control, action: Action) -> Result<Value, String> {
        let r = self
            .runtimes
            .get_mut(&p.calc_id)
            .ok_or("unknown calculation")?;
        if r.pending.len() >= 16 || r.seen.len() >= 4096 || r.audit.len() >= 4096 {
            return Err("control capacity exhausted".into());
        }
        if r.seen.contains(&p.request_id) {
            return Err("duplicate request_id".into());
        }
        r.seen.insert(p.request_id.clone());
        let result = r.validate(&p);
        if result.is_err() {
            r.record(&p, action, r.snapshot(), &result, None);
        }
        let disposition = match result {
            Ok(()) => "pending_host_approval".to_string(),
            Err(e) => e,
        };
        let id = p.request_id.clone();
        if disposition == "pending_host_approval" {
            r.pending.insert(id.clone(), Pending { control: p, action });
        }
        Ok(json!({"request_id":id,"disposition":disposition,"snapshot":Self::snapshot(r)}))
    }
    fn call(&mut self, name: &str, args: Value) -> Result<Value, String> {
        match name {
            "create_calculation" => {
                let p: schema::Create = serde_json::from_value(args).map_err(error)?;
                if self.runtimes.contains_key(&p.calc_id)
                    || self.runtimes.len() >= self.quotas.calculations
                {
                    return Err("calculation identity or capacity".into());
                }
                let problem = self
                    .inputs
                    .remove(&p.input_id)
                    .ok_or("unknown or consumed validated input")?;
                let worker = Worker::new(problem)?;
                let r = Calculation {
                    id: p.calc_id.clone(),
                    mode: p.mode,
                    worker,
                    pending: BTreeMap::new(),
                    seen: Default::default(),
                    audit: vec![],
                    events: Default::default(),
                    sequence: 0,
                    dropped: 0,
                };
                let snapshot = Self::snapshot(&r);
                self.runtimes.insert(p.calc_id, r);
                Ok(snapshot)
            }
            "propose_policy_change" => {
                let p: schema::Policy = serde_json::from_value(args).map_err(error)?;
                let (p, patch) = p.into();
                self.propose(p, Action::Policy(patch))
            }
            "pause" | "resume" | "cancel" => {
                let p = serde_json::from_value(args).map_err(error)?;
                self.propose(
                    p,
                    Action::Control(match name {
                        "pause" => Control::Pause,
                        "resume" => Control::Resume,
                        _ => Control::Cancel,
                    }),
                )
            }
            "run" | "step" => {
                let (calc, steps) = if name == "run" {
                    let p: schema::Run = serde_json::from_value(args).map_err(error)?;
                    (p.calc_id, p.steps)
                } else {
                    let p: schema::Calculation = serde_json::from_value(args).map_err(error)?;
                    (p.calc_id, 1)
                };
                if steps > self.quotas.total_steps.saturating_sub(self.steps) {
                    return Err("step quota".into());
                }
                let r = self.runtimes.get_mut(&calc).ok_or("unknown calculation")?;
                for _ in 0..steps {
                    if r.state().terminal() || r.state() == RuntimeState::Paused {
                        break;
                    }
                    self.steps += 1;
                    let result = r.worker.apply(Action::Step);
                    r.sequence += 1;
                    if r.events.len() < 256 {
                        r.events.push_back(json!({"sequence":r.sequence,"snapshot":Self::snapshot(r),"error":result.as_ref().err()}));
                    } else {
                        r.dropped += 1;
                    }
                    result?;
                }
                Ok(Self::snapshot(r))
            }
            "get_events" => {
                let p: schema::Events = serde_json::from_value(args).map_err(error)?;
                let r = self
                    .runtimes
                    .get_mut(&p.calc_id)
                    .ok_or("unknown calculation")?;
                let events: Vec<_> = r.events.drain(..p.limit.min(r.events.len())).collect();
                Ok(json!({"events":events,"dropped_events":r.dropped,"audit":r.audit}))
            }
            _ => {
                let p: schema::Calculation = serde_json::from_value(args).map_err(error)?;
                let r = self.runtimes.get(&p.calc_id).ok_or("unknown calculation")?;
                let s = Self::snapshot(r);
                if name == "export_result" {
                    if !r.state().terminal() {
                        return Err("result requires terminal state".into());
                    }
                    let energy = r.snapshot().last_step.and_then(|r| r.energy_ha);
                    Ok(
                        json!({"snapshot":s,"audit":r.audit,"scope":"scalar boundary result; no density, orbitals or restart data","csv":format!("iteration,energy_hartree\n{},{}\n",r.snapshot().last_step.map_or(0, |r|r.iteration),energy.map(|e|e.to_string()).unwrap_or_default())}),
                    )
                } else {
                    Ok(s)
                }
            }
        }
    }
}

fn step_json(r: StepReport) -> Value {
    json!({"iteration":r.iteration,"kind":format!("{:?}",r.kind),"stage":r.stage,"estimator_ha":r.estimator_ha,"dr2_ry":r.dr2_ry,"energy_ha":r.energy_ha,"ethr_ry":r.ethr_ry,"h_columns":r.h_columns,"reused":r.reused,"safe_point":r.safe_point})
}
struct Pending {
    control: schema::Control,
    action: Action,
}
/// Transport approval and bounded audit state; numerical policy validation lives
/// exclusively in the reference runtime on the worker thread.
pub struct Calculation {
    id: String,
    mode: schema::Mode,
    worker: Worker,
    pending: BTreeMap<String, Pending>,
    seen: std::collections::BTreeSet<String>,
    audit: Vec<Value>,
    events: std::collections::VecDeque<Value>,
    sequence: u64,
    dropped: usize,
}
impl Calculation {
    pub fn snapshot(&self) -> Snapshot {
        self.worker.snapshot
    }
    pub fn state(&self) -> RuntimeState {
        self.snapshot().state
    }
    pub fn policy_revision(&self) -> u64 {
        self.snapshot().policy_revision
    }
    pub fn audit(&self) -> &[Value] {
        &self.audit
    }
    fn validate(&self, p: &schema::Control) -> Result<(), String> {
        let s = self.snapshot();
        if !matches!(self.mode, schema::Mode::Supervised) {
            return Err("ModeDisallowsControl".into());
        }
        if s.state.terminal() {
            return Err("terminal calculation".into());
        }
        if p.observed_state_version != s.state_version
            || p.observed_policy_revision != s.policy_revision
        {
            return Err("StaleState".into());
        }
        if s.last_step.is_some_and(|s| s.iteration > p.max_iteration) {
            return Err("request expired at max_iteration".into());
        }
        Ok(())
    }
    fn record(
        &mut self,
        p: &schema::Control,
        action: Action,
        before: Snapshot,
        result: &Result<(), String>,
        approval: Option<&Approval>,
    ) {
        let outcome = match result {
            Err(error) => error.as_str(),
            Ok(())
                if matches!(action, Action::Control(Control::Pause | Control::Cancel))
                    && !before.safe_point =>
            {
                "QueuedUntilSafePoint"
            }
            Ok(()) => "Applied",
        };
        self.audit.push(json!({"schema_version":2,"calc_id":self.id,"sequence":self.audit.len()+1,"request_id":p.request_id,"reason":p.reason,"action":format!("{action:?}"),"requester":"mcp-client","approval_id":approval.map(|a|&a.id),"approver":approval.map(|a|&a.approver),"outcome":outcome,"before_state_version":before.state_version,"after_state_version":self.snapshot().state_version,"before_policy_revision":before.policy_revision,"after_policy_revision":self.snapshot().policy_revision,"completed_maps":self.snapshot().completed_maps}));
    }
}
