//! Minimal MCP adapter over the existing safe-point runtime. No filesystem tools.
mod event;
pub mod schema;
mod transport;
use dft_engine::{
    CompletionMode, ExecutionConfig, ExecutionMode, ScfProblem, ScfSession, policy::PolicyChange,
    runtime::*,
};
use serde_json::{Value, json};
use std::{
    collections::BTreeMap,
    path::{Path, PathBuf},
    time::UNIX_EPOCH,
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
    inputs: BTreeMap<String, ScfProblem>,
    runtimes: BTreeMap<String, Runtime>,
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
    pub fn register_input(&mut self, id: &str, problem: ScfProblem) -> Result<(), String> {
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
    pub fn runtime(&self, id: &str) -> Option<&Runtime> {
        self.runtimes.get(id)
    }
    /// Only the authenticated embedding host calls this, outside MCP dispatch.
    pub fn approve(&mut self, calc: &str, request: &str, approval: Approval) -> Result<(), String> {
        let runtime = self.runtimes.get_mut(calc).ok_or("unknown calculation")?;
        runtime.approve(request, approval).map_err(error)?;
        runtime.safe_point().map_err(error)
    }
    pub fn deny(&mut self, calc: &str, request: &str) -> Result<(), String> {
        self.runtimes
            .get_mut(calc)
            .ok_or("unknown calculation")?
            .deny(request)
            .map_err(error)
    }
    fn snapshot(runtime: &Runtime) -> Value {
        let s = runtime.session().snapshot();
        json!({"schema_version":1,"calc_id":s.calc_id,"state":format!("{:?}",runtime.state()),"mode":format!("{:?}",s.mode),"iteration":s.iteration,"state_version":s.state_version,"policy_revision":runtime.policy_revision(),"energy_hartree":s.energy.map(|e|e.total_hartree()),"density_residual_bohr_minus_3_over_2":s.density_residual,"charge_error_relative":s.charge_error_relative,"max_eigen_residual_hartree":s.max_eigen_residual_hartree})
    }
    fn audit(runtime: &Runtime) -> Value {
        json!(runtime.audit().iter().map(|a|json!({"schema_version":a.schema_version,"calc_id":a.calc_id,"sequence":a.sequence,"timestamp_unix_ms":a.timestamp.duration_since(UNIX_EPOCH).unwrap_or_default().as_millis().to_string(),"request_id":a.request.id,"requester":a.requester,"reason":a.request.reason,"outcome":format!("{:?}",a.outcome),"before_state_version":a.before.state_version,"after_state_version":a.after.state_version,"record_debug":format!("{a:?}")})).collect::<Vec<_>>())
    }
    fn propose(&mut self, p: schema::Control, action: ControlAction) -> Result<Value, String> {
        let r = self
            .runtimes
            .get_mut(&p.calc_id)
            .ok_or("unknown calculation")?;
        let id = p.request_id.clone();
        r.submit(
            ControlRequest {
                schema_version: 1,
                calc_id: p.calc_id,
                id: p.request_id,
                observed_state_version: p.observed_state_version,
                observed_policy_revision: p.observed_policy_revision,
                action,
                preconditions: vec![],
                expiry: Expiry {
                    max_iteration: Some(p.max_iteration),
                    deadline: None,
                },
                stale_policy: StalePolicy::Reject,
                reason: p.reason,
            },
            "mcp-client".into(),
            Authority::PendingApproval,
        )
        .map_err(error)?;
        r.safe_point().map_err(error)?;
        let disposition = r
            .audit()
            .iter()
            .rev()
            .find(|a| a.request.id == id)
            .map(|a| format!("{:?}", a.outcome))
            .unwrap_or_else(|| "pending_host_approval".into());
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
                let session = ScfSession::new(
                    problem,
                    CompletionMode::RetainLast,
                    ExecutionConfig {
                        mode: match p.mode {
                            schema::Mode::Observe => ExecutionMode::Observed,
                            schema::Mode::Supervised => ExecutionMode::AgentControlled,
                        },
                        calc_id: p.calc_id.clone(),
                        routine_capacity: 256,
                        critical_capacity: 128,
                    },
                )
                .map_err(error)?;
                let mut r = Runtime::new(session, RuntimeConfig::default()).map_err(error)?;
                let initialized = r.initialize().map_err(error);
                let snapshot = Self::snapshot(&r);
                self.runtimes.insert(p.calc_id, r);
                initialized?;
                Ok(snapshot)
            }
            "propose_policy_change" => {
                let p: schema::Policy = serde_json::from_value(args).map_err(error)?;
                let (p, beta) = p.into();
                self.propose(
                    p,
                    ControlAction::Policy(vec![PolicyChange::MixingBeta(beta)]),
                )
            }
            "pause" | "resume" | "cancel" => {
                let p = serde_json::from_value(args).map_err(error)?;
                self.propose(
                    p,
                    match name {
                        "pause" => ControlAction::Pause,
                        "resume" => ControlAction::Resume,
                        _ => ControlAction::Cancel,
                    },
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
                    r.step().map_err(error)?;
                }
                Ok(Self::snapshot(r))
            }
            "get_events" => {
                let p: schema::Events = serde_json::from_value(args).map_err(error)?;
                let r = self
                    .runtimes
                    .get_mut(&p.calc_id)
                    .ok_or("unknown calculation")?;
                let mut events = vec![];
                if let Some(queue) = r.events() {
                    for _ in 0..p.limit {
                        let Some(e) = queue.pop() else { break };
                        events.push(json!({"schema_version":e.schema_version,"calc_id":e.calc_id.as_ref(),"sequence":e.sequence,"iteration":e.iteration,"state_version":e.state_version,"timestamp_unix_ms":e.timestamp.duration_since(UNIX_EPOCH).unwrap_or_default().as_millis().to_string(),"kind":event::payload(e.kind)}));
                    }
                }
                Ok(json!({"events":events,"audit":Self::audit(r)}))
            }
            _ => {
                let p: schema::Calculation = serde_json::from_value(args).map_err(error)?;
                let r = self.runtimes.get(&p.calc_id).ok_or("unknown calculation")?;
                let s = Self::snapshot(r);
                if name == "export_result" {
                    if !r.state().terminal() {
                        return Err("result requires terminal state".into());
                    }
                    let energy = r.session().snapshot().energy.map(|e| e.total_hartree());
                    Ok(
                        json!({"snapshot":s,"audit":Self::audit(r),"scope":"scalar boundary result; no density, orbitals or restart data","csv":format!("iteration,energy_hartree\n{},{}\n",r.session().snapshot().iteration,energy.map(|e|e.to_string()).unwrap_or_default())}),
                    )
                } else {
                    Ok(s)
                }
            }
        }
    }
}
