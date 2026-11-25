//! Unix proposal endpoint; only the launch terminal grants operator approval.
use super::*;
use dft_engine::{policy::PolicyChange, runtime::*};
use serde::{Deserialize, Serialize};
use std::{
    io::{BufRead, IsTerminal},
    os::unix::{fs::PermissionsExt, net::UnixDatagram},
    sync::mpsc::Receiver,
    time::{Duration, SystemTime},
};

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub(super) struct Proposal {
    id: String,
    calc_id: String,
    observed_state_version: u64,
    observed_policy_revision: u64,
    max_iteration: usize,
    action: String,
    beta: Option<f64>,
    history: Option<usize>,
    threshold_ry: Option<f64>,
    verification: Option<bool>,
}
impl Proposal {
    fn request(self) -> Result<ControlRequest, String> {
        if self.history.is_some() || self.threshold_ry.is_some() || self.verification.is_some() {
            return Err(
                "legacy control channel supports only pause/resume/cancel/mixing_beta".into(),
            );
        }
        let action = match self.action.as_str() {
            "pause" if self.beta.is_none() => ControlAction::Pause,
            "resume" if self.beta.is_none() => ControlAction::Resume,
            "cancel" if self.beta.is_none() => ControlAction::Cancel,
            "mixing_beta" => {
                let beta = self.beta.filter(|x|x.is_finite() && *x > 0.0 && *x <= 1.0)
                    .ok_or("beta: valid range is finite (0, 1]")?;
                ControlAction::Policy(vec![PolicyChange::MixingBeta(beta)])
            }
            _ => return Err("action: valid values are pause, resume, cancel, mixing_beta; beta is only valid for mixing_beta".into()),
        };
        Ok(ControlRequest {
            schema_version: 1,
            calc_id: self.calc_id,
            id: self.id,
            observed_state_version: self.observed_state_version,
            observed_policy_revision: self.observed_policy_revision,
            action,
            preconditions: vec![],
            expiry: Expiry {
                max_iteration: Some(self.max_iteration),
                deadline: None,
            },
            stale_policy: StalePolicy::Reject,
            reason: "local supervised CLI proposal".into(),
        })
    }
}
pub struct Channel {
    socket: UnixDatagram,
    path: PathBuf,
    approvals: Receiver<String>,
    pending: std::collections::BTreeSet<String>,
}
impl Channel {
    pub fn new(path: &Path) -> Result<Self, RunError> {
        // A private directory closes the bind/chmod permission race.
        let parent = path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        if std::fs::metadata(parent)?.permissions().mode() & 0o077 != 0 {
            return Err(fail(format!(
                "control_socket: parent {} must have permissions 0700; use --control-socket in a private directory",
                parent.display()
            )));
        }
        let socket = UnixDatagram::bind(path)?;
        let (sender, approvals) = mpsc::sync_channel(16);
        let channel = Self {
            socket,
            path: path.into(),
            approvals,
            pending: Default::default(),
        };
        std::fs::set_permissions(path, std::fs::Permissions::from_mode(0o600))?;
        channel.socket.set_nonblocking(true)?;
        if io::stdin().is_terminal() {
            thread::spawn(move || {
                for line in io::stdin().lock().lines().take_while(Result::is_ok) {
                    if sender.send(line.unwrap()).is_err() {
                        break;
                    }
                }
            });
        }
        eprintln!(
            "Control channel: {}. Supervised proposals only. Operator terminal: approve ID or deny ID. Non-terminal stdin cannot approve; absent approval continues current policy.",
            path.display()
        );
        Ok(channel)
    }
    pub fn poll(&mut self, runtime: &mut Runtime) -> Result<(), RunError> {
        let mut bytes = [0; 4097];
        // Bound intake per boundary, including malformed traffic.
        for _ in 0..16 {
            match self.socket.recv_from(&mut bytes) {
                Ok((n, peer)) => {
                    let response = (|| -> Result<Value,String> {
                        if n > 4096 { return Err("proposal: maximum 4096 bytes".into()); }
                        let p: Proposal = serde_json::from_slice(&bytes[..n]).map_err(|e|format!("{e:?}"))?;
                        let request = p.request()?;
                        let id = request.id.clone();
                        let display = format!("{request:?}");
                        runtime.submit(request, "local-agent".into(), Authority::PendingApproval).map_err(|e|format!("{e:?}"))?;
                        self.pending.insert(id.clone());
                        // Do not render untrusted raw strings as terminal escape sequences.
                        eprintln!("Pending proposal {id:?}: {display}. Type approve ID or deny ID; no automatic approval");
                        Ok(json!({"id":id,"status":"pending_approval"}))
                    })().unwrap_or_else(|e|json!({"error":e}));
                    if let Some(path) = peer.as_pathname() {
                        // Replies are best effort. Never wait for the external agent.
                        let _ = self
                            .socket
                            .send_to(&serde_json::to_vec(&response).map_err(fail)?, path);
                    }
                }
                Err(e) if e.kind() == io::ErrorKind::WouldBlock => break,
                Err(e) => return Err(fail(e)),
            }
        }
        for _ in 0..16 {
            let Ok(line) = self.approvals.try_recv() else {
                break;
            };
            let words: Vec<_> = line.split_whitespace().collect();
            if words.len() != 2 || !self.pending.contains(words[1]) {
                eprintln!("approval: expected approve ID or deny ID for a pending proposal");
                continue;
            }
            let id = words[1];
            let result = match words[0] {
                "approve" => runtime.approve(
                    id,
                    Approval {
                        id: format!("terminal-{id}"),
                        approver: "launch-terminal".into(),
                        expires_at: SystemTime::now() + Duration::from_secs(60),
                        allow_revalidation: false,
                    },
                ),
                "deny" => runtime.deny(id),
                _ => {
                    eprintln!("approval: valid commands are approve ID, deny ID");
                    continue;
                }
            };
            if let Err(e) = result {
                eprintln!("approval {id:?}: {e:?}");
            }
            self.pending.remove(id);
        }
        Ok(())
    }
}
impl Drop for Channel {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.path);
    }
}

pub struct Journal(pub File);
impl AuditSink for Journal {
    fn append(&mut self, records: &[AuditRecord]) -> Result<(), JournalError> {
        use std::io::{Seek, SeekFrom};
        let start = self.0.stream_position().map_err(|_| JournalError)?;
        let batch = records.iter().map(|r| json!({"schema_version":1,"calc_id":r.calc_id,"sequence":r.sequence,
            "request_id":r.request.id,"action":format!("{:?}",r.request.action),"requester":r.requester,
            "observed_state_version":r.request.observed_state_version,"observed_policy_revision":r.request.observed_policy_revision,
            "max_iteration":r.request.expiry.max_iteration,"outcome":format!("{:?}",r.outcome),
            "approver":r.approval.as_ref().map(|a|&a.approver),"iteration":r.after.iteration,
            "state":format!("{:?}",r.after.state),"state_version":r.after.state_version,"policy_revision":r.after.policy_revision}));
        let result = (|| -> io::Result<()> {
            for value in batch {
                serde_json::to_writer(&mut self.0, &value)?;
                self.0.write_all(b"\n")?;
            }
            self.0.sync_data()
        })();
        if result.is_err() {
            let _ = self.0.set_len(start);
            let _ = self.0.seek(SeekFrom::Start(start));
            let _ = self.0.sync_data();
            return Err(JournalError);
        }
        Ok(())
    }
}

#[cfg(test)]
use crate::test_fixtures as fixture;

#[cfg(test)]
mod tests {
    use super::*;
    use dft_engine::{
        SpeciesPotential,
        input::{Atom, parse_upf},
    };

    #[test]
    fn proposal_requires_separate_approval_and_rejects_stale_version() {
        let config = crate::config::Config::parse("[system]\nname='test'\nlattice_unit='bohr'\nlattice_vectors=[[4.0,0.0,0.0],[0.0,4.0,0.0],[0.0,0.0,4.0]]\n[[system.atoms]]\nspecies='Si'\nfractional=[0.13,0.21,0.32]\n[pseudopotentials]\nSi='synthetic.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n").unwrap();
        let prepared = config.prepare().unwrap();
        let upf = parse_upf(&fixture::fixture()).unwrap();
        let problem = ScfProblem::new(
            &prepared.cell,
            &[Atom::new("Si", [0.13, 0.21, 0.32])],
            &[SpeciesPotential {
                label: "Si",
                potential: &upf,
            }],
            crate::run::options(&config).unwrap(),
        )
        .unwrap();
        let session = ScfSession::new(
            problem,
            CompletionMode::RetainLast,
            ExecutionConfig {
                mode: ExecutionMode::AgentControlled,
                calc_id: "test".into(),
                ..ExecutionConfig::default()
            },
        )
        .unwrap();
        let root = PathBuf::from(format!("/tmp/p11h-control-{}", std::process::id()));
        std::fs::create_dir(&root).unwrap();
        std::fs::set_permissions(&root, std::fs::Permissions::from_mode(0o700)).unwrap();
        let (socket, endpoint) = UnixDatagram::pair().unwrap();
        endpoint.set_nonblocking(true).unwrap();
        // Inject operator decisions via the separate trusted terminal side.
        let (operator, receiver) = mpsc::sync_channel(16);
        let mut channel = Channel {
            socket: endpoint,
            path: root.join("control.sock"),
            approvals: receiver,
            pending: Default::default(),
        };
        let audit = File::create(root.join("audit.jsonl")).unwrap();
        let mut runtime = Runtime::new(session, RuntimeConfig::default())
            .unwrap()
            .with_audit_sink(Box::new(Journal(audit)))
            .unwrap();
        runtime.initialize().unwrap();

        let proposal = |id: &str, action: &str, version: u64| json!({"id":id,"calc_id":"test","observed_state_version":version,"observed_policy_revision":0,"max_iteration":20,"action":action});
        socket
            .send(&serde_json::to_vec(&proposal("pause-1", "pause", 1)).unwrap())
            .unwrap();
        channel.poll(&mut runtime).unwrap();
        runtime.safe_point().unwrap();
        assert_eq!(runtime.state(), RuntimeState::Initialized);
        assert!(runtime.audit().is_empty());
        operator.send("approve pause-1".into()).unwrap();
        channel.poll(&mut runtime).unwrap();
        runtime.safe_point().unwrap();
        assert_eq!(runtime.state(), RuntimeState::Paused);
        let version = runtime.snapshot().state_version;
        socket
            .send(&serde_json::to_vec(&proposal("resume-1", "resume", version)).unwrap())
            .unwrap();
        channel.poll(&mut runtime).unwrap();
        operator.send("approve resume-1".into()).unwrap();
        channel.poll(&mut runtime).unwrap();
        runtime.safe_point().unwrap();
        assert_eq!(runtime.state(), RuntimeState::Initialized);
        socket
            .send(&serde_json::to_vec(&proposal("stale", "pause", 1)).unwrap())
            .unwrap();
        channel.poll(&mut runtime).unwrap();
        operator.send("approve stale".into()).unwrap();
        channel.poll(&mut runtime).unwrap();
        runtime.safe_point().unwrap();
        assert_eq!(runtime.state(), RuntimeState::Initialized);
        assert!(matches!(
            runtime.audit().last().unwrap().outcome,
            AuditOutcome::Rejected(_)
        ));
        let audit = std::fs::read_to_string(root.join("audit.jsonl")).unwrap();
        assert!(audit.contains("launch-terminal") && audit.contains("Rejected"));
        drop(channel);
        assert!(!root.join("control.sock").exists());
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn proposal_cannot_self_approve_or_hide_invalid_beta() {
        let value = json!({"id":"x","calc_id":"test","observed_state_version":1,"observed_policy_revision":0,"max_iteration":2,"action":"mixing_beta","beta":2.0});
        let p: Proposal = serde_json::from_value(value.clone()).unwrap();
        assert!(p.request().unwrap_err().contains("(0, 1]"));
        let mut forged = value;
        forged["approved"] = json!(true);
        assert!(serde_json::from_value::<Proposal>(forged).is_err());
    }
}

/// Reference transport keeps the same supervised socket protocol. Admission
/// never grants authority; only lines from the launch terminal approve.
/// The bounded ID ledger prevents a retry from applying twice.
pub(super) struct ReferenceChannel {
    channel: Channel,
    requests: std::collections::BTreeMap<String, (Proposal, &'static str, std::time::Instant)>,
    journal: File,
}
impl Proposal {
    fn reference_action(&self) -> Result<ReferenceAction, String> {
        use dft_engine::refscf::runtime::{Control, PolicyPatch};
        let fields = usize::from(self.beta.is_some())
            + usize::from(self.history.is_some())
            + usize::from(self.threshold_ry.is_some())
            + usize::from(self.verification.is_some());
        let patch = match self.action.as_str() {
            "pause" | "resume" | "cancel" if fields == 0 => return Ok(ReferenceAction::Control(match self.action.as_str() {
                "pause" => Control::Pause, "resume" => Control::Resume, _ => Control::Cancel,
            })),
            "mixing_beta" if fields == 1 && self.beta.is_some_and(|v| v.is_finite() && v > 0. && v <= 1.) => PolicyPatch { beta: self.beta, ..Default::default() },
            "mixing_history" if fields == 1 && self.history.is_some_and(|v| (1..=25).contains(&v)) => PolicyPatch { history: self.history, ..Default::default() },
            "threshold_ry" if fields == 1 && self.threshold_ry.is_some_and(|v| v.is_finite() && v > 0.) => PolicyPatch { threshold_ry: self.threshold_ry, ..Default::default() },
            "verification" if fields == 1 && self.verification.is_some() => PolicyPatch { verification: self.verification, ..Default::default() },
            _ => return Err("reference control: supported actions are pause/resume/cancel, mixing_beta (0,1], mixing_history 1..=25, threshold_ry > 0 (tightening only), verification boolean; legacy knobs and extra parameters are rejected".into()),
        };
        Ok(ReferenceAction::Policy(patch))
    }
    fn validate_reference(
        &self,
        calc_id: &str,
        snapshot: dft_engine::refscf::runtime::Snapshot,
    ) -> Result<(), String> {
        if self.id.is_empty()
            || self.id.len() > 128
            || !self.id.bytes().all(|b| b.is_ascii_graphic())
        {
            return Err("invalid request id".into());
        }
        if self.calc_id != calc_id {
            return Err("wrong calculation".into());
        }
        if self.observed_state_version != snapshot.state_version {
            return Err("stale/future state version".into());
        }
        if self.observed_policy_revision != snapshot.policy_revision {
            return Err("stale policy revision".into());
        }
        if snapshot.last_step.map_or(0, |r| r.iteration) > self.max_iteration {
            return Err("expired iteration".into());
        }
        if snapshot.state.terminal() || !snapshot.safe_point {
            return Err("requires live safe point".into());
        }
        self.reference_action()?;
        Ok(())
    }
}
enum ReferenceAction {
    Control(dft_engine::refscf::runtime::Control),
    Policy(dft_engine::refscf::runtime::PolicyPatch),
}
impl ReferenceChannel {
    pub(super) fn new(path: &Path, out: &Path) -> Result<Self, RunError> {
        Ok(Self {
            channel: Channel::new(path)?,
            requests: Default::default(),
            journal: OpenOptions::new()
                .create_new(true)
                .write(true)
                .open(out.join("control-audit.jsonl"))?,
        })
    }
    fn record(&mut self, proposal: &Proposal, status: &str) -> Result<(), RunError> {
        serde_json::to_writer(
            &mut self.journal,
            &json!({"schema_version":1,"driver":"reference","request":proposal,"outcome":status}),
        )
        .map_err(fail)?;
        self.journal.write_all(b"\n")?;
        self.journal.sync_data()?;
        Ok(())
    }
    pub(super) fn poll(
        &mut self,
        calc_id: &str,
        snapshot: dft_engine::refscf::runtime::Snapshot,
    ) -> Result<Vec<Proposal>, RunError> {
        let mut bytes = [0; 4097];
        for _ in 0..16 {
            match self.channel.socket.recv_from(&mut bytes) {
                Ok((n, peer)) => {
                    let response = (|| -> Result<Value, String> {
                        if n > 4096 {
                            return Err("proposal exceeds 4096 bytes".into());
                        }
                        let proposal: Proposal =
                            serde_json::from_slice(&bytes[..n]).map_err(|e| e.to_string())?;
                        if let Some((old, status, _)) = self.requests.get(&proposal.id) {
                            return if old == &proposal {
                                Ok(json!({"id":proposal.id,"status":status}))
                            } else {
                                Err("request id conflict".into())
                            };
                        }
                        proposal.validate_reference(calc_id, snapshot)?;
                        if self.requests.len() >= 1024 {
                            return Err("reference control request quota exhausted".into());
                        }
                        self.record(&proposal, "pending_approval")
                            .map_err(|e| e.cause)?;
                        let response = json!({"id":proposal.id,"status":"pending_approval"});
                        eprintln!(
                            "Pending reference proposal {:?}: {:?}; approve ID or deny ID",
                            proposal.id, proposal.action
                        );
                        self.requests.insert(
                            proposal.id.clone(),
                            (proposal, "pending_approval", std::time::Instant::now()),
                        );
                        Ok(response)
                    })()
                    .unwrap_or_else(|e| json!({"error":e}));
                    if let Some(path) = peer.as_pathname() {
                        let _ = self
                            .channel
                            .socket
                            .send_to(&serde_json::to_vec(&response).map_err(fail)?, path);
                    }
                }
                Err(e) if e.kind() == io::ErrorKind::WouldBlock => break,
                Err(e) => return Err(fail(e)),
            }
        }
        let mut approved = Vec::new();
        for _ in 0..16 {
            let Ok(line) = self.channel.approvals.try_recv() else {
                break;
            };
            let words: Vec<_> = line.split_whitespace().collect();
            if words.len() != 2 || !matches!(words[0], "approve" | "deny") {
                continue;
            }
            let Some((proposal, status, admitted)) = self.requests.get(words[1]).cloned() else {
                continue;
            };
            if status != "pending_approval" {
                continue;
            }
            let disposition = if words[0] == "deny" {
                "denied"
            } else if admitted.elapsed() > Duration::from_secs(60) {
                "expired"
            } else if proposal.validate_reference(calc_id, snapshot).is_err() {
                "rejected_stale"
            } else {
                "approved"
            };
            // Durable intent precedes the mutation on every rank.
            self.record(&proposal, disposition)?;
            self.requests.get_mut(words[1]).unwrap().1 = disposition;
            if disposition == "approved" {
                approved.push(proposal);
            }
        }
        Ok(approved)
    }
    pub(super) fn disposition(
        &mut self,
        proposal: &Proposal,
        result: &Result<Value, String>,
    ) -> Result<(), RunError> {
        let status = if result.is_ok() {
            "applied"
        } else {
            "rejected"
        };
        self.requests.get_mut(&proposal.id).expect("admitted").1 = status;
        serde_json::to_writer(
            &mut self.journal,
            &json!({"schema_version":1,"driver":"reference",
            "request":proposal,"outcome":status,"result":result,"approver":"launch-terminal"}),
        )
        .map_err(fail)?;
        self.journal.write_all(b"\n")?;
        self.journal.sync_data()?;
        Ok(())
    }
    pub(super) fn close(&mut self) -> Result<(), RunError> {
        let pending: Vec<_> = self
            .requests
            .values()
            .filter(|(_, status, _)| *status == "pending_approval")
            .map(|(p, _, _)| p.clone())
            .collect();
        for proposal in pending {
            self.record(&proposal, "terminal_without_approval")?;
        }
        Ok(())
    }
}
pub(super) fn apply_reference<T: dft_engine::refscf::runtime::Completion>(
    proposal: &Proposal,
    calc_id: &str,
    runtime: &mut dft_engine::refscf::runtime::Runtime<'_, T>,
    verification_required: bool,
) -> Result<Value, String> {
    proposal.validate_reference(calc_id, runtime.snapshot())?;
    match proposal.reference_action()? {
        ReferenceAction::Control(control) => runtime.control(control).map_err(|e| e.to_string())?,
        ReferenceAction::Policy(patch) => {
            if verification_required && patch.verification == Some(false) {
                return Err("verification is required for requested derivatives".into());
            }
            runtime.apply_policy(patch).map_err(|e| e.to_string())?;
        }
    }
    let after = runtime.snapshot();
    Ok(
        json!({"state":format!("{:?}",after.state),"state_version":after.state_version,"policy_revision":after.policy_revision}),
    )
}

#[cfg(test)]
mod reference_channel_tests {
    use super::*;
    use dft_engine::{
        GammaReal, ScfOptions, SpeciesPotential,
        refscf::{driver, profile::Profile, runtime::Runtime as ReferenceRuntime},
        refstart::StartInput,
    };
    #[test]
    fn reference_socket_requires_approval_and_rejects_replay_and_stale_requests() {
        let _scope = dft_numerics::switches::RunConfig::from_pairs(&[])
            .unwrap()
            .enter();
        let cell =
            dft_geometry::Cell::from_user_rows_bohr([[5., 0., 0.], [0., 5., 0.], [0., 0., 5.]])
                .unwrap();
        let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
        let mut text = fixture::gaussian_with(2., 1., true, 161);
        let parsed = dft_upf::parse_upf(&text).unwrap();
        let radial = parsed
            .mesh()
            .r()
            .iter()
            .map(|r| format!("{:.17e}", r * (-r * r).exp()))
            .collect::<Vec<_>>()
            .join(" ");
        let start = text.find("<PP_RHOATOM>").unwrap() + "<PP_RHOATOM>".len();
        let end = text.find("</PP_RHOATOM>").unwrap();
        text.replace_range(start..end, &radial);
        text=text.replace("number_of_wfc=\"0\"","number_of_wfc=\"1\"").replace("</UPF>",&format!("<PP_PSWFC><PP_CHI.1 l=\"0\" label=\"s\" occupation=\"2\">{radial}</PP_CHI.1></PP_PSWFC></UPF>"));
        let pp = dft_upf::parse_upf(&text).unwrap();
        let potentials = [SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let options = ScfOptions {
            cutoff_hartree: 1.,
            bands: Some(4),
            gamma_real: GammaReal::On,
            ..Default::default()
        };
        let mut runtime = ReferenceRuntime::new(driver::stepper(
            StartInput {
                cell: &cell,
                atoms: &atoms,
                potentials: &potentials,
                options: &options,
            },
            Profile::default(),
        ));
        runtime.initialize().unwrap();
        let root = std::env::temp_dir().join(format!("r8c-channel-{}", std::process::id()));
        std::fs::create_dir_all(&root).unwrap();
        let (peer, socket) = UnixDatagram::pair().unwrap();
        socket.set_nonblocking(true).unwrap();
        let (operator, approvals) = mpsc::sync_channel(16);
        let mut channel = ReferenceChannel {
            channel: Channel {
                socket,
                path: root.join("unused.sock"),
                approvals,
                pending: Default::default(),
            },
            requests: Default::default(),
            journal: File::create(root.join("audit.jsonl")).unwrap(),
        };
        let proposal = json!({"id":"pause","calc_id":"test","observed_state_version":1,"observed_policy_revision":0,"max_iteration":100,"action":"pause"});
        peer.send(&serde_json::to_vec(&proposal).unwrap()).unwrap();
        assert!(channel.poll("test", runtime.snapshot()).unwrap().is_empty());
        assert_eq!(runtime.snapshot().state, RuntimeState::Initialized);
        operator.send("approve pause".into()).unwrap();
        let approved = channel.poll("test", runtime.snapshot()).unwrap();
        assert_eq!(approved.len(), 1);
        let result = apply_reference(&approved[0], "test", &mut runtime, false);
        assert!(result.is_ok());
        channel.disposition(&approved[0], &result).unwrap();
        assert_eq!(runtime.snapshot().state, RuntimeState::Paused);
        peer.send(&serde_json::to_vec(&proposal).unwrap()).unwrap();
        operator.send("approve pause".into()).unwrap();
        assert!(channel.poll("test", runtime.snapshot()).unwrap().is_empty());
        let mut stale = proposal.clone();
        stale["id"] = json!("stale");
        peer.send(&serde_json::to_vec(&stale).unwrap()).unwrap();
        assert!(channel.poll("test", runtime.snapshot()).unwrap().is_empty());
        assert!(!channel.requests.contains_key("stale"));
        for (id, action, param, value) in [
            ("beta", "mixing_beta", "beta", json!(0.3)),
            ("verify", "verification", "verification", json!(false)),
            ("bad", "max_iterations", "history", json!(20)),
        ] {
            let mut p = proposal.clone();
            p["id"] = json!(id);
            p["action"] = json!(action);
            p[param] = value;
            p["observed_state_version"] = json!(runtime.snapshot().state_version);
            p["observed_policy_revision"] = json!(runtime.snapshot().policy_revision);
            peer.send(&serde_json::to_vec(&p).unwrap()).unwrap();
            assert!(channel.poll("test", runtime.snapshot()).unwrap().is_empty());
            if id == "bad" {
                assert!(!channel.requests.contains_key(id));
                continue;
            }
            operator.send(format!("approve {id}")).unwrap();
            let approved = channel.poll("test", runtime.snapshot()).unwrap();
            assert_eq!(approved.len(), 1);
            let result = apply_reference(&approved[0], "test", &mut runtime, id == "verify");
            assert_eq!(result.is_ok(), id == "beta");
            channel.disposition(&approved[0], &result).unwrap();
        }
        assert_eq!(runtime.snapshot().policy_revision, 1);
        runtime
            .control(dft_engine::refscf::runtime::Control::Cancel)
            .unwrap();
        channel.close().unwrap();
        drop(channel);
        let journal = std::fs::read_to_string(root.join("audit.jsonl")).unwrap();
        assert!(
            journal.contains("pending_approval")
                && journal.contains("approved")
                && journal.contains("applied")
        );
        std::fs::remove_dir_all(root).unwrap();
    }
}
