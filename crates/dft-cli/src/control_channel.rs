//! Unix proposal endpoint; only the launch terminal grants operator approval.
use super::*;
use dft_engine::{policy::PolicyChange, runtime::*};
use serde::Deserialize;
use std::{
    io::{BufRead, IsTerminal},
    os::unix::{fs::PermissionsExt, net::UnixDatagram},
    sync::mpsc::Receiver,
    time::{Duration, SystemTime},
};

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Proposal {
    id: String,
    calc_id: String,
    observed_state_version: u64,
    observed_policy_revision: u64,
    max_iteration: usize,
    action: String,
    beta: Option<f64>,
}
impl Proposal {
    fn request(self) -> Result<ControlRequest, String> {
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
