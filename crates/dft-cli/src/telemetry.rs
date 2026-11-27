//! CLI boundary adapters. Bounded publication keeps sink I/O out of SCF steps.
use crate::run::{RunError, engine_error};
use dft_engine::{
    ScfOutcome, ScfProblem,
    events::{EventKind, EventQueue, ExecutionConfig, ExecutionMode, ScfEvent},
    runtime::{Runtime, RuntimeConfig, RuntimeState},
    session::{CompletionMode, ScfSession, StepStatus},
};
use serde_json::{Value, json};
use std::{
    fs::{File, OpenOptions},
    io::{self, Write},
    path::{Path, PathBuf},
    sync::mpsc::{self, SyncSender},
    thread::{self, JoinHandle},
    time::UNIX_EPOCH,
};

#[derive(Default)]
pub struct RunMode {
    pub observe: bool,
    pub agent: bool,
    pub events: Option<String>,
    pub control_socket: Option<PathBuf>,
}
fn fail(cause: impl std::fmt::Display) -> RunError {
    RunError {
        code: 5,
        cause: cause.to_string(),
    }
}
impl RunMode {
    pub fn validate(&self) -> Result<(), RunError> {
        if self.observe && self.agent {
            return Err(fail("mode: choose exactly one of --observe or --agent"));
        }
        if self.events.is_some() && !self.observe && !self.agent {
            return Err(fail("events: requires --observe or --agent"));
        }
        if self.control_socket.is_some() && !self.agent {
            return Err(fail("control_socket: requires --agent"));
        }
        Ok(())
    }
}

/// Stable CLI schema, independent of Rust Debug formatting.
pub fn event_json(e: &ScfEvent) -> Value {
    let (kind, data) = match e.kind {
        EventKind::Initialized => ("initialized", json!({})),
        EventKind::IterationFinished {
            energy_hartree,
            density_residual,
            charge_error_relative,
            mixing_applied,
        } => (
            "iteration_finished",
            json!({"energy_hartree":energy_hartree,"density_residual_bohr_minus_3_over_2":density_residual,"charge_error_relative":charge_error_relative,"mixing_applied":mixing_applied,"solver_status":"iteration_completed"}),
        ),
        EventKind::CertificateEvaluated {
            passed,
            density_residual,
            occupied_eigenvalue_change_hartree,
        } => (
            "certificate_evaluated",
            json!({"passed":passed,"density_residual_bohr_minus_3_over_2":density_residual,"occupied_eigenvalue_change_hartree":occupied_eigenvalue_change_hartree}),
        ),
        EventKind::MixerReset => ("mixer_reset", json!({})),
        EventKind::Diverged => ("diverged", json!({})),
        EventKind::Converged => ("converged", json!({})),
        EventKind::NotConverged => ("not_converged", json!({})),
        EventKind::Failed { code } => ("failed", json!({"code":format!("{code:?}")})),
        EventKind::EventsDropped {
            count,
            first_missing_seq,
            last_missing_seq,
        } => (
            "events_dropped",
            json!({"count":count,"first_missing_seq":first_missing_seq,"last_missing_seq":last_missing_seq}),
        ),
    };
    json!({"schema_version":e.schema_version,"calc_id":e.calc_id.as_ref(),"sequence":e.sequence,
        "timestamp_unix_seconds":e.timestamp.duration_since(UNIX_EPOCH).map(|d|d.as_secs_f64()).unwrap_or(0.0),
        "iteration":e.iteration,"state_version":e.state_version,"type":kind,"data":data})
}

fn human(v: &Value) -> String {
    let d = &v["data"];
    if v["type"] == "iteration_finished" && d["driver"] == "reference" {
        format!(
            "iteration {} ({}): energy={} Ha; estimator={} Ha; dr2={} Ry; ethr={} Ry; H columns={}",
            v["iteration"],
            d["kind"].as_str().unwrap_or("map"),
            d["energy_hartree"],
            d["estimator_hartree"],
            d["dr2_ry"],
            d["ethr_ry"],
            d["h_columns"]
        )
    } else if v["type"] == "iteration_finished" {
        format!(
            "iteration {}: energy={} Ha; residual={} Bohr^-3/2; solver=iteration_completed",
            v["iteration"], d["energy_hartree"], d["density_residual_bohr_minus_3_over_2"]
        )
    } else {
        format!(
            "iteration {}: {} {}",
            v["iteration"],
            v["type"].as_str().unwrap_or("unknown"),
            d
        )
    }
}

/// Reference map payload in the existing event envelope. Missing energy and
/// legacy residual diagnostics remain null; dr2 is an energy metric in Ry,
/// never relabelled as the legacy real-space density residual.
pub fn reference_step_json(
    calc_id: &str,
    sequence: u64,
    state_version: u64,
    timestamp: std::time::SystemTime,
    report: dft_engine::refscf::stepper::StepReport,
) -> Value {
    use dft_engine::refscf::stepper::StepKind;
    let kind = match report.kind {
        StepKind::Ordinary => "ordinary",
        StepKind::Verification => "verification",
        StepKind::Retry => "retry",
    };
    json!({
        "schema_version":"1.0", "calc_id":calc_id, "sequence":sequence,
        "timestamp_unix_seconds":timestamp.duration_since(UNIX_EPOCH).map(|d|d.as_secs_f64()).unwrap_or(0.0),
        "iteration":report.iteration, "state_version":state_version,
        "type":"iteration_finished",
        "data":{
            "driver":"reference", "kind":kind, "stage":report.stage,
            "energy_hartree":report.energy_ha, "estimator_hartree":report.estimator_ha,
            "dr2_ry":report.dr2_ry, "ethr_ry":report.ethr_ry,
            "h_columns":report.h_columns, "reused":report.reused,
            "safe_point":report.safe_point,
            "density_residual_bohr_minus_3_over_2":null,
            "charge_error_relative":null, "mixing_applied":null,
            "solver_status":"map_completed"
        }
    })
}

/// Publication is bounded and nonblocking. Capacity exhaustion fails at a safe
/// boundary; queued records are drained before returning the error.
struct Sink {
    sender: SyncSender<Value>,
    worker: JoinHandle<io::Result<()>>,
}
impl Sink {
    fn new(out: &Path, mirror: Option<&str>) -> Result<Self, RunError> {
        // Never overwrite a previous log or an input via a symlink.
        let mut journal = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(out.join("events.jsonl"))?;
        let mut output: Option<Box<dyn Write + Send>> = match mirror {
            None => None,
            Some("-") => Some(Box::new(io::stdout())),
            Some(s) if s.starts_with("unix:") => None,
            Some(path) => Some(Box::new(
                OpenOptions::new().write(true).create_new(true).open(path)?,
            )),
        };
        #[cfg(unix)]
        let socket = if let Some(path) = mirror.and_then(|s| s.strip_prefix("unix:")) {
            let socket = std::os::unix::net::UnixDatagram::unbound()?;
            socket.connect(path)?;
            socket.set_nonblocking(true)?;
            Some(socket)
        } else {
            None
        };
        #[cfg(not(unix))]
        if mirror.is_some_and(|s| s.starts_with("unix:")) {
            return Err(fail("events: unix sockets require Unix"));
        }
        let (sender, receiver) = mpsc::sync_channel::<Value>(256);
        let worker = thread::spawn(move || {
            for value in receiver {
                let mut line = serde_json::to_vec(&value)?;
                line.push(b'\n');
                journal.write_all(&line)?;
                journal.flush()?;
                writeln!(io::stderr(), "{}", human(&value))?;
                if let Some(w) = &mut output {
                    w.write_all(&line)?;
                    w.flush()?;
                }
                #[cfg(unix)]
                if let Some(socket) = &socket {
                    socket.send(&line)?;
                }
            }
            journal.sync_data()
        });
        Ok(Self { sender, worker })
    }
    fn drain(&self, events: Option<&mut EventQueue>) -> Result<(), RunError> {
        if let Some(events) = events {
            while let Some(event) = events.pop() {
                self.sender.try_send(event_json(&event)).map_err(|e| {
                    fail(format!(
                        "telemetry: bounded sink unavailable ({e}); stopped at iteration boundary"
                    ))
                })?;
            }
        }
        Ok(())
    }
    fn finish(self) -> Result<(), RunError> {
        drop(self.sender);
        self.worker
            .join()
            .map_err(|_| fail("telemetry worker panicked"))?
            .map_err(fail)
    }
}

pub fn drive(problem: ScfProblem, mode: &RunMode, out: &Path) -> Result<ScfOutcome, RunError> {
    let config = ExecutionConfig {
        mode: if mode.agent {
            ExecutionMode::AgentControlled
        } else {
            ExecutionMode::Observed
        },
        // One log per output directory; each new calculation requires a fresh log.
        calc_id: format!(
            "cli-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .map_err(fail)?
                .as_nanos()
        ),
        ..ExecutionConfig::default()
    };
    let session =
        ScfSession::new(problem, CompletionMode::RetainLast, config).map_err(engine_error)?;
    let sink = Sink::new(out, mode.events.as_deref())?;
    let result = if mode.agent {
        drive_agent(session, &sink, mode, out)
    } else {
        drive_observed(session, &sink)
    };
    let flushed = sink.finish();
    flushed?;
    result
}
fn drive_observed(mut session: ScfSession, sink: &Sink) -> Result<ScfOutcome, RunError> {
    let initialized = session.initialize();
    sink.drain(session.events())?;
    initialized.map_err(engine_error)?;
    loop {
        let result = session.step();
        sink.drain(session.events())?;
        if result.map_err(engine_error)?.status != StepStatus::Ready {
            break;
        }
    }
    session
        .into_result()
        .ok_or_else(|| fail("session: no final result"))
}

#[cfg(unix)]
#[path = "control_channel.rs"]
mod control_channel;
#[cfg(unix)]
fn drive_agent(
    session: ScfSession,
    sink: &Sink,
    mode: &RunMode,
    out: &Path,
) -> Result<ScfOutcome, RunError> {
    let path = if let Some(path) = &mode.control_socket {
        path.clone()
    } else {
        use std::os::unix::fs::DirBuilderExt;
        let directory = out.join("control");
        match std::fs::DirBuilder::new().mode(0o700).create(&directory) {
            Ok(()) => {}
            Err(e) if e.kind() == io::ErrorKind::AlreadyExists => {}
            Err(e) => return Err(fail(e)),
        }
        directory.join("control.sock")
    };
    let mut channel = control_channel::Channel::new(&path)?;
    let audit = OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(out.join("control-audit.jsonl"))?;
    let mut runtime = Runtime::new(session, RuntimeConfig::default())
        .map_err(|e| fail(format!("{e:?}")))?
        .with_audit_sink(Box::new(control_channel::Journal(audit)))
        .map_err(|e| fail(format!("{e:?}")))?;
    let initialized = runtime.initialize();
    sink.drain(runtime.events())?;
    initialized.map_err(|e| fail(format!("{e:?}")))?;
    loop {
        channel.poll(&mut runtime)?;
        runtime.safe_point().map_err(|e| fail(format!("{e:?}")))?;
        if runtime.state().terminal() {
            break;
        }
        if runtime.state() == RuntimeState::Paused {
            thread::sleep(std::time::Duration::from_millis(20));
            continue;
        }
        let result = runtime.step();
        sink.drain(runtime.events())?;
        result.map_err(|e| fail(format!("{e:?}")))?;
        if runtime.state().terminal() {
            break;
        }
    }
    let result = runtime.into_result();
    result.outcome.ok_or_else(|| {
        fail(format!(
            "control: terminal {:?}; no numerical result",
            result.state
        ))
    })
}
#[cfg(not(unix))]
fn drive_agent(_: ScfSession, _: &Sink, _: &RunMode, _: &Path) -> Result<ScfOutcome, RunError> {
    Err(fail("agent: local control channel currently requires Unix"))
}

/// The reference facade supplies all numerical work. I/O and supervision occur
/// only between its returns, and only rank zero accepts external controls.
pub fn drive_reference(
    input: dft_engine::refstart::StartInput<'_>,
    request: dft_engine::refscf::api::Request,
    mode: &RunMode,
    out: &Path,
) -> Result<dft_engine::refscf::api::Delivered, RunError> {
    use dft_engine::refscf::runtime::Runtime as ReferenceRuntime;
    use dft_mpi::Collective;
    let comm = dft_engine::parallel::SharedComm::active();
    let comm = comm.as_ref().map(|c| c as &dyn Collective);
    let root = comm.is_none_or(|c| c.rank() == 0);
    let identity = if root {
        json!(format!(
            "reference-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .map_err(fail)?
                .as_nanos()
        ))
    } else {
        Value::Null
    };
    let identity = reference_broadcast(Ok(identity), comm)?;
    let calc_id = identity.as_str().ok_or_else(|| fail("runtime identity"))?;
    let required = request.derivatives.forces || request.derivatives.stress;
    let threshold_ry = request.profile.threshold_ry;
    // Only root mirrors events to a caller-specified external endpoint.
    let sink = reference_agree(
        Sink::new(out, if root { mode.events.as_deref() } else { None }),
        comm,
    )?;
    let mut runtime = ReferenceRuntime::new(dft_engine::refscf::api::stepper(input, request));
    #[cfg(unix)]
    let mut channel = {
        let create = || -> Result<Option<control_channel::ReferenceChannel>, RunError> {
            if !mode.agent || !root {
                return Ok(None);
            }
            let path = if let Some(path) = &mode.control_socket {
                path.clone()
            } else {
                use std::os::unix::fs::DirBuilderExt;
                let dir = out.join("control");
                std::fs::DirBuilder::new().mode(0o700).create(&dir)?;
                dir.join("control.sock")
            };
            Ok(Some(control_channel::ReferenceChannel::new(&path, out)?))
        };
        reference_agree(create(), comm)?
    };
    #[cfg(not(unix))]
    if mode.agent {
        return Err(fail("reference agent control requires Unix"));
    }
    let mut sequence = 0;
    let publish = |value: Value| -> Result<(), RunError> {
        reference_agree(
            sink.sender
                .try_send(value)
                .map_err(|e| fail(format!("reference telemetry: {e}"))),
            comm,
        )
    };
    let run = (|| -> Result<(), RunError> {
        crate::progress::scf_start(threshold_ry);
        reference_agree(runtime.initialize().map_err(engine_error), comm)?;
        crate::progress::phase("scf");
        sequence += 1;
        publish(reference_state_json(
            calc_id,
            sequence,
            runtime.snapshot(),
            "initialized",
        ))?;
        while !runtime.snapshot().state.terminal() {
            #[cfg(unix)]
            if mode.agent && runtime.snapshot().safe_point {
                let proposals = if let Some(channel) = &mut channel {
                    channel
                        .poll(calc_id, runtime.snapshot())
                        .and_then(|p| serde_json::to_value(p).map_err(fail))
                } else {
                    Ok(json!([]))
                };
                let proposals: Vec<control_channel::Proposal> =
                    serde_json::from_value(reference_broadcast(proposals, comm)?).map_err(fail)?;
                for proposal in proposals {
                    let applied = control_channel::apply_reference(
                        &proposal,
                        calc_id,
                        &mut runtime,
                        required,
                    );
                    let journal = if let Some(channel) = &mut channel {
                        channel.disposition(&proposal, &applied)
                    } else {
                        Ok(())
                    };
                    reference_agree(journal, comm)?;
                    if let Some(comm) = comm {
                        let mut applied_count = [0.];
                        comm.tree_sum_f64(&[f64::from(applied.is_ok())], &mut applied_count)
                            .map_err(fail)?;
                        if applied_count[0] != 0. && applied_count[0] != comm.size() as f64 {
                            return Err(fail(
                                "reference control application disagreed across ranks; stopped before next map",
                            ));
                        }
                    }
                    sequence += 1;
                    let mut event = reference_state_json(
                        calc_id,
                        sequence,
                        runtime.snapshot(),
                        "control_disposition",
                    );
                    event["data"]["request"] = serde_json::to_value(&proposal).map_err(fail)?;
                    event["data"]["result"] = json!(applied);
                    publish(event)?;
                }
            }
            if runtime.snapshot().state.terminal() {
                break;
            }
            if runtime.snapshot().state == RuntimeState::Paused {
                thread::sleep(std::time::Duration::from_millis(20));
                continue;
            }
            let report = reference_agree(runtime.step().map_err(engine_error), comm)?;
            crate::progress::scf_map(&report);
            sequence += 1;
            publish(reference_step_json(
                calc_id,
                sequence,
                runtime.snapshot().state_version,
                std::time::SystemTime::now(),
                report,
            ))?;
        }
        sequence += 1;
        let kind = match runtime.snapshot().state {
            RuntimeState::Converged => "converged",
            RuntimeState::Exhausted => "not_converged",
            RuntimeState::Cancelled => "cancelled",
            _ => "failed",
        };
        crate::progress::scf_end(Some(kind == "converged"), None);
        crate::progress::phase("post-scf");
        publish(reference_state_json(
            calc_id,
            sequence,
            runtime.snapshot(),
            kind,
        ))?;
        Ok(())
    })();
    let failed_event = if let Err(error) = &run {
        sequence += 1;
        let mut event = reference_state_json(calc_id, sequence, runtime.snapshot(), "failed");
        event["data"]["error"] = json!(error.cause);
        publish(event)
    } else {
        Ok(())
    };
    let snapshot = runtime.snapshot();
    let summary = json!({"schema_version":1,"driver":"reference","calc_id":calc_id,
        "mode":if mode.agent {"agent"} else {"observed"},
        "state":format!("{:?}",snapshot.state),"state_version":snapshot.state_version,
        "completed_maps":snapshot.completed_maps,"safe_point":snapshot.safe_point,
        "policy_revision":snapshot.policy_revision,"delivered":runtime.result().is_some(),
        "error":run.as_ref().err().map(|e|&e.cause),
        "last_step":snapshot.last_step.map(|r|reference_step_json(calc_id, sequence, snapshot.state_version, UNIX_EPOCH,r)),
        "policies":runtime.policy_records().iter().map(|r|json!({"revision":r.revision,"completed_maps":r.completed_maps,
            "iteration":r.iteration,"state_version":r.state_version,"history_reset":r.history_reset,
            "beta":r.patch.beta,"history":r.patch.history,"threshold_ry":r.patch.threshold_ry,"verification":r.patch.verification})).collect::<Vec<_>>()});
    let saved = (|| -> Result<(), RunError> {
        let mut file = OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(out.join("runtime.json"))?;
        serde_json::to_writer_pretty(&mut file, &summary).map_err(fail)?;
        file.write_all(b"\n")?;
        file.sync_data()?;
        Ok(())
    })();
    #[cfg(unix)]
    let closed = if let Some(channel) = &mut channel {
        channel.close()
    } else {
        Ok(())
    };
    #[cfg(not(unix))]
    let closed = Ok(());
    // Drain durable event output even when the driver or control channel fails.
    let flushed = sink.finish();
    reference_agree(saved.and(closed).and(flushed).and(failed_event), comm)?;
    run?;
    if snapshot.state == RuntimeState::Cancelled {
        return Err(fail(
            "reference control: cancelled at a safe point; partial record is in runtime.json",
        ));
    }
    runtime.finish().map_err(engine_error)
}
fn reference_state_json(
    calc_id: &str,
    sequence: u64,
    snapshot: dft_engine::refscf::runtime::Snapshot,
    kind: &str,
) -> Value {
    json!({"schema_version":"1.0","calc_id":calc_id,"sequence":sequence,
        "timestamp_unix_seconds":std::time::SystemTime::now().duration_since(UNIX_EPOCH).map(|d|d.as_secs_f64()).unwrap_or(0.),
        "iteration":snapshot.last_step.map_or(0, |s|s.iteration),"state_version":snapshot.state_version,"type":kind,
        "data":{"driver":"reference","state":format!("{:?}",snapshot.state),"policy_revision":snapshot.policy_revision,
            "safe_point":snapshot.safe_point,"completed_maps":snapshot.completed_maps}})
}
fn reference_agree<T>(
    result: Result<T, RunError>,
    comm: Option<&dyn dft_mpi::Collective>,
) -> Result<T, RunError> {
    if let Some(comm) = comm {
        let mut failures = [0.];
        comm.tree_sum_f64(&[f64::from(result.is_err())], &mut failures)
            .map_err(fail)?;
        if failures[0] > 0. {
            return Err(result
                .err()
                .unwrap_or_else(|| fail("reference runtime: another rank failed")));
        }
    }
    result
}
fn reference_broadcast(
    value: Result<Value, RunError>,
    comm: Option<&dyn dft_mpi::Collective>,
) -> Result<Value, RunError> {
    let Some(comm) = comm else {
        return value;
    };
    let bytes = if comm.rank() == 0 {
        serde_json::to_vec(&value.map_err(|e| e.cause)).map_err(fail)?
    } else {
        Vec::new()
    };
    let mut size = [bytes.len() as f64];
    comm.bcast_f64(&mut size, 0).map_err(fail)?;
    if size[0] > 131072. {
        return Err(fail("reference control broadcast exceeds bound"));
    }
    let mut payload = vec![0.; size[0] as usize];
    for (v, b) in payload.iter_mut().zip(bytes) {
        *v = f64::from(b);
    }
    comm.bcast_f64(&mut payload, 0).map_err(fail)?;
    let bytes: Vec<_> = payload.into_iter().map(|v| v as u8).collect();
    serde_json::from_slice::<Result<Value, String>>(&bytes)
        .map_err(fail)?
        .map_err(fail)
}

/// Attach the audit and report the effective settings after live policy changes.
/// Preserve launch settings explicitly instead of labelling them as effective.
pub fn annotate_reference_result(result: &mut Value, out: &Path) -> Result<(), RunError> {
    let mut runtime: Value =
        serde_json::from_slice(&std::fs::read(out.join("runtime.json"))?).map_err(fail)?;
    runtime["launch_policy"] = json!({"mixing":result["scf"]["mixing"],
        "threshold_ha":result["single_path"]["standard_threshold_ha"],
        "verification":result["single_path"]["verification"]});
    if let Some(policies) = runtime["policies"].as_array() {
        for policy in policies {
            for key in ["beta", "history"] {
                if !policy[key].is_null() {
                    result["scf"]["mixing"][key] = policy[key].clone();
                }
            }
            if let Some(ry) = policy["threshold_ry"].as_f64() {
                result["reference_scf"]["threshold_ry"] = json!(ry);
                result["reference_scf"]["threshold_ha"] = json!(ry / 2.);
                result["standard_threshold_ha"] = json!(ry / 2.);
                for key in ["scf", "single_path", "execution"] {
                    result[key]["standard_threshold_ha"] = json!(ry / 2.);
                }
                result["single_path"]["threshold_source"] = json!("runtime_policy");
            }
            if let Some(enabled) = policy["verification"].as_bool() {
                let name = if enabled {
                    dft_engine::Verification::FixedPoint
                } else {
                    dft_engine::Verification::Off
                }
                .name();
                for key in ["single_path", "execution"] {
                    result[key]["verification"] = json!(name);
                }
            }
        }
    }
    result["execution"]["runtime_mode"] = runtime["mode"].clone();
    result["runtime"] = runtime;
    Ok(())
}
#[cfg(test)]
mod reference_collective_tests {
    use super::*;
    use dft_mpi::Collective;
    #[test]
    fn root_controls_and_rank_errors_are_collective() {
        dft_mpi::simulate(3, |comm| {
            let value = if comm.rank() == 0 {
                json!([{"action":"pause","id":"p"}])
            } else {
                Value::Null
            };
            let received = reference_broadcast(Ok(value), Some(comm)).unwrap();
            assert_eq!(received, json!([{"action":"pause","id":"p"}]));
            let result = if comm.rank() == 0 {
                Err(fail("socket failed"))
            } else {
                Ok(Value::Null)
            };
            assert!(
                reference_broadcast(result, Some(comm))
                    .unwrap_err()
                    .cause
                    .contains("socket failed")
            );
            let result = if comm.rank() == 2 {
                Err(fail("sink failed"))
            } else {
                Ok(())
            };
            assert!(reference_agree(result, Some(comm)).is_err());
        });
    }
}
