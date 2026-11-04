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
    if v["type"] == "iteration_finished" {
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
