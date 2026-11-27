//! UX1 progress stream and job lock (ADR 0053, schema `baandi-progress-1`).
//!
//! Every job started through the CLI (bare `baandi` or an alias subcommand)
//! writes `progress.jsonl` and holds `baandi.lock` in its output directory;
//! `baandi watch` (`crate::monitor`) tails the file. Records are written only
//! at job, phase, SCF-map and relaxation-step boundaries, never inside a map:
//! one small JSON line per boundary, written and flushed directly. Only MPI
//! rank zero writes. The stream is advisory: a failed write disables it with
//! one warning and never changes or stops the calculation, and nothing here
//! feeds back into the numerics (the per-map driver loop is the same sequence
//! of calls as `refscf::api::solve`).
use crate::config::Config;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::{
    cell::RefCell,
    fs::{self, File, OpenOptions},
    io::Write,
    path::{Path, PathBuf},
    sync::OnceLock,
    time::{Instant, SystemTime, UNIX_EPOCH},
};

pub const FILE: &str = "progress.jsonl";
/// The stream of the previous job in the same directory (one generation).
pub const PREVIOUS: &str = "progress.prev.jsonl";
pub const LOCK: &str = "baandi.lock";
pub const SCHEMA: &str = "baandi-progress-1";
pub const LOCK_SCHEMA: &str = "baandi-lock-1";

fn unix_now() -> f64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0., |d| d.as_secs_f64())
}

/// Host name used for lock liveness (a lock from another host is never
/// declared stale by this host).
pub fn host() -> &'static str {
    static HOST: OnceLock<String> = OnceLock::new();
    HOST.get_or_init(|| {
        std::process::Command::new("hostname")
            .output()
            .ok()
            .and_then(|o| String::from_utf8(o.stdout).ok())
            .map(|s| s.trim().to_owned())
            .filter(|s| !s.is_empty())
            .unwrap_or_else(|| "unknown".into())
    })
}

/// Contents of `baandi.lock`.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct LockInfo {
    pub schema: String,
    pub pid: u32,
    pub host: String,
    pub started_unix: f64,
    pub task: String,
    pub input: String,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Liveness {
    /// The process is running on this host.
    Alive,
    /// This host has no such process: the lock is stale.
    Dead,
    /// The lock is from another host; it cannot be checked from here.
    Remote(String),
}

/// Whether a process id runs on this host (`ps -p`, no signals, no FFI). A
/// zombie (exited, not yet reaped) is not running.
pub fn pid_alive(pid: u32) -> bool {
    if pid == std::process::id() {
        return true;
    }
    std::process::Command::new("ps")
        .args(["-p", &pid.to_string(), "-o", "stat="])
        .output()
        .is_ok_and(|o| {
            let stat = o.stdout.trim_ascii();
            o.status.success() && !stat.is_empty() && !stat.starts_with(b"Z")
        })
}

pub fn liveness(lock: &LockInfo) -> Liveness {
    if lock.host != host() {
        Liveness::Remote(lock.host.clone())
    } else if pid_alive(lock.pid) {
        Liveness::Alive
    } else {
        Liveness::Dead
    }
}

pub fn read_lock(dir: &Path) -> Option<LockInfo> {
    serde_json::from_slice(&fs::read(dir.join(LOCK)).ok()?).ok()
}

struct Writer {
    file: File,
    sequence: u64,
    start: Instant,
    cycle: u64,
    maps: u64,
    last_energy: Option<f64>,
    map_clock: Instant,
    relax_frames: usize,
    relax_failures: usize,
}

// Per thread: the job's driver thread owns the stream (the reference SCF and
// the relax drivers run on the caller's thread; worker threads never emit),
// so in-process callers and parallel tests cannot write into another job.
thread_local! {
    static SINK: RefCell<Option<Writer>> = const { RefCell::new(None) };
}

fn with_sink<R>(f: impl FnOnce(&mut Writer) -> R) -> Option<R> {
    SINK.with_borrow_mut(|s| s.as_mut().map(f))
}

/// True when this thread writes a progress stream (rank zero of a CLI job).
pub fn active() -> bool {
    SINK.with_borrow(Option::is_some)
}

/// Append one record. The envelope (`v`, `seq`, `t`, `elapsed_s`, `event`) is
/// added here; `fields` must be a JSON object.
pub fn emit(event: &str, fields: Value) {
    SINK.with_borrow_mut(|guard| {
        let Some(w) = guard.as_mut() else {
            return;
        };
        w.sequence += 1;
        let mut record = json!({"v":1,"seq":w.sequence,"t":unix_now(),
        "elapsed_s":w.start.elapsed().as_secs_f64(),"event":event});
        if let (Some(r), Value::Object(f)) = (record.as_object_mut(), fields) {
            r.extend(f);
        }
        let mut line = serde_json::to_vec(&record).unwrap_or_default();
        line.push(b'\n');
        if let Err(e) = w.file.write_all(&line).and_then(|()| w.file.flush()) {
            eprintln!("warning: progress stream disabled after a write error: {e}");
            *guard = None;
        }
    });
}

pub fn phase(name: &str) {
    emit("phase", json!({"phase":name}));
}
pub fn warning(source: &str, message: &str) {
    emit("warning", json!({"source":source,"message":message}));
}

/// Start of one SCF (a relaxation runs one per ionic step).
pub fn scf_start(threshold_ry: f64) {
    let cycle = with_sink(|w| {
        w.cycle += 1;
        w.maps = 0;
        w.last_energy = None;
        w.map_clock = Instant::now();
        w.cycle
    });
    if let Some(cycle) = cycle {
        emit(
            "scf_start",
            json!({"cycle":cycle,"threshold_ry":threshold_ry,"threshold_ha":threshold_ry/2.}),
        );
    }
}

/// One completed SCF map (ordinary, verification or retry).
pub fn scf_map(report: &dft_engine::refscf::stepper::StepReport) {
    use dft_engine::refscf::stepper::StepKind;
    let Some((cycle, seconds, de)) = with_sink(|w| {
        let seconds = w.map_clock.elapsed().as_secs_f64();
        w.map_clock = Instant::now();
        w.maps += 1;
        let de = match (report.energy_ha, w.last_energy) {
            (Some(e), Some(p)) => Some(e - p),
            _ => None,
        };
        if report.energy_ha.is_some() {
            w.last_energy = report.energy_ha;
        }
        (w.cycle, seconds, de)
    }) else {
        return;
    };
    let kind = match report.kind {
        StepKind::Ordinary => "ordinary",
        StepKind::Verification => "verification",
        StepKind::Retry => "retry",
    };
    emit(
        "scf_map",
        json!({"cycle":cycle,"iteration":report.iteration,"kind":kind,"stage":report.stage,
            "energy_ha":report.energy_ha,"de_ha":de,"dr2_ry":report.dr2_ry,
            "estimator_ha":report.estimator_ha,"ethr_ry":report.ethr_ry,
            "h_columns":report.h_columns,"reused":report.reused,"seconds":seconds}),
    );
    if let Some(stall) = &report.stall {
        warning(
            "scf",
            &format!(
                "stall diagnostic at map {}: no decrease of the stop estimate for three maps with an ill-conditioned history",
                stall.iteration
            ),
        );
    }
}

pub fn scf_end(converged: Option<bool>, error: Option<&str>) {
    let Some((cycle, maps)) = with_sink(|w| (w.cycle, w.maps)) else {
        return;
    };
    emit(
        "scf_end",
        json!({"cycle":cycle,"maps":maps,"converged":converged,"error":error}),
    );
}

/// The reference SCF of `refscf::api::solve`, driven map by map when a
/// progress stream is active. `api::solve` is `stepper(..).run_to_end()`, and
/// this loop makes the same calls in the same order (initialize, step until
/// done, finish), so the numerical result is identical; only boundary
/// records are added.
pub fn solve(
    input: dft_engine::refstart::StartInput<'_>,
    request: dft_engine::refscf::api::Request,
) -> Result<dft_engine::refscf::api::Delivered, dft_engine::ScfError> {
    use dft_engine::refscf::api;
    if !active() {
        return api::solve(input, request);
    }
    scf_start(request.profile.threshold_ry);
    let result = (|| {
        let mut stepper = api::stepper(input, request);
        stepper.initialize()?;
        phase("scf");
        while !stepper.is_done() {
            let report = stepper.step()?;
            scf_map(&report);
        }
        phase("post-scf");
        stepper.finish()
    })();
    match &result {
        Ok(delivered) => {
            let converged = match &delivered.run {
                api::Run::Scalar(r) => Some(r.outcome.converged),
                _ => None,
            };
            scf_end(converged, None);
        }
        Err(e) => scf_end(None, Some(&e.to_string())),
    }
    result
}

/// One relaxation frame written by the relax drivers' output step. Repeated
/// calls with the same frame and failure counts are ignored.
#[allow(clippy::too_many_arguments)]
pub fn relax_frame(
    frames: usize,
    energy_ha: f64,
    max_force_ha_bohr: f64,
    step_length_bohr: f64,
    scf_iterations: usize,
    provenance: &str,
    stress_ha_bohr3: Option<[[f64; 3]; 3]>,
    status: &str,
    failures: &[Value],
) {
    let Some((new_frame, new_failures)) = with_sink(|w| {
        let r = (frames != w.relax_frames, failures.len() > w.relax_failures);
        let first_failure = w.relax_failures;
        w.relax_frames = frames;
        w.relax_failures = failures.len();
        (r.0, r.1.then_some(first_failure))
    }) else {
        return;
    };
    if let Some(first) = new_failures {
        for f in &failures[first..] {
            let reason = f["reason"].as_str().unwrap_or("unknown");
            warning(
                "relax",
                &format!(
                    "step {} attempt {} failed: {reason}",
                    f["step"], f["attempt"]
                ),
            );
        }
    }
    if new_frame {
        let pressure = stress_ha_bohr3.map(|t| -(t[0][0] + t[1][1] + t[2][2]) / 3.);
        emit(
            "relax_step",
            json!({"step":frames.saturating_sub(1),"energy_ha":energy_ha,
                "max_force_ha_bohr":max_force_ha_bohr,"step_length_bohr":step_length_bohr,
                "scf_iterations":scf_iterations,"provenance":provenance,
                "pressure_ha_bohr3":pressure,"status":status}),
        );
    }
}

/// System summary of the `start` record (read from the parsed input only).
pub fn summary(config: &Config) -> Value {
    let mut species: Vec<(String, usize)> = Vec::new();
    for a in &config.system.atoms {
        match species.iter_mut().find(|(s, _)| *s == a.species) {
            Some((_, n)) => *n += 1,
            None => species.push((a.species.clone(), 1)),
        }
    }
    let formula: String = species
        .iter()
        .map(|(s, n)| {
            if *n == 1 {
                s.clone()
            } else {
                format!("{s}{n}")
            }
        })
        .collect();
    let scf = config.scf.as_ref();
    let relax = config.relax.as_ref().map(|r| {
        json!({"mode":format!("{:?}",r.mode).to_lowercase(),"force_tol_ha_bohr":r.force_tol_ha_bohr,
            "energy_tol_ha":r.energy_tol_ha,"stress_tol_ha_bohr3":r.stress_tol_ha_bohr3,"max_steps":r.max_steps})
    });
    json!({
        "name":config.system.name,
        "formula":formula,
        "atoms":config.system.atoms.len(),
        "kpoints":config.resolved_points().ok().map(|p|p.len()),
        "nbands":scf.and_then(|s|s.nbands),
        "ecutwfc_ry":config.basis.ecutwfc_ry,
        "volume_bohr3":config.cell().ok().map(|c|c.volume_bohr3()),
        "spin":config.collinear(),
        "symmetry":format!("{:?}",config.symmetry.mode).to_lowercase(),
        "threads":scf.and_then(|s|s.threads),
        "ranks":crate::parallel::size(),
        "relax":relax,
    })
}

/// A started job: holds the lock and the installed writer until `finish` or
/// drop (a drop without `finish` records `aborted`).
pub struct Session {
    dir: Option<PathBuf>,
    finished: bool,
    /// The session created the directory: a job that wrote nothing else there
    /// (rejected before any output) leaves no directory behind, as before.
    created: bool,
}

impl Session {
    /// A session that writes nothing (MPI ranks above zero, `--no-progress`).
    pub fn inert() -> Self {
        Self {
            dir: None,
            finished: true,
            created: false,
        }
    }

    /// Take the lock of `dir` (replacing a stale lock), move a previous
    /// stream to `progress.prev.jsonl`, open a fresh stream and write `start`.
    /// A live lock (same host, process running) or a lock from another host
    /// refuses the job.
    pub fn begin(dir: &Path, task: &str, input: &Path, system: Value) -> Result<Self, String> {
        if crate::parallel::rank() > 0 {
            return Ok(Self::inert());
        }
        let ctx = |e: std::io::Error| format!("progress: {}: {e}", dir.display());
        let created = !dir.exists();
        fs::create_dir_all(dir).map_err(ctx)?;
        let info = LockInfo {
            schema: LOCK_SCHEMA.into(),
            pid: std::process::id(),
            host: host().into(),
            started_unix: unix_now(),
            task: task.into(),
            input: input.display().to_string(),
        };
        let mut stale = None;
        let lock_path = dir.join(LOCK);
        for _ in 0..2 {
            match OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&lock_path)
            {
                Ok(mut f) => {
                    let bytes = serde_json::to_vec(&info).map_err(|e| e.to_string())?;
                    f.write_all(&bytes).map_err(ctx)?;
                    f.flush().map_err(ctx)?;
                    stale = stale.or(Some(None));
                    break;
                }
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {
                    let previous = read_lock(dir);
                    match previous.as_ref().map(liveness) {
                        Some(Liveness::Alive) => {
                            let p = previous.unwrap();
                            return Err(format!(
                                "a BAANDI job (pid {}, task {}) is running in {}; watch it with `baandi watch {}`",
                                p.pid,
                                p.task,
                                dir.display(),
                                dir.display()
                            ));
                        }
                        Some(Liveness::Remote(h)) => {
                            return Err(format!(
                                "{} is locked by a job on host {h}; remove {} only if that job has ended",
                                dir.display(),
                                lock_path.display()
                            ));
                        }
                        Some(Liveness::Dead) | None => {
                            fs::remove_file(&lock_path).map_err(ctx)?;
                            stale = Some(Some(previous.map_or_else(
                                || "unreadable lock replaced".to_owned(),
                                |p| format!("stale lock of pid {} ({}) replaced", p.pid, p.task),
                            )));
                        }
                    }
                }
                Err(e) => return Err(ctx(e)),
            }
        }
        let Some(stale) = stale else {
            return Err(format!("progress: could not take {}", lock_path.display()));
        };
        let session = Self {
            dir: Some(dir.to_owned()),
            finished: false,
            created,
        };
        let path = dir.join(FILE);
        if let Ok(meta) = fs::symlink_metadata(&path) {
            if !meta.is_file() {
                return Err(format!("progress: refusing nonregular {}", path.display()));
            }
            fs::rename(&path, dir.join(PREVIOUS)).map_err(ctx)?;
        }
        let file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .map_err(ctx)?;
        SINK.set(Some(Writer {
            file,
            sequence: 0,
            start: Instant::now(),
            cycle: 0,
            maps: 0,
            last_energy: None,
            map_clock: Instant::now(),
            relax_frames: 0,
            relax_failures: 0,
        }));
        emit(
            "start",
            json!({"schema":SCHEMA,"task":task,"input":input.display().to_string(),
                "output_dir":dir.display().to_string(),"pid":info.pid,"host":info.host,
                "version":env!("CARGO_PKG_VERSION"),"system":system}),
        );
        if let Some(message) = stale {
            warning("lock", &message);
        }
        phase("setup");
        Ok(session)
    }

    /// Record the final status (`completed`, `not_converged`, `failed`) and
    /// release the lock.
    pub fn finish(mut self, status: &str, exit_code: u8, message: Option<&str>) {
        if self.dir.is_some() {
            emit(
                "end",
                json!({"status":status,"exit_code":exit_code,"message":message}),
            );
        }
        self.finished = true;
    }
}

impl Drop for Session {
    fn drop(&mut self) {
        let Some(dir) = self.dir.take() else {
            return;
        };
        if !self.finished {
            emit(
                "end",
                json!({"status":"aborted","exit_code":Value::Null,"message":"job ended without a final status"}),
            );
        }
        if let Some(w) = SINK.take() {
            let _ = w.file.sync_data();
        }
        let _ = fs::remove_file(dir.join(LOCK));
        let only_stream = fs::read_dir(&dir).is_ok_and(|entries| {
            entries
                .flatten()
                .all(|e| e.file_name() == FILE || e.file_name() == PREVIOUS)
        });
        if self.created && only_stream {
            let _ = fs::remove_file(dir.join(FILE));
            let _ = fs::remove_dir(&dir);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scratch(name: &str) -> PathBuf {
        let dir =
            std::env::temp_dir().join(format!("baandi-progress-{name}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        dir
    }
    fn records(dir: &Path) -> Vec<Value> {
        fs::read_to_string(dir.join(FILE))
            .unwrap()
            .lines()
            .map(|l| serde_json::from_str(l).unwrap())
            .collect()
    }

    #[test]
    fn session_writes_envelope_and_releases_the_lock() {
        let dir = scratch("session");
        fs::create_dir_all(&dir).unwrap();
        let s =
            Session::begin(&dir, "scf", Path::new("in.toml"), json!({"formula":"Si2"})).unwrap();
        assert!(active());
        let lock = read_lock(&dir).unwrap();
        assert_eq!(lock.pid, std::process::id());
        assert_eq!(liveness(&lock), Liveness::Alive);
        scf_start(2e-8);
        relax_frame(1, -7.9, 0.01, 0., 5, "initial", None, "running", &[]);
        relax_frame(1, -7.9, 0.01, 0., 5, "initial", None, "running", &[]);
        s.finish("completed", 0, None);
        assert!(!active());
        assert!(!dir.join(LOCK).exists());
        let r = records(&dir);
        let events: Vec<_> = r.iter().map(|v| v["event"].as_str().unwrap()).collect();
        assert_eq!(
            events,
            ["start", "phase", "scf_start", "relax_step", "end"],
            "duplicate frames are not repeated"
        );
        for (i, v) in r.iter().enumerate() {
            assert_eq!(v["v"], 1);
            assert_eq!(v["seq"], i as u64 + 1);
        }
        assert_eq!(r[0]["schema"], SCHEMA);
        assert_eq!(r[2]["threshold_ha"], 1e-8);
        // A second job in the same directory keeps one previous generation.
        let s = Session::begin(&dir, "scf", Path::new("in.toml"), Value::Null).unwrap();
        drop(s);
        assert_eq!(records(&dir).last().unwrap()["status"], "aborted");
        assert!(dir.join(PREVIOUS).is_file());
        fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn a_job_rejected_before_any_output_leaves_no_directory() {
        let dir = scratch("rejected");
        let s = Session::begin(&dir, "scf", Path::new("in.toml"), Value::Null).unwrap();
        s.finish("failed", 1, Some("input error"));
        assert!(!dir.exists());
        // An existing directory, or one with other outputs, is kept.
        let s = Session::begin(&dir, "scf", Path::new("in.toml"), Value::Null).unwrap();
        fs::write(dir.join("results.json"), "{}").unwrap();
        s.finish("failed", 4, None);
        assert!(dir.join(FILE).is_file() && !dir.join(LOCK).exists());
        fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn live_locks_refuse_and_stale_locks_are_replaced() {
        let dir = scratch("lock");
        fs::create_dir_all(&dir).unwrap();
        let mut info = LockInfo {
            schema: LOCK_SCHEMA.into(),
            pid: std::process::id(),
            host: host().into(),
            started_unix: 0.,
            task: "relax".into(),
            input: "x.toml".into(),
        };
        fs::write(dir.join(LOCK), serde_json::to_vec(&info).unwrap()).unwrap();
        let e = Session::begin(&dir, "scf", Path::new("in.toml"), Value::Null)
            .err()
            .unwrap();
        assert!(
            e.contains("is running") && e.contains("baandi watch"),
            "{e}"
        );
        info.host = "elsewhere.invalid".into();
        fs::write(dir.join(LOCK), serde_json::to_vec(&info).unwrap()).unwrap();
        let e = Session::begin(&dir, "scf", Path::new("in.toml"), Value::Null)
            .err()
            .unwrap();
        assert!(e.contains("elsewhere.invalid"), "{e}");
        // A finished child's pid is a dead process on this host.
        let mut child = std::process::Command::new("true").spawn().unwrap();
        child.wait().unwrap();
        info.host = host().into();
        info.pid = child.id();
        assert!(!pid_alive(info.pid));
        fs::write(dir.join(LOCK), serde_json::to_vec(&info).unwrap()).unwrap();
        let s = Session::begin(&dir, "scf", Path::new("in.toml"), Value::Null).unwrap();
        assert_eq!(read_lock(&dir).unwrap().pid, std::process::id());
        s.finish("completed", 0, None);
        let r = records(&dir);
        assert_eq!(r[1]["event"], "warning");
        assert!(r[1]["message"].as_str().unwrap().contains("stale lock"));
        fs::remove_dir_all(&dir).unwrap();
    }
}
