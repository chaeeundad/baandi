//! UX1 job monitor (ADR 0053): folds a `progress.jsonl` stream
//! (`crate::progress`, schema `baandi-progress-1`) into a view state and
//! renders it as a full-screen terminal UI (`monitor::tui`, ratatui) or as a
//! plain, refreshing text view (`--plain`, or when stdout is not a terminal).
//! It only reads files: no IPC with the running job.
use crate::progress;
use serde_json::Value;
use std::{
    fmt::Write as _,
    fs,
    io::{Read, Seek, SeekFrom},
    path::{Path, PathBuf},
};

mod tui;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Map {
    pub iteration: u64,
    pub kind: String,
    pub energy_ha: Option<f64>,
    pub de_ha: Option<f64>,
    pub dr2_ry: f64,
    pub seconds: f64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Cycle {
    pub index: u64,
    pub threshold_ry: f64,
    pub maps: Vec<Map>,
    /// `(maps, converged)` once `scf_end` arrived.
    pub end: Option<(u64, Option<bool>)>,
    pub error: Option<String>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Step {
    pub step: u64,
    pub energy_ha: f64,
    pub max_force_ha_bohr: f64,
    pub step_length_bohr: f64,
    pub scf_iterations: u64,
    pub provenance: String,
    pub pressure_ha_bohr3: Option<f64>,
    pub seconds: f64,
}

/// Everything the views show, folded from the stream in order.
#[derive(Clone, Debug, Default)]
pub struct State {
    pub start: Option<Value>,
    pub phase: String,
    pub cycles: Vec<Cycle>,
    pub steps: Vec<Step>,
    /// `(elapsed seconds, text)`, oldest first.
    pub events: Vec<(f64, String)>,
    pub end: Option<Value>,
    pub elapsed_s: f64,
    pub start_unix: Option<f64>,
    pub records: usize,
    /// Lines that did not parse (a torn last line is not counted: it is kept
    /// until its newline arrives).
    pub invalid: usize,
    last_step_elapsed: f64,
}

const MAX_EVENTS: usize = 500;

fn num(v: &Value) -> Option<f64> {
    v.as_f64()
}
fn sci(v: Option<f64>) -> String {
    v.map_or_else(|| "-".into(), |x| format!("{x:.2e}"))
}

impl State {
    pub fn parse(text: &str) -> Self {
        let mut s = Self::default();
        for line in text.lines() {
            s.line(line);
        }
        s
    }
    fn line(&mut self, line: &str) {
        if line.trim().is_empty() {
            return;
        }
        match serde_json::from_str::<Value>(line) {
            Ok(v) if v["v"] == 1 => self.apply(&v),
            _ => self.invalid += 1,
        }
    }
    fn event(&mut self, text: String) {
        self.events.push((self.elapsed_s, text));
        if self.events.len() > MAX_EVENTS {
            self.events.remove(0);
        }
    }
    pub fn apply(&mut self, v: &Value) {
        self.records += 1;
        if let Some(e) = num(&v["elapsed_s"]) {
            self.elapsed_s = e;
        }
        let event = v["event"].as_str().unwrap_or("");
        match event {
            "start" => {
                self.start_unix = num(&v["t"]);
                self.start = Some(v.clone());
                self.phase = "setup".into();
            }
            "phase" => {
                self.phase = v["phase"].as_str().unwrap_or("?").into();
            }
            "scf_start" => self.cycles.push(Cycle {
                index: v["cycle"].as_u64().unwrap_or(0),
                threshold_ry: num(&v["threshold_ry"]).unwrap_or(f64::NAN),
                ..Cycle::default()
            }),
            "scf_map" => {
                let map = Map {
                    iteration: v["iteration"].as_u64().unwrap_or(0),
                    kind: v["kind"].as_str().unwrap_or("?").into(),
                    energy_ha: num(&v["energy_ha"]),
                    de_ha: num(&v["de_ha"]),
                    dr2_ry: num(&v["dr2_ry"]).unwrap_or(f64::NAN),
                    seconds: num(&v["seconds"]).unwrap_or(0.),
                };
                if map.kind != "ordinary" {
                    self.event(format!(
                        "{} map {}: dr2 {} Ry",
                        map.kind,
                        map.iteration,
                        sci(Some(map.dr2_ry))
                    ));
                }
                if self.cycles.is_empty() {
                    self.cycles.push(Cycle::default());
                }
                self.cycles.last_mut().unwrap().maps.push(map);
            }
            "scf_end" => {
                let maps = v["maps"].as_u64().unwrap_or(0);
                let converged = v["converged"].as_bool();
                let error = v["error"].as_str().map(str::to_owned);
                let text = match (&error, converged) {
                    (Some(e), _) => format!("SCF {} failed: {e}", v["cycle"]),
                    (None, Some(false)) => {
                        format!("SCF {} not converged after {maps} maps", v["cycle"])
                    }
                    _ => format!("SCF {} finished in {maps} maps", v["cycle"]),
                };
                self.event(text);
                if let Some(c) = self.cycles.last_mut() {
                    c.end = Some((maps, converged));
                    c.error = error;
                }
            }
            "relax_step" => {
                let step = Step {
                    step: v["step"].as_u64().unwrap_or(0),
                    energy_ha: num(&v["energy_ha"]).unwrap_or(f64::NAN),
                    max_force_ha_bohr: num(&v["max_force_ha_bohr"]).unwrap_or(f64::NAN),
                    step_length_bohr: num(&v["step_length_bohr"]).unwrap_or(0.),
                    scf_iterations: v["scf_iterations"].as_u64().unwrap_or(0),
                    provenance: v["provenance"].as_str().unwrap_or("").into(),
                    pressure_ha_bohr3: num(&v["pressure_ha_bohr3"]),
                    seconds: self.elapsed_s - self.last_step_elapsed,
                };
                self.last_step_elapsed = self.elapsed_s;
                self.event(format!(
                    "relax step {}: E {:.10} Ha, max|F| {} Ha/Bohr",
                    step.step,
                    step.energy_ha,
                    sci(Some(step.max_force_ha_bohr))
                ));
                self.steps.push(step);
            }
            "warning" => {
                let text = format!(
                    "warning ({}): {}",
                    v["source"].as_str().unwrap_or("?"),
                    v["message"].as_str().unwrap_or("")
                );
                self.event(text);
            }
            "end" => {
                let text = format!(
                    "end: {}{}",
                    v["status"].as_str().unwrap_or("?"),
                    v["message"]
                        .as_str()
                        .map(|m| format!(" ({m})"))
                        .unwrap_or_default()
                );
                self.event(text);
                self.end = Some(v.clone());
                self.phase = "finished".into();
            }
            _ => {}
        }
    }
    pub fn task(&self) -> &str {
        self.start
            .as_ref()
            .and_then(|s| s["task"].as_str())
            .unwrap_or("?")
    }
    fn system(&self) -> &Value {
        static NULL: Value = Value::Null;
        self.start.as_ref().map_or(&NULL, |s| &s["system"])
    }
    /// One-line system summary.
    pub fn system_line(&self) -> String {
        let s = self.system();
        if s.is_null() {
            return "(input not parsed)".into();
        }
        let mut parts = vec![format!(
            "{} ({})",
            s["formula"].as_str().unwrap_or("?"),
            s["name"].as_str().unwrap_or("")
        )];
        parts.push(match s["atoms"].as_u64() {
            Some(1) => "1 atom".into(),
            Some(n) => format!("{n} atoms"),
            None => "atoms ?".into(),
        });
        parts.push(match s["kpoints"].as_u64() {
            Some(1) => "1 k point".into(),
            Some(k) => format!("{k} k points"),
            None => "k ?".into(),
        });
        parts.push(match s["nbands"].as_u64() {
            Some(n) => format!("{n} bands"),
            None => "bands auto".into(),
        });
        if let Some(e) = num(&s["ecutwfc_ry"]) {
            parts.push(format!("{e} Ry"));
        }
        let ranks = s["ranks"].as_u64().unwrap_or(1);
        parts.push(match s["threads"].as_u64() {
            Some(t) => format!("{ranks} rank(s) x {t} threads"),
            None => format!("{ranks} rank(s)"),
        });
        if s["spin"] == true {
            parts.push("collinear spin".into());
        }
        if s["symmetry"] == "auto" {
            parts.push("symmetry auto".into());
        }
        parts.join(" | ")
    }
    pub fn relax_settings(&self) -> Option<&Value> {
        let r = &self.system()["relax"];
        (!r.is_null()).then_some(r)
    }
    /// Human phase line (`SCF 2 (relax step 1), map 5 [verification]`).
    pub fn phase_line(&self) -> String {
        if let Some(end) = &self.end {
            return format!("finished: {}", end["status"].as_str().unwrap_or("?"));
        }
        let relax = (self.task().contains("relax")).then(|| {
            format!(
                " (relax step {})",
                self.steps.last().map_or(0, |s| s.step + 1)
            )
        });
        match (self.phase.as_str(), self.cycles.last()) {
            ("scf", Some(c)) => {
                let last = c.maps.last();
                format!(
                    "SCF {}{}, map {}{}",
                    c.index,
                    relax.unwrap_or_default(),
                    last.map_or(0, |m| m.iteration),
                    last.map(|m| format!(" [{}]", m.kind)).unwrap_or_default()
                )
            }
            ("post-scf", Some(c)) => format!(
                "after SCF {}{}: forces/stress/output",
                c.index,
                relax.unwrap_or_default()
            ),
            (p, _) => p.into(),
        }
    }
}

/// Whether the job of a directory is still running.
#[derive(Clone, Debug, PartialEq)]
pub enum JobStatus {
    Running {
        pid: u32,
    },
    /// Lock from another host; liveness cannot be checked here.
    Remote {
        host: String,
        pid: u32,
    },
    Finished {
        status: String,
        exit_code: Option<u64>,
    },
    /// No final record and no live process (killed, or a stale lock).
    Stopped,
}
impl JobStatus {
    pub fn label(&self) -> String {
        match self {
            Self::Running { pid } => format!("RUNNING pid {pid}"),
            Self::Remote { host, pid } => format!("RUNNING? pid {pid} on {host}"),
            Self::Finished {
                status,
                exit_code: Some(c),
            } => format!("{} (exit {c})", status.to_uppercase()),
            Self::Finished { status, .. } => status.to_uppercase(),
            Self::Stopped => "STOPPED (no final record, process gone)".into(),
        }
    }
    pub fn live(&self) -> bool {
        matches!(self, Self::Running { .. } | Self::Remote { .. })
    }
}

pub fn job_status(dir: &Path, state: &State) -> JobStatus {
    if let Some(end) = &state.end {
        return JobStatus::Finished {
            status: end["status"].as_str().unwrap_or("?").into(),
            exit_code: end["exit_code"].as_u64(),
        };
    }
    match progress::read_lock(dir) {
        Some(lock) => match progress::liveness(&lock) {
            progress::Liveness::Alive => JobStatus::Running { pid: lock.pid },
            progress::Liveness::Remote(host) => JobStatus::Remote {
                host,
                pid: lock.pid,
            },
            progress::Liveness::Dead => JobStatus::Stopped,
        },
        None => JobStatus::Stopped,
    }
}

/// Whether `dir` holds the lock of a live job (bare `baandi` then watches).
pub fn live_job(dir: &Path) -> bool {
    progress::read_lock(dir).is_some_and(|l| {
        matches!(
            progress::liveness(&l),
            progress::Liveness::Alive | progress::Liveness::Remote(_)
        )
    })
}

/// Unicode sparkline of finite values (lowest = ' ' bar).
pub fn sparkline(values: &[f64]) -> String {
    const BARS: [char; 8] = ['▁', '▂', '▃', '▄', '▅', '▆', '▇', '█'];
    let finite: Vec<f64> = values.iter().copied().filter(|v| v.is_finite()).collect();
    let (lo, hi) = finite
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &v| {
            (a.min(v), b.max(v))
        });
    values
        .iter()
        .map(|v| {
            if !v.is_finite() {
                ' '
            } else if hi > lo {
                BARS[(((v - lo) / (hi - lo)) * 7.).round() as usize]
            } else {
                BARS[3]
            }
        })
        .collect()
}

/// Plain text view. `now_unix` adds the live elapsed time; `None` uses the
/// last record's elapsed time (deterministic, for finished jobs and tests).
pub fn render_plain(
    dir: &Path,
    state: &State,
    status: &JobStatus,
    now_unix: Option<f64>,
) -> String {
    let mut o = String::new();
    let elapsed = match (now_unix, state.start_unix, status.live()) {
        (Some(now), Some(start), true) => now - start,
        _ => state.elapsed_s,
    };
    let start = state.start.as_ref();
    let _ = writeln!(
        o,
        "BAANDI {}  [{}]  elapsed {:.1} s",
        dir.display(),
        status.label(),
        elapsed
    );
    let _ = writeln!(
        o,
        "Task: {}   Input: {}",
        state.task(),
        start.and_then(|s| s["input"].as_str()).unwrap_or("?")
    );
    let _ = writeln!(o, "System: {}", state.system_line());
    let _ = writeln!(o, "Phase: {}", state.phase_line());
    if let Some(c) = state.cycles.last() {
        let _ = writeln!(
            o,
            "\nSCF {} of {}  stop: dr2 < {} Ry",
            c.index,
            state.cycles.len(),
            sci(Some(c.threshold_ry))
        );
        let _ = writeln!(
            o,
            "  {:>4}  {:<12} {:>20} {:>10} {:>10} {:>8}",
            "map", "kind", "energy (Ha)", "dE (Ha)", "dr2 (Ry)", "time (s)"
        );
        let skip = c.maps.len().saturating_sub(15);
        if skip > 0 {
            let _ = writeln!(o, "  ... {skip} earlier maps");
        }
        for m in &c.maps[skip..] {
            let _ = writeln!(
                o,
                "  {:>4}  {:<12} {:>20} {:>10} {:>10} {:>8.3}",
                m.iteration,
                m.kind,
                m.energy_ha
                    .map_or_else(|| "-".into(), |e| format!("{e:.12}")),
                sci(m.de_ha),
                sci(Some(m.dr2_ry)),
                m.seconds
            );
        }
        let logs: Vec<f64> = c.maps.iter().map(|m| m.dr2_ry.log10()).collect();
        let _ = writeln!(
            o,
            "  log10 dr2: {}  (threshold {:.1})",
            sparkline(&logs),
            c.threshold_ry.log10()
        );
        if let Some((maps, converged)) = c.end {
            let _ = writeln!(
                o,
                "  SCF {} ended after {maps} maps{}",
                c.index,
                match converged {
                    Some(true) => ", converged",
                    Some(false) => ", NOT converged",
                    None => "",
                }
            );
        }
    }
    if !state.steps.is_empty() || state.relax_settings().is_some() {
        let r = state.relax_settings();
        let _ = writeln!(
            o,
            "\nRelaxation  force tol {} Ha/Bohr, energy tol {} Ha, max steps {}",
            sci(r.and_then(|r| num(&r["force_tol_ha_bohr"]))),
            sci(r.and_then(|r| num(&r["energy_tol_ha"]))),
            r.map_or(Value::Null, |r| r["max_steps"].clone())
        );
        let _ = writeln!(
            o,
            "  {:>4}  {:>20} {:>12} {:>10} {:>5} {:>8}  provenance",
            "step", "energy (Ha)", "max|F| Ha/B", "step Bohr", "maps", "time (s)"
        );
        let skip = state.steps.len().saturating_sub(15);
        if skip > 0 {
            let _ = writeln!(o, "  ... {skip} earlier steps");
        }
        for s in &state.steps[skip..] {
            let _ = writeln!(
                o,
                "  {:>4}  {:>20.12} {:>12} {:>10} {:>5} {:>8.3}  {}{}",
                s.step,
                s.energy_ha,
                sci(Some(s.max_force_ha_bohr)),
                sci(Some(s.step_length_bohr)),
                s.scf_iterations,
                s.seconds,
                s.provenance,
                s.pressure_ha_bohr3
                    .map(|p| format!(" P={p:.3e} Ha/Bohr^3"))
                    .unwrap_or_default()
            );
        }
        let forces: Vec<f64> = state
            .steps
            .iter()
            .map(|s| s.max_force_ha_bohr.log10())
            .collect();
        if !forces.is_empty() {
            let _ = writeln!(o, "  log10 max|F|: {}", sparkline(&forces));
        }
    }
    let _ = writeln!(o, "\nEvents (latest {}):", state.events.len().min(10));
    for (t, text) in state.events.iter().rev().take(10).rev() {
        let _ = writeln!(o, "  {t:>9.2} s  {text}");
    }
    if let Some(end) = &state.end {
        let _ = writeln!(
            o,
            "\nFinal: {} (exit {}) after {:.1} s",
            end["status"].as_str().unwrap_or("?"),
            end["exit_code"],
            state.elapsed_s
        );
    }
    if state.invalid > 0 {
        let _ = writeln!(o, "({} unreadable lines skipped)", state.invalid);
    }
    o
}

/// Incremental reader of a growing stream: complete lines only, restart from
/// the top when the file is replaced or truncated (a new job).
pub struct Tail {
    path: PathBuf,
    offset: u64,
    partial: String,
    pub state: State,
}
impl Tail {
    pub fn new(dir: &Path) -> Self {
        Self {
            path: dir.join(progress::FILE),
            offset: 0,
            partial: String::new(),
            state: State::default(),
        }
    }
    /// Read new complete lines; returns whether the state changed.
    pub fn poll(&mut self) -> bool {
        let Ok(mut file) = fs::File::open(&self.path) else {
            return false;
        };
        let len = file.metadata().map_or(0, |m| m.len());
        if len < self.offset {
            *self = Self::new(self.path.parent().unwrap_or(Path::new(".")));
        }
        if len == self.offset {
            return false;
        }
        let mut bytes = Vec::new();
        if file.seek(SeekFrom::Start(self.offset)).is_err() || file.read_to_end(&mut bytes).is_err()
        {
            return false;
        }
        self.offset += bytes.len() as u64;
        self.partial.push_str(&String::from_utf8_lossy(&bytes));
        let mut changed = false;
        while let Some(i) = self.partial.find('\n') {
            let line: String = self.partial.drain(..=i).collect();
            self.state.line(&line);
            changed = true;
        }
        changed
    }
}

/// The directory whose stream to show: `dir` itself when it holds a stream or
/// lock; else the output directory of `dir/input.toml`; else the most
/// recently modified immediate subdirectory holding a stream.
pub fn locate(dir: &Path) -> Result<PathBuf, String> {
    let has = |d: &Path| d.join(progress::FILE).is_file() || d.join(progress::LOCK).is_file();
    if has(dir) {
        return Ok(dir.to_owned());
    }
    let input = dir.join(crate::task::DEFAULT_INPUT);
    if let Ok(text) = fs::read_to_string(&input)
        && let Ok(value) = toml::from_str::<toml::Value>(&text)
    {
        let out = dir.join(crate::task::output_dir(&value));
        if has(&out) {
            return Ok(out);
        }
    }
    let mut best: Option<(std::time::SystemTime, PathBuf)> = None;
    for entry in fs::read_dir(dir)
        .map_err(|e| format!("{}: {e}", dir.display()))?
        .flatten()
    {
        let p = entry.path();
        if let Ok(m) = fs::metadata(p.join(progress::FILE))
            && let Ok(t) = m.modified()
            && best.as_ref().is_none_or(|(b, _)| t > *b)
        {
            best = Some((t, p));
        }
    }
    best.map(|(_, p)| p).ok_or_else(|| {
        format!(
            "no BAANDI job found in {} (no {} here, in the output directory of {}, or in a subdirectory)",
            dir.display(),
            progress::FILE,
            crate::task::DEFAULT_INPUT
        )
    })
}

/// `baandi watch`: TUI on a terminal, else (or with `plain`) a text view that
/// is reprinted when the stream changes until the job ends. `once` prints the
/// current view and returns.
pub fn watch(dir: &Path, plain: bool, once: bool) -> Result<(), String> {
    use std::io::IsTerminal;
    let dir = locate(dir)?;
    let mut tail = Tail::new(&dir);
    tail.poll();
    if !(plain || once) && std::io::stdout().is_terminal() && std::io::stdin().is_terminal() {
        return tui::run(&dir, tail).map_err(|e| format!("terminal UI: {e}"));
    }
    let now = || {
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0., |d| d.as_secs_f64())
    };
    let mut status = job_status(&dir, &tail.state);
    print!("{}", render_plain(&dir, &tail.state, &status, Some(now())));
    if once {
        return Ok(());
    }
    while status.live() {
        std::thread::sleep(std::time::Duration::from_secs(2));
        let changed = tail.poll();
        let next = job_status(&dir, &tail.state);
        if changed || next != status {
            status = next;
            println!("\n{}", "-".repeat(72));
            print!("{}", render_plain(&dir, &tail.state, &status, Some(now())));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const RELAX: &str = include_str!("monitor/fixtures/relax.jsonl");
    const RELAX_VIEW: &str = include_str!("monitor/fixtures/relax.plain.txt");

    #[test]
    fn recorded_relax_stream_folds_and_renders_the_snapshot() {
        let state = State::parse(RELAX);
        assert_eq!(state.invalid, 0);
        assert_eq!(state.records, 48);
        assert_eq!(state.cycles.len(), 4);
        assert_eq!(state.cycles[0].maps.len(), 7);
        assert_eq!(state.steps.len(), 4);
        let status = job_status(Path::new("/nonexistent"), &state);
        assert_eq!(
            status,
            JobStatus::Finished {
                status: "not_converged".into(),
                exit_code: Some(3)
            }
        );
        let view = render_plain(Path::new("results"), &state, &status, None);
        if std::env::var_os("UPDATE_SNAPSHOTS").is_some() {
            fs::write(
                concat!(
                    env!("CARGO_MANIFEST_DIR"),
                    "/src/monitor/fixtures/relax.plain.txt"
                ),
                &view,
            )
            .unwrap();
        }
        assert_eq!(view, RELAX_VIEW, "UPDATE_SNAPSHOTS=1 regenerates");
    }

    #[test]
    fn torn_lines_wait_and_a_running_stream_shows_the_live_phase() {
        let dir = std::env::temp_dir().join(format!("baandi-monitor-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        let lines: Vec<&str> = RELAX.lines().collect();
        // First nine records and half of the tenth: the torn line is pending.
        let mut text = lines[..9].join("\n");
        text.push('\n');
        text.push_str(&lines[9][..lines[9].len() / 2]);
        fs::write(dir.join(progress::FILE), &text).unwrap();
        let mut tail = Tail::new(&dir);
        assert!(tail.poll());
        assert_eq!(tail.state.records, 9);
        assert_eq!(tail.state.invalid, 0);
        assert_eq!(job_status(&dir, &tail.state), JobStatus::Stopped);
        assert!(
            tail.state
                .phase_line()
                .starts_with("SCF 1 (relax step 0), map")
        );
        text.push_str(&lines[9][lines[9].len() / 2..]);
        text.push('\n');
        fs::write(dir.join(progress::FILE), &text).unwrap();
        assert!(tail.poll());
        assert_eq!(tail.state.records, 10);
        assert_eq!(locate(&dir).unwrap(), dir);
        let parent = dir
            .parent()
            .unwrap()
            .join(format!("baandi-monitor-p-{}", std::process::id()));
        let _ = fs::remove_dir_all(&parent);
        fs::create_dir_all(parent.join("out")).unwrap();
        fs::write(parent.join("out").join(progress::FILE), "").unwrap();
        fs::write(parent.join("input.toml"), "[output]\ndirectory='out'\n").unwrap();
        assert_eq!(locate(&parent).unwrap(), parent.join("out"));
        fs::remove_file(parent.join("input.toml")).unwrap();
        assert_eq!(locate(&parent).unwrap(), parent.join("out"));
        fs::remove_dir_all(&parent).unwrap();
        fs::remove_dir_all(&dir).unwrap();
    }
}
