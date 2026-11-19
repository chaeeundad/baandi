//! RUN-CONFIG (2026-11-18): the run's switch configuration from the input
//! file's `[performance]` and `[profile]` tables and the environment
//! (environment > input > default), resolved once per run.
//!
//! `main` installs it as the process configuration before the first switch is
//! read (so code that still reads the legacy accessors, on any thread, sees the
//! input file's values) and checks that every MPI rank resolved the same
//! effective values before any optional collective. `run::execute_mode`
//! resolves it again from the parsed input and uses the process object when
//! the values agree (the CLI binary), or enters a thread scope otherwise (an
//! in-process caller that runs several inputs).
use dft_numerics::switches::{self, RunConfig, ScopeGuard, Table};

/// `(table, key, text)` triples of the `[performance]` and `[profile]` tables
/// of a TOML document. Booleans become `on`/`off`, numbers their text.
pub fn input_values(
    performance: &std::collections::BTreeMap<String, toml::Value>,
    profile: &std::collections::BTreeMap<String, toml::Value>,
) -> Result<Vec<(Table, String, String)>, String> {
    let mut out = Vec::new();
    for (table, values) in [(Table::Performance, performance), (Table::Profile, profile)] {
        for (key, value) in values {
            let text = match value {
                toml::Value::String(s) => s.clone(),
                toml::Value::Boolean(b) => (if *b { "on" } else { "off" }).to_string(),
                toml::Value::Integer(i) => i.to_string(),
                toml::Value::Float(f) => f.to_string(),
                other => {
                    return Err(format!(
                        "[{}] {key}: expected a string, boolean or number, found {}",
                        table.as_str(),
                        other.type_str()
                    ));
                }
            };
            out.push((table, key.clone(), text));
        }
    }
    Ok(out)
}

type SwitchTable = std::collections::BTreeMap<String, toml::Value>;

fn tables_of(doc: &toml::Table) -> Result<(SwitchTable, SwitchTable), String> {
    let table = |name: &str| -> Result<SwitchTable, String> {
        match doc.get(name) {
            None => Ok(Default::default()),
            Some(toml::Value::Table(t)) => {
                Ok(t.iter().map(|(k, v)| (k.clone(), v.clone())).collect())
            }
            Some(_) => Err(format!("[{name}] must be a table")),
        }
    };
    Ok((table("performance")?, table("profile")?))
}

/// Resolve the configuration of the input file at `path` (environment >
/// input > default). An unreadable or unparsable file gives `Ok(None)`: the
/// command reports that error itself later; an invalid switch table is an
/// error.
pub fn resolve_file(path: &std::path::Path) -> Result<Option<RunConfig>, String> {
    let Ok(text) = std::fs::read_to_string(path) else {
        return Ok(None);
    };
    let Ok(doc) = toml::from_str::<toml::Table>(&text) else {
        return Ok(None);
    };
    let (performance, profile) = tables_of(&doc)?;
    let input = input_values(&performance, &profile)?;
    RunConfig::resolve_env_and_input(&input).map(Some)
}

/// `main`: install the input file's configuration as the process
/// configuration (before the first switch is read).
pub fn install_for_input(path: Option<&std::path::Path>) -> Result<(), String> {
    let Some(path) = path else {
        return Ok(());
    };
    match resolve_file(path)? {
        Some(config) => switches::install_process(config),
        None => Ok(()),
    }
}

/// The configuration of one run: the process object when its effective values
/// equal the input's (the CLI binary installed it), else `config` entered as
/// the calling thread's scope (returned guard) so that this run's code paths,
/// including its worker-pool regions, read it.
pub fn for_run(config: &crate::config::Config) -> Result<(RunConfig, Option<ScopeGuard>), String> {
    let input = input_values(&config.performance, &config.profile)?;
    let resolved = RunConfig::resolve_env_and_input(&input)?;
    let current = RunConfig::current();
    if current.agreement_hash() == resolved.agreement_hash()
        && switches::TABLE
            .iter()
            .all(|s| current.source(s.id) == resolved.source(s.id) || !switches::in_agreement(s.id))
    {
        Ok((current, None))
    } else {
        let guard = resolved.enter();
        Ok((resolved, Some(guard)))
    }
}

/// Cross-rank agreement (one allgatherv of three doubles per rank, before any
/// optional collective): every rank must have resolved the same effective
/// values of the agreement set. On a mismatch every rank returns the same
/// error, listing the hashes by rank. Every rank adopts rank 0's run ID.
pub fn agree<C: dft_mpi::Collective + ?Sized>(comm: &C, config: RunConfig) -> Result<(), String> {
    let ranks = comm.size();
    if ranks <= 1 {
        return Ok(());
    }
    let split = |x: u64| [(x >> 32) as f64, (x & 0xffff_ffff) as f64];
    let h = split(config.agreement_hash());
    let id = config.run_id();
    let local = [h[0], h[1], ((id >> 32) as f64)];
    let mut all = vec![0.0; 3 * ranks];
    comm.allgatherv_f64(&local, &vec![3; ranks], &mut all)
        .map_err(|e| format!("switch agreement exchange: {e}"))?;
    let hashes: Vec<u64> = all
        .chunks(3)
        .map(|c| ((c[0] as u64) << 32) | (c[1] as u64))
        .collect();
    // Rank 0's run ID: the high half from the exchange, the low half from a
    // second word (the exchange carries 32 bits per double exactly).
    let mut low = [(id & 0xffff_ffff) as f64];
    comm.bcast_f64(&mut low, 0)
        .map_err(|e| format!("run ID broadcast: {e}"))?;
    config.adopt_run_id(((all[2] as u64) << 32) | (low[0] as u64));
    if hashes.iter().any(|h| *h != hashes[0]) {
        let list: Vec<String> = hashes
            .iter()
            .enumerate()
            .map(|(r, h)| format!("rank {r}: {h:016x}"))
            .collect();
        return Err(format!(
            "switch configuration differs between MPI ranks (environment or input file); every rank must resolve the same [performance]/[profile] values: {}",
            list.join(", ")
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use dft_numerics::switches::{Id, Source};

    #[test]
    fn toml_values_map_to_switch_texts() {
        let doc: toml::Table = toml::from_str(
            "[performance]\nm2_grid = false\neigen_factor = 0.05\nladder_secant = 4\n[profile]\nphase_profile = \"on\"\n",
        )
        .unwrap();
        let (p, q) = tables_of(&doc).unwrap();
        let input = input_values(&p, &q).unwrap();
        let c = RunConfig::from_input(&input).unwrap();
        assert_eq!(c.choice(Id::M2_GRID), "off");
        assert_eq!(c.number(Id::EIGEN_FACTOR), Some(0.05));
        assert_eq!(c.choice(Id::LADDER_SECANT), "4");
        assert!(c.on(Id::PHASE_PROFILE));
        assert_eq!(c.source(Id::PHASE_PROFILE), Source::Input);
        let bad: toml::Table = toml::from_str("[performance]\nm2_grid = [1]\n").unwrap();
        let (p, q) = tables_of(&bad).unwrap();
        assert!(
            input_values(&p, &q)
                .unwrap_err()
                .contains("expected a string")
        );
    }

    /// Simulated ranks (threads with a shared communicator): equal
    /// configurations agree and adopt rank 0's run ID; a mismatch is the same
    /// error on every rank.
    #[test]
    fn ranks_agree_or_fail_together() {
        use dft_mpi::Collective;
        for mismatch in [false, true] {
            let out = dft_mpi::simulate(3, |comm| {
                let value = if mismatch && comm.rank() == 2 {
                    "off"
                } else {
                    "on"
                };
                let c = RunConfig::from_pairs(&[("m2_grid", value)]).unwrap();
                (agree(comm, c), c.run_id())
            });
            if mismatch {
                let first = out[0].0.as_ref().unwrap_err();
                assert!(first.contains("differs between MPI ranks"));
                assert!(out.iter().all(|(r, _)| r.as_ref().unwrap_err() == first));
            } else {
                assert!(out.iter().all(|(r, _)| r.is_ok()));
            }
            assert!(out.iter().all(|(_, id)| *id == out[0].1), "run ID adopted");
        }
    }
}

/// CLI lifecycle clocks (F-LAUNCH-IO). No numerical clocks or collectives are
/// added. Each rank appends its own record AFTER MPI_Finalize to its existing
/// results file; this avoids another gather/barrier or a second serialization
/// of the large scientific result. Library callers do not own process lifetime.
pub mod launch {
    use serde_json::{Value, json};
    use std::{
        cell::RefCell,
        collections::BTreeMap,
        fs::OpenOptions,
        io::{self, Read, Seek, SeekFrom, Write},
        path::{Path, PathBuf},
        time::Instant,
    };

    struct Launch {
        start: Instant,
        seconds: BTreeMap<&'static str, f64>,
        milestones: BTreeMap<&'static str, f64>,
        covered: Option<f64>,
        prepare_reused: Option<bool>,
        output: Option<PathBuf>,
    }
    thread_local! {
        static STATE: RefCell<Option<Launch>> = const { RefCell::new(None) };
    }

    /// Installed only by main, after resolving input and environment switches.
    /// The four bootstrap stamps are unavoidable even with profiling off: an
    /// input [profile] table is not available at process entry. Thereafter off
    /// means no clock reads, allocation, file writes or collective calls.
    pub fn begin(start: Instant, mpi_ready: Instant, args_ready: Instant, config_ready: Instant) {
        if !dft_engine::profile::enabled() {
            return;
        }
        STATE.with(|s| {
            *s.borrow_mut() = Some(Launch {
                start,
                seconds: BTreeMap::from([
                    ("mpi_init", (mpi_ready - start).as_secs_f64()),
                    ("arguments", (args_ready - mpi_ready).as_secs_f64()),
                    ("switch_input", (config_ready - args_ready).as_secs_f64()),
                ]),
                milestones: BTreeMap::new(),
                covered: None,
                prepare_reused: None,
                output: None,
            })
        });
    }
    pub(crate) fn clock() -> Option<Instant> {
        STATE.with(|s| s.borrow().as_ref().map(|_| Instant::now()))
    }
    pub(crate) fn elapsed(key: &'static str, clock: Option<Instant>) {
        if let Some(clock) = clock {
            let seconds = clock.elapsed().as_secs_f64();
            STATE.with(|s| {
                if let Some(s) = s.borrow_mut().as_mut() {
                    *s.seconds.entry(key).or_default() += seconds;
                }
            });
        }
    }
    pub fn mark(key: &'static str) {
        STATE.with(|s| {
            if let Some(s) = s.borrow_mut().as_mut() {
                s.milestones.insert(key, s.start.elapsed().as_secs_f64());
            }
        });
    }
    pub(crate) fn since(key: &'static str, milestone: &'static str) {
        STATE.with(|s| {
            if let Some(s) = s.borrow_mut().as_mut()
                && let Some(t) = s.milestones.get(milestone)
            {
                s.seconds.insert(key, s.start.elapsed().as_secs_f64() - t);
            }
        });
    }
    pub(crate) fn covered(seconds: f64) {
        STATE.with(|s| {
            if let Some(s) = s.borrow_mut().as_mut() {
                s.covered = Some(seconds);
            }
        });
    }
    pub(crate) fn output(path: &Path) {
        STATE.with(|s| {
            if let Some(s) = s.borrow_mut().as_mut() {
                s.output = Some(path.to_owned());
            }
        });
    }
    pub(crate) fn prepared(reused: bool) {
        STATE.with(|s| {
            if let Some(s) = s.borrow_mut().as_mut() {
                s.prepare_reused = Some(reused);
            }
        });
    }
    /// Time an operation from the binary without exposing the private clock.
    pub fn timed<T>(key: &'static str, f: impl FnOnce() -> T) -> T {
        let start = clock();
        let result = f();
        elapsed(key, start);
        result
    }
    fn record(s: &Launch) -> Value {
        let wall = s.start.elapsed().as_secs_f64();
        let entry = s.milestones.get("run_entry").copied();
        let finish = s.milestones.get("profile_finish").copied();
        json!({"schema":"launch-1", "rank":crate::parallel::rank(), "ranks":crate::parallel::size(),
            "mpi_backend":if cfg!(feature="mpi") { "mpi" } else { "single_rank_stub" },
            "prepare_reused":s.prepare_reused,
            "seconds":s.seconds, "milestones_s":s.milestones,
            "main_wall_s":wall, "covered_s":s.covered,
            "outside_covered_s":s.covered.map(|covered| wall-covered),
            "closure": {"before_init_s":entry, "after_profile_s":finish.map(|t| wall-t),
                "within_profile_residual_s":entry.zip(finish).zip(s.covered).map(|((a,b),c)| b-a-c)},
            "note":"Local monotonic seconds from main entry; no cross-rank clock subtraction. init.cli already includes input/preflight/UPF; engine_construction and materialize are also inside covered_s. Detail timers are nested, not additive. preflight includes Config::prepare, its UPF parse and per-k basis checks; prepare_reuse only takes that prepared value. Binary SHA is inside run_block. Each rank retains its existing output files; read rank*/results.json for other ranks. main_wall_s ends before this small launch append; external wall minus main_wall_s includes pre-main launcher/loader, final append and post-main runtime/OS teardown. No fsync or new barrier. MPI finalize time includes waiting for slower ranks."})
    }
    fn append(path: &Path, launch: &Value) -> io::Result<()> {
        // Preserve every original byte except the final top-level brace. The
        // file was emitted by serde_json::to_vec_pretty and has no newline.
        let mut file = OpenOptions::new().read(true).write(true).open(path)?;
        file.seek(SeekFrom::End(-1))?;
        let mut last = [0];
        file.read_exact(&mut last)?;
        if last != *b"}" {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "results.json must end with an object",
            ));
        }
        let mut bytes = b",\n  \"launch\": ".to_vec();
        serde_json::to_writer_pretty(&mut bytes, launch)?;
        bytes.extend_from_slice(b"\n}");
        file.seek(SeekFrom::End(-1))?;
        file.write_all(&bytes)?;
        file.flush()
    }
    /// Called after dropping the engine communicator and the MPI universe.
    /// A missing output (unsupported driver or early failure) adds no file.
    pub fn finish() -> io::Result<()> {
        STATE.with(|s| {
            if let Some(s) = s.borrow_mut().take()
                && let Some(path) = &s.output {
                let record = record(&s);
                append(path, &record)?;
                eprintln!("launch rank {}: main {:.6} s; outside covered {:.6} s (local; see results.json)",
                    crate::parallel::rank(), record["main_wall_s"].as_f64().unwrap(), record["outside_covered_s"].as_f64().unwrap_or(0.));
            }
            Ok(())
        })
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        #[test]
        fn append_preserves_science_and_bytes_and_rejects_truncation() {
            let path = std::env::temp_dir().join(format!("launch-{}.json", std::process::id()));
            let science = json!({"density":"abc", "run":{"completed":true}, "total_energy":{"ha":-1.234567890123456}});
            let before = serde_json::to_vec_pretty(&science).unwrap();
            std::fs::write(&path, &before).unwrap();
            let profile = json!({"schema":"launch-1", "seconds":{"mpi_finalize":0.25}});
            append(&path, &profile).unwrap();
            let after = std::fs::read(&path).unwrap();
            assert_eq!(&after[..before.len() - 1], &before[..before.len() - 1]);
            let mut parsed: Value = serde_json::from_slice(&after).unwrap();
            assert_eq!(
                parsed.as_object_mut().unwrap().remove("launch"),
                Some(profile)
            );
            assert_eq!(parsed, science);
            std::fs::write(&path, b"{\"truncated\":").unwrap();
            assert!(append(&path, &json!({})).is_err());
            std::fs::remove_file(path).unwrap();
        }
        #[test]
        fn inactive_has_no_clock_or_output() {
            assert!(clock().is_none());
            mark("unused");
            output(Path::new("/not/an/output"));
            finish().unwrap();
            assert!(clock().is_none());
        }
        #[test]
        fn closure_uses_this_ranks_covered_wall_and_clears_state() {
            let start = Instant::now() - std::time::Duration::from_secs(1);
            let state = Launch {
                start,
                seconds: BTreeMap::new(),
                milestones: BTreeMap::from([("run_entry", 0.1), ("profile_finish", 0.8)]),
                covered: Some(0.6),
                prepare_reused: Some(true),
                output: None,
            };
            let r = record(&state);
            let closure = &r["closure"];
            let sum = closure["before_init_s"].as_f64().unwrap()
                + closure["after_profile_s"].as_f64().unwrap()
                + closure["within_profile_residual_s"].as_f64().unwrap();
            assert!((sum - r["outside_covered_s"].as_f64().unwrap()).abs() < 1e-14);
            STATE.with(|s| *s.borrow_mut() = Some(state));
            assert!(clock().is_some());
            finish().unwrap();
            assert!(clock().is_none());
        }
    }
}
