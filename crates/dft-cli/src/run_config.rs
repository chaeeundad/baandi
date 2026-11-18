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
