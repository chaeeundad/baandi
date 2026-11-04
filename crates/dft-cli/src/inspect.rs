//! Read-only explanation of recorded evidence. No inferred convergence or causes.
use crate::{config::Config, grid_report};
use serde_json::{Value, json};
use std::{
    fs::File,
    io::{BufRead, BufReader, Read},
    path::{Path, PathBuf},
};
type Error = Box<dyn std::error::Error>;

fn source(path: &Path) -> PathBuf {
    if path.is_dir() {
        // Prefer the event stream: a previous result must not hide a newer run.
        let log = path.join("events.jsonl");
        if log.exists() {
            log
        } else {
            path.join("results.json")
        }
    } else {
        path.into()
    }
}
fn evidence(path: &Path) -> Result<Value, Error> {
    let path = source(path);
    if path.extension().is_some_and(|e| e == "jsonl") {
        let mut reader = BufReader::new(File::open(&path)?);
        let mut last = None;
        let mut sequence = 0;
        let mut calc = None;
        let mut partial = false;
        let mut last_iteration = None;
        let mut last_certificate = None;
        loop {
            let mut bytes = Vec::new();
            let n = reader.by_ref().take(65537).read_until(b'\n', &mut bytes)?;
            if n == 0 {
                break;
            }
            if n > 65536 {
                return Err("event: maximum record size is 65536 bytes".into());
            }
            if bytes.last() != Some(&b'\n') {
                partial = true;
                break;
            }
            let value: Value = serde_json::from_slice(&bytes)?;
            if value["schema_version"] != "1.0" {
                return Err("event.schema_version: supported value is 1.0".into());
            }
            let seq = value["sequence"]
                .as_u64()
                .ok_or("event.sequence: expected positive integer")?;
            let id = value["calc_id"]
                .as_str()
                .ok_or("event.calc_id: expected string")?
                .to_owned();
            if seq <= sequence || calc.as_ref().is_some_and(|c| c != &id) {
                return Err(
                    "event log: sequence must increase and calc_id must remain unchanged".into(),
                );
            }
            if value["iteration"].as_u64().is_none()
                || value["state_version"].as_u64().is_none()
                || !value["data"].is_object()
            {
                return Err("event: iteration/state_version must be nonnegative integers and data must be an object".into());
            }
            match value["type"].as_str() {
                Some("iteration_finished") => last_iteration = Some(value.clone()),
                Some("certificate_evaluated") => last_certificate = Some(value.clone()),
                Some(
                    "initialized" | "mixer_reset" | "diverged" | "converged" | "not_converged"
                    | "failed" | "events_dropped",
                ) => {}
                _ => return Err("event.type: unsupported event type".into()),
            }
            sequence = seq;
            calc = Some(id);
            last = Some(value);
        }
        let state = match last.as_ref().and_then(|v| v["type"].as_str()) {
            Some("converged") => "converged",
            Some("not_converged") => "not_converged",
            Some("failed" | "diverged") => "failed",
            _ => "nonterminal: process liveness unknown",
        };
        let mut result = json!({"source":path,"state":state,"partial_tail_ignored":partial,"last_event":last,
            "last_iteration":last_iteration,"last_certificate":last_certificate});
        let audit_path = path.with_file_name("control-audit.jsonl");
        if audit_path.exists() {
            let mut reader = BufReader::new(File::open(&audit_path)?);
            let mut last_control = None;
            loop {
                let mut bytes = Vec::new();
                let n = reader.by_ref().take(65537).read_until(b'\n', &mut bytes)?;
                if n == 0 {
                    break;
                }
                if n > 65536 {
                    return Err("audit: maximum record size is 65536 bytes".into());
                }
                if bytes.last() != Some(&b'\n') {
                    result["partial_audit_tail_ignored"] = json!(true);
                    break;
                }
                let record: Value = serde_json::from_slice(&bytes)?;
                if record["schema_version"] != 1 || record["calc_id"].as_str() != calc.as_deref() {
                    return Err("audit: schema or calculation identity mismatch".into());
                }
                last_control = Some(record);
            }
            if let Some(record) = &last_control
                && record["state_version"].as_u64() > result["last_event"]["state_version"].as_u64()
                && record["outcome"] == "Applied"
            {
                result["state"] = if record["state"] == "Cancelled" {
                    json!("cancelled")
                } else {
                    json!("nonterminal: process liveness unknown")
                };
            }
            result["last_control"] = json!(last_control);
        }
        Ok(result)
    } else {
        let value: Value = serde_json::from_reader(File::open(&path)?)?;
        if value["schema_version"] != 1 || !value["convergence"]["converged"].is_boolean() {
            return Err("source: expected CLI results.json (schema_version=1) or events.jsonl; checkpoint inspection is not available in this adapter".into());
        }
        Ok(
            json!({"source":path,"state":if value["convergence"]["converged"]==true {"converged"} else {"not_converged"},
            "convergence":value["convergence"],"energies":value["energies"],"preflight":value["preflight"]}),
        )
    }
}
pub fn status(path: &Path) -> Result<(), Error> {
    println!("{}", serde_json::to_string_pretty(&evidence(path)?)?);
    Ok(())
}
pub fn explain(path: &Path) -> Result<(), Error> {
    if path.extension().is_some_and(|e| e == "toml") {
        let config = Config::parse(&std::fs::read_to_string(path)?)?;
        let preflight = crate::config::preflight::check(&config, path)?;
        let prepared = config.prepare()?;
        print!("{}", grid_report(&config, &prepared));
        println!("XC: {}", config.xc.functional.engine().name());
        preflight.print();
        println!(
            "Resolved SCF controls (energies in Ha): {:#?}",
            crate::run::options(&config)?
        );
        println!(
            "Preflight validates supported input and resource estimates; it does not establish numerical convergence or physical accuracy."
        );
    } else {
        println!(
            "Recorded evidence (energies in Ha; normalized density residual in Bohr^-3/2):\n{}",
            serde_json::to_string_pretty(&evidence(path)?)?
        );
        println!(
            "Diagnostics: convergence is reported only from recorded engine evidence. A nonterminal event cannot establish whether a process is still running. Hypotheses: none inferred; a residual or exhausted budget alone does not establish a physical cause."
        );
    }
    Ok(())
}
