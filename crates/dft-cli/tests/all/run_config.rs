//! RUN-CONFIG (2026-11-18): `[performance]` and `[profile]` tables of the input
//! file, precedence environment > input > default, strict input errors, the
//! extended switch manifest (requested, source, consulted, path, fallback) and
//! the `run` block of results.json (run ID, exit reason, binary SHA-256,
//! completion marker). Every run is its own `dft-cli` process.
#[allow(dead_code)]
#[path = "../common/mod.rs"]
mod common;
use serde_json::Value;
use sha2::{Digest, Sha256};
use std::{fs, path::Path, process::Command};

const BASE: &str = "[system]\nname='synthetic run config'\nlattice_unit='bohr'\nlattice_vectors=[[8.0,0.0,0.0],[0.0,8.0,0.0],[0.0,0.0,8.0]]\n[[system.atoms]]\nspecies='H'\nfractional=[0.0,0.0,0.0]\n[[system.atoms]]\nspecies='H'\nfractional=[0.27,0.23,0.19]\n[pseudopotentials]\nH='syn.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[scf]\nnbands=4\n";

struct Run {
    code: Option<i32>,
    stderr: String,
    result: Option<Value>,
    density_sha: Option<String>,
}

fn run(dir: &Path, name: &str, tables: &str, env: &[(&str, &str)]) -> Run {
    let text = format!("{BASE}[output]\ndirectory='{name}'\n{tables}");
    let input = format!("{name}.toml");
    fs::write(dir.join(&input), text).unwrap();
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_dft-cli"));
    cmd.current_dir(dir).args(["run", &input]);
    for name in [
        "BAANDI_STOP_METRIC",
        "BAANDI_PHASE_PROFILE",
        "BAANDI_M2_GRID",
    ] {
        cmd.env_remove(name);
    }
    for (k, v) in env {
        cmd.env(k, v);
    }
    let out = cmd.output().unwrap();
    let read = |f: &str| fs::read(dir.join(name).join(f)).ok();
    Run {
        code: out.status.code(),
        stderr: String::from_utf8_lossy(&out.stderr).into_owned(),
        result: read("results.json").map(|b| serde_json::from_slice(&b).unwrap()),
        density_sha: read("density.bin").map(|b| format!("{:x}", Sha256::digest(b))),
    }
}

fn ok(r: &Run) -> &Value {
    assert_eq!(r.code, Some(0), "{}", r.stderr);
    r.result.as_ref().unwrap()
}

#[test]
fn input_tables_precedence_manifest_and_run_block() {
    let dir = std::env::temp_dir().join(format!("run-config-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("syn.xml"),
        common::gaussian_with(1., 1.5, false, 401).replace("element=\"Xx\"", "element=\"H\""),
    )
    .unwrap();

    // Default: no tables.
    let base = run(&dir, "base", "", &[]);
    let b = ok(&base);
    let r = &b["run"];
    assert_eq!(r["schema"], "run-report-1");
    assert_eq!(r["completed"], true);
    assert_eq!(r["exit_reason"], "converged");
    assert_eq!(r["run_id"].as_str().unwrap().len(), 16);
    let exe = fs::read(env!("CARGO_BIN_EXE_dft-cli")).unwrap();
    assert_eq!(
        r["binary"]["sha256"].as_str().unwrap(),
        format!("{:x}", Sha256::digest(&exe))
    );
    assert!(r["counters"]["session_steps"].as_u64().unwrap() > 0);
    let sw = &b["scf"]["switches"]["BAANDI_STOP_METRIC"];
    assert_eq!(sw["source"], "default");
    assert_eq!(sw["consulted"], true);
    assert_eq!(sw["requested"], Value::Null);
    assert_eq!(sw["path"], "kept");
    // MPI-only switch: not consulted by a one-rank run.
    assert_eq!(b["scf"]["switches"]["BAANDI_M2_STICKS"]["consulted"], false);

    // Explicit defaults in the input: bit-identical science output.
    let same = run(
        &dir,
        "explicit_defaults",
        "[performance]\nstop_metric='kept'\nm2_grid=true\n",
        &[],
    );
    let s = ok(&same);
    assert_eq!(same.density_sha, base.density_sha);
    assert_eq!(s["total_energy"]["ha"], b["total_energy"]["ha"]);
    assert_eq!(
        s["scf"]["switches"]["BAANDI_STOP_METRIC"]["source"],
        "input"
    );
    assert_eq!(
        s["run"]["switch_agreement_hash"], b["run"]["switch_agreement_hash"],
        "same effective values, same hash"
    );
    assert_ne!(s["run"]["run_id"], b["run"]["run_id"]);

    // Input value acts; profile table; path recorded.
    let fft = run(
        &dir,
        "fft",
        "[performance]\nstop_metric='fft'\n[profile]\nphase_profile=true\n",
        &[],
    );
    let f = ok(&fft);
    let sw = &f["scf"]["switches"]["BAANDI_STOP_METRIC"];
    assert_eq!(sw["value"], "fft");
    assert_eq!(sw["source"], "input");
    assert_eq!(sw["requested"], "fft");
    assert_eq!(sw["path"], "fft");
    assert_eq!(
        f["scf"]["switches"]["BAANDI_PHASE_PROFILE"]["source"],
        "input"
    );
    assert!(f.get("phase_profile").is_some(), "profile table acts");
    assert_eq!(f["run"]["counters"]["stop_kept"], 0);
    assert!(f["run"]["counters"]["stop_fft"].as_u64().unwrap() > 0);
    assert_ne!(
        f["run"]["switch_agreement_hash"],
        b["run"]["switch_agreement_hash"]
    );
    let de = f["total_energy"]["ha"].as_f64().unwrap() - b["total_energy"]["ha"].as_f64().unwrap();
    assert!(de.abs() <= 5e-9, "{de:e}");

    // scripts/parse_results.py on real output: accepted with the binary that
    // wrote it; rejected when truncated (interrupted write) or stale.
    let parser = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../scripts/parse_results.py");
    let results = dir.join("fft/results.json");
    let parse = |file: &Path, extra: &[&str]| {
        let out = Command::new("python3")
            .arg("-I")
            .arg(&parser)
            .arg(file)
            .args(extra)
            .output()
            .unwrap();
        (
            out.status.code(),
            String::from_utf8_lossy(&out.stdout).into_owned(),
        )
    };
    let (code, text) = parse(&results, &["--binary", env!("CARGO_BIN_EXE_dft-cli")]);
    assert_eq!(code, Some(0), "{text}");
    assert!(
        text.starts_with("ACCEPT") && text.contains("imbalance="),
        "{text}"
    );
    let full = fs::read(&results).unwrap();
    let cut = dir.join("truncated.json");
    fs::write(&cut, &full[..full.len() * 2 / 3]).unwrap();
    let (code, text) = parse(&cut, &[]);
    assert_eq!(code, Some(2));
    assert!(text.contains("REJECT: truncated"), "{text}");
    let (code, text) = parse(&results, &["--expect-sha", &"0".repeat(64)]);
    assert_eq!(code, Some(2));
    assert!(text.contains("REJECT: stale binary"), "{text}");
    let closure = &f["phase_profile"]["closure"];
    assert!(
        closure["residual_min_s"].as_f64().unwrap() >= -1e-3,
        "{closure}"
    );

    // Environment over input.
    let env = run(
        &dir,
        "env",
        "[performance]\nstop_metric='fft'\n",
        &[("BAANDI_STOP_METRIC", "kept")],
    );
    let e = ok(&env);
    let sw = &e["scf"]["switches"]["BAANDI_STOP_METRIC"];
    assert_eq!(
        (sw["value"].as_str(), sw["source"].as_str()),
        (Some("kept"), Some("env"))
    );
    assert_eq!(env.density_sha, base.density_sha);

    // Strict input: errors before any work, exit 1, no results file.
    for (name, tables, message) in [
        (
            "bad_value",
            "[performance]\nstop_metric='maybe'\n",
            "allowed values fft, kept",
        ),
        (
            "bad_key",
            "[performance]\nno_such_switch='on'\n",
            "unknown switch",
        ),
        (
            "wrong_table",
            "[performance]\nphase_profile='on'\n",
            "belongs to the [profile] table",
        ),
        (
            "env_only",
            "[performance]\npseudo_dir='/tmp'\n",
            "environment only",
        ),
    ] {
        let r = run(&dir, name, tables, &[]);
        assert_eq!(r.code, Some(1), "{name}: {}", r.stderr);
        assert!(r.stderr.contains(message), "{name}: {}", r.stderr);
        assert!(r.result.is_none(), "{name}");
    }
    let _ = fs::remove_dir_all(&dir);
}
