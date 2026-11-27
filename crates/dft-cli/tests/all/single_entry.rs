//! UX1 (ADR 0053): the bare `baandi` entry, the progress stream, the job lock
//! and the plain monitor view, through the process boundary.
#[allow(dead_code)]
#[path = "../common/mod.rs"]
mod common;
use serde_json::Value;
use sha2::{Digest, Sha256};
use std::{
    fs,
    path::{Path, PathBuf},
    process::{Command, Output},
};

const SCF: &str = "[system]\nname='synthetic small'\nlattice_unit='bohr'\nlattice_vectors=[[4.0,0.0,0.0],[0.0,4.0,0.0],[0.0,0.0,4.0]]\n[[system.atoms]]\nspecies='Si'\nfractional=[0.13,0.21,0.32]\n[pseudopotentials]\nSi='synthetic.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n";

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("baandi-ux1-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("synthetic.xml"),
        common::fixture().replace("element=\"Xx\"", "element=\"Si\""),
    )
    .unwrap();
    dir
}
fn baandi(dir: &Path) -> Command {
    let mut c = Command::new(env!("CARGO_BIN_EXE_baandi"));
    c.current_dir(dir).env("BAANDI_BINARY_HASH_CACHE", "off");
    c
}
fn ok(o: &Output, code: i32) {
    assert_eq!(
        o.status.code(),
        Some(code),
        "stdout={}\nstderr={}",
        String::from_utf8_lossy(&o.stdout),
        String::from_utf8_lossy(&o.stderr)
    );
}
fn events(dir: &Path) -> Vec<Value> {
    fs::read_to_string(dir.join("progress.jsonl"))
        .unwrap()
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect()
}
fn hash(p: &Path) -> String {
    format!("{:x}", Sha256::digest(fs::read(p).unwrap()))
}
fn energy(dir: &Path) -> u64 {
    let v: Value = serde_json::from_slice(&fs::read(dir.join("results.json")).unwrap()).unwrap();
    v["energies"]["total"]["ha"]
        .as_f64()
        .unwrap_or_else(|| panic!("no total energy in {v}"))
        .to_bits()
}

/// Bare `baandi` runs `[task] type = "scf"` like `baandi run`; the stream
/// does not change the numbers (density bytes and total energy bits equal
/// with and without it); the alias writes a stream too; the plain monitor
/// shows the finished job.
#[test]
fn bare_entry_runs_the_task_and_the_stream_leaves_results_identical() {
    let dir = scratch("scf");
    fs::write(dir.join("input.toml"), format!("[task]\ntype='scf'\n{SCF}")).unwrap();
    let bare = baandi(&dir).output().unwrap();
    ok(&bare, 0);
    let out = dir.join("results");
    assert!(!out.join("baandi.lock").exists(), "lock released");
    let e = events(&out);
    let kinds: Vec<&str> = e.iter().map(|v| v["event"].as_str().unwrap()).collect();
    assert_eq!(kinds.first(), Some(&"start"));
    assert_eq!(kinds.last(), Some(&"end"));
    for k in ["scf_start", "scf_map", "scf_end"] {
        assert!(kinds.contains(&k), "{kinds:?}");
    }
    assert_eq!(e[0]["schema"], "baandi-progress-1");
    assert_eq!(e[0]["task"], "scf");
    assert_eq!(e[0]["system"]["formula"], "Si");
    assert_eq!(e.last().unwrap()["status"], "completed");
    assert!(e.iter().any(|v| v["kind"] == "verification"));
    for (i, v) in e.iter().enumerate() {
        assert_eq!(v["seq"], i as u64 + 1);
    }

    let silent = baandi(&dir)
        .args([
            "run",
            "input.toml",
            "--output-dir",
            "silent",
            "--no-progress",
        ])
        .output()
        .unwrap();
    ok(&silent, 0);
    assert!(!dir.join("silent/progress.jsonl").exists());
    let alias = baandi(&dir)
        .args(["run", "input.toml", "--output-dir", "alias"])
        .output()
        .unwrap();
    ok(&alias, 0);
    assert_eq!(
        events(&dir.join("alias")).last().unwrap()["status"],
        "completed"
    );
    for other in ["silent", "alias"] {
        assert_eq!(
            hash(&out.join("density.bin")),
            hash(&dir.join(other).join("density.bin")),
            "{other}"
        );
        assert_eq!(energy(&out), energy(&dir.join(other)), "{other}");
    }

    let view = baandi(&dir)
        .args(["watch", "--plain", "--once"])
        .output()
        .unwrap();
    ok(&view, 0);
    let text = String::from_utf8(view.stdout).unwrap();
    assert!(text.contains("[COMPLETED (exit 0)]"), "{text}");
    assert!(
        text.contains("Task: scf") && text.contains("SCF 1 of 1"),
        "{text}"
    );
    // From inside the output directory a bare `baandi` (no input.toml there)
    // shows the finished job instead of failing.
    let inside = baandi(&out).arg("--plain").output().unwrap();
    ok(&inside, 0);
    assert!(String::from_utf8_lossy(&inside.stdout).contains("COMPLETED"));
    fs::remove_dir_all(&dir).unwrap();
}

/// Input and lock errors of the single entry, and the watch of a live job.
#[test]
fn bare_entry_rejects_ambiguous_input_and_watches_a_running_job() {
    let dir = scratch("errors");
    let o = baandi(&dir).output().unwrap();
    ok(&o, 1);
    assert!(String::from_utf8_lossy(&o.stderr).contains("no input.toml"));
    fs::write(
        dir.join("input.toml"),
        format!("{SCF}[relax]\nmax_steps=2\n"),
    )
    .unwrap();
    let o = baandi(&dir).output().unwrap();
    ok(&o, 1);
    assert!(String::from_utf8_lossy(&o.stderr).contains("does not select relaxation"));
    let o = baandi(&dir)
        .args(["-i", "input.toml", "grid", "input.toml"])
        .output()
        .unwrap();
    assert!(!o.status.success());
    fs::write(
        dir.join("job.toml"),
        format!("[task]\ntype='dos'\n{SCF}[output]\ndirectory='busy'\n"),
    )
    .unwrap();

    // A live lock (a running `sleep`) refuses a second job in the directory and
    // turns a bare `baandi -i` into the monitor, which returns once the
    // process is gone.
    let mut sleeper = Command::new("sleep").arg("3").spawn().unwrap();
    let busy = dir.join("busy");
    fs::create_dir_all(&busy).unwrap();
    let host = String::from_utf8(Command::new("hostname").output().unwrap().stdout).unwrap();
    fs::write(
        busy.join("baandi.lock"),
        serde_json::json!({"schema":"baandi-lock-1","pid":sleeper.id(),"host":host.trim(),
            "started_unix":0.0,"task":"dos","input":"job.toml"})
        .to_string(),
    )
    .unwrap();
    fs::write(
        busy.join("progress.jsonl"),
        "{\"v\":1,\"seq\":1,\"t\":0.0,\"elapsed_s\":0.0,\"event\":\"start\",\"task\":\"dos\",\"input\":\"job.toml\",\"system\":null}\n",
    )
    .unwrap();
    let refused = baandi(&dir).args(["dos", "job.toml"]).output().unwrap();
    ok(&refused, 5);
    assert!(String::from_utf8_lossy(&refused.stderr).contains("is running"));
    let reaper = std::thread::spawn(move || sleeper.wait().unwrap());
    let watched = baandi(&dir).args(["-i", "job.toml"]).output().unwrap();
    reaper.join().unwrap();
    ok(&watched, 0);
    let text = String::from_utf8_lossy(&watched.stdout);
    assert!(text.contains("RUNNING pid"), "{text}");
    assert!(
        text.contains("STOPPED"),
        "a stopped job is reported: {text}"
    );
    fs::remove_dir_all(&dir).unwrap();
}
