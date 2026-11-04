#[allow(dead_code)]
#[path = "../../dft-engine/tests/common/mod.rs"]
mod common;
#[allow(dead_code)]
#[path = "../../dft-engine/tests/common/checkpoint.rs"]
mod fixture;
use dft_cli::checkpoint_io;
use dft_engine::{checkpoint::*, *};
use serde_json::Value;
use std::{
    fs,
    path::Path,
    process::Command,
    sync::atomic::{AtomicUsize, Ordering},
};

static SERIAL: AtomicUsize = AtomicUsize::new(0);
fn inspect(command: &str, path: &Path) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_dft-cli"))
        .arg(command)
        .arg(path)
        .output()
        .unwrap()
}
fn check(cp: &Checkpoint, state: &str, iteration: Option<usize>) {
    let dir = std::env::temp_dir().join(format!(
        "p11h-inspect-{}-{}",
        std::process::id(),
        SERIAL.fetch_add(1, Ordering::SeqCst)
    ));
    fs::create_dir(&dir).unwrap();
    // Content detection must work even with a TOML extension.
    let path = dir.join("saved.toml");
    checkpoint_io::write(&path, cp).unwrap();
    let before = fs::read(&path).unwrap();
    let output = inspect("status", &path);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let v: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(v["kind"], format!("{:?}", cp.kind));
    assert_eq!(v["state"], state);
    assert_eq!(v["iteration"], serde_json::json!(iteration));
    assert_eq!(v["calc_id"], cp.calc_id);
    assert_eq!(
        v["converged"],
        serde_json::json!(cp.result.as_ref().map(|r| r.converged))
    );
    if let Some(r) = &cp.result {
        // Printed JSON is for inspection; compare with the existing scalar display
        // tolerance, not a restart guarantee (the container itself preserves bits).
        let shown = v["result"]["energy_hartree"].as_f64().unwrap();
        assert!((shown - r.energy.total_hartree()).abs() < 1e-14);
    }
    if let Some(c) = &cp.continuation {
        assert_eq!(v["state_version"], c.state_version);
        assert_eq!(v["last_iteration"].is_null(), c.history.is_empty());
    }
    let explanation = inspect("explain", &path);
    assert!(
        explanation.status.success(),
        "{}",
        String::from_utf8_lossy(&explanation.stderr)
    );
    let text = String::from_utf8(explanation.stdout).unwrap();
    assert!(
        text.contains("Hypotheses: none inferred")
            && text.contains("restart compatibility")
            && text.contains("No SCF executed")
    );
    assert_eq!(before, fs::read(&path).unwrap());
    assert_eq!(fs::read_dir(&dir).unwrap().count(), 1);
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn all_checkpoint_kinds_report_recorded_state_without_running_or_writing() {
    let mut session = ScfSession::new(
        fixture::fixture(0, Mixing::default(), 1, 1),
        CompletionMode::RetainLast,
        ExecutionConfig::default(),
    )
    .unwrap();
    session.initialize().unwrap();
    let capture = |s: &ScfSession, kind| {
        let (m, e) = fixture::identity(0, 1);
        s.checkpoint(kind, m, e).unwrap()
    };
    check(
        &capture(&session, CheckpointKind::ExactContinuation),
        "nonterminal: process liveness unknown",
        Some(0),
    );
    check(
        &capture(&session, CheckpointKind::WarmStart),
        "warm_start: convergence not recorded",
        None,
    );
    session.step().unwrap();
    let cp = capture(&session, CheckpointKind::ResultExport);
    assert!(!cp.result.as_ref().unwrap().converged);
    check(&cp, "not_converged", Some(1));

    let mut session = ScfSession::new(
        fixture::fixture(0, Mixing::default(), 1, 200),
        CompletionMode::RetainLast,
        ExecutionConfig::default(),
    )
    .unwrap();
    session.initialize().unwrap();
    session.step().unwrap();
    check(
        &capture(&session, CheckpointKind::ExactContinuation),
        "nonterminal: process liveness unknown",
        Some(1),
    );
    while session.step().unwrap().status == StepStatus::Ready {}
    let cp = capture(&session, CheckpointKind::ResultExport);
    assert!(cp.result.as_ref().unwrap().converged);
    let iterations = cp.result.as_ref().unwrap().history.len();
    check(&cp, "converged", Some(iterations));
    println!(
        "P11H checkpoint inspection: three kinds, boundaries 0/1, unconverged and converged exports; scalar energy agreement < 1e-14 Ha; files unchanged."
    );
}

#[test]
fn damaged_and_inconsistent_checkpoints_fail_closed() {
    let (m, e) = fixture::identity(0, 1);
    let mut s = ScfSession::new(
        fixture::fixture(0, Mixing::default(), 1, 1),
        CompletionMode::RetainLast,
        ExecutionConfig::default(),
    )
    .unwrap();
    s.initialize().unwrap();
    let cp = s
        .checkpoint(CheckpointKind::ExactContinuation, m, e)
        .unwrap();
    let bytes = checkpoint_io::encode(&cp).unwrap();
    let mut checksum = bytes.clone();
    *checksum.last_mut().unwrap() ^= 1;
    let mut schema = cp.clone();
    schema.schema_version = 2;
    let mut missing = cp.clone();
    missing.continuation = None;
    let mut wrong_kind = cp;
    wrong_kind.kind = CheckpointKind::WarmStart;
    let path = std::env::temp_dir().join(format!("p11h-bad-{}.cp", std::process::id()));
    for bytes in [
        checksum,
        bytes[..20].to_vec(),
        checkpoint_io::encode(&schema).unwrap(),
        checkpoint_io::encode(&missing).unwrap(),
        checkpoint_io::encode(&wrong_kind).unwrap(),
        b"wrong magic".to_vec(),
    ] {
        fs::write(&path, &bytes).unwrap();
        for command in ["status", "explain"] {
            let output = inspect(command, &path);
            assert!(!output.status.success());
            assert!(output.stdout.is_empty());
            assert!(!output.stderr.is_empty());
        }
        assert_eq!(fs::read(&path).unwrap(), bytes);
    }
    fs::remove_file(path).unwrap();
}
