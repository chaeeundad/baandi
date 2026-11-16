//! ADR 0043 amendment ONEPATH-4 (user decision 2026-11-16): the last strict
//! fallbacks are gone. The observe/agent runtime modes, `bands`/`dos`,
//! `checkpoint`/`restart` run the single standard path with the fixed-point
//! verification; `scf.locking = "auto"` is an input error on the default path
//! (tests/single_path.rs); strict stays reachable only through
//! `BAANDI_REFERENCE_STRICT=on`.
#[allow(dead_code)]
mod common;
use dft_cli::checkpoint_io as io;
use serde_json::Value;
use sha2::{Digest, Sha256};
use std::{
    fs,
    path::{Path, PathBuf},
    process::{Command, Output},
    sync::atomic::{AtomicUsize, Ordering},
};

static SERIAL: AtomicUsize = AtomicUsize::new(0);
const BASE: &str = "[system]\nname='onepath4'\nlattice_unit='bohr'\nlattice_vectors=[[8.0,0.0,0.0],[0.0,8.0,0.0],[0.0,0.0,8.0]]\n[[system.atoms]]\nspecies='H'\nfractional=[0.0,0.0,0.0]\n[[system.atoms]]\nspecies='H'\nfractional=[0.27,0.23,0.19]\n[pseudopotentials]\nH='syn.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n";

struct Case(PathBuf);
impl Case {
    fn new(scf: &str) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "onepath4-{}-{}",
            std::process::id(),
            SERIAL.fetch_add(1, Ordering::SeqCst)
        ));
        fs::create_dir_all(&dir).unwrap();
        fs::write(
            dir.join("syn.xml"),
            common::gaussian_with(1., 1.5, false, 401).replace("element=\"Xx\"", "element=\"H\""),
        )
        .unwrap();
        fs::write(
            dir.join("input.toml"),
            format!("{BASE}[scf]\nnbands=4\n{scf}"),
        )
        .unwrap();
        Self(dir)
    }
    /// The default path: the reference switch is removed from the child.
    fn cli(&self, args: &[&str], env: &[(&str, &str)]) -> Output {
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_dft-cli"));
        cmd.current_dir(&self.0)
            .env_remove("BAANDI_REFERENCE_STRICT")
            .args(args);
        for (k, v) in env {
            cmd.env(k, v);
        }
        cmd.output().unwrap()
    }
    fn json(&self, path: &str) -> Value {
        serde_json::from_slice(&fs::read(self.0.join(path)).unwrap()).unwrap()
    }
    fn sha(&self, path: &str) -> String {
        format!("{:x}", Sha256::digest(fs::read(self.0.join(path)).unwrap()))
    }
}
impl Drop for Case {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn ok(o: &Output) {
    assert!(
        o.status.success(),
        "{}\n{}",
        String::from_utf8_lossy(&o.stdout),
        String::from_utf8_lossy(&o.stderr)
    );
}
fn verified_standard(v: &Value) {
    assert_eq!(v["single_path"]["status"], "applied");
    assert_eq!(v["single_path"]["effective_profile"], "standard");
    assert_eq!(v["single_path"]["verification"], "fixed_point");
    assert_eq!(v["convergence"]["converged"], true);
    assert_eq!(v["verified"], true);
    let attempts = v["verification"]["attempts"].as_array().unwrap();
    assert_eq!(attempts.last().unwrap()["passed"], true);
}

/// `run --observe` runs the same standard path with the verification as a
/// classic run (cutoff ladder off: the runtime modes have no coarse stage):
/// same energy bits, density, attempts; the event log reports the passed
/// verification as `certificate_evaluated`.
#[test]
fn observe_mode_runs_the_verified_standard_path() {
    let c = Case::new("cutoff_ladder='off'\n");
    ok(&c.cli(&["run", "input.toml", "--output-dir", "classic"], &[]));
    ok(&c.cli(
        &[
            "run",
            "input.toml",
            "--output-dir",
            "observed",
            "--observe",
            "--events",
            "events.jsonl",
        ],
        &[],
    ));
    let classic = c.json("classic/results.json");
    let observed = c.json("observed/results.json");
    verified_standard(&classic);
    verified_standard(&observed);
    assert_eq!(
        classic["total_energy"]["ha"].as_f64().unwrap().to_bits(),
        observed["total_energy"]["ha"].as_f64().unwrap().to_bits()
    );
    assert_eq!(c.sha("classic/density.bin"), c.sha("observed/density.bin"));
    assert_eq!(
        classic["verification"]["attempts"]
            .as_array()
            .unwrap()
            .len(),
        observed["verification"]["attempts"]
            .as_array()
            .unwrap()
            .len()
    );
    let events = fs::read_to_string(c.0.join("events.jsonl")).unwrap();
    let last_certificate = events
        .lines()
        .map(|l| serde_json::from_str::<Value>(l).unwrap())
        .rfind(|e| e["type"] == "certificate_evaluated")
        .expect("verification reported as certificate_evaluated");
    assert_eq!(last_certificate["data"]["passed"], true);
    // The reference switch keeps the strict path in the runtime modes too.
    ok(&c.cli(
        &["run", "input.toml", "--output-dir", "strict", "--observe"],
        &[("BAANDI_REFERENCE_STRICT", "on")],
    ));
    let strict = c.json("strict/results.json");
    assert_eq!(strict["single_path"]["status"], "reference_strict");
    assert!(strict.get("verification").is_none());
    let de = (strict["total_energy"]["ha"].as_f64().unwrap()
        - observed["total_energy"]["ha"].as_f64().unwrap())
    .abs();
    println!("observe: standard vs strict reference |dE| = {de:e} Ha");
    assert!(de <= 5e-9, "{de:e} Ha");
}

/// `run --agent` (Unix socket binding; ignored like `run.rs::agent_socket_mode`):
/// `cargo test -p dft-cli --test onepath_last_fallbacks agent_mode -- --ignored --nocapture`.
#[test]
#[ignore = "requires filesystem Unix socket binding"]
fn agent_mode_runs_the_verified_standard_path() {
    let c = Case::new("cutoff_ladder='off'\n");
    let out = c.0.join("agent");
    fs::create_dir(&out).unwrap();
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        fs::set_permissions(&out, fs::Permissions::from_mode(0o700)).unwrap();
    }
    ok(&c.cli(
        &["run", "input.toml", "--output-dir", "agent", "--agent"],
        &[],
    ));
    verified_standard(&c.json("agent/results.json"));
}

/// `bands`/`dos` run their SCF on the verified standard path and record it.
#[test]
fn bands_and_dos_scf_run_the_verified_standard_path() {
    // The spectrum fixture of tests/bands_dos.rs (one synthetic Si-labelled ion).
    let c = Case::new("");
    fs::write(c.0.join("synthetic.xml"), common::fixture()).unwrap();
    fs::write(c.0.join("input.toml"), "[system]\nname='synthetic spectrum'\nlattice_unit='bohr'\nlattice_vectors=[[4.0,0.0,0.0],[0.0,4.0,0.0],[0.0,0.0,4.0]]\n[[system.atoms]]\nspecies='Si'\nfractional=[0.13,0.21,0.32]\n[pseudopotentials]\nSi='synthetic.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[bands]\npoints_per_segment=2\npoints=[{label='Gamma',fractional=[0.0,0.0,0.0]},{label='A',fractional=[0.1,0.0,0.0]},{label='B',fractional=[0.0,0.1,0.0],break_before=true},{label='Gamma',fractional=[0.0,0.0,0.0]}]\n[dos]\ndos_kgrid=[1,1,1]\n").unwrap();
    for (command, file) in [("bands", "bands.json"), ("dos", "dos.json")] {
        let o = c.cli(&[command, "input.toml", "--output-dir", "out"], &[]);
        ok(&o);
        let stdout = String::from_utf8_lossy(&o.stdout);
        assert!(
            stdout.contains("SCF path: standard") && stdout.contains("passed: true"),
            "{stdout}"
        );
        if command == "bands" {
            let meta = c.json(&format!("out/{file}"));
            assert_eq!(meta["scf_path"]["effective_profile"], "standard");
            assert_eq!(meta["scf_path"]["verified"], true);
            assert_eq!(meta["scf_path"]["certified"], false);
        }
    }
    let o = c.cli(
        &["bands", "input.toml", "--output-dir", "strict"],
        &[("BAANDI_REFERENCE_STRICT", "on")],
    );
    ok(&o);
    assert_eq!(
        c.json("strict/bands.json")["scf_path"]["effective_profile"],
        "strict"
    );
}

fn stripped(path: &Path) -> String {
    let mut cp = io::read(path).unwrap();
    let r = cp.result.as_mut().unwrap();
    for h in &mut r.history {
        h.timings = Default::default();
        if let Some(c) = &mut h.fixed_point_certificate {
            c.timings = Default::default();
        }
    }
    for a in &mut r.verification_attempts {
        a.wall = std::time::Duration::ZERO;
    }
    format!("{:?}|{:?}", cp.electronic_state, cp.result)
}

/// `checkpoint`/`restart` carry the verification. With REUSE_PSI off an exact
/// restart reproduces the uninterrupted run bit for bit (state, energies,
/// history, verification attempts without wall time). With REUSE_PSI on (the
/// default) the checkpointed session forgets its kept state, so the restart
/// matches the uninterrupted run within the REUSE_PSI rounding class
/// (pre-registered bound 1e-10 Ha, the bound of the REUSE_PSI adoption tests).
#[test]
fn checkpoint_restart_carries_the_verification() {
    let c = Case::new("");
    let off = [("BAANDI_REUSE_PSI", "off")];
    let exact = c.cli(
        &["checkpoint", "input.toml", "exact.cp", "--iteration", "2"],
        &off,
    );
    ok(&exact);
    assert!(
        String::from_utf8_lossy(&exact.stdout).contains("verification: fixed_point"),
        "{}",
        String::from_utf8_lossy(&exact.stdout)
    );
    let saved = io::read(&c.0.join("exact.cp")).unwrap();
    assert_eq!(saved.verification, dft_engine::Verification::FixedPoint);
    assert_eq!(
        saved.continuation.as_ref().unwrap().options.verification,
        dft_engine::Verification::FixedPoint
    );
    ok(&c.cli(&["restart", "input.toml", "exact.cp", "restarted.cp"], &off));
    ok(&c.cli(
        &[
            "checkpoint",
            "input.toml",
            "baseline.cp",
            "--kind",
            "result",
        ],
        &off,
    ));
    let restarted = io::read(&c.0.join("restarted.cp")).unwrap();
    let r = restarted.result.as_ref().unwrap();
    assert!(r.converged);
    assert!(r.verification_attempts.last().is_some_and(|a| a.passed));
    assert_eq!(
        stripped(&c.0.join("baseline.cp")),
        stripped(&c.0.join("restarted.cp"))
    );
    // Inspection shows the carried option and attempts.
    let status = c.cli(&["status", "restarted.cp"], &[]);
    ok(&status);
    let v: Value = serde_json::from_slice(&status.stdout).unwrap();
    assert_eq!(v["verification"], "fixed_point");
    assert_eq!(v["verified"], true);
    // REUSE_PSI on (default): rounding class against the uninterrupted run.
    ok(&c.cli(
        &[
            "checkpoint",
            "input.toml",
            "exact-reuse.cp",
            "--iteration",
            "2",
        ],
        &[],
    ));
    ok(&c.cli(
        &[
            "restart",
            "input.toml",
            "exact-reuse.cp",
            "restarted-reuse.cp",
        ],
        &[],
    ));
    ok(&c.cli(
        &[
            "checkpoint",
            "input.toml",
            "baseline-reuse.cp",
            "--kind",
            "result",
        ],
        &[],
    ));
    let energy = |p: &str| {
        let cp = io::read(&c.0.join(p)).unwrap();
        let r = cp.result.unwrap();
        assert!(r.converged && r.verification_attempts.last().unwrap().passed);
        r.energy.total_hartree()
    };
    let de = (energy("restarted-reuse.cp") - energy("baseline-reuse.cp")).abs();
    println!("checkpoint/restart with REUSE_PSI on: |dE| against uninterrupted = {de:e} Ha");
    assert!(de <= 1e-10, "{de:e} Ha");
}

/// Wire: the verification option and attempts round-trip bit for bit; a file
/// without them (verification off) carries no new key, so files written
/// before ONEPATH-4 decode unchanged.
#[test]
fn v1_codec_round_trips_verification_and_omits_it_when_off() {
    let c = Case::new("max_iterations=3\n");
    ok(&c.cli(
        &["checkpoint", "input.toml", "plain.cp", "--iteration", "1"],
        &[("BAANDI_REFERENCE_STRICT", "on")],
    ));
    let mut cp = io::read(&c.0.join("plain.cp")).unwrap();
    let bytes = io::encode(&cp).unwrap();
    let header = String::from_utf8_lossy(&bytes[24..]).into_owned();
    assert!(!header.contains("verification"), "strict file gained a key");
    cp.verification = dft_engine::Verification::FixedPoint;
    let cont = cp.continuation.as_mut().unwrap();
    cont.options.verification = dft_engine::Verification::FixedPoint;
    cont.verification_attempts
        .push(dft_engine::VerificationAttempt {
            iteration: 1,
            passed: false,
            coulomb_error_ha: 3.25e-9,
            threshold_ha: 5e-9,
            eigen_tolerance_ha: 1e-6,
            residual_bound_ha: Some(1.5e-10),
            density_residual: 2e-5,
            occupied_eigenvalue_change_ha: 4e-8,
            strict_criteria_met: false,
            davidson_iterations: 3,
            operator_applications: 17,
            wall: std::time::Duration::new(0, 123_456),
            spin_density_residual: None,
            occupation_change_electrons: None,
            chemical_potential_change_ha: None,
            moment_change: None,
        });
    let bytes = io::encode(&cp).unwrap();
    let back = io::decode(&bytes).unwrap();
    assert_eq!(io::encode(&back).unwrap(), bytes);
    assert_eq!(back.verification, dft_engine::Verification::FixedPoint);
    let attempts = &back.continuation.as_ref().unwrap().verification_attempts;
    assert_eq!(
        format!("{attempts:?}"),
        format!(
            "{:?}",
            cp.continuation.as_ref().unwrap().verification_attempts
        )
    );
}
