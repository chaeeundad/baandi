//! ADR 0043 single SCF path, the default since the user decision of
//! 2026-11-16 (`scf.convergence_profile` removed). The strict path stays as a
//! reference tool behind `BAANDI_REFERENCE_STRICT=on`.
#[allow(dead_code)]
mod common;
use serde_json::Value;
use sha2::{Digest, Sha256};
use std::{fs, path::PathBuf, process::Command};

const BASE: &str = "[system]\nname='synthetic single path'\nlattice_unit='bohr'\nlattice_vectors=[[8.0,0.0,0.0],[0.0,8.0,0.0],[0.0,0.0,8.0]]\n[[system.atoms]]\nspecies='H'\nfractional=[0.0,0.0,0.0]\n[[system.atoms]]\nspecies='H'\nfractional=[0.27,0.23,0.19]\n[pseudopotentials]\nH='syn.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n";

struct Run {
    code: Option<i32>,
    stderr: String,
    result: Option<Value>,
    density_sha: Option<String>,
}

/// `reference`: value of `BAANDI_REFERENCE_STRICT` (None = unset).
fn run(dir: &PathBuf, name: &str, reference: Option<&str>, scf: &str, output: &str) -> Run {
    let text = format!("{BASE}[scf]\nnbands=4\n{scf}[output]\ndirectory='{name}'\n{output}");
    let input = format!("{name}.toml");
    fs::write(dir.join(&input), text).unwrap();
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_dft-cli"));
    cmd.current_dir(dir).args(["run", &input]);
    match reference {
        Some(v) => cmd.env("BAANDI_REFERENCE_STRICT", v),
        None => cmd.env_remove("BAANDI_REFERENCE_STRICT"),
    };
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

fn energy_bits(v: &Value) -> u64 {
    v["total_energy"]["ha"].as_f64().unwrap().to_bits()
}

#[test]
fn single_path_default_contract() {
    let dir = std::env::temp_dir().join(format!("adr43-single-path-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("syn.xml"),
        common::gaussian_with(1., 1.5, false, 401).replace("element=\"Xx\"", "element=\"H\""),
    )
    .unwrap();

    // Reference switch on: the strict path with its certificate (same input).
    let strict = run(&dir, "strict", Some("on"), "", "forces=true\n");
    let u = ok(&strict);
    assert!(u.get("verification").is_none() && u.get("convergence_profile").is_none());
    assert!(u["timings"].get("verification_s").is_none());
    assert_eq!(u["single_path"]["status"], "reference_strict");
    assert_eq!(u["single_path"]["effective_profile"], "strict");
    assert_eq!(u["single_path"]["verification"], "strict_certificate");
    assert_eq!(u["execution"]["effective_profile"], "strict");
    assert_eq!(u["execution"]["requested"]["reference_strict"], true);
    assert_eq!(
        u["execution"]["non_default_switches"]["BAANDI_REFERENCE_STRICT"]["value"],
        "on"
    );
    assert_eq!(
        u["scf"]["switches"]["BAANDI_REFERENCE_STRICT"]["value"],
        "on"
    );
    assert_eq!(
        u["forces"]["preconditions"]["fixed_point_certificate_passed"],
        true
    );

    // Default (switch unset and explicit off): standard path, verification, forces.
    let on = run(&dir, "on", None, "", "forces=true\n");
    let off = run(&dir, "off", Some("off"), "", "forces=true\n");
    let v = ok(&on);
    let o = ok(&off);
    assert_eq!(on.density_sha, off.density_sha);
    assert_eq!(energy_bits(v), energy_bits(o));
    assert_eq!(v["single_path"]["status"], "applied");
    assert_eq!(v["single_path"]["verification"], "fixed_point");
    assert_eq!(v["convergence_profile"], "standard");
    assert_eq!(v["standard_threshold_ha"], 5e-9);
    assert_eq!(v["certified"], false);
    assert_eq!(v["verified"], true);
    assert_eq!(v["convergence"]["converged"], true);
    // F01: requested and effective path, solver and backend are recorded.
    let e = &v["execution"];
    assert_eq!(e["effective_profile"], "standard");
    assert_eq!(e["standard_threshold_ha"], 5e-9);
    assert_eq!(e["verification"], "fixed_point");
    assert_eq!(e["requested"]["standard_threshold_ha"], Value::Null);
    assert_eq!(e["requested"]["reference_strict"], false);
    assert_eq!(e["parallel_model"], "threads");
    assert_eq!(e["real_gamma"], v["gamma_real_used"]);
    let solver = e["effective_solver"].as_str().unwrap();
    if v["gamma_real_used"] == true {
        assert!(solver.starts_with("real_gamma_"), "{solver}");
    } else {
        assert_eq!(solver, "complex_davidson");
    }
    assert!(e["backend"]["linear_algebra"].is_string());
    assert!(
        e["non_default_switches"]
            .get("BAANDI_REFERENCE_STRICT")
            .is_none()
    );
    let attempts = v["verification"]["attempts"].as_array().unwrap();
    assert!(!attempts.is_empty());
    let last = attempts.last().unwrap();
    assert_eq!(last["passed"], true);
    assert!(last["coulomb_error_ha"].as_f64().unwrap() < 5e-9);
    assert!(last["h_vectors"].as_u64().unwrap() > 0);
    let seconds = v["timings"]["verification_s"].as_f64().unwrap();
    assert!(seconds > 0.0 && seconds <= v["timings"]["total_s"].as_f64().unwrap());
    assert_eq!(
        v["work"]["certificate_solves"].as_u64().unwrap(),
        attempts.len() as u64
    );
    assert_eq!(v["forces"]["status"], "available");
    assert_eq!(
        v["forces"]["preconditions"]["convergence_profile"],
        "standard"
    );
    assert_eq!(
        v["forces"]["preconditions"]["fixed_point_verification_passed"],
        true
    );
    assert!(
        Command::new("python3")
            .arg("-I")
            .arg(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/validate_force_schema.py"))
            .arg(dir.join("on/results.json"))
            .status()
            .unwrap()
            .success()
    );
    let worst = v["forces"]["values"]
        .as_array()
        .unwrap()
        .iter()
        .zip(u["forces"]["values"].as_array().unwrap())
        .flat_map(|(a, b)| {
            (0..3).map(move |i| (a[i].as_f64().unwrap() - b[i].as_f64().unwrap()).abs())
        })
        .fold(0.0_f64, f64::max);
    println!(
        "single path: attempts {}, verification_s {seconds:.6}, max |F_standard - F_strict| {worst:e} Ha/Bohr, dE {:e} Ha",
        attempts.len(),
        v["total_energy"]["ha"].as_f64().unwrap() - u["total_energy"]["ha"].as_f64().unwrap()
    );
    assert!(worst < 1e-4, "{worst:e} Ha/Bohr");

    // verify = off: same stop step and delivered state (V4); no verification.
    let unverified = run(&dir, "unverified", None, "verify='off'\n", "");
    let n = ok(&unverified);
    assert_eq!(n["single_path"]["verification"], "off");
    assert_eq!(n["execution"]["requested"]["verify"], "off");
    assert!(n.get("verification").is_none() && n.get("verified").is_none());
    if attempts.len() == 1 {
        assert_eq!(energy_bits(n), energy_bits(v));
        assert_eq!(unverified.density_sha, on.density_sha);
    }
    // Forces need the verification.
    let refused = run(&dir, "refused", None, "verify='off'\n", "forces=true\n");
    assert_eq!(refused.code, Some(1), "{}", refused.stderr);
    assert!(
        refused.stderr.contains("verification"),
        "{}",
        refused.stderr
    );
    // The removed profile key is an input error that names the replacement.
    for profile in ["strict", "standard"] {
        let bad = run(
            &dir,
            &format!("removed-{profile}"),
            None,
            &format!("convergence_profile='{profile}'\n"),
            "",
        );
        assert_eq!(bad.code, Some(1), "{}", bad.stderr);
        for needle in [
            "convergence_profile was removed",
            "scf.standard_threshold_ha",
            "scf.verify",
            "BAANDI_REFERENCE_STRICT",
        ] {
            assert!(bad.stderr.contains(needle), "{}", bad.stderr);
        }
    }
    // Inputs the standard path cannot run keep the strict path, with a reason.
    let locking = run(&dir, "locking", None, "locking='auto'\n", "");
    let l = ok(&locking);
    assert_eq!(l["single_path"]["status"], "not_applicable");
    assert_eq!(l["single_path"]["reason"], "scf.locking = auto");
    assert_eq!(l["execution"]["effective_profile"], "strict");
    assert!(l.get("convergence_profile").is_none());
    // ... and reject standard-only controls rather than dropping them.
    let dropped = run(
        &dir,
        "dropped",
        None,
        "locking='auto'\nstandard_threshold_ha=5e-11\n",
        "",
    );
    assert_eq!(dropped.code, Some(1), "{}", dropped.stderr);
    assert!(
        dropped.stderr.contains("scf.locking = auto"),
        "{}",
        dropped.stderr
    );
    // A tighter threshold on the default path.
    let tight = run(&dir, "tight", None, "standard_threshold_ha=5e-11\n", "");
    let t = ok(&tight);
    assert_eq!(t["standard_threshold_ha"], 5e-11);
    assert_eq!(t["execution"]["requested"]["standard_threshold_ha"], 5e-11);
    assert_eq!(t["verified"], true);
    fs::remove_dir_all(dir).unwrap();
}
