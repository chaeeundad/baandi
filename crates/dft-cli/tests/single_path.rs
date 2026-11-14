//! ADR 0043 single SCF path behind `BAANDI_SINGLE_PATH` (default off).
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

fn run(dir: &PathBuf, name: &str, switch: Option<&str>, scf: &str, output: &str) -> Run {
    let text = format!("{BASE}[scf]\nnbands=4\n{scf}[output]\ndirectory='{name}'\n{output}");
    let input = format!("{name}.toml");
    fs::write(dir.join(&input), text).unwrap();
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_dft-cli"));
    cmd.current_dir(dir).args(["run", &input]);
    match switch {
        Some(v) => cmd.env("BAANDI_SINGLE_PATH", v),
        None => cmd.env_remove("BAANDI_SINGLE_PATH"),
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
fn single_path_switch_contract() {
    let dir = std::env::temp_dir().join(format!("adr43-single-path-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("syn.xml"),
        common::gaussian_with(1., 1.5, false, 401).replace("element=\"Xx\"", "element=\"H\""),
    )
    .unwrap();

    // Switch off (unset and explicit) is the strict path, bit for bit.
    let unset = run(&dir, "unset", None, "", "forces=true\n");
    let off = run(&dir, "off", Some("off"), "", "forces=true\n");
    let (u, o) = (ok(&unset), ok(&off));
    assert_eq!(unset.density_sha, off.density_sha);
    assert_eq!(energy_bits(u), energy_bits(o));
    assert_eq!(u["forces"]["values"], o["forces"]["values"]);
    for v in [u, o] {
        assert!(v.get("single_path").is_none() && v.get("verification").is_none());
        assert!(v.get("convergence_profile").is_none(), "strict path");
        assert!(v["timings"].get("verification_s").is_none());
    }
    assert_eq!(o["scf"]["switches"]["BAANDI_SINGLE_PATH"]["value"], "off");

    // Switch on, explicit strict: today's strict behaviour.
    let strict = run(
        &dir,
        "strict",
        Some("on"),
        "convergence_profile='strict'\n",
        "",
    );
    let s = ok(&strict);
    assert_eq!(strict.density_sha, unset.density_sha);
    assert_eq!(energy_bits(s), energy_bits(u));
    assert_eq!(s["single_path"]["status"], "explicit_profile");
    assert_eq!(s["single_path"]["verification"], "off");

    // Switch on, plain input: standard profile, verification, forces.
    let on = run(&dir, "on", Some("on"), "", "forces=true\n");
    let v = ok(&on);
    assert_eq!(v["single_path"]["status"], "applied");
    assert_eq!(v["single_path"]["verification"], "fixed_point");
    assert_eq!(v["convergence_profile"], "standard");
    assert_eq!(v["standard_threshold_ha"], 5e-9);
    assert_eq!(v["certified"], false);
    assert_eq!(v["verified"], true);
    assert_eq!(v["convergence"]["converged"], true);
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
    let unverified = run(&dir, "unverified", Some("on"), "verify='off'\n", "");
    let n = ok(&unverified);
    assert_eq!(n["single_path"]["verification"], "off");
    assert!(n.get("verification").is_none() && n.get("verified").is_none());
    if attempts.len() == 1 {
        assert_eq!(energy_bits(n), energy_bits(v));
        assert_eq!(unverified.density_sha, on.density_sha);
    }
    // Forces need the verification on the single path.
    let refused = run(
        &dir,
        "refused",
        Some("on"),
        "verify='off'\n",
        "forces=true\n",
    );
    assert_eq!(refused.code, Some(1), "{}", refused.stderr);
    assert!(
        refused.stderr.contains("verification"),
        "{}",
        refused.stderr
    );
    // verify with an explicit strict profile is an input error.
    let bad = run(
        &dir,
        "bad",
        None,
        "convergence_profile='strict'\nverify='off'\n",
        "",
    );
    assert_eq!(bad.code, Some(1), "{}", bad.stderr);
    assert!(bad.stderr.contains("scf.verify"), "{}", bad.stderr);
    // Inputs the standard profile cannot run keep the strict path, with a reason.
    let locking = run(&dir, "locking", Some("on"), "locking='auto'\n", "");
    let l = ok(&locking);
    assert_eq!(l["single_path"]["status"], "not_applicable");
    assert_eq!(l["single_path"]["reason"], "scf.locking = auto");
    assert!(l.get("convergence_profile").is_none());
    // A threshold without a profile belongs to the single path only.
    let tight = run(
        &dir,
        "tight",
        Some("on"),
        "standard_threshold_ha=5e-11\n",
        "",
    );
    let t = ok(&tight);
    assert_eq!(t["standard_threshold_ha"], 5e-11);
    assert_eq!(t["verified"], true);
    let orphan = run(&dir, "orphan", None, "standard_threshold_ha=5e-11\n", "");
    assert_eq!(orphan.code, Some(1), "{}", orphan.stderr);
    fs::remove_dir_all(dir).unwrap();
}
