//! ONEPATH-3 (user decision 2026-11-16): collinear spin and Fermi-Dirac
//! occupations run on the single standard path with the fixed-point
//! verification (docs/design/onepath3-spin-metals.md, ADR 0043 amendment
//! ONEPATH-3). The strict path (`BAANDI_REFERENCE_STRICT=on`) is the reference.
//! Synthetic fixtures only (no real pseudopotential data).
#[allow(dead_code)]
mod common;
use serde_json::Value;
use std::{fs, path::Path, process::Command};

/// Ferromagnetic fixed sector: He (s projector, 2 e) and H (local, 1 e) in a
/// skewed cell, N_up = 2, N_down = 1 (M = 1 mu_B), Gamma.
const FERRO: &str = "[system]\nname='onepath3 ferro'\nlattice_unit='bohr'\nlattice_vectors=[[6.0,0.0,0.0],[1.1,5.8,0.0],[0.7,0.9,6.3]]\n[[system.atoms]]\nspecies='He'\nfractional=[0.0,0.0,0.0]\n[[system.atoms]]\nspecies='H'\nfractional=[0.41,0.37,0.33]\n[pseudopotentials]\nHe='nl.xml'\nH='loc.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pw'\n[spin]\nmode='collinear'\nn_up=2\nn_down=1\n";
/// Smeared metal: one H (local, 1 e) in a skewed cell, 2x2x2 mesh, FD
/// theta = 0.01 Ha, 12 retained bands.
const METAL: &str = "[system]\nname='onepath3 metal'\nlattice_unit='bohr'\nlattice_vectors=[[5.0,0.0,0.0],[0.6,5.2,0.0],[0.4,0.5,5.4]]\n[[system.atoms]]\nspecies='H'\nfractional=[0.0,0.0,0.0]\n[pseudopotentials]\nH='loc.xml'\n[basis]\necutwfc_ry=5.0\nkgrid=[2,2,2]\n[xc]\nfunctional='lda-pz'\n[occupations]\nmode='fermi_dirac'\nsmearing_ha=0.01\nband_cap=12\n";
/// Spin-polarized smeared fixture: one H (local, 1 e), unconstrained moment
/// (single Fermi level), seed 0.5, Gamma.
const SPIN_METAL: &str = "[system]\nname='onepath3 spin metal'\nlattice_unit='bohr'\nlattice_vectors=[[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0]]\n[[system.atoms]]\nspecies='H'\nfractional=[0.0,0.0,0.0]\n[pseudopotentials]\nH='loc.xml'\n[basis]\necutwfc_ry=5.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pw'\n[occupations]\nmode='fermi_dirac'\nsmearing_ha=0.01\nband_cap=12\n[spin]\nmode='collinear'\ninitial_moment_seed=0.5\n";

fn setup(tag: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("onepath3-{tag}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    for (file, z, rc, nonlocal, element) in [
        ("loc.xml", 1., 1.5, false, "H"),
        ("nl.xml", 2., 1.2, true, "He"),
    ] {
        let upf = common::gaussian_with(z, rc, nonlocal, 401)
            .replace("element=\"Xx\"", &format!("element=\"{element}\""));
        fs::write(dir.join(file), &upf).unwrap();
        // Spin-capable functional label for the lda-pw fixtures.
        fs::write(
            dir.join(file.replace(".xml", "-pw.xml")),
            upf.replace("SLA PZ", "SLA PW"),
        )
        .unwrap();
    }
    dir
}

struct Out {
    code: Option<i32>,
    text: String,
    result: Option<Value>,
}

/// `command`: run, validate or relax; `strict`: BAANDI_REFERENCE_STRICT=on.
fn cli(dir: &Path, command: &str, name: &str, input: &str, strict: bool) -> Out {
    let input = if input.contains("functional='lda-pw'") {
        input
            .replace("'loc.xml'", "'loc-pw.xml'")
            .replace("'nl.xml'", "'nl-pw.xml'")
    } else {
        input.to_owned()
    };
    let _ = fs::remove_dir_all(dir.join(name));
    fs::write(dir.join(format!("{name}.toml")), input).unwrap();
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_dft-cli"));
    cmd.current_dir(dir);
    if command == "relax" {
        cmd.args([command, &format!("{name}.toml"), "--output-dir", name]);
    } else {
        cmd.args([command, &format!("{name}.toml")]);
    }
    if strict {
        cmd.env("BAANDI_REFERENCE_STRICT", "on");
    } else {
        cmd.env_remove("BAANDI_REFERENCE_STRICT");
    }
    let out = cmd.output().unwrap();
    Out {
        code: out.status.code(),
        text: format!(
            "{}\n{}",
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        ),
        result: fs::read(dir.join(name).join("results.json"))
            .ok()
            .map(|b| serde_json::from_slice(&b).unwrap()),
    }
}

fn ok(o: &Out) -> &Value {
    assert_eq!(o.code, Some(0), "{}", o.text);
    o.result.as_ref().unwrap()
}

fn with(base: &str, name: &str, scf: &str, output: &str) -> String {
    format!("{base}[scf]\n{scf}[output]\ndirectory='{name}'\n{output}")
}

/// Energy compared on the variational quantity: F under Fermi-Dirac, E else.
fn energy(v: &Value) -> f64 {
    v["thermodynamics"]["free_energy_ha"]
        .as_f64()
        .unwrap_or_else(|| v["total_energy"]["ha"].as_f64().unwrap())
}

fn numbers(v: &Value) -> Vec<f64> {
    v.as_array()
        .map(|a| a.iter().map(|x| x.as_f64().unwrap()).collect())
        .unwrap_or_default()
}

/// sum_k w_k sum_{s,b} |f_a - f_b| (electrons); occupations exclude k weights.
fn occupation_difference(a: &Value, b: &Value) -> f64 {
    let block = |v: &Value| {
        let o = if v["occupations"]["occupations_electrons"].is_object() {
            &v["occupations"]["occupations_electrons"]
        } else {
            &v["state_contracts"]["occupations"]["occupations_electrons"]
        };
        let shape = numbers(&o["shape"]);
        (numbers(&o["values"]), shape[1] as usize, shape[2] as usize)
    };
    let weights = numbers(&a["state_contracts"]["occupations"]["k_weights"]);
    let ((fa, nk, nb), (fb, _, _)) = (block(a), block(b));
    fa.iter()
        .zip(&fb)
        .enumerate()
        .map(|(i, (x, y))| weights[(i / nb) % nk] * (x - y).abs())
        .sum()
}

/// (b), (c): standard path (default 5e-9 Ha, verification on) against the
/// strict reference on three fixtures. Bounds fixed before measuring
/// (docs/design/onepath3-spin-metals.md section 7): |dE| (F under FD) <= 5e-9
/// Ha; |dM| <= 1e-4 mu_B; |dmu| <= 1e-4 Ha; sum w |df| <= 1e-4 electrons.
/// cargo test --release -p dft-cli --test onepath_spin_metals -- --ignored --nocapture
#[test]
#[ignore = "release-only: three fixtures, standard against strict reference"]
fn standard_spin_and_metals_agree_with_strict_reference() {
    let dir = setup("compare");
    for (fixture, base) in [
        ("ferro", FERRO),
        ("metal", METAL),
        ("spin-metal", SPIN_METAL),
    ] {
        let bands = if base == FERRO { "nbands=4\n" } else { "" };
        let strict_name = format!("{fixture}-strict");
        let strict = cli(
            &dir,
            "run",
            &strict_name,
            &with(base, &strict_name, bands, ""),
            true,
        );
        let s = ok(&strict);
        assert_eq!(s["single_path"]["status"], "reference_strict");
        assert_eq!(s["convergence"]["converged"], true);
        let name = format!("{fixture}-standard");
        let standard = cli(&dir, "run", &name, &with(base, &name, bands, ""), false);
        let r = ok(&standard);
        assert_eq!(r["single_path"]["status"], "applied", "{fixture}");
        assert_eq!(r["single_path"]["reason"], Value::Null);
        assert_eq!(r["execution"]["effective_solver"], "complex_davidson");
        assert_eq!(r["convergence"]["converged"], true);
        assert_eq!(r["verified"], true, "{}", r["verification"]);
        let last = r["verification"]["attempts"]
            .as_array()
            .unwrap()
            .last()
            .unwrap()
            .clone();
        assert_eq!(last["passed"], true);
        assert!(last["coulomb_error_ha"].as_f64().unwrap() < 5e-9);
        let de = (energy(r) - energy(s)).abs();
        let dm = (r["spin"]["moment_mu_b"].as_f64().unwrap()
            - s["spin"]["moment_mu_b"].as_f64().unwrap())
        .abs();
        let mu = |v: &Value| numbers(&v["occupations"]["chemical_potential"]["values_ha"]);
        let dmu = mu(r)
            .iter()
            .zip(mu(s))
            .map(|(a, b)| (a - b).abs())
            .fold(0., f64::max);
        let df = occupation_difference(r, s);
        let iterations = |v: &Value| v["convergence"]["iterations"].as_u64().unwrap();
        println!(
            "ONEPATH-3 {fixture}: E_std={:.15e} E_strict={:.15e} |dE|={de:.3e} Ha |dM|={dm:.3e} mu_B |dmu|={dmu:.3e} Ha sum w|df|={df:.3e} e; SCF std={} strict={} attempts={} last E_H={:.3e} Ha residual estimate={}",
            energy(r),
            energy(s),
            iterations(r),
            iterations(s),
            r["verification"]["attempts"].as_array().unwrap().len(),
            last["coulomb_error_ha"].as_f64().unwrap(),
            last["residual_bound_ha"]
        );
        assert!(de <= 5e-9, "{fixture}: |dE| {de:e} Ha");
        assert!(dm <= 1e-4, "{fixture}: |dM| {dm:e} mu_B");
        assert!(dmu <= 1e-4, "{fixture}: |dmu| {dmu:e} Ha");
        assert!(df <= 1e-4, "{fixture}: sum w|df| {df:e} electrons");
    }
    let _ = fs::remove_dir_all(dir);
}

/// (c) forces of a verified spin state, (d) path records and input errors.
#[test]
fn spin_and_metal_path_records_and_input_errors() {
    let dir = setup("records");
    // Ferromagnetic fixture on the default path, with forces from the
    // verified state (standard preconditions recorded).
    let r = cli(
        &dir,
        "run",
        "ferro",
        &with(FERRO, "ferro", "nbands=4\n", "forces=true\n"),
        false,
    );
    let v = ok(&r);
    assert_eq!(
        (&v["single_path"]["status"], &v["single_path"]["reason"]),
        (&Value::from("applied"), &Value::Null)
    );
    assert_eq!(v["single_path"]["verification"], "fixed_point");
    assert_eq!(v["execution"]["path_status"], "applied");
    assert_eq!(v["execution"]["effective_solver"], "complex_davidson");
    assert_eq!(v["convergence_profile"], "standard");
    assert_eq!(v["verified"], true);
    assert!(v["timings"]["verification_s"].is_number());
    assert!(v["convergence"]["history"][0]["standard_error_ha"].is_number());
    let p = &v["forces"]["preconditions"];
    assert_eq!(p["convergence_profile"], "standard");
    assert_eq!(p["fixed_point_verification_passed"], true);
    assert!(
        Command::new("python3")
            .arg("-I")
            .arg(Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/validate_force_schema.py"))
            .arg(dir.join("ferro/results.json"))
            .status()
            .unwrap()
            .success()
    );
    // verify = off: an unverified standard run (converged, no record).
    let off = cli(
        &dir,
        "run",
        "ferro-off",
        &with(FERRO, "ferro-off", "nbands=4\nverify='off'\n", ""),
        false,
    );
    let v = ok(&off);
    assert_eq!(v["single_path"]["verification"], "off");
    assert_eq!(v["convergence"]["converged"], true);
    assert!(v.get("verification").is_none() && v.get("verified").is_none());
    // Standard-only keys are accepted on spin and FD inputs.
    let tight = cli(
        &dir,
        "validate",
        "metal-tight",
        &with(METAL, "metal-tight", "standard_threshold_ha=5e-11\n", ""),
        false,
    );
    assert_eq!(tight.code, Some(0), "{}", tight.text);
    // Input errors (exit 1, no results.json): forces with verify = off, a
    // spin relaxation, the stall guard (closed-shell experiment only).
    for (name, input, command, needle) in [
        (
            "forces-off",
            with(
                FERRO,
                "forces-off",
                "nbands=4\nverify='off'\n",
                "forces=true\n",
            ),
            "run",
            "fixed-point verification",
        ),
        (
            "relax-spin",
            with(FERRO, "relax-spin", "nbands=4\n", "") + "[relax]\n",
            "relax",
            "unpolarized fixed occupations",
        ),
        (
            "stall",
            with(METAL, "stall", "standard_stall_guard='on'\n", ""),
            "run",
            "closed-shell",
        ),
    ] {
        let o = cli(&dir, command, name, &input, false);
        assert_eq!(o.code, Some(1), "{name}: {}", o.text);
        assert!(o.text.contains(needle), "{name}: {}", o.text);
        assert!(o.result.is_none(), "{name}");
    }
    let _ = fs::remove_dir_all(dir);
}
