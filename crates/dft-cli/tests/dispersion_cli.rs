//! D3-DISPERSION (ADR 0048): `[dispersion]` in the CLI.
//!
//! - Default build: the table is rejected at preflight (no code path runs).
//! - `--features dftd3` with `BAANDI_DFTD3_LIB` set: preflight rejections, and a
//!   synthetic end-to-end run where total energy = SCF energy + E_disp exactly
//!   and the D3 forces/stress are added to the engine's; a relax step uses them.
//! - Ignored: a real PBE graphite run with the PseudoDojo C file (`DFT_NC_DIR`).
//!
//! `BAANDI_DFTD3_LIB=<dest>/lib/libs-dftd3.dylib cargo test -p dft-cli
//! --features dftd3 --test dispersion_cli -- --nocapture`
#[allow(dead_code)]
mod common;
use std::{fs, path::Path, process::Command};

fn scratch(name: &str) -> std::path::PathBuf {
    let root = std::env::var_os("BAANDI_TEST_TMP")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(std::env::temp_dir);
    let dir = root.join(format!("d3-cli-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

/// Two-layer carbon cell with the synthetic LDA Gaussian ion labelled C.
fn synthetic_input(dir: &Path, extra: &str) -> String {
    fs::write(
        dir.join("c.xml"),
        common::gaussian_with(4., 1.0, true, 401).replace("element=\"Xx\"", "element=\"C\""),
    )
    .unwrap();
    format!(
        "[system]\nname='synthetic D3'\nlattice_unit='bohr'\nlattice_vectors=[[4.6563,0.0,0.0],[-2.32815,4.0324,0.0],[0.3,0.1,9.5]]\n\
[[system.atoms]]\nspecies='C'\nfractional=[0.0,0.0,0.25]\n\
[[system.atoms]]\nspecies='C'\nfractional=[0.36,0.64,0.27]\n\
[[system.atoms]]\nspecies='C'\nfractional=[0.02,0.01,0.74]\n\
[[system.atoms]]\nspecies='C'\nfractional=[0.66,0.31,0.77]\n\
[pseudopotentials]\nC='c.xml'\n[basis]\necutwfc_ry=6.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[scf]\nnbands=10\n\
[output]\ndirectory='out'\nforces=true\nstress=true\n{extra}"
    )
}

fn run(dir: &Path, args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_baandi"))
        .current_dir(dir)
        .args(args)
        .output()
        .unwrap()
}

#[cfg(feature = "dftd3")]
fn results(dir: &Path, sub: &str) -> serde_json::Value {
    serde_json::from_slice(&fs::read(dir.join(sub).join("results.json")).unwrap()).unwrap()
}

#[cfg(not(feature = "dftd3"))]
#[test]
fn default_build_rejects_dispersion_at_preflight() {
    let dir = scratch("default");
    fs::write(
        dir.join("input.toml"),
        synthetic_input(&dir, "[dispersion]\nmethod='d3bj'\n"),
    )
    .unwrap();
    let out = run(&dir, &["validate", "input.toml"]);
    let text =
        String::from_utf8_lossy(&out.stderr).to_string() + &String::from_utf8_lossy(&out.stdout);
    assert!(!out.status.success(), "{text}");
    assert!(text.contains("--features dftd3"), "{text}");
    let out = run(&dir, &["run", "input.toml"]);
    assert!(!out.status.success());
    assert!(!dir.join("out/results.json").exists());
    fs::remove_dir_all(dir).ok();
}

#[cfg(feature = "dftd3")]
fn library() -> Option<String> {
    match std::env::var("BAANDI_DFTD3_LIB") {
        Ok(p) if !p.is_empty() => Some(p),
        _ => {
            eprintln!("SKIP: BAANDI_DFTD3_LIB is not set (scripts/fetch_dftd3.py)");
            None
        }
    }
}

#[cfg(feature = "dftd3")]
#[test]
fn preflight_rejections() {
    if library().is_none() {
        return;
    }
    let dir = scratch("reject");
    for (table, needle) in [
        ("method='d3bj'\n", "no D3 parameter set"),
        ("method='d3mbj'\n", "unknown method"),
        (
            "method='d3bj'\ns6=1.0\na1=0.4\n",
            "partial damping parameters",
        ),
        (
            "method='d3zero'\ns6=1.0\ns8=0.7\nrs6=1.2\na1=0.3\n",
            "does not belong",
        ),
        (
            "method='d3bj'\ns6=1.0\ns8=0.7875\na1=0.4289\na2=4.4407\nlibrary='/nonexistent/libs-dftd3.dylib'\n",
            "cannot load",
        ),
        ("method='d3bj'\nunknown=1\n", "unknown field"),
    ] {
        fs::write(
            dir.join("input.toml"),
            synthetic_input(&dir, &format!("[dispersion]\n{table}")),
        )
        .unwrap();
        let out = run(&dir, &["validate", "input.toml"]);
        let text = String::from_utf8_lossy(&out.stderr).to_string()
            + &String::from_utf8_lossy(&out.stdout);
        assert!(!out.status.success(), "{table}: {text}");
        assert!(
            text.contains(needle),
            "{table}: expected {needle:?} in {text}"
        );
    }
    fs::remove_dir_all(dir).ok();
}

#[cfg(feature = "dftd3")]
fn max_abs_diff(a: &serde_json::Value, b: &serde_json::Value) -> f64 {
    match (a, b) {
        (serde_json::Value::Array(x), serde_json::Value::Array(y)) => {
            assert_eq!(x.len(), y.len());
            x.iter()
                .zip(y)
                .map(|(x, y)| max_abs_diff(x, y))
                .fold(0., f64::max)
        }
        _ => (a.as_f64().unwrap() - b.as_f64().unwrap()).abs(),
    }
}

#[cfg(feature = "dftd3")]
fn max_abs(a: &serde_json::Value) -> f64 {
    match a {
        serde_json::Value::Array(x) => x.iter().map(max_abs).fold(0., f64::max),
        _ => a.as_f64().unwrap().abs(),
    }
}

#[cfg(feature = "dftd3")]
fn sum(a: &serde_json::Value, b: &serde_json::Value) -> serde_json::Value {
    match (a, b) {
        (serde_json::Value::Array(x), serde_json::Value::Array(y)) => {
            serde_json::Value::Array(x.iter().zip(y).map(|(x, y)| sum(x, y)).collect())
        }
        _ => serde_json::json!(a.as_f64().unwrap() + b.as_f64().unwrap()),
    }
}

/// Run with and without `[dispersion]` and check the decomposition.
#[cfg(feature = "dftd3")]
fn check_decomposition(dir: &Path, base: &str, table: &str, label: &str) -> serde_json::Value {
    fs::write(dir.join("plain.toml"), base.replace("'out'", "'plain'")).unwrap();
    fs::write(
        dir.join("d3.toml"),
        base.replace("'out'", "'d3'") + "[dispersion]\n" + table,
    )
    .unwrap();
    for name in ["plain.toml", "d3.toml"] {
        let out = run(dir, &["run", name]);
        assert!(
            out.status.success(),
            "{label} {name}: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        if name == "d3.toml" {
            assert!(String::from_utf8_lossy(&out.stdout).contains("Dispersion (DFT-D3"));
        }
    }
    let plain = results(dir, "plain");
    let d3 = results(dir, "d3");
    let disp = &d3["dispersion"];
    let e_disp = disp["energy_ha"].as_f64().unwrap();
    let get = |v: &serde_json::Value, path: &[&str]| -> Option<f64> {
        path.iter().fold(v, |v, k| &v[*k]).as_f64()
    };
    // Every reported total: (path in the plain run, path in the D3 run). The
    // scalar path has total_energy; spin/FD runs have thermodynamics instead.
    let mut pairs = vec![(
        vec!["energies", "total", "ha"],
        vec!["energies", "total_without_dispersion", "ha"],
    )];
    if plain["total_energy"].is_object() {
        pairs.push((
            vec!["total_energy", "ha"],
            vec!["total_energy", "without_dispersion_ha"],
        ));
    }
    if plain["thermodynamics"].is_object() {
        pairs.push((vec!["thermodynamics", "free_energy_ha"], vec![]));
        pairs.push((vec!["thermodynamics", "internal_energy_ha"], vec![]));
    }
    let (mut e_scf, mut e_total) = (f64::NAN, f64::NAN);
    for (path, split) in &pairs {
        let scf = get(&plain, path).unwrap_or_else(|| panic!("{label}: {path:?}"));
        let total = get(&d3, path).unwrap();
        // The SCF does not see the D3 term: same SCF energy bit for bit.
        if !split.is_empty() {
            assert_eq!(
                get(&d3, split).unwrap().to_bits(),
                scf.to_bits(),
                "{label}: {split:?}"
            );
        }
        assert_eq!(
            total.to_bits(),
            (scf + e_disp).to_bits(),
            "{label}: {path:?}"
        );
        (e_scf, e_total) = (scf, total);
    }
    assert_eq!(d3["energies"]["dispersion"]["ha"].as_f64().unwrap(), e_disp);
    assert!(e_disp < 0.);
    let df = max_abs_diff(
        &d3["forces"]["values"],
        &sum(&plain["forces"]["values"], &disp["forces_ha_bohr"]),
    );
    let ds = max_abs_diff(
        &d3["stress"]["values"],
        &sum(&plain["stress"]["values"], &disp["stress_ha_bohr3"]),
    );
    assert_eq!(
        d3["forces"]["per_term"]["dispersion"],
        disp["forces_ha_bohr"]
    );
    assert_eq!(
        d3["stress"]["per_term"]["dispersion"],
        disp["stress_ha_bohr3"]
    );
    let pressure = -(0..3)
        .map(|a| d3["stress"]["values"][a][a].as_f64().unwrap())
        .sum::<f64>()
        / 3.;
    assert!((pressure - d3["stress"]["pressure_ha_bohr3"].as_f64().unwrap()).abs() < 1e-15);
    println!(
        "{label}: E_scf {e_scf:.12} Ha + E_disp {e_disp:.12e} Ha = E_total {e_total:.12} Ha; \
max|F_d3 - (F_scf + F_disp)| {df:.2e} Ha/Bohr; max|sigma_d3 - (sigma_scf + sigma_disp)| {ds:.2e} Ha/Bohr^3; \
max|F_disp| {:.3e} Ha/Bohr; P_disp {:.3e} Ha/Bohr^3",
        max_abs(&disp["forces_ha_bohr"]),
        disp["pressure_ha_bohr3"].as_f64().unwrap()
    );
    // Same SCF (deterministic), so the only difference is the added term.
    assert!(df <= 1e-12, "{label}: force decomposition {df:e}");
    assert!(ds <= 1e-12, "{label}: stress decomposition {ds:e}");
    d3
}

#[cfg(feature = "dftd3")]
#[test]
fn synthetic_end_to_end_run_and_relax_step() {
    if library().is_none() {
        return;
    }
    let dir = scratch("e2e");
    let base = synthetic_input(&dir, "");
    let explicit = "method='d3bj'\ns6=1.0\ns8=0.7875\na1=0.4289\na2=4.4407\n";
    let d3 = check_decomposition(&dir, &base, explicit, "synthetic LDA + explicit d3bj");
    assert_eq!(d3["dispersion"]["parameters"]["source"], "explicit_input");
    assert_eq!(d3["dispersion"]["library"]["version"], "1.7.0");
    // One relaxation step with and without D3: the same SCF per frame, so the
    // first frame differs exactly by the D3 term and the step follows it.
    let relax_run = |name: &str, d3: bool| -> serde_json::Value {
        let mut text = base
            .replace("'out'", &format!("'{name}'"))
            .replace("stress=true\n", "");
        if d3 {
            text = text + "[dispersion]\n" + explicit;
        }
        text += "[relax]\nmax_steps=1\nforce_tol_ha_bohr=1e-12\n";
        fs::write(dir.join(format!("{name}.toml")), text).unwrap();
        let out = run(
            &dir,
            &["relax", &format!("{name}.toml"), "--output-dir", name],
        );
        serde_json::from_slice(
            &fs::read(dir.join(name).join("results.json")).unwrap_or_else(|_| {
                panic!("relax {name}: {}", String::from_utf8_lossy(&out.stderr))
            }),
        )
        .unwrap()
    };
    let relax = relax_run("relax_d3", true);
    let plain = relax_run("relax_plain", false);
    let frames = relax["relaxation"]["frames"].as_array().unwrap();
    let plain_frames = plain["relaxation"]["frames"].as_array().unwrap();
    assert!(
        frames.len() >= 2 && plain_frames.len() >= 2,
        "one step expected"
    );
    let e_disp = d3["dispersion"]["energy_ha"].as_f64().unwrap();
    let de = frames[0]["energy_ha"].as_f64().unwrap()
        - plain_frames[0]["energy_ha"].as_f64().unwrap()
        - e_disp;
    let df = max_abs_diff(
        &frames[0]["forces_ha_bohr"],
        &sum(
            &plain_frames[0]["forces_ha_bohr"],
            &d3["dispersion"]["forces_ha_bohr"],
        ),
    );
    let moved = max_abs_diff(&frames[1]["positions"], &plain_frames[1]["positions"]);
    println!(
        "relax frame 0: E_d3 - E_plain - E_disp = {de:.2e} Ha; max|F_d3 - (F_plain + F_disp)| {df:.2e} Ha/Bohr; \
step 1 position difference with/without D3 {moved:.3e} Bohr"
    );
    assert!(de.abs() <= 1e-12, "relax energy {de:e}");
    assert!(df <= 1e-12, "relax forces {df:e}");
    assert!(moved > 1e-6, "the step did not follow the D3 forces");
    assert_eq!(relax["dispersion"]["method"], "d3bj");
    fs::remove_dir_all(dir).ok();
}

/// Real PBE graphite (AB, 4 C, PseudoDojo nc-sr-04 C.upf from `DFT_NC_DIR`),
/// library D3(BJ) and D3(0) parameters for PBE.
/// `DFT_NC_DIR=... BAANDI_DFTD3_LIB=... cargo test --release -p dft-cli --features dftd3
/// --test dispersion_cli real_pbe_graphite -- --ignored --nocapture`
#[cfg(feature = "dftd3")]
#[test]
#[ignore = "requires DFT_NC_DIR (PseudoDojo C.upf) and the D3 library; release mode"]
fn real_pbe_graphite() {
    if library().is_none() {
        return;
    }
    let nc = std::path::PathBuf::from(std::env::var_os("DFT_NC_DIR").expect("set DFT_NC_DIR"));
    let dir = scratch("graphite");
    fs::copy(nc.join("C.upf"), dir.join("C.upf")).unwrap();
    let base = "[system]\nname='graphite AB PBE-D3'\nlattice_unit='bohr'\nlattice_vectors=[[4.6563,0.0,0.0],[-2.32815,4.032531,0.0],[0.0,0.0,12.682]]\n\
[[system.atoms]]\nspecies='C'\nfractional=[0.0,0.0,0.25]\n\
[[system.atoms]]\nspecies='C'\nfractional=[0.0,0.0,0.75]\n\
[[system.atoms]]\nspecies='C'\nfractional=[0.3333333333333333,0.6666666666666666,0.25]\n\
[[system.atoms]]\nspecies='C'\nfractional=[0.6666666666666666,0.3333333333333333,0.76]\n\
[pseudopotentials]\nC='C.upf'\n[basis]\necutwfc_ry=40.0\nkgrid=[2,2,1]\n[xc]\nfunctional='pbe'\n\
[occupations]\nmode='fermi_dirac'\nsmearing_ha=0.01\n\
[output]\ndirectory='out'\nforces=true\nstress=true\n"
        .to_string();
    for (table, label) in [
        ("method='d3bj'\n", "PBE-D3(BJ) graphite"),
        ("method='d3zero'\n", "PBE-D3(0) graphite"),
    ] {
        let d3 = check_decomposition(&dir, &base, table, label);
        assert_eq!(
            d3["dispersion"]["parameters"]["source"],
            "library_parameter_set"
        );
    }
    fs::remove_dir_all(dir).ok();
}
