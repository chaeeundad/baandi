//! ONEPATH-1 (user decision 2026-11-16): forces and fixed-cell relaxation run
//! on the single standard path with the fixed-point verification (ADR 0043
//! item 5). The strict path (`BAANDI_REFERENCE_STRICT=on`) is the reference.
//! Synthetic fixtures only (no real pseudopotential data).
#[allow(dead_code)]
#[path = "../common/mod.rs"]
mod common;
use serde_json::Value;
use std::{fs, path::Path, process::Command};

/// Two-atom H fixture in a cubic 8 Bohr cell (local-only Gaussian ions).
const CUBIC: &str = "[system]\nname='onepath cubic'\nlattice_unit='bohr'\nlattice_vectors=[[8.0,0.0,0.0],[0.0,8.0,0.0],[0.0,0.0,8.0]]\n[[system.atoms]]\nspecies='H'\nfractional=[0.0,0.0,0.0]\n[[system.atoms]]\nspecies='H'\nfractional=[0.27,0.23,0.19]\n[pseudopotentials]\nH='loc.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n";
/// Four-atom fixture (two species, six electrons) in a skewed cell with s, s and p projectors (nonlocal
/// Gaussian ions): a different cell shape and projector structure.
const SKEW: &str = "[system]\nname='onepath skew'\nlattice_unit='bohr'\nlattice_vectors=[[7.0,0.0,0.0],[1.4,6.6,0.0],[0.9,1.2,7.4]]\n[[system.atoms]]\nspecies='He'\nfractional=[0.0,0.0,0.0]\n[[system.atoms]]\nspecies='He'\nfractional=[0.31,0.22,0.17]\n[[system.atoms]]\nspecies='H'\nfractional=[0.62,0.55,0.48]\n[[system.atoms]]\nspecies='H'\nfractional=[0.71,0.83,0.66]\n[pseudopotentials]\nHe='nl.xml'\nH='loc.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n";

fn setup(tag: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("onepath1-{tag}-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("loc.xml"),
        common::gaussian_with(1., 1.5, false, 401).replace("element=\"Xx\"", "element=\"H\""),
    )
    .unwrap();
    fs::write(
        dir.join("nl.xml"),
        common::gaussian_with(2., 1.2, true, 401).replace("element=\"Xx\"", "element=\"He\""),
    )
    .unwrap();
    dir
}

struct Out {
    code: Option<i32>,
    text: String,
    result: Option<Value>,
}

/// `command`: run or relax; `strict`: BAANDI_REFERENCE_STRICT=on.
fn cli(dir: &Path, command: &str, name: &str, input: &str, strict: bool) -> Out {
    // R6C: this suite validates the legacy standard/strict contracts.
    let input = input.replace("[scf]", "[scf]\ndriver='legacy'");
    fs::write(dir.join(format!("{name}.toml")), input).unwrap();
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_dft-cli"));
    cmd.current_dir(dir);
    if command == "run" {
        cmd.args(["run", &format!("{name}.toml")]);
    } else {
        cmd.args([command, &format!("{name}.toml"), "--output-dir", name]);
    }
    if strict {
        cmd.env("BAANDI_REFERENCE_STRICT", "on");
    } else {
        cmd.env_remove("BAANDI_REFERENCE_STRICT");
    }
    cmd.env_remove("BAANDI_RELAX_EXTRAPOLATION");
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

fn forces(v: &Value) -> Vec<f64> {
    v["forces"]["values"]
        .as_array()
        .unwrap()
        .iter()
        .flat_map(|a| (0..3).map(move |i| a[i].as_f64().unwrap()))
        .collect()
}

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0., f64::max)
}

/// (a) and (c): forces of the default standard path with verification against
/// the strict reference, at the default threshold 5e-9 Ha and at 5e-11 Ha, on
/// two fixtures (cubic local-only, skewed with nonlocal projectors).
///
/// Tolerance (fixed before measuring, not tuned): the largest component
/// difference must stay below 1e-3 of the largest strict component at 5e-9 Ha
/// (the 0.1 % recorded on Si in the ledger) and below 1e-4 of it at 5e-11 Ha
/// (the strict class is two orders tighter). Measured values are printed.
#[test]
fn standard_forces_agree_with_strict_reference() {
    let dir = setup("forces");
    for (fixture, base) in [("cubic", CUBIC), ("skew", SKEW)] {
        let input = |name: &str, scf: &str| {
            format!("{base}[scf]\nnbands=4\n{scf}[output]\ndirectory='{name}'\nforces=true\n")
        };
        let strict = cli(
            &dir,
            "run",
            &format!("{fixture}-strict"),
            &input(&format!("{fixture}-strict"), ""),
            true,
        );
        let s = ok(&strict);
        assert_eq!(s["single_path"]["status"], "reference_strict");
        assert_eq!(
            s["forces"]["preconditions"]["fixed_point_certificate_passed"],
            true
        );
        let fs_strict = forces(s);
        let scale = fs_strict.iter().map(|v| v.abs()).fold(0., f64::max);
        assert!(scale > 1e-4, "fixture has no force to compare: {scale:e}");
        for (threshold, relative) in [(None, 1e-3), (Some(5e-11), 1e-4)] {
            let name = format!(
                "{fixture}-{}",
                threshold.map_or("default".into(), |t: f64| format!("{t:e}"))
            );
            let scf = threshold.map_or(String::new(), |t| format!("standard_threshold_ha={t:e}\n"));
            let standard = cli(&dir, "run", &name, &input(&name, &scf), false);
            let v = ok(&standard);
            let t = threshold.unwrap_or(5e-9);
            // (c) the single path applied and the verification is recorded next to the forces.
            assert_eq!(v["single_path"]["status"], "applied");
            assert_eq!(v["single_path"]["verification"], "fixed_point");
            assert_eq!(v["verified"], true);
            let p = &v["forces"]["preconditions"];
            assert_eq!(p["convergence_profile"], "standard");
            assert_eq!(p["fixed_point_verification_passed"], true);
            assert_eq!(p["standard_threshold_ha"].as_f64().unwrap(), t);
            assert!(p["verification_attempts"].as_u64().unwrap() >= 1);
            assert!(p["verification_coulomb_error_ha"].as_f64().unwrap() < t);
            assert!(
                Command::new("python3")
                    .arg("-I")
                    .arg(
                        Path::new(env!("CARGO_MANIFEST_DIR"))
                            .join("tests/validate_force_schema.py")
                    )
                    .arg(dir.join(&name).join("results.json"))
                    .status()
                    .unwrap()
                    .success()
            );
            let diff = max_abs_diff(&forces(v), &fs_strict);
            let de = v["total_energy"]["ha"].as_f64().unwrap()
                - s["total_energy"]["ha"].as_f64().unwrap();
            println!(
                "ONEPATH-1 forces {fixture} threshold {t:e} Ha: max|F_std - F_strict| {diff:.3e} Ha/Bohr, max|F_strict| {scale:.6e} Ha/Bohr, relative {:.3e}, dE {de:.3e} Ha, attempts {}, coulomb {:.3e} Ha",
                diff / scale,
                p["verification_attempts"],
                p["verification_coulomb_error_ha"].as_f64().unwrap()
            );
            assert!(
                diff <= relative * scale,
                "{fixture} {t:e}: {diff:e} vs {scale:e}"
            );
            assert!(de.abs() <= t, "{fixture} {t:e}: dE {de:e} Ha");
        }
    }
    fs::remove_dir_all(dir).unwrap();
}

/// (d) verify = off with forces or with a relaxation (fixed-cell or, since
/// ONEPATH-2, variable-cell) is an input error.
#[test]
fn verification_off_and_variable_cell_inputs_are_rejected() {
    let dir = setup("reject");
    let scf = "[scf]\nnbands=4\nverify='off'\n";
    let forces = cli(
        &dir,
        "run",
        "forces-off",
        &format!("{CUBIC}{scf}[output]\ndirectory='forces-off'\nforces=true\n"),
        false,
    );
    assert_eq!(forces.code, Some(1), "{}", forces.text);
    assert!(forces.text.contains("verification"), "{}", forces.text);
    let relax = cli(
        &dir,
        "relax",
        "relax-off",
        &format!("{CUBIC}{scf}[relax]\nmax_steps=1\nfixed=[{{atom=0}}]\n"),
        false,
    );
    assert_eq!(relax.code, Some(1), "{}", relax.text);
    assert!(relax.text.contains("verification"), "{}", relax.text);
    assert!(relax.result.is_none());
    // ONEPATH-2: variable-cell relaxation runs the standard path too, so
    // verify = off is an input error there as well.
    let cell = cli(
        &dir,
        "relax",
        "cell-off",
        &format!("{CUBIC}{scf}[relax]\nmode='cell'\nmax_steps=1\n"),
        false,
    );
    assert_eq!(cell.code, Some(1), "{}", cell.text);
    assert!(cell.text.contains("verification"), "{}", cell.text);
    assert!(cell.result.is_none());
    fs::remove_dir_all(dir).unwrap();
}

/// Two H ions in a compact skewed cell (edges about 4 Bohr): the free atom has
/// eight near images of the fixed one, so its minimum is stiff in every
/// direction (curvature of order 0.1 Ha/Bohr^2, estimated from a force scan of
/// this ion: |F| about 0.07 Ha/Bohr at 2.8 Bohr and 0.15 Ha/Bohr at 2.2 Bohr).
const COMPACT: &str = "[system]\nname='onepath compact skew'\nlattice_unit='bohr'\nlattice_vectors=[[4.0,0.0,0.0],[0.4,4.0,0.0],[0.3,0.2,4.2]]\n[[system.atoms]]\nspecies='H'\nfractional=[0.0,0.0,0.0]\n[[system.atoms]]\nspecies='H'\nfractional=[0.4,0.42,0.37]\n[pseudopotentials]\nH='loc.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n";

/// (b) and (c): a short fixed-cell relaxation (one atom fixed, the other
/// displaced in all three directions) on the standard path with verification
/// reaches the same structure as the strict reference. Every ionic step is
/// verified.
///
/// Tolerance: final Cartesian positions agree within 1e-3 Bohr per component
/// for force_tol 1e-4 Ha/Bohr. With a curvature of order 0.1 Ha/Bohr^2, a
/// residual force below 1e-4 Ha/Bohr leaves each run within about 1e-3 Bohr of
/// the minimum. (A first attempt on the cubic two-atom fixture,
/// whose only minimum is the flat cell-symmetric point with curvature about
/// 2e-4 Ha/Bohr^2, ended 5.4e-3 Bohr apart with both runs converged: there the
/// stopping force, not the SCF path, sets the position; recorded in the ledger
/// row ONEPATH-1.)
#[test]
#[ignore = "two multi-step relaxations; run cargo test --release -p dft-cli --test all onepath_forces_relax:: -- --ignored --nocapture"]
fn standard_relaxation_matches_strict_reference() {
    let dir = setup("relax");
    let input = COMPACT.to_string()
        + "[scf]\nnbands=4\n[relax]\nmax_steps=40\nforce_tol_ha_bohr=1e-4\nfixed=[{atom=0}]\n";
    let mut finals = Vec::new();
    for strict in [true, false] {
        let name = if strict { "strict" } else { "standard" };
        let started = std::time::Instant::now();
        let out = cli(&dir, "relax", name, &input, strict);
        let wall = started.elapsed().as_secs_f64();
        let v = ok(&out);
        let r = &v["relaxation"];
        assert_eq!(r["status"], "converged");
        let frames = r["frames"].as_array().unwrap();
        if strict {
            assert_eq!(v["single_path"]["status"], "reference_strict");
            assert!(frames.iter().all(|f| f.get("verification").is_none()));
        } else {
            assert_eq!(v["single_path"]["status"], "applied");
            assert_eq!(v["single_path"]["verification"], "fixed_point");
            assert_eq!(v["scf"]["standard_threshold_ha"], 5e-9);
            for f in frames {
                let ver = &f["verification"];
                assert_eq!(ver["passed"], true, "{f}");
                assert!(ver["coulomb_error_ha"].as_f64().unwrap() < 5e-9);
            }
        }
        let last = frames.last().unwrap();
        let x: Vec<f64> = last["positions"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_f64().unwrap())
            .collect();
        let scf: u64 = frames
            .iter()
            .map(|f| f["scf_iterations"].as_u64().unwrap())
            .sum();
        let h: u64 = frames
            .iter()
            .map(|f| f["scf_h_vectors"].as_u64().unwrap())
            .sum();
        let attempts: u64 = frames
            .iter()
            .filter_map(|f| f["verification"]["attempts"].as_u64())
            .sum();
        println!(
            "ONEPATH-1 relax {name}: ionic steps {}, total SCF iterations {scf}, H vectors {h}, verification attempts {attempts}, final free-atom position {:?} Bohr, final max force {:.3e} Ha/Bohr, wall {wall:.2} s",
            frames.len() - 1,
            &x[3..],
            last["max_force_ha_bohr"].as_f64().unwrap()
        );
        finals.push(x);
    }
    let gap = max_abs_diff(&finals[0], &finals[1]);
    println!("ONEPATH-1 relax: max |x_standard - x_strict| {gap:.3e} Bohr");
    assert!(gap < 1e-3, "{gap:e} Bohr");
    fs::remove_dir_all(dir).unwrap();
}
