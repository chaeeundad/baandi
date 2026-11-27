#[allow(dead_code)]
mod common;
use std::{
    fs,
    path::PathBuf,
    process::Command,
    sync::atomic::{AtomicUsize, Ordering},
};
static NEXT: AtomicUsize = AtomicUsize::new(0);
struct Case {
    dir: PathBuf,
}
impl Case {
    fn new() -> Self {
        let dir = std::env::temp_dir().join(format!(
            "mf03-cli-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&dir).unwrap();
        fs::write(
            dir.join("syn.xml"),
            common::gaussian_with(1.0, 1.0, false, 401)
                .replace("SLA PZ", "SLA PW")
                .replace("element=\"Xx\"", "element=\"H\""),
        )
        .unwrap();
        Self { dir }
    }
    fn input(&self, spin: &str) -> String {
        format!(
            "[system]\nname='synthetic spin'\nlattice_unit='bohr'\nlattice_vectors=[[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0]]\n[[system.atoms]]\nspecies='H'\nfractional=[0.0,0.0,0.0]\n[pseudopotentials]\nH='syn.xml'\n[basis]\necutwfc_ry=5.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pw'\n[scf]\nnbands=4\n[spin]\n{spin}\n"
        )
    }
    fn run(&self, text: &str, command: &str) -> std::process::Output {
        self.run_env(text, command, &[])
    }
    fn run_env(&self, text: &str, command: &str, env: &[(&str, &str)]) -> std::process::Output {
        fs::write(self.dir.join("input.toml"), text).unwrap();
        // These tests check the strict reference (certificate, digests): since
        // ONEPATH-3 spin runs take the standard path unless the reference switch
        // is on (tests/all/onepath_spin_metals.rs covers the standard path).
        Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .env("BAANDI_REFERENCE_STRICT", "on")
            .envs(env.iter().copied())
            .current_dir(&self.dir)
            .args([command, "input.toml"])
            .output()
            .unwrap()
    }
}
impl Drop for Case {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.dir);
    }
}
#[test]
fn spin_preflight_rejections() {
    let c = Case::new();
    for (spin, key) in [
        ("mode='collinear'", "spin.constraint"),
        ("mode='collinear'\nn_up=1", "spin.constraint"),
        (
            "mode='collinear'\nn_up=1\nn_down=0\ntotal_moment=1.0",
            "spin.constraint",
        ),
        ("mode='collinear'\ntotal_moment=0.0", "spin.constraint"),
        ("mode='collinear'\ntotal_moment=3.0", "spin.constraint"),
        ("mode='collinear'\nn_up=1\nn_down=0\naxis='q'", "spin.axis"),
        (
            "mode='collinear'\nn_up=1\nn_down=0\naxis_cartesian=[0.0,0.0,0.0]",
            "spin.axis",
        ),
        ("mode='unpolarized'\nn_up=1", "spin.mode"),
    ] {
        for command in ["validate", "run"] {
            let out = c.run(&c.input(spin), command);
            assert_eq!(out.status.code(), Some(1));
            assert!(
                String::from_utf8_lossy(&out.stderr).contains(key),
                "{:?}",
                out
            );
        }
    }
    let base = c.input("mode='collinear'\nn_up=1\nn_down=0");
    let cases = [
        (
            base.replace(
                "functional='lda-pw'",
                "functional='lda-pz'\nallow_functional_override=true",
            ),
            "xc.spin: required polarized capability unavailable for lda-pz; select pbe or lda-pw",
        ),
        (base.replace("nbands=4", "nbands=1"), "empty band"),
    ];
    for (text, key) in cases {
        for command in ["validate", "run"] {
            let out = c.run(&text, command);
            assert_eq!(out.status.code(), Some(1));
            assert!(
                String::from_utf8_lossy(&out.stderr).contains(key),
                "{:?}",
                out
            );
        }
    } // GT1 (ADR 0050): collinear time reversal is accepted by default; the
    // switch restores the former rejection.
    let tr = base.replace("[basis]", "[basis]\nktime_reversal=true");
    let out = c.run(&tr, "validate");
    assert_eq!(out.status.code(), Some(0), "{out:?}");
    for command in ["validate", "run"] {
        let out = c.run_env(&tr, command, &[("BAANDI_COLLINEAR_TIME_REVERSAL", "off")]);
        assert_eq!(out.status.code(), Some(1));
        assert!(
            String::from_utf8_lossy(&out.stderr).contains("symmetry.operations"),
            "{out:?}"
        );
    }
}
#[test]
fn collinear_cli_certificate_schemas_and_cubes() {
    let c = Case::new();
    let text = c.input("mode='collinear'\ntotal_moment=1.0\naxis='z'\ninitial_moment_seed=0.5");
    let out = c.run(&text, "run");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let path = c.dir.join("results/results.json");
    let a: serde_json::Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
    assert_eq!(a["schema_version"], 2);
    assert_eq!(a["convergence"]["certificate"]["passed"], true);
    assert_eq!(
        a["spin"]["channel_capacity_electrons"],
        serde_json::json!([1, 1])
    );
    let status = Command::new("python3")
        .args(["-I"])
        .arg(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/validate_spin_results.py"))
        .arg(&path)
        .status()
        .unwrap();
    assert!(status.success());
    let out = c.run(
        &text.replace("initial_moment_seed=0.5", "initial_moment_seed=-0.5"),
        "run",
    );
    assert!(out.status.success());
    let b: serde_json::Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
    assert_eq!(
        a["state_contracts"]["model_digest"],
        b["state_contracts"]["model_digest"]
    );
    assert_eq!(
        a["state_contracts"]["state_digest"],
        b["state_contracts"]["state_digest"]
    );
    assert_ne!(a["provenance"], b["provenance"]);
    let out = c.run(
        &text.replace("total_moment=1.0", "total_moment=-1.0"),
        "run",
    );
    assert!(out.status.success());
    let b: serde_json::Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
    assert_ne!(
        a["state_contracts"]["model_digest"],
        b["state_contracts"]["model_digest"]
    );
    assert_ne!(
        a["state_contracts"]["state_digest"],
        b["state_contracts"]["state_digest"]
    );
    assert_eq!(a["total_energy"], b["total_energy"]);
}

#[test]
fn spin_force_extension_v2_validates() {
    let c = Case::new();
    let text = c.input("mode='collinear'\nn_up=1\nn_down=0") + "\n[output]\nforces=true\n";
    let out = c.run(&text, "run");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let path = c.dir.join("results/results.json");
    let r: serde_json::Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
    assert_eq!(r["schema_version"], 2);
    assert_eq!(r["state_contracts"]["extension_version"], 2);
    assert!(
        Command::new("python3")
            .arg("-I")
            .arg(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/validate_force_schema.py"))
            .arg(path)
            .status()
            .unwrap()
            .success()
    );
}

/// Run with P36B_BASELINE_BINARY=/absolute/baseline and optionally
/// P36B_CLI_EVIDENCE=/absolute/evidence cargo test --release -p dft-cli
/// --test spin shared_symmetry_budget_matrix -- --ignored --nocapture.
#[test]
#[ignore = "release-only before/after symmetry and nine-budget matrix"]
fn shared_symmetry_budget_matrix() {
    fn scientific(mut value: serde_json::Value) -> String {
        let object = value.as_object_mut().unwrap();
        for key in [
            "timings",
            "preflight",
            "tool",
            "input_sha256",
            "scf",
            "threading",
        ] {
            object.remove(key);
        }
        // JSON serialization preserves round-trip floats, including signed zero.
        serde_json::to_string(&value).unwrap()
    }
    for mode in ["spin", "fd", "spin_fd", "spin_fd_moment"] {
        let c = Case::new();
        let sector = match mode {
            "fd" => "mode='unpolarized'",
            "spin_fd" => "mode='collinear'",
            "spin_fd_moment" => "mode='collinear'\ntotal_moment=0.3",
            _ => "mode='collinear'\nn_up=1\nn_down=0",
        };
        let mut base = c.input(sector).replace(
            "kpoints=[[0.0,0.0,0.0]]",
            "kgrid=[2,2,2]\nkshift=[0.5,0.5,0.5]\nktime_reversal=false",
        );
        if mode == "fd" {
            base = base.replace("ktime_reversal=false", "ktime_reversal=true");
        }
        if mode != "spin" {
            base = base.replace("nbands=4", "nbands=12");
            base.push_str("\n[occupations]\nmode='fermi_dirac'\nsmearing_ha=0.01\nband_cap=12\n");
        }
        for symmetry in ["off", "auto"] {
            let mut reference = None;
            for threads in [1, 2, 3, 4, 6, 8, 10, 12, 16] {
                let input = base.replace("[scf]", &format!("[scf]\nthreads={threads}"))
                    + &format!("\n[symmetry]\nmode='{symmetry}'\n");
                fs::write(c.dir.join("input.toml"), input).unwrap();
                let mut binaries =
                    vec![("candidate", PathBuf::from(env!("CARGO_BIN_EXE_dft-cli")))];
                if let Some(path) = std::env::var_os("P36B_BASELINE_BINARY") {
                    binaries.push(("baseline", PathBuf::from(path)));
                }
                for (arm, binary) in binaries {
                    let label = format!("{mode}-{symmetry}-T{threads}-{arm}");
                    let out = if let Some(path) = std::env::var_os("P36B_CLI_EVIDENCE") {
                        PathBuf::from(path).join(&label)
                    } else {
                        c.dir.join(&label)
                    };
                    let result = Command::new(binary)
                        .current_dir(&c.dir)
                        .args(["run", "input.toml", "--output-dir"])
                        .arg(&out)
                        .output()
                        .unwrap();
                    assert!(
                        result.status.success(),
                        "{label}: {}",
                        String::from_utf8_lossy(&result.stderr)
                    );
                    let result: serde_json::Value =
                        serde_json::from_slice(&fs::read(out.join("results.json")).unwrap())
                            .unwrap();
                    let signature = (
                        scientific(result),
                        fs::read(out.join("density.bin")).unwrap(),
                    );
                    if let Some(expected) = &reference {
                        assert_eq!(&signature, expected, "{label}");
                    } else {
                        reference = Some(signature);
                    }
                    println!("PASS {label}");
                }
            }
        }
    }
}

/// GT1 (ADR 0050): time reversal in the k reduction against the full mesh.
/// Four synthetic Z=1/2 atoms (two electrons) on a regular tetrahedron in a cubic cell (point group -43m, no
/// inversion), nonlocal synthetic potential with an l = 1 projector, Gamma
/// centred 3x3x1 mesh (TR folds 9 -> 5). For each model (fixed and FD,
/// unpolarized and collinear) the full mesh without symmetry is the oracle;
/// TR only (5 k) and crystal + TR (3 k) must agree in energy (free energy
/// for FD) within the strict class 1e-10 Ha and in forces within 1e-6
/// Ha/Bohr (the bound of the existing path comparisons; strict SCF leaves
/// ~1e-8 Ha/Bohr of force noise). The geometry test time_reversal_folds_only_without_inversion
/// covers the 3D mesh where TR folds beyond the crystal group.
#[test]
fn time_reversal_matches_full_mesh() {
    time_reversal_case("on", &["fixed", "fd", "spin", "spin_fd"], 1e-10);
}

/// The same comparison on the default (standard) reference spin/FD driver,
/// whose own guard used to reject collinear time reversal; the standard
/// threshold is tightened to 1e-12 Ha so that the strict 1e-10 Ha class
/// applies between the full and the reduced run.
#[test]
fn time_reversal_standard_spin_driver() {
    time_reversal_case("off", &["spin", "spin_fd"], 1e-10);
}

fn time_reversal_case(strict: &str, selected: &[&str], tolerance: f64) {
    let c = Case::new();
    fs::write(
        c.dir.join("nl.xml"),
        common::gaussian_with(0.5, 1.0, true, 401)
            .replace("SLA PZ", "SLA PW")
            .replace("element=\"Xx\"", "element=\"H\""),
    )
    .unwrap();
    let a = 0.1;
    let atoms: String = [[a, a, a], [a, -a, -a], [-a, a, -a], [-a, -a, a]]
        .iter()
        .map(|p| format!("[[system.atoms]]\nspecies='H'\nfractional={p:?}\n"))
        .collect();
    let models = [
        ("fixed", "[spin]\nmode='unpolarized'\n", "nbands=4", ""),
        (
            "fd",
            "[spin]\nmode='unpolarized'\n",
            "nbands=8",
            "[occupations]\nmode='fermi_dirac'\nsmearing_ha=0.02\nband_cap=8\n",
        ),
        (
            "spin",
            "[spin]\nmode='collinear'\nn_up=1\nn_down=1\n",
            "nbands=4",
            "",
        ),
        (
            "spin_fd",
            "[spin]\nmode='collinear'\ntotal_moment=0.5\n",
            "nbands=8",
            "[occupations]\nmode='fermi_dirac'\nsmearing_ha=0.02\nband_cap=8\n",
        ),
    ];
    let energy = |r: &serde_json::Value| {
        r["thermodynamics"]["free_energy_ha"]
            .as_f64()
            .or_else(|| r["total_energy"]["ha"].as_f64())
            .or_else(|| r["energies"]["total_ha"].as_f64())
            .unwrap_or_else(|| {
                panic!(
                    "no energy in {:?}",
                    r.as_object().map(|o| o.keys().collect::<Vec<_>>())
                )
            })
    };
    let threshold = if strict == "off" {
        "standard_threshold_ha=1e-12\n"
    } else {
        ""
    };
    for (model, spin, bands, occupations) in models {
        if !selected.contains(&model) {
            continue;
        }
        let mut oracle: Option<(f64, Vec<f64>)> = None;
        for (arm, symmetry, tr, nks) in [
            ("full", "off", false, 9),
            ("tr", "off", true, 5),
            ("crystal_tr", "auto", true, 3),
        ] {
            let text = format!(
                "[system]\nname='tr {model} {arm}'\nlattice_unit='bohr'\nlattice_vectors=[[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0]]\n{atoms}[pseudopotentials]\nH='nl.xml'\n[basis]\necutwfc_ry=6.0\nkgrid=[3,3,1]\nktime_reversal={tr}\nfft_grid=[18,18,18]\n[xc]\nfunctional='lda-pw'\n[scf]\n{bands}\n{threshold}{spin}{occupations}[symmetry]\nmode='{symmetry}'\n[output]\nforces=true\n"
            );
            let out = c.run_env(&text, "run", &[("BAANDI_REFERENCE_STRICT", strict)]);
            assert!(
                out.status.success(),
                "{model}/{arm}: {}",
                String::from_utf8_lossy(&out.stderr)
            );
            let r: serde_json::Value =
                serde_json::from_slice(&fs::read(c.dir.join("results/results.json")).unwrap())
                    .unwrap();
            assert_eq!(r["time_reversal"]["used"], tr, "{model}/{arm}");
            let e = energy(&r);
            assert_eq!(r["forces"]["status"], "available", "{model}/{arm}");
            let f: Vec<f64> = r["forces"]["values"]
                .as_array()
                .unwrap()
                .iter()
                .flat_map(|v| v.as_array().unwrap().iter().map(|x| x.as_f64().unwrap()))
                .collect();
            let k = r["symmetry"]["irreducible_k_count"]
                .as_u64()
                .or_else(|| {
                    r["band_data"]["k_points"]
                        .as_array()
                        .map(|v| v.len() as u64)
                })
                .or_else(|| {
                    r["identity"]["basis"]["k_points"]
                        .as_array()
                        .map(|v| v.len() as u64)
                })
                .or_else(|| r["kpoints"]["count"].as_u64());
            println!("GT1 strict={strict} {model}/{arm}: nks={k:?} E={e:.15e} Ha");
            assert_eq!(k, Some(nks), "{model}/{arm}");
            match &oracle {
                None => oracle = Some((e, f)),
                Some((e0, f0)) => {
                    assert!(
                        (e - e0).abs() <= tolerance,
                        "{model}/{arm}: dE={:e} Ha",
                        e - e0
                    );
                    assert_eq!(f.len(), f0.len(), "{model}/{arm}");
                    let df = f
                        .iter()
                        .zip(f0)
                        .map(|(x, y)| (x - y).abs())
                        .fold(0.0, f64::max);
                    println!(
                        "GT1 {model}/{arm}: dE={:e} Ha max|dF|={df:e} Ha/Bohr",
                        e - e0
                    );
                    assert!(df <= 1e-6, "{model}/{arm}: max|dF|={df:e} Ha/Bohr");
                }
            }
        }
    }
}
