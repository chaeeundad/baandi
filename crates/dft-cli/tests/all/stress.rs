// R6C: pin legacy-specific contracts explicitly; optional atomic data no longer route there.
#[allow(dead_code)]
#[path = "../common/mod.rs"]
mod common;
use std::{fs, process::Command};
#[test]
fn stress_output_identity_schema_and_paths() {
    let dir = std::env::temp_dir().join(format!("p22-stress-cli-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("syn.xml"),
        common::gaussian_with(1., 1.5, false, 401).replace("element=\"Xx\"", "element=\"H\""),
    )
    .unwrap();
    let base = "[system]\nname='synthetic stress'\nlattice_unit='bohr'\nlattice_vectors=[[8.0,0.0,0.0],[0.0,8.0,0.0],[0.0,0.0,8.0]]\n[[system.atoms]]\nspecies='H'\nfractional=[0.0,0.0,0.0]\n[[system.atoms]]\nspecies='H'\nfractional=[0.27,0.23,0.19]\n[pseudopotentials]\nH='syn.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[scf]\nnbands=4\n[output]\ndirectory='out'\nstress=true\n";
    // ONEPATH-2 (ADR 0043 item 5): stress runs the standard path with the
    // fixed-point verification; the stress-off identity run takes the same
    // path. The strict reference (BAANDI_REFERENCE_STRICT=on) is checked too.
    let run_with = |text: &str, reference: &str| {
        fs::write(dir.join("input.toml"), text).unwrap();
        Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .env("BAANDI_REFERENCE_STRICT", "on")
            .current_dir(&dir)
            .env("BAANDI_REFERENCE_STRICT", reference)
            .args(["run", "input.toml"])
            .output()
            .unwrap()
    };
    let run = |text: &str| run_with(text, "on");
    let output = run(base);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(String::from_utf8_lossy(&output.stdout).contains("tensile-positive"));
    let result: serde_json::Value =
        serde_json::from_slice(&fs::read(dir.join("out/results.json")).unwrap()).unwrap();
    let density = fs::read(dir.join("out/density.bin")).unwrap();
    let schema_ok = || {
        Command::new("python3")
            .arg("-I")
            .arg(
                std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                    .join("tests/validate_stress_schema.py"),
            )
            .arg(dir.join("out/results.json"))
            .status()
            .unwrap()
            .success()
    };
    assert!(schema_ok());
    let off = run(&base.replace("stress=true", "stress=false"));
    assert!(off.status.success());
    let no_stress: serde_json::Value =
        serde_json::from_slice(&fs::read(dir.join("out/results.json")).unwrap()).unwrap();
    assert_eq!(result["total_energy"], no_stress["total_energy"]);
    assert_eq!(density, fs::read(dir.join("out/density.bin")).unwrap());
    // PULAY-OPTIN (2026-11-26): the cutoff-scan assessment is off by default.
    // Off: no scan SCF and a not_evaluated record; on: the three scans. The
    // tensor is bit for bit the same either way.
    assert_eq!(
        result["stress"]["pulay_assessment"]["status"],
        "not_evaluated"
    );
    assert_eq!(result["stress"]["cutoff_scf_runs"], 0);
    let with_assessment = base.replace(
        "stress=true\n",
        "stress=true\nstress_pulay_assessment=true\n",
    );
    let started = std::time::Instant::now();
    let on = run(&with_assessment);
    let on_seconds = started.elapsed().as_secs_f64();
    assert!(on.status.success());
    assert!(schema_ok());
    let on: serde_json::Value =
        serde_json::from_slice(&fs::read(dir.join("out/results.json")).unwrap()).unwrap();
    assert_eq!(on["stress"]["cutoff_scf_runs"], 3);
    assert_eq!(
        on["stress"]["pulay_assessment"]["method"],
        "isotropic_cutoff_derivative_5_10_percent"
    );
    let started = std::time::Instant::now();
    assert!(run(base).status.success());
    let off_seconds = started.elapsed().as_secs_f64();
    let bits = |v: &serde_json::Value| -> Vec<u64> {
        v.as_array()
            .unwrap()
            .iter()
            .flat_map(|row| row.as_array().unwrap().iter())
            .map(|x| x.as_f64().unwrap().to_bits())
            .collect()
    };
    for key in ["values", "symmetric", "antisymmetric"] {
        assert_eq!(
            bits(&result["stress"][key]),
            bits(&on["stress"][key]),
            "{key}"
        );
    }
    assert_eq!(
        result["stress"]["pressure_ha_bohr3"]
            .as_f64()
            .unwrap()
            .to_bits(),
        on["stress"]["pressure_ha_bohr3"]
            .as_f64()
            .unwrap()
            .to_bits()
    );
    println!("PULAY-OPTIN wall: assessment off {off_seconds:.3} s, on {on_seconds:.3} s");
    // The strict reference keeps its certificate preconditions.
    let strict = run_with(&with_assessment, "on");
    assert!(strict.status.success());
    assert!(schema_ok());
    let strict: serde_json::Value =
        serde_json::from_slice(&fs::read(dir.join("out/results.json")).unwrap()).unwrap();
    assert_eq!(strict["single_path"]["status"], "reference_strict");
    assert_eq!(
        strict["stress"]["preconditions"]["convergence_profile"],
        "strict"
    );
    assert_eq!(
        strict["stress"]["pulay_assessment"]["scan_profile"],
        "standard_then_strict"
    );
    fs::remove_dir_all(dir).unwrap();
}

#[test]
#[ignore = "release spin/FD envelope scans: cargo test --release -p dft-cli --test all stress:: -- --include-ignored --nocapture"]
fn spin_and_finite_stress_envelopes() {
    let dir = std::env::temp_dir().join(format!("p22-spin-stress-cli-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("syn.xml"),
        common::gaussian_with(1., 1., false, 401)
            .replace("element=\"Xx\"", "element=\"H\"")
            .replace("SLA PZ", "SLA PW"),
    )
    .unwrap();
    for (name, policy, spin) in [
        ("spin", "", "[spin]\nmode='collinear'\nn_up=1\nn_down=1\n"),
        (
            "fd",
            "[occupations]\nmode='fermi_dirac'\nsmearing_ha=0.01\nband_cap=12\n",
            "",
        ),
        (
            "spin-fd",
            "[occupations]\nmode='fermi_dirac'\nsmearing_ha=0.01\nband_cap=12\n",
            "[spin]\nmode='collinear'\ntotal_moment=0.3\n",
        ),
    ] {
        // FD-FORCES: scf.nbands sets the FD band count; FD cases keep their
        // historical count, the cap (12).
        let bands = if policy.is_empty() { "nbands=4\n" } else { "" };
        let text = format!(
            "[system]\nname='synthetic stress'\nlattice_unit='bohr'\nlattice_vectors=[[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0]]\n[[system.atoms]]\nspecies='H'\nfractional=[0.0,0.0,0.0]\n[[system.atoms]]\nspecies='H'\nfractional=[0.27,0.23,0.19]\n[pseudopotentials]\nH='syn.xml'\n[basis]\necutwfc_ry=5.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pw'\n[scf]\n{bands}{policy}{spin}[output]\ndirectory='{name}'\nstress=true\n"
        );
        fs::write(dir.join("input.toml"), text).unwrap();
        let out = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .env("BAANDI_REFERENCE_STRICT", "on")
            .current_dir(&dir)
            .args(["run", "input.toml"])
            .output()
            .unwrap();
        assert!(
            out.status.success(),
            "{name}: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        assert!(
            Command::new("python3")
                .arg("-I")
                .arg(
                    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                        .join("tests/validate_stress_schema.py")
                )
                .arg(dir.join(name).join("results.json"))
                .status()
                .unwrap()
                .success()
        );
        let r: serde_json::Value =
            serde_json::from_slice(&fs::read(dir.join(name).join("results.json")).unwrap())
                .unwrap();
        assert_eq!(
            r["stress"]["quantity"],
            if name == "spin" { "E" } else { "F" }
        );
        assert!(dir.join(name).join("stress.txt").is_file());
    }
    fs::remove_dir_all(dir).unwrap();
}
