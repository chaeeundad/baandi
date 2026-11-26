// R6C: pin legacy-specific contracts explicitly; optional atomic data no longer route there.
#[allow(dead_code)]
mod common;
use std::{fs, process::Command};

#[test]
fn force_output_and_standard_rejection() {
    let dir = std::env::temp_dir().join(format!("p20-force-cli-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("syn.xml"),
        common::gaussian_with(1., 1.5, false, 401).replace("element=\"Xx\"", "element=\"H\""),
    )
    .unwrap();
    let base = "[system]\nname='synthetic forces'\nlattice_unit='bohr'\nlattice_vectors=[[8.0,0.0,0.0],[0.0,8.0,0.0],[0.0,0.0,8.0]]\n[[system.atoms]]\nspecies='H'\nfractional=[0.0,0.0,0.0]\n[[system.atoms]]\nspecies='H'\nfractional=[0.27,0.23,0.19]\n[pseudopotentials]\nH='syn.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[scf]\nnbands=4\n[output]\ndirectory='out'\nforces=true\n";
    fs::write(dir.join("input.toml"), base).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
        .env("BAANDI_REFERENCE_STRICT", "on")
        .current_dir(&dir)
        .args(["run", "input.toml"])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let result: serde_json::Value =
        serde_json::from_slice(&fs::read(dir.join("out/results.json")).unwrap()).unwrap();
    assert!(
        Command::new("python3")
            .arg("-I")
            .arg(
                std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                    .join("tests/validate_force_schema.py")
            )
            .arg(dir.join("out/results.json"))
            .status()
            .unwrap()
            .success()
    );
    let f = &result["forces"];
    assert_eq!(f["status"], "available");
    assert_eq!(f["units"], "Ha/Bohr");
    assert_eq!(f["values"].as_array().unwrap().len(), 2);
    for atom in 0..2 {
        for axis in 0..3 {
            let sum: f64 = ["ewald", "local", "nonlocal", "core"]
                .iter()
                .map(|t| f["per_term"][t][atom][axis].as_f64().unwrap())
                .sum();
            assert!((sum - f["values"][atom][axis].as_f64().unwrap()).abs() < 1e-14);
        }
    }
    assert!(String::from_utf8_lossy(&output.stdout).contains("Cartesian forces"));
    let symmetric =
        base.replace("[0.27,0.23,0.19]", "[0.25,0.0,0.0]") + "\n[symmetry]\nmode='auto'\n";
    fs::write(dir.join("input.toml"), symmetric).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
        .env("BAANDI_REFERENCE_STRICT", "on")
        .current_dir(&dir)
        .args(["run", "input.toml"])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let result: serde_json::Value =
        serde_json::from_slice(&fs::read(dir.join("out/results.json")).unwrap()).unwrap();
    assert_eq!(result["forces"]["status"], "available");
    for value in result["forces"]["net_force"].as_array().unwrap() {
        assert!(value.as_f64().unwrap().abs() < 1e-12);
    }
    fs::remove_dir_all(dir).unwrap();
}

/// FD-FORCES: the default (reference) driver writes forces and stress for
/// Fermi-Dirac and collinear deliveries, honours scf.nbands below the band cap
/// and rejects scf.nbands above it at preflight.
#[test]
fn reference_fermi_dirac_and_spin_write_forces_and_stress() {
    let dir = std::env::temp_dir().join(format!("fd-force-cli-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let pp = common::gaussian_with(1., 1.5, true, 161).replace("element=\"Xx\"", "element=\"H\"");
    fs::write(dir.join("syn.xml"), &pp).unwrap();
    fs::write(dir.join("syn-pw.xml"), pp.replace("SLA PZ", "SLA PW")).unwrap();
    let base = "[system]\nname='synthetic FD forces'\nlattice_unit='bohr'\nlattice_vectors=[[7.0,0.0,0.0],[0.3,7.0,0.0],[0.0,0.2,7.0]]\n[[system.atoms]]\nspecies='H'\nfractional=[0.1,0.2,0.3]\n[[system.atoms]]\nspecies='H'\nfractional=[0.37,0.43,0.49]\n[pseudopotentials]\nH='syn.xml'\n[basis]\necutwfc_ry=4.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[scf]\nNBANDS\n[output]\ndirectory='out'\nforces=true\nstress=true\n";
    let run = |text: String| {
        let _ = fs::remove_dir_all(dir.join("out"));
        fs::write(dir.join("input.toml"), text).unwrap();
        Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .current_dir(&dir)
            .args(["run", "input.toml"])
            .output()
            .unwrap()
    };
    let check = |output: std::process::Output, quantity: &str| -> serde_json::Value {
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let result: serde_json::Value =
            serde_json::from_slice(&fs::read(dir.join("out/results.json")).unwrap()).unwrap();
        assert_eq!(result["output_profile"], "reference_spin_fd_v1");
        assert_eq!(result["verified"], true);
        for key in ["forces", "stress"] {
            let block = &result[key];
            assert_eq!(block["status"], "available", "{key}: {block}");
            assert_eq!(block["quantity"], quantity);
            assert_eq!(block["verification"], "verified");
            assert_eq!(block["preconditions"]["convergence_profile"], "standard");
            assert_eq!(
                block["preconditions"]["fixed_point_verification_passed"],
                true
            );
        }
        let f = &result["forces"];
        assert_eq!(f["values"].as_array().unwrap().len(), 2);
        for atom in 0..2 {
            for axis in 0..3 {
                let sum: f64 = ["ewald", "local", "nonlocal", "core"]
                    .iter()
                    .map(|t| f["per_term"][t][atom][axis].as_f64().unwrap())
                    .sum();
                assert!((sum - f["values"][atom][axis].as_f64().unwrap()).abs() < 1e-14);
            }
        }
        assert!(
            result["stress"]["pressure_ha_bohr3"]
                .as_f64()
                .unwrap()
                .is_finite()
        );
        // PULAY-OPTIN: the assessment is opt-in; its record and the scan count
        // agree either way.
        let pulay = &result["stress"]["pulay_assessment"];
        if pulay["status"] == "not_evaluated" {
            assert_eq!(result["stress"]["cutoff_scf_runs"], 0);
        } else {
            assert_eq!(pulay["method"], "isotropic_cutoff_derivative_5_10_percent");
            assert_eq!(result["stress"]["cutoff_scf_runs"], 3);
        }
        let stdout = String::from_utf8_lossy(&output.stdout);
        assert!(stdout.contains("Cartesian forces") && stdout.contains("stress"));
        result
    };
    let fd = |theta: f64| {
        format!("[occupations]\nmode='fermi_dirac'\nsmearing_ha={theta}\nband_cap=16\n")
    };
    // Without scf.nbands the count grows automatically up to the cap and the
    // derivatives use the grown state (TS = 9.9e-4 Ha).
    // The assessment is enabled here (spin/FD scans) and must leave the tensor
    // bit for bit unchanged against the default-off run.
    let result = check(
        run(base.replace("NBANDS", "").replace(
            "stress=true\n",
            "stress=true\nstress_pulay_assessment=true\n",
        ) + &fd(0.03)),
        "F",
    );
    assert_eq!(result["stress"]["cutoff_scf_runs"], 3);
    let off = check(run(base.replace("NBANDS", "") + &fd(0.03)), "F");
    assert_eq!(off["stress"]["pulay_assessment"]["status"], "not_evaluated");
    assert_eq!(off["stress"]["values"], result["stress"]["values"]);
    let growth = &result["reference_scf"]["band_count_growth"];
    assert_eq!(growth["mode"], "automatic", "{growth}");
    assert_eq!(growth["band_cap"], 16);
    assert!(growth["final_bands"].as_u64().unwrap() <= 16);
    assert_eq!(result["band_data"]["shape"][2], growth["final_bands"]);
    assert!(
        result["thermodynamics"]["entropy_term_ha"]
            .as_f64()
            .unwrap()
            > 1e-4
    );
    // scf.nbands below the cap sets the band count of the FD run.
    let result = check(run(base.replace("NBANDS", "nbands=12") + &fd(0.02)), "F");
    assert_eq!(result["band_data"]["shape"][2], 12);
    assert_eq!(
        result["reference_scf"]["band_count_growth"]["mode"],
        "fixed"
    );
    // scf.nbands above the explicit cap is rejected before any SCF.
    let output = run(base.replace("NBANDS", "nbands=20") + &fd(0.03));
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("scf.nbands = 20 exceeds occupations.band_cap = 16"),
        "{stderr}"
    );
    assert!(!dir.join("out/results.json").exists());
    // Collinear fixed sector (E) and collinear FD at fixed moment (F).
    let spin = base
        .replace("H='syn.xml'", "H='syn-pw.xml'")
        .replace("functional='lda-pz'", "functional='lda-pw'");
    check(
        run(spin.replace("NBANDS", "nbands=4") + "[spin]\nmode='collinear'\nn_up=2\nn_down=0\n"),
        "E",
    );
    check(
        run(spin.replace("NBANDS", "")
            + &fd(0.02)
            + "[spin]\nmode='collinear'\ntotal_moment=0.5\n"),
        "F",
    );
    fs::remove_dir_all(dir).unwrap();
}
