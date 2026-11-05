#[allow(dead_code)]
mod common;
use std::{fs, process::Command};
#[test]
fn stress_output_identity_schema_and_standard_rejection() {
    let dir = std::env::temp_dir().join(format!("p22-stress-cli-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("syn.xml"),
        common::gaussian_with(1., 1.5, false, 401).replace("element=\"Xx\"", "element=\"H\""),
    )
    .unwrap();
    let base = "[system]\nname='synthetic stress'\nlattice_unit='bohr'\nlattice_vectors=[[8.0,0.0,0.0],[0.0,8.0,0.0],[0.0,0.0,8.0]]\n[[system.atoms]]\nspecies='H'\nfractional=[0.0,0.0,0.0]\n[[system.atoms]]\nspecies='H'\nfractional=[0.27,0.23,0.19]\n[pseudopotentials]\nH='syn.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[scf]\nnbands=4\n[output]\ndirectory='out'\nstress=true\n";
    let run = |text: &str| {
        fs::write(dir.join("input.toml"), text).unwrap();
        Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .current_dir(&dir)
            .args(["run", "input.toml"])
            .output()
            .unwrap()
    };
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
    assert!(
        Command::new("python3")
            .arg("-I")
            .arg(
                std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                    .join("tests/validate_stress_schema.py")
            )
            .arg(dir.join("out/results.json"))
            .status()
            .unwrap()
            .success()
    );
    let off = run(&base.replace("stress=true", "stress=false"));
    assert!(off.status.success());
    let no_stress: serde_json::Value =
        serde_json::from_slice(&fs::read(dir.join("out/results.json")).unwrap()).unwrap();
    assert_eq!(result["total_energy"], no_stress["total_energy"]);
    assert_eq!(density, fs::read(dir.join("out/density.bin")).unwrap());
    let standard = run(&base.replace("[scf]", "[scf]\nconvergence_profile='standard'"));
    assert_eq!(standard.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&standard.stderr).contains("strict"));
    fs::remove_dir_all(dir).unwrap();
}

#[test]
#[ignore = "release spin/FD envelope scans: cargo test --release -p dft-cli --test stress -- --include-ignored --nocapture"]
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
        let text = format!(
            "[system]\nname='synthetic stress'\nlattice_unit='bohr'\nlattice_vectors=[[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0]]\n[[system.atoms]]\nspecies='H'\nfractional=[0.0,0.0,0.0]\n[[system.atoms]]\nspecies='H'\nfractional=[0.27,0.23,0.19]\n[pseudopotentials]\nH='syn.xml'\n[basis]\necutwfc_ry=5.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pw'\n[scf]\nnbands=4\n{policy}{spin}[output]\ndirectory='{name}'\nstress=true\n"
        );
        fs::write(dir.join("input.toml"), text).unwrap();
        let out = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
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
