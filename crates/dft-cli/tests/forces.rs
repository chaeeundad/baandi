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
    let standard = base.replace("[scf]", "[scf]\nconvergence_profile='standard'");
    fs::write(dir.join("input.toml"), standard).unwrap();
    for action in ["validate", "run"] {
        let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .current_dir(&dir)
            .args([action, "input.toml"])
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(1));
        assert!(String::from_utf8_lossy(&output.stderr).contains("output.forces"));
    }
    fs::remove_dir_all(dir).unwrap();
}
