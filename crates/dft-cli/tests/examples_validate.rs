//! Validate-only coverage of shipped workflows; never run SCF in debug tests.
#[allow(dead_code)]
mod common;

use std::{fs, path::PathBuf, process::Command};

struct Scratch(PathBuf);
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

#[test]
fn every_workflow_passes_preflight_with_synthetic_potentials() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let mut inputs: Vec<_> = fs::read_dir(root.join("examples"))
        .unwrap()
        .map(|entry| entry.unwrap().path().join("input.toml"))
        .filter(|path| path.is_file())
        .collect();
    inputs.sort();
    assert!(inputs.len() >= 13, "a tutorial input is missing");
    let scratch =
        Scratch(std::env::temp_dir().join(format!("baandi-examples-{}", std::process::id())));
    fs::create_dir(&scratch.0).unwrap();
    for input in inputs {
        let workflow = input.parent().unwrap().file_name().unwrap();
        let dir = scratch.0.join(workflow);
        fs::create_dir(&dir).unwrap();
        assert!(input.with_file_name("README.md").is_file());
        let mut config: toml::Value = toml::from_str(&fs::read_to_string(&input).unwrap()).unwrap();
        let functional = match config["xc"]["functional"].as_str().unwrap() {
            "lda-pz" => "SLA PZ",
            "lda-pw" => "SLA PW",
            "pbe" => "PBE",
            other => panic!("add synthetic functional coverage for {other}"),
        };
        for (species, entry) in config["pseudopotentials"].as_table_mut().unwrap() {
            // Match real-file valence/element/XC without bypassing preflight's
            // identity or occupation checks. These are not material potentials.
            let valence = match species.as_str() {
                "Si" => 4.0,
                "Al" => 3.0,
                "Fe" => 16.0,
                other => panic!("add synthetic valence coverage for {other}"),
            };
            let xml = common::gaussian_with(valence, 1.0, true, 101)
                .replace("element=\"Xx\"", &format!("element=\"{species}\""))
                .replace(
                    "functional=\"SLA PZ\"",
                    &format!("functional=\"{functional}\""),
                );
            let filename = format!("{species}.xml");
            fs::write(dir.join(&filename), xml).unwrap();
            *entry = toml::Value::String(filename);
        }
        let prepared = dir.join("input.toml");
        fs::write(&prepared, toml::to_string(&config).unwrap()).unwrap();
        let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .env("BAANDI_REFERENCE_STRICT", "off")
            .current_dir(&dir)
            .arg("validate")
            .arg(&prepared)
            .arg("--json")
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}: {}\n{}",
            input.display(),
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
        assert_eq!(report["valid"], true, "{}", input.display());
        assert_eq!(report["overrides"], serde_json::json!([]));
        assert!(!dir.join("results").exists(), "validate must not run SCF");
    }
}
