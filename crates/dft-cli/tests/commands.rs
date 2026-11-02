//! Exercise the public process boundary: output and nonzero error exits.
use std::{path::PathBuf, process::Command};

fn example() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples/si_lda.toml")
}

#[test]
fn example_commands_and_failures() {
    for command in ["validate", "grid"] {
        let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .arg(command)
            .arg(example())
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let text = String::from_utf8(output.stdout).unwrap();
        if command == "validate" {
            assert!(text.contains("Valid: Si diamond"));
            assert!(text.contains("file existence not checked"));
        } else {
            for label in [
                "Bohr^3",
                "Bohr^-1",
                "Ry",
                "Ha",
                "Plane-wave count:",
                "Density-sphere count:",
                "FFT dims:",
                "automatic",
            ] {
                assert!(text.contains(label), "missing {label}");
            }
        }
    }
    let missing = example().with_file_name("does-not-exist.toml");
    let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
        .arg("validate")
        .arg(&missing)
        .output()
        .unwrap();
    assert!(!output.status.success());
    let error = String::from_utf8(output.stderr).unwrap();
    assert!(error.contains("error:"));
    assert!(error.contains("does-not-exist.toml"));
    let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
        .arg("unsupported-command")
        .output()
        .unwrap();
    assert!(!output.status.success());
}
