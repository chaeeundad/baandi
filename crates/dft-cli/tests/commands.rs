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

#[test]
fn convert_units_and_usage_errors() {
    for (args, expected) in [
        (["30", "Ry", "Ha"], 15.),
        (["1", "GPa", "kbar"], 10.),
        (["-2", "Ha", "Ry"], -4.),
    ] {
        let out = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .arg("convert")
            .args(args)
            .output()
            .unwrap();
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        let text = String::from_utf8(out.stdout).unwrap();
        assert_eq!(text.trim().parse::<f64>().unwrap(), expected);
        assert_eq!(
            text.trim()
                .split('e')
                .next()
                .unwrap()
                .chars()
                .filter(char::is_ascii_digit)
                .count(),
            15
        );
    }
    for args in [
        ["30", "watts", "Ha"],
        ["30", "Ry", "nm"],
        ["NaN", "Ha", "Ry"],
        ["inf", "Ha", "Ry"],
        ["1e308", "J", "Ha"],
    ] {
        let out = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .arg("convert")
            .args(args)
            .output()
            .unwrap();
        assert_eq!(out.status.code(), Some(2));
    }
}

#[test]
fn unit_strings_resolve_to_legacy_input() {
    use dft_cli::config::Config;
    let old = std::fs::read_to_string(example()).unwrap();
    let new = old.replace("ecutwfc_ry = 30.0", "cutoff = '30 Ry'");
    let legacy = Config::parse(&old).unwrap();
    let unit = Config::parse(&new).unwrap();
    assert_eq!(unit.basis.ecutwfc_ry, legacy.basis.ecutwfc_ry);
    assert_eq!(
        unit.prepare().unwrap().wavefunctions.cutoff(),
        legacy.prepare().unwrap().wavefunctions.cutoff()
    );
    for change in [
        "cutoff = '30 unknown'",
        "ecutwfc_ry = 30.0\ncutoff = '30 Ry'",
    ] {
        assert!(Config::parse(&old.replace("ecutwfc_ry = 30.0", change)).is_err());
    }
    let new =
        format!("{new}\n[scf]\nenergy_tolerance = '2e-8 Ry'\neigen_tolerance_energy = '1e-8 Ha'\n");
    let config = Config::parse(&new).unwrap();
    let scf = config.scf.unwrap();
    assert_eq!(scf.energy_tolerance_ha, Some(1e-8));
    assert_eq!(scf.eigen_tolerance, Some(1e-8));
    for field in ["energy_tolerance_ha = 1e-8", "eigen_tolerance = 1e-8"] {
        assert!(Config::parse(&format!("{new}{field}\n")).is_err());
    }
    let unit_lattice = old.replace("lattice_unit = \"angstrom\"", "").replace(
        "lattice_vectors = [[0.0, 2.715, 2.715], [2.715, 0.0, 2.715], [2.715, 2.715, 0.0]]",
        "lattice_vectors_with_units = [['0 A', '2.715 A', '2.715 A'], ['2.715 A', '0 A', '2.715 A'], ['2.715 A', '2.715 A', '0 A']]");
    let cell = Config::parse(&unit_lattice).unwrap().cell().unwrap();
    let expected = legacy.cell().unwrap();
    assert!((cell.volume_bohr3() / expected.volume_bohr3() - 1.).abs() < 1e-14);
    assert!(
        Config::parse(&unit_lattice.replace("[system]", "[system]\nlattice_unit = 'bohr'"))
            .is_err()
    );
}
