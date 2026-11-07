#[allow(dead_code)]
#[path = "common/mod.rs"]
mod common;
use serde_json::Value;
use std::{fs, process::Command};

#[test]
fn small_synthetic_si_mode_matches_manual_energy_density_and_counts() {
    let dir = std::env::temp_dir().join(format!("reference-mode-si-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("syn.xml"),
        common::gaussian_with(4., 1., false, 400)
            .replace("SLA PZ", "PBE")
            .replace("element=\"Xx\"", "element=\"Si\""),
    )
    .unwrap();
    let base = "[system]\nname='synthetic Si'\nlattice_unit='bohr'\nlattice_vectors=[[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0]]\n[[system.atoms]]\nspecies='Si'\nfractional=[0.0,0.0,0.0]\n[pseudopotentials]\nSi='syn.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='pbe'\n";
    let manual = "pw92_prefactor='original'\n[pseudo]\nlocal_mesh='qe_like'\nradial_even_rule='simpson_quadratic_end'\n";
    for symmetry in [false, true] {
        let suffix = if symmetry {
            "[symmetry]\nmode='auto'\n"
        } else {
            ""
        };
        for (name, extra) in [
            ("mode", "[compat]\nmode='reference'\n".to_string()),
            (
                "manual",
                format!(
                    "{manual}{}",
                    if symmetry {
                        "[symmetry]\nmode='auto'\ndensity_group='crystal'\n"
                    } else {
                        ""
                    }
                ),
            ),
        ] {
            let inp = dir.join(format!("{name}.toml"));
            fs::write(
                &inp,
                format!("{base}{extra}{}", if name == "mode" { suffix } else { "" }),
            )
            .unwrap();
            let run = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
                .args(["run"])
                .arg(&inp)
                .arg("--output-dir")
                .arg(dir.join(name))
                .output()
                .unwrap();
            assert!(
                run.status.success(),
                "{}",
                String::from_utf8_lossy(&run.stderr)
            );
            if name == "mode" {
                assert!(String::from_utf8_lossy(&run.stderr).contains("compat_mode"));
                let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
                    .arg("validate")
                    .arg(&inp)
                    .arg("--json")
                    .output()
                    .unwrap();
                assert!(
                    output.status.success(),
                    "{}",
                    String::from_utf8_lossy(&output.stderr)
                );
                let report: Value = serde_json::from_slice(&output.stdout).unwrap();
                assert_eq!(report["compat_mode"], "reference");
                assert_eq!(
                    report["density_group"],
                    if symmetry { "crystal" } else { "mesh" }
                );
                assert!(
                    report["scope"]
                        .as_str()
                        .unwrap()
                        .contains("unpolarized PBE energies and densities")
                );
                let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
                    .arg("explain")
                    .arg(&inp)
                    .output()
                    .unwrap();
                assert!(output.status.success());
                assert!(String::from_utf8_lossy(&output.stdout).contains("compat_mode"));
                assert!(
                    String::from_utf8_lossy(&output.stderr)
                        .contains("unpolarized PBE energies and densities")
                );
            }
        }
        let a: Value =
            serde_json::from_slice(&fs::read(dir.join("mode/results.json")).unwrap()).unwrap();
        let b: Value =
            serde_json::from_slice(&fs::read(dir.join("manual/results.json")).unwrap()).unwrap();
        assert_eq!(a["compat_mode"], "reference");
        for key in [
            "pw92_prefactor",
            "local_mesh",
            "radial_even_rule",
            "density_group",
        ] {
            assert!(a["compat"][key].is_string(), "{key}");
        }
        assert_eq!(a["total_energy"], b["total_energy"]);
        assert_eq!(
            a["convergence"]["iterations"],
            b["convergence"]["iterations"]
        );
        assert_eq!(a["convergence"]["converged"], true);
        assert_eq!(
            fs::read(dir.join("mode/density.bin")).unwrap(),
            fs::read(dir.join("manual/density.bin")).unwrap()
        );
    }
    fs::remove_dir_all(dir).unwrap();
}

/// cargo test --release -p dft-cli --test reference_mode_equivalence material_mode_matches_manual -- --ignored --nocapture
#[test]
#[ignore = "external PseudoDojo SCFs; run in release mode"]
fn material_mode_matches_manual() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let status = Command::new("python3")
        .args(["-I"])
        .arg(root.join("crates/dft-engine/tests/reference_mode_materials.py"))
        .arg("--binary")
        .arg(env!("CARGO_BIN_EXE_dft-cli"))
        .arg("--evidence")
        .arg(root.join("P44A-evidence/materials"))
        .status()
        .unwrap();
    assert!(status.success());
}
