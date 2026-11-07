#[allow(dead_code)]
mod common;
use serde_json::Value;
use std::{fs, process::Command};

#[test]
fn config_fallback_results_and_explain() {
    let dir = std::env::temp_dir().join(format!("p46a-cli-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("syn.xml"),
        common::gaussian_with(2., 1., false, 101)
            .replace("SLA PZ", "SLA PW")
            .replace("element=\"Xx\"", "element=\"He\""),
    )
    .unwrap();
    let base = "[system]\nname='synthetic'\nlattice_unit='bohr'\nlattice_vectors=[[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0]]\n[[system.atoms]]\nspecies='He'\nfractional=[0.0,0.0,0.0]\n[pseudopotentials]\nHe='syn.xml'\n[basis]\necutwfc_ry=4.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pw'\n[scf]\nnbands=5\n";
    for (name, suffix, reason) in [
        ("missing", "", "PP_PSWFC"),
        (
            "spin",
            "[spin]\nmode='collinear'\nn_up=1\nn_down=1\n",
            "spin-polarized",
        ),
        (
            "fd",
            "[occupations]\nmode='fermi_dirac'\nsmearing_ha=0.01\nband_cap=5\n",
            "Fermi-Dirac",
        ),
    ] {
        let mut energies = Vec::new();
        for mode in ["free_electron", "atomic"] {
            let input = dir.join("input.toml");
            fs::write(
                &input,
                format!("{base}start_wavefunctions='{mode}'\n{suffix}"),
            )
            .unwrap();
            let out = dir.join(format!("{name}-{mode}"));
            let run = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
                .arg("run")
                .arg(&input)
                .arg("--output-dir")
                .arg(&out)
                .output()
                .unwrap();
            assert!(
                run.status.success(),
                "{name} {mode}: {}",
                String::from_utf8_lossy(&run.stderr)
            );
            let result: Value =
                serde_json::from_slice(&fs::read(out.join("results.json")).unwrap()).unwrap();
            energies.push(result["total_energy"]["ha"].clone());
            if mode == "atomic" {
                assert_eq!(result["scf"]["start_wavefunctions"], "atomic");
                assert!(
                    result["preflight"]["atomic_start"]["fallback_reason"]
                        .as_str()
                        .unwrap()
                        .contains(reason)
                );
                assert_eq!(
                    result["preflight"]["atomic_start"]["atomic_orbitals_used"],
                    0
                );
                let explain = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
                    .arg("explain")
                    .arg(&input)
                    .output()
                    .unwrap();
                assert!(explain.status.success());
                assert!(String::from_utf8_lossy(&explain.stderr).contains(reason));
            }
        }
        assert_eq!(energies[0], energies[1], "{name}");
    }
    let input = dir.join("input.toml");
    fs::write(&input, format!("{base}start_wavefunctions='typo'\n")).unwrap();
    let bad = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
        .arg("validate")
        .arg(input)
        .output()
        .unwrap();
    assert!(!bad.status.success());
    assert!(String::from_utf8_lossy(&bad.stderr).contains("start_wavefunctions"));
    fs::remove_dir_all(dir).unwrap();
}
