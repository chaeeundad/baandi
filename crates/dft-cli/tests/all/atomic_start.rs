#[allow(dead_code)]
#[path = "../common/mod.rs"]
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
        // ONEPATH-3B: a collinear fixed sector no longer excludes the atomic
        // start; this fixture falls back only for its missing orbitals.
        (
            "spin",
            "[spin]\nmode='collinear'\nn_up=1\nn_down=1\n",
            "PP_PSWFC",
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

/// ONEPATH-3B: a collinear fixed sector with PP_PSWFC orbitals starts from the
/// atomic subspace; preflight and results.json record it (no fallback), the
/// span counted once, Rayleigh-Ritz once per channel.
#[test]
fn spin_fixed_sector_records_the_atomic_start() {
    let dir = std::env::temp_dir().join(format!("onepath3b-cli-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let text = common::gaussian_with(2., 1., false, 101)
        .replace("SLA PZ", "SLA PW")
        .replace("element=\"Xx\"", "element=\"He\"");
    let r = dft_upf::parse_upf(&text).unwrap().mesh().r().to_vec();
    let chi = r
        .iter()
        .map(|r| format!("{}", r * (-0.7 * r * r).exp()))
        .collect::<Vec<_>>()
        .join(" ");
    fs::write(
        dir.join("syn.xml"),
        text.replace("number_of_wfc=\"0\"", "number_of_wfc=\"1\"")
            .replace(
                "</UPF>",
                &format!("<PP_PSWFC><PP_CHI.1 l=\"0\" label=\"s\" occupation=\"2\">{chi}</PP_CHI.1></PP_PSWFC></UPF>"),
            ),
    )
    .unwrap();
    let input = dir.join("input.toml");
    fs::write(
        &input,
        "[system]\nname='synthetic'\nlattice_unit='bohr'\nlattice_vectors=[[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0]]\n[[system.atoms]]\nspecies='He'\nfractional=[0.0,0.0,0.0]\n[pseudopotentials]\nHe='syn.xml'\n[basis]\necutwfc_ry=4.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pw'\n[scf]\nnbands=5\n[spin]\nmode='collinear'\nn_up=1\nn_down=1\n",
    )
    .unwrap();
    let out = dir.join("spin-atomic");
    let run = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
        .arg("run")
        .arg(&input)
        .arg("--output-dir")
        .arg(&out)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let result: Value =
        serde_json::from_slice(&fs::read(out.join("results.json")).unwrap()).unwrap();
    assert_eq!(result["convergence"]["converged"], true);
    assert_eq!(
        result["preflight"]["atomic_start"]["fallback_reason"],
        Value::Null
    );
    let a = &result["atomic_start"];
    assert_eq!(a["requested"], "atomic");
    assert_eq!(a["fallback_reason"], Value::Null);
    assert_eq!(a["atomic_orbitals_used"], 1);
    // n_up == n_down: with BAANDI_SPIN_SHARE on (default since 2026-11-17) the two identical
    // channels share one atomic-start EVD (bit-identical to solving them separately).
    assert_eq!(a["evd_orders"].as_array().unwrap().len(), 1);
    fs::remove_dir_all(dir).unwrap();
}
