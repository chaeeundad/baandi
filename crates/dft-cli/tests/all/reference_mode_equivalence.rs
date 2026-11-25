#[allow(dead_code)]
#[path = "../common/mod.rs"]
mod common;
use serde_json::Value;
use std::{fs, process::Command};

#[test]
fn small_synthetic_si_mode_matches_manual_energy_density_and_counts() {
    let dir = std::env::temp_dir().join(format!("reference-mode-si-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    // Smooth closed-shell atomic start avoids cutoff ringing and the resulting
    // positivity repair outside the symmetry projector's density sphere.
    let rc: f64 = 1.8;
    let mut xml = common::gaussian_with(2., rc, false, 400)
        .replace("SLA PZ", "PBE")
        .replace("element=\"Xx\"", "element=\"Si\"");
    let pp = dft_upf::parse_upf(&xml).unwrap();
    let radial = pp
        .mesh()
        .r()
        .iter()
        .map(|r| format!("{:.17e}", r * (-r * r / (rc * rc)).exp()))
        .collect::<Vec<_>>()
        .join(" ");
    let begin = xml.find("<PP_RHOATOM>").unwrap() + "<PP_RHOATOM>".len();
    let end = xml.find("</PP_RHOATOM>").unwrap();
    let charge = pp
        .mesh()
        .r()
        .iter()
        .map(|r| {
            format!(
                "{:.17e}",
                8. / (std::f64::consts::PI.sqrt() * rc.powi(3))
                    * r
                    * r
                    * (-r * r / (rc * rc)).exp()
            )
        })
        .collect::<Vec<_>>()
        .join(" ");
    xml.replace_range(begin..end, &charge);
    xml = xml.replace("number_of_wfc=\"0\"", "number_of_wfc=\"1\"")
        .replace("</UPF>", &format!("<PP_PSWFC><PP_CHI.1 l=\"0\" label=\"s\" occupation=\"2\">{radial}</PP_CHI.1></PP_PSWFC></UPF>"));
    fs::write(dir.join("syn.xml"), xml).unwrap();
    // This tiny cubic basis has only a few independent density modes; avoid
    // asking the unregularized mixer for a rank-deficient long history.
    let base = "[system]\nname='synthetic Si'\nlattice_unit='bohr'\nlattice_vectors=[[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0]]\n[[system.atoms]]\nspecies='Si'\nfractional=[0.0,0.0,0.0]\n[pseudopotentials]\nSi='syn.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[scf]\ndriver='reference'\nmixing_history=2\n[xc]\nfunctional='pbe'\n";
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
        assert_eq!(a["scf"]["driver"], "reference");
        assert_eq!(b["scf"]["driver"], "reference");
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

/// cargo test --release -p dft-cli --test all reference_mode_equivalence::material_mode_matches_manual -- --ignored --nocapture
#[test]
#[ignore = "external PseudoDojo SCFs; run in release mode"]
fn material_mode_matches_manual() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let nc = std::env::var_os("DFT_NC_DIR")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| root.join("../pseudo/pseudodojo/ex"));
    assert!(
        ["Si", "Mg", "O"]
            .iter()
            .all(|s| nc.join(format!("{s}.upf")).is_file()),
        "set DFT_NC_DIR to the external PseudoDojo directory"
    );
    let dir = std::env::temp_dir().join(format!("reference-mode-materials-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    for case in ["si-20", "si-30", "si-40", "mgo-40", "mgo-50", "mgo-60"] {
        // Historical material inputs request legacy-only floor/mixer/tolerances.
        // Keep their physical model, and explicitly exercise reference routing.
        let text = fs::read_to_string(root.join(format!("tests/reference/materials/{case}.toml")))
            .unwrap();
        let mut input: toml::Value = toml::from_str(&text).unwrap();
        input["xc"].as_table_mut().unwrap().remove("density_floor");
        let scf = input["scf"].as_table_mut().unwrap();
        for key in [
            "mixing_q0_inverse_bohr",
            "energy_tolerance_ha",
            "density_residual_tolerance",
        ] {
            scf.remove(key);
        }
        scf.insert("driver".into(), "reference".into());
        scf.insert("mixing".into(), "pulay".into());
        for (species, path) in input["pseudopotentials"].as_table_mut().unwrap() {
            *path = nc
                .join(format!("{species}.upf"))
                .canonicalize()
                .unwrap()
                .to_str()
                .unwrap()
                .into();
        }
        let base = toml::to_string(&input).unwrap();
        let variants = if matches!(case, "si-20" | "mgo-40") {
            2
        } else {
            1
        };
        let mut results = Vec::new();
        for variant in 0..variants {
            let text = if variant == 0 {
                format!("{base}\n[compat]\nmode='reference'\n")
            } else {
                format!(
                    "{}\n[pseudo]\nlocal_mesh='qe_like'\nradial_even_rule='simpson_quadratic_end'\n",
                    base.replace("[xc]", "[xc]\npw92_prefactor='original'")
                )
            };
            let path = dir.join(format!("{case}-{variant}.toml"));
            let out = dir.join(format!("{case}-{variant}"));
            fs::write(&path, text).unwrap();
            let run = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
                .env("BAANDI_REFERENCE_STRICT", "off")
                .arg("run")
                .arg(&path)
                .arg("--output-dir")
                .arg(&out)
                .output()
                .unwrap();
            assert!(
                run.status.success(),
                "{case}: {}",
                String::from_utf8_lossy(&run.stderr)
            );
            let report: Value =
                serde_json::from_slice(&fs::read(out.join("results.json")).unwrap()).unwrap();
            assert_eq!(report["scf"]["driver"], "reference");
            assert_eq!(report["convergence"]["converged"], true);
            if variant == 0 {
                assert_eq!(report["compat_mode"], "reference");
            }
            results.push((
                report["total_energy"].clone(),
                report["convergence"]["iterations"].clone(),
                fs::read(out.join("density.bin")).unwrap(),
            ));
        }
        if results.len() == 2 {
            assert_eq!(results[0], results[1], "{case}");
        }
    }
    fs::remove_dir_all(dir).unwrap();
}
