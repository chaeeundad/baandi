//! G5B: executable output contracts and the fixed-cell k-mesh trajectory oracle.
#[allow(dead_code)]
#[path = "common/mod.rs"]
mod common;
#[allow(dead_code)]
#[path = "../src/refscf_check.rs"]
mod reader;
use serde_json::{Value, json};
use std::{
    fs,
    path::{Path, PathBuf},
};

fn directory(name: &str) -> PathBuf {
    let path = std::env::temp_dir().join(format!("g5b-{name}-{}", std::process::id()));
    fs::create_dir_all(&path).unwrap();
    path
}
fn run(root: &Path, name: &str, text: &str) -> Value {
    let input = root.join(format!("{name}.toml"));
    let out = root.join(name);
    fs::write(&input, text).unwrap();
    dft_cli::run::execute(&input, Some(&out)).unwrap();
    serde_json::from_slice(&fs::read(out.join("results.json")).unwrap()).unwrap()
}
fn synthetic(root: &Path) -> String {
    let mut pp =
        common::gaussian_with(2., 1.8, true, 161).replace("element=\"Xx\"", "element=\"He\"");
    let parsed = dft_upf::parse_upf(&pp).unwrap();
    let array = |f: &dyn Fn(f64) -> f64| {
        parsed
            .mesh()
            .r()
            .iter()
            .map(|&r| format!("{:.17e}", f(r)))
            .collect::<Vec<_>>()
            .join(" ")
    };
    let start = pp.find("<PP_RHOATOM>").unwrap() + "<PP_RHOATOM>".len();
    let end = pp.find("</PP_RHOATOM>").unwrap();
    pp.replace_range(
        start..end,
        &array(&|r| {
            8. / (std::f64::consts::PI.sqrt() * 1.8_f64.powi(3))
                * r
                * r
                * (-(r / 1.8).powi(2)).exp()
        }),
    );
    pp = pp.replace("number_of_wfc=\"0\"", "number_of_wfc=\"1\"").replace("</UPF>",
        &format!("<PP_PSWFC><PP_CHI.1 l=\"0\" label=\"s\" occupation=\"2\">{}</PP_CHI.1></PP_PSWFC></UPF>", array(&|r| r * (-0.6 * r * r).exp())));
    fs::write(root.join("syn.xml"), pp).unwrap();
    "[system]\nname='G5B'\nlattice_unit='bohr'\nlattice_vectors=[[7.0,0.0,0.0],[0.2,7.0,0.0],[0.0,0.1,7.0]]\n[[system.atoms]]\nspecies='He'\nfractional=[0.13,0.21,0.32]\n[pseudopotentials]\nHe='syn.xml'\n[basis]\necutwfc_ry=2.0\nkgrid=[2,2,2]\nkshift=[0.5,0.5,0.5]\n[xc]\nfunctional='lda-pz'\n[scf]\nnbands=3\nmax_iterations=150\nstandard_threshold_ha=5e-13\ndriver='reference'\n".into()
}
#[test]
fn reference_verify_off_labels_and_stop_to_verified_k_fd_spin_outputs() {
    let root = directory("outputs");
    let base = synthetic(&root);
    let mut docs = Vec::new();
    for mode in ["fixed_point", "off"] {
        docs.push(run(
            &root,
            mode,
            &format!(
                "{base}verify='{mode}'\n[output]\nforces=true\nstress=true\n[output.volumetric]\nquantities=['density']\nformats=['cube']\n"
            ),
        ));
    }
    let on = &docs[0];
    let off = &docs[1];
    assert_eq!(on["convergence"]["converged"], true);
    assert_eq!(off["convergence"]["converged"], true);
    for key in ["total_energy", "eigenvalues", "occupations", "kpoints"] {
        assert_eq!(on[key], off[key], "{key}");
    }
    for name in ["forces", "stress"] {
        assert_eq!(off[name]["verification"], "unverified");
        assert_eq!(
            off[name]["preconditions"]["fixed_point_verification_passed"],
            false
        );
        assert_eq!(on[name]["values"], off[name]["values"]);
    }
    assert_eq!(
        fs::read(root.join("fixed_point/density.cube")).unwrap(),
        fs::read(root.join("off/density.cube")).unwrap()
    );
    assert!(off["verification"].is_null());
    assert!(
        off["reference_scf"]["run"]["verification_calls"]
            .as_array()
            .unwrap()
            .is_empty()
    );
    assert_eq!(
        on["reference_scf"]["run"]["verification_calls"]
            .as_array()
            .unwrap()
            .len(),
        on["kpoints"]["count"].as_u64().unwrap() as usize
    );
    assert_eq!(on["verification"]["attempts"].as_array().unwrap().len(), 1);
    assert_eq!(
        on["convergence"]["iterations"],
        off["convergence"]["iterations"]
    );
    for (doc, expected) in [(&docs[0], "VERIFIED"), (&docs[1], "UNVERIFIED")] {
        assert_eq!(
            reader::stop_to_verified_summary(doc).unwrap()["status"],
            expected
        );
    }
    // Real FD and fixed-sector spin results have verdicts and durations, but
    // do not have the Gamma/k driver's monotonic map clocks.
    for (name, extra) in [
        (
            "fd",
            "[occupations]\nmode='fermi_dirac'\nsmearing_ha=0.01\nband_cap=12\n",
        ),
        ("spin", "[spin]\nmode='collinear'\nn_up=1\nn_down=1\n"),
    ] {
        // FD-FORCES: scf.nbands now sets the FD band count; this fixture keeps
        // its historical FD count, the cap (12).
        let base = if name == "fd" {
            base.replace("nbands=3\n", "")
        } else {
            base.clone()
        };
        let text = base
            .replace(
                "functional='lda-pz'",
                "functional='lda-pw'\nallow_functional_override=true",
            )
            .replace(
                "kshift=[0.5,0.5,0.5]",
                "kshift=[0.5,0.5,0.5]\nktime_reversal=false",
            );
        let doc = run(&root, name, &(text + extra));
        let summary = reader::stop_to_verified_summary(&doc).unwrap();
        assert_eq!(summary["status"], "VERIFIED");
        assert!(summary["interval_wall_s"].is_null());
        assert_eq!(summary["parts_sum_checked"], false);
        let mut invalid = doc.clone();
        invalid["verification"]["attempts"][0]["iteration"] = json!(0);
        assert!(reader::stop_to_verified_summary(&invalid).is_err());
        invalid = doc.clone();
        invalid["verified"] = json!(false);
        assert!(reader::stop_to_verified_summary(&invalid).is_err());
        let output = std::process::Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .arg("refscf-check")
            .arg("--stop-to-verified")
            .arg(root.join(name).join("results.json"))
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
    // Ionic relaxation retains its verification precondition.
    assert!(
        dft_cli::config::Config::parse(&(base + "verify='off'\n[relax]\nmax_steps=3\n")).is_err()
    );
}

#[test]
#[ignore = "external Si2 trajectory; DFT_SI_UPF=/path/Si.pz-vbc.UPF cargo test --release -p dft-cli --test reference_gaps reference_k_relax -- --ignored --nocapture"]
fn reference_k_relax_three_steps_matches_strict_oracle() {
    let root = directory("si2-relax");
    let pp = std::env::var("DFT_SI_UPF").unwrap();
    let mut frames = Vec::new();
    for driver in ["strict_oracle", "reference"] {
        let switches = dft_numerics::switches::RunConfig::from_pairs(&[(
            "reference_strict",
            if driver == "strict_oracle" {
                "on"
            } else {
                "off"
            },
        )])
        .unwrap();
        let _scope = switches.enter();
        let input = root.join(format!("{driver}.toml"));
        fs::write(
            &input,
            format!(
                r#"
[system]
name='G5B Si2 mesh relaxation'
lattice_unit='bohr'
lattice_vectors=[[0.0,5.13,5.13],[5.13,0.0,5.13],[5.13,5.13,0.0]]
[[system.atoms]]
species='Si'
fractional=[0.0,0.0,0.0]
[[system.atoms]]
species='Si'
fractional=[0.253,0.25,0.25]
[pseudopotentials]
Si={pp:?}
[basis]
ecutwfc_ry=8.0
kgrid=[2,2,2]
ktime_reversal=false
[xc]
functional='lda-pz'
[scf]
driver='reference'
nbands=5
max_iterations=200
standard_threshold_ha=5e-13
start_density='atomic'
start_wavefunctions='atomic'
[relax]
max_steps=3
force_tol_ha_bohr=1e-12
energy_tol_ha=1e-12
"#
            ),
        )
        .unwrap();
        let out = root.join(driver);
        let result = dft_cli::relax::execute(&input, &out, false);
        assert!(
            result.is_ok() || result.as_ref().is_err_and(|e| e.code == 3),
            "{result:?}"
        );
        let doc: Value =
            serde_json::from_slice(&fs::read(out.join("results.json")).unwrap()).unwrap();
        let trajectory = doc["relaxation"]["frames"].as_array().unwrap().clone();
        assert_eq!(trajectory.len(), 4);
        if driver == "reference" {
            for frame in &trajectory[1..] {
                assert_eq!(frame["provenance"], "warm_transfer");
                assert!(frame["start_fallback"].is_null());
            }
        }
        frames.push(trajectory);
    }
    for (step, (a, b)) in frames[0].iter().zip(&frames[1]).enumerate() {
        let de = (a["energy_ha"].as_f64().unwrap() - b["energy_ha"].as_f64().unwrap()).abs();
        let flatten = |v: &Value| -> Vec<f64> {
            v.as_array()
                .unwrap()
                .iter()
                .flat_map(|row| row.as_array().unwrap().iter().map(|x| x.as_f64().unwrap()))
                .collect()
        };
        let df = flatten(&a["forces_ha_bohr"])
            .iter()
            .zip(flatten(&b["forces_ha_bohr"]))
            .map(|(x, y)| (x - y).abs())
            .fold(0., f64::max);
        let dr = a["positions"]
            .as_array()
            .unwrap()
            .iter()
            .zip(b["positions"].as_array().unwrap())
            .map(|(x, y)| (x.as_f64().unwrap() - y.as_f64().unwrap()).abs())
            .fold(0., f64::max);
        println!("G5B k222 step={step} dE={de:e} Ha dF={df:e} Ha/Bohr dR={dr:e} Bohr");
        assert!(de <= 5e-9 && df <= 1e-6 && dr <= 1e-6);
    }
}
