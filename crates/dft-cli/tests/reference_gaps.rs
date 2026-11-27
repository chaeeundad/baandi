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

/// Automatic FD band count on the reference driver: scf.nbands omitted grows
/// the count on coverage rejections (one stderr line each) and records it in
/// reference_scf.band_count_growth; explicit scf.nbands keeps the fixed count
/// (the cap) and reaches the same free energy.
#[test]
fn reference_fd_automatic_band_count() {
    let root = directory("fd-bands");
    let base = synthetic(&root)
        .replace("nbands=3\n", "")
        .replace(
            "functional='lda-pz'",
            "functional='lda-pw'\nallow_functional_override=true",
        )
        .replace(
            "kshift=[0.5,0.5,0.5]",
            "kshift=[0.5,0.5,0.5]\nktime_reversal=false",
        )
        + FD_GROWTH;
    let invoke = |name: &str, text: &str| {
        fs::write(root.join(format!("{name}.toml")), text).unwrap();
        let mut c = std::process::Command::new(env!("CARGO_BIN_EXE_dft-cli"));
        for (k, _) in std::env::vars().filter(|(k, _)| k.starts_with("BAANDI_")) {
            c.env_remove(k);
        }
        let output = c
            .current_dir(&root)
            .args(["run", &format!("{name}.toml"), "--output-dir", name])
            .output()
            .unwrap();
        let doc = fs::read(root.join(name).join("results.json"))
            .ok()
            .map(|b| serde_json::from_slice::<Value>(&b).unwrap());
        (output, doc)
    };
    let (output, doc) = invoke("auto", &base);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    let doc_auto = doc.unwrap();
    let growth = &doc_auto["reference_scf"]["band_count_growth"];
    assert_eq!(growth["mode"], "automatic", "{growth}");
    // GB2: the reference FD start is the atomic one (density and orbital
    // span); preflight no longer reports an FD fallback for this driver.
    assert_eq!(doc_auto["reference_scf"]["start"]["density"], "atomic");
    assert!(
        doc_auto["reference_scf"]["start"]["orbitals"]
            .as_str()
            .unwrap()
            .starts_with("atomic"),
        "{}",
        doc_auto["reference_scf"]["start"]
    );
    assert!(!stderr.contains("Fermi-Dirac occupations"), "{stderr}");
    assert_eq!(
        growth["growth_trigger"],
        "stopping, budget-final or verification map (GB2)"
    );
    assert_eq!(growth["initial_bands"], 5);
    assert_eq!(growth["band_cap"], 12);
    let events = growth["events"].as_array().unwrap();
    assert!(!events.is_empty(), "{growth}");
    assert_eq!(
        stderr
            .lines()
            .filter(|l| l.starts_with("band growth:"))
            .count(),
        events.len(),
        "{stderr}"
    );
    assert_eq!(growth["final_bands"], events.last().unwrap()["to_bands"]);
    assert_eq!(doc_auto["occupations"]["shape"][2], growth["final_bands"]);
    assert_eq!(doc_auto["convergence"]["converged"], true);
    // Explicit scf.nbands: fixed count = scf.nbands, no growth.
    let (output, doc) = invoke("explicit", &base.replace("[scf]\n", "[scf]\nnbands=12\n"));
    assert!(output.status.success());
    let doc_explicit = doc.unwrap();
    let fixed = &doc_explicit["reference_scf"]["band_count_growth"];
    assert_eq!(fixed["mode"], "fixed");
    assert_eq!(fixed["final_bands"], 12);
    assert!(!String::from_utf8_lossy(&output.stderr).contains("band growth:"));
    let f = |d: &Value| d["thermodynamics"]["free_energy_ha"].as_f64().unwrap();
    let mu = |d: &Value| {
        d["occupations"]["chemical_potential"][0]["mu_ha"]
            .as_f64()
            .unwrap()
    };
    eprintln!(
        "auto {growth}\nF auto {:.15} explicit {:.15} diff {:.3e}; mu diff {:.3e}",
        f(&doc_auto),
        f(&doc_explicit),
        f(&doc_auto) - f(&doc_explicit),
        mu(&doc_auto) - mu(&doc_explicit)
    );
    assert!((f(&doc_auto) - f(&doc_explicit)).abs() <= 5e-9);
    assert!((mu(&doc_auto) - mu(&doc_explicit)).abs() <= 1e-7);
    // Cap reached: the existing message plus the count reached.
    let (output, _) = invoke("capped", &base.replace("band_cap=12", "band_cap=6"));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(!output.status.success());
    assert!(
        stderr.contains("insufficient empty-state coverage at explicit cap")
            && stderr.contains("automatic band count reached 6 bands; occupations.band_cap = 6"),
        "{stderr}"
    );
    // Explicit scf.nbands below the needed count: the existing failure,
    // without growth (the automatic run needed more than 5 bands).
    let (output, _) = invoke(
        "explicit-short",
        &base.replace("[scf]\n", "[scf]\nnbands=5\n"),
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(!output.status.success());
    assert!(
        stderr.contains("insufficient empty-state coverage at explicit cap; increase cap and repeat the convergence study\")")
            && !stderr.contains("band growth:"),
        "{stderr}"
    );
    fs::remove_dir_all(root).unwrap();
}
const FD_GROWTH: &str = "[occupations]\nmode='fermi_dirac'\nsmearing_ha=0.02\nband_cap=12\n";

/// GB1: odd electron count (one synthetic atom with Z = 3, no atomic start
/// data, so the start is non-atomic), spin-unpolarized FD, neither
/// `scf.nbands` nor `band_cap`: the default cap max(32, start + two growth
/// steps) gives empty states, the automatic count runs; an explicit count
/// equal to the delivered automatic count runs to the same free energy; an
/// explicit count too small at the fixed point is refused at delivery with
/// the unchanged message.
#[test]
fn reference_fd_odd_electrons_default_cap() {
    let root = directory("fd-odd");
    fs::write(
        root.join("odd.xml"),
        common::gaussian_with(3., 1.8, true, 161).replace("element=\"Xx\"", "element=\"Li\""),
    )
    .unwrap();
    let base = "[system]\nname='GB1'\nlattice_unit='bohr'\nlattice_vectors=[[7.0,0.0,0.0],[0.2,7.0,0.0],[0.0,0.1,7.0]]\n[[system.atoms]]\nspecies='Li'\nfractional=[0.13,0.21,0.32]\n[pseudopotentials]\nLi='odd.xml'\n[basis]\necutwfc_ry=6.0\nkgrid=[2,2,2]\nkshift=[0.5,0.5,0.5]\nktime_reversal=false\n[xc]\nfunctional='lda-pw'\nallow_functional_override=true\n[scf]\nmax_iterations=150\nstandard_threshold_ha=5e-13\ndriver='reference'\n[occupations]\nmode='fermi_dirac'\nsmearing_ha=0.02\n";
    let invoke = |name: &str, text: &str| {
        fs::write(root.join(format!("{name}.toml")), text).unwrap();
        let mut c = std::process::Command::new(env!("CARGO_BIN_EXE_dft-cli"));
        for (k, _) in std::env::vars().filter(|(k, _)| k.starts_with("BAANDI_")) {
            c.env_remove(k);
        }
        let output = c
            .current_dir(&root)
            .args(["run", &format!("{name}.toml"), "--output-dir", name])
            .output()
            .unwrap();
        let doc = fs::read(root.join(name).join("results.json"))
            .ok()
            .map(|b| serde_json::from_slice::<Value>(&b).unwrap());
        (output, doc)
    };
    let (output, doc) = invoke("auto", base);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    let auto = doc.unwrap();
    let growth = &auto["reference_scf"]["band_count_growth"];
    assert_eq!(growth["mode"], "automatic", "{growth}");
    // n = 1.5: start max(2 + 4, ceil(1.8)) = 6; default cap max(32, 6->10->14).
    assert_eq!(growth["initial_bands"], 6);
    assert_eq!(growth["band_cap"], 32);
    assert_eq!(auto["convergence"]["converged"], true);
    let bands = growth["final_bands"].as_u64().unwrap();
    let explicit_input = base.replace("[scf]\n", &format!("[scf]\nnbands={bands}\n"));
    let (output, doc) = invoke("explicit", &explicit_input);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    let explicit = doc.unwrap();
    let fixed = &explicit["reference_scf"]["band_count_growth"];
    assert_eq!(fixed["mode"], "fixed");
    assert_eq!(fixed["final_bands"], bands);
    assert!(fixed["transient_coverage_misses"].is_array());
    let f = |d: &Value| d["thermodynamics"]["free_energy_ha"].as_f64().unwrap();
    eprintln!(
        "GB1 odd: auto {growth}\nF auto {:.15} explicit {:.15} diff {:.3e}; misses {}",
        f(&auto),
        f(&explicit),
        f(&auto) - f(&explicit),
        fixed["transient_coverage_misses"]
    );
    assert!((f(&auto) - f(&explicit)).abs() <= 5e-9);
    let (output, _) = invoke("short", &base.replace("[scf]\n", "[scf]\nnbands=2\n"));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(!output.status.success());
    assert!(
        stderr.contains("insufficient empty-state coverage at explicit cap; increase cap and repeat the convergence study\")")
            && stderr.contains("fd coverage: delivered")
            && !stderr.contains("band growth:"),
        "{stderr}"
    );
    fs::remove_dir_all(root).unwrap();
}
