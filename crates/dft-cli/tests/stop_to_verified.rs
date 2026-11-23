//! Independent wall-interval fixtures; seconds here are synthetic, not timings.
#[allow(dead_code)]
#[path = "../src/refscf_check.rs"]
mod reader;
use serde_json::{Value, json};

fn record() -> Value {
    let maps: Vec<_> = (1..=12)
        .map(|i| {
            let start = 10. * (i - 1) as f64;
            json!({"iteration":i,"dr2_ry":if i>=10 {1e-9} else {1.},
            "map_start_s":start,"estimator_s":start+2.,
            "verification_start_s":if i>=10 {Some(start+3.)} else {None},
            "verification_end_s":if i>=10 {Some(start+5.)} else {None},
            "map_end_s":start+8.})
        })
        .collect();
    json!({"run":{"completed":true},"reference_scf":{
        "stop_to_verified_clock":"reference-scf-wall-1","threshold_ry":1e-8,
        "iterations":maps,"verification":{"mode":"fixed_point","passed":true,
            "attempts":[{"iteration":10,"passed":false},{"iteration":11,"passed":false},
                {"iteration":12,"passed":true}]}}})
}

#[test]
fn stop_to_verified_three_attempt_partition() {
    let summary = reader::stop_to_verified_summary(&record()).unwrap();
    assert_eq!(summary["first_stop_map"], 10);
    assert_eq!(summary["passing_map"], 12);
    assert_eq!(summary["attempt_count"], 3);
    // [92,115]: attempts [93,95],[103,105],[113,115] = 6;
    // recovery [95,98],[100,103],[105,108],[110,113] = 12;
    // initial glue [92,93] and inter-map gaps [98,100],[108,110] = 5.
    assert_eq!(summary["interval_wall_s"], 23.);
    assert_eq!(summary["verification_attempts_s"], 6.);
    assert_eq!(summary["ordinary_scf_after_failure_s"], 12.);
    assert_eq!(summary["unassigned_s"], 5.);
    assert_eq!(summary["closure_error_s"], 0.);
    assert_eq!(summary["parts_sum_checked"], true);
}

#[test]
fn stop_to_verified_first_pass_has_zero_recovery() {
    let mut doc = record();
    doc["reference_scf"]["iterations"]
        .as_array_mut()
        .unwrap()
        .truncate(10);
    let attempts = doc["reference_scf"]["verification"]["attempts"]
        .as_array_mut()
        .unwrap();
    attempts.truncate(1);
    attempts[0]["passed"] = json!(true);
    let summary = reader::stop_to_verified_summary(&doc).unwrap();
    assert_eq!(summary["passing_map"], 10);
    assert_eq!(summary["attempt_count"], 1);
    assert_eq!(summary["interval_wall_s"], 3.);
    assert_eq!(summary["verification_attempts_s"], 2.);
    assert_eq!(summary["ordinary_scf_after_failure_s"], 0.);
    assert_eq!(summary["unassigned_s"], 1.);
}

#[test]
fn stop_to_verified_unverified_and_no_stop_are_explicit() {
    let mut doc = record();
    doc["reference_scf"]["verification"]["passed"] = json!(false);
    doc["reference_scf"]["verification"]["attempts"][2]["passed"] = json!(false);
    let summary = reader::stop_to_verified_summary(&doc).unwrap();
    assert_eq!(summary["status"], "UNVERIFIED");
    assert!(summary["passing_map"].is_null());
    assert_eq!(summary["interval_wall_s"], 26.);
    assert_eq!(summary["ordinary_scf_after_failure_s"], 15.);
    assert_eq!(summary["unassigned_s"], 5.);
    // A non-stopping map after the failed attempt is ordinary recovery too.
    doc["reference_scf"]["iterations"][11]["dr2_ry"] = json!(1.);
    doc["reference_scf"]["iterations"][11]["verification_start_s"] = Value::Null;
    doc["reference_scf"]["iterations"][11]["verification_end_s"] = Value::Null;
    doc["reference_scf"]["verification"]["attempts"]
        .as_array_mut()
        .unwrap()
        .truncate(2);
    let summary = reader::stop_to_verified_summary(&doc).unwrap();
    assert_eq!(summary["attempt_count"], 2);
    assert_eq!(summary["ordinary_scf_after_failure_s"], 17.);
    assert_eq!(summary["verification_attempts_s"], 4.);
    assert_eq!(summary["unassigned_s"], 5.);
    doc["reference_scf"]["iterations"]
        .as_array_mut()
        .unwrap()
        .truncate(9);
    doc["reference_scf"]["verification"]["attempts"] = json!([]);
    let summary = reader::stop_to_verified_summary(&doc).unwrap();
    assert_eq!(summary["status"], "UNVERIFIED");
    assert!(summary["first_stop_map"].is_null());
    assert!(summary["interval_wall_s"].is_null());
    assert_eq!(summary["attempt_count"], 0);
    assert_eq!(summary["parts_sum_checked"], false);
}

#[test]
fn stop_to_verified_verification_off_never_implies_verified() {
    let mut doc = record();
    doc["reference_scf"]["verification"] = Value::Null;
    let maps = doc["reference_scf"]["iterations"].as_array_mut().unwrap();
    maps.truncate(10);
    for map in maps {
        map["verification_start_s"] = Value::Null;
        map["verification_end_s"] = Value::Null;
    }
    let summary = reader::stop_to_verified_summary(&doc).unwrap();
    assert_eq!(summary["status"], "UNVERIFIED");
    assert!(summary["passing_map"].is_null());
    assert_eq!(summary["attempt_count"], 0);
    assert_eq!(summary["interval_wall_s"], 6.);
    assert_eq!(summary["verification_attempts_s"], 0.);
    assert_eq!(summary["ordinary_scf_after_failure_s"], 0.);
    assert_eq!(summary["unassigned_s"], 6.);
}

#[test]
fn stop_to_verified_rejects_bad_evidence() {
    let mut bad = Vec::new();
    let mut doc = record();
    doc["run"]["completed"] = json!(false);
    bad.push(doc);
    let mut doc = record();
    doc["reference_scf"]["stop_to_verified_clock"] = json!("old");
    bad.push(doc);
    let mut doc = record();
    doc["reference_scf"]["iterations"][10]["map_start_s"] = json!(97.);
    bad.push(doc);
    let mut doc = record();
    doc["reference_scf"]["iterations"][9]["verification_end_s"] = json!(99.);
    bad.push(doc);
    let mut doc = record();
    doc["reference_scf"]["iterations"][9]["verification_start_s"] = Value::Null;
    bad.push(doc);
    let mut doc = record();
    doc["reference_scf"]["iterations"][9]["estimator_s"] = json!(-1.);
    bad.push(doc);
    let mut doc = record();
    doc["reference_scf"]["iterations"][9]["map_end_s"] = Value::Null;
    bad.push(doc);
    let mut doc = record();
    doc["reference_scf"]["verification"]["attempts"][0]["iteration"] = json!(11);
    bad.push(doc);
    let mut doc = record();
    doc["reference_scf"]["verification"]["attempts"][0]["passed"] = json!(true);
    bad.push(doc);
    let mut doc = record();
    doc["reference_scf"]["verification"] = Value::Null;
    bad.push(doc);
    for (i, doc) in bad.iter().enumerate() {
        assert!(
            reader::stop_to_verified_summary(doc).is_err(),
            "bad evidence {i}"
        );
    }
}

#[allow(dead_code)]
mod common;

#[test]
fn stop_to_verified_cli_records_at_both_diagnostic_levels() {
    use std::{fs, process::Command};
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(format!(
        "../../target/tmp/stop-to-verified-cli-{}",
        std::process::id()
    ));
    fs::create_dir_all(&dir).unwrap();
    let mut pp =
        common::gaussian_with(2., 1.8, false, 161).replace("element=\"Xx\"", "element=\"He\"");
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
    pp = pp.replace("number_of_wfc=\"0\"","number_of_wfc=\"1\"").replace("</UPF>",&format!("<PP_PSWFC><PP_CHI.1 l=\"0\" label=\"s\" occupation=\"2\">{}</PP_CHI.1></PP_PSWFC></UPF>",array(&|r|r*(-0.6*r*r).exp())));
    fs::write(dir.join("syn.xml"), pp).unwrap();
    let text = "[system]\nname='reference test'\nlattice_unit='bohr'\nlattice_vectors=[[7.0,0.0,0.0],[0.2,7.0,0.0],[0.0,0.1,7.0]]\n[[system.atoms]]\nspecies='He'\nfractional=[0.13,0.21,0.32]\n[pseudopotentials]\nHe='syn.xml'\n[basis]\necutwfc_ry=4.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[scf]\nnbands=3\nmax_iterations=100\nmixing='broyden'\nmixing_beta=0.7\nmixing_history=8\n";

    for (level, profile) in [
        ("lightweight", ""),
        ("full", "[profile]\ndavidson_stages=true\n"),
    ] {
        fs::write(
            dir.join("in.toml"),
            format!("{text}driver='reference'\n{profile}"),
        )
        .unwrap();
        let mut command = Command::new(env!("CARGO_BIN_EXE_dft-cli"));
        for (key, _) in std::env::vars_os() {
            if key.to_string_lossy().starts_with("BAANDI_") {
                command.env_remove(key);
            }
        }
        let output = command
            .current_dir(&dir)
            .args(["run", "in.toml", "--output-dir", level])
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let path = dir.join(level).join("results.json");
        let doc: Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
        let r = &doc["reference_scf"];
        assert_eq!(r["diagnostics"]["level"], level);
        assert_eq!(r["stop_to_verified_clock"], "reference-scf-wall-1");
        let expected = reader::stop_to_verified_summary(&doc).unwrap();
        assert_eq!(expected["status"], "VERIFIED");
        for (map, history) in r["iterations"]
            .as_array()
            .unwrap()
            .iter()
            .zip(doc["convergence"]["history"].as_array().unwrap())
        {
            for key in [
                "map_start_s",
                "estimator_s",
                "verification_start_s",
                "verification_end_s",
                "map_end_s",
            ] {
                assert!(map.get(key).is_some(), "{level}: missing {key}");
                assert_eq!(map[key], history["reference_scf"][key]);
            }
        }
        let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .args(["refscf-check", "--stop-to-verified"])
            .arg(&path)
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let actual: Value = serde_json::from_slice(&output.stdout).unwrap();
        assert_eq!(actual, expected);
        // Completed JSON must bypass TOML input parsing, while QE replay keeps its required flag.
        let missing_mode = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .arg("refscf-check")
            .arg(&path)
            .output()
            .unwrap();
        assert!(!missing_mode.status.success());
    }
    fs::remove_dir_all(dir).unwrap();
}
