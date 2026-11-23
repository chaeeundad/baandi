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
    cli_roundtrip(false);
}

#[test]
fn stop_to_verified_cli_skips_and_reuses_at_both_diagnostic_levels() {
    cli_roundtrip(true);
}

fn cli_roundtrip(active: bool) {
    use std::{fs, process::Command};
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(format!(
        "../../target/tmp/stop-to-verified-cli-{}-{active}",
        std::process::id()
    ));
    fs::create_dir_all(&dir).unwrap();
    let (z, rc) = if active { (4., 1.2_f64) } else { (2., 1.8_f64) };
    let mut pp =
        common::gaussian_with(z, rc, active, 161).replace("element=\"Xx\"", "element=\"He\"");
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
            4. * z / (std::f64::consts::PI.sqrt() * rc.powi(3)) * r * r * (-(r / rc).powi(2)).exp()
        }),
    );
    pp = pp.replace("number_of_wfc=\"0\"","number_of_wfc=\"1\"").replace("</UPF>",&format!("<PP_PSWFC><PP_CHI.1 l=\"0\" label=\"s\" occupation=\"2\">{}</PP_CHI.1></PP_PSWFC></UPF>",array(&|r|r*(-0.6*r*r).exp())));
    if active {
        pp = pp
            .replace("number_of_wfc=\"1\"", "number_of_wfc=\"2\"")
            .replace(
                "</PP_PSWFC>",
                &format!(
                    "<PP_CHI.2 l=\"1\" label=\"p\" occupation=\"2\">{}</PP_CHI.2></PP_PSWFC>",
                    array(&|r| r * r * (-0.5 * r * r).exp())
                ),
            );
    }
    fs::write(dir.join("syn.xml"), pp).unwrap();
    let text = "[system]\nname='reference test'\nlattice_unit='bohr'\nlattice_vectors=[[7.0,0.0,0.0],[0.2,7.0,0.0],[0.0,0.1,7.0]]\n[[system.atoms]]\nspecies='He'\nfractional=[0.13,0.21,0.32]\n[pseudopotentials]\nHe='syn.xml'\n[basis]\necutwfc_ry=4.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[scf]\nnbands=3\nmax_iterations=100\nmixing='broyden'\nmixing_beta=0.7\nmixing_history=8\n";

    let text = if active {
        text.replace("7.0", "11.0")
            .replace("mixing_history=8", "mixing_history=4")
    } else {
        text.to_string()
    };
    for (variant, reuse, gate, first) in [
        ("off", "off", "off", "off"),
        ("default", "off", "on", "off"),
        ("first", "off", "on", "on"),
        ("first_inert", "off", "off", "on"),
        ("reuse", "on", "off", "off"),
        ("gate", "off", "on", "off"),
        ("both", "on", "on", "off"),
    ] {
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
            if variant != "default" {
                command
                    .env("BAANDI_REFSCF_VERIFY_REUSE", reuse)
                    .env("BAANDI_REFSCF_VERIFY_GATE", gate)
                    .env("BAANDI_REFSCF_VERIFY_GATE_FIRST", first);
            }
            let output_dir = format!("{level}-{variant}");
            let output = command
                .current_dir(&dir)
                .args(["run", "in.toml", "--output-dir", &output_dir])
                .output()
                .unwrap();
            assert!(
                output.status.success(),
                "{}",
                String::from_utf8_lossy(&output.stderr)
            );
            let path = dir.join(&output_dir).join("results.json");
            let doc: Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
            let r = &doc["reference_scf"];
            assert_eq!(r["diagnostics"]["level"], level);
            assert_eq!(r["stop_to_verified_clock"], "reference-scf-wall-1");
            let expected = reader::stop_to_verified_summary(&doc).unwrap();
            assert_eq!(expected["status"], "VERIFIED");
            if variant == "first" {
                let first_stop = expected["first_stop_map"].as_u64().unwrap() as usize - 1;
                for field in ["first_prediction_ha", "contraction_prediction_ha"] {
                    let mut invalid = doc.clone();
                    invalid["reference_scf"]["iterations"][first_stop]["verification_cost"]
                        [field] = json!(-1.);
                    assert!(reader::stop_to_verified_summary(&invalid).is_err());
                }
            }
            let costs = &r["verification_cost"];
            if variant == "off" || variant == "first_inert" {
                assert!(r.get("verification_cost").is_none());
                assert!(expected.get("skipped").is_none());
            } else {
                let maps = r["iterations"].as_array().unwrap();
                let reused = maps
                    .iter()
                    .filter(|m| !m["verification_cost"]["reused_from"].is_null())
                    .count();
                let skipped = maps
                    .iter()
                    .filter(|m| m["verification_cost"]["skipped"] == true)
                    .count();
                assert_eq!(costs["reused"], reused);
                assert_eq!(costs["skipped"], skipped);
                assert_eq!(costs["attempts"], expected["attempt_count"]);
                assert_eq!(costs["logical_scf_iterations"], maps.len());
                let ordinary = maps.len() - reused
                    + maps
                        .iter()
                        .filter(|m| m["first_iteration_redo"] == true)
                        .count();
                assert_eq!(costs["ordinary_map_requests"], ordinary);
                assert_eq!(
                    costs["map_requests"],
                    ordinary as u64 + expected["attempt_count"].as_u64().unwrap()
                );
                let calls = maps
                    .iter()
                    .map(|m| m["calls"].as_array().unwrap().len())
                    .sum::<usize>()
                    + r["run"]["verification_calls"].as_array().unwrap().len();
                assert_eq!(costs["solver_calls"], calls);
                assert_eq!(costs["evaluated_maps"], calls);
                assert_eq!(
                    costs["verification_maps"],
                    r["run"]["verification_calls"].as_array().unwrap().len()
                );
                assert_eq!(
                    costs["evaluated_ordinary_maps"],
                    maps.iter()
                        .map(|m| m["calls"].as_array().unwrap().len())
                        .sum::<usize>()
                );
                if gate == "on" {
                    assert_eq!(expected["skipped"], skipped);
                }
                if active && variant == "gate" {
                    assert!(skipped > 0);
                }
                if active && variant == "reuse" {
                    assert!(reused > 0);
                }
                eprintln!(
                    "CLI_VERIFY_COST active={active} level={level} variant={variant} counts={costs}"
                );
            }
            for (map, history) in r["iterations"]
                .as_array()
                .unwrap()
                .iter()
                .zip(doc["convergence"]["history"].as_array().unwrap())
            {
                assert_eq!(
                    map.get("verification_cost"),
                    history["reference_scf"].get("verification_cost")
                );
                if variant == "off" || variant == "first_inert" {
                    assert!(map.get("verification_cost").is_none());
                } else {
                    let cost = &map["verification_cost"];
                    assert_eq!(
                        cost.as_object().unwrap().len(),
                        if variant == "first" { 8 } else { 5 }
                    );
                    assert_eq!(cost["gate_enabled"], gate == "on");
                    let attempt = r["verification"]["attempts"]
                        .as_array()
                        .unwrap()
                        .iter()
                        .find(|a| a["iteration"] == map["iteration"]);
                    if let Some(attempt) = attempt {
                        assert_eq!(cost["actual_error_ha"], attempt["coulomb_error_ha"]);
                        if variant == "first" {
                            assert!(cost["first_prediction_ha"].as_f64().is_some());
                            let prediction = expected["verification_predictions"]
                                .as_array()
                                .unwrap()
                                .iter()
                                .find(|p| p["iteration"] == map["iteration"])
                                .unwrap();
                            assert_eq!(prediction["actual_error_ha"], cost["actual_error_ha"]);
                            assert_eq!(
                                prediction["first_prediction_ha"],
                                cost["first_prediction_ha"]
                            );
                            assert_eq!(
                                prediction["contraction_prediction_ha"],
                                cost["contraction_prediction_ha"]
                            );
                            assert_eq!(
                                prediction["predicted_error_ha"],
                                cost["predicted_error_ha"]
                            );
                        }
                    } else {
                        assert!(cost["actual_error_ha"].is_null());
                    }
                    if cost["skipped"] == true {
                        assert!(cost["predicted_error_ha"].as_f64().unwrap() >= 5e-9);
                        assert!(attempt.is_none());
                    }
                }
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
    }
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn stop_to_verified_skipped_attempt_partition() {
    let mut doc = record();
    // First stop is skipped: its work remains in the interval, not an attempt.
    let map = &mut doc["reference_scf"]["iterations"][9];
    map["verification_start_s"] = Value::Null;
    map["verification_end_s"] = Value::Null;
    map["verification_cost"] = json!({"gate_enabled":true,"skipped":true});
    doc["reference_scf"]["verification"]["attempts"]
        .as_array_mut()
        .unwrap()
        .remove(0);
    let summary = reader::stop_to_verified_summary(&doc).unwrap();
    assert_eq!(summary["status"], "VERIFIED");
    assert_eq!(summary["first_stop_map"], 10);
    assert_eq!(summary["skipped"], 1);
    assert_eq!(summary["attempt_count"], 2);
    assert_eq!(summary["interval_wall_s"], 23.);
    assert_eq!(summary["verification_attempts_s"], 4.);
    assert_eq!(summary["ordinary_scf_after_failure_s"], 6.);
    assert_eq!(summary["unassigned_s"], 13.);
    assert_eq!(summary["closure_error_s"], 0.);
    let mut bad = doc.clone();
    bad["reference_scf"]["iterations"][9]["verification_cost"]["gate_enabled"] = json!(false);
    assert!(reader::stop_to_verified_summary(&bad).is_err());
    let mut bad = doc.clone();
    bad["reference_scf"]["iterations"][9]["verification_start_s"] = json!(93.);
    bad["reference_scf"]["iterations"][9]["verification_end_s"] = json!(95.);
    assert!(reader::stop_to_verified_summary(&bad).is_err());
    doc["reference_scf"]["iterations"]
        .as_array_mut()
        .unwrap()
        .truncate(10);
    doc["reference_scf"]["verification"]["attempts"] = json!([]);
    doc["reference_scf"]["verification"]["passed"] = json!(false);
    assert_eq!(
        reader::stop_to_verified_summary(&doc).unwrap()["status"],
        "UNVERIFIED"
    );
}
