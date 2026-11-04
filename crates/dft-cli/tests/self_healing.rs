#[allow(dead_code)]
#[path = "../../dft-engine/tests/common/mod.rs"]
mod common;
#[path = "../src/heal.rs"]
mod heal;
use dft_engine::{ScfOptions, ScfProblem, SpeciesPotential};
use serde_json::json;
use std::{fs, path::PathBuf};

fn problem(budget: usize) -> ScfProblem {
    let p = common::potential(4., 1., true, 801);
    ScfProblem::new(
        &common::cell(),
        &common::atoms(),
        &[SpeciesPotential {
            label: "Xx",
            potential: &p,
        }],
        ScfOptions {
            max_iterations: budget,
            threads: 1,
            seed: 2026,
            ..Default::default()
        },
    )
    .unwrap()
}

#[test]
fn recorded_schema_rejects_unknown_authority_and_missing_version() {
    let text = include_str!("../../../tests/reference/runtime/agent-decision.json");
    let mut v: serde_json::Value = serde_json::from_str(text).unwrap();
    assert!(serde_json::from_value::<heal::Decision>(v.clone()).is_ok());
    v["authority"] = json!("operator");
    assert!(serde_json::from_value::<heal::Decision>(v.clone()).is_err());
    v.as_object_mut().unwrap().remove("authority");
    v.as_object_mut().unwrap().remove("observed_state_version");
    assert!(serde_json::from_value::<heal::Decision>(v).is_err());
}

#[test]
fn recorded_stale_or_future_observations_cannot_mutate_policy() {
    use dft_engine::{CompletionMode, ExecutionConfig, ExecutionMode, ScfSession, runtime::*};
    use std::time::{Duration, SystemTime};
    for (version, expected) in [
        (0, ControlError::StaleState),
        (99, ControlError::FutureState),
    ] {
        let mut decision: heal::Decision = serde_json::from_str(include_str!(
            "../../../tests/reference/runtime/agent-decision.json"
        ))
        .unwrap();
        decision.observed_state_version = version;
        let mut runtime = Runtime::new(
            ScfSession::new(
                problem(1),
                CompletionMode::RetainLast,
                ExecutionConfig {
                    mode: ExecutionMode::AgentControlled,
                    calc_id: "self-healing-demo".into(),
                    ..Default::default()
                },
            )
            .unwrap(),
            RuntimeConfig::default(),
        )
        .unwrap();
        runtime.initialize().unwrap();
        let before = format!("{:?}", runtime.snapshot());
        runtime
            .submit(
                decision.request(),
                "recorded-agent".into(),
                Authority::Approved(Approval {
                    id: "host-review".into(),
                    approver: "host".into(),
                    expires_at: SystemTime::now() + Duration::from_secs(60),
                    allow_revalidation: false,
                }),
            )
            .unwrap();
        runtime.safe_point().unwrap();
        assert_eq!(
            runtime.audit().last().unwrap().outcome,
            AuditOutcome::Rejected(expected)
        );
        assert_eq!(before, format!("{:?}", runtime.snapshot()));
    }
}

#[test]
#[ignore = "real SCF demo; scripts/demo_self_healing.sh (release, offline)"]
fn three_arm_demo() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let out = std::env::var_os("DFT_HEAL_OUTPUT")
        .map(PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join(format!("dft-heal-{}", std::process::id())));
    fs::create_dir_all(&out).unwrap();
    // A failed rerun must not leave an old successful aggregate report.
    if out.join("results.json").exists() {
        fs::remove_file(out.join("results.json")).unwrap();
    }
    let decision = heal::read_decision(&root.join("tests/reference/runtime/agent-decision.json"));
    let rule = heal::Decision {
        schema_version: 1,
        calc_id: "self-healing-demo".into(),
        request_id: "demo-budget-1".into(),
        observed_state_version: 1,
        observed_policy_revision: 0,
        apply_after_iteration: 0,
        max_iterations: 200,
    };
    let mut summaries = vec![];
    let mut initial_density = None;
    let mut identities = std::collections::BTreeMap::new();
    let mut energies = std::collections::BTreeMap::new();
    for repeat in 0..2 {
        for (arm, budget, action, probe) in [
            ("reference", 1000, None, false),
            ("baseline", 1, None, false),
            ("rules", 1, Some(&rule), false),
            ("agent", 1, Some(&decision), true),
        ] {
            let name = format!("{arm}-{repeat}");
            let mut run = heal::run(problem(budget), &out.join(&name), action, arm, probe);
            assert_eq!(run.outcome.converged, arm != "baseline");
            let initial = run.summary["initial_density_sha256"]
                .as_str()
                .unwrap()
                .to_owned();
            if let Some(prior) = &initial_density {
                assert_eq!(&initial, prior);
            }
            initial_density = Some(initial);
            let audit = fs::read_to_string(out.join(&name).join("audit.jsonl")).unwrap();
            let mut previous = serde_json::Value::Null;
            let mut outcomes = vec![];
            for line in audit.lines() {
                let mut record: serde_json::Value = serde_json::from_str(line).unwrap();
                assert_eq!(record["previous_record_sha256"], previous);
                let digest = record
                    .as_object_mut()
                    .unwrap()
                    .remove("record_sha256")
                    .unwrap();
                assert_eq!(
                    digest,
                    json!(heal::hash(serde_json::to_vec(&record).unwrap()))
                );
                previous = digest;
                outcomes.push(record["outcome"].as_str().unwrap().to_owned());
            }
            let expected: Vec<&str> = match arm {
                "agent" => vec!["Rejected(Policy(InvalidBeta))", "Accepted", "Applied"],
                "rules" => vec!["Accepted", "Applied"],
                _ => vec![],
            };
            assert_eq!(outcomes, expected);
            let energy = run.outcome.result.total_hartree();
            let reference = *energies.entry("reference").or_insert(energy);
            let delta = (energy - reference).abs();
            if arm != "baseline" {
                assert!(delta <= 1e-10, "{name}: energy mismatch {delta} Ha");
                assert!(
                    run.outcome
                        .result
                        .history
                        .last()
                        .unwrap()
                        .fixed_point_certificate
                        .is_some()
                );
            }
            // Hash EVERY array element's IEEE bits, including all complex coefficients.
            let mut bytes = Vec::new();
            for v in &run.outcome.result.density {
                bytes.extend(v.to_bits().to_le_bytes());
            }
            for k in &run.outcome.result.k_points {
                for v in &k.solution.eigenvalues {
                    bytes.extend(v.to_bits().to_le_bytes());
                }
                for col in 0..k.solution.eigenvectors.cols() {
                    for v in k.solution.eigenvectors.column(col) {
                        bytes.extend(v.re.to_bits().to_le_bytes());
                        bytes.extend(v.im.to_bits().to_le_bytes());
                    }
                }
            }
            let arrays = heal::hash(&bytes);
            let normalized = heal::normalized_result(&mut run.outcome);
            let identity = (arrays.clone(), heal::hash(&normalized));
            if let Some(prior) = identities.get(arm) {
                assert_eq!(&identity, prior, "repeat {name}");
            }
            if arm == "agent" || arm == "rules" {
                assert_eq!(&identity, identities.get("reference").unwrap(), "{name}");
            }
            identities.insert(arm, identity);
            // Full normalized history/certificate with array hash; Debug text is build-local, not a wire format.
            fs::write(out.join(&name).join("normalized-result.txt"), normalized).unwrap();
            heal::write_json(
                &out.join(&name).join("arrays.json"),
                &json!({"sha256":arrays,"bytes":bytes.len(),
                "encoding":"little-endian f64 bits: density iterator, then each k eigenvalues and column-major complex (re,im)",
                "density_shape":run.outcome.result.density.shape()}),
            );
            run.summary["arm"] = json!(arm);
            run.summary["repeat"] = json!(repeat);
            run.summary["energy_error_hartree"] = if arm == "baseline" {
                json!(null)
            } else {
                json!(delta)
            };
            println!(
                "{name}: converged={} iterations={} energy={energy:.17} Ha delta={delta:.3e} Ha",
                run.outcome.converged,
                run.outcome.result.history.len()
            );
            summaries.push(run.summary);
        }
    }
    let sources: Vec<_> = [
        "Cargo.lock",
        "rust-toolchain.toml",
        "crates/dft-cli/src/heal.rs",
        "crates/dft-cli/tests/self_healing.rs",
        "crates/dft-engine/tests/common/mod.rs",
        "crates/dft-engine/src/scf.rs",
        "crates/dft-engine/src/session.rs",
        "crates/dft-engine/src/runtime.rs",
        "crates/dft-engine/src/policy.rs",
        "crates/dft-engine/src/events.rs",
        "scripts/demo_self_healing.sh",
        "crates/dft-engine/src/diagnostics.rs",
        "tests/reference/runtime/agent-decision.json",
    ]
    .iter()
    .map(|p| json!({"path":p,"sha256":heal::hash(fs::read(root.join(p)).unwrap())}))
    .collect();
    heal::write_json(
        &out.join("results.json"),
        &json!({"schema_version":1,"task_id":"P11E-SELF-HEALING-DEMO",
        "scope":"synthetic two-ion Gamma budget recovery, not instability healing or QE comparison",
        "command":"scripts/demo_self_healing.sh","reference":"fresh same-build BAANDI, max_iterations=1000",
        "energy_tolerance_hartree":1e-10,"repeated_arms":4,"repeats":2,"full_result_and_array_identity":true,
        "sources":sources,"runs":summaries}),
    );
}
