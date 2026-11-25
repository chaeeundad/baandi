//! R1 routing matrix: real parsed synthetic UPFs, no numerical SCF required.
#[allow(dead_code)]
mod common;
use dft_cli::{
    config::Config,
    driver_select::{self as route, Command, Driver, Mode, Source},
};
use dft_numerics::switches::RunConfig;
use serde_json::Value;
use std::{fs, path::Path};

fn atomic(name: &str) -> String {
    let mut pp =
        common::gaussian_with(4., 1., true, 161).replace("element=\"Xx\"", "element=\"Si\"");
    let p = dft_upf::parse_upf(&pp).unwrap();
    let values = p
        .mesh()
        .r()
        .iter()
        .map(|r| format!("{:.17e}", r * (-r * r).exp()))
        .collect::<Vec<_>>()
        .join(" ");
    if name != "l11-zero-rhoatom" {
        let start = pp.find("<PP_RHOATOM>").unwrap() + "<PP_RHOATOM>".len();
        let end = pp.find("</PP_RHOATOM>").unwrap();
        pp.replace_range(start..end, &values);
    }
    if name == "l11-missing-rhoatom" {
        let start = pp.find("<PP_RHOATOM>").unwrap();
        let end = pp.find("</PP_RHOATOM>").unwrap() + "</PP_RHOATOM>".len();
        pp.replace_range(start..end, "");
    }
    if name != "l11-missing-pswfc" {
        let l = if name == "l11-high-l" { 4 } else { 0 };
        pp = pp.replace("number_of_wfc=\"0\"", "number_of_wfc=\"1\"").replace("</UPF>", &format!("<PP_PSWFC><PP_CHI.1 l=\"{l}\" label=\"s\" occupation=\"4\">{values}</PP_CHI.1></PP_PSWFC></UPF>"));
    }
    pp
}
fn reason(key: &str) -> &'static str {
    match key {
        "G4D_FALLBACK" => route::G4D_FALLBACK,
        "G4B_FALLBACK" => route::G4B_FALLBACK,
        "NSCF" => route::NSCF,
        "VARIABLE_CELL" => route::VARIABLE_CELL,
        "EXTRAPOLATION" => route::EXTRAPOLATION,
        "CHECKPOINT" => route::CHECKPOINT,
        "RUNTIME" => route::RUNTIME,
        "SESSION" => route::SESSION,
        "LEGACY_INPUT" => route::LEGACY_INPUT,
        "COMPAT" => route::COMPAT,
        "ATOMIC_DATA" => route::ATOMIC_DATA,
        _ => panic!("unknown reason {key}"),
    }
}
#[test]
fn corpus_precedence_and_fallback_matrix() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/routing");
    let expectations: Value =
        serde_json::from_str(&fs::read_to_string(root.join("expectations.json")).unwrap()).unwrap();
    assert_eq!(
        expectations.as_object().unwrap().len(),
        fs::read_dir(&root)
            .unwrap()
            .filter(|e| e
                .as_ref()
                .unwrap()
                .path()
                .extension()
                .is_some_and(|ext| ext == "toml"))
            .count()
    );
    let mut count = 0;
    let mut parse_rejections = 0;
    for (name, expected) in expectations.as_object().unwrap() {
        let text = fs::read_to_string(root.join(format!("{name}.toml"))).unwrap();
        let mode = match expected["mode"].as_str().unwrap() {
            "observe" => Mode::Observe,
            "agent" => Mode::Agent,
            _ => Mode::Classic,
        };
        let command = match expected["command"].as_str().unwrap() {
            "relax" => Command::Relax,
            "bands" => Command::Bands,
            "dos" => Command::Dos,
            "checkpoint" => Command::Checkpoint,
            "session" => Command::Session,
            _ => Command::Run,
        };
        for strict in [false, true] {
            for (env, input, source) in [
                (None, None, Source::Default),
                (Some("on"), None, Source::Env),
                (None, Some("reference"), Source::Input),
                (Some("on"), Some("legacy"), Source::Input),
            ] {
                let _scope = RunConfig::resolve_with(
                    |key| match key {
                        "BAANDI_REFERENCE_SCF" => env.map(Into::into),
                        "BAANDI_REFERENCE_STRICT" if strict => Some("on".into()),
                        _ => None,
                    },
                    &[],
                )
                .unwrap()
                .enter();
                let input_text = if let Some(input) = input {
                    text.replace("[scf]", &format!("[scf]\ndriver = '{input}'"))
                } else {
                    text.clone()
                };
                let config = Config::parse_unvalidated(&input_text)
                    .unwrap_or_else(|e| panic!("{name}: {e}"));
                let parsed = dft_upf::parse_upf(&atomic(name));
                if let Some(fragment) = expected["upf_error"].as_str() {
                    assert!(parsed.unwrap_err().to_string().contains(fragment), "{name}");
                    parse_rejections += 1;
                    continue;
                }
                let pp = parsed.unwrap_or_else(|e| panic!("{name}: {e}"));
                let selection = route::select(&config, mode, command, &[&pp]);
                let wants_reference =
                    input == Some("reference") || (input.is_none() && env != Some("off"));
                let fallback = expected["reference_reason"].as_str().map(reason);
                let driver = if strict {
                    Driver::StrictOracle
                } else if wants_reference && fallback.is_none() {
                    Driver::Reference
                } else {
                    Driver::LegacyStandard
                };
                assert_eq!(
                    selection.driver, driver,
                    "{name}, {strict}, {env:?}, {input:?}"
                );
                assert_eq!(
                    selection.reason,
                    if strict {
                        Some(route::STRICT)
                    } else if wants_reference {
                        fallback
                    } else {
                        None
                    },
                    "{name}"
                );
                assert_eq!(
                    selection.source,
                    if strict { Source::Env } else { source },
                    "{name}"
                );
                // Invalid input stays invalid; choosing a fallback never repairs it.
                if !strict && let Some(fragment) = expected["parse_error"].as_str() {
                    assert!(
                        config
                            .validate()
                            .unwrap_err()
                            .to_string()
                            .contains(fragment)
                    );
                } else {
                    config.validate().unwrap_or_else(|e| panic!("{name}: {e}"));
                }
                if let Some(fragment) = expected["options_error"].as_str() {
                    assert!(
                        dft_cli::run::options_for(&config, selection)
                            .unwrap_err()
                            .to_string()
                            .contains(fragment),
                        "{name}"
                    );
                }
                let mut record = serde_json::json!({"scf":{"preserved":true}});
                selection.record(&mut record);
                assert_eq!(
                    record["scf"]["driver"],
                    serde_json::to_value(driver).unwrap()
                );
                assert_eq!(
                    record["scf"]["driver_reason"],
                    serde_json::json!(selection.reason)
                );
                assert_eq!(
                    record["scf"]["driver_source"],
                    serde_json::to_value(selection.source).unwrap()
                );
                assert_eq!(record["scf"]["preserved"], true);
                count += 1;
            }
        }
    }
    println!("{count} routing selections and {parse_rejections} UPF parse rejections passed");
}
#[test]
fn switch_reads_are_owned_by_selector() {
    fn visit(path: &Path) {
        for entry in fs::read_dir(path).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                visit(&path);
            } else if path.extension().is_some_and(|e| e == "rs")
                && path.file_name().unwrap() != "driver_select.rs"
            {
                let text = fs::read_to_string(&path).unwrap();
                let compact: String = text
                    .lines()
                    .filter(|l| !l.trim_start().starts_with("//"))
                    .collect::<String>()
                    .chars()
                    .filter(|c| !c.is_whitespace())
                    .collect();
                for id in ["REFERENCE_SCF", "REFERENCE_STRICT"] {
                    assert!(
                        !compact.contains(&format!("Id::{id}")),
                        "direct registry access in {}",
                        path.display()
                    );
                    assert!(
                        !compact.contains(&format!("var(\"BAANDI_{id}\""))
                            && !compact.contains(&format!("var_os(\"BAANDI_{id}\"")),
                        "direct environment read in {}",
                        path.display()
                    );
                }
            }
        }
    }
    visit(&Path::new(env!("CARGO_MANIFEST_DIR")).join("src"));
}

#[test]
fn cli_records_executed_driver_and_strict_precedence() {
    let dir = std::env::temp_dir().join(format!("r1-routing-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(dir.join("syn.xml"), atomic("gamma-lda")).unwrap();
    let base = "[system]\nname='routing'\nlattice_unit='bohr'\nlattice_vectors=[[7.0,0.0,0.0],[0.2,7.0,0.0],[0.0,0.1,7.0]]\n[[system.atoms]]\nspecies='Si'\nfractional=[0.13,0.21,0.32]\n[pseudopotentials]\nSi='syn.xml'\n[basis]\necutwfc_ry=2.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[scf]\nnbands=3\nmax_iterations=2\n";
    for strict in [false, true] {
        for (index, (env, input, driver, source)) in [
            (None, None, "reference", "default"),
            (Some("on"), None, "reference", "env"),
            (None, Some("reference"), "reference", "input"),
            (Some("on"), Some("legacy"), "legacy", "input"),
        ]
        .into_iter()
        .enumerate()
        {
            let text = if let Some(input) = input {
                format!("{base}driver='{input}'\n")
            } else {
                base.into()
            };
            let input_path = dir.join("input.toml");
            fs::write(&input_path, text).unwrap();
            let out = dir.join(format!("out-{strict}-{index}"));
            let mut cmd = std::process::Command::new(env!("CARGO_BIN_EXE_dft-cli"));
            for (key, _) in std::env::vars().filter(|(k, _)| k.starts_with("BAANDI_")) {
                cmd.env_remove(key);
            }
            cmd.args([
                "run",
                input_path.to_str().unwrap(),
                "--output-dir",
                out.to_str().unwrap(),
            ]);
            if let Some(env) = env {
                cmd.env("BAANDI_REFERENCE_SCF", env);
            }
            if strict {
                cmd.env("BAANDI_REFERENCE_STRICT", "on");
            }
            let output = cmd.output().unwrap();
            assert!(
                matches!(output.status.code(), Some(0 | 3)),
                "{}\n{}",
                String::from_utf8_lossy(&output.stdout),
                String::from_utf8_lossy(&output.stderr)
            );
            let result: Value =
                serde_json::from_slice(&fs::read(out.join("results.json")).unwrap()).unwrap();
            assert_eq!(
                result["scf"]["driver"],
                if strict { "strict_oracle" } else { driver }
            );
            assert_eq!(
                result["scf"]["driver_source"],
                if strict { "env" } else { source }
            );
            assert_eq!(
                result["scf"]["driver_reason"],
                if strict {
                    serde_json::json!(route::STRICT)
                } else {
                    Value::Null
                }
            );
            assert_eq!(
                String::from_utf8_lossy(&output.stderr)
                    .matches("SCF driver:")
                    .count(),
                1
            );
            if strict {
                assert!(result.get("reference_scf").is_none());
            }
        }
    }
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn command_boundaries_and_input_switch_provenance() {
    let text = include_str!("fixtures/routing/gamma-lda.toml");
    let _scope = RunConfig::from_pairs(&[("reference_scf", "on")])
        .unwrap()
        .enter();
    let config = Config::parse(text).unwrap();
    let pp = dft_upf::parse_upf(&atomic("gamma-lda")).unwrap();
    for command in [Command::Run, Command::Relax, Command::Explain] {
        let selection = route::select(&config, Mode::Classic, command, &[&pp]);
        assert_eq!(selection.driver, Driver::Reference);
        assert_eq!(selection.source, Source::Input);
    }
    for command in [Command::Checkpoint, Command::Restart] {
        let selection = route::select(&config, Mode::Classic, command, &[&pp]);
        assert_eq!(selection.driver, Driver::LegacyStandard);
        assert_eq!(selection.reason, Some(route::CHECKPOINT));
        assert!(
            dft_cli::run::options_for(&config, selection)
                .unwrap()
                .standard
                .is_some()
        );
    }
    let selection = route::select(
        &config,
        Mode::Classic,
        Command::RefscfMap { reference: None },
        &[&pp],
    );
    assert!(selection.reference());
    assert_eq!(selection.source, Source::Input);
    for reference in [false, true] {
        let selection = route::select(
            &config,
            Mode::Classic,
            Command::RefscfMap {
                reference: Some(reference),
            },
            &[&pp],
        );
        assert_eq!(selection.reference(), reference);
        let _strict = RunConfig::from_pairs(&[("reference_strict", "on")])
            .unwrap()
            .enter();
        assert!(
            route::select(
                &config,
                Mode::Classic,
                Command::RefscfMap {
                    reference: Some(reference)
                },
                &[&pp]
            )
            .strict()
        );
    }
}

#[test]
fn checkpoint_restart_and_explain_record_selection() {
    let dir = std::env::temp_dir().join(format!("r1-checkpoint-routing-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(dir.join("synthetic-atomic.upf"), atomic("gamma-lda")).unwrap();
    let input = dir.join("input.toml");
    let text = include_str!("fixtures/routing/gamma-lda.toml")
        .replace("max_iterations = 100", "max_iterations = 2");
    fs::write(&input, text).unwrap();
    let invoke = |args: &[&str], strict: bool| {
        let mut command = std::process::Command::new(env!("CARGO_BIN_EXE_dft-cli"));
        for (key, _) in std::env::vars().filter(|(k, _)| k.starts_with("BAANDI_")) {
            command.env_remove(key);
        }
        command
            .current_dir(&dir)
            .args(args)
            .env("BAANDI_REFERENCE_SCF", "on")
            .env("BAANDI_REFERENCE_STRICT", if strict { "on" } else { "off" })
            .output()
            .unwrap()
    };
    for strict in [false, true] {
        let cp = format!("state-{strict}.cp");
        let end = format!("end-{strict}.cp");
        let output = invoke(
            &["checkpoint", "input.toml", &cp, "--iteration", "1"],
            strict,
        );
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let output = invoke(&["restart", "input.toml", &cp, &end], strict);
        // An exhausted two-step SCF still exports its result and route.
        assert!(
            output.status.success()
                || String::from_utf8_lossy(&output.stderr).contains("NOT CONVERGED"),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        for name in [cp, end] {
            let result: Value = serde_json::from_slice(
                &fs::read(dir.join(format!("{name}.results.json"))).unwrap(),
            )
            .unwrap();
            assert_eq!(
                result["scf"]["driver"],
                if strict { "strict_oracle" } else { "legacy" }
            );
            assert_eq!(
                result["scf"]["driver_reason"],
                if strict {
                    route::STRICT
                } else {
                    route::CHECKPOINT
                }
            );
            assert_eq!(result["scf"]["driver_source"], "env");
        }
        let output = invoke(&["explain", "input.toml"], strict);
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let text = String::from_utf8(output.stdout).unwrap();
        let result: Value =
            serde_json::from_str(text.lines().find(|l| l.starts_with("{\"scf\":")).unwrap())
                .unwrap();
        assert_eq!(
            result["scf"]["driver"],
            if strict { "strict_oracle" } else { "reference" }
        );
    }
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn preexisting_spin_fd_fallback_mixer_defaults_are_preserved() {
    let _scope = RunConfig::from_pairs(&[("reference_scf", "on")])
        .unwrap()
        .enter();
    for text in [
        include_str!("fixtures/routing/l01-free-moment-pbe-fd.toml"),
        include_str!("fixtures/routing/l02-moment-fd.toml"),
    ] {
        let config = Config::parse(text).unwrap();
        let options = dft_cli::run::options(&config).unwrap();
        assert!(matches!(
            options.mixing,
            dft_engine::Mixing::Pulay {
                beta: 0.5,
                history: 8
            }
        ));
        let config = Config::parse(&text.replace("[scf]", "[scf]\ndriver='legacy'\n")).unwrap();
        let options = dft_cli::run::options(&config).unwrap();
        assert!(matches!(
            options.mixing,
            dft_engine::Mixing::Pulay {
                beta: 0.5,
                history: 6
            }
        ));
    }
}
