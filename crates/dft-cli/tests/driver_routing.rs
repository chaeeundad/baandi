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
    let z = if name == "l06-checkpoint" { 2. } else { 4. };
    let mut pp =
        common::gaussian_with(z, 1., true, 161).replace("element=\"Xx\"", "element=\"Si\"");
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
        pp = pp.replace("number_of_wfc=\"0\"", "number_of_wfc=\"1\"").replace("</UPF>", &format!("<PP_PSWFC><PP_CHI.1 l=\"{l}\" label=\"s\" occupation=\"{z}\">{values}</PP_CHI.1></PP_PSWFC></UPF>"));
    }
    pp
}

#[test]
fn corpus_precedence_and_fallback_matrix() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/routing");
    let expectations: Value =
        serde_json::from_str(&fs::read_to_string(root.join("expectations.json")).unwrap()).unwrap();
    for (name, _) in expectations.as_object().unwrap() {
        let text = fs::read_to_string(root.join(format!("{name}.toml"))).unwrap();
        for strict in [false, true] {
            let _scope =
                RunConfig::from_pairs(&[("reference_strict", if strict { "on" } else { "off" })])
                    .unwrap()
                    .enter();
            for explicit in [false, true] {
                let text = if explicit {
                    text.replace("[scf]", "[scf]\ndriver='reference'")
                } else {
                    text.clone()
                };
                if name.starts_with("l09-") {
                    let error = Config::parse(&text).unwrap_err().to_string();
                    assert!(error.contains("unknown"), "{name}: {error}");
                    continue;
                }
                let config = Config::parse(&text).unwrap_or_else(|e| panic!("{name}: {e}"));
                let pp = dft_upf::parse_upf(&atomic(name)).unwrap();
                for mode in [Mode::Classic, Mode::Observe, Mode::Agent] {
                    for command in [
                        Command::Run,
                        Command::Relax,
                        Command::Bands,
                        Command::Dos,
                        Command::Checkpoint,
                        Command::Restart,
                        Command::Explain,
                        Command::Session,
                    ] {
                        let selection = route::select(&config, mode, command, &[&pp]);
                        let g4d = name == "l01-free-moment-pbe-fd";
                        assert_eq!(selection.strict(), strict || g4d, "{name}");
                        assert_eq!(
                            selection.reason,
                            if strict {
                                Some(route::STRICT)
                            } else if g4d {
                                Some(route::G4D_FALLBACK)
                            } else {
                                None
                            }
                        );
                        let mut record = serde_json::json!({});
                        selection.record(&mut record);
                        assert_eq!(
                            record["scf"]["driver"],
                            if strict || g4d {
                                "strict_oracle"
                            } else {
                                "reference"
                            }
                        );
                        if g4d && !strict {
                            assert_eq!(
                                record["reference_driver_fallback_reason"],
                                route::G4D_FALLBACK
                            );
                        }
                    }
                }
            }
        }
    }
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
            (Some("on"), None, "reference", "default"),
            (None, Some("reference"), "reference", "input"),
            (Some("off"), Some("reference"), "reference", "input"),
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
            // This coarse synthetic fixture tests routing, not the XC charge gate.
            cmd.env("BAANDI_XC_SCREEN_NEGATIVE", "off");
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
    let _scope = RunConfig::resolve_with(|_| None, &[]).unwrap().enter();
    let config = Config::parse(text).unwrap();
    let pp = dft_upf::parse_upf(&atomic("gamma-lda")).unwrap();
    for command in [
        Command::Run,
        Command::Relax,
        Command::Explain,
        Command::Bands,
        Command::Dos,
    ] {
        let selection = route::select(&config, Mode::Classic, command, &[&pp]);
        assert_eq!(selection.driver, Driver::Reference);
        assert_eq!(selection.source, Source::Default);
    }
    for command in [Command::Checkpoint, Command::Restart] {
        let selection = route::select(&config, Mode::Classic, command, &[&pp]);
        assert_eq!(selection.driver, Driver::Reference);
        assert_eq!(selection.reason, None);
        assert!(dft_cli::run::options_for(&config, selection).is_ok());
    }
    let selection = route::select(
        &config,
        Mode::Classic,
        Command::RefscfMap { reference: None },
        &[&pp],
    );
    assert!(selection.reference());
    assert_eq!(selection.source, Source::Default);
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
    fs::write(dir.join("synthetic-atomic.upf"), atomic("l06-checkpoint")).unwrap();
    let input = dir.join("input.toml");
    // L6 now executes the reference driver; use a numerical fixture, rather
    // than the parse-only routing corpus's intentionally minimal Si2 grid.
    let text = r#"[system]
name = "L6 checkpoint routing"
lattice_unit = "bohr"
lattice_vectors = [[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0]]
[[system.atoms]]
species = "Si"
fractional = [0.13,0.21,0.32]
[pseudopotentials]
Si = "synthetic-atomic.upf"
[basis]
ecutwfc_ry = 5.0
kpoints = [[0.0,0.0,0.0]]
[xc]
functional = "lda-pz"
[scf]
nbands = 4
max_iterations = 2
"#;
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
                if strict { "strict_oracle" } else { "reference" }
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
                result["scf"]["driver_source"],
                if strict { "env" } else { "default" }
            );
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
fn moment_fd_reference_and_g4d_fallback_mixer_defaults() {
    // Fallback compatibility is an implicit/environment request contract.
    let _scope = RunConfig::resolve_with(
        |key| (key == "BAANDI_REFERENCE_SCF").then(|| "on".into()),
        &[],
    )
    .unwrap()
    .enter();
    let moment = include_str!("fixtures/routing/l02-moment-fd.toml");
    let free_lda = moment.replace("total_moment = 2.0\n", "");
    let fixed_pbe = moment.replace("lda-pw", "pbe");
    // R6a uses reference beta for moment FD; only L1/G4D retains fallback beta.
    for (text, beta) in [
        (
            include_str!("fixtures/routing/l01-free-moment-pbe-fd.toml"),
            0.3,
        ),
        (moment, 0.7),
        (free_lda.as_str(), 0.7),
        (fixed_pbe.as_str(), 0.7),
    ] {
        let config = Config::parse(text).unwrap();
        let selection = route::preliminary(&config);
        assert_eq!(selection.reference(), beta == 0.7);
        let options = dft_cli::run::options(&config).unwrap();
        assert!(matches!(
            options.mixing,
            dft_engine::Mixing::Pulay { beta: actual, history: 8 } if actual == beta
        ));
    }
}

#[test]
fn explicit_reference_legacy_keys_reject_before_mixer_validation() {
    let root =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/routing/explicit-reference");
    for entry in fs::read_dir(root).unwrap() {
        let text = fs::read_to_string(entry.unwrap().path()).unwrap();
        for strict in ["off", "on"] {
            let _scope = RunConfig::from_pairs(&[("reference_strict", strict)])
                .unwrap()
                .enter();
            assert!(
                Config::parse(&text)
                    .unwrap_err()
                    .to_string()
                    .contains("unknown")
            );
        }
    }
    let base = include_str!("fixtures/routing/gamma-lda.toml");
    assert!(
        Config::parse(&base.replace("[scf]", "[scf]\ndriver='legacy'"))
            .unwrap_err()
            .to_string()
            .contains("unknown")
    );
}

#[test]
fn explain_prints_selection_before_validation_errors() {
    let dir = std::env::temp_dir().join(format!("r1b-explain-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let base = include_str!("fixtures/routing/explicit-reference/l09-linear.toml");
    for (name, text, fragment) in [
        ("conflict", base.to_string(), "scf.mixing"),
        (
            "multiple",
            base.replace(
                "mixing = \"linear\"",
                "mixing = \"broyden\"\nenergy_tolerance_ha = 1e-8\neigen_tolerance = 1e-9",
            ),
            "unknown field",
        ),
        (
            "preflight",
            base.replace("mixing = \"linear\"", "mixing = \"broyden\""),
            "synthetic-atomic.upf",
        ),
    ] {
        let input = dir.join(format!("{name}.toml"));
        let log = dir.join(format!("{name}.log"));
        fs::write(&input, text).unwrap();
        let output = fs::File::create(&log).unwrap();
        let mut command = std::process::Command::new(env!("CARGO_BIN_EXE_dft-cli"));
        for (key, _) in std::env::vars().filter(|(k, _)| k.starts_with("BAANDI_")) {
            command.env_remove(key);
        }
        let status = command
            .arg("explain")
            .arg(&input)
            .stdout(output.try_clone().unwrap())
            .stderr(output)
            .status()
            .unwrap();
        assert!(!status.success());
        let text = fs::read_to_string(log).unwrap();

        assert!(text.contains(fragment), "{text}");
        assert!(!text.contains("unknown mixer"), "{text}");
    }
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn optional_atomic_data_runs_on_reference() {
    let dir = std::env::temp_dir().join(format!("r6c-start-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    for variant in ["no-rhoatom", "no-atomic", "zero-rhoatom", "no-pswfc", "l4"] {
        let upf = dir.join("syn.xml");
        let generated = std::process::Command::new("python3")
            .arg("-I")
            .arg(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/../../scripts/gen_synthetic_upf.py"
            ))
            .arg(&upf)
            .args([
                "--valence",
                "2",
                "--radius",
                "1.8",
                "--points",
                "161",
                "--start-variant",
                variant,
            ])
            .output()
            .unwrap();
        assert!(
            generated.status.success(),
            "{}",
            String::from_utf8_lossy(&generated.stderr)
        );
        // Only identity metadata changes; this remains a non-physical fixture.
        let text = fs::read_to_string(&upf)
            .unwrap()
            .replace("element=\"Xx\"", "element=\"He\"");
        let mut cases = vec![(false, "gamma"), (true, "gamma")];
        if variant == "no-atomic" {
            cases.extend([(true, "k"), (true, "fd"), (true, "spin")]);
        }
        if variant == "no-pswfc" {
            cases.push((true, "spin"));
        }
        for (explicit, mode) in cases {
            fs::write(
                &upf,
                if mode == "spin" {
                    text.replace("SLA PZ", "SLA PW")
                } else {
                    text.clone()
                },
            )
            .unwrap();
            let input = dir.join("input.toml");
            let mut input_text = format!(
                "[system]\nname='R6C synthetic'\nlattice_unit='bohr'\nlattice_vectors=[[7.0,0.0,0.0],[0.2,7.0,0.0],[0.0,0.1,7.0]]\n[[system.atoms]]\nspecies='Xx'\nfractional=[0.13,0.21,0.32]\n[pseudopotentials]\nXx={{path='syn.xml',element='He'}}\n[basis]\necutwfc_ry=2.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[scf]\nnbands=3\nmax_iterations=160\n{}",
                if explicit { "driver='reference'\n" } else { "" }
            );
            if mode != "gamma" {
                input_text = input_text.replace("nbands=3", "nbands=5");
            }
            match mode {
                "k" => {
                    input_text = input_text
                        .replace("kpoints=[[0.0,0.0,0.0]]", "kpoints=[[0.13,-0.21,0.07,1.0]]")
                }
                "fd" => input_text
                    .push_str("[occupations]\nmode='fermi_dirac'\nsmearing_ha=0.01\nband_cap=5\n"),
                "spin" => {
                    input_text = input_text.replace("lda-pz", "lda-pw");
                    input_text.push_str("[spin]\nmode='collinear'\nn_up=1\nn_down=1\n");
                }
                _ => (),
            }
            fs::write(&input, input_text).unwrap();
            let out = dir.join(format!("{variant}-{explicit}-{mode}"));
            let mut cmd = std::process::Command::new(env!("CARGO_BIN_EXE_dft-cli"));
            for (key, _) in std::env::vars().filter(|(k, _)| k.starts_with("BAANDI_")) {
                cmd.env_remove(key);
            }
            let output = cmd
                .arg("run")
                .arg(&input)
                .arg("--output-dir")
                .arg(&out)
                .output()
                .unwrap();
            assert!(
                output.status.success(),
                "{variant}/{explicit}/{mode}: {}\n{}",
                String::from_utf8_lossy(&output.stdout),
                String::from_utf8_lossy(&output.stderr)
            );
            let result: Value =
                serde_json::from_slice(&fs::read(out.join("results.json")).unwrap()).unwrap();
            assert_eq!(result["scf"]["driver"], "reference");
            assert_eq!(
                result["reference_scf"]["start"]["density"],
                if mode == "spin" || matches!(variant, "no-rhoatom" | "no-atomic" | "zero-rhoatom")
                {
                    "uniform"
                } else {
                    "atomic"
                }
            );
            assert_eq!(
                result["reference_scf"]["start"]["orbitals"],
                if matches!(variant, "no-rhoatom" | "zero-rhoatom") {
                    "atomic"
                } else {
                    "kinetic"
                }
            );

            assert_eq!(result["scf"]["driver_reason"], Value::Null);
            assert_eq!(
                result["scf"]["driver_source"],
                if explicit { "input" } else { "default" }
            );
        }
    }
    fs::remove_dir_all(dir).unwrap();
}
