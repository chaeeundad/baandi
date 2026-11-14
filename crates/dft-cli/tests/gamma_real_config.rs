#[allow(dead_code)]
#[path = "../../dft-engine/tests/common/mod.rs"]
mod common;
#[allow(dead_code)]
#[path = "../../dft-engine/tests/common/checkpoint.rs"]
mod fixture;
use dft_cli::{checkpoint_io, config::Config, run};
use dft_engine::{
    CompletionMode, ExecutionConfig, GammaReal, Mixing, ScfSession, checkpoint::CheckpointKind,
};

#[test]
fn key_validation_and_engine_plumbing() {
    let text = include_str!("../../../examples/si_lda.toml");
    let config = Config::parse(text).unwrap();
    assert_eq!(run::options(&config).unwrap().gamma_real, GammaReal::Auto);
    for (key, value) in [
        ("off", GammaReal::Off),
        ("auto", GammaReal::Auto),
        ("on", GammaReal::On),
    ] {
        let config = Config::parse(&format!("{text}\n[scf]\ngamma_real = \"{key}\"\n")).unwrap();
        assert_eq!(run::options(&config).unwrap().gamma_real, value);
    }
    for value in ["\"ON\"", "\"invalid\"", "true", "1"] {
        assert!(Config::parse(&format!("{text}\n[scf]\ngamma_real = {value}\n")).is_err());
    }
}

#[test]
fn checkpoint_optional_option_and_roundtrip() {
    let (model, execution) = fixture::identity(0, 1);
    let mut s = ScfSession::new(
        fixture::fixture(0, Mixing::default(), 1, 200),
        CompletionMode::Strict,
        ExecutionConfig::default(),
    )
    .unwrap();
    s.initialize().unwrap();
    let mut cp = s
        .checkpoint(CheckpointKind::ExactContinuation, model, execution)
        .unwrap();
    // The default is now Auto; the legacy wire format is the Off header.
    cp.continuation.as_mut().unwrap().options.gamma_real = GammaReal::Off;
    let old = checkpoint_io::encode(&cp).unwrap();
    let h = u64::from_le_bytes(old[8..16].try_into().unwrap()) as usize;
    let header: serde_json::Value = serde_json::from_slice(&old[24..24 + h]).unwrap();
    // Independent wire check: legacy off headers have no new key or version.
    assert_eq!(header["schema_version"], 1);
    assert!(
        header["continuation"]["options"]
            .get("gamma_real")
            .is_none()
    );
    assert_eq!(
        checkpoint_io::decode(&old)
            .unwrap()
            .continuation
            .unwrap()
            .options
            .gamma_real,
        GammaReal::Off
    );
    for mode in [GammaReal::Auto, GammaReal::On] {
        cp.continuation.as_mut().unwrap().options.gamma_real = mode;
        let bytes = checkpoint_io::encode(&cp).unwrap();
        assert_eq!(
            checkpoint_io::decode(&bytes)
                .unwrap()
                .continuation
                .unwrap()
                .options
                .gamma_real,
            mode
        );
    }
}

#[allow(dead_code)]
#[path = "common/mod.rs"]
mod cli_fixture;

#[test]
fn standard_combination_reaches_real_engine_and_labels_output() {
    use std::fs;
    let dir = std::env::temp_dir().join(format!("p26g-cli-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("syn.xml"),
        cli_fixture::gaussian_with(2., 1., true, 401).replace("element=\"Xx\"", "element=\"He\""),
    )
    .unwrap();
    for mode in ["auto", "on"] {
        let text = format!(
            "[system]\nname='R6 synthetic'\nlattice_unit='bohr'\nlattice_vectors=[[4.0,0.0,0.0],[0.0,4.0,0.0],[0.0,0.0,4.0]]\n[[system.atoms]]\nspecies='He'\nfractional=[0.13,0.21,0.32]\n[pseudopotentials]\nHe='syn.xml'\n[basis]\necutwfc_ry=5.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[scf]\nnbands=5\nconvergence_profile='standard'\ngamma_real='{mode}'\n[output]\ndirectory='out-{mode}'\n"
        );
        let config = Config::parse(&text).unwrap();
        let opts = run::options(&config).unwrap();
        assert_eq!(opts.standard, Some(dft_engine::StandardProfile::default()));
        assert_ne!(opts.gamma_real, GammaReal::Off);
        fs::write(dir.join("input.toml"), &text).unwrap();
        let output = std::process::Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .current_dir(&dir)
            .args(["run", "input.toml"])
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let result: serde_json::Value = serde_json::from_slice(
            &fs::read(dir.join(format!("out-{mode}/results.json"))).unwrap(),
        )
        .unwrap();
        assert_eq!(result["convergence_profile"], "standard");
        assert_eq!(result["certified"], false);
        assert_eq!(result["gamma_real_used"], true);
        assert!(result["gamma_real_fallback_reason"].is_null());
        assert_eq!(result["standard_threshold_ha"], 5e-9);
        assert_eq!(result["work"]["certificate_solves"], 0);
        assert!(String::from_utf8_lossy(&output.stdout).contains("certified: false"));
    }
    fs::remove_dir_all(dir).unwrap();
}

/// P0.2 work timers (`BAANDI_DAVIDSON_STAGES=on`): results.json carries the FFT
/// split by kind, the h_apply parts and every block per run stage; the pooled
/// keys are the sums of the run stages; turning the timers on changes no
/// energy, eigenvalue, density bit or iteration count.
#[test]
fn work_timers_report_by_kind_and_run_stage() {
    use std::fs;
    let dir = std::env::temp_dir().join(format!("p02-cli-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("syn.xml"),
        cli_fixture::gaussian_with(2., 1., true, 401).replace("element=\"Xx\"", "element=\"He\""),
    )
    .unwrap();
    let text = "[system]\nname='P0.2 synthetic'\nlattice_unit='bohr'\nlattice_vectors=[[4.0,0.0,0.0],[0.0,4.0,0.0],[0.0,0.0,4.0]]\n[[system.atoms]]\nspecies='He'\nfractional=[0.13,0.21,0.32]\n[pseudopotentials]\nHe='syn.xml'\n[basis]\necutwfc_ry=5.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[scf]\nnbands=5\n";
    fs::write(dir.join("input.toml"), text).unwrap();
    let run = |stages: &str| {
        let out = format!("out-{stages}");
        let output = std::process::Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .current_dir(&dir)
            .env("BAANDI_DAVIDSON_STAGES", stages)
            .args(["run", "input.toml", "--output-dir", &out])
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        serde_json::from_slice::<serde_json::Value>(
            &fs::read(dir.join(out).join("results.json")).unwrap(),
        )
        .unwrap()
    };
    let off = run("off");
    let on = run("on");
    assert_eq!(on["gamma_real_used"], true);
    for key in ["total_energy", "eigenvalues"] {
        assert_eq!(on[key], off[key], "{key}");
    }
    assert_eq!(
        on["convergence"]["iterations"],
        off["convergence"]["iterations"]
    );
    let density = |out: &str| fs::read(dir.join(out).join("density.bin")).ok();
    assert!(density("out-on").is_some());
    assert_eq!(density("out-on"), density("out-off"));
    let work = &on["convergence"]["work"];
    let stages = &work["davidson_stages"];
    let zero_off = &off["convergence"]["work"]["davidson_stages"];
    assert_eq!(zero_off["fft"]["pair_transforms"], 0);
    assert_eq!(zero_off["h_apply_parts"]["apply_calls"], 0);
    // Existing keys keep their place and meaning.
    for key in [
        "timed_solves",
        "iteration_s_by_correction_width",
        "iterations_by_correction_width",
        "stage_s",
        "note",
    ] {
        assert!(stages.get(key).is_some(), "{key}");
    }
    let fft = &stages["fft"];
    assert_eq!(fft.as_object().unwrap().len(), 14);
    let parts = &stages["h_apply_parts"];
    assert_eq!(parts.as_object().unwrap().len(), 7);
    let n = |v: &serde_json::Value| v.as_u64().unwrap();
    let f = |v: &serde_json::Value| v.as_f64().unwrap();
    assert!(n(&fft["pair_transforms"]) > 0);
    assert!(n(&fft["density_transforms"]) > 0);
    assert!(n(&parts["apply_calls"]) > 0);
    // One pair carries at most two columns; every h_apply column is in a pair.
    assert!(n(&fft["pair_columns"]) <= 2 * n(&fft["pair_transforms"]));
    assert_eq!(fft["pair_columns"], parts["apply_columns"]);
    let by = stages["by_run_stage"].as_object().unwrap();
    let mut labels: Vec<_> = by.keys().map(String::as_str).collect();
    labels.sort_unstable();
    assert_eq!(labels, ["fine", "ladder_coarse", "verification"]);
    for key in [
        "pair_transforms",
        "pair_columns",
        "density_transforms",
        "other_transforms",
    ] {
        let sum: u64 = by.values().map(|b| n(&b["fft"][key])).sum();
        assert_eq!(sum, n(&fft[key]), "{key}");
    }
    let sum: u64 = by.values().map(|b| n(&b["timed_solves"])).sum();
    assert_eq!(sum, n(&stages["timed_solves"]));
    let h: f64 = by.values().map(|b| f(&b["stage_s"]["h_apply"])).sum();
    assert!((h - f(&stages["stage_s"]["h_apply"])).abs() <= 1e-9 * (1. + h));
    if n(&work["certificate_solves"]) > 0 {
        assert!(n(&by["verification"]["timed_solves"]) > 0);
        assert!(n(&by["verification"]["fft"]["pair_transforms"]) > 0);
    }
    fs::remove_dir_all(dir).unwrap();
}
