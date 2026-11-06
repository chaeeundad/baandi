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
    assert_eq!(run::options(&config).unwrap().gamma_real, GammaReal::Off);
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
