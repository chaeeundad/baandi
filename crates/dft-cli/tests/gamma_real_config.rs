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
