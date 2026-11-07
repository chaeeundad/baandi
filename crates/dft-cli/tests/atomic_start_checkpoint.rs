#[allow(dead_code)]
#[path = "../../dft-engine/tests/common/mod.rs"]
mod common;
#[allow(dead_code)]
#[path = "../../dft-engine/tests/common/checkpoint.rs"]
mod fixture;
use dft_cli::checkpoint_io as io;
use dft_engine::{checkpoint::*, *};
use sha2::{Digest, Sha256};

#[test]
fn old_checkpoint_defaults_and_atomic_wire_roundtrip() {
    let (m, e) = fixture::identity(0, 1);
    let mut s = ScfSession::new(
        fixture::fixture(0, Mixing::default(), 1, 200),
        CompletionMode::Strict,
        ExecutionConfig::default(),
    )
    .unwrap();
    s.initialize().unwrap();
    let mut cp = s
        .checkpoint(CheckpointKind::ExactContinuation, m, e)
        .unwrap();
    let bytes = io::encode(&cp).unwrap();
    let len = u64::from_le_bytes(bytes[8..16].try_into().unwrap()) as usize;
    let header = std::str::from_utf8(&bytes[24..24 + len]).unwrap();
    assert!(!header.contains("start_wavefunctions"));
    assert!(!header.contains("atomic_start"));
    let old = io::decode(&bytes).unwrap();
    assert_eq!(
        old.continuation.unwrap().options.start_wavefunctions,
        StartWavefunctions::FreeElectron
    );
    let c = cp.continuation.as_mut().unwrap();
    c.options.start_wavefunctions = StartWavefunctions::Atomic;
    c.atomic_start = Some(AtomicStartReport {
        atomic_orbitals: 4,
        initial_h_applications: 5,
        gram_rotation_macs: 100,
        evd_orders: vec![5],
        ..Default::default()
    });
    let bytes = io::encode(&cp).unwrap();
    let decoded = io::decode(&bytes).unwrap();
    let c = decoded.continuation.as_ref().unwrap();
    assert_eq!(c.options.start_wavefunctions, StartWavefunctions::Atomic);
    assert_eq!(c.atomic_start.as_ref().unwrap().initial_h_applications, 5);
    assert_eq!(io::encode(&decoded).unwrap(), bytes);
    // Unknown option values fail even with a valid container checksum.
    let len = u64::from_le_bytes(bytes[8..16].try_into().unwrap()) as usize;
    let mut h: serde_json::Value = serde_json::from_slice(&bytes[24..24 + len]).unwrap();
    h["continuation"]["options"]["start_wavefunctions"] = serde_json::json!("unknown");
    let h = serde_json::to_vec(&h).unwrap();
    let mut bad = bytes[..8].to_vec();
    bad.extend_from_slice(&(h.len() as u64).to_le_bytes());
    bad.extend_from_slice(&bytes[16..24]);
    bad.extend_from_slice(&h);
    bad.extend_from_slice(&bytes[24 + len..bytes.len() - 32]);
    bad.extend_from_slice(&Sha256::digest(&bad));
    assert!(io::decode(&bad).is_err());
}
