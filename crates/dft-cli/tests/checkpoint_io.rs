#[allow(dead_code)]
#[path = "../../dft-engine/tests/common/mod.rs"]
mod common;
#[allow(dead_code)]
#[path = "../../dft-engine/tests/common/checkpoint.rs"]
mod fixture;
use dft_cli::checkpoint_io as io;
use dft_engine::{checkpoint::*, *};
use fixture::*;
use sha2::{Digest, Sha256};
fn cp() -> Checkpoint {
    let (m, e) = identity(0, 1);
    let mut s = ScfSession::new(
        fixture(0, Mixing::default(), 1, 200),
        CompletionMode::Strict,
        ExecutionConfig::default(),
    )
    .unwrap();
    s.initialize().unwrap();
    s.checkpoint(CheckpointKind::ExactContinuation, m, e)
        .unwrap()
}
fn rewrite_header(bytes: &[u8], f: impl FnOnce(&mut serde_json::Value)) -> Vec<u8> {
    let len = u64::from_le_bytes(bytes[8..16].try_into().unwrap()) as usize;
    let mut header: serde_json::Value = serde_json::from_slice(&bytes[24..24 + len]).unwrap();
    f(&mut header);
    let header = serde_json::to_vec(&header).unwrap();
    let mut out = bytes[..8].to_vec();
    out.extend_from_slice(&(header.len() as u64).to_le_bytes());
    out.extend_from_slice(&bytes[16..24]);
    out.extend_from_slice(&header);
    out.extend_from_slice(&bytes[24 + len..bytes.len() - 32]);
    let hash = Sha256::digest(&out);
    out.extend_from_slice(&hash);
    out
}
#[test]
fn layout_checksum_and_strict_header() {
    let cp = cp();
    let bytes = io::encode(&cp).unwrap();
    assert_eq!(&bytes[..8], b"BAANDICP");
    let h = u64::from_le_bytes(bytes[8..16].try_into().unwrap()) as usize;
    let b = u64::from_le_bytes(bytes[16..24].try_into().unwrap()) as usize;
    assert_eq!(24 + h + b + 32, bytes.len());
    let v: serde_json::Value = serde_json::from_slice(&bytes[24..24 + h]).unwrap();
    let offset = v["electronic_state"]["density"]["offset"].as_u64().unwrap() as usize;
    // Independent specification check, not the production reader.
    for (i, x) in cp.electronic_state.density.iter().enumerate() {
        assert_eq!(
            &bytes[24 + h + offset + 8 * i..24 + h + offset + 8 * i + 8],
            &x.to_le_bytes()
        );
    }
    assert_eq!(
        Sha256::digest(&bytes[..bytes.len() - 32])[..],
        bytes[bytes.len() - 32..]
    );
    assert_eq!(io::encode(&io::decode(&bytes).unwrap()).unwrap(), bytes);
    for index in [0, 8, 16, 24, 24 + h, bytes.len() - 1] {
        let mut bad = bytes.clone();
        bad[index] ^= 1;
        assert!(io::decode(&bad).is_err());
    }
    for size in [0, 8, 24, bytes.len() - 1] {
        assert!(io::decode(&bytes[..size]).is_err());
    }
    let mut trailing = bytes.clone();
    trailing.push(0);
    assert!(io::decode(&trailing).is_err());
    for bad in [
        rewrite_header(&bytes, |v| v["schema_version"] = 2.into()),
        rewrite_header(&bytes, |v| v["unknown"] = true.into()),
        rewrite_header(&bytes, |v| {
            v["electronic_state"]["density"]["offset"] = 8.into()
        }),
        rewrite_header(&bytes, |v| {
            v["electronic_state"]["density"]["count"] = u64::MAX.into()
        }),
        rewrite_header(&bytes, |v| {
            v["model"]["cutoff_hartree"] = f64::NAN.to_bits().into()
        }),
    ] {
        assert!(io::decode(&bad).is_err());
    }
    // Duplicate keys with a valid checksum must also reject.
    let mut head = String::from_utf8(bytes[24..24 + h].to_vec()).unwrap();
    head.insert_str(1, "\"schema_version\":1,");
    let mut dup = bytes[..8].to_vec();
    dup.extend_from_slice(&(head.len() as u64).to_le_bytes());
    dup.extend_from_slice(&bytes[16..24]);
    dup.extend_from_slice(head.as_bytes());
    dup.extend_from_slice(&bytes[24 + h..bytes.len() - 32]);
    let digest = Sha256::digest(&dup);
    dup.extend_from_slice(&digest);
    assert!(io::decode(&dup).is_err());
}
#[test]
fn atomic_no_clobber_publication() {
    let path = std::env::temp_dir().join(format!("p11-checkpoint-{}", std::process::id()));
    let cp = cp();
    io::write(&path, &cp).unwrap();
    assert!(io::write(&path, &cp).is_err());
    assert_eq!(
        io::encode(&io::read(&path).unwrap()).unwrap(),
        io::encode(&cp).unwrap()
    );
    std::fs::remove_file(path).unwrap();
}
#[test]
#[ignore = "release checkpoint file restart matrix"]
fn binary_roundtrip_restart() {
    for kind in 0..3 {
        for threads in [1, 4] {
            for mixing in [
                Mixing::Linear { beta: 0.3 },
                Mixing::default(),
                Mixing::KerkerPulay {
                    beta: 0.3,
                    history: 6,
                    q0_inverse_bohr: 1.0,
                },
            ] {
                let (m, e) = identity(kind, threads);
                let mut s = ScfSession::new(
                    fixture(kind, mixing.clone(), threads, 200),
                    CompletionMode::Strict,
                    ExecutionConfig::default(),
                )
                .unwrap();
                s.initialize().unwrap();
                for _ in 0..3 {
                    s.step().unwrap();
                }
                let cp = s
                    .checkpoint(CheckpointKind::ExactContinuation, m.clone(), e.clone())
                    .unwrap();
                let bytes = io::encode(&cp).unwrap();
                let cp = io::decode(&bytes).unwrap();
                let restored = ScfSession::restore_exact(
                    fixture(kind, mixing.clone(), threads, 200),
                    &cp,
                    &m,
                    &e,
                )
                .unwrap();
                let baseline = finish(s);
                assert_identical(&baseline, &finish(restored));
                println!(
                    "file restart kind={kind} threads={threads} mixing={mixing:?} bytes={} energy={:.17} Ha",
                    bytes.len(),
                    baseline.total_hartree()
                );
            }
        }
    }
}
#[test]
#[ignore = "release runtime checkpoint with committed policy and audit"]
fn runtime_policy_roundtrip() {
    use dft_engine::policy::PolicyChange;
    use dft_engine::runtime::*;
    use std::time::{Duration, SystemTime};
    for threads in [1, 4] {
        let (m, e) = identity(1, threads);
        let s = ScfSession::new(
            fixture(1, Mixing::default(), threads, 200),
            CompletionMode::Strict,
            ExecutionConfig {
                mode: ExecutionMode::AgentControlled,
                ..Default::default()
            },
        )
        .unwrap();
        let mut runtime = Runtime::new(s, RuntimeConfig::default()).unwrap();
        runtime.initialize().unwrap();
        runtime.step().unwrap();
        let snap = runtime.snapshot();
        let request = ControlRequest {
            schema_version: 1,
            calc_id: "scf".into(),
            id: "beta".into(),
            observed_state_version: snap.state_version,
            observed_policy_revision: snap.policy_revision,
            action: ControlAction::Policy(vec![PolicyChange::MixingBeta(0.2)]),
            preconditions: vec![],
            expiry: Expiry {
                max_iteration: Some(100),
                deadline: None,
            },
            stale_policy: StalePolicy::Reject,
            reason: "checkpoint fixture".into(),
        };
        runtime
            .submit(
                request.clone(),
                "test-host".into(),
                Authority::Approved(Approval {
                    id: "approval".into(),
                    approver: "test-host".into(),
                    expires_at: SystemTime::now() + Duration::from_secs(600),
                    allow_revalidation: false,
                }),
            )
            .unwrap();
        runtime.safe_point().unwrap();
        for _ in 0..2 {
            runtime.step().unwrap();
        }
        let mut pause = request.clone();
        pause.id = "pause".into();
        pause.action = ControlAction::Pause;
        pause.observed_state_version = runtime.snapshot().state_version;
        pause.observed_policy_revision = runtime.policy_revision();
        runtime
            .submit(
                pause,
                "test-host".into(),
                Authority::Approved(Approval {
                    id: "pause-approval".into(),
                    approver: "test-host".into(),
                    expires_at: SystemTime::now() + Duration::from_secs(600),
                    allow_revalidation: false,
                }),
            )
            .unwrap();
        runtime.safe_point().unwrap();
        assert_eq!(runtime.state(), RuntimeState::Paused);
        while runtime.events().unwrap().pop().is_some() {}
        let saved = runtime.checkpoint(m.clone(), e.clone()).unwrap();
        let bytes = io::encode(&saved).unwrap();
        let saved = io::decode(&bytes).unwrap();
        let mut restored = Runtime::restore_exact(
            fixture(
                1,
                Mixing::Pulay {
                    beta: 0.2,
                    history: 6,
                },
                threads,
                200,
            ),
            &saved,
            &m,
            &e,
        )
        .unwrap();
        assert_eq!(restored.policy_revision(), 1);
        assert_eq!(
            format!("{:?}", runtime.audit()),
            format!("{:?}", restored.audit())
        );
        for r in [&mut runtime, &mut restored] {
            assert_eq!(r.state(), RuntimeState::Paused);
            let mut resume = request.clone();
            resume.id = "resume".into();
            resume.action = ControlAction::Resume;
            resume.observed_state_version = r.snapshot().state_version;
            resume.observed_policy_revision = r.policy_revision();
            r.submit(
                resume,
                "test-host".into(),
                Authority::Approved(Approval {
                    id: "resume-approval".into(),
                    approver: "test-host".into(),
                    expires_at: SystemTime::now() + Duration::from_secs(600),
                    allow_revalidation: false,
                }),
            )
            .unwrap();
            r.safe_point().unwrap();
            while !r.state().terminal() {
                r.step().unwrap();
                while r.events().unwrap().pop().is_some() {}
            }
        }
        let a = runtime.into_result();
        let b = restored.into_result();
        assert_eq!(a.policy_changes, b.policy_changes);
        assert_identical(
            &normalized(a.outcome.unwrap().result),
            &normalized(b.outcome.unwrap().result),
        );
        println!(
            "runtime file restart threads={threads} audit_records={} bytes={}",
            saved.runtime.as_ref().unwrap().audit.len(),
            bytes.len()
        );
    }
}

#[test]
fn result_export_does_not_confer_restart_or_convergence() {
    let (m, e) = identity(0, 1);
    let mut s = ScfSession::new(
        fixture(0, Mixing::default(), 1, 1),
        CompletionMode::RetainLast,
        ExecutionConfig::default(),
    )
    .unwrap();
    s.initialize().unwrap();
    s.step().unwrap();
    let cp = s
        .checkpoint(CheckpointKind::ResultExport, m.clone(), e.clone())
        .unwrap();
    let copy = io::decode(&io::encode(&cp).unwrap()).unwrap();
    assert!(!copy.result.unwrap().converged);
    assert!(cp.continuation.is_none());
    assert_eq!(
        ScfSession::restore_exact(fixture(0, Mixing::default(), 1, 200), &cp, &m, &e)
            .err()
            .unwrap(),
        Rejection::Kind
    );
    assert_eq!(
        ScfSession::warm_start(
            fixture(0, Mixing::default(), 1, 200),
            &cp,
            &m,
            &e,
            CompletionMode::Strict,
            ExecutionConfig {
                calc_id: "new".into(),
                ..Default::default()
            }
        )
        .err()
        .unwrap(),
        Rejection::Kind
    );
}
#[test]
#[ignore = "release warm-start file convergence"]
fn warm_start_file() {
    for kind in 0..3 {
        let (m, e) = identity(kind, 4);
        let mut source = ScfSession::new(
            fixture(kind, Mixing::Linear { beta: 0.3 }, 4, 200),
            CompletionMode::Strict,
            ExecutionConfig::default(),
        )
        .unwrap();
        source.initialize().unwrap();
        for _ in 0..3 {
            source.step().unwrap();
        }
        let cp = source
            .checkpoint(CheckpointKind::WarmStart, m.clone(), e.clone())
            .unwrap();
        let cp = io::decode(&io::encode(&cp).unwrap()).unwrap();
        let mut baseline = ScfSession::new(
            fixture(kind, Mixing::default(), 4, 200),
            CompletionMode::Strict,
            ExecutionConfig::default(),
        )
        .unwrap();
        baseline.initialize().unwrap();
        let baseline = finish(baseline);
        let restored = ScfSession::warm_start(
            fixture(kind, Mixing::default(), 4, 200),
            &cp,
            &m,
            &e,
            CompletionMode::Strict,
            ExecutionConfig {
                calc_id: "warm".into(),
                ..Default::default()
            },
        )
        .unwrap();
        let result = finish(restored);
        let delta = (baseline.total_hartree() - result.total_hartree()).abs();
        assert!(delta < 1e-10);
        println!("warm file kind={kind} threads=4 delta={delta:.17e} Ha");
    }
}

#[allow(dead_code)]
#[path = "common/mod.rs"]
mod cli_common;
#[test]
#[ignore = "release cross-process CLI restart"]
fn command_line_cross_process() {
    use std::{fs, process::Command};
    let dir = std::env::temp_dir().join(format!("p11-cli-process-{}", std::process::id()));
    fs::create_dir(&dir).unwrap();
    fs::write(
        dir.join("synthetic.xml"),
        cli_common::fixture().replace("element=\"Xx\"", "element=\"Si\""),
    )
    .unwrap();
    let text = "[system]\nname='checkpoint fixture'\nlattice_unit='bohr'\nlattice_vectors=[[4.0,0.0,0.0],[0.0,4.0,0.0],[0.0,0.0,4.0]]\n[[system.atoms]]\nspecies='Si'\nfractional=[0.13,0.21,0.32]\n[pseudopotentials]\nSi='synthetic.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n";
    fs::write(dir.join("input.toml"), text).unwrap();
    let invoke = |args: &[&str]| {
        let o = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .current_dir(&dir)
            .args(args)
            .output()
            .unwrap();
        assert!(
            o.status.success(),
            "{:?}: {} {}",
            args,
            String::from_utf8_lossy(&o.stdout),
            String::from_utf8_lossy(&o.stderr)
        );
    };
    invoke(&["checkpoint", "input.toml", "exact.cp", "--iteration", "3"]);
    invoke(&["restart", "input.toml", "exact.cp", "restarted.cp"]);
    invoke(&[
        "checkpoint",
        "input.toml",
        "baseline.cp",
        "--kind",
        "result",
    ]);
    let mut a = io::read(&dir.join("baseline.cp")).unwrap();
    let mut b = io::read(&dir.join("restarted.cp")).unwrap();
    for cp in [&mut a, &mut b] {
        for h in &mut cp.result.as_mut().unwrap().history {
            h.timings = PhaseTimings::default();
            if let Some(c) = &mut h.fixed_point_certificate {
                c.timings = PhaseTimings::default();
            }
        }
    }
    assert_eq!(
        format!("{:?}", a.electronic_state),
        format!("{:?}", b.electronic_state)
    );
    assert_eq!(format!("{:?}", a.result), format!("{:?}", b.result));
    assert_eq!(a.execution, b.execution);
    invoke(&[
        "checkpoint",
        "input.toml",
        "warm.cp",
        "--kind",
        "warm",
        "--iteration",
        "3",
    ]);
    invoke(&[
        "restart",
        "input.toml",
        "warm.cp",
        "warm-result.cp",
        "--warm",
    ]);
    let mut corrupt = fs::read(dir.join("exact.cp")).unwrap();
    corrupt[24] ^= 1;
    fs::write(dir.join("bad.cp"), corrupt).unwrap();
    let failed = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
        .current_dir(&dir)
        .args(["restart", "input.toml", "bad.cp", "bad-result.cp"])
        .output()
        .unwrap();
    assert!(!failed.status.success());
    assert!(!dir.join("bad-result.cp").exists());
    println!(
        "cross-process CLI energy={:.17} Ha",
        a.result.unwrap().energy.total_hartree()
    );
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn v1_codec_rejects_nlcc_instead_of_dropping_its_ledger() {
    let (m, e) = identity(0, 1);
    let mut s = ScfSession::new(
        fixture(0, Mixing::default(), 1, 200),
        CompletionMode::Strict,
        ExecutionConfig::default(),
    )
    .unwrap();
    s.initialize().unwrap();
    s.step().unwrap();
    let mut checkpoint = s
        .checkpoint(CheckpointKind::ExactContinuation, m, e)
        .unwrap();
    let bytes = io::encode(&checkpoint).unwrap();
    assert_eq!(io::encode(&io::decode(&bytes).unwrap()).unwrap(), bytes);
    checkpoint.continuation.as_mut().unwrap().history[0]
        .energy
        .nlcc = Some(dft_engine::energy::CoreXcLedger {
        core_electrons: 1.,
        core_potential_integral_ha: -0.1,
        valence_potential_integral_ha: -0.2,
    });
    assert!(
        io::encode(&checkpoint)
            .unwrap_err()
            .to_string()
            .contains("NLCC")
    );
}
