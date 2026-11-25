//! Reference v2 file identity against an uninterrupted trajectory.
//! Full matrix: cargo test --release -p dft-cli --test checkpoint_v2 -- --include-ignored --nocapture
#[allow(dead_code)]
#[path = "../../dft-engine/tests/common/mod.rs"]
mod common;
use dft_cli::checkpoint_io as io;
use dft_engine::{
    GammaReal, ScfOptions, SpeciesPotential, Verification,
    refscf::{checkpoint::Identity, driver::Run, k_driver, profile::Profile, stepper::Stepper},
    refstart::StartInput,
};
use dft_geometry::{
    Atom, Cell,
    kpoints::{KPoint, monkhorst_pack},
};
use dft_mpi::Collective;
use sha2::{Digest, Sha256};
use std::fs;

fn identity(ranks: usize) -> Identity {
    Identity {
        input_sha256: "1".repeat(64),
        upf_sha256: vec!["2".repeat(64)],
        binary_sha256: "3".repeat(64),
        platform: "synthetic-test".into(),
        threads: 1,
        ranks,
    }
}
fn bits(run: &Run) -> (u64, Vec<u64>, Vec<u64>) {
    let r = &run.outcome.result;
    (
        r.energy.total_hartree().to_bits(),
        r.density.iter().map(|v| v.to_bits()).collect(),
        r.k_points
            .iter()
            .flat_map(|k| k.solution.eigenvectors.as_slice())
            .flat_map(|v| [v.re.to_bits(), v.im.to_bits()])
            .collect(),
    )
}
fn hash(v: &[u64]) -> String {
    format!(
        "{:016x}",
        v.iter()
            .flat_map(|v| v.to_le_bytes())
            .fold(0xcbf29ce484222325u64, |h, b| (h ^ u64::from(b))
                .wrapping_mul(0x100000001b3))
    )
}
fn boundary(s: &mut Stepper<'_, Run>, iteration: usize) {
    s.initialize().unwrap();
    if iteration == 0 {
        return;
    }
    loop {
        let r = s.step().unwrap();
        assert!(!s.is_done(), "terminal before requested checkpoint");
        if r.iteration >= iteration && r.safe_point {
            break;
        }
    }
}
fn fixture(mode: &str, ranks: usize, budget: usize, points_to_save: &[usize], reuse: bool) {
    let skew = mode == "skew-two-species";
    let cell = Cell::from_user_rows_bohr(if skew {
        [[6., 0., 0.], [0.6, 5.5, 0.], [0.1, 0.3, 6.5]]
    } else {
        [[5., 0., 0.], [0., 5., 0.], [0., 0., 5.]]
    })
    .unwrap();
    let mut atoms = vec![Atom::new("Xx", [0.13, 0.21, 0.32])];
    if skew {
        atoms.push(Atom::new("Yy", [0.63, 0.71, 0.72]));
    }
    let pp = common::potential_with_atomic_start(2., 1., true, 161, "SLA PZ");
    let pp2 = common::potential_with_atomic_start(2., 1.2, true, 161, "SLA PZ");
    let potentials = [
        SpeciesPotential {
            label: "Xx",
            potential: &pp,
        },
        SpeciesPotential {
            label: "Yy",
            potential: &pp2,
        },
    ];
    let options = ScfOptions {
        cutoff_hartree: 2.5,
        bands: Some(if skew { 6 } else { 4 }),
        gamma_real: if mode == "k222s" {
            GammaReal::Off
        } else {
            GammaReal::On
        },
        max_iterations: budget,
        verification: Verification::FixedPoint,
        ..Default::default()
    };
    let points = if mode == "k222s" {
        monkhorst_pack([2; 3], [0.5; 3]).unwrap()
    } else {
        vec![KPoint {
            fractional: [0.; 3],
            weight: 1.,
        }]
    };
    dft_mpi::simulate(ranks, |comm| {
        let _scope = dft_numerics::switches::RunConfig::from_pairs(&[])
            .unwrap()
            .enter();
        let build = || {
            k_driver::stepper_with(
                StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &options,
                },
                Profile::default(),
                &points,
                comm,
            )
        };
        let uninterrupted = build().run_to_end().unwrap();
        let expected = bits(&uninterrupted);
        if budget > 10 {
            assert!(
                uninterrupted.outcome.converged,
                "{mode}: {}",
                uninterrupted.exit_reason
            );
        }
        for &at in points_to_save {
            let mut s = build();
            boundary(&mut s, at);
            let cp = s.checkpoint(identity(ranks), comm.rank()).unwrap();
            drop(s);
            let file = io::gather_reference(cp, comm).unwrap();
            let path = std::env::temp_dir().join(format!(
                "r8b-v2-{}-{mode}-{ranks}-{budget}-{at}-{reuse}.cp",
                std::process::id()
            ));
            if comm.rank() == 0 {
                io::write_reference(&path, &file).unwrap();
                assert!(
                    io::write_reference(&path, &file).is_err(),
                    "clobbered immutable checkpoint"
                );
                let bytes = fs::read(&path).unwrap();
                let h = u64::from_le_bytes(bytes[8..16].try_into().unwrap()) as usize;
                let b = u64::from_le_bytes(bytes[16..24].try_into().unwrap()) as usize;
                assert_eq!(bytes.len(), 56 + h + b);
                assert_eq!(
                    Sha256::digest(&bytes[..bytes.len() - 32])[..],
                    bytes[bytes.len() - 32..]
                );
                let header: serde_json::Value = serde_json::from_slice(&bytes[24..24 + h]).unwrap();
                assert_eq!(header["schema_version"], 2);
                let mut end = 0;
                for (rank, r) in header["states"].as_array().unwrap().iter().enumerate() {
                    let offset = r["offset"].as_u64().unwrap() as usize;
                    let count = r["count"].as_u64().unwrap() as usize;
                    assert_eq!(offset, end);
                    assert_eq!(
                        &bytes[24 + h + offset..24 + h + offset + count],
                        file.states[rank]
                    );
                    end += count;
                }
                assert_eq!(end, b);
            }
            comm.barrier().unwrap();
            let decoded = io::read_reference(&path).unwrap();
            let cp = decoded.checkpoint(comm.rank()).unwrap();
            let mut restored = build();
            restored
                .restore(&cp, &identity(ranks), comm.rank())
                .unwrap();
            comm.barrier().unwrap();
            if comm.rank() == 0 {
                fs::remove_file(&path).unwrap();
            }
            let resumed = restored.run_to_end().unwrap();
            assert_eq!(
                bits(&resumed),
                expected,
                "restart at {at}, {mode}, ranks {ranks}"
            );
            assert_eq!(
                resumed.outcome.result.history.len(),
                uninterrupted.outcome.result.history.len()
            );
            if comm.rank() == 0 {
                let row = format!(
                    "{{\"fixture\":\"{mode}\",\"ranks\":{ranks},\"budget\":{budget},\"converged\":{},\"safe_point\":{},\"reuse\":{reuse},\"energy_ha\":{:.17e},\"energy_bits\":\"{:016x}\",\"density_fnv1a\":\"{}\",\"orbitals_fnv1a\":\"{}\",\"identity\":true}}",
                    uninterrupted.outcome.converged,
                    cp.iteration,
                    f64::from_bits(expected.0),
                    expected.0,
                    hash(&expected.1),
                    hash(&expected.2)
                );
                println!("R8B_FILE {row}");
                if let Ok(path) = std::env::var("BAANDI_CHECKPOINT_FILE_EVIDENCE") {
                    use std::io::Write;
                    writeln!(
                        std::fs::OpenOptions::new()
                            .create(true)
                            .append(true)
                            .open(path)
                            .unwrap(),
                        "{row}"
                    )
                    .unwrap();
                }
            }
        }
    });
}
#[test]
fn short_file_identity() {
    fixture("gamma", 1, 6, &[0, 1, 3], false);
}
#[test]
#[ignore = "full convergence matrix; run in release (module header)"]
fn file_gamma_k222s_skew_ranks_1_3_8() {
    for mode in ["gamma", "k222s", "skew-two-species"] {
        for ranks in [1, 3, 8] {
            fixture(mode, ranks, 400, &[0, 1, 3], false);
        }
    }
}

#[test]
fn rejects_identity_and_corruption() {
    let cp = dft_engine::refscf::checkpoint::Checkpoint {
        schema_version: 2,
        identity: identity(1),
        rank: 0,
        iteration: 0,
        state: vec![0; 8],
    };
    let mut wrong = identity(1);
    wrong.input_sha256 = "4".repeat(64);
    assert!(
        cp.validate(&wrong, 0)
            .unwrap_err()
            .to_string()
            .contains("identity mismatch")
    );
    assert!(cp.validate(&identity(1), 1).is_err());
    let mut bad = cp.clone();
    bad.schema_version = 3;
    assert!(bad.validate(&identity(1), 0).is_err());
}

#[test]
fn container_corruption_and_strict_schema() {
    let file = io::ReferenceFile {
        identity: identity(1),
        iteration: 0,
        states: vec![vec![0; 32]],
        result: None,
    };
    let bytes = io::encode_reference(&file).unwrap();
    for index in [0, 8, 16, 24, bytes.len() - 33, bytes.len() - 1] {
        let mut bad = bytes.clone();
        bad[index] ^= 1;
        assert!(io::decode_reference(&bad).is_err());
    }
    for len in [0, 8, 24, bytes.len() - 1] {
        assert!(io::decode_reference(&bytes[..len]).is_err());
    }
    let h = u64::from_le_bytes(bytes[8..16].try_into().unwrap()) as usize;
    for field in ["version", "offset", "count", "rank", "unknown", "identity"] {
        let mut header: serde_json::Value = serde_json::from_slice(&bytes[24..24 + h]).unwrap();
        match field {
            "version" => header["schema_version"] = 99.into(),
            "offset" => header["states"][0]["offset"] = 8.into(),
            "count" => header["states"][0]["count"] = u64::MAX.into(),
            "rank" => header["identity"]["ranks"] = 2.into(),
            "unknown" => header["unknown"] = true.into(),
            _ => header["identity"]["binary_sha256"] = "invalid".into(),
        }
        let header = serde_json::to_vec(&header).unwrap();
        let mut bad = bytes[..8].to_vec();
        bad.extend_from_slice(&(header.len() as u64).to_le_bytes());
        bad.extend_from_slice(&bytes[16..24]);
        bad.extend(header);
        bad.extend_from_slice(&bytes[24 + h..bytes.len() - 32]);
        let hash = Sha256::digest(&bad);
        bad.extend_from_slice(&hash);
        assert!(io::decode_reference(&bad).is_err(), "{field}");
    }
    let mut header = String::from_utf8(bytes[24..24 + h].to_vec()).unwrap();
    header.insert_str(1, "\"schema_version\":2,");
    let mut bad = bytes[..8].to_vec();
    bad.extend_from_slice(&(header.len() as u64).to_le_bytes());
    bad.extend_from_slice(&bytes[16..24]);
    bad.extend(header.as_bytes());
    bad.extend_from_slice(&bytes[24 + h..bytes.len() - 32]);
    let hash = Sha256::digest(&bad);
    bad.extend_from_slice(&hash);
    assert!(io::decode_reference(&bad).is_err());
}
#[allow(dead_code)]
#[path = "common/mod.rs"]
mod cli_common;
fn atomic_text() -> String {
    let mut text = cli_common::gaussian_with(2., 1., true, 161);
    let p = dft_upf::parse_upf(&text).unwrap();
    let values = p
        .mesh()
        .r()
        .iter()
        .map(|r| format!("{:.17e}", r * (-r * r).exp()))
        .collect::<Vec<_>>()
        .join(" ");
    let a = text.find("<PP_RHOATOM>").unwrap() + "<PP_RHOATOM>".len();
    let b = text.find("</PP_RHOATOM>").unwrap();
    text.replace_range(a..b, &values);
    text.replace("number_of_wfc=\"0\"","number_of_wfc=\"1\"").replace("</UPF>",&format!("<PP_PSWFC><PP_CHI.1 l=\"0\" label=\"s\" occupation=\"2\">{values}</PP_CHI.1></PP_PSWFC></UPF>"))
}
#[test]
fn cross_process_restart_inspect_and_input_rejection() {
    let dir = std::env::temp_dir().join(format!("r8b-process-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(dir.join("synthetic.upf"), atomic_text()).unwrap();
    let input = r#"[system]
name="R8b synthetic"
lattice_unit="bohr"
lattice_vectors=[[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0]]
[[system.atoms]]
species="Xx"
fractional=[0.13,0.21,0.32]
[pseudopotentials]
Xx="synthetic.upf"
[basis]
ecutwfc_ry=5.0
kpoints=[[0.0,0.0,0.0]]
[xc]
functional="lda-pz"
[scf]
max_iterations=200
nbands=4
"#;
    fs::write(dir.join("input.toml"), input).unwrap();
    let invoke = |args: &[&str]| {
        let mut c = std::process::Command::new(env!("CARGO_BIN_EXE_dft-cli"));
        for (k, _) in std::env::vars().filter(|(k, _)| k.starts_with("BAANDI_")) {
            c.env_remove(k);
        }
        c.current_dir(&dir).args(args).output().unwrap()
    };
    for args in [
        vec!["checkpoint", "input.toml", "full.cp", "--kind", "result"],
        vec!["checkpoint", "input.toml", "part.cp", "--iteration", "3"],
        vec!["restart", "input.toml", "part.cp", "resumed.cp"],
    ] {
        let output = invoke(&args);
        assert!(
            output.status.success(),
            "{:?}: {}",
            args,
            String::from_utf8_lossy(&output.stderr)
        );
    }
    let a = io::read_reference(&dir.join("full.cp"))
        .unwrap()
        .result
        .unwrap();
    let b = io::read_reference(&dir.join("resumed.cp"))
        .unwrap()
        .result
        .unwrap();
    assert_eq!(a.energy_ha_bits, b.energy_ha_bits);
    assert_eq!(a.density_sha256, b.density_sha256);
    assert_eq!(a.orbitals_sha256, b.orbitals_sha256);
    assert_eq!(a.verification_attempts, b.verification_attempts);
    fs::write(dir.join("changed.toml"), input.replace("0.13", "0.14")).unwrap();
    let out = invoke(&["restart", "changed.toml", "part.cp", "wrong.cp"]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("identity mismatch"));
    assert!(!dir.join("wrong.cp").exists());
    let upf = fs::read_to_string(dir.join("synthetic.upf")).unwrap();
    fs::write(
        dir.join("synthetic.upf"),
        format!("{upf}\n<!-- identity change -->\n"),
    )
    .unwrap();
    let out = invoke(&["restart", "input.toml", "part.cp", "wrong-upf.cp"]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("identity mismatch"));
    fs::write(dir.join("synthetic.upf"), upf).unwrap();
    let mut wrong_binary = io::read_reference(&dir.join("part.cp")).unwrap();
    wrong_binary.identity.binary_sha256 = "0".repeat(64);
    io::write_reference(&dir.join("wrong-binary.cp"), &wrong_binary).unwrap();
    let out = invoke(&["restart", "input.toml", "wrong-binary.cp", "wrong-build.cp"]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("identity mismatch"));
    let mut data = fs::read(dir.join("part.cp")).unwrap();
    let n = data.len();
    data[n - 33] ^= 1;
    fs::write(dir.join("corrupt.cp"), data).unwrap();
    let out = invoke(&["restart", "input.toml", "corrupt.cp", "bad.cp"]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("checksum"));
    println!("R8B_CROSS_PROCESS {}", serde_json::to_string(&a).unwrap());
    let out = invoke(&["explain", "part.cp"]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(String::from_utf8_lossy(&out.stdout).contains('2'));
    fs::remove_dir_all(dir).unwrap();
}
