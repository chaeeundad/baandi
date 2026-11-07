#[allow(dead_code)]
#[path = "common/mod.rs"]
mod common;
use std::{fs, path::PathBuf, process::Command};

/// DFT_REFERENCE_BASELINE=$PWD/P44A-evidence/baseline-cli cargo test --release
/// -p dft-cli --test reference_mode_relaxation -- --ignored --nocapture
#[test]
#[ignore = "release-only baseline versus candidate ionic and cell trajectories"]
fn no_mode_relaxation_outputs_are_byte_identical() {
    let Some(baseline) = std::env::var_os("DFT_REFERENCE_BASELINE").map(PathBuf::from) else {
        eprintln!("SKIP: set DFT_REFERENCE_BASELINE to the pre-change executable");
        return;
    };
    assert!(baseline.is_file());
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let dir = root
        .join("P44A-evidence/relaxation-identity")
        .join(std::process::id().to_string());
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("syn.xml"),
        common::gaussian_with(1., 1.5, false, 801).replace("element=\"Xx\"", "element=\"H\""),
    )
    .unwrap();
    let base = "[system]\nname='relaxation oracle'\nlattice_unit='bohr'\nlattice_vectors=[[8.0,0.0,0.0],[0.0,8.0,0.0],[0.0,0.0,8.0]]\n[[system.atoms]]\nspecies='H'\nfractional=[0.0,0.0,0.0]\n[[system.atoms]]\nspecies='H'\nfractional=[0.27,0.23,0.19]\n[pseudopotentials]\nH='syn.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[scf]\nnbands=4\nmax_iterations=300\nenergy_tolerance_ha=1e-12\ndensity_residual_tolerance=1e-10\neigen_tolerance=1e-11\n[relax]\nmax_steps=4\nforce_tol_ha_bohr=1e-8\n";
    for (case, extra) in [
        (
            "ions",
            "fixed=[{atom=0},{atom=1,components=[false,true,false]}]\n",
        ),
        ("cell", "mode='ions_cell'\ncell_dofs='isotropic'\n"),
    ] {
        let inp = dir.join(format!("{case}.toml"));
        fs::write(&inp, format!("{base}{extra}")).unwrap();
        for (label, binary) in [
            ("baseline", baseline.clone()),
            ("candidate", PathBuf::from(env!("CARGO_BIN_EXE_dft-cli"))),
        ] {
            let out = dir.join(format!("{case}-{label}"));
            let output = Command::new(&binary)
                .arg("relax")
                .arg(&inp)
                .arg("--output-dir")
                .arg(&out)
                .output()
                .unwrap();
            fs::write(dir.join(format!("{case}-{label}.stdout")), &output.stdout).unwrap();
            fs::write(dir.join(format!("{case}-{label}.stderr")), &output.stderr).unwrap();
            assert_eq!(
                output.status.code(),
                Some(3),
                "{}",
                String::from_utf8_lossy(&output.stderr)
            );
        }
        // Build provenance is expected to differ across executables. Verify its
        // exact derivation before comparing all remaining bytes, including every
        // floating-point state bit and optimizer history in the restart payload.
        let mut results = Vec::new();
        let mut restarts = Vec::new();
        for (label, binary) in [
            ("baseline", baseline.clone()),
            ("candidate", PathBuf::from(env!("CARGO_BIN_EXE_dft-cli"))),
        ] {
            let digest = dft_cli::run::sha256(&fs::read(binary).unwrap());
            let mut identity = fs::read(&inp).unwrap();
            identity.extend(digest.as_bytes());
            identity.extend(std::env::consts::ARCH.as_bytes());
            identity.extend(std::env::consts::OS.as_bytes());
            identity
                .extend(dft_cli::run::sha256(&fs::read(dir.join("syn.xml")).unwrap()).as_bytes());
            let identity = dft_cli::run::sha256(&identity);
            let out = dir.join(format!("{case}-{label}"));
            let result = fs::read_to_string(out.join("results.json")).unwrap();
            let parsed: serde_json::Value = serde_json::from_str(&result).unwrap();
            assert_eq!(parsed["provenance"]["executable_sha256"], digest);
            results.push(result.replace(&digest, "VERIFIED_EXECUTABLE_SHA256"));
            let envelope: serde_json::Value =
                serde_json::from_slice(&fs::read(out.join("restart.json")).unwrap()).unwrap();
            let payload = envelope["payload"].as_str().unwrap();
            assert_eq!(envelope["sha256"], dft_cli::run::sha256(payload.as_bytes()));
            let parsed: serde_json::Value = serde_json::from_str(payload).unwrap();
            assert_eq!(
                if case == "cell" {
                    &parsed["ionic"]["identity"]
                } else {
                    &parsed["identity"]
                },
                &identity
            );
            restarts.push(payload.replace(&identity, "VERIFIED_BUILD_INPUT_IDENTITY"));
        }
        assert_eq!(
            results[0], results[1],
            "{case}: result bytes excluding verified executable hash"
        );
        assert_eq!(
            restarts[0], restarts[1],
            "{case}: restart payload bytes excluding verified build/input identity"
        );
        println!(
            "{case}: result and restart payload byte-identical after verified build provenance normalization"
        );
        let mut files = vec!["trajectory.extxyz", "POSCAR", "final.toml", "relax.log"];
        if case == "cell" {
            files.push("final.cif");
        }
        for file in files {
            assert_eq!(
                fs::read(dir.join(format!("{case}-baseline")).join(file)).unwrap(),
                fs::read(dir.join(format!("{case}-candidate")).join(file)).unwrap(),
                "{case}: {file}"
            );
            println!("{case}: {file} byte-identical");
        }
    }
}
