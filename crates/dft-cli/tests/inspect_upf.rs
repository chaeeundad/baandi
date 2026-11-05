//! Process-boundary tests; synthetic XML is generated at runtime only.
use std::process::Command;
fn fixture() -> String {
    let n = 401;
    let r = (0..n)
        .map(|i| format!("{}", i as f64 * 0.025))
        .collect::<Vec<_>>()
        .join(" ");
    let rab = vec!["0.025"; n].join(" ");
    // Smooth Coulomb-tailed synthetic potential, in Ry.
    let local = (0..n)
        .map(|i| format!("{}", -8.0 / ((i as f64 * 0.025).powi(2) + 1.0).sqrt()))
        .collect::<Vec<_>>()
        .join(" ");
    let beta = (0..n)
        .map(|i| {
            let r = i as f64 * 0.025;
            format!("{}", r * (-r * r).exp())
        })
        .collect::<Vec<_>>()
        .join(" ");
    let zero = vec!["0"; n].join(" ");
    format!(
        r#"<UPF version="2.0.1"><PP_INFO>SYNTHETIC CLI TEST ONLY</PP_INFO>
<PP_HEADER generated="synthetic" author="tests" date="2026" comment="synthetic" element="Xx"
pseudo_type="NC" relativistic="nonrelativistic" is_ultrasoft="F" is_paw="F" is_coulomb="F"
has_so="F" has_wfc="F" has_gipaw="F" core_correction="F" functional="SLA PZ"
z_valence="4" total_psenergy="0" wfc_cutoff="0" rho_cutoff="0" l_max="0" l_max_rho="0"
l_local="-1" mesh_size="{n}" number_of_wfc="0" number_of_proj="1"/>
<PP_MESH dx="0.025" mesh="{n}" xmin="0" rmax="10" zmesh="1"><PP_R>{r}</PP_R><PP_RAB>{rab}</PP_RAB></PP_MESH>
<PP_LOCAL>{local}</PP_LOCAL><PP_NONLOCAL><PP_BETA.1 index="1" angular_momentum="0" cutoff_radius="10" cutoff_radius_index="401">{beta}</PP_BETA.1><PP_DIJ>2</PP_DIJ></PP_NONLOCAL><PP_RHOATOM>{zero}</PP_RHOATOM></UPF>"#
    )
}
#[test]
fn inspect_synthetic_success_and_failure() {
    let dir = std::env::temp_dir().join(format!("dft-inspect-synthetic-{}", std::process::id()));
    std::fs::create_dir(&dir).unwrap();
    let path = dir.join("synthetic.xml");
    std::fs::write(&path, fixture()).unwrap();
    let run = || {
        Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .arg("inspect-upf")
            .arg(&path)
            .output()
            .unwrap()
    };
    let result = run();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let output = String::from_utf8(result.stdout).unwrap();
    for label in [
        "Element: Xx",
        "z_valence: 4",
        "Functional: SLA PZ",
        "Type: NC",
        "401 points",
        "Bohr",
        "Projectors: 1",
        "l=0, cutoff index=401",
        "D matrix: 1 x 1, Ha, row-major [1.0]",
        "Support: SUPPORTED",
    ] {
        assert!(output.contains(label), "missing {label}");
    }
    std::fs::write(
        &path,
        fixture().replace("core_correction=\"F\"", "core_correction=\"T\""),
    )
    .unwrap();
    let result = run();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("invalid field PP_NLCC"));
    std::fs::write(&path, "<UPF").unwrap();
    let result = run();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("malformed UPF"));
    std::fs::remove_file(&path).unwrap();
    assert!(!run().status.success());
    std::fs::remove_dir(dir).unwrap();
}

#[test]
fn json_inspection_lists_independent_blockers_and_ordered_channels() {
    let dir = std::env::temp_dir().join(format!("dft-inspect-blockers-{}", std::process::id()));
    std::fs::create_dir(&dir).unwrap();
    let path = dir.join("synthetic.upf");
    for (text, expected) in [
        (fixture().replace("is_paw=\"F\"", "is_paw=\"T\""), "PAW"),
        (
            fixture().replace("is_ultrasoft=\"F\"", "is_ultrasoft=\"T\""),
            "USPP",
        ),
        (
            fixture().replace("relativistic=\"nonrelativistic\"", "relativistic=\"full\""),
            "full-relativistic/SOC",
        ),
        (fixture().replace("has_so=\"F\"", "has_so=\"T\""), "SOC"),
        (fixture().replace("index=\"1\"", "index=\"2\""), "malformed"),
        (
            fixture().replace("angular_momentum=\"0\"", "angular_momentum=\"bad\""),
            "invalid index/angular_momentum",
        ),
        (
            fixture().replace("<PP_DIJ>2</PP_DIJ>", ""),
            "missing coupling matrix",
        ),
    ] {
        std::fs::write(&path, text).unwrap();
        let result = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .args(["inspect-upf", path.to_str().unwrap(), "--json"])
            .output()
            .unwrap();
        assert!(result.status.success()); // Inspection succeeded; SCF is blocked.
        let report: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
        assert_eq!(report["supported_by_engine"], false);
        assert!(
            report["blockers"].to_string().contains(expected),
            "{report}"
        );
        assert_eq!(report["header_facts"]["mesh_size"], 401);
        assert_eq!(report["sha256"].as_str().unwrap().len(), 64);
        assert!(report["license_reference"].is_null());
    }
    std::fs::write(&path, fixture()).unwrap();
    let report = inspection_report(&path, None).unwrap();
    assert_eq!(
        report["header_facts"]["projector_channels"],
        serde_json::json!([{"index":1,"l":0}])
    );
    let manifest = dir.join("manifest.json");
    let entry = serde_json::json!({"datasets":[{"element":"Xx", "sha256":report["sha256"], "license":{"identifier":"synthetic"},"qualification_status":"unqualified"}]});
    std::fs::write(&manifest, entry.to_string()).unwrap();
    let matched = inspection_report(&path, Some(&manifest)).unwrap();
    assert_eq!(matched["license_reference"]["identifier"], "synthetic");
    std::fs::write(&path, fixture().replace("SLA PZ", "PBE")).unwrap();
    assert!(inspection_report(&path, Some(&manifest)).is_err());
    std::fs::remove_dir_all(dir).unwrap();
}

/// Header inspection only; SCF smoke is in the engine real_nc test.
/// DFT_NC_DIR=/external/path cargo test --release -p dft-cli --test inspect_upf real_nc_headers -- --ignored --nocapture
#[test]
#[ignore = "requires external licensed UPFs via DFT_NC_DIR"]
fn real_nc_headers() {
    let dir = std::path::PathBuf::from(std::env::var_os("DFT_NC_DIR").expect("set DFT_NC_DIR"));
    let manifest =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/data/nc-datasets.json");
    for element in ["O", "Ni", "Al", "Mg", "Si"] {
        let start = std::time::Instant::now();
        let report =
            inspection_report(&dir.join(format!("{element}.upf")), Some(&manifest)).unwrap();
        assert_eq!(report["header_facts"]["element"], element);
        assert_eq!(report["qualification_status"], "unqualified");
        assert!(report["header_facts"]["number_of_proj"].as_u64().unwrap() > 3);
        assert_eq!(report["supported_by_engine"], true);
        assert_eq!(report["blockers"], serde_json::json!([]));
        println!(
            "{element}: inspection_s={:.6}, inspection blockers: {}",
            start.elapsed().as_secs_f64(),
            report["blockers"]
        );
    }
}

fn inspection_report(
    path: &std::path::Path,
    manifest: Option<&std::path::Path>,
) -> Result<serde_json::Value, String> {
    let mut command = Command::new(env!("CARGO_BIN_EXE_dft-cli"));
    command.arg("inspect-upf").arg(path).arg("--json");
    if let Some(manifest) = manifest {
        command.arg("--manifest").arg(manifest);
    }
    let output = command.output().unwrap();
    if !output.status.success() {
        return Err(String::from_utf8_lossy(&output.stderr).into_owned());
    }
    serde_json::from_slice(&output.stdout).map_err(|e| e.to_string())
}
