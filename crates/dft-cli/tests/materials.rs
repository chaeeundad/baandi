//! Real UPFs are external and never copied. Run only in release mode.
#[test]
#[ignore = "requires DFT_NC_DIR real PseudoDojo files; release bulk smoke"]
fn real_nc_bulk_insulators() {
    assert!(std::env::var_os("DFT_NC_DIR").is_some(), "set DFT_NC_DIR");
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let status = std::process::Command::new("python3")
        .arg("-I")
        .arg(root.join("scripts/run_materials_smoke.py"))
        .arg("--binary")
        .arg(env!("CARGO_BIN_EXE_dft-cli"))
        .current_dir(&root)
        .status()
        .unwrap();
    assert!(status.success(), "bulk smoke exit {status}");
    let explicit: serde_json::Value = serde_json::from_slice(
        &std::fs::read(root.join("tests/reference/materials/results.json")).unwrap(),
    )
    .unwrap();
    let status = std::process::Command::new("python3")
        .arg("-I")
        .arg(root.join("scripts/run_materials_smoke.py"))
        .arg("--binary")
        .arg(env!("CARGO_BIN_EXE_dft-cli"))
        .arg("--implicit-spin")
        .arg("--default-floor")
        .arg("--output")
        .arg(root.join("tests/reference/materials/implicit-results.json"))
        .current_dir(&root)
        .status()
        .unwrap();
    assert!(status.success());
    let implicit: serde_json::Value = serde_json::from_slice(
        &std::fs::read(root.join("tests/reference/materials/implicit-results.json")).unwrap(),
    )
    .unwrap();
    for (a, b) in explicit["cases"]
        .as_array()
        .unwrap()
        .iter()
        .zip(implicit["cases"].as_array().unwrap())
    {
        assert_eq!(a["density_sha256"], b["density_sha256"]);
        assert_eq!(a["total_energy_ha"], b["total_energy_ha"]);
        assert_eq!(a["iterations"], b["iterations"]);
    }
}
