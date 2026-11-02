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
    assert!(
        String::from_utf8_lossy(&result.stderr).contains("unsupported feature: core correction")
    );
    std::fs::write(&path, "<UPF").unwrap();
    let result = run();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("malformed UPF"));
    std::fs::remove_file(&path).unwrap();
    assert!(!run().status.success());
    std::fs::remove_dir(dir).unwrap();
}
