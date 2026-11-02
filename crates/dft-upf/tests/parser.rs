mod common;
use dft_upf::{UpfError, parse_upf};
fn replace_body(text: &str, tag: &str, value: &str) -> String {
    let start = text.find(&format!("<{tag}>")).unwrap() + tag.len() + 2;
    let end = text.find(&format!("</{tag}>")).unwrap();
    format!("{}{value}{}", &text[..start], &text[end..])
}
#[test]
fn synthetic_round_trip() {
    let p = parse_upf(&common::fixture()).unwrap();
    assert_eq!(p.element(), "Xx");
    assert_eq!(p.z_valence(), 4.0);
    assert_eq!(p.functional(), "SLA PZ");
    let (r, rab) = common::mesh(601, 12.0);
    assert_eq!(p.mesh().r(), r);
    assert_eq!(p.mesh().rab(), rab);
    assert_eq!(
        p.dij_hartree(),
        &[1.0, 0.2, 0.0, 0.2, 2.0, 0.0, 0.0, 0.0, 3.0]
    );
    for (i, b) in p.projectors().iter().enumerate() {
        assert_eq!(b.l(), if i == 2 { 1 } else { 0 });
        assert_eq!(b.cutoff_index(), 601);
        for (j, r) in r.iter().enumerate() {
            assert_eq!(
                b.radial_times_r()[j],
                r.powi(b.l() as i32 + 1) * (-r * r / (i + 1) as f64).exp()
            );
        }
    }
    assert_eq!(p.local_hartree()[0], -8.0 / std::f64::consts::PI.sqrt());
    for (i, r) in r.iter().enumerate().skip(1) {
        assert_eq!(p.local_hartree()[i], -4.0 * libm::erf(*r) / r);
    }
}
#[test]
fn negative_variants() {
    let s = common::fixture();
    for (old, new, feature) in [
        ("is_paw=\"F\"", "is_paw=\"T\"", "PAW"),
        ("is_ultrasoft=\"F\"", "is_ultrasoft=\"T\"", "ultrasoft"),
        ("has_so=\"F\"", "has_so=\"T\"", "spin-orbit"),
        (
            "core_correction=\"F\"",
            "core_correction=\"T\"",
            "core correction",
        ),
        ("pseudo_type=\"NC\"", "pseudo_type=\"SL\"", "pseudo_type=SL"),
    ] {
        assert!(
            matches!(parse_upf(&s.replace(old,new)),Err(UpfError::Unsupported(f)) if f==feature)
        );
    }
    assert!(matches!(
        parse_upf(&s.replace("version=\"2.0.1\"", "version=\"1.0\"")),
        Err(UpfError::Version(_))
    ));
    assert!(matches!(
        parse_upf(&s[..s.len() - 10]),
        Err(UpfError::Parse(_))
    ));
    for tag in ["PP_LOCAL", "PP_BETA.1"] {
        let start = s.find(&format!("<{tag}")).unwrap();
        let end = s.find(&format!("</{tag}>")).unwrap() + tag.len() + 3;
        let removed = format!("{}{}", &s[..start], &s[end..]);
        assert!(matches!(parse_upf(&removed),Err(UpfError::Missing(f)) if f==tag));
    }
    let (mut r, _) = common::mesh(601, 12.0);
    r[20] = r[19];
    for (text, field) in [
        (replace_body(&s, "PP_R", &common::array(&r)), "mesh"),
        (replace_body(&s, "PP_RAB", "1 2"), "PP_RAB"),
        (
            replace_body(&s, "PP_LOCAL", &vec!["NaN"; 601].join(" ")),
            "PP_LOCAL",
        ),
        (replace_body(&s, "PP_DIJ", "1 2"), "PP_DIJ"),
        (
            s.replace("z_valence=\"4\"", "z_valence=\"NaN\""),
            "z_valence",
        ),
        (
            s.replace("cutoff_radius_index=\"601\"", "cutoff_radius_index=\"602\""),
            "cutoff_radius_index",
        ),
    ] {
        let result = parse_upf(&text);
        assert!(
            matches!(result,Err(UpfError::Invalid{field:f,..}) if f==field),
            "expected {field}"
        );
    }
    assert!(matches!(
        parse_upf(&replace_body(&s, "PP_DIJ", "2 0.5 0 0.4 4 0 0 0 6")),
        Err(UpfError::NonSymmetricD)
    ));
    assert!(
        matches!(parse_upf(&s.replace("</UPF>","<PP_SPIN_ORB/></UPF>")),Err(UpfError::Unsupported(f)) if f=="PP_SPIN_ORB")
    );
    // Comments containing fake metadata cannot change the support verdict.
    assert!(
        parse_upf(&s.replace("<PP_HEADER", "<!-- <PP_HEADER is_paw='T'/> --><PP_HEADER")).is_ok()
    );
}

#[test]
fn supplemental_metadata_and_storage_boundaries() {
    let s = common::fixture();
    let p = parse_upf(&s.replace("SLA PZ", "SLA&#32;PZ")).unwrap();
    assert_eq!(p.functional(), "SLA PZ");
    assert_eq!(p.raw().local.values[0], 2.0 * p.local_hartree()[0]);
    assert_eq!(p.raw().nonlocal.dij.values[0], 2.0 * p.dij_hartree()[0]);
    for (key, feature) in [
        ("core_correction", "core correction"),
        ("has_so", "spin-orbit"),
        ("has_gipaw", "GIPAW"),
        ("has_wfc", "all-electron wavefunctions"),
    ] {
        let text = s.replace(&format!("{key}=\"F\""), &format!("{key}=\".T.\""));
        assert!(matches!(parse_upf(&text), Err(UpfError::Unsupported(f)) if f == feature));
    }
    for (old, new, field) in [
        ("<PP_LOCAL>", "<PP_LOCAL size=\"600\">", "PP_LOCAL"),
        ("<PP_DIJ>", "<PP_DIJ size=\"8\">", "PP_DIJ"),
        (
            "cutoff_radius_index=\"601\"",
            "cutoff_radius_index=\"600\"",
            "PP_BETA",
        ),
        ("index=\"1\"", "index=\"2\"", "PP_BETA/index"),
        (
            "number_of_proj=\"3\"",
            "number_of_proj=\"2\"",
            "number_of_proj",
        ),
    ] {
        assert!(
            matches!(parse_upf(&s.replacen(old, new, 1)), Err(UpfError::Invalid{field:f,..}) if f == field)
        );
    }
    assert!(
        matches!(parse_upf(&s.replace("<PP_LOCAL>", "<PP_LOCAL units=\"Ha\">")), Err(UpfError::Unsupported(f)) if f.contains("units=Ha"))
    );
    assert!(
        matches!(parse_upf(&s.replace("number_of_proj=\"3\"", "number_of_proj=\"4\"")), Err(UpfError::Unsupported(f)) if f.contains("three projectors"))
    );
    assert!(
        matches!(parse_upf(&s.replace("angular_momentum=\"1\"", "angular_momentum=\"4\"")), Err(UpfError::Unsupported(f)) if f == "projector l=4")
    );
    assert!(
        matches!(parse_upf(&s.replace("<PP_DIJ>", "<PP_AUGMENTATION/><PP_DIJ>")), Err(UpfError::Unsupported(f)) if f == "PP_AUGMENTATION")
    );
    assert!(matches!(
        parse_upf(&format!("{s}{s}")),
        Err(UpfError::Parse(_))
    ));
    let (r, mut rab) = common::mesh(601, 12.0);
    rab[8] = 0.0;
    assert!(
        matches!(parse_upf(&replace_body(&s,"PP_RAB",&common::array(&rab))),Err(UpfError::Invalid{field,..}) if field=="mesh")
    );
    // A genuinely shortened projector is accepted and transformed only through its cutoff.
    let start = s.find("<PP_BETA.1").unwrap();
    let body = s[start..].find('>').unwrap() + start + 1;
    let end = s.find("</PP_BETA.1>").unwrap();
    let u: Vec<_> = r[..401].iter().map(|r| r * (-r * r).exp()).collect();
    let shorter = format!(
        "{}{}{}",
        s[..body].replacen(
            "cutoff_radius_index=\"601\"",
            "cutoff_radius_index=\"401\"",
            1
        ),
        common::array(&u),
        &s[end..]
    );
    let p = parse_upf(&shorter).unwrap();
    assert_eq!(p.projectors()[0].cutoff_index(), 401);
    assert_eq!(p.projectors()[0].radial_times_r(), u);
    assert!(p.projector_transform(0, 1.0).unwrap().is_finite());
}

#[test]
fn real_world_layouts_are_bit_exact() {
    use common::{Options, fixture_with_options};
    // Zero tolerance in the stored UPF units and converted atomic units.
    let bits = |xs: &[f64]| xs.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    let baseline = Options {
        wavefunctions: true,
        ..Options::default()
    };
    let expected = parse_upf(&fixture_with_options(601, &baseline)).unwrap();
    for (relativistic, wrap, separator, crlf) in [
        ("no", false, " ", false),
        ("nonrelativistic", false, " ", false),
        ("scalar", true, " ", false),
        ("scalar", true, " ", true),
        ("scalar", false, "\t", false),
        ("scalar", false, "  ", false),
        ("scalar", false, "\r", false),
        ("scalar", false, "\x0c", false),
        ("scalar", false, "\x0b", false),
        ("no", true, " \t  \r\x0c", true),
    ] {
        let options = Options {
            relativistic,
            wrap,
            separator,
            crlf,
            wavefunctions: true,
        };
        let text = fixture_with_options(601, &options);
        if wrap {
            // Both the mesh (601) and D (9) end with a one-value row.
            for tag in ["PP_R", "PP_DIJ"] {
                let start = text.find(&format!("<{tag} columns=\"4\">")).unwrap();
                let end = text.find(&format!("</{tag}>")).unwrap();
                assert_eq!(
                    text[start..end]
                        .lines()
                        .last()
                        .unwrap()
                        .split_ascii_whitespace()
                        .count(),
                    1
                );
            }
        }
        let actual = parse_upf(&text).unwrap();
        for (a, b) in [
            (actual.mesh().r(), expected.mesh().r()),
            (actual.mesh().rab(), expected.mesh().rab()),
            (actual.local_hartree(), expected.local_hartree()),
            (actual.dij_hartree(), expected.dij_hartree()),
            (&actual.raw().local.values, &expected.raw().local.values),
            (
                &actual.raw().nonlocal.dij.values,
                &expected.raw().nonlocal.dij.values,
            ),
            (&actual.raw().rhoatom.values, &expected.raw().rhoatom.values),
        ] {
            assert_eq!(bits(a), bits(b));
        }
        for (a, b) in actual.projectors().iter().zip(expected.projectors()) {
            assert_eq!(bits(a.radial_times_r()), bits(b.radial_times_r()));
        }
        use upf::model::PpChiNode;
        let orbital_bits = |p: &dft_upf::Pseudopotential| {
            p.raw()
                .pswfc
                .as_ref()
                .unwrap()
                .orbitals
                .iter()
                .map(|node| {
                    let (PpChiNode::Chi1(w) | PpChiNode::Chi2(w) | PpChiNode::Chi3(w)) = node;
                    bits(&w.values)
                })
                .collect::<Vec<_>>()
        };
        assert_eq!(orbital_bits(&actual), orbital_bits(&expected));
        assert_eq!(actual.z_valence().to_bits(), expected.z_valence().to_bits());
        // Also compare every raw header/mesh/projector/orbital attribute. Only
        // the requested relativistic spelling differs from the baseline.
        let mut raw = actual.raw().clone();
        assert_eq!(raw.header.relativistic, relativistic);
        raw.header.relativistic = expected.raw().header.relativistic.clone();
        assert_eq!(&raw, expected.raw());
    }
}

#[test]
fn unsupported_relativistic_values_are_named() {
    for value in ["full", "foo"] {
        let text = common::fixture_with_options(
            601,
            &common::Options {
                relativistic: value,
                ..common::Options::default()
            },
        );
        let error = parse_upf(&text).unwrap_err();
        assert_eq!(
            error.to_string(),
            format!("unsupported feature: relativistic={value}")
        );
        assert!(matches!(error, UpfError::Unsupported(f) if f == format!("relativistic={value}")));
    }
}

#[test]
fn inspect_wrapped_nonrelativistic_synthetic_file() {
    // This test lives in dft-upf to respect ticket ownership. Cargo builds the
    // CLI explicitly because CARGO_BIN_EXE is only set for a package's own bins.
    let path = std::env::temp_dir().join(format!("dft-upf-wrapped-{}.xml", std::process::id()));
    let text = common::fixture_with_options(
        601,
        &common::Options {
            relativistic: "no",
            wrap: true,
            crlf: true,
            ..common::Options::default()
        },
    );
    std::fs::write(&path, text).unwrap();
    let result = std::process::Command::new(env!("CARGO"))
        .args([
            "run",
            "--offline",
            "--quiet",
            "-p",
            "dft-cli",
            "--",
            "inspect-upf",
        ])
        .arg(&path)
        .current_dir(concat!(env!("CARGO_MANIFEST_DIR"), "/../.."))
        .output()
        .unwrap();
    std::fs::remove_file(path).unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let stdout = String::from_utf8(result.stdout).unwrap();
    assert!(stdout.contains("Support: SUPPORTED"), "{stdout}");
    assert!(stdout.contains("601 points"), "{stdout}");
    assert!(stdout.contains("z_valence: 4"), "{stdout}");
}

#[test]
fn functional_label_preserves_whitespace() {
    let label = " SLA  PZ   NOGX NOGC ";
    let p = parse_upf(&common::fixture().replace("SLA PZ", label)).unwrap();
    assert_eq!(p.functional(), label);
    assert!(parse_upf(&common::fixture().replace("SLA PZ", "   ")).is_err());
}
