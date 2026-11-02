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
