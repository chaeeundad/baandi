#[allow(dead_code)]
#[path = "../common/mod.rs"]
mod common;
use dft_upf::{parse_upf, radial::RadialEvenRule};

fn fixture(metadata: &str, prefix: &str, tail: &str) -> String {
    let mut text = common::fixture_with_points(7);
    let start = text.find("<PP_BETA.1 ").unwrap();
    let end = text.find("</PP_BETA.1>").unwrap() + "</PP_BETA.1>".len();
    text.replace_range(
        start..end,
        &format!(
            "<PP_BETA.1 index=\"1\" angular_momentum=\"0\" {metadata}>{prefix} {tail}</PP_BETA.1>"
        ),
    );
    text
}

const PREFIX: &str = "0 1.000000000000e0 2.000000000000e0 1.000000000000e0";

#[test]
fn full_mesh_roundoff_tails_have_identical_supported_integrals() {
    for cutoff in [4, 5] {
        let prefix = if cutoff == 4 {
            PREFIX.to_string()
        } else {
            format!("{PREFIX} 0")
        };
        let zeros = vec!["0"; 7 - cutoff].join(" ");
        let tails = vec!["-1.000000000000e-12"; 7 - cutoff].join(" ");
        let indexed = format!("cutoff_radius_index=\"{cutoff}\"");
        let r = common::mesh(7, 12.0).0;
        for metadata in [
            indexed.clone(),
            // Index wins over a conflicting informational radius.
            format!("{indexed} cutoff_radius=\"0\""),
            format!("cutoff_radius=\"{}\"", r[cutoff - 1]),
            format!("cutoff_radius=\"{}\"", (r[cutoff - 1] + r[cutoff]) / 2.0),
        ] {
            let clean = parse_upf(&fixture(&metadata, &prefix, &zeros)).unwrap();
            let tailed = parse_upf(&fixture(&metadata, &prefix, &tails)).unwrap();
            let short = parse_upf(&fixture(&metadata, &prefix, "")).unwrap();
            match &tailed.raw().nonlocal.betas[0] {
                upf::model::PpBetaNode::Beta1(beta) => {
                    assert_eq!(beta.values.len(), 7);
                    assert_eq!(beta.values[cutoff].to_bits(), (-1e-12_f64).to_bits());
                }
                _ => panic!("first raw projector has wrong tag"),
            }
            for pp in [&clean, &tailed, &short] {
                assert_eq!(pp.projectors()[0].cutoff_index(), cutoff);
                assert_eq!(pp.projectors()[0].radial_times_r().len(), cutoff);
            }
            for rule in [
                RadialEvenRule::Simpson3_8,
                RadialEvenRule::SimpsonQuadraticEnd,
            ] {
                for q in [0.0, 0.1, 1.0, 8.0] {
                    // q is in Bohr^-1; the radial integral is in Bohr^3/2.
                    let expected = clean.projector_transform_with_rule(0, q, rule).unwrap();
                    for pp in [&tailed, &short] {
                        let actual = pp.projector_transform_with_rule(0, q, rule).unwrap();
                        assert!((actual - expected).abs() <= 1e-12);
                        assert_eq!(actual.to_bits(), expected.to_bits());
                    }
                }
            }
        }
    }
}

#[test]
fn real_tail_diagnostic_and_nonfinite_samples() {
    let metadata = "cutoff_radius_index=\"4\"";
    let error = parse_upf(&fixture(metadata, PREFIX, "0 -0.01 0.001"))
        .unwrap_err()
        .to_string();
    assert_eq!(
        error,
        format!(
            "invalid field PP_BETA.1: cutoff index 4: largest tail magnitude {:.17e} Bohr^-1/2 at index 6; relative tolerance {:.17e}, projector maximum {:.17e} Bohr^-1/2",
            0.01, 1e-12, 2.0
        )
    );
    for value in ["NaN", "inf", "-inf"] {
        let error = parse_upf(&fixture(metadata, PREFIX, &format!("0 {value} 0")))
            .unwrap_err()
            .to_string();
        assert_eq!(
            error,
            "invalid field PP_BETA.1: non-finite sample at index 6 (cutoff index 4)"
        );
        assert!(parse_upf(&fixture(metadata, &format!("0 {value} 2 1"), "0 0 0")).is_err());
    }
}

#[test]
fn precision_scaling_zero_support_and_metadata_boundaries() {
    let metadata = "cutoff_radius_index=\"4\"";
    // Same relative tail budget under amplitude and sign changes.
    for exponent in [-80, 0, 80] {
        let prefix = format!("0 -1.000000000000e{exponent} 2.000000000000e{exponent} 0");
        assert!(
            parse_upf(&fixture(
                metadata,
                &prefix,
                &format!("0 1e{} 0", exponent - 12)
            ))
            .is_ok()
        );
        assert!(
            parse_upf(&fixture(
                metadata,
                &prefix,
                &format!("0 3e{} 0", exponent - 12)
            ))
            .is_err()
        );
    }
    let precise = PREFIX.replace(".000000000000", ".0000000000000000");
    assert!(parse_upf(&fixture(metadata, &precise, "0 1e-12 0")).is_err());
    assert!(parse_upf(&fixture(metadata, "0 1 2 1", "0 1e-12 0")).is_err());
    // An underflowed relative ratio cannot create a budget for plain integers.
    assert!(parse_upf(&fixture(metadata, "0 100 200 100", "0 5e-324 0")).is_err());
    assert!(parse_upf(&fixture(metadata, "0 1.0 2.0 1.0", "0 0.01 0")).is_err());
    assert!(parse_upf(&fixture(metadata, "0 0 0 0", "0 1e-99 0")).is_err());
    for metadata in [
        "",
        "cutoff_radius_index=\"0\"",
        "cutoff_radius_index=\"8\"",
        "cutoff_radius=\"NaN\"",
        "cutoff_radius=\"-1\"",
    ] {
        assert!(parse_upf(&fixture(metadata, PREFIX, "0 0 0")).is_err());
    }
}
