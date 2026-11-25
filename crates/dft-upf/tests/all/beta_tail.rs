#[allow(dead_code)]
#[path = "../common/mod.rs"]
mod common;
use dft_upf::{parse_upf, radial::RadialEvenRule};
use std::f64::consts::PI;

fn fixture(n: usize, metadata: &str, samples: &str) -> String {
    let mut text = common::fixture_with_points(n);
    let start = text.find("<PP_BETA.1 ").unwrap();
    let end = text.find("</PP_BETA.1>").unwrap() + "</PP_BETA.1>".len();
    text.replace_range(
        start..end,
        &format!("<PP_BETA.1 index=\"1\" angular_momentum=\"0\" {metadata}>{samples}</PP_BETA.1>"),
    );
    text
}

const RULES: [RadialEvenRule; 2] = [
    RadialEvenRule::Simpson3_8,
    RadialEvenRule::SimpsonQuadraticEnd,
];

// Independent full-mesh weight-vector quadrature, with j0 from sin(x)/x.
// Does not call production integration, Bessel, or projector transforms.
fn direct_full_mesh(r: &[f64], rab: &[f64], u: &[f64], q: f64, rule: RadialEvenRule) -> f64 {
    let n = r.len();
    let mut weights = vec![0.0; n];
    let end = if !n.is_multiple_of(2) {
        n - 1
    } else if rule == RadialEvenRule::Simpson3_8 {
        n - 4
    } else {
        n - 2
    };
    for i in (0..end).step_by(2) {
        for (j, w) in [1.0 / 3.0, 4.0 / 3.0, 1.0 / 3.0].into_iter().enumerate() {
            weights[i + j] += w;
        }
    }
    if n.is_multiple_of(2) {
        let panel: &[f64] = if rule == RadialEvenRule::Simpson3_8 {
            &[3.0 / 8.0, 9.0 / 8.0, 9.0 / 8.0, 3.0 / 8.0]
        } else {
            &[-1.0 / 12.0, 8.0 / 12.0, 5.0 / 12.0]
        };
        for (j, w) in panel.iter().enumerate() {
            weights[n - panel.len() + j] += w;
        }
    }
    4.0 * PI
        * (0..n)
            .map(|i| {
                let x = q * r[i];
                let j0 = if x == 0.0 { 1.0 } else { x.sin() / x };
                weights[i] * r[i] * u[i] * j0 * rab[i]
            })
            .sum::<f64>()
}

#[test]
fn oncv_and_hgh_tails_match_direct_full_mesh_quadrature() {
    let mut maximum_error = 0.0_f64;
    for n in [367, 368, 898, 899] {
        let (r, rab) = common::mesh(n, 12.0);
        let oncv = n < 400;
        let cutoff = if oncv { n - 6 } else { n - 125 };
        // Smooth sign-changing continuation, or an analytic Gaussian tail.
        // These are synthetic functions, not third-party pseudopotential data.
        let u: Vec<_> = r
            .iter()
            .map(|&x| {
                if oncv {
                    0.01 * x * (r[n - 4] - x) * (-0.1 * x * x).exp()
                } else {
                    0.1 * x * (-x * x).exp()
                }
            })
            .collect();
        assert!(u[cutoff..].iter().any(|&x| x != 0.0));
        if oncv {
            assert!(u[cutoff] > 0.0 && u[n - 1] < 0.0);
        }
        let pp = parse_upf(&fixture(
            n,
            &format!("cutoff_radius_index=\"{cutoff}\""),
            &common::array(&u),
        ))
        .unwrap();
        let p = &pp.projectors()[0];
        assert_eq!(p.declared_cutoff_index(), cutoff);
        assert_eq!(p.effective_support_index(), n);
        assert_eq!(p.cutoff_index(), n);
        assert_eq!(
            p.largest_tail_magnitude(),
            u[cutoff..].iter().map(|x| x.abs()).fold(0.0, f64::max)
        );
        assert_eq!(p.radial_times_r(), u);
        for rule in RULES {
            for q in [0.0, 0.1, 1.0, 8.0] {
                // q: Bohr^-1; transform and error: Bohr^3/2.
                let expected = direct_full_mesh(&r, &rab, &u, q, rule);
                let actual = pp.projector_transform_with_rule(0, q, rule).unwrap();
                let error = (actual - expected).abs();
                maximum_error = maximum_error.max(error);
                assert!(
                    error <= 1e-14,
                    "n={n} q={q} rule={rule:?}: {error:e} Bohr^3/2"
                );
            }
        }
    }
    println!("maximum full-mesh quadrature error: {maximum_error:.17e} Bohr^3/2");
}

#[test]
fn zero_tails_and_shortened_files_remain_bit_identical() {
    let r = common::mesh(7, 12.0).0;
    for cutoff in [4, 5] {
        let mut u = vec![0.0, 1.0, 2.0, 0.0];
        u.resize(cutoff, 0.0);
        let indexed = format!("cutoff_radius_index=\"{cutoff}\"");
        for metadata in [
            indexed.clone(),
            format!("{indexed} cutoff_radius=\"0\""),
            format!("cutoff_radius=\"{}\"", r[cutoff - 1]),
            format!("cutoff_radius=\"{}\"", (r[cutoff - 1] + r[cutoff]) / 2.0),
        ] {
            let short = parse_upf(&fixture(7, &metadata, &common::array(&u))).unwrap();
            let mut padded = u.clone();
            padded.resize(7, -0.0);
            let full = parse_upf(&fixture(7, &metadata, &common::array(&padded))).unwrap();
            for pp in [&short, &full] {
                let p = &pp.projectors()[0];
                assert_eq!(p.declared_cutoff_index(), cutoff);
                assert_eq!(p.effective_support_index(), cutoff);
                assert_eq!(p.largest_tail_magnitude(), 0.0);
                assert_eq!(p.radial_times_r(), u);
            }
            for rule in RULES {
                for q in [0.0, 0.1, 1.0, 8.0] {
                    assert_eq!(
                        short
                            .projector_transform_with_rule(0, q, rule)
                            .unwrap()
                            .to_bits(),
                        full.projector_transform_with_rule(0, q, rule)
                            .unwrap()
                            .to_bits()
                    );
                }
            }
        }
    }
}

#[test]
fn finite_tails_have_no_threshold_and_only_exact_zero_padding_is_removed() {
    for tail in ["5e-324", "-1e-99", "1e-12", "-0.01", "1e80"] {
        for prefix in ["0 0 0 0", "0 1 2 1", "0 1.0000000000000000 2.0 1.0"] {
            let pp = parse_upf(&fixture(
                7,
                "cutoff_radius_index=\"4\"",
                &format!("{prefix} 0 {tail} -0"),
            ))
            .unwrap();
            let p = &pp.projectors()[0];
            assert_eq!(p.effective_support_index(), 6);
            assert_eq!(p.radial_times_r().len(), 6);
            assert_eq!(
                p.radial_times_r()[5].to_bits(),
                tail.parse::<f64>().unwrap().to_bits()
            );
            assert_eq!(
                p.largest_tail_magnitude(),
                tail.parse::<f64>().unwrap().abs()
            );
            match &pp.raw().nonlocal.betas[0] {
                upf::model::PpBetaNode::Beta1(b) => assert_eq!(b.values.len(), 7),
                _ => panic!("wrong raw projector tag"),
            }
        }
    }
}

#[test]
fn nonfinite_samples_and_invalid_metadata_reject() {
    for value in ["NaN", "inf", "-inf"] {
        for i in 0..7 {
            let mut u = ["0"; 7];
            u[i] = value;
            let error = parse_upf(&fixture(7, "cutoff_radius_index=\"4\"", &u.join(" ")))
                .unwrap_err()
                .to_string();
            assert_eq!(
                error,
                format!(
                    "invalid field PP_BETA.1: non-finite sample at index {} (cutoff index 4)",
                    i + 1
                )
            );
        }
    }
    for metadata in [
        "",
        "cutoff_radius_index=\"0\"",
        "cutoff_radius_index=\"8\"",
        "cutoff_radius=\"NaN\"",
        "cutoff_radius=\"-1\"",
    ] {
        assert!(parse_upf(&fixture(7, metadata, "0 1 2 1 0 0 0")).is_err());
    }
}
