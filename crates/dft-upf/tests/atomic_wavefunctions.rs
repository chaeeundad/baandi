mod common;
use dft_upf::parse_upf;

#[test]
fn optional_atomic_wavefunctions_preserve_metadata_and_samples() {
    assert!(
        parse_upf(&common::fixture())
            .unwrap()
            .atomic_wavefunctions()
            .is_none()
    );
    let text = common::fixture_with_options(
        31,
        &common::Options {
            wavefunctions: true,
            ..Default::default()
        },
    )
    .replace("occupation=\"1\"", "occupation=\" 0.0 \"");
    let pp = parse_upf(&text).unwrap();
    let orbitals = pp.atomic_wavefunctions().unwrap();
    assert_eq!(orbitals.len(), 3);
    for (i, chi) in orbitals.iter().enumerate() {
        assert_eq!(chi.l, 0);
        assert_eq!(chi.occupation, Some(0.0));
        assert_eq!(
            chi.label.as_deref(),
            Some(format!("synthetic{}", i + 1).as_str())
        );
        assert_eq!(chi.radial_times_r, pp.mesh().r());
    }
    for bad in ["NaN", "-1"] {
        assert!(parse_upf(&text.replace(" 0.0 ", bad)).is_err());
    }
    assert!(parse_upf(&text.replace("l=\"0\"", "l=\"-1\"")).is_err());
}
