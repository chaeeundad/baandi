use dft_units::*;
#[test]
fn codata_2018_and_definitions() {
    // NIST CODATA 2018: https://physics.nist.gov/cuu/pdf/all_2018.pdf
    // Absolute tolerance is one last printed decimal digit of the MVP constants.
    assert!((HARTREE_IN_EV - 27.211_386_245_988).abs() <= 1e-12);
    assert!((BOHR_IN_ANGSTROM - 0.529_177_210_903).abs() <= 1e-12);
    assert_eq!(from_hartree(1., EnergyUnit::Rydberg), 2.);
    assert_eq!(to_hartree(30., EnergyUnit::Rydberg), 15.);
    let ry_ev = convert_energy(1., EnergyUnit::Rydberg, EnergyUnit::ElectronVolt);
    assert!((ry_ev - 13.605_693_122_994).abs() <= 1e-12);
    let pressure = convert_pressure(
        1.,
        PressureUnit::HartreePerBohrCubed,
        PressureUnit::Gigapascal,
    );
    let independent = 4.359_744_722_207_1e-18 / (5.291_772_109_03e-11_f64).powi(3) / 1e9;
    assert!((pressure / independent - 1.).abs() < 1e-9);
    assert_eq!(
        convert_pressure(1., PressureUnit::Gigapascal, PressureUnit::Kilobar),
        10.
    );
    println!(
        "1 Ha = {HARTREE_IN_EV:.15} eV; 1 Bohr = {BOHR_IN_ANGSTROM:.15} Angstrom; 1 Ry = {ry_ev:.15} eV; 1 Ha/Bohr^3 = {pressure:.15} GPa"
    );
    assert_eq!(Angstrom::from(Bohr(1.)).0, BOHR_IN_ANGSTROM);
    assert_eq!(ElectronVolt::from(Hartree(1.)).0, HARTREE_IN_EV);
    assert_eq!(Hartree::from(Rydberg(30.)).0, 15.);
}
#[test]
fn all_pairs_round_trip() {
    fn check<U: Copy>(units: &[U], convert: impl Fn(f64, U, U) -> f64) {
        let mut worst = 0_f64;
        for &a in units {
            for &b in units {
                for x in [-1e6, -3.7, 1e-9, 0.1, 1., 1e6] {
                    let y = convert(convert(x, a, b), b, a);
                    worst = worst.max(((y - x) / x).abs());
                    assert_eq!(convert(0., a, b), 0.);
                }
            }
        }
        println!("all-pair worst relative round-trip error = {worst:e}");
        assert!(worst < 1e-14);
    }
    check(EnergyUnit::ALL, convert_energy);
    check(LengthUnit::ALL, convert_length);
    check(PressureUnit::ALL, convert_pressure);
    for x in [-1e6, -3.7, 1e-9, 0.1, 1., 1e6] {
        for y in [
            Bohr::from(Angstrom::from(Bohr(x))).0,
            Angstrom::from(Bohr::from(Angstrom(x))).0,
            Hartree::from(ElectronVolt::from(Hartree(x))).0,
            ElectronVolt::from(Hartree::from(ElectronVolt(x))).0,
            Hartree::from(Rydberg::from(Hartree(x))).0,
            Rydberg::from(Hartree::from(Rydberg(x))).0,
        ] {
            assert!(((y - x) / x).abs() < 1e-15); // Preserve the original wrapper tolerance.
        }
    }
}
#[test]
fn strict_parser() {
    for text in ["30 Ry", "30Ry", "30 rY", " 3e1 rydberg "] {
        assert_eq!(parse_energy(text).unwrap(), 15.);
    }
    assert!(parse_energy("408.2 eV").unwrap().is_finite());
    assert_eq!(
        parse_quantity::<EnergyUnit>("1eV").unwrap(),
        (1., EnergyUnit::ElectronVolt)
    );
    assert!(parse_length("5.43 Angstrom").unwrap().is_finite());
    for text in [
        "",
        "30",
        "nan Ha",
        "inf Ha",
        "-infRy",
        "30 unknown",
        "30 Ry extra",
        "x Ry",
        "1e Ry",
        "1e309Ha",
    ] {
        let err = parse_energy(text).unwrap_err();
        assert_eq!(err.text, text);
    }
    assert_eq!(
        parse_energy("NaN Ha").unwrap_err().kind,
        UnitErrorKind::NonFinite
    );
    assert_eq!(
        parse_energy("30").unwrap_err().kind,
        UnitErrorKind::MissingUnit
    );
    assert_eq!(
        parse_energy("30 Ha Ry").unwrap_err().kind,
        UnitErrorKind::ExtraTokens
    );
    for alias in ["bohr", "a0", "angstrom", "ang", "A", "nm", "pm", "m"] {
        assert!(parse_length(&format!("1{alias}")).is_ok());
    }
    assert!(parse_length("1 Ha").is_err());
    assert!(parse_energy("1e308 J").is_err());
}
