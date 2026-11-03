use dft_units::*;

fn ulp(text: &str, expected: f64) {
    let actual = eval_expr(text).unwrap();
    let distance = actual.to_bits().abs_diff(expected.to_bits());
    println!("{text} = {actual:.17e}; std = {expected:.17e}; ulps = {distance}");
    assert!(distance <= 1, "{text}: {actual} vs {expected}");
}

#[test]
fn arithmetic_and_precedence() {
    for (text, expected) in [
        ("1/4", 0.25_f64),
        ("2+1/3", 2.0 + 1.0 / 3.0),
        ("2+3*4", 14.0),
        ("2\t+\n3", 5.0),
        ("4+1e-13", 4.0 + 1e-13),
        ("2*1E+3", 2000.0),
        ("(2+3)*4", 20.0),
        ("2^3^2", 64.0),
        ("-2^2", 4.0),
        ("-(2^2)", -4.0),
        ("8/4/2", 1.0),
        ("8-4-2", 2.0),
        (" \t2 + 1 / 3\n", 2.0 + 1.0 / 3.0),
        (" ( 2 + 3 ) * 4 ", 20.0),
        ("pi", std::f64::consts::PI),
        ("e", std::f64::consts::E),
        ("e+pi", std::f64::consts::E + std::f64::consts::PI),
    ] {
        let actual = eval_expr(text).unwrap();
        println!("{text:?} = {actual:.17e} (exact bits)");
        assert_eq!(actual.to_bits(), expected.to_bits(), "{text}");
    }
    ulp("atan(0.5*pi)", (0.5 * std::f64::consts::PI).atan());
    ulp("sqrt(3)/4", 3.0_f64.sqrt() / 4.0);
    for s in ["-0", "+2", "1e-3", " 2.5 ", "1.7976931348623157e308"] {
        assert_eq!(
            eval_expr(s).unwrap().to_bits(),
            s.trim().parse::<f64>().unwrap().to_bits()
        );
    }
    let resolved = resolve_expr(" 1/4 ").unwrap();
    assert_eq!(resolved.text, " 1/4 ");
    assert_eq!(resolved.value, 0.25);
}

#[test]
fn functions_match_std_within_one_ulp() {
    type Function = (&'static str, fn(f64) -> f64);
    let functions: &[Function] = &[
        ("sqrt", f64::sqrt),
        ("abs", f64::abs),
        ("sin", f64::sin),
        ("cos", f64::cos),
        ("tan", f64::tan),
        ("asin", f64::asin),
        ("acos", f64::acos),
        ("atan", f64::atan),
        ("exp", f64::exp),
        ("ln", f64::ln),
        ("log10", f64::log10),
    ];
    for &(name, function) in functions {
        for x in [0.125, 0.5, 0.875] {
            ulp(&format!("{name}({x})"), function(x));
        }
    }
    ulp("abs(-2)", 2.0);
    ulp("sin(-0.5)", (-0.5_f64).sin());
}

#[test]
fn typed_errors_quote_input() {
    for (text, kind) in [
        ("", ExprErrorKind::Empty),
        ("  ", ExprErrorKind::Empty),
        ("2+", ExprErrorKind::Syntax),
        ("sqrt(", ExprErrorKind::Syntax),
        ("mystery(2)", ExprErrorKind::Syntax),
        ("x=2", ExprErrorKind::Syntax),
        ("x", ExprErrorKind::UnknownName),
        ("pi+x", ExprErrorKind::UnknownName),
        ("1/0", ExprErrorKind::NonFinite),
        ("sqrt(-1)", ExprErrorKind::NonFinite),
        ("asin(2)", ExprErrorKind::NonFinite),
        ("NaN", ExprErrorKind::NonFinite),
        ("inf", ExprErrorKind::NonFinite),
        ("exp(1000)", ExprErrorKind::NonFinite),
    ] {
        let err = eval_expr(text).unwrap_err();
        assert_eq!(err.kind, kind, "{err}");
        assert_eq!(err.text, text);
        assert!(err.to_string().contains(&format!("{text:?}")));
    }
    for text in [
        "log(2)", "log2(2)", "floor(2)", "sinh(2)", "PI", "E", "tau", "π", "{pi}", "2%1", "2**3",
        "1,2", "2<3", "2 & 3", "pi(2)",
    ] {
        assert!(eval_expr(text).is_err(), "unexpectedly accepted {text}");
    }
}

#[test]
fn integers() {
    for (s, n) in [
        ("8/2", 4),
        ("-8/2", -4),
        ("0", 0),
        ("4+1e-13", 4),
        ("-9223372036854775808", i64::MIN),
    ] {
        assert_eq!(eval_expr_integer(s).unwrap(), n);
    }
    for s in ["7/2", "4+1e-10", "1e-15"] {
        assert_eq!(
            eval_expr_integer(s).unwrap_err().kind,
            ExprErrorKind::NotInteger
        );
    }
    for s in ["9223372036854775808", "-1e20"] {
        assert_eq!(
            eval_expr_integer(s).unwrap_err().kind,
            ExprErrorKind::IntegerOutOfRange
        );
    }
}

#[test]
fn quantities_and_legacy() {
    let (v, u) = parse_quantity_expr::<LengthUnit>("5.43/2 Angstrom").unwrap();
    assert_eq!(v, 2.715);
    assert_eq!(u, LengthUnit::Angstrom);
    println!("5.43/2 Angstrom = {v:.17e} Angstrom");
    let (v, u) = parse_quantity_expr::<EnergyUnit>("0.5*30 Ry").unwrap();
    assert_eq!(v, 15.0);
    assert_eq!(to_hartree(v, u), 7.5);
    println!("0.5*30 Ry = {v} Ry = {} Ha", to_hartree(v, u));
    assert_eq!(
        parse_quantity_expr::<LengthUnit>("2 * 5.43 Angstrom")
            .unwrap()
            .0,
        10.86
    );
    assert_eq!(
        parse_quantity_expr::<LengthUnit>("(2 + 3) / 2\t a0 ")
            .unwrap()
            .0,
        2.5
    );
    assert_eq!(
        parse_quantity_expr::<PressureUnit>("2/4 Ha/Bohr^3")
            .unwrap()
            .0,
        0.5
    );
    for s in ["30 Ry", "30Ry", "1eV", "-0 Ha", "1e-3 eV", "  30  Ry  "] {
        let old = parse_quantity::<EnergyUnit>(s).unwrap();
        let new = parse_quantity_expr::<EnergyUnit>(s).unwrap();
        assert_eq!(old.0.to_bits(), new.0.to_bits());
        assert_eq!(old.1, new.1);
    }
    for s in ["", "30", "NaN Ry", "30 nope"] {
        assert_eq!(
            parse_quantity_expr::<EnergyUnit>(s).unwrap_err().kind,
            ExprErrorKind::Unit(parse_quantity::<EnergyUnit>(s).unwrap_err().kind)
        );
    }
    for s in [
        "30 Ry extra",
        "1/0 Ry",
        "2+ Ry",
        "2 * Ry",
        "2 Ry + 3 Ry",
        "2/4Ry",
    ] {
        let e = parse_quantity_expr::<EnergyUnit>(s).unwrap_err();
        assert_eq!(e.text, s);
    }
}

#[test]
fn no_typed_math_constants_in_sources() {
    // Derive forbidden digit sequences from std, never duplicate the constants.
    let prefixes: Vec<String> = [std::f64::consts::PI, std::f64::consts::E]
        .iter()
        .map(|v| {
            v.to_string()
                .chars()
                .filter(char::is_ascii_digit)
                .take(7)
                .collect()
        })
        .collect();
    fn scan(path: &std::path::Path, prefixes: &[String]) {
        for entry in std::fs::read_dir(path).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                scan(&path, prefixes);
            } else if path.extension().is_some_and(|e| e == "rs") {
                let source = std::fs::read_to_string(&path).unwrap();
                let digits: String = source.chars().filter(char::is_ascii_digit).collect();
                for prefix in prefixes {
                    assert!(
                        !digits.contains(prefix),
                        "typed math constant in {}",
                        path.display()
                    );
                }
            }
        }
    }
    scan(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src"),
        &prefixes,
    );
}
