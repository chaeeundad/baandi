//! Adapter against the independent `s-dftd3` command-line program from the same
//! conda-forge package (ticket D3-DISPERSION acceptance): the program reads the
//! structure from a Turbomole `coord` file in Bohr and writes energy, gradient
//! and virial to JSON; the adapter calls the C API in-process.
//!
//! Run: `BAANDI_DFTD3_LIB=<dest>/lib/libs-dftd3.dylib cargo test -p dft-dispersion
//! --features dftd3 --test oracle_sdftd3 -- --nocapture`. Skips when the
//! variable is unset.
#![cfg(feature = "dftd3")]

mod common;

use common::Fixture;
use dft_dispersion::{Damping, Overrides, Parameters, compute, symbol};
use std::path::Path;
use std::process::Command;

const E_TOL: f64 = 1e-10;
const G_TOL: f64 = 1e-9;
const S_TOL: f64 = 1e-9;

fn scratch(name: &str) -> std::path::PathBuf {
    let root = std::env::var_os("BAANDI_TEST_TMP")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(std::env::temp_dir);
    let dir = root.join(format!("dftd3-oracle-{}-{name}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn write_coord(f: &Fixture, dir: &Path) -> std::path::PathBuf {
    let mut text = String::from("$coord\n");
    for (z, r) in f.numbers.iter().zip(&f.positions) {
        text += &format!(
            " {:.17e} {:.17e} {:.17e} {}\n",
            r[0],
            r[1],
            r[2],
            symbol(*z).unwrap().to_ascii_lowercase()
        );
    }
    text += "$periodic 3\n$lattice bohr\n";
    for v in &f.lattice_rows {
        text += &format!(" {:.17e} {:.17e} {:.17e}\n", v[0], v[1], v[2]);
    }
    text += "$end\n";
    let path = dir.join("coord");
    std::fs::write(&path, text).unwrap();
    path
}

struct Oracle {
    energy: f64,
    gradient: Vec<f64>,
    virial: Vec<f64>,
}

fn run_program(program: &Path, dir: &Path, args: &[&str]) -> Oracle {
    let out = Command::new(program)
        .current_dir(dir)
        .args(["run"])
        .args(args)
        .args([
            "--noedisp",
            "--json",
            "out.json",
            "--grad",
            "grad.txt",
            "coord",
        ])
        .output()
        .expect("run s-dftd3");
    assert!(
        out.status.success(),
        "s-dftd3 {args:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let json: serde_json::Value =
        serde_json::from_slice(&std::fs::read(dir.join("out.json")).unwrap()).unwrap();
    let list = |k: &str| -> Vec<f64> {
        json[k]
            .as_array()
            .unwrap_or_else(|| panic!("s-dftd3 JSON has no {k}"))
            .iter()
            .map(|v| v.as_f64().unwrap())
            .collect()
    };
    Oracle {
        energy: json["energy"].as_f64().unwrap(),
        gradient: list("gradient"),
        virial: list("virial"),
    }
}

#[test]
fn adapter_matches_s_dftd3_program() {
    let Some(lib) = common::library() else { return };
    let Some(program) = common::program(&lib) else {
        return;
    };
    let explicit_bj = Overrides {
        s6: Some(1.0),
        s8: Some(0.7875),
        a1: Some(0.4289),
        a2: Some(4.4407),
        ..Overrides::default()
    };
    let explicit_zero = Overrides {
        s6: Some(1.0),
        s8: Some(0.722),
        rs6: Some(1.217),
        ..Overrides::default()
    };
    let none = Overrides::default();
    let no_atm = Overrides {
        s9: Some(0.0),
        ..none
    };
    // (label, program arguments, damping, overrides)
    let cases: Vec<(&str, Vec<&str>, Damping, Overrides)> = vec![
        (
            "pbe-d3bj-atm",
            vec!["--bj", "pbe", "--atm"],
            Damping::Rational,
            none,
        ),
        (
            "pbe-d3zero-atm",
            vec!["--zero", "pbe", "--atm"],
            Damping::Zero,
            none,
        ),
        ("pbe-d3bj", vec!["--bj", "pbe"], Damping::Rational, no_atm),
        ("pbe-d3zero", vec!["--zero", "pbe"], Damping::Zero, no_atm),
        (
            "explicit-bj-atm",
            vec!["--bj-param", "1.0", "0.7875", "0.4289", "4.4407", "--atm"],
            Damping::Rational,
            explicit_bj,
        ),
        (
            "explicit-zero-atm",
            vec!["--zero-param", "1.0", "0.722", "1.217", "--atm"],
            Damping::Zero,
            explicit_zero,
        ),
    ];
    let mut worst = (0.0_f64, 0.0_f64, 0.0_f64);
    println!(
        "{:18} {:18} {:>22} {:>10} {:>10} {:>10}",
        "structure", "case", "E_disp (Ha)", "|dE|", "max|dg|", "max|dsig|"
    );
    for f in common::all() {
        let dir = scratch(f.name);
        write_coord(&f, &dir);
        for (label, args, damping, overrides) in &cases {
            let p = Parameters::resolve(*damping, Some("pbe"), overrides).unwrap();
            let ours = compute(&lib, &f.structure(), &p).unwrap();
            let theirs = run_program(&program, &dir, args);
            let de = (ours.energy_ha - theirs.energy).abs();
            let dg = ours
                .gradient_ha_bohr
                .iter()
                .flatten()
                .zip(&theirs.gradient)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max);
            let ds = ours
                .virial_ha
                .iter()
                .flatten()
                .zip(&theirs.virial)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max);
            assert_eq!(theirs.gradient.len(), 3 * f.numbers.len());
            assert_eq!(theirs.virial.len(), 9);
            assert!(
                ours.energy_ha < 0.0,
                "{} {label}: E_disp must be negative",
                f.name
            );
            println!(
                "{:18} {:18} {:22.15e} {:10.2e} {:10.2e} {:10.2e}",
                f.name, label, ours.energy_ha, de, dg, ds
            );
            assert!(de <= E_TOL, "{} {label}: |dE| = {de:e}", f.name);
            assert!(dg <= G_TOL, "{} {label}: max |d grad| = {dg:e}", f.name);
            assert!(ds <= S_TOL, "{} {label}: max |d sigma| = {ds:e}", f.name);
            worst = (worst.0.max(de), worst.1.max(dg), worst.2.max(ds));
        }
        std::fs::remove_dir_all(&dir).ok();
    }
    println!(
        "worst: |dE| {:.2e} Ha, |d grad| {:.2e} Ha/Bohr, |d sigma| {:.2e} Ha",
        worst.0, worst.1, worst.2
    );
}

#[test]
fn library_errors_are_surfaced() {
    let Some(lib) = common::library() else { return };
    let f = common::graphite_ab();
    let p = Parameters::Library {
        damping: Damping::Rational,
        functional: "no-such-functional".into(),
        three_body: true,
    };
    let e = compute(&lib, &f.structure(), &p).unwrap_err();
    println!("library error: {e}");
    assert!(matches!(e, dft_dispersion::Error::Library(_)), "{e:?}");
    let missing = dft_dispersion::library_version(Path::new("/nonexistent/libs-dftd3.dylib"));
    assert!(matches!(missing, Err(dft_dispersion::Error::Load(_))));
    assert_eq!(dft_dispersion::library_version(&lib).unwrap(), "1.7.0");
}
