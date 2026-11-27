//! A `[system]` block written by the workstation converter
//! (`python/baandi-convert`, committed fixture produced offline from a
//! synthetic cached record) is an ordinary explicit input: it validates and
//! gives bit-identical energies to the hand-written equivalent, and its
//! `[system.provenance]` table is copied into results.json unchanged.
#[allow(dead_code)]
#[path = "../common/mod.rs"]
mod common;
use serde_json::Value;
use std::{fs, path::PathBuf, process::Command};

const CONVERTED: &str = include_str!("../fixtures/structure_convert/si_mp0_converted.toml");
const HAND_WRITTEN: &str = "[system]
name = 'Si2'
lattice_unit = 'angstrom'
lattice_vectors = [[0.0, 2.715, 2.715], [2.715, 0.0, 2.715], [2.715, 2.715, 0.0]]
[[system.atoms]]
species = 'Si'
fractional = [0.0, 0.0, 0.0]
[[system.atoms]]
species = 'Si'
fractional = [0.25, 0.25, 0.25]
";
const REST: &str = "
[pseudopotentials]
Si = 'synthetic.xml'
[basis]
ecutwfc_ry = 8.0
kpoints = [[0.0, 0.0, 0.0]]
[xc]
functional = 'lda-pz'
";

struct Scratch(PathBuf);
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn baandi(dir: &std::path::Path, args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_dft-cli"))
        .env("BAANDI_BINARY_HASH_CACHE", "off")
        // The synthetic potential is not a material potential; its density dips
        // below the XC screening limit, so screening is off as in tests/all/run.rs.
        .env("BAANDI_XC_SCREEN_NEGATIVE", "off")
        // Same converging synthetic case as the large case of tests/all/run.rs.
        .env("BAANDI_REFERENCE_STRICT", "on")
        .current_dir(dir)
        .args(args)
        .output()
        .unwrap()
}

#[test]
fn converted_input_matches_hand_written_bits_and_records_provenance() {
    let scratch = Scratch(
        std::env::temp_dir().join(format!("baandi-structure-convert-{}", std::process::id())),
    );
    fs::create_dir_all(&scratch.0).unwrap();
    let dir = &scratch.0;
    fs::write(
        dir.join("synthetic.xml"),
        common::gaussian().replace("element=\"Xx\"", "element=\"Si\""),
    )
    .unwrap();
    fs::write(dir.join("converted.toml"), format!("{CONVERTED}{REST}")).unwrap();
    fs::write(dir.join("hand.toml"), format!("{HAND_WRITTEN}{REST}")).unwrap();

    // Same parsed structure, bit for bit.
    let parse = |name: &str| {
        dft_cli::config::Config::parse(&fs::read_to_string(dir.join(name)).unwrap()).unwrap()
    };
    let (converted, hand) = (parse("converted.toml"), parse("hand.toml"));
    assert_eq!(
        converted.cell().unwrap().lattice_columns_bohr(),
        hand.cell().unwrap().lattice_columns_bohr()
    );
    for (a, b) in converted.system.atoms.iter().zip(&hand.system.atoms) {
        assert_eq!(a.species, b.species);
        assert_eq!(
            a.fractional.map(f64::to_bits),
            b.fractional.map(f64::to_bits)
        );
    }
    assert!(hand.system.provenance_json().is_none());

    let mut energies = Vec::new();
    for name in ["converted", "hand"] {
        let input = format!("{name}.toml");
        let out = baandi(dir, &["validate", &input]);
        assert!(
            out.status.success(),
            "{name} validate: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let results = format!("results-{name}");
        let out = baandi(dir, &["run", &input, "--output-dir", &results]);
        assert!(
            out.status.success(),
            "{name} run: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let json: Value =
            serde_json::from_slice(&fs::read(dir.join(&results).join("results.json")).unwrap())
                .unwrap();
        energies.push(json["total_energy"]["ha"].as_f64().unwrap().to_bits());
        if name == "converted" {
            let p = &json["structure_provenance"];
            assert_eq!(p["source_kind"], "materials_project");
            assert_eq!(p["mp_id"], "mp-0");
            assert_eq!(p["converter"], "baandi-convert");
            assert_eq!(p["source_sha256"].as_str().unwrap().len(), 64);
        } else {
            assert!(json.get("structure_provenance").is_none());
        }
    }
    assert_eq!(energies[0], energies[1], "energy bits differ");
}

#[test]
fn provenance_rejects_nested_tables() {
    let base = format!("{CONVERTED}{REST}").replace("synthetic.xml", "/nonexistent.xml");
    for extra in ["nested = { a = 1 }", "list = [{ a = 1 }]"] {
        let text = base.replace(
            "[system.provenance]\n",
            &format!("[system.provenance]\n{extra}\n"),
        );
        let err = dft_cli::config::Config::parse(&text)
            .unwrap_err()
            .to_string();
        assert!(err.contains("system.provenance"), "{err}");
    }
    // Unknown keys elsewhere in [system] stay errors.
    assert!(
        dft_cli::config::Config::parse(&base.replace("[system]\n", "[system]\ntypo = 1\n"))
            .is_err()
    );
}

/// SI1: `baandi-convert --compact` writes the `[structure]` form of the same
/// record (committed fixture); it parses to the same bits as the `[system]`
/// output and keeps the provenance table.
#[test]
fn compact_converter_output_matches_system_form_bits() {
    const COMPACT: &str = include_str!("../fixtures/structure_convert/si_mp0_compact.toml");
    let rest = REST.replace("synthetic.xml", "/nonexistent.xml");
    let compact = dft_cli::config::Config::parse(&format!("{COMPACT}{rest}")).unwrap();
    let system = dft_cli::config::Config::parse(&format!("{CONVERTED}{rest}")).unwrap();
    assert_eq!(
        compact.cell().unwrap().lattice_columns_bohr(),
        system.cell().unwrap().lattice_columns_bohr()
    );
    assert_eq!(compact.system.name, system.system.name);
    for (a, b) in compact.system.atoms.iter().zip(&system.system.atoms) {
        assert_eq!(a.species, b.species);
        assert_eq!(
            a.fractional.map(f64::to_bits),
            b.fractional.map(f64::to_bits)
        );
    }
    let (p, q) = (
        compact.system.provenance_json().unwrap(),
        system.system.provenance_json().unwrap(),
    );
    assert_eq!(p, q);
    assert_eq!(compact.structure_source.kind, "compact");
}
