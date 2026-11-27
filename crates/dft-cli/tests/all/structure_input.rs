//! SI1 structure input: compact `[structure]` block, POSCAR/extended XYZ
//! files and `[pseudopotentials] dir`.
//!
//! * Readers against an independent reader: `fixtures/structure_input/` holds
//!   files written by ASE 3.26 (or hand-written layouts) and `reference.json`,
//!   ASE's reading of each (`make_fixtures.py`, test-only, offline). BAANDI's
//!   native readers must give the same lattice (relative 1e-12), species,
//!   fractional positions (absolute 1e-12), Cartesian freeze masks and moments.
//! * Bit identity: the same structure written as `[system]`, compact
//!   fractional, compact Cartesian and a POSCAR file gives the same parsed
//!   bits and the same total-energy bits (synthetic potential, strict path).
//! * Relax constraints from compact flags and POSCAR selective dynamics equal
//!   `[relax] fixed` (one ionic step; identical constraints and energies).
//! * Moments: per-atom magmom maps onto `spin.initial_moment_seed`.
#[allow(dead_code)]
#[path = "../common/mod.rs"]
mod common;
use dft_cli::config::Config;
use serde_json::Value;
use std::{
    fs,
    path::{Path, PathBuf},
    process::Command,
};

struct Scratch(PathBuf);
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn scratch(tag: &str) -> Scratch {
    let dir = std::env::temp_dir().join(format!("baandi-si1-{tag}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    Scratch(dir)
}
fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/structure_input")
}
const REST: &str = "
[pseudopotentials]
dir = '.'
[basis]
ecutwfc_ry = 8.0
kpoints = [[0.0, 0.0, 0.0]]
[xc]
functional = 'lda-pz'
";
fn to_bohr(x: f64) -> f64 {
    dft_units::Bohr::from(dft_units::Angstrom(x)).0
}

#[test]
fn readers_match_ase() {
    let reference: Value =
        serde_json::from_slice(&fs::read(fixtures().join("reference.json")).unwrap()).unwrap();
    assert_eq!(reference["reader"], "ase");
    let cases = reference["cases"].as_array().unwrap();
    assert!(cases.len() >= 10);
    for case in cases {
        let file = case["file"].as_str().unwrap();
        let frame = case["frame"].as_i64().unwrap();
        let magnetic = case["magmoms"]
            .as_array()
            .is_some_and(|m| m.iter().any(|x| x.as_f64() != Some(0.0)));
        let frame_key = if case["format"] == "extxyz" {
            format!("frame = {frame}\n")
        } else {
            String::new()
        };
        let spin = if magnetic {
            "[spin]\nmode = 'collinear'\n"
        } else {
            ""
        };
        let rest = REST.replace("dir = '.'\n", "");
        let text = format!("[structure]\nfile = '{file}'\n{frame_key}{rest}{spin}");
        let config = Config::parse_unvalidated_at(&text, &fixtures().join("input.toml"))
            .unwrap_or_else(|e| panic!("{file}: {e}"));
        let tag = format!("{file} frame {frame}");
        // Lattice: rows in Angstrom (ASE) against the parsed cell in Bohr.
        let cols = config.cell().unwrap().lattice_columns_bohr();
        let rows = case["lattice_rows_angstrom"].as_array().unwrap();
        for v in 0..3 {
            for c in 0..3 {
                let ase = to_bohr(rows[v][c].as_f64().unwrap());
                let ours = cols[c][v];
                assert!(
                    (ase - ours).abs() <= 1e-12 * ase.abs().max(1.0),
                    "{tag}: lattice[{v}][{c}] {ours} vs ASE {ase}"
                );
            }
        }
        let species: Vec<&str> = config
            .system
            .atoms
            .iter()
            .map(|a| a.species.as_str())
            .collect();
        let ase_species: Vec<&str> = case["species"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap())
            .collect();
        assert_eq!(species, ase_species, "{tag}");
        for (i, (atom, ase)) in config
            .system
            .atoms
            .iter()
            .zip(case["fractional"].as_array().unwrap())
            .enumerate()
        {
            for c in 0..3 {
                let ase = ase[c].as_f64().unwrap();
                assert!(
                    (atom.fractional[c] - ase).abs() <= 1e-12,
                    "{tag}: atom {i} fractional[{c}] {} vs ASE {ase}",
                    atom.fractional[c]
                );
            }
        }
        let mut fixed = vec![None; species.len()];
        for f in &config.structure_fixed {
            fixed[f.atom] = Some(f.components);
        }
        let ase_fixed: Vec<Option<[bool; 3]>> = case["fixed"]
            .as_array()
            .unwrap()
            .iter()
            .map(|m| {
                m.as_array()
                    .map(|m| std::array::from_fn(|c| m[c].as_bool().unwrap()))
            })
            .collect();
        assert_eq!(fixed, ase_fixed, "{tag}: constraints");
        let ours = &config.structure_source.magmoms_mu_b;
        match case["magmoms"].as_array() {
            None => assert!(ours.is_empty(), "{tag}"),
            Some(m) => {
                let ase: Vec<Option<f64>> = m.iter().map(Value::as_f64).collect();
                assert_eq!(*ours, ase, "{tag}: moments");
                let seed = config.spin.as_ref().unwrap().initial_moment_seed;
                assert_eq!(
                    seed,
                    Some(ase.iter().map(|x| x.unwrap()).sum::<f64>()),
                    "{tag}"
                );
            }
        }
        let source = serde_json::to_value(&config.structure_source).unwrap();
        assert_eq!(source["kind"], "file");
        assert_eq!(source["file"]["format"], case["format"]);
        assert_eq!(source["file"]["sha256"].as_str().unwrap().len(), 64);
    }
}

const SYSTEM: &str = "[system]
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
const COMPACT: &str = "[structure]
name = 'Si2'
lattice_unit = 'angstrom'
lattice = [[0.0, 2.715, 2.715], [2.715, 0.0, 2.715], [2.715, 2.715, 0.0]]
coordinates = 'fractional'
atoms = '''
Si 0.0  0.0  0.0
Si 0.25 0.25 0.25   # second atom
'''
";
const CARTESIAN_SYSTEM: &str = "[system]
name = 'Si2'
lattice_unit = 'angstrom'
lattice_vectors = [[0.0, 2.715, 2.715], [2.715, 0.0, 2.715], [2.715, 2.715, 0.0]]
[[system.atoms]]
species = 'Si'
cartesian = [0.0, 0.0, 0.0]
[[system.atoms]]
species = 'Si'
cartesian = [1.3575, 1.3575, 1.3575]
";
const CARTESIAN_COMPACT: &str = "[structure]
name = 'Si2'
lattice_unit = 'angstrom'
lattice = [[0.0, 2.715, 2.715], [2.715, 0.0, 2.715], [2.715, 2.715, 0.0]]
coordinates = 'cartesian'
atoms = '''
Si 0 0 0
Si 1.3575 1.3575 1.3575
'''
";
const POSCAR: &str = "Si2
1.0
0.0 2.715 2.715
2.715 0.0 2.715
2.715 2.715 0.0
Si
2
Direct
0.0 0.0 0.0
0.25 0.25 0.25
";

fn baandi(dir: &Path, args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_dft-cli"))
        .env("BAANDI_BINARY_HASH_CACHE", "off")
        // Synthetic potential, as in tests/all/structure_convert.rs.
        .env("BAANDI_XC_SCREEN_NEGATIVE", "off")
        .env("BAANDI_REFERENCE_STRICT", "on")
        .current_dir(dir)
        .args(args)
        .output()
        .unwrap()
}

fn atoms_bits(c: &Config) -> Vec<(String, [u64; 3])> {
    c.system
        .atoms
        .iter()
        .map(|a| (a.species.clone(), a.fractional.map(f64::to_bits)))
        .collect()
}

/// Old and new input forms give the same parsed bits and the same energy bits.
#[test]
fn input_forms_are_bit_identical() {
    let s = scratch("forms");
    let dir = &s.0;
    let pp = dir.join("pp");
    fs::create_dir(&pp).unwrap();
    fs::write(
        pp.join("Si.pz-test.UPF"),
        common::gaussian().replace("element=\"Xx\"", "element=\"Si\""),
    )
    .unwrap();
    fs::write(dir.join("POSCAR"), POSCAR).unwrap();
    let explicit_pp = "[pseudopotentials]\nSi = 'pp/Si.pz-test.UPF'\n";
    let rest = REST.replace("[pseudopotentials]\ndir = '.'\n", "");
    let inputs = [
        ("system", format!("{SYSTEM}{explicit_pp}{rest}")),
        (
            "compact",
            format!("{COMPACT}[pseudopotentials]\ndir = 'pp'\n{rest}"),
        ),
        (
            "file",
            format!(
                "[structure]\nfile = 'POSCAR'\nname = 'Si2'\n[pseudopotentials]\ndir = 'pp'\n{rest}"
            ),
        ),
        (
            "cartesian_system",
            format!("{CARTESIAN_SYSTEM}{explicit_pp}{rest}"),
        ),
        (
            "cartesian_compact",
            format!("{CARTESIAN_COMPACT}{explicit_pp}{rest}"),
        ),
    ];
    let mut parsed = Vec::new();
    for (name, text) in &inputs {
        let path = dir.join(format!("{name}.toml"));
        fs::write(&path, text).unwrap();
        let c = Config::parse_at(text, &path).unwrap_or_else(|e| panic!("{name}: {e}"));
        parsed.push((
            c.cell()
                .unwrap()
                .lattice_columns_bohr()
                .map(|r| r.map(f64::to_bits)),
            atoms_bits(&c),
            c.system.name.clone(),
        ));
        assert_eq!(
            c.pseudopotentials["Si"],
            Path::new("pp/Si.pz-test.UPF"),
            "{name}"
        );
    }
    // Fractional forms agree with each other, Cartesian forms with each other.
    assert_eq!(parsed[0], parsed[1]);
    assert_eq!(parsed[0], parsed[2]);
    assert_eq!(parsed[3], parsed[4]);

    let mut energies = Vec::new();
    for (name, _) in &inputs[..3] {
        let input = format!("{name}.toml");
        let out = baandi(dir, &["validate", &input]);
        assert!(
            out.status.success(),
            "{name}: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let results = format!("results-{name}");
        let out = baandi(dir, &["run", &input, "--output-dir", &results]);
        assert!(
            out.status.success(),
            "{name}: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let json: Value =
            serde_json::from_slice(&fs::read(dir.join(&results).join("results.json")).unwrap())
                .unwrap();
        energies.push(json["total_energy"]["ha"].as_f64().unwrap().to_bits());
        if *name == "system" {
            // [system] inputs keep their results schema: no structure key.
            assert!(json.get("structure").is_none());
            assert!(json.get("pseudopotential_lookup").is_none());
            continue;
        }
        let structure = &json["structure"];
        assert_eq!(structure["schema"], "baandi-structure-1");
        assert_eq!(structure["source"]["kind"], *name);
        assert_eq!(structure["species"], serde_json::json!(["Si", "Si"]));
        assert_eq!(
            structure["fractional"][1],
            serde_json::json!([0.25, 0.25, 0.25])
        );
        let a = structure["lattice_rows_bohr"][0][1].as_f64().unwrap();
        assert_eq!(a.to_bits(), to_bohr(2.715).to_bits());
        if *name == "file" {
            let f = &structure["source"]["file"];
            assert_eq!(f["format"], "poscar");
            assert_eq!(f["path"], "POSCAR");
            assert_eq!(
                f["sha256"].as_str().unwrap(),
                format!(
                    "{:x}",
                    <sha2::Sha256 as sha2::Digest>::digest(POSCAR.as_bytes())
                )
            );
        }
        let l = &json["pseudopotential_lookup"];
        assert_eq!(l["dir"], "pp");
        assert_eq!(l["resolved"]["Si"], "pp/Si.pz-test.UPF");
        assert_eq!(l["rule"]["Si"], "prefix:Si");
    }
    assert_eq!(energies[0], energies[1], "compact energy bits differ");
    assert_eq!(energies[0], energies[2], "POSCAR energy bits differ");
}

/// Pseudopotential directory lookup order and errors.
#[test]
fn pseudopotential_dir_rules() {
    let s = scratch("ppdir");
    let dir = &s.0;
    let input = dir.join("input.toml");
    let parse = |atoms: &str, table: &str| {
        let text = format!(
            "[structure]\nlattice_unit='bohr'\nlattice=[[8.0,0,0],[0,8.0,0],[0,0,8.0]]\ncoordinates='fractional'\natoms='''\n{atoms}\n'''\n[pseudopotentials]\n{table}\n[basis]\necutwfc_ry=4.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n"
        );
        Config::parse_unvalidated_at(&text, &input).map(|c| {
            (
                c.pseudopotentials
                    .iter()
                    .map(|(k, v)| (k.clone(), v.to_string_lossy().into_owned()))
                    .collect::<Vec<_>>(),
                serde_json::to_value(c.pseudo_lookup).unwrap(),
            )
        })
    };
    let d = dir.join("lib");
    fs::create_dir(&d).unwrap();
    for f in [
        "Si.upf",
        "Si.UPF.bak",
        "fe.UPF",
        "O_pbe_v1.upf",
        "O.pbe-alt.UPF",
        "Mg-sp.upf",
        "readme.txt",
    ] {
        fs::write(d.join(f), "x").unwrap();
    }
    // Exact name.
    let (map, rec) = parse("Si 0 0 0", "dir='lib'").unwrap();
    assert_eq!(map, vec![("Si".into(), "lib/Si.upf".into())]);
    assert_eq!(rec["rule"]["Si"], "exact:Si");
    // Case-insensitive name, element prefix of a label.
    let (map, rec) = parse("Fe1 0 0 0\nMg 0.5 0.5 0.5", "dir='lib'").unwrap();
    assert_eq!(
        map,
        vec![
            ("Fe1".into(), "lib/fe.UPF".into()),
            ("Mg".into(), "lib/Mg-sp.upf".into())
        ]
    );
    assert_eq!(rec["rule"]["Fe1"], "case_insensitive:Fe");
    assert_eq!(rec["rule"]["Mg"], "prefix:Mg");
    // Explicit entries win and need no file.
    let (map, rec) = parse("O 0 0 0\nSi 0.5 0 0", "dir='lib'\nO='other/O.upf'").unwrap();
    assert_eq!(
        map,
        vec![
            ("O".into(), "other/O.upf".into()),
            ("Si".into(), "lib/Si.upf".into())
        ]
    );
    assert_eq!(rec["rule"]["O"], "explicit");
    // Two prefixed candidates: ambiguous; no candidate: missing.
    let err = parse("O 0 0 0", "dir='lib'").unwrap_err().to_string();
    assert!(
        err.contains("several files") && err.contains("explicit entry"),
        "{err}"
    );
    let err = parse("C 0 0 0", "dir='lib'").unwrap_err().to_string();
    assert!(
        err.contains("no pseudopotential for species \"C\""),
        "{err}"
    );
    let err = parse("Si 0 0 0", "dir='missing'").unwrap_err().to_string();
    assert!(err.contains("pseudopotentials.dir"), "{err}");
    let err = parse("Si 0 0 0", "dir=1").unwrap_err().to_string();
    assert!(err.contains("must be a string"), "{err}");
}

/// Input-level conflicts and errors are rejected with a message naming them.
#[test]
fn rejects_conflicts_and_bad_structures() {
    let s = scratch("errors");
    let dir = &s.0;
    let input = dir.join("input.toml");
    fs::write(dir.join("POSCAR"), POSCAR).unwrap();
    fs::write(dir.join("vasp4.POSCAR"), POSCAR.replace("Si\n2\n", "2\n")).unwrap();
    fs::write(
        dir.join("slab.extxyz"),
        "1\nLattice=\"5 0 0 0 5 0 0 0 20\" Properties=species:S:1:pos:R:3 pbc=\"T T F\"\nH 0 0 0\n",
    )
    .unwrap();
    fs::write(
        dir.join("nolattice.xyz"),
        "1\nProperties=species:S:1:pos:R:3\nH 0 0 0\n",
    )
    .unwrap();
    fs::write(dir.join("structure.cif"), "data_x\n").unwrap();
    let rest = REST.replace("dir = '.'", "Si = 'x.upf'\nH = 'x.upf'");
    let cases: Vec<(String, &str)> = vec![
        (
            format!("{SYSTEM}{COMPACT}{rest}"),
            "[structure] and [system] both given",
        ),
        (
            format!(
                "{}{rest}[relax]\nfixed=[{{atom=0}}]\n",
                COMPACT.replace("Si 0.0  0.0  0.0", "Si 0 0 0 fixed=xyz")
            ),
            "relax.fixed and per-atom fixed flags",
        ),
        (
            format!(
                "{}{rest}",
                COMPACT.replace("Si 0.0  0.0  0.0", "Si 0 0 0 magmom=1")
            ),
            "need [spin] mode = \"collinear\"",
        ),
        (
            format!(
                "{}{rest}[spin]\nmode='collinear'\ninitial_moment_seed=1.0\n",
                COMPACT.replace("Si 0.0  0.0  0.0", "Si 0 0 0 magmom=1")
            ),
            "spin.initial_moment_seed and per-atom magmom both given",
        ),
        (
            format!("{}file = 'POSCAR'\n{rest}", COMPACT),
            "structure.file conflicts with",
        ),
        (
            format!("{}{rest}", COMPACT.replace("name = 'Si2'", "typo = 1")),
            "unknown field",
        ),
        (
            format!(
                "{}{rest}",
                COMPACT.replace("coordinates = 'fractional'\n", "")
            ),
            "structure.coordinates is required",
        ),
        (
            format!("{}{rest}", COMPACT.replace("'fractional'", "'crystal'")),
            "use \"fractional\" or \"cartesian\"",
        ),
        (
            format!("{}{rest}", COMPACT.replace("'angstrom'", "'pm'")),
            "use \"angstrom\" or \"bohr\"",
        ),
        (
            format!(
                "{}{rest}",
                COMPACT.replace("Si 0.25 0.25 0.25", "Si 0.25 0.25")
            ),
            "expected `SPECIES c1 c2 c3",
        ),
        (
            format!(
                "{}{rest}",
                COMPACT.replace("Si 0.25 0.25 0.25", "Si 0.25 0.25 0.25 fixed=xq")
            ),
            "unknown component",
        ),
        (
            format!(
                "{}{rest}",
                COMPACT.replace("Si 0.25 0.25 0.25", "Si 0.25 0.25 0.25 spin=up")
            ),
            "unknown flag",
        ),
        (
            format!("[structure]\nfile='vasp4.POSCAR'\n{rest}"),
            "VASP 4 layout",
        ),
        (
            format!("[structure]\nfile='slab.extxyz'\n{rest}"),
            "not periodic in every direction",
        ),
        (
            format!("[structure]\nfile='nolattice.xyz'\n{rest}"),
            "no Lattice",
        ),
        (
            format!("[structure]\nfile='structure.cif'\n{rest}"),
            "baandi-convert",
        ),
        (
            format!("[structure]\nfile='missing.vasp'\n{rest}"),
            "missing.vasp",
        ),
        (
            format!("[structure]\nfile='POSCAR'\nframe=0\n{rest}"),
            "extended XYZ files only",
        ),
        (
            format!("[structure]\nfile='slab.extxyz'\nframe=3\n{rest}"),
            "frame(s)",
        ),
    ];
    for (text, needle) in cases {
        let err = Config::parse_unvalidated_at(&text, &input)
            .map(|_| ())
            .unwrap_err()
            .to_string();
        assert!(err.contains(needle), "expected {needle:?} in {err:?}");
    }
}

/// Per-atom magmom flags give the same spin seed as the explicit key.
#[test]
fn magmom_maps_to_initial_moment_seed() {
    let spin = "[spin]\nmode='collinear'\nn_up=5\nn_down=3\n";
    let rest = REST
        .replace("dir = '.'", "Si = 'x.upf'")
        .replace("'lda-pz'", "'lda-pw'");
    let compact = COMPACT
        .replace("Si 0.0  0.0  0.0", "Si 0 0 0 magmom=1.5")
        .replace("Si 0.25 0.25 0.25 ", "Si 0.25 0.25 0.25 magmom=0.5");
    let a = Config::parse_unvalidated(&format!("{compact}{rest}{spin}")).unwrap();
    let b = Config::parse_unvalidated(&format!("{SYSTEM}{rest}{spin}initial_moment_seed=2.0\n"))
        .unwrap();
    assert_eq!(a.spin.as_ref().unwrap().initial_moment_seed, Some(2.0));
    assert_eq!(
        a.spin.as_ref().unwrap().initial_moment_seed,
        b.spin.as_ref().unwrap().initial_moment_seed
    );
    assert_eq!(a.structure_source.magmoms_mu_b, vec![Some(1.5), Some(0.5)]);
    assert_eq!(
        a.structure_source.initial_moment_seed_from_magmoms,
        Some(2.0)
    );
    assert_eq!(atoms_bits(&a), atoms_bits(&b));
    // All-zero moments need no [spin] table and add no seed.
    let zero = COMPACT.replace("Si 0.0  0.0  0.0", "Si 0 0 0 magmom=0");
    let c = Config::parse_unvalidated(&format!("{zero}{rest}")).unwrap();
    assert!(c.spin.is_none());
    assert_eq!(c.structure_source.magmoms_mu_b, vec![Some(0.0), None]);
}

fn relax_input(structure: &str, relax_extra: &str) -> String {
    format!(
        "{structure}[pseudopotentials]\nH='syn.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[scf]\nnbands=4\nmax_iterations=300\n[relax]\nmax_steps=1\nforce_tol_ha_bohr=1e-8\n{relax_extra}"
    )
}

/// Relax constraints from compact flags and from POSCAR selective dynamics
/// equal `[relax] fixed`: the same recorded constraints and the same
/// frames (one ionic step, exit 3 = step budget exhausted).
#[test]
fn relax_constraints_from_flags_and_selective_dynamics() {
    let s = scratch("relax");
    let dir = &s.0;
    fs::write(
        dir.join("syn.xml"),
        common::gaussian_with(1., 1.5, false, 401).replace("element=\"Xx\"", "element=\"H\""),
    )
    .unwrap();
    let system = "[system]\nname='H2'\nlattice_unit='bohr'\nlattice_vectors=[[8.0,0.0,0.0],[0.0,8.0,0.0],[0.0,0.0,8.0]]\n[[system.atoms]]\nspecies='H'\nfractional=[0.0,0.0,0.0]\n[[system.atoms]]\nspecies='H'\nfractional=[0.27,0.23,0.19]\n";
    let compact = "[structure]\nname='H2'\nlattice_unit='bohr'\nlattice=[[8.0,0.0,0.0],[0.0,8.0,0.0],[0.0,0.0,8.0]]\ncoordinates='fractional'\natoms='''\nH 0.0 0.0 0.0 fixed=true\nH 0.27 0.23 0.19 fixed=y\n'''\n";
    // Bohr cell written in Angstrom-scaled POSCAR: scale = 1 Bohr in Angstrom.
    let poscar = format!(
        "H2\n{}\n8.0 0.0 0.0\n0.0 8.0 0.0\n0.0 0.0 8.0\nH\n2\nSelective dynamics\nDirect\n0.0 0.0 0.0 F F F\n0.27 0.23 0.19 T F T\n",
        1.0 / to_bohr(1.0)
    );
    fs::write(dir.join("POSCAR"), &poscar).unwrap();
    let inputs = [
        (
            "legacy",
            relax_input(
                system,
                "fixed=[{atom=0},{atom=1,components=[false,true,false]}]\n",
            ),
        ),
        ("compact", relax_input(compact, "")),
        (
            "poscar",
            relax_input("[structure]\nname='H2'\nfile='POSCAR'\n", ""),
        ),
    ];
    let mut frames = Vec::new();
    for (name, text) in &inputs {
        let input = format!("{name}.toml");
        fs::write(dir.join(&input), text).unwrap();
        let out = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .env("BAANDI_BINARY_HASH_CACHE", "off")
            .current_dir(dir)
            .args(["relax", &input, "--output-dir", &format!("out-{name}")])
            .output()
            .unwrap();
        assert_eq!(
            out.status.code(),
            Some(3),
            "{name}: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let r: Value = serde_json::from_slice(
            &fs::read(dir.join(format!("out-{name}")).join("results.json")).unwrap(),
        )
        .unwrap();
        assert_eq!(
            r["relaxation"]["constraints"],
            serde_json::json!([{"atom":0,"components":[true,true,true]},{"atom":1,"components":[false,true,false]}]),
            "{name}"
        );
        frames.push((
            r["relaxation"]["energy_history_ha"].clone(),
            r["relaxation"]["frames"]
                .as_array()
                .unwrap()
                .iter()
                .map(|f| f["fractional"].clone())
                .collect::<Vec<_>>(),
        ));
        if *name != "legacy" {
            assert_eq!(
                r["structure"]["fixed"],
                serde_json::json!([[true, true, true], [false, true, false]]),
                "{name}"
            );
        }
    }
    assert_eq!(
        frames[0], frames[1],
        "compact flags differ from relax.fixed"
    );
    // The POSCAR cell is 8 Bohr only to roundoff of the Angstrom scale, so
    // the energies agree to roundoff, not bits.
    let e = |i: usize| {
        frames[i]
            .0
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_f64().unwrap())
            .collect::<Vec<_>>()
    };
    for (a, b) in e(0).iter().zip(e(2)) {
        assert!((a - b).abs() < 1e-9, "{a} vs {b}");
    }
    // The constrained component of atom 1 (y) and atom 0 never move.
    let f = &frames[2].1;
    assert!(f.len() >= 2, "no ionic step was taken");
    assert_ne!(f[0][1], f.last().unwrap()[1], "atom 1 did not move");
    assert_eq!(f[0][0], f.last().unwrap()[0]);
    assert_eq!(f[0][1][1], f.last().unwrap()[1][1]);
}
