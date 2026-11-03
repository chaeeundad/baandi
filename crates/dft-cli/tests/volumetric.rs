use dft_cli::run::volumetric_out;
mod config {
    pub use dft_cli::config::*;
}
use dft_geometry::Cell;
use dft_physics::data::*;
use std::{fs, process::Command};
use volumetric_out::*;

#[test]
fn cube_independent_round_trip_and_xsf() {
    for handedness in [1., -1.] {
        let cell =
            Cell::from_user_rows_bohr([[4., 0.2, 0.], [1., 3., 0.], [0., 0.4, handedness * 2.]])
                .unwrap();
        let atoms = vec![OutputAtom {
            atomic_number: 14,
            valence_charge: 4.,
            position_bohr: [1., 2., 3.],
        }];
        let mut data =
            Array3::from_shape_fn([2, 3, 4], |(i, j, k)| (1 + i * 12 + j * 4 + k) as f64);
        let factor = 8.0 / (data.sum() * cell.volume_bohr3() / data.len() as f64);
        data.mapv_inplace(|x| x * factor);
        let mut bytes = Vec::new();
        write_cube(&mut bytes, &cell, &atoms, &data, "density").unwrap();
        let path =
            std::env::temp_dir().join(format!("p8b-{}-{handedness}.cube", std::process::id()));
        fs::write(&path, bytes).unwrap();
        let reader =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../scripts/read_cube.py");
        let output = Command::new("python3")
            .args(["-I"])
            .arg(reader)
            .arg(&path)
            .arg("--values")
            .output()
            .unwrap();
        fs::remove_file(path).unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let parsed: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
        for (v, json) in data.iter().zip(parsed["values"].as_array().unwrap()) {
            assert_eq!(v.to_bits(), json.as_f64().unwrap().to_bits());
        }
        assert_eq!(parsed["dims"], serde_json::json!([2, 3, 4]));
        assert_eq!(parsed["atoms"][0]["atomic_number"], 14);
        assert_eq!(parsed["atoms"][0]["charge"], 4.);
        assert_eq!(
            parsed["atoms"][0]["position_bohr"],
            serde_json::json!([1., 2., 3.])
        );
        let a = cell.lattice_columns_bohr();
        for axis in 0..3 {
            for (cart, row) in a.iter().enumerate() {
                assert!(
                    (parsed["lattice_vectors_bohr"][axis][cart].as_f64().unwrap() - row[axis])
                        .abs()
                        < 1e-14
                );
            }
        }
        let integral = parsed["integral"].as_f64().unwrap();
        assert!((integral - 8.).abs() < 1e-10);
        println!(
            "cube handedness={handedness}: integral={integral:.16} electrons; samples round trip bitwise"
        );
        let mut xsf = Vec::new();
        write_xsf(&mut xsf, &cell, &atoms, &data, "density").unwrap();
        let text = String::from_utf8(xsf).unwrap();
        let lines: Vec<_> = text.lines().collect();
        let start = lines
            .iter()
            .position(|l| *l == "BEGIN_DATAGRID_3D_density")
            .unwrap();
        assert_eq!(lines[start + 1], "3 4 5");
        let values: Vec<f64> = lines[start + 6..]
            .iter()
            .take_while(|l| **l != "END_DATAGRID_3D")
            .map(|l| l.parse().unwrap())
            .collect();
        assert_eq!(values.len(), 3 * 4 * 5);
        for k in 0..5 {
            for j in 0..4 {
                for i in 0..3 {
                    assert_eq!(values[(k * 4 + j) * 3 + i], data[[i % 2, j % 3, k % 4]]);
                }
            }
        }
        assert!(text.ends_with("END_DATAGRID_3D\nEND_BLOCK_DATAGRID_3D\n"));
    }
}

#[test]
fn symbols_and_selection() {
    assert_eq!(atomic_number("H").unwrap(), 1);
    assert_eq!(atomic_number("Si").unwrap(), 14);
    assert_eq!(atomic_number("Og").unwrap(), 118);
    for bad in ["si", "Si1", "Xx", ""] {
        assert!(atomic_number(bad).is_err());
    }
    let base = include_str!("../../../examples/si_lda.toml");
    let valid = "\n[output.volumetric]\nquantities=['density','elf','partial_charge']\nformats=['cube','xsf','vasp','vti']\nkpoint=1\nbands=[1,2]\n";
    assert!(config::Config::parse(&format!("{base}{valid}")).is_ok());
    for (from, to) in [
        ("kpoint=1", "kpoint=0"),
        ("bands=[1,2]", "bands=[2,1]"),
        ("'cube','xsf'", "'cube','cube'"),
        ("'density','elf','partial_charge'", "'unknown'"),
        ("formats=['cube','xsf','vasp','vti']", "formats=[]"),
    ] {
        assert!(config::Config::parse(&format!("{base}{}", valid.replace(from, to))).is_err());
    }
}

#[allow(dead_code)]
mod common;

fn synthetic_config(threads: usize) -> config::Config {
    config::Config::parse(&format!("[system]\nname='volumetric synthetic'\nlattice_unit='bohr'\nlattice_vectors=[[4.0,0.0,0.0],[0.0,4.0,0.0],[0.0,0.0,4.0]]\n[[system.atoms]]\nspecies='Si'\nfractional=[0.13,0.21,0.32]\n[pseudopotentials]\nSi='unused.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.125,0.0,0.0,0.3],[-0.125,0.0,0.0,0.7]]\nktime_reversal=false\n[xc]\nfunctional='lda-pz'\n[scf]\nthreads={threads}\n[output.volumetric]\nquantities=['density','elf','hartree','local_potential','partial_charge']\nformats=['cube','xsf','vasp','vti']\nbands=[1,2]\nkpoint=2\n")).unwrap()
}

#[test]
fn final_scf_exports_identical_across_thread_counts() {
    let root = std::env::temp_dir().join(format!("p8b-scf-{}", std::process::id()));
    let pp = dft_engine::input::parse_upf(
        &common::fixture()
            .replace("z_valence=\"4\"", "z_valence=\"2\"")
            .replace("Xx", "Si"),
    )
    .unwrap();
    let refs = [dft_engine::SpeciesPotential {
        label: "Si",
        potential: &pp,
    }];
    for threads in [1, 4] {
        let config = synthetic_config(threads);
        let prepared = config.prepare().unwrap();
        let atoms = config
            .system
            .atoms
            .iter()
            .map(|a| dft_engine::input::Atom::new(&a.species, a.fractional))
            .collect::<Vec<_>>();
        let points = config.basis.resolve_kpoints().unwrap();
        let result = dft_engine::ScfProblem::with_kpoints(
            &prepared.cell,
            &atoms,
            &refs,
            dft_cli::run::options(&config).unwrap(),
            &points,
            dft_engine::OccupationPolicy::Insulator {
                minimum_gap_hartree: 1e-7,
            },
        )
        .unwrap()
        .run()
        .unwrap();
        let out = root.join(threads.to_string());
        fs::create_dir_all(&out).unwrap();
        preflight(&config, &out).unwrap();
        write_outputs(&config, &prepared, &refs, &result, &out).unwrap();
        for (q, expected) in [("density", 2.0), ("partial_charge", 2.0)] {
            let output = Command::new("python3")
                .arg("-I")
                .arg(
                    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                        .join("../../scripts/read_cube.py"),
                )
                .arg(out.join(format!("{q}.cube")))
                .output()
                .unwrap();
            assert!(output.status.success());
            let parsed: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
            assert!((parsed["integral"].as_f64().unwrap() - expected).abs() < 1e-10);
        }
        assert!(preflight(&config, &out).is_err());
        assert!(write_outputs(&config, &prepared, &refs, &result, &out).is_err());
        for q in [
            "density",
            "elf",
            "hartree",
            "local_potential",
            "partial_charge",
        ] {
            for f in ["cube", "xsf", "vasp", "vti"] {
                let name = output_name(q, f);
                if threads == 4 {
                    assert_eq!(
                        fs::read(root.join("1").join(&name)).unwrap(),
                        fs::read(out.join(&name)).unwrap(),
                        "{name}"
                    );
                }
            }
        }
    }
    println!("SCF 2 k points: all 20 volumetric files bit-identical for threads=1,4");
    fs::remove_dir_all(root).unwrap();
}

#[test]
#[ignore = "real Si UPF required; run release with DFT_SI_UPF and DFT_VOLUMETRIC_OUTPUT"]
fn real_si_example() {
    let upf = std::env::var("DFT_SI_UPF").expect("set DFT_SI_UPF to a local NC LDA Si UPF");
    let out = std::env::var("DFT_VOLUMETRIC_OUTPUT")
        .expect("set DFT_VOLUMETRIC_OUTPUT to a fresh directory");
    let out = std::path::Path::new(&out);
    fs::create_dir_all(out).unwrap();
    let text = format!(
        "{}\n[output.volumetric]\nquantities=['density','elf','hartree','local_potential','partial_charge']\nformats=['cube','xsf','vasp','vti']\n",
        include_str!("../../../examples/si_lda.toml").replace("/absolute/path/to/Si.UPF", &upf)
    );
    let input = out.join("input.toml");
    use std::io::Write;
    let mut file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&input)
        .unwrap();
    file.write_all(text.as_bytes()).unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
        .arg("run")
        .arg(&input)
        .arg("--output-dir")
        .arg(out)
        .output()
        .unwrap();
    fs::write(out.join("run.stdout"), &result.stdout).unwrap();
    fs::write(out.join("run.stderr"), &result.stderr).unwrap();
    assert!(
        result.status.success(),
        "{}\n{}",
        String::from_utf8_lossy(&result.stdout),
        String::from_utf8_lossy(&result.stderr)
    );
    for q in [
        "density",
        "elf",
        "hartree",
        "local_potential",
        "partial_charge",
    ] {
        validate_formats(out, q);
        for f in ["cube", "xsf", "vasp", "vti"] {
            let path = out.join(output_name(q, f));
            assert!(path.is_file());
            println!("written: {}", path.display());
        }
    }
}

fn validate_formats(out: &std::path::Path, quantity: &str) {
    let output = Command::new("python3")
        .arg("-I")
        .arg(std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/validate_volumetric.py"))
        .arg(out.join(output_name(quantity, "cube")))
        .arg(out.join(output_name(quantity, "vasp")))
        .arg(out.join(output_name(quantity, "vti")))
        .arg(quantity)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    println!("{}", String::from_utf8_lossy(&output.stdout));
}

#[test]
fn vasp_vti_skew_handedness_and_atom_grouping() {
    let root = std::env::temp_dir().join(format!("p8b-formats-{}", std::process::id()));
    fs::create_dir_all(&root).unwrap();
    for handedness in [1.0, -1.0] {
        let cell =
            Cell::from_user_rows_bohr([[4., 0.2, 0.], [1., 3., 0.], [0., 0.4, handedness * 2.]])
                .unwrap();
        let atoms = vec![
            OutputAtom {
                atomic_number: 14,
                valence_charge: 4.,
                position_bohr: [1., 2., 3.],
            },
            OutputAtom {
                atomic_number: 8,
                valence_charge: 6.,
                position_bohr: [0.2, 0.3, 0.4],
            },
            OutputAtom {
                atomic_number: 14,
                valence_charge: 4.,
                position_bohr: [0.5, 0.6, 0.7],
            },
        ];
        // Index ramp catches transpose and periodic-endpoint mistakes.
        let data = Array3::from_shape_fn([2, 3, 4], |(i, j, k)| (i * 12 + j * 4 + k) as f64 / 24.0);
        for q in [
            "density",
            "elf",
            "hartree",
            "local_potential",
            "partial_charge",
        ] {
            let mut cube = Vec::new();
            write_cube(&mut cube, &cell, &atoms, &data, q).unwrap();
            let mut vasp = Vec::new();
            write_vasp(&mut vasp, &cell, &atoms, &data, q).unwrap();
            let mut vti = Vec::new();
            write_vti(&mut vti, &cell, &data, q).unwrap();
            fs::write(root.join(output_name(q, "cube")), cube).unwrap();
            fs::write(root.join(output_name(q, "vasp")), vasp).unwrap();
            fs::write(root.join(output_name(q, "vti")), vti).unwrap();
            validate_formats(&root, q);
        }
    }
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn normal_cli_emits_all_formats_and_protects_existing_files() {
    let root = std::env::temp_dir().join(format!("p8b-cli-{}", std::process::id()));
    fs::create_dir_all(&root).unwrap();
    fs::write(
        root.join("synthetic.xml"),
        common::fixture().replace("Xx", "Si"),
    )
    .unwrap();
    let input = root.join("input.toml");
    let text = "[system]\nname='synthetic'\nlattice_unit='bohr'\nlattice_vectors=[[4.0,0.0,0.0],[0.0,4.0,0.0],[0.0,0.0,4.0]]\n[[system.atoms]]\nspecies='Si'\nfractional=[0.13,0.21,0.32]\n[pseudopotentials]\nSi='synthetic.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[output.volumetric]\nquantities=['density','elf','hartree','local_potential','partial_charge']\nformats=['cube','xsf','vasp','vti']\n";
    fs::write(&input, text).unwrap();
    let out = root.join("out");
    let run = || {
        Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .arg("run")
            .arg(&input)
            .arg("--output-dir")
            .arg(&out)
            .output()
            .unwrap()
    };
    let result = run();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    for q in [
        "density",
        "elf",
        "hartree",
        "local_potential",
        "partial_charge",
    ] {
        validate_formats(&out, q);
    }
    let before = fs::read(out.join("density.cube")).unwrap();
    assert_eq!(run().status.code(), Some(5));
    assert_eq!(fs::read(out.join("density.cube")).unwrap(), before);
    // Legacy density.bin is exactly the same final grid as the cube.
    let bytes = fs::read(out.join("density.bin")).unwrap();
    let cube = fs::read_to_string(out.join("density.cube")).unwrap();
    let values: Vec<f64> = cube
        .lines()
        .skip(7)
        .flat_map(str::split_whitespace)
        .map(|x| x.parse().unwrap())
        .collect();
    assert_eq!(values.len() * 8, bytes.len());
    for (value, raw) in values.iter().zip(bytes.as_chunks::<8>().0) {
        assert_eq!(value.to_bits(), f64::from_le_bytes(*raw).to_bits());
    }
    fs::remove_dir_all(root).unwrap();
}
