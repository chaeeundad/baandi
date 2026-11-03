mod common;
use dft_cli::{config::Config, run::options};
use dft_engine::{
    ScfProblem, SpeciesPotential,
    input::{Atom, parse_upf},
};
use serde_json::Value;
use sha2::{Digest, Sha256};
use std::{
    fs,
    path::{Path, PathBuf},
    process::{Command, Output},
    sync::atomic::{AtomicUsize, Ordering},
};
static SERIAL: AtomicUsize = AtomicUsize::new(0);
struct Case {
    dir: PathBuf,
    text: String,
    fixture: String,
}
impl Case {
    fn new(large: bool) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "dft-run-{}-{}",
            std::process::id(),
            SERIAL.fetch_add(1, Ordering::SeqCst)
        ));
        fs::create_dir_all(&dir).unwrap();
        let fixture = if large {
            common::gaussian()
        } else {
            common::fixture()
        };
        fs::write(dir.join("synthetic.xml"), &fixture).unwrap();
        let text = if large {
            include_str!("../../../examples/si_lda.toml")
                .replace("/absolute/path/to/Si.UPF", "synthetic.xml")
                .replace("30.0", "8.0")
        } else {
            "[system]\nname='synthetic small'\nlattice_unit='bohr'\nlattice_vectors=[[4.0,0.0,0.0],[0.0,4.0,0.0],[0.0,0.0,4.0]]\n[[system.atoms]]\nspecies='Si'\nfractional=[0.13,0.21,0.32]\n[pseudopotentials]\nSi='synthetic.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n".into()
        };
        Self { dir, text, fixture }
    }
    fn run(&self, suffix: &str, override_dir: Option<&Path>) -> Output {
        fs::write(
            self.dir.join("input.toml"),
            format!("{}{suffix}", self.text),
        )
        .unwrap();
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_dft-cli"));
        cmd.current_dir(&self.dir).args(["run", "input.toml"]);
        if let Some(dir) = override_dir {
            cmd.arg("--output-dir").arg(dir);
        }
        cmd.output().unwrap()
    }
    fn json(&self, dir: &str, file: &str) -> Value {
        serde_json::from_slice(&fs::read(self.dir.join(dir).join(file)).unwrap()).unwrap()
    }
}
impl Drop for Case {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.dir).unwrap();
    }
}
fn status(o: &Output, code: i32) {
    assert_eq!(
        o.status.code(),
        Some(code),
        "stdout={}\nstderr={}",
        String::from_utf8_lossy(&o.stdout),
        String::from_utf8_lossy(&o.stderr)
    );
    if code != 0 {
        assert!(!o.stderr.is_empty());
    }
}
fn hash(p: &Path) -> String {
    format!("{:x}", Sha256::digest(fs::read(p).unwrap()))
}
fn strip_times(v: &mut Value) {
    match v {
        Value::Object(m) => {
            m.remove("timestamp");
            m.remove("timings");
            for x in m.values_mut() {
                strip_times(x);
            }
        }
        Value::Array(a) => {
            for x in a {
                strip_times(x);
            }
        }
        _ => {}
    }
}
fn verify(large: bool) {
    let c = Case::new(large);
    let run = c.run("", None);
    status(&run, 0);
    let mut first = c.json("results", "results.json");
    // Exact legacy top-level shape, checked on actual CLI output.
    assert_eq!(first.as_object().unwrap().len(), 22);
    assert!(first.get("kpoints").is_none());
    for key in [
        "schema_version",
        "tool",
        "timestamp",
        "input_sha256",
        "species",
        "cell",
        "atoms",
        "cutoff",
        "basis",
        "electrons",
        "bands",
        "occupations",
        "xc",
        "scf",
        "convergence",
        "energies",
        "eigenvalues",
        "timings",
        "threading",
        "memory_estimate",
        "units",
    ] {
        assert!(!first[key].is_null(), "{key}");
    }
    assert_eq!(first["convergence"]["converged"], true);
    let energy = first["energies"]["total"]["ha"].as_f64().unwrap();
    assert_eq!(
        first["total_energy"]["ha"],
        first["energies"]["total"]["ha"]
    );
    assert_eq!(
        first["total_energy"]["ry"],
        first["energies"]["total"]["ry"]
    );
    assert_eq!(
        first["total_energy"]["ev"].as_f64().unwrap(),
        dft_units::from_hartree(energy, dft_units::EnergyUnit::ElectronVolt)
    );
    let sum: f64 = ["kinetic", "local", "nonlocal", "hartree", "xc", "ewald"]
        .iter()
        .map(|t| first["energies"][t]["ha"].as_f64().unwrap())
        .sum();
    assert!((sum - energy).abs() < 1e-12); // Ha
    assert_eq!(first["input_sha256"], hash(&c.dir.join("input.toml")));
    assert_eq!(
        first["species"][0]["sha256"],
        hash(&c.dir.join("synthetic.xml"))
    );
    let config = Config::parse(&c.text).unwrap();
    let pp = parse_upf(&c.fixture).unwrap();
    let atoms: Vec<_> = config
        .system
        .atoms
        .iter()
        .map(|a| Atom::new(&a.species, a.fractional))
        .collect();
    let direct = ScfProblem::new(
        &config.cell().unwrap(),
        &atoms,
        &[SpeciesPotential {
            label: "Si",
            potential: &pp,
        }],
        options(&config).unwrap(),
    )
    .unwrap()
    .run()
    .unwrap();
    assert!((direct.total_hartree() - energy).abs() < 1e-12); // Ha, independent engine invocation
    let meta = c.json("results", "density.json");
    let data = fs::read(c.dir.join("results/density.bin")).unwrap();
    let density: Vec<f64> = data
        .as_chunks::<8>()
        .0
        .iter()
        .map(|b| f64::from_le_bytes(*b))
        .collect();
    let dims: Vec<usize> = meta["fft_grid"]
        .as_array()
        .unwrap()
        .iter()
        .map(|n| n.as_u64().unwrap() as usize)
        .collect();
    assert_eq!(density.len(), dims.iter().product::<usize>());
    assert_eq!(data.len(), density.len() * 8);
    assert_eq!(meta["endianness"], "little");
    assert_eq!(meta["dtype"], "f64");
    assert_eq!(meta["fft_grid"], first["basis"]["fft_grid"]);
    assert_eq!(
        meta["lattice_columns_bohr"],
        first["cell"]["lattice_columns_bohr"]
    );
    assert_eq!(meta["sha256"], hash(&c.dir.join("results/density.bin")));
    let integral =
        density.iter().sum::<f64>() * meta["volume_bohr3"].as_f64().unwrap() / density.len() as f64;
    let electrons = first["electrons"].as_f64().unwrap();
    assert_eq!(meta["normalization"]["expected_electrons"], electrons);
    assert!(
        (meta["normalization"]["integral_electrons"]
            .as_f64()
            .unwrap()
            - integral)
            .abs()
            < 1e-12
    );
    assert!((integral - electrons).abs() / electrons < 1e-10);
    assert!(
        first["convergence"]["last_iteration"]["relative_charge_error"]
            .as_f64()
            .unwrap()
            < 1e-8
    );
    for [i, j, k] in [
        [0, 0, 0],
        [1, 2, 1],
        [dims[0] - 1, dims[1] - 1, dims[2] - 1],
    ] {
        assert_eq!(
            density[(i * dims[1] + j) * dims[2] + k],
            direct.density[[i, j, k]]
        );
    }
    println!(
        "large={large}; energy={energy:.16} Ha; energy sum error={:.3e} Ha; direct error={:.3e} Ha; charge error={:.3e} electrons; iterations={}",
        (sum - energy).abs(),
        (direct.total_hartree() - energy).abs(),
        (integral - electrons).abs(),
        first["convergence"]["iterations"]
    );
    status(&c.run("", Some(Path::new("second"))), 0);
    let mut second = c.json("second", "results.json");
    strip_times(&mut first);
    strip_times(&mut second);
    assert_eq!(
        serde_json::to_vec_pretty(&first).unwrap(),
        serde_json::to_vec_pretty(&second).unwrap()
    );
    assert_eq!(data, fs::read(c.dir.join("second/density.bin")).unwrap());
}
#[test]
fn small_end_to_end() {
    verify(false);
}
#[test]
#[ignore = "two-atom 8 Ry SCF: run in release; see docs/cli.md"]
fn fcc_end_to_end() {
    verify(true);
}
#[test]
fn failures_and_paths() {
    let mut c = Case::new(false);
    status(&c.run("\n[scf]\nmax_iterations=2", None), 3);
    assert_eq!(
        c.json("results", "results.json")["convergence"]["converged"],
        false
    );
    let o = c.run("\n[scf]\nmax_iterations=2", None);
    assert!(!String::from_utf8_lossy(&o.stdout).contains("converged"));
    status(
        &c.run("\n[scf]\nmixing='linear'\nmixing_beta=100.0", None),
        4,
    );
    status(&c.run("\n[scf]\ntypo=1", None), 1);
    c.text = c
        .text
        .replace("kpoints=[[0.0,0.0,0.0]]", "kpoints=[[0.1,0.0,0.0]]");
    status(&c.run("", None), 1);
    c.text = c
        .text
        .replace("kpoints=[[0.1,0.0,0.0]]", "kpoints=[[0.0,0.0,0.0]]");
    for (from, to) in [
        ("is_paw=\"F\"", "is_paw=\"T\""),
        ("is_ultrasoft=\"F\"", "is_ultrasoft=\"T\""),
        ("SLA PZ", "PBE"),
    ] {
        fs::write(
            c.dir.join("synthetic.xml"),
            common::fixture().replace(from, to),
        )
        .unwrap();
        status(&c.run("", None), 2);
    }
    fs::remove_file(c.dir.join("synthetic.xml")).unwrap();
    status(&c.run("", None), 5);
    fs::write(c.dir.join("synthetic.xml"), common::fixture()).unwrap();
    status(&c.run("\n[output]\ndirectory='relative/nested'", None), 0);
    assert!(c.dir.join("relative/nested/results.json").exists());
    let configured_absolute = c.dir.join("configured-absolute");
    status(
        &c.run(
            &format!("\n[output]\ndirectory='{}'", configured_absolute.display()),
            None,
        ),
        0,
    );
    assert!(configured_absolute.join("results.json").exists());
    let absolute = c.dir.join("absolute");
    status(&c.run("\n[output]\ndirectory='unused'", Some(&absolute)), 0);
    assert!(!c.dir.join("unused").exists());
    status(&c.run("", Some(&absolute)), 0); // prior outputs may be overwritten
    fs::write(c.dir.join("blocked"), "file").unwrap();
    status(&c.run("", Some(Path::new("blocked"))), 5);
}
#[test]
#[cfg(unix)]
fn protects_inputs() {
    use std::os::unix::fs::symlink;
    let c = Case::new(false);
    fs::create_dir(c.dir.join("out")).unwrap();
    fs::write(c.dir.join("input.toml"), &c.text).unwrap();
    for link in [false, true] {
        let dest = c.dir.join("out/density.bin");
        if link {
            symlink(c.dir.join("synthetic.xml"), &dest).unwrap();
        } else {
            fs::hard_link(c.dir.join("synthetic.xml"), &dest).unwrap();
        }
        status(&c.run("", Some(Path::new("out"))), 5);
        assert_eq!(
            fs::read_to_string(c.dir.join("synthetic.xml")).unwrap(),
            common::fixture()
        );
        fs::remove_file(dest).unwrap();
    }
}

#[test]
fn functional_compatibility_and_original_label() {
    let c = Case::new(false);
    for label in ["SLA  PZ   NOGX NOGC", " SLA  PZ   NOGX NOGC "] {
        fs::write(
            c.dir.join("synthetic.xml"),
            c.fixture.replace("SLA PZ", label),
        )
        .unwrap();
        status(&c.run("", None), 0);
        let results = c.json("results", "results.json");
        assert_eq!(results["convergence"]["converged"], true);
        assert_eq!(results["species"][0]["functional"], label);
    }
    fs::write(
        c.dir.join("synthetic.xml"),
        c.fixture.replace("SLA PZ", "SLA PW"),
    )
    .unwrap();
    let output = c.run("", None);
    status(&output, 2);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr
            .contains("unsupported UPF functional \"SLA PW\": requires SLA PZ or SLA PZ NOGX NOGC"),
        "{stderr}"
    );
}

#[test]
fn unit_cutoff_preserves_scf_and_legacy_outputs() {
    let mut c = Case::new(false);
    status(&c.run("", None), 0);
    let mut legacy = c.json("results", "results.json");
    let density = fs::read(c.dir.join("results/density.bin")).unwrap();
    c.text = c.text.replace("ecutwfc_ry=3.0", "cutoff='3 Ry'");
    status(&c.run("", None), 0);
    let mut units = c.json("results", "results.json");
    assert_ne!(legacy["input_sha256"], units["input_sha256"]);
    for report in [&mut legacy, &mut units] {
        strip_times(report);
        report.as_object_mut().unwrap().remove("input_sha256");
    }
    assert_eq!(legacy, units); // All old and new numerical fields, exactly.
    assert_eq!(
        density,
        fs::read(c.dir.join("results/density.bin")).unwrap()
    );
    c.text = c.text.replace("cutoff='3 Ry'", "cutoff='3 unknown'");
    status(&c.run("", None), 1);
    c.text = c
        .text
        .replace("cutoff='3 unknown'", "cutoff='3 Ry'\necutwfc_ry=3.0");
    status(&c.run("", None), 1);
}

#[test]
fn kpoint_inputs_and_errors() {
    let mut c = Case::new(false);
    let original = c.text.clone();
    for (input, message) in [
        ("kgrid=[0,1,1]", "grid"),
        ("kgrid=[2,1,1]\nkshift=[0.1,0.0,0.0]", "shift"),
        ("kpoints=[[0.0,0.0,0.0,0.9]]", "weight"),
        ("kpoints=[[nan,0.0,0.0,1.0]]", "finite"),
        ("kpoints=[[0,0,0,inf]]", "weight"),
        ("kpoints=[[0,0,0,-1]]", "weight"),
        ("kpoints=[]", "empty"),
        ("kgrid=[1,1,1]\nkshift=[nan,0,0]", "shift"),
        ("kpoints=[[0.0,0.0,0.0,1.0]]\nkgrid=[1,1,1]", "both"),
        ("kpoints=[[0,0,0]]\nkshift=[0,0,0]", "requires"),
    ] {
        c.text = original.replace("kpoints=[[0.0,0.0,0.0]]", input);
        let out = c.run("", None);
        status(&out, 1);
        assert!(
            String::from_utf8_lossy(&out.stderr).contains(message),
            "{input}: {}",
            String::from_utf8_lossy(&out.stderr)
        );
    }
    let config = Config::parse(&original.replace(
        "kpoints=[[0.0,0.0,0.0]]",
        "kgrid=[2,2,2]\nkshift=[0.5,0.5,0.5]",
    ))
    .unwrap();
    let points = config.basis.resolve_kpoints().unwrap();
    assert_eq!(points.len(), 4);
    assert!(points.iter().all(|k| k.fractional != [0.0; 3]));
}

#[test]
#[ignore = "multiple SCF solves; cargo test --release -p dft-cli --test run kpoint_scf -- --ignored --nocapture"]
fn kpoint_scf() {
    let mut c = Case::new(false);
    c.fixture = c.fixture.replace("z_valence=\"4\"", "z_valence=\"2\"");
    fs::write(c.dir.join("synthetic.xml"), &c.fixture).unwrap();
    let original = c.text.clone();
    let controls = "\n[scf]\nnbands=3\ndensity_residual_tolerance=1e-12\nenergy_tolerance_ha=1e-13\neigen_tolerance=1e-12\nmax_iterations=300\n";
    let forms = [
        "kgrid=[2,1,1]\nkshift=[0.5,0,0]\nktime_reversal=false",
        "kpoints=[[0.25,0,0,0.5],[-0.25,0,0,0.5]]",
        "kgrid=[2,1,1]\nkshift=[0.5,0,0]",
        "kpoints=[[0.25,0,0,0.5],[-0.25,0,0,0.5]]\nktime_reversal=true",
    ];
    let mut energies = Vec::new();
    for (i, form) in forms.iter().enumerate() {
        c.text = original.replace("kpoints=[[0.0,0.0,0.0]]", form);
        let out = c.run(controls, None);
        status(&out, 0);
        let r = c.json("results", "results.json");
        assert_eq!(r["kpoints"]["count"], if i < 2 { 2 } else { 1 });
        assert!(String::from_utf8_lossy(&out.stdout).contains("indirect"));
        let e = r["energies"]["total"]["ha"].as_f64().unwrap();
        println!("k form {i}: E={e:.16e} Ha");
        energies.push(e);
    }
    for e in &energies {
        assert!((e - energies[0]).abs() < 1e-10);
    }
    c.text = original.clone();
    status(&c.run("", None), 0);
    let r = c.json("results", "results.json");
    let mut keys: Vec<_> = r.as_object().unwrap().keys().map(String::as_str).collect();
    keys.sort();
    let mut expected = vec![
        "schema_version",
        "tool",
        "timestamp",
        "input_sha256",
        "species",
        "cell",
        "atoms",
        "cutoff",
        "basis",
        "electrons",
        "bands",
        "occupations",
        "xc",
        "scf",
        "convergence",
        "energies",
        "total_energy",
        "eigenvalues",
        "timings",
        "threading",
        "memory_estimate",
        "units",
    ];
    expected.sort();
    assert_eq!(keys, expected);
    let legacy_energy = r["energies"]["total"]["ha"].as_f64().unwrap();
    for form in [
        "kgrid=[1,1,1]",
        "kpoints=[[0,0,0,1]]",
        "kpoints=[[0,0,0]]\nkgrid=[1,1,1]",
    ] {
        c.text = original.replace("kpoints=[[0.0,0.0,0.0]]", form);
        status(&c.run("", None), 0);
        let gamma = c.json("results", "results.json");
        let energy = gamma["energies"]["total"]["ha"].as_f64().unwrap();
        println!(
            "Gamma legacy/new energy difference: {:.16e} Ha",
            (energy - legacy_energy).abs()
        );
        assert_eq!(energy.to_bits(), legacy_energy.to_bits());
        assert_eq!(gamma["kpoints"]["count"], 1);
    }

    c.text = original
        .replace(
            "kpoints=[[0.0,0.0,0.0]]",
            "kgrid=[4,4,1]\nktime_reversal=false",
        )
        .replace("ecutwfc_ry=3.0", "ecutwfc_ry=5.0");
    let out = c.run("\n[scf]\nnbands=3\nmax_iterations=1", None);
    status(&out, 3);
    let console = String::from_utf8_lossy(&out.stdout);
    assert!(console.contains("remaining k points: see results.json"));
    assert_eq!(
        console
            .lines()
            .filter(|line| line.starts_with('k') && line.contains("weight"))
            .count(),
        6
    );
    assert_eq!(
        c.json("results", "results.json")["kpoints"]["points"]
            .as_array()
            .unwrap()
            .len(),
        16
    );
    fs::write(
        c.dir.join("a.xml"),
        common::gaussian_with(1.0, 0.55, true, 401),
    )
    .unwrap();
    fs::write(
        c.dir.join("b.xml"),
        common::gaussian_with(1.0, 1.15, true, 401),
    )
    .unwrap();
    c.text = original.replace("species='Si'", "species='A'")
        .replace("fractional=[0.13,0.21,0.32]", "fractional=[0.11,0.21,0.17]\n[[system.atoms]]\nspecies='B'\nfractional=[0.39,0.57,0.68]")
        .replace("Si='synthetic.xml'", "A='a.xml'\nB='b.xml'")
        .replace("ecutwfc_ry=3.0", "ecutwfc_ry=5.0")
        .replace("kpoints=[[0.0,0.0,0.0]]", "kgrid=[3,1,1]");
    let out = c.run(controls, None);
    status(&out, 4);
    let message = String::from_utf8_lossy(&out.stderr);
    println!("gap-less sampling: {message}");
    assert!(message.contains("gap"));
}
