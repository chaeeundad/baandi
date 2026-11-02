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
