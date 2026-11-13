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
// Change only the identity metadata; the synthetic radial arrays remain unchanged.
fn silicon_fixture() -> String {
    common::fixture().replace("element=\"Xx\"", "element=\"Si\"")
}
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
            common::gaussian().replace("element=\"Xx\"", "element=\"Si\"")
        } else {
            silicon_fixture()
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
        // These legacy-output tests pin the pre-ladder path (the C1 cutoff ladder is
        // on by default since 2026-11-10 and has its own tests).
        cmd.env("BAANDI_CUTOFF_LADDER", "off");
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
    // Legacy fields plus the additive preflight report, checked on CLI output.
    // Legacy fields plus additive per-band accuracy, checked on actual output.
    // 24 legacy keys plus gamma_real_used and gamma_real_fallback_reason (default Auto
    // mode, ADR 0040) and atomic_start (default atomic start, ADR 0042).
    assert_eq!(first.as_object().unwrap().len(), 27);
    assert!(first.get("kpoints").is_none());
    for key in [
        "preflight",
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
    // P48F: the closure of the time budget adds up to total_s.
    let closure = &first["timings"]["closure"];
    let parts: f64 = [
        "setup_s",
        "phases_s",
        "certificate_outside_phases_s",
        "step_overhead_s",
        "ladder_coarse_s",
        "remainder_s",
    ]
    .iter()
    .map(|k| closure[*k].as_f64().unwrap())
    .sum();
    let total = first["timings"]["total_s"].as_f64().unwrap();
    assert!(
        (parts - total).abs() <= 1e-9 * total.max(1.0),
        "{parts} {total}"
    );
    assert!(first["timings"]["certificate"]["blocks"].as_u64().unwrap() >= 1);
    assert!(first["timings"]["steps_s"].as_f64().unwrap() > 0.0);
    let occupied = first["occupations"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|f| f.as_f64().unwrap() > 0.0)
        .count();
    let accuracy = first["band_accuracy"].as_array().unwrap();
    assert_eq!(accuracy.len(), first["bands"].as_u64().unwrap() as usize);
    for (j, band) in accuracy.iter().enumerate() {
        // The first empty band (j == occupied) takes the buffer threshold in steps above
        // the tolerance floor (BAANDI_EMPTY_LOOSE, default on); the certificate
        // re-solve checks it strictly.
        if j != occupied {
            assert_eq!(band["buffer"], j > occupied);
        }
        assert!(band["residual_ha"].as_f64().unwrap().is_finite());
        if j <= occupied {
            assert_eq!(band["converged"], true);
        }
    }

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
        1, // Invalid beta now fails preflight instead of diverging in the solver.
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
            silicon_fixture().replace(from, to),
        )
        .unwrap();
        status(&c.run("", None), 1);
    }
    fs::remove_file(c.dir.join("synthetic.xml")).unwrap();
    status(&c.run("", None), 1);
    fs::write(c.dir.join("synthetic.xml"), silicon_fixture()).unwrap();
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
            silicon_fixture()
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
    status(&output, 1);
    let stderr = String::from_utf8_lossy(&output.stderr);
    for message in [
        "xc.functional",
        "SLA PW",
        "lda-pz",
        "select a matching UPF/functional",
    ] {
        assert!(stderr.contains(message), "{stderr}");
    }
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
    assert!(legacy.get("input_expressions").is_none());
    assert_eq!(
        units
            .as_object_mut()
            .unwrap()
            .remove("input_expressions")
            .unwrap(),
        serde_json::json!([{ "key_path": "basis.cutoff", "text": "3 Ry", "value": 3.0 }])
    );
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
        "band_accuracy",
        "preflight",
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
        common::gaussian_with(1.0, 0.55, true, 401).replace("element=\"Xx\"", "element=\"H\""),
    )
    .unwrap();
    fs::write(
        c.dir.join("b.xml"),
        common::gaussian_with(1.0, 1.15, true, 401).replace("element=\"Xx\"", "element=\"H\""),
    )
    .unwrap();
    c.text = original.replace("species='Si'", "species='A'")
        .replace("fractional=[0.13,0.21,0.32]", "fractional=[0.11,0.21,0.17]\n[[system.atoms]]\nspecies='B'\nfractional=[0.39,0.57,0.68]")
        .replace("Si='synthetic.xml'", "A={path='a.xml', element='H'}\nB={path='b.xml', element='H'}")
        .replace("ecutwfc_ry=3.0", "ecutwfc_ry=5.0")
        .replace("kpoints=[[0.0,0.0,0.0]]", "kgrid=[3,1,1]");
    let out = c.run(controls, None);
    status(&out, 4);
    let message = String::from_utf8_lossy(&out.stderr);
    println!("gap-less sampling: {message}");
    assert!(message.contains("gap"));
}

#[test]
fn expressions_preserve_energy_and_record_resolutions() {
    let mut c = Case::new(false);
    c.text = c.text.replace("0.13,0.21,0.32", "0.25,0.21,0.32");
    status(&c.run("", Some(Path::new("plain"))), 0);
    let plain = c.json("plain", "results.json");
    // Switch manifest (work order A1): the whole table, effective value and source.
    let switches = plain["scf"]["switches"].as_object().unwrap();
    assert_eq!(switches.len(), 47);
    let ladder = &switches["BAANDI_CUTOFF_LADDER"];
    assert_eq!(
        (&ladder["value"], &ladder["source"]),
        (&serde_json::json!("off"), &serde_json::json!("env"))
    );
    if std::env::var_os("BAANDI_BUFFER_TOL").is_none() {
        let tol = &switches["BAANDI_BUFFER_TOL"];
        assert_eq!(
            (&tol["value"], &tol["source"]),
            (&serde_json::json!(1e-3), &serde_json::json!("default"))
        );
    }
    let expected_keys = [
        "band_accuracy",
        "preflight",
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
        // Recorded by the default Auto mode (ADR 0040).
        "gamma_real_used",
        "gamma_real_fallback_reason",
        // Recorded by the default atomic start (ADR 0042).
        "atomic_start",
    ];
    let keys: std::collections::BTreeSet<_> = plain
        .as_object()
        .unwrap()
        .keys()
        .map(String::as_str)
        .collect();
    assert_eq!(keys, expected_keys.into_iter().collect());
    c.text = c
        .text
        .replace("4.0", "'8/2'")
        .replace("0.25,0.21,0.32", "'1/4',0.21,0.32")
        .replace("ecutwfc_ry=3.0", "cutoff='0.5*6 Ry'");
    status(&c.run("", Some(Path::new("expr"))), 0);
    let expr = c.json("expr", "results.json");
    for key in [
        "cell",
        "atoms",
        "cutoff",
        "basis",
        "scf",
        "energies",
        "eigenvalues",
    ] {
        assert_eq!(plain[key], expr[key], "{key}");
    }
    let energy = plain["total_energy"]["ha"].as_f64().unwrap();
    let delta = (energy - expr["total_energy"]["ha"].as_f64().unwrap()).abs();
    assert!(delta <= 1e-12);
    let records = expr["input_expressions"].as_array().unwrap();
    assert_eq!(records.len(), 5);
    assert!(records.contains(
        &serde_json::json!({"key_path":"system.atoms[0].fractional[0]", "text":"1/4", "value":0.25})
    ));
    assert!(
        records.contains(
            &serde_json::json!({"key_path":"basis.cutoff", "text":"0.5*6 Ry", "value":3.0})
        )
    );
    let out = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
        .arg("validate")
        .arg(c.dir.join("input.toml"))
        .output()
        .unwrap();
    status(&out, 0);
    let stdout = String::from_utf8(out.stdout).unwrap();
    assert!(stdout.contains("system.atoms[0].fractional[0]: \"1/4\" = 0.25"));
    println!(
        "expression equivalence: total={energy:.17e} Ha, delta={delta:.3e} Ha, records={}",
        records.len()
    );

    c.text = c.text.replace("'1/4'", "'1/0'");
    let invalid = c.run("", Some(Path::new("invalid")));
    status(&invalid, 1);
    let stderr = String::from_utf8(invalid.stderr).unwrap();
    assert!(stderr.contains("system.atoms[0].fractional[0]") && stderr.contains("\"1/0\""));
    let invalid = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
        .arg("validate")
        .arg(c.dir.join("input.toml"))
        .output()
        .unwrap();
    status(&invalid, 1);

    // The quantity expression is resolved identically to its numeric quantity.
    let old = include_str!("../../../examples/si_lda.toml");
    let unit = old.replace("lattice_unit = \"angstrom\"", "").replace("lattice_vectors = [[0.0, 2.715, 2.715], [2.715, 0.0, 2.715], [2.715, 2.715, 0.0]]", "lattice_vectors_with_units = [['0 Bohr', '5.43/2 Angstrom', '5.43/2 Angstrom'], ['5.43/2 Angstrom', '0 Bohr', '5.43/2 Angstrom'], ['5.43/2 Angstrom', '5.43/2 Angstrom', '0 Bohr']]");
    let a = Config::parse(&unit).unwrap();
    let b = Config::parse(&unit.replace("5.43/2", "2.715")).unwrap();
    assert_eq!(a.system.lattice_vectors, b.system.lattice_vectors);
    println!(
        "5.43/2 Angstrom = {:.17e} Bohr",
        a.system.lattice_vectors[0][1]
    );
}

#[test]
fn quantity_expression_energy_equivalence() {
    let mut c = Case::new(false);
    let edge = dft_units::parse_length("2.715 Angstrom").unwrap();
    c.text = c
        .text
        .replace("4.0", &format!("{edge:.17e}"))
        .replace("0.13,0.21,0.32", "0.25,0.21,0.32");
    status(&c.run("", Some(Path::new("numeric"))), 0);
    let numeric = c.json("numeric", "results.json");
    let lattice_line = c
        .text
        .lines()
        .find(|line| line.starts_with("lattice_vectors="))
        .unwrap()
        .to_owned();
    c.text = c.text.replace("lattice_unit='bohr'", "")
        .replace(&lattice_line, "lattice_vectors_with_units=[['5.43/2 Angstrom','0 Bohr','0 Bohr'],['0 Bohr','5.43/2 Angstrom','0 Bohr'],['0 Bohr','0 Bohr','5.43/2 Angstrom']]")
        .replace("0.25,0.21,0.32", "'1/4',0.21,0.32");
    status(&c.run("", Some(Path::new("quantity"))), 0);
    let quantity = c.json("quantity", "results.json");
    assert_eq!(numeric["cell"], quantity["cell"]);
    assert_eq!(numeric["atoms"], quantity["atoms"]);
    let energy = numeric["total_energy"]["ha"].as_f64().unwrap();
    let delta = (energy - quantity["total_energy"]["ha"].as_f64().unwrap()).abs();
    assert!(delta <= 1e-12);
    println!("quantity expression equivalence: total={energy:.17e} Ha, delta={delta:.3e} Ha");
}

#[test]
fn pbe_selection_override_and_spin_rejection() {
    let mut c = Case::new(false);
    c.text = c.text.replace("functional='lda-pz'", "functional='pbe'");
    status(&c.run("", None), 1);
    let output = c.run("allow_functional_override=true\n", None);
    status(&output, 0);
    assert!(String::from_utf8_lossy(&output.stdout).contains("override_applied=true"));
    let r = c.json("results", "results.json");
    assert_eq!(r["xc"], "pbe");
    assert_eq!(r["allow_functional_override"], true);
    assert_eq!(r["functional_override_applied"], true);
    fs::write(
        c.dir.join("synthetic.xml"),
        c.fixture.replace("SLA PZ", "PBE"),
    )
    .unwrap();
    status(&c.run("", None), 0);
    let r = c.json("results", "results.json");
    assert_eq!(r["functional_override_applied"], false);
    assert!(Config::parse(&format!("{}nspin=2\n", c.text)).is_err());
    assert!(Config::parse(&c.text.replace("[system]", "[system]\nnspin=2")).is_err());
    c.text = c.text.replace("functional='pbe'", "functional='lda-pz'");
    status(&c.run("allow_functional_override=true\n", None), 0);
    let r = c.json("results", "results.json");
    assert_eq!(r["functional_override_applied"], true);
    assert_eq!(r["xc"], "lda-pz");
}

#[test]
fn threading_metadata_and_identical_outputs() {
    for nk in [1, 2] {
        let mut c = Case::new(false);
        if nk == 2 {
            c.text = c.text.replace(
                "kpoints=[[0.0,0.0,0.0]]",
                "kpoints=[[0.0,0.0,0.0,0.5],[0.05,0.0,0.0,0.5]]",
            );
        }
        let mut reference = None;
        for requested in [None, Some(0), Some(1), Some(2), Some(3), Some(8), Some(16)] {
            let name = format!("threads-{requested:?}");
            let suffix = requested.map_or(String::new(), |n| format!("\n[scf]\nthreads={n}"));
            let output = c.run(&suffix, Some(Path::new(&name)));
            if requested == Some(0) {
                // Zero is now an invalid request, not silently normalized to one.
                status(&output, 1);
                assert!(String::from_utf8_lossy(&output.stderr).contains("scf.threads"));
                assert!(!c.dir.join(&name).exists());
                continue;
            }
            status(&output, 0);
            let mut json = c.json(&name, "results.json");
            let requested = requested.unwrap_or(1);
            assert_eq!(json["threading"]["threads_requested"], requested);
            assert_eq!(json["threading"]["threads_used"], requested.max(1));
            assert_eq!(json["threading"]["threads"], requested.max(1));
            let outer = requested.max(1).min(nk);
            let inner = requested.max(1);
            assert_eq!(json["threading"]["outer_threads"], outer);
            assert_eq!(json["threading"]["inner_threads"], inner);
            // Outer and inner tasks share one pool, not separate worker teams.
            assert!(outer <= inner);
            assert!(
                json["threading"]["backend"]
                    .as_str()
                    .unwrap()
                    .contains("persistent Rayon pool")
            );
            assert_eq!(json["scf"]["threads"], requested);
            strip_times(&mut json);
            let object = json.as_object_mut().unwrap();
            // The summary describes the requested threads; compare its other
            // preflight fields as well as every numerical result below.
            let summary = object["preflight"]["summary"].as_str().unwrap();
            assert!(summary.contains(&format!("{requested} threads")));
            object["preflight"]
                .as_object_mut()
                .unwrap()
                .remove("summary");
            object.remove("threading");
            object.remove("input_sha256");
            object["scf"].as_object_mut().unwrap().remove("threads");
            let density = fs::read(c.dir.join(&name).join("density.bin")).unwrap();
            if let Some((ref_json, ref_density)) = &reference {
                assert_eq!(ref_json, &json);
                assert_eq!(ref_density, &density);
            } else {
                reference = Some((json, density));
            }
        }
    }
}

fn baandi() -> Command {
    let mut command = Command::new(concat!(env!("CARGO_MANIFEST_DIR"), "/src/baandi"));
    let bin_dir = Path::new(env!("CARGO_BIN_EXE_dft-cli")).parent().unwrap();
    let mut paths = vec![bin_dir.to_path_buf()];
    if let Some(path) = std::env::var_os("PATH") {
        paths.extend(std::env::split_paths(&path));
    }
    command.env("PATH", std::env::join_paths(paths).unwrap());
    // Legacy-output comparison across CLI modes: pin the pre-ladder path (the
    // observe/agent modes do not run the C1 ladder).
    command.env("BAANDI_CUTOFF_LADDER", "off");
    command
}

/// CLI modes reuse the same retained-result SCF: exact density bytes and all
/// numerical JSON values must agree (zero tolerance, energies in Ha).
#[test]
fn execution_modes_are_bit_identical_and_alias_preserves_classic() {
    compare_execution_modes(false);
}

#[test]
#[ignore = "requires filesystem Unix socket binding; cargo test -p dft-cli --test run agent_socket_mode -- --ignored --nocapture"]
fn agent_socket_mode() {
    compare_execution_modes(true);
}

fn compare_execution_modes(include_agent: bool) {
    let c = Case::new(false);
    let suffix = "[scf]\nmax_iterations=3\n";
    let classic = c.run(suffix, Some(&c.dir.join("classic")));
    status(&classic, 3);
    let mut expected = c.json("classic", "results.json");
    strip_times(&mut expected);
    for (name, flags) in [
        ("alias", vec![]),
        ("observed", vec!["--observe", "--events", "-"]),
        ("agent", vec!["--agent"]),
    ] {
        if name == "agent" && !include_agent {
            continue;
        }
        let out = c.dir.join(name);
        fs::create_dir(&out).unwrap();
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            fs::set_permissions(&out, fs::Permissions::from_mode(0o700)).unwrap();
        }
        let result = baandi()
            .current_dir(&c.dir)
            .args(["run", "input.toml", "--output-dir", name])
            .args(flags)
            .output()
            .unwrap();
        status(&result, 3);
        assert_eq!(
            hash(&out.join("density.bin")),
            hash(&c.dir.join("classic/density.bin"))
        );
        let mut actual = c.json(name, "results.json");
        strip_times(&mut actual);
        assert_eq!(expected, actual, "mode {name}");
        if name == "alias" {
            let stable_console = |bytes: &[u8]| {
                String::from_utf8_lossy(bytes)
                    .lines()
                    .filter(|line| !line.starts_with("Total wall:"))
                    .map(|line| {
                        if line
                            .split_whitespace()
                            .next()
                            .is_some_and(|first| first.parse::<usize>().is_ok())
                        {
                            let (values, _) = line.rsplit_once(' ').unwrap();
                            values.to_owned()
                        } else {
                            line.to_owned()
                        }
                    })
                    .collect::<Vec<_>>()
            };
            assert_eq!(
                stable_console(&classic.stdout),
                stable_console(&result.stdout)
            );
            assert!(!out.join("events.jsonl").exists());
        } else {
            let log = fs::read_to_string(out.join("events.jsonl")).unwrap();
            let events: Vec<Value> = log
                .lines()
                .map(|l| serde_json::from_str(l).unwrap())
                .collect();
            assert_eq!(events.first().unwrap()["type"], "initialized");
            assert_eq!(events.last().unwrap()["type"], "not_converged");
            assert_eq!(
                events
                    .iter()
                    .filter(|v| v["type"] == "iteration_finished")
                    .count(),
                3
            );
            assert!(String::from_utf8_lossy(&result.stderr).contains("solver=iteration_completed"));
            if name == "observed" {
                assert_eq!(String::from_utf8(result.stdout).unwrap(), log);
            }
            let status = baandi().arg("status").arg(&out).output().unwrap();
            assert!(status.status.success());
            let status: Value = serde_json::from_slice(&status.stdout).unwrap();
            assert_eq!(status["state"], "not_converged");
        }
    }
    println!(
        "P11H synthetic fixture: Classic/alias/Observed (agent={include_agent}), 3 iterations; density bytes and numerical JSON identical (0 Ha energy difference)."
    );
}

#[test]
fn explain_and_status_report_evidence_and_reject_invalid_input() {
    let c = Case::new(false);
    fs::write(c.dir.join("input.toml"), &c.text).unwrap();
    let explain = baandi()
        .current_dir(&c.dir)
        .args(["explain", "input.toml"])
        .output()
        .unwrap();
    assert!(
        explain.status.success(),
        "{}",
        String::from_utf8_lossy(&explain.stderr)
    );
    assert!(
        String::from_utf8_lossy(&explain.stdout)
            .contains("does not establish numerical convergence")
    );
    assert!(!c.dir.join("results").exists());
    fs::write(
        c.dir.join("input.toml"),
        format!("{}\n[scf]\nmixing_beta=2.0\n", c.text),
    )
    .unwrap();
    let invalid = baandi()
        .current_dir(&c.dir)
        .args(["explain", "input.toml"])
        .output()
        .unwrap();
    assert!(!invalid.status.success());
    let message = String::from_utf8_lossy(&invalid.stderr);
    assert!(
        message.contains("scf.mixing_beta") && message.contains("(0, 1]"),
        "{message}"
    );
    let path = c.dir.join("events.jsonl");
    let event = r#"{"schema_version":"1.0","calc_id":"test","sequence":1,"iteration":1,"state_version":1,"type":"iteration_finished","data":{}}"#;
    fs::write(&path, format!("{event}\n{{\"partial\":")).unwrap();
    let out = baandi().arg("status").arg(&path).output().unwrap();
    assert!(out.status.success());
    let v: Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(v["partial_tail_ignored"], true);
    assert_eq!(v["state"], "nonterminal: process liveness unknown");
    fs::write(&path, format!("{event}\n{event}\n")).unwrap();
    let out = baandi().arg("status").arg(&path).output().unwrap();
    assert!(!out.status.success());
}

#[cfg(unix)]
#[test]
#[ignore = "requires filesystem Unix socket binding; cargo test -p dft-cli --test run observed_socket -- --ignored --nocapture"]
fn observed_socket_mirror_and_existing_log_protection() {
    use std::os::unix::net::UnixDatagram;
    let c = Case::new(false);
    fs::write(
        c.dir.join("input.toml"),
        format!("{}\n[scf]\nmax_iterations=1\n", c.text),
    )
    .unwrap();
    // Darwin socket paths are limited to 104 bytes; /tmp avoids long temp roots.
    let socket_path = PathBuf::from(format!(
        "/tmp/p11h-{}-{}.sock",
        std::process::id(),
        SERIAL.fetch_add(1, Ordering::SeqCst)
    ));
    let socket = UnixDatagram::bind(&socket_path).unwrap();
    socket
        .set_read_timeout(Some(std::time::Duration::from_secs(5)))
        .unwrap();
    let receiver = std::thread::spawn(move || {
        let mut events = Vec::new();
        let mut bytes = [0; 8192];
        loop {
            let n = socket.recv(&mut bytes).unwrap();
            let v: Value = serde_json::from_slice(&bytes[..n]).unwrap();
            let terminal = v["type"] == "not_converged";
            events.push(v);
            if terminal {
                return events;
            }
        }
    });
    let mirror = format!("unix:{}", socket_path.display());
    let mut command = baandi();
    command
        .current_dir(&c.dir)
        .args(["run", "input.toml", "--observe", "--events", &mirror]);
    let first = command.output().unwrap();
    status(&first, 3);
    let events = receiver.join().unwrap();
    fs::remove_file(socket_path).unwrap();
    assert_eq!(events.len(), 3);
    let before = fs::read(c.dir.join("results/events.jsonl")).unwrap();
    let repeat = command.output().unwrap();
    status(&repeat, 5);
    assert_eq!(
        before,
        fs::read(c.dir.join("results/events.jsonl")).unwrap()
    );
}

#[test]
fn observed_file_mirror_is_exact_and_logs_are_never_overwritten() {
    let c = Case::new(false);
    fs::write(
        c.dir.join("input.toml"),
        format!("{}\n[scf]\nmax_iterations=1\n", c.text),
    )
    .unwrap();
    let mut command = baandi();
    command.current_dir(&c.dir).args([
        "run",
        "input.toml",
        "--observe",
        "--events",
        "mirror.jsonl",
    ]);
    let first = command.output().unwrap();
    status(&first, 3);
    let before = fs::read(c.dir.join("results/events.jsonl")).unwrap();
    assert_eq!(before, fs::read(c.dir.join("mirror.jsonl")).unwrap());
    let repeat = command.output().unwrap();
    status(&repeat, 5);
    assert_eq!(
        before,
        fs::read(c.dir.join("results/events.jsonl")).unwrap()
    );
    assert_eq!(before, fs::read(c.dir.join("mirror.jsonl")).unwrap());
}

#[test]
fn final_spectral_guard_is_additive_and_audits_only_reported_state() {
    let c = Case::new(false);
    // Pinned to the complex path: the guard forces it, and bit-identical energies are the claim.
    status(
        &c.run("[scf]\ngamma_real='off'\n", Some(Path::new("off"))),
        0,
    );
    let off = c.json("off", "results.json");
    assert_eq!(off["scf"]["spectral_guard"], "off");
    assert!(off.get("spectral_completeness").is_none());
    let run = c.run(
        "[scf]\nspectral_guard='final'\ngamma_real='off'\n",
        Some(Path::new("final")),
    );
    status(&run, 0);
    let guarded = c.json("final", "results.json");
    assert_eq!(guarded["schema_version"], off["schema_version"]);
    assert_eq!(guarded["total_energy"], off["total_energy"]);
    assert_eq!(guarded["eigenvalues"], off["eigenvalues"]);
    assert_eq!(
        hash(&c.dir.join("off/density.bin")),
        hash(&c.dir.join("final/density.bin"))
    );
    let report = &guarded["spectral_completeness"];
    assert_eq!(report["passed"], true);
    assert_eq!(report["retries"], 0);
    let events = report["events"].as_array().unwrap();
    assert_eq!(events.len(), 1);
    assert_eq!(events[0]["iteration"], guarded["convergence"]["iterations"]);
    assert_eq!(events[0]["k_point"], 0);
    assert_eq!(events[0]["method"], "full_space_dense");
    assert_eq!(
        events[0]["extra_operator_applications"],
        report["extra_operator_applications"]
    );
    assert!(report["extra_operator_applications"].as_u64().unwrap() > 0);
    assert_eq!(
        guarded.as_object().unwrap().len(),
        off.as_object().unwrap().len() + 1
    );
    let bad = c.run(
        "[scf]\nspectral_guard='sometimes'\n",
        Some(Path::new("bad")),
    );
    status(&bad, 1);
    assert!(String::from_utf8_lossy(&bad.stderr).contains("scf.spectral_guard"));
}
