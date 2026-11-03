//! Process-boundary preflight checks using generated synthetic silicon only.
#[allow(dead_code)]
mod common;
use std::{
    fs,
    path::PathBuf,
    process::{Command, Output},
    sync::atomic::{AtomicUsize, Ordering},
};
static NEXT: AtomicUsize = AtomicUsize::new(0);
struct Case {
    dir: PathBuf,
}
impl Case {
    fn new() -> Self {
        let dir = std::env::temp_dir().join(format!(
            "preflight-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&dir).unwrap();
        fs::write(
            dir.join("synthetic.xml"),
            common::fixture()
                .replace("element=\"Xx\"", "element=\"Si\"")
                .replace("wfc_cutoff=\"0\"", "wfc_cutoff=\"40\""),
        )
        .unwrap();
        Self { dir }
    }
    fn text(&self, label: &str, entry: &str, scf: &str) -> String {
        format!(
            "[system]\nname='synthetic'\nlattice_unit='bohr'\nlattice_vectors=[[4.0,0.0,0.0],[0.0,4.0,0.0],[0.0,0.0,4.0]]\n[[system.atoms]]\nspecies='{label}'\nfractional=[0.0,0.0,0.0]\n[pseudopotentials]\n{label}={entry}\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[scf]\n{scf}\n"
        )
    }
    fn command(&self, text: &str, command: &str, json: bool) -> Output {
        fs::write(self.dir.join("input.toml"), text).unwrap();
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_dft-cli"));
        cmd.current_dir(&self.dir).args([command, "input.toml"]);
        if json {
            cmd.arg("--json");
        }
        cmd.output().unwrap()
    }
}
impl Drop for Case {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.dir).unwrap();
    }
}
fn reject(out: Output, keys: &[&str]) {
    assert_eq!(out.status.code(), Some(1));
    let err = String::from_utf8(out.stderr).unwrap();
    for key in keys {
        assert!(err.contains(key), "missing {key}: {err}");
    }
}
#[test]
fn identity_and_schema() {
    let c = Case::new();
    let wrong = c.text("Fe", "'synthetic.xml'", "");
    for command in ["validate", "run"] {
        reject(
            c.command(&wrong, command, false),
            &[
                "pseudopotentials.Fe",
                "Fe",
                "Si",
                "synthetic.xml",
                "correct",
            ],
        );
    }
    reject(
        c.command(&c.text("Si_up", "'synthetic.xml'", ""), "validate", false),
        &["pseudopotentials.Si_up.element", "explicit"],
    );
    reject(
        c.command(
            &c.text("Si_up", "{path='synthetic.xml', element='Fe'}", ""),
            "validate",
            false,
        ),
        &["Si_up", "Fe", "Si"],
    );
    let good = c.text("Si_up", "{path='synthetic.xml', element='Si'}", "");
    let out = c.command(&good, "validate", true);
    assert_eq!(out.status.code(), Some(0), "{:?}", out);
    let report: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(report["schema_version"], 1);
    assert_eq!(report["valid"], true);
    assert_eq!(report["electrons"], 4);
    assert_eq!(report["bands"], 6);
    assert_eq!(report["warnings"].as_array().unwrap().len(), 1);
    assert!(report["warnings"][0].as_str().unwrap().contains("40 Ry"));
    assert!(report["overrides"].as_array().unwrap().is_empty());
    assert_eq!(report["fft_grid"].as_array().unwrap().len(), 3);
    assert!(report["estimated_array_bytes"].as_u64().unwrap() > 0);
    let summary = report["summary"].as_str().unwrap();
    for part in [
        "Si_up:1",
        "1 atoms",
        "4 electrons",
        "lda-pz",
        "k mesh",
        "3 Ry",
        "6 bands",
        "1 threads",
    ] {
        assert!(summary.contains(part), "{summary}");
    }
    let out = c.command(&good, "validate", false);
    assert!(String::from_utf8_lossy(&out.stderr).contains("warning:"));
}
#[test]
fn scf_controls_and_electrons() {
    let c = Case::new();
    for setting in [
        "max_iterations=0",
        "max_iterations=-1",
        "threads=0",
        "energy_tolerance_ha=0.0",
        "density_residual_tolerance=-1.0",
        "eigen_tolerance=nan",
        "mixing_beta=0.0",
        "mixing_beta=1.1",
        "mixing='unknown'",
        "nbands=1",
    ] {
        let key = format!("scf.{}", setting.split('=').next().unwrap());
        reject(
            c.command(&c.text("Si", "'synthetic.xml'", setting), "validate", false),
            &[&key],
        );
    }
    let text = c.text("Si", "'synthetic.xml'", "");
    reject(
        c.command(
            &text.replace("ecutwfc_ry=3.0", "ecutwfc_ry=0.001"),
            "validate",
            false,
        ),
        &["basis.ecutwfc_ry"],
    );
    reject(
        c.command(
            &text.replace("kpoints=[[0.0,0.0,0.0]]", "kgrid=[0,1,1]"),
            "validate",
            false,
        ),
        &["basis"],
    );
    let path = c.dir.join("synthetic.xml");
    fs::write(
        &path,
        fs::read_to_string(&path)
            .unwrap()
            .replace("z_valence=\"4\"", "z_valence=\"3\""),
    )
    .unwrap();
    reject(
        c.command(&text, "validate", false),
        &["odd electron count 3", "system.atoms"],
    );
}
#[test]
fn override_is_printed_and_persisted() {
    let c = Case::new();
    let text = c.text(
        "Fe",
        "{path='synthetic.xml', allow_element_mismatch=true}",
        "max_iterations=1",
    );
    let out = c.command(&text, "validate", true);
    assert_eq!(out.status.code(), Some(0));
    let report: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(report["overrides"].as_array().unwrap().len(), 1);
    let out = c.command(&text, "run", false);
    assert_eq!(out.status.code(), Some(3), "{:?}", out);
    assert!(String::from_utf8_lossy(&out.stderr).contains("allow_element_mismatch = true"));
    let result: serde_json::Value =
        serde_json::from_slice(&fs::read(c.dir.join("results/results.json")).unwrap()).unwrap();
    assert_eq!(result["preflight"], report);
}

#[test]
fn functional_files_and_kpoint_gap() {
    let c = Case::new();
    let text = c.text("Si", "'synthetic.xml'", "extra_bands=0");
    assert_eq!(c.command(&text, "validate", false).status.code(), Some(0));
    reject(
        c.command(
            &text.replace("kpoints=[[0.0,0.0,0.0]]", "kgrid=[1,1,1]"),
            "validate",
            false,
        ),
        &["scf.extra_bands", "empty band"],
    );
    let text = c.text("Si", "'missing.xml'", "");
    reject(
        c.command(&text, "validate", false),
        &["pseudopotentials.Si", "missing.xml", "readable"],
    );
    let text = c.text("Si", "'synthetic.xml'", "");
    let wrong = text.replace("lda-pz", "pbe");
    reject(
        c.command(&wrong, "validate", false),
        &["xc.functional", "SLA PZ", "allow_functional_override"],
    );
    let out = c.command(
        &wrong.replace(
            "functional='pbe'",
            "functional='pbe'\nallow_functional_override=true",
        ),
        "validate",
        true,
    );
    assert_eq!(out.status.code(), Some(0));
    let report: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert!(
        report["overrides"][0]
            .as_str()
            .unwrap()
            .contains("xc.allow_functional_override = true")
    );
    fs::write(c.dir.join("synthetic.xml"), "<UPF>").unwrap();
    reject(
        c.command(&text, "validate", false),
        &["pseudopotentials.Si", "synthetic.xml"],
    );
}

#[test]
fn valid_preflight_preserves_engine_energy_and_density() {
    use dft_cli::{config::Config, run::options};
    use dft_engine::{ScfProblem, SpeciesPotential, input::Atom};
    let c = Case::new();
    let text = c.text("Si", "'synthetic.xml'", "max_iterations=1");
    let config = Config::parse(&text).unwrap();
    let potential = dft_upf::from_file(c.dir.join("synthetic.xml")).unwrap();
    let direct = ScfProblem::new(
        &config.cell().unwrap(),
        &[Atom::new("Si", [0.; 3])],
        &[SpeciesPotential {
            label: "Si",
            potential: &potential,
        }],
        options(&config).unwrap(),
    )
    .unwrap()
    .run_outcome()
    .unwrap()
    .result;
    let out = c.command(&text, "run", false);
    assert_eq!(out.status.code(), Some(3));
    let result: serde_json::Value =
        serde_json::from_slice(&fs::read(c.dir.join("results/results.json")).unwrap()).unwrap();
    let energy = result["energies"]["total"]["ha"].as_f64().unwrap();
    assert!((energy - direct.total_hartree()).abs() < 1e-12);
    let expected: Vec<_> = direct
        .density
        .iter()
        .flat_map(|x| x.to_le_bytes())
        .collect();
    assert_eq!(
        fs::read(c.dir.join("results/density.bin")).unwrap(),
        expected
    );
    println!(
        "Synthetic 1-step preflight/direct engine energy difference: {} Ha; density bit-identical",
        energy - direct.total_hartree()
    );
}

#[test]
fn example_templates_validate_with_available_synthetic_si() {
    let c = Case::new();
    for template in [
        include_str!("../../../examples/si_lda.toml"),
        include_str!("../../../examples/si_lda_expr.toml"),
        include_str!("../../../examples/si_lda_kgrid.toml"),
    ] {
        let text = template.replace("/absolute/path/to/Si.UPF", "synthetic.xml");
        let out = c.command(&text, "validate", true);
        assert_eq!(out.status.code(), Some(0), "{:?}", out);
        let report: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
        assert_eq!(report["electrons"], 8);
        assert_eq!(report["bands"], 8);
    }
}
