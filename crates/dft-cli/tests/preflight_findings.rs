//! Stable additive findings from existing rules; synthetic UPF only.
#[allow(dead_code)]
mod common;
use dft_cli::config::{Config, preflight};
use serde_json::{Value, json};
use std::{
    fs,
    path::PathBuf,
    process::Command,
    sync::atomic::{AtomicUsize, Ordering},
};

static NEXT: AtomicUsize = AtomicUsize::new(0);
struct Fixture(PathBuf);
impl Fixture {
    fn new(suggestion: f64) -> Self {
        let path = std::env::temp_dir().join(format!(
            "findings-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&path).unwrap();
        fs::write(
            path.join("synthetic.xml"),
            common::fixture()
                .replace("element=\"Xx\"", "element=\"Si\"")
                .replace("wfc_cutoff=\"0\"", &format!("wfc_cutoff=\"{suggestion}\"")),
        )
        .unwrap();
        Self(path)
    }
    fn input(&self, element: &str, allow: bool, functional: &str, override_xc: bool) -> String {
        format!(
            r#"[system]
name='synthetic'
lattice_unit='bohr'
lattice_vectors=[[4.0,0.0,0.0],[0.0,4.0,0.0],[0.0,0.0,4.0]]
[[system.atoms]]
species='site.1'
fractional=[0.0,0.0,0.0]
[pseudopotentials]
'site.1'={{path='synthetic.xml', element='{element}', allow_element_mismatch={allow}}}
[basis]
ecutwfc_ry=3.0
kpoints=[[0.0,0.0,0.0]]
[xc]
functional='{functional}'
allow_functional_override={override_xc}
"#
        )
    }
    fn report(&self, text: &str) -> Value {
        let config = Config::parse(text).unwrap();
        serde_json::to_value(preflight::check(&config, &self.0.join("input.toml")).unwrap())
            .unwrap()
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.0).unwrap();
    }
}

#[test]
fn cutoff_code_boundaries_schema_and_determinism() {
    for (suggestion, count) in [(0., 0), (2., 0), (3., 0), (40., 1)] {
        let f = Fixture::new(suggestion);
        let text = f.input("Si", false, "lda-pz", false);
        let report = f.report(&text);
        assert_eq!(report, f.report(&text));
        assert_eq!(report["schema_version"], 1);
        assert_eq!(report["findings_schema_version"], 1);
        assert_eq!(report["valid"], true);
        assert_eq!(report["findings"].as_array().unwrap().len(), count);
        if count != 0 {
            assert_eq!(
                report["findings"][0],
                json!({
                    "code": "PF_CUTOFF_BELOW_SUGGESTED",
                    "severity": "warning",
                    "key_path": "basis.ecutwfc_ry",
                    "evidence": {"species": "site.1", "cutoff_ry": 3.0, "suggested_cutoff_ry": 40.0},
                    "fix_hint": "Increase cutoff and test convergence; the UPF suggestion is not a convergence certificate."
                })
            );
            let typed: preflight::Finding =
                serde_json::from_value(report["findings"][0].clone()).unwrap();
            assert_eq!(serde_json::to_value(typed).unwrap(), report["findings"][0]);
        }
        assert_eq!(
            report["structure"],
            json!({
                "lattice_vectors_bohr": [[4.,0.,0.],[0.,4.,0.],[0.,0.,4.]],
                "atoms": [{"species":"site.1", "fractional":[0.,0.,0.]}]
            })
        );
    }
}

#[test]
fn override_codes_have_positive_and_negative_fixtures() {
    let f = Fixture::new(0.);
    assert_eq!(
        f.report(&f.input("Si", false, "lda-pz", false))["findings"],
        json!([])
    );
    for (element, mismatch) in [("Fe", true), ("Si", false)] {
        let report = f.report(&f.input(element, true, "lda-pz", false));
        assert_eq!(report["findings"][0]["code"], "PF_ELEMENT_OVERRIDE");
        assert_eq!(report["findings"][0]["severity"], "info");
        assert_eq!(
            report["findings"][0]["key_path"],
            "pseudopotentials.\"site.1\".allow_element_mismatch"
        );
        assert_eq!(
            report["findings"][0]["evidence"],
            json!({
                "species":"site.1", "declared_element":element, "upf_element":"Si", "mismatch":mismatch
            })
        );
    }
    let report = f.report(&f.input("Si", false, "pbe", true));
    assert_eq!(report["findings"][0]["code"], "PF_FUNCTIONAL_OVERRIDE");
    assert_eq!(report["findings"][0]["severity"], "info");
    assert_eq!(
        report["findings"][0]["key_path"],
        "xc.allow_functional_override"
    );
    assert_eq!(
        report["findings"][0]["evidence"],
        json!({
            "species":"site.1", "requested_functional":"pbe", "upf_functional":"SLA PZ"
        })
    );
    assert_eq!(
        f.report(&f.input("Si", false, "lda-pz", true))["findings"],
        json!([])
    );
    for text in [
        f.input("Fe", false, "lda-pz", false),
        f.input("Si", false, "pbe", false),
    ] {
        assert!(preflight::check(&Config::parse(&text).unwrap(), &f.0.join("input.toml")).is_err());
    }
}

#[test]
fn external_client_receives_findings_without_agent_or_input_mutation() {
    let f = Fixture::new(40.);
    let text = f.input("Si", false, "lda-pz", false);
    let path = f.0.join("input.toml");
    fs::write(&path, &text).unwrap();
    let run = |json: bool| {
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_dft-cli"));
        cmd.args(["validate", path.to_str().unwrap()]);
        if json {
            cmd.arg("--json");
        }
        cmd.output().unwrap()
    };
    let out = run(true);
    assert_eq!(out.status.code(), Some(0));
    let report: Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(report, f.report(&text));
    let plain = run(false);
    assert_eq!(plain.status.code(), Some(0));
    assert!(String::from_utf8_lossy(&plain.stderr).contains("40 Ry"));
    assert_eq!(fs::read_to_string(path).unwrap(), text);
    assert!(!f.0.join("results").exists());
}

#[test]
fn structure_unit_conversion_preserves_input_order_and_skew() {
    let f = Fixture::new(0.);
    let text = f
        .input("Si", false, "lda-pz", false)
        .replace("lattice_unit='bohr'", "lattice_unit='angstrom'")
        .replace("[4.0,0.0,0.0]", "[4.0,1.0,0.0]");
    let report = f.report(&text);
    let factor = dft_units::Bohr::from(dft_units::Angstrom(1.)).0;
    assert_eq!(
        report["structure"]["lattice_vectors_bohr"],
        json!([
            [4. * factor, factor, 0.],
            [0., 4. * factor, 0.],
            [0., 0., 4. * factor]
        ])
    );
}
