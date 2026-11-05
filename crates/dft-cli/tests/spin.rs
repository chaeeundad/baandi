#[allow(dead_code)]
mod common;
use std::{
    fs,
    path::PathBuf,
    process::Command,
    sync::atomic::{AtomicUsize, Ordering},
};
static NEXT: AtomicUsize = AtomicUsize::new(0);
struct Case {
    dir: PathBuf,
}
impl Case {
    fn new() -> Self {
        let dir = std::env::temp_dir().join(format!(
            "mf03-cli-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&dir).unwrap();
        fs::write(
            dir.join("syn.xml"),
            common::gaussian_with(1.0, 1.0, false, 401)
                .replace("SLA PZ", "SLA PW")
                .replace("element=\"Xx\"", "element=\"H\""),
        )
        .unwrap();
        Self { dir }
    }
    fn input(&self, spin: &str) -> String {
        format!(
            "[system]\nname='synthetic spin'\nlattice_unit='bohr'\nlattice_vectors=[[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0]]\n[[system.atoms]]\nspecies='H'\nfractional=[0.0,0.0,0.0]\n[pseudopotentials]\nH='syn.xml'\n[basis]\necutwfc_ry=5.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pw'\n[scf]\nnbands=4\neigen_tolerance=1e-10\ndensity_residual_tolerance=1e-10\nenergy_tolerance_ha=1e-11\n[spin]\n{spin}\n"
        )
    }
    fn run(&self, text: &str, command: &str) -> std::process::Output {
        fs::write(self.dir.join("input.toml"), text).unwrap();
        Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .current_dir(&self.dir)
            .args([command, "input.toml"])
            .output()
            .unwrap()
    }
}
impl Drop for Case {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.dir);
    }
}
#[test]
fn spin_preflight_rejections() {
    let c = Case::new();
    for (spin, key) in [
        ("mode='collinear'", "spin.constraint"),
        ("mode='collinear'\nn_up=1", "spin.constraint"),
        (
            "mode='collinear'\nn_up=1\nn_down=0\ntotal_moment=1.0",
            "spin.constraint",
        ),
        ("mode='collinear'\ntotal_moment=0.0", "spin.constraint"),
        ("mode='collinear'\ntotal_moment=3.0", "spin.constraint"),
        ("mode='collinear'\nn_up=1\nn_down=0\naxis='q'", "spin.axis"),
        (
            "mode='collinear'\nn_up=1\nn_down=0\naxis_cartesian=[0.0,0.0,0.0]",
            "spin.axis",
        ),
        ("mode='unpolarized'\nn_up=1", "spin.mode"),
    ] {
        for command in ["validate", "run"] {
            let out = c.run(&c.input(spin), command);
            assert_eq!(out.status.code(), Some(1));
            assert!(
                String::from_utf8_lossy(&out.stderr).contains(key),
                "{:?}",
                out
            );
        }
    }
    let base = c.input("mode='collinear'\nn_up=1\nn_down=0");
    let cases = [
        (
            base.replace(
                "functional='lda-pw'",
                "functional='lda-pz'\nallow_functional_override=true",
            ),
            "xc.spin: required polarized capability unavailable for lda-pz; select pbe or lda-pw",
        ),
        (
            base.replace("[basis]", "[basis]\nktime_reversal=true"),
            "symmetry.operations",
        ),
        (base.replace("nbands=4", "nbands=1"), "empty band"),
    ];
    for (text, key) in cases {
        for command in ["validate", "run"] {
            let out = c.run(&text, command);
            assert_eq!(out.status.code(), Some(1));
            assert!(
                String::from_utf8_lossy(&out.stderr).contains(key),
                "{:?}",
                out
            );
        }
    }
}
#[test]
fn collinear_cli_certificate_schemas_and_cubes() {
    let c = Case::new();
    let text = c.input("mode='collinear'\ntotal_moment=1.0\naxis='z'\ninitial_moment_seed=0.5");
    let out = c.run(&text, "run");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let path = c.dir.join("results/results.json");
    let a: serde_json::Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
    assert_eq!(a["schema_version"], 2);
    assert_eq!(a["convergence"]["certificate"]["passed"], true);
    assert_eq!(
        a["spin"]["channel_capacity_electrons"],
        serde_json::json!([1, 1])
    );
    let status = Command::new("python3")
        .args(["-I"])
        .arg(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/validate_spin_results.py"))
        .arg(&path)
        .status()
        .unwrap();
    assert!(status.success());
    let out = c.run(
        &text.replace("initial_moment_seed=0.5", "initial_moment_seed=-0.5"),
        "run",
    );
    assert!(out.status.success());
    let b: serde_json::Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
    assert_eq!(
        a["state_contracts"]["model_digest"],
        b["state_contracts"]["model_digest"]
    );
    assert_eq!(
        a["state_contracts"]["state_digest"],
        b["state_contracts"]["state_digest"]
    );
    assert_ne!(a["provenance"], b["provenance"]);
    let out = c.run(
        &text.replace("total_moment=1.0", "total_moment=-1.0"),
        "run",
    );
    assert!(out.status.success());
    let b: serde_json::Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
    assert_ne!(
        a["state_contracts"]["model_digest"],
        b["state_contracts"]["model_digest"]
    );
    assert_ne!(
        a["state_contracts"]["state_digest"],
        b["state_contracts"]["state_digest"]
    );
    assert_eq!(a["total_energy"], b["total_energy"]);
}

#[test]
fn final_guard_reports_both_spin_channels_without_changing_schema_version() {
    let mut reference = None;
    for threads in [1, 2, 4, 8] {
        let c = Case::new();
        let text = c
            .input("mode='collinear'\ntotal_moment=1.0\naxis='z'\ninitial_moment_seed=0.5")
            .replace(
                "[scf]",
                &format!("[scf]\nspectral_guard='final'\nthreads={threads}"),
            );
        let out = c.run(&text, "run");
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        let result: serde_json::Value =
            serde_json::from_slice(&fs::read(c.dir.join("results/results.json")).unwrap()).unwrap();
        assert_eq!(result["schema_version"], 2);
        let report = &result["spectral_completeness"];
        assert_eq!(report["passed"], true);
        assert_eq!(report["retries"], 0);
        let events = report["events"].as_array().unwrap();
        assert_eq!(events.len(), 2);
        for (s, e) in events.iter().enumerate() {
            assert_eq!(e["spin"], s);
            assert_eq!(e["k_point"], 0);
            assert_eq!(e["iteration"], result["convergence"]["iterations"]);
            assert_eq!(e["passed"], true);
        }
        let mut audit = report.clone();
        audit.as_object_mut().unwrap().remove("timings");
        let scientific = serde_json::json!([
            result["total_energy"],
            result["band_data"],
            result["density"],
            audit
        ]);
        if let Some(previous) = &reference {
            assert_eq!(previous, &scientific);
        } else {
            reference = Some(scientific);
        }
    }
}

#[test]
fn spin_force_extension_v2_validates() {
    let c = Case::new();
    let text = c.input("mode='collinear'\nn_up=1\nn_down=0") + "\n[output]\nforces=true\n";
    let out = c.run(&text, "run");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let path = c.dir.join("results/results.json");
    let r: serde_json::Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
    assert_eq!(r["schema_version"], 2);
    assert_eq!(r["state_contracts"]["extension_version"], 2);
    assert!(
        Command::new("python3")
            .arg("-I")
            .arg(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/validate_force_schema.py"))
            .arg(path)
            .status()
            .unwrap()
            .success()
    );
}
