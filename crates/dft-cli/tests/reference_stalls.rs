//! Observational G4G records: no tolerance, stop-rule or production-default changes.
use serde_json::Value;
use std::{fs, path::PathBuf, process::Command};

struct Fixture(PathBuf);
impl Fixture {
    fn new() -> Self {
        let dir = std::env::temp_dir().join(format!("g4g-cli-records-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        // Smaller-box version of the G4G Gaussian atom, generated without
        // committing UPF data. beta/history below are explicit test inputs only.
        let generated = Command::new("python3")
            .args(["-I", "-c", include_str!("g4g_fixture.py")])
            .output()
            .unwrap();
        assert!(
            generated.status.success(),
            "{}",
            String::from_utf8_lossy(&generated.stderr)
        );
        fs::write(dir.join("syn.xml"), generated.stdout).unwrap();
        Self(dir)
    }
    fn run(&self, max_iterations: usize) -> (Value, String) {
        fs::write(
            self.0.join("input.toml"),
            format!(
                r#"
[system]
name='G4G synthetic diagnostic fixture'
lattice_unit='bohr'
lattice_vectors=[[8.0,0.0,0.0],[0.0,8.0,0.0],[0.0,0.0,8.0]]
[[system.atoms]]
species='O'
fractional=[0.13,0.21,0.32]
[pseudopotentials]
O='syn.xml'
[basis]
ecutwfc_ry=12.0
kpoints=[[0.0,0.0,0.0]]
[xc]
functional='pbe'
[scf]
threads=1
nbands=6
max_iterations={max_iterations}
mixing_beta=0.6
mixing_history=15
[spin]
mode='collinear'
n_up=4
n_down=2
[performance]
reference_strict='off'
# The stall is a rounding-dependent trajectory (see the test); GK10 border
# Gram rounding moves this fixture past it on macOS, so pin the full Gram.
cdav_border_gram='off'
"#
            ),
        )
        .unwrap();
        self.execute()
    }

    fn execute(&self) -> (Value, String) {
        let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .current_dir(&self.0)
            .env("BAANDI_REFERENCE_STRICT", "off")
            .args(["run", "input.toml"])
            .output()
            .unwrap();
        let stderr = String::from_utf8(output.stderr).unwrap();
        assert_eq!(output.status.code(), Some(3), "{stderr}");
        let result =
            serde_json::from_slice(&fs::read(self.0.join("results/results.json")).unwrap())
                .unwrap();
        (result, stderr)
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

#[test]
fn fixed_sector_stall_records_and_channel_gaps() {
    let fixture = Fixture::new();
    let (result, stderr) = fixture.run(20);
    assert_eq!(result["scf"]["driver"], "reference");
    assert_eq!(result["convergence"]["converged"], false);
    let reference = &result["reference_scf"];
    let rows = reference["iterations"].as_array().unwrap();
    assert_eq!(rows.len(), 20);
    let stalls = reference["stalls"].as_array().unwrap();
    // Whether the 20-map fixture reaches the stall predicate depends on the
    // platform's rounding (it does on macOS/aarch64, not on dev08 x86_64), so
    // require it only where it is known to occur; everywhere the records must
    // match the predicate recomputed below.
    if cfg!(all(target_os = "macos", target_arch = "aarch64")) {
        assert!(
            !stalls.is_empty(),
            "fixture must exercise a real stall: {stderr}"
        );
    }
    let mut expected = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        for field in ["channel_gap_ha", "channel_edge_residual_ha"] {
            let channels = row[field].as_array().unwrap();
            assert_eq!(channels.len(), 2);
            assert!(
                channels
                    .iter()
                    .all(|v| v.as_f64().is_some_and(f64::is_finite))
            );
        }
        assert_eq!(
            row["occupied_subspace_loss"],
            serde_json::json!([null, null])
        );
        // Recompute the documented predicate independently from serialized history.
        let triggered = i >= 3
            && (i - 2..=i).all(|j| {
                rows[j]["mixing"]["condition"]
                    .as_f64()
                    .is_some_and(|c| c > 1e12)
                    && rows[j]["dr2_ry"].as_f64().unwrap()
                        >= rows[j - 1]["dr2_ry"].as_f64().unwrap()
            });
        assert_eq!(row.get("stall").is_some(), triggered);
        if triggered {
            let stall = &row["stall"];
            assert_eq!(stall["iteration"], row["iteration"]);
            assert_eq!(stall["channel_gap_ha"], row["channel_gap_ha"]);
            assert_eq!(stall["condition"], row["mixing"]["condition"]);
            assert_eq!(stall["estimate_ha"], row["dr2_ha"]);
            assert!(stall["orbital_rotation_fraction"].is_null());
            let recommendation = stall["recommendation"].as_str().unwrap();
            assert!(recommendation.contains("different fractional (spherical) state"));
            assert!(recommendation.contains("BAANDI_REFERENCE_STRICT=on"));
            let prefix = format!("SCF stall: map={} ", i + 1);
            let lines: Vec<_> = stderr.lines().filter(|l| l.starts_with(&prefix)).collect();
            assert_eq!(lines.len(), 1);
            assert!(lines[0].contains("condition="));
            assert!(lines[0].contains("estimate="));
            assert!(lines[0].contains(" Ha; "));
            assert!(lines[0].ends_with(recommendation));
            expected.push(stall.clone());
        }
    }
    assert_eq!(*stalls, expected);
    assert_eq!(
        stderr
            .lines()
            .filter(|l| l.starts_with("SCF stall:"))
            .count(),
        stalls.len()
    );
    let (early, stderr) = fixture.run(3);
    assert_eq!(early["reference_scf"]["stalls"], serde_json::json!([]));
    assert!(!stderr.contains("SCF stall:"));
    assert!(
        early["reference_scf"]["iterations"]
            .as_array()
            .unwrap()
            .iter()
            .all(|r| r.get("stall").is_none())
    );
    // Both scalar serializers expose the same fields and preserve the history alias.
    let scalar = fs::read_to_string(fixture.0.join("input.toml"))
        .unwrap()
        .replace("[spin]\nmode='collinear'\nn_up=4\nn_down=2\n", "")
        .replace("max_iterations=3", "max_iterations=2");
    for points in ["[[0.0,0.0,0.0]]", "[[0.125,0.0,0.0,1.0]]"] {
        fs::write(
            fixture.0.join("input.toml"),
            scalar.replace("[[0.0,0.0,0.0]]", points),
        )
        .unwrap();
        let (result, stderr) = fixture.execute();
        assert!(!stderr.contains("SCF stall:"));
        assert_eq!(result["reference_scf"]["stalls"], serde_json::json!([]));
        let rows = result["reference_scf"]["iterations"].as_array().unwrap();
        assert_eq!(rows.len(), 2);
        for (i, row) in rows.iter().enumerate() {
            for field in ["channel_gap_ha", "channel_edge_residual_ha"] {
                assert!(row[field][0].as_f64().is_some_and(f64::is_finite));
                assert!(row[field][1].is_null());
            }
            assert_eq!(*row, result["convergence"]["history"][i]["reference_scf"]);
        }
    }
}
