//! Small Mac-only protocol and physics regression tests; no ML dependencies.
#[allow(dead_code)]
mod common;
use serde_json::Value;
use std::{
    fs,
    path::PathBuf,
    process::Command,
    sync::atomic::{AtomicUsize, Ordering},
};
static SERIAL: AtomicUsize = AtomicUsize::new(0);
struct Case {
    dir: PathBuf,
}
impl Case {
    fn new() -> Self {
        let dir = std::env::temp_dir().join(format!(
            "exp-scf-{}-{}",
            std::process::id(),
            SERIAL.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&dir).unwrap();
        let mut pp =
            common::gaussian_with(2., 1.8, false, 161).replace("element=\"Xx\"", "element=\"He\"");
        let parsed = dft_upf::parse_upf(&pp).unwrap();
        let array = |f: &dyn Fn(f64) -> f64| {
            parsed
                .mesh()
                .r()
                .iter()
                .map(|&r| format!("{:.17e}", f(r)))
                .collect::<Vec<_>>()
                .join(" ")
        };
        let start = pp.find("<PP_RHOATOM>").unwrap() + "<PP_RHOATOM>".len();
        let end = pp.find("</PP_RHOATOM>").unwrap();
        pp.replace_range(
            start..end,
            &array(&|r| {
                8. / (std::f64::consts::PI.sqrt() * 1.8_f64.powi(3))
                    * r
                    * r
                    * (-(r / 1.8).powi(2)).exp()
            }),
        );
        pp = pp.replace("number_of_wfc=\"0\"","number_of_wfc=\"1\"").replace("</UPF>",&format!("<PP_PSWFC><PP_CHI.1 l=\"0\" label=\"s\" occupation=\"2\">{}</PP_CHI.1></PP_PSWFC></UPF>",array(&|r|r*(-0.6*r*r).exp())));
        fs::write(dir.join("syn.xml"), pp).unwrap();
        let text = "[system]\nname='reference test'\nlattice_unit='bohr'\nlattice_vectors=[[7.0,0.0,0.0],[0.2,7.0,0.0],[0.0,0.1,7.0]]\n[[system.atoms]]\nspecies='He'\nfractional=[0.13,0.21,0.32]\n[pseudopotentials]\nHe='syn.xml'\n[basis]\necutwfc_ry=4.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[scf]\nnbands=3\nmax_iterations=100\nmixing='broyden'\nmixing_beta=0.7\nmixing_history=8\n";

        fs::write(dir.join("in.toml"), format!("{text}driver='reference'\n")).unwrap();
        Self { dir }
    }
    fn controller(&self, mode: &str) -> PathBuf {
        use std::os::unix::fs::PermissionsExt;
        let script =
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/exp_echo_controller.py");
        let path = self.dir.join(format!("controller-{mode}"));
        // Shell quoting remains valid for paths containing apostrophes.
        let quoted = script.to_str().unwrap().replace('\'', "'\\''");
        fs::write(
            &path,
            format!("#!/bin/sh\nexec python3 -I '{quoted}' {mode}\n"),
        )
        .unwrap();
        fs::set_permissions(&path, fs::Permissions::from_mode(0o700)).unwrap();
        path
    }
    fn run(&self, name: &str, env: &[(&str, &str)]) -> Value {
        let mut command = Command::new(env!("CARGO_BIN_EXE_dft-cli"));
        // Isolate every registered switch from the invoking shell.
        for switch in dft_numerics::switches::TABLE {
            command.env_remove(switch.name);
        }
        command
            .current_dir(&self.dir)
            .args(["run", "in.toml", "--output-dir", name]);
        command.envs(env.iter().copied());
        let out = command.output().unwrap();
        assert!(
            out.status.success(),
            "{}: {}",
            name,
            String::from_utf8_lossy(&out.stderr)
        );
        let result = self.result(name);
        assert_eq!(result["verified"], true);
        result
    }
    fn result(&self, name: &str) -> Value {
        serde_json::from_slice(&fs::read(self.dir.join(name).join("results.json")).unwrap())
            .unwrap()
    }
    fn trace(&self, name: &str) -> Vec<Value> {
        fs::read_to_string(self.dir.join(name).join("scf_trace.jsonl"))
            .unwrap()
            .lines()
            .map(|s| serde_json::from_str(s).unwrap())
            .collect()
    }
    fn identical(&self, a: &str, b: &str) {
        let a_result = self.result(a);
        let b_result = self.result(b);
        for key in ["total_energy", "verification"] {
            assert_eq!(
                scientific(&a_result[key]),
                scientific(&b_result[key]),
                "{key}"
            );
        }
        assert_eq!(
            scientific(&a_result["convergence"]["history"]),
            scientific(&b_result["convergence"]["history"])
        );
        assert_eq!(
            scientific(&a_result["reference_scf"]["iterations"]),
            scientific(&b_result["reference_scf"]["iterations"])
        );
        assert_eq!(
            fs::read(self.dir.join(a).join("density.bin")).unwrap(),
            fs::read(self.dir.join(b).join("density.bin")).unwrap()
        );
    }
}
impl Drop for Case {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.dir);
    }
}
fn scientific(v: &Value) -> Value {
    match v {
        Value::Array(a) => Value::Array(a.iter().map(scientific).collect()),
        Value::Object(o) => Value::Object(
            o.iter()
                .filter(|(k, _)| {
                    !k.ends_with("_s")
                        && !matches!(
                            k.as_str(),
                            "timings"
                                | "wall_s"
                                | "ritz_changes_ha"
                                | "ritz_values_ha"
                                | "inner_clocks"
                        )
                })
                .map(|(k, v)| (k.clone(), scientific(v)))
                .collect(),
        ),
        _ => v.clone(),
    }
}
#[test]
fn echo_and_trace_preserve_default_bits_at_both_diagnostic_levels() {
    let c = Case::new();
    let controller = c.controller("echo");
    c.run("default", &[]);
    assert!(!c.dir.join("default/scf_trace.jsonl").exists());
    c.run(
        "explicit-off",
        &[("EXP_SCF_TRACE", "off"), ("EXP_SCF_CONTROLLER", "")],
    );
    c.identical("default", "explicit-off");
    for level in ["off", "on"] {
        c.run(
            level,
            &[
                ("EXP_SCF_TRACE", "on"),
                ("EXP_SCF_CONTROLLER", controller.to_str().unwrap()),
                ("BAANDI_DAVIDSON_STAGES", level),
            ],
        );
        c.identical("default", level);
        let trace = c.trace(level);
        assert_eq!(
            trace.len(),
            c.result(level)["reference_scf"]["iterations"]
                .as_array()
                .unwrap()
                .len()
        );
        for (i, row) in trace.iter().enumerate() {
            assert_eq!(row["iteration"], i + 1);
            assert_eq!(row["controller_status"], "accepted");
            assert_eq!(row["observation"]["dr2_ry"], row["dr2_ry"]);
            assert!(row["controller_wall_s"].as_f64().unwrap() > 0.);
            assert!(
                row["map_wall_s"].as_f64().unwrap() >= row["controller_wall_s"].as_f64().unwrap()
            );
            assert!(row["ritz_last_max_change_ha"].is_number());
            if i > 0 {
                let delta = row["total_energy_ha"].as_f64().unwrap()
                    - trace[i - 1]["total_energy_ha"].as_f64().unwrap();
                assert!((delta - row["delta_total_energy_ha"].as_f64().unwrap()).abs() < 1e-14);
            }
        }
        assert_eq!(trace.last().unwrap()["verification_status"], "passed");
        assert_eq!(trace.last().unwrap()["action_applied"], false);
        let columns: u64 = trace.iter().map(|r| r["h_columns"].as_u64().unwrap()).sum();
        let verification_columns: u64 = trace
            .iter()
            .map(|r| r["verification_h_columns"].as_u64().unwrap())
            .sum();
        let iterations: u64 = trace
            .iter()
            .map(|r| r["davidson_iterations"].as_u64().unwrap())
            .sum();
        eprintln!(
            "echo {level}: maps={}, E={} Ha, ordinary H={columns}, verification H={verification_columns}, Davidson={iterations}, verification error={} Ha, density bit-identical",
            trace.len(),
            trace.last().unwrap()["total_energy_ha"],
            trace.last().unwrap()["verification_error_ha"]
        );
    }
}
#[test]
fn trace_only_is_inert() {
    let c = Case::new();
    c.run("default", &[]);
    c.run("trace", &[("EXP_SCF_TRACE", "on")]);
    c.identical("default", "trace");
    assert!(
        c.trace("trace")
            .iter()
            .all(|r| r["controller_status"] == "off")
    );
}
#[test]
fn malformed_missing_and_out_of_range_replies_fall_back() {
    let c = Case::new();
    c.run("default", &[]);
    for mode in ["garbage", "eof", "range", "missing"] {
        let path = c.controller(mode);
        c.run(mode, &[("EXP_SCF_CONTROLLER", path.to_str().unwrap())]);
        c.identical("default", mode);
        assert!(c.trace(mode).iter().all(|r| {
            r["controller_status"]
                .as_str()
                .unwrap()
                .starts_with("fallback:")
        }));
    }
}
#[test]
fn spawn_failure_falls_back() {
    let c = Case::new();
    c.run("default", &[]);
    c.run(
        "missing",
        &[("EXP_SCF_CONTROLLER", c.dir.join("absent").to_str().unwrap())],
    );
    c.identical("default", "missing");
    assert!(
        c.trace("missing")[0]["controller_status"]
            .as_str()
            .unwrap()
            .contains("spawn:")
    );
}
#[test]
fn seeded_exploration_repeats_full_trajectory() {
    let c = Case::new();
    for name in ["one", "two"] {
        c.run(name, &[("EXP_SCF_BETA_EXPLORE", "42")]);
    }
    c.identical("one", "two");
    assert_eq!(
        scientific(&Value::from(c.trace("one"))),
        scientific(&Value::from(c.trace("two")))
    );
    let first = c.trace("one");
    assert_eq!(first[0]["exploration_seed"], "42");
    assert!(
        first
            .iter()
            .all(|r| (0.2..=1.).contains(&r["beta_selected"].as_f64().unwrap()))
    );
    c.run("other", &[("EXP_SCF_BETA_EXPLORE", "43")]);
    assert_ne!(
        first[0]["exploration_draw_u64"],
        c.trace("other")[0]["exploration_draw_u64"]
    );
}
#[test]
fn reset_history_and_beta_above_one_are_applied() {
    let c = Case::new();
    let controller = c.controller("reset");
    c.run(
        "reset",
        &[("EXP_SCF_CONTROLLER", controller.to_str().unwrap())],
    );
    let rows = c.trace("reset");
    assert_eq!(rows[1]["beta_used"], 1.2);
    assert_eq!(rows[1]["reset_history"], true);
    assert_eq!(rows[1]["history_occupancy"], 0);
    assert_eq!(rows[2]["history_occupancy"], 1);
}
#[test]
fn experiments_reject_disabled_verification_and_conflicting_actions() {
    let c = Case::new();
    let text = fs::read_to_string(c.dir.join("in.toml")).unwrap();
    fs::write(c.dir.join("in.toml"), format!("{text}verify='off'\n")).unwrap();
    let run = |env: &[(&str, &str)]| {
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_dft-cli"));
        for switch in dft_numerics::switches::TABLE {
            cmd.env_remove(switch.name);
        }
        cmd.current_dir(&c.dir)
            .args(["run", "in.toml", "--output-dir", "reject"])
            .envs(env.iter().copied())
            .output()
            .unwrap()
    };
    let out = run(&[("EXP_SCF_TRACE", "on")]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("fixed_point"));
    fs::write(c.dir.join("in.toml"), text).unwrap();
    let out = run(&[
        ("EXP_SCF_CONTROLLER", "absent"),
        ("EXP_SCF_BETA_EXPLORE", "42"),
    ]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("mutually exclusive"));
}

#[test]
fn exploration_known_draw_and_broyden_reset_oracle() {
    use dft_engine::refscf::{
        broyden::Broyden,
        experiment::Exploration,
        metric::{Metric, Serial},
        profile::Profile,
    };
    use dft_physics::data::Complex64 as C;
    let mut rng = Exploration::new(42);
    assert_eq!(rng.draw().exploration_draw, Some(13679457532755275413));
    let profile = Profile::default();
    let mut mixer = Broyden::new(profile, &[1., 1.], &[true; 2], 1., 2.).unwrap();
    let mut metric = Metric::new(1., &[1., 1.], &[true; 2], 1, 350).unwrap();
    let x = [C::new(1., 0.), C::new(2., 0.)];
    let f = [C::new(0.5, 0.), C::new(-0.5, 0.)];
    let mut next = x;
    mixer
        .update(&Serial, &mut metric, &x, &f, &mut next)
        .unwrap();
    for bad in [0., -0.1, 1.5001, f64::NAN, f64::INFINITY] {
        assert!(mixer.action(bad, true).is_err());
    }
    mixer.action(1.5, true).unwrap();
    let diagnostics = mixer
        .update(&Serial, &mut metric, &x, &f, &mut next)
        .unwrap();
    assert_eq!(diagnostics.occupancy, 0);
    assert_eq!(next, [C::new(1.75, 0.), C::new(1.25, 0.)]);
}

#[test]
fn exhausted_budget_retains_trace_without_verified_delivery() {
    let c = Case::new();
    let text = fs::read_to_string(c.dir.join("in.toml")).unwrap();
    fs::write(
        c.dir.join("in.toml"),
        text.replace("max_iterations=100", "max_iterations=1"),
    )
    .unwrap();
    let path = c.controller("echo");
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_dft-cli"));
    for switch in dft_numerics::switches::TABLE {
        cmd.env_remove(switch.name);
    }
    let out = cmd
        .current_dir(&c.dir)
        .env("EXP_SCF_CONTROLLER", &path)
        .args(["run", "in.toml", "--output-dir", "exhausted"])
        .output()
        .unwrap();
    assert!(!out.status.success());
    assert_eq!(c.result("exhausted")["verified"], false);
    assert_eq!(c.trace("exhausted").len(), 1);
    assert_ne!(c.trace("exhausted")[0]["verification_status"], "passed");
}
