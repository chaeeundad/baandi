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
        Self(fs::canonicalize(path).unwrap())
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
        let mut structure = report["structure"].clone();
        structure.as_object_mut().unwrap().remove("review");
        assert_eq!(
            structure,
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

fn codes(report: &Value) -> Vec<&str> {
    report["findings"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v["code"].as_str().unwrap())
        .collect()
}
impl Fixture {
    fn inspect(&self, text: &str) -> Value {
        fs::write(self.0.join("input.toml"), text).unwrap();
        preflight::inspect(&self.0.join("input.toml"))
    }
}

#[test]
fn geometry_policies_positive_negative_and_boundaries() {
    let f = Fixture::new(0.);
    let base = f.input("Si", false, "lda-pz", false);
    for (policy, code) in [
        ("minimum_distance_bohr=4.01", "PF_MIN_DISTANCE"),
        ("minimum_volume_bohr3=65.0", "PF_CELL_VOLUME"),
        ("slab_axis=2\nminimum_vacuum_bohr=4.01", "PF_SLAB_VACUUM"),
    ] {
        let report = f.inspect(&format!("{base}\n[preflight]\n{policy}\n"));
        assert_eq!(report["valid"], true);
        assert!(codes(&report).contains(&code), "{report}");
        let baseline = f.inspect(&base);
        assert!(!codes(&baseline).contains(&code));
    }
    let boundary = f.inspect(&format!("{base}\n[preflight]\nminimum_distance_bohr=4.0\nminimum_volume_bohr3=64.0\nslab_axis=2\nminimum_vacuum_bohr=4.0\nminimum_volume_ratio=1.0\n"));
    assert!(codes(&boundary).is_empty());
    let skew = base.replace("[0.0,4.0,0.0]", "[8.0,4.0,0.0]");
    let skew = f.inspect(&format!("{skew}\n[preflight]\nminimum_volume_ratio=0.5\n"));
    assert_eq!(skew["valid"], true);
    assert!(codes(&skew).contains(&"PF_CELL_SKEW"));
    assert_eq!(skew["structure"]["review"]["minimum_distance_bohr"], 4.0);
    // Skew slab uses perpendicular height, not the length of the slanted vector.
    let slab = base.replace("[0.0,0.0,4.0]", "[8.0,0.0,4.0]");
    let slab = f.inspect(&format!(
        "{slab}\n[preflight]\nslab_axis=2\nminimum_vacuum_bohr=4.1\n"
    ));
    assert_eq!(slab["structure"]["review"]["vacuum_gap_bohr"], 4.0);
    let planes = base.replace("[0.0,0.0,4.0]", "[0.0,0.0,20.0]").replace(
        "[pseudopotentials]",
        "[[system.atoms]]\nspecies='site.1'\nfractional=[0.0,0.0,0.2]\n[pseudopotentials]",
    );
    let planes = format!("{planes}\n[preflight]\nslab_axis=2\n");
    let (_, metrics) =
        preflight::structure::review(&Config::parse_unvalidated(&planes).unwrap()).unwrap();
    assert_eq!(metrics.vacuum_gap_bohr, Some(16.0));
    let mesh = base
        .replace("kpoints=[[0.0,0.0,0.0]]", "kgrid=[1,1,1]")
        .replace("ecutwfc_ry=3.0", "ecutwfc_ry=12.0");
    let policy = "\n[preflight]\nmaximum_k_spacing_inverse_bohr=1.0\n";
    assert!(codes(&f.inspect(&format!("{mesh}{policy}"))).contains(&"PF_KMESH_CELL"));
    let dense = f.inspect(&format!(
        "{}{policy}",
        mesh.replace("kgrid=[1,1,1]", "kgrid=[2,2,2]")
    ));
    assert_eq!(dense["valid"], true);
    assert!(!codes(&dense).contains(&"PF_KMESH_CELL"));
    assert_eq!(
        f.inspect(&format!("{base}{policy}"))["structure"]["review"]["k_spacing_inverse_bohr"],
        Value::Null
    );
}

#[test]
fn failures_are_typed_json_with_unchanged_exit_status() {
    let f = Fixture::new(0.);
    let base = f.input("Si", false, "lda-pz", false);
    for (text, code) in [
        ("invalid = [".into(), "PF_INPUT_INVALID"),
        (
            base.replace("[0.0,0.0,4.0]", "[0.0,0.0,0.0]"),
            "PF_CELL_VOLUME",
        ),
        (
            base.replace(
                "[pseudopotentials]",
                "[[system.atoms]]\nspecies='site.1'\nfractional=[1.0,-2.0,3.0]\n[pseudopotentials]",
            ),
            "PF_ATOM_OVERLAP",
        ),
        (f.input("Fe", false, "lda-pz", false), "PF_SPECIES_IDENTITY"),
    ] {
        let report = f.inspect(&text);
        assert_eq!(report["valid"], false, "{report}");
        assert_eq!(codes(&report), [code]);
        assert_eq!(report["coverage"], "stopped_at_first_error");
        let finding: preflight::Finding =
            serde_json::from_value(report["findings"][0].clone()).unwrap();
        assert_eq!(
            serde_json::to_value(finding).unwrap(),
            report["findings"][0]
        );
        let out = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .args([
                "check",
                f.0.join("input.toml").to_str().unwrap(),
                "--json",
                "--explain",
            ])
            .output()
            .unwrap();
        assert_eq!(out.status.code(), Some(1));
        assert_eq!(
            serde_json::from_slice::<Value>(&out.stdout).unwrap(),
            report
        );
    }
    assert!(codes(&f.inspect(&base)).is_empty());
    let path = f.0.join("synthetic.xml");
    let original = fs::read_to_string(&path).unwrap();
    fs::write(
        &path,
        original.replace("z_valence=\"4\"", "z_valence=\"3\""),
    )
    .unwrap();
    assert_eq!(codes(&f.inspect(&base)), ["PF_ELECTRON_PARITY"]);
    fs::write(&path, original).unwrap();
    assert!(codes(&f.inspect(&base)).is_empty());
    let missing = preflight::inspect(&f.0.join("missing.toml"));
    assert_eq!(codes(&missing), ["PF_INPUT_INVALID"]);
}

#[test]
fn periodic_minimum_matches_independent_brute_force_for_skew_and_slab() {
    let f = Fixture::new(0.);
    for rows in [
        [[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]],
        [[4., 0., 0.], [8., 4., 0.], [0., 0., 4.]],
        [[4., 0., 0.], [0., 4., 0.], [0., 0., 20.]],
        [[0., 4., 0.], [4., 0., 0.], [1., 2., 7.]],
    ] {
        let text = f.input("Si", false, "lda-pz", false)
            .replace("[[4.0,0.0,0.0],[0.0,4.0,0.0],[0.0,0.0,4.0]]", &format!("{rows:?}"))
            .replace("[pseudopotentials]", "[[system.atoms]]\nspecies='site.1'\nfractional=[0.49,0.49,0.1]\n[pseudopotentials]");
        // Review only; the oracle uses Cartesian row combinations over a wider cube.
        let config = Config::parse_unvalidated(&text).unwrap();
        let (_, metrics) = preflight::structure::review(&config).unwrap();
        let mut brute = f64::INFINITY;
        for delta in [[0., 0., 0.], [0.49, 0.49, 0.1]] {
            for i in -6..=6 {
                for j in -6..=6 {
                    for k in -6..=6 {
                        if delta == [0.; 3] && [i, j, k] == [0; 3] {
                            continue;
                        }
                        let frac = [
                            delta[0] + f64::from(i),
                            delta[1] + f64::from(j),
                            delta[2] + f64::from(k),
                        ];
                        let r: [f64; 3] =
                            std::array::from_fn(|c| (0..3).map(|v| rows[v][c] * frac[v]).sum());
                        brute = brute.min(r.iter().map(|x| x * x).sum::<f64>().sqrt());
                    }
                }
            }
        }
        assert!((metrics.minimum_distance_bohr - brute).abs() < 1e-12);
    }
}

#[test]
fn invalid_review_policy_and_search_limits_fail_explicitly() {
    let f = Fixture::new(0.);
    let base = f.input("Si", false, "lda-pz", false);
    for policy in [
        "minimum_distance_bohr=nan",
        "minimum_volume_ratio=1.1",
        "slab_axis=3",
        "minimum_vacuum_bohr=2.0",
    ] {
        assert_eq!(
            codes(&f.inspect(&format!("{base}\n[preflight]\n{policy}"))),
            ["PF_INPUT_INVALID"]
        );
    }
    let extreme = base.replace("[0.0,4.0,0.0]", "[4.0,0.000001,0.0]");
    let report = f.inspect(&extreme);
    assert_eq!(report["valid"], false);
    assert!(report.to_string().contains("budget exceeded"));
}

fn proposal(f: &Fixture, text: &str, candidate: &str) -> preflight::proposal::Proposal {
    use preflight::proposal::{digest, findings_bytes, potential_hashes};
    let input = f.0.join("input.toml");
    fs::write(&input, text).unwrap();
    let diff = format!(
        "--- a/input.toml\n+++ b/input.toml\n@@ -1,{} +1,{} @@\n{}{}",
        text.lines().count(),
        candidate.lines().count(),
        text.lines().map(|l| format!("-{l}\n")).collect::<String>(),
        candidate
            .lines()
            .map(|l| format!("+{l}\n"))
            .collect::<String>()
    );
    let before: toml::Value = toml::from_str(text).unwrap();
    let after: toml::Value = toml::from_str(candidate).unwrap();
    serde_json::from_value(json!({
        "schema_version":1,
        "input_sha256":digest(text.as_bytes()),
        "findings_sha256":digest(&findings_bytes(&input).unwrap()),
        "pseudopotentials":potential_hashes(&Config::parse_unvalidated(text).unwrap(), &input).unwrap(),
        "candidate_pseudopotentials":potential_hashes(&Config::parse_unvalidated(candidate).unwrap(), &input).unwrap(),
        "advisory":"Review cutoff and test convergence.",
        "diff":diff,
        "edits":[{"key_path":["basis","ecutwfc_ry"],"before":before["basis"]["ecutwfc_ry"],"after":after["basis"]["ecutwfc_ry"],"reason":"UPF suggests a higher cutoff.","finding_codes":["PF_CUTOFF_BELOW_SUGGESTED"]}],
        "candidate_input_sha256":digest(candidate.as_bytes()),
        "agent":"synthetic-test-client"
    })).unwrap()
}

#[test]
fn proposal_validation_rejects_stale_undeclared_and_malformed_edits() {
    use preflight::proposal::validate;
    let f = Fixture::new(40.);
    let text = f.input("Si", false, "lda-pz", false);
    let changed = text.replace("ecutwfc_ry=3.0", "ecutwfc_ry=4.0");
    let good = proposal(&f, &text, &changed);
    let input = f.0.join("input.toml");
    let (candidate, result) = validate(&input, &good).unwrap();
    assert_eq!(candidate, changed);
    assert_eq!(result["valid"], true);
    let value = serde_json::to_value(&good).unwrap();
    for field in [
        "schema_version",
        "input_sha256",
        "findings_sha256",
        "candidate_input_sha256",
        "agent",
        "advisory",
        "diff",
        "edits",
        "candidate_pseudopotentials",
    ] {
        let mut bad = value.clone();
        match field {
            "schema_version" => bad[field] = json!(2),
            "edits" => bad[field] = json!([]),
            "candidate_pseudopotentials" => bad[field] = json!({}),
            _ => bad[field] = json!(""),
        }
        let p = serde_json::from_value(bad).unwrap();
        assert!(validate(&input, &p).is_err(), "accepted invalid {field}");
    }
    for diff in [
        good.diff.replace("b/input.toml", "../../input.toml"),
        good.diff.replace("@@ -1,", "@@ -2,"),
        format!("{}--- a/second\n+++ b/second\n", good.diff),
        good.diff.replace("-name='synthetic'", "-name='stale'"),
    ] {
        let mut bad = value.clone();
        bad["diff"] = json!(diff);
        assert!(validate(&input, &serde_json::from_value(bad).unwrap()).is_err());
    }
    for field in ["reason", "finding_codes", "before", "key_path"] {
        let mut bad = value.clone();
        bad["edits"][0][field] = match field {
            "finding_codes" => json!(["PF_FAKE"]),
            "before" => json!(99.),
            "key_path" => json!(["xc", "functional"]),
            _ => json!(""),
        };
        assert!(validate(&input, &serde_json::from_value(bad).unwrap()).is_err());
    }
    let mut unknown = value.clone();
    unknown["auto_apply"] = json!(true);
    assert!(serde_json::from_value::<preflight::proposal::Proposal>(unknown).is_err());
    let extra = changed.replace("fractional=[0.0,0.0,0.0]", "fractional=[0.1,0.0,0.0]");
    let undeclared = proposal(&f, &text, &extra);
    assert!(validate(&input, &undeclared).is_err());
    fs::write(&input, format!("{text}\n")).unwrap();
    assert!(validate(&input, &good).is_err());
    fs::write(&input, &text).unwrap();
    let upf = f.0.join("synthetic.xml");
    fs::write(&upf, format!("{}\n", fs::read_to_string(&upf).unwrap())).unwrap();
    assert!(validate(&input, &good).is_err());
}

#[test]
fn explicit_review_is_additive_durable_and_never_applies_edits() {
    use preflight::proposal::{Decision, record};
    let f = Fixture::new(40.);
    let text = f.input("Si", false, "lda-pz", false);
    let changed = text.replace("ecutwfc_ry=3.0", "ecutwfc_ry=4.0");
    let p = proposal(&f, &text, &changed);
    let path = f.0.join("proposal.json");
    fs::write(&path, serde_json::to_vec_pretty(&p).unwrap()).unwrap();
    let output = f.0.join("results");
    fs::create_dir(&output).unwrap();
    fs::write(output.join("results.json"), "existing numerical results").unwrap();
    let input = f.0.join("input.toml");
    let sidecar = record(
        &input,
        &path,
        &output,
        Decision::Approved,
        "test operator",
        "reviewed exact diff",
    )
    .unwrap();
    let bytes = fs::read(&sidecar).unwrap();
    let audit: Value = serde_json::from_slice(&bytes).unwrap();
    assert_eq!(
        preflight::proposal::digest(audit["proposal_source"].as_str().unwrap().as_bytes()),
        audit["proposal_sha256"].as_str().unwrap()
    );
    assert_eq!(audit["applied"], false);
    assert_eq!(audit["decision"]["value"], "approved");
    assert_eq!(audit["candidate_input"], changed);
    assert_eq!(audit["candidate_preflight"]["valid"], true);
    assert_eq!(audit["proposal"], serde_json::to_value(&p).unwrap());
    assert_eq!(fs::read_to_string(&input).unwrap(), text);
    assert_eq!(
        fs::read_to_string(output.join("results.json")).unwrap(),
        "existing numerical results"
    );
    assert!(
        record(
            &input,
            &path,
            &output,
            Decision::Rejected,
            "test",
            "reconsidered"
        )
        .is_err()
    );
    assert_eq!(fs::read(&sidecar).unwrap(), bytes);
    assert!(record(&input, &path, &output, Decision::Approved, "", "reason").is_err());
    assert!(
        !fs::read_dir(&output).unwrap().any(|e| e
            .unwrap()
            .path()
            .extension()
            .is_some_and(|x| x == "tmp"))
    );
    let invalid_candidate = text.replace("ecutwfc_ry=3.0", "ecutwfc_ry=0.0");
    let invalid = proposal(&f, &text, &invalid_candidate);
    fs::write(&path, serde_json::to_vec(&invalid).unwrap()).unwrap();
    assert!(record(&input, &path, &output, Decision::Approved, "test", "reason").is_err());
    let rejected = record(
        &input,
        &path,
        &output,
        Decision::Rejected,
        "test",
        "invalid cutoff",
    )
    .unwrap();
    let rejected: Value = serde_json::from_slice(&fs::read(rejected).unwrap()).unwrap();
    assert_eq!(rejected["decision"]["value"], "rejected");
    assert_eq!(rejected["candidate_preflight"]["valid"], false);
    assert_eq!(
        rejected["candidate_preflight"]["findings_schema_version"],
        1
    );
}

#[test]
fn review_cli_requires_explicit_decision_and_rejects_symlink_destination() {
    let f = Fixture::new(40.);
    let text = f.input("Si", false, "lda-pz", false);
    let p = proposal(&f, &text, &text.replace("ecutwfc_ry=3.0", "ecutwfc_ry=4.0"));
    fs::write(f.0.join("proposal.json"), serde_json::to_vec(&p).unwrap()).unwrap();
    let preview = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
        .current_dir(&f.0)
        .args(["preflight-review", "input.toml", "proposal.json"])
        .output()
        .unwrap();
    assert_eq!(preview.status.code(), Some(0));
    assert!(String::from_utf8_lossy(&preview.stdout).contains("no decision recorded"));
    assert!(!f.0.join("results").exists());
    let mut command = Command::new(env!("CARGO_BIN_EXE_dft-cli"));
    command.current_dir(&f.0).args([
        "preflight-review",
        "input.toml",
        "proposal.json",
        "--output-dir",
        "results",
        "--actor",
        "operator",
        "--reason",
        "reviewed",
    ]);
    assert_eq!(command.output().unwrap().status.code(), Some(2));
    assert!(!f.0.join("results").exists());
    assert_eq!(
        command
            .args(["--decision", "approved"])
            .output()
            .unwrap()
            .status
            .code(),
        Some(0)
    );
    assert_eq!(fs::read_to_string(f.0.join("input.toml")).unwrap(), text);
    #[cfg(unix)]
    {
        std::os::unix::fs::symlink(f.0.join("results"), f.0.join("link")).unwrap();
        assert!(
            preflight::proposal::record(
                &f.0.join("input.toml"),
                &f.0.join("proposal.json"),
                &f.0.join("link"),
                preflight::proposal::Decision::Approved,
                "operator",
                "reviewed"
            )
            .is_err()
        );
    }
}
