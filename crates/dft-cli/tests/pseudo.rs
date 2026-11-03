//! Registry process tests use only synthetic XML generated in memory.
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{
    fs,
    path::PathBuf,
    process::{Command, Output},
    sync::atomic::{AtomicUsize, Ordering},
};
fn fixture() -> String {
    let n = 401;
    let r = (0..n)
        .map(|i| format!("{}", i as f64 * 0.025))
        .collect::<Vec<_>>()
        .join(" ");
    let rab = vec!["0.025"; n].join(" ");
    // Smooth Coulomb-tailed synthetic potential, in Ry.
    let local = (0..n)
        .map(|i| format!("{}", -8.0 / ((i as f64 * 0.025).powi(2) + 1.0).sqrt()))
        .collect::<Vec<_>>()
        .join(" ");
    let beta = (0..n)
        .map(|i| {
            let r = i as f64 * 0.025;
            format!("{}", r * (-r * r).exp())
        })
        .collect::<Vec<_>>()
        .join(" ");
    let zero = vec!["0"; n].join(" ");
    format!(
        r#"<UPF version="2.0.1"><PP_INFO>SYNTHETIC CLI TEST ONLY</PP_INFO>
<PP_HEADER generated="synthetic" author="tests" date="2026" comment="synthetic" element="Xx"
pseudo_type="NC" relativistic="nonrelativistic" is_ultrasoft="F" is_paw="F" is_coulomb="F"
has_so="F" has_wfc="F" has_gipaw="F" core_correction="F" functional="SLA PZ"
z_valence="4" total_psenergy="0" wfc_cutoff="0" rho_cutoff="0" l_max="0" l_max_rho="0"
l_local="-1" mesh_size="{n}" number_of_wfc="0" number_of_proj="1"/>
<PP_MESH dx="0.025" mesh="{n}" xmin="0" rmax="10" zmesh="1"><PP_R>{r}</PP_R><PP_RAB>{rab}</PP_RAB></PP_MESH>
<PP_LOCAL>{local}</PP_LOCAL><PP_NONLOCAL><PP_BETA.1 index="1" angular_momentum="0" cutoff_radius="10" cutoff_radius_index="401">{beta}</PP_BETA.1><PP_DIJ>2</PP_DIJ></PP_NONLOCAL><PP_RHOATOM>{zero}</PP_RHOATOM></UPF>"#
    )
}
static NEXT: AtomicUsize = AtomicUsize::new(0);
struct Test {
    dir: PathBuf,
    registry: PathBuf,
    entry: Value,
}
impl Test {
    fn new() -> Self {
        let dir = std::env::temp_dir().join(format!(
            "baandi-registry-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&dir).unwrap();
        let registry = dir.join("registry.json");
        let data = fixture();
        let entry = json!({"id":"synthetic-si", "element":"Si", "source":{"family":"Synthetic", "version":"1", "file_url":"https://example.invalid/synthetic.upf", "project_url":"https://example.invalid/"}, "license":{"identifier":"CC0-1.0","url":"https://example.invalid/license", "verified":true},"attribution":"Synthetic test only", "citations":["Not a scientific reference"], "formats":["UPF2"], "xc":"SLA PZ", "valence_electrons":4, "relativistic":"nonrelativistic", "pseudo_type":"NC", "nlcc":false, "spin_orbit":false,"paw":false,"ultrasoft":false, "projector_counts":{"0":1}, "suggested_cutoffs":{"wavefunction_ry":30,"density_ry":120,"source":"https://example.invalid/cutoffs"}, "sha256":format!("{:x}",Sha256::digest(data.as_bytes())),"size_bytes":data.len(),"date_checked":"2026-11-03","supported_by_engine":true,"support_reasons":[],"status":"verified"});
        fs::write(dir.join("synthetic.xml"), data).unwrap();
        let t = Self {
            dir,
            registry,
            entry,
        };
        t.save(vec![t.entry.clone()]);
        t
    }
    fn save(&self, entries: Vec<Value>) {
        fs::write(
            &self.registry,
            serde_json::to_vec(&json!({"schema_version":1,"entries":entries})).unwrap(),
        )
        .unwrap();
    }
    fn command(&self) -> Command {
        let mut c = Command::new(env!("CARGO_BIN_EXE_dft-cli"));
        c.arg("pseudo")
            .arg("--registry")
            .arg(&self.registry)
            .env("BAANDI_PSEUDO_DIR", self.dir.join("cache"));
        c
    }
    fn run(&self, args: &[&str]) -> Output {
        self.command().args(args).output().unwrap()
    }
}
impl Drop for Test {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.dir).unwrap();
    }
}
fn out(o: &Output) -> String {
    format!(
        "{}{}",
        String::from_utf8_lossy(&o.stdout),
        String::from_utf8_lossy(&o.stderr)
    )
}
#[test]
fn offline_commands_and_unknown_identity() {
    let t = Test::new();
    for args in [
        vec!["list", "--element", "Si", "--supported"],
        vec!["show", "synthetic-si"],
        vec!["suggest", "Si2"],
    ] {
        let o = t.run(&args);
        assert!(o.status.success(), "{}", out(&o));
        assert!(out(&o).contains("synthetic-si"));
    }
    let o = t.run(&["show", "synthetic-si"]);
    for s in [
        "CC0-1.0",
        "Synthetic test only",
        "Not a scientific reference",
        "https://example.invalid/",
    ] {
        assert!(out(&o).contains(s));
    }
    assert!(!t.run(&["show", "missing"]).status.success());
    assert!(!t.run(&["suggest", "si2"]).status.success());
    assert!(!t.run(&["suggest", "Si0"]).status.success());
    assert!(t.run(&["list", "--element", "O"]).stdout.is_empty());
    let file = t.dir.join("synthetic.xml");
    let o = t.run(&["verify", file.to_str().unwrap()]);
    assert!(o.status.success(), "{}", out(&o));
    assert!(out(&o).contains("MATCH synthetic-si"));
    fs::write(&file, b"different synthetic bytes").unwrap();
    let o = t.run(&["verify", file.to_str().unwrap(), "--id", "synthetic-si"]);
    assert!(!o.status.success());
    assert!(out(&o).contains("mismatch"));
    let o = t.run(&["verify", file.to_str().unwrap()]);
    assert!(!o.status.success());
    assert!(out(&o).contains("UNKNOWN"));
}
#[test]
fn suggestions_rank_cutoffs_and_refuse_capability_gaps() {
    let t = Test::new();
    let mut nlcc = t.entry.clone();
    nlcc["id"] = json!("synthetic-fe");
    nlcc["element"] = json!("Fe");
    nlcc["nlcc"] = json!(true);
    nlcc["supported_by_engine"] = json!(false);
    nlcc["support_reasons"] = json!(["nonlinear core correction unsupported"]);
    let mut low = t.entry.clone();
    low["id"] = json!("low-si");
    low["suggested_cutoffs"]["wavefunction_ry"] = json!(20);
    let mut unknown = t.entry.clone();
    unknown["id"] = json!("unknown-cutoff");
    unknown["suggested_cutoffs"] = Value::Null;
    t.save(vec![unknown, t.entry.clone(), nlcc.clone(), low]);
    let o = t.run(&["suggest", "SiFe2O3"]);
    assert!(o.status.success());
    let text = out(&o);
    assert!(text.find("low-si").unwrap() < text.find("synthetic-si").unwrap());
    assert!(text.find("synthetic-si").unwrap() < text.find("unknown-cutoff").unwrap());
    assert!(text.contains("Fe: no supported entry"));
    assert!(text.contains("nonlinear core correction"));
    assert!(text.contains("O: no supported entry"));
    nlcc["supported_by_engine"] = json!(true);
    t.save(vec![nlcc]);
    assert!(!t.run(&["list"]).status.success());
    let mut four = t.entry.clone();
    four["projector_counts"] = json!({"0":2,"1":2});
    t.save(vec![four]);
    assert!(!t.run(&["suggest", "Si"]).status.success());
}
#[test]
fn cache_resolution_and_no_network_side_effects() {
    let t = Test::new();
    let o = t.run(&["path"]);
    assert_eq!(out(&o).trim(), t.dir.join("cache").to_str().unwrap());
    let o = t
        .command()
        .env("BAANDI_PSEUDO_DIR", "relative")
        .arg("path")
        .output()
        .unwrap();
    assert!(!o.status.success());
    let mut c = t.command();
    c.env_remove("BAANDI_PSEUDO_DIR")
        .env("HOME", &t.dir)
        .env_remove("XDG_CACHE_HOME")
        .env("LOCALAPPDATA", &t.dir);
    let o = c.arg("path").output().unwrap();
    assert!(o.status.success());
    let suffix = if cfg!(target_os = "macos") {
        "Library/Caches/baandi/pseudopotentials"
    } else if cfg!(windows) {
        "baandi/pseudopotentials"
    } else {
        ".cache/baandi/pseudopotentials"
    };
    assert_eq!(out(&o).trim(), t.dir.join(suffix).to_str().unwrap());
    if !cfg!(target_os = "macos") && !cfg!(windows) {
        let o = t
            .command()
            .env_remove("BAANDI_PSEUDO_DIR")
            .env("XDG_CACHE_HOME", &t.dir)
            .arg("path")
            .output()
            .unwrap();
        assert_eq!(
            out(&o).trim(),
            t.dir.join("baandi/pseudopotentials").to_str().unwrap()
        );
    }
    for (args, expected) in [
        (vec!["fetch", "synthetic-si"], "--yes required"),
        (
            vec!["fetch", "synthetic-si", "--offline"],
            "offline cache miss",
        ),
        (
            vec!["fetch", "synthetic-si", "--yes"],
            "not implemented (P9A)",
        ),
    ] {
        let o = t.run(&args);
        assert!(!o.status.success());
        assert!(out(&o).contains(expected), "{}", out(&o));
    }
    assert!(!t.dir.join("cache").exists());
    let o = t
        .command()
        .env(
            "BAANDI_PSEUDO_DIR",
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("bad-cache"),
        )
        .args(["fetch", "synthetic-si", "--yes"])
        .output()
        .unwrap();
    assert!(out(&o).contains("outside repositories"));
    let mut e = t.entry.clone();
    e["license"]["verified"] = json!(false);
    t.save(vec![e]);
    assert!(out(&t.run(&["fetch", "synthetic-si", "--yes"])).contains("license unverified"));
}
#[test]
fn cached_fetch_checks_hash_and_provenance() {
    let t = Test::new();
    let dir = t.dir.join("cache");
    fs::create_dir(&dir).unwrap();
    let file = dir.join(format!(
        "synthetic-si-{}.upf",
        t.entry["sha256"].as_str().unwrap()
    ));
    fs::write(&file, fixture()).unwrap();
    let manifest = file.with_extension("json");
    fs::write(&manifest,serde_json::to_vec(&json!({"id":"synthetic-si","url":t.entry["source"]["file_url"],"sha256":t.entry["sha256"],"downloaded_at_unix_seconds":1})).unwrap()).unwrap();
    let o = t.run(&["fetch", "synthetic-si", "--offline"]);
    assert!(o.status.success(), "{}", out(&o));
    assert!(out(&o).contains(file.to_str().unwrap()));
    fs::write(&manifest, "{}").unwrap();
    assert!(
        !t.run(&["fetch", "synthetic-si", "--offline"])
            .status
            .success()
    );
    fs::write(&file, "corrupted").unwrap();
    let o = t.run(&["fetch", "synthetic-si", "--yes"]);
    assert!(!o.status.success());
    assert!(out(&o).contains("mismatch"));
    assert_eq!(fs::read_to_string(&file).unwrap(), "corrupted");
}
#[test]
fn inspect_json_and_builder_use_strict_reader() {
    let t = Test::new();
    let file = t.dir.join("synthetic.xml");
    for (text, supported) in [
        (fixture(), true),
        (
            fixture().replace("core_correction=\"F\"", "core_correction=\"T\""),
            false,
        ),
        (fixture().replace("SLA PZ", "PBE"), true),
    ] {
        fs::write(&file, &text).unwrap();
        let o = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .args(["inspect-upf", file.to_str().unwrap(), "--json"])
            .output()
            .unwrap();
        assert!(o.status.success(), "{}", out(&o));
        let v: Value = serde_json::from_slice(&o.stdout).unwrap();
        assert_eq!(v["supported_by_engine"], supported);
        assert_eq!(v["element"], "Xx");
        assert_eq!(v["projector_counts"], json!({"0":1}));
        assert_eq!(v["nlcc"], !supported);
    }
    fs::write(&file, fixture()).unwrap();
    let metadata = t.dir.join("metadata.json");
    fs::write(
        &metadata,
        serde_json::to_vec(&json!({"synthetic.xml":t.entry})).unwrap(),
    )
    .unwrap();
    let output = t.dir.join("built.json");
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let o = Command::new("python3")
        .arg("-I")
        .arg(root.join("scripts/pseudo_registry_build.py"))
        .arg(&file)
        .arg("--metadata")
        .arg(metadata)
        .arg("--cli")
        .arg(env!("CARGO_BIN_EXE_dft-cli"))
        .arg("--output")
        .arg(&output)
        .output()
        .unwrap();
    assert!(o.status.success(), "{}", out(&o));
    let built: Value = serde_json::from_slice(&fs::read(output).unwrap()).unwrap();
    assert_eq!(built["entries"][0]["supported_by_engine"], true);
    assert_eq!(built["entries"][0]["sha256"], t.entry["sha256"]);
    assert_eq!(built["entries"][0]["element"], "Xx");
}
