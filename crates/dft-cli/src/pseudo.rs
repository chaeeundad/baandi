//! Metadata-only registry and cache adapter. No network client is linked here.
use clap::Subcommand;
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{
    collections::BTreeSet,
    fs,
    io::{Read, Write},
    path::{Path, PathBuf},
    time::{SystemTime, UNIX_EPOCH},
};
type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const REGISTRY: &str = include_str!("../../../registry/pseudopotentials.json");

#[derive(Subcommand)]
pub enum Command {
    List {
        #[arg(long)]
        element: Option<String>,
        #[arg(long)]
        supported: bool,
    },
    Show {
        id: String,
    },
    Verify {
        file: PathBuf,
        #[arg(long)]
        id: Option<String>,
    },
    Path,
    Suggest {
        #[arg(required = true, num_args = 1..)]
        elements: Vec<String>,
    },
    Fetch {
        id: String,
        #[arg(long)]
        yes: bool,
        #[arg(long)]
        offline: bool,
    },
}
fn str_field<'a>(e: &'a Value, key: &str) -> &'a str {
    e[key].as_str().unwrap_or("")
}
fn digest(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}
fn valid_hash(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn xc_supported(s: &str) -> bool {
    use dft_engine::Functional;
    Functional::LdaPz.matches_upf(s)
        || Functional::LdaPw.matches_upf(s)
        || Functional::Pbe.matches_upf(s)
}
fn reasons(e: &Value) -> Vec<String> {
    let mut r = Vec::new();
    if e["status"] != "verified" || !valid_hash(str_field(e, "sha256")) {
        r.push("File identity unverified or deprecated".into());
    }
    if !e["formats"]
        .as_array()
        .is_some_and(|a| a.iter().any(|v| v == "UPF2"))
    {
        r.push("Reader requires UPF v2".into());
    }
    if e["pseudo_type"] != "NC" {
        r.push("Only NC supported; PAW/ultrasoft need augmentation".into());
    }
    for (key, label) in [
        ("nlcc", "nonlinear core correction"),
        ("spin_orbit", "spin-orbit"),
        ("paw", "PAW"),
        ("ultrasoft", "ultrasoft"),
    ] {
        if e[key] != false {
            r.push(format!("{label} unsupported or unverified"));
        }
    }
    if !matches!(
        str_field(e, "relativistic"),
        "scalar" | "nonrelativistic" | "no"
    ) {
        r.push("Relativistic treatment unsupported or unverified".into());
    }
    let counts = e["projector_counts"].as_object();
    if !counts.is_some_and(|c| {
        c.iter().all(|(l, n)| {
            l.parse::<usize>().is_ok_and(|l| l <= 3) && n.as_u64().is_some_and(|n| n <= 3)
        }) && c.values().filter_map(Value::as_u64).sum::<u64>() <= 3
    }) {
        r.push("Reader supports three projectors TOTAL and l <= 3".into());
    }
    if !xc_supported(str_field(e, "xc")) {
        r.push("XC requires recognized LDA-PZ or PBE label".into());
    }
    if !e["valence_electrons"].as_f64().is_some_and(|v| v > 0.0) {
        r.push("Valence unverified".into());
    }
    r
}
fn load(path: Option<&Path>) -> Result<Vec<Value>> {
    let text = path
        .map(fs::read_to_string)
        .transpose()?
        .unwrap_or_else(|| REGISTRY.into());
    let root: Value = serde_json::from_str(&text)?;
    if root["schema_version"] != 1 {
        return Err("unknown registry schema version".into());
    }
    let entries = root["entries"].as_array().ok_or("missing entries")?.clone();
    let mut ids = BTreeSet::new();
    for e in &entries {
        let id = str_field(e, "id");
        if id.is_empty()
            || !id
                .bytes()
                .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b == b'-' || b == b'.')
            || id.starts_with('.')
            || !ids.insert(id)
        {
            return Err("invalid/duplicate registry id".into());
        }
        if !e["sha256"].is_null() && !valid_hash(str_field(e, "sha256")) {
            return Err("invalid SHA-256".into());
        }
        if e["sha256"].is_null() && e["status"] != "unverified" {
            return Err("null hash must be unverified".into());
        }
        if e["status"] == "verified"
            && (e["size_bytes"].as_u64().unwrap_or(0) == 0 || !e["source"]["file_url"].is_string())
        {
            return Err("verified entry missing size or URL".into());
        }
        if e["supported_by_engine"] == true
            && (!reasons(e).is_empty()
                || e["support_reasons"]
                    .as_array()
                    .is_none_or(|a| !a.is_empty()))
        {
            return Err(format!("{id}: contradictory engine support").into());
        }
    }
    Ok(entries)
}
fn entry<'a>(entries: &'a [Value], id: &str) -> Result<&'a Value> {
    entries
        .iter()
        .find(|e| e["id"] == id)
        .ok_or_else(|| format!("unknown registry id: {id}").into())
}
fn env_path(name: &str) -> Result<Option<PathBuf>> {
    match std::env::var_os(name) {
        Some(v) => {
            let p = PathBuf::from(v);
            if !p.is_absolute() {
                return Err(format!("{name} must be a nonempty absolute path").into());
            }
            Ok(Some(p))
        }
        None => Ok(None),
    }
}
fn cache_dir() -> Result<PathBuf> {
    if let Some(p) = env_path("BAANDI_PSEUDO_DIR")? {
        return Ok(p);
    }
    let base = if cfg!(target_os = "macos") {
        env_path("HOME")?
            .ok_or("HOME is unset")?
            .join("Library/Caches")
    } else if cfg!(windows) {
        env_path("LOCALAPPDATA")?.ok_or("LOCALAPPDATA is unset")?
    } else if let Some(p) = env_path("XDG_CACHE_HOME")? {
        p
    } else {
        env_path("HOME")?.ok_or("HOME is unset")?.join(".cache")
    };
    Ok(base.join("baandi/pseudopotentials"))
}
// Resolve existing ancestors before creating anything, including symlink parents.
fn outside_repository(path: &Path) -> Result<()> {
    if path
        .components()
        .any(|c| matches!(c, std::path::Component::ParentDir))
    {
        return Err("cache path must not contain ..".into());
    }
    let existing = path
        .ancestors()
        .find(|p| p.exists())
        .ok_or("cache has no existing ancestor")?
        .canonicalize()?;
    let repo = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()?;
    if existing.starts_with(repo) || existing.ancestors().any(|p| p.join(".git").exists()) {
        return Err("pseudopotential cache must be outside repositories".into());
    }
    Ok(())
}
fn elements(args: &[String]) -> Result<BTreeSet<String>> {
    let mut out = BTreeSet::new();
    for arg in args {
        let b = arg.as_bytes();
        let mut i = 0;
        while i < b.len() {
            if b[i].is_ascii_whitespace() || b[i] == b',' {
                i += 1;
                continue;
            }
            if !b[i].is_ascii_uppercase() {
                return Err(format!("invalid composition: {arg}; use SiO2 or Si O").into());
            }
            let start = i;
            i += 1;
            if i < b.len() && b[i].is_ascii_lowercase() {
                i += 1;
            }
            out.insert(arg[start..i].to_string());
            let n = i;
            while i < b.len() && b[i].is_ascii_digit() {
                i += 1;
            }
            if n < i && arg[n..i].parse::<u64>().ok().is_none_or(|n| n == 0) {
                return Err("stoichiometry must be a positive integer".into());
            }
        }
    }
    if out.is_empty() {
        return Err("no elements supplied".into());
    }
    Ok(out)
}
fn eligible(e: &Value) -> bool {
    e["supported_by_engine"] == true && reasons(e).is_empty()
}
fn check_bytes(e: &Value, bytes: &[u8]) -> Result<()> {
    if digest(bytes) != str_field(e, "sha256")
        || e["size_bytes"].as_u64() != Some(bytes.len() as u64)
    {
        return Err("SHA-256 or size mismatch; file not accepted".into());
    }
    Ok(())
}
/// P9A integration boundary: HTTPS only, no redirects, normal TLS verification,
/// bounded response size, timeout and errors on every non-200 response.
trait Fetcher {
    fn download(&self, url: &str, max_bytes: u64) -> Result<Vec<u8>>;
}
struct NetworkUnavailable;
impl Fetcher for NetworkUnavailable {
    fn download(&self, _url: &str, _max_bytes: u64) -> Result<Vec<u8>> {
        Err("network fetch is not implemented (P9A); use a verified cache with --offline".into())
    }
}
fn read_bounded(path: &Path, size: u64) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    fs::File::open(path)?
        .take(size.saturating_add(1))
        .read_to_end(&mut bytes)?;
    Ok(bytes)
}
fn fetch(e: &Value, yes: bool, offline: bool, transport: &dyn Fetcher) -> Result<PathBuf> {
    fetch_into(e, yes, offline, transport, &cache_dir()?)
}
fn fetch_into(
    e: &Value,
    yes: bool,
    offline: bool,
    transport: &dyn Fetcher,
    dir: &Path,
) -> Result<PathBuf> {
    println!(
        "License: {}\nAttribution: {}\nCitations: {}\nProject: {}",
        e["license"], e["attribution"], e["citations"], e["source"]["project_url"]
    );
    if e["license"]["verified"] != true || !e["license"]["identifier"].is_string() {
        return Err("license unverified: read the source page; automatic download refused".into());
    }
    if !eligible(e) {
        return Err(format!(
            "entry not usable: {} {:?}",
            e["support_reasons"],
            reasons(e)
        )
        .into());
    }
    outside_repository(dir)?;
    let name = format!("{}-{}.upf", str_field(e, "id"), str_field(e, "sha256"));
    let path = dir.join(name);
    let manifest = path.with_extension("json");
    let size = e["size_bytes"].as_u64().ok_or("missing size")?;
    if fs::symlink_metadata(&path).is_ok() {
        if !fs::symlink_metadata(&path)?.file_type().is_file() {
            return Err("cache file is not a regular file".into());
        }
        check_bytes(e, &read_bounded(&path, size)?)?;
        let m: Value = serde_json::from_str(&fs::read_to_string(&manifest)?)?;
        if m["id"] != e["id"]
            || m["url"] != e["source"]["file_url"]
            || m["sha256"] != e["sha256"]
            || m["downloaded_at_unix_seconds"].as_u64().is_none()
        {
            return Err("cache provenance mismatch".into());
        }
        return Ok(path);
    }
    if offline {
        return Err("offline cache miss".into());
    }
    if !yes {
        return Err(
            "explicit --yes required before network access; review license and citations above"
                .into(),
        );
    }
    let url = e["source"]["file_url"].as_str().ok_or("missing file URL")?;
    if !url.starts_with("https://") {
        return Err("download requires HTTPS; redirects are forbidden".into());
    }
    // Fetch to memory: mismatch drops all bytes and creates no download file.
    let bytes = transport.download(url, size)?;
    check_bytes(e, &bytes)?;
    fs::create_dir_all(dir)?;
    outside_repository(dir)?;
    let now = SystemTime::now().duration_since(UNIX_EPOCH)?;
    let staging = dir.join(format!(
        ".download-{}-{}",
        std::process::id(),
        now.as_nanos()
    ));
    let result = (|| -> Result<()> {
        let mut f = fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&staging)?;
        f.write_all(&bytes)?;
        f.sync_all()?;
        // Hard link is atomic and never replaces another process's cache entry.
        fs::hard_link(&staging, &path)?;
        let provenance = json!({"id":e["id"], "url":url, "sha256":e["sha256"], "downloaded_at_unix_seconds":now.as_secs()});
        let mut m = fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&manifest)?;
        m.write_all(serde_json::to_string_pretty(&provenance)?.as_bytes())?;
        m.sync_all()?;
        Ok(())
    })();
    let _ = fs::remove_file(staging);
    result?;
    Ok(path)
}
pub fn execute(command: &Command, registry: Option<&Path>) -> Result<()> {
    if matches!(command, Command::Path) {
        println!("{}", cache_dir()?.display());
        return Ok(());
    }
    let entries = load(registry)?;
    match command {
        Command::List { element, supported } => {
            for e in &entries {
                if element.as_ref().is_none_or(|x| e["element"] == *x)
                    && (!supported || eligible(e))
                {
                    println!(
                        "{} {} {} supported={} {}",
                        str_field(e, "id"),
                        str_field(e, "element"),
                        str_field(e, "status"),
                        eligible(e),
                        e["support_reasons"]
                    );
                }
            }
        }
        Command::Show { id } => println!("{}", serde_json::to_string_pretty(entry(&entries, id)?)?),
        Command::Verify { file, id } => {
            let bytes = fs::read(file)?;
            let hash = digest(&bytes);
            println!("SHA-256: {hash}");
            if let Some(id) = id {
                let e = entry(&entries, id)?;
                check_bytes(e, &bytes)?;
                println!(
                    "MATCH {id} status={} supported={}",
                    e["status"],
                    eligible(e)
                );
            } else {
                let matches: Vec<_> = entries.iter().filter(|e| e["sha256"] == hash).collect();
                if matches.is_empty() {
                    return Err(
                        "UNKNOWN: no registry SHA-256 matches; identity and support unverified"
                            .into(),
                    );
                }
                for e in matches {
                    check_bytes(e, &bytes)?;
                    println!(
                        "MATCH {} status={} supported={}",
                        str_field(e, "id"),
                        e["status"],
                        eligible(e)
                    );
                }
            }
        }
        Command::Suggest { elements: args } => {
            for element in elements(args)? {
                let candidates: Vec<_> =
                    entries.iter().filter(|e| e["element"] == element).collect();
                let mut usable: Vec<_> =
                    candidates.iter().copied().filter(|e| eligible(e)).collect();
                usable.sort_by(|a, b| {
                    let cutoff = |v: &Value| {
                        v["suggested_cutoffs"]["wavefunction_ry"]
                            .as_f64()
                            .unwrap_or(f64::INFINITY)
                    };
                    cutoff(a)
                        .total_cmp(&cutoff(b))
                        .then_with(|| str_field(a, "id").cmp(str_field(b, "id")))
                });
                if usable.is_empty() {
                    println!("{element}: no supported entry");
                    if candidates.is_empty() {
                        println!("  No inspected file registered; engine compatibility unverified");
                    }
                    for e in candidates {
                        println!(
                            "  {}: {} {:?}",
                            str_field(e, "id"),
                            e["support_reasons"],
                            reasons(e)
                        );
                    }
                }
                for e in usable {
                    println!(
                        "{element}: {} cutoff={} (Ry; ascending source cutoff, unknown last, then id; XC={})",
                        str_field(e, "id"),
                        e["suggested_cutoffs"],
                        e["xc"]
                    );
                }
            }
        }
        Command::Fetch { id, yes, offline } => println!(
            "{}",
            fetch(entry(&entries, id)?, *yes, *offline, &NetworkUnavailable)?.display()
        ),
        Command::Path => unreachable!(),
    }
    Ok(())
}

pub fn inspection_report(path: &Path, manifest: Option<&Path>) -> Result<Value> {
    let text = fs::read_to_string(path)?;
    let hash = digest(text.as_bytes());
    let mut report = match dft_upf::parse_upf(&text) {
        Ok(p) => {
            let mut counts = serde_json::Map::new();
            for b in p.projectors() {
                let n = counts.entry(b.l().to_string()).or_insert(json!(0));
                *n = json!(n.as_u64().unwrap() + 1);
            }
            let supported = xc_supported(p.functional());
            json!({"element":p.element(),"valence_electrons":p.z_valence(),"xc":p.functional(),"relativistic":p.raw().header.relativistic,"pseudo_type":"NC","nlcc":p.core_density().is_some(),"spin_orbit":false,"paw":false,"ultrasoft":false,"projector_counts":counts,"formats":["UPF2"],"supported_by_engine":supported,"support_reasons":if supported {vec![]} else {vec!["XC label unsupported"]}})
        }
        Err(e) => json!({"supported_by_engine":false,"support_reasons":[e.to_string()]}),
    };
    if let Ok(facts) = dft_upf::header_facts(&text) {
        let a = &facts.attributes;
        let flag = |key: &str| -> Option<bool> {
            match a.get(key)?.trim().to_ascii_lowercase().as_str() {
                "t" | ".t." | "true" | ".true." | "1" => Some(true),
                "f" | ".f." | "false" | ".false." | "0" => Some(false),
                _ => None,
            }
        };
        let fields = json!({"element":a.get("element").map(|s|s.trim()), "valence_electrons":a.get("z_valence").and_then(|s|s.trim().parse::<f64>().ok()), "xc":a.get("functional"), "relativistic":a.get("relativistic"), "pseudo_type":a.get("pseudo_type"), "nlcc":flag("core_correction"), "spin_orbit":flag("has_so"), "paw":flag("is_paw"), "ultrasoft":flag("is_ultrasoft"), "projector_counts":facts.projector_counts, "formats": if matches!(facts.version.as_str(),"2.0"|"2.0.0"|"2.0.1") {vec!["UPF2"]} else {vec!["UPF"]}});
        let integer = |key| a.get(key).and_then(|s| s.trim().parse::<usize>().ok());
        let real = |key| a.get(key).and_then(|s| s.trim().parse::<f64>().ok());
        report["header_facts"] = json!({
            "element": a.get("element").map(|s| s.trim()),
            "z_valence": real("z_valence"), "pseudo_type": a.get("pseudo_type"),
            "relativistic": a.get("relativistic"), "functional": a.get("functional"),
            "nlcc": flag("core_correction"), "number_of_wfc": integer("number_of_wfc"),
            "mesh_size": integer("mesh_size"), "number_of_proj": integer("number_of_proj"),
            "projector_channels": facts.projector_channels.iter().map(|&(index,l)| json!({"index":index,"l":l})).collect::<Vec<_>>(),
            "cutoff_hints_ry": {"wavefunction":real("wfc_cutoff"),"density":real("rho_cutoff")}
        });
        let mut blockers = facts.blockers;
        if let Some(reasons) = report["support_reasons"].as_array() {
            blockers.extend(
                reasons
                    .iter()
                    .filter_map(Value::as_str)
                    .map(|s| format!("upf.validation: {s}; supply a supported, valid dataset")),
            );
        }
        report["blockers"] = json!(blockers);
        if !blockers.is_empty() {
            report["supported_by_engine"] = json!(false);
        }
        for (key, value) in fields.as_object().unwrap() {
            report[key] = value.clone();
        }
    }
    report["sha256"] = json!(hash);
    report["schema_version"] = json!(1);
    report["license_reference"] = Value::Null;
    report["qualification_status"] = json!("unqualified");
    if report.get("blockers").is_none() {
        report["blockers"] = report["support_reasons"].clone();
    }
    if let Some(path) = manifest {
        let m: Value = serde_json::from_slice(&fs::read(path)?)?;
        let entries = m["datasets"]
            .as_array()
            .ok_or("manifest.datasets: expected array")?;
        let matches: Vec<_> = entries.iter().filter(|e| e["sha256"] == hash).collect();
        if matches.len() != 1 {
            return Err("manifest.sha256: expected exactly one matching file identity; verify manifest and file".into());
        }
        let entry = matches[0];
        if entry["element"] != report["header_facts"]["element"] {
            return Err("manifest.element: differs from file; repair manifest".into());
        }
        report["license_reference"] = entry["license"].clone();
        report["qualification_status"] = entry["qualification_status"].clone();
    }
    Ok(report)
}

pub fn inspect_json(path: &Path, manifest: Option<&Path>) -> Result<()> {
    println!(
        "{}",
        serde_json::to_string_pretty(&inspection_report(path, manifest)?)?
    );
    Ok(())
}
#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;
    struct Fake {
        calls: Cell<usize>,
        data: Vec<u8>,
    }
    impl Fetcher for Fake {
        fn download(&self, _url: &str, _max_bytes: u64) -> Result<Vec<u8>> {
            self.calls.set(self.calls.get() + 1);
            Ok(self.data.clone())
        }
    }
    #[test]
    fn fetch_boundary_verifies_before_publication_and_never_overwrites() {
        let dir = std::env::temp_dir().join(format!(
            "baandi-fetch-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let bytes = b"synthetic transport test";
        let mut e: Value = serde_json::from_str::<Value>(REGISTRY).unwrap()["entries"][0].clone();
        e["id"] = json!("synthetic");
        e["status"] = json!("verified");
        e["sha256"] = json!(digest(bytes));
        e["size_bytes"] = json!(bytes.len());
        e["source"]["file_url"] = json!("https://example.invalid/file");
        e["license"]["verified"] = json!(true);
        e["formats"] = json!(["UPF2"]);
        e["xc"] = json!("SLA PZ");
        e["valence_electrons"] = json!(4);
        e["relativistic"] = json!("scalar");
        e["projector_counts"] = json!({"0":1});
        e["supported_by_engine"] = json!(true);
        e["support_reasons"] = json!([]);
        for key in ["nlcc", "spin_orbit", "paw", "ultrasoft"] {
            e[key] = json!(false);
        }
        let fake = Fake {
            calls: Cell::new(0),
            data: bytes.to_vec(),
        };
        assert!(fetch_into(&e, false, false, &fake, &dir).is_err());
        assert!(fetch_into(&e, true, true, &fake, &dir).is_err());
        assert_eq!(fake.calls.get(), 0);
        let bad = Fake {
            calls: Cell::new(0),
            data: b"wrong".to_vec(),
        };
        assert!(fetch_into(&e, true, false, &bad, &dir).is_err());
        assert!(!dir.exists());
        let path = fetch_into(&e, true, false, &fake, &dir).unwrap();
        assert_eq!(fake.calls.get(), 1);
        assert_eq!(fs::read(&path).unwrap(), bytes);
        assert!(fetch_into(&e, false, true, &fake, &dir).is_ok());
        assert_eq!(fake.calls.get(), 1);
        assert_eq!(fs::read_dir(&dir).unwrap().count(), 2);
        fs::write(&path, b"different").unwrap();
        assert!(fetch_into(&e, true, false, &fake, &dir).is_err());
        assert_eq!(fake.calls.get(), 1);
        assert_eq!(fs::read(&path).unwrap(), b"different");
        fs::remove_dir_all(dir).unwrap();
    }
}
