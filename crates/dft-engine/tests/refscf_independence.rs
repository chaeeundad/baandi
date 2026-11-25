//! R9b boundary guard. No legacy-module exceptions, including unit-test fixtures.
//! `refstart`, shared force identities, and the synthetic `pool_fixture` are
//! intentional shared dependencies. Public compatibility re-exports live outside
//! the reference implementation and do not waive this check.
use std::path::{Path, PathBuf};

const LEGACY: &[&str] = &[
    "scf",
    "session",
    "convergence",
    "collinear",
    "finite_temperature",
    "cutoff_ladder",
    "gamma_real",
    "m2",
    "mix_dist",
    "band_density",
];

// A conservative grep-style scan: whitespace cannot hide a path, and grouped
// imports (including aliases and globs) are checked through their semicolon.
// Skip full-line comments so prose such as "use the collinear result" is not
// mistaken for an import. Strings and inline comments remain conservative.
fn legacy_paths(source: &str) -> Vec<&str> {
    let tokens: Vec<_> = source
        .lines()
        .filter(|line| !line.trim_start().starts_with("//"))
        .flat_map(|line| line.split_inclusive(|c: char| !c.is_ascii_alphanumeric() && c != '_'))
        .flat_map(|part| {
            let end = part.trim_end_matches(|c: char| !c.is_ascii_alphanumeric() && c != '_');
            [end, &part[end.len()..]]
        })
        .filter(|t| !t.is_empty() && !t.chars().all(char::is_whitespace))
        .collect();
    let mut in_use = false;
    let mut found = Vec::new();
    for (i, token) in tokens.iter().enumerate() {
        if *token == "use" {
            in_use = true;
        } else if *token == ";" {
            in_use = false;
        }
        if LEGACY.contains(token) && (in_use || tokens.get(i + 1..i + 3) == Some(&[":", ":"])) {
            found.push(*token);
        }
    }
    found
}

fn rust_sources(dir: &Path, files: &mut Vec<PathBuf>) {
    for entry in std::fs::read_dir(dir).unwrap() {
        let path = entry.unwrap().path();
        if path.is_dir() {
            rust_sources(&path, files);
        } else if path.extension().is_some_and(|ext| ext == "rs") {
            files.push(path);
        }
    }
}

#[test]
fn reference_modules_are_independent_of_legacy_drivers() {
    // Check the guard itself against common ways of spelling a dependency.
    for module in LEGACY {
        for source in [
            format!("crate::{module}::Thing"),
            format!("use crate::{{model::ScfError, {module} as old}};"),
            format!("use super :: {module} :: *;"),
            format!("use crate::{{{module}}};"),
        ] {
            assert!(!legacy_paths(&source).is_empty(), "missed {source}");
        }
    }
    assert!(legacy_paths("options.gamma_real == crate::model::GammaReal::On").is_empty());
    assert!(legacy_paths("use crate::{model::ScfError, refscf::driver};").is_empty());

    let src = Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
    let lib = std::fs::read_to_string(src.join("lib.rs")).unwrap();
    assert!(lib.contains("pub mod refscf;"));
    let scf = std::fs::read_to_string(src.join("scf.rs")).unwrap();
    assert!(!scf.contains("pub mod refscf;"));
    assert!(scf.contains("pub use crate::refscf;"));

    let mut files = vec![src.join("refstart.rs")];
    rust_sources(&src.join("refscf"), &mut files);
    files.sort();
    let mut violations = Vec::new();
    for path in files {
        let source = std::fs::read_to_string(&path).unwrap();
        for module in legacy_paths(&source) {
            violations.push(format!("{}: {module}", path.display()));
        }
    }
    assert!(
        violations.is_empty(),
        "legacy dependencies:\n{}",
        violations.join("\n")
    );
}
