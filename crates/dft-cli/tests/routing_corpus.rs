//! Corpus integrity and current parse expectations. Historical R0 routing
//! comments remain for provenance; driver_routing executes the R1 expectations.
use dft_cli::{config::Config, run_config::input_values};
use dft_numerics::switches::RunConfig;
use std::{collections::BTreeSet, fs, path::Path};

fn comment<'a>(text: &'a str, key: &str) -> &'a str {
    text.lines()
        .take_while(|line| line.starts_with('#'))
        .find_map(|line| line.strip_prefix(&format!("# {key} = ")))
        .filter(|value| !value.is_empty())
        .unwrap_or_else(|| panic!("missing leading {key} comment"))
}

#[test]
fn routing_fixtures_parse_with_current_reference_switch() {
    let _scope = RunConfig::from_pairs(&[("reference_scf", "on"), ("reference_strict", "off")])
        .unwrap()
        .enter();
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/routing");
    let mut files: Vec<_> = fs::read_dir(root)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.extension().is_some_and(|ext| ext == "toml"))
        .collect();
    files.sort();
    let mut names = BTreeSet::new();
    for path in &files {
        names.insert(path.file_stem().unwrap().to_str().unwrap().to_owned());
        let text = fs::read_to_string(path).unwrap();
        assert!(
            text.starts_with("# expected_driver = "),
            "{}",
            path.display()
        );
        assert!(matches!(
            comment(&text, "expected_driver"),
            "reference" | "legacy" | "rejected"
        ));
        for key in ["expected_reason", "command", "switches", "upf_requirements"] {
            comment(&text, key);
        }
        let expected = comment(&text, "expected_parse");
        match Config::parse(&text) {
            Ok(config) => {
                assert_eq!(expected, "ok", "{} unexpectedly accepted", path.display());
                // Type-check registry entries too; Config itself keeps these as TOML values.
                let input = input_values(&config.performance, &config.profile).unwrap();
                RunConfig::from_input(&input).unwrap();
            }
            Err(error) => {
                let fragment = expected
                    .strip_prefix("error:")
                    .unwrap_or_else(|| panic!("{} unexpectedly rejected: {error}", path.display()));
                assert!(
                    error.to_string().contains(fragment),
                    "{}: {error}",
                    path.display()
                );
                assert_eq!(comment(&text, "expected_driver"), "rejected");
            }
        }
    }
    for row in 1..=11 {
        assert!(
            names
                .iter()
                .any(|name| name.starts_with(&format!("l{row:02}-"))),
            "missing L{row}"
        );
    }
    for name in [
        "gamma-lda",
        "k-mesh",
        "fd",
        "fixed-spin",
        "symmetry",
        "pbe-nlcc",
        "relax-ions",
    ] {
        assert!(names.contains(name), "missing positive route {name}");
    }
    println!(
        "validated {} routing fixtures (configuration only)",
        files.len()
    );
}
