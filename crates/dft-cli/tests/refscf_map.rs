//! External QE maps: slow material tests run explicitly in release. They test
//! the diagnostic's disagreement reporting, NOT roundoff agreement with QE.
use serde_json::Value;
use std::{fs, path::PathBuf, process::Command};
fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}
fn fixture() -> PathBuf {
    root().join("tests/reference/qe/i25-si8/si8-np1")
}
fn input() -> PathBuf {
    let upf = std::env::var_os("DFT_UPF")
        .expect("set DFT_UPF to the external Si.pz-vbc.UPF (see map-dev08-request.md)");
    let upf = fs::canonicalize(upf).unwrap();
    assert_eq!(
        dft_cli::run::sha256(&fs::read(&upf).unwrap()),
        "e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217"
    );
    let dir = root().join(format!(
        "target/scratch/refscf-map-test-{}",
        std::process::id()
    ));
    fs::create_dir_all(&dir).unwrap();
    let mut text =
        fs::read_to_string(root().join("docs/performance/lanes/RSCF/si8-qe-replay.toml")).unwrap();
    text = text.replace("/path/to/Si.pz-vbc.UPF", upf.to_str().unwrap());
    text.push_str("\nnbands = 20\n");
    let path = dir.join("input.toml");
    fs::write(&path, text).unwrap();
    path
}
fn maps(
    input: &std::path::Path,
    solver: &str,
    ethr: Option<&str>,
) -> (std::process::Output, Vec<Value>) {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_dft-cli"));
    cmd.args(["refscf-map"])
        .arg(input)
        .arg("--qe")
        .arg(fixture())
        .args(["--steps", "3", "--solver", solver]);
    if let Some(t) = ethr {
        cmd.args(["--ethr", t]);
    } else {
        cmd.arg("--ethr-from-stdout");
    }
    let out = cmd.output().unwrap();
    let rows = String::from_utf8_lossy(&out.stdout)
        .lines()
        .map(|s| serde_json::from_str(s).unwrap())
        .collect();
    (out, rows)
}
#[test]
fn map_cli_exposes_units_and_conflicts() {
    let out = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
        .args(["refscf-map", "--help"])
        .output()
        .unwrap();
    assert!(out.status.success());
    let help = String::from_utf8(out.stdout).unwrap();
    assert!(help.contains("threshold in Ry"));
    assert!(help.contains("legacy, reference"));
    let out = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
        .args([
            "refscf-map",
            "unused.toml",
            "--qe",
            "unused",
            "--ethr",
            "1e-4",
            "--ethr-from-stdout",
        ])
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&out.stderr).contains("cannot be used with"));
}
#[test]
#[ignore = "external UPF and material solves; DFT_UPF=/path/Si.pz-vbc.UPF cargo test --release -p dft-cli --test refscf_map qe_si8_maps -- --ignored --nocapture"]
fn qe_si8_maps() {
    let input = input();
    let evidence: Value =
        serde_json::from_str(include_str!("fixtures/refscf-map-si8.json")).unwrap();
    for solver in ["legacy", "reference"] {
        let (out, rows) = maps(&input, solver, None);
        assert_eq!(
            out.status.code(),
            Some(4),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        assert_eq!(rows.len(), 5, "{}", String::from_utf8_lossy(&out.stderr));
        assert_eq!(rows[0]["g_vectors"], 11936);
        assert_eq!(rows[0]["grid"], serde_json::json!([36, 36, 36]));
        for (k, ethr) in [1e-2, 1.37e-3, 6.11e-5].into_iter().enumerate() {
            let r = &rows[k + 1];
            // Preserve the measured discrepancy, without turning it into an
            // allowed QE error. The CLI still applies the 2.72e-16 equality gate.
            let recorded = &evidence["runs"][solver]["records"][k + 1];
            for statistic in ["max", "rms"] {
                assert!(
                    (r["all_g"][statistic].as_f64().unwrap()
                        - recorded["all_g"][statistic].as_f64().unwrap())
                    .abs()
                        < 1e-12,
                    "{solver} k={} measured {statistic} changed",
                    k + 1
                );
            }
            assert_eq!(r["iteration"], k + 1);
            assert_eq!(r["ethr_ry"], ethr);
            assert_eq!(r["ethr_ha"], ethr / 2.);
            assert_eq!(r["first_failing_term"], "rho_out all G");
            assert_eq!(r["lowest_10_nonzero_shells"].as_array().unwrap().len(), 10);
            assert!(r["input_roundtrip"]["max"].as_f64().unwrap() < 2.72e-16);
            assert!(
                r["potential"]["input_hartree_vs_analytic"]["max"]
                    .as_f64()
                    .unwrap()
                    < 1e-12
            );
            for key in ["input_e", "actual_e", "expected_e"] {
                assert!((r["total_charge"][key].as_f64().unwrap() - 32.).abs() < 1e-10);
            }
            assert!(r["all_g"]["max"].as_f64().unwrap() > 2.72e-16);
            assert!(r["density_distance_ry"].as_f64().unwrap() > 0.);
            assert_eq!(r["bands"]["qe_ev"].as_array().unwrap().len(), 20);
            assert_eq!(r["occupations"]["baandi_e"], r["occupations"]["qe_e"]);
            assert_eq!(r["energy"]["terms"].as_array().unwrap().len(), 5);
            assert!(
                r["inner_solve"]["converged"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .all(|v| v == true)
            );
            println!(
                "{solver} k={} max={} RMS={} metric={} Ry",
                k + 1,
                r["all_g"]["max"],
                r["all_g"]["rms"],
                r["density_distance_ry"]
            );
        }
        assert_eq!(rows[4]["first_failing_term"], "k=1 rho_out all G");
        assert_eq!(rows[4]["agreement"], false);
    }
}
#[test]
#[ignore = "external UPF and tight material solves; DFT_UPF=/path/Si.pz-vbc.UPF cargo test --release -p dft-cli --test refscf_map tight_solvers -- --ignored --nocapture"]
fn tight_solvers() {
    let input = input();
    let (a, legacy) = maps(&input, "legacy", Some("1e-12"));
    let (b, reference) = maps(&input, "reference", Some("1e-12"));
    assert_eq!(
        a.status.code(),
        Some(4),
        "{}",
        String::from_utf8_lossy(&a.stderr)
    );
    assert_eq!(
        b.status.code(),
        Some(4),
        "{}",
        String::from_utf8_lossy(&b.stderr)
    );
    assert_eq!(legacy.len(), 5);
    assert_eq!(reference.len(), 5);
    for k in 1..=3 {
        // Separate solvers at a common tight threshold: occupied eigenvalues,
        // charge and explicit orbital energies. No comparison to rounded QE
        // eigenvalues is used as an oracle for this internal regression.
        let a = &legacy[k];
        let b = &reference[k];
        let ea = a["bands"]["baandi_ev"].as_array().unwrap();
        let eb = b["bands"]["baandi_ev"].as_array().unwrap();
        let max = ea
            .iter()
            .zip(eb)
            .take(16)
            .map(|(a, b)| (a.as_f64().unwrap() - b.as_f64().unwrap()).abs())
            .fold(0., f64::max);
        assert!(max < 1e-8, "k={k} occupied band difference {max} eV");
        let de = (a["energy"]["baandi_direct_orbital_total_ha"]
            .as_f64()
            .unwrap()
            - b["energy"]["baandi_direct_orbital_total_ha"]
                .as_f64()
                .unwrap())
        .abs();
        // Fixed non-self-consistent input: orbital functional is not stationary
        // in these orbitals, so it need not have eigenvalue-level agreement.
        assert!(de < 1e-6, "k={k} direct functional difference {de} Ha");
        println!("tight k={k} occupied bands {max:.9e} eV; direct functional {de:.9e} Ha");
    }
}

#[test]
#[ignore = "external UPF and simulated MPI; DFT_UPF=/path/Si.pz-vbc.UPF cargo test --release -p dft-cli --test refscf_map qe_si8_reference_ranks -- --ignored --nocapture"]
fn qe_si8_reference_ranks() {
    use dft_engine::refscf::{density::DensityPlan, metric::Serial, qe_reference::QeDensity};
    use dft_engine::{ScfProblem, SpeciesPotential, parallel};
    use dft_mpi::{SimRank, simulate};
    let path = input();
    let config = dft_cli::config::Config::parse(&fs::read_to_string(&path).unwrap()).unwrap();
    let prep = config.prepare().unwrap();
    let pp =
        dft_upf::parse_upf(&fs::read_to_string(&config.pseudopotentials["Si"]).unwrap()).unwrap();
    let opts = dft_cli::run::options(&config).unwrap();
    let atoms: Vec<_> = config
        .system
        .atoms
        .iter()
        .map(|a| dft_geometry::Atom::new(&a.species, a.fractional))
        .collect();
    let mut plan = DensityPlan::new(&Serial, &prep.cell, prep.grid.dims(), 60., 60.).unwrap();
    let read = |name: &str| {
        QeDensity::decode(&fs::read(fixture().join(name).join("charge-density.dat")).unwrap())
            .unwrap()
            .aligned(&plan.miller)
            .unwrap()
    };
    let out = read("k1-out");
    let mix = read("k1-mix");
    let input: Vec<_> = mix
        .iter()
        .zip(&out)
        .map(|(m, y)| (*m - 0.7 * y) / 0.3)
        .collect();
    let mut planes = vec![0.; plan.plane_len()];
    plan.inverse(&Serial, &input, &mut planes).unwrap();
    let mut density = dft_physics::data::Array3::zeros(prep.grid.dims());
    plan.write_serial_full(&planes, &mut density).unwrap();
    let run = || {
        let mut problem = ScfProblem::new(
            &prep.cell,
            &atoms,
            &[SpeciesPotential {
                label: "Si",
                potential: &pp,
            }],
            opts.clone(),
        )
        .unwrap();
        problem
            .reference_fixed_map(&density, None, 0.005, true)
            .unwrap()
    };
    let serial = run();
    for ranks in [2, 3] {
        let maps = simulate(ranks, |comm: &SimRank| {
            parallel::set_communicator(Some(std::sync::Arc::new(comm.clone())));
            parallel::set_m2_override(Some(true));
            let result = run();
            parallel::set_communicator(None);
            parallel::set_m2_override(None);
            result
        });
        for r in maps {
            assert!(r.converged.iter().all(|c| *c));
            let max = r
                .density
                .iter()
                .zip(&serial.density)
                .map(|(a, b)| (a - b).abs())
                .fold(0., f64::max);
            assert!(max < 1e-12, "ranks={ranks} density {max} e/bohr^3");
            assert!((r.energy.total_hartree() - serial.energy.total_hartree()).abs() < 1e-10);
            assert_eq!(r.record.iterations.len(), serial.record.iterations.len());
            assert_eq!(r.record.h_columns, serial.record.h_columns);
            println!("simulated ranks={ranks} max density {max:.9e} e/bohr^3");
        }
    }
}
