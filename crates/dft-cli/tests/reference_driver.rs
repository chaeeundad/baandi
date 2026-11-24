//! Reference-driver boundary: actual CLI runs, including validation and JSON.
#[allow(dead_code)]
mod common;
use serde_json::Value;
use std::{fs, process::Command};

fn scientific(value: &Value) -> Value {
    match value {
        Value::Array(a) => Value::Array(a.iter().map(scientific).collect()),
        Value::Object(o) => Value::Object(
            o.iter()
                .filter(|(k, _)| {
                    !matches!(
                        k.as_str(),
                        "timings" | "wall_s" | "ritz_changes_ha" | "inner_clocks"
                    ) && !k.ends_with("_s")
                })
                .map(|(k, v)| (k.clone(), scientific(v)))
                .collect(),
        ),
        _ => value.clone(),
    }
}
#[test]
fn broyden_reference_inputs_and_verified_records() {
    let dir = std::env::temp_dir().join(format!("rscf2-cli-{}", std::process::id()));
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
    let mut previous: Option<Value> = None;
    for (name, driver, env, input_profile, env_profile, full) in [
        ("input", "reference", "off", "", None, false),
        ("env", "legacy", "on", "", None, false),
        (
            "full",
            "reference",
            "off",
            "[profile]\ndavidson_stages=true\n",
            None,
            true,
        ),
        (
            "override",
            "reference",
            "off",
            "[profile]\ndavidson_stages=true\n",
            Some("off"),
            false,
        ),
        (
            "env-full",
            "reference",
            "off",
            "[profile]\ndavidson_stages=false\n",
            Some("on"),
            true,
        ),
    ] {
        fs::write(
            dir.join("in.toml"),
            format!("{text}driver='{driver}'\n{input_profile}"),
        )
        .unwrap();
        let mut command = Command::new(env!("CARGO_BIN_EXE_dft-cli"));
        command
            .current_dir(&dir)
            .env("BAANDI_REFERENCE_SCF", env)
            .env("BAANDI_REFERENCE_STRICT", "off")
            .env_remove("BAANDI_DAVIDSON_STAGES")
            .args(["run", "in.toml", "--output-dir", name]);
        if let Some(value) = env_profile {
            command.env("BAANDI_DAVIDSON_STAGES", value);
        }
        let out = command.output().unwrap();
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        let v: Value =
            serde_json::from_slice(&fs::read(dir.join(name).join("results.json")).unwrap())
                .unwrap();
        assert_eq!(v["scf"]["driver"], "reference");
        assert_eq!(v["scf"]["mixing"]["method"], "broyden");
        assert_eq!(v["scf"]["mixing_metric"], "hartree");
        assert_eq!(v["execution"]["davidson_mode"], "reference_generalized");
        assert_eq!(v["verified"], true);
        assert_eq!(v["run"]["completed"], true);
        assert_eq!(
            v["reference_scf"]["diagnostics"]["schema"],
            "reference-scf-diagnostics-1"
        );
        assert_eq!(
            v["reference_scf"]["diagnostics"]["level"],
            if full { "full" } else { "lightweight" }
        );
        assert_eq!(
            v["execution"]["effective_solver"],
            "reference_block_davidson"
        );
        let rows = v["reference_scf"]["iterations"].as_array().unwrap();
        assert!(!rows.is_empty());
        for row in rows {
            assert_eq!(
                row["dr2_ry"].as_f64().unwrap(),
                2. * row["dr2_ha"].as_f64().unwrap()
            );
            assert!(row["h_columns"].as_u64().unwrap() > 0);
            assert!(row["active_bands_sum"].as_u64().unwrap() > 0);
            assert!(!row["reduced_problem_orders"].as_array().unwrap().is_empty());
            // REFSCF_LEAN_ENERGY (default on since 2026-11-26): intermediate maps may
            // carry a null audit; the delivered (last) map always has the full audit.
            let total = &row["energies_ha"]["total"];
            assert!(total.is_null() || total.as_f64().unwrap().is_finite());
            // Preserve both historical trace-reader paths, including full Ritz
            // matrices. No reader must infer missing work from an empty trace.
            let i = row["iteration"].as_u64().unwrap() as usize - 1;
            assert_eq!(*row, v["convergence"]["history"][i]["reference_scf"]);
        }
        assert!(
            rows.last().unwrap()["energies_ha"]["total"]
                .as_f64()
                .unwrap()
                .is_finite()
        );
        for call in rows
            .iter()
            .flat_map(|r| r["calls"].as_array().unwrap())
            .chain(
                v["reference_scf"]["run"]["verification_calls"]
                    .as_array()
                    .unwrap(),
            )
        {
            let trace = call["ritz_changes_ha"].as_array().unwrap();
            if full {
                assert_eq!(
                    trace.len(),
                    call["davidson_iterations"].as_u64().unwrap() as usize
                );
                assert!(trace.iter().all(|row| row.as_array().unwrap().len() == 3));
                assert!(call["inner_clocks"]["total_s"].as_f64().unwrap() > 0.);
                assert_eq!(
                    call["inner_clocks"]["stage_s"].as_object().unwrap().len(),
                    8
                );
            } else {
                assert!(trace.is_empty());
                assert!(call["inner_clocks"].is_null());
            }
        }
        for stage in dft_engine::refscf::record::INNER_STAGES {
            let expected: f64 = rows
                .iter()
                .flat_map(|r| r["calls"].as_array().unwrap())
                .chain(
                    v["reference_scf"]["run"]["verification_calls"]
                        .as_array()
                        .unwrap(),
                )
                .map(|c| c["inner_clocks"]["stage_s"][stage].as_f64().unwrap_or(0.))
                .sum();
            let actual = v["convergence"]["work"]["davidson_stages"]["stage_s"][stage]
                .as_f64()
                .unwrap();
            assert!(
                (expected - actual).abs() < 1e-12,
                "stage {stage} must include verification"
            );
        }
        assert!(
            v["reference_scf"]["wall_s"]["verification"]
                .as_f64()
                .unwrap()
                > 0.
        );
        if let Some(p) = previous {
            assert_eq!(v["total_energy"], p["total_energy"]);
            assert_eq!(v["eigenvalues"], p["eigenvalues"]);
            assert_eq!(
                scientific(&v["verification"]),
                scientific(&p["verification"])
            );
            assert_eq!(
                scientific(&v["convergence"]["history"]),
                scientific(&p["convergence"]["history"])
            );
            assert_eq!(
                scientific(&v["reference_scf"]["run"]["verification_calls"]),
                scientific(&p["reference_scf"]["run"]["verification_calls"])
            );
            assert_eq!(
                fs::read(dir.join("input/density.bin")).unwrap(),
                fs::read(dir.join(name).join("density.bin")).unwrap()
            );
        }
        previous = Some(v);
    }
    fs::write(dir.join("in.toml"), format!("{text}driver='legacy'\n")).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
        .current_dir(&dir)
        .env("BAANDI_REFERENCE_SCF", "off")
        .env("BAANDI_REFERENCE_STRICT", "off")
        .args(["run", "in.toml", "--output-dir", "rejected"])
        .output()
        .unwrap();
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("mixing"));
    fs::remove_dir_all(dir).unwrap();
}

/// Diagnosis against genuine, committed QE output (never a QE source import).
#[test]
#[ignore = "external Si UPF required; DFT_SI_UPF=... cargo test --release -p dft-cli --test reference_driver qe_start_density -- --ignored --nocapture"]
fn qe_start_density() {
    use dft_engine::{
        SpeciesPotential,
        refscf::{density::DensityPlan, metric::Serial, qe_reference::QeDensity},
        refstart::{self, StartInput},
    };
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let text =
        fs::read_to_string(root.join("docs/performance/lanes/RSCF/si8-qe-replay.toml")).unwrap();
    let config = dft_cli::config::Config::parse(&text).unwrap();
    let cell = config.cell().unwrap();
    let atoms: Vec<_> = config
        .system
        .atoms
        .iter()
        .map(|a| dft_engine::input::Atom::new(&a.species, a.fractional))
        .collect();
    let pp = dft_upf::parse_upf(
        &fs::read_to_string(std::env::var("DFT_SI_UPF").expect("DFT_SI_UPF required")).unwrap(),
    )
    .unwrap();
    let pots = [SpeciesPotential {
        label: "Si",
        potential: &pp,
    }];
    let options = dft_cli::run::options(&config).unwrap();
    let state = refstart::start(
        StartInput {
            cell: &cell,
            atoms: &atoms,
            potentials: &pots,
            options: &options,
        },
        &Serial,
    )
    .unwrap();
    let mut plan = DensityPlan::new(
        &Serial,
        &cell,
        state.grid.dims(),
        4. * options.cutoff_hartree,
        4. * options.cutoff_hartree,
    )
    .unwrap();
    let mut ours = vec![dft_physics::data::Complex64::default(); plan.miller.len()];
    plan.forward(&Serial, &state.density, &mut ours).unwrap();
    let read = |mode: &str| {
        QeDensity::decode(
            &fs::read(root.join(format!(
                "tests/reference/qe/i25-si8/si8-np1/k1-{mode}/charge-density.dat"
            )))
            .unwrap(),
        )
        .unwrap()
        .aligned(&plan.miller)
        .unwrap()
    };
    let output = read("out");
    let mixed = read("mix");
    let initial: Vec<_> = mixed
        .iter()
        .zip(output)
        .map(|(m, y)| (*m - 0.7 * y) / 0.3)
        .collect();
    let max = ours
        .iter()
        .zip(&initial)
        .map(|(a, b)| (*a - *b).norm())
        .fold(0., f64::max);
    let rms = (ours
        .iter()
        .zip(&initial)
        .map(|(a, b)| (*a - *b).norm_sqr())
        .sum::<f64>()
        / ours.len() as f64)
        .sqrt();
    let zero = plan.miller.iter().position(|m| *m == [0; 3]).unwrap();
    println!(
        "QE_START max_e_bohr3={max:.17e} rms_e_bohr3={rms:.17e} zero_delta_e_bohr3={:.17e}",
        (ours[zero] - initial[zero]).norm()
    );
    for i in 0..ours.len() {
        if (ours[i] - initial[i]).norm() == max {
            println!(
                "QE_START_WORST miller={:?} ours={:?} qe={:?}",
                plan.miller[i], ours[i], initial[i]
            );
        }
    }
    assert!((ours[zero] - initial[zero]).norm() < 1e-15);
    assert!(max.is_finite());
}
