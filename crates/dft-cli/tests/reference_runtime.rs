//! Reference runtime lifecycle oracles; synthetic fixtures, no QE evidence.
#[allow(dead_code)]
#[path = "../../dft-engine/tests/common/mod.rs"]
mod common;
use dft_engine::{
    GammaReal, ScfOptions, SpeciesPotential, Verification,
    refscf::{
        driver,
        profile::Profile,
        runtime::{Control, Runtime},
    },
    refstart::StartInput,
    runtime::RuntimeState,
};
use dft_geometry::{Atom, Cell};
use dft_numerics::switches::RunConfig;
use sha2::{Digest, Sha256};

fn fingerprint(run: &driver::Run) -> (u64, String) {
    let result = &run.outcome.result;
    let mut hash = Sha256::new();
    for x in &result.density {
        hash.update(x.to_bits().to_le_bytes());
    }
    (
        result.total_hartree().to_bits(),
        format!("{:x}", hash.finalize()),
    )
}

#[test]
fn observation_and_pause_resume_preserve_bits_on_one_and_three_ranks() {
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0.2, 5., 0.], [0., 0.1, 5.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    let pp = common::potential_with_atomic_start(2., 1., true, 161, "SLA PZ");
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: 1.,
        bands: Some(4),
        gamma_real: GammaReal::On,
        max_iterations: 100,
        verification: Verification::Off,
        ..Default::default()
    };
    let mut observations = Vec::new();
    for verification in [Verification::Off, Verification::FixedPoint] {
        let options = ScfOptions {
            verification,
            ..options.clone()
        };
        for ranks in [1, 3] {
            let rows = dft_mpi::simulate(ranks, |comm| {
                let mut rows = Vec::new();
                let _scope = RunConfig::from_pairs(&[]).unwrap().enter();
                let input = || StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &options,
                };
                let plain = driver::run_with(input(), Profile::default(), comm).unwrap();
                assert!(plain.outcome.converged);
                for pause in [false, true] {
                    let mut runtime =
                        Runtime::new(driver::stepper_with(input(), Profile::default(), comm));
                    assert!(runtime.step().is_err());
                    runtime.initialize().unwrap();
                    assert!(runtime.initialize().is_err());
                    let mut reports = Vec::new();
                    while !runtime.snapshot().state.terminal() {
                        if pause && runtime.snapshot().safe_point {
                            let before = runtime.snapshot();
                            runtime.control(Control::Pause).unwrap();
                            assert_eq!(runtime.snapshot().state, RuntimeState::Paused);
                            assert!(runtime.step().is_err());
                            runtime.control(Control::Resume).unwrap();
                            assert_eq!(runtime.snapshot().completed_maps, before.completed_maps);
                            assert_eq!(runtime.snapshot().last_step, before.last_step);
                        }
                        let report = runtime.step().unwrap();
                        let event = dft_cli::telemetry::reference_step_json(
                            "reference-runtime-test",
                            reports.len() as u64 + 1,
                            runtime.snapshot().state_version,
                            std::time::UNIX_EPOCH,
                            report,
                        );
                        assert_eq!(event["data"]["h_columns"], report.h_columns);
                        assert_eq!(
                            event["data"]["energy_hartree"],
                            serde_json::json!(report.energy_ha)
                        );
                        assert!(event["data"]["density_residual_bohr_minus_3_over_2"].is_null());
                        reports.push(event);
                    }
                    let snapshot = runtime.snapshot();
                    assert_eq!(snapshot.completed_maps, reports.len());
                    assert!(snapshot.safe_point);
                    assert!(runtime.control(Control::Cancel).is_err());
                    let run = runtime.finish().unwrap();
                    assert_eq!(fingerprint(&plain), fingerprint(&run));
                    assert_eq!(
                        plain
                            .outcome
                            .result
                            .density
                            .iter()
                            .map(|x| x.to_bits())
                            .collect::<Vec<_>>(),
                        run.outcome
                            .result
                            .density
                            .iter()
                            .map(|x| x.to_bits())
                            .collect::<Vec<_>>()
                    );
                    rows.push(serde_json::json!({
                        "ranks":ranks, "rank":dft_mpi::Collective::rank(comm),
                        "pause_resume":pause, "verification":format!("{verification:?}"),
                        "energy_hartree":run.outcome.result.total_hartree(),
                        "energy_bits":fingerprint(&run).0.to_string(),
                        "density_sha256":fingerprint(&run).1,
                        "converged":run.outcome.converged,
                        "identity":true, "events":reports
                    }));
                    eprintln!(
                        "ranks={ranks} pause={pause} energy_bits={} density_sha256={} maps={}",
                        fingerprint(&run).0,
                        fingerprint(&run).1,
                        reports.len()
                    );
                }
                rows
            });
            observations.extend(rows.into_iter().flatten());
        }
    }
    evidence("observation_runs", serde_json::json!(observations));
}

#[test]
fn cancellation_preserves_partial_record_and_releases_continuation() {
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0., 5., 0.], [0., 0., 5.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    let pp = common::potential_with_atomic_start(2., 1., true, 161, "SLA PZ");
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: 1.,
        bands: Some(4),
        gamma_real: GammaReal::On,
        max_iterations: 20,
        ..Default::default()
    };
    for ranks in [1, 3] {
        dft_mpi::simulate(ranks, |comm| {
            let _scope = RunConfig::from_pairs(&[]).unwrap().enter();
            let input = || StartInput {
                cell: &cell,
                atoms: &atoms,
                potentials: &potentials,
                options: &options,
            };
            let mut runtime = Runtime::new(driver::stepper_with(input(), Profile::default(), comm));
            runtime.initialize().unwrap();
            runtime.step().unwrap();
            let before = runtime.snapshot();
            runtime.control(Control::Pause).unwrap();
            runtime.control(Control::Cancel).unwrap();
            while !runtime.snapshot().state.terminal() {
                assert_ne!(runtime.snapshot().state, RuntimeState::Paused);
                runtime.step().unwrap();
            }
            let after = runtime.snapshot();
            assert_eq!(after.state, RuntimeState::Cancelled);
            assert!(after.safe_point);
            assert!(after.completed_maps >= before.completed_maps);
            assert!(after.last_step.is_some());
            evidence_row(
                "cancellation",
                serde_json::json!({"ranks":ranks,"rank":dft_mpi::Collective::rank(comm),"state":"Cancelled","completed_maps":after.completed_maps,"safe_point":after.safe_point,"delivered":false}),
            );
            assert!(runtime.result().is_none());
            assert!(runtime.step().is_err());
            assert!(runtime.finish().is_err());
            let mut next = Runtime::new(driver::stepper_with(input(), Profile::default(), comm));
            next.initialize().unwrap();
            next.control(Control::Cancel).unwrap();
            assert_eq!(next.snapshot().completed_maps, 0);
        });
    }
}

#[allow(dead_code)]
#[path = "../src/heal.rs"]
mod heal;

#[test]
fn policy_changes_are_real_atomic_and_replayable() {
    use dft_engine::refscf::runtime::PolicyPatch;
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0.2, 5., 0.], [0., 0.1, 5.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    let pp = common::potential_with_atomic_start(2., 1., true, 161, "SLA PZ");
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: 1.,
        bands: Some(4),
        gamma_real: GammaReal::On,
        max_iterations: 100,
        ..Default::default()
    };
    for ranks in [1, 3] {
        dft_mpi::simulate(ranks, |comm| {
            let _scope = RunConfig::from_pairs(&[]).unwrap().enter();
            let input = || StartInput {
                cell: &cell,
                atoms: &atoms,
                potentials: &potentials,
                options: &options,
            };
            let mut previous = None;
            for _ in 0..2 {
                let mut runtime =
                    Runtime::new(driver::stepper_with(input(), Profile::default(), comm));
                runtime.initialize().unwrap();
                let mut changed = false;
                while !runtime.snapshot().state.terminal() {
                    let report = runtime.step().unwrap();
                    if !changed
                        && runtime.snapshot().safe_point
                        && !runtime.snapshot().state.terminal()
                    {
                        runtime.control(Control::Pause).unwrap();
                        let before = runtime.snapshot();
                        assert!(
                            runtime
                                .apply_policy(PolicyPatch {
                                    beta: Some(0.),
                                    ..Default::default()
                                })
                                .is_err()
                        );
                        assert_eq!(runtime.snapshot(), before);
                        assert!(
                            runtime
                                .apply_policy(PolicyPatch {
                                    threshold_ry: Some(1.),
                                    ..Default::default()
                                })
                                .is_err()
                        );
                        assert_eq!(runtime.snapshot(), before);
                        let record = runtime
                            .apply_policy(PolicyPatch {
                                beta: Some(0.3),
                                history: Some(4),
                                threshold_ry: Some(5e-9),
                                verification: Some(true),
                            })
                            .unwrap();
                        assert!(record.history_reset);
                        assert_eq!(record.iteration, report.iteration);
                        assert!(
                            runtime
                                .apply_policy(PolicyPatch {
                                    beta: Some(0.4),
                                    ..Default::default()
                                })
                                .is_err()
                        );
                        runtime.control(Control::Resume).unwrap();
                        changed = true;
                    } else if !report.safe_point {
                        assert!(
                            runtime
                                .apply_policy(PolicyPatch {
                                    beta: Some(0.4),
                                    ..Default::default()
                                })
                                .is_err()
                        );
                    }
                }
                assert!(changed);
                assert_eq!(runtime.policy_records().len(), 1);
                let run = runtime.finish().unwrap();
                assert!(run.outcome.converged);
                assert!(
                    run.records
                        .iter()
                        .any(|r| r.beta.to_bits() == 0.3f64.to_bits())
                );
                assert!(run.outcome.result.verification.as_ref().unwrap().passed());
                let bits = fingerprint(&run);
                if let Some(previous) = previous.replace(bits.clone()) {
                    assert_eq!(bits, previous);
                }
                evidence_row(
                    "policy_replays",
                    serde_json::json!({"ranks":ranks,"rank":dft_mpi::Collective::rank(comm),"energy_bits":bits.0.to_string(),"energy_ha":run.outcome.result.total_hartree(),"density_sha256":bits.1,"beta":0.3,"history":4,"threshold_ry":5e-9,"verification":true,"converged":true}),
                );
            }
            // The historical healing operation rejects before any state change.
            let short_options = ScfOptions {
                max_iterations: 1,
                ..options.clone()
            };
            let (run, record) = heal::reference_budget_demo(driver::stepper_with(
                StartInput {
                    options: &short_options,
                    ..input()
                },
                Profile::default(),
                comm,
            ))
            .unwrap();
            assert!(!run.outcome.converged);
            assert_eq!(record["status"], "not_applicable");
            assert_eq!(record["policy_revision"], 0);
            evidence_row("self_healing", record);
        });
    }
}

#[allow(dead_code)]
#[path = "common/mod.rs"]
mod cli_fixture;

#[test]
fn cli_plain_observe_and_agent_have_identical_outputs() {
    use std::{fs, process::Command};
    let root = std::env::var_os("BAANDI_TEST_TMP")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(std::env::temp_dir)
        .join(format!("r8c-cli-{}", std::process::id()));
    fs::create_dir_all(&root).unwrap();
    let mut pp =
        cli_fixture::gaussian_with(2., 1., true, 161).replace("element=\"Xx\"", "element=\"He\"");
    let parsed = dft_upf::parse_upf(&pp).unwrap();
    let radial = parsed
        .mesh()
        .r()
        .iter()
        .map(|r| format!("{:.17e}", r * (-r * r).exp()))
        .collect::<Vec<_>>()
        .join(" ");
    let start = pp.find("<PP_RHOATOM>").unwrap() + "<PP_RHOATOM>".len();
    let end = pp.find("</PP_RHOATOM>").unwrap();
    pp.replace_range(start..end, &radial);
    pp=pp.replace("number_of_wfc=\"0\"","number_of_wfc=\"1\"").replace("</UPF>",&format!("<PP_PSWFC><PP_CHI.1 l=\"0\" label=\"s\" occupation=\"2\">{radial}</PP_CHI.1></PP_PSWFC></UPF>"));
    fs::write(root.join("synthetic.upf"), pp).unwrap();
    let input = "[system]\nname='runtime'\nlattice_unit='bohr'\nlattice_vectors=[[5.0,0.0,0.0],[0.2,5.0,0.0],[0.0,0.1,5.0]]\n[[system.atoms]]\nspecies='Xx'\nfractional=[0.13,0.21,0.32]\n[pseudopotentials]\nXx='synthetic.upf'\n[basis]\necutwfc_ry=2.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[scf]\ndriver='reference'\nnbands=4\nmax_iterations=100\n";
    fs::write(root.join("input.toml"), input.replace("Xx", "He")).unwrap();
    let mut plain = None;
    // Some restricted runners deny filesystem Unix sockets while allowing pairs.
    // Report that capability explicitly; the channel unit test uses a real pair.
    #[cfg(unix)]
    let agent_available = match std::os::unix::net::UnixDatagram::bind(root.join("probe.sock")) {
        Ok(socket) => {
            drop(socket);
            fs::remove_file(root.join("probe.sock")).unwrap();
            true
        }
        Err(e) if e.kind() == std::io::ErrorKind::PermissionDenied => {
            eprintln!("CLI agent endpoint not exercised: sandbox denies Unix socket bind: {e}");
            evidence(
                "agent_endpoint",
                serde_json::json!({"status":"unavailable","reason":e.to_string(),"channel_pair_unit_test":"separate cargo lib test"}),
            );
            false
        }
        Err(e) => panic!("socket capability probe: {e}"),
    };
    #[cfg(not(unix))]
    let agent_available = false;
    let mut cli_rows = Vec::new();
    for mode in ["plain", "observe", "agent"] {
        if mode == "agent" && !agent_available {
            continue;
        }
        let output = root.join(mode);
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_dft-cli"));
        cmd.arg("run")
            .arg(root.join("input.toml"))
            .arg("--output-dir")
            .arg(&output)
            .env("BAANDI_REFERENCE_STRICT", "off")
            .stdin(std::process::Stdio::null());
        if mode != "plain" {
            cmd.arg(format!("--{mode}"));
        }
        let result = cmd.output().unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let result: serde_json::Value =
            serde_json::from_slice(&fs::read(output.join("results.json")).unwrap()).unwrap();
        assert_eq!(result["scf"]["driver"], "reference");
        let density = fs::read(output.join("density.bin")).unwrap();
        let energy = result["energies"]["total"]["ha"]
            .as_f64()
            .unwrap()
            .to_bits();
        cli_rows.push(serde_json::json!({"mode":mode,"driver":result["scf"]["driver"],"energy_bits":energy.to_string(),"energy_ha":result["energies"]["total"]["ha"],"density_sha256":format!("{:x}",Sha256::digest(&density))}));
        if let Some((bits, ref bytes)) = plain {
            assert_eq!(energy, bits);
            assert_eq!(&density, bytes);
        } else {
            plain = Some((energy, density));
        }
        if mode != "plain" {
            let events: Vec<serde_json::Value> = fs::read_to_string(output.join("events.jsonl"))
                .unwrap()
                .lines()
                .map(|l| serde_json::from_str(l).unwrap())
                .collect();
            assert_eq!(events.first().unwrap()["type"], "initialized");
            assert_eq!(events.last().unwrap()["type"], "converged");
            for (i, e) in events.iter().enumerate() {
                assert_eq!(e["sequence"], i + 1);
            }
            assert!(events.iter().any(|e| e["data"]["kind"] == "verification"));
            let summary: serde_json::Value =
                serde_json::from_slice(&fs::read(output.join("runtime.json")).unwrap()).unwrap();
            assert_eq!(summary["delivered"], true);
            assert_eq!(summary["policy_revision"], 0);
            if mode == "agent" {
                assert!(output.join("control-audit.jsonl").exists());
                assert!(!output.join("control/control.sock").exists());
            }
        }
    }
    evidence("cli_runs", serde_json::json!(cli_rows));
    fs::remove_dir_all(root).unwrap();
}

fn evidence(key: &str, value: serde_json::Value) {
    save_evidence(key, value, false);
}
fn evidence_row(key: &str, value: serde_json::Value) {
    save_evidence(key, value, true);
}
fn save_evidence(key: &str, value: serde_json::Value, append: bool) {
    if std::env::var_os("R8C_WRITE_EVIDENCE").is_none() {
        return;
    }
    static EVIDENCE: std::sync::OnceLock<std::sync::Mutex<serde_json::Value>> =
        std::sync::OnceLock::new();
    let mut report=EVIDENCE.get_or_init(|| std::sync::Mutex::new(serde_json::json!({
        "schema":"reference-runtime-1",
        "command":"R8C_WRITE_EVIDENCE=1 cargo test -p dft-cli --test reference_runtime -- --nocapture",
        "reference":"same-build plain reference driver; synthetic fixtures, no QE reference",
        "rank_scope":"engine tests use simulated ranks 1/3; CLI subprocess tests use one rank",
        "limits":["No native MPI timing claim","No budget-extension or numerical-instability healing claim","Event timestamps in engine fixtures use UNIX_EPOCH"]
    }))).lock().unwrap();
    if append {
        report
            .as_object_mut()
            .unwrap()
            .entry(key)
            .or_insert_with(|| serde_json::json!([]))
            .as_array_mut()
            .unwrap()
            .push(value);
    } else {
        report[key] = value;
    }
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/reference/runtime/generated/reference-results.json");
    std::fs::write(
        path,
        format!("{}\n", serde_json::to_string_pretty(&*report).unwrap()),
    )
    .unwrap();
}

#[test]
fn complex_k_and_spin_policy_hooks_match_plain_trajectories() {
    use dft_engine::{
        Functional,
        collinear::SpinSector,
        refscf::{
            electronic::{self, Policy},
            k_driver,
            runtime::{Completion, PolicyPatch},
        },
    };
    use dft_geometry::kpoints::KPoint;
    fn complete<T: Completion>(mut runtime: Runtime<'_, T>) -> Result<T, dft_engine::ScfError> {
        while !runtime.snapshot().state.terminal() {
            if let Err(e) = runtime.step() {
                assert_eq!(runtime.snapshot().state, RuntimeState::Failed);
                return Err(e);
            }
        }
        runtime.finish()
    }
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0.2, 5., 0.], [0., 0.1, 5.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    let points = [KPoint {
        fractional: [0., 0., 0.],
        weight: 1.,
    }];
    for spin in [false, true] {
        let pp = common::potential_with_atomic_start(
            if spin { 1. } else { 2. },
            1.8,
            true,
            161,
            if spin { "SLA PW" } else { "SLA PZ" },
        );
        let potentials = [SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let options = ScfOptions {
            cutoff_hartree: 2.5,
            bands: Some(4),
            gamma_real: GammaReal::Off,
            max_iterations: 100,
            functional: if spin {
                Functional::LdaPw
            } else {
                Functional::LdaPz
            },
            ..Default::default()
        };
        for ranks in [1, 3] {
            dft_mpi::simulate(ranks, |comm| {
                let _scope = RunConfig::from_pairs(&[]).unwrap().enter();
                let input = || StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &options,
                };
                let patch = PolicyPatch {
                    beta: Some(0.4),
                    history: Some(3),
                    threshold_ry: Some(5e-9),
                    verification: Some(false),
                };
                let final_profile = Profile {
                    beta: 0.4,
                    history: 3,
                    threshold_ry: 5e-9,
                    ..Default::default()
                };
                let final_options = ScfOptions {
                    verification: Verification::Off,
                    ..options.clone()
                };
                let outcome = if spin {
                    let policy = || {
                        Policy::FixedSpin(SpinSector {
                            n_up: 1,
                            n_down: 0,
                            axis_cartesian: [0., 0., 1.],
                            initial_moment_seed: Some(0.3),
                        })
                    };
                    let plain = electronic::run_with(
                        StartInput {
                            options: &final_options,
                            ..input()
                        },
                        final_profile,
                        &points,
                        policy(),
                        comm,
                    );
                    let mut runtime = Runtime::new(electronic::stepper_with(
                        input(),
                        Profile::default(),
                        &points,
                        policy(),
                        comm,
                    ));
                    runtime.initialize().unwrap();
                    runtime.apply_policy(patch).unwrap();
                    match (plain, complete(runtime)) {
                        (Ok(plain), Ok(run)) => {
                            assert_eq!(run.converged, plain.converged);
                            assert!(run.verification.is_none());
                            assert_eq!(
                                plain.energy.total_hartree().to_bits(),
                                run.energy.total_hartree().to_bits()
                            );
                            assert_eq!(
                                plain
                                    .density
                                    .iter()
                                    .flatten()
                                    .map(|x| x.to_bits())
                                    .collect::<Vec<_>>(),
                                run.density
                                    .iter()
                                    .flatten()
                                    .map(|x| x.to_bits())
                                    .collect::<Vec<_>>()
                            );
                            serde_json::json!({"converged":run.converged,"energy_bits_identical":true,"density_bits_identical":true})
                        }
                        (Err(a), Err(b)) => {
                            assert_eq!(a.to_string(), b.to_string());
                            serde_json::json!({"matching_numerical_rejection":a.to_string()})
                        }
                        _ => panic!("reference spin runtime changed success/failure"),
                    }
                } else {
                    let plain = k_driver::run_with(
                        StartInput {
                            options: &final_options,
                            ..input()
                        },
                        final_profile,
                        &points,
                        comm,
                    );
                    let mut runtime = Runtime::new(k_driver::stepper_with(
                        input(),
                        Profile::default(),
                        &points,
                        comm,
                    ));
                    runtime.initialize().unwrap();
                    runtime.apply_policy(patch).unwrap();
                    match (plain, complete(runtime)) {
                        (Ok(plain), Ok(run)) => {
                            assert_eq!(run.outcome.converged, plain.outcome.converged);
                            assert!(run.outcome.result.verification.is_none());
                            assert_eq!(fingerprint(&plain), fingerprint(&run));
                            serde_json::json!({"converged":run.outcome.converged,"energy_bits_identical":true,"density_bits_identical":true})
                        }
                        (Err(a), Err(b)) => {
                            assert_eq!(a.to_string(), b.to_string());
                            serde_json::json!({"matching_numerical_rejection":a.to_string()})
                        }
                        _ => panic!("reference k runtime changed success/failure"),
                    }
                };
                evidence_row(
                    "complex_policy_outcomes",
                    serde_json::json!({"spin":spin,"ranks":ranks,"rank":dft_mpi::Collective::rank(comm),"outcome":outcome,"verification":"off","beta":0.4,"history":3,"threshold_ry":5e-9}),
                );
            });
        }
    }
}

#[test]
fn output_policy_metadata_matches_applied_controls() {
    let root = std::env::temp_dir().join(format!("r8c-metadata-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    std::fs::write(root.join("runtime.json"),serde_json::to_vec(&serde_json::json!({"mode":"agent","policies":[{"beta":0.3,"history":4,"threshold_ry":1e-9,"verification":false}]})).unwrap()).unwrap();
    let mut result = serde_json::json!({"scf":{"mixing":{"method":"broyden","beta":0.7,"history":8},"standard_threshold_ha":5e-9},"single_path":{"verification":"fixed_point","standard_threshold_ha":5e-9},"execution":{"verification":"fixed_point","standard_threshold_ha":5e-9},"reference_scf":{}});
    dft_cli::telemetry::annotate_reference_result(&mut result, &root).unwrap();
    assert_eq!(result["scf"]["mixing"]["beta"], 0.3);
    assert_eq!(result["runtime"]["launch_policy"]["mixing"]["beta"], 0.7);
    assert_eq!(result["execution"]["verification"], "off");
    assert_eq!(result["single_path"]["standard_threshold_ha"], 5e-10);
    assert_eq!(result["reference_scf"]["threshold_ry"], 1e-9);
    std::fs::remove_dir_all(root).unwrap();
}
