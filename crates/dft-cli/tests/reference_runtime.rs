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
    let mut evidence = Vec::new();
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
            evidence.extend(rows.into_iter().flatten());
        }
    }
    if std::env::var_os("R8C_WRITE_EVIDENCE").is_some() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference/runtime/generated/reference-results.json");
        let value = serde_json::json!({
            "schema":"reference-runtime-foundation-1", "status":"partial",
            "command":"R8C_WRITE_EVIDENCE=1 cargo test -p dft-cli --test reference_runtime -- --nocapture",
            "reference":"same-build plain reference driver; synthetic skew-cell fixture, 1 Ha cutoff",
            "scope":"engine lifecycle adapter and CLI event serialization only; not CLI observe/agent dispatch",
            "ranks":"simulated collectives, not native MPI",
            "self_healing":{"status":"not_applicable","reason":"The stored legacy scenario extends max_iterations; budget extension is excluded by the R8c policy allow-list."},
            "limitations":["Policy mutation and transport integration pending owned-file scope approval", "L7 remains until full runtime integration", "Event timestamps fixed to UNIX_EPOCH for deterministic test capture"],
            "runs":evidence
        });
        std::fs::write(
            path,
            format!("{}\n", serde_json::to_string_pretty(&value).unwrap()),
        )
        .unwrap();
    }
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
    let _scope = RunConfig::from_pairs(&[]).unwrap().enter();
    let input = || StartInput {
        cell: &cell,
        atoms: &atoms,
        potentials: &potentials,
        options: &options,
    };
    let mut runtime = Runtime::new(driver::stepper(input(), Profile::default()));
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
    assert!(runtime.result().is_none());
    assert!(runtime.step().is_err());
    assert!(runtime.finish().is_err());
    let mut next = Runtime::new(driver::stepper(input(), Profile::default()));
    next.initialize().unwrap();
    next.control(Control::Cancel).unwrap();
    assert_eq!(next.snapshot().completed_maps, 0);
}
