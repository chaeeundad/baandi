//! Boundary-only measurement adapter; no solver or diagnostic reimplementation.
use dft_engine::{
    diagnostics::{Rule, Sample, Status, Thresholds, diagnose},
    policy::PolicyChange,
    runtime::*,
    *,
};
use dft_geometry::{Atom, Cell, kpoints::monkhorst_pack};
use serde::Deserialize;
use serde_json::{Value, json};
use std::{
    fs,
    time::{Duration, Instant, SystemTime},
};

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Input {
    id: String,
    lattice_bohr: [[f64; 3]; 3],
    atoms_fractional: Vec<[f64; 3]>,
    kgrid: [usize; 3],
    bands: usize,
    cutoff_hartree: f64,
    mixer: String,
    beta: f64,
    solver: String,
    max_iterations: usize,
    intervention: bool,
}
fn mixing(kind: &str, beta: f64) -> Mixing {
    match kind {
        "linear" => Mixing::Linear { beta },
        "pulay" => Mixing::Pulay { beta, history: 6 },
        "kerker" => Mixing::KerkerPulay {
            beta,
            history: 6,
            q0_inverse_bohr: 1.0,
        },
        _ => panic!("unknown mixer"),
    }
}
fn run(input: &Input, upf: &str) -> Result<Value, Box<dyn std::error::Error>> {
    let start = Instant::now();
    let cell = Cell::from_user_rows_bohr(input.lattice_bohr)?;
    let atoms: Vec<_> = input
        .atoms_fractional
        .iter()
        .map(|p| Atom::new("Xx", *p))
        .collect();
    let pp = dft_upf::parse_upf(upf)?;
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: input.cutoff_hartree,
        bands: Some(input.bands),
        mixing: mixing(&input.mixer, input.beta),
        max_iterations: input.max_iterations,
        eigensolver: match input.solver.as_str() {
            "davidson" => dft_solvers::SolverChoice::Davidson,
            "hybrid" => dft_solvers::SolverChoice::Hybrid(Default::default()),
            _ => panic!("unknown solver"),
        },
        ..Default::default()
    };
    let expanded = format!("{options:?}");
    let problem = if input.kgrid == [1, 1, 1] {
        ScfProblem::new(&cell, &atoms, &potentials, options)?
    } else {
        ScfProblem::with_kpoints(
            &cell,
            &atoms,
            &potentials,
            options,
            &monkhorst_pack(input.kgrid, [0.; 3])?,
            OccupationPolicy::Insulator {
                minimum_gap_hartree: 1e-6,
            },
        )?
    };
    let mut r = Runtime::new(
        ScfSession::new(
            problem,
            CompletionMode::RetainLast,
            ExecutionConfig {
                mode: ExecutionMode::AgentControlled,
                calc_id: input.id.clone(),
                ..Default::default()
            },
        )?,
        RuntimeConfig::default(),
    )
    .map_err(|e| format!("{e:?}"))?;
    r.initialize().map_err(|e| format!("{e:?}"))?;
    let mut samples = Vec::new();
    let mut findings = Vec::new();
    let mut history = Vec::new();
    let mut error = None;
    let mut terminal_evidence = Value::Null;
    let mut applied = false;
    while !r.state().terminal() {
        let step = r.step();
        // Drain every boundary; authoritative records below are not reconstructed from events.
        while let Some(e) = r.events().unwrap().pop() {
            assert!(!matches!(
                e.kind,
                dft_engine::events::EventKind::EventsDropped { .. }
            ));
        }
        if let Some(h) = r
            .session()
            .history()
            .last()
            .filter(|h| h.iteration > samples.len())
        {
            samples.push(Sample::from_record(h, r.policy_revision()));
            history.push(json!({"iteration":h.iteration,"revision":r.policy_revision(),"energy_hartree":h.energy.total_hartree(),"residual":h.density_residual,"certificate_debug":h.fixed_point_certificate.as_ref().map(|c|format!("{c:?}"))}));
            let f = diagnose(&samples, &[], Thresholds::from(&r.snapshot().options))
                .map_err(|e| format!("{e:?}"))?;
            let trigger = f.iter().any(|f| {
                f.status == Status::Triggered
                    && matches!(f.rule, Rule::Stagnation | Rule::Oscillation)
            });
            for finding in f {
                if finding.status == Status::Triggered {
                    findings.push(json!({"rule":format!("{:?}",finding.rule),"version":finding.rule_version,"first":finding.iteration_first,"last":finding.iteration_last,"revision":finding.policy_revision,"evidence_debug":format!("{:?}",finding.evidence)}));
                }
            }
            if input.intervention && trigger && !applied && !r.state().terminal() {
                let s = r.snapshot();
                let q=ControlRequest{schema_version:1,calc_id:input.id.clone(),id:"switch-pulay-01".into(),observed_state_version:s.state_version,observed_policy_revision:s.policy_revision,
                    action:ControlAction::Policy(vec![PolicyChange::Mixing(mixing("pulay",0.1)),PolicyChange::ResetHistory]),
                    preconditions:vec![Precondition::PolicyRevision(s.policy_revision),Precondition::Mode(ExecutionMode::AgentControlled)],
                    expiry:Expiry{max_iteration:Some(s.iteration),deadline:None},stale_policy:StalePolicy::Reject,reason:"P11D triggered stagnation or oscillation; switch to Pulay beta=0.1 and reset".into()};
                r.submit(
                    q,
                    "stagnation-harness".into(),
                    Authority::Approved(Approval {
                        id: "experiment-host".into(),
                        approver: "experiment-host".into(),
                        expires_at: SystemTime::now() + Duration::from_secs(120),
                        allow_revalidation: false,
                    }),
                )
                .map_err(|e| format!("{e:?}"))?;
                r.safe_point().map_err(|e| format!("{e:?}"))?;
                assert_eq!(r.audit().last().unwrap().outcome, AuditOutcome::Applied);
                applied = true;
            }
        }
        if let Err(e) = step {
            terminal_evidence = match &e {
                RuntimeError::Scf(ScfError::NonInsulating {
                    gap_hartree,
                    required_hartree,
                }) => {
                    json!({"code":"NonInsulating","gap_hartree":gap_hartree,"required_hartree":required_hartree})
                }
                RuntimeError::Scf(ScfError::Diverged { .. }) => json!({"code":"Diverged"}),
                RuntimeError::Scf(ScfError::NonFinite) => json!({"code":"NonFinite"}),
                _ => json!({"code":"Other","details_available":"debug error only"}),
            };
            error = Some(format!("{e:?}"));
            break;
        }
    }
    let audit:Vec<_>=r.audit().iter().map(|a|json!({"outcome":format!("{:?}",a.outcome),"request_debug":format!("{:?}",a.request),"iteration":a.before.iteration,"before_revision":a.before.policy_revision,"after_revision":a.after.policy_revision,"history_reset":a.history_reset})).collect();
    let result = r.into_result();
    let converged = result.outcome.as_ref().is_some_and(|o| o.converged);
    let energy = result
        .outcome
        .as_ref()
        .filter(|o| o.converged)
        .map(|o| o.result.total_hartree());
    // Sampled indirect HOMO/LUMO gap, also for the legacy Gamma path.
    // Gamma does not enforce the multi-k insulating gate; report that distinction.
    let gaps = result.outcome.as_ref().map(|o| {
        let occupied = o.result.occupations.iter().filter(|n| **n > 0.0).count();
        let homo = o
            .result
            .k_points
            .iter()
            .map(|k| k.solution.eigenvalues[occupied - 1])
            .fold(f64::NEG_INFINITY, f64::max);
        let lumo = o
            .result
            .k_points
            .iter()
            .map(|k| k.solution.eigenvalues[occupied])
            .fold(f64::INFINITY, f64::min);
        lumo - homo
    });
    Ok(
        json!({"state":format!("{:?}",result.state),"converged":converged,"certified_energy_hartree":energy,"sampled_gap_hartree":gaps,"gap_from_converged_state":converged,"insulating_gap_gate":input.kgrid!=[1,1,1],"iterations":history.len(),"final_residual":history.last().map(|v|&v["residual"]),"history":history,"findings":findings,"charge_sloshing_hint":"not_computed","error":error,"terminal_evidence":terminal_evidence,"audit":audit,"intervention_applied":applied,"expanded_options_debug":expanded,"wall_seconds":start.elapsed().as_secs_f64()}),
    )
}
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let input: Input = serde_json::from_str(&fs::read_to_string(&args[1]).unwrap()).unwrap();
    let value = match run(&input, &fs::read_to_string(&args[2]).unwrap()) {
        Ok(v) => v,
        Err(e) => json!({"state":"SetupError","converged":false,"error":format!("{e:?}")}),
    };
    fs::write(&args[3], serde_json::to_string_pretty(&value).unwrap()).unwrap();
}
