#[allow(dead_code)]
mod common;
use dft_engine::diagnostics::*;
use dft_engine::events::{EventKind, FailureCode, ScfEvent};
use dft_engine::*;
use dft_geometry::{Atom, Cell, kpoints::KPoint};

fn thresholds() -> Thresholds {
    Thresholds::from(&ScfOptions::default())
}
fn samples(residuals: &[f64]) -> Vec<Sample> {
    residuals
        .iter()
        .enumerate()
        .map(|(i, &r)| Sample {
            iteration: i + 1,
            policy_revision: 0,
            energy_hartree: -(i as f64),
            density_residual: r,
            charge_error_relative: 0.0,
            min_density: None,
            eigen_residual_hartree: None,
            input_charge_error_relative: None,
        })
        .collect()
}
fn evaluate(h: &[Sample]) -> Vec<Finding> {
    diagnose(h, &[], thresholds()).unwrap()
}
fn status(h: &[Sample], rule: Rule) -> Status {
    evaluate(h)
        .into_iter()
        .find(|f| f.rule == rule)
        .unwrap()
        .status
}
#[test]
fn patterns_and_evidence() {
    let plateau = samples(&[1.0; 7]);
    let f = evaluate(&plateau);
    assert_eq!(f, evaluate(&plateau));
    assert_eq!(f[2].status, Status::Triggered);
    assert_eq!(
        f[2].evidence,
        Evidence::Stagnation {
            residuals: [1.0; 6],
            range_ratio: 1.0,
            endpoint_ratio: 1.0
        }
    );
    assert_eq!(f[2].hypothesis, Some(Hypothesis::PoorMixing));
    assert_eq!(f[2].candidates, [Candidate::ReduceBetaAndResetHistory]);
    assert_eq!(f[4].status, Status::NotComputed);
    let mut oscillation = plateau.clone();
    for (i, s) in oscillation.iter_mut().enumerate() {
        s.energy_hartree = (i % 2) as f64;
    }
    let f = evaluate(&oscillation);
    assert_eq!(f[3].status, Status::Triggered);
    assert_eq!(
        f[3].evidence,
        Evidence::Oscillation {
            delta_energy_hartree: [1., -1., 1., -1., 1., -1.],
            sign_changes: 5,
            residual_ratio: 1.
        }
    );
    assert_eq!(
        status(&samples(&[1., 21., 22.]), Rule::Divergence),
        Status::Clear
    );
    let f = evaluate(&samples(&[1., 21., 22., 23.]));
    assert_eq!(f[1].status, Status::Triggered);
    assert_eq!(
        f[1].evidence,
        Evidence::Divergence {
            residual: 23.,
            prior_best: 1.,
            factor: 20.,
            count: 3
        }
    );
    assert_eq!(
        status(&samples(&[1., 21., 1., 22., 23.]), Rule::Divergence),
        Status::Clear
    );
    let healthy = samples(&[1., 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625]);
    assert!(
        evaluate(&healthy)
            .iter()
            .all(|f| f.status != Status::Triggered)
    );
}
#[test]
fn window_boundaries_zero_and_tolerance() {
    assert!(evaluate(&[]).iter().all(|f| f.status != Status::Triggered));
    assert_eq!(
        status(&samples(&[1.; 5]), Rule::Stagnation),
        Status::InsufficientEvidence
    );
    assert_eq!(
        status(&samples(&[1.; 6]), Rule::Oscillation),
        Status::InsufficientEvidence
    );
    for r in [0.0, thresholds().density_tolerance] {
        assert!(
            evaluate(&samples(&[r; 7]))
                .iter()
                .all(|f| f.status != Status::Triggered)
        );
    }
    assert_eq!(
        status(&samples(&[1., 1., 1., 1., 1., 1.21]), Rule::Stagnation),
        Status::Clear
    );
    assert_eq!(
        status(&samples(&[1., 1., 1., 1., 1., 0.94]), Rule::Stagnation),
        Status::Clear
    );
    let mut h = samples(&[1.; 7]);
    for (i, s) in h.iter_mut().enumerate() {
        s.energy_hartree = (i % 2) as f64 * thresholds().energy_tolerance_hartree;
    }
    assert_eq!(status(&h, Rule::Oscillation), Status::Clear);
    h[6].energy_hartree = f64::MAX;
    h[5].energy_hartree = -f64::MAX;
    assert_eq!(status(&h, Rule::Oscillation), Status::InsufficientEvidence);
}
#[test]
fn revisions_gaps_and_invalid_thresholds() {
    let mut h = samples(&[1.; 7]);
    h[6].policy_revision = 1;
    assert_eq!(status(&h, Rule::Stagnation), Status::InsufficientEvidence);
    h[6].policy_revision = 0;
    h[3].iteration += 1;
    assert_eq!(status(&h, Rule::Stagnation), Status::InsufficientEvidence);
    h[3].iteration -= 1;
    h[3].density_residual = -1.;
    assert_eq!(status(&h, Rule::Stagnation), Status::InsufficientEvidence);
    let mut t = thresholds();
    t.divergence_patience = 0;
    assert_eq!(diagnose(&h, &[], t), Err(InvalidThresholds));
    t = thresholds();
    t.density_tolerance = f64::NAN;
    assert_eq!(diagnose(&h, &[], t), Err(InvalidThresholds));
    let partial = samples(&[1., 21., 22., 23.]);
    assert_eq!(
        status(&partial[1..], Rule::Divergence),
        Status::InsufficientEvidence
    );
}
#[test]
fn nonfinite_is_classified_and_suppresses_advice() {
    for (v, class) in [
        (f64::NAN, NonFiniteClass::Nan),
        (f64::INFINITY, NonFiniteClass::PositiveInfinity),
        (f64::NEG_INFINITY, NonFiniteClass::NegativeInfinity),
    ] {
        let mut h = samples(&[1.; 7]);
        h[3].charge_error_relative = v;
        let f = evaluate(&h);
        assert_eq!(f, evaluate(&h));
        assert_eq!(f[0].status, Status::Triggered);
        assert_eq!(
            f[0].evidence,
            Evidence::NonFinite {
                iteration: 4,
                field: Field::ChargeErrorRelative,
                classification: class
            }
        );
        assert!(f[1..].iter().all(|f| f.status != Status::Triggered));
    }
    let mut h = samples(&[1.; 7]);
    h[6].min_density = Some(f64::NAN);
    assert_eq!(status(&h, Rule::NonFinite), Status::Triggered);
}
#[test]
fn sloshing_requires_explicit_consistent_probes_and_stagnation() {
    let h = samples(&[1.; 7]);
    let mut p: Vec<_> = (5..=7)
        .map(|iteration| FourierProbe {
            iteration,
            policy_revision: 0,
            low_q_power: 8.,
            total_power: 10.,
            q_min_inverse_bohr: 0.1,
            basis_id: "fixture-basis".into(),
        })
        .collect();
    let f = diagnose(&h, &p, thresholds()).unwrap();
    assert_eq!(f[4].status, Status::Triggered);
    assert_eq!(
        f[4].evidence,
        Evidence::Sloshing {
            fractions: [0.8; 3],
            q_min_inverse_bohr: 0.1,
            basis_id: "fixture-basis".into(),
            iterations: [5, 6, 7]
        }
    );
    p[1].low_q_power = 7.9;
    assert_eq!(
        diagnose(&h, &p, thresholds()).unwrap()[4].status,
        Status::Clear
    );
    p[1].total_power = 0.;
    assert_eq!(
        diagnose(&h, &p, thresholds()).unwrap()[4].status,
        Status::InsufficientEvidence
    );
    p[1].total_power = 10.;
    p[1].basis_id = "different".into();
    assert_eq!(
        diagnose(&h, &p, thresholds()).unwrap()[4].status,
        Status::InsufficientEvidence
    );
    p[1].basis_id = "fixture-basis".into();
    assert_eq!(
        diagnose(
            &samples(&[1., 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625]),
            &p,
            thresholds()
        )
        .unwrap()[4]
            .status,
        Status::Clear
    );
}
fn event(kind: EventKind) -> ScfEvent {
    ScfEvent {
        schema_version: "1",
        calc_id: "fixture".into(),
        sequence: 1,
        timestamp: std::time::SystemTime::UNIX_EPOCH,
        iteration: 1,
        state_version: 2,
        kind,
    }
}
#[test]
fn event_adapters_do_not_infer_from_untyped_errors_or_time() {
    let e = event(EventKind::IterationFinished {
        energy_hartree: -1.,
        density_residual: 1.,
        charge_error_relative: 0.,
        mixing_applied: true,
    });
    let s = Sample::from_event(&e, 3).unwrap();
    assert_eq!(s.policy_revision, 3);
    assert_eq!(s.energy_hartree, -1.);
    assert!(diagnose_failure(&e, 3).is_none());
    let mut e = event(EventKind::Failed {
        code: FailureCode::NonFinite,
    });
    let f = diagnose_failure(&e, 3).unwrap();
    e.timestamp += std::time::Duration::from_secs(10);
    assert_eq!(f, diagnose_failure(&e, 3).unwrap());
    assert_eq!(f.rule, Rule::NonFinite);
    assert!(Sample::from_event(&e, 3).is_none());
    assert!(
        diagnose_failure(
            &event(EventKind::Failed {
                code: FailureCode::Numerical
            }),
            0
        )
        .is_none()
    );
}

fn fixture(kind: usize, mixing: Mixing) -> (ScfProblem, ScfOptions) {
    let p = common::potential(2., 1., kind != 0, 401);
    let q = common::potential(2., 1.3, false, 401);
    let cell = Cell::from_user_rows_bohr(match kind {
        0 => [[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]],
        1 => [[5., 0., 0.], [0.5, 4., 0.], [0.2, 0.4, 6.]],
        2 => [[4., 0., 0.], [0., 4., 0.], [0., 0., 7.]],
        _ => [[4., 0., 0.], [0., 5., 0.], [0., 0., 8.]],
    })
    .unwrap();
    let mut atoms = vec![Atom::new("Xx", [0.13, 0.21, 0.32])];
    if kind == 1 {
        atoms.push(Atom::new("Yy", [0.6, 0.7, 0.8]));
    }
    let species = [
        SpeciesPotential {
            label: "Xx",
            potential: &p,
        },
        SpeciesPotential {
            label: "Yy",
            potential: &q,
        },
    ];
    let points: Vec<_> = (0..if kind == 2 { 4 } else { 1 })
        .map(|i| KPoint {
            fractional: [0.025 * i as f64, 0.015 * i as f64, 0.],
            weight: if kind == 2 { 0.25 } else { 1. },
        })
        .collect();
    let o = ScfOptions {
        cutoff_hartree: 2.5,
        bands: Some(if kind == 1 { 4 } else { 3 }),
        mixing,
        ..Default::default()
    };
    let problem = ScfProblem::with_kpoints(
        &cell,
        &atoms,
        &species[..if kind == 1 { 2 } else { 1 }],
        o.clone(),
        &points,
        OccupationPolicy::Insulator {
            minimum_gap_hartree: 1e-7,
        },
    )
    .unwrap();
    (problem, o)
}
/// Actual SCF runs on synthetic ions, not fabricated histories or QE references.
/// Cubic local, skew two-species, elongated four-k, and slab-like Gamma cells.
/// Small bands only; this does not calibrate a general default or large-band use.
#[test]
#[ignore = "SCF matrix; cargo test --release -p dft-engine --test diagnostics -- --include-ignored --nocapture"]
fn four_systems_have_no_false_positives_at_any_boundary() {
    for kind in 0..4 {
        let (problem, o) = fixture(kind, Mixing::default());
        let result = problem.run().unwrap();
        let h: Vec<_> = result
            .history
            .iter()
            .map(|r| Sample::from_record(r, 0))
            .collect();
        for n in 1..=h.len() {
            let f = diagnose(&h[..n], &[], Thresholds::from(&o)).unwrap();
            assert!(
                f.iter().all(|f| f.status != Status::Triggered),
                "fixture={kind} iteration={n}: {f:?}"
            );
            assert_eq!(f, diagnose(&h[..n], &[], Thresholds::from(&o)).unwrap());
        }
        println!(
            "fixture={kind} boundaries={} energy={:.17} Ha final_residual={:.9e} zero_false_positives",
            h.len(),
            result.total_hartree(),
            h.last().unwrap().density_residual
        );
    }
}
#[test]
#[ignore = "SCF failure; cargo test --release -p dft-engine --test diagnostics -- --include-ignored --nocapture"]
fn large_beta_failure_is_not_fabricated_as_residual_growth() {
    let (problem, o) = fixture(0, Mixing::Linear { beta: 100. });
    let mut s = ScfSession::new(
        problem,
        CompletionMode::Strict,
        ExecutionConfig {
            mode: ExecutionMode::Observed,
            ..Default::default()
        },
    )
    .unwrap();
    s.initialize().unwrap();
    let error = loop {
        if let Err(e) = s.step() {
            break e;
        }
    };
    assert!(matches!(error, ScfError::Diverged { .. }));
    let h: Vec<_> = s
        .history()
        .iter()
        .map(|r| Sample::from_record(r, 0))
        .collect();
    let findings = diagnose(&h, &[], Thresholds::from(&o)).unwrap();
    let mut failure = None;
    while let Some(e) = s.events().unwrap().pop() {
        if let Some(f) = diagnose_failure(&e, 0) {
            failure = Some(f);
        }
    }
    assert_eq!(failure.unwrap().rule, Rule::Divergence);
    println!(
        "beta=100 terminal divergence; completed_records={}; residual_growth_status={:?}; error={error}",
        h.len(),
        findings[1].status
    );
}

/// Re-run the P7PC eight-atom input with today's repaired buffer-band solver.
/// The historical 400-step/4.238824601895987e-6 Ha inner-solver failure is NOT a
/// completed density history. It cannot honestly test the SCF stagnation rule.
/// No real UPF is embedded; pass DFT_DIAGNOSTICS_SI_UPF to the read-only fixture.
#[test]
#[ignore = "real Si8; DFT_DIAGNOSTICS_SI_UPF=/path/to/Si.pz-vbc.UPF cargo test --release -p dft-engine --test diagnostics p7pc -- --ignored --nocapture"]
fn p7pc_current_si8_run() {
    let path =
        std::env::var("DFT_DIAGNOSTICS_SI_UPF").expect("read-only Si validation UPF required");
    let pp = dft_upf::parse_upf(&std::fs::read_to_string(path).unwrap()).unwrap();
    let a = 10.261_212_856_717_933;
    let cell = Cell::from_user_rows_bohr([[a, 0., 0.], [0., a, 0.], [0., 0., a]]).unwrap();
    let atoms: Vec<_> = [
        [0., 0., 0.],
        [0., 0.5, 0.5],
        [0.5, 0., 0.5],
        [0.5, 0.5, 0.],
        [0.25, 0.25, 0.25],
        [0.25, 0.75, 0.75],
        [0.75, 0.25, 0.75],
        [0.75, 0.75, 0.25],
    ]
    .into_iter()
    .map(|r| Atom::new("Si", r))
    .collect();
    let o = ScfOptions {
        cutoff_hartree: 15.,
        max_iterations: 100,
        mixing: Mixing::Pulay {
            beta: 0.5,
            history: 6,
        },
        ..Default::default()
    };
    let result = ScfProblem::new(
        &cell,
        &atoms,
        &[SpeciesPotential {
            label: "Si",
            potential: &pp,
        }],
        o.clone(),
    )
    .unwrap()
    .run()
    .unwrap();
    let h: Vec<_> = result
        .history
        .iter()
        .map(|r| Sample::from_record(r, 0))
        .collect();
    let mut triggers = 0;
    for n in 1..=h.len() {
        let f = diagnose(&h[..n], &[], Thresholds::from(&o)).unwrap();
        triggers += f.iter().filter(|f| f.status == Status::Triggered).count();
    }
    println!(
        "P7PC current solver: boundaries={} energy={:.17} Ha residual={:.9e} diagnostic_triggers={triggers}; historical inner-solver stagnation is not SCF stagnation",
        h.len(),
        result.total_hartree(),
        h.last().unwrap().density_residual
    );
}
