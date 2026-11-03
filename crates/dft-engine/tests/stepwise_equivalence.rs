#[allow(dead_code)]
mod common;
use dft_engine::events::{EventKind, FailureCode};
use dft_engine::*;
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use std::time::{Duration, Instant};
fn normalized(mut r: ScfResult) -> ScfResult {
    r.setup_time = Duration::ZERO;
    r.total_time = Duration::ZERO;
    for h in &mut r.history {
        h.timings = PhaseTimings::default();
        if let Some(c) = &mut h.fixed_point_certificate {
            c.timings = PhaseTimings::default();
        }
    }
    r
}
fn assert_identical(a: &ScfResult, b: &ScfResult) {
    assert_eq!(format!("{a:?}"), format!("{b:?}"));
    assert_eq!(a.total_hartree().to_bits(), b.total_hartree().to_bits());
    for (x, y) in a.history.iter().zip(&b.history) {
        assert_eq!(
            x.energy.total_hartree().to_bits(),
            y.energy.total_hartree().to_bits()
        );
        assert_eq!(x.density_residual.to_bits(), y.density_residual.to_bits());
    }
    assert_eq!(a.density.shape(), b.density.shape());
    for (x, y) in a.density.iter().zip(&b.density) {
        assert_eq!(x.to_bits(), y.to_bits());
    }
    for (x, y) in a.k_points.iter().zip(&b.k_points) {
        assert_eq!(x.solution, y.solution);
        for j in 0..x.solution.eigenvectors.cols() {
            for (v, w) in x
                .solution
                .eigenvectors
                .column(j)
                .iter()
                .zip(y.solution.eigenvectors.column(j))
            {
                assert_eq!(v.re.to_bits(), w.re.to_bits());
                assert_eq!(v.im.to_bits(), w.im.to_bits());
            }
        }
    }
}

fn fixture(kind: usize, mixing: Mixing, threads: usize, max_iterations: usize) -> ScfProblem {
    let p = common::potential(2.0, 1.0, true, 401);
    let q = common::potential(2.0, 1.3, false, 401);
    let cell = Cell::from_user_rows_bohr(if kind == 2 {
        [[5., 0., 0.], [0.5, 4., 0.], [0.2, 0.4, 6.]]
    } else {
        [[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]
    })
    .unwrap();
    let mut atoms = vec![Atom::new("Xx", [0.13, 0.21, 0.32])];
    if kind == 2 {
        atoms.push(Atom::new("Yy", [0.6, 0.7, 0.8]));
    }
    let points: Vec<_> = (0..if kind == 1 { 4 } else { 1 })
        .map(|i| KPoint {
            fractional: [0.025 * i as f64, 0.015 * i as f64, 0.],
            weight: if kind == 1 { 0.25 } else { 1. },
        })
        .collect();
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
    ScfProblem::with_kpoints(
        &cell,
        &atoms,
        &species[..if kind == 2 { 2 } else { 1 }],
        ScfOptions {
            threads,
            mixing,
            max_iterations,
            cutoff_hartree: 2.5,
            bands: Some(if kind == 2 { 4 } else { 3 }),
            ..Default::default()
        },
        &points,
        OccupationPolicy::Insulator {
            minimum_gap_hartree: 1e-7,
        },
    )
    .unwrap()
}
fn drive(problem: ScfProblem, mode: ExecutionMode) -> ScfResult {
    let mut session = ScfSession::new(
        problem,
        CompletionMode::Strict,
        ExecutionConfig {
            mode,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(session.snapshot().state_version, 0);
    assert!(session.result().is_none());
    session.initialize().unwrap();
    assert_eq!(session.snapshot().state_version, 1);
    let mut last_seq = 0;
    let mut steps = 0;
    loop {
        let report = session.step().unwrap();
        steps += 1;
        assert_eq!(report.iteration, steps);
        assert_eq!(report.state_version, steps as u64 + 1);
        let snap = session.snapshot();
        assert_eq!(snap.iteration, steps);
        assert_eq!(
            snap.energy.unwrap().total_hartree().to_bits(),
            session
                .history()
                .last()
                .unwrap()
                .energy
                .total_hartree()
                .to_bits()
        );
        if let Some(queue) = session.events() {
            while let Some(event) = queue.pop() {
                assert!(event.sequence > last_seq);
                last_seq = event.sequence;
                assert!(event.state_version <= report.state_version);
            }
        } else {
            assert_eq!(mode, ExecutionMode::Classic);
        }
        if report.status != StepStatus::Ready {
            assert_eq!(report.status, StepStatus::Converged);
            assert!(!report.mixing_applied);
            break;
        }
        assert!(report.mixing_applied);
    }
    assert!(session.step().is_err());
    session.into_result().unwrap().result
}
#[test]
#[ignore = "full matrix; cargo test --release -p dft-engine --test stepwise_equivalence -- --ignored --nocapture"]
fn modes_mixers_threads_and_species_are_bit_identical() {
    for kind in 0..3 {
        for mixing in [
            Mixing::Linear { beta: 0.3 },
            Mixing::Pulay {
                beta: 0.3,
                history: 6,
            },
            Mixing::KerkerPulay {
                beta: 0.3,
                history: 6,
                q0_inverse_bohr: 1.0,
            },
        ] {
            for threads in [1, 4] {
                let reference =
                    normalized(fixture(kind, mixing.clone(), threads, 200).run().unwrap());
                println!(
                    "CURRENT {kind} {mixing:?} {threads} {:016x}",
                    fingerprint(&reference)
                );
                for mode in [
                    ExecutionMode::Classic,
                    ExecutionMode::Observed,
                    ExecutionMode::AgentControlled,
                ] {
                    let result =
                        normalized(drive(fixture(kind, mixing.clone(), threads, 200), mode));
                    assert_identical(&reference, &result);
                }
                println!(
                    "fixture={kind} mixing={mixing:?} threads={threads}: all modes bit-identical; E={:.17} Ha; steps={}",
                    reference.total_hartree(),
                    reference.history.len()
                );
            }
        }
    }
}
#[test]
fn retained_and_strict_budget_semantics() {
    let mixing = Mixing::Linear { beta: 0.3 };
    let reference = normalized(
        fixture(0, mixing.clone(), 1, 1)
            .run_outcome()
            .unwrap()
            .result,
    );
    for mode in [
        ExecutionMode::Classic,
        ExecutionMode::Observed,
        ExecutionMode::AgentControlled,
    ] {
        let mut session = ScfSession::new(
            fixture(0, mixing.clone(), 1, 1),
            CompletionMode::RetainLast,
            ExecutionConfig {
                mode,
                ..Default::default()
            },
        )
        .unwrap();
        session.initialize().unwrap();
        let report = session.step().unwrap();
        assert_eq!(report.status, StepStatus::NotConverged);
        assert!(!report.mixing_applied);
        assert_identical(
            &reference,
            &normalized(session.into_result().unwrap().result),
        );
        let mut strict = ScfSession::new(
            fixture(0, mixing.clone(), 1, 1),
            CompletionMode::Strict,
            ExecutionConfig {
                mode,
                ..Default::default()
            },
        )
        .unwrap();
        strict.initialize().unwrap();
        assert_eq!(strict.step().unwrap().status, StepStatus::NotConverged);
        assert!(strict.result().is_none());
        assert!(strict.history()[0].mixing_history_used.is_some());
        assert_eq!(strict.snapshot().iteration, 1);
    }
}
#[test]
fn routine_loss_is_reported_and_critical_exhaustion_prevents_work() {
    let mut session = ScfSession::new(
        fixture(0, Mixing::Linear { beta: 0.3 }, 1, 200),
        CompletionMode::Strict,
        ExecutionConfig {
            mode: ExecutionMode::Observed,
            routine_capacity: 0,
            critical_capacity: 6,
            ..Default::default()
        },
    )
    .unwrap();
    session.initialize().unwrap();
    // Initialization consumed one credit; six must be reserved before work.
    assert!(session.step().is_err());
    assert_eq!(session.snapshot().iteration, 0);
    assert!(session.history().is_empty());
    let events = session.events().unwrap();
    assert!(matches!(events.pop().unwrap().kind, EventKind::Initialized));
    assert!(matches!(
        events.pop().unwrap().kind,
        EventKind::Failed {
            code: FailureCode::TelemetryExhausted
        }
    ));
    let mut session = ScfSession::new(
        fixture(0, Mixing::Linear { beta: 0.3 }, 1, 1),
        CompletionMode::RetainLast,
        ExecutionConfig {
            mode: ExecutionMode::Observed,
            routine_capacity: 0,
            ..Default::default()
        },
    )
    .unwrap();
    session.initialize().unwrap();
    session.step().unwrap();
    let mut seq = 0;
    let mut loss = 0;
    while let Some(e) = session.events().unwrap().pop() {
        assert!(e.sequence > seq);
        seq = e.sequence;
        if let EventKind::EventsDropped {
            count,
            first_missing_seq,
            last_missing_seq,
        } = e.kind
        {
            assert_eq!(count, 1);
            assert_eq!(first_missing_seq, last_missing_seq);
            loss += count;
        }
    }
    assert_eq!(loss, 1);
}
#[test]
#[ignore = "release timing only; cargo test --release -p dft-engine --test stepwise_equivalence mode_timings -- --ignored --nocapture"]
fn mode_timings() {
    let modes = [
        ExecutionMode::Classic,
        ExecutionMode::Observed,
        ExecutionMode::AgentControlled,
    ];
    let mut samples = [Vec::new(), Vec::new(), Vec::new()];
    // Warm up first, rotate order to reduce shared-host drift bias.
    for repeat in 0..111 {
        for offset in 0..3 {
            let index = (repeat + offset) % 3;
            let problem = fixture(
                0,
                Mixing::Pulay {
                    beta: 0.3,
                    history: 6,
                },
                1,
                200,
            );
            let start = Instant::now();
            drive(problem, modes[index]);
            let seconds = start.elapsed().as_secs_f64();
            if repeat >= 10 {
                samples[index].push(seconds);
            }
        }
    }
    for (mode, values) in modes.into_iter().zip(&mut samples) {
        values.sort_by(f64::total_cmp);
        println!(
            "{mode:?}: median={} s min={} s max={} s",
            values[50], values[0], values[100]
        );
    }
}

#[test]
fn divergence_is_terminal_and_observable() {
    for mode in [
        ExecutionMode::Classic,
        ExecutionMode::Observed,
        ExecutionMode::AgentControlled,
    ] {
        let mut session = ScfSession::new(
            fixture(0, Mixing::Linear { beta: 100. }, 1, 200),
            CompletionMode::Strict,
            ExecutionConfig {
                mode,
                ..Default::default()
            },
        )
        .unwrap();
        session.initialize().unwrap();
        let error = loop {
            if let Err(error) = session.step() {
                break error;
            }
        };
        assert!(matches!(error, ScfError::Diverged { .. }));
        assert!(session.result().is_none());
        let version = session.snapshot().state_version;
        assert_eq!(session.snapshot().status, StepStatus::Failed);
        assert!(session.step().is_err());
        assert_eq!(session.snapshot().state_version, version);
        if let Some(queue) = session.events() {
            let mut diverged = false;
            while let Some(e) = queue.pop() {
                diverged |= matches!(e.kind, EventKind::Diverged);
            }
            assert!(diverged);
        }
    }
}

// FNV-1a over all scalar Debug fields (timings zeroed) and full array bits.
// Explicit array traversal prevents pretty-print truncation hiding a change.
fn fingerprint(result: &ScfResult) -> u64 {
    let mut hash = 0xcbf29ce484222325u64;
    let mut feed = |bytes: &[u8]| {
        for b in bytes {
            hash = (hash ^ u64::from(*b)).wrapping_mul(0x100000001b3);
        }
    };
    feed(format!("{result:?}").as_bytes());
    for v in &result.density {
        feed(&v.to_bits().to_le_bytes());
    }
    for k in &result.k_points {
        for j in 0..k.solution.eigenvectors.cols() {
            for v in k.solution.eigenvectors.column(j) {
                feed(&v.re.to_bits().to_le_bytes());
                feed(&v.im.to_bits().to_le_bytes());
            }
        }
    }
    hash
}

#[test]
fn slow_consumer_preserves_certificates_and_reports_all_routine_loss() {
    let mut session = ScfSession::new(
        fixture(
            0,
            Mixing::Pulay {
                beta: 0.3,
                history: 6,
            },
            1,
            200,
        ),
        CompletionMode::Strict,
        ExecutionConfig {
            mode: ExecutionMode::Observed,
            routine_capacity: 1,
            ..Default::default()
        },
    )
    .unwrap();
    session.initialize().unwrap();
    while session.step().unwrap().status == StepStatus::Ready {}
    let iterations = session.snapshot().iteration as u64;
    let mut last_seq = 0;
    let mut lost = 0;
    let mut routine = 0;
    let mut certified = 0;
    let mut terminal = 0;
    while let Some(e) = session.events().unwrap().pop() {
        assert!(e.sequence > last_seq);
        last_seq = e.sequence;
        match e.kind {
            EventKind::EventsDropped { count, .. } => lost += count,
            EventKind::IterationFinished { .. } => routine += 1,
            EventKind::CertificateEvaluated { passed: true, .. } => certified += 1,
            EventKind::Converged => terminal += 1,
            _ => {}
        }
    }
    assert_eq!(routine, 1);
    assert_eq!(routine + lost, iterations);
    assert_eq!(certified, 1);
    assert_eq!(terminal, 1);
}

#[test]
fn invalid_mixer_initialization_emits_failure_without_a_record() {
    let mut session = ScfSession::new(
        fixture(0, Mixing::Linear { beta: f64::NAN }, 1, 200),
        CompletionMode::Strict,
        ExecutionConfig {
            mode: ExecutionMode::Observed,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(session.step().is_err());
    assert_eq!(session.snapshot().state_version, 0);
    assert!(session.initialize().is_err());
    assert_eq!(session.snapshot().state_version, 1);
    assert!(session.history().is_empty());
    assert!(matches!(
        session.events().unwrap().pop().unwrap().kind,
        EventKind::Failed { .. }
    ));
    assert!(session.initialize().is_err());
    assert_eq!(session.snapshot().state_version, 1);
}
