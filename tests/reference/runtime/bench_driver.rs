//! Standalone release benchmark adapter; no instrumentation in engine kernels.
#[cfg(not(pre_refactor))]
#[allow(dead_code)]
#[path = "../../../crates/dft-engine/tests/common/checkpoint.rs"]
mod checkpoint_fixture;
#[allow(dead_code)]
#[path = "../../../crates/dft-engine/tests/common/mod.rs"]
mod common;
use dft_engine::*;
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::time::{Duration, Instant};
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

fn digest(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}
fn identity(mut out: ScfOutcome) -> Value {
    let r = &mut out.result;
    let energy = r.total_hartree();
    let iterations = r.history.len();
    let certificates = r
        .history
        .iter()
        .filter(|h| h.fixed_point_certificate.is_some())
        .count();
    r.setup_time = Duration::ZERO;
    r.total_time = Duration::ZERO;
    for h in &mut r.history {
        h.timings = PhaseTimings::default();
        if let Some(c) = &mut h.fixed_point_certificate {
            c.timings = PhaseTimings::default();
        }
    }
    let mut bits = Vec::new();
    bits.extend(energy.to_bits().to_le_bytes());
    for v in &r.density {
        bits.extend(v.to_bits().to_le_bytes());
    }
    for k in &r.k_points {
        for v in &k.solution.eigenvalues {
            bits.extend(v.to_bits().to_le_bytes());
        }
        for col in 0..k.solution.eigenvectors.cols() {
            for v in k.solution.eigenvectors.column(col) {
                bits.extend(v.re.to_bits().to_le_bytes());
                bits.extend(v.im.to_bits().to_le_bytes());
            }
        }
    }
    json!({"energy_hartree":energy,"iterations":iterations,"certificates":certificates,
        "converged":out.converged,"arrays_sha256":digest(&bits),
        "normalized_result_sha256":digest(format!("{r:?}").as_bytes())})
}
#[cfg(not(pre_refactor))]
#[derive(Default)]
struct Sink {
    events: u64,
    dropped: u64,
    high_water: usize,
    missing: Vec<(u64, u64)>,
    latency_ns: Vec<u128>,
    serialization_ns: u128,
    bytes: usize,
}
#[cfg(not(pre_refactor))]
fn drain(r: &mut runtime::Runtime, sink: &mut Sink, instrument: bool) {
    use events::EventKind;
    if let Some(q) = r.events() {
        sink.high_water = sink.high_water.max(q.len());
        while let Some(e) = q.pop() {
            sink.events += 1;
            if let EventKind::EventsDropped {
                count,
                first_missing_seq,
                last_missing_seq,
            } = e.kind
            {
                sink.dropped += count;
                sink.missing.push((first_missing_seq, last_missing_seq));
            }
            if instrument {
                if let Ok(elapsed) = e.timestamp.elapsed() {
                    sink.latency_ns.push(elapsed.as_nanos());
                }
                let t = Instant::now();
                let bytes = serde_json::to_vec(&dft_cli::telemetry::event_json(&e)).unwrap();
                sink.serialization_ns += t.elapsed().as_nanos();
                sink.bytes += bytes.len();
                std::hint::black_box(bytes);
            }
            std::hint::black_box(e);
        }
    }
}
#[cfg(not(pre_refactor))]
fn drive(
    p: ScfProblem,
    mode: &str,
    probe: &str,
    kind: usize,
    threads: usize,
    budget: usize,
) -> (ScfOutcome, Value) {
    use runtime::*;
    // Normal modes use a typed scalar sink without probe-report allocations.
    if probe == "normal" {
        if mode == "classic" {
            return (p.run_outcome().unwrap(), json!({}));
        }
        let mut runtime = Runtime::new(
            ScfSession::new(
                p,
                CompletionMode::RetainLast,
                ExecutionConfig {
                    mode: if mode == "observed" {
                        ExecutionMode::Observed
                    } else {
                        ExecutionMode::AgentControlled
                    },
                    ..Default::default()
                },
            )
            .unwrap(),
            RuntimeConfig::default(),
        )
        .unwrap();
        let mut sink = Sink::default();
        runtime.initialize().unwrap();
        while !runtime.state().terminal() {
            runtime.step().unwrap();
            drain(&mut runtime, &mut sink, false);
        }
        drain(&mut runtime, &mut sink, false);
        std::hint::black_box(sink);
        return (runtime.into_result().outcome.unwrap(), json!({}));
    }
    let mut stats =
        json!({"events":0,"dropped":0,"queue_high_water":0,"missing_intervals":[],"steps_ns":[]});
    // Classic uses the production conventional entry point, with no observer.
    if mode == "classic" {
        return (p.run_outcome().unwrap(), stats);
    }
    let mut r = Runtime::new(
        ScfSession::new(
            p,
            CompletionMode::RetainLast,
            ExecutionConfig {
                mode: if mode == "observed" {
                    ExecutionMode::Observed
                } else {
                    ExecutionMode::AgentControlled
                },
                routine_capacity: if probe == "saturation" { 1 } else { 256 },
                ..Default::default()
            },
        )
        .unwrap(),
        RuntimeConfig::default(),
    )
    .unwrap();
    let mut sink = Sink::default();
    r.initialize().unwrap();
    if probe == "control" {
        for (i, action) in [ControlAction::Pause, ControlAction::Resume]
            .into_iter()
            .enumerate()
        {
            let req = ControlRequest {
                schema_version: 1,
                calc_id: "scf".into(),
                id: format!("request-{i}"),
                observed_state_version: r.snapshot().state_version,
                observed_policy_revision: r.policy_revision(),
                action,
                preconditions: vec![],
                expiry: Expiry {
                    max_iteration: Some(200),
                    deadline: None,
                },
                stale_policy: StalePolicy::Reject,
                reason: "benchmark lifecycle latency".into(),
            };
            let authority = Authority::Approved(Approval {
                id: "host".into(),
                approver: "host".into(),
                expires_at: std::time::SystemTime::now() + Duration::from_secs(60),
                allow_revalidation: false,
            });
            let t = Instant::now();
            r.submit(req, "host".into(), authority).unwrap();
            let submit = t.elapsed().as_nanos();
            r.safe_point().unwrap();
            stats[format!("control_{i}")] = json!({"submit_ns":submit,"receipt_to_applied_ns":t.elapsed().as_nanos(),"boundaries":0,"approval_wait_ns":0});
        }
        assert_eq!(r.audit().len(), 4);
    }
    while !r.state().terminal() {
        let t = (probe != "normal").then(Instant::now);
        r.step().unwrap();
        let step_ns = t.map_or(0, |v| v.elapsed().as_nanos());
        // Probe instrumentation only; normal runs avoid growing a trace per step.
        if probe != "normal" {
            stats["steps_ns"]
                .as_array_mut()
                .unwrap()
                .push(json!(step_ns));
        }
        if probe != "saturation" {
            drain(&mut r, &mut sink, probe == "events");
        }
        if probe == "checkpoint" && r.snapshot().iteration == 1 {
            let (model, mut execution) = checkpoint_fixture::identity(kind, threads);
            execution.build_sha256 =
                digest(&std::fs::read(std::env::current_exe().unwrap()).unwrap());
            let t = Instant::now();
            let cp = r.checkpoint(model, execution).unwrap();
            stats["checkpoint_capture_ns"] = json!(t.elapsed().as_nanos());
            let t = Instant::now();
            let bytes = dft_cli::checkpoint_io::encode(&cp).unwrap();
            stats["checkpoint_encode_ns"] = json!(t.elapsed().as_nanos());
            stats["checkpoint_bytes"] = json!(bytes.len());
            let path = std::env::temp_dir()
                .join(format!("runtime-bench-{}.checkpoint", std::process::id()));
            let t = Instant::now();
            dft_cli::checkpoint_io::write(&path, &cp).unwrap();
            stats["checkpoint_sync_persist_ns"] = json!(t.elapsed().as_nanos());
            let restored = dft_cli::checkpoint_io::read(&path).unwrap();
            std::fs::remove_file(path).unwrap();
            assert_eq!(dft_cli::checkpoint_io::encode(&restored).unwrap(), bytes);
            let t = Instant::now();
            r = Runtime::restore_exact(
                fixture(kind, Mixing::default(), threads, budget),
                &restored,
                &restored.model,
                &restored.execution,
            )
            .unwrap();
            stats["checkpoint_restore_including_setup_ns"] = json!(t.elapsed().as_nanos());
        }
    }
    drain(&mut r, &mut sink, probe == "events");
    if probe == "events" {
        use events::{EventKind, ScfEvent};
        let calc_id: std::sync::Arc<str> = "microbench".into();
        let count = 10000;
        let start = Instant::now();
        for i in 0..count {
            std::hint::black_box(ScfEvent {
                schema_version: "1.0",
                calc_id: calc_id.clone(),
                sequence: i,
                timestamp: std::time::SystemTime::now(),
                iteration: i as usize,
                state_version: i,
                kind: EventKind::IterationFinished {
                    energy_hartree: 0.,
                    density_residual: 0.,
                    charge_error_relative: 0.,
                    mixing_applied: false,
                },
            });
        }
        stats["typed_construction_ns_per_event"] =
            json!(start.elapsed().as_nanos() as f64 / count as f64);
        let start = Instant::now();
        for _ in 0..count {
            std::hint::black_box(Instant::now().elapsed());
        }
        stats["timer_pair_ns"] = json!(start.elapsed().as_nanos() as f64 / count as f64);
        stats["construction_count"] = json!(count);
        stats["event_size_bytes"] = json!(std::mem::size_of::<ScfEvent>());
    }
    stats["events"] = json!(sink.events);
    stats["dropped"] = json!(sink.dropped);
    stats["queue_high_water"] = json!(sink.high_water);
    stats["missing_intervals"] = json!(sink.missing);
    stats["delivery_wall_latency_ns"] = json!(sink.latency_ns);
    stats["serialization_ns"] = json!(sink.serialization_ns);
    stats["serialized_bytes"] = json!(sink.bytes);
    (r.into_result().outcome.unwrap(), stats)
}
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let kind: usize = args[2].parse().unwrap();
    let budget: usize = args[3].parse().unwrap();
    let threads: usize = args[4].parse().unwrap();
    let start = Instant::now();
    let p = fixture(kind, Mixing::default(), threads, budget);
    let setup_ns = start.elapsed().as_nanos();
    let run = Instant::now();
    #[cfg(pre_refactor)]
    let (out, stats) = (p.run_outcome().unwrap(), json!({}));
    #[cfg(not(pre_refactor))]
    let (out, stats) = drive(p, &args[1], &args[5], kind, threads, budget);
    let scf_ns = run.elapsed().as_nanos();
    let total_ns = start.elapsed().as_nanos();
    // Fingerprint and stdout are intentionally outside the in-process timer.
    println!(
        "{}",
        json!({"setup_ns":setup_ns,"scf_ns":scf_ns,"total_ns":total_ns,
        "identity":identity(out),"probe":stats})
    );
}
