//! G-SCHED-OVERHEAD: measure schedule-only allocation deltas around real
//! single-thread SCF steps. This binary has one test, so no concurrent test
//! harness can contaminate the process-wide allocator counters. Fixture setup,
//! assertions and collection are outside the measured regions.
#[allow(dead_code)]
#[path = "common/mod.rs"]
mod common;
use dft_engine::*;
use dft_geometry::{Atom, Cell};
use dft_numerics::switches::RunConfig;
use std::sync::atomic::Ordering;

// Test-only System allocator adapter (as in pair_fused_plan.rs).
#[allow(unsafe_code)]
mod probe {
    use std::{
        alloc::{GlobalAlloc, Layout, System},
        sync::atomic::{AtomicBool, AtomicUsize, Ordering},
    };
    pub static TRACK: AtomicBool = AtomicBool::new(false);
    pub static COUNT: AtomicUsize = AtomicUsize::new(0);
    pub struct Allocator;
    fn record() {
        if TRACK.load(Ordering::Relaxed) {
            COUNT.fetch_add(1, Ordering::Relaxed);
        }
    }
    // SAFETY: every pointer and layout is forwarded unchanged to the System
    // allocator; the counters are atomics and never allocate.
    unsafe impl GlobalAlloc for Allocator {
        unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
            record();
            // SAFETY: the caller supplies a valid layout.
            unsafe { System.alloc(layout) }
        }
        unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
            record();
            // SAFETY: the caller supplies a valid layout.
            unsafe { System.alloc_zeroed(layout) }
        }
        unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
            // SAFETY: pointer and layout come from this forwarding allocator.
            unsafe { System.dealloc(ptr, layout) }
        }
        unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, size: usize) -> *mut u8 {
            record();
            // SAFETY: pointer, layout and size keep the GlobalAlloc contract.
            unsafe { System.realloc(ptr, layout, size) }
        }
    }
}
#[global_allocator]
static ALLOCATOR: probe::Allocator = probe::Allocator;

fn measured<T>(f: impl FnOnce() -> T) -> (T, usize) {
    probe::COUNT.store(0, Ordering::SeqCst);
    probe::TRACK.store(true, Ordering::SeqCst);
    let value = f();
    probe::TRACK.store(false, Ordering::SeqCst);
    (value, probe::COUNT.load(Ordering::SeqCst))
}

fn drive(mode: &str) -> (ScfResult, Vec<usize>) {
    let run = RunConfig::from_pairs(&[("verify_schedule", mode)]).unwrap();
    let _scope = run.enter();
    let potential = common::potential(2.0, 1.0, true, 401);
    let cell = Cell::from_user_rows_bohr([[4., 0., 0.], [0.3, 4.2, 0.], [0., 0.2, 4.4]]).unwrap();
    let problem = ScfProblem::new(
        &cell,
        &[Atom::new("Xx", [0.13, 0.21, 0.32])],
        &[SpeciesPotential {
            label: "Xx",
            potential: &potential,
        }],
        ScfOptions {
            cutoff_hartree: 2.5,
            bands: Some(3),
            threads: 1,
            standard: Some(StandardProfile { threshold_ha: 5e-9 }),
            verification: Verification::FixedPoint,
            ..ScfOptions::strict_reference()
        },
    )
    .unwrap();
    let mut session =
        ScfSession::new(problem, CompletionMode::Strict, ExecutionConfig::default()).unwrap();
    session.initialize().unwrap();
    let mut counts = Vec::with_capacity(100);
    loop {
        let (step, count) = measured(|| session.step().unwrap());
        counts.push(count);
        if step.status == StepStatus::Converged {
            break;
        }
        assert_eq!(step.status, StepStatus::Ready);
    }
    (session.into_result().unwrap().result, counts)
}

#[test]
fn no_failed_attempt_adds_zero_allocations_per_ordinary_map() {
    // Prove the counter observes allocations; prevent allocation elision.
    let (allocation, count) = measured(|| std::hint::black_box(vec![0u8; 123]));
    assert_eq!(count, 1);
    drop(allocation);
    // Warm all library/thread-local paths before the paired measurement.
    drop(drive("off"));
    let normalize = |r: &mut ScfResult| {
        r.setup_time = std::time::Duration::ZERO;
        r.total_time = std::time::Duration::ZERO;
        for h in &mut r.history {
            h.timings = PhaseTimings::default();
            if let Some(c) = &mut h.fixed_point_certificate {
                c.timings = PhaseTimings::default();
            }
        }
        for a in &mut r.verification.as_mut().unwrap().attempts {
            a.wall = std::time::Duration::ZERO;
        }
        format!("{r:?}")
    };
    let (mut off, baseline) = drive("off");
    let expected = normalize(&mut off);
    for mode in ["on", "trace"] {
        let (mut result, counts) = drive(mode);
        assert_eq!(baseline.len(), counts.len());
        let report = result.verification.as_mut().unwrap();
        assert_eq!(report.attempts.len(), 1);
        assert!(report.attempts[0].passed);
        let schedule = report.schedule.take().unwrap();
        assert_eq!(schedule.prediction_evaluations, 0);
        assert_eq!(schedule.records.len(), 1);
        assert_eq!(schedule.records.capacity(), 1);
        for (i, (a, b)) in baseline.iter().zip(&counts).enumerate() {
            // The final step also constructs the public result's diagnostic
            // Vec. All preceding ordinary maps have exactly zero added allocations.
            let output_allocations = if i + 1 == counts.len() { 1 } else { 0 };
            println!(
                "{mode} map {}: off={a}, schedule={b}, extra={}",
                i + 1,
                b - a
            );
            assert_eq!(*b, a + output_allocations, "{mode} map {}", i + 1);
        }
        // Compare complete scientific state, including orbitals, history and
        // attempts, excluding wall clocks and the optional diagnostic journal.
        assert_eq!(expected, normalize(&mut result));
    }
}
