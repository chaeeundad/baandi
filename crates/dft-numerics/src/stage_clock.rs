//! F21 (review 2026-11-16): time mode of the MPI stage and collective counters
//! (`BAANDI_STAGE_TIMERS=on|off`, default on = the behaviour before the
//! switch). With `off` the counters run count-only: every stage and collective
//! still counts its calls (and collectives their payload bytes), but no clock
//! is read (`start` returns `None`, no `Instant::now` or `elapsed`), so an
//! acceptance run can be timed without the instrumentation and the same binary
//! measures the instrumentation overhead (on/off A/B). Stage seconds then read
//! 0. The serial work timers keep their own switches (`BAANDI_DAVIDSON_STAGES`,
//! work_timers).
use std::sync::atomic::{AtomicU8, Ordering::Relaxed};
use std::time::Instant;

/// 0: from the switch; 1: forced on; 2: forced off (tests).
static FORCED: AtomicU8 = AtomicU8::new(0);

/// Whether stage and collective times are measured.
#[inline]
pub fn on() -> bool {
    match FORCED.load(Relaxed) {
        1 => true,
        2 => false,
        _ => {
            static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
            *ON.get_or_init(|| crate::switches::on(crate::switches::Id::STAGE_TIMERS))
        }
    }
}

/// Tests: force the time mode for the whole process (`None`: the switch).
pub fn force(mode: Option<bool>) {
    FORCED.store(
        match mode {
            None => 0,
            Some(true) => 1,
            Some(false) => 2,
        },
        Relaxed,
    );
}

/// A stage clock: `Some(now)` in time mode, `None` count-only.
#[inline]
pub fn start() -> Option<Instant> {
    on().then(Instant::now)
}

/// A clock read for another purpose, kept for a stage only in time mode.
#[inline]
pub fn keep(clock: Instant) -> Option<Instant> {
    on().then_some(clock)
}

/// Nanoseconds since `start` (0 count-only).
#[inline]
pub fn elapsed_ns(start: Option<Instant>) -> u64 {
    start.map_or(0, |s| s.elapsed().as_nanos() as u64)
}

/// Per-thread mirror of the stage, collective and operator counters
/// (ticket BREAKTHROUGH-0-INSTRUMENT, the exclusive phase profile of
/// `dft_engine::profile`). Every counted site adds its call (and, in time mode,
/// its nanoseconds) to a fixed slot of the calling thread when the mirror is
/// enabled (one relaxed load otherwise). The engine takes snapshots of the
/// calling thread at phase boundaries and attributes the deltas to the phase in
/// progress; with simulated ranks (threads of one process) every rank thread
/// therefore has its own counts, which the process-wide atomics cannot give.
/// No allocation, no I/O; nothing numerical reads these values.
pub mod local {
    use std::cell::Cell;
    use std::sync::atomic::{AtomicBool, Ordering::Relaxed};

    /// Distributed solver stages (`dft_solvers::distributed_davidson::STAGES` order).
    pub const DIST_STAGE: usize = 0;
    /// Engine M2 stages (`dft_engine::m2` stage order).
    pub const M2_STAGE: usize = 16;
    /// Collectives: allgatherv, tree_sum, bcast, alltoall, barrier.
    pub const COLLECTIVE: usize = 32;
    /// Distributed operator counts (calls only): applied columns, potential
    /// change columns, forward pair transforms, inverse pair transforms, kept
    /// plane pairs.
    pub const OP_COUNT: usize = 40;
    /// Serial Davidson stages (`dft_solvers::davidson_stages::Stage` order).
    pub const SERIAL_STAGE: usize = 48;
    /// Number of slots.
    pub const SLOTS: usize = 64;

    static ENABLED: AtomicBool = AtomicBool::new(false);
    thread_local! {
        static NANOS: [Cell<u64>; SLOTS] = const { [const { Cell::new(0) }; SLOTS] };
        static CALLS: [Cell<u64>; SLOTS] = const { [const { Cell::new(0) }; SLOTS] };
    }

    /// Turn the mirror on for the process (the engine does it when its phase
    /// profile is requested). Off by default: the sites cost one relaxed load.
    pub fn enable(on: bool) {
        ENABLED.store(on, Relaxed);
    }
    #[inline]
    pub fn enabled() -> bool {
        ENABLED.load(Relaxed)
    }
    /// Add one call and `nanos` to `slot` on this thread.
    #[inline]
    pub fn add(slot: usize, nanos: u64) {
        if ENABLED.load(Relaxed) {
            add_n(slot, nanos, 1);
        }
    }
    /// Add `calls` calls and `nanos` to `slot` on this thread.
    #[inline]
    pub fn add_n(slot: usize, nanos: u64, calls: u64) {
        if ENABLED.load(Relaxed) && slot < SLOTS {
            NANOS.with(|n| n[slot].set(n[slot].get().wrapping_add(nanos)));
            CALLS.with(|c| c[slot].set(c[slot].get().wrapping_add(calls)));
        }
    }

    /// Cumulative values of the calling thread.
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub struct Snapshot {
        pub nanos: [u64; SLOTS],
        pub calls: [u64; SLOTS],
        pub fft: crate::real::fft_profile::Snapshot,
    }
    impl Default for Snapshot {
        fn default() -> Self {
            Self {
                nanos: [0; SLOTS],
                calls: [0; SLOTS],
                fft: crate::real::fft_profile::Snapshot::default(),
            }
        }
    }
    impl Snapshot {
        /// Values accumulated after `earlier` was taken.
        pub fn since(&self, earlier: &Snapshot) -> Snapshot {
            let mut out = Snapshot::default();
            for i in 0..SLOTS {
                out.nanos[i] = self.nanos[i].wrapping_sub(earlier.nanos[i]);
                out.calls[i] = self.calls[i].wrapping_sub(earlier.calls[i]);
            }
            out.fft = self.fft.since(earlier.fft);
            out
        }
    }
    /// The calling thread's cumulative counters.
    pub fn snapshot() -> Snapshot {
        let mut out = Snapshot::default();
        if crate::switches::RunConfig::current().on(crate::switches::Id::FFT_PROFILE) {
            out.fft = crate::real::fft_profile::snapshot();
        }
        NANOS.with(|n| {
            for (o, c) in out.nanos.iter_mut().zip(n.iter()) {
                *o = c.get();
            }
        });
        CALLS.with(|n| {
            for (o, c) in out.calls.iter_mut().zip(n.iter()) {
                *o = c.get();
            }
        });
        out
    }
}
