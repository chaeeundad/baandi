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
