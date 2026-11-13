//! P48F diagnostics: wall clocks of the work that the per-iteration phase
//! timings (`PhaseTimings`) do not cover, so that a result file can close its
//! time budget (total = setup + phases + ladder coarse stage + certificate
//! work outside the phases + step overhead + remainder).
//!
//! Process-global, cumulative, relaxed atomics (the pattern of the Davidson
//! counters): no numerical value depends on them. One CLI run is one
//! calculation; callers that run several calculations in one process take a
//! [`snapshot`] before and after and subtract ([`PhaseClock::since`]).
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{Duration, Instant};

#[derive(Clone, Copy)]
enum Slot {
    /// Wall of `ScfSession::step` (every session, coarse ladder stage included).
    Steps = 0,
    /// Wall of the fixed-point certificate block (screen, re-solve, checks).
    Certificate = 1,
    /// Number of certificate blocks entered (screened ones included).
    CertificateCount = 2,
    /// Wall of the certificate pre-screen alone (inside `Certificate`).
    Screen = 3,
    /// Wall of the atomic-orbital start (span, H Q, Ritz); it is charged to
    /// the potential phase of the first map.
    AtomicStart = 4,
}
const SLOTS: usize = 5;
static CLOCK: [AtomicU64; SLOTS] = [const { AtomicU64::new(0) }; SLOTS];

fn add(slot: Slot, value: u64) {
    CLOCK[slot as usize].fetch_add(value, Ordering::Relaxed);
}
fn nanos(d: Duration) -> u64 {
    u64::try_from(d.as_nanos()).unwrap_or(u64::MAX)
}
pub(crate) fn add_steps(since: Instant) {
    add(Slot::Steps, nanos(since.elapsed()));
}
pub(crate) fn add_certificate(since: Instant) {
    add(Slot::Certificate, nanos(since.elapsed()));
    add(Slot::CertificateCount, 1);
}
pub(crate) fn add_screen(seconds: f64) {
    add(
        Slot::Screen,
        nanos(Duration::from_secs_f64(seconds.max(0.0))),
    );
}
pub(crate) fn add_atomic_start(since: Instant) {
    add(Slot::AtomicStart, nanos(since.elapsed()));
}

/// Cumulative clocks of this process.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PhaseClock {
    pub steps: Duration,
    pub certificate: Duration,
    pub certificate_blocks: u64,
    pub screen: Duration,
    pub atomic_start: Duration,
}
impl PhaseClock {
    /// Clocks accumulated after `earlier` was taken.
    pub fn since(&self, earlier: &PhaseClock) -> PhaseClock {
        PhaseClock {
            steps: self.steps.saturating_sub(earlier.steps),
            certificate: self.certificate.saturating_sub(earlier.certificate),
            certificate_blocks: self
                .certificate_blocks
                .saturating_sub(earlier.certificate_blocks),
            screen: self.screen.saturating_sub(earlier.screen),
            atomic_start: self.atomic_start.saturating_sub(earlier.atomic_start),
        }
    }
}

/// Current cumulative clocks.
pub fn snapshot() -> PhaseClock {
    let get = |s: Slot| CLOCK[s as usize].load(Ordering::Relaxed);
    PhaseClock {
        steps: Duration::from_nanos(get(Slot::Steps)),
        certificate: Duration::from_nanos(get(Slot::Certificate)),
        certificate_blocks: get(Slot::CertificateCount),
        screen: Duration::from_nanos(get(Slot::Screen)),
        atomic_start: Duration::from_nanos(get(Slot::AtomicStart)),
    }
}
