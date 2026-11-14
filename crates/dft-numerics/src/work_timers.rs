//! Work timers of the operator application (work order beat-QE P0.2), active
//! only under `BAANDI_DAVIDSON_STAGES=on` (the Davidson stage-timer switch).
//!
//! * Wave-function FFT time split by kind: the paired orbital transforms of the
//!   real Gamma `h_apply` (inverse, potential multiply, forward), the density
//!   transforms, and other wave transforms (complex-path local applies, the
//!   folded odd tail). Each kind carries its transform and column counts, so the
//!   FFT time can be set beside QE's `fftw` timer (wave transforms only).
//! * The real Gamma `h_apply` parts: whole application, local part (pack, pair
//!   FFTs, unpack), nonlocal projection P^T X and nonlocal back product
//!   (coupling, P W, Y += image).
//! * Every slot is kept per run stage ([`RunStage`]: fine, ladder coarse,
//!   verification), selected by a process-wide index that the engine sets at
//!   serial boundaries.
//!
//! Cost model: when off, each instrumented site does one relaxed atomic load
//! and a predicted branch; no clock is read and no arithmetic changes (the
//! timed and untimed paths run the same operations in the same order). When on,
//! clocks are read per pair transform (per slab of the mirror-slab plan), per
//! density band or pair, and per block application; the sums go to relaxed
//! atomics once per call. No allocation, no I/O. Times are summed over worker
//! threads (thread-seconds); counts are independent of the thread count.
use std::sync::atomic::{AtomicU8, AtomicU64, AtomicUsize, Ordering::Relaxed};
use std::time::Instant;

/// Run stage of the work being timed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunStage {
    /// Full-cutoff SCF iterations (the default stage).
    Fine = 0,
    /// Cutoff-ladder coarse stage.
    Coarse = 1,
    /// Fixed-point certificate / verification solves.
    Verification = 2,
}
/// Number of run stages.
pub const RUN_STAGES: usize = 3;
/// Labels in index order of [`RunStage`].
pub const RUN_STAGE_LABELS: [&str; RUN_STAGES] = ["fine", "ladder_coarse", "verification"];

static RUN_STAGE: AtomicUsize = AtomicUsize::new(RunStage::Fine as usize);

/// Index of the current run stage (the default is `Fine`).
#[inline]
pub fn run_stage() -> usize {
    RUN_STAGE.load(Relaxed)
}
/// Whether the current run stage is `Fine`.
pub fn in_fine_stage() -> bool {
    run_stage() == RunStage::Fine as usize
}
/// Restores the previous run stage on drop.
#[must_use]
pub struct RunStageGuard {
    previous: usize,
}
impl Drop for RunStageGuard {
    fn drop(&mut self) {
        RUN_STAGE.store(self.previous, Relaxed);
    }
}
/// Enter `stage` until the guard is dropped. Call at serial boundaries only
/// (between solves), never from a worker.
pub fn enter(stage: RunStage) -> RunStageGuard {
    RunStageGuard {
        previous: RUN_STAGE.swap(stage as usize, Relaxed),
    }
}

// 0 = not read yet, 1 = off, 2 = on.
static STATE: AtomicU8 = AtomicU8::new(0);

/// Whether the work timers are on (`BAANDI_DAVIDSON_STAGES=on`, read once).
#[inline]
pub fn enabled() -> bool {
    match STATE.load(Relaxed) {
        2 => true,
        1 => false,
        _ => init(),
    }
}
#[cold]
fn init() -> bool {
    let on = crate::switches::on(crate::switches::Id::DAVIDSON_STAGES);
    let _ = STATE.compare_exchange(0, if on { 2 } else { 1 }, Relaxed, Relaxed);
    STATE.load(Relaxed) == 2
}
/// Test hook: force the timers on or off for this process.
pub fn set_enabled(on: bool) {
    STATE.store(if on { 2 } else { 1 }, Relaxed);
}

/// Counter slots (nanoseconds or counts).
#[derive(Debug, Clone, Copy)]
pub enum Slot {
    /// Paired orbital transforms in `h_apply`: inverse (scatter, three axes).
    PairInverseNs = 0,
    /// Paired orbital transforms: local potential multiply on the grid.
    PairPotentialNs = 1,
    /// Paired orbital transforms: forward (three axes, gather, 1/N).
    PairForwardNs = 2,
    /// Paired applications of a plan without the split (whole call).
    PairUnsplitNs = 3,
    /// Pair applications (one inverse and one forward transform each).
    PairTransforms = 4,
    /// Orbital columns carried by the pair applications (2 per pair, 1 for a
    /// zero-partner tail).
    PairColumns = 5,
    /// Density transforms (inverse transform and |psi|^2 accumulation).
    DensityNs = 6,
    /// Density inverse transforms.
    DensityTransforms = 7,
    /// Bands carried by the density transforms.
    DensityColumns = 8,
    /// Other wave transforms (complex local apply, folded odd tail): whole
    /// call (one inverse, potential multiply, one forward).
    OtherNs = 9,
    /// Other applications (one inverse and one forward transform each).
    OtherTransforms = 10,
    /// Columns carried by other applications.
    OtherColumns = 11,
    /// Real Gamma `h_apply` block applications: whole call.
    ApplyNs = 12,
    /// Real Gamma `h_apply` block applications.
    ApplyCalls = 13,
    /// Columns of those applications.
    ApplyColumns = 14,
    /// Local part (pack, pair FFTs and multiply, unpack), wall of the caller.
    LocalNs = 15,
    /// Nonlocal projection P^T X.
    NonlocalProjectionNs = 16,
    /// Nonlocal back product (coupling D, P W, Y += image, output scan).
    NonlocalBackNs = 17,
}
/// Number of slots.
pub const SLOTS: usize = 18;

static COUNTS: [[AtomicU64; SLOTS]; RUN_STAGES] =
    [const { [const { AtomicU64::new(0) }; SLOTS] }; RUN_STAGES];

/// Add `value` to `slot` of the current run stage.
#[inline]
pub fn add(slot: Slot, value: u64) {
    COUNTS[run_stage()][slot as usize].fetch_add(value, Relaxed);
}
/// `Some(now)` when the timers are on.
#[inline]
pub fn start() -> Option<Instant> {
    enabled().then(Instant::now)
}
/// Add the time since `start` to `slot` (no-op for `None`).
#[inline]
pub fn stop(start: Option<Instant>, slot: Slot) {
    if let Some(t) = start {
        add(slot, nanos(t));
    }
}
/// Nanoseconds since `t`, saturating.
pub fn nanos(t: Instant) -> u64 {
    u64::try_from(t.elapsed().as_nanos()).unwrap_or(u64::MAX)
}

/// Phase of one fused local application.
#[derive(Debug, Clone, Copy)]
pub enum Phase {
    Inverse = 0,
    Potential = 1,
    Forward = 2,
}
/// Lap recorder of a fused local application; `NoLaps` compiles to nothing.
pub trait Laps {
    /// Charge the time since the previous lap to `phase`.
    fn lap(&mut self, phase: Phase);
}
/// The untimed recorder (all methods empty).
pub struct NoLaps;
impl Laps for NoLaps {
    #[inline(always)]
    fn lap(&mut self, _: Phase) {}
}
/// Timed recorder: per-phase nanoseconds of one call, accumulated locally and
/// added to the shared counters once by [`LocalLaps::record_pair`].
pub struct LocalLaps {
    last: Instant,
    nanos: [u64; 3],
    unsplit: u64,
}
impl LocalLaps {
    pub fn start() -> Self {
        Self {
            last: Instant::now(),
            nanos: [0; 3],
            unsplit: 0,
        }
    }
    /// Charge the whole time since start to the unsplit slot (plans without
    /// a phase split).
    pub fn lap_unsplit(&mut self) {
        let now = Instant::now();
        self.unsplit += u64::try_from((now - self.last).as_nanos()).unwrap_or(u64::MAX);
        self.last = now;
    }
    /// Add this call's phases to the pair slots, with one pair and `columns`.
    pub fn record_pair(&self, columns: u64) {
        add(Slot::PairInverseNs, self.nanos[0]);
        add(Slot::PairPotentialNs, self.nanos[1]);
        add(Slot::PairForwardNs, self.nanos[2]);
        add(Slot::PairUnsplitNs, self.unsplit);
        add(Slot::PairTransforms, 1);
        add(Slot::PairColumns, columns);
    }
}
impl Laps for LocalLaps {
    #[inline]
    fn lap(&mut self, phase: Phase) {
        let now = Instant::now();
        self.nanos[phase as usize] +=
            u64::try_from((now - self.last).as_nanos()).unwrap_or(u64::MAX);
        self.last = now;
    }
}

/// Snapshot of one run stage, seconds and counts.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct WorkTimes {
    pub pair_inverse: f64,
    pub pair_potential: f64,
    pub pair_forward: f64,
    pub pair_unsplit: f64,
    pub pair_transforms: u64,
    pub pair_columns: u64,
    pub density: f64,
    pub density_transforms: u64,
    pub density_columns: u64,
    pub other: f64,
    pub other_transforms: u64,
    pub other_columns: u64,
    pub apply: f64,
    pub apply_calls: u64,
    pub apply_columns: u64,
    pub local: f64,
    pub nonlocal_projection: f64,
    pub nonlocal_back: f64,
}
impl WorkTimes {
    /// Element-wise sum.
    pub fn plus(&self, o: &WorkTimes) -> WorkTimes {
        WorkTimes {
            pair_inverse: self.pair_inverse + o.pair_inverse,
            pair_potential: self.pair_potential + o.pair_potential,
            pair_forward: self.pair_forward + o.pair_forward,
            pair_unsplit: self.pair_unsplit + o.pair_unsplit,
            pair_transforms: self.pair_transforms + o.pair_transforms,
            pair_columns: self.pair_columns + o.pair_columns,
            density: self.density + o.density,
            density_transforms: self.density_transforms + o.density_transforms,
            density_columns: self.density_columns + o.density_columns,
            other: self.other + o.other,
            other_transforms: self.other_transforms + o.other_transforms,
            other_columns: self.other_columns + o.other_columns,
            apply: self.apply + o.apply,
            apply_calls: self.apply_calls + o.apply_calls,
            apply_columns: self.apply_columns + o.apply_columns,
            local: self.local + o.local,
            nonlocal_projection: self.nonlocal_projection + o.nonlocal_projection,
            nonlocal_back: self.nonlocal_back + o.nonlocal_back,
        }
    }
}

/// Cumulative counters of this process per run stage (index order of
/// [`RunStage`]). All zero unless the timers were on.
pub fn snapshot() -> [WorkTimes; RUN_STAGES] {
    std::array::from_fn(|stage| {
        let g = |s: Slot| COUNTS[stage][s as usize].load(Relaxed);
        let s = |slot: Slot| g(slot) as f64 * 1e-9;
        WorkTimes {
            pair_inverse: s(Slot::PairInverseNs),
            pair_potential: s(Slot::PairPotentialNs),
            pair_forward: s(Slot::PairForwardNs),
            pair_unsplit: s(Slot::PairUnsplitNs),
            pair_transforms: g(Slot::PairTransforms),
            pair_columns: g(Slot::PairColumns),
            density: s(Slot::DensityNs),
            density_transforms: g(Slot::DensityTransforms),
            density_columns: g(Slot::DensityColumns),
            other: s(Slot::OtherNs),
            other_transforms: g(Slot::OtherTransforms),
            other_columns: g(Slot::OtherColumns),
            apply: s(Slot::ApplyNs),
            apply_calls: g(Slot::ApplyCalls),
            apply_columns: g(Slot::ApplyColumns),
            local: s(Slot::LocalNs),
            nonlocal_projection: s(Slot::NonlocalProjectionNs),
            nonlocal_back: s(Slot::NonlocalBackNs),
        }
    })
}
