//! Research switches and process-wide counters of the serial real Davidson.
//!
//! * Stage timers (`BAANDI_DAVIDSON_STAGES=on`, default off): wall time per
//!   Davidson iteration bucketed by correction width, and the dense-algebra
//!   stages. Clock reads and relaxed atomic adds only; no I/O, no allocation.
//!   When off, no clock is read. Timers never influence the arithmetic.
//! * B1 (`BAANDI_B1=standard|strict|on`, default off): Ritz vectors, H images and
//!   residual vectors only for the roots selected for correction in narrow
//!   iterations (proposal B1, ledger row B1-SELECTED-ROOTS). `standard` acts on
//!   Ritz-change solves only, `strict` on residual solves only, `on` on both.
//!   `BAANDI_B1_SKETCH=count` replaces the dense strict sketch by a CountSketch.
//! * B5 (`BAANDI_B5=on`, default off): when exactly one column was appended,
//!   the reduced eigenproblem is updated from the previous decomposition by the
//!   secular equation (module `bordered`), with a full solve on any doubt and
//!   after `bordered::MAX_STREAK` updates in a row.
//!
//! * Algebra-stream switches (`BAANDI_RITZ_FUSED`, `BAANDI_NARROW_PROJ`,
//!   `BAANDI_NARROW_REDUCED`, `BAANDI_B5=par`; default off): counters in
//!   [`algebra_counters`] show how often each one acted.
//!
//! The counters are cumulative for the process (all solves, all threads) and are
//! read once at the reporting boundary through [`stage_seconds`] and
//! [`b1_counters`].
use std::sync::atomic::{AtomicU64, Ordering::Relaxed};
use std::time::Instant;

/// B1 selected-root products mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum B1Mode {
    /// Every iteration rotates all bands (the default path, bit-identical).
    #[default]
    Off,
    /// Ritz-change (standard-profile) solves only; exact flags, no sketch.
    Standard,
    /// Residual (strict-profile) solves only; sketched flags, explicit
    /// residuals for every root that is not clearly converged.
    Strict,
    /// Both.
    On,
}
impl B1Mode {
    pub fn from_env() -> Self {
        match dft_numerics::switches::choice(dft_numerics::switches::Id::B1) {
            "standard" => Self::Standard,
            "strict" => Self::Strict,
            "off" => Self::Off,
            // Default ON with the count sketch since the adoption rule (dev08
            // B1-PART-B: work identical to off in all 240 runs, si32 strict -4.2 / -2.0
            // percent at T=1 / T=8, si32 standard -1.7 / -0.8; si16 strict +1.1 / +0.9
            // percent is unresolved at 4 and 3 of 10 seeds). `BAANDI_B1=off` disables.
            _ => Self::On,
        }
    }
    /// Whether the mode acts on a solve with this convergence semantics.
    pub fn acts_on(self, ritz: bool) -> bool {
        match self {
            Self::Off => false,
            Self::Standard => ritz,
            Self::Strict => !ritz,
            Self::On => true,
        }
    }
}

/// Form of the B1 strict sketch (`BAANDI_B1_SKETCH`, default `dense`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SketchKind {
    /// Dense sign matrix, entries +-1/sqrt(s) (proposal B1): s N multiply-adds
    /// per sketched column.
    #[default]
    Dense,
    /// CountSketch (one +-1 per coordinate in a seeded bucket): N adds per
    /// sketched column; same variance bound of ||S x||^2.
    Count,
}

/// Options of the research entry point `RealDavidson::solve_research`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ResearchOptions {
    pub b1: B1Mode,
    pub sketch: SketchKind,
    /// B5 bordered update of the reduced eigenproblem (`BAANDI_B5=on`).
    pub b5: bool,
    pub stage_timers: bool,
    /// B6 cached projector products for the return certificate (`BAANDI_B6=on`).
    pub b6: bool,
    /// B5 bordered update with the column-parallel coefficient product and
    /// root-parallel secular solve (`BAANDI_B5=par`; implies `b5`).
    pub b5_parallel: bool,
    /// RITZ-FUSED (`BAANDI_RITZ_FUSED=on`): Ritz vectors and residual vectors
    /// in one pass over V and HV; H X is formed only for a contraction restart.
    /// Bit-identical to the default path.
    pub ritz_fused: bool,
    /// NARROW-PROJ (`BAANDI_NARROW_PROJ=on`): blocks of at most
    /// `real_davidson::NARROW_WIDTH` columns are projected against the basis by
    /// streaming dot products and one fused update pass (rounding class).
    pub narrow_proj: bool,
    /// NARROW-REDUCED (`BAANDI_NARROW_REDUCED=on`): the border of the reduced
    /// matrix for at most `NARROW_WIDTH` new columns by the same streaming dot
    /// kernel (rounding class).
    pub narrow_reduced: bool,
}
impl ResearchOptions {
    /// `BAANDI_B1` and `BAANDI_DAVIDSON_STAGES=on`, read once per solve.
    pub fn from_env() -> Self {
        Self {
            b1: B1Mode::from_env(),
            sketch: match dft_numerics::switches::choice(dft_numerics::switches::Id::B1_SKETCH) {
                "dense" => SketchKind::Dense,
                _ => SketchKind::Count,
            },
            // B5 default ON since the adoption rule (dev08 B5-PART-C: 160 runs, zero
            // violations, mechanism fires on 53 percent of strict Davidson eigenproblems,
            // si32 strict T=1 -0.77 percent 10/10 seeds, nothing regresses);
            // `BAANDI_B5=off` disables, `par` threads the update.
            b5: dft_numerics::switches::choice(dft_numerics::switches::Id::B5) != "off",
            b5_parallel: dft_numerics::switches::choice(dft_numerics::switches::Id::B5) == "par",
            ritz_fused: dft_numerics::switches::on(dft_numerics::switches::Id::RITZ_FUSED),
            narrow_proj: dft_numerics::switches::on(dft_numerics::switches::Id::NARROW_PROJ),
            narrow_reduced: dft_numerics::switches::on(dft_numerics::switches::Id::NARROW_REDUCED),
            stage_timers: dft_numerics::switches::on(dft_numerics::switches::Id::DAVIDSON_STAGES),
            b6: dft_numerics::switches::on(dft_numerics::switches::Id::B6),
        }
    }
}

/// Dense-algebra stages timed inside one Davidson iteration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Stage {
    /// Incremental reduced matrix V^T H V.
    Reduced = 0,
    /// Reduced dense eigensolve.
    Eigen = 1,
    /// Ritz vectors and H images (all bands, or the B1 subset).
    Rotate = 2,
    /// Residual vectors and norms.
    Residual = 3,
    /// Projection against the basis and Gram-Schmidt (append and restart).
    Orthonormalize = 4,
    /// Operator applications (corrections, restarts and certificates).
    HApply = 5,
    /// B1 strict sketch: estimates and sketch maintenance.
    Sketch = 6,
    /// Preconditioned corrections.
    Precondition = 7,
}
const STAGES: usize = 8;
/// Correction-width buckets: 0 (returning or no correction), 1, 2-8, 9-64, 65+.
const WIDTHS: usize = 5;

static STAGE_NANOS: [AtomicU64; STAGES] = [const { AtomicU64::new(0) }; STAGES];
static WIDTH_NANOS: [AtomicU64; WIDTHS] = [const { AtomicU64::new(0) }; WIDTHS];
static WIDTH_COUNTS: [AtomicU64; WIDTHS] = [const { AtomicU64::new(0) }; WIDTHS];
static TIMED_SOLVES: AtomicU64 = AtomicU64::new(0);

fn width_bucket(width: usize) -> usize {
    match width {
        0 => 0,
        1 => 1,
        2..=8 => 2,
        9..=64 => 3,
        _ => 4,
    }
}

/// Copyable timer handle; every method is a no-op when off.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Timers {
    on: bool,
}
impl Timers {
    pub(crate) fn new(on: bool) -> Self {
        if on {
            TIMED_SOLVES.fetch_add(1, Relaxed);
        }
        Self { on }
    }
    #[inline]
    pub(crate) fn start(self) -> Option<Instant> {
        self.on.then(Instant::now)
    }
    #[inline]
    pub(crate) fn stop(self, start: Option<Instant>, stage: Stage) {
        if let Some(t) = start {
            STAGE_NANOS[stage as usize].fetch_add(nanos(t), Relaxed);
        }
    }
    /// Whole-iteration wall time, bucketed by the correction width it appended.
    #[inline]
    pub(crate) fn iteration(self, start: Option<Instant>, width: usize) {
        if let Some(t) = start {
            let b = width_bucket(width);
            WIDTH_NANOS[b].fetch_add(nanos(t), Relaxed);
            WIDTH_COUNTS[b].fetch_add(1, Relaxed);
        }
    }
}
fn nanos(t: Instant) -> u64 {
    u64::try_from(t.elapsed().as_nanos()).unwrap_or(u64::MAX)
}

/// Cumulative stage timers of this process, seconds. All zero unless some solve
/// ran with `BAANDI_DAVIDSON_STAGES=on`.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct StageSeconds {
    /// Solves that ran with timers on.
    pub timed_solves: u64,
    /// Iteration wall time per correction-width bucket:
    /// [returning/none, 1, 2-8, 9-64, 65+].
    pub width_seconds: [f64; WIDTHS],
    /// Iterations per bucket (same order).
    pub width_iterations: [u64; WIDTHS],
    pub reduced: f64,
    pub eigen: f64,
    pub rotate: f64,
    pub residual: f64,
    pub orthonormalize: f64,
    pub h_apply: f64,
    pub sketch: f64,
    pub precondition: f64,
}
/// Width bucket labels in the order of `StageSeconds::width_seconds`.
pub const WIDTH_LABELS: [&str; WIDTHS] = ["0", "1", "2-8", "9-64", "65+"];

pub fn stage_seconds() -> StageSeconds {
    let s = |i: usize| STAGE_NANOS[i].load(Relaxed) as f64 * 1e-9;
    let mut out = StageSeconds {
        timed_solves: TIMED_SOLVES.load(Relaxed),
        reduced: s(Stage::Reduced as usize),
        eigen: s(Stage::Eigen as usize),
        rotate: s(Stage::Rotate as usize),
        residual: s(Stage::Residual as usize),
        orthonormalize: s(Stage::Orthonormalize as usize),
        h_apply: s(Stage::HApply as usize),
        sketch: s(Stage::Sketch as usize),
        precondition: s(Stage::Precondition as usize),
        ..StageSeconds::default()
    };
    for b in 0..WIDTHS {
        out.width_seconds[b] = WIDTH_NANOS[b].load(Relaxed) as f64 * 1e-9;
        out.width_iterations[b] = WIDTH_COUNTS[b].load(Relaxed);
    }
    out
}

/// Cumulative B1 work counters of this process (all zero when B1 is off).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct B1Counters {
    /// Solves in which B1 was active.
    pub solves: u64,
    /// Iterations that rotated only the selected (standard) or not clearly
    /// converged (strict) roots.
    pub narrow_iterations: u64,
    /// Iterations of B1-active solves that took the full all-band path
    /// (wide, restart, return, iteration limit, or fallback).
    pub full_iterations: u64,
    /// Ritz columns formed in narrow iterations (each column is one V and one
    /// HV column).
    pub narrow_columns: u64,
    /// Ritz columns the narrow iterations did not form (bands minus formed).
    pub avoided_columns: u64,
    /// Strict: narrow attempts whose explicit residuals all passed, so the
    /// iteration was redone on the full path (return candidates).
    pub fallbacks: u64,
    /// Strict: roots whose sketch estimate was below half their threshold
    /// while the explicit residual of the same iteration exceeded it (checked
    /// in every full iteration of a sketched solve).
    pub sketch_disagreements: u64,
    /// Strict: sketch estimates evaluated (roots x iterations).
    pub sketch_estimates: u64,
}
static B1_COUNTS: [AtomicU64; 8] = [const { AtomicU64::new(0) }; 8];
#[derive(Debug, Clone, Copy)]
pub(crate) enum B1Count {
    Solves = 0,
    Narrow = 1,
    Full = 2,
    NarrowColumns = 3,
    Avoided = 4,
    Fallbacks = 5,
    Disagreements = 6,
    Estimates = 7,
}
#[inline]
pub(crate) fn b1_add(counter: B1Count, value: usize) {
    B1_COUNTS[counter as usize].fetch_add(value as u64, Relaxed);
}
pub fn b1_counters() -> B1Counters {
    let g = |c: B1Count| B1_COUNTS[c as usize].load(Relaxed);
    B1Counters {
        solves: g(B1Count::Solves),
        narrow_iterations: g(B1Count::Narrow),
        full_iterations: g(B1Count::Full),
        narrow_columns: g(B1Count::NarrowColumns),
        avoided_columns: g(B1Count::Avoided),
        fallbacks: g(B1Count::Fallbacks),
        sketch_disagreements: g(B1Count::Disagreements),
        sketch_estimates: g(B1Count::Estimates),
    }
}

/// Cumulative B5 counters (all zero when B5 is off).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct B5Counters {
    /// Reduced eigenproblems solved by the bordered update.
    pub bordered: u64,
    /// One-column iterations where the update declined (conditioning doubt)
    /// and the full solve ran.
    pub declined: u64,
    /// Full solves in B5-active solves (first iteration, wider appends,
    /// restarts, streak refresh, declines).
    pub full: u64,
}
static B5_COUNTS: [AtomicU64; 3] = [const { AtomicU64::new(0) }; 3];
#[derive(Debug, Clone, Copy)]
pub(crate) enum B5Count {
    Bordered = 0,
    Declined = 1,
    Full = 2,
}
#[inline]
pub(crate) fn b5_add(counter: B5Count) {
    B5_COUNTS[counter as usize].fetch_add(1, Relaxed);
}
pub fn b5_counters() -> B5Counters {
    let g = |c: B5Count| B5_COUNTS[c as usize].load(Relaxed);
    B5Counters {
        bordered: g(B5Count::Bordered),
        declined: g(B5Count::Declined),
        full: g(B5Count::Full),
    }
}

/// Cumulative algebra-stream counters (all zero when the switches are off).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct AlgebraCounters {
    /// Iterations whose Ritz vectors and residuals came from the fused pass.
    pub fused_iterations: u64,
    /// Contraction restarts that formed H X late (fused iterations only).
    pub deferred_hx: u64,
    /// Projection passes taken by the narrow kernel.
    pub narrow_projections: u64,
    /// Reduced-matrix borders taken by the narrow kernel.
    pub narrow_reduced: u64,
    /// Bordered updates whose coefficient product ran column-parallel.
    pub b5_parallel: u64,
}
static ALGEBRA_COUNTS: [AtomicU64; 5] = [const { AtomicU64::new(0) }; 5];
#[derive(Debug, Clone, Copy)]
pub(crate) enum AlgebraCount {
    Fused = 0,
    DeferredHx = 1,
    NarrowProjection = 2,
    NarrowReduced = 3,
    B5Parallel = 4,
}
#[inline]
pub(crate) fn algebra_add(counter: AlgebraCount) {
    ALGEBRA_COUNTS[counter as usize].fetch_add(1, Relaxed);
}
pub fn algebra_counters() -> AlgebraCounters {
    let g = |c: AlgebraCount| ALGEBRA_COUNTS[c as usize].load(Relaxed);
    AlgebraCounters {
        fused_iterations: g(AlgebraCount::Fused),
        deferred_hx: g(AlgebraCount::DeferredHx),
        narrow_projections: g(AlgebraCount::NarrowProjection),
        narrow_reduced: g(AlgebraCount::NarrowReduced),
        b5_parallel: g(AlgebraCount::B5Parallel),
    }
}
