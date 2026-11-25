//! Frozen legacy Davidson arithmetic and diagnostic controls.
//! Strict solves use B1 CountSketch and serial B5 bordered updates; Ritz-change
//! solves use lean algebra. Fused/narrow alternatives and restart-history/cache
//! experiments are retired. Stage timers and process counters remain diagnostic.
use std::sync::atomic::{AtomicU64, Ordering::Relaxed};
use std::time::Instant;

/// Frozen B1 strict sketch margin (R11b).
pub const B1_DEFAULT_MARGIN: f64 = 0.5;
/// Miss bound of the automatic B1 sketch row rule, in standard deviations of
/// the relative sketch error. A sketched norm has relative standard error about
/// 1/sqrt(2 s) for s rows; a root whose true residual exceeds its threshold is
/// skipped only when its estimate falls below `margin` times the threshold, a
/// relative error of at least (1 - margin), i.e. (1 - margin) sqrt(2 s) standard
/// deviations. The default (margin 0.5, s = 128) gives 8; the rule keeps 8 at
/// any margin. A statistical rule, independent of the material.
const B1_SKETCH_SIGMAS: f64 = 8.0;

/// Diagnostic controls of the frozen legacy solver. Arithmetic is fixed.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ResearchOptions {
    pub stage_timers: bool,
}
impl ResearchOptions {
    /// Fixed eight-sigma rule at margin 0.5, capped by the caller.
    pub fn b1_sketch_rows(&self) -> usize {
        let z = B1_SKETCH_SIGMAS / (1. - B1_DEFAULT_MARGIN);
        (0.5 * z * z).ceil() as usize
    }
    pub fn from_env() -> Self {
        Self {
            stage_timers: dft_numerics::switches::on(dft_numerics::switches::Id::DAVIDSON_STAGES),
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

// Per run stage (`dft_numerics::work_timers::RunStage`: fine, ladder coarse,
// verification), selected by the process-wide index the engine sets between solves.
const RUN_STAGES: usize = dft_numerics::work_timers::RUN_STAGES;
static STAGE_NANOS: [[AtomicU64; STAGES]; RUN_STAGES] =
    [const { [const { AtomicU64::new(0) }; STAGES] }; RUN_STAGES];
static WIDTH_NANOS: [[AtomicU64; WIDTHS]; RUN_STAGES] =
    [const { [const { AtomicU64::new(0) }; WIDTHS] }; RUN_STAGES];
static WIDTH_COUNTS: [[AtomicU64; WIDTHS]; RUN_STAGES] =
    [const { [const { AtomicU64::new(0) }; WIDTHS] }; RUN_STAGES];
static TIMED_SOLVES: [AtomicU64; RUN_STAGES] = [const { AtomicU64::new(0) }; RUN_STAGES];
fn run_stage() -> usize {
    dft_numerics::work_timers::run_stage()
}

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
            TIMED_SOLVES[run_stage()].fetch_add(1, Relaxed);
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
            let ns = nanos(t);
            STAGE_NANOS[run_stage()][stage as usize].fetch_add(ns, Relaxed);
            // Per-thread mirror for the engine's phase profile.
            dft_numerics::stage_clock::local::add(
                dft_numerics::stage_clock::local::SERIAL_STAGE + stage as usize,
                ns,
            );
        }
    }
    /// Whole-iteration wall time, bucketed by the correction width it appended.
    #[inline]
    pub(crate) fn iteration(self, start: Option<Instant>, width: usize) {
        if let Some(t) = start {
            let b = width_bucket(width);
            let r = run_stage();
            WIDTH_NANOS[r][b].fetch_add(nanos(t), Relaxed);
            WIDTH_COUNTS[r][b].fetch_add(1, Relaxed);
        }
    }
}
fn nanos(t: Instant) -> u64 {
    u64::try_from(t.elapsed().as_nanos()).unwrap_or(u64::MAX)
}

/// Distributed solver (`distributed_davidson::STAGES` index `dist`, `ns` of one
/// stage entry, 0 count-only): charge the matching serial stage of the current
/// run stage, so that `stage_seconds_by_run_stage` covers the MPI path too.
/// Mapping: reduced_gram -> reduced, reduced_eigen -> eigen, rotate, residual,
/// orthonormalize, apply_total -> h_apply, precondition; solver_total counts
/// one timed solve. The other distributed labels have no serial counterpart
/// (they stay in the phase profile). Called only with the work timers on
/// (`BAANDI_DAVIDSON_STAGES=on`); not mirrored to the per-thread profile slots
/// (the profile reads the distributed slots already).
pub(crate) fn add_distributed(dist: usize, ns: u64) {
    let stage = match dist {
        0 => Stage::Reduced,
        1 => Stage::Eigen,
        2 => Stage::Rotate,
        3 => Stage::Residual,
        4 => Stage::Orthonormalize,
        5 => Stage::HApply,
        12 => Stage::Precondition,
        15 => {
            TIMED_SOLVES[run_stage()].fetch_add(1, Relaxed);
            return;
        }
        _ => return,
    };
    STAGE_NANOS[run_stage()][stage as usize].fetch_add(ns, Relaxed);
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

/// Cumulative stage timers of this process, summed over the run stages.
pub fn stage_seconds() -> StageSeconds {
    let all = stage_seconds_by_run_stage();
    let mut out = StageSeconds::default();
    for r in &all {
        out.timed_solves += r.timed_solves;
        for b in 0..WIDTHS {
            out.width_seconds[b] += r.width_seconds[b];
            out.width_iterations[b] += r.width_iterations[b];
        }
        out.reduced += r.reduced;
        out.eigen += r.eigen;
        out.rotate += r.rotate;
        out.residual += r.residual;
        out.orthonormalize += r.orthonormalize;
        out.h_apply += r.h_apply;
        out.sketch += r.sketch;
        out.precondition += r.precondition;
    }
    out
}

/// Cumulative stage timers per run stage, in the index order of
/// `dft_numerics::work_timers::RunStage` (fine, ladder coarse, verification).
pub fn stage_seconds_by_run_stage() -> [StageSeconds; RUN_STAGES] {
    std::array::from_fn(stage_seconds_of)
}
fn stage_seconds_of(r: usize) -> StageSeconds {
    let s = |i: usize| STAGE_NANOS[r][i].load(Relaxed) as f64 * 1e-9;
    let mut out = StageSeconds {
        timed_solves: TIMED_SOLVES[r].load(Relaxed),
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
        out.width_seconds[b] = WIDTH_NANOS[r][b].load(Relaxed) as f64 * 1e-9;
        out.width_iterations[b] = WIDTH_COUNTS[r][b].load(Relaxed);
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

/// Historical algebra-stream counters (all zero after retirement).
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
/// Retired alternatives retain zero-valued diagnostic schema fields.
pub fn algebra_counters() -> AlgebraCounters {
    AlgebraCounters::default()
}

/// GEMM sites of the serial real Davidson for the work counters (flops are
/// 2 m n k per product, counted from the shapes at the call; no effect on the
/// arithmetic). The in-block Gram-Schmidt loops are counted as 4 n flops per
/// projected (column, earlier column) pair and repetition.
#[derive(Debug, Clone, Copy)]
pub(crate) enum FlopSite {
    /// Border of the reduced matrix V^T (H V).
    Reduced = 0,
    /// Lean: border of the basis overlap V^T V.
    Overlap = 1,
    /// Ritz vectors, H images and residuals inside the loop (classic: every
    /// iteration incl. the returning one; lean: the unconverged roots only).
    RitzLoop = 2,
    /// Lean: all-band Ritz vectors and residuals at the return.
    Return = 3,
    /// Contraction restarts (lean: X and H X of all bands; classic: the late
    /// H X of fused iterations).
    Contraction = 4,
    /// Projection of a block against the basis (Q^T X and Q S, every pass).
    Projection = 5,
    /// In-block Gram-Schmidt (panel GEMMs and the scalar sweep).
    GramSchmidt = 6,
    /// Lean fallback: rotation of the whole basis to its Ritz basis.
    Fallback = 7,
    /// Lean small dense products of order m (Cholesky border, inverse factor,
    /// reduced border of T H T^T, coefficients T^T Y).
    Small = 8,
}
const FLOP_SITES: usize = 9;
static FLOPS: [AtomicU64; FLOP_SITES] = [const { AtomicU64::new(0) }; FLOP_SITES];
/// Full dense reduced eigensolves: count and sum of m^3.
static EIGEN_FULL: [AtomicU64; 2] = [const { AtomicU64::new(0) }; 2];
#[inline]
pub(crate) fn flops_add(site: FlopSite, rows: usize, inner: usize, cols: usize) {
    let f = 2_u64
        .saturating_mul(rows as u64)
        .saturating_mul(inner as u64)
        .saturating_mul(cols as u64);
    FLOPS[site as usize].fetch_add(f, Relaxed);
}
#[inline]
pub(crate) fn flops_add_raw(site: FlopSite, flops: u64) {
    FLOPS[site as usize].fetch_add(flops, Relaxed);
}
#[inline]
pub(crate) fn eigen_full_add(order: usize) {
    EIGEN_FULL[0].fetch_add(1, Relaxed);
    EIGEN_FULL[1].fetch_add((order as u64).saturating_pow(3), Relaxed);
}
/// Cumulative GEMM flops of the serial real Davidson by site (both paths).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct FlopCounters {
    pub reduced: u64,
    pub overlap: u64,
    pub ritz_loop: u64,
    pub ret: u64,
    pub contraction: u64,
    pub projection: u64,
    pub gram_schmidt: u64,
    pub fallback: u64,
    pub small: u64,
    /// Full dense reduced eigensolves (B5 bordered updates excluded).
    pub eigen_full: u64,
    /// Sum of m^3 over those solves.
    pub eigen_full_m3: u64,
}
impl FlopCounters {
    /// Sum of the GEMM sites (the eigensolves are not GEMM work).
    pub fn total(&self) -> u64 {
        self.reduced
            + self.overlap
            + self.ritz_loop
            + self.ret
            + self.contraction
            + self.projection
            + self.gram_schmidt
            + self.fallback
            + self.small
    }
}
pub fn flop_counters() -> FlopCounters {
    let g = |s: FlopSite| FLOPS[s as usize].load(Relaxed);
    FlopCounters {
        reduced: g(FlopSite::Reduced),
        overlap: g(FlopSite::Overlap),
        ritz_loop: g(FlopSite::RitzLoop),
        ret: g(FlopSite::Return),
        contraction: g(FlopSite::Contraction),
        projection: g(FlopSite::Projection),
        gram_schmidt: g(FlopSite::GramSchmidt),
        fallback: g(FlopSite::Fallback),
        small: g(FlopSite::Small),
        eigen_full: EIGEN_FULL[0].load(Relaxed),
        eigen_full_m3: EIGEN_FULL[1].load(Relaxed),
    }
}

/// Cumulative counters of the lean Davidson (all zero when it never ran).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct LeanCounters {
    pub solves: u64,
    pub iterations: u64,
    /// Ritz and residual columns formed inside the loop (unconverged roots).
    pub residual_columns: u64,
    /// Correction blocks appended (each one is a pivot test).
    pub appends: u64,
    pub appended_columns: u64,
    /// Appends whose overlap pivot ratio fell below the bound (fallback stage
    /// 1: the block alone was projected against the basis through the factor,
    /// passed through the in-block Gram-Schmidt and tested again).
    pub reprojections: u64,
    /// Appends that failed stage 1 too (stage 2: the basis was rotated to its
    /// Ritz basis once and the block projected classically).
    pub fallbacks: u64,
    /// Correction columns dropped by the classic rank decisions of either stage.
    pub fallback_dropped_columns: u64,
    /// Bases rebuilt by projection and H application (factor failure right
    /// after a start or a contraction; expected zero).
    pub rebuilds: u64,
    pub rebuild_columns: u64,
    pub contractions: u64,
    pub returns: u64,
    /// Appends whose smallest pivot ratio was below 1e-1 and 1e-2 (headroom
    /// above the bound).
    pub pivot_below_1e_1: u64,
    pub pivot_below_1e_2: u64,
    /// Smallest pivot ratio seen in any append (1 when none).
    pub pivot_min: f64,
    /// REUSE-PSI: lean solves that started from the kept X and H X (no
    /// initial H application), armed starts declined because `initial` was
    /// not the kept X, returns that kept X and H X, and the H columns the
    /// kept starts did not apply.
    pub kept_starts: u64,
    pub kept_declined: u64,
    pub kept_stored: u64,
    pub kept_columns: u64,
    /// Retired projector-cache diagnostics: returns that left P^T X in the
    /// workspace, returns whose projection cache was incomplete (the caller
    /// computes P^T X directly), kept starts that also took the kept P^T X,
    /// and projector-product columns copied from the operator after applies.
    pub proj_cache_returns: u64,
    pub proj_cache_misses: u64,
    pub proj_cache_kept_loads: u64,
    pub proj_cache_columns: u64,
    /// Retired fused-border counters (always zero), retained for CLI reporting.
    pub fused_borders: u64,
    pub fused_wasted_columns: u64,
    /// Retired restart-history diagnostics: contractions that kept at least one
    /// previous direction, directions kept, candidates dropped by the pivot
    /// rule in coefficient space, and contractions whose new basis failed the
    /// factor test and fell back to X alone.
    pub keep_contractions: u64,
    pub kept_directions: u64,
    pub keep_dropped: u64,
    pub keep_fallbacks: u64,
}
impl LeanCounters {
    /// Stage-2 fallbacks (basis orthonormalizations) per append (the
    /// work-order kill rule: below 0.05).
    pub fn fallback_rate(&self) -> f64 {
        if self.appends == 0 {
            0.
        } else {
            self.fallbacks as f64 / self.appends as f64
        }
    }
    /// Appends that needed any fallback stage (pivot test failed once).
    pub fn reprojection_rate(&self) -> f64 {
        if self.appends == 0 {
            0.
        } else {
            self.reprojections as f64 / self.appends as f64
        }
    }
}
#[derive(Debug, Clone, Copy)]
pub(crate) enum LeanCount {
    Solves = 0,
    Iterations = 1,
    ResidualColumns = 2,
    Appends = 3,
    AppendedColumns = 4,
    Fallbacks = 5,
    FallbackDropped = 6,
    Rebuilds = 7,
    RebuildColumns = 8,
    Contractions = 9,
    Returns = 10,
    PivotBelow1 = 11,
    PivotBelow2 = 12,
    Reprojections = 13,
    KeptStarts = 14,
    KeptDeclined = 15,
    KeptStored = 16,
    KeptColumns = 17,
    ProjCacheReturns = 18,
    ProjCacheMisses = 19,
    ProjCacheKeptLoads = 20,
    ProjCacheColumns = 21,
    FusedBorders = 22,
    FusedWasted = 23,
    KeepContractions = 24,
    KeptDirections = 25,
    KeepDropped = 26,
    KeepFallbacks = 27,
}
static LEAN_COUNTS: [AtomicU64; 28] = [const { AtomicU64::new(0) }; 28];
/// Smallest pivot ratio, as f64 bits (positive floats order like their bits).
static LEAN_PIVOT_MIN: AtomicU64 = AtomicU64::new(0x3ff0_0000_0000_0000);
#[inline]
pub(crate) fn lean_add(counter: LeanCount, value: usize) {
    LEAN_COUNTS[counter as usize].fetch_add(value as u64, Relaxed);
}
pub(crate) fn lean_pivot(ratio: f64) {
    if ratio < 1e-1 {
        lean_add(LeanCount::PivotBelow1, 1);
    }
    if ratio < 1e-2 {
        lean_add(LeanCount::PivotBelow2, 1);
    }
    let bits = if ratio.is_finite() && ratio > 0. {
        ratio.to_bits()
    } else {
        0
    };
    LEAN_PIVOT_MIN.fetch_min(bits, Relaxed);
}
pub fn lean_counters() -> LeanCounters {
    let g = |c: LeanCount| LEAN_COUNTS[c as usize].load(Relaxed);
    LeanCounters {
        solves: g(LeanCount::Solves),
        iterations: g(LeanCount::Iterations),
        residual_columns: g(LeanCount::ResidualColumns),
        appends: g(LeanCount::Appends),
        appended_columns: g(LeanCount::AppendedColumns),
        reprojections: g(LeanCount::Reprojections),
        fallbacks: g(LeanCount::Fallbacks),
        fallback_dropped_columns: g(LeanCount::FallbackDropped),
        rebuilds: g(LeanCount::Rebuilds),
        rebuild_columns: g(LeanCount::RebuildColumns),
        contractions: g(LeanCount::Contractions),
        returns: g(LeanCount::Returns),
        pivot_below_1e_1: g(LeanCount::PivotBelow1),
        pivot_below_1e_2: g(LeanCount::PivotBelow2),
        pivot_min: f64::from_bits(LEAN_PIVOT_MIN.load(Relaxed)),
        kept_starts: g(LeanCount::KeptStarts),
        kept_declined: g(LeanCount::KeptDeclined),
        kept_stored: g(LeanCount::KeptStored),
        kept_columns: g(LeanCount::KeptColumns),
        proj_cache_returns: g(LeanCount::ProjCacheReturns),
        proj_cache_misses: g(LeanCount::ProjCacheMisses),
        proj_cache_kept_loads: g(LeanCount::ProjCacheKeptLoads),
        proj_cache_columns: g(LeanCount::ProjCacheColumns),
        fused_borders: g(LeanCount::FusedBorders),
        fused_wasted_columns: g(LeanCount::FusedWasted),
        keep_contractions: g(LeanCount::KeepContractions),
        kept_directions: g(LeanCount::KeptDirections),
        keep_dropped: g(LeanCount::KeepDropped),
        keep_fallbacks: g(LeanCount::KeepFallbacks),
    }
}
