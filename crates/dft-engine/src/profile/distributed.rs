//! Shared distributed-stage instrumentation and counters.
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::time::Instant;

#[derive(Clone, Copy)]
#[cfg_attr(not(feature = "strict-oracle"), allow(dead_code))]
pub(crate) enum Stage {
    Density,
    Xc,
    NonlocalEnergy,
    ScreenIntegrals,
    ScreenClusters,
    AtomicStart,
    /// Replicated solve interface before the distributed solver (operator
    /// update, kinetic diagonal, ORC packing of the N x B start; MPI-SCALING-GAP).
    SolvePrepare,
    /// Replicated solve interface after it (ORC unpacking of the N x B result).
    SolveUnpack,
    /// MPI-LAZY-ORBITALS: on-demand gather and unpacking of a lazy block.
    LazyMaterialize,
    /// MPI-LAZY-ORBITALS: kinetic energy of the audit from the local rows.
    LazyKinetic,
    /// BREAKTHROUGH-4: MPI-KEPT-START update of the kept H X to the new
    /// potential (the delta-V application of `prepare_m2_kept_start`).
    KeptUpdate,
    /// BREAKTHROUGH-4: the solve interface between the preparation and the
    /// solver call (local diagonal, start rows, kept-start bookkeeping).
    SolveEntry,
    /// BREAKTHROUGH-4: the solve interface after the solver call (lazy handle
    /// copy of the local rows, kept-start flags, the old potential copy).
    SolveFinish,
    /// BREAKTHROUGH-4 (`BAANDI_PHASE_WAIT=on` only): timed barrier at the entry
    /// of every M2 solve (arrival skew of the ranks after the potential phase).
    EntryWait,
    /// BREAKTHROUGH-4: `ensure_state` of a real Gamma solve (the replicated
    /// serial real Hamiltonian: reused with a potential copy, or rebuilt from
    /// the complex operator after a new problem or a ladder stage).
    RealState,
}
#[cfg(feature = "strict-oracle")]
const STAGE_NAMES: [&str; 15] = [
    "density",
    "xc",
    "nonlocal_energy",
    "screen_integrals",
    "screen_clusters",
    "atomic_start",
    "solve_prepare",
    "solve_unpack",
    "lazy_materialize",
    "lazy_kinetic",
    "kept_update",
    "solve_entry",
    "solve_finish",
    "entry_wait",
    "real_state",
];
#[cfg(feature = "strict-oracle")]
const STAGE_CALL_NAMES: [&str; 15] = [
    "calls.density",
    "calls.xc",
    "calls.nonlocal_energy",
    "calls.screen_integrals",
    "calls.screen_clusters",
    "calls.atomic_start",
    "calls.solve_prepare",
    "calls.solve_unpack",
    "calls.lazy_materialize",
    "calls.lazy_kinetic",
    "calls.kept_update",
    "calls.solve_entry",
    "calls.solve_finish",
    "calls.entry_wait",
    "calls.real_state",
];
static STAGE_NANOS: [AtomicU64; 15] = [const { AtomicU64::new(0) }; 15];
static STAGE_CALLS: [AtomicUsize; 15] = [const { AtomicUsize::new(0) }; 15];
/// F21: `start` from `dft_numerics::stage_clock::start` (or `keep`); `None`
/// with `BAANDI_STAGE_TIMERS=off` counts the entry without reading a clock.
pub(crate) fn add(stage: Stage, start: Option<Instant>) {
    STAGE_CALLS[stage as usize].fetch_add(1, Ordering::Relaxed);
    let ns = dft_numerics::stage_clock::elapsed_ns(start);
    if start.is_some() {
        STAGE_NANOS[stage as usize].fetch_add(ns, Ordering::Relaxed);
    }
    // Per-thread mirror for the phase profile (`crate::profile`).
    dft_numerics::stage_clock::local::add(
        dft_numerics::stage_clock::local::M2_STAGE + stage as usize,
        ns,
    );
}
#[derive(Clone, Copy)]
#[cfg_attr(not(feature = "strict-oracle"), allow(dead_code))]
pub(crate) enum Count {
    Density,
    Xc,
    NonlocalEnergy,
    ScreenIntegrals,
    ScreenClusterBands,
    ScreenDistClusters,
    AtomicSplitStarts,
    /// MPI-KEPT-START (`BAANDI_M2_REUSE=on`): solves armed with the kept X and
    /// H X (the solver may still decline).
    ReuseArmed,
    /// Armed solves the solver started from the kept X and H X.
    ReuseTaken,
    /// Armed solves the solver declined (initial differs, Gram defect).
    ReuseDeclined,
    /// Armed solves whose potential did not change (no transform at all).
    ReuseUnchanged,
    /// Eligible solves without a kept return (first solve, rebuilt operator,
    /// failed or non-keeping previous solve).
    ReuseNotKept,
    /// H columns not applied thanks to a taken kept start.
    ReuseColumnsSaved,
    /// Density builds that kept plane slabs for the next start.
    ReusePlaneKeeps,
    /// States whose plane slabs exceeded the memory gate.
    ReuseGated,
    /// Largest kept plane-slab storage of a rank, bytes (maximum, not a sum).
    ReusePlaneBytes,
    /// MPI-LAZY-ORBITALS (`BAANDI_M2_LAZY=on`): solves that returned only the
    /// local rows (no N x B allgatherv, no complex unpack).
    LazyReturns,
    /// Solves started from a lazy block's local rows (no ORC pack, no row
    /// extraction).
    LazyStarts,
    /// Lazy blocks gathered on demand (result, checkpoint, ladder, audits).
    LazyMaterializations,
    /// Lazy starts whose rows did not match the operator on some rank (then
    /// gathered and packed; expected 0).
    LazyStartMismatches,
    /// `BAANDI_M2_ATOMIC_BLOCK=on`: atomic starts orthonormalized by the
    /// row-distributed CholQR2.
    AtomicBlockStarts,
    /// Block attempts that fell back to the split Gram-Schmidt (a rank could
    /// not fill its rows, a pivot below the dependence guard, nonfinite data).
    AtomicBlockFallbacks,
    /// INIT-SCALE-2: constructions whose ionic local coefficients were split
    /// over the ranks (one allgatherv; bit-identical).
    IonicLocalSplits,
    /// INIT-SCALE-3: Ewald evaluations split over the ranks (chained sums;
    /// bit-identical; only from `EWALD_SPLIT_MIN_RANKS` ranks on).
    EwaldSplits,
    /// INIT-DENSITY-SPLIT: atomic-superposition start densities whose
    /// coefficients were split over the ranks (one allgatherv; bit-identical).
    DensitySplits,
    /// LAZY-RESULT: delivered results that kept their lazy blocks (the caller
    /// reads no orbitals; no gather at delivery).
    LazyResults,
}
#[cfg(feature = "strict-oracle")]
const COUNT_NAMES: [&str; 26] = [
    "density_builds",
    "xc_evaluations",
    "nonlocal_energy_audits",
    "screen_integrals",
    "screen_cluster_bands",
    "screen_dist_clusters",
    "atomic_split_starts",
    "reuse_armed",
    "reuse_taken",
    "reuse_declined",
    "reuse_unchanged_potential",
    "reuse_not_kept",
    "reuse_columns_saved",
    "reuse_plane_keeps",
    "reuse_gated",
    "reuse_plane_bytes_max",
    "lazy_returns",
    "lazy_starts",
    "lazy_materializations",
    "lazy_start_mismatches",
    "atomic_block_starts",
    "atomic_block_fallbacks",
    "ionic_local_splits",
    "ewald_splits",
    "density_splits",
    "lazy_results",
];
static COUNTS: [AtomicUsize; 26] = [const { AtomicUsize::new(0) }; 26];
#[cfg(feature = "strict-oracle")]
pub(crate) fn count_atomic_block_fallback() {
    count(Count::AtomicBlockFallbacks);
}
#[cfg(feature = "strict-oracle")]
pub(crate) fn count_max(c: Count, n: usize) {
    COUNTS[c as usize].fetch_max(n, Ordering::Relaxed);
}
#[cfg(feature = "strict-oracle")]
pub(crate) fn count_n(c: Count, n: usize) {
    COUNTS[c as usize].fetch_add(n, Ordering::Relaxed);
}
pub(crate) fn count(c: Count) {
    COUNTS[c as usize].fetch_add(1, Ordering::Relaxed);
}

/// Engine-level M2d counters and stage seconds of this process (counters end
/// with one `calls.<stage>` entry per stage; seconds 0 with
/// `BAANDI_STAGE_TIMERS=off`).
#[cfg(feature = "strict-oracle")]
pub fn stats() -> dft_physics::distributed_real_gamma::M2Counters {
    let handoff_stages = [(
        "handoff.orbital_copy",
        crate::refscf::prolong::orbital_copy_seconds(),
    )];
    (
        COUNT_NAMES
            .iter()
            .zip(&COUNTS)
            .map(|(n, c)| (*n, c.load(Ordering::Relaxed)))
            .chain(
                STAGE_CALL_NAMES
                    .iter()
                    .zip(&STAGE_CALLS)
                    .map(|(n, c)| (*n, c.load(Ordering::Relaxed))),
            )
            .chain(dft_physics::hamiltonian::owned_state_stats())
            .collect(),
        STAGE_NAMES
            .iter()
            .zip(&STAGE_NANOS)
            .map(|(n, c)| (*n, c.load(Ordering::Relaxed) as f64 * 1e-9))
            .chain(handoff_stages)
            .collect(),
    )
}
