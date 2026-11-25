//! Shared distributed-stage instrumentation and counters.
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::time::Instant;

#[derive(Clone, Copy)]
#[cfg_attr(not(feature = "strict-oracle"), allow(dead_code))]
pub(crate) enum Stage {
    Density = 0,
    Xc = 1,
    NonlocalEnergy = 2,
    ScreenIntegrals = 3,
    ScreenClusters = 4,
    AtomicStart = 5,
    /// Replicated solve interface before the distributed solver (operator
    /// update, kinetic diagonal, ORC packing of the N x B start; MPI-SCALING-GAP).
    SolvePrepare = 6,
    /// Replicated solve interface after it (ORC unpacking of the N x B result).
    SolveUnpack = 7,
    /// MPI-LAZY-ORBITALS: on-demand gather and unpacking of a lazy block.
    LazyMaterialize = 8,
    /// BREAKTHROUGH-4 (`BAANDI_PHASE_WAIT=on` only): timed barrier at the entry
    /// of every M2 solve (arrival skew of the ranks after the potential phase).
    EntryWait = 13,
    /// BREAKTHROUGH-4: `ensure_state` of a real Gamma solve (the replicated
    /// serial real Hamiltonian: reused with a potential copy, or rebuilt from
    /// the complex operator after a new problem or a ladder stage).
    RealState = 14,
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
    Density = 0,
    Xc = 1,
    NonlocalEnergy = 2,
    ScreenIntegrals = 3,
    ScreenClusterBands = 4,
    ScreenDistClusters = 5,
    AtomicSplitStarts = 6,
    /// Lazy blocks gathered on demand (result, checkpoint, ladder, audits).
    LazyMaterializations = 18,
    /// `BAANDI_M2_ATOMIC_BLOCK=on`: atomic starts orthonormalized by the
    /// row-distributed CholQR2.
    AtomicBlockStarts = 20,
    /// Block attempts that fell back to the split Gram-Schmidt (a rank could
    /// not fill its rows, a pivot below the dependence guard, nonfinite data).
    AtomicBlockFallbacks = 21,
    /// INIT-SCALE-2: constructions whose ionic local coefficients were split
    /// over the ranks (one allgatherv; bit-identical).
    IonicLocalSplits = 22,
    /// INIT-SCALE-3: Ewald evaluations split over the ranks (chained sums;
    /// bit-identical; only from `EWALD_SPLIT_MIN_RANKS` ranks on).
    EwaldSplits = 23,
    /// INIT-DENSITY-SPLIT: atomic-superposition start densities whose
    /// coefficients were split over the ranks (one allgatherv; bit-identical).
    DensitySplits = 24,
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
