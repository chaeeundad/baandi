//! MPI stage M2d: distributed SCF stages on the M2 path (`BAANDI_MPI_M2=on` with a
//! communicator of two or more ranks; the default and the one-rank paths never
//! enter this module).
//!
//! * Density: occupation-weighted |psi_b|^2 of the distributed ORC eigenvectors on
//!   each rank's FFT planes (stick FFT, one batched transpose), summed over bands
//!   in ascending order on the planes and assembled into the replicated grid with
//!   one allgatherv (a pure permutation: the full density does not depend on how
//!   the planes are split). The per-band plane contributions are kept for the
//!   certificate screen (`weighted_integrals`: local sums plus one tree_sum).
//! * Potential: LDA (PZ, PW) exchange-correlation is pointwise; each rank evaluates
//!   a contiguous block of the row-major grid, the potential is assembled with one
//!   allgatherv (bit-identical values to the replicated evaluation) and the energy
//!   and potential integral are reduced with one deterministic tree_sum. Hartree
//!   (one 3D FFT pair) stays replicated.
//! * Energies: the nonlocal energy of the audit from the distributed projections
//!   (one tree_sum); the output-density XC energy as above.
//! * Mixing stays replicated (0.3 percent of the si32 wall at one rank, ledger
//!   MPI-M1-GATE4; distributing it would add a grid allgatherv per step that costs
//!   about what it saves).
//!
//! Determinism (ADR 0041): every reduction is a fixed-order local sum followed by
//! the binomial tree_sum, so a fixed rank count is bit-reproducible; different
//! rank counts differ by rounding only. Every collective here is entered by all
//! ranks under the same replicated conditions (no data-dependent branch).
use crate::parallel::SharedComm;
use crate::{Functional, Pw92Prefactor, ScfError, ScfOptions};
use dft_mpi::{Collective, block_counts, block_displacements};
use dft_physics::data::{Array3, Hartree};
use dft_physics::{PhysicsError, xc};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::time::Instant;

/// The M2 communicator when the M2 path is requested and two or more ranks run.
pub(crate) fn active() -> Option<SharedComm> {
    SharedComm::active().filter(|_| crate::parallel::m2_enabled())
}

/// Engine-side M2 switches of the coarse stage and the certificate
/// (docs/performance/mpi-coarse-certificate-analysis.md; both default off).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct EngineTuning {
    /// `BAANDI_M2_SCREEN=dist`: certificate-screen clusters from the distributed
    /// band values (plane-local sums, one tree_sum per cluster).
    pub screen_dist: bool,
    /// `BAANDI_M2_ATOMIC=split`: atomic-start span, Gram-Schmidt first passes,
    /// H Q and Gram columns split over ranks.
    pub atomic_split: bool,
    /// `BAANDI_M2_ATOMIC_BLOCK=on` (default off, rounding class): with
    /// `atomic_split`, the span is filled by row blocks and orthonormalized by
    /// a row-distributed CholQR2 instead of the column-cyclic two-pass
    /// Gram-Schmidt (whose second pass is a serial chain over the columns).
    pub atomic_block: bool,
}
impl EngineTuning {
    pub fn from_env() -> Self {
        use dft_numerics::switches::{Id, choice};
        Self {
            screen_dist: choice(Id::M2_SCREEN) == "dist",
            atomic_split: choice(Id::M2_ATOMIC) == "split",
            atomic_block: choice(Id::M2_ATOMIC_BLOCK) == "on",
        }
    }
}
thread_local! {
    static TUNING_OVERRIDE: std::cell::Cell<Option<EngineTuning>> = const { std::cell::Cell::new(None) };
}
/// Override the engine M2 switches for the calling thread (simulated ranks in
/// tests); `None` restores the environment.
pub fn set_tuning_override(t: Option<EngineTuning>) {
    TUNING_OVERRIDE.with(|o| o.set(t));
}
/// The engine M2 switches in force on this thread.
pub fn tuning() -> EngineTuning {
    TUNING_OVERRIDE.with(|o| o.get()).unwrap_or_else(|| {
        static ENV: std::sync::OnceLock<EngineTuning> = std::sync::OnceLock::new();
        *ENV.get_or_init(EngineTuning::from_env)
    })
}

thread_local! {
    // Diagnostic: the last screen estimates of this thread (tests compare the
    // replicated and distributed screens); bounded, no numerical use.
    static SCREEN_BOUNDS: std::cell::RefCell<Vec<f64>> = const { std::cell::RefCell::new(Vec::new()) };
}
pub(crate) fn record_screen_bound(bound: f64) {
    SCREEN_BOUNDS.with(|b| {
        let mut b = b.borrow_mut();
        if b.len() < 64 {
            b.push(bound);
        }
    });
}
/// Certificate-screen estimates (Ha) recorded on this thread since the last
/// call, oldest first (at most 64; diagnostics and tests).
pub fn take_screen_bounds() -> Vec<f64> {
    SCREEN_BOUNDS.with(|b| std::mem::take(&mut *b.borrow_mut()))
}

/// `BAANDI_M2_BATCH=off` restores the per-pair FFT transposes of M2b/M2c (A/B).
pub(crate) fn batched_transposes() -> bool {
    dft_numerics::switches::on(dft_numerics::switches::Id::M2_BATCH)
}

/// Per-band plane-local density contributions of the latest M2 density build.
pub(crate) struct M2Density {
    /// bands x plane points, [band][plane][row][z].
    pub(crate) contributions: Vec<f64>,
    pub(crate) bands: usize,
    pub(crate) points: usize,
    /// Row-major grid offset of every plane-local point.
    pub(crate) offsets: Vec<usize>,
    /// True when the latest SCF map built its density here (else the serial
    /// `BandDensity` holds the current contributions).
    pub(crate) valid: bool,
    pub(crate) partial: Vec<f64>,
}

impl M2Density {
    pub(crate) fn new(offsets: Vec<usize>) -> Self {
        Self {
            contributions: Vec::new(),
            bands: 0,
            points: offsets.len(),
            offsets,
            valid: false,
            partial: Vec::new(),
        }
    }

    /// out[b] = sum over grid points of band b's contribution times `field`
    /// (full row-major grid), as `BandDensity::weighted_integrals`; local sums in
    /// plane order plus one tree_sum. Collective.
    pub(crate) fn weighted_integrals<K: Collective>(
        &mut self,
        comm: &K,
        field: &[f64],
        out: &mut [f64],
    ) -> Result<(), ScfError> {
        let start = dft_numerics::stage_clock::start();
        self.partial.clear();
        self.partial.resize(out.len(), 0.);
        for (b, p) in self.partial.iter_mut().enumerate().take(self.bands) {
            let c = &self.contributions[b * self.points..(b + 1) * self.points];
            *p = c
                .iter()
                .zip(&self.offsets)
                .map(|(c, &o)| c * field[o])
                .sum();
        }
        comm.tree_sum_f64(&self.partial, out)
            .map_err(|e| ScfError::Numerical(e.to_string()))?;
        add(Stage::ScreenIntegrals, start);
        count(Count::ScreenIntegrals);
        Ok(())
    }
}

/// Distributed LDA XC of the replicated density `n` when the functional and
/// options allow it (LDA PZ or PW, no core density, no density floor, modified
/// PW92 prefactor); None means: use the replicated evaluation. Collective.
pub(crate) fn xc<K: Collective>(
    comm: &K,
    options: &ScfOptions,
    n: &Array3<f64>,
    core: Option<&Array3<f64>>,
    floor: Option<f64>,
    volume: f64,
) -> Option<Result<xc::XcGrid, PhysicsError>> {
    if core.is_some()
        || floor.is_some()
        || options.pw92_prefactor != Pw92Prefactor::Modified
        || !matches!(options.functional, Functional::LdaPz | Functional::LdaPw)
        || n.as_slice().is_none()
    {
        return None;
    }
    Some(xc_lda(comm, options.functional, n, volume))
}

fn xc_lda<K: Collective>(
    comm: &K,
    functional: Functional,
    n: &Array3<f64>,
    volume: f64,
) -> Result<xc::XcGrid, PhysicsError> {
    let start = dft_numerics::stage_clock::start();
    if !(volume.is_finite() && volume > 0.) || n.is_empty() {
        return Err(PhysicsError::Shape);
    }
    let all = n.as_slice().ok_or(PhysicsError::Shape)?;
    let counts = block_counts(all.len(), comm.size());
    let displs = block_displacements(&counts);
    let (lo, len) = (displs[comm.rank()], counts[comm.rank()]);
    let local = &all[lo..lo + len];
    let mut potential_local = vec![0.; len];
    // [energy sum, integral sum, failure flag] of this rank's points.
    let mut sums = [0.; 3];
    let ok = (|| -> Result<(), PhysicsError> {
        if len == 0 {
            return Ok(());
        }
        match functional {
            Functional::LdaPz => {
                for (v, &n) in potential_local.iter_mut().zip(local) {
                    let p = xc::lda_pz(n)?;
                    let n = n.max(0.0);
                    *v = p.v_xc();
                    sums[0] += n * p.eps_xc();
                    sums[1] += n * *v;
                }
            }
            _ => {
                // Volume = point count makes the library weight exactly 1.
                let block = Array3::from_shape_vec((1, 1, len), local.to_vec())
                    .map_err(|_| PhysicsError::Shape)?;
                let r = dft_physics::spin::unpolarized_pw(&block, len as f64)?;
                potential_local.copy_from_slice(r.potential.as_slice().unwrap());
                sums[0] = r.energy.0;
                sums[1] = r.potential_integral.0;
            }
        }
        Ok(())
    })();
    if ok.is_err() {
        sums = [0., 0., 1.];
    }
    let mut total = [0.; 3];
    comm.tree_sum_f64(&sums, &mut total)
        .map_err(|e| PhysicsError::Xc(e.to_string()))?;
    if total[2] != 0. {
        // Same error on every rank (the local one where it occurred).
        return Err(ok
            .err()
            .unwrap_or_else(|| PhysicsError::Xc("invalid density on another rank".into())));
    }
    let mut potential = Array3::zeros(n.raw_dim());
    comm.allgatherv_f64(&potential_local, &counts, potential.as_slice_mut().unwrap())
        .map_err(|e| PhysicsError::Xc(e.to_string()))?;
    let weight = volume / all.len() as f64;
    let (energy, integral) = (total[0] * weight, total[1] * weight);
    if !energy.is_finite() || !integral.is_finite() || potential.iter().any(|v| !v.is_finite()) {
        return Err(PhysicsError::Unrepresentable);
    }
    add(Stage::Xc, start);
    count(Count::Xc);
    Ok(xc::XcGrid {
        potential,
        energy: Hartree(energy),
        potential_integral: Hartree(integral),
    })
}

#[derive(Clone, Copy)]
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
}
const COUNT_NAMES: [&str; 23] = [
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
];
static COUNTS: [AtomicUsize; 23] = [const { AtomicUsize::new(0) }; 23];
pub(crate) fn count_atomic_block_fallback() {
    count(Count::AtomicBlockFallbacks);
}
pub(crate) fn count_max(c: Count, n: usize) {
    COUNTS[c as usize].fetch_max(n, Ordering::Relaxed);
}
pub(crate) fn count_n(c: Count, n: usize) {
    COUNTS[c as usize].fetch_add(n, Ordering::Relaxed);
}
pub(crate) fn count(c: Count) {
    COUNTS[c as usize].fetch_add(1, Ordering::Relaxed);
}

/// Engine-level M2d counters and stage seconds of this process (counters end
/// with one `calls.<stage>` entry per stage; seconds 0 with
/// `BAANDI_STAGE_TIMERS=off`).
pub fn stats() -> dft_physics::distributed_real_gamma::M2Counters {
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
            .collect(),
        STAGE_NAMES
            .iter()
            .zip(&STAGE_NANOS)
            .map(|(n, c)| (*n, c.load(Ordering::Relaxed) as f64 * 1e-9))
            .collect(),
    )
}
