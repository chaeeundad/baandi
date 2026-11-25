//! Row-block distributed real Davidson (ADR 0041, M1c). Interface in and out is
//! replicated (full initial vectors, full diagonal, full eigenvectors on every
//! rank); inside, the basis and the products are plane-wave row blocks, the dense
//! algebra is local GEMM plus the deterministic tree sum, and H is applied to
//! column pairs owned by one rank after gathering the correction block. With one
//! rank the caller should use the threads-only `RealDavidson` (bit-identical by
//! construction); this solver is used for two or more ranks.
use crate::real_davidson::lane_profile;
use crate::{
    BandTolerances, RealLinearOperator, SolverError,
    real_davidson::{RealEigenSolution, free_electron_initial_guess},
};
use dft_mpi::{Collective, block_counts, block_displacements};
use dft_numerics::{
    linalg::DEPENDENCE_TOLERANCE,
    pool::WorkerPool,
    real::{RealEighWorkspace, RealGemmWorkspace, RealMatrix},
};

/// A real symmetric operator whose basis rows are distributed over the ranks of a
/// communicator (M2c). Local row `l` of this rank is global (ORC) row
/// `local_rows()[l]`; the rank partitions must cover `0..global_len()` exactly once.
pub trait DistributedRealOperator {
    fn global_len(&self) -> usize;
    fn local_rows(&self) -> &[usize];
    fn worker_pool(&self) -> Option<std::sync::Arc<WorkerPool>> {
        None
    }
    fn inner_threads(&self) -> usize {
        1
    }
    /// y_local = H x_local for the local rows; every rank calls with the same column
    /// count (the call contains collectives).
    fn apply_block<C: Collective>(
        &mut self,
        comm: &C,
        x_local: &RealMatrix,
        y_local: &mut RealMatrix,
    ) -> Result<(), SolverError>;
}

mod workspace;
pub use workspace::{
    DistributedDavidsonWorkspace, LeanAlgebraStats, ScratchAllocationStats, WorkspaceStats,
    lean_algebra_stats, release_workspace, scratch_allocation_stats, set_iteration_probe,
    set_scratch_allocation_profile, start_gram_stats, workspace_stats,
};

fn comm_error(e: dft_mpi::MpiError) -> SolverError {
    SolverError::Application(e.to_string())
}

use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

/// `solve_start`: workspace path, from entry to the first iteration (row
/// extraction of the start, kept-start check, initial orthonormalization; the
/// initial apply is in `apply_total`). `solve_return`: the kept return and the
/// N x B gather of the eigenvectors at a converged return (MPI-SCALING-GAP).
/// `restart` (BREAKTHROUGH-0-INSTRUMENT): the restart block of the workspace
/// paths (Gram defect of X, contraction copy or re-orthonormalization and
/// re-application of the basis; the lean contraction's basis replacement).
/// Its rotation is counted in `rotate`.
///
/// BREAKTHROUGH-4-OTHER-AND-INIT (finer exclusive labels of the code the
/// earlier stages did not cover): `setup` (the per-solve row-index exchange of
/// `solve_rows_impl` and the workspace preparation before `solve_start`),
/// `precondition` (convergence flags, selected-root bookkeeping, coefficient
/// copies and the preconditioned correction block), `append` (copies of the new
/// columns into the basis and product storage and the finiteness scan of every
/// applied block). `wait` is NOT exclusive: with `BAANDI_PHASE_WAIT=on` a timed
/// barrier before every collective of the solver measures the time spent
/// waiting for the last rank (inside whichever stage issued the collective).
/// `solver_total` is the whole solver call (the container of every other
/// stage; solve time outside it belongs to the engine's interface).
pub const STAGES: [&str; 16] = [
    "reduced_gram",
    "reduced_eigen",
    "rotate",
    "residual",
    "orthonormalize",
    "apply_total",
    "apply_operator",
    "other_solve",
    "solve_start",
    "solve_return",
    "restart",
    "setup",
    "precondition",
    "append",
    "wait",
    "solver_total",
];
/// Stage indices of the BREAKTHROUGH-4 labels (see `STAGES`).
pub(crate) const STAGE_SETUP: usize = 11;
pub(crate) const STAGE_PRECONDITION: usize = 12;
pub(crate) const STAGE_APPEND: usize = 13;
const STAGE_WAIT: usize = 14;
const STAGE_TOTAL: usize = 15;
static STAGE_NANOS: [AtomicU64; 16] = [const { AtomicU64::new(0) }; 16];
static STAGE_CALLS: [AtomicU64; 16] = [const { AtomicU64::new(0) }; 16];

// Solves that reused the cached layout, solves that exchanged and stored it.
static LAYOUT_COUNTS: [AtomicU64; 2] = [const { AtomicU64::new(0) }; 2];
/// (solves that reused the cached layout, solves that exchanged and stored it)
/// of the retired layout-cache experiment (always zero).
pub fn layout_cache_stats() -> (u64, u64) {
    (
        LAYOUT_COUNTS[0].load(Ordering::Relaxed),
        LAYOUT_COUNTS[1].load(Ordering::Relaxed),
    )
}

/// `BAANDI_PHASE_WAIT=on` (read once): a timed barrier before every collective
/// of the distributed solver (diagnostic; timing only, the arithmetic is
/// unchanged).
pub fn wait_probe() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| dft_numerics::switches::on(dft_numerics::switches::Id::PHASE_WAIT))
}

/// Wall seconds of this process in each stage of the distributed solver (cheap
/// per-stage accumulation, no I/O in the loop; read once for the run report).
/// All 0 with `BAANDI_STAGE_TIMERS=off` (count-only; `stage_calls` still counts).
pub fn stage_seconds() -> Vec<(&'static str, f64)> {
    (0..STAGES.len())
        .map(|i| {
            (
                STAGES[i],
                STAGE_NANOS[i].load(Ordering::Relaxed) as f64 * 1e-9,
            )
        })
        .collect()
}
/// Entries into each stage (both time modes), in `STAGES` order.
pub fn stage_calls() -> Vec<(&'static str, u64)> {
    (0..STAGES.len())
        .map(|i| (STAGES[i], STAGE_CALLS[i].load(Ordering::Relaxed)))
        .collect()
}
/// F21: `start` from `dft_numerics::stage_clock::start` (`None` count-only).
pub(crate) fn stage_add(i: usize, start: Option<Instant>) {
    STAGE_CALLS[i].fetch_add(1, Ordering::Relaxed);
    let ns = dft_numerics::stage_clock::elapsed_ns(start);
    if start.is_some() {
        STAGE_NANOS[i].fetch_add(ns, Ordering::Relaxed);
    }
    // Per-thread mirror for the engine's phase profile (off: one relaxed load).
    dft_numerics::stage_clock::local::add(dft_numerics::stage_clock::local::DIST_STAGE + i, ns);
    // Serial-stage accumulators of `BAANDI_DAVIDSON_STAGES=on` (results.json
    // davidson_stages.by_run_stage.*.stage_s; off: one relaxed load).
    if dft_numerics::work_timers::enabled() {
        crate::davidson_stages::add_distributed(i, ns);
    }
}

// Reduced eigenproblem accounting (docs/performance/mpi-reduced-eigen-analysis.md):
// one histogram slot per matrix order (orders above the last slot are counted in
// it), plus the number of solves by kind. Process-wide, no allocation.
const ORDER_SLOTS: usize = 4096;
static SOLVE_ORDERS: [AtomicU64; ORDER_SLOTS] = [const { AtomicU64::new(0) }; ORDER_SLOTS];
// Sum of m^3 (f64 bits), single-column appends, bordered updates, root broadcasts.
static SOLVE_COUNTS: [AtomicU64; 4] = [const { AtomicU64::new(0) }; 4];

/// Distribution of the reduced eigenproblems solved by this process: orders,
/// cubic work and how each was solved. `p90` and `max` are exact orders.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ReducedSolveStats {
    pub bordered: u64,
    /// Reduced eigenproblems (one per Davidson iteration).
    pub solves: u64,
    pub mean_order: f64,
    pub p90_order: usize,
    pub max_order: usize,
    /// Sum of m^3 over all solves (the dense EVD cost model).
    pub sum_cubed: f64,
    /// Iterations whose reduced matrix grew by exactly one column.
    pub single_column: u64,
    /// Retired root-broadcast counter (always zero).
    pub broadcast: u64,
}

/// The reduced-solve distribution accumulated so far by this process.
pub fn reduced_solve_stats() -> ReducedSolveStats {
    let counts: Vec<u64> = SOLVE_ORDERS
        .iter()
        .map(|c| c.load(Ordering::Relaxed))
        .collect();
    let solves: u64 = counts.iter().sum();
    let mut sum = 0.0;
    let mut max_order = 0;
    for (m, &c) in counts.iter().enumerate() {
        sum += m as f64 * c as f64;
        if c > 0 {
            max_order = m;
        }
    }
    let target = (solves as f64 * 0.9).ceil() as u64;
    let mut seen = 0;
    let mut p90_order = 0;
    for (m, &c) in counts.iter().enumerate() {
        seen += c;
        if c > 0 && seen >= target {
            p90_order = m;
            break;
        }
    }
    let load = |i: usize| SOLVE_COUNTS[i].load(Ordering::Relaxed);
    ReducedSolveStats {
        solves,
        mean_order: if solves > 0 { sum / solves as f64 } else { 0.0 },
        p90_order,
        max_order,
        sum_cubed: f64::from_bits(load(0)),
        single_column: load(1),
        bordered: load(2),
        broadcast: load(3),
    }
}

fn record_solve(m: usize, fresh: usize) {
    SOLVE_ORDERS[m.min(ORDER_SLOTS - 1)].fetch_add(1, Ordering::Relaxed);
    let cubed = (m as f64).powi(3);
    // f64 accumulation through the bit pattern (compare-and-swap loop).
    let _ = SOLVE_COUNTS[0].try_update(Ordering::Relaxed, Ordering::Relaxed, |b| {
        Some((f64::from_bits(b) + cubed).to_bits())
    });
    if fresh == 1 {
        SOLVE_COUNTS[1].fetch_add(1, Ordering::Relaxed);
    }
}

// Ritz rotations V C and (H V) C: rotation calls (pairs), output columns,
// flops (2 n m c per product, both products counted). Process-wide.
static ROTATION_COUNTS: [AtomicU64; 3] = [const { AtomicU64::new(0) }; 3];

/// Ritz rotation work of the distributed solver (review F05 accounting).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct RotationStats {
    /// Rotation pairs (V C and (H V) C together).
    pub rotations: u64,
    /// Output columns of those pairs (per product).
    pub columns: u64,
    /// GEMM flops of both products, 2 n_local m c each.
    pub flops: u64,
}

/// The rotation work accumulated so far by this process.
pub fn rotation_stats() -> RotationStats {
    let load = |i: usize| ROTATION_COUNTS[i].load(Ordering::Relaxed);
    RotationStats {
        rotations: load(0),
        columns: load(1),
        flops: load(2),
    }
}

fn rotation_add(nl: usize, m: usize, cols: usize) {
    ROTATION_COUNTS[0].fetch_add(1, Ordering::Relaxed);
    ROTATION_COUNTS[1].fetch_add(cols as u64, Ordering::Relaxed);
    ROTATION_COUNTS[2].fetch_add(4 * (nl * m * cols) as u64, Ordering::Relaxed);
}

// ---- F08: reduction counters (count only; order and semantics unchanged) ----
// Every tree sum issued inside a distributed solve passes through `Counted`,
// which adds one call and the payload length to the histogram of the site the
// solver last declared on this thread (`site`). No timing, no allocation.
const REDUCTION_SITES: [&str; 11] = [
    "reduced_gram",
    "residual_norms",
    "ortho_normalize",
    "ortho_projection",
    "panel_cross",
    "panel_dot",
    "panel_norm",
    "restart_defect",
    "apply",
    "setup",
    "lean_border",
];
/// Payload buckets (doubles per call): 1, 2-8, 9-64, 65-512, 513-4096, above.
pub const REDUCTION_BUCKETS: [&str; 6] = ["1", "2-8", "9-64", "65-512", "513-4096", ">4096"];
static REDUCTION_CALLS: [[AtomicU64; 6]; 11] = [const { [const { AtomicU64::new(0) }; 6] }; 11];
static REDUCTION_PAYLOAD: [AtomicU64; 11] = [const { AtomicU64::new(0) }; 11];
static REDUCTION_SOLVES: AtomicU64 = AtomicU64::new(0);
thread_local! {
    static SITE: std::cell::Cell<usize> = const { std::cell::Cell::new(9) };
}

/// Declare the reduction site of the following tree sums on this thread.
pub(crate) fn site(s: Site) {
    SITE.with(|c| c.set(s as usize));
}

#[derive(Clone, Copy)]
pub(crate) enum Site {
    ReducedGram = 0,
    ResidualNorms = 1,
    OrthoNormalize = 2,
    OrthoProjection = 3,
    PanelCross = 4,
    PanelDot = 5,
    PanelNorm = 6,
    RestartDefect = 7,
    Apply = 8,
    Setup = 9,
    /// MPI-OVERLAP-DAVIDSON: the one tree sum per iteration of the lean
    /// algebra (S border, H border and the selected residual norms).
    LeanBorder = 10,
}

fn bucket(len: usize) -> usize {
    match len {
        0 | 1 => 0,
        2..=8 => 1,
        9..=64 => 2,
        65..=512 => 3,
        513..=4096 => 4,
        _ => 5,
    }
}

/// Tree-sum calls and payload of one site of the distributed solver.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReductionSite {
    pub site: &'static str,
    /// Calls per payload bucket (`REDUCTION_BUCKETS`).
    pub histogram: [u64; 6],
    pub calls: u64,
    /// Doubles reduced in total.
    pub payload: u64,
}

/// Reduction counts of this process (summed over the simulated ranks of a test
/// process) by site, and the number of distributed solves they cover.
pub fn reduction_stats() -> (u64, Vec<ReductionSite>) {
    let sites = REDUCTION_SITES
        .iter()
        .enumerate()
        .map(|(i, &name)| {
            let histogram: [u64; 6] =
                std::array::from_fn(|b| REDUCTION_CALLS[i][b].load(Ordering::Relaxed));
            ReductionSite {
                site: name,
                histogram,
                calls: histogram.iter().sum(),
                payload: REDUCTION_PAYLOAD[i].load(Ordering::Relaxed),
            }
        })
        .collect();
    (REDUCTION_SOLVES.load(Ordering::Relaxed), sites)
}

/// The solver's view of the communicator: forwards every call unchanged and
/// counts the tree sums by site.
struct Counted<'a, C: Collective> {
    inner: &'a C,
}
impl<C: Collective> Counted<'_, C> {
    /// `BAANDI_PHASE_WAIT=on`: a timed barrier before the collective, so the
    /// time spent waiting for the last rank is reported as `wait` and the
    /// collective's own clock measures the transfer (one load when off).
    fn probe(&self) -> Result<(), dft_mpi::MpiError> {
        if wait_probe() && self.inner.size() > 1 {
            let t = dft_numerics::stage_clock::start();
            self.inner.barrier()?;
            stage_add(STAGE_WAIT, t);
        }
        Ok(())
    }
}
impl<C: Collective> Collective for Counted<'_, C> {
    fn rank(&self) -> usize {
        self.inner.rank()
    }
    fn size(&self) -> usize {
        self.inner.size()
    }
    fn barrier(&self) -> Result<(), dft_mpi::MpiError> {
        self.inner.barrier()
    }
    fn allgatherv_f64(
        &self,
        local: &[f64],
        counts: &[usize],
        out: &mut [f64],
    ) -> Result<(), dft_mpi::MpiError> {
        self.probe()?;
        self.inner.allgatherv_f64(local, counts, out)
    }
    fn tree_sum_f64(&self, local: &[f64], out: &mut [f64]) -> Result<(), dft_mpi::MpiError> {
        let s = SITE.with(|c| c.get());
        let _lane = lane_profile::span(
            lane_profile::Op::TreeSum,
            local.len(),
            self.inner.size(),
            0,
            s,
        );
        REDUCTION_CALLS[s][bucket(local.len())].fetch_add(1, Ordering::Relaxed);
        REDUCTION_PAYLOAD[s].fetch_add(local.len() as u64, Ordering::Relaxed);
        self.probe()?;
        dft_numerics::stage_clock::arrival::site(401 + s);
        self.inner.tree_sum_f64(local, out)
    }
    fn bcast_f64(&self, data: &mut [f64], root: usize) -> Result<(), dft_mpi::MpiError> {
        self.probe()?;
        self.inner.bcast_f64(data, root)
    }
    fn alltoall_f64(
        &self,
        send: &[f64],
        block: usize,
        recv: &mut [f64],
    ) -> Result<(), dft_mpi::MpiError> {
        self.probe()?;
        self.inner.alltoall_f64(send, block, recv)
    }
}

/// Solver-only agreement: the operator uses the underlying communicator, then
/// records a local failure. On failure only, a quiet NaN replaces the first
/// scalar of the next nonempty solver sum. IEEE addition propagates it to
/// every rank; no failed numerical payload is consumed. Ordinary inputs and
/// payload sizes are unchanged. The operator does not use this wrapper.
struct ApplyAgreement<'a, C: Collective> {
    inner: &'a C,
    failed: std::cell::Cell<bool>,
    scratch: workspace::AgreementScratch,
}
impl<C: Collective> Collective for ApplyAgreement<'_, C> {
    fn rank(&self) -> usize {
        self.inner.rank()
    }
    fn size(&self) -> usize {
        self.inner.size()
    }
    fn barrier(&self) -> Result<(), dft_mpi::MpiError> {
        self.inner.barrier()
    }
    fn allgatherv_f64(
        &self,
        local: &[f64],
        counts: &[usize],
        out: &mut [f64],
    ) -> Result<(), dft_mpi::MpiError> {
        self.inner.allgatherv_f64(local, counts, out)
    }
    fn bcast_f64(&self, data: &mut [f64], root: usize) -> Result<(), dft_mpi::MpiError> {
        self.inner.bcast_f64(data, root)
    }
    fn alltoall_f64(
        &self,
        send: &[f64],
        block: usize,
        recv: &mut [f64],
    ) -> Result<(), dft_mpi::MpiError> {
        self.inner.alltoall_f64(send, block, recv)
    }
    fn tree_sum_f64(&self, local: &[f64], out: &mut [f64]) -> Result<(), dft_mpi::MpiError> {
        if self.failed.get() && !local.is_empty() {
            let mut input = self.scratch.buffer.borrow_mut();
            input[..local.len()].copy_from_slice(local);
            input[0] = f64::NAN;
            self.inner.tree_sum_f64(&input[..local.len()], out)?;
        } else {
            self.inner.tree_sum_f64(local, out)?;
        }
        if out.first().is_some_and(|x| x.is_nan()) {
            return Err(dft_mpi::MpiError::State(
                "distributed operator failed on at least one rank",
            ));
        }
        Ok(())
    }
}

/// Row partition: `index[r]` lists the global row of every local row of rank `r`
/// (contiguous blocks for M1, the stick partition for M2c). Algebra helpers only
/// use the local row count; the partition matters at the replicated interface.
struct Layout {
    n: usize,
    rank: usize,
    counts: Vec<usize>,
    index: Vec<Vec<usize>>,
}
impl Layout {
    fn nl(&self) -> usize {
        self.counts[self.rank]
    }
    fn contiguous(n: usize, rank: usize, ranks: usize) -> Self {
        let counts = block_counts(n, ranks);
        let displs = block_displacements(&counts);
        let index = (0..ranks)
            .map(|r| (displs[r]..displs[r] + counts[r]).collect())
            .collect();
        Self {
            n,
            rank,
            counts,
            index,
        }
    }
}

/// Rows of a replicated matrix owned by this rank.
fn local_rows(layout: &Layout, full: &RealMatrix) -> Result<RealMatrix, SolverError> {
    let mine = &layout.index[layout.rank];
    let mut data = Vec::with_capacity(mine.len() * full.cols());
    for column in full.as_slice().chunks(full.rows().max(1)).take(full.cols()) {
        data.extend(mine.iter().map(|&r| column[r]));
    }
    Ok(RealMatrix::new(mine.len(), full.cols(), data)?)
}

/// Gather the full n x w matrix on every rank from the row blocks (one allgatherv).
fn gather_full<C: Collective>(
    comm: &C,
    layout: &Layout,
    local: &RealMatrix,
) -> Result<RealMatrix, SolverError> {
    let w = local.cols();
    let counts: Vec<usize> = layout.counts.iter().map(|c| c * w).collect();
    let mut flat = vec![0.0; layout.n * w];
    comm.allgatherv_f64(local.as_slice(), &counts, &mut flat)
        .map_err(comm_error)?;
    let mut data = vec![0.0; layout.n * w];
    let mut at = 0;
    for (r, &rc) in layout.counts.iter().enumerate() {
        for j in 0..w {
            let src = &flat[at + j * rc..at + (j + 1) * rc];
            for (&g, &v) in layout.index[r].iter().zip(src) {
                data[j * layout.n + g] = v;
            }
        }
        at += rc * w;
    }
    Ok(RealMatrix::new(layout.n, w, data)?)
}

/// H applied to the columns of a row-block matrix: gather the block, apply H to
/// the column pairs this rank owns (the pair partition, and therefore every
/// column's result, does not depend on the rank count), gather the products and
/// keep the local rows.
fn apply_distributed<C: Collective>(
    comm: &C,
    layout: &Layout,
    op: &dyn RealLinearOperator,
    x_local: &RealMatrix,
) -> Result<RealMatrix, SolverError> {
    let apply_start = dft_numerics::stage_clock::start();
    let w = x_local.cols();
    let n = layout.n;
    let full = gather_full(comm, layout, x_local)?;
    let pairs = w.div_ceil(2);
    let pair_counts = block_counts(pairs, comm.size());
    let pair_displs = block_displacements(&pair_counts);
    let first = (2 * pair_displs[comm.rank()]).min(w);
    let last = (2 * (pair_displs[comm.rank()] + pair_counts[comm.rank()])).min(w);
    let owned = last - first;
    let mut mine = RealMatrix::zeros(n, owned)?;
    if owned > 0 {
        let x = RealMatrix::new(n, owned, full.as_slice()[first * n..last * n].to_vec())?;
        let op_start = dft_numerics::stage_clock::start();
        op.apply_block(&x, &mut mine)?;
        stage_add(6, op_start);
        if mine.rows() != n
            || mine.cols() != owned
            || mine.as_slice().iter().any(|v| !v.is_finite())
        {
            return Err(SolverError::Application(
                "operator returned an invalid block".into(),
            ));
        }
    }
    let col_counts: Vec<usize> = (0..comm.size())
        .map(|r| {
            let f = (2 * pair_displs[r]).min(w);
            let l = (2 * (pair_displs[r] + pair_counts[r])).min(w);
            (l - f) * n
        })
        .collect();
    let mut products = vec![0.0; n * w];
    comm.allgatherv_f64(mine.as_slice(), &col_counts, &mut products)
        .map_err(comm_error)?;
    let out = local_rows(layout, &RealMatrix::new(n, w, products)?);
    stage_add(5, apply_start);
    out
}

/// Largest absolute difference between the solver's `initial` and the kept X
/// (local rows) for which a kept start is taken (as the serial
/// `real_davidson::KEPT_START_TOLERANCE`: the pack/unpack round trip of the
/// complex warm start changes X by rounding only).
pub const KEPT_START_TOLERANCE: f64 = 1e-14;

/// Start block of a stick-row solve (MPI-LAZY-ORBITALS): replicated global
/// rows (N x w, the same on every rank) or this rank's rows only (n_loc x w, in
/// the operator's `local_rows` order).
#[derive(Clone, Copy, Debug)]
pub enum DistributedStart<'a> {
    Global(&'a RealMatrix),
    Local(&'a RealMatrix),
}

/// Kept start of the stick-row solver (MPI-KEPT-START, `BAANDI_M2_REUSE`): the
/// H images of the latest Ritz-change return on this rank's rows, owned by the
/// caller across solves (the X itself is the caller's `x_local`). All flags
/// are replicated: every rank sets and reads them under the same conditions.
#[derive(Debug)]
pub struct DistributedKeptStart {
    hx: RealMatrix,
    keep: bool,
    stored: bool,
    armed: bool,
    taken: bool,
    declined: bool,
}

impl Default for DistributedKeptStart {
    fn default() -> Self {
        Self::new()
    }
}

impl DistributedKeptStart {
    pub fn new() -> Self {
        Self {
            hx: RealMatrix::new(0, 0, Vec::new()).expect("empty matrix"),
            keep: false,
            stored: false,
            armed: false,
            taken: false,
            declined: false,
        }
    }
    /// Keep H X at the next Ritz-change return (`stored` afterwards).
    pub fn set_keep(&mut self, keep: bool) {
        self.keep = keep;
    }
    /// H X of the latest keeping return (local rows x bands), if stored.
    pub fn stored_mut(&mut self) -> Option<&mut RealMatrix> {
        self.stored.then_some(&mut self.hx)
    }
    pub fn is_stored(&self) -> bool {
        self.stored
    }
    /// Let the next solve try the kept start (the caller has brought H X to the
    /// current operator). Without a stored return this does nothing.
    pub fn arm(&mut self) {
        self.armed = self.stored;
    }
    /// Forget the stored return (rebuilt operator, failed solve, checkpoint).
    pub fn invalidate(&mut self) {
        self.stored = false;
        self.armed = false;
    }
    /// Whether the latest solve started from the kept X and H X.
    pub fn taken(&self) -> bool {
        self.taken
    }
    /// Whether the latest solve was armed but declined the kept start
    /// (initial differs from X, or the Gram defect of X exceeds 1e-12).
    pub fn declined(&self) -> bool {
        self.declined
    }
    /// Whether the next solve is armed to try the kept start.
    pub fn armed(&self) -> bool {
        self.armed
    }
    /// H X of the latest keeping return, if stored (read).
    pub fn stored(&self) -> Option<&RealMatrix> {
        self.stored.then_some(&self.hx)
    }
    /// A kept start holding `hx` as its stored H X (snapshot replay,
    /// BREAKTHROUGH-0-INSTRUMENT): stored, not armed, keeping off.
    pub fn with_stored(hx: RealMatrix) -> Self {
        Self {
            hx,
            keep: false,
            stored: true,
            armed: false,
            taken: false,
            declined: false,
        }
    }
}

/// Run-time choices of the distributed solver (MPI entry cost,
/// docs/performance/mpi-entry-cost-analysis.md). `from_env` reads the switches;
/// Oracle algorithm choices are frozen at their defaults (R11b). Every remaining choice uses
/// replicated data only, so all ranks take the same branches.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct DistributedTuning {
    /// `BAANDI_M2_DGEMM=on`: rotations and basis projections through
    /// `multiply_slices_wide_into` (per-panel dgemm at any width on MKL builds;
    /// bit-identical on pure builds).
    pub wide_gemm: bool,
    /// DGKS keep ratio of the serial order (`BAANDI_ORTHO`); `None` = two passes.
    pub dgks: Option<f64>,
    /// Panel rule of the serial order: `Some(0)` = ORTHO-BLOCK shape rule on the
    /// widest local row count, `Some(w)` = fixed width w (MGS-B), `None` = the
    /// whole block as one panel (the classic sweep).
    pub panel: Option<usize>,
    /// `BAANDI_RESTART_KEEP` (ticket BREAKTHROUGH-1-RESTART, R1-b): previous
    /// directions kept at a contraction of the lean algebra path
    /// (`RestartKeep`; default off = the plain restart, bit for bit).
    pub restart_keep: crate::real_davidson::RestartKeep,
}

/// In-panel reduction scheme of the distributed panel Gram-Schmidt
/// (shared kernel; legacy oracle always uses Fused). The
/// in-panel dot and norm reductions are about 88 percent of all reductions of the
/// distributed solver and carry at most 64 doubles each (latency bound, ledger
/// M2-REDUCTION-COUNTS-F08).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PanelGs {
    /// Two dot reductions and one norm reduction per column (the M1/M2 code).
    #[default]
    Classic,
    /// Two passes; the norm of the once-projected column rides with the second
    /// pass's dot products and the final norm is the Pythagorean update
    /// `|v1|^2 - |s2|^2` (explicit norm reduction when that update would cancel
    /// more than a factor 4). Two reductions per column; rounding class.
    Pythag,
    /// The input norm rides with the first pass's dot products; the second pass
    /// is skipped when the DGKS rule holds for the column (`|v1| >= keep |v|`,
    /// keep = the `BAANDI_ORTHO` ratio the serial solver applies against the
    /// basis), else `Pythag`'s second pass. About one reduction per column;
    /// rounding class.
    Fused,
}
impl DistributedTuning {
    /// Shared switches and the frozen backend-dependent serial panel rule.
    pub fn from_env() -> Self {
        use dft_numerics::switches::{Id, on};
        let panel = Some(if dft_numerics::MKL_BACKEND {
            crate::real_davidson::BLOCK_COLUMNS
        } else {
            0
        });
        Self {
            wide_gemm: on(Id::M2_DGEMM),
            dgks: crate::real_davidson::dgks_keep(),
            panel,
            restart_keep: crate::real_davidson::RestartKeep::from_env(),
        }
    }
    /// The tuning of the calling thread: the test override when set
    /// (`set_tuning_override`), else the switches (`from_env`).
    pub fn current() -> Self {
        TUNING_OVERRIDE
            .with(|t| t.get())
            .unwrap_or_else(Self::from_env)
    }
    /// Panel width for local blocks of at most `max_rows` rows (same on every rank).
    fn panel_width(&self, max_rows: usize) -> usize {
        match self.panel {
            Some(0) => crate::real_davidson::ortho_panel_width(max_rows),
            Some(w) => w,
            None => usize::MAX,
        }
    }
}

thread_local! {
    // Test-only per-thread override of the switches (simulated ranks run on their
    // own threads, so each rank sets it inside its closure).
    static TUNING_OVERRIDE: std::cell::Cell<Option<DistributedTuning>> =
        const { std::cell::Cell::new(None) };
}

/// Override the distributed solver's switches for the calling thread (tests);
/// `None` restores the environment.
pub fn set_tuning_override(tuning: Option<DistributedTuning>) {
    TUNING_OVERRIDE.with(|t| t.set(tuning));
}

/// Orthogonalise `v` (local rows) against the `q.len() / rows` orthonormal
/// in-panel columns `q` (column-major) and return its remaining global norm,
/// with the reductions of `mode` (`PanelGs`). `Classic` is the arithmetic of the
/// earlier panel code bit for bit. Every branch depends on reduced values only,
/// so all ranks take the same branches and enter the same reductions.
#[allow(clippy::too_many_arguments)]
pub(crate) fn panel_column<C: Collective>(
    comm: &C,
    q: &[f64],
    rows: usize,
    k: usize,
    v: &mut [f64],
    local: &mut Vec<f64>,
    sums: &mut Vec<f64>,
    mode: PanelGs,
    keep: Option<f64>,
) -> Result<f64, SolverError> {
    // The column count is replicated even when this rank owns zero rows.
    debug_assert_eq!(q.len(), rows * k);
    // Local dot products with the k columns, optionally followed by |v|^2.
    let dots = |v: &[f64], local: &mut Vec<f64>, norm: bool| {
        local.clear();
        for c in 0..k {
            let qc = &q[c * rows..(c + 1) * rows];
            local.push(qc.iter().zip(v).map(|(a, b)| a * b).sum::<f64>());
        }
        if norm {
            local.push(v.iter().map(|a| a * a).sum::<f64>());
        }
    };
    let subtract = |v: &mut [f64], sums: &[f64]| {
        for (c, s) in sums.iter().enumerate().take(k) {
            let qc = &q[c * rows..(c + 1) * rows];
            for (vi, qi) in v.iter_mut().zip(qc) {
                *vi -= s * qi;
            }
        }
    };
    let explicit_norm = |v: &[f64]| -> Result<f64, SolverError> {
        let local_norm = [v.iter().map(|a| a * a).sum::<f64>()];
        let mut total = [0.0];
        site(Site::PanelNorm);
        comm.tree_sum_f64(&local_norm, &mut total)
            .map_err(comm_error)?;
        Ok(total[0].sqrt())
    };
    let reduce = |local: &[f64], sums: &mut Vec<f64>| -> Result<(), SolverError> {
        sums.clear();
        sums.resize(local.len(), 0.0);
        site(Site::PanelDot);
        comm.tree_sum_f64(local, sums).map_err(comm_error)
    };
    if mode == PanelGs::Classic || k == 0 {
        if k > 0 {
            for _ in 0..2 {
                dots(v, local, false);
                reduce(local, sums)?;
                subtract(v, sums);
            }
        }
        return explicit_norm(v);
    }
    // First pass: the dots, plus the input norm when the column may stop here.
    let fused = mode == PanelGs::Fused && keep.is_some();
    dots(v, local, fused);
    reduce(local, sums)?;
    subtract(v, sums);
    if fused {
        let keep = keep.unwrap_or(0.0);
        let (n0, removed) = (sums[k], sum_squares(&sums[..k]));
        if n0 > 0.0 && n0 - removed >= keep * keep * n0 {
            return match pythagorean(n0, removed) {
                Some(len) => Ok(len),
                None => explicit_norm(v),
            };
        }
    }
    // Second pass with the norm of the once-projected column.
    dots(v, local, true);
    reduce(local, sums)?;
    let (n1, removed) = (sums[k], sum_squares(&sums[..k]));
    subtract(v, sums);
    match pythagorean(n1, removed) {
        Some(len) => Ok(len),
        None => explicit_norm(v),
    }
}

fn sum_squares(s: &[f64]) -> f64 {
    s.iter().map(|x| x * x).sum()
}

/// Norm of u - Q s (s = Q^T u, Q orthonormal) from |u|^2 and |s|^2, only when the
/// update cancels at most a factor 4 (relative rounding below about 4 eps);
/// `None` asks for an explicit norm reduction.
fn pythagorean(norm2: f64, removed: f64) -> Option<f64> {
    let left = norm2 - removed;
    (left >= 0.25 * norm2).then(|| left.sqrt())
}

fn column_of(x: &RealMatrix, j: usize) -> &[f64] {
    &x.as_slice()[j * x.rows()..(j + 1) * x.rows()]
}

/// Frozen replicated dense reduced eigenproblem.
fn solve_reduced(
    eig: &mut RealEighWorkspace,
    reduced: &RealMatrix,
    prefix: usize,
    values: &mut [f64],
    coefficients: &mut RealMatrix,
) -> Result<(), SolverError> {
    let _lane = lane_profile::span(
        lane_profile::Op::EigenStep,
        reduced.rows(),
        prefix,
        reduced.rows(),
        0,
    );

    let m = reduced.rows();
    let _full = lane_profile::span(lane_profile::Op::Eigen, m, m, m, 0);
    eig.solve_into(reduced, values, coefficients)?;
    Ok(())
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        (((z ^ (z >> 31)) >> 11) as f64) * (1. / ((1_u64 << 53) as f64)) - 0.5
    }
}

impl crate::RealDavidson {
    /// Distributed solve over the ranks of `comm` (two or more; with one rank use
    /// `solve_policy`/`solve_ritz_changes`). Every rank must call it with the same
    /// replicated arguments and returns the same replicated solution.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_distributed<C: Collective>(
        &self,
        comm: &C,
        op: &dyn RealLinearOperator,
        initial: Option<&RealMatrix>,
        diagonal: &[f64],
        policy: BandTolerances,
        ritz: bool,
        certify: bool,
        capacity_per_band: usize,
    ) -> Result<RealEigenSolution, SolverError> {
        let counted = Counted { inner: comm };
        let comm = &counted;
        REDUCTION_SOLVES.fetch_add(1, Ordering::Relaxed);
        site(Site::Setup);
        let n = op.dimension();
        if diagonal.len() != n {
            return Err(SolverError::InvalidParameter("bands or diagonal length"));
        }
        let layout = Layout::contiguous(n, comm.rank(), comm.size());
        let nl = layout.nl();
        let offset = layout.index[comm.rank()].first().copied().unwrap_or(0);
        let diag_local = &diagonal[offset..offset + nl];
        let tuning = DistributedTuning::current();
        {
            let mut apply = |x: &RealMatrix, y: &mut RealMatrix| -> Result<(), SolverError> {
                site(Site::Apply);
                let r = apply_distributed(comm, &layout, op, x)?;
                y.as_mut_slice().copy_from_slice(r.as_slice());
                Ok(())
            };
            workspace::solve_core_ws(
                self,
                comm,
                &layout,
                op.worker_pool(),
                op.inner_threads(),
                &mut apply,
                initial.map(DistributedStart::Global),
                diag_local,
                Some(diagonal),
                policy,
                ritz,
                certify,
                capacity_per_band,
                tuning,
                None,
                true,
            )
        }
    }

    /// Distributed solve directly on the row partition of a distributed operator
    /// (M2c): the basis and the products stay row-distributed, H is the operator's
    /// own distributed apply (no gather of the block, no product allgather), and
    /// the eigenvectors are gathered to the replicated interface only at the end,
    /// in global (ORC) row order. `initial` is the replicated start (global rows),
    /// `diagonal_local` the preconditioner diagonal on the operator's local rows.
    /// All ranks must call it together with the same replicated arguments.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_distributed_rows<C: Collective, O: DistributedRealOperator>(
        &self,
        comm: &C,
        op: &mut O,
        initial: Option<&RealMatrix>,
        diagonal_local: &[f64],
        policy: BandTolerances,
        ritz: bool,
        certify: bool,
        capacity_per_band: usize,
    ) -> Result<RealEigenSolution, SolverError> {
        self.solve_distributed_rows_tuned(
            comm,
            op,
            initial,
            diagonal_local,
            policy,
            ritz,
            certify,
            capacity_per_band,
            DistributedTuning::current(),
        )
    }

    /// `solve_distributed_rows` with explicit run-time choices (tests and A/B);
    /// the plain entry point reads them from the switches.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_distributed_rows_tuned<C: Collective, O: DistributedRealOperator>(
        &self,
        comm: &C,
        op: &mut O,
        initial: Option<&RealMatrix>,
        diagonal_local: &[f64],
        policy: BandTolerances,
        ritz: bool,
        certify: bool,
        capacity_per_band: usize,
        tuning: DistributedTuning,
    ) -> Result<RealEigenSolution, SolverError> {
        self.solve_rows_impl(
            comm,
            op,
            initial.map(DistributedStart::Global),
            diagonal_local,
            policy,
            ritz,
            certify,
            capacity_per_band,
            tuning,
            None,
            true,
        )
    }

    /// MPI-KEPT-START (`BAANDI_M2_REUSE=on`): `solve_distributed_rows_tuned` with
    /// a kept start. `x_local` receives the local rows of the returned
    /// eigenvectors (this rank's rows of `eigenvectors`, in local row order; no
    /// extraction from the gathered block needed). When `kept` is armed and holds
    /// a return, `x_local` must hold that return's X and `kept` its H X at the
    /// CURRENT operator (the caller adds (V_new - V_old) X); the solve then starts
    /// from basis X and products H X without orthonormalization or an initial
    /// apply, if (replicated decision) `initial` equals X on every rank to
    /// `KEPT_START_TOLERANCE` and the Gram defect of X is at most 1e-12 (the
    /// contraction rule of a restart). Otherwise it takes the plain start. A
    /// Ritz-change return with `kept.keep` stores H X (the return rotation of the
    /// H images). Needs the workspace path (`true` or
    /// `selected_roots`); otherwise the kept data is dropped and the plain solver
    /// runs.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_distributed_rows_kept<C: Collective, O: DistributedRealOperator>(
        &self,
        comm: &C,
        op: &mut O,
        initial: Option<&RealMatrix>,
        diagonal_local: &[f64],
        policy: BandTolerances,
        ritz: bool,
        certify: bool,
        capacity_per_band: usize,
        tuning: DistributedTuning,
        x_local: &mut RealMatrix,
        kept: &mut DistributedKeptStart,
    ) -> Result<RealEigenSolution, SolverError> {
        self.solve_rows_impl(
            comm,
            op,
            initial.map(DistributedStart::Global),
            diagonal_local,
            policy,
            ritz,
            certify,
            capacity_per_band,
            tuning,
            Some((x_local, kept)),
            true,
        )
    }

    /// MPI-LAZY-ORBITALS (`BAANDI_M2_LAZY=on`): `solve_distributed_rows_kept`
    /// with the start given either replicated (`DistributedStart::Global`) or as
    /// this rank's rows in the operator's local row order
    /// (`DistributedStart::Local`, no replicated block needed), and with the
    /// final N x B allgatherv optional. With `gather = false` on the workspace
    /// path the returned `eigenvectors` is a 0 x bands placeholder and the
    /// eigenvectors exist only as the local rows in `x_local` (every rank holds
    /// its own rows). On the non-workspace path a local start is gathered first
    /// and the full block is always returned; callers check
    /// `eigenvectors.rows()`. All ranks must pass the same variant and `gather`.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_distributed_rows_local<C: Collective, O: DistributedRealOperator>(
        &self,
        comm: &C,
        op: &mut O,
        initial: Option<DistributedStart<'_>>,
        diagonal_local: &[f64],
        policy: BandTolerances,
        ritz: bool,
        certify: bool,
        capacity_per_band: usize,
        tuning: DistributedTuning,
        x_local: &mut RealMatrix,
        kept: &mut DistributedKeptStart,
        gather: bool,
    ) -> Result<RealEigenSolution, SolverError> {
        self.solve_rows_impl(
            comm,
            op,
            initial,
            diagonal_local,
            policy,
            ritz,
            certify,
            capacity_per_band,
            tuning,
            Some((x_local, kept)),
            gather,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_rows_impl<C: Collective, O: DistributedRealOperator>(
        &self,
        comm: &C,
        op: &mut O,
        initial: Option<DistributedStart<'_>>,
        diagonal_local: &[f64],
        policy: BandTolerances,
        ritz: bool,
        certify: bool,
        capacity_per_band: usize,
        tuning: DistributedTuning,
        mut kept: Option<(&mut RealMatrix, &mut DistributedKeptStart)>,
        gather: bool,
    ) -> Result<RealEigenSolution, SolverError> {
        if let Some((_, k)) = kept.as_mut() {
            k.taken = false;
            k.declined = false;
        }
        let t_total = dft_numerics::stage_clock::start();
        let t_setup = dft_numerics::stage_clock::start();
        let counted = Counted { inner: comm };
        let comm = &counted;
        REDUCTION_SOLVES.fetch_add(1, Ordering::Relaxed);
        site(Site::Setup);
        let n = op.global_len();
        let mine: &[usize] = op.local_rows();
        let _lane_solve =
            lane_profile::solve(mine.len(), self.bands, capacity_per_band, self.tolerance);
        let invalid = diagonal_local.len() != mine.len()
            || mine.iter().any(|&r| r >= n)
            || matches!(initial, Some(DistributedStart::Local(x)) if x.rows() != mine.len());
        let ranks = comm.size();
        let layout: std::sync::Arc<Layout> = {
            // Setup (once per solve): every rank learns the row lists of all ranks.
            let mut counts_f = vec![0.0; ranks];
            comm.allgatherv_f64(
                &[if invalid { f64::NAN } else { mine.len() as f64 }],
                &vec![1; ranks],
                &mut counts_f,
            )
            .map_err(comm_error)?;
            if counts_f.iter().any(|v| v.is_nan()) {
                return Err(SolverError::InvalidParameter(
                    "invalid local rows, diagonal or start",
                ));
            }
            let counts: Vec<usize> = counts_f.iter().map(|&v| v as usize).collect();
            if counts.iter().sum::<usize>() != n {
                return Err(SolverError::InvalidParameter(
                    "row partition does not cover the basis",
                ));
            }
            let mut flat = vec![0.0; n];
            let mine_f: Vec<f64> = mine.iter().map(|&r| r as f64).collect();
            comm.allgatherv_f64(&mine_f, &counts, &mut flat)
                .map_err(comm_error)?;
            let mut at = 0;
            let index: Vec<Vec<usize>> = counts
                .iter()
                .map(|&c| {
                    let v = flat[at..at + c].iter().map(|&r| r as usize).collect();
                    at += c;
                    v
                })
                .collect();
            let layout = std::sync::Arc::new(Layout {
                n,
                rank: comm.rank(),
                counts,
                index,
            });
            {}
            layout
        };
        let cap = self
            .max_subspace
            .max(self.bands.saturating_mul(capacity_per_band))
            .min(n);
        let payload = 2 * cap * cap + 2 * cap + 1;
        let persistent = true;
        let agreement = ApplyAgreement {
            inner: comm,
            failed: std::cell::Cell::new(false),
            scratch: workspace::AgreementScratch::acquire(payload, persistent),
        };
        stage_add(STAGE_SETUP, t_setup);
        {
            let supplied = op.worker_pool();
            let threads = op.inner_threads();
            let mut apply = |x: &RealMatrix, y: &mut RealMatrix| -> Result<(), SolverError> {
                site(Site::Apply);
                let start = dft_numerics::stage_clock::start();
                y.as_mut_slice().fill(0.0);
                let failed = op.apply_block(comm, x, y).is_err();
                stage_add(5, start);
                stage_add(6, start);
                let scan = dft_numerics::stage_clock::start();
                let invalid = y.as_slice().iter().any(|v| !v.is_finite());
                stage_add(STAGE_APPEND, scan);
                if failed || invalid {
                    agreement.failed.set(true);
                    y.as_mut_slice().fill(0.0);
                }
                Ok(())
            };
            let result = workspace::solve_core_ws(
                self,
                &agreement,
                &layout,
                supplied,
                threads,
                &mut apply,
                initial,
                diagonal_local,
                None,
                policy,
                ritz,
                certify,
                capacity_per_band,
                tuning,
                kept,
                gather,
            );
            stage_add(STAGE_TOTAL, t_total);
            result
        }
    }
}

#[cfg(test)]
mod panel_gs_tests {
    use super::*;
    use dft_mpi::simulate;

    /// Rows of a global vector owned by `rank` (cyclic split).
    fn local(global: &[f64], rank: usize, ranks: usize) -> Vec<f64> {
        global.iter().skip(rank).step_by(ranks).copied().collect()
    }

    /// `panel_column` in every mode: an exactly dependent column is dropped (norm
    /// at the dependence floor); a nearly dependent one (relative 1e-7) and a
    /// generic one keep their norm to 1e-9 relative of the classic mode and the
    /// orthogonalised vectors agree to 1e-12; every rank returns the same bits.
    /// The fused mode takes one reduction for the generic column and two for the
    /// (nearly) dependent ones (the per-column DGKS rule sends them to the second
    /// pass).
    #[test]
    fn panel_column_modes_agree_and_drop_dependent_columns() {
        let n = 97;
        let k = 3;
        let mut q = vec![0.0; n * k];
        for c in 0..k {
            for i in (c * 30)..(c * 30 + 30) {
                q[c * n + i] = (1.0 + (i % 7) as f64).sqrt();
            }
            let s: f64 = q[c * n..(c + 1) * n]
                .iter()
                .map(|x| x * x)
                .sum::<f64>()
                .sqrt();
            q[c * n..(c + 1) * n].iter_mut().for_each(|x| *x /= s);
        }
        let generic: Vec<f64> = (0..n).map(|i| ((i * 37 % 11) as f64 - 5.0) * 0.1).collect();
        let span: Vec<f64> = (0..n)
            .map(|i| 0.3 * q[i] - 0.7 * q[n + i] + 0.2 * q[2 * n + i])
            .collect();
        let near: Vec<f64> = span
            .iter()
            .zip(&generic)
            .map(|(s, g)| s + 1e-7 * g)
            .collect();
        for ranks in [1usize, 2, 3] {
            let out = simulate(ranks, |comm| {
                let rows = local(&generic, comm.rank(), ranks).len();
                let mut ql = Vec::new();
                for c in 0..k {
                    ql.extend(local(&q[c * n..(c + 1) * n], comm.rank(), ranks));
                }
                let mut res = Vec::new();
                for mode in [PanelGs::Classic, PanelGs::Pythag, PanelGs::Fused] {
                    for vector in [&generic, &span, &near] {
                        let mut v = local(vector, comm.rank(), ranks);
                        let (mut l, mut s) = (Vec::new(), Vec::new());
                        let calls = std::cell::Cell::new(0u64);
                        let counting = CountingComm {
                            inner: comm,
                            calls: &calls,
                        };
                        let len = panel_column(
                            &counting,
                            &ql,
                            rows,
                            k,
                            &mut v,
                            &mut l,
                            &mut s,
                            mode,
                            Some(0.1),
                        )
                        .unwrap();
                        res.push((mode, len, v, calls.get()));
                    }
                }
                res
            });
            for rank_res in &out {
                for (a, b) in rank_res.iter().zip(&out[0]) {
                    assert_eq!(a.1.to_bits(), b.1.to_bits(), "ranks {ranks}");
                }
            }
            let r = &out[0];
            for (i, (mode, len, v, calls)) in r.iter().enumerate() {
                let c = &r[i % 3];
                if i % 3 == 1 {
                    assert!(*len <= DEPENDENCE_TOLERANCE, "{mode:?} dependent {len}");
                } else {
                    assert!((len - c.1).abs() <= 1e-9 * c.1, "{mode:?} {len} vs {}", c.1);
                    for (a, b) in v.iter().zip(&c.2) {
                        assert!((a - b).abs() <= 1e-12 * c.1, "{mode:?}");
                    }
                }
                let expect = match (mode, i % 3) {
                    (PanelGs::Classic, _) => 3,
                    (PanelGs::Pythag, _) => 2,
                    (PanelGs::Fused, 0) => 1,
                    (PanelGs::Fused, _) => 2,
                };
                // A cancelling Pythagorean update adds one explicit norm.
                assert!(
                    *calls == expect || (*calls == expect + 1 && i % 3 != 0),
                    "{mode:?} case {i}: {calls} reductions"
                );
            }
        }
    }

    struct CountingComm<'a, C> {
        inner: &'a C,
        calls: &'a std::cell::Cell<u64>,
    }
    impl<C: Collective> Collective for CountingComm<'_, C> {
        fn rank(&self) -> usize {
            self.inner.rank()
        }
        fn size(&self) -> usize {
            self.inner.size()
        }
        fn barrier(&self) -> Result<(), dft_mpi::MpiError> {
            self.inner.barrier()
        }
        fn allgatherv_f64(
            &self,
            local: &[f64],
            counts: &[usize],
            out: &mut [f64],
        ) -> Result<(), dft_mpi::MpiError> {
            self.inner.allgatherv_f64(local, counts, out)
        }
        fn tree_sum_f64(&self, local: &[f64], out: &mut [f64]) -> Result<(), dft_mpi::MpiError> {
            self.calls.set(self.calls.get() + 1);
            self.inner.tree_sum_f64(local, out)
        }
        fn bcast_f64(&self, data: &mut [f64], root: usize) -> Result<(), dft_mpi::MpiError> {
            self.inner.bcast_f64(data, root)
        }
        fn alltoall_f64(
            &self,
            send: &[f64],
            block: usize,
            recv: &mut [f64],
        ) -> Result<(), dft_mpi::MpiError> {
            self.inner.alltoall_f64(send, block, recv)
        }
    }
}
