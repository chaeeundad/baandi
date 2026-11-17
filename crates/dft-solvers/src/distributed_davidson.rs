//! Row-block distributed real Davidson (ADR 0041, M1c). Interface in and out is
//! replicated (full initial vectors, full diagonal, full eigenvectors on every
//! rank); inside, the basis and the products are plane-wave row blocks, the dense
//! algebra is local GEMM plus the deterministic tree sum, and H is applied to
//! column pairs owned by one rank after gathering the correction block. With one
//! rank the caller should use the threads-only `RealDavidson` (bit-identical by
//! construction); this solver is used for two or more ranks.
use crate::{
    BandTolerances, RealLinearOperator, SolverError,
    bordered::{Bordered, MAX_STREAK},
    distributed::{column_norms_squared, gram, project_out_with, rotate_local_with},
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
    DistributedDavidsonWorkspace, LeanAlgebraStats, WorkspaceStats, lean_algebra_stats,
    release_workspace, set_iteration_probe, start_gram_stats, workspace_stats,
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

thread_local! {
    /// This rank's local rows and the row layout of its previous solve
    /// (`BAANDI_M2_LAYOUT_CACHE=on`).
    static LAYOUT_CACHE: std::cell::RefCell<Option<(Vec<usize>, std::sync::Arc<Layout>)>> =
        const { std::cell::RefCell::new(None) };
    static LAYOUT_CACHE_OVERRIDE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}
// Solves that reused the cached layout, solves that exchanged and stored it.
static LAYOUT_COUNTS: [AtomicU64; 2] = [const { AtomicU64::new(0) }; 2];
/// `BAANDI_M2_LAYOUT_CACHE=on|off` (default off; per-thread override for
/// simulated ranks): reuse the row-index layout of the previous solve when
/// every rank's local rows are unchanged (bit-identical: same layout).
fn layout_cache() -> bool {
    LAYOUT_CACHE_OVERRIDE.with(|o| o.get()).unwrap_or_else(|| {
        static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
        *ON.get_or_init(|| dft_numerics::switches::on(dft_numerics::switches::Id::M2_LAYOUT_CACHE))
    })
}
/// Override `BAANDI_M2_LAYOUT_CACHE` on the calling thread (`None`: the switch).
pub fn set_layout_cache_override(on: Option<bool>) {
    LAYOUT_CACHE_OVERRIDE.with(|o| o.set(on));
    if on == Some(false) {
        LAYOUT_CACHE.with(|l| *l.borrow_mut() = None);
    }
}
/// (solves that reused the cached layout, solves that exchanged and stored it)
/// of `BAANDI_M2_LAYOUT_CACHE=on`, summed over the threads of this process.
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
    /// Reduced eigenproblems (one per Davidson iteration).
    pub solves: u64,
    pub mean_order: f64,
    pub p90_order: usize,
    pub max_order: usize,
    /// Sum of m^3 over all solves (the dense EVD cost model).
    pub sum_cubed: f64,
    /// Iterations whose reduced matrix grew by exactly one column.
    pub single_column: u64,
    /// Solves done by the bordered (secular) update (`BAANDI_M2_B5=on`).
    pub bordered: u64,
    /// Solves done on rank 0 and broadcast (`BAANDI_M2_EIGEN=root`).
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
        REDUCTION_CALLS[s][bucket(local.len())].fetch_add(1, Ordering::Relaxed);
        REDUCTION_PAYLOAD[s].fetch_add(local.len() as u64, Ordering::Relaxed);
        self.probe()?;
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
/// `Default` is the classic M1/M2 code (all off). Every choice is a function of
/// replicated data only, so all ranks take the same branches.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct DistributedTuning {
    /// `BAANDI_M2_DGEMM=on`: rotations and basis projections through
    /// `multiply_slices_wide_into` (per-panel dgemm at any width on MKL builds;
    /// bit-identical on pure builds).
    pub wide_gemm: bool,
    /// `BAANDI_M2_ORTHO=serial`: orthonormalize in the order of the serial solver
    /// (one basis projection per pass, the DGKS second-pass rule, panel
    /// Gram-Schmidt with the serial panel rule) instead of three projections and
    /// two whole-block column sweeps. Rounding class.
    pub serial_ortho: bool,
    /// DGKS keep ratio of the serial order (`BAANDI_ORTHO`); `None` = two passes.
    pub dgks: Option<f64>,
    /// Panel rule of the serial order: `Some(0)` = ORTHO-BLOCK shape rule on the
    /// widest local row count, `Some(w)` = fixed width w (MGS-B), `None` = the
    /// whole block as one panel (the classic sweep).
    pub panel: Option<usize>,
    /// `BAANDI_M2_LEAN=on`: residual and correction blocks by slice loops
    /// (bit-identical).
    pub lean: bool,
    /// `BAANDI_M2_EIGEN=root`: rank 0 solves the reduced eigenproblem and
    /// broadcasts the lowest `bands` pairs (bit-identical: the same solve on the
    /// same replicated matrix; removes the concurrent replicated solves).
    pub eigen_root: bool,
    /// `BAANDI_M2_B5=on`: bordered (secular) update of the reduced
    /// eigenproblem when exactly one column was appended (the serial solver's B5;
    /// rounding class).
    pub bordered: bool,
    /// `BAANDI_M2_WORKSPACE=on`: run on the persistent per-rank workspace
    /// (`workspace.rs`; storage only, bit-identical).
    pub workspace: bool,
    /// `BAANDI_M2_SELECTED_ROOTS=on`: on Ritz-change (standard) solves, decide
    /// the unconverged roots from the Ritz-value change before any rotation and
    /// form Ritz vectors, H images and residuals for those roots only; all bands
    /// are rotated only at return and at restarts (review F05, the serial lean
    /// solver's rule). Runs on the workspace path. Strict (residual) solves are
    /// unchanged. Rounding class (column subsets of the same GEMM).
    pub selected_roots: bool,
    /// `BAANDI_M2_PANEL_GS`: reduction scheme of the in-panel Gram-Schmidt
    /// (`PanelGs`; default classic, three reductions per column).
    pub panel_gs: PanelGs,
    /// `BAANDI_M2_LEAN_ALGEBRA=on` (MPI-OVERLAP-DAVIDSON): Ritz-change
    /// (standard) solves run the serial lean solver's algebra on the stick
    /// rows: unit-norm corrections appended without Gram-Schmidt, overlap S and
    /// its inverse Cholesky factor kept beside H_r, generalized reduced problem;
    /// one tree sum per iteration carries the S border, the H border and the
    /// selected residual norms (`workspace/lean.rs`,
    /// docs/design/mpi-overlap-davidson.md). Implies the workspace path and
    /// selected-root rotation on those solves; strict solves are unchanged.
    /// Rounding class.
    pub lean_algebra: bool,
    /// `BAANDI_RESTART_KEEP` (ticket BREAKTHROUGH-1-RESTART, R1-b): previous
    /// directions kept at a contraction of the lean algebra path
    /// (`RestartKeep`; default off = the plain restart, bit for bit).
    pub restart_keep: crate::real_davidson::RestartKeep,
}

/// In-panel reduction scheme of the distributed panel Gram-Schmidt
/// (`BAANDI_M2_PANEL_GS`, docs/performance/profile-2026-11-16-mpi.md). The
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
    /// The switches of this process; the serial-order panel and DGKS rules are the
    /// serial solver's (`BAANDI_ORTHO_BLOCK`, `BAANDI_MGS`, `BAANDI_ORTHO`).
    pub fn from_env() -> Self {
        use dft_numerics::switches::{Id, choice, on};
        let panel = if on(Id::ORTHO_BLOCK) {
            Some(0)
        } else {
            match choice(Id::MGS) {
                "blocked" => Some(crate::real_davidson::BLOCK_COLUMNS),
                _ => None,
            }
        };
        Self {
            wide_gemm: on(Id::M2_DGEMM),
            serial_ortho: choice(Id::M2_ORTHO) == "serial",
            dgks: crate::real_davidson::dgks_keep(),
            panel,
            lean: on(Id::M2_LEAN),
            eigen_root: choice(Id::M2_EIGEN) == "root",
            bordered: on(Id::M2_B5),
            workspace: on(Id::M2_WORKSPACE),
            selected_roots: on(Id::M2_SELECTED_ROOTS),
            panel_gs: match choice(Id::M2_PANEL_GS) {
                "pythag" => PanelGs::Pythag,
                "fused" => PanelGs::Fused,
                _ => PanelGs::Classic,
            },
            lean_algebra: on(Id::M2_LEAN_ALGEBRA),
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

thread_local! {
    static START_GRAM_OVERRIDE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}
/// Force `BAANDI_M2_START` for the calling thread (tests with simulated ranks;
/// every rank thread sets the same value); `None` restores the switch.
pub fn set_start_gram_override(gram: Option<bool>) {
    START_GRAM_OVERRIDE.with(|o| o.set(gram));
}
/// `BAANDI_M2_START=gram` (MPI-SCALING-GAP): the workspace solver takes an
/// orthonormal start block (Gram defect at most 1e-12, the restart contraction
/// rule) as its basis after one Gram reduction instead of the panel Gram-Schmidt.
pub fn start_gram() -> bool {
    START_GRAM_OVERRIDE.with(|o| o.get()).unwrap_or_else(|| {
        dft_numerics::switches::choice(dft_numerics::switches::Id::M2_START) == "gram"
    })
}

/// Serial-order orthonormalization (`BAANDI_M2_ORTHO=serial`), the distributed
/// counterpart of the serial solver's projection: normalise (one reduction);
/// per pass project once against `against` (one Gram reduction, local rotation)
/// and run the panel Gram-Schmidt; the second pass is skipped when the DGKS rule
/// holds (every column kept and the smallest kept norm at least `dgks`). Without
/// `against` one pass. Rank decisions use the same relative floor as the serial
/// solver and depend on reduced values only (identical on every rank).
fn orthonormalize_serial<C: Collective>(
    comm: &C,
    ws: &mut RealGemmWorkspace,
    against: Option<&RealMatrix>,
    block: &RealMatrix,
    tuning: &DistributedTuning,
    panel: usize,
) -> Result<RealMatrix, SolverError> {
    let rows = block.rows();
    let mut x = block.clone();
    site(Site::OrthoNormalize);
    let norms = column_norms_squared(comm, &x)?;
    for (j, nn) in norms.iter().enumerate() {
        let len = nn.sqrt();
        for v in x.as_mut_slice()[j * rows..(j + 1) * rows].iter_mut() {
            *v = if len == 0.0 { 0.0 } else { *v / len };
        }
    }
    let against = against.filter(|q| q.cols() > 0);
    let mut skip_second = false;
    for pass in 0..if against.is_some() { 2 } else { 1 } {
        if (pass > 0 && skip_second) || x.cols() == 0 {
            break;
        }
        let count = x.cols();
        if let Some(q) = against {
            let mut overlap = RealMatrix::zeros(q.cols(), count)?;
            let mut scratch = RealMatrix::zeros(rows, count)?;
            site(Site::OrthoProjection);
            project_out_with(
                comm,
                ws,
                q,
                &mut x,
                &mut overlap,
                &mut scratch,
                tuning.wide_gemm,
            )?;
        }
        let (y, min_len) = panel_gram_schmidt(comm, ws, &x, panel, tuning)?;
        if pass == 0
            && against.is_some()
            && let Some(keep) = tuning.dgks
        {
            skip_second = y.cols() == count && min_len >= keep;
        }
        x = y;
    }
    Ok(x)
}

/// Panel Gram-Schmidt of a row-distributed block: panels of `panel` columns; each
/// panel is projected twice against the columns already kept (Gram reduction plus
/// local rotation), then every column of the panel is orthogonalised twice
/// against the panel's earlier kept columns (one reduction of the dot products
/// per pass, as `mgs_panel`) and normalised (one reduction); a column whose
/// remaining norm is not above the dependence tolerance is dropped. With one
/// panel covering the block this is exactly `mgs_panel`. Returns the kept columns
/// and the smallest remaining norm seen (dropped columns included).
fn panel_gram_schmidt<C: Collective>(
    comm: &C,
    ws: &mut RealGemmWorkspace,
    x: &RealMatrix,
    panel: usize,
    tuning: &DistributedTuning,
) -> Result<(RealMatrix, f64), SolverError> {
    let wide = tuning.wide_gemm;
    let (rows, cols) = (x.rows(), x.cols());
    let mut out = RealMatrix::zeros(rows, cols)?;
    let mut kept = 0usize;
    let mut min_len = f64::INFINITY;
    let mut start = 0;
    while start < cols {
        let w = panel.max(1).min(cols - start);
        let mut p = RealMatrix::new(
            rows,
            w,
            x.as_slice()[start * rows..(start + w) * rows].to_vec(),
        )?;
        if kept > 0 {
            let q = RealMatrix::new(rows, kept, out.as_slice()[..kept * rows].to_vec())?;
            let mut overlap = RealMatrix::zeros(kept, w)?;
            let mut scratch = RealMatrix::zeros(rows, w)?;
            site(Site::PanelCross);
            for _ in 0..2 {
                project_out_with(comm, ws, &q, &mut p, &mut overlap, &mut scratch, wide)?;
            }
        }
        let first = kept;
        let mut local = Vec::with_capacity(w + 1);
        let mut sums = Vec::with_capacity(w + 1);
        let mut v = vec![0.0; rows];
        for jj in 0..w {
            v.copy_from_slice(column_of(&p, jj));
            let len = panel_column(
                comm,
                &out.as_slice()[first * rows..kept * rows],
                rows,
                &mut v,
                &mut local,
                &mut sums,
                tuning.panel_gs,
                tuning.dgks,
            )?;
            if len.is_finite() {
                min_len = min_len.min(len);
            }
            if !len.is_finite() || len <= DEPENDENCE_TOLERANCE {
                continue;
            }
            for (o, a) in out.as_mut_slice()[kept * rows..(kept + 1) * rows]
                .iter_mut()
                .zip(&v)
            {
                *o = a / len;
            }
            kept += 1;
        }
        start += w;
    }
    out.reshape_within_storage(rows, kept)?;
    Ok((out, min_len))
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
    v: &mut [f64],
    local: &mut Vec<f64>,
    sums: &mut Vec<f64>,
    mode: PanelGs,
    keep: Option<f64>,
) -> Result<f64, SolverError> {
    let k = q.len().checked_div(rows).unwrap_or(0);
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

/// Sequential modified Gram-Schmidt inside the panel (two passes per column), the
/// distributed counterpart of the serial solver's in-block orthogonalisation:
/// a column is dropped only when its remaining relative norm is at most the
/// dependence tolerance, so near-dependent corrections are kept exactly as in the
/// serial solver (a Gram-matrix Cholesky would discard them). Each reduction is a
/// deterministic tree sum of a short vector.
fn mgs_panel<C: Collective>(comm: &C, x: &RealMatrix) -> Result<RealMatrix, SolverError> {
    let rows = x.rows();
    let mut kept: Vec<f64> = Vec::new(); // column-major, `count` columns of `rows`
    let mut count = 0usize;
    for j in 0..x.cols() {
        let mut v = column_of(x, j).to_vec();
        for _ in 0..2 {
            if count > 0 {
                let mut local = vec![0.0; count];
                for (c, q) in local.iter_mut().zip(kept.chunks(rows)) {
                    *c = q.iter().zip(&v).map(|(a, b)| a * b).sum();
                }
                let mut sums = vec![0.0; count];
                comm.tree_sum_f64(&local, &mut sums).map_err(comm_error)?;
                for (c, q) in sums.iter().zip(kept.chunks(rows)) {
                    for (vi, qi) in v.iter_mut().zip(q) {
                        *vi -= c * qi;
                    }
                }
            }
        }
        let local = [v.iter().map(|a| a * a).sum::<f64>()];
        let mut total = [0.0];
        comm.tree_sum_f64(&local, &mut total).map_err(comm_error)?;
        let len = total[0].sqrt();
        if !len.is_finite() || len <= DEPENDENCE_TOLERANCE {
            continue;
        }
        kept.extend(v.iter().map(|a| a / len));
        count += 1;
    }
    Ok(RealMatrix::new(rows, count, kept)?)
}

fn column_of(x: &RealMatrix, j: usize) -> &[f64] {
    &x.as_slice()[j * x.rows()..(j + 1) * x.rows()]
}

/// Normalise the columns, project twice against `against` (block Gram-Schmidt),
/// drop columns that vanish, then orthonormalise the rest with two CholeskyQR
/// passes. Returns the kept block (possibly with zero columns).
fn orthonormalize<C: Collective>(
    comm: &C,
    ws: &mut RealGemmWorkspace,
    against: Option<&RealMatrix>,
    block: &RealMatrix,
    wide: bool,
) -> Result<RealMatrix, SolverError> {
    let rows = block.rows();
    let mut x = block.clone();
    let norms = column_norms_squared(comm, &x)?;
    for (j, nn) in norms.iter().enumerate() {
        let len = nn.sqrt();
        let col = &mut x.as_mut_slice()[j * rows..(j + 1) * rows];
        for v in col.iter_mut() {
            *v = if len == 0.0 { 0.0 } else { *v / len };
        }
    }
    if let Some(q) = against.filter(|q| q.cols() > 0) {
        for _ in 0..2 {
            let mut overlap = RealMatrix::zeros(q.cols(), x.cols())?;
            let mut scratch = RealMatrix::zeros(rows, x.cols())?;
            project_out_with(comm, ws, q, &mut x, &mut overlap, &mut scratch, wide)?;
        }
    }
    // Drop columns that vanished after the projection, renormalise the others.
    let norms = column_norms_squared(comm, &x)?;
    let survivors: Vec<usize> = (0..x.cols())
        .filter(|&j| norms[j].sqrt() > DEPENDENCE_TOLERANCE)
        .collect();
    if survivors.is_empty() {
        return Ok(RealMatrix::zeros(rows, 0)?);
    }
    let mut y = RealMatrix::zeros(rows, survivors.len())?;
    for (c, &j) in survivors.iter().enumerate() {
        let len = norms[j].sqrt();
        for i in 0..rows {
            y[(i, c)] = x[(i, j)] / len;
        }
    }
    let mut y = mgs_panel(comm, &y)?;
    // Renormalising a nearly dependent column amplifies rounding error along the
    // basis; one more projection plus a panel pass (nothing drops now) restores
    // orthogonality to working precision.
    if let Some(q) = against.filter(|q| q.cols() > 0 && y.cols() > 0) {
        let mut overlap = RealMatrix::zeros(q.cols(), y.cols())?;
        let mut scratch = RealMatrix::zeros(rows, y.cols())?;
        project_out_with(comm, ws, q, &mut y, &mut overlap, &mut scratch, wide)?;
        y = mgs_panel(comm, &y)?;
    }
    Ok(y)
}

/// The orthonormalization selected by `tuning` (classic or serial order).
fn orthonormalize_tuned<C: Collective>(
    comm: &C,
    ws: &mut RealGemmWorkspace,
    against: Option<&RealMatrix>,
    block: &RealMatrix,
    tuning: &DistributedTuning,
    panel: usize,
) -> Result<RealMatrix, SolverError> {
    if tuning.serial_ortho {
        orthonormalize_serial(comm, ws, against, block, tuning, panel)
    } else {
        orthonormalize(comm, ws, against, block, tuning.wide_gemm)
    }
}

/// One reduced eigenproblem: the bordered (secular) update when exactly one
/// column was appended to an unbroken base (`BAANDI_M2_B5=on`, same rule and
/// streak bound as the serial solver), otherwise the dense solve. With B5 the
/// result becomes the next base. All pairs are written.
fn solve_reduced(
    eig: &mut RealEighWorkspace,
    b5: Option<(&mut Bordered, &mut RealGemmWorkspace)>,
    reduced: &RealMatrix,
    prefix: usize,
    values: &mut [f64],
    coefficients: &mut RealMatrix,
) -> Result<(), SolverError> {
    let m = reduced.rows();
    let Some((b, ws)) = b5 else {
        eig.solve_into(reduced, values, coefficients)?;
        return Ok(());
    };
    if prefix == 0 {
        b.invalidate();
    }
    let bordered = prefix + 1 == m
        && b.order == prefix
        && b.streak < MAX_STREAK
        && b.update(ws, reduced, values, coefficients, None)?;
    if bordered {
        b.streak += 1;
        SOLVE_COUNTS[2].fetch_add(1, Ordering::Relaxed);
    } else {
        eig.solve_into(reduced, values, coefficients)?;
        b.streak = 0;
    }
    b.store(values, coefficients);
    Ok(())
}

fn append_columns(a: &mut RealMatrix, b: &RealMatrix) -> Result<(), SolverError> {
    let mut data = a.as_slice().to_vec();
    data.extend_from_slice(b.as_slice());
    *a = RealMatrix::new(a.rows(), a.cols() + b.cols(), data)?;
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
        if tuning.workspace || tuning.selected_roots || tuning.lean_algebra {
            let mut apply = |x: &RealMatrix, y: &mut RealMatrix| -> Result<(), SolverError> {
                site(Site::Apply);
                let r = apply_distributed(comm, &layout, op, x)?;
                y.as_mut_slice().copy_from_slice(r.as_slice());
                Ok(())
            };
            return workspace::solve_core_ws(
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
            );
        }
        let pool = match op.worker_pool() {
            Some(p) => p,
            None => WorkerPool::new(op.inner_threads())
                .map_err(|e| SolverError::Application(e.to_string()))?,
        };
        let mut apply = |x: &RealMatrix| {
            site(Site::Apply);
            apply_distributed(comm, &layout, op, x)
        };
        self.solve_core(
            comm,
            &layout,
            pool,
            &mut apply,
            initial,
            diag_local,
            Some(diagonal),
            policy,
            ritz,
            certify,
            capacity_per_band,
            tuning,
        )
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
    /// H images). Needs the workspace path (`tuning.workspace` or
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
        if diagonal_local.len() != mine.len() || mine.iter().any(|&r| r >= n) {
            return Err(SolverError::InvalidParameter("bands or diagonal length"));
        }
        let ranks = comm.size();
        // BREAKTHROUGH-4 (`BAANDI_M2_LAYOUT_CACHE=on`): the row lists of all
        // ranks are reused from this thread's previous solve when every rank
        // holds the same local rows (one 1-double reduction decides, so the
        // decision is replicated); otherwise exchanged as before.
        let cache = layout_cache();
        let reuse = if cache {
            let hit = LAYOUT_CACHE.with(|l| {
                l.borrow().as_ref().is_some_and(|(rows, lay)| {
                    lay.n == n
                        && lay.rank == comm.rank()
                        && lay.counts.len() == ranks
                        && rows.as_slice() == mine
                })
            });
            let mut flag = [0.0];
            comm.tree_sum_f64(&[f64::from(u8::from(!hit))], &mut flag)
                .map_err(comm_error)?;
            flag[0] == 0.0
        } else {
            false
        };
        let layout: std::sync::Arc<Layout> = if reuse {
            LAYOUT_COUNTS[0].fetch_add(1, Ordering::Relaxed);
            LAYOUT_CACHE.with(|l| l.borrow().as_ref().expect("checked above").1.clone())
        } else {
            // Setup (once per solve): every rank learns the row lists of all ranks.
            let mut counts_f = vec![0.0; ranks];
            comm.allgatherv_f64(&[mine.len() as f64], &vec![1; ranks], &mut counts_f)
                .map_err(comm_error)?;
            let counts: Vec<usize> = counts_f.iter().map(|&c| c as usize).collect();
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
            if cache {
                LAYOUT_COUNTS[1].fetch_add(1, Ordering::Relaxed);
                LAYOUT_CACHE.with(|l| *l.borrow_mut() = Some((mine.to_vec(), layout.clone())));
            }
            layout
        };
        let nl = layout.nl();
        stage_add(STAGE_SETUP, t_setup);
        if tuning.workspace || tuning.selected_roots || tuning.lean_algebra {
            let supplied = op.worker_pool();
            let threads = op.inner_threads();
            let mut apply = |x: &RealMatrix, y: &mut RealMatrix| -> Result<(), SolverError> {
                site(Site::Apply);
                let start = dft_numerics::stage_clock::start();
                y.as_mut_slice().fill(0.0);
                op.apply_block(comm, x, y)?;
                stage_add(5, start);
                stage_add(6, start);
                let scan = dft_numerics::stage_clock::start();
                let invalid = y.as_slice().iter().any(|v| !v.is_finite());
                stage_add(STAGE_APPEND, scan);
                if invalid {
                    return Err(SolverError::Application(
                        "operator returned an invalid block".into(),
                    ));
                }
                Ok(())
            };
            let result = workspace::solve_core_ws(
                self,
                comm,
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
            return result;
        }
        // The allocating path has no kept start: drop it, then give the local
        // rows of the result below.
        if let Some((_, k)) = kept.as_mut() {
            k.invalidate();
        }
        // It also needs a replicated start: gather a local one (MPI-LAZY-ORBITALS
        // with the workspace off), and it always returns the gathered block.
        let gathered;
        let initial = match initial {
            Some(DistributedStart::Global(x)) => Some(x),
            Some(DistributedStart::Local(x)) => {
                if x.rows() != nl {
                    return Err(SolverError::InvalidParameter(
                        "local start rows differ from the operator's rows",
                    ));
                }
                gathered = gather_full(comm, &layout, x)?;
                Some(&gathered)
            }
            None => None,
        };
        let pool = match op.worker_pool() {
            Some(p) => p,
            None => WorkerPool::new(op.inner_threads())
                .map_err(|e| SolverError::Application(e.to_string()))?,
        };
        let mut apply = |x: &RealMatrix| -> Result<RealMatrix, SolverError> {
            site(Site::Apply);
            let start = dft_numerics::stage_clock::start();
            let mut y = RealMatrix::zeros(nl, x.cols())?;
            op.apply_block(comm, x, &mut y)?;
            stage_add(5, start);
            stage_add(6, start);
            if y.as_slice().iter().any(|v| !v.is_finite()) {
                return Err(SolverError::Application(
                    "operator returned an invalid block".into(),
                ));
            }
            Ok(y)
        };
        let result = self.solve_core(
            comm,
            &layout,
            pool,
            &mut apply,
            initial,
            diagonal_local,
            None,
            policy,
            ritz,
            certify,
            capacity_per_band,
            tuning,
        )?;
        if let Some((x, _)) = kept {
            *x = local_rows(&layout, &result.eigenvectors)?;
        }
        stage_add(STAGE_TOTAL, t_total);
        Ok(result)
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_core<C: Collective>(
        &self,
        comm: &C,
        layout: &Layout,
        pool: std::sync::Arc<WorkerPool>,
        apply: &mut dyn FnMut(&RealMatrix) -> Result<RealMatrix, SolverError>,
        initial: Option<&RealMatrix>,
        diag_local: &[f64],
        full_diagonal: Option<&[f64]>,
        policy: BandTolerances,
        ritz: bool,
        certify: bool,
        capacity_per_band: usize,
        tuning: DistributedTuning,
    ) -> Result<RealEigenSolution, SolverError> {
        let n = layout.n;
        if !(2..=12).contains(&capacity_per_band) {
            return Err(SolverError::InvalidParameter(
                "davidson_capacity_per_band must be 2 to 12",
            ));
        }
        if self.bands == 0 || self.bands > n {
            return Err(SolverError::InvalidParameter("bands or diagonal length"));
        }
        if !self.tolerance.is_finite()
            || self.tolerance <= 0.0
            || self.max_iterations == 0
            || policy.strict_bands == 0
            || policy.strict_bands > self.bands
            || !policy.buffer_tolerance.is_finite()
            || policy.buffer_tolerance < self.tolerance
        {
            return Err(SolverError::InvalidParameter(
                "invalid tolerance, iteration limit or strict prefix",
            ));
        }
        let bands = self.bands;
        let capacity = self.max_subspace.max(bands * capacity_per_band).min(n);
        let nl = layout.nl();
        let max_rows = layout.counts.iter().copied().max().unwrap_or(1);
        let b5_pool = pool.clone();
        let mut ws = RealGemmWorkspace::with_pool(max_rows.max(1), capacity.max(1), pool)?;
        // Reduced eigenproblem: replicated on every rank, or on rank 0 with a
        // broadcast of the lowest `bands` pairs (status, values, coefficients).
        let root_mode = tuning.eigen_root && comm.size() > 1;
        let solves_here = !root_mode || comm.rank() == 0;
        let mut eig = if solves_here {
            Some(RealEighWorkspace::new(capacity)?)
        } else {
            None
        };
        // B5 state and its own GEMM workspace (the coefficient product is
        // capacity x capacity, independent of the local row count).
        let mut b5 = if tuning.bordered && solves_here {
            Some((
                Bordered::new(capacity),
                RealGemmWorkspace::with_pool(capacity.max(1), capacity.max(1), b5_pool)?,
            ))
        } else {
            None
        };
        let mut eigen_buffer = if root_mode {
            vec![0.0; 1 + bands + capacity * bands]
        } else {
            Vec::new()
        };
        let mut rng = Rng(self.seed);
        let panel = tuning.panel_width(max_rows.max(1));
        let wide = tuning.wide_gemm;

        // Initial basis.
        let start = match initial {
            Some(x) => x.clone(),
            None => {
                let full: Vec<f64> = match full_diagonal {
                    Some(d) => d.to_vec(),
                    None => {
                        let mut flat = vec![0.0; n];
                        comm.allgatherv_f64(diag_local, &layout.counts, &mut flat)
                            .map_err(comm_error)?;
                        let mut d = vec![0.0; n];
                        let mut at = 0;
                        for (r, &c) in layout.counts.iter().enumerate() {
                            for (&g, &v) in layout.index[r].iter().zip(&flat[at..at + c]) {
                                d[g] = v;
                            }
                            at += c;
                        }
                        d
                    }
                };
                free_electron_initial_guess(&full, bands)?
            }
        };
        let mut basis = orthonormalize_tuned(
            comm,
            &mut ws,
            None,
            &local_rows(layout, &start)?,
            &tuning,
            panel,
        )?;
        while basis.cols() < bands {
            let w = bands - basis.cols();
            let random = RealMatrix::new(n, w, (0..n * w).map(|_| rng.next()).collect())?;
            let extra = orthonormalize_tuned(
                comm,
                &mut ws,
                Some(&basis),
                &local_rows(layout, &random)?,
                &tuning,
                panel,
            )?;
            if extra.cols() == 0 {
                return Err(SolverError::InvalidParameter(
                    "unable to complete initial basis",
                ));
            }
            append_columns(&mut basis, &extra)?;
        }
        let mut applications = 0usize;
        let mut products = apply(&basis)?;
        applications += basis.cols();

        let threshold = |j: usize| {
            if j < policy.strict_bands - usize::from(ritz) {
                self.tolerance
            } else {
                policy.buffer_tolerance
            }
        };
        let meets_policy = |r: &[f64]| r.iter().enumerate().all(|(j, &x)| x <= threshold(j));
        let meets_strict = |r: &[f64]| {
            r[..policy.strict_bands]
                .iter()
                .all(|&x| x <= self.tolerance)
        };
        let mut history = Vec::with_capacity(self.max_iterations);
        let (mut restarts, mut contractions, mut certificates) = (0usize, 0usize, 0usize);
        let mut previous = vec![f64::NAN; bands];
        let mut reduced_width = 0usize;
        let mut reduced_prev = RealMatrix::zeros(1, 1)?;
        let mut max_basis = 0usize;
        let buffer: Vec<bool> = (0..bands).map(|j| j >= policy.strict_bands).collect();

        for iteration in 1..=self.max_iterations {
            workspace::call_probe(iteration);
            let m = basis.cols();
            max_basis = max_basis.max(m);
            // Reduced matrix V^T HV (replicated), symmetrised from the upper triangle.
            // Incremental, as in the serial solver: only the columns added since the last
            // iteration are formed (m x new products); the rest is reused.
            let t0 = dft_numerics::stage_clock::start();
            let mut reduced = RealMatrix::zeros(m, m)?;
            for c in 0..reduced_width {
                for r in 0..reduced_width {
                    reduced[(r, c)] = reduced_prev[(r, c)];
                }
            }
            let fresh = m - reduced_width;
            let prefix = reduced_width;
            let mut fresh_products = RealMatrix::zeros(nl, fresh)?;
            fresh_products
                .as_mut_slice()
                .copy_from_slice(&products.as_slice()[reduced_width * nl..m * nl]);
            let mut block = RealMatrix::zeros(m, fresh)?;
            site(Site::ReducedGram);
            gram(comm, &mut ws, &basis, &fresh_products, &mut block)?;
            for c in reduced_width..m {
                for r in 0..=c {
                    let v = block[(r, c - reduced_width)];
                    reduced[(r, c)] = v;
                    reduced[(c, r)] = v;
                }
            }
            reduced_width = m;
            stage_add(0, t0);
            reduced_prev = reduced.clone();
            let mut values = vec![0.0; m];
            let mut coefficients = RealMatrix::zeros(m, m)?;
            let t0 = dft_numerics::stage_clock::start();
            record_solve(m, fresh);
            let mut status: Result<(), SolverError> = Ok(());
            if let Some(eig) = eig.as_mut() {
                status = solve_reduced(
                    eig,
                    b5.as_mut().map(|(b, g)| (b, g)),
                    &reduced,
                    prefix,
                    &mut values,
                    &mut coefficients,
                );
            }
            if root_mode {
                let buffer = &mut eigen_buffer[..1 + bands + m * bands];
                if comm.rank() == 0 {
                    buffer[0] = if status.is_ok() { 0.0 } else { 1.0 };
                    buffer[1..1 + bands].copy_from_slice(&values[..bands]);
                    buffer[1 + bands..].copy_from_slice(&coefficients.as_slice()[..m * bands]);
                }
                comm.bcast_f64(buffer, 0).map_err(comm_error)?;
                if comm.rank() != 0 {
                    if buffer[0] != 0.0 {
                        status = Err(SolverError::Application(
                            "reduced eigenproblem failed on rank 0".into(),
                        ));
                    } else {
                        values[..bands].copy_from_slice(&buffer[1..1 + bands]);
                        coefficients.as_mut_slice()[..m * bands]
                            .copy_from_slice(&buffer[1 + bands..]);
                    }
                }
                SOLVE_COUNTS[3].fetch_add(1, Ordering::Relaxed);
            }
            status?;
            stage_add(1, t0);
            let t0 = dft_numerics::stage_clock::start();
            let lowest = RealMatrix::new(m, bands, coefficients.as_slice()[..m * bands].to_vec())?;
            let mut vectors = RealMatrix::zeros(nl, bands)?;
            let mut rotated = RealMatrix::zeros(nl, bands)?;
            rotate_local_with(&mut ws, &basis, &lowest, &mut vectors, wide)?;
            rotate_local_with(&mut ws, &products, &lowest, &mut rotated, wide)?;
            rotation_add(nl, m, bands);
            stage_add(2, t0);
            let residual_of = |hv: &RealMatrix,
                               v: &RealMatrix,
                               comm: &C|
             -> Result<(RealMatrix, Vec<f64>), SolverError> {
                let mut r = RealMatrix::zeros(nl, bands)?;
                if tuning.lean {
                    // M2_LEAN: the same element expression by column slices.
                    for (j, &lambda) in values.iter().enumerate().take(bands) {
                        let (h, x) = (column_of(hv, j), column_of(v, j));
                        for ((r, &h), &x) in r.as_mut_slice()[j * nl..(j + 1) * nl]
                            .iter_mut()
                            .zip(h)
                            .zip(x)
                        {
                            *r = h - lambda * x;
                        }
                    }
                } else {
                    for j in 0..bands {
                        for i in 0..nl {
                            r[(i, j)] = hv[(i, j)] - values[j] * v[(i, j)];
                        }
                    }
                }
                site(Site::ResidualNorms);
                let norms = column_norms_squared(comm, &r)?
                    .into_iter()
                    .map(f64::sqrt)
                    .collect();
                Ok((r, norms))
            };
            let t0 = dft_numerics::stage_clock::start();
            let (mut residual_block, mut residuals) = residual_of(&rotated, &vectors, comm)?;
            stage_add(3, t0);
            let mut certified = false;
            let mut trusted = true;
            if !ritz
                && (certify || iteration == self.max_iterations)
                && (meets_policy(&residuals)
                    || (iteration == self.max_iterations && meets_strict(&residuals)))
            {
                certificates += 1;
                let actual = apply(&vectors)?;
                applications += bands;
                (residual_block, residuals) = residual_of(&actual, &vectors, comm)?;
                certified = true;
                trusted = meets_policy(&residuals);
            }
            let maximum = residuals.iter().copied().fold(0.0, f64::max);
            history.push(maximum);
            let converged: Vec<bool> = (0..bands)
                .map(|j| {
                    if ritz {
                        iteration > 1 && (values[j] - previous[j]).abs() < threshold(j)
                    } else {
                        residuals[j] <= threshold(j)
                    }
                })
                .collect();
            previous.copy_from_slice(&values[..bands]);
            let finished = converged.iter().all(|&b| b)
                || (iteration == self.max_iterations
                    && converged[..policy.strict_bands].iter().all(|&b| b));
            if finished {
                let full = gather_full(comm, layout, &vectors)?;
                return Ok(RealEigenSolution {
                    eigenvalues: values[..bands].to_vec(),
                    eigenvectors: full,
                    residuals,
                    converged,
                    buffer,
                    residuals_certified: certified,
                    iterations: iteration,
                    history,
                    operator_applications: applications,
                    restarts,
                    contracted_restarts: contractions,
                    certification_attempts: certificates,
                    max_basis_dimension: max_basis,
                });
            }
            if iteration == self.max_iterations {
                return Err(SolverError::NotConverged {
                    iterations: iteration,
                    max_residual: maximum,
                    residuals,
                });
            }
            if ritz && iteration == 1 && basis.cols() == n {
                continue;
            }
            // Corrections for the unconverged roots (row-local preconditioning).
            let selected: Vec<usize> = (0..bands)
                .filter(|&j| !converged[j])
                .take(bands.min(capacity - bands))
                .collect();
            let mut corrections = RealMatrix::zeros(nl, selected.len())?;
            let precondition = |d: f64, lambda: f64, r: f64| {
                let delta = d - lambda;
                let denominator = if delta < 0.0 {
                    -delta.abs().max(0.1)
                } else {
                    delta.max(0.1)
                };
                -r / denominator
            };
            for (k, &j) in selected.iter().enumerate() {
                if tuning.lean {
                    let r = column_of(&residual_block, j);
                    for ((c, &d), &r) in corrections.as_mut_slice()[k * nl..(k + 1) * nl]
                        .iter_mut()
                        .zip(diag_local)
                        .zip(r)
                    {
                        *c = precondition(d, values[j], r);
                    }
                } else {
                    for i in 0..nl {
                        corrections[(i, k)] =
                            precondition(diag_local[i], values[j], residual_block[(i, j)]);
                    }
                }
            }
            if !trusted || basis.cols() + corrections.cols() > capacity {
                restarts += 1;
                reduced_width = 0;
                let mut defect = RealMatrix::zeros(bands, bands)?;
                site(Site::RestartDefect);
                gram(comm, &mut ws, &vectors, &vectors, &mut defect)?;
                let mut worst = 0.0_f64;
                for j in 0..bands {
                    for i in 0..bands {
                        worst = worst.max((defect[(i, j)] - if i == j { 1.0 } else { 0.0 }).abs());
                    }
                }
                if trusted && worst <= 1e-12 {
                    contractions += 1;
                    basis = vectors.clone();
                    products = rotated.clone();
                } else {
                    basis = orthonormalize_tuned(comm, &mut ws, None, &vectors, &tuning, panel)?;
                    products = apply(&basis)?;
                    applications += basis.cols();
                }
            }
            let t0 = dft_numerics::stage_clock::start();
            let mut extra =
                orthonormalize_tuned(comm, &mut ws, Some(&basis), &corrections, &tuning, panel)?;
            stage_add(4, t0);
            if extra.cols() == 0 {
                let w = bands.min(capacity - basis.cols());
                let random = RealMatrix::new(n, w, (0..n * w).map(|_| rng.next()).collect())?;
                extra = orthonormalize_tuned(
                    comm,
                    &mut ws,
                    Some(&basis),
                    &local_rows(layout, &random)?,
                    &tuning,
                    panel,
                )?;
            }
            if extra.cols() == 0 {
                return Err(SolverError::NotConverged {
                    iterations: iteration,
                    max_residual: maximum,
                    residuals,
                });
            }
            let extra_products = apply(&extra)?;
            applications += extra.cols();
            append_columns(&mut products, &extra_products)?;
            append_columns(&mut basis, &extra)?;
        }
        unreachable!("validated positive iteration limit")
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
