//! F06 (`BAANDI_M2_WORKSPACE=on`): the distributed Davidson on a persistent
//! per-rank workspace. Same arithmetic, in the same order, as `solve_core`
//! (bit-identical at a fixed rank count); what changes is storage only:
//! - basis and H products live in buffers of `capacity` columns, new columns are
//!   written in place (no copy of the whole basis per iteration);
//! - the reduced matrix is kept with the fixed leading dimension `capacity` and
//!   only its new border is written (no per-iteration zero matrix and no clone of
//!   the previous one); the dense solver reads a packed m x m copy;
//! - the fresh products, the lowest coefficient columns and the Gram operands
//!   are passed as views (slices), not copies;
//! - Ritz vectors, residuals, corrections, orthonormalization scratch, the Gram
//!   partial sums, the GEMM, EVD and B5 workspaces and the root broadcast buffer
//!   are reused across iterations and across solves (one workspace per thread,
//!   that is per MPI rank or per simulated rank).
//!
//! Buffers grow only; every growth is counted (`workspace_stats`). The classic
//! (non-serial) orthonormalization order still allocates inside (it is not the
//! default, `BAANDI_M2_ORTHO=serial` is), and so does the M1 apply.
use super::*;
use crate::distributed::{
    column_norms_squared_into, gram_slices, project_out_slices, rotate_slices_with,
};
use std::cell::{Cell, RefCell};
use std::sync::Arc;

/// Growth and reuse counters of the persistent workspace (process-wide).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct WorkspaceStats {
    /// Solves that ran on the workspace.
    pub solves: u64,
    /// Solves during which no buffer grew and no inner workspace was rebuilt.
    pub warm_solves: u64,
    /// Buffer growth events (each one allocation).
    pub growths: u64,
    /// Bytes requested by those growths (f64 entries x 8, or the inner
    /// workspace estimate for a GEMM/EVD/B5 rebuild).
    pub growth_bytes: u64,
}

// solves, warm solves, growths, growth bytes
static WS_COUNTS: [AtomicU64; 4] = [const { AtomicU64::new(0) }; 4];

/// The workspace counters accumulated so far by this process.
pub fn workspace_stats() -> WorkspaceStats {
    let load = |i: usize| WS_COUNTS[i].load(Ordering::Relaxed);
    WorkspaceStats {
        solves: load(0),
        warm_solves: load(1),
        growths: load(2),
        growth_bytes: load(3),
    }
}

thread_local! {
    static WORKSPACE: RefCell<Option<Box<DistributedDavidsonWorkspace>>> =
        const { RefCell::new(None) };
    // Growth events of the solve running on this thread (warm-solve detection).
    static SOLVE_GROWTHS: Cell<u64> = const { Cell::new(0) };
    static PROBE: Cell<Option<fn(usize)>> = const { Cell::new(None) };
}

/// Drop the calling thread's (rank's) persistent workspace.
pub fn release_workspace() {
    WORKSPACE.with(|w| w.borrow_mut().take());
}

/// Test hook: `probe(iteration)` is called at the start of every Davidson
/// iteration of a workspace solve on the calling thread (allocation tests).
#[doc(hidden)]
pub fn set_iteration_probe(probe: Option<fn(usize)>) {
    PROBE.with(|p| p.set(probe));
}

/// Call the iteration probe of the calling thread, if any (both solver paths).
pub(super) fn call_probe(iteration: usize) {
    if let Some(probe) = PROBE.with(|p| p.get()) {
        probe(iteration);
    }
}

fn grown(bytes: usize) {
    WS_COUNTS[2].fetch_add(1, Ordering::Relaxed);
    WS_COUNTS[3].fetch_add(bytes as u64, Ordering::Relaxed);
    SOLVE_GROWTHS.with(|g| g.set(g.get() + 1));
}

/// Give `m` the active shape rows x cols, growing its storage if needed.
fn fit(m: &mut RealMatrix, rows: usize, cols: usize) -> Result<(), SolverError> {
    if m.reshape_within_storage(rows, cols).is_err() {
        let entries = rows
            .checked_mul(cols)
            .ok_or(SolverError::InvalidParameter("workspace size"))?;
        grown(entries * 8);
        m.reserve_storage(entries);
        m.reshape_within_storage(rows, cols)?;
    }
    Ok(())
}

/// Make sure `m` can hold `entries` values without growing later.
fn reserve(m: &mut RealMatrix, entries: usize) -> Result<(), SolverError> {
    let (rows, cols) = (m.rows(), m.cols());
    fit(m, entries, 1)?;
    m.reshape_within_storage(rows, cols)?;
    Ok(())
}

fn fit_vec(v: &mut Vec<f64>, len: usize) {
    if v.len() < len {
        grown((len - v.len()) * 8);
        v.resize(len, 0.0);
    }
}

fn empty() -> RealMatrix {
    RealMatrix::new(0, 0, Vec::new()).expect("empty matrix")
}

/// Scratch of the serial-order orthonormalization.
struct OrthoScratch {
    x: RealMatrix,
    y: RealMatrix,
    p: RealMatrix,
    overlap: RealMatrix,
    scratch: RealMatrix,
    partial: RealMatrix,
    v: Vec<f64>,
    local: Vec<f64>,
    sums: Vec<f64>,
    norms: Vec<f64>,
    norm_partial: Vec<f64>,
}
impl OrthoScratch {
    fn new() -> Self {
        Self {
            x: empty(),
            y: empty(),
            p: empty(),
            overlap: empty(),
            scratch: empty(),
            partial: empty(),
            v: Vec::new(),
            local: Vec::new(),
            sums: Vec::new(),
            norms: Vec::new(),
            norm_partial: Vec::new(),
        }
    }
}

/// Persistent buffers of the distributed Davidson for one rank.
pub struct DistributedDavidsonWorkspace {
    /// (max local rows, capacity, pool address) of `gemm`.
    gemm: Option<((usize, usize, usize), RealGemmWorkspace)>,
    /// Pool used when the operator supplies none, with its thread count.
    own_pool: Option<(usize, Arc<WorkerPool>)>,
    eig: Option<(usize, RealEighWorkspace)>,
    b5: Option<((usize, usize), Bordered, RealGemmWorkspace)>,
    eigen_buffer: Vec<f64>,
    basis: RealMatrix,
    products: RealMatrix,
    /// Reduced matrix, leading dimension `held_ld`.
    held: Vec<f64>,
    reduced: RealMatrix,
    coefficients: RealMatrix,
    values: Vec<f64>,
    block: RealMatrix,
    partial: RealMatrix,
    vectors: RealMatrix,
    rotated: RealMatrix,
    residual: RealMatrix,
    actual: RealMatrix,
    corrections: RealMatrix,
    extra: RealMatrix,
    extra_products: RealMatrix,
    defect: RealMatrix,
    norms: Vec<f64>,
    norm_partial: Vec<f64>,
    converged: Vec<bool>,
    selected: Vec<usize>,
    ortho: OrthoScratch,
}

impl DistributedDavidsonWorkspace {
    fn new() -> Self {
        Self {
            gemm: None,
            own_pool: None,
            eig: None,
            b5: None,
            eigen_buffer: Vec::new(),
            basis: empty(),
            products: empty(),
            held: Vec::new(),
            reduced: empty(),
            coefficients: empty(),
            values: Vec::new(),
            block: empty(),
            partial: empty(),
            vectors: empty(),
            rotated: empty(),
            residual: empty(),
            actual: empty(),
            corrections: empty(),
            extra: empty(),
            extra_products: empty(),
            defect: empty(),
            norms: Vec::new(),
            norm_partial: Vec::new(),
            converged: Vec::new(),
            selected: Vec::new(),
            ortho: OrthoScratch::new(),
        }
    }
}

/// The operator's pool, or a pool of `inner_threads` kept in the workspace.
fn pool_of(
    ws: &mut DistributedDavidsonWorkspace,
    supplied: Option<Arc<WorkerPool>>,
    inner_threads: usize,
) -> Result<Arc<WorkerPool>, SolverError> {
    if let Some(p) = supplied {
        return Ok(p);
    }
    if let Some((t, p)) = &ws.own_pool
        && *t == inner_threads
    {
        return Ok(p.clone());
    }
    let p = WorkerPool::new(inner_threads).map_err(|e| SolverError::Application(e.to_string()))?;
    ws.own_pool = Some((inner_threads, p.clone()));
    Ok(p)
}

/// Residuals hv - lambda v of the first `cols` columns (both element orders of
/// `solve_core`) and their norms (one reduction).
#[allow(clippy::too_many_arguments)]
fn residual_into<C: Collective>(
    comm: &C,
    hv: &RealMatrix,
    v: &RealMatrix,
    lambdas: &[f64],
    lean: bool,
    out: &mut RealMatrix,
    norms: &mut Vec<f64>,
    partial: &mut Vec<f64>,
) -> Result<(), SolverError> {
    let (nl, cols) = (v.rows(), lambdas.len());
    fit(out, nl, cols)?;
    if lean {
        for (j, &lambda) in lambdas.iter().enumerate() {
            let (h, x) = (column_of(hv, j), column_of(v, j));
            for ((r, &h), &x) in out.as_mut_slice()[j * nl..(j + 1) * nl]
                .iter_mut()
                .zip(h)
                .zip(x)
            {
                *r = h - lambda * x;
            }
        }
    } else {
        for (j, &lambda) in lambdas.iter().enumerate() {
            for i in 0..nl {
                out[(i, j)] = hv[(i, j)] - lambda * v[(i, j)];
            }
        }
    }
    column_norms_squared_into(comm, out, partial, norms)?;
    for x in norms.iter_mut() {
        *x = x.sqrt();
    }
    Ok(())
}

/// Panel Gram-Schmidt of `o.x` into `o.y` (the arithmetic of
/// `panel_gram_schmidt`, on the scratch); returns the smallest remaining norm.
fn panel_gs_into<C: Collective>(
    comm: &C,
    gemm: &mut RealGemmWorkspace,
    panel: usize,
    wide: bool,
    o: &mut OrthoScratch,
) -> Result<f64, SolverError> {
    let OrthoScratch {
        x,
        y,
        p,
        overlap,
        scratch,
        partial,
        v,
        local,
        sums,
        ..
    } = o;
    let (rows, cols) = (x.rows(), x.cols());
    fit(y, rows, cols)?;
    let mut kept = 0usize;
    let mut min_len = f64::INFINITY;
    let mut start = 0;
    if v.len() < rows {
        grown((rows - v.len()) * 8);
    }
    v.clear();
    v.resize(rows, 0.0);
    while start < cols {
        let w = panel.max(1).min(cols - start);
        fit(p, rows, w)?;
        p.as_mut_slice()
            .copy_from_slice(&x.as_slice()[start * rows..(start + w) * rows]);
        if kept > 0 {
            fit(overlap, kept, w)?;
            fit(scratch, rows, w)?;
            fit(partial, kept, w)?;
            for _ in 0..2 {
                project_out_slices(
                    comm,
                    gemm,
                    &y.as_slice()[..kept * rows],
                    kept,
                    p,
                    overlap,
                    scratch,
                    partial,
                    wide,
                )?;
            }
        }
        let first = kept;
        for jj in 0..w {
            v.copy_from_slice(column_of(p, jj));
            for _ in 0..2 {
                let k = kept - first;
                if k > 0 {
                    local.clear();
                    for c in first..kept {
                        let q = &y.as_slice()[c * rows..(c + 1) * rows];
                        local.push(q.iter().zip(v.iter()).map(|(a, b)| a * b).sum::<f64>());
                    }
                    sums.clear();
                    sums.resize(k, 0.0);
                    comm.tree_sum_f64(local, sums).map_err(comm_error)?;
                    for (c, s) in (first..kept).zip(sums.iter()) {
                        let q = &y.as_slice()[c * rows..(c + 1) * rows];
                        for (vi, qi) in v.iter_mut().zip(q) {
                            *vi -= s * qi;
                        }
                    }
                }
            }
            let local_norm = [v.iter().map(|a| a * a).sum::<f64>()];
            let mut total = [0.0];
            comm.tree_sum_f64(&local_norm, &mut total)
                .map_err(comm_error)?;
            let len = total[0].sqrt();
            if len.is_finite() {
                min_len = min_len.min(len);
            }
            if !len.is_finite() || len <= DEPENDENCE_TOLERANCE {
                continue;
            }
            for (o, a) in y.as_mut_slice()[kept * rows..(kept + 1) * rows]
                .iter_mut()
                .zip(v.iter())
            {
                *o = a / len;
            }
            kept += 1;
        }
        start += w;
    }
    y.reshape_within_storage(rows, kept)?;
    Ok(min_len)
}

/// `orthonormalize_tuned` into `out`. The serial order runs on the scratch
/// (same arithmetic as `orthonormalize_serial`); the classic order calls the
/// allocating classic routine and copies its result.
#[allow(clippy::too_many_arguments)]
fn ortho_into<C: Collective>(
    comm: &C,
    gemm: &mut RealGemmWorkspace,
    against: Option<&RealMatrix>,
    block: &RealMatrix,
    tuning: &DistributedTuning,
    panel: usize,
    o: &mut OrthoScratch,
    out: &mut RealMatrix,
) -> Result<(), SolverError> {
    if !tuning.serial_ortho {
        let r = orthonormalize(comm, gemm, against, block, tuning.wide_gemm)?;
        fit(out, r.rows(), r.cols())?;
        out.as_mut_slice().copy_from_slice(r.as_slice());
        return Ok(());
    }
    let rows = block.rows();
    fit(&mut o.x, rows, block.cols())?;
    o.x.as_mut_slice().copy_from_slice(block.as_slice());
    column_norms_squared_into(comm, &o.x, &mut o.norm_partial, &mut o.norms)?;
    for (j, nn) in o.norms.iter().enumerate() {
        let len = nn.sqrt();
        for v in o.x.as_mut_slice()[j * rows..(j + 1) * rows].iter_mut() {
            *v = if len == 0.0 { 0.0 } else { *v / len };
        }
    }
    let against = against.filter(|q| q.cols() > 0);
    let mut skip_second = false;
    for pass in 0..if against.is_some() { 2 } else { 1 } {
        if (pass > 0 && skip_second) || o.x.cols() == 0 {
            break;
        }
        let count = o.x.cols();
        if let Some(q) = against {
            fit(&mut o.overlap, q.cols(), count)?;
            fit(&mut o.scratch, rows, count)?;
            fit(&mut o.partial, q.cols(), count)?;
            project_out_slices(
                comm,
                gemm,
                q.as_slice(),
                q.cols(),
                &mut o.x,
                &mut o.overlap,
                &mut o.scratch,
                &mut o.partial,
                tuning.wide_gemm,
            )?;
        }
        let min_len = panel_gs_into(comm, gemm, panel, tuning.wide_gemm, o)?;
        if pass == 0
            && against.is_some()
            && let Some(keep) = tuning.dgks
        {
            skip_second = o.y.cols() == count && min_len >= keep;
        }
        std::mem::swap(&mut o.x, &mut o.y);
    }
    fit(out, rows, o.x.cols())?;
    out.as_mut_slice().copy_from_slice(o.x.as_slice());
    Ok(())
}

/// Append the columns of `extra` to `a` in place (storage permitting).
fn append_in_place(a: &mut RealMatrix, extra: &RealMatrix) -> Result<(), SolverError> {
    let (rows, m) = (a.rows(), a.cols());
    if extra.rows() != rows {
        return Err(SolverError::InvalidParameter("append dimensions"));
    }
    fit(a, rows, m + extra.cols())?;
    a.as_mut_slice()[m * rows..].copy_from_slice(extra.as_slice());
    Ok(())
}

/// Workspace form of `solve_core` (see the module documentation). `apply`
/// writes H x into a caller-shaped output (local rows x columns of x).
#[allow(clippy::too_many_arguments)]
pub(super) fn solve_core_ws<C: Collective>(
    solver: &crate::RealDavidson,
    comm: &C,
    layout: &Layout,
    supplied_pool: Option<Arc<WorkerPool>>,
    inner_threads: usize,
    apply: &mut dyn FnMut(&RealMatrix, &mut RealMatrix) -> Result<(), SolverError>,
    initial: Option<&RealMatrix>,
    diag_local: &[f64],
    full_diagonal: Option<&[f64]>,
    policy: BandTolerances,
    ritz: bool,
    certify: bool,
    capacity_per_band: usize,
    tuning: DistributedTuning,
) -> Result<RealEigenSolution, SolverError> {
    let mut ws = WORKSPACE
        .with(|w| w.borrow_mut().take())
        .unwrap_or_else(|| Box::new(DistributedDavidsonWorkspace::new()));
    SOLVE_GROWTHS.with(|g| g.set(0));
    let result = solve_on(
        &mut ws,
        solver,
        comm,
        layout,
        supplied_pool,
        inner_threads,
        apply,
        initial,
        diag_local,
        full_diagonal,
        policy,
        ritz,
        certify,
        capacity_per_band,
        tuning,
    );
    WS_COUNTS[0].fetch_add(1, Ordering::Relaxed);
    if SOLVE_GROWTHS.with(|g| g.get()) == 0 {
        WS_COUNTS[1].fetch_add(1, Ordering::Relaxed);
    }
    WORKSPACE.with(|w| *w.borrow_mut() = Some(ws));
    result
}

#[allow(clippy::too_many_arguments)]
fn solve_on<C: Collective>(
    ws: &mut DistributedDavidsonWorkspace,
    solver: &crate::RealDavidson,
    comm: &C,
    layout: &Layout,
    supplied_pool: Option<Arc<WorkerPool>>,
    inner_threads: usize,
    apply: &mut dyn FnMut(&RealMatrix, &mut RealMatrix) -> Result<(), SolverError>,
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
    if solver.bands == 0 || solver.bands > n {
        return Err(SolverError::InvalidParameter("bands or diagonal length"));
    }
    if !solver.tolerance.is_finite()
        || solver.tolerance <= 0.0
        || solver.max_iterations == 0
        || policy.strict_bands == 0
        || policy.strict_bands > solver.bands
        || !policy.buffer_tolerance.is_finite()
        || policy.buffer_tolerance < solver.tolerance
    {
        return Err(SolverError::InvalidParameter(
            "invalid tolerance, iteration limit or strict prefix",
        ));
    }
    let bands = solver.bands;
    let capacity = solver.max_subspace.max(bands * capacity_per_band).min(n);
    let nl = layout.nl();
    let max_rows = layout.counts.iter().copied().max().unwrap_or(1);
    let pool = pool_of(ws, supplied_pool, inner_threads)?;
    let key = (
        max_rows.max(1),
        capacity.max(1),
        Arc::as_ptr(&pool) as usize,
    );
    if ws.gemm.as_ref().map(|(k, _)| *k) != Some(key) {
        grown(max_rows.max(1).div_ceil(256) * capacity * capacity * 8);
        ws.gemm = Some((
            key,
            RealGemmWorkspace::with_pool(max_rows.max(1), capacity.max(1), pool.clone())?,
        ));
    }
    let root_mode = tuning.eigen_root && comm.size() > 1;
    let solves_here = !root_mode || comm.rank() == 0;
    if solves_here && ws.eig.as_ref().map(|(c, _)| *c) != Some(capacity) {
        grown(3 * capacity * capacity * 8);
        ws.eig = Some((capacity, RealEighWorkspace::new(capacity)?));
    }
    let use_b5 = tuning.bordered && solves_here;
    let b5_key = (capacity, Arc::as_ptr(&pool) as usize);
    if use_b5 && ws.b5.as_ref().map(|(k, _, _)| *k) != Some(b5_key) {
        grown(6 * capacity * capacity * 8);
        ws.b5 = Some((
            b5_key,
            Bordered::new(capacity),
            RealGemmWorkspace::with_pool(capacity.max(1), capacity.max(1), pool.clone())?,
        ));
    }
    if let Some((_, b, _)) = ws.b5.as_mut() {
        b.invalidate();
    }
    if root_mode {
        fit_vec(&mut ws.eigen_buffer, 1 + bands + capacity * bands);
    }
    if ws.held.len() < capacity * capacity {
        grown((capacity * capacity - ws.held.len()) * 8);
        ws.held.resize(capacity * capacity, 0.0);
    }
    let ld = capacity;
    fit_vec(&mut ws.values, capacity);
    for m in [&mut ws.basis, &mut ws.products] {
        reserve(m, nl * capacity)?;
    }
    for m in [&mut ws.reduced, &mut ws.coefficients] {
        reserve(m, capacity * capacity)?;
    }
    for m in [&mut ws.block, &mut ws.partial] {
        reserve(m, capacity * bands)?;
    }
    for m in [
        &mut ws.vectors,
        &mut ws.rotated,
        &mut ws.residual,
        &mut ws.corrections,
        &mut ws.extra,
        &mut ws.extra_products,
    ] {
        reserve(m, nl * bands)?;
    }
    let mut rng = Rng(solver.seed);
    let panel = tuning.panel_width(max_rows.max(1));
    let wide = tuning.wide_gemm;

    let DistributedDavidsonWorkspace {
        gemm,
        eig,
        b5,
        eigen_buffer,
        basis,
        products,
        held,
        reduced,
        coefficients,
        values,
        block,
        partial,
        vectors,
        rotated,
        residual,
        actual,
        corrections,
        extra,
        extra_products,
        defect,
        norms,
        norm_partial,
        converged,
        selected,
        ortho,
        ..
    } = ws;
    let gemm = &mut gemm.as_mut().expect("GEMM workspace").1;

    // Initial basis.
    let cold;
    let start = match initial {
        Some(x) => x,
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
            cold = free_electron_initial_guess(&full, bands)?;
            &cold
        }
    };
    ortho_into(
        comm,
        gemm,
        None,
        &local_rows(layout, start)?,
        &tuning,
        panel,
        ortho,
        basis,
    )?;
    while basis.cols() < bands {
        let w = bands - basis.cols();
        let random = RealMatrix::new(n, w, (0..n * w).map(|_| rng.next()).collect())?;
        ortho_into(
            comm,
            gemm,
            Some(basis),
            &local_rows(layout, &random)?,
            &tuning,
            panel,
            ortho,
            extra,
        )?;
        if extra.cols() == 0 {
            return Err(SolverError::InvalidParameter(
                "unable to complete initial basis",
            ));
        }
        append_in_place(basis, extra)?;
    }
    let mut applications = 0usize;
    fit(products, nl, basis.cols())?;
    apply(basis, products)?;
    applications += basis.cols();

    let threshold = |j: usize| {
        if j < policy.strict_bands - usize::from(ritz) {
            solver.tolerance
        } else {
            policy.buffer_tolerance
        }
    };
    let meets_policy = |r: &[f64]| r.iter().enumerate().all(|(j, &x)| x <= threshold(j));
    let meets_strict = |r: &[f64]| {
        r[..policy.strict_bands]
            .iter()
            .all(|&x| x <= solver.tolerance)
    };
    let mut history = Vec::with_capacity(solver.max_iterations);
    let (mut restarts, mut contractions, mut certificates) = (0usize, 0usize, 0usize);
    let mut previous = vec![f64::NAN; bands];
    let mut reduced_width = 0usize;
    let mut max_basis = 0usize;
    let buffer: Vec<bool> = (0..bands).map(|j| j >= policy.strict_bands).collect();

    for iteration in 1..=solver.max_iterations {
        call_probe(iteration);
        let m = basis.cols();
        max_basis = max_basis.max(m);
        // Reduced matrix: only the border of the columns added since the last
        // iteration is formed (Gram of the basis with the fresh products) and
        // written, symmetrised, into the fixed-leading-dimension storage.
        let t0 = Instant::now();
        let fresh = m - reduced_width;
        let prefix = reduced_width;
        fit(block, m, fresh)?;
        fit(partial, m, fresh)?;
        gram_slices(
            comm,
            gemm,
            basis.as_slice(),
            m,
            &products.as_slice()[reduced_width * nl..m * nl],
            fresh,
            nl,
            partial,
            block,
        )?;
        for c in reduced_width..m {
            for r in 0..=c {
                let v = block[(r, c - reduced_width)];
                held[r + c * ld] = v;
                held[c + r * ld] = v;
            }
        }
        reduced_width = m;
        fit(reduced, m, m)?;
        for c in 0..m {
            reduced.as_mut_slice()[c * m..(c + 1) * m].copy_from_slice(&held[c * ld..c * ld + m]);
        }
        stage_add(0, t0);
        fit(coefficients, m, m)?;
        let values = &mut values[..m];
        let t0 = Instant::now();
        record_solve(m, fresh);
        let mut status: Result<(), SolverError> = Ok(());
        if solves_here && let Some((_, eig)) = eig.as_mut() {
            status = solve_reduced(
                eig,
                if use_b5 {
                    b5.as_mut().map(|(_, b, g)| (b, g))
                } else {
                    None
                },
                reduced,
                prefix,
                values,
                coefficients,
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
                    coefficients.as_mut_slice()[..m * bands].copy_from_slice(&buffer[1 + bands..]);
                }
            }
            SOLVE_COUNTS[3].fetch_add(1, Ordering::Relaxed);
        }
        status?;
        stage_add(1, t0);
        let t0 = Instant::now();
        fit(vectors, nl, bands)?;
        fit(rotated, nl, bands)?;
        let lowest = &coefficients.as_slice()[..m * bands];
        rotate_slices_with(gemm, basis.as_slice(), nl, m, lowest, bands, vectors, wide)?;
        rotate_slices_with(
            gemm,
            products.as_slice(),
            nl,
            m,
            lowest,
            bands,
            rotated,
            wide,
        )?;
        stage_add(2, t0);
        let t0 = Instant::now();
        residual_into(
            comm,
            rotated,
            vectors,
            &values[..bands],
            tuning.lean,
            residual,
            norms,
            norm_partial,
        )?;
        stage_add(3, t0);
        let mut certified = false;
        let mut trusted = true;
        if !ritz
            && (certify || iteration == solver.max_iterations)
            && (meets_policy(norms) || (iteration == solver.max_iterations && meets_strict(norms)))
        {
            certificates += 1;
            fit(actual, nl, bands)?;
            apply(vectors, actual)?;
            applications += bands;
            residual_into(
                comm,
                actual,
                vectors,
                &values[..bands],
                tuning.lean,
                residual,
                norms,
                norm_partial,
            )?;
            certified = true;
            trusted = meets_policy(norms);
        }
        let maximum = norms.iter().copied().fold(0.0, f64::max);
        history.push(maximum);
        converged.clear();
        converged.extend((0..bands).map(|j| {
            if ritz {
                iteration > 1 && (values[j] - previous[j]).abs() < threshold(j)
            } else {
                norms[j] <= threshold(j)
            }
        }));
        previous.copy_from_slice(&values[..bands]);
        let finished = converged.iter().all(|&b| b)
            || (iteration == solver.max_iterations
                && converged[..policy.strict_bands].iter().all(|&b| b));
        if finished {
            let full = gather_full(comm, layout, vectors)?;
            return Ok(RealEigenSolution {
                eigenvalues: values[..bands].to_vec(),
                eigenvectors: full,
                residuals: norms.clone(),
                converged: converged.clone(),
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
        if iteration == solver.max_iterations {
            return Err(SolverError::NotConverged {
                iterations: iteration,
                max_residual: maximum,
                residuals: norms.clone(),
            });
        }
        if ritz && iteration == 1 && basis.cols() == n {
            continue;
        }
        // Corrections for the unconverged roots (row-local preconditioning).
        selected.clear();
        selected.extend(
            (0..bands)
                .filter(|&j| !converged[j])
                .take(bands.min(capacity - bands)),
        );
        fit(corrections, nl, selected.len())?;
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
                let r = column_of(residual, j);
                for ((c, &d), &r) in corrections.as_mut_slice()[k * nl..(k + 1) * nl]
                    .iter_mut()
                    .zip(diag_local)
                    .zip(r)
                {
                    *c = precondition(d, values[j], r);
                }
            } else {
                for i in 0..nl {
                    corrections[(i, k)] = precondition(diag_local[i], values[j], residual[(i, j)]);
                }
            }
        }
        if !trusted || basis.cols() + corrections.cols() > capacity {
            restarts += 1;
            reduced_width = 0;
            fit(defect, bands, bands)?;
            fit(partial, bands, bands)?;
            gram_slices(
                comm,
                gemm,
                vectors.as_slice(),
                bands,
                vectors.as_slice(),
                bands,
                nl,
                partial,
                defect,
            )?;
            let mut worst = 0.0_f64;
            for j in 0..bands {
                for i in 0..bands {
                    worst = worst.max((defect[(i, j)] - if i == j { 1.0 } else { 0.0 }).abs());
                }
            }
            if trusted && worst <= 1e-12 {
                contractions += 1;
                fit(basis, nl, bands)?;
                basis.as_mut_slice().copy_from_slice(vectors.as_slice());
                fit(products, nl, bands)?;
                products.as_mut_slice().copy_from_slice(rotated.as_slice());
            } else {
                ortho_into(comm, gemm, None, vectors, &tuning, panel, ortho, basis)?;
                fit(products, nl, basis.cols())?;
                apply(basis, products)?;
                applications += basis.cols();
            }
        }
        let t0 = Instant::now();
        ortho_into(
            comm,
            gemm,
            Some(basis),
            corrections,
            &tuning,
            panel,
            ortho,
            extra,
        )?;
        stage_add(4, t0);
        if extra.cols() == 0 {
            let w = bands.min(capacity - basis.cols());
            let random = RealMatrix::new(n, w, (0..n * w).map(|_| rng.next()).collect())?;
            ortho_into(
                comm,
                gemm,
                Some(basis),
                &local_rows(layout, &random)?,
                &tuning,
                panel,
                ortho,
                extra,
            )?;
        }
        if extra.cols() == 0 {
            return Err(SolverError::NotConverged {
                iterations: iteration,
                max_residual: maximum,
                residuals: norms.clone(),
            });
        }
        fit(extra_products, nl, extra.cols())?;
        apply(extra, extra_products)?;
        applications += extra.cols();
        append_in_place(products, extra_products)?;
        append_in_place(basis, extra)?;
    }
    unreachable!("validated positive iteration limit")
}
