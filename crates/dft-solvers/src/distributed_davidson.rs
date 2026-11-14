//! Row-block distributed real Davidson (ADR 0041, M1c). Interface in and out is
//! replicated (full initial vectors, full diagonal, full eigenvectors on every
//! rank); inside, the basis and the products are plane-wave row blocks, the dense
//! algebra is local GEMM plus the deterministic tree sum, and H is applied to
//! column pairs owned by one rank after gathering the correction block. With one
//! rank the caller should use the threads-only `RealDavidson` (bit-identical by
//! construction); this solver is used for two or more ranks.
use crate::{
    BandTolerances, RealLinearOperator, SolverError,
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

fn comm_error(e: dft_mpi::MpiError) -> SolverError {
    SolverError::Application(e.to_string())
}

use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

const STAGES: [&str; 8] = [
    "reduced_gram",
    "reduced_eigen",
    "rotate",
    "residual",
    "orthonormalize",
    "apply_total",
    "apply_operator",
    "other_solve",
];
static STAGE_NANOS: [AtomicU64; 8] = [const { AtomicU64::new(0) }; 8];

/// Wall seconds of this process in each stage of the distributed solver (cheap
/// per-stage accumulation, no I/O in the loop; read once for the run report).
pub fn stage_seconds() -> Vec<(&'static str, f64)> {
    (0..8)
        .map(|i| {
            (
                STAGES[i],
                STAGE_NANOS[i].load(Ordering::Relaxed) as f64 * 1e-9,
            )
        })
        .collect()
}
fn stage_add(i: usize, start: Instant) {
    STAGE_NANOS[i].fetch_add(start.elapsed().as_nanos() as u64, Ordering::Relaxed);
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
    let apply_start = Instant::now();
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
        let op_start = Instant::now();
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
        let (y, min_len) = panel_gram_schmidt(comm, ws, &x, panel, tuning.wide_gemm)?;
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
    wide: bool,
) -> Result<(RealMatrix, f64), SolverError> {
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
            for _ in 0..2 {
                project_out_with(comm, ws, &q, &mut p, &mut overlap, &mut scratch, wide)?;
            }
        }
        let first = kept;
        let mut local = Vec::with_capacity(w);
        let mut sums = Vec::with_capacity(w);
        let mut v = vec![0.0; rows];
        for jj in 0..w {
            v.copy_from_slice(column_of(&p, jj));
            for _ in 0..2 {
                let k = kept - first;
                if k > 0 {
                    local.clear();
                    for c in first..kept {
                        let q = &out.as_slice()[c * rows..(c + 1) * rows];
                        local.push(q.iter().zip(&v).map(|(a, b)| a * b).sum::<f64>());
                    }
                    sums.clear();
                    sums.resize(k, 0.0);
                    comm.tree_sum_f64(&local, &mut sums).map_err(comm_error)?;
                    for (c, s) in (first..kept).zip(&sums) {
                        let q = &out.as_slice()[c * rows..(c + 1) * rows];
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
        let n = op.dimension();
        if diagonal.len() != n {
            return Err(SolverError::InvalidParameter("bands or diagonal length"));
        }
        let layout = Layout::contiguous(n, comm.rank(), comm.size());
        let nl = layout.nl();
        let offset = layout.index[comm.rank()].first().copied().unwrap_or(0);
        let diag_local = &diagonal[offset..offset + nl];
        let pool = match op.worker_pool() {
            Some(p) => p,
            None => WorkerPool::new(op.inner_threads())
                .map_err(|e| SolverError::Application(e.to_string()))?,
        };
        let mut apply = |x: &RealMatrix| apply_distributed(comm, &layout, op, x);
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
            DistributedTuning::current(),
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
        let n = op.global_len();
        let mine: Vec<usize> = op.local_rows().to_vec();
        if diagonal_local.len() != mine.len() || mine.iter().any(|&r| r >= n) {
            return Err(SolverError::InvalidParameter("bands or diagonal length"));
        }
        // Setup (once per solve): every rank learns the row lists of all ranks.
        let ranks = comm.size();
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
        let layout = Layout {
            n,
            rank: comm.rank(),
            counts,
            index,
        };
        let pool = match op.worker_pool() {
            Some(p) => p,
            None => WorkerPool::new(op.inner_threads())
                .map_err(|e| SolverError::Application(e.to_string()))?,
        };
        let nl = layout.nl();
        let mut apply = |x: &RealMatrix| -> Result<RealMatrix, SolverError> {
            let start = Instant::now();
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
        self.solve_core(
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
        )
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
        let mut ws = RealGemmWorkspace::with_pool(max_rows.max(1), capacity.max(1), pool)?;
        let mut eig = RealEighWorkspace::new(capacity)?;
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
            let m = basis.cols();
            max_basis = max_basis.max(m);
            // Reduced matrix V^T HV (replicated), symmetrised from the upper triangle.
            // Incremental, as in the serial solver: only the columns added since the last
            // iteration are formed (m x new products); the rest is reused.
            let t0 = Instant::now();
            let mut reduced = RealMatrix::zeros(m, m)?;
            for c in 0..reduced_width {
                for r in 0..reduced_width {
                    reduced[(r, c)] = reduced_prev[(r, c)];
                }
            }
            let fresh = m - reduced_width;
            let mut fresh_products = RealMatrix::zeros(nl, fresh)?;
            fresh_products
                .as_mut_slice()
                .copy_from_slice(&products.as_slice()[reduced_width * nl..m * nl]);
            let mut block = RealMatrix::zeros(m, fresh)?;
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
            let t0 = Instant::now();
            eig.solve_into(&reduced, &mut values, &mut coefficients)?;
            stage_add(1, t0);
            let t0 = Instant::now();
            let lowest = RealMatrix::new(m, bands, coefficients.as_slice()[..m * bands].to_vec())?;
            let mut vectors = RealMatrix::zeros(nl, bands)?;
            let mut rotated = RealMatrix::zeros(nl, bands)?;
            rotate_local_with(&mut ws, &basis, &lowest, &mut vectors, wide)?;
            rotate_local_with(&mut ws, &products, &lowest, &mut rotated, wide)?;
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
                let norms = column_norms_squared(comm, &r)?
                    .into_iter()
                    .map(f64::sqrt)
                    .collect();
                Ok((r, norms))
            };
            let t0 = Instant::now();
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
            let t0 = Instant::now();
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
