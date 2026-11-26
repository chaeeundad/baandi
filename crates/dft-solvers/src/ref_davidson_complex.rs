//! Complex, owned-row reference block Davidson. Gamma keeps its real solver.
//! Uses the reference Ritz-change/empty-band rules and Hermitian overlap metric.
//! All ranks process the same k point; there are no k-parallel pools here.
//!
//! GK3 (`BAANDI_CDAV_REUSE`, default on, bit-identical): the workspace is kept
//! per calling thread across calls (k points, SCF maps) like the real path's
//! row cache; error flags ride as one extra element of the next numerical
//! reduction (or the reduced-eigenpair broadcast) instead of their own
//! collectives; GEMM operands already known finite are not rescanned. Off
//! restores per-call allocation, one agreement per product and full scans.
use crate::ref_davidson::accounting::{Counted, clock, elapsed};
use crate::{
    SolverError,
    ref_davidson::{self, Iteration, Options, Record, ReducedPlacement},
};
use dft_mpi::Collective;
use dft_numerics::{
    linalg::{
        ComplexMatrix as Matrix, ComplexOverlapFactor, DEPENDENCE_TOLERANCE, GemmWorkspace,
        MatrixOp,
    },
    pool::WorkerPool,
};
use num_complex::Complex64 as C;
use std::{cell::RefCell, sync::Arc};

pub trait DistributedComplexOperator {
    fn global_len(&self) -> usize;
    fn local_rows(&self) -> &[usize];
    fn worker_pool(&self) -> Option<Arc<WorkerPool>> {
        None
    }
    /// Implementations agree local failures before entering later internal collectives.
    /// `out` has the shape of `x` and must be overwritten completely: its
    /// previous contents are unspecified (reused solver storage).
    fn apply_block<K: Collective>(
        &mut self,
        comm: &K,
        x: &Matrix,
        out: &mut Matrix,
    ) -> Result<(), SolverError>;
}
#[derive(Debug)]
pub struct Solution {
    pub eigenvalues: Vec<f64>,
    pub eigenvectors: Matrix,
    pub hx: Option<Matrix>,
    pub converged: Vec<bool>,
    pub ritz_changes_ha: Vec<Vec<f64>>,
    pub ritz_values_ha: Vec<Vec<[f64; 2]>>,
    pub band_thresholds_ha: Vec<f64>,
    pub record: Record,
}
fn error(e: impl std::fmt::Display) -> SolverError {
    SolverError::Application(e.to_string())
}
const DISAGREED: &str = "complex reference Davidson invalid input on a rank";
fn agree<K: Collective>(comm: &K, bad: bool) -> Result<(), SolverError> {
    let mut out = [0.];
    comm.tree_sum_f64(&[f64::from(bad)], &mut out)
        .map_err(error)?;
    if out[0] != 0. {
        return Err(error(DISAGREED));
    }
    Ok(())
}
/// Error-flag agreement. Unmerged (switch off) every local check is agreed at
/// once by its own collective. Merged, local failures accumulate and travel as
/// one extra element of the next numerical reduction or broadcast, which every
/// rank reaches with the same call sequence; no decision that differs between
/// ranks reads failed data before that collective. The data elements keep
/// their reduction association: the tree sum is element-wise.
struct Agreement {
    merge: bool,
    pending: bool,
}
impl Agreement {
    fn note<K: Collective>(&mut self, comm: &K, bad: bool) -> Result<(), SolverError> {
        if self.merge {
            self.pending |= bad;
            Ok(())
        } else {
            agree(comm, bad)
        }
    }
    /// A standalone agreement where no numerical collective follows.
    fn flush<K: Collective>(&mut self, comm: &K) -> Result<(), SolverError> {
        if self.merge {
            agree(comm, std::mem::take(&mut self.pending))
        } else {
            Ok(())
        }
    }
    /// Tree sum of `local[..len]` carrying the pending flag (merged only) in
    /// element `len`; `local` and `all` hold at least `len + 1` elements.
    fn sum<K: Collective>(
        &mut self,
        comm: &K,
        local: &mut [f64],
        all: &mut [f64],
        len: usize,
    ) -> Result<(), SolverError> {
        let total = if self.merge {
            local[len] = f64::from(std::mem::take(&mut self.pending));
            len + 1
        } else {
            len
        };
        comm.tree_sum_f64(&local[..total], &mut all[..total])
            .map_err(error)?;
        if total > len && all[len] != 0. {
            return Err(error(DISAGREED));
        }
        Ok(())
    }
}
fn reduce_into<K: Collective>(
    comm: &K,
    a: &mut Matrix,
    local: &mut [f64],
    all: &mut [f64],
    ag: &mut Agreement,
) -> Result<(), SolverError> {
    let len = 2 * a.as_slice().len();
    for (p, z) in local[..len]
        .as_chunks_mut::<2>()
        .0
        .iter_mut()
        .zip(a.as_slice())
    {
        *p = [z.re, z.im];
    }
    ag.sum(comm, local, all, len)?;
    for (z, p) in a
        .as_mut_slice()
        .iter_mut()
        .zip(all[..len].as_chunks::<2>().0)
    {
        *z = C::new(p[0], p[1]);
    }
    Ok(())
}
#[cfg(test)]
fn reduce<K: Collective>(comm: &K, a: &mut Matrix) -> Result<(), SolverError> {
    let len = 2 * a.as_slice().len();
    let mut ag = Agreement {
        merge: false,
        pending: false,
    };
    reduce_into(comm, a, &mut vec![0.; len], &mut vec![0.; len], &mut ag)
}
struct Products {
    gemm: GemmWorkspace,
    local: Vec<f64>,
    all: Vec<f64>,
    /// Skip finite scans of operands flagged known finite (GK3 arm).
    skip_scans: bool,
}
impl Products {
    /// `a_finite`/`b_finite`: the operand was checked finite since its last
    /// mutation (see the invariant in `solve_with_initial_hx`).
    #[allow(clippy::too_many_arguments)]
    fn multiply(
        &mut self,
        a: &Matrix,
        oa: MatrixOp,
        a_finite: bool,
        b: &Matrix,
        ob: MatrixOp,
        b_finite: bool,
        out: &mut Matrix,
    ) -> Result<(), dft_numerics::linalg::LinalgError> {
        if self.skip_scans {
            self.gemm
                .multiply_into_known_finite(a, oa, a_finite, b, ob, b_finite, out)
        } else {
            self.gemm.multiply_into(a, oa, b, ob, out)
        }
    }
    #[allow(clippy::too_many_arguments)]
    fn gram<K: Collective>(
        &mut self,
        comm: &K,
        ag: &mut Agreement,
        (a, a_finite): (&Matrix, bool),
        (b, b_finite): (&Matrix, bool),
        out: &mut Matrix,
        r: &mut Iteration,
    ) -> Result<(), SolverError> {
        let t = clock(r.profile);
        out.reshape_within_capacity(a.cols(), b.cols())?;
        let product = self.multiply(
            a,
            MatrixOp::Adjoint,
            a_finite,
            b,
            MatrixOp::Identity,
            b_finite,
            out,
        );
        ag.note(comm, product.is_err())?;
        r.gemm("complex_gram", a.cols(), b.cols(), a.rows());
        if let Some(g) = r.gemms.last_mut() {
            g.flops *= 4;
        }
        reduce_into(comm, out, &mut self.local, &mut self.all, ag)?;
        // Agreed above: a failed product on any rank already returned.
        product?;
        r.seconds[1] += elapsed(t);
        Ok(())
    }
    /// Local rotation; its failure flag is agreed by `ag` (at once when
    /// unmerged, else with the next collective).
    #[allow(clippy::too_many_arguments)]
    fn rotate<K: Collective>(
        &mut self,
        comm: &K,
        ag: &mut Agreement,
        (a, a_finite): (&Matrix, bool),
        c: &Matrix,
        out: &mut Matrix,
        r: &mut Iteration,
    ) -> Result<(), SolverError> {
        out.resize_columns(c.cols())?;
        let product = self.multiply(
            a,
            MatrixOp::Identity,
            a_finite,
            c,
            MatrixOp::Identity,
            false,
            out,
        );
        ag.note(comm, product.is_err())?;
        r.gemm("complex_rotate", a.rows(), c.cols(), a.cols());
        if let Some(g) = r.gemms.last_mut() {
            g.flops *= 4;
        }
        Ok(())
    }
}
fn hermitize(a: &mut Matrix) {
    for j in 0..a.cols() {
        a[(j, j)].im = 0.;
        for i in 0..j {
            let z = (a[(i, j)] + a[(j, i)].conj()) * 0.5;
            a[(i, j)] = z;
            a[(j, i)] = z.conj();
        }
    }
}
fn identity(a: &mut Matrix, n: usize) -> Result<(), SolverError> {
    a.reshape_within_capacity(n, n)?;
    a.as_mut_slice().fill(C::default());
    for i in 0..n {
        a[(i, i)] = C::new(1., 0.);
    }
    Ok(())
}
fn factor(
    f: &mut ComplexOverlapFactor,
    s: &Matrix,
    r: &mut Iteration,
) -> Result<bool, SolverError> {
    let t = clock(r.profile);
    let ok = f.factor(s, crate::real_davidson::PIVOT_RATIO_MIN)?;
    r.seconds[2] += elapsed(t);
    Ok(ok)
}
#[allow(clippy::too_many_arguments)]
fn eigen<K: Collective>(
    comm: &K,
    ag: &mut Agreement,
    f: &mut ComplexOverlapFactor,
    h: &Matrix,
    o: Options,
    r: &mut Iteration,
    values: &mut [f64],
    c: &mut Matrix,
    packed: &mut [f64],
    canonical: &mut [C],
) -> Result<(), SolverError> {
    let t = clock(o.profile);
    let root = o.placement == ReducedPlacement::Root;
    c.reshape_within_capacity(h.rows(), o.bands)?;
    let result = if !root || comm.rank() == 0 {
        r.backend_solves += 1;
        f.solve_lowest_into(h, values, c)
    } else {
        Ok(())
    };
    if root && ag.merge {
        // Every caller reaches here right after a fused reduction, so no
        // flag is pending; only rank 0 can fail and its flag leads the one
        // broadcast of [flag, values, coefficients] (copies, identical bits).
        debug_assert!(!ag.pending);
        let n = values.len();
        let len = 2 * c.as_slice().len();
        if comm.rank() == 0 {
            packed[0] = f64::from(result.is_err());
            packed[1..=n].copy_from_slice(values);
            for (p, z) in packed[1 + n..1 + n + len]
                .as_chunks_mut::<2>()
                .0
                .iter_mut()
                .zip(c.as_slice())
            {
                *p = [z.re, z.im];
            }
        }
        comm.bcast_f64(&mut packed[..1 + n + len], 0)
            .map_err(error)?;
        if packed[0] != 0. {
            return Err(error(DISAGREED));
        }
        values.copy_from_slice(&packed[1..=n]);
        for (z, p) in c
            .as_mut_slice()
            .iter_mut()
            .zip(packed[1 + n..1 + n + len].as_chunks::<2>().0)
        {
            *z = C::new(p[0], p[1]);
        }
    } else {
        agree(comm, result.is_err())?;
        result?;
        if root {
            comm.bcast_f64(values, 0).map_err(error)?;
            let len = 2 * c.as_slice().len();
            for (p, z) in packed[..len]
                .as_chunks_mut::<2>()
                .0
                .iter_mut()
                .zip(c.as_slice())
            {
                *p = [z.re, z.im];
            }
            comm.bcast_f64(&mut packed[..len], 0).map_err(error)?;
            for (z, p) in c
                .as_mut_slice()
                .iter_mut()
                .zip(packed[..len].as_chunks::<2>().0)
            {
                *z = C::new(p[0], p[1]);
            }
        }
    }
    if !canonical.is_empty() {
        canonicalize_clusters(values, c, canonical);
    }
    r.reduced_solves += 1;
    r.reduced_order = h.rows();
    r.seconds[3] += elapsed(t);
    Ok(())
}
/// Choose a deterministic basis of a roundoff-degenerate Ritz cluster by
/// projecting ordered coordinate axes into its coefficient subspace. The
/// cluster transformation is unitary even for a generalized overlap metric.
/// Scratch is allocated by the caller before the Davidson iteration loop.
pub fn canonicalize_clusters(values: &[f64], c: &mut Matrix, scratch: &mut [C]) {
    let scale = values.iter().map(|v| v.abs()).fold(0., f64::max);
    let tolerance = 64. * f64::EPSILON * scale;
    let mut first = 0;
    while first < values.len() {
        let mut end = first + 1;
        while end < values.len() && values[end] - values[first] <= tolerance {
            end += 1;
        }
        let k = end - first;
        if k > 1 {
            let rows = c.rows();
            let (old, rest) = scratch.split_at_mut(rows * k);
            let (rotation, work) = rest.split_at_mut(k * k);
            for j in 0..k {
                old[j * rows..(j + 1) * rows].copy_from_slice(c.column(first + j));
            }
            let mut kept = 0;
            for row in 0..rows {
                for j in 0..k {
                    work[j] = old[j * rows + row].conj();
                }
                for _ in 0..2 {
                    for j in 0..kept {
                        let dot: C = (0..k).map(|i| rotation[j * k + i].conj() * work[i]).sum();
                        for i in 0..k {
                            work[i] -= dot * rotation[j * k + i];
                        }
                    }
                }
                let norm = work[..k].iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
                if norm <= 64. * f64::EPSILON {
                    continue;
                }
                for i in 0..k {
                    rotation[kept * k + i] = work[i] / norm;
                }
                kept += 1;
                if kept == k {
                    break;
                }
            }
            if kept == k {
                for j in 0..k {
                    for row in 0..rows {
                        c[(row, first + j)] = (0..k)
                            .map(|i| old[i * rows + row] * rotation[j * k + i])
                            .sum();
                    }
                }
            }
        }
        first = end;
    }
}
/// All per-call storage of one solve. Numerical content never carries over
/// between calls: every matrix is (re)shaped and completely written before it
/// is read, the factors are refactored, and the operator overwrites its
/// output. Only allocations (and their pages) are reused.
struct Workspace {
    cap: usize,
    n: usize,
    /// Rows reserved in the row matrices (largest local row count seen).
    rows: usize,
    long: usize,
    threads: usize,
    /// `None`: the solver's own one-worker pool (operator without a pool).
    pool: Option<Arc<WorkerPool>>,
    canonical_on: bool,
    products: Products,
    f: ComplexOverlapFactor,
    trial: ComplexOverlapFactor,
    s: Matrix,
    h: Matrix,
    c: Matrix,
    selected: Matrix,
    transform: Matrix,
    proj: Matrix,
    q: Matrix,
    hq: Matrix,
    v: Matrix,
    hv: Matrix,
    x: Matrix,
    hx: Matrix,
    e: Matrix,
    he: Matrix,
    canonical: Vec<C>,
}
thread_local! {
    static CACHE: RefCell<Option<Workspace>> = const { RefCell::new(None) };
}
/// Drop this thread's kept complex workspace (run end, before delivery).
pub(crate) fn release_workspace() {
    CACHE.with(|slot| *slot.borrow_mut() = None);
}
fn rows_matrix(rows: usize, capacity: usize) -> Result<Matrix, SolverError> {
    // Reserved, not zeroed: callers shape and write before reading.
    let mut a = Matrix::zeros(rows, 0)?;
    a.reserve_columns(capacity)?;
    Ok(a)
}
/// Product workspace with its identity key: (workspace, long, threads, pool).
type GemmParts = (GemmWorkspace, usize, usize, Option<Arc<WorkerPool>>);
fn gemm_for(
    rows: usize,
    cap: usize,
    pool: Option<Arc<WorkerPool>>,
) -> Result<GemmParts, SolverError> {
    let run = match &pool {
        Some(p) => p.clone(),
        None => WorkerPool::new(1).map_err(error)?,
    };
    let long = rows.max(cap);
    let threads = run.workers_created().max(1);
    Ok((
        GemmWorkspace::with_pool(long, cap, threads, run)?,
        long,
        threads,
        pool,
    ))
}
impl Workspace {
    fn new(
        rows: usize,
        n: usize,
        cap: usize,
        pool: Option<Arc<WorkerPool>>,
        canonical_on: bool,
    ) -> Result<Self, SolverError> {
        let (gemm, long, threads, pool) = gemm_for(rows, cap, pool)?;
        let s = Matrix::zeros(cap, cap)?;
        Ok(Self {
            cap,
            n,
            rows,
            long,
            threads,
            pool,
            canonical_on,
            // Gram (2 cap^2 + flag) and root broadcast (flag + n + 2 cap n).
            products: Products {
                gemm,
                local: vec![0.; 2 * cap * cap + cap + 1],
                all: vec![0.; 2 * cap * cap + cap + 1],
                skip_scans: false,
            },
            f: ComplexOverlapFactor::with_capacity(cap)?,
            trial: ComplexOverlapFactor::with_capacity(cap)?,
            h: s.clone(),
            c: s.clone(),
            selected: s.clone(),
            transform: s.clone(),
            proj: s.clone(),
            s,
            q: rows_matrix(rows, cap)?,
            hq: rows_matrix(rows, cap)?,
            v: rows_matrix(rows, cap)?,
            hv: rows_matrix(rows, cap)?,
            x: rows_matrix(rows, cap)?,
            hx: rows_matrix(rows, cap)?,
            e: rows_matrix(rows, n)?,
            he: rows_matrix(rows, n)?,
            canonical: if canonical_on {
                vec![C::default(); cap * n + n * n + n]
            } else {
                Vec::new()
            },
        })
    }
    /// Take this thread's kept workspace when it serves the request, else
    /// build one (rows grow monotonically so a k loop settles after one map).
    fn take(
        rows: usize,
        n: usize,
        cap: usize,
        pool: Option<Arc<WorkerPool>>,
        canonical_on: bool,
    ) -> Result<Self, SolverError> {
        let kept = CACHE.with(|slot| slot.borrow_mut().take());
        let Some(mut w) = kept else {
            return Self::new(rows, n, cap, pool, canonical_on);
        };
        if w.cap != cap || w.n != n || w.canonical_on != canonical_on || rows > w.rows {
            let rows = rows.max(if w.cap == cap && w.n == n { w.rows } else { 0 });
            drop(w);
            return Self::new(rows, n, cap, pool, canonical_on);
        }
        let same_pool = match (&w.pool, &pool) {
            (None, None) => true,
            (Some(a), Some(b)) => Arc::ptr_eq(a, b),
            _ => false,
        };
        let threads = pool
            .as_ref()
            .map_or(w.threads, |p| p.workers_created().max(1));
        if !same_pool || threads != w.threads || w.long < rows.max(cap) {
            // Tile/scheduling identity: rebuild only the product workspace.
            let (gemm, long, threads, pool) = gemm_for(w.rows, cap, pool)?;
            w.products.gemm = gemm;
            w.long = long;
            w.threads = threads;
            w.pool = pool;
        }
        Ok(w)
    }
}

pub fn solve<K: Collective, O: DistributedComplexOperator>(
    comm: &K,
    operator: &mut O,
    start: &Matrix,
    diagonal: &[f64],
    options: Options,
) -> Result<Solution, SolverError> {
    solve_with_initial_hx(comm, operator, start, diagonal, options, None)
}
pub fn solve_with_initial_hx<K: Collective, O: DistributedComplexOperator>(
    comm: &K,
    operator: &mut O,
    start: &Matrix,
    diagonal: &[f64],
    o: Options,
    initial_hx: Option<&Matrix>,
) -> Result<Solution, SolverError> {
    let total = clock(o.profile);
    let comm = &Counted::new(comm, o.profile);
    let reuse = dft_numerics::switches::on(dft_numerics::switches::Id::CDAV_REUSE);
    let mut ag = Agreement {
        merge: reuse,
        pending: false,
    };
    let n = o.bands;
    let dim = operator.global_len();
    let rows = operator.local_rows().len();
    let cap = ref_davidson::workspace::capacity(dim, n, o.capacity);
    let bad = n == 0
        || n > dim
        || cap < n.saturating_mul(2).min(dim)
        || o.occupied > n
        || o.max_iterations == 0
        || !o.threshold.is_finite()
        || o.threshold <= 0.
        || !o.band_threshold(n.saturating_sub(1)).is_finite()
        || start.rows() != rows
        || start.cols() != n
        || diagonal.len() != rows
        || start.validate_finite().is_err()
        || diagonal.iter().any(|x| !x.is_finite())
        || initial_hx
            .is_some_and(|x| x.rows() != rows || x.cols() != n || x.validate_finite().is_err());
    let has = f64::from(initial_hx.is_some());
    if reuse {
        // Input agreement and the initial-HX presence count in one sum.
        let mut out = [0.; 2];
        comm.tree_sum_f64(&[f64::from(bad), has], &mut out)
            .map_err(error)?;
        if out[0] != 0. {
            return Err(error(DISAGREED));
        }
        if out[1] != 0. && out[1] != comm.size() as f64 {
            return Err(error("complex initial HX differs across ranks"));
        }
    } else {
        agree(comm, bad)?;
        let mut out = [0.];
        comm.tree_sum_f64(&[has], &mut out).map_err(error)?;
        if out[0] != 0. && out[0] != comm.size() as f64 {
            return Err(error("complex initial HX differs across ranks"));
        }
    }
    let canonical_on =
        dft_numerics::switches::on(dft_numerics::switches::Id::REFSCF_CANONICAL_CLUSTER);
    let mut w = if reuse {
        Workspace::take(rows, n, cap, operator.worker_pool(), canonical_on)?
    } else {
        Workspace::new(rows, n, cap, operator.worker_pool(), canonical_on)?
    };
    w.products.skip_scans = reuse;
    let Workspace {
        products,
        f,
        trial,
        s,
        h,
        c,
        selected,
        transform,
        proj,
        q,
        hq,
        v,
        hv,
        x,
        hx,
        e,
        he,
        canonical,
        ..
    } = &mut w;
    for a in [&mut *q, &mut *hq, &mut *x, &mut *hx] {
        a.reshape_within_capacity(rows, 0)?;
    }
    for a in [&mut *v, &mut *hv, &mut *e, &mut *he] {
        a.reshape_within_capacity(rows, 0)?;
    }
    // Finite-operand invariant (GK3): at every product below, V and HV were
    // checked finite since their last mutation: START and HV by the input
    // agreement or the operator check, appended corrections by their finite
    // agreed norms (|e_ij| <= norm), appended HE by its explicit check,
    // rotated/repaired bases as validated product outputs or by a scanned
    // Gram. Small coefficient operands are always scanned.
    v.replace_columns(start)?;
    let mut values = vec![0.; n];
    let mut old = vec![0.; n];
    let mut record = Record::new(o.max_iterations, cap, o.profile, o.placement);
    let t = clock(o.profile);
    if let Some(h) = initial_hx {
        hv.replace_columns(h)?;
    } else {
        hv.resize_columns(n)?;
        let applied = operator.apply_block(comm, v, hv);
        ag.note(comm, applied.is_err() || hv.validate_finite().is_err())?;
        if !ag.merge {
            applied?;
        }
        record.setup.h_calls += 1;
        record.setup.h_columns += n;
    }
    record.setup.seconds[0] += elapsed(t);
    products.gram(comm, &mut ag, (v, true), (v, true), s, &mut record.setup)?;
    hermitize(s);
    if !factor(f, s, &mut record.setup)? {
        return Err(error(
            "complex reference Davidson start is numerically dependent",
        ));
    }
    products.gram(comm, &mut ag, (v, true), (hv, true), h, &mut record.setup)?;
    hermitize(h);
    eigen(
        comm,
        &mut ag,
        f,
        h,
        o,
        &mut record.setup,
        &mut values,
        c,
        &mut products.local,
        canonical,
    )?;
    record.setup.communication = comm.snapshot();
    let mut active: Vec<_> = (0..n).collect();
    let mut converged = vec![false; n];
    let mut changes = if o.profile {
        vec![vec![0.; n]; o.max_iterations]
    } else {
        Vec::new()
    };
    let mut history = if o.profile {
        vec![vec![[0.; 2]; n]; o.max_iterations]
    } else {
        Vec::new()
    };
    let bounds = diagonal
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &v| {
            (lo.min(v), hi.max(v))
        });
    let mut contracted = false;
    let mut contractions = 0usize;
    let mut reset_disabled = false;
    // One spare element for a fused agreement flag.
    let mut norms = vec![0.; n + 1];
    let mut sums = vec![0.; n + 1];
    let mut iterations = 0;
    for (step, r) in record.iterations.iter_mut().enumerate() {
        let before = comm.snapshot();
        let m = v.cols();
        r.basis_before = m;
        r.active_bands.extend_from_slice(&active);
        let t = clock(o.profile);
        selected.reshape_within_capacity(m, active.len())?;
        for (j, &b) in active.iter().enumerate() {
            selected.column_mut(j).copy_from_slice(c.column(b));
        }
        if contracted {
            x.resize_columns(active.len())?;
            hx.resize_columns(active.len())?;
            for (j, &band) in active.iter().enumerate() {
                x.column_mut(j).copy_from_slice(v.column(band));
                hx.column_mut(j).copy_from_slice(hv.column(band));
            }
        } else {
            products.rotate(comm, &mut ag, (v, true), selected, x, r)?;
            products.rotate(comm, &mut ag, (hv, true), selected, hx, r)?;
        }
        contracted = false;
        e.resize_columns(active.len())?;
        norms[..active.len()].fill(0.);
        for (j, &b) in active.iter().enumerate() {
            for (i, &d) in diagonal.iter().enumerate() {
                let denominator = if ref_davidson::stream::sqrt_safe(bounds, values[b]) {
                    ref_davidson::stream::sqrt_denominator(d, values[b])
                } else {
                    ref_davidson::preconditioner_denominator(d, values[b])
                };
                let z = (hx[(i, j)] - values[b] * x[(i, j)]) / denominator;
                e[(i, j)] = z;
                norms[j] += z.norm_sqr();
            }
        }
        ag.sum(comm, &mut norms, &mut sums, active.len())?;
        let mut kept = 0;
        for (j, &norm2) in sums[..active.len()].iter().enumerate() {
            let norm = norm2.sqrt();
            if !norm.is_finite() {
                return Err(error("nonfinite complex correction"));
            }
            if norm <= DEPENDENCE_TOLERANCE {
                r.dropped_corrections += 1;
                continue;
            }
            for i in 0..rows {
                e[(i, kept)] = e[(i, j)] / norm;
            }
            kept += 1;
        }
        e.resize_columns(kept)?;
        r.seconds[4] += elapsed(t);
        // Try the complete normalized correction block, then repair only if
        // its metric is dependent or the full physical space has been reached.
        let mut accepted = false;
        if m + kept <= cap && kept > 0 {
            v.append_columns(e)?;
            products.gram(comm, &mut ag, (v, true), (v, true), s, r)?;
            hermitize(s);
            accepted = factor(trial, s, r)?;
        }
        if !accepted {
            let t = clock(o.profile);
            let child_seconds = r.seconds[1] + r.seconds[2];
            r.repairs += 1;
            v.resize_columns(m)?;
            // Whiten the OLD basis and HV without reapplying H. Preserve the
            // physical Ritz block through a freshly measured coordinate map.
            // The rare repair path scans every operand.
            products.rotate(comm, &mut ag, (v, false), c, x, r)?;
            f.basis_transform_into(transform)?;
            products.rotate(comm, &mut ag, (v, false), transform, q, r)?;
            products.rotate(comm, &mut ag, (hv, false), transform, hq, r)?;
            products.gram(comm, &mut ag, (q, false), (x, false), c, r)?;
            std::mem::swap(v, q);
            std::mem::swap(hv, hq);
            for _ in 0..2 {
                products.gram(comm, &mut ag, (v, false), (e, false), proj, r)?;
                products.rotate(comm, &mut ag, (v, false), proj, x, r)?;
                for (a, b) in e.as_mut_slice().iter_mut().zip(x.as_slice()) {
                    *a -= b;
                }
            }
            ag.flush(comm)?;
            let input = e.cols();
            kept = 0;
            for j in 0..input {
                for _ in 0..2 {
                    for a in 0..kept {
                        let z: C = (0..rows).map(|i| e[(i, a)].conj() * e[(i, j)]).sum();
                        let mut sum = [0.; 2];
                        comm.tree_sum_f64(&[z.re, z.im], &mut sum).map_err(error)?;
                        let z = C::new(sum[0], sum[1]);
                        for i in 0..rows {
                            let q = e[(i, a)];
                            e[(i, j)] -= z * q;
                        }
                    }
                }
                let z: f64 = e.column(j).iter().map(|z| z.norm_sqr()).sum();
                let mut sum = [0.];
                comm.tree_sum_f64(&[z], &mut sum).map_err(error)?;
                let norm = sum[0].sqrt();
                if norm <= DEPENDENCE_TOLERANCE || m + kept == cap {
                    r.dropped_corrections += 1;
                    continue;
                }
                for i in 0..rows {
                    e[(i, kept)] = e[(i, j)] / norm;
                }
                kept += 1;
            }
            e.resize_columns(kept)?;
            v.append_columns(e)?;
            // Scanned: establishes the invariant for the repaired basis.
            products.gram(comm, &mut ag, (v, false), (v, false), s, r)?;
            hermitize(s);
            accepted = factor(trial, s, r)?;
            if !accepted {
                return Err(error("complex overlap remains dependent after repair"));
            }
            r.seconds[5] += (elapsed(t) - (r.seconds[1] + r.seconds[2] - child_seconds)).max(0.);
        }
        std::mem::swap(f, trial);
        old.copy_from_slice(&values);
        if kept > 0 {
            let t = clock(o.profile);
            he.resize_columns(kept)?;
            let applied = operator.apply_block(comm, e, he);
            ag.note(comm, applied.is_err() || he.validate_finite().is_err())?;
            if !ag.merge {
                applied?;
            }
            r.h_calls += 1;
            r.h_columns += kept;
            r.seconds[0] += elapsed(t);
            hv.append_columns(he)?;
            products.gram(comm, &mut ag, (v, true), (hv, true), h, r)?;
            hermitize(h);
            eigen(
                comm,
                &mut ag,
                f,
                h,
                o,
                r,
                &mut values,
                c,
                &mut products.local,
                canonical,
            )?;
        } else {
            r.reduced_order = m;
        }
        active.clear();
        for i in 0..n {
            let delta = (values[i] - old[i]).abs();
            converged[i] = delta < o.band_threshold(i);
            if !converged[i] {
                active.push(i);
            }
            if o.profile {
                changes[step][i] = delta;
                history[step][i] = [old[i], values[i]];
            }
        }
        r.unconverged = active.len();
        iterations = step + 1;
        if !active.is_empty() && iterations < o.max_iterations && v.cols() + active.len() > cap {
            let t = clock(o.profile);
            let child_seconds = r.seconds[1] + r.seconds[2];
            products.rotate(comm, &mut ag, (v, true), c, x, r)?;
            products.rotate(comm, &mut ag, (hv, true), c, hx, r)?;
            {
                std::mem::swap(v, x);
                std::mem::swap(hv, hx);
            }
            contractions += 1;
            contracted = true;
            identity(c, n)?;
            let mut ideal = false;
            {
                r.resets.attempts += 1;
                if reset_disabled {
                    r.resets.disabled_refreshes += 1;
                } else if contractions.is_multiple_of(ref_davidson::ideal_reset::REFRESH_INTERVAL) {
                    r.resets.periodic_refreshes += 1;
                } else {
                    for (j, norm) in norms[..n].iter_mut().enumerate() {
                        *norm = v.column(j).iter().map(|z| z.norm_sqr()).sum();
                    }
                    ag.sum(comm, &mut norms, &mut sums, n)?;
                    ideal = sums[..n].iter().all(|v| {
                        v.is_finite()
                            && (v - 1.).abs() <= ref_davidson::ideal_reset::ORTHOGONALITY_LIMIT
                    });
                    if ideal {
                        r.resets.ideal += 1;
                    } else {
                        reset_disabled = true;
                        r.resets.guard_fallbacks += 1;
                    }
                }
            }
            if ideal {
                identity(s, n)?;
                identity(h, n)?;
                for j in 0..n {
                    h[(j, j)] = C::new(values[j], 0.);
                }
            } else {
                products.gram(comm, &mut ag, (v, true), (v, true), s, r)?;
                hermitize(s);
                products.gram(comm, &mut ag, (v, true), (hv, true), h, r)?;
                hermitize(h);
                if !reset_disabled {
                    let mut defect = 0_f64;
                    for j in 0..n {
                        for i in 0..n {
                            defect = defect.max((s[(i, j)] - C::new(f64::from(i == j), 0.)).norm());
                        }
                    }
                    r.resets.max_measured_defect = defect;
                    if !defect.is_finite()
                        || defect > ref_davidson::ideal_reset::ORTHOGONALITY_LIMIT
                    {
                        reset_disabled = true;
                        r.resets.guard_fallbacks += 1;
                    }
                }
            }
            if !factor(f, s, r)? {
                return Err(error("complex restart metric is dependent"));
            }
            r.restarts += 1;
            r.seconds[6] += (elapsed(t) - (r.seconds[1] + r.seconds[2] - child_seconds)).max(0.);
        }
        r.communication = comm.snapshot().difference(before);
        if active.is_empty() || iterations == o.max_iterations {
            break;
        }
    }
    record.iterations.truncate(iterations);
    changes.truncate(iterations);
    history.truncate(iterations);
    let t = clock(o.profile);
    products.rotate(comm, &mut ag, (v, true), c, x, &mut record.finish)?;
    if o.return_hx {
        products.rotate(comm, &mut ag, (hv, true), c, hx, &mut record.finish)?;
    }
    ag.flush(comm)?;
    record.finish.seconds[7] += elapsed(t);
    record.communication = comm.snapshot();
    record.total_seconds = elapsed(total);
    record.summarize();
    let (eigenvectors, returned_hx) = if reuse {
        // Tight copies; the capacity-wide storage stays with the thread.
        let tight = |a: &Matrix| Matrix::new(a.rows(), a.cols(), a.as_slice().to_vec());
        let pair = (tight(x)?, o.return_hx.then(|| tight(hx)).transpose()?);
        CACHE.with(|slot| *slot.borrow_mut() = Some(w));
        pair
    } else {
        let hx = std::mem::replace(&mut w.hx, Matrix::zeros(0, 0)?);
        (
            std::mem::replace(&mut w.x, Matrix::zeros(0, 0)?),
            o.return_hx.then_some(hx),
        )
    };
    Ok(Solution {
        eigenvalues: values,
        eigenvectors,
        hx: returned_hx,
        converged,
        ritz_changes_ha: changes,
        ritz_values_ha: history,
        band_thresholds_ha: (0..n).map(|i| o.band_threshold(i)).collect(),
        record,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    struct Dense {
        h: Matrix,
        rows: Vec<usize>,
    }
    impl DistributedComplexOperator for Dense {
        fn global_len(&self) -> usize {
            self.h.rows()
        }
        fn local_rows(&self) -> &[usize] {
            &self.rows
        }
        fn apply_block<K: Collective>(
            &mut self,
            comm: &K,
            x: &Matrix,
            out: &mut Matrix,
        ) -> Result<(), SolverError> {
            let mut full = Matrix::zeros(self.h.rows(), x.cols())?;
            for j in 0..x.cols() {
                for (i, &g) in self.rows.iter().enumerate() {
                    full[(g, j)] = x[(i, j)];
                }
            }
            reduce(comm, &mut full)?;
            for j in 0..x.cols() {
                for (i, &g) in self.rows.iter().enumerate() {
                    out[(i, j)] = (0..self.h.rows())
                        .map(|a| self.h[(g, a)] * full[(a, j)])
                        .sum();
                }
            }
            Ok(())
        }
    }
    #[test]
    fn complex_reference_matches_hermitian_oracle_and_supplied_hx() {
        for ranks in [1, 3, 8] {
            dft_mpi::simulate(ranks, |comm| {
                let dim = 11;
                let n = 3;
                let mut h = Matrix::zeros(dim, dim).unwrap();
                for i in 0..dim {
                    h[(i, i)] = C::new(i as f64 * 0.3 - 0.7, 0.);
                    if i > 0 {
                        h[(i - 1, i)] = C::new(0.09, 0.07);
                        h[(i, i - 1)] = h[(i - 1, i)].conj();
                    }
                }
                let (exact, _) = dft_numerics::linalg::hermitian_eigh(&h).unwrap();
                let rows: Vec<_> = (0..dim).filter(|i| i % ranks == comm.rank()).collect();
                let mut start = Matrix::zeros(rows.len(), n).unwrap();
                for j in 0..n {
                    for (i, &g) in rows.iter().enumerate() {
                        start[(i, j)] = C::new(f64::from(g == j), 0.)
                            + C::new(
                                ((g + 1) * (j + 2)) as f64 * 0.001,
                                ((g + 2) * (j + 1)) as f64 * 0.002,
                            );
                    }
                }
                let diagonal: Vec<_> = rows.iter().map(|&i| h[(i, i)].re).collect();
                let mut op = Dense { h, rows };
                for (placement, flags) in [
                    (ReducedPlacement::Root, "off"),
                    (ReducedPlacement::Replicated, "on"),
                ] {
                    let config = dft_numerics::switches::RunConfig::from_pairs(&[]).unwrap();
                    let _scope = config.enter();
                    let mut o = Options::new(n, 1e-13, n);
                    o.max_iterations = 100;
                    o.return_hx = true;
                    o.profile = true;
                    o.placement = placement;
                    let a = solve(comm, &mut op, &start, &diagonal, o).unwrap();
                    assert!(a.converged.iter().all(|b| *b));
                    if flags == "on" {
                        assert!(a.record.resets.attempts > 0);
                    }
                    for (e, want) in a.eigenvalues.iter().zip(&exact) {
                        assert!((e - want).abs() < 1e-11, "{e} {want}");
                    }
                    let mut hs = Matrix::zeros(start.rows(), n).unwrap();
                    op.apply_block(comm, &start, &mut hs).unwrap();
                    let b = solve_with_initial_hx(comm, &mut op, &start, &diagonal, o, Some(&hs))
                        .unwrap();
                    assert_eq!(a.eigenvalues, b.eigenvalues);
                    assert_eq!(a.eigenvectors, b.eigenvectors);
                    assert_eq!(a.record.h_columns, b.record.h_columns + n);
                }
                let mut bad = start.clone();
                bad.column_mut(1).copy_from_slice(start.column(0));
                assert!(solve(comm, &mut op, &bad, &diagonal, Options::new(n, 1e-9, n)).is_err());
            });
        }
    }
    /// FD buffer bands (engine `electronic_map::buffer_bands`): a strongly
    /// coupled 16-state cluster cut in the middle by the retained count n = 20.
    /// The top retained Ritz pairs converge at the rate set by the next
    /// unretained cluster state, so requiring them without extra columns
    /// exhausts six 20-iteration calls at 5e-14 Ha (the reference ethr floor
    /// 1e-13 Ry / 2). Solving n + max(4, ceil(n/10)) = 24 columns and requiring
    /// only the first 20 converges in at most two calls to the dense answer.
    #[test]
    fn cut_cluster_converges_with_buffer_columns() {
        let (dim, n, cluster) = (400, 20, 16);
        let mut seed = 12345_u64;
        let mut rnd = || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((seed >> 11) as f64 / (1_u64 << 53) as f64) - 0.5
        };
        let mut h = Matrix::zeros(dim, dim).unwrap();
        let lo = n - cluster / 2;
        for i in 0..dim {
            h[(i, i)] = C::new(
                if i < lo {
                    -1. + i as f64 / lo as f64
                } else if i < lo + cluster {
                    1. + 1e-2 * (i - lo) as f64 / cluster as f64
                } else {
                    1.5 + 50. * ((i - lo - cluster) as f64 / dim as f64).powi(2)
                },
                0.,
            );
        }
        for i in 0..dim {
            for j in i + 1..dim {
                let z = C::new(rnd(), rnd()) * 0.015 / (1. + (j - i) as f64);
                h[(i, j)] += z;
                h[(j, i)] += z.conj();
            }
        }
        let (exact, _) = dft_numerics::linalg::hermitian_eigh(&h).unwrap();
        let width = n + 4_usize.max(n.div_ceil(10));
        let mut start = Matrix::zeros(dim, width).unwrap();
        for j in 0..width {
            for i in 0..dim {
                start[(i, j)] = C::new(f64::from(i == j), 0.) + C::new(rnd(), rnd()) * 0.2;
            }
        }
        let diagonal: Vec<_> = (0..dim).map(|i| h[(i, i)].re).collect();
        // Calls needed until the first n flags converge (None: six exhausted).
        let calls = |columns: usize| -> Option<(usize, f64)> {
            dft_mpi::simulate(1, |comm| {
                let mut op = Dense {
                    h: h.clone(),
                    rows: (0..dim).collect(),
                };
                let mut block = Matrix::zeros(dim, columns).unwrap();
                for j in 0..columns {
                    block.column_mut(j).copy_from_slice(start.column(j));
                }
                let o = Options::new(columns, 5e-14, n);
                for call in 1..=6 {
                    let s = solve(comm, &mut op, &block, &diagonal, o).unwrap();
                    if s.converged[..n].iter().all(|&v| v) {
                        let error = (0..n)
                            .map(|i| (s.eigenvalues[i] - exact[i]).abs())
                            .fold(0., f64::max);
                        return Some((call, error));
                    }
                    block = s.eigenvectors;
                }
                None
            })[0]
        };
        let unbuffered = calls(n);
        let (used, error) = calls(width).expect("buffer columns converge");
        assert!(used <= 2, "{used} calls");
        // Negative control: without buffer columns the cut cluster exhausts the
        // six calls on the native backend; with MKL rounding it may scrape
        // through, but only with clearly more calls than the buffered solve.
        if let Some((unbuffered_calls, _)) = unbuffered {
            assert!(
                unbuffered_calls > used + 2,
                "cut cluster without buffer: {unbuffered_calls} calls vs {used} buffered"
            );
        }
        assert!(error < 1e-11, "{error:e} Ha");
    }
    /// A k-like sequence: Hermitian operators of different dimensions (local
    /// row counts change between calls, as between k points), repeated twice
    /// (two "SCF maps"), both reduced placements, with fixtures that restart
    /// (capacity 2n < dim) and repair (n > dim/2). Returns a bit hash of every
    /// eigenpair, HX and decision count, and the per-iteration collectives.
    fn k_sequence(comm: &dft_mpi::SimRank, arm: &str) -> (u64, u64, usize, usize, [f64; 2]) {
        let config = dft_numerics::switches::RunConfig::from_pairs(&[("cdav_reuse", arm)]).unwrap();
        let _scope = config.enter();
        let mut hash = 0xcbf2_9ce4_8422_2325_u64;
        let mut mix = |bits: u64| {
            for b in bits.to_le_bytes() {
                hash = (hash ^ u64::from(b)).wrapping_mul(0x100_0000_01b3);
            }
        };
        let (mut collectives, mut iterations, mut repairs, mut restarts) = (0, 0, 0, 0);
        // Ordinary iterations (one H block, no repair/restart), per placement.
        let mut plain = [(0_u64, 0_u64); 2];
        for _map in 0..2 {
            for (k, &(dim, n)) in [(37, 6), (41, 6), (29, 6), (9, 5), (11, 4)]
                .iter()
                .enumerate()
            {
                let mut h = Matrix::zeros(dim, dim).unwrap();
                for i in 0..dim {
                    h[(i, i)] = C::new(i as f64 * 0.21 - 0.9 + 0.05 * k as f64, 0.);
                    for j in 0..i {
                        let z = C::new(
                            0.03 / (1. + (i - j) as f64) * ((i * 3 + j + k) as f64).sin(),
                            0.02 / (1. + (i - j) as f64) * ((i + 5 * j + k) as f64).cos(),
                        );
                        h[(i, j)] = z;
                        h[(j, i)] = z.conj();
                    }
                }
                let rows: Vec<_> = (0..dim)
                    .filter(|i| i % comm.size() == comm.rank())
                    .collect();
                let mut start = Matrix::zeros(rows.len(), n).unwrap();
                for j in 0..n {
                    for (i, &g) in rows.iter().enumerate() {
                        start[(i, j)] = C::new(f64::from(g == j), 0.)
                            + C::new(
                                ((g + 1) * (j + 3)) as f64 * 0.003,
                                ((g + 2) * (j + 1)) as f64 * -0.002,
                            );
                    }
                }
                let diagonal: Vec<_> = rows.iter().map(|&i| h[(i, i)].re).collect();
                let mut op = Dense { h, rows };
                for placement in [ReducedPlacement::Root, ReducedPlacement::Replicated] {
                    let mut o = Options::new(n, 1e-12, n);
                    o.max_iterations = 40;
                    o.return_hx = true;
                    o.profile = true;
                    o.placement = placement;
                    let s = solve(comm, &mut op, &start, &diagonal, o).unwrap();
                    for v in &s.eigenvalues {
                        mix(v.to_bits());
                    }
                    for z in s.eigenvectors.as_slice() {
                        mix(z.re.to_bits());
                        mix(z.im.to_bits());
                    }
                    for z in s.hx.as_ref().unwrap().as_slice() {
                        mix(z.re.to_bits());
                        mix(z.im.to_bits());
                    }
                    for c in &s.converged {
                        mix(u64::from(*c));
                    }
                    mix(s.record.h_columns as u64);
                    mix(s.record.iterations.len() as u64);
                    for r in &s.record.iterations {
                        let count = r.communication.reductions + r.communication.broadcasts;
                        collectives += count;
                        repairs += r.repairs;
                        restarts += r.restarts;
                        if r.repairs == 0 && r.restarts == 0 && r.h_calls == 1 {
                            let slot = usize::from(placement == ReducedPlacement::Replicated);
                            plain[slot].0 += count;
                            plain[slot].1 += 1;
                        }
                    }
                    iterations += s.record.iterations.len();
                }
            }
        }
        assert!(
            repairs > 0 && restarts > 0,
            "{repairs} repairs {restarts} restarts"
        );
        let ratio = |(c, i): (u64, u64)| c as f64 / i.max(1) as f64;
        (
            hash,
            collectives,
            iterations,
            repairs,
            [ratio(plain[0]), ratio(plain[1])],
        )
    }
    #[test]
    fn workspace_reuse_and_merged_agreements_are_bit_identical() {
        for ranks in [1, 3, 8] {
            let out = dft_mpi::simulate(ranks, |comm| {
                let off = k_sequence(comm, "off");
                let on = k_sequence(comm, "on");
                let again = k_sequence(comm, "on");
                (off, on, again)
            });
            for (off, on, again) in &out {
                assert_eq!(off.0, on.0, "ranks {ranks}");
                assert_eq!(on.0, again.0, "ranks {ranks}");
                assert_eq!(off.2, on.2);
            }
            let (off, on, _) = out[0];
            // The Dense test operator adds one reduction per H block.
            eprintln!(
                "GK3 ranks {ranks}: hash {:016x}, iterations {}, repairs {}, collectives per iteration (all) off {:.2} on {:.2}; ordinary iteration root off {:.2} on {:.2}, replicated off {:.2} on {:.2}",
                on.0,
                on.2,
                on.3,
                off.1 as f64 / off.2 as f64,
                on.1 as f64 / on.2 as f64,
                off.4[0],
                on.4[0],
                off.4[1],
                on.4[1]
            );
            if ranks > 1 {
                assert!(on.4[0] < off.4[0] && on.4[1] < off.4[1]);
            }
        }
    }
    /// A fault on ONE rank in a later H block (non-finite output or an
    /// error) must become an error on every rank in both GK3 arms and both
    /// placements, without a hang; the kept workspace then serves a clean
    /// solve bit-identically.
    #[test]
    fn complex_reference_agrees_one_rank_fault_mid_solve() {
        struct Late {
            inner: Dense,
            calls: usize,
            fault_call: usize,
            fault_rank: usize,
            nan: bool,
        }
        impl DistributedComplexOperator for Late {
            fn global_len(&self) -> usize {
                self.inner.global_len()
            }
            fn local_rows(&self) -> &[usize] {
                self.inner.local_rows()
            }
            fn apply_block<K: Collective>(
                &mut self,
                comm: &K,
                x: &Matrix,
                out: &mut Matrix,
            ) -> Result<(), SolverError> {
                self.inner.apply_block(comm, x, out)?;
                self.calls += 1;
                if self.calls == self.fault_call && comm.rank() == self.fault_rank {
                    if self.nan {
                        if let Some(z) = out.as_mut_slice().first_mut() {
                            *z = C::new(f64::NAN, 0.);
                        }
                    } else {
                        return Err(error("injected local failure"));
                    }
                }
                Ok(())
            }
        }
        let dim = 23;
        let n = 4;
        let mut h = Matrix::zeros(dim, dim).unwrap();
        for i in 0..dim {
            h[(i, i)] = C::new(i as f64 * 0.17 - 0.5, 0.);
            if i > 0 {
                h[(i - 1, i)] = C::new(0.05, -0.03);
                h[(i, i - 1)] = h[(i - 1, i)].conj();
            }
        }
        for ranks in [2, 3, 8] {
            for arm in ["off", "on"] {
                for placement in [ReducedPlacement::Root, ReducedPlacement::Replicated] {
                    for (fault_call, nan) in [(1, true), (3, true), (3, false), (5, true)] {
                        dft_mpi::simulate(ranks, |comm| {
                            let config = dft_numerics::switches::RunConfig::from_pairs(&[(
                                "cdav_reuse",
                                arm,
                            )])
                            .unwrap();
                            let _scope = config.enter();
                            let rows: Vec<_> =
                                (0..dim).filter(|i| i % ranks == comm.rank()).collect();
                            let mut start = Matrix::zeros(rows.len(), n).unwrap();
                            for j in 0..n {
                                for (i, &g) in rows.iter().enumerate() {
                                    start[(i, j)] = C::new(f64::from(g == j), 0.)
                                        + C::new(0.01 * (g + j) as f64, 0.002 * g as f64);
                                }
                            }
                            let diagonal: Vec<_> = rows.iter().map(|&i| h[(i, i)].re).collect();
                            let mut o = Options::new(n, 1e-12, n);
                            o.max_iterations = 30;
                            o.placement = placement;
                            let mut clean = Dense {
                                h: h.clone(),
                                rows: rows.clone(),
                            };
                            let want = solve(comm, &mut clean, &start, &diagonal, o).unwrap();
                            let mut op = Late {
                                inner: Dense { h: h.clone(), rows },
                                calls: 0,
                                fault_call,
                                fault_rank: ranks - 1,
                                nan,
                            };
                            assert!(
                                solve(comm, &mut op, &start, &diagonal, o).is_err(),
                                "rank {} ranks {ranks} arm {arm} call {fault_call}",
                                comm.rank()
                            );
                            let again = solve(comm, &mut clean, &start, &diagonal, o).unwrap();
                            assert_eq!(want.eigenvalues, again.eigenvalues);
                            assert_eq!(want.eigenvectors, again.eigenvectors);
                        });
                    }
                }
            }
        }
    }
    #[test]
    fn complex_reference_agrees_nonfinite_operator_output() {
        struct Fault {
            rows: Vec<usize>,
        }
        impl DistributedComplexOperator for Fault {
            fn global_len(&self) -> usize {
                1
            }
            fn local_rows(&self) -> &[usize] {
                &self.rows
            }
            fn apply_block<K: Collective>(
                &mut self,
                _comm: &K,
                _x: &Matrix,
                out: &mut Matrix,
            ) -> Result<(), SolverError> {
                out.as_mut_slice().fill(C::new(f64::NAN, 0.));
                Ok(())
            }
        }
        for ranks in [1, 3, 8] {
            dft_mpi::simulate(ranks, |comm| {
                let rows = if comm.rank() == 0 {
                    vec![0]
                } else {
                    Vec::new()
                };
                let mut start = Matrix::zeros(rows.len(), 1).unwrap();
                start.as_mut_slice().fill(C::new(1., 0.));
                let diagonal = vec![1.; rows.len()];
                let mut op = Fault { rows };
                assert!(solve(comm, &mut op, &start, &diagonal, Options::new(1, 1e-8, 1)).is_err());
            });
        }
    }
}
