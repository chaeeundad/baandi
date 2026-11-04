//! Fixed backend decomposition, independent of the requested thread budget.
use super::*;
use crate::pool::WorkerPool;
use crate::pool::analysis::{Piece, Span};
use std::sync::Arc;

// A cache-sized plane-wave panel, independent of machine and worker count.
#[cfg(not(feature = "mkl"))]
const PANEL: usize = 256;

// Projection-only partition: 2969 rows yield 47 independent pieces at every T.
const PROJECTION_PANEL: usize = 64;
// The candidate failed the multi-system iteration gate. Preserve it for
// controlled experiments without changing the accepted default or MKL work.
#[allow(unexpected_cfgs)]
const BLOCK_QR: bool = cfg!(feature = "block-qr");

#[cfg(not(feature = "mkl"))]
struct Tile {
    c: ComplexMatrix,
}

/// Reusable storage for tall-skinny products and their Gram matrices.
/// Numerical scratch is allocated at construction. SCF callers share their
/// owned worker pool with Hamiltonian and projection phases.
pub struct GemmWorkspace {
    #[cfg(not(feature = "mkl"))]
    tiles: Vec<Tile>,
    long: usize,
    short: usize,
    threads: usize,
    pool: Arc<WorkerPool>,
}
impl GemmWorkspace {
    pub fn new(long: usize, short: usize, threads: usize) -> Result<Self, LinalgError> {
        let pool = WorkerPool::new(threads).map_err(|_| LinalgError::Decomposition)?;
        Self::with_pool(long, short, threads, pool)
    }

    pub fn with_pool(
        long: usize,
        short: usize,
        threads: usize,
        pool: Arc<WorkerPool>,
    ) -> Result<Self, LinalgError> {
        #[cfg(not(feature = "mkl"))]
        let tiles = {
            let size = short;
            let mut tiles = Vec::with_capacity(
                long.div_ceil(if BLOCK_QR || cfg!(test) {
                    PROJECTION_PANEL
                } else {
                    PANEL
                })
                .max(1),
            );
            for _ in 0..long
                .div_ceil(if BLOCK_QR || cfg!(test) {
                    PROJECTION_PANEL
                } else {
                    PANEL
                })
                .max(1)
            {
                tiles.push(Tile {
                    c: ComplexMatrix::zeros(size, short)?,
                });
            }
            tiles
        };
        Ok(Self {
            #[cfg(not(feature = "mkl"))]
            tiles,
            long,
            short,
            threads: threads.max(1),
            pool,
        })
    }

    /// Allocating convenience boundary for legacy Davidson result ownership.
    /// Use `multiply_into` when the caller already owns reusable output storage.
    pub fn multiply(
        &mut self,
        a: &ComplexMatrix,
        oa: MatrixOp,
        b: &ComplexMatrix,
        ob: MatrixOp,
    ) -> Result<ComplexMatrix, LinalgError> {
        let (m, _) = shape(a, oa);
        let (_, n) = shape(b, ob);
        let mut out = ComplexMatrix::zeros(m, n)?;
        self.multiply_into(a, oa, b, ob, &mut out)?;
        Ok(out)
    }

    /// Rotate a basis and its Hamiltonian image in one scheduling region.
    /// Each fixed panel still performs the same two sequential backend products.
    /// Result allocation matches the two legacy allocating multiply calls.
    pub fn rotate_pair(
        &mut self,
        a: &ComplexMatrix,
        h: &ComplexMatrix,
        b: &ComplexMatrix,
    ) -> Result<(ComplexMatrix, ComplexMatrix), LinalgError> {
        if a.rows != h.rows || a.cols != h.cols {
            return Err(LinalgError::Dimension);
        }
        #[cfg(not(feature = "mkl"))]
        if self.threads > 1 && a.rows >= a.cols {
            let (m, k, n) = (a.rows, a.cols, b.cols);
            if k != b.rows || m > self.long || k.max(n) > self.short {
                return Err(LinalgError::Dimension);
            }
            a.validate_finite()?;
            h.validate_finite()?;
            b.validate_finite()?;
            let mut out = ComplexMatrix::zeros(m, n)?;
            let mut hout = ComplexMatrix::zeros(m, n)?;
            let job = Job {
                a,
                oa: MatrixOp::Identity,
                b,
                ob: MatrixOp::Identity,
                m,
                k,
                n,
                reduce: false,
                panel: PANEL,
            };
            let hjob = Job { a: h, ..job };
            let count = m.div_ceil(PANEL).max(1);
            // Keep the existing scheduling threshold for each logical product.
            let threads = if m.saturating_mul(n).saturating_mul(k) < 262_144 || count < 4 {
                1
            } else {
                self.threads.min(count)
            };
            crate::pool::analysis::logical_gemm(m, n, k);
            crate::pool::analysis::logical_gemm(m, n, k);
            let output = faer::MatMut::from_column_major_slice_mut(&mut out.data, m, n);
            let houtput = faer::MatMut::from_column_major_slice_mut(&mut hout.data, m, n);
            let _region = Span::region("rotation_pair", threads);
            if threads <= 1 {
                rotate_pair(output, houtput, 0, 1, &job, &hjob);
            } else {
                self.pool
                    .install(|| rotate_pair(output, houtput, 0, threads, &job, &hjob));
            }
            out.validate_finite()?;
            hout.validate_finite()?;
            return Ok((out, hout));
        }
        Ok((
            self.multiply(a, MatrixOp::Identity, b, MatrixOp::Identity)?,
            self.multiply(h, MatrixOp::Identity, b, MatrixOp::Identity)?,
        ))
    }

    /// Pure Rust uses fixed 256-row panels for rotations, or 256-plane-wave partial
    /// products for Gram matrices. Partial products are summed in panel order.
    /// Identical tile shapes and backend calls are used even with one worker.
    /// MKL uses one sequential call for the entire logical product at every
    /// budget, avoiding packed operand copies and small BLAS calls.
    pub fn multiply_into(
        &mut self,
        a: &ComplexMatrix,
        oa: MatrixOp,
        b: &ComplexMatrix,
        ob: MatrixOp,
        out: &mut ComplexMatrix,
    ) -> Result<(), LinalgError> {
        a.validate_finite()?;
        b.validate_finite()?;
        self.multiply_into_finite_inputs(a, oa, b, ob, out)
    }

    // Both operands must have been checked since their last mutation. Private
    // so only the projection workspace can reuse its validation of the basis.
    fn multiply_into_finite_inputs(
        &mut self,
        a: &ComplexMatrix,
        oa: MatrixOp,
        b: &ComplexMatrix,
        ob: MatrixOp,
        out: &mut ComplexMatrix,
    ) -> Result<(), LinalgError> {
        self.multiply_panels(a, oa, b, ob, out, 256)
    }

    #[allow(clippy::too_many_arguments)] // GEMM operands plus a fixed data partition.
    fn multiply_panels(
        &mut self,
        a: &ComplexMatrix,
        oa: MatrixOp,
        b: &ComplexMatrix,
        ob: MatrixOp,
        out: &mut ComplexMatrix,
        _panel: usize,
    ) -> Result<(), LinalgError> {
        #[cfg(not(feature = "mkl"))]
        let panel = _panel;
        let (m, k) = shape(a, oa);
        let (l, n) = shape(b, ob);
        if k != l
            || out.rows != m
            || out.cols != n
            || n > self.short
            || m.max(k) > self.long
            || m.min(k) > self.short
        {
            return Err(LinalgError::Dimension);
        }
        crate::pool::analysis::logical_gemm(m, n, k);
        #[cfg(feature = "mkl")]
        {
            crate::pool::analysis::executed_gemm(m, n, k);
            gemm_into_finite_inputs(a, oa, b, ob, out)
        }
        #[cfg(not(feature = "mkl"))]
        {
            let reduce = k > m;
            let count = if reduce { k } else { m }.div_ceil(panel).max(1);
            let job = Job {
                a,
                oa,
                b,
                ob,
                m,
                k,
                n,
                reduce,
                panel,
            };
            // Fixed 2,097,152-real-flop cutoff (8 flops per complex multiply-add).
            // Require four panels to avoid a single full panel plus a tiny tail.
            // This affects scheduling only, never tile shapes or arithmetic order.
            let products = m.saturating_mul(n).saturating_mul(k);
            let threads = if products < 262_144 || count < 4 {
                1
            } else {
                self.threads.min(count)
            };
            let _region = Span::region(if reduce { "gram" } else { "rotation" }, threads);
            if !reduce {
                let output = faer::MatMut::from_column_major_slice_mut(&mut out.data, m, n);
                if threads <= 1 {
                    rotate(output, 0, 1, &job);
                } else {
                    self.pool.install(|| rotate(output, 0, threads, &job));
                }
                return out.validate_finite();
            }
            if threads <= 1 {
                run(&mut self.tiles[..count], 0, 1, &job)?;
            } else {
                self.pool
                    .install(|| run(&mut self.tiles[..count], 0, threads, &job))?;
            }
            let _assembly = Span::serial("gram_assembly");
            out.data.fill(Complex64::default());
            for tile in &self.tiles[..count] {
                for (y, x) in out.data.iter_mut().zip(&tile.c.data) {
                    *y += x;
                }
            }
            out.validate_finite()
        }
    }
}
fn shape(a: &ComplexMatrix, op: MatrixOp) -> (usize, usize) {
    match op {
        MatrixOp::Identity => (a.rows, a.cols),
        MatrixOp::Adjoint => (a.cols, a.rows),
    }
}
#[cfg(test)]
fn value(a: &ComplexMatrix, op: MatrixOp, i: usize, j: usize) -> Complex64 {
    match op {
        MatrixOp::Identity => a[(i, j)],
        MatrixOp::Adjoint => a[(j, i)].conj(),
    }
}
fn reshape(a: &mut ComplexMatrix, rows: usize, cols: usize) {
    assert!(rows * cols <= a.data.capacity());
    a.data.resize(rows * cols, Complex64::default());
    a.rows = rows;
    a.cols = cols;
}
#[cfg(not(feature = "mkl"))]
struct Job<'a> {
    a: &'a ComplexMatrix,
    oa: MatrixOp,
    b: &'a ComplexMatrix,
    ob: MatrixOp,
    m: usize,
    k: usize,
    n: usize,
    reduce: bool,
    panel: usize,
}
#[cfg(not(feature = "mkl"))]
fn run(tiles: &mut [Tile], first: usize, threads: usize, job: &Job<'_>) -> Result<(), LinalgError> {
    if threads <= 1 || tiles.len() == 1 {
        for (offset, tile) in tiles.iter_mut().enumerate() {
            let _piece = Piece::new();
            let start = (first + offset) * job.panel;
            let (rows, inner) = if job.reduce {
                (job.m, (job.k - start).min(job.panel))
            } else {
                ((job.m - start).min(job.panel), job.k)
            };
            reshape(&mut tile.c, rows, job.n);
            #[cfg(not(feature = "mkl"))]
            borrowed_product(
                job,
                start,
                rows,
                inner,
                faer::MatMut::from_column_major_slice_mut(&mut tile.c.data, rows, job.n),
            );
        }
        return Ok(());
    }
    let middle = tiles.len() / 2;
    let (left, right) = tiles.split_at_mut(middle);
    let (left, right) = WorkerPool::join(
        || run(left, first, threads, job),
        || run(right, first + middle, threads, job),
    );
    left?;
    right
}

// Each output row belongs to exactly one fixed panel. Strided views write
// directly into the caller's matrix; there is no rotation packing/assembly.
#[cfg(not(feature = "mkl"))]
fn rotate(mut out: faer::MatMut<'_, Complex64>, first: usize, threads: usize, job: &Job<'_>) {
    if threads <= 1 || out.nrows() <= job.panel {
        let mut start = first;
        while out.nrows() > 0 {
            let _piece = Piece::new();
            let rows = out.nrows().min(job.panel);
            let (panel, rest) = out.split_at_row_mut(rows);
            borrowed_product(job, start, rows, job.k, panel);
            start += rows;
            out = rest;
        }
        return;
    }
    let panels = out.nrows().div_ceil(job.panel);
    let middle = panels / 2 * job.panel;
    let (left, right) = out.split_at_row_mut(middle);
    WorkerPool::join(
        || rotate(left, first, threads, job),
        || rotate(right, first + middle, threads, job),
    );
}

#[cfg(not(feature = "mkl"))]
fn rotate_pair(
    mut out: faer::MatMut<'_, Complex64>,
    mut hout: faer::MatMut<'_, Complex64>,
    first: usize,
    threads: usize,
    job: &Job<'_>,
    hjob: &Job<'_>,
) {
    if threads <= 1 || out.nrows() <= job.panel {
        let mut start = first;
        while out.nrows() > 0 {
            let _piece = Piece::new();
            let rows = out.nrows().min(job.panel);
            let (panel, rest) = out.split_at_row_mut(rows);
            let (hpanel, hrest) = hout.split_at_row_mut(rows);
            borrowed_product(job, start, rows, job.k, panel);
            borrowed_product(hjob, start, rows, hjob.k, hpanel);
            start += rows;
            out = rest;
            hout = hrest;
        }
        return;
    }
    let panels = out.nrows().div_ceil(job.panel);
    let middle = panels / 2 * job.panel;
    let (left, right) = out.split_at_row_mut(middle);
    let (hleft, hright) = hout.split_at_row_mut(middle);
    WorkerPool::join(
        || rotate_pair(left, hleft, first, threads, job, hjob),
        || rotate_pair(right, hright, first + middle, threads, job, hjob),
    );
}

#[cfg(not(feature = "mkl"))]
fn borrowed_product(
    job: &Job<'_>,
    start: usize,
    rows: usize,
    inner: usize,
    mut c: faer::MatMut<'_, Complex64>,
) {
    crate::pool::analysis::executed_gemm(rows, job.n, inner);
    let row = if job.reduce { 0 } else { start };
    let reduction = if job.reduce { start } else { 0 };
    let a = faer::MatRef::from_column_major_slice(&job.a.data, job.a.rows, job.a.cols);
    let b = faer::MatRef::from_column_major_slice(&job.b.data, job.b.rows, job.b.cols);
    let a = match job.oa {
        MatrixOp::Identity => a.submatrix(row, reduction, rows, inner),
        MatrixOp::Adjoint => a.submatrix(reduction, row, inner, rows),
    };
    let b = match job.ob {
        MatrixOp::Identity => b.submatrix(reduction, 0, inner, job.n),
        MatrixOp::Adjoint => b.submatrix(0, reduction, job.n, inner),
    };
    let one = Complex64::new(1.0, 0.0);
    match (job.oa, job.ob) {
        (MatrixOp::Identity, MatrixOp::Identity) => {
            matmul(&mut c, Accum::Replace, a, b, one, Par::Seq)
        }
        (MatrixOp::Adjoint, MatrixOp::Identity) => {
            matmul(&mut c, Accum::Replace, a.adjoint(), b, one, Par::Seq)
        }
        (MatrixOp::Identity, MatrixOp::Adjoint) => {
            matmul(&mut c, Accum::Replace, a, b.adjoint(), one, Par::Seq)
        }
        (MatrixOp::Adjoint, MatrixOp::Adjoint) => matmul(
            &mut c,
            Accum::Replace,
            a.adjoint(),
            b.adjoint(),
            one,
            Par::Seq,
        ),
    }
}

/// Block CGS2 against the established basis, interleaved with CholeskyQR2
/// inside the incoming block. Ill-conditioned blocks retain two-pass MGS and
/// its rank decisions relative to original norms.
/// All projection buffers are reserved once, outside Davidson iterations.
pub struct ProjectionWorkspace {
    normalized: ComplexMatrix,
    coefficients: ComplexMatrix,
    projected: ComplexMatrix,
    column: Vec<Complex64>,
    indices: Vec<usize>,
    gram: ComplexMatrix,
    inverse: ComplexMatrix,
    backup: ComplexMatrix,
    cholesky_scratch: faer::dyn_stack::MemBuffer,
}
impl ProjectionWorkspace {
    pub fn new(rows: usize, capacity: usize) -> Result<Self, LinalgError> {
        let qr_capacity = if BLOCK_QR || cfg!(test) { capacity } else { 0 };
        Ok(Self {
            normalized: ComplexMatrix::zeros(rows, capacity)?,
            coefficients: ComplexMatrix::zeros(capacity, capacity)?,
            projected: ComplexMatrix::zeros(rows, capacity)?,
            column: vec![Complex64::default(); rows],
            indices: Vec::with_capacity(capacity),
            gram: ComplexMatrix::zeros(qr_capacity, qr_capacity)?,
            inverse: ComplexMatrix::zeros(qr_capacity, qr_capacity)?,
            backup: ComplexMatrix::zeros(rows, qr_capacity)?,
            cholesky_scratch: faer::dyn_stack::MemBuffer::new(
                faer::linalg::cholesky::llt::factor::cholesky_in_place_scratch::<Complex64>(
                    qr_capacity,
                    Par::Seq,
                    Default::default(),
                ),
            ),
        })
    }
    // Gram conditioning is a suitability check, not a rank decision. If the
    // normal equations lose too much accuracy, restore the input and use the
    // existing relative-norm rank policy. All buffers and backend scratch are
    // reserved at construction, including the rollback copy.
    fn try_cholesky_qr(&mut self, gemm: &mut GemmWorkspace) -> Result<bool, LinalgError> {
        let (rows, n) = (self.normalized.rows, self.normalized.cols);
        if n < 2 {
            return Ok(false);
        }
        reshape(&mut self.backup, rows, n);
        crate::pool::analysis::copied(self.normalized.data.len() * size_of::<Complex64>());
        self.backup.data.copy_from_slice(&self.normalized.data);
        reshape(&mut self.gram, n, n);
        reshape(&mut self.inverse, n, n);
        reshape(&mut self.projected, rows, n);
        let _phase = Span::new("block_qr");
        for _ in 0..2 {
            gemm.multiply_panels(
                &self.normalized,
                MatrixOp::Adjoint,
                &self.normalized,
                MatrixOp::Identity,
                &mut self.gram,
                PROJECTION_PANEL,
            )?;
            let factor = Span::serial("block_qr_factor");
            let scale = (0..n).map(|j| self.gram[(j, j)].re).fold(0.0, f64::max);
            let ok = faer::linalg::cholesky::llt::factor::cholesky_in_place(
                faer::MatMut::from_column_major_slice_mut(&mut self.gram.data, n, n),
                Default::default(),
                Par::Seq,
                faer::dyn_stack::MemStack::new(&mut self.cholesky_scratch),
                Default::default(),
            )
            .is_ok();
            if !ok
                || (0..n).any(|j| {
                    let pivot = self.gram[(j, j)].re;
                    !pivot.is_finite()
                        || pivot <= DEPENDENCE_TOLERANCE
                        || pivot * pivot <= f64::EPSILON.sqrt() * scale
                })
            {
                crate::pool::analysis::copied(self.normalized.data.len() * size_of::<Complex64>());
                self.normalized.data.copy_from_slice(&self.backup.data);
                return Ok(false);
            }
            // Invert only the small triangular factor through the library;
            // the tall update X L^{-H} is a row-partitioned GEMM.
            self.inverse.data.fill(Complex64::default());
            for j in 0..n {
                self.inverse[(j, j)] = Complex64::new(1.0, 0.0);
            }
            faer::linalg::triangular_solve::solve_lower_triangular_in_place(
                faer::MatRef::from_column_major_slice(&self.gram.data, n, n),
                faer::MatMut::from_column_major_slice_mut(&mut self.inverse.data, n, n),
                Par::Seq,
            );
            drop(factor);
            gemm.multiply_panels(
                &self.normalized,
                MatrixOp::Identity,
                &self.inverse,
                MatrixOp::Adjoint,
                &mut self.projected,
                PROJECTION_PANEL,
            )?;
            crate::pool::analysis::copied(self.normalized.data.len() * size_of::<Complex64>());
            self.normalized.data.copy_from_slice(&self.projected.data);
        }
        Ok(true)
    }

    pub fn orthonormalize(
        &mut self,
        block: &ComplexMatrix,
        against: Option<&ComplexMatrix>,
        gemm: &mut GemmWorkspace,
    ) -> Result<Orthonormalized, LinalgError> {
        let _phase = Span::new("projection");
        let prepare = Span::serial("projection_prepare");
        block.validate_finite()?;
        if block.rows != self.column.len() || block.data.len() > self.normalized.data.capacity() {
            return Err(LinalgError::Dimension);
        }
        if let Some(q) = against {
            q.validate_finite()?;
            if q.rows != block.rows || q.cols * block.cols > self.coefficients.data.capacity() {
                return Err(LinalgError::Dimension);
            }
        }
        reshape(&mut self.normalized, block.rows, block.cols);
        for j in 0..block.cols {
            let length = norm(block.column(j));
            if !length.is_finite() {
                return Err(LinalgError::NonFinite);
            }
            for (out, x) in self
                .normalized
                .column_mut(j)
                .iter_mut()
                .zip(block.column(j))
            {
                *out = if length == 0.0 {
                    Complex64::default()
                } else {
                    *x / length
                };
            }
        }
        drop(prepare);
        // Interleave projection and normalization: a nearly dependent incoming
        // block can amplify roundoff during its internal orthogonalization.
        // Reprojecting the normalized block prevents that amplification from
        // leaving a component in the established basis.
        let against = against.filter(|q| q.cols > 0);
        let mut data = Vec::with_capacity(block.data.len());
        let mut dropped = Vec::with_capacity(block.cols);
        self.indices.clear();
        self.indices.extend(0..block.cols);
        let mut count = block.cols;
        for pass in 0..if against.is_some() { 2 } else { 1 } {
            if pass > 0 {
                reshape(&mut self.normalized, block.rows, count);
                crate::pool::analysis::copied(self.normalized.data.len() * size_of::<Complex64>());
                self.normalized.data.copy_from_slice(&data);
            }
            if let Some(q) = against {
                reshape(&mut self.coefficients, q.cols, count);
                reshape(&mut self.projected, block.rows, count);
                // q was checked at entry and is immutable. Normalization (or
                // the preceding MGS pass) changed the incoming block; check it
                // once. The first product checks coefficients on output.
                self.normalized.validate_finite()?;
                gemm.multiply_panels(
                    q,
                    MatrixOp::Adjoint,
                    &self.normalized,
                    MatrixOp::Identity,
                    &mut self.coefficients,
                    if BLOCK_QR && block.rows >= 2048 {
                        PROJECTION_PANEL
                    } else {
                        256
                    },
                )?;
                gemm.multiply_panels(
                    q,
                    MatrixOp::Identity,
                    &self.coefficients,
                    MatrixOp::Identity,
                    &mut self.projected,
                    if BLOCK_QR && block.rows >= 2048 {
                        PROJECTION_PANEL
                    } else {
                        256
                    },
                )?;
                let _update = Span::serial("projection_update");
                for (y, p) in self.normalized.data.iter_mut().zip(&self.projected.data) {
                    *y -= p;
                }
            }
            data.clear();
            if BLOCK_QR && block.rows >= 2048 && self.try_cholesky_qr(gemm)? {
                crate::pool::analysis::copied(self.normalized.data.len() * size_of::<Complex64>());
                data.extend_from_slice(&self.normalized.data);
                continue;
            }
            count = 0;
            let mut prepared = None;
            for j in 0..self.normalized.cols {
                let second = Span::serial("mgs_second_pass");
                let remaining = if let Some(remaining) = prepared.take() {
                    remaining
                } else {
                    crate::pool::analysis::copied(block.rows * size_of::<Complex64>());
                    self.column.copy_from_slice(self.normalized.column(j));
                    // Tall blocks receive the first pass as earlier vectors are
                    // finalized. Small blocks retain the original two-pass loop.
                    for _ in 0..if block.rows >= 2048 { 1 } else { 2 } {
                        for k in 0..count {
                            subtract_projection(
                                &mut self.column,
                                &data[k * block.rows..(k + 1) * block.rows],
                            );
                        }
                    }
                    norm(&self.column)
                };
                drop(second);
                if !remaining.is_finite() {
                    return Err(LinalgError::NonFinite);
                }
                if remaining <= DEPENDENCE_TOLERANCE {
                    dropped.push(self.indices[j]);
                } else {
                    data.extend(self.column.iter().map(|x| x / remaining));
                    self.indices[count] = self.indices[j];
                    count += 1;
                    if block.rows < 2048 {
                        continue;
                    }
                    let q = &data[(count - 1) * block.rows..count * block.rows];
                    let future = &mut self.normalized.data[(j + 1) * block.rows..];
                    let threads = if future.len() >= 16_384 {
                        gemm.threads
                    } else {
                        1
                    };
                    let _region = Span::region("mgs_first_pass", threads);
                    if threads > 1 {
                        let (next, rest) = future.split_at_mut(block.rows);
                        let column = &mut self.column;
                        // The next column's dependent pass overlaps independent
                        // first-pass work on later columns. No basis entry is
                        // appended until both branches finish. Every column sees
                        // exactly the historical ascending projection sequence.
                        let (remaining, ()) = gemm.pool.install(|| {
                            WorkerPool::join(
                                || {
                                    crate::pool::analysis::copied(
                                        block.rows * size_of::<Complex64>(),
                                    );
                                    column.copy_from_slice(next);
                                    {
                                        let _piece = Piece::new();
                                        subtract_projection(column, q);
                                    }
                                    let _second = Span::serial("mgs_second_pass");
                                    for basis in data.chunks_exact(block.rows) {
                                        subtract_projection(column, basis);
                                    }
                                    norm(column)
                                },
                                || project_columns(rest, block.rows, q, threads),
                            )
                        });
                        prepared = Some(remaining);
                    } else {
                        project_columns(future, block.rows, q, 1);
                    }
                }
            }
            self.indices.truncate(count);
        }
        dropped.sort_unstable();
        Ok(Orthonormalized {
            vectors: ComplexMatrix::new(block.rows, count, data)?,
            dropped,
        })
    }
}

// Fixed one-column pieces, recursively stealable; each dot/update uses the
// unchanged scalar kernel and ascending basis order.
fn project_columns(columns: &mut [Complex64], rows: usize, q: &[Complex64], threads: usize) {
    let count = columns.len() / rows;
    if threads <= 1 || count <= 1 {
        for column in columns.chunks_exact_mut(rows) {
            let _piece = Piece::new();
            subtract_projection(column, q);
        }
    } else {
        let (left, right) = columns.split_at_mut(count / 2 * rows);
        WorkerPool::join(
            || project_columns(left, rows, q, threads),
            || project_columns(right, rows, q, threads),
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn matrix(rows: usize, cols: usize) -> ComplexMatrix {
        ComplexMatrix::new(
            rows,
            cols,
            (0..rows * cols)
                .map(|i| {
                    Complex64::new(
                        ((i * i * 17 + cols * 31 + 3) as f64).sin(),
                        ((i * i * 7 + cols * 13 + 1) as f64).cos(),
                    )
                })
                .collect(),
        )
        .unwrap()
    }
    fn bits(a: &ComplexMatrix, b: &ComplexMatrix) {
        assert_eq!((a.rows, a.cols), (b.rows, b.cols));
        for (a, b) in a.data.iter().zip(&b.data) {
            assert_eq!(
                (a.re.to_bits(), a.im.to_bits()),
                (b.re.to_bits(), b.im.to_bits())
            );
        }
    }
    #[test]
    fn products_all_adjoint_shapes_partial_tiles_and_budgets() {
        use MatrixOp::{Adjoint, Identity};
        for (m, k, n) in [
            (513, 11, 7),
            (11, 513, 7),
            (513, 23, 19),
            (23, 513, 19),
            (41, 1025, 33),
            (1025, 41, 33),
            (7, 11, 3),
            (0, 0, 0),
            (7, 0, 3),
        ] {
            for oa in [Identity, Adjoint] {
                for ob in [Identity, Adjoint] {
                    let a = if matches!(oa, Identity) {
                        matrix(m, k)
                    } else {
                        matrix(k, m)
                    };
                    let b = if matches!(ob, Identity) {
                        matrix(k, n)
                    } else {
                        matrix(n, k)
                    };
                    let mut reference = ComplexMatrix::zeros(m, n).unwrap();
                    // Independent scalar complex multiplication oracle.
                    for j in 0..n {
                        for i in 0..m {
                            for l in 0..k {
                                reference[(i, j)] += value(&a, oa, i, l) * value(&b, ob, l, j);
                            }
                        }
                    }
                    let mut first = None;
                    for threads in [0, 1, 2, 3, 4, 8, 16] {
                        let mut work =
                            GemmWorkspace::new(m.max(k), m.min(k).max(n), threads).unwrap();
                        #[cfg(not(feature = "mkl"))]
                        let pointers: Vec<_> =
                            work.tiles.iter().map(|t| t.c.data.as_ptr()).collect();
                        let result = work.multiply(&a, oa, &b, ob).unwrap();
                        for (a, b) in result.data.iter().zip(&reference.data) {
                            assert!((*a - *b).norm() < 2e-12 * (1.0 + b.norm()));
                        }
                        if let Some(prev) = &first {
                            bits(prev, &result)
                        } else {
                            first = Some(result)
                        }
                        let mut repeat = ComplexMatrix::zeros(m, n).unwrap();
                        work.multiply_into(&a, oa, &b, ob, &mut repeat).unwrap();
                        bits(first.as_ref().unwrap(), &repeat);
                        #[cfg(not(feature = "mkl"))]
                        for (t, p) in work.tiles.iter().zip(pointers) {
                            assert_eq!(t.c.data.as_ptr(), p);
                        }
                    }
                }
            }
        }
    }
    #[test]
    fn paired_rotation_matches_separate_products() {
        for (m, k, n) in [(0, 0, 0), (7, 11, 3), (513, 11, 7), (4097, 41, 33)] {
            let a = matrix(m, k);
            let mut h = matrix(m, k);
            for z in &mut h.data {
                *z *= Complex64::new(0.7, -0.2);
            }
            let b = matrix(k, n);
            let mut first = None;
            for threads in [1, 2, 3, 4, 8, 12] {
                let mut work = GemmWorkspace::new(m.max(k), k.max(n), threads).unwrap();
                let expected = work
                    .multiply(&a, MatrixOp::Identity, &b, MatrixOp::Identity)
                    .unwrap();
                let hexpected = work
                    .multiply(&h, MatrixOp::Identity, &b, MatrixOp::Identity)
                    .unwrap();
                let (actual, hactual) = work.rotate_pair(&a, &h, &b).unwrap();
                bits(&expected, &actual);
                bits(&hexpected, &hactual);
                if let Some((a, h)) = &first {
                    bits(a, &actual);
                    bits(h, &hactual);
                } else {
                    first = Some((actual, hactual));
                }
            }
        }
    }

    #[test]
    fn right_looking_first_pass_matches_scalar_two_pass_mgs_bits() {
        for (rows, width) in [(257, 9), (2053, 17)] {
            let mut block = matrix(rows, width);
            block.column_mut(1).fill(Complex64::default());
            for i in 0..rows {
                block[(i, 3)] = block[(i, 2)];
            }
            for z in block.column_mut(4) {
                *z *= 1e-200;
            }
            for z in block.column_mut(5) {
                *z *= 1e200;
            }
            let expected = orthonormalize(&block, None).unwrap();
            for threads in [1, 2, 4, 8] {
                let mut gemm = GemmWorkspace::new(rows, width, threads).unwrap();
                let mut projection = ProjectionWorkspace::new(rows, width).unwrap();
                let actual = projection.orthonormalize(&block, None, &mut gemm).unwrap();
                assert_eq!(expected.dropped, actual.dropped);
                bits(&expected.vectors, &actual.vectors);
            }
        }
    }

    fn check_block_qr(rows: usize, width: usize) {
        let block = matrix(rows, width);
        let mut first = None;
        for threads in [1, 2, 4, 8, 12] {
            let mut work = GemmWorkspace::new(rows, width, threads).unwrap();
            let mut projection = ProjectionWorkspace::new(rows, width).unwrap();
            reshape(&mut projection.normalized, rows, width);
            projection.normalized.data.copy_from_slice(&block.data);
            assert!(projection.try_cholesky_qr(&mut work).unwrap());
            let q = &projection.normalized;
            let gram = gemm(q, MatrixOp::Adjoint, q, MatrixOp::Identity).unwrap();
            let mut defect = 0.0_f64;
            for j in 0..width {
                for i in 0..width {
                    defect = defect.max(
                        (gram[(i, j)] - Complex64::new(if i == j { 1.0 } else { 0.0 }, 0.0)).norm(),
                    );
                }
            }
            assert!(
                defect < 1e-12,
                "dimensionless orthogonality defect {defect}"
            );
            // A positive triangular QR factor preserves the ordered column
            // span: compare against independently accumulated scalar MGS.
            let scalar = orthonormalize(&block, None).unwrap().vectors;
            assert!(
                q.data
                    .iter()
                    .zip(&scalar.data)
                    .all(|(a, b)| (*a - *b).norm() < 1e-12)
            );
            if let Some(previous) = &first {
                bits(previous, q);
            } else {
                first = Some(q.clone());
            }
        }
    }

    #[test]
    fn block_cholesky_qr_preserves_span_and_thread_identity() {
        check_block_qr(257, 7);
    }

    /// cargo test --release -p dft-numerics block_qr_si8_shapes -- --ignored --nocapture
    #[test]
    #[ignore]
    fn block_qr_si8_shapes() {
        assert!(2969_usize.div_ceil(PROJECTION_PANEL) >= 4 * 8);
        for width in [20, 40, 60] {
            check_block_qr(2969, width);
        }
    }

    #[test]
    fn blocked_projection_rank_scale_orthogonality_and_identity() {
        let rows = 517;
        let against = orthonormalize(&matrix(rows, 9), None).unwrap().vectors;
        let mut block = matrix(rows, 6);
        // Exact dependency on the established basis, a zero and duplicate.
        block.column_mut(0).copy_from_slice(against.column(0));
        block.column_mut(1).fill(Complex64::default());
        for i in 0..rows {
            block[(i, 3)] = block[(i, 2)];
        }
        for z in block.column_mut(4) {
            *z *= 1e-200;
        }
        for z in block.column_mut(5) {
            *z *= 1e200;
        }
        let mut first = None;
        for threads in [1, 2, 4, 8] {
            let mut work = GemmWorkspace::new(rows, 16, threads).unwrap();
            let mut projection = ProjectionWorkspace::new(rows, 16).unwrap();
            let result = projection
                .orthonormalize(&block, Some(&against), &mut work)
                .unwrap();
            assert_eq!(result.dropped, [0, 1, 3]);
            let gram = gemm(
                &result.vectors,
                MatrixOp::Adjoint,
                &result.vectors,
                MatrixOp::Identity,
            )
            .unwrap();
            for j in 0..gram.cols {
                for i in 0..gram.rows {
                    assert!(
                        (gram[(i, j)] - Complex64::new(if i == j { 1.0 } else { 0.0 }, 0.0)).norm()
                            < 1e-12
                    );
                }
            }
            let cross = gemm(
                &against,
                MatrixOp::Adjoint,
                &result.vectors,
                MatrixOp::Identity,
            )
            .unwrap();
            assert!(cross.data.iter().all(|x| x.norm() < 1e-12));
            if let Some(prev) = &first {
                bits(prev, &result.vectors)
            } else {
                first = Some(result.vectors)
            }
        }
    }
    #[test]
    fn reproject_after_normalizing_nearly_dependent_block() {
        let rows = 517;
        let q = orthonormalize(&matrix(rows, 5), None).unwrap().vectors;
        let r = matrix(rows, 2);
        let mut block = ComplexMatrix::zeros(rows, 2).unwrap();
        for i in 0..rows {
            block[(i, 0)] = q[(i, 0)] + 1e-6 * r[(i, 0)];
            block[(i, 1)] = q[(i, 1)] + 1e-6 * (r[(i, 0)] + 1e-4 * r[(i, 1)]);
        }
        let mut first = None;
        for threads in [1, 2, 4, 8] {
            let mut work = GemmWorkspace::new(rows, 8, threads).unwrap();
            let mut projection = ProjectionWorkspace::new(rows, 8).unwrap();
            let result = projection
                .orthonormalize(&block, Some(&q), &mut work)
                .unwrap();
            assert!(result.dropped.is_empty());
            let cross = gemm(&q, MatrixOp::Adjoint, &result.vectors, MatrixOp::Identity).unwrap();
            assert!(cross.data.iter().all(|z| z.norm() < 1e-12));
            let gram = gemm(
                &result.vectors,
                MatrixOp::Adjoint,
                &result.vectors,
                MatrixOp::Identity,
            )
            .unwrap();
            for j in 0..2 {
                for i in 0..2 {
                    assert!(
                        (gram[(i, j)] - Complex64::new(if i == j { 1.0 } else { 0.0 }, 0.0)).norm()
                            < 1e-12
                    );
                }
            }
            if let Some(previous) = &first {
                bits(previous, &result.vectors);
            } else {
                first = Some(result.vectors);
            }
        }
    }

    /// Release-only scheduling experiment, not a performance assertion.
    /// cargo test --release -p dft-numerics granularity_table -- --ignored --nocapture
    #[test]
    #[ignore]
    #[cfg(not(feature = "mkl"))]
    fn granularity_table() {
        use std::time::Instant;
        let pool = WorkerPool::new(8).unwrap();
        println!("rows,width,complex_madds,inline_us,pool_us");
        for rows in [257, 1025, 4097] {
            for width in [8, 16, 32, 64] {
                let a = matrix(rows, width);
                let mut work = GemmWorkspace::with_pool(rows, width, 8, pool.clone()).unwrap();
                let mut out = ComplexMatrix::zeros(width, width).unwrap();
                let job = Job {
                    a: &a,
                    oa: MatrixOp::Adjoint,
                    b: &a,
                    ob: MatrixOp::Identity,
                    m: width,
                    k: rows,
                    n: width,
                    reduce: true,
                    panel: PANEL,
                };
                let count = rows.div_ceil(PANEL);
                let mut times = [0.0; 2];
                for (mode, time) in times.iter_mut().enumerate() {
                    let mut samples = Vec::new();
                    for _ in 0..5 {
                        let start = Instant::now();
                        for _ in 0..50 {
                            if mode == 0 {
                                run(&mut work.tiles[..count], 0, 1, &job).unwrap();
                            } else {
                                pool.install(|| {
                                    run(&mut work.tiles[..count], 0, 8.min(count), &job)
                                })
                                .unwrap();
                            }
                            out.data.fill(Complex64::default());
                            for tile in &work.tiles[..count] {
                                for (y, x) in out.data.iter_mut().zip(&tile.c.data) {
                                    *y += x;
                                }
                            }
                            std::hint::black_box(&out);
                        }
                        samples.push(start.elapsed().as_secs_f64() * 1e6 / 50.0);
                    }
                    samples.sort_by(f64::total_cmp);
                    *time = samples[2];
                }
                println!(
                    "{rows},{width},{},{:.3},{:.3}",
                    rows * width * width,
                    times[0],
                    times[1]
                );
            }
        }
    }

    #[test]
    fn projection_revalidates_mutated_inputs_and_rejects_overflow() {
        let mut work = GemmWorkspace::new(517, 8, 4).unwrap();
        let mut projection = ProjectionWorkspace::new(517, 8).unwrap();
        let block = matrix(517, 2);
        let mut q = matrix(517, 3);
        projection
            .orthonormalize(&block, Some(&q), &mut work)
            .unwrap();
        q[(516, 2)] = Complex64::new(0.0, f64::NAN);
        assert!(matches!(
            projection.orthonormalize(&block, Some(&q), &mut work),
            Err(LinalgError::NonFinite)
        ));
        q.data.fill(Complex64::new(f64::MAX, 0.0));
        let mut block = block;
        block.data.fill(Complex64::new(1.0, 0.0));
        assert!(matches!(
            projection.orthonormalize(&block, Some(&q), &mut work),
            Err(LinalgError::NonFinite)
        ));
    }

    #[test]
    fn dimensions_and_nonfinite_reject_before_arithmetic() {
        let mut work = GemmWorkspace::new(8, 4, 8).unwrap();
        let a = matrix(8, 3);
        assert!(matches!(
            work.multiply(&a, MatrixOp::Identity, &matrix(4, 2), MatrixOp::Identity),
            Err(LinalgError::Dimension)
        ));
        let mut a = a;
        a[(0, 0)] = Complex64::new(f64::NAN, 0.0);
        assert!(matches!(
            work.multiply(&a, MatrixOp::Identity, &matrix(3, 2), MatrixOp::Identity),
            Err(LinalgError::NonFinite)
        ));
    }
}
