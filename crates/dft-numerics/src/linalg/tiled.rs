//! Fixed backend decomposition, independent of the requested thread budget.
use super::*;
use crate::pool::WorkerPool;
use std::sync::Arc;

// A cache-sized plane-wave panel, independent of machine and worker count.
#[cfg(not(feature = "mkl"))]
const PANEL: usize = 256;

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
    #[cfg(not(feature = "mkl"))]
    threads: usize,
    #[cfg(not(feature = "mkl"))]
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
            let mut tiles = Vec::with_capacity(long.div_ceil(PANEL).max(1));
            for _ in 0..long.div_ceil(PANEL).max(1) {
                tiles.push(Tile {
                    c: ComplexMatrix::zeros(size, short)?,
                });
            }
            tiles
        };
        #[cfg(feature = "mkl")]
        let _ = (threads, pool);
        Ok(Self {
            #[cfg(not(feature = "mkl"))]
            tiles,
            long,
            short,
            #[cfg(not(feature = "mkl"))]
            threads: threads.max(1),
            #[cfg(not(feature = "mkl"))]
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
        #[cfg(feature = "mkl")]
        {
            gemm_into(a, oa, b, ob, out)
        }
        #[cfg(not(feature = "mkl"))]
        {
            a.validate_finite()?;
            b.validate_finite()?;
            let reduce = k > m;
            let count = if reduce { k } else { m }.div_ceil(PANEL).max(1);
            let job = Job {
                a,
                oa,
                b,
                ob,
                m,
                k,
                n,
                reduce,
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
}
#[cfg(not(feature = "mkl"))]
fn run(tiles: &mut [Tile], first: usize, threads: usize, job: &Job<'_>) -> Result<(), LinalgError> {
    if threads <= 1 || tiles.len() == 1 {
        for (offset, tile) in tiles.iter_mut().enumerate() {
            let start = (first + offset) * PANEL;
            let (rows, inner) = if job.reduce {
                (job.m, (job.k - start).min(PANEL))
            } else {
                ((job.m - start).min(PANEL), job.k)
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
    if threads <= 1 || out.nrows() <= PANEL {
        let mut start = first;
        while out.nrows() > 0 {
            let rows = out.nrows().min(PANEL);
            let (panel, rest) = out.split_at_row_mut(rows);
            borrowed_product(job, start, rows, job.k, panel);
            start += rows;
            out = rest;
        }
        return;
    }
    let panels = out.nrows().div_ceil(PANEL);
    let middle = panels / 2 * PANEL;
    let (left, right) = out.split_at_row_mut(middle);
    WorkerPool::join(
        || rotate(left, first, threads, job),
        || rotate(right, first + middle, threads, job),
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

/// Block CGS2 against the established basis, interleaved with two-pass MGS
/// within the small incoming block. Rank decisions remain relative to original norms.
/// All projection buffers are reserved once, outside Davidson iterations.
pub struct ProjectionWorkspace {
    normalized: ComplexMatrix,
    coefficients: ComplexMatrix,
    projected: ComplexMatrix,
    column: Vec<Complex64>,
    indices: Vec<usize>,
}
impl ProjectionWorkspace {
    pub fn new(rows: usize, capacity: usize) -> Result<Self, LinalgError> {
        Ok(Self {
            normalized: ComplexMatrix::zeros(rows, capacity)?,
            coefficients: ComplexMatrix::zeros(capacity, capacity)?,
            projected: ComplexMatrix::zeros(rows, capacity)?,
            column: vec![Complex64::default(); rows],
            indices: Vec::with_capacity(capacity),
        })
    }
    pub fn orthonormalize(
        &mut self,
        block: &ComplexMatrix,
        against: Option<&ComplexMatrix>,
        gemm: &mut GemmWorkspace,
    ) -> Result<Orthonormalized, LinalgError> {
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
                self.normalized.data.copy_from_slice(&data);
            }
            if let Some(q) = against {
                reshape(&mut self.coefficients, q.cols, count);
                reshape(&mut self.projected, block.rows, count);
                gemm.multiply_into(
                    q,
                    MatrixOp::Adjoint,
                    &self.normalized,
                    MatrixOp::Identity,
                    &mut self.coefficients,
                )?;
                gemm.multiply_into(
                    q,
                    MatrixOp::Identity,
                    &self.coefficients,
                    MatrixOp::Identity,
                    &mut self.projected,
                )?;
                for (y, p) in self.normalized.data.iter_mut().zip(&self.projected.data) {
                    *y -= p;
                }
            }
            data.clear();
            count = 0;
            for j in 0..self.normalized.cols {
                self.column.copy_from_slice(self.normalized.column(j));
                for _ in 0..2 {
                    for k in 0..count {
                        subtract_projection(
                            &mut self.column,
                            &data[k * block.rows..(k + 1) * block.rows],
                        );
                    }
                }
                let remaining = norm(&self.column);
                if !remaining.is_finite() {
                    return Err(LinalgError::NonFinite);
                }
                if remaining <= DEPENDENCE_TOLERANCE {
                    dropped.push(self.indices[j]);
                } else {
                    data.extend(self.column.iter().map(|x| x / remaining));
                    self.indices[count] = self.indices[j];
                    count += 1;
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
