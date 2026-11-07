//! Real column-major algebra for the isolated Gamma library path.
use crate::{linalg::LinalgError, pool::WorkerPool};
use faer::{Accum, MatMut, MatRef, Par, linalg::matmul::matmul};
use std::{
    ops::{Index, IndexMut},
    sync::Arc,
};

/// Borrowed real matrix; faer views preserve strides for submatrices.
pub type RealMatrixRef<'a> = MatRef<'a, f64>;
/// Exclusive borrowed real matrix, including strided submatrices.
pub type RealMatrixMut<'a> = MatMut<'a, f64>;

/// Owned column-major storage; entry (i,j) is at i+j*rows.
#[derive(Debug)]
pub struct RealMatrix {
    rows: usize,
    cols: usize,
    data: Vec<f64>,
}
fn length(rows: usize, cols: usize) -> Result<usize, LinalgError> {
    rows.checked_mul(cols)
        .filter(|&n| n <= isize::MAX as usize / size_of::<f64>())
        .ok_or(LinalgError::Dimension)
}
impl RealMatrix {
    /// Adopt storage after checking dimensions and addressable size.
    pub fn new(rows: usize, cols: usize, data: Vec<f64>) -> Result<Self, LinalgError> {
        if length(rows, cols)? != data.len() {
            return Err(LinalgError::Dimension);
        }
        Ok(Self { rows, cols, data })
    }
    /// Allocate a zero matrix with checked dimensions.
    pub fn zeros(rows: usize, cols: usize) -> Result<Self, LinalgError> {
        Self::new(rows, cols, vec![0.; length(rows, cols)?])
    }
    /// Change the active contiguous shape within initialized storage. No allocation,
    /// clearing or copying occurs; the caller must overwrite newly exposed entries.
    pub fn reshape_within_storage(&mut self, rows: usize, cols: usize) -> Result<(), LinalgError> {
        if length(rows, cols)? > self.data.len() {
            return Err(LinalgError::Dimension);
        }
        self.rows = rows;
        self.cols = cols;
        Ok(())
    }
    /// Reserve initialized backing storage before a solve starts, preserving entries.
    pub fn reserve_storage(&mut self, entries: usize) {
        if entries > self.data.len() {
            self.data.resize(entries, 0.);
        }
    }
    /// Row count.
    pub fn rows(&self) -> usize {
        self.rows
    }
    /// Column count.
    pub fn cols(&self) -> usize {
        self.cols
    }
    /// Contiguous column-major storage.
    pub fn as_slice(&self) -> &[f64] {
        &self.data[..self.rows * self.cols]
    }
    /// Mutable contiguous column-major storage.
    pub fn as_mut_slice(&mut self) -> &mut [f64] {
        &mut self.data[..self.rows * self.cols]
    }
    /// Borrow without copying or allocating.
    pub fn as_ref(&self) -> RealMatrixRef<'_> {
        MatRef::from_column_major_slice(self.as_slice(), self.rows, self.cols)
    }
    /// Borrow exclusively without copying or allocating.
    pub fn as_mut(&mut self) -> RealMatrixMut<'_> {
        MatMut::from_column_major_slice_mut(
            &mut self.data[..self.rows * self.cols],
            self.rows,
            self.cols,
        )
    }
}
impl Clone for RealMatrix {
    fn clone(&self) -> Self {
        Self {
            rows: self.rows,
            cols: self.cols,
            data: self.as_slice().to_vec(),
        }
    }
}
impl PartialEq for RealMatrix {
    fn eq(&self, other: &Self) -> bool {
        self.rows == other.rows && self.cols == other.cols && self.as_slice() == other.as_slice()
    }
}
impl Index<(usize, usize)> for RealMatrix {
    type Output = f64;
    fn index(&self, (i, j): (usize, usize)) -> &f64 {
        assert!(i < self.rows && j < self.cols);
        &self.data[i + j * self.rows]
    }
}
impl IndexMut<(usize, usize)> for RealMatrix {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut f64 {
        assert!(i < self.rows && j < self.cols);
        &mut self.data[i + j * self.rows]
    }
}
fn finite(a: RealMatrixRef<'_>) -> Result<(), LinalgError> {
    if let Some(a) = a.try_as_col_major() {
        return if a
            .col_iter()
            .all(|col| col.as_slice().iter().all(|x| x.is_finite()))
        {
            Ok(())
        } else {
            Err(LinalgError::NonFinite)
        };
    }
    for j in 0..a.ncols() {
        for i in 0..a.nrows() {
            if !a[(i, j)].is_finite() {
                return Err(LinalgError::NonFinite);
            }
        }
    }
    Ok(())
}
/// Left operand operation; the right operand is untransposed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RealOp {
    /// A B, partitioned along output rows.
    Identity,
    /// A^T B, partitioned along the reduction dimension.
    Transpose,
}
const PANEL: usize = 256;
// Match the complex scheduling guard's 2,097,152 real flops: a real
// multiply-add costs two flops, hence 1,048,576 real multiply-adds.
// This is a scheduling heuristic, not a measured real-path speed claim.
const PARALLEL_WORK: usize = 1_048_576;

/// Reusable fixed-panel Gram scratch and shared executor. Construction allocates;
/// multiply_into neither allocates numerical storage nor creates worker threads.
/// The faer backend stays sequential within each panel, including MKL builds.
pub struct RealGemmWorkspace {
    long: usize,
    short: usize,
    tiles: Vec<RealMatrix>,
    pool: Arc<WorkerPool>,
}
impl RealGemmWorkspace {
    /// Reserve up to `long` full rows and `short` subspace columns.
    pub fn with_pool(
        long: usize,
        short: usize,
        pool: Arc<WorkerPool>,
    ) -> Result<Self, LinalgError> {
        length(short, short)?;
        length(
            long.div_ceil(PANEL),
            short.checked_mul(short).ok_or(LinalgError::Dimension)?,
        )?;
        let tiles = (0..long.div_ceil(PANEL))
            .map(|_| RealMatrix::zeros(short, short))
            .collect::<Result<_, _>>()?;
        Ok(Self {
            long,
            short,
            tiles,
            pool,
        })
    }
    /// Rotate V and HV by the same Ritz coefficients in one parallel region.
    /// Each side retains the original 256-row pieces and bounded kernel calls
    /// (24 independent pieces for a 2969-row pair). The guard is unchanged.
    pub fn rotate_pair(
        &mut self,
        v: RealMatrixRef<'_>,
        hv: RealMatrixRef<'_>,
        c: RealMatrixRef<'_>,
        mut out_v: RealMatrixMut<'_>,
        mut out_hv: RealMatrixMut<'_>,
    ) -> Result<(), LinalgError> {
        let (m, k, n) = (v.nrows(), v.ncols(), c.ncols());
        if hv.nrows() != m
            || hv.ncols() != k
            || c.nrows() != k
            || out_v.nrows() != m
            || out_hv.nrows() != m
            || out_v.ncols() != n
            || out_hv.ncols() != n
            || m > self.long
            || k > self.short
            || n > self.short
        {
            return Err(LinalgError::Dimension);
        }
        finite(v)?;
        finite(hv)?;
        finite(c)?;
        if m == 0 || n == 0 {
            return Ok(());
        }
        if k == 0 {
            out_v.fill(0.);
            out_hv.fill(0.);
            return Ok(());
        }
        let parallel =
            m.div_ceil(PANEL) >= 4 && m.saturating_mul(n).saturating_mul(k) >= PARALLEL_WORK;
        if parallel {
            self.pool.install(|| {
                WorkerPool::join(
                    || rotate(out_v.as_mut(), v, c, true),
                    || rotate(out_hv.as_mut(), hv, c, true),
                );
            });
        } else {
            rotate(out_v.as_mut(), v, c, false);
            rotate(out_hv.as_mut(), hv, c, false);
        }
        finite(out_v.as_ref())?;
        finite(out_hv.as_ref())
    }
    /// Overwrite output with A B or A^T B. Dimensions, capacity, and finiteness
    /// are checked. Zero reduction length produces zeros. Inputs may be views.
    pub fn multiply_into(
        &mut self,
        a: RealMatrixRef<'_>,
        op: RealOp,
        b: RealMatrixRef<'_>,
        mut out: RealMatrixMut<'_>,
    ) -> Result<(), LinalgError> {
        let (m, k) = if op == RealOp::Transpose {
            (a.ncols(), a.nrows())
        } else {
            (a.nrows(), a.ncols())
        };
        let n = b.ncols();
        if k != b.nrows()
            || out.nrows() != m
            || out.ncols() != n
            || a.nrows() > self.long
            || a.ncols() > self.short
            || n > self.short
        {
            return Err(LinalgError::Dimension);
        }
        finite(a)?;
        finite(b)?;
        if m == 0 || n == 0 {
            return Ok(());
        }
        if k == 0 {
            out.fill(0.);
            return Ok(());
        }
        let count = if op == RealOp::Transpose { k } else { m }.div_ceil(PANEL);
        let parallel = count >= 4 && m.saturating_mul(n).saturating_mul(k) >= PARALLEL_WORK;
        if op == RealOp::Transpose {
            let tiles = &mut self.tiles[..count];
            let mut work = || gram(tiles, 0, a, b, m, n, parallel);
            if parallel {
                self.pool.install(work);
            } else {
                work();
            }
            let mut sum = || assemble(out.as_mut(), tiles, 0, parallel);
            if parallel {
                self.pool.install(sum);
            } else {
                sum();
            }
        } else {
            let mut work = || rotate(out.as_mut(), a, b, parallel);
            if parallel {
                self.pool.install(work);
            } else {
                work();
            }
        }
        finite(out.as_ref())
    }
}
// General faer GEMM allocates opaque packing buffers. These call shapes stay
// within its allocation-free 4096-MAC nano-GEMM branch. Rotation uses the full
// reduction panel instead of repeated 16-term updates; Gram retains the shape
// that won its separate panel benchmark. Neither grouping depends on workers.
fn bounded_matmul(
    mut out: RealMatrixMut<'_>,
    a: RealMatrixRef<'_>,
    b: RealMatrixRef<'_>,
    transpose: bool,
) {
    let (rows, cols, reduction) = if transpose { (16, 16, 16) } else { (4, 4, 256) };
    for j in (0..out.ncols()).step_by(cols) {
        let n = cols.min(out.ncols() - j);
        for i in (0..out.nrows()).step_by(rows) {
            let m = rows.min(out.nrows() - i);
            for k in (0..a.ncols()).step_by(reduction) {
                let len = reduction.min(a.ncols() - k);
                matmul(
                    out.as_mut().submatrix_mut(i, j, m, n),
                    if k == 0 { Accum::Replace } else { Accum::Add },
                    a.submatrix(i, k, m, len),
                    b.submatrix(k, j, len, n),
                    1.,
                    Par::Seq,
                );
            }
        }
    }
}

fn gram(
    tiles: &mut [RealMatrix],
    start: usize,
    a: RealMatrixRef<'_>,
    b: RealMatrixRef<'_>,
    m: usize,
    n: usize,
    parallel: bool,
) {
    if parallel && tiles.len() > 1 {
        let half = tiles.len() / 2;
        let (left, right) = tiles.split_at_mut(half);
        WorkerPool::join(
            || gram(left, start, a, b, m, n, true),
            || gram(right, start + half, a, b, m, n, true),
        );
    } else {
        for (t, tile) in tiles.iter_mut().enumerate() {
            let row = (start + t) * PANEL;
            let len = PANEL.min(a.nrows() - row);
            bounded_matmul(
                tile.as_mut().submatrix_mut(0, 0, m, n),
                a.subrows(row, len).transpose(),
                b.subrows(row, len),
                true,
            );
        }
    }
}
// Output columns are independent. Four-column pieces and ascending panel sums
// are fixed for every thread budget. Retain the initial +0 addition (including
// its signed-zero behavior) while overwriting, rather than clearing, output.
fn assemble(mut out: RealMatrixMut<'_>, tiles: &[RealMatrix], first: usize, parallel: bool) {
    let count = out.ncols();
    if parallel && count > 4 {
        let split = (count.div_ceil(4) / 2) * 4;
        let (left, right) = out.split_at_col_mut(split);
        WorkerPool::join(
            || assemble(left, tiles, first, true),
            || assemble(right, tiles, first + split, true),
        );
    } else {
        let m = out.nrows();
        for j in 0..count {
            for (t, tile) in tiles.iter().enumerate() {
                let src =
                    &tile.as_slice()[(first + j) * tile.rows()..(first + j) * tile.rows() + m];
                if let Some(col) = out.as_mut().col_mut(j).try_as_col_major_mut() {
                    for (dst, &value) in col.as_slice_mut().iter_mut().zip(src) {
                        *dst = if t == 0 { 0. + value } else { *dst + value };
                    }
                } else {
                    for (i, &value) in src.iter().enumerate() {
                        out[(i, j)] = if t == 0 {
                            0. + value
                        } else {
                            out[(i, j)] + value
                        };
                    }
                }
            }
        }
    }
}
fn rotate(mut out: RealMatrixMut<'_>, a: RealMatrixRef<'_>, b: RealMatrixRef<'_>, parallel: bool) {
    let count = a.nrows().div_ceil(PANEL);
    if parallel && count > 1 {
        let split = (count / 2) * PANEL;
        let (left, right) = out.split_at_row_mut(split);
        WorkerPool::join(
            || rotate(left, a.subrows(0, split), b, true),
            || rotate(right, a.subrows(split, a.nrows() - split), b, true),
        );
    } else {
        for row in (0..a.nrows()).step_by(PANEL) {
            let len = PANEL.min(a.nrows() - row);
            bounded_matmul(
                out.as_mut().subrows_mut(row, len),
                a.subrows(row, len),
                b,
                false,
            );
        }
    }
}

/// Real symmetric eigensolve, ascending eigenvalues, sequential faer backend.
/// Largest-magnitude component of each eigenvector is positive (first index
/// wins ties). Degenerate eigenspaces are not assigned a unique basis.
/// Symmetry uses the existing adapter's 1e-12 relative entrywise allowance,
/// with unit absolute floor; faer consumes the lower triangle after validation.
pub fn symmetric_eigh(a: &RealMatrix) -> Result<(Vec<f64>, RealMatrix), LinalgError> {
    let mut ws = RealEighWorkspace::new(a.rows())?;
    let mut values = vec![0.; a.rows()];
    let mut vectors = RealMatrix::zeros(a.rows(), a.rows())?;
    ws.solve_into(a, &mut values, &mut vectors)?;
    Ok((values, vectors))
}
/// Caller-owned EVD storage. Every active size is covered at construction;
/// backend algorithm, strides and sequential arithmetic match symmetric_eigh.
pub struct RealEighWorkspace {
    capacity: usize,
    diagonal: Vec<f64>,
    buffer: faer::dyn_stack::MemBuffer,
}
impl RealEighWorkspace {
    /// Reserve scratch for every possible active subspace size.
    pub fn new(capacity: usize) -> Result<Self, LinalgError> {
        length(capacity, capacity)?;
        let requirements: Vec<_> = (1..=capacity)
            .map(|n| {
                faer::linalg::evd::self_adjoint_evd_scratch::<f64>(
                    n,
                    faer::linalg::evd::ComputeEigenvectors::Yes,
                    Par::Seq,
                    Default::default(),
                )
            })
            .collect();
        Ok(Self {
            capacity,
            diagonal: vec![0.; capacity],
            buffer: faer::dyn_stack::MemBuffer::new(faer::dyn_stack::StackReq::any_of(
                &requirements,
            )),
        })
    }
    /// Reuse output and scratch storage. The faer backend can still allocate
    /// internal GEMM packing/bookkeeping; this is not a whole-call no-allocation
    /// guarantee (see the P33A counting-allocator diagnostic).
    pub fn solve_into(
        &mut self,
        a: &RealMatrix,
        values: &mut [f64],
        vectors: &mut RealMatrix,
    ) -> Result<(), LinalgError> {
        if a.rows != a.cols {
            return Err(LinalgError::NotSquare);
        }
        finite(a.as_ref())?;
        let allowed = crate::linalg::HERMITIAN_TOLERANCE
            * a.as_slice().iter().map(|v| v.abs()).fold(1., f64::max);
        let mut defect: f64 = 0.;
        for j in 0..a.cols {
            for i in 0..j {
                defect = defect.max((a[(i, j)] - a[(j, i)]).abs());
            }
        }
        if defect > allowed {
            return Err(LinalgError::NotHermitian { defect, allowed });
        }
        let n = a.rows;
        if n > self.capacity || values.len() != n || vectors.rows() != n || vectors.cols() != n {
            return Err(LinalgError::Dimension);
        }
        if n == 0 {
            return Ok(());
        }
        faer::linalg::evd::self_adjoint_evd(
            a.as_ref(),
            faer::diag::DiagMut::from_slice_mut(&mut self.diagonal[..n]),
            Some(vectors.as_mut()),
            Par::Seq,
            faer::dyn_stack::MemStack::new(&mut self.buffer),
            Default::default(),
        )
        .map_err(|_| LinalgError::Decomposition)?;
        values.copy_from_slice(&self.diagonal[..n]);
        if values.iter().any(|v| !v.is_finite()) {
            return Err(LinalgError::NonFinite);
        }
        finite(vectors.as_ref())?;
        for j in 0..n {
            let mut pivot = 0;
            for i in 1..n {
                if vectors[(i, j)].abs() > vectors[(pivot, j)].abs() {
                    pivot = i;
                }
            }
            if vectors[(pivot, j)].is_sign_negative() {
                for i in 0..n {
                    vectors[(i, j)] = -vectors[(i, j)];
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::{ComplexMatrix, GemmWorkspace, MatrixOp};
    use num_complex::Complex64;
    fn fixture(r: usize, c: usize) -> RealMatrix {
        RealMatrix::new(
            r,
            c,
            (0..r * c)
                .map(|i| ((i * 13 + 7) as f64).sin() / 32.)
                .collect(),
        )
        .unwrap()
    }
    #[test]
    fn products_reference_and_pool_identity() {
        // Includes empty, partial, exact panels and parallel-dispatch shapes.
        for (long, short) in [(0, 3), (1, 1), (255, 5), (256, 4), (257, 7), (1025, 33)] {
            for op in [RealOp::Identity, RealOp::Transpose] {
                let a = fixture(long, short);
                let b = fixture(if op == RealOp::Identity { short } else { long }, short);
                let m = if op == RealOp::Identity { long } else { short };
                let k = b.rows;
                let mut expected = RealMatrix::zeros(m, short).unwrap();
                for j in 0..short {
                    for i in 0..m {
                        for p in 0..k {
                            expected[(i, j)] += if op == RealOp::Identity {
                                a[(i, p)]
                            } else {
                                a[(p, i)]
                            } * b[(p, j)];
                        }
                    }
                }
                let ca = ComplexMatrix::new(
                    a.rows,
                    a.cols,
                    a.data.iter().map(|&x| Complex64::new(x, 0.)).collect(),
                )
                .unwrap();
                let cb = ComplexMatrix::new(
                    b.rows,
                    b.cols,
                    b.data.iter().map(|&x| Complex64::new(x, 0.)).collect(),
                )
                .unwrap();
                let mut cw = GemmWorkspace::new(long.max(short), short, 1).unwrap();
                let complex = cw
                    .multiply(
                        &ca,
                        if op == RealOp::Identity {
                            MatrixOp::Identity
                        } else {
                            MatrixOp::Adjoint
                        },
                        &cb,
                        MatrixOp::Identity,
                    )
                    .unwrap();
                let mut first = None;
                for threads in [1, 2, 4, 8] {
                    let pool = WorkerPool::new(threads).unwrap();
                    let mut ws = RealGemmWorkspace::with_pool(long, short, pool.clone()).unwrap();
                    let mut out = RealMatrix::zeros(m, short).unwrap();
                    ws.multiply_into(a.as_ref(), op, b.as_ref(), out.as_mut())
                        .unwrap();
                    // Unit-scale products: 1e-14 allows backend FMA/reduction rounding.
                    let defect = out
                        .data
                        .iter()
                        .zip(&expected.data)
                        .map(|(a, b)| (a - b).abs())
                        .fold(0., f64::max);
                    assert!(defect < 1e-14, "{defect:e}");
                    assert!(
                        out.data
                            .iter()
                            .zip(complex.as_slice())
                            .all(|(a, b)| (a - b.re).abs() < 1e-14 && b.im == 0.)
                    );
                    let bits: Vec<_> = out.data.iter().map(|x| x.to_bits()).collect();
                    if let Some(ref first) = first {
                        assert_eq!(&bits, first);
                    } else {
                        first = Some(bits);
                    }
                    // Resident/nested dispatch must retain the same arithmetic too.
                    pool.install(|| ws.multiply_into(a.as_ref(), op, b.as_ref(), out.as_mut()))
                        .unwrap();
                    assert_eq!(
                        out.data.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                        *first.as_ref().unwrap()
                    );
                    if threads == 1 {
                        println!("GEMM {op:?} {long}x{short}: naive defect={defect:e}");
                    }
                }
            }
        }
    }
    #[test]
    fn eigen_residual_sign_and_errors() {
        for n in [0, 1, 2, 9, 24] {
            let mut a = fixture(n, n);
            for j in 0..n {
                for i in 0..j {
                    a[(i, j)] = a[(j, i)];
                }
                a[(j, j)] += (j / 3) as f64;
            }
            let (values, v) = symmetric_eigh(&a).unwrap();
            assert!(values.windows(2).all(|w| w[0] <= w[1]));
            let mut residual: f64 = 0.;
            let mut orth: f64 = 0.;
            for j in 0..n {
                let mut pivot = 0;
                for i in 1..n {
                    if v[(i, j)].abs() > v[(pivot, j)].abs() {
                        pivot = i;
                    }
                }
                assert!(v[(pivot, j)] >= 0.);
                for i in 0..n {
                    residual = residual.max(
                        ((0..n).map(|k| a[(i, k)] * v[(k, j)]).sum::<f64>()
                            - values[j] * v[(i, j)])
                            .abs(),
                    );
                    orth = orth.max(
                        ((0..n).map(|k| v[(k, i)] * v[(k, j)]).sum::<f64>()
                            - if i == j { 1. } else { 0. })
                        .abs(),
                    );
                }
            }
            println!("EVD n={n}: residual={residual:e}, orthogonality={orth:e}");
            // Small matrices with O(1) spectrum: tight backward error allowance.
            assert!(residual < 1e-13 && orth < 1e-13);
            assert_eq!(symmetric_eigh(&a).unwrap(), (values, v));
        }
        assert!(RealMatrix::zeros(usize::MAX, 2).is_err());
        assert!(RealMatrix::new(2, 2, vec![0.]).is_err());
        assert!(symmetric_eigh(&fixture(2, 3)).is_err());
        let mut a = fixture(2, 2);
        a[(0, 1)] = 9.;
        assert!(symmetric_eigh(&a).is_err());
        a[(0, 0)] = f64::NAN;
        assert!(symmetric_eigh(&a).is_err());
        let pool = WorkerPool::new(1).unwrap();
        let mut ws = RealGemmWorkspace::with_pool(2, 2, pool).unwrap();
        let b = fixture(2, 2);
        let mut out = fixture(2, 2);
        assert!(
            ws.multiply_into(a.as_ref(), RealOp::Identity, b.as_ref(), out.as_mut())
                .is_err()
        );
        assert!(
            ws.multiply_into(
                b.as_ref().subrows(0, 1),
                RealOp::Transpose,
                b.as_ref(),
                out.as_mut()
            )
            .is_err()
        );
    }
}
