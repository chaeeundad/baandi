//! Dense backend boundary; all matrices are column-major and backend independent.
#[cfg(not(feature = "mkl"))]
use faer::{Accum, linalg::matmul::matmul};
use faer::{Mat, Par};
use num_complex::Complex64;
use std::ops::{Index, IndexMut};
use thiserror::Error;

/// Relative entrywise Hermiticity tolerance, with an absolute floor of one.
pub const HERMITIAN_TOLERANCE: f64 = 1e-12;
/// A column is dropped if its remaining norm is <= this times its original norm.
pub const DEPENDENCE_TOLERANCE: f64 = 1e-12;

#[derive(Debug, Clone, PartialEq, Error)]
pub enum LinalgError {
    #[error("incompatible matrix dimensions or storage length")]
    Dimension,
    #[error("matrix must be square")]
    NotSquare,
    #[error("matrix contains a non-finite value or arithmetic overflowed")]
    NonFinite,
    #[error("matrix is not Hermitian: defect {defect:e}, allowed {allowed:e}")]
    NotHermitian { defect: f64, allowed: f64 },
    #[error("backend eigendecomposition failed")]
    Decomposition,
}

/// Entry (i,j) is at i + j*rows. Dimensions and storage length cannot diverge.
#[derive(Debug, Clone, PartialEq)]
pub struct ComplexMatrix {
    rows: usize,
    cols: usize,
    data: Vec<Complex64>,
}
impl ComplexMatrix {
    pub fn new(rows: usize, cols: usize, data: Vec<Complex64>) -> Result<Self, LinalgError> {
        if rows.checked_mul(cols) != Some(data.len())
            || data.len() > isize::MAX as usize / size_of::<Complex64>()
        {
            return Err(LinalgError::Dimension);
        }
        Ok(Self { rows, cols, data })
    }
    pub fn zeros(rows: usize, cols: usize) -> Result<Self, LinalgError> {
        let len = rows
            .checked_mul(cols)
            .filter(|&n| n <= isize::MAX as usize / size_of::<Complex64>())
            .ok_or(LinalgError::Dimension)?;
        Self::new(rows, cols, vec![Complex64::new(0.0, 0.0); len])
    }
    /// Reserve a fixed column capacity during workspace construction. Live
    /// columns remain contiguous with leading dimension `rows`; spare columns
    /// are not visible to validation or backend matrix views.
    pub fn reserve_columns(&mut self, capacity: usize) -> Result<(), LinalgError> {
        let len = self
            .rows
            .checked_mul(capacity)
            .filter(|&n| n <= isize::MAX as usize / size_of::<Complex64>())
            .ok_or(LinalgError::Dimension)?;
        if capacity < self.cols {
            return Err(LinalgError::Dimension);
        }
        self.data.reserve_exact(len.saturating_sub(self.data.len()));
        Ok(())
    }

    /// Change the live prefix without allocation. New columns are zeroed;
    /// shrinking never copies retained columns or releases their capacity.
    pub fn resize_columns(&mut self, cols: usize) -> Result<(), LinalgError> {
        let len = self.rows.checked_mul(cols).ok_or(LinalgError::Dimension)?;
        if len > self.data.capacity() {
            return Err(LinalgError::Dimension);
        }
        self.data.resize(len, Complex64::default());
        self.cols = cols;
        Ok(())
    }

    /// Append only the new column range into reserved subspace storage.
    pub fn append_columns(&mut self, extra: &Self) -> Result<(), LinalgError> {
        if self.rows != extra.rows {
            return Err(LinalgError::Dimension);
        }
        let cols = self
            .cols
            .checked_add(extra.cols)
            .ok_or(LinalgError::Dimension)?;
        let len = self.rows.checked_mul(cols).ok_or(LinalgError::Dimension)?;
        if len > self.data.capacity() {
            return Err(LinalgError::Dimension);
        }
        crate::pool::analysis::copied(extra.data.len() * size_of::<Complex64>());
        self.data.extend_from_slice(&extra.data);
        self.cols = cols;
        Ok(())
    }

    /// Replace a live prefix while retaining the workspace allocation.
    pub fn replace_columns(&mut self, source: &Self) -> Result<(), LinalgError> {
        if self.rows != source.rows || source.data.len() > self.data.capacity() {
            return Err(LinalgError::Dimension);
        }
        crate::pool::analysis::copied(source.data.len() * size_of::<Complex64>());
        self.data.clear();
        self.data.extend_from_slice(&source.data);
        self.cols = source.cols;
        Ok(())
    }

    pub fn rows(&self) -> usize {
        self.rows
    }
    pub fn cols(&self) -> usize {
        self.cols
    }
    pub fn as_slice(&self) -> &[Complex64] {
        &self.data
    }
    /// Mutable column-major storage; dimensions and allocation remain unchanged.
    pub fn as_mut_slice(&mut self) -> &mut [Complex64] {
        &mut self.data
    }
    pub fn column(&self, j: usize) -> &[Complex64] {
        assert!(j < self.cols);
        &self.data[j * self.rows..(j + 1) * self.rows]
    }
    pub fn column_mut(&mut self, j: usize) -> &mut [Complex64] {
        assert!(j < self.cols);
        &mut self.data[j * self.rows..(j + 1) * self.rows]
    }
    pub fn validate_finite(&self) -> Result<(), LinalgError> {
        // Reduce without an early-exit branch per element, permitting SIMD.
        // Preserve rejection of either nonfinite component, including NaNs.
        let finite = self.data.iter().fold(true, |finite, z| {
            finite & z.re.is_finite() & z.im.is_finite()
        });
        if finite {
            Ok(())
        } else {
            Err(LinalgError::NonFinite)
        }
    }
}
impl Index<(usize, usize)> for ComplexMatrix {
    type Output = Complex64;
    fn index(&self, (i, j): (usize, usize)) -> &Self::Output {
        assert!(i < self.rows && j < self.cols);
        &self.data[i + j * self.rows]
    }
}
impl IndexMut<(usize, usize)> for ComplexMatrix {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut Self::Output {
        assert!(i < self.rows && j < self.cols);
        &mut self.data[i + j * self.rows]
    }
}
fn backend(a: &ComplexMatrix) -> Mat<Complex64> {
    crate::pool::analysis::copied(a.data.len() * size_of::<Complex64>());
    Mat::from_fn(a.rows, a.cols, |i, j| a[(i, j)])
}
fn owned(a: faer::MatRef<'_, Complex64>) -> Result<ComplexMatrix, LinalgError> {
    crate::pool::analysis::copied(a.nrows() * a.ncols() * size_of::<Complex64>());
    let mut result = ComplexMatrix::zeros(a.nrows(), a.ncols())?;
    for j in 0..a.ncols() {
        for i in 0..a.nrows() {
            result[(i, j)] = a[(i, j)];
        }
    }
    result.validate_finite()?;
    Ok(result)
}

/// Eigenvalues ascending, eigenvectors in corresponding columns. Empty square matrices are allowed.
pub fn hermitian_eigh(a: &ComplexMatrix) -> Result<(Vec<f64>, ComplexMatrix), LinalgError> {
    if a.rows != a.cols {
        return Err(LinalgError::NotSquare);
    }
    a.validate_finite()?;
    let scale = a.data.iter().map(|z| z.norm()).fold(1.0, f64::max);
    let allowed = HERMITIAN_TOLERANCE * scale;
    let mut defect: f64 = 0.0;
    for j in 0..a.cols {
        for i in 0..=j {
            defect = defect.max((a[(i, j)] - a[(j, i)].conj()).norm());
        }
    }
    if !scale.is_finite() || !defect.is_finite() {
        return Err(LinalgError::NonFinite);
    }
    if defect > allowed {
        return Err(LinalgError::NotHermitian { defect, allowed });
    }
    if a.rows == 0 {
        return Ok((vec![], a.clone()));
    }
    let a = backend(a);
    // Explicit sequential execution avoids process-global faer configuration.
    let mut vectors = Mat::zeros(a.nrows(), a.nrows());
    let mut diagonal = faer::diag::Diag::zeros(a.nrows());
    let scratch = faer::linalg::evd::self_adjoint_evd_scratch::<Complex64>(
        a.nrows(),
        faer::linalg::evd::ComputeEigenvectors::Yes,
        Par::Seq,
        Default::default(),
    );
    let mut buffer = faer::dyn_stack::MemBuffer::new(scratch);
    faer::linalg::evd::self_adjoint_evd(
        a.as_ref(),
        diagonal.as_mut(),
        Some(vectors.as_mut()),
        Par::Seq,
        faer::dyn_stack::MemStack::new(&mut buffer),
        Default::default(),
    )
    .map_err(|_| LinalgError::Decomposition)?;
    let values: Vec<_> = (0..a.nrows()).map(|i| diagonal[i].re).collect();
    if !values.iter().all(|x| x.is_finite()) {
        return Err(LinalgError::NonFinite);
    }
    Ok((values, owned(vectors.as_ref())?))
}

#[derive(Debug, Clone, Copy)]
pub enum MatrixOp {
    Identity,
    Adjoint,
}
/// Returns op(a)*op(b); backend GEMM is explicitly sequential.
#[cfg(not(feature = "mkl"))]
pub fn gemm(
    a: &ComplexMatrix,
    op_a: MatrixOp,
    b: &ComplexMatrix,
    op_b: MatrixOp,
) -> Result<ComplexMatrix, LinalgError> {
    a.validate_finite()?;
    b.validate_finite()?;
    let a = backend(a);
    let b = backend(b);
    let (m, k) = match op_a {
        MatrixOp::Identity => (a.nrows(), a.ncols()),
        MatrixOp::Adjoint => (a.ncols(), a.nrows()),
    };
    let (l, n) = match op_b {
        MatrixOp::Identity => (b.nrows(), b.ncols()),
        MatrixOp::Adjoint => (b.ncols(), b.nrows()),
    };
    if k != l {
        return Err(LinalgError::Dimension);
    }
    let mut c = Mat::zeros(m, n);
    let one = Complex64::new(1.0, 0.0);
    match (op_a, op_b) {
        (MatrixOp::Identity, MatrixOp::Identity) => {
            matmul(&mut c, Accum::Replace, &a, &b, one, Par::Seq)
        }
        (MatrixOp::Adjoint, MatrixOp::Identity) => {
            matmul(&mut c, Accum::Replace, a.adjoint(), &b, one, Par::Seq)
        }
        (MatrixOp::Identity, MatrixOp::Adjoint) => {
            matmul(&mut c, Accum::Replace, &a, b.adjoint(), one, Par::Seq)
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
    owned(c.as_ref())
}

/// MKL LP64 GEMM; eigendecomposition remains on faer.
#[cfg(feature = "mkl")]
pub fn gemm(
    a: &ComplexMatrix,
    op_a: MatrixOp,
    b: &ComplexMatrix,
    op_b: MatrixOp,
) -> Result<ComplexMatrix, LinalgError> {
    a.validate_finite()?;
    b.validate_finite()?;
    let adj_a = matches!(op_a, MatrixOp::Adjoint);
    let adj_b = matches!(op_b, MatrixOp::Adjoint);
    let (m, k) = if adj_a {
        (a.cols, a.rows)
    } else {
        (a.rows, a.cols)
    };
    let (l, n) = if adj_b {
        (b.cols, b.rows)
    } else {
        (b.rows, b.cols)
    };
    if k != l {
        return Err(LinalgError::Dimension);
    }
    let ma =
        dft_mkl::Matrix::new(a.as_slice(), a.rows, a.cols).map_err(|_| LinalgError::Dimension)?;
    let mb =
        dft_mkl::Matrix::new(b.as_slice(), b.rows, b.cols).map_err(|_| LinalgError::Dimension)?;
    let mut result = ComplexMatrix::zeros(m, n)?;
    dft_mkl::gemm(ma, adj_a, mb, adj_b, &mut result.data).map_err(|_| LinalgError::Dimension)?;
    result.validate_finite()?;
    Ok(result)
}

#[derive(Debug)]
pub struct Orthonormalized {
    pub vectors: ComplexMatrix,
    /// Original input column indices removed by rank filtering.
    pub dropped: Vec<usize>,
}
/// Stable Euclidean norm (avoids squaring very large/small components).
pub fn norm(x: &[Complex64]) -> f64 {
    let mut sums = [0.0; 4];
    let (chunks, remainder) = x.as_chunks::<2>();
    for pair in chunks {
        sums[0] += pair[0].re * pair[0].re;
        sums[1] += pair[0].im * pair[0].im;
        sums[2] += pair[1].re * pair[1].re;
        sums[3] += pair[1].im * pair[1].im;
    }
    for z in remainder {
        sums[0] += z.re * z.re;
        sums[1] += z.im * z.im;
    }
    let sum = (sums[0] + sums[1]) + (sums[2] + sums[3]);
    // At most 2*n subnormal products lose <= MIN_POSITIVE*EPSILON each.
    // This conservative threshold bounds that loss relative to sum by EPSILON.
    if sum.is_finite() && sum > f64::MIN_POSITIVE * (2.0 * x.len() as f64).max(1.0) {
        sum.sqrt()
    } else {
        x.iter().fold(0.0_f64, |s, z| s.hypot(z.re).hypot(z.im))
    }
}
/// Two-pass modified Gram-Schmidt. `against` must already be orthonormal.
/// Filtering is relative to each original column norm, hence scale independent.
pub fn orthonormalize(
    block: &ComplexMatrix,
    against: Option<&ComplexMatrix>,
) -> Result<Orthonormalized, LinalgError> {
    block.validate_finite()?;
    if let Some(q) = against {
        q.validate_finite()?;
        if q.rows != block.rows {
            return Err(LinalgError::Dimension);
        }
    }
    let mut data = Vec::new();
    let mut dropped = Vec::new();
    let mut count = 0;
    for j in 0..block.cols {
        let mut v = block.column(j).to_vec();
        let original = norm(&v);
        if !original.is_finite() {
            return Err(LinalgError::NonFinite);
        }
        // Normalize first so rank decisions also work for very small input columns.
        if original == 0.0 {
            dropped.push(j);
            continue;
        }
        for z in &mut v {
            *z /= original;
        }
        for _ in 0..2 {
            if let Some(q) = against {
                for k in 0..q.cols {
                    subtract_projection(&mut v, q.column(k));
                }
            }
            for k in 0..count {
                subtract_projection(&mut v, &data[k * block.rows..(k + 1) * block.rows]);
            }
        }
        let remaining = norm(&v);
        if !remaining.is_finite() {
            return Err(LinalgError::NonFinite);
        }
        if remaining <= DEPENDENCE_TOLERANCE {
            dropped.push(j);
        } else {
            data.extend(v.iter().map(|z| z / remaining));
            count += 1;
        }
    }
    Ok(Orthonormalized {
        vectors: ComplexMatrix::new(block.rows, count, data)?,
        dropped,
    })
}
fn subtract_projection(v: &mut [Complex64], q: &[Complex64]) {
    let mut sums = [Complex64::new(0.0, 0.0); 4];
    let (qc, qr) = q.as_chunks::<4>();
    let (vc, vr) = v.as_chunks::<4>();
    for (q, v) in qc.iter().zip(vc) {
        for k in 0..4 {
            sums[k] += q[k].conj() * v[k];
        }
    }
    let mut dot = (sums[0] + sums[1]) + (sums[2] + sums[3]);
    for (q, v) in qr.iter().zip(vr) {
        dot += q.conj() * v;
    }
    for (v, q) in v.iter_mut().zip(q) {
        *v -= dot * q;
    }
}

#[cfg(test)]
mod tests;

/// MKL variant of `gemm_into`: writes `op(a) * op(b)` into `output` (beta = 0).
#[cfg(feature = "mkl")]
pub fn gemm_into(
    a: &ComplexMatrix,
    op_a: MatrixOp,
    b: &ComplexMatrix,
    op_b: MatrixOp,
    output: &mut ComplexMatrix,
) -> Result<(), LinalgError> {
    a.validate_finite()?;
    b.validate_finite()?;
    gemm_into_finite_inputs(a, op_a, b, op_b, output)
}

// Internal boundary: callers have checked both inputs since their last mutation.
// Dimensions and output finiteness are still checked here.
#[cfg(feature = "mkl")]
fn gemm_into_finite_inputs(
    a: &ComplexMatrix,
    op_a: MatrixOp,
    b: &ComplexMatrix,
    op_b: MatrixOp,
    output: &mut ComplexMatrix,
) -> Result<(), LinalgError> {
    let adj_a = matches!(op_a, MatrixOp::Adjoint);
    let adj_b = matches!(op_b, MatrixOp::Adjoint);
    let (m, k) = if adj_a {
        (a.cols, a.rows)
    } else {
        (a.rows, a.cols)
    };
    let (l, n) = if adj_b {
        (b.cols, b.rows)
    } else {
        (b.rows, b.cols)
    };
    if k != l || output.rows != m || output.cols != n {
        return Err(LinalgError::Dimension);
    }
    let ma =
        dft_mkl::Matrix::new(a.as_slice(), a.rows, a.cols).map_err(|_| LinalgError::Dimension)?;
    let mb =
        dft_mkl::Matrix::new(b.as_slice(), b.rows, b.cols).map_err(|_| LinalgError::Dimension)?;
    dft_mkl::gemm(ma, adj_a, mb, adj_b, &mut output.data).map_err(|_| LinalgError::Dimension)?;
    output.validate_finite()
}

/// Allocation-free adapter into an existing column-major output. Uses sequential
/// faer multiplication and borrowed views; operands must not alias the output.
#[cfg(not(feature = "mkl"))]
pub fn gemm_into(
    a: &ComplexMatrix,
    op_a: MatrixOp,
    b: &ComplexMatrix,
    op_b: MatrixOp,
    output: &mut ComplexMatrix,
) -> Result<(), LinalgError> {
    a.validate_finite()?;
    b.validate_finite()?;
    let a = faer::MatRef::from_column_major_slice(&a.data, a.rows, a.cols);
    let b = faer::MatRef::from_column_major_slice(&b.data, b.rows, b.cols);
    let (m, k) = match op_a {
        MatrixOp::Identity => (a.nrows(), a.ncols()),
        MatrixOp::Adjoint => (a.ncols(), a.nrows()),
    };
    let (l, n) = match op_b {
        MatrixOp::Identity => (b.nrows(), b.ncols()),
        MatrixOp::Adjoint => (b.ncols(), b.nrows()),
    };
    if k != l || output.rows != m || output.cols != n {
        return Err(LinalgError::Dimension);
    }
    let mut c = faer::MatMut::from_column_major_slice_mut(&mut output.data, m, n);
    let one = Complex64::new(1.0, 0.0);
    match (op_a, op_b) {
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
    output.validate_finite()
}

/// Fixed-size sequential Hermitian eigensolve with persistent backend scratch.
/// Construction allocates; `solve` borrows all inputs and outputs.
pub struct HermitianWorkspace {
    dimension: usize,
    diagonal: faer::diag::Diag<Complex64>,
    scratch: faer::dyn_stack::MemBuffer,
}
impl HermitianWorkspace {
    pub fn new(dimension: usize) -> Result<Self, LinalgError> {
        // Validate sizes before asking the backend to size its scratch.
        dimension
            .checked_mul(dimension)
            .filter(|n| *n <= isize::MAX as usize / size_of::<Complex64>())
            .ok_or(LinalgError::Dimension)?;
        // The backend scratch query assumes n >= 1; preserve empty-matrix
        // semantics without calling it for zero dimensions.
        let scratch = if dimension == 0 {
            faer::dyn_stack::StackReq::EMPTY
        } else {
            faer::linalg::evd::self_adjoint_evd_scratch::<Complex64>(
                dimension,
                faer::linalg::evd::ComputeEigenvectors::Yes,
                Par::Seq,
                Default::default(),
            )
        };
        Ok(Self {
            dimension,
            diagonal: faer::diag::Diag::zeros(dimension),
            scratch: faer::dyn_stack::MemBuffer::new(scratch),
        })
    }
    pub fn solve(
        &mut self,
        a: &ComplexMatrix,
        values: &mut [f64],
        vectors: &mut ComplexMatrix,
    ) -> Result<(), LinalgError> {
        let n = self.dimension;
        if a.rows != n || a.cols != n || values.len() != n || vectors.rows != n || vectors.cols != n
        {
            return Err(LinalgError::Dimension);
        }
        self.solve_prefix(a, values, vectors)
    }

    /// Variable-size solve into storage reserved for this workspace's maximum
    /// dimension. Shapes change without allocation; the fixed-size API retains
    /// its dimension checks.
    pub fn solve_into(
        &mut self,
        a: &ComplexMatrix,
        values: &mut Vec<f64>,
        vectors: &mut ComplexMatrix,
    ) -> Result<(), LinalgError> {
        let n = a.rows;
        if a.cols != n
            || n > self.dimension
            || values.capacity() < n
            || vectors.data.capacity() < n * n
        {
            return Err(LinalgError::Dimension);
        }
        vectors.data.resize(n * n, Complex64::default());
        vectors.rows = n;
        vectors.cols = n;
        values.resize(n, 0.0);
        self.solve_prefix(a, values, vectors)
    }

    fn solve_prefix(
        &mut self,
        a: &ComplexMatrix,
        values: &mut [f64],
        vectors: &mut ComplexMatrix,
    ) -> Result<(), LinalgError> {
        let n = a.rows;
        a.validate_finite()?;
        let scale = a.data.iter().map(|z| z.norm()).fold(1.0, f64::max);
        let allowed = HERMITIAN_TOLERANCE * scale;
        let mut defect: f64 = 0.0;
        for j in 0..n {
            for i in 0..=j {
                defect = defect.max((a[(i, j)] - a[(j, i)].conj()).norm());
            }
        }
        if !scale.is_finite() || !defect.is_finite() {
            return Err(LinalgError::NonFinite);
        }
        if defect > allowed {
            return Err(LinalgError::NotHermitian { defect, allowed });
        }
        if n == 0 {
            return Ok(());
        }
        faer::linalg::evd::self_adjoint_evd(
            faer::MatRef::from_column_major_slice(&a.data, n, n),
            self.diagonal
                .as_mut()
                .column_vector_mut()
                .subrows_mut(0, n)
                .as_diagonal_mut(),
            Some(faer::MatMut::from_column_major_slice_mut(
                &mut vectors.data,
                n,
                n,
            )),
            Par::Seq,
            faer::dyn_stack::MemStack::new(&mut self.scratch),
            Default::default(),
        )
        .map_err(|_| LinalgError::Decomposition)?;
        for (i, value) in values.iter_mut().enumerate() {
            *value = self.diagonal[i].re;
        }
        if values.iter().any(|v| !v.is_finite()) {
            return Err(LinalgError::NonFinite);
        }
        vectors.validate_finite()
    }
}

/// Two-pass MGS into preallocated storage, returning retained rank. Trailing
/// columns are zeroed. `scratch` has exactly block.rows() entries. Same relative
/// rank policy and projection arithmetic as `orthonormalize`.
pub fn orthonormalize_into(
    block: &ComplexMatrix,
    output: &mut ComplexMatrix,
    scratch: &mut [Complex64],
) -> Result<usize, LinalgError> {
    if block.rows != output.rows || block.cols != output.cols || scratch.len() != block.rows {
        return Err(LinalgError::Dimension);
    }
    block.validate_finite()?;
    let mut count = 0;
    for j in 0..block.cols {
        scratch.copy_from_slice(block.column(j));
        let original = norm(scratch);
        if !original.is_finite() {
            return Err(LinalgError::NonFinite);
        }
        if original == 0.0 {
            continue;
        }
        for z in scratch.iter_mut() {
            *z /= original;
        }
        for _ in 0..2 {
            for k in 0..count {
                subtract_projection(scratch, output.column(k));
            }
        }
        let remaining = norm(scratch);
        if !remaining.is_finite() {
            return Err(LinalgError::NonFinite);
        }
        if remaining > DEPENDENCE_TOLERANCE {
            for (z, v) in output.column_mut(count).iter_mut().zip(scratch.iter()) {
                *z = *v / remaining;
            }
            count += 1;
        }
    }
    for j in count..output.cols {
        output.column_mut(j).fill(Complex64::default());
    }
    Ok(count)
}

pub mod accounting;
mod tiled;
pub use tiled::{GemmWorkspace, ProjectionWorkspace};

#[cfg(test)]
mod column_storage_tests {
    use super::*;
    #[test]
    fn live_column_ranges_preserve_storage_and_reject_overflow() {
        let mut a = ComplexMatrix::zeros(3, 1).unwrap();
        a[(1, 0)] = Complex64::new(2.0, -3.0);
        a.reserve_columns(4).unwrap();
        let address = a.as_slice().as_ptr();
        let b = a.clone();
        a.append_columns(&b).unwrap();
        assert_eq!(a.column(0), a.column(1));
        a.resize_columns(1).unwrap();
        a.resize_columns(4).unwrap();
        assert_eq!(a.column(0), b.column(0));
        assert!(a.as_slice()[3..].iter().all(|z| *z == Complex64::default()));
        assert_eq!(a.as_slice().as_ptr(), address);
        let before = a.clone();
        assert_eq!(a.resize_columns(5), Err(LinalgError::Dimension));
        assert_eq!(a.append_columns(&b), Err(LinalgError::Dimension));
        assert_eq!(a, before);
        a.replace_columns(&b).unwrap();
        assert_eq!(a, b);
        assert_eq!(a.as_slice().as_ptr(), address);
    }
}
