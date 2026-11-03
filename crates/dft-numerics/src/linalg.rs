//! Dense backend boundary; all matrices are column-major and backend independent.
use faer::{Accum, Mat, Par, linalg::matmul::matmul};
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
    pub fn rows(&self) -> usize {
        self.rows
    }
    pub fn cols(&self) -> usize {
        self.cols
    }
    pub fn as_slice(&self) -> &[Complex64] {
        &self.data
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
        if self
            .data
            .iter()
            .all(|z| z.re.is_finite() && z.im.is_finite())
        {
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
    Mat::from_fn(a.rows, a.cols, |i, j| a[(i, j)])
}
fn owned(a: faer::MatRef<'_, Complex64>) -> Result<ComplexMatrix, LinalgError> {
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
