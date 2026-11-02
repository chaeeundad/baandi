//! Backend-independent dense and matrix-free Hermitian eigensolvers.
use dft_numerics::linalg::{
    ComplexMatrix, LinalgError,
    MatrixOp::{Adjoint, Identity},
    gemm, hermitian_eigh, norm, orthonormalize,
};
use num_complex::Complex64;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum SolverError {
    #[error(transparent)]
    Linalg(#[from] LinalgError),
    #[error("invalid solver configuration: {0}")]
    InvalidParameter(&'static str),
    #[error("operator/preconditioner failure: {0}")]
    Application(String),
    #[error("dense dimension {dimension} exceeds limit {limit}")]
    DenseLimit { dimension: usize, limit: usize },
    #[error("not converged after {iterations} iterations; max residual {max_residual:e}")]
    NotConverged {
        iterations: usize,
        max_residual: f64,
        residuals: Vec<f64>,
    },
}

/// Hermitian by contract. Implementations must overwrite all of y, without changing x.
/// Both slices have length dimension(); errors must be returned, not hidden.
pub trait LinearOperator {
    fn dimension(&self) -> usize;
    fn apply(&self, x: &[Complex64], y: &mut [Complex64]) -> Result<(), SolverError>;
    fn apply_block(&self, x: &ComplexMatrix, y: &mut ComplexMatrix) -> Result<(), SolverError> {
        if x.rows() != self.dimension() || y.rows() != x.rows() || y.cols() != x.cols() {
            return Err(LinalgError::Dimension.into());
        }
        for j in 0..x.cols() {
            self.apply(x.column(j), y.column_mut(j))?;
        }
        y.validate_finite()?;
        Ok(())
    }
}
/// Transform residual columns, using the matching Ritz values; overwrite output.
pub trait Preconditioner {
    fn apply(
        &self,
        eigenvalues: &[f64],
        residuals: &ComplexMatrix,
        output: &mut ComplexMatrix,
    ) -> Result<(), SolverError>;
}

/// Correction -r/(d-lambda), with sign-preserving denominator floor (zero takes +).
#[derive(Debug, Clone)]
pub struct DiagonalPreconditioner {
    diagonal: Vec<f64>,
    denominator_floor: f64,
}
impl DiagonalPreconditioner {
    /// The positive floor has the same units as the diagonal/eigenvalues.
    pub fn new(diagonal: Vec<f64>, denominator_floor: f64) -> Result<Self, SolverError> {
        if diagonal.iter().any(|d| !d.is_finite())
            || !denominator_floor.is_finite()
            || denominator_floor <= 0.0
        {
            return Err(SolverError::InvalidParameter(
                "finite diagonal and positive finite denominator floor required",
            ));
        }
        Ok(Self {
            diagonal,
            denominator_floor,
        })
    }
}
impl Preconditioner for DiagonalPreconditioner {
    fn apply(
        &self,
        eigenvalues: &[f64],
        residuals: &ComplexMatrix,
        output: &mut ComplexMatrix,
    ) -> Result<(), SolverError> {
        if residuals.rows() != self.diagonal.len()
            || eigenvalues.len() != residuals.cols()
            || output.rows() != residuals.rows()
            || output.cols() != residuals.cols()
        {
            return Err(LinalgError::Dimension.into());
        }
        residuals.validate_finite()?;
        if eigenvalues.iter().any(|x| !x.is_finite()) {
            return Err(LinalgError::NonFinite.into());
        }
        for (j, &lambda) in eigenvalues.iter().enumerate() {
            for (i, &d) in self.diagonal.iter().enumerate() {
                let delta = d - lambda;
                if !delta.is_finite() {
                    return Err(LinalgError::NonFinite.into());
                }
                let denom = if delta < 0.0 {
                    -delta.abs().max(self.denominator_floor)
                } else {
                    delta.max(self.denominator_floor)
                };
                output[(i, j)] = -residuals[(i, j)] / denom;
            }
        }
        output.validate_finite()?;
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct EigenSolution {
    pub eigenvalues: Vec<f64>,
    pub eigenvectors: ComplexMatrix,
    pub residuals: Vec<f64>,
    pub iterations: usize,
    /// Maximum absolute residual 2-norm at each Rayleigh-Ritz iteration.
    pub history: Vec<f64>,
}

/// Oracle only: rejects dimensions above 1024 before allocating or applying H.
#[derive(Debug, Default, Clone, Copy)]
pub struct DenseEigensolver;
impl DenseEigensolver {
    pub const MAX_DIMENSION: usize = 1024;
    pub fn solve(
        &self,
        operator: &dyn LinearOperator,
        bands: usize,
    ) -> Result<EigenSolution, SolverError> {
        let n = operator.dimension();
        if n > Self::MAX_DIMENSION {
            return Err(SolverError::DenseLimit {
                dimension: n,
                limit: Self::MAX_DIMENSION,
            });
        }
        validate_bands(n, bands)?;
        let mut unit = ComplexMatrix::zeros(n, n)?;
        for i in 0..n {
            unit[(i, i)] = Complex64::new(1.0, 0.0);
        }
        let matrix = apply(operator, &unit)?;
        let (mut values, vectors) = hermitian_eigh(&matrix)?;
        values.truncate(bands);
        let vectors = first_columns(&vectors, bands)?;
        let (_, residuals) = residual(operator, &vectors, &values)?;
        let maximum = max_residual(&residuals);
        Ok(EigenSolution {
            eigenvalues: values,
            eigenvectors: vectors,
            residuals,
            iterations: 1,
            history: vec![maximum],
        })
    }
}

#[derive(Debug, Clone)]
pub struct Davidson {
    pub bands: usize,
    /// Absolute residual norm, in the operator's eigenvalue units.
    pub tolerance: f64,
    pub max_iterations: usize,
    pub max_subspace: usize,
    /// Maximum number of residual corrections appended per iteration.
    pub block_size: usize,
    pub seed: u64,
}
impl Davidson {
    pub fn solve(
        &self,
        operator: &dyn LinearOperator,
        initial: Option<&ComplexMatrix>,
        preconditioner: Option<&dyn Preconditioner>,
    ) -> Result<EigenSolution, SolverError> {
        let n = operator.dimension();
        validate_bands(n, self.bands)?;
        let capacity = self.max_subspace.min(n);
        if !self.tolerance.is_finite()
            || self.tolerance <= 0.0
            || self.max_iterations == 0
            || self.block_size == 0
            || capacity < self.bands
            || (capacity == self.bands && self.bands < n)
        {
            return Err(SolverError::InvalidParameter(
                "positive tolerance/iterations/block size and subspace larger than bands required (unless bands=n)",
            ));
        }
        let mut rng = SplitMix64(self.seed);
        let mut basis = if let Some(initial) = initial {
            if initial.rows() != n || initial.cols() > capacity {
                return Err(LinalgError::Dimension.into());
            }
            orthonormalize(initial, None)?.vectors
        } else {
            ComplexMatrix::zeros(n, 0)?
        };
        let initial_width = self.bands.max(self.block_size).min(capacity);
        // Fill missing rank with independent seeded columns, including deficient user guesses.
        while basis.cols() < initial_width {
            let random = rng.block(n, initial_width - basis.cols())?;
            let extra = orthonormalize(&random, Some(&basis))?.vectors;
            if extra.cols() == 0 {
                return Err(SolverError::InvalidParameter(
                    "unable to complete initial basis",
                ));
            }
            basis = concatenate(&basis, &extra)?;
        }
        let mut history = Vec::new();
        for iteration in 1..=self.max_iterations {
            let h_basis = apply(operator, &basis)?;
            let projected = gemm(&basis, Adjoint, &h_basis, Identity)?;
            let (mut values, coefficients) = hermitian_eigh(&projected)?;
            values.truncate(self.bands);
            let vectors = gemm(
                &basis,
                Identity,
                &first_columns(&coefficients, self.bands)?,
                Identity,
            )?;
            // Reapply H: convergence is based on actual operator residuals, not a recurrence.
            let (residual_block, residuals) = residual(operator, &vectors, &values)?;
            let maximum = max_residual(&residuals);
            history.push(maximum);
            if residuals.iter().all(|&r| r <= self.tolerance) {
                return Ok(EigenSolution {
                    eigenvalues: values,
                    eigenvectors: vectors,
                    residuals,
                    iterations: iteration,
                    history,
                });
            }
            if iteration == self.max_iterations {
                return Err(SolverError::NotConverged {
                    iterations: iteration,
                    max_residual: maximum,
                    residuals,
                });
            }
            let selected: Vec<_> = residuals
                .iter()
                .enumerate()
                .filter(|(_, r)| **r > self.tolerance)
                .map(|(j, _)| j)
                .take(self.block_size.min(capacity - self.bands))
                .collect();
            let mut corrections = ComplexMatrix::zeros(n, selected.len())?;
            let mut shifts = Vec::new();
            for (k, &j) in selected.iter().enumerate() {
                corrections
                    .column_mut(k)
                    .copy_from_slice(residual_block.column(j));
                shifts.push(values[j]);
            }
            if let Some(p) = preconditioner {
                let mut output = ComplexMatrix::zeros(n, corrections.cols())?;
                p.apply(&shifts, &corrections, &mut output)?;
                if output.rows() != n || output.cols() != corrections.cols() {
                    return Err(LinalgError::Dimension.into());
                }
                output.validate_finite()?;
                corrections = output;
            }
            if basis.cols() + corrections.cols() > capacity {
                basis = orthonormalize(&vectors, None)?.vectors;
            }
            let mut extra = orthonormalize(&corrections, Some(&basis))?.vectors;
            // A stalled direction can occur with an exact preconditioner or an invariant
            // initial space. Try fresh independent directions, never infer convergence.
            if extra.cols() == 0 {
                let width = self.block_size.min(capacity - basis.cols());
                let random = rng.block(n, width)?;
                extra = orthonormalize(&random, Some(&basis))?.vectors;
            }
            if extra.cols() == 0 {
                return Err(SolverError::NotConverged {
                    iterations: iteration,
                    max_residual: maximum,
                    residuals,
                });
            }
            basis = concatenate(&basis, &extra)?;
        }
        unreachable!("positive iteration limit validated above")
    }
}
fn validate_bands(n: usize, bands: usize) -> Result<(), SolverError> {
    if bands == 0 || bands > n {
        Err(SolverError::InvalidParameter(
            "require 1 <= bands <= dimension",
        ))
    } else {
        Ok(())
    }
}
fn apply(op: &dyn LinearOperator, x: &ComplexMatrix) -> Result<ComplexMatrix, SolverError> {
    let mut y = ComplexMatrix::zeros(op.dimension(), x.cols())?;
    op.apply_block(x, &mut y)?;
    if y.rows() != op.dimension() || y.cols() != x.cols() {
        return Err(LinalgError::Dimension.into());
    }
    y.validate_finite()?;
    Ok(y)
}
fn residual(
    op: &dyn LinearOperator,
    vectors: &ComplexMatrix,
    values: &[f64],
) -> Result<(ComplexMatrix, Vec<f64>), SolverError> {
    let mut r = apply(op, vectors)?;
    let mut norms = Vec::new();
    for (j, &lambda) in values.iter().enumerate() {
        for i in 0..r.rows() {
            r[(i, j)] -= vectors[(i, j)] * lambda;
        }
        let magnitude = norm(r.column(j));
        if !magnitude.is_finite() {
            return Err(LinalgError::NonFinite.into());
        }
        norms.push(magnitude);
    }
    Ok((r, norms))
}
fn max_residual(r: &[f64]) -> f64 {
    r.iter().copied().fold(0.0, f64::max)
}
fn first_columns(a: &ComplexMatrix, count: usize) -> Result<ComplexMatrix, LinalgError> {
    ComplexMatrix::new(a.rows(), count, a.as_slice()[..a.rows() * count].to_vec())
}
fn concatenate(a: &ComplexMatrix, b: &ComplexMatrix) -> Result<ComplexMatrix, LinalgError> {
    if a.rows() != b.rows() {
        return Err(LinalgError::Dimension);
    }
    let mut data = a.as_slice().to_vec();
    data.extend_from_slice(b.as_slice());
    ComplexMatrix::new(a.rows(), a.cols() + b.cols(), data)
}
/// SplitMix64, deterministic and local: no process-global random state.
struct SplitMix64(u64);
impl SplitMix64 {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        ((z ^ (z >> 31)) >> 11) as f64 * (1.0 / ((1_u64 << 53) as f64)) - 0.5
    }
    fn block(&mut self, rows: usize, cols: usize) -> Result<ComplexMatrix, LinalgError> {
        let mut block = ComplexMatrix::zeros(rows, cols)?;
        for j in 0..cols {
            for i in 0..rows {
                block[(i, j)] = Complex64::new(self.next(), self.next());
            }
        }
        Ok(block)
    }
}
#[cfg(test)]
mod tests;
