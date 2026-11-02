//! Reusable complex FFTs with explicit plane-wave normalization.
pub mod linalg;

use ndarray::Array3;
use ndrustfft::{FftHandler, ndfft, ndifft};
use num_complex::Complex64;
use thiserror::Error;

#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum FftError {
    #[error("grid dimensions must be nonzero and their storage size must fit isize")]
    InvalidGrid,
    #[error("input and output must have the plan's shape and standard C layout")]
    IncompatibleArray,
}

/// Narrow backend boundary. Arrays use (i,j,k), with k contiguous.
pub trait Fft3d {
    fn shape(&self) -> [usize; 3];
    /// f(r) = sum_G f(G) exp(+i G.r), without scaling.
    fn to_real(
        &mut self,
        input: &Array3<Complex64>,
        output: &mut Array3<Complex64>,
    ) -> Result<(), FftError>;
    /// f(G) = sum_r f(r) exp(-i G.r) / N.
    fn to_reciprocal(
        &mut self,
        input: &Array3<Complex64>,
        output: &mut Array3<Complex64>,
    ) -> Result<(), FftError>;
}

/// Three persistent handlers and two persistent full-grid scratch buffers.
/// ndrustfft may still allocate internal lane scratch; this is not allocation-free.
pub struct FftPlan3d {
    shape: [usize; 3],
    handlers: [FftHandler<f64>; 3],
    first: Array3<Complex64>,
    second: Array3<Complex64>,
    size: f64,
}
impl FftPlan3d {
    pub fn new(shape: [usize; 3]) -> Result<Self, FftError> {
        let n = shape
            .iter()
            .try_fold(1_usize, |n, &d| n.checked_mul(d))
            .ok_or(FftError::InvalidGrid)?;
        if n == 0 || n > isize::MAX as usize / std::mem::size_of::<Complex64>() {
            return Err(FftError::InvalidGrid);
        }
        let dims = (shape[0], shape[1], shape[2]);
        Ok(Self {
            shape,
            handlers: shape.map(FftHandler::new),
            first: Array3::zeros(dims),
            second: Array3::zeros(dims),
            size: n as f64,
        })
    }
    fn transform(
        &mut self,
        input: &Array3<Complex64>,
        output: &mut Array3<Complex64>,
        inverse: bool,
    ) -> Result<(), FftError> {
        if input.shape() != self.shape
            || output.shape() != self.shape
            || !input.is_standard_layout()
            || !output.is_standard_layout()
        {
            return Err(FftError::IncompatibleArray);
        }
        // Measured ndrustfft default: forward unscaled, inverse divided by N.
        if inverse {
            ndifft(input, &mut self.first, &self.handlers[0], 0);
            ndifft(&self.first, &mut self.second, &self.handlers[1], 1);
            ndifft(&self.second, output, &self.handlers[2], 2);
            output.mapv_inplace(|z| z * self.size);
        } else {
            ndfft(input, &mut self.first, &self.handlers[0], 0);
            ndfft(&self.first, &mut self.second, &self.handlers[1], 1);
            ndfft(&self.second, output, &self.handlers[2], 2);
            output.mapv_inplace(|z| z / self.size);
        }
        Ok(())
    }
}
impl Fft3d for FftPlan3d {
    fn shape(&self) -> [usize; 3] {
        self.shape
    }
    fn to_real(
        &mut self,
        input: &Array3<Complex64>,
        output: &mut Array3<Complex64>,
    ) -> Result<(), FftError> {
        self.transform(input, output, true)
    }
    fn to_reciprocal(
        &mut self,
        input: &Array3<Complex64>,
        output: &mut Array3<Complex64>,
    ) -> Result<(), FftError> {
        self.transform(input, output, false)
    }
}

#[cfg(test)]
mod tests;
