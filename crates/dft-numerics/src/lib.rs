//! Reusable complex FFTs with explicit plane-wave normalization.
pub mod linalg;

use ndarray::Array3;
use num_complex::Complex64;
use rustfft::{Fft, FftPlanner};
use std::sync::Arc;
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

/// Persistent rustfft plans, shared scratch, and one reusable strided lane.
/// After construction, transforms only copy slices, gather/scatter lanes, scale,
/// and call rustfft's allocation-free `process_with_scratch`. No buffers grow.
/// Scratch requirements are checked before every axis pass in debug builds.
pub struct FftPlan3d {
    shape: [usize; 3],
    forward: [Arc<dyn Fft<f64>>; 3],
    inverse: [Arc<dyn Fft<f64>>; 3],
    scratch: Vec<Complex64>,
    lane: Vec<Complex64>,
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
        let mut planner = FftPlanner::new();
        let forward = shape.map(|n| planner.plan_fft_forward(n));
        let inverse = shape.map(|n| planner.plan_fft_inverse(n));
        let scratch_len = forward
            .iter()
            .chain(&inverse)
            .map(|plan| plan.get_inplace_scratch_len())
            .max()
            .unwrap();
        Ok(Self {
            shape,
            forward,
            inverse,
            scratch: vec![Complex64::default(); scratch_len],
            lane: vec![Complex64::default(); *shape.iter().max().unwrap()],
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
        let data = output.as_slice_mut().unwrap();
        data.copy_from_slice(input.as_slice().unwrap());
        let plans = if inverse {
            &self.inverse
        } else {
            &self.forward
        };
        for (axis, plan) in plans.iter().enumerate() {
            let len = self.shape[axis];
            let stride: usize = self.shape[axis + 1..].iter().product();
            debug_assert!(self.scratch.len() >= plan.get_inplace_scratch_len());
            debug_assert!(self.lane.len() >= len);
            let scale = 1.0 / len as f64;
            for block in data.chunks_exact_mut(len * stride) {
                for offset in 0..stride {
                    // The last axis is already contiguous; other axes reuse one lane.
                    let lane = if stride == 1 {
                        &mut *block
                    } else {
                        for (i, z) in self.lane[..len].iter_mut().enumerate() {
                            *z = block[i * stride + offset];
                        }
                        &mut self.lane[..len]
                    };
                    plan.process_with_scratch(lane, &mut self.scratch);
                    if inverse {
                        // Match ndrustfft's per-axis multiplication exactly.
                        for z in lane.iter_mut() {
                            *z *= scale;
                        }
                    }
                    if stride != 1 {
                        for (i, z) in self.lane[..len].iter().enumerate() {
                            block[i * stride + offset] = *z;
                        }
                    }
                }
            }
        }
        for z in data {
            if inverse {
                *z *= self.size;
            } else {
                *z /= self.size;
            }
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
