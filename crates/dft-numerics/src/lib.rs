//! Reusable complex FFTs with explicit plane-wave normalization.
pub mod fused;
pub mod fused_pair;
pub mod fused_pair_plan;
pub use fused_pair::{PairDensityOperator, PairDensityWorkspace, PairMap, PairPlan};
pub use fused_pair_plan::PairFusedPlan;
pub mod linalg;
pub mod real;
pub use real::{
    RealGemmWorkspace, RealMatrix, RealMatrixMut, RealMatrixRef, RealOp, symmetric_eigh,
};
#[cfg(feature = "mkl")]
pub mod mkl_backend;
pub mod pool;
#[cfg(not(feature = "mkl"))]
pub use RustFftPlan3d as FftPlan3d;
#[cfg(feature = "mkl")]
pub use mkl_backend::MklFftPlan3d as FftPlan3d;

use ndarray::Array3;
use num_complex::Complex64;
use rustfft::{Fft, FftPlanner};
use std::sync::Arc;
use thiserror::Error;

#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum FftError {
    #[error("grid dimensions must be nonzero and their storage size must fit isize")]
    InvalidGrid,
    #[cfg(feature = "mkl")]
    #[error("MKL FFT: {0}")]
    Mkl(#[from] dft_mkl::Error),
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

/// Precomputed lane dependencies for a reciprocal-space support, in axis order
/// 0, 1, 2. Indices may be unordered or duplicated; empty support is valid.
/// Construction allocates; transforms reuse these lists without allocation.
pub struct SpherePattern {
    shape: [usize; 3],
    present: Vec<bool>,
    inverse_lanes: [Vec<usize>; 3],
    forward_lanes: [Vec<usize>; 3],
}
impl SpherePattern {
    pub fn from_indices(shape: [usize; 3], indices: &[[usize; 3]]) -> Result<Self, FftError> {
        let n = shape
            .iter()
            .try_fold(1_usize, |n, &d| n.checked_mul(d))
            .ok_or(FftError::InvalidGrid)?;
        if n == 0 || n > isize::MAX as usize / std::mem::size_of::<Complex64>() {
            return Err(FftError::InvalidGrid);
        }
        let [nx, ny, nz] = shape;
        let mut present = vec![false; n];
        let mut jk = vec![false; ny * nz];
        let mut ks = vec![false; nz];
        let mut is = vec![false; nx];
        let mut ij = vec![false; nx * ny];
        for &[i, j, k] in indices {
            if i >= nx || j >= ny || k >= nz {
                return Err(FftError::IncompatibleArray);
            }
            present[(i * ny + j) * nz + k] = true;
            jk[j * nz + k] = true;
            ks[k] = true;
            is[i] = true;
            ij[i * ny + j] = true;
        }
        let inverse_lanes = [
            (0..ny * nz).filter(|&q| jk[q]).collect(),
            (0..nx * nz)
                .filter(|&q| ks[q % nz])
                .map(|q| q / nz * ny * nz + q % nz)
                .collect(),
            (0..nx * ny).map(|q| q * nz).collect(),
        ];
        let forward_lanes = [
            (0..ny * nz).collect(),
            (0..nx * nz)
                .filter(|&q| is[q / nz])
                .map(|q| q / nz * ny * nz + q % nz)
                .collect(),
            (0..nx * ny).filter(|&q| ij[q]).map(|q| q * nz).collect(),
        ];
        Ok(Self {
            shape,
            present,
            inverse_lanes,
            forward_lanes,
        })
    }
}

/// Persistent rustfft plans, shared scratch, and one reusable strided lane.
/// After construction, transforms only copy slices, gather/scatter lanes, scale,
/// and call rustfft's allocation-free `process_with_scratch`. No buffers grow.
/// Scratch requirements are checked before every axis pass in debug builds.
pub struct RustFftPlan3d {
    shape: [usize; 3],
    forward: [Arc<dyn Fft<f64>>; 3],
    inverse: [Arc<dyn Fft<f64>>; 3],
    scratch: Vec<Complex64>,
    lane: Vec<Complex64>,
    size: f64,
    #[cfg(test)]
    pruned_lane_count: usize,
}
impl RustFftPlan3d {
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
            #[cfg(test)]
            pruned_lane_count: 0,
        })
    }
    /// Inverse FFT of input supported on `pattern`. Input outside the pattern
    /// must be canonical complex zero (+0,+0), as produced by grid scatter.
    /// Debug assertions check this while copying, without an extra grid scan.
    /// Axis 0 cannot populate an absent (j,k) lane; axis 1 cannot populate an
    /// absent k plane. Those skipped lanes remain exactly zero before axis 2.
    pub fn to_real_pruned(
        &mut self,
        input: &Array3<Complex64>,
        pattern: &SpherePattern,
        output: &mut Array3<Complex64>,
    ) -> Result<(), FftError> {
        self.transform_pruned(input, pattern, output, true)
    }

    /// Forward FFT, bit-identical to the full transform on pattern entries.
    /// Entries outside the pattern are set to zero, NOT the full FFT values.
    /// Axis 1 does not mix i, and axis 2 does not mix (i,j), so omitted
    /// lanes cannot contribute to a requested coefficient.
    pub fn to_reciprocal_pruned(
        &mut self,
        input: &Array3<Complex64>,
        pattern: &SpherePattern,
        output: &mut Array3<Complex64>,
    ) -> Result<(), FftError> {
        self.transform_pruned(input, pattern, output, false)
    }

    fn transform_pruned(
        &mut self,
        input: &Array3<Complex64>,
        pattern: &SpherePattern,
        output: &mut Array3<Complex64>,
        inverse: bool,
    ) -> Result<(), FftError> {
        if pattern.shape != self.shape
            || input.shape() != self.shape
            || output.shape() != self.shape
            || !input.is_standard_layout()
            || !output.is_standard_layout()
        {
            return Err(FftError::IncompatibleArray);
        }
        let data = output.as_slice_mut().unwrap();
        for (q, (dst, src)) in data.iter_mut().zip(input.as_slice().unwrap()).enumerate() {
            debug_assert!(
                !inverse || pattern.present[q] || (src.re.to_bits() == 0 && src.im.to_bits() == 0),
                "inverse input outside pattern must be canonical zero"
            );
            *dst = *src;
        }
        let (plans, lanes) = if inverse {
            (&self.inverse, &pattern.inverse_lanes)
        } else {
            (&self.forward, &pattern.forward_lanes)
        };
        for (axis, plan) in plans.iter().enumerate() {
            let len = self.shape[axis];
            let stride: usize = self.shape[axis + 1..].iter().product();
            let scale = 1.0 / len as f64;
            debug_assert!(self.scratch.len() >= plan.get_inplace_scratch_len());
            for &start in &lanes[axis] {
                let lane = if stride == 1 {
                    &mut data[start..start + len]
                } else {
                    for (i, z) in self.lane[..len].iter_mut().enumerate() {
                        *z = data[start + i * stride];
                    }
                    &mut self.lane[..len]
                };
                plan.process_with_scratch(lane, &mut self.scratch);
                #[cfg(test)]
                {
                    self.pruned_lane_count += 1;
                }
                if inverse {
                    for z in lane.iter_mut() {
                        *z *= scale;
                    }
                }
                if stride != 1 {
                    for (i, z) in self.lane[..len].iter().enumerate() {
                        data[start + i * stride] = *z;
                    }
                }
            }
        }
        for (q, z) in data.iter_mut().enumerate() {
            if inverse {
                *z *= self.size;
            } else if pattern.present[q] {
                *z /= self.size;
            } else {
                *z = Complex64::default();
            }
        }
        Ok(())
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
impl Fft3d for RustFftPlan3d {
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
mod tests {
    // Existing tests include implementation-specific lane-count assertions.
    // They always exercise Rust; MKL contracts have their own tests.
    use super::RustFftPlan3d as FftPlan3d;
    include!("tests.rs");
}
