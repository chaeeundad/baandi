//! Sphere-to-sphere local operators with reusable, worker-owned storage.
//!
//! The ordered implementation preserves the historical x/y/z arithmetic in
//! both directions. The plane implementation instead pairs z/y/x with x/y/z;
//! it is an explicit experiment, never the default. See batched-fft-study.md.
use crate::{FftError, SpherePattern};
use num_complex::Complex64 as C;
use rustfft::{Fft, FftPlanner};
use std::sync::Arc;

/// Backend boundary for periodic Bloch coefficients (the k phase is absent).
/// Plans are constructed from grid dimensions and *ordered* sphere indices.
/// Grid slices use C order: `(x * ny + y) * nz + z`. Repeated indices have
/// scatter's last-write-wins semantics; outputs retain the supplied ordering.
///
/// Implementations own descriptors and scratch, allocate only at construction,
/// and use no internal threads. One plan belongs to one outer band/k worker.
/// Vendor implementations may batch 1D sticks with arbitrary stride/distance,
/// prune plane columns, and fuse multiplication without exposing those layouts.
/// Normalization is inverse unscaled, forward divided by grid size. Backends
/// must state their arithmetic order and bitwise/normwise equivalence contract.
/// Slice/volume errors are detected before changing outputs. Finite coefficient
/// and potential validation belongs to the physics caller, outside this kernel.
pub trait FusedOperator {
    fn shape(&self) -> [usize; 3];
    fn coefficient_count(&self) -> usize;
    /// Write P V psi in the input sphere order. Potential values are in Ha;
    /// kinetic and nonlocal terms are not included.
    fn apply_local(&mut self, c: &[C], potential: &[f64], out: &mut [C]) -> Result<(), FftError>;
    /// Add `weight * |psi|^2 / volume` in that arithmetic order, in e/Bohr^3.
    /// Weight is separate from volume (Bohr^3) to avoid changing rounding.
    /// Exactly zero weight leaves the accumulator untouched and skips FFTs.
    fn accumulate_density(
        &mut self,
        c: &[C],
        weight: f64,
        volume: f64,
        density: &mut [f64],
    ) -> Result<(), FftError>;
}

struct Kernels {
    forward: [Arc<dyn Fft<f64>>; 3],
    inverse: [Arc<dyn Fft<f64>>; 3],
    scratch: Vec<C>,
    lane: Vec<C>,
}
impl Kernels {
    fn new(shape: [usize; 3]) -> Self {
        let mut planner = FftPlanner::new();
        let forward = shape.map(|n| planner.plan_fft_forward(n));
        let inverse = shape.map(|n| planner.plan_fft_inverse(n));
        let scratch = forward
            .iter()
            .chain(&inverse)
            .map(|p| p.get_inplace_scratch_len())
            .max()
            .unwrap();
        Self {
            forward,
            inverse,
            scratch: vec![C::default(); scratch],
            lane: vec![C::default(); *shape.iter().max().unwrap()],
        }
    }
    fn lane(&mut self, axis: usize, inverse: bool, n: usize) {
        let plan = if inverse {
            &self.inverse[axis]
        } else {
            &self.forward[axis]
        };
        debug_assert!(self.scratch.len() >= plan.get_inplace_scratch_len());
        plan.process_with_scratch(&mut self.lane[..n], &mut self.scratch);
        if inverse {
            for z in &mut self.lane[..n] {
                *z *= 1.0 / n as f64;
            }
        }
    }
}

/// Same x/y/z order and per-axis inverse scaling as `crate::RustFftPlan3d`.
/// Scatter uses compact x sticks; inverse y results use only active z planes.
/// The final inverse z pass fuses scaling with multiplication or density.
/// For each real-space y, an xz plane holds the inverse z outputs, multiplies
/// V, and transforms x forward. Only requested x columns are gathered into
/// the compact forward intermediate, followed by forward y/z. This asymmetric
/// fusion preserves both historical axis orders without forming a full real
/// grid. There is no full-grid zero-fill, scatter or transform copy. Buffers
/// become dense only when the support itself occupies every column/plane.
pub struct OrderedFusedPlan {
    shape: [usize; 3],
    size: usize,
    kernels: Kernels,
    indices: Vec<[usize; 3]>,
    scatter: Vec<usize>,
    sticks: Vec<C>,
    // For each active z plane, one optional x stick for every y coordinate.
    planes: Vec<usize>,
    plane_sticks: Vec<usize>,
    partial: Vec<C>,
    active_x: Vec<usize>,
    forward_z: Vec<usize>,
    gather: Vec<usize>,
    forward: Vec<C>,
    real_plane: Vec<C>,
}
impl OrderedFusedPlan {
    pub fn new(shape: [usize; 3], indices: &[[usize; 3]]) -> Result<Self, FftError> {
        let pattern = SpherePattern::from_indices(shape, indices)?;
        let [nx, ny, nz] = shape;
        let mut slots = vec![usize::MAX; ny * nz];
        for (slot, &jk) in pattern.inverse_lanes[0].iter().enumerate() {
            slots[jk] = slot;
        }
        let scatter = indices
            .iter()
            .map(|&[x, y, z]| slots[y * nz + z] * nx + x)
            .collect();
        let planes: Vec<_> = (0..nz)
            .filter(|&z| (0..ny).any(|y| slots[y * nz + z] != usize::MAX))
            .collect();
        let plane_sticks = planes
            .iter()
            .flat_map(|&z| (0..ny).map(move |y| y * nz + z))
            .map(|jk| slots[jk])
            .collect();
        let partial = vec![C::default(); nx * ny * planes.len()];
        let sticks = vec![C::default(); nx * pattern.inverse_lanes[0].len()];
        let active_x: Vec<_> = (0..nx)
            .filter(|&x| indices.iter().any(|p| p[0] == x))
            .collect();
        let mut x_slots = vec![usize::MAX; nx];
        for (slot, &x) in active_x.iter().enumerate() {
            x_slots[x] = slot;
        }
        let forward_z = pattern.forward_lanes[2]
            .iter()
            .map(|&start| {
                let xy = start / nz;
                (x_slots[xy / ny] * ny + xy % ny) * nz
            })
            .collect();
        let gather = indices
            .iter()
            .map(|&[x, y, z]| (x_slots[x] * ny + y) * nz + z)
            .collect();
        let forward = vec![C::default(); active_x.len() * ny * nz];
        Ok(Self {
            shape,
            size: nx * ny * nz,
            kernels: Kernels::new(shape),
            indices: indices.to_vec(),
            scatter,
            sticks,
            planes,
            plane_sticks,
            partial,
            active_x,
            forward_z,
            gather,
            forward,
            real_plane: vec![C::default(); nx * nz],
        })
    }
    fn inverse_xy(&mut self, c: &[C]) {
        let [nx, ny, _] = self.shape;
        self.sticks.fill(C::default());
        for (&slot, &value) in self.scatter.iter().zip(c) {
            self.sticks[slot] = value;
        }
        // Preserve the oracle's per-lane call (batch kernels may round differently).
        for stick in self.sticks.chunks_exact_mut(nx) {
            self.kernels.inverse[0].process_with_scratch(stick, &mut self.kernels.scratch);
            for z in stick {
                *z *= 1.0 / nx as f64;
            }
        }
        let nk = self.planes.len();
        for (k, columns) in self.plane_sticks.chunks_exact(ny).enumerate() {
            for x in 0..nx {
                for (y, &slot) in columns.iter().enumerate() {
                    self.kernels.lane[y] = if slot == usize::MAX {
                        C::default()
                    } else {
                        self.sticks[slot * nx + x]
                    };
                }
                self.kernels.lane(1, true, ny);
                for y in 0..ny {
                    self.partial[(x * ny + y) * nk + k] = self.kernels.lane[y];
                }
            }
        }
    }
    fn inverse_z(&mut self, xy: usize) {
        let nz = self.shape[2];
        self.kernels.lane[..nz].fill(C::default());
        for (k, &z) in self.planes.iter().enumerate() {
            self.kernels.lane[z] = self.partial[xy * self.planes.len() + k];
        }
        self.kernels.lane(2, true, nz);
        for z in &mut self.kernels.lane[..nz] {
            *z *= self.size as f64;
        }
    }
    fn validate(&self, c: &[C], grid_len: usize) -> Result<(), FftError> {
        if c.len() != self.indices.len() || grid_len != self.size {
            Err(FftError::IncompatibleArray)
        } else {
            Ok(())
        }
    }
}
impl FusedOperator for OrderedFusedPlan {
    fn shape(&self) -> [usize; 3] {
        self.shape
    }
    fn coefficient_count(&self) -> usize {
        self.indices.len()
    }
    fn apply_local(&mut self, c: &[C], potential: &[f64], out: &mut [C]) -> Result<(), FftError> {
        self.validate(c, potential.len())?;
        if out.len() != c.len() {
            return Err(FftError::IncompatibleArray);
        }
        let [nx, ny, nz] = self.shape;
        self.inverse_xy(c);
        // Independent y planes may be completed in any order without changing
        // any lane's floating-point operation sequence.
        for y in 0..ny {
            for x in 0..nx {
                self.inverse_z(x * ny + y);
                for z in 0..nz {
                    self.real_plane[x * nz + z] =
                        self.kernels.lane[z] * potential[(x * ny + y) * nz + z];
                }
            }
            for z in 0..nz {
                for x in 0..nx {
                    self.kernels.lane[x] = self.real_plane[x * nz + z];
                }
                self.kernels.lane(0, false, nx);
                for (slot, &x) in self.active_x.iter().enumerate() {
                    self.forward[(slot * ny + y) * nz + z] = self.kernels.lane[x];
                }
            }
        }
        for slot in 0..self.active_x.len() {
            for z in 0..nz {
                for y in 0..ny {
                    self.kernels.lane[y] = self.forward[(slot * ny + y) * nz + z];
                }
                self.kernels.lane(1, false, ny);
                for y in 0..ny {
                    self.forward[(slot * ny + y) * nz + z] = self.kernels.lane[y];
                }
            }
        }
        for &start in &self.forward_z {
            self.kernels.forward[2].process_with_scratch(
                &mut self.forward[start..start + nz],
                &mut self.kernels.scratch,
            );
        }
        for (dst, &slot) in out.iter_mut().zip(&self.gather) {
            *dst = self.forward[slot] / self.size as f64;
        }
        Ok(())
    }
    fn accumulate_density(
        &mut self,
        c: &[C],
        weight: f64,
        volume: f64,
        density: &mut [f64],
    ) -> Result<(), FftError> {
        self.validate(c, density.len())?;
        if !volume.is_finite() || volume <= 0.0 || !weight.is_finite() {
            return Err(FftError::IncompatibleArray);
        }
        if weight == 0.0 {
            return Ok(());
        }
        self.inverse_xy(c);
        let [nx, ny, nz] = self.shape;
        for xy in 0..nx * ny {
            self.inverse_z(xy);
            for z in 0..nz {
                density[xy * nz + z] += weight * self.kernels.lane[z].norm_sqr() / volume;
            }
        }
        Ok(())
    }
}

/// Experimental z-stick / xy-plane implementation. Inverse z/y/x differs from
/// the legacy x/y/z order; this backend does NOT promise bitwise equivalence.
/// Storage is compact z sticks and two xy planes, never a full transform grid.
/// Inactive x columns of `plane` stay canonical zero for the plan's lifetime.
/// Pure Rust has no native 2D plan: dense planes still use separable 1D calls.
/// Pruning y columns is enabled only at <=60% occupancy; otherwise all y lanes
/// are evaluated. Both policies keep the same x/y arithmetic within a plane.
pub struct PlaneFusedPlan {
    shape: [usize; 3],
    size: usize,
    kernels: Kernels,
    indices: Vec<[usize; 3]>,
    scatter: Vec<usize>,
    slots: Vec<usize>,
    columns: Vec<usize>,
    gather_columns: Vec<usize>,
    sticks: Vec<C>,
    plane: Vec<C>,
    real: Vec<C>,
}
impl PlaneFusedPlan {
    pub fn new(shape: [usize; 3], indices: &[[usize; 3]]) -> Result<Self, FftError> {
        Self::with_pruning(shape, indices, true)
    }
    /// Disable pruning to measure the dense separable plane control. `true`
    /// requests pruning only when the active-column share is at most 60%.
    pub fn with_pruning(
        shape: [usize; 3],
        indices: &[[usize; 3]],
        prune: bool,
    ) -> Result<Self, FftError> {
        let pattern = SpherePattern::from_indices(shape, indices)?;
        let [nx, ny, nz] = shape;
        let mut slots = vec![usize::MAX; nx * ny];
        for (slot, &start) in pattern.forward_lanes[2].iter().enumerate() {
            slots[start / nz] = slot;
        }
        let scatter = indices
            .iter()
            .map(|&[x, y, z]| slots[x * ny + y] * nz + z)
            .collect();
        let active: Vec<_> = (0..nx)
            .filter(|&x| slots[x * ny..(x + 1) * ny].iter().any(|&s| s != usize::MAX))
            .collect();
        let columns = if prune && active.len() * 5 <= nx * 3 {
            active.clone()
        } else {
            (0..nx).collect()
        };
        Ok(Self {
            shape,
            size: nx * ny * nz,
            kernels: Kernels::new(shape),
            indices: indices.to_vec(),
            scatter,
            slots,
            columns,
            gather_columns: active,
            sticks: vec![C::default(); pattern.forward_lanes[2].len() * nz],
            plane: vec![C::default(); nx * ny],
            real: vec![C::default(); nx * ny],
        })
    }
    fn validate(&self, c: &[C], n: usize) -> Result<(), FftError> {
        if c.len() != self.indices.len() || n != self.size {
            Err(FftError::IncompatibleArray)
        } else {
            Ok(())
        }
    }
    fn inverse_sticks(&mut self, c: &[C]) {
        self.sticks.fill(C::default());
        for (&slot, &value) in self.scatter.iter().zip(c) {
            self.sticks[slot] = value;
        }
        // Packed contiguous sticks permit one batched rustfft invocation.
        if !self.sticks.is_empty() {
            self.kernels.inverse[2]
                .process_with_scratch(&mut self.sticks, &mut self.kernels.scratch);
            for z in &mut self.sticks {
                *z *= 1.0 / self.shape[2] as f64;
            }
        }
    }
    fn inverse_plane(&mut self, z: usize) {
        let [nx, ny, nz] = self.shape;
        for &x in &self.columns {
            let lane = &mut self.plane[x * ny..(x + 1) * ny];
            for (y, v) in lane.iter_mut().enumerate() {
                let slot = self.slots[x * ny + y];
                *v = if slot == usize::MAX {
                    C::default()
                } else {
                    self.sticks[slot * nz + z]
                };
            }
            self.kernels.inverse[1].process_with_scratch(lane, &mut self.kernels.scratch);
            for v in lane {
                *v *= 1.0 / ny as f64;
            }
        }
        for y in 0..ny {
            for x in 0..nx {
                self.kernels.lane[x] = self.plane[x * ny + y];
            }
            self.kernels.lane(0, true, nx);
            for x in 0..nx {
                self.real[x * ny + y] = self.kernels.lane[x] * self.size as f64;
            }
        }
    }
    fn forward_plane(&mut self, z: usize) {
        let [nx, ny, nz] = self.shape;
        for y in 0..ny {
            for x in 0..nx {
                self.kernels.lane[x] = self.real[x * ny + y];
            }
            self.kernels.lane(0, false, nx);
            for &x in &self.gather_columns {
                self.plane[x * ny + y] = self.kernels.lane[x];
            }
        }
        for &x in &self.columns {
            let lane = &mut self.plane[x * ny..(x + 1) * ny];
            self.kernels.forward[1].process_with_scratch(lane, &mut self.kernels.scratch);
            for (y, &v) in lane.iter().enumerate() {
                let slot = self.slots[x * ny + y];
                if slot != usize::MAX {
                    self.sticks[slot * nz + z] = v;
                }
            }
        }
    }
}
impl FusedOperator for PlaneFusedPlan {
    fn shape(&self) -> [usize; 3] {
        self.shape
    }
    fn coefficient_count(&self) -> usize {
        self.indices.len()
    }
    fn apply_local(&mut self, c: &[C], potential: &[f64], out: &mut [C]) -> Result<(), FftError> {
        self.validate(c, potential.len())?;
        if out.len() != c.len() {
            return Err(FftError::IncompatibleArray);
        }
        self.inverse_sticks(c);
        let [nx, ny, nz] = self.shape;
        for z in 0..nz {
            self.inverse_plane(z);
            for xy in 0..nx * ny {
                self.real[xy] *= potential[xy * nz + z];
            }
            self.forward_plane(z);
        }
        if !self.sticks.is_empty() {
            self.kernels.forward[2]
                .process_with_scratch(&mut self.sticks, &mut self.kernels.scratch);
        }
        for (dst, &slot) in out.iter_mut().zip(&self.scatter) {
            *dst = self.sticks[slot] / self.size as f64;
        }
        Ok(())
    }
    fn accumulate_density(
        &mut self,
        c: &[C],
        weight: f64,
        volume: f64,
        density: &mut [f64],
    ) -> Result<(), FftError> {
        self.validate(c, density.len())?;
        if !volume.is_finite() || volume <= 0.0 || !weight.is_finite() {
            return Err(FftError::IncompatibleArray);
        }
        if weight == 0.0 {
            return Ok(());
        }
        self.inverse_sticks(c);
        let [nx, ny, nz] = self.shape;
        for z in 0..nz {
            self.inverse_plane(z);
            for xy in 0..nx * ny {
                density[xy * nz + z] += weight * self.real[xy].norm_sqr() / volume;
            }
        }
        Ok(())
    }
}

/// Selected backend for physics callers. MKL keeps its existing full-grid
/// backend until a vendor fused implementation is measured; it never silently
/// falls back to Rust FFTs when the `mkl` feature is selected.
#[cfg(not(feature = "mkl"))]
pub type FusedPlan = OrderedFusedPlan;
#[cfg(feature = "mkl")]
pub type FusedPlan = GridFusedPlan<crate::FftPlan3d>;

/// Allocation-free compatibility adapter for an existing `Fft3d` backend.
/// This preserves backend selection while vendor-specific fused kernels are
/// developed. Unlike the optimized Rust plans it retains full-grid scratch.
/// It is also a testable fallback for any conforming caller-supplied FFT.
pub struct GridFusedPlan<F: crate::Fft3d> {
    fft: F,
    indices: Vec<[usize; 3]>,
    reciprocal: ndarray::Array3<C>,
    real: ndarray::Array3<C>,
}
impl<F: crate::Fft3d> GridFusedPlan<F> {
    pub fn from_fft(fft: F, indices: &[[usize; 3]]) -> Result<Self, FftError> {
        let shape = fft.shape();
        SpherePattern::from_indices(shape, indices)?;
        Ok(Self {
            fft,
            indices: indices.to_vec(),
            reciprocal: ndarray::Array3::zeros(shape),
            real: ndarray::Array3::zeros(shape),
        })
    }
    fn validate(&self, c: &[C], n: usize) -> Result<(), FftError> {
        if c.len() != self.indices.len() || n != self.real.len() {
            Err(FftError::IncompatibleArray)
        } else {
            Ok(())
        }
    }
    fn inverse(&mut self, c: &[C]) -> Result<(), FftError> {
        self.real.fill(C::default());
        for (&index, &value) in self.indices.iter().zip(c) {
            self.real[index] = value;
        }
        self.fft
            .transform_in_place(&mut self.real, &mut self.reciprocal, true)
    }
}
impl GridFusedPlan<crate::FftPlan3d> {
    pub fn new(shape: [usize; 3], indices: &[[usize; 3]]) -> Result<Self, FftError> {
        Self::from_fft(crate::FftPlan3d::new(shape)?, indices)
    }
}
impl<F: crate::Fft3d> FusedOperator for GridFusedPlan<F> {
    fn shape(&self) -> [usize; 3] {
        self.fft.shape()
    }
    fn coefficient_count(&self) -> usize {
        self.indices.len()
    }
    fn apply_local(&mut self, c: &[C], potential: &[f64], out: &mut [C]) -> Result<(), FftError> {
        self.validate(c, potential.len())?;
        if out.len() != c.len() {
            return Err(FftError::IncompatibleArray);
        }
        self.inverse(c)?;
        for (z, &v) in self.real.iter_mut().zip(potential) {
            *z *= v;
        }
        self.fft
            .transform_in_place(&mut self.real, &mut self.reciprocal, false)?;
        for (dst, &index) in out.iter_mut().zip(&self.indices) {
            *dst = self.real[index];
        }
        Ok(())
    }
    fn accumulate_density(
        &mut self,
        c: &[C],
        weight: f64,
        volume: f64,
        density: &mut [f64],
    ) -> Result<(), FftError> {
        self.validate(c, density.len())?;
        if !weight.is_finite() || !volume.is_finite() || volume <= 0.0 {
            return Err(FftError::IncompatibleArray);
        }
        if weight == 0.0 {
            return Ok(());
        }
        self.inverse(c)?;
        for (n, z) in density.iter_mut().zip(&self.real) {
            *n += weight * z.norm_sqr() / volume;
        }
        Ok(())
    }
}

// Additive real-Gamma density entries. The existing single-band kernels above
// remain untouched; reuse exactly their inverse stages and storage.
impl crate::fused_pair::PairDensityOperator for OrderedFusedPlan {
    fn accumulate_density_pair(
        &mut self,
        phi: &[C],
        weights: [f64; 2],
        volume: f64,
        density: &mut [f64],
    ) -> Result<(), FftError> {
        self.validate(phi, density.len())?;
        crate::fused_pair::validate_density(weights, volume)?;
        if weights == [0.; 2] {
            return Ok(());
        }
        self.inverse_xy(phi);
        let [nx, ny, nz] = self.shape;
        for xy in 0..nx * ny {
            self.inverse_z(xy);
            for z in 0..nz {
                let value = self.kernels.lane[z];
                density[xy * nz + z] +=
                    (weights[0] * value.re * value.re + weights[1] * value.im * value.im) / volume;
            }
        }
        Ok(())
    }
}
impl crate::fused_pair::SplitPairDensity for OrderedFusedPlan {
    fn accumulate_density_pair_split(
        &mut self,
        phi: &[C],
        weights: [f64; 2],
        volume: f64,
        a: &mut [f64],
        b: &mut [f64],
    ) -> Result<(), FftError> {
        self.validate(phi, a.len())?;
        crate::fused_pair::validate_density(weights, volume)?;
        if b.len() != a.len() {
            return Err(FftError::IncompatibleArray);
        }
        if weights == [0.; 2] {
            return Ok(());
        }
        self.inverse_xy(phi);
        let [nx, ny, nz] = self.shape;
        for xy in 0..nx * ny {
            self.inverse_z(xy);
            let range = xy * nz..(xy + 1) * nz;
            for ((value, a), b) in self.kernels.lane[..nz]
                .iter()
                .zip(&mut a[range.clone()])
                .zip(&mut b[range])
            {
                *a += weights[0] * (value.re * value.re) / volume;
                *b += weights[1] * (value.im * value.im) / volume;
            }
        }
        Ok(())
    }
}
impl<F: crate::Fft3d> crate::fused_pair::SplitPairDensity for GridFusedPlan<F> {
    fn accumulate_density_pair_split(
        &mut self,
        phi: &[C],
        weights: [f64; 2],
        volume: f64,
        a: &mut [f64],
        b: &mut [f64],
    ) -> Result<(), FftError> {
        self.validate(phi, a.len())?;
        crate::fused_pair::validate_density(weights, volume)?;
        if b.len() != a.len() {
            return Err(FftError::IncompatibleArray);
        }
        if weights == [0.; 2] {
            return Ok(());
        }
        self.inverse(phi)?;
        for ((value, a), b) in self.real.iter().zip(a.iter_mut()).zip(b.iter_mut()) {
            *a += weights[0] * (value.re * value.re) / volume;
            *b += weights[1] * (value.im * value.im) / volume;
        }
        Ok(())
    }
}
impl<F: crate::Fft3d> crate::fused_pair::PairDensityOperator for GridFusedPlan<F> {
    fn accumulate_density_pair(
        &mut self,
        phi: &[C],
        weights: [f64; 2],
        volume: f64,
        density: &mut [f64],
    ) -> Result<(), FftError> {
        self.validate(phi, density.len())?;
        crate::fused_pair::validate_density(weights, volume)?;
        if weights == [0.; 2] {
            return Ok(());
        }
        self.inverse(phi)?;
        for (dst, value) in density.iter_mut().zip(&self.real) {
            *dst += (weights[0] * value.re * value.re + weights[1] * value.im * value.im) / volume;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn vec_state<T>(v: &Vec<T>) -> (usize, usize, usize) {
        (v.as_ptr() as usize, v.len(), v.capacity())
    }
    #[test]
    fn compatibility_adapter_reuses_scratch_and_preserves_density() {
        let indices = [[0, 0, 0], [1, 1, 1], [2, 2, 2]];
        let mut plan =
            GridFusedPlan::from_fft(crate::RustFftPlan3d::new([5, 7, 9]).unwrap(), &indices)
                .unwrap();
        let pointers = (plan.real.as_ptr(), plan.reciprocal.as_ptr());
        let mut ordered = OrderedFusedPlan::new([5, 7, 9], &indices).unwrap();
        let c = [C::new(0.3, -0.7); 3];
        let mut actual = vec![0.; 315];
        let mut expected = actual.clone();
        let mut out = c;
        for weight in [1., 0., 0.3] {
            plan.apply_local(&c, &[0.7; 315], &mut out).unwrap();
            plan.accumulate_density(&c, weight, 20., &mut actual)
                .unwrap();
            ordered
                .accumulate_density(&c, weight, 20., &mut expected)
                .unwrap();
            assert_eq!(pointers, (plan.real.as_ptr(), plan.reciprocal.as_ptr()));
            assert_eq!(
                actual.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                expected.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
            );
        }
        let saved = actual.clone();
        assert!(plan.accumulate_density(&c, 1., 0., &mut actual).is_err());
        assert!(plan.apply_local(&[], &[0.7; 315], &mut out).is_err());
        assert_eq!(actual, saved);
        assert_eq!(pointers, (plan.real.as_ptr(), plan.reciprocal.as_ptr()));
    }
    #[test]
    fn dense_plane_policy_keeps_unrequested_columns_zero() {
        let indices = [[0, 0, 0], [1, 1, 1], [2, 2, 2], [3, 3, 3]];
        for prune in [false, true] {
            let mut plan = PlaneFusedPlan::with_pruning([5, 7, 9], &indices, prune).unwrap();
            let c = [C::new(0.3, -0.7); 4];
            let mut out = c;
            for _ in 0..3 {
                plan.apply_local(&c, &[0.7; 315], &mut out).unwrap();
                for z in &plan.plane[4 * 7..] {
                    assert_eq!((z.re.to_bits(), z.im.to_bits()), (0, 0));
                }
            }
        }
    }
    #[test]
    fn scratch_stability_and_inactive_plane_columns() {
        let indices = [[0, 0, 0], [1, 1, 1], [0, 1, 2]];
        let mut ordered = OrderedFusedPlan::new([5, 7, 9], &indices).unwrap();
        let mut plane = PlaneFusedPlan::new([5, 7, 9], &indices).unwrap();
        let ordered_state = |p: &OrderedFusedPlan| {
            [
                vec_state(&p.sticks),
                vec_state(&p.partial),
                vec_state(&p.forward),
                vec_state(&p.real_plane),
                vec_state(&p.kernels.scratch),
                vec_state(&p.kernels.lane),
            ]
        };
        let plane_state = |p: &PlaneFusedPlan| {
            [
                vec_state(&p.sticks),
                vec_state(&p.plane),
                vec_state(&p.real),
                vec_state(&p.kernels.scratch),
                vec_state(&p.kernels.lane),
            ]
        };
        let a = ordered_state(&ordered);
        let b = plane_state(&plane);
        let c = [C::new(0.3, -0.7); 3];
        let v = vec![0.7; 315];
        let mut out = [C::default(); 3];
        let mut rho = vec![0.; 315];
        for _ in 0..4 {
            ordered.apply_local(&c, &v, &mut out).unwrap();
            ordered.accumulate_density(&c, 1., 20., &mut rho).unwrap();
            plane.apply_local(&c, &v, &mut out).unwrap();
            plane.accumulate_density(&c, 1., 20., &mut rho).unwrap();
            assert!(ordered.apply_local(&[], &v, &mut out).is_err());
            assert!(plane.apply_local(&[], &v, &mut out).is_err());
            assert_eq!(a, ordered_state(&ordered));
            assert_eq!(b, plane_state(&plane));
            for z in &plane.plane[2 * 7..] {
                assert_eq!((z.re.to_bits(), z.im.to_bits()), (0, 0));
            }
        }
    }
}

#[cfg(test)]
mod batch_study {
    use super::*;
    use std::{hint::black_box, time::Instant};
    // Test-only alternative: all bands advance through each stage together,
    // and the multiply loads one potential value per plane point per group.
    // No new production batching policy is introduced by this experiment.
    fn grouped(
        plans: &mut [PlaneFusedPlan],
        coefficients: &[Vec<C>],
        potential: &[f64],
        out: &mut [Vec<C>],
    ) {
        let [nx, ny, nz] = plans[0].shape;
        for (p, c) in plans.iter_mut().zip(coefficients) {
            p.inverse_sticks(c);
        }
        for z in 0..nz {
            for p in plans.iter_mut() {
                p.inverse_plane(z);
            }
            for xy in 0..nx * ny {
                let v = potential[xy * nz + z];
                for p in plans.iter_mut() {
                    p.real[xy] *= v;
                }
            }
            for p in plans.iter_mut() {
                p.forward_plane(z);
            }
        }
        for (p, out) in plans.iter_mut().zip(out) {
            if !p.sticks.is_empty() {
                p.kernels.forward[2].process_with_scratch(&mut p.sticks, &mut p.kernels.scratch);
            }
            for (dst, &slot) in out.iter_mut().zip(&p.scatter) {
                *dst = p.sticks[slot] / p.size as f64;
            }
        }
    }
    #[test]
    #[ignore = "cargo test --release -p dft-numerics grouped_band_study -- --ignored --nocapture"]
    fn grouped_band_study() {
        for shape in [[25, 25, 25], [36, 36, 36], [15, 20, 24]] {
            let [nx, ny, nz] = shape;
            let mut indices = Vec::new();
            let signed = |i: usize, n: usize| {
                if i <= n / 2 {
                    i as f64
                } else {
                    i as f64 - n as f64
                }
            };
            for x in 0..nx {
                for y in 0..ny {
                    for z in 0..nz {
                        if signed(x, nx).powi(2) + signed(y, ny).powi(2) + signed(z, nz).powi(2)
                            < (nx as f64 / 4.).powi(2)
                        {
                            indices.push([x, y, z]);
                        }
                    }
                }
            }
            let potential: Vec<_> = (0..nx * ny * nz).map(|i| (i as f64).sin()).collect();
            let coefficients: Vec<Vec<C>> = (0..32)
                .map(|b| {
                    (0..indices.len())
                        .map(|i| C::new(((i + b * 3) as f64).sin(), ((i + b * 7) as f64).cos()))
                        .collect()
                })
                .collect();
            let mut oracle = PlaneFusedPlan::new(shape, &indices).unwrap();
            let mut expected = vec![vec![C::default(); indices.len()]; 32];
            for (c, o) in coefficients.iter().zip(&mut expected) {
                oracle.apply_local(c, &potential, o).unwrap();
            }
            for batch in [1, 2, 4, 8, 16, 32] {
                let mut plans: Vec<_> = (0..batch)
                    .map(|_| PlaneFusedPlan::new(shape, &indices).unwrap())
                    .collect();
                let mut out = expected.clone();
                let mut samples = Vec::new();
                for _ in 0..5 {
                    let t = Instant::now();
                    for _ in 0..4 {
                        for (c, o) in coefficients.chunks(batch).zip(out.chunks_mut(batch)) {
                            grouped(&mut plans, black_box(c), &potential, o);
                        }
                        black_box(&out);
                    }
                    samples.push(t.elapsed().as_secs_f64() * 1e6 / (4. * 32.));
                    for (a, b) in out.iter().flatten().zip(expected.iter().flatten()) {
                        assert_eq!(
                            (a.re.to_bits(), a.im.to_bits()),
                            (b.re.to_bits(), b.im.to_bits())
                        );
                    }
                }
                println!("grouped shape={shape:?} batch={batch} us_per_band={samples:?}");
            }
        }
    }
}
