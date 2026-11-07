//! Mirror-slab FFT for real Gamma pairs (P39A, pair-fft-efficiency section 5).
//! Unscaled inverse, mirrored forward, and one final 1/N factor. Rounding
//! differs from OrderedFusedPlan; this type never replaces the complex alias.
use crate::{FftError, fused::FusedOperator, fused_pair::PairDensityOperator};
use num_complex::Complex64 as C;
use rustfft::{Fft, FftPlanner};
use std::{collections::BTreeSet, sync::Arc};

/// Backend selection lives in numerics, which owns the unified `mkl` feature.
/// Physics has no independent backend feature. MKL retains its existing plan.
#[cfg(not(feature = "mkl"))]
pub type RealPairPlan = PairFusedPlan;
#[cfg(feature = "mkl")]
pub type RealPairPlan = crate::fused::FusedPlan;

/// Serial, construction-sized pair workspace. The last axis is always z;
/// choose the fewest occupied sticks among x/y (ties choose x). Keeping z
/// last takes precedence over minimizing sticks over all three axes.
pub struct PairFusedPlan {
    shape: [usize; 3],
    first: usize,
    n: [usize; 3],
    size: usize,
    scatter: Vec<usize>,
    planes: Vec<usize>,
    // (index in the active-plane middle batch, stick slot), plane-major.
    occupied: Vec<(usize, usize)>,
    sticks: Vec<C>,
    middle: Vec<C>,
    slab: Vec<C>,
    scratch: Vec<C>,
    inverse: [Arc<dyn Fft<f64>>; 3],
    forward: [Arc<dyn Fft<f64>>; 3],
}
impl PairFusedPlan {
    pub fn new(shape: [usize; 3], indices: &[[usize; 3]]) -> Result<Self, FftError> {
        let size = shape
            .iter()
            .try_fold(1usize, |a, &b| a.checked_mul(b))
            .ok_or(FftError::InvalidGrid)?;
        if size == 0 || size > isize::MAX as usize / size_of::<C>() {
            return Err(FftError::InvalidGrid);
        }
        if indices.iter().any(|p| (0..3).any(|d| p[d] >= shape[d])) {
            return Err(FftError::IncompatibleArray);
        }
        let supports: [_; 2] = std::array::from_fn(|a| {
            indices
                .iter()
                .map(|p| (p[2], p[1 - a]))
                .collect::<BTreeSet<_>>()
        });
        let first = usize::from(supports[1].len() < supports[0].len());
        let n = [shape[first], shape[1 - first], shape[2]];
        let planes: Vec<_> = supports[first]
            .iter()
            .map(|p| p.0)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        let mut slots = vec![usize::MAX; n[1] * n[2]];
        let mut occupied = Vec::new();
        for (slot, &(z, y)) in supports[first].iter().enumerate() {
            slots[z * n[1] + y] = slot;
            occupied.push((planes.binary_search(&z).unwrap() * n[1] + y, slot));
        }
        let scatter = indices
            .iter()
            .map(|p| slots[p[2] * n[1] + p[1 - first]] * n[0] + p[first])
            .collect();
        let mut planner = FftPlanner::new();
        let inverse = n.map(|d| planner.plan_fft_inverse(d));
        let forward = n.map(|d| planner.plan_fft_forward(d));
        let scratch = inverse
            .iter()
            .chain(&forward)
            .map(|p| p.get_inplace_scratch_len())
            .max()
            .unwrap();
        Ok(Self {
            shape,
            first,
            n,
            size,
            scatter,
            middle: vec![C::default(); planes.len() * n[1]],
            sticks: vec![C::default(); supports[first].len() * n[0]],
            slab: vec![C::default(); n[1] * n[2]],
            scratch: vec![C::default(); scratch],
            planes,
            occupied,
            inverse,
            forward,
        })
    }
    fn validate(&self, c: &[C], grid: usize) -> Result<(), FftError> {
        if c.len() != self.scatter.len() || grid != self.size {
            Err(FftError::IncompatibleArray)
        } else {
            Ok(())
        }
    }
    fn inverse_sticks(&mut self, c: &[C]) {
        self.sticks.fill(C::default());
        for (&i, &v) in self.scatter.iter().zip(c) {
            self.sticks[i] = v;
        }
        if !self.sticks.is_empty() {
            self.inverse[0].process_with_scratch(&mut self.sticks, &mut self.scratch);
        }
    }
    fn inverse_slab(&mut self, x: usize) {
        let [nx, ny, nz] = self.n;
        self.middle.fill(C::default());
        for &(i, slot) in &self.occupied {
            self.middle[i] = self.sticks[slot * nx + x];
        }
        if !self.middle.is_empty() {
            self.inverse[1].process_with_scratch(&mut self.middle, &mut self.scratch);
        }
        self.slab.fill(C::default());
        for (k, &z) in self.planes.iter().enumerate() {
            for y in 0..ny {
                self.slab[y * nz + z] = self.middle[k * ny + y];
            }
        }
        self.inverse[2].process_with_scratch(&mut self.slab, &mut self.scratch);
    }
    fn grid_start(&self, x: usize, y: usize) -> usize {
        if self.first == 0 {
            (x * self.shape[1] + y) * self.n[2]
        } else {
            (y * self.shape[1] + x) * self.n[2]
        }
    }
    fn forward_slab(&mut self, x: usize) {
        let [nx, ny, nz] = self.n;
        self.forward[2].process_with_scratch(&mut self.slab, &mut self.scratch);
        for (k, &z) in self.planes.iter().enumerate() {
            for y in 0..ny {
                self.middle[k * ny + y] = self.slab[y * nz + z];
            }
        }
        if !self.middle.is_empty() {
            self.forward[1].process_with_scratch(&mut self.middle, &mut self.scratch);
        }
        for &(i, slot) in &self.occupied {
            self.sticks[slot * nx + x] = self.middle[i];
        }
    }
}
impl FusedOperator for PairFusedPlan {
    fn shape(&self) -> [usize; 3] {
        self.shape
    }
    fn coefficient_count(&self) -> usize {
        self.scatter.len()
    }
    fn apply_local(&mut self, c: &[C], potential: &[f64], out: &mut [C]) -> Result<(), FftError> {
        self.validate(c, potential.len())?;
        if out.len() != c.len() {
            return Err(FftError::IncompatibleArray);
        }
        self.inverse_sticks(c);
        for x in 0..self.n[0] {
            self.inverse_slab(x);
            for y in 0..self.n[1] {
                let start = self.grid_start(x, y);
                for (a, &v) in self.slab[y * self.n[2]..(y + 1) * self.n[2]]
                    .iter_mut()
                    .zip(&potential[start..start + self.n[2]])
                {
                    *a *= v;
                }
            }
            self.forward_slab(x);
        }
        if !self.sticks.is_empty() {
            self.forward[0].process_with_scratch(&mut self.sticks, &mut self.scratch);
        }
        for (o, &i) in out.iter_mut().zip(&self.scatter) {
            *o = self.sticks[i] / self.size as f64;
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
        crate::fused_pair::validate_density([weight; 2], volume)?;
        if weight == 0. {
            return Ok(());
        }
        self.inverse_sticks(c);
        for x in 0..self.n[0] {
            self.inverse_slab(x);
            for y in 0..self.n[1] {
                let start = self.grid_start(x, y);
                for z in 0..self.n[2] {
                    density[start + z] += weight * self.slab[y * self.n[2] + z].norm_sqr() / volume;
                }
            }
        }
        Ok(())
    }
}
impl PairDensityOperator for PairFusedPlan {
    fn accumulate_density_pair(
        &mut self,
        c: &[C],
        weights: [f64; 2],
        volume: f64,
        density: &mut [f64],
    ) -> Result<(), FftError> {
        self.validate(c, density.len())?;
        crate::fused_pair::validate_density(weights, volume)?;
        if weights == [0.; 2] {
            return Ok(());
        }
        self.inverse_sticks(c);
        for x in 0..self.n[0] {
            self.inverse_slab(x);
            for y in 0..self.n[1] {
                let start = self.grid_start(x, y);
                for z in 0..self.n[2] {
                    let v = self.slab[y * self.n[2] + z];
                    density[start + z] +=
                        (weights[0] * v.re * v.re + weights[1] * v.im * v.im) / volume;
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn retained_storage() {
        let mut p =
            PairFusedPlan::new([15, 20, 24], &[[0, 0, 0], [1, 2, 3], [14, 18, 21]]).unwrap();
        fn state(p: &PairFusedPlan) -> [(usize, usize, usize); 4] {
            [&p.sticks, &p.middle, &p.slab, &p.scratch]
                .map(|v| (v.as_ptr() as usize, v.len(), v.capacity()))
        }
        let before = state(&p);
        let mut out = [C::default(); 3];
        let mut rho = vec![0.; 7200];
        for _ in 0..4 {
            p.apply_local(&[C::new(0.3, 0.7); 3], &[0.2; 7200], &mut out)
                .unwrap();
            p.accumulate_density_pair(&out, [0.3, 0.7], 10., &mut rho)
                .unwrap();
            assert_eq!(before, state(&p));
        }
    }
}
