//! Mirror-slab FFT for real Gamma pairs (P39A, pair-fft-efficiency section 5).
//! Unscaled inverse, mirrored forward, and one final 1/N factor. Rounding
//! differs from OrderedFusedPlan; this type never replaces the complex alias.
use crate::{FftError, fused::FusedOperator, fused_pair::PairDensityOperator};
use num_complex::Complex64 as C;
use rustfft::{Fft, FftPlanner};
use std::{collections::BTreeSet, sync::Arc};

/// Backend selection lives in numerics, which owns the unified `mkl` feature.
/// Physics has no independent backend feature. MKL retains its existing plan
/// unless `BAANDI_PAIR_FFT=pruned` selects the mirror-slab layout on MKL lanes.
#[cfg(not(feature = "mkl"))]
pub type RealPairPlan = PairFusedPlan;
#[cfg(feature = "mkl")]
pub type RealPairPlan = MklPairPlan;

/// Unscaled in-place batches of 1D transforms for the three lane buffers of
/// `PairFusedPlan`: axis `a` always acts on whole buffers of `counts[a]`
/// contiguous lanes of length `n[a]` (sticks, middle, slab).
pub trait PairLanes: Sized {
    fn plan(n: [usize; 3], counts: [usize; 3]) -> Result<Self, FftError>;
    fn run(&mut self, axis: usize, inverse: bool, data: &mut [C]) -> Result<(), FftError>;
}

/// rustfft lanes (the historical kernels; bit-identical to the pre-trait code).
pub struct RustLanes {
    scratch: Vec<C>,
    inverse: [Arc<dyn Fft<f64>>; 3],
    forward: [Arc<dyn Fft<f64>>; 3],
}
impl PairLanes for RustLanes {
    fn plan(n: [usize; 3], _: [usize; 3]) -> Result<Self, FftError> {
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
            scratch: vec![C::default(); scratch],
            inverse,
            forward,
        })
    }
    fn run(&mut self, axis: usize, inverse: bool, data: &mut [C]) -> Result<(), FftError> {
        let plan = if inverse {
            &self.inverse[axis]
        } else {
            &self.forward[axis]
        };
        plan.process_with_scratch(data, &mut self.scratch);
        Ok(())
    }
}

/// MKL DFTI lanes: one batch descriptor per axis buffer (rounding differs from
/// rustfft; results stay independent of the thread budget, one plan per piece).
#[cfg(feature = "mkl")]
pub struct MklLanes {
    plans: [Option<dft_mkl::Plan1dBatch>; 3],
}
#[cfg(feature = "mkl")]
impl PairLanes for MklLanes {
    fn plan(n: [usize; 3], counts: [usize; 3]) -> Result<Self, FftError> {
        let mut plans = [None, None, None];
        for axis in 0..3 {
            if counts[axis] > 0 {
                plans[axis] = Some(dft_mkl::Plan1dBatch::new(n[axis], counts[axis])?);
            }
        }
        Ok(Self { plans })
    }
    fn run(&mut self, axis: usize, inverse: bool, data: &mut [C]) -> Result<(), FftError> {
        match &mut self.plans[axis] {
            Some(plan) => Ok(plan.transform_in_place(data, inverse)?),
            None if data.is_empty() => Ok(()),
            None => Err(FftError::IncompatibleArray),
        }
    }
}

/// MKL builds: the existing full-grid plan (default, unchanged bits) or, with
/// `BAANDI_PAIR_FFT=pruned` read at construction, the mirror-slab pair plan on
/// MKL 1D lanes (rounding class).
#[cfg(feature = "mkl")]
pub enum MklPairPlan {
    Grid(crate::fused::FusedPlan),
    Pruned(PairFusedPlan<MklLanes>),
}
#[cfg(feature = "mkl")]
impl MklPairPlan {
    pub fn new(shape: [usize; 3], indices: &[[usize; 3]]) -> Result<Self, FftError> {
        if crate::switches::choice(crate::switches::Id::PAIR_FFT) == "pruned" {
            Ok(Self::Pruned(PairFusedPlan::with_lanes(shape, indices)?))
        } else {
            Ok(Self::Grid(crate::fused::FusedPlan::new(shape, indices)?))
        }
    }
}
#[cfg(feature = "mkl")]
impl FusedOperator for MklPairPlan {
    fn shape(&self) -> [usize; 3] {
        match self {
            Self::Grid(p) => p.shape(),
            Self::Pruned(p) => p.shape(),
        }
    }
    fn coefficient_count(&self) -> usize {
        match self {
            Self::Grid(p) => p.coefficient_count(),
            Self::Pruned(p) => p.coefficient_count(),
        }
    }
    fn apply_local(&mut self, c: &[C], potential: &[f64], out: &mut [C]) -> Result<(), FftError> {
        match self {
            Self::Grid(p) => p.apply_local(c, potential, out),
            Self::Pruned(p) => p.apply_local(c, potential, out),
        }
    }
    fn apply_local_laps(
        &mut self,
        c: &[C],
        potential: &[f64],
        out: &mut [C],
        laps: &mut crate::work_timers::LocalLaps,
    ) -> Result<(), FftError> {
        match self {
            Self::Grid(p) => p.apply_local_laps(c, potential, out, laps),
            Self::Pruned(p) => p.apply_local_laps(c, potential, out, laps),
        }
    }
    fn accumulate_density(
        &mut self,
        c: &[C],
        weight: f64,
        volume: f64,
        density: &mut [f64],
    ) -> Result<(), FftError> {
        match self {
            Self::Grid(p) => p.accumulate_density(c, weight, volume, density),
            Self::Pruned(p) => p.accumulate_density(c, weight, volume, density),
        }
    }
}
#[cfg(feature = "mkl")]
impl PairDensityOperator for MklPairPlan {
    fn accumulate_density_pair(
        &mut self,
        c: &[C],
        weights: [f64; 2],
        volume: f64,
        density: &mut [f64],
    ) -> Result<(), FftError> {
        match self {
            Self::Grid(p) => p.accumulate_density_pair(c, weights, volume, density),
            Self::Pruned(p) => p.accumulate_density_pair(c, weights, volume, density),
        }
    }
}

/// Serial, construction-sized pair workspace. The last axis is always z;
/// choose the fewest occupied sticks among x/y (ties choose x). Keeping z
/// last takes precedence over minimizing sticks over all three axes.
pub struct PairFusedPlan<L: PairLanes = RustLanes> {
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
    lanes: L,
}
impl PairFusedPlan<RustLanes> {
    pub fn new(shape: [usize; 3], indices: &[[usize; 3]]) -> Result<Self, FftError> {
        Self::with_lanes(shape, indices)
    }
}
impl<L: PairLanes> PairFusedPlan<L> {
    /// Same layout with the selected 1D lane kernels.
    pub fn with_lanes(shape: [usize; 3], indices: &[[usize; 3]]) -> Result<Self, FftError> {
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
        let lanes = L::plan(n, [supports[first].len(), planes.len(), n[1]])?;
        Ok(Self {
            shape,
            first,
            n,
            size,
            scatter,
            middle: vec![C::default(); planes.len() * n[1]],
            sticks: vec![C::default(); supports[first].len() * n[0]],
            slab: vec![C::default(); n[1] * n[2]],
            planes,
            occupied,
            lanes,
        })
    }
    fn validate(&self, c: &[C], grid: usize) -> Result<(), FftError> {
        if c.len() != self.scatter.len() || grid != self.size {
            Err(FftError::IncompatibleArray)
        } else {
            Ok(())
        }
    }
    fn inverse_sticks(&mut self, c: &[C]) -> Result<(), FftError> {
        self.sticks.fill(C::default());
        for (&i, &v) in self.scatter.iter().zip(c) {
            self.sticks[i] = v;
        }
        if !self.sticks.is_empty() {
            self.lanes.run(0, true, &mut self.sticks)?;
        }
        Ok(())
    }
    fn inverse_slab(&mut self, x: usize) -> Result<(), FftError> {
        let [nx, ny, nz] = self.n;
        self.middle.fill(C::default());
        for &(i, slot) in &self.occupied {
            self.middle[i] = self.sticks[slot * nx + x];
        }
        if !self.middle.is_empty() {
            self.lanes.run(1, true, &mut self.middle)?;
        }
        self.slab.fill(C::default());
        for (k, &z) in self.planes.iter().enumerate() {
            for y in 0..ny {
                self.slab[y * nz + z] = self.middle[k * ny + y];
            }
        }
        self.lanes.run(2, true, &mut self.slab)
    }
    fn grid_start(&self, x: usize, y: usize) -> usize {
        if self.first == 0 {
            (x * self.shape[1] + y) * self.n[2]
        } else {
            (y * self.shape[1] + x) * self.n[2]
        }
    }
    fn forward_slab(&mut self, x: usize) -> Result<(), FftError> {
        let [nx, ny, nz] = self.n;
        self.lanes.run(2, false, &mut self.slab)?;
        for (k, &z) in self.planes.iter().enumerate() {
            for y in 0..ny {
                self.middle[k * ny + y] = self.slab[y * nz + z];
            }
        }
        if !self.middle.is_empty() {
            self.lanes.run(1, false, &mut self.middle)?;
        }
        for &(i, slot) in &self.occupied {
            self.sticks[slot * nx + x] = self.middle[i];
        }
        Ok(())
    }
    /// The single body of `apply_local`; `NoLaps` makes the laps vanish. Laps
    /// are taken per slab (inverse, multiply, forward), not per element.
    fn apply_local_with<K: crate::work_timers::Laps>(
        &mut self,
        c: &[C],
        potential: &[f64],
        out: &mut [C],
        laps: &mut K,
    ) -> Result<(), FftError> {
        use crate::work_timers::Phase;
        self.validate(c, potential.len())?;
        if out.len() != c.len() {
            return Err(FftError::IncompatibleArray);
        }
        self.inverse_sticks(c)?;
        laps.lap(Phase::Inverse);
        for x in 0..self.n[0] {
            self.inverse_slab(x)?;
            laps.lap(Phase::Inverse);
            for y in 0..self.n[1] {
                let start = self.grid_start(x, y);
                for (a, &v) in self.slab[y * self.n[2]..(y + 1) * self.n[2]]
                    .iter_mut()
                    .zip(&potential[start..start + self.n[2]])
                {
                    *a *= v;
                }
            }
            laps.lap(Phase::Potential);
            self.forward_slab(x)?;
            laps.lap(Phase::Forward);
        }
        if !self.sticks.is_empty() {
            self.lanes.run(0, false, &mut self.sticks)?;
        }
        for (o, &i) in out.iter_mut().zip(&self.scatter) {
            *o = self.sticks[i] / self.size as f64;
        }
        laps.lap(Phase::Forward);
        Ok(())
    }
}
impl<L: PairLanes> FusedOperator for PairFusedPlan<L> {
    fn shape(&self) -> [usize; 3] {
        self.shape
    }
    fn coefficient_count(&self) -> usize {
        self.scatter.len()
    }
    fn apply_local(&mut self, c: &[C], potential: &[f64], out: &mut [C]) -> Result<(), FftError> {
        self.apply_local_with(c, potential, out, &mut crate::work_timers::NoLaps)
    }
    fn apply_local_laps(
        &mut self,
        c: &[C],
        potential: &[f64],
        out: &mut [C],
        laps: &mut crate::work_timers::LocalLaps,
    ) -> Result<(), FftError> {
        self.apply_local_with(c, potential, out, laps)
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
        self.inverse_sticks(c)?;
        for x in 0..self.n[0] {
            self.inverse_slab(x)?;
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
impl<L: PairLanes> PairDensityOperator for PairFusedPlan<L> {
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
        self.inverse_sticks(c)?;
        for x in 0..self.n[0] {
            self.inverse_slab(x)?;
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
            [&p.sticks, &p.middle, &p.slab, &p.lanes.scratch]
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

    /// dev08 (MKL builds only): the mirror-slab plan on MKL DFTI lanes equals
    /// the rustfft lanes to relative 1e-14 (local and pair density), is
    /// repeatable bit for bit, and the default `MklPairPlan` stays the grid plan.
    /// cargo test --features mkl -p dft-numerics --lib mkl_lanes_match_rust_lanes
    #[cfg(feature = "mkl")]
    #[test]
    fn mkl_lanes_match_rust_lanes() {
        for shape in [[15, 20, 24], [6, 6, 6], [25, 9, 36], [36, 36, 36]] {
            let h = |n: usize| (n as i64 - 1) / 4;
            let wrap = |v: i64, n: usize| v.rem_euclid(n as i64) as usize;
            let mut indices = Vec::new();
            for x in -h(shape[0])..=h(shape[0]) {
                for y in -h(shape[1])..=h(shape[1]) {
                    for z in -h(shape[2])..=h(shape[2]) {
                        if x * x + y * y + z * z <= h(shape[0]).pow(2) + 1 {
                            indices.push([wrap(x, shape[0]), wrap(y, shape[1]), wrap(z, shape[2])]);
                        }
                    }
                }
            }
            let size: usize = shape.iter().product();
            let v: Vec<f64> = (0..size)
                .map(|i| 0.3 * ((i * 7 + 1) as f64).sin())
                .collect();
            let c: Vec<C> = (0..indices.len())
                .map(|i| C::new(((i * 3 + 1) as f64).cos(), ((i * 5 + 2) as f64).sin()))
                .collect();
            let mut rust = PairFusedPlan::new(shape, &indices).unwrap();
            let mut mkl = PairFusedPlan::<MklLanes>::with_lanes(shape, &indices).unwrap();
            let (mut a, mut b, mut b2) = (vec![C::default(); c.len()], c.clone(), c.clone());
            rust.apply_local(&c, &v, &mut a).unwrap();
            mkl.apply_local(&c, &v, &mut b).unwrap();
            mkl.apply_local(&c, &v, &mut b2).unwrap();
            assert_eq!(b, b2);
            let norm = a.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
            let diff = a
                .iter()
                .zip(&b)
                .map(|(x, y)| (x - y).norm_sqr())
                .sum::<f64>()
                .sqrt();
            assert!(diff <= 1e-14 * norm, "{shape:?} local {diff:e}");
            let (mut ra, mut rb) = (vec![0.; size], vec![0.; size]);
            rust.accumulate_density_pair(&c, [0.3, 0.7], 10., &mut ra)
                .unwrap();
            mkl.accumulate_density_pair(&c, [0.3, 0.7], 10., &mut rb)
                .unwrap();
            let norm = ra.iter().map(|x| x * x).sum::<f64>().sqrt();
            let diff = ra
                .iter()
                .zip(&rb)
                .map(|(x, y)| (x - y) * (x - y))
                .sum::<f64>()
                .sqrt();
            assert!(diff <= 1e-14 * norm, "{shape:?} density {diff:e}");
        }
        // The switch is read once per process; tests in this module do not set it.
        if crate::switches::text(crate::switches::Id::PAIR_FFT).is_none() {
            assert!(matches!(
                MklPairPlan::new([6, 6, 6], &[[0, 0, 0]]).unwrap(),
                MklPairPlan::Grid(_)
            ));
        }
    }
}
