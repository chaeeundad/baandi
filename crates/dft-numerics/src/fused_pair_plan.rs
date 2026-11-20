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

/// Lane kernels chosen at run time (MPI M2 distributed FFT, `BAANDI_M2_FFT`):
/// rustfft, or MKL DFTI batch lanes when `vendor` is requested on an MKL build.
/// A pure build has no vendor lanes and always uses rustfft (bit-identical to
/// `RustLanes`).
pub enum BackendLanes {
    Rust(RustLanes),
    #[cfg(feature = "mkl")]
    Mkl(MklLanes),
}
impl BackendLanes {
    /// Plan lanes of lengths `n` with `counts[a]` lanes per buffer of axis `a`.
    pub fn plan(n: [usize; 3], counts: [usize; 3], vendor: bool) -> Result<Self, FftError> {
        #[cfg(feature = "mkl")]
        if vendor {
            return Ok(Self::Mkl(MklLanes::plan(n, counts)?));
        }
        let _ = vendor;
        Ok(Self::Rust(RustLanes::plan(n, counts)?))
    }
    /// True when the MKL lanes are active.
    pub fn is_vendor(&self) -> bool {
        !matches!(self, Self::Rust(_))
    }
    /// Unscaled in-place transforms of all lanes of the axis buffer `data`.
    pub fn run(&mut self, axis: usize, inverse: bool, data: &mut [C]) -> Result<(), FftError> {
        match self {
            Self::Rust(l) => l.run(axis, inverse, data),
            #[cfg(feature = "mkl")]
            Self::Mkl(l) => l.run(axis, inverse, data),
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
    profile: crate::real::fft_profile::Probe,
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
            profile: crate::real::fft_profile::Probe::new(n),
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
        let mark = self.profile.start();
        self.sticks.fill(C::default());
        for (&i, &v) in self.scatter.iter().zip(c) {
            self.sticks[i] = v;
        }
        self.profile.finish(
            mark,
            0,
            self.scatter.len(),
            16 * self.sticks.len() + 32 * self.scatter.len(),
            None,
        );
        if !self.sticks.is_empty() {
            let mark = self.profile.start();
            self.lanes.run(0, true, &mut self.sticks)?;
            self.profile
                .finish(mark, 1, self.sticks.len(), 32 * self.sticks.len(), Some(0));
        }
        Ok(())
    }
    fn inverse_slab(&mut self, x: usize) -> Result<(), FftError> {
        let [nx, ny, nz] = self.n;
        let mark = self.profile.start();
        self.middle.fill(C::default());
        for &(i, slot) in &self.occupied {
            self.middle[i] = self.sticks[slot * nx + x];
        }
        self.profile.finish(
            mark,
            4,
            self.occupied.len(),
            16 * self.middle.len() + 32 * self.occupied.len(),
            None,
        );
        if !self.middle.is_empty() {
            let mark = self.profile.start();
            self.lanes.run(1, true, &mut self.middle)?;
            self.profile
                .finish(mark, 5, self.middle.len(), 32 * self.middle.len(), Some(1));
        }
        let mark = self.profile.start();
        self.slab.fill(C::default());
        for (k, &z) in self.planes.iter().enumerate() {
            for y in 0..ny {
                self.slab[y * nz + z] = self.middle[k * ny + y];
            }
        }
        self.profile.finish(
            mark,
            6,
            self.middle.len(),
            16 * self.slab.len() + 32 * self.middle.len(),
            None,
        );
        let mark = self.profile.start();
        let result = self.lanes.run(2, true, &mut self.slab);
        self.profile
            .finish(mark, 7, self.slab.len(), 32 * self.slab.len(), Some(2));
        result
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
        let mark = self.profile.start();
        self.lanes.run(2, false, &mut self.slab)?;
        self.profile
            .finish(mark, 9, self.slab.len(), 32 * self.slab.len(), Some(2));
        let mark = self.profile.start();
        for (k, &z) in self.planes.iter().enumerate() {
            for y in 0..ny {
                self.middle[k * ny + y] = self.slab[y * nz + z];
            }
        }
        self.profile
            .finish(mark, 10, self.middle.len(), 32 * self.middle.len(), None);
        if !self.middle.is_empty() {
            let mark = self.profile.start();
            self.lanes.run(1, false, &mut self.middle)?;
            self.profile
                .finish(mark, 11, self.middle.len(), 32 * self.middle.len(), Some(1));
        }
        let mark = self.profile.start();
        for &(i, slot) in &self.occupied {
            self.sticks[slot * nx + x] = self.middle[i];
        }
        self.profile.finish(
            mark,
            12,
            self.occupied.len(),
            32 * self.occupied.len(),
            None,
        );
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
        let _profile_scope = self.profile.apply_scope();
        self.validate(c, potential.len())?;
        if out.len() != c.len() {
            return Err(FftError::IncompatibleArray);
        }
        self.inverse_sticks(c)?;
        laps.lap(Phase::Inverse);
        for x in 0..self.n[0] {
            self.inverse_slab(x)?;
            laps.lap(Phase::Inverse);
            let mark = self.profile.start();
            for y in 0..self.n[1] {
                let start = self.grid_start(x, y);
                for (a, &v) in self.slab[y * self.n[2]..(y + 1) * self.n[2]]
                    .iter_mut()
                    .zip(&potential[start..start + self.n[2]])
                {
                    *a *= v;
                }
            }
            self.profile
                .finish(mark, 8, self.slab.len(), 40 * self.slab.len(), None);
            laps.lap(Phase::Potential);
            self.forward_slab(x)?;
            laps.lap(Phase::Forward);
        }
        if !self.sticks.is_empty() {
            let mark = self.profile.start();
            self.lanes.run(0, false, &mut self.sticks)?;
            self.profile
                .finish(mark, 15, self.sticks.len(), 32 * self.sticks.len(), Some(0));
        }
        let mark = self.profile.start();
        for (o, &i) in out.iter_mut().zip(&self.scatter) {
            *o = self.sticks[i] / self.size as f64;
        }
        self.profile
            .finish(mark, 16, self.scatter.len(), 32 * self.scatter.len(), None);
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
impl<L: PairLanes> crate::fused_pair::ForwardProduct for PairFusedPlan<L> {
    /// The forward half of `apply_local_with` on a kept C-order grid: slab by
    /// slab, `grid * potential`, then the same forward lanes, gather and 1/N.
    fn forward_product(
        &mut self,
        grid: &[C],
        potential: &[f64],
        out: &mut [C],
    ) -> Result<(), FftError> {
        let _profile_scope = self.profile.apply_scope();
        if out.len() != self.scatter.len()
            || grid.len() != self.size
            || potential.len() != self.size
        {
            return Err(FftError::IncompatibleArray);
        }
        let nz = self.n[2];
        for x in 0..self.n[0] {
            let mark = self.profile.start();
            for y in 0..self.n[1] {
                let start = self.grid_start(x, y);
                for ((a, &g), &v) in self.slab[y * nz..(y + 1) * nz]
                    .iter_mut()
                    .zip(&grid[start..start + nz])
                    .zip(&potential[start..start + nz])
                {
                    *a = g * v;
                }
            }
            self.profile
                .finish(mark, 8, self.slab.len(), 40 * self.slab.len(), None);
            self.forward_slab(x)?;
        }
        if !self.sticks.is_empty() {
            let mark = self.profile.start();
            self.lanes.run(0, false, &mut self.sticks)?;
            self.profile
                .finish(mark, 15, self.sticks.len(), 32 * self.sticks.len(), Some(0));
        }
        let mark = self.profile.start();
        for (o, &i) in out.iter_mut().zip(&self.scatter) {
            *o = self.sticks[i] / self.size as f64;
        }
        self.profile
            .finish(mark, 16, self.scatter.len(), 32 * self.scatter.len(), None);
        Ok(())
    }
}
#[cfg(feature = "mkl")]
impl crate::fused_pair::ForwardProduct for MklPairPlan {
    fn forward_product(
        &mut self,
        grid: &[C],
        potential: &[f64],
        out: &mut [C],
    ) -> Result<(), FftError> {
        match self {
            Self::Grid(p) => p.forward_product(grid, potential, out),
            Self::Pruned(p) => p.forward_product(grid, potential, out),
        }
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
    fn lane_fft_serial_bits_and_density_exclusion() {
        use crate::fused_pair::ForwardProduct;
        use crate::real::fft_profile;
        use crate::switches::RunConfig;
        for shape in [[4, 4, 4], [3, 5, 7], [2, 3, 12]] {
            let indices = vec![[0, 0, 0], [1, 0, 0], [0, 1, 1]];
            let c = vec![C::new(0.3, -0.2); indices.len()];
            let total = shape.iter().product();
            let v: Vec<_> = (0..total).map(|i| (i as f64).cos()).collect();
            let grid = vec![C::new(-0.2, 0.7); total];
            let mut reference = Vec::new();
            for arm in ["off", "on"] {
                let cfg = RunConfig::from_pairs(&[
                    ("fft_profile", arm),
                    ("phase_profile", "on"),
                    ("stage_timers", "off"),
                ])
                .unwrap();
                let _scope = cfg.enter();
                let mut plan = PairFusedPlan::new(shape, &indices).unwrap();
                let before = fft_profile::snapshot();
                let mut out = c.clone();
                plan.apply_local(&c, &v, &mut out).unwrap();
                let mut kept = c.clone();
                plan.forward_product(&grid, &v, &mut kept).unwrap();
                let delta = fft_profile::snapshot().since(before);
                let bits: Vec<_> = out
                    .iter()
                    .chain(&kept)
                    .map(|z| (z.re.to_bits(), z.im.to_bits()))
                    .collect();
                if arm == "off" {
                    reference = bits;
                    assert_eq!(delta, fft_profile::Snapshot::default());
                } else {
                    assert_eq!(bits, reference);
                    assert_eq!(delta.0[8][2], 2 * total as u64);
                    assert_eq!(delta.0[8][3], 80 * total as u64);
                    assert!(delta.0.iter().all(|v| v[0] == 0));
                }
                let before = fft_profile::snapshot();
                let mut density = vec![0.; total];
                plan.accumulate_density(&c, 1., 10., &mut density).unwrap();
                assert_eq!(
                    fft_profile::snapshot(),
                    before,
                    "density is outside the apply scope"
                );
            }
        }
    }

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

    /// REUSE-PSI: the kept inverse transform of the density plans (ordered and
    /// grid adapters) followed by `forward_product` equals `apply_local` of the
    /// pair plan and of the grid adapter (relative 1e-14), and keeping the
    /// grid leaves the split pair density bit-identical.
    #[test]
    fn kept_grid_forward_matches_apply_local() {
        use crate::fused::{GridFusedPlan, OrderedFusedPlan};
        use crate::fused_pair::{ForwardProduct, SplitPairDensity};
        for shape in [[15, 20, 24], [9, 12, 10]] {
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
            let relative = |a: &[C], b: &[C]| {
                let d = a
                    .iter()
                    .zip(b)
                    .map(|(x, y)| (x - y).norm_sqr())
                    .sum::<f64>();
                (d / b.iter().map(|y| y.norm_sqr()).sum::<f64>()).sqrt()
            };
            let mut pair = PairFusedPlan::new(shape, &indices).unwrap();
            let mut grid = GridFusedPlan::new(shape, &indices).unwrap();
            let mut expected = vec![C::default(); c.len()];
            pair.apply_local(&c, &v, &mut expected).unwrap();
            let mut ordered = OrderedFusedPlan::new(shape, &indices).unwrap();
            let mut keeps = Vec::new();
            for which in 0..2 {
                let (mut a, mut b) = (vec![0.; size], vec![0.; size]);
                let (mut a0, mut b0) = (vec![0.; size], vec![0.; size]);
                let mut keep = vec![C::default(); size];
                if which == 0 {
                    ordered
                        .accumulate_density_pair_split_keep(
                            &c,
                            [0.3, 0.7],
                            10.,
                            &mut a,
                            &mut b,
                            &mut keep,
                        )
                        .unwrap();
                    ordered
                        .accumulate_density_pair_split(&c, [0.3, 0.7], 10., &mut a0, &mut b0)
                        .unwrap();
                } else {
                    grid.accumulate_density_pair_split_keep(
                        &c,
                        [0.3, 0.7],
                        10.,
                        &mut a,
                        &mut b,
                        &mut keep,
                    )
                    .unwrap();
                    grid.accumulate_density_pair_split(&c, [0.3, 0.7], 10., &mut a0, &mut b0)
                        .unwrap();
                }
                assert_eq!((a, b), (a0, b0), "{shape:?} plan {which}: density changed");
                keeps.push(keep);
            }
            for keep in &keeps {
                let mut out = vec![C::default(); c.len()];
                pair.forward_product(keep, &v, &mut out).unwrap();
                let e = relative(&out, &expected);
                assert!(e <= 1e-14, "{shape:?} pair forward {e:e}");
                grid.forward_product(keep, &v, &mut out).unwrap();
                let e = relative(&out, &expected);
                assert!(e <= 1e-14, "{shape:?} grid forward {e:e}");
            }
            // The kept grid is untouched by the forward product.
            let copy = keeps[0].clone();
            let mut out = vec![C::default(); c.len()];
            pair.forward_product(&keeps[0], &v, &mut out).unwrap();
            assert_eq!(copy, keeps[0]);
        }
    }

    /// dev08 (MKL builds only): the mirror-slab plan on MKL DFTI lanes equals
    /// the rustfft lanes to relative 1e-14 (local and pair density), is
    /// repeatable bit for bit, and `MklPairPlan` selects the configured variant
    /// (MKL defaults to pruned; the explicit grid override remains available).
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
        // The registry default is pruned with MKL, independent of grid shape.
        // Scope the configuration so the caller's environment cannot skip or
        // redirect this check; exercise both explicit variants as well.
        for (pairs, pruned) in [
            (&[][..], true),
            (&[("pair_fft", "grid")][..], false),
            (&[("pair_fft", "pruned")][..], true),
        ] {
            let _config = crate::switches::RunConfig::from_pairs(pairs)
                .unwrap()
                .enter();
            let plan = MklPairPlan::new([6, 6, 6], &[[0, 0, 0]]).unwrap();
            assert_eq!(matches!(&plan, MklPairPlan::Pruned(_)), pruned);
            assert_eq!(matches!(&plan, MklPairPlan::Grid(_)), !pruned);
        }
    }
}
