//! Real Gamma pairs on unchanged full-sphere FFT plans (design section 6).
//! The caller supplies `HalfSphere::plus/minus`; this layer cannot depend on
//! dft-basis, which already depends on dft-numerics. No engine dispatch lives here.
use crate::{FftError, fused::FusedOperator, pool::WorkerPool};
use num_complex::Complex64 as C;
use std::sync::Arc;

/// Additive backend operation: one inverse transform, no forward transform.
pub trait PairDensityOperator: FusedOperator {
    /// Add `(weights[0]*Re(psi)^2 + weights[1]*Im(psi)^2)/volume`.
    /// Coefficients describe the combined real functions `psi_a + i psi_b`.
    /// Grid order and units match `FusedOperator::accumulate_density`.
    fn accumulate_density_pair(
        &mut self,
        phi: &[C],
        weights: [f64; 2],
        volume: f64,
        density: &mut [f64],
    ) -> Result<(), FftError>;
}

/// DENSITY-PAIR (engine switch `BAANDI_DENSITY_PAIR=on`): one inverse transform
/// of the combined coefficients `c_a + i c_b` of two real functions, with each
/// band's density added into its own grid: `a += weights[0] Re(psi)^2 / volume`,
/// `b += weights[1] Im(psi)^2 / volume`. Grid order and units match
/// `FusedOperator::accumulate_density`; per band the value equals the
/// single-band kernel to rounding (the dropped part is rounding noise).
pub trait SplitPairDensity: FusedOperator {
    fn accumulate_density_pair_split(
        &mut self,
        phi: &[C],
        weights: [f64; 2],
        volume: f64,
        a: &mut [f64],
        b: &mut [f64],
    ) -> Result<(), FftError>;
}

/// Validated full-sphere mirror table, independent of the basis crate.
#[derive(Clone)]
pub struct PairMap {
    shape: [usize; 3],
    plus: Vec<u32>,
    minus: Vec<u32>,
    len: usize,
}
impl PairMap {
    /// Validate coverage, unique grid locations, G=0 and FFT inversion.
    pub fn new(
        shape: [usize; 3],
        indices: &[[usize; 3]],
        plus: &[u32],
        minus: &[u32],
    ) -> Result<Self, FftError> {
        crate::SpherePattern::from_indices(shape, indices)?;
        let n = indices.len();
        if plus.is_empty()
            || plus.len() != minus.len()
            || plus.len().checked_mul(2).and_then(|n| n.checked_sub(1)) != Some(n)
            || plus[0] != 0
            || minus[0] != 0
            || indices[0] != [0; 3]
        {
            return Err(FftError::IncompatibleArray);
        }
        let mut seen = vec![false; n];
        seen[0] = true;
        let mut grid = std::collections::HashSet::new();
        for &index in indices {
            if !grid.insert(index) {
                return Err(FftError::IncompatibleArray);
            }
        }
        for (&p, &m) in plus.iter().zip(minus).skip(1) {
            let (p, m) = (p as usize, m as usize);
            if p >= n
                || m >= n
                || p == m
                || seen[p]
                || seen[m]
                || (0..3).any(|d| indices[m][d] != (shape[d] - indices[p][d]) % shape[d])
            {
                return Err(FftError::IncompatibleArray);
            }
            seen[p] = true;
            seen[m] = true;
        }
        Ok(Self {
            shape,
            plus: plus.to_vec(),
            minus: minus.to_vec(),
            len: n,
        })
    }
    /// Real coordinate count, equal to the full-sphere complex count.
    pub fn len(&self) -> usize {
        self.len
    }
    /// Valid maps always contain G=0.
    pub fn is_empty(&self) -> bool {
        false
    }
    fn valid_bands(&self, a: &[f64], b: Option<&[f64]>) -> bool {
        a.len() == self.len && b.is_none_or(|b| b.len() == self.len)
    }
    /// Pack ORC bands; `None` is the exact zero partner for an odd tail.
    /// Finite input validation belongs to the physics caller. Zero vectors are
    /// exempt from the norm-ratio assertion; all nonzero partners must be balanced.
    pub fn pack(&self, a: &[f64], b: Option<&[f64]>, phi: &mut [C]) -> Result<(), FftError> {
        if !self.valid_bands(a, b) || phi.len() != self.len {
            return Err(FftError::IncompatibleArray);
        }
        #[cfg(debug_assertions)]
        if let Some(b) = b {
            let na = a.iter().fold(0.0_f64, |n, &x| n.hypot(x));
            let nb = b.iter().fold(0.0_f64, |n, &x| n.hypot(x));
            debug_assert!(
                na == 0. || nb == 0. || na.max(nb) / na.min(nb) <= 1e3,
                "paired nonzero band norms must be within a factor 1e3"
            );
        }
        let b_at = |i| b.map_or(0., |b| b[i]);
        phi[0] = C::new(a[0], b_at(0));
        for (i, (&p, &m)) in self.plus.iter().zip(&self.minus).enumerate().skip(1) {
            let ca = C::new(a[2 * i - 1], a[2 * i]) / std::f64::consts::SQRT_2;
            let cb = C::new(b_at(2 * i - 1), b_at(2 * i)) / std::f64::consts::SQRT_2;
            phi[p as usize] = ca + C::new(-cb.im, cb.re);
            phi[m as usize] = ca.conj() + C::new(cb.im, cb.re);
        }
        Ok(())
    }
    /// Separate by parity. G=0 stores only Re(h) for a and Im(h) for b;
    /// each recovered band's imaginary G=0 leftover is discarded.
    pub fn unpack(
        &self,
        h: &[C],
        a: &mut [f64],
        mut b: Option<&mut [f64]>,
    ) -> Result<(), FftError> {
        if h.len() != self.len || !self.valid_bands(a, b.as_deref()) {
            return Err(FftError::IncompatibleArray);
        }
        a[0] = h[0].re;
        if let Some(b) = b.as_deref_mut() {
            b[0] = h[0].im;
        }
        for (i, (&p, &m)) in self.plus.iter().zip(&self.minus).enumerate().skip(1) {
            let even = (h[p as usize] + h[m as usize].conj()) * 0.5;
            let odd = (h[p as usize] - h[m as usize].conj()) * 0.5;
            a[2 * i - 1] = even.re * std::f64::consts::SQRT_2;
            a[2 * i] = even.im * std::f64::consts::SQRT_2;
            if let Some(b) = b.as_deref_mut() {
                b[2 * i - 1] = odd.im * std::f64::consts::SQRT_2;
                b[2 * i] = -odd.re * std::f64::consts::SQRT_2;
            }
        }
        Ok(())
    }
}

/// Persistent pair scratch. All numerical allocation occurs at construction.
pub struct PairPlan<P> {
    plan: P,
    map: PairMap,
    phi: Vec<C>,
    image: Vec<C>,
}
impl<P: FusedOperator> PairPlan<P> {
    /// The supplied plan must use the same ordered indices as `map`.
    pub fn new(plan: P, map: PairMap) -> Result<Self, FftError> {
        if plan.coefficient_count() != map.len || plan.shape() != map.shape {
            return Err(FftError::IncompatibleArray);
        }
        let n = map.len;
        Ok(Self {
            plan,
            map,
            phi: vec![C::default(); n],
            image: vec![C::default(); n],
        })
    }
    /// One inverse and one forward transform for two real ORC columns.
    pub fn apply_local_pair(
        &mut self,
        a: &[f64],
        b: Option<&[f64]>,
        potential: &[f64],
        out_a: &mut [f64],
        out_b: Option<&mut [f64]>,
    ) -> Result<(), FftError> {
        if !self.map.valid_bands(out_a, out_b.as_deref()) {
            return Err(FftError::IncompatibleArray);
        }
        self.map.pack(a, b, &mut self.phi)?;
        self.plan
            .apply_local(&self.phi, potential, &mut self.image)?;
        self.map.unpack(&self.image, out_a, out_b)
    }
    /// Column-major block, pairs (0,1), (2,3), ...; odd tail uses a zero partner.
    pub fn apply_block(
        &mut self,
        bands: &[f64],
        potential: &[f64],
        out: &mut [f64],
    ) -> Result<(), FftError> {
        let n = self.map.len;
        if !bands.len().is_multiple_of(n)
            || out.len() != bands.len()
            || potential.len() != self.plan.shape().iter().product::<usize>()
        {
            return Err(FftError::IncompatibleArray);
        }
        for (pair, dst) in bands.chunks(2 * n).zip(out.chunks_mut(2 * n)) {
            let (a, b) = pair.split_at(n);
            let (oa, ob) = dst.split_at_mut(n);
            self.apply_local_pair(
                a,
                (!b.is_empty()).then_some(b),
                potential,
                oa,
                (!ob.is_empty()).then_some(ob),
            )?;
        }
        Ok(())
    }
}
impl<P: PairDensityOperator> PairPlan<P> {
    /// Add density from one inverse transform; skip only two exactly zero weights.
    pub fn accumulate_density_pair(
        &mut self,
        a: &[f64],
        b: Option<&[f64]>,
        weights: [f64; 2],
        volume: f64,
        density: &mut [f64],
    ) -> Result<(), FftError> {
        if !self.map.valid_bands(a, b) {
            return Err(FftError::IncompatibleArray);
        }
        validate_density(weights, volume)?;
        if density.len() != self.plan.shape().iter().product::<usize>() {
            return Err(FftError::IncompatibleArray);
        }
        if weights == [0.; 2] {
            return Ok(());
        }
        self.map.pack(a, b, &mut self.phi)?;
        self.plan
            .accumulate_density_pair(&self.phi, weights, volume, density)
    }
}

pub(crate) fn validate_density(weights: [f64; 2], volume: f64) -> Result<(), FftError> {
    if !volume.is_finite() || volume <= 0. || weights.iter().any(|w| !w.is_finite()) {
        Err(FftError::IncompatibleArray)
    } else {
        Ok(())
    }
}

// Four pairs per piece is a fixed arithmetic grouping, not a tuned dispatch
// threshold. Storage scales with ceil(bands/8), independently of worker count.
const BANDS_PER_PIECE: usize = 8;
struct DensityPiece<P> {
    pair: PairPlan<P>,
    density: Vec<f64>,
}
/// Reusable deterministic density workspace. Plans and partial grids are owned
/// per fixed piece; reduction is ascending piece order for every thread budget.
pub struct PairDensityWorkspace<P> {
    pieces: Vec<DensityPiece<P>>,
    pool: Arc<WorkerPool>,
    bands: usize,
    n: usize,
    grid_len: usize,
}
impl<P: PairDensityOperator + Send> PairDensityWorkspace<P> {
    /// `make_plan` constructs identical ordered full-sphere plans.
    pub fn new(
        map: PairMap,
        bands: usize,
        pool: Arc<WorkerPool>,
        mut make_plan: impl FnMut() -> Result<P, FftError>,
    ) -> Result<Self, FftError> {
        let first = make_plan()?;
        let shape = first.shape();
        let grid_len = shape.iter().product();
        let n = map.len;
        let mut pieces = Vec::with_capacity(bands.div_ceil(BANDS_PER_PIECE));
        let mut first = Some(first);
        for _ in 0..bands.div_ceil(BANDS_PER_PIECE) {
            let plan = match first.take() {
                Some(p) => p,
                None => make_plan()?,
            };
            if plan.shape() != shape {
                return Err(FftError::IncompatibleArray);
            }
            pieces.push(DensityPiece {
                pair: PairPlan::new(plan, map.clone())?,
                density: vec![0.; grid_len],
            });
        }
        Ok(Self {
            pieces,
            pool,
            bands,
            n,
            grid_len,
        })
    }
    /// Overwrite density in e/Bohr^3. No allocations or worker creation here.
    pub fn density(
        &mut self,
        bands: &[f64],
        weights: &[f64],
        volume: f64,
        out: &mut [f64],
    ) -> Result<(), FftError> {
        if self.bands.checked_mul(self.n) != Some(bands.len())
            || weights.len() != self.bands
            || out.len() != self.grid_len
            || weights.iter().any(|w| !w.is_finite())
        {
            return Err(FftError::IncompatibleArray);
        }
        validate_density([0.; 2], volume)?;
        self.pool
            .install(|| density_pieces(&mut self.pieces, bands, weights, self.n, volume))?;
        out.fill(0.);
        for piece in &self.pieces {
            for (dst, &src) in out.iter_mut().zip(&piece.density) {
                *dst += src;
            }
        }
        Ok(())
    }
}
fn density_pieces<P: PairDensityOperator + Send>(
    pieces: &mut [DensityPiece<P>],
    bands: &[f64],
    weights: &[f64],
    n: usize,
    volume: f64,
) -> Result<(), FftError> {
    if pieces.len() > 1 {
        let mid = pieces.len() / 2;
        let (left, right) = pieces.split_at_mut(mid);
        let split = mid * BANDS_PER_PIECE;
        let (lb, rb) = bands.split_at(split * n);
        let (lw, rw) = weights.split_at(split);
        let (a, b) = WorkerPool::join(
            || density_pieces(left, lb, lw, n, volume),
            || density_pieces(right, rb, rw, n, volume),
        );
        a?;
        b
    } else if let Some(piece) = pieces.first_mut() {
        piece.density.fill(0.);
        for (pair, w) in bands.chunks(2 * n).zip(weights.chunks(2)) {
            let (a, b) = pair.split_at(n);
            piece.pair.accumulate_density_pair(
                a,
                (!b.is_empty()).then_some(b),
                [w[0], *w.get(1).unwrap_or(&0.)],
                volume,
                &mut piece.density,
            )?;
        }
        Ok(())
    } else {
        Ok(())
    }
}

/// Pool-aware local application with one fixed band pair per piece. Each piece
/// owns its FFT and pack/unpack scratch, independent of the worker count.
/// Construction reserves the maximum block; applications allocate no storage.
pub struct PairLocalWorkspace<P> {
    pieces: Vec<PairPlan<P>>,
    pool: Arc<WorkerPool>,
    capacity: usize,
    n: usize,
    grid_len: usize,
}
impl<P: FusedOperator + Send> PairLocalWorkspace<P> {
    /// Plans must use the ordered full-sphere indices of `map`.
    pub fn new(
        map: PairMap,
        capacity: usize,
        pool: Arc<WorkerPool>,
        mut make_plan: impl FnMut() -> Result<P, FftError>,
    ) -> Result<Self, FftError> {
        let n = map.len;
        if capacity.checked_mul(n).is_none() {
            return Err(FftError::IncompatibleArray);
        }
        let grid_len = map.shape.iter().product();
        let mut pieces = Vec::with_capacity(capacity.div_ceil(2));
        for _ in 0..capacity.div_ceil(2) {
            pieces.push(PairPlan::new(make_plan()?, map.clone())?);
        }
        Ok(Self {
            pieces,
            pool,
            capacity,
            n,
            grid_len,
        })
    }

    /// Preserve serial pair arithmetic, including the zero partner of an odd
    /// tail. Disjoint output slices need neither locks nor a reduction.
    pub fn apply_block(
        &mut self,
        bands: &[f64],
        potential: &[f64],
        out: &mut [f64],
    ) -> Result<(), FftError> {
        if !bands.len().is_multiple_of(self.n)
            || bands.len() / self.n > self.capacity
            || out.len() != bands.len()
            || potential.len() != self.grid_len
        {
            return Err(FftError::IncompatibleArray);
        }
        let count = (bands.len() / self.n).div_ceil(2);
        let pieces = &mut self.pieces[..count];
        self.pool
            .install(|| local_pieces(pieces, bands, potential, out, self.n))
    }
}
fn local_pieces<P: FusedOperator + Send>(
    pieces: &mut [PairPlan<P>],
    bands: &[f64],
    potential: &[f64],
    out: &mut [f64],
    n: usize,
) -> Result<(), FftError> {
    if pieces.len() > 1 {
        let mid = pieces.len() / 2;
        let (left, right) = pieces.split_at_mut(mid);
        let (lb, rb) = bands.split_at(mid * 2 * n);
        let (lo, ro) = out.split_at_mut(mid * 2 * n);
        let (a, b) = WorkerPool::join(
            || local_pieces(left, lb, potential, lo, n),
            || local_pieces(right, rb, potential, ro, n),
        );
        a?;
        b
    } else if let Some(piece) = pieces.first_mut() {
        piece.apply_block(bands, potential, out)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fused::OrderedFusedPlan;

    #[test]
    fn reusable_pair_storage() {
        let indices = [[0, 0, 0], [1, 0, 0], [4, 0, 0]];
        let map = PairMap::new([5, 1, 1], &indices, &[0, 1], &[0, 2]).unwrap();
        let mut ws = PairDensityWorkspace::new(map, 17, WorkerPool::new(4).unwrap(), || {
            OrderedFusedPlan::new([5, 1, 1], &indices)
        })
        .unwrap();
        // Stable pointer/length/capacity checks for every new numerical buffer.
        // Backend process_with_scratch allocation contracts are inherited from
        // the existing fused plans; no allocator instrumentation/unsafe added.
        fn state(ws: &PairDensityWorkspace<OrderedFusedPlan>) -> Vec<(usize, usize, usize)> {
            ws.pieces
                .iter()
                .flat_map(|p| {
                    [
                        (
                            p.pair.phi.as_ptr() as usize,
                            p.pair.phi.len(),
                            p.pair.phi.capacity(),
                        ),
                        (
                            p.pair.image.as_ptr() as usize,
                            p.pair.image.len(),
                            p.pair.image.capacity(),
                        ),
                        (
                            p.density.as_ptr() as usize,
                            p.density.len(),
                            p.density.capacity(),
                        ),
                        (
                            p.pair.map.plus.as_ptr() as usize,
                            p.pair.map.plus.len(),
                            p.pair.map.plus.capacity(),
                        ),
                        (
                            p.pair.map.minus.as_ptr() as usize,
                            p.pair.map.minus.len(),
                            p.pair.map.minus.capacity(),
                        ),
                    ]
                })
                .collect()
        }
        let initial = state(&ws);
        let bands = [0.2; 51];
        let weights = [0.5; 17];
        let mut density = [0.; 5];
        let mut out = [0.; 51];
        for _ in 0..10 {
            ws.density(&bands, &weights, 10., &mut density).unwrap();
            ws.pieces[0]
                .pair
                .apply_block(&bands, &[0.7; 5], &mut out)
                .unwrap();
            assert_eq!(state(&ws), initial);
        }
    }
}
