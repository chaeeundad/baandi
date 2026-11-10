//! HAPPLY-FOLD (physics switch `BAANDI_HAPPLY_FOLD=on`): the local potential
//! of ONE real Gamma column through a half-length complex 3D transform.
//!
//! A real field psi on an `nx x ny x nz` grid with even `nz` is folded along z
//! into the complex half grid `Z(x,y,m) = psi(x,y,2m) + i psi(x,y,2m+1)`,
//! `M = nz/2`. With `w = exp(2 pi i / nz)` and P the full-grid coefficients,
//! the unscaled inverse transform of the half grid reproduces Z when its
//! spectrum is `Zhat(kx,ky,kz mod M) = sum over the two full kz of
//! P(k) (1 + i w^kz)`. After the multiplication by the real potential
//! (`V(2m)` on the real part, `V(2m+1)` on the imaginary part; in C order these
//! are the potential entries `2j` and `2j+1` of half-grid point `j`), the
//! forward transform (scaled by 1/(N/2)) gives `A = Zhat'(k')` and with the
//! mirror `B = conj(Zhat'(-k'))` the full coefficient of the product is
//! `F(k) = ((A + B) - i w^-kz (A - B)) / 4`. This is the textbook
//! two-real-halves identity (even/odd samples), written here from the DFT
//! definition. The full-grid pair path spends a whole complex transform on a
//! single real column with an exact zero partner; the fold spends half of it.
//!
//! Rounding class (different transform lengths and summation order), serial
//! per call, so results do not depend on the thread budget. Only representative
//! (plus) outputs are formed; the product of real fields is Hermitian, which is
//! what the pair path's even part extracts as well.
use crate::{Fft3d, FftError, fused_pair::PairMap};
use ndarray::Array3;
use num_complex::Complex64 as C;

/// Construction-sized folded single-column plan; no allocation on application.
pub struct FoldedRealPlan<F: Fft3d> {
    fft: F,
    full: usize,
    grid: Array3<C>,
    scratch: Array3<C>,
    plus_index: Vec<usize>,
    minus_index: Vec<usize>,
    plus_in: Vec<C>,
    minus_in: Vec<C>,
    plus_out: Vec<C>,
    len: usize,
}

/// Half-grid shape for a full grid, or `None` when the last axis is odd.
pub fn folded_shape(shape: [usize; 3]) -> Option<[usize; 3]> {
    (shape[2] >= 2 && shape[2].is_multiple_of(2)).then(|| [shape[0], shape[1], shape[2] / 2])
}

impl<F: Fft3d> FoldedRealPlan<F> {
    /// `fft` must have shape `folded_shape(shape)`; `indices`, `plus`, `minus`
    /// are the ordered full-sphere indices and the HalfSphere mirror tables
    /// (validated through `PairMap::new`).
    pub fn from_fft(
        fft: F,
        shape: [usize; 3],
        indices: &[[usize; 3]],
        plus: &[u32],
        minus: &[u32],
    ) -> Result<Self, FftError> {
        let half = folded_shape(shape).ok_or(FftError::InvalidGrid)?;
        if fft.shape() != half {
            return Err(FftError::IncompatibleArray);
        }
        let map = PairMap::new(shape, indices, plus, minus)?;
        let [_, ny, nz] = shape;
        let m = half[2];
        let at = |p: [usize; 3]| (p[0] * ny + p[1]) * m + p[2] % m;
        let angle = |p: [usize; 3]| 2. * std::f64::consts::PI * p[2] as f64 / nz as f64;
        let inward = |p: [usize; 3]| {
            let (s, c) = angle(p).sin_cos();
            // 1 + i w^z
            C::new(1. - s, c)
        };
        let outward = |p: [usize; 3]| {
            let (s, c) = angle(p).sin_cos();
            // -i w^-z = -i (cos - i sin)
            C::new(-s, -c)
        };
        let point = |i: u32| indices[i as usize];
        Ok(Self {
            fft,
            full: shape.iter().product(),
            grid: Array3::zeros(half),
            scratch: Array3::zeros(half),
            plus_index: plus.iter().map(|&i| at(point(i))).collect(),
            minus_index: minus.iter().map(|&i| at(point(i))).collect(),
            plus_in: plus.iter().map(|&i| inward(point(i))).collect(),
            minus_in: minus.iter().map(|&i| inward(point(i))).collect(),
            plus_out: plus.iter().map(|&i| outward(point(i))).collect(),
            len: map.len(),
        })
    }

    /// Overwrite `out` with the ORC coordinates of P V psi for one real ORC
    /// column `a` (same convention as `PairPlan::apply_local_pair` with a zero
    /// partner). Potential in Ha on the full grid, C order.
    pub fn apply_local_real(
        &mut self,
        a: &[f64],
        potential: &[f64],
        out: &mut [f64],
    ) -> Result<(), FftError> {
        if a.len() != self.len || out.len() != self.len || potential.len() != self.full {
            return Err(FftError::IncompatibleArray);
        }
        let sqrt2 = std::f64::consts::SQRT_2;
        {
            let g = self
                .grid
                .as_slice_mut()
                .ok_or(FftError::IncompatibleArray)?;
            g.fill(C::default());
            g[self.plus_index[0]] += C::new(a[0], 0.) * self.plus_in[0];
            for i in 1..self.plus_index.len() {
                let c = C::new(a[2 * i - 1], a[2 * i]) / sqrt2;
                g[self.plus_index[i]] += c * self.plus_in[i];
                g[self.minus_index[i]] += c.conj() * self.minus_in[i];
            }
        }
        self.fft
            .transform_in_place(&mut self.grid, &mut self.scratch, true)?;
        {
            let g = self
                .grid
                .as_slice_mut()
                .ok_or(FftError::IncompatibleArray)?;
            for (z, v) in g.iter_mut().zip(potential.as_chunks::<2>().0) {
                *z = C::new(z.re * v[0], z.im * v[1]);
            }
        }
        self.fft
            .transform_in_place(&mut self.grid, &mut self.scratch, false)?;
        let g = self.grid.as_slice().ok_or(FftError::IncompatibleArray)?;
        let value = |i: usize| {
            let a = g[self.plus_index[i]];
            let b = g[self.minus_index[i]].conj();
            ((a + b) + self.plus_out[i] * (a - b)) * 0.25
        };
        out[0] = value(0).re;
        for i in 1..self.plus_index.len() {
            let f = value(i);
            out[2 * i - 1] = f.re * sqrt2;
            out[2 * i] = f.im * sqrt2;
        }
        Ok(())
    }
}

impl FoldedRealPlan<crate::FftPlan3d> {
    /// Selected backend (MKL DFTI in `mkl` builds, rustfft otherwise);
    /// `Ok(None)` when the last grid axis is odd (caller keeps the pair path).
    pub fn new(
        shape: [usize; 3],
        indices: &[[usize; 3]],
        plus: &[u32],
        minus: &[u32],
    ) -> Result<Option<Self>, FftError> {
        match folded_shape(shape) {
            None => Ok(None),
            Some(half) => {
                Self::from_fft(crate::FftPlan3d::new(half)?, shape, indices, plus, minus).map(Some)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{PairPlan, fused::OrderedFusedPlan};

    /// Synthetic sphere with G=0 first, then (plus, minus) mirror pairs.
    pub(crate) fn sphere(shape: [usize; 3], r2: i64) -> (Vec<[usize; 3]>, Vec<u32>, Vec<u32>) {
        let wrap = |v: i64, n: usize| v.rem_euclid(n as i64) as usize;
        let mut indices = vec![[0; 3]];
        let (mut plus, mut minus) = (vec![0u32], vec![0u32]);
        let h = |n: usize| (n as i64 - 1) / 2;
        for x in -h(shape[0])..=h(shape[0]) {
            for y in -h(shape[1])..=h(shape[1]) {
                for z in -h(shape[2])..=h(shape[2]) {
                    if x * x + y * y + z * z > r2 || (x, y, z) <= (0, 0, 0) {
                        continue;
                    }
                    plus.push(indices.len() as u32);
                    indices.push([wrap(x, shape[0]), wrap(y, shape[1]), wrap(z, shape[2])]);
                    minus.push(indices.len() as u32);
                    indices.push([wrap(-x, shape[0]), wrap(-y, shape[1]), wrap(-z, shape[2])]);
                }
            }
        }
        (indices, plus, minus)
    }

    /// Mac proxy for the MKL path (full-grid `GridFusedPlan` with a zero
    /// partner against the fold on the same 3D FFT implementation, rustfft
    /// here, MKL DFTI on dev08). Observation only.
    /// cargo test --release -p dft-numerics --lib fold_full_grid_proxy -- --ignored --nocapture
    #[test]
    #[ignore]
    fn fold_full_grid_proxy() {
        use crate::fused::GridFusedPlan;
        use std::time::Instant;
        for (shape, r2) in [
            ([48, 48, 48], 121),
            ([60, 60, 60], 196),
            ([72, 72, 72], 289),
        ] {
            let (indices, plus, minus) = sphere(shape, r2);
            let n = indices.len();
            let size: usize = shape.iter().product();
            let v: Vec<f64> = (0..size)
                .map(|i| 0.3 * ((i * 7 + 1) as f64).sin())
                .collect();
            let a: Vec<f64> = (0..n).map(|i| ((i * 13 + 5) as f64).cos()).collect();
            let map = PairMap::new(shape, &indices, &plus, &minus).unwrap();
            let mut pair = PairPlan::new(
                GridFusedPlan::from_fft(crate::RustFftPlan3d::new(shape).unwrap(), &indices)
                    .unwrap(),
                map,
            )
            .unwrap();
            let mut fold = FoldedRealPlan::<crate::RustFftPlan3d>::from_fft(
                crate::RustFftPlan3d::new(folded_shape(shape).unwrap()).unwrap(),
                shape,
                &indices,
                &plus,
                &minus,
            )
            .unwrap();
            let (mut x, mut y) = (vec![0.; n], vec![0.; n]);
            let (mut tp, mut tf) = (Vec::new(), Vec::new());
            for round in 0..15 {
                for k in 0..2 {
                    let start = Instant::now();
                    if (k + round) % 2 == 0 {
                        pair.apply_local_pair(&a, None, &v, &mut x, None).unwrap();
                        tp.push(start.elapsed().as_secs_f64() * 1e6);
                    } else {
                        fold.apply_local_real(&a, &v, &mut y).unwrap();
                        tf.push(start.elapsed().as_secs_f64() * 1e6);
                    }
                }
            }
            tp.sort_by(f64::total_cmp);
            tf.sort_by(f64::total_cmp);
            println!(
                "{shape:?} rows={n} full-grid pair={:.0}us fold={:.0}us ratio={:.3}",
                tp[7],
                tf[7],
                tf[7] / tp[7]
            );
        }
    }

    #[test]
    fn fold_matches_pair_with_zero_partner() {
        for (shape, r2) in [
            ([6, 6, 6], 4),
            ([5, 7, 8], 6),
            ([8, 6, 10], 9),
            ([9, 9, 4], 3),
            ([4, 6, 12], 16),
        ] {
            let (indices, plus, minus) = sphere(shape, r2);
            let n = indices.len();
            let size: usize = shape.iter().product();
            let v: Vec<f64> = (0..size)
                .map(|i| 0.3 * ((i * 7 + 1) as f64).sin())
                .collect();
            let a: Vec<f64> = (0..n).map(|i| ((i * 13 + 5) as f64).cos()).collect();
            let map = PairMap::new(shape, &indices, &plus, &minus).unwrap();
            let mut pair =
                PairPlan::new(OrderedFusedPlan::new(shape, &indices).unwrap(), map).unwrap();
            let mut expected = vec![0.; n];
            pair.apply_local_pair(&a, None, &v, &mut expected, None)
                .unwrap();
            let mut fold = FoldedRealPlan::<crate::RustFftPlan3d>::from_fft(
                crate::RustFftPlan3d::new(folded_shape(shape).unwrap()).unwrap(),
                shape,
                &indices,
                &plus,
                &minus,
            )
            .unwrap();
            let mut actual = vec![0.; n];
            for _ in 0..2 {
                fold.apply_local_real(&a, &v, &mut actual).unwrap();
                let diff = actual
                    .iter()
                    .zip(&expected)
                    .map(|(x, y)| (x - y) * (x - y))
                    .sum::<f64>()
                    .sqrt();
                let norm = expected.iter().map(|x| x * x).sum::<f64>().sqrt();
                assert!(norm > 1e-3 && diff <= 1e-14 * norm, "{shape:?}: {diff:e}");
            }
            assert!(fold.apply_local_real(&a[1..], &v, &mut actual).is_err());
            assert!(fold.apply_local_real(&a, &v[1..], &mut actual).is_err());
        }
        assert!(folded_shape([4, 4, 5]).is_none());
        let (indices, plus, minus) = sphere([4, 4, 6], 2);
        assert!(
            FoldedRealPlan::<crate::RustFftPlan3d>::from_fft(
                crate::RustFftPlan3d::new([4, 4, 6]).unwrap(),
                [4, 4, 6],
                &indices,
                &plus,
                &minus
            )
            .is_err()
        );
    }
}
