//! Gamma mirror maps and orthonormal real coordinates (ORC).
use crate::{FftGridSpec, PlaneWaveBasis};
use num_complex::Complex64;
use std::collections::HashMap;
use thiserror::Error;

/// Absolute coefficient tolerance: allows roundoff at conversion boundaries,
/// while rejecting symmetry defects above 1e-14 in normalized orbitals.
pub const ORC_TOLERANCE: f64 = 1e-14;

/// Invalid Gamma basis, grid, or coefficient conversion.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum HalfSphereError {
    /// Only exactly Gamma bases have the required real structure.
    #[error("basis is not Gamma")]
    NotGamma,
    /// Missing zero, duplicate Miller index, missing mirror, or index overflow.
    #[error("basis is not an inversion-symmetric Gamma sphere")]
    InvalidBasis,
    /// The supplied FFT mapping does not describe this basis.
    #[error("FFT mapping does not match the mirror table")]
    InvalidGrid,
    /// Wrong vector length.
    #[error("incorrect coefficient count")]
    Dimension,
    /// Nonfinite coefficients, overflow, or broken Hermitian symmetry.
    #[error("coefficients are not finite and Hermitian within ORC_TOLERANCE")]
    InvalidCoefficients,
}

/// Representatives start with zero, then lexicographically positive Miller
/// triples in full-list order. Indices refer to the original full sphere.
#[derive(Debug, Clone)]
pub struct HalfSphere {
    plus: Vec<u32>,
    minus: Vec<u32>,
    mirror: Vec<u32>,
}
impl HalfSphere {
    /// Validate Gamma, inversion symmetry, and the supplied FFT mapping.
    pub fn new(basis: &PlaneWaveBasis, grid: &FftGridSpec) -> Result<Self, HalfSphereError> {
        if basis.k_fractional() != [0.; 3] {
            return Err(HalfSphereError::NotGamma);
        }
        let gs = basis.vectors();
        if gs.is_empty()
            || gs.len() > u32::MAX as usize
            || gs[0].miller != [0; 3]
            || gs[0].norm_squared != 0.0
        {
            return Err(HalfSphereError::InvalidBasis);
        }
        let map: HashMap<_, _> = gs
            .iter()
            .enumerate()
            .map(|(i, g)| (g.miller, i as u32))
            .collect();
        if map.len() != gs.len() {
            return Err(HalfSphereError::InvalidBasis);
        }
        let mut out = Self {
            plus: vec![0],
            minus: vec![0],
            mirror: Vec::with_capacity(gs.len()),
        };
        for (i, g) in gs.iter().enumerate() {
            let mut neg = [0; 3];
            for (n, m) in neg.iter_mut().zip(g.miller) {
                *n = m.checked_neg().ok_or(HalfSphereError::InvalidBasis)?;
            }
            let &j = map.get(&neg).ok_or(HalfSphereError::InvalidBasis)?;
            out.mirror.push(j);
            if g.miller > [0; 3] {
                out.plus.push(i as u32);
                out.minus.push(j);
            }
        }
        if grid.g_to_fft().len() != gs.len() {
            return Err(HalfSphereError::InvalidGrid);
        }
        for (i, g) in gs.iter().enumerate() {
            let expected = FftGridSpec::map_miller(g.miller, grid.dims())
                .map_err(|_| HalfSphereError::InvalidGrid)?;
            if grid.g_to_fft()[i] != expected
                || (0..3).any(|d| {
                    grid.g_to_fft()[out.mirror[i] as usize][d]
                        != (grid.dims()[d] - expected[d]) % grid.dims()[d]
                })
            {
                return Err(HalfSphereError::InvalidGrid);
            }
        }
        Ok(out)
    }
    /// Full-list representative indices, including zero.
    pub fn plus(&self) -> &[u32] {
        &self.plus
    }
    /// Full-list conjugate partner indices, including zero.
    pub fn minus(&self) -> &[u32] {
        &self.minus
    }
    /// Involution mapping each full index to its negative Miller index.
    pub fn mirror(&self) -> &[u32] {
        &self.mirror
    }
    /// Number of real coordinates (equal to the full complex basis length).
    pub fn full_len(&self) -> usize {
        self.mirror.len()
    }
    /// Pack into caller-owned storage. Validate before writing; accepted small
    /// symmetry defects are projected by retaining the positive representative.
    pub fn pack(&self, full: &[Complex64], real: &mut [f64]) -> Result<(), HalfSphereError> {
        if full.len() != self.full_len() || real.len() != self.full_len() {
            return Err(HalfSphereError::Dimension);
        }
        if full.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
            || full[0].im.abs() > ORC_TOLERANCE
        {
            return Err(HalfSphereError::InvalidCoefficients);
        }
        for (&p, &m) in self.plus.iter().zip(&self.minus).skip(1) {
            let z = full[p as usize];
            if (z - full[m as usize].conj()).norm() > ORC_TOLERANCE
                || !(z.re * std::f64::consts::SQRT_2).is_finite()
                || !(z.im * std::f64::consts::SQRT_2).is_finite()
            {
                return Err(HalfSphereError::InvalidCoefficients);
            }
        }
        real[0] = full[0].re;
        for (i, &p) in self.plus.iter().enumerate().skip(1) {
            real[2 * i - 1] = std::f64::consts::SQRT_2 * full[p as usize].re;
            real[2 * i] = std::f64::consts::SQRT_2 * full[p as usize].im;
        }
        Ok(())
    }
    /// Hermitian completion into caller-owned full-sphere storage, without allocation.
    pub fn unpack(&self, real: &[f64], full: &mut [Complex64]) -> Result<(), HalfSphereError> {
        if full.len() != self.full_len() || real.len() != self.full_len() {
            return Err(HalfSphereError::Dimension);
        }
        if real.iter().any(|x| !x.is_finite()) {
            return Err(HalfSphereError::InvalidCoefficients);
        }
        full[0] = Complex64::new(real[0], 0.);
        for (i, (&p, &m)) in self.plus.iter().zip(&self.minus).enumerate().skip(1) {
            let z = Complex64::new(real[2 * i - 1], real[2 * i]) / std::f64::consts::SQRT_2;
            full[p as usize] = z;
            full[m as usize] = z.conj();
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use dft_geometry::Cell;
    use dft_units::Hartree;
    #[test]
    fn mirrors_and_orc_four_shapes() {
        for rows in [
            [[8., 0., 0.], [0., 8., 0.], [0., 0., 8.]],
            [[8., 1., 0.], [0., 7., 2.], [1., 0., 9.]],
            [[4., 0., 0.], [0., 5., 0.], [0., 0., 20.]],
            [[7., 0., 0.], [0., 8., 0.], [0., 0., 25.]],
        ] {
            let cell = Cell::from_user_rows_bohr(rows).unwrap();
            let basis = PlaneWaveBasis::new(&cell, Hartree(1.)).unwrap();
            let grid = FftGridSpec::new(&basis, &basis, None).unwrap();
            let half = HalfSphere::new(&basis, &grid).unwrap();
            let n = half.full_len();
            assert_eq!(2 * half.plus().len() - 1, n);
            for (i, &j) in half.mirror().iter().enumerate() {
                assert_eq!(half.mirror()[j as usize] as usize, i);
            }
            assert!(half.plus()[1..].windows(2).all(|w| w[0] < w[1]));
            let x: Vec<_> = (0..n)
                .map(|i| ((i * 17 + 3) as f64).sin() / (n as f64).sqrt())
                .collect();
            let y: Vec<_> = (0..n)
                .map(|i| ((i * 7 + 1) as f64).cos() / (n as f64).sqrt())
                .collect();
            let mut a = vec![Complex64::default(); n];
            let mut b = a.clone();
            half.unpack(&x, &mut a).unwrap();
            half.unpack(&y, &mut b).unwrap();
            let mut back = vec![0.; n];
            half.pack(&a, &mut back).unwrap();
            let dot: f64 = x.iter().zip(&y).map(|(x, y)| x * y).sum();
            let cdot: Complex64 = a.iter().zip(&b).map(|(a, b)| a.conj() * b).sum();
            let norm: f64 = x.iter().map(|x| x * x).sum();
            let cnorm: f64 = a.iter().map(|z| z.norm_sqr()).sum();
            let defect = (dot - cdot.re)
                .abs()
                .max(cdot.im.abs())
                .max((norm - cnorm).abs());
            println!("ORC n={n}: isometry defect={defect:e}");
            // Unit-scale vectors: 1e-15 admits only sqrt(2) and sum rounding.
            assert!(defect < 1e-15);
            assert!(back.iter().zip(&x).all(|(a, b)| (a - b).abs() < 1e-15));
            a[0].im = 2. * ORC_TOLERANCE;
            assert_eq!(
                half.pack(&a, &mut back),
                Err(HalfSphereError::InvalidCoefficients)
            );
            a[0].im = 0.;
            a[half.minus()[1] as usize].re += 1e-10;
            assert!(half.pack(&a, &mut back).is_err());
            a[0].re = f64::NAN;
            assert!(half.pack(&a, &mut back).is_err());
            assert!(half.unpack(&x[..n - 1], &mut a).is_err());
            let mut broken = basis.clone();
            broken.vectors.pop();
            assert_eq!(
                HalfSphere::new(&broken, &grid).unwrap_err(),
                HalfSphereError::InvalidBasis
            );
            let shifted = PlaneWaveBasis::at_k(&cell, Hartree(1.), [0.1, 0., 0.]).unwrap();
            assert_eq!(
                HalfSphere::new(&shifted, &grid).unwrap_err(),
                HalfSphereError::NotGamma
            );
            let mut wrong = grid.clone();
            wrong.g_to_fft[1] = [0; 3];
            assert_eq!(
                HalfSphere::new(&basis, &wrong).unwrap_err(),
                HalfSphereError::InvalidGrid
            );
        }
    }
}
