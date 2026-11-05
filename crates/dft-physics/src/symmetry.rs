//! Allocation-free application of a prepared density group projector.
use crate::PhysicsError;
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_geometry::symmetry::{SpinAction, SymmetryOperation};
use dft_numerics::Fft3d;
use ndarray::Array3;
use num_complex::Complex64;
use std::{collections::BTreeMap, f64::consts::TAU};

pub struct DensitySymmetrizer {
    /// For each output G: input R^T G and exp(-2*pi*i*G.tau)/|group|.
    terms: Vec<(usize, usize, Complex64)>,
    real: Array3<Complex64>,
    input: Array3<Complex64>,
    output: Array3<Complex64>,
}
impl DensitySymmetrizer {
    /// The density sphere must be closed under every operation. Grid padding
    /// outside that sphere is zero; no modulo alias of a rotated G is allowed.
    pub fn new(
        basis: &PlaneWaveBasis,
        grid: &FftGridSpec,
        operations: &[SymmetryOperation],
    ) -> Result<Self, PhysicsError> {
        if operations.is_empty()
            || operations
                .iter()
                .any(|o| o.spin_action != SpinAction::Identity || o.antiunitary)
        {
            return Err(PhysicsError::Shape);
        }
        let dims = grid.dims();
        let index =
            |m| FftGridSpec::map_miller(m, dims).map(|p| (p[0] * dims[1] + p[1]) * dims[2] + p[2]);
        let lookup: BTreeMap<_, _> = basis
            .vectors()
            .iter()
            .map(|g| Ok((g.miller, index(g.miller)?)))
            .collect::<Result<_, dft_basis::BasisError>>()?;
        let mut terms = Vec::new();
        for g in basis.vectors() {
            for op in operations {
                let source =
                    std::array::from_fn(|i| (0..3).map(|j| op.rotation[j][i] * g.miller[j]).sum());
                let src = *lookup.get(&source).ok_or(PhysicsError::Shape)?;
                let phase = -TAU
                    * (0..3)
                        .map(|i| g.miller[i] as f64 * op.translation[i])
                        .sum::<f64>();
                terms.push((
                    lookup[&g.miller],
                    src,
                    Complex64::from_polar(1.0 / operations.len() as f64, phase),
                ));
            }
        }
        Ok(Self {
            terms,
            real: Array3::zeros(dims),
            input: Array3::zeros(dims),
            output: Array3::zeros(dims),
        })
    }
    /// Fail closed for externally supplied densities outside the invariant
    /// subspace. Absolute allowance is 1e-12 electrons/Bohr^3 for FFT roundoff.
    pub fn check_invariant(
        &mut self,
        density: &Array3<f64>,
        fft: &mut impl Fft3d,
    ) -> Result<bool, PhysicsError> {
        self.transform(density, fft)?;
        Ok(density
            .iter()
            .zip(&self.real)
            .all(|(n, z)| n.is_finite() && (n - z.re).abs() <= 1e-12 && z.im.abs() <= 1e-12))
    }
    /// Each collinear channel uses the same unitary identity-spin projector.
    pub fn apply(
        &mut self,
        density: &mut Array3<f64>,
        fft: &mut impl Fft3d,
    ) -> Result<(), PhysicsError> {
        self.transform(density, fft)?;
        for (n, z) in density.iter_mut().zip(&self.real) {
            *n = z.re;
        }
        Ok(())
    }
    fn transform(
        &mut self,
        density: &Array3<f64>,
        fft: &mut impl Fft3d,
    ) -> Result<(), PhysicsError> {
        if density.shape() != self.real.shape() {
            return Err(PhysicsError::Shape);
        }
        for (z, &n) in self.real.iter_mut().zip(density.iter()) {
            *z = Complex64::new(n, 0.0);
        }
        fft.to_reciprocal(&self.real, &mut self.input)?;
        self.output.fill(Complex64::default());
        let input = self.input.as_slice().unwrap();
        let output = self.output.as_slice_mut().unwrap();
        for &(dst, src, phase) in &self.terms {
            output[dst] += phase * input[src];
        }
        fft.to_real(&self.output, &mut self.real)?;
        if self
            .real
            .iter()
            .any(|z| !z.re.is_finite() || z.im.abs() > 1e-12 * z.re.abs().max(1.0))
        {
            return Err(PhysicsError::Unrepresentable);
        }
        Ok(())
    }
}
