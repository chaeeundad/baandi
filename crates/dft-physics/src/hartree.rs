//! Full-grid Coulomb kernel, with the G=0 coefficient explicitly zero.
use crate::{PhysicsError, xc::checked_density};
use dft_geometry::Cell;
use dft_numerics::Fft3d;
use dft_units::Hartree;
use ndarray::Array3;
use num_complex::Complex64;
use std::f64::consts::PI;

pub struct HartreeGrid {
    /// Ha; complex to preserve the canonical Nyquist convention on skew cells.
    pub potential: Array3<Complex64>,
    pub reciprocal_potential: Array3<Complex64>,
    pub energy: Hartree,
}

/// Squared reciprocal lengths in Bohr^-2, including negative even Nyquist.
pub fn squared_frequencies(cell: &Cell, dims: [usize; 3]) -> Result<Array3<f64>, PhysicsError> {
    if dims.contains(&0) {
        return Err(PhysicsError::Shape);
    }
    let count = dims
        .iter()
        .try_fold(1_usize, |count, &n| count.checked_mul(n))
        .ok_or(PhysicsError::Unrepresentable)?;
    if count > isize::MAX as usize / std::mem::size_of::<f64>() {
        return Err(PhysicsError::Unrepresentable);
    }
    let b = cell.reciprocal_columns_inverse_bohr();
    let values = Array3::from_shape_fn(dims, |(i, j, k)| {
        let indices = [i, j, k];
        let m: [f64; 3] = std::array::from_fn(|a| {
            if indices[a] <= (dims[a] - 1) / 2 {
                indices[a] as f64
            } else {
                indices[a] as f64 - dims[a] as f64
            }
        });
        (0..3)
            .map(|a| (0..3).map(|j| b[a][j] * m[j]).sum::<f64>().powi(2))
            .sum::<f64>()
    });
    if values.iter().any(|g| !g.is_finite()) {
        return Err(PhysicsError::Unrepresentable);
    }
    Ok(values)
}

/// Frequencies outside the density sphere are the caller's responsibility.
/// No radial filtering or Nyquist symmetrization is applied.
pub fn evaluate(
    n: &Array3<f64>,
    cell: &Cell,
    fft: &mut impl Fft3d,
) -> Result<HartreeGrid, PhysicsError> {
    if n.shape() != fft.shape() {
        return Err(PhysicsError::Shape);
    }
    let mut real = Array3::zeros(n.raw_dim());
    for (z, &n) in real.iter_mut().zip(n) {
        *z = Complex64::new(checked_density(n)?, 0.0);
    }
    let mut reciprocal = Array3::zeros(n.raw_dim());
    fft.to_reciprocal(&real, &mut reciprocal)?;
    let g2 = squared_frequencies(cell, fft.shape())?;
    let mut energy = 0.0;
    for ((index, z), &g2) in reciprocal.indexed_iter_mut().zip(&g2) {
        if index == (0, 0, 0) {
            *z = Complex64::default();
        } else {
            let kernel = 4.0 * PI / g2;
            energy += kernel * z.norm_sqr();
            *z *= kernel;
        }
    }
    energy *= 0.5 * cell.volume_bohr3();
    fft.to_real(&reciprocal, &mut real)?;
    if !energy.is_finite() || real.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
        return Err(PhysicsError::Unrepresentable);
    }
    Ok(HartreeGrid {
        potential: real,
        reciprocal_potential: reciprocal,
        energy: Hartree(energy),
    })
}
