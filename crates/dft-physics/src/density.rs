//! Real number density, in electrons/Bohr^3.
use crate::{PhysicsError, orbitals, volume};
use dft_basis::FftGridSpec;
use dft_numerics::Fft3d;
use ndarray::Array3;
use num_complex::Complex64;

/// Caller supplies normalized orbitals in the grid's basis order and a grid
/// large enough for their density sphere. No renormalization is performed.
pub fn from_orbitals(
    coefficients: &[Vec<Complex64>],
    occupations: &[f64],
    grid: &FftGridSpec,
    volume_bohr3: f64,
    fft: &mut impl Fft3d,
) -> Result<Array3<f64>, PhysicsError> {
    volume(volume_bohr3)?;
    orbitals(coefficients, occupations, grid.g_to_fft().len())?;
    if grid.dims() != fft.shape() {
        return Err(PhysicsError::Shape);
    }
    let mut density = Array3::zeros(grid.dims());
    let mut psi = Array3::zeros(grid.dims());
    for (c, &f) in coefficients.iter().zip(occupations) {
        fft.to_real(&grid.scatter(c)?, &mut psi)?;
        density.zip_mut_with(&psi, |n, z| *n += f * z.norm_sqr() / volume_bohr3);
    }
    if density.iter().any(|n: &f64| !n.is_finite()) {
        return Err(PhysicsError::Unrepresentable);
    }
    Ok(density)
}
