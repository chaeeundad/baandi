//! Real number density, in electrons/Bohr^3.
use crate::{PhysicsError, orbitals, volume};
use dft_basis::FftGridSpec;
use dft_numerics::Fft3d;
// Re-export through the physics adapter for orchestration without a new dependency.
pub use dft_numerics::fused::{FusedOperator, FusedPlan};
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
        if f == 0.0 {
            continue;
        }
        fft.to_real(&grid.scatter(c)?, &mut psi)?;
        density.zip_mut_with(&psi, |n, z| *n += f * z.norm_sqr() / volume_bohr3);
    }
    if density.iter().any(|n: &f64| !n.is_finite()) {
        return Err(PhysicsError::Unrepresentable);
    }
    Ok(density)
}

/// Sphere-based density adapter. The caller constructs the plan with this
/// grid's ordered `g_to_fft()` indices and retains it across calls. The original
/// `from_orbitals` remains available for arbitrary `Fft3d` backends and oracles.
pub fn from_orbitals_fused(
    coefficients: &[Vec<Complex64>],
    occupations: &[f64],
    grid: &FftGridSpec,
    volume_bohr3: f64,
    fft: &mut impl FusedOperator,
) -> Result<Array3<f64>, PhysicsError> {
    volume(volume_bohr3)?;
    orbitals(coefficients, occupations, grid.g_to_fft().len())?;
    if grid.dims() != fft.shape() || grid.g_to_fft().len() != fft.coefficient_count() {
        return Err(PhysicsError::Shape);
    }
    let mut density = Array3::zeros(grid.dims());
    for (c, &weight) in coefficients.iter().zip(occupations) {
        fft.accumulate_density(c, weight, volume_bohr3, density.as_slice_mut().unwrap())?;
    }
    if density.iter().any(|n: &f64| !n.is_finite()) {
        return Err(PhysicsError::Unrepresentable);
    }
    Ok(density)
}
