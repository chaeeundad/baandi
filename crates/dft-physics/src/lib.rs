//! Gamma-only functionals and pseudopotential operators in Hartree atomic units.
pub mod density;
pub mod ewald;
pub mod hamiltonian;
pub mod hartree;
pub mod kinetic;
pub mod xc;

use num_complex::Complex64;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum PhysicsError {
    #[error("invalid density {0}: must be finite and >= -1e-14 Bohr^-3")]
    Density(f64),
    #[error("volume must be finite and positive")]
    Volume,
    #[error("incompatible orbital, occupation, coefficient or grid dimensions")]
    Shape,
    #[error("occupations must be finite and nonnegative; coefficients must be finite")]
    Orbitals,
    #[error("invalid Ewald charges, positions, eta or tail target")]
    EwaldInput,
    #[error("distinct ions occupy the same periodic position")]
    CoincidentIons,
    #[error("numerical result or enumeration bound is not representable")]
    Unrepresentable,
    #[error(transparent)]
    Basis(#[from] dft_basis::BasisError),
    #[error(transparent)]
    Fft(#[from] dft_numerics::FftError),
}

fn volume(volume: f64) -> Result<(), PhysicsError> {
    if volume.is_finite() && volume > 0.0 {
        Ok(())
    } else {
        Err(PhysicsError::Volume)
    }
}

fn orbitals(
    coefficients: &[Vec<Complex64>],
    occupations: &[f64],
    count: usize,
) -> Result<(), PhysicsError> {
    if coefficients.len() != occupations.len() || coefficients.iter().any(|c| c.len() != count) {
        return Err(PhysicsError::Shape);
    }
    if occupations.iter().any(|f| !f.is_finite() || *f < 0.0)
        || coefficients
            .iter()
            .flatten()
            .any(|c| !c.re.is_finite() || !c.im.is_finite())
    {
        return Err(PhysicsError::Orbitals);
    }
    Ok(())
}

#[cfg(test)]
mod tests;
