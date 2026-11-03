//! Diagonal plane-wave kinetic energy.
use crate::{PhysicsError, orbitals};
use dft_basis::PlaneWaveBasis;
use dft_units::Hartree;
use num_complex::Complex64;

/// T = sum_b f_b sum_G |c_b(G)|^2 |k+G|^2/2, in Ha.
pub fn energy(
    basis: &PlaneWaveBasis,
    coefficients: &[Vec<Complex64>],
    occupations: &[f64],
) -> Result<Hartree, PhysicsError> {
    orbitals(coefficients, occupations, basis.vectors().len())?;
    let energy = coefficients
        .iter()
        .zip(occupations)
        .map(|(c, f)| {
            f * c
                .iter()
                .zip(basis.momenta())
                .map(|(c, g)| 0.5 * g.norm_squared * c.norm_sqr())
                .sum::<f64>()
        })
        .sum::<f64>();
    if !energy.is_finite() {
        return Err(PhysicsError::Unrepresentable);
    }
    Ok(Hartree(energy))
}
