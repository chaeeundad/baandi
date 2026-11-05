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

/// Frozen-Miller-set tensile-positive kinetic stress, in Ha/Bohr^3.
/// Coefficients and fractional k are fixed; occupations include k weights
/// supplied by the caller. This is one contribution, not a total stress.
pub fn stress(
    basis: &PlaneWaveBasis,
    coefficients: &[Vec<Complex64>],
    occupations: &[f64],
    volume_bohr3: f64,
) -> Result<[[f64; 3]; 3], PhysicsError> {
    crate::volume(volume_bohr3)?;
    orbitals(coefficients, occupations, basis.vectors().len())?;
    let mut result = [[0.0; 3]; 3];
    // q(eta) = (I+eta)^-T q, hence d(q^2/2)/d eta_ab = -q_a q_b.
    for (c, f) in coefficients.iter().zip(occupations) {
        for (c, q) in c.iter().zip(basis.momenta()) {
            let weight = f * c.norm_sqr() / volume_bohr3;
            for (a, row) in result.iter_mut().enumerate() {
                for (b, value) in row.iter_mut().enumerate() {
                    *value -= weight * q.cartesian[a] * q.cartesian[b];
                }
            }
        }
    }
    if result.iter().flatten().any(|v| !v.is_finite()) {
        return Err(PhysicsError::Unrepresentable);
    }
    Ok(result)
}
