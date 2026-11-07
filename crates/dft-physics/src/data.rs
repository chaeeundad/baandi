//! Types used by the physics APIs, re-exported for orchestration clients.
//! This avoids requiring clients to add backend dependencies merely to construct
//! the arrays, FFT plans, and unit values accepted by the existing physics API.
pub use dft_numerics::{
    Fft3d, FftError, FftPlan3d,
    linalg::{ComplexMatrix, LinalgError},
    pool::WorkerPool,
};
pub use dft_units::Hartree;
pub use ndarray::Array3;
pub use num_complex::Complex64;

/// Gaussian Coulomb screening functions, provided by the adopted libm backend.
/// Arguments are dimensionless. These are also useful for analytic fixtures.
pub use libm::{erf, erfc};

/// Research core quadrature dispatch; the core mesh extent is unchanged.
pub fn core_density_with_rule(
    density_basis: &dft_basis::PlaneWaveBasis,
    grid: &dft_basis::FftGridSpec,
    volume_bohr3: f64,
    species: &[crate::hamiltonian::Species<'_>],
    fft: &mut impl Fft3d,
    rule: dft_upf::radial::RadialEvenRule,
) -> Result<Option<Array3<f64>>, crate::hamiltonian::HamiltonianError> {
    if rule == dft_upf::radial::RadialEvenRule::Simpson3_8 {
        return crate::density::core_density(density_basis, grid, volume_bohr3, species, fft);
    }

    use crate::hamiltonian::{HamiltonianError, structure_factor};
    use dft_upf::radial::spherical_density_transform_with_rule;
    if !species.iter().any(|s| s.potential.core_density().is_some()) {
        return Ok(None);
    }
    if !volume_bohr3.is_finite() || volume_bohr3 <= 0.0 {
        return Err(HamiltonianError::Volume);
    }
    let mut coefficients = vec![Complex64::default(); density_basis.vectors().len()];
    for s in species {
        if let Some(core) = s.potential.core_density() {
            let mut shells = std::collections::BTreeMap::new();
            for (out, g) in coefficients.iter_mut().zip(density_basis.vectors()) {
                let radial = if let Some(&v) = shells.get(&g.norm_squared.to_bits()) {
                    v
                } else {
                    let v = spherical_density_transform_with_rule(
                        s.potential.mesh(),
                        core,
                        g.norm_squared.sqrt(),
                        rule,
                    )?;
                    shells.insert(g.norm_squared.to_bits(), v);
                    v
                };
                *out += radial / volume_bohr3 * structure_factor(g.cartesian, s.positions_bohr)?;
            }
        }
    }
    let mut real = Array3::zeros(grid.dims());
    fft.to_real(&grid.scatter(&coefficients)?, &mut real)?;
    let scale = real.iter().map(|z| z.re.abs()).fold(1.0, f64::max);
    if real
        .iter()
        .any(|z| !z.re.is_finite() || !z.im.is_finite() || z.im.abs() > 1e-12 * scale)
    {
        return Err(HamiltonianError::NonFinite);
    }
    Ok(Some(real.mapv(|z| z.re)))
}
