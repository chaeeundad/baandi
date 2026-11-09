//! Starting density from the superposition of the UPF atomic charges (the
//! reference code's default starting potential), as an opt-in alternative to the
//! uniform density. Unpolarized fixed-occupation problems only.
use super::*;

fn numerical(e: impl std::fmt::Display) -> ScfError {
    ScfError::Numerical(e.to_string())
}

/// Starting density of the SCF. Default `Uniform` is bit-identical to before.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum StartDensity {
    #[default]
    Uniform,
    /// Sum of the atomic charge densities of the pseudopotentials (UPF
    /// `PP_RHOATOM`, 4 pi r^2 rho), negative ringing clipped, charge renormalised.
    Atomic,
}
impl StartDensity {
    pub fn name(self) -> &'static str {
        match self {
            Self::Uniform => "uniform",
            Self::Atomic => "atomic",
        }
    }
}

/// rho(G) = (1/Omega) sum_atoms exp(-i G.tau) int 4 pi r^2 rho_at(r) sinc(G r) dr,
/// on the density sphere, transformed to the grid, clipped at zero and rescaled
/// to the electron count. None when a species carries no atomic charge samples.
pub(super) fn atomic_superposition(
    cell: &Cell,
    density_basis: &PlaneWaveBasis,
    dg: &FftGridSpec,
    species: &[Species<'_>],
    fft: &mut FftPlan3d,
    electrons: f64,
) -> Result<Option<Array3<f64>>, ScfError> {
    let volume = cell.volume_bohr3();
    let vectors = density_basis.vectors();
    let mut coefficients = vec![Complex64::default(); vectors.len()];
    for s in species {
        let rho = &s.potential.raw().rhoatom.values;
        let mesh = s.potential.mesh();
        if rho.len() != mesh.r().len() || rho.iter().all(|v| *v == 0.0) {
            return Ok(None);
        }
        let mut cache = std::collections::HashMap::new();
        let mut radial = Vec::with_capacity(vectors.len());
        for g in vectors {
            let key = g.norm_squared.to_bits();
            let value = if let Some(v) = cache.get(&key) {
                *v
            } else {
                let q = g.norm_squared.sqrt();
                let samples: Vec<f64> = mesh
                    .r()
                    .iter()
                    .zip(rho)
                    .map(|(&r, &f)| {
                        if q == 0.0 {
                            f
                        } else {
                            f * dft_upf::radial::spherical_bessel(0, q * r).unwrap_or(0.0)
                        }
                    })
                    .collect();
                let v = mesh.integrate(&samples).map_err(numerical)? / volume;
                cache.insert(key, v);
                v
            };
            radial.push(value);
        }
        for pos in s.positions_bohr {
            for ((c, g), f) in coefficients.iter_mut().zip(vectors).zip(&radial) {
                let phase =
                    -(g.cartesian[0] * pos[0] + g.cartesian[1] * pos[1] + g.cartesian[2] * pos[2]);
                *c += Complex64::new(phase.cos(), phase.sin()) * f;
            }
        }
    }
    let reciprocal = dg.scatter(&coefficients).map_err(numerical)?;
    let mut real = Array3::zeros(reciprocal.raw_dim());
    fft.to_real(&reciprocal, &mut real)?;
    let mut density = real.mapv(|z: Complex64| z.re.max(0.0));
    let charge = density.sum() * volume / density.len() as f64;
    if !(charge.is_finite() && charge > 0.0) {
        return Ok(None);
    }
    density.mapv_inplace(|v| v * electrons / charge);
    Ok(Some(density))
}
