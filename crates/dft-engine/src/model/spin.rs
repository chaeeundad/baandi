//! Shared spin sectors and finite-temperature contracts.
use super::{ScfError, ScfOptions};

#[derive(Clone, Debug)]
pub struct SpinSector {
    pub n_up: usize,
    pub n_down: usize,
    /// Channel order is up/down along this unit Cartesian axis (no SOC).
    pub axis_cartesian: [f64; 3],
    /// Provenance only; the fixed sector, not this seed, constrains the moment.
    pub initial_moment_seed: Option<f64>,
}
impl SpinSector {
    pub fn validate(&self, electrons: f64) -> Result<(), ScfError> {
        let n = self
            .n_up
            .checked_add(self.n_down)
            .ok_or(ScfError::InvalidInput(
                "spin.constraint: electron count overflow; reduce counts",
            ))?;
        if n == 0 || (n as f64 - electrons).abs() > 1e-10 {
            return Err(ScfError::InvalidInput(
                "spin.constraint: channel counts must sum to valence electrons; correct n_up and n_down",
            ));
        }
        if self.axis_cartesian.iter().any(|v| !v.is_finite())
            || (self.axis_cartesian.iter().map(|v| v * v).sum::<f64>() - 1.0).abs() > 1e-12
        {
            return Err(ScfError::InvalidInput(
                "spin.axis: requires a finite unit Cartesian vector; normalize the axis",
            ));
        }
        if self
            .initial_moment_seed
            .is_some_and(|v| !v.is_finite() || v.abs() > electrons)
        {
            return Err(ScfError::InvalidInput(
                "spin.initial_moment_seed: outside electron capacity; use a finite seed with abs(seed)<=N",
            ));
        }
        Ok(())
    }
    pub fn moment(&self) -> f64 {
        self.n_up as f64 - self.n_down as f64
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ThermodynamicLedger {
    pub internal_energy_ha: f64,
    /// Positive TS; the correction in the free-energy ledger is its negative.
    pub entropy_term_ha: f64,
    pub free_energy_ha: f64,
}
impl ThermodynamicLedger {
    pub fn new(energy: f64, entropy_term: f64) -> Self {
        Self {
            internal_energy_ha: energy,
            entropy_term_ha: entropy_term,
            free_energy_ha: energy - entropy_term,
        }
    }
}
pub(crate) fn validate_options(
    options: &ScfOptions,
    theta_ha: f64,
    band_cap: usize,
) -> Result<(), ScfError> {
    if !theta_ha.is_finite() || theta_ha <= 0.0 {
        return Err(ScfError::InvalidInput(
            "occupations.smearing: require positive finite theta; use fixed for zero temperature",
        ));
    }
    if band_cap == 0 || options.bands.is_some_and(|b| b > band_cap) {
        return Err(ScfError::InvalidInput(
            "occupations.band_cap: require a positive cap at least as large as requested bands",
        ));
    }
    if options.density_tolerance > 1e-7
        || options.energy_tolerance > 1e-8
        || options.eigen_tolerance > 1e-8
    {
        return Err(ScfError::InvalidInput(
            "occupations.certificate: tolerances exceed approved caps; use density<=1e-7 and energy/eigen<=1e-8 Ha",
        ));
    }
    if options.spectral_guard != crate::SpectralGuard::Off {
        return Err(ScfError::InvalidInput(
            "occupations.spectral_guard: finite-T guard unavailable; use off",
        ));
    }
    Ok(())
}
