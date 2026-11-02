//! Slater exchange (alpha=2/3) and unpolarized PZ81 correlation.
use crate::{PhysicsError, volume};
use dft_units::Hartree;
use ndarray::Array3;
use std::f64::consts::PI;

pub const DENSITY_THRESHOLD: f64 = 1e-10;
pub const NEGATIVE_ROUNDOFF: f64 = 1e-14;

pub(crate) fn checked_density(n: f64) -> Result<f64, PhysicsError> {
    if !n.is_finite() || n < -NEGATIVE_ROUNDOFF {
        Err(PhysicsError::Density(n))
    } else {
        Ok(n.max(0.0))
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct XcPoint {
    /// Exchange energy per electron (Ha).
    pub eps_x: f64,
    /// Correlation energy per electron (Ha).
    pub eps_c: f64,
    /// Exchange potential (Ha).
    pub v_x: f64,
    /// Correlation potential (Ha).
    pub v_c: f64,
}
impl XcPoint {
    pub fn eps_xc(self) -> f64 {
        self.eps_x + self.eps_c
    }
    pub fn v_xc(self) -> f64 {
        self.v_x + self.v_c
    }
}

/// Allocation-free point evaluation. At/below threshold all four values vanish.
pub fn lda_pz(n: f64) -> Result<XcPoint, PhysicsError> {
    let n = checked_density(n)?;
    if n <= DENSITY_THRESHOLD {
        return Ok(XcPoint::default());
    }
    let rs = (3.0 / (4.0 * PI)).cbrt() / n.cbrt();
    let eps_x = -0.75 * (3.0 / PI).cbrt() * n.cbrt();
    let (eps_c, derivative) = if rs >= 1.0 {
        let denominator = 1.0 + 1.0529 * rs.sqrt() + 0.3334 * rs;
        (
            -0.1423 / denominator,
            0.1423 * (1.0529 / (2.0 * rs.sqrt()) + 0.3334) / denominator.powi(2),
        )
    } else {
        (
            0.0311 * rs.ln() - 0.048 + 0.0020 * rs * rs.ln() - 0.0116 * rs,
            0.0311 / rs + 0.0020 * (rs.ln() + 1.0) - 0.0116,
        )
    };
    // drs/dn = -rs/(3n); d(n eps_c)/dn = eps_c - rs eps_c'/3.
    Ok(XcPoint {
        eps_x,
        eps_c,
        v_x: 4.0 * eps_x / 3.0,
        v_c: eps_c - rs * derivative / 3.0,
    })
}

pub struct XcGrid {
    pub potential: Array3<f64>,
    pub energy: Hartree,
    /// Integral n(r) v_xc(r) d^3r (Ha), for energy double counting.
    pub potential_integral: Hartree,
}

pub fn evaluate(n: &Array3<f64>, volume_bohr3: f64) -> Result<XcGrid, PhysicsError> {
    volume(volume_bohr3)?;
    if n.is_empty() {
        return Err(PhysicsError::Shape);
    }
    let mut potential = Array3::zeros(n.raw_dim());
    let mut energy = 0.0;
    let mut integral = 0.0;
    for (v, &n) in potential.iter_mut().zip(n) {
        let p = lda_pz(n)?;
        let n = checked_density(n)?;
        *v = p.v_xc();
        energy += n * p.eps_xc();
        integral += n * *v;
    }
    let weight = volume_bohr3 / n.len() as f64;
    energy *= weight;
    integral *= weight;
    if !energy.is_finite() || !integral.is_finite() {
        return Err(PhysicsError::Unrepresentable);
    }
    Ok(XcGrid {
        potential,
        energy: Hartree(energy),
        potential_integral: Hartree(integral),
    })
}
