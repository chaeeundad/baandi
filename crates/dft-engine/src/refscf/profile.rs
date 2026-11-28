//! Resolved once at the CLI boundary; no environment reads in the driver.
use super::{Error, invalid};
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Preconditioning {
    Plain,
    ThomasFermi,
}
#[derive(Clone, Copy, Debug)]
pub struct Profile {
    pub beta: f64,
    pub history: usize,
    /// Total residual Hartree energy expressed in Ry, not per electron.
    pub threshold_ry: f64,
    pub preconditioning: Preconditioning,
}
impl Default for Profile {
    fn default() -> Self {
        Self {
            beta: 0.7,
            history: 8,
            threshold_ry: 1e-8,
            preconditioning: Preconditioning::Plain,
        }
    }
}
impl Profile {
    pub fn validate(self) -> Result<Self, Error> {
        if !self.beta.is_finite() || !(0.0 < self.beta && self.beta <= 1.0) {
            return Err(invalid("beta must be finite and in (0,1]"));
        }
        if !(1..=25).contains(&self.history) {
            return Err(invalid("history must be in 1..=25"));
        }
        if !self.threshold_ry.is_finite() || self.threshold_ry <= 0.0 {
            return Err(invalid("threshold_ry must be finite and positive"));
        }
        Ok(self)
    }
    /// Standard class 5e-9 Ha maps to conv_thr=1e-8 Ry. This is an estimator,
    /// not a proof of delivered total-energy error; verification is separate.
    pub fn from_standard_ha(threshold_ha: f64) -> Result<Self, Error> {
        Self {
            threshold_ry: 2.0 * threshold_ha,
            ..Self::default()
        }
        .validate()
    }
}
pub fn next_ethr_ry(previous: f64, dr2_ry: f64, electrons: f64) -> f64 {
    previous.min(0.1 * dr2_ry / electrons.max(1.0)).max(1e-13)
}

/// Adaptive inner PDE tolerance ratchet for implicit solvation (ADR 0061):
/// clamps `1e-2 * min(sqrt(max(dr2_ry, 0)), previous)` into `[target, max(target, 1e-4)]`.
pub fn next_solvation_tol(previous: f64, dr2_ry: f64, electrons: f64, target: f64) -> f64 {
    let _ = electrons;
    let target = if target.is_finite() && target > 0.0 {
        target
    } else {
        1.0e-10
    };
    let upper = 1.0e-4_f64.max(target);
    if !dr2_ry.is_finite() {
        return upper;
    }
    let prev = if previous.is_finite() && previous > 0.0 {
        previous
    } else {
        upper
    };
    (1.0e-2 * dr2_ry.max(0.0).sqrt().min(prev)).clamp(target, upper)
}
