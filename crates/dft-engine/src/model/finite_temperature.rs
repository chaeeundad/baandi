//! Shared finite-temperature result and occupation data.
use super::*;
use dft_geometry::kpoints::KPoint;
#[cfg(feature = "strict-oracle")]
use dft_physics::occupations::Layout;
use dft_physics::occupations::{Constraint, Evaluation, FermiDirac};

pub struct OccupationWorkspace {
    pub(crate) bands: usize,
    #[cfg(feature = "strict-oracle")]
    pub(crate) eigenvalues: Vec<f64>,
    pub(crate) occupations: Vec<f64>,
    pub(crate) weights: Vec<f64>,
    #[cfg(feature = "strict-oracle")]
    capacities: Vec<f64>,
    policy: FermiDirac,
    constraint: Constraint,
    pub evaluation: Option<Evaluation>,
}
impl OccupationWorkspace {
    pub fn new(
        bands: usize,
        points: &[KPoint],
        collinear: bool,
        theta_ha: f64,
        constraint: Constraint,
    ) -> Result<Self, ScfError> {
        let ns = if collinear { 2 } else { 1 };
        let len = bands
            .checked_mul(points.len())
            .and_then(|n| n.checked_mul(ns))
            .ok_or(ScfError::InvalidInput("occupations.layout: size overflow"))?;
        Ok(Self {
            bands,
            #[cfg(feature = "strict-oracle")]
            eigenvalues: vec![0.0; len],
            occupations: vec![0.0; len],
            weights: points.iter().map(|k| k.weight).collect(),
            #[cfg(feature = "strict-oracle")]
            capacities: if collinear { vec![1.0, 1.0] } else { vec![2.0] },
            policy: FermiDirac::new(theta_ha),
            constraint,
            evaluation: None,
        })
    }
    #[cfg(feature = "strict-oracle")]
    pub(crate) fn evaluate(&mut self, electrons: f64) -> Result<(), ScfError> {
        self.evaluation = None;
        self.evaluation = Some(
            self.policy
                .evaluate(
                    Layout {
                        #[cfg(feature = "strict-oracle")]
                        capacities: &self.capacities,
                        weights: &self.weights,
                        bands: self.bands,
                    },
                    &self.eigenvalues,
                    electrons,
                    self.constraint,
                    &mut self.occupations,
                )
                .map_err(|e| ScfError::InvalidInput(e.0))?,
        );
        Ok(())
    }
    pub(crate) fn block(&self, s: usize, k: usize) -> &[f64] {
        let start = (s * self.weights.len() + k) * self.bands;
        &self.occupations[start..start + self.bands]
    }
    #[cfg(feature = "strict-oracle")]
    pub(crate) fn restore(&mut self, occupations: &[f64], evaluation: Evaluation) {
        self.occupations.copy_from_slice(occupations);
        self.evaluation = Some(evaluation);
    }
    pub fn occupations(&self) -> &[f64] {
        &self.occupations
    }
    pub fn constraint(&self) -> Constraint {
        self.constraint
    }
    pub fn theta_ha(&self) -> f64 {
        self.policy.theta_ha
    }
}
pub use crate::model::ThermodynamicLedger;

#[derive(Clone, Debug)]
pub struct FiniteTCertificate {
    pub passed: bool,
    pub charge_density_residual: f64,
    pub spin_density_residual: Option<f64>,
    pub moment_residual: Option<f64>,
    pub retained_eigen_residual_ha: f64,
    pub occupation_change_electrons: f64,
    pub electron_residual: f64,
    pub free_energy_change_ha: f64,
    pub edge_max_fraction: f64,
    pub edge_electrons: f64,
}
pub struct FiniteTResult {
    pub gamma_real_used: bool,
    pub gamma_real_fallback_reason: Option<String>,
    pub(crate) force_identity: crate::force_identity::ForceIdentity,
    pub converged: bool,
    pub iterations: usize,
    pub state: ScfStep,
    pub occupations: OccupationWorkspace,
    pub thermodynamics: ThermodynamicLedger,
    pub certificate: Option<FiniteTCertificate>,
    /// ONEPATH-3: the standard profile of the run (None: strict reference).
    pub standard: Option<StandardProfile>,
    /// ADR 0043 verification attempts (Some only on the standard path with
    /// the verification on).
    pub verification: Option<VerificationReport>,
    /// Standard path: Coulomb estimate E_H[n_out - n_in] of every step, Ha.
    pub standard_errors_ha: Vec<f64>,
    /// Eigen tolerance of every step's solve, Ha.
    pub requested_eigen_tolerances: Vec<f64>,
}
