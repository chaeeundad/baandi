//! Shared collinear result contracts.
pub use super::SpinSector;
use super::finite_temperature::OccupationWorkspace;
use super::*;
use crate::energy::EnergyAudit;
use dft_physics::spin::SpinDensity;

#[derive(Clone, Debug)]
pub struct SpinCertificate {
    pub passed: bool,
    pub charge_density_residual: f64,
    pub spin_density_residual: f64,
    pub energy_change_ha: f64,
    pub occupied_eigenvalue_change_ha: f64,
    pub strict_eigen_residual_ha: [f64; 2],
    pub channel_electrons: [f64; 2],
    pub gap_ha: [Option<f64>; 2],
    pub density_tolerance: f64,
    pub energy_tolerance_ha: f64,
    pub eigen_tolerance_ha: f64,
}
#[derive(Clone, Debug)]
pub struct SpinIteration {
    pub iteration: usize,
    pub energy_ha: f64,
    pub charge_density_residual: f64,
    pub spin_density_residual: f64,
    pub mixing_fraction: f64,
    /// Standard profile only (ONEPATH-3): the collinear Coulomb estimate
    /// E_H[d_n] + E_H[d_m] + E_0[d_m] in Ha (`CoulombMetric::spin_error`).
    pub standard_error_ha: Option<f64>,
    /// Eigen tolerance of this step's solve, Ha.
    pub requested_eigen_tolerance: f64,
    /// Standard profile only: the first map was repeated at the ratcheted
    /// tolerance (same rule as the closed-shell standard path).
    pub standard_redo: bool,
}
/// Count-only diagnostics across ordinary, redo and verification maps.
#[derive(Clone, Debug, Default)]
pub struct SpinShareReport {
    pub shared_maps: usize,
    pub fallback_counts: [usize; 9],
    /// Actual Davidson H columns (atomic-start columns are reported separately).
    pub operator_applications: usize,
    /// Logical two-channel iterations, including reused solutions.
    pub davidson_iterations: usize,
}
impl SpinShareReport {
    pub const REASONS: [&'static str; 9] = [
        "disabled",
        "strict_reference",
        "mpi",
        "finite_temperature",
        "nonzero_seed",
        "unequal_sector",
        "unequal_occupations_or_bands",
        "different_potentials",
        "different_starts",
    ];
    pub fn fallbacks(&self) -> impl Iterator<Item = (&'static str, usize)> + '_ {
        Self::REASONS.into_iter().zip(self.fallback_counts)
    }
}

pub struct CollinearResult {
    pub spin_share: SpinShareReport,
    /// ONEPATH-3B: atomic start of both channels (ADR 0042). Span work is
    /// counted once per k point, Rayleigh-Ritz work once per channel and k.
    pub atomic_start: AtomicStartReport,
    /// ONEPATH-3: the standard profile of the run (None: strict reference).
    pub standard: Option<StandardProfile>,
    /// ADR 0043 verification attempts (Some only on the standard path with
    /// the verification on).
    pub verification: Option<VerificationReport>,
    pub gamma_real_used: bool,
    pub gamma_real_fallback_reason: Option<String>,
    pub(crate) force_identity: crate::force_identity::ForceIdentity,
    pub finite_certificate: Option<super::finite_temperature::FiniteTCertificate>,
    pub finite_occupations: Option<OccupationWorkspace>,
    pub thermodynamics: Option<ThermodynamicLedger>,
    pub spectral_completeness: Option<crate::SpectralReport>,
    pub converged: bool,
    pub sector: SpinSector,
    pub density: SpinDensity,
    /// Spin-major, then input k-point order. Occupations exclude k weights.
    pub k_points: [Vec<KPointSolution>; 2],
    /// Fixed-policy band-only view. Empty for finite T: use finite_occupations,
    /// whose explicit spin/k/band layout cannot be represented by this field.
    pub occupations: [Vec<f64>; 2],
    pub energy: EnergyAudit,
    pub certificate: Option<SpinCertificate>,
    pub history: Vec<SpinIteration>,
    pub channel_electrons: [f64; 2],
    pub moment_mu_b: f64,
}
impl CollinearResult {
    pub fn legacy_scalar_result(&self) -> Result<&ScfResult, ScfError> {
        Err(ScfError::InvalidInput(
            "result.spin: scalar legacy view is ambiguous; request spin-resolved data",
        ))
    }
}
