//! Shared convergence profiles and fixed-point assessment.
use super::{FixedPointCertificate, ScfOptions};
use std::time::Duration;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StandardProfile {
    /// Total Coulomb residual energy, Ha; never normalized by electron count.
    pub threshold_ha: f64,
}
impl Default for StandardProfile {
    fn default() -> Self {
        Self { threshold_ha: 5e-9 }
    }
}
impl StandardProfile {
    pub const CAP_HA: f64 = 5e-3;
    pub const FLOOR_HA: f64 = 5e-14;
    /// Threshold, Ha, that stress needs to meet the 1e-7 Ha/Bohr^3 stress class
    /// (ADR 0043 amendment ONEPATH-2 item 5, provisional until dev08): the CLI
    /// default for inputs that ask for stress, and the loosest threshold the
    /// Python session's `stress()` accepts. A rule by observable, not material.
    pub const STRESS_THRESHOLD_HA: f64 = 5e-11;
    pub fn ratchet(previous: f64, error_ha: f64, electrons: f64) -> f64 {
        previous
            .min(0.1 * error_ha / electrons.max(1.0))
            .clamp(Self::FLOOR_HA, Self::CAP_HA)
    }
    pub fn empty_threshold(tolerance: f64) -> f64 {
        (5.0 * tolerance).max(5e-6)
    }

    /// Residual guard (`BAANDI_STANDARD_RESIDUAL_GUARD`, default off since 2026-11-27): an
    /// EMPIRICAL estimator of the band-energy error left by loose eigensolves,
    /// from the per-band residual norms of the latest eigensolve, in Ha:
    ///
    ///   B = c * sum_k w_k sum_i f_i ||r_ik||^2,   c = 1 Ha^-1.
    ///
    /// Motivation (not a proof): for a normalized Ritz vector x with Ritz value
    /// theta and residual r = (H - theta) x, the Ritz value and the orbital
    /// energy functional are stationary at the exact eigenvector, so their
    /// error is of second order in ||r||; the Kato-Temple inequality bounds it by
    /// ||r||^2 / gap, where gap is the distance to the rest of the spectrum. The
    /// density-based Coulomb estimator does not see this term because it
    /// measures input/output density mismatch only.
    ///
    /// This is NOT a bound (review F03, 2026-11-16). c = 1 Ha^-1 stands for
    /// 1 / gap with the gap scale set to the atomic unit; nothing guarantees
    /// that the actual spectral separation is at least 1 Ha, and an SCF
    /// fixed-point residual is not the total-energy error either. Counter-example
    /// (2 x 2, arithmetic): H = diag(0, 0.01 Ha), x = (sqrt(1 - 1e-6), 1e-3)
    /// gives a Rayleigh-quotient error of 1e-8 Ha while B = ||r||^2 * 1 Ha^-1 =
    /// 9.99999e-11 Ha, below the 5e-9 Ha class threshold. The guard is therefore
    /// an empirical estimator whose adequacy rests on measured delivered errors
    /// against strict references (ADR 0043, "Residual guard is an estimator");
    /// it only delays a stop and no tolerance depends on it being a bound.
    /// Empty bands (f_i = 0) do not contribute.
    pub fn residual_energy_estimate<'a>(
        points: impl IntoIterator<Item = (f64, &'a [f64])>,
        occupations: &[f64],
    ) -> f64 {
        Self::residual_energy_estimate_blocks(
            points
                .into_iter()
                .map(|(weight, residuals)| (weight, residuals, occupations)),
        )
    }

    /// The same estimator over (spin, k) blocks with their own occupations
    /// (ONEPATH-3): sum over blocks of w_k sum_i f_i ||r_i||^2, f in electrons
    /// including the channel capacity (1 per spin channel, 2 closed shell;
    /// fractional under Fermi-Dirac), k weights excluded from f. Exact zeros do
    /// not contribute. Two spin channels with f = 1 reproduce the closed-shell
    /// 2 ||r||^2.
    pub fn residual_energy_estimate_blocks<'a>(
        blocks: impl IntoIterator<Item = (f64, &'a [f64], &'a [f64])>,
    ) -> f64 {
        const C_PER_HA: f64 = 1.0;
        C_PER_HA
            * blocks
                .into_iter()
                .map(|(weight, residuals, occupations)| {
                    weight
                        * residuals
                            .iter()
                            .zip(occupations)
                            .filter(|(_, f)| **f > 0.0)
                            .map(|(r, f)| f * r * r)
                            .sum::<f64>()
                })
                .sum::<f64>()
    }

    /// The guarded stop: the ordinary stop decision AND the residual estimate
    /// below the class threshold. It can only delay a stop, never cause one.
    pub fn residual_guard(stop: bool, bound_ha: f64, threshold_ha: f64) -> bool {
        stop && bound_ha < threshold_ha
    }
}

/// Post-convergence fixed-point verification of the standard profile (ADR 0043).
/// `Off` (default) keeps the standard path exactly as before; the strict
/// profile ignores it (its certificate always runs).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Verification {
    #[default]
    Off,
    /// One more Kohn-Sham map at the delivered density; pass when the Coulomb
    /// fixed-point residual E_H[F(rho_out) - rho_out] is below the threshold.
    FixedPoint,
}
impl Verification {
    pub fn name(self) -> &'static str {
        match self {
            Self::Off => "off",
            Self::FixedPoint => "fixed_point",
        }
    }
}

/// One verification attempt (ADR 0043, V1 to V6). Energies in Ha.
#[derive(Clone, Debug)]
pub struct VerificationAttempt {
    /// SCF iteration whose output density was verified.
    pub iteration: usize,
    pub passed: bool,
    /// E_H[F(rho_out) - rho_out], the pass criterion (V2).
    pub coulomb_error_ha: f64,
    pub threshold_ha: f64,
    /// Eigen tolerance of the verification solve (the next ratchet value).
    pub eigen_tolerance_ha: f64,
    /// Residual estimate (`StandardProfile::residual_energy_estimate`, an
    /// empirical estimator, not a bound; the name is kept for result-file
    /// compatibility) of the verification solve when the residual guard is on.
    pub residual_bound_ha: Option<f64>,
    /// Information only (V3): density L2 change and occupied eigenvalue change,
    /// and whether the strict certificate criteria would also have passed.
    pub density_residual: f64,
    pub occupied_eigenvalue_change_ha: f64,
    pub strict_criteria_met: bool,
    pub davidson_iterations: usize,
    pub operator_applications: usize,
    pub wall: std::time::Duration,
    /// ONEPATH-3, information only (None on the closed-shell path): spin
    /// density L2 change of a collinear run (`density_residual` is then the
    /// charge density change), weighted occupation change sum w |Delta f|
    /// (electrons), largest change of a chemical potential (Ha) and total
    /// moment change (mu_B) of the verification map against the stop step.
    pub spin_density_residual: Option<f64>,
    pub occupation_change_electrons: Option<f64>,
    pub chemical_potential_change_ha: Option<f64>,
    pub moment_change: Option<f64>,
}

/// Diagnostic decisions at eligible verification boundaries, never certificates.
#[derive(Clone, Debug)]
pub struct VerificationScheduleEntry {
    pub iteration: usize,
    pub predicted_additional_maps: usize,
    /// Actual action: trace always attempts; on may wait.
    pub action: &'static str,
    pub would_wait: bool,
    /// Absent on a wait, which evaluates no verification map.
    pub coulomb_error_ha: Option<f64>,
    pub ratio_to_previous_attempt: Option<f64>,
    pub passed: Option<bool>,
}

/// Lazily allocated at the first attempt, outside the ordinary-map path.
#[derive(Clone, Debug)]
pub struct VerificationSchedule {
    pub mode: &'static str,
    pub prediction_evaluations: usize,
    pub records: Vec<VerificationScheduleEntry>,
}

/// All verification attempts of one run; passed when the last one passed.
#[derive(Clone, Debug, Default)]
pub struct VerificationReport {
    pub attempts: Vec<VerificationAttempt>,
    /// None in off mode. Diagnostic history only; never consulted by `passed`.
    pub schedule: Option<VerificationSchedule>,
}
impl VerificationReport {
    pub fn passed(&self) -> bool {
        self.attempts.last().is_some_and(|a| a.passed)
    }
    pub fn wall(&self) -> std::time::Duration {
        self.attempts.iter().map(|a| a.wall).sum()
    }
}

/// Shared ADR 0043 fixed-point assessment, independent of a session or solver.
/// The caller must first evaluate a fresh map at the delivered output density,
/// at the next ratchet tolerance, and supply its measured certificate metrics.
/// This function does not replace that map with an estimator or alter the
/// delivered state. The reference driver uses it without legacy switches.
#[allow(clippy::too_many_arguments)]
pub fn assess_fixed_point(
    options: &ScfOptions,
    iteration: usize,
    eigen_tolerance_ha: f64,
    coulomb_error_ha: f64,
    residual_bound_ha: Option<f64>,
    record: &FixedPointCertificate,
    wall: Duration,
) -> VerificationAttempt {
    let threshold_ha = options
        .standard
        .expect("standard verification")
        .threshold_ha;
    let strict_criteria_met = record.density_residual <= options.density_tolerance
        && record.occupied_eigenvalue_change <= CERT_EIGENVALUE_LIMIT_HA
        && record
            .strict_residuals_per_k
            .iter()
            .flatten()
            .all(|r| *r <= 0.1 * options.eigen_tolerance);
    VerificationAttempt {
        iteration,
        passed: coulomb_error_ha.is_finite()
            && coulomb_error_ha < threshold_ha
            && residual_bound_ha.is_none_or(|b| b.is_finite() && b < threshold_ha),
        coulomb_error_ha,
        threshold_ha,
        eigen_tolerance_ha,
        residual_bound_ha,
        density_residual: record.density_residual,
        occupied_eigenvalue_change_ha: record.occupied_eigenvalue_change,
        strict_criteria_met,
        davidson_iterations: record.davidson_iterations,
        operator_applications: record.operator_applications,
        wall,
        spin_density_residual: None,
        occupation_change_electrons: None,
        chemical_potential_change_ha: None,
        moment_change: None,
    }
}

/// Largest occupied-eigenvalue change (Ha) the fixed-point certificate accepts.
pub(crate) const CERT_EIGENVALUE_LIMIT_HA: f64 = 1e-7;
