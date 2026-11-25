//! Shared SCF outputs and fixed-point records.
use super::{CompatMode, ScfError, StandardProfile, VerificationReport};
pub use crate::AtomicStartReport;
use crate::energy::{self, EnergyAudit};
use dft_geometry::kpoints::KPoint;
use dft_physics::{
    data::{Array3, ComplexMatrix},
    ewald,
};
use dft_solvers::EigenSolution;
use std::time::Duration;

#[derive(Clone, Debug, Default)]
pub struct PhaseTimings {
    pub potential: Duration,
    pub diagonalization: Duration,
    pub density: Duration,
    pub energy: Duration,
    pub mixing: Duration,
}
#[derive(Clone, Debug)]
pub struct IterationRecord {
    pub standard_error_ha: Option<f64>,
    pub standard_redo: bool,
    pub iteration: usize,
    pub energy: EnergyAudit,
    pub energy_change: Option<f64>,
    pub density_residual: f64,
    pub max_eigen_residual: f64,
    pub requested_eigen_tolerance: f64,
    pub davidson_iterations: usize,
    /// Hamiltonian vector applications in this solve; certificate work is separate.
    pub operator_applications: usize,
    pub davidson_iterations_per_k: Vec<usize>,
    pub operator_applications_per_k: Vec<usize>,
    pub charge_error: f64,
    pub input_charge_error: f64,
    pub min_density: f64,
    /// Additional Pulay feasibility damping; one means the full proposal was
    /// accepted. This is NOT the configured beta: beta already multiplies the
    /// residual in x + beta R - (Delta x + beta Delta R) gamma. A Pulay step
    /// need not be parallel to R, so it has no single effective scalar beta.
    pub mixing_fraction: f64,
    pub mixer_reset: bool,
    /// Raw SCF residual ratio to the previous iteration; None initially/after zero.
    pub density_contraction: Option<f64>,
    /// Differences used by this proposal, or None if no mixing was needed.
    pub mixing_history_used: Option<usize>,
    /// Cumulative feasibility resets; excludes parameter-validation reset.
    pub mixing_resets: usize,
    /// New strategy only: a raw convex step rescued an infeasible first correction.
    pub mixing_raw_recovery: bool,
    pub fixed_point_certificate: Option<FixedPointCertificate>,
    pub timings: PhaseTimings,
    /// SCF-TRAJECTORY: what the mixer did with this map's pair (coefficients,
    /// Gram conditioning, predicted residual, step length). Diagnostic only;
    /// None when no mixing followed or the mixer reports none. Not carried by
    /// checkpoint v1.
    pub pulay: Option<dft_solvers::mixing::PulayStepDiagnostics>,
}
#[derive(Clone, Debug)]
pub struct FixedPointCertificate {
    /// True residuals (Ha) for the occupied bands and first empty band.
    pub strict_residuals_per_k: Vec<Vec<f64>>,
    /// Informational only: these higher empty bands are not certified.
    pub buffer_residuals_per_k: Vec<Vec<f64>>,
    pub density_residual: f64,
    pub occupied_eigenvalue_change: f64,
    pub davidson_iterations: usize,
    /// Hamiltonian vector applications in this certificate solve.
    pub operator_applications: usize,
    pub davidson_iterations_per_k: Vec<usize>,
    pub operator_applications_per_k: Vec<usize>,
    pub timings: PhaseTimings,
}
/// Eigenpairs at one canonical fractional reciprocal point, in its basis order.
#[derive(Clone, Debug)]
pub struct KPointSolution {
    /// Canonical fractional k and its normalized integration weight.
    pub point: KPoint,
    /// Eigenvalues/residuals in Ha and normalized coefficients in this k basis.
    pub solution: EigenSolution,
    /// MPI-LAZY-ORBITALS (`BAANDI_M2_LAZY=on`, session maps only): the block
    /// held as distributed stick rows; `solution.eigenvectors` is then a
    /// 0 x bands placeholder until `ScfProblem::materialize_orbitals`. Never set
    /// on public `step_kpoints*` maps or with the switch off; set on a result
    /// only after `gamma_real::set_delivered_orbitals(false)` (LAZY-RESULT,
    /// `ScfResult::materialize_orbitals`).
    pub(crate) lazy: Option<crate::refscf::lazy::LazyOrbitals>,
}

/// Warm-start matrices and lazy rows of earlier k points (MPI-LAZY-ORBITALS).
#[doc(hidden)]
pub fn split_warm(
    points: Option<&[KPointSolution]>,
) -> (
    Option<Vec<ComplexMatrix>>,
    Vec<Option<crate::refscf::lazy::LazyOrbitals>>,
) {
    match points {
        None => (None, Vec::new()),
        Some(p) => (
            Some(p.iter().map(|k| k.solution.eigenvectors.clone()).collect()),
            p.iter().map(|k| k.lazy.clone()).collect(),
        ),
    }
}

/// Replace every lazy block of `k_points` by its gathered complex form
/// (collective on the M2 path) and drop the handles.
pub(crate) fn materialize_kpoints(k_points: &mut [KPointSolution]) -> Result<(), ScfError> {
    for k in k_points {
        if let Some(lazy) = k.lazy.take() {
            k.solution.eigenvectors = lazy.into_full()?;
        }
    }
    Ok(())
}

impl ScfResult {
    /// LAZY-RESULT: whether every k block is held as a full complex block
    /// (false only for an M2 result delivered lazily after
    /// `gamma_real::set_delivered_orbitals(false)`).
    pub fn orbitals_delivered(&self) -> bool {
        self.k_points.iter().all(|k| k.lazy.is_none())
    }
    /// Gather lazily delivered orbital blocks (collective on the M2 path:
    /// every rank must call it at the same point); a no-op otherwise.
    pub fn materialize_orbitals(&mut self) -> Result<(), ScfError> {
        if self.orbitals_delivered() {
            return Ok(());
        }
        materialize_kpoints(&mut self.k_points)?;
        if let Some(first) = self.k_points.first() {
            self.orbitals = first.solution.eigenvectors.clone();
        }
        Ok(())
    }
}

#[derive(Debug)]
pub struct ScfResult {
    pub atomic_start: AtomicStartReport,
    pub compat_mode: CompatMode,
    /// Density projection policy used to produce this state; derivatives are unsupported.
    pub crystal_density_group: bool,
    pub standard_stall_guard: bool,
    /// Experimental density-residual forcing cap; default off, standard only.
    pub standard_density_forcing: bool,
    pub pw92_prefactor: crate::Pw92Prefactor,
    pub local_mesh: crate::LocalMesh,
    pub radial_even_rule: crate::RadialEvenRule,
    pub gamma_real_used: bool,
    pub gamma_real_fallback_reason: Option<String>,
    pub(crate) force_identity: crate::force_identity::ForceIdentity,
    pub standard: Option<StandardProfile>,
    pub certified: bool,
    /// ADR 0043: verification attempts (Some only when verification was on).
    pub verification: Option<VerificationReport>,
    pub spectral_completeness: Option<crate::SpectralReport>,
    pub energy: EnergyAudit,
    /// B2 (`BAANDI_ENERGY_IDENTITY=on`): identity against explicit nonlocal
    /// energy at the final state; `energy` always holds the explicit audit.
    pub energy_identity: Option<energy::EnergyIdentityCheck>,
    /// First k point, retained for Gamma API compatibility; see `k_points`.
    pub eigenvalues: Vec<f64>,
    pub occupations: Vec<f64>,
    /// First k point, retained for Gamma API compatibility.
    pub orbitals: ComplexMatrix,
    /// All eigenpairs with normalized weights and canonical fractional k.
    pub k_points: Vec<KPointSolution>,
    /// Sampled indirect gap in Ha; None for the legacy Gamma constructor.
    pub gap_hartree: Option<f64>,
    pub density: Array3<f64>,
    pub charge: f64,
    pub history: Vec<IterationRecord>,
    pub ewald_breakdown: ewald::EwaldEnergy,
    pub setup_time: Duration,
    pub total_time: Duration,
}
impl ScfResult {
    pub fn total_hartree(&self) -> f64 {
        self.energy.total_hartree()
    }
    pub fn total_rydberg(&self) -> f64 {
        self.energy.total_rydberg()
    }
}
/// Output-oriented completion, including a nonconverged final evaluated state.
/// Numerical failures still return `ScfError`; this is not a convergence bypass.
pub struct ScfOutcome {
    pub converged: bool,
    pub result: ScfResult,
}
