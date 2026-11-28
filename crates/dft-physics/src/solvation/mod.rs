//! Track B Gate 2 prototype: MPI-distributed generalized Poisson / linearized
//! Poisson-Boltzmann (LPB) implicit solvation solver in Hartree atomic units.
//!
//! Conventions follow `docs/physics-conventions.md` and
//! `docs/research/extensible-dft/implicit-solvation-theory.md`:
//! - Electron number density `n(r) >= 0` in Bohr^-3 (`rho_el = -n`).
//! - Electrostatic potential `phi(r)` (energy per +1e, Ha); electron potential
//!   energy `U_el(r) = -phi(r)` (Ha).
//! - Discrete SPD operator `A = sum_a D_a^T eps(r) D_a + eps_b * kappa^2 * lambda(r)`
//!   using BAANDI's skew-adjoint Nyquist-zeroed spectral derivative `D_a`.
//! - Polarization-potential split `Delta phi_pol = phi_solv - phi_vac` plus
//!   exact `G = 0` background correction `E_bg,corr = (Q_ion - N_e) * sum_I alpha_z,I^sm`,
//!   reducing identically to zero in vacuum (`eps_b = 1, kappa = 0, tau = 0, beta = 0`).
pub mod cavity;
pub mod poisson;

pub use cavity::{
    build_ionic_pseudo_charge_g, evaluate_cavity_fields, ionic_solvation_forces,
    skew_adjoint_g_cart,
};
pub use poisson::SolvationWorkspace;

use dft_mpi::MpiError;
use dft_numerics::linalg::LinalgError;
use dft_solvers::distributed_fft::DistFftError;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum SolvationError {
    #[error("invalid solvation parameter: {0}")]
    InvalidParameter(&'static str),
    #[error("incompatible solvation grid or slice shape")]
    Shape,
    #[error("non-finite value encountered in solvation evaluation")]
    NonFinite,
    #[error("solvation solver breakdown: {0}")]
    Breakdown(&'static str),
    #[error(
        "solvation solver did not converge in {iterations} iterations (residual {residual_norm:.3e})"
    )]
    NotConverged {
        iterations: usize,
        residual_norm: f64,
    },
    #[error(transparent)]
    DistFft(#[from] DistFftError),
    #[error(transparent)]
    Mpi(#[from] MpiError),
    #[error(transparent)]
    Linalg(#[from] LinalgError),
}

/// Numerical algorithm for the generalized Poisson / LPB PDE.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SolverKind {
    /// Strictly-SPD direct G-space PCG with `np = 2` real-space complex pairing and
    /// mirror-stick unpacking (0 setup FFTs, 4 FFTs/iter on mirror). Exact discrete
    /// variational operator baseline and guaranteed-SPD fallback.
    DirectPcg,
    /// High-speed continuum-equivalent Fisicaro et al. (2016) `phi_tilde = sqrt(eps) * phi`
    /// PCG (2 FFTs/iter; 3.4x-12.1x faster than `DirectPcg` on well-resolved grids).
    FisicaroPcg,
    /// Benchmark arm: polarization-charge fixed-point iteration with Anderson/DIIS (3 FFTs/iter).
    PolarizationDiis,
    /// Gate 2 experimental arm (negative result on sharp cavities): warm-started Fisicaro PCG
    /// followed by `DirectPcg` polish. Retained for benchmark reproducibility; high-G aliasing
    /// between the product-rule and discrete divergence operators prevents Stage 2 from saving
    /// iterations over standalone `DirectPcg`.
    HybridFisicaroDirect,
}

/// Smooth density-dependent cavity switching function `S(n) in [0, 1]`.
/// `S -> 1` in bulk solvent (`n -> 0`), `S -> 0` inside solute (`n` large).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CavityModel {
    /// VASPsol `erfc` switching (Mathew et al. 2014, DOI `10.1063/1.4865107`).
    VaspSol { nc_bohr3: f64, sigma: f64 },
    /// SCCS logarithmic-sine switching (Andreussi et al. 2012, DOI `10.1063/1.3676407`).
    Sccs {
        rho_min_bohr3: f64,
        rho_max_bohr3: f64,
        /// When true, `eps(n) = exp(ln(eps_b) * S(n))`; when false, linear `1 + (eps_b - 1) * S(n)`.
        log_epsilon: bool,
    },
}

impl CavityModel {
    /// Standard VASPsol water cavity parameters (`n_c = 0.0025 Bohr^-3, sigma = 0.6`).
    pub fn vaspsol_water() -> Self {
        Self::VaspSol {
            nc_bohr3: 0.0025,
            sigma: 0.6,
        }
    }

    /// Standard SCCS water cavity parameters (`rho_min = 1e-4, rho_max = 1.5e-3 Bohr^-3`).
    pub fn sccs_water() -> Self {
        Self::Sccs {
            rho_min_bohr3: 1.0e-4,
            rho_max_bohr3: 1.5e-3,
            log_epsilon: true,
        }
    }
}

/// Physical parameters of the implicit solvent and optional LPB electrolyte.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolvationParams {
    /// Bulk relative permittivity `eps_b >= 1.0` (`1.0` = vacuum).
    pub epsilon_bulk: f64,
    /// Inverse bulk Debye screening length `kappa >= 0.0` in Bohr^-1.
    /// The LPB screening operator uses `eps_b * kappa^2 * lambda(r)`.
    pub kappa_inv_bohr: f64,
    /// Effective cavity surface tension `tau >= 0.0` in Ha/Bohr^2.
    pub tau_ha_bohr2: f64,
    /// Optional cavity volume coefficient `beta` in Ha/Bohr^3.
    pub beta_ha_bohr3: f64,
    /// Regularization `delta_S > 0` (Bohr^-1) in `sqrt(|D S|^2 + delta_S^2) - delta_S`.
    pub cavity_grad_reg: f64,
    /// Density floor (Bohr^-3) below which `S = 1, S' = 0` smoothly without `ln(<= 0)`.
    pub density_floor_bohr3: f64,
    /// Cavity switching model.
    pub cavity: CavityModel,
}

impl SolvationParams {
    /// Exact vacuum limit (`eps_b = 1, kappa = 0, tau = 0, beta = 0`).
    pub fn vacuum() -> Self {
        Self {
            epsilon_bulk: 1.0,
            kappa_inv_bohr: 0.0,
            tau_ha_bohr2: 0.0,
            beta_ha_bohr3: 0.0,
            cavity_grad_reg: 1.0e-14,
            density_floor_bohr3: 1.0e-16,
            cavity: CavityModel::vaspsol_water(),
        }
    }

    /// Standard VASPsol aqueous parameters (`eps_b = 78.4`, neutral solvent `kappa = 0`).
    pub fn vaspsol_water() -> Self {
        Self {
            epsilon_bulk: 78.4,
            kappa_inv_bohr: 0.0,
            tau_ha_bohr2: 0.0,
            beta_ha_bohr3: 0.0,
            cavity_grad_reg: 1.0e-14,
            density_floor_bohr3: 1.0e-16,
            cavity: CavityModel::vaspsol_water(),
        }
    }

    /// Standard SCCS aqueous parameters (`eps_b = 78.4`, neutral solvent `kappa = 0`).
    pub fn sccs_water() -> Self {
        Self {
            epsilon_bulk: 78.4,
            kappa_inv_bohr: 0.0,
            tau_ha_bohr2: 0.0,
            beta_ha_bohr3: 0.0,
            cavity_grad_reg: 1.0e-14,
            density_floor_bohr3: 1.0e-16,
            cavity: CavityModel::sccs_water(),
        }
    }

    /// True when every physical solvation coupling is identically off.
    pub fn is_vacuum(&self) -> bool {
        self.epsilon_bulk == 1.0
            && self.kappa_inv_bohr == 0.0
            && self.tau_ha_bohr2 == 0.0
            && self.beta_ha_bohr3 == 0.0
    }

    /// Validate all physical parameters before evaluation.
    pub fn validate(&self) -> Result<(), SolvationError> {
        if !self.epsilon_bulk.is_finite() || self.epsilon_bulk < 1.0 {
            return Err(SolvationError::InvalidParameter(
                "epsilon_bulk must be finite and >= 1.0",
            ));
        }
        if !self.kappa_inv_bohr.is_finite() || self.kappa_inv_bohr < 0.0 {
            return Err(SolvationError::InvalidParameter(
                "kappa_inv_bohr must be finite and >= 0.0",
            ));
        }
        if !self.tau_ha_bohr2.is_finite() || self.tau_ha_bohr2 < 0.0 {
            return Err(SolvationError::InvalidParameter(
                "tau_ha_bohr2 must be finite and >= 0.0",
            ));
        }
        if !self.beta_ha_bohr3.is_finite() {
            return Err(SolvationError::InvalidParameter(
                "beta_ha_bohr3 must be finite",
            ));
        }
        if !self.cavity_grad_reg.is_finite() || self.cavity_grad_reg <= 0.0 {
            return Err(SolvationError::InvalidParameter(
                "cavity_grad_reg must be finite and > 0.0",
            ));
        }
        if !self.density_floor_bohr3.is_finite() || self.density_floor_bohr3 <= 0.0 {
            return Err(SolvationError::InvalidParameter(
                "density_floor_bohr3 must be finite and > 0.0",
            ));
        }
        match self.cavity {
            CavityModel::VaspSol { nc_bohr3, sigma } => {
                if !nc_bohr3.is_finite() || nc_bohr3 <= 0.0 || !sigma.is_finite() || sigma <= 0.0 {
                    return Err(SolvationError::InvalidParameter(
                        "VASPsol nc_bohr3 and sigma must be finite and > 0.0",
                    ));
                }
            }
            CavityModel::Sccs {
                rho_min_bohr3,
                rho_max_bohr3,
                ..
            } => {
                if !rho_min_bohr3.is_finite()
                    || !rho_max_bohr3.is_finite()
                    || rho_min_bohr3 <= 0.0
                    || rho_max_bohr3 <= rho_min_bohr3
                {
                    return Err(SolvationError::InvalidParameter(
                        "SCCS requires 0 < rho_min_bohr3 < rho_max_bohr3",
                    ));
                }
            }
        }
        Ok(())
    }
}

/// Numerical tolerances and solver options.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolverConfig {
    pub kind: SolverKind,
    /// Relative L2 residual tolerance: `||r||_2 <= tolerance * max(1.0, ||b||_2)`.
    pub tolerance: f64,
    pub max_iterations: usize,
    /// History depth for `SolverKind::PolarizationDiis`.
    pub diis_history: usize,
    /// Reuse previous solution in `SolvationWorkspace` with 0-FFT bulk-delta predictor.
    pub warm_start: bool,
}

impl Default for SolverConfig {
    fn default() -> Self {
        Self {
            kind: SolverKind::DirectPcg,
            tolerance: 1.0e-12,
            max_iterations: 200,
            diis_history: 8,
            warm_start: true,
        }
    }
}

impl SolverConfig {
    pub fn validate(&self) -> Result<(), SolvationError> {
        if !self.tolerance.is_finite() || self.tolerance <= 0.0 {
            return Err(SolvationError::InvalidParameter(
                "solver tolerance must be finite and > 0.0",
            ));
        }
        if self.max_iterations == 0 {
            return Err(SolvationError::InvalidParameter(
                "solver max_iterations must be >= 1",
            ));
        }
        if self.kind == SolverKind::PolarizationDiis && self.diis_history == 0 {
            return Err(SolvationError::InvalidParameter(
                "diis_history must be >= 1",
            ));
        }
        Ok(())
    }
}

/// Smooth Gaussian ionic pseudo-charge representing one valence ion core.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IonicPseudoCharge {
    pub charge: f64,
    pub fractional_position: [f64; 3],
    pub width_bohr: f64,
}

/// Diagnostics from one generalized Poisson / LPB solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolveReport {
    pub kind: SolverKind,
    pub iterations: usize,
    pub converged: bool,
    pub residual_norm: f64,
    pub relative_residual: f64,
    /// Directional 3D complex FFTs executed during setup + iterations + post-processing.
    pub fft_count: usize,
    /// True when the rank-local mirror-stick pairing path (`np = 2` on `r -> G`) was used.
    pub mirror_paired: bool,
    /// True when the solve warm-started from a retained previous workspace state.
    pub warm_started: bool,
}

/// Complete variational output of one implicit solvation evaluation on owned planes.
#[derive(Clone, Debug)]
pub struct SolvationEvaluation {
    /// Dual electrostatic solvation energy `Delta E_es = 0.5 * <rho_sol, Delta phi_pol>` (Ha).
    pub electrostatic_ha: f64,
    /// Exact `G = 0` background-uncompensation correction `(Q_ion - N_e) * sum_I alpha_z,I^sm` (Ha).
    pub background_correction_ha: f64,
    /// Regularized non-electrostatic cavitation (+ volume) energy `E_cav` (Ha).
    pub cavitation_ha: f64,
    /// Total solvation energy `Delta E_es + E_bg,corr + E_cav` (Ha).
    pub total_ha: f64,
    /// Electrostatic reaction potential `V_es^solv = -Delta phi_pol - delta_{kappa>0} sum alpha_z^sm` (Ha).
    pub v_es_planes: Vec<f64>,
    /// Dielectric cavity potential `V_diel = -1/(8*pi) * (d eps/dn) * |D phi_solv|^2` (Ha).
    pub v_diel_planes: Vec<f64>,
    /// LPB ionic-screening potential `V_scr = -eps_b*kappa^2/(8*pi) * (d lambda/dn) * phi_solv^2` (Ha).
    pub v_scr_planes: Vec<f64>,
    /// Cavitation + volume potential `V_cav` (Ha).
    pub v_cav_planes: Vec<f64>,
    /// Total solvation Kohn-Sham potential `V_solv = V_es^solv + V_diel + V_scr + V_cav` (Ha).
    pub v_total_planes: Vec<f64>,
    /// Polarization electrostatic potential `Delta phi_pol = phi_solv - phi_vac` on owned planes (Ha/e).
    pub phi_pol_planes: Vec<f64>,
    /// Total electrostatic potential `phi_solv = phi_vac + Delta phi_pol` on owned planes (Ha/e).
    pub phi_solv_planes: Vec<f64>,
    /// Explicit Gaussian ionic pseudo-charge solvation forces `F_I^solv,ion` (Ha/Bohr).
    pub ionic_forces_ha_bohr: Vec<[f64; 3]>,
    /// Inner solver telemetry.
    pub report: SolveReport,
}
