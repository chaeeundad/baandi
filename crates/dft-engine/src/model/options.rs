//! Shared SCF options and occupation contracts.
use super::{StandardProfile, Verification};
use crate::Functional;
pub use crate::{StartDensity, StartWavefunctions};
use dft_basis::PlaneWaveBasis;
use dft_geometry::Cell;
use dft_physics::data::Hartree;

#[derive(Clone, Debug)]
pub enum Mixing {
    Linear {
        beta: f64,
    },
    Pulay {
        beta: f64,
        history: usize,
    },
    /// Optional filtering; q0 is a positive inverse screening length in Bohr^-1.
    KerkerPulay {
        beta: f64,
        history: usize,
        q0_inverse_bohr: f64,
    },
}
impl Default for Mixing {
    fn default() -> Self {
        Self::Pulay {
            beta: 0.3,
            history: 6,
        }
    }
}
/// Provenance label only; the CLI resolves the existing numerical options.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum CompatMode {
    #[default]
    None,
    Reference,
}
impl CompatMode {
    pub fn name(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::Reference => "reference",
        }
    }
}

#[derive(Clone, Debug)]
pub struct ScfOptions {
    pub start_wavefunctions: StartWavefunctions,
    /// Starting density (default uniform, bit-identical to the original start).
    pub start_density: StartDensity,
    pub compat_mode: CompatMode,
    /// Research density projector uses all crystal operations; requires enable_symmetry.
    pub crystal_density_group: bool,
    pub gamma_real: crate::GammaReal,
    pub standard: Option<StandardProfile>,
    /// Experimental standard-profile stagnation guard; never enabled by default.
    pub standard_stall_guard: bool,
    /// Experimental density-residual forcing cap; default off, standard only.
    pub standard_density_forcing: bool,
    /// ADR 0043: fixed-point verification after a standard stop (default
    /// `FixedPoint` since ONEPATH-4, `Off` in `strict_reference()`; requires
    /// the standard profile).
    pub verification: Verification,
    pub hartree_mixing: bool,
    pub spectral_guard: crate::SpectralGuard,
    pub locking: dft_solvers::locking::Locking,
    /// Experimental opt-in hybrid; Davidson remains the default.
    pub eigensolver: dft_solvers::SolverChoice,
    /// Independent k-point workers; zero means one, capped at the k-point count.
    pub threads: usize,
    pub functional: Functional,
    pub pw92_prefactor: crate::Pw92Prefactor,
    pub local_mesh: crate::LocalMesh,
    pub radial_even_rule: crate::RadialEvenRule,
    pub allow_functional_override: bool,
    pub cutoff_hartree: f64,
    pub bands: Option<usize>,
    pub extra_bands: usize,
    pub grid_override: Option<[usize; 3]>,
    /// Opt-in frozen Miller protocol. The reference cell selects all reciprocal
    /// lists; an explicit reference FFT grid is required. None is legacy.
    pub strain_reference: Option<Box<Cell>>,
    pub mixing: Mixing,
    pub density_tolerance: f64,
    pub energy_tolerance: f64,
    pub eigen_tolerance: f64,
    /// Legacy schedule (residual_eigen_schedule=false):
    /// tau_k = max(0.1*eigen_tolerance, initial_eigen_tolerance*decay^(k-1)).
    pub initial_eigen_tolerance: f64,
    pub eigen_tolerance_decay: f64,
    /// Use density-residual scheduling; false selects the legacy geometric schedule.
    pub residual_eigen_schedule: bool,
    /// Start cold solves from deterministic lowest-kinetic plane waves.
    pub free_electron_initial_guess: bool,
    pub max_iterations: usize,
    pub davidson_max_iterations: usize,
    /// Minimum working capacity, combined with the configured per-band floor.
    pub davidson_subspace: usize,
    /// Research capacity floor per band, inclusive range 2..=12; default 6.
    pub davidson_capacity_per_band: usize,
    /// Experimental kinetic Davidson damping; default preserves legacy arithmetic.
    pub preconditioner: dft_solvers::DavidsonPreconditioner,
    pub seed: u64,
    pub divergence_factor: f64,
    pub divergence_patience: usize,
}
impl ScfOptions {
    pub(crate) fn basis(
        &self,
        cell: &Cell,
        k: [f64; 3],
        scale: f64,
    ) -> Result<PlaneWaveBasis, dft_basis::BasisError> {
        let basis = PlaneWaveBasis::at_k(
            self.strain_reference.as_deref().unwrap_or(cell),
            Hartree(scale * self.cutoff_hartree),
            k,
        )?;
        if let Some(reference) = &self.strain_reference {
            basis.in_cell_with_frozen_density(cell, reference)
        } else {
            Ok(basis)
        }
    }
}
/// ADR 0043 amendment ONEPATH-4 (user decision 2026-11-16): the engine default
/// is the single standard SCF path with the fixed-point verification after
/// convergence (`StandardProfile::default()`, 5e-9 Ha; Pulay beta 0.5, the
/// standard-class mixing default of the CLI). The strict path is a reference
/// tool reached only through `ScfOptions::strict_reference()` (mirrored by the
/// CLI and Python switch `BAANDI_REFERENCE_STRICT=on`, see
/// `ScfOptions::from_reference_switch`). Every other field is shared.
impl Default for ScfOptions {
    fn default() -> Self {
        Self {
            standard: Some(StandardProfile::default()),
            verification: Verification::FixedPoint,
            mixing: Mixing::Pulay {
                beta: 0.5,
                history: 6,
            },
            ..Self::strict_reference()
        }
    }
}
impl ScfOptions {
    /// The strict reference path (the engine default before ONEPATH-4, bit for
    /// bit): residual stop rule and the strict fixed-point certificate, Pulay
    /// beta 0.3. For accuracy gates, references and tests that need the strict
    /// arithmetic; production callers use `ScfOptions::default()`.
    pub fn strict_reference() -> Self {
        Self {
            start_wavefunctions: StartWavefunctions::default(),
            start_density: StartDensity::Uniform,
            compat_mode: CompatMode::None,
            crystal_density_group: false,
            gamma_real: crate::GammaReal::default(),
            standard: None,
            standard_stall_guard: false,
            standard_density_forcing: false,
            verification: Verification::Off,
            hartree_mixing: false,
            spectral_guard: crate::SpectralGuard::Off,
            locking: dft_solvers::locking::Locking::Off,
            eigensolver: dft_solvers::SolverChoice::default(),
            threads: 1,
            functional: Functional::LdaPz,
            pw92_prefactor: crate::Pw92Prefactor::Modified,
            local_mesh: crate::LocalMesh::Baandi,
            radial_even_rule: crate::RadialEvenRule::Simpson3_8,
            allow_functional_override: false,
            cutoff_hartree: 5.0,
            bands: None,
            extra_bands: 4,
            grid_override: None,
            strain_reference: None,
            mixing: Mixing::default(),
            density_tolerance: 1e-7,
            energy_tolerance: 1e-8,
            eigen_tolerance: 1e-8,
            initial_eigen_tolerance: 1e-5,
            eigen_tolerance_decay: 0.5,
            residual_eigen_schedule: true,
            free_electron_initial_guess: true,
            max_iterations: 200,
            davidson_max_iterations: 400,
            davidson_subspace: 48,
            davidson_capacity_per_band: 6,
            preconditioner: dft_solvers::DavidsonPreconditioner::Default,
            seed: 2026,
            divergence_factor: 20.0,
            divergence_patience: 3,
        }
    }
    /// Host default that mirrors the reference switch of the CLI:
    /// `BAANDI_REFERENCE_STRICT=on` gives `strict_reference()`, otherwise
    /// `default()` (the standard path with verification). Hosts without an
    /// input file (Python session) use this; engine tests and the numerical
    /// core never read the environment through `default()`.
    pub fn from_reference_switch() -> Self {
        if dft_numerics::switches::on(dft_numerics::switches::Id::REFERENCE_STRICT) {
            Self::strict_reference()
        } else {
            Self::default()
        }
    }
    /// True when these options select the strict reference path.
    pub fn is_strict_reference(&self) -> bool {
        self.standard.is_none()
    }
}
/// Whether `bands` fit a basis of `plane_waves` functions. The solvers raise the
/// Davidson capacity to at least `capacity_per_band * bands` (every solve path), so the check uses
/// that effective capacity rather than the raw `davidson_subspace` option: a raw
/// limit of 48 must not reject a 68-band calculation the solver can run.
pub fn bands_fit_subspace(
    davidson_subspace: usize,
    capacity_per_band: usize,
    bands: usize,
    plane_waves: usize,
) -> bool {
    if !(2..=12).contains(&capacity_per_band) {
        return false;
    }
    let capacity = davidson_subspace
        .max(capacity_per_band.saturating_mul(bands))
        .min(plane_waves);
    bands <= plane_waves && (capacity > bands || bands == plane_waves)
}

/// Explicit occupation contract for the k-point API. Partial occupations and
/// smearing are unsupported. At least one empty band per k is required.
#[derive(Clone, Copy, Debug)]
pub enum OccupationPolicy {
    /// Require a sampled indirect gap larger than this positive threshold (Ha)
    /// and twice the requested eigen residual tolerance.
    Insulator { minimum_gap_hartree: f64 },
    /// Explicitly rejected: this engine cannot compute metallic occupations.
    PartialOccupations,
}

/// On fails closed; Auto records why a complex solve was selected.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum GammaReal {
    Off,
    /// Default (ADR 0040): the real path when the problem qualifies, otherwise
    /// the complex path with the reason recorded.
    #[default]
    Auto,
    On,
}
impl GammaReal {
    pub fn name(self) -> &'static str {
        match self {
            Self::Off => "off",
            Self::Auto => "auto",
            Self::On => "on",
        }
    }
}
