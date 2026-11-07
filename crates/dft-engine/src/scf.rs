#[path = "convergence.rs"]
mod convergence;
use convergence::CoulombMetric;
pub use convergence::StandardProfile;
#[path = "collinear.rs"]
pub mod collinear;
#[path = "finite_temperature.rs"]
pub mod finite_temperature;
#[path = "replay_capture.rs"]
mod replay_capture;
#[path = "session.rs"]
pub mod session;
use crate::Functional;
use crate::energy::{self, AuditInput, EnergyAudit};
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_geometry::{
    Atom, Cell,
    kpoints::{self, KPoint},
};
use dft_physics::{
    data::*,
    ewald::{self, EwaldOptions, Ion},
    hamiltonian::{self, Hamiltonian, NonlocalOperator, Species},
    hartree, xc,
};
#[doc(hidden)]
pub use dft_solvers::work_accounting;
use dft_solvers::{
    BandTolerances, Davidson, DiagonalPreconditioner, EigenSolution,
    kerker::{Kerker, MixerHistory, PreconditionedPulay},
    mixing::{LinearMixer, PulayMixer},
};
use dft_upf::Pseudopotential;
use session::{CompletionMode, ScfSession, StepStatus};
use std::{
    fmt,
    time::{Duration, Instant},
};

#[derive(Debug)]
pub enum ScfError {
    InvalidInput(&'static str),
    UnsupportedFunctional(String),
    OddElectronCount(f64),
    /// The sampled indirect gap does not certify fixed integer occupations.
    NonInsulating {
        gap_hartree: f64,
        required_hartree: f64,
    },
    InsufficientBands {
        bands: usize,
        occupied: usize,
    },
    Numerical(String),
    NonFinite,
    /// Bounded observation storage cannot reserve another complete step.
    TelemetryExhausted,
    NotConverged {
        history: Vec<IterationRecord>,
    },
    Diverged {
        reason: String,
        history: Vec<IterationRecord>,
    },
}
impl fmt::Display for ScfError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Self::UnsupportedFunctional(label) = self {
            write!(
                f,
                "unsupported UPF functional {label:?}: requires SLA PZ or SLA PZ NOGX NOGC"
            )
        } else {
            write!(f, "{self:?}")
        }
    }
}
impl std::error::Error for ScfError {}
macro_rules! numerical_error { ($($t:ty),*) => { $(impl From<$t> for ScfError { fn from(e:$t)->Self { Self::Numerical(e.to_string()) } })* }; }
numerical_error!(dft_basis::BasisError, FftError);
impl From<dft_physics::PhysicsError> for ScfError {
    fn from(e: dft_physics::PhysicsError) -> Self {
        match e {
            dft_physics::PhysicsError::Unrepresentable => Self::NonFinite,
            dft_physics::PhysicsError::Density(n) if !n.is_finite() => Self::NonFinite,
            _ => Self::Numerical(e.to_string()),
        }
    }
}
impl From<hamiltonian::HamiltonianError> for ScfError {
    fn from(e: hamiltonian::HamiltonianError) -> Self {
        match e {
            hamiltonian::HamiltonianError::NonFinite => Self::NonFinite,
            _ => Self::Numerical(e.to_string()),
        }
    }
}
impl From<dft_solvers::SolverError> for ScfError {
    fn from(e: dft_solvers::SolverError) -> Self {
        match e {
            dft_solvers::SolverError::Linalg(LinalgError::NonFinite) => Self::NonFinite,
            _ => Self::Numerical(e.to_string()),
        }
    }
}

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
#[derive(Clone, Debug)]
pub struct ScfOptions {
    /// Research density projector uses all crystal operations; requires enable_symmetry.
    pub crystal_density_group: bool,
    pub gamma_real: crate::GammaReal,
    pub standard: Option<StandardProfile>,
    /// Experimental standard-profile stagnation guard; never enabled by default.
    pub standard_stall_guard: bool,
    /// Experimental density-residual forcing cap; default off, standard only.
    pub standard_density_forcing: bool,
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
impl Default for ScfOptions {
    fn default() -> Self {
        Self {
            crystal_density_group: false,
            gamma_real: crate::GammaReal::Off,
            standard: None,
            standard_stall_guard: false,
            standard_density_forcing: false,
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
            seed: 2026,
            divergence_factor: 20.0,
            divergence_patience: 3,
        }
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
    /// Pulay feasibility damping; one means the full proposal was accepted.
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

/// Eigenpairs at one canonical fractional reciprocal point, in its basis order.
#[derive(Clone, Debug)]
pub struct KPointSolution {
    /// Canonical fractional k and its normalized integration weight.
    pub point: KPoint,
    /// Eigenvalues/residuals in Ha and normalized coefficients in this k basis.
    pub solution: EigenSolution,
}

#[derive(Debug)]
pub struct ScfResult {
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
    pub spectral_completeness: Option<crate::SpectralReport>,
    pub energy: EnergyAudit,
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
/// Species labels resolve explicitly against Atom::species; the PP must declare SLA PZ or SLA PZ NOGX NOGC.
pub struct SpeciesPotential<'a> {
    pub label: &'a str,
    pub potential: &'a Pseudopotential,
}

/// Immutable geometry-verified symmetry preparation shared with input adapters.
/// Private fields prevent callers from supplying an unverified orbit map.
#[derive(Debug)]
pub struct PreparedSymmetry {
    lattice: [[f64; 3]; 3],
    atoms: Vec<Atom>,
    full: Vec<KPoint>,
    tolerance_bohr: f64,
    time_reversal: bool,
    map: dft_geometry::symmetry::KOrbitMap,
    crystal_operations: Vec<dft_geometry::symmetry::SymmetryOperation>,
}
impl PreparedSymmetry {
    pub fn new(
        cell: &Cell,
        atoms: &[Atom],
        full: &[KPoint],
        tolerance_bohr: f64,
        time_reversal: bool,
    ) -> Result<Self, ScfError> {
        let crystal_operations = dft_geometry::symmetry::detect(cell, atoms, tolerance_bohr)
            .map_err(|e| ScfError::Numerical(e.to_string()))?;
        let map =
            dft_geometry::symmetry::reduce_mesh(full, &crystal_operations, time_reversal, false)
                .map_err(|e| ScfError::Numerical(e.to_string()))?;
        Ok(Self {
            lattice: cell.lattice_columns_bohr(),
            atoms: atoms.to_vec(),
            full: full.to_vec(),
            tolerance_bohr,
            time_reversal,
            map,
            crystal_operations,
        })
    }
    /// Public configurations may change during relaxation. Reuse only for the
    /// exact geometry, mesh and detection settings that produced this record.
    pub fn matches(
        &self,
        cell: &Cell,
        atoms: &[Atom],
        full: &[KPoint],
        tolerance_bohr: f64,
        time_reversal: bool,
    ) -> bool {
        self.lattice == cell.lattice_columns_bohr()
            && self.atoms == atoms
            && self.full == full
            && self.tolerance_bohr.to_bits() == tolerance_bohr.to_bits()
            && self.time_reversal == time_reversal
    }
    pub fn map(&self) -> &dft_geometry::symmetry::KOrbitMap {
        &self.map
    }
    pub fn density_operations(
        &self,
        crystal: bool,
    ) -> Vec<dft_geometry::symmetry::SymmetryOperation> {
        if crystal {
            self.crystal_operations.clone()
        } else {
            self.map
                .operations
                .iter()
                .filter(|o| !o.antiunitary)
                .cloned()
                .collect()
        }
    }
}

/// Cached ionic terms, Hamiltonian application plan, and density FFT plan.
/// `step` is exposed for independent post-SCF fixed-point verification.
pub struct ScfProblem {
    force_model: crate::force_identity::ModelDescription,
    atoms: Vec<Atom>,
    symmetry: Option<dft_physics::symmetry::DensitySymmetrizer>,
    coulomb: Option<CoulombMetric>,
    cell: Cell,
    states: Vec<KState>,
    grid: FftGridSpec,
    ionic_g: Array3<Complex64>,
    ionic_real: Array3<f64>,
    core_density: Option<Array3<f64>>,
    xc_density: Array3<f64>,
    density_floor: Option<f64>,
    fft: FftPlan3d,
    // One reusable FFT workspace per active k worker, shared across phases.
    worker_ffts: Vec<FftPlan3d>,
    pool: std::sync::Arc<WorkerPool>,
    ewald: ewald::EwaldEnergy,
    occupations: Vec<f64>,
    electrons: f64,
    options: ScfOptions,
    setup_time: Duration,
    minimum_gap: Option<f64>,
    gamma_real_exclusion: Option<&'static str>,
}
struct KState {
    gamma_real: crate::GammaRealDiagnostics,
    point: KPoint,
    basis: PlaneWaveBasis,
    grid: FftGridSpec,
    hamiltonian: Hamiltonian,
    band_density: crate::band_density::BandDensity,
}
/// Disjoint state ownership needs Send, never Sync for Hamiltonian scratch.
/// Join every worker before inspecting results in ascending k order. Catch each
/// solve panic so later indices in the same group are still evaluated.
fn ordered_solves<T: Send, R: Send>(
    states: &mut [T],
    requested: usize,
    pool: &WorkerPool,
    solve: impl Fn(usize, &mut T) -> Result<R, ScfError> + Sync,
) -> Result<Vec<R>, ScfError> {
    let solve = &solve;
    pool.install(|| ordered_solves_at(states, requested, 0, solve))
}

fn ordered_solves_at<T: Send, R: Send>(
    states: &mut [T],
    requested: usize,
    first: usize,
    solve: &(impl Fn(usize, &mut T) -> Result<R, ScfError> + Sync),
) -> Result<Vec<R>, ScfError> {
    if states.is_empty() {
        return Ok(Vec::new());
    }
    let workers = requested.max(1).min(states.len());
    if workers == 1 {
        // Evaluate every item before propagating the first error, like the
        // scoped path. Default execution creates no operating-system thread.
        let results: Vec<_> = states
            .iter_mut()
            .enumerate()
            .map(|(index, state)| {
                let index = first + index;
                std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| solve(index, state)))
                    .unwrap_or_else(|_| {
                        Err(ScfError::Numerical(format!(
                            "k index {index}: solve panicked"
                        )))
                    })
                    .map_err(|e| ScfError::Numerical(format!("k index {index}: {e}")))
            })
            .collect();
        return results.into_iter().collect();
    }
    let left_workers = workers / 2;
    let middle = states.len() * left_workers / workers;
    let (left, right) = states.split_at_mut(middle);
    let (left, right) = WorkerPool::join(
        || ordered_solves_at(left, left_workers, first, solve),
        || ordered_solves_at(right, workers - left_workers, first + middle, solve),
    );
    // Both halves have completed before the earliest indexed error is selected.
    let mut results = left?;
    results.extend(right?);
    Ok(results)
}

/// Bound live partial results by the workspace count. Fold each completed
/// chunk in ascending global k order before starting the next chunk; neither
/// chunk boundaries nor completion order can change floating-point association.
fn ordered_chunks<T: Send, W: Send, R: Send>(
    states: &mut [T],
    workspaces: &mut [W],
    pool: &WorkerPool,
    work: impl Fn(usize, &mut T, &mut W) -> Result<R, ScfError> + Sync,
    mut fold: impl FnMut(R) + Send,
) -> Result<(), ScfError> {
    assert!(!workspaces.is_empty());
    let width = workspaces.len();
    let work = &work;
    // Keep the bounded chunk loop resident too: one external synchronization
    // for the phase, rather than waking the pool again for every k chunk.
    pool.install(|| {
        for (chunk, states) in states.chunks_mut(width).enumerate() {
            let mut jobs: Vec<_> = states.iter_mut().zip(workspaces.iter_mut()).collect();
            let partials =
                ordered_solves_at(&mut jobs, width, chunk * width, &|i, (state, workspace)| {
                    work(i, state, workspace)
                })?;
            for partial in partials {
                fold(partial);
            }
        }
        Ok(())
    })
}

/// Output-oriented completion, including a nonconverged final evaluated state.
/// Numerical failures still return `ScfError`; this is not a convergence bypass.
pub struct ScfOutcome {
    pub converged: bool,
    pub result: ScfResult,
}

pub struct ScfStep {
    /// Density projection policy carried with finite-temperature results as well.
    pub crystal_density_group: bool,
    /// First-k compatibility view; `k_points` contains every solution.
    pub solution: EigenSolution,
    /// Every k block in input quadrature order.
    pub k_points: Vec<KPointSolution>,
    pub density: Array3<f64>,
    pub energy: EnergyAudit,
    pub timings: PhaseTimings,
}
impl ScfProblem {
    /// Dispatch evidence in input k order for the latest solve.
    pub fn gamma_real_diagnostics(&self) -> impl Iterator<Item = &crate::GammaRealDiagnostics> {
        self.states.iter().map(|s| &s.gamma_real)
    }

    pub fn new(
        cell: &Cell,
        atoms: &[Atom],
        potentials: &[SpeciesPotential<'_>],
        options: ScfOptions,
    ) -> Result<Self, ScfError> {
        Self::construct(
            cell,
            atoms,
            potentials,
            options,
            &[KPoint {
                fractional: [0.0; 3],
                weight: 1.0,
            }],
            None,
            None,
        )
    }
    /// Construct a weighted Bloch SCF problem. Coordinates are fractional
    /// reciprocal; weights must be positive and sum to one within 1e-14.
    /// Integer occupations are accepted only with an explicit insulating policy.
    /// The sampled indirect gap is checked on the converged state and certificate.
    pub fn with_kpoints(
        cell: &Cell,
        atoms: &[Atom],
        potentials: &[SpeciesPotential<'_>],
        options: ScfOptions,
        points: &[KPoint],
        policy: OccupationPolicy,
    ) -> Result<Self, ScfError> {
        let OccupationPolicy::Insulator {
            minimum_gap_hartree,
        } = policy
        else {
            return Err(ScfError::InvalidInput(
                "partial occupations/smearing unsupported",
            ));
        };
        if !minimum_gap_hartree.is_finite() || minimum_gap_hartree <= 0.0 {
            return Err(ScfError::InvalidInput(
                "positive finite insulating gap required",
            ));
        }
        Self::construct(
            cell,
            atoms,
            potentials,
            options,
            points,
            Some(minimum_gap_hartree),
            None,
        )
    }
    #[allow(clippy::too_many_arguments)]
    fn construct(
        cell: &Cell,
        atoms: &[Atom],
        potentials: &[SpeciesPotential<'_>],
        options: ScfOptions,
        points: &[KPoint],
        minimum_gap: Option<f64>,
        channel: Option<usize>,
    ) -> Result<Self, ScfError> {
        if options.pw92_prefactor == crate::Pw92Prefactor::Original
            && (options.functional != Functional::Pbe || channel.is_some())
        {
            return Err(ScfError::InvalidInput(
                "pw92_prefactor=original requires unpolarized PBE",
            ));
        }
        let gamma_real_exclusion = if channel.is_some() {
            Some("spin-polarized solve")
        } else if points.len() != 1 {
            Some("multiple k points or shifted mesh")
        } else {
            None
        };
        let mut gamma_real = crate::GammaRealDiagnostics::default();
        if options.gamma_real != crate::GammaReal::Off
            && let Some(reason) = gamma_real_exclusion
        {
            gamma_real.fallback(options.gamma_real, reason)?;
        }
        if options.standard_density_forcing
            && (options.standard.is_none() || options.standard_stall_guard)
        {
            return Err(ScfError::InvalidInput(
                "standard_density_forcing requires standard and stall guard off",
            ));
        }
        if options.standard_stall_guard && options.standard.is_none() {
            return Err(ScfError::InvalidInput(
                "standard_stall_guard requires the standard profile",
            ));
        }
        if options
            .standard
            .is_some_and(|p| !p.threshold_ha.is_finite() || p.threshold_ha <= 0.0)
        {
            return Err(ScfError::InvalidInput(
                "positive finite standard_threshold_ha required",
            ));
        }
        if channel.is_some() && (options.standard.is_some() || options.hartree_mixing) {
            return Err(ScfError::InvalidInput(
                "standard profile and Hartree metric currently require unpolarized SCF",
            ));
        }
        if options.standard.is_some()
            && !matches!(options.eigensolver, dft_solvers::SolverChoice::Davidson)
        {
            return Err(ScfError::InvalidInput("standard profile requires Davidson"));
        }
        if options.locking == dft_solvers::locking::Locking::Auto
            && (options.standard.is_some()
                || !matches!(options.eigensolver, dft_solvers::SolverChoice::Davidson))
        {
            return Err(ScfError::InvalidInput("locking requires strict Davidson"));
        }
        kpoints::validate(points)
            .map_err(|_| ScfError::InvalidInput("invalid k-point quadrature"))?;
        let start = Instant::now();
        if atoms.is_empty() || potentials.is_empty() {
            return Err(ScfError::InvalidInput("empty system"));
        }
        for value in [
            options.cutoff_hartree,
            options.density_tolerance,
            options.energy_tolerance,
            options.eigen_tolerance,
            options.initial_eigen_tolerance,
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(ScfError::InvalidInput(
                    "positive finite cutoffs/tolerances required",
                ));
            }
        }
        if !(2..=12).contains(&options.davidson_capacity_per_band) {
            return Err(ScfError::InvalidInput(
                "davidson_capacity_per_band must be an integer from 2 to 12 (minimum 2 for a full correction block)",
            ));
        }
        if options.max_iterations == 0
            || options.davidson_max_iterations == 0
            || options.divergence_patience == 0
            || !options.divergence_factor.is_finite()
            || options.divergence_factor <= 1.0
            || !options.eigen_tolerance_decay.is_finite()
            || options.eigen_tolerance_decay <= 0.0
            || options.eigen_tolerance_decay > 1.0
        {
            return Err(ScfError::InvalidInput(
                "iteration/schedule/divergence configuration",
            ));
        }
        let mut positions = vec![Vec::new(); potentials.len()];
        let mut ions = Vec::new();
        for (i, s) in potentials.iter().enumerate() {
            if potentials[..i].iter().any(|p| p.label == s.label) {
                return Err(ScfError::InvalidInput("duplicate species"));
            }
            if !options.functional.matches_upf(s.potential.functional())
                && !options.allow_functional_override
            {
                if options.functional == Functional::Pbe {
                    return Err(ScfError::InvalidInput(
                        "UPF label mismatches PBE; explicit allow_functional_override required",
                    ));
                }
                if options.functional == Functional::LdaPw {
                    return Err(ScfError::InvalidInput(
                        "UPF label mismatches LDA-PW; explicit allow_functional_override required",
                    ));
                }
                return Err(ScfError::UnsupportedFunctional(
                    s.potential.functional().to_string(),
                ));
            }
        }
        for atom in atoms {
            let index = potentials
                .iter()
                .position(|p| p.label == atom.species)
                .ok_or(ScfError::InvalidInput("unknown species"))?;
            if atom.fractional_position.iter().any(|r| !r.is_finite()) {
                return Err(ScfError::InvalidInput("nonfinite atomic position"));
            }
            positions[index].push(
                cell.fractional_to_cartesian_bohr(
                    atom.fractional_position.map(|r| r.rem_euclid(1.0)),
                ),
            );
            ions.push(Ion {
                charge: potentials[index].potential.z_valence(),
                fractional_position: atom.fractional_position,
            });
        }
        let electrons = ions.iter().map(|i| i.charge).sum::<f64>();
        if !electrons.is_finite()
            || electrons <= 0.0
            || (electrons - electrons.round()).abs() > 1e-10
            || electrons > usize::MAX as f64
        {
            return Err(ScfError::InvalidInput(
                "positive integral valence count required",
            ));
        }
        let count = electrons.round() as usize;
        if channel.is_none() && !count.is_multiple_of(2) {
            return Err(ScfError::OddElectronCount(electrons));
        }
        let occupied = channel.unwrap_or(count / 2);
        let bands = options.bands.unwrap_or(
            occupied
                .checked_add(options.extra_bands)
                .ok_or(ScfError::InvalidInput("band count overflow"))?,
        );
        if bands < occupied {
            return Err(ScfError::InsufficientBands { bands, occupied });
        }
        if minimum_gap.is_some() && bands <= occupied {
            return Err(ScfError::InvalidInput(
                "gap certification requires an empty band at every k",
            ));
        }
        if options.strain_reference.is_some() && options.grid_override.is_none() {
            return Err(ScfError::InvalidInput(
                "strain reference requires a frozen FFT grid",
            ));
        }
        let gamma = options.basis(cell, [0.; 3], 1.0)?;
        let density_basis = options.basis(cell, [0.; 3], 4.0)?;
        let grid = FftGridSpec::new(&gamma, &density_basis, options.grid_override)?;
        let dg = FftGridSpec::new(&density_basis, &density_basis, Some(grid.dims()))?;
        let species: Vec<_> = potentials
            .iter()
            .zip(&positions)
            .map(|(p, r)| Species {
                potential: p.potential,
                positions_bohr: r,
            })
            .collect();
        let ionic_g = dg.scatter(&hamiltonian::ionic_local_with_conventions(
            &density_basis,
            cell.volume_bohr3(),
            &species,
            options.local_mesh,
            options.radial_even_rule,
        )?)?;
        let mut fft = FftPlan3d::new(grid.dims())?;
        let core_density = dft_physics::data::core_density_with_rule(
            &density_basis,
            &dg,
            cell.volume_bohr3(),
            &species,
            &mut fft,
            options.radial_even_rule,
        )?;
        let xc_density = Array3::zeros(
            if core_density.is_some() || options.pw92_prefactor == crate::Pw92Prefactor::Original {
                grid.dims()
            } else {
                [0; 3]
            },
        );
        let mut ionic = Array3::zeros(grid.dims());
        fft.to_real(&ionic_g, &mut ionic)?;
        let scale = ionic.iter().map(|z| z.re.abs()).fold(1.0, f64::max);
        if ionic.iter().any(|z| z.im.abs() > 1e-12 * scale) {
            return Err(ScfError::Numerical("ionic imaginary potential".into()));
        }
        let ionic_real = ionic.mapv(|z| z.re);
        let pool =
            WorkerPool::new(options.threads).map_err(|e| ScfError::Numerical(e.to_string()))?;
        // Every k solve shares the persistent pool. A full-pool budget lets
        // idle outer workers steal fixed inner pieces and gives each worker
        // its own Hamiltonian scratch slot (no modulo-slot collisions).
        let inner = options.threads.max(1);
        let mut construction_points = points.to_vec();
        let states = ordered_solves(
            &mut construction_points,
            options.threads,
            &pool,
            |_, point| {
                let basis = options.basis(cell, point.fractional, 1.0)?;
                // The solvers apply the configured per-band capacity floor, so use
                // the same effective capacity instead of the raw option.
                if !bands_fit_subspace(
                    options.davidson_subspace,
                    options.davidson_capacity_per_band,
                    bands,
                    basis.vectors().len(),
                ) {
                    return Err(ScfError::InvalidInput(
                        "bands/subspace incompatible with k basis",
                    ));
                }
                let k_grid = FftGridSpec::new(&basis, &density_basis, Some(grid.dims()))?;
                let hamiltonian = Hamiltonian::new(
                    cell,
                    &basis,
                    ionic_real.clone(),
                    NonlocalOperator::new_with_rule(
                        &basis,
                        cell.volume_bohr3(),
                        &species,
                        options.radial_even_rule,
                    )?,
                )?
                .with_band_pool(inner, pool.clone())?;
                let band_density =
                    crate::band_density::BandDensity::new(&k_grid, bands, inner, pool.clone())?;
                Ok(KState {
                    gamma_real: gamma_real.clone(),
                    point: KPoint {
                        fractional: basis.k_fractional(),
                        weight: point.weight,
                    },
                    basis,
                    grid: k_grid,
                    hamiltonian,
                    band_density,
                })
            },
        )?;
        let worker_ffts = (0..options.threads.max(1).min(states.len()))
            .map(|_| FftPlan3d::new(grid.dims()))
            .collect::<Result<Vec<_>, _>>()?;
        let ewald = ewald::energy(cell, &ions, EwaldOptions::default())?;
        let occupations = (0..bands)
            .map(|b| {
                if b < occupied {
                    if channel.is_some() { 1.0 } else { 2.0 }
                } else {
                    0.0
                }
            })
            .collect();
        let coulomb = options
            .standard
            .map(|_| CoulombMetric::new(cell, options.cutoff_hartree, grid.dims()))
            .transpose()?;
        let force_model = crate::force_identity::ModelDescription::new(
            cell,
            atoms,
            potentials,
            &options,
            grid.dims(),
            states.iter().map(|s| s.point).collect(),
        );
        Ok(Self {
            gamma_real_exclusion,
            force_model,
            coulomb,
            cell: cell.clone(),
            states,
            grid,
            ionic_g,
            ionic_real,
            core_density,
            xc_density,
            density_floor: None,
            atoms: atoms.to_vec(),
            symmetry: None,
            fft,
            worker_ffts,
            pool,
            ewald,
            occupations,
            electrons,
            options,
            setup_time: start.elapsed(),
            minimum_gap,
        })
    }
    /// Enable only a geometry-verified reduction of the supplied full mesh.
    /// No fields, local constraints or Hubbard terms exist in this SCF model.
    pub fn enable_symmetry(
        &mut self,
        full: &[KPoint],
        tolerance_bohr: f64,
        time_reversal: bool,
    ) -> Result<(), ScfError> {
        let prepared =
            PreparedSymmetry::new(&self.cell, &self.atoms, full, tolerance_bohr, time_reversal)?;
        self.enable_prepared_symmetry(&prepared)
    }
    /// Reuse an immutable verified map, rejecting a different engine geometry.
    pub fn enable_prepared_symmetry(
        &mut self,
        prepared: &PreparedSymmetry,
    ) -> Result<(), ScfError> {
        if prepared.lattice != self.cell.lattice_columns_bohr() || prepared.atoms != self.atoms {
            return Err(ScfError::InvalidInput(
                "symmetry: prepared geometry differs from engine",
            ));
        }
        let map = prepared.map();
        if self.states.len() != map.points.len()
            || self
                .states
                .iter()
                .zip(&map.points)
                .any(|(s, p)| s.point != *p)
        {
            return Err(ScfError::InvalidInput(
                "symmetry: engine k points differ from verified irreducible mesh",
            ));
        }
        let ops = prepared.density_operations(self.options.crystal_density_group);
        if !dft_geometry::symmetry::preserves_fft_grid(self.grid.dims(), &ops) {
            return Err(ScfError::InvalidInput(
                "symmetry.fft_grid: rotations/translations must preserve nonlinear XC quadrature; use compatible_fft_grid at setup",
            ));
        }
        let basis = self.options.basis(&self.cell, [0.; 3], 4.0)?;
        self.force_model.symmetry = ops.clone();
        self.symmetry = Some(dft_physics::symmetry::DensitySymmetrizer::new(
            &basis, &self.grid, &ops,
        )?);
        Ok(())
    }

    /// Opt-in XC-only floor; the legacy construction path is unchanged.
    pub fn set_density_floor(&mut self, floor: f64) -> Result<(), ScfError> {
        if !floor.is_finite() || floor <= 0.0 {
            return Err(ScfError::InvalidInput(
                "xc.density_floor: require finite positive electrons/Bohr^3",
            ));
        }
        self.density_floor = Some(floor);
        self.force_model.floor = Some(floor);
        self.xc_density = Array3::zeros(self.grid.dims());
        Ok(())
    }
    pub(crate) fn validate_clipped_charge(clipped: f64, electrons: f64) -> Result<(), String> {
        let limit = 1e-3_f64.max(1e-4 * electrons);
        if !clipped.is_finite() || clipped > limit {
            return Err(format!(
                "xc.density_floor: converged clipped negative charge {clipped} electrons exceeds {limit}; increase cutoff/grid and inspect the dataset; no certified result"
            ));
        }
        Ok(())
    }
    fn negative_charge(&self, n: &Array3<f64>) -> f64 {
        n.iter()
            .enumerate()
            .map(|(i, &v)| {
                let c = self
                    .core_density
                    .as_ref()
                    .map_or(0.0, |c| c.as_slice().unwrap()[i]);
                (-(v + c)).max(0.0)
            })
            .sum::<f64>()
            * self.cell.volume_bohr3()
            / n.len() as f64
    }
    pub fn core_density(&self) -> Option<&Array3<f64>> {
        self.core_density.as_ref()
    }
    pub fn grid_dimensions(&self) -> [usize; 3] {
        self.grid.dims()
    }
    pub fn electron_count(&self) -> f64 {
        self.electrons
    }
    pub fn uniform_density(&self) -> Array3<f64> {
        Array3::from_elem(self.grid.dims(), self.electrons / self.cell.volume_bohr3())
    }
    pub fn step(
        &mut self,
        n: &Array3<f64>,
        initial: Option<&ComplexMatrix>,
        tolerance: f64,
    ) -> Result<ScfStep, ScfError> {
        if self.states.len() != 1 {
            return Err(ScfError::InvalidInput(
                "use step_kpoints for multiple k points",
            ));
        }
        let initial = initial.map(|v| vec![v.clone()]);
        self.step_kpoints(n, initial.as_deref(), tolerance)
    }
    /// One fixed-point map, with warm-start matrices in quadrature order.
    /// Intermediate maps use integer occupations; gap acceptance occurs in run.
    pub fn step_kpoints(
        &mut self,
        n: &Array3<f64>,
        initial: Option<&[ComplexMatrix]>,
        tolerance: f64,
    ) -> Result<ScfStep, ScfError> {
        self.step_kpoints_with_certification(n, initial, tolerance, true)
    }

    fn step_kpoints_with_certification(
        &mut self,
        n: &Array3<f64>,
        initial: Option<&[ComplexMatrix]>,
        tolerance: f64,
        certify: bool,
    ) -> Result<ScfStep, ScfError> {
        self.step_with_occupations(n, initial, tolerance, certify, None)
    }

    fn step_with_occupations(
        &mut self,
        n: &Array3<f64>,
        initial: Option<&[ComplexMatrix]>,
        tolerance: f64,
        certify: bool,
        mut finite: Option<&mut finite_temperature::OccupationWorkspace>,
    ) -> Result<ScfStep, ScfError> {
        if self.options.crystal_density_group && self.symmetry.is_none() {
            return Err(ScfError::InvalidInput(
                "symmetry.density_group=crystal requires enabled auto symmetry",
            ));
        }
        if initial.is_some_and(|v| v.len() != self.states.len()) {
            return Err(ScfError::InvalidInput("one warm start required per k"));
        }
        if n.shape() != self.grid.dims() {
            return Err(ScfError::InvalidInput("input density grid"));
        }
        let start = Instant::now();
        if let Some(symmetry) = &mut self.symmetry
            && !symmetry.check_invariant(n, &mut self.fft)?
        {
            return Err(ScfError::InvalidInput(
                "symmetry.density: input breaks the selected symmetry; use full mesh",
            ));
        }
        let h = hartree::evaluate(n, &self.cell, &mut self.fft)?;
        let x = self.options.evaluate_xc(
            n,
            self.core_density.as_ref(),
            &mut self.xc_density,
            &self.cell,
            &mut self.fft,
            self.density_floor,
        )?;
        let mut effective = self.ionic_real.clone();
        for ((v, h), x) in effective.iter_mut().zip(&h.potential).zip(&x.potential) {
            *v += h.re + x;
        }
        ordered_solves(
            &mut self.states,
            self.options.threads,
            &self.pool,
            |_, state| {
                state.hamiltonian.set_local_potential(effective.clone())?;
                Ok(())
            },
        )?;
        let mut timings = PhaseTimings {
            potential: start.elapsed(),
            ..Default::default()
        };
        let start = Instant::now();
        let solver = Davidson {
            bands: self.occupations.len(),
            tolerance,
            max_iterations: self.options.davidson_max_iterations,
            max_subspace: self.options.davidson_subspace.max(
                self.options
                    .davidson_capacity_per_band
                    .saturating_mul(self.occupations.len()),
            ),
            block_size: self.occupations.len(),
            seed: self.options.seed,
        };
        // Occupied states determine density; the first empty state determines
        // the sampled gap. Higher empty states are only subspace buffers.
        let policy = BandTolerances {
            strict_bands: (self.occupations.iter().filter(|&&f| f > 0.0).count() + 1)
                .min(self.occupations.len()),
            buffer_tolerance: if self.options.standard.is_some() {
                StandardProfile::empty_threshold(tolerance)
            } else {
                tolerance.max(1e-5)
            },
        };
        // A single k point still has an intra-k budget. Keep its complete solve
        // resident on that pool, so each band/GEMM phase does not enter from
        // the external caller and wait on a LockLatch again.
        let k_points = if self.options.threads <= 1 {
            let mut k_points = Vec::new();
            for (i, state) in self.states.iter_mut().enumerate() {
                let preconditioner = DiagonalPreconditioner::new(
                    state.hamiltonian.kinetic_diagonal().to_vec(),
                    0.1,
                )?;
                let guess = if initial.is_none() && self.options.free_electron_initial_guess {
                    Some(free_electron_guess(
                        state.hamiltonian.kinetic_diagonal(),
                        self.occupations.len(),
                    )?)
                } else {
                    None
                };
                let real_solution = crate::gamma_real::solve_real_gamma(
                    &state.hamiltonian,
                    &state.basis,
                    initial.map(|v| &v[i]),
                    state.hamiltonian.kinetic_diagonal(),
                    self.occupations.len(),
                    tolerance,
                    policy,
                    certify,
                    &self.options,
                    self.pool.clone(),
                    self.gamma_real_exclusion
                        .or(finite.as_ref().map(|_| "Fermi-Dirac occupations")),
                    &mut state.gamma_real,
                )?;
                let solution = if let Some(solution) = real_solution {
                    solution
                } else {
                    std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        if self.options.standard.is_some() {
                            solver.solve_ritz_changes(
                                &state.hamiltonian,
                                initial.map(|v| &v[i]).or(guess.as_ref()),
                                Some(&preconditioner),
                                policy,
                            )
                        } else if self.options.locking == dft_solvers::locking::Locking::Auto {
                            solver.solve_locked(
                                &state.hamiltonian,
                                initial.map(|v| &v[i]).or(guess.as_ref()),
                                Some(&preconditioner),
                                policy,
                            )
                        } else {
                            self.options
                                .eigensolver
                                .solve_with_band_tolerances_and_certification(
                                    &solver,
                                    &state.hamiltonian,
                                    initial.map(|v| &v[i]).or(guess.as_ref()),
                                    Some(&preconditioner),
                                    policy,
                                    certify,
                                )
                        }
                    }))
                    .map_err(|_| ScfError::Numerical(format!("k index {i}: solve panicked")))?
                    .map_err(|e| ScfError::Numerical(format!("k index {i}: {e}")))?
                };
                k_points.push(KPointSolution {
                    point: state.point,
                    solution,
                });
            }
            k_points
        } else {
            let bands = self.occupations.len();
            let free_guess = self.options.free_electron_initial_guess;
            ordered_solves(
                &mut self.states,
                self.options.threads,
                &self.pool,
                |i, state| {
                    let preconditioner = DiagonalPreconditioner::new(
                        state.hamiltonian.kinetic_diagonal().to_vec(),
                        0.1,
                    )?;
                    let guess = if initial.is_none() && free_guess {
                        Some(free_electron_guess(
                            state.hamiltonian.kinetic_diagonal(),
                            bands,
                        )?)
                    } else {
                        None
                    };
                    let real_solution = crate::gamma_real::solve_real_gamma(
                        &state.hamiltonian,
                        &state.basis,
                        initial.map(|v| &v[i]),
                        state.hamiltonian.kinetic_diagonal(),
                        bands,
                        tolerance,
                        policy,
                        certify,
                        &self.options,
                        self.pool.clone(),
                        self.gamma_real_exclusion
                            .or(finite.as_ref().map(|_| "Fermi-Dirac occupations")),
                        &mut state.gamma_real,
                    )?;
                    let solution = if let Some(solution) = real_solution {
                        solution
                    } else {
                        if self.options.standard.is_some() {
                            solver.solve_ritz_changes(
                                &state.hamiltonian,
                                initial.map(|v| &v[i]).or(guess.as_ref()),
                                Some(&preconditioner),
                                policy,
                            )
                        } else if self.options.locking == dft_solvers::locking::Locking::Auto {
                            solver.solve_locked(
                                &state.hamiltonian,
                                initial.map(|v| &v[i]).or(guess.as_ref()),
                                Some(&preconditioner),
                                policy,
                            )
                        } else {
                            self.options
                                .eigensolver
                                .solve_with_band_tolerances_and_certification(
                                    &solver,
                                    &state.hamiltonian,
                                    initial.map(|v| &v[i]).or(guess.as_ref()),
                                    Some(&preconditioner),
                                    policy,
                                    certify,
                                )
                        }?
                    };
                    Ok(KPointSolution {
                        point: state.point,
                        solution,
                    })
                },
            )?
        };
        timings.diagonalization = start.elapsed();
        if let Some(w) = finite.as_mut() {
            for (k, point) in k_points.iter().enumerate() {
                w.eigenvalues[k * w.bands..(k + 1) * w.bands]
                    .copy_from_slice(&point.solution.eigenvalues);
            }
            w.evaluate(self.electrons)?;
        }
        let start = Instant::now();
        let occupations = &self.occupations;
        let finite = finite.as_deref();
        let volume = self.cell.volume_bohr3();
        let orbitals = ordered_solves(
            &mut self.states,
            self.options.threads,
            &self.pool,
            |i, state| {
                let vectors = &k_points[i].solution.eigenvectors;
                let orbitals: Vec<_> = (0..vectors.cols())
                    .map(|b| vectors.column(b).to_vec())
                    .collect();
                state.band_density.compute(
                    &orbitals,
                    finite.map_or(occupations.as_slice(), |w| w.block(0, i)),
                    &state.grid,
                    volume,
                )?;
                Ok(orbitals)
            },
        )?;
        let mut output = Array3::zeros(self.grid.dims());
        // Every k workspace owns its density output. Fold in fixed k order after
        // all jobs finish; no allocated partial arrays or worker locks.
        for state in &self.states {
            for (out, n) in output.iter_mut().zip(state.band_density.density()) {
                *out += state.point.weight * n;
            }
        }
        if let Some(symmetry) = &mut self.symmetry {
            symmetry.apply(&mut output, &mut self.fft)?;
        }
        timings.density = start.elapsed();
        let start = Instant::now();
        let output_xc = if self.core_density.is_some()
            || self.density_floor.is_some()
            || self.options.pw92_prefactor == crate::Pw92Prefactor::Original
        {
            Some(self.options.evaluate_xc(
                &output,
                self.core_density.as_ref(),
                &mut self.xc_density,
                &self.cell,
                &mut self.fft,
                self.density_floor,
            )?)
        } else {
            None
        };
        let shared_audit = energy::prepare_shared(
            &AuditInput {
                cell: &self.cell,
                basis: &self.states[0].basis,
                grid: &self.states[0].grid,
                ionic_g: &self.ionic_g,
                nonlocal: self.states[0].hamiltonian.nonlocal(),
                orbitals: &orbitals[0],
                occupations: finite.map_or(self.occupations.as_slice(), |w| w.block(0, 0)),
                eigenvalues: &k_points[0].solution.eigenvalues,
                eigen_residuals: &k_points[0].solution.residuals,
                input_density: n,
                output_density: &output,
                input_hartree: &h,
                input_xc: &x,
                ewald: self.ewald.energy.0,
            },
            &mut self.fft,
            self.options.functional,
            output_xc.as_ref(),
            self.core_density.as_ref(),
        )?;
        let mut audits = Vec::new();
        ordered_chunks(
            &mut self.states,
            &mut self.worker_ffts,
            &self.pool,
            |i, state, _fft| {
                let k = &k_points[i];
                let audit = energy::evaluate_shared(
                    AuditInput {
                        cell: &self.cell,
                        basis: &state.basis,
                        grid: &state.grid,
                        ionic_g: &self.ionic_g,
                        nonlocal: state.hamiltonian.nonlocal(),
                        orbitals: &orbitals[i],
                        occupations: finite.map_or(self.occupations.as_slice(), |w| w.block(0, i)),
                        eigenvalues: &k.solution.eigenvalues,
                        eigen_residuals: &k.solution.residuals,
                        input_density: n,
                        output_density: &output,
                        input_hartree: &h,
                        input_xc: &x,
                        ewald: self.ewald.energy.0,
                    },
                    &shared_audit,
                )?;
                Ok((state.point.weight, audit))
            },
            |audit| audits.push(audit),
        )?;
        let mut energy = energy::weighted_audits(&audits);
        energy.xc_floor = self.density_floor.map(|floor| energy::XcFloorAudit {
            density_floor: floor,
            input_negative_charge: self.negative_charge(n),
            output_negative_charge: self.negative_charge(&output),
        });
        timings.energy = start.elapsed();
        Ok(ScfStep {
            crystal_density_group: self.options.crystal_density_group,
            solution: k_points[0].solution.clone(),
            k_points,
            density: output,
            energy,
            timings,
        })
    }
    /// Original success-only API; exhausted budgets remain errors.
    pub fn run(self) -> Result<ScfResult, ScfError> {
        let outcome = self.run_impl(false)?;
        if outcome.converged {
            Ok(outcome.result)
        } else {
            Err(ScfError::NotConverged {
                history: outcome.result.history,
            })
        }
    }
    /// Retain the last evaluated state on budget exhaustion for output clients.
    /// Only the existing fixed-point certificate can set `converged` to true.
    pub fn run_outcome(self) -> Result<ScfOutcome, ScfError> {
        self.run_impl(true)
    }
    fn run_impl(self, retain_last: bool) -> Result<ScfOutcome, ScfError> {
        let _analysis = dft_solvers::analysis::Span::new("scf");
        let completion = if retain_last {
            CompletionMode::RetainLast
        } else {
            CompletionMode::Strict
        };
        let mut session =
            ScfSession::new(self, completion, crate::events::ExecutionConfig::default())?;
        session.initialize()?;
        while session.step()?.status == StepStatus::Ready {}
        if !retain_last && session.result().is_none() {
            return Err(ScfError::NotConverged {
                history: session.history().to_vec(),
            });
        }
        session
            .into_result()
            .ok_or(ScfError::InvalidInput("result unavailable"))
    }
    fn sampled_gap(&self, points: &[KPointSolution]) -> f64 {
        let occupied = self.occupations.iter().filter(|f| **f > 0.0).count();
        let valence = points
            .iter()
            .map(|k| k.solution.eigenvalues[occupied - 1])
            .fold(f64::NEG_INFINITY, f64::max);
        let conduction = points
            .iter()
            .map(|k| k.solution.eigenvalues[occupied])
            .fold(f64::INFINITY, f64::min);
        conduction - valence
    }
    fn check_gap(&self, points: &[KPointSolution]) -> Result<Option<f64>, ScfError> {
        let Some(minimum) = self.minimum_gap else {
            return Ok(None);
        };
        let gap = self.sampled_gap(points);
        let required = minimum.max(2.0 * self.options.eigen_tolerance);
        if !gap.is_finite() || gap <= required {
            return Err(ScfError::NonInsulating {
                gap_hartree: gap,
                required_hartree: required,
            });
        }
        Ok(Some(gap))
    }
}
pub fn density_distance(a: &Array3<f64>, b: &Array3<f64>, volume: f64, electrons: f64) -> f64 {
    if a.shape() != b.shape()
        || a.is_empty()
        || !volume.is_finite()
        || volume <= 0.0
        || !electrons.is_finite()
        || electrons <= 0.0
    {
        return f64::NAN;
    }
    (volume / a.len() as f64 * a.iter().zip(b).map(|(a, b)| (a - b).powi(2)).sum::<f64>()).sqrt()
        / electrons
}
/// Coordinate vectors are independent; Davidson supplies seeded columns if rank is missing.
/// Equal kinetic energies retain the original basis index order.
fn free_electron_guess(kinetic: &[f64], bands: usize) -> Result<ComplexMatrix, ScfError> {
    let mut order: Vec<_> = (0..kinetic.len()).collect();
    order.sort_by(|&a, &b| kinetic[a].total_cmp(&kinetic[b]).then(a.cmp(&b)));
    let mut guess = ComplexMatrix::zeros(kinetic.len(), bands.min(kinetic.len()))
        .map_err(dft_solvers::SolverError::from)?;
    for (j, &i) in order.iter().take(guess.cols()).enumerate() {
        guess.column_mut(j)[i] = Complex64::new(1.0, 0.0);
    }
    Ok(guess)
}

fn converged(r: &IterationRecord, o: &ScfOptions) -> bool {
    r.density_residual <= o.density_tolerance
        && r.energy_change
            .is_some_and(|de| de.abs() <= o.energy_tolerance)
        && r.max_eigen_residual <= o.eigen_tolerance
        && r.requested_eigen_tolerance <= o.eigen_tolerance
        && r.charge_error <= 1e-8
        && r.input_charge_error <= 1e-8
}

#[cfg(test)]
#[allow(dead_code)]
#[path = "../tests/common/mod.rs"]
mod pool_fixture;

#[cfg(test)]
mod tests {
    use super::*;

    /// Bounded-work shape/projector benchmark, not a converged material reference.
    /// cargo test --release -p dft-engine worker_utilization_synthetic_twenty -- --ignored --nocapture --test-threads=1
    #[test]
    #[ignore]
    fn worker_utilization_synthetic_twenty() {
        let cell =
            Cell::from_user_rows_bohr([[16.0, 0.0, 0.0], [1.5, 12.0, 0.0], [0.5, 0.7, 10.0]])
                .unwrap();
        let atoms: Vec<_> = (0..20)
            .map(|i| {
                Atom::new(
                    "Xx",
                    [
                        (i % 4) as f64 / 4.0 + 0.07,
                        (i / 4) as f64 / 5.0 + 0.03,
                        if i % 2 == 0 { 0.35 } else { 0.65 },
                    ],
                )
            })
            .collect();
        for nonlocal in [false, true] {
            let potential = pool_fixture::potential(2.0, 1.0, nonlocal, 101);
            let mut reference = None;
            for threads in [1, 2, 4, 8, 12] {
                let start = std::time::Instant::now();
                let outcome = ScfProblem::new(
                    &cell,
                    &atoms,
                    &[SpeciesPotential {
                        label: "Xx",
                        potential: &potential,
                    }],
                    ScfOptions {
                        threads,
                        cutoff_hartree: 10.0,
                        bands: Some(24),
                        max_iterations: 6,
                        ..Default::default()
                    },
                )
                .unwrap()
                .run_outcome()
                .unwrap();
                let elapsed = start.elapsed().as_secs_f64();
                let result = outcome.result;
                let bits: Vec<_> = result.density.iter().map(|x| x.to_bits()).collect();
                let values = (result.total_hartree().to_bits(), bits);
                if let Some(expected) = &reference {
                    assert_eq!(&values, expected);
                } else {
                    reference = Some(values);
                }
                println!(
                    "synthetic20 nonlocal={nonlocal} T={threads} wall_s={elapsed:.9} energy_ha={:.17} iterations={} converged={}",
                    result.total_hartree(),
                    result.history.len(),
                    outcome.converged
                );
            }
        }
    }

    #[test]
    fn scf_reuses_workers_across_iterations() {
        let potential = pool_fixture::potential(2.0, 1.0, true, 101);
        let cell = Cell::from_user_rows_bohr([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]).unwrap();
        for threads in [1, 4] {
            let problem = ScfProblem::new(
                &cell,
                &[Atom::new("Xx", [0.13, 0.21, 0.32])],
                &[SpeciesPotential {
                    label: "Xx",
                    potential: &potential,
                }],
                ScfOptions {
                    threads,
                    cutoff_hartree: 2.5,
                    bands: Some(3),
                    ..Default::default()
                },
            )
            .unwrap();
            let pool = problem.pool.clone();
            let expected = if threads == 1 { 0 } else { threads };
            assert_eq!(pool.workers_created(), expected);
            let result = problem.run().unwrap();
            assert!(result.history.len() > 1);
            assert_eq!(pool.workers_created(), expected);
            println!(
                "SCF T={threads}: workers created={expected}, iterations={}",
                result.history.len()
            );
        }
    }

    #[test]
    fn free_electron_order_and_determinism() {
        let kinetic = [3.0, 0.0, 1.0, 1.0, 2.0];
        let a = free_electron_guess(&kinetic, 4).unwrap();
        assert_eq!(a, free_electron_guess(&kinetic, 4).unwrap());
        for (j, row) in [1, 2, 3, 4].into_iter().enumerate() {
            for i in 0..kinetic.len() {
                assert_eq!(
                    a.column(j)[i],
                    Complex64::new(if i == row { 1.0 } else { 0.0 }, 0.0)
                );
            }
        }
    }

    #[test]
    fn convergence_requires_every_certificate_on_same_iteration() {
        let o = ScfOptions::default();
        let energy = EnergyAudit {
            xc_floor: None,
            kinetic: 0.0,
            local: 0.0,
            local_zero_mode: 0.0,
            nonlocal: 0.0,
            hartree: 0.0,
            xc: 0.0,
            nlcc: None,
            ewald: 0.0,
            band_sum: 0.0,
            double_counting: 0.0,
            band_form: 0.0,
            harris_foulkes: 0.0,
            functional_difference_bound: 0.0,
        };
        let good = IterationRecord {
            standard_error_ha: None,
            standard_redo: false,
            iteration: 2,
            energy,
            energy_change: Some(0.0),
            density_residual: 0.0,
            max_eigen_residual: 0.0,
            requested_eigen_tolerance: o.eigen_tolerance,
            davidson_iterations: 1,
            operator_applications: 0,
            davidson_iterations_per_k: vec![1],
            operator_applications_per_k: vec![0],
            charge_error: 0.0,
            input_charge_error: 0.0,
            min_density: 1.0,
            mixing_fraction: 1.0,
            mixer_reset: false,
            density_contraction: None,
            mixing_history_used: None,
            mixing_resets: 0,
            mixing_raw_recovery: false,
            fixed_point_certificate: None,
            timings: PhaseTimings::default(),
        };
        assert!(converged(&good, &o));
        let mut bad = good.clone();
        bad.energy_change = None;
        assert!(!converged(&bad, &o));
        for failure in 0..6 {
            let mut bad = good.clone();
            match failure {
                0 => bad.energy_change = Some(2.0 * o.energy_tolerance),
                1 => bad.density_residual = 2.0 * o.density_tolerance,
                2 => bad.max_eigen_residual = 2.0 * o.eigen_tolerance,
                3 => bad.charge_error = 2e-8,
                4 => bad.requested_eigen_tolerance = 2.0 * o.eigen_tolerance,
                _ => bad.input_charge_error = 2e-8,
            }
            assert!(!converged(&bad, &o));
        }
        let mut bad = good;
        bad.density_residual = f64::NAN;
        assert!(!converged(&bad, &o));
    }
}

#[cfg(test)]
mod threading_tests {
    use super::*;

    #[test]
    fn chunked_reduction_preserves_global_k_order() {
        // Regrouping the first three terms gives 0 instead of 1. Include an
        // incomplete chunk and more workers than k points.
        let values = [1e16_f64, -1e16, 1.0, 2.0, 4.0];
        let expected = values.iter().fold(0.0_f64, |sum, value| sum + value);
        assert_eq!(expected, 7.0);
        for threads in [1, 2, 3, 4, 8] {
            let mut states = values;
            let mut workspaces = vec![(); threads.min(states.len())];
            let mut sum = 0.0_f64;
            let mut indices = Vec::new();
            ordered_chunks(
                &mut states,
                &mut workspaces,
                &WorkerPool::new(threads).unwrap(),
                |i, value, _| Ok((i, *value)),
                |(i, value)| {
                    indices.push(i);
                    sum += value;
                },
            )
            .unwrap();
            assert_eq!(indices, vec![0, 1, 2, 3, 4]);
            assert_eq!(sum.to_bits(), expected.to_bits());
            let error = ordered_chunks(
                &mut states,
                &mut workspaces,
                &WorkerPool::new(threads).unwrap(),
                |i, _, _| {
                    if i == 4 {
                        return Err(ScfError::InvalidInput("injected chunk failure"));
                    }
                    Ok(())
                },
                |_| {},
            )
            .unwrap_err();
            assert!(error.to_string().contains("k index 4"), "{error}");
        }
    }

    #[test]
    fn ordered_errors_and_panics_join_every_solve() {
        for threads in [1, 2, 3, 4, 8] {
            for panic_at_two in [false, true] {
                let mut visited = [false; 5];
                let error = ordered_solves(
                    &mut visited,
                    threads,
                    &WorkerPool::new(threads).unwrap(),
                    |i, visited| {
                        *visited = true;
                        if i == 2 {
                            // Make the lower-index failure finish later than index 3.
                            std::thread::sleep(Duration::from_millis(10));
                            assert!(!panic_at_two, "injected solve panic");
                            return Err(ScfError::Numerical("injected failure two".into()));
                        }
                        if i == 3 {
                            return Err(ScfError::Numerical("injected failure three".into()));
                        }
                        Ok(i)
                    },
                )
                .unwrap_err();
                assert!(error.to_string().contains("k index 2"), "{error}");
                assert!(visited.into_iter().all(|v| v));
            }
            let mut states = [0; 5];
            assert_eq!(
                ordered_solves(
                    &mut states,
                    threads,
                    &WorkerPool::new(threads).unwrap(),
                    |i, state| {
                        *state = i;
                        Ok(i)
                    }
                )
                .unwrap(),
                vec![0, 1, 2, 3, 4]
            );
            assert_eq!(states, [0, 1, 2, 3, 4]);
        }
    }
}

#[path = "spectral_scf.rs"]
mod spectral_scf;

#[cfg(test)]
mod bands_fit_tests {
    use super::bands_fit_subspace;
    #[test]
    fn effective_capacity_is_at_least_six_times_bands() {
        // 68 bands, 11877 plane waves, raw subspace 48: runnable (si32 benchmark).
        assert!(bands_fit_subspace(48, 6, 68, 11877));
        assert!(bands_fit_subspace(48, 6, 20, 2969));
        // More bands than plane waves is still rejected.
        assert!(!bands_fit_subspace(48, 6, 100, 99));
        // bands equal to the basis size stays allowed (full-space solve).
        assert!(bands_fit_subspace(48, 6, 50, 50));
        // One plane wave more than the bands: capacity is the basis size, allowed.
        assert!(bands_fit_subspace(48, 6, 50, 51));
    }
}
