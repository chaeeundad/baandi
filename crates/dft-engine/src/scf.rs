#[path = "atomic_start.rs"]
mod atomic_start;
#[path = "start_density.rs"]
mod start_density;
pub use atomic_start::{AtomicStartReport, StartWavefunctions};
pub use dft_solvers::DavidsonPreconditioner;
pub use start_density::StartDensity;
#[path = "convergence.rs"]
mod convergence;
use convergence::CoulombMetric;
pub use convergence::{StandardProfile, Verification, VerificationAttempt, VerificationReport};
#[path = "collinear.rs"]
pub mod collinear;
/// C1 cutoff ladder (coarse-cutoff start, full-cutoff finish); opt-in.
#[path = "cutoff_ladder.rs"]
pub mod cutoff_ladder;
/// Track A pre-converger; a child of `scf` for read access to problem internals.
#[path = "direct_min.rs"]
pub mod direct_min;
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
    /// ADR 0043: fixed-point verification after a standard stop (default off;
    /// requires the standard profile).
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
impl Default for ScfOptions {
    fn default() -> Self {
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
    start_density: Option<Array3<f64>>,
    atomic_species: Option<Vec<atomic_start::AtomicSpecies>>,
    atomic_start: AtomicStartReport,
    #[cfg(test)]
    force_real_rebuild: bool,
    #[cfg(test)]
    force_energy_identity: bool,
    /// Test hook: the next this-many verifications fail regardless of the criterion.
    #[cfg(test)]
    pub(crate) force_verification_failures: usize,
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
    /// POTENTIAL-REUSE: Hartree/XC fields of a density built by the certificate
    /// pre-screen, consumed by the next map only if its input density is
    /// bit-identical (a deterministic function of the same inputs, so reuse is
    /// exact). Always taken (and dropped on mismatch) by the next map.
    reused_fields: Option<ReusedFields>,
    /// Per-worker FFT plans of the certificate pre-screen, created on first use.
    screen_ffts: Vec<FftPlan3d>,
    /// B2: set by the session for ordinary (non-certificate) maps when
    /// `BAANDI_ENERGY_IDENTITY=on`; integer occupations only.
    energy_identity: bool,
}
struct ReusedFields {
    density: Array3<f64>,
    hartree: hartree::HartreeGrid,
    xc: xc::XcGrid,
    effective: Array3<f64>,
}
struct KState {
    real_solve: Option<crate::gamma_real::RealSolveState>,
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

/// Cluster matrix sums sum_r conj(psi_i) dV psi_j of the certificate pre-screen,
/// one serial sum per entry, entries split over the pool (order-independent).
fn cluster_entries(psi: &[&[Complex64]], delta: &[f64], entries: &mut [(usize, usize, Complex64)]) {
    if entries.len() <= 2 {
        for (i, j, out) in entries.iter_mut() {
            let mut sum = Complex64::default();
            for ((a, b), d) in psi[*i].iter().zip(psi[*j]).zip(delta) {
                sum += a.conj() * b * *d;
            }
            *out = sum;
        }
        return;
    }
    let (left, right) = entries.split_at_mut(entries.len() / 2);
    WorkerPool::join(
        || cluster_entries(psi, delta, left),
        || cluster_entries(psi, delta, right),
    );
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
    /// True when `energy.nonlocal` comes from the band-sum identity (B2,
    /// `BAANDI_ENERGY_IDENTITY=on`, non-final iterations only).
    pub nonlocal_from_identity: bool,
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
        if options.preconditioner != dft_solvers::DavidsonPreconditioner::Default
            && (channel.is_some()
                || !matches!(options.eigensolver, dft_solvers::SolverChoice::Davidson))
        {
            return Err(ScfError::InvalidInput(
                "nondefault preconditioner requires unpolarized Davidson",
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
        if options.verification != Verification::Off && options.standard.is_none() {
            return Err(ScfError::InvalidInput(
                "verification requires the standard profile; the strict profile always runs its certificate",
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
                "standard profile and Hartree metric currently require closed-shell fixed occupations (not collinear spin or Fermi-Dirac smearing runs); use the strict profile with smearing",
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
        let start_density = if options.start_density == StartDensity::Atomic {
            start_density::atomic_superposition(
                cell,
                &density_basis,
                &dg,
                &species,
                &mut fft,
                electrons,
            )?
        } else {
            None
        };
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
                    real_solve: None,
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
        let mut atomic_start = AtomicStartReport::default();
        let atomic_species = if options.start_wavefunctions == StartWavefunctions::Atomic {
            let reason = if channel.is_some() {
                Some("spin-polarized occupations")
            } else if species.iter().any(|s| {
                s.potential
                    .atomic_wavefunctions()
                    .is_none_or(|w| w.is_empty())
            }) {
                Some("missing PP_PSWFC atomic orbitals")
            } else if species.iter().any(|s| {
                s.potential
                    .atomic_wavefunctions()
                    .unwrap()
                    .iter()
                    .any(|w| w.l > 3)
            }) {
                Some("atomic angular momentum above l=3 is unsupported")
            } else {
                None
            };
            atomic_start.fallback_reason = reason.map(str::to_owned);
            if reason.is_none() {
                Some(atomic_start::prepare(&species)?)
            } else {
                None
            }
        } else {
            None
        };
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
            #[cfg(test)]
            force_real_rebuild: false,
            #[cfg(test)]
            force_energy_identity: false,
            #[cfg(test)]
            force_verification_failures: 0,
            gamma_real_exclusion,
            reused_fields: None,
            screen_ffts: Vec::new(),
            energy_identity: false,
            start_density,
            atomic_species,
            atomic_start,
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
    pub fn atomic_start_report(&self) -> &AtomicStartReport {
        &self.atomic_start
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
    /// First SCF input density: the atomic superposition when requested and
    /// available, otherwise the uniform density.
    /// The superposition of atomic charges at this problem's positions (built
    /// only with `start_density = Atomic`; None otherwise or when a species has
    /// no atomic charge samples). Charge equals the electron count.
    pub fn atomic_superposition_density(&self) -> Option<&Array3<f64>> {
        self.start_density.as_ref()
    }
    pub fn starting_density(&self) -> Array3<f64> {
        self.start_density
            .clone()
            .unwrap_or_else(|| self.uniform_density())
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

    /// Measurement access to a fixed-potential map with explicit residual certification.
    /// This does not change the SCF schedule or the ordinary session dispatcher.
    #[doc(hidden)]
    pub fn step_kpoints_with_certification(
        &mut self,
        n: &Array3<f64>,
        initial: Option<&[ComplexMatrix]>,
        tolerance: f64,
        certify: bool,
    ) -> Result<ScfStep, ScfError> {
        self.step_with_occupations(n, initial, tolerance, certify, None)
    }

    /// Hartree and XC fields of an input density and the effective local
    /// potential built from them (ionic + Hartree + XC), in Ha.
    fn input_fields(
        &mut self,
        n: &Array3<f64>,
    ) -> Result<(hartree::HartreeGrid, xc::XcGrid, Array3<f64>), ScfError> {
        let h = hartree::evaluate(n, &self.cell, &mut self.fft)?;
        // M2d: pointwise LDA XC split over ranks (bit-identical potential).
        let m2_xc = crate::m2::active().and_then(|comm| {
            crate::m2::xc(
                &comm,
                &self.options,
                n,
                self.core_density.as_ref(),
                self.density_floor,
                self.cell.volume_bohr3(),
            )
        });
        let x = match m2_xc {
            Some(x) => x?,
            None => self.options.evaluate_xc(
                n,
                self.core_density.as_ref(),
                &mut self.xc_density,
                &self.cell,
                &mut self.fft,
                self.density_floor,
            )?,
        };
        let mut effective = self.ionic_real.clone();
        for ((v, h), x) in effective.iter_mut().zip(&h.potential).zip(&x.potential) {
            *v += h.re + x;
        }
        Ok((h, x, effective))
    }

    /// M2d: occupation-weighted output density of the latest solve from the
    /// distributed ORC eigenvectors (stick FFT on each rank's planes, planes
    /// assembled with one allgatherv). None when the M2 path is not active for this
    /// map (then the replicated `BandDensity` runs). Every rank takes the same
    /// branch: the conditions are replicated data.
    fn m2_density(
        &mut self,
        k_points: &[KPointSolution],
        integer_occupations: bool,
    ) -> Result<Option<Array3<f64>>, ScfError> {
        for state in &mut self.states {
            if let Some(d) = state
                .real_solve
                .as_mut()
                .and_then(|r| r.m2_density.as_mut())
            {
                d.valid = false;
            }
        }
        let Some(comm) = crate::m2::active() else {
            return Ok(None);
        };
        if !integer_occupations || self.states.len() != 1 || k_points.len() != 1 {
            return Ok(None);
        }
        let start = Instant::now();
        let dims = self.grid.dims();
        let volume = self.cell.volume_bohr3();
        let occupations = &self.occupations;
        let state = &mut self.states[0];
        let weight = state.point.weight;
        let grid_dims = state.grid.dims();
        let Some(rs) = state.real_solve.as_mut() else {
            return Ok(None);
        };
        let (Some(dist), Some(x)) = (rs.dist.as_mut(), rs.m2_local.as_ref()) else {
            return Ok(None);
        };
        let bands = x.cols();
        if dist.grid_shape() != dims
            || grid_dims != dims
            || bands != occupations.len()
            || k_points[0].solution.eigenvectors.cols() != bands
        {
            return Ok(None);
        }
        let err = |e: dft_physics::distributed_real_gamma::DistributedRealError| {
            ScfError::Numerical(format!("M2 density: {e}"))
        };
        // Occupied prefix: later bands contribute exactly zero.
        let m = occupations
            .iter()
            .rposition(|&f| f != 0.0)
            .map_or(0, |i| i + 1);
        let n_loc = x.rows();
        let mut xm =
            dft_numerics::RealMatrix::zeros(n_loc, m).map_err(dft_solvers::SolverError::from)?;
        xm.as_mut_slice()
            .copy_from_slice(&x.as_slice()[..n_loc * m]);
        let md = rs
            .m2_density
            .get_or_insert_with(|| crate::m2::M2Density::new(dist.plane_grid_offsets()));
        let points = md.points;
        md.bands = bands;
        md.contributions.clear();
        md.contributions.resize(bands * points, 0.);
        dist.band_densities(
            &comm,
            &xm,
            &occupations[..m],
            volume,
            &mut md.contributions[..m * points],
        )
        .map_err(err)?;
        // Ascending band order on the planes, as the replicated fold.
        let mut local = vec![0.0; points];
        for b in 0..m {
            for (l, c) in local
                .iter_mut()
                .zip(&md.contributions[b * points..(b + 1) * points])
            {
                *l += c;
            }
        }
        let mut full = Array3::zeros(dims);
        dist.gather_planes(&comm, &local, full.as_slice_mut().unwrap())
            .map_err(err)?;
        if full.iter().any(|n: &f64| !n.is_finite()) {
            return Err(ScfError::NonFinite);
        }
        let mut output = Array3::zeros(dims);
        for (out, n) in output.iter_mut().zip(&full) {
            *out += weight * n;
        }
        md.valid = true;
        crate::m2::add(crate::m2::Stage::Density, start);
        crate::m2::count(crate::m2::Count::Density);
        Ok(Some(output))
    }

    /// M2d: occupation-weighted nonlocal energy sum_b f_b x_b^T P D P^T x_b of the
    /// latest M2 solve (one tree_sum); None when no M2 eigenvectors are held.
    fn m2_nonlocal_energy(&mut self) -> Result<Option<f64>, ScfError> {
        let Some(comm) = crate::m2::active() else {
            return Ok(None);
        };
        if self.states.len() != 1 {
            return Ok(None);
        }
        let start = Instant::now();
        let occupations = &self.occupations;
        let Some(rs) = self.states[0].real_solve.as_mut() else {
            return Ok(None);
        };
        let (Some(dist), Some(x)) = (rs.dist.as_mut(), rs.m2_local.as_ref()) else {
            return Ok(None);
        };
        if x.cols() != occupations.len() {
            return Ok(None);
        }
        let m = occupations
            .iter()
            .rposition(|&f| f != 0.0)
            .map_or(0, |i| i + 1);
        let n_loc = x.rows();
        let mut xm =
            dft_numerics::RealMatrix::zeros(n_loc, m).map_err(dft_solvers::SolverError::from)?;
        xm.as_mut_slice()
            .copy_from_slice(&x.as_slice()[..n_loc * m]);
        let mut values = vec![0.0; m];
        dist.nonlocal_band_energies(&comm, &xm, &mut values)
            .map_err(|e| ScfError::Numerical(format!("M2 nonlocal energy: {e}")))?;
        let mut sum = 0.0;
        for (v, f) in values.iter().zip(occupations) {
            if *f != 0.0 {
                sum += f * v;
            }
        }
        crate::m2::add(crate::m2::Stage::NonlocalEnergy, start);
        crate::m2::count(crate::m2::Count::NonlocalEnergy);
        Ok(Some(sum))
    }

    /// Explicit k-weighted nonlocal energy of a map's orbitals with the audit's
    /// arithmetic (B2 final-state audit): equals the `nonlocal` field an
    /// explicit audit of the same map would report, bit for bit.
    pub(crate) fn explicit_nonlocal_energy(
        &self,
        k_points: &[KPointSolution],
    ) -> Result<f64, ScfError> {
        if k_points.len() != self.states.len() {
            return Err(ScfError::InvalidInput("one solution per k point"));
        }
        let mut values = Vec::with_capacity(k_points.len());
        for (state, k) in self.states.iter().zip(k_points) {
            let vectors = &k.solution.eigenvectors;
            let orbitals: Vec<_> = (0..vectors.cols())
                .map(|b| vectors.column(b).to_vec())
                .collect();
            let nonlocal = state.hamiltonian.nonlocal();
            let occupations = self.occupations.as_slice();
            let orbitals = orbitals.as_slice();
            let value = self.pool.install(move || {
                energy::explicit_nonlocal_energy(nonlocal, orbitals, occupations)
            })?;
            values.push((state.point.weight, value));
        }
        Ok(values.iter().map(|(w, v)| w * v).sum())
    }

    /// Drop fields kept by `certificate_shift_estimate` (POTENTIAL-REUSE).
    pub(crate) fn clear_reused_fields(&mut self) {
        self.reused_fields = None;
    }

    /// First-order estimate (Ha) of the occupied-eigenvalue change that the
    /// fixed-point certificate would measure when it re-solves at `density`.
    ///
    /// The certificate compares the occupied eigenvalues of the Hamiltonian built
    /// from `density` with those of the step just taken (potential of its input
    /// density). With dV the change of the effective local potential, the shifts
    /// of a group of bands that is separated from all other bands by more than
    /// `CLUSTER_GAP_HA` are, to first order, the eigenvalues of
    /// diag(eps) + <psi_i|dV|psi_j> over that group (degenerate perturbation
    /// theory); a single band needs only its diagonal element, which the stored
    /// per-band densities of the step provide. Coupling between groups is second
    /// order, about |dV|^2 / gap (1e-9 Ha or less for shifts near the limit), so
    /// the caller applies a margin. Groups larger than `MAX_CLUSTER` or holding an
    /// unoccupied band are skipped (their shifts are then not bounded). The value
    /// is therefore a lower estimate of the true maximum, never an upper bound.
    /// It never certifies anything; it can only predict a certain failure.
    ///
    /// With `keep_fields` the Hartree/XC fields built here are kept for the next
    /// map (`step_*`) at the same density (POTENTIAL-REUSE), so a certificate
    /// solve that follows a passing screen does not build them again.
    pub(crate) fn certificate_shift_estimate(
        &mut self,
        density: &Array3<f64>,
        k_points: &[KPointSolution],
        keep_fields: bool,
    ) -> Result<Option<f64>, ScfError> {
        const CLUSTER_GAP_HA: f64 = 1e-3;
        const MAX_CLUSTER: usize = 24;
        self.reused_fields = None;
        if k_points.len() != self.states.len() {
            return Ok(None);
        }
        let (h, x, effective) = self.input_fields(density)?;
        let standard = effective.as_standard_layout().into_owned();
        if keep_fields {
            self.reused_fields = Some(ReusedFields {
                density: density.clone(),
                hartree: h,
                xc: x,
                effective,
            });
        }
        let effective = standard;
        let target = effective.as_slice().expect("standard layout");
        let points = target.len();
        let weight = self.cell.volume_bohr3() / points as f64;
        let bands = self.occupations.len();
        let mut bound = 0.0_f64;
        let mut delta = vec![0.0; points];
        let mut integrals = vec![0.0; bands];
        let mut psi: Vec<Array3<Complex64>> = Vec::new();
        let mut entries: Vec<(usize, usize, Complex64)> = Vec::new();
        for (state, k) in self.states.iter_mut().zip(k_points) {
            let current = state.hamiltonian.local_potential_values();
            let eigenvalues = &k.solution.eigenvalues;
            if current.len() != points || eigenvalues.len() != bands {
                return Ok(None);
            }
            for ((d, t), c) in delta.iter_mut().zip(target).zip(current) {
                *d = t - c;
            }
            match (
                crate::m2::active(),
                state
                    .real_solve
                    .as_mut()
                    .and_then(|r| r.m2_density.as_mut())
                    .filter(|d| d.valid),
            ) {
                (Some(comm), Some(d)) => d.weighted_integrals(&comm, &delta, &mut integrals)?,
                _ => state
                    .band_density
                    .weighted_integrals(&delta, &mut integrals),
            }
            let mut start = 0;
            while start < bands {
                let mut end = start + 1;
                while end < bands && eigenvalues[end] - eigenvalues[end - 1] < CLUSTER_GAP_HA {
                    end += 1;
                }
                let size = end - start;
                if self.occupations[start..end].iter().any(|f| *f <= 0.0) || size > MAX_CLUSTER {
                    start = end;
                    continue;
                }
                if size == 1 {
                    bound = bound.max((integrals[start] * weight / self.occupations[start]).abs());
                    start = end;
                    continue;
                }
                // Cluster: build diag(eps) + M with M_ij = mean_r conj(psi_i) dV psi_j.
                psi.resize_with(size, || Array3::zeros(state.grid.dims()));
                // One band transform per pool worker (own FFT plan each) and the
                // cluster entries split over the pool; every transform and every
                // entry keeps its serial arithmetic, so the estimate is bit-identical
                // for any thread count.
                if self.screen_ffts.is_empty() {
                    self.screen_ffts = (0..self.options.threads.max(1))
                        .map(|_| FftPlan3d::new(self.grid.dims()))
                        .collect::<Result<Vec<_>, _>>()?;
                }
                let vectors = &k.solution.eigenvectors;
                let grid = &state.grid;
                let mut jobs: Vec<_> = psi.iter_mut().zip(start..end).collect();
                ordered_chunks(
                    &mut jobs,
                    &mut self.screen_ffts,
                    &self.pool,
                    |_, (slot, b), fft| {
                        fft.to_real(&grid.scatter(vectors.column(*b))?, slot)?;
                        Ok(())
                    },
                    |()| {},
                )?;
                drop(jobs);
                entries.clear();
                for j in 0..size {
                    for i in 0..=j {
                        entries.push((i, j, Complex64::default()));
                    }
                }
                let slices: Vec<&[Complex64]> = psi
                    .iter()
                    .map(|p| p.as_slice().expect("standard layout"))
                    .collect();
                let (slices, delta) = (slices.as_slice(), delta.as_slice());
                let entries_slice = entries.as_mut_slice();
                self.pool
                    .install(move || cluster_entries(slices, delta, entries_slice));
                let mut matrix = ComplexMatrix::zeros(size, size)
                    .map_err(|e| ScfError::Numerical(e.to_string()))?;
                for &(i, j, sum) in &entries {
                    let mut value = sum / points as f64;
                    if i == j {
                        value = Complex64::new(value.re + eigenvalues[start + i], 0.0);
                    }
                    matrix.column_mut(j)[i] = value;
                    matrix.column_mut(i)[j] = value.conj();
                }
                let (shifted, _) = dft_numerics::linalg::hermitian_eigh(&matrix)
                    .map_err(|e| ScfError::Numerical(e.to_string()))?;
                for (k, lambda) in shifted.iter().enumerate() {
                    bound = bound.max((lambda - eigenvalues[start + k]).abs());
                }
                start = end;
            }
        }
        if bound.is_finite() {
            Ok(Some(bound))
        } else {
            Ok(None)
        }
    }

    fn step_with_occupations(
        &mut self,
        n: &Array3<f64>,
        initial: Option<&[ComplexMatrix]>,
        tolerance: f64,
        certify: bool,
        mut finite: Option<&mut finite_temperature::OccupationWorkspace>,
    ) -> Result<ScfStep, ScfError> {
        #[cfg(test)]
        if self.force_real_rebuild {
            for state in &mut self.states {
                state.real_solve = None;
            }
        }
        if finite.is_some()
            && self.options.preconditioner != dft_solvers::DavidsonPreconditioner::Default
        {
            return Err(ScfError::InvalidInput(
                "nondefault preconditioner does not support finite occupations",
            ));
        }
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
        let (h, x, effective) = match self.reused_fields.take() {
            Some(r)
                if r.density.shape() == n.shape()
                    && r.density
                        .iter()
                        .zip(n.iter())
                        .all(|(a, b)| a.to_bits() == b.to_bits()) =>
            {
                (r.hartree, r.xc, r.effective)
            }
            _ => self.input_fields(n)?,
        };
        ordered_solves(
            &mut self.states,
            self.options.threads,
            &self.pool,
            |_, state| {
                state.hamiltonian.set_local_potential(effective.clone())?;
                Ok(())
            },
        )?;
        // Only a cold map may create atomic vectors. Existing restart/warm
        // vectors pass unchanged through both solver bridges.
        let atomic_initial = if self.options.start_wavefunctions == StartWavefunctions::Atomic {
            if initial.is_some() {
                if self.atomic_start.initial_h_applications == 0
                    && self.atomic_start.fallback_reason.is_none()
                {
                    self.atomic_start.fallback_reason =
                        Some("supplied initial vectors retained".into());
                }
                None
            } else if finite.is_some() {
                self.atomic_start.fallback_reason = Some("Fermi-Dirac occupations".into());
                None
            } else if let Some(species) = &self.atomic_species {
                let atomic_clock = Instant::now();
                let mut guesses = Vec::with_capacity(self.states.len());
                for state in &mut self.states {
                    let span = atomic_start::span(
                        &self.cell,
                        &state.basis,
                        &state.grid,
                        species,
                        &self.options,
                        self.occupations.len(),
                        dft_numerics::switches::on(dft_numerics::switches::Id::PHASE_THREADS)
                            .then_some(&*self.pool),
                    )?;
                    // ATOMIC-REAL: H Q on the real Gamma path when eligible.
                    let hq = if crate::gamma_real::atomic_real() {
                        crate::gamma_real::atomic_apply_real(
                            &state.hamiltonian,
                            &state.basis,
                            span.vectors(),
                            self.occupations.len(),
                            &self.options,
                            &self.pool,
                            self.gamma_real_exclusion,
                            &mut state.real_solve,
                        )?
                    } else {
                        None
                    };
                    let (guess, report) = atomic_start::ritz(
                        &span,
                        &state.hamiltonian,
                        self.occupations.len(),
                        state.point.fractional == [0.; 3],
                        hq,
                    )?;
                    self.atomic_start.atomic_orbitals += report.atomic_orbitals;
                    self.atomic_start.initial_h_applications += report.initial_h_applications;
                    self.atomic_start.gram_rotation_macs += report.gram_rotation_macs;
                    self.atomic_start.mgs_projection_macs += report.mgs_projection_macs;
                    self.atomic_start.radial_samples += report.radial_samples;
                    self.atomic_start.evd_orders.extend(report.evd_orders);
                    guesses.push(guess);
                }
                crate::phase_clock::add_atomic_start(atomic_clock);
                Some(guesses)
            } else {
                None
            }
        } else {
            None
        };
        let initial = initial.or(atomic_initial.as_deref());
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
        // BAANDI_EMPTY_LOOSE (default on, `off` disables): in SCF steps whose
        // eigen tolerance is above the floor 0.1 * eigen_tolerance (every step except
        // the last ones and the certificate re-solve) only the occupied bands are
        // converged strictly; the first empty band takes the buffer threshold.
        let policy = band_tolerances(
            self.occupations.iter().filter(|&&f| f > 0.0).count(),
            self.occupations.len(),
            tolerance,
            &self.options,
        );
        // A single k point still has an intra-k budget. Keep its complete solve
        // resident on that pool, so each band/GEMM phase does not enter from
        // the external caller and wait on a LockLatch again.
        // The distributed solver issues MPI collectives, which must come from the
        // thread that owns the communicator: solve the k points sequentially on the
        // caller thread when one is active (intra-rank threads still serve GEMM/FFT).
        let k_points =
            if self.options.threads <= 1 || crate::parallel::SharedComm::active().is_some() {
                let mut k_points = Vec::new();
                for (i, state) in self.states.iter_mut().enumerate() {
                    let preconditioner = DiagonalPreconditioner::new(
                        state.hamiltonian.kinetic_diagonal().to_vec(),
                        0.1,
                    )?
                    .with_variant(self.options.preconditioner);
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
                        &mut state.real_solve,
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
                        )?
                        .with_variant(self.options.preconditioner);
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
                            &mut state.real_solve,
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
        // M2d: distributed density on the M2 path (None: replicated band density).
        let m2_output = self.m2_density(&k_points, finite.is_none())?;
        let m2 = m2_output.is_some();
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
                // DENSITY-PAIR applies only to Hermitian-completed vectors of
                // this step's real Gamma solve; the M2 path builds its density
                // from the distributed eigenvectors.
                if !m2 {
                    state.band_density.compute_with(
                        &orbitals,
                        finite.map_or(occupations.as_slice(), |w| w.block(0, i)),
                        &state.grid,
                        volume,
                        state.gamma_real.gamma_real_used,
                    )?;
                }
                Ok(orbitals)
            },
        )?;
        let mut output = if let Some(output) = m2_output {
            output
        } else {
            let mut output = Array3::zeros(self.grid.dims());
            // Every k workspace owns its density output. Fold in fixed k order after
            // all jobs finish; no allocated partial arrays or worker locks.
            for state in &self.states {
                for (out, n) in output.iter_mut().zip(state.band_density.density()) {
                    *out += state.point.weight * n;
                }
            }
            output
        };
        if let Some(symmetry) = &mut self.symmetry {
            symmetry.apply(&mut output, &mut self.fft)?;
        }
        timings.density = start.elapsed();
        let start = Instant::now();
        let m2_xc = crate::m2::active().and_then(|comm| {
            crate::m2::xc(
                &comm,
                &self.options,
                &output,
                self.core_density.as_ref(),
                self.density_floor,
                volume,
            )
        });
        let output_xc = if let Some(x) = m2_xc {
            Some(x?)
        } else if self.core_density.is_some()
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
        // B2: Fermi-Dirac occupations fall back to the explicit audit.
        let identity = self.energy_identity && finite.is_none();
        // M2d: nonlocal energy from the distributed projections (one tree_sum).
        let m2_nonlocal = if m2 && !identity {
            self.m2_nonlocal_energy()?
        } else {
            None
        };
        let mut audits = Vec::new();
        ordered_chunks(
            &mut self.states,
            &mut self.worker_ffts,
            &self.pool,
            |i, state, _fft| {
                let k = &k_points[i];
                // NL-ENERGY-REAL: values of this step's real solve, if computed.
                let partials = state
                    .real_solve
                    .as_ref()
                    .filter(|_| state.gamma_real.gamma_real_used)
                    .and_then(|r| r.nonlocal_partials());
                let audit = energy::evaluate_shared_with_partials(
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
                    !identity && m2_nonlocal.is_none(),
                    partials,
                )?;
                Ok((state.point.weight, audit))
            },
            |audit| audits.push(audit),
        )?;
        if let (Some(value), Some((_, audit))) = (m2_nonlocal, audits.first_mut()) {
            audit.nonlocal = value;
        }
        let mut energy = energy::weighted_audits(&audits);
        if identity {
            energy::apply_band_sum_identity(&mut energy);
        }
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
            nonlocal_from_identity: identity,
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
/// Band-tolerance policy shared by every Davidson driver (unpolarized, collinear
/// spin, finite temperature): occupied bands strict; in steps whose eigen
/// tolerance is above its floor (0.1 x eigen_tolerance) the first empty band takes
/// the buffer threshold (EMPTY_LOOSE) and the buffer floor is BUFFER_TOL; the
/// standard profile uses its own empty-state threshold.
pub(crate) fn band_tolerances(
    occupied_count: usize,
    bands: usize,
    tolerance: f64,
    options: &ScfOptions,
) -> BandTolerances {
    let intermediate = empty_loose()
        && options.standard.is_none()
        && tolerance > 0.1 * options.eigen_tolerance * (1.0 + 1e-9);
    BandTolerances {
        strict_bands: (occupied_count + usize::from(!intermediate))
            .min(bands)
            .max(1),
        buffer_tolerance: if options.standard.is_some() {
            StandardProfile::empty_threshold(tolerance)
        } else if intermediate {
            tolerance.max(buffer_loose_floor())
        } else {
            tolerance.max(1e-5)
        },
    }
}
/// Empty-band looseness: default ON since the adoption rule (dev08 EMPTY-LOOSE:
/// strict perturbed si32 -13.1 / -13.7 percent wall at T=1 / T=8, 10 of 10 seeds,
/// Davidson work -60.7 percent; certificate verdicts and certified states
/// unchanged). `BAANDI_EMPTY_LOOSE=off` restores the strict first empty band in
/// every step; see `step_with_occupations`.
fn empty_loose() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| dft_numerics::switches::on(dft_numerics::switches::Id::EMPTY_LOOSE))
}
/// Residual floor of the empty buffer bands in SCF steps whose eigen tolerance is
/// above its floor: default 1e-3 since the adoption rule (dev08 BUFFER-TOL: si32
/// strict perturbed -7.54 percent wall, 10 of 10 seeds, optimum interior: 1e-2 is
/// past the turn; certificates stay strict). `BAANDI_BUFFER_TOL=<x>` (1e-5 to
/// 1e-1) overrides; 1e-5 restores the old floor.
fn buffer_loose_floor() -> f64 {
    static FLOOR: std::sync::OnceLock<f64> = std::sync::OnceLock::new();
    *FLOOR.get_or_init(|| {
        dft_numerics::switches::number(dft_numerics::switches::Id::BUFFER_TOL).unwrap_or(1e-3)
    })
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

    mod persistent_real_state {
        include!("../tests/common/persistent_real_state_unit.rs");
    }

    /// The first-order estimate of the certificate's occupied-eigenvalue change
    /// (singles and a threefold-degenerate occupied cluster at an s + p atom in a
    /// cubic cell) agrees with the real re-solve once the SCF is nearly converged.
    #[test]
    fn certificate_shift_estimate_matches_resolve() {
        let potential = pool_fixture::potential(8.0, 1.0, true, 101);
        let cell = Cell::from_user_rows_bohr([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]]).unwrap();
        let mut problem = ScfProblem::new(
            &cell,
            &[Atom::new("Xx", [0.0, 0.0, 0.0])],
            &[SpeciesPotential {
                label: "Xx",
                potential: &potential,
            }],
            ScfOptions {
                threads: 1,
                cutoff_hartree: 3.0,
                bands: Some(6),
                ..Default::default()
            },
        )
        .unwrap();
        let occupied = problem.occupations.iter().filter(|f| **f > 0.0).count();
        assert!(occupied >= 3, "need an s + p occupied set");
        let mut n = problem.uniform_density();
        let mut warm: Option<Vec<ComplexMatrix>> = None;
        for _ in 0..14 {
            let step = problem.step_kpoints(&n, warm.as_deref(), 1e-11).unwrap();
            warm = Some(
                step.k_points
                    .iter()
                    .map(|k| k.solution.eigenvectors.clone())
                    .collect(),
            );
            let out = step.density;
            n.zip_mut_with(&out, |a, b| *a += 0.5 * (b - *a));
        }
        let step = problem.step_kpoints(&n, warm.as_deref(), 1e-11).unwrap();
        let estimate = problem
            .certificate_shift_estimate(&step.density, &step.k_points, false)
            .unwrap()
            .unwrap();
        let warm: Vec<_> = step
            .k_points
            .iter()
            .map(|k| k.solution.eigenvectors.clone())
            .collect();
        let check = problem
            .step_kpoints(&step.density, Some(&warm), 1e-11)
            .unwrap();
        let actual = check
            .k_points
            .iter()
            .zip(&step.k_points)
            .flat_map(|(a, b)| a.solution.eigenvalues.iter().zip(&b.solution.eigenvalues))
            .take(occupied)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        println!("estimate {estimate:.6e} actual {actual:.6e}");
        assert!(actual > 1e-12, "fixture must still move: {actual}");
        assert!(
            (estimate - actual).abs() <= 2e-3 * actual + 1e-13,
            "estimate {estimate:e} actual {actual:e}"
        );
    }

    /// POTENTIAL-REUSE: the certificate solve that follows a screen at the same
    /// density gives bit-identical eigenvalues, density and energy whether or not
    /// it reuses the screen's Hartree/XC fields; the kept fields are consumed by
    /// that one map and never survive it. The estimate itself (pool-parallel
    /// cluster transforms and entries) is bit-identical for 1 and 4 threads.
    #[test]
    fn certificate_potential_reuse_is_bit_identical() {
        let potential = pool_fixture::potential(8.0, 1.0, true, 101);
        // Cubic: the occupied p shell is a threefold cluster (matrix path).
        let cell = Cell::from_user_rows_bohr([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]]).unwrap();
        let run = |keep: bool, threads: usize| {
            let mut problem = ScfProblem::new(
                &cell,
                &[Atom::new("Xx", [0.0, 0.0, 0.0])],
                &[SpeciesPotential {
                    label: "Xx",
                    potential: &potential,
                }],
                ScfOptions {
                    threads,
                    cutoff_hartree: 3.0,
                    bands: Some(6),
                    ..Default::default()
                },
            )
            .unwrap();
            let mut n = problem.uniform_density();
            let mut warm: Option<Vec<ComplexMatrix>> = None;
            for _ in 0..4 {
                let step = problem.step_kpoints(&n, warm.as_deref(), 1e-9).unwrap();
                warm = Some(
                    step.k_points
                        .iter()
                        .map(|k| k.solution.eigenvectors.clone())
                        .collect(),
                );
                n.zip_mut_with(&step.density, |a, b| *a += 0.5 * (b - *a));
            }
            let step = problem.step_kpoints(&n, warm.as_deref(), 1e-9).unwrap();
            let estimate = problem
                .certificate_shift_estimate(&step.density, &step.k_points, keep)
                .unwrap()
                .unwrap();
            assert_eq!(problem.reused_fields.is_some(), keep);
            let warm: Vec<_> = step
                .k_points
                .iter()
                .map(|k| k.solution.eigenvectors.clone())
                .collect();
            let check = problem
                .step_kpoints(&step.density, Some(&warm), 1e-10)
                .unwrap();
            assert!(
                problem.reused_fields.is_none(),
                "kept fields must be consumed"
            );
            let eig: Vec<u64> = check.k_points[0]
                .solution
                .eigenvalues
                .iter()
                .map(|e| e.to_bits())
                .collect();
            let density: Vec<u64> = check.density.iter().map(|v| v.to_bits()).collect();
            (
                estimate.to_bits(),
                eig,
                density,
                check.energy.total_hartree().to_bits(),
            )
        };
        let reference = run(true, 1);
        assert_eq!(reference, run(false, 1));
        assert_eq!(reference, run(true, 4));
    }

    /// B2: one map with the band-sum identity differs from the explicit audit
    /// only in `nonlocal` (to rounding of the eigenpairs), every other audit field
    /// is bit-identical, and the explicit final-state audit reproduces the
    /// explicit map's nonlocal energy bit for bit. A full SCF with the identity
    /// in non-final iterations reaches the same certified state and energy.
    #[test]
    fn energy_identity_matches_explicit_audit() {
        let potential = pool_fixture::potential(8.0, 1.0, true, 101);
        let cell =
            Cell::from_user_rows_bohr([[6., 0., 0.], [0.4, 6.2, 0.], [0.3, 0., 5.8]]).unwrap();
        let make = || {
            ScfProblem::new(
                &cell,
                &[Atom::new("Xx", [0.0, 0.0, 0.0])],
                &[SpeciesPotential {
                    label: "Xx",
                    potential: &potential,
                }],
                ScfOptions {
                    threads: 1,
                    cutoff_hartree: 3.0,
                    bands: Some(6),
                    ..Default::default()
                },
            )
            .unwrap()
        };
        let mut problem = make();
        let n = problem.uniform_density();
        let explicit = problem.step_kpoints(&n, None, 1e-10).unwrap();
        let mut problem = make();
        problem.energy_identity = true;
        let identity = problem.step_kpoints(&n, None, 1e-10).unwrap();
        assert!(identity.nonlocal_from_identity && !explicit.nonlocal_from_identity);
        let (a, b) = (&explicit.energy, &identity.energy);
        for (x, y) in [
            (a.kinetic, b.kinetic),
            (a.local, b.local),
            (a.hartree, b.hartree),
            (a.xc, b.xc),
            (a.band_sum, b.band_sum),
            (a.double_counting, b.double_counting),
            (a.band_form, b.band_form),
        ] {
            assert_eq!(x.to_bits(), y.to_bits());
        }
        let difference = b.nonlocal - a.nonlocal;
        println!(
            "nonlocal explicit {:.15e} identity {:.15e} difference {difference:.3e}",
            a.nonlocal, b.nonlocal
        );
        assert!(a.nonlocal.abs() > 1e-3, "fixture needs a nonlocal term");
        assert!(
            difference.abs() < 1e-9,
            "identity difference {difference:e}"
        );
        // The identity total is the band-structure form by construction.
        assert!((b.total_hartree() - b.band_form).abs() < 1e-11);
        let audit = problem
            .explicit_nonlocal_energy(&identity.k_points)
            .unwrap();
        // NL_ENERGY_REAL (default on since the adoption rule) computes the step audit
        // from a real projector GEMM; the explicit function sums complex projections:
        // the two agree to rounding, not bit for bit.
        assert!(
            (audit - a.nonlocal).abs() < 1e-13,
            "{audit:e} {:e}",
            a.nonlocal
        );

        let off = make().run().unwrap();
        let mut on = make();
        on.force_energy_identity = true;
        let on = on.run().unwrap();
        assert!(off.energy_identity.is_none());
        let check = on.energy_identity.clone().expect("identity report");
        println!("final {check:?}");
        assert!(check.difference_ha.abs() < 1e-9);
        assert_eq!(on.history.len(), off.history.len());
        // With the rounding-class defaults (DGKS skipping, real-GEMM audits) the two
        // maps can differ by an ulp; the certified final state must agree to 1e-12.
        assert!(
            (on.total_hartree() - off.total_hartree()).abs() < 1e-12,
            "final explicit audit must equal the off run"
        );
        let worst = on
            .density
            .iter()
            .zip(off.density.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        assert!(worst < 1e-10, "density differs by {worst:e}");
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
#[path = "verification_tests.rs"]
mod verification_tests;

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
