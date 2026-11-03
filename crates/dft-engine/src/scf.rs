use crate::Functional;
use crate::energy::{self, AuditInput, EnergyAudit};
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_geometry::{
    Atom, Cell,
    kpoints::{self, KPoint},
};
use dft_physics::{
    data::*,
    density,
    ewald::{self, EwaldOptions, Ion},
    hamiltonian::{self, Hamiltonian, NonlocalOperator, Species},
    hartree, xc,
};
use dft_solvers::{
    BandTolerances, Davidson, DiagonalPreconditioner, EigenSolution,
    kerker::{Kerker, MixerHistory, PreconditionedPulay},
    mixing::{LinearMixer, PulayMixer},
};
use dft_upf::Pseudopotential;
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
    /// Experimental opt-in hybrid; Davidson remains the default.
    pub eigensolver: dft_solvers::SolverChoice,
    /// Independent k-point workers; zero means one, capped at the k-point count.
    pub threads: usize,
    pub functional: Functional,
    pub allow_functional_override: bool,
    pub cutoff_hartree: f64,
    pub bands: Option<usize>,
    pub extra_bands: usize,
    pub grid_override: Option<[usize; 3]>,
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
    /// Minimum working capacity; the engine also reserves six columns per band.
    pub davidson_subspace: usize,
    pub seed: u64,
    pub divergence_factor: f64,
    pub divergence_patience: usize,
}
impl Default for ScfOptions {
    fn default() -> Self {
        Self {
            eigensolver: dft_solvers::SolverChoice::default(),
            threads: 1,
            functional: Functional::LdaPz,
            allow_functional_override: false,
            cutoff_hartree: 5.0,
            bands: None,
            extra_bands: 4,
            grid_override: None,
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
            seed: 2026,
            divergence_factor: 20.0,
            divergence_patience: 3,
        }
    }
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

/// Cached ionic terms, Hamiltonian application plan, and density FFT plan.
/// `step` is exposed for independent post-SCF fixed-point verification.
pub struct ScfProblem {
    cell: Cell,
    states: Vec<KState>,
    grid: FftGridSpec,
    ionic_g: Array3<Complex64>,
    ionic_real: Array3<f64>,
    fft: FftPlan3d,
    // One reusable FFT workspace per active k worker, shared across phases.
    worker_ffts: Vec<FftPlan3d>,
    ewald: ewald::EwaldEnergy,
    occupations: Vec<f64>,
    electrons: f64,
    options: ScfOptions,
    setup_time: Duration,
    minimum_gap: Option<f64>,
}
struct KState {
    point: KPoint,
    basis: PlaneWaveBasis,
    grid: FftGridSpec,
    hamiltonian: Hamiltonian,
    band_density: Option<crate::band_density::BandDensity>,
}
/// Disjoint state ownership needs Send, never Sync for Hamiltonian scratch.
/// Join every worker before inspecting results in ascending k order. Catch each
/// solve panic so later indices in the same group are still evaluated.
fn ordered_solves<T: Send, R: Send>(
    states: &mut [T],
    requested: usize,
    solve: impl Fn(usize, &mut T) -> Result<R, ScfError> + Sync,
) -> Result<Vec<R>, ScfError> {
    ordered_solves_at(states, requested, 0, solve)
}

fn ordered_solves_at<T: Send, R: Send>(
    states: &mut [T],
    requested: usize,
    first: usize,
    solve: impl Fn(usize, &mut T) -> Result<R, ScfError> + Sync,
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
    std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(workers);
        let mut remaining = states;
        let mut offset = first;
        for worker in 0..workers {
            let count = remaining.len().div_ceil(workers - worker);
            let (group, rest) = remaining.split_at_mut(count);
            remaining = rest;
            let start = offset;
            offset += count;
            let solve = &solve;
            handles.push(scope.spawn(move || {
                group
                    .iter_mut()
                    .enumerate()
                    .map(|(local, state)| {
                        let index = start + local;
                        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                            solve(index, state)
                        }))
                        .unwrap_or_else(|_| {
                            Err(ScfError::Numerical(format!(
                                "k index {index}: solve panicked"
                            )))
                        })
                        .map_err(|e| ScfError::Numerical(format!("k index {index}: {e}")))
                    })
                    .collect::<Vec<_>>()
            }));
        }
        // Do not short circuit joins: scope must never rethrow an unjoined panic.
        let joined: Vec<_> = handles.into_iter().map(|h| h.join()).collect();
        let mut results = Vec::new();
        for group in joined {
            let group =
                group.map_err(|_| ScfError::Numerical("k worker panicked outside solve".into()))?;
            for result in group {
                results.push(result?);
            }
        }
        Ok(results)
    })
}

/// Bound live partial results by the workspace count. Fold each completed
/// chunk in ascending global k order before starting the next chunk; neither
/// chunk boundaries nor completion order can change floating-point association.
fn ordered_chunks<T: Send, W: Send, R: Send>(
    states: &mut [T],
    workspaces: &mut [W],
    work: impl Fn(usize, &mut T, &mut W) -> Result<R, ScfError> + Sync,
    mut fold: impl FnMut(R),
) -> Result<(), ScfError> {
    assert!(!workspaces.is_empty());
    let width = workspaces.len();
    for (chunk, states) in states.chunks_mut(width).enumerate() {
        let mut jobs: Vec<_> = states.iter_mut().zip(workspaces.iter_mut()).collect();
        let partials =
            ordered_solves_at(&mut jobs, width, chunk * width, |i, (state, workspace)| {
                work(i, state, workspace)
            })?;
        for partial in partials {
            fold(partial);
        }
    }
    Ok(())
}

/// Output-oriented completion, including a nonconverged final evaluated state.
/// Numerical failures still return `ScfError`; this is not a convergence bypass.
pub struct ScfOutcome {
    pub converged: bool,
    pub result: ScfResult,
}

pub struct ScfStep {
    /// First-k compatibility view; `k_points` contains every solution.
    pub solution: EigenSolution,
    /// Every k block in input quadrature order.
    pub k_points: Vec<KPointSolution>,
    pub density: Array3<f64>,
    pub energy: EnergyAudit,
    pub timings: PhaseTimings,
}
impl ScfProblem {
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
        )
    }
    fn construct(
        cell: &Cell,
        atoms: &[Atom],
        potentials: &[SpeciesPotential<'_>],
        options: ScfOptions,
        points: &[KPoint],
        minimum_gap: Option<f64>,
    ) -> Result<Self, ScfError> {
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
        if !count.is_multiple_of(2) {
            return Err(ScfError::OddElectronCount(electrons));
        }
        let occupied = count / 2;
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
        let gamma = PlaneWaveBasis::new(cell, Hartree(options.cutoff_hartree))?;
        let density_basis = gamma.density_basis(cell)?;
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
        let ionic_g = dg.scatter(&hamiltonian::ionic_local(
            &density_basis,
            cell.volume_bohr3(),
            &species,
        )?)?;
        let mut fft = FftPlan3d::new(grid.dims())?;
        let mut ionic = Array3::zeros(grid.dims());
        fft.to_real(&ionic_g, &mut ionic)?;
        let scale = ionic.iter().map(|z| z.re.abs()).fold(1.0, f64::max);
        if ionic.iter().any(|z| z.im.abs() > 1e-12 * scale) {
            return Err(ScfError::Numerical("ionic imaginary potential".into()));
        }
        let ionic_real = ionic.mapv(|z| z.re);
        let outer = options.threads.max(1).min(points.len());
        let inner = (options.threads.max(1) / outer).max(1);
        let mut construction_points = points.to_vec();
        let states = ordered_solves(&mut construction_points, options.threads, |_, point| {
            let basis =
                PlaneWaveBasis::at_k(cell, Hartree(options.cutoff_hartree), point.fractional)?;
            if bands > basis.vectors().len()
                || (options.davidson_subspace.min(basis.vectors().len()) <= bands
                    && bands != basis.vectors().len())
            {
                return Err(ScfError::InvalidInput(
                    "bands/subspace incompatible with k basis",
                ));
            }
            let k_grid = FftGridSpec::new(&basis, &density_basis, Some(grid.dims()))?;
            let hamiltonian = Hamiltonian::new(
                cell,
                &basis,
                ionic_real.clone(),
                NonlocalOperator::new(&basis, cell.volume_bohr3(), &species)?,
            )?
            .with_band_threads(inner)?;
            let band_density = Some(crate::band_density::BandDensity::new(
                &k_grid, bands, inner,
            )?);
            Ok(KState {
                point: KPoint {
                    fractional: basis.k_fractional(),
                    weight: point.weight,
                },
                basis,
                grid: k_grid,
                hamiltonian,
                band_density,
            })
        })?;
        let worker_ffts = (0..options.threads.max(1).min(states.len()))
            .map(|_| FftPlan3d::new(grid.dims()))
            .collect::<Result<Vec<_>, _>>()?;
        let ewald = ewald::energy(cell, &ions, EwaldOptions::default())?;
        let occupations = (0..bands)
            .map(|b| if b < occupied { 2.0 } else { 0.0 })
            .collect();
        Ok(Self {
            cell: cell.clone(),
            states,
            grid,
            ionic_g,
            ionic_real,
            fft,
            worker_ffts,
            ewald,
            occupations,
            electrons,
            options,
            setup_time: start.elapsed(),
            minimum_gap,
        })
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
        if initial.is_some_and(|v| v.len() != self.states.len()) {
            return Err(ScfError::InvalidInput("one warm start required per k"));
        }
        if n.shape() != self.grid.dims() {
            return Err(ScfError::InvalidInput("input density grid"));
        }
        let start = Instant::now();
        let h = hartree::evaluate(n, &self.cell, &mut self.fft)?;
        let x = self
            .options
            .functional
            .evaluate(n, &self.cell, &mut self.fft)?;
        let mut effective = self.ionic_real.clone();
        for ((v, h), x) in effective.iter_mut().zip(&h.potential).zip(&x.potential) {
            *v += h.re + x;
        }
        ordered_solves(&mut self.states, self.options.threads, |_, state| {
            state.hamiltonian.set_local_potential(effective.clone())?;
            Ok(())
        })?;
        let mut timings = PhaseTimings {
            potential: start.elapsed(),
            ..Default::default()
        };
        let start = Instant::now();
        let solver = Davidson {
            bands: self.occupations.len(),
            tolerance,
            max_iterations: self.options.davidson_max_iterations,
            max_subspace: self
                .options
                .davidson_subspace
                .max(6 * self.occupations.len()),
            block_size: self.occupations.len(),
            seed: self.options.seed,
        };
        // Occupied states determine density; the first empty state determines
        // the sampled gap. Higher empty states are only subspace buffers.
        let policy = BandTolerances {
            strict_bands: (self.occupations.iter().filter(|&&f| f > 0.0).count() + 1)
                .min(self.occupations.len()),
            buffer_tolerance: tolerance.max(1e-5),
        };
        let k_points = if self.options.threads.max(1).min(self.states.len()) == 1 {
            let mut k_points = Vec::new();
            for (i, state) in self.states.iter().enumerate() {
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
                let solution = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
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
                }))
                .map_err(|_| ScfError::Numerical(format!("k index {i}: solve panicked")))?
                .map_err(|e| ScfError::Numerical(format!("k index {i}: {e}")))?;
                k_points.push(KPointSolution {
                    point: state.point,
                    solution,
                });
            }
            k_points
        } else {
            let bands = self.occupations.len();
            let free_guess = self.options.free_electron_initial_guess;
            ordered_solves(&mut self.states, self.options.threads, |i, state| {
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
                let solution = self
                    .options
                    .eigensolver
                    .solve_with_band_tolerances_and_certification(
                        &solver,
                        &state.hamiltonian,
                        initial.map(|v| &v[i]).or(guess.as_ref()),
                        Some(&preconditioner),
                        policy,
                        certify,
                    )?;
                Ok(KPointSolution {
                    point: state.point,
                    solution,
                })
            })?
        };
        timings.diagonalization = start.elapsed();
        let start = Instant::now();
        let mut orbitals = Vec::with_capacity(k_points.len());
        let mut output = Array3::zeros(self.grid.dims());
        ordered_chunks(
            &mut self.states,
            &mut self.worker_ffts,
            |i, state, fft| {
                let vectors = &k_points[i].solution.eigenvectors;
                let orbitals: Vec<_> = (0..vectors.cols())
                    .map(|b| vectors.column(b).to_vec())
                    .collect();
                let nk = if let Some(workspace) = &mut state.band_density {
                    workspace.compute(
                        &orbitals,
                        &self.occupations,
                        &state.grid,
                        self.cell.volume_bohr3(),
                    )?
                } else {
                    density::from_orbitals(
                        &orbitals,
                        &self.occupations,
                        &state.grid,
                        self.cell.volume_bohr3(),
                        fft,
                    )?
                };
                Ok((state.point.weight, nk, orbitals))
            },
            |(weight, nk, k_orbitals)| {
                for (out, n) in output.iter_mut().zip(nk) {
                    *out += weight * n;
                }
                orbitals.push(k_orbitals);
            },
        )?;
        timings.density = start.elapsed();
        let start = Instant::now();
        let mut audits = Vec::new();
        ordered_chunks(
            &mut self.states,
            &mut self.worker_ffts,
            |i, state, fft| {
                let k = &k_points[i];
                let audit = energy::evaluate_with_functional(
                    AuditInput {
                        cell: &self.cell,
                        basis: &state.basis,
                        grid: &state.grid,
                        ionic_g: &self.ionic_g,
                        nonlocal: state.hamiltonian.nonlocal(),
                        orbitals: &orbitals[i],
                        occupations: &self.occupations,
                        eigenvalues: &k.solution.eigenvalues,
                        eigen_residuals: &k.solution.residuals,
                        input_density: n,
                        output_density: &output,
                        input_hartree: &h,
                        input_xc: &x,
                        ewald: self.ewald.energy.0,
                    },
                    fft,
                    self.options.functional,
                )?;
                Ok((state.point.weight, audit))
            },
            |audit| audits.push(audit),
        )?;
        let energy = energy::weighted_audits(&audits);
        timings.energy = start.elapsed();
        Ok(ScfStep {
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
    fn run_impl(mut self, retain_last: bool) -> Result<ScfOutcome, ScfError> {
        let start = Instant::now();
        let mut mixer: Box<dyn MixerHistory> = match self.options.mixing {
            Mixing::Linear { beta } => Box::new(LinearMixer { beta }),
            Mixing::Pulay { beta, history } => Box::new(PulayMixer::new(beta, history)?),
            Mixing::KerkerPulay {
                beta,
                history,
                q0_inverse_bohr,
            } => {
                let dims = self.grid.dims();
                let mut fft = FftPlan3d::new(dims)?;
                let mut source = Array3::<Complex64>::zeros(dims);
                let mut target = source.clone();
                let transform = move |values: &mut [Complex64], inverse: bool| {
                    source.as_slice_mut().unwrap().copy_from_slice(values);
                    let result = if inverse {
                        fft.to_real(&source, &mut target)
                    } else {
                        fft.to_reciprocal(&source, &mut target)
                    };
                    result.map_err(|e| dft_solvers::SolverError::Application(e.to_string()))?;
                    values.copy_from_slice(target.as_slice().unwrap());
                    Ok(())
                };
                let filter = Kerker::new(
                    dims,
                    self.cell.reciprocal_columns_inverse_bohr(),
                    q0_inverse_bohr,
                    transform,
                )?;
                Box::new(PreconditionedPulay::new(
                    beta,
                    history,
                    dims.iter().product(),
                    filter,
                )?)
            }
        };
        // Validate mixer parameters even if a future system converges immediately.
        let mut n = self.uniform_density();
        mixer.mix(
            n.as_slice().unwrap(),
            n.as_slice().unwrap(),
            self.electrons / self.cell.volume_bohr3(),
        )?;
        mixer.reset();
        let mut history: Vec<IterationRecord> = Vec::new();
        let mut previous: Option<Vec<ComplexMatrix>> = None;
        let weight = self.cell.volume_bohr3() / n.len() as f64;
        let mut best = f64::INFINITY;
        let mut growing = 0;
        let mut mixing_resets = 0;
        for iteration in 1..=self.options.max_iterations {
            if n.iter()
                .any(|n| !n.is_finite() || *n < -xc::NEGATIVE_ROUNDOFF)
            {
                return Err(ScfError::Diverged {
                    reason: "mixed density is nonfinite or negative beyond 1e-14 Bohr^-3".into(),
                    history,
                });
            }
            let tolerance = if self.options.residual_eigen_schedule {
                // Residual is the normalized density L2 norm in atomic units.
                // 0.1 is an empirical scheduling coefficient in atomic units.
                history
                    .last()
                    .map_or(1e-2, |r| 0.1 * r.density_residual)
                    .min(1e-2)
                    .max(0.1 * self.options.eigen_tolerance)
            } else {
                (self.options.initial_eigen_tolerance
                    * self
                        .options
                        .eigen_tolerance_decay
                        .powf((iteration - 1) as f64))
                .max(0.1 * self.options.eigen_tolerance)
            };
            let step = self
                .step_kpoints_with_certification(
                    &n,
                    previous.as_deref(),
                    tolerance,
                    tolerance <= self.options.eigen_tolerance,
                )
                .map_err(|e| match e {
                    ScfError::NonFinite => ScfError::Diverged {
                        reason: "nonfinite SCF arithmetic".into(),
                        history: history.clone(),
                    },
                    other => other,
                })?;
            let residual =
                density_distance(&step.density, &n, self.cell.volume_bohr3(), self.electrons);
            let charge = step.density.sum() * weight;
            let charge_error = (charge - self.electrons).abs() / self.electrons;
            let input_charge_error = (n.sum() * weight - self.electrons).abs() / self.electrons;
            let energy_change = history
                .last()
                .map(|r| step.energy.total_hartree() - r.energy.total_hartree());
            let max_eigen_residual = step
                .k_points
                .iter()
                .flat_map(|k| k.solution.residuals.iter().zip(&self.occupations))
                .filter(|(_, f)| **f > 0.0)
                .map(|(r, _)| *r)
                .fold(0.0, f64::max);
            let record = IterationRecord {
                fixed_point_certificate: None,
                mixing_fraction: 1.0,
                mixer_reset: false,
                density_contraction: history
                    .last()
                    .filter(|r| r.density_residual > 0.0)
                    .map(|r| residual / r.density_residual),
                mixing_history_used: None,
                mixing_resets,
                mixing_raw_recovery: false,
                iteration,
                energy: step.energy.clone(),
                energy_change,
                density_residual: residual,
                max_eigen_residual,
                requested_eigen_tolerance: tolerance,
                davidson_iterations: step.k_points.iter().map(|k| k.solution.iterations).sum(),
                operator_applications: step
                    .k_points
                    .iter()
                    .map(|k| k.solution.operator_applications)
                    .sum(),
                davidson_iterations_per_k: step
                    .k_points
                    .iter()
                    .map(|k| k.solution.iterations)
                    .collect(),
                operator_applications_per_k: step
                    .k_points
                    .iter()
                    .map(|k| k.solution.operator_applications)
                    .collect(),
                charge_error,
                input_charge_error,
                min_density: step.density.iter().copied().fold(f64::INFINITY, f64::min),
                timings: step.timings,
            };
            let converged = converged(&record, &self.options);
            history.push(record);
            if !residual.is_finite()
                || !step.energy.total_hartree().is_finite()
                || !charge_error.is_finite()
                || !max_eigen_residual.is_finite()
            {
                return Err(ScfError::Diverged {
                    reason: "nonfinite SCF observable".into(),
                    history,
                });
            }
            if residual > self.options.divergence_factor * best.max(self.options.density_tolerance)
            {
                growing += 1;
            } else {
                growing = 0;
            }
            best = best.min(residual);
            if growing >= self.options.divergence_patience {
                return Err(ScfError::Diverged {
                    reason: "density residual exceeded growth limit".into(),
                    history,
                });
            }
            // Certify the returned output density itself, not just the input
            // fixed-point residual. A small-gap response can amplify that residual.
            let certified = if converged {
                let warm: Vec<_> = step
                    .k_points
                    .iter()
                    .map(|k| k.solution.eigenvectors.clone())
                    .collect();
                let check = self
                    .step_kpoints(
                        &step.density,
                        Some(&warm),
                        0.1 * self.options.eigen_tolerance,
                    )
                    .map_err(|e| match e {
                        ScfError::NonFinite => ScfError::Diverged {
                            reason: "nonfinite certification arithmetic".into(),
                            history: history.clone(),
                        },
                        other => other,
                    })?;
                let density_residual = density_distance(
                    &check.density,
                    &step.density,
                    self.cell.volume_bohr3(),
                    self.electrons,
                );
                let occupied_eigenvalue_change = check
                    .k_points
                    .iter()
                    .zip(&step.k_points)
                    .flat_map(|(a, b)| {
                        a.solution
                            .eigenvalues
                            .iter()
                            .zip(&b.solution.eigenvalues)
                            .zip(&self.occupations)
                    })
                    .filter(|(_, f)| **f > 0.0)
                    .map(|((a, b), _)| (a - b).abs())
                    .fold(0.0, f64::max);
                let passed = density_residual <= self.options.density_tolerance
                    && occupied_eigenvalue_change <= 1e-7;
                history.last_mut().unwrap().fixed_point_certificate = Some(FixedPointCertificate {
                    strict_residuals_per_k: check
                        .k_points
                        .iter()
                        .map(|k| {
                            k.solution
                                .residuals
                                .iter()
                                .zip(&k.solution.buffer)
                                .filter_map(|(&r, &buffer)| (!buffer).then_some(r))
                                .collect()
                        })
                        .collect(),
                    buffer_residuals_per_k: check
                        .k_points
                        .iter()
                        .map(|k| {
                            k.solution
                                .residuals
                                .iter()
                                .zip(&k.solution.buffer)
                                .filter_map(|(&r, &buffer)| buffer.then_some(r))
                                .collect()
                        })
                        .collect(),
                    density_residual,
                    occupied_eigenvalue_change,
                    davidson_iterations: check.k_points.iter().map(|k| k.solution.iterations).sum(),
                    operator_applications: check
                        .k_points
                        .iter()
                        .map(|k| k.solution.operator_applications)
                        .sum(),
                    davidson_iterations_per_k: check
                        .k_points
                        .iter()
                        .map(|k| k.solution.iterations)
                        .collect(),
                    operator_applications_per_k: check
                        .k_points
                        .iter()
                        .map(|k| k.solution.operator_applications)
                        .collect(),
                    timings: check.timings,
                });
                if passed {
                    self.check_gap(&check.k_points)?;
                }
                passed
            } else {
                false
            };
            if certified || (retain_last && iteration == self.options.max_iterations) {
                let gap_hartree = if certified {
                    self.check_gap(&step.k_points)?
                } else {
                    self.minimum_gap.map(|_| self.sampled_gap(&step.k_points))
                };
                return Ok(ScfOutcome {
                    converged: certified,
                    result: ScfResult {
                        k_points: step.k_points,
                        gap_hartree,
                        energy: step.energy,
                        eigenvalues: step.solution.eigenvalues,
                        occupations: self.occupations,
                        orbitals: step.solution.eigenvectors,
                        density: step.density,
                        charge,
                        history,
                        ewald_breakdown: self.ewald,
                        setup_time: self.setup_time,
                        total_time: start.elapsed() + self.setup_time,
                    },
                });
            }
            let mix_start = Instant::now();
            let mut mixed = mixer
                .mix(
                    n.as_slice().unwrap(),
                    step.density.as_slice().unwrap(),
                    self.electrons / self.cell.volume_bohr3(),
                )
                .map_err(|e| ScfError::Diverged {
                    reason: e.to_string(),
                    history: history.clone(),
                })?;
            history.last_mut().unwrap().mixing_history_used =
                Some(mixer.history_len().saturating_sub(1));
            history.last_mut().unwrap().mixing_raw_recovery = mixer.used_raw_recovery();
            // Anderson extrapolation is unconstrained. A convex step toward its
            // proposal preserves charge and Fourier support while keeping XC in
            // its density domain. Never clip individual spatial samples.
            if matches!(
                self.options.mixing,
                Mixing::Pulay { .. } | Mixing::KerkerPulay { .. }
            ) {
                let fraction = n
                    .iter()
                    .zip(&mixed)
                    .filter(|(_, b)| **b < 0.0)
                    .map(|(a, b)| 0.9 * a / (a - b))
                    .fold(1.0, f64::min);
                if fraction < 1.0 {
                    for (b, a) in mixed.iter_mut().zip(&n) {
                        *b = a + fraction * (*b - a);
                    }
                    mixer.reset();
                    let record = history.last_mut().unwrap();
                    record.mixing_fraction = fraction;
                    record.mixer_reset = true;
                    mixing_resets += 1;
                    record.mixing_resets = mixing_resets;
                }
            }
            n.as_slice_mut().unwrap().copy_from_slice(&mixed);
            history.last_mut().unwrap().timings.mixing = mix_start.elapsed();
            previous = Some(
                step.k_points
                    .into_iter()
                    .map(|k| k.solution.eigenvectors)
                    .collect(),
            );
        }
        Err(ScfError::NotConverged { history })
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
mod tests {
    use super::*;

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
            kinetic: 0.0,
            local: 0.0,
            local_zero_mode: 0.0,
            nonlocal: 0.0,
            hartree: 0.0,
            xc: 0.0,
            ewald: 0.0,
            band_sum: 0.0,
            double_counting: 0.0,
            band_form: 0.0,
            harris_foulkes: 0.0,
            functional_difference_bound: 0.0,
        };
        let good = IterationRecord {
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
                let error = ordered_solves(&mut visited, threads, |i, visited| {
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
                })
                .unwrap_err();
                assert!(error.to_string().contains("k index 2"), "{error}");
                assert!(visited.into_iter().all(|v| v));
            }
            let mut states = [0; 5];
            assert_eq!(
                ordered_solves(&mut states, threads, |i, state| {
                    *state = i;
                    Ok(i)
                })
                .unwrap(),
                vec![0, 1, 2, 3, 4]
            );
            assert_eq!(states, [0, 1, 2, 3, 4]);
        }
    }
}
