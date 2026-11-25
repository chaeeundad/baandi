use super::atomic_start;
pub use crate::refscf;
#[path = "start_density.rs"]
mod start_density;
pub use atomic_start::{AtomicStartReport, StartWavefunctions};
pub use dft_solvers::DavidsonPreconditioner;
pub use start_density::StartDensity;
#[path = "convergence.rs"]
pub(crate) mod convergence;
pub use crate::model::{StandardProfile, Verification, VerificationAttempt, VerificationReport};
#[path = "collinear.rs"]
pub mod collinear;
/// Compatibility exports for reference prolongation.
#[path = "cutoff_ladder.rs"]
pub mod cutoff_ladder;
#[path = "finite_temperature.rs"]
pub mod finite_temperature;
#[path = "session.rs"]
pub mod session;
use crate::Functional;
use crate::energy::{self, AuditInput, EnergyAudit};
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_geometry::{
    Atom, Cell,
    kpoints::{self, KPoint},
};
use dft_mpi::Collective;
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
use session::{CompletionMode, ScfSession, StepStatus};
use std::time::{Duration, Instant};

pub use crate::model::{PreparedSymmetry, ScfError, SpeciesPotential};

pub use crate::model::{CompatMode, Mixing, ScfOptions, bands_fit_subspace};

pub use crate::model::{FixedPointCertificate, IterationRecord, PhaseTimings};

pub use crate::model::OccupationPolicy;

pub(crate) use crate::model::result::materialize_kpoints;
pub use crate::model::{KPointSolution, ScfResult, split_warm};

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
    force_model: crate::force_identity::ModelDescription,
    /// RUN-CONFIG: the run's switch configuration (the calling thread's at
    /// construction; `set_run_config` replaces it) and the session switches
    /// resolved from it once.
    run: dft_numerics::switches::RunConfig,
    pub(crate) sw: session::SessionSwitches,
    atoms: Vec<Atom>,
    symmetry: Option<dft_physics::symmetry::DensitySymmetrizer>,
    /// `BAANDI_M2_GRID=on` (ticket M2-GRID, default off): slab plan of the
    /// distributed full-grid stages, built at the first eligible map.
    grid_dist: Option<crate::grid_dist::GridPlan>,
    /// M2-OWNED (`BAANDI_M2_LOCAL_FIELDS` / `_LOCAL_DENSITY`, default off):
    /// epoch-tagged owned fields of the M2 path.
    pub(crate) owned: crate::grid_dist::OwnedFields,
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

// BREAKTHROUGH-4 (`BAANDI_EWALD_CACHE=on`, default off): constructions on the
// same thread with bit-equal lattice, charges and fractional positions reuse
// the Ewald energy (the cutoff-ladder coarse problem after the full-cutoff
// one). The energy is a pure function of exactly these inputs and the default
// options, so a hit is bit-identical to recomputing it.
static EWALD_CACHE_COUNTS: [std::sync::atomic::AtomicU64; 2] =
    [const { std::sync::atomic::AtomicU64::new(0) }; 2];
thread_local! {
    static EWALD_CACHE_OVERRIDE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
    static EWALD_LAST: std::cell::RefCell<Option<(Vec<u64>, ewald::EwaldEnergy)>> =
        const { std::cell::RefCell::new(None) };
}
/// (hits, computed) Ewald energies of problem constructions since process
/// start (hits only with `BAANDI_EWALD_CACHE=on`).
pub fn ewald_cache_counters() -> [u64; 2] {
    use std::sync::atomic::Ordering::Relaxed;
    [
        EWALD_CACHE_COUNTS[0].load(Relaxed),
        EWALD_CACHE_COUNTS[1].load(Relaxed),
    ]
}
/// Override `BAANDI_EWALD_CACHE` for constructions on the calling thread
/// (tests); `None` restores the environment.
pub fn set_ewald_cache_override(on: Option<bool>) {
    EWALD_CACHE_OVERRIDE.with(|o| o.set(on));
}
/// Ionic local coefficients over the density sphere (INIT-SCALE-2, always on).
/// With a registered communicator of two or more ranks every rank evaluates
/// its contiguous block of G vectors with the per-G arithmetic of the full call
/// (`ionic_local_with_conventions_rows`) and one allgatherv assembles them:
/// bit-identical to the replicated evaluation for any rank count. When any
/// rank failed, every rank repeats the replicated call, which returns the
/// replicated error. Collective with a communicator.
fn ionic_local_coefficients(
    density: &PlaneWaveBasis,
    volume_bohr3: f64,
    species: &[Species<'_>],
    local_mesh: dft_upf::radial::LocalMesh,
    rule: dft_upf::radial::RadialEvenRule,
) -> Result<Vec<Complex64>, ScfError> {
    use dft_mpi::Collective;
    let full = || {
        hamiltonian::ionic_local_with_conventions(density, volume_bohr3, species, local_mesh, rule)
    };
    let Some(comm) = crate::parallel::SharedComm::active() else {
        return Ok(full()?);
    };
    let n = density.vectors().len();
    let counts = dft_mpi::block_counts(n, comm.size());
    let first: usize = counts[..comm.rank()].iter().sum();
    let local = hamiltonian::ionic_local_with_conventions_rows(
        density,
        volume_bohr3,
        species,
        local_mesh,
        rule,
        first..first + counts[comm.rank()],
    );
    // Per rank: 2 values per coefficient, then a failure flag.
    let mut send = vec![0.0; 2 * counts[comm.rank()] + 1];
    match &local {
        Ok(values) => {
            for (b, z) in send.as_chunks_mut::<2>().0.iter_mut().zip(values) {
                b[0] = z.re;
                b[1] = z.im;
            }
        }
        Err(_) => *send.last_mut().unwrap() = 1.0,
    }
    let gathered: Vec<usize> = counts.iter().map(|c| 2 * c + 1).collect();
    let mut all = vec![0.0; gathered.iter().sum()];
    comm.allgatherv_f64(&send, &gathered, &mut all)
        .map_err(|e| ScfError::Numerical(e.to_string()))?;
    let mut result = Vec::with_capacity(n);
    let mut at = 0;
    for &c in &counts {
        if all[at + 2 * c] != 0.0 {
            return Ok(full()?);
        }
        result.extend(
            all[at..at + 2 * c]
                .as_chunks::<2>()
                .0
                .iter()
                .map(|b| Complex64::new(b[0], b[1])),
        );
        at += 2 * c + 1;
    }
    crate::m2::count(crate::m2::Count::IonicLocalSplits);
    Ok(result)
}

/// Smallest communicator for which the Ewald sum is split over the ranks
/// (INIT-EWALD-GATE). dev08 (master e7c94c8, 10-step screens, max over ranks):
/// at np=2 the split was 83 to 102 % slower than the replicated sum (si32
/// 0.0595 -> 0.1089 s, si64 0.2451 -> 0.4949 s: per-rank candidate setup plus
/// the chained broadcasts outweigh halving the terms), at np=8 it was 48 %
/// faster. The result is bit-identical on both sides of the gate.
/// Re-measured on 6be8f09 (master, si64 np=2 `init.ewald`): split 0.4977 s, the
/// gated (replicated) path 0.8178 s, i.e. the gate made np=2 worse (si32 0.1084
/// vs 0.1813 s): the split is now the faster branch from 2 ranks on.
pub const EWALD_SPLIT_MIN_RANKS: usize = 2;
/// Whether `ewald_evaluate` splits the sum on a communicator of `ranks` ranks.
pub fn ewald_split_applies(ranks: usize) -> bool {
    ranks >= EWALD_SPLIT_MIN_RANKS
}
/// Ewald energy (INIT-SCALE-3, always on). With a registered communicator of
/// two or more ranks every rank evaluates the terms of its contiguous block of
/// the real-space and the reciprocal candidates (`ewald::EwaldPlan`, the
/// expensive erfc and structure-factor work; INIT-EWALD-FAST: a real-space
/// candidate is an ion pair and enumeration row in the reference scheme, an
/// ion pair i <= j in the balanced scheme), then the two running sums pass
/// down the ranks: rank r adds its terms in order to the sums received from
/// rank r - 1 and broadcasts them (P broadcasts of three values). This is the
/// serial loop's association exactly, so the result is bit-identical to
/// `ewald::energy` for any rank count. A failure on any rank (coincident ions)
/// is carried in the third value and makes every rank repeat the serial call,
/// which returns the serial error. Collective with a communicator.
///
/// The split runs only for communicators of at least
/// [`EWALD_SPLIT_MIN_RANKS`] ranks; smaller ones evaluate the serial sum on
/// every rank (bit-identical either way, no collective).
fn ewald_evaluate(cell: &Cell, ions: &[Ion]) -> Result<ewald::EwaldEnergy, ScfError> {
    use dft_mpi::Collective;
    let serial = || Ok(ewald::energy(cell, ions, EwaldOptions::default())?);
    let Some(comm) = crate::parallel::SharedComm::active() else {
        return serial();
    };
    if !ewald_split_applies(comm.size()) {
        return serial();
    }
    // The plan depends only on replicated inputs: every rank succeeds or fails alike.
    let Ok(plan) = ewald::EwaldPlan::new(cell, ions, EwaldOptions::default()) else {
        return serial();
    };
    let (rank, ranks) = (comm.rank(), comm.size());
    let block = |n: usize| {
        let counts = dft_mpi::block_counts(n, ranks);
        let first: usize = counts[..rank].iter().sum();
        first..first + counts[rank]
    };
    let (mut real, mut reciprocal) = (Vec::new(), Vec::new());
    let failed = plan
        .real_terms(block(plan.real_candidates()), &mut real)
        .is_err();
    if !failed {
        plan.reciprocal_terms(block(plan.reciprocal_candidates()), &mut reciprocal);
    }
    // [real sum, reciprocal sum, failure flag], passed from rank 0 to P - 1.
    let mut state = [0.0; 3];
    for root in 0..ranks {
        if root == rank {
            if failed {
                state[2] = 1.0;
            } else if state[2] == 0.0 {
                state[0] = ewald::accumulate(state[0], &real);
                state[1] = ewald::accumulate(state[1], &reciprocal);
            }
        }
        comm.bcast_f64(&mut state, root)
            .map_err(|e| ScfError::Numerical(e.to_string()))?;
    }
    if state[2] != 0.0 {
        return serial();
    }
    crate::m2::count(crate::m2::Count::EwaldSplits);
    Ok(plan.finish_sums(state[0], state[1])?)
}

fn ewald_energy(
    cell: &Cell,
    ions: &[Ion],
    run: dft_numerics::switches::RunConfig,
) -> Result<ewald::EwaldEnergy, ScfError> {
    use dft_numerics::switches::RunCounter;
    use std::sync::atomic::Ordering::Relaxed;
    let on = EWALD_CACHE_OVERRIDE
        .with(|o| o.get())
        .unwrap_or_else(|| run.on(dft_numerics::switches::Id::EWALD_CACHE));
    if !on {
        EWALD_CACHE_COUNTS[1].fetch_add(1, Relaxed);
        run.count(RunCounter::EwaldComputed, 1);
        return ewald_evaluate(cell, ions);
    }
    let mut key = Vec::with_capacity(10 + 4 * ions.len());
    // The scheme (BAANDI_EWALD_FAST or a per-thread override) selects a
    // different pure function: part of the key.
    key.push(u64::from(ewald::scheme() == ewald::Scheme::Balanced));
    key.extend(
        cell.lattice_columns_bohr()
            .iter()
            .flatten()
            .map(|v| v.to_bits()),
    );
    for ion in ions {
        key.push(ion.charge.to_bits());
        key.extend(ion.fractional_position.iter().map(|v| v.to_bits()));
    }
    let hit = EWALD_LAST.with(|l| {
        l.borrow()
            .as_ref()
            .filter(|(k, _)| *k == key)
            .map(|(_, e)| e.clone())
    });
    if let Some(e) = hit {
        EWALD_CACHE_COUNTS[0].fetch_add(1, Relaxed);
        run.count(RunCounter::EwaldHits, 1);
        return Ok(e);
    }
    EWALD_CACHE_COUNTS[1].fetch_add(1, Relaxed);
    run.count(RunCounter::EwaldComputed, 1);
    let e = ewald_evaluate(cell, ions)?;
    EWALD_LAST.with(|l| *l.borrow_mut() = Some((key, e.clone())));
    Ok(e)
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

pub use crate::model::ScfOutcome;

pub use crate::model::ScfStep;
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
        // RUN-CONFIG: construction resolves its switches from the calling
        // thread's configuration (scope or process) once.
        if options.standard.is_some() || options.verification != crate::Verification::Off {
            return Err(ScfError::InvalidInput(
                "the strict oracle requires ScfOptions::strict_reference()",
            ));
        }
        let run = dft_numerics::switches::RunConfig::current();
        run.count(dft_numerics::switches::RunCounter::Problems, 1);
        use dft_numerics::switches::{InitProfile, RunCounter as IC};
        let init_profile = InitProfile::current();
        let _init_total = init_profile.map(|p| p.span(IC::InitProblemCalls, IC::InitProblemNanos));
        if let Some(p) = init_profile {
            p.count(IC::InitAtoms, atoms.len());
            p.count(IC::InitKpoints, points.len());
        }
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

        if options.locking == dft_solvers::locking::Locking::Auto
            && !matches!(options.eigensolver, dft_solvers::SolverChoice::Davidson)
        {
            return Err(ScfError::InvalidInput("locking requires strict Davidson"));
        }
        // BREAKTHROUGH-4: construction parts for the phase profile (no clock is
        // read when the profile is off).
        let part = crate::profile::clock();
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
        crate::profile::init_part(crate::profile::InitPart::Rest, part);
        let part = crate::profile::clock();
        let gamma = options.basis(cell, [0.; 3], 1.0)?;
        let density_basis = options.basis(cell, [0.; 3], 4.0)?;

        let grid = FftGridSpec::new(&gamma, &density_basis, options.grid_override)?;
        let dg = FftGridSpec::new(&density_basis, &density_basis, Some(grid.dims()))?;
        if let Some(p) = init_profile {
            p.count(IC::InitBands, bands);
            p.count(IC::InitGridPoints, grid.dims().iter().product());
            p.count(IC::InitWaveRows, gamma.vectors().len());
            p.count(IC::InitDensityRows, density_basis.vectors().len());
        }
        crate::profile::init_part(crate::profile::InitPart::Basis, part);
        let part = crate::profile::clock();
        let species: Vec<_> = potentials
            .iter()
            .zip(&positions)
            .map(|(p, r)| Species {
                potential: p.potential,
                positions_bohr: r,
            })
            .collect();
        let ionic_g = dg.scatter(&ionic_local_coefficients(
            &density_basis,
            cell.volume_bohr3(),
            &species,
            options.local_mesh,
            options.radial_even_rule,
        )?)?;
        crate::profile::init_part(crate::profile::InitPart::IonicLocal, part);
        let part = crate::profile::clock();
        let mut fft = FftPlan3d::new(grid.dims())?;
        let span = init_profile.map(|p| p.span(IC::InitCoreDensityCalls, IC::InitCoreDensityNanos));
        let core_density = dft_physics::data::core_density_with_rule(
            &density_basis,
            &dg,
            cell.volume_bohr3(),
            &species,
            &mut fft,
            options.radial_even_rule,
        )?;
        drop(span);
        let span =
            init_profile.map(|p| p.span(IC::InitStartDensityCalls, IC::InitStartDensityNanos));
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
        drop(span);
        let xc_density = Array3::zeros(
            if core_density.is_some() || options.pw92_prefactor == crate::Pw92Prefactor::Original {
                grid.dims()
            } else {
                [0; 3]
            },
        );
        let span = init_profile.map(|p| p.span(IC::InitIonicFftCalls, IC::InitIonicFftNanos));
        let ionic_real = if let Some(comm) = crate::m2::active() {
            crate::m2::initial_ionic_real(&comm, &mut fft, &ionic_g)?
        } else {
            let mut ionic = Array3::zeros(grid.dims());
            fft.to_real(&ionic_g, &mut ionic)?;
            let scale = ionic.iter().map(|z| z.re.abs()).fold(1.0, f64::max);
            if ionic.iter().any(|z| z.im.abs() > 1e-12 * scale) {
                return Err(ScfError::Numerical("ionic imaginary potential".into()));
            }
            ionic.mapv(|z| z.re)
        };
        drop(span);
        crate::profile::init_part(crate::profile::InitPart::Densities, part);
        let part = crate::profile::clock();
        let pool =
            WorkerPool::new(options.threads).map_err(|e| ScfError::Numerical(e.to_string()))?;
        // Every k solve shares the persistent pool. A full-pool budget lets
        // idle outer workers steal fixed inner pieces and gives each worker
        // its own Hamiltonian scratch slot (no modulo-slot collisions).
        let inner = options.threads.max(1);
        let mut construction_points = points.to_vec();
        // Capture MPI identity on the caller before entering pool workers;
        // the communicator registry itself is thread-local.
        let owned_layout = crate::m2::active()
            .filter(|_| {
                gamma_real_exclusion.is_none() && options.gamma_real != crate::GammaReal::Off
            })
            .map(|comm| (comm.rank(), comm.size()));
        let states = ordered_solves(
            &mut construction_points,
            options.threads,
            &pool,
            |_, point| {
                // INIT-REUSE: the exact Gamma point reuses the Gamma basis built
                // above (the same call with the same arguments).
                let basis = if point.fractional.map(f64::to_bits) == [0; 3] {
                    gamma.clone()
                } else {
                    options.basis(cell, point.fractional, 1.0)?
                };
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
                let span =
                    init_profile.map(|p| p.span(IC::InitStateGridCalls, IC::InitStateGridNanos));
                let k_grid = FftGridSpec::new(&basis, &density_basis, Some(grid.dims()))?;
                drop(span);
                // INIT-REUSE: the k grid is the Hamiltonian's grid (same basis,
                // same density sphere vectors, same dims); passing it avoids
                // building the density sphere again per k state.
                let nonlocal = if let Some((rank, ranks)) =
                    owned_layout.filter(|_| point.fractional == [0.; 3])
                {
                    let indices = dft_physics::distributed_real_gamma::owned_projector_indices(
                        &basis, &k_grid, rank, ranks,
                    )
                    .map_err(|e| ScfError::Numerical(e.to_string()))?;
                    // PHASE_VALIDATE instruments only this owned constructor;
                    // the serial/eager branch below performs no such scan.
                    // Validation stays before phase materialization/collectives.
                    NonlocalOperator::new_owned_with_rule(
                        &basis,
                        cell.volume_bohr3(),
                        &species,
                        options.radial_even_rule,
                        &indices,
                    )?
                } else {
                    NonlocalOperator::new_with_rule(
                        &basis,
                        cell.volume_bohr3(),
                        &species,
                        options.radial_even_rule,
                    )?
                };
                let span = init_profile
                    .map(|p| p.span(IC::InitHamiltonianCalls, IC::InitHamiltonianNanos));
                let hamiltonian = Hamiltonian::new_with_grid(
                    &basis,
                    k_grid.clone(),
                    ionic_real.clone(),
                    nonlocal,
                )?
                .with_band_pool(inner, pool.clone())?;
                drop(span);
                let span = init_profile
                    .map(|p| p.span(IC::InitBandDensityCalls, IC::InitBandDensityNanos));
                let band_density =
                    crate::band_density::BandDensity::new(&k_grid, bands, inner, pool.clone())?;
                drop(span);
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
        crate::profile::init_part(crate::profile::InitPart::States, part);
        let part = crate::profile::clock();
        let mut atomic_start = AtomicStartReport::default();
        let atomic_species = if options.start_wavefunctions == StartWavefunctions::Atomic {
            // ONEPATH-3B: spin channels (collinear fixed sector) prepare the
            // atomic start too; the Fermi-Dirac constructors (scalar and
            // collinear) clear it with their own reason.
            let reason = if species.iter().any(|s| {
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
        crate::profile::init_part(crate::profile::InitPart::AtomicPrepare, part);
        let part = crate::profile::clock();
        let span = init_profile.map(|p| p.span(IC::InitWorkerFftsCalls, IC::InitWorkerFftsNanos));
        let worker_ffts = (0..options.threads.max(1).min(states.len()))
            .map(|_| FftPlan3d::new(grid.dims()))
            .collect::<Result<Vec<_>, _>>()?;
        drop(span);
        crate::profile::init_part(crate::profile::InitPart::Rest, part);
        let part = crate::profile::clock();
        let ewald = ewald_energy(cell, &ions, run)?;
        crate::profile::init_part(crate::profile::InitPart::Ewald, part);
        let part = crate::profile::clock();
        let occupations = (0..bands)
            .map(|b| {
                if b < occupied {
                    if channel.is_some() { 1.0 } else { 2.0 }
                } else {
                    0.0
                }
            })
            .collect();
        let force_model = crate::force_identity::ModelDescription::new(
            cell,
            atoms,
            potentials,
            &options,
            grid.dims(),
            states.iter().map(|s| s.point).collect(),
        );
        crate::profile::init_part(crate::profile::InitPart::Rest, part);
        Ok(Self {
            #[cfg(test)]
            force_real_rebuild: false,
            #[cfg(test)]
            force_energy_identity: false,
            gamma_real_exclusion,
            reused_fields: None,
            screen_ffts: Vec::new(),
            energy_identity: false,
            start_density,
            atomic_species,
            atomic_start,
            force_model,
            run,
            sw: session::SessionSwitches::resolve(run),
            grid_dist: None,
            owned: Default::default(),
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
    /// Drivers whose start does not use the atomic subspace (Fermi-Dirac
    /// occupations) drop the prepared species and record the reason
    /// (ADR 0042; ONEPATH-3B moved this out of the spin-channel construction).
    pub(crate) fn clear_atomic_start(&mut self, reason: &str) {
        if self.atomic_species.take().is_some() {
            self.atomic_start.fallback_reason = Some(reason.to_owned());
        }
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
    /// The real-space ionic local potential built at construction (Hartree,
    /// FFT grid). Read-only; identity tests of the distributed construction.
    #[doc(hidden)]
    pub fn ionic_potential(&self) -> &Array3<f64> {
        &self.ionic_real
    }
    /// The Ewald energy terms computed at construction. Read-only; identity
    /// tests of the distributed construction.
    #[doc(hidden)]
    pub fn ewald_terms(&self) -> &ewald::EwaldEnergy {
        &self.ewald
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

    /// The pointwise LDA flavour when the distributed grid stages apply to this
    /// problem (`BAANDI_M2_GRID=on`, M2 path, LDA PZ/PW with the modified PW92
    /// prefactor, no core density, density floor or symmetry); replicated.
    fn grid_dist_lda(&self) -> Option<crate::grid_dist::Lda> {
        if !crate::grid_dist::requested()
            || self.core_density.is_some()
            || self.density_floor.is_some()
            || self.symmetry.is_some()
            || self.options.pw92_prefactor != crate::Pw92Prefactor::Modified
        {
            return None;
        }
        match self.options.functional {
            Functional::LdaPz => Some(crate::grid_dist::Lda::Pz),
            Functional::LdaPw => Some(crate::grid_dist::Lda::Pw),
            _ => None,
        }
    }

    /// Fields of an SCF map's input density: the distributed grid stages when
    /// eligible (`true`: V_H and V_xc stay on the plan's planes, the returned
    /// grids carry only the energies and empty potential arrays), otherwise
    /// `input_fields`. Collective on the M2 path.
    /// M2-OWNED (`BAANDI_M2_LOCAL_FIELDS=on`): whether this map's effective
    /// potential may stay on the owned planes (replicated conditions only: the
    /// distributed grid path, the standard profile, one exact-Gamma state whose
    /// previous map used the distributed operator, a warm start so that no
    /// atomic start applies the complex Hamiltonian, no snapshot capture, no
    /// finite temperature). Counts the reason of a decline.
    fn local_fields_eligible(&self, _warm: bool, _finite: bool) -> bool {
        use crate::grid_dist as g;
        if !g::local_fields_requested() || crate::m2::active().is_none() {
            return false;
        }
        let reason = if let Some(d) = self.owned.declined {
            Some(d.index())
        } else if self.grid_dist_lda().is_none() {
            Some(self.grid_path_reason())
        } else {
            Some(1)
        };
        if let Some(r) = reason {
            g::owned_declined(r);
            return false;
        }
        true
    }

    /// Strict maps always keep the full complex potential current.
    pub(crate) fn ensure_complex_potential(&mut self) -> Result<(), ScfError> {
        Ok(())
    }

    /// M2-OWNED (`BAANDI_M2_LOCAL_DENSITY`): static eligibility of the owned
    /// density loop (replicated): both owned switches and `BAANDI_M2_MIX`, the
    /// M2 path, the distributed grid path, the standard profile, one k state,
    /// plain Pulay without a Hartree metric, no symmetry, no snapshot capture,
    /// no replay capture. A requested but ineligible loop is counted.
    pub(crate) fn record_owned_density_decline(&self) {
        use crate::grid_dist as g;
        if !g::local_density_requested() {
            return;
        }
        let reason = if let Some(d) = self.owned.declined {
            Some(d.index())
        } else if crate::m2::active().is_none() || !crate::mix_dist::requested() {
            Some(g::R_GRID_PATH)
        } else if self.grid_dist_lda().is_none() {
            Some(self.grid_path_reason())
        } else {
            Some(1)
        };
        if let Some(r) = reason {
            g::density_declined(r);
        }
    }

    /// Why `grid_dist_lda` is None: `functional` for a non-LDA functional
    /// (PBE and other GGAs) with `BAANDI_M2_GRID` on, else `grid_path`
    /// (switch off, core density, density floor, symmetry, PW92 prefactor).
    fn grid_path_reason(&self) -> usize {
        use crate::grid_dist as g;
        if g::requested()
            && !matches!(
                self.options.functional,
                Functional::LdaPz | Functional::LdaPw
            )
        {
            g::R_FUNCTIONAL
        } else {
            g::R_GRID_PATH
        }
    }

    /// M2-OWNED: declare a consumer of active-boundary state before the
    /// session is initialized (`grid_dist::OwnedDecline`). Both owned modes
    /// then take the gathered path for this problem with the reason counted;
    /// `Checkpoint` also selects the replicated Pulay mixer (whose state an
    /// exact checkpoint carries). Every rank must declare the same reason.
    /// The first declaration is kept.
    pub fn decline_owned(&mut self, reason: crate::grid_dist::OwnedDecline) {
        if self.owned.declined.is_none() {
            self.owned.declined = Some(reason);
        }
    }
    /// The declared owned-mode decline, if any.
    pub fn owned_declined(&self) -> Option<crate::grid_dist::OwnedDecline> {
        self.owned.declined
    }

    /// Session metrics of a map: (residual `density_distance(output, input)`,
    /// output sum, input sum, output minimum); the full-array arithmetic of
    /// the session, or `owned_metrics` for an owned pair. Collective when owned.
    pub(crate) fn step_metrics(
        &mut self,
        input: &Array3<f64>,
        output: &Array3<f64>,
    ) -> Result<(f64, f64, f64, f64), ScfError> {
        Ok((
            density_distance(output, input, self.cell.volume_bohr3(), self.electrons),
            output.sum(),
            input.sum(),
            output.iter().copied().fold(f64::INFINITY, f64::min),
        ))
    }

    fn map_fields(
        &mut self,
        n: &Array3<f64>,
    ) -> Result<(hartree::HartreeGrid, xc::XcGrid, Array3<f64>, bool), ScfError> {
        let comm = crate::m2::active();
        let lda = comm.as_ref().and(self.grid_dist_lda());
        let (Some(comm), Some(lda), Some(values), Some(ionic)) =
            (comm, lda, n.as_slice(), self.ionic_real.as_slice())
        else {
            if crate::grid_dist::requested() && crate::m2::active().is_some() {
                crate::grid_dist::count_declined();
            }
            let (h, x, effective) = self.input_fields(n)?;
            return Ok((h, x, effective, false));
        };
        let field = crate::grid_dist::Field::Full(values);
        let dims = self.grid.dims();
        if self.grid_dist.as_ref().map(|p| p.dims()) != Some(dims) {
            let weights = None;
            self.grid_dist = Some(crate::grid_dist::GridPlan::new(
                &comm,
                dims,
                self.cell.reciprocal_columns_inverse_bohr(),
                self.cell.volume_bohr3(),
                lda,
                weights,
            )?);
        }
        let plan = self.grid_dist.as_mut().expect("plan built above");
        let (effective, (hartree_energy, xc_energy, xc_integral)) = {
            let mut effective = Array3::zeros(dims);
            let energies = plan.potential_field(
                &comm,
                field,
                ionic,
                Some(effective.as_slice_mut().expect("standard layout")),
            )?;
            (effective, energies)
        };

        let empty = (0, 0, 0);
        Ok((
            hartree::HartreeGrid {
                potential: Array3::zeros(empty),
                reciprocal_potential: Array3::zeros(empty),
                energy: Hartree(hartree_energy),
            },
            xc::XcGrid {
                potential: Array3::zeros(empty),
                energy: Hartree(xc_energy),
                potential_integral: Hartree(xc_integral),
            },
            effective,
            true,
        ))
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
        let start = dft_numerics::stage_clock::start();
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
        let mut local = vec![0.0; points];
        let keep: &mut [Complex64] = &mut [];
        {
            md.contributions.clear();
            md.contributions.resize(bands * points, 0.);
            dist.band_densities_keeping(
                &comm,
                &xm,
                &occupations[..m],
                volume,
                &mut md.contributions[..m * points],
                keep,
            )
            .map_err(err)?;
            // Ascending band order on the planes, as the replicated fold.
            for b in 0..m {
                for (l, c) in local
                    .iter_mut()
                    .zip(&md.contributions[b * points..(b + 1) * points])
                {
                    *l += c;
                }
            }
        }
        let output = {
            let mut full = Array3::zeros(dims);
            let gather_clock = dft_numerics::stage_clock::start();
            dist.gather_planes(&comm, &local, full.as_slice_mut().unwrap())
                .map_err(err)?;
            {
                use crate::grid_dist as g;
                g::owned_add(g::O_DENSITY_GATHERS, 1);
                g::owned_add(g::O_DENSITY_GATHER_VALUES, full.len() as u64);
                g::owned_time(g::T_DENSITY_GATHER, gather_clock);
            }
            if full.iter().any(|n: &f64| !n.is_finite()) {
                return Err(ScfError::NonFinite);
            }
            let mut output = Array3::zeros(dims);
            for (out, n) in output.iter_mut().zip(&full) {
                *out += weight * n;
            }
            output
        };
        // Without per-band contributions the screen data is not held.
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
        let start = dft_numerics::stage_clock::start();
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
    /// The Hartree/XC fields built here are kept for the next
    /// map (`step_*`) at the same density (POTENTIAL-REUSE), so a certificate
    /// solve that follows a passing screen does not build them again.
    pub(crate) fn certificate_shift_estimate(
        &mut self,
        density: &Array3<f64>,
        k_points: &[KPointSolution],
    ) -> Result<Option<f64>, ScfError> {
        const CLUSTER_GAP_HA: f64 = 1e-3;
        const MAX_CLUSTER: usize = 24;
        self.reused_fields = None;
        // M2-OWNED: the screen reads the complex Hamiltonian's full potential.
        self.ensure_complex_potential()?;
        if k_points.len() != self.states.len() {
            return Ok(None);
        }
        let (h, x, effective) = self.input_fields(density)?;
        let standard = effective.as_standard_layout().into_owned();
        {
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
        let m2_comm = crate::m2::active();
        let m2_on = m2_comm.is_some();
        let screen_comm = m2_comm.filter(|_| crate::m2::tuning().screen_dist);
        // Plane-local delta, band values and entry sums of the distributed screen.
        let (mut delta_local, mut values, mut partial, mut reduced) =
            (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for (state, k) in self.states.iter_mut().zip(k_points) {
            let current = state.hamiltonian.local_potential_values();
            let eigenvalues = &k.solution.eigenvalues;
            if current.len() != points || eigenvalues.len() != bands {
                return Ok(None);
            }
            for ((d, t), c) in delta.iter_mut().zip(target).zip(current) {
                *d = t - c;
            }
            delta_local.clear();
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
                entries.clear();
                for j in 0..size {
                    for i in 0..=j {
                        entries.push((i, j, Complex64::default()));
                    }
                }
                // BAANDI_M2_SCREEN=dist: the band values on this rank's FFT planes
                // (distributed transform), plane-local sums and one tree_sum per
                // cluster instead of full-grid transforms on every rank. All
                // conditions are replicated, so every rank takes the same branch.
                let cluster_clock = dft_numerics::stage_clock::start();
                if let Some(comm) = screen_comm.as_ref()
                    && let Some(rs) = state.real_solve.as_mut()
                    && let (Some(dist), Some(x), Some(md)) = (
                        rs.dist.as_mut(),
                        rs.m2_local.as_ref(),
                        rs.m2_density.as_ref().filter(|d| d.valid),
                    )
                    && x.cols() == bands
                    && md.points == dist.plane_points()
                    && md.offsets.len() == md.points
                {
                    let pts = md.points;
                    if delta_local.len() != pts {
                        delta_local.clear();
                        delta_local.extend(md.offsets.iter().map(|&o| delta[o]));
                    }
                    let n_loc = x.rows();
                    let mut xm = dft_numerics::RealMatrix::zeros(n_loc, size)
                        .map_err(dft_solvers::SolverError::from)?;
                    xm.as_mut_slice()
                        .copy_from_slice(&x.as_slice()[start * n_loc..end * n_loc]);
                    values.clear();
                    values.resize(size * pts, 0.);
                    dist.band_values(comm, &xm, &mut values)
                        .map_err(|e| ScfError::Numerical(format!("M2 screen: {e}")))?;
                    partial.clear();
                    for &(i, j, _) in &entries {
                        let (a, b) = (
                            &values[i * pts..(i + 1) * pts],
                            &values[j * pts..(j + 1) * pts],
                        );
                        let mut sum = 0.0;
                        for ((a, b), d) in a.iter().zip(b).zip(&delta_local) {
                            sum += a * b * d;
                        }
                        partial.push(sum);
                    }
                    reduced.clear();
                    reduced.resize(partial.len(), 0.);
                    dft_mpi::Collective::tree_sum_f64(comm, &partial, &mut reduced)
                        .map_err(|e| ScfError::Numerical(e.to_string()))?;
                    for (e, v) in entries.iter_mut().zip(&reduced) {
                        e.2 = Complex64::new(*v, 0.0);
                    }
                    crate::m2::count(crate::m2::Count::ScreenDistClusters);
                } else {
                    psi.resize_with(size, || Array3::zeros(state.grid.dims()));
                    // One band transform per pool worker (own FFT plan each) and the
                    // cluster entries split over the pool; every transform and every
                    // entry keeps its serial arithmetic, so the estimate is
                    // bit-identical for any thread count.
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
                    let slices: Vec<&[Complex64]> = psi
                        .iter()
                        .map(|p| p.as_slice().expect("standard layout"))
                        .collect();
                    let (slices, delta) = (slices.as_slice(), delta.as_slice());
                    let entries_slice = entries.as_mut_slice();
                    self.pool
                        .install(move || cluster_entries(slices, delta, entries_slice));
                }
                if m2_on {
                    crate::m2::add(crate::m2::Stage::ScreenClusters, cluster_clock);
                    crate::m2::count_n(crate::m2::Count::ScreenClusterBands, size);
                }
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
        crate::m2::record_screen_bound(bound);
        if bound.is_finite() {
            Ok(Some(bound))
        } else {
            Ok(None)
        }
    }

    /// Atomic-start vectors of a cold map (one per k state) for the local
    /// potential currently set on the states; None without atomic species.
    /// The single implementation used by every cold map.
    fn atomic_guesses(&mut self) -> Result<Option<Vec<ComplexMatrix>>, ScfError> {
        if self.atomic_species.is_none() {
            return Ok(None);
        }
        // M2-OWNED: the atomic start applies the complex Hamiltonian.
        self.ensure_complex_potential()?;
        let Some(species) = &self.atomic_species else {
            return Ok(None);
        };
        let atomic_clock = Instant::now();
        // BAANDI_M2_ATOMIC=split: span, Gram-Schmidt and H Q split over ranks.
        let split = crate::m2::active().filter(|_| crate::m2::tuning().atomic_split);
        let mut guesses = Vec::with_capacity(self.states.len());
        for state in &mut self.states {
            let part = crate::profile::clock();
            let span = atomic_start::span(
                &self.cell,
                &state.basis,
                &state.grid,
                species,
                &self.options,
                self.occupations.len(),
                Some(&*self.pool),
                split.as_ref(),
            )?;
            crate::profile::atomic_part(crate::profile::AtomicPart::Span, part);
            let part = crate::profile::clock();
            // Owned init prepares the small D B^H Q once. The existing split
            // Ritz still applies H and forms Gram columns on their owners.
            let atomic_weights = if state.hamiltonian.nonlocal().has_owned_rows() {
                match crate::m2::active() {
                    Some(comm) => {
                        let weights = crate::m2::owned_atomic_weights(
                            &comm,
                            &state.hamiltonian,
                            span.vectors(),
                        )?;
                        Some(if let Some(split) = split.as_ref() {
                            let counts = dft_mpi::block_counts(weights.cols(), split.size());
                            let starts = dft_mpi::block_displacements(&counts);
                            let lo = starts[split.rank()] * weights.rows();
                            let len = counts[split.rank()] * weights.rows();
                            ComplexMatrix::new(
                                weights.rows(),
                                counts[split.rank()],
                                weights.as_slice()[lo..lo + len].to_vec(),
                            )
                            .map_err(dft_solvers::SolverError::from)?
                        } else {
                            weights
                        })
                    }
                    None => None,
                }
            } else {
                None
            };
            // ATOMIC-REAL remains the unchanged serial alternative.
            let hq = if atomic_weights.is_none() {
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
            crate::profile::atomic_part(crate::profile::AtomicPart::Hq, part);
            let part = crate::profile::clock();
            let ritz = || {
                atomic_start::ritz(
                    &span,
                    &state.hamiltonian,
                    self.occupations.len(),
                    state.point.fractional == [0.; 3],
                    hq,
                    split.as_ref(),
                )
            };
            let (guess, report) = match atomic_weights {
                Some(weights) => state.hamiltonian.with_atomic_weights(weights, ritz)?,
                None => ritz(),
            }?;
            crate::profile::atomic_part(crate::profile::AtomicPart::Ritz, part);
            self.atomic_start.atomic_orbitals += report.atomic_orbitals;
            self.atomic_start.initial_h_applications += report.initial_h_applications;
            self.atomic_start.gram_rotation_macs += report.gram_rotation_macs;
            self.atomic_start.mgs_projection_macs += report.mgs_projection_macs;
            self.atomic_start.radial_samples += report.radial_samples;
            self.atomic_start.evd_orders.extend(report.evd_orders);
            guesses.push(guess);
        }
        if split.is_some() {
            crate::m2::count(crate::m2::Count::AtomicSplitStarts);
        }
        if crate::m2::active().is_some() {
            crate::m2::add(
                crate::m2::Stage::AtomicStart,
                dft_numerics::stage_clock::keep(atomic_clock),
            );
        }
        crate::phase_clock::add_atomic_start(atomic_clock);
        crate::profile::add_atomic_start(atomic_clock.elapsed().as_secs_f64());
        Ok(Some(guesses))
    }

    /// Measurement access (stage-0 benchmark `examples/davidson_bench.rs` of
    /// dft-cli): sets the local potential of density `n` exactly as a map does
    /// and returns the atomic-start vectors a cold map at `n` would pass to the
    /// eigensolver (None when the problem has no atomic start). A map at `n`
    /// given these vectors as `initial` runs the identical solve without
    /// recomputing the start. Not used by the SCF.
    #[doc(hidden)]
    pub fn atomic_start_vectors(
        &mut self,
        n: &Array3<f64>,
    ) -> Result<Option<Vec<ComplexMatrix>>, ScfError> {
        if n.shape() != self.grid.dims() {
            return Err(ScfError::InvalidInput("input density grid"));
        }
        if self.options.start_wavefunctions != StartWavefunctions::Atomic {
            return Ok(None);
        }
        let (_, _, effective) = self.input_fields(n)?;
        ordered_solves(
            &mut self.states,
            self.options.threads,
            &self.pool,
            |_, state| {
                state.hamiltonian.set_local_potential(effective.clone())?;
                Ok(())
            },
        )?;
        // M2-OWNED: the complex Hamiltonians now hold the potential of `n`;
        // an earlier owned map's planes must not overwrite it.
        self.atomic_guesses()
    }

    fn step_with_occupations(
        &mut self,
        n: &Array3<f64>,
        initial: Option<&[ComplexMatrix]>,
        tolerance: f64,
        certify: bool,
        finite: Option<&mut finite_temperature::OccupationWorkspace>,
    ) -> Result<ScfStep, ScfError> {
        self.step_map(n, initial, &[], tolerance, certify, finite)
    }

    /// Fixed-potential map accepting materializable warm-start blocks.
    pub(crate) fn step_session(
        &mut self,
        n: &Array3<f64>,
        initial: Option<&[ComplexMatrix]>,
        initial_lazy: &[Option<crate::gamma_real::LazyOrbitals>],
        tolerance: f64,
        certify: bool,
    ) -> Result<ScfStep, ScfError> {
        self.step_map(n, initial, initial_lazy, tolerance, certify, None)
    }

    /// Compatibility map entry; delivered orbitals are always materialized.
    #[doc(hidden)]
    pub fn step_kpoints_lazy(
        &mut self,
        n: &Array3<f64>,
        initial: Option<&[ComplexMatrix]>,
        initial_lazy: &[Option<crate::gamma_real::LazyOrbitals>],
        tolerance: f64,
        certify: bool,
    ) -> Result<ScfStep, ScfError> {
        self.step_session(n, initial, initial_lazy, tolerance, certify)
    }

    /// Gather every lazy block of `k_points` into its full complex form
    /// (collective on the M2 path: every rank calls it at the same point) and
    /// drop the handles. No-op for full blocks.
    #[doc(hidden)]
    pub fn materialize_orbitals(&self, k_points: &mut [KPointSolution]) -> Result<(), ScfError> {
        materialize_kpoints(k_points)
    }

    /// Snapshot replay (BREAKTHROUGH-0-INSTRUMENT, `davidson_bench --snapshot`;
    /// measurement only, never used by the SCF): one session map at density
    /// `n` whose real Gamma solve starts from `start` (this rank's rows on the
    /// M2 path, the full ORC block on the threads path) with its kept H X
    /// armed when given, exactly as the captured solve did. One k point only.
    #[doc(hidden)]
    pub fn replay_solve(
        &mut self,
        n: &Array3<f64>,
        start: crate::snapshot::ReplayStart,
        tolerance: f64,
    ) -> Result<ScfStep, ScfError> {
        if self.states.len() != 1 {
            return Err(ScfError::InvalidInput("replay needs one k point"));
        }
        eprintln!(
            "snapshot replay: legacy-unchecked (no operator/row identity check); use replay_checked"
        );
        let bands = self.occupations.len();
        crate::snapshot::set_contract(None);
        crate::snapshot::inject(crate::profile::rank_key(), start);
        let placeholder =
            vec![ComplexMatrix::zeros(0, bands).map_err(dft_solvers::SolverError::from)?];
        self.step_session(n, Some(&placeholder), &[None], tolerance, false)
    }

    /// Checked snapshot replay (ISSUE-18, measurement only). Before anything
    /// is injected, the current operator's fingerprint (`operator_fingerprint`
    /// of `n`) must equal the captured one and the density bytes must equal
    /// the captured density (exact and A/B modes); the result is agreed over
    /// all ranks (one tree sum) so that every rank returns the same error
    /// before any solve. The solver entry then checks this rank's rows,
    /// dimensions, potential, stopping profile and kept-product status
    /// against `contract` (agreed again). The captured tolerance is restored.
    /// `Legacy` mode skips the identity checks with a printed warning.
    pub fn replay_checked(
        &mut self,
        n: &Array3<f64>,
        start: crate::snapshot::ReplayStart,
        contract: &crate::snapshot::ReplayContract,
    ) -> Result<(ScfStep, crate::snapshot::ReplayReport), ScfError> {
        use crate::snapshot::ReplayMode;
        if self.states.len() != 1 {
            return Err(ScfError::InvalidInput("replay needs one k point"));
        }
        let mut problems = Vec::new();
        if contract.mode == ReplayMode::Legacy {
            eprintln!(
                "snapshot replay {}: LEGACY manifest without operator fingerprint: operator, rows, potential and profile are NOT checked; not an exact replay",
                contract.label
            );
        } else {
            let current = self.operator_fingerprint(n);
            match &contract.operator {
                Some(captured) => problems.extend(captured.differences(&current)),
                None => problems.push("operator fingerprint: absent".to_string()),
            }
            let density = crate::snapshot::fnv(n.as_slice().unwrap_or(&[]));
            if contract.density_fnv != Some(density) {
                problems.push(format!(
                    "density bytes: manifest {} file {}",
                    contract
                        .density_fnv
                        .map_or("absent".into(), crate::snapshot::hex),
                    crate::snapshot::hex(density)
                ));
            }
        }
        // Agreement before any rank injects or solves.
        let mut any = f64::from(u8::from(!problems.is_empty()));
        if let Some(c) = crate::parallel::SharedComm::active() {
            use dft_mpi::Collective;
            let mut out = [0.0];
            c.tree_sum_f64(&[any], &mut out)
                .map_err(|e| ScfError::Numerical(format!("replay identity agreement: {e}")))?;
            any = out[0];
        }
        if any > 0.0 {
            return Err(ScfError::Numerical(format!(
                "snapshot {} ({} replay): operator identity mismatch: {}",
                contract.label,
                contract.mode.name(),
                if problems.is_empty() {
                    "on another rank".to_string()
                } else {
                    problems.join("; ")
                }
            )));
        }
        let bands = self.occupations.len();
        crate::snapshot::set_report(None);
        crate::snapshot::set_contract(Some(contract.clone()));
        crate::snapshot::inject(crate::profile::rank_key(), start);
        let placeholder =
            vec![ComplexMatrix::zeros(0, bands).map_err(dft_solvers::SolverError::from)?];
        let step = self.step_session(n, Some(&placeholder), &[None], contract.tolerance, false);
        // A solve that never reached the real Gamma entry leaves its contract.
        crate::snapshot::set_contract(None);
        let step = step?;
        let report = crate::snapshot::take_report().ok_or(ScfError::InvalidInput(
            "snapshot replay: the injected start did not reach the real Gamma solver",
        ))?;
        Ok((step, report))
    }

    /// The operator fingerprint of a map at input density `n` (ISSUE-18;
    /// computed only for snapshot capture and replay, never in a plain SCF).
    pub fn operator_fingerprint(&self, n: &Array3<f64>) -> crate::snapshot::OperatorFingerprint {
        use crate::snapshot::{fnv, hex};
        let words = |h: u64, v: &[u64]| {
            v.iter()
                .fold(h, |h, x| (h ^ x).wrapping_mul(0x0100_0000_01b3))
        };
        const SEED: u64 = 0xcbf2_9ce4_8422_2325;
        let mut fp = crate::snapshot::OperatorFingerprint::default();
        fp.push("version", crate::snapshot::FINGERPRINT_VERSION.to_string());
        fp.push("convention", crate::snapshot::CONVENTION);
        let cell = self.cell.lattice_columns_bohr();
        fp.push("cell", hex(fnv(cell.as_flattened())));
        let species: Vec<&str> = self.atoms.iter().map(|a| a.species.as_str()).collect();
        fp.push("species", species.join("+"));
        let positions: Vec<f64> = self
            .atoms
            .iter()
            .flat_map(|a| a.fractional_position)
            .collect();
        fp.push("positions", hex(fnv(&positions)));
        // Pseudopotential content as the operator sees it: local ionic field
        // (G space), core density, and per state the nonlocal projector
        // columns and couplings.
        let mut pseudo = SEED;
        for z in self.ionic_g.iter() {
            pseudo = words(pseudo, &[z.re.to_bits(), z.im.to_bits()]);
        }
        if let Some(core) = &self.core_density {
            pseudo = words(pseudo, &[core.len() as u64]);
            for x in core.iter() {
                pseudo = words(pseudo, &[x.to_bits()]);
            }
        }
        let mut basis = SEED;
        let mut basis_len = Vec::new();
        for state in &self.states {
            for atom in state.hamiltonian.nonlocal().atoms() {
                pseudo = words(pseudo, &[fnv(atom.coupling_hartree())]);
                for column in atom.columns() {
                    for z in column {
                        pseudo = words(pseudo, &[z.re.to_bits(), z.im.to_bits()]);
                    }
                }
            }
            let momenta = state.hamiltonian.basis_momenta();
            basis_len.push(momenta.len().to_string());
            for g in momenta {
                basis = words(basis, &g.miller.map(|m| m as i64 as u64));
            }
            basis = words(basis, &state.point.fractional.map(f64::to_bits));
        }
        fp.push("pseudo", hex(pseudo));
        fp.push(
            "xc",
            format!(
                "{}/{:?}/{:?}/{:?}",
                self.options.functional.name(),
                self.options.pw92_prefactor,
                self.options.local_mesh,
                self.options.radial_even_rule
            ),
        );
        fp.push("spin", "unpolarized");
        fp.push("occupations", hex(fnv(&self.occupations)));
        fp.push("electrons", hex(self.electrons.to_bits()));
        fp.push("cutoff", hex(self.options.cutoff_hartree.to_bits()));
        fp.push("basis", hex(basis));
        fp.push("basis_len", basis_len.join("+"));
        let d = self.grid.dims();
        fp.push("grid", format!("{}x{}x{}", d[0], d[1], d[2]));
        fp.push(
            "density_floor",
            self.density_floor
                .map_or("none".into(), |f| hex(f.to_bits())),
        );
        fp.push("symmetry_ops", self.force_model.symmetry.len().to_string());
        fp.push("density", hex(fnv(n.as_slice().unwrap_or(&[]))));
        fp
    }

    #[allow(clippy::too_many_arguments)]
    fn step_map(
        &mut self,
        n: &Array3<f64>,
        initial: Option<&[ComplexMatrix]>,
        initial_lazy: &[Option<crate::gamma_real::LazyOrbitals>],
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
        // The strict residual solve always materializes its result.
        let gathered_initial: Option<Vec<ComplexMatrix>>;
        let lazy_of = |i: usize| initial_lazy.get(i).and_then(Option::as_ref);
        let initial = match initial {
            Some(v) if v.iter().any(|m| m.rows() == 0) => {
                let mut full = Vec::with_capacity(v.len());
                for (i, m) in v.iter().enumerate() {
                    full.push(match lazy_of(i) {
                        Some(h) if m.rows() == 0 => h.materialize()?.clone(),
                        _ => m.clone(),
                    });
                }
                gathered_initial = Some(full);
                gathered_initial.as_deref()
            }
            other => other,
        };
        if n.shape() != self.grid.dims() {
            return Err(ScfError::InvalidInput("input density grid"));
        }
        // BREAKTHROUGH-0-INSTRUMENT: a new map for the phase profile; the
        // snapshot capture takes a copy of the input density (both default off).
        crate::profile::map_begin();
        {
            crate::snapshot::set_density(n);
        }
        // ISSUE-18: the operator fingerprint of the captured map (capture only).
        if crate::snapshot::dir().is_some() {
            crate::snapshot::set_operator(&self.operator_fingerprint(n));
        }
        let start = Instant::now();
        if let Some(symmetry) = &mut self.symmetry
            && !symmetry.check_invariant(n, &mut self.fft)?
        {
            return Err(ScfError::InvalidInput(
                "symmetry.density: input breaks the selected symmetry; use full mesh",
            ));
        }
        let (h, x, effective, dist_fields) = match self.reused_fields.take() {
            Some(r)
                if r.density.shape() == n.shape()
                    && r.density
                        .iter()
                        .zip(n.iter())
                        .all(|(a, b)| a.to_bits() == b.to_bits()) =>
            {
                (r.hartree, r.xc, r.effective, false)
            }
            _ => {
                self.local_fields_eligible(initial.is_some(), finite.is_some());
                self.map_fields(n)?
            }
        };
        {
            ordered_solves(
                &mut self.states,
                self.options.threads,
                &self.pool,
                |_, state| {
                    state.hamiltonian.set_local_potential(effective.clone())?;
                    Ok(())
                },
            )?;
            if crate::m2::active().is_some() {
                use crate::grid_dist as g;
                g::owned_add(g::O_FULL_COPIES, self.states.len() as u64);
                g::owned_add(
                    g::O_FULL_COPY_VALUES,
                    (self.states.len() * effective.len()) as u64,
                );
            }
        }
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
            } else {
                self.atomic_guesses()?
            }
        } else {
            None
        };
        let initial = initial.or(atomic_initial.as_deref());
        let mut timings = PhaseTimings {
            potential: start.elapsed(),
            ..Default::default()
        };
        crate::profile::map_phase(
            crate::profile::Phase::Potential,
            timings.potential.as_secs_f64(),
        );
        crate::profile::solve_begin();
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
        // Frozen empty-band policy: in SCF steps whose
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
                        lazy_of(i),
                        state.hamiltonian.kinetic_diagonal(),
                        self.occupations.len(),
                        tolerance,
                        policy,
                        &self.occupations,
                        self.electrons,
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
                        // M2-OWNED: the complex solver applies the full local
                        // potential; gather the owned one first (the same
                        // gather as `ensure_complex_potential`, replicated
                        // fallback decision).

                        // A lazy warm start of a solve that left the real path
                        // (replicated options) is gathered first.
                        let gathered;
                        let warm = match (initial.map(|v| &v[i]), lazy_of(i)) {
                            (Some(m), Some(h)) if m.rows() == 0 => {
                                gathered = h.materialize()?.clone();
                                Some(&gathered)
                            }
                            (w, _) => w,
                        };
                        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                            if self.options.locking == dft_solvers::locking::Locking::Auto {
                                solver.solve_locked(
                                    &state.hamiltonian,
                                    warm.or(guess.as_ref()),
                                    Some(&preconditioner),
                                    policy,
                                )
                            } else {
                                self.options
                                    .eigensolver
                                    .solve_with_band_tolerances_and_certification(
                                        &solver,
                                        &state.hamiltonian,
                                        warm.or(guess.as_ref()),
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
                        lazy: None,
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
                            None,
                            state.hamiltonian.kinetic_diagonal(),
                            bands,
                            tolerance,
                            policy,
                            &self.occupations,
                            self.electrons,
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
                            if self.options.locking == dft_solvers::locking::Locking::Auto {
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
                            lazy: None,
                        })
                    },
                )?
            };
        timings.diagonalization = start.elapsed();

        // Snapshot file writing (inside the solve window) is not solve time.
        let io = crate::snapshot::take_io_seconds();
        if io > 0.0 {
            timings.diagonalization = timings
                .diagonalization
                .saturating_sub(std::time::Duration::from_secs_f64(io));
        }
        crate::profile::map_phase(
            crate::profile::Phase::Solve,
            timings.diagonalization.as_secs_f64(),
        );
        crate::profile::solve_end(
            k_points.iter().map(|k| k.solution.iterations).sum(),
            k_points
                .iter()
                .map(|k| k.solution.operator_applications)
                .sum(),
        );
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
        let mut k_points = k_points;
        if !m2 && k_points.iter().any(|k| k.lazy.is_some()) {
            // MPI-LAZY-ORBITALS: the replicated band density needs the block.
            materialize_kpoints(&mut k_points)?;
        }
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
                    let weights = finite.map_or(occupations.as_slice(), |w| w.block(0, i));
                    let real = state.gamma_real.gamma_real_used;
                    state.band_density.compute_with(
                        &orbitals,
                        weights,
                        &state.grid,
                        volume,
                        real,
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
        crate::profile::map_phase(
            crate::profile::Phase::Density,
            timings.density.as_secs_f64(),
        );
        let start = Instant::now();
        // M2-GRID: the audit terms from the distributed output density (the
        // plan holds V_H and V_xc of this map's input on its planes).
        let dist_terms = match (
            dist_fields,
            crate::m2::active(),
            self.grid_dist.as_mut().filter(|p| p.has_fields()),
            n.as_slice(),
            output.as_slice(),
            self.ionic_real.as_slice(),
        ) {
            (true, Some(comm), Some(plan), Some(input), Some(out), Some(ionic)) => {
                use crate::grid_dist::Field;
                let in_field = Field::Full(input);
                let out_field = Field::Full(out);
                Some(plan.audit(&comm, in_field, out_field, ionic)?)
            }
            _ => None,
        };

        let m2_xc = if dist_terms.is_some() {
            None
        } else {
            crate::m2::active().and_then(|comm| {
                crate::m2::xc(
                    &comm,
                    &self.options,
                    &output,
                    self.core_density.as_ref(),
                    self.density_floor,
                    volume,
                )
            })
        };
        let output_xc = if let Some(x) = m2_xc {
            Some(x?)
        } else if dist_terms.is_some() {
            None
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
        let shared_audit = if let Some(terms) = dist_terms {
            energy::SharedAudit::from_distributed(terms.values(), self.ionic_g[[0, 0, 0]].re)
        } else {
            energy::prepare_shared(
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
            )?
        };

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
                let audit = energy::evaluate_shared_with_kinetic(
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
                    None,
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
        crate::profile::map_phase(crate::profile::Phase::Energy, timings.energy.as_secs_f64());
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
    /// RUN-CONFIG: the run configuration of this problem.
    pub fn run_config(&self) -> dft_numerics::switches::RunConfig {
        self.run
    }
    /// RUN-CONFIG: run this problem under `run` (the explicit per-call way;
    /// construction took the calling thread's configuration). The session
    /// switches and the construction-time choices that depend on run switches
    /// (`BAANDI_STOP_METRIC`) are resolved again from it; sessions enter it as
    /// their thread scope while they initialize and step, so code that still
    /// reads the legacy accessors sees it too. Construction-only switches
    /// (`BAANDI_EWALD_CACHE`, the Ewald and ionic schemes) keep the values the
    /// construction used: to construct under `run`, enter `run.enter()` first.
    pub fn set_run_config(&mut self, run: dft_numerics::switches::RunConfig) {
        self.run = run;
        self.sw = session::SessionSwitches::resolve(run);
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
    let intermediate = tolerance > 0.1 * options.eigen_tolerance * (1.0 + 1e-9);
    BandTolerances {
        strict_bands: (occupied_count + usize::from(!intermediate))
            .min(bands)
            .max(1),
        buffer_tolerance: if intermediate {
            tolerance.max(1e-3)
        } else {
            tolerance.max(1e-5)
        },
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
use crate::pool_fixture;

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
                        ..ScfOptions::strict_reference()
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
        include!("../../tests/common/persistent_real_state_unit.rs");
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
                ..ScfOptions::strict_reference()
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
            .certificate_shift_estimate(&step.density, &step.k_points)
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
        let run = |threads: usize| {
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
                    ..ScfOptions::strict_reference()
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
                .certificate_shift_estimate(&step.density, &step.k_points)
                .unwrap()
                .unwrap();
            assert!(problem.reused_fields.is_some());
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
        let reference = run(1);
        assert_eq!(reference, run(4));
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
                    ..ScfOptions::strict_reference()
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
                    ..ScfOptions::strict_reference()
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
        let o = ScfOptions::strict_reference();
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
            pulay: None,
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

/// Diagnostic fixed-potential maps, opt-in only.
#[path = "fixed_map.rs"]
pub mod fixed_map;
