//! Track A, stages 0 and 1 (speed-restructure proposal, section 3): all-band
//! direct minimization of the Kohn-Sham energy over the occupied subspace, used
//! only as a PRE-CONVERGER for the existing strict SCF session.
//!
//! Scope: one exact Gamma point on the real (ORC) path, fixed integer
//! occupations, norm-conserving, the engine's own functional (LDA first). Every
//! other problem falls back to the ordinary SCF path with a recorded reason.
//!
//! Placement: a child module of `scf` (not `dft-solvers`), because the energy
//! needs the problem's private ionic potential, Ewald term, Hartree/XC setup and
//! Hamiltonian state; `dft-solvers` sits below `dft-physics`/`dft-engine` in the
//! dependency graph and cannot see them.
//!
//! Variant A0: energy and gradient from the orbitals (density from the real pair
//! path, potential from the existing Hartree/XC code, G = H C from
//! `RealHamiltonian::apply_block`), projected gradient R = G - C (C^T G),
//! kinetic (TPA) preconditioning, Polak-Ribiere conjugate gradient with reset,
//! Loewdin retraction C(t) = orth(C + t D) and a parabolic line search that
//! accepts only a strictly lower evaluated energy. The minimizer never declares
//! convergence of the calculation: its output is handed to `ScfSession`, whose
//! unchanged strict steps and fixed-point certificate produce the result.
//!
//! Known caveat (docs/validation/scf-fixed-point-multiplicity.md): a minimizer
//! descends to a local minimum of the energy, while a mixed SCF iteration can
//! stop at another self-consistent solution. On systems with several fixed
//! points the two paths can certify different states; `StartPath` records which
//! path produced the start state.
//!
//! Hot-loop rules: all orbital-sized buffers are allocated before the loop. The
//! Hartree/XC evaluation and `Hamiltonian::set_local_potential` are the
//! existing, allocating engine APIs (the same calls the SCF step makes); the
//! dense symmetric EVD may allocate faer bookkeeping (see `RealEighWorkspace`).
//! No logging or I/O. Every reduction has a fixed order, so results do not
//! depend on the thread count.
use super::*;
use crate::checkpoint::{ElectronicState, Orbitals};
use dft_basis::HalfSphere;
use dft_numerics::{
    PairDensityWorkspace, PairMap, RealGemmWorkspace, RealMatrix, RealOp,
    fused_pair_plan::RealPairPlan, real::RealEighWorkspace,
};
use dft_physics::real_gamma::RealHamiltonian;

fn num(e: impl fmt::Display) -> ScfError {
    ScfError::Numerical(format!("direct_min: {e}"))
}

/// Engine-level method switch. The default keeps the ordinary SCF path; there is
/// deliberately no CLI key yet (`[scf] method` is a later, separate decision).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ScfMethod {
    #[default]
    Scf,
    /// Direct minimization pre-converger, then the strict SCF session.
    DirectThenScf,
}

/// Minimizer controls. Defaults are stage-0 choices, not tuned on any system.
#[derive(Clone, Debug)]
pub struct DirectMinOptions {
    /// Stop when the largest column norm of the projected gradient (Ha) is below this.
    pub gradient_tolerance: f64,
    pub max_iterations: usize,
    /// Steepest-descent restart period of the conjugate-gradient recursion.
    pub cg_reset_period: usize,
    /// Stop as stagnated if the gradient norm did not fall below 0.9 times its
    /// value this many iterations earlier.
    pub stagnation_window: usize,
    /// Energy evaluations allowed in one line search before it fails.
    pub max_line_search_evaluations: usize,
    /// First trial step length (dimensionless, along the preconditioned direction).
    pub initial_step: f64,
}
impl Default for DirectMinOptions {
    fn default() -> Self {
        Self {
            // Stage-0 sweep 1e-5/1e-6/1e-7 on six fixtures (Si and synthetic,
            // docs/performance/track-a-stage0.md): 1e-5 leaves a density error
            // that makes the strict certificate fail repeatedly; 1e-7 reaches the
            // energy rounding floor (failed line searches). Not tuned on Si alone.
            gradient_tolerance: 1e-6,
            max_iterations: 300,
            cg_reset_period: 50,
            stagnation_window: 30,
            max_line_search_evaluations: 8,
            initial_step: 1.0,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum StopReason {
    /// Gradient tolerance reached (the minimizer's own stop; never a certificate).
    GradientTolerance,
    MaxIterations,
    Stagnated,
    /// No evaluated point in the line search had a lower energy. The state is the
    /// last accepted one; the strict session starts from it.
    LineSearchFailed,
    /// The orthonormality defect of an accepted C exceeded 1e-10.
    OrthonormalityDefect,
}

/// Work and outcome counters of one minimization. H columns are real ORC
/// columns applied by `RealHamiltonian::apply_block` (plus the complex atomic
/// start columns, counted one each).
#[derive(Clone, Debug)]
pub struct DirectMinReport {
    pub occupied_bands: usize,
    pub start: &'static str,
    pub start_h_columns: usize,
    pub iterations: usize,
    /// Energy evaluations (each one apply of all occupied columns), start included.
    pub energy_evaluations: usize,
    pub h_columns: usize,
    pub cg_resets: usize,
    pub energy_hartree: f64,
    pub projected_gradient: f64,
    pub max_orthonormality_defect: f64,
    pub stop: StopReason,
    /// Accepted energies, start first (Ha); non-increasing by construction.
    pub energies: Vec<f64>,
    /// Max column norm of the projected gradient per iteration (Ha).
    pub gradients: Vec<f64>,
}
impl DirectMinReport {
    /// H applications per occupied band (H-equivalents per band).
    pub fn h_equivalents_per_band(&self) -> f64 {
        self.h_columns as f64 / self.occupied_bands as f64
    }
}

/// Which path produced the start state of the strict session.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum StartPath {
    /// Ordinary SCF start (default method).
    Scf,
    /// Density and Ritz-rotated occupied orbitals from the direct minimizer.
    DirectMinimization,
    /// DirectThenScf was requested but the problem is ineligible; ordinary SCF ran.
    FallbackScf(String),
}

pub struct MethodOutcome {
    pub outcome: ScfOutcome,
    pub start_path: StartPath,
    pub minimizer: Option<DirectMinReport>,
}

/// Minimizer output: occupied orbitals (real ORC columns, Ritz rotated,
/// ascending Ritz values), their density and the report.
pub struct DirectMinState {
    pub orbitals: RealMatrix,
    pub ritz_values: Vec<f64>,
    pub density: Array3<f64>,
    /// Buffer-band guesses for the hand-off (real columns, may be empty): the
    /// unoccupied atomic Ritz vectors of the start, when the start was atomic.
    pub buffer_guess: Vec<f64>,
    pub report: DirectMinReport,
}

/// Why a problem cannot use the minimizer (None means eligible).
pub fn ineligibility(problem: &ScfProblem) -> Option<String> {
    if problem.states.len() != 1 {
        return Some("multiple k points".into());
    }
    if let Some(reason) = problem.gamma_real_exclusion {
        return Some(reason.into());
    }
    if problem.states[0].basis.k_fractional() != [0.; 3] {
        return Some("k is not exactly Gamma".into());
    }
    if problem.options.gamma_real == crate::GammaReal::Off {
        return Some("gamma_real off".into());
    }
    if problem.symmetry.is_some() {
        return Some("density symmetrization enabled".into());
    }
    if problem.density_floor.is_some() {
        return Some("XC density floor enabled".into());
    }
    let occupied = problem.occupations.iter().filter(|&&f| f > 0.).count();
    if occupied == 0
        || problem
            .occupations
            .iter()
            .any(|&f| f != 0. && f != problem.occupations[0])
        || problem.occupations[..occupied].iter().any(|&f| f <= 0.)
    {
        return Some("occupations are not fixed and uniform".into());
    }
    None
}

/// Run the selected method. `ScfMethod::Scf` is exactly `run_outcome`/`run`.
pub fn run_with_method(
    problem: ScfProblem,
    method: ScfMethod,
    options: &DirectMinOptions,
    retain_last: bool,
) -> Result<MethodOutcome, ScfError> {
    let fallback = match method {
        ScfMethod::Scf => None,
        ScfMethod::DirectThenScf => Some(ineligibility(&problem)),
    };
    let reason = match fallback {
        None => {
            let outcome = finish(problem, None, retain_last)?;
            return Ok(MethodOutcome {
                outcome,
                start_path: StartPath::Scf,
                minimizer: None,
            });
        }
        Some(reason) => reason,
    };
    if let Some(reason) = reason {
        let outcome = finish(problem, None, retain_last)?;
        return Ok(MethodOutcome {
            outcome,
            start_path: StartPath::FallbackScf(reason),
            minimizer: None,
        });
    }
    let mut problem = problem;
    let state = minimize(&mut problem, options)?;
    let electronic = hand_off_state(&problem, &state)?;
    let report = state.report;
    let outcome = finish(problem, Some(&electronic), retain_last)?;
    Ok(MethodOutcome {
        outcome,
        start_path: StartPath::DirectMinimization,
        minimizer: Some(report),
    })
}

/// Drive the unchanged strict session, optionally from an imported state
/// (`hand_off_state`). Only the session's fixed-point certificate converges.
pub fn finish(
    problem: ScfProblem,
    start: Option<&ElectronicState>,
    retain_last: bool,
) -> Result<ScfOutcome, ScfError> {
    let completion = if retain_last {
        CompletionMode::RetainLast
    } else {
        CompletionMode::Strict
    };
    let mut session = ScfSession::new(
        problem,
        completion,
        crate::events::ExecutionConfig::default(),
    )?;
    match start {
        Some(state) => session.initialize_electronic(state)?,
        None => session.initialize()?,
    }
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

/// Build the hand-off state: minimizer density, Ritz-rotated occupied orbitals
/// and buffer columns (the start's unoccupied atomic Ritz vectors, then
/// lowest-kinetic ORC unit vectors, each orthogonalized twice against all
/// previous columns), Hermitian-completed into the complex basis.
pub fn hand_off_state(
    problem: &ScfProblem,
    state: &DirectMinState,
) -> Result<ElectronicState, ScfError> {
    let k = &problem.states[0];
    let half = HalfSphere::new(&k.basis, k.hamiltonian.fft_grid()).map_err(num)?;
    let n = half.full_len();
    let m = state.orbitals.cols();
    let bands = problem.occupations.len();
    let mut real = vec![0.; n * bands];
    real[..n * m].copy_from_slice(state.orbitals.as_slice());
    let kinetic = real_kinetic(&half, k.hamiltonian.kinetic_diagonal());
    let mut order: Vec<_> = (0..n).collect();
    order.sort_by(|&a, &b| kinetic[a].total_cmp(&kinetic[b]).then(a.cmp(&b)));
    let mut filled = m;
    let mut v = vec![0.; n];
    let guesses = state.buffer_guess.len() / n;
    for candidate in 0..guesses + n {
        if filled == bands {
            break;
        }
        if candidate < guesses {
            v.copy_from_slice(&state.buffer_guess[candidate * n..(candidate + 1) * n]);
        } else {
            v.fill(0.);
            v[order[candidate - guesses]] = 1.;
        }
        for _ in 0..2 {
            for j in 0..filled {
                let col = &real[j * n..(j + 1) * n];
                let dot: f64 = col.iter().zip(&v).map(|(a, b)| a * b).sum();
                for (x, a) in v.iter_mut().zip(col) {
                    *x -= dot * a;
                }
            }
        }
        let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm > 0.1 {
            for (dst, x) in real[filled * n..(filled + 1) * n].iter_mut().zip(&v) {
                *dst = x / norm;
            }
            filled += 1;
        }
    }
    if filled != bands {
        return Err(num("cannot complete buffer bands"));
    }
    let mut full = vec![Complex64::default(); n];
    let mut coefficients = Vec::with_capacity(2 * n * bands);
    for j in 0..bands {
        half.unpack(&real[j * n..(j + 1) * n], &mut full)
            .map_err(num)?;
        coefficients.extend(full.iter().flat_map(|z| [z.re, z.im]));
    }
    Ok(ElectronicState {
        grid: problem.grid.dims(),
        density: state.density.iter().copied().collect(),
        orbitals: vec![Orbitals {
            rows: n,
            cols: bands,
            coefficients,
        }],
        occupations: problem.occupations.clone(),
    })
}

fn real_kinetic(half: &HalfSphere, complex: &[f64]) -> Vec<f64> {
    let mut kinetic = vec![0.; half.full_len()];
    kinetic[0] = complex[0];
    for (i, &p) in half.plus().iter().enumerate().skip(1) {
        kinetic[2 * i - 1] = complex[p as usize];
        kinetic[2 * i] = complex[p as usize];
    }
    kinetic
}

/// Energy evaluation resources (borrowed disjointly from the workspace).
struct Evaluator {
    h: RealHamiltonian,
    density: PairDensityWorkspace<RealPairPlan>,
    weights: Vec<f64>,
    volume: f64,
    applied_columns: usize,
    evaluations: usize,
}
impl Evaluator {
    /// E[C] = sum_j f_j c_j.(H c_j) - E_H[n] + E_xc[n] - int n v_xc + E_ewald,
    /// with H built from n = n(C); writes G = H C and n. Equals the audit's
    /// direct functional (kinetic+local+nonlocal+Hartree+XC+Ewald) for this C.
    fn energy(
        &mut self,
        problem: &mut ScfProblem,
        c: &RealMatrix,
        g: &mut RealMatrix,
        n: &mut Array3<f64>,
    ) -> Result<f64, ScfError> {
        self.density
            .density(
                c.as_slice(),
                &self.weights,
                self.volume,
                n.as_slice_mut().expect("standard layout"),
            )
            .map_err(num)?;
        if n.iter().any(|v| !v.is_finite()) {
            return Err(ScfError::NonFinite);
        }
        let (hartree, xc, effective) = problem.input_fields(n)?;
        let state = &mut problem.states[0];
        state.hamiltonian.set_local_potential(effective)?;
        if !self.h.update_from_complex(&state.hamiltonian, &state.basis) {
            return Err(num("real Hamiltonian lost its source identity"));
        }
        self.h.apply_block(c, g).map_err(num)?;
        self.applied_columns += c.cols();
        self.evaluations += 1;
        let rows = c.rows();
        let mut band = 0.;
        for (j, f) in self.weights.iter().enumerate() {
            let cj = &c.as_slice()[j * rows..(j + 1) * rows];
            let gj = &g.as_slice()[j * rows..(j + 1) * rows];
            band += f * cj.iter().zip(gj).map(|(a, b)| a * b).sum::<f64>();
        }
        let weight = self.volume / n.len() as f64;
        let xc_double = weight * n.iter().zip(&xc.potential).map(|(n, v)| n * v).sum::<f64>();
        let energy = band - hartree.energy.0 + xc.energy.0 - xc_double + problem.ewald.energy.0;
        if !energy.is_finite() {
            return Err(ScfError::NonFinite);
        }
        Ok(energy)
    }
}

/// Small dense scratch for Gram/rotation work (m x m), reused every iteration.
struct Dense {
    gemm: RealGemmWorkspace,
    eigh: RealEighWorkspace,
    a: RealMatrix,
    b: RealMatrix,
    v: RealMatrix,
    values: Vec<f64>,
}
impl Dense {
    /// out = y (y^T y)^{-1/2} (Loewdin). Returns the smallest Gram eigenvalue.
    fn lowdin(&mut self, y: &RealMatrix, out: &mut RealMatrix) -> Result<f64, ScfError> {
        let m = y.cols();
        self.gemm
            .multiply_into(y.as_ref(), RealOp::Transpose, y.as_ref(), self.a.as_mut())
            .map_err(num)?;
        symmetrize(&mut self.a);
        self.eigh
            .solve_into(&self.a, &mut self.values, &mut self.v)
            .map_err(num)?;
        let smallest = self.values[0];
        if smallest.is_nan() || smallest <= 1e-12 {
            return Err(num("rank-deficient retraction"));
        }
        for j in 0..m {
            for i in 0..m {
                let mut s = 0.;
                for k in 0..m {
                    s += self.v[(i, k)] * self.v[(j, k)] / self.values[k].sqrt();
                }
                self.b[(i, j)] = s;
            }
        }
        self.gemm
            .multiply_into(y.as_ref(), RealOp::Identity, self.b.as_ref(), out.as_mut())
            .map_err(num)?;
        Ok(smallest)
    }
    /// a = x^T y (symmetrized when `symmetric`).
    fn gram(&mut self, x: &RealMatrix, y: &RealMatrix, symmetric: bool) -> Result<(), ScfError> {
        self.gemm
            .multiply_into(x.as_ref(), RealOp::Transpose, y.as_ref(), self.a.as_mut())
            .map_err(num)?;
        if symmetric {
            symmetrize(&mut self.a);
        }
        Ok(())
    }
    /// out = x a.
    fn rotate(&mut self, x: &RealMatrix, out: &mut RealMatrix) -> Result<(), ScfError> {
        self.gemm
            .multiply_into(x.as_ref(), RealOp::Identity, self.a.as_ref(), out.as_mut())
            .map_err(num)
    }
}
fn symmetrize(a: &mut RealMatrix) {
    let m = a.rows();
    for j in 0..m {
        for i in 0..j {
            let v = 0.5 * (a[(i, j)] + a[(j, i)]);
            a[(i, j)] = v;
            a[(j, i)] = v;
        }
    }
}
fn dot(a: &RealMatrix, b: &RealMatrix) -> f64 {
    a.as_slice()
        .iter()
        .zip(b.as_slice())
        .map(|(x, y)| x * y)
        .sum()
}

struct Trial {
    c: RealMatrix,
    g: RealMatrix,
    n: Array3<f64>,
    energy: f64,
    step: f64,
}

/// Initial occupied block: atomic Rayleigh-Ritz at the starting density when the
/// problem's start is atomic and available, else lowest-kinetic ORC unit vectors.
fn initial_block(
    problem: &mut ScfProblem,
    half: &HalfSphere,
    kinetic: &[f64],
    m: usize,
    c: &mut RealMatrix,
    buffer: &mut Vec<f64>,
) -> Result<(&'static str, usize), ScfError> {
    let n = half.full_len();
    if problem.options.start_wavefunctions == StartWavefunctions::Atomic
        && problem.atomic_species.is_some()
    {
        let start = problem.starting_density();
        let (_, _, effective) = problem.input_fields(&start)?;
        let species = problem.atomic_species.as_ref().expect("checked above");
        let state = &mut problem.states[0];
        state.hamiltonian.set_local_potential(effective)?;
        let span = atomic_start::span(
            &problem.cell,
            &state.basis,
            &state.grid,
            species,
            &problem.options,
            problem.occupations.len(),
            None,
        )?;
        let (guess, report) = atomic_start::ritz(
            &span,
            &state.hamiltonian,
            problem.occupations.len(),
            true,
            None,
        )?;
        let mut ok = true;
        for j in 0..m {
            if half
                .pack(guess.column(j), &mut c.as_mut_slice()[j * n..(j + 1) * n])
                .is_err()
            {
                ok = false;
                break;
            }
        }
        if ok {
            let mut column = vec![0.; n];
            for j in m..guess.cols() {
                if half.pack(guess.column(j), &mut column).is_ok() {
                    buffer.extend_from_slice(&column);
                }
            }
            return Ok(("atomic", report.initial_h_applications));
        }
    }
    let mut order: Vec<_> = (0..n).collect();
    order.sort_by(|&a, &b| kinetic[a].total_cmp(&kinetic[b]).then(a.cmp(&b)));
    c.as_mut_slice().fill(0.);
    for (j, &i) in order.iter().take(m).enumerate() {
        c[(i, j)] = 1.;
    }
    // A pure unit-vector block can lie in one symmetry sector; the gradient flow
    // then never leaves it (observed: syn-fcc reached a state 0.049 Ha above the
    // certified SCF one). Add a seeded, kinetically damped component on every
    // row so that no sector is empty. Deterministic (xorshift on options.seed).
    let mut x = problem.options.seed ^ 0x9e37_79b9_7f4a_7c15;
    for j in 0..m {
        for (i, t) in kinetic.iter().enumerate() {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            let u = (x >> 11) as f64 / (1u64 << 53) as f64 - 0.5;
            c[(i, j)] += 0.1 * u / (1. + t * t);
        }
    }
    Ok(("free_electron_perturbed", 0))
}

/// Run variant A0 on an eligible problem. The problem's Hamiltonian potential is
/// left at the last evaluated density; the strict session rebuilds it.
pub fn minimize(
    problem: &mut ScfProblem,
    options: &DirectMinOptions,
) -> Result<DirectMinState, ScfError> {
    if let Some(reason) = ineligibility(problem) {
        return Err(ScfError::Numerical(format!(
            "direct_min: ineligible problem: {reason}"
        )));
    }
    if options.gradient_tolerance.is_nan()
        || options.gradient_tolerance <= 0.
        || options.max_iterations == 0
        || options.cg_reset_period == 0
        || options.stagnation_window == 0
        || options.max_line_search_evaluations < 2
        || options.initial_step.is_nan()
        || options.initial_step <= 0.
    {
        return Err(ScfError::InvalidInput("direct_min options"));
    }
    let m = problem.occupations.iter().filter(|&&f| f > 0.).count();
    let f = problem.occupations[0];
    let pool = problem.pool.clone();
    let volume = problem.cell.volume_bohr3();
    let dims = problem.grid.dims();
    let (half, kinetic, h, density_ws) = {
        let k = &problem.states[0];
        let grid = k.hamiltonian.fft_grid();
        let half = HalfSphere::new(&k.basis, grid).map_err(num)?;
        let kinetic = real_kinetic(&half, k.hamiltonian.kinetic_diagonal());
        let h = RealHamiltonian::from_complex(&k.hamiltonian, &k.basis, m, pool.clone())
            .map_err(num)?;
        let map =
            PairMap::new(grid.dims(), grid.g_to_fft(), half.plus(), half.minus()).map_err(num)?;
        let density_ws = PairDensityWorkspace::new(map, m, pool.clone(), || {
            RealPairPlan::new(grid.dims(), grid.g_to_fft())
        })
        .map_err(num)?;
        (half, kinetic, h, density_ws)
    };
    let n = half.full_len();
    if m > n {
        return Err(num("more occupied bands than real coordinates"));
    }
    let zeros = || RealMatrix::zeros(n, m).map_err(num);
    let mut eval = Evaluator {
        h,
        density: density_ws,
        weights: vec![f; m],
        volume,
        applied_columns: 0,
        evaluations: 0,
    };
    let mut dense = Dense {
        gemm: RealGemmWorkspace::with_pool(n, m, pool).map_err(num)?,
        eigh: RealEighWorkspace::new(m).map_err(num)?,
        a: RealMatrix::zeros(m, m).map_err(num)?,
        b: RealMatrix::zeros(m, m).map_err(num)?,
        v: RealMatrix::zeros(m, m).map_err(num)?,
        values: vec![0.; m],
    };
    let (mut c, mut g, mut r, mut z, mut d, mut r_old, mut y) = (
        zeros()?,
        zeros()?,
        zeros()?,
        zeros()?,
        zeros()?,
        zeros()?,
        zeros()?,
    );
    let trial = || -> Result<Trial, ScfError> {
        Ok(Trial {
            c: zeros()?,
            g: zeros()?,
            n: Array3::zeros(dims),
            energy: f64::INFINITY,
            step: 0.,
        })
    };
    let mut trials = [trial()?, trial()?];
    let mut density = Array3::zeros(dims);
    let mut buffer_guess = Vec::new();
    let mut report = DirectMinReport {
        occupied_bands: m,
        start: "",
        start_h_columns: 0,
        iterations: 0,
        energy_evaluations: 0,
        h_columns: 0,
        cg_resets: 0,
        energy_hartree: f64::NAN,
        projected_gradient: f64::NAN,
        max_orthonormality_defect: 0.,
        stop: StopReason::MaxIterations,
        energies: Vec::with_capacity(options.max_iterations + 1),
        gradients: Vec::with_capacity(options.max_iterations + 1),
    };
    // Start: initial block, Loewdin orthonormalization, first energy.
    let (start, start_columns) =
        initial_block(problem, &half, &kinetic, m, &mut y, &mut buffer_guess)?;
    report.start = start;
    report.start_h_columns = start_columns;
    dense.lowdin(&y, &mut c)?;
    let mut energy = eval.energy(problem, &c, &mut g, &mut density)?;
    report.energies.push(energy);
    let mut step = options.initial_step;
    let mut zr_old = 0.;
    let mut have_direction = false;
    let mut since_reset = 0usize;
    let mut scale = vec![0.; m];
    for iteration in 0..options.max_iterations {
        // Orthonormality defect of the accepted block.
        dense.gram(&c, &c, false)?;
        let mut defect: f64 = 0.;
        for j in 0..m {
            for i in 0..m {
                let target = if i == j { 1. } else { 0. };
                defect = defect.max((dense.a[(i, j)] - target).abs());
            }
        }
        report.max_orthonormality_defect = report.max_orthonormality_defect.max(defect);
        if defect > 1e-10 {
            report.stop = StopReason::OrthonormalityDefect;
            break;
        }
        // Projected gradient R = G - C (C^T G).
        dense.gram(&c, &g, true)?;
        dense.rotate(&c, &mut r)?;
        for (r, g) in r.as_mut_slice().iter_mut().zip(g.as_slice()) {
            *r = g - *r;
        }
        let gnorm = (0..m)
            .map(|j| {
                r.as_slice()[j * n..(j + 1) * n]
                    .iter()
                    .map(|x| x * x)
                    .sum::<f64>()
                    .sqrt()
            })
            .fold(0., f64::max);
        report.gradients.push(gnorm);
        report.projected_gradient = gnorm;
        if gnorm <= options.gradient_tolerance {
            report.stop = StopReason::GradientTolerance;
            break;
        }
        let w = options.stagnation_window;
        if iteration >= w && gnorm > 0.9 * report.gradients[iteration - w] {
            report.stop = StopReason::Stagnated;
            break;
        }
        // TPA preconditioning per band, then projection out of span(C).
        for (j, s) in scale.iter_mut().enumerate() {
            let cj = &c.as_slice()[j * n..(j + 1) * n];
            *s = dft_solvers::kinetic_scale(cj.iter().zip(&kinetic).map(|(x, t)| (x * x, *t)))
                .map_err(num)?;
        }
        for j in 0..m {
            for i in 0..n {
                z[(i, j)] = dft_solvers::tpa_weight(kinetic[i], scale[j]) * r[(i, j)];
            }
        }
        dense.gram(&c, &z, false)?;
        dense.rotate(&c, &mut y)?;
        for (z, y) in z.as_mut_slice().iter_mut().zip(y.as_slice()) {
            *z -= y;
        }
        let zr = dot(&z, &r);
        // Polak-Ribiere (non-negative) with periodic reset; old direction is
        // transported by projection onto the new tangent space.
        let mut beta = 0.;
        if have_direction && since_reset + 1 < options.cg_reset_period && zr_old > 0. {
            beta = ((zr - dot(&z, &r_old)) / zr_old).max(0.);
        }
        if beta == 0. {
            report.cg_resets += 1;
            since_reset = 0;
        } else {
            since_reset += 1;
        }
        if beta > 0. {
            dense.gram(&c, &d, false)?;
            dense.rotate(&c, &mut y)?;
            for ((d, y), z) in d
                .as_mut_slice()
                .iter_mut()
                .zip(y.as_slice())
                .zip(z.as_slice())
            {
                *d = -z + beta * (*d - y);
            }
        } else {
            for (d, z) in d.as_mut_slice().iter_mut().zip(z.as_slice()) {
                *d = -z;
            }
        }
        let mut slope = 2. * f * dot(&d, &g);
        if slope.is_nan() || slope >= 0. {
            for (d, z) in d.as_mut_slice().iter_mut().zip(z.as_slice()) {
                *d = -z;
            }
            report.cg_resets += 1;
            since_reset = 0;
            slope = 2. * f * dot(&d, &g);
        }
        std::mem::swap(&mut r, &mut r_old);
        zr_old = zr;
        have_direction = true;
        // Line search: parabola through E(0), E'(0) and E(t); accept only a
        // strictly lower evaluated energy, otherwise shrink.
        let mut accepted: Option<usize> = None;
        let mut evaluations = 0;
        let mut t = step;
        while evaluations + 1 < options.max_line_search_evaluations && accepted.is_none() {
            retract(&mut dense, &c, &d, t, &mut y, &mut trials[0])?;
            trials[0].energy =
                eval.energy(problem, &trials[0].c, &mut trials[0].g, &mut trials[0].n)?;
            evaluations += 1;
            let e1 = trials[0].energy;
            let curvature = (e1 - energy - slope * t) / (t * t);
            if curvature > 0. {
                let fitted = (-slope / (2. * curvature)).clamp(0.02 * t, 20. * t);
                if (fitted - t).abs() <= 0.1 * t && e1 < energy {
                    accepted = Some(0);
                    break;
                }
                retract(&mut dense, &c, &d, fitted, &mut y, &mut trials[1])?;
                trials[1].energy =
                    eval.energy(problem, &trials[1].c, &mut trials[1].g, &mut trials[1].n)?;
                evaluations += 1;
                let best = if trials[1].energy <= trials[0].energy {
                    1
                } else {
                    0
                };
                if trials[best].energy < energy {
                    accepted = Some(best);
                } else {
                    t = 0.25 * fitted.min(t);
                }
            } else if e1 < energy {
                accepted = Some(0);
            } else {
                t *= 0.25;
            }
        }
        report.iterations = iteration + 1;
        let Some(best) = accepted else {
            report.stop = StopReason::LineSearchFailed;
            break;
        };
        let trial = &mut trials[best];
        std::mem::swap(&mut c, &mut trial.c);
        std::mem::swap(&mut g, &mut trial.g);
        std::mem::swap(&mut density, &mut trial.n);
        energy = trial.energy;
        step = trial.step;
        report.energies.push(energy);
        if iteration + 1 == options.max_iterations {
            report.stop = StopReason::MaxIterations;
        }
    }
    // Final Ritz rotation: C <- C V, G <- G V with V the eigenvectors of C^T G.
    dense.gram(&c, &g, true)?;
    let mut ritz = vec![0.; m];
    dense
        .eigh
        .solve_into(&dense.a, &mut ritz, &mut dense.v)
        .map_err(num)?;
    std::mem::swap(&mut dense.a, &mut dense.v);
    dense.rotate(&c, &mut y)?;
    std::mem::swap(&mut c, &mut y);
    report.energy_hartree = energy;
    report.energy_evaluations = eval.evaluations;
    report.h_columns = eval.applied_columns + report.start_h_columns;
    Ok(DirectMinState {
        orbitals: c,
        ritz_values: ritz,
        density,
        buffer_guess,
        report,
    })
}

/// trial.c = orth(C + t D), recording t.
fn retract(
    dense: &mut Dense,
    c: &RealMatrix,
    d: &RealMatrix,
    t: f64,
    y: &mut RealMatrix,
    trial: &mut Trial,
) -> Result<(), ScfError> {
    for ((y, c), d) in y
        .as_mut_slice()
        .iter_mut()
        .zip(c.as_slice())
        .zip(d.as_slice())
    {
        *y = c + t * d;
    }
    dense.lowdin(y, &mut trial.c)?;
    trial.step = t;
    Ok(())
}
