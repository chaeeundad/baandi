//! Default-off SCF bridge to ORC. Derived from real-gamma-design sections 4–10.
//! Completion is the only output boundary; downstream physics stays complex.
use crate::{ScfError, ScfOptions};
use dft_basis::PlaneWaveBasis;
use dft_mpi::Collective;
use dft_numerics::{RealMatrix, linalg::ComplexMatrix, pool::WorkerPool};
#[cfg(test)]
use dft_physics::data::Complex64;
use dft_physics::{
    distributed_real_gamma::DistributedRealHamiltonian, hamiltonian::Hamiltonian,
    real_gamma::RealHamiltonian,
};
use dft_solvers::{
    BandTolerances, EigenSolution, RealLinearOperator, SolverError, real_davidson::RealDavidson,
};
use std::{cell::RefCell, sync::Arc};

pub use crate::model::GammaReal;

/// Per-k dispatch evidence for the latest solve, also available before completion.
#[derive(Clone, Debug, Default)]
pub struct GammaRealDiagnostics {
    pub gamma_real_used: bool,
    pub fallback_reason: Option<String>,
}
impl GammaRealDiagnostics {
    pub(crate) fn fallback(
        &mut self,
        option: GammaReal,
        reason: impl Into<String>,
    ) -> Result<(), ScfError> {
        self.gamma_real_used = false;
        let reason = reason.into();
        if option == GammaReal::On {
            return Err(ScfError::Numerical(format!(
                "gamma_real=on is not eligible: {reason}"
            )));
        }
        self.fallback_reason = Some(reason);
        Ok(())
    }
}

/// Owned by one KState. The immutable Hamiltonian token guards physical data;
/// the key guards scratch sizing and solver options. Per-iteration tolerances
/// and certification policy are deliberately supplied afresh, not cached.
pub(crate) struct RealSolveState {
    /// Serial real Gamma operator. ISSUE-20 (frozen lazy construction): None
    /// on the M2 path until a fallback applies it (`ensure_serial_op`);
    /// always Some otherwise.
    op: Option<Operator>,
    /// ISSUE-20: reciprocal-layout descriptor (ORC half sphere) of the state,
    /// equal to the serial operator's when that exists.
    half: Arc<dft_basis::HalfSphere>,
    /// Pool, operator identity and grid the state was built for (reuse check
    /// of a state without a serial operator).
    pool: Arc<WorkerPool>,
    source_identity: Arc<()>,
    grid_dims: [usize; 3],
    key: RealSolveKey,
    workspace: Option<dft_solvers::real_davidson::RealDavidsonWorkspace>,
    /// Stick-distributed Hamiltonian of this rank (MPI stage M2c), built only when
    /// `parallel::m2_enabled()` and a communicator with two or more ranks is active.
    pub(crate) dist: Option<DistributedRealHamiltonian>,
    /// Local rows (stick partition) of the eigenvectors of the latest solve when
    /// it ran on the M2 path; None after any other solve (M2d density/energy).
    pub(crate) m2_local: Option<RealMatrix>,
    pub(crate) m2_density: Option<crate::m2::M2Density>,
    /// NL-ENERGY-REAL: per-band <x|V_nl|x> (Ha) of the latest real solve.
    nl_values: Vec<f64>,
    nl_valid: bool,
}
impl RealSolveState {
    /// NL-ENERGY-REAL values of the latest solve of this state, if computed.
    pub(crate) fn nonlocal_partials(&self) -> Option<&[f64]> {
        self.nl_valid.then_some(self.nl_values.as_slice())
    }
}

thread_local! {
    static M2_DIST_DISABLED: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}
/// ISSUE-20 test hook: treat the distributed operator of the calling thread
/// as unavailable (its construction "fails"), so that the M1 fallback, and
/// with it the lazy serial construction, is exercised. Never set in
/// production code.
#[doc(hidden)]
pub fn set_m2_dist_disabled_for_tests(value: bool) {
    M2_DIST_DISABLED.with(|o| o.set(value));
}

#[cfg(test)]
use crate::refscf::lazy::assemble_lazy_block;
pub use crate::refscf::lazy::{LazyOrbitals, set_delivered_orbitals};

/// C3a: audits of the real nonlocal energy that read the solver's kept
/// projections (index 0) or formed P^T X directly (index 1).
static NL_AUDITS: [std::sync::atomic::AtomicU64; 2] =
    [const { std::sync::atomic::AtomicU64::new(0) }; 2];
fn count_nl_audit(cached: bool) {
    NL_AUDITS[usize::from(!cached)].fetch_add(1, std::sync::atomic::Ordering::Relaxed);
}
/// Historical cache-audit diagnostics: (audits from the kept projections, audits
/// with a direct P^T X product) of the threads-only real Gamma path.
pub fn nl_audit_counters() -> (u64, u64) {
    let g = |i: usize| NL_AUDITS[i].load(std::sync::atomic::Ordering::Relaxed);
    (g(0), g(1))
}

#[derive(PartialEq, Eq)]
struct RealSolveKey {
    bands: usize,
    capacity: usize,
    threads: usize,
    max_iterations: usize,
    subspace: usize,
    capacity_per_band: usize,
    seed: u64,
    option: GammaReal,
}

struct Operator {
    h: RefCell<RealHamiltonian>,
    n: usize,
    pool: Arc<WorkerPool>,
}
impl RealLinearOperator for Operator {
    fn dimension(&self) -> usize {
        self.n
    }
    fn worker_pool(&self) -> Option<Arc<WorkerPool>> {
        Some(self.pool.clone())
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), SolverError> {
        let mut input = RealMatrix::zeros(self.n, 1)?;
        input.as_mut_slice().copy_from_slice(x);
        let mut out = RealMatrix::zeros(self.n, 1)?;
        self.apply_block(&input, &mut out)?;
        y.copy_from_slice(out.as_slice());
        Ok(())
    }
    fn apply_block(&self, x: &RealMatrix, y: &mut RealMatrix) -> Result<(), SolverError> {
        self.h
            .borrow_mut()
            .apply_block(x, y)
            .map_err(|e| SolverError::Application(e.to_string()))
    }
    // Cached-projector certificate hooks.
    fn projector_count(&self) -> usize {
        self.h.borrow().projector_count()
    }
    fn copy_last_projections(&self, m: usize, out: &mut [f64]) -> bool {
        match self.h.borrow().last_projections(m) {
            Ok(p) if p.len() == out.len() => {
                out.copy_from_slice(p);
                true
            }
            _ => false,
        }
    }
    fn apply_block_with_projections(
        &self,
        x: &RealMatrix,
        projections: &[f64],
        y: &mut RealMatrix,
    ) -> Result<(), SolverError> {
        self.h
            .borrow_mut()
            .apply_block_with_projections(x, projections, y)
            .map_err(|e| SolverError::Application(e.to_string()))
    }
}

fn real_capacity(options: &ScfOptions, bands: usize, n: usize) -> usize {
    options
        .davidson_subspace
        .max(options.davidson_capacity_per_band.saturating_mul(bands))
        .min(n)
}
/// Reuse the state when its key, pool and immutable operator data match (the
/// local potential is refreshed), otherwise build it. Unchanged construction
/// logic of `solve_real_gamma`, shared with ATOMIC-REAL.
fn ensure_state(
    hamiltonian: &Hamiltonian,
    basis: &PlaneWaveBasis,
    bands: usize,
    capacity: usize,
    options: &ScfOptions,
    pool: &Arc<WorkerPool>,
    state: &mut Option<RealSolveState>,
) -> Result<(), dft_physics::RealGammaError> {
    ensure_state_with(
        hamiltonian,
        basis,
        bands,
        capacity,
        options,
        pool,
        state,
        true,
    )
}

/// ISSUE-20: serial real Hamiltonian constructions of this process (test and
/// diagnostic counter; `RealHamiltonian::from_complex` calls of solve states).
static SERIAL_BUILDS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
/// ISSUE-20: number of serial real Gamma operator constructions of solve
/// states so far in this process (counts only; never read by the numerics).
#[doc(hidden)]
pub fn serial_operator_builds() -> u64 {
    SERIAL_BUILDS.load(std::sync::atomic::Ordering::Relaxed)
}
fn build_operator(
    hamiltonian: &Hamiltonian,
    basis: &PlaneWaveBasis,
    capacity: usize,
    pool: &Arc<WorkerPool>,
) -> Result<Operator, dft_physics::RealGammaError> {
    SERIAL_BUILDS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let real = RealHamiltonian::from_complex(hamiltonian, basis, capacity, pool.clone())?;
    Ok(Operator {
        n: real.half_sphere().full_len(),
        h: RefCell::new(real),
        pool: pool.clone(),
    })
}

/// ISSUE-20: build the serial operator of a state that has none (lazy M2
/// state reached by a fallback), at `capacity`; a no-op when it exists. The
/// rest of the state (distributed operator, kept data) is untouched.
fn ensure_serial_op(
    hamiltonian: &Hamiltonian,
    basis: &PlaneWaveBasis,
    capacity: usize,
    state: &mut RealSolveState,
) -> Result<(), dft_physics::RealGammaError> {
    if state.op.is_none() {
        state.op = Some(build_operator(hamiltonian, basis, capacity, &state.pool)?);
    }
    Ok(())
}

/// `ensure_state` with the serial operator optional: `build_op = false`
/// (ISSUE-20 lazy M2 state) builds or reuses a state holding only the
/// half-sphere descriptor; a reused state keeps a serial operator it already
/// has (its potential is refreshed).
#[allow(clippy::too_many_arguments)]
fn ensure_state_with(
    hamiltonian: &Hamiltonian,
    basis: &PlaneWaveBasis,
    bands: usize,
    capacity: usize,
    options: &ScfOptions,
    pool: &Arc<WorkerPool>,
    state: &mut Option<RealSolveState>,
    build_op: bool,
) -> Result<(), dft_physics::RealGammaError> {
    let key = RealSolveKey {
        bands,
        capacity,
        threads: options.threads,
        max_iterations: options.davidson_max_iterations,
        subspace: options.davidson_subspace,
        capacity_per_band: options.davidson_capacity_per_band,
        seed: options.seed,
        option: options.gamma_real,
    };
    let reusable = state.as_mut().is_some_and(|state| {
        state.key == key
            && Arc::ptr_eq(&state.pool, pool)
            && match state.op.as_mut() {
                Some(op) => op.h.get_mut().update_from_complex(hamiltonian, basis),
                // The checks of `RealHamiltonian::update_from_complex` without
                // the potential copy (a descriptor-only state has none).
                None => {
                    Arc::ptr_eq(&state.source_identity, &hamiltonian.operator_identity())
                        && hamiltonian.basis_momenta() == basis.momenta()
                        && basis.k_fractional() == [0.; 3]
                        && state.grid_dims == hamiltonian.fft_grid().dims()
                }
            }
    });
    if reusable {
        if build_op {
            let state = state.as_mut().expect("reusable state");
            ensure_serial_op(hamiltonian, basis, capacity, state)?;
        }
        return Ok(());
    }
    let (op, half) = if build_op {
        let op = build_operator(hamiltonian, basis, capacity, pool)?;
        let half = Arc::new(op.h.borrow().half_sphere().clone());
        (Some(op), half)
    } else {
        // The construction checks of `RealHamiltonian::from_complex` that
        // precede its plans (shape, half sphere).
        if hamiltonian.basis_momenta() != basis.momenta() {
            return Err(dft_physics::RealGammaError::Shape);
        }
        let half = dft_basis::HalfSphere::new(basis, hamiltonian.fft_grid())?;
        (None, Arc::new(half))
    };
    *state = Some(RealSolveState {
        op,
        half,
        pool: pool.clone(),
        source_identity: hamiltonian.operator_identity(),
        grid_dims: hamiltonian.fft_grid().dims(),
        key,
        workspace: None,
        dist: None,
        m2_local: None,
        m2_density: None,
        nl_values: Vec::new(),
        nl_valid: false,
    });
    Ok(())
}

/// ATOMIC-REAL: H Q for the atomic-start Ritz step through the real Gamma
/// Hamiltonian (pair FFTs, real projector GEMM) instead of the complex apply.
/// The state is built exactly as `solve_real_gamma` builds it, so the first
/// solve reuses it. `None` (caller keeps the complex apply) unless the k state
/// is exactly Gamma, the real path is selected and eligible on the threads-only
/// path, Q packs to ORC within the existing tolerance and fits the capacity.
#[allow(clippy::too_many_arguments)]
pub(crate) fn atomic_apply_real(
    hamiltonian: &Hamiltonian,
    basis: &PlaneWaveBasis,
    q: &ComplexMatrix,
    bands: usize,
    options: &ScfOptions,
    pool: &Arc<WorkerPool>,
    exclusion: Option<&str>,
    state: &mut Option<RealSolveState>,
) -> Result<Option<ComplexMatrix>, ScfError> {
    if options.gamma_real == GammaReal::Off
        || exclusion.is_some()
        || basis.k_fractional() != [0.; 3]
        || options.locking != dft_solvers::locking::Locking::Off
        || options.spectral_guard != crate::SpectralGuard::Off
        || !matches!(options.eigensolver, dft_solvers::SolverChoice::Davidson)
        || crate::parallel::SharedComm::active().is_some()
    {
        return Ok(None);
    }
    let n = q.rows();
    let capacity = real_capacity(options, bands, hamiltonian.kinetic_diagonal().len());
    if q.cols() > capacity
        || ensure_state(hamiltonian, basis, bands, capacity, options, pool, state).is_err()
    {
        return Ok(None);
    }
    let state = state.as_mut().unwrap();
    let Some(op) = state.op.as_ref() else {
        return Ok(None);
    };
    let mut h = op.h.borrow_mut();
    if h.half_sphere().full_len() != n {
        return Ok(None);
    }
    let mut packed = RealMatrix::zeros(n, q.cols()).map_err(SolverError::from)?;
    let half = h.half_sphere();
    if pool
        .install(|| {
            pair_columns(q.as_slice(), packed.as_mut_slice(), n, &|src, dst| {
                half.pack(src, dst)
            })
        })
        .is_err()
    {
        return Ok(None);
    }
    let mut image = RealMatrix::zeros(n, q.cols()).map_err(SolverError::from)?;
    h.apply_block(&packed, &mut image)
        .map_err(|e| ScfError::Numerical(e.to_string()))?;
    let mut hq = ComplexMatrix::zeros(n, q.cols()).map_err(SolverError::from)?;
    let half = h.half_sphere();
    pool.install(|| {
        pair_columns(image.as_slice(), hq.as_mut_slice(), n, &|src, dst| {
            half.unpack(src, dst)
        })
    })
    .map_err(|e| ScfError::Numerical(e.to_string()))?;
    Ok(Some(hq))
}

/// Try the real path without modifying an incompatible warm start. Construction
/// failures are eligibility failures; solver failures propagate, never fall back.
/// Certification follows the same caller policy as complex Davidson.
#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_real_gamma(
    hamiltonian: &Hamiltonian,
    basis: &PlaneWaveBasis,
    initial: Option<&ComplexMatrix>,
    initial_lazy: Option<&LazyOrbitals>,
    kinetic: &[f64],
    bands: usize,
    tolerance: f64,
    policy: BandTolerances,
    occupations: &[f64],
    electrons: f64,
    certify: bool,
    options: &ScfOptions,
    pool: Arc<WorkerPool>,
    exclusion: Option<&str>,
    diagnostics: &mut GammaRealDiagnostics,
    state: &mut Option<RealSolveState>,
) -> Result<Option<EigenSolution>, ScfError> {
    let option = options.gamma_real;
    let comm = crate::parallel::SharedComm::active();
    let reason = if option == GammaReal::Off {
        Some("option off")
    } else if let Some(reason) = exclusion {
        Some(reason)
    } else if basis.k_fractional() != [0.; 3] {
        Some("k is not exactly Gamma")
    } else if options.locking != dft_solvers::locking::Locking::Off {
        Some("locking enabled")
    } else if options.spectral_guard != crate::SpectralGuard::Off {
        Some("spectral guard enabled")
    } else if !matches!(options.eigensolver, dft_solvers::SolverChoice::Davidson) {
        Some("hybrid/RMM eigensolver")
    } else if comm.is_some()
        && options.preconditioner != dft_solvers::DavidsonPreconditioner::Default
    {
        Some("distributed Davidson supports the default preconditioner only")
    } else {
        None
    };
    if let Some(reason) = reason {
        *state = None;
        if diagnostics.fallback_reason.as_deref() != Some(reason) {
            diagnostics.fallback(option, reason)?;
        }
        return Ok(None);
    }
    let capacity = real_capacity(options, bands, kinetic.len());
    let t_state = dft_numerics::stage_clock::start();
    // ISSUE-20 (frozen lazy construction): on the M2 path the state holds
    // only the half-sphere descriptor; the serial operator is built below
    // only when the distributed operator is unavailable. Replicated condition.
    let lazy_serial = comm.is_some() && crate::parallel::m2_enabled();
    let ensured = ensure_state_with(
        hamiltonian,
        basis,
        bands,
        capacity,
        options,
        &pool,
        state,
        !lazy_serial,
    );
    crate::m2::add(crate::m2::Stage::RealState, t_state);
    if let Err(e) = ensured {
        *state = None;
        diagnostics.fallback(option, format!("real Hamiltonian construction: {e}"))?;
        return Ok(None);
    }
    let state = state.as_mut().unwrap();
    // BREAKTHROUGH-4 (`BAANDI_PHASE_WAIT=on`, diagnostic): the arrival skew of
    // the ranks at the solve entry, apart from the solve's own work. The
    // condition is replicated (every rank reaches this point for this map).
    if let Some(c) = comm.as_ref().filter(|_| {
        crate::parallel::m2_enabled() && dft_solvers::distributed_davidson::wait_probe()
    }) {
        let t = dft_numerics::stage_clock::start();
        c.barrier()
            .map_err(|e| ScfError::Numerical(format!("entry wait: {e}")))?;
        crate::m2::add(crate::m2::Stage::EntryWait, t);
    }
    state.nl_valid = false;
    // The strict solve constructs a fresh distributed start.
    state.m2_local = None;
    let t_prepare = dft_numerics::stage_clock::start();
    // M2c: the distributed Hamiltonian follows the SCF potential like the serial one.
    match comm.as_ref().filter(|_| crate::parallel::m2_enabled()) {
        Some(c) => {
            let fresh = state.dist.as_mut().is_some_and(|d| {
                let fresh = d.update_from_complex(hamiltonian, basis);
                if fresh {
                    use crate::grid_dist as g;
                    g::owned_add(g::O_FULL_COPIES, 1);
                    g::owned_add(g::O_FULL_COPY_VALUES, d.potential().len() as u64);
                }
                fresh
            });
            if !fresh {
                state.dist = if M2_DIST_DISABLED.with(|o| o.get()) {
                    None
                } else {
                    DistributedRealHamiltonian::new(
                        hamiltonian,
                        basis,
                        capacity,
                        c.rank(),
                        c.size(),
                        pool.clone(),
                    )
                    .ok()
                };
                state.m2_density = None;
                // A rebuilt operator (atoms, projectors, basis, grid): nothing kept.
            }
            if let Some(d) = state.dist.as_mut() {
                d.set_batched(crate::m2::batched_transposes());
            }
        }
        None => {
            state.dist = None;
            state.m2_density = None;
        }
    }

    if state.dist.is_none()
        && let Err(e) = ensure_serial_op(hamiltonian, basis, capacity, state)
    {
        // ISSUE-20: a lazy M2 state reached the M1 fallback; the serial
        // operator is built now (same construction, same failure policy).
        diagnostics.fallback(option, format!("real Hamiltonian construction: {e}"))?;
        return Ok(None);
    }
    let real = state.op.as_ref().map(|op| op.h.borrow());
    let half: &dft_basis::HalfSphere = &state.half;
    let n = half.full_len();
    let mut diagonal = vec![0.; n];
    diagonal[0] = kinetic[0];
    for (i, &p) in half.plus().iter().enumerate().skip(1) {
        diagonal[2 * i - 1] = kinetic[p as usize];
        diagonal[2 * i] = kinetic[p as usize];
    }
    // Snapshot replay (BREAKTHROUGH-0-INSTRUMENT): an injected start (this
    // rank's rows on the M2 path, the ORC block on the threads path) and its
    // kept H X replace the placeholder `initial` of `ScfProblem::replay_solve`.
    let mut injected = crate::snapshot::take_injection();
    // ISSUE-18: check the injected start against the captured contract
    // (rows, dimensions, potential, stopping profile, kept-product status)
    // before the solver sees it; the verdict is agreed over all ranks (one
    // tree sum of the mismatch classes; every replaying rank reaches this
    // point) so that no rank enters a divergent solve.
    if let Some(inj) = injected.as_mut() {
        let contract = crate::snapshot::take_contract();
        let mode = contract
            .as_ref()
            .map_or(crate::snapshot::ReplayMode::Legacy, |c| c.mode);
        let stick = comm.is_some() && state.dist.is_some();
        let (verdict, notes) = match &contract {
            Some(c) => {
                let d = state.dist.as_ref().filter(|_| stick);
                let v = crate::snapshot::check_solve(
                    c,
                    &crate::snapshot::SolveIdentity {
                        layout: if stick {
                            crate::snapshot::Layout::StickRows
                        } else {
                            crate::snapshot::Layout::OrcFull
                        },
                        rows: d.map(|d| d.local_rows()),
                        n_global: n,
                        bands,
                        x_rows: inj.x.rows(),
                        x_cols: inj.x.cols(),
                        has_hx: inj.hx.is_some(),
                        potential: d.map_or_else(
                            || {
                                real.as_ref()
                                    .expect("serial operator off the M2 path")
                                    .potential()
                            },
                            |d| d.potential(),
                        ),
                        tolerance,
                        strict_bands: policy.strict_bands,
                        buffer_tolerance: policy.buffer_tolerance,
                        ritz: false,
                    },
                );
                let notes = v.notes.clone();
                (Some(v), notes)
            }
            None => (None, vec!["legacy replay_solve: unchecked".to_string()]),
        };
        let mut flags: Vec<f64> = verdict.as_ref().map_or(vec![0.0; 5], |v| {
            v.errors.iter().map(|&e| f64::from(u8::from(e))).collect()
        });
        if let Some(c) = comm.as_ref() {
            let mut out = vec![0.0; flags.len()];
            c.tree_sum_f64(&flags, &mut out)
                .map_err(|e| ScfError::Numerical(format!("replay identity agreement: {e}")))?;
            flags = out;
        }
        if flags.iter().any(|&f| f > 0.0) {
            let classes: Vec<&str> = crate::snapshot::CHECK_CLASSES
                .iter()
                .zip(&flags)
                .filter(|(_, f)| **f > 0.0)
                .map(|(c, _)| *c)
                .collect();
            let local = verdict.map(|v| v.messages).unwrap_or_default();
            return Err(ScfError::Numerical(format!(
                "snapshot replay ({}): identity mismatch [{}] (classes agreed over ranks): {}",
                mode.name(),
                classes.join(", "),
                if local.is_empty() {
                    "on another rank".to_string()
                } else {
                    local.join("; ")
                }
            )));
        }
        let recompute = verdict.as_ref().is_some_and(|v| v.drop_hx);
        if recompute {
            inj.hx = None;
        }
        crate::snapshot::set_report(Some(crate::snapshot::ReplayReport {
            mode,
            hx_injected: inj.hx.is_some(),
            hx_recomputed: recompute,
            notes,
        }));
    }
    let mut injected_packed: Option<RealMatrix> = None;
    let mut injected_hx: Option<RealMatrix> = None;
    let initial = match (initial, initial_lazy, injected) {
        (Some(m), _, Some(inj)) if m.rows() == 0 => {
            injected_packed = Some(inj.x);
            injected_hx = inj.hx;
            None
        }
        (Some(m), Some(h), _) if m.rows() == 0 => Some(h.materialize()?),
        (Some(m), None, _) if m.rows() == 0 => {
            return Err(ScfError::InvalidInput(
                "lazy warm start without its distributed rows",
            ));
        }
        (m, _, _) => m,
    };
    let packed = if let Some(p) = injected_packed.take() {
        Some(p)
    } else if let Some(initial) = initial {
        let mut packed = RealMatrix::zeros(n, initial.cols()).map_err(SolverError::from)?;
        let converted = pool.install(|| {
            pair_columns(initial.as_slice(), packed.as_mut_slice(), n, &|src, dst| {
                half.pack(src, dst)
            })
        });
        if let Err(e) = converted {
            diagnostics.fallback(
                option,
                format!("warm start is not ORC-compatible (1e-14 absolute): {e}"),
            )?;
            return Ok(None);
        }
        Some(packed)
    } else {
        None
    };
    if state.dist.is_some() {
        crate::m2::add(crate::m2::Stage::SolvePrepare, t_prepare);
    }
    let solver = RealDavidson {
        bands,
        tolerance,
        max_iterations: options.davidson_max_iterations,
        max_subspace: capacity,
        seed: options.seed,
    };
    drop(real);
    // Snapshot replay on the threads path: the captured kept X and H X.
    if comm.is_none()
        && let (Some(x), Some(hx)) = (packed.as_ref(), injected_hx.take())
    {
        match state.workspace.as_mut() {
            Some(ws) => ws.load_kept_start(x, &hx)?,
            None => {
                // ISSUE-18: a checked replay must not silently drop the
                // captured kept-product status (threads path: no collective).
                if let Some(r) = crate::snapshot::take_report() {
                    if r.mode == crate::snapshot::ReplayMode::Exact {
                        return Err(ScfError::Numerical(
                            "snapshot replay (exact): kept-product status: no solver workspace (REUSE_PSI off), captured H X cannot be used".into(),
                        ));
                    }
                    let mut r = r;
                    r.hx_injected = false;
                    r.hx_recomputed = true;
                    r.notes.push("no solver workspace: H X recomputed".into());
                    crate::snapshot::set_report(Some(r));
                }
                eprintln!("snapshot replay: no solver workspace (REUSE_PSI off): kept H X not used")
            }
        }
    }
    // Snapshot capture of a threads-only solve (M2 solves are handled below, where
    // their rows and kept H X are known).
    let occupied = occupations.iter().filter(|&&f| f > 0.0).count();
    let criterion = crate::snapshot::Criterion::ResidualNorm;
    let capacity_per_band = options.davidson_capacity_per_band;
    let threads_capture = if comm.is_none() {
        let ws = state.workspace.as_ref();
        let armed = ws.is_some_and(|w| w.kept_start_armed());
        let kept = ws.and_then(|w| w.kept_start());
        let hx = kept.map(|(_, hx)| hx);
        // The basis of a taken kept start is the kept X itself, not the
        // round-tripped `packed` start: capture and analyse that block.
        let x_block: Option<&RealMatrix> = if armed {
            kept.map(|(x, _)| x)
        } else {
            packed.as_ref()
        };
        crate::snapshot::capture_entry(crate::snapshot::Entry {
            rank: 0,
            ranks: 1,
            layout: crate::snapshot::Layout::OrcFull,
            x: x_block,
            hx,
            armed,
            rows_fnv: 0,
            n_global: n,
            bands,
            occupied,
            occupations,
            criterion,
            tolerance,
            strict_bands: policy.strict_bands,
            buffer_tolerance: policy.buffer_tolerance,
            ritz: false,
            potential: state
                .op
                .as_ref()
                .expect("serial operator on the threads path")
                .h
                .borrow()
                .potential(),
            grid: state
                .op
                .as_ref()
                .expect("serial operator on the threads path")
                .h
                .borrow()
                .grid_dims(),
            cutoff_hartree: options.cutoff_hartree,
            electrons,
            capacity_per_band,
        })
    } else {
        None
    };
    // None only on the M2 path of a lazy state (ISSUE-20); every branch that
    // applies the serial operator runs with `dist` None, where it was built.
    let op_opt = state.op.as_ref();
    let dist = state.dist.as_mut();
    let result = if let (Some(comm), Some(dist)) = (&comm, dist) {
        crate::parallel::count_distributed_solve();
        crate::parallel::count_distributed_rows_solve();
        let local_diagonal: Vec<f64> = dist.local_rows().iter().map(|&r| diagonal[r]).collect();
        {
            // Strict (residual) solves and the allocating path: capture only
            // (no kept start exists here).
            let start_rows = match packed.as_ref() {
                Some(p) => Some(
                    dist.extract_local(p)
                        .map_err(|e| ScfError::Numerical(e.to_string()))?,
                ),
                None => None,
            };
            let capture = crate::snapshot::capture_entry(crate::snapshot::Entry {
                rank: comm.rank(),
                ranks: comm.size(),
                layout: crate::snapshot::Layout::StickRows,
                x: start_rows.as_ref(),
                hx: None,
                armed: false,
                rows_fnv: crate::snapshot::fnv_usize(dist.local_rows()),
                n_global: n,
                bands,
                occupied,
                occupations,
                criterion,
                tolerance,
                strict_bands: policy.strict_bands,
                buffer_tolerance: policy.buffer_tolerance,
                ritz: false,
                potential: dist.potential(),
                grid: dist.grid_shape(),
                cutoff_hartree: options.cutoff_hartree,
                electrons,
                capacity_per_band,
            });
            let result = solver.solve_distributed_rows(
                comm,
                dist,
                packed.as_ref(),
                &local_diagonal,
                policy,
                false,
                certify,
                options.davidson_capacity_per_band,
            )?;
            if let Some(c) = capture {
                crate::snapshot::capture_exit(c, &result, false);
            }
            crate::profile::note_solve(result.restarts, result.contracted_restarts, false);
            // M2d: keep this rank's rows for the distributed density and energy.
            state.m2_local = Some(
                dist.extract_local(&result.eigenvectors)
                    .map_err(|e| ScfError::Numerical(e.to_string()))?,
            );
            result
        }
    } else if let Some(comm) = &comm {
        let op = op_opt.ok_or(ScfError::InvalidInput("serial real operator missing"))?;
        crate::parallel::count_distributed_solve();
        solver.solve_distributed(
            comm,
            op,
            packed.as_ref(),
            &diagonal,
            policy,
            false,
            certify,
            options.davidson_capacity_per_band,
        )?
    } else {
        let op = op_opt.ok_or(ScfError::InvalidInput("serial real operator missing"))?;
        let armed_before = state
            .workspace
            .as_ref()
            .is_some_and(|w| w.kept_start_armed());
        let result = solver.solve_reusing_with_preconditioner(
            op,
            packed.as_ref(),
            &diagonal,
            policy,
            false,
            certify,
            options.davidson_capacity_per_band,
            &mut state.workspace,
            options.preconditioner,
        )?;
        // An armed kept start that is no longer armed after the solve was
        // consumed (taken, or declined and dropped): the lean counters tell
        // which; here the taken flag is the armed-before flag.
        let kept_taken = armed_before
            && !state
                .workspace
                .as_ref()
                .is_some_and(|w| w.kept_start_armed());
        if let Some(c) = threads_capture {
            crate::snapshot::capture_exit(c, &result, kept_taken);
        }
        crate::profile::note_solve(result.restarts, result.contracted_restarts, kept_taken);
        result
    };
    // NL-ENERGY-REAL: nonlocal expectation values of the returned real vectors
    // for the energy audit of this map (threads-only path).
    if comm.is_none()
        && let Some(op) = op_opt
    {
        state.nl_values.resize(bands, 0.0);
        count_nl_audit(false);
        state.nl_valid =
            op.h.borrow_mut()
                .nonlocal_expectation_values(&result.eigenvectors, &mut state.nl_values)
                .is_ok();
    }

    let t_unpack = dft_numerics::stage_clock::start();
    let mut eigenvectors = ComplexMatrix::zeros(n, bands).map_err(SolverError::from)?;
    {
        let half: &dft_basis::HalfSphere = &state.half;
        pool.install(|| {
            pair_columns(
                result.eigenvectors.as_slice(),
                eigenvectors.as_mut_slice(),
                n,
                &|src, dst| half.unpack(src, dst),
            )
        })
        .map_err(|e| ScfError::Numerical(e.to_string()))?;
    }
    if state.dist.is_some() {
        crate::m2::add(crate::m2::Stage::SolveUnpack, t_unpack);
    }
    diagnostics.gamma_real_used = true;
    diagnostics.fallback_reason = None;
    Ok(Some(EigenSolution {
        locking: None,
        spectral_completeness: None,
        eigenvalues: result.eigenvalues,
        eigenvectors,
        residuals: result.residuals,
        converged: result.converged,
        buffer: result.buffer,
        residuals_certified: result.residuals_certified,
        iterations: result.iterations,
        history: result.history,
        operator_applications: result.operator_applications,
    }))
}

/// Test adapter for packing selected real-Gamma rows.
#[cfg(test)]
fn pack_orbital_rows(
    half: &dft_basis::HalfSphere,
    full: &[Complex64],
    rows: &[usize],
    out: &mut [f64],
) -> Result<(), dft_basis::HalfSphereError> {
    use dft_basis::{HalfSphereError as E, ORC_TOLERANCE};
    use std::f64::consts::SQRT_2;
    if full.len() != half.full_len()
        || out.len() != rows.len()
        || rows.iter().any(|&r| r >= half.full_len())
    {
        return Err(E::Dimension);
    }
    if full.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
        || full[0].im.abs() > ORC_TOLERANCE
    {
        return Err(E::InvalidCoefficients);
    }
    for (&p, &m) in half.plus().iter().zip(half.minus()).skip(1) {
        let z = full[p as usize];
        if (z - full[m as usize].conj()).norm() > ORC_TOLERANCE
            || !(z.re * SQRT_2).is_finite()
            || !(z.im * SQRT_2).is_finite()
        {
            return Err(E::InvalidCoefficients);
        }
    }
    for (out, &r) in out.iter_mut().zip(rows) {
        *out = if r == 0 {
            full[0].re
        } else {
            let z = full[half.plus()[r.div_ceil(2)] as usize];
            SQRT_2 * if r % 2 == 1 { z.re } else { z.im }
        };
    }
    Ok(())
}

fn pair_columns<T: Sync, U: Send, E: Send>(
    input: &[T],
    output: &mut [U],
    n: usize,
    convert: &(impl Fn(&[T], &mut [U]) -> Result<(), E> + Sync),
) -> Result<(), E> {
    let pairs = (input.len() / n).div_ceil(2);
    if pairs > 1 {
        let split = (pairs / 2) * 2 * n;
        let (li, ri) = input.split_at(split);
        let (lo, ro) = output.split_at_mut(split);
        let (a, b) = WorkerPool::join(
            || pair_columns(li, lo, n, convert),
            || pair_columns(ri, ro, n, convert),
        );
        a?;
        b
    } else {
        for (src, dst) in input.chunks(n).zip(output.chunks_mut(n)) {
            convert(src, dst)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod atomic_real_tests {
    use super::*;
    use dft_geometry::Cell;
    use dft_physics::data::Array3;
    use dft_physics::data::Complex64;
    use dft_physics::hamiltonian::NonlocalOperator;

    /// ATOMIC-REAL: H Q through the real Hamiltonian equals the complex apply to
    /// rounding for Hermitian-completed Q, leaves a state the solve reuses, and
    /// declines non-Hermitian input and the off option.
    #[test]
    fn atomic_apply_real_matches_complex() {
        let cell =
            Cell::from_user_rows_bohr([[4.2, 0.3, 0.], [0., 3.9, 0.2], [0.1, 0., 4.4]]).unwrap();
        let options = ScfOptions {
            gamma_real: GammaReal::Auto,
            cutoff_hartree: 2.5,
            ..ScfOptions::strict_reference()
        };
        let basis = options.basis(&cell, [0.; 3], 1.).unwrap();
        let grid = dft_basis::FftGridSpec::new(&basis, &basis.density_basis(&cell).unwrap(), None)
            .unwrap();
        let potential = Array3::from_shape_fn(grid.dims(), |(i, j, k)| {
            0.1 * ((i + 2 * j + 3 * k) as f64).sin()
        });
        let h = Hamiltonian::new(
            &cell,
            &basis,
            potential,
            NonlocalOperator::new(&basis, cell.volume_bohr3(), &[]).unwrap(),
        )
        .unwrap();
        let half = dft_basis::HalfSphere::new(&basis, &grid).unwrap();
        let n = basis.vectors().len();
        let cols = 5;
        let mut q = ComplexMatrix::zeros(n, cols).unwrap();
        for j in 0..cols {
            let real: Vec<f64> = (0..n).map(|i| ((i * 3 + j * 7 + 1) as f64).sin()).collect();
            half.unpack(&real, q.column_mut(j)).unwrap();
        }
        let mut expected = ComplexMatrix::zeros(n, cols).unwrap();
        dft_solvers::LinearOperator::apply_block(&h, &q, &mut expected).unwrap();
        let pool = WorkerPool::new(2).unwrap();
        let mut state = None;
        let hq = atomic_apply_real(&h, &basis, &q, 4, &options, &pool, None, &mut state)
            .unwrap()
            .expect("eligible");
        let scale = expected
            .as_slice()
            .iter()
            .map(|z| z.norm())
            .fold(0., f64::max);
        let worst = hq
            .as_slice()
            .iter()
            .zip(expected.as_slice())
            .map(|(a, b)| (a - b).norm())
            .fold(0., f64::max);
        assert!(worst <= 1e-13 * scale, "{worst:e} of {scale:e}");
        assert!(state.is_some());
        // The solve's own construction check reuses the state built here.
        let before = state
            .as_ref()
            .map(|s| s.op.as_ref().unwrap().h.as_ptr() as usize);
        let capacity = real_capacity(&options, 4, h.kinetic_diagonal().len());
        ensure_state(&h, &basis, 4, capacity, &options, &pool, &mut state).unwrap();
        assert_eq!(
            before,
            state
                .as_ref()
                .map(|s| s.op.as_ref().unwrap().h.as_ptr() as usize)
        );
        // Non-Hermitian input and the off option fall back to the complex apply.
        let mut bad = q.clone();
        bad.column_mut(0)[1] += Complex64::new(0., 1e-6);
        assert!(
            atomic_apply_real(&h, &basis, &bad, 4, &options, &pool, None, &mut state)
                .unwrap()
                .is_none()
        );
        let off = ScfOptions {
            gamma_real: GammaReal::Off,
            ..options.clone()
        };
        assert!(
            atomic_apply_real(&h, &basis, &q, 4, &off, &pool, None, &mut state)
                .unwrap()
                .is_none()
        );
    }
}

#[cfg(test)]
mod lazy_serial_tests {
    use super::*;
    use dft_geometry::Cell;
    use dft_physics::data::Array3;
    use dft_physics::hamiltonian::NonlocalOperator;

    /// ISSUE-20: a descriptor-only state (no serial operator) is reused across
    /// a potential update of the same operator, rebuilt when the operator,
    /// bands, options, capacity or pool change, and gets its serial operator
    /// lazily (in place, descriptor kept) when a caller needs it; the lazily
    /// built operator equals an eagerly built one.
    #[test]
    fn descriptor_state_reuse_invalidation_and_lazy_operator() {
        let cell =
            Cell::from_user_rows_bohr([[4.2, 0.3, 0.], [0., 3.9, 0.2], [0.1, 0., 4.4]]).unwrap();
        let mut options = ScfOptions {
            gamma_real: GammaReal::Auto,
            cutoff_hartree: 2.5,
            ..ScfOptions::strict_reference()
        };
        let basis = options.basis(&cell, [0.; 3], 1.).unwrap();
        let grid = dft_basis::FftGridSpec::new(&basis, &basis.density_basis(&cell).unwrap(), None)
            .unwrap();
        let make = || {
            Hamiltonian::new(
                &cell,
                &basis,
                Array3::zeros(grid.dims()),
                NonlocalOperator::new(&basis, cell.volume_bohr3(), &[]).unwrap(),
            )
            .unwrap()
        };
        let mut h = make();
        let mut pool = WorkerPool::new(1).unwrap();
        let mut state = None;
        let half_ptr = |s: &Option<RealSolveState>| s.as_ref().unwrap().half.plus().as_ptr();
        let has_op = |s: &Option<RealSolveState>| s.as_ref().unwrap().op.is_some();
        let n_pw = h.kinetic_diagonal().len();
        let cap = |o: &ScfOptions, bands| real_capacity(o, bands, n_pw);
        // (The process-wide build counter is checked by tests/m2_lazy_serial.rs;
        // other unit tests of this binary build operators concurrently.)
        let c = cap(&options, 3);
        ensure_state_with(&h, &basis, 3, c, &options, &pool, &mut state, false).unwrap();
        assert!(!has_op(&state), "descriptor only");
        let reference = dft_basis::HalfSphere::new(&basis, h.fft_grid()).unwrap();
        assert_eq!(state.as_ref().unwrap().half.plus(), reference.plus());
        assert_eq!(state.as_ref().unwrap().half.minus(), reference.minus());
        let mut previous = half_ptr(&state);
        // Potential update of the same operator: reused.
        h.set_local_potential(Array3::from_elem(grid.dims(), 0.25))
            .unwrap();
        ensure_state_with(&h, &basis, 3, c, &options, &pool, &mut state, false).unwrap();
        assert_eq!(half_ptr(&state), previous, "potential update reuses");
        // Lazy serial operator: built in place, descriptor kept, potential current.
        ensure_serial_op(&h, &basis, c, state.as_mut().unwrap()).unwrap();
        assert!(has_op(&state));
        assert_eq!(half_ptr(&state), previous, "descriptor kept");
        {
            let s = state.as_ref().unwrap();
            let lazy = s.op.as_ref().unwrap().h.borrow();
            assert_eq!(lazy.potential(), h.local_potential_values());
            assert_eq!(lazy.half_sphere().plus(), s.half.plus());
            // Same operator as an eager construction: identical images.
            let mut eager = None;
            ensure_state(&h, &basis, 3, c, &options, &pool, &mut eager).unwrap();
            let eager = eager.unwrap();
            let n = s.half.full_len();
            let x = RealMatrix::new(
                n,
                2,
                (0..2 * n).map(|i| ((i * 7 + 3) as f64).sin()).collect(),
            )
            .unwrap();
            let mut y_lazy = RealMatrix::zeros(n, 2).unwrap();
            let mut y_eager = RealMatrix::zeros(n, 2).unwrap();
            drop(lazy);
            s.op.as_ref().unwrap().apply_block(&x, &mut y_lazy).unwrap();
            eager
                .op
                .as_ref()
                .unwrap()
                .apply_block(&x, &mut y_eager)
                .unwrap();
            let bits =
                |m: &RealMatrix| m.as_slice().iter().map(|v| v.to_bits()).collect::<Vec<_>>();
            assert_eq!(bits(&y_lazy), bits(&y_eager));
        }
        // A state with a serial operator stays reusable for descriptor-only calls.
        ensure_state_with(&h, &basis, 3, c, &options, &pool, &mut state, false).unwrap();
        assert!(has_op(&state));
        assert_eq!(half_ptr(&state), previous);
        // Invalidation: each change rebuilds a descriptor-only state.
        for step in 0..4 {
            let bands = if step == 0 { 4 } else { 3 };
            match step {
                1 => options.seed += 1,
                2 => pool = WorkerPool::new(2).unwrap(),
                3 => h = make(),
                _ => {}
            }
            let c = cap(&options, bands);
            ensure_state_with(&h, &basis, bands, c, &options, &pool, &mut state, false).unwrap();
            assert_ne!(half_ptr(&state), previous, "step {step} rebuilds");
            assert!(!has_op(&state), "step {step}: descriptor only");
            previous = half_ptr(&state);
        }
        // Capacity change (key) rebuilds as well.
        let c2 = c - 1;
        ensure_state_with(&h, &basis, 3, c2, &options, &pool, &mut state, false).unwrap();
        assert_ne!(half_ptr(&state), previous);
        assert!(!has_op(&state));
    }
}

#[cfg(test)]
mod lazy_assemble_tests {
    use super::*;
    use dft_geometry::Cell;

    /// The former assembly: scatter into a global ORC block, then unpack column
    /// by column.
    fn assemble_reference(
        half: &dft_basis::HalfSphere,
        counts: &[usize],
        index: &[f64],
        flat: &[f64],
        n: usize,
        bands: usize,
    ) -> Result<ComplexMatrix, ScfError> {
        let mut real = vec![0.0; n * bands];
        let (mut at, mut offset) = (0, 0);
        for &c in counts {
            let rows = &index[offset..offset + c];
            for j in 0..bands {
                let src = &flat[at + j * c..at + (j + 1) * c];
                for (&g, &v) in rows.iter().zip(src) {
                    real[j * n + g as usize] = v;
                }
            }
            at += c * bands;
            offset += c;
        }
        let mut full = ComplexMatrix::zeros(n, bands).map_err(SolverError::from)?;
        for j in 0..bands {
            half.unpack(&real[j * n..(j + 1) * n], full.column_mut(j))
                .map_err(|e| ScfError::Numerical(format!("lazy orbitals: {e}")))?;
        }
        Ok(full)
    }

    /// LAZY-ASSEMBLE-FUSED: bit-identical blocks for cubic, skewed and
    /// elongated cells, 1 to 5 rank segments of uneven sizes with permuted
    /// rows, signed zeros and tiny values; the same errors for a nonfinite
    /// value and a wrong length.
    #[test]
    fn fused_assembly_matches_scatter_then_unpack() {
        for rows in [
            [[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]],
            [[6., 1., 0.], [0., 5.5, 1.2], [0.7, 0., 6.5]],
            [[4., 0., 0.], [0., 4.5, 0.], [0., 0., 14.]],
        ] {
            let cell = Cell::from_user_rows_bohr(rows).unwrap();
            let basis = dft_basis::PlaneWaveBasis::new(&cell, dft_units::Hartree(2.)).unwrap();
            let grid = dft_basis::FftGridSpec::new(&basis, &basis, None).unwrap();
            let half = dft_basis::HalfSphere::new(&basis, &grid).unwrap();
            let n = half.full_len();
            let bands = 5;
            // A fixed permutation of the ORC coordinates (multiplicative step).
            let step = (1..n).rev().find(|s| gcd(*s, n) == 1).unwrap_or(1);
            let order: Vec<usize> = (0..n).map(|k| (k * step + 3) % n).collect();
            for parts in 1..=5usize {
                let mut counts = vec![n / parts; parts];
                counts[0] += n % parts;
                if parts > 2 {
                    // Uneven: move rows from the last segment to the second.
                    let d = counts[parts - 1] / 3;
                    counts[parts - 1] -= d;
                    counts[1] += d;
                }
                let index: Vec<f64> = order.iter().map(|&g| g as f64).collect();
                let mut flat = Vec::with_capacity(n * bands);
                let mut offset = 0;
                for &c in &counts {
                    for j in 0..bands {
                        for (k, &g) in order[offset..offset + c].iter().enumerate() {
                            let v = ((g * 13 + j * 7 + k) as f64).sin() * 1e-3;
                            flat.push(match (g + j) % 11 {
                                0 => -0.0,
                                1 => 1e-310,
                                _ => v,
                            });
                        }
                    }
                    offset += c;
                }
                let a = assemble_lazy_block(&half, &counts, &index, &flat, n, bands).unwrap();
                let b = assemble_reference(&half, &counts, &index, &flat, n, bands).unwrap();
                for (x, y) in a.as_slice().iter().zip(b.as_slice()) {
                    assert_eq!(x.re.to_bits(), y.re.to_bits());
                    assert_eq!(x.im.to_bits(), y.im.to_bits());
                }
                let mut bad = flat.clone();
                bad[flat.len() / 2] = f64::NAN;
                let ea = assemble_lazy_block(&half, &counts, &index, &bad, n, bands).unwrap_err();
                let eb = assemble_reference(&half, &counts, &index, &bad, n, bands).unwrap_err();
                assert_eq!(ea.to_string(), eb.to_string());
            }
            let e = assemble_lazy_block(&half, &[n - 1], &vec![0.0; n - 1], &[], n - 1, 0);
            assert!(e.is_err());
        }
    }

    fn gcd(a: usize, b: usize) -> usize {
        if b == 0 { a } else { gcd(b, a % b) }
    }
}

// The existing scf::pool_fixture is private; keep this test-only import local
// rather than change that unrelated module's visibility outside the ticket.
#[cfg(test)]
#[allow(dead_code, clippy::duplicate_mod)]
#[path = "../../tests/common/mod.rs"]
mod owned_input_fixtures;

#[cfg(test)]
mod owned_input_tests {
    use super::*;
    use dft_basis::{FftGridSpec, HalfSphere};
    use dft_geometry::Cell;

    #[test]
    fn owned_pack_matches_full_validation_and_bits() {
        for lattice in [
            [[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]],
            [[5., 0.3, 0.], [0., 7., 0.4], [0.2, 0., 9.]],
            [[4., 0., 0.], [0., 4., 0.], [0., 0., 16.]],
        ] {
            let cell = Cell::from_user_rows_bohr(lattice).unwrap();
            for cutoff_hartree in [1.2, 2.5] {
                let basis = ScfOptions {
                    cutoff_hartree,
                    ..ScfOptions::strict_reference()
                }
                .basis(&cell, [0.; 3], 1.)
                .unwrap();
                let grid =
                    FftGridSpec::new(&basis, &basis.density_basis(&cell).unwrap(), None).unwrap();
                let half = HalfSphere::new(&basis, &grid).unwrap();
                let n = half.full_len();
                let real: Vec<_> = (0..n)
                    .map(|i| match i % 5 {
                        0 => -0.0,
                        1 => 1e-310,
                        _ => (i as f64).sin() * 0.01,
                    })
                    .collect();
                let mut full = vec![Complex64::default(); n];
                half.unpack(&real, &mut full).unwrap();
                for defect in 0..7 {
                    let mut x = full.clone();
                    match defect {
                        1 => x[0].im = dft_basis::ORC_TOLERANCE,
                        2 => x[0].im = 2. * dft_basis::ORC_TOLERANCE,
                        3 => x[n - 1].re = f64::NAN,
                        4 => x[n - 1].im = f64::INFINITY,
                        5 => x[half.minus()[1] as usize].re += 1e-10,
                        6 => {
                            x[half.plus()[1] as usize] = Complex64::new(f64::MAX, 0.);
                            x[half.minus()[1] as usize] = Complex64::new(f64::MAX, -0.);
                        }
                        _ => {}
                    }
                    let mut reference = vec![42.; n];
                    let verdict = half.pack(&x, &mut reference);
                    for ranks in [1, 2, 3, 8] {
                        for rank in 0..ranks {
                            let owners = dft_solvers::distributed_fft::position_owners(
                                grid.dims(),
                                grid.g_to_fft(),
                                ranks,
                            )
                            .unwrap();
                            let rows: Vec<_> = (0..n)
                                .filter(|&r| {
                                    let g = if r == 0 {
                                        0
                                    } else {
                                        half.plus()[r.div_ceil(2)] as usize
                                    };
                                    owners[g] == rank
                                })
                                .collect();
                            let mut got = vec![42.; rows.len()];
                            assert_eq!(pack_orbital_rows(&half, &x, &rows, &mut got), verdict);
                            for (&r, v) in rows.iter().zip(got) {
                                assert_eq!(v.to_bits(), reference[r].to_bits());
                            }
                        }
                    }
                    // An empty owner still checks every coefficient.
                    assert_eq!(pack_orbital_rows(&half, &x, &[], &mut []), verdict);
                }
                assert_eq!(
                    pack_orbital_rows(&half, &full, &[n], &mut [0.]),
                    Err(dft_basis::HalfSphereError::Dimension)
                );
            }
        }
    }

    #[test]
    fn owned_input_descriptor_invalidation() {
        use dft_physics::data::Array3;
        use dft_physics::hamiltonian::{NonlocalOperator, Species};
        let pool = WorkerPool::new(1).unwrap();
        let mut state = None;
        let mut old_half = None;
        // Each immutable model/layout change supplies a new Hamiltonian. No
        // content cache may accept it merely because its dimensions match.
        for change in 0..6 {
            let edge = if change == 1 { 6.3 } else { 6. };
            let cell =
                Cell::from_user_rows_bohr([[edge, 0., 0.], [0., 6., 0.], [0., 0., 6.]]).unwrap();
            let options = ScfOptions {
                cutoff_hartree: if change == 4 { 2. } else { 1.2 },
                ..ScfOptions::strict_reference()
            };
            let basis = options.basis(&cell, [0.; 3], 1.).unwrap();
            let density = basis.density_basis(&cell).unwrap();
            let initial_grid = FftGridSpec::new(&basis, &density, None).unwrap();
            let override_dims = (change == 5).then(|| initial_grid.dims().map(|n| n * 2));
            let grid = FftGridSpec::new(&basis, &density, override_dims).unwrap();
            let pp = super::owned_input_fixtures::potential(
                if change == 3 { 2. } else { 4. },
                1.,
                true,
                101,
            );
            let positions = [[if change == 2 { 0.4 } else { 0.3 }, 0.7, 1.2]];
            let nl = NonlocalOperator::new(
                &basis,
                cell.volume_bohr3(),
                &[Species {
                    potential: &pp,
                    positions_bohr: &positions,
                }],
            )
            .unwrap();
            let mut h =
                Hamiltonian::new_with_grid(&basis, grid.clone(), Array3::zeros(grid.dims()), nl)
                    .unwrap();
            ensure_state_with(&h, &basis, 3, 6, &options, &pool, &mut state, false).unwrap();
            let half = state.as_ref().unwrap().half.clone();
            if let Some(old) = &old_half {
                assert!(!Arc::ptr_eq(old, &half));
            }
            assert!(state.as_ref().unwrap().op.is_none());
            // Density-dependent potential must NOT invalidate geometry.
            h.set_local_potential(Array3::from_elem(grid.dims(), 0.125))
                .unwrap();
            ensure_state_with(&h, &basis, 3, 6, &options, &pool, &mut state, false).unwrap();
            assert!(Arc::ptr_eq(&half, &state.as_ref().unwrap().half));
            // New rank layouts are constructed, never cached by this bridge.
            for ranks in [1, 2, 3, 8] {
                let mut seen = vec![false; half.full_len()];
                for rank in 0..ranks {
                    let d =
                        DistributedRealHamiltonian::new(&h, &basis, 2, rank, ranks, pool.clone())
                            .unwrap();
                    assert!(d.same_operator(&h, &basis));
                    for &row in d.local_rows() {
                        assert!(!seen[row]);
                        seen[row] = true;
                    }
                }
                assert!(seen.into_iter().all(|x| x));
            }
            old_half = Some(half);
        }
    }
}
