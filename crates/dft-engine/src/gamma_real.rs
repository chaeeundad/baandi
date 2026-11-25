//! Default-off SCF bridge to ORC. Derived from real-gamma-design sections 4–10.
//! Completion is the only output boundary; downstream physics stays complex.
use crate::{ScfError, ScfOptions};
use dft_basis::PlaneWaveBasis;
use dft_mpi::Collective;
use dft_numerics::{RealMatrix, linalg::ComplexMatrix, pool::WorkerPool};
use dft_physics::{
    data::Complex64, distributed_real_gamma::DistributedRealHamiltonian, hamiltonian::Hamiltonian,
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
    /// Serial real Gamma operator. ISSUE-20 (`BAANDI_M2_LAZY_SERIAL=on`): None
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
    /// REUSE-PSI state (work order C2), created at the first eligible solve.
    reuse: Option<PsiReuse>,
    /// MPI-KEPT-START state (`BAANDI_M2_REUSE=on`) of the M2 stick-row path.
    pub(crate) m2_reuse: Option<M2Reuse>,
    /// MPI-LAZY-ORBITALS: half sphere shared by this state's lazy blocks, the
    /// operator's local rows and their ORC kinetic diagonal (Ha) of the latest
    /// lazy solve, and the block that solve returned (taken by the caller).
    lazy_half: Option<Arc<dft_basis::HalfSphere>>,
    lazy_rows: Arc<[usize]>,
    lazy_diag: Vec<f64>,
    pub(crate) lazy_latest: Option<LazyOrbitals>,
    /// Set by an exact-continuation checkpoint: the next solve forgets the kept
    /// X, H X and psi(r) and takes the plain start, exactly as a restored session
    /// (which has no kept state) does, so a continued and a restored run agree bit
    /// for bit.
    pub(crate) forget_reuse: std::sync::atomic::AtomicBool,
}
impl RealSolveState {
    /// NL-ENERGY-REAL values of the latest solve of this state, if computed.
    pub(crate) fn nonlocal_partials(&self) -> Option<&[f64]> {
        self.nl_valid.then_some(self.nl_values.as_slice())
    }
    /// REUSE-PSI: slots for the psi(r) grids of the leading `pairs` band pairs
    /// of the latest solve's orbitals, for the density build to fill. `None`
    /// when the latest solve kept no X and H X, or when the grids would exceed
    /// `PSI_KEEP_MAX_BYTES` (then the next start transforms every pair). The
    /// grids are allocated here once per state (outside numerical loops).
    pub(crate) fn psi_keep(
        &mut self,
        pairs: usize,
        dims: [usize; 3],
    ) -> Option<(&mut [Vec<Complex64>], &mut [bool])> {
        let reuse = self.reuse.as_mut().filter(|r| r.collect)?;
        let grid_len = dims.iter().product::<usize>();
        if dims != reuse.dims || grid_len != reuse.v_old.len() {
            return None;
        }
        if reuse.grids.len() != pairs {
            let bytes = pairs
                .saturating_mul(grid_len)
                .saturating_mul(std::mem::size_of::<Complex64>());
            if bytes > PSI_KEEP_MAX_BYTES {
                if !reuse.gated {
                    reuse.gated = true;
                    reuse_add(ReuseCount::Gated, 1);
                }
                reuse.grids = Vec::new();
                reuse.kept = Vec::new();
                return None;
            }
            reuse.grids = (0..pairs)
                .map(|_| vec![Complex64::default(); grid_len])
                .collect();
            reuse.kept = vec![false; pairs];
            reuse_max(ReuseCount::BytesKept, bytes);
        }
        Some((&mut reuse.grids, &mut reuse.kept))
    }
}

/// REUSE-PSI (`BAANDI_REUSE_PSI=on`, default off; work order beat-QE C2;
/// rounding class). Every SCF map transforms the returned orbitals X to real
/// space for the density, and the next lean Davidson solve starts by applying
/// H to the same X. With the switch on, the lean return keeps X and H X (its
/// return rotation of the stored H images), the density build keeps psi(r) of
/// each occupied band pair, and the next solve starts from
/// H_new X = H_old X + (V_new - V_old) X: the potential change multiplies the
/// kept psi(r), one forward transform per pair; kinetic and nonlocal parts are
/// unchanged and not applied. Pairs without a kept grid (empty bands, odd
/// tails) take the full pair transform with the potential change. Fallback to
/// the plain start whenever the solver's `initial` is not the kept X, the
/// operator was rebuilt (new atoms or projectors, a different basis or grid,
/// the other ladder grid), the solve is not a threads-only lean Ritz-value
/// solve, or the grids exceed the memory gate.
pub(crate) fn reuse_psi() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| dft_numerics::switches::on(dft_numerics::switches::Id::REUSE_PSI))
}
/// REUSE-PSI memory gate for the kept psi(r) grids of one k state (pair
/// complex grids, 16 bytes per point and pair). A resource bound, not a
/// physical parameter: above it the kept H X is still used and every pair is
/// transformed in both directions.
pub const PSI_KEEP_MAX_BYTES: usize = 2 << 30;

pub(crate) struct PsiReuse {
    /// Whether the latest solve kept X and H X (the density build may fill
    /// the grids for them).
    collect: bool,
    /// Grid of the local potential (the density grid must match).
    dims: [usize; 3],
    /// Potential (Ha, grid C order) the kept H X belongs to.
    v_old: Vec<f64>,
    /// Potential change scratch.
    delta: Vec<f64>,
    /// (V_new - V_old) X, n x bands.
    image: RealMatrix,
    /// psi(r) of band pairs (0,1), (2,3), ... of the kept X (unscaled inverse
    /// transforms), valid where `kept` is set.
    grids: Vec<Vec<Complex64>>,
    kept: Vec<bool>,
    /// The memory gate declined the grids for this state.
    gated: bool,
}

#[derive(Clone, Copy)]
enum ReuseCount {
    /// Solves started from the kept X and H X (armed; the solver may still
    /// decline when `initial` is not the kept X).
    Armed = 0,
    /// Armed starts whose potential did not change (no transform at all).
    Unchanged = 1,
    /// Pairs whose potential change came from a kept grid (one forward
    /// transform; the inverse transform is saved).
    KeptPairs = 2,
    /// Pairs transformed in both directions for the potential change.
    FullPairs = 3,
    /// Columns whose nonlocal products (P^T X and the back product) and kinetic
    /// pass were not formed.
    ColumnsSaved = 4,
    /// Eligible solves with no kept start (first solve, rebuilt operator,
    /// changed grid, failed previous solve).
    NotKept = 5,
    /// States whose grids exceeded the memory gate.
    Gated = 6,
    /// Largest kept grid storage of a state, bytes.
    BytesKept = 7,
}
static REUSE_COUNTS: [std::sync::atomic::AtomicU64; 8] =
    [const { std::sync::atomic::AtomicU64::new(0) }; 8];
fn reuse_add(c: ReuseCount, v: usize) {
    REUSE_COUNTS[c as usize].fetch_add(v as u64, std::sync::atomic::Ordering::Relaxed);
}
fn reuse_max(c: ReuseCount, v: usize) {
    REUSE_COUNTS[c as usize].fetch_max(v as u64, std::sync::atomic::Ordering::Relaxed);
}
/// REUSE-PSI counters of this process (all zero with the switch off).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ReusePsiCounters {
    pub armed: u64,
    pub unchanged_potential: u64,
    pub kept_pairs: u64,
    pub full_pairs: u64,
    pub columns_saved: u64,
    pub not_kept: u64,
    pub gated: u64,
    pub bytes_kept: u64,
}
impl ReusePsiCounters {
    /// Inverse transforms saved against the plain start (one per kept pair).
    pub fn transforms_saved(&self) -> u64 {
        self.kept_pairs
    }
}
pub fn reuse_psi_counters() -> ReusePsiCounters {
    let g = |c: ReuseCount| REUSE_COUNTS[c as usize].load(std::sync::atomic::Ordering::Relaxed);
    ReusePsiCounters {
        armed: g(ReuseCount::Armed),
        unchanged_potential: g(ReuseCount::Unchanged),
        kept_pairs: g(ReuseCount::KeptPairs),
        full_pairs: g(ReuseCount::FullPairs),
        columns_saved: g(ReuseCount::ColumnsSaved),
        not_kept: g(ReuseCount::NotKept),
        gated: g(ReuseCount::Gated),
        bytes_kept: g(ReuseCount::BytesKept),
    }
}
/// MPI-KEPT-START (`BAANDI_M2_REUSE=on`, default off; rounding class;
/// docs/design/mpi-kept-start.md): the REUSE-PSI kept start on the M2
/// stick-row path. A standard (Ritz-change) M2 solve keeps the local rows of X
/// and of H X (the return rotation of the H images); the M2 density build keeps
/// the real-space plane slabs of the leading fully occupied pairs; the next
/// standard solve of the same operator starts from X and
/// H_new X = H_old X + (V_new - V_old) X (one forward transform per kept pair,
/// both transforms for the other pairs; no kinetic or nonlocal work, no
/// orthonormalization, no initial apply) when the solver finds `initial` equal
/// to X to 1e-14 and the Gram defect of X at most 1e-12. The local rows of X
/// come from the solver, not from the gathered block. Per-thread override for
/// simulated ranks: `set_m2_reuse_override`.
pub(crate) fn m2_reuse() -> bool {
    M2_REUSE_OVERRIDE.with(|o| o.get()).unwrap_or_else(|| {
        static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
        *ON.get_or_init(|| dft_numerics::switches::on(dft_numerics::switches::Id::M2_REUSE))
    })
}
thread_local! {
    static M2_REUSE_OVERRIDE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}
/// Override `BAANDI_M2_REUSE` for the calling thread (simulated ranks in
/// tests); `None` restores the environment.
pub fn set_m2_reuse_override(value: Option<bool>) {
    M2_REUSE_OVERRIDE.with(|o| o.set(value));
}

/// MPI-LAZY-ORBITALS (`BAANDI_M2_LAZY=on`, default off;
/// docs/design/mpi-lazy-orbitals.md): standard M2 solves of the SCF session
/// return only this rank's stick rows; the full complex block is gathered on
/// demand (`LazyOrbitals::materialize`). Per-thread override for simulated
/// ranks: `set_m2_lazy_override`.
pub(crate) fn m2_lazy() -> bool {
    M2_LAZY_OVERRIDE.with(|o| o.get()).unwrap_or_else(|| {
        static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
        *ON.get_or_init(|| dft_numerics::switches::on(dft_numerics::switches::Id::M2_LAZY))
    })
}
thread_local! {
    static M2_LAZY_OVERRIDE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}
/// Override `BAANDI_M2_LAZY` for the calling thread (simulated ranks in tests);
/// `None` restores the environment.
pub fn set_m2_lazy_override(value: Option<bool>) {
    M2_LAZY_OVERRIDE.with(|o| o.set(value));
}

/// BREAKTHROUGH-4: block capacity of the serial real Hamiltonian on the M2
/// path with `BAANDI_M2_LIGHT_STATE=on` (one pair plan; it is never applied).
const LIGHT_CAPACITY: usize = 2;
/// `BAANDI_M2_LIGHT_STATE=on|off` (default off; per-thread override
/// `set_m2_light_state_override`): build the replicated serial real
/// Hamiltonian of an M2 solve state with the minimal block capacity.
/// Bit-identical (that operator only supplies the ORC half sphere there).
pub(crate) fn m2_light_state() -> bool {
    M2_LIGHT_OVERRIDE.with(|o| o.get()).unwrap_or_else(|| {
        static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
        *ON.get_or_init(|| dft_numerics::switches::on(dft_numerics::switches::Id::M2_LIGHT_STATE))
    })
}
thread_local! {
    static M2_LIGHT_OVERRIDE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}
/// Override `BAANDI_M2_LIGHT_STATE` for the calling thread (simulated ranks
/// in tests); `None` restores the environment.
pub fn set_m2_light_state_override(value: Option<bool>) {
    M2_LIGHT_OVERRIDE.with(|o| o.set(value));
}

/// ISSUE-20 `BAANDI_M2_LAZY_SERIAL=on|off` (default off; per-thread override
/// `set_m2_lazy_serial_override`): an M2 solve state owns only the ORC
/// half-sphere descriptor and the distributed operator; the serial real
/// Hamiltonian is built only if a fallback applies it. Bit-identical.
pub(crate) fn m2_lazy_serial() -> bool {
    M2_LAZY_SERIAL_OVERRIDE
        .with(|o| o.get())
        .unwrap_or_else(|| {
            static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
            *ON.get_or_init(|| {
                dft_numerics::switches::on(dft_numerics::switches::Id::M2_LAZY_SERIAL)
            })
        })
}
thread_local! {
    static M2_LAZY_SERIAL_OVERRIDE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
    static M2_DIST_DISABLED: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}
/// Override `BAANDI_M2_LAZY_SERIAL` for the calling thread (simulated ranks in
/// tests); `None` restores the environment.
pub fn set_m2_lazy_serial_override(value: Option<bool>) {
    M2_LAZY_SERIAL_OVERRIDE.with(|o| o.set(value));
}
/// ISSUE-20 test hook: treat the distributed operator of the calling thread
/// as unavailable (its construction "fails"), so that the M1 fallback, and
/// with it the lazy serial construction, is exercised. Never set in
/// production code.
#[doc(hidden)]
pub fn set_m2_dist_disabled_for_tests(value: bool) {
    M2_DIST_DISABLED.with(|o| o.set(value));
}

use crate::refscf::lazy::LazyInner;
#[cfg(test)]
use crate::refscf::lazy::assemble_lazy_block;
pub(crate) use crate::refscf::lazy::lazy_result_allowed;
pub use crate::refscf::lazy::{LazyOrbitals, set_delivered_orbitals};

/// MPI-LAZY-ORBITALS: the audit's kinetic energy of a lazy block, with the ORC
/// kinetic diagonal of the operator rows of `state` (kept by the solve that
/// produced the block). One tree sum; collective.
pub(crate) fn lazy_kinetic_energy(
    state: &RealSolveState,
    lazy: &LazyOrbitals,
    occupations: &[f64],
) -> Result<f64, ScfError> {
    let start = dft_numerics::stage_clock::start();
    if !lazy.matches(&state.lazy_rows, lazy.0.n) {
        return Err(ScfError::Numerical(
            "lazy kinetic energy: rows differ from the operator".into(),
        ));
    }
    let value = lazy.kinetic_energy(&state.lazy_diag, occupations)?;
    crate::m2::add(crate::m2::Stage::LazyKinetic, start);
    Ok(value)
}

/// MPI-KEPT-START state of one k state on one rank. Every flag is replicated
/// (set under conditions every rank evaluates alike); the matrices hold this
/// rank's rows and the slabs this rank's FFT planes.
pub(crate) struct M2Reuse {
    /// H X of the latest keeping return and the solver flags.
    pub(crate) kept: dft_solvers::distributed_davidson::DistributedKeptStart,
    /// Potential (Ha, full grid) the stored H X belongs to.
    v_old: Vec<f64>,
    /// M2-OWNED: `v_old` holds only this rank's operator planes (replicated
    /// flag: set by every rank for the same solve).
    v_old_owned: bool,
    delta: Vec<f64>,
    /// (V_new - V_old) X on the local rows.
    image: RealMatrix,
    /// Plane slabs of the leading pairs of the kept X (pair t at
    /// `planes[t * plane_points..]`), valid for the first `kept_pairs` pairs.
    pub(crate) planes: Vec<Complex64>,
    pub(crate) kept_pairs: usize,
    /// The latest solve stored a return: the density build may keep slabs.
    pub(crate) collect: bool,
    pub(crate) gated: bool,
}

impl M2Reuse {
    /// Slab storage for `pairs` pairs of `plane_points` points, or None when
    /// nothing was stored or the replicated estimate `pairs * ceil(grid /
    /// ranks)` points exceeds `PSI_KEEP_MAX_BYTES` (the decision must not
    /// depend on the rank's own plane count, so that every rank takes the same
    /// branch). Grows once (setup of the state), never in a numerical loop.
    pub(crate) fn plane_slots(
        &mut self,
        pairs: usize,
        plane_points: usize,
        grid_points: usize,
        ranks: usize,
    ) -> Option<&mut [Complex64]> {
        let estimate = pairs
            .saturating_mul(grid_points.div_ceil(ranks.max(1)))
            .saturating_mul(std::mem::size_of::<Complex64>());
        if !self.collect || pairs == 0 {
            return None;
        }
        if estimate > PSI_KEEP_MAX_BYTES {
            if !self.gated {
                self.gated = true;
                crate::m2::count(crate::m2::Count::ReuseGated);
            }
            return None;
        }
        let len = pairs * plane_points;
        if self.planes.len() < len {
            self.planes.resize(len, Complex64::default());
            crate::m2::count_max(
                crate::m2::Count::ReusePlaneBytes,
                len * std::mem::size_of::<Complex64>(),
            );
        }
        Some(&mut self.planes[..len])
    }
}

/// NL-ENERGY-REAL (opt-in `BAANDI_NL_ENERGY_REAL=on`, default off; rounding
/// class): on the real Gamma path the audit's per-band nonlocal energies come
/// from one real projector product of the returned ORC vectors (P^T X by the
/// nonlocal GEMM) instead of the complex projector sums over the full sphere.
pub(crate) fn nl_energy_real() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| dft_numerics::switches::on(dft_numerics::switches::Id::NL_ENERGY_REAL))
}
/// C3a: audits of the real nonlocal energy that read the solver's kept
/// projections (index 0) or formed P^T X directly (index 1).
static NL_AUDITS: [std::sync::atomic::AtomicU64; 2] =
    [const { std::sync::atomic::AtomicU64::new(0) }; 2];
fn count_nl_audit(cached: bool) {
    NL_AUDITS[usize::from(!cached)].fetch_add(1, std::sync::atomic::Ordering::Relaxed);
}
/// C3a (`BAANDI_LEAN_PROJ_CACHE`): (audits from the kept projections, audits
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
    standard_threshold: Option<u64>,
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
    // B6 hooks (used by the solver only with `BAANDI_B6=on`).
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
        standard_threshold: options.standard.map(|p| p.threshold_ha.to_bits()),
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
        reuse: None,
        m2_reuse: None,
        lazy_half: None,
        lazy_rows: Arc::from(Vec::new()),
        lazy_diag: Vec::new(),
        lazy_latest: None,
        forget_reuse: std::sync::atomic::AtomicBool::new(false),
    });
    Ok(())
}

/// ATOMIC-REAL (opt-in `BAANDI_ATOMIC_REAL=on`, default off; rounding class).
pub(crate) fn atomic_real() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| dft_numerics::switches::on(dft_numerics::switches::Id::ATOMIC_REAL))
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
    lazy: bool,
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
    // M2-OWNED (`BAANDI_M2_LOCAL_FIELDS`): this map's effective potential on
    // the distributed operator's planes and its epoch; the complex
    // Hamiltonian's potential is then stale and must not reach the operator.
    owned: Option<(&[f64], u64)>,
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
    // BREAKTHROUGH-4 (`BAANDI_M2_LIGHT_STATE=on`): on the M2 path the serial
    // real Hamiltonian only supplies its half sphere (ORC packing) and is
    // never applied, so it is built with the minimal block capacity (one pair
    // plan instead of capacity / 2 full-grid pair plans and scratch per rank).
    // Replicated condition; the full operator is rebuilt below if the
    // distributed operator cannot be built (the M1 fallback applies it).
    // ISSUE-20 (`BAANDI_M2_LAZY_SERIAL=on`): on the M2 path the state holds
    // only the half-sphere descriptor; the serial operator is built below
    // only when the distributed operator is unavailable. Replicated condition.
    let lazy_serial = comm.is_some() && crate::parallel::m2_enabled() && m2_lazy_serial();
    let light = !lazy_serial && comm.is_some() && crate::parallel::m2_enabled() && m2_light_state();
    let state_capacity = if light {
        LIGHT_CAPACITY.min(capacity)
    } else {
        capacity
    };
    let ensured = ensure_state_with(
        hamiltonian,
        basis,
        bands,
        state_capacity,
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
    state.lazy_latest = None;
    // MPI-KEPT-START: the previous solve's local rows are the kept X (and the
    // storage of this solve's local rows).
    let previous_local = state.m2_local.take();
    let t_prepare = dft_numerics::stage_clock::start();
    // M2c: the distributed Hamiltonian follows the SCF potential like the serial one.
    match comm.as_ref().filter(|_| crate::parallel::m2_enabled()) {
        Some(c) => {
            // M2-OWNED: identity check only; the owned planes are set below
            // (the complex potential is stale and must not be copied).
            let fresh = match owned {
                Some(_) => state
                    .dist
                    .as_ref()
                    .is_some_and(|d| d.same_operator(hamiltonian, basis)),
                None => state.dist.as_mut().is_some_and(|d| {
                    let fresh = d.update_from_complex(hamiltonian, basis);
                    if fresh {
                        // The full-grid potential copy of the gathered path.
                        use crate::grid_dist as g;
                        g::owned_add(g::O_FULL_COPIES, 1);
                        g::owned_add(g::O_FULL_COPY_VALUES, d.potential().len() as u64);
                    }
                    fresh
                }),
            };
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
                state.m2_reuse = None;
            }
            if let Some(d) = state.dist.as_mut() {
                d.set_batched(crate::m2::batched_transposes());
            }
            match (owned, state.dist.as_mut()) {
                (Some((planes, epoch)), Some(d)) => {
                    if d.owned_epoch().is_none() {
                        // Entering the owned mode fills the full buffer with NaN once.
                        use crate::grid_dist as g;
                        g::owned_add(g::O_FULL_COPIES, 1);
                        g::owned_add(g::O_FULL_COPY_VALUES, d.potential().len() as u64);
                    }
                    d.set_owned_potential(planes, epoch)
                        .map_err(|e| ScfError::Numerical(format!("M2 owned potential: {e}")))?;
                }
                // The serial (M1) fallback would apply the stale complex
                // potential: refuse (replicated: the operator exists on every
                // rank or on none).
                (Some(_), None) => {
                    return Err(ScfError::Numerical(
                        "M2 owned potential: distributed operator unavailable".into(),
                    ));
                }
                _ => {}
            }
        }
        None => {
            state.dist = None;
            state.m2_density = None;
            state.m2_reuse = None;
        }
    }
    if light && state.dist.is_none() {
        // The M1 fallback applies the serial operator: build it in full.
        let mut full = None;
        if let Err(e) = ensure_state(
            hamiltonian,
            basis,
            bands,
            capacity,
            options,
            &pool,
            &mut full,
        ) {
            diagnostics.fallback(option, format!("real Hamiltonian construction: {e}"))?;
            return Ok(None);
        }
        let full = full.expect("built above");
        state.op = full.op;
        state.half = full.half;
        state.key = full.key;
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
    // MPI-LAZY-ORBITALS: standard M2 solves on the workspace path return only
    // the local rows when the caller asks for it (replicated conditions).
    let tuning = dft_solvers::distributed_davidson::DistributedTuning::current();
    let lazy_path = lazy
        && comm.is_some()
        && state.dist.is_some()
        && options.standard.is_some()
        && (tuning.workspace || tuning.selected_roots || tuning.lean_algebra);
    if lazy_path && state.lazy_half.is_none() {
        state.lazy_half = Some(state.half.clone());
    }
    // A lazy warm start (0 x B placeholder plus the rows): this rank's rows
    // enter the solver directly when every rank's rows match its operator
    // (one 1-double tree sum); otherwise the block is gathered and packed.
    let mut lazy_start: Option<RealMatrix> = None;
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
                        ritz: options.standard.is_some(),
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
            if comm.is_some() && state.dist.is_some() {
                lazy_start = Some(inj.x);
            } else {
                injected_packed = Some(inj.x);
            }
            injected_hx = inj.hx;
            None
        }
        (Some(m), Some(h), _) if m.rows() == 0 => {
            let mut local = false;
            if lazy_path && let (Some(c), Some(d)) = (comm.as_ref(), state.dist.as_ref()) {
                let mismatch = !h.matches(d.local_rows(), n);
                let mut flag = [0.0];
                c.tree_sum_f64(&[f64::from(u8::from(mismatch))], &mut flag)
                    .map_err(|e| ScfError::Numerical(format!("lazy start: {e}")))?;
                local = flag[0] == 0.0;
                if !local {
                    crate::m2::count(crate::m2::Count::LazyStartMismatches);
                }
            }
            if local {
                lazy_start = Some(h.round_trip_start()?);
                crate::m2::count(crate::m2::Count::LazyStarts);
                None
            } else {
                Some(h.materialize()?)
            }
        }
        (Some(m), None, _) if m.rows() == 0 => {
            return Err(ScfError::InvalidInput(
                "lazy warm start without its distributed rows",
            ));
        }
        (m, _, _) => m,
    };
    // Only the workspace branch below accepts DistributedStart::Local.
    // Resolve once at solve entry; no cache of geometry/layout is introduced.
    let owned_input_eligible = comm.is_some()
        && state.dist.is_some()
        && options.standard.is_some()
        && (m2_reuse() || lazy_path)
        && (tuning.workspace || tuning.selected_roots || tuning.lean_algebra);
    let owned_input = owned_input_eligible
        && dft_numerics::switches::on(dft_numerics::switches::Id::M2_OWNED_INPUT);
    let input_profile = dft_numerics::switches::InitProfile::current();
    let packed = if let Some(p) = injected_packed.take() {
        Some(p)
    } else if let Some(initial) = initial {
        let rows = state.dist.as_ref().map(|d| d.local_rows());
        let out_rows = if owned_input { rows.unwrap().len() } else { n };
        let mut packed = RealMatrix::zeros(out_rows, initial.cols()).map_err(SolverError::from)?;
        if owned_input_eligible && let Some(p) = input_profile {
            use dft_numerics::switches::RunCounter as C;
            let local_bytes = 8 * rows.unwrap().len() * initial.cols();
            let full_bytes = 8 * n * initial.cols();
            p.count(C::OwnedInputCalls, 1);
            p.count(C::OwnedInputDirectCalls, usize::from(owned_input));
            p.count(C::OwnedInputBaselineBytes, full_bytes + local_bytes);
            p.count(
                C::OwnedInputSavedBytes,
                if owned_input { full_bytes } else { 0 },
            );
            p.count(C::OwnedInputAllocatedBytes, 8 * packed.as_slice().len());
        }
        let converted = if owned_input {
            // Validate every full column, even on ranks owning no rows, before
            // exposing the local start. No full-length ORC temporary.
            (0..initial.cols()).try_for_each(|j| {
                pack_orbital_rows(
                    half,
                    initial.column(j),
                    rows.unwrap(),
                    &mut packed.as_mut_slice()[j * out_rows..(j + 1) * out_rows],
                )
            })
        } else {
            pool.install(|| {
                pair_columns(initial.as_slice(), packed.as_mut_slice(), n, &|src, dst| {
                    half.pack(src, dst)
                })
            })
        };
        if let Err(e) = converted {
            diagnostics.fallback(
                option,
                format!("warm start is not ORC-compatible (1e-14 absolute): {e}"),
            )?;
            return Ok(None);
        }
        if owned_input {
            lazy_start = Some(packed);
            None
        } else {
            Some(packed)
        }
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
    // REUSE-PSI: threads-only Ritz-value (standard) solves of this state.
    let reuse = comm.is_none() && options.standard.is_some() && reuse_psi();
    if state
        .forget_reuse
        .swap(false, std::sync::atomic::Ordering::Relaxed)
    {
        state.reuse = None;
        state.m2_reuse = None;
        if let Some(ws) = state.workspace.as_mut() {
            ws.keep_returned_products(false)?;
        }
    }
    if reuse {
        prepare_kept_start(state, &solver, options, packed.is_some())?;
    } else if state.reuse.take().is_some()
        && let Some(ws) = state.workspace.as_mut()
    {
        ws.keep_returned_products(false)?;
    }
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
    let criterion = if options.standard.is_some() {
        crate::snapshot::Criterion::RitzChange
    } else {
        crate::snapshot::Criterion::ResidualNorm
    };
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
            ritz: options.standard.is_some(),
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
    // The SCF driver supplies the same ratcheted tolerance and band policy to
    // both representations. Standard returns remain explicitly uncertified.
    let m2_kept = comm.is_some()
        && dist.is_some()
        && options.standard.is_some()
        && m2_reuse()
        && (tuning.workspace || tuning.selected_roots || tuning.lean_algebra);
    if !m2_kept {
        state.m2_reuse = None;
    }
    let result = if let (Some(comm), Some(dist)) = (&comm, dist) {
        // M2-OWNED: the operator must hold exactly this map's owned field (or
        // the full complex potential when nothing is owned); anything else is
        // a stale potential and is rejected here, at the consumer.
        let expected = owned.map(|(_, e)| e);
        if dist.owned_epoch() != expected {
            return Err(ScfError::Numerical(format!(
                "M2 owned potential: operator epoch {:?}, map epoch {expected:?}",
                dist.owned_epoch()
            )));
        }
        if expected.is_some() {
            crate::grid_dist::owned_add(crate::grid_dist::O_EPOCH_CHECKS, 1);
        }
        let t_entry = dft_numerics::stage_clock::start();
        crate::parallel::count_distributed_solve();
        crate::parallel::count_distributed_rows_solve();
        let local_diagonal: Vec<f64> = dist.local_rows().iter().map(|&r| diagonal[r]).collect();
        if m2_kept || lazy_path {
            use dft_solvers::distributed_davidson::{DistributedKeptStart, DistributedStart};
            let mut x_local = previous_local
                .unwrap_or_else(|| RealMatrix::new(0, 0, Vec::new()).expect("empty matrix"));
            // Snapshot replay: the injected rows are the kept X; the injected
            // H X is the stored return brought to this potential (no delta).
            if let (Some(hx), Some(l)) = (injected_hx.take(), lazy_start.as_ref()) {
                x_local = RealMatrix::new(l.rows(), l.cols(), l.as_slice().to_vec())
                    .map_err(SolverError::from)?;
                state.m2_reuse = Some(M2Reuse {
                    kept: DistributedKeptStart::with_stored(hx),
                    v_old: dist.potential().to_vec(),
                    v_old_owned: false,
                    delta: Vec::new(),
                    image: RealMatrix::new(0, 0, Vec::new()).expect("empty matrix"),
                    planes: Vec::new(),
                    kept_pairs: 0,
                    collect: true,
                    gated: false,
                });
            }
            let start = match (&lazy_start, &packed) {
                (Some(l), _) => Some(DistributedStart::Local(l)),
                (None, Some(p)) => Some(DistributedStart::Global(p)),
                (None, None) => None,
            };
            let gather = !lazy_path;
            // This rank's rows of the start (for the capture and the diagnostic).
            let start_rows: Option<RealMatrix> = match (&lazy_start, &packed) {
                (Some(_), _) => None,
                (None, Some(p)) => Some(
                    dist.extract_local(p)
                        .map_err(|e| ScfError::Numerical(e.to_string()))?,
                ),
                (None, None) => None,
            };
            if owned_input_eligible
                && let (Some(p), Some(start_rows)) = (input_profile, start_rows.as_ref())
            {
                p.count(
                    dft_numerics::switches::RunCounter::OwnedInputAllocatedBytes,
                    8 * start_rows.as_slice().len(),
                );
            }
            let x_rows: Option<&RealMatrix> = lazy_start.as_ref().or(start_rows.as_ref());
            // BREAKTHROUGH-4: clock of the interface after the solver call.
            let t_finish;
            let result = if m2_kept {
                let r = state.m2_reuse.get_or_insert_with(|| M2Reuse {
                    kept: DistributedKeptStart::new(),
                    v_old: Vec::new(),
                    v_old_owned: false,
                    delta: Vec::new(),
                    image: RealMatrix::new(0, 0, Vec::new()).expect("empty matrix"),
                    planes: Vec::new(),
                    kept_pairs: 0,
                    collect: false,
                    gated: false,
                });
                crate::m2::add(crate::m2::Stage::SolveEntry, t_entry);
                let t_kept = dft_numerics::stage_clock::start();
                prepare_m2_kept_start(
                    comm,
                    dist,
                    r,
                    &x_local,
                    start.is_some(),
                    bands,
                    owned.is_some(),
                )?;
                crate::m2::add(crate::m2::Stage::KeptUpdate, t_kept);
                let armed = r.kept.armed();
                let hx = r.kept.stored().filter(|_| armed);
                // The solver's basis of a taken kept start is the kept X
                // itself (`x_local`), not the round-tripped start: capture
                // and analyse that block when the start is armed.
                let x_rows: Option<&RealMatrix> = if armed { Some(&x_local) } else { x_rows };
                let capture = crate::snapshot::capture_entry(crate::snapshot::Entry {
                    rank: comm.rank(),
                    ranks: comm.size(),
                    layout: crate::snapshot::Layout::StickRows,
                    x: x_rows,
                    hx,
                    armed,
                    rows_fnv: crate::snapshot::fnv_usize(dist.local_rows()),
                    n_global: n,
                    bands,
                    occupied,
                    occupations,
                    criterion,
                    tolerance,
                    strict_bands: policy.strict_bands,
                    buffer_tolerance: policy.buffer_tolerance,
                    ritz: true,
                    potential: dist.potential(),
                    grid: dist.grid_shape(),
                    cutoff_hartree: options.cutoff_hartree,
                    electrons,
                    capacity_per_band,
                });
                let solved = solver.solve_distributed_rows_local(
                    comm,
                    dist,
                    start,
                    &local_diagonal,
                    policy,
                    true,
                    false,
                    options.davidson_capacity_per_band,
                    tuning,
                    &mut x_local,
                    &mut r.kept,
                    gather,
                );
                t_finish = dft_numerics::stage_clock::start();
                if r.kept.taken() {
                    crate::m2::count(crate::m2::Count::ReuseTaken);
                    crate::m2::count_n(crate::m2::Count::ReuseColumnsSaved, bands);
                } else if r.kept.declined() {
                    crate::m2::count(crate::m2::Count::ReuseDeclined);
                }
                let result = solved?;
                if let Some(c) = capture {
                    crate::snapshot::capture_exit(c, &result, r.kept.taken());
                }
                crate::profile::note_solve(
                    result.restarts,
                    result.contracted_restarts,
                    r.kept.taken(),
                );
                // The density build of this map may keep slabs for the stored X.
                r.collect = r.kept.is_stored();
                r.kept_pairs = 0;
                r.v_old.clear();
                if owned.is_some() {
                    // M2-OWNED: the potential of the stored H X on this
                    // rank's planes only (the other entries are not held).
                    dist.read_planes(dist.potential(), &mut r.v_old);
                    r.v_old_owned = true;
                } else {
                    r.v_old.extend_from_slice(dist.potential());
                    r.v_old_owned = false;
                }
                result
            } else {
                // MPI-LAZY-ORBITALS without the kept start: a kept-start object
                // that never keeps (only the local rows are returned in x_local).
                let mut plain = DistributedKeptStart::new();
                let capture = crate::snapshot::capture_entry(crate::snapshot::Entry {
                    rank: comm.rank(),
                    ranks: comm.size(),
                    layout: crate::snapshot::Layout::StickRows,
                    x: x_rows,
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
                    ritz: true,
                    potential: dist.potential(),
                    grid: dist.grid_shape(),
                    cutoff_hartree: options.cutoff_hartree,
                    electrons,
                    capacity_per_band,
                });
                crate::m2::add(crate::m2::Stage::SolveEntry, t_entry);
                let result = solver.solve_distributed_rows_local(
                    comm,
                    dist,
                    start,
                    &local_diagonal,
                    policy,
                    true,
                    false,
                    options.davidson_capacity_per_band,
                    tuning,
                    &mut x_local,
                    &mut plain,
                    gather,
                )?;
                t_finish = dft_numerics::stage_clock::start();
                if let Some(c) = capture {
                    crate::snapshot::capture_exit(c, &result, false);
                }
                crate::profile::note_solve(result.restarts, result.contracted_restarts, false);
                result
            };
            if lazy_path && result.eigenvectors.rows() == 0 {
                // The block stays on the ranks: this rank's rows are its state.
                if *state.lazy_rows != *dist.local_rows() {
                    state.lazy_rows = Arc::from(dist.local_rows());
                }
                state.lazy_diag = local_diagonal;
                let x =
                    RealMatrix::new(x_local.rows(), x_local.cols(), x_local.as_slice().to_vec())
                        .map_err(SolverError::from)?;
                state.lazy_latest = Some(LazyOrbitals(Arc::new(LazyInner {
                    x,
                    rows: state.lazy_rows.clone(),
                    n,
                    half: state.lazy_half.clone().expect("set on the lazy path"),
                    comm: comm.clone(),
                    full: std::sync::OnceLock::new(),
                    handoff_zero_pairs: None,
                })));
                crate::m2::count(crate::m2::Count::LazyReturns);
            }
            // M2d: this rank's rows for the distributed density and energy,
            // straight from the solver (F07: no extraction from the gathered block).
            state.m2_local = Some(x_local);
            crate::m2::add(crate::m2::Stage::SolveFinish, t_finish);
            result
        } else {
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
                ritz: options.standard.is_some(),
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
                options.standard.is_some(),
                certify && options.standard.is_none(),
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
            options.standard.is_some(),
            certify && options.standard.is_none(),
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
            options.standard.is_some(),
            certify && options.standard.is_none(),
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
    if reuse && let (Some(r), Some(ws)) = (state.reuse.as_mut(), state.workspace.as_ref()) {
        // The density build of this map may keep psi(r) for the kept X.
        r.collect = ws.kept_start().is_some();
        r.kept.fill(false);
        r.v_old.copy_from_slice(
            state
                .op
                .as_ref()
                .expect("serial operator on the threads path")
                .h
                .borrow()
                .potential(),
        );
    }
    // NL-ENERGY-REAL: nonlocal expectation values of the returned real vectors
    // for the energy audit of this map (threads-only path).
    // C3a (`BAANDI_LEAN_PROJ_CACHE=on`): a lean return leaves P^T X = (P^T V) C
    // of exactly these vectors in the workspace (same operator, this solve);
    // otherwise the projector product is formed here.
    if comm.is_none()
        && nl_energy_real()
        && let Some(op) = op_opt
    {
        state.nl_values.resize(bands, 0.0);
        let cached = state
            .workspace
            .as_ref()
            .and_then(|ws| ws.returned_projections());
        state.nl_valid = match cached {
            Some(p) => {
                count_nl_audit(true);
                op.h.borrow()
                    .nonlocal_expectation_values_from_projections(p, &mut state.nl_values)
                    .is_ok()
            }
            None => {
                count_nl_audit(false);
                op.h.borrow_mut()
                    .nonlocal_expectation_values(&result.eigenvectors, &mut state.nl_values)
                    .is_ok()
            }
        };
    }
    // MPI-LAZY-ORBITALS: a lazy return has no block to unpack; the caller takes
    // `lazy_latest` and keeps the 0 x B placeholder.
    if result.eigenvectors.rows() == 0 && state.lazy_latest.is_some() {
        diagnostics.gamma_real_used = true;
        diagnostics.fallback_reason = None;
        return Ok(Some(EigenSolution {
            locking: None,
            spectral_completeness: None,
            eigenvalues: result.eigenvalues,
            eigenvectors: ComplexMatrix::zeros(0, bands).map_err(SolverError::from)?,
            residuals: result.residuals,
            converged: result.converged,
            buffer: result.buffer,
            residuals_certified: result.residuals_certified,
            iterations: result.iterations,
            history: result.history,
            operator_applications: result.operator_applications,
        }));
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

/// REUSE-PSI: ask the workspace to keep X and H X at the coming return, and,
/// when the previous solve of this state kept them, bring H X to the current
/// potential, H_new X = H_old X + (V_new - V_old) X (kinetic and nonlocal parts
/// are unchanged because the state was not rebuilt), and arm the kept start.
/// The solver takes it only if its `initial` is that X.
fn prepare_kept_start(
    state: &mut RealSolveState,
    solver: &RealDavidson,
    options: &ScfOptions,
    has_initial: bool,
) -> Result<(), ScfError> {
    let RealSolveState {
        op,
        workspace,
        reuse,
        ..
    } = state;
    let op = op
        .as_ref()
        .ok_or(ScfError::InvalidInput("serial real operator missing"))?;
    let ws = solver.prepare_workspace(op, options.davidson_capacity_per_band, workspace)?;
    ws.keep_returned_products(true)?;
    let mut h = op.h.borrow_mut();
    let dims = h.grid_dims();
    let grid_len = h.potential().len();
    let (n, bands) = (op.n, solver.bands);
    if reuse.as_ref().is_none_or(|r| {
        r.dims != dims
            || r.v_old.len() != grid_len
            || r.image.rows() != n
            || r.image.cols() != bands
    }) {
        *reuse = Some(PsiReuse {
            collect: false,
            dims,
            v_old: vec![0.; grid_len],
            delta: vec![0.; grid_len],
            image: RealMatrix::zeros(n, bands).map_err(SolverError::from)?,
            grids: Vec::new(),
            kept: Vec::new(),
            gated: false,
        });
    }
    let r = reuse.as_mut().unwrap();
    let collected = std::mem::replace(&mut r.collect, false);
    let armed = match ws.kept_start_mut() {
        Some((x, hx)) if collected && has_initial => {
            let mut changed = false;
            for ((d, &new), &old) in r.delta.iter_mut().zip(h.potential()).zip(&r.v_old) {
                *d = new - old;
                changed |= *d != 0.;
            }
            if changed {
                h.apply_local_kept(x, &r.grids, &r.kept, &r.delta, &mut r.image)
                    .map_err(|e| ScfError::Numerical(e.to_string()))?;
                for (y, &d) in hx.as_mut_slice().iter_mut().zip(r.image.as_slice()) {
                    *y += d;
                }
                let pairs = bands.div_ceil(2);
                let kept = r.kept.iter().take(pairs).filter(|&&k| k).count();
                reuse_add(ReuseCount::KeptPairs, kept);
                reuse_add(ReuseCount::FullPairs, pairs - kept);
            } else {
                reuse_add(ReuseCount::Unchanged, 1);
            }
            true
        }
        _ => false,
    };
    if armed {
        ws.arm_kept_start();
        reuse_add(ReuseCount::Armed, 1);
        reuse_add(ReuseCount::ColumnsSaved, bands);
    } else {
        ws.invalidate_kept_start();
        reuse_add(ReuseCount::NotKept, 1);
    }
    Ok(())
}

/// MPI-KEPT-START: keep H X at the coming return, and, when the previous solve
/// of this operator stored X (`x`, local rows) and H X, bring H X to the current
/// potential, H_new X = H_old X + (V_new - V_old) X, and arm the kept start
/// (the solver takes it only if `initial` is that X). Every condition is
/// replicated, so all ranks enter the same collectives.
fn prepare_m2_kept_start(
    comm: &crate::parallel::SharedComm,
    dist: &mut DistributedRealHamiltonian,
    r: &mut M2Reuse,
    x: &RealMatrix,
    has_initial: bool,
    bands: usize,
    owned: bool,
) -> Result<(), ScfError> {
    r.kept.set_keep(true);
    let collected = std::mem::replace(&mut r.collect, false);
    let grid = dist.potential().len();
    let stored = r.kept.is_stored();
    // `v_old_owned` is replicated; an owned v_old has this rank's plane count.
    let v_old_valid = if r.v_old_owned {
        r.v_old.len() == dist.plane_points()
    } else {
        r.v_old.len() == grid
    };
    if !(collected && stored && has_initial && v_old_valid && x.cols() == bands) {
        r.kept.invalidate();
        crate::m2::count(crate::m2::Count::ReuseNotKept);
        return Ok(());
    }
    let changed = if owned {
        // M2-OWNED: V_new - V_old on this rank's planes only (the apply reads
        // no other entry; the rest of `delta` stays 0); the "changed" decision
        // is agreed over the ranks by one 1-double tree sum (the full-grid
        // test of the replicated path, same verdict).
        let mut planes = Vec::new();
        dist.read_planes(dist.potential(), &mut planes);
        let mut old = Vec::new();
        let old: &[f64] = if r.v_old_owned {
            &r.v_old
        } else {
            dist.read_planes(&r.v_old, &mut old);
            &old
        };
        let mut local = false;
        for (p, &o) in planes.iter_mut().zip(old) {
            *p -= o;
            local |= *p != 0.;
        }
        if r.delta.len() != grid {
            r.delta.clear();
            r.delta.resize(grid, 0.);
        }
        dist.write_planes(&planes, &mut r.delta);
        let mut any = [0.0];
        comm.tree_sum_f64(&[f64::from(u8::from(local))], &mut any)
            .map_err(|e| ScfError::Numerical(format!("M2 kept start: {e}")))?;
        any[0] != 0.
    } else {
        r.delta.clear();
        let mut changed = false;
        for (&new, &old) in dist.potential().iter().zip(&r.v_old) {
            let d = new - old;
            changed |= d != 0.;
            r.delta.push(d);
        }
        changed
    };
    if changed {
        if r.image.rows() != x.rows() || r.image.cols() != bands {
            r.image = RealMatrix::zeros(x.rows(), bands).map_err(SolverError::from)?;
        }
        let kept_pairs = r.kept_pairs.min(bands / 2);
        dist.apply_local_change(comm, x, &r.planes, kept_pairs, &r.delta, &mut r.image)
            .map_err(|e| ScfError::Numerical(format!("M2 kept start: {e}")))?;
        let hx = r.kept.stored_mut().expect("stored");
        for (y, &d) in hx.as_mut_slice().iter_mut().zip(r.image.as_slice()) {
            *y += d;
        }
    } else {
        crate::m2::count(crate::m2::Count::ReuseUnchanged);
    }
    r.kept.arm();
    crate::m2::count(crate::m2::Count::ReuseArmed);
    Ok(())
}

// Boundary conversions use the same fixed pair partition as local application.
// They need only disjoint output slices; all storage belongs to the caller.
/// Owned-row form of HalfSphere::pack. Keep its full-column validation and
/// arithmetic exactly: projector packing has a different tolerance contract
/// and must not be substituted here. Rows come from the current operator.
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
mod reuse_psi_tests {
    use super::*;
    use dft_geometry::Cell;
    use dft_physics::data::Array3;
    use dft_physics::hamiltonian::NonlocalOperator;

    /// REUSE-PSI memory gate and validity: grids are offered only after a
    /// solve kept X and H X, only for the local-potential grid, and never
    /// above `PSI_KEEP_MAX_BYTES` (then nothing is allocated and the state is
    /// counted once as gated).
    #[test]
    fn psi_keep_gate_and_validity() {
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
        let h = Hamiltonian::new(
            &cell,
            &basis,
            Array3::zeros(grid.dims()),
            NonlocalOperator::new(&basis, cell.volume_bohr3(), &[]).unwrap(),
        )
        .unwrap();
        let pool = WorkerPool::new(1).unwrap();
        let mut state = None;
        let capacity = real_capacity(&options, 4, h.kinetic_diagonal().len());
        ensure_state(&h, &basis, 4, capacity, &options, &pool, &mut state).unwrap();
        let state = state.as_mut().unwrap();
        let dims = grid.dims();
        let len: usize = dims.iter().product();
        assert!(state.psi_keep(1, dims).is_none(), "no reuse state yet");
        state.reuse = Some(PsiReuse {
            collect: false,
            dims,
            v_old: vec![0.; len],
            delta: vec![0.; len],
            image: RealMatrix::zeros(state.op.as_ref().unwrap().n, 4).unwrap(),
            grids: Vec::new(),
            kept: Vec::new(),
            gated: false,
        });
        assert!(state.psi_keep(1, dims).is_none(), "nothing kept at return");
        state.reuse.as_mut().unwrap().collect = true;
        assert!(state.psi_keep(1, [dims[0] + 1, dims[1], dims[2]]).is_none());
        let (grids, kept) = state.psi_keep(2, dims).expect("two pairs fit");
        assert_eq!((grids.len(), kept.len(), grids[1].len()), (2, 2, len));
        let gated_before = reuse_psi_counters().gated;
        let too_many = PSI_KEEP_MAX_BYTES / (len * std::mem::size_of::<Complex64>()) + 1;
        assert!(state.psi_keep(too_many, dims).is_none());
        let reuse = state.reuse.as_ref().unwrap();
        assert!(reuse.gated && reuse.grids.is_empty());
        assert!(reuse_psi_counters().gated > gated_before);
    }
}

#[cfg(test)]
mod persistence_tests {
    include!("../tests/common/persistent_real_invalidation_unit.rs");
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
#[path = "../tests/common/mod.rs"]
mod owned_input_fixtures;

#[cfg(test)]
mod owned_input_tests {
    use super::*;
    use crate::{ScfProblem, ScfResult, SpeciesPotential};
    use dft_basis::{FftGridSpec, HalfSphere};
    use dft_geometry::{Atom, Cell};
    use dft_mpi::{SimRank, simulate};
    use dft_numerics::switches::{RunConfig, RunCounter as C};

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
                    ..ScfOptions::default()
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
                ..ScfOptions::default()
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

    fn history(r: &ScfResult) -> String {
        let mut history = r.history.clone();
        for h in &mut history {
            h.timings = Default::default();
            if let Some(c) = h.fixed_point_certificate.as_mut() {
                c.timings = Default::default();
            }
        }
        format!("{history:?}")
    }

    fn run_pair(
        cell: &Cell,
        atoms: &[Atom],
        pp: &dft_upf::Pseudopotential,
        ranks: usize,
        options: &ScfOptions,
        fallback: bool,
    ) {
        let run = |direct| {
            simulate(ranks, |comm: &SimRank| {
                let config = RunConfig::from_pairs(&[
                    ("m2_owned_input", if direct { "on" } else { "off" }),
                    ("init_profile", "on"),
                    ("stage_timers", "off"),
                    // Fallback may not consume an owned potential.
                    ("m2_local_fields", if fallback { "off" } else { "on" }),
                    ("m2_local_density", if fallback { "off" } else { "on" }),
                ])
                .unwrap();
                let _scope = config.enter();
                crate::parallel::set_m2_override(Some(true));
                crate::parallel::set_communicator(Some(Arc::new(comm.clone())));
                set_m2_dist_disabled_for_tests(fallback);
                // Synthetic fixtures lack atomic orbitals. Exercise their
                // materialized warm inputs instead of the already-local bypass.
                if pp.atomic_wavefunctions().is_none_or(|w| w.is_empty()) {
                    set_m2_lazy_override(Some(false));
                }
                let problem = ScfProblem::new(
                    cell,
                    atoms,
                    &[SpeciesPotential {
                        label: "Xx",
                        potential: pp,
                    }],
                    options.clone(),
                )
                .unwrap();
                let result = problem.run().unwrap();
                set_m2_dist_disabled_for_tests(false);
                set_m2_lazy_override(None);
                crate::parallel::set_communicator(None);
                crate::parallel::set_m2_override(None);
                (result, config.counters())
            })
        };
        let a = run(false);
        let b = run(true);
        for (rank, ((a, ca), (b, cb))) in a.iter().zip(&b).enumerate() {
            assert_eq!(a.total_hartree().to_bits(), b.total_hartree().to_bits());
            assert_eq!(
                a.density.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                b.density.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
            );
            assert_eq!(history(a), history(b));
            assert_eq!(
                format!("{:?}", a.k_points[0].solution),
                format!("{:?}", b.k_points[0].solution)
            );
            if options.standard.is_some() {
                assert!(b.verification.as_ref().unwrap().passed());
            }
            let g = |c: &[u64], k: C| c[k as usize];
            assert_eq!(g(ca, C::OwnedInputDirectCalls), 0);
            assert_eq!(g(ca, C::OwnedInputSavedBytes), 0);
            assert_eq!(
                g(ca, C::OwnedInputAllocatedBytes),
                g(ca, C::OwnedInputBaselineBytes)
            );
            assert_eq!(
                g(cb, C::OwnedInputAllocatedBytes) + g(cb, C::OwnedInputSavedBytes),
                g(cb, C::OwnedInputBaselineBytes)
            );
            assert_eq!(
                g(ca, C::OwnedInputBaselineBytes),
                g(cb, C::OwnedInputBaselineBytes)
            );
            if ranks > 1 && !fallback && options.standard.is_some() {
                assert!(g(cb, C::OwnedInputCalls) > 0);
                assert_eq!(g(cb, C::OwnedInputCalls), g(cb, C::OwnedInputDirectCalls));
                let basis = options.basis(cell, [0.; 3], 1.).unwrap();
                let grid =
                    FftGridSpec::new(&basis, &basis.density_basis(cell).unwrap(), None).unwrap();
                let half = HalfSphere::new(&basis, &grid).unwrap();
                let owners = dft_solvers::distributed_fft::position_owners(
                    grid.dims(),
                    grid.g_to_fft(),
                    ranks,
                )
                .unwrap();
                let local = (0..half.full_len())
                    .filter(|&r| {
                        let p = if r == 0 {
                            0
                        } else {
                            half.plus()[r.div_ceil(2)] as usize
                        };
                        owners[p] == rank
                    })
                    .count();
                let bands = b.k_points[0].solution.eigenvalues.len();
                let calls = g(cb, C::OwnedInputCalls);
                assert_eq!(
                    g(cb, C::OwnedInputSavedBytes),
                    calls * (8 * half.full_len() * bands) as u64
                );
                assert_eq!(
                    g(cb, C::OwnedInputAllocatedBytes),
                    calls * (8 * local * bands) as u64
                );
            } else {
                assert_eq!(g(cb, C::OwnedInputCalls), 0);
            }
            println!(
                "OWNED ranks={ranks} rank={rank} atoms={} fallback={fallback} energy_bits={:016x} scf={} calls={} baseline_bytes={} allocated_bytes={} saved_bytes={}",
                atoms.len(),
                b.total_hartree().to_bits(),
                b.history.len(),
                g(cb, C::OwnedInputCalls),
                g(cb, C::OwnedInputBaselineBytes),
                g(cb, C::OwnedInputAllocatedBytes),
                g(cb, C::OwnedInputSavedBytes)
            );
        }
    }

    #[test]
    #[ignore = "verified synthetic SCF, fallback and strict controls; run in release"]
    fn owned_input_synthetic_fallback_and_strict() {
        let cell =
            Cell::from_user_rows_bohr([[6.2, 0., 0.], [1.1, 5.7, 0.], [0.4, 0.7, 6.6]]).unwrap();
        let atoms = [
            Atom::new("Xx", [0.13, 0.21, 0.32]),
            Atom::new("Xx", [0.61, 0.57, 0.74]),
        ];
        for nonlocal in [false, true] {
            let pp = super::owned_input_fixtures::potential(4., 1., nonlocal, 101);
            for mode in 0..3 {
                let options = ScfOptions {
                    cutoff_hartree: 3.,
                    bands: Some(6),
                    ..if mode == 2 {
                        ScfOptions::strict_reference()
                    } else {
                        ScfOptions::default()
                    }
                };
                run_pair(&cell, &atoms, &pp, 2, &options, mode == 1);
            }
        }
    }

    #[test]
    #[ignore = "full verified SCF matrix; run release with DFT_SI_UPF (see OWNED audit)"]
    fn owned_input_si_matrix() {
        let pp = dft_upf::parse_upf(
            &std::fs::read_to_string(
                std::env::var("DFT_SI_UPF").expect("external Si.pz-vbc.UPF required"),
            )
            .unwrap(),
        )
        .unwrap();
        let base = [
            [0., 0., 0.],
            [0., 0.5, 0.5],
            [0.5, 0., 0.5],
            [0.5, 0.5, 0.],
            [0.25, 0.25, 0.25],
            [0.25, 0.75, 0.75],
            [0.75, 0.25, 0.75],
            [0.75, 0.75, 0.25],
        ];
        for (rx, ry, skew) in [(1, 1, false), (2, 1, false), (2, 2, false), (1, 1, true)] {
            let a = 10.261_212_856_717_933;
            let mut cell = Cell::from_user_rows_bohr([
                [a * rx as f64, 0., 0.],
                [if skew { 0.7 } else { 0. }, a * ry as f64, 0.],
                [0., if skew { 0.3 } else { 0. }, a],
            ])
            .unwrap();
            let mut atoms = Vec::new();
            for x in 0..rx {
                for y in 0..ry {
                    for p in base {
                        atoms.push(Atom::new(
                            "Xx",
                            [
                                (p[0] + x as f64) / rx as f64,
                                (p[1] + y as f64) / ry as f64,
                                p[2],
                            ],
                        ));
                    }
                }
            }
            let options = ScfOptions {
                cutoff_hartree: 15.,
                start_wavefunctions: crate::StartWavefunctions::Atomic,
                davidson_capacity_per_band: 3,
                ..ScfOptions::default()
            };
            if skew {
                // Existing benchmarks/cases/si-skew.toml, not a newly strained
                // conventional cell with an unqualified convergence trajectory.
                let h = 5.130_606_428_358_966;
                cell = Cell::from_user_rows_bohr([[0.7, h, h], [h, 0., h], [h, h, 0.3]]).unwrap();
                atoms = vec![Atom::new("Xx", [0.; 3]), Atom::new("Xx", [0.25; 3])];
            }
            for ranks in [1, 2, 3, 8] {
                run_pair(&cell, &atoms, &pp, ranks, &options, false);
            }
        }
    }
}
