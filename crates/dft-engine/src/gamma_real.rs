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

/// Per-k dispatch evidence for the latest solve, also available before completion.
#[derive(Clone, Debug, Default)]
pub struct GammaRealDiagnostics {
    pub gamma_real_used: bool,
    pub fallback_reason: Option<String>,
    /// BREAKTHROUGH-3: the session's `Advance` request for the coming map and
    /// its outcome (taken by `step_map` after the solve; never set otherwise).
    pub(crate) advance: Option<crate::advance::AdvanceSlot>,
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
    op: Operator,
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

/// MPI-LAZY-ORBITALS: the eigenvectors of one standard M2 solve held as this
/// rank's stick rows (ORC real coordinates, the operator's `local_rows` order),
/// the authoritative SCF state of the M2 path. The complex N x B block is
/// formed only by `materialize` (three allgathervs and the Hermitian completion
/// of the gathered path, so its bits equal the gathered block's) and cached
/// with the handle. Cloning shares the data. Every rank must call
/// `materialize` at the same point of the (replicated) control flow.
#[derive(Clone)]
pub struct LazyOrbitals(Arc<LazyInner>);
struct LazyInner {
    /// n_loc x bands, this rank's rows.
    x: RealMatrix,
    /// Global ORC row index of each local row.
    rows: Arc<[usize]>,
    /// Global ORC length N.
    n: usize,
    half: Arc<dft_basis::HalfSphere>,
    comm: crate::parallel::SharedComm,
    full: std::sync::OnceLock<ComplexMatrix>,
}
impl std::fmt::Debug for LazyOrbitals {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("LazyOrbitals")
            .field("n", &self.0.n)
            .field("local_rows", &self.0.x.rows())
            .field("bands", &self.0.x.cols())
            .field("materialized", &self.0.full.get().is_some())
            .finish()
    }
}
impl LazyOrbitals {
    pub fn bands(&self) -> usize {
        self.0.x.cols()
    }
    /// Global plane-wave count (rows of the materialized block).
    pub fn global_len(&self) -> usize {
        self.0.n
    }
    /// Whether this handle holds the rows of `rows` (the operator's local
    /// rows) of an N-row block.
    fn matches(&self, rows: &[usize], n: usize) -> bool {
        self.0.n == n && *self.0.rows == *rows
    }
    /// The local rows after the elementwise ORC round trip sqrt2 * (x / sqrt2)
    /// that pack(unpack(x)) applies to every coordinate but the G = 0 row:
    /// bit-identical to the local rows of the packed gathered warm start.
    fn round_trip_start(&self) -> Result<RealMatrix, ScfError> {
        let x = &self.0.x;
        let nl = x.rows();
        let mut out = RealMatrix::zeros(nl, x.cols()).map_err(SolverError::from)?;
        let s = std::f64::consts::SQRT_2;
        for (dst, src) in out
            .as_mut_slice()
            .chunks_mut(nl.max(1))
            .zip(x.as_slice().chunks(nl.max(1)))
        {
            for ((d, &v), &g) in dst.iter_mut().zip(src).zip(self.0.rows.iter()) {
                *d = if g == 0 { v } else { s * (v / s) };
            }
        }
        Ok(out)
    }
    /// The full complex N x B block (gathered once, then cached). Collective:
    /// every rank must call it at the same point.
    pub fn materialize(&self) -> Result<&ComplexMatrix, ScfError> {
        if let Some(m) = self.0.full.get() {
            return Ok(m);
        }
        let start = dft_numerics::stage_clock::start();
        let inner = &*self.0;
        let comm = &inner.comm;
        let size = comm.size();
        let (n, nl, bands) = (inner.n, inner.x.rows(), inner.x.cols());
        let err = |e: dft_mpi::MpiError| ScfError::Numerical(format!("lazy orbitals: {e}"));
        let mut counts_f = vec![0.0; size];
        comm.allgatherv_f64(&[nl as f64], &vec![1; size], &mut counts_f)
            .map_err(err)?;
        let counts: Vec<usize> = counts_f.iter().map(|&c| c as usize).collect();
        if counts.iter().sum::<usize>() != n {
            return Err(ScfError::Numerical(
                "lazy orbitals: row partition does not cover the basis".into(),
            ));
        }
        let mine: Vec<f64> = inner.rows.iter().map(|&r| r as f64).collect();
        let mut index = vec![0.0; n];
        comm.allgatherv_f64(&mine, &counts, &mut index)
            .map_err(err)?;
        let block: Vec<usize> = counts.iter().map(|c| c * bands).collect();
        let mut flat = vec![0.0; n * bands];
        comm.allgatherv_f64(inner.x.as_slice(), &block, &mut flat)
            .map_err(err)?;
        // Scatter to global ORC order, column by column (exact data movement).
        let mut real = vec![0.0; n * bands];
        let (mut at, mut offset) = (0, 0);
        for &c in &counts {
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
            inner
                .half
                .unpack(&real[j * n..(j + 1) * n], full.column_mut(j))
                .map_err(|e| ScfError::Numerical(format!("lazy orbitals: {e}")))?;
        }
        let _ = inner.full.set(full);
        crate::m2::count(crate::m2::Count::LazyMaterializations);
        crate::m2::add(crate::m2::Stage::LazyMaterialize, start);
        Ok(inner.full.get().expect("set above"))
    }
    /// The full block, owned (gathered if not yet; moved out of the cache when
    /// this is the last handle, cloned otherwise). Collective like `materialize`.
    pub fn into_full(self) -> Result<ComplexMatrix, ScfError> {
        self.materialize()?;
        match Arc::try_unwrap(self.0) {
            Ok(inner) => Ok(inner.full.into_inner().expect("materialized")),
            Err(shared) => Ok(shared.full.get().expect("materialized").clone()),
        }
    }
    /// Occupation-weighted kinetic energy sum_b f_b sum_G |G|^2/2 |c_b(G)|^2 from
    /// the local rows: sum_r d_r x_rb^2 with the ORC kinetic diagonal d (the
    /// representative's |G|^2/2 on both of its coordinates, |G|^2/2 at G = 0),
    /// band values reduced with one tree sum, then summed in band order.
    /// Rounding class against the complex sum (another summation order).
    fn kinetic_energy(&self, diagonal: &[f64], occupations: &[f64]) -> Result<f64, ScfError> {
        let x = &self.0.x;
        let nl = x.rows();
        let m = occupations
            .iter()
            .rposition(|&f| f != 0.0)
            .map_or(0, |i| i + 1);
        if diagonal.len() != nl || occupations.len() != x.cols() {
            return Err(ScfError::InvalidInput("lazy kinetic energy dimensions"));
        }
        let mut partial = vec![0.0; m];
        for (p, column) in partial.iter_mut().zip(x.as_slice().chunks(nl.max(1))) {
            *p = column
                .iter()
                .zip(diagonal)
                .map(|(&v, &d)| d * v * v)
                .sum::<f64>();
        }
        let mut values = vec![0.0; m];
        self.0
            .comm
            .tree_sum_f64(&partial, &mut values)
            .map_err(|e| ScfError::Numerical(format!("lazy kinetic energy: {e}")))?;
        let mut sum = 0.0;
        for (v, f) in values.iter().zip(occupations) {
            if *f != 0.0 {
                sum += f * v;
            }
        }
        Ok(sum)
    }
}

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
            && Arc::ptr_eq(&state.op.pool, pool)
            && state.op.h.get_mut().update_from_complex(hamiltonian, basis)
    });
    if !reusable {
        let real = RealHamiltonian::from_complex(hamiltonian, basis, capacity, pool.clone())?;
        *state = Some(RealSolveState {
            op: Operator {
                n: real.half_sphere().full_len(),
                h: RefCell::new(real),
                pool: pool.clone(),
            },
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
    }
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
    let mut h = state.op.h.borrow_mut();
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
    if let Err(e) = ensure_state(hamiltonian, basis, bands, capacity, options, &pool, state) {
        *state = None;
        diagnostics.fallback(option, format!("real Hamiltonian construction: {e}"))?;
        return Ok(None);
    }
    let state = state.as_mut().unwrap();
    state.nl_valid = false;
    state.lazy_latest = None;
    // MPI-KEPT-START: the previous solve's local rows are the kept X (and the
    // storage of this solve's local rows).
    let previous_local = state.m2_local.take();
    let t_prepare = dft_numerics::stage_clock::start();
    // M2c: the distributed Hamiltonian follows the SCF potential like the serial one.
    match comm.as_ref().filter(|_| crate::parallel::m2_enabled()) {
        Some(c) => {
            let fresh = state
                .dist
                .as_mut()
                .is_some_and(|d| d.update_from_complex(hamiltonian, basis));
            if !fresh {
                state.dist = DistributedRealHamiltonian::new(
                    hamiltonian,
                    basis,
                    capacity,
                    c.rank(),
                    c.size(),
                    pool.clone(),
                )
                .ok();
                state.m2_density = None;
                // A rebuilt operator (atoms, projectors, basis, grid): nothing kept.
                state.m2_reuse = None;
            }
            if let Some(d) = state.dist.as_mut() {
                d.set_batched(crate::m2::batched_transposes());
            }
        }
        None => {
            state.dist = None;
            state.m2_density = None;
            state.m2_reuse = None;
        }
    }
    let real = state.op.h.borrow();
    let half = real.half_sphere();
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
        state.lazy_half = Some(Arc::new(half.clone()));
    }
    // A lazy warm start (0 x B placeholder plus the rows): this rank's rows
    // enter the solver directly when every rank's rows match its operator
    // (one 1-double tree sum); otherwise the block is gathered and packed.
    let mut lazy_start: Option<RealMatrix> = None;
    let initial = match (initial, initial_lazy) {
        (Some(m), Some(h)) if m.rows() == 0 => {
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
        (Some(m), None) if m.rows() == 0 => {
            return Err(ScfError::InvalidInput(
                "lazy warm start without its distributed rows",
            ));
        }
        (m, _) => m,
    };
    let packed = if let Some(initial) = initial {
        let mut packed = RealMatrix::zeros(n, initial.cols()).map_err(SolverError::from)?;
        if let Err(e) = pool.install(|| {
            pair_columns(initial.as_slice(), packed.as_mut_slice(), n, &|src, dst| {
                half.pack(src, dst)
            })
        }) {
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
    let op = &state.op;
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
    if comm.is_some()
        && let Some(slot) = diagnostics.advance.as_mut()
    {
        // BREAKTHROUGH-3: the distributed solvers have no Advance contract yet.
        slot.declined = Some("distributed solver");
    }
    let result = if let (Some(comm), Some(dist)) = (&comm, dist) {
        crate::parallel::count_distributed_solve();
        crate::parallel::count_distributed_rows_solve();
        let local_diagonal: Vec<f64> = dist.local_rows().iter().map(|&r| diagonal[r]).collect();
        if m2_kept || lazy_path {
            use dft_solvers::distributed_davidson::{DistributedKeptStart, DistributedStart};
            let mut x_local = previous_local
                .unwrap_or_else(|| RealMatrix::new(0, 0, Vec::new()).expect("empty matrix"));
            let start = match (&lazy_start, &packed) {
                (Some(l), _) => Some(DistributedStart::Local(l)),
                (None, Some(p)) => Some(DistributedStart::Global(p)),
                (None, None) => None,
            };
            let gather = !lazy_path;
            let result = if m2_kept {
                let r = state.m2_reuse.get_or_insert_with(|| M2Reuse {
                    kept: DistributedKeptStart::new(),
                    v_old: Vec::new(),
                    delta: Vec::new(),
                    image: RealMatrix::new(0, 0, Vec::new()).expect("empty matrix"),
                    planes: Vec::new(),
                    kept_pairs: 0,
                    collect: false,
                    gated: false,
                });
                prepare_m2_kept_start(comm, dist, r, &x_local, start.is_some(), bands)?;
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
                if r.kept.taken() {
                    crate::m2::count(crate::m2::Count::ReuseTaken);
                    crate::m2::count_n(crate::m2::Count::ReuseColumnsSaved, bands);
                } else if r.kept.declined() {
                    crate::m2::count(crate::m2::Count::ReuseDeclined);
                }
                let result = solved?;
                // The density build of this map may keep slabs for the stored X.
                r.collect = r.kept.is_stored();
                r.kept_pairs = 0;
                r.v_old.clear();
                r.v_old.extend_from_slice(dist.potential());
                result
            } else {
                // MPI-LAZY-ORBITALS without the kept start: a kept-start object
                // that never keeps (only the local rows are returned in x_local).
                let mut plain = DistributedKeptStart::new();
                solver.solve_distributed_rows_local(
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
                )?
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
                })));
                crate::m2::count(crate::m2::Count::LazyReturns);
            }
            // M2d: this rank's rows for the distributed density and energy,
            // straight from the solver (F07: no extraction from the gathered block).
            state.m2_local = Some(x_local);
            result
        } else {
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
            // M2d: keep this rank's rows for the distributed density and energy.
            state.m2_local = Some(
                dist.extract_local(&result.eigenvectors)
                    .map_err(|e| ScfError::Numerical(e.to_string()))?,
            );
            result
        }
    } else if let Some(comm) = &comm {
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
    } else if let Some(slot) = diagnostics
        .advance
        .as_mut()
        .filter(|_| options.standard.is_some())
    {
        // BREAKTHROUGH-3: bounded subspace update of this map. A solver
        // configuration that cannot advance (classic Davidson, second
        // vector) or a solver error declines with its reason and the ordinary
        // converging solve runs instead.
        match solver.advance_reusing(
            op,
            packed.as_ref(),
            &diagonal,
            policy,
            slot.budget,
            options.davidson_capacity_per_band,
            &mut state.workspace,
            options.preconditioner,
        ) {
            Ok(advanced) => {
                slot.outcome = Some(advanced.diagnostics);
                advanced.solution
            }
            Err(e) => {
                slot.declined = Some(match e {
                    SolverError::InvalidParameter(_) => "lean Ritz-change solver not active",
                    SolverError::NotConverged { .. } => "advance exhausted its iteration limit",
                    _ => "advance solver error",
                });
                solver.solve_reusing_with_preconditioner(
                    op,
                    packed.as_ref(),
                    &diagonal,
                    policy,
                    true,
                    false,
                    options.davidson_capacity_per_band,
                    &mut state.workspace,
                    options.preconditioner,
                )?
            }
        }
    } else {
        if let Some(slot) = diagnostics.advance.as_mut() {
            slot.declined = Some("strict solve");
        }
        solver.solve_reusing_with_preconditioner(
            op,
            packed.as_ref(),
            &diagonal,
            policy,
            options.standard.is_some(),
            certify && options.standard.is_none(),
            options.davidson_capacity_per_band,
            &mut state.workspace,
            options.preconditioner,
        )?
    };
    if reuse && let (Some(r), Some(ws)) = (state.reuse.as_mut(), state.workspace.as_ref()) {
        // The density build of this map may keep psi(r) for the kept X.
        r.collect = ws.kept_start().is_some();
        r.kept.fill(false);
        r.v_old.copy_from_slice(state.op.h.borrow().potential());
    }
    // NL-ENERGY-REAL: nonlocal expectation values of the returned real vectors
    // for the energy audit of this map (threads-only path).
    // C3a (`BAANDI_LEAN_PROJ_CACHE=on`): a lean return leaves P^T X = (P^T V) C
    // of exactly these vectors in the workspace (same operator, this solve);
    // otherwise the projector product is formed here.
    if comm.is_none() && nl_energy_real() {
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
        let h = op.h.borrow();
        let half = h.half_sphere();
        op.pool
            .install(|| {
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
) -> Result<(), ScfError> {
    r.kept.set_keep(true);
    let collected = std::mem::replace(&mut r.collect, false);
    let grid = dist.potential().len();
    let stored = r.kept.is_stored();
    if !(collected && stored && has_initial && r.v_old.len() == grid && x.cols() == bands) {
        r.kept.invalidate();
        crate::m2::count(crate::m2::Count::ReuseNotKept);
        return Ok(());
    }
    r.delta.clear();
    let mut changed = false;
    for (&new, &old) in dist.potential().iter().zip(&r.v_old) {
        let d = new - old;
        changed |= d != 0.;
        r.delta.push(d);
    }
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
        let before = state.as_ref().map(|s| s.op.h.as_ptr() as usize);
        let capacity = real_capacity(&options, 4, h.kinetic_diagonal().len());
        ensure_state(&h, &basis, 4, capacity, &options, &pool, &mut state).unwrap();
        assert_eq!(before, state.as_ref().map(|s| s.op.h.as_ptr() as usize));
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
            image: RealMatrix::zeros(state.op.n, 4).unwrap(),
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
