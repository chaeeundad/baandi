//! Exclusive nested phase profile of one SCF run (ticket
//! BREAKTHROUGH-0-INSTRUMENT, review 2026-11-17 section 13.1).
//!
//! `BAANDI_PHASE_PROFILE=on` (default off: no clock read, no record kept, the
//! arithmetic is never touched either way). Every rank thread keeps its own
//! accumulators (simulated ranks are threads of one process, so process-wide
//! atomics could not separate them); a report takes the MAXIMUM over ranks of
//! every entry (the critical path) and the minimum (the imbalance), never a sum
//! of CPU seconds.
//!
//! Top level, mutually exclusive buckets of the run wall:
//! `init` (problem construction and the atomic start, which the first map's
//! potential phase would otherwise carry), `coarse` (ladder coarse SCF maps and
//! their mixing), `fine` (ordinary full-cutoff maps and their mixing;
//! `fine_redo` lists the discarded first-step solves separately), `verify_pass`
//! and `verify_fail` (fixed-point verification maps by their outcome),
//! `materialize` (lazy orbital gathers and, when the caller reports it, output
//! writing). What is not covered by any bucket is the caller's remainder
//! (`wall - sum`, printed by the harnesses).
//!
//! Inside each map bucket: `potential`, `solve`, `density`, `energy`, `mixing`
//! (the SCF phase clocks; `potential` excludes the atomic start). Inside
//! `solve`, exclusive stages of the eigensolver: `solve.start`, `solve.apply`
//! (full H columns and potential-change columns, with the pair transform
//! counts), `solve.reduced_gram` (overlap / reduced border), `solve.reduced_eigen`
//! (factor / EVD), `solve.rotate`, `solve.residual`, `solve.orthonormalize`,
//! `solve.restart`, `solve.return`, and `solve.other` (the rest of the solve
//! wall: the correction loop, copies, waits outside stages). Collective seconds
//! and calls are listed beside them (they are INSIDE the stages, not exclusive).
//! Counts (maps, Davidson iterations, H columns, restarts, FFT pairs, collective
//! calls) are kept in count-only mode too (`BAANDI_STAGE_TIMERS=off`).
//!
//! The serial (threads-only) path reports its Davidson stage clocks under the
//! same `solve.*` names when `BAANDI_DAVIDSON_STAGES=on`; without it the serial
//! solve is one number.
use dft_mpi::Collective;
use dft_numerics::stage_clock::local::{
    self, COLLECTIVE, DIST_STAGE, M2_STAGE, OP_COUNT, SERIAL_STAGE, Snapshot,
};
/// Davidson lane diagnostic transport; serialization remains in the CLI.
pub use dft_solvers::real_davidson::lane_profile as davidson_lane;
use std::cell::{Cell, RefCell};
use std::sync::{Mutex, OnceLock};

/// Top-level phase in progress on the calling thread.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Top {
    Init = 0,
    Coarse = 1,
    Fine = 2,
    Verification = 3,
    Materialize = 4,
}

/// Exclusive buckets of the report.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Bucket {
    Init = 0,
    Coarse = 1,
    Fine = 2,
    FineRedo = 3,
    VerifyPass = 4,
    VerifyFail = 5,
    Materialize = 6,
    /// BREAKTHROUGH-4: named work between the maps that no other bucket covers
    /// (ladder hand-off, session set-up, stop rule and iteration record).
    Between = 7,
}
pub const BUCKETS: usize = 8;
pub const BUCKET_LABELS: [&str; BUCKETS] = [
    "init",
    "coarse",
    "fine",
    "fine_redo",
    "verify_pass",
    "verify_fail",
    "materialize",
    "between",
];

/// How a recorded map is filed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MapKind {
    Ordinary,
    Redo,
    VerifyPass,
    VerifyFail,
}

/// Map phase of the SCF clocks.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Phase {
    Potential,
    Solve,
    Density,
    Energy,
}

// Seconds slots. Slots 32 and above (BREAKTHROUGH-4-OTHER-AND-INIT): the solve
// time the earlier stages left in `solve.other` (solver `setup`,
// `precondition`, `append`; engine interface `kept_update`, `entry`, `finish`,
// `entry_wait`), the split of the rest into the solver's loop body
// (`other_solver` = solver call minus its stages) and the engine interface
// (`other_engine`), the wait-probe time inside the solver's collectives
// (`solve.wait`, not exclusive), the parts of `init` (problem construction and
// the atomic start) and the named parts of the `between` bucket.
pub const SEC_LABELS: [&str; NSEC] = [
    "wall",
    "potential",
    "solve",
    "density",
    "energy",
    "mixing",
    "atomic_start",
    "io",
    "solve.start",
    "solve.apply",
    "solve.reduced_gram",
    "solve.reduced_eigen",
    "solve.rotate",
    "solve.residual",
    "solve.orthonormalize",
    "solve.restart",
    "solve.return",
    "solve.other",
    "coll.allgatherv",
    "coll.tree_sum",
    "coll.bcast",
    "coll.alltoall",
    "coll.barrier",
    "m2.density",
    "m2.xc",
    "m2.nonlocal_energy",
    "m2.solve_prepare",
    "m2.solve_unpack",
    "m2.lazy_kinetic",
    "m2.lazy_materialize",
    "m2.atomic_start",
    "m2.screen",
    "solve.setup",
    "solve.precondition",
    "solve.append",
    "solve.kept_update",
    "solve.entry",
    "solve.finish",
    "solve.entry_wait",
    "solve.wait",
    "solve.solver_total",
    "solve.other_solver",
    "solve.other_engine",
    "init.cli",
    "init.basis",
    "init.ionic_local",
    "init.densities",
    "init.states",
    "init.atomic_prepare",
    "init.ewald",
    "init.rest",
    "atomic.span",
    "atomic.hq",
    "atomic.ritz",
    "between.ladder_handoff",
    "between.session_setup",
    "between.stop_rule",
    "between.verify_glue",
    "solve.real_state",
    "atomic.radial",
    "atomic.fill",
    "atomic.orthonormalize",
    "atomic.finish",
    "atomic.gram",
    "atomic.eigh",
    "atomic.rotate",
    "between.ladder_truncation",
    "nonlocal.projection",
    "nonlocal.reduce",
    "nonlocal.coupling",
    "nonlocal.back",
    "nonlocal.add",
    "nonlocal.audit",
];
const S_WALL: usize = 0;
const S_POTENTIAL: usize = 1;
const S_SOLVE: usize = 2;
const S_DENSITY: usize = 3;
const S_ENERGY: usize = 4;
const S_MIXING: usize = 5;
const S_ATOMIC: usize = 6;
const S_IO: usize = 7;
const S_SOLVE_START: usize = 8;
const S_SOLVE_APPLY: usize = 9;
const S_SOLVE_GRAM: usize = 10;
const S_SOLVE_EIGEN: usize = 11;
const S_SOLVE_ROTATE: usize = 12;
const S_SOLVE_RESIDUAL: usize = 13;
const S_SOLVE_ORTHO: usize = 14;
const S_SOLVE_RESTART: usize = 15;
const S_SOLVE_RETURN: usize = 16;
const S_SOLVE_OTHER: usize = 17;
const S_COLL: usize = 18;
const S_M2: usize = 23;
const S_M2_LAZY_MATERIALIZE: usize = 29;
const S_M2_SOLVE_PREPARE: usize = 26;
const S_M2_SOLVE_UNPACK: usize = 27;
const S_SOLVE_SETUP: usize = 32;
const S_SOLVE_PRECONDITION: usize = 33;
const S_SOLVE_APPEND: usize = 34;
const S_SOLVE_KEPT_UPDATE: usize = 35;
const S_SOLVE_ENTRY: usize = 36;
const S_SOLVE_FINISH: usize = 37;
const S_SOLVE_ENTRY_WAIT: usize = 38;
const S_SOLVE_WAIT: usize = 39;
const S_SOLVER_TOTAL: usize = 40;
const S_OTHER_SOLVER: usize = 41;
const S_OTHER_ENGINE: usize = 42;
const S_INIT_PART: usize = 43;
const S_ATOMIC_PART: usize = 51;
const S_BETWEEN_PART: usize = 54;
const BETWEEN_PARTS: usize = 4;
const S_SOLVE_REAL_STATE: usize = 58;
const S_ATOMIC_DETAIL: usize = 59;
const S_BETWEEN_TRUNCATION: usize = 66;
const S_NONLOCAL: usize = 67;
pub const NSEC: usize = 73;

/// Parts of problem construction (`init` bucket; seconds only, the bucket wall
/// is the caller's construction clock).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InitPart {
    /// Caller work before the engine (input parsing, preflight, basis
    /// preparation and pseudopotential reading of the CLI).
    Cli = 0,
    /// Wavefunction and density bases and FFT grid specifications.
    Basis = 1,
    /// Ionic local potential on the density sphere.
    IonicLocal = 2,
    /// FFT plan, core density, atomic start density and the real-space ionic
    /// potential.
    Densities = 3,
    /// Worker pool and the k states (basis, nonlocal projectors, Hamiltonian,
    /// band-density workspace).
    States = 4,
    /// Atomic-orbital radial tables of the atomic start.
    AtomicPrepare = 5,
    /// Ewald energy.
    Ewald = 6,
    /// Validation, worker FFT plans, Coulomb metric, force model.
    Rest = 7,
}
/// Parts of the atomic start (`init` bucket; seconds only). `Span`, `Hq` and
/// `Ritz` are the three steps; the others split them (span: radial tables,
/// column fill, orthonormalization, plane-wave completion; Ritz: H and the
/// Gram matrix, small eigenproblem, rotation to the start vectors).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AtomicPart {
    /// Atomic-orbital span on the plane-wave basis.
    Span,
    /// H applied to the span (threads path with ATOMIC-REAL only).
    Hq,
    /// Gram, small eigenproblem and rotation (Rayleigh-Ritz).
    Ritz,
    Radial,
    Fill,
    Orthonormalize,
    Finish,
    Gram,
    Eigh,
    Rotate,
}
impl AtomicPart {
    fn slot(self) -> usize {
        match self {
            Self::Span => S_ATOMIC_PART,
            Self::Hq => S_ATOMIC_PART + 1,
            Self::Ritz => S_ATOMIC_PART + 2,
            Self::Radial => S_ATOMIC_DETAIL,
            Self::Fill => S_ATOMIC_DETAIL + 1,
            Self::Orthonormalize => S_ATOMIC_DETAIL + 2,
            Self::Finish => S_ATOMIC_DETAIL + 3,
            Self::Gram => S_ATOMIC_DETAIL + 4,
            Self::Eigh => S_ATOMIC_DETAIL + 5,
            Self::Rotate => S_ATOMIC_DETAIL + 6,
        }
    }
}
/// Named parts of the `between` bucket (each adds to its wall).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BetweenPart {
    /// Cutoff ladder: orbital and density prolongation of the hand-off.
    LadderHandoff,
    /// Session construction and initialization (mixer validation, imported
    /// state, seeds).
    SessionSetup,
    /// Stop-rule estimate, density residual and the iteration record of an
    /// SCF step (outside its map).
    StopRule,
    /// Verification block outside its map (set-up and bookkeeping).
    VerifyGlue,
    /// Cutoff ladder: the report's atomic truncation distance (prolongation
    /// of the coarse atomic start density).
    LadderTruncation,
}
impl BetweenPart {
    fn slot(self) -> usize {
        match self {
            Self::LadderHandoff => S_BETWEEN_PART,
            Self::SessionSetup => S_BETWEEN_PART + 1,
            Self::StopRule => S_BETWEEN_PART + 2,
            Self::VerifyGlue => S_BETWEEN_PART + 3,
            Self::LadderTruncation => S_BETWEEN_TRUNCATION,
        }
    }
}

// Count slots.
pub const COUNT_LABELS: [&str; NCOUNT] = [
    "maps",
    "davidson_iterations",
    "h_columns",
    "h_delta_columns",
    "restarts",
    "contractions",
    "kept_starts",
    "fft_forward_pairs",
    "fft_inverse_pairs",
    "kept_plane_pairs",
    "op_columns",
    "coll.allgatherv.calls",
    "coll.tree_sum.calls",
    "coll.bcast.calls",
    "coll.alltoall.calls",
    "coll.barrier.calls",
    "solve.apply.calls",
    "solve.reduced_eigen.calls",
    "solve.restart.calls",
    "solve.orthonormalize.calls",
    "nonlocal.projection.calls",
    "nonlocal.reduce.calls",
    "nonlocal.coupling.calls",
    "nonlocal.back.calls",
    "nonlocal.add.calls",
    "nonlocal.audit.calls",
    "nonlocal.projection.flops",
    "nonlocal.reduce.flops",
    "nonlocal.coupling.flops",
    "nonlocal.back.flops",
    "nonlocal.add.flops",
    "nonlocal.audit.flops",
    "nonlocal.projection.bytes",
    "nonlocal.reduce.bytes",
    "nonlocal.coupling.bytes",
    "nonlocal.back.bytes",
    "nonlocal.add.bytes",
    "nonlocal.audit.bytes",
    "nonlocal.rows_sum",
    "nonlocal.projectors_sum",
    "nonlocal.bands_sum",
    "nonlocal.width1",
    "nonlocal.width2",
    "nonlocal.width3_8",
    "nonlocal.width9_32",
    "nonlocal.width33_plus",
];
const C_MAPS: usize = 0;
const C_ITER: usize = 1;
const C_HCOLS: usize = 2;
const C_HDELTA: usize = 3;
const C_RESTARTS: usize = 4;
const C_CONTRACTIONS: usize = 5;
const C_KEPT: usize = 6;
const C_FFT_FWD: usize = 7;
const C_FFT_INV: usize = 8;
const C_KEPT_PAIRS: usize = 9;
const C_OP_COLS: usize = 10;
const C_COLL: usize = 11;
const C_APPLY_CALLS: usize = 16;
const C_EIGEN_CALLS: usize = 17;
const C_RESTART_CALLS: usize = 18;
const C_ORTHO_CALLS: usize = 19;
const C_NONLOCAL: usize = 20;
pub const NCOUNT: usize = 46;

/// Accumulated seconds and counts of one bucket.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Accum {
    pub seconds: [f64; NSEC],
    pub counts: [u64; NCOUNT],
    /// LANE-FFT: separate optional lane accounting; not part of wall closure.
    pub fft: dft_numerics::real::fft_profile::Snapshot,
}
impl Default for Accum {
    fn default() -> Self {
        Self {
            seconds: [0.0; NSEC],
            counts: [0; NCOUNT],
            fft: dft_numerics::real::fft_profile::Snapshot::default(),
        }
    }
}
impl Accum {
    fn add(&mut self, other: &Accum) {
        self.fft.add(other.fft);
        for (a, b) in self.seconds.iter_mut().zip(&other.seconds) {
            *a += b;
        }
        for (a, b) in self.counts.iter_mut().zip(&other.counts) {
            *a += b;
        }
    }
    /// Attribute a per-thread counter delta: solver stages, collectives, M2
    /// stages, operator counts. Lazy materializations are returned separately
    /// (they belong to the `materialize` bucket wherever they happen).
    fn add_delta(&mut self, d: &Snapshot) -> f64 {
        self.fft.add(d.fft);
        let ns = |i: usize| d.nanos[i] as f64 * 1e-9;
        // Distributed solver stages (STAGES order: reduced_gram, reduced_eigen,
        // rotate, residual, orthonormalize, apply_total, apply_operator,
        // other_solve, solve_start, solve_return, restart).
        self.seconds[S_SOLVE_GRAM] += ns(DIST_STAGE);
        self.seconds[S_SOLVE_EIGEN] += ns(DIST_STAGE + 1);
        self.seconds[S_SOLVE_ROTATE] += ns(DIST_STAGE + 2);
        self.seconds[S_SOLVE_RESIDUAL] += ns(DIST_STAGE + 3);
        self.seconds[S_SOLVE_ORTHO] += ns(DIST_STAGE + 4);
        self.seconds[S_SOLVE_APPLY] += ns(DIST_STAGE + 5);
        self.seconds[S_SOLVE_START] += ns(DIST_STAGE + 8);
        self.seconds[S_SOLVE_RETURN] += ns(DIST_STAGE + 9);
        self.seconds[S_SOLVE_RESTART] += ns(DIST_STAGE + 10);
        self.counts[C_APPLY_CALLS] += d.calls[DIST_STAGE + 5];
        self.counts[C_EIGEN_CALLS] += d.calls[DIST_STAGE + 1];
        self.counts[C_RESTART_CALLS] += d.calls[DIST_STAGE + 10];
        self.counts[C_ORTHO_CALLS] += d.calls[DIST_STAGE + 4];
        // BREAKTHROUGH-4 solver labels (setup, precondition, append, wait,
        // solver_total) and engine interface stages (m2 kept_update,
        // solve_entry, solve_finish, entry_wait).
        self.seconds[S_SOLVE_SETUP] += ns(DIST_STAGE + 11);
        self.seconds[S_SOLVE_PRECONDITION] += ns(DIST_STAGE + 12);
        self.seconds[S_SOLVE_APPEND] += ns(DIST_STAGE + 13);
        self.seconds[S_SOLVE_WAIT] += ns(DIST_STAGE + 14);
        self.seconds[S_SOLVER_TOTAL] += ns(DIST_STAGE + 15);
        self.seconds[S_SOLVE_KEPT_UPDATE] += ns(M2_STAGE + 10);
        self.seconds[S_SOLVE_ENTRY] += ns(M2_STAGE + 11);
        self.seconds[S_SOLVE_FINISH] += ns(M2_STAGE + 12);
        self.seconds[S_SOLVE_ENTRY_WAIT] += ns(M2_STAGE + 13);
        self.seconds[S_SOLVE_REAL_STATE] += ns(M2_STAGE + 14);
        // Serial Davidson stages (Stage order: reduced, eigen, rotate, residual,
        // orthonormalize, h_apply, sketch, precondition), only with
        // BAANDI_DAVIDSON_STAGES=on.
        self.seconds[S_SOLVE_GRAM] += ns(SERIAL_STAGE);
        self.seconds[S_SOLVE_EIGEN] += ns(SERIAL_STAGE + 1);
        self.seconds[S_SOLVE_ROTATE] += ns(SERIAL_STAGE + 2);
        self.seconds[S_SOLVE_RESIDUAL] += ns(SERIAL_STAGE + 3);
        self.seconds[S_SOLVE_ORTHO] += ns(SERIAL_STAGE + 4);
        self.seconds[S_SOLVE_APPLY] += ns(SERIAL_STAGE + 5);
        self.seconds[S_SOLVE_OTHER] += ns(SERIAL_STAGE + 6) + ns(SERIAL_STAGE + 7);
        self.counts[C_APPLY_CALLS] += d.calls[SERIAL_STAGE + 5];
        self.counts[C_EIGEN_CALLS] += d.calls[SERIAL_STAGE + 1];
        // Collectives.
        for i in 0..5 {
            self.seconds[S_COLL + i] += ns(COLLECTIVE + i);
            self.counts[C_COLL + i] += d.calls[COLLECTIVE + i];
        }
        // M2 engine stages (m2::Stage order: density, xc, nonlocal_energy,
        // screen_integrals, screen_clusters, atomic_start, solve_prepare,
        // solve_unpack, lazy_materialize, lazy_kinetic).
        self.seconds[S_M2] += ns(M2_STAGE);
        self.seconds[S_M2 + 1] += ns(M2_STAGE + 1);
        self.seconds[S_M2 + 2] += ns(M2_STAGE + 2);
        self.seconds[S_M2 + 3] += ns(M2_STAGE + 6);
        self.seconds[S_M2 + 4] += ns(M2_STAGE + 7);
        self.seconds[S_M2 + 5] += ns(M2_STAGE + 9);
        self.seconds[S_M2 + 7] += ns(M2_STAGE + 5);
        self.seconds[S_M2 + 8] += ns(M2_STAGE + 3) + ns(M2_STAGE + 4);
        // Operator counts.
        self.counts[C_OP_COLS] += d.calls[OP_COUNT];
        self.counts[C_HDELTA] += d.calls[OP_COUNT + 1];
        self.counts[C_FFT_FWD] += d.calls[OP_COUNT + 2];
        self.counts[C_FFT_INV] += d.calls[OP_COUNT + 3];
        self.counts[C_KEPT_PAIRS] += d.calls[OP_COUNT + 4];
        ns(M2_STAGE + 8)
    }
}

/// One rank's buckets.
#[derive(Clone, Debug, PartialEq)]
pub struct RankProfile {
    pub rank: usize,
    pub buckets: [Accum; BUCKETS],
    /// RUN-REPORT: this rank's own external run wall (seconds), set by the
    /// caller before `gather`; 0 = not reported (no closure residual).
    pub wall: f64,
}

// LANE-NONLOCAL instrumentation wiring. Kept separate from the shared solver
// snapshot slots so other lanes can append their own fields independently.
struct NonlocalMark(dft_physics::distributed_real_gamma::lane_profile::Snapshot);
impl NonlocalMark {
    fn new() -> Self {
        Self(dft_physics::distributed_real_gamma::lane_profile::snapshot())
    }
    fn flush_into(&mut self, target: &mut Accum) {
        use dft_physics::distributed_real_gamma::lane_profile::{self, SHAPES, STAGES};
        let now = lane_profile::snapshot();
        let delta = now.since(self.0);
        self.0 = now;
        for i in 0..STAGES {
            target.seconds[S_NONLOCAL + i] += delta.nanos[i] as f64 * 1e-9;
            target.counts[C_NONLOCAL + i] += delta.calls[i];
            target.counts[C_NONLOCAL + STAGES + i] += delta.flops[i];
            target.counts[C_NONLOCAL + 2 * STAGES + i] += delta.bytes[i];
        }
        for i in 0..SHAPES {
            target.counts[C_NONLOCAL + 3 * STAGES + i] += delta.shapes[i];
        }
    }
}

struct State {
    top: Top,
    buckets: [Accum; BUCKETS],
    /// The map in progress (committed by the session or by the next map).
    pending: Option<Accum>,
    /// Pending map's atomic-start seconds (excluded from its potential).
    pending_atomic: f64,
    mark: Snapshot,
    nonlocal_mark: NonlocalMark,
}
impl State {
    fn new() -> Self {
        Self {
            top: Top::Fine,
            buckets: [Accum::default(); BUCKETS],
            pending: None,
            pending_atomic: 0.0,
            mark: local::snapshot(),
            nonlocal_mark: NonlocalMark::new(),
        }
    }
    fn ordinary_bucket(&self) -> Bucket {
        match self.top {
            Top::Init => Bucket::Init,
            Top::Coarse => Bucket::Coarse,
            Top::Fine => Bucket::Fine,
            Top::Verification => Bucket::VerifyPass,
            Top::Materialize => Bucket::Materialize,
        }
    }
    /// Attribute the thread's counter delta since the last mark to `bucket`.
    fn flush(&mut self, bucket: Bucket) {
        let now = local::snapshot();
        let d = now.since(&self.mark);
        self.mark = now;
        let target = &mut self.buckets[bucket as usize];
        self.nonlocal_mark.flush_into(target);
        let lazy = target.add_delta(&d);
        let m = &mut self.buckets[Bucket::Materialize as usize];
        m.seconds[S_M2_LAZY_MATERIALIZE] += lazy;
        m.seconds[S_IO] += lazy;
        m.seconds[S_WALL] += lazy;
    }
    fn commit(&mut self, kind: MapKind) {
        let Some(mut map) = self.pending.take() else {
            return;
        };
        let bucket = match kind {
            MapKind::Ordinary => self.ordinary_bucket(),
            MapKind::Redo => Bucket::FineRedo,
            MapKind::VerifyPass => Bucket::VerifyPass,
            MapKind::VerifyFail => Bucket::VerifyFail,
        };
        map.seconds[S_POTENTIAL] = (map.seconds[S_POTENTIAL] - self.pending_atomic).max(0.0);
        self.pending_atomic = 0.0;
        map.seconds[S_WALL] = map.seconds[S_POTENTIAL]
            + map.seconds[S_SOLVE]
            + map.seconds[S_DENSITY]
            + map.seconds[S_ENERGY];
        map.counts[C_MAPS] = 1;
        self.buckets[bucket as usize].add(&map);
        self.flush(bucket);
    }
}

thread_local! {
    static STATE: RefCell<Option<State>> = const { RefCell::new(None) };
    static SOLVE_MARK: Cell<Option<Snapshot>> = const { Cell::new(None) };
}
static REGISTRY: Mutex<Vec<RankProfile>> = Mutex::new(Vec::new());

// 0: the switch; 1: forced on; 2: forced off (tests).
static FORCED: std::sync::atomic::AtomicU8 = std::sync::atomic::AtomicU8::new(0);
/// Tests: force the profile on or off for the process (`None`: the switch).
pub fn force(mode: Option<bool>) {
    FORCED.store(
        match mode {
            None => 0,
            Some(true) => 1,
            Some(false) => 2,
        },
        std::sync::atomic::Ordering::Relaxed,
    );
    if mode == Some(true) {
        local::enable(true);
    }
}
/// `BAANDI_PHASE_PROFILE=on`, read once; turns the per-thread counter mirror on.
pub fn enabled() -> bool {
    match FORCED.load(std::sync::atomic::Ordering::Relaxed) {
        1 => return true,
        2 => return false,
        _ => {}
    }
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| {
        let on = dft_numerics::switches::on(dft_numerics::switches::Id::PHASE_PROFILE);
        if on {
            local::enable(true);
        }
        on
    })
}

fn with<R>(f: impl FnOnce(&mut State) -> R) -> Option<R> {
    if !enabled() {
        return None;
    }
    STATE.with(|s| {
        let mut s = s.borrow_mut();
        let s = s.get_or_insert_with(State::new);
        Some(f(s))
    })
}

/// Restores the previous top-level phase on drop.
pub struct Guard {
    previous: Option<Top>,
}
impl Drop for Guard {
    fn drop(&mut self) {
        if let Some(p) = self.previous {
            with(|s| {
                let leaving = s.ordinary_bucket();
                s.flush(leaving);
                s.top = p;
            });
        }
    }
}
/// Enter a top-level phase on this thread until the guard is dropped.
#[must_use]
pub fn enter(top: Top) -> Guard {
    let previous = with(|s| {
        let leaving = s.ordinary_bucket();
        s.flush(leaving);
        std::mem::replace(&mut s.top, top)
    });
    Guard { previous }
}

/// A map starts: any uncommitted map is filed as ordinary.
pub(crate) fn map_begin() {
    with(|s| {
        s.commit(MapKind::Ordinary);
        s.pending = Some(Accum::default());
    });
}
/// Seconds of one phase of the pending map.
pub(crate) fn map_phase(phase: Phase, seconds: f64) {
    with(|s| {
        if let Some(m) = s.pending.as_mut() {
            let i = match phase {
                Phase::Potential => S_POTENTIAL,
                Phase::Solve => S_SOLVE,
                Phase::Density => S_DENSITY,
                Phase::Energy => S_ENERGY,
            };
            m.seconds[i] += seconds;
        }
    });
}
/// File the pending map.
pub fn map_commit(kind: MapKind) {
    with(|s| s.commit(kind));
}
/// The eigensolve of the pending map begins (the thread's counters are marked).
pub(crate) fn solve_begin() {
    if dft_solvers::real_davidson::lane_profile::requested() {
        let c = context();
        dft_solvers::real_davidson::lane_profile::map(c.stage, c.kind, c.iteration);
    }
    if enabled() {
        SOLVE_MARK.with(|m| m.set(Some(local::snapshot())));
    }
}
/// Restarts, contractions and the kept-start flag of the latest real Gamma
/// solve, keyed by rank (the solver entry notes them, the map reads them).
static SOLVE_NOTES: Mutex<Vec<[u64; 3]>> = Mutex::new(Vec::new());
pub(crate) fn note_solve(restarts: usize, contractions: usize, kept_taken: bool) {
    if !enabled() {
        return;
    }
    let r = rank_key();
    let mut v = SOLVE_NOTES.lock().expect("solve notes");
    if v.len() <= r {
        v.resize(r + 1, [0; 3]);
    }
    let n = &mut v[r];
    n[0] += restarts as u64;
    n[1] += contractions as u64;
    n[2] += u64::from(kept_taken);
}
fn take_solve_note() -> [u64; 3] {
    let r = rank_key();
    let mut v = SOLVE_NOTES.lock().expect("solve notes");
    v.get_mut(r).map_or([0; 3], std::mem::take)
}
/// The eigensolve of the pending map ended: its stage deltas and work.
pub(crate) fn solve_end(iterations: usize, h_columns: usize) {
    let Some(mark) = SOLVE_MARK.with(|m| m.take()) else {
        return;
    };
    let [restarts, contractions, kept] = take_solve_note();
    with(|s| {
        let now = local::snapshot();
        let d = now.since(&mark);
        s.mark = now;
        if let Some(m) = s.pending.as_mut() {
            let lazy = m.add_delta(&d);
            let mat = &mut s.buckets[Bucket::Materialize as usize];
            mat.seconds[S_M2_LAZY_MATERIALIZE] += lazy;
            mat.seconds[S_IO] += lazy;
            mat.seconds[S_WALL] += lazy;
            m.counts[C_ITER] += iterations as u64;
            m.counts[C_HCOLS] += h_columns as u64;
            m.counts[C_RESTARTS] += restarts;
            m.counts[C_CONTRACTIONS] += contractions;
            m.counts[C_KEPT] += kept;
        }
    });
}
/// Mixing seconds of the current ordinary bucket.
pub(crate) fn add_mixing(seconds: f64) {
    with(|s| {
        let b = s.ordinary_bucket();
        let a = &mut s.buckets[b as usize];
        a.seconds[S_MIXING] += seconds;
        a.seconds[S_WALL] += seconds;
    });
}
/// Atomic start seconds: filed under `init`, excluded from the pending map's
/// potential phase.
pub(crate) fn add_atomic_start(seconds: f64) {
    with(|s| {
        let a = &mut s.buckets[Bucket::Init as usize];
        a.seconds[S_ATOMIC] += seconds;
        a.seconds[S_WALL] += seconds;
        if s.pending.is_some() {
            s.pending_atomic += seconds;
        }
    });
}
/// Problem construction seconds (caller-reported), under `init`.
pub fn add_init(seconds: f64) {
    with(|s| {
        let a = &mut s.buckets[Bucket::Init as usize];
        a.seconds[S_WALL] += seconds;
    });
}
/// A clock for a profile part: `Some(now)` only when the profile is on (no
/// clock read otherwise).
pub fn clock() -> Option<std::time::Instant> {
    enabled().then(std::time::Instant::now)
}
fn since(clock: Option<std::time::Instant>) -> f64 {
    clock.map_or(0.0, |c| c.elapsed().as_secs_f64())
}
/// Seconds of one part of problem construction since `clock` (under `init`,
/// not added to its wall: the caller's construction clock is the wall).
pub fn init_part(part: InitPart, clock: Option<std::time::Instant>) {
    if clock.is_none() {
        return;
    }
    let s = since(clock);
    with(|st| st.buckets[Bucket::Init as usize].seconds[S_INIT_PART + part as usize] += s);
}
/// Seconds of one part of the atomic start since `clock` (under `init`; the
/// atomic start total is its wall).
pub(crate) fn atomic_part(part: AtomicPart, clock: Option<std::time::Instant>) {
    if clock.is_none() {
        return;
    }
    let s = since(clock);
    with(|st| st.buckets[Bucket::Init as usize].seconds[part.slot()] += s);
}
/// Seconds of one named part of the `between` bucket since `clock` (added to
/// that bucket's wall).
pub fn between(part: BetweenPart, clock: Option<std::time::Instant>) {
    if clock.is_none() {
        return;
    }
    let s = since(clock);
    with(|st| {
        let a = &mut st.buckets[Bucket::Between as usize];
        a.seconds[part.slot()] += s;
        a.seconds[S_WALL] += s;
    });
}
/// Output / materialization seconds (caller-reported), under `materialize`.
pub fn add_io(seconds: f64) {
    with(|s| {
        let a = &mut s.buckets[Bucket::Materialize as usize];
        a.seconds[S_IO] += seconds;
        a.seconds[S_WALL] += seconds;
    });
}

/// End of this rank's run: commit the pending map, file the remaining counter
/// delta under the current phase, publish the profile to the process registry
/// (simulated ranks) and return it. `None` when the profile is off.
pub fn finish_rank(rank: usize) -> Option<RankProfile> {
    let p = with(|s| {
        s.commit(MapKind::Ordinary);
        let b = s.ordinary_bucket();
        s.flush(b);
        let p = RankProfile {
            rank,
            buckets: s.buckets,
            wall: 0.0,
        };
        *s = State::new();
        p
    })?;
    REGISTRY.lock().expect("profile registry").push(p.clone());
    Some(p)
}
/// Profiles published by `finish_rank` since the last call (simulated ranks of
/// one process), sorted by rank.
pub fn take_registry() -> Vec<RankProfile> {
    let mut v = std::mem::take(&mut *REGISTRY.lock().expect("profile registry"));
    v.sort_by_key(|p| p.rank);
    v
}
/// Real MPI: every rank's profile (one allgatherv), sorted by rank.
pub fn gather<C: Collective + ?Sized>(
    comm: &C,
    mine: &RankProfile,
) -> Result<Vec<RankProfile>, String> {
    let fft_words = if fft_lane_enabled() {
        FFT_WIRE_WORDS
    } else {
        0
    };
    let per = 2 + BUCKETS * (NSEC + NCOUNT + fft_words);
    let mut local = Vec::with_capacity(per);
    local.push(mine.rank as f64);
    local.push(mine.wall);
    for b in &mine.buckets {
        local.extend_from_slice(&b.seconds);
        local.extend(b.counts.iter().map(|&c| c as f64));
        if fft_words != 0 {
            for &v in b.fft.0.iter().flatten() {
                local.extend([u64::from(v as u32) as f64, (v >> 32) as f64]);
            }
        }
    }
    let ranks = comm.size();
    let mut all = vec![0.0; per * ranks];
    comm.allgatherv_f64(&local, &vec![per; ranks], &mut all)
        .map_err(|e| e.to_string())?;
    let mut out: Vec<RankProfile> = all
        .chunks(per)
        .map(|c| {
            let mut p = RankProfile {
                rank: c[0] as usize,
                buckets: [Accum::default(); BUCKETS],
                wall: c[1],
            };
            let mut at = 2;
            for b in &mut p.buckets {
                b.seconds.copy_from_slice(&c[at..at + NSEC]);
                at += NSEC;
                for (k, v) in b.counts.iter_mut().zip(&c[at..at + NCOUNT]) {
                    *k = *v as u64;
                }
                at += NCOUNT;
                if fft_words != 0 {
                    for v in b.fft.0.iter_mut().flatten() {
                        *v = c[at] as u64 | ((c[at + 1] as u64) << 32);
                        at += 2;
                    }
                }
            }
            p
        })
        .collect();
    out.sort_by_key(|p| p.rank);
    Ok(out)
}

/// `solve.other` of one rank's bucket: the solve wall not covered by the
/// exclusive stages (solver stages, the BREAKTHROUGH-4 solver labels and the
/// engine's solve interface stages), split into the solver's uncovered loop
/// body (`other_solver`: the solver call minus its stages) and the rest
/// (`other_engine`). The serial path's sketch/precondition seconds already in
/// `solve.other` count as covered.
fn derive_other(a: &mut Accum) {
    let solver_parts: f64 = (S_SOLVE_START..=S_SOLVE_RETURN)
        .chain([S_SOLVE_SETUP, S_SOLVE_PRECONDITION, S_SOLVE_APPEND])
        .map(|i| a.seconds[i])
        .sum();
    let engine_parts: f64 = [
        S_M2_SOLVE_PREPARE,
        S_M2_SOLVE_UNPACK,
        S_SOLVE_KEPT_UPDATE,
        S_SOLVE_ENTRY,
        S_SOLVE_FINISH,
        S_SOLVE_ENTRY_WAIT,
        S_SOLVE_REAL_STATE,
    ]
    .iter()
    .map(|&i| a.seconds[i])
    .sum();
    let parts = solver_parts + engine_parts + a.seconds[S_SOLVE_OTHER];
    let other = (a.seconds[S_SOLVE] - parts).max(0.0);
    a.seconds[S_SOLVE_OTHER] = other;
    let in_solver = if a.seconds[S_SOLVER_TOTAL] > 0.0 {
        (a.seconds[S_SOLVER_TOTAL] - solver_parts).clamp(0.0, other)
    } else {
        0.0
    };
    a.seconds[S_OTHER_SOLVER] = in_solver;
    a.seconds[S_OTHER_ENGINE] = other - in_solver;
}

// LANE-FFT Phase 0: independent section to avoid changes to shared slot lists.
const FFT_WIRE_WORDS: usize =
    2 * dft_numerics::real::fft_profile::STAGES * dft_numerics::real::fft_profile::METRICS;
fn fft_lane_enabled() -> bool {
    dft_numerics::switches::RunConfig::current().on(dft_numerics::switches::Id::FFT_PROFILE)
}
#[derive(Clone, Debug, PartialEq)]
pub struct FftRankProfile {
    pub rank: usize,
    pub maps: [u64; BUCKETS],
    pub buckets: [dft_numerics::real::fft_profile::Snapshot; BUCKETS],
}

/// One report entry: (label, maximum over ranks, minimum over ranks).
pub type Entry = (&'static str, f64, f64);

/// Maximum and minimum over ranks of every entry of every bucket.
#[derive(Clone, Debug, PartialEq)]
pub struct Report {
    pub ranks: usize,
    pub max: [Accum; BUCKETS],
    pub min: [Accum; BUCKETS],
    /// RUN-REPORT: each rank's exclusive tree closed on that rank alone.
    pub per_rank: Vec<RankClosure>,
    pub fft_ranks: Vec<FftRankProfile>,
}

/// One rank's closed exclusive tree: the sum of its own bucket walls
/// (`covered`, every bucket of the SAME rank) against its own external wall.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RankClosure {
    pub rank: usize,
    /// External wall of the rank (None when the caller did not report it).
    pub wall: Option<f64>,
    pub covered: f64,
    /// `wall - covered` (None without a wall): the rank's unattributed time.
    /// Negative beyond clock granularity means double counting.
    pub residual: Option<f64>,
}

/// Closure summary over ranks (each figure from per-rank closed trees, never a
/// sum of maxima taken on different ranks).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Closure {
    pub covered_max: f64,
    pub covered_min: f64,
    /// Rank with the largest covered time (the critical rank).
    pub critical_rank: usize,
    /// (max - min) / max of the covered time, 0 for one rank.
    pub imbalance: f64,
    pub wall_max: Option<f64>,
    pub wall_min: Option<f64>,
    pub residual_max: Option<f64>,
    pub residual_min: Option<f64>,
}
impl Report {
    pub fn from_ranks(profiles: &[RankProfile]) -> Self {
        let mut max = [Accum::default(); BUCKETS];
        let mut min = [Accum::default(); BUCKETS];
        // BREAKTHROUGH-4: the derived `solve.other` entries are formed per rank
        // first (a sum of per-stage maxima over different ranks can exceed any
        // rank's solve), then reduced like every other entry.
        let mut derived: Vec<RankProfile> = profiles.to_vec();
        for p in &mut derived {
            for a in &mut p.buckets {
                derive_other(a);
            }
        }
        let profiles = &derived[..];
        for b in 0..BUCKETS {
            for i in 0..NSEC {
                let it = profiles.iter().map(|p| p.buckets[b].seconds[i]);
                max[b].seconds[i] = it.clone().fold(0.0, f64::max);
                min[b].seconds[i] = it.fold(f64::INFINITY, f64::min);
                if !min[b].seconds[i].is_finite() {
                    min[b].seconds[i] = 0.0;
                }
            }
            for i in 0..NCOUNT {
                let it = profiles.iter().map(|p| p.buckets[b].counts[i]);
                max[b].counts[i] = it.clone().max().unwrap_or(0);
                min[b].counts[i] = it.min().unwrap_or(0);
            }
        }
        let per_rank = profiles
            .iter()
            .map(|p| {
                let covered: f64 = p.buckets.iter().map(|a| a.seconds[S_WALL]).sum();
                let wall = (p.wall > 0.0).then_some(p.wall);
                RankClosure {
                    rank: p.rank,
                    wall,
                    covered,
                    residual: wall.map(|w| w - covered),
                }
            })
            .collect();
        Self {
            ranks: profiles.len(),
            max,
            min,
            per_rank,
            fft_ranks: if fft_lane_enabled() {
                profiles
                    .iter()
                    .map(|p| FftRankProfile {
                        rank: p.rank,
                        maps: std::array::from_fn(|b| p.buckets[b].counts[C_MAPS]),
                        buckets: std::array::from_fn(|b| p.buckets[b].fft),
                    })
                    .collect()
            } else {
                Vec::new()
            },
        }
    }
    /// Sum of the bucket walls (maximum over ranks per bucket). Legacy figure:
    /// the maxima of different buckets can belong to different ranks, so this
    /// is not a closed tree; use [`Report::closure`].
    pub fn covered_wall(&self) -> f64 {
        self.max.iter().map(|a| a.seconds[S_WALL]).sum()
    }
    /// Per-rank closure summary (None without ranks).
    pub fn closure(&self) -> Option<Closure> {
        let first = self.per_rank.first()?;
        let mut c = Closure {
            covered_max: first.covered,
            covered_min: first.covered,
            critical_rank: first.rank,
            imbalance: 0.0,
            wall_max: first.wall,
            wall_min: first.wall,
            residual_max: first.residual,
            residual_min: first.residual,
        };
        let opt = |a: Option<f64>, b: Option<f64>, f: fn(f64, f64) -> f64| match (a, b) {
            (Some(a), Some(b)) => Some(f(a, b)),
            _ => None,
        };
        for r in &self.per_rank[1..] {
            if r.covered > c.covered_max {
                c.covered_max = r.covered;
                c.critical_rank = r.rank;
            }
            c.covered_min = c.covered_min.min(r.covered);
            c.wall_max = opt(c.wall_max, r.wall, f64::max);
            c.wall_min = opt(c.wall_min, r.wall, f64::min);
            c.residual_max = opt(c.residual_max, r.residual, f64::max);
            c.residual_min = opt(c.residual_min, r.residual, f64::min);
        }
        if c.covered_max > 0.0 {
            c.imbalance = (c.covered_max - c.covered_min) / c.covered_max;
        }
        Some(c)
    }
    /// Compact text table (seconds are the maximum over ranks, `[min]` the
    /// minimum; `run_wall` gives the remainder line).
    pub fn table(&self, run_wall: Option<f64>) -> String {
        use std::fmt::Write;
        let mut out = String::new();
        let _ = writeln!(
            out,
            "phase profile ({} rank{}; seconds = max over ranks [min over ranks]; counts = max over ranks)",
            self.ranks,
            if self.ranks == 1 { "" } else { "s" }
        );
        let _ = writeln!(
            out,
            "  {:<12} {:>10} {:>10} {:>5} {:>9} {:>9} {:>9} {:>9} {:>9}",
            "bucket", "wall", "[min]", "maps", "potential", "solve", "density", "energy", "mixing"
        );
        for (b, label) in BUCKET_LABELS.iter().enumerate() {
            let (mx, mn) = (&self.max[b], &self.min[b]);
            if mx.seconds[S_WALL] == 0.0 && mx.counts[C_MAPS] == 0 {
                continue;
            }
            let _ = writeln!(
                out,
                "  {:<12} {:>10.4} {:>10.4} {:>5} {:>9.4} {:>9.4} {:>9.4} {:>9.4} {:>9.4}",
                label,
                mx.seconds[S_WALL],
                mn.seconds[S_WALL],
                mx.counts[C_MAPS],
                mx.seconds[S_POTENTIAL],
                mx.seconds[S_SOLVE],
                mx.seconds[S_DENSITY],
                mx.seconds[S_ENERGY],
                mx.seconds[S_MIXING]
            );
            if b == Bucket::Init as usize {
                let _ = writeln!(
                    out,
                    "  {:<12} atomic start {:.4} s (span {:.4} = radial {:.4} + fill {:.4} + orthonormalize {:.4} + finish {:.4}; H span {:.4}; Ritz {:.4} = gram {:.4} + eigh {:.4} + rotate {:.4})",
                    "",
                    mx.seconds[S_ATOMIC],
                    mx.seconds[S_ATOMIC_PART],
                    mx.seconds[S_ATOMIC_DETAIL],
                    mx.seconds[S_ATOMIC_DETAIL + 1],
                    mx.seconds[S_ATOMIC_DETAIL + 2],
                    mx.seconds[S_ATOMIC_DETAIL + 3],
                    mx.seconds[S_ATOMIC_PART + 1],
                    mx.seconds[S_ATOMIC_PART + 2],
                    mx.seconds[S_ATOMIC_DETAIL + 4],
                    mx.seconds[S_ATOMIC_DETAIL + 5],
                    mx.seconds[S_ATOMIC_DETAIL + 6]
                );
                let parts: Vec<String> = (S_INIT_PART..S_ATOMIC_PART)
                    .map(|i| {
                        format!(
                            "{} {:.4} [{:.4}]",
                            &SEC_LABELS[i][5..],
                            mx.seconds[i],
                            mn.seconds[i]
                        )
                    })
                    .collect();
                let _ = writeln!(out, "  {:<12} construction: {}", "", parts.join(", "));
            }
            if b == Bucket::Between as usize {
                let parts: Vec<String> = (S_BETWEEN_PART..S_BETWEEN_PART + BETWEEN_PARTS)
                    .chain([S_BETWEEN_TRUNCATION])
                    .map(|i| format!("{} {:.4}", &SEC_LABELS[i][8..], mx.seconds[i]))
                    .collect();
                let _ = writeln!(out, "  {:<12} {}", "", parts.join(", "));
            }
            if b == Bucket::Materialize as usize {
                let _ = writeln!(
                    out,
                    "  {:<12} lazy materialize {:.4} s, io {:.4} s",
                    "",
                    mx.seconds[S_M2_LAZY_MATERIALIZE],
                    mx.seconds[S_IO] - mx.seconds[S_M2_LAZY_MATERIALIZE]
                );
            }
        }
        let covered = self.covered_wall();
        match run_wall {
            Some(w) => {
                let _ = writeln!(
                    out,
                    "  covered {:.4} s of run wall {:.4} s; remainder {:.4} s ({:.1} %)",
                    covered,
                    w,
                    w - covered,
                    100.0 * (w - covered) / w.max(1e-300)
                );
            }
            None => {
                let _ = writeln!(out, "  covered {:.4} s", covered);
            }
        }
        let _ = writeln!(
            out,
            "  solve stages (exclusive; seconds max [min]; counts max): bucket | start apply reduced_gram reduced_eigen rotate residual orthonormalize restart return other | Davidson, H columns (full, delta-V), restarts, contractions, kept starts, FFT pairs (fwd, inv), tree sums, alltoalls"
        );
        for (b, label) in BUCKET_LABELS.iter().enumerate() {
            let (mx, mn) = (&self.max[b], &self.min[b]);
            if mx.seconds[S_SOLVE] == 0.0 && mx.counts[C_ITER] == 0 {
                continue;
            }
            let stages: Vec<String> = (S_SOLVE_START..=S_SOLVE_OTHER)
                .map(|i| format!("{:.4} [{:.4}]", mx.seconds[i], mn.seconds[i]))
                .collect();
            let _ = writeln!(
                out,
                "  {:<12} | {} | {} iter, H {} ({} full, {} delta), restarts {}, contractions {}, kept {}, FFT pairs {} fwd {} inv, tree sums {} ({:.4} s), alltoalls {} ({:.4} s)",
                label,
                stages.join(" "),
                mx.counts[C_ITER],
                mx.counts[C_HCOLS] + mx.counts[C_HDELTA],
                mx.counts[C_HCOLS],
                mx.counts[C_HDELTA],
                mx.counts[C_RESTARTS],
                mx.counts[C_CONTRACTIONS],
                mx.counts[C_KEPT],
                mx.counts[C_FFT_FWD],
                mx.counts[C_FFT_INV],
                mx.counts[C_COLL + 1],
                mx.seconds[S_COLL + 1],
                mx.counts[C_COLL + 3],
                mx.seconds[S_COLL + 3]
            );
        }
        let _ = writeln!(
            out,
            "  solve labels (BREAKTHROUGH-4; exclusive except wait; seconds max [min]): bucket | setup precondition append | real_state prepare entry kept_update finish unpack entry_wait | other = other_solver + other_engine | wait (inside stages), barrier"
        );
        for (b, label) in BUCKET_LABELS.iter().enumerate() {
            let (mx, mn) = (&self.max[b], &self.min[b]);
            if mx.seconds[S_SOLVE] == 0.0 {
                continue;
            }
            let f = |idx: &[usize]| -> String {
                idx.iter()
                    .map(|&i| format!("{:.4} [{:.4}]", mx.seconds[i], mn.seconds[i]))
                    .collect::<Vec<_>>()
                    .join(" ")
            };
            let _ = writeln!(
                out,
                "  {:<12} | {} | {} | {:.4} = {:.4} + {:.4} | {:.4}, {:.4} ({} calls)",
                label,
                f(&[S_SOLVE_SETUP, S_SOLVE_PRECONDITION, S_SOLVE_APPEND]),
                f(&[
                    S_SOLVE_REAL_STATE,
                    S_M2_SOLVE_PREPARE,
                    S_SOLVE_ENTRY,
                    S_SOLVE_KEPT_UPDATE,
                    S_SOLVE_FINISH,
                    S_M2_SOLVE_UNPACK,
                    S_SOLVE_ENTRY_WAIT
                ]),
                mx.seconds[S_SOLVE_OTHER],
                mx.seconds[S_OTHER_SOLVER],
                mx.seconds[S_OTHER_ENGINE],
                mx.seconds[S_SOLVE_WAIT],
                mx.seconds[S_COLL + 4],
                mx.counts[C_COLL + 4]
            );
        }
        out
    }
    /// Bucket entries as (bucket label, [(entry label, max, min)]) for JSON
    /// writers: seconds first, counts after.
    pub fn entries(&self) -> Vec<(&'static str, Vec<Entry>)> {
        (0..BUCKETS)
            .map(|b| {
                let mut v: Vec<(&'static str, f64, f64)> = SEC_LABELS
                    .iter()
                    .enumerate()
                    .map(|(i, l)| (*l, self.max[b].seconds[i], self.min[b].seconds[i]))
                    .collect();
                v.extend(COUNT_LABELS.iter().enumerate().map(|(i, l)| {
                    (
                        *l,
                        self.max[b].counts[i] as f64,
                        self.min[b].counts[i] as f64,
                    )
                }));
                (BUCKET_LABELS[b], v)
            })
            .collect()
    }
}

/// Context of the solve in progress (stage, map kind and SCF iteration), set
/// by the session and the ladder; read by the snapshot and first-block
/// diagnostics for their labels. Keyed by the rank of the registered
/// communicator (0 on the threads-only path), so that a solve running on a
/// worker thread of the pool still finds the context its session set.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Context {
    pub stage: &'static str,
    pub kind: &'static str,
    pub iteration: usize,
}
impl Default for Context {
    fn default() -> Self {
        Self {
            stage: "fine",
            kind: "ordinary",
            iteration: 0,
        }
    }
}
static CONTEXTS: Mutex<Vec<Context>> = Mutex::new(Vec::new());
/// Rank key of the calling thread: the registered communicator's rank, else 0.
pub(crate) fn rank_key() -> usize {
    crate::parallel::SharedComm::active().map_or(0, |c| c.rank())
}
pub fn set_context(c: Context) {
    let r = rank_key();
    let mut v = CONTEXTS.lock().expect("context");
    if v.len() <= r {
        v.resize(r + 1, Context::default());
    }
    v[r] = c;
}
pub fn context() -> Context {
    let r = rank_key();
    CONTEXTS
        .lock()
        .expect("context")
        .get(r)
        .copied()
        .unwrap_or_default()
}
/// Set the stage of the context until the guard drops.
pub struct ContextGuard(Context);
impl Drop for ContextGuard {
    fn drop(&mut self) {
        set_context(self.0);
    }
}
#[must_use]
pub fn context_stage(stage: &'static str) -> ContextGuard {
    let previous = context();
    set_context(Context { stage, ..previous });
    ContextGuard(previous)
}

#[cfg(test)]
mod closure_tests {
    use super::*;

    fn rank(rank: usize, init: f64, fine: f64, wall: f64) -> RankProfile {
        let mut buckets = [Accum::default(); BUCKETS];
        buckets[Bucket::Init as usize].seconds[S_WALL] = init;
        buckets[Bucket::Fine as usize].seconds[S_WALL] = fine;
        RankProfile {
            rank,
            buckets,
            wall,
        }
    }

    /// Rank imbalance with the bucket maxima on different ranks: the legacy
    /// sum of maxima exceeds every rank's wall (an impossible row), the
    /// per-rank closure does not.
    #[test]
    fn per_rank_closure_differs_from_sum_of_maxima() {
        let r = Report::from_ranks(&[rank(0, 1.0, 3.0, 4.5), rank(1, 2.0, 1.5, 4.0)]);
        assert_eq!(r.covered_wall(), 5.0);
        let c = r.closure().unwrap();
        assert_eq!(
            (c.covered_max, c.covered_min, c.critical_rank),
            (4.0, 3.5, 0)
        );
        assert_eq!(c.imbalance, 0.125);
        assert_eq!((c.wall_max, c.wall_min), (Some(4.5), Some(4.0)));
        assert_eq!((c.residual_max, c.residual_min), (Some(0.5), Some(0.5)));
        assert!(r.covered_wall() > c.wall_max.unwrap());
        assert!(r.per_rank.iter().all(|p| p.residual.unwrap() >= 0.0));
        // No wall reported: no residual.
        let r = Report::from_ranks(&[rank(0, 1.0, 1.0, 0.0)]);
        assert_eq!(r.closure().unwrap().residual_max, None);
    }

    /// The wall travels with the profile through the real-MPI gather path
    /// (simulated ranks).
    #[test]
    fn gather_carries_the_rank_wall() {
        let out = dft_mpi::simulate(3, |comm| {
            let mine = rank(
                comm.rank(),
                0.5,
                comm.rank() as f64,
                1.0 + comm.rank() as f64,
            );
            gather(comm, &mine).unwrap()
        });
        for all in out {
            assert_eq!(all.len(), 3);
            for (r, p) in all.iter().enumerate() {
                assert_eq!(p.rank, r);
                assert_eq!(p.wall, 1.0 + r as f64);
                assert_eq!(p.buckets[Bucket::Fine as usize].seconds[S_WALL], r as f64);
            }
            let c = Report::from_ranks(&all).closure().unwrap();
            assert_eq!(c.critical_rank, 2);
            assert_eq!(c.residual_min, Some(0.5));
        }
    }
}

// LANE-FFT wire-format regression: real MPI uses the same final allgatherv.
#[cfg(test)]
mod fft_lane_tests {
    use super::*;
    #[test]
    fn lane_fft_gather_keeps_rank_buckets_and_u64_bits() {
        for enabled in [false, true] {
            dft_mpi::simulate(3, |comm| {
                let cfg = dft_numerics::switches::RunConfig::from_pairs(&[
                    ("fft_profile", if enabled { "on" } else { "off" }),
                    ("phase_profile", "on"),
                ])
                .unwrap();
                let _scope = cfg.enter();
                let mut mine = RankProfile {
                    rank: comm.rank(),
                    buckets: [Accum::default(); BUCKETS],
                    wall: 0.,
                };
                mine.buckets[2].counts[C_MAPS] = 11 + comm.rank() as u64;
                if enabled {
                    mine.buckets[2].fft.0[8] = [1, 2, 3, (1_u64 << 60) + comm.rank() as u64, 5];
                }
                let all = gather(comm, &mine).unwrap();
                assert_eq!(all[comm.rank()], mine);
                let report = Report::from_ranks(&all);
                assert_eq!(report.fft_ranks.len(), if enabled { 3 } else { 0 });
                for r in &report.fft_ranks {
                    assert_eq!(r.maps[2], 11 + r.rank as u64);
                    assert_eq!(r.buckets[2].0[8][3], (1_u64 << 60) + r.rank as u64);
                }
            });
        }
    }
}

/// Per-attempt verification lane instrumentation (Phase 0).
#[path = "verification_profile.rs"]
pub mod verification;
