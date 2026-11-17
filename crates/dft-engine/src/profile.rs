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
}
pub const BUCKETS: usize = 7;
pub const BUCKET_LABELS: [&str; BUCKETS] = [
    "init",
    "coarse",
    "fine",
    "fine_redo",
    "verify_pass",
    "verify_fail",
    "materialize",
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

// Seconds slots.
pub const SEC_LABELS: [&str; 32] = [
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
pub const NSEC: usize = 32;

// Count slots.
pub const COUNT_LABELS: [&str; 20] = [
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
pub const NCOUNT: usize = 20;

/// Accumulated seconds and counts of one bucket.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Accum {
    pub seconds: [f64; NSEC],
    pub counts: [u64; NCOUNT],
}
impl Default for Accum {
    fn default() -> Self {
        Self {
            seconds: [0.0; NSEC],
            counts: [0; NCOUNT],
        }
    }
}
impl Accum {
    fn add(&mut self, other: &Accum) {
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
}

struct State {
    top: Top,
    buckets: [Accum; BUCKETS],
    /// The map in progress (committed by the session or by the next map).
    pending: Option<Accum>,
    /// Pending map's atomic-start seconds (excluded from its potential).
    pending_atomic: f64,
    mark: Snapshot,
}
impl State {
    fn new() -> Self {
        Self {
            top: Top::Fine,
            buckets: [Accum::default(); BUCKETS],
            pending: None,
            pending_atomic: 0.0,
            mark: local::snapshot(),
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
        let lazy = self.buckets[bucket as usize].add_delta(&d);
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
    let per = 1 + BUCKETS * (NSEC + NCOUNT);
    let mut local = Vec::with_capacity(per);
    local.push(mine.rank as f64);
    for b in &mine.buckets {
        local.extend_from_slice(&b.seconds);
        local.extend(b.counts.iter().map(|&c| c as f64));
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
            };
            let mut at = 1;
            for b in &mut p.buckets {
                b.seconds.copy_from_slice(&c[at..at + NSEC]);
                at += NSEC;
                for (k, v) in b.counts.iter_mut().zip(&c[at..at + NCOUNT]) {
                    *k = *v as u64;
                }
                at += NCOUNT;
            }
            p
        })
        .collect();
    out.sort_by_key(|p| p.rank);
    Ok(out)
}

/// One report entry: (label, maximum over ranks, minimum over ranks).
pub type Entry = (&'static str, f64, f64);

/// Maximum and minimum over ranks of every entry of every bucket.
#[derive(Clone, Debug, PartialEq)]
pub struct Report {
    pub ranks: usize,
    pub max: [Accum; BUCKETS],
    pub min: [Accum; BUCKETS],
}
impl Report {
    pub fn from_ranks(profiles: &[RankProfile]) -> Self {
        let mut max = [Accum::default(); BUCKETS];
        let mut min = [Accum::default(); BUCKETS];
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
            // `solve.other`: the solve wall not covered by the exclusive stages.
            for a in [&mut max[b], &mut min[b]] {
                let parts: f64 = (S_SOLVE_START..=S_SOLVE_RETURN)
                    .map(|i| a.seconds[i])
                    .sum::<f64>()
                    + a.seconds[S_SOLVE_OTHER];
                a.seconds[S_SOLVE_OTHER] = (a.seconds[S_SOLVE] - parts).max(0.0);
            }
        }
        Self {
            ranks: profiles.len(),
            max,
            min,
        }
    }
    /// Sum of the bucket walls (maximum over ranks per bucket).
    pub fn covered_wall(&self) -> f64 {
        self.max.iter().map(|a| a.seconds[S_WALL]).sum()
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
                    "  {:<12} atomic start {:.4} s",
                    "", mx.seconds[S_ATOMIC]
                );
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
