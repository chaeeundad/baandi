//! Distributed full-grid stages of the M2 path (ticket M2-GRID, 2026-11-18;
//! `BAANDI_M2_GRID=on|off`, default off, rounding class).
//!
//! Without the switch every rank of the M2 path repeats the same full-grid work
//! per SCF map: the Hartree FFT pair of the input density, the effective
//! potential sum, the Hartree FFT of the output density in the energy audit and
//! the audit's grid sums. None of it shrinks with the rank count. With the switch
//! each rank works on its own real-space planes and reciprocal sticks of a
//! slab-decomposed full-grid FFT (the M2a `DistributedPairFft` with every grid
//! point as an index: sticks along the first FFT axis, contiguous planes of the
//! first axis per rank, two transposes per transform):
//!
//! * potential: forward transform of the (replicated) input density on the
//!   owned planes, 4 pi / G^2 on the owned sticks, inverse transform, LDA XC on
//!   the owned plane points, effective potential ionic + (V_H + V_xc) on the
//!   planes; ONE allgatherv brings the effective potential (all planes) and the
//!   partial sums (E_H per stick, XC energy and integral per plane, error flags)
//!   to every rank. V_H and V_xc stay on the planes for the audit.
//! * energy audit: forward transform of the output density on the planes,
//!   E_H per stick; charge, double counting, the two Bregman linear terms, the
//!   local energy (real-space form w sum n V_ion, Parseval-equal to the
//!   reciprocal form) and the output XC energy per plane; ONE allgatherv.
//! * stop rule (`BAANDI_STOP_METRIC=kept`): the owned stick coefficients of the
//!   input and output densities are kept; the Coulomb estimate is a per-stick
//!   partial sum plus one small allgatherv.
//!
//! Determinism and numeric class. The transforms move data between ranks but
//! never change the arithmetic of a plane (2D lines) or of a stick (1D line), so
//! every transformed value is the same for every rank count. Every reduction is
//! a canonical sum: partial sums per plane (row-major inside the plane) or per
//! stick (ascending first-axis index), gathered and added in ascending plane or
//! stick order. The results therefore do not depend on the rank count (P >= 2),
//! and a fixed rank count is bit-reproducible. Against the replicated path
//! (`FftPlan3d`, row-major full-grid sums) the potential and the energies differ
//! by rounding only: rounding class, hence the switch.
//!
//! Every branch condition is replicated data (switch, functional, communicator,
//! grid) and every rank enters the same collectives. An invalid density on one
//! rank does not leave the collectives early: it travels as a flag in the
//! gathered sums and every rank returns the same error afterwards.
use crate::ScfError;
use dft_mpi::Collective;
use std::sync::atomic::{AtomicU64, Ordering};

thread_local! {
    static OVERRIDE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}
/// Override `BAANDI_M2_GRID` for problems used on the calling thread (simulated
/// ranks in tests; every rank thread must set the same value).
pub fn set_override(on: Option<bool>) {
    OVERRIDE.with(|o| o.set(on));
}
/// True when the distributed grid stages are requested on this thread.
pub fn requested() -> bool {
    OVERRIDE
        .with(|o| o.get())
        .unwrap_or_else(|| dft_numerics::switches::on(dft_numerics::switches::Id::M2_GRID))
}

static COUNTS: [AtomicU64; 5] = [const { AtomicU64::new(0) }; 5];
#[cfg(feature = "strict-oracle")]
const C_POTENTIAL: usize = 0;
#[cfg(feature = "strict-oracle")]
const C_AUDIT: usize = 1;
#[cfg(feature = "strict-oracle")]
const C_STOP: usize = 2;
#[cfg(feature = "strict-oracle")]
const C_PLANS: usize = 3;
#[cfg(feature = "strict-oracle")]
const C_DECLINED: usize = 4;
/// (distributed potential maps, distributed audits, distributed stop-rule
/// estimates, plans built, maps that declined: not eligible) since process
/// start, summed over the threads of the process (simulated ranks add up).
pub fn counters() -> [u64; 5] {
    std::array::from_fn(|i| COUNTS[i].load(Ordering::Relaxed))
}
#[cfg(feature = "strict-oracle")]
fn bump(i: usize) {
    COUNTS[i].fetch_add(1, Ordering::Relaxed);
}
#[cfg(feature = "strict-oracle")]
pub(crate) fn count_declined() {
    bump(C_DECLINED);
}

// ---------------------------------------------------------------------------
// M2-OWNED (`BAANDI_M2_LOCAL_FIELDS`, `BAANDI_M2_LOCAL_DENSITY`, default off):
// switches, counters and the plane redistribution of the owned fields.

thread_local! {
    static OWNED_OVERRIDE: std::cell::Cell<Option<(bool, bool)>> =
        const { std::cell::Cell::new(None) };
}
/// Override (local fields, local density) for the calling thread (simulated
/// ranks in tests; every rank thread must set the same value).
pub fn set_owned_override(v: Option<(bool, bool)>) {
    OWNED_OVERRIDE.with(|o| o.set(v));
}
/// `BAANDI_M2_LOCAL_FIELDS=on` on this thread.
pub fn local_fields_requested() -> bool {
    OWNED_OVERRIDE.with(|o| o.get()).map_or_else(
        || dft_numerics::switches::on(dft_numerics::switches::Id::M2_LOCAL_FIELDS),
        |v| v.0,
    )
}
/// `BAANDI_M2_LOCAL_DENSITY=on` on this thread.
pub fn local_density_requested() -> bool {
    OWNED_OVERRIDE.with(|o| o.get()).map_or_else(
        || dft_numerics::switches::on(dft_numerics::switches::Id::M2_LOCAL_DENSITY),
        |v| v.1,
    )
}

#[cfg(feature = "strict-oracle")]
pub(crate) const O_FIELDS_OWNED_MAPS: usize = 0;
#[cfg(feature = "strict-oracle")]
pub(crate) const O_FIELDS_DECLINED_MAPS: usize = 1;
#[cfg(feature = "strict-oracle")]
pub(crate) const O_POTENTIAL_GATHERS: usize = 2;
#[cfg(feature = "strict-oracle")]
pub(crate) const O_POTENTIAL_GATHER_VALUES: usize = 3;
#[cfg(feature = "strict-oracle")]
pub(crate) const O_CONSUMER_GATHERS: usize = 4;
#[cfg(feature = "strict-oracle")]
pub(crate) const O_CONSUMER_GATHER_VALUES: usize = 5;
pub(crate) const O_REDISTRIBUTIONS: usize = 6;
pub(crate) const O_REDISTRIBUTION_VALUES: usize = 7;
#[cfg(feature = "strict-oracle")]
pub(crate) const O_FULL_COPIES: usize = 8;
#[cfg(feature = "strict-oracle")]
pub(crate) const O_FULL_COPY_VALUES: usize = 9;
#[cfg(feature = "strict-oracle")]
pub(crate) const O_DENSITY_OWNED_MAPS: usize = 10;
#[cfg(feature = "strict-oracle")]
pub(crate) const O_DENSITY_DECLINED: usize = 11;
#[cfg(feature = "strict-oracle")]
pub(crate) const O_DENSITY_GATHERS: usize = 12;
#[cfg(feature = "strict-oracle")]
pub(crate) const O_DENSITY_GATHER_VALUES: usize = 13;
#[cfg(feature = "strict-oracle")]
pub(crate) const O_MIXER_GATHERS: usize = 14;
#[cfg(feature = "strict-oracle")]
pub(crate) const O_MIXER_GATHER_VALUES: usize = 15;
#[cfg(feature = "strict-oracle")]
pub(crate) const O_MATERIALIZATIONS: usize = 16;
#[cfg(feature = "strict-oracle")]
pub(crate) const O_MATERIALIZATION_VALUES: usize = 17;
#[cfg(feature = "strict-oracle")]
pub(crate) const O_EPOCH_CHECKS: usize = 18;
#[cfg(feature = "strict-oracle")]
pub(crate) const O_OWNED_REDUCTIONS: usize = 19;
#[cfg(feature = "strict-oracle")]
pub(crate) const O_PLANE_VALUES: usize = 20;
const OWNED_NAMES: [&str; 21] = [
    "fields_owned_maps",
    "fields_declined_maps",
    "potential_gathers",
    "potential_gather_values",
    "consumer_potential_gathers",
    "consumer_potential_gather_values",
    "redistributions",
    "redistribution_values",
    "full_copies",
    "full_copy_values",
    "density_owned_maps",
    "density_declined_maps",
    "density_gathers",
    "density_gather_values",
    "mixer_gathers",
    "mixer_gather_values",
    "density_materializations",
    "density_materialization_values",
    "epoch_checks",
    "owned_reductions",
    "plane_values_per_rank_max",
];
static OWNED: [AtomicU64; 21] = [const { AtomicU64::new(0) }; 21];
#[cfg(feature = "strict-oracle")]
pub(crate) const T_POTENTIAL_GATHER: usize = 0;
#[cfg(feature = "strict-oracle")]
pub(crate) const T_CONSUMER_GATHER: usize = 1;
pub(crate) const T_REDISTRIBUTE: usize = 2;
#[cfg(feature = "strict-oracle")]
pub(crate) const T_DENSITY_GATHER: usize = 3;
#[cfg(feature = "strict-oracle")]
pub(crate) const T_MIXER_GATHER: usize = 4;
#[cfg(feature = "strict-oracle")]
pub(crate) const T_MATERIALIZE: usize = 5;
#[cfg(feature = "strict-oracle")]
pub(crate) const T_OWNED_REDUCE: usize = 6;
const OWNED_TIME_NAMES: [&str; 7] = [
    "potential_gather",
    "consumer_potential_gather",
    "redistribute",
    "density_gather",
    "mixer_gather",
    "density_materialize",
    "owned_reduce",
];
static OWNED_NANOS: [AtomicU64; 7] = [const { AtomicU64::new(0) }; 7];
/// Reasons a map (fields) or a session (density) that asked for the owned
/// mode took the gathered path. Indices 7.. are the declared or automatic
/// consumer declines (`OwnedDecline`) and the model reasons spin/functional.
pub(crate) const DECLINE_REASONS: [&str; 11] = [
    "grid_path",
    "standard",
    "kpoints",
    "cold_map",
    "operator",
    "snapshot",
    "finite_temperature",
    "checkpoint",
    "runtime_control",
    "spin",
    "functional",
];
#[cfg(feature = "strict-oracle")]
pub(crate) const R_GRID_PATH: usize = 0;
#[cfg(feature = "strict-oracle")]
pub(crate) const R_SPIN: usize = 9;
#[cfg(feature = "strict-oracle")]
pub(crate) const R_FUNCTIONAL: usize = 10;
static DECLINED: [AtomicU64; 11] = [const { AtomicU64::new(0) }; 11];
static DENSITY_DECLINED: [AtomicU64; 11] = [const { AtomicU64::new(0) }; 11];

/// A consumer that needs the full density, potential or mixer state at an
/// ACTIVE SCF boundary, declared before the session is initialized
/// (`ScfProblem::decline_owned`; the session declares `RuntimeControl` itself
/// for agent-controlled runs and `Checkpoint` for an exact restore). The run
/// then takes the gathered path for both owned switches (and, for
/// `Checkpoint`, the replicated Pulay mixer, whose state checkpoint v1 can
/// carry), counted as `owned.fields_declined.<reason>` /
/// `owned.density_declined.<reason>`. Consumers of the DELIVERED result
/// (forces, stress, volumetric output, density files) need no decline: the
/// delivered density and orbitals are gathered at delivery.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OwnedDecline {
    /// Exact or warm checkpoint of an active boundary, or an exact restore.
    Checkpoint,
    /// Runtime policy control (mixer replacement, session density at a
    /// boundary): `ExecutionMode::AgentControlled`.
    RuntimeControl,
}
impl OwnedDecline {
    pub(crate) fn index(self) -> usize {
        match self {
            Self::Checkpoint => 7,
            Self::RuntimeControl => 8,
        }
    }
    pub fn name(self) -> &'static str {
        DECLINE_REASONS[self.index()]
    }
}
/// Count a decline of the owned density loop with its reason (once per
/// session initialization that asked for it).
#[cfg(feature = "strict-oracle")]
pub(crate) fn density_declined(reason: usize) {
    DENSITY_DECLINED[reason].fetch_add(1, Ordering::Relaxed);
    owned_add(O_DENSITY_DECLINED, 1);
}
/// A spin (collinear) SCF asked for the owned mode: it has its own driver
/// and always takes the gathered path (counted once per run, both switches).
#[cfg(feature = "strict-oracle")]
pub(crate) fn count_spin_decline() {
    if crate::m2::active().is_none() {
        return;
    }
    if local_fields_requested() {
        owned_declined(R_SPIN);
    }
    if local_density_requested() {
        density_declined(R_SPIN);
    }
}
pub(crate) fn owned_add(i: usize, n: u64) {
    OWNED[i].fetch_add(n, Ordering::Relaxed);
}
#[cfg(feature = "strict-oracle")]
pub(crate) fn owned_max(i: usize, n: u64) {
    OWNED[i].fetch_max(n, Ordering::Relaxed);
}
/// Charge the time since `start` (None with `BAANDI_STAGE_TIMERS=off`: no clock).
pub(crate) fn owned_time(i: usize, start: Option<std::time::Instant>) {
    if start.is_some() {
        OWNED_NANOS[i].fetch_add(
            dft_numerics::stage_clock::elapsed_ns(start),
            Ordering::Relaxed,
        );
    }
}
#[cfg(feature = "strict-oracle")]
pub(crate) fn owned_declined(reason: usize) {
    DECLINED[reason].fetch_add(1, Ordering::Relaxed);
    owned_add(O_FIELDS_DECLINED_MAPS, 1);
}
/// M2-OWNED counters (`owned.<name>`, process-wide, simulated ranks add up)
/// and exposed seconds (`owned.<stage>`, 0 with `BAANDI_STAGE_TIMERS=off`).
/// Values are f64 doubles (8 bytes each). results.json `parallel.m2`.
pub fn owned_stats() -> crate::parallel::M2Stats {
    let mut counts: Vec<(String, usize)> = OWNED_NAMES
        .iter()
        .zip(&OWNED)
        .map(|(n, c)| (format!("owned.{n}"), c.load(Ordering::Relaxed) as usize))
        .collect();
    counts.extend(DECLINE_REASONS.iter().zip(&DECLINED).map(|(n, c)| {
        (
            format!("owned.fields_declined.{n}"),
            c.load(Ordering::Relaxed) as usize,
        )
    }));
    counts.extend(DECLINE_REASONS.iter().zip(&DENSITY_DECLINED).map(|(n, c)| {
        (
            format!("owned.density_declined.{n}"),
            c.load(Ordering::Relaxed) as usize,
        )
    }));
    let stages = OWNED_TIME_NAMES
        .iter()
        .zip(&OWNED_NANOS)
        .map(|(n, c)| {
            (
                format!("owned.{n}"),
                c.load(Ordering::Relaxed) as f64 * 1e-9,
            )
        })
        .collect();
    (counts, stages)
}

/// Epoch-tagged plane-local field of the slab layout (M2-OWNED): the values
/// of this rank's planes and the epoch of the map that produced them. A
/// consumer compares the epoch it expects with the tag, so a stale field is
/// rejected instead of used.
#[derive(Clone, Debug)]
#[cfg(feature = "strict-oracle")]
pub(crate) struct Owned {
    pub(crate) epoch: u64,
    pub(crate) values: Vec<f64>,
}

/// Owned-field state of one problem (M2-OWNED). `potential` is the effective
/// potential of the latest owned map on the distributed operator's planes;
/// `complex_stale` says that the complex Hamiltonians still hold an older full
/// potential (a consumer of the full potential gathers it first).
#[derive(Default)]
#[cfg(feature = "strict-oracle")]
pub(crate) struct OwnedFields {
    epoch: u64,
    pub(crate) potential: Option<Owned>,
    pub(crate) complex_stale: bool,
    pub(crate) map: Option<PlaneMap>,
    pub(crate) scratch: Vec<f64>,
    /// `BAANDI_M2_LOCAL_DENSITY` loop active (set by the session): maps keep
    /// their output density on the slab planes (`output`), and a map given an
    /// empty input array reads its input from `input`.
    pub(crate) density_mode: bool,
    pub(crate) input: Option<Owned>,
    pub(crate) output: Option<Owned>,
    /// Declared consumer of active-boundary state (`OwnedDecline`): both
    /// owned modes decline for this problem.
    pub(crate) declined: Option<OwnedDecline>,
}
#[cfg(feature = "strict-oracle")]
impl OwnedFields {
    /// A new epoch for a field produced now (monotonic per problem).
    pub(crate) fn next_epoch(&mut self) -> u64 {
        self.epoch += 1;
        self.epoch
    }
}

/// Plane start of rank r for n0 planes over `ranks` ranks (the rule of
/// `DistributedPairFft`: r * (n0 / P) + min(r, n0 % P)).
fn plane_starts(n0: usize, ranks: usize) -> Vec<usize> {
    (0..=ranks)
        .map(|r| r * (n0 / ranks) + r.min(n0 % ranks))
        .collect()
}

/// Redistribution of plane-local fields between two slab layouts of one grid
/// (layout A: the slab plan, planes along `first_a`; layout B: the distributed
/// operator, planes along `first_b`). Both use the same plane rule, so equal
/// axes mean equal layouts and the move is a local copy; otherwise one
/// alltoall of equal padded blocks (every rank sends its A planes cut by the
/// B plane ranges). A pure permutation: the values keep their bits.
pub(crate) struct PlaneMap {
    dims: [usize; 3],
    rank: usize,
    first_a: usize,
    first_b: usize,
    starts_a: Vec<usize>,
    starts_b: Vec<usize>,
    block: usize,
    send: Vec<f64>,
    recv: Vec<f64>,
}

impl PlaneMap {
    pub(crate) fn new(
        dims: [usize; 3],
        rank: usize,
        ranks: usize,
        first_a: usize,
        first_b: usize,
    ) -> Self {
        let starts_a = plane_starts(dims[first_a], ranks);
        let starts_b = plane_starts(dims[first_b], ranks);
        let max = |s: &[usize]| s.windows(2).map(|w| w[1] - w[0]).max().unwrap_or(0);
        let block = max(&starts_a) * max(&starts_b) * dims[2];
        Self {
            dims,
            rank,
            first_a,
            first_b,
            starts_a,
            starts_b,
            block,
            send: Vec::new(),
            recv: Vec::new(),
        }
    }
    #[cfg(feature = "strict-oracle")]
    pub(crate) fn identity(&self) -> bool {
        self.first_a == self.first_b
    }
    #[cfg(feature = "strict-oracle")]
    pub(crate) fn matches(&self, dims: [usize; 3], first_a: usize, first_b: usize) -> bool {
        self.dims == dims && self.first_a == first_a && self.first_b == first_b
    }
    /// Layout A -> layout B. Collective unless `identity()`.
    pub(crate) fn a_to_b<K: Collective>(
        &mut self,
        comm: &K,
        src: &[f64],
        dst: &mut Vec<f64>,
    ) -> Result<(), ScfError> {
        let (sa, sb) = (
            std::mem::take(&mut self.starts_a),
            std::mem::take(&mut self.starts_b),
        );
        let r = self.moved(comm, (self.first_a, &sa), (self.first_b, &sb), src, dst);
        (self.starts_a, self.starts_b) = (sa, sb);
        r
    }
    /// Layout B -> layout A. Collective unless `identity()`.
    pub(crate) fn b_to_a<K: Collective>(
        &mut self,
        comm: &K,
        src: &[f64],
        dst: &mut Vec<f64>,
    ) -> Result<(), ScfError> {
        let (sa, sb) = (
            std::mem::take(&mut self.starts_a),
            std::mem::take(&mut self.starts_b),
        );
        let r = self.moved(comm, (self.first_b, &sb), (self.first_a, &sa), src, dst);
        (self.starts_a, self.starts_b) = (sa, sb);
        r
    }
    fn moved<K: Collective>(
        &mut self,
        comm: &K,
        (fs, ss): (usize, &[usize]),
        (ft, st): (usize, &[usize]),
        src: &[f64],
        dst: &mut Vec<f64>,
    ) -> Result<(), ScfError> {
        let nz = self.dims[2];
        let rank = self.rank;
        let ranks = ss.len() - 1;
        if fs == ft {
            dst.clear();
            dst.extend_from_slice(src);
            return Ok(());
        }
        let clock = dft_numerics::stage_clock::start();
        // Rows of the source layout run along the target's plane axis and
        // vice versa (the two axes are 0 and 1).
        let (rows_s, rows_t) = (self.dims[ft], self.dims[fs]);
        let (s0, s1) = (ss[rank], ss[rank + 1]);
        if src.len() != (s1 - s0) * rows_s * nz {
            return Err(ScfError::InvalidInput("M2 owned redistribution: length"));
        }
        self.send.clear();
        self.send.resize(self.block * ranks, 0.0);
        for q in 0..ranks {
            let mut k = q * self.block;
            for x in s0..s1 {
                for y in st[q]..st[q + 1] {
                    let at = ((x - s0) * rows_s + y) * nz;
                    self.send[k..k + nz].copy_from_slice(&src[at..at + nz]);
                    k += nz;
                }
            }
        }
        self.recv.resize(self.block * ranks, 0.0);
        if self.block > 0 {
            comm.alltoall_f64(&self.send, self.block, &mut self.recv)
                .map_err(|e| ScfError::Numerical(format!("M2 owned redistribution: {e}")))?;
        }
        let (t0, t1) = (st[rank], st[rank + 1]);
        dst.clear();
        dst.resize((t1 - t0) * rows_t * nz, 0.0);
        for r in 0..ranks {
            let mut k = r * self.block;
            for x in ss[r]..ss[r + 1] {
                for y in t0..t1 {
                    let at = ((y - t0) * rows_t + x) * nz;
                    dst[at..at + nz].copy_from_slice(&self.recv[k..k + nz]);
                    k += nz;
                }
            }
        }
        owned_add(O_REDISTRIBUTIONS, 1);
        owned_add(O_REDISTRIBUTION_VALUES, (self.block * ranks) as u64);
        owned_time(T_REDISTRIBUTE, clock);
        Ok(())
    }
}

#[cfg(feature = "strict-oracle")]
#[path = "oracle/grid_dist.rs"]
mod legacy;
#[cfg(feature = "strict-oracle")]
pub(crate) use legacy::*;
