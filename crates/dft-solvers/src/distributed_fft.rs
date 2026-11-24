//! G-vector (stick) distributed pruned 3D FFT for packed real pairs (MPI stage M2a).
//!
//! Placement: this module lives in dft-solvers because that crate already depends
//! on both `dft_mpi::Collective` and dft-numerics, and dft-numerics must stay free
//! of the MPI layer. Pairing of two real columns into one complex vector is done
//! upstream (`PairMap`); this plan transforms one complex coefficient vector.
//!
//! Conventions follow `PairFusedPlan`: the first transformed axis is the one with
//! fewer occupied sticks (ties: x), z is last. Sticks are columns along the first
//! axis, keyed by (z, other axis).
//!
//! Partition rule (deterministic in shape, indices, ranks): sticks sorted by
//! length descending (ties: lower (z, other) slot) are assigned one by one to the
//! rank with the smallest accumulated length (ties: lower rank); with
//! `BAANDI_M2_STICKS=mirror` (MPI-SCALING-GAP) a stick and its Gamma mirror stick
//! are assigned together as one unit (no pair representative is split, so the
//! mirror exchanges of the real Gamma operator vanish). Each rank owns
//! the sticks assigned to it; its local coefficients are the entries of the global
//! coefficient vector that fall on those sticks, in ascending global position
//! (`local_positions`). Real-space planes are contiguous ranges of the first axis
//! (sizes differ by at most one; ranks beyond n0 own none).
//!
//! Transpose: one `alltoall_f64` per direction with equal blocks of
//! 2*max_sticks*max_planes f64 per transformed vector; shorter (rank, plane-range)
//! pairs are zero-padded and the padding is never read. Batches (M2d): the `_batch`
//! entry points transform np vectors with ONE alltoall per direction (block
//! 2*max_sticks*max_planes*np, destination-major, vector-minor); every vector's
//! arithmetic is exactly that of the single-vector path, so a batch is bit-identical
//! to np single calls. Transpose buffers grow to the widest batch seen and are then
//! reused (no allocation once the widest batch has occurred). Per-plane work reuses the serial order of
//! operations (clear, fill, batched 1D FFT along the second axis, 1D FFT along
//! z), so each plane's slab is expected to be bit-identical to the serial plan's
//! for any rank count (tests assert bit identity for one rank only). The 1/N
//! factor is applied at the end of the reciprocal direction.
//! 1D kernels (MPI entry cost, `BAANDI_M2_FFT`): rustfft by default; `mkl` selects
//! MKL DFTI batch lanes (one descriptor per lane buffer, as `PairFusedPlan<MklLanes>`)
//! on MKL builds, a rounding change; pure builds always use rustfft.
//! `BAANDI_FFT_FUSED_ROWS=on`: ordinary local-potential applications stream
//! the last-axis inverse/multiply/forward through two reusable rows. Immutable
//! inverse input retains zero holes; the standalone plane and kept-plane APIs
//! retain the slab route. RustFFT only (vendor lanes fall back), conservatively
//! rounding class; off retains the original transforms and data movement.
//! No allocation after batch warmup; no unsafe; no logging.
//! `BAANDI_FFT_A2A_BALANCE=on` caps sticks per rank (mirror units stay intact)
//! only when it reduces the equal alltoall block. Transform arithmetic stays
//! unchanged; downstream reductions over redistributed rows are rounding class.
use dft_mpi::{Collective, MpiError};
use dft_numerics::FftError;
use dft_numerics::fused_pair_plan::BackendLanes;
use num_complex::Complex64 as C;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::{Arc, Mutex};
use thiserror::Error;

#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum DistFftError {
    #[error("FFT: {0}")]
    Fft(#[from] FftError),
    #[error("MPI: {0}")]
    Mpi(#[from] MpiError),
}

/// A constructor-only cache owned by one run/rank. rustfft caches by length and
/// direction; MKL descriptors additionally depend on batch count. No global
/// cache survives a run. Work arrays remain private to each FFT workspace.
#[derive(Clone, Default)]
pub struct PlanCache(Arc<Mutex<CachedPlans>>);
struct CachedPlans {
    share_kernels: bool,
    cache_partition: bool,
    partition: Option<CachedPartition>,
    partition_hits: usize,
    rust: rustfft::FftPlanner<f64>,
    lengths: BTreeSet<usize>,
    vendor: BTreeMap<(usize, usize), Arc<Mutex<BackendLanes>>>,
}
impl Default for CachedPlans {
    fn default() -> Self {
        Self {
            share_kernels: true,
            cache_partition: false,
            partition: None,
            partition_hits: 0,
            rust: rustfft::FftPlanner::new(),
            lengths: BTreeSet::new(),
            vendor: BTreeMap::new(),
        }
    }
}
thread_local! {
    static PLAN_CACHE: std::cell::RefCell<Option<PlanCache>> = const { std::cell::RefCell::new(None) };
}
/// Constructor scope only; intentionally cannot move to a worker thread.
pub struct PlanScope {
    previous: Option<PlanCache>,
    _thread: std::marker::PhantomData<std::rc::Rc<()>>,
}
impl Drop for PlanScope {
    fn drop(&mut self) {
        PLAN_CACHE.with(|c| *c.borrow_mut() = self.previous.take());
    }
}
impl PlanCache {
    /// Keep at most the most recent exact partition during constructors.
    /// Kernel sharing and layout reuse are independent experiments.
    pub fn with_partition_cache(share_kernels: bool, cache_partition: bool) -> Self {
        Self(Arc::new(Mutex::new(CachedPlans {
            share_kernels,
            cache_partition,
            ..Default::default()
        })))
    }
    pub fn partition_hits(&self) -> usize {
        self.0.lock().unwrap().partition_hits
    }
    pub fn enter(&self) -> PlanScope {
        PlanScope {
            previous: PLAN_CACHE.with(|c| c.replace(Some(self.clone()))),
            _thread: std::marker::PhantomData,
        }
    }
    /// Distinct rustfft (length, direction) requests and cached batch dispatchers
    /// (explicit MKL or calibrated lanes, including the pure-build fallback).
    /// Counts are construction evidence, not hot-loop instrumentation.
    pub fn counts(&self) -> (usize, usize) {
        let c = self.0.lock().unwrap();
        (2 * c.lengths.len(), c.vendor.len())
    }
    fn rust_pair(&self, n: usize) -> [Arc<dyn rustfft::Fft<f64>>; 2] {
        let mut c = self.0.lock().unwrap();
        c.lengths.insert(n);
        [c.rust.plan_fft_forward(n), c.rust.plan_fft_inverse(n)]
    }
}

enum Lanes {
    Private(BackendLanes),
    Shared(SharedLanes),
}
enum SharedLanes {
    Rust {
        plans: [[Arc<dyn rustfft::Fft<f64>>; 2]; 3],
        scratch: Vec<C>,
    },
    Vendor([Option<Arc<Mutex<BackendLanes>>>; 3]),
}
impl Clone for SharedLanes {
    fn clone(&self) -> Self {
        match self {
            Self::Rust { plans, scratch } => Self::Rust {
                plans: plans.clone(),
                scratch: vec![C::default(); scratch.len()],
            },
            Self::Vendor(plans) => Self::Vendor(plans.clone()),
        }
    }
}
impl Lanes {
    fn plan(n: [usize; 3], counts: [usize; 3], vendor: bool) -> Result<Self, FftError> {
        let cache = PLAN_CACHE.with(|c| c.borrow().clone());
        let Some(cache) = cache else {
            return BackendLanes::plan(n, counts, vendor).map(Self::Private);
        };
        if !cache.0.lock().unwrap().share_kernels {
            return BackendLanes::plan(n, counts, vendor).map(Self::Private);
        }
        // P3-PERMAP: the shared rustfft fast path used to bypass BackendLanes
        // entirely, including its run-scoped selector. Reuse the existing
        // per-(length,batch) shared dispatch entries for calibrated lanes too;
        // workspace sharing and every transform call remain unchanged.
        let calibrated =
            dft_numerics::switches::on(dft_numerics::switches::Id::FFT_KERNEL_BY_LENGTH);
        if calibrated || (vendor && dft_numerics::MKL_BACKEND) {
            let mut c = cache.0.lock().unwrap();
            let mut plans = [None, None, None];
            for a in 0..3 {
                if counts[a] != 0 {
                    let key = (n[a], counts[a]);
                    if let std::collections::btree_map::Entry::Vacant(e) = c.vendor.entry(key) {
                        // Only axis zero has a descriptor; zero-count axes do
                        // not create MKL plans. Use the existing FFI adapter.
                        e.insert(Arc::new(Mutex::new(BackendLanes::plan(
                            [n[a], 1, 1],
                            [counts[a], 0, 0],
                            true,
                        )?)));
                    }
                    plans[a] = Some(c.vendor[&key].clone());
                }
            }
            Ok(Self::Shared(SharedLanes::Vendor(plans)))
        } else {
            let plans = n.map(|d| cache.rust_pair(d));
            let scratch = plans
                .iter()
                .flatten()
                .map(|p| p.get_inplace_scratch_len())
                .max()
                .unwrap();
            Ok(Self::Shared(SharedLanes::Rust {
                plans,
                scratch: vec![C::default(); scratch],
            }))
        }
    }
    fn is_vendor(&self) -> bool {
        match self {
            Self::Private(p) => p.is_vendor(),
            Self::Shared(SharedLanes::Vendor(_)) => true,
            Self::Shared(SharedLanes::Rust { .. }) => false,
        }
    }
    fn run(&mut self, axis: usize, inverse: bool, data: &mut [C]) -> Result<(), FftError> {
        match self {
            Self::Private(p) => p.run(axis, inverse, data),
            Self::Shared(SharedLanes::Rust { plans, scratch }) => {
                plans[axis][usize::from(inverse)].process_with_scratch(data, scratch);
                Ok(())
            }
            Self::Shared(SharedLanes::Vendor(plans)) => match &plans[axis] {
                Some(p) => p.lock().unwrap().run(0, inverse, data),
                None if data.is_empty() => Ok(()),
                None => Err(FftError::IncompatibleArray),
            },
        }
    }
}

#[derive(Clone)]
struct Partition {
    first: usize,
    n: [usize; 3],
    sticks: Vec<(usize, usize)>,
    planes: Vec<usize>,
    slots: Vec<usize>,
    owner: Vec<usize>,
    local_index: Vec<usize>,
    counts: Vec<usize>,
}

struct CachedPartition {
    shape: [usize; 3],
    indices: Vec<[usize; 3]>,
    ranks: usize,
    mirror: bool,
    balance: bool,
    value: Partition,
}

thread_local! {
    static MIRROR_OVERRIDE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}
/// Force the stick rule of the calling thread (tests with simulated ranks; every
/// rank thread must set the same value); `None` restores `BAANDI_M2_STICKS`.
pub fn set_mirror_sticks_override(mirror: Option<bool>) {
    MIRROR_OVERRIDE.with(|o| o.set(mirror));
}
/// True when the mirror-paired stick partition is in force on this thread.
pub fn mirror_sticks() -> bool {
    MIRROR_OVERRIDE.with(|o| o.get()).unwrap_or_else(|| {
        dft_numerics::switches::choice(dft_numerics::switches::Id::M2_STICKS) == "mirror"
    })
}

/// Stick partition shared by `new` and `position_owners` (inputs already validated).
fn partition(shape: [usize; 3], indices: &[[usize; 3]], ranks: usize) -> Partition {
    let cache = PLAN_CACHE.with(|c| c.borrow().clone());
    let Some(cache) = cache else {
        return partition_uncached(shape, indices, ranks);
    };
    let mut cached = cache.0.lock().unwrap();
    if !cached.cache_partition {
        return partition_uncached(shape, indices, ranks);
    }
    let mirror = mirror_sticks();
    let balance = dft_numerics::switches::on(dft_numerics::switches::Id::FFT_A2A_BALANCE);
    if let Some(p) = &cached.partition
        && p.shape == shape
        && p.ranks == ranks
        && p.mirror == mirror
        && p.balance == balance
        && p.indices == indices
    {
        let value = p.value.clone();
        cached.partition_hits += 1;
        return value;
    }
    let value = partition_uncached(shape, indices, ranks);
    cached.partition = Some(CachedPartition {
        shape,
        indices: indices.to_vec(),
        ranks,
        mirror,
        balance,
        value: value.clone(),
    });
    value
}

fn partition_uncached(shape: [usize; 3], indices: &[[usize; 3]], ranks: usize) -> Partition {
    let supports: [_; 2] = std::array::from_fn(|a| {
        indices
            .iter()
            .map(|p| (p[2], p[1 - a]))
            .collect::<BTreeSet<_>>()
    });
    let first = usize::from(supports[1].len() < supports[0].len());
    let n = [shape[first], shape[1 - first], shape[2]];
    let sticks: Vec<(usize, usize)> = supports[first].iter().copied().collect();
    let planes: Vec<usize> = sticks
        .iter()
        .map(|p| p.0)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let mut len = vec![0usize; sticks.len()];
    let slots: Vec<usize> = indices
        .iter()
        .map(|p| {
            let s = sticks.binary_search(&(p[2], p[1 - first])).unwrap();
            len[s] += 1;
            s
        })
        .collect();
    // Balanced greedy assignment by stick length.
    let mut order: Vec<usize> = (0..sticks.len()).collect();
    order.sort_by(|&a, &b| len[b].cmp(&len[a]).then(a.cmp(&b)));
    let mut load = vec![0usize; ranks];
    let mut owner = vec![usize::MAX; sticks.len()];
    // `BAANDI_M2_STICKS=mirror` (MPI-SCALING-GAP): a stick (z, o) and its mirror
    // ((-z) mod n2, (-o) mod n1) are assigned to one rank as a unit, so the two
    // positions +G and -G of every Gamma pair representative live on one rank.
    // A stick without a mirror in the set (non-symmetric index sets) or its own
    // mirror is assigned alone, exactly as in the greedy rule.
    let mirror = mirror_sticks();
    let n_other = n[1];
    let n_z = n[2];
    for &s in &order {
        if owner[s] != usize::MAX {
            continue;
        }
        let r = (0..ranks).min_by_key(|&r| (load[r], r)).unwrap();
        owner[s] = r;
        load[r] += len[s];
        if mirror {
            let (z, o) = sticks[s];
            let key = ((n_z - z) % n_z, (n_other - o) % n_other);
            if let Ok(m) = sticks.binary_search(&key)
                && owner[m] == usize::MAX
            {
                owner[m] = r;
                load[r] += len[m];
            }
        }
    }
    // P3-A2A: equal-block alltoall pays for the largest *stick count*, not
    // coefficient load. Keep the historical assignment unless a count-capped
    // assignment strictly reduces that block. No changes to FFT arithmetic.
    if dft_numerics::switches::on(dft_numerics::switches::Id::FFT_A2A_BALANCE) {
        balance_transpose_sticks(&sticks, &len, n, ranks, mirror, &mut owner);
    }
    let mut local_index = vec![0usize; sticks.len()];
    let mut counts = vec![0usize; ranks];
    for s in 0..sticks.len() {
        local_index[s] = counts[owner[s]];
        counts[owner[s]] += 1;
    }
    Partition {
        first,
        n,
        sticks,
        planes,
        slots,
        owner,
        local_index,
        counts,
    }
}

/// Minimize the maximum stick count without splitting a mirror unit. Units
/// have size one or two: capacity c is feasible iff R*c >= S and
/// R*floor(c/2) >= number_of_pairs. Assign pairs before singletons, longest
/// coefficient work first within each class, to the least-loaded feasible rank.
/// This is a layout change (reductions over coefficient rows can round
/// differently), not compression or a change to the mathematical FFT support.
fn balance_transpose_sticks(
    sticks: &[(usize, usize)],
    lengths: &[usize],
    n: [usize; 3],
    ranks: usize,
    mirror: bool,
    owner: &mut Vec<usize>,
) {
    let mut units = Vec::new();
    for (s, &(z, y)) in sticks.iter().enumerate() {
        let partner = if mirror {
            sticks
                .binary_search(&((n[2] - z) % n[2], (n[1] - y) % n[1]))
                .unwrap_or(s)
        } else {
            s
        };
        if partner < s {
            continue;
        }
        let size = 1 + usize::from(partner != s);
        let work = lengths[s] + if partner != s { lengths[partner] } else { 0 };
        units.push((s, partner, size, work));
    }
    let pairs = units.iter().filter(|u| u.2 == 2).count();
    let mut cap = sticks.len().div_ceil(ranks);
    if ranks * (cap / 2) < pairs {
        cap += 1;
    }
    let mut old_counts = vec![0usize; ranks];
    for &r in owner.iter() {
        old_counts[r] += 1;
    }
    if cap >= old_counts.into_iter().max().unwrap_or(0) {
        return;
    }
    units.sort_by(|a, b| b.2.cmp(&a.2).then(b.3.cmp(&a.3)).then(a.0.cmp(&b.0)));
    let mut counts = vec![0usize; ranks];
    let mut load = vec![0usize; ranks];
    let mut candidate = vec![0; sticks.len()];
    for (s, partner, size, work) in units {
        let r = (0..ranks)
            .filter(|&r| counts[r] + size <= cap)
            .min_by_key(|&r| (load[r], counts[r], r))
            .expect("feasible one/two-stick capacity");
        candidate[s] = r;
        candidate[partner] = r;
        counts[r] += size;
        load[r] += work;
    }
    *owner = candidate;
}

/// Owning rank of every global coefficient position (same partition as `new`).
pub fn position_owners(
    shape: [usize; 3],
    indices: &[[usize; 3]],
    ranks: usize,
) -> Result<Vec<usize>, DistFftError> {
    if ranks == 0 || shape.contains(&0) {
        return Err(FftError::InvalidGrid.into());
    }
    if indices.iter().any(|p| (0..3).any(|d| p[d] >= shape[d])) {
        return Err(FftError::IncompatibleArray.into());
    }
    let part = partition(shape, indices, ranks);
    Ok(part.slots.iter().map(|&s| part.owner[s]).collect())
}

/// Complement of scatter destinations, coalesced once at construction. Clearing
/// spans avoids both redundant zero stores and a per-point hole-index stream.
#[derive(Clone)]
struct HoleClear {
    spans: Vec<std::ops::Range<usize>>,
    points: usize,
}
impl HoleClear {
    fn new(len: usize, occupied: impl Iterator<Item = usize>) -> Self {
        let mut mask = vec![false; len];
        for i in occupied {
            mask[i] = true;
        }
        let mut spans = Vec::new();
        let mut start = 0;
        while start < len {
            if mask[start] {
                start += 1;
                continue;
            }
            let mut end = start + 1;
            while end < len && !mask[end] {
                end += 1;
            }
            spans.push(start..end);
            start = end;
        }
        let points = spans.iter().map(|s| s.len()).sum();
        Self { spans, points }
    }
    fn clear(&self, data: &mut [C]) {
        for span in &self.spans {
            data[span.clone()].fill(C::default());
        }
    }
}

/// Last-axis streaming workspace. The spectral input is immutable during the
/// inverse FFT: inactive z entries stay zero for the lifetime of the plan.
/// Only one row is materialized, multiplied and transformed back at a time.
/// RustFFT's immutable interface may choose different arithmetic from its
/// in-place interface, so this experiment is conservatively rounding class.
struct FusedRows {
    spectral: Vec<C>,
    real: Vec<C>,
    scratch: Vec<C>,
    inverse: std::sync::Arc<dyn rustfft::Fft<f64>>,
    forward: std::sync::Arc<dyn rustfft::Fft<f64>>,
}
impl FusedRows {
    fn new(nz: usize) -> Self {
        let cache = PLAN_CACHE.with(|c| c.borrow().clone());
        let [forward, inverse] = if let Some(cache) = cache {
            cache.rust_pair(nz)
        } else {
            let mut planner = rustfft::FftPlanner::new();
            let inverse = planner.plan_fft_inverse(nz);
            let forward = planner.plan_fft_forward(nz);
            [forward, inverse]
        };
        let scratch = inverse
            .get_immutable_scratch_len()
            .max(forward.get_inplace_scratch_len());
        Self {
            spectral: vec![C::default(); nz],
            real: vec![C::default(); nz],
            scratch: vec![C::default(); scratch],
            inverse,
            forward,
        }
    }
}

pub struct DistributedPairFft {
    shape: [usize; 3],
    first: usize,
    n: [usize; 3],
    total: usize,
    ranks: usize,
    local_positions: Vec<usize>,
    // local stick slot * n0 + first-axis index, one per local coefficient.
    scatter: Vec<usize>,
    local_sticks: usize,
    max_planes: usize,
    block: usize, // complex elements per rank block
    plane_starts: Vec<usize>,
    rank: usize,
    planes: Vec<usize>,
    // (index in the middle batch, owner rank, complex offset of the stick inside
    // one vector's rank block).
    occupied: Vec<(usize, usize, usize)>,
    sticks: Vec<C>,
    middle: Vec<C>,
    slab: Vec<C>,
    send: Vec<f64>,
    recv: Vec<f64>,
    /// 1D kernels: rustfft (default) or MKL DFTI batch lanes (`BAANDI_M2_FFT=mkl`
    /// on MKL builds). Lane buffers: local sticks, middle batch, plane slab.
    lanes: Lanes,
    profile: dft_numerics::real::fft_profile::Probe,
    kept_fuse: bool,
    holes: Option<[HoleClear; 3]>,
    fused_rows: Option<FusedRows>,
}

impl DistributedPairFft {
    /// Plan with the 1D kernels of `BAANDI_M2_FFT` (read once per process).
    pub fn new(
        shape: [usize; 3],
        indices: &[[usize; 3]],
        rank: usize,
        ranks: usize,
    ) -> Result<Self, DistFftError> {
        let vendor = dft_numerics::switches::choice(dft_numerics::switches::Id::M2_FFT) == "mkl";
        Self::with_lanes(shape, indices, rank, ranks, vendor)
    }

    /// Plan with rustfft lanes (`vendor` false) or MKL DFTI batch lanes (`vendor`
    /// true on an MKL build; a pure build has none and uses rustfft, bit-identical).
    pub fn with_lanes(
        shape: [usize; 3],
        indices: &[[usize; 3]],
        rank: usize,
        ranks: usize,
        vendor: bool,
    ) -> Result<Self, DistFftError> {
        let total = shape
            .iter()
            .try_fold(1usize, |a, &b| a.checked_mul(b))
            .ok_or(FftError::InvalidGrid)?;
        if total == 0 || total > isize::MAX as usize / size_of::<C>() || ranks == 0 || rank >= ranks
        {
            return Err(FftError::InvalidGrid.into());
        }
        if indices.iter().any(|p| (0..3).any(|d| p[d] >= shape[d])) {
            return Err(FftError::IncompatibleArray.into());
        }
        let Partition {
            first,
            n,
            sticks,
            planes,
            slots,
            owner,
            local_index,
            counts,
        } = partition(shape, indices, ranks);
        let max_sticks = counts.iter().copied().max().unwrap_or(0);
        let plane_starts: Vec<usize> = (0..=ranks)
            .map(|r| r * (n[0] / ranks) + r.min(n[0] % ranks))
            .collect();
        let max_planes = plane_starts
            .windows(2)
            .map(|w| w[1] - w[0])
            .max()
            .unwrap_or(0);
        let block = max_sticks * max_planes;
        let mut local_positions = Vec::new();
        let mut scatter = Vec::new();
        for (g, (p, &s)) in indices.iter().zip(&slots).enumerate() {
            if owner[s] == rank {
                local_positions.push(g);
                scatter.push(local_index[s] * n[0] + p[first]);
            }
        }
        let occupied = sticks
            .iter()
            .enumerate()
            .map(|(s, &(z, y))| {
                (
                    planes.binary_search(&z).unwrap() * n[1] + y,
                    owner[s],
                    local_index[s] * max_planes,
                )
            })
            .collect::<Vec<_>>();
        let local_sticks = counts[rank];
        let lanes = Lanes::plan(n, [local_sticks, planes.len(), n[1]], vendor)?;
        let holes =
            dft_numerics::switches::on(dft_numerics::switches::Id::FFT_SPARSE_CLEAR).then(|| {
                [
                    HoleClear::new(local_sticks * n[0], scatter.iter().copied()),
                    HoleClear::new(planes.len() * n[1], occupied.iter().map(|p| p.0)),
                    HoleClear::new(
                        n[1] * n[2],
                        (0..n[1]).flat_map(|y| planes.iter().map(move |&z| y * n[2] + z)),
                    ),
                ]
            });
        let config = dft_numerics::switches::RunConfig::current();
        let fused_rows = if config.on(dft_numerics::switches::Id::FFT_FUSED_ROWS) {
            if lanes.is_vendor() {
                config.note_path(
                    dft_numerics::switches::Id::FFT_FUSED_ROWS,
                    "slab",
                    Some("immutable row FFT requires rustfft lanes"),
                );
                None
            } else {
                config.note_path(dft_numerics::switches::Id::FFT_FUSED_ROWS, "rows", None);
                Some(FusedRows::new(n[2]))
            }
        } else {
            None
        };
        Ok(Self {
            fused_rows,
            holes,
            shape,
            first,
            n,
            total,
            ranks,
            local_positions,
            scatter,
            local_sticks,
            max_planes,
            block,
            rank,
            plane_starts,
            middle: vec![C::default(); planes.len() * n[1]],
            sticks: vec![C::default(); local_sticks * n[0]],
            slab: vec![C::default(); n[1] * n[2]],
            planes,
            occupied,
            send: vec![0.; 2 * block * ranks],
            recv: vec![0.; 2 * block * ranks],
            lanes,
            profile: dft_numerics::real::fft_profile::Probe::new(n),
            kept_fuse: dft_numerics::switches::on(dft_numerics::switches::Id::FFT_KEPT_FUSE),
        })
    }

    /// Independent zeroed work buffers with the same layout and shared kernels.
    /// Returns None for the historical private-plan path. No planning or
    /// partitioning is repeated, and no mutable FFT buffer is shared.
    pub fn shared_workspace(&self) -> Option<Self> {
        let Lanes::Shared(lanes) = &self.lanes else {
            return None;
        };
        Some(Self {
            shape: self.shape,
            first: self.first,
            n: self.n,
            total: self.total,
            ranks: self.ranks,
            rank: self.rank,
            local_positions: self.local_positions.clone(),
            scatter: self.scatter.clone(),
            local_sticks: self.local_sticks,
            max_planes: self.max_planes,
            block: self.block,
            plane_starts: self.plane_starts.clone(),
            planes: self.planes.clone(),
            occupied: self.occupied.clone(),
            sticks: vec![C::default(); self.sticks.len()],
            middle: vec![C::default(); self.middle.len()],
            slab: vec![C::default(); self.slab.len()],
            send: vec![0.; self.send.len()],
            recv: vec![0.; self.recv.len()],
            lanes: Lanes::Shared(lanes.clone()),
            profile: dft_numerics::real::fft_profile::Probe::new(self.n),
            kept_fuse: self.kept_fuse,
            holes: self.holes.clone(),
            fused_rows: self.fused_rows.as_ref().map(|r| FusedRows {
                spectral: vec![C::default(); r.spectral.len()],
                real: vec![C::default(); r.real.len()],
                scratch: vec![C::default(); r.scratch.len()],
                inverse: r.inverse.clone(),
                forward: r.forward.clone(),
            }),
        })
    }

    /// Positions in the global coefficient vector owned by this rank (ascending);
    /// local coefficient vectors follow this order.
    pub fn local_positions(&self) -> &[usize] {
        &self.local_positions
    }
    pub fn local_len(&self) -> usize {
        self.scatter.len()
    }
    pub fn shape(&self) -> [usize; 3] {
        self.shape
    }
    /// Number of owned real-space planes (contiguous along the first FFT axis).
    pub fn plane_count(&self) -> usize {
        self.plane_starts[self.rank + 1] - self.plane_starts[self.rank]
    }
    /// Grid index of owned plane i along grid axis `first_axis()`.
    pub fn plane_x(&self, i: usize) -> usize {
        self.plane_starts[self.rank] + i
    }
    pub fn first_axis(&self) -> usize {
        self.first
    }
    /// True when the MKL DFTI lanes are active (`BAANDI_M2_FFT=mkl` on an MKL build).
    pub fn vendor_lanes(&self) -> bool {
        self.lanes.is_vendor()
    }
    /// Offset into a row-major (shape) grid of the z-line at (plane i, row y); the
    /// line has `shape[2]` entries and matches `plane()[y*nz..]`.
    pub fn grid_start(&self, i: usize, y: usize) -> usize {
        let x = self.plane_x(i);
        if self.first == 0 {
            (x * self.shape[1] + y) * self.n[2]
        } else {
            (y * self.shape[1] + x) * self.n[2]
        }
    }
    /// Extent along the second FFT axis (rows of `plane`).
    pub fn rows(&self) -> usize {
        self.n[1]
    }
    /// Current plane slab, row-major [y][z] (set by `to_real_plane`).
    pub fn plane(&self) -> &[C] {
        &self.slab
    }
    pub fn plane_mut(&mut self) -> &mut [C] {
        &mut self.slab
    }

    fn all_to_all<K: Collective>(&mut self, comm: &K, np: usize) -> Result<(), DistFftError> {
        if comm.size() != self.ranks || comm.rank() != self.rank {
            return Err(MpiError::Shape("communicator differs from the plan").into());
        }
        if self.block > 0 {
            let len = 2 * self.block * self.ranks * np;
            comm.alltoall_f64(
                &self.send[..len],
                2 * self.block * np,
                &mut self.recv[..len],
            )?;
        }
        Ok(())
    }

    /// Grow the transpose buffers to batches of `np` vectors (no-op once grown).
    fn reserve_batch(&mut self, np: usize) {
        let len = 2 * self.block * self.ranks * np;
        if self.send.len() < len {
            self.send.resize(len, 0.);
            self.recv.resize(len, 0.);
        }
    }

    /// Complex offset of vector t's block for rank q in a batch of `np`.
    fn base(&self, q: usize, t: usize, np: usize) -> usize {
        (q * np + t) * self.block
    }

    /// Step 1 of G->r: 1D transforms on local sticks and the transpose.
    pub fn begin_to_real<K: Collective>(&mut self, comm: &K, c: &[C]) -> Result<(), DistFftError> {
        self.begin_to_real_batch(comm, c, 1)
    }

    /// Step 1 of G->r for `np` vectors stored one after another in `c` (each
    /// `local_len` long): ONE transpose for all of them.
    pub fn begin_to_real_batch<K: Collective>(
        &mut self,
        comm: &K,
        c: &[C],
        np: usize,
    ) -> Result<(), DistFftError> {
        let ll = self.scatter.len();
        if c.len() != ll * np {
            return Err(FftError::IncompatibleArray.into());
        }
        self.begin_to_real_batch_from(comm, np, |t, g| c[t * ll + g])
    }

    /// Scatter coefficients supplied in local-position order directly into the
    /// sticks. `coefficient(t, g)` is called once per local position per vector,
    /// in the same order as the slice path. No intermediate coefficient array.
    /// All ranks must supply the same batch width; the callback must not perform collectives.
    fn begin_to_real_batch_from<K: Collective>(
        &mut self,
        comm: &K,
        np: usize,
        mut coefficient: impl FnMut(usize, usize) -> C,
    ) -> Result<(), DistFftError> {
        self.reserve_batch(np);
        let (nx, mp) = (self.n[0], self.max_planes);
        for t in 0..np {
            let mark = self.profile.start();
            let clear_points = if let Some(holes) = &self.holes {
                holes[0].clear(&mut self.sticks);
                holes[0].points
            } else {
                self.sticks.fill(C::default());
                self.sticks.len()
            };
            for (g, &i) in self.scatter.iter().enumerate() {
                self.sticks[i] = coefficient(t, g);
            }
            self.profile.finish(
                mark,
                0,
                self.scatter.len(),
                16 * clear_points + 32 * self.scatter.len(),
                None,
            );
            if !self.sticks.is_empty() {
                let mark = self.profile.start();
                self.lanes.run(0, true, &mut self.sticks)?;
                self.profile
                    .finish(mark, 1, self.sticks.len(), 32 * self.sticks.len(), Some(0));
            }
            let mark = self.profile.start();
            for q in 0..self.ranks {
                let (x0, x1) = (self.plane_starts[q], self.plane_starts[q + 1]);
                let qb = self.base(q, t, np);
                for sp in 0..self.local_sticks {
                    let base = 2 * (qb + sp * mp);
                    for (xi, x) in (x0..x1).enumerate() {
                        let v = self.sticks[sp * nx + x];
                        self.send[base + 2 * xi] = v.re;
                        self.send[base + 2 * xi + 1] = v.im;
                    }
                }
            }
            self.profile
                .finish(mark, 2, self.sticks.len(), 32 * self.sticks.len(), None);
        }
        let mark = self.profile.start();
        if self.block > 0 {
            dft_numerics::stage_clock::arrival::site(104);
        }
        let result = self.all_to_all(comm, np);
        if self.block > 0 {
            self.profile.finish(
                mark,
                3,
                self.block * self.ranks * np,
                16 * self.block * self.ranks * np,
                None,
            );
        }
        result
    }

    /// Step 2 of G->r for owned plane i: 2D transform; result in `plane()`.
    pub fn to_real_plane(&mut self, i: usize) -> Result<(), DistFftError> {
        self.to_real_plane_in(0, 1, i)
    }

    /// Step 2 of G->r for vector t of a batch of `np`, owned plane i.
    pub fn to_real_plane_in(&mut self, t: usize, np: usize, i: usize) -> Result<(), DistFftError> {
        let [_, ny, nz] = self.n;
        self.inverse_middle(t, np, i)?;
        let mark = self.profile.start();
        let clear_points = if let Some(holes) = &self.holes {
            holes[2].clear(&mut self.slab);
            holes[2].points
        } else {
            self.slab.fill(C::default());
            self.slab.len()
        };
        for (k, &z) in self.planes.iter().enumerate() {
            for y in 0..ny {
                self.slab[y * nz + z] = self.middle[k * ny + y];
            }
        }
        self.profile.finish(
            mark,
            6,
            self.middle.len(),
            16 * clear_points + 32 * self.middle.len(),
            None,
        );
        let mark = self.profile.start();
        self.lanes.run(2, true, &mut self.slab)?;
        self.profile
            .finish(mark, 7, self.slab.len(), 32 * self.slab.len(), Some(2));
        Ok(())
    }

    /// r->G for owned plane i from `plane()`; results are staged for the transpose.
    pub fn to_reciprocal_plane(&mut self, i: usize) -> Result<(), DistFftError> {
        self.to_reciprocal_plane_in(0, 1, i)
    }

    /// r->G of vector t of a batch of `np`, owned plane i (staged for the transpose).
    pub fn to_reciprocal_plane_in(
        &mut self,
        t: usize,
        np: usize,
        i: usize,
    ) -> Result<(), DistFftError> {
        let [_, ny, nz] = self.n;
        let mark = self.profile.start();
        self.lanes.run(2, false, &mut self.slab)?;
        self.profile
            .finish(mark, 9, self.slab.len(), 32 * self.slab.len(), Some(2));
        let mark = self.profile.start();
        for (k, &z) in self.planes.iter().enumerate() {
            for y in 0..ny {
                self.middle[k * ny + y] = self.slab[y * nz + z];
            }
        }
        self.profile
            .finish(mark, 10, self.middle.len(), 32 * self.middle.len(), None);
        self.forward_middle(t, np, i)?;
        Ok(())
    }

    fn inverse_middle(&mut self, t: usize, np: usize, i: usize) -> Result<(), DistFftError> {
        let mark = self.profile.start();
        let clear_points = if let Some(holes) = &self.holes {
            holes[1].clear(&mut self.middle);
            holes[1].points
        } else {
            self.middle.fill(C::default());
            self.middle.len()
        };
        for &(m, owner, off) in &self.occupied {
            let k = 2 * ((owner * np + t) * self.block + off + i);
            self.middle[m] = C::new(self.recv[k], self.recv[k + 1]);
        }
        self.profile.finish(
            mark,
            4,
            self.occupied.len(),
            16 * clear_points + 32 * self.occupied.len(),
            None,
        );
        if !self.middle.is_empty() {
            let mark = self.profile.start();
            self.lanes.run(1, true, &mut self.middle)?;
            self.profile
                .finish(mark, 5, self.middle.len(), 32 * self.middle.len(), Some(1));
        }
        Ok(())
    }

    fn forward_middle(&mut self, t: usize, np: usize, i: usize) -> Result<(), DistFftError> {
        if !self.middle.is_empty() {
            let mark = self.profile.start();
            self.lanes.run(1, false, &mut self.middle)?;
            self.profile
                .finish(mark, 11, self.middle.len(), 32 * self.middle.len(), Some(1));
        }
        let mark = self.profile.start();
        for &(m, owner, off) in &self.occupied {
            let k = 2 * ((owner * np + t) * self.block + off + i);
            self.send[k] = self.middle[m].re;
            self.send[k + 1] = self.middle[m].im;
        }
        self.profile.finish(
            mark,
            12,
            self.occupied.len(),
            32 * self.occupied.len(),
            None,
        );
        Ok(())
    }

    /// Fused last-axis round trip for an ordinary local-potential application.
    /// Middle entries at row y are dead after the inverse gather and can be
    /// overwritten by the forward scatter without affecting any later row.
    fn multiply_plane_rows(
        &mut self,
        potential: &[f64],
        t: usize,
        np: usize,
        i: usize,
    ) -> Result<(), DistFftError> {
        self.inverse_middle(t, np, i)?;
        let [_, ny, nz] = self.n;
        for y in 0..ny {
            let start = self.grid_start(i, y);
            let row = self.fused_rows.as_mut().expect("fused row plan");
            let mark = self.profile.start();
            for (k, &z) in self.planes.iter().enumerate() {
                row.spectral[z] = self.middle[k * ny + y];
            }
            self.profile
                .finish(mark, 6, self.planes.len(), 32 * self.planes.len(), None);
            let mark = self.profile.start();
            row.inverse.process_immutable_with_scratch(
                &row.spectral,
                &mut row.real,
                &mut row.scratch,
            );
            self.profile.finish(mark, 7, nz, 32 * nz, Some(2));
            let mark = self.profile.start();
            for (a, &v) in row.real.iter_mut().zip(&potential[start..start + nz]) {
                *a *= v;
            }
            self.profile.finish(mark, 8, nz, 40 * nz, None);
            let mark = self.profile.start();
            row.forward
                .process_with_scratch(&mut row.real, &mut row.scratch);
            self.profile.finish(mark, 9, nz, 32 * nz, Some(2));
            let mark = self.profile.start();
            for (k, &z) in self.planes.iter().enumerate() {
                self.middle[k * ny + y] = row.real[z];
            }
            self.profile
                .finish(mark, 10, self.planes.len(), 32 * self.planes.len(), None);
        }
        self.forward_middle(t, np, i)
    }

    /// Transpose back, 1D transforms on local sticks, 1/N scaling; writes local coefficients.
    pub fn finish_to_reciprocal<K: Collective>(
        &mut self,
        comm: &K,
        out: &mut [C],
    ) -> Result<(), DistFftError> {
        self.finish_to_reciprocal_batch(comm, out, 1)
    }

    /// Batched `finish_to_reciprocal`: ONE transpose back for `np` vectors; `out`
    /// receives the np local coefficient vectors one after another.
    pub fn finish_to_reciprocal_batch<K: Collective>(
        &mut self,
        comm: &K,
        out: &mut [C],
        np: usize,
    ) -> Result<(), DistFftError> {
        let ll = self.scatter.len();
        if out.len() != ll * np || self.send.len() < 2 * self.block * self.ranks * np {
            return Err(FftError::IncompatibleArray.into());
        }
        let mark = self.profile.start();
        if self.block > 0 {
            dft_numerics::stage_clock::arrival::site(114);
        }
        self.all_to_all(comm, np)?;
        if self.block > 0 {
            self.profile.finish(
                mark,
                13,
                self.block * self.ranks * np,
                16 * self.block * self.ranks * np,
                None,
            );
        }
        let (nx, mp) = (self.n[0], self.max_planes);
        for t in 0..np {
            let mark = self.profile.start();
            for q in 0..self.ranks {
                let x0 = self.plane_starts[q];
                let qb = self.base(q, t, np);
                for sp in 0..self.local_sticks {
                    let base = 2 * (qb + sp * mp);
                    for xi in 0..self.plane_starts[q + 1] - x0 {
                        self.sticks[sp * nx + x0 + xi] =
                            C::new(self.recv[base + 2 * xi], self.recv[base + 2 * xi + 1]);
                    }
                }
            }
            self.profile
                .finish(mark, 14, self.sticks.len(), 32 * self.sticks.len(), None);
            if !self.sticks.is_empty() {
                let mark = self.profile.start();
                self.lanes.run(0, false, &mut self.sticks)?;
                self.profile
                    .finish(mark, 15, self.sticks.len(), 32 * self.sticks.len(), Some(0));
            }
            let mark = self.profile.start();
            for (o, &i) in out[t * ll..(t + 1) * ll].iter_mut().zip(&self.scatter) {
                *o = self.sticks[i] / self.total as f64;
            }
            self.profile
                .finish(mark, 16, self.scatter.len(), 32 * self.scatter.len(), None);
        }
        Ok(())
    }

    /// Local potential application: G -> r, multiply by `potential` (full row-major
    /// grid, read on owned planes only), r -> G. Equals `PairFusedPlan::apply_local`.
    pub fn apply_local<K: Collective>(
        &mut self,
        comm: &K,
        c: &[C],
        potential: &[f64],
        out: &mut [C],
    ) -> Result<(), DistFftError> {
        let _profile_scope = self.profile.apply_scope();
        if potential.len() != self.total {
            return Err(FftError::IncompatibleArray.into());
        }
        self.begin_to_real(comm, c)?;
        self.multiply_planes(potential, 1)?;
        self.finish_to_reciprocal(comm, out)
    }

    /// `apply_local` for `np` vectors stored one after another in `data`
    /// (overwritten by the result) with one transpose per direction; bit-identical
    /// to np calls of `apply_local`.
    pub fn apply_local_batch<K: Collective>(
        &mut self,
        comm: &K,
        data: &mut [C],
        potential: &[f64],
        np: usize,
    ) -> Result<(), DistFftError> {
        let _profile_scope = self.profile.apply_scope();
        if potential.len() != self.total {
            return Err(FftError::IncompatibleArray.into());
        }
        self.begin_to_real_batch(comm, data, np)?;
        self.multiply_planes(potential, np)?;
        self.finish_to_reciprocal_batch(comm, data, np)
    }

    /// Local-potential application with input packing fused into the stick
    /// scatter. The callback supplies vector t, local position g; output retains
    /// the ordinary local coefficient layout and final 1/N scaling. Callback
    /// source traffic is caller-owned; FFT profiling counts a logical complex
    /// coefficient read, as on the slice path. No allocation after batch warmup.
    pub fn apply_local_batch_from<K: Collective>(
        &mut self,
        comm: &K,
        coefficient: impl FnMut(usize, usize) -> C,
        out: &mut [C],
        potential: &[f64],
        np: usize,
    ) -> Result<(), DistFftError> {
        let _profile_scope = self.profile.apply_scope();
        if potential.len() != self.total || out.len() != self.local_len() * np {
            return Err(FftError::IncompatibleArray.into());
        }
        self.begin_to_real_batch_from(comm, np, coefficient)?;
        self.multiply_planes(potential, np)?;
        self.finish_to_reciprocal_batch(comm, out, np)
    }

    fn multiply_planes(&mut self, potential: &[f64], np: usize) -> Result<(), DistFftError> {
        for t in 0..np {
            for i in 0..self.plane_count() {
                if self.fused_rows.is_some() {
                    self.multiply_plane_rows(potential, t, np, i)?;
                } else {
                    self.to_real_plane_in(t, np, i)?;
                    self.multiply_slab(potential, i);
                    self.to_reciprocal_plane_in(t, np, i)?;
                }
            }
        }
        Ok(())
    }

    /// Multiply the current slab (owned plane i) by `potential` (full row-major grid).
    fn multiply_slab(&mut self, potential: &[f64], i: usize) {
        let mark = self.profile.start();
        let (ny, nz) = (self.n[1], self.n[2]);
        for y in 0..ny {
            let start = self.grid_start(i, y);
            for (a, &v) in self.slab[y * nz..(y + 1) * nz]
                .iter_mut()
                .zip(&potential[start..start + nz])
            {
                *a *= v;
            }
        }
        self.profile
            .finish(mark, 8, self.slab.len(), 40 * self.slab.len(), None);
    }

    /// MPI-KEPT-START: `apply_local_batch` for `np` vectors of which the first
    /// `kept` have their real-space planes supplied in `planes` (vector t, owned
    /// plane i at `planes[t * plane_points() + i * rows() * shape[2]..]`, the
    /// slab `plane()` that `to_real_plane_in` leaves for that vector), so that
    /// only the r -> G transform runs for them; the remaining `np - kept` vectors
    /// are read from `data[kept * local_len()..]` and transformed both ways. All
    /// np results are written to `data` with ONE transpose back. A kept plane
    /// equal to the slab the inverse transform would produce gives the bits of
    /// `apply_local_batch` (the transforms of one vector do not depend on the
    /// batch). Collectives: one transpose to r when `kept < np`, one back; every
    /// rank must call with the same np and kept.
    pub fn apply_local_batch_kept<K: Collective>(
        &mut self,
        comm: &K,
        data: &mut [C],
        potential: &[f64],
        np: usize,
        kept: usize,
        planes: &[C],
    ) -> Result<(), DistFftError> {
        let _profile_scope = self.profile.apply_scope();
        let ll = self.scatter.len();
        let pp = self.plane_points();
        if potential.len() != self.total
            || data.len() != ll * np
            || kept > np
            || planes.len() < kept * pp
        {
            return Err(FftError::IncompatibleArray.into());
        }
        let fresh = np - kept;
        if fresh > 0 {
            self.begin_to_real_batch(comm, &data[kept * ll..], fresh)?;
        }
        // Keeps the received fresh data (resize preserves the contents).
        self.reserve_batch(np);
        let line = self.n[1] * self.n[2];
        for t in 0..kept {
            for i in 0..self.plane_count() {
                let source = &planes[t * pp + i * line..t * pp + (i + 1) * line];
                if self.kept_fuse {
                    // Same scalar operations as copy then multiply_slab. Read the
                    // immutable kept plane directly, avoiding a slab read/write pass.
                    let mark = self.profile.start();
                    let nz = self.n[2];
                    for y in 0..self.n[1] {
                        let start = self.grid_start(i, y);
                        for ((dst, &src), &v) in self.slab[y * nz..(y + 1) * nz]
                            .iter_mut()
                            .zip(&source[y * nz..(y + 1) * nz])
                            .zip(&potential[start..start + nz])
                        {
                            let mut value = src;
                            value *= v;
                            *dst = value;
                        }
                    }
                    self.profile.finish(mark, 8, line, 40 * line, None);
                } else {
                    let mark = self.profile.start();
                    self.slab.copy_from_slice(source);
                    self.profile.finish(mark, 17, line, 32 * line, None);
                    self.multiply_slab(potential, i);
                }
                self.to_reciprocal_plane_in(t, np, i)?;
            }
        }
        for t in 0..fresh {
            for i in 0..self.plane_count() {
                self.to_real_plane_in(t, fresh, i)?;
                self.multiply_slab(potential, i);
                self.to_reciprocal_plane_in(kept + t, np, i)?;
            }
        }
        self.finish_to_reciprocal_batch(comm, data, np)
    }

    /// Grid points on this rank's planes (`plane_count * rows * shape[2]`): the
    /// length of a plane-local field ordered [plane][row][z].
    pub fn plane_points(&self) -> usize {
        self.plane_count() * self.n[1] * self.n[2]
    }

    /// Assemble a full row-major grid field from every rank's plane-local field
    /// (one allgatherv; the copy into the grid is a pure permutation, so the
    /// result does not depend on the rank count). `scratch` grows to the grid size
    /// on first use.
    pub fn gather_planes<K: Collective>(
        &self,
        comm: &K,
        local: &[f64],
        scratch: &mut Vec<f64>,
        full: &mut [f64],
    ) -> Result<(), DistFftError> {
        if comm.size() != self.ranks || comm.rank() != self.rank {
            return Err(MpiError::Shape("communicator differs from the plan").into());
        }
        if local.len() != self.plane_points() || full.len() != self.total {
            return Err(FftError::IncompatibleArray.into());
        }
        let (ny, nz) = (self.n[1], self.n[2]);
        let counts: Vec<usize> = self
            .plane_starts
            .windows(2)
            .map(|w| (w[1] - w[0]) * ny * nz)
            .collect();
        scratch.resize(self.total, 0.);
        comm.allgatherv_f64(local, &counts, scratch)?;
        // Plane x of the concatenation starts at x*ny*nz for every rank count.
        for x in 0..self.n[0] {
            for y in 0..ny {
                let start = if self.first == 0 {
                    (x * self.shape[1] + y) * nz
                } else {
                    (y * self.shape[1] + x) * nz
                };
                full[start..start + nz].copy_from_slice(&scratch[(x * ny + y) * nz..][..nz]);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use dft_mpi::simulate;
    use dft_numerics::fused::FusedOperator;
    use dft_numerics::{Fft3d, FftPlan3d, PairFusedPlan};
    use ndarray::Array3;

    #[test]
    fn a2a_catalog_payload_accounting() {
        use dft_numerics::switches::RunConfig;
        // Geometry/cutoff from benchmarks/cases/si{64,128}-gamma.toml.
        // No SCF or pseudopotential and no inference about remote timing.
        for (shape, edges, npw, sticks, payload) in [
            ([72, 72, 72], [20.522425713435865; 3], 23871, 1005, 145152),
            (
                [72, 72, 144],
                [20.522425713435865, 20.522425713435865, 41.04485142687173],
                47851,
                2005,
                290304,
            ),
        ] {
            let signed = |i: usize, n: usize| {
                if i <= n / 2 {
                    i as f64
                } else {
                    i as f64 - n as f64
                }
            };
            let mut idx = Vec::new();
            for x in 0..shape[0] {
                for y in 0..shape[1] {
                    for z in 0..shape[2] {
                        let pos = [x, y, z];
                        let g2 = (0..3)
                            .map(|a| {
                                (std::f64::consts::TAU * signed(pos[a], shape[a]) / edges[a])
                                    .powi(2)
                            })
                            .sum::<f64>();
                        if g2 <= 30. {
                            idx.push(pos);
                        }
                    }
                }
            }
            assert_eq!(idx.len(), npw);
            let mut original = None;
            for arm in ["off", "on"] {
                let config =
                    RunConfig::from_pairs(&[("fft_a2a_balance", arm), ("m2_sticks", "mirror")])
                        .unwrap();
                let _scope = config.enter();
                let p = partition(shape, &idx, 8);
                assert_eq!(p.sticks.len(), sticks);
                let bytes = 16 * 8 * p.counts.iter().max().unwrap() * p.n[0].div_ceil(8);
                assert_eq!(bytes, payload);
                if let Some(owner) = &original {
                    assert_eq!(&p.owner, owner);
                } else {
                    original = Some(p.owner);
                }
                eprintln!(
                    "P3-A2A catalog {shape:?} {arm}: Npw={npw} S={sticks} counts={:?} payload={bytes} B/pair/rank/direction; useful mean={} B",
                    p.counts,
                    16 * sticks * p.n[0] / 8
                );
            }
        }
    }

    #[test]
    fn a2a_balance_preserves_fft_bits_and_minimizes_padding() {
        type RankResult = (Vec<usize>, usize, Vec<(Vec<C>, Vec<C>)>);
        use dft_numerics::real::fft_profile;
        use dft_numerics::switches::{Id, RunConfig};
        assert!(!RunConfig::from_pairs(&[]).unwrap().on(Id::FFT_A2A_BALANCE));
        let mut fixtures = cases();
        fixtures.extend([
            ([5, 7, 11], sphere([5, 7, 11], 4., [0.; 3])),
            ([7, 11, 13], sphere([7, 11, 13], 9., [0.; 3])),
            ([2, 3, 5], vec![[0, 0, 0]]),
            ([2, 3, 5], vec![]),
        ]);
        let bits = |v: &[C]| {
            v.iter()
                .map(|z| (z.re.to_bits(), z.im.to_bits()))
                .collect::<Vec<_>>()
        };
        let mut reductions = 0;
        for (shape, idx) in fixtures {
            let potential: Vec<_> = (0..shape.iter().product())
                .map(|i| 0.5 + (i as f64).sin())
                .collect();
            for ranks in [1, 2, 3, 8] {
                for mirror in ["greedy", "mirror"] {
                    let run = |arm| {
                        simulate(ranks, |comm| {
                            let config = RunConfig::from_pairs(&[
                                ("fft_a2a_balance", arm),
                                ("m2_sticks", mirror),
                                ("fft_profile", "on"),
                                ("phase_profile", "on"),
                                ("stage_timers", "off"),
                            ])
                            .unwrap();
                            let _scope = config.enter();
                            let owners = position_owners(shape, &idx, ranks).unwrap();
                            let mut p =
                                DistributedPairFft::new(shape, &idx, comm.rank(), ranks).unwrap();
                            for &g in p.local_positions() {
                                assert_eq!(owners[g], comm.rank());
                            }
                            if mirror == "mirror" {
                                for (g, pos) in idx.iter().enumerate() {
                                    let partner =
                                        std::array::from_fn(|a| (shape[a] - pos[a]) % shape[a]);
                                    if let Some(m) = idx.iter().position(|v| *v == partner) {
                                        assert_eq!(owners[g], owners[m]);
                                    }
                                }
                            }
                            let mut outputs = Vec::new();
                            let global_len = idx.len();
                            // Grow, shrink, empty, and reuse the widest storage.
                            p.reserve_batch(5);
                            let pointers = (p.send.as_ptr(), p.recv.as_ptr());
                            for np in [5, 1, 0, 3, 5] {
                                let input = random_c(idx.len() * np, 99);
                                let local: Vec<_> = (0..np)
                                    .flat_map(|t| {
                                        p.local_positions().iter().map({
                                            let input = &input;
                                            move |&g| input[t * global_len + g]
                                        })
                                    })
                                    .collect();
                                let before = fft_profile::snapshot();
                                let mut out = local.clone();
                                p.apply_local_batch(comm, &mut out, &potential, np).unwrap();
                                let delta = fft_profile::snapshot().since(before);
                                for stage in [3, 13] {
                                    assert_eq!(
                                        delta.0[stage][3],
                                        (16 * ranks * p.block * np) as u64
                                    );
                                }
                                // Standalone inverse and kept-plane API must
                                // use precisely the same new ownership layout.
                                p.begin_to_real_batch(comm, &local, np).unwrap();
                                let mut planes = Vec::new();
                                for t in 0..np {
                                    for i in 0..p.plane_count() {
                                        p.to_real_plane_in(t, np, i).unwrap();
                                        planes.extend_from_slice(p.plane());
                                    }
                                }
                                let mut kept = local;
                                p.apply_local_batch_kept(
                                    comm, &mut kept, &potential, np, np, &planes,
                                )
                                .unwrap();
                                assert_eq!(bits(&out), bits(&kept));
                                outputs.push((out, planes));
                            }
                            assert_eq!(pointers, (p.send.as_ptr(), p.recv.as_ptr()));
                            (p.local_positions().to_vec(), p.block, outputs)
                        })
                    };
                    let off = run("off");
                    let on = run("on");
                    assert!(on[0].1 <= off[0].1);
                    if on[0].1 < off[0].1 {
                        reductions += 1;
                        eprintln!(
                            "P3-A2A {shape:?} P={ranks} {mirror}: bytes/pair/rank/direction {} -> {}",
                            16 * ranks * off[0].1,
                            16 * ranks * on[0].1
                        );
                    } else {
                        for (a, b) in off.iter().zip(&on) {
                            assert_eq!(a.0, b.0, "no payload saving: retain original ownership");
                        }
                    }
                    for (j, np) in [5, 1, 0, 3, 5].into_iter().enumerate() {
                        let assemble = |parts: &[RankResult]| {
                            let mut out = vec![C::default(); idx.len() * np];
                            for (positions, _, values) in parts {
                                for t in 0..np {
                                    for (k, &g) in positions.iter().enumerate() {
                                        out[t * idx.len() + g] =
                                            values[j].0[t * positions.len() + k];
                                    }
                                }
                            }
                            bits(&out)
                        };
                        assert_eq!(assemble(&off), assemble(&on));
                        for (a, b) in off.iter().zip(&on) {
                            assert_eq!(bits(&a.2[j].1), bits(&b.2[j].1));
                        }
                    }
                    // Independently compute the optimal mirror-safe capacity.
                    let config = RunConfig::from_pairs(&[("m2_sticks", mirror)]).unwrap();
                    let _scope = config.enter();
                    let part = partition(shape, &idx, ranks);
                    let pairs = if mirror == "mirror" {
                        part.sticks
                            .iter()
                            .enumerate()
                            .filter(|(s, (z, y))| {
                                part.sticks
                                    .binary_search(&(
                                        (part.n[2] - z) % part.n[2],
                                        (part.n[1] - y) % part.n[1],
                                    ))
                                    .is_ok_and(|m| m > *s)
                            })
                            .count()
                    } else {
                        0
                    };
                    let optimum = (0..=part.sticks.len())
                        .find(|&c| ranks * c >= part.sticks.len() && ranks * (c / 2) >= pairs)
                        .unwrap();
                    assert_eq!(on[0].1, optimum * part.n[0].div_ceil(ranks));
                }
            }
        }
        assert!(reductions > 0);
        eprintln!("P3-A2A: {reductions} layouts reduced; FFT max difference 0 (bitwise)");
    }

    /// P3-PERMAP: actual kernel dispatch, including empty owners, odd/prime
    /// lengths and a sheared reciprocal support. MKL: relative 2e-14;
    /// pure: bit-identical on/off. No material-dependent choice is asserted.
    #[test]
    fn kernel_by_length_preserves_transforms() {
        use dft_numerics::switches::{Id, RunConfig};
        let bits = |v: &[C]| {
            v.iter()
                .map(|z| (z.re.to_bits(), z.im.to_bits()))
                .collect::<Vec<_>>()
        };
        assert!(
            !RunConfig::from_pairs(&[])
                .unwrap()
                .on(Id::FFT_KERNEL_BY_LENGTH)
        );
        for (shape, shear) in [
            ([6, 6, 6], [0.; 3]),
            ([5, 7, 11], [0.5, 0.3, 0.]),
            ([3, 5, 24], [0.; 3]),
            ([9, 9, 15], [0.; 3]),
            ([5, 7, 13], [0.; 3]),
        ] {
            let idx = sphere(shape, 3., shear);
            for ranks in [1, 3, 8] {
                simulate(ranks, |comm| {
                    let mut reference = Vec::new();
                    for arm in ["off", "on"] {
                        let cfg = RunConfig::from_pairs(&[
                            ("fft_kernel_by_length", arm),
                            ("fft_fused_rows", "off"),
                        ])
                        .unwrap();
                        let _scope = cfg.enter();
                        let mut p =
                            DistributedPairFft::with_lanes(shape, &idx, comm.rank(), ranks, false)
                                .unwrap();
                        let input = random_c(p.local_len() * 2, 137 + comm.rank() as u64);
                        let mut actual = input.clone();
                        let size = shape.iter().product();
                        p.apply_local_batch(comm, &mut actual, &vec![1.; size], 2)
                            .unwrap();
                        assert!(
                            input.is_empty() || rel(&actual, &input) < 2e-14,
                            "round trip {shape:?}"
                        );
                        let potential: Vec<_> =
                            (0..size).map(|i| (i as f64 * 0.13).cos()).collect();
                        actual.copy_from_slice(&input);
                        p.apply_local_batch(comm, &mut actual, &potential, 2)
                            .unwrap();
                        if arm == "off" {
                            reference = actual;
                        } else {
                            assert!(
                                reference.is_empty() || rel(&actual, &reference) < 2e-14,
                                "apply {shape:?}"
                            );
                            if !p.vendor_lanes() {
                                assert_eq!(bits(&actual), bits(&reference));
                            }
                            let mut again = input.clone();
                            p.apply_local_batch(comm, &mut again, &potential, 2)
                                .unwrap();
                            assert_eq!(bits(&again), bits(&actual));
                        }
                    }
                });
            }
        }
    }

    #[test]
    fn fused_rows_respects_backend_selection() {
        use dft_numerics::switches::{Id, RunConfig};
        let config = RunConfig::from_pairs(&[]).unwrap();
        assert!(!config.on(Id::FFT_FUSED_ROWS));
        let config = RunConfig::from_pairs(&[("fft_fused_rows", "on")]).unwrap();
        let _guard = config.enter();
        let shape = [3, 5, 7];
        let idx = sphere(shape, 4., [0.; 3]);
        simulate(2, |comm| {
            let _guard = config.enter();
            let mut plan =
                DistributedPairFft::with_lanes(shape, &idx, comm.rank(), 2, true).unwrap();
            assert_eq!(plan.fused_rows.is_none(), plan.vendor_lanes());
            let local = random_c(plan.local_len(), 135);
            let mut result = local.clone();
            plan.apply_local_batch(comm, &mut result, &[1.; 105], 1)
                .unwrap();
            assert!(rel(&result, &local) < 2e-14);
        });
    }

    #[test]
    fn fused_rows_preserve_transforms_and_remove_slab_passes() {
        use dft_numerics::{real::fft_profile, switches::RunConfig};
        let mut fixtures = vec![
            ([6, 6, 6], sphere([6, 6, 6], 4., [0.; 3])),
            ([5, 7, 15], sphere([5, 7, 15], 4., [0.; 3])),
            ([5, 7, 11], sphere([5, 7, 11], 5., [0.5, 0.3, 0.])),
            ([3, 5, 7], vec![[0, 0, 0], [0, 1, 0], [0, 3, 0]]),
            ([1, 2, 3], vec![]),
            ([3, 3, 17], sphere([3, 3, 17], 2., [0.; 3])),
            ([3, 3, 36], sphere([3, 3, 36], 2., [0.; 3])),
            ([3, 3, 72], sphere([3, 3, 72], 2., [0.; 3])),
        ];
        fixtures.push((
            [3, 3, 3],
            (0..3)
                .flat_map(|x| (0..3).flat_map(move |y| (0..3).map(move |z| [x, y, z])))
                .collect(),
        ));
        let bits = |a: &[C]| {
            a.iter()
                .map(|v| (v.re.to_bits(), v.im.to_bits()))
                .collect::<Vec<_>>()
        };
        let rel = |a: &[C], b: &[C]| {
            assert_eq!(a.len(), b.len());
            let numerator = a
                .iter()
                .zip(b)
                .map(|(x, y)| (x - y).norm())
                .fold(0_f64, f64::max);
            let denominator = b.iter().map(|x| x.norm()).fold(0_f64, f64::max);
            if denominator == 0. {
                numerator
            } else {
                numerator / denominator
            }
        };
        let mut max_relative = 0_f64;
        let mut changed_bits = 0;
        let mut comparisons = 0;
        for (shape, idx) in fixtures {
            for ranks in [1, 2, 3, 8] {
                for sparse in ["off", "on"] {
                    let results = simulate(ranks, |comm| {
                        let make = |on| {
                            let config = RunConfig::from_pairs(&[
                                ("fft_fused_rows", on),
                                ("fft_sparse_clear", sparse),
                                ("fft_profile", "on"),
                                ("phase_profile", "on"),
                                ("stage_timers", "off"),
                            ])
                            .unwrap();
                            let _guard = config.enter();
                            DistributedPairFft::with_lanes(shape, &idx, comm.rank(), ranks, false)
                                .unwrap()
                        };
                        let (mut old, mut new) = (make("off"), make("on"));
                        let ll = old.local_len();
                        let mut error = 0_f64;
                        let mut changed = 0;
                        // Independent inverse-only check against the public plane API.
                        let c = random_c(ll, 473);
                        old.begin_to_real(comm, &c).unwrap();
                        new.begin_to_real(comm, &c).unwrap();
                        for i in 0..old.plane_count() {
                            old.to_real_plane(i).unwrap();
                            new.inverse_middle(0, 1, i).unwrap();
                            let [_, ny, nz] = new.n;
                            let row = new.fused_rows.as_mut().unwrap();
                            for y in 0..ny {
                                for (k, &z) in new.planes.iter().enumerate() {
                                    row.spectral[z] = new.middle[k * ny + y];
                                }
                                row.inverse.process_immutable_with_scratch(
                                    &row.spectral,
                                    &mut row.real,
                                    &mut row.scratch,
                                );
                                let want = &old.plane()[y * nz..(y + 1) * nz];
                                error = error.max(rel(&row.real, want));
                                changed += usize::from(bits(&row.real) != bits(want));
                            }
                        }
                        // Dirty reuse, growing/shrinking batches, direct callback and
                        // zero width; V=1 also checks the paired transform round trip.
                        for (iteration, np) in [3, 1, 5, 0, 2].into_iter().enumerate() {
                            let potential: Vec<_> = (0..shape.iter().product())
                                .map(|j| {
                                    if iteration == 0 {
                                        1.
                                    } else if j % 5 == 0 {
                                        -0.0
                                    } else {
                                        (j as f64).cos()
                                    }
                                })
                                .collect();
                            let input = random_c(ll * np, 321 + np as u64);
                            let mut expected = input.clone();
                            let before = fft_profile::snapshot();
                            old.apply_local_batch(comm, &mut expected, &potential, np)
                                .unwrap();
                            let old_counts = fft_profile::snapshot().since(before);
                            let row = new.fused_rows.as_ref().unwrap();
                            let pointers = (
                                row.spectral.as_ptr(),
                                row.real.as_ptr(),
                                row.scratch.as_ptr(),
                                new.middle.as_ptr(),
                            );
                            new.slab.fill(C::new(123., -456.));
                            let mut actual = input.clone();
                            let before = fft_profile::snapshot();
                            new.apply_local_batch_from(
                                comm,
                                |t, g| input[t * ll + g],
                                &mut actual,
                                &potential,
                                np,
                            )
                            .unwrap();
                            let counts = fft_profile::snapshot().since(before);
                            error = error.max(rel(&actual, &expected));
                            changed += usize::from(bits(&actual) != bits(&expected));
                            if iteration == 0 {
                                assert!(rel(&actual, &input) < 2e-14, "round trip {shape:?}");
                            }
                            assert!(
                                new.slab.iter().all(|&v| v == C::new(123., -456.)),
                                "no slab materialization"
                            );
                            let row = new.fused_rows.as_ref().unwrap();
                            assert_eq!(
                                pointers,
                                (
                                    row.spectral.as_ptr(),
                                    row.real.as_ptr(),
                                    row.scratch.as_ptr(),
                                    new.middle.as_ptr()
                                )
                            );
                            for z in 0..shape[2] {
                                if !new.planes.contains(&z) {
                                    assert_eq!(row.spectral[z], C::default(), "immutable holes");
                                }
                            }
                            let clear = if sparse == "on" {
                                new.slab.len() - new.middle.len()
                            } else {
                                new.slab.len()
                            };
                            assert_eq!(
                                old_counts.0[6][3] - counts.0[6][3],
                                (16 * np * new.plane_count() * clear) as u64
                            );
                            for stage in 0..fft_profile::STAGES {
                                assert_eq!(
                                    old_counts.0[stage][2], counts.0[stage][2],
                                    "point count stage {stage}"
                                );
                                if stage != 6 {
                                    assert_eq!(
                                        old_counts.0[stage][3], counts.0[stage][3],
                                        "bytes stage {stage}"
                                    );
                                }
                            }
                        }
                        assert!(error < 2e-14, "{shape:?} ranks {ranks}: {error:e}");
                        (error, changed)
                    });
                    for (error, changed) in results {
                        max_relative = max_relative.max(error);
                        changed_bits += changed;
                        comparisons += 5;
                    }
                }
            }
        }
        eprintln!(
            "P3-FFT2: {comparisons} rank-local paired transforms, max relative difference {max_relative:e}, {changed_bits} non-bitwise arrays/rows; slab untouched, zero-hole stores eliminated"
        );
    }

    #[test]
    fn sparse_clear_and_direct_input_preserve_transform_bits() {
        use dft_numerics::{real::fft_profile, switches::RunConfig};
        let mut fixtures = vec![
            ([6, 6, 6], sphere([6, 6, 6], 4., [0.; 3])),
            ([5, 7, 15], sphere([5, 7, 15], 4., [0.; 3])),
            ([5, 7, 9], sphere([5, 7, 9], 5., [0.5, 0.3, 0.])),
            // Force the other first axis and include ranks owning no planes.
            ([3, 5, 7], vec![[0, 0, 0], [0, 1, 0], [0, 3, 0]]),
            ([1, 2, 3], vec![]),
        ];
        fixtures.push((
            [3, 3, 3],
            (0..3)
                .flat_map(|x| (0..3).flat_map(move |y| (0..3).map(move |z| [x, y, z])))
                .collect(),
        ));
        let mut comparisons = 0;
        for (shape, idx) in fixtures {
            let potential: Vec<_> = (0..shape.iter().product())
                .map(|i| if i % 5 == 0 { -0.0 } else { (i as f64).cos() })
                .collect();
            for vendor in [false, true] {
                for ranks in [1, 2, 3, 8] {
                    simulate(ranks, |comm| {
                        let mut plans: Vec<_> = ["off", "on"]
                            .into_iter()
                            .map(|sparse| {
                                let config = RunConfig::from_pairs(&[
                                    ("fft_sparse_clear", sparse),
                                    ("fft_profile", "on"),
                                    ("phase_profile", "on"),
                                    ("stage_timers", "off"),
                                ])
                                .unwrap();
                                let _guard = config.enter();
                                DistributedPairFft::with_lanes(
                                    shape,
                                    &idx,
                                    comm.rank(),
                                    ranks,
                                    vendor,
                                )
                                .unwrap()
                            })
                            .collect();
                        // Reuse dirty buffers, shrink/grow batches, and zero-width batches.
                        for np in [3, 1, 5, 0, 2] {
                            let ll = plans[0].local_len();
                            let input = random_c(ll * np, 321 + np as u64);
                            let mut reference = Vec::new();
                            let mut base_bytes = 0;
                            for (sparse, plan) in plans.iter_mut().enumerate() {
                                for direct in [false, true] {
                                    let before = fft_profile::snapshot();
                                    let mut data = input.clone();
                                    if direct {
                                        plan.apply_local_batch_from(
                                            comm,
                                            |t, g| input[t * ll + g],
                                            &mut data,
                                            &potential,
                                            np,
                                        )
                                        .unwrap();
                                    } else {
                                        plan.apply_local_batch(comm, &mut data, &potential, np)
                                            .unwrap();
                                    }
                                    let delta = fft_profile::snapshot().since(before);
                                    let bits: Vec<_> = data
                                        .iter()
                                        .map(|v| (v.re.to_bits(), v.im.to_bits()))
                                        .collect();
                                    let bytes =
                                        [0, 4, 6].iter().map(|&i| delta.0[i][3]).sum::<u64>();
                                    if sparse == 0 && !direct {
                                        reference = bits;
                                        base_bytes = bytes;
                                    } else {
                                        assert_eq!(
                                            bits, reference,
                                            "{shape:?} ranks {ranks} np {np}"
                                        );
                                        let saved = if sparse == 0 {
                                            0
                                        } else {
                                            16 * np
                                                * (ll
                                                    + plan.plane_count()
                                                        * (plan.occupied.len() + plan.middle.len()))
                                        };
                                        assert_eq!(base_bytes - bytes, saved as u64);
                                    }
                                }
                            }
                        }
                    });
                    comparisons += 5 * 3 * ranks;
                }
            }
        }
        eprintln!(
            "P3-FFT: {comparisons} rank-local on/off transform comparisons, max difference 0 (bitwise); exact sparse-clear byte savings verified"
        );
    }

    #[test]
    fn kept_fusion_bits_counters_and_rank_independence() {
        use dft_numerics::real::fft_profile;
        use dft_numerics::switches::RunConfig;
        let mut fixtures = cases();
        // First axis y, ranks with no planes/sticks, and empty reciprocal support.
        fixtures.push(([3, 5, 7], vec![[0, 0, 0], [0, 1, 0], [0, 2, 0], [0, 3, 0]]));
        fixtures.push(([1, 2, 3], vec![]));
        for (shape, idx) in fixtures {
            let potential: Vec<_> = (0..shape.iter().product())
                .map(|i| if i % 5 == 0 { -0.0 } else { (i as f64).cos() })
                .collect();
            for np in [1, 3, 8] {
                let input = random_c(idx.len() * np, 509);
                let mut rank_one = Vec::new();
                for ranks in [1, 2, 3, 8] {
                    let results = simulate(ranks, |comm| {
                        let cfg = RunConfig::from_pairs(&[
                            ("fft_kept_fuse", "off"),
                            ("fft_profile", "on"),
                            ("phase_profile", "on"),
                            ("stage_timers", "off"),
                        ])
                        .unwrap();
                        let _scope = cfg.enter();
                        let mut p =
                            DistributedPairFft::new(shape, &idx, comm.rank(), ranks).unwrap();
                        let local: Vec<_> = (0..np)
                            .flat_map(|t| {
                                p.local_positions()
                                    .iter()
                                    .map(|&g| input[t * idx.len() + g])
                                    .collect::<Vec<_>>()
                            })
                            .collect();
                        let mut plain = local.clone();
                        p.apply_local_batch(comm, &mut plain, &potential, np)
                            .unwrap();
                        p.begin_to_real_batch(comm, &local, np).unwrap();
                        let pp = p.plane_points();
                        let line = p.rows() * shape[2];
                        let mut planes = vec![C::default(); np * pp];
                        for t in 0..np {
                            for i in 0..p.plane_count() {
                                p.to_real_plane_in(t, np, i).unwrap();
                                planes[t * pp + i * line..t * pp + (i + 1) * line]
                                    .copy_from_slice(p.plane());
                            }
                        }
                        let bits = |xs: &[C]| {
                            xs.iter()
                                .map(|v| (v.re.to_bits(), v.im.to_bits()))
                                .collect::<Vec<_>>()
                        };
                        let plane_bits = bits(&planes);
                        for kept in [0, 1, np] {
                            let mut off_counts = fft_profile::Snapshot::default();
                            for enabled in [false, true] {
                                p.kept_fuse = enabled;
                                let pointers = (p.slab.as_ptr(), p.send.as_ptr(), p.recv.as_ptr());
                                let mut out = local.clone();
                                let before = fft_profile::snapshot();
                                p.apply_local_batch_kept(
                                    comm, &mut out, &potential, np, kept, &planes,
                                )
                                .unwrap();
                                let counts = fft_profile::snapshot().since(before);
                                assert_eq!(bits(&out), bits(&plain));
                                assert_eq!(bits(&planes), plane_bits);
                                assert_eq!(
                                    pointers,
                                    (p.slab.as_ptr(), p.send.as_ptr(), p.recv.as_ptr())
                                );
                                if enabled {
                                    assert_eq!(counts.0[17], [0; 5]);
                                    for stage in 0..17 {
                                        assert_eq!(counts.0[stage], off_counts.0[stage]);
                                    }
                                } else {
                                    assert_eq!(counts.0[17][1], (kept * p.plane_count()) as u64);
                                    assert_eq!(counts.0[17][3], (32 * kept * pp) as u64);
                                    off_counts = counts;
                                }
                            }
                        }
                        (p.local_positions.clone(), plain)
                    });
                    let mut global = vec![C::default(); input.len()];
                    for (positions, values) in results {
                        for t in 0..np {
                            for (l, &g) in positions.iter().enumerate() {
                                global[t * idx.len() + g] = values[t * positions.len() + l];
                            }
                        }
                    }
                    let bits: Vec<_> = global
                        .iter()
                        .map(|v| (v.re.to_bits(), v.im.to_bits()))
                        .collect();
                    if ranks == 1 {
                        rank_one = bits;
                    } else {
                        assert_eq!(bits, rank_one);
                    }
                }
            }
        }
    }

    #[test]
    fn lane_fft_instrumentation_preserves_bits_and_counts() {
        use dft_numerics::real::fft_profile;
        use dft_numerics::switches::RunConfig;
        let mut fixtures = cases();
        fixtures.push(([2, 3, 4], vec![]));
        for (shape, idx) in fixtures {
            let pot: Vec<_> = (0..shape.iter().product())
                .map(|i| (i as f64).cos())
                .collect();
            for ranks in [1, 2, 3, 8] {
                simulate(ranks, |comm| {
                    let mut reference = Vec::new();
                    for (arm, timers) in [("off", "off"), ("on", "off"), ("on", "on")] {
                        let cfg = RunConfig::from_pairs(&[
                            ("fft_profile", arm),
                            ("phase_profile", "on"),
                            ("stage_timers", timers),
                        ])
                        .unwrap();
                        let _scope = cfg.enter();
                        let mut p =
                            DistributedPairFft::new(shape, &idx, comm.rank(), ranks).unwrap();
                        let np = 3;
                        let mut data = random_c(p.local_len() * np, 103);
                        let before = fft_profile::snapshot();
                        p.apply_local_batch(comm, &mut data, &pot, np).unwrap();
                        let delta = fft_profile::snapshot().since(before);
                        if arm == "off" {
                            assert_eq!(delta, fft_profile::Snapshot::default());
                            reference = data;
                        } else {
                            for (a, b) in data.iter().zip(&reference) {
                                assert_eq!(
                                    (a.re.to_bits(), a.im.to_bits()),
                                    (b.re.to_bits(), b.im.to_bits())
                                );
                            }
                            assert_eq!(delta.0[0][1], np as u64);
                            assert_eq!(delta.0[3][1], u64::from(p.block > 0));
                            assert_eq!(delta.0[13][1], u64::from(p.block > 0));
                            assert_eq!(delta.0[8][1], (np * p.plane_count()) as u64);
                            assert_eq!(delta.0[8][2], (np * p.plane_points()) as u64);
                            assert_eq!(delta.0[8][3], (40 * np * p.plane_points()) as u64);
                            if timers == "off" {
                                assert!(delta.0.iter().all(|row| row[0] == 0));
                            }
                        }
                    }
                });
            }
        }
    }

    fn lcg(s: &mut u64) -> f64 {
        *s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*s >> 11) as f64 / (1u64 << 53) as f64) - 0.5
    }

    fn sphere(shape: [usize; 3], radius2: f64, shift: [f64; 3]) -> Vec<[usize; 3]> {
        let f = |i: usize, n: usize| {
            if i <= n / 2 {
                i as f64
            } else {
                i as f64 - n as f64
            }
        };
        let mut v = Vec::new();
        for i in 0..shape[0] {
            for j in 0..shape[1] {
                for k in 0..shape[2] {
                    let g = [
                        f(i, shape[0]) + shift[0],
                        f(j, shape[1]) + shift[1],
                        f(k, shape[2]) + shift[2],
                    ];
                    if g.iter().map(|x| x * x).sum::<f64>() <= radius2 {
                        v.push([i, j, k]);
                    }
                }
            }
        }
        v
    }

    type Case = ([usize; 3], Vec<[usize; 3]>);
    fn cases() -> Vec<Case> {
        let skew: Vec<[usize; 3]> = sphere([12, 14, 16], 20., [0.5, 0.3, 0.])
            .into_iter()
            .filter(|p| (p[0] + 2 * p[1] + p[2]) % 3 != 0)
            .collect();
        vec![
            ([16, 16, 16], sphere([16, 16, 16], 16., [0.; 3])),
            ([15, 15, 45], sphere([15, 15, 45], 49., [0.; 3])),
            ([12, 14, 16], skew),
        ]
    }

    fn random_c(n: usize, seed: u64) -> Vec<C> {
        let mut s = seed;
        (0..n).map(|_| C::new(lcg(&mut s), lcg(&mut s))).collect()
    }

    fn rel(a: &[C], b: &[C]) -> f64 {
        let num = a
            .iter()
            .zip(b)
            .map(|(x, y)| (x - y).norm())
            .fold(0., f64::max);
        let den = b.iter().map(|x| x.norm()).fold(0., f64::max);
        num / den
    }

    /// Distributed G->r field gathered into a full grid.
    fn dist_field(shape: [usize; 3], idx: &[[usize; 3]], c: &[C], ranks: usize) -> Vec<C> {
        let res = simulate(ranks, |comm| {
            let mut p = DistributedPairFft::new(shape, idx, comm.rank(), ranks).unwrap();
            let local: Vec<C> = p.local_positions().iter().map(|&g| c[g]).collect();
            p.begin_to_real(comm, &local).unwrap();
            let mut lines = Vec::new();
            for i in 0..p.plane_count() {
                p.to_real_plane(i).unwrap();
                for y in 0..p.rows() {
                    lines.push((
                        p.grid_start(i, y),
                        p.plane()[y * shape[2]..(y + 1) * shape[2]].to_vec(),
                    ));
                }
            }
            lines
        });
        let mut grid = vec![C::default(); shape.iter().product()];
        let mut covered = 0;
        for lines in &res {
            for (s, l) in lines {
                grid[*s..*s + l.len()].copy_from_slice(l);
                covered += l.len();
            }
        }
        assert_eq!(covered, grid.len());
        grid
    }

    #[test]
    fn partition_is_exact_and_balanced() {
        for (shape, idx) in cases() {
            for ranks in 1..=8 {
                let mut seen = vec![0; idx.len()];
                let mut lens = Vec::new();
                for r in 0..ranks {
                    let p = DistributedPairFft::new(shape, &idx, r, ranks).unwrap();
                    lens.push(p.local_len());
                    for &g in p.local_positions() {
                        seen[g] += 1;
                    }
                }
                assert!(seen.iter().all(|&s| s == 1));
                let max = *lens.iter().max().unwrap();
                let min = *lens.iter().min().unwrap();
                assert!(max - min <= shape[0].max(shape[1]), "{lens:?}");
            }
        }
    }

    #[test]
    fn mirror_partition_keeps_gamma_pairs_on_one_rank() {
        // BAANDI_M2_STICKS=mirror: exact partition for every case; on the
        // Gamma-symmetric spheres +G and -G always share an owner, while the greedy
        // rule splits some pairs for two or more ranks.
        for (case, (shape, idx)) in cases().into_iter().enumerate() {
            let mirror_of = |p: &[usize; 3]| {
                [
                    (shape[0] - p[0]) % shape[0],
                    (shape[1] - p[1]) % shape[1],
                    (shape[2] - p[2]) % shape[2],
                ]
            };
            let symmetric = idx.iter().all(|p| idx.contains(&mirror_of(p)));
            for ranks in 1..=8 {
                for mirror in [false, true] {
                    set_mirror_sticks_override(Some(mirror));
                    let owners = position_owners(shape, &idx, ranks).unwrap();
                    let mut seen = vec![0; idx.len()];
                    for r in 0..ranks {
                        let p = DistributedPairFft::new(shape, &idx, r, ranks).unwrap();
                        for &g in p.local_positions() {
                            seen[g] += 1;
                            assert_eq!(owners[g], r);
                        }
                    }
                    set_mirror_sticks_override(None);
                    assert!(seen.iter().all(|&s| s == 1), "case {case} ranks {ranks}");
                    if !symmetric {
                        continue;
                    }
                    let split = idx
                        .iter()
                        .enumerate()
                        .filter(|(g, p)| {
                            let m = idx.iter().position(|q| *q == mirror_of(p)).unwrap();
                            owners[*g] != owners[m]
                        })
                        .count();
                    if mirror {
                        assert_eq!(split, 0, "case {case} ranks {ranks}");
                    } else if ranks >= 2 {
                        assert!(split > 0, "case {case} ranks {ranks}");
                    }
                }
            }
        }
    }

    #[test]
    fn forward_matches_serial_field() {
        for (shape, idx) in cases() {
            let n = shape.iter().product::<usize>() as f64;
            let c = random_c(idx.len(), 7);
            let mut full = Array3::<C>::zeros(shape);
            for (p, v) in idx.iter().zip(&c) {
                full[[p[0], p[1], p[2]]] = *v;
            }
            let mut serial = FftPlan3d::new(shape).unwrap();
            let mut field = Array3::<C>::zeros(shape);
            serial.to_real(&full, &mut field).unwrap();
            let reference: Vec<C> = field.iter().copied().collect();
            let mut worst = 0f64;
            for ranks in 1..=8 {
                let grid = dist_field(shape, &idx, &c, ranks);
                let e = rel(&grid, &reference);
                worst = worst.max(e);
                assert!(e < 8. * n * f64::EPSILON, "{shape:?} ranks {ranks}: {e:e}");
                assert_eq!(
                    grid,
                    dist_field(shape, &idx, &c, ranks),
                    "bit reproducibility at fixed rank count"
                );
            }
            eprintln!(
                "forward {shape:?} worst rel {worst:e} (limit {:e})",
                8. * n * f64::EPSILON
            );
        }
    }

    #[test]
    fn roundtrip_and_local_potential_match_serial() {
        for (shape, idx) in cases() {
            let n = shape.iter().product::<usize>();
            let c = random_c(idx.len(), 11);
            let mut s = 99;
            let pot: Vec<f64> = (0..n).map(|_| 1.0 + lcg(&mut s)).collect();
            let ones = vec![1.0; n];
            let mut serial = PairFusedPlan::new(shape, &idx).unwrap();
            let mut want = vec![C::default(); idx.len()];
            serial.apply_local(&c, &pot, &mut want).unwrap();
            for ranks in 1..=8 {
                let parts = simulate(ranks, |comm| {
                    let mut p = DistributedPairFft::new(shape, &idx, comm.rank(), ranks).unwrap();
                    let local: Vec<C> = p.local_positions().iter().map(|&g| c[g]).collect();
                    let mut rt = vec![C::default(); local.len()];
                    p.apply_local(comm, &local, &ones, &mut rt).unwrap();
                    let mut out = vec![C::default(); local.len()];
                    p.apply_local(comm, &local, &pot, &mut out).unwrap();
                    (p.local_positions().to_vec(), rt, out)
                });
                let mut got = vec![C::default(); idx.len()];
                let mut rt = vec![C::default(); idx.len()];
                for (pos, r, o) in &parts {
                    for (k, &g) in pos.iter().enumerate() {
                        got[g] = o[k];
                        rt[g] = r[k];
                    }
                }
                let e_rt = rel(&rt, &c);
                let e_pot = rel(&got, &want);
                eprintln!(
                    "{shape:?} ranks {ranks}: roundtrip {e_rt:e} potential {e_pot:e} bitwise {}",
                    got == want
                );
                assert!(e_rt < 4. * n as f64 * f64::EPSILON);
                assert!(e_pot < 1e-12);
                if ranks == 1 {
                    assert_eq!(got, want, "one rank is bit-identical to PairFusedPlan");
                }
            }
        }
    }

    /// `BAANDI_M2_FFT=mkl` lanes: within 1e-14 relative of the rustfft lanes on
    /// every rank count (MKL builds), bit-identical on pure builds (no vendor lanes).
    #[test]
    fn vendor_lanes_match_rust_lanes() {
        for (shape, idx) in cases() {
            let n = shape.iter().product::<usize>();
            let c = random_c(idx.len(), 21);
            let mut s = 17;
            let pot: Vec<f64> = (0..n).map(|_| 1.0 + lcg(&mut s)).collect();
            for ranks in [1usize, 2, 3, 5] {
                let parts = simulate(ranks, |comm| {
                    let mut out = Vec::new();
                    for vendor in [false, true] {
                        let mut p =
                            DistributedPairFft::with_lanes(shape, &idx, comm.rank(), ranks, vendor)
                                .unwrap();
                        assert_eq!(p.vendor_lanes(), vendor && dft_numerics::MKL_BACKEND);
                        let mut data: Vec<C> = p.local_positions().iter().map(|&g| c[g]).collect();
                        let np = 1;
                        p.apply_local_batch(comm, &mut data, &pot, np).unwrap();
                        out.push(data);
                    }
                    out
                });
                for part in &parts {
                    let e = rel(&part[1], &part[0]);
                    assert!(e < 1e-14, "{shape:?} P={ranks}: vendor lanes {e:e}");
                    if !dft_numerics::MKL_BACKEND {
                        assert_eq!(part[0], part[1], "pure build: identical lanes");
                    }
                }
            }
        }
    }

    #[test]
    fn batched_apply_equals_single_calls_bitwise_and_gather_assembles_the_grid() {
        for (shape, idx) in cases() {
            let n = shape.iter().product::<usize>();
            let mut s = 5;
            let pot: Vec<f64> = (0..n).map(|_| 1.0 + lcg(&mut s)).collect();
            let field: Vec<f64> = (0..n).map(|_| lcg(&mut s)).collect();
            let np = 3;
            let vecs: Vec<Vec<C>> = (0..np)
                .map(|t| random_c(idx.len(), 40 + t as u64))
                .collect();
            for ranks in 1..=5 {
                let res = simulate(ranks, |comm| {
                    let mut p = DistributedPairFft::new(shape, &idx, comm.rank(), ranks).unwrap();
                    let ll = p.local_len();
                    let mut single = vec![C::default(); ll * np];
                    for t in 0..np {
                        let local: Vec<C> =
                            p.local_positions().iter().map(|&g| vecs[t][g]).collect();
                        p.apply_local(comm, &local, &pot, &mut single[t * ll..(t + 1) * ll])
                            .unwrap();
                    }
                    let vecs = &vecs;
                    let mut batch: Vec<C> = (0..np)
                        .flat_map(|t| p.local_positions().iter().map(move |&g| vecs[t][g]))
                        .collect();
                    p.apply_local_batch(comm, &mut batch, &pot, np).unwrap();
                    // Plane-local copy of a known field, gathered back.
                    let mut local = Vec::with_capacity(p.plane_points());
                    for i in 0..p.plane_count() {
                        for y in 0..p.rows() {
                            let st = p.grid_start(i, y);
                            local.extend_from_slice(&field[st..st + shape[2]]);
                        }
                    }
                    let mut scratch = Vec::new();
                    let mut full = vec![0.; n];
                    p.gather_planes(comm, &local, &mut scratch, &mut full)
                        .unwrap();
                    (single == batch, full == field)
                });
                assert!(
                    res.iter().all(|r| r.0),
                    "{shape:?} P={ranks}: batch differs"
                );
                assert!(
                    res.iter().all(|r| r.1),
                    "{shape:?} P={ranks}: gather differs"
                );
            }
        }
    }

    /// MPI-KEPT-START: with the planes of the first `kept` vectors taken from
    /// their own inverse transforms, `apply_local_batch_kept` gives the bits of
    /// `apply_local_batch` for every kept count, 1 to 5 ranks.
    #[test]
    fn kept_planes_batch_equals_plain_batch_bitwise() {
        for (shape, idx) in cases() {
            let n = shape.iter().product::<usize>();
            let mut s = 11;
            let pot: Vec<f64> = (0..n).map(|_| lcg(&mut s) - 0.5).collect();
            let np = 3;
            let vecs: Vec<Vec<C>> = (0..np)
                .map(|t| random_c(idx.len(), 70 + t as u64))
                .collect();
            for ranks in 1..=5 {
                let res = simulate(ranks, |comm| {
                    let mut p = DistributedPairFft::new(shape, &idx, comm.rank(), ranks).unwrap();
                    let vecs = &vecs;
                    let input: Vec<C> = (0..np)
                        .flat_map(|t| p.local_positions().iter().map(move |&g| vecs[t][g]))
                        .collect();
                    let mut plain = input.clone();
                    p.apply_local_batch(comm, &mut plain, &pot, np).unwrap();
                    // Planes of every vector from the inverse transform alone.
                    let pp = p.plane_points();
                    let mut planes = vec![C::default(); np * pp];
                    p.begin_to_real_batch(comm, &input, np).unwrap();
                    let line = p.rows() * shape[2];
                    for t in 0..np {
                        for i in 0..p.plane_count() {
                            p.to_real_plane_in(t, np, i).unwrap();
                            planes[t * pp + i * line..t * pp + (i + 1) * line]
                                .copy_from_slice(p.plane());
                        }
                    }
                    (0..=np)
                        .map(|kept| {
                            let mut out = input.clone();
                            p.apply_local_batch_kept(comm, &mut out, &pot, np, kept, &planes)
                                .unwrap();
                            out == plain
                        })
                        .collect::<Vec<_>>()
                });
                for r in &res {
                    assert!(r.iter().all(|&b| b), "{shape:?} P={ranks}: {r:?}");
                }
            }
        }
    }

    #[test]
    fn no_reallocation_after_construction() {
        let (shape, idx) = cases().remove(0);
        simulate(1, |comm| {
            let mut p = DistributedPairFft::new(shape, &idx, 0, 1).unwrap();
            let c = random_c(idx.len(), 3);
            let pot = vec![1.; 4096];
            let mut out = vec![C::default(); idx.len()];
            let ptrs = |p: &DistributedPairFft| {
                [
                    p.send.as_ptr() as usize,
                    p.recv.as_ptr() as usize,
                    p.sticks.as_ptr() as usize,
                    p.middle.as_ptr() as usize,
                ]
            };
            let before = ptrs(&p);
            for _ in 0..3 {
                p.apply_local(comm, &c, &pot, &mut out).unwrap();
            }
            assert_eq!(before, ptrs(&p));
        });
    }
}

#[cfg(test)]
mod partition_cache_tests {
    use super::*;
    use dft_numerics::switches::RunConfig;

    #[test]
    fn exact_partition_key_and_private_workspaces() {
        let shape = [4, 5, 6];
        let indices = vec![[0, 0, 0], [1, 2, 3], [3, 3, 3], [2, 1, 4]];
        let cache = PlanCache::with_partition_cache(true, true);
        let _scope = cache.enter();
        for ranks in [1, 3, 8] {
            for mirror in ["greedy", "mirror"] {
                for balance in ["off", "on"] {
                    let config = RunConfig::from_pairs(&[
                        ("m2_sticks", mirror),
                        ("fft_a2a_balance", balance),
                    ])
                    .unwrap();
                    let _run = config.enter();
                    let expected = partition_uncached(shape, &indices, ranks);
                    let hits = cache.partition_hits();
                    let owners = position_owners(shape, &indices, ranks).unwrap();
                    assert_eq!(
                        cache.partition_hits(),
                        hits,
                        "changed policy/ranks must miss"
                    );
                    assert_eq!(
                        owners,
                        expected
                            .slots
                            .iter()
                            .map(|&s| expected.owner[s])
                            .collect::<Vec<_>>()
                    );
                    let a = DistributedPairFft::new(shape, &indices, 0, ranks).unwrap();
                    assert_eq!(cache.partition_hits(), hits + 1);
                    let mut b = a.shared_workspace().unwrap();
                    b.plane_mut().fill(C::new(1., 2.));
                    assert!(a.plane().iter().all(|v| *v == C::default()));
                }
            }
        }
        let mut changed = indices.clone();
        changed.reverse();
        let hits = cache.partition_hits();
        position_owners(shape, &changed, 8).unwrap();
        assert_eq!(
            cache.partition_hits(),
            hits,
            "ordered indices are part of key"
        );
        position_owners([5, 5, 6], &changed, 8).unwrap();
        assert_eq!(cache.partition_hits(), hits, "shape is part of key");
        let cache_private = PlanCache::with_partition_cache(false, true);
        {
            let _inner = cache_private.enter();
            let p = DistributedPairFft::new(shape, &indices, 0, 3).unwrap();
            assert!(p.shared_workspace().is_none());
            assert_eq!(cache_private.counts(), (0, 0));
        }
        // Nested constructor scopes restore the original cache.
        position_owners([5, 5, 6], &changed, 8).unwrap();
        assert_eq!(cache.partition_hits(), hits + 1);
    }
}
