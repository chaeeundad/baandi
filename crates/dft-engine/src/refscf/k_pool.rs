//! GK9 k-point pools (ADR 0041 option D; `BAANDI_REFSCF_K_POOLS=off|auto|on`,
//! `BAANDI_REFSCF_K_NPOOL=<n>`).
//!
//! Without pools every k point is G/stick-distributed over all ranks and the
//! k points run one after another, so a small-N_pw many-k case pays the
//! collective latency of every Davidson iteration of every k. With pools the
//! ranks form `npool` pools (contiguous rank blocks, sizes
//! `dft_mpi::block_counts(ranks, npool)`, one sub-communicator each by
//! `Collective::split`); pool j owns the contiguous k block
//! `dft_mpi::block_counts(nk, npool)[j]` and distributes each of its k points'
//! plane-wave rows over its ranks exactly as the distributed layout does over
//! all ranks (the same operator, start and Davidson code on the pool
//! communicator). Stage 1 is the special case of one rank per pool (serial
//! collective, full-row operators). Per map the ranks exchange on the world
//! communicator:
//!
//! * the effective potential: one allgatherv of the field planes, then a
//!   local transpose/slice to each operator's own planes (a permutation);
//! * the per-k records (eigenvalues, residual norms, solve records and the
//!   per-band energy terms, reduced inside the pool so they are replicated
//!   there): one allgatherv of the counts (it also carries a failure flag, so
//!   a rank-local error never leaves the others in a collective) and one
//!   allgatherv of the records, sent by the pool leaders;
//! * the per-k densities: one alltoall in which every rank sends the part of
//!   its operator planes of every owned k that each field-plane owner holds;
//!   the owner adds them in global k order.
//!
//! Determinism. Every sum over k (density, band sum, kinetic, nonlocal,
//! residual bound) is replayed in global k order exactly as the distributed
//! loop forms it; per-k quantities are computed by a pool-sized operator. One
//! pool is the distributed layout itself (never routed through this module),
//! one rank never pools; results are in the rounding class across rank and
//! pool counts and bit-reproducible at a fixed (ranks, npool).
//!
//! `auto` (GK9B rule). Candidates: npool dividing the rank count (equal pool
//! sizes) and not above the k count. Predicted k-loop time per map, in units
//! of the work of one plane-wave row:
//!   T(npool) = max over pools j of sum_{k in j} (N_pw,k / p + L ceil(log2 p)),
//! with p = ranks / npool ranks per pool and L = `POOL_LATENCY_ROWS`, the
//! collective latency of one distributed Davidson step expressed in rows of
//! local work (a machine constant, not a material one). The first term is
//! the work split over the pool plus the k load imbalance (the max over
//! pools), the second the cost of distributing a k over p ranks (tree
//! reductions and transposes grow with log2 p). The smallest T wins, ties go
//! to more pools.
use super::{
    checkpoint::{Codec, Reader},
    metric::Serial,
    numerical,
};
use crate::model::ScfError;
use dft_geometry::kpoints::KPoint;
use dft_mpi::{Collective, MpiError};
use dft_numerics::switches::{Id, RunConfig};
use std::ops::Range;
use std::sync::Arc;

/// `L` of the `auto` rule: collective latency of one distributed Davidson step
/// in plane-wave rows of local work. Fitted on the dev08 MPI8 GK9 stage-1 A/B
/// (one pool per rank against the distributed layout: GeTe 260 k, N_pw 4.6k,
/// -11 %; B12 172 k, 6.6k, -6 %): solving T(8)/T(1) for L gives 27 and 24 rows;
/// as the whole SCF (not only the k loop) entered the ratios this is a lower
/// bound, rounded up to 32. Machine-level; recalibrate with the
/// `REFSCF_K_NPOOL` scan on the target machine.
pub const POOL_LATENCY_ROWS: f64 = 32.;

/// Pool sub-communicator of this rank: `None` for a one-rank pool (serial
/// collective), shared so that the start's channels and maps reuse it.
pub(crate) type PoolComm = Option<Option<Arc<dft_mpi::Comm>>>;

/// Ownership of whole k points by pools of ranks.
#[derive(Clone, Debug)]
pub struct KPool {
    /// Global k count.
    pub nk: usize,
    /// Pool of this rank and its rank inside the pool.
    pool: usize,
    rank_in_pool: usize,
    /// Ranks per pool and first world rank of each pool.
    pub pool_sizes: Vec<usize>,
    pub pool_first: Vec<usize>,
    /// k points per pool (`dft_mpi::block_counts(nk, npool)`).
    pub counts: Vec<usize>,
    /// First global k of each pool.
    pub first: Vec<usize>,
    /// The global k list (weights for the k-order replays).
    pub points: Vec<KPoint>,
    /// Plane-wave basis size per global k.
    pub basis_len: Vec<usize>,
    /// Plane axis of every global k operator (set by the start).
    pub op_axis: Vec<usize>,
    comm: Option<Arc<dft_mpi::Comm>>,
}

/// Divisors of `ranks` not above `limit`, ascending.
fn divisors(ranks: usize, limit: usize) -> Vec<usize> {
    (1..=ranks.min(limit))
        .filter(|&d| ranks.is_multiple_of(d))
        .collect()
}

/// Predicted k-loop time of `npool` pools (rows; see the module notes).
pub fn predicted_rows(basis_len: &[usize], ranks: usize, npool: usize) -> f64 {
    let sizes = dft_mpi::block_counts(ranks, npool);
    let counts = dft_mpi::block_counts(basis_len.len(), npool);
    let firsts = dft_mpi::block_displacements(&counts);
    (0..npool)
        .map(|j| {
            let p = sizes[j].max(1);
            let latency = POOL_LATENCY_ROWS * (usize::BITS - (p - 1).leading_zeros()) as f64;
            basis_len[firsts[j]..firsts[j] + counts[j]]
                .iter()
                .map(|&n| n as f64 / p as f64 + if p > 1 { latency } else { 0. })
                .sum::<f64>()
        })
        .fold(0., f64::max)
}

/// The `auto` pool count (module notes).
pub fn auto_npool(basis_len: &[usize], ranks: usize) -> usize {
    let mut best = (f64::INFINITY, 1);
    for npool in divisors(ranks, basis_len.len()) {
        let t = predicted_rows(basis_len, ranks, npool);
        // Ascending candidates: `<=` sends ties to more pools.
        if t <= best.0 {
            best = (t, npool);
        }
    }
    best.1
}

/// Pool count of a run: `REFSCF_K_NPOOL` forces (clamped to ranks and k
/// count), `on` is one rank per pool (largest divisor of the rank count not
/// above the k count when there are fewer k than ranks), `auto` the rule.
pub fn npool_for(choice: &str, forced: Option<f64>, basis_len: &[usize], ranks: usize) -> usize {
    let nk = basis_len.len().max(1);
    if let Some(n) = forced {
        return (n as usize).clamp(1, ranks.min(nk).max(1));
    }
    match choice {
        "on" => divisors(ranks, nk).last().copied().unwrap_or(1),
        "auto" => auto_npool(basis_len, ranks),
        _ => 1,
    }
}

impl KPool {
    /// Decide the layout of one run (collective: every rank calls it with the
    /// same arguments). Pools need the switch not `off`, more than one rank,
    /// a caller that supports them (`allowed`) and a pool count above one;
    /// otherwise the distributed layout stays and the reason is recorded.
    pub(crate) fn decide<K: Collective>(
        points: &[KPoint],
        basis_len: Vec<usize>,
        comm: &K,
        allowed: bool,
    ) -> Result<Option<Self>, ScfError> {
        let config = RunConfig::current();
        let choice = config.choice(Id::REFSCF_K_POOLS);
        if choice == "off" {
            return Ok(None);
        }
        let forced = config.number(Id::REFSCF_K_NPOOL);
        let ranks = comm.size();
        let npool = npool_for(choice, forced, &basis_len, ranks);
        let fallback = if !allowed {
            Some("caller requires the distributed layout (external trace or rotation probe)")
        } else if ranks < 2 {
            Some("one rank")
        } else if npool < 2 {
            Some("one pool")
        } else {
            None
        };
        if let Some(reason) = fallback {
            config.note_path(Id::REFSCF_K_POOLS, "distributed", Some(reason));
            return Ok(None);
        }
        let mut pool = Self::new(points, basis_len, comm.rank(), ranks, npool);
        if pool.pool_sizes.iter().any(|&p| p > 1) {
            // Collective on every rank (MPI_Comm_split), also one-rank pools.
            // A communicator adapter without sub-communicators refuses on
            // every rank alike (the trait default): keep the distributed layout.
            let sub = match comm.split(pool.pool, comm.rank()) {
                Ok(sub) => sub,
                Err(MpiError::State(reason)) => {
                    config.note_path(Id::REFSCF_K_POOLS, "distributed", Some(reason));
                    return Ok(None);
                }
                Err(e) => return Err(numerical(e).into()),
            };
            if sub.size() != pool.pool_size() || sub.rank() != pool.rank_in_pool {
                return Err(ScfError::InvalidInput("k pools: sub-communicator layout"));
            }
            if sub.size() > 1 {
                pool.comm = Some(Arc::new(sub));
            }
        }
        // "pools": one rank per pool; "pools:npool=N": pools of several ranks.
        let path = if pool.pool_sizes.iter().all(|&p| p == 1) {
            "pools".to_string()
        } else {
            format!("pools:npool={npool}")
        };
        config.note_path(Id::REFSCF_K_POOLS, &path, None);
        Ok(Some(pool))
    }
    /// Layout without a sub-communicator (unit tests and one-rank pools).
    pub(crate) fn new(
        points: &[KPoint],
        basis_len: Vec<usize>,
        rank: usize,
        ranks: usize,
        npool: usize,
    ) -> Self {
        let pool_sizes = dft_mpi::block_counts(ranks, npool);
        let pool_first = dft_mpi::block_displacements(&pool_sizes);
        let pool = (0..npool)
            .find(|&j| pool_first[j] <= rank && rank < pool_first[j] + pool_sizes[j])
            .expect("rank within a pool");
        let counts = dft_mpi::block_counts(points.len(), npool);
        let first = dft_mpi::block_displacements(&counts);
        Self {
            nk: points.len(),
            pool,
            rank_in_pool: rank - pool_first[pool],
            pool_sizes,
            pool_first,
            counts,
            first,
            points: points.to_vec(),
            basis_len,
            op_axis: Vec::new(),
            comm: None,
        }
    }
    pub fn npool(&self) -> usize {
        self.counts.len()
    }
    /// Pool of this rank, its rank inside the pool and the pool size.
    pub fn position(&self) -> (usize, usize, usize) {
        (self.pool, self.rank_in_pool, self.pool_size())
    }
    pub fn pool_size(&self) -> usize {
        self.pool_sizes[self.pool]
    }
    /// The pool's rank 0 sends the pool's records.
    pub fn leader(&self) -> bool {
        self.rank_in_pool == 0
    }
    /// Global k indices owned by this rank's pool.
    pub fn local(&self) -> Range<usize> {
        self.first[self.pool]..self.first[self.pool] + self.counts[self.pool]
    }
    /// Pool that owns global k point `k`.
    pub fn owner(&self, k: usize) -> usize {
        (0..self.counts.len())
            .find(|&q| self.first[q] <= k && k < self.first[q] + self.counts[q])
            .expect("k index within the pool")
    }
    pub fn weights(&self) -> Vec<f64> {
        self.points.iter().map(|p| p.weight).collect()
    }
    /// The per-k collective layout of a start (`KComm::new`).
    pub(crate) fn comm_of(pool: Option<&Self>) -> PoolComm {
        pool.map(|p| p.comm.clone())
    }
}

/// The collective of the per-k work: the world communicator without pools,
/// the serial collective for one-rank pools, the pool communicator
/// otherwise. One type, so the per-k code is shared.
pub(crate) enum KComm<'a, K> {
    World(&'a K),
    Own,
    Pool(&'a dft_mpi::Comm),
}
impl<'a, K> KComm<'a, K> {
    pub(crate) fn new(comm: &'a K, layout: &'a PoolComm) -> Self {
        match layout {
            None => Self::World(comm),
            Some(None) => Self::Own,
            Some(Some(c)) => Self::Pool(c),
        }
    }
}
macro_rules! kdispatch {
    ($self:ident, $c:ident => $e:expr) => {
        match $self {
            Self::World($c) => $e,
            Self::Own => {
                let $c = &Serial;
                $e
            }
            Self::Pool($c) => $e,
        }
    };
}
impl<K: Collective> Collective for KComm<'_, K> {
    fn rank(&self) -> usize {
        kdispatch!(self, c => c.rank())
    }
    fn size(&self) -> usize {
        kdispatch!(self, c => c.size())
    }
    fn barrier(&self) -> Result<(), MpiError> {
        kdispatch!(self, c => c.barrier())
    }
    fn allgatherv_f64(
        &self,
        local: &[f64],
        counts: &[usize],
        out: &mut [f64],
    ) -> Result<(), MpiError> {
        kdispatch!(self, c => c.allgatherv_f64(local, counts, out))
    }
    fn tree_sum_f64(&self, local: &[f64], out: &mut [f64]) -> Result<(), MpiError> {
        kdispatch!(self, c => c.tree_sum_f64(local, out))
    }
    fn bcast_f64(&self, data: &mut [f64], root: usize) -> Result<(), MpiError> {
        kdispatch!(self, c => c.bcast_f64(data, root))
    }
    fn alltoall_f64(&self, send: &[f64], block: usize, recv: &mut [f64]) -> Result<(), MpiError> {
        kdispatch!(self, c => c.alltoall_f64(send, block, recv))
    }
}

/// Point count of each rank's field planes (planes along `axis` in {0, 1},
/// the plane rule of `DistributedPairFft`).
fn field_counts(dims: [usize; 3], axis: usize, ranks: usize) -> Vec<usize> {
    let starts = crate::grid_dist::plane_starts(dims[axis], ranks);
    let plane = dims[1 - axis] * dims[2];
    starts.windows(2).map(|w| (w[1] - w[0]) * plane).collect()
}

/// A field distributed on the field planes, gathered once to the full grid
/// (planes along the field axis) and transposed locally on demand to the
/// other plane axis; an operator reads its own planes. Pure copies: values
/// keep their bits.
pub(crate) struct FullField {
    dims: [usize; 3],
    axis: usize,
    full: Vec<f64>,
    other: Option<Vec<f64>>,
}
impl FullField {
    pub(crate) fn gather<K: Collective>(
        comm: &K,
        dims: [usize; 3],
        axis: usize,
        local: &[f64],
    ) -> Result<Self, ScfError> {
        let counts = field_counts(dims, axis, comm.size());
        if local.len() != counts[comm.rank()] {
            return Err(ScfError::InvalidInput("k pools: field plane length"));
        }
        let mut full = vec![0.; dims.iter().product()];
        comm.allgatherv_f64(local, &counts, &mut full)
            .map_err(numerical)?;
        Ok(Self {
            dims,
            axis,
            full,
            other: None,
        })
    }
    /// The full grid in the plane layout of an operator with plane axis `axis`.
    pub(crate) fn for_axis(&mut self, axis: usize) -> Result<&[f64], ScfError> {
        if axis == self.axis {
            return Ok(&self.full);
        }
        if self.other.is_none() {
            let mut map = crate::grid_dist::PlaneMap::new(self.dims, 0, 1, self.axis, axis);
            let mut out = Vec::new();
            map.a_to_b(&Serial, &self.full, &mut out)?;
            self.other = Some(out);
        }
        Ok(self.other.as_deref().expect("transposed"))
    }
    /// The planes of rank `rank` of a `size`-rank operator with plane axis
    /// `axis` (the whole grid for one rank).
    pub(crate) fn for_operator(
        &mut self,
        axis: usize,
        rank: usize,
        size: usize,
    ) -> Result<&[f64], ScfError> {
        let plane = self.dims[1 - axis] * self.dims[2];
        let starts = crate::grid_dist::plane_starts(self.dims[axis], size);
        let full = self.for_axis(axis)?;
        Ok(&full[starts[rank] * plane..starts[rank + 1] * plane])
    }
}

/// The points that source planes `b` (range along axis `b_axis`, operator
/// layout) hold of destination planes `a` (range along `a_axis`, field
/// layout), as runs (destination offset, source offset, length) in
/// destination order. Both layouts: planes, then the other of axes 0/1, then
/// the last axis.
fn runs(
    dims: [usize; 3],
    (a_axis, a0, a1): (usize, usize, usize),
    (b_axis, b0, b1): (usize, usize, usize),
    mut f: impl FnMut(usize, usize, usize),
) {
    let nz = dims[2];
    if a_axis == b_axis {
        let (lo, hi) = (a0.max(b0), a1.min(b1));
        if lo < hi {
            let plane = dims[1 - a_axis] * nz;
            f((lo - a0) * plane, (lo - b0) * plane, (hi - lo) * plane);
        }
    } else {
        for x in a0..a1 {
            for y in b0..b1 {
                f(
                    ((x - a0) * dims[b_axis] + y) * nz,
                    ((y - b0) * dims[a_axis] + x) * nz,
                    nz,
                );
            }
        }
    }
}
fn run_len(dims: [usize; 3], a: (usize, usize, usize), b: (usize, usize, usize)) -> usize {
    let mut n = 0;
    if a.0 == b.0 {
        runs(dims, a, b, |_, _, len| n += len);
    } else {
        n = (a.2 - a.1) * (b.2 - b.1) * dims[2];
    }
    n
}

/// One alltoall of the per-k densities of every owned k (and channel): each
/// rank receives the part of its field planes that every pool rank holds of
/// every k and adds them in global (channel, k) order, the order of the
/// distributed loop.
pub(crate) struct DensityExchange {
    dims: [usize; 3],
    axis: usize,
    channels: usize,
    /// Field plane starts along `axis` over the world ranks.
    starts: Vec<usize>,
    block: usize,
    send: Vec<f64>,
    /// Fill position per destination rank and the next item expected.
    fill: Vec<usize>,
    next: usize,
}
impl DensityExchange {
    /// Operator planes of world rank `r` for global k `k`.
    fn source(pool: &KPool, dims: [usize; 3], r: usize, k: usize) -> (usize, usize, usize) {
        let j = pool.owner(k);
        let b = pool.op_axis[k];
        let s = crate::grid_dist::plane_starts(dims[b], pool.pool_sizes[j]);
        let i = r - pool.pool_first[j];
        (b, s[i], s[i + 1])
    }
    fn dest(&self, q: usize) -> (usize, usize, usize) {
        (self.axis, self.starts[q], self.starts[q + 1])
    }
    pub(crate) fn new(pool: &KPool, dims: [usize; 3], axis: usize, channels: usize) -> Self {
        let ranks: usize = pool.pool_sizes.iter().sum();
        let starts = crate::grid_dist::plane_starts(dims[axis], ranks);
        let mut ex = Self {
            dims,
            axis,
            channels,
            starts,
            block: 0,
            send: Vec::new(),
            fill: Vec::new(),
            next: 0,
        };
        let mut block = 0;
        for j in 0..pool.npool() {
            for r in pool.pool_first[j]..pool.pool_first[j] + pool.pool_sizes[j] {
                for q in 0..ranks {
                    let per_channel: usize = (pool.first[j]..pool.first[j] + pool.counts[j])
                        .map(|k| run_len(dims, ex.dest(q), Self::source(pool, dims, r, k)))
                        .sum();
                    block = block.max(channels * per_channel);
                }
            }
        }
        ex.block = block;
        ex.send = vec![0.; block * ranks];
        ex.fill = (0..ranks).map(|q| q * block).collect();
        ex
    }
    /// Add the density of owned k `local` (index within this pool's block)
    /// and `channel`, on this rank's planes of its operator (plane axis
    /// `op_axis`). Calls come in (channel, local k) order.
    pub(crate) fn add(
        &mut self,
        pool: &KPool,
        channel: usize,
        local: usize,
        op_axis: usize,
        density: &[f64],
    ) -> Result<(), ScfError> {
        let k = pool.local().start + local;
        if channel * pool.counts[pool.pool] + local != self.next
            || pool.op_axis.get(k) != Some(&op_axis)
        {
            return Err(ScfError::InvalidInput(
                "k pools: density exchange order/axis",
            ));
        }
        self.next += 1;
        let me: usize = pool.pool_first[pool.pool] + pool.rank_in_pool;
        let src = Self::source(pool, self.dims, me, k);
        let plane = self.dims[1 - src.0] * self.dims[2];
        if density.len() != (src.2 - src.1) * plane {
            return Err(ScfError::InvalidInput("k pools: density length"));
        }
        for q in 0..self.fill.len() {
            let dest = self.dest(q);
            let mut at = self.fill[q];
            runs(self.dims, dest, src, |_, s, len| {
                self.send[at..at + len].copy_from_slice(&density[s..s + len]);
                at += len;
            });
            self.fill[q] = at;
        }
        Ok(())
    }
    /// Exchange and add into `out[channel]` (this rank's field planes,
    /// already holding the starting values, normally zeros).
    pub(crate) fn finish<K: Collective>(
        self,
        comm: &K,
        pool: &KPool,
        out: &mut [Vec<f64>],
    ) -> Result<(), ScfError> {
        let ranks = self.fill.len();
        let me = comm.rank();
        let dest = self.dest(me);
        let mine = (dest.2 - dest.1) * self.dims[1 - self.axis] * self.dims[2];
        let mut recv = vec![0.; self.block * ranks];
        if self.block > 0 {
            comm.alltoall_f64(&self.send, self.block, &mut recv)
                .map_err(numerical)?;
        }
        if out.len() != self.channels || out.iter().any(|o| o.len() != mine) {
            return Err(ScfError::InvalidInput("k pools: density planes"));
        }
        // Each source rank's items arrive in its own (channel, k) order,
        // which is the visiting order below.
        let mut cursor: Vec<usize> = (0..ranks).map(|r| r * self.block).collect();
        for (channel, density) in out.iter_mut().enumerate() {
            let _ = channel;
            for k in 0..pool.nk {
                let j = pool.owner(k);
                let members = pool.pool_first[j]..pool.pool_first[j] + pool.pool_sizes[j];
                for (r, slot) in cursor
                    .iter_mut()
                    .enumerate()
                    .take(members.end)
                    .skip(members.start)
                {
                    let src = Self::source(pool, self.dims, r, k);
                    let mut at = *slot;
                    runs(self.dims, dest, src, |d, _, len| {
                        for (a, b) in density[d..d + len].iter_mut().zip(&recv[at..at + len]) {
                            *a += b;
                        }
                        at += len;
                    });
                    *slot = at;
                }
            }
        }
        Ok(())
    }
}

/// Agree on a rank-local outcome before the next collective: one allgatherv
/// of a flag; on failure the lowest failing rank broadcasts its message.
/// The failing rank returns its own error, the others a numerical error
/// naming it.
pub(crate) fn agree<K: Collective, T>(comm: &K, local: Result<T, ScfError>) -> Result<T, ScfError> {
    let mut flags = vec![0.; comm.size()];
    comm.allgatherv_f64(
        &[f64::from(u8::from(local.is_err()))],
        &vec![1; comm.size()],
        &mut flags,
    )
    .map_err(numerical)?;
    fail_from(comm, &flags, local)
}
fn fail_from<K: Collective, T>(
    comm: &K,
    failed: &[f64],
    local: Result<T, ScfError>,
) -> Result<T, ScfError> {
    let Some(root) = failed.iter().position(|&f| f != 0.) else {
        return local;
    };
    let message = match &local {
        Err(e) if comm.rank() == root => e.to_string(),
        _ => String::new(),
    };
    let mut len = [message.len() as f64];
    comm.bcast_f64(&mut len, root).map_err(numerical)?;
    let mut bytes: Vec<f64> = message.bytes().map(f64::from).collect();
    bytes.resize(len[0] as usize, 0.);
    comm.bcast_f64(&mut bytes, root).map_err(numerical)?;
    match local {
        Err(e) => Err(e),
        Ok(_) => {
            let text: String = bytes.iter().map(|&b| char::from(b as u8)).collect();
            Err(ScfError::Numerical(format!(
                "k pools: rank {root} failed: {text}"
            )))
        }
    }
}

/// Allgather the records of every pool's owned k points (in rank order of the
/// pool leaders, which is global k order). The count exchange carries the failure flag: a rank
/// whose local work failed sends -1, and every rank returns an error.
/// Ranks with `contribute == false` (pool members other than the leader,
/// whose records duplicate the leader's) send none.
pub(super) fn allgather<K: Collective, T: Codec>(
    comm: &K,
    local: Result<Vec<T>, ScfError>,
    contribute: bool,
) -> Result<Vec<T>, ScfError> {
    let ranks = comm.size();
    let local = local.map(|items| if contribute { items } else { Vec::new() });
    let encoded = local.and_then(|items| {
        let mut bytes = Vec::new();
        items.put(&mut bytes)?;
        Ok(bytes)
    });
    let words = encoded.as_ref().map_or(-1., |b| b.len().div_ceil(8) as f64);
    let mut lengths = vec![0.; ranks];
    comm.allgatherv_f64(&[words], &vec![1; ranks], &mut lengths)
        .map_err(numerical)?;
    let failed: Vec<f64> = lengths
        .iter()
        .map(|&w| f64::from(u8::from(w < 0.)))
        .collect();
    let bytes = fail_from(comm, &failed, encoded)?;
    let counts: Vec<usize> = lengths.iter().map(|&w| w as usize).collect();
    let mut packed: Vec<f64> = bytes
        .chunks(8)
        .map(|c| {
            let mut word = [0u8; 8];
            word[..c.len()].copy_from_slice(c);
            f64::from_bits(u64::from_le_bytes(word))
        })
        .collect();
    packed.resize(counts[comm.rank()], 0.);
    let mut all = vec![0.; counts.iter().sum()];
    comm.allgatherv_f64(&packed, &counts, &mut all)
        .map_err(numerical)?;
    let mut out = Vec::new();
    let mut at = 0;
    for &count in &counts {
        let bytes: Vec<u8> = all[at..at + count]
            .iter()
            .flat_map(|w| w.to_bits().to_le_bytes())
            .collect();
        at += count;
        let mut reader = Reader::new(&bytes);
        let items: Vec<T> = Vec::get(&mut reader)?;
        if reader.remaining() >= 8 {
            return Err(ScfError::InvalidInput("k pools: record exchange length"));
        }
        out.extend(items);
    }
    Ok(out)
}

/// Delivery: the full orbital matrix (rows in basis order) of every global k
/// on every rank. Without pools each k is gathered from its distributed rows
/// as before; with pools the owning pool gathers it on its communicator and
/// its leader broadcasts it.
pub(super) fn deliver_orbitals<K: Collective>(
    state: &crate::refstart::KStartState,
    comm: &K,
    owned: &[&dft_numerics::linalg::ComplexMatrix],
) -> Result<Vec<dft_numerics::linalg::ComplexMatrix>, ScfError> {
    use dft_numerics::linalg::ComplexMatrix;
    use dft_physics::data::Complex64;
    let Some(pool) = &state.pool else {
        return state
            .blocks
            .iter()
            .zip(owned)
            .map(|(block, x)| super::k_driver::gather_orbitals(block, comm, x))
            .collect();
    };
    let layout = KPool::comm_of(Some(pool));
    let kc = KComm::new(comm, &layout);
    // Every pool owns at least one k and every block has the run's width.
    let width = owned.first().map_or(0, |x| x.cols());
    let mut out = Vec::with_capacity(pool.nk);
    let local = pool.local();
    for k in 0..pool.nk {
        let owner = pool.owner(k);
        let n = pool.basis_len[k];
        let mut data = vec![0.; 2 * n * width];
        if local.contains(&k) {
            let li = k - local.start;
            let (block, x) = (&state.blocks[li], owned[li]);
            if x.cols() != width || block.basis.vectors().len() != n {
                return Err(ScfError::InvalidInput("k pools: delivered orbital shape"));
            }
            for j in 0..width {
                for (i, &g) in block.operator.local_rows().iter().enumerate() {
                    data[2 * (g + j * n)] = x[(i, j)].re;
                    data[2 * (g + j * n) + 1] = x[(i, j)].im;
                }
            }
            if pool.pool_size() > 1 {
                let mut all = vec![0.; data.len()];
                kc.tree_sum_f64(&data, &mut all).map_err(numerical)?;
                data = all;
            }
        }
        comm.bcast_f64(&mut data, pool.pool_first[owner])
            .map_err(numerical)?;
        out.push(
            ComplexMatrix::new(
                n,
                width,
                data.as_chunks::<2>()
                    .0
                    .iter()
                    .map(|p| Complex64::new(p[0], p[1]))
                    .collect(),
            )
            .map_err(numerical)?,
        );
    }
    Ok(out)
}

#[cfg(test)]
#[allow(clippy::useless_vec, clippy::needless_range_loop)]
mod tests {
    use super::*;

    fn points(n: usize) -> Vec<KPoint> {
        (0..n)
            .map(|i| KPoint {
                fractional: [0.1 * i as f64, 0., 0.],
                weight: 1. / n as f64,
            })
            .collect()
    }

    #[test]
    fn ownership_is_contiguous_and_exhaustive() {
        for nk in [1usize, 3, 8, 13, 260] {
            for ranks in [1usize, 2, 3, 8] {
                for npool in 1..=ranks.min(nk) {
                    let pts = points(nk);
                    let mut seen = vec![0; nk];
                    let mut leaders = 0;
                    for rank in 0..ranks {
                        let pool = KPool::new(&pts, vec![1; nk], rank, ranks, npool);
                        leaders += usize::from(pool.leader());
                        let (j, i, p) = pool.position();
                        assert_eq!(pool.pool_first[j] + i, rank);
                        assert!(i < p);
                        for k in pool.local() {
                            assert_eq!(pool.owner(k), j, "nk={nk} P={ranks} k={k}");
                            if pool.leader() {
                                seen[k] += 1;
                            }
                        }
                    }
                    assert_eq!(leaders, npool);
                    assert!(seen.iter().all(|&s| s == 1));
                }
            }
        }
    }

    #[test]
    fn auto_rule_trades_imbalance_against_pool_width() {
        // Many k, small N_pw: one rank per pool.
        assert_eq!(auto_npool(&vec![4600; 260], 8), 4);
        assert_eq!(npool_for("on", None, &vec![4600; 260], 8), 8);
        // Few k: pools of several ranks (never more pools than k); three k
        // on eight ranks stay distributed (two pools would hold 2 and 1 k).
        assert_eq!(auto_npool(&vec![2000; 4], 8), 4);
        assert_eq!(auto_npool(&vec![2000; 3], 8), 1);
        assert_eq!(npool_for("on", None, &vec![2000; 3], 8), 2);
        assert_eq!(auto_npool(&vec![2000; 1], 8), 1);
        // Large N_pw with a remainder (51 k on 8 ranks): pools of two ranks
        // remove most of the 7-of-51 imbalance.
        assert_eq!(auto_npool(&vec![15500; 51], 8), 4);
        // Tiny N_pw: latency dominates, one rank per pool.
        assert_eq!(auto_npool(&vec![300; 51], 8), 8);
        // Forced counts clamp to ranks and k.
        assert_eq!(npool_for("auto", Some(16.), &vec![1; 51], 8), 8);
        assert_eq!(npool_for("auto", Some(5.), &vec![1; 3], 8), 3);
        assert_eq!(npool_for("off", None, &vec![1; 51], 8), 1);
        for (n, ranks) in [(51usize, 8usize), (7, 6), (260, 12)] {
            let b = vec![5000; n];
            let best = auto_npool(&b, ranks);
            for d in divisors(ranks, n) {
                assert!(predicted_rows(&b, ranks, best) <= predicted_rows(&b, ranks, d));
            }
        }
    }

    #[test]
    fn density_exchange_sums_in_global_k_order_on_both_axes() {
        let dims = [5, 4, 3];
        let nk = 7;
        let pts = points(nk);
        let full = |k: usize, c: usize| -> Vec<f64> {
            (0..60)
                .map(|i| {
                    0.1 * (k as f64 + 1.)
                        + 1e-17 * (i * (k + 3) + c) as f64
                        + 1e16 * f64::from(u8::from(k % 3 == 1))
                        - 1e16 * f64::from(u8::from(k % 3 == 2))
                })
                .collect()
        };
        for axis in [0usize, 1] {
            // Operator axes vary per k (both kinds of runs in one exchange).
            for op_pattern in 0..3usize {
                let op_axis: Vec<usize> = (0..nk)
                    .map(|k| match op_pattern {
                        0 => axis,
                        1 => 1 - axis,
                        _ => k % 2,
                    })
                    .collect();
                // Reference: the one-rank loop in global k order.
                let mut reference = vec![vec![0.; 60]; 2];
                for (c, r) in reference.iter_mut().enumerate() {
                    for k in 0..nk {
                        let mut field = Vec::new();
                        let mut map = crate::grid_dist::PlaneMap::new(dims, 0, 1, axis, op_axis[k]);
                        map.b_to_a(&Serial, &full(k, c), &mut field).unwrap();
                        for (a, b) in r.iter_mut().zip(field) {
                            *a += b;
                        }
                    }
                }
                for (ranks, npool) in [
                    (2usize, 2usize),
                    (3, 3),
                    (4, 4),
                    (4, 2),
                    (5, 2),
                    (6, 3),
                    (3, 2),
                ] {
                    let got = dft_mpi::simulate(ranks, |comm| {
                        let mut pool = KPool::new(&pts, vec![1; nk], comm.rank(), ranks, npool);
                        pool.op_axis = op_axis.clone();
                        let (_, i, p) = pool.position();
                        let mut ex = DensityExchange::new(&pool, dims, axis, 2);
                        for c in 0..2 {
                            for (li, k) in pool.local().enumerate() {
                                let b = op_axis[k];
                                let plane = dims[1 - b] * dims[2];
                                let s = crate::grid_dist::plane_starts(dims[b], p);
                                // `full` is the k density in layout b (as
                                // the reference reads it): this rank's planes.
                                let mine = &full(k, c)[s[i] * plane..s[i + 1] * plane];
                                ex.add(&pool, c, li, b, mine).unwrap();
                            }
                        }
                        let counts = field_counts(dims, axis, ranks);
                        let mut out = vec![vec![0.; counts[comm.rank()]]; 2];
                        ex.finish(comm, &pool, &mut out).unwrap();
                        out
                    });
                    let counts = field_counts(dims, axis, ranks);
                    let offsets = dft_mpi::block_displacements(&counts);
                    for (r, out) in got.iter().enumerate() {
                        for c in 0..2 {
                            let want = &reference[c][offsets[r]..offsets[r] + counts[r]];
                            assert_eq!(
                                out[c].iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                                want.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                                "axis {axis} ops {op_pattern} P={ranks} npool={npool} rank {r} channel {c}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn full_field_gather_and_transpose_are_permutations() {
        let dims = [4, 5, 2];
        for axis in [0usize, 1] {
            let full: Vec<f64> = (0..40).map(|i| i as f64 + 0.25).collect();
            for ranks in [1usize, 2, 3] {
                let got = dft_mpi::simulate(ranks, |comm| {
                    let counts = field_counts(dims, axis, ranks);
                    let offsets = dft_mpi::block_displacements(&counts);
                    let r = comm.rank();
                    let local = &full[offsets[r]..offsets[r] + counts[r]];
                    let mut f = FullField::gather(comm, dims, axis, local).unwrap();
                    (
                        f.for_axis(axis).unwrap().to_vec(),
                        f.for_axis(1 - axis).unwrap().to_vec(),
                    )
                });
                let mut map = crate::grid_dist::PlaneMap::new(dims, 0, 1, axis, 1 - axis);
                let mut other = Vec::new();
                map.a_to_b(&Serial, &full, &mut other).unwrap();
                for (same, transposed) in got {
                    assert_eq!(same, full);
                    assert_eq!(transposed, other);
                }
            }
        }
    }

    #[test]
    fn records_and_failures_reach_every_rank() {
        let got = dft_mpi::simulate(3, |comm| {
            let mine: Vec<(Vec<f64>, usize)> = (0..comm.rank() + 1)
                .map(|i| (vec![comm.rank() as f64, -0.0, f64::MIN_POSITIVE], i))
                .collect();
            allgather(comm, Ok(mine), true).unwrap()
        });
        for all in &got {
            assert_eq!(all.len(), 6);
            assert_eq!(all[5].0[0], 2.);
            assert_eq!(all[1].0[1].to_bits(), (-0.0f64).to_bits());
        }
        let errors = dft_mpi::simulate(3, |comm| {
            let local: Result<Vec<f64>, _> = if comm.rank() == 1 {
                Err(ScfError::InvalidInput("synthetic failure"))
            } else {
                Ok(vec![1.])
            };
            allgather(comm, local, comm.rank() != 2)
                .unwrap_err()
                .to_string()
        });
        for e in errors {
            assert!(e.contains("synthetic failure"), "{e}");
        }
        // Non-leaders send nothing: only the contributing ranks' records arrive.
        let got = dft_mpi::simulate(4, |comm| {
            allgather(comm, Ok(vec![comm.rank() as f64]), comm.rank() % 2 == 0).unwrap()
        });
        for all in got {
            assert_eq!(all, [0., 2.]);
        }
        let agreed = dft_mpi::simulate(2, |comm| {
            agree(
                comm,
                if comm.rank() == 0 {
                    Ok(())
                } else {
                    Err(ScfError::NonFinite)
                },
            )
            .is_err()
        });
        assert_eq!(agreed, [true, true]);
    }
}
