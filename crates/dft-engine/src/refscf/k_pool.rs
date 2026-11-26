//! GK9 stage 1: k-point pools of size one (ADR 0041 option D,
//! `BAANDI_REFSCF_K_POOLS=off|auto|on`, default off).
//!
//! Without pools every k point is G/stick-distributed over all ranks and the
//! k points run one after another, so a small-N_pw many-k case pays the
//! collective latency of every Davidson iteration of every k. With pools each
//! rank owns a contiguous block of whole k points (`dft_mpi::block_counts`,
//! the rule of the MPI-DERIV derivative ranges), builds only their operators
//! with all plane-wave rows (`DistributedComplexHamiltonian` of one rank) and
//! solves them with the serial collective. Per map the ranks exchange:
//!
//! * the effective potential: one allgatherv of the field planes, then a
//!   local transpose to each operator's plane axis (a permutation);
//! * the per-k records (eigenvalues, residual norms, solve records and the
//!   per-band energy terms): one allgatherv of the counts (it also carries a
//!   failure flag, so a rank-local error never leaves the others in a
//!   collective) and one allgatherv of the encoded records;
//! * the per-k densities: one alltoall of the field planes of every owned
//!   k, summed on the owner of each plane in global k order.
//!
//! Determinism. Every per-k quantity is computed by a one-rank operator, so it
//! does not depend on the rank count; every sum over k (density, band sum,
//! kinetic, nonlocal, residual bound) is replayed in global k order exactly as
//! the one-rank loop forms it. What still depends on the rank count is the
//! field layout of the density/potential stages (as without pools): results
//! are in the rounding class across rank counts and bit-reproducible at a
//! fixed rank count. One rank never uses pools (bit-identical to today).
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

/// Ownership of whole k points by ranks (pool size one).
#[derive(Clone, Debug)]
pub struct KPool {
    /// Global k count.
    pub nk: usize,
    rank: usize,
    /// k points per rank (`dft_mpi::block_counts(nk, ranks)`).
    pub counts: Vec<usize>,
    /// First global k of each rank.
    pub first: Vec<usize>,
    /// The global k list (weights for the k-order replays).
    pub points: Vec<KPoint>,
    /// Plane-wave basis size per global k.
    pub basis_len: Vec<usize>,
}
impl KPool {
    /// Decide the layout of one run. Pools need the switch (`auto` or `on`),
    /// more than one rank, at least one k point per rank and a caller that
    /// supports them (`allowed`); otherwise the distributed layout stays and
    /// the reason is recorded on the switch.
    pub(crate) fn decide<K: Collective>(
        points: &[KPoint],
        basis_len: Vec<usize>,
        comm: &K,
        allowed: bool,
    ) -> Option<Self> {
        let config = RunConfig::current();
        let choice = config.choice(Id::REFSCF_K_POOLS);
        if choice == "off" {
            return None;
        }
        let ranks = comm.size();
        let fallback = if !allowed {
            Some("caller requires the distributed layout (external trace or rotation probe)")
        } else if ranks < 2 {
            Some("one rank")
        } else if points.len() < ranks {
            Some("fewer k points than ranks (stage 2: sub-communicators)")
        } else {
            None
        };
        if let Some(reason) = fallback {
            config.note_path(Id::REFSCF_K_POOLS, "distributed", Some(reason));
            return None;
        }
        config.note_path(Id::REFSCF_K_POOLS, "pools", None);
        Some(Self::new(points, basis_len, comm.rank(), ranks))
    }
    pub(crate) fn new(points: &[KPoint], basis_len: Vec<usize>, rank: usize, ranks: usize) -> Self {
        let counts = dft_mpi::block_counts(points.len(), ranks);
        let first = dft_mpi::block_displacements(&counts);
        Self {
            nk: points.len(),
            rank,
            counts,
            first,
            points: points.to_vec(),
            basis_len,
        }
    }
    /// Global k indices owned by this rank.
    pub fn local(&self) -> Range<usize> {
        self.first[self.rank]..self.first[self.rank] + self.counts[self.rank]
    }
    /// Rank that owns global k point `k`.
    pub fn owner(&self, k: usize) -> usize {
        (0..self.counts.len())
            .find(|&q| self.first[q] <= k && k < self.first[q] + self.counts[q])
            .expect("k index within the pool")
    }
    pub fn weights(&self) -> Vec<f64> {
        self.points.iter().map(|p| p.weight).collect()
    }
}

/// The collective of the per-k work: the world communicator without pools,
/// the serial collective with them. One type, so the per-k code is shared.
pub(crate) enum KComm<'a, K> {
    World(&'a K),
    Own,
}
impl<'a, K> KComm<'a, K> {
    pub(crate) fn new(comm: &'a K, pooled: bool) -> Self {
        if pooled { Self::Own } else { Self::World(comm) }
    }
}
impl<K: Collective> Collective for KComm<'_, K> {
    fn rank(&self) -> usize {
        match self {
            Self::World(c) => c.rank(),
            Self::Own => 0,
        }
    }
    fn size(&self) -> usize {
        match self {
            Self::World(c) => c.size(),
            Self::Own => 1,
        }
    }
    fn barrier(&self) -> Result<(), MpiError> {
        match self {
            Self::World(c) => c.barrier(),
            Self::Own => Serial.barrier(),
        }
    }
    fn allgatherv_f64(
        &self,
        local: &[f64],
        counts: &[usize],
        out: &mut [f64],
    ) -> Result<(), MpiError> {
        match self {
            Self::World(c) => c.allgatherv_f64(local, counts, out),
            Self::Own => Serial.allgatherv_f64(local, counts, out),
        }
    }
    fn tree_sum_f64(&self, local: &[f64], out: &mut [f64]) -> Result<(), MpiError> {
        match self {
            Self::World(c) => c.tree_sum_f64(local, out),
            Self::Own => Serial.tree_sum_f64(local, out),
        }
    }
    fn bcast_f64(&self, data: &mut [f64], root: usize) -> Result<(), MpiError> {
        match self {
            Self::World(c) => c.bcast_f64(data, root),
            Self::Own => Serial.bcast_f64(data, root),
        }
    }
    fn alltoall_f64(&self, send: &[f64], block: usize, recv: &mut [f64]) -> Result<(), MpiError> {
        match self {
            Self::World(c) => c.alltoall_f64(send, block, recv),
            Self::Own => Serial.alltoall_f64(send, block, recv),
        }
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
/// other plane axis of a one-rank operator. Pure copies: values keep their bits.
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
}

/// One alltoall of the per-k densities of every owned k (and channel): each
/// rank receives the field planes it owns of every k and adds them in global
/// (channel, k) order, the order of the one-rank loop.
pub(crate) struct DensityExchange {
    dims: [usize; 3],
    axis: usize,
    channels: usize,
    /// Point count of each rank's field planes.
    counts: Vec<usize>,
    offsets: Vec<usize>,
    block: usize,
    send: Vec<f64>,
    full: Vec<f64>,
}
impl DensityExchange {
    pub(crate) fn new(pool: &KPool, dims: [usize; 3], axis: usize, channels: usize) -> Self {
        let ranks = pool.counts.len();
        let counts = field_counts(dims, axis, ranks);
        let offsets = dft_mpi::block_displacements(&counts);
        let items = channels * pool.counts.iter().copied().max().unwrap_or(0);
        let block = items * counts.iter().copied().max().unwrap_or(0);
        Self {
            dims,
            axis,
            channels,
            counts,
            offsets,
            block,
            send: vec![0.; block * ranks],
            full: Vec::new(),
        }
    }
    /// Add the density of owned k `local` (index within this rank's block)
    /// and `channel`, given in the plane layout of its one-rank operator
    /// (plane axis `op_axis`).
    pub(crate) fn add(
        &mut self,
        pool: &KPool,
        channel: usize,
        local: usize,
        op_axis: usize,
        density: &[f64],
    ) -> Result<(), ScfError> {
        let field = if op_axis == self.axis {
            density
        } else {
            let mut map = crate::grid_dist::PlaneMap::new(self.dims, 0, 1, self.axis, op_axis);
            map.b_to_a(&Serial, density, &mut self.full)?;
            &self.full
        };
        if field.len() != self.dims.iter().product::<usize>() {
            return Err(ScfError::InvalidInput("k pools: density length"));
        }
        let item = channel * pool.counts[pool.rank] + local;
        for (q, (&count, &offset)) in self.counts.iter().zip(&self.offsets).enumerate() {
            let at = q * self.block + item * count;
            self.send[at..at + count].copy_from_slice(&field[offset..offset + count]);
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
        let ranks = self.counts.len();
        let mine = self.counts[comm.rank()];
        let mut recv = vec![0.; self.block * ranks];
        if self.block > 0 {
            comm.alltoall_f64(&self.send, self.block, &mut recv)
                .map_err(numerical)?;
        }
        if out.len() != self.channels || out.iter().any(|o| o.len() != mine) {
            return Err(ScfError::InvalidInput("k pools: density planes"));
        }
        for (channel, density) in out.iter_mut().enumerate() {
            for k in 0..pool.nk {
                let r = pool.owner(k);
                let item = channel * pool.counts[r] + (k - pool.first[r]);
                let at = r * self.block + item * mine;
                for (a, b) in density.iter_mut().zip(&recv[at..at + mine]) {
                    *a += b;
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

/// Allgather the records of every rank's owned k points (in rank order, which
/// is global k order). The count exchange carries the failure flag: a rank
/// whose local work failed sends -1, and every rank returns an error.
pub(super) fn allgather<K: Collective, T: Codec>(
    comm: &K,
    local: Result<Vec<T>, ScfError>,
) -> Result<Vec<T>, ScfError> {
    let ranks = comm.size();
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
/// as before; with pools the owner broadcasts its full-row block.
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
    // Every rank owns at least one k and every block has the run's width.
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
        }
        comm.bcast_f64(&mut data, owner).map_err(numerical)?;
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
                let pts = points(nk);
                let mut seen = vec![0; nk];
                for rank in 0..ranks {
                    let pool = KPool::new(&pts, vec![1; nk], rank, ranks);
                    for k in pool.local() {
                        assert_eq!(pool.owner(k), rank, "nk={nk} P={ranks} k={k}");
                        seen[k] += 1;
                    }
                }
                assert!(seen.iter().all(|&s| s == 1));
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
            for op_axis in [0usize, 1] {
                // Serial reference: the one-rank loop in global k order.
                let mut reference = vec![vec![0.; 60]; 2];
                for (c, r) in reference.iter_mut().enumerate() {
                    for k in 0..nk {
                        let mut field = Vec::new();
                        let mut map = crate::grid_dist::PlaneMap::new(dims, 0, 1, axis, op_axis);
                        map.b_to_a(&Serial, &full(k, c), &mut field).unwrap();
                        for (a, b) in r.iter_mut().zip(field) {
                            *a += b;
                        }
                    }
                }
                for ranks in [2usize, 3, 4] {
                    let got = dft_mpi::simulate(ranks, |comm| {
                        let pool = KPool::new(&pts, vec![1; nk], comm.rank(), ranks);
                        let mut ex = DensityExchange::new(&pool, dims, axis, 2);
                        for c in 0..2 {
                            for (li, k) in pool.local().enumerate() {
                                ex.add(&pool, c, li, op_axis, &full(k, c)).unwrap();
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
                                "axis {axis} op {op_axis} P={ranks} rank {r} channel {c}"
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
            allgather(comm, Ok(mine)).unwrap()
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
            allgather(comm, local).unwrap_err().to_string()
        });
        for e in errors {
            assert!(e.contains("synthetic failure"), "{e}");
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
