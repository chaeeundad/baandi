//! Backend-independent collective interface (ADR 0041) and an in-process
//! simulator that runs P "ranks" as threads, so the distributed algorithms can be
//! unit-tested on a machine without MPI. All implementations follow the same
//! fixed association for sums, so results agree bit for bit for a fixed rank count.
use crate::MpiError;
use std::sync::{Arc, Barrier, Mutex};

pub trait Collective {
    fn rank(&self) -> usize;
    fn size(&self) -> usize;
    fn barrier(&self) -> Result<(), MpiError>;
    /// Concatenate every rank's `local` (counts[r] items from rank r) into `out`.
    fn allgatherv_f64(
        &self,
        local: &[f64],
        counts: &[usize],
        out: &mut [f64],
    ) -> Result<(), MpiError>;
    /// Elementwise sum over ranks with the fixed binomial-tree association;
    /// every rank receives identical bits.
    fn tree_sum_f64(&self, local: &[f64], out: &mut [f64]) -> Result<(), MpiError>;
    /// Broadcast from `root`.
    fn bcast_f64(&self, data: &mut [f64], root: usize) -> Result<(), MpiError>;
    /// Equal-block all-to-all: block r of `send` goes to rank r, block r of
    /// `recv` comes from rank r.
    fn alltoall_f64(&self, send: &[f64], block: usize, recv: &mut [f64]) -> Result<(), MpiError>;
    /// `MPI_Comm_split` semantics: ranks passing the same `color` form one
    /// sub-communicator, ranked by `key` and then by their rank here.
    /// Collective. The default serves one-rank collectives only; adapters
    /// that wrap a multi-rank communicator forward it.
    fn split(&self, color: usize, key: usize) -> Result<Comm, MpiError> {
        let _ = (color, key);
        if self.size() == 1 {
            Ok(Comm::Single)
        } else {
            Err(MpiError::State(
                "this collective adapter does not support sub-communicators",
            ))
        }
    }
}

/// A sub-communicator made by [`Collective::split`] (GK9 k-point pools).
pub enum Comm {
    /// A group of one rank.
    Single,
    /// A group of simulated ranks ([`simulate`]).
    Sim(SimRank),
    /// An MPI communicator from `MPI_Comm_split`.
    #[cfg(feature = "mpi")]
    Mpi(crate::ffi::SubComm),
}
impl std::fmt::Debug for Comm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Comm(rank {} of {})", self.rank(), self.size())
    }
}
macro_rules! dispatch {
    ($self:ident, $c:ident => $e:expr, $single:expr) => {
        match $self {
            Comm::Single => $single,
            Comm::Sim($c) => $e,
            #[cfg(feature = "mpi")]
            Comm::Mpi($c) => $e,
        }
    };
}
impl Collective for Comm {
    fn rank(&self) -> usize {
        dispatch!(self, c => c.rank(), 0)
    }
    fn size(&self) -> usize {
        dispatch!(self, c => c.size(), 1)
    }
    fn barrier(&self) -> Result<(), MpiError> {
        dispatch!(self, c => c.barrier(), Ok(()))
    }
    fn allgatherv_f64(
        &self,
        local: &[f64],
        counts: &[usize],
        out: &mut [f64],
    ) -> Result<(), MpiError> {
        dispatch!(self, c => c.allgatherv_f64(local, counts, out), {
            if counts != [local.len()] || out.len() != local.len() {
                return Err(MpiError::Shape("allgatherv counts"));
            }
            out.copy_from_slice(local);
            Ok(())
        })
    }
    fn tree_sum_f64(&self, local: &[f64], out: &mut [f64]) -> Result<(), MpiError> {
        dispatch!(self, c => c.tree_sum_f64(local, out), {
            if out.len() != local.len() {
                return Err(MpiError::Shape("tree_sum lengths"));
            }
            out.copy_from_slice(local);
            Ok(())
        })
    }
    fn bcast_f64(&self, data: &mut [f64], root: usize) -> Result<(), MpiError> {
        dispatch!(self, c => c.bcast_f64(data, root), if root == 0 {
            Ok(())
        } else {
            Err(MpiError::Shape("root out of range"))
        })
    }
    fn alltoall_f64(&self, send: &[f64], block: usize, recv: &mut [f64]) -> Result<(), MpiError> {
        dispatch!(self, c => c.alltoall_f64(send, block, recv), {
            if send.len() != block || recv.len() != block {
                return Err(MpiError::Shape("alltoall lengths"));
            }
            recv.copy_from_slice(send);
            Ok(())
        })
    }
    fn split(&self, color: usize, key: usize) -> Result<Comm, MpiError> {
        match self {
            Comm::Single => Ok(Comm::Single),
            Comm::Sim(c) => Collective::split(c, color, key),
            #[cfg(feature = "mpi")]
            Comm::Mpi(c) => Ok(Comm::Mpi(c.split(color, key)?)),
        }
    }
}

/// Members of `color` ordered by (key, rank) from the gathered (color, key)
/// pairs of every rank; `MPI_Comm_split` ordering.
pub fn split_members(pairs: &[(usize, usize)], color: usize) -> Vec<usize> {
    let mut members: Vec<usize> = (0..pairs.len()).filter(|&r| pairs[r].0 == color).collect();
    members.sort_by_key(|&r| (pairs[r].1, r));
    members
}

impl Collective for crate::World {
    fn rank(&self) -> usize {
        crate::World::rank(self)
    }
    fn size(&self) -> usize {
        crate::World::size(self)
    }
    fn barrier(&self) -> Result<(), MpiError> {
        crate::World::barrier(self)
    }
    fn allgatherv_f64(
        &self,
        local: &[f64],
        counts: &[usize],
        out: &mut [f64],
    ) -> Result<(), MpiError> {
        crate::World::allgatherv_f64(self, local, counts, out)
    }
    fn tree_sum_f64(&self, local: &[f64], out: &mut [f64]) -> Result<(), MpiError> {
        crate::World::tree_sum_f64(self, local, out)
    }
    fn bcast_f64(&self, data: &mut [f64], root: usize) -> Result<(), MpiError> {
        crate::World::bcast_f64(self, data, root)
    }
    fn alltoall_f64(&self, send: &[f64], block: usize, recv: &mut [f64]) -> Result<(), MpiError> {
        crate::World::alltoall_f64(self, send, block, recv)
    }
    fn split(&self, color: usize, key: usize) -> Result<Comm, MpiError> {
        #[cfg(feature = "mpi")]
        {
            Ok(Comm::Mpi(crate::World::split(self, color, key)?))
        }
        #[cfg(not(feature = "mpi"))]
        {
            let _ = (color, key);
            Ok(Comm::Single)
        }
    }
}

/// The binomial-tree association over contributions held in rank order
/// (the arithmetic of `World::tree_sum_f64`, computed locally).
pub fn tree_associate(contributions: &mut [Vec<f64>]) -> Vec<f64> {
    let size = contributions.len();
    let mut step = 1;
    while step < size {
        for r in (0..size).step_by(2 * step) {
            if r + step < size {
                let (left, right) = contributions.split_at_mut(r + step);
                for (a, b) in left[r].iter_mut().zip(&right[0]) {
                    *a += *b;
                }
            }
        }
        step <<= 1;
    }
    contributions[0].clone()
}

struct Shared {
    size: usize,
    barrier: Barrier,
    slots: Mutex<Vec<Vec<f64>>>,
    /// Hand-over of sub-group state during `split` (slot of the group leader).
    groups: Mutex<Vec<Option<Arc<Shared>>>>,
}
impl Shared {
    fn new(size: usize) -> Self {
        Self {
            size,
            barrier: Barrier::new(size),
            slots: Mutex::new(vec![Vec::new(); size]),
            groups: Mutex::new(vec![None; size]),
        }
    }
}
/// One simulated rank (see [`simulate`]).
#[derive(Clone)]
pub struct SimRank {
    rank: usize,
    shared: Arc<Shared>,
}
impl SimRank {
    fn publish(&self, data: &[f64]) {
        {
            // Reuse the slot's storage (grow-only), so a steady-state collective of
            // the simulator allocates nothing.
            let mut slots = self.shared.slots.lock().unwrap();
            let slot = &mut slots[self.rank];
            slot.clear();
            slot.extend_from_slice(data);
        }
        self.shared.barrier.wait();
    }
    fn release(&self) {
        self.shared.barrier.wait();
    }
}
impl Collective for SimRank {
    fn rank(&self) -> usize {
        self.rank
    }
    fn size(&self) -> usize {
        self.shared.size
    }
    fn barrier(&self) -> Result<(), MpiError> {
        self.shared.barrier.wait();
        Ok(())
    }
    fn allgatherv_f64(
        &self,
        local: &[f64],
        counts: &[usize],
        out: &mut [f64],
    ) -> Result<(), MpiError> {
        if counts.len() != self.size()
            || counts[self.rank] != local.len()
            || counts.iter().sum::<usize>() != out.len()
        {
            return Err(MpiError::Shape("allgatherv counts"));
        }
        self.publish(local);
        {
            let slots = self.shared.slots.lock().unwrap();
            let mut at = 0;
            for (r, &c) in counts.iter().enumerate() {
                if slots[r].len() != c {
                    return Err(MpiError::Shape("allgatherv peer count"));
                }
                out[at..at + c].copy_from_slice(&slots[r]);
                at += c;
            }
        }
        self.release();
        Ok(())
    }
    fn tree_sum_f64(&self, local: &[f64], out: &mut [f64]) -> Result<(), MpiError> {
        if out.len() != local.len() {
            return Err(MpiError::Shape("tree_sum lengths"));
        }
        self.publish(local);
        {
            // The association of `tree_associate`, element by element on a
            // stack copy of the P contributions (no clone of the slots).
            let slots = self.shared.slots.lock().unwrap();
            let size = slots.len();
            if size > 64 {
                let mut all = slots.clone();
                out.copy_from_slice(&tree_associate(&mut all));
            } else {
                let mut acc = [0.0f64; 64];
                for (i, o) in out.iter_mut().enumerate() {
                    for (a, slot) in acc.iter_mut().zip(slots.iter()) {
                        *a = slot[i];
                    }
                    let mut step = 1;
                    while step < size {
                        for r in (0..size).step_by(2 * step) {
                            if r + step < size {
                                acc[r] += acc[r + step];
                            }
                        }
                        step <<= 1;
                    }
                    *o = acc[0];
                }
            }
        }
        self.release();
        Ok(())
    }
    fn bcast_f64(&self, data: &mut [f64], root: usize) -> Result<(), MpiError> {
        if root >= self.size() {
            return Err(MpiError::Shape("root out of range"));
        }
        self.publish(data);
        {
            let slots = self.shared.slots.lock().unwrap();
            data.copy_from_slice(&slots[root]);
        }
        self.release();
        Ok(())
    }
    fn split(&self, color: usize, key: usize) -> Result<Comm, MpiError> {
        let size = self.size();
        let mut all = vec![0.; 2 * size];
        self.allgatherv_f64(&[color as f64, key as f64], &vec![2; size], &mut all)?;
        let pairs: Vec<(usize, usize)> = all
            .as_chunks::<2>()
            .0
            .iter()
            .map(|p| (p[0] as usize, p[1] as usize))
            .collect();
        let members = split_members(&pairs, color);
        let leader = members[0];
        if leader == self.rank {
            self.shared.groups.lock().unwrap()[leader] = Some(Arc::new(Shared::new(members.len())));
        }
        self.shared.barrier.wait();
        let shared = self.shared.groups.lock().unwrap()[leader]
            .clone()
            .expect("group published by its leader");
        self.shared.barrier.wait();
        Ok(Comm::Sim(SimRank {
            rank: members
                .iter()
                .position(|&r| r == self.rank)
                .expect("member"),
            shared,
        }))
    }
    fn alltoall_f64(&self, send: &[f64], block: usize, recv: &mut [f64]) -> Result<(), MpiError> {
        let size = self.size();
        if send.len() != block * size || recv.len() != block * size {
            return Err(MpiError::Shape("alltoall lengths"));
        }
        self.publish(send);
        {
            let slots = self.shared.slots.lock().unwrap();
            for r in 0..size {
                recv[r * block..(r + 1) * block]
                    .copy_from_slice(&slots[r][self.rank * block..(self.rank + 1) * block]);
            }
        }
        self.release();
        Ok(())
    }
}

/// Run `work` on `ranks` simulated ranks (threads) and return the per-rank results
/// in rank order.
pub fn simulate<R: Send>(ranks: usize, work: impl Fn(&SimRank) -> R + Sync) -> Vec<R> {
    let shared = Arc::new(Shared::new(ranks));
    let work = &work;
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..ranks)
            .map(|rank| {
                let comm = SimRank {
                    rank,
                    shared: shared.clone(),
                };
                scope.spawn(move || work(&comm))
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{block_counts, block_displacements};

    /// Sub-groups behave as independent communicators with MPI_Comm_split
    /// ordering (key, then rank); their collectives agree with the serial
    /// association of their members, and groups can be split again.
    #[test]
    fn split_groups_are_independent_communicators() {
        for ranks in [1usize, 2, 3, 5, 8] {
            let out = simulate(ranks, |comm| {
                let color = comm.rank() % 2;
                let key = ranks - comm.rank();
                let sub = comm.split(color, key).unwrap();
                let v = [comm.rank() as f64 + 0.25, 1e16 * (comm.rank() % 3) as f64];
                let mut s = [0.; 2];
                sub.tree_sum_f64(&v, &mut s).unwrap();
                let mut all = vec![0.; sub.size()];
                sub.allgatherv_f64(&[comm.rank() as f64], &vec![1; sub.size()], &mut all)
                    .unwrap();
                let inner = sub.split(0, sub.rank()).unwrap();
                (sub.rank(), sub.size(), s, all, inner.size())
            });
            for (r, (sub_rank, sub_size, s, all, inner)) in out.iter().enumerate() {
                let mut members: Vec<usize> = (0..ranks).filter(|q| q % 2 == r % 2).collect();
                members.reverse(); // key = ranks - rank: descending rank
                assert_eq!(*sub_size, members.len());
                assert_eq!(*inner, members.len());
                assert_eq!(members[*sub_rank], r);
                assert_eq!(all, &members.iter().map(|&q| q as f64).collect::<Vec<_>>());
                let mut c: Vec<Vec<f64>> = members
                    .iter()
                    .map(|&q| vec![q as f64 + 0.25, 1e16 * (q % 3) as f64])
                    .collect();
                assert_eq!(s.to_vec(), tree_associate(&mut c));
            }
        }
    }

    #[test]
    fn collectives_agree_between_ranks_and_with_the_serial_association() {
        for ranks in [1usize, 2, 3, 4, 5, 8] {
            let results = simulate(ranks, |comm| {
                let v: Vec<f64> = (0..5)
                    .map(|i| 0.1 * (comm.rank() + 1) as f64 + 1e-9 * i as f64)
                    .collect();
                let mut sum = vec![0.0; 5];
                comm.tree_sum_f64(&v, &mut sum).unwrap();
                let counts = block_counts(23, ranks);
                let local: Vec<f64> = (0..counts[comm.rank()])
                    .map(|i| (comm.rank() * 100 + i) as f64)
                    .collect();
                let mut all = vec![0.0; 23];
                comm.allgatherv_f64(&local, &counts, &mut all).unwrap();
                let send: Vec<f64> = (0..2 * ranks)
                    .map(|i| (comm.rank() * 1000 + i) as f64)
                    .collect();
                let mut recv = vec![0.0; 2 * ranks];
                comm.alltoall_f64(&send, 2, &mut recv).unwrap();
                (sum, all, recv)
            });
            let mut contributions: Vec<Vec<f64>> = (0..ranks)
                .map(|r| {
                    (0..5)
                        .map(|i| 0.1 * (r + 1) as f64 + 1e-9 * i as f64)
                        .collect()
                })
                .collect();
            let expect = tree_associate(&mut contributions);
            let counts = block_counts(23, ranks);
            let displs = block_displacements(&counts);
            for (rank, (sum, all, recv)) in results.iter().enumerate() {
                assert_eq!(
                    sum.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                    expect.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
                );
                for r in 0..ranks {
                    for i in 0..counts[r] {
                        assert_eq!(all[displs[r] + i], (r * 100 + i) as f64);
                    }
                    for k in 0..2 {
                        assert_eq!(recv[2 * r + k], (r * 1000 + 2 * rank + k) as f64);
                    }
                }
            }
        }
    }
}
