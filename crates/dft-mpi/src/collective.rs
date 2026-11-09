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
}
/// One simulated rank (see [`simulate`]).
#[derive(Clone)]
pub struct SimRank {
    rank: usize,
    shared: Arc<Shared>,
}
impl SimRank {
    fn publish(&self, data: &[f64]) {
        self.shared.slots.lock().unwrap()[self.rank] = data.to_vec();
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
            let mut all = self.shared.slots.lock().unwrap().clone();
            out.copy_from_slice(&tree_associate(&mut all));
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
    let shared = Arc::new(Shared {
        size: ranks,
        barrier: Barrier::new(ranks),
        slots: Mutex::new(vec![Vec::new(); ranks]),
    });
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
