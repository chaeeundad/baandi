//! Single-rank stand-in used when feature `mpi` is off.
use crate::MpiError;

pub struct Universe;
#[derive(Clone, Copy)]
pub struct World;
impl Universe {
    pub fn init() -> Result<Self, MpiError> {
        Ok(Self)
    }
    pub fn world(&self) -> World {
        World
    }
}
impl World {
    pub fn rank(&self) -> usize {
        0
    }
    pub fn size(&self) -> usize {
        1
    }
    pub fn barrier(&self) -> Result<(), MpiError> {
        Ok(())
    }
    pub fn wtime(&self) -> f64 {
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0.0, |d| d.as_secs_f64())
    }
    /// Broadcast from `root` (must be 0).
    pub fn bcast_f64(&self, _data: &mut [f64], root: usize) -> Result<(), MpiError> {
        if root != 0 {
            return Err(MpiError::Shape("root out of range"));
        }
        Ok(())
    }
    /// Gather every rank's `local` (counts[r] items from rank r) into `out`.
    pub fn allgatherv_f64(
        &self,
        local: &[f64],
        counts: &[usize],
        out: &mut [f64],
    ) -> Result<(), MpiError> {
        if counts.len() != 1 || counts[0] != local.len() || out.len() != local.len() {
            return Err(MpiError::Shape("allgatherv counts"));
        }
        out.copy_from_slice(local);
        Ok(())
    }
    /// Elementwise sum over ranks in rank order (bitwise independent of the
    /// MPI library's reduction tree).
    pub fn ordered_sum_f64(&self, local: &[f64], out: &mut [f64]) -> Result<(), MpiError> {
        if out.len() != local.len() {
            return Err(MpiError::Shape("ordered_sum lengths"));
        }
        out.copy_from_slice(local);
        Ok(())
    }
    /// Fixed-order tree sum (single rank: the input).
    pub fn tree_sum_f64(&self, local: &[f64], out: &mut [f64]) -> Result<(), MpiError> {
        self.ordered_sum_f64(local, out)
    }
    /// All-to-all with equal blocks of `block` items per rank pair.
    pub fn alltoall_f64(
        &self,
        send: &[f64],
        block: usize,
        recv: &mut [f64],
    ) -> Result<(), MpiError> {
        if send.len() != block || recv.len() != block {
            return Err(MpiError::Shape("alltoall lengths"));
        }
        recv.copy_from_slice(send);
        Ok(())
    }
}
