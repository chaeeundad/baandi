//! Thin MPI adapter (ADR 0041, milestone M0). The FFI follows the Open MPI ABI
//! (pointer handles, `ompi_mpi_*` predefined objects) and is only compiled with
//! feature `mpi`; without it `World` is a single-rank stand-in whose collectives
//! have the same semantics (rank 0 of 1), so dependent code runs unchanged.
//!
//! Reductions of floating-point data are deliberately NOT exposed as
//! `MPI_Allreduce`: [`World::ordered_sum`] gathers the contributions and adds them
//! locally in rank order, which makes the result independent of the library's
//! internal reduction tree (determinism contract of ADR 0041).
use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MpiError {
    /// An MPI call returned this nonzero error code.
    Call(&'static str, i32),
    /// Buffer lengths, counts or displacements do not match.
    Shape(&'static str),
    /// MPI was already initialised or finalised elsewhere.
    State(&'static str),
}
impl fmt::Display for MpiError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Call(name, code) => write!(f, "{name} failed with MPI error code {code}"),
            Self::Shape(what) => write!(f, "invalid MPI buffer shape: {what}"),
            Self::State(what) => write!(f, "invalid MPI state: {what}"),
        }
    }
}
impl std::error::Error for MpiError {}

#[cfg(feature = "mpi")]
mod ffi;
#[cfg(feature = "mpi")]
pub use ffi::{Universe, World};
#[cfg(not(feature = "mpi"))]
mod single;
#[cfg(not(feature = "mpi"))]
pub use single::{Universe, World};

/// Row counts of a block distribution: `total` items over `ranks`, the first
/// `total % ranks` ranks get one extra. A fixed function of (total, ranks).
pub fn block_counts(total: usize, ranks: usize) -> Vec<usize> {
    let base = total / ranks.max(1);
    let extra = total % ranks.max(1);
    (0..ranks).map(|r| base + usize::from(r < extra)).collect()
}
/// Start offsets matching [`block_counts`].
pub fn block_displacements(counts: &[usize]) -> Vec<usize> {
    let mut at = 0;
    counts
        .iter()
        .map(|c| {
            let d = at;
            at += c;
            d
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn block_distribution_is_fixed_and_exhaustive() {
        for total in [0usize, 1, 7, 68, 11877] {
            for ranks in [1usize, 2, 3, 8] {
                let c = block_counts(total, ranks);
                assert_eq!(c.len(), ranks);
                assert_eq!(c.iter().sum::<usize>(), total);
                assert!(c.iter().max().unwrap() - c.iter().min().unwrap() <= 1);
                let d = block_displacements(&c);
                assert_eq!(d[0], 0);
            }
        }
    }
    #[test]
    fn stand_in_is_a_single_rank() {
        #[cfg(not(feature = "mpi"))]
        {
            let universe = Universe::init().unwrap();
            let world = universe.world();
            assert_eq!((world.rank(), world.size()), (0, 1));
            let mut out = vec![0.0; 3];
            world
                .allgatherv_f64(&[1.0, 2.0, 3.0], &[3], &mut out)
                .unwrap();
            assert_eq!(out, [1.0, 2.0, 3.0]);
            let mut sum = vec![0.0; 2];
            world.ordered_sum_f64(&[4.0, 5.0], &mut sum).unwrap();
            assert_eq!(sum, [4.0, 5.0]);
            world.barrier().unwrap();
        }
    }
}
