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

mod collective;
pub use collective::{Collective, SimRank, simulate, tree_associate};
#[cfg(feature = "mpi")]
mod ffi;
#[cfg(feature = "mpi")]
pub use ffi::{Universe, World};
#[cfg(not(feature = "mpi"))]
mod single;
#[cfg(not(feature = "mpi"))]
pub use single::{Universe, World};

/// Message pattern of `World::tree_sum_f64` (`BAANDI_MPI_TREE`, set once by the
/// host after reading the switch). Both give the binomial-tree association bit
/// for bit when the rank count is a power of two (IEEE addition is commutative
/// and each round adds the lower block's partial before the higher block's);
/// other rank counts always use `Binomial`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TreeSum {
    /// Binomial reduce to rank 0, then a broadcast: 2 log2(P) message rounds.
    Binomial,
    /// Recursive doubling (pairwise exchange with rank XOR 2^s): log2(P)
    /// rounds, every rank ends with the sum.
    RecursiveDoubling,
}
static TREE_SUM: std::sync::atomic::AtomicU8 = std::sync::atomic::AtomicU8::new(0);
/// Select the message pattern of `tree_sum_f64` for this process.
pub fn set_tree_sum(pattern: TreeSum) {
    TREE_SUM.store(
        u8::from(pattern == TreeSum::RecursiveDoubling),
        std::sync::atomic::Ordering::Relaxed,
    );
}
/// The message pattern selected by `set_tree_sum` (default `Binomial`).
pub fn tree_sum_pattern() -> TreeSum {
    if TREE_SUM.load(std::sync::atomic::Ordering::Relaxed) == 1 {
        TreeSum::RecursiveDoubling
    } else {
        TreeSum::Binomial
    }
}
/// The arithmetic of the recursive-doubling exchange over contributions in rank
/// order (what every rank holds at the end), for P a power of two: in round s
/// rank r and r XOR s both form (lower partner's partial) + (higher partner's
/// partial). Equals [`tree_associate`] bit for bit (test below).
pub fn recursive_doubling_associate(contributions: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let size = contributions.len();
    assert!(size.is_power_of_two(), "recursive doubling needs 2^k ranks");
    let mut held: Vec<Vec<f64>> = contributions.to_vec();
    let mut step = 1;
    while step < size {
        let before = held.clone();
        for (r, h) in held.iter_mut().enumerate() {
            let partner = r ^ step;
            let (lo, hi) = if r < partner {
                (&before[r], &before[partner])
            } else {
                (&before[partner], &before[r])
            };
            for ((o, a), b) in h.iter_mut().zip(lo).zip(hi) {
                *o = a + b;
            }
        }
        step <<= 1;
    }
    held
}

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
    /// Recursive doubling reproduces the binomial association bit for bit on every
    /// rank for P = 1, 2, 4, 8, 16 (values chosen so that association matters).
    #[test]
    fn recursive_doubling_equals_the_binomial_tree() {
        for size in [1usize, 2, 4, 8, 16] {
            let contributions: Vec<Vec<f64>> = (0..size)
                .map(|r| {
                    (0..5)
                        .map(|i| {
                            let x = (r * 7 + i * 3) as f64;
                            0.1 * x + 1e-17 * x * x - (if r % 3 == 0 { 1e16 } else { 0.0 })
                                + (if r % 3 == 1 { 1e16 } else { 0.0 })
                        })
                        .collect()
                })
                .collect();
            let tree = tree_associate(&mut contributions.clone());
            for held in recursive_doubling_associate(&contributions) {
                assert_eq!(
                    held.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                    tree.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                    "P={size}"
                );
            }
        }
        set_tree_sum(TreeSum::RecursiveDoubling);
        assert_eq!(tree_sum_pattern(), TreeSum::RecursiveDoubling);
        set_tree_sum(TreeSum::Binomial);
        assert_eq!(tree_sum_pattern(), TreeSum::Binomial);
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
