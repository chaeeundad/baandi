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

/// Whether this binary links the real MPI backend (including rank-one runs).
pub const MPI_ENABLED: bool = cfg!(feature = "mpi");

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
    /// REL-TREE-RHD: for payloads of at least [`halving_min_doubles`] doubles a
    /// reduce-scatter by recursive halving (partner rank XOR 2^s in the same
    /// order 1, 2, 4, ... as recursive doubling, each rank keeps the lower or
    /// upper half of the shared segment by bit s of its rank) followed by an
    /// allgather by recursive doubling (pure copies, reverse order): 2 log2(P)
    /// rounds moving 2 (P - 1) / P of the payload instead of log2(P) times the
    /// payload. Every element is associated exactly as in recursive doubling
    /// ([`halving_doubling_associate`]), so the bits are those of the binomial
    /// tree. Smaller payloads use recursive doubling.
    HalvingDoubling,
}
/// Smallest payload (doubles) for which `HalvingDoubling` uses the halving
/// pattern (below it the doubling pattern, latency-bound). A message-size
/// crossover of the machine, not of the physics: calibrate with
/// `examples/mpi_tree_bench.rs` on the target node.
pub const HALVING_MIN_DOUBLES: usize = 2048;
/// Sentinel: choose the crossover from the communicator size (see
/// [`halving_min_doubles_for`]).
pub const HALVING_AUTO: usize = usize::MAX;
static HALVING_MIN: std::sync::atomic::AtomicUsize =
    std::sync::atomic::AtomicUsize::new(HALVING_AUTO);
/// Override the halving crossover of this process (measurement only: the
/// microbenchmark forces 0 to time and bit-check the halving pattern on every
/// payload; any value gives the same bits).
pub fn set_halving_min_doubles(doubles: usize) {
    HALVING_MIN.store(doubles, std::sync::atomic::Ordering::Relaxed);
}
/// The halving crossover override in use ([`HALVING_AUTO`] = size dependent).
pub fn halving_min_doubles() -> usize {
    HALVING_MIN.load(std::sync::atomic::Ordering::Relaxed)
}
/// Crossover for a communicator of `size` ranks. dev08 microbenchmark
/// (mpi_tree_bench, Open MPI 5.0.8, 2026-11-18; us per call, worst rank): halving
/// beats doubling from about 55000 doubles at 2 ranks, 4096 at 4 ranks and 8192 at
/// 8 ranks; at the former fixed 2048 it was 28-59 % slower than doubling. Any
/// value gives the same bits.
pub fn halving_min_doubles_for(size: usize) -> usize {
    let forced = halving_min_doubles();
    if forced != HALVING_AUTO {
        return forced;
    }
    match size {
        0..=2 => 55_000,
        3..=4 => 4_096,
        _ => 8_192,
    }
}
static TREE_SUM: std::sync::atomic::AtomicU8 = std::sync::atomic::AtomicU8::new(0);
/// Select the message pattern of `tree_sum_f64` for this process.
pub fn set_tree_sum(pattern: TreeSum) {
    TREE_SUM.store(
        match pattern {
            TreeSum::Binomial => 0,
            TreeSum::RecursiveDoubling => 1,
            TreeSum::HalvingDoubling => 2,
        },
        std::sync::atomic::Ordering::Relaxed,
    );
}
/// The message pattern selected by `set_tree_sum` (default `Binomial`).
pub fn tree_sum_pattern() -> TreeSum {
    match TREE_SUM.load(std::sync::atomic::Ordering::Relaxed) {
        1 => TreeSum::RecursiveDoubling,
        2 => TreeSum::HalvingDoubling,
        _ => TreeSum::Binomial,
    }
}
/// The arithmetic of the halving/doubling exchange (`TreeSum::HalvingDoubling`
/// on payloads it splits) over contributions in rank order, for P a power of
/// two: the segment bookkeeping and the additions of the reduce-scatter, then
/// the copies of the allgather, exactly as `World::tree_sum_f64` performs them.
pub fn halving_doubling_associate(contributions: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let size = contributions.len();
    assert!(size.is_power_of_two(), "halving needs 2^k ranks");
    let n = contributions.first().map_or(0, Vec::len);
    let mut held: Vec<Vec<f64>> = contributions.to_vec();
    // Reduce-scatter: (segment before the round, kept segment) per rank.
    let mut seg: Vec<(usize, usize)> = vec![(0, n); size];
    let mut history: Vec<Vec<(usize, usize)>> = Vec::new();
    let mut step = 1;
    while step < size {
        history.push(seg.clone());
        let before = held.clone();
        for r in 0..size {
            let (lo, hi) = seg[r];
            let (keep, _) = halving_split(lo, hi, r & step != 0);
            let partner = r ^ step;
            assert_eq!(seg[partner], (lo, hi), "partners share the segment");
            for i in keep.0..keep.1 {
                held[r][i] = before[r][i] + before[partner][i];
            }
        }
        for (r, s) in seg.iter_mut().enumerate() {
            *s = halving_split(s.0, s.1, r & step != 0).0;
        }
        step <<= 1;
    }
    // Allgather: reverse rounds, the partner's kept segment is copied in.
    while let Some(parent) = history.pop() {
        step >>= 1;
        let before = held.clone();
        for (r, h) in held.iter_mut().enumerate() {
            let partner = r ^ step;
            let (lo, hi) = seg[partner];
            h[lo..hi].copy_from_slice(&before[partner][lo..hi]);
        }
        seg = parent;
    }
    held
}
/// (kept, sent) halves of the shared segment [lo, hi) in a halving round: the
/// rank with the round's bit clear keeps the lower half.
pub(crate) fn halving_split(lo: usize, hi: usize, upper: bool) -> ((usize, usize), (usize, usize)) {
    let mid = lo + (hi - lo) / 2;
    if upper {
        ((mid, hi), (lo, mid))
    } else {
        ((lo, mid), (mid, hi))
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
        set_tree_sum(TreeSum::HalvingDoubling);
        assert_eq!(tree_sum_pattern(), TreeSum::HalvingDoubling);
        set_tree_sum(TreeSum::Binomial);
        assert_eq!(tree_sum_pattern(), TreeSum::Binomial);
    }
    /// REL-TREE-RHD: halving/doubling reproduces the binomial association bit
    /// for bit on every rank for P = 1..16 (powers of two) and payload lengths
    /// shorter than, equal to and longer than P, odd and even (segments of
    /// length zero included), with values where the association matters.
    #[test]
    fn halving_doubling_equals_the_binomial_tree() {
        for size in [1usize, 2, 4, 8, 16] {
            for n in [0usize, 1, 2, 3, 5, 7, 8, 15, 16, 17, 33, 100] {
                let contributions: Vec<Vec<f64>> = (0..size)
                    .map(|r| {
                        (0..n)
                            .map(|i| {
                                let x = (r * 7 + i * 3) as f64;
                                0.1 * x + 1e-17 * x * x - (if r % 3 == 0 { 1e16 } else { 0.0 })
                                    + (if (r + i) % 3 == 1 { 1e16 } else { 0.0 })
                            })
                            .collect()
                    })
                    .collect();
                let tree = tree_associate(&mut contributions.clone());
                for (r, held) in halving_doubling_associate(&contributions)
                    .iter()
                    .enumerate()
                {
                    assert_eq!(
                        held.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                        tree.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                        "P={size} n={n} rank={r}"
                    );
                }
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
