//! Process-level MPI state of the CLI (ADR 0041). Rank 0 writes the canonical
//! outputs; other ranks write the same files under `<output>/rank<r>/`, which is
//! also the evidence that ranks agree.
use std::sync::atomic::{AtomicUsize, Ordering};

static RANK: AtomicUsize = AtomicUsize::new(0);
static SIZE: AtomicUsize = AtomicUsize::new(1);

pub fn set(rank: usize, size: usize) {
    RANK.store(rank, Ordering::Relaxed);
    SIZE.store(size, Ordering::Relaxed);
}
pub fn rank() -> usize {
    RANK.load(Ordering::Relaxed)
}
pub fn size() -> usize {
    SIZE.load(Ordering::Relaxed)
}
