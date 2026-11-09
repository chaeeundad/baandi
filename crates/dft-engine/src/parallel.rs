//! Process-wide parallel communicator (ADR 0041). MPI is itself process-global
//! (one world per process), so the engine reads a registered communicator instead of
//! threading one through every option struct. Unset or one rank means the
//! threads-only path, which is bit-identical to a build without MPI.
use dft_mpi::{Collective, MpiError};
use std::{
    cell::RefCell,
    sync::{
        Arc,
        atomic::{AtomicU64, AtomicUsize, Ordering},
    },
};

type Shared = Arc<dyn Collective + Send + Sync>;
thread_local! {
    // Per thread: real MPI registers it on the calling (main) thread; the rank
    // simulator registers one per simulated rank.
    static COMM: RefCell<Option<Shared>> = const { RefCell::new(None) };
}
static DISTRIBUTED_SOLVES: AtomicUsize = AtomicUsize::new(0);
// Per-operation call counts and wall nanoseconds spent inside collectives
// (allgatherv, tree_sum, bcast, alltoall, barrier), process-wide.
static CALLS: [AtomicUsize; 5] = [const { AtomicUsize::new(0) }; 5];
static NANOS: [AtomicU64; 5] = [const { AtomicU64::new(0) }; 5];
const OPS: [&str; 5] = ["allgatherv", "tree_sum", "bcast", "alltoall", "barrier"];

/// (operation, calls, seconds) for every collective used so far.
pub fn collective_stats() -> Vec<(&'static str, usize, f64)> {
    (0..5)
        .map(|i| {
            (
                OPS[i],
                CALLS[i].load(Ordering::Relaxed),
                NANOS[i].load(Ordering::Relaxed) as f64 * 1e-9,
            )
        })
        .collect()
}

fn timed<R>(op: usize, f: impl FnOnce() -> R) -> R {
    let start = std::time::Instant::now();
    let r = f();
    CALLS[op].fetch_add(1, Ordering::Relaxed);
    NANOS[op].fetch_add(start.elapsed().as_nanos() as u64, Ordering::Relaxed);
    r
}

/// Register (or clear) the communicator of the calling thread.
pub fn set_communicator(comm: Option<Shared>) {
    COMM.with(|c| *c.borrow_mut() = comm);
}

/// Number of eigensolves that ran on the distributed path in this process.
pub fn distributed_solves() -> usize {
    DISTRIBUTED_SOLVES.load(Ordering::Relaxed)
}

pub(crate) fn count_distributed_solve() {
    DISTRIBUTED_SOLVES.fetch_add(1, Ordering::Relaxed);
}

/// A cloneable handle that forwards every collective to the registered communicator.
#[derive(Clone)]
pub struct SharedComm(Shared);

impl SharedComm {
    /// The registered communicator when it has two or more ranks.
    pub fn active() -> Option<Self> {
        COMM.with(|c| {
            c.borrow()
                .as_ref()
                .filter(|c| c.size() > 1)
                .map(|c| Self(c.clone()))
        })
    }
}

impl Collective for SharedComm {
    fn rank(&self) -> usize {
        self.0.rank()
    }
    fn size(&self) -> usize {
        self.0.size()
    }
    fn barrier(&self) -> Result<(), MpiError> {
        timed(4, || self.0.barrier())
    }
    fn allgatherv_f64(
        &self,
        local: &[f64],
        counts: &[usize],
        out: &mut [f64],
    ) -> Result<(), MpiError> {
        timed(0, || self.0.allgatherv_f64(local, counts, out))
    }
    fn tree_sum_f64(&self, local: &[f64], out: &mut [f64]) -> Result<(), MpiError> {
        timed(1, || self.0.tree_sum_f64(local, out))
    }
    fn bcast_f64(&self, data: &mut [f64], root: usize) -> Result<(), MpiError> {
        timed(2, || self.0.bcast_f64(data, root))
    }
    fn alltoall_f64(&self, send: &[f64], block: usize, recv: &mut [f64]) -> Result<(), MpiError> {
        timed(3, || self.0.alltoall_f64(send, block, recv))
    }
}
