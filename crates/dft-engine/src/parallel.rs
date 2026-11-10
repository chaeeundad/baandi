//! Process-wide parallel communicator (ADR 0041). MPI is itself process-global
//! (one world per process), so the engine reads a registered communicator instead of
//! threading one through every option struct. Unset or one rank means the
//! threads-only path, which is bit-identical to a build without MPI.
use dft_mpi::{Collective, MpiError};
pub use dft_solvers::distributed_davidson::stage_seconds;
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
thread_local! {
    // Test-only per-thread override of BAANDI_MPI_M2 (simulated ranks run on their own threads).
    static M2_OVERRIDE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}
static DISTRIBUTED_ROWS_SOLVES: AtomicUsize = AtomicUsize::new(0);
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

/// Stick-row distributed Davidson (MPI stage M2c) requested: `BAANDI_MPI_M2=on`
/// (default off) unless overridden for the calling thread by `set_m2_override`.
pub fn m2_enabled() -> bool {
    M2_OVERRIDE
        .with(|o| o.get())
        .unwrap_or_else(|| std::env::var("BAANDI_MPI_M2").is_ok_and(|v| v == "on"))
}

/// (counter name, value) list and (stage name, seconds) list.
pub type M2Stats = (Vec<(String, usize)>, Vec<(String, f64)>);

/// Proof-of-branch counters and stage seconds of the M2 path in this process:
/// (counters, stages). Counters: `rows_solves` (M2c eigensolves), the distributed
/// operator's `batched_applies` / `unbatched_applies` / `applied_columns` /
/// `band_density_builds` / `nonlocal_energy_evaluations`, and the engine's M2d
/// `density_builds` / `xc_evaluations` / `nonlocal_energy_audits` /
/// `screen_integrals`. Stages: operator `apply_*`, `band_densities`,
/// `nonlocal_energy`, `gather_planes` and engine `density`, `xc`,
/// `nonlocal_energy`, `screen_integrals` (engine stages include the operator
/// work they call; operator stages are named `op.*`).
pub fn m2_stats() -> M2Stats {
    let (op_counts, op_stages) = dft_physics::distributed_real_gamma::m2_counters();
    let (engine_counts, engine_stages) = crate::m2::stats();
    let mut counts = vec![("rows_solves".to_string(), distributed_rows_solves())];
    counts.extend(op_counts.into_iter().map(|(n, v)| (n.to_string(), v)));
    counts.extend(engine_counts.into_iter().map(|(n, v)| (n.to_string(), v)));
    let mut stages: Vec<(String, f64)> = engine_stages
        .into_iter()
        .map(|(n, v)| (n.to_string(), v))
        .collect();
    stages.extend(op_stages.into_iter().map(|(n, v)| (format!("op.{n}"), v)));
    (counts, stages)
}

/// Override `BAANDI_MPI_M2` for the calling thread (tests); `None` restores the environment.
pub fn set_m2_override(value: Option<bool>) {
    M2_OVERRIDE.with(|o| o.set(value));
}

/// Number of eigensolves that ran on the M2c stick-row path in this process.
pub fn distributed_rows_solves() -> usize {
    DISTRIBUTED_ROWS_SOLVES.load(Ordering::Relaxed)
}

pub(crate) fn count_distributed_rows_solve() {
    DISTRIBUTED_ROWS_SOLVES.fetch_add(1, Ordering::Relaxed);
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
