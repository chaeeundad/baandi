//! Full-diagnostic, rank-local work attribution. No serialization or numerical
//! decisions here. Bytes are logical traffic, never measured DRAM/wire bytes.
use std::{cell::RefCell, time::Instant};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(usize)]
pub(crate) enum Stage {
    Setup,
    ApplyH,
    ApplyCheck,
    Append,
    BorderS,
    BorderH,
    Factor,
    Transform,
    Eigen,
    Coefficients,
    Corrections,
    Precondition,
    Normalize,
    Repair,
    RepairPanel,
    Restart,
    Return,
    Decisions,
    Glue,
}
pub const NAMES: [&str; 19] = [
    "setup",
    "apply_h",
    "apply_check",
    "append",
    "border_s",
    "border_h",
    "factor",
    "transform",
    "reduced_eigen",
    "coefficients",
    "corrections",
    "precondition",
    "normalize",
    "repair",
    "repair_panel",
    "restart",
    "return",
    "decisions",
    "glue",
];
#[derive(Clone, Copy, Debug, Default)]
pub struct Counters {
    pub seconds: f64,
    pub calls: u64,
    pub gemm_calls: u64,
    pub gemm_flops: u64,
    pub gemm_operand_bytes: u64,
    /// Scalar nominal FLOPs. sqrt/hypot are counted separately.
    pub scalar_flops: u64,
    pub sqrt_calls: u64,
    pub hypot_calls: u64,
    /// Destination bytes of explicit copies, excluding arithmetic writes.
    pub copied_bytes: u64,
    pub zeroed_bytes: u64,
    /// Logical non-GEMM reads + writes, including both sides of copies.
    pub streamed_bytes: u64,
    /// barrier, tree_sum, broadcast, allgatherv, alltoall (logical input).
    pub collective_calls: [u64; 5],
    pub collective_bytes: [u64; 5],
    pub collective_seconds: f64,
    /// Panels whose rejected scalar factor work cannot be counted exactly
    /// without instrumenting the shared adapter. Nominal full panel is charged.
    pub rejected_factor_panels: u64,
}
#[derive(Clone, Copy, Debug)]
pub struct Product {
    pub stage: usize,
    pub site: &'static str,
    pub m: usize,
    pub n: usize,
    pub k: usize,
}
#[derive(Clone, Debug)]
pub struct Work {
    pub stages: [Counters; NAMES.len()],
    /// Logical adapter invocations, not backend microtiles/BLAS-internal calls.
    pub gemms: Vec<Product>,
}
struct Active {
    work: Work,
    stage: Stage,
    since: Instant,
}
impl Active {
    fn switch(&mut self, stage: Stage) -> Stage {
        let now = Instant::now();
        self.work.stages[self.stage as usize].seconds +=
            now.duration_since(self.since).as_secs_f64();
        self.since = now;
        std::mem::replace(&mut self.stage, stage)
    }
}
thread_local! {
    static ACTIVE: RefCell<Option<Active>> = const { RefCell::new(None) };
    static COMPLETED: RefCell<Vec<Work>> = const { RefCell::new(Vec::new()) };
}
/// Drained at the CLI reporting boundary. Failed solves do not publish partial work.
pub fn take_solves() -> Vec<Work> {
    COMPLETED.with(|s| std::mem::take(&mut *s.borrow_mut()))
}
pub(crate) struct Scope(bool);
impl Scope {
    pub fn new(on: bool, iterations: usize, cap: usize) -> Self {
        if on {
            // Each iteration can repair and restart, including repeated panel
            // factorizations. This conservative bound avoids hot-loop growth.
            let capacity = (iterations + 2) * (96 + 32 * cap.div_ceil(64));
            let work = Work {
                stages: [Counters::default(); NAMES.len()],
                gemms: Vec::with_capacity(capacity),
            };
            ACTIVE.with(|a| {
                assert!(a.borrow().is_none(), "nested Davidson work scope");
                *a.borrow_mut() = Some(Active {
                    work,
                    stage: Stage::Setup,
                    since: Instant::now(),
                });
            });
            add(true, |s| s.calls += 1);
        }
        Self(on)
    }
    pub fn finish(mut self) {
        if self.0 {
            self.0 = false;
            let mut a = ACTIVE.with(|a| a.borrow_mut().take().expect("work scope"));
            a.switch(Stage::Setup);
            COMPLETED.with(|s| s.borrow_mut().push(a.work));
        }
    }
}
impl Drop for Scope {
    fn drop(&mut self) {
        if self.0 {
            ACTIVE.with(|a| {
                a.borrow_mut().take();
            });
        }
    }
}
/// Nested spans pause their parent, so stage seconds are exclusive. An off
/// span has no clock reads, TLS access, allocations, or counter updates.
pub(crate) struct Span(Option<Stage>);
pub(crate) fn span(on: bool, stage: Stage) -> Span {
    if !on {
        return Span(None);
    }
    Span(ACTIVE.with(|a| {
        a.borrow_mut().as_mut().map(|a| {
            a.work.stages[stage as usize].calls += 1;
            a.switch(stage)
        })
    }))
}
impl Drop for Span {
    fn drop(&mut self) {
        if let Some(previous) = self.0 {
            ACTIVE.with(|a| {
                if let Some(a) = a.borrow_mut().as_mut() {
                    a.switch(previous);
                }
            });
        }
    }
}
#[inline]
pub(crate) fn add(on: bool, f: impl FnOnce(&mut Counters)) {
    if on {
        ACTIVE.with(|a| {
            if let Some(a) = a.borrow_mut().as_mut() {
                f(&mut a.work.stages[a.stage as usize]);
            }
        });
    }
}
pub(crate) fn copy(on: bool, elements: usize) {
    add(on, |s| {
        s.copied_bytes += 8 * elements as u64;
        s.streamed_bytes += 16 * elements as u64;
    });
}
pub(crate) fn gemm(on: bool, site: &'static str, m: usize, n: usize, k: usize) {
    if on {
        ACTIVE.with(|a| {
            if let Some(a) = a.borrow_mut().as_mut() {
                let s = &mut a.work.stages[a.stage as usize];
                s.gemm_calls += 1;
                s.gemm_flops += 2 * m as u64 * n as u64 * k as u64;
                s.gemm_operand_bytes +=
                    8 * (m as u64 * k as u64 + k as u64 * n as u64 + m as u64 * n as u64);
                assert!(
                    a.work.gemms.len() < a.work.gemms.capacity(),
                    "Davidson work trace bound"
                );
                a.work.gemms.push(Product {
                    stage: a.stage as usize,
                    site,
                    m,
                    n,
                    k,
                });
            }
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ref_davidson::accounting::Counted;
    use dft_mpi::Collective;

    #[test]
    fn all_collective_kinds_are_forwarded_and_attributed_once() {
        dft_mpi::simulate(3, |comm| {
            take_solves();
            let scope = Scope::new(true, 1, 2);
            let c = Counted::new(comm, true);
            {
                let _stage = span(true, Stage::Repair);
                c.barrier().unwrap();
                let mut sum = [0.];
                c.tree_sum_f64(&[comm.rank() as f64], &mut sum).unwrap();
                assert_eq!(sum, [3.]);
                let mut broadcast = [comm.rank() as f64; 2];
                c.bcast_f64(&mut broadcast, 1).unwrap();
                assert_eq!(broadcast, [1.; 2]);
                let mut gathered = [0.; 3];
                c.allgatherv_f64(&[comm.rank() as f64], &[1; 3], &mut gathered)
                    .unwrap();
                assert_eq!(gathered, [0., 1., 2.]);
                let mut exchanged = [0.; 3];
                c.alltoall_f64(&[comm.rank() as f64; 3], 1, &mut exchanged)
                    .unwrap();
                assert_eq!(exchanged, [0., 1., 2.]);
            }
            scope.finish();
            let record = take_solves().pop().unwrap();
            let repair = &record.stages[Stage::Repair as usize];
            assert_eq!(repair.collective_calls, [1; 5]);
            assert_eq!(repair.collective_bytes, [0, 8, 16, 8, 24]);
            assert!(repair.collective_seconds > 0.);
            assert!(repair.collective_seconds <= repair.seconds);
            assert_eq!(
                record.stages[Stage::Setup as usize].collective_calls,
                [0; 5]
            );
        });
    }
}

/// Initialized backing storage observed at allocation/ownership boundaries.
/// Component maxima are independent, not a simultaneous process peak. Backend
/// GEMM/EVD and private LeanScratch allocations are deliberately not inferred.
#[derive(Clone, Copy, Debug, Default)]
pub struct WorkspaceMemory {
    pub rows: [usize; 6], // V, HV, X, HX, E, HE at workspace construction
    pub reduced: usize,
    pub vectors: usize,
    pub returned: usize,
    pub cached: usize,
}
thread_local! {
    static MEMORY: std::cell::Cell<Option<WorkspaceMemory>> = const { std::cell::Cell::new(None) };
}
/// Run-local, thread-bound census; nested scopes restore their caller on drop.
pub struct MemoryScope {
    previous: Option<WorkspaceMemory>,
    _thread: std::marker::PhantomData<std::rc::Rc<()>>,
}
impl MemoryScope {
    pub fn begin() -> Self {
        Self {
            previous: MEMORY.with(|m| m.replace(Some(WorkspaceMemory::default()))),
            _thread: std::marker::PhantomData,
        }
    }
    pub fn snapshot(&self) -> WorkspaceMemory {
        MEMORY.with(|m| m.get().expect("active workspace memory scope"))
    }
}
impl Drop for MemoryScope {
    fn drop(&mut self) {
        MEMORY.with(|m| m.set(self.previous));
    }
}
pub(crate) fn memory_observe(value: WorkspaceMemory) {
    MEMORY.with(|m| {
        if let Some(mut peak) = m.get() {
            for (p, v) in peak.rows.iter_mut().zip(value.rows) {
                *p = (*p).max(v);
            }
            peak.reduced = peak.reduced.max(value.reduced);
            peak.vectors = peak.vectors.max(value.vectors);
            peak.returned = peak.returned.max(value.returned);
            peak.cached = peak.cached.max(value.cached);
            m.set(Some(peak));
        }
    });
}

/// Release a completed reference run's private row cache before delivery,
/// including the complex (k-point) Davidson workspace kept by GK3.
pub fn release_row_cache() {
    super::workspace::release_row_cache();
    crate::ref_davidson_complex::release_workspace();
}
