use super::ReducedPlacement;
use dft_mpi::{Collective, MpiError};
use std::{cell::Cell, time::Instant};

/// Exclusive stage clocks, seconds. Total includes setup/allocation and glue.
pub const STAGES: [&str; 8] = [
    "apply_h",
    "border",
    "factor",
    "reduced_eigen",
    "corrections",
    "repair",
    "restart",
    "return",
];
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Communication {
    /// Logical collective invocations and input payload bytes PER RANK, not
    /// wire traffic. Operator-internal collectives are included, tagged by stage.
    pub reductions: u64,
    pub reduction_bytes: u64,
    pub broadcasts: u64,
    pub broadcast_bytes: u64,
    pub gathers: u64,
    pub gather_bytes: u64,
    pub alltoalls: u64,
    pub alltoall_bytes: u64,
}
impl Communication {
    pub(crate) fn difference(self, a: Self) -> Self {
        Self {
            reductions: self.reductions - a.reductions,
            reduction_bytes: self.reduction_bytes - a.reduction_bytes,
            broadcasts: self.broadcasts - a.broadcasts,
            broadcast_bytes: self.broadcast_bytes - a.broadcast_bytes,
            gathers: self.gathers - a.gathers,
            gather_bytes: self.gather_bytes - a.gather_bytes,
            alltoalls: self.alltoalls - a.alltoalls,
            alltoall_bytes: self.alltoall_bytes - a.alltoall_bytes,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Gemm {
    pub stage: &'static str,
    pub m: usize,
    pub n: usize,
    pub k: usize,
    pub flops: u64,
}
#[derive(Clone, Debug)]
pub struct Iteration {
    /// Incoming band indices, zero-based. Flags are recomputed, never locked.
    pub active_bands: Vec<usize>,
    pub basis_before: usize,
    pub h_calls: usize,
    pub h_columns: usize,
    pub reduced_order: usize,
    /// Logical reduced problems, identical on all ranks.
    pub reduced_solves: usize,
    /// Actual backend calls on this rank (zero on non-root in Root mode).
    pub backend_solves: usize,
    pub restarts: usize,
    pub repairs: usize,
    pub dropped_corrections: usize,
    pub unconverged: usize,
    pub communication: Communication,
    pub gemms: Vec<Gemm>,
    pub seconds: [f64; 8],
}
impl Iteration {
    fn new(cap: usize) -> Self {
        Self {
            active_bands: Vec::with_capacity(cap),
            basis_before: 0,
            h_calls: 0,
            h_columns: 0,
            reduced_order: 0,
            reduced_solves: 0,
            backend_solves: 0,
            restarts: 0,
            repairs: 0,
            dropped_corrections: 0,
            unconverged: 0,
            communication: Communication::default(),
            gemms: Vec::with_capacity(64 + 24 * cap.div_ceil(64)),
            seconds: [0.0; 8],
        }
    }
    pub(crate) fn gemm(&mut self, stage: &'static str, m: usize, n: usize, k: usize) {
        self.gemms.push(Gemm {
            stage,
            m,
            n,
            k,
            flops: 2 * m as u64 * n as u64 * k as u64,
        });
    }
}
#[derive(Debug)]
pub struct Record {
    pub apply_h: Option<ApplyH>,
    pub setup: Iteration,
    pub iterations: Vec<Iteration>,
    pub finish: Iteration,
    pub communication: Communication,
    pub h_calls: usize,
    pub h_columns: usize,
    pub reduced_solves: usize,
    pub backend_solves: usize,
    pub restarts: usize,
    pub gemm_flops: u64,
    pub profile: bool,
    pub placement: ReducedPlacement,
    pub total_seconds: f64,
}
impl Record {
    pub(crate) fn new(it: usize, cap: usize, profile: bool, placement: ReducedPlacement) -> Self {
        Self {
            apply_h: None,
            setup: Iteration::new(cap),
            iterations: (0..it).map(|_| Iteration::new(cap)).collect(),
            finish: Iteration::new(cap),
            communication: Communication::default(),
            h_calls: 0,
            h_columns: 0,
            reduced_solves: 0,
            backend_solves: 0,
            restarts: 0,
            gemm_flops: 0,
            profile,
            placement,
            total_seconds: 0.0,
        }
    }
    pub(crate) fn summarize(&mut self) {
        for r in std::iter::once(&self.setup)
            .chain(&self.iterations)
            .chain(std::iter::once(&self.finish))
        {
            self.h_calls += r.h_calls;
            self.h_columns += r.h_columns;
            self.reduced_solves += r.reduced_solves;
            self.backend_solves += r.backend_solves;
            self.restarts += r.restarts;
            self.gemm_flops += r.gemms.iter().map(|g| g.flops).sum::<u64>();
        }
    }
}
pub(crate) fn clock(on: bool) -> Option<Instant> {
    on.then(Instant::now)
}
pub(crate) fn elapsed(t: Option<Instant>) -> f64 {
    t.map_or(0.0, |t| t.elapsed().as_secs_f64())
}
pub(crate) struct Counted<'a, C> {
    inner: &'a C,
    counts: Cell<Communication>,
}
impl<'a, C> Counted<'a, C> {
    pub fn new(inner: &'a C) -> Self {
        Self {
            inner,
            counts: Cell::new(Communication::default()),
        }
    }
    pub fn snapshot(&self) -> Communication {
        self.counts.get()
    }
}
impl<C: Collective> Collective for Counted<'_, C> {
    fn rank(&self) -> usize {
        self.inner.rank()
    }
    fn size(&self) -> usize {
        self.inner.size()
    }
    fn barrier(&self) -> Result<(), MpiError> {
        self.inner.barrier()
    }
    fn tree_sum_f64(&self, a: &[f64], b: &mut [f64]) -> Result<(), MpiError> {
        let mut s = self.counts.get();
        s.reductions += 1;
        s.reduction_bytes += 8 * a.len() as u64;
        self.counts.set(s);
        self.inner.tree_sum_f64(a, b)
    }
    fn bcast_f64(&self, a: &mut [f64], root: usize) -> Result<(), MpiError> {
        let mut s = self.counts.get();
        s.broadcasts += 1;
        s.broadcast_bytes += 8 * a.len() as u64;
        self.counts.set(s);
        self.inner.bcast_f64(a, root)
    }
    fn allgatherv_f64(&self, a: &[f64], counts: &[usize], b: &mut [f64]) -> Result<(), MpiError> {
        let mut s = self.counts.get();
        s.gathers += 1;
        s.gather_bytes += 8 * a.len() as u64;
        self.counts.set(s);
        self.inner.allgatherv_f64(a, counts, b)
    }
    fn alltoall_f64(&self, a: &[f64], block: usize, b: &mut [f64]) -> Result<(), MpiError> {
        let mut s = self.counts.get();
        s.alltoalls += 1;
        s.alltoall_bytes += 8 * a.len() as u64;
        self.counts.set(s);
        self.inner.alltoall_f64(a, block, b)
    }
}

/// Rank-local H-application totals for one reference solve. No numerical code
/// reads these counters. FFT entries follow `real::fft_profile::NAMES`.
#[derive(Clone, Debug, Default)]
pub struct ApplyH {
    pub seconds: [f64; 9],
    pub fft: dft_numerics::real::fft_profile::Snapshot,
    pub columns_per_call: Vec<usize>,
    pub pairs: usize,
    pub odd_tails: usize,
    pub wide_calls: usize,
    pub direct_calls: usize,
}
pub const APPLY_H_STAGES: [&str; 9] = [
    "total",
    "mirror",
    "pack",
    "unpack",
    "kinetic",
    "projection",
    "projection_reduce",
    "coupling",
    "back_and_add",
];
thread_local! {
    static APPLY_H: std::cell::RefCell<Option<ApplyH>> = const { std::cell::RefCell::new(None) };
    static APPLY_H_SOLVES: std::cell::RefCell<Vec<ApplyH>> = const { std::cell::RefCell::new(Vec::new()) };
}
/// Scope is owned by the solve, including early-error cleanup. Storage is
/// reserved before entering the Davidson loop; no per-H-call allocation.
pub(crate) struct ApplyHScope(bool);
impl ApplyHScope {
    pub fn new(enabled: bool, max_calls: usize) -> Self {
        if enabled {
            APPLY_H.with(|s| {
                assert!(s.borrow().is_none(), "nested reference H profiling");
                *s.borrow_mut() = Some(ApplyH {
                    columns_per_call: Vec::with_capacity(max_calls),
                    ..ApplyH::default()
                });
            });
        }
        Self(enabled)
    }
    pub fn finish(mut self) -> Option<ApplyH> {
        if !self.0 {
            return None;
        }
        self.0 = false;
        let record = APPLY_H
            .with(|s| s.borrow_mut().take())
            .expect("active H profile");
        APPLY_H_SOLVES.with(|s| s.borrow_mut().push(record.clone()));
        Some(record)
    }
}
impl Drop for ApplyHScope {
    fn drop(&mut self) {
        if self.0 {
            APPLY_H.with(|s| {
                s.borrow_mut().take();
            });
        }
    }
}
/// Consume completed records at the reporting boundary, never in a hot loop.
pub fn take_apply_h_solves() -> Vec<ApplyH> {
    APPLY_H_SOLVES.with(|s| std::mem::take(&mut *s.borrow_mut()))
}
/// True only inside a profiled reference solve on this rank/thread.
pub fn apply_h_active() -> bool {
    APPLY_H.with(|s| s.borrow().is_some())
}
pub fn record_apply_h(
    seconds: [f64; 9],
    fft: dft_numerics::real::fft_profile::Snapshot,
    columns: usize,
    wide: bool,
    nonlocal: bool,
) {
    APPLY_H.with(|s| {
        if let Some(s) = s.borrow_mut().as_mut() {
            for (a, b) in s.seconds.iter_mut().zip(seconds) {
                *a += b;
            }
            s.fft.add(fft);
            // Capacity is max_iterations + setup; exceeding it is a contract bug.
            assert!(s.columns_per_call.len() < s.columns_per_call.capacity());
            s.columns_per_call.push(columns);
            s.pairs += columns.div_ceil(2);
            s.odd_tails += columns % 2;
            if nonlocal {
                if wide {
                    s.wide_calls += 2;
                } else {
                    s.direct_calls += 2;
                }
            }
        }
    });
}
