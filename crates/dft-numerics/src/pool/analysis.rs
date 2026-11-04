//! Compile-time analysis only. Default builds expose inline zero-sized no-ops.
//! Records are buffered and printed only when the outer SCF scope finishes.
#[cfg(not(feature = "worker-analysis"))]
mod implementation {
    #[inline(always)]
    pub fn copied(_: usize) {}
    #[inline(always)]
    pub fn logical_gemm(_: usize, _: usize, _: usize) {}
    #[inline(always)]
    pub fn executed_gemm(_: usize, _: usize, _: usize) {}
    #[derive(Clone)]
    pub struct Context;
    pub struct Span;
    pub struct Piece;
    impl Drop for Span {
        #[inline(always)]
        fn drop(&mut self) {}
    }
    #[inline(always)]
    pub fn context() -> Context {
        Context
    }
    #[inline(always)]
    pub fn with_context<R>(_: Context, work: impl FnOnce() -> R) -> R {
        work()
    }
    impl Span {
        #[inline(always)]
        pub fn new(_: &'static str) -> Self {
            Self
        }
        #[inline(always)]
        pub fn region(_: &'static str, _: usize) -> Self {
            Self
        }
        #[inline(always)]
        pub fn serial(_: &'static str) -> Self {
            Self
        }
    }
    impl Piece {
        #[inline(always)]
        pub fn new() -> Self {
            Self
        }
    }
    impl Default for Piece {
        fn default() -> Self {
            Self::new()
        }
    }
}
#[cfg(feature = "worker-analysis")]
mod implementation {
    use std::{
        cell::RefCell,
        sync::{
            Arc, Mutex, OnceLock,
            atomic::{AtomicBool, AtomicU64, Ordering},
        },
        time::Instant,
    };
    const SLOTS: usize = 256;
    struct Work {
        busy: [AtomicU64; SLOTS],
        pieces: AtomicU64,
    }
    #[derive(Clone, Default)]
    pub struct Context(Option<Arc<Work>>, Option<Arc<Cost>>);
    #[derive(Default)]
    struct Cost {
        copied: AtomicU64,
        logical: AtomicU64,
        executed: AtomicU64,
    }
    pub fn copied(bytes: usize) {
        if let Some(cost) = context().1 {
            cost.copied.fetch_add(bytes as u64, Ordering::Relaxed);
        }
    }
    pub fn logical_gemm(m: usize, n: usize, k: usize) {
        if let Some(cost) = context().1 {
            cost.logical
                .fetch_add(2 * m as u64 * n as u64 * k as u64, Ordering::Relaxed);
        }
    }
    pub fn executed_gemm(m: usize, n: usize, k: usize) {
        if let Some(cost) = context().1 {
            cost.executed
                .fetch_add(2 * m as u64 * n as u64 * k as u64, Ordering::Relaxed);
        }
    }
    thread_local! { static CURRENT: RefCell<Context> = RefCell::new(Context::default()); }
    static EPOCH: OnceLock<Instant> = OnceLock::new();
    static EVENTS: Mutex<Vec<Event>> = Mutex::new(Vec::new());
    static ACTIVE_SCF: AtomicBool = AtomicBool::new(false);
    fn now() -> u64 {
        EPOCH.get_or_init(Instant::now).elapsed().as_nanos() as u64
    }
    pub fn context() -> Context {
        CURRENT.with(|c| c.borrow().clone())
    }
    struct Restore(Context);
    impl Drop for Restore {
        fn drop(&mut self) {
            CURRENT.with(|c| *c.borrow_mut() = self.0.clone());
        }
    }
    pub fn with_context<R>(context: Context, work: impl FnOnce() -> R) -> R {
        let _restore = Restore(CURRENT.with(|c| c.replace(context)));
        work()
    }
    struct Event {
        name: &'static str,
        start: u64,
        end: u64,
        workers: usize,
        pieces: u64,
        busy: Vec<u64>,
        nested_busy: bool,
        cost: Option<[u64; 3]>,
    }
    pub struct Span {
        name: &'static str,
        start: u64,
        workers: usize,
        work: Option<Arc<Work>>,
        previous: Option<Context>,
        serial_slot: Option<usize>,
        serial_context: Option<Context>,
        cost: Option<(Arc<Cost>, Context)>,
    }
    impl Span {
        pub fn new(name: &'static str) -> Self {
            if name == "scf" {
                assert!(
                    !ACTIVE_SCF.swap(true, Ordering::SeqCst),
                    "analysis supports one active SCF per process"
                );
                EVENTS.lock().unwrap().clear();
            }
            let cost = if name == "iteration" {
                let cost = Arc::new(Cost::default());
                let previous = CURRENT.with(|c| {
                    let mut current = c.borrow_mut();
                    let old = current.clone();
                    current.1 = Some(cost.clone());
                    old
                });
                Some((cost, previous))
            } else {
                None
            };
            Self {
                name,
                start: now(),
                workers: 0,
                work: None,
                previous: None,
                serial_slot: None,
                serial_context: None,
                cost,
            }
        }
        pub fn serial(name: &'static str) -> Self {
            let mut span = Self::new(name);
            span.serial_slot = Some(rayon::current_thread_index().unwrap_or(SLOTS - 1));
            span.serial_context = Some(context());
            span
        }
        pub fn region(name: &'static str, workers: usize) -> Self {
            assert!(workers < SLOTS, "analysis supports at most 255 workers");
            let work = Arc::new(Work {
                busy: std::array::from_fn(|_| AtomicU64::new(0)),
                pieces: AtomicU64::new(0),
            });
            let previous = CURRENT.with(|c| {
                let inherited_cost = c.borrow().1.clone();
                c.replace(Context(Some(work.clone()), inherited_cost))
            });
            Self {
                name,
                start: now(),
                workers,
                work: Some(work),
                previous: Some(previous),
                serial_slot: None,
                serial_context: None,
                cost: None,
            }
        }
    }
    impl Drop for Span {
        fn drop(&mut self) {
            let end = now();
            if let Some(previous) = self.previous.take() {
                CURRENT.with(|c| *c.borrow_mut() = previous);
            }
            let parent = self.serial_context.as_ref().and_then(|c| c.0.as_ref());
            if let (Some(slot), Some(parent)) = (self.serial_slot, parent) {
                parent.busy[slot].fetch_add(end - self.start, Ordering::Relaxed);
            }
            let (pieces, busy) = match &self.work {
                Some(w) => (
                    w.pieces.load(Ordering::Relaxed),
                    w.busy.iter().map(|x| x.load(Ordering::Relaxed)).collect(),
                ),
                None => {
                    if let Some(slot) = self.serial_slot {
                        let mut busy = vec![0; SLOTS];
                        busy[slot] = end - self.start;
                        (0, busy)
                    } else {
                        (0, Vec::new())
                    }
                }
            };
            let cost = self.cost.take().map(|(cost, previous)| {
                CURRENT.with(|c| *c.borrow_mut() = previous);
                [
                    cost.logical.load(Ordering::Relaxed),
                    cost.executed.load(Ordering::Relaxed),
                    cost.copied.load(Ordering::Relaxed),
                ]
            });
            let mut events = EVENTS.lock().unwrap();
            events.push(Event {
                name: self.name,
                start: self.start,
                end,
                workers: self.workers,
                pieces,
                busy,
                nested_busy: parent.is_some(),
                cost,
            });
            if self.name == "scf" {
                ACTIVE_SCF.store(false, Ordering::SeqCst);
                // Diagnostic build only: no output inside a numerical loop.
                for e in events.drain(..) {
                    if let Some([logical, executed, copied]) = e.cost {
                        eprintln!(
                            "P10M_COST,{},{},{},{},{}",
                            e.start, e.end, logical, executed, copied
                        );
                    }
                    eprintln!(
                        "P10L,{},{},{},{},{},{},{}",
                        e.name,
                        e.start,
                        e.end,
                        e.workers,
                        e.pieces,
                        e.busy
                            .iter()
                            .map(u64::to_string)
                            .collect::<Vec<_>>()
                            .join(";"),
                        u8::from(e.nested_busy)
                    );
                }
            }
        }
    }
    pub struct Piece {
        start: u64,
        context: Context,
        slot: usize,
    }
    impl Piece {
        pub fn new() -> Self {
            Self {
                start: now(),
                context: context(),
                slot: rayon::current_thread_index().unwrap_or(SLOTS - 1),
            }
        }
    }
    impl Default for Piece {
        fn default() -> Self {
            Self::new()
        }
    }
    impl Drop for Piece {
        fn drop(&mut self) {
            let elapsed = now() - self.start;
            if let Some(w) = &self.context.0 {
                w.busy[self.slot].fetch_add(elapsed, Ordering::Relaxed);
                w.pieces.fetch_add(1, Ordering::Relaxed);
            }
        }
    }
    #[cfg(test)]
    mod tests {
        use super::*;
        #[test]
        fn iteration_costs_follow_workers_without_counting_nested_regions_twice() {
            let pool = crate::pool::WorkerPool::new(4).unwrap();
            let iteration = Span::new("iteration");
            let cost = context().1.unwrap();
            logical_gemm(13, 7, 65);
            pool.install(|| {
                let _region = Span::region("test_cost", 4);
                crate::pool::WorkerPool::join(
                    || {
                        executed_gemm(13, 7, 32);
                        copied(1024);
                    },
                    || {
                        executed_gemm(13, 7, 33);
                        copied(2048);
                    },
                );
            });
            assert_eq!(cost.logical.load(Ordering::Relaxed), 2 * 13 * 7 * 65);
            assert_eq!(cost.executed.load(Ordering::Relaxed), 2 * 13 * 7 * 65);
            assert_eq!(cost.copied.load(Ordering::Relaxed), 3072);
            drop(iteration);
            assert!(context().1.is_none());
        }

        #[test]
        fn context_follows_stolen_work_and_restores_after_panic() {
            let pool = crate::pool::WorkerPool::new(4).unwrap();
            let span = Span::region("test", 4);
            pool.install(|| {
                crate::pool::WorkerPool::join(
                    || {
                        let _piece = Piece::new();
                        for value in 0..1024_u64 {
                            std::hint::black_box(value.wrapping_mul(17));
                        }
                    },
                    || {
                        let _piece = Piece::new();
                        for value in 0..1024_u64 {
                            std::hint::black_box(value.wrapping_mul(19));
                        }
                    },
                )
            });
            let work = span.work.as_ref().unwrap();
            assert_eq!(work.pieces.load(Ordering::Relaxed), 2);
            assert!(
                work.busy[..4]
                    .iter()
                    .map(|x| x.load(Ordering::Relaxed))
                    .sum::<u64>()
                    > 0
            );
            assert_eq!(work.busy[SLOTS - 1].load(Ordering::Relaxed), 0);
            {
                let _serial = Span::serial("nested_test_serial");
                for value in 0..1024_u64 {
                    std::hint::black_box(value.wrapping_mul(23));
                }
            }
            assert!(work.busy[SLOTS - 1].load(Ordering::Relaxed) > 0);
            assert_eq!(work.pieces.load(Ordering::Relaxed), 2);
            assert!(
                EVENTS
                    .lock()
                    .unwrap()
                    .iter()
                    .any(|e| e.name == "nested_test_serial" && e.nested_busy)
            );
            let result =
                std::panic::catch_unwind(|| with_context(Context::default(), || panic!("test")));
            assert!(result.is_err());
            assert!(Arc::ptr_eq(context().0.as_ref().unwrap(), work));
            drop(span);
            assert!(context().0.is_none());
        }
    }
}
pub use implementation::*;
