//! An owned executor shared by an SCF's sequential and nested parallel phases.
//! No global Rayon configuration; the one-thread path never creates a worker.
use std::sync::Arc;

// Analysis builds opt in with RUSTFLAGS='--cfg feature="worker-analysis"'.
// Kept source-only because this ticket does not own Cargo manifests.
#[allow(unexpected_cfgs)]
#[doc(hidden)]
pub mod analysis;

pub struct WorkerPool {
    pool: Option<rayon::ThreadPool>,
    created: Arc<std::sync::atomic::AtomicUsize>,
}
impl WorkerPool {
    pub fn new(threads: usize) -> Result<Arc<Self>, rayon::ThreadPoolBuildError> {
        let created = Arc::new(std::sync::atomic::AtomicUsize::new(0));
        let pool = if threads <= 1 {
            None
        } else {
            let builder = rayon::ThreadPoolBuilder::new().num_threads(threads);
            let builder = {
                let created = created.clone();
                builder.start_handler(move |_| {
                    created.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
                })
            };
            let pool = builder.build()?;
            pool.broadcast(|_| ());
            Some(pool)
        };
        Ok(Arc::new(Self { pool, created }))
    }

    /// Workers actually started at construction. Only startup writes this
    /// diagnostic counter; submitting a job does not touch it.
    pub fn workers_created(&self) -> usize {
        self.created.load(std::sync::atomic::Ordering::SeqCst)
    }

    /// Install only at a parallel region boundary. Nested joins use this same
    /// executor, so outer k and inner band/tile tasks share one thread budget.
    pub fn install<R: Send>(&self, work: impl FnOnce() -> R + Send) -> R {
        let _span = analysis::Span::new("pool_install");
        let context = analysis::context();
        let work = || analysis::with_context(context, work);
        match &self.pool {
            Some(pool) => pool.install(work),
            None => work(),
        }
    }

    /// Index in the currently executing pool, for bounded reusable scratch.
    pub fn current_worker_index() -> Option<usize> {
        rayon::current_thread_index()
    }

    /// Call inside `install`. Both branches finish before returning or unwinding.
    pub fn join<A: Send, B: Send>(
        left: impl FnOnce() -> A + Send,
        right: impl FnOnce() -> B + Send,
    ) -> (A, B) {
        let context = analysis::context();
        let other = context.clone();
        let left = || analysis::with_context(context, left);
        let right = || analysis::with_context(other, right);
        if rayon::current_thread_index().is_some() {
            rayon::join(left, right)
        } else {
            // Serial installations must never initialize Rayon's global pool.
            let left = std::panic::catch_unwind(std::panic::AssertUnwindSafe(left));
            let right = std::panic::catch_unwind(std::panic::AssertUnwindSafe(right));
            match (left, right) {
                (Ok(a), Ok(b)) => (a, b),
                (Err(panic), _) | (_, Err(panic)) => std::panic::resume_unwind(panic),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::Ordering;

    fn empty_region(tasks: usize) {
        if tasks <= 1 {
            std::hint::black_box(());
        } else {
            WorkerPool::join(
                || empty_region(tasks / 2),
                || empty_region(tasks - tasks / 2),
            );
        }
    }

    /// cargo test --release -p dft-numerics dispatch_latency -- --ignored --nocapture
    #[test]
    #[ignore]
    fn dispatch_latency() {
        use std::time::Instant;
        println!("threads,ingress_us,external_region_us,resident_region_us");
        for threads in [1, 2, 4, 8, 12, 16] {
            let pool = WorkerPool::new(threads).unwrap();
            for _ in 0..100 {
                pool.install(|| empty_region(threads));
            }
            let mut times = [0.0; 3];
            for (mode, time) in times.iter_mut().enumerate() {
                let mut samples = Vec::new();
                for _ in 0..9 {
                    let batch = || {
                        let start = Instant::now();
                        for _ in 0..1000 {
                            match mode {
                                0 => pool.install(|| std::hint::black_box(())),
                                1 => pool.install(|| empty_region(threads)),
                                _ => empty_region(threads),
                            }
                        }
                        start.elapsed().as_secs_f64() * 1e6 / 1000.0
                    };
                    samples.push(if mode == 2 {
                        pool.install(batch)
                    } else {
                        batch()
                    });
                }
                samples.sort_by(f64::total_cmp);
                *time = samples[4];
            }
            assert_eq!(
                pool.workers_created(),
                if threads == 1 { 0 } else { threads }
            );
            println!("{threads},{:.3},{:.3},{:.3}", times[0], times[1], times[2]);
        }
    }

    #[test]
    fn workers_created_once_across_nested_phases() {
        let serial = WorkerPool::new(1).unwrap();
        assert_eq!(serial.install(|| 17), 17);
        assert_eq!(serial.created.load(Ordering::SeqCst), 0);
        let pool = WorkerPool::new(4).unwrap();
        // Broadcast waits for every worker, making the creation counter exact.
        pool.pool.as_ref().unwrap().broadcast(|_| ());
        assert_eq!(pool.created.load(Ordering::SeqCst), 4);
        for _ in 0..100 {
            let (a, b) = pool.install(|| {
                WorkerPool::join(
                    || 13,
                    || {
                        let (a, b) = WorkerPool::join(|| 3, || 7);
                        a + b
                    },
                )
            });
            assert_eq!((a, b), (13, 10));
        }
        assert_eq!(pool.created.load(Ordering::SeqCst), 4);
        let panic = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            pool.install(|| WorkerPool::join(|| panic!("test"), || 0))
        }));
        assert!(panic.is_err());
        assert_eq!(pool.install(|| 11), 11);
        assert_eq!(pool.created.load(Ordering::SeqCst), 4);
    }
}
