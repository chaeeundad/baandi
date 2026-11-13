//! P48B: no steady-state allocation in the packed pure-Rust GEMM. A counting
//! global allocator (this test binary only; one test function, so no other test
//! thread allocates while counting) records every allocation of every thread,
//! pool workers included, while the workspace entry points run on the
//! benchmark shapes after a warm-up. `Packed` resolves to the kernel of the
//! target (x86_64: faer's private x86 GEMM, whose thread-local packing buffer
//! is created on a worker's first call; elsewhere the register-tile kernel,
//! stack storage only), so this test proves the property of the kernel the
//! machine running it would use.
use dft_numerics::{
    pool::WorkerPool,
    real::{GemmSite, PureGemm, RealGemmWorkspace, RealMatrix, RealOp},
};

// Test-only System allocator adapter (as in pair_fused_plan.rs).
#[allow(unsafe_code)]
mod probe {
    use std::{
        alloc::{GlobalAlloc, Layout, System},
        sync::atomic::{AtomicBool, AtomicUsize, Ordering},
    };
    pub static TRACK: AtomicBool = AtomicBool::new(false);
    pub static COUNT: AtomicUsize = AtomicUsize::new(0);
    pub struct Allocator;
    fn record() {
        if TRACK.load(Ordering::Relaxed) {
            COUNT.fetch_add(1, Ordering::Relaxed);
        }
    }
    // SAFETY: every pointer and layout is forwarded unchanged to the System
    // allocator; the counters are atomics and never allocate.
    unsafe impl GlobalAlloc for Allocator {
        unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
            record();
            // SAFETY: the caller supplies a valid layout.
            unsafe { System.alloc(layout) }
        }
        unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
            record();
            // SAFETY: the caller supplies a valid layout.
            unsafe { System.alloc_zeroed(layout) }
        }
        unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
            // SAFETY: pointer and layout come from this forwarding allocator.
            unsafe { System.dealloc(ptr, layout) }
        }
        unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, size: usize) -> *mut u8 {
            record();
            // SAFETY: pointer, layout and size keep the GlobalAlloc contract.
            unsafe { System.realloc(ptr, layout, size) }
        }
    }
}
#[global_allocator]
static ALLOCATOR: probe::Allocator = probe::Allocator;

fn matrix(rows: usize, cols: usize, phase: usize) -> RealMatrix {
    RealMatrix::new(
        rows,
        cols,
        (0..rows * cols)
            .map(|i| ((i * 7 + phase) as f64).sin())
            .collect(),
    )
    .unwrap()
}

struct Case {
    rows: usize,
    k: usize,
    n: usize,
    v: RealMatrix,
    hv: RealMatrix,
    c: RealMatrix,
    values: Vec<f64>,
    x: RealMatrix,
    y: RealMatrix,
    g: RealMatrix,
}
impl Case {
    fn new(rows: usize, k: usize, n: usize) -> Self {
        Self {
            rows,
            k,
            n,
            v: matrix(rows, k, 1),
            hv: matrix(rows, k, 2),
            c: matrix(k, n, 3),
            values: (0..n).map(|j| j as f64).collect(),
            x: RealMatrix::zeros(rows, n).unwrap(),
            y: RealMatrix::zeros(rows, n).unwrap(),
            g: RealMatrix::zeros(k, n).unwrap(),
        }
    }
    fn run(&mut self, ws: &mut RealGemmWorkspace) {
        let (rows, k, n) = (self.rows, self.k, self.n);
        ws.rotate_pair_slices(
            self.v.as_slice(),
            self.hv.as_slice(),
            k,
            self.c.as_slice(),
            n,
            self.x.as_mut_slice(),
            self.y.as_mut_slice(),
        )
        .unwrap();
        ws.rotate_residual_slices(
            self.v.as_slice(),
            self.hv.as_slice(),
            k,
            self.c.as_slice(),
            n,
            &self.values,
            self.x.as_mut_slice(),
            self.y.as_mut_slice(),
        )
        .unwrap();
        ws.multiply_slices_into(
            GemmSite::Gram,
            self.v.as_slice(),
            rows,
            k,
            rows,
            RealOp::Transpose,
            self.hv.as_slice(),
            n,
            rows,
            self.g.as_mut(),
        )
        .unwrap();
        ws.multiply_into(
            self.v.as_ref(),
            RealOp::Identity,
            self.c.as_ref(),
            self.x.as_mut(),
        )
        .unwrap();
        ws.multiply_column_groups(
            self.v.as_slice(),
            rows,
            k,
            rows,
            self.c.as_slice(),
            n,
            k,
            self.y.as_mut_slice(),
        )
        .unwrap();
    }
}

// Release only: debug faer builds allocate inside the nano branch too (3
// allocations per measured phase on this Mac) and the run takes about a minute.
#[test]
#[cfg_attr(
    debug_assertions,
    ignore = "release allocation gate: cargo test --release -p dft-numerics --test pure_gemm_alloc -- --nocapture"
)]
fn packed_gemm_has_no_steady_state_allocation() {
    use std::sync::atomic::Ordering;
    let shapes = [
        (2969, 216, 36),
        (5939, 408, 68),
        (1025, 40, 1),
        (300, 270, 9),
    ];
    let kernels: &[PureGemm] = if cfg!(target_arch = "x86_64") {
        &[PureGemm::Nano, PureGemm::Packed]
    } else {
        &[PureGemm::Nano, PureGemm::Packed, PureGemm::Simd]
    };
    for &kernel in kernels {
        for threads in [1, 4] {
            let pool = WorkerPool::new(threads).unwrap();
            let mut cases: Vec<Case> = shapes.iter().map(|&(r, k, n)| Case::new(r, k, n)).collect();
            let mut ws = RealGemmWorkspace::with_pool(5939, 408, pool.clone()).unwrap();
            ws.set_pure_gemm(kernel);
            // Everything runs inside one installation, as in the SCF: the
            // workspace's own installations are then nested and inject no job
            // into the executor's global queue (whose block allocation every 63
            // outside injections is an executor property, not the kernel's).
            let count = pool.install(|| {
                // Warm-up: every worker's first call (thread-local packing
                // buffer of the x86 library kernel) happens here.
                for _ in 0..8 {
                    for case in &mut cases {
                        case.run(&mut ws);
                    }
                }
                probe::COUNT.store(0, Ordering::SeqCst);
                probe::TRACK.store(true, Ordering::SeqCst);
                for _ in 0..10 {
                    for case in &mut cases {
                        case.run(&mut ws);
                    }
                }
                probe::TRACK.store(false, Ordering::SeqCst);
                probe::COUNT.load(Ordering::SeqCst)
            });
            println!("{kernel:?} threads {threads}: {count} allocations in the measured calls");
            assert_eq!(count, 0, "{kernel:?} allocated at {threads} threads");
        }
    }
}
