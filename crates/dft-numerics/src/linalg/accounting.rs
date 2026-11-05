//! Source-only logical GEMM measurement. Counts use the historical 2*m*n*k
//! convention (multiply by four for complex real-operation equivalents).
//! Storage is reserved at workspace construction; no per-call allocation/I/O.
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{LazyLock, Mutex};
use std::time::Instant;
static NEXT_SOLVE: AtomicUsize = AtomicUsize::new(0);
pub fn solve_id() -> usize {
    if ENABLED {
        NEXT_SOLVE.fetch_add(1, Ordering::Relaxed)
    } else {
        0
    }
}

pub const ENABLED: bool = option_env!("P14B_GEMM_ACCOUNTING").is_some();
const CAPACITY: usize = 65536;
static RECORDS: LazyLock<Mutex<Vec<Record>>> =
    LazyLock::new(|| Mutex::new(Vec::with_capacity(CAPACITY)));

pub struct Record {
    pub solve: usize,
    pub step: usize,
    pub minimum: u64,
    pub site: &'static str,
    pub m: usize,
    pub n: usize,
    pub k: usize,
    pub calls: u64,
    pub flops: u64,
    pub nanoseconds: u128,
}

pub fn prepare() {
    if ENABLED {
        LazyLock::force(&RECORDS);
    }
}

pub struct Span(Option<(Instant, Record)>);
impl Span {
    #[inline]
    pub fn new(
        solve: usize,
        step: usize,
        site: &'static str,
        shape: [usize; 3],
        calls: u64,
        minimum: u64,
    ) -> Self {
        let [m, n, k] = shape;
        Self(ENABLED.then(|| {
            (
                Instant::now(),
                Record {
                    solve,
                    step,
                    minimum,
                    site,
                    m,
                    n,
                    k,
                    calls,
                    flops: calls * 2 * m as u64 * n as u64 * k as u64,
                    nanoseconds: 0,
                },
            )
        }))
    }
}
impl Drop for Span {
    fn drop(&mut self) {
        if let Some((start, mut record)) = self.0.take() {
            record.nanoseconds = start.elapsed().as_nanos();
            let mut records = RECORDS.lock().expect("GEMM accounting lock");
            if let Some(old) = records.iter_mut().find(|r| {
                (r.solve, r.step, r.site, r.m, r.n, r.k)
                    == (
                        record.solve,
                        record.step,
                        record.site,
                        record.m,
                        record.n,
                        record.k,
                    )
            }) {
                old.minimum += record.minimum;
                old.calls += record.calls;
                old.flops += record.flops;
                old.nanoseconds += record.nanoseconds;
            } else {
                assert!(records.len() < CAPACITY, "GEMM shape capacity exceeded");
                records.push(record);
            }
        }
    }
}

/// Drain only at an outer reporting boundary, after all workers have joined.
pub fn take() -> Vec<Record> {
    if ENABLED {
        RECORDS
            .lock()
            .expect("GEMM accounting lock")
            .drain(..)
            .collect()
    } else {
        Vec::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::{ComplexMatrix, GemmWorkspace};

    /// P14B_GEMM_ACCOUNTING=1 cargo test --release -p dft-numerics accounting_shapes -- --ignored
    #[test]
    #[ignore = "requires source-only GEMM measurement build"]
    fn accounting_shapes() {
        assert!(std::hint::black_box(ENABLED));
        take();
        for threads in [1, 4] {
            let a = ComplexMatrix::zeros(513, 7).unwrap();
            let b = ComplexMatrix::zeros(7, 3).unwrap();
            let mut work = GemmWorkspace::new(513, 16, threads).unwrap();
            work.accounting_site("oracle_pair");
            work.rotate_pair(&a, &a, &b).unwrap();
            work.reduced_hamiltonian(&a, &a, None).unwrap();
            let rows = take();
            let pair = rows.iter().find(|r| r.site == "oracle_pair").unwrap();
            assert_eq!((pair.m, pair.n, pair.k, pair.calls), (513, 3, 7, 2));
            assert_eq!(pair.flops, 4 * 513 * 3 * 7);
            let gram = rows
                .iter()
                .find(|r| r.site == "reduced_hamiltonian")
                .unwrap();
            assert_eq!(gram.flops, 2 * 7 * 7 * 513);
            assert_eq!(gram.minimum, 2 * 7 * 7 * 513);
            assert_eq!(rows.len(), 2);
            let old = ComplexMatrix::zeros(7, 7).unwrap();
            let extended = ComplexMatrix::zeros(513, 11).unwrap();
            let result = work
                .reduced_hamiltonian(&extended, &extended, Some(&old))
                .unwrap();
            let rows = take();
            assert_eq!(rows.len(), 2);
            assert!(
                rows.iter()
                    .all(|r| r.site == "reduced_extension" && r.calls == 1)
            );
            assert_eq!(
                rows.iter().map(|r| r.flops).sum::<u64>(),
                2 * 513 * (11 * 11 - 7 * 7)
            );
            work.reduced_hamiltonian(&extended, &extended, Some(&result))
                .unwrap();
            assert!(take().is_empty(), "unchanged width requires no GEMM");
        }
    }
}
