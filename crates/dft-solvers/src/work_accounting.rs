//! Source-only P14A measurement build (`P14A_WORK_ACCOUNTING=1 cargo build`).
//! No clocks, serialization, or output in the solver. Root arrays are allocated
//! once per solve; the bounded collector is reserved before entering iterations.
use std::sync::{LazyLock, Mutex};

pub const ENABLED: bool = option_env!("P14A_WORK_ACCOUNTING").is_some();
const CAPACITY: usize = 16384;
static RECORDS: LazyLock<Mutex<Vec<Record>>> =
    LazyLock::new(|| Mutex::new(Vec::with_capacity(CAPACITY)));

#[derive(Default)]
pub struct Record {
    pub bands: usize,
    pub strict_bands: usize,
    pub tolerance: f64,
    pub buffer_tolerance: f64,
    pub iterations: usize,
    pub block_calls: usize,
    pub vector_applications: usize,
    pub buffer_origin_applications: usize,
    pub unattributed_applications: usize,
    pub restarts: usize,
    pub first_below: Vec<usize>,
    pub last_above: Vec<usize>,
    pub returned: bool,
}

pub struct Trace(Option<Record>);
impl Trace {
    pub fn new(bands: usize, strict_bands: usize, tolerance: f64, buffer_tolerance: f64) -> Self {
        Self(ENABLED.then(|| {
            LazyLock::force(&RECORDS);
            Record {
                bands,
                strict_bands,
                tolerance,
                buffer_tolerance,
                first_below: vec![0; bands],
                last_above: vec![0; bands],
                ..Record::default()
            }
        }))
    }

    pub fn apply(&mut self, columns: usize, buffer: Option<usize>) {
        if let Some(r) = &mut self.0 {
            r.block_calls += 1;
            r.vector_applications += columns;
            if let Some(buffer) = buffer {
                r.buffer_origin_applications += buffer;
            } else {
                r.unattributed_applications += columns;
            }
        }
    }

    pub fn residuals(&mut self, iteration: usize, residuals: &[f64]) {
        if let Some(r) = &mut self.0 {
            r.iterations = iteration;
            for (j, &value) in residuals.iter().enumerate() {
                let tolerance = if j < r.strict_bands {
                    r.tolerance
                } else {
                    r.buffer_tolerance
                };
                if value <= tolerance {
                    if r.first_below[j] == 0 {
                        r.first_below[j] = iteration;
                    }
                } else {
                    r.last_above[j] = iteration;
                }
            }
        }
    }

    pub fn restart(&mut self) {
        if let Some(r) = &mut self.0 {
            r.restarts += 1;
        }
    }

    pub fn returned(&mut self) {
        if let Some(r) = &mut self.0 {
            r.returned = true;
        }
    }
}

impl Drop for Trace {
    fn drop(&mut self) {
        if let Some(record) = self.0.take() {
            let mut records = RECORDS.lock().expect("work accounting collector");
            assert!(
                records.len() < CAPACITY,
                "work accounting capacity exceeded"
            );
            records.push(record);
        }
    }
}

/// Called by the CLI after a complete SCF, never by a numerical kernel.
pub fn take() -> Vec<Record> {
    if ENABLED {
        RECORDS
            .lock()
            .expect("work accounting collector")
            .drain(..)
            .collect()
    } else {
        Vec::new()
    }
}
