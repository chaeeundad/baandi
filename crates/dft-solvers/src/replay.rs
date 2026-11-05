//! P15A source-only frozen-operator measurement. Enable at compile time with
//! `P15A_REPLAY=1`. All output is drained at the existing CLI reporting boundary.
//! Replays borrow the same immutable operator (including potential/nonlocal data)
//! and initial columns, never feed their result back to SCF, and use fixed seeds.
use std::sync::{LazyLock, Mutex};

pub const ENABLED: bool = option_env!("P15A_REPLAY").is_some();
static RECORDS: LazyLock<Mutex<Vec<Record>>> = LazyLock::new(|| Mutex::new(Vec::new()));
static METRICS: Mutex<Vec<(usize, f64, f64)>> = Mutex::new(Vec::new());

#[derive(Default)]
pub(crate) struct Trace {
    pub bands: usize,
    pub rows: Vec<(usize, usize, usize, usize)>,
    pub values: Vec<f64>,
    pub residuals: Vec<f64>,
    pub blocks: usize,
    pub restarts: usize,
}
impl Trace {
    pub fn new(bands: usize, iterations: usize) -> Self {
        Self {
            bands,
            rows: Vec::with_capacity(iterations),
            values: Vec::with_capacity(iterations * bands),
            residuals: Vec::with_capacity(iterations * bands),
            ..Self::default()
        }
    }
    pub fn step(&mut self, basis: usize, vectors: usize, values: &[f64], residuals: &[f64]) {
        assert!(
            self.rows.len() < self.rows.capacity(),
            "replay trace capacity"
        );
        self.rows.push((basis, vectors, self.blocks, self.restarts));
        self.values.extend_from_slice(values);
        self.residuals.extend_from_slice(residuals);
    }
}
pub(crate) struct Record {
    pub traces: Vec<(f64, usize, f64, bool, Trace)>,
    pub dense_values: Option<Vec<f64>>,
}
pub(crate) fn store(record: Record) {
    RECORDS.lock().expect("replay records").push(record);
}
/// Capture-only engine hook. dr2 is in Ry; electron count is dimensionless.
pub fn metric(iteration: usize, dr2_ry: f64, electrons: f64) {
    METRICS
        .lock()
        .expect("replay metrics")
        .push((iteration, dr2_ry, electrons));
}
pub(crate) fn drain() {
    if !ENABLED {
        return;
    }
    for (id, record) in RECORDS
        .lock()
        .expect("replay records")
        .drain(..)
        .enumerate()
    {
        if let Some(values) = record.dense_values {
            eprintln!("P15A_DENSE,{id},{values:?}");
        }
        for (variant, (tolerance, strict, buffer, success, trace)) in
            record.traces.into_iter().enumerate()
        {
            eprintln!("P15A_SOLVE,{id},{variant},{tolerance:?},{strict},{buffer:?},{success}");
            for (step, &(basis, vectors, blocks, restarts)) in trace.rows.iter().enumerate() {
                let range = step * trace.bands..(step + 1) * trace.bands;
                eprintln!(
                    "P15A_STEP,{id},{variant},{},{basis},{vectors},{blocks},{restarts},{:?},{:?}",
                    step + 1,
                    &trace.values[range.clone()],
                    &trace.residuals[range]
                );
            }
        }
    }
    for (iteration, dr2, electrons) in METRICS.lock().expect("replay metrics").drain(..) {
        eprintln!("P15A_METRIC,{iteration},{dr2:?},{electrons:?}");
    }
}

/// Optional independent full-space spectrum for selected sequential solve IDs.
/// This is bounded measurement code, not a change to DenseEigensolver's public
/// dimension guard. All true residuals are checked before using the oracle.
pub(crate) fn dense_oracle(
    operator: &dyn crate::LinearOperator,
    bands: usize,
) -> Result<Option<Vec<f64>>, crate::SolverError> {
    let id = RECORDS.lock().expect("replay records").len();
    let selected = option_env!("P15A_DENSE_IDS")
        .is_some_and(|ids| ids.split(',').any(|s| s.parse::<usize>().ok() == Some(id)));
    let n = operator.dimension();
    if n > 128 && !selected {
        return Ok(None);
    }
    if n > 4096 {
        return Err(crate::SolverError::DenseLimit {
            dimension: n,
            limit: 4096,
        });
    }
    let mut identity = crate::ComplexMatrix::zeros(n, n)?;
    for i in 0..n {
        identity[(i, i)] = crate::Complex64::new(1.0, 0.0);
    }
    let matrix = crate::apply(operator, &identity)?;
    drop(identity);
    let (mut values, mut vectors) = crate::hermitian_eigh(&matrix)?;
    values.truncate(bands);
    vectors.resize_columns(bands)?;
    let (_, residuals) = crate::residual(operator, &vectors, &values)?;
    if residuals.iter().any(|&r| r > 1e-10) {
        return Err(crate::SolverError::Application(
            "dense replay oracle residual exceeds 1e-10 Ha".into(),
        ));
    }
    Ok(Some(values))
}
