//! LANE-VERIFY Phase 0. Per-attempt/rank records; no numerical decisions.
//! Enabled only by VERIFY_PROFILE plus PHASE_PROFILE. STAGE_TIMERS controls
//! all new clocks. Records are gathered by the adapter after the run, never
//! by a numerical loop. Counts fit exactly in the f64 MPI transport (< 2^53).
use std::{cell::RefCell, time::Instant};

pub const PARTS: [&str; 9] = [
    "owned_input",
    "warm_start",
    "map",
    "coulomb",
    "metrics",
    "decision",
    "record",
    "restore",
    "schedule",
];
pub const SIZES: [&str; 12] = [
    "grid_x",
    "grid_y",
    "grid_z",
    "full_density_elements",
    "owned_density_elements",
    "warm_complex_bytes",
    "k_blocks",
    "lazy_blocks",
    "bands_sum",
    "plane_waves_sum",
    "history_records",
    "previous_attempts",
];
pub const MAP_PARTS: [&str; 4] = ["potential", "solve", "density", "energy"];

#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    pub rank: usize,
    pub iteration: usize,
    pub cutoff_ha: f64,
    /// 0: verification attempt, 1: scheduler evaluation.
    pub event: u8,
    /// Attempt: -1 error, 0 rejected/clipped, 1 passed. Schedule: 0 proceed, 1 defer.
    pub status: i8,
    pub timed: bool,
    pub sizes: [u64; 12],
    pub seconds: [f64; 9],
    /// Nested INSIDE `seconds[2]`, never added to its total.
    pub map_seconds: [f64; 4],
    pub h_vectors: usize,
    pub davidson_iterations: usize,
}
thread_local! {
    static RECORDS: RefCell<Vec<Record>> = const { RefCell::new(Vec::new()) };
}

pub fn enabled() -> bool {
    dft_numerics::switches::on(dft_numerics::switches::Id::VERIFY_PROFILE) && super::enabled()
}

pub(crate) struct Trace {
    record: Record,
    mark: Option<Instant>,
}
impl Trace {
    pub(crate) fn new(iteration: usize, cutoff_ha: f64, sizes: [u64; 12], timed: bool) -> Self {
        Self {
            record: Record {
                rank: super::rank_key(),
                iteration,
                cutoff_ha,
                event: 0,
                status: -1,
                timed,
                sizes,
                seconds: [0.; 9],
                map_seconds: [0.; 4],
                h_vectors: 0,
                davidson_iterations: 0,
            },
            mark: timed.then(Instant::now),
        }
    }
    pub(crate) fn map(&mut self, seconds: [f64; 4]) {
        if self.record.timed {
            self.record.map_seconds = seconds;
        }
    }
    pub(crate) fn finish(mut self, result: Option<(bool, usize, usize)>) {
        if let Some((passed, h, iterations)) = result {
            self.record.status = i8::from(passed);
            self.record.h_vectors = h;
            self.record.davidson_iterations = iterations;
        }
        RECORDS.with(|r| r.borrow_mut().push(self.record));
    }
}
#[inline]
pub(crate) fn mark(trace: &mut Option<Trace>, part: usize) {
    if let Some(t) = trace.as_mut()
        && let Some(previous) = t.mark
    {
        let now = Instant::now();
        t.record.seconds[part] += now.duration_since(previous).as_secs_f64();
        t.mark = Some(now);
    }
}
/// Drain at a run boundary; each OS thread represents one real/simulated rank.
pub fn take() -> Vec<Record> {
    RECORDS.with(|r| std::mem::take(&mut *r.borrow_mut()))
}

/// Two adapter-only allgathers: lengths then records. Different numbers of
/// records per rank are supported; no assumption about rank-local map success.
pub fn gather<C: dft_mpi::Collective + ?Sized>(
    comm: &C,
    records: &[Record],
) -> Result<Vec<Record>, String> {
    const WIDTH: usize = 8 + SIZES.len() + PARTS.len() + MAP_PARTS.len();
    let mut lengths = vec![0.; comm.size()];
    comm.allgatherv_f64(
        &[(records.len() * WIDTH) as f64],
        &vec![1; comm.size()],
        &mut lengths,
    )
    .map_err(|e| e.to_string())?;
    let counts: Vec<_> = lengths.iter().map(|&n| n as usize).collect();
    let mut local = Vec::with_capacity(records.len() * WIDTH);
    for r in records {
        local.extend([
            r.rank as f64,
            r.iteration as f64,
            r.cutoff_ha,
            r.status as f64,
            u8::from(r.timed) as f64,
            r.h_vectors as f64,
            r.davidson_iterations as f64,
            r.event as f64,
        ]);
        local.extend(r.sizes.iter().map(|&n| n as f64));
        local.extend(r.seconds);
        local.extend(r.map_seconds);
    }
    let mut all = vec![0.; counts.iter().sum()];
    if all.is_empty() {
        return Ok(Vec::new());
    }
    comm.allgatherv_f64(&local, &counts, &mut all)
        .map_err(|e| e.to_string())?;
    Ok(all
        .as_chunks::<WIDTH>()
        .0
        .iter()
        .map(|v| Record {
            rank: v[0] as usize,
            iteration: v[1] as usize,
            cutoff_ha: v[2],
            status: v[3] as i8,
            timed: v[4] != 0.,
            h_vectors: v[5] as usize,
            davidson_iterations: v[6] as usize,
            event: v[7] as u8,
            sizes: std::array::from_fn(|i| v[8 + i] as u64),
            seconds: std::array::from_fn(|i| v[20 + i]),
            map_seconds: std::array::from_fn(|i| v[29 + i]),
        })
        .collect())
}

/// Scheduler-only event; called only for an actual predictor evaluation.
pub(crate) fn schedule(
    iteration: usize,
    cutoff_ha: f64,
    history: usize,
    attempts: usize,
    clock: Option<Instant>,
    deferred: bool,
) {
    let mut sizes = [0; 12];
    sizes[10] = history as u64;
    sizes[11] = attempts as u64;
    let mut seconds = [0.; 9];
    seconds[8] = clock.map_or(0., |c| c.elapsed().as_secs_f64());
    RECORDS.with(|r| {
        r.borrow_mut().push(Record {
            rank: super::rank_key(),
            iteration,
            cutoff_ha,
            event: 1,
            status: i8::from(deferred),
            timed: clock.is_some(),
            sizes,
            seconds,
            map_seconds: [0.; 4],
            h_vectors: 0,
            davidson_iterations: 0,
        })
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn gather_preserves_unequal_rank_records_and_count_only_payloads() {
        for empty in [false, true] {
            let out = dft_mpi::simulate(3, |comm: &dft_mpi::SimRank| {
                use dft_mpi::Collective;
                let records: Vec<_> = (0..if empty { 0 } else { comm.rank() })
                    .map(|i| Record {
                        rank: comm.rank(),
                        iteration: i + 1,
                        cutoff_ha: 2.5,
                        event: i as u8,
                        status: 0,
                        timed: false,
                        sizes: [comm.rank() as u64; 12],
                        seconds: [0.; 9],
                        map_seconds: [0.; 4],
                        h_vectors: 7,
                        davidson_iterations: 2,
                    })
                    .collect();
                let all = gather(comm, &records).unwrap();
                assert_eq!(
                    all.iter()
                        .filter(|r| r.rank == comm.rank())
                        .cloned()
                        .collect::<Vec<_>>(),
                    records
                );
                all
            });
            assert_eq!(out[0], out[1]);
            assert_eq!(out[1], out[2]);
            assert_eq!(out[0].len(), if empty { 0 } else { 3 });
        }
    }
}
