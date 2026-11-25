//! Bounded boundary telemetry. No serialization or external callbacks in SCF.
use crate::ScfError;
use std::{collections::VecDeque, sync::Arc, time::SystemTime};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ExecutionMode {
    #[default]
    Classic,
    Observed,
    /// Observation with no interventions; control transactions belong to P11B.
    AgentControlled,
}
#[derive(Clone, Debug)]
pub struct ExecutionConfig {
    pub mode: ExecutionMode,
    pub calc_id: String,
    pub routine_capacity: usize,
    pub critical_capacity: usize,
}
impl Default for ExecutionConfig {
    fn default() -> Self {
        Self {
            mode: ExecutionMode::Classic,
            calc_id: "scf".into(),
            routine_capacity: 256,
            critical_capacity: 128,
        }
    }
}
impl ExecutionConfig {
    pub(crate) fn validate(&self) -> Result<(), ScfError> {
        if self.calc_id.is_empty()
            || self.calc_id.len() > 128
            || !self.calc_id.is_ascii()
            || self.critical_capacity < 6
            || self.critical_capacity > 65536
            || self.routine_capacity > 65536
        {
            return Err(ScfError::InvalidInput(
                "invalid execution identity or queue capacity",
            ));
        }
        Ok(())
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FailureCode {
    NonFinite,
    Diverged,
    Numerical,
    InvalidInput,
    NonInsulating,
    TelemetryExhausted,
}
#[derive(Clone, Copy, Debug)]
pub enum EventKind {
    Initialized,
    IterationFinished {
        energy_hartree: f64,
        density_residual: f64,
        charge_error_relative: f64,
        mixing_applied: bool,
    },
    CertificateEvaluated {
        passed: bool,
        density_residual: Option<f64>,
        occupied_eigenvalue_change_hartree: Option<f64>,
    },
    MixerReset,
    Diverged,
    Converged,
    NotConverged,
    Failed {
        code: FailureCode,
    },
    EventsDropped {
        count: u64,
        first_missing_seq: u64,
        last_missing_seq: u64,
    },
}
#[derive(Clone, Debug)]
pub struct ScfEvent {
    pub schema_version: &'static str,
    pub calc_id: Arc<str>,
    pub sequence: u64,
    /// Wall clock at the completed boundary, never used for numerical decisions.
    pub timestamp: SystemTime,
    pub iteration: usize,
    pub state_version: u64,
    pub kind: EventKind,
}
/// Called explicitly by a consumer after stepping, never from a numerical step.
/// External I/O should consume popped events in an adapter worker.
pub trait Observer {
    fn observe(&mut self, event: &ScfEvent);
}
#[derive(Default)]
pub struct NoopObserver;
impl Observer for NoopObserver {
    fn observe(&mut self, _: &ScfEvent) {}
}

/// Single-owner nonblocking producer queue; adapters own synchronization/fan-out.
/// Routine events are dropped when their partition is full. Critical events and
/// loss summaries use separate credits, plus two emergency slots. No journal is
/// implemented here: consumers must drain/persist between steps. Exhaustion is a
/// terminal resource error before another eigensolve, never an implicit wait.
/// Pending loss is flushed before critical events and every 64 iterations.
pub struct EventQueue {
    calc_id: Arc<str>,
    routine: VecDeque<ScfEvent>,
    critical: VecDeque<ScfEvent>,
    routine_capacity: usize,
    critical_capacity: usize,
    sequence: u64,
    dropped: u64,
    first_missing: u64,
    last_missing: u64,
}
impl EventQueue {
    pub(crate) fn new(config: &ExecutionConfig) -> Self {
        Self {
            calc_id: Arc::from(config.calc_id.as_str()),
            routine: VecDeque::with_capacity(config.routine_capacity),
            critical: VecDeque::with_capacity(config.critical_capacity + 2),
            routine_capacity: config.routine_capacity,
            critical_capacity: config.critical_capacity,
            sequence: 0,
            dropped: 0,
            first_missing: 0,
            last_missing: 0,
        }
    }
    pub(crate) fn reserve_step(&self) -> bool {
        self.critical_capacity.saturating_sub(self.critical.len()) >= 6
            && self.sequence < u64::MAX - 16
    }
    pub fn sequence(&self) -> u64 {
        self.sequence
    }
    pub fn len(&self) -> usize {
        self.routine.len() + self.critical.len()
    }
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
    pub fn pop(&mut self) -> Option<ScfEvent> {
        match (self.routine.front(), self.critical.front()) {
            (Some(r), Some(c)) if r.sequence < c.sequence => self.routine.pop_front(),
            (_, Some(_)) => self.critical.pop_front(),
            _ => self.routine.pop_front(),
        }
    }
    pub fn drain_to(&mut self, observer: &mut impl Observer) {
        while let Some(event) = self.pop() {
            observer.observe(&event);
        }
    }
    pub(crate) fn push(&mut self, kind: EventKind, iteration: usize, state_version: u64) {
        let routine = matches!(kind, EventKind::IterationFinished { .. });
        if self.dropped > 0 && (!routine || iteration.is_multiple_of(64)) {
            let summary = EventKind::EventsDropped {
                count: self.dropped,
                first_missing_seq: self.first_missing,
                last_missing_seq: self.last_missing,
            };
            self.dropped = 0;
            self.push(summary, iteration, state_version);
        }
        self.sequence += 1;
        if routine && self.routine.len() == self.routine_capacity {
            if self.dropped == 0 {
                self.first_missing = self.sequence;
            }
            self.dropped += 1;
            self.last_missing = self.sequence;
            return;
        }
        let event = ScfEvent {
            schema_version: "1.0",
            calc_id: self.calc_id.clone(),
            sequence: self.sequence,
            timestamp: SystemTime::now(),
            iteration,
            state_version,
            kind,
        };
        if routine {
            self.routine.push_back(event);
        } else {
            assert!(self.critical.len() < self.critical_capacity + 2);
            self.critical.push_back(event);
        }
    }
}

impl EventQueue {
    /// The adapter must drain and durably publish telemetry before checkpointing.
    pub(crate) fn checkpoint_sequence(&self) -> Result<u64, crate::checkpoint::Rejection> {
        if !self.routine.is_empty() || !self.critical.is_empty() || self.dropped != 0 {
            return Err(crate::checkpoint::Rejection::UndrainedEvents);
        }
        Ok(self.sequence)
    }
    pub(crate) fn restore_sequence(&mut self, sequence: u64) {
        self.sequence = sequence;
    }
}
