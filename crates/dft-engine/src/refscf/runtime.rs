//! Single-owner lifecycle adapter for a reference continuation.
//!
//! Controls are trusted host operations, not deserializable agent authority.
//! Transport approval, audit persistence and collective command distribution
//! belong to the caller. Every rank must receive the same control sequence.
use super::stepper::{StepReport, Stepper};
use crate::{ScfError, runtime::RuntimeState};

/// The numerical driver alone determines convergence.
pub trait Completion {
    fn converged(&self) -> bool;
}
impl Completion for super::driver::Run {
    fn converged(&self) -> bool {
        self.outcome.converged
    }
}
impl Completion for super::electronic::Run {
    fn converged(&self) -> bool {
        self.converged
    }
}
impl Completion for super::api::Delivered {
    fn converged(&self) -> bool {
        match &self.run {
            super::api::Run::Scalar(run) => run.converged(),
            super::api::Run::Electronic(run) => run.converged(),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Control {
    Pause,
    Resume,
    Cancel,
}

/// Scalar partial record, including the last completed map on cancellation or
/// failure. It never claims that a cancelled continuation delivered a density.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Snapshot {
    pub state: RuntimeState,
    pub state_version: u64,
    pub completed_maps: usize,
    pub last_step: Option<StepReport>,
    pub safe_point: bool,
}

/// Pausing retains the live continuation; cancellation drops it only at a safe
/// boundary. No map is replayed and no numerical arrays are copied.
pub struct Runtime<'a, T: Completion> {
    stepper: Option<Stepper<'a, T>>,
    snapshot: Snapshot,
    pending: Option<Control>,
}
impl<'a, T: Completion> Runtime<'a, T> {
    pub fn new(stepper: Stepper<'a, T>) -> Self {
        Self {
            stepper: Some(stepper),
            snapshot: Snapshot {
                state: RuntimeState::Created,
                state_version: 0,
                completed_maps: 0,
                last_step: None,
                safe_point: false,
            },
            pending: None,
        }
    }
    pub fn snapshot(&self) -> Snapshot {
        self.snapshot
    }
    fn transition(&mut self, state: RuntimeState) {
        self.snapshot.state = state;
        self.snapshot.state_version += 1;
    }
    fn failure(&mut self) {
        self.stepper = None;
        self.pending = None;
        self.snapshot.safe_point = false;
        self.transition(RuntimeState::Failed);
    }
    pub fn initialize(&mut self) -> Result<(), ScfError> {
        if self.snapshot.state != RuntimeState::Created {
            return Err(invalid("reference runtime: initialize exactly once"));
        }
        if let Err(e) = self.stepper.as_mut().expect("created").initialize() {
            self.failure();
            return Err(e);
        }
        self.snapshot.safe_point = true;
        self.transition(RuntimeState::Initialized);
        Ok(())
    }
    /// Queue an approved host control. At an unsafe boundary, pause/cancel wait
    /// for the following safe return. Cancellation takes precedence over pause.
    pub fn control(&mut self, control: Control) -> Result<(), ScfError> {
        let state = self.snapshot.state;
        if state == RuntimeState::Created || state.terminal() {
            return Err(invalid(
                "reference runtime: control requires a live initialized run",
            ));
        }
        match control {
            Control::Resume => {
                if state != RuntimeState::Paused {
                    return Err(invalid("reference runtime: resume requires paused state"));
                }
                self.transition(RuntimeState::Initialized);
            }
            Control::Pause if state == RuntimeState::Paused => {}
            Control::Pause | Control::Cancel => {
                if self.pending != Some(Control::Cancel) {
                    self.pending = Some(control);
                }
                self.apply_pending();
            }
        }
        Ok(())
    }
    fn apply_pending(&mut self) {
        if self.snapshot.state.terminal() {
            // Terminal numerical delivery wins over pending host controls.
            self.pending = None;
        } else if self.snapshot.safe_point {
            match self.pending.take() {
                Some(Control::Cancel) => {
                    self.stepper = None;
                    self.transition(RuntimeState::Cancelled);
                }
                Some(Control::Pause) => self.transition(RuntimeState::Paused),
                _ => {}
            }
        }
    }
    pub fn step(&mut self) -> Result<StepReport, ScfError> {
        if self.snapshot.state != RuntimeState::Initialized {
            return Err(invalid(
                "reference runtime: step requires an unpaused live run",
            ));
        }
        let stepper = self.stepper.as_mut().expect("initialized");
        let report = match stepper.step() {
            Ok(report) => report,
            Err(e) => {
                self.failure();
                return Err(e);
            }
        };
        let state = if stepper.is_done() {
            if stepper.result().expect("successful delivery").converged() {
                RuntimeState::Converged
            } else {
                RuntimeState::Exhausted
            }
        } else {
            RuntimeState::Initialized
        };
        self.snapshot.last_step = Some(report);
        self.snapshot.completed_maps += 1;
        self.snapshot.safe_point = report.safe_point;
        self.transition(state);
        self.apply_pending();
        Ok(report)
    }
    pub fn result(&self) -> Option<&T> {
        self.stepper.as_ref().and_then(Stepper::result)
    }
    pub fn finish(mut self) -> Result<T, ScfError> {
        self.stepper
            .take()
            .ok_or_else(|| invalid("reference runtime: no delivered result"))?
            .finish()
    }
}
fn invalid(message: &'static str) -> ScfError {
    ScfError::InvalidInput(message)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::refscf::stepper::StepKind;

    struct Delivered(bool);
    impl Completion for Delivered {
        fn converged(&self) -> bool {
            self.0
        }
    }
    fn report(safe_point: bool) -> StepReport {
        StepReport {
            iteration: 1,
            kind: StepKind::Ordinary,
            stage: "fine",
            estimator_ha: 1.,
            dr2_ry: 2.,
            energy_ha: None,
            ethr_ry: 0.1,
            h_columns: 4,
            reused: false,
            safe_point,
        }
    }
    #[test]
    fn unsafe_pause_defers_and_cancel_supersedes_it() {
        let stepper = Stepper::new(|boundary| async move {
            boundary.ready().await;
            boundary.emit(report(false)).await;
            boundary.emit(report(true)).await;
            boundary.publish(report(true));
            Ok(Delivered(true))
        });
        let mut runtime = Runtime::new(stepper);
        runtime.initialize().unwrap();
        runtime.step().unwrap();
        runtime.control(Control::Pause).unwrap();
        assert_eq!(runtime.snapshot().state, RuntimeState::Initialized);
        runtime.control(Control::Cancel).unwrap();
        assert_eq!(runtime.snapshot().completed_maps, 1);
        runtime.step().unwrap();
        assert_eq!(runtime.snapshot().state, RuntimeState::Cancelled);
        assert_eq!(runtime.snapshot().completed_maps, 2);
        assert!(runtime.result().is_none());
    }
    #[test]
    fn terminal_numerics_win_over_pending_cancel() {
        for converged in [false, true] {
            let stepper = Stepper::new(|boundary| async move {
                boundary.ready().await;
                boundary.emit(report(false)).await;
                boundary.publish(report(true));
                Ok(Delivered(converged))
            });
            let mut runtime = Runtime::new(stepper);
            runtime.initialize().unwrap();
            runtime.step().unwrap();
            runtime.control(Control::Cancel).unwrap();
            runtime.step().unwrap();
            assert_eq!(
                runtime.snapshot().state,
                if converged {
                    RuntimeState::Converged
                } else {
                    RuntimeState::Exhausted
                }
            );
            assert_eq!(runtime.finish().unwrap().converged(), converged);
        }
    }
    #[test]
    fn numerical_failure_keeps_last_completed_map() {
        let stepper: Stepper<'_, Delivered> = Stepper::new(|boundary| async move {
            boundary.ready().await;
            boundary.emit(report(true)).await;
            Err(invalid("test numerical failure"))
        });
        let mut runtime = Runtime::new(stepper);
        runtime.initialize().unwrap();
        runtime.step().unwrap();
        assert!(runtime.step().is_err());
        let snapshot = runtime.snapshot();
        assert_eq!(snapshot.state, RuntimeState::Failed);
        assert_eq!(snapshot.last_step, Some(report(true)));
        assert_eq!(snapshot.completed_maps, 1);
        assert!(!snapshot.safe_point);
        assert!(runtime.finish().is_err());
    }
}
