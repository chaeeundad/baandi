//! Synchronous, single-owner reference SCF continuation.
//!
//! Rust's pinned future stores the existing loop locals without copying numerical
//! arrays. It is polled only here, once per map: no executor, worker, wakeup,
//! callback, serialization or allocation is involved at a step boundary.
use crate::ScfError;
use std::{
    cell::{Cell, RefCell},
    future::Future,
    pin::Pin,
    rc::Rc,
    task::{Context, Poll, Waker},
};

thread_local! {
    // Existing solver/census scopes span yields. Reject interleaved suspended
    // runs on one thread; nested synchronous derivative solves remain legal.
    static SUSPENDED: Cell<Option<usize>> = const { Cell::new(None) };
}
type Continuation<'a, T> = Pin<Box<dyn Future<Output = Result<T, ScfError>> + 'a>>;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StepKind {
    Ordinary,
    Verification,
    Retry,
}

/// All energies/estimators have explicit units. H columns count only this map,
/// including its inner Davidson retries, never previous SCF maps.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StepReport {
    pub iteration: usize,
    pub kind: StepKind,
    pub stage: &'static str,
    pub estimator_ha: f64,
    pub dr2_ry: f64,
    pub energy_ha: Option<f64>,
    pub ethr_ry: f64,
    pub h_columns: usize,
    /// A reused failed verification is an ordinary map with zero new H work.
    pub reused: bool,
    /// False while an ordinary candidate awaits verification or a first-map
    /// retry. No mixing/control intervention is legal at such a boundary.
    pub safe_point: bool,
}

impl StepReport {
    pub(super) fn new(
        iteration: usize,
        kind: StepKind,
        dr2_ry: f64,
        energy: f64,
        ethr_ry: f64,
        h_columns: usize,
    ) -> Self {
        Self {
            iteration,
            kind,
            stage: "fine",
            estimator_ha: dr2_ry / 2.,
            dr2_ry,
            energy_ha: energy.is_finite().then_some(energy),
            ethr_ry,
            h_columns,
            reused: false,
            safe_point: false,
        }
    }
}

type StateResponse = Result<(Vec<u8>, usize), ScfError>;

#[derive(Default)]
pub(super) struct Boundary {
    report: Cell<Option<StepReport>>,
    initialized: Cell<bool>,
    request: RefCell<Option<Option<Vec<u8>>>>,
    response: RefCell<Option<StateResponse>>,
    state_available: Cell<bool>,
}
impl Boundary {
    pub fn publish(&self, report: StepReport) {
        assert!(self.report.replace(Some(report)).is_none());
    }
    pub async fn emit(&self, report: StepReport) {
        self.publish(report);
        self.pause().await;
    }
    pub async fn ready_state(
        &self,
        mut visit: impl FnMut(&mut super::checkpoint::Exchange<'_>) -> Result<usize, ScfError>,
    ) {
        if !self.initialized.replace(true) {
            self.pause_state(&mut visit).await;
        }
    }
    pub async fn emit_state(
        &self,
        report: StepReport,
        mut visit: impl FnMut(&mut super::checkpoint::Exchange<'_>) -> Result<usize, ScfError>,
    ) {
        self.publish(report);
        self.pause_state(&mut visit).await;
    }
    async fn pause_state(
        &self,
        visit: &mut impl FnMut(&mut super::checkpoint::Exchange<'_>) -> Result<usize, ScfError>,
    ) {
        let mut yielded = false;
        std::future::poll_fn(|_| {
            self.state_available.set(true);
            if let Some(request) = self.request.borrow_mut().take() {
                let mut exchange = match &request {
                    Some(bytes) => super::checkpoint::Exchange::load(bytes),
                    None => super::checkpoint::Exchange::save(),
                };
                let result = visit(&mut exchange)
                    .and_then(|iteration| exchange.finish().map(|bytes| (bytes, iteration)));
                *self.response.borrow_mut() = Some(result);
                yielded = true;
                return Poll::Pending;
            }
            if std::mem::replace(&mut yielded, true) {
                self.state_available.set(false);
                Poll::Ready(())
            } else {
                Poll::Pending
            }
        })
        .await;
    }
    async fn pause(&self) {
        let mut yielded = false;
        std::future::poll_fn(|_| {
            if std::mem::replace(&mut yielded, true) {
                Poll::Ready(())
            } else {
                Poll::Pending
            }
        })
        .await;
    }
}

/// Construct through `driver::stepper[_with]`, `k_driver::stepper_with`,
/// `electronic::stepper_with` or `api::stepper`. Initialize exactly once before
/// stepping. Finish consumes only an already delivered result, never runs a map.
/// Keep the run configuration scope installed and drive all ranks in lockstep.
/// This object is deliberately !Send/!Sync: solver scopes are thread-local.
pub struct Stepper<'a, T> {
    continuation: Option<Continuation<'a, T>>,
    boundary: Rc<Boundary>,
    result: Option<T>,
    initialized: bool,
    done: bool,
    safe: bool,
    iteration: usize,
}
impl<'a, T> Stepper<'a, T> {
    pub(super) fn new<F: Future<Output = Result<T, ScfError>> + 'a>(
        build: impl FnOnce(Rc<Boundary>) -> F,
    ) -> Self {
        let boundary = Rc::new(Boundary::default());
        Self {
            continuation: Some(Box::pin(build(boundary.clone()))),
            boundary,
            result: None,
            initialized: false,
            done: false,
            safe: true,
            iteration: 0,
        }
    }
    fn advance(&mut self) -> Result<(), ScfError> {
        SUSPENDED.set(None);
        let mut cx = Context::from_waker(Waker::noop());
        let poll = self
            .continuation
            .as_mut()
            .expect("live continuation")
            .as_mut()
            .poll(&mut cx);
        if let Poll::Ready(result) = poll {
            self.done = true;
            // Release all thread-local scopes even when numerical execution fails.
            self.continuation = None;
            self.result = Some(result?);
        } else {
            SUSPENDED.set(Some(self.owner()));
        }
        Ok(())
    }
    fn owner(&self) -> usize {
        Rc::as_ptr(&self.boundary) as usize
    }
    pub fn initialize(&mut self) -> Result<(), ScfError> {
        if self.initialized || self.done {
            return Err(ScfError::InvalidInput(
                "reference stepper: already initialized or terminal",
            ));
        }
        if SUSPENDED.get().is_some() {
            return Err(ScfError::InvalidInput(
                "reference stepper: another run is suspended on this thread",
            ));
        }
        self.initialized = true;
        self.advance()
    }
    pub fn step(&mut self) -> Result<StepReport, ScfError> {
        if !self.initialized || self.done {
            return Err(ScfError::InvalidInput(
                "reference stepper: initialize before stepping; terminal runs cannot step",
            ));
        }
        self.advance()?;
        let report = self.boundary.report.take().ok_or(ScfError::InvalidInput(
            "reference stepper: missing map boundary",
        ))?;
        self.safe = report.safe_point;
        self.iteration = report.iteration;
        Ok(report)
    }
    /// Capture without running another map. All ranks call at the same safe point.
    pub fn checkpoint(
        &mut self,
        identity: super::checkpoint::Identity,
        rank: usize,
    ) -> Result<super::checkpoint::Checkpoint, ScfError> {
        identity.validate()?;
        if !self.initialized || self.done || !self.safe || !self.boundary.state_available.get() {
            return Err(super::checkpoint::invalid(
                "capture requires a live safe point",
            ));
        }
        if rank >= identity.ranks {
            return Err(super::checkpoint::invalid("rank"));
        }
        *self.boundary.request.borrow_mut() = Some(None);
        self.advance()?;
        let (state, iteration) = self
            .boundary
            .response
            .borrow_mut()
            .take()
            .ok_or_else(|| super::checkpoint::invalid("missing state response"))??;
        if iteration != self.iteration {
            return Err(super::checkpoint::invalid("iteration boundary mismatch"));
        }
        Ok(super::checkpoint::Checkpoint {
            schema_version: super::checkpoint::SCHEMA_VERSION,
            identity,
            rank,
            iteration,
            state,
        })
    }
    /// Restore into a fresh stepper with the same input/configuration. Initialization
    /// rebuilds immutable plans; no completed SCF maps are replayed.
    pub fn restore(
        &mut self,
        checkpoint: &super::checkpoint::Checkpoint,
        expected: &super::checkpoint::Identity,
        rank: usize,
    ) -> Result<(), ScfError> {
        if self.initialized || self.done {
            return Err(super::checkpoint::invalid(
                "restore requires an uninitialized stepper",
            ));
        }
        checkpoint.validate(expected, rank)?;
        *self.boundary.request.borrow_mut() = Some(Some(checkpoint.state.clone()));
        self.initialize()?;
        let response = self
            .boundary
            .response
            .borrow_mut()
            .take()
            .ok_or_else(|| super::checkpoint::invalid("driver has no checkpoint boundary"))
            .and_then(|r| r);
        match response {
            Ok((_, iteration)) if iteration == checkpoint.iteration => {
                self.iteration = iteration;
                self.safe = true;
                Ok(())
            }
            result => {
                self.continuation = None;
                self.done = true;
                SUSPENDED.set(None);
                Err(result
                    .err()
                    .unwrap_or_else(|| super::checkpoint::invalid("iteration metadata mismatch")))
            }
        }
    }
    pub fn is_done(&self) -> bool {
        self.done
    }
    pub fn result(&self) -> Option<&T> {
        self.result.as_ref()
    }
    pub fn finish(mut self) -> Result<T, ScfError> {
        self.result.take().ok_or(ScfError::InvalidInput(
            "reference stepper: no delivered result",
        ))
    }
    pub fn run_to_end(mut self) -> Result<T, ScfError> {
        if !self.initialized {
            self.initialize()?;
        }
        while !self.is_done() {
            self.step()?;
        }
        self.finish()
    }
}

impl<T> Drop for Stepper<'_, T> {
    fn drop(&mut self) {
        // Drop numerical scopes before allowing another run on this thread.
        self.continuation = None;
        if SUSPENDED.get() == Some(self.owner()) {
            SUSPENDED.set(None);
        }
    }
}
