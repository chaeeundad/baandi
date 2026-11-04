//! Incremental execution of the existing SCF iteration, with immutable boundaries.
use super::*;
use crate::events::{EventKind, EventQueue, ExecutionConfig, ExecutionMode, FailureCode};
use crate::policy::{PolicyChange, ToleranceFloors, prepare_policy};
use crate::runtime::ControlError;
use crate::state_view::StateSnapshot;

/// Prepared under the runtime's exclusive session ownership. All fallible work
/// precedes journal publication; commit performs only owned replacements.
pub(crate) struct BoundaryPolicy {
    options: ScfOptions,
    mixer: Option<Box<dyn MixerHistory>>,
}
pub(crate) struct BoundaryMode {
    mode: ExecutionMode,
    events: Option<EventQueue>,
    replace_events: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CompletionMode {
    Strict,
    RetainLast,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StepStatus {
    Created,
    Ready,
    Converged,
    NotConverged,
    Failed,
}
#[derive(Clone, Copy, Debug)]
pub struct StepReport {
    pub iteration: usize,
    pub state_version: u64,
    pub status: StepStatus,
    pub mixing_applied: bool,
}
struct Continuation {
    n: Array3<f64>,
    mixer: Box<dyn MixerHistory>,
    previous: Option<Vec<ComplexMatrix>>,
    history: Vec<IterationRecord>,
    weight: f64,
    best: f64,
    growing: usize,
    mixing_resets: usize,
    start: Instant,
}
/// Owns all continuation state. No background worker or hidden SCF loop.
/// Use `runtime::Runtime` for audited control; direct sessions retain the lean driver.
pub struct ScfSession {
    problem: ScfProblem,
    completion: CompletionMode,
    config: ExecutionConfig,
    continuation: Option<Continuation>,
    outcome: Option<ScfOutcome>,
    status: StepStatus,
    iteration: usize,
    state_version: u64,
    events: Option<EventQueue>,
    certificate_attempted: bool,
}
impl ScfSession {
    pub(crate) fn prepare_boundary_policy(
        &self,
        changes: &[PolicyChange],
        floors: Option<ToleranceFloors>,
    ) -> Result<BoundaryPolicy, ControlError> {
        if self.status != StepStatus::Ready || self.config.mode != ExecutionMode::AgentControlled {
            return Err(ControlError::ModeDisallowsControl);
        }
        self.check_control_version()?;
        let prepared = prepare_policy(&self.problem.options, self.iteration, floors, changes)
            .map_err(ControlError::Policy)?;
        let mixer = if prepared.resets_history() {
            let c = self
                .continuation
                .as_ref()
                .ok_or(ControlError::PrepareFailed)?;
            if c.mixing_resets == usize::MAX {
                return Err(ControlError::Capacity);
            }
            let mut mixer = self
                .make_mixer(&prepared.options().mixing)
                .map_err(|_| ControlError::PrepareFailed)?;
            mixer
                .mix(
                    c.n.as_slice().unwrap(),
                    c.n.as_slice().unwrap(),
                    self.problem.electrons / self.problem.cell.volume_bohr3(),
                )
                .map_err(|_| ControlError::PrepareFailed)?;
            mixer.reset();
            Some(mixer)
        } else {
            None
        };
        Ok(BoundaryPolicy {
            options: prepared.options().clone(),
            mixer,
        })
    }
    pub(crate) fn commit_boundary_policy(&mut self, prepared: BoundaryPolicy) {
        self.problem.options = prepared.options;
        if let Some(mixer) = prepared.mixer {
            let c = self.continuation.as_mut().expect("prepared continuation");
            c.mixer = mixer;
            c.mixing_resets += 1;
        }
        self.commit_control_version();
    }
    pub(crate) fn check_control_version(&self) -> Result<(), ControlError> {
        if self.state_version >= u64::MAX - 1 {
            Err(ControlError::Capacity)
        } else {
            Ok(())
        }
    }
    pub(crate) fn commit_control_version(&mut self) {
        self.state_version += 1;
    }
    pub(crate) fn prepare_boundary_mode(
        &self,
        mode: ExecutionMode,
    ) -> Result<BoundaryMode, ControlError> {
        self.check_control_version()?;
        if self.status != StepStatus::Created && self.status != StepStatus::Ready {
            return Err(ControlError::Terminal);
        }
        if mode == self.config.mode {
            return Err(ControlError::NoChange);
        }
        let created = self.status == StepStatus::Created;
        if !created
            && (mode == ExecutionMode::Classic || self.config.mode == ExecutionMode::Classic)
        {
            return Err(ControlError::ModeUnavailableAfterStart);
        }
        let events = if created && mode != ExecutionMode::Classic {
            Some(EventQueue::new(&self.config))
        } else {
            None
        };
        Ok(BoundaryMode {
            mode,
            events,
            replace_events: created,
        })
    }
    pub(crate) fn commit_boundary_mode(&mut self, prepared: BoundaryMode) {
        self.config.mode = prepared.mode;
        if prepared.replace_events {
            self.events = prepared.events;
        }
        self.commit_control_version();
    }
    pub fn new(
        problem: ScfProblem,
        completion: CompletionMode,
        config: ExecutionConfig,
    ) -> Result<Self, ScfError> {
        config.validate()?;
        let events = (config.mode != ExecutionMode::Classic).then(|| EventQueue::new(&config));
        Ok(Self {
            problem,
            completion,
            config,
            continuation: None,
            outcome: None,
            status: StepStatus::Created,
            iteration: 0,
            state_version: 0,
            events,
            certificate_attempted: false,
        })
    }
    pub fn initialize(&mut self) -> Result<(), ScfError> {
        if self.continuation.is_some() || self.status != StepStatus::Created {
            return Err(ScfError::InvalidInput(
                "session already initialized or failed",
            ));
        }
        match self.initialize_inner() {
            Ok(()) => {
                self.state_version += 1;
                self.status = StepStatus::Ready;
                if self.events.is_some() {
                    self.emit(EventKind::Initialized);
                }
                Ok(())
            }
            Err(e) => {
                self.fail(&e);
                Err(e)
            }
        }
    }
    fn make_mixer(&self, mixing: &Mixing) -> Result<Box<dyn MixerHistory>, ScfError> {
        let mixer: Box<dyn MixerHistory> = match *mixing {
            Mixing::Linear { beta } => Box::new(LinearMixer { beta }),
            Mixing::Pulay { beta, history } => Box::new(PulayMixer::new(beta, history)?),
            Mixing::KerkerPulay {
                beta,
                history,
                q0_inverse_bohr,
            } => {
                let dims = self.problem.grid.dims();
                let mut fft = FftPlan3d::new(dims)?;
                let mut source = Array3::<Complex64>::zeros(dims);
                let mut target = source.clone();
                let transform = move |values: &mut [Complex64], inverse: bool| {
                    source.as_slice_mut().unwrap().copy_from_slice(values);
                    let result = if inverse {
                        fft.to_real(&source, &mut target)
                    } else {
                        fft.to_reciprocal(&source, &mut target)
                    };
                    result.map_err(|e| dft_solvers::SolverError::Application(e.to_string()))?;
                    values.copy_from_slice(target.as_slice().unwrap());
                    Ok(())
                };
                let filter = Kerker::new(
                    dims,
                    self.problem.cell.reciprocal_columns_inverse_bohr(),
                    q0_inverse_bohr,
                    transform,
                )?;
                Box::new(PreconditionedPulay::new(
                    beta,
                    history,
                    dims.iter().product(),
                    filter,
                )?)
            }
        };
        Ok(mixer)
    }
    fn initialize_inner(&mut self) -> Result<(), ScfError> {
        let start = Instant::now();
        let mut mixer = self.make_mixer(&self.problem.options.mixing)?;
        // Validate mixer parameters even if a future system converges immediately.
        let n = self.problem.uniform_density();
        mixer.mix(
            n.as_slice().unwrap(),
            n.as_slice().unwrap(),
            self.problem.electrons / self.problem.cell.volume_bohr3(),
        )?;
        mixer.reset();
        let weight = self.problem.cell.volume_bohr3() / n.len() as f64;
        self.continuation = Some(Continuation {
            n,
            mixer,
            previous: None,
            history: Vec::new(),
            weight,
            best: f64::INFINITY,
            growing: 0,
            mixing_resets: 0,
            start,
        });
        Ok(())
    }
    fn advance(&mut self) -> Result<StepStatus, ScfError> {
        let retain_last = self.completion == CompletionMode::RetainLast;
        let iteration = self.iteration + 1;
        let Continuation {
            n,
            mixer,
            previous,
            history,
            weight,
            best,
            growing,
            mixing_resets,
            start,
        } = self.continuation.as_mut().unwrap();
        let weight = *weight;
        if n.iter()
            .any(|n| !n.is_finite() || *n < -xc::NEGATIVE_ROUNDOFF)
        {
            return Err(ScfError::Diverged {
                reason: "mixed density is nonfinite or negative beyond 1e-14 Bohr^-3".into(),
                history: std::mem::take(history),
            });
        }
        let tolerance = if self.problem.options.residual_eigen_schedule {
            // Residual is the normalized density L2 norm in atomic units.
            // 0.1 is an empirical scheduling coefficient in atomic units.
            history
                .last()
                .map_or(1e-2, |r| 0.1 * r.density_residual)
                .min(1e-2)
                .max(0.1 * self.problem.options.eigen_tolerance)
        } else {
            (self.problem.options.initial_eigen_tolerance
                * self
                    .problem
                    .options
                    .eigen_tolerance_decay
                    .powf((iteration - 1) as f64))
            .max(0.1 * self.problem.options.eigen_tolerance)
        };
        let step = self
            .problem
            .step_kpoints_with_certification(
                n,
                previous.as_deref(),
                tolerance,
                tolerance <= self.problem.options.eigen_tolerance,
            )
            .map_err(|e| match e {
                ScfError::NonFinite => ScfError::Diverged {
                    reason: "nonfinite SCF arithmetic".into(),
                    history: history.clone(),
                },
                other => other,
            })?;
        let residual = density_distance(
            &step.density,
            n,
            self.problem.cell.volume_bohr3(),
            self.problem.electrons,
        );
        let charge = step.density.sum() * weight;
        let charge_error = (charge - self.problem.electrons).abs() / self.problem.electrons;
        let input_charge_error =
            (n.sum() * weight - self.problem.electrons).abs() / self.problem.electrons;
        let energy_change = history
            .last()
            .map(|r| step.energy.total_hartree() - r.energy.total_hartree());
        let max_eigen_residual = step
            .k_points
            .iter()
            .flat_map(|k| k.solution.residuals.iter().zip(&self.problem.occupations))
            .filter(|(_, f)| **f > 0.0)
            .map(|(r, _)| *r)
            .fold(0.0, f64::max);
        let record = IterationRecord {
            fixed_point_certificate: None,
            mixing_fraction: 1.0,
            mixer_reset: false,
            density_contraction: history
                .last()
                .filter(|r| r.density_residual > 0.0)
                .map(|r| residual / r.density_residual),
            mixing_history_used: None,
            mixing_resets: *mixing_resets,
            mixing_raw_recovery: false,
            iteration,
            energy: step.energy.clone(),
            energy_change,
            density_residual: residual,
            max_eigen_residual,
            requested_eigen_tolerance: tolerance,
            davidson_iterations: step.k_points.iter().map(|k| k.solution.iterations).sum(),
            operator_applications: step
                .k_points
                .iter()
                .map(|k| k.solution.operator_applications)
                .sum(),
            davidson_iterations_per_k: step
                .k_points
                .iter()
                .map(|k| k.solution.iterations)
                .collect(),
            operator_applications_per_k: step
                .k_points
                .iter()
                .map(|k| k.solution.operator_applications)
                .collect(),
            charge_error,
            input_charge_error,
            min_density: step.density.iter().copied().fold(f64::INFINITY, f64::min),
            timings: step.timings,
        };
        let converged = converged(&record, &self.problem.options);
        history.push(record);
        if !residual.is_finite()
            || !step.energy.total_hartree().is_finite()
            || !charge_error.is_finite()
            || !max_eigen_residual.is_finite()
        {
            return Err(ScfError::Diverged {
                reason: "nonfinite SCF observable".into(),
                history: std::mem::take(history),
            });
        }
        if residual
            > self.problem.options.divergence_factor
                * best.max(self.problem.options.density_tolerance)
        {
            *growing += 1;
        } else {
            *growing = 0;
        }
        *best = best.min(residual);
        if *growing >= self.problem.options.divergence_patience {
            return Err(ScfError::Diverged {
                reason: "density residual exceeded growth limit".into(),
                history: std::mem::take(history),
            });
        }
        // Certify the returned output density itself, not just the input
        // fixed-point residual. A small-gap response can amplify that residual.
        let certified = if converged {
            self.certificate_attempted = true;
            let warm: Vec<_> = step
                .k_points
                .iter()
                .map(|k| k.solution.eigenvectors.clone())
                .collect();
            let check = self
                .problem
                .step_kpoints(
                    &step.density,
                    Some(&warm),
                    0.1 * self.problem.options.eigen_tolerance,
                )
                .map_err(|e| match e {
                    ScfError::NonFinite => ScfError::Diverged {
                        reason: "nonfinite certification arithmetic".into(),
                        history: history.clone(),
                    },
                    other => other,
                })?;
            let density_residual = density_distance(
                &check.density,
                &step.density,
                self.problem.cell.volume_bohr3(),
                self.problem.electrons,
            );
            let occupied_eigenvalue_change = check
                .k_points
                .iter()
                .zip(&step.k_points)
                .flat_map(|(a, b)| {
                    a.solution
                        .eigenvalues
                        .iter()
                        .zip(&b.solution.eigenvalues)
                        .zip(&self.problem.occupations)
                })
                .filter(|(_, f)| **f > 0.0)
                .map(|((a, b), _)| (a - b).abs())
                .fold(0.0, f64::max);
            let passed = density_residual <= self.problem.options.density_tolerance
                && occupied_eigenvalue_change <= 1e-7;
            history.last_mut().unwrap().fixed_point_certificate = Some(FixedPointCertificate {
                strict_residuals_per_k: check
                    .k_points
                    .iter()
                    .map(|k| {
                        k.solution
                            .residuals
                            .iter()
                            .zip(&k.solution.buffer)
                            .filter_map(|(&r, &buffer)| (!buffer).then_some(r))
                            .collect()
                    })
                    .collect(),
                buffer_residuals_per_k: check
                    .k_points
                    .iter()
                    .map(|k| {
                        k.solution
                            .residuals
                            .iter()
                            .zip(&k.solution.buffer)
                            .filter_map(|(&r, &buffer)| buffer.then_some(r))
                            .collect()
                    })
                    .collect(),
                density_residual,
                occupied_eigenvalue_change,
                davidson_iterations: check.k_points.iter().map(|k| k.solution.iterations).sum(),
                operator_applications: check
                    .k_points
                    .iter()
                    .map(|k| k.solution.operator_applications)
                    .sum(),
                davidson_iterations_per_k: check
                    .k_points
                    .iter()
                    .map(|k| k.solution.iterations)
                    .collect(),
                operator_applications_per_k: check
                    .k_points
                    .iter()
                    .map(|k| k.solution.operator_applications)
                    .collect(),
                timings: check.timings,
            });
            if passed {
                self.problem.check_gap(&check.k_points)?;
            }
            passed
        } else {
            false
        };
        if certified || (retain_last && iteration == self.problem.options.max_iterations) {
            let gap_hartree = if certified {
                self.problem.check_gap(&step.k_points)?
            } else {
                self.problem
                    .minimum_gap
                    .map(|_| self.problem.sampled_gap(&step.k_points))
            };
            self.outcome = Some(ScfOutcome {
                converged: certified,
                result: ScfResult {
                    k_points: step.k_points,
                    gap_hartree,
                    energy: step.energy,
                    eigenvalues: step.solution.eigenvalues,
                    occupations: self.problem.occupations.clone(),
                    orbitals: step.solution.eigenvectors,
                    density: step.density,
                    charge,
                    history: std::mem::take(history),
                    ewald_breakdown: ewald::EwaldEnergy {
                        ..self.problem.ewald
                    },
                    setup_time: self.problem.setup_time,
                    total_time: start.elapsed() + self.problem.setup_time,
                },
            });
            return Ok(if certified {
                StepStatus::Converged
            } else {
                StepStatus::NotConverged
            });
        }
        let mix_start = Instant::now();
        let mut mixed = mixer
            .mix(
                n.as_slice().unwrap(),
                step.density.as_slice().unwrap(),
                self.problem.electrons / self.problem.cell.volume_bohr3(),
            )
            .map_err(|e| ScfError::Diverged {
                reason: e.to_string(),
                history: history.clone(),
            })?;
        history.last_mut().unwrap().mixing_history_used =
            Some(mixer.history_len().saturating_sub(1));
        history.last_mut().unwrap().mixing_raw_recovery = mixer.used_raw_recovery();
        // Anderson extrapolation is unconstrained. A convex step toward its
        // proposal preserves charge and Fourier support while keeping XC in
        // its density domain. Never clip individual spatial samples.
        if matches!(
            self.problem.options.mixing,
            Mixing::Pulay { .. } | Mixing::KerkerPulay { .. }
        ) {
            let fraction = n
                .iter()
                .zip(&mixed)
                .filter(|(_, b)| **b < 0.0)
                .map(|(a, b)| 0.9 * a / (a - b))
                .fold(1.0, f64::min);
            if fraction < 1.0 {
                for (b, a) in mixed.iter_mut().zip(n.iter()) {
                    *b = a + fraction * (*b - a);
                }
                mixer.reset();
                let record = history.last_mut().unwrap();
                record.mixing_fraction = fraction;
                record.mixer_reset = true;
                *mixing_resets += 1;
                record.mixing_resets = *mixing_resets;
            }
        }
        n.as_slice_mut().unwrap().copy_from_slice(&mixed);
        history.last_mut().unwrap().timings.mixing = mix_start.elapsed();
        *previous = Some(
            step.k_points
                .into_iter()
                .map(|k| k.solution.eigenvectors)
                .collect(),
        );
        if iteration == self.problem.options.max_iterations {
            Ok(StepStatus::NotConverged)
        } else {
            Ok(StepStatus::Ready)
        }
    }
}

impl ScfSession {
    /// Exactly one complete iteration. Terminal/created calls reject without mutation.
    pub fn step(&mut self) -> Result<StepReport, ScfError> {
        if self.continuation.is_none() || self.status != StepStatus::Ready {
            return Err(ScfError::InvalidInput("session is not ready"));
        }
        if self.state_version == u64::MAX - 1 {
            let error = ScfError::Numerical("state version exhausted".into());
            self.fail(&error);
            return Err(error);
        }
        if self.events.as_ref().is_some_and(|q| !q.reserve_step()) {
            let error = ScfError::TelemetryExhausted;
            self.status = StepStatus::Failed;
            self.state_version += 1;
            self.emit(EventKind::Failed {
                code: FailureCode::TelemetryExhausted,
            });
            return Err(error);
        }
        self.certificate_attempted = false;
        let result = self.advance();
        match &result {
            Ok(status) => {
                self.status = *status;
                self.iteration += 1;
                self.state_version += 1;
            }
            Err(e) => {
                if let ScfError::Diverged { history, .. } = e {
                    self.continuation.as_mut().unwrap().history = history.clone();
                }
                self.fail(e);
                return result.map(|_| unreachable!());
            }
        }
        let mixing_applied = self.outcome.is_none();
        // Classic constructs no payloads and reads no telemetry clock.
        if self.events.is_some() {
            if let Some(c) = self
                .history()
                .last()
                .and_then(|r| r.fixed_point_certificate.as_ref())
            {
                self.emit(EventKind::CertificateEvaluated {
                    passed: self.status == StepStatus::Converged,
                    density_residual: c.density_residual.is_finite().then_some(c.density_residual),
                    occupied_eigenvalue_change_hartree: c
                        .occupied_eigenvalue_change
                        .is_finite()
                        .then_some(c.occupied_eigenvalue_change),
                });
            }
            let record = self.history().last().unwrap();
            let reset = record.mixer_reset;
            self.emit(EventKind::IterationFinished {
                energy_hartree: record.energy.total_hartree(),
                density_residual: record.density_residual,
                charge_error_relative: record.charge_error,
                mixing_applied,
            });
            if reset {
                self.emit(EventKind::MixerReset);
            }
            match self.status {
                StepStatus::Converged => self.emit(EventKind::Converged),
                StepStatus::NotConverged => self.emit(EventKind::NotConverged),
                _ => {}
            }
        }
        result.map(|_| StepReport {
            iteration: self.iteration,
            state_version: self.state_version,
            status: self.status,
            mixing_applied,
        })
    }
    fn emit(&mut self, kind: EventKind) {
        if let Some(queue) = &mut self.events {
            queue.push(kind, self.iteration, self.state_version);
        }
    }
    fn fail(&mut self, error: &ScfError) {
        self.status = StepStatus::Failed;
        self.state_version += 1;
        if self.events.is_some() {
            if self.certificate_attempted {
                let c = self
                    .history()
                    .last()
                    .and_then(|r| r.fixed_point_certificate.as_ref());
                self.emit(EventKind::CertificateEvaluated {
                    passed: false,
                    density_residual: c.map(|c| c.density_residual).filter(|v| v.is_finite()),
                    occupied_eigenvalue_change_hartree: c
                        .map(|c| c.occupied_eigenvalue_change)
                        .filter(|v| v.is_finite()),
                });
            }
            let code = match error {
                ScfError::NonFinite => FailureCode::NonFinite,
                ScfError::Diverged { .. } => FailureCode::Diverged,
                ScfError::NonInsulating { .. } => FailureCode::NonInsulating,
                ScfError::InvalidInput(_) => FailureCode::InvalidInput,
                _ => FailureCode::Numerical,
            };
            if code == FailureCode::Diverged {
                self.emit(EventKind::Diverged);
            }
            self.emit(EventKind::Failed { code });
        }
    }
    pub fn result(&self) -> Option<&ScfOutcome> {
        self.outcome.as_ref()
    }
    pub fn into_result(self) -> Option<ScfOutcome> {
        self.outcome
    }
    /// Explicit potentially large history access, excluded from events/snapshots.
    pub fn history(&self) -> &[IterationRecord] {
        self.outcome
            .as_ref()
            .map(|o| o.result.history.as_slice())
            .or_else(|| self.continuation.as_ref().map(|c| c.history.as_slice()))
            .unwrap_or(&[])
    }
    /// Next input density at an active boundary, returned output at a terminal one.
    pub fn density(&self) -> Option<&Array3<f64>> {
        self.outcome
            .as_ref()
            .map(|o| &o.result.density)
            .or_else(|| self.continuation.as_ref().map(|c| &c.n))
    }
    pub fn events(&mut self) -> Option<&mut EventQueue> {
        self.events.as_mut()
    }
    pub fn snapshot(&self) -> StateSnapshot<'_> {
        let record = self.history().last();
        StateSnapshot {
            schema_version: "1.0",
            calc_id: &self.config.calc_id,
            mode: self.config.mode,
            sequence: self.events.as_ref().map_or(0, EventQueue::sequence),
            state_version: self.state_version,
            iteration: self.iteration,
            status: self.status,
            record_iteration: record.map(|r| r.iteration),
            energy_change_hartree: record.and_then(|r| r.energy_change),
            requested_eigen_tolerance_hartree: record.map(|r| r.requested_eigen_tolerance),
            input_charge_error_relative: record.map(|r| r.input_charge_error),
            energy: record.map(|r| &r.energy),
            density_residual: record.map(|r| r.density_residual),
            max_eigen_residual_hartree: record.map(|r| r.max_eigen_residual),
            charge_error_relative: record.map(|r| r.charge_error),
            options: &self.problem.options,
        }
    }
}
