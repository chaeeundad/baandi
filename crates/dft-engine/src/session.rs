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
    spectral_report: Option<crate::SpectralReport>,
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
        if prepared.options().preconditioner != dft_solvers::DavidsonPreconditioner::Default
            && !matches!(
                prepared.options().eigensolver,
                dft_solvers::SolverChoice::Davidson
            )
        {
            return Err(ControlError::PrepareFailed);
        }
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
        let spectral_report = (problem.options.spectral_guard == crate::SpectralGuard::Final)
            .then(crate::SpectralReport::default);
        Ok(Self {
            spectral_report,
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
    /// Start a new, uncertified SCF with a validated fixed-cell density transfer.
    /// This resets orbitals, mixer and history and is never exact continuation.
    pub fn initialize_density_transfer(
        &mut self,
        source: &ScfProblem,
        density: &Array3<f64>,
    ) -> Result<(), ScfError> {
        self.initialize_transferred_density(source, density, false)
    }
    /// Affine cell warm transfer on the identical fractional FFT grid and frozen
    /// Miller support. Samples scale by V_old/V_new, preserving quadrature charge.
    /// No interpolation, orbitals, mixer history or certificate are transferred.
    pub fn initialize_cell_density_transfer(
        &mut self,
        source: &ScfProblem,
        density: &Array3<f64>,
    ) -> Result<(), ScfError> {
        let ratio = source.cell.volume_bohr3() / self.problem.cell.volume_bohr3();
        let transferred = density.mapv(|v| v * ratio);
        self.initialize_transferred_density(source, &transferred, true)
    }
    fn initialize_transferred_density(
        &mut self,
        source: &ScfProblem,
        density: &Array3<f64>,
        changed_cell: bool,
    ) -> Result<(), ScfError> {
        let matches = if changed_cell {
            self.problem
                .force_model
                .cell_transfer_matches(&source.force_model)
        } else {
            self.problem
                .force_model
                .transfer_matches(&source.force_model)
        };
        if !matches
            || self.problem.occupations != source.occupations
            || self.problem.electrons.to_bits() != source.electrons.to_bits()
            || self.problem.minimum_gap != source.minimum_gap
        {
            return Err(ScfError::InvalidInput(
                "warm_transfer: model differs beyond atom positions",
            ));
        }
        if density.shape() != self.problem.grid.dims()
            || !density.is_standard_layout()
            || density
                .iter()
                .any(|v| !v.is_finite() || *v < -xc::NEGATIVE_ROUNDOFF)
        {
            return Err(ScfError::InvalidInput(
                "warm_transfer: invalid density shape, C-order layout or values",
            ));
        }
        let charge = density.sum() * self.problem.cell.volume_bohr3() / density.len() as f64;
        if (charge - self.problem.electrons).abs() / self.problem.electrons > 1e-8 {
            return Err(ScfError::InvalidInput(
                "warm_transfer: relative charge error exceeds 1e-8",
            ));
        }
        self.initialize()?;
        self.continuation.as_mut().unwrap().n = density.clone();
        Ok(())
    }
    fn make_mixer(&self, mixing: &Mixing) -> Result<Box<dyn MixerHistory>, ScfError> {
        let mixer: Box<dyn MixerHistory> = match *mixing {
            Mixing::Linear { beta } => Box::new(LinearMixer { beta }),
            Mixing::Pulay { beta, history } => {
                let mut p = PulayMixer::new(beta, history)?;
                if self.problem.options.hartree_mixing {
                    p.set_metric(Box::new(CoulombMetric::new(
                        &self.problem.cell,
                        self.problem.options.cutoff_hartree,
                        self.problem.grid.dims(),
                    )?));
                }
                Box::new(p)
            }
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
                let mut p = PreconditionedPulay::new(beta, history, dims.iter().product(), filter)?;
                if self.problem.options.hartree_mixing {
                    p.set_metric(Box::new(CoulombMetric::new(
                        &self.problem.cell,
                        self.problem.options.cutoff_hartree,
                        dims,
                    )?));
                }
                Box::new(p)
            }
        };
        Ok(mixer)
    }
    fn initialize_inner(&mut self) -> Result<(), ScfError> {
        let start = Instant::now();
        let mut mixer = self.make_mixer(&self.problem.options.mixing)?;
        // Validate mixer parameters even if a future system converges immediately.
        let n = self.problem.starting_density();
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
        let standard = self.problem.options.standard;
        let mut tolerance = if standard.is_some() {
            let ordinary = history.last().map_or(StandardProfile::CAP_HA, |r| {
                StandardProfile::ratchet(
                    r.requested_eigen_tolerance,
                    r.standard_error_ha.expect("standard history metric"),
                    self.problem.electrons,
                )
            });
            if self.problem.options.standard_density_forcing {
                StandardProfile::density_forcing(ordinary, history)
            } else if self.problem.options.standard_stall_guard {
                StandardProfile::stall_guard(ordinary, history)
            } else {
                ordinary
            }
        } else if self.problem.options.residual_eigen_schedule {
            // Residual is the normalized density L2 norm in atomic units.
            // 0.1 is an empirical scheduling coefficient in atomic units.
            history
                .last()
                .map_or(1e-2, |r| schedule_factor() * r.density_residual)
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
        // B2: ordinary maps may take the nonlocal energy from the band-sum
        // identity; the flag is cleared again before any certificate work.
        let identity = crate::energy::energy_identity_requested();
        #[cfg(test)]
        let identity = identity || self.problem.force_energy_identity;
        self.problem.energy_identity = identity;
        let mut step = self
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
        let mut standard_error_ha = self
            .problem
            .coulomb
            .as_mut()
            .map(|m| m.error(n, &step.density))
            .transpose()?;
        let standard_redo = iteration == 1
            && standard_error_ha.is_some_and(|e| e < tolerance * self.problem.electrons.max(1.0));
        if standard_redo {
            tolerance = StandardProfile::ratchet(
                tolerance,
                standard_error_ha.unwrap(),
                self.problem.electrons,
            );
            let warm: Vec<_> = step
                .k_points
                .iter()
                .map(|k| k.solution.eigenvectors.clone())
                .collect();
            let mut redo =
                self.problem
                    .step_kpoints_with_certification(n, Some(&warm), tolerance, false)?;
            // Charge every discarded first solve, including its potential,
            // density and audit work. No mixer history exists before the redo.
            for (new, old) in redo.k_points.iter_mut().zip(&step.k_points) {
                new.solution.iterations += old.solution.iterations;
                new.solution.operator_applications += old.solution.operator_applications;
            }
            redo.timings.potential += step.timings.potential;
            redo.timings.diagonalization += step.timings.diagonalization;
            redo.timings.density += step.timings.density;
            redo.timings.energy += step.timings.energy;
            step = redo;
            standard_error_ha = Some(
                self.problem
                    .coulomb
                    .as_mut()
                    .unwrap()
                    .error(n, &step.density)?,
            );
        }
        if dft_solvers::replay::ENABLED {
            // Measurement only: signed density residual, normalized FFT, density
            // sphere and G=0 omission. Do not pass a signed field to the physical
            // Hartree API, which correctly rejects negative electron densities.
            assert_eq!(
                self.problem.options.threads, 1,
                "P15A requires ordered one-thread captures"
            );
            let dr2_ry = super::replay_capture::density_error_ry(
                n,
                &step.density,
                &self.problem.cell,
                self.problem.options.cutoff_hartree,
                &mut self.problem.fft,
            )?;
            dft_solvers::replay::metric(iteration, dr2_ry, self.problem.electrons);
        }
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
            standard_error_ha,
            standard_redo,
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
        let converged = if let Some(p) = standard {
            standard_error_ha.is_some_and(|e| e < p.threshold_ha)
        } else {
            converged(&record, &self.problem.options)
        };
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
        self.problem.energy_identity = false;
        let mut rescued = None;
        let mut unable_to_verify = false;
        // Opt-in pre-screen (BAANDI_CERT_SCREEN=on): predict a certain failure of the
        // occupied-eigenvalue criterion from first-order perturbation theory before
        // paying for the tight re-solve. It can only reject; passing still needs
        // the full solve below. `trace` prints the bound and runs the solve anyway.
        let screen_start = Instant::now();
        let screen = if converged && standard.is_none() {
            match cert_screen_mode() {
                CertScreen::Off => None,
                mode => self
                    .problem
                    .certificate_shift_estimate(&step.density, &step.k_points, potential_reuse())?
                    .map(|bound| (mode, bound)),
            }
        } else {
            None
        };
        let screen_seconds = screen_start.elapsed().as_secs_f64();
        let screened = screen.and_then(|(mode, bound)| {
            (mode == CertScreen::On && bound > CERT_SCREEN_MARGIN * CERT_EIGENVALUE_LIMIT_HA)
                .then_some(bound)
        });
        let certified = if let Some(bound) = screened {
            // Certain failure: the next SCF iteration continues exactly as after a
            // failed solve (that solve's result is discarded either way). The record
            // carries the step's own residuals, zero solver work, the first-order
            // eigenvalue-shift bound and the SCF density residual; zero Davidson
            // iterations marks it as screened, not solved.
            self.certificate_attempted = true;
            let k = step.k_points.len();
            history.last_mut().unwrap().fixed_point_certificate = Some(FixedPointCertificate {
                strict_residuals_per_k: step
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
                buffer_residuals_per_k: step
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
                density_residual: residual,
                occupied_eigenvalue_change: bound,
                davidson_iterations: 0,
                operator_applications: 0,
                davidson_iterations_per_k: vec![0; k],
                operator_applications_per_k: vec![0; k],
                timings: PhaseTimings::default(),
            });
            false
        } else if converged && standard.is_some() {
            // Standard convergence is explicitly uncertified. An independently
            // requested spectral audit remains a separate, charged operation.
            let mut passed = true;
            if let Some(report) = &mut self.spectral_report {
                (rescued, unable_to_verify) = self.problem.audit_reported_spectrum(
                    n,
                    &mut step.k_points,
                    iteration,
                    report,
                )?;
                passed = report.passed;
            }
            if passed {
                self.problem.check_gap(&step.k_points)?;
            }
            passed
        } else if converged {
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
            if let Some((_, bound)) = screen {
                eprintln!(
                    "CERT-SCREEN iter={iteration} bound={bound:.4e} actual={occupied_eigenvalue_change:.4e} density={density_residual:.4e} screen_s={screen_seconds:.6} davidson={}",
                    check
                        .k_points
                        .iter()
                        .map(|k| k.solution.iterations)
                        .sum::<usize>()
                );
            }
            let mut passed = density_residual <= self.problem.options.density_tolerance
                && occupied_eigenvalue_change <= CERT_EIGENVALUE_LIMIT_HA;
            if passed && self.problem.density_floor.is_some() {
                let clipped = self
                    .problem
                    .negative_charge(&step.density)
                    .max(self.problem.negative_charge(&check.density));
                if let Err(reason) =
                    ScfProblem::validate_clipped_charge(clipped, self.problem.electrons)
                {
                    return Err(ScfError::Diverged {
                        reason,
                        history: std::mem::take(history),
                    });
                }
            }
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
            if passed && let Some(report) = &mut self.spectral_report {
                (rescued, unable_to_verify) = self.problem.audit_reported_spectrum(
                    n,
                    &mut step.k_points,
                    iteration,
                    report,
                )?;
                passed = report.passed;
            }
            if passed {
                self.problem.check_gap(&check.k_points)?;
            }
            passed
        } else {
            false
        };
        // Fields kept by the screen belong to this certificate only.
        self.problem.clear_reused_fields();
        if certified
            || unable_to_verify
            || (retain_last && iteration == self.problem.options.max_iterations)
        {
            // B2: the final state always carries the explicit nonlocal audit.
            let energy_identity = if step.nonlocal_from_identity {
                let audit = self.problem.explicit_nonlocal_energy(&step.k_points)?;
                let identity = step.energy.nonlocal;
                step.energy.nonlocal = audit;
                step.nonlocal_from_identity = false;
                Some(crate::energy::EnergyIdentityCheck {
                    identity_nonlocal_ha: identity,
                    audit_nonlocal_ha: audit,
                    difference_ha: identity - audit,
                })
            } else {
                None
            };
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
                    atomic_start: self.problem.atomic_start.clone(),
                    compat_mode: self.problem.options.compat_mode,
                    crystal_density_group: self.problem.options.crystal_density_group,
                    standard_stall_guard: self.problem.options.standard_stall_guard,
                    standard_density_forcing: self.problem.options.standard_density_forcing,
                    pw92_prefactor: self.problem.options.pw92_prefactor,
                    local_mesh: self.problem.options.local_mesh,
                    radial_even_rule: self.problem.options.radial_even_rule,
                    gamma_real_used: self
                        .problem
                        .states
                        .iter()
                        .any(|s| s.gamma_real.gamma_real_used),
                    gamma_real_fallback_reason: self
                        .problem
                        .states
                        .iter()
                        .find_map(|s| s.gamma_real.fallback_reason.clone()),
                    force_identity: self.problem.force_model.with_policy(
                        crate::force_identity::OccupationIdentity::new(
                            &[&self.problem.occupations],
                            None,
                            None,
                        ),
                    ),
                    standard,
                    certified: certified && standard.is_none(),
                    spectral_completeness: self.spectral_report.clone().map(|mut r| {
                        r.passed &= certified;
                        r
                    }),
                    k_points: step.k_points,
                    gap_hartree,
                    energy: step.energy,
                    energy_identity,
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
        if rescued.is_some() {
            // The old history converged in an incomplete invariant subspace.
            // Its near-zero residual is not a meaningful growth baseline after
            // inserting the missing occupied states. Start a fresh mixing history.
            mixer.reset();
            *best = f64::INFINITY;
            *growing = 0;
            *mixing_resets += 1;
            let record = history.last_mut().unwrap();
            record.mixer_reset = true;
            record.mixing_resets = *mixing_resets;
        }
        // Experimental guard (BAANDI_MIX_GROWTH_RESET=<factor>, default off): a density
        // residual that grew by more than the factor since the previous iteration means
        // the Pulay history produced a bad step and is ill-conditioned; drop it.
        if let Some(limit) = growth_reset_limit() {
            let record = history.last_mut().unwrap();
            if record.density_contraction.is_some_and(|c| c > limit) && !record.mixer_reset {
                mixer.reset();
                *mixing_resets += 1;
                record.mixer_reset = true;
                record.mixing_resets = *mixing_resets;
            }
        }
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
        *previous = rescued.or_else(|| {
            Some(
                step.k_points
                    .into_iter()
                    .map(|k| k.solution.eigenvectors)
                    .collect(),
            )
        });
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
    /// C1 ladder hand-off: next input density and warm orbitals at an active boundary.
    pub(crate) fn warm_state(&self) -> Option<(&Array3<f64>, Option<&[ComplexMatrix]>)> {
        self.continuation
            .as_ref()
            .map(|c| (&c.n, c.previous.as_deref()))
    }
    pub(crate) fn problem_ref(&self) -> &ScfProblem {
        &self.problem
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

// Persistence is opt-in at completed boundaries. No call or branch is added to
// initialize/advance/step, including the Classic execution path.
impl ScfSession {
    fn check_checkpoint_model(
        &self,
        model: &crate::checkpoint::ModelIdentity,
    ) -> Result<(), crate::checkpoint::Rejection> {
        use crate::checkpoint::Rejection as R;
        model.validate()?;
        if self.problem.symmetry.is_some() || self.problem.options.crystal_density_group {
            return Err(R::Boundary);
        }
        // Checkpoint v1 has neither cumulative spectral audit records nor the
        // locking policy. Fail closed instead of dropping either on resume.
        if self.problem.options.spectral_guard == crate::SpectralGuard::Final
            || self.problem.options.locking != dft_solvers::locking::Locking::Off
        {
            return Err(R::Boundary);
        }
        if self.problem.core_density.is_some() || self.problem.density_floor.is_some() {
            return Err(R::Xc);
        }
        if model.cell_bohr != self.problem.cell.lattice_columns_bohr() {
            return Err(R::Cell);
        }
        if model.cutoff_hartree.to_bits() != self.problem.options.cutoff_hartree.to_bits() {
            return Err(R::Cutoff);
        }
        if model.electrons.to_bits() != self.problem.electrons.to_bits() {
            return Err(R::ElectronCount);
        }
        if model.xc != self.problem.options.functional.name() {
            return Err(R::Xc);
        }
        if model.spin != "unpolarized" {
            return Err(R::Spin);
        }
        let points: Vec<_> = self
            .problem
            .states
            .iter()
            .map(|s| {
                let k = s.point;
                [k.fractional[0], k.fractional[1], k.fractional[2], k.weight]
            })
            .collect();
        if model.kpoints != points {
            return Err(R::KMesh);
        }
        Ok(())
    }
    pub fn checkpoint(
        &self,
        kind: crate::checkpoint::CheckpointKind,
        model: crate::checkpoint::ModelIdentity,
        execution: crate::checkpoint::ExecutionIdentity,
    ) -> Result<crate::checkpoint::Checkpoint, crate::checkpoint::Rejection> {
        use crate::checkpoint::*;
        self.check_checkpoint_model(&model)?;
        if execution.threads != self.problem.options.threads || !is_digest(&execution.build_sha256)
        {
            return Err(Rejection::Identity);
        }
        let orbitals = |m: &ComplexMatrix| Orbitals {
            rows: m.rows(),
            cols: m.cols(),
            coefficients: m.as_slice().iter().flat_map(|z| [z.re, z.im]).collect(),
        };
        let density = self.density().ok_or(Rejection::Boundary)?;
        let previous = if let Some(o) = &self.outcome {
            o.result
                .k_points
                .iter()
                .map(|k| orbitals(&k.solution.eigenvectors))
                .collect()
        } else {
            self.continuation
                .as_ref()
                .and_then(|c| c.previous.as_ref())
                .map(|v| v.iter().map(orbitals).collect())
                .unwrap_or_default()
        };
        let sequence = if kind == CheckpointKind::ExactContinuation {
            self.events
                .as_ref()
                .map(EventQueue::checkpoint_sequence)
                .transpose()?
                .unwrap_or(0)
        } else {
            self.snapshot().sequence
        };
        let continuation = if kind == CheckpointKind::ExactContinuation {
            if self.status != StepStatus::Ready {
                return Err(Rejection::Boundary);
            }
            let c = self.continuation.as_ref().ok_or(Rejection::Boundary)?;
            Some(ContinuationState {
                atomic_start: (self.problem.options.start_wavefunctions
                    == StartWavefunctions::Atomic)
                    .then(|| self.problem.atomic_start.clone()),
                options: self.problem.options.clone(),
                completion: self.completion,
                config: self.config.clone(),
                iteration: self.iteration,
                state_version: self.state_version,
                sequence,
                mixer: c.mixer.checkpoint().map_err(|_| Rejection::Mixer)?,
                history: c.history.clone(),
                best: c.best.is_finite().then_some(c.best),
                growing: c.growing,
                mixing_resets: c.mixing_resets,
                minimum_gap: self.problem.minimum_gap,
            })
        } else {
            None
        };
        let result = if kind == CheckpointKind::ResultExport {
            let o = self.outcome.as_ref().ok_or(Rejection::Boundary)?;
            let r = &o.result;
            Some(ExportResult {
                converged: o.converged,
                energy: r.energy.clone(),
                charge: r.charge,
                gap_hartree: r.gap_hartree,
                eigenvalues: r
                    .k_points
                    .iter()
                    .map(|k| k.solution.eigenvalues.clone())
                    .collect(),
                residuals: r
                    .k_points
                    .iter()
                    .map(|k| k.solution.residuals.clone())
                    .collect(),
                history: r.history.clone(),
            })
        } else {
            None
        };
        Ok(Checkpoint {
            standard_stall_guard: self.problem.options.standard_stall_guard,
            standard_density_forcing: self.problem.options.standard_density_forcing,
            standard: self.problem.options.standard,
            hartree_mixing: self.problem.options.hartree_mixing,
            schema_version: SCHEMA_VERSION,
            kind,
            calc_id: self.config.calc_id.clone(),
            sequence,
            timestamp: std::time::SystemTime::now(),
            model,
            execution,
            electronic_state: ElectronicState {
                grid: self.problem.grid.dims(),
                density: density.iter().copied().collect(),
                orbitals: previous,
                occupations: self.problem.occupations.clone(),
            },
            continuation,
            runtime: None,
            result,
        })
    }
    /// Track A hand-off hook (direct minimization pre-converger, `scf::direct_min`):
    /// ordinary initialization, then the same density/orbital validation and import
    /// as `warm_start`. Needed because `warm_start` requires host identity digests
    /// and rejects NLCC, and `initialize_density_transfer` drops the orbitals.
    /// Nothing is certified here; only the unchanged strict steps can converge.
    pub(crate) fn initialize_electronic(
        &mut self,
        state: &crate::checkpoint::ElectronicState,
    ) -> Result<(), ScfError> {
        self.initialize()?;
        self.import_electronic(state)
            .map_err(|e| ScfError::Numerical(format!("direct_min hand-off state rejected: {e}")))
    }
    fn import_electronic(
        &mut self,
        state: &crate::checkpoint::ElectronicState,
    ) -> Result<(), crate::checkpoint::Rejection> {
        use crate::checkpoint::Rejection as R;
        let dims = self.problem.grid.dims();
        if state.grid != dims || state.density.len() != dims.iter().product::<usize>() {
            return Err(R::Grid);
        }
        if state.occupations != self.problem.occupations {
            return Err(R::Occupations);
        }
        if state
            .density
            .iter()
            .any(|n| !n.is_finite() || *n < -xc::NEGATIVE_ROUNDOFF)
        {
            return Err(R::Density);
        }
        let n = Array3::from_shape_vec(dims, state.density.clone()).map_err(|_| R::Grid)?;
        let charge = n.sum() * (self.problem.cell.volume_bohr3() / n.len() as f64);
        if (charge - self.problem.electrons).abs() / self.problem.electrons > 1e-8 {
            return Err(R::Density);
        }
        let previous = if state.orbitals.is_empty() {
            None
        } else {
            if state.orbitals.len() != self.problem.states.len() {
                return Err(R::Orbitals);
            }
            let mut previous = Vec::new();
            for (o, k) in state.orbitals.iter().zip(&self.problem.states) {
                if o.rows != k.basis.vectors().len()
                    || o.cols != self.problem.occupations.len()
                    || o.rows.checked_mul(o.cols).and_then(|n| n.checked_mul(2))
                        != Some(o.coefficients.len())
                    || o.coefficients.iter().any(|x| !x.is_finite())
                {
                    return Err(R::Orbitals);
                }
                let m = ComplexMatrix::new(
                    o.rows,
                    o.cols,
                    o.coefficients
                        .as_chunks::<2>()
                        .0
                        .iter()
                        .map(|z| Complex64::new(z[0], z[1]))
                        .collect(),
                )
                .map_err(|_| R::Orbitals)?;
                for j in 0..m.cols() {
                    for i in 0..=j {
                        let overlap: Complex64 = m
                            .column(i)
                            .iter()
                            .zip(m.column(j))
                            .map(|(a, b)| a.conj() * b)
                            .sum();
                        if (overlap - Complex64::new(if i == j { 1.0 } else { 0.0 }, 0.0)).norm()
                            > 1e-8
                        {
                            return Err(R::Orbitals);
                        }
                    }
                }
                previous.push(m);
            }
            Some(previous)
        };
        let c = self.continuation.as_mut().ok_or(R::Boundary)?;
        c.n = n;
        c.previous = previous;
        Ok(())
    }
    pub fn restore_exact(
        problem: ScfProblem,
        checkpoint: &crate::checkpoint::Checkpoint,
        model: &crate::checkpoint::ModelIdentity,
        execution: &crate::checkpoint::ExecutionIdentity,
    ) -> Result<Self, crate::checkpoint::Rejection> {
        use crate::checkpoint::*;
        if checkpoint.runtime.is_some() {
            return Err(Rejection::Runtime);
        }
        Self::restore_exact_inner(problem, checkpoint, model, execution)
    }
    pub(crate) fn restore_exact_inner(
        problem: ScfProblem,
        checkpoint: &crate::checkpoint::Checkpoint,
        model: &crate::checkpoint::ModelIdentity,
        execution: &crate::checkpoint::ExecutionIdentity,
    ) -> Result<Self, crate::checkpoint::Rejection> {
        use crate::checkpoint::*;
        checkpoint.compatible(CheckpointKind::ExactContinuation, model, execution)?;
        let c = checkpoint.continuation.as_ref().ok_or(Rejection::Kind)?;
        if checkpoint.result.is_some() {
            return Err(Rejection::Kind);
        }
        if checkpoint.standard != c.options.standard
            || checkpoint.standard_stall_guard != c.options.standard_stall_guard
            || checkpoint.standard_density_forcing != c.options.standard_density_forcing
            || checkpoint.hartree_mixing != c.options.hartree_mixing
        {
            return Err(Rejection::Options);
        }
        if !options_identical(&problem.options, &c.options) || problem.minimum_gap != c.minimum_gap
        {
            return Err(Rejection::Options);
        }
        if c.config.calc_id != checkpoint.calc_id
            || c.sequence != checkpoint.sequence
            || c.iteration
                .checked_add(1)
                .is_none_or(|n| c.state_version < n as u64)
            || c.state_version >= u64::MAX - 1
            || c.sequence >= u64::MAX - 16
            || execution.threads != c.options.threads
        {
            return Err(Rejection::Identity);
        }
        if c.iteration >= c.options.max_iterations
            || c.history.len() != c.iteration
            || c.history.iter().enumerate().any(|(i, r)| {
                r.iteration != i + 1
                    || !r.energy.total_hartree().is_finite()
                    || !r.density_residual.is_finite()
                    || r.standard_error_ha.is_some() != c.options.standard.is_some()
                    || r.standard_error_ha
                        .is_some_and(|v| !v.is_finite() || v < 0.0)
                    || (c.options.standard.is_some()
                        && (!r.requested_eigen_tolerance.is_finite()
                            || !(StandardProfile::FLOOR_HA..=StandardProfile::CAP_HA)
                                .contains(&r.requested_eigen_tolerance)
                            || r.fixed_point_certificate.is_some()))
            })
            || (c.iteration > 0 && checkpoint.electronic_state.orbitals.is_empty())
            || c.best.is_some_and(|v| !v.is_finite() || v < 0.0)
            || (c.iteration > 0 && c.best.is_none())
            || c.growing > c.iteration
            || c.growing >= c.options.divergence_patience
            || c.mixing_resets as u64 > c.state_version
            || c.best.map(f64::to_bits)
                != c.history
                    .iter()
                    .map(|r| r.density_residual)
                    .reduce(f64::min)
                    .map(f64::to_bits)
        {
            return Err(Rejection::History);
        }
        if let MixerState::Kerker { diagnostics, .. } = &c.mixer
            && (diagnostics.raw_recoveries > c.iteration
                || diagnostics.resets as u64 > c.state_version)
        {
            return Err(Rejection::Mixer);
        }
        let mut session =
            Self::new(problem, c.completion, c.config.clone()).map_err(|_| Rejection::Options)?;
        session.check_checkpoint_model(model)?;
        session.initialize().map_err(|_| Rejection::Mixer)?;
        session.import_electronic(&checkpoint.electronic_state)?;
        session.problem.atomic_start = c.atomic_start.clone().unwrap_or_default();
        let state = session.continuation.as_mut().unwrap();
        state
            .mixer
            .restore(&c.mixer, state.n.len())
            .map_err(|_| Rejection::Mixer)?;
        state.history = c.history.clone();
        state.best = c.best.unwrap_or(f64::INFINITY);
        state.growing = c.growing;
        state.mixing_resets = c.mixing_resets;
        session.iteration = c.iteration;
        session.state_version = c.state_version;
        // initialize emitted an event; restored streams start after the saved,
        // durably drained cursor and must not fabricate a second initialization.
        session.events =
            (c.config.mode != ExecutionMode::Classic).then(|| EventQueue::new(&c.config));
        if let Some(q) = &mut session.events {
            q.restore_sequence(c.sequence);
        }
        Ok(session)
    }
    pub fn warm_start(
        problem: ScfProblem,
        checkpoint: &crate::checkpoint::Checkpoint,
        model: &crate::checkpoint::ModelIdentity,
        execution: &crate::checkpoint::ExecutionIdentity,
        completion: CompletionMode,
        config: ExecutionConfig,
    ) -> Result<Self, crate::checkpoint::Rejection> {
        use crate::checkpoint::*;
        checkpoint.compatible(CheckpointKind::WarmStart, model, execution)?;
        if checkpoint.continuation.is_some()
            || checkpoint.runtime.is_some()
            || checkpoint.result.is_some()
        {
            return Err(Rejection::Kind);
        }
        if config.calc_id == checkpoint.calc_id {
            return Err(Rejection::Identity);
        }
        let mut session = Self::new(problem, completion, config).map_err(|_| Rejection::Options)?;
        session.check_checkpoint_model(model)?;
        session.initialize().map_err(|_| Rejection::Mixer)?;
        session.import_electronic(&checkpoint.electronic_state)?;
        Ok(session)
    }
}

/// Largest occupied-eigenvalue change (Ha) the fixed-point certificate accepts.
const CERT_EIGENVALUE_LIMIT_HA: f64 = 1e-7;
/// The pre-screen rejects only above this multiple of the limit, so its
/// second-order modelling error (about 1e-9 Ha) is far inside the margin.
const CERT_SCREEN_MARGIN: f64 = 1.02;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum CertScreen {
    Off,
    On,
    Trace,
}

/// Certificate screen: default ON since the user's adopt-what-improves rule
/// (2026-11-10; dev08 ledger SINGLE-NODE-AB1: bit-identical density and energy in
/// 40 of 40 pairs, -19 percent wall on si32 strict). `BAANDI_CERT_SCREEN=off`
/// restores the unscreened path, `trace` prints the estimates.
fn cert_screen_mode() -> CertScreen {
    static MODE: std::sync::OnceLock<CertScreen> = std::sync::OnceLock::new();
    *MODE.get_or_init(|| match std::env::var("BAANDI_CERT_SCREEN").as_deref() {
        Ok("off") => CertScreen::Off,
        Ok("trace") => CertScreen::Trace,
        _ => CertScreen::On,
    })
}

/// POTENTIAL-REUSE: with the screen active (`on` or `trace`), the Hartree/XC
/// fields it builds for the output density are reused by the certificate solve
/// at that same density. Bit-identical either way; `BAANDI_POTENTIAL_REUSE=off`
/// disables the reuse for A/B checks. Without the screen nothing is kept.
fn potential_reuse() -> bool {
    static MODE: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *MODE.get_or_init(|| std::env::var("BAANDI_POTENTIAL_REUSE").as_deref() != Ok("off"))
}

fn growth_reset_limit() -> Option<f64> {
    static LIMIT: std::sync::OnceLock<Option<f64>> = std::sync::OnceLock::new();
    *LIMIT.get_or_init(|| {
        std::env::var("BAANDI_MIX_GROWTH_RESET")
            .ok()
            .and_then(|v| v.parse::<f64>().ok())
            .filter(|v| v.is_finite() && *v > 1.0)
    })
}

/// Coefficient of the residual-driven eigen tolerance schedule (tolerance =
/// factor x density residual). Default 0.1 (the established empirical value);
/// `BAANDI_EIGEN_FACTOR=<x>` (0.001 to 1) is an experiment switch for dev08 A/B.
fn schedule_factor() -> f64 {
    static FACTOR: std::sync::OnceLock<f64> = std::sync::OnceLock::new();
    *FACTOR.get_or_init(|| {
        std::env::var("BAANDI_EIGEN_FACTOR")
            .ok()
            .and_then(|v| v.parse::<f64>().ok())
            .filter(|v| v.is_finite() && (1e-3..=1.0).contains(v))
            .unwrap_or(0.1)
    })
}
