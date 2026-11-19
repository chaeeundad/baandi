//! Incremental execution of the existing SCF iteration, with immutable boundaries.
use super::*;
use crate::events::{EventKind, EventQueue, ExecutionConfig, ExecutionMode, FailureCode};
use crate::policy::{PolicyChange, ToleranceFloors, prepare_policy};
use crate::runtime::ControlError;
use crate::state_view::StateSnapshot;
use dft_solvers::mixing::DensityMixer as _;

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
    /// Next input density; an EMPTY array in the owned density loop
    /// (`owned`), whose input lives in `problem.owned.input`.
    n: Array3<f64>,
    mixer: Box<dyn MixerHistory>,
    /// M2-OWNED (`BAANDI_M2_LOCAL_DENSITY`): the plane-mode distributed Pulay
    /// mixer of the owned density loop (then `mixer` is not used).
    owned: Option<crate::mix_dist::DistributedPulay<crate::parallel::SharedComm>>,
    previous: Option<Vec<ComplexMatrix>>,
    /// MPI-LAZY-ORBITALS: distributed rows of the blocks of `previous` that are
    /// 0 x B placeholders (same order; empty when none is lazy).
    previous_lazy: Vec<Option<crate::gamma_real::LazyOrbitals>>,
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
    /// Standard profile only: eigen tolerance of the first step instead of the
    /// cap (cutoff-ladder hand-off; set before the first step, crate-private).
    standard_seed: Option<f64>,
    /// Strict residual schedule only: eigen tolerance of the first step instead
    /// of the schedule's 1e-2 start (cutoff-ladder hand-off, D1
    /// `BAANDI_LADDER_FIRST_TOL`; set before the first step, crate-private).
    strict_seed: Option<f64>,
    /// ADR 0043: attempts of the post-convergence verification.
    verification: VerificationReport,
    /// BREAKTHROUGH-3 (`BAANDI_SCF_ADVANCE=on`): correction blocks of the next
    /// advanced map (adapted from the previous advanced map's inner/SCF ratio).
    advance_blocks: usize,
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
            standard_seed: None,
            strict_seed: None,
            verification: VerificationReport::default(),
            advance_blocks: crate::advance::starting_blocks(),
        })
    }

    /// BREAKTHROUGH-3: whether the coming ordinary map may run as a bounded
    /// subspace update (design `docs/design/bounded-subspace-scf.md`, section
    /// 5): switch on, standard path, a history exists (iteration 1 and its
    /// redo initialize the subspace exactly), the previous map's SCF residual
    /// is at or above the threshold (otherwise the stop rule is in reach and
    /// the map must be exact), one k point, closed-shell fixed occupations
    /// (equal occupation of every occupied band, so P = X X^T is the density
    /// matrix) and no runtime solver override.
    fn advance_eligible(problem: &ScfProblem, history: &[IterationRecord]) -> bool {
        let Some(profile) = problem.options.standard else {
            return false;
        };
        if !crate::advance::enabled() || problem.states.len() != 1 {
            return false;
        }
        if crate::advance::fine_only() && !dft_numerics::work_timers::in_fine_stage() {
            return false;
        }
        let Some(last) = history.last() else {
            return false;
        };
        if last
            .standard_error_ha
            .is_none_or(|e| e < profile.threshold_ha)
        {
            return false;
        }
        let occupied: Vec<f64> = problem
            .occupations
            .iter()
            .copied()
            .filter(|&f| f > 0.0)
            .collect();
        !occupied.is_empty()
            && occupied.iter().all(|&f| f == occupied[0])
            && matches!(
                problem.options.eigensolver,
                dft_solvers::SolverChoice::Davidson
            )
    }
    pub fn initialize(&mut self) -> Result<(), ScfError> {
        self.initialize_mode(false)
    }
    fn initialize_mode(&mut self, handoff: bool) -> Result<(), ScfError> {
        // RUN-CONFIG: legacy switch reads below see the problem's configuration.
        let _run = self.problem.run.enter();
        if self.continuation.is_some() || self.status != StepStatus::Created {
            return Err(ScfError::InvalidInput(
                "session already initialized or failed",
            ));
        }
        match self.initialize_inner(handoff) {
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
        self.set_next_input(density.clone())
    }
    /// SCF-TRAJECTORY Coulomb guard of the distributed mixer
    /// (`BAANDI_MIX_COULOMB_GUARD=trace|on`): a slab plan of the M2 grid with
    /// the Coulomb-metric weights (the plane layout of the owned density; only
    /// its Coulomb part is used) and the refit flag. None when off.
    /// Collective (plan construction).
    fn guard_plan(&self) -> Result<Option<(crate::grid_dist::GridPlan, bool)>, ScfError> {
        let guard = dft_numerics::switches::choice(dft_numerics::switches::Id::MIX_COULOMB_GUARD);
        if guard == "off" {
            return Ok(None);
        }
        let comm = crate::m2::active().ok_or(ScfError::InvalidInput(
            "distributed Coulomb guard without the M2 path",
        ))?;
        let dims = self.problem.grid.dims();
        let built;
        let weights = match self.problem.coulomb.as_ref() {
            Some(c) => c.weights(),
            None => {
                built = CoulombMetric::new(&self.problem.cell, &self.problem.options, dims)?;
                built.weights()
            }
        };
        let plan = crate::grid_dist::GridPlan::new(
            &comm,
            dims,
            self.problem.cell.reciprocal_columns_inverse_bohr(),
            self.problem.cell.volume_bohr3(),
            crate::grid_dist::Lda::Pz,
            Some(weights),
        )?;
        Ok(Some((plan, guard == "on")))
    }
    fn make_mixer(&self, mixing: &Mixing) -> Result<Box<dyn MixerHistory>, ScfError> {
        let mixer: Box<dyn MixerHistory> = match *mixing {
            Mixing::Linear { beta } => Box::new(LinearMixer { beta }),
            // M2-GRID (`BAANDI_M2_MIX=on`, default off): plain Pulay on the M2
            // path keeps its window on the rank's grid block.
            Mixing::Pulay { beta, history }
                if !self.problem.options.hartree_mixing
                    && self.problem.owned.declined
                        != Some(crate::grid_dist::OwnedDecline::Checkpoint)
                    && crate::mix_dist::requested()
                    && crate::m2::active().is_some() =>
            {
                let mut m = crate::mix_dist::DistributedPulay::new(
                    crate::m2::active().expect("checked"),
                    beta,
                    history,
                    self.problem.grid.dims().iter().product(),
                )?;
                if let Some((plan, refit)) = self.guard_plan()? {
                    m.set_coulomb_guard(plan, refit);
                }
                Box::new(m)
            }
            Mixing::Pulay { beta, history } => {
                let mut p = PulayMixer::new(beta, history)?;
                if self.problem.options.hartree_mixing {
                    p.set_metric(Box::new(CoulombMetric::new(
                        &self.problem.cell,
                        &self.problem.options,
                        self.problem.grid.dims(),
                    )?));
                } else {
                    // SCF-TRAJECTORY (`BAANDI_MIX_COULOMB_GUARD`, default off).
                    let guard = dft_numerics::switches::choice(
                        dft_numerics::switches::Id::MIX_COULOMB_GUARD,
                    );
                    if guard != "off" {
                        p.set_coulomb_guard(
                            Box::new(CoulombMetric::new(
                                &self.problem.cell,
                                &self.problem.options,
                                self.problem.grid.dims(),
                            )?),
                            guard == "on",
                        );
                    }
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
                        &self.problem.options,
                        dims,
                    )?));
                }
                Box::new(p)
            }
        };
        Ok(mixer)
    }
    fn initialize_inner(&mut self, handoff: bool) -> Result<(), ScfError> {
        let start = Instant::now();
        // M2-OWNED: an agent-controlled run may replace the mixer or read the
        // density at an active boundary; the owned modes decline for it.
        if self.config.mode == crate::events::ExecutionMode::AgentControlled {
            self.problem
                .decline_owned(crate::grid_dist::OwnedDecline::RuntimeControl);
        }
        let mut mixer = self.make_mixer(&self.problem.options.mixing)?;
        // Validate mixer parameters even if a future system converges immediately.
        let n = self.problem.starting_density();
        if handoff {
            // The initial trial is discarded. Validate the same density and
            // Pulay parameters locally instead of communicating three times.
            // Eligibility is restricted to plain Pulay for this constructor.
            let Mixing::Pulay { beta, history } = self.problem.options.mixing else {
                return Err(ScfError::InvalidInput("owned handoff requires Pulay"));
            };
            PulayMixer::new(beta, history)?.mix(
                n.as_slice().unwrap(),
                n.as_slice().unwrap(),
                self.problem.electrons / self.problem.cell.volume_bohr3(),
            )?;
        } else {
            mixer.mix(
                n.as_slice().unwrap(),
                n.as_slice().unwrap(),
                self.problem.electrons / self.problem.cell.volume_bohr3(),
            )?;
        }
        mixer.reset();
        let weight = self.problem.cell.volume_bohr3() / n.len() as f64;
        // M2-OWNED (`BAANDI_M2_LOCAL_DENSITY`, default off): the density stays
        // on the slab planes through the audit, the mixer and the next
        // potential; static replicated eligibility (no spectral audit here).
        let owned = if self.spectral_report.is_none() && self.problem.owned_density_eligible() {
            let Mixing::Pulay { beta, history } = self.problem.options.mixing else {
                unreachable!("eligibility requires plain Pulay");
            };
            let input = self.problem.owned_from_full(&n)?;
            let total = n.len();
            let (comm, plan) = self.problem.grid_plan()?;
            let (len, counts) = (plan.plane_len(), plan.plane_counts().to_vec());
            let mut m = crate::mix_dist::DistributedPulay::new_planes(
                comm, beta, history, total, len, &counts,
            )?;
            if let Some((plan, refit)) = self.guard_plan()? {
                m.set_coulomb_guard(plan, refit);
            }
            self.problem.owned.input = Some(input);
            self.problem.owned.density_mode = true;
            Some(m)
        } else {
            self.problem.owned.input = None;
            self.problem.owned.density_mode = false;
            None
        };
        let n = if owned.is_some() {
            Array3::zeros((0, 0, 0))
        } else {
            n
        };
        self.continuation = Some(Continuation {
            n,
            mixer,
            owned,
            previous: None,
            previous_lazy: Vec::new(),
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
            owned,
            previous,
            previous_lazy,
            history,
            weight,
            best,
            growing,
            mixing_resets,
            start,
        } = self.continuation.as_mut().unwrap();
        let weight = *weight;
        if owned.is_some() && !n.is_empty() {
            // A materialized copy (`materialize_density`) of the owned input:
            // the owned planes stay authoritative; the copy is dropped so that
            // it cannot go stale.
            *n = Array3::zeros((0, 0, 0));
        }
        let invalid = if owned.is_some() {
            self.problem.owned_input_invalid()?
        } else {
            n.iter()
                .any(|n| !n.is_finite() || *n < -xc::NEGATIVE_ROUNDOFF)
        };
        if invalid {
            return Err(ScfError::Diverged {
                reason: "mixed density is nonfinite or negative beyond 1e-14 Bohr^-3".into(),
                history: std::mem::take(history),
            });
        }
        let standard = self.problem.options.standard;
        let mut tolerance = if standard.is_some() {
            let first = self.standard_seed.unwrap_or(StandardProfile::CAP_HA);
            let ordinary = history.last().map_or(first, |r| {
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
            residual_schedule(
                &self.problem.options,
                self.problem.sw.eigen_factor,
                history.last().map(|r| r.density_residual),
                self.strict_seed.unwrap_or(1e-2),
            )
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
        // BREAKTHROUGH-0-INSTRUMENT: context of the map (snapshot labels,
        // first-block records) and its phase-profile bucket.
        let context = crate::profile::context();
        crate::profile::set_context(crate::profile::Context {
            kind: if previous.is_none() {
                "cold"
            } else {
                "ordinary"
            },
            iteration,
            ..context
        });
        // BREAKTHROUGH-3: request a bounded subspace update of this map.
        let advance_requested = Self::advance_eligible(&self.problem, history);
        if advance_requested {
            // The block count decides; the column cap only bounds a start
            // that is not the kept one (its initial application is one block).
            let bands = self.problem.occupations.len();
            self.problem.states[0].gamma_real.advance = Some(crate::advance::AdvanceSlot::new(
                dft_solvers::real_davidson::SubspaceBudget {
                    correction_blocks: self.advance_blocks,
                    h_columns: (self.advance_blocks + 1).saturating_mul(bands).max(1),
                },
            ));
        }
        let mut step = self
            .problem
            .step_session(
                n,
                previous.as_deref(),
                previous_lazy,
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
        let stop_clock = crate::profile::clock();
        let mut standard_error_ha =
            self.problem
                .stop_error_step(n, &step.density, step.owned_density.as_ref())?;
        crate::profile::between(crate::profile::BetweenPart::StopRule, stop_clock);
        crate::profile::map_commit(crate::profile::MapKind::Ordinary);
        // BREAKTHROUGH-3: completion rules of an advanced map (design section 5).
        let mut advance_record = None;
        match step.advance.take() {
            Some(slot) => match slot.outcome {
                Some(diagnostics) => {
                    let error = standard_error_ha.expect("standard history metric");
                    let inner = StandardProfile::residual_energy_estimate(
                        step.k_points
                            .iter()
                            .map(|k| (k.point.weight, k.solution.residuals.as_slice())),
                        &self.problem.occupations,
                    );
                    let occupied = self
                        .problem
                        .occupations
                        .iter()
                        .filter(|&&f| f > 0.0)
                        .count();
                    let values = &step.k_points[0].solution.eigenvalues;
                    let gap_stable = occupied == 0
                        || occupied >= values.len()
                        || values[occupied - 1] < values[occupied];
                    let threshold = standard.expect("standard profile").threshold_ha;
                    let previous_error = history.last().and_then(|r| r.standard_error_ha);
                    let trigger = crate::advance::completion_trigger(
                        error,
                        inner,
                        threshold,
                        previous_error,
                        gap_stable,
                    );
                    let mut record = crate::advance::AdvanceRecord {
                        budget_blocks: slot.budget.correction_blocks,
                        blocks_used: diagnostics.correction_blocks,
                        h_columns: diagnostics.h_columns,
                        kept_start: diagnostics.kept_start,
                        unconverged_roots: diagnostics.unconverged_roots,
                        inner_estimate_ha: inner,
                        scf_error_ha: error,
                        completion: trigger,
                        completion_h_columns: 0,
                    };
                    if trigger.is_some() {
                        // Completion: the exact solve at the same density from
                        // the advanced orbitals (the kept start applies, the
                        // potential is unchanged). The advance is charged to
                        // the step like a first-iteration redo.
                        let (warm, warm_lazy) = split_warm(Some(&step.k_points));
                        let mut redo = self
                            .problem
                            .step_session(n, warm.as_deref(), &warm_lazy, tolerance, false)
                            .map_err(|e| match e {
                                ScfError::NonFinite => ScfError::Diverged {
                                    reason: "nonfinite SCF arithmetic".into(),
                                    history: history.clone(),
                                },
                                other => other,
                            })?;
                        record.completion_h_columns = redo
                            .k_points
                            .iter()
                            .map(|k| k.solution.operator_applications)
                            .sum();
                        for (new, old) in redo.k_points.iter_mut().zip(&step.k_points) {
                            new.solution.iterations += old.solution.iterations;
                            new.solution.operator_applications +=
                                old.solution.operator_applications;
                        }
                        redo.timings.potential += step.timings.potential;
                        redo.timings.diagonalization += step.timings.diagonalization;
                        redo.timings.density += step.timings.density;
                        redo.timings.energy += step.timings.energy;
                        step = redo;
                        standard_error_ha = self.problem.stop_error_step(
                            n,
                            &step.density,
                            step.owned_density.as_ref(),
                        )?;
                        // Phase profile: the completion solve is filed as a redo.
                        crate::profile::map_commit(crate::profile::MapKind::Redo);
                    }
                    let ratio = if error > 0.0 {
                        inner / error
                    } else {
                        f64::INFINITY
                    };
                    self.advance_blocks =
                        crate::advance::adapt_blocks(self.advance_blocks, ratio, trigger);
                    crate::advance::count_record(&record, &diagnostics);
                    advance_record = Some(record);
                }
                None => crate::advance::count_declined(),
            },
            None if crate::advance::enabled() && standard.is_some() => {
                crate::advance::count_exact_map();
            }
            None => {}
        }
        let standard_redo = iteration == 1
            && standard_error_ha.is_some_and(|e| e < tolerance * self.problem.electrons.max(1.0));
        if standard_redo {
            tolerance = StandardProfile::ratchet(
                tolerance,
                standard_error_ha.unwrap(),
                self.problem.electrons,
            );
            let (warm, warm_lazy) = split_warm(Some(&step.k_points));
            crate::profile::set_context(crate::profile::Context {
                kind: "redo",
                iteration,
                ..context
            });
            let mut redo =
                self.problem
                    .step_session(n, warm.as_deref(), &warm_lazy, tolerance, false)?;
            crate::profile::map_commit(crate::profile::MapKind::Redo);
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
            standard_error_ha =
                self.problem
                    .stop_error_step(n, &step.density, step.owned_density.as_ref())?;
        }
        let stop_clock = crate::profile::clock();
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
        // M2-OWNED: owned pairs reduce their per-plane partials (one
        // allgatherv); full pairs keep the row-major arithmetic.
        let (residual, output_sum, input_sum, min_density) =
            self.problem
                .step_metrics(n, &step.density, step.owned_density.as_ref())?;
        let charge = output_sum * weight;
        let charge_error = (charge - self.problem.electrons).abs() / self.problem.electrons;
        let input_charge_error =
            (input_sum * weight - self.problem.electrons).abs() / self.problem.electrons;
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
            min_density,
            timings: step.timings,
            advance: advance_record,
            pulay: None,
        };
        let converged = if let Some(p) = standard {
            let stop = standard_error_ha.is_some_and(|e| {
                e < p.threshold_ha
                    || (self.problem.sw.standard_stop_confirm
                        && history
                            .last()
                            .and_then(|r| r.standard_error_ha)
                            .is_some_and(|previous| {
                                // Confirmed stop (offline fit on 60 dev08 trajectories, ledger
                                // STANDARD-STOP): within 2x the threshold, the previous step
                                // within 20x, and a contraction of at least 3.3x.
                                e < 2.0 * p.threshold_ha
                                    && previous <= 20.0 * p.threshold_ha
                                    && e <= 0.3 * previous
                            }))
            });
            // Residual guard (default on; StandardProfile::residual_energy_estimate,
            // an empirical estimator, not a bound: review F03): the Coulomb
            // estimator does not see the second-order band-energy error of loose
            // eigensolves, so also require that estimate below the class
            // threshold. Only evaluated when the density rule already stops.
            if stop && self.problem.sw.residual_guard {
                let bound = StandardProfile::residual_energy_estimate(
                    step.k_points
                        .iter()
                        .map(|k| (k.point.weight, k.solution.residuals.as_slice())),
                    &self.problem.occupations,
                );
                StandardProfile::residual_guard(stop, bound, p.threshold_ha)
            } else {
                stop
            }
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
        crate::profile::between(crate::profile::BetweenPart::StopRule, stop_clock);
        // Certify the returned output density itself, not just the input
        // fixed-point residual. A small-gap response can amplify that residual.
        self.problem.energy_identity = false;
        let mut rescued = None;
        // VERIFY-REUSE (`BAANDI_VERIFY_REUSE=warm|feed`): orbitals and output
        // density F(rho_out) of a failed verification solve.
        let mut verification_reuse: Option<(VerifyReuse, VerificationMap)> = None;
        let mut unable_to_verify = false;
        // Opt-in pre-screen (BAANDI_CERT_SCREEN=on): predict a certain failure of the
        // occupied-eigenvalue criterion from first-order perturbation theory before
        // paying for the tight re-solve. It can only reject; passing still needs
        // the full solve below. `trace` prints the bound and runs the solve anyway.
        // Scheduling only vetoes an attempt. A deferred boundary continues
        // ordinary mixing; only the unchanged verification below can pass.
        let defer_verification = if converged
            && standard.is_some()
            && self.problem.options.verification == Verification::FixedPoint
            && self.problem.sw.verify_schedule != "off"
            && self.spectral_report.is_none()
        {
            let prediction = verification_schedule(history, &self.verification.attempts);
            let wait = prediction.extra_maps > 0 && iteration < self.problem.options.max_iterations;
            if self.problem.sw.verify_schedule == "trace" {
                eprintln!(
                    "VERIFY-SCHEDULE map={iteration} samples={} failure_vote={:.6} extra_maps={} would_wait={wait}",
                    prediction.samples, prediction.failure_vote, prediction.extra_maps,
                );
            }
            self.problem.sw.verify_schedule == "on" && wait
        } else {
            false
        };
        let certificate_clock = Instant::now();
        // P0.2: work of the certificate block is reported as the verification
        // run stage (a coarse-stage certificate stays in the coarse stage).
        let verification = dft_numerics::work_timers::in_fine_stage().then(|| {
            dft_numerics::work_timers::enter(dft_numerics::work_timers::RunStage::Verification)
        });
        // Phase profile: the certificate block is its own top-level bucket.
        let profile_verification = converged.then(|| {
            (
                crate::profile::enter(crate::profile::Top::Verification),
                crate::profile::context_stage("verification"),
            )
        });
        let screen_start = Instant::now();
        let screen = if converged && standard.is_none() {
            match self.problem.sw.cert_screen {
                CertScreen::Off => None,
                mode => self
                    .problem
                    .certificate_shift_estimate(
                        &step.density,
                        &step.k_points,
                        self.problem.sw.potential_reuse,
                    )?
                    .map(|bound| (mode, bound)),
            }
        } else {
            None
        };
        let screen_seconds = screen_start.elapsed().as_secs_f64();
        if screen.is_some() {
            crate::phase_clock::add_screen(screen_seconds);
        }
        let screened = screen.and_then(|(mode, bound)| {
            (mode == CertScreen::On && bound > CERT_SCREEN_MARGIN * CERT_EIGENVALUE_LIMIT_HA)
                .then_some(bound)
        });
        let certified = if defer_verification {
            false
        } else if let Some(bound) = screened {
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
                // MPI-LAZY-ORBITALS: the spectral audit reads the full blocks.
                materialize_kpoints(&mut step.k_points)?;
                (rescued, unable_to_verify) = self.problem.audit_reported_spectrum(
                    n,
                    &mut step.k_points,
                    iteration,
                    report,
                )?;
                passed = report.passed;
            }
            if passed && self.problem.options.verification == Verification::FixedPoint {
                // ADR 0043: one more map at the delivered density; the delivered
                // state stays the converged step's own (V4).
                self.certificate_attempted = true;
                let (attempt, record, clipped, reuse) = verify_fixed_point(
                    &mut self.problem,
                    &step.density,
                    step.owned_density.as_ref(),
                    &step.k_points,
                    tolerance,
                    standard_error_ha.expect("standard history metric"),
                    iteration,
                    (n, step.energy.total_hartree()),
                )
                .map_err(|e| match e {
                    ScfError::NonFinite => ScfError::Diverged {
                        reason: "nonfinite verification arithmetic".into(),
                        history: history.clone(),
                    },
                    other => other,
                })?;
                if let Some(reason) = clipped {
                    return Err(ScfError::Diverged {
                        reason,
                        history: std::mem::take(history),
                    });
                }
                passed = attempt.passed;
                if !passed {
                    // feed2: the first failed attempt keeps the default
                    // continuation; repeated failures feed.
                    let earlier_failures = self
                        .verification
                        .attempts
                        .iter()
                        .filter(|a| !a.passed)
                        .count();
                    verification_reuse = match verify_reuse(&self.problem) {
                        VerifyReuse::Off => None,
                        VerifyReuse::Feed2 if earlier_failures == 0 => None,
                        mode => Some((mode, reuse)),
                    };
                }
                self.verification.attempts.push(attempt);
                history.last_mut().unwrap().fixed_point_certificate = Some(record);
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
        drop(verification);
        drop(profile_verification);
        if self.certificate_attempted {
            crate::phase_clock::add_certificate(certificate_clock);
        }
        if certified
            || unable_to_verify
            || (retain_last && iteration == self.problem.options.max_iterations)
        {
            // MPI-LAZY-ORBITALS: the delivered state carries full blocks (one
            // gather per lazy k block; every rank reaches this point together).
            // LAZY-RESULT: kept distributed when the caller declared it reads
            // no delivered orbitals (replicated decision) and the identity
            // audit below does not need them; `ScfResult::materialize_orbitals`
            // gathers them later. `orbitals` is then the 0 x B placeholder.
            if step.k_points.iter().any(|k| k.lazy.is_some()) {
                if step.nonlocal_from_identity || !crate::gamma_real::lazy_result_allowed() {
                    materialize_kpoints(&mut step.k_points)?;
                } else {
                    crate::m2::count(crate::m2::Count::LazyResults);
                }
                if let Some(first) = step.k_points.first() {
                    step.solution.eigenvectors = first.solution.eigenvectors.clone();
                }
            }
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
            // M2-OWNED: the delivered density is gathered once here (output).
            if let Some(o) = step.owned_density.take() {
                step.density = self.problem.materialize_owned(&o)?;
            }
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
                    verification: (standard.is_some()
                        && self.problem.options.verification != Verification::Off)
                        .then(|| self.verification.clone()),
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
        if let Some(limit) = self.problem.sw.growth_reset {
            let record = history.last_mut().unwrap();
            if record.density_contraction.is_some_and(|c| c > limit) && !record.mixer_reset {
                mixer.reset();
                *mixing_resets += 1;
                record.mixer_reset = true;
                record.mixing_resets = *mixing_resets;
            }
        }
        // VERIFY-REUSE: the failed verification evaluated the map at rho_out.
        // feed/feed2/half: its pair enters the mixer after the step's own and the
        // next input is mixed from it (half: midpoint with the step's proposal);
        // replace: only its pair enters (the step's pair is not recorded).
        let feed = verification_reuse
            .as_ref()
            .filter(|(mode, _)| *mode != VerifyReuse::Warm)
            .map(|(mode, (_, d, _))| (*mode, d));
        let mean = self.problem.electrons / self.problem.cell.volume_bohr3();
        let diverged = |e: dft_solvers::SolverError| ScfError::Diverged {
            reason: e.to_string(),
            history: history.clone(),
        };
        if let Some(om) = owned.as_mut() {
            // M2-OWNED: the Pulay update, the feed of a failed verification and
            // the convex safeguard on this rank's planes; the next input stays
            // owned (no gather). The same branch structure as below.
            if rescued.is_some() || history.last().is_some_and(|r| r.mixer_reset) {
                om.reset();
            }
            let input = self
                .problem
                .owned
                .input
                .take()
                .ok_or(ScfError::InvalidInput("M2 owned mixing without input"))?;
            let out = step
                .owned_density
                .as_ref()
                .ok_or(ScfError::InvalidInput("M2 owned mixing without output"))?;
            let fed = verification_reuse
                .as_ref()
                .filter(|(mode, _)| *mode != VerifyReuse::Warm)
                .map(|(mode, (_, _, d))| (*mode, d.as_ref()));
            let mut mixed = match fed {
                Some((VerifyReuse::Replace, _)) => Vec::new(),
                _ => om
                    .mix_block(&input.values, &out.values, mean)
                    .map_err(diverged)?,
            };
            if let Some((mode, fed)) = fed {
                let fed = fed.ok_or(ScfError::InvalidInput("M2 owned feed without density"))?;
                let step_proposal = std::mem::take(&mut mixed);
                mixed = om
                    .mix_block(&out.values, &fed.values, mean)
                    .map_err(diverged)?;
                if mode == VerifyReuse::Half {
                    for (m, s) in mixed.iter_mut().zip(&step_proposal) {
                        *m = 0.5 * (*m + s);
                    }
                }
            }
            let base: &[f64] = if fed.is_some() {
                &out.values
            } else {
                &input.values
            };
            history.last_mut().unwrap().mixing_history_used =
                Some(om.history_len().saturating_sub(1));
            history.last_mut().unwrap().mixing_raw_recovery = om.used_raw_recovery();
            history.last_mut().unwrap().pulay = om.pulay_step();
            let fraction = self.problem.owned_convex_fraction(base, &mixed)?;
            if fraction < 1.0 {
                for (b, a) in mixed.iter_mut().zip(base) {
                    *b = a + fraction * (*b - a);
                }
                om.reset();
                let record = history.last_mut().unwrap();
                record.mixing_fraction = fraction;
                record.mixer_reset = true;
                *mixing_resets += 1;
                record.mixing_resets = *mixing_resets;
            }
            let epoch = self.problem.owned.next_epoch();
            self.problem.owned.input = Some(crate::grid_dist::Owned {
                epoch,
                values: mixed,
            });
            history.last_mut().unwrap().timings.mixing = mix_start.elapsed();
            crate::profile::add_mixing(mix_start.elapsed().as_secs_f64());
            let verification_orbitals = verification_reuse.map(|(_, (o, _, _))| o);
            let (warm, lazy) = match (rescued, verification_orbitals) {
                (Some(r), _) => (r, Vec::new()),
                (None, Some(points)) => split_points(points),
                (None, None) => split_points(step.k_points),
            };
            *previous = Some(warm);
            *previous_lazy = lazy;
            return if iteration == self.problem.options.max_iterations {
                Ok(StepStatus::NotConverged)
            } else {
                Ok(StepStatus::Ready)
            };
        }
        let mut mixed = match feed {
            Some((VerifyReuse::Replace, _)) => Vec::new(),
            _ => mixer
                .mix(
                    n.as_slice().unwrap(),
                    step.density.as_slice().unwrap(),
                    mean,
                )
                .map_err(diverged)?,
        };
        if let Some((mode, fed)) = feed {
            let step_proposal = std::mem::take(&mut mixed);
            mixed = mixer
                .mix(
                    step.density.as_slice().unwrap(),
                    fed.as_slice().unwrap(),
                    mean,
                )
                .map_err(diverged)?;
            if mode == VerifyReuse::Half {
                for (m, s) in mixed.iter_mut().zip(&step_proposal) {
                    *m = 0.5 * (*m + s);
                }
            }
        }
        let base: &Array3<f64> = if feed.is_some() { &step.density } else { n };
        history.last_mut().unwrap().mixing_history_used =
            Some(mixer.history_len().saturating_sub(1));
        history.last_mut().unwrap().mixing_raw_recovery = mixer.used_raw_recovery();
        history.last_mut().unwrap().pulay = mixer.pulay_step();
        // Anderson extrapolation is unconstrained. A convex step toward its
        // proposal preserves charge and Fourier support while keeping XC in
        // its density domain. Never clip individual spatial samples.
        if matches!(
            self.problem.options.mixing,
            Mixing::Pulay { .. } | Mixing::KerkerPulay { .. }
        ) {
            let fraction = base
                .iter()
                .zip(&mixed)
                .filter(|(_, b)| **b < 0.0)
                .map(|(a, b)| 0.9 * a / (a - b))
                .fold(1.0, f64::min);
            if fraction < 1.0 {
                for (b, a) in mixed.iter_mut().zip(base.iter()) {
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
        crate::profile::add_mixing(mix_start.elapsed().as_secs_f64());
        // MPI-LAZY-ORBITALS: lazy blocks travel with their distributed rows.
        let verification_orbitals = verification_reuse.map(|(_, (o, _, _))| o);
        let (warm, lazy) = match (rescued, verification_orbitals) {
            (Some(r), _) => (r, Vec::new()),
            (None, Some(points)) => split_points(points),
            (None, None) => split_points(step.k_points),
        };
        *previous = Some(warm);
        *previous_lazy = lazy;
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
        // RUN-CONFIG: legacy switch reads below (and on the problem's worker
        // pool) see the problem's configuration.
        let _run = self.problem.run.enter();
        self.problem
            .run
            .count(dft_numerics::switches::RunCounter::SessionSteps, 1);
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
        let clock = Instant::now();
        let result = self.advance();
        crate::phase_clock::add_steps(clock);
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
    /// None at an active boundary of the owned density loop
    /// (`BAANDI_M2_LOCAL_DENSITY`): the input is held on the ranks' planes;
    /// `materialize_density` gathers it (collective).
    pub fn density(&self) -> Option<&Array3<f64>> {
        self.outcome
            .as_ref()
            .map(|o| &o.result.density)
            .or_else(|| {
                self.continuation
                    .as_ref()
                    .filter(|c| c.owned.is_none() || !c.n.is_empty())
                    .map(|c| &c.n)
            })
    }
    /// M2-OWNED: gather the owned next input density into the continuation so
    /// that `density()` returns it (one allgatherv, charged as a
    /// materialization; a no-op otherwise). Collective: every rank calls it at
    /// the same boundary. The owned loop continues from its owned input.
    pub fn materialize_density(&mut self) -> Result<(), ScfError> {
        let Some(c) = self.continuation.as_mut() else {
            return Ok(());
        };
        if c.owned.is_none() || !c.n.is_empty() {
            return Ok(());
        }
        let input = self
            .problem
            .owned
            .input
            .take()
            .ok_or(ScfError::InvalidInput("M2 owned boundary without input"))?;
        let full = self.problem.materialize_owned(&input);
        self.problem.owned.input = Some(input);
        c.n = full?;
        Ok(())
    }
    /// The next input density of an active boundary: kept on the planes in
    /// the owned loop (the full array is cleared so it cannot go stale).
    fn set_next_input(&mut self, n: Array3<f64>) -> Result<(), ScfError> {
        let owned = self
            .continuation
            .as_ref()
            .is_some_and(|c| c.owned.is_some());
        let n = if owned {
            let o = self.problem.owned_from_full(&n)?;
            self.problem.owned.input = Some(o);
            Array3::zeros((0, 0, 0))
        } else {
            n
        };
        if let Some(c) = self.continuation.as_mut() {
            c.n = n;
        }
        Ok(())
    }
    /// C1 ladder hand-off: next input density and warm orbitals at an active boundary.
    /// Lazy blocks (MPI-LAZY-ORBITALS) and an owned density (M2-OWNED) are
    /// gathered here (collective).
    pub(crate) fn warm_state(&mut self) -> Result<Option<WarmState<'_>>, ScfError> {
        let density_span = crate::cutoff_ladder::HandoffSpan::enter(
            crate::cutoff_ladder::HandoffPart::ExportDensity,
        );
        self.materialize_density()?;
        drop(density_span);
        let Some(c) = self.continuation.as_ref() else {
            return Ok(None);
        };
        let _span = crate::cutoff_ladder::HandoffSpan::enter(
            crate::cutoff_ladder::HandoffPart::Materialize,
        );
        Ok(Some((&c.n, full_previous(c)?)))
    }
    /// Export local planes and a lazy block. The descriptor agrees
    /// availability, consumer identity, overlaps and payload sizes in one collective.
    pub(crate) fn warm_owned_handoff(
        &mut self,
        fine: &ScfProblem,
    ) -> Result<Option<crate::cutoff_ladder::OwnedExport>, ScfError> {
        use dft_mpi::Collective;
        let Some(comm) = crate::parallel::SharedComm::active() else {
            return Ok(None);
        };
        let lazy = self
            .continuation
            .as_ref()
            .and_then(|c| c.previous_lazy.first())
            .and_then(Option::as_ref)
            .or_else(|| {
                self.result()
                    .and_then(|r| r.result.k_points.first())
                    .and_then(|k| k.lazy.as_ref())
            })
            .filter(|h| !h.handoff_materialized())
            .cloned();
        let dims = self.problem.grid.dims();
        let (density, axis) = if self.problem.owned.density_mode {
            let values = self
                .problem
                .owned
                .input
                .as_ref()
                .map(|n| n.values.clone())
                .unwrap_or_default();
            let axis = self.problem.grid_plan()?.1.first_axis();
            (values, axis)
        } else {
            let full = self
                .continuation
                .as_ref()
                .map(|c| &c.n)
                .or_else(|| self.result().map(|r| &r.result.density));
            let total = dims.iter().product::<usize>();
            let start = comm.rank() * (total / comm.size()) + comm.rank().min(total % comm.size());
            let len = total / comm.size() + usize::from(comm.rank() < total % comm.size());
            (
                full.and_then(|n| n.as_slice())
                    .and_then(|n| n.get(start..start + len))
                    .unwrap_or_default()
                    .to_vec(),
                2,
            )
        };
        // A deterministic fingerprint covers the fine consumer's validation
        // inputs. A local options/layout mismatch rejects before payload exchange.
        let identity = format!(
            "{:?}/{:?}/{:?}/{:?}",
            fine.options,
            fine.occupations,
            fine.electrons,
            fine.cell.lattice_columns_bohr()
        )
        .bytes()
        .fold(0xcbf29ce484222325u64, |h, b| {
            (h ^ u64::from(b)).wrapping_mul(0x100000001b3)
        });
        let eligible = fine.options.standard.is_some()
            && fine.owned_declined().is_none()
            && matches!(fine.options.mixing, Mixing::Pulay { .. })
            && !fine.options.hartree_mixing
            && crate::snapshot::dir().is_none()
            && crate::gamma_real::owned_handoff_supported();
        crate::cutoff_ladder::OwnedExport::prepare(
            lazy,
            density,
            dims,
            axis,
            &self.problem.states[0].basis,
            &fine.states[0].basis,
            fine.states[0].hamiltonian.fft_grid(),
            eligible,
            identity,
            self.problem.occupations.len(),
        )
    }
    pub(crate) fn problem_ref(&self) -> &ScfProblem {
        &self.problem
    }
    /// Mixer memory of the current continuation (None before initialization or
    /// for opaque mixers). Read-only; used by the cutoff ladder.
    pub(crate) fn mixer_state(&self) -> Option<crate::checkpoint::MixerState> {
        self.continuation
            .as_ref()
            .and_then(|c| c.mixer.checkpoint().ok())
    }
    /// Cutoff-ladder hand-off: seed the empty mixer of a just initialized session
    /// with difference pairs (dx, dr) of a related map, oldest first.
    pub(crate) fn seed_mixer_differences(
        &mut self,
        dx: &[Vec<f64>],
        dr: &[Vec<f64>],
    ) -> Result<(), ScfError> {
        if self.status != StepStatus::Ready || self.iteration != 0 {
            return Err(ScfError::InvalidInput(
                "seed_mixer_differences: only before the first step",
            ));
        }
        let mean = self.problem.electrons / self.problem.cell.volume_bohr3();
        if self
            .continuation
            .as_ref()
            .is_some_and(|c| c.owned.is_some())
        {
            // M2-OWNED: the seeds enter the plane-mode mixer as this rank's
            // planes (no collective).
            let (_, plan) = self.problem.grid_plan()?;
            let block = |v: &Vec<f64>| {
                let mut out = Vec::new();
                plan.planes_of(v, &mut out);
                out
            };
            let bx: Vec<Vec<f64>> = dx.iter().map(block).collect();
            let br: Vec<Vec<f64>> = dr.iter().map(block).collect();
            let c = self.continuation.as_mut().expect("checked");
            return c
                .owned
                .as_mut()
                .expect("checked")
                .seed_block(&bx, &br, mean)
                .map_err(|e| ScfError::Numerical(format!("seed_mixer_differences: {e}")));
        }
        let c = self
            .continuation
            .as_mut()
            .ok_or(ScfError::InvalidInput("seed_mixer_differences: no session"))?;
        c.mixer
            .seed_differences(dx, dr, mean)
            .map_err(|e| ScfError::Numerical(format!("seed_mixer_differences: {e}")))
    }
    /// Cutoff-ladder hand-off: first-step eigen tolerance of the standard
    /// profile (clamped to the ratchet's [floor, cap]). Ignored by strict.
    pub(crate) fn seed_standard_tolerance(&mut self, tolerance: f64) -> Result<(), ScfError> {
        if self.iteration != 0 || !tolerance.is_finite() || tolerance <= 0.0 {
            return Err(ScfError::InvalidInput(
                "seed_standard_tolerance: positive value before the first step",
            ));
        }
        self.standard_seed =
            Some(tolerance.clamp(StandardProfile::FLOOR_HA, StandardProfile::CAP_HA));
        Ok(())
    }
    /// Cutoff-ladder hand-off (D1): first-step eigen tolerance of the strict
    /// residual schedule, clamped by the schedule's own bounds. Ignored by the
    /// standard profile and by the legacy decay schedule.
    pub(crate) fn seed_strict_tolerance(&mut self, tolerance: f64) -> Result<(), ScfError> {
        if self.iteration != 0 || !tolerance.is_finite() || tolerance <= 0.0 {
            return Err(ScfError::InvalidInput(
                "seed_strict_tolerance: positive value before the first step",
            ));
        }
        self.strict_seed = Some(tolerance);
        Ok(())
    }
    /// Strict residual schedule: the eigen tolerance the next step would use
    /// (None for the standard profile, the legacy schedule or an empty history).
    pub(crate) fn next_strict_tolerance(&self) -> Option<f64> {
        let o = &self.problem.options;
        if o.standard.is_some() || !o.residual_eigen_schedule {
            return None;
        }
        self.history().last().map(|r| {
            residual_schedule(
                o,
                self.problem.sw.eigen_factor,
                Some(r.density_residual),
                1e-2,
            )
        })
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
            // MPI-LAZY-ORBITALS: lazy blocks are gathered (every rank runs the
            // same capture, so the collective is entered by all of them).
            match self.continuation.as_ref() {
                Some(c) if c.previous_lazy.iter().any(Option::is_some) => full_previous(c)
                    .map_err(|_| Rejection::Orbitals)?
                    .map(|v| v.iter().map(orbitals).collect())
                    .unwrap_or_default(),
                Some(c) => c
                    .previous
                    .as_ref()
                    .map(|v| v.iter().map(orbitals).collect())
                    .unwrap_or_default(),
                None => Vec::new(),
            }
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
            self.problem.forget_real_reuse();
            let c = self.continuation.as_ref().ok_or(Rejection::Boundary)?;
            // A pending first-step standard seed is not part of checkpoint v1.
            if (self.standard_seed.is_some() || self.strict_seed.is_some()) && c.history.is_empty()
            {
                return Err(Rejection::Boundary);
            }
            Some(ContinuationState {
                verification_attempts: self.verification.attempts.clone(),
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
                verification_attempts: r
                    .verification
                    .as_ref()
                    .map(|v| v.attempts.clone())
                    .unwrap_or_default(),
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
            verification: self.problem.options.verification,
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
    /// The ladder scopes initialization and import separately; ordinary
    /// checkpoint/direct-min consumers retain their existing entry point.
    pub(crate) fn initialize_ladder_electronic(
        &mut self,
        state: &crate::checkpoint::ElectronicState,
    ) -> Result<(), ScfError> {
        self.initialize()?;
        let _span =
            crate::cutoff_ladder::HandoffSpan::enter(crate::cutoff_ladder::HandoffPart::Import);
        self.import_electronic(state)
            .map_err(|e| ScfError::Numerical(format!("direct_min hand-off state rejected: {e}")))
    }
    pub(crate) fn initialize_owned_ladder(
        &mut self,
        state: &crate::checkpoint::ElectronicState,
        owned: crate::gamma_real::LazyOrbitals,
    ) -> Result<(), ScfError> {
        let _setup = crate::cutoff_ladder::OwnedSetup::enter();
        self.initialize_mode(true)?;
        let _import =
            crate::cutoff_ladder::HandoffSpan::enter(crate::cutoff_ladder::HandoffPart::Import);
        // The transcript already validated the global shape, overlaps and
        // identical consumer configuration. Density is replicated identically;
        // these unchanged import checks therefore have the same result everywhere.
        if !state.orbitals.is_empty()
            || owned.global_len() != self.problem.states[0].basis.vectors().len()
            || owned.bands() != self.problem.occupations.len()
        {
            return Err(ScfError::InvalidInput("owned handoff consumer shape"));
        }
        self.import_electronic(state)
            .map_err(|e| ScfError::Numerical(format!("owned handoff import: {e}")))?;
        let c = self.continuation.as_mut().expect("initialized");
        c.previous = Some(vec![
            ComplexMatrix::zeros(0, owned.bands())
                .map_err(|e| ScfError::Numerical(e.to_string()))?,
        ]);
        c.previous_lazy = vec![Some(owned)];
        Ok(())
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
                // LADDER-IMPORT-SPLIT (`BAANDI_IMPORT_CHECK=split`, default):
                // the same per-entry overlaps, lockstep over column blocks and
                // split over the ranks of a registered communicator (every
                // rank imports the same replicated state), one verdict
                // reduction; `serial` is the former replicated scalar loop.
                let split = self.problem.sw.import_split;
                let check_span = crate::cutoff_ladder::HandoffSpan::check();
                let orthonormal = if split {
                    let comm = crate::parallel::SharedComm::active();
                    orthonormal_split(&m, comm.as_ref().map(|c| c as &dyn dft_mpi::Collective))
                        .map_err(|_| R::Orbitals)?
                } else {
                    orthonormal_serial(&m)
                };
                drop(check_span);
                if !orthonormal {
                    return Err(R::Orbitals);
                }
                previous.push(m);
            }
            Some(previous)
        };
        if self.continuation.is_none() {
            return Err(R::Boundary);
        }
        // M2-OWNED: the imported density becomes the owned input of the loop.
        self.set_next_input(n).map_err(|_| R::Grid)?;
        let c = self.continuation.as_mut().ok_or(R::Boundary)?;
        c.previous = previous;
        c.previous_lazy = Vec::new();
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
        mut problem: ScfProblem,
        checkpoint: &crate::checkpoint::Checkpoint,
        model: &crate::checkpoint::ModelIdentity,
        execution: &crate::checkpoint::ExecutionIdentity,
    ) -> Result<Self, crate::checkpoint::Rejection> {
        use crate::checkpoint::*;
        // M2-OWNED: the restored mixer history and density are full arrays.
        problem.decline_owned(crate::grid_dist::OwnedDecline::Checkpoint);
        checkpoint.compatible(CheckpointKind::ExactContinuation, model, execution)?;
        let c = checkpoint.continuation.as_ref().ok_or(Rejection::Kind)?;
        if checkpoint.result.is_some() {
            return Err(Rejection::Kind);
        }
        if checkpoint.standard != c.options.standard
            || checkpoint.verification != c.options.verification
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
                            // A standard history holds a certificate record
                            // only for a verification attempt (ONEPATH-4).
                            || (r.fixed_point_certificate.is_some()
                                && !c
                                    .verification_attempts
                                    .iter()
                                    .any(|a| a.iteration == r.iteration))))
            })
            || !verification_attempts_consistent(c)
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
        session.verification.attempts = c.verification_attempts.clone();
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

/// ONEPATH-4: the saved verification attempts of an exact continuation are
/// the failed attempts of its history, one per history entry that carries a
/// verification record, in iteration order, each with finite numbers and the
/// session's threshold. Without the verification there are none.
fn verification_attempts_consistent(c: &crate::checkpoint::ContinuationState) -> bool {
    let attempts = &c.verification_attempts;
    if c.options.verification != Verification::FixedPoint || c.options.standard.is_none() {
        return attempts.is_empty();
    }
    let threshold = c.options.standard.map(|p| p.threshold_ha);
    let recorded: Vec<usize> = c
        .history
        .iter()
        .filter(|r| r.fixed_point_certificate.is_some())
        .map(|r| r.iteration)
        .collect();
    recorded.len() == attempts.len()
        && attempts.iter().zip(&recorded).all(|(a, &iteration)| {
            a.iteration == iteration
                && !a.passed
                && Some(a.threshold_ha.to_bits()) == threshold.map(f64::to_bits)
                && [
                    a.coulomb_error_ha,
                    a.eigen_tolerance_ha,
                    a.density_residual,
                    a.occupied_eigenvalue_change_ha,
                ]
                .iter()
                .chain(a.residual_bound_ha.iter())
                .all(|v| v.is_finite() && *v >= 0.0)
        })
}

/// Attempt, its history record, a clipped-charge failure, and the solve's
/// orbitals and output density F(rho_out) (VERIFY-REUSE).
/// Next input density and warm orbitals at an active boundary (ladder hand-off).
pub(crate) type WarmState<'a> = (&'a Array3<f64>, Option<Vec<ComplexMatrix>>);

/// The continuation's warm-start blocks with every lazy one gathered
/// (MPI-LAZY-ORBITALS; collective when any is lazy).
fn full_previous(c: &Continuation) -> Result<Option<Vec<ComplexMatrix>>, ScfError> {
    let Some(previous) = c.previous.as_ref() else {
        return Ok(None);
    };
    let mut out = Vec::with_capacity(previous.len());
    for (i, m) in previous.iter().enumerate() {
        out.push(match c.previous_lazy.get(i).and_then(Option::as_ref) {
            Some(lazy) if m.rows() == 0 => lazy.materialize()?.clone(),
            _ => m.clone(),
        });
    }
    Ok(Some(out))
}

/// Cost proxy from the measured attempt/map ratio: m / (m + 1.3 m).
const VERIFY_WAIT_BREAK_EVEN: f64 = 1.0 / 2.3;

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct VerificationPrediction {
    pub(crate) samples: usize,
    pub(crate) failure_vote: f64,
    pub(crate) extra_maps: usize,
}

/// A stateless retry forecast. Only existing scalar history is inspected; this
/// does not evaluate a density, alter a tolerance or certify a state. Keeping
/// all inputs in history also preserves exact checkpoint continuation.
pub(crate) fn verification_schedule(
    history: &[IterationRecord],
    attempts: &[VerificationAttempt],
) -> VerificationPrediction {
    let fallback = VerificationPrediction::default();
    let Some(last) = attempts.last().filter(|a| !a.passed) else {
        return fallback;
    };
    let Some(current) = history.last() else {
        return fallback;
    };
    let elapsed = current.iteration.saturating_sub(last.iteration);
    // Three maps is a bounded experiment, never an unbounded convergence veto.
    if elapsed >= 3
        || !last.coulomb_error_ha.is_finite()
        || last.coulomb_error_ha <= last.threshold_ha
    {
        return fallback;
    }
    let Some(end) = history.iter().position(|r| r.iteration == last.iteration) else {
        return fallback;
    };
    if history[end..]
        .iter()
        .any(|r| r.mixer_reset || r.mixing_raw_recovery)
    {
        return fallback;
    }
    let Some(anchor) = history[end].standard_error_ha else {
        return fallback;
    };
    let Some(now) = current.standard_error_ha else {
        return fallback;
    };
    if !anchor.is_finite() || anchor <= 0.0 || !now.is_finite() || now <= 0.0 {
        return fallback;
    }
    // Re-anchor after the ordinary maps actually performed. Pre-failure
    // contraction alone is poorly calibrated when a failure feeds the mixer.
    let estimate = last.coulomb_error_ha * (now / anchor);
    let start = history.len().saturating_sub(2);
    let Some(previous) = history.get(start).and_then(|r| r.standard_error_ha) else {
        return fallback;
    };
    let mut q = now / previous;
    if !estimate.is_finite() || !q.is_finite() || q <= 0.0 || q >= 1.0 {
        return fallback;
    }
    // A second failed attempt supplies an independent per-map rate. Combine
    // rates geometrically, without fitting a material-dependent scale.
    if attempts.len() >= 2 {
        let previous = &attempts[attempts.len() - 2];
        let distance = last.iteration.saturating_sub(previous.iteration);
        if distance == 0 {
            return fallback;
        }
        let measured =
            (last.coulomb_error_ha / previous.coulomb_error_ha).powf(1.0 / distance as f64);
        if !measured.is_finite() || measured <= 0.0 || measured >= 1.0 {
            return fallback;
        }
        q = (q * measured).sqrt();
    }
    // This point forecast is a 0/1 vote, not a probabilistic bound. Its
    // empirical calibration against real attempts is reported in the design.
    let vote = |extra: usize| f64::from(estimate * q.powi(extra as i32) >= last.threshold_ha);
    VerificationPrediction {
        samples: 1,
        failure_vote: vote(0),
        extra_maps: (0..=3 - elapsed)
            .find(|&extra| vote(extra) <= VERIFY_WAIT_BREAK_EVEN)
            .unwrap_or(3 - elapsed),
    }
}

/// Orbitals (k points, possibly lazy) and output density F(rho_out) of a
/// verification solve.
/// The third entry is the owned F(rho_out) of the owned density loop (the
/// array is then empty).
type VerificationMap = (
    Vec<KPointSolution>,
    Array3<f64>,
    Option<crate::grid_dist::Owned>,
);

/// Owned warm-start blocks and their lazy rows (MPI-LAZY-ORBITALS), moved out
/// of k points; the lazy list is empty when no block is lazy.
fn split_points(
    points: Vec<KPointSolution>,
) -> (
    Vec<ComplexMatrix>,
    Vec<Option<crate::gamma_real::LazyOrbitals>>,
) {
    let any_lazy = points.iter().any(|k| k.lazy.is_some());
    let mut warm = Vec::with_capacity(points.len());
    let mut lazy = Vec::new();
    for k in points {
        warm.push(k.solution.eigenvectors);
        if any_lazy {
            lazy.push(k.lazy);
        }
    }
    (warm, lazy)
}

type VerificationOutcome = (
    VerificationAttempt,
    FixedPointCertificate,
    Option<String>,
    VerificationMap,
);

/// ADR 0043 verification: one more ordinary map at the converged output density
/// (V1), pass on the Coulomb fixed-point residual at the standard threshold and,
/// with the residual guard on, the residual estimate of that solve (V2). Strict
/// certificate quantities are computed as information only (V3). Returns the
/// attempt, its record for the history entry, and a clipped-charge failure.
/// `stop_map` (input density, energy) feeds only the SCF-TAIL cross estimate.
#[allow(clippy::too_many_arguments)]
fn verify_fixed_point(
    problem: &mut ScfProblem,
    density: &Array3<f64>,
    owned: Option<&crate::grid_dist::Owned>,
    points: &[KPointSolution],
    tolerance: f64,
    standard_error_ha: f64,
    iteration: usize,
    stop_map: (&Array3<f64>, f64),
) -> Result<VerificationOutcome, ScfError> {
    // M2-OWNED: the verification map's input is the delivered (owned) output;
    // the loop's own input is restored afterwards on every path.
    let saved = owned.map(|o| problem.owned.input.replace(o.clone()));
    let result = verify_fixed_point_inner(
        problem,
        density,
        points,
        tolerance,
        standard_error_ha,
        iteration,
        stop_map,
    );
    if let Some(s) = saved {
        problem.owned.input = s;
    }
    result
}

fn verify_fixed_point_inner(
    problem: &mut ScfProblem,
    density: &Array3<f64>,
    points: &[KPointSolution],
    tolerance: f64,
    standard_error_ha: f64,
    iteration: usize,
    (stop_input, stop_energy_ha): (&Array3<f64>, f64),
) -> Result<VerificationOutcome, ScfError> {
    let clock = Instant::now();
    let threshold = problem
        .options
        .standard
        .expect("verification requires the standard profile")
        .threshold_ha;
    let eigen_tolerance = StandardProfile::ratchet(tolerance, standard_error_ha, problem.electrons);
    let (warm, warm_lazy) = split_warm(Some(points));
    crate::profile::set_context(crate::profile::Context {
        stage: "verification",
        kind: "verify",
        iteration,
    });
    problem
        .run
        .count(dft_numerics::switches::RunCounter::Verifications, 1);
    let check = problem.step_session(
        density,
        warm.as_deref(),
        &warm_lazy,
        eigen_tolerance,
        // Same flag as an ordinary step (the standard path ignores it).
        eigen_tolerance <= problem.options.eigen_tolerance,
    )?;
    let glue_clock = crate::profile::clock();
    let coulomb_error_ha = problem
        .stop_error_step(density, &check.density, check.owned_density.as_ref())?
        .expect("standard Coulomb metric");
    let residual_bound_ha = problem.sw.residual_guard.then(|| {
        StandardProfile::residual_energy_estimate(
            check
                .k_points
                .iter()
                .map(|k| (k.point.weight, k.solution.residuals.as_slice())),
            &problem.occupations,
        )
    });
    let (density_residual, _, _, _) =
        problem.step_metrics(density, &check.density, check.owned_density.as_ref())?;
    let occupied_eigenvalue_change_ha = check
        .k_points
        .iter()
        .zip(points)
        .flat_map(|(a, b)| {
            a.solution
                .eigenvalues
                .iter()
                .zip(&b.solution.eigenvalues)
                .zip(&problem.occupations)
        })
        .filter(|(_, f)| **f > 0.0)
        .map(|((a, b), _)| (a - b).abs())
        .fold(0.0, f64::max);
    if !coulomb_error_ha.is_finite()
        || !density_residual.is_finite()
        || !occupied_eigenvalue_change_ha.is_finite()
    {
        return Err(ScfError::NonFinite);
    }
    let guard_ok = residual_bound_ha.is_none_or(|b| b < threshold);
    #[allow(unused_mut)]
    let mut passed = coulomb_error_ha < threshold && guard_ok;
    // SCF-TAIL (`BAANDI_STOP_PREDICTOR`): the cross estimate needs the stop
    // map's input and both outputs as full arrays (the owned-plane loop keeps
    // them on the planes: old criterion only).
    if problem.sw.stop_predictor != StopPredictor::Off
        && !stop_input.is_empty()
        && !density.is_empty()
        && !check.density.is_empty()
        && let Some(metric) = problem.coulomb.as_mut()
    {
        let sum_ha = metric.error(stop_input, &check.density)?;
        let cross_ha = verification_cross(standard_error_ha, coulomb_error_ha, sum_ha);
        let cross_pass = cross_ha > 0.0
            && cross_ha < threshold / CROSS_SAFETY
            && standard_error_ha < threshold
            && guard_ok;
        eprintln!(
            "STOP-PREDICTOR iter={iteration} stop_ha={standard_error_ha:.6e} verify_ha={coulomb_error_ha:.6e} cross_ha={cross_ha:.6e} verify_energy_change_ha={:.6e} stop_energy_ha={stop_energy_ha:.15e} threshold_ha={threshold:.3e} old_pass={passed} cross_pass={cross_pass}",
            check.energy.total_hartree() - stop_energy_ha,
        );
        if problem.sw.stop_predictor == StopPredictor::On && !passed && cross_pass {
            passed = true;
        }
    }
    #[cfg(test)]
    if problem.force_verification_failures > 0 {
        problem.force_verification_failures -= 1;
        passed = false;
    }
    let mut clipped = None;
    if passed && problem.density_floor.is_some() {
        let charge = problem
            .negative_charge(density)
            .max(problem.negative_charge(&check.density));
        clipped = ScfProblem::validate_clipped_charge(charge, problem.electrons).err();
    }
    // Phase profile: file the verification map by its outcome; the snapshot
    // manifest of this map records it too (outside any timed phase).
    crate::profile::map_commit(if passed && clipped.is_none() {
        crate::profile::MapKind::VerifyPass
    } else {
        crate::profile::MapKind::VerifyFail
    });
    crate::snapshot::annotate_verification(passed && clipped.is_none());
    if passed && clipped.is_none() {
        problem.check_gap(&check.k_points)?;
    }
    crate::profile::between(crate::profile::BetweenPart::VerifyGlue, glue_clock);
    let residuals = |buffered: bool| -> Vec<Vec<f64>> {
        check
            .k_points
            .iter()
            .map(|k| {
                k.solution
                    .residuals
                    .iter()
                    .zip(&k.solution.buffer)
                    .filter_map(|(&r, &buffer)| (buffer == buffered).then_some(r))
                    .collect()
            })
            .collect()
    };
    let strict_residuals_per_k = residuals(false);
    let strict_criteria_met = density_residual <= problem.options.density_tolerance
        && occupied_eigenvalue_change_ha <= CERT_EIGENVALUE_LIMIT_HA
        && strict_residuals_per_k
            .iter()
            .flatten()
            .all(|r| *r <= 0.1 * problem.options.eigen_tolerance);
    let davidson_iterations = check.k_points.iter().map(|k| k.solution.iterations).sum();
    let operator_applications = check
        .k_points
        .iter()
        .map(|k| k.solution.operator_applications)
        .sum();
    let record = FixedPointCertificate {
        strict_residuals_per_k,
        buffer_residuals_per_k: residuals(true),
        density_residual,
        occupied_eigenvalue_change: occupied_eigenvalue_change_ha,
        davidson_iterations,
        operator_applications,
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
    };
    let attempt = VerificationAttempt {
        iteration,
        passed: passed && clipped.is_none(),
        coulomb_error_ha,
        threshold_ha: threshold,
        eigen_tolerance_ha: eigen_tolerance,
        residual_bound_ha,
        density_residual,
        occupied_eigenvalue_change_ha,
        strict_criteria_met,
        davidson_iterations,
        operator_applications,
        wall: clock.elapsed(),
        spin_density_residual: None,
        occupation_change_electrons: None,
        chemical_potential_change_ha: None,
        moment_change: None,
    };
    Ok((
        attempt,
        record,
        clipped,
        (check.k_points, check.density, check.owned_density),
    ))
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum VerifyReuse {
    Off,
    Warm,
    Feed,
    Feed2,
    Replace,
    Half,
}
/// `BAANDI_VERIFY_REUSE=off|warm|feed` (default off; ledger VERIFY-REUSE).
fn verify_reuse(#[allow(unused_variables)] problem: &ScfProblem) -> VerifyReuse {
    #[cfg(test)]
    if let Some(mode) = problem.force_verify_reuse {
        return mode;
    }
    problem.sw.verify_reuse
}

/// RUN-CONFIG (2026-11-18): the switches of the SCF session, resolved once
/// from the problem's run configuration at construction (and again by
/// `ScfProblem::set_run_config`), never looked up inside the iteration. The
/// parsing rules are those of the former per-process accessors.
#[derive(Clone, Copy, Debug)]
pub(crate) struct SessionSwitches {
    /// `BAANDI_STANDARD_STOP=confirm`.
    pub(crate) standard_stop_confirm: bool,
    /// `BAANDI_STANDARD_RESIDUAL_GUARD=on`.
    pub(crate) residual_guard: bool,
    /// `BAANDI_IMPORT_CHECK=split`.
    pub(crate) import_split: bool,
    /// `BAANDI_VERIFY_REUSE`.
    pub(crate) verify_reuse: VerifyReuse,
    /// `BAANDI_CERT_SCREEN`.
    pub(crate) cert_screen: CertScreen,
    /// `BAANDI_POTENTIAL_REUSE=on`.
    pub(crate) potential_reuse: bool,
    /// `BAANDI_MIX_GROWTH_RESET`.
    pub(crate) growth_reset: Option<f64>,
    /// `BAANDI_EIGEN_FACTOR` (default 0.1).
    pub(crate) eigen_factor: f64,
    /// `BAANDI_STOP_PREDICTOR` (SCF-TAIL, default off).
    pub(crate) stop_predictor: StopPredictor,
    /// F-VERIFY-SCHEDULE: resolved once; no environment reads in the loop.
    pub(crate) verify_schedule: &'static str,
}

/// SCF-TAIL (`BAANDI_STOP_PREDICTOR=off|trace|on`, default off): the cross
/// estimate of the fixed-point verification (`verification_cross`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum StopPredictor {
    Off,
    Trace,
    On,
}

/// Safety factor of the cross estimate (SCF-TAIL calibration, ledger
/// SCF-TAIL-PREDICTOR): an attempt passes on it only below threshold / factor;
/// the smallest measured -X / T (calibration systems) is 0.654, so 4 leaves a
/// margin of 2.6 over it.
pub(crate) const CROSS_SAFETY: f64 = 4.0;

/// Cross estimate of a verification (SCF-TAIL), returned as -X >= 0 in the
/// model. With r0 = rho_out - rho_in (the stop map's residual, E_H[r0] =
/// `stop_ha`) and r1 = F(rho_out) - rho_out (the verification's, E_H[r1] =
/// `verify_ha`), `sum_ha` = E_H[r0 + r1] = E_H[F(rho_out) - rho_in] gives the
/// Coulomb form X = r0 K r1 by polarization. Linear Hartree model
/// (rho_out = rho* + J delta, J = chi0 v = -v^-1/2 A v^1/2, A = v^1/2 (-chi0)
/// v^1/2 >= 0, x = v^1/2 delta): r0 = -(1 + A) x, r1 = A (1 + A) x; true
/// energy error of the delivered output T = x A (1 + A) x / 2, E_H[r0] =
/// x (1 + A)^2 x / 2, E_H[r1] = x A^2 (1 + A)^2 x / 2 and -X =
/// x A (1 + A)^2 x / 2, so -X / T = 1 + A >= 1 per mode while E_H[r1] / T =
/// A (1 + A) (below 1 for A < 0.618, eps (eps - 1) for screened long-wavelength
/// modes). Not a proof for the real map (XC kernel, inexact solves,
/// nonlinearity): measured -X / T reaches 0.654 (ledger SCF-TAIL-PREDICTOR), so
/// the pass rule carries the safety factor `CROSS_SAFETY`.
pub(crate) fn verification_cross(stop_ha: f64, verify_ha: f64, sum_ha: f64) -> f64 {
    -0.5 * (sum_ha - stop_ha - verify_ha)
}
impl SessionSwitches {
    pub(crate) fn resolve(run: dft_numerics::switches::RunConfig) -> Self {
        use dft_numerics::switches::Id;
        Self {
            standard_stop_confirm: run.choice(Id::STANDARD_STOP) == "confirm",
            residual_guard: run.on(Id::STANDARD_RESIDUAL_GUARD),
            import_split: run.choice(Id::IMPORT_CHECK) == "split",
            verify_reuse: match run.choice(Id::VERIFY_REUSE) {
                "warm" => VerifyReuse::Warm,
                "feed" => VerifyReuse::Feed,
                "feed2" => VerifyReuse::Feed2,
                "replace" => VerifyReuse::Replace,
                "half" => VerifyReuse::Half,
                _ => VerifyReuse::Off,
            },
            cert_screen: match run.choice(Id::CERT_SCREEN) {
                "off" => CertScreen::Off,
                "trace" => CertScreen::Trace,
                _ => CertScreen::On,
            },
            potential_reuse: run.on(Id::POTENTIAL_REUSE),
            growth_reset: run.number(Id::MIX_GROWTH_RESET),
            eigen_factor: run.number(Id::EIGEN_FACTOR).unwrap_or(0.1),
            verify_schedule: run.choice(Id::VERIFY_SCHEDULE),
            stop_predictor: match run.choice(Id::STOP_PREDICTOR) {
                "trace" => StopPredictor::Trace,
                "on" => StopPredictor::On,
                _ => StopPredictor::Off,
            },
        }
    }
}

/// Largest occupied-eigenvalue change (Ha) the fixed-point certificate accepts.
const CERT_EIGENVALUE_LIMIT_HA: f64 = 1e-7;
/// The pre-screen rejects only above this multiple of the limit, so its
/// second-order modelling error (about 1e-9 Ha) is far inside the margin.
const CERT_SCREEN_MARGIN: f64 = 1.02;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum CertScreen {
    Off,
    On,
    Trace,
}

// Certificate screen (`SessionSwitches::cert_screen`): default ON since the
// user's adopt-what-improves rule (2026-11-10; dev08 ledger SINGLE-NODE-AB1:
// bit-identical density and energy in 40 of 40 pairs, -19 percent wall on si32
// strict). `BAANDI_CERT_SCREEN=off` restores the unscreened path, `trace`
// prints the estimates. POTENTIAL-REUSE (`SessionSwitches::potential_reuse`):
// with the screen active, the Hartree/XC fields it builds for the output
// density are reused by the certificate solve at that same density
// (bit-identical either way). MIX_GROWTH_RESET: `SessionSwitches::growth_reset`.

/// Coefficient of the residual-driven eigen tolerance schedule (tolerance =
/// factor x density residual). Default 0.1 (the established empirical value);
/// `BAANDI_EIGEN_FACTOR=<x>` (0.001 to 1) is an experiment switch for dev08 A/B.
/// Strict residual eigen schedule: factor x last density residual (the
/// normalized density L2 norm in atomic units; 0.1 is an empirical scheduling
/// coefficient), `first` without history, clamped to [0.1 eigen_tolerance, 1e-2].
fn residual_schedule(
    options: &ScfOptions,
    factor: f64,
    last_residual: Option<f64>,
    first: f64,
) -> f64 {
    last_residual
        .map_or(first, |r| factor * r)
        .min(1e-2)
        .max(0.1 * options.eigen_tolerance)
}

/// Import check of a warm orbital block: every overlap <c_i|c_j> (i <= j) within
/// 1e-8 of the identity, each overlap the sequential sum over rows of
/// conj(c_i) c_j from zero (the former replicated loop; kept for
/// `BAANDI_IMPORT_CHECK=serial` and as the reference of the tests).
pub(crate) fn orthonormal_serial(m: &ComplexMatrix) -> bool {
    for j in 0..m.cols() {
        for i in 0..=j {
            let overlap: Complex64 = m
                .column(i)
                .iter()
                .zip(m.column(j))
                .map(|(a, b)| a.conj() * b)
                .sum();
            if (overlap - Complex64::new(if i == j { 1.0 } else { 0.0 }, 0.0)).norm() > 1e-8 {
                return false;
            }
        }
    }
    true
}

/// Lockstep width of `orthonormal_split` (independent accumulators per row pass).
const IMPORT_LOCKSTEP: usize = 8;

/// LADDER-IMPORT-SPLIT: the verdict of `orthonormal_serial` with the same
/// per-entry arithmetic (each overlap the sequential row sum from zero of
/// conj(c_i) c_j, so every overlap is bit-identical), computed IMPORT_LOCKSTEP
/// columns i at a time against one column j (independent accumulators instead
/// of one latency-bound chain) and, with a communicator, for the columns j of
/// this rank only (j mod P == rank; the triangle work balances cyclically).
/// One 1-double tree sum of the failure flags gives every rank the same
/// verdict. Collective when `comm` is given: every rank must call it with the
/// same block.
pub(crate) fn orthonormal_split(
    m: &ComplexMatrix,
    comm: Option<&dyn dft_mpi::Collective>,
) -> Result<bool, dft_mpi::MpiError> {
    let (rank, size) = comm.map_or((0, 1), |c| (c.rank(), c.size()));
    let mut bad = false;
    'columns: for j in (rank..m.cols()).step_by(size.max(1)) {
        let cj = m.column(j);
        let mut i0 = 0;
        while i0 <= j {
            let w = IMPORT_LOCKSTEP.min(j + 1 - i0);
            let mut acc = [Complex64::new(0.0, 0.0); IMPORT_LOCKSTEP];
            let cols: [&[Complex64]; IMPORT_LOCKSTEP] =
                std::array::from_fn(|t| m.column(if t < w { i0 + t } else { i0 }));
            for (k, &b) in cj.iter().enumerate() {
                for t in 0..w {
                    // AddAssign adds re and im separately, as `acc + x` does.
                    acc[t] += cols[t][k].conj() * b;
                }
            }
            for (t, overlap) in acc.iter().take(w).enumerate() {
                let i = i0 + t;
                if (overlap - Complex64::new(if i == j { 1.0 } else { 0.0 }, 0.0)).norm() > 1e-8 {
                    bad = true;
                    break 'columns;
                }
            }
            i0 += w;
        }
    }
    let Some(c) = comm else {
        return Ok(!bad);
    };
    let mut flag = [0.0];
    c.tree_sum_f64(&[f64::from(u8::from(bad))], &mut flag)?;
    Ok(flag[0] == 0.0)
}

#[cfg(test)]
mod import_check_tests {
    use super::*;

    fn block(rows: usize, cols: usize, seed: u64) -> ComplexMatrix {
        // Orthonormal columns by twice-repeated Gram-Schmidt of a pseudo-random block.
        let mut s = seed;
        let mut next = || {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((s >> 11) as f64 / (1u64 << 53) as f64) - 0.5
        };
        let mut m = ComplexMatrix::zeros(rows, cols).unwrap();
        for j in 0..cols {
            for v in m.column_mut(j).iter_mut() {
                *v = Complex64::new(next(), next());
            }
            for _ in 0..2 {
                for i in 0..j {
                    let p: Complex64 = m
                        .column(i)
                        .iter()
                        .zip(m.column(j))
                        .map(|(a, b)| a.conj() * b)
                        .sum();
                    let ci = m.column(i).to_vec();
                    for (v, a) in m.column_mut(j).iter_mut().zip(ci) {
                        *v -= p * a;
                    }
                }
            }
            let n = m.column(j).iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
            for v in m.column_mut(j).iter_mut() {
                *v /= n;
            }
        }
        m
    }

    /// Verdict identical to the serial loop: orthonormal blocks, defects just
    /// above and below 1e-8 at varied (i, j), NaN and infinity, for 1 to 5
    /// simulated ranks and widths around the lockstep width.
    #[test]
    fn split_import_check_matches_serial_verdict() {
        let mut cases = Vec::new();
        for (rows, cols) in [(1, 1), (7, 3), (40, 8), (61, 9), (90, 17), (33, 33)] {
            let m = block(rows, cols, rows as u64 * 31 + cols as u64);
            cases.push(m.clone());
            for (i, j, scale) in [(0, 0, 2e-8), (0, 0, 2e-9), (cols / 2, cols - 1, 3e-8)] {
                let mut d = m.clone();
                // Add a little of column i to column j: an overlap (or norm)
                // defect of about `scale`.
                let ci = d.column(i).to_vec();
                for (v, a) in d.column_mut(j).iter_mut().zip(ci) {
                    *v += scale * a;
                }
                cases.push(d);
            }
            let mut nan = m.clone();
            nan.column_mut(cols - 1)[rows - 1] = Complex64::new(f64::NAN, 0.0);
            cases.push(nan);
            let mut inf = m.clone();
            inf.column_mut(0)[0] = Complex64::new(f64::INFINITY, 0.0);
            cases.push(inf);
        }
        let mut verdicts = [0usize; 2];
        for m in &cases {
            let want = orthonormal_serial(m);
            verdicts[usize::from(want)] += 1;
            assert_eq!(orthonormal_split(m, None).unwrap(), want);
            for ranks in 2..=5 {
                let got = dft_mpi::simulate(ranks, |r| {
                    orthonormal_split(m, Some(r as &dyn dft_mpi::Collective)).unwrap()
                });
                assert!(got.iter().all(|&g| g == want), "ranks {ranks}");
            }
        }
        assert!(verdicts[0] > 0 && verdicts[1] > 0, "{verdicts:?}");
    }
}
