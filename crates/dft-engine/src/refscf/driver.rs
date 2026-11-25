//! Independent reference SCF with an optional coarse start. No ScfProblem/session or
//! legacy solver dispatch is constructed, including during verification.
// Reporting bridge for the CLI, whose production dependency is the engine.
pub use dft_solvers::ref_rmm::{
    STAGES as RMM_STAGES, scheduled as rmm_scheduled, take_records as take_rmm_records,
};

use super::memory::{self, Component as Mem, bytes, observe};
use super::stepper::{Boundary, StepKind, StepReport, Stepper};
use super::{
    broyden::{Broyden, Diagnostics},
    density::DensityPlan,
    map,
    metric::{Metric, Serial},
    numerical,
    profile::{Profile, next_ethr_ry},
    record::{CorrectionRecord, Record, SolveRecord, VerificationCost},
};
use crate::Functional;
use crate::{
    model::*,
    refstart::{self, StartInput, StartState},
};
use dft_geometry::{Cell, kpoints::KPoint};
use dft_mpi::Collective;
use dft_physics::data::{Complex64, ComplexMatrix};
use dft_solvers::EigenSolution;
use std::time::{Duration, Instant};

/// Empirical verification-error/SCF-estimator scaling, never a certificate.
/// A skipped candidate is followed by an actual observation at the next
/// eligible candidate, bounding how long stale calibration can defer a check.
#[derive(Default)]
pub struct VerificationPredictor {
    pub(super) observed: Option<(f64, f64)>,
    pub(super) skipped_last: bool,
}
impl VerificationPredictor {
    pub fn observe(&mut self, dr2_ry: f64, error_ha: f64) {
        self.observed = Some((dr2_ry, error_ha));
        self.skipped_last = false;
    }
    pub fn predict(&self, dr2_ry: f64) -> Option<f64> {
        let (previous, error) = self.observed?;
        if !previous.is_finite()
            || previous <= 0.
            || !error.is_finite()
            || error < 0.
            || !dr2_ry.is_finite()
            || dr2_ry <= 0.
            || dr2_ry >= previous
        {
            return None;
        }
        let prediction = error * (dr2_ry / previous);
        prediction.is_finite().then_some(prediction)
    }
    pub fn skip(&mut self, prediction: Option<f64>, threshold_ha: f64) -> bool {
        let skip =
            !self.skipped_last && prediction.is_some_and(|p| p.is_finite() && p >= threshold_ha);
        self.skipped_last = skip;
        skip
    }
}

pub struct Run {
    pub start_report: crate::refstart::StartReport,
    pub symmetry: Option<super::symmetry::Report>,
    pub verification_k_call_indices: Vec<usize>,
    pub memory_breakdown: memory::Breakdown,
    pub kept_hx: super::kept_hx::Report,
    pub ladder: Option<super::ladder::Report>,
    pub qe_clock: Vec<super::record::QeClock>,
    pub full_diagnostics: bool,
    pub outcome: ScfOutcome,
    pub records: Vec<Record>,
    pub reference_stopped: bool,
    pub first_reference_stop: Option<usize>,
    pub exit_reason: String,
    pub setup: super::record::SetupTimings,
    pub output: super::record::SetupTimings,
    pub start_s: f64,
    pub output_s: f64,
    pub verification_calls: Vec<SolveRecord>,
    /// Conservative named-array payload envelope, not allocator high water/RSS.
    pub peak_payload_bytes_per_rank: Vec<usize>,
}
pub fn run(input: StartInput<'_>, profile: Profile) -> Result<Run, ScfError> {
    stepper(input, profile).run_to_end()
}
/// Select the active communicator once and retain it through delivery.
pub fn stepper(input: StartInput<'_>, profile: Profile) -> Stepper<'_, Run> {
    let comm = crate::parallel::SharedComm::active();
    Stepper::new(move |boundary| async move {
        if let Some(comm) = comm {
            run_inner(input, profile, &comm, Some(comm.clone()), &boundary).await
        } else {
            run_inner(input, profile, &Serial, None, &boundary).await
        }
    })
}
pub fn stepper_with<'a, K: Collective>(
    input: StartInput<'a>,
    profile: Profile,
    comm: &'a K,
) -> Stepper<'a, Run> {
    Stepper::new(
        move |boundary| async move { run_inner(input, profile, comm, None, &boundary).await },
    )
}
/// Identical parsed input on every rank; also supports simulated collectives.
pub fn run_with<K: Collective>(
    input: StartInput<'_>,
    profile: Profile,
    comm: &K,
) -> Result<Run, ScfError> {
    stepper_with(input, profile, comm).run_to_end()
}

/// Continue an owned coarse seed through fresh fine-stage setup and verified
/// delivery. No mixer or kept-HX cache is transferred. This shares the exact
/// handoff used by `REFSCF_LADDER` and permits inspecting injection beforehand.
pub fn run_seeded_with<K: Collective>(
    input: StartInput<'_>,
    profile: Profile,
    comm: &K,
    seed: super::ladder::Seed,
) -> Result<Run, ScfError> {
    Stepper::new(move |boundary| async move {
        run_seeded_continuation(input, profile, comm, seed, &boundary).await
    })
    .run_to_end()
}
pub(super) async fn run_seeded_continuation<K: Collective>(
    input: StartInput<'_>,
    profile: Profile,
    comm: &K,
    seed: super::ladder::Seed,
    boundary: &Boundary,
) -> Result<Run, ScfError> {
    let StageRun::Fine(run) = run_stage(
        input,
        profile,
        comm,
        None,
        Some(seed),
        false,
        Instant::now(),
        boundary,
    )
    .await?
    else {
        unreachable!()
    };
    Ok(*run)
}

// Only run() can attach the exact owned communicator to a lazy result.
// Generic borrowed run_with() retains its eager public delivery contract.
pub(super) async fn run_inner<K: Collective>(
    input: StartInput<'_>,
    profile: Profile,
    comm: &K,
    owned_comm: Option<crate::parallel::SharedComm>,
    boundary: &Boundary,
) -> Result<Run, ScfError> {
    let start = Instant::now();
    use super::ladder::{FRACTION, MAX_COARSE_MAPS, THRESHOLD_RATIO};
    let ladder = dft_numerics::switches::on(dft_numerics::switches::Id::REFSCF_LADDER);
    let (seed, report) = if ladder {
        let mut coarse = input.options.clone();
        coarse.cutoff_hartree *= FRACTION;
        coarse.grid_override = None;
        coarse.verification = Verification::Off;
        coarse.max_iterations = coarse.max_iterations.min(MAX_COARSE_MAPS);
        let coarse_profile = Profile {
            threshold_ry: profile.threshold_ry
                * dft_numerics::switches::number(dft_numerics::switches::Id::REFSCF_LADDER_RATIO)
                    .unwrap_or(THRESHOLD_RATIO),
            ..profile
        }
        .validate()?;
        let StageRun::Coarse(seed, mut report) = run_stage(
            StartInput {
                options: &coarse,
                ..input
            },
            coarse_profile,
            comm,
            None,
            None,
            true,
            start,
            boundary,
        )
        .await?
        else {
            unreachable!()
        };
        report.coarse_seconds = start.elapsed().as_secs_f64();
        (Some(*seed), Some(report))
    } else {
        (None, None)
    };
    let StageRun::Fine(mut run) = run_stage(
        input, profile, comm, owned_comm, seed, false, start, boundary,
    )
    .await?
    else {
        unreachable!()
    };
    if let Some(coarse) = &report {
        run.kept_hx.merge(&coarse.kept_hx);
    }
    run.ladder = report;
    Ok(*run)
}

enum StageRun {
    Fine(Box<Run>),
    Coarse(Box<super::ladder::Seed>, super::ladder::Report),
}

#[allow(clippy::too_many_arguments)]
async fn run_stage<K: Collective>(
    input: StartInput<'_>,
    profile: Profile,
    comm: &K,
    owned_comm: Option<crate::parallel::SharedComm>,
    seed: Option<super::ladder::Seed>,
    coarse: bool,
    start: Instant,
    boundary: &Boundary,
) -> Result<StageRun, ScfError> {
    let stage_start = Instant::now();
    let memory_scope = memory::Scope::begin();
    // Resolve once from this run's scope, never a process-cached profile flag.
    use dft_numerics::switches::{Id, RunConfig};
    let config = RunConfig::current();
    let full_diagnostics = config.on(Id::DAVIDSON_STAGES);
    let trim_memory = config.on(Id::REFSCF_MEMORY_TRIM);
    let rmm = config.on(Id::REFSCF_RMM).then(|| {
        let mut options = dft_solvers::ref_rmm::Options::new(0, 0., 0);
        options.steps = config.number(Id::REFSCF_RMM_STEPS).unwrap_or(4.) as usize;
        options.block = config.number(Id::REFSCF_RMM_BLOCK).unwrap_or(16.) as usize;
        options
    });
    let qe_scope = full_diagnostics.then(dft_solvers::ref_davidson::QeClockScope::begin);
    let mut qe_verification = dft_solvers::ref_davidson::QeCounters::default();
    let mut qe_mixing_calls = 0;

    let mut profile = profile.validate()?;
    let o = input.options;
    if o.standard.is_none()
        || o.gamma_real != crate::model::GammaReal::On
        || o.crystal_density_group
        || o.max_iterations == 0
    {
        return Err(ScfError::InvalidInput(
            "reference SCF requires real Gamma, standard class, no symmetry and positive iteration budget",
        ));
    }
    let mut options = o.clone();
    options.standard = Some(StandardProfile {
        threshold_ha: profile.threshold_ry / 2.,
    });
    let cell = input.cell;
    let prolonged_start = seed.is_some();
    let mut state = refstart::start_prolonged(
        StartInput {
            options: &options,
            ..input
        },
        comm,
        seed,
    )?;
    let mut setup = std::mem::take(&mut state.setup);
    let prep_comm = refstart::TimedCollective::new(comm);
    let clock = Instant::now();
    let dims = state.grid.dims();
    let volume = cell.volume_bohr3();
    let dv = volume / dims.iter().product::<usize>() as f64;
    let bands = state.orbitals.cols();
    let occupied = state.electrons as usize / 2;
    let occupations: Vec<_> = (0..bands)
        .map(|b| if b < occupied { 2. } else { 0. })
        .collect();
    let mut plan = DensityPlan::with_fft(
        &prep_comm,
        cell,
        dims,
        4. * o.cutoff_hartree,
        4. * o.cutoff_hartree,
        state.fields.shared_fft(),
    )?;
    if plan.first_axis() != state.fields.first_axis() || plan.plane_len() != state.density.len() {
        return Err(ScfError::InvalidInput("reference field layouts differ"));
    }
    setup.record("density_fft_plan", clock);
    let clock = Instant::now();
    let mut metric = Metric::new(volume, &plan.g2, &plan.smooth, comm.size(), 350)?;
    let mut mixer = Broyden::new(profile, &plan.g2, &plan.smooth, volume, state.electrons)?;
    let mut n = std::mem::take(&mut state.density);
    let mut x = vec![Complex64::default(); plan.miller.len()];
    let mut y = x.clone();
    let mut f = x.clone();
    let mut next = x.clone();
    let mut next_planes = vec![0.; n.len()];
    plan.forward(&prep_comm, &n, &mut x)?;
    set_zero(&plan, &mut x, state.electrons / volume);
    observe(
        Mem::Start,
        bytes(&state.density)
            + bytes(&state.density_g)
            + bytes(&state.potential)
            + bytes(&state.ionic)
            + bytes(&state.diagonal_static)
            + bytes(&state.ritz_values_hartree),
    );
    if trim_memory {
        // These start-only outputs have no consumer in the reference map. Drop
        // backing storage (clear would retain capacity) before the first solve.
        observe(
            Mem::ReleasedStart,
            bytes(&state.density_g) + bytes(&state.potential) + bytes(&state.ritz_values_hartree),
        );
        state.density_g = Vec::new();
        state.potential = Vec::new();
        state.ritz_values_hartree = Vec::new();
    }
    observe(
        Mem::DriverReciprocal,
        bytes(&x) + bytes(&y) + bytes(&f) + bytes(&next),
    );
    observe(Mem::DriverPlanes, bytes(&n) + bytes(&next_planes));
    // Peak envelope: startup, owned Davidson arrays/factors plus current field,
    // mixer/history buffers. Excludes backend/FFT metadata, stacks and allocator.
    let cap = (2 * bands).min(state.basis.vectors().len());
    let peak = state
        .memory
        .principal_peak_bound_bytes
        .max(
            state.memory.retained_field_bytes
                + state.memory.real_projector_bytes
                + 8 * state.operator.local_len() * (5 * cap + 3 * bands)
                + 8 * (20 * cap * cap + 4 * cap)
                + 8 * state.operator.local_len().div_ceil(256) * cap * cap
                + 16 * x.len() * (2 * profile.history + 12)
                + 8 * n.len() * 16,
        )
        .max(
            // Delivery currently materializes both public orbital copies and the
            // reduction's full real input/output buffers on every rank.
            state.memory.retained_field_bytes
                + state.memory.real_projector_bytes
                + 64 * state.basis.vectors().len() * bands
                + 16 * dims.iter().product::<usize>()
                + 16 * x.len() * (2 * profile.history + 12)
                + 8 * state.operator.local_len() * (5 * cap + 3 * bands),
        )
        // Exact field cache (density + HXC + VH), retained candidate/check
        // input fields and a temporary cloned field pair.
        + 8 * n.len() * 9
        + {
            // Retained X/HX, supplied HX and delta workspace, owned potentials,
            // bounded slabs and the existing full-index delta-kernel adapter.
            32 * state.operator.local_len() * bands
                + 24 * state.operator.plane_points()
                + (16 * (bands / 2) * state.operator.plane_points())
                    .min(super::kept_hx::SLAB_LIMIT_BYTES)
                + 8 * dims.iter().product::<usize>()
        };
    let mut peaks = vec![0.; comm.size()];
    let mut local_peak = peaks.clone();
    local_peak[comm.rank()] = peak as f64;
    prep_comm
        .tree_sum_f64(&local_peak, &mut peaks)
        .map_err(numerical)?;
    let peak_payload_bytes_per_rank = peaks.into_iter().map(|v| v as usize).collect();
    setup.record("mixer_workspace_initial_transform", clock);
    setup.collective_s += prep_comm.seconds();
    let start_time = start.elapsed();
    setup.total_s = stage_start.elapsed().as_secs_f64();
    let mut history = Vec::<IterationRecord>::with_capacity(o.max_iterations);
    let mut records = Vec::with_capacity(o.max_iterations);
    let mut ethr = 1e-2;
    let mut last_dr2 = f64::INFINITY;
    let mut verification =
        (o.verification == Verification::FixedPoint).then(VerificationReport::default);
    let mut verification_calls = Vec::new();
    let mut first_reference_stop = None;
    let mut pending: Option<(map::Map, f64, usize)> = None;
    let mut predictor = VerificationPredictor::default();
    let permap = config.on(Id::REFSCF_PERMAP);
    let mut field_cache = None;

    let mut kept = super::kept_hx::KeptHx::new(true);
    let mut completed = 0usize;
    // R8b: closure creation borrows locals only. State visits occur exclusively
    // on an explicit checkpoint poll; ordinary maps allocate nothing here.
    macro_rules! checkpoint_state {
        () => {
            |request: super::stepper::Visit<'_, '_>| -> Result<usize, ScfError> {
                let io = match request {
                    super::stepper::Visit::Checkpoint(io) => io,
                    super::stepper::Visit::Policy(patch) => {
                        if coarse || config.on(Id::REFSCF_LADDER) {
                            return Err(ScfError::InvalidInput(
                                "reference runtime policy changes require a fixed cutoff",
                            ));
                        }
                        let next = patch.prepare(profile)?;
                        if patch.beta.is_some() || patch.history.is_some() {
                            let replacement = Broyden::new(
                                next,
                                &plan.g2,
                                &plan.smooth,
                                volume,
                                state.electrons,
                            )?;

                            mixer = replacement;
                        }
                        if let Some(enabled) = patch.verification {
                            if enabled != verification.is_some() {
                                verification = enabled.then(VerificationReport::default);
                            }
                        }
                        options.standard = Some(StandardProfile {
                            threshold_ha: next.threshold_ry / 2.,
                        });
                        profile = next;
                        pending = None;
                        predictor = VerificationPredictor::default();
                        return Ok(completed);
                    }
                };
                if coarse || prolonged_start || config.on(Id::REFSCF_LADDER) {
                    return Err(super::checkpoint::invalid(
                        "checkpoint requires an unseeded fixed-cutoff run",
                    ));
                }
                super::checkpoint::header(io, &input, profile, comm, "gamma", dims)?;
                io.same(&plan.g2)?;
                io.value(&mut completed)?;
                io.planes(&mut n)?;
                io.reciprocal(&mut x)?;
                mixer.checkpoint(io, &metric)?;
                super::checkpoint::gamma_operator(io, &mut state)?;
                kept.checkpoint(io, &state.operator, bands)?;
                io.value(&mut ethr)?;
                io.value(&mut last_dr2)?;
                io.value(&mut history)?;
                io.value(&mut records)?;
                io.value(&mut verification)?;
                io.value(&mut verification_calls)?;
                io.value(&mut first_reference_stop)?;
                io.value(&mut predictor.observed)?;
                io.value(&mut predictor.skipped_last)?;
                io.value(&mut field_cache)?;
                io.value(&mut pending)?;
                super::checkpoint::cached_fields(&field_cache, n.len())?;
                if let Some((map, _, _)) = &pending {
                    map.validate_checkpoint(state.operator.local_len(), bands, n.len())?;
                    return Err(super::checkpoint::invalid(
                        "retired failed-verification reuse",
                    ));
                }
                super::checkpoint::counters(
                    completed,
                    o.max_iterations,
                    history.len(),
                    ethr,
                    last_dr2,
                )?;
                if records.len() != completed
                    || verification.as_ref().is_some_and(|v| {
                        v.attempts
                            .iter()
                            .any(|a| a.passed || a.iteration > completed)
                    })
                {
                    return Err(super::checkpoint::invalid("record/verification history"));
                }
                Ok(completed)
            }
        };
    }
    boundary.ready_state(checkpoint_state!()).await;
    for iteration in completed + 1..=o.max_iterations {
        let map_start_s = start.elapsed().as_secs_f64();
        let verification_call_start = verification_calls.len();
        let reused_from = None;
        {
            if iteration > 1 {
                ethr = next_ethr_ry(
                    if iteration == 2 { 1e-2 } else { ethr },
                    last_dr2,
                    state.electrons,
                );
            }
        }
        let requested = ethr;
        let mut redo = false;
        let mut reported_columns = 0;
        let mut calls = Vec::new();
        let mut timings = PhaseTimings::default();
        let (mut map_s, mut transform_s, mut metric_s) = (0., 0., 0.);
        let (mut step, dr2) = loop {
            let clock = Instant::now();
            let step = {
                let step = solve_map(
                    &mut state,
                    comm,
                    cell,
                    o.functional,
                    &n,
                    &occupations,
                    ethr,
                    iteration == 1 && !redo,
                    full_diagnostics,
                    rmm.filter(|_| dft_solvers::ref_rmm::scheduled(iteration))
                        .map(|mut o| {
                            o.scf_iteration = iteration;
                            o
                        }),
                    &mut calls,
                    &mut timings,
                    &mut field_cache,
                    true,
                    &mut kept,
                )?;
                map_s += clock.elapsed().as_secs_f64();
                step
            };
            let clock = Instant::now();
            plan.forward(comm, &step.density, &mut y)?;
            set_zero(&plan, &mut y, state.electrons / volume);
            for ((v, y), x) in f.iter_mut().zip(&y).zip(&x) {
                *v = y - x;
            }
            transform_s += clock.elapsed().as_secs_f64();
            let clock = Instant::now();
            let dr2 = metric.norm(comm, &f)?;
            metric_s += clock.elapsed().as_secs_f64();
            if iteration == 1 && !redo && dr2 < ethr * state.electrons.max(1.) {
                let mut report = StepReport::new(
                    iteration,
                    StepKind::Ordinary,
                    dr2,
                    step.energy.total_hartree(),
                    ethr,
                    calls.iter().map(|c| c.h_columns).sum(),
                );
                report.stage = if coarse { "coarse" } else { "fine" };
                reported_columns = report.h_columns;
                boundary.emit(report).await;
                ethr = next_ethr_ry(ethr, dr2, state.electrons);
                redo = true;
                kept.invalidate();
                state.orbitals = step.solve.eigenvectors;
                continue;
            }
            break (step, dr2);
        };
        last_dr2 = dr2;
        let stopped = dr2 < profile.threshold_ry;
        let estimator_s = start.elapsed().as_secs_f64();
        if stopped && first_reference_stop.is_none() {
            first_reference_stop = Some(iteration);
        }
        let mut certificate = None;
        // Compute before the attempt and before inserting this iterate into
        // the mixer. Later attempts retain the adopted empirical gate.
        let first_prediction_ha = None;
        let contraction_prediction_ha = None;
        let predicted_error_ha = if stopped && verification.is_some() {
            predictor.predict(dr2)
        } else {
            None
        };
        let skipped = stopped
            && verification.is_some()
            && predictor.skip(predicted_error_ha, profile.threshold_ry / 2.);
        if ((stopped && !skipped) || iteration == o.max_iterations) && step.energy.kinetic.is_nan()
        {
            let clock = Instant::now();
            map::complete_audit(
                &mut step,
                &mut state,
                comm,
                cell,
                o.functional,
                &n,
                &occupations,
                &mut field_cache,
            )?;
            let elapsed = clock.elapsed();
            timings.energy += elapsed;
            map_s += elapsed.as_secs_f64();
        }
        let mut map_report = StepReport::new(
            iteration,
            if redo {
                StepKind::Retry
            } else {
                StepKind::Ordinary
            },
            dr2,
            step.energy.total_hartree(),
            ethr,
            calls.iter().map(|c| c.h_columns).sum::<usize>() - reported_columns,
        );
        map_report.stage = if coarse { "coarse" } else { "fine" };
        map_report.reused = reused_from.is_some();
        let mut verification_warm = None;
        let mut actual_error_ha = None;
        let mut verification_s = 0.;
        let (mut verification_start_s, mut verification_end_s) = (None, None);
        if stopped
            && !skipped
            && let Some(report) = &mut verification
        {
            boundary.emit(map_report).await;
            let qe_before = qe_scope.as_ref().map(|s| s.snapshot());
            verification_start_s = Some(start.elapsed().as_secs_f64());
            let clock = Instant::now();
            state.orbitals = step.solve.eigenvectors.clone();
            let tol = next_ethr_ry(ethr, dr2, state.electrons);
            let mut t = PhaseTimings::default();
            let first_call = verification_calls.len();

            kept.preserve_input(true);
            let check = solve_map(
                &mut state,
                comm,
                cell,
                o.functional,
                &step.density,
                &occupations,
                tol,
                false,
                full_diagnostics,
                None, // Verification always uses a converged Davidson solve.
                &mut verification_calls,
                &mut t,
                &mut field_cache,
                false,
                &mut kept,
            )?;
            kept.preserve_input(false);
            let current_calls = &verification_calls[first_call..];
            verification_warm = Some(check.warm_start);
            let mut checked = y.clone();
            observe(Mem::Verification, bytes(&checked));
            plan.forward(comm, &check.density, &mut checked)?;
            for (v, y) in checked.iter_mut().zip(&y) {
                *v -= y;
            }
            let error = metric.norm(comm, &checked)? / 2.;
            actual_error_ha = Some(error);
            {
                predictor.observe(dr2, error);
            }
            let mut l2 = [0.];
            comm.tree_sum_f64(
                &[check
                    .density
                    .iter()
                    .zip(&step.density)
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()],
                &mut l2,
            )
            .map_err(numerical)?;
            let record = FixedPointCertificate {
                strict_residuals_per_k: vec![check.residuals[..occupied].to_vec()],
                buffer_residuals_per_k: vec![check.residuals[occupied..].to_vec()],
                density_residual: (dv * l2[0]).sqrt(),
                occupied_eigenvalue_change: check
                    .solve
                    .eigenvalues
                    .iter()
                    .zip(&step.solve.eigenvalues)
                    .take(occupied)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0., f64::max),
                davidson_iterations: current_calls.iter().map(|c| c.iterations).sum(),
                operator_applications: current_calls.iter().map(|c| c.h_columns).sum(),
                davidson_iterations_per_k: vec![current_calls.iter().map(|c| c.iterations).sum()],
                operator_applications_per_k: vec![current_calls.iter().map(|c| c.h_columns).sum()],
                timings: t,
            };
            let bound = StandardProfile::residual_energy_estimate(
                std::iter::once((1., check.residuals.as_slice())),
                &occupations,
            );
            let attempt = crate::model::assess_fixed_point(
                &options,
                iteration,
                tol / 2.,
                error,
                Some(bound),
                &record,
                clock.elapsed(),
            );

            map_report = StepReport::new(
                iteration,
                StepKind::Verification,
                error * 2.,
                check.energy.total_hartree(),
                tol,
                current_calls.iter().map(|c| c.h_columns).sum(),
            );
            map_report.stage = if coarse { "coarse" } else { "fine" };

            report.attempts.push(attempt);
            certificate = Some(record);
            verification_s = clock.elapsed().as_secs_f64();
            verification_end_s = Some(start.elapsed().as_secs_f64());
            if let (Some(scope), Some(before)) = (&qe_scope, qe_before) {
                qe_verification.add(scope.snapshot().without(before));
            }
        }
        let finish = stopped && verification.as_ref().is_none_or(VerificationReport::passed);
        let mut diagnostics = Diagnostics {
            iteration,
            occupancy: mixer.occupancy(),
            ..Default::default()
        };
        let (mut mixing_s, mut correction_s, mut descf) = (0., 0., 0.);
        let mut variational = step.energy.band_form;
        let (mut vh, mut vx) = (step.energy.hartree, step.energy.xc);
        let mut failure = None;
        let mut fused_local = None;
        if !finish {
            if full_diagnostics {
                qe_mixing_calls += 1;
            }
            let clock = Instant::now();
            match mixer.update(comm, &mut metric, &x, &f, &mut next) {
                Ok(d) => diagnostics = d,
                Err(e) => failure = Some(e.to_string()),
            }
            if failure.is_none() {
                {
                    plan.inverse(comm, &next, &mut next_planes)?;
                }
            }
            mixing_s = clock.elapsed().as_secs_f64();
            if failure.is_none() {
                let clock = Instant::now();
                let fields = map::cached_fields(
                    &mut state,
                    comm,
                    cell,
                    o.functional,
                    &next_planes,
                    &mut field_cache,
                )?;
                if permap {
                    // Reuse the exact LEAN_ENERGY cache entry. Each accumulator
                    // retains its original point order; one grid pass and one
                    // four-scalar tree sum replace four passes and two sums.
                    let mut local = [0.; 4];
                    for (((n, y), next), v) in n
                        .iter()
                        .zip(&step.density)
                        .zip(&next_planes)
                        .zip(&fields.hxc)
                    {
                        local[0] += (n - y).powi(2);
                        local[1] += n;
                        local[2] += y;
                        local[3] += (next - y) * v;
                    }
                    local[3] *= -dv;
                    fused_local = Some(local);
                } else {
                    let local = -dv
                        * next_planes
                            .iter()
                            .zip(&step.density)
                            .zip(&fields.hxc)
                            .map(|((n, y), v)| (n - y) * v)
                            .sum::<f64>();
                    let mut total = [0.];
                    comm.tree_sum_f64(&[local], &mut total).map_err(numerical)?;
                    descf = total[0];
                }
                vh = fields.hartree;
                vx = fields.xc;
                variational = step.energy.band_sum - step.energy.double_counting
                    + vh
                    + vx
                    + step.energy.ewald
                    + descf;
                correction_s = clock.elapsed().as_secs_f64();
            }
        }
        if failure.is_some() && step.energy.kinetic.is_nan() {
            let clock = Instant::now();
            map::complete_audit(
                &mut step,
                &mut state,
                comm,
                cell,
                o.functional,
                &n,
                &occupations,
                &mut field_cache,
            )?;
            let elapsed = clock.elapsed();
            timings.energy += elapsed;
            map_s += elapsed.as_secs_f64();
            variational = step.energy.band_form;
            vh = step.energy.hartree;
            vx = step.energy.xc;
            if map_report.kind != StepKind::Verification {
                let energy = step.energy.total_hartree();
                map_report.energy_ha = energy.is_finite().then_some(energy);
            }
        }
        let mut sums = [0.; 3];
        if let Some(local) = fused_local {
            let clock = Instant::now();
            let mut total = [0.; 4];
            comm.tree_sum_f64(&local, &mut total).map_err(numerical)?;
            sums.copy_from_slice(&total[..3]);
            descf = total[3];
            variational = step.energy.band_sum - step.energy.double_counting
                + vh
                + vx
                + step.energy.ewald
                + descf;
            correction_s += clock.elapsed().as_secs_f64();
        } else {
            comm.tree_sum_f64(
                &[
                    n.iter()
                        .zip(&step.density)
                        .map(|(n, y)| (n - y).powi(2))
                        .sum::<f64>(),
                    n.iter().sum::<f64>(),
                    step.density.iter().sum::<f64>(),
                ],
                &mut sums,
            )
            .map_err(numerical)?;
        }
        let dav = calls.iter().map(|c| c.iterations).sum();
        let columns = calls.iter().map(|c| c.h_columns).sum();
        timings.mixing = Duration::from_secs_f64(mixing_s);
        history.push(IterationRecord {
            iteration,
            energy: step.energy.clone(),
            energy_change: history
                .last()
                .map(|r| step.energy.total_hartree() - r.energy.total_hartree())
                .filter(|e| e.is_finite()),
            standard_error_ha: Some(dr2 / 2.),
            standard_redo: redo,
            density_residual: (dv * sums[0]).sqrt(),
            max_eigen_residual: step.residuals.iter().copied().fold(0., f64::max),
            requested_eigen_tolerance: ethr / 2.,
            davidson_iterations: dav,
            operator_applications: columns,
            davidson_iterations_per_k: vec![dav],
            operator_applications_per_k: vec![columns],
            charge_error: dv * sums[2] - state.electrons,
            input_charge_error: dv * sums[1] - state.electrons,
            min_density: global_min(comm, &step.density)?,
            mixing_fraction: 1.,
            mixer_reset: false,
            density_contraction: None,
            mixing_history_used: Some(diagnostics.occupancy),
            mixing_resets: 0,
            mixing_raw_recovery: false,
            fixed_point_certificate: certificate,
            timings: timings.clone(),
            pulay: None,
        });
        records.push(Record {
            k_call_indices: Vec::new(),
            warm_start: step.warm_start && !calls.is_empty(),
            stage: if coarse { "coarse" } else { "fine" },
            iteration,
            verification_cost: Some(VerificationCost {
                gate_enabled: true,
                first_gate_enabled: false,
                first_prediction_ha,
                contraction_prediction_ha,
                skipped,
                predicted_error_ha,
                actual_error_ha,
                reused_from,
            }),
            map_start_s,
            estimator_s,
            verification_energy: None,
            verification_warm,
            verification_start_s,
            verification_end_s,
            map_end_s: start.elapsed().as_secs_f64(),
            dr2_ry: dr2,
            ethr_requested_ry: requested,
            ethr_ry: ethr,
            first_iteration_redo: redo,
            beta: profile.beta,
            mixing: diagnostics,
            calls,
            verification_calls: verification_call_start..verification_calls.len(),
            davidson_iterations: dav,
            h_columns: columns,
            map_s,
            solve_s: timings.diagonalization.as_secs_f64(),
            density_s: timings.density.as_secs_f64(),
            potential_s: timings.potential.as_secs_f64(),
            energy_s: timings.energy.as_secs_f64(),
            density_transform_s: transform_s,
            metric_s,
            mixing_s,
            correction_s,
            verification_s,
            descf_ha: descf,
            variational_energy_ha: variational,
            energy: step.energy.clone(),
            variational_hartree_ha: vh,
            variational_xc_ha: vx,
        });
        map_report.safe_point = failure.is_none();
        if finish || failure.is_some() || iteration == o.max_iterations {
            if !coarse {
                boundary.publish(map_report);
            }
            let verified = verification.as_ref().is_none_or(VerificationReport::passed);
            let converged = stopped && verified && failure.is_none();
            let failed = failure.is_some();
            let exit_reason = failure.unwrap_or_else(|| {
                if !stopped {
                    "reference SCF maximum iterations reached".into()
                } else if !verified {
                    "reference estimator stop reached; independent fixed-point verification failed"
                        .into()
                } else if verification.is_none() {
                    "reference estimator stop reached; verification explicitly off".into()
                } else {
                    "reference estimator stop and independent fixed-point verification passed"
                        .into()
                }
            });
            if coarse {
                if trim_memory {
                    dft_solvers::ref_davidson::work::release_row_cache();
                }
                if failed {
                    return Err(ScfError::Numerical(exit_reason));
                }
                let density_g = state.fields.forward(comm, &step.density)?;
                let seed = super::ladder::Seed {
                    basis: state.basis,
                    grid: state.grid,
                    half: state.half,
                    rows: state.operator.local_rows().to_vec(),
                    orbitals: step.solve.eigenvectors,
                    indices: state.fields.reciprocal_indices().to_vec(),
                    density_g,
                };
                let delivered = StageRun::Coarse(
                    Box::new(seed),
                    super::ladder::Report {
                        kept_hx: kept.report(comm)?,
                        coarse_records: records,
                        coarse_setup: setup,
                        coarse_seconds: stage_start.elapsed().as_secs_f64(),
                        coarse_threshold_ry: profile.threshold_ry,
                        coarse_stopped: stopped,
                    },
                );
                boundary.emit(map_report).await;
                return Ok(delivered);
            }
            if trim_memory {
                kept.release_storage();
                dft_solvers::ref_davidson::work::release_row_cache();
                state.orbitals = dft_numerics::RealMatrix::zeros(0, bands).map_err(numerical)?;
                step.solve.hx = None;
                drop(field_cache);
                drop(mixer);
                drop(metric);
                drop(plan);
                drop(x);
                drop(y);
                drop(f);
                drop(next);
                drop(n);
                drop(next_planes);
            }
            let clock = Instant::now();
            records.last_mut().expect("current map").map_end_s = start.elapsed().as_secs_f64();
            let mut output = super::record::SetupTimings::default();
            let delivery_comm = refstart::TimedCollective::new(comm);
            let phase = Instant::now();
            let density = state.fields.gather(&delivery_comm, &step.density)?;
            output.record("density_gather", phase);
            let phase = Instant::now();
            let lazy_comm = owned_comm
                .clone()
                .filter(|_| crate::refscf::lazy::lazy_result_allowed());
            let (orbitals, lazy) = if let Some(comm) = lazy_comm {
                let placeholder = ComplexMatrix::zeros(0, bands).map_err(numerical)?;
                let lazy = crate::refscf::lazy::LazyOrbitals::from_owned_handoff(
                    std::mem::replace(
                        &mut step.solve.eigenvectors,
                        dft_numerics::RealMatrix::zeros(0, bands).map_err(numerical)?,
                    ),
                    state.operator.local_rows().to_vec(),
                    state.half.clone(),
                    comm,
                    Vec::new(),
                );
                (placeholder, Some(lazy))
            } else {
                (
                    gather_orbitals(&state, &delivery_comm, &step.solve.eigenvectors)?,
                    None,
                )
            };
            output.record("orbital_gather_unpack", phase);
            output.collective_s = delivery_comm.seconds();
            let phase = Instant::now();
            let point = KPoint {
                fractional: [0.; 3],
                weight: 1.,
            };
            let k = KPointSolution {
                point,
                lazy,
                solution: EigenSolution {
                    locking: None,
                    spectral_completeness: None,
                    eigenvalues: std::mem::take(&mut step.solve.eigenvalues),
                    eigenvectors: orbitals.clone(),
                    residuals: step.residuals,
                    converged: step.solve.converged,
                    buffer: occupations.iter().map(|&f| f == 0.).collect(),
                    residuals_certified: false,
                    iterations: dav,
                    history: Vec::new(),
                    operator_applications: columns,
                },
            };
            let identity = crate::force_identity::ModelDescription::new(
                cell,
                input.atoms,
                input.potentials,
                &options,
                dims,
                vec![point],
            );
            let result = ScfResult {
                atomic_start: AtomicStartReport {
                    atomic_orbitals: if prolonged_start {
                        0
                    } else {
                        state.memory.atomic_columns
                    },
                    initial_h_applications: if prolonged_start {
                        0
                    } else {
                        state.memory.atomic_columns.max(bands)
                    },
                    ..Default::default()
                },
                compat_mode: options.compat_mode,
                crystal_density_group: false,
                standard_stall_guard: false,
                standard_density_forcing: false,
                pw92_prefactor: options.pw92_prefactor,
                local_mesh: options.local_mesh,
                radial_even_rule: options.radial_even_rule,
                gamma_real_used: true,
                gamma_real_fallback_reason: None,
                force_identity: identity.with_policy(
                    crate::force_identity::OccupationIdentity::new(&[&occupations], None, None),
                ),
                standard: options.standard,
                certified: false,
                verification,
                spectral_completeness: None,
                eigenvalues: k.solution.eigenvalues.clone(),
                orbitals,
                occupations,
                gap_hartree: (bands > occupied).then(|| {
                    k.solution.eigenvalues[occupied] - k.solution.eigenvalues[occupied - 1]
                }),
                k_points: vec![k],
                energy: step.energy,
                energy_identity: None,
                charge: density.sum() * dv,
                density,
                history,
                ewald_breakdown: state.ewald,
                setup_time: start_time,
                total_time: start.elapsed(),
            };
            observe(
                Mem::Delivered,
                std::mem::size_of_val(result.density.as_slice().unwrap())
                    + std::mem::size_of_val(result.orbitals.as_slice())
                    + std::mem::size_of_val(result.k_points[0].solution.eigenvectors.as_slice()),
            );
            output.record("result_assembly", phase);
            output.total_s = clock.elapsed().as_secs_f64();
            let qe_clock = qe_scope.as_ref().map_or_else(Vec::new, |scope| {
                qe_clock_table(
                    comm.rank(),
                    &setup,
                    &output,
                    &records,
                    qe_mixing_calls,
                    scope.snapshot().without(qe_verification),
                )
            });
            return Ok(StageRun::Fine(Box::new(Run {
                start_report: state.report,
                symmetry: None,
                verification_k_call_indices: Vec::new(),
                memory_breakdown: memory_scope.report(comm)?,
                kept_hx: kept.report(comm)?,
                ladder: None,
                qe_clock,
                setup,
                output,
                full_diagnostics,
                outcome: ScfOutcome { converged, result },
                records,
                reference_stopped: first_reference_stop.is_some(),
                first_reference_stop,
                exit_reason,
                start_s: start_time.as_secs_f64(),
                output_s: clock.elapsed().as_secs_f64(),
                verification_calls,
                peak_payload_bytes_per_rank,
            })));
        }
        std::mem::swap(&mut x, &mut next);
        n.copy_from_slice(&next_planes);
        state.orbitals = step.solve.eigenvectors;
        records.last_mut().expect("current map").map_end_s = start.elapsed().as_secs_f64();
        completed = iteration;
        boundary.emit_state(map_report, checkpoint_state!()).await;
    }
    unreachable!("positive budget checked")
}
fn set_zero(plan: &DensityPlan, values: &mut [Complex64], charge: f64) {
    for (v, m) in values.iter_mut().zip(&plan.miller) {
        if *m == [0; 3] {
            *v = Complex64::new(charge, 0.);
        }
    }
}
fn global_min<K: Collective>(comm: &K, values: &[f64]) -> Result<f64, ScfError> {
    let mut local = vec![0.; comm.size()];
    local[comm.rank()] = values.iter().copied().fold(f64::MAX, f64::min);
    let mut all = local.clone();
    comm.tree_sum_f64(&local, &mut all).map_err(numerical)?;
    Ok(all.into_iter().fold(f64::MAX, f64::min))
}
#[allow(clippy::too_many_arguments)]
fn solve_map<K: Collective>(
    state: &mut StartState,
    comm: &K,
    cell: &Cell,
    functional: Functional,
    n: &[f64],
    occupations: &[f64],
    ethr: f64,
    first: bool,
    diagnostics: bool,
    rmm: Option<dft_solvers::ref_rmm::Options>,
    calls: &mut Vec<SolveRecord>,
    timings: &mut PhaseTimings,
    cache: &mut Option<(Vec<f64>, map::Fields)>,
    lean: bool,
    kept: &mut super::kept_hx::KeptHx,
) -> Result<map::Map, ScfError> {
    // QE's documented bounded re-entry: at most six calls, retry when more
    // than max(5, bands/4) roots remain. Keep every flag/count in the record;
    // SCF and the independent certificate still decide delivered convergence.
    let mut warm_start = false;
    for retry in 0..6 {
        if retry > 0 {
            kept.invalidate();
        }
        let mut step = map::evaluate(
            state,
            comm,
            cell,
            functional,
            n,
            occupations,
            ethr,
            first && retry == 0,
            diagnostics,
            cache,
            lean,
            rmm,
            kept,
        )?;
        warm_start |= step.warm_start;
        step.warm_start = warm_start;
        let r = step.solve.record.as_ref().expect("fresh solve accounting");
        let unconverged = step.solve.converged.iter().filter(|&&v| !v).count();
        calls.push(SolveRecord {
            ethr_ry: ethr,
            first_step: first && retry == 0,
            iterations: r.iterations.len(),
            active_bands: r.iterations.iter().map(|i| i.active_bands.len()).collect(),
            h_columns: r.h_columns,
            reduced_orders: std::iter::once(&r.setup)
                .chain(&r.iterations)
                .chain(std::iter::once(&r.finish))
                .filter(|i| i.reduced_solves > 0)
                .map(|i| i.reduced_order)
                .collect(),
            restarts: r.restarts,
            unconverged,
            ritz_changes_ha: std::mem::take(&mut step.solve.ritz_changes_ha),
            inner_clocks: (diagnostics && rmm.is_none()).then(|| {
                let mut stages = [0.; 8];
                for row in std::iter::once(&r.setup)
                    .chain(&r.iterations)
                    .chain(std::iter::once(&r.finish))
                {
                    for (sum, seconds) in stages.iter_mut().zip(row.seconds) {
                        *sum += seconds;
                    }
                }
                (r.total_seconds, stages)
            }),
            ritz_values_ha: std::mem::take(&mut step.solve.ritz_values_ha),
            band_thresholds_ha: step.solve.band_thresholds_ha.clone(),
            corrections: r
                .iterations
                .iter()
                .map(|i| CorrectionRecord {
                    active_bands: i.active_bands.clone(),
                    basis_before: i.basis_before,
                    h_columns: i.h_columns,
                    reduced_order: i.reduced_order,
                    repairs: i.repairs,
                    dropped_corrections: i.dropped_corrections,
                    restarts: i.restarts,
                })
                .collect(),
            repairs: r.iterations.iter().map(|i| i.repairs).sum(),
            dropped_corrections: r.iterations.iter().map(|i| i.dropped_corrections).sum(),
        });
        timings.potential += step.timings.potential;
        timings.diagonalization += step.timings.diagonalization;
        timings.density += step.timings.density;
        timings.energy += step.timings.energy;
        if rmm.is_some() || unconverged <= 5.max(occupations.len() / 4) {
            return Ok(step);
        }
        if retry == 5 {
            return Err(ScfError::Numerical(format!(
                "reference Davidson exhausted six calls: {unconverged} unconverged bands"
            )));
        }
        state.orbitals = step.solve.eigenvectors;
    }
    unreachable!()
}
fn gather_orbitals<K: Collective>(
    state: &StartState,
    comm: &K,
    x: &dft_numerics::RealMatrix,
) -> Result<ComplexMatrix, ScfError> {
    let n = state.basis.vectors().len();
    let mut local = vec![0.; n * x.cols()];
    for j in 0..x.cols() {
        for (i, &r) in state.operator.local_rows().iter().enumerate() {
            local[j * n + r] = x[(i, j)];
        }
    }
    let mut full = local.clone();
    comm.tree_sum_f64(&local, &mut full).map_err(numerical)?;
    observe(Mem::Gather, bytes(&local) + bytes(&full));
    if dft_numerics::switches::on(dft_numerics::switches::Id::REFSCF_MEMORY_TRIM) {
        drop(local);
    }
    let mut out = ComplexMatrix::zeros(n, x.cols()).map_err(numerical)?;
    for j in 0..x.cols() {
        state
            .half
            .unpack(&full[j * n..(j + 1) * n], out.column_mut(j))
            .map_err(numerical)?;
    }
    Ok(out)
}

/// Assemble outside numerical loops. Existing max/min reports remain independent.
fn qe_clock_table(
    rank: usize,
    setup: &super::record::SetupTimings,
    output: &super::record::SetupTimings,
    records: &[Record],
    mixing_calls: usize,
    counters: dft_solvers::ref_davidson::QeCounters,
) -> Vec<super::record::QeClock> {
    use super::record::QeClock;
    let mut rows = Vec::with_capacity(20);
    // These boundaries cannot be recovered from the owned files without
    // confusing unlike work. Keep explicit holes instead of proxies.
    for name in [
        "init_run",
        "electrons",
        "c_bands",
        "fftw",
        "rdiaghg",
        "davidson_remainder",
        "sum_band",
        "v_of_rho",
        "mix_rho",
    ] {
        rows.push(QeClock::unavailable(name, rank));
    }
    for (name, slot) in [
        ("h_psi", 0),
        ("calbec", 3),
        ("add_vuspsi", 4),
        ("regterg", 5),
    ] {
        rows.push(QeClock::measured(
            name,
            counters.seconds[slot],
            counters.calls[slot],
            rank,
            "exact",
        ));
    }
    rows.push(QeClock::measured(
        "vloc_psi",
        counters.seconds[1] - counters.seconds[2],
        counters.calls[1],
        rank,
        "exact",
    ));
    let sum = |get: fn(&Record) -> f64| records.iter().map(get).sum::<f64>();
    let maps = records.iter().map(|r| r.calls.len()).sum();
    let verify = sum(|r| r.verification_s);
    for (name, seconds, calls) in [
        ("baandi_setup", setup.total_s, 1),
        (
            "baandi_scf",
            records.last().map_or(0., |r| r.map_end_s)
                - records.first().map_or(0., |r| r.map_start_s)
                - verify,
            1,
        ),
        ("baandi_solve", sum(|r| r.solve_s), maps),
        ("baandi_density", sum(|r| r.density_s), maps),
        ("baandi_input_potential", sum(|r| r.potential_s), maps),
        ("baandi_mix", sum(|r| r.mixing_s), mixing_calls),
        (
            "baandi_reduced_standard",
            counters.seconds[7],
            counters.calls[7],
        ),
        (
            "baandi_davidson_remainder",
            counters.seconds[5] - counters.seconds[6] - counters.seconds[7],
            counters.calls[5],
        ),
        (
            "verification",
            verify,
            records
                .iter()
                .filter(|r| r.verification_start_s.is_some())
                .count(),
        ),
        ("output", output.total_s, 1),
    ] {
        rows.push(QeClock::measured(name, seconds, calls, rank, "baandi_only"));
    }
    rows
}

#[cfg(test)]
mod delivery_tests {
    use super::*;
    use crate::model::SpeciesPotential;
    use dft_numerics::switches::RunConfig;
    use std::sync::Arc;

    #[test]
    fn reference_lazy_delivery_preserves_science_and_materializes() {
        for ranks in [2, 3] {
            dft_mpi::simulate(ranks, |comm| {
                crate::parallel::set_communicator(Some(Arc::new(comm.clone())));
                let pp =
                    crate::pool_fixture::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
                let cell = dft_geometry::Cell::from_user_rows_bohr([
                    [7., 0., 0.],
                    [0.2, 7., 0.],
                    [0., 0.1, 7.],
                ])
                .unwrap();
                let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
                let potentials = [SpeciesPotential {
                    label: "Xx",
                    potential: &pp,
                }];
                let options = ScfOptions {
                    cutoff_hartree: 1.,
                    bands: Some(3),
                    max_iterations: 100,
                    gamma_real: crate::model::GammaReal::On,
                    ..Default::default()
                };
                let mut baseline: Option<ScfResult> = None;
                for (needed, trim) in [(false, "off"), (true, "off"), (false, "on"), (true, "on")] {
                    let config = RunConfig::from_pairs(&[
                        ("refscf_memory_trim", trim),
                        ("m2_lazy_result", "on"),
                    ])
                    .unwrap();
                    let _guard = config.enter();
                    crate::refscf::lazy::set_delivered_orbitals(needed);
                    let run = run(
                        StartInput {
                            cell: &cell,
                            atoms: &atoms,
                            potentials: &potentials,
                            options: &options,
                        },
                        Profile::default(),
                    )
                    .unwrap();
                    assert!(run.outcome.converged);
                    let mut r = run.outcome.result;
                    assert!(r.verification.as_ref().unwrap().passed());
                    assert_eq!(r.orbitals_delivered(), needed);
                    if !r.orbitals_delivered() {
                        assert_eq!(r.orbitals.rows(), 0);
                    }
                    r.materialize_orbitals().unwrap();
                    assert!(r.orbitals_delivered());
                    if let Some(base) = &baseline {
                        assert_eq!(
                            base.energy.total_hartree().to_bits(),
                            r.energy.total_hartree().to_bits()
                        );
                        assert_eq!(base.charge.to_bits(), r.charge.to_bits());
                        assert_eq!(base.eigenvalues, r.eigenvalues);
                        assert_eq!(base.occupations, r.occupations);
                        assert_eq!(base.history.len(), r.history.len());
                        assert_eq!(
                            base.density
                                .iter()
                                .flat_map(|v| v.to_le_bytes())
                                .collect::<Vec<_>>(),
                            r.density
                                .iter()
                                .flat_map(|v| v.to_le_bytes())
                                .collect::<Vec<_>>()
                        );
                        assert_eq!(
                            (base.orbitals.rows(), base.orbitals.cols()),
                            (r.orbitals.rows(), r.orbitals.cols())
                        );
                        for (a, b) in base.orbitals.as_slice().iter().zip(r.orbitals.as_slice()) {
                            assert_eq!(
                                (a.re.to_bits(), a.im.to_bits()),
                                (b.re.to_bits(), b.im.to_bits())
                            );
                        }
                    } else {
                        baseline = Some(r);
                    }
                }
                crate::refscf::lazy::set_delivered_orbitals(true);
                crate::parallel::set_communicator(None);
            });
        }
    }
}
