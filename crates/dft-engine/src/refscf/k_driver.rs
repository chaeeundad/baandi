//! Sequential-k reference driver. Density mixing, ethr, stopping and independent
//! verification retain the real reference loop's ordering. Gamma dispatch stays
//! in driver::run; this sibling owns only complex k lists.
use super::memory::{self, Component as Mem, bytes, observe};
use super::stepper::{Boundary, StepKind, StepReport, Stepper};
use super::{
    broyden::{Broyden, Diagnostics},
    density::DensityPlan,
    k_map as map,
    metric::{Metric, Serial},
    numerical,
    profile::{Profile, next_ethr_ry},
    record::{CorrectionRecord, Record, SolveRecord, VerificationCost},
};
use crate::Functional;
use crate::{
    model::*,
    refstart::{self, KStartState as StartState, StartInput},
};
use dft_geometry::Cell;
use dft_mpi::Collective;
use dft_physics::data::{Complex64, ComplexMatrix};
use dft_solvers::EigenSolution;
use std::time::{Duration, Instant};

use super::driver::{Run, VerificationPredictor, VerifyFeed, feed_attempt};
use dft_geometry::kpoints::{self, KPoint};
pub fn run(input: StartInput<'_>, profile: Profile, points: &[KPoint]) -> Result<Run, ScfError> {
    if points.len() == 1
        && points[0].fractional == [0.; 3]
        && input.options.gamma_real == crate::model::GammaReal::On
    {
        kpoints::validate(points).map_err(numerical)?;
        return super::driver::run(input, profile);
    }
    if let Some(comm) = crate::parallel::SharedComm::active() {
        run_with(input, profile, points, &comm)
    } else {
        run_with(input, profile, points, &Serial)
    }
}
pub fn run_with<K: Collective>(
    input: StartInput<'_>,
    profile: Profile,
    points: &[KPoint],
    comm: &K,
) -> Result<Run, ScfError> {
    run_inner(input, profile, points, comm, None, None)
}
/// Run the verified irreducible list and weights from legacy preparation.
pub fn run_prepared(
    input: StartInput<'_>,
    profile: Profile,
    prepared: &super::symmetry::Prepared,
) -> Result<Run, ScfError> {
    if let Some(comm) = crate::parallel::SharedComm::active() {
        run_prepared_with(input, profile, prepared, &comm)
    } else {
        run_prepared_with(input, profile, prepared, &Serial)
    }
}
pub fn run_prepared_with<K: Collective>(
    input: StartInput<'_>,
    profile: Profile,
    prepared: &super::symmetry::Prepared,
    comm: &K,
) -> Result<Run, ScfError> {
    run_inner(
        input,
        profile,
        prepared.points(),
        comm,
        Some(prepared),
        None,
    )
}
/// Fixed-cell seed; mixer and verification history are always fresh.
pub fn run_seeded_with<K: Collective>(
    input: StartInput<'_>,
    profile: Profile,
    points: &[KPoint],
    comm: &K,
    seed: super::delivered::KSeed,
) -> Result<Run, ScfError> {
    run_inner(input, profile, points, comm, None, Some(seed))
}
pub fn stepper_with<'a, K: Collective>(
    input: StartInput<'a>,
    profile: Profile,
    points: &'a [KPoint],
    comm: &'a K,
) -> Stepper<'a, Run> {
    Stepper::new(move |boundary| async move {
        run_continuation(input, profile, points, comm, None, None, &boundary).await
    })
}
fn run_inner<K: Collective>(
    input: StartInput<'_>,
    profile: Profile,
    points: &[KPoint],
    comm: &K,
    prepared: Option<&super::symmetry::Prepared>,
    seed: Option<super::delivered::KSeed>,
) -> Result<Run, ScfError> {
    Stepper::new(move |boundary| async move {
        run_continuation(input, profile, points, comm, prepared, seed, &boundary).await
    })
    .run_to_end()
}
#[allow(clippy::too_many_arguments)]
pub(super) async fn run_continuation<K: Collective>(
    input: StartInput<'_>,
    profile: Profile,
    points: &[KPoint],
    comm: &K,
    prepared: Option<&super::symmetry::Prepared>,
    seed: Option<super::delivered::KSeed>,
    boundary: &Boundary,
) -> Result<Run, ScfError> {
    kpoints::validate(points).map_err(numerical)?;
    if seed.is_none()
        && prepared.is_none()
        && points.len() == 1
        && points[0].fractional == [0.; 3]
        && input.options.gamma_real == crate::model::GammaReal::On
    {
        return super::driver::run_inner(input, profile, comm, None, boundary).await;
    }
    let start = Instant::now();
    let stage_start = Instant::now();
    let memory_scope = memory::Scope::begin();
    // Resolve once from this run's scope, never a process-cached profile flag.
    use dft_numerics::switches::{Id, RunConfig};
    let config = RunConfig::current();
    let full_diagnostics = config.on(Id::DAVIDSON_STAGES);
    let trim_memory = config.on(Id::REFSCF_MEMORY_TRIM);
    let rmm: Option<dft_solvers::ref_rmm::Options> = None;
    if config.on(Id::REFSCF_RMM) || config.on(Id::REFSCF_LADDER) {
        return Err(ScfError::InvalidInput(
            "reference k lists currently require Davidson and final-cutoff start (RMM/ladder unavailable)",
        ));
    }
    let qe_scope = full_diagnostics.then(dft_solvers::ref_davidson::QeClockScope::begin);
    let mut qe_verification = dft_solvers::ref_davidson::QeCounters::default();

    let mut profile = profile.validate()?;
    let o = input.options;
    // Prepared symmetry must not supply an implicit grid for a frozen segment.
    refstart::validate_strain_grid(o)?;
    if o.standard.is_none()
        || (o.crystal_density_group && prepared.is_none())
        || o.max_iterations == 0
    {
        return Err(ScfError::InvalidInput(
            "reference k SCF requires standard class, prepared symmetry when requested, and positive iteration budget",
        ));
    }
    let mut options = o.clone();
    options.standard = Some(StandardProfile {
        threshold_ha: profile.threshold_ry / 2.,
    });
    if let Some(p) = prepared {
        p.validate(&input, points)?;
        options.grid_override = Some(p.grid(&input)?);
    }
    let cell = input.cell;
    let mut state = refstart::start_kpoints(
        StartInput {
            options: &options,
            ..input
        },
        points,
        comm,
    )?;
    if let Some(seed) = seed {
        seed.inject(&input, &mut state)?;
        state.common.report = refstart::StartReport {
            density: refstart::DensityStartKind::Supplied,
            orbitals: refstart::OrbitalStartKind::Supplied,
        };
    }
    let mut symmetry = prepared
        .map(|p| {
            super::symmetry::Workspace::new(&input, p, state.grid.dims(), &state.fields, 1, comm)
        })
        .transpose()?;
    if let Some(s) = &mut symmetry {
        let initial = [state.density.clone()];
        s.check_start(&mut state.fields, comm, &initial)?;
        state.potential_symmetry = s.take_potential_projector();
    }
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
        o.strain_reference.as_deref(),
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
    // Named-array envelope: retained common setup plus all k states and one
    // active complex workspace. Backend packing, allocator and FFT metadata
    // remain excluded, as in the real reference memory contract.
    let mut retained = state.memory.principal_peak_bound_bytes;
    let mut active = 0;
    for block in &state.blocks {
        let rows = block.operator.local_rows().len();
        let cap = (2 * bands).min(block.basis.vectors().len());
        retained += 8 * block.operator.plane_points()
            + block.operator.payload_bytes()
            + 16 * rows * bands * 8
            + 16 * block.basis.vectors().len() * bands * 2;
        active = active.max(16 * rows * (6 * cap + 2 * bands) + 16 * cap * cap * 24);
    }
    let peak = retained + active + 16 * x.len() * (2 * profile.history + 12) + 8 * n.len() * 28;
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
    let mut stalls = Vec::with_capacity(o.max_iterations);
    let mut map_diagnostics = Vec::with_capacity(o.max_iterations);
    let mut history = Vec::<IterationRecord>::with_capacity(o.max_iterations);
    let mut records = Vec::with_capacity(o.max_iterations);
    let mut ethr = 1e-2;
    let mut last_dr2 = f64::INFINITY;
    let mut verification =
        (o.verification == Verification::FixedPoint).then(VerificationReport::default);
    let mut verification_calls = Vec::new();
    let mut verification_k_call_indices = Vec::new();
    let mut first_reference_stop = None;
    let mut pending: Option<(map::Map, f64, usize)> = None;
    let mut predictor = VerificationPredictor::default();
    let feed = VerifyFeed::resolve(config);
    let skip_run = super::driver::verify_skip_run(config);
    let response_predictor = super::driver::verify_response_predictor(config);
    let permap = config.on(Id::REFSCF_PERMAP);
    let mut field_cache = None;

    let mut kept = super::k_kept::KeptHx::new(true);
    let mut completed = 0usize;
    // R8b: closure creation borrows locals only. State visits occur exclusively
    // on an explicit checkpoint poll; ordinary maps allocate nothing here.
    macro_rules! checkpoint_state {
        () => {
            |request: super::stepper::Visit<'_, '_>| -> Result<usize, ScfError> {
                let io = match request {
                    super::stepper::Visit::Checkpoint(io) => io,
                    super::stepper::Visit::Policy(patch) => {
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
                super::checkpoint::header(io, &input, profile, comm, "complex-k", dims)?;
                super::checkpoint::symmetry(io, &mut symmetry)?;
                io.same(&plan.g2)?;
                io.value(&mut completed)?;
                io.planes(&mut n)?;
                io.reciprocal(&mut x)?;
                mixer.checkpoint(io, &metric)?;
                super::checkpoint::k_state(io, &mut state)?;
                kept.checkpoint(io, &state)?;
                io.value(&mut verification_k_call_indices)?;
                io.value(&mut ethr)?;
                io.value(&mut last_dr2)?;
                io.value(&mut history)?;
                io.value(&mut records)?;
                io.value(&mut verification)?;
                io.value(&mut verification_calls)?;
                io.value(&mut first_reference_stop)?;
                io.value(&mut predictor.observed)?;
                io.value(&mut predictor.skipped_last)?;
                if skip_run {
                    io.value(&mut predictor.last_skip_dr2)?;
                }
                io.value(&mut field_cache)?;
                io.value(&mut pending)?;
                super::checkpoint::cached_fields(&field_cache, n.len())?;
                if let Some((map, _, _)) = &pending {
                    map.validate_checkpoint(&state, n.len())?;
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
        if let Some(s) = &mut symmetry {
            s.iteration = iteration;
            s.verification = false;
        }
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
        if iteration == 1 {
            // GK10B (`K_START_KEPT_HX`): seed the kept HX with the start's
            // H X of the start Ritz orbitals; `prepare` hits only while the
            // orbitals are still those exact vectors (not after a seed) and
            // corrects the potential change exactly.
            for (k, (x, hx, potential)) in
                std::mem::take(&mut state.start_hx).into_iter().enumerate()
            {
                kept.retain(k, &state, &x, &hx, &potential);
            }
        }
        let requested = ethr;
        let mut redo = false;
        let mut reported_columns = 0;
        let mut calls = Vec::new();
        let mut call_indices = Vec::new();
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
                    &mut call_indices,
                    &mut timings,
                    &mut field_cache,
                    true,
                    &mut kept,
                    &mut symmetry,
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
                report.stage = "fine";
                reported_columns = report.h_columns;
                boundary.emit(report).await;
                ethr = next_ethr_ry(ethr, dr2, state.electrons);
                redo = true;
                // GK10B (`K_REDO_KEPT_HX`): the redo solves the same input
                // potential from the coarse solve's X, whose HX is kept
                // exactly (zero potential change); off applies H again.
                if !config.on(Id::K_REDO_KEPT_HX) {
                    kept.invalidate();
                }
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
            if response_predictor && predictor.observed.is_none() {
                // GV1: no attempt yet; the secant output-response estimate.
                mixer
                    .predict_output_residual(comm, &mut metric, &x, &f)?
                    .map(|v| v / 2.)
            } else {
                predictor.predict(dr2)
            }
        } else {
            None
        };
        let skipped = stopped
            && verification.is_some()
            && predictor.skip_with(predicted_error_ha, profile.threshold_ry / 2., dr2, skip_run);
        if ((stopped && !skipped) || iteration == o.max_iterations) && step.energy.hartree.is_nan()
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
        map_report.stage = "fine";
        map_report.reused = reused_from.is_some();
        let mut verification_warm = None;
        let mut actual_error_ha = None;
        let mut verification_s = 0.;
        let (mut verification_start_s, mut verification_end_s) = (None, None);
        let mut fed = None;
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
            if let Some(s) = &mut symmetry {
                s.verification = true;
            }
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
                &mut verification_k_call_indices,
                &mut t,
                &mut field_cache,
                false,
                &mut kept,
                &mut symmetry,
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
                strict_residuals_per_k: check
                    .solve
                    .per_k
                    .iter()
                    .map(|k| k.residuals[..occupied].to_vec())
                    .collect(),
                buffer_residuals_per_k: check
                    .solve
                    .per_k
                    .iter()
                    .map(|k| k.residuals[occupied..].to_vec())
                    .collect(),
                density_residual: (dv * l2[0]).sqrt(),
                occupied_eigenvalue_change: check
                    .solve
                    .per_k
                    .iter()
                    .zip(&step.solve.per_k)
                    .flat_map(|(a, b)| a.eigenvalues.iter().zip(&b.eigenvalues).take(occupied))
                    .map(|(a, b)| (a - b).abs())
                    .fold(0., f64::max),
                davidson_iterations: current_calls.iter().map(|c| c.iterations).sum(),
                operator_applications: current_calls.iter().map(|c| c.h_columns).sum(),
                davidson_iterations_per_k: (0..points.len())
                    .map(|k| {
                        current_calls
                            .iter()
                            .enumerate()
                            .filter(|(i, _)| verification_k_call_indices[first_call + *i] == k)
                            .map(|(_, c)| c)
                            .map(|c| c.iterations)
                            .sum()
                    })
                    .collect(),
                operator_applications_per_k: (0..points.len())
                    .map(|k| {
                        current_calls
                            .iter()
                            .enumerate()
                            .filter(|(i, _)| verification_k_call_indices[first_call + *i] == k)
                            .map(|(_, c)| c)
                            .map(|c| c.h_columns)
                            .sum()
                    })
                    .collect(),
                timings: t,
            };
            let bound = StandardProfile::residual_energy_estimate(
                points
                    .iter()
                    .zip(&check.solve.per_k)
                    .map(|(p, k)| (p.weight, k.residuals.as_slice())),
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
            map_report.stage = "fine";
            if !attempt.passed && feed != VerifyFeed::Off {
                // GV1: keep F(y) for the mixer (G = 0 carries no residual).
                set_zero(&plan, &mut checked, 0.);
                fed = Some(checked);
            }

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
            let clock = Instant::now();
            match mixer.update(comm, &mut metric, &x, &f, &mut next) {
                Ok(d) => diagnostics = d,
                Err(e) => failure = Some(e.to_string()),
            }
            if failure.is_none()
                && let Some(g) = &fed
            {
                match feed_attempt(feed, &mut mixer, comm, &mut metric, &y, g, &mut next) {
                    Ok(d) => diagnostics = d,
                    Err(e) => failure = Some(e.to_string()),
                }
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
        if failure.is_some() && step.energy.hartree.is_nan() {
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
            k_call_indices: call_indices,
            warm_start: step.warm_start && !calls.is_empty(),
            stage: "fine",
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
        let mut sampled = super::stall::MapDiagnostic {
            iteration,
            ..Default::default()
        };
        if occupied > 0 && occupied < bands {
            let high = step
                .solve
                .per_k
                .iter()
                .map(|b| b.eigenvalues[occupied - 1])
                .fold(f64::NEG_INFINITY, f64::max);
            let low = step
                .solve
                .per_k
                .iter()
                .map(|b| b.eigenvalues[occupied])
                .fold(f64::INFINITY, f64::min);
            sampled = super::stall::MapDiagnostic {
                iteration,
                channel_gap_ha: [Some(low - high), None],
                channel_edge_residual_ha: [
                    Some(
                        step.solve
                            .per_k
                            .iter()
                            .map(|b| b.residuals[occupied - 1] + b.residuals[occupied])
                            .fold(0., f64::max),
                    ),
                    None,
                ],
                ..Default::default()
            };
        }
        map_diagnostics.push(sampled);
        if super::stall::detected(records.iter().rev().map(|r| (r.dr2_ry, r.mixing.condition))) {
            let diagnostic = super::stall::Diagnostic {
                iteration,
                channel_gap_ha: sampled.channel_gap_ha,
                orbital_rotation_fraction: None,
            };
            map_report.stall = Some(diagnostic);
            stalls.push(diagnostic);
        }
        map_report.safe_point = failure.is_none();
        if finish || failure.is_some() || iteration == o.max_iterations {
            boundary.publish(map_report);
            let verified = verification.as_ref().is_none_or(VerificationReport::passed);
            let converged = stopped && verified && failure.is_none();
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
            if trim_memory {
                kept.release_storage();
                dft_solvers::ref_davidson::work::release_row_cache();
                state.orbitals = refstart::KOrbitals(Vec::new());
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
            let mut k_points = Vec::with_capacity(points.len());
            // GK18: a pooled run whose caller does not read the orbitals keeps
            // them with their owners (materialized on demand, same bits).
            let lazy = if crate::refscf::lazy::lazy_result_allowed() && state.pool.is_some() {
                super::k_pool::lazy_orbitals(
                    &state,
                    comm,
                    std::mem::take(&mut step.solve.eigenvectors.0),
                    bands,
                )?
            } else {
                None
            };
            // GK9 pools: the owner of each k broadcasts its full-row block.
            let (delivered, lazy): (Vec<_>, Vec<_>) = match lazy {
                Some(handles) => handles
                    .into_iter()
                    .map(|h| ComplexMatrix::zeros(0, bands).map(|z| (z, Some(h))))
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(numerical)?
                    .into_iter()
                    .unzip(),
                None => super::k_pool::deliver_orbitals(
                    &state,
                    &delivery_comm,
                    &step.solve.eigenvectors.0.iter().collect::<Vec<_>>(),
                )?
                .into_iter()
                .map(|x| (x, None))
                .unzip(),
            };
            for (((meta, orbitals), lazy), point) in
                step.solve.per_k.iter().zip(delivered).zip(lazy).zip(points)
            {
                k_points.push(KPointSolution {
                    point: *point,
                    lazy,
                    solution: EigenSolution {
                        locking: None,
                        spectral_completeness: None,
                        eigenvalues: meta.eigenvalues.clone(),
                        eigenvectors: orbitals,
                        residuals: meta.residuals.clone(),
                        converged: meta.converged.clone(),
                        buffer: occupations.iter().map(|&f| f == 0.).collect(),
                        residuals_certified: false,
                        iterations: meta.calls.iter().map(|c| c.iterations).sum(),
                        history: Vec::new(),
                        operator_applications: meta.calls.iter().map(|c| c.h_columns).sum(),
                    },
                });
            }
            let orbitals = k_points[0].solution.eigenvectors.clone();
            // GK18B: a lazy delivery forms handles only (no collective).
            let lazy_delivery = k_points.iter().any(|k| k.lazy.is_some());
            output.record(
                if lazy_delivery {
                    "complex_orbital_lazy"
                } else {
                    "complex_orbital_gather"
                },
                phase,
            );
            output.collective_s = delivery_comm.seconds();
            let phase = Instant::now();
            let mut identity = crate::force_identity::ModelDescription::new(
                cell,
                input.atoms,
                input.potentials,
                &options,
                dims,
                points.to_vec(),
            );
            if let Some(p) = prepared {
                identity.symmetry = p.operations(o.crystal_density_group).to_vec();
            }
            let result = ScfResult {
                atomic_start: AtomicStartReport {
                    atomic_orbitals: state.memory.atomic_columns,
                    initial_h_applications: state.initial_h_columns,
                    ..Default::default()
                },
                compat_mode: options.compat_mode,
                crystal_density_group: o.crystal_density_group,
                standard_stall_guard: false,
                standard_density_forcing: false,
                pw92_prefactor: options.pw92_prefactor,
                local_mesh: options.local_mesh,
                radial_even_rule: options.radial_even_rule,
                gamma_real_used: false,
                gamma_real_fallback_reason: None,
                force_identity: identity.with_policy(
                    crate::force_identity::OccupationIdentity::new(
                        &vec![occupations.as_slice(); points.len()],
                        None,
                        None,
                    ),
                ),
                standard: options.standard,
                certified: false,
                verification,
                spectral_completeness: None,
                eigenvalues: k_points[0].solution.eigenvalues.clone(),
                orbitals,
                occupations,
                gap_hartree: (bands > occupied).then(|| {
                    let valence = k_points
                        .iter()
                        .map(|k| k.solution.eigenvalues[occupied - 1])
                        .fold(f64::NEG_INFINITY, f64::max);
                    let conduction = k_points
                        .iter()
                        .map(|k| k.solution.eigenvalues[occupied])
                        .fold(f64::INFINITY, f64::min);
                    conduction - valence
                }),
                k_points,
                energy: step.energy,
                energy_identity: None,
                charge: density.sum() * dv,
                density,
                history,
                ewald_breakdown: state.common.ewald,
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
            let qe_clock = Vec::new();
            return Ok(Run {
                stalls,
                map_diagnostics,
                start_report: state.common.report,
                symmetry: symmetry.map(|s| s.report),
                verification_k_call_indices,
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
            });
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
    _rmm: Option<dft_solvers::ref_rmm::Options>,
    calls: &mut Vec<SolveRecord>,
    indices: &mut Vec<usize>,
    timings: &mut PhaseTimings,
    cache: &mut Option<(Vec<f64>, super::map::Fields)>,
    lean: bool,
    kept: &mut super::k_kept::KeptHx,
    symmetry: &mut Option<super::symmetry::Workspace>,
) -> Result<map::Map, ScfError> {
    let step = map::evaluate(
        state,
        comm,
        cell,
        functional,
        n,
        occupations,
        ethr,
        first,
        diagnostics,
        cache,
        lean,
        kept,
        symmetry,
    )?;
    for (index, k) in step.solve.per_k.iter().enumerate() {
        indices.extend(std::iter::repeat_n(index, k.calls.len()));
        calls.extend(k.calls.iter().cloned());
    }
    timings.potential += step.timings.potential;
    timings.diagonalization += step.timings.diagonalization;
    timings.density += step.timings.density;
    timings.energy += step.timings.energy;
    Ok(step)
}
pub(super) fn gather_orbitals<K: Collective>(
    block: &refstart::KBlock,
    comm: &K,
    x: &ComplexMatrix,
) -> Result<ComplexMatrix, ScfError> {
    let n = block.basis.vectors().len();
    let mut local = vec![0.; 2 * n * x.cols()];
    for j in 0..x.cols() {
        for (i, &g) in block.operator.local_rows().iter().enumerate() {
            local[2 * (g + j * n)] = x[(i, j)].re;
            local[2 * (g + j * n) + 1] = x[(i, j)].im;
        }
    }
    let mut all = local.clone();
    comm.tree_sum_f64(&local, &mut all).map_err(numerical)?;
    ComplexMatrix::new(
        n,
        x.cols(),
        all.as_chunks::<2>()
            .0
            .iter()
            .map(|p| Complex64::new(p[0], p[1]))
            .collect(),
    )
    .map_err(|e| numerical(e).into())
}
pub(super) fn solve_record(
    solve: &mut dft_solvers::ref_davidson_complex::Solution,
    _k: usize,
    ethr: f64,
    first: bool,
    diagnostics: bool,
) -> SolveRecord {
    let r = &solve.record;
    let unconverged = solve.converged.iter().filter(|&&v| !v).count();
    SolveRecord {
        ethr_ry: ethr,
        first_step: first,
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
        ritz_changes_ha: std::mem::take(&mut solve.ritz_changes_ha),
        inner_clocks: diagnostics.then(|| {
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
        ritz_values_ha: std::mem::take(&mut solve.ritz_values_ha),
        band_thresholds_ha: solve.band_thresholds_ha.clone(),
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
    }
}
