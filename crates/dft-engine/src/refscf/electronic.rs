//! Reference SCF for scalar FD and collinear fixed/FD occupations.
//! The fixed unpolarized entry points never call this module.
use super::{
    broyden::{Broyden, ConditioningDiagnostics, Diagnostics},
    density::DensityPlan,
    electronic_map::{self as map, Map},
    k_kept::KeptHx,
    metric::{Metric, Serial},
    profile::{Profile, next_ethr_ry},
    record::SolveRecord,
};
use crate::{
    Functional, ScfError, Verification, VerificationAttempt, VerificationReport,
    collinear::SpinSector,
    energy::EnergyAudit,
    finite_temperature::ThermodynamicLedger,
    refstart::{self, StartInput},
};
use dft_geometry::kpoints::KPoint;
use dft_mpi::Collective;
use dft_numerics::{
    linalg::ComplexMatrix,
    switches::{Id, RunConfig},
};
use dft_physics::{
    data::{Array3, Complex64 as C},
    occupations::{Constraint, Evaluation},
};
use std::time::Instant;

#[derive(Clone, Debug)]
pub enum Policy {
    FermiDirac {
        theta_ha: f64,
        band_cap: usize,
    },
    FixedSpin(SpinSector),
    SpinFermiDirac {
        sector: SpinSector,
        theta_ha: f64,
        constraint: Constraint,
        band_cap: usize,
    },
}
impl Policy {
    pub fn channels(&self) -> usize {
        if matches!(self, Self::FermiDirac { .. }) {
            1
        } else {
            2
        }
    }
    pub(super) fn finite(&self) -> Option<(f64, Constraint)> {
        match self {
            Self::FermiDirac { theta_ha, .. } => Some((*theta_ha, Constraint::Total)),
            Self::SpinFermiDirac {
                theta_ha,
                constraint,
                ..
            } => Some((*theta_ha, *constraint)),
            _ => None,
        }
    }
    pub(super) fn fixed_count(&self, s: usize) -> usize {
        match self {
            Self::FixedSpin(sector) => [sector.n_up, sector.n_down][s],
            _ => 0,
        }
    }
    pub(super) fn strict_bands(&self, s: usize, bands: usize) -> usize {
        if self.finite().is_some() {
            bands
        } else {
            (self.fixed_count(s) + 1).min(bands)
        }
    }
    fn target_moment(&self) -> Option<f64> {
        match self {
            Self::FixedSpin(s) => Some(s.moment()),
            Self::SpinFermiDirac {
                constraint: Constraint::FixedMoment(m),
                ..
            } => Some(*m),
            _ => None,
        }
    }
}
#[derive(Debug)]
pub struct BandBlock {
    pub point: KPoint,
    pub eigenvalues_ha: Vec<f64>,
    pub occupations: Vec<f64>,
    pub residuals_ha: Vec<f64>,
    pub orbitals: ComplexMatrix,
}
pub struct Iteration {
    pub iteration: usize,
    pub dr2_ry: f64,
    pub ethr_ry: f64,
    pub first_iteration_redo: bool,
    pub energy_ha: f64,
    pub entropy_term_ha: f64,
    pub free_energy_ha: f64,
    pub moment_mu_b: f64,
    pub residual_guard_ha: f64,
    pub mixing: Diagnostics,
    pub history_conditioning: ConditioningDiagnostics,
    pub mixing_fraction: f64,
    /// Channel/k/call order; includes first-step redo and bounded retries.
    pub calls: Vec<(usize, usize, SolveRecord)>,
}
pub struct Run {
    pub converged: bool,
    pub density: Vec<Array3<f64>>,
    pub blocks: Vec<Vec<BandBlock>>,
    pub energy: EnergyAudit,
    pub thermodynamics: ThermodynamicLedger,
    pub occupations: Option<Evaluation>,
    pub channel_electrons: Vec<f64>,
    pub moment_mu_b: f64,
    pub absolute_moment_mu_b: f64,
    pub history: Vec<Iteration>,
    pub verification: Option<VerificationReport>,
    pub verification_calls: Vec<(usize, usize, SolveRecord)>,
    pub kept_hx: Vec<super::kept_hx::Report>,
    pub setup_s: f64,
    pub total_s: f64,
}
/// Opt-in diagnostic records retained even when the run rejects a history.
/// Storage is reserved at setup; adapters inspect/serialize only after return.
#[derive(Clone, Debug, Default)]
pub struct MapTrace {
    pub iteration: usize,
    pub verification: bool,
    pub redo: bool,
    pub estimator_terms_ry: [f64; 3],
    pub free_energy_ha: f64,
    pub moment_mu_b: f64,
    pub chemical_potentials_ha: [Option<f64>; 2],
    pub electron_residual: f64,
    pub moment_residual: Option<f64>,
    pub occupations: Vec<f64>,
    pub condition: Option<f64>,
    pub history_pair_norm2_ry: [f64; 25],
    pub history_max_cosine: f64,
    pub stiffness_ha: f64,
    pub dropped_pairs: usize,
    pub dropped_conditions: [f64; 25],
}
#[derive(Default)]
pub struct Trace {
    records: Vec<MapTrace>,
    used: usize,
}
impl Trace {
    pub fn records(&self) -> &[MapTrace] {
        &self.records[..self.used]
    }
    fn prepare(&mut self, maps: usize, occupations: usize) {
        self.used = 0;
        self.records = (0..maps)
            .map(|_| MapTrace {
                occupations: vec![0.; occupations],
                ..Default::default()
            })
            .collect();
    }
    #[allow(clippy::too_many_arguments)]
    fn record<K: Collective>(
        &mut self,
        comm: &K,
        metric: &mut Metric,
        g2: &[f64],
        f: &[C],
        step: &Map,
        iteration: usize,
        verification: bool,
        redo: bool,
    ) -> Result<(), ScfError> {
        let record = &mut self.records[self.used];
        self.used += 1;
        record.iteration = iteration;
        record.verification = verification;
        record.redo = redo;
        record.estimator_terms_ry = if step.density.len() == 2 {
            metric.spin_terms(comm, f, g2)?
        } else {
            [0.; 3]
        };
        record.free_energy_ha = step.free_energy();
        record.moment_mu_b = step.moment;
        if let Some(e) = step.evaluation {
            record.chemical_potentials_ha = e.roots.map(|r| r.map(|r| r.mu_ha));
            record.electron_residual = e.electron_residual;
            record.moment_residual = e.moment_residual;
        }
        for (a, b) in record
            .occupations
            .iter_mut()
            .zip(step.blocks.iter().flatten().flat_map(|b| &b.occupations))
        {
            *a = *b;
        }
        Ok(())
    }
}
/// Gamma uses the complex solver as required by ONEPATH-3.
pub fn run(
    input: StartInput<'_>,
    profile: Profile,
    points: &[KPoint],
    policy: Policy,
) -> Result<Run, ScfError> {
    if let Some(comm) = crate::parallel::SharedComm::active() {
        run_with(input, profile, points, policy, &comm)
    } else {
        run_with(input, profile, points, policy, &Serial)
    }
}
pub fn run_with<K: Collective>(
    input: StartInput<'_>,
    profile: Profile,
    points: &[KPoint],
    policy: Policy,
    comm: &K,
) -> Result<Run, ScfError> {
    run_with_trace(input, profile, points, policy, comm, None)
}
/// Diagnostic entry; `None` retains the ordinary driver arithmetic and storage.
pub fn run_with_trace<K: Collective>(
    input: StartInput<'_>,
    profile: Profile,
    points: &[KPoint],
    policy: Policy,
    comm: &K,
    mut trace: Option<&mut Trace>,
) -> Result<Run, ScfError> {
    let clock = Instant::now();
    let profile = profile.validate()?;
    let opts = input.options;
    let config = RunConfig::current();
    if opts.standard.is_none() || opts.max_iterations == 0 || opts.crystal_density_group {
        return Err(ScfError::InvalidInput(
            "reference spin/FD requires standard profile, positive budget, no crystal symmetry",
        ));
    }
    if config.on(Id::REFSCF_RMM)
        || config.on(Id::REFSCF_LADDER)
        || config.on(Id::REFSCF_VERIFY_ENERGY)
        || config.on(Id::REFSCF_VERIFY_REUSE)
    {
        return Err(ScfError::InvalidInput(
            "reference spin/FD supports Davidson and fixed-point verification only; RMM, ladder, energy verification and failed-map reuse unavailable",
        ));
    }
    for id in [
        Id::REFSCF_VERIFY_GATE,
        Id::REFSCF_VERIFY_GATE_FIRST,
        Id::REFSCF_LEAN_ENERGY,
    ] {
        config.note_path(
            id,
            "inactive",
            Some("reference spin/FD uses complete map audits and verifies each stop"),
        );
    }
    if policy.channels() == 2 && opts.functional != Functional::LdaPw {
        return Err(ScfError::InvalidInput(
            "reference spin XC requires lda-pw: legacy polarized lda-pz unavailable; distributed PBE has no two-channel evaluator",
        ));
    }
    let electrons = input.atoms.iter().try_fold(0., |n, a| {
        input
            .potentials
            .iter()
            .find(|p| p.label == a.species)
            .map(|p| n + p.potential.z_valence())
            .ok_or(ScfError::InvalidInput("unknown species"))
    })?;
    let mut options = opts.clone();
    options.gamma_real = crate::GammaReal::Off;
    let bands = match &policy {
        Policy::FermiDirac { band_cap, .. } | Policy::SpinFermiDirac { band_cap, .. } => *band_cap,
        Policy::FixedSpin(sector) => options.bands.unwrap_or(
            sector
                .n_up
                .max(sector.n_down)
                .checked_add(options.extra_bands)
                .ok_or(ScfError::InvalidInput("band overflow"))?,
        ),
    };
    if let Some((theta, constraint)) = policy.finite() {
        crate::scf::finite_temperature::validate_options(&options, theta, bands)?;
        if !theta.is_finite()
            || theta <= 0.
            || bands == 0
            || options.bands.is_some_and(|b| b > bands)
            || electrons >= 2. * bands as f64
        {
            return Err(ScfError::InvalidInput(
                "occupations: positive theta and band cap with empty-state capacity required",
            ));
        }
        if let Constraint::FixedMoment(m) = constraint
            && (!m.is_finite()
                || m.abs() >= electrons
                || (electrons + m.abs()) / 2. >= bands as f64)
        {
            return Err(ScfError::InvalidInput(
                "occupations.electron_target: fixed moment requires interior channel targets",
            ));
        }
    }
    let seed_moment = match &policy {
        Policy::FixedSpin(sector) | Policy::SpinFermiDirac { sector, .. } => {
            sector.validate(electrons)?;
            if bands <= sector.n_up.max(sector.n_down) {
                return Err(ScfError::InvalidInput(
                    "spin.constraint: each channel requires an empty band",
                ));
            }
            policy
                .target_moment()
                .unwrap_or(sector.initial_moment_seed.unwrap_or(0.))
        }
        _ => 0.,
    };
    options.bands = Some(bands);
    let mut states = (0..policy.channels())
        .map(|_| {
            refstart::start_electronic(
                StartInput {
                    options: &options,
                    ..input
                },
                points,
                comm,
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let volume = input.cell.volume_bohr3();
    let dims = states[0].grid.dims();
    let dv = volume / dims.iter().product::<usize>() as f64;
    // Legacy collinear starts from uniform channel counts. A positive seed
    // permits its common feasibility step without clipping atomic-density zeros.
    let mut n: Vec<Vec<f64>> = states
        .iter()
        .enumerate()
        .map(|(s, state)| {
            if policy.channels() == 1 {
                state.density.clone()
            } else {
                vec![
                    (electrons + if s == 0 { seed_moment } else { -seed_moment }) / (2. * volume);
                    state.density.len()
                ]
            }
        })
        .collect();
    let mut plan = DensityPlan::with_fft(
        comm,
        input.cell,
        dims,
        4. * opts.cutoff_hartree,
        4. * opts.cutoff_hartree,
        states[0].fields.shared_fft(),
    )?;
    let mut metric = if policy.channels() == 2 {
        Metric::spin(comm, volume, &plan.g2, &plan.smooth, 350)?
    } else {
        Metric::new(volume, &plan.g2, &plan.smooth, comm.size(), 350)?
    };
    let mut mixer = if policy.channels() == 2 {
        Broyden::spin(profile, &plan.g2, &plan.smooth, volume, electrons)?
    } else {
        Broyden::new(profile, &plan.g2, &plan.smooth, volume, electrons)?
    };
    let spin_fd = matches!(policy, Policy::SpinFermiDirac { .. });
    let stiffness = spin_fd && config.on(Id::REFSCF_SPIN_STIFFNESS);
    if spin_fd && config.on(Id::REFSCF_HISTORY_CONDITIONING) {
        mixer.condition_history();
    }
    let mut kept: Vec<_> = (0..policy.channels())
        .map(|_| KeptHx::new(config.on(Id::REFSCF_KEPT_HX)))
        .collect();
    let mut x = transform(&mut plan, comm, &n, &policy, electrons, volume)?;
    let mut workspace = map::Workspace::new(&states[0], comm, input.cell, policy.channels() == 2)?;
    if let Some(trace) = trace.as_deref_mut() {
        trace.prepare(
            2 * opts.max_iterations + 1,
            policy.channels() * points.len() * bands,
        );
    }
    let setup_s = clock.elapsed().as_secs_f64();
    let mut history = Vec::new();
    let mut verification =
        (opts.verification == Verification::FixedPoint).then(VerificationReport::default);
    let mut verification_calls = Vec::new();
    let mut ethr = 1e-2;
    let mut previous = f64::INFINITY;
    let guard = config.on(Id::STANDARD_RESIDUAL_GUARD);
    let diagnostic = config.on(Id::DAVIDSON_STAGES);
    for iteration in 1..=opts.max_iterations {
        if iteration > 1 {
            ethr = next_ethr_ry(
                if iteration == 2 { 1e-2 } else { ethr },
                previous,
                electrons,
            );
        }
        let mut calls = Vec::new();
        let mut redo = false;
        let (step, y, f, dr2) = loop {
            let step = map::evaluate(
                &mut workspace,
                &mut states,
                &mut kept,
                comm,
                input.cell,
                opts.functional,
                &n,
                &policy,
                ethr,
                iteration == 1 && !redo,
                diagnostic,
            )?;
            collect_calls(&step, &mut calls);
            let y = transform(&mut plan, comm, &step.density, &policy, electrons, volume)?;
            let f: Vec<_> = y.iter().zip(&x).map(|(a, b)| a - b).collect();
            let base_dr2 = metric.norm(comm, &f)?;
            let stiffness_ha = if stiffness {
                workspace.spin_stiffness(comm, &n, &step.density, dv)?
            } else {
                0.
            };
            let dr2 = if stiffness {
                base_dr2 + 2. * stiffness_ha
            } else {
                base_dr2
            };
            if let Some(trace) = trace.as_deref_mut() {
                trace.record(
                    comm,
                    &mut metric,
                    &plan.g2,
                    &f,
                    &step,
                    iteration,
                    false,
                    redo,
                )?;
                trace.records[trace.used - 1].stiffness_ha = stiffness_ha;
            }
            if iteration == 1 && !redo && dr2 < ethr * electrons.max(1.) {
                ethr = next_ethr_ry(ethr, dr2, electrons);
                redo = true;
                warm(&mut states, &step);
                for cache in &mut kept {
                    cache.invalidate();
                }
                continue;
            }
            break (step, y, f, dr2);
        };
        let ordinary_trace = trace.as_ref().map(|t| t.used - 1);
        previous = dr2;
        let valid = valid(&step, &policy, electrons);
        let stopped = dr2 < profile.threshold_ry
            && (!guard || step.guard < profile.threshold_ry / 2.)
            && valid;
        if stopped {
            check_gaps(&step, &policy, opts.eigen_tolerance)?;
        }
        let mut passed = stopped;
        if stopped && let Some(report) = &mut verification {
            let start = Instant::now();
            warm(&mut states, &step);
            let warm_check = config.on(Id::REFSCF_VERIFY_WARM);
            for cache in &mut kept {
                if !warm_check {
                    cache.invalidate();
                }
                cache.preserve_input(true);
            }
            let tolerance = next_ethr_ry(ethr, dr2, electrons);
            let fresh = map::evaluate(
                &mut workspace,
                &mut states,
                &mut kept,
                comm,
                input.cell,
                opts.functional,
                &step.density,
                &policy,
                tolerance,
                false,
                diagnostic,
            )?;
            for cache in &mut kept {
                cache.preserve_input(false);
            }
            collect_calls(&fresh, &mut verification_calls);
            check_gaps(&fresh, &policy, opts.eigen_tolerance)?;
            let z = transform(&mut plan, comm, &fresh.density, &policy, electrons, volume)?;
            let residual: Vec<_> = z.iter().zip(&y).map(|(a, b)| a - b).collect();
            if let Some(trace) = trace.as_deref_mut() {
                trace.record(
                    comm,
                    &mut metric,
                    &plan.g2,
                    &residual,
                    &fresh,
                    iteration,
                    true,
                    false,
                )?;
            }
            let base_error = metric.norm(comm, &residual)? / 2.;
            let stiffness_ha = if stiffness {
                workspace.spin_stiffness(comm, &step.density, &fresh.density, dv)?
            } else {
                0.
            };
            let error = if stiffness {
                base_error + stiffness_ha
            } else {
                base_error
            };
            if let Some(trace) = trace.as_deref_mut() {
                trace.records[trace.used - 1].stiffness_ha = stiffness_ha;
            }
            let distance = distances(comm, &step, &fresh, dv, electrons)?;
            let mut eig_change = 0_f64;
            let mut occ_change = 0.;
            let mut eigen_residual = 0_f64;
            for (s, (a, b)) in step.blocks.iter().zip(&fresh.blocks).enumerate() {
                for (k, (a, b)) in a.iter().zip(b).enumerate() {
                    for i in 0..bands {
                        occ_change +=
                            points[k].weight * (a.occupations[i] - b.occupations[i]).abs();
                        if policy.finite().is_some() || i < policy.fixed_count(s) {
                            eig_change =
                                eig_change.max((a.eigenvalues[i] - b.eigenvalues[i]).abs());
                            eigen_residual = eigen_residual.max(b.residuals[i]);
                        }
                    }
                }
            }
            let mu_change = step.evaluation.zip(fresh.evaluation).map(|(a, b)| {
                a.roots
                    .iter()
                    .zip(b.roots)
                    .filter_map(|(a, b)| a.zip(b))
                    .map(|(a, b)| (a.mu_ha - b.mu_ha).abs())
                    .fold(0., f64::max)
            });
            passed = error < profile.threshold_ry / 2.
                && (!guard || fresh.guard < profile.threshold_ry / 2.)
                && valid_map(&fresh, &policy, electrons);
            let all_calls = fresh.blocks.iter().flatten().flat_map(|b| &b.calls);
            report.attempts.push(VerificationAttempt {
                iteration,
                passed,
                coulomb_error_ha: error,
                threshold_ha: profile.threshold_ry / 2.,
                eigen_tolerance_ha: tolerance / 2.,
                residual_bound_ha: guard.then_some(fresh.guard),
                density_residual: distance.0,
                spin_density_residual: (policy.channels() == 2).then_some(distance.1),
                occupied_eigenvalue_change_ha: eig_change,
                strict_criteria_met: distance.0 <= opts.density_tolerance
                    && distance.1 <= opts.density_tolerance
                    && eigen_residual <= opts.eigen_tolerance
                    && (step.free_energy() - fresh.free_energy()).abs() <= opts.energy_tolerance
                    && if policy.finite().is_some() {
                        occ_change <= 1e-8 * electrons.max(1.)
                    } else {
                        eig_change <= 1e-7
                    },
                davidson_iterations: all_calls.clone().map(|c| c.iterations).sum(),
                operator_applications: all_calls.map(|c| c.h_columns).sum(),
                wall: start.elapsed(),
                occupation_change_electrons: policy.finite().map(|_| occ_change),
                chemical_potential_change_ha: mu_change,
                moment_change: (policy.channels() == 2)
                    .then_some((step.moment - fresh.moment).abs()),
            });
        }
        let mut mixing = Diagnostics::default();
        let mut history_conditioning = ConditioningDiagnostics::default();
        let mut mixing_fraction = 1.;
        if !passed && iteration < opts.max_iterations {
            let mut next = x.clone();
            let update = mixer.update(comm, &mut metric, &x, &f, &mut next);
            if let Some(trace) = trace.as_deref_mut() {
                let condition = match &update {
                    Ok(d) => d.condition,
                    Err(super::Error::Singular { condition, .. }) => Some(*condition),
                    _ => None,
                };
                let record = &mut trace.records[ordinary_trace.unwrap()];
                record.condition = condition;
                let d = mixer.conditioning_diagnostics();
                record.dropped_pairs = d.dropped_pairs;
                record.dropped_conditions = d.dropped_conditions;
                (record.history_pair_norm2_ry, record.history_max_cosine) =
                    mixer.history_geometry();
            }
            mixing = update?;
            history_conditioning = mixer.conditioning_diagnostics();
            let mut proposal = inverse(&mut plan, comm, &next, policy.channels())?;
            if policy.channels() == 2 {
                for (s, channel) in proposal.iter_mut().enumerate() {
                    if policy.finite().is_none() && policy.fixed_count(s) == 0 {
                        channel.fill(0.);
                    }
                }
                // Same common feasibility fraction as the legacy spin driver;
                // preserve both number constraints and never clip a channel.
                let mut fraction = 1_f64;
                for (old, new) in n.iter().zip(&proposal) {
                    for (&a, &b) in old.iter().zip(new) {
                        if b < 0. {
                            fraction = fraction.min(0.9 * a / (a - b));
                        }
                    }
                }
                let mut fractions = vec![0.; comm.size()];
                comm.allgatherv_f64(&[fraction], &vec![1; comm.size()], &mut fractions)
                    .map_err(super::numerical)?;
                mixing_fraction = fractions.into_iter().fold(1_f64, f64::min);
                if !mixing_fraction.is_finite() || mixing_fraction <= 0. {
                    return Err(ScfError::InvalidInput(
                        "spin mixing has no feasible common step",
                    ));
                }
                if mixing_fraction < 1. {
                    for (old, new) in n.iter().zip(&mut proposal) {
                        for (&a, b) in old.iter().zip(new) {
                            *b = a + mixing_fraction * (*b - a);
                        }
                    }
                }
                x = transform(&mut plan, comm, &proposal, &policy, electrons, volume)?;
            } else {
                x = next;
            }
            n = proposal;
            warm(&mut states, &step);
        }
        history.push(Iteration {
            iteration,
            dr2_ry: dr2,
            ethr_ry: ethr,
            first_iteration_redo: redo,
            energy_ha: step.energy.total_hartree(),
            entropy_term_ha: step.entropy(),
            free_energy_ha: step.free_energy(),
            moment_mu_b: step.moment,
            residual_guard_ha: step.guard,
            mixing,
            history_conditioning,
            mixing_fraction,
            calls,
        });
        if passed || iteration == opts.max_iterations {
            let mut densities = Vec::new();
            let mut blocks = Vec::new();
            for (s, channel) in step.blocks.iter().enumerate() {
                densities.push(states[s].fields.gather(comm, &step.density[s])?);
                let mut delivered = Vec::new();
                for (k, b) in channel.iter().enumerate() {
                    delivered.push(BandBlock {
                        point: points[k],
                        eigenvalues_ha: b.eigenvalues.clone(),
                        occupations: b.occupations.clone(),
                        residuals_ha: b.residuals.clone(),
                        orbitals: super::k_driver::gather_orbitals(
                            &states[s].blocks[k],
                            comm,
                            &b.orbitals,
                        )?,
                    });
                }
                blocks.push(delivered);
            }
            return Ok(Run {
                converged: passed,
                density: densities,
                blocks,
                thermodynamics: ThermodynamicLedger::new(
                    step.energy.total_hartree(),
                    step.entropy(),
                ),
                energy: step.energy,
                occupations: step.evaluation,
                channel_electrons: step.charges,
                moment_mu_b: step.moment,
                absolute_moment_mu_b: step.absolute_moment,
                history,
                verification,
                verification_calls,
                kept_hx: kept
                    .iter()
                    .map(|k| k.report(comm))
                    .collect::<Result<_, _>>()?,
                setup_s,
                total_s: clock.elapsed().as_secs_f64(),
            });
        }
    }
    unreachable!("positive budget checked")
}
fn collect_calls(step: &Map, out: &mut Vec<(usize, usize, SolveRecord)>) {
    for (s, c) in step.blocks.iter().enumerate() {
        for (k, b) in c.iter().enumerate() {
            out.extend(b.calls.iter().cloned().map(|call| (s, k, call)));
        }
    }
}
fn warm(states: &mut [refstart::KStartState], step: &Map) {
    for (state, channel) in states.iter_mut().zip(&step.blocks) {
        for (x, b) in state.orbitals.0.iter_mut().zip(channel) {
            x.clone_from(&b.orbitals);
        }
    }
}
fn valid(step: &Map, policy: &Policy, electrons: f64) -> bool {
    valid_map(step, policy, electrons)
}
fn valid_map(step: &Map, policy: &Policy, electrons: f64) -> bool {
    let tolerance = 1e-10 * electrons.max(1.);
    (step.charges.iter().sum::<f64>() - electrons).abs() <= tolerance
        && policy
            .target_moment()
            .is_none_or(|m| (step.moment - m).abs() <= tolerance)
        && step.evaluation.is_none_or(|e| {
            e.electron_residual.abs() <= tolerance
                && e.moment_residual.unwrap_or(0.).abs() <= tolerance
        })
}
fn check_gaps(step: &Map, policy: &Policy, eigen_tolerance: f64) -> Result<(), ScfError> {
    if policy.finite().is_some() {
        return Ok(());
    }
    for (s, channel) in step.blocks.iter().enumerate() {
        let n = policy.fixed_count(s);
        if n == 0 {
            continue;
        }
        let high = channel
            .iter()
            .map(|b| b.eigenvalues[n - 1])
            .fold(f64::NEG_INFINITY, f64::max);
        let low = channel
            .iter()
            .map(|b| b.eigenvalues[n])
            .fold(f64::INFINITY, f64::min);
        if low - high <= 1e-7_f64.max(2. * eigen_tolerance) {
            return Err(ScfError::NonInsulating {
                gap_hartree: low - high,
                required_hartree: 1e-7_f64.max(2. * eigen_tolerance),
            });
        }
    }
    Ok(())
}
fn transform<K: Collective>(
    plan: &mut DensityPlan,
    comm: &K,
    density: &[Vec<f64>],
    policy: &Policy,
    electrons: f64,
    volume: f64,
) -> Result<Vec<C>, ScfError> {
    let size = plan.miller.len();
    let mut values = vec![C::default(); size * density.len()];
    let total: Vec<_> = (0..density[0].len())
        .map(|i| density.iter().map(|n| n[i]).sum())
        .collect();
    plan.forward(comm, &total, &mut values[..size])?;
    if density.len() == 2 {
        let spin: Vec<_> = density[0]
            .iter()
            .zip(&density[1])
            .map(|(u, d)| u - d)
            .collect();
        plan.forward(comm, &spin, &mut values[size..])?;
    }
    for (i, m) in plan.miller.iter().enumerate() {
        if *m == [0; 3] {
            values[i] = C::new(electrons / volume, 0.);
            if let Some(moment) = policy.target_moment() {
                values[size + i] = C::new(moment / volume, 0.);
            }
        }
    }
    Ok(values)
}
fn inverse<K: Collective>(
    plan: &mut DensityPlan,
    comm: &K,
    x: &[C],
    ns: usize,
) -> Result<Vec<Vec<f64>>, ScfError> {
    let size = plan.miller.len();
    let mut total = vec![0.; plan.plane_len()];
    plan.inverse(comm, &x[..size], &mut total)?;
    if ns == 1 {
        return Ok(vec![total]);
    }
    let mut spin = vec![0.; total.len()];
    plan.inverse(comm, &x[size..], &mut spin)?;
    Ok(vec![
        total.iter().zip(&spin).map(|(n, m)| (n + m) / 2.).collect(),
        total.iter().zip(&spin).map(|(n, m)| (n - m) / 2.).collect(),
    ])
}
fn distances<K: Collective>(
    comm: &K,
    a: &Map,
    b: &Map,
    dv: f64,
    electrons: f64,
) -> Result<(f64, f64), ScfError> {
    let mut local = [0.; 2];
    for i in 0..a.density[0].len() {
        let up = a.density[0][i] - b.density[0][i];
        let down = if a.density.len() == 2 {
            a.density[1][i] - b.density[1][i]
        } else {
            0.
        };
        local[0] += (up + down).powi(2);
        if a.density.len() == 2 {
            local[1] += (up - down).powi(2);
        }
    }
    let sums = map::sum(comm, &local)?;
    Ok((
        (dv * sums[0]).sqrt() / electrons,
        (dv * sums[1]).sqrt() / electrons,
    ))
}
