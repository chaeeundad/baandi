//! Reference SCF for scalar FD and collinear fixed/FD occupations.
//! The fixed unpolarized entry points never call this module.
use super::stepper::{Boundary, StepKind, StepReport, Stepper};
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
    Functional,
    energy::EnergyAudit,
    model::{
        ScfError, SpinSector, ThermodynamicLedger, Verification, VerificationAttempt,
        VerificationReport,
    },
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
    /// Stall events from this invocation (also emitted at step boundaries).
    pub stalls: Vec<super::stall::Diagnostic>,
    pub start_report: crate::refstart::StartReport,
    pub symmetry: Option<super::symmetry::Report>,
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
    /// Diagnostic only: complete PBE secant and nonlinear energy remainder.
    /// Zero unless tracing free-moment PBE; neither changes the stop rule.
    pub xc_response_ha: f64,
    pub xc_remainder_ha: f64,
    /// Frozen-spectrum FD number derivative, including rounded-away holes.
    pub fd_number_slope_e_per_ha: f64,
    /// Number of combined valence plus half-core channel samples changing mask.
    pub xc_mask_crossings: usize,
    /// Input mask compared with the preceding ordinary input (not output/verification).
    pub input_screened_points: usize,
    pub input_mask_changes: usize,
    /// Absolute screened input charge by channel, including any immutable core.
    pub input_screened_charge_electrons: [f64; 2],
    pub mixing_fraction: f64,
    /// Sampled indirect gap and maximum edge Ritz residual sum (Ha), fixed sectors only.
    pub channel_gap_ha: [Option<f64>; 2],
    pub channel_edge_residual_ha: [Option<f64>; 2],
    /// max_k (N_occ - ||C_previous^* C_current||_F^2) / N_occ.
    /// Invariant under phases and rotations within the occupied subspace.
    pub occupied_subspace_loss: [Option<f64>; 2],
    /// Fraction of squared residual norm in the explicitly requested rotation span.
    pub rotation_fraction: Option<f64>,
    pub rotation_rank: usize,
    pub dropped_pairs: usize,
    pub dropped_conditions: [f64; 25],
}
#[derive(Default)]
pub struct Trace {
    records: Vec<MapTrace>,
    used: usize,
    previous_mask: Vec<bool>,
    mask_seen: bool,
    mask_counts: Vec<f64>,
    mask_partition: Vec<usize>,
    previous_orbitals: Vec<Vec<ComplexMatrix>>,
    orbitals_seen: bool,
    overlap_local: Vec<f64>,
    overlap_global: Vec<f64>,
    rotation_request: Option<(Vec<[usize; 4]>, f64)>,
    rotation_probe: Option<super::rotation::Probe>,
}
impl Trace {
    /// Select (channel, k point, occupied band, empty band) pairs, zero based.
    /// Measures both phase quadratures; it does not assert these are zero modes.
    pub fn with_rotation_probe(pairs: Vec<[usize; 4]>, angle_radians: f64) -> Self {
        Self {
            rotation_request: Some((pairs, angle_radians)),
            ..Self::default()
        }
    }
    pub fn records(&self) -> &[MapTrace] {
        &self.records[..self.used]
    }
    fn prepare(&mut self, maps: usize, occupations: usize, samples: usize, ranks: usize) {
        self.previous_mask = vec![false; samples];
        self.mask_seen = false;
        self.mask_counts = vec![0.; 4 * ranks];
        self.mask_partition = vec![4; ranks];
        self.used = 0;
        self.records = (0..maps)
            .map(|_| MapTrace {
                occupations: vec![0.; occupations],
                ..Default::default()
            })
            .collect();
    }
    fn prepare_orbitals(&mut self, states: &[crate::refstart::KStartState], policy: &Policy) {
        if policy.finite().is_none() {
            self.previous_orbitals = states.iter().map(|s| s.orbitals.0.clone()).collect();
            let n = (0..policy.channels())
                .map(|s| policy.fixed_count(s))
                .max()
                .unwrap_or(0);
            self.overlap_local = vec![0.; 2 * n * n];
            self.overlap_global = vec![0.; 2 * n * n];
        }
        self.orbitals_seen = false;
    }
    fn spectrum<K: Collective>(
        &mut self,
        comm: &K,
        step: &Map,
        policy: &Policy,
    ) -> Result<(), ScfError> {
        if policy.finite().is_some() {
            return Ok(());
        }
        let record = &mut self.records[self.used - 1];
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
            record.channel_gap_ha[s] = Some(low - high);
            record.channel_edge_residual_ha[s] = Some(
                channel
                    .iter()
                    .map(|b| b.residuals[n - 1] + b.residuals[n])
                    .fold(0., f64::max),
            );
            if record.verification {
                continue;
            }
            let mut loss = 0_f64;
            for (k, b) in channel.iter().enumerate() {
                let previous = &mut self.previous_orbitals[s][k];
                if self.orbitals_seen {
                    loss = loss.max(subspace_loss(
                        comm,
                        previous,
                        &b.orbitals,
                        n,
                        &mut self.overlap_local,
                        &mut self.overlap_global,
                    )?);
                }
                previous
                    .as_mut_slice()
                    .copy_from_slice(b.orbitals.as_slice());
            }
            if self.orbitals_seen {
                record.occupied_subspace_loss[s] = Some(loss);
            }
        }
        if !record.verification {
            self.orbitals_seen = true;
        }
        Ok(())
    }
    fn input_mask<K: Collective>(
        &mut self,
        comm: &K,
        density: &[Vec<f64>],
        core: Option<&[f64]>,
        dv: f64,
    ) -> Result<(), ScfError> {
        if !dft_physics::xc::screen_negative() {
            return Ok(());
        }
        let mut counts = [0.; 4];
        for (j, &n) in density.iter().flatten().enumerate() {
            let c = core.map_or(0., |c| c[j % density[0].len()] / density.len() as f64);
            let screened = n + c <= dft_physics::xc::DENSITY_THRESHOLD;
            counts[0] += f64::from(screened);
            if screened {
                counts[2 + j / density[0].len()] += (n + c).abs() * dv;
            }
            counts[1] += f64::from(self.mask_seen && self.previous_mask[j] != screened);
            self.previous_mask[j] = screened;
        }
        comm.allgatherv_f64(&counts, &self.mask_partition, &mut self.mask_counts)
            .map_err(super::numerical)?;
        let record = &mut self.records[self.used - 1];
        record.input_screened_points = self
            .mask_counts
            .as_chunks::<4>()
            .0
            .iter()
            .map(|c| c[0] as usize)
            .sum();
        record.input_mask_changes = self
            .mask_counts
            .as_chunks::<4>()
            .0
            .iter()
            .map(|c| c[1] as usize)
            .sum();
        for s in 0..density.len() {
            record.input_screened_charge_electrons[s] = self
                .mask_counts
                .as_chunks::<4>()
                .0
                .iter()
                .map(|c| c[2 + s])
                .sum();
        }
        self.mask_seen = true;
        Ok(())
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
        record.xc_response_ha = step.xc_response_ha;
        record.xc_remainder_ha = step.xc_remainder_ha;
        record.fd_number_slope_e_per_ha = step.fd_number_slope_e_per_ha;
        record.xc_mask_crossings = step.xc_mask_crossings;
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
/// Fraction of occupied projector weight lost between two maps. Reducing the
/// complex overlaps before squaring is essential on distributed plane-wave rows.
fn subspace_loss<K: Collective>(
    comm: &K,
    previous: &ComplexMatrix,
    current: &ComplexMatrix,
    occupied: usize,
    local: &mut [f64],
    global: &mut [f64],
) -> Result<f64, ScfError> {
    for i in 0..occupied {
        for j in 0..occupied {
            let mut value = C::default();
            for row in 0..previous.rows() {
                value += previous[(row, i)].conj() * current[(row, j)];
            }
            local[2 * (i * occupied + j)] = value.re;
            local[2 * (i * occupied + j) + 1] = value.im;
        }
    }
    comm.tree_sum_f64(
        &local[..2 * occupied * occupied],
        &mut global[..2 * occupied * occupied],
    )
    .map_err(super::numerical)?;
    let retained: f64 = global[..2 * occupied * occupied]
        .iter()
        .map(|v| v * v)
        .sum();
    Ok(((occupied as f64 - retained) / occupied as f64).clamp(0., 1.))
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
    run_inner(input, profile, points, policy, comm, None, None)
}
/// Diagnostic entry; `None` retains the ordinary driver arithmetic and storage.
pub fn run_with_trace<K: Collective>(
    input: StartInput<'_>,
    profile: Profile,
    points: &[KPoint],
    policy: Policy,
    comm: &K,
    trace: Option<&mut Trace>,
) -> Result<Run, ScfError> {
    run_inner(input, profile, points, policy, comm, None, trace)
}
pub fn run_prepared(
    input: StartInput<'_>,
    profile: Profile,
    prepared: &super::symmetry::Prepared,
    policy: Policy,
) -> Result<Run, ScfError> {
    if let Some(comm) = crate::parallel::SharedComm::active() {
        run_prepared_with(input, profile, prepared, policy, &comm)
    } else {
        run_prepared_with(input, profile, prepared, policy, &Serial)
    }
}
pub fn run_prepared_with<K: Collective>(
    input: StartInput<'_>,
    profile: Profile,
    prepared: &super::symmetry::Prepared,
    policy: Policy,
    comm: &K,
) -> Result<Run, ScfError> {
    run_inner(
        input,
        profile,
        prepared.points(),
        policy,
        comm,
        Some(prepared),
        None,
    )
}
pub fn stepper_with<'a, K: Collective>(
    input: StartInput<'a>,
    profile: Profile,
    points: &'a [KPoint],
    policy: Policy,
    comm: &'a K,
) -> Stepper<'a, Run> {
    Stepper::new(move |boundary| async move {
        run_continuation(input, profile, points, policy, comm, None, None, &boundary).await
    })
}
fn run_inner<K: Collective>(
    input: StartInput<'_>,
    profile: Profile,
    points: &[KPoint],
    policy: Policy,
    comm: &K,
    prepared: Option<&super::symmetry::Prepared>,
    trace: Option<&mut Trace>,
) -> Result<Run, ScfError> {
    Stepper::new(move |boundary| async move {
        run_continuation(
            input, profile, points, policy, comm, prepared, trace, &boundary,
        )
        .await
    })
    .run_to_end()
}
#[allow(clippy::too_many_arguments)]
pub(super) async fn run_continuation<K: Collective>(
    input: StartInput<'_>,
    profile: Profile,
    points: &[KPoint],
    policy: Policy,
    comm: &K,
    prepared: Option<&super::symmetry::Prepared>,
    mut trace: Option<&mut Trace>,
    boundary: &Boundary,
) -> Result<Run, ScfError> {
    if policy.channels() == 2 && prepared.is_some_and(|p| p.time_reversal) {
        return Err(ScfError::InvalidInput(
            "symmetry: collinear preparation must disable time reversal",
        ));
    }
    let clock = Instant::now();
    let mut profile = profile.validate()?;
    let opts = input.options;
    let config = RunConfig::current();
    if opts.standard.is_none()
        || opts.max_iterations == 0
        || (opts.crystal_density_group && prepared.is_none())
    {
        return Err(ScfError::InvalidInput(
            "reference spin/FD requires standard profile, positive budget, prepared symmetry when requested",
        ));
    }
    if config.on(Id::REFSCF_RMM) || config.on(Id::REFSCF_LADDER) {
        return Err(ScfError::InvalidInput(
            "reference spin/FD supports Davidson and fixed-point verification only; RMM, ladder, energy verification and failed-map reuse unavailable",
        ));
    }
    if policy.channels() == 2 && opts.functional == Functional::LdaPz {
        return Err(ScfError::InvalidInput(
            "reference spin XC requires lda-pw or pbe: legacy polarized lda-pz unavailable",
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
    options.gamma_real = crate::model::GammaReal::Off;
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
        crate::model::spin::validate_options(&options, theta, bands)?;
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
    if let Some(p) = prepared {
        p.validate(&input, points)?;
        options.grid_override = Some(p.grid(&input)?);
    }
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
    let mut symmetry = prepared
        .map(|p| {
            super::symmetry::Workspace::new(
                &input,
                p,
                dims,
                &states[0].fields,
                policy.channels(),
                comm,
            )
        })
        .transpose()?;
    if let Some(s) = &mut symmetry {
        s.check_start(&mut states[0].fields, comm, &n)?;
    }
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
    let stiffness = spin_fd;
    if spin_fd {
        mixer.condition_history();
    }
    let mut kept: Vec<_> = (0..policy.channels()).map(|_| KeptHx::new(true)).collect();
    let mut x = transform(&mut plan, comm, &n, &policy, electrons, volume)?;
    let mut workspace = map::Workspace::new(
        &states[0],
        comm,
        input.cell,
        policy.channels() == 2,
        opts.functional,
    )?;
    if let Some(trace) = trace.as_deref_mut() {
        trace.prepare(
            2 * opts.max_iterations + 1,
            policy.channels() * points.len() * bands,
            policy.channels() * n[0].len(),
            comm.size(),
        );
    }
    if let Some(trace) = trace.as_deref_mut() {
        trace.prepare_orbitals(&states, &policy);
        trace.rotation_probe = trace
            .rotation_request
            .as_ref()
            .map(|(pairs, angle)| {
                super::rotation::Probe::new(pairs, *angle, &states, &plan, comm, volume)
            })
            .transpose()?;
    }
    workspace.symmetry = symmetry;
    let setup_s = clock.elapsed().as_secs_f64();
    let mut stalls = Vec::with_capacity(opts.max_iterations);
    let mut history = Vec::new();
    let mut verification =
        (opts.verification == Verification::FixedPoint).then(VerificationReport::default);
    let mut verification_calls = Vec::new();
    let mut ethr = 1e-2;
    let mut previous = f64::INFINITY;
    let guard = config.on(Id::STANDARD_RESIDUAL_GUARD);
    let diagnostic = config.on(Id::DAVIDSON_STAGES);
    let mut completed = 0usize;
    macro_rules! checkpoint_state {
        () => {
            |request: super::stepper::Visit<'_, '_>| -> Result<usize, ScfError> {
                let io = match request {
                    super::stepper::Visit::Checkpoint(io) => io,
                    super::stepper::Visit::Policy(patch) => {
                        let next = patch.prepare(profile)?;
                        if patch.beta.is_some() || patch.history.is_some() {
                            let mut replacement = if policy.channels() == 2 {
                                Broyden::spin(next, &plan.g2, &plan.smooth, volume, electrons)?
                            } else {
                                Broyden::new(next, &plan.g2, &plan.smooth, volume, electrons)?
                            };
                            if spin_fd {
                                replacement.condition_history();
                            }
                            mixer = replacement;
                        }
                        if let Some(enabled) = patch.verification {
                            if enabled != verification.is_some() {
                                verification = enabled.then(VerificationReport::default);
                            }
                        }

                        profile = next;

                        return Ok(completed);
                    }
                };
                if trace.is_some() {
                    return Err(super::checkpoint::invalid(
                        "electronic checkpoint with external trace is unsupported",
                    ));
                }
                super::checkpoint::header(io, &input, profile, comm, "electronic", dims)?;
                io.same(&format!("{policy:?}"))?;
                super::checkpoint::symmetry(io, &mut workspace.symmetry)?;
                io.same(&plan.g2)?;
                io.value(&mut completed)?;
                io.same(&n.len())?;
                for channel in &mut n {
                    io.planes(channel)?;
                }
                io.reciprocal(&mut x)?;
                mixer.checkpoint(io, &metric)?;
                for (state, cache) in states.iter_mut().zip(&mut kept) {
                    super::checkpoint::k_state(io, state)?;
                    cache.checkpoint(io, state)?;
                }
                io.value(&mut ethr)?;
                io.value(&mut previous)?;
                io.value(&mut history)?;
                io.value(&mut verification)?;
                io.value(&mut verification_calls)?;
                super::checkpoint::counters(
                    completed,
                    opts.max_iterations,
                    history.len(),
                    ethr,
                    previous,
                )?;
                Ok(completed)
            }
        };
    }
    boundary.ready_state(checkpoint_state!()).await;
    for iteration in completed + 1..=opts.max_iterations {
        if let Some(s) = &mut workspace.symmetry {
            s.iteration = iteration;
            s.verification = false;
        }
        if iteration > 1 {
            ethr = next_ethr_ry(
                if iteration == 2 { 1e-2 } else { ethr },
                previous,
                electrons,
            );
        }
        let mut calls = Vec::new();
        let mut redo = false;
        let mut reported_columns = 0;
        let (step, y, f, dr2) = loop {
            let xc_stage = dft_physics::xc::screening_scope("map", false, Some(electrons));
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
                trace.is_some(),
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
            // Restore thread-local diagnostics before any yielded boundary.
            drop(xc_stage);
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
                trace.input_mask(comm, &n, states[0].core_density.as_deref(), dv)?;
                trace.spectrum(comm, &step, &policy)?;
                if let Some(probe) = &mut trace.rotation_probe {
                    let (fraction, rank) =
                        probe.measure(&mut states, &step, &mut plan, &mut metric, comm, &f)?;
                    trace.records[trace.used - 1].rotation_fraction = Some(fraction);
                    trace.records[trace.used - 1].rotation_rank = rank;
                }
            }
            if iteration == 1 && !redo && dr2 < ethr * electrons.max(1.) {
                let report = StepReport::new(
                    iteration,
                    StepKind::Ordinary,
                    dr2,
                    step.energy.total_hartree(),
                    ethr,
                    calls.iter().map(|(_, _, c)| c.h_columns).sum(),
                );
                reported_columns = report.h_columns;
                boundary.emit(report).await;
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
            calls.iter().map(|(_, _, c)| c.h_columns).sum::<usize>() - reported_columns,
        );
        let mut passed = stopped;
        if stopped && let Some(report) = &mut verification {
            boundary.emit(map_report).await;
            let start = Instant::now();
            warm(&mut states, &step);

            for cache in &mut kept {
                cache.preserve_input(true);
            }
            let tolerance = next_ethr_ry(ethr, dr2, electrons);
            if let Some(s) = &mut workspace.symmetry {
                s.verification = true;
            }
            let _xc_stage =
                dft_physics::xc::screening_scope("verification", false, Some(electrons));
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
                trace.is_some(),
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
                trace.spectrum(comm, &fresh, &policy)?;
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
            map_report = StepReport::new(
                iteration,
                StepKind::Verification,
                error * 2.,
                fresh.energy.total_hartree(),
                tolerance,
                all_calls.clone().map(|c| c.h_columns).sum(),
            );
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
        if let Some(trace) = trace.as_deref_mut() {
            trace.records[ordinary_trace.unwrap()].mixing_fraction = mixing_fraction;
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
        if super::stall::detected(history.iter().rev().map(|r| (r.dr2_ry, r.mixing.condition))) {
            let mut gaps = [None; 2];
            for (s, channel) in step.blocks.iter().enumerate() {
                let n = policy.fixed_count(s);
                if n > 0 {
                    let high = channel
                        .iter()
                        .map(|b| b.eigenvalues[n - 1])
                        .fold(f64::NEG_INFINITY, f64::max);
                    let low = channel
                        .iter()
                        .map(|b| b.eigenvalues[n])
                        .fold(f64::INFINITY, f64::min);
                    gaps[s] = Some(low - high);
                }
            }
            let diagnostic = super::stall::Diagnostic {
                iteration,
                channel_gap_ha: gaps,
                orbital_rotation_fraction: trace
                    .as_ref()
                    .and_then(|t| t.records[ordinary_trace.unwrap()].rotation_fraction),
            };
            map_report.stall = Some(diagnostic);
            stalls.push(diagnostic);
        }
        map_report.safe_point = true;
        if passed || iteration == opts.max_iterations {
            boundary.publish(map_report);
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
                stalls,
                start_report: refstart::StartReport {
                    // Collinear channels use uniform counts, irrespective of the
                    // atomic density used to prepare their initial orbital span.
                    density: if states.len() == 2 {
                        refstart::DensityStartKind::Uniform
                    } else {
                        states[0].report.density
                    },
                    ..states[0].report
                },
                symmetry: workspace.symmetry.map(|s| s.report),
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
        completed = iteration;
        boundary.emit_state(map_report, checkpoint_state!()).await;
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

#[cfg(test)]
mod mask_trace_tests {
    use super::*;

    #[test]
    fn mask_trace_detects_location_changes_with_unchanged_counts_and_empty_owner() {
        let config = RunConfig::from_pairs(&[("xc_screen_negative", "on")]).unwrap();
        dft_mpi::simulate(3, |comm| {
            let _scope = config.enter();
            let t = dft_physics::xc::DENSITY_THRESHOLD;
            let mut trace = Trace::default();
            let density = if comm.rank() == 0 {
                vec![vec![0., 2. * t], vec![-t, 4. * t]]
            } else {
                vec![vec![], vec![]]
            };
            let core = if comm.rank() == 0 {
                vec![t, 0.]
            } else {
                vec![]
            };
            trace.prepare(2, 0, density.iter().map(Vec::len).sum(), comm.size());
            trace.used = 1;
            trace.input_mask(comm, &density, Some(&core), 2.).unwrap();
            assert_eq!(trace.records()[0].input_screened_points, 2);
            assert_eq!(trace.records()[0].input_mask_changes, 0);
            assert_eq!(trace.records()[0].input_screened_charge_electrons, [t, t]);
            let next = if comm.rank() == 0 {
                vec![vec![2. * t, 0.], vec![-t, 4. * t]]
            } else {
                vec![vec![], vec![]]
            };
            trace.used = 2;
            trace.input_mask(comm, &next, Some(&core), 2.).unwrap();
            assert_eq!(trace.records()[1].input_screened_points, 2);
            assert_eq!(trace.records()[1].input_mask_changes, 2);
            assert_eq!(trace.records()[1].input_screened_charge_electrons, [0., t]);
        });
    }
}

#[cfg(test)]
mod gap_tests {
    use super::*;
    #[test]
    fn occupied_projector_overlap_ignores_gauge_but_detects_empty_band_rotation() {
        for ranks in [1, 3, 8] {
            dft_mpi::simulate(ranks, |comm| {
                let rows: Vec<_> = (0..4).filter(|i| i % ranks == comm.rank()).collect();
                let mut previous = ComplexMatrix::zeros(rows.len(), 2).unwrap();
                let mut current = previous.clone();
                for (i, &row) in rows.iter().enumerate() {
                    for j in 0..2 {
                        previous[(i, j)] = C::new(f64::from(row == j), 0.);
                    }
                    // Swap the occupied basis and apply complex phases: same projector.
                    current[(i, 0)] = C::new(0., f64::from(row == 1));
                    current[(i, 1)] = C::new(-f64::from(row == 0), 0.);
                }
                let (mut local, mut global) = (vec![0.; 8], vec![0.; 8]);
                assert_eq!(
                    subspace_loss(comm, &previous, &current, 2, &mut local, &mut global).unwrap(),
                    0.
                );
                for (i, &row) in rows.iter().enumerate() {
                    current[(i, 0)] = C::new(f64::from(row == 0), 0.);
                    current[(i, 1)] = C::new(f64::from(row == 2), 0.);
                }
                assert_eq!(
                    subspace_loss(comm, &previous, &current, 2, &mut local, &mut global).unwrap(),
                    0.5
                );
            });
        }
    }
}
