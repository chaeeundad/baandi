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
    /// Fermi-Dirac occupations (the variational quantity is F = E - TS).
    pub fn is_finite(&self) -> bool {
        self.finite().is_some()
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
/// Automatic Fermi-Dirac band count (reference driver). The policy's
/// `band_cap` is the ceiling; the FD coverage criterion is unchanged.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BandGrowth {
    /// Start at `max(ceil(n) + 4, ceil(1.2 n))`, n the largest channel
    /// electron target (N/2 unpolarized; max(N/2, n_up, n_down) at free
    /// moment; (N+|m|)/2 at fixed moment), bounded by the cap.
    Automatic,
    /// Start at an explicit count (tests and studies), bounded by the cap.
    From(usize),
}
impl BandGrowth {
    /// The initial count of [`BandGrowth::Automatic`] for a channel target.
    pub fn automatic_initial(channel_electrons: f64) -> usize {
        let n = channel_electrons.max(0.);
        (n.ceil() as usize + 4).max((1.2 * n).ceil() as usize)
    }
    /// GB1 (2026-11-27): the band cap implied for a channel target `n` when
    /// the caller gives none: the automatic start followed by two growth
    /// steps, `next(next(automatic_initial(n)))`. It always exceeds `n`, so a
    /// finite chemical potential exists at the retained-capacity endpoint for
    /// every electron count, odd or even (the CLI applies a floor of 32).
    pub fn automatic_cap(channel_electrons: f64) -> usize {
        let start = Self::automatic_initial(channel_electrons);
        Self::next(Self::next(start, usize::MAX), usize::MAX)
    }
    /// One growth step: `max(bands + 4, ceil(1.2 bands))`, bounded by `limit`.
    pub fn next(bands: usize, limit: usize) -> usize {
        (bands + 4)
            .max((1.2 * bands as f64).ceil() as usize)
            .min(limit)
    }
}
/// One automatic band growth: the map (ordinary or verification) of
/// `iteration` failed the coverage criterion at `from` bands and was solved
/// again from the same input density at `to` bands.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GrowthEvent {
    pub iteration: usize,
    pub verification: bool,
    pub from: usize,
    pub to: usize,
    pub edge_max_fraction: f64,
    pub edge_electrons: f64,
}
/// GB1: a map of a fixed-count FD run that missed the empty-state coverage
/// criterion before delivery (diagnostic only; the delivered state, its
/// verification map and the last map of an exhausted budget are checked hard).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoverageMiss {
    pub iteration: usize,
    /// The first-iteration redo map (same input density, tighter ethr).
    pub redo: bool,
    pub edge_max_fraction: f64,
    pub edge_electrons: f64,
}
#[derive(Debug)]
pub struct BandBlock {
    pub point: KPoint,
    pub eigenvalues_ha: Vec<f64>,
    pub occupations: Vec<f64>,
    pub residuals_ha: Vec<f64>,
    /// 0 x bands placeholder while `lazy` holds the block (GK18B; see
    /// [`Run::materialize_orbitals`]).
    pub orbitals: ComplexMatrix,
    pub(crate) lazy: Option<super::lazy::LazyOrbitals>,
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
    pub map_diagnostics: Vec<super::stall::MapDiagnostic>,
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
    /// Automatic band count: initial count, cap and growth events (`None`
    /// when the count was fixed by the caller).
    pub band_growth: Option<BandGrowthReport>,
    /// Fixed-count FD: transient maps that missed the coverage criterion
    /// (empty with automatic growth, which grows instead; since a restart
    /// only the maps solved by this invocation).
    pub coverage_misses: Vec<CoverageMiss>,
    pub setup_s: f64,
    pub total_s: f64,
    /// Occupation policy of this delivery (derivatives read theta, cap, sector).
    pub policy: Policy,
    /// Setup-time model record checked by post-SCF derivatives.
    pub(crate) force_identity: crate::force_identity::ForceIdentity,
}
impl Run {
    /// GK18B: whether every k block is held as a full complex block (false
    /// only after `lazy::set_delivered_orbitals(false)` on a pooled run).
    pub fn orbitals_delivered(&self) -> bool {
        self.blocks.iter().flatten().all(|b| b.lazy.is_none())
    }
    /// Gather lazily delivered blocks (collective: every rank calls it at the
    /// same point; channels then k in order); a no-op otherwise.
    pub fn materialize_orbitals(&mut self) -> Result<(), ScfError> {
        for b in self.blocks.iter_mut().flatten() {
            if let Some(lazy) = b.lazy.take() {
                b.orbitals = lazy.into_full()?;
            }
        }
        Ok(())
    }
    /// Derivatives evaluate the Mermin free energy F = E - TS for finite
    /// temperature and the internal energy for a fixed spin sector.
    pub fn quantity(&self) -> &'static str {
        if self.policy.is_finite() { "F" } else { "E" }
    }
    pub(crate) fn occupation_identity(&self) -> crate::force_identity::OccupationIdentity {
        occupation_identity(&self.policy)
    }
}
fn occupation_identity(policy: &Policy) -> crate::force_identity::OccupationIdentity {
    let (finite, sector) = match policy {
        Policy::FermiDirac { theta_ha, band_cap } => {
            (Some((*theta_ha, Constraint::Total, *band_cap)), None)
        }
        Policy::FixedSpin(sector) => (None, Some(sector)),
        Policy::SpinFermiDirac {
            sector,
            theta_ha,
            constraint,
            band_cap,
        } => (Some((*theta_ha, *constraint, *band_cap)), Some(sector)),
    };
    crate::force_identity::OccupationIdentity::reference(finite, sector)
}
#[derive(Clone, Debug, PartialEq)]
pub struct BandGrowthReport {
    pub initial_bands: usize,
    pub final_bands: usize,
    pub band_cap: usize,
    pub events: Vec<GrowthEvent>,
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
    /// Spin-major, then k/band order, Ha; allocated only for an external trace.
    pub eigenvalues_ha: Vec<f64>,
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
                eigenvalues_ha: vec![0.; occupations],
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
        for (a, b) in record
            .eigenvalues_ha
            .iter_mut()
            .zip(step.blocks.iter().flatten().flat_map(|b| &b.eigenvalues))
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
        run_continuation(
            input, profile, points, policy, None, comm, None, None, &boundary,
        )
        .await
    })
}
/// [`stepper_with`] with automatic FD band growth (as `api::Request::band_growth`).
pub fn stepper_growing<'a, K: Collective>(
    input: StartInput<'a>,
    profile: Profile,
    points: &'a [KPoint],
    policy: Policy,
    growth: BandGrowth,
    comm: &'a K,
) -> Stepper<'a, Run> {
    Stepper::new(move |boundary| async move {
        run_continuation(
            input,
            profile,
            points,
            policy,
            Some(growth),
            comm,
            None,
            None,
            &boundary,
        )
        .await
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
            input, profile, points, policy, None, comm, prepared, trace, &boundary,
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
    growth: Option<BandGrowth>,
    comm: &K,
    prepared: Option<&super::symmetry::Prepared>,
    mut trace: Option<&mut Trace>,
    boundary: &Boundary,
) -> Result<Run, ScfError> {
    let config = RunConfig::current();
    // GT1 (ADR 0050): per-channel time reversal is valid for collinear
    // non-SOC models; off restores the former rejection.
    if policy.channels() == 2
        && prepared.is_some_and(|p| p.time_reversal)
        && !config.on(Id::COLLINEAR_TIME_REVERSAL)
    {
        return Err(ScfError::InvalidInput(
            "symmetry: collinear preparation must disable time reversal",
        ));
    }
    let clock = Instant::now();
    let mut profile = profile.validate()?;
    let opts = input.options;
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
    let cap = match &policy {
        Policy::FermiDirac { band_cap, .. } | Policy::SpinFermiDirac { band_cap, .. } => *band_cap,
        Policy::FixedSpin(sector) => options.bands.unwrap_or(
            sector
                .n_up
                .max(sector.n_down)
                .checked_add(options.extra_bands)
                .ok_or(ScfError::InvalidInput("band overflow"))?,
        ),
    };
    if growth.is_some() && (policy.finite().is_none() || options.bands.is_some() || trace.is_some())
    {
        return Err(ScfError::InvalidInput(
            "occupations.bands: automatic band growth requires Fermi-Dirac occupations, no explicit band count and no trace",
        ));
    }
    let mut bands = match growth {
        None => cap,
        Some(BandGrowth::From(n)) => n.min(cap),
        Some(BandGrowth::Automatic) => {
            let channel = match &policy {
                Policy::SpinFermiDirac {
                    constraint: Constraint::FixedMoment(m),
                    ..
                } => (electrons + m.abs()) / 2.,
                Policy::SpinFermiDirac { sector, .. } => {
                    (electrons / 2.).max(sector.n_up.max(sector.n_down) as f64)
                }
                _ => electrons / 2.,
            };
            BandGrowth::automatic_initial(channel).min(cap)
        }
    };
    let initial_bands = bands;
    let mut growth_events: Vec<GrowthEvent> = Vec::new();
    let mut coverage_misses: Vec<CoverageMiss> = Vec::new();
    if let Some((theta, constraint)) = policy.finite() {
        crate::model::spin::validate_options(&options, theta, cap)?;
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
    let mut states = refstart::start_electronic_channels(
        StartInput {
            options: &options,
            ..input
        },
        points,
        comm,
        policy.channels(),
        // GK9 pools: the external trace (subspace loss, rotation probe)
        // reads every k block's rows and keeps the distributed layout.
        trace.is_none(),
    )?;
    // FD buffer bands: the solved block holds `bands` plus the buffer; the
    // start's columns are kept bit for bit and the buffer is appended.
    let width = map::solve_width(&policy, bands, min_basis(&states[0]));
    if width > bands {
        for state in states.iter_mut() {
            refstart::extend_k_orbitals(state, width, comm)?;
        }
    }
    let volume = input.cell.volume_bohr3();
    let dims = states[0].grid.dims();
    let dv = volume / dims.iter().product::<usize>() as f64;
    // Legacy collinear starts from uniform channel counts. A positive seed
    // permits its common feasibility step without clipping atomic-density zeros.
    // GB2 (`SPIN_FD_ATOMIC_DENSITY`, collinear FD only, default off): when
    // the start has an atomic density, each channel starts from it scaled to
    // its count (N +- m) / 2, the spin analogue of the unpolarized atomic
    // start (zeros of the clipped atomic density are kept).
    let atomic_channels = matches!(policy, Policy::SpinFermiDirac { .. })
        && states[0].report.density == refstart::DensityStartKind::Atomic
        && config.on(Id::SPIN_FD_ATOMIC_DENSITY);
    let mut n: Vec<Vec<f64>> = states
        .iter()
        .enumerate()
        .map(|(s, state)| {
            if policy.channels() == 1 {
                state.density.clone()
            } else if atomic_channels {
                let scale = (electrons + if s == 0 { seed_moment } else { -seed_moment })
                    / (2. * electrons);
                state.density.iter().map(|v| v * scale).collect()
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
        states[0].potential_symmetry = s.take_potential_projector();
    }
    let mut plan = DensityPlan::with_fft(
        comm,
        input.cell,
        opts.strain_reference.as_deref(),
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
    let mut pulay = if dft_numerics::switches::on(dft_numerics::switches::Id::REFSCF_SPIN_PULAY) {
        if comm.size() != 1 || policy.finite().is_some() {
            return Err(ScfError::InvalidInput("serial fixed research only"));
        }
        Some(
            (0..policy.channels())
                .map(|_| {
                    dft_solvers::mixing::PulayMixer::new(profile.beta, profile.history).unwrap()
                })
                .collect::<Vec<_>>(),
        )
    } else {
        None
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
    let mut map_diagnostics = Vec::with_capacity(opts.max_iterations);
    let mut history = Vec::new();
    let mut verification =
        (opts.verification == Verification::FixedPoint).then(VerificationReport::default);
    let mut verification_calls = Vec::new();
    let mut ethr = 1e-2;
    let mut previous = f64::INFINITY;
    let guard = config.on(Id::STANDARD_RESIDUAL_GUARD);
    let diagnostic = config.on(Id::DAVIDSON_STAGES);
    let mut completed = 0usize;
    // GF2 (`FD_EMPTY_THRESHOLD`): occupations of the last SCF map per
    // (channel, global k), the empty-band rule's input for the next map.
    let empty_rule = policy.finite().is_some() && config.choice(Id::FD_EMPTY_THRESHOLD) != "off";
    let mut previous_occupations: Option<Vec<Vec<Vec<f64>>>> = None;
    let empty_verify = empty_rule && config.on(Id::FD_EMPTY_VERIFY);
    // GB2 (`FD_VERIFY_PREDICTOR`): the closed-shell verification predictor.
    // After an attempt with error e at residual r0, a later stop at r skips
    // its attempt (one ordinary map instead) when e r / r0 >= the threshold,
    // never twice in a row. The pass criterion is unchanged.
    let use_predictor = config.on(Id::FD_VERIFY_PREDICTOR);
    // GM1 (`SCF_ETHR_GUARD`): map accuracy relative to the SCF residual.
    let ethr_guard = config.choice(Id::SCF_ETHR_GUARD);
    let redo_every_map = matches!(ethr_guard, "redo" | "both");
    let forecast_ethr = matches!(ethr_guard, "forecast" | "both");
    let mut predictor = super::driver::VerificationPredictor::default();
    // GV1 (`VERIFY_FEED`): a failed attempt's pair (y, F(y)) may enter the
    // Broyden history (reference Broyden only; the Pulay variant ignores it).
    let feed = super::driver::VerifyFeed::resolve(config);
    let skip_run = super::driver::verify_skip_run(config);
    let response_predictor = super::driver::verify_response_predictor(config);
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
                // Automatic band count only: the fixed-count format is unchanged.
                if growth.is_some() {
                    io.same(&format!("{growth:?}"))?;
                    let mut count = bands;
                    io.value(&mut count)?;
                    if count > cap || count < bands {
                        return Err(super::checkpoint::invalid("band growth count"));
                    }
                    if count != bands {
                        // Shape only: the orbital block is read next. The
                        // solved block includes the FD buffer.
                        let width = map::fd_width(count, min_basis(&states[0]));
                        for state in states.iter_mut() {
                            for (block, x) in state.blocks.iter_mut().zip(&mut state.orbitals.0) {
                                *x = dft_numerics::linalg::ComplexMatrix::zeros(x.rows(), width)
                                    .map_err(super::numerical)?;
                                block
                                    .operator
                                    .ensure_capacity((2 * width).min(block.basis.vectors().len()))
                                    .map_err(super::numerical)?;
                            }
                        }
                        bands = count;
                    }
                    let mut events: Vec<_> = growth_events
                        .iter()
                        .map(|e| {
                            (
                                (e.iteration, e.verification),
                                (e.from, e.to),
                                (e.edge_max_fraction, e.edge_electrons),
                            )
                        })
                        .collect();
                    io.value(&mut events)?;
                    growth_events = events
                        .into_iter()
                        .map(
                            |((iteration, verification), (from, to), (f, e))| GrowthEvent {
                                iteration,
                                verification,
                                from,
                                to,
                                edge_max_fraction: f,
                                edge_electrons: e,
                            },
                        )
                        .collect();
                }
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
                // GF2: only with the rule on (the switch is in the header's
                // agreement hash); with `FD_EMPTY_THRESHOLD=off` and for
                // fixed-sector policies the format is the former one.
                if empty_rule {
                    io.value(&mut previous_occupations)?;
                }
                // GB2: only with the predictor on (switch in the header hash).
                if use_predictor {
                    io.value(&mut predictor.observed)?;
                    io.value(&mut predictor.skipped_last)?;
                    if skip_run {
                        io.value(&mut predictor.last_skip_dr2)?;
                    }
                }
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
    // GF1: seed every channel's kept HX with the
    // start's H X of the start Ritz orbitals (incl. the FD buffer columns);
    // `prepare` hits only while the orbitals are still those exact vectors
    // and corrects the potential change exactly. A restored run drops it.
    for (s, state) in states.iter_mut().enumerate() {
        let start_hx = std::mem::take(&mut state.start_hx);
        if completed == 0 {
            for (k, (x, hx, potential)) in start_hx.into_iter().enumerate() {
                kept[s].retain(k, state, &x, &hx, &potential);
            }
        }
    }
    // GB2 (`FD_GROWTH_AT_STOP`, default on): automatic growth reacts to a
    // coverage miss only at a stopping, budget-final or verification map;
    // off restores the pre-GB2 growth on every map.
    let legacy_growth = !config.on(Id::FD_GROWTH_AT_STOP);
    let mut transient_shortfall = map::Shortfall::default();

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
        let mut verification_grew = false;
        let (mut step, y, f, dr2) = loop {
            let xc_stage = dft_physics::xc::screening_scope("map", false, Some(electrons));
            // GB2: the ordinary map is never rejected for coverage; a miss is
            // judged below once its residual is known (growth only at a
            // stopping or budget-final map, a recorded transient otherwise).
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
                (growth.is_some() && legacy_growth).then_some(&mut transient_shortfall),
                bands,
                previous_occupations.as_deref(),
            );
            let step = match step {
                Err(ScfError::InvalidInput(m)) if growth.is_some() && m == map::COVERAGE => {
                    // GB2_TRANSIENT_GROWTH=off: the pre-GB2 rule (every map).
                    drop(xc_stage);
                    grow(
                        &mut states,
                        &mut kept,
                        comm,
                        &mut bands,
                        cap,
                        (iteration, false),
                        transient_shortfall,
                        &mut growth_events,
                    )?;
                    continue;
                }
                other => other?,
            };
            collect_calls(&step, &mut calls);
            if empty_rule {
                previous_occupations = Some(occupations_of(&step));
            }
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
            if let Some(miss) = step.coverage {
                let stopping = dr2 < profile.threshold_ry
                    && (!guard || step.guard < profile.threshold_ry / 2.)
                    && valid(&step, &policy, electrons);
                if growth.is_some() && (stopping || iteration == opts.max_iterations) {
                    // GB2: a map that would be delivered misses the criterion:
                    // grow from its solved orbitals and solve the same input
                    // density again (mixing history untouched).
                    warm(&mut states, &step);
                    grow(
                        &mut states,
                        &mut kept,
                        comm,
                        &mut bands,
                        cap,
                        (iteration, false),
                        miss,
                        &mut growth_events,
                    )?;
                    continue;
                }
                // A transient map may miss the criterion (an early density far
                // from the fixed point); the delivered state is checked.
                if comm.rank() == 0 && coverage_misses.is_empty() {
                    eprintln!(
                        "fd coverage: iteration {iteration}: transient map misses the empty-state criterion at {bands} bands (edge fraction {:.3e}, edge electrons {:.3e}); checked again at delivery",
                        miss.edge_max_fraction, miss.edge_electrons,
                    );
                }
                coverage_misses.push(CoverageMiss {
                    iteration,
                    redo,
                    edge_max_fraction: miss.edge_max_fraction,
                    edge_electrons: miss.edge_electrons,
                });
            }
            // The first-map rule; with GM1 `SCF_ETHR_GUARD` redo on every
            // ordinary map: a residual below the map's own diagonalization
            // error estimate is not trusted for mixing.
            if (iteration == 1 || redo_every_map) && !redo && dr2 < ethr * electrons.max(1.) {
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
                // GF1: the redo solves the same input potential from the
                // coarse solve's X, whose HX is kept exactly (zero potential
                // change).
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
            delivered_coverage(&step, comm, bands, false)?;
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
        let skipped = stopped && use_predictor && verification.is_some() && {
            let prediction = if response_predictor && predictor.observed.is_none() {
                // GV1: no attempt yet; the secant output-response estimate.
                mixer
                    .predict_output_residual(comm, &mut metric, &x, &f)?
                    .map(|v| v / 2.)
            } else {
                predictor.predict(dr2)
            };
            predictor.skip_with(prediction, profile.threshold_ry / 2., dr2, skip_run)
        };
        let mut passed = stopped && !skipped;
        let mut fed = None;
        if passed && let Some(report) = &mut verification {
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
            let before = bands;
            let mut shortfall = map::Shortfall::default();
            let fresh = loop {
                match map::evaluate(
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
                    growth.is_some().then_some(&mut shortfall),
                    bands,
                    // GF2: verification maps keep every retained band at the
                    // occupied threshold unless FD_EMPTY_VERIFY=on (the
                    // verification criterion itself is unchanged).
                    if empty_verify {
                        previous_occupations.as_deref()
                    } else {
                        None
                    },
                ) {
                    Err(ScfError::InvalidInput(m)) if growth.is_some() && m == map::COVERAGE => {
                        grow(
                            &mut states,
                            &mut kept,
                            comm,
                            &mut bands,
                            cap,
                            (iteration, true),
                            shortfall,
                            &mut growth_events,
                        )?;
                    }
                    other => break other?,
                }
            };
            verification_grew = bands != before;
            if verification_grew {
                warm(&mut states, &fresh);
            }
            for cache in &mut kept {
                cache.preserve_input(false);
            }
            collect_calls(&fresh, &mut verification_calls);
            delivered_coverage(&fresh, comm, bands, true)?;
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
            predictor.observe(dr2, error);
            if let Some(trace) = trace.as_deref_mut() {
                trace.records[trace.used - 1].stiffness_ha = stiffness_ha;
            }
            let distance = distances(comm, &step, &fresh, dv, electrons)?;
            let mut eig_change = 0_f64;
            let mut occ_change = 0.;
            let mut eigen_residual = 0_f64;
            for (s, (a, b)) in step.blocks.iter().zip(&fresh.blocks).enumerate() {
                for (k, (a, b)) in a.iter().zip(b).enumerate() {
                    // A verification-map growth adds empty bands absent
                    // from the ordinary map: their occupation counts in full.
                    for i in 0..bands {
                        let before = a.occupations.get(i).copied().unwrap_or(0.);
                        occ_change += points[k].weight * (before - b.occupations[i]).abs();
                        if i >= a.occupations.len() {
                            continue;
                        }
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
            // GV1: a band-count change inside the attempt leaves the pair
            // consistent (both densities are full densities), so it is kept.
            if !passed && feed != super::driver::VerifyFeed::Off && pulay.is_none() {
                fed = Some(residual);
            }
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
            if let Some(g) = &fed {
                mixing = super::driver::feed_attempt(
                    feed,
                    &mut mixer,
                    comm,
                    &mut metric,
                    &y,
                    g,
                    &mut next,
                )?;
            }
            history_conditioning = mixer.conditioning_diagnostics();
            // GM1 (`SCF_ETHR_GUARD` forecast): the next map's eigenvalue threshold
            // follows min(dr2, |F~|^2), the residual the secant model leaves
            // unexplained, so a map is solved finely enough for the residual
            // it is expected to produce, not only for the previous one.
            if forecast_ethr {
                let reduced = metric.norm(comm, mixer.reduced_residual())?;
                if reduced.is_finite() {
                    previous = previous.min(reduced);
                }
            }
            let mut proposal = inverse(&mut plan, comm, &next, policy.channels())?;
            if let Some(mixers) = &mut pulay {
                use dft_solvers::mixing::DensityMixer;
                for s in 0..policy.channels() {
                    proposal[s] = mixers[s]
                        .mix(
                            &n[s],
                            &step.density[s],
                            policy.fixed_count(s) as f64 / volume,
                        )
                        .map_err(super::numerical)?;
                }
            }
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
                    if let Some(mixers) = &mut pulay {
                        use dft_solvers::mixing::DensityMixer;
                        for mixer in mixers {
                            mixer.reset();
                        }
                    }
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
            // After a verification-map growth the states hold the enlarged
            // verification solution; the smaller ordinary block is not restored.
            if !verification_grew {
                warm(&mut states, &step);
            }
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
        let mut sampled = super::stall::MapDiagnostic {
            iteration,
            ..Default::default()
        };
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
                sampled.channel_gap_ha[s] = Some(low - high);
                sampled.channel_edge_residual_ha[s] = Some(
                    channel
                        .iter()
                        .map(|b| b.residuals[n - 1] + b.residuals[n])
                        .fold(0., f64::max),
                );
            }
        }
        if let Some(t) = trace.as_ref() {
            sampled.occupied_subspace_loss =
                t.records[ordinary_trace.unwrap()].occupied_subspace_loss;
        }
        map_diagnostics.push(sampled);
        if super::stall::detected(history.iter().rev().map(|r| (r.dr2_ry, r.mixing.condition))) {
            let diagnostic = super::stall::Diagnostic {
                iteration,
                channel_gap_ha: sampled.channel_gap_ha,
                orbital_rotation_fraction: trace
                    .as_ref()
                    .and_then(|t| t.records[ordinary_trace.unwrap()].rotation_fraction),
            };
            map_report.stall = Some(diagnostic);
            stalls.push(diagnostic);
        }
        map_report.safe_point = true;
        if passed || iteration == opts.max_iterations {
            // The last map of an exhausted budget is delivered too.
            delivered_coverage(&step, comm, bands, false)?;
            boundary.publish(map_report);
            let mut densities = Vec::new();
            let mut blocks = Vec::new();
            for (s, channel) in step.blocks.iter_mut().enumerate() {
                densities.push(states[s].fields.gather(comm, &step.density[s])?);
                let first = states[s].pool.as_ref().map_or(0, |p| p.local().start);
                let local = states[s].blocks.len();
                // GK18B: a pooled run whose caller does not read the orbitals
                // keeps each k block with its owner (materialized on demand by
                // `Run::materialize_orbitals`, same bits as the eager delivery).
                let lazy = if states[s].pool.is_some() && super::lazy::lazy_result_allowed() {
                    let width = channel.get(first).map_or(0, |b| b.orbitals.cols());
                    let owned = channel
                        .iter_mut()
                        .skip(first)
                        .take(local)
                        .map(|b| {
                            ComplexMatrix::zeros(0, 0)
                                .map(|empty| std::mem::replace(&mut b.orbitals, empty))
                        })
                        .collect::<Result<Vec<_>, _>>()
                        .map_err(super::numerical)?;
                    super::k_pool::lazy_orbitals(&states[s], comm, owned, width)?
                } else {
                    None
                };
                let (orbitals, lazy): (Vec<_>, Vec<_>) = match lazy {
                    Some(handles) => handles
                        .into_iter()
                        .map(|h| ComplexMatrix::zeros(0, h.bands()).map(|z| (z, Some(h))))
                        .collect::<Result<Vec<_>, _>>()
                        .map_err(super::numerical)?
                        .into_iter()
                        .unzip(),
                    None => {
                        // GK9 pools: the owner of each k broadcasts its full-row block.
                        let owned: Vec<_> = channel
                            .iter()
                            .skip(first)
                            .take(local)
                            .map(|b| &b.orbitals)
                            .collect();
                        super::k_pool::deliver_orbitals(&states[s], comm, &owned)?
                            .into_iter()
                            .map(|x| (x, None))
                            .unzip()
                    }
                };
                let mut delivered = Vec::new();
                for (k, ((b, orbitals), lazy)) in channel.iter().zip(orbitals).zip(lazy).enumerate()
                {
                    delivered.push(BandBlock {
                        point: points[k],
                        eigenvalues_ha: b.eigenvalues.clone(),
                        occupations: b.occupations.clone(),
                        residuals_ha: b.residuals.clone(),
                        orbitals,
                        lazy,
                    });
                }
                blocks.push(delivered);
            }
            let mut identity = crate::force_identity::ModelDescription::new(
                input.cell,
                input.atoms,
                input.potentials,
                &options,
                dims,
                points.to_vec(),
            );
            if let Some(p) = prepared {
                identity.symmetry = p.operations(opts.crystal_density_group).to_vec();
            }
            return Ok(Run {
                stalls,
                map_diagnostics,
                start_report: refstart::StartReport {
                    // Collinear channels use uniform counts, irrespective of the
                    // atomic density used to prepare their initial orbital span.
                    density: if states.len() == 2 && !atomic_channels {
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
                band_growth: growth.map(|_| BandGrowthReport {
                    initial_bands,
                    final_bands: bands,
                    band_cap: cap,
                    events: growth_events.clone(),
                }),
                coverage_misses: coverage_misses.clone(),
                setup_s,
                total_s: clock.elapsed().as_secs_f64(),
                force_identity: identity.with_policy(occupation_identity(&policy)),
                policy: policy.clone(),
            });
        }
        completed = iteration;
        boundary.emit_state(map_report, checkpoint_state!()).await;
    }
    unreachable!("positive budget checked")
}
/// One automatic growth step after a coverage rejection: extend every channel
/// and k block with the start's kinetic-ordered complements and drop the kept
/// H X entries (their column count changed). Density, Broyden history, ethr
/// and the iteration counter are untouched: the same map is solved again.
#[allow(clippy::too_many_arguments)]
fn grow<K: Collective>(
    states: &mut [refstart::KStartState],
    kept: &mut [KeptHx],
    comm: &K,
    bands: &mut usize,
    cap: usize,
    (iteration, verification): (usize, bool),
    shortfall: map::Shortfall,
    events: &mut Vec<GrowthEvent>,
) -> Result<(), ScfError> {
    let basis = min_basis(&states[0]);
    let next = BandGrowth::next(*bands, cap.min(basis));
    if next <= *bands {
        // ScfError carries a static message; this terminal error is built once.
        let message = format!(
            "{} (automatic band count reached {} bands; occupations.band_cap = {cap}{})",
            map::COVERAGE,
            *bands,
            if basis < cap {
                ", smallest k basis reached"
            } else {
                ""
            }
        );
        return Err(ScfError::InvalidInput(Box::leak(message.into_boxed_str())));
    }
    // Growth applies to the retained count; the FD buffer follows it.
    let width = map::fd_width(next, basis);
    for state in states.iter_mut() {
        refstart::extend_k_orbitals(state, width, comm)?;
    }
    for cache in kept.iter_mut() {
        cache.invalidate();
    }
    if comm.rank() == 0 {
        eprintln!(
            "band growth: iteration {iteration}{}: FD coverage insufficient at {} bands (edge fraction {:.3e}, edge electrons {:.3e}); continuing with {next} bands (cap {cap})",
            if verification {
                " (verification map)"
            } else {
                ""
            },
            *bands,
            shortfall.edge_max_fraction,
            shortfall.edge_electrons,
        );
    }
    events.push(GrowthEvent {
        iteration,
        verification,
        from: *bands,
        to: next,
        edge_max_fraction: shortfall.edge_max_fraction,
        edge_electrons: shortfall.edge_electrons,
    });
    *bands = next;
    Ok(())
}
/// GB1: the hard FD coverage check of a delivered map (fixed band count; with
/// automatic growth such a map never reaches delivery). The refusal message is
/// the unchanged `map::COVERAGE`; the numbers go to stderr on rank 0.
fn delivered_coverage<K: Collective>(
    step: &Map,
    comm: &K,
    bands: usize,
    verification: bool,
) -> Result<(), ScfError> {
    let Some(miss) = step.coverage else {
        return Ok(());
    };
    if comm.rank() == 0 {
        eprintln!(
            "fd coverage: delivered {}map at {bands} bands: top-band occupation fraction {:.3e} (limit {:.0e}), edge electrons {:.3e} (limit {:.0e} N)",
            if verification { "verification " } else { "" },
            miss.edge_max_fraction,
            map::EDGE_FRACTION,
            miss.edge_electrons,
            map::EDGE_ELECTRONS,
        );
    }
    Err(ScfError::InvalidInput(map::COVERAGE))
}
/// Smallest plane-wave basis over the k points of one channel (every global
/// k, also with GK9 pools).
fn min_basis(state: &refstart::KStartState) -> usize {
    state.basis_sizes().into_iter().min().unwrap_or(0)
}
fn collect_calls(step: &Map, out: &mut Vec<(usize, usize, SolveRecord)>) {
    for (s, c) in step.blocks.iter().enumerate() {
        for (k, b) in c.iter().enumerate() {
            out.extend(b.calls.iter().cloned().map(|call| (s, k, call)));
        }
    }
}
/// GF2: per (channel, global k) occupations of the retained bands of a map.
fn occupations_of(step: &Map) -> Vec<Vec<Vec<f64>>> {
    step.blocks
        .iter()
        .map(|channel| channel.iter().map(|b| b.occupations.clone()).collect())
        .collect()
}
fn warm(states: &mut [refstart::KStartState], step: &Map) {
    for (state, channel) in states.iter_mut().zip(&step.blocks) {
        // GK9 pools: this rank's blocks are the global k `local()`.
        let first = state.pool.as_ref().map_or(0, |p| p.local().start);
        for (x, b) in state.orbitals.0.iter_mut().zip(channel.iter().skip(first)) {
            x.clone_from(b.warm_start());
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
