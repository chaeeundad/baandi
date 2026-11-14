//! C1 cutoff ladder (design note docs/roadmap/c1-cutoff-ladder.md): the first
//! part of an SCF runs at a reduced cutoff E1 = fraction * E2 on the smaller FFT
//! grid that E1 implies; the state is then prolonged EXACTLY into the full basis
//! and the unchanged strict session finishes at the full cutoff E2.
//!
//! Exactness of the prolongation: plane-wave spheres are nested, so every coarse
//! coefficient c(G) is copied to the same Miller index of the fine basis and the
//! new G vectors get zero (an injection: orthonormality and the Hermitian Gamma
//! symmetry c(-G) = conj c(G) are preserved bit for bit). The density is a
//! trigonometric polynomial on the coarse density sphere; its Fourier
//! coefficients are copied to the fine grid and zero padded (Fourier
//! interpolation, exact for a band-limited field; the G = 0 term and so the
//! charge are unchanged).
//!
//! The result of the calculation is the full-cutoff session's: its energy audit,
//! eigenpairs and fixed-point certificate are all at the full cutoff; the coarse
//! stage supplies only the start state (density and warm orbitals). By default
//! the mixer history starts empty at the full cutoff: absolute (input, residual)
//! pairs do not translate between the two maps, whose fixed points differ.
//!
//! S2 (`secant_pairs > 0`, `BAANDI_LADDER_SECANT`, default off): the DIFFERENCE
//! pairs of the coarse Pulay history do translate to first order (dr_i is the
//! response of the residual map to dx_i, dominated by the low-G dielectric
//! response that the coarse grid holds). The last h of them are prolonged with
//! the same exact Fourier embedding as the density and seeded into the fine
//! mixer as a closed block of differences (never differenced against a fine
//! pair). The first fine step is a plain evaluation at the prolonged density;
//! its mix uses the seeded least squares with the unchanged spectral cut, and
//! the seeds age out first as fine differences arrive.
//!
//! Scope: one exact Gamma point on the real (ORC) path, fixed uniform integer
//! occupations, no symmetry projector, no frozen-Miller strain reference. The
//! strict profile always; the standard profile only with `standard = true`
//! (`BAANDI_CUTOFF_LADDER_STANDARD=on`, S1, default off): its coarse stage runs
//! the standard profile at E1 (uncertified, stopped by the same residual rule)
//! and the full-cutoff standard session is unchanged except for the first-step
//! eigen tolerance (`StandardSeed`). Everything else runs the ordinary SCF with
//! a recorded reason.
//!
//! General rule for the coarse level (no per-system constant): E1 is a stated
//! fraction of the full cutoff (default 0.5, the design's pre-registered arm);
//! the coarse stage stops when its density residual falls to a fixed level
//! (default 1e-4 Bohr^-3/2, pre-registered in the design note) or after a fixed
//! iteration budget. Neither value is tuned on any system.
use super::*;
use crate::checkpoint::{ElectronicState, Orbitals};
use std::collections::HashMap;

/// First-step eigen tolerance of the full-cutoff standard session after a
/// ladder hand-off.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StandardSeed {
    /// The ratchet cap, exactly as a cold standard session (the first-iteration
    /// redo rule then tightens it if the start is already good).
    Cap,
    /// The tolerance the coarse standard session would have used next:
    /// ratchet(last coarse tolerance, last coarse Coulomb error, N_e). The fine
    /// session's own ratchet takes over from its first record.
    Coarse,
}

/// Ladder controls. Off unless a caller runs `run_ladder`.
#[derive(Clone, Debug)]
pub struct CutoffLadderOptions {
    /// Coarse cutoff as a fraction of the full cutoff, in (0, 1).
    pub fraction: f64,
    /// Hand off when the coarse density residual (Bohr^-3/2) is at or below this.
    pub coarse_density_residual: f64,
    /// Budget of coarse iterations; the state reached is handed off either way.
    pub coarse_max_iterations: usize,
    /// S1: the standard profile is eligible (`BAANDI_CUTOFF_LADDER_STANDARD=on`;
    /// default off).
    pub standard: bool,
    /// S1: first-step tolerance rule of the fine standard session
    /// (`BAANDI_LADDER_STANDARD_SEED=cap|coarse`; default coarse).
    pub standard_seed: StandardSeed,
    /// S2: number h of coarse Pulay difference pairs prolonged into the fine
    /// mixer; 0 is off (`BAANDI_LADDER_SECANT=off|on|<h>`, on means 4; default off).
    pub secant_pairs: usize,
}
impl Default for CutoffLadderOptions {
    /// Pre-registered values (fraction 0.5, residual 1e-4, budget 30) and the
    /// S1/S2 experiment switches read from the environment (an invalid value
    /// means off / the default).
    fn default() -> Self {
        Self {
            standard: dft_numerics::switches::on(
                dft_numerics::switches::Id::CUTOFF_LADDER_STANDARD,
            ),
            standard_seed: match dft_numerics::switches::choice(
                dft_numerics::switches::Id::LADDER_STANDARD_SEED,
            ) {
                "cap" => StandardSeed::Cap,
                _ => StandardSeed::Coarse,
            },
            secant_pairs: match dft_numerics::switches::choice(
                dft_numerics::switches::Id::LADDER_SECANT,
            ) {
                "off" => 0,
                "on" => 4,
                v => v.parse().unwrap_or(0),
            },
            ..Self::baseline()
        }
    }
}
impl CutoffLadderOptions {
    /// The pre-registered ladder with both experiment switches off (no
    /// environment read): the strict-only ladder as adopted.
    pub fn baseline() -> Self {
        Self {
            fraction: 0.5,
            coarse_density_residual: 1e-4,
            coarse_max_iterations: 30,
            standard: false,
            standard_seed: StandardSeed::Coarse,
            secant_pairs: 0,
        }
    }
    pub fn validate(&self) -> Result<(), ScfError> {
        if !(self.fraction.is_finite() && self.fraction > 0.0 && self.fraction < 1.0)
            || !(self.coarse_density_residual.is_finite() && self.coarse_density_residual > 0.0)
            || self.coarse_max_iterations == 0
            || self.secant_pairs > 64
        {
            return Err(ScfError::InvalidInput(
                "cutoff_ladder: fraction in (0,1), positive residual and budget, secant <= 64 required",
            ));
        }
        Ok(())
    }
}

/// Work and outcome of the coarse stage. H columns are counted in the coarse
/// basis; `coarse_grid_points / fine_grid_points` converts them to full-grid
/// equivalents for an H-cost comparison.
#[derive(Clone, Debug, Default)]
pub struct CutoffLadderReport {
    /// True when the full-cutoff session started from the prolonged state.
    pub used: bool,
    pub fallback_reason: Option<String>,
    pub fraction: f64,
    pub coarse_cutoff_hartree: f64,
    pub coarse_grid: [usize; 3],
    pub fine_grid: [usize; 3],
    pub coarse_plane_waves: usize,
    pub fine_plane_waves: usize,
    pub coarse_iterations: usize,
    pub coarse_final_density_residual: f64,
    pub coarse_davidson_iterations: usize,
    /// Ordinary solve H columns of the coarse history (no certificate ran).
    pub coarse_h_vectors: usize,
    /// H columns of the coarse atomic start (counted in the coarse basis).
    pub coarse_start_h_vectors: usize,
    pub coarse_setup: Duration,
    /// Setup, coarse session and prolongation together (unchanged meaning).
    pub coarse_time: Duration,
    /// P48F diagnostics, parts of `coarse_time`: wall of the coarse session
    /// (initialization and steps) and of the hand-off (orbital and density
    /// prolongation with its checks; S2 secant pairs are not included).
    pub coarse_session: Duration,
    pub coarse_prolong: Duration,
    /// Phase timings summed over the coarse history (certificates included).
    pub coarse_phases: PhaseTimings,
    /// Process clocks accumulated inside the coarse session (atomic start,
    /// step wall, certificate blocks).
    pub coarse_clock: crate::phase_clock::PhaseClock,
    /// The full-cutoff session runs the standard profile (S1).
    pub standard: bool,
    /// First-step eigen tolerance given to the fine standard session (Ha);
    /// None means the cap (cold-start rule).
    pub standard_seed_ha: Option<f64>,
    /// S2: requested and seeded difference pairs, and why none were seeded.
    pub secant_requested: usize,
    pub secant_seeded: usize,
    pub secant_fallback: Option<String>,
}
impl CutoffLadderReport {
    /// Coarse H columns scaled by the FFT grid ratio (H cost per column is
    /// dominated by the grid transforms): full-grid H-column equivalents.
    pub fn coarse_h_equivalents(&self) -> f64 {
        let ratio = self.coarse_grid.iter().product::<usize>() as f64
            / self.fine_grid.iter().product::<usize>().max(1) as f64;
        (self.coarse_h_vectors + self.coarse_start_h_vectors) as f64 * ratio
    }
}

/// Why a problem cannot use the ladder (None means eligible).
pub fn ineligibility(problem: &ScfProblem, ladder: &CutoffLadderOptions) -> Option<String> {
    if ladder.validate().is_err() {
        return Some("invalid ladder options".into());
    }
    if problem.states.len() != 1 {
        return Some("multiple k points".into());
    }
    if let Some(reason) = problem.gamma_real_exclusion {
        return Some(reason.into());
    }
    if problem.states[0].basis.k_fractional() != [0.; 3] {
        return Some("k is not exactly Gamma".into());
    }
    let o = &problem.options;
    if o.gamma_real == crate::GammaReal::Off {
        return Some("gamma_real off".into());
    }
    if o.standard.is_some() && !ladder.standard {
        return Some("standard profile".into());
    }
    if o.locking != dft_solvers::locking::Locking::Off
        || o.spectral_guard != crate::SpectralGuard::Off
        || !matches!(o.eigensolver, dft_solvers::SolverChoice::Davidson)
    {
        return Some("solver options outside the real Gamma path".into());
    }
    if o.strain_reference.is_some() {
        return Some("frozen-Miller strain reference".into());
    }
    if o.compat_mode != CompatMode::None {
        return Some("compatibility mode".into());
    }
    if problem.symmetry.is_some() {
        return Some("density symmetrization enabled".into());
    }
    let occupied = problem.occupations.iter().filter(|&&f| f > 0.).count();
    if occupied == 0
        || problem
            .occupations
            .iter()
            .any(|&f| f != 0. && f != problem.occupations[0])
        || problem.occupations[..occupied].iter().any(|&f| f <= 0.)
    {
        return Some("occupations are not fixed and uniform".into());
    }
    None
}

/// Copy coefficients by Miller index into the fine basis (zero for new G).
/// Fails if a coarse G is missing from the fine basis (bases not nested).
pub fn prolong_orbitals(
    coarse: &PlaneWaveBasis,
    fine: &PlaneWaveBasis,
    orbitals: &ComplexMatrix,
) -> Result<ComplexMatrix, ScfError> {
    if orbitals.rows() != coarse.vectors().len() {
        return Err(ScfError::InvalidInput(
            "cutoff_ladder: orbital rows differ from the coarse basis",
        ));
    }
    let index: HashMap<[i32; 3], usize> = fine
        .vectors()
        .iter()
        .enumerate()
        .map(|(i, g)| (g.miller, i))
        .collect();
    let map = coarse
        .vectors()
        .iter()
        .map(|g| index.get(&g.miller).copied())
        .collect::<Option<Vec<_>>>()
        .ok_or(ScfError::InvalidInput(
            "cutoff_ladder: coarse basis is not a subset of the fine basis",
        ))?;
    let rows = fine.vectors().len();
    let mut out = ComplexMatrix::zeros(rows, orbitals.cols())
        .map_err(|e| ScfError::Numerical(format!("cutoff_ladder: {e}")))?;
    for j in 0..orbitals.cols() {
        let dst = out.column_mut(j);
        for (src, &i) in orbitals.column(j).iter().zip(&map) {
            dst[i] = *src;
        }
    }
    Ok(out)
}

/// Fourier interpolation of a real periodic field from the coarse FFT grid to a
/// finer one (zero padding). The Nyquist plane of an even coarse dimension is
/// dropped: a field band-limited to the coarse density sphere has no weight
/// there because the grid holds that sphere (dims >= 2 max|m| + 1).
pub fn prolong_density(
    coarse: &Array3<f64>,
    fine_dims: [usize; 3],
) -> Result<Array3<f64>, ScfError> {
    let cd = [coarse.dim().0, coarse.dim().1, coarse.dim().2];
    let mut fft = FftPlan3d::new(cd)?;
    let input = coarse.mapv(|v| Complex64::new(v, 0.0));
    let mut reciprocal = Array3::<Complex64>::zeros(cd);
    fft.to_reciprocal(&input, &mut reciprocal)?;
    let mut padded = Array3::<Complex64>::zeros(fine_dims);
    let signed = |i: usize, n: usize| -> Option<i32> {
        if n.is_multiple_of(2) && i == n / 2 {
            None
        } else if i <= (n - 1) / 2 {
            Some(i as i32)
        } else {
            Some(i as i32 - n as i32)
        }
    };
    for ((i, j, k), z) in reciprocal.indexed_iter() {
        let (Some(a), Some(b), Some(c)) = (signed(i, cd[0]), signed(j, cd[1]), signed(k, cd[2]))
        else {
            continue;
        };
        let [p, q, r] = FftGridSpec::map_miller([a, b, c], fine_dims).map_err(|_| {
            ScfError::InvalidInput("cutoff_ladder: coarse frequencies exceed the fine grid")
        })?;
        padded[[p, q, r]] = *z;
    }
    let mut real = Array3::<Complex64>::zeros(fine_dims);
    let mut fine_fft = FftPlan3d::new(fine_dims)?;
    fine_fft.to_real(&padded, &mut real)?;
    Ok(real.mapv(|z| z.re))
}

fn coarse_problem(
    fine: &ScfProblem,
    potentials: &[SpeciesPotential<'_>],
    ladder: &CutoffLadderOptions,
) -> Result<ScfProblem, ScfError> {
    let mut options = fine.options.clone();
    options.cutoff_hartree *= ladder.fraction;
    // A fixed grid belongs to the full cutoff; the coarse grid follows E1.
    options.grid_override = None;
    // ADR 0043: the coarse stage only supplies a start; verification belongs to the fine run.
    options.verification = crate::Verification::Off;
    let mut coarse = ScfProblem::new(&fine.cell, &fine.atoms, potentials, options)?;
    if let Some(floor) = fine.density_floor {
        coarse.set_density_floor(floor)?;
    }
    Ok(coarse)
}

/// What the coarse stage hands to the full cutoff.
struct CoarseState {
    density: Array3<f64>,
    orbitals: ComplexMatrix,
    basis: PlaneWaveBasis,
    /// Pulay memory of the coarse session (only when S2 asks for it).
    mixer: Option<crate::checkpoint::MixerState>,
    /// Standard profile: the tolerance the coarse ratchet would use next.
    standard_next: Option<f64>,
}

/// Coarse stage: drive the ordinary session at E1 and return its next input
/// density and the warm orbitals of its last solve.
fn coarse_stage(
    coarse: ScfProblem,
    ladder: &CutoffLadderOptions,
    report: &mut CutoffLadderReport,
) -> Result<CoarseState, ScfError> {
    // P0.2: Davidson and operator work of this stage is reported as ladder coarse.
    let _run_stage = dft_numerics::work_timers::enter(dft_numerics::work_timers::RunStage::Coarse);
    let basis = coarse.states[0].basis.clone();
    let session_start = Instant::now();
    let clock = crate::phase_clock::snapshot();
    let mut session = ScfSession::new(
        coarse,
        CompletionMode::RetainLast,
        crate::events::ExecutionConfig::default(),
    )?;
    session.initialize()?;
    loop {
        let status = session.step()?.status;
        let last = session
            .history()
            .last()
            .map_or(f64::INFINITY, |r| r.density_residual);
        if status != StepStatus::Ready
            || last <= ladder.coarse_density_residual
            || session.history().len() >= ladder.coarse_max_iterations
        {
            break;
        }
    }
    report.coarse_session = session_start.elapsed();
    report.coarse_clock = crate::phase_clock::snapshot().since(&clock);
    let history = session.history();
    for t in history.iter().flat_map(|r| {
        std::iter::once(&r.timings).chain(r.fixed_point_certificate.as_ref().map(|c| &c.timings))
    }) {
        report.coarse_phases.potential += t.potential;
        report.coarse_phases.diagonalization += t.diagonalization;
        report.coarse_phases.density += t.density;
        report.coarse_phases.energy += t.energy;
        report.coarse_phases.mixing += t.mixing;
    }
    report.coarse_iterations = history.len();
    report.coarse_final_density_residual = history.last().map_or(f64::NAN, |r| r.density_residual);
    report.coarse_davidson_iterations = history.iter().map(|r| r.davidson_iterations).sum();
    report.coarse_h_vectors = history
        .iter()
        .map(|r| {
            r.operator_applications
                + r.fixed_point_certificate
                    .as_ref()
                    .map_or(0, |c| c.operator_applications)
        })
        .sum();
    report.coarse_start_h_vectors = session.problem_ref().atomic_start.initial_h_applications;
    let standard_next = history.last().and_then(|r| {
        r.standard_error_ha.map(|e| {
            StandardProfile::ratchet(
                r.requested_eigen_tolerance,
                e,
                session.problem_ref().electrons,
            )
        })
    });
    let mixer = if ladder.secant_pairs > 0 {
        session.mixer_state()
    } else {
        None
    };
    let (density, orbitals) = match session.warm_state() {
        Some((n, Some(orbitals))) if orbitals.len() == 1 => (n.clone(), orbitals[0].clone()),
        _ => {
            // Terminal coarse session (converged or exhausted): use its result.
            let outcome = session.result().ok_or(ScfError::InvalidInput(
                "cutoff_ladder: coarse state unavailable",
            ))?;
            (
                outcome.result.density.clone(),
                outcome.result.k_points[0].solution.eigenvectors.clone(),
            )
        }
    };
    Ok(CoarseState {
        density,
        orbitals,
        basis,
        mixer,
        standard_next,
    })
}

/// Prolonged difference pairs (dx_i, dr_i), oldest first.
type SecantPairs = (Vec<Vec<f64>>, Vec<Vec<f64>>);

/// S2: the last `h` difference pairs (dx, dr) of the coarse Pulay memory,
/// oldest first, each prolonged to the fine grid by the exact Fourier
/// embedding. Differences are formed inside the coarse history only.
fn prolonged_secant_pairs(
    mixer: Option<&crate::checkpoint::MixerState>,
    coarse_dims: [usize; 3],
    fine_dims: [usize; 3],
    h: usize,
) -> Result<SecantPairs, String> {
    let Some(crate::checkpoint::MixerState::Pulay(state)) = mixer else {
        return Err("coarse mixer is not plain Pulay".into());
    };
    let size = coarse_dims.iter().product::<usize>();
    let pairs = state.inputs.len().min(state.residuals.len());
    if pairs < 2 {
        return Err("coarse history has no difference pair".into());
    }
    let first = (pairs - 1).saturating_sub(h);
    let prolong = |a: &[f64], b: &[f64]| -> Result<Vec<f64>, String> {
        if a.len() != size || b.len() != size {
            return Err("coarse history size differs from the coarse grid".into());
        }
        let d = Array3::from_shape_vec(coarse_dims, a.iter().zip(b).map(|(x, y)| x - y).collect())
            .map_err(|e| e.to_string())?;
        Ok(prolong_density(&d, fine_dims)
            .map_err(|e| e.to_string())?
            .into_iter()
            .collect())
    };
    let mut dx = Vec::with_capacity(pairs - 1 - first);
    let mut dr = Vec::with_capacity(pairs - 1 - first);
    for i in first..pairs - 1 {
        dx.push(prolong(&state.inputs[i + 1], &state.inputs[i])?);
        dr.push(prolong(&state.residuals[i + 1], &state.residuals[i])?);
    }
    Ok((dx, dr))
}

/// Result of a ladder run: the full-cutoff outcome and the coarse-stage report.
pub struct LadderOutcome {
    pub outcome: ScfOutcome,
    pub report: CutoffLadderReport,
}

/// Run the ladder on an eligible full-cutoff problem, or the ordinary SCF with
/// the reason recorded. `potentials` must be the ones `fine` was built from.
/// Only the full-cutoff session's certificate can set convergence.
pub fn run_ladder(
    fine: ScfProblem,
    potentials: &[SpeciesPotential<'_>],
    ladder: &CutoffLadderOptions,
    retain_last: bool,
) -> Result<LadderOutcome, ScfError> {
    let mut report = CutoffLadderReport {
        standard: fine.options.standard.is_some(),
        secant_requested: ladder.secant_pairs,
        fraction: ladder.fraction,
        coarse_cutoff_hartree: fine.options.cutoff_hartree * ladder.fraction,
        fine_grid: fine.grid.dims(),
        fine_plane_waves: fine.states[0].basis.vectors().len(),
        ..Default::default()
    };
    if let Some(reason) = ineligibility(&fine, ladder) {
        report.fallback_reason = Some(reason);
        let outcome = direct_min::finish(fine, None, retain_last)?;
        return Ok(LadderOutcome { outcome, report });
    }
    let setup = Instant::now();
    let start = match coarse_problem(&fine, potentials, ladder) {
        Ok(coarse) => {
            report.coarse_setup = setup.elapsed();
            report.coarse_grid = coarse.grid.dims();
            report.coarse_plane_waves = coarse.states[0].basis.vectors().len();
            coarse_stage(coarse, ladder, &mut report).and_then(|c| {
                let prolong = Instant::now();
                let state = hand_off(&fine, &c.basis, &c.density, &c.orbitals);
                report.coarse_prolong = prolong.elapsed();
                state.map(|state| (state, c))
            })
        }
        Err(e) => Err(e),
    };
    let (electronic, seeds) = match start {
        Ok((state, coarse)) => {
            report.used = true;
            if fine.options.standard.is_some() && ladder.standard_seed == StandardSeed::Coarse {
                report.standard_seed_ha = coarse.standard_next;
            }
            let seeds = if ladder.secant_pairs > 0 {
                match prolonged_secant_pairs(
                    coarse.mixer.as_ref(),
                    report.coarse_grid,
                    fine.grid.dims(),
                    ladder.secant_pairs,
                ) {
                    Ok(pairs) => Some(pairs),
                    Err(reason) => {
                        report.secant_fallback = Some(reason);
                        None
                    }
                }
            } else {
                None
            };
            (Some(state), seeds)
        }
        Err(e) => {
            report.fallback_reason = Some(format!("coarse stage: {e}"));
            (None, None)
        }
    };
    report.coarse_time = setup.elapsed();
    let Some(electronic) = electronic else {
        let outcome = direct_min::finish(fine, None, retain_last)?;
        return Ok(LadderOutcome { outcome, report });
    };
    let outcome = finish_from(fine, &electronic, seeds, retain_last, &mut report)?;
    Ok(LadderOutcome { outcome, report })
}

/// `direct_min::finish` from the prolonged state, plus the optional S1
/// standard first-step tolerance and S2 mixer seeds (each recorded; a rejected
/// seed is a recorded fallback, never an error). Only the session's own
/// convergence rule ends the run.
fn finish_from(
    fine: ScfProblem,
    start: &ElectronicState,
    seeds: Option<SecantPairs>,
    retain_last: bool,
    report: &mut CutoffLadderReport,
) -> Result<ScfOutcome, ScfError> {
    if report.standard_seed_ha.is_none() && seeds.is_none() {
        return direct_min::finish(fine, Some(start), retain_last);
    }
    let completion = if retain_last {
        CompletionMode::RetainLast
    } else {
        CompletionMode::Strict
    };
    let mut session = ScfSession::new(fine, completion, crate::events::ExecutionConfig::default())?;
    session.initialize_electronic(start)?;
    if let Some(t) = report.standard_seed_ha {
        session.seed_standard_tolerance(t)?;
    }
    if let Some((dx, dr)) = seeds {
        match session.seed_mixer_differences(&dx, &dr) {
            Ok(()) => report.secant_seeded = dx.len(),
            Err(e) => report.secant_fallback = Some(e.to_string()),
        }
    }
    while session.step()?.status == StepStatus::Ready {}
    if !retain_last && session.result().is_none() {
        return Err(ScfError::NotConverged {
            history: session.history().to_vec(),
        });
    }
    session
        .into_result()
        .ok_or(ScfError::InvalidInput("result unavailable"))
}

/// Prolonged density and orbitals as the start state of the full-cutoff session,
/// with the session's own import checks repeated here so that a rejection is a
/// recorded fallback, not an error.
fn hand_off(
    fine: &ScfProblem,
    coarse_basis: &PlaneWaveBasis,
    density: &Array3<f64>,
    orbitals: &ComplexMatrix,
) -> Result<ElectronicState, ScfError> {
    let fine_basis = &fine.states[0].basis;
    let prolonged = prolong_orbitals(coarse_basis, fine_basis, orbitals)?;
    let n = prolong_density(density, fine.grid.dims())?;
    if n.iter()
        .any(|v| !v.is_finite() || *v < -xc::NEGATIVE_ROUNDOFF)
    {
        return Err(ScfError::InvalidInput(
            "cutoff_ladder: prolonged density is negative or nonfinite",
        ));
    }
    let charge = n.sum() * fine.cell.volume_bohr3() / n.len() as f64;
    if (charge - fine.electrons).abs() / fine.electrons > 1e-8 {
        return Err(ScfError::InvalidInput(
            "cutoff_ladder: prolonged density charge error exceeds 1e-8",
        ));
    }
    let mut coefficients = Vec::with_capacity(2 * prolonged.rows() * prolonged.cols());
    coefficients.extend(prolonged.as_slice().iter().flat_map(|z| [z.re, z.im]));
    Ok(ElectronicState {
        grid: fine.grid.dims(),
        density: n.iter().copied().collect(),
        orbitals: vec![Orbitals {
            rows: prolonged.rows(),
            cols: prolonged.cols(),
            coefficients,
        }],
        occupations: fine.occupations.clone(),
    })
}
