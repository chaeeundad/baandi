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
//! stage supplies only the start state (density and warm orbitals). The mixer
//! history starts empty at the full cutoff (histories do not translate between
//! bases).
//!
//! Scope: one exact Gamma point on the real (ORC) path, fixed uniform integer
//! occupations, no symmetry projector, no frozen-Miller strain reference, strict
//! profile only. Everything else runs the ordinary SCF with a recorded reason.
//!
//! General rule for the coarse level (no per-system constant): E1 is a stated
//! fraction of the full cutoff (default 0.5, the design's pre-registered arm);
//! the coarse stage stops when its density residual falls to a fixed level
//! (default 1e-4 Bohr^-3/2, pre-registered in the design note) or after a fixed
//! iteration budget. Neither value is tuned on any system.
use super::*;
use crate::checkpoint::{ElectronicState, Orbitals};
use std::collections::HashMap;

/// Ladder controls. Off unless a caller runs `run_ladder`.
#[derive(Clone, Debug)]
pub struct CutoffLadderOptions {
    /// Coarse cutoff as a fraction of the full cutoff, in (0, 1).
    pub fraction: f64,
    /// Hand off when the coarse density residual (Bohr^-3/2) is at or below this.
    pub coarse_density_residual: f64,
    /// Budget of coarse iterations; the state reached is handed off either way.
    pub coarse_max_iterations: usize,
}
impl Default for CutoffLadderOptions {
    fn default() -> Self {
        Self {
            fraction: 0.5,
            coarse_density_residual: 1e-4,
            coarse_max_iterations: 30,
        }
    }
}
impl CutoffLadderOptions {
    pub fn validate(&self) -> Result<(), ScfError> {
        if !(self.fraction.is_finite() && self.fraction > 0.0 && self.fraction < 1.0)
            || !(self.coarse_density_residual.is_finite() && self.coarse_density_residual > 0.0)
            || self.coarse_max_iterations == 0
        {
            return Err(ScfError::InvalidInput(
                "cutoff_ladder: fraction in (0,1), positive residual and budget required",
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
    pub coarse_time: Duration,
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
    if o.standard.is_some() {
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
    let mut coarse = ScfProblem::new(&fine.cell, &fine.atoms, potentials, options)?;
    if let Some(floor) = fine.density_floor {
        coarse.set_density_floor(floor)?;
    }
    Ok(coarse)
}

/// Coarse stage: drive the ordinary session at E1 and return its next input
/// density and the warm orbitals of its last solve.
fn coarse_stage(
    coarse: ScfProblem,
    ladder: &CutoffLadderOptions,
    report: &mut CutoffLadderReport,
) -> Result<(Array3<f64>, ComplexMatrix, PlaneWaveBasis), ScfError> {
    let basis = coarse.states[0].basis.clone();
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
    let history = session.history();
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
    Ok((density, orbitals, basis))
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
            coarse_stage(coarse, ladder, &mut report)
                .and_then(|(n, c, basis)| hand_off(&fine, &basis, &n, &c))
        }
        Err(e) => Err(e),
    };
    report.coarse_time = setup.elapsed();
    let electronic = match start {
        Ok(state) => {
            report.used = true;
            Some(state)
        }
        Err(e) => {
            report.fallback_reason = Some(format!("coarse stage: {e}"));
            None
        }
    };
    let outcome = direct_min::finish(fine, electronic.as_ref(), retain_last)?;
    Ok(LadderOutcome { outcome, report })
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
