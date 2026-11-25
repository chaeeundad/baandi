//! Fixed-default initialization of the strict oracle (D2).
//! The coarse half-cutoff stage supplies only density and orbitals. Its mixer
//! history and tolerances are discarded. Experiment variants and switch reads
//! are intentionally absent; the reference driver uses its own ladder.
use super::*;
use crate::checkpoint::{ElectronicState, Orbitals};
pub use crate::refscf::prolong::{prolong_density, prolong_orbitals};

fn ineligibility(problem: &ScfProblem) -> Option<String> {
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

fn coarse_problem(
    fine: &ScfProblem,
    potentials: &[SpeciesPotential<'_>],
) -> Result<ScfProblem, ScfError> {
    let mut options = fine.options.clone();
    options.cutoff_hartree *= 0.5;
    // A fixed grid belongs to the full cutoff; the coarse grid follows E1.
    options.grid_override = None;
    // ADR 0043: the coarse stage only supplies a start; verification belongs to the fine run.
    options.verification = crate::Verification::Off;
    // RUN-CONFIG: the coarse problem is constructed and run under the fine
    // problem's configuration.
    let mut coarse = {
        let _run = fine.run_config().enter();
        ScfProblem::new(&fine.cell, &fine.atoms, potentials, options)?
    };
    coarse.set_run_config(fine.run_config());
    if let Some(floor) = fine.density_floor {
        coarse.set_density_floor(floor)?;
    }
    Ok(coarse)
}

fn coarse_start(coarse: ScfProblem, fine: &ScfProblem) -> Result<ElectronicState, ScfError> {
    let _stage = dft_numerics::work_timers::enter(dft_numerics::work_timers::RunStage::Coarse);
    let _profile = crate::profile::enter(crate::profile::Top::Coarse);
    let _context = crate::profile::context_stage("coarse");
    let basis = coarse.states[0].basis.clone();
    let mut session = ScfSession::new(
        coarse,
        CompletionMode::RetainLast,
        crate::events::ExecutionConfig::default(),
    )?;
    session.initialize()?;
    loop {
        let status = session.step()?.status;
        let residual = session
            .history()
            .last()
            .map_or(f64::INFINITY, |r| r.density_residual);
        if status != StepStatus::Ready || residual <= 1e-4 || session.history().len() >= 30 {
            break;
        }
    }
    let (density, orbitals) = match session.warm_state()? {
        Some((n, Some(mut orbitals))) if orbitals.len() == 1 => {
            (n.clone(), orbitals.pop().expect("one block"))
        }
        _ => {
            let outcome = session.result().ok_or(ScfError::InvalidInput(
                "cutoff_ladder: coarse state unavailable",
            ))?;
            let k = &outcome.result.k_points[0];
            let orbitals = match &k.lazy {
                Some(lazy) => lazy.materialize()?.clone(),
                None => k.solution.eigenvectors.clone(),
            };
            (outcome.result.density.clone(), orbitals)
        }
    };
    hand_off(fine, &basis, &density, &orbitals)
}

/// Run the strict oracle with its historical default half-cutoff start.
/// Ineligible or unsuccessful coarse stages retain the ordinary strict start.
pub fn run_oracle(
    fine: ScfProblem,
    potentials: &[SpeciesPotential<'_>],
) -> Result<ScfOutcome, ScfError> {
    let _run = fine.run_config().enter();
    if ineligibility(&fine).is_some() {
        return fine.run_outcome();
    }
    let start = coarse_problem(&fine, potentials).and_then(|coarse| coarse_start(coarse, &fine));
    let Ok(start) = start else {
        return fine.run_outcome();
    };
    let mut session = ScfSession::new(
        fine,
        CompletionMode::RetainLast,
        crate::events::ExecutionConfig::default(),
    )?;
    session.initialize_oracle_ladder(&start)?;
    while session.step()?.status == StepStatus::Ready {}
    session
        .into_result()
        .ok_or(ScfError::InvalidInput("result unavailable"))
}

fn hand_off(
    fine: &ScfProblem,
    coarse_basis: &PlaneWaveBasis,
    density: &Array3<f64>,
    orbitals: &ComplexMatrix,
) -> Result<ElectronicState, ScfError> {
    let fine_basis = &fine.states[0].basis;
    let prolonged = prolong_orbitals(coarse_basis, fine_basis, orbitals)?;
    let state = handoff_density(fine, density)?;
    let mut coefficients = Vec::with_capacity(2 * prolonged.rows() * prolonged.cols());
    coefficients.extend(prolonged.as_slice().iter().flat_map(|z| [z.re, z.im]));
    Ok(ElectronicState {
        grid: fine.grid.dims(),
        density: state.density,
        orbitals: vec![Orbitals {
            rows: prolonged.rows(),
            cols: prolonged.cols(),
            coefficients,
        }],
        occupations: fine.occupations.clone(),
    })
}

fn handoff_density(fine: &ScfProblem, density: &Array3<f64>) -> Result<ElectronicState, ScfError> {
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
    Ok(ElectronicState {
        grid: fine.grid.dims(),
        density: n.iter().copied().collect(),
        orbitals: Vec::new(),
        occupations: fine.occupations.clone(),
    })
}
