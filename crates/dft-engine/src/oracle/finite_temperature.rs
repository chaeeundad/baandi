//! Explicit finite-T driver. Legacy sessions and checkpoints retain their policy.
use super::*;
use dft_physics::occupations::Constraint;

pub use crate::model::finite_temperature::*;
pub(crate) use crate::model::spin::validate_options;
pub struct FiniteTProblem {
    problem: ScfProblem,
    occupations: OccupationWorkspace,
}
impl FiniteTProblem {
    /// Preparation chooses the explicit cap, including automatic empty states.
    /// All retained states are solved strictly. This conservative initial policy
    /// avoids reallocating band workspaces during the SCF loop.
    pub fn new(
        cell: &Cell,
        atoms: &[Atom],
        potentials: &[SpeciesPotential<'_>],
        mut options: ScfOptions,
        points: &[KPoint],
        theta_ha: f64,
        band_cap: usize,
    ) -> Result<Self, ScfError> {
        validate_options(&options, theta_ha, band_cap)?;
        let mut gamma_real = crate::GammaRealDiagnostics::default();
        if options.gamma_real != crate::GammaReal::Off {
            gamma_real.fallback(options.gamma_real, "Fermi-Dirac occupations")?;
        }
        options.bands = Some(band_cap);
        let mut problem =
            ScfProblem::construct(cell, atoms, potentials, options, points, None, Some(0))?;
        problem.gamma_real_exclusion = Some("Fermi-Dirac occupations");
        problem.clear_atomic_start("Fermi-Dirac occupations");
        for state in &mut problem.states {
            state.gamma_real = gamma_real.clone();
        }
        // Used only to mark every retained eigenpair strict, never for density.
        problem.occupations.fill(1.0);
        let occupations =
            OccupationWorkspace::new(band_cap, points, false, theta_ha, Constraint::Total)?;
        if problem.electrons >= 2.0 * band_cap as f64 {
            return Err(ScfError::InvalidInput(
                "occupations.bands: cap has no empty-state capacity; increase cap",
            ));
        }
        Ok(Self {
            problem,
            occupations,
        })
    }
    pub fn enable_symmetry(
        &mut self,
        full: &[KPoint],
        tolerance_bohr: f64,
        time_reversal: bool,
    ) -> Result<(), ScfError> {
        self.problem
            .enable_symmetry(full, tolerance_bohr, time_reversal)
    }
    /// Reuse the same verified preparation as the input adapter.
    pub fn enable_prepared_symmetry(
        &mut self,
        prepared: &PreparedSymmetry,
    ) -> Result<(), ScfError> {
        self.problem.enable_prepared_symmetry(prepared)
    }
    pub fn run(mut self) -> Result<FiniteTResult, ScfError> {
        let p = &mut self.problem;
        let opts = p.options.clone();
        let volume = p.cell.volume_bohr3();
        let n = p.electrons;
        let mut density = Array3::from_elem(p.grid.dims(), n / volume);
        let mut mixer = super::collinear::make_spin_mixer(p)?;
        let mut previous: Option<Vec<ComplexMatrix>> = None;
        let mut last_free: Option<f64> = None;
        let mut best: f64 = f64::INFINITY;
        let mut growing = 0;
        let mut saved = vec![0.0; self.occupations.occupations.len()];
        let mut requested_eigen_tolerances = Vec::new();
        for iteration in 1..=opts.max_iterations {
            let tolerance = { 0.1 * opts.eigen_tolerance };
            let state = p.step_with_occupations(
                &density,
                previous.as_deref(),
                tolerance,
                true,
                Some(&mut self.occupations),
            )?;

            requested_eigen_tolerances.push(tolerance);
            let evaluation = self.occupations.evaluation.unwrap();
            if evaluation.edge_max_fraction > 1e-8 || evaluation.edge_electrons > 1e-10 * n.max(1.0)
            {
                return Err(ScfError::InvalidInput(
                    "occupations.bands: insufficient empty-state coverage at explicit cap; increase cap and repeat the convergence study",
                ));
            }
            let ledger =
                ThermodynamicLedger::new(state.energy.total_hartree(), evaluation.entropy_term_ha);
            let residual = density_distance(&density, &state.density, volume, n);
            if residual > opts.divergence_factor * best.max(opts.density_tolerance) {
                growing += 1;
            } else {
                growing = 0;
            }
            best = best.min(residual);
            if growing >= opts.divergence_patience {
                return Err(ScfError::Numerical(
                    "finite-T SCF density residual exceeded configured growth limit".into(),
                ));
            }
            let candidate = {
                residual <= opts.density_tolerance
                    && last_free
                        .is_some_and(|f| (ledger.free_energy_ha - f).abs() <= opts.energy_tolerance)
            };
            let mut certificate = None;
            let mut passed = false;
            if candidate {
                saved.copy_from_slice(self.occupations.occupations());
                let fresh = p.step_with_occupations(
                    &state.density,
                    None,
                    0.1 * opts.eigen_tolerance,
                    true,
                    Some(&mut self.occupations),
                )?;
                let fresh_eval = self.occupations.evaluation.unwrap();
                let nr = density_distance(&state.density, &fresh.density, volume, n);
                let er = fresh
                    .k_points
                    .iter()
                    .flat_map(|k| k.solution.residuals.iter())
                    .copied()
                    .fold(0.0, f64::max);
                let df = (fresh.energy.total_hartree()
                    - fresh_eval.entropy_term_ha
                    - ledger.free_energy_ha)
                    .abs();
                let change = saved
                    .iter()
                    .zip(self.occupations.occupations())
                    .enumerate()
                    .map(|(i, (a, b))| {
                        self.occupations.weights[i / self.occupations.bands] * (a - b).abs()
                    })
                    .sum::<f64>();
                passed = nr <= opts.density_tolerance
                    && er <= opts.eigen_tolerance
                    && df <= opts.energy_tolerance
                    && change <= 1e-8 * n.max(1.0)
                    && fresh_eval.electron_residual.abs() <= 1e-10 * n.max(1.0)
                    && (fresh.density.sum() * volume / fresh.density.len() as f64 - n).abs()
                        <= 1e-10 * n.max(1.0)
                    && fresh_eval.edge_max_fraction <= 1e-8
                    && fresh_eval.edge_electrons <= 1e-10 * n.max(1.0);
                certificate = Some(FiniteTCertificate {
                    passed,
                    charge_density_residual: nr,
                    spin_density_residual: None,
                    moment_residual: None,
                    retained_eigen_residual_ha: er,
                    occupation_change_electrons: change,
                    electron_residual: fresh_eval.electron_residual,
                    free_energy_change_ha: df,
                    edge_max_fraction: fresh_eval.edge_max_fraction,
                    edge_electrons: fresh_eval.edge_electrons,
                });
                self.occupations.occupations.copy_from_slice(&saved);
                self.occupations.evaluation = Some(evaluation);
            }
            if passed || iteration == opts.max_iterations {
                return Ok(FiniteTResult {
                    gamma_real_used: false,
                    gamma_real_fallback_reason: self
                        .problem
                        .states
                        .iter()
                        .find_map(|s| s.gamma_real.fallback_reason.clone()),
                    force_identity: self.problem.force_model.with_policy(
                        crate::force_identity::OccupationIdentity::new(
                            &[],
                            Some(&self.occupations),
                            None,
                        ),
                    ),
                    converged: passed,
                    iterations: iteration,
                    state,
                    occupations: self.occupations,
                    thermodynamics: ledger,
                    certificate,
                    standard: None,
                    verification: None,
                    standard_errors_ha: Vec::new(),
                    requested_eigen_tolerances,
                });
            }
            let mut proposal = mixer.mix(
                density.as_slice().unwrap(),
                state.density.as_slice().unwrap(),
                n / volume,
            )?;
            if matches!(
                opts.mixing,
                Mixing::Pulay { .. } | Mixing::KerkerPulay { .. }
            ) {
                let fraction = density
                    .iter()
                    .zip(&proposal)
                    .filter(|(_, b)| **b < 0.0)
                    .map(|(a, b)| 0.9 * a / (a - b))
                    .fold(1.0, f64::min);
                if fraction < 1.0 {
                    for (b, a) in proposal.iter_mut().zip(&density) {
                        *b = a + fraction * (*b - a);
                    }
                    mixer.reset();
                }
            }
            density.as_slice_mut().unwrap().copy_from_slice(&proposal);
            previous = Some(
                state
                    .k_points
                    .into_iter()
                    .map(|k| k.solution.eigenvectors)
                    .collect(),
            );
            last_free = Some(ledger.free_energy_ha);
        }
        unreachable!()
    }
}

#[cfg(test)]
#[path = "finite_temperature_tests.rs"]
mod tests;
