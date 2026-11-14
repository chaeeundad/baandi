//! Additive standard-profile controls and reusable signed Coulomb metric.
use crate::scf::*;
use dft_solvers::mixing::ResidualMetric;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StandardProfile {
    /// Total Coulomb residual energy, Ha; never normalized by electron count.
    pub threshold_ha: f64,
}
impl Default for StandardProfile {
    fn default() -> Self {
        Self { threshold_ha: 5e-9 }
    }
}
impl StandardProfile {
    pub const CAP_HA: f64 = 5e-3;
    pub const FLOOR_HA: f64 = 5e-14;
    pub fn ratchet(previous: f64, error_ha: f64, electrons: f64) -> f64 {
        previous
            .min(0.1 * error_ha / electrons.max(1.0))
            .clamp(Self::FLOOR_HA, Self::CAP_HA)
    }
    /// Experimental inexact-solve cap: square the previous normalized density
    /// residual (in atomic units) to map a norm scale to a Ritz-energy scale.
    /// This heuristic is not an energy bound or an inexact-Newton guarantee.
    /// Never loosen the ordinary ratchet; retain its existing numerical floor.
    pub(crate) fn density_forcing(tolerance: f64, history: &[IterationRecord]) -> f64 {
        history.last().map_or(tolerance, |r| {
            tolerance
                .min(r.density_residual.powi(2))
                .max(Self::FLOOR_HA)
        })
    }

    pub fn empty_threshold(tolerance: f64) -> f64 {
        (5.0 * tolerance).max(5e-6)
    }

    /// Five completed intervals without a factor-two error reduction.
    /// Cap against the window-entry tolerance, not the latest tolerance: this
    /// prevents repeated cuts on every step of the same plateau. History alone
    /// determines the action, including after an exact checkpoint continuation.
    pub(crate) fn stall_guard(tolerance: f64, history: &[IterationRecord]) -> f64 {
        const INTERVALS: usize = 5;
        const TIGHTEN: f64 = 0.1;
        if history.len() <= INTERVALS {
            return tolerance;
        }
        let window = &history[history.len() - INTERVALS - 1..];
        let Some(entry_error) = window[0].standard_error_ha else {
            return tolerance;
        };
        if window[1..]
            .iter()
            .all(|r| r.standard_error_ha.is_some_and(|e| e > 0.5 * entry_error))
        {
            tolerance
                .min(TIGHTEN * window[0].requested_eigen_tolerance)
                .max(Self::FLOOR_HA)
        } else {
            tolerance
        }
    }

    /// Residual guard (`BAANDI_STANDARD_RESIDUAL_GUARD`, default off): a
    /// second-order energy-error estimate from the per-band residual norms of
    /// the latest eigensolve, in Ha:
    ///
    ///   B = c * sum_k w_k sum_i f_i ||r_ik||^2,   c = 1 Ha^-1.
    ///
    /// Derivation (general, no material input): for a normalized Ritz vector x
    /// with Ritz value theta and residual r = (H - theta) x, the eigenvector
    /// error is delta ~ (H - lambda)^-1 r and both the Ritz value and the
    /// orbital (Kohn-Sham) energy functional, being stationary at the exact
    /// eigenvector, carry an error of second order,
    /// |theta - lambda| <= ||r||^2 / delta_gap (Kato-Temple), where delta_gap
    /// is the distance to the rest of the spectrum. The band-energy error is the
    /// occupation- and k-weighted sum of these terms. The density-based Coulomb
    /// estimator does not see this term because it measures input/output
    /// density mismatch only. `c` stands for 1 / delta_gap with the gap scale
    /// set to the atomic unit (1 Ha); it is a units choice, not a fit. Empty
    /// bands (f_i = 0) do not contribute.
    pub fn residual_energy_bound<'a>(
        points: impl IntoIterator<Item = (f64, &'a [f64])>,
        occupations: &[f64],
    ) -> f64 {
        const C_PER_HA: f64 = 1.0;
        C_PER_HA
            * points
                .into_iter()
                .map(|(weight, residuals)| {
                    weight
                        * residuals
                            .iter()
                            .zip(occupations)
                            .filter(|(_, f)| **f > 0.0)
                            .map(|(r, f)| f * r * r)
                            .sum::<f64>()
                })
                .sum::<f64>()
    }

    /// The guarded stop: the ordinary stop decision AND the residual bound
    /// below the class threshold. It can only delay a stop, never cause one.
    pub fn residual_guard(stop: bool, bound_ha: f64, threshold_ha: f64) -> bool {
        stop && bound_ha < threshold_ha
    }
}

/// Post-convergence fixed-point verification of the standard profile (ADR 0043).
/// `Off` (default) keeps the standard path exactly as before; the strict
/// profile ignores it (its certificate always runs).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Verification {
    #[default]
    Off,
    /// One more Kohn-Sham map at the delivered density; pass when the Coulomb
    /// fixed-point residual E_H[F(rho_out) - rho_out] is below the threshold.
    FixedPoint,
}
impl Verification {
    pub fn name(self) -> &'static str {
        match self {
            Self::Off => "off",
            Self::FixedPoint => "fixed_point",
        }
    }
}

/// One verification attempt (ADR 0043, V1 to V6). Energies in Ha.
#[derive(Clone, Debug)]
pub struct VerificationAttempt {
    /// SCF iteration whose output density was verified.
    pub iteration: usize,
    pub passed: bool,
    /// E_H[F(rho_out) - rho_out], the pass criterion (V2).
    pub coulomb_error_ha: f64,
    pub threshold_ha: f64,
    /// Eigen tolerance of the verification solve (the next ratchet value).
    pub eigen_tolerance_ha: f64,
    /// Residual bound of the verification solve when the residual guard is on.
    pub residual_bound_ha: Option<f64>,
    /// Information only (V3): density L2 change and occupied eigenvalue change,
    /// and whether the strict certificate criteria would also have passed.
    pub density_residual: f64,
    pub occupied_eigenvalue_change_ha: f64,
    pub strict_criteria_met: bool,
    pub davidson_iterations: usize,
    pub operator_applications: usize,
    pub wall: std::time::Duration,
}

/// All verification attempts of one run; passed when the last one passed.
#[derive(Clone, Debug, Default)]
pub struct VerificationReport {
    pub attempts: Vec<VerificationAttempt>,
}
impl VerificationReport {
    pub fn passed(&self) -> bool {
        self.attempts.last().is_some_and(|a| a.passed)
    }
    pub fn wall(&self) -> std::time::Duration {
        self.attempts.iter().map(|a| a.wall).sum()
    }
}

/// K maps a signed density to half its Hartree potential times voxel volume.
/// Thus r.K.r = E_H[r] in Ha. G=0 and modes outside the density sphere vanish.
/// All FFT buffers and sphere weights are created once, outside SCF steps.
pub(crate) struct CoulombMetric {
    fft: FftPlan3d,
    source: Array3<Complex64>,
    reciprocal: Array3<Complex64>,
    target: Array3<Complex64>,
    weights: Vec<f64>,
    residual: Vec<f64>,
    mapped: Vec<f64>,
}
impl std::fmt::Debug for CoulombMetric {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CoulombMetric")
            .field("grid", &self.fft.shape())
            .finish()
    }
}
impl CoulombMetric {
    pub(crate) fn new(cell: &Cell, cutoff: f64, dims: [usize; 3]) -> Result<Self, ScfError> {
        let fft = FftPlan3d::new(dims)?;
        let mut weights = vec![0.0; dims.iter().product()];
        let sphere = PlaneWaveBasis::new(cell, Hartree(cutoff))?.density_basis(cell)?;
        let scale = cell.volume_bohr3() / (2.0 * weights.len() as f64);
        for g in sphere.vectors() {
            if g.norm_squared > 0.0 {
                let [i, j, k] = FftGridSpec::map_miller(g.miller, dims)?;
                weights[(i * dims[1] + j) * dims[2] + k] =
                    scale * 4.0 * std::f64::consts::PI / g.norm_squared;
            }
        }
        Ok(Self {
            fft,
            source: Array3::zeros(dims),
            reciprocal: Array3::zeros(dims),
            target: Array3::zeros(dims),
            residual: vec![0.0; weights.len()],
            mapped: vec![0.0; weights.len()],
            weights,
        })
    }
    fn transform(&mut self) -> Result<(), dft_solvers::SolverError> {
        for (z, r) in self.source.iter_mut().zip(&self.residual) {
            *z = (*r).into();
        }
        self.fft
            .to_reciprocal(&self.source, &mut self.reciprocal)
            .map_err(|e| dft_solvers::SolverError::Application(e.to_string()))?;
        for (z, w) in self.reciprocal.iter_mut().zip(&self.weights) {
            *z *= w;
        }
        self.fft
            .to_real(&self.reciprocal, &mut self.target)
            .map_err(|e| dft_solvers::SolverError::Application(e.to_string()))?;
        for (out, z) in self.mapped.iter_mut().zip(&self.target) {
            *out = z.re;
        }
        Ok(())
    }
    pub(crate) fn error(
        &mut self,
        input: &Array3<f64>,
        output: &Array3<f64>,
    ) -> Result<f64, ScfError> {
        for ((r, a), b) in self.residual.iter_mut().zip(input).zip(output) {
            *r = b - a;
        }
        self.transform()?;
        let energy: f64 = self
            .residual
            .iter()
            .zip(&self.mapped)
            .map(|(a, b)| a * b)
            .sum();
        if !energy.is_finite() {
            return Err(ScfError::NonFinite);
        }
        Ok(energy.max(0.0))
    }
}
impl ResidualMetric for CoulombMetric {
    fn duplicate(&self) -> Box<dyn ResidualMetric> {
        Box::new(Self {
            fft: FftPlan3d::new(self.fft.shape()).expect("validated FFT grid"),
            source: self.source.clone(),
            reciprocal: self.reciprocal.clone(),
            target: self.target.clone(),
            weights: self.weights.clone(),
            residual: self.residual.clone(),
            mapped: self.mapped.clone(),
        })
    }
    fn apply(
        &mut self,
        residual: &[f64],
        output: &mut [f64],
    ) -> Result<(), dft_solvers::SolverError> {
        if residual.len() != self.residual.len() || output.len() != residual.len() {
            return Err(dft_solvers::SolverError::InvalidParameter(
                "Coulomb metric grid mismatch",
            ));
        }
        self.residual.copy_from_slice(residual);
        self.transform()?;
        output.copy_from_slice(&self.mapped);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn plateau_record() -> IterationRecord {
        IterationRecord {
            standard_error_ha: Some(1.0),
            standard_redo: false,
            iteration: 1,
            energy: EnergyAudit {
                xc_floor: None,
                nlcc: None,
                kinetic: 0.0,
                local: 0.0,
                local_zero_mode: 0.0,
                nonlocal: 0.0,
                hartree: 0.0,
                xc: 0.0,
                ewald: 0.0,
                band_sum: 0.0,
                double_counting: 0.0,
                band_form: 0.0,
                harris_foulkes: 0.0,
                functional_difference_bound: 0.0,
            },
            energy_change: None,
            density_residual: 1.0,
            max_eigen_residual: 1.0,
            requested_eigen_tolerance: 1e-4,
            davidson_iterations: 2,
            operator_applications: 0,
            davidson_iterations_per_k: vec![],
            operator_applications_per_k: vec![],
            charge_error: 0.0,
            input_charge_error: 0.0,
            min_density: 0.0,
            mixing_fraction: 1.0,
            mixer_reset: false,
            density_contraction: None,
            mixing_history_used: None,
            mixing_resets: 0,
            mixing_raw_recovery: false,
            fixed_point_certificate: None,
            timings: PhaseTimings::default(),
        }
    }

    #[test]
    fn stall_window_progress_floor_and_no_repeated_compounding() {
        let mut history = vec![plateau_record(); 6];
        let tolerance = 1e-4;
        assert_eq!(
            StandardProfile::stall_guard(tolerance, &history[..5]),
            tolerance
        );
        let tightened = StandardProfile::stall_guard(tolerance, &history);
        assert!(tightened < tolerance);
        // The same plateau cannot compound the cut before the window advances.
        assert_eq!(StandardProfile::stall_guard(tightened, &history), tightened);
        assert_eq!(
            StandardProfile::stall_guard(StandardProfile::FLOOR_HA, &history),
            StandardProfile::FLOOR_HA
        );
        history[3].standard_error_ha = Some(0.5);
        assert_eq!(StandardProfile::stall_guard(tolerance, &history), tolerance);
        history[3].standard_error_ha = None;
        assert_eq!(StandardProfile::stall_guard(tolerance, &history), tolerance);
        history[3].standard_error_ha = Some(0.49);
        history[5].standard_error_ha = Some(2.0);
        assert_eq!(StandardProfile::stall_guard(tolerance, &history), tolerance);
    }

    #[test]
    fn density_forcing_never_loosens_and_respects_floor() {
        assert_eq!(StandardProfile::density_forcing(1e-4, &[]), 1e-4);
        let mut r = plateau_record();
        for residual in [0.0, 1e-10, 1e-5, 0.01, 1.0] {
            r.density_residual = residual;
            let t = StandardProfile::density_forcing(1e-4, &[r.clone()]);
            assert!((StandardProfile::FLOOR_HA..=1e-4).contains(&t));
        }
        r.density_residual = 1e-4;
        assert_eq!(StandardProfile::density_forcing(1e-4, &[r]), 1e-8);
    }

    #[test]
    fn signed_coulomb_energy_analytic_and_ratchet() {
        let side: f64 = 8.0;
        let cell =
            Cell::from_user_rows_bohr([[side, 0., 0.], [0., side, 0.], [0., 0., side]]).unwrap();
        let mut metric = CoulombMetric::new(&cell, 0.2, [12, 12, 12]).unwrap();
        let a = 0.03;
        let pi = std::f64::consts::PI;
        let input = Array3::zeros([12, 12, 12]);
        let output = Array3::from_shape_fn([12, 12, 12], |(i, _, _)| {
            0.2 + a * (2.0 * pi * i as f64 / 12.0).cos() + 0.07 * (8.0 * pi * i as f64 / 12.0).cos()
        });
        let expected = pi * side.powi(3) * a * a / (2.0 * pi / side).powi(2);
        let e = metric.error(&input, &output).unwrap();
        assert!((e - expected).abs() < 1e-12, "Ha: {e} vs {expected}");
        assert_eq!(e, metric.error(&output, &input).unwrap());
        assert_eq!(metric.error(&output, &output).unwrap(), 0.0);
        assert_eq!(StandardProfile::ratchet(1e-5, 1.0, 32.0), 1e-5);
        assert_eq!(
            StandardProfile::ratchet(1e-5, 0.0, 32.0),
            StandardProfile::FLOOR_HA
        );
        assert_eq!(StandardProfile::ratchet(0.005, 0.032, 32.0), 0.0001);
    }

    #[test]
    fn residual_energy_bound_is_weighted_sum_of_occupied_squares() {
        // Two k points (weights 0.25, 0.75), bands f = [2, 2, 0]; the empty
        // band's large residual must not count. Units: residual Ha, bound Ha.
        let occ = [2.0, 2.0, 0.0];
        let k0 = [1e-5, 2e-5, 1.0];
        let k1 = [3e-5, 0.0, 1.0];
        let b = StandardProfile::residual_energy_bound([(0.25, &k0[..]), (0.75, &k1[..])], &occ);
        let expected = 0.25 * 2.0 * (1e-10 + 4e-10) + 0.75 * 2.0 * 9e-10;
        assert!((b - expected).abs() <= 1e-24, "{b} vs {expected}");
        assert_eq!(
            StandardProfile::residual_energy_bound(std::iter::empty(), &occ),
            0.0
        );
        // The dev08 case that motivated the guard: 16 doubly occupied bands
        // with residual 3.65e-5 Ha give 4.3e-8 Ha, above the 5e-9 class.
        let r = [3.65e-5; 16];
        let b = StandardProfile::residual_energy_bound([(1.0, &r[..])], &[2.0; 16]);
        assert!((b - 32.0 * 3.65e-5 * 3.65e-5).abs() < 1e-20);
        assert!(!StandardProfile::residual_guard(true, b, 5e-9));
    }

    #[test]
    fn residual_guard_can_only_delay_a_stop() {
        for stop in [false, true] {
            for bound in [0.0, 1e-12, 4.9e-9, 5e-9, 1e-6, f64::INFINITY, f64::NAN] {
                let guarded = StandardProfile::residual_guard(stop, bound, 5e-9);
                // Never stops where the ordinary rule would continue.
                assert!(!guarded || stop, "stop={stop} bound={bound}");
                assert_eq!(guarded, stop && bound < 5e-9);
            }
        }
    }
}
