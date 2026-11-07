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
}
