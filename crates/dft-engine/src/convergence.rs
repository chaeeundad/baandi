//! Additive standard-profile controls and reusable signed Coulomb metric.
use crate::scf::*;
use dft_solvers::mixing::ResidualMetric;

/// Geometry-only hypothesis, not a model of the material's dielectric response.
/// Restore the equal-volume cube's lowest squared frequency in the denominator
/// of the longest represented density mode. The existing Kerker adapter owns
/// FFT normalization, conjugacy on skew/Nyquist grids and the zero-charge mode.
pub(crate) fn cell_kerker_scale(cell: &Cell, basis: &PlaneWaveBasis) -> Option<f64> {
    let g_min2 = basis
        .vectors()
        .iter()
        .map(|g| g.norm_squared)
        .filter(|g2| *g2 > 0.0)
        .fold(f64::INFINITY, f64::min);
    let cube2 = (2.0 * std::f64::consts::PI / cell.volume_bohr3().cbrt()).powi(2);
    let q02 = cube2 - g_min2;
    // Numerical equality guard only (cell inversion / volume / cbrt roundoff).
    (q02 > 64.0 * f64::EPSILON * cube2).then(|| q02.sqrt())
}

pub use crate::model::verify::{
    StandardProfile, Verification, VerificationAttempt, VerificationReport, VerificationSchedule,
    VerificationScheduleEntry,
};

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
    /// The density sphere is the problem's own density list: with a frozen
    /// strain reference (`options.strain_reference`, variable-cell relaxation)
    /// the reference cell's Miller set expressed in `cell`, otherwise the
    /// `cell` sphere at 4 x cutoff (unchanged). ONEPATH-2B: building it from
    /// `cell` alone ignored the frozen reference, so after a cell expansion
    /// the sphere left the frozen FFT grid ("Miller index outside the
    /// canonical frequencies") on the standard path, whose stop rule always
    /// builds this metric.
    pub(crate) fn new(
        cell: &Cell,
        options: &ScfOptions,
        dims: [usize; 3],
    ) -> Result<Self, ScfError> {
        let sphere = options.basis(cell, [0.0; 3], 1.0)?.density_basis(cell)?;
        Self::with_sphere(cell, &sphere, dims)
    }
    /// As `new` with the density sphere supplied by the caller (INIT-REUSE:
    /// problem construction passes its own density list, whose vectors equal
    /// `options.basis(cell, [0; 3], 1.0)?.density_basis(cell)?` element for
    /// element, frozen strain reference included, so the weights are the same
    /// bits without building the sphere again).
    pub(crate) fn with_sphere(
        cell: &Cell,
        sphere: &dft_basis::PlaneWaveBasis,
        dims: [usize; 3],
    ) -> Result<Self, ScfError> {
        let fft = FftPlan3d::new(dims)?;
        let mut weights = vec![0.0; dims.iter().product()];
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
    /// Full-grid weights (row-major; zero outside the density sphere and at G=0).
    pub(crate) fn weights(&self) -> &[f64] {
        &self.weights
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
    /// `error` from the Fourier coefficients of the two densities (ticket
    /// MPI-REMAINING, `BAANDI_STOP_METRIC=kept`): with the forward transform
    /// scaled by 1/N and the unscaled inverse, sum_r d(r) (W d)(r) =
    /// N sum_G w(G) |d(G)|^2 for the real residual d and the symmetric weights
    /// w, so no transform is needed. Rounding class against `error` (the
    /// residual is transformed as two densities; Parseval order of summation).
    pub(crate) fn error_reciprocal(
        &self,
        input: &Array3<Complex64>,
        output: &Array3<Complex64>,
    ) -> Result<f64, ScfError> {
        let (Some(a), Some(b)) = (input.as_slice(), output.as_slice()) else {
            return Err(ScfError::InvalidInput("Coulomb metric layout"));
        };
        if a.len() != self.weights.len() || b.len() != self.weights.len() {
            return Err(ScfError::InvalidInput("Coulomb metric grid mismatch"));
        }
        let sum: f64 = self
            .weights
            .iter()
            .zip(a.iter().zip(b))
            .filter(|(w, _)| **w != 0.0)
            .map(|(w, (a, b))| w * (b - a).norm_sqr())
            .sum();
        let energy = sum * self.weights.len() as f64;
        if !energy.is_finite() {
            return Err(ScfError::NonFinite);
        }
        Ok(energy.max(0.0))
    }
    /// ONEPATH-3 collinear estimator (docs/design/onepath3-spin-metals.md):
    /// E_H[d_n] + E_H[d_m] + E_0[d_m], d_s = output_s - input_s, d_n the charge
    /// and d_m the magnetization residual. E_0 weights the uniform part of d_m
    /// (a change of the total moment, possible under an unconstrained
    /// Fermi-Dirac moment) like the longest wavelength of the density sphere:
    /// E_0 = N_grid * max weight * <d_m>^2 = (Omega/2) (4 pi / G_min^2) <d_m>^2.
    /// An estimator, not an energy bound. Unpolarized limit (d_up = d_down):
    /// exactly `error` of the total densities.
    pub(crate) fn spin_error(
        &mut self,
        input: &[Array3<f64>; 2],
        output: &[Array3<f64>; 2],
    ) -> Result<f64, ScfError> {
        let mut energy = 0.0;
        let mut mean = 0.0;
        for sign in [1.0, -1.0] {
            for (((r, a0), a1), (b0, b1)) in self
                .residual
                .iter_mut()
                .zip(&input[0])
                .zip(&input[1])
                .zip(output[0].iter().zip(&output[1]))
            {
                *r = (b0 - a0) + sign * (b1 - a1);
            }
            if sign < 0.0 {
                mean = self.residual.iter().sum::<f64>() / self.residual.len() as f64;
            }
            self.transform()?;
            energy += self
                .residual
                .iter()
                .zip(&self.mapped)
                .map(|(a, b)| a * b)
                .sum::<f64>()
                .max(0.0);
        }
        let largest = self.weights.iter().copied().fold(0.0, f64::max);
        energy += self.weights.len() as f64 * largest * mean * mean;
        if !energy.is_finite() {
            return Err(ScfError::NonFinite);
        }
        Ok(energy)
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

    #[test]
    fn cell_kerker_scale_is_geometric_and_cubic_is_inert() {
        let scale = |rows| {
            let c = Cell::from_user_rows_bohr(rows).unwrap();
            let b = ScfOptions {
                cutoff_hartree: 1.0,
                ..ScfOptions::default()
            }
            .basis(&c, [0.0; 3], 4.0)
            .unwrap();
            cell_kerker_scale(&c, &b)
        };
        assert_eq!(scale([[8., 0., 0.], [0., 8., 0.], [0., 0., 8.]]), None);
        let q = scale([[8., 0., 0.], [0., 8., 0.], [0., 0., 32.]]).unwrap();
        let expected = ((2.0 * std::f64::consts::PI / 2048.0_f64.cbrt()).powi(2)
            - (2.0 * std::f64::consts::PI / 32.0).powi(2))
        .sqrt();
        assert!((q - expected).abs() < 1e-14, "q0 in Bohr^-1");
        let doubled = scale([[16., 0., 0.], [0., 16., 0.], [0., 0., 64.]]).unwrap();
        assert!((2.0 * doubled - q).abs() < 1e-14);
        // Unimodular re-basing: the shortest represented G need not be a column.
        let rebased = scale([[8., 0., 0.], [8., 8., 0.], [0., 0., 32.]]).unwrap();
        assert!((rebased - q).abs() < 1e-14);
    }

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
            pulay: None,
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
        let mut metric = CoulombMetric::new(
            &cell,
            &ScfOptions {
                cutoff_hartree: 0.2,
                ..ScfOptions::default()
            },
            [12, 12, 12],
        )
        .unwrap();
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

    /// MPI-REMAINING `BAANDI_STOP_METRIC=kept`: the Parseval form from the
    /// forward coefficients of the two densities equals the transform form
    /// to rounding (relative 1e-12) on a skewed cell with a generic density
    /// pair, and to the analytic value on the plane wave above.
    #[test]
    fn coulomb_estimate_from_kept_coefficients_matches_the_transform() {
        use dft_numerics::Fft3d;
        let cell =
            Cell::from_user_rows_bohr([[7.0, 0., 0.], [0.9, 6.4, 0.], [0.4, 0.7, 8.1]]).unwrap();
        let dims = [12, 10, 14];
        let mut metric = CoulombMetric::new(
            &cell,
            &ScfOptions {
                cutoff_hartree: 0.6,
                ..ScfOptions::default()
            },
            dims,
        )
        .unwrap();
        let field = |s: f64| {
            Array3::from_shape_fn(dims, |(i, j, k)| {
                let x = (i * 31 + j * 17 + k * 7) as f64;
                0.05 + 0.01 * ((x * s).sin() + 0.5 * (x * 0.37 * s + 1.0).cos())
            })
        };
        let (input, output) = (field(0.11), field(0.13));
        let mut fft = FftPlan3d::new(dims).unwrap();
        let hat = |n: &Array3<f64>, fft: &mut FftPlan3d| {
            let mut out = Array3::zeros(dims);
            fft.to_reciprocal(&n.mapv(|v| Complex64::new(v, 0.0)), &mut out)
                .unwrap();
            out
        };
        let (a, b) = (hat(&input, &mut fft), hat(&output, &mut fft));
        let transform = metric.error(&input, &output).unwrap();
        let kept = metric.error_reciprocal(&a, &b).unwrap();
        println!("Coulomb estimate transform {transform:e} kept {kept:e} Ha");
        assert!(transform > 0.0);
        assert!((kept - transform).abs() <= 1e-12 * transform);
        assert_eq!(metric.error_reciprocal(&b, &b).unwrap(), 0.0);
    }

    /// ONEPATH-3 collinear estimator: equal channel residuals reduce to the
    /// closed-shell E_H of the total; a pure magnetization wave counts with the
    /// same Coulomb weight; a uniform moment change counts like the longest
    /// wavelength of the sphere (analytic values, Ha).
    #[test]
    fn spin_coulomb_estimator_limits() {
        let side: f64 = 8.0;
        let cell =
            Cell::from_user_rows_bohr([[side, 0., 0.], [0., side, 0.], [0., 0., side]]).unwrap();
        let mut metric = CoulombMetric::new(
            &cell,
            &ScfOptions {
                cutoff_hartree: 0.2,
                ..ScfOptions::default()
            },
            [12, 12, 12],
        )
        .unwrap();
        let pi = std::f64::consts::PI;
        let a = 0.03;
        let wave = |amp: f64, mean: f64| {
            Array3::from_shape_fn([12, 12, 12], |(i, _, _)| {
                mean + amp * (2.0 * pi * i as f64 / 12.0).cos()
            })
        };
        let zero = [Array3::zeros([12, 12, 12]), Array3::zeros([12, 12, 12])];
        let g2 = (2.0 * pi / side).powi(2);
        // Unpolarized limit: d_up = d_down = d/2 gives E_H[d].
        let half = [wave(0.5 * a, 0.1), wave(0.5 * a, 0.1)];
        let total = wave(a, 0.2);
        let closed = metric.error(&Array3::zeros([12, 12, 12]), &total).unwrap();
        let spin = metric.spin_error(&zero, &half).unwrap();
        assert!((spin - closed).abs() <= 1e-14, "{spin:e} vs {closed:e}");
        // Magnetization wave only (d_n = 0): E_H[d_m] = pi Omega a^2 / G^2.
        let split = [wave(0.5 * a, 0.1), wave(-0.5 * a, 0.1)];
        let expected = pi * side.powi(3) * a * a / g2;
        let e = metric.spin_error(&zero, &split).unwrap();
        assert!(
            (e - expected).abs() <= 1e-12 * expected,
            "{e:e} vs {expected:e}"
        );
        // Uniform moment change dM: (Omega/2)(4 pi/G_min^2)(dM/Omega)^2.
        let dm: f64 = 0.01;
        let uniform = [
            wave(0.0, 0.5 * dm / side.powi(3)),
            wave(0.0, -0.5 * dm / side.powi(3)),
        ];
        let expected = 2.0 * pi * dm * dm / (side.powi(3) * g2);
        let e = metric.spin_error(&zero, &uniform).unwrap();
        assert!(
            (e - expected).abs() <= 1e-12 * expected,
            "{e:e} vs {expected:e}"
        );
    }

    #[test]
    fn residual_energy_estimate_is_weighted_sum_of_occupied_squares() {
        // Two k points (weights 0.25, 0.75), bands f = [2, 2, 0]; the empty
        // band's large residual must not count. Units: residual Ha, bound Ha.
        let occ = [2.0, 2.0, 0.0];
        let k0 = [1e-5, 2e-5, 1.0];
        let k1 = [3e-5, 0.0, 1.0];
        let b = StandardProfile::residual_energy_estimate([(0.25, &k0[..]), (0.75, &k1[..])], &occ);
        let expected = 0.25 * 2.0 * (1e-10 + 4e-10) + 0.75 * 2.0 * 9e-10;
        assert!((b - expected).abs() <= 1e-24, "{b} vs {expected}");
        assert_eq!(
            StandardProfile::residual_energy_estimate(std::iter::empty(), &occ),
            0.0
        );
        // The dev08 case that motivated the guard: 16 doubly occupied bands
        // with residual 3.65e-5 Ha give 4.3e-8 Ha, above the 5e-9 class.
        let r = [3.65e-5; 16];
        let b = StandardProfile::residual_energy_estimate([(1.0, &r[..])], &[2.0; 16]);
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
