//! Explicit fixed-sector collinear SCF. Legacy scalar result/session APIs remain separate.
use super::finite_temperature::{OccupationWorkspace, ThermodynamicLedger};
use super::*;
use dft_physics::occupations::Constraint;
use dft_physics::spin::{SpinDensity, SpinFunctional, SpinHartree, SpinXc, SpinXcResponse};

#[derive(Clone, Debug)]
pub struct SpinSector {
    pub n_up: usize,
    pub n_down: usize,
    /// Channel order is up/down along this unit Cartesian axis (no SOC).
    pub axis_cartesian: [f64; 3],
    /// Provenance only; the fixed sector, not this seed, constrains the moment.
    pub initial_moment_seed: Option<f64>,
}
impl SpinSector {
    pub fn validate(&self, electrons: f64) -> Result<(), ScfError> {
        let n = self
            .n_up
            .checked_add(self.n_down)
            .ok_or(ScfError::InvalidInput(
                "spin.constraint: electron count overflow; reduce counts",
            ))?;
        if n == 0 || (n as f64 - electrons).abs() > 1e-10 {
            return Err(ScfError::InvalidInput(
                "spin.constraint: channel counts must sum to valence electrons; correct n_up and n_down",
            ));
        }
        if self.axis_cartesian.iter().any(|v| !v.is_finite())
            || (self.axis_cartesian.iter().map(|v| v * v).sum::<f64>() - 1.0).abs() > 1e-12
        {
            return Err(ScfError::InvalidInput(
                "spin.axis: requires a finite unit Cartesian vector; normalize the axis",
            ));
        }
        if self
            .initial_moment_seed
            .is_some_and(|v| !v.is_finite() || v.abs() > electrons)
        {
            return Err(ScfError::InvalidInput(
                "spin.initial_moment_seed: outside electron capacity; use a finite seed with abs(seed)<=N",
            ));
        }
        Ok(())
    }
    pub fn moment(&self) -> f64 {
        self.n_up as f64 - self.n_down as f64
    }
}

#[derive(Clone, Debug)]
pub struct SpinCertificate {
    pub passed: bool,
    pub charge_density_residual: f64,
    pub spin_density_residual: f64,
    pub energy_change_ha: f64,
    pub occupied_eigenvalue_change_ha: f64,
    pub strict_eigen_residual_ha: [f64; 2],
    pub channel_electrons: [f64; 2],
    pub gap_ha: [Option<f64>; 2],
    pub density_tolerance: f64,
    pub energy_tolerance_ha: f64,
    pub eigen_tolerance_ha: f64,
}
#[derive(Clone, Debug)]
pub struct SpinIteration {
    pub iteration: usize,
    pub energy_ha: f64,
    pub charge_density_residual: f64,
    pub spin_density_residual: f64,
    pub mixing_fraction: f64,
}
pub struct CollinearResult {
    pub finite_certificate: Option<super::finite_temperature::FiniteTCertificate>,
    pub finite_occupations: Option<OccupationWorkspace>,
    pub thermodynamics: Option<ThermodynamicLedger>,
    pub spectral_completeness: Option<crate::SpectralReport>,
    pub converged: bool,
    pub sector: SpinSector,
    pub density: SpinDensity,
    /// Spin-major, then input k-point order. Occupations exclude k weights.
    pub k_points: [Vec<KPointSolution>; 2],
    /// Fixed-policy band-only view. Empty for finite T: use finite_occupations,
    /// whose explicit spin/k/band layout cannot be represented by this field.
    pub occupations: [Vec<f64>; 2],
    pub energy: EnergyAudit,
    pub certificate: Option<SpinCertificate>,
    pub history: Vec<SpinIteration>,
    pub channel_electrons: [f64; 2],
    pub moment_mu_b: f64,
}
impl CollinearResult {
    pub fn legacy_scalar_result(&self) -> Result<&ScfResult, ScfError> {
        Err(ScfError::InvalidInput(
            "result.spin: scalar legacy view is ambiguous; request spin-resolved data",
        ))
    }
}
struct BlockWork {
    preconditioner: DiagonalPreconditioner,
    initial: ComplexMatrix,
    orbitals: Vec<Vec<Complex64>>,
    nonlocal: Vec<Complex64>,
    projections: Vec<Complex64>,
    solution: Option<EigenSolution>,
    saved: Option<EigenSolution>,
}
pub struct CollinearProblem {
    finite: Option<OccupationWorkspace>,
    channels: [ScfProblem; 2],
    sector: SpinSector,
    blocks: [Vec<BlockWork>; 2],
    xc: SpinXc,
    input_xc: SpinXcResponse,
    output_xc: SpinXcResponse,
    input_h: SpinHartree,
    output_h: SpinHartree,
    total: Array3<f64>,
    effective: Array3<f64>,
    mixers: [Box<dyn MixerHistory>; 2],
}
impl CollinearProblem {
    /// `points` must be a full quadrature: no magnetic/time-reversal reduction.
    /// Explicit sector; no search over moments or unconstrained ground-state claim.
    pub fn new(
        cell: &Cell,
        atoms: &[Atom],
        potentials: &[SpeciesPotential<'_>],
        mut options: ScfOptions,
        points: &[KPoint],
        sector: SpinSector,
    ) -> Result<Self, ScfError> {
        let functional = match options.functional {
            Functional::LdaPz => {
                return Err(ScfError::InvalidInput(
                    "xc.spin: required polarized capability unavailable for lda-pz; select pbe or lda-pw",
                ));
            }
            Functional::LdaPw => SpinFunctional::LdaPw,
            Functional::Pbe => SpinFunctional::Pbe,
        };
        let electrons = atoms.iter().try_fold(0.0, |sum, a| {
            potentials
                .iter()
                .find(|p| p.label == a.species)
                .map(|p| sum + p.potential.z_valence())
                .ok_or(ScfError::InvalidInput("unknown species"))
        })?;
        sector.validate(electrons)?;
        let max_occupied = sector.n_up.max(sector.n_down);
        options.bands = Some(
            options.bands.unwrap_or(
                max_occupied
                    .checked_add(options.extra_bands)
                    .ok_or(ScfError::InvalidInput("band count overflow"))?,
            ),
        );
        if options.bands.unwrap() <= max_occupied {
            return Err(ScfError::InvalidInput(
                "spin.constraint: fixed channels require an empty band; increase bands",
            ));
        }
        let up = ScfProblem::construct(
            cell,
            atoms,
            potentials,
            options.clone(),
            points,
            Some(1e-7),
            Some(sector.n_up),
        )?;
        let down = ScfProblem::construct(
            cell,
            atoms,
            potentials,
            options,
            points,
            Some(1e-7),
            Some(sector.n_down),
        )?;
        let dims = up.grid.dims();
        let mut blocks = [
            Vec::with_capacity(points.len()),
            Vec::with_capacity(points.len()),
        ];
        for (s, p) in [&up, &down].into_iter().enumerate() {
            for state in &p.states {
                let bands = p.occupations.len();
                let rows = state.basis.vectors().len();
                blocks[s].push(BlockWork {
                    preconditioner: DiagonalPreconditioner::new(
                        state.hamiltonian.kinetic_diagonal().to_vec(),
                        0.1,
                    )?,
                    initial: free_electron_guess(state.hamiltonian.kinetic_diagonal(), bands)?,
                    orbitals: vec![vec![Complex64::default(); rows]; bands],
                    nonlocal: vec![Complex64::default(); rows],
                    projections: vec![
                        Complex64::default();
                        state.hamiltonian.nonlocal().scratch_len()
                    ],
                    solution: None,
                    saved: None,
                });
            }
        }
        let mixers = [make_spin_mixer(&up)?, make_spin_mixer(&down)?];
        Ok(Self {
            finite: None,
            channels: [up, down],
            sector,
            blocks,
            xc: SpinXc::new(functional, cell, dims)?,
            input_xc: SpinXcResponse::new(dims),
            output_xc: SpinXcResponse::new(dims),
            input_h: SpinHartree::new(cell, dims)?,
            output_h: SpinHartree::new(cell, dims)?,
            total: Array3::zeros(dims),
            effective: Array3::zeros(dims),
            mixers,
        })
    }
    #[allow(clippy::too_many_arguments)]
    pub fn new_fermi_dirac(
        cell: &Cell,
        atoms: &[Atom],
        potentials: &[SpeciesPotential<'_>],
        mut options: ScfOptions,
        points: &[KPoint],
        sector: SpinSector,
        theta_ha: f64,
        constraint: Constraint,
        band_cap: usize,
    ) -> Result<Self, ScfError> {
        super::finite_temperature::validate_options(&options, theta_ha, band_cap)?;
        options.bands = Some(band_cap);
        let mut problem = Self::new(cell, atoms, potentials, options, points, sector)?;
        let finite = OccupationWorkspace::new(band_cap, points, true, theta_ha, constraint)?;
        for p in &mut problem.channels {
            p.occupations.fill(1.0);
        }
        problem.finite = Some(finite);
        Ok(problem)
    }
    pub fn grid_dimensions(&self) -> [usize; 3] {
        self.channels[0].grid.dims()
    }
    fn map(
        &mut self,
        density: &[Array3<f64>; 2],
        output: &mut [Array3<f64>; 2],
        tolerance: f64,
    ) -> Result<EnergyAudit, ScfError> {
        self.total.assign(&density[0]);
        self.total += &density[1];
        self.input_h.evaluate(&self.total)?;
        self.xc.evaluate_with_core(
            density,
            self.channels[0].core_density.as_ref(),
            &mut self.input_xc,
        )?;
        let mut kinetic = 0.0;
        let mut nonlocal = 0.0;
        let mut band_sum = 0.0;
        let mut eig_bound = 0.0;
        for (s, channel_output) in output.iter_mut().enumerate() {
            let p = &mut self.channels[s];
            let bands = p.occupations.len();
            let solver = Davidson {
                bands,
                tolerance,
                max_iterations: p.options.davidson_max_iterations,
                max_subspace: p.options.davidson_subspace.max(6 * bands),
                block_size: bands,
                seed: p.options.seed,
            };
            let policy = BandTolerances {
                strict_bands: (p.occupations.iter().filter(|&&f| f > 0.0).count() + 1).min(bands),
                buffer_tolerance: tolerance.max(1e-5),
            };
            for (((v, &ion), &h), &x) in self
                .effective
                .iter_mut()
                .zip(&p.ionic_real)
                .zip(&self.input_h.potential)
                .zip(&self.input_xc.potentials[s])
            {
                *v = ion + h + x;
            }
            channel_output.fill(0.0);
            for (state, work) in p.states.iter_mut().zip(&mut self.blocks[s]) {
                state.hamiltonian.copy_local_potential(&self.effective)?;
                // Reuse the established solver. Its owned result/workspace lifecycle
                // is unchanged; all new band/grid/audit buffers live in BlockWork.
                let guess = work.solution.as_ref().map(|v| &v.eigenvectors).or_else(|| {
                    (work.saved.is_some() || p.options.free_electron_initial_guess)
                        .then_some(&work.initial)
                });
                let solution = p
                    .options
                    .eigensolver
                    .solve_with_band_tolerances_and_certification(
                        &solver,
                        &state.hamiltonian,
                        guess,
                        Some(&work.preconditioner),
                        policy,
                        tolerance <= p.options.eigen_tolerance,
                    )?;
                for (b, orbital) in work.orbitals.iter_mut().enumerate() {
                    orbital.copy_from_slice(solution.eigenvectors.column(b));
                }
                work.solution = Some(solution);
            }
        }
        if let Some(finite) = &mut self.finite {
            for s in 0..2 {
                for (k, work) in self.blocks[s].iter().enumerate() {
                    let first = (s * self.blocks[s].len() + k) * finite.bands;
                    finite.eigenvalues[first..first + finite.bands]
                        .copy_from_slice(&work.solution.as_ref().unwrap().eigenvalues);
                }
            }
            finite.evaluate(self.channels[0].electrons)?;
        }
        for (s, channel_output) in output.iter_mut().enumerate() {
            let p = &mut self.channels[s];
            let volume = p.cell.volume_bohr3();
            for (k, (state, work)) in p.states.iter_mut().zip(&mut self.blocks[s]).enumerate() {
                let occupations = self
                    .finite
                    .as_ref()
                    .map_or(p.occupations.as_slice(), |f| f.block(s, k));
                let solution = work.solution.as_ref().unwrap();
                state
                    .band_density
                    .compute(&work.orbitals, occupations, &state.grid, volume)?;
                for (out, &n) in channel_output.iter_mut().zip(state.band_density.density()) {
                    *out += state.point.weight * n;
                }
                for (b, (&f, orbital)) in occupations.iter().zip(&work.orbitals).enumerate() {
                    if f == 0.0 {
                        continue;
                    }
                    let weight = f * state.point.weight;
                    kinetic += weight
                        * orbital
                            .iter()
                            .zip(state.hamiltonian.kinetic_diagonal())
                            .map(|(c, t)| c.norm_sqr() * t)
                            .sum::<f64>();
                    state.hamiltonian.nonlocal().apply_into(
                        orbital,
                        &mut work.nonlocal,
                        &mut work.projections,
                    )?;
                    nonlocal += weight
                        * orbital
                            .iter()
                            .zip(&work.nonlocal)
                            .map(|(c, v)| (c.conj() * v).re)
                            .sum::<f64>();
                    band_sum += weight * solution.eigenvalues[b];
                    eig_bound += weight * solution.residuals[b];
                }
            }
        }
        self.total.assign(&output[0]);
        self.total += &output[1];
        self.output_h.evaluate(&self.total)?;
        self.xc.evaluate_with_core(
            output,
            self.channels[0].core_density.as_ref(),
            &mut self.output_xc,
        )?;
        let p = &self.channels[0];
        let weight = p.cell.volume_bohr3() / self.total.len() as f64;
        let local = weight
            * self
                .total
                .iter()
                .zip(&p.ionic_real)
                .map(|(n, v)| n * v)
                .sum::<f64>();
        let dc = weight
            * (0..2)
                .map(|s| {
                    output[s]
                        .iter()
                        .zip(&self.input_h.potential)
                        .zip(&self.input_xc.potentials[s])
                        .map(|((n, h), x)| n * (h + x))
                        .sum::<f64>()
                })
                .sum::<f64>();
        let linear_h = weight
            * (0..2)
                .map(|s| {
                    output[s]
                        .iter()
                        .zip(&density[s])
                        .zip(&self.input_h.potential)
                        .map(|((o, i), h)| (o - i) * h)
                        .sum::<f64>()
                })
                .sum::<f64>();
        let linear_x = weight
            * (0..2)
                .map(|s| {
                    output[s]
                        .iter()
                        .zip(&density[s])
                        .zip(&self.input_xc.potentials[s])
                        .map(|((o, i), x)| (o - i) * x)
                        .sum::<f64>()
                })
                .sum::<f64>();
        let ewald = p.ewald.energy.0;
        let audit = EnergyAudit {
            xc_floor: None,
            kinetic,
            local,
            local_zero_mode: self.total.sum() * weight * p.ionic_g[[0, 0, 0]].re,
            nonlocal,
            hartree: self.output_h.energy_ha,
            xc: self.output_xc.energy_ha,
            nlcc: self.channels[0]
                .core_density
                .as_ref()
                .map(|core| crate::energy::CoreXcLedger {
                    core_electrons: core.sum() * weight,
                    core_potential_integral_ha: self.output_xc.core_potential_integral_ha,
                    valence_potential_integral_ha: self.output_xc.valence_potential_integral_ha,
                }),
            ewald,
            band_sum,
            double_counting: dc,
            band_form: band_sum - dc + self.output_h.energy_ha + self.output_xc.energy_ha + ewald,
            harris_foulkes: band_sum - self.input_h.energy_ha + self.input_xc.energy_ha
                - self.input_xc.valence_potential_integral_ha
                + ewald,
            functional_difference_bound: (self.output_h.energy_ha
                - self.input_h.energy_ha
                - linear_h)
                .abs()
                + (self.output_xc.energy_ha - self.input_xc.energy_ha - linear_x).abs()
                + eig_bound,
        };
        if !audit.total_hartree().is_finite() || !audit.band_form.is_finite() {
            return Err(ScfError::NonFinite);
        }
        Ok(audit)
    }
    fn strict_residuals(&self) -> [f64; 2] {
        std::array::from_fn(|s| {
            self.blocks[s]
                .iter()
                .flat_map(|w| {
                    let a = w.solution.as_ref().unwrap();
                    a.residuals
                        .iter()
                        .zip(&a.buffer)
                        .filter_map(|(&r, &buffer)| (!buffer).then_some(r))
                })
                .fold(0.0, f64::max)
        })
    }
    fn gaps(&self) -> [Option<f64>; 2] {
        std::array::from_fn(|s| {
            let occupied = if s == 0 {
                self.sector.n_up
            } else {
                self.sector.n_down
            };
            if occupied == 0 {
                return None;
            }
            let high = self.blocks[s]
                .iter()
                .map(|w| w.solution.as_ref().unwrap().eigenvalues[occupied - 1])
                .fold(f64::NEG_INFINITY, f64::max);
            let low = self.blocks[s]
                .iter()
                .map(|w| w.solution.as_ref().unwrap().eigenvalues[occupied])
                .fold(f64::INFINITY, f64::min);
            Some(low - high)
        })
    }
    pub fn run(mut self) -> Result<CollinearResult, ScfError> {
        let o = self.channels[0].options.clone();
        let dims = self.grid_dimensions();
        let volume = self.channels[0].cell.volume_bohr3();
        let electrons = self.channels[0].electrons;
        let counts = if let Some(f) = &self.finite {
            match f.constraint() {
                Constraint::FixedMoment(m) => [(electrons + m) / 2.0, (electrons - m) / 2.0],
                Constraint::Total => {
                    let m = self.sector.initial_moment_seed.unwrap_or(0.0);
                    [(electrons + m) / 2.0, (electrons - m) / 2.0]
                }
            }
        } else {
            [self.sector.n_up as f64, self.sector.n_down as f64]
        };
        let weight = volume / dims.iter().product::<usize>() as f64;
        let mut density = std::array::from_fn(|s| Array3::from_elem(dims, counts[s] / volume));
        let mut output = std::array::from_fn(|_| Array3::zeros(dims));
        let mut check = std::array::from_fn(|_| Array3::zeros(dims));
        let mut proposals = std::array::from_fn::<_, 2, _>(|_| Array3::zeros(dims));
        let mut spectral_report =
            (o.spectral_guard == crate::SpectralGuard::Final).then(crate::SpectralReport::default);
        let mut history = Vec::with_capacity(o.max_iterations);
        let mut saved_occupations = self
            .finite
            .as_ref()
            .map(|f| vec![0.0; f.occupations().len()]);
        let mut previous_energy: Option<f64> = None;
        let mut previous_residual: f64 = 0.1;
        let mut best: f64 = f64::INFINITY;
        let mut growing = 0;
        for iteration in 1..=o.max_iterations {
            let mut finite_certificate = None;
            let mut final_certificate = None;
            let tolerance = if o.residual_eigen_schedule {
                (0.1 * previous_residual)
                    .min(1e-2)
                    .max(0.1 * o.eigen_tolerance)
            } else {
                (o.initial_eigen_tolerance * o.eigen_tolerance_decay.powf((iteration - 1) as f64))
                    .max(0.1 * o.eigen_tolerance)
            };
            let energy = self
                .map(&density, &mut output, tolerance)
                .map_err(|error| {
                    ScfError::Numerical(format!("collinear SCF iteration {iteration}: {error}"))
                })?;
            let finite_evaluation = self.finite.as_ref().and_then(|f| f.evaluation);
            let free_energy =
                energy.total_hartree() - finite_evaluation.map_or(0.0, |f| f.entropy_term_ha);
            if finite_evaluation.is_some_and(|f| {
                f.edge_max_fraction > 1e-8 || f.edge_electrons > 1e-10 * electrons.max(1.0)
            }) {
                return Err(ScfError::InvalidInput(
                    "occupations.bands: insufficient empty-state coverage at cap; increase cap",
                ));
            }
            let (charge_residual, spin_residual) =
                spin_distances(&density, &output, weight, electrons);
            let strict = self.strict_residuals();
            let charge: [f64; 2] = std::array::from_fn(|s| output[s].sum() * weight);
            let numbers_ok = charge
                .iter()
                .zip(counts)
                .all(|(a, b)| (a - b).abs() <= 1e-10 * electrons.max(1.0));
            if !charge_residual.is_finite() || !spin_residual.is_finite() {
                return Err(ScfError::NonFinite);
            }
            let residual = charge_residual.max(spin_residual);
            if residual > o.divergence_factor * best.max(o.density_tolerance) {
                growing += 1;
            } else {
                growing = 0;
            }
            best = best.min(residual);
            if growing >= o.divergence_patience {
                return Err(ScfError::Numerical(format!(
                    "collinear SCF iteration {iteration}: spin density residual exceeded configured growth limit; last evaluated energy_ha={:.16e}, electrons={:.16e}",
                    energy.total_hartree(),
                    charge.iter().sum::<f64>()
                )));
            }
            let numbers_ok = if let Some(f) = finite_evaluation {
                f.electron_residual.abs() <= 1e-10 * electrons.max(1.0)
                    && f.moment_residual.unwrap_or(0.0).abs() <= 1e-10 * electrons.max(1.0)
            } else {
                numbers_ok
                    && ((charge[0] + charge[1]) - electrons).abs() <= 1e-10 * electrons.max(1.0)
                    && ((charge[0] - charge[1]) - self.sector.moment()).abs()
                        <= 1e-10 * electrons.max(1.0)
            };
            let candidate = charge_residual <= o.density_tolerance
                && spin_residual <= o.density_tolerance
                && strict.iter().all(|r| *r <= o.eigen_tolerance)
                && numbers_ok
                && previous_energy.is_some_and(|e| (free_energy - e).abs() <= o.energy_tolerance)
                && tolerance <= o.eigen_tolerance;
            history.push(SpinIteration {
                iteration,
                energy_ha: energy.total_hartree(),
                charge_density_residual: charge_residual,
                spin_density_residual: spin_residual,
                mixing_fraction: 1.0,
            });
            let mut passed = false;
            let mut unable_to_verify = false;
            let mut rescued = false;
            if candidate {
                let candidate_gaps = if self.finite.is_some() {
                    [None; 2]
                } else {
                    self.gaps()
                };
                if let (Some(f), Some(saved)) = (&self.finite, &mut saved_occupations) {
                    saved.copy_from_slice(f.occupations());
                }
                // Retain the candidate orbitals while checking a fresh map.
                for channel in &mut self.blocks {
                    for w in channel {
                        w.initial
                            .as_mut_slice()
                            .copy_from_slice(w.solution.as_ref().unwrap().eigenvectors.as_slice());
                        w.saved = w.solution.take();
                    }
                }
                let fresh = self
                    .map(&output, &mut check, 0.1 * o.eigen_tolerance)
                    .map_err(|error| {
                        ScfError::Numerical(format!(
                            "collinear certificate iteration {iteration}: {error}"
                        ))
                    })?;
                let (nr, mr) = spin_distances(&output, &check, weight, electrons);
                let residuals = self.strict_residuals();
                let gaps = if self.finite.is_some() {
                    [None; 2]
                } else {
                    self.gaps()
                };
                let fresh_charge = std::array::from_fn(|s| check[s].sum() * weight);
                let mut eigen_change: f64 = 0.0;
                for s in 0..2 {
                    for w in &self.blocks[s] {
                        for (b, &f) in self.channels[s].occupations.iter().enumerate() {
                            if f > 0.0 {
                                eigen_change = eigen_change.max(
                                    (w.solution.as_ref().unwrap().eigenvalues[b]
                                        - w.saved.as_ref().unwrap().eigenvalues[b])
                                        .abs(),
                                );
                            }
                        }
                    }
                }
                let de = (fresh.total_hartree() - energy.total_hartree()).abs();
                passed = nr <= o.density_tolerance
                    && mr <= o.density_tolerance
                    && de <= o.energy_tolerance
                    && eigen_change <= 1e-7
                    && residuals.iter().all(|r| *r <= o.eigen_tolerance)
                    && fresh_charge
                        .iter()
                        .zip(counts)
                        .all(|(a, b)| (a - b).abs() <= 1e-10 * electrons.max(1.0));
                passed = passed
                    && ((fresh_charge[0] + fresh_charge[1]) - electrons).abs()
                        <= 1e-10 * electrons.max(1.0)
                    && ((fresh_charge[0] - fresh_charge[1]) - self.sector.moment()).abs()
                        <= 1e-10 * electrons.max(1.0);
                if let Some(f) = &self.finite {
                    let fresh_eval = f.evaluation.unwrap();
                    let change = saved_occupations
                        .as_ref()
                        .unwrap()
                        .iter()
                        .zip(f.occupations())
                        .enumerate()
                        .map(|(i, (a, b))| {
                            self.channels[0].states[(i / f.bands) % self.channels[0].states.len()]
                                .point
                                .weight
                                * (a - b).abs()
                        })
                        .sum::<f64>();
                    passed = nr <= o.density_tolerance
                        && mr <= o.density_tolerance
                        && (fresh.total_hartree() - fresh_eval.entropy_term_ha - free_energy).abs()
                            <= o.energy_tolerance
                        && residuals.iter().all(|r| *r <= o.eigen_tolerance)
                        && change <= 1e-8 * electrons.max(1.0)
                        && fresh_eval.electron_residual.abs() <= 1e-10 * electrons.max(1.0)
                        && fresh_eval.moment_residual.unwrap_or(0.0).abs()
                            <= 1e-10 * electrons.max(1.0)
                        && (fresh_charge[0] + fresh_charge[1] - electrons).abs()
                            <= 1e-10 * electrons.max(1.0)
                        && match f.constraint() {
                            Constraint::Total => true,
                            Constraint::FixedMoment(m) => {
                                (fresh_charge[0] - fresh_charge[1] - m).abs()
                                    <= 1e-10 * electrons.max(1.0)
                            }
                        }
                        && fresh_eval.edge_max_fraction <= 1e-8
                        && fresh_eval.edge_electrons <= 1e-10 * electrons.max(1.0);
                    finite_certificate = Some(super::finite_temperature::FiniteTCertificate {
                        passed,
                        charge_density_residual: nr,
                        spin_density_residual: Some(mr),
                        retained_eigen_residual_ha: residuals[0].max(residuals[1]),
                        occupation_change_electrons: change,
                        electron_residual: fresh_eval.electron_residual,
                        moment_residual: fresh_eval.moment_residual,
                        free_energy_change_ha: (fresh.total_hartree()
                            - fresh_eval.entropy_term_ha
                            - free_energy)
                            .abs(),
                        edge_max_fraction: fresh_eval.edge_max_fraction,
                        edge_electrons: fresh_eval.edge_electrons,
                    });
                }
                if passed {
                    for gap in gaps.into_iter().chain(candidate_gaps).flatten() {
                        let required = 1e-7_f64.max(2.0 * o.eigen_tolerance);
                        if gap <= required {
                            return Err(ScfError::NonInsulating {
                                gap_hartree: gap,
                                required_hartree: required,
                            });
                        }
                    }
                }
                final_certificate = Some(SpinCertificate {
                    passed,
                    charge_density_residual: nr,
                    spin_density_residual: mr,
                    energy_change_ha: de,
                    occupied_eigenvalue_change_ha: eigen_change,
                    strict_eigen_residual_ha: residuals,
                    channel_electrons: fresh_charge,
                    gap_ha: gaps,
                    density_tolerance: o.density_tolerance,
                    energy_tolerance_ha: o.energy_tolerance,
                    eigen_tolerance_ha: o.eigen_tolerance,
                });
                if passed && let Some(report) = &mut spectral_report {
                    (rescued, unable_to_verify) = self.audit_reported_spectrum(
                        &density,
                        iteration,
                        iteration < o.max_iterations,
                        report,
                    )?;
                    passed = report.passed;
                }
                if let (Some(f), Some(saved)) = (&mut self.finite, &saved_occupations) {
                    f.restore(saved, finite_evaluation.unwrap());
                }
                // The ordinary candidate remains the reported state. Repaired
                // starts may enter a subsequent ordinary map, never this result.
                for channel in &mut self.blocks {
                    for w in channel {
                        w.solution = w.saved.take();
                    }
                }
            }

            if passed || unable_to_verify || iteration == o.max_iterations {
                let k_points = std::array::from_fn(|s| {
                    self.blocks[s]
                        .iter_mut()
                        .zip(&self.channels[s].states)
                        .map(|(w, k)| KPointSolution {
                            point: k.point,
                            solution: w.solution.take().unwrap(),
                        })
                        .collect()
                });
                return Ok(CollinearResult {
                    finite_certificate,
                    thermodynamics: finite_evaluation.map(|f| {
                        ThermodynamicLedger::new(energy.total_hartree(), f.entropy_term_ha)
                    }),
                    finite_occupations: self.finite,
                    spectral_completeness: spectral_report.map(|mut r| {
                        r.passed &= passed;
                        r
                    }),
                    converged: passed,
                    sector: self.sector,
                    density: SpinDensity::Collinear {
                        grid_id: "density".into(),
                        channels: output,
                    },
                    k_points,
                    occupations: if finite_evaluation.is_some() {
                        [Vec::new(), Vec::new()]
                    } else {
                        [
                            self.channels[0].occupations.clone(),
                            self.channels[1].occupations.clone(),
                        ]
                    },
                    energy,
                    certificate: if finite_evaluation.is_some() {
                        None
                    } else {
                        final_certificate
                    },
                    history,
                    channel_electrons: charge,
                    moment_mu_b: charge[0] - charge[1],
                });
            }
            if rescued {
                for mixer in &mut self.mixers {
                    mixer.reset();
                }
                best = f64::INFINITY;
                growing = 0;
            }
            for s in 0..2 {
                if counts[s] == 0.0 && self.finite.is_none() {
                    proposals[s].fill(0.0);
                    continue;
                }
                // Existing mixer owns its return/history allocations; no extra
                // channel buffer is allocated by this adapter.
                let proposal = self.mixers[s].mix(
                    density[s].as_slice().unwrap(),
                    output[s].as_slice().unwrap(),
                    if self.finite.is_some() {
                        charge[s] / volume
                    } else {
                        counts[s] / volume
                    },
                )?;
                proposals[s]
                    .as_slice_mut()
                    .unwrap()
                    .copy_from_slice(&proposal);
            }
            let mut fraction: f64 = 1.0;
            for s in 0..2 {
                for (&a, &b) in density[s].iter().zip(&proposals[s]) {
                    if b < 0.0 {
                        fraction = fraction.min(0.9 * a / (a - b));
                    }
                }
            }
            if !fraction.is_finite() || fraction <= 0.0 {
                return Err(ScfError::Numerical(
                    "spin mixing has no feasible common step".into(),
                ));
            }
            for s in 0..2 {
                for (a, &b) in density[s].iter_mut().zip(&proposals[s]) {
                    *a += fraction * (b - *a);
                }
            }
            if fraction < 1.0 {
                for mixer in &mut self.mixers {
                    mixer.reset();
                }
            }
            history.last_mut().unwrap().mixing_fraction = fraction;
            previous_energy = Some(free_energy);
            previous_residual = charge_residual.max(spin_residual);
        }
        unreachable!()
    }
}
fn spin_distances(a: &[Array3<f64>; 2], b: &[Array3<f64>; 2], weight: f64, n: f64) -> (f64, f64) {
    let mut charge = 0.0;
    let mut spin = 0.0;
    for (((a0, a1), b0), b1) in a[0].iter().zip(&a[1]).zip(&b[0]).zip(&b[1]) {
        let u = b0 - a0;
        let d = b1 - a1;
        charge += (u + d).powi(2);
        spin += (u - d).powi(2);
    }
    ((weight * charge).sqrt() / n, (weight * spin).sqrt() / n)
}
pub(super) fn make_spin_mixer(p: &ScfProblem) -> Result<Box<dyn MixerHistory>, ScfError> {
    Ok(match p.options.mixing {
        Mixing::Linear { beta } => Box::new(LinearMixer { beta }),
        Mixing::Pulay { beta, history } => Box::new(PulayMixer::new(beta, history)?),
        Mixing::KerkerPulay {
            beta,
            history,
            q0_inverse_bohr,
        } => {
            let dims = p.grid.dims();
            let mut fft = FftPlan3d::new(dims)?;
            let mut source = Array3::zeros(dims);
            let mut target = Array3::zeros(dims);
            let transform = move |values: &mut [Complex64], inverse: bool| {
                source.as_slice_mut().unwrap().copy_from_slice(values);
                if inverse {
                    fft.to_real(&source, &mut target)
                } else {
                    fft.to_reciprocal(&source, &mut target)
                }
                .map_err(|e| dft_solvers::SolverError::Application(e.to_string()))?;
                values.copy_from_slice(target.as_slice().unwrap());
                Ok(())
            };
            Box::new(PreconditionedPulay::new(
                beta,
                history,
                dims.iter().product(),
                Kerker::new(
                    dims,
                    p.cell.reciprocal_columns_inverse_bohr(),
                    q0_inverse_bohr,
                    transform,
                )?,
            )?)
        }
    })
}

impl CollinearProblem {
    fn audit_reported_spectrum(
        &mut self,
        input: &[Array3<f64>; 2],
        iteration: usize,
        allow_rescue: bool,
        report: &mut crate::SpectralReport,
    ) -> Result<(bool, bool), ScfError> {
        let start = Instant::now();
        self.total.assign(&input[0]);
        self.total += &input[1];
        self.input_h.evaluate(&self.total)?;
        self.xc.evaluate_with_core(
            input,
            self.channels[0].core_density.as_ref(),
            &mut self.input_xc,
        )?;
        let first = report.events.len();
        let mut failed = false;
        for s in 0..2 {
            let p = &mut self.channels[s];
            for (((v, ion), h), x) in self
                .effective
                .iter_mut()
                .zip(&p.ionic_real)
                .zip(&self.input_h.potential)
                .zip(&self.input_xc.potentials[s])
            {
                *v = ion + h + x;
            }
            let bands = p.occupations.len();
            let solver = Davidson {
                bands,
                tolerance: 1e-10,
                max_iterations: p.options.davidson_max_iterations,
                max_subspace: p.options.davidson_subspace.max(6 * bands),
                block_size: bands,
                seed: p.options.seed,
            };
            for (k, (state, w)) in p.states.iter_mut().zip(&mut self.blocks[s]).enumerate() {
                state.hamiltonian.copy_local_potential(&self.effective)?;
                match report.audit(
                    iteration,
                    Some(s),
                    k,
                    &solver,
                    &state.hamiltonian,
                    &w.preconditioner,
                    w.saved.as_ref().unwrap(),
                ) {
                    Some(checked) => {
                        if checked.spectral_completeness.as_ref().unwrap().retries == 0 {
                            w.saved.as_mut().unwrap().spectral_completeness =
                                checked.spectral_completeness;
                        } else if allow_rescue {
                            w.saved = Some(checked);
                        }
                    }
                    None => failed = true,
                }
            }
        }
        report.passed = report.events[first..].iter().all(|e| e.passed);
        report.elapsed += start.elapsed();
        Ok((report.events[first..].iter().any(|e| e.retries > 0), failed))
    }
}
