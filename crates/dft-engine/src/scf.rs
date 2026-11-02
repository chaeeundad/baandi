use crate::energy::{self, AuditInput, EnergyAudit};
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_geometry::{Atom, Cell};
use dft_physics::{
    data::*,
    density,
    ewald::{self, EwaldOptions, Ion},
    hamiltonian::{self, Hamiltonian, NonlocalOperator, Species},
    hartree, xc,
};
use dft_solvers::{
    Davidson, DiagonalPreconditioner, EigenSolution,
    mixing::{DensityMixer, LinearMixer, PulayMixer},
};
use dft_upf::Pseudopotential;
use std::{
    fmt,
    time::{Duration, Instant},
};

#[derive(Debug)]
pub enum ScfError {
    InvalidInput(&'static str),
    OddElectronCount(f64),
    InsufficientBands {
        bands: usize,
        occupied: usize,
    },
    Numerical(String),
    NonFinite,
    NotConverged {
        history: Vec<IterationRecord>,
    },
    Diverged {
        reason: String,
        history: Vec<IterationRecord>,
    },
}
impl fmt::Display for ScfError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}
impl std::error::Error for ScfError {}
macro_rules! numerical_error { ($($t:ty),*) => { $(impl From<$t> for ScfError { fn from(e:$t)->Self { Self::Numerical(e.to_string()) } })* }; }
numerical_error!(dft_basis::BasisError, FftError);
impl From<dft_physics::PhysicsError> for ScfError {
    fn from(e: dft_physics::PhysicsError) -> Self {
        match e {
            dft_physics::PhysicsError::Unrepresentable => Self::NonFinite,
            dft_physics::PhysicsError::Density(n) if !n.is_finite() => Self::NonFinite,
            _ => Self::Numerical(e.to_string()),
        }
    }
}
impl From<hamiltonian::HamiltonianError> for ScfError {
    fn from(e: hamiltonian::HamiltonianError) -> Self {
        match e {
            hamiltonian::HamiltonianError::NonFinite => Self::NonFinite,
            _ => Self::Numerical(e.to_string()),
        }
    }
}
impl From<dft_solvers::SolverError> for ScfError {
    fn from(e: dft_solvers::SolverError) -> Self {
        match e {
            dft_solvers::SolverError::Linalg(LinalgError::NonFinite) => Self::NonFinite,
            _ => Self::Numerical(e.to_string()),
        }
    }
}

#[derive(Clone, Debug)]
pub enum Mixing {
    Linear { beta: f64 },
    Pulay { beta: f64, history: usize },
}
impl Default for Mixing {
    fn default() -> Self {
        Self::Pulay {
            beta: 0.3,
            history: 6,
        }
    }
}
#[derive(Clone, Debug)]
pub struct ScfOptions {
    pub cutoff_hartree: f64,
    pub bands: Option<usize>,
    pub extra_bands: usize,
    pub grid_override: Option<[usize; 3]>,
    pub mixing: Mixing,
    pub density_tolerance: f64,
    pub energy_tolerance: f64,
    pub eigen_tolerance: f64,
    /// tau_k = max(0.1*eigen_tolerance, initial_eigen_tolerance*decay^(k-1)).
    pub initial_eigen_tolerance: f64,
    pub eigen_tolerance_decay: f64,
    pub max_iterations: usize,
    pub davidson_max_iterations: usize,
    pub davidson_subspace: usize,
    pub seed: u64,
    pub divergence_factor: f64,
    pub divergence_patience: usize,
}
impl Default for ScfOptions {
    fn default() -> Self {
        Self {
            cutoff_hartree: 5.0,
            bands: None,
            extra_bands: 4,
            grid_override: None,
            mixing: Mixing::default(),
            density_tolerance: 1e-7,
            energy_tolerance: 1e-8,
            eigen_tolerance: 1e-8,
            initial_eigen_tolerance: 1e-5,
            eigen_tolerance_decay: 0.5,
            max_iterations: 200,
            davidson_max_iterations: 400,
            davidson_subspace: 48,
            seed: 2026,
            divergence_factor: 20.0,
            divergence_patience: 3,
        }
    }
}
#[derive(Clone, Debug, Default)]
pub struct PhaseTimings {
    pub potential: Duration,
    pub diagonalization: Duration,
    pub density: Duration,
    pub energy: Duration,
    pub mixing: Duration,
}
#[derive(Clone, Debug)]
pub struct IterationRecord {
    pub iteration: usize,
    pub energy: EnergyAudit,
    pub energy_change: Option<f64>,
    pub density_residual: f64,
    pub max_eigen_residual: f64,
    pub requested_eigen_tolerance: f64,
    pub davidson_iterations: usize,
    pub charge_error: f64,
    pub input_charge_error: f64,
    pub min_density: f64,
    /// Pulay feasibility damping; one means the full proposal was accepted.
    pub mixing_fraction: f64,
    pub mixer_reset: bool,
    pub fixed_point_certificate: Option<FixedPointCertificate>,
    pub timings: PhaseTimings,
}
#[derive(Clone, Debug)]
pub struct FixedPointCertificate {
    pub density_residual: f64,
    pub occupied_eigenvalue_change: f64,
    pub davidson_iterations: usize,
    pub timings: PhaseTimings,
}
#[derive(Debug)]
pub struct ScfResult {
    pub energy: EnergyAudit,
    pub eigenvalues: Vec<f64>,
    pub occupations: Vec<f64>,
    pub orbitals: ComplexMatrix,
    pub density: Array3<f64>,
    pub charge: f64,
    pub history: Vec<IterationRecord>,
    pub ewald_breakdown: ewald::EwaldEnergy,
    pub setup_time: Duration,
    pub total_time: Duration,
}
impl ScfResult {
    pub fn total_hartree(&self) -> f64 {
        self.energy.total_hartree()
    }
    pub fn total_rydberg(&self) -> f64 {
        self.energy.total_rydberg()
    }
}
/// Species labels resolve explicitly against Atom::species; the PP must declare SLA PZ.
pub struct SpeciesPotential<'a> {
    pub label: &'a str,
    pub potential: &'a Pseudopotential,
}

/// Cached ionic terms, Hamiltonian application plan, and density FFT plan.
/// `step` is exposed for independent post-SCF fixed-point verification.
pub struct ScfProblem {
    cell: Cell,
    basis: PlaneWaveBasis,
    grid: FftGridSpec,
    ionic_g: Array3<Complex64>,
    ionic_real: Array3<f64>,
    hamiltonian: Hamiltonian,
    fft: FftPlan3d,
    ewald: ewald::EwaldEnergy,
    occupations: Vec<f64>,
    electrons: f64,
    options: ScfOptions,
    setup_time: Duration,
}
pub struct ScfStep {
    pub solution: EigenSolution,
    pub density: Array3<f64>,
    pub energy: EnergyAudit,
    pub timings: PhaseTimings,
}
impl ScfProblem {
    pub fn new(
        cell: &Cell,
        atoms: &[Atom],
        potentials: &[SpeciesPotential<'_>],
        options: ScfOptions,
    ) -> Result<Self, ScfError> {
        let start = Instant::now();
        if atoms.is_empty() || potentials.is_empty() {
            return Err(ScfError::InvalidInput("empty system"));
        }
        for value in [
            options.cutoff_hartree,
            options.density_tolerance,
            options.energy_tolerance,
            options.eigen_tolerance,
            options.initial_eigen_tolerance,
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(ScfError::InvalidInput(
                    "positive finite cutoffs/tolerances required",
                ));
            }
        }
        if options.max_iterations == 0
            || options.davidson_max_iterations == 0
            || options.divergence_patience == 0
            || !options.divergence_factor.is_finite()
            || options.divergence_factor <= 1.0
            || !options.eigen_tolerance_decay.is_finite()
            || options.eigen_tolerance_decay <= 0.0
            || options.eigen_tolerance_decay > 1.0
        {
            return Err(ScfError::InvalidInput(
                "iteration/schedule/divergence configuration",
            ));
        }
        let mut positions = vec![Vec::new(); potentials.len()];
        let mut ions = Vec::new();
        for (i, s) in potentials.iter().enumerate() {
            if potentials[..i].iter().any(|p| p.label == s.label) {
                return Err(ScfError::InvalidInput("duplicate species"));
            }
            if s.potential
                .functional()
                .split_whitespace()
                .collect::<Vec<_>>()
                != ["SLA", "PZ"]
            {
                return Err(ScfError::InvalidInput(
                    "only explicitly labelled SLA PZ potentials are supported",
                ));
            }
        }
        for atom in atoms {
            let index = potentials
                .iter()
                .position(|p| p.label == atom.species)
                .ok_or(ScfError::InvalidInput("unknown species"))?;
            if atom.fractional_position.iter().any(|r| !r.is_finite()) {
                return Err(ScfError::InvalidInput("nonfinite atomic position"));
            }
            positions[index].push(
                cell.fractional_to_cartesian_bohr(
                    atom.fractional_position.map(|r| r.rem_euclid(1.0)),
                ),
            );
            ions.push(Ion {
                charge: potentials[index].potential.z_valence(),
                fractional_position: atom.fractional_position,
            });
        }
        let electrons = ions.iter().map(|i| i.charge).sum::<f64>();
        if !electrons.is_finite()
            || electrons <= 0.0
            || (electrons - electrons.round()).abs() > 1e-10
            || electrons > usize::MAX as f64
        {
            return Err(ScfError::InvalidInput(
                "positive integral valence count required",
            ));
        }
        let count = electrons.round() as usize;
        if !count.is_multiple_of(2) {
            return Err(ScfError::OddElectronCount(electrons));
        }
        let occupied = count / 2;
        let bands = options.bands.unwrap_or(
            occupied
                .checked_add(options.extra_bands)
                .ok_or(ScfError::InvalidInput("band count overflow"))?,
        );
        if bands < occupied {
            return Err(ScfError::InsufficientBands { bands, occupied });
        }
        let basis = PlaneWaveBasis::new(cell, Hartree(options.cutoff_hartree))?;
        if bands > basis.vectors().len()
            || (options.davidson_subspace.min(basis.vectors().len()) <= bands
                && bands != basis.vectors().len())
        {
            return Err(ScfError::InvalidInput(
                "bands/subspace incompatible with basis",
            ));
        }
        let density_basis = basis.density_basis(cell)?;
        let grid = FftGridSpec::new(&basis, &density_basis, options.grid_override)?;
        let dg = FftGridSpec::new(&density_basis, &density_basis, Some(grid.dims()))?;
        let species: Vec<_> = potentials
            .iter()
            .zip(&positions)
            .map(|(p, r)| Species {
                potential: p.potential,
                positions_bohr: r,
            })
            .collect();
        let ionic_g = dg.scatter(&hamiltonian::ionic_local(
            &density_basis,
            cell.volume_bohr3(),
            &species,
        )?)?;
        let mut fft = FftPlan3d::new(grid.dims())?;
        let mut ionic = Array3::zeros(grid.dims());
        fft.to_real(&ionic_g, &mut ionic)?;
        let scale = ionic.iter().map(|z| z.re.abs()).fold(1.0, f64::max);
        if ionic.iter().any(|z| z.im.abs() > 1e-12 * scale) {
            return Err(ScfError::Numerical("ionic imaginary potential".into()));
        }
        let ionic_real = ionic.mapv(|z| z.re);
        let hamiltonian = Hamiltonian::new(
            cell,
            &basis,
            ionic_real.clone(),
            NonlocalOperator::new(&basis, cell.volume_bohr3(), &species)?,
        )?;
        let ewald = ewald::energy(cell, &ions, EwaldOptions::default())?;
        let occupations = (0..bands)
            .map(|b| if b < occupied { 2.0 } else { 0.0 })
            .collect();
        Ok(Self {
            cell: cell.clone(),
            basis,
            grid,
            ionic_g,
            ionic_real,
            hamiltonian,
            fft,
            ewald,
            occupations,
            electrons,
            options,
            setup_time: start.elapsed(),
        })
    }
    pub fn grid_dimensions(&self) -> [usize; 3] {
        self.grid.dims()
    }
    pub fn electron_count(&self) -> f64 {
        self.electrons
    }
    pub fn uniform_density(&self) -> Array3<f64> {
        Array3::from_elem(self.grid.dims(), self.electrons / self.cell.volume_bohr3())
    }
    pub fn step(
        &mut self,
        n: &Array3<f64>,
        initial: Option<&ComplexMatrix>,
        tolerance: f64,
    ) -> Result<ScfStep, ScfError> {
        if n.shape() != self.grid.dims() {
            return Err(ScfError::InvalidInput("input density grid"));
        }
        let start = Instant::now();
        let h = hartree::evaluate(n, &self.cell, &mut self.fft)?;
        let x = xc::evaluate(n, self.cell.volume_bohr3())?;
        let mut effective = self.ionic_real.clone();
        for ((v, h), x) in effective.iter_mut().zip(&h.potential).zip(&x.potential) {
            *v += h.re + x;
        }
        self.hamiltonian.set_local_potential(effective)?;
        let mut timings = PhaseTimings {
            potential: start.elapsed(),
            ..Default::default()
        };
        let start = Instant::now();
        let solver = Davidson {
            bands: self.occupations.len(),
            tolerance,
            max_iterations: self.options.davidson_max_iterations,
            max_subspace: self.options.davidson_subspace,
            block_size: self.occupations.len(),
            seed: self.options.seed,
        };
        let preconditioner =
            DiagonalPreconditioner::new(self.hamiltonian.kinetic_diagonal().to_vec(), 0.1)?;
        let solution = solver.solve(&self.hamiltonian, initial, Some(&preconditioner))?;
        timings.diagonalization = start.elapsed();
        let start = Instant::now();
        let orbitals: Vec<_> = (0..solution.eigenvectors.cols())
            .map(|b| solution.eigenvectors.column(b).to_vec())
            .collect();
        let output = density::from_orbitals(
            &orbitals,
            &self.occupations,
            &self.grid,
            self.cell.volume_bohr3(),
            &mut self.fft,
        )?;
        timings.density = start.elapsed();
        let start = Instant::now();
        let energy = energy::evaluate(
            AuditInput {
                cell: &self.cell,
                basis: &self.basis,
                grid: &self.grid,
                ionic_g: &self.ionic_g,
                nonlocal: self.hamiltonian.nonlocal(),
                orbitals: &orbitals,
                occupations: &self.occupations,
                eigenvalues: &solution.eigenvalues,
                eigen_residuals: &solution.residuals,
                input_density: n,
                output_density: &output,
                input_hartree: &h,
                input_xc: &x,
                ewald: self.ewald.energy.0,
            },
            &mut self.fft,
        )?;
        timings.energy = start.elapsed();
        Ok(ScfStep {
            solution,
            density: output,
            energy,
            timings,
        })
    }
    pub fn run(mut self) -> Result<ScfResult, ScfError> {
        let start = Instant::now();
        let mut mixer: Box<dyn DensityMixer> = match self.options.mixing {
            Mixing::Linear { beta } => Box::new(LinearMixer { beta }),
            Mixing::Pulay { beta, history } => Box::new(PulayMixer::new(beta, history)?),
        };
        // Validate mixer parameters even if a future system converges immediately.
        let mut n = self.uniform_density();
        mixer.mix(
            n.as_slice().unwrap(),
            n.as_slice().unwrap(),
            self.electrons / self.cell.volume_bohr3(),
        )?;
        mixer.reset();
        let mut history: Vec<IterationRecord> = Vec::new();
        let mut previous = None;
        let weight = self.cell.volume_bohr3() / n.len() as f64;
        let mut best = f64::INFINITY;
        let mut growing = 0;
        for iteration in 1..=self.options.max_iterations {
            if n.iter()
                .any(|n| !n.is_finite() || *n < -xc::NEGATIVE_ROUNDOFF)
            {
                return Err(ScfError::Diverged {
                    reason: "mixed density is nonfinite or negative beyond 1e-14 Bohr^-3".into(),
                    history,
                });
            }
            let tolerance = (self.options.initial_eigen_tolerance
                * self
                    .options
                    .eigen_tolerance_decay
                    .powf((iteration - 1) as f64))
            .max(0.1 * self.options.eigen_tolerance);
            let step = self
                .step(&n, previous.as_ref(), tolerance)
                .map_err(|e| match e {
                    ScfError::NonFinite => ScfError::Diverged {
                        reason: "nonfinite SCF arithmetic".into(),
                        history: history.clone(),
                    },
                    other => other,
                })?;
            let residual =
                density_distance(&step.density, &n, self.cell.volume_bohr3(), self.electrons);
            let charge = step.density.sum() * weight;
            let charge_error = (charge - self.electrons).abs() / self.electrons;
            let input_charge_error = (n.sum() * weight - self.electrons).abs() / self.electrons;
            let energy_change = history
                .last()
                .map(|r| step.energy.total_hartree() - r.energy.total_hartree());
            let max_eigen_residual = step
                .solution
                .residuals
                .iter()
                .zip(&self.occupations)
                .filter(|(_, f)| **f > 0.0)
                .map(|(r, _)| *r)
                .fold(0.0, f64::max);
            let record = IterationRecord {
                fixed_point_certificate: None,
                mixing_fraction: 1.0,
                mixer_reset: false,
                iteration,
                energy: step.energy.clone(),
                energy_change,
                density_residual: residual,
                max_eigen_residual,
                requested_eigen_tolerance: tolerance,
                davidson_iterations: step.solution.iterations,
                charge_error,
                input_charge_error,
                min_density: step.density.iter().copied().fold(f64::INFINITY, f64::min),
                timings: step.timings,
            };
            let converged = converged(&record, &self.options);
            history.push(record);
            if !residual.is_finite()
                || !step.energy.total_hartree().is_finite()
                || !charge_error.is_finite()
                || !max_eigen_residual.is_finite()
            {
                return Err(ScfError::Diverged {
                    reason: "nonfinite SCF observable".into(),
                    history,
                });
            }
            if residual > self.options.divergence_factor * best.max(self.options.density_tolerance)
            {
                growing += 1;
            } else {
                growing = 0;
            }
            best = best.min(residual);
            if growing >= self.options.divergence_patience {
                return Err(ScfError::Diverged {
                    reason: "density residual exceeded growth limit".into(),
                    history,
                });
            }
            // Certify the returned output density itself, not just the input
            // fixed-point residual. A small-gap response can amplify that residual.
            let certified = if converged {
                let check = self
                    .step(
                        &step.density,
                        Some(&step.solution.eigenvectors),
                        0.1 * self.options.eigen_tolerance,
                    )
                    .map_err(|e| match e {
                        ScfError::NonFinite => ScfError::Diverged {
                            reason: "nonfinite certification arithmetic".into(),
                            history: history.clone(),
                        },
                        other => other,
                    })?;
                let density_residual = density_distance(
                    &check.density,
                    &step.density,
                    self.cell.volume_bohr3(),
                    self.electrons,
                );
                let occupied_eigenvalue_change = check
                    .solution
                    .eigenvalues
                    .iter()
                    .zip(&step.solution.eigenvalues)
                    .zip(&self.occupations)
                    .filter(|(_, f)| **f > 0.0)
                    .map(|((a, b), _)| (a - b).abs())
                    .fold(0.0, f64::max);
                let passed = density_residual <= self.options.density_tolerance
                    && occupied_eigenvalue_change <= 1e-7;
                history.last_mut().unwrap().fixed_point_certificate = Some(FixedPointCertificate {
                    density_residual,
                    occupied_eigenvalue_change,
                    davidson_iterations: check.solution.iterations,
                    timings: check.timings,
                });
                passed
            } else {
                false
            };
            if certified {
                return Ok(ScfResult {
                    energy: step.energy,
                    eigenvalues: step.solution.eigenvalues,
                    occupations: self.occupations,
                    orbitals: step.solution.eigenvectors,
                    density: step.density,
                    charge,
                    history,
                    ewald_breakdown: self.ewald,
                    setup_time: self.setup_time,
                    total_time: start.elapsed() + self.setup_time,
                });
            }
            let mix_start = Instant::now();
            let mut mixed = mixer
                .mix(
                    n.as_slice().unwrap(),
                    step.density.as_slice().unwrap(),
                    self.electrons / self.cell.volume_bohr3(),
                )
                .map_err(|e| ScfError::Diverged {
                    reason: e.to_string(),
                    history: history.clone(),
                })?;
            // Anderson extrapolation is unconstrained. A convex step toward its
            // proposal preserves charge and Fourier support while keeping XC in
            // its density domain. Never clip individual spatial samples.
            if matches!(self.options.mixing, Mixing::Pulay { .. }) {
                let fraction = n
                    .iter()
                    .zip(&mixed)
                    .filter(|(_, b)| **b < 0.0)
                    .map(|(a, b)| 0.9 * a / (a - b))
                    .fold(1.0, f64::min);
                if fraction < 1.0 {
                    for (b, a) in mixed.iter_mut().zip(&n) {
                        *b = a + fraction * (*b - a);
                    }
                    mixer.reset();
                    let record = history.last_mut().unwrap();
                    record.mixing_fraction = fraction;
                    record.mixer_reset = true;
                }
            }
            n.as_slice_mut().unwrap().copy_from_slice(&mixed);
            history.last_mut().unwrap().timings.mixing = mix_start.elapsed();
            previous = Some(step.solution.eigenvectors);
        }
        Err(ScfError::NotConverged { history })
    }
}
pub fn density_distance(a: &Array3<f64>, b: &Array3<f64>, volume: f64, electrons: f64) -> f64 {
    if a.shape() != b.shape()
        || a.is_empty()
        || !volume.is_finite()
        || volume <= 0.0
        || !electrons.is_finite()
        || electrons <= 0.0
    {
        return f64::NAN;
    }
    (volume / a.len() as f64 * a.iter().zip(b).map(|(a, b)| (a - b).powi(2)).sum::<f64>()).sqrt()
        / electrons
}
fn converged(r: &IterationRecord, o: &ScfOptions) -> bool {
    r.density_residual <= o.density_tolerance
        && r.energy_change
            .is_some_and(|de| de.abs() <= o.energy_tolerance)
        && r.max_eigen_residual <= o.eigen_tolerance
        && r.requested_eigen_tolerance <= o.eigen_tolerance
        && r.charge_error <= 1e-8
        && r.input_charge_error <= 1e-8
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn convergence_requires_every_certificate_on_same_iteration() {
        let o = ScfOptions::default();
        let energy = EnergyAudit {
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
        };
        let good = IterationRecord {
            iteration: 2,
            energy,
            energy_change: Some(0.0),
            density_residual: 0.0,
            max_eigen_residual: 0.0,
            requested_eigen_tolerance: o.eigen_tolerance,
            davidson_iterations: 1,
            charge_error: 0.0,
            input_charge_error: 0.0,
            min_density: 1.0,
            mixing_fraction: 1.0,
            mixer_reset: false,
            fixed_point_certificate: None,
            timings: PhaseTimings::default(),
        };
        assert!(converged(&good, &o));
        let mut bad = good.clone();
        bad.energy_change = None;
        assert!(!converged(&bad, &o));
        for failure in 0..6 {
            let mut bad = good.clone();
            match failure {
                0 => bad.energy_change = Some(2.0 * o.energy_tolerance),
                1 => bad.density_residual = 2.0 * o.density_tolerance,
                2 => bad.max_eigen_residual = 2.0 * o.eigen_tolerance,
                3 => bad.charge_error = 2e-8,
                4 => bad.requested_eigen_tolerance = 2.0 * o.eigen_tolerance,
                _ => bad.input_charge_error = 2e-8,
            }
            assert!(!converged(&bad, &o));
        }
        let mut bad = good;
        bad.density_residual = f64::NAN;
        assert!(!converged(&bad, &o));
    }
}
