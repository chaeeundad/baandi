//! Fixed-density spectra. Reuses the SCF fixed-point map without feeding its
//! output density back into the Hamiltonian. No potential formula is duplicated.
use crate::{KPointSolution, OccupationPolicy, ScfError, ScfOptions, ScfProblem, SpeciesPotential};
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use dft_physics::data::{Array3, erf};
use std::fmt;

#[derive(Debug)]
pub enum NscfError {
    InvalidInput(&'static str),
    Setup(ScfError),
    KPoint { index: usize, cause: ScfError },
}
impl fmt::Display for NscfError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInput(s) => write!(f, "invalid NSCF input: {s}"),
            Self::Setup(e) => write!(f, "NSCF setup/potential: {e}"),
            Self::KPoint { index, cause } => write!(f, "NSCF k index {index}: {cause}"),
        }
    }
}
impl std::error::Error for NscfError {}
impl From<ScfError> for NscfError {
    fn from(cause: ScfError) -> Self {
        // The existing ordered SCF worker boundary encodes the index in this
        // prefix. Translate it once at the adapter boundary, never in a kernel.
        if let ScfError::Numerical(message) = &cause
            && let Some(index) = message
                .strip_prefix("k index ")
                .and_then(|s| s.split(':').next())
                .and_then(|s| s.parse().ok())
        {
            return Self::KPoint { index, cause };
        }
        Self::Setup(cause)
    }
}

#[derive(Debug)]
pub struct NscfResult {
    /// Original coordinates and weights (path coordinates are not wrapped).
    pub points: Vec<KPoint>,
    /// Eigenvalues and per-band residual norms in Hartree.
    pub solutions: Vec<KPointSolution>,
    pub converged: bool,
    pub tolerance_hartree: f64,
    pub electrons: f64,
}

/// Evaluate the fixed input density at independent k points using the same
/// deterministic worker machinery and Davidson as SCF. `options.eigen_tolerance`
/// is the residual-norm threshold (default 1e-8 Ha), strict for every band.
/// The caller supplies the original cell, atoms, pseudopotentials and cutoff.
/// The grid is fixed to the supplied density; no interpolation is performed.
pub fn nscf(
    cell: &Cell,
    atoms: &[Atom],
    potentials: &[SpeciesPotential<'_>],
    density: &Array3<f64>,
    points: &[KPoint],
    mut options: ScfOptions,
) -> Result<NscfResult, NscfError> {
    if density.iter().any(|v| !v.is_finite() || *v < -1e-14) {
        return Err(NscfError::InvalidInput(
            "density must be finite and nonnegative",
        ));
    }
    options.extra_bands = options.extra_bands.max(4);
    options.grid_override = Some([density.shape()[0], density.shape()[1], density.shape()[2]]);
    let tolerance = options.eigen_tolerance;
    let requested_bands = options.bands;
    let mut problem = ScfProblem::with_kpoints(
        cell,
        atoms,
        potentials,
        options,
        points,
        OccupationPolicy::Insulator {
            minimum_gap_hartree: 1e-7,
        },
    )?;
    let electrons = problem.electron_count();
    if requested_bands.is_some_and(|n| n < electrons as usize / 2 + 4) {
        return Err(NscfError::InvalidInput(
            "at least four empty bands required",
        ));
    }
    let charge = density.sum() * cell.volume_bohr3() / density.len() as f64;
    if (charge - electrons).abs() > 1e-8 * electrons {
        return Err(NscfError::InvalidInput(
            "density charge does not match system",
        ));
    }
    let step = problem.step_kpoints(density, None, tolerance)?;
    for (index, point) in step.k_points.iter().enumerate() {
        if point
            .solution
            .residuals
            .iter()
            .any(|r| !r.is_finite() || *r > tolerance)
        {
            return Err(NscfError::KPoint {
                index,
                cause: ScfError::Numerical("unconverged band residual".into()),
            });
        }
    }
    Ok(NscfResult {
        points: points.to_vec(),
        solutions: step.k_points,
        converged: true,
        tolerance_hartree: tolerance,
        electrons,
    })
}

/// Spin-unpolarized Gaussian DOS. Each eigenvalue carries two states times its
/// normalized k weight; sigma is the standard deviation, in the same units as
/// energies. Integrated DOS is the analytic integral from minus infinity.
pub fn gaussian_dos(
    spectra: &[Vec<f64>],
    weights: &[f64],
    energies: &[f64],
    sigma: f64,
) -> Result<Vec<[f64; 3]>, NscfError> {
    if spectra.is_empty()
        || spectra.len() != weights.len()
        || spectra[0].is_empty()
        || spectra
            .iter()
            .any(|s| s.len() != spectra[0].len() || s.iter().any(|e| !e.is_finite()))
        || weights.iter().any(|w| !w.is_finite() || *w <= 0.0)
        || (weights.iter().sum::<f64>() - 1.0).abs() > 1e-12
        || !sigma.is_finite()
        || sigma <= 0.0
        || energies.is_empty()
        || energies.iter().any(|e| !e.is_finite())
        || energies.windows(2).any(|p| p[1] <= p[0])
    {
        return Err(NscfError::InvalidInput(
            "invalid Gaussian spectrum/grid/weights/sigma",
        ));
    }
    let norm = 1.0 / (sigma * (2.0 * std::f64::consts::PI).sqrt());
    if !norm.is_finite() || norm <= 0.0 {
        return Err(NscfError::InvalidInput(
            "Gaussian width is not representable",
        ));
    }
    Ok(energies
        .iter()
        .map(|&energy| {
            let mut dos = 0.0;
            let mut integrated = 0.0;
            for (spectrum, &weight) in spectra.iter().zip(weights) {
                for &eigenvalue in spectrum {
                    let t = (energy - eigenvalue) / sigma;
                    dos += 2.0 * weight * norm * (-0.5 * t * t).exp();
                    integrated += weight * (1.0 + erf(t / std::f64::consts::SQRT_2));
                }
            }
            [energy, dos, integrated]
        })
        .collect())
}
