//! Default-off SCF bridge to ORC. Derived from real-gamma-design sections 4–10.
//! Completion is the only output boundary; downstream physics stays complex.
use crate::{ScfError, ScfOptions};
use dft_basis::PlaneWaveBasis;
use dft_numerics::{RealMatrix, linalg::ComplexMatrix, pool::WorkerPool};
use dft_physics::{hamiltonian::Hamiltonian, real_gamma::RealHamiltonian};
use dft_solvers::{
    BandTolerances, EigenSolution, RealLinearOperator, SolverError, real_davidson::RealDavidson,
};
use std::{cell::RefCell, sync::Arc};

/// On fails closed; Auto records why a complex solve was selected.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum GammaReal {
    #[default]
    Off,
    Auto,
    On,
}
impl GammaReal {
    pub fn name(self) -> &'static str {
        match self {
            Self::Off => "off",
            Self::Auto => "auto",
            Self::On => "on",
        }
    }
}

/// Per-k dispatch evidence for the latest solve, also available before completion.
#[derive(Clone, Debug, Default)]
pub struct GammaRealDiagnostics {
    pub gamma_real_used: bool,
    pub fallback_reason: Option<String>,
}
impl GammaRealDiagnostics {
    pub(crate) fn fallback(
        &mut self,
        option: GammaReal,
        reason: impl Into<String>,
    ) -> Result<(), ScfError> {
        self.gamma_real_used = false;
        let reason = reason.into();
        if option == GammaReal::On {
            return Err(ScfError::Numerical(format!(
                "gamma_real=on is not eligible: {reason}"
            )));
        }
        self.fallback_reason = Some(reason);
        Ok(())
    }
}

struct Operator {
    h: RefCell<RealHamiltonian>,
    n: usize,
    pool: Arc<WorkerPool>,
}
impl RealLinearOperator for Operator {
    fn dimension(&self) -> usize {
        self.n
    }
    fn worker_pool(&self) -> Option<Arc<WorkerPool>> {
        Some(self.pool.clone())
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), SolverError> {
        let mut input = RealMatrix::zeros(self.n, 1)?;
        input.as_mut_slice().copy_from_slice(x);
        let mut out = RealMatrix::zeros(self.n, 1)?;
        self.apply_block(&input, &mut out)?;
        y.copy_from_slice(out.as_slice());
        Ok(())
    }
    fn apply_block(&self, x: &RealMatrix, y: &mut RealMatrix) -> Result<(), SolverError> {
        self.h
            .borrow_mut()
            .apply_block(x, y)
            .map_err(|e| SolverError::Application(e.to_string()))
    }
}

/// Try the real path without modifying an incompatible warm start. Construction
/// failures are eligibility failures; solver failures propagate, never fall back.
/// Certification follows the same caller policy as complex Davidson.
#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_real_gamma(
    hamiltonian: &Hamiltonian,
    basis: &PlaneWaveBasis,
    initial: Option<&ComplexMatrix>,
    kinetic: &[f64],
    bands: usize,
    tolerance: f64,
    policy: BandTolerances,
    certify: bool,
    options: &ScfOptions,
    pool: Arc<WorkerPool>,
    exclusion: Option<&str>,
    diagnostics: &mut GammaRealDiagnostics,
) -> Result<Option<EigenSolution>, ScfError> {
    let option = options.gamma_real;
    let reason = if option == GammaReal::Off {
        Some("option off")
    } else if let Some(reason) = exclusion {
        Some(reason)
    } else if basis.k_fractional() != [0.; 3] {
        Some("k is not exactly Gamma")
    } else if options.locking != dft_solvers::locking::Locking::Off {
        Some("locking enabled")
    } else if options.spectral_guard != crate::SpectralGuard::Off {
        Some("spectral guard enabled")
    } else if !matches!(options.eigensolver, dft_solvers::SolverChoice::Davidson) {
        Some("hybrid/RMM eigensolver")
    } else {
        None
    };
    if let Some(reason) = reason {
        if diagnostics.fallback_reason.as_deref() != Some(reason) {
            diagnostics.fallback(option, reason)?;
        }
        return Ok(None);
    }
    let capacity = options
        .davidson_subspace
        .max(options.davidson_capacity_per_band.saturating_mul(bands))
        .min(kinetic.len());
    let real = match RealHamiltonian::from_complex(hamiltonian, basis, capacity, pool.clone()) {
        Ok(real) => real,
        Err(e) => {
            diagnostics.fallback(option, format!("real Hamiltonian construction: {e}"))?;
            return Ok(None);
        }
    };
    let half = real.half_sphere();
    let n = half.full_len();
    let mut diagonal = vec![0.; n];
    diagonal[0] = kinetic[0];
    for (i, &p) in half.plus().iter().enumerate().skip(1) {
        diagonal[2 * i - 1] = kinetic[p as usize];
        diagonal[2 * i] = kinetic[p as usize];
    }
    let packed = if let Some(initial) = initial {
        let mut packed = RealMatrix::zeros(n, initial.cols()).map_err(SolverError::from)?;
        if let Err(e) = pool.install(|| {
            pair_columns(initial.as_slice(), packed.as_mut_slice(), n, &|src, dst| {
                half.pack(src, dst)
            })
        }) {
            diagnostics.fallback(
                option,
                format!("warm start is not ORC-compatible (1e-14 absolute): {e}"),
            )?;
            return Ok(None);
        }
        Some(packed)
    } else {
        None
    };
    let solver = RealDavidson {
        bands,
        tolerance,
        max_iterations: options.davidson_max_iterations,
        max_subspace: capacity,
        seed: options.seed,
    };
    let op = Operator {
        h: RefCell::new(real),
        n,
        pool,
    };
    // The SCF driver supplies the same ratcheted tolerance and band policy to
    // both representations. Standard returns remain explicitly uncertified.
    let result = solver.solve_with_capacity_per_band(
        &op,
        packed.as_ref(),
        &diagonal,
        policy,
        options.standard.is_some(),
        certify && options.standard.is_none(),
        options.davidson_capacity_per_band,
    )?;
    let mut eigenvectors = ComplexMatrix::zeros(n, bands).map_err(SolverError::from)?;
    {
        let h = op.h.borrow();
        let half = h.half_sphere();
        op.pool
            .install(|| {
                pair_columns(
                    result.eigenvectors.as_slice(),
                    eigenvectors.as_mut_slice(),
                    n,
                    &|src, dst| half.unpack(src, dst),
                )
            })
            .map_err(|e| ScfError::Numerical(e.to_string()))?;
    }
    diagnostics.gamma_real_used = true;
    diagnostics.fallback_reason = None;
    Ok(Some(EigenSolution {
        locking: None,
        spectral_completeness: None,
        eigenvalues: result.eigenvalues,
        eigenvectors,
        residuals: result.residuals,
        converged: result.converged,
        buffer: result.buffer,
        residuals_certified: result.residuals_certified,
        iterations: result.iterations,
        history: result.history,
        operator_applications: result.operator_applications,
    }))
}

// Boundary conversions use the same fixed pair partition as local application.
// They need only disjoint output slices; all storage belongs to the caller.
fn pair_columns<T: Sync, U: Send, E: Send>(
    input: &[T],
    output: &mut [U],
    n: usize,
    convert: &(impl Fn(&[T], &mut [U]) -> Result<(), E> + Sync),
) -> Result<(), E> {
    let pairs = (input.len() / n).div_ceil(2);
    if pairs > 1 {
        let split = (pairs / 2) * 2 * n;
        let (li, ri) = input.split_at(split);
        let (lo, ro) = output.split_at_mut(split);
        let (a, b) = WorkerPool::join(
            || pair_columns(li, lo, n, convert),
            || pair_columns(ri, ro, n, convert),
        );
        a?;
        b
    } else {
        for (src, dst) in input.chunks(n).zip(output.chunks_mut(n)) {
            convert(src, dst)?;
        }
        Ok(())
    }
}
