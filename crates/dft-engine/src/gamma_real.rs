//! Default-off SCF bridge to ORC. Derived from real-gamma-design sections 4–10.
//! Completion is the only output boundary; downstream physics stays complex.
use crate::{ScfError, ScfOptions};
use dft_basis::PlaneWaveBasis;
use dft_mpi::Collective;
use dft_numerics::{RealMatrix, linalg::ComplexMatrix, pool::WorkerPool};
use dft_physics::{
    distributed_real_gamma::DistributedRealHamiltonian, hamiltonian::Hamiltonian,
    real_gamma::RealHamiltonian,
};
use dft_solvers::{
    BandTolerances, EigenSolution, RealLinearOperator, SolverError, real_davidson::RealDavidson,
};
use std::{cell::RefCell, sync::Arc};

/// On fails closed; Auto records why a complex solve was selected.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum GammaReal {
    Off,
    /// Default (ADR 0040): the real path when the problem qualifies, otherwise
    /// the complex path with the reason recorded.
    #[default]
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

/// Owned by one KState. The immutable Hamiltonian token guards physical data;
/// the key guards scratch sizing and solver options. Per-iteration tolerances
/// and certification policy are deliberately supplied afresh, not cached.
pub(crate) struct RealSolveState {
    op: Operator,
    key: RealSolveKey,
    workspace: Option<dft_solvers::real_davidson::RealDavidsonWorkspace>,
    /// Stick-distributed Hamiltonian of this rank (MPI stage M2c), built only when
    /// `parallel::m2_enabled()` and a communicator with two or more ranks is active.
    pub(crate) dist: Option<DistributedRealHamiltonian>,
    /// Local rows (stick partition) of the eigenvectors of the latest solve when
    /// it ran on the M2 path; None after any other solve (M2d density/energy).
    pub(crate) m2_local: Option<RealMatrix>,
    pub(crate) m2_density: Option<crate::m2::M2Density>,
}

#[derive(PartialEq, Eq)]
struct RealSolveKey {
    bands: usize,
    capacity: usize,
    threads: usize,
    max_iterations: usize,
    subspace: usize,
    capacity_per_band: usize,
    seed: u64,
    standard_threshold: Option<u64>,
    option: GammaReal,
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
    state: &mut Option<RealSolveState>,
) -> Result<Option<EigenSolution>, ScfError> {
    let option = options.gamma_real;
    let comm = crate::parallel::SharedComm::active();
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
    } else if comm.is_some()
        && options.preconditioner != dft_solvers::DavidsonPreconditioner::Default
    {
        Some("distributed Davidson supports the default preconditioner only")
    } else {
        None
    };
    if let Some(reason) = reason {
        *state = None;
        if diagnostics.fallback_reason.as_deref() != Some(reason) {
            diagnostics.fallback(option, reason)?;
        }
        return Ok(None);
    }
    let capacity = options
        .davidson_subspace
        .max(options.davidson_capacity_per_band.saturating_mul(bands))
        .min(kinetic.len());
    let key = RealSolveKey {
        bands,
        capacity,
        threads: options.threads,
        max_iterations: options.davidson_max_iterations,
        subspace: options.davidson_subspace,
        capacity_per_band: options.davidson_capacity_per_band,
        seed: options.seed,
        standard_threshold: options.standard.map(|p| p.threshold_ha.to_bits()),
        option,
    };
    let reusable = state.as_mut().is_some_and(|state| {
        state.key == key
            && Arc::ptr_eq(&state.op.pool, &pool)
            && state.op.h.get_mut().update_from_complex(hamiltonian, basis)
    });
    if !reusable {
        let real = match RealHamiltonian::from_complex(hamiltonian, basis, capacity, pool.clone()) {
            Ok(real) => real,
            Err(e) => {
                *state = None;
                diagnostics.fallback(option, format!("real Hamiltonian construction: {e}"))?;
                return Ok(None);
            }
        };
        *state = Some(RealSolveState {
            op: Operator {
                n: real.half_sphere().full_len(),
                h: RefCell::new(real),
                pool: pool.clone(),
            },
            key,
            workspace: None,
            dist: None,
            m2_local: None,
            m2_density: None,
        });
    }
    let state = state.as_mut().unwrap();
    state.m2_local = None;
    // M2c: the distributed Hamiltonian follows the SCF potential like the serial one.
    match comm.as_ref().filter(|_| crate::parallel::m2_enabled()) {
        Some(c) => {
            let fresh = state
                .dist
                .as_mut()
                .is_some_and(|d| d.update_from_complex(hamiltonian, basis));
            if !fresh {
                state.dist = DistributedRealHamiltonian::new(
                    hamiltonian,
                    basis,
                    capacity,
                    c.rank(),
                    c.size(),
                    pool.clone(),
                )
                .ok();
                state.m2_density = None;
            }
            if let Some(d) = state.dist.as_mut() {
                d.set_batched(crate::m2::batched_transposes());
            }
        }
        None => {
            state.dist = None;
            state.m2_density = None;
        }
    }
    let real = state.op.h.borrow();
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
    drop(real);
    let op = &state.op;
    let dist = state.dist.as_mut();
    // The SCF driver supplies the same ratcheted tolerance and band policy to
    // both representations. Standard returns remain explicitly uncertified.
    let result = if let (Some(comm), Some(dist)) = (&comm, dist) {
        crate::parallel::count_distributed_solve();
        crate::parallel::count_distributed_rows_solve();
        let local_diagonal: Vec<f64> = dist.local_rows().iter().map(|&r| diagonal[r]).collect();
        let result = solver.solve_distributed_rows(
            comm,
            dist,
            packed.as_ref(),
            &local_diagonal,
            policy,
            options.standard.is_some(),
            certify && options.standard.is_none(),
            options.davidson_capacity_per_band,
        )?;
        // M2d: keep this rank's rows for the distributed density and energy.
        state.m2_local = Some(
            dist.extract_local(&result.eigenvectors)
                .map_err(|e| ScfError::Numerical(e.to_string()))?,
        );
        result
    } else if let Some(comm) = &comm {
        crate::parallel::count_distributed_solve();
        solver.solve_distributed(
            comm,
            op,
            packed.as_ref(),
            &diagonal,
            policy,
            options.standard.is_some(),
            certify && options.standard.is_none(),
            options.davidson_capacity_per_band,
        )?
    } else {
        solver.solve_reusing_with_preconditioner(
            op,
            packed.as_ref(),
            &diagonal,
            policy,
            options.standard.is_some(),
            certify && options.standard.is_none(),
            options.davidson_capacity_per_band,
            &mut state.workspace,
            options.preconditioner,
        )?
    };
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

#[cfg(test)]
mod persistence_tests {
    include!("../tests/common/persistent_real_invalidation_unit.rs");
}
