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
    dist: Option<DistributedRealHamiltonian>,
    /// NL-ENERGY-REAL: per-band <x|V_nl|x> (Ha) of the latest real solve.
    nl_values: Vec<f64>,
    nl_valid: bool,
}
impl RealSolveState {
    /// NL-ENERGY-REAL values of the latest solve of this state, if computed.
    pub(crate) fn nonlocal_partials(&self) -> Option<&[f64]> {
        self.nl_valid.then_some(self.nl_values.as_slice())
    }
}
/// NL-ENERGY-REAL (opt-in `BAANDI_NL_ENERGY_REAL=on`, default off; rounding
/// class): on the real Gamma path the audit's per-band nonlocal energies come
/// from one real projector product of the returned ORC vectors (P^T X by the
/// nonlocal GEMM) instead of the complex projector sums over the full sphere.
pub(crate) fn nl_energy_real() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("BAANDI_NL_ENERGY_REAL").as_deref() == Ok("on"))
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
    // B6 hooks (used by the solver only with `BAANDI_B6=on`).
    fn projector_count(&self) -> usize {
        self.h.borrow().projector_count()
    }
    fn copy_last_projections(&self, m: usize, out: &mut [f64]) -> bool {
        match self.h.borrow().last_projections(m) {
            Ok(p) if p.len() == out.len() => {
                out.copy_from_slice(p);
                true
            }
            _ => false,
        }
    }
    fn apply_block_with_projections(
        &self,
        x: &RealMatrix,
        projections: &[f64],
        y: &mut RealMatrix,
    ) -> Result<(), SolverError> {
        self.h
            .borrow_mut()
            .apply_block_with_projections(x, projections, y)
            .map_err(|e| SolverError::Application(e.to_string()))
    }
}

fn real_capacity(options: &ScfOptions, bands: usize, n: usize) -> usize {
    options
        .davidson_subspace
        .max(options.davidson_capacity_per_band.saturating_mul(bands))
        .min(n)
}
/// Reuse the state when its key, pool and immutable operator data match (the
/// local potential is refreshed), otherwise build it. Unchanged construction
/// logic of `solve_real_gamma`, shared with ATOMIC-REAL.
fn ensure_state(
    hamiltonian: &Hamiltonian,
    basis: &PlaneWaveBasis,
    bands: usize,
    capacity: usize,
    options: &ScfOptions,
    pool: &Arc<WorkerPool>,
    state: &mut Option<RealSolveState>,
) -> Result<(), dft_physics::RealGammaError> {
    let key = RealSolveKey {
        bands,
        capacity,
        threads: options.threads,
        max_iterations: options.davidson_max_iterations,
        subspace: options.davidson_subspace,
        capacity_per_band: options.davidson_capacity_per_band,
        seed: options.seed,
        standard_threshold: options.standard.map(|p| p.threshold_ha.to_bits()),
        option: options.gamma_real,
    };
    let reusable = state.as_mut().is_some_and(|state| {
        state.key == key
            && Arc::ptr_eq(&state.op.pool, pool)
            && state.op.h.get_mut().update_from_complex(hamiltonian, basis)
    });
    if !reusable {
        let real = RealHamiltonian::from_complex(hamiltonian, basis, capacity, pool.clone())?;
        *state = Some(RealSolveState {
            op: Operator {
                n: real.half_sphere().full_len(),
                h: RefCell::new(real),
                pool: pool.clone(),
            },
            key,
            workspace: None,
            dist: None,
            nl_values: Vec::new(),
            nl_valid: false,
        });
    }
    Ok(())
}

/// ATOMIC-REAL (opt-in `BAANDI_ATOMIC_REAL=on`, default off; rounding class).
pub(crate) fn atomic_real() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("BAANDI_ATOMIC_REAL").as_deref() == Ok("on"))
}

/// ATOMIC-REAL: H Q for the atomic-start Ritz step through the real Gamma
/// Hamiltonian (pair FFTs, real projector GEMM) instead of the complex apply.
/// The state is built exactly as `solve_real_gamma` builds it, so the first
/// solve reuses it. `None` (caller keeps the complex apply) unless the k state
/// is exactly Gamma, the real path is selected and eligible on the threads-only
/// path, Q packs to ORC within the existing tolerance and fits the capacity.
#[allow(clippy::too_many_arguments)]
pub(crate) fn atomic_apply_real(
    hamiltonian: &Hamiltonian,
    basis: &PlaneWaveBasis,
    q: &ComplexMatrix,
    bands: usize,
    options: &ScfOptions,
    pool: &Arc<WorkerPool>,
    exclusion: Option<&str>,
    state: &mut Option<RealSolveState>,
) -> Result<Option<ComplexMatrix>, ScfError> {
    if options.gamma_real == GammaReal::Off
        || exclusion.is_some()
        || basis.k_fractional() != [0.; 3]
        || options.locking != dft_solvers::locking::Locking::Off
        || options.spectral_guard != crate::SpectralGuard::Off
        || !matches!(options.eigensolver, dft_solvers::SolverChoice::Davidson)
        || crate::parallel::SharedComm::active().is_some()
    {
        return Ok(None);
    }
    let n = q.rows();
    let capacity = real_capacity(options, bands, hamiltonian.kinetic_diagonal().len());
    if q.cols() > capacity
        || ensure_state(hamiltonian, basis, bands, capacity, options, pool, state).is_err()
    {
        return Ok(None);
    }
    let state = state.as_mut().unwrap();
    let mut h = state.op.h.borrow_mut();
    if h.half_sphere().full_len() != n {
        return Ok(None);
    }
    let mut packed = RealMatrix::zeros(n, q.cols()).map_err(SolverError::from)?;
    let half = h.half_sphere();
    if pool
        .install(|| {
            pair_columns(q.as_slice(), packed.as_mut_slice(), n, &|src, dst| {
                half.pack(src, dst)
            })
        })
        .is_err()
    {
        return Ok(None);
    }
    let mut image = RealMatrix::zeros(n, q.cols()).map_err(SolverError::from)?;
    h.apply_block(&packed, &mut image)
        .map_err(|e| ScfError::Numerical(e.to_string()))?;
    let mut hq = ComplexMatrix::zeros(n, q.cols()).map_err(SolverError::from)?;
    let half = h.half_sphere();
    pool.install(|| {
        pair_columns(image.as_slice(), hq.as_mut_slice(), n, &|src, dst| {
            half.unpack(src, dst)
        })
    })
    .map_err(|e| ScfError::Numerical(e.to_string()))?;
    Ok(Some(hq))
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
    let capacity = real_capacity(options, bands, kinetic.len());
    if let Err(e) = ensure_state(hamiltonian, basis, bands, capacity, options, &pool, state) {
        *state = None;
        diagnostics.fallback(option, format!("real Hamiltonian construction: {e}"))?;
        return Ok(None);
    }
    let state = state.as_mut().unwrap();
    state.nl_valid = false;
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
            }
        }
        None => state.dist = None,
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
        solver.solve_distributed_rows(
            comm,
            dist,
            packed.as_ref(),
            &local_diagonal,
            policy,
            options.standard.is_some(),
            certify && options.standard.is_none(),
            options.davidson_capacity_per_band,
        )?
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
    // NL-ENERGY-REAL: nonlocal expectation values of the returned real vectors
    // for the energy audit of this map (threads-only path).
    if comm.is_none() && nl_energy_real() {
        state.nl_values.resize(bands, 0.0);
        state.nl_valid = op
            .h
            .borrow_mut()
            .nonlocal_expectation_values(&result.eigenvectors, &mut state.nl_values)
            .is_ok();
    }
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
mod atomic_real_tests {
    use super::*;
    use dft_geometry::Cell;
    use dft_physics::data::Array3;
    use dft_physics::hamiltonian::NonlocalOperator;
    use dft_physics::data::Complex64;

    /// ATOMIC-REAL: H Q through the real Hamiltonian equals the complex apply to
    /// rounding for Hermitian-completed Q, leaves a state the solve reuses, and
    /// declines non-Hermitian input and the off option.
    #[test]
    fn atomic_apply_real_matches_complex() {
        let cell =
            Cell::from_user_rows_bohr([[4.2, 0.3, 0.], [0., 3.9, 0.2], [0.1, 0., 4.4]]).unwrap();
        let options = ScfOptions {
            gamma_real: GammaReal::Auto,
            cutoff_hartree: 2.5,
            ..Default::default()
        };
        let basis = options.basis(&cell, [0.; 3], 1.).unwrap();
        let grid = dft_basis::FftGridSpec::new(&basis, &basis.density_basis(&cell).unwrap(), None)
            .unwrap();
        let potential = Array3::from_shape_fn(grid.dims(), |(i, j, k)| {
            0.1 * ((i + 2 * j + 3 * k) as f64).sin()
        });
        let h = Hamiltonian::new(
            &cell,
            &basis,
            potential,
            NonlocalOperator::new(&basis, cell.volume_bohr3(), &[]).unwrap(),
        )
        .unwrap();
        let half = dft_basis::HalfSphere::new(&basis, &grid).unwrap();
        let n = basis.vectors().len();
        let cols = 5;
        let mut q = ComplexMatrix::zeros(n, cols).unwrap();
        for j in 0..cols {
            let real: Vec<f64> = (0..n).map(|i| ((i * 3 + j * 7 + 1) as f64).sin()).collect();
            half.unpack(&real, q.column_mut(j)).unwrap();
        }
        let mut expected = ComplexMatrix::zeros(n, cols).unwrap();
        dft_solvers::LinearOperator::apply_block(&h, &q, &mut expected).unwrap();
        let pool = WorkerPool::new(2).unwrap();
        let mut state = None;
        let hq = atomic_apply_real(&h, &basis, &q, 4, &options, &pool, None, &mut state)
            .unwrap()
            .expect("eligible");
        let scale = expected.as_slice().iter().map(|z| z.norm()).fold(0., f64::max);
        let worst = hq
            .as_slice()
            .iter()
            .zip(expected.as_slice())
            .map(|(a, b)| (a - b).norm())
            .fold(0., f64::max);
        assert!(worst <= 1e-13 * scale, "{worst:e} of {scale:e}");
        assert!(state.is_some());
        // The solve's own construction check reuses the state built here.
        let before = state.as_ref().map(|s| s.op.h.as_ptr() as usize);
        let capacity = real_capacity(&options, 4, h.kinetic_diagonal().len());
        ensure_state(&h, &basis, 4, capacity, &options, &pool, &mut state).unwrap();
        assert_eq!(before, state.as_ref().map(|s| s.op.h.as_ptr() as usize));
        // Non-Hermitian input and the off option fall back to the complex apply.
        let mut bad = q.clone();
        bad.column_mut(0)[1] += Complex64::new(0., 1e-6);
        assert!(
            atomic_apply_real(&h, &basis, &bad, 4, &options, &pool, None, &mut state)
                .unwrap()
                .is_none()
        );
        let off = ScfOptions {
            gamma_real: GammaReal::Off,
            ..options.clone()
        };
        assert!(
            atomic_apply_real(&h, &basis, &q, 4, &off, &pool, None, &mut state)
                .unwrap()
                .is_none()
        );
    }
}

#[cfg(test)]
mod persistence_tests {
    include!("../tests/common/persistent_real_invalidation_unit.rs");
}
