//! Reference-SCF block Davidson (R-SCF-DAV), opt-in by direct call only.
//!
//! Derived from `docs/research/qe-study/davidson-iteration-structure-gamma.md`
//! (own-words algorithm note); no QE source used. All energies are Hartree.
//! Local rows use the operator's orthonormal real coordinates (ORC), including
//! the Gamma pair weights. Never apply an additional G=0 or factor-two weight.
//!
//! Ritz-change convergence is NOT a residual certificate or a proof of spectral
//! completeness. The caller must supply a start spanning the desired states and
//! independently verify the final SCF. An iteration limit returns an explicitly
//! unconverged result, so the driver can retry with the returned vectors.

mod accounting;
mod workspace;
pub use accounting::{Communication, Gemm, Iteration, Record, STAGES};
use accounting::{Counted, clock, elapsed};
use workspace::Workspace;

use crate::{SolverError, distributed_davidson::DistributedRealOperator};
use dft_mpi::Collective;
use dft_numerics::RealMatrix;

/// Placement is explicit so dev08 can compare identical arithmetic and counts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReducedPlacement {
    Replicated,
    Root,
}

#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub bands: usize,
    /// Occupied Ritz-change threshold, Ha (strict inequality).
    pub threshold: f64,
    /// Leading occupied bands; all bands use `threshold` if `first_step`.
    pub occupied: usize,
    pub first_step: bool,
    /// Zero selects min(global dimension, 2 * bands).
    pub capacity: usize,
    /// Default 20 correction iterations per call. Exhaustion is not success.
    pub max_iterations: usize,
    pub placement: ReducedPlacement,
    pub return_hx: bool,
    /// Stage clocks only; work/communication accounting is always available.
    pub profile: bool,
}
impl Options {
    pub fn new(bands: usize, threshold: f64, occupied: usize) -> Self {
        Self {
            bands,
            threshold,
            occupied,
            first_step: false,
            capacity: 0,
            max_iterations: 20,
            placement: ReducedPlacement::Replicated,
            return_hx: false,
            profile: false,
        }
    }
    pub fn band_threshold(&self, band: usize) -> f64 {
        if self.first_step || band < self.occupied {
            self.threshold
        } else {
            (5.0 * self.threshold).max(5e-6)
        }
    }
}

#[derive(Debug)]
pub struct Solution {
    pub eigenvalues: Vec<f64>,
    /// Local rows, in exactly `operator.local_rows()` order. No final gather.
    pub eigenvectors: RealMatrix,
    pub hx: Option<RealMatrix>,
    pub converged: Vec<bool>,
    /// Per correction, absolute Ritz changes in Ha, retained for trajectory audits.
    pub ritz_changes_ha: Vec<Vec<f64>>,
    /// Previous/current Ritz values per correction and band (Ha), for rounding audits.
    pub ritz_values_ha: Vec<Vec<[f64; 2]>>,
    /// Actual strict-inequality thresholds, including the empty-band rule (Ha).
    pub band_thresholds_ha: Vec<f64>,
    pub record: Record,
}

/// Solve H X = X E. `start` and `hamiltonian_diagonal` contain ONLY this rank's ORC rows;
/// `hamiltonian_diagonal` is the full plane-wave diagonal in Ha: kinetic plus
/// the local average plus the nonlocal diagonal. Share each Gamma-pair value
/// between its real/imaginary rows, without pair weights. Options are replicated.
/// The row partitions must cover the global dimension once (operator contract).
/// The communicator supplies fixed-tree sums; only rounding changes with ranks.
/// No retained H X is accepted: every call applies H to all starting columns.
/// Scratch and trace storage are reserved once, before the iteration loop.
pub fn solve<C: Collective, O: DistributedRealOperator>(
    comm: &C,
    operator: &mut O,
    start: &RealMatrix,
    hamiltonian_diagonal: &[f64],
    options: Options,
) -> Result<Solution, SolverError> {
    let total_clock = clock(options.profile);
    let c = Counted::new(comm);
    let n = options.bands;
    let dim = operator.global_len();
    let cap = if options.capacity == 0 {
        n.saturating_mul(2)
    } else {
        options.capacity
    }
    .min(dim);
    let bad = n == 0
        || n > dim
        || cap < n.saturating_mul(2).min(dim)
        || options.occupied > n
        || options.max_iterations == 0
        || !options.threshold.is_finite()
        || options.threshold <= 0.0
        || !options.band_threshold(n.saturating_sub(1)).is_finite()
        || start.rows() != operator.local_rows().len()
        || start.cols() != n
        || hamiltonian_diagonal.len() != start.rows()
        || !start
            .as_slice()
            .iter()
            .chain(hamiltonian_diagonal)
            .all(|v| v.is_finite());
    // Agree on local invalid inputs before any rank enters operator collectives.
    let mut bad_sum = [0.0];
    c.tree_sum_f64(&[f64::from(bad)], &mut bad_sum)
        .map_err(error)?;
    if bad_sum[0] != 0.0 {
        return Err(SolverError::InvalidParameter(
            "reference Davidson dimensions/threshold/start",
        ));
    }
    let mut w = Workspace::new(operator, n, cap)?;
    let mut record = Record::new(
        options.max_iterations,
        cap,
        options.profile,
        options.placement,
    );
    record.setup.communication = c.snapshot();
    w.v.as_mut_slice()[..start.rows() * n].copy_from_slice(start.as_slice());
    w.v.reshape_within_storage(start.rows(), n)?;
    w.hv.reshape_within_storage(start.rows(), n)?;
    let t = clock(options.profile);
    w.apply(&c, operator, start, true, &mut record.setup)?;
    record.setup.seconds[0] += elapsed(t);
    w.border(&c, 0, n, false, &mut record.setup, options.profile)?;
    if !w.factor(0, n, &mut record.setup, options.profile)? {
        return Err(SolverError::InvalidParameter(
            "reference Davidson start is numerically dependent",
        ));
    }
    w.border(&c, 0, n, true, &mut record.setup, options.profile)?;
    w.eigen(&c, 0, n, &mut record.setup, &options)?;
    record.setup.communication = c.snapshot();
    let mut m = n;
    let mut active = Vec::with_capacity(n);
    active.extend(0..n);
    let mut converged = vec![false; n];
    let mut iterations = 0;
    let mut ritz_changes_ha = vec![vec![0.; n]; options.max_iterations];
    let mut ritz_values_ha = vec![vec![[0.; 2]; n]; options.max_iterations];
    for (step, ritz_changes) in ritz_changes_ha.iter_mut().enumerate() {
        let before = c.snapshot();
        let r = &mut record.iterations[step];
        r.active_bands.extend_from_slice(&active);
        r.basis_before = m;
        // A smaller explicitly requested capacity may need a contraction before
        // the first correction. Never silently drop requested active roots.
        if m + active.len() > cap && m > n {
            w.contract(&c, m, r, options.profile)?;
            m = n;
        }
        if m + active.len() > cap && cap < dim {
            return Err(SolverError::InvalidParameter(
                "capacity cannot hold the active correction block; use at least 2*bands",
            ));
        }
        w.corrections(&c, hamiltonian_diagonal, &active, m, r, options.profile)?;
        let mut k = w.e.cols();
        // Exact eigenvectors have zero corrections; appending them would make S
        // singular. The general rank-repair path handles this too.
        if m + k > cap || k == 0 {
            k = w.repair(&c, m, r, options.profile)?;
        } else {
            w.append_e(m)?;
            w.border(&c, m, m + k, false, r, options.profile)?;
            if !w.factor(m, m + k, r, options.profile)? {
                k = w.repair(&c, m, r, options.profile)?;
            }
        }
        if k > 0 {
            let t = clock(options.profile);
            w.apply_corrections(&c, operator, m, r)?;
            r.seconds[0] += elapsed(t);
            w.border(&c, m, m + k, true, r, options.profile)?;
            w.old_values.copy_from_slice(&w.values[..n]);
            w.eigen(&c, m, m + k, r, &options)?;
            m += k;
        } else {
            w.old_values.copy_from_slice(&w.values[..n]);
            r.reduced_order = m;
        }
        active.clear();
        for (i, flag) in converged.iter_mut().enumerate() {
            let change = (w.values[i] - w.old_values[i]).abs();
            ritz_changes[i] = change;
            ritz_values_ha[step][i] = [w.old_values[i], w.values[i]];
            *flag = change < options.band_threshold(i);
            if !*flag {
                active.push(i);
            }
        }
        r.unconverged = active.len();
        r.communication = c.snapshot().difference(before);
        iterations = step + 1;
        if active.is_empty() || iterations == options.max_iterations {
            break;
        }
        if m + active.len() > cap {
            w.contract(&c, m, r, options.profile)?;
            m = n;
            r.communication = c.snapshot().difference(before);
        }
    }
    record.iterations.truncate(iterations);
    ritz_changes_ha.truncate(iterations);
    ritz_values_ha.truncate(iterations);
    w.finish(m, options.return_hx, &mut record.finish, options.profile)?;
    record.communication = c.snapshot();
    record.total_seconds = elapsed(total_clock);
    record.summarize();
    Ok(Solution {
        ritz_values_ha,
        band_thresholds_ha: (0..n).map(|i| options.band_threshold(i)).collect(),
        eigenvalues: w.values[..n].to_vec(),
        eigenvectors: w.x,
        hx: options.return_hx.then_some(w.hx),
        converged,
        ritz_changes_ha,
        record,
    })
}

fn error(e: dft_mpi::MpiError) -> SolverError {
    SolverError::Application(e.to_string())
}

/// Smooth positive diagonal regularizer. The dimensionless argument uses
/// 1 Ry = 0.5 Ha: d_Ry = (1 + x + sqrt(1 + (x-1)^2))/2, x=2(h-e).
/// The algebra below avoids cancellation for negative x. No material parameter.
pub fn preconditioner_denominator(hamiltonian_diagonal_ha: f64, eigenvalue_ha: f64) -> f64 {
    let x = 2.0 * (hamiltonian_diagonal_ha - eigenvalue_ha);
    let z = (x - 1.0).abs();
    0.25 * (if x >= 1.0 { 2.0 * x } else { 2.0 } + 1.0 / (z.hypot(1.0) + z))
}
