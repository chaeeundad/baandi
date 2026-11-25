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

pub(crate) mod accounting;
pub(crate) mod ideal_reset;
pub use ideal_reset::ResetStats;
pub(crate) mod stream;
pub mod work;
pub(crate) mod workspace;
pub use accounting::{
    APPLY_H_STAGES, ApplyH, Communication, Gemm, Iteration, QeClockScope, QeCounters, Record,
    STAGES, apply_h_active, qe_clock_end, qe_clock_start, record_apply_h, take_apply_h_solves,
};
use accounting::{Counted, clock, elapsed};
use work::{Stage, span};
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
    /// Stage clocks and per-band Ritz history. Work/communication accounting
    /// is always available; false allocates no Ritz history.
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
    /// Per correction, absolute Ritz changes in Ha; empty when profile is false.
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
/// This entry always applies H to all starting columns.
/// Scratch and trace storage are reserved once, before the iteration loop.
pub fn solve<C: Collective, O: DistributedRealOperator>(
    comm: &C,
    operator: &mut O,
    start: &RealMatrix,
    hamiltonian_diagonal: &[f64],
    options: Options,
) -> Result<Solution, SolverError> {
    solve_with_initial_hx(comm, operator, start, hamiltonian_diagonal, options, None)
}

/// As `solve`, with an optional H(start) supplied by the caller. The caller
/// guarantees the exact start, fixed operator and current potential identity.
/// All ranks must supply a product or none. Shape, finiteness, overlap, rank and
/// convergence checks remain active; the start is not rotated before use.
pub fn solve_with_initial_hx<C: Collective, O: DistributedRealOperator>(
    comm: &C,
    operator: &mut O,
    start: &RealMatrix,
    hamiltonian_diagonal: &[f64],
    options: Options,
    initial_hx: Option<&RealMatrix>,
) -> Result<Solution, SolverError> {
    let apply_h_scope = accounting::ApplyHScope::new(options.profile, options.max_iterations + 1);
    let total_clock = clock(options.profile);
    let c = Counted::new(comm, options.profile);
    let n = options.bands;
    let dim = operator.global_len();
    let cap = workspace::capacity(dim, n, options.capacity);
    let work_scope = work::Scope::new(options.profile, options.max_iterations, cap);

    let mut bad = initial_hx.is_some_and(|hx| hx.rows() != start.rows() || hx.cols() != n)
        || n == 0
        || n > dim
        || cap < n.saturating_mul(2).min(dim)
        || options.occupied > n
        || options.max_iterations == 0
        || !options.threshold.is_finite()
        || options.threshold <= 0.0
        || !options.band_threshold(n.saturating_sub(1)).is_finite()
        || start.rows() != operator.local_rows().len()
        || start.cols() != n
        || hamiltonian_diagonal.len() != start.rows();
    // Allocate only after local shapes are valid. Copy/finite classification
    // share loads; collective input agreement still precedes every H call.

    bad |= initial_hx.is_some_and(|hx| hx.as_slice().iter().any(|v| !v.is_finite()))
        || !start
            .as_slice()
            .iter()
            .chain(hamiltonian_diagonal)
            .all(|v| v.is_finite());
    work::add(options.profile, |s| {
        s.streamed_bytes += 8
            * (start.rows() * n + start.rows() + initial_hx.map_or(0, |hx| hx.as_slice().len()))
                as u64;
    });
    // Agree on local invalid inputs before any rank enters operator collectives.
    let mut bad_sum = [0.0; 2];
    c.tree_sum_f64(
        &[f64::from(bad), f64::from(initial_hx.is_some())],
        &mut bad_sum,
    )
    .map_err(error)?;
    if bad_sum[0] != 0.0 || (bad_sum[1] != 0.0 && bad_sum[1] != comm.size() as f64) {
        return Err(SolverError::InvalidParameter(
            "reference Davidson dimensions/threshold/start",
        ));
    }
    let mut w = Workspace::new(operator, n, cap)?;
    work::add(options.profile, |s| {
        // Explicit arrays of Workspace::new; shared adapter scratch excluded.
        let doubles =
            4 * start.rows() * cap + 2 * start.rows() * n + 6 * cap * cap + 3 * cap + n + 2
                - w.reused_elements;
        s.zeroed_bytes += 8 * doubles as u64;
        s.streamed_bytes += 8 * doubles as u64;
    });
    {
        w.diagonal_bounds = hamiltonian_diagonal
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &h| {
                (lo.min(h), hi.max(h))
            });
        work::add(options.profile, |s| {
            s.streamed_bytes += 8 * start.rows() as u64
        });
    }
    // Resolve once per solve. The reference workspace resets C to identity
    // on contraction, even though it remeasures S/H and their factor afterward.
    // Only the immediately following correction can use these stored Ritz rows.

    work::add(options.profile, |s| {
        s.zeroed_bytes += 16 * n as u64;
        s.streamed_bytes += 16 * n as u64;
    });
    let mut restart_scratch = (vec![0.0; n], vec![0.0; n]);
    let mut contracted = false;
    let mut record = Record::new(
        options.max_iterations,
        cap,
        options.profile,
        options.placement,
    );
    record.setup.communication = c.snapshot();
    {
        work::copy(options.profile, start.rows() * n);
        w.v.as_mut_slice()[..start.rows() * n].copy_from_slice(start.as_slice());
    }
    w.v.reshape_within_storage(start.rows(), n)?;
    w.hv.reshape_within_storage(start.rows(), n)?;
    let t = clock(options.profile);
    if let Some(hx) = initial_hx {
        {
            w.hv.as_mut_slice().copy_from_slice(hx.as_slice());
            work::copy(options.profile, start.rows() * n);
        }
    } else {
        w.apply(&c, operator, start, true, &mut record.setup)?;
    }
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
    let mut ritz_changes_ha = if options.profile {
        vec![vec![0.; n]; options.max_iterations]
    } else {
        Vec::new()
    };
    let mut ritz_values_ha = if options.profile {
        vec![vec![[0.; 2]; n]; options.max_iterations]
    } else {
        Vec::new()
    };
    for (step, r) in record.iterations.iter_mut().enumerate() {
        let _glue = span(options.profile, Stage::Glue);
        let before = c.snapshot();
        r.active_bands.extend_from_slice(&active);
        r.basis_before = m;
        // A smaller explicitly requested capacity may need a contraction before
        // the first correction. Never silently drop requested active roots.
        if m + active.len() > cap && m > n {
            w.contract(&c, m, r, options.profile)?;
            m = n;
            contracted = true;
        }
        if m + active.len() > cap && cap < dim {
            return Err(SolverError::InvalidParameter(
                "capacity cannot hold the active correction block; use at least 2*bands",
            ));
        }
        if contracted {
            let (norms, sums) = &mut restart_scratch;
            restart_corrections(
                &c,
                &mut w,
                hamiltonian_diagonal,
                &active,
                norms,
                sums,
                r,
                options.profile,
            )?;
        } else {
            w.corrections(&c, hamiltonian_diagonal, &active, m, r, options.profile)?;
        }
        // Neither a new reduced solve nor rank repair preserves identity C.
        contracted = false;
        let mut k = w.e.cols();
        // Exact eigenvectors have zero corrections; appending them would make S
        // singular. The general rank-repair path handles this too.
        if m + k > cap || k == 0 {
            k = w.repair(&c, m, r, options.profile)?;
        } else {
            w.append_e(m, options.profile)?;
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
            work::copy(options.profile, n);
            w.old_values.copy_from_slice(&w.values[..n]);
            w.eigen(&c, m, m + k, r, &options)?;
            m += k;
        } else {
            work::copy(options.profile, n);
            w.old_values.copy_from_slice(&w.values[..n]);
            r.reduced_order = m;
        }
        let decision_span = span(options.profile, Stage::Decisions);
        active.clear();
        for (i, flag) in converged.iter_mut().enumerate() {
            let change = (w.values[i] - w.old_values[i]).abs();
            if options.profile {
                ritz_changes_ha[step][i] = change;
                ritz_values_ha[step][i] = [w.old_values[i], w.values[i]];
            }
            *flag = change < options.band_threshold(i);
            if !*flag {
                active.push(i);
            }
        }
        r.unconverged = active.len();
        work::add(options.profile, |s| s.scalar_flops += n as u64);
        drop(decision_span);
        r.communication = c.snapshot().difference(before);
        iterations = step + 1;
        if active.is_empty() || iterations == options.max_iterations {
            break;
        }
        if m + active.len() > cap {
            w.contract(&c, m, r, options.profile)?;
            m = n;
            contracted = true;
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
    if options.profile {
        accounting::qe_clock_add(5, record.total_seconds, 1);
        for r in std::iter::once(&record.setup)
            .chain(&record.iterations)
            .chain(std::iter::once(&record.finish))
        {
            accounting::qe_clock_add(6, r.seconds[0], r.h_calls);
            accounting::qe_clock_add(7, r.seconds[3], r.reduced_solves);
            accounting::qe_clock_add(8, r.seconds[2], 0);
        }
    }
    record.apply_h = apply_h_scope.finish();
    work_scope.finish();
    w.recycle_rows();
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

/// C is exactly identity immediately after contract; V/HV already are the
/// retained Ritz vectors/products. Keep the baseline residual, summation,
/// normalization, compaction and drop tests in their original order. This is
/// not stale-band locking and does not skip the measured restart refresh.
#[allow(clippy::too_many_arguments)]
fn restart_corrections<C: Collective>(
    c: &C,
    w: &mut Workspace,
    diagonal: &[f64],
    active: &[usize],
    norms: &mut [f64],
    sums: &mut [f64],
    r: &mut Iteration,
    profile: bool,
) -> Result<(), SolverError> {
    let _stage = span(profile, Stage::Corrections);
    let t = clock(profile);
    let precondition = span(profile, Stage::Precondition);
    let rows = w.v.rows();
    let k = active.len();
    w.e.reshape_within_storage(rows, k)?;
    {
        let mut sqrt_columns = 0;
        for (j, &band) in active.iter().enumerate() {
            let x = &w.v.as_slice()[band * rows..(band + 1) * rows];
            let hx = &w.hv.as_slice()[band * rows..(band + 1) * rows];
            let out = &mut w.e.as_mut_slice()[j * rows..(j + 1) * rows];
            let safe = stream::sqrt_safe(w.diagonal_bounds, w.values[band]);
            norms[j] = if safe {
                stream::column::<true>(diagonal, x, hx, out, w.values[band])
            } else {
                stream::column::<false>(diagonal, x, hx, out, w.values[band])
            };
            sqrt_columns += usize::from(safe);
        }
        stream::count(profile, rows, sqrt_columns, k - sqrt_columns);
    }
    drop(precondition);
    w.normalize_restart(c, norms, sums, r, profile)?;
    r.seconds[4] += elapsed(t);
    Ok(())
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

#[cfg(test)]
mod diagnostic_tests {
    use super::*;

    struct Tridiagonal(Vec<usize>);
    impl DistributedRealOperator for Tridiagonal {
        fn global_len(&self) -> usize {
            self.0.len()
        }
        fn local_rows(&self) -> &[usize] {
            &self.0
        }
        fn apply_block<C: Collective>(
            &mut self,
            _: &C,
            x: &RealMatrix,
            y: &mut RealMatrix,
        ) -> Result<(), SolverError> {
            for j in 0..x.cols() {
                for i in 0..x.rows() {
                    y[(i, j)] = (i + 1) as f64 * x[(i, j)];
                    if i > 0 {
                        y[(i, j)] += 0.1 * x[(i - 1, j)];
                    }
                    if i + 1 < x.rows() {
                        y[(i, j)] += 0.1 * x[(i + 1, j)];
                    }
                }
            }
            Ok(())
        }
    }

    struct Partitioned {
        dim: usize,
        rows: Vec<usize>,
        clustered: bool,
    }
    impl Partitioned {
        fn diagonal(&self, row: usize) -> f64 {
            if self.clustered && row < 2 {
                1.0
            } else {
                1.0 + row as f64 * 0.4
            }
        }
    }
    impl DistributedRealOperator for Partitioned {
        fn global_len(&self) -> usize {
            self.dim
        }
        fn local_rows(&self) -> &[usize] {
            &self.rows
        }
        fn apply_block<C: Collective>(
            &mut self,
            c: &C,
            x: &RealMatrix,
            y: &mut RealMatrix,
        ) -> Result<(), SolverError> {
            for j in 0..x.cols() {
                let dot: f64 = self
                    .rows
                    .iter()
                    .enumerate()
                    .map(|(i, &row)| ((row + 1) as f64).sin() * x[(i, j)])
                    .sum();
                let mut sum = [0.0];
                c.tree_sum_f64(&[dot], &mut sum).map_err(error)?;
                for (i, &row) in self.rows.iter().enumerate() {
                    y[(i, j)] =
                        self.diagonal(row) * x[(i, j)] + 0.03 * ((row + 1) as f64).sin() * sum[0];
                }
            }
            Ok(())
        }
    }

    fn same_bits(a: &[f64], b: &[f64]) {
        assert_eq!(
            a.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            b.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
        );
    }

    #[test]
    fn restart_corrections_preserve_decisions_columns_and_remove_identity_flops() {
        use dft_numerics::switches::RunConfig;
        let off = RunConfig::from_pairs(&[]).unwrap();
        let on = RunConfig::from_pairs(&[]).unwrap();
        for (dim, n, clustered) in [(12, 3, false), (17, 5, true), (7, 4, true)] {
            for ranks in [1, 2, 3, 8] {
                dft_mpi::simulate(ranks, |comm| {
                    // Cyclic partition deliberately differs from the usual block rows.
                    let rows: Vec<_> = (comm.rank()..dim).step_by(ranks).collect();
                    let mut op = Partitioned {
                        dim,
                        rows,
                        clustered,
                    };
                    let mut x = RealMatrix::zeros(op.rows.len(), n).unwrap();
                    for j in 0..n {
                        for (i, &row) in op.rows.iter().enumerate() {
                            x[(i, j)] = if row == j {
                                1.0
                            } else {
                                0.03 * ((row * 7 + j * 13 + 1) as f64).sin()
                            };
                        }
                    }
                    let diag: Vec<_> = op
                        .rows
                        .iter()
                        .map(|&row| op.diagonal(row) + 0.03 * ((row + 1) as f64).sin().powi(2))
                        .collect();
                    for placement in [ReducedPlacement::Root, ReducedPlacement::Replicated] {
                        for budget in [1, 40] {
                            let mut options = Options::new(n, 1e-10, n - 1);
                            options.max_iterations = budget;
                            options.placement = placement;
                            options.profile = true;
                            options.return_hx = true;
                            let a = {
                                let _scope = off.enter();
                                solve(comm, &mut op, &x, &diag, options).unwrap()
                            };
                            let b = {
                                let _scope = on.enter();
                                solve(comm, &mut op, &x, &diag, options).unwrap()
                            };
                            same_bits(&a.eigenvalues, &b.eigenvalues);
                            same_bits(a.eigenvectors.as_slice(), b.eigenvectors.as_slice());
                            same_bits(
                                a.hx.as_ref().unwrap().as_slice(),
                                b.hx.as_ref().unwrap().as_slice(),
                            );
                            assert_eq!(a.converged, b.converged);
                            assert_eq!(a.ritz_values_ha, b.ritz_values_ha);
                            assert_eq!(a.ritz_changes_ha, b.ritz_changes_ha);
                            assert_eq!(a.record.communication, b.record.communication);
                            assert_eq!(a.record.h_columns, b.record.h_columns);
                            assert_eq!(a.record.h_calls, b.record.h_calls);
                            assert_eq!(a.record.restarts, b.record.restarts);
                            assert_eq!(a.record.iterations.len(), b.record.iterations.len());
                            let saved = 0;
                            let mut restart = false;
                            for (ar, br) in a.record.iterations.iter().zip(&b.record.iterations) {
                                assert_eq!(ar.active_bands, br.active_bands);
                                assert_eq!(ar.h_columns, br.h_columns);
                                assert_eq!(ar.reduced_order, br.reduced_order);
                                assert_eq!(ar.restarts, br.restarts);
                                assert_eq!(ar.repairs, br.repairs);
                                assert_eq!(ar.dropped_corrections, br.dropped_corrections);
                                assert_eq!(ar.unconverged, br.unconverged);
                                assert_eq!(ar.gemms, br.gemms);
                                if restart {
                                    assert!(
                                        !br.gemms.iter().any(|g| matches!(
                                            g.stage,
                                            "residual_x" | "residual_hx"
                                        ))
                                    );
                                }
                                restart = ar.restarts > 0;
                            }
                            assert_eq!(a.record.gemm_flops, b.record.gemm_flops);
                            if budget > 1 {
                                assert!(a.converged.iter().all(|v| *v));
                                if dim < 2 * n {
                                    assert!(a.record.iterations.iter().any(|it| it.repairs > 0));
                                    assert!(
                                        a.record
                                            .iterations
                                            .iter()
                                            .any(|it| it.dropped_corrections > 0)
                                    );
                                }
                            } else {
                                assert_eq!(saved, 0);
                            }
                            if comm.rank() == 0 && placement == ReducedPlacement::Root && budget > 1
                            {
                                eprintln!(
                                    "DAV_DENSE dim={dim} bands={n} ranks={ranks} iterations={} h_columns={} restarts={} saved_flops_rank0={saved} delta_e_ha=0 bits=identical",
                                    a.record.iterations.len(),
                                    a.record.h_columns,
                                    a.record.restarts
                                );
                            }
                        }
                    }
                });
            }
        }
    }

    #[test]
    fn optional_ritz_storage_and_clocks_preserve_science_and_work() {
        dft_mpi::simulate(1, |comm| {
            let mut op = Tridiagonal((0..12).collect());
            let mut initial = RealMatrix::zeros(12, 3).unwrap();
            for j in 0..3 {
                for i in 0..12 {
                    initial[(i, j)] = if i == j {
                        1.
                    } else {
                        0.01 * ((i * 7 + j) as f64).sin()
                    };
                }
            }
            let diagonal: Vec<_> = (1..=12).map(f64::from).collect();
            for budget in [1, 20] {
                let mut options = Options::new(3, 1e-10, 2);
                options.max_iterations = budget;
                options.return_hx = true;
                let qe_scope = QeClockScope::begin();
                let light = solve(comm, &mut op, &initial, &diagonal, options).unwrap();
                assert_eq!(qe_scope.snapshot().calls, [0; 9]);
                options.profile = true;
                let full = solve(comm, &mut op, &initial, &diagonal, options).unwrap();
                let clocks = qe_scope.snapshot();
                assert_eq!(clocks.calls[5], 1);
                assert_eq!(clocks.calls[6], full.record.h_calls);
                assert_eq!(clocks.calls[7], full.record.reduced_solves);
                assert_eq!(clocks.seconds[5], full.record.total_seconds);
                assert!(light.ritz_changes_ha.is_empty());
                assert_eq!(light.ritz_changes_ha.capacity(), 0);
                assert_eq!(light.record.total_seconds, 0.);
                assert!(full.record.total_seconds > 0.);
                assert_eq!(full.ritz_changes_ha.len(), full.record.iterations.len());
                assert!(full.ritz_changes_ha.iter().all(|row| row.len() == 3));
                assert_eq!(light.eigenvalues, full.eigenvalues);
                assert_eq!(light.eigenvectors.as_slice(), full.eigenvectors.as_slice());
                assert_eq!(light.hx.unwrap().as_slice(), full.hx.unwrap().as_slice());
                assert_eq!(light.converged, full.converged);
                assert_eq!(light.record.communication, full.record.communication);
                assert_eq!(light.record.h_columns, full.record.h_columns);
                assert_eq!(light.record.gemm_flops, full.record.gemm_flops);
                assert_eq!(light.record.iterations.len(), full.record.iterations.len());
                for (a, b) in std::iter::once(&light.record.setup)
                    .chain(&light.record.iterations)
                    .chain(std::iter::once(&light.record.finish))
                    .zip(
                        std::iter::once(&full.record.setup)
                            .chain(&full.record.iterations)
                            .chain(std::iter::once(&full.record.finish)),
                    )
                {
                    assert_eq!(a.seconds, [0.; 8]);
                    assert_eq!(a.active_bands, b.active_bands);
                    assert_eq!(a.gemms, b.gemms);
                    assert_eq!(a.reduced_order, b.reduced_order);
                    assert_eq!(a.restarts, b.restarts);
                    assert_eq!(a.repairs, b.repairs);
                    assert_eq!(a.unconverged, b.unconverged);
                }
                if budget == 1 {
                    assert!(full.converged.iter().any(|v| !v));
                } else {
                    assert!(full.converged.iter().all(|v| *v));
                }
            }
        });
    }
    #[test]
    fn davrem_switches_preserve_bits_and_account_for_removed_work() {
        use dft_numerics::switches::RunConfig;
        let base = RunConfig::from_pairs(&[]).unwrap();
        let arms: Vec<_> = [("on", "off"), ("off", "on"), ("on", "on")]
            .into_iter()
            .map(|(swap, reset)| (swap, reset, RunConfig::from_pairs(&[]).unwrap()))
            .collect();
        for (dim, n, capacity, clustered) in [(12, 3, 6, false), (17, 5, 15, true), (7, 4, 7, true)]
        {
            for ranks in [1, 3, 8] {
                dft_mpi::simulate(ranks, |comm| {
                    let mut op = Partitioned {
                        dim,
                        rows: (comm.rank()..dim).step_by(ranks).collect(),
                        clustered,
                    };
                    let rows = op.rows.len();
                    let mut x = RealMatrix::zeros(rows, n).unwrap();
                    for j in 0..n {
                        for (i, &row) in op.rows.iter().enumerate() {
                            x[(i, j)] = if row == j {
                                1.
                            } else {
                                0.03 * ((row * 7 + j * 13 + 1) as f64).sin()
                            };
                        }
                    }
                    let diag: Vec<_> = op
                        .rows
                        .iter()
                        .map(|&row| op.diagonal(row) + 0.03 * ((row + 1) as f64).sin().powi(2))
                        .collect();
                    for placement in [ReducedPlacement::Root, ReducedPlacement::Replicated] {
                        let mut options = Options::new(n, 1e-10, n - 1);
                        options.max_iterations = 40;
                        options.capacity = capacity;
                        options.placement = placement;
                        options.return_hx = true;
                        // Full vs off: no retained work, clocks or trace storage off.
                        work::take_solves();
                        let light = {
                            let _g = base.enter();
                            solve(comm, &mut op, &x, &diag, options).unwrap()
                        };
                        assert!(work::take_solves().is_empty());
                        options.profile = true;
                        let a = {
                            workspace::release_row_cache();
                            let _g = base.enter();
                            solve(comm, &mut op, &x, &diag, options).unwrap()
                        };
                        let aw = work::take_solves().pop().unwrap();
                        same_bits(&light.eigenvalues, &a.eigenvalues);
                        same_bits(light.eigenvectors.as_slice(), a.eigenvectors.as_slice());
                        assert_eq!(light.record.h_columns, a.record.h_columns);
                        assert_eq!(
                            aw.stages.iter().map(|s| s.gemm_flops).sum::<u64>(),
                            a.record.gemm_flops
                        );
                        let counts = a.record.communication;
                        for (kind, calls, bytes) in [
                            (1, counts.reductions, counts.reduction_bytes),
                            (2, counts.broadcasts, counts.broadcast_bytes),
                            (3, counts.gathers, counts.gather_bytes),
                            (4, counts.alltoalls, counts.alltoall_bytes),
                        ] {
                            assert_eq!(
                                aw.stages
                                    .iter()
                                    .map(|s| s.collective_calls[kind])
                                    .sum::<u64>(),
                                calls
                            );
                            assert_eq!(
                                aw.stages
                                    .iter()
                                    .map(|s| s.collective_bytes[kind])
                                    .sum::<u64>(),
                                bytes
                            );
                        }
                        let work_s: f64 = aw.stages.iter().map(|s| s.seconds).sum();
                        assert!(work_s > 0.);
                        assert!(aw.stages.iter().all(|s| s.collective_seconds <= s.seconds));
                        assert!(a.converged.iter().all(|&v| v));
                        for &(swap, reset, config) in &arms {
                            let b = {
                                let _g = config.enter();
                                solve(comm, &mut op, &x, &diag, options).unwrap()
                            };
                            let bw = work::take_solves().pop().unwrap();
                            same_bits(&a.eigenvalues, &b.eigenvalues);
                            same_bits(a.eigenvectors.as_slice(), b.eigenvectors.as_slice());
                            same_bits(
                                a.hx.as_ref().unwrap().as_slice(),
                                b.hx.as_ref().unwrap().as_slice(),
                            );
                            assert_eq!(a.converged, b.converged);
                            assert_eq!(a.ritz_values_ha, b.ritz_values_ha);
                            assert_eq!(a.record.communication, b.record.communication);
                            assert_eq!(a.record.h_columns, b.record.h_columns);
                            assert_eq!(a.record.h_calls, b.record.h_calls);
                            assert_eq!(a.record.iterations.len(), b.record.iterations.len());
                            for (a, b) in a.record.iterations.iter().zip(&b.record.iterations) {
                                assert_eq!(a.active_bands, b.active_bands);
                                assert_eq!(a.reduced_order, b.reduced_order);
                                assert_eq!(a.repairs, b.repairs);
                                assert_eq!(a.dropped_corrections, b.dropped_corrections);
                                assert_eq!(a.gemms, b.gemms);
                            }
                            assert_eq!(
                                aw.stages[Stage::Restart as usize].copied_bytes,
                                bw.stages[Stage::Restart as usize].copied_bytes
                            );
                            let expected = 0;
                            let reset_saved = 0;
                            if comm.rank() == 0 && placement == ReducedPlacement::Root {
                                eprintln!(
                                    "DAVREM dim={dim} bands={n} ranks={ranks} swap={swap} reset={reset} H={} iterations={} restarts={} restart_copy_saved_B={expected} reset_write_saved_B={reset_saved} delta_e_Ha=0 bits=identical",
                                    a.record.h_columns,
                                    a.record.iterations.len(),
                                    a.record.restarts
                                );
                            }
                        }
                    }
                });
            }
        }
    }

    #[test]
    fn davrem_one_correction_has_analytic_row_work_and_error_cleanup() {
        dft_mpi::simulate(1, |comm| {
            // The analytic counts below describe the P3-DAVREM baseline passes; the
            // P3-STREAM switches (default on since 2026-11-25) change the counters.
            let pinned = dft_numerics::switches::RunConfig::from_pairs(&[]).unwrap();
            let _scope = pinned.enter();
            let mut op = Tridiagonal((0..12).collect());
            let mut x = RealMatrix::zeros(12, 3).unwrap();
            for j in 0..3 {
                x[(j, j)] = 1.;
            }
            let diag: Vec<_> = (1..=12).map(f64::from).collect();
            let mut options = Options::new(3, 1e-10, 3);
            options.profile = true;
            options.max_iterations = 1;
            work::take_solves();
            let a = solve(comm, &mut op, &x, &diag, options).unwrap();
            let w = work::take_solves().pop().unwrap();
            let pre = &w.stages[Stage::Precondition as usize];
            assert_eq!(pre.hypot_calls, 0);
            assert_eq!(pre.sqrt_calls, 36);
            assert_eq!(pre.scalar_flops, 15 * 36);
            assert_eq!(pre.streamed_bytes, 32 * 36 + 24);
            assert_eq!(w.stages[Stage::Normalize as usize].sqrt_calls, 3);
            assert_eq!(w.stages[Stage::Normalize as usize].collective_calls[1], 1);
            assert_eq!(w.stages[Stage::Normalize as usize].collective_bytes[1], 24);
            assert_eq!(w.stages[Stage::Return as usize].gemm_calls, 1);
            assert_eq!(
                w.gemms.len() as u64,
                w.stages.iter().map(|s| s.gemm_calls).sum::<u64>()
            );
            assert_eq!(
                w.stages.iter().map(|s| s.gemm_flops).sum::<u64>(),
                a.record.gemm_flops
            );
            x[(0, 0)] = f64::NAN;
            assert!(solve(comm, &mut op, &x, &diag, options).is_err());
            assert!(work::take_solves().is_empty());
            x[(0, 0)] = 1.;
            solve(comm, &mut op, &x, &diag, options).unwrap();
            assert_eq!(work::take_solves().len(), 1);
        });
    }
    #[test]
    fn stream_switches_preserve_decisions_and_columns() {
        use dft_numerics::switches::RunConfig;
        for (dim, n, capacity, clustered) in [(12, 3, 6, false), (17, 5, 15, true), (7, 4, 7, true)]
        {
            for ranks in [1, 3, 8] {
                dft_mpi::simulate(ranks, |comm| {
                    let mut op = Partitioned {
                        dim,
                        rows: (comm.rank()..dim).step_by(ranks).collect(),
                        clustered,
                    };
                    let rows = op.rows.len();
                    let mut x = RealMatrix::zeros(rows, n).unwrap();
                    for j in 0..n {
                        for (i, &row) in op.rows.iter().enumerate() {
                            x[(i, j)] = if row == j {
                                1.0
                            } else {
                                0.03 * ((row * 7 + j * 13 + 1) as f64).sin()
                            };
                        }
                    }
                    let diag: Vec<_> = op
                        .rows
                        .iter()
                        .map(|&row| op.diagonal(row) + 0.03 * ((row + 1) as f64).sin().powi(2))
                        .collect();
                    for restart in ["off", "on"] {
                        let config = |_mask: usize| RunConfig::from_pairs(&[]).unwrap();
                        let mut options = Options::new(n, 1e-10, n - 1);
                        options.capacity = capacity;
                        options.max_iterations = 40;
                        options.profile = true;
                        options.return_hx = true;
                        work::take_solves();
                        let a = {
                            let _g = config(0).enter();
                            solve(comm, &mut op, &x, &diag, options).unwrap()
                        };
                        let _aw = work::take_solves().pop().unwrap();
                        assert!(a.converged.iter().all(|&v| v));
                        for mask in [1, 2, 4, 8, 14, 15] {
                            // Two solves exercise both cold and reused buffers.
                            for repeat in 0..2 {
                                let b = {
                                    let _g = config(mask).enter();
                                    solve(comm, &mut op, &x, &diag, options).unwrap()
                                };
                                let bw = work::take_solves().pop().unwrap();
                                assert_eq!(a.converged, b.converged);
                                assert_eq!(
                                    a.record.h_columns, b.record.h_columns,
                                    "mask={mask} dim={dim} ranks={ranks}"
                                );
                                assert_eq!(a.record.h_calls, b.record.h_calls);
                                assert_eq!(a.record.communication, b.record.communication);
                                assert_eq!(a.record.iterations.len(), b.record.iterations.len());
                                for (a, b) in a.record.iterations.iter().zip(&b.record.iterations) {
                                    assert_eq!(a.active_bands, b.active_bands);
                                    assert_eq!(a.reduced_order, b.reduced_order);
                                    assert_eq!(a.repairs, b.repairs);
                                    assert_eq!(a.dropped_corrections, b.dropped_corrections);
                                    assert_eq!(a.gemms, b.gemms);
                                }
                                let delta = a
                                    .eigenvalues
                                    .iter()
                                    .zip(&b.eigenvalues)
                                    .map(|(a, b)| (a - b).abs())
                                    .fold(0.0, f64::max);
                                assert!(delta < 1e-12, "delta={delta:e} Ha mask={mask} dim={dim}");
                                if mask & 1 == 0 {
                                    same_bits(&a.eigenvalues, &b.eigenvalues);
                                    same_bits(a.eigenvectors.as_slice(), b.eigenvectors.as_slice());
                                    same_bits(
                                        a.hx.as_ref().unwrap().as_slice(),
                                        b.hx.as_ref().unwrap().as_slice(),
                                    );
                                    assert_eq!(a.ritz_values_ha, b.ritz_values_ha);
                                } else {
                                    assert_eq!(
                                        bw.stages[Stage::Precondition as usize].hypot_calls,
                                        0
                                    );
                                }

                                if comm.rank() == 0 && repeat == 1 {
                                    eprintln!(
                                        "STREAM dim={dim} bands={n} ranks={ranks} restart={restart} mask={mask} H={} iterations={} delta_e_Ha={delta:e} decisions=identical",
                                        a.record.h_columns,
                                        a.record.iterations.len()
                                    );
                                }
                            }
                        }
                    }
                });
            }
        }
    }
    #[test]
    fn stream_check_append_agrees_errors_and_nonfinite_on_one_rank() {
        use dft_numerics::switches::RunConfig;
        struct Failing {
            op: Partitioned,
            calls: usize,
            nan: bool,
        }
        impl DistributedRealOperator for Failing {
            fn global_len(&self) -> usize {
                self.op.global_len()
            }
            fn local_rows(&self) -> &[usize] {
                self.op.local_rows()
            }
            fn apply_block<C: Collective>(
                &mut self,
                c: &C,
                x: &RealMatrix,
                y: &mut RealMatrix,
            ) -> Result<(), SolverError> {
                self.op.apply_block(c, x, y)?;
                self.calls += 1;
                if self.calls == 2 && c.rank() == 1 {
                    if self.nan {
                        y[(0, 0)] = f64::NAN;
                    } else {
                        return Err(SolverError::Application("injected".into()));
                    }
                }
                Ok(())
            }
        }
        for nan in [false, true] {
            dft_mpi::simulate(3, |comm| {
                let _g = RunConfig::from_pairs(&[]).unwrap().enter();
                let rows: Vec<_> = (comm.rank()..12).step_by(3).collect();
                let mut x = RealMatrix::zeros(rows.len(), 3).unwrap();
                for (i, &row) in rows.iter().enumerate() {
                    if row < 3 {
                        x[(i, row)] = 1.0;
                    }
                }
                let diag: Vec<_> = rows.iter().map(|&i| 1.0 + i as f64 * 0.4).collect();
                let mut op = Failing {
                    op: Partitioned {
                        dim: 12,
                        rows,
                        clustered: false,
                    },
                    calls: 0,
                    nan,
                };
                let mut options = Options::new(3, 1e-10, 3);
                options.profile = true;
                work::take_solves();
                assert!(solve(comm, &mut op, &x, &diag, options).is_err());
                assert!(work::take_solves().is_empty());
                // A failed solve must neither publish work nor retain dirty buffers.
                solve(comm, &mut op.op, &x, &diag, options).unwrap();
                assert_eq!(work::take_solves().len(), 1);
            });
        }
    }
    #[test]
    fn stream_extreme_columns_use_hypot_and_fusions_save_counted_reads() {
        use dft_numerics::switches::RunConfig;
        dft_mpi::simulate(1, |comm| {
            let mut op = Tridiagonal((0..12).collect());
            let mut x = RealMatrix::zeros(12, 3).unwrap();
            for j in 0..3 {
                x[(j, j)] = 1.0;
            }
            let mut diag: Vec<_> = (1..=12).map(f64::from).collect();
            let mut options = Options::new(3, 1e-10, 3);
            options.profile = true;
            options.max_iterations = 1;
            for extreme in [false, true] {
                if extreme {
                    diag[0] = -f64::MAX / 8.0;
                }
                work::take_solves();
                let a = {
                    let _g = RunConfig::from_pairs(&[]).unwrap().enter();
                    solve(comm, &mut op, &x, &diag, options).unwrap()
                };
                let aw = work::take_solves().pop().unwrap();
                let b = {
                    let _g = RunConfig::from_pairs(&[]).unwrap().enter();
                    solve(comm, &mut op, &x, &diag, options).unwrap()
                };
                let bw = work::take_solves().pop().unwrap();
                assert_eq!(a.record.h_columns, b.record.h_columns);
                assert_eq!(
                    bw.stages[Stage::Precondition as usize].hypot_calls,
                    if extreme { 36 } else { 0 }
                );
                assert_eq!(
                    bw.stages[Stage::Precondition as usize].sqrt_calls,
                    if extreme { 0 } else { 36 }
                );
                if extreme {
                    same_bits(&a.eigenvalues, &b.eigenvalues);
                }
                let traffic = |w: &work::Work| {
                    [Stage::Normalize, Stage::Append, Stage::ApplyCheck]
                        .iter()
                        .map(|&s| w.stages[s as usize].streamed_bytes)
                        .sum::<u64>()
                };
                assert_eq!(traffic(&aw), traffic(&bw));
            }
        });
    }
}
