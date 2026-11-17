//! MPI-OVERLAP-DAVIDSON (`BAANDI_M2_LEAN_ALGEBRA=on`, default off): the serial
//! lean solver's algebra (`lean_davidson.rs`) on the distributed stick rows, for
//! Ritz-change (standard) solves of the workspace path. Design and reduction
//! counts: docs/design/mpi-overlap-davidson.md.
//!
//! Distributed (local rows only): the basis V and its H images HV (corrections
//! appended without any projection), Ritz vectors and residuals of the selected
//! roots. Replicated (identical on every rank, computed from reduced data only):
//! the overlap S = V^T V, its inverse Cholesky factor T (bordered, the serial
//! `extend_factor` with the pivot rule `PIVOT_RATIO_MIN`), H_r = V^T H V, the
//! reduced matrix A = T H_r T^T, its eigenpairs and the coefficients C = T^T Y.
//!
//! Two tree sums per iteration in the normal case (no per-column chain):
//! - the S border (site `lean_border`): the raw corrections E are appended to
//!   V without projection and ONE reduction gives V^T E (m x k, the whole
//!   fresh basis after a contraction), whose diagonal is |E_j|^2, together
//!   with the squared norms of the selected residuals (last column). The new
//!   columns of V and S are scaled to unit norm locally (the serial solver
//!   normalizes before the overlap: rounding class), the factor is extended
//!   with the pivot test, and only the accepted block is applied by H;
//! - the H border (site `reduced_gram`) V^T (H E) at the next iteration, as
//!   on the orthonormal path.
//!
//! A first version also delayed the S border to the H border (one tree sum
//! per iteration) and tested the pivot after H; at the warm si32 solve 15 of
//! 22 appends fail the pivot test (stage 1), so it applied H to 15 percent
//! more columns: rejected (docs/design/mpi-overlap-davidson.md, section 6).
//!
//! Fallbacks (the serial rules, same constants): stage 1 projects the
//! normalized block once against V through the factor (E' = E - V T^T T V^T E,
//! V^T E is the S border already reduced), drops columns with |E'_j| <=
//! `DEPENDENCE_TOLERANCE` (one norm reduction), runs the in-block Gram-Schmidt
//! of the serial-order orthonormalization (its rank floor) and tests again
//! with the next S border; stage 2 rotates V to its full Ritz basis (span
//! unchanged, no H application) unless the basis already is a fresh one
//! (contraction or rebuild) and projects the surviving columns with the
//! classic projection against the basis (random block if nothing survives);
//! a fresh basis that fails the test is rebuilt by Gram-Schmidt and H
//! application (counted; expected zero). Every branch depends on reduced
//! values only, so all ranks take the same branches.
use super::*;
use crate::real_davidson::lean::{LeanScratch, extend_factor, reduced_border, ritz_coefficients};
use dft_numerics::GemmSite;
use dft_numerics::real::RealOp;

// solves, iterations, border reductions, appends, appended columns,
// reprojections, Ritz-basis fallbacks, rebuilds, wasted H columns,
// kept/gram starts taken by the factor test, contractions, R1-b contractions
// that kept directions, kept directions, dropped candidates, R1-b fallbacks
static LEAN_COUNTS: [AtomicU64; 15] = [const { AtomicU64::new(0) }; 15];
// Smallest accepted or rejected pivot ratio seen (f64 bits; starts at +inf).
static LEAN_MIN_PIVOT: AtomicU64 = AtomicU64::new(0x7ff0_0000_0000_0000);

/// Counters of the distributed lean algebra (process-wide, summed over the
/// simulated ranks of a test process).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct LeanAlgebraStats {
    pub solves: u64,
    pub iterations: u64,
    /// `lean_border` tree sums (one per iteration that has new columns).
    pub border_reductions: u64,
    pub appends: u64,
    pub appended_columns: u64,
    /// Stage-1 fallbacks (block re-projected through the factor).
    pub reprojections: u64,
    /// Stage-2 fallbacks (Ritz-basis rotation and classic projection).
    pub fallbacks: u64,
    /// Fresh bases rebuilt by Gram-Schmidt and H application.
    pub rebuilds: u64,
    /// H columns applied to blocks a fallback rejected or a rebuild replaced.
    pub wasted_columns: u64,
    /// Kept or Gram-checked starts accepted by the pivot test.
    pub factored_starts: u64,
    /// Contraction restarts.
    pub contractions: u64,
    /// R1-b (`BAANDI_RESTART_KEEP`): contractions that kept previous
    /// directions, directions kept, candidates dropped by the pivot rule in
    /// coefficient space, and contractions whose new basis failed the factor
    /// test and fell back to X alone (one extra border reduction each).
    pub keep_contractions: u64,
    pub kept_directions: u64,
    pub keep_dropped: u64,
    pub keep_fallbacks: u64,
    /// Smallest pivot ratio (squared pivot over the overlap diagonal) seen.
    pub min_pivot: f64,
}

/// The distributed lean algebra counters of this process.
pub fn lean_algebra_stats() -> LeanAlgebraStats {
    let l = |i: usize| LEAN_COUNTS[i].load(Ordering::Relaxed);
    LeanAlgebraStats {
        solves: l(0),
        iterations: l(1),
        border_reductions: l(2),
        appends: l(3),
        appended_columns: l(4),
        reprojections: l(5),
        fallbacks: l(6),
        rebuilds: l(7),
        wasted_columns: l(8),
        factored_starts: l(9),
        contractions: l(10),
        keep_contractions: l(11),
        kept_directions: l(12),
        keep_dropped: l(13),
        keep_fallbacks: l(14),
        min_pivot: f64::from_bits(LEAN_MIN_PIVOT.load(Ordering::Relaxed)),
    }
}

fn count(i: usize, v: usize) {
    LEAN_COUNTS[i].fetch_add(v as u64, Ordering::Relaxed);
}

fn pivot(ratio: f64) {
    if ratio.is_finite() {
        let _ = LEAN_MIN_PIVOT.try_update(Ordering::Relaxed, Ordering::Relaxed, |b| {
            (ratio < f64::from_bits(b)).then_some(ratio.to_bits())
        });
    }
}

/// Replicated lean scratch of one rank: the serial `LeanScratch` (S, T, A, C),
/// a small GEMM workspace for the replicated products (capacity x capacity
/// shapes, fixed panels: thread-count independent and identical on every
/// rank), and the border reduction buffers.
pub(in crate::distributed_davidson) struct LeanState {
    key: (usize, usize, usize),
    scratch: LeanScratch,
    small: RealGemmWorkspace,
    partial: RealMatrix,
    sum: RealMatrix,
    /// Local squared norms of the pending selected residuals.
    pending: Vec<f64>,
    /// Column scales of the border normalization.
    lens: Vec<f64>,
}

/// Make the rank's lean state fit (capacity, bands, pool); counted growth.
pub(super) fn prepare(
    slot: &mut Option<Box<LeanState>>,
    capacity: usize,
    bands: usize,
    pool: &Arc<WorkerPool>,
) -> Result<(), SolverError> {
    let key = (capacity, bands, Arc::as_ptr(pool) as usize);
    if slot.as_ref().map(|s| s.key) == Some(key) {
        return Ok(());
    }
    let cap = capacity.max(1);
    grown(12 * cap * cap * 8);
    let mut partial = empty();
    let mut sum = empty();
    reserve(&mut partial, cap * (2 * cap + 1))?;
    reserve(&mut sum, cap * (2 * cap + 1))?;
    *slot = Some(Box::new(LeanState {
        key,
        scratch: LeanScratch::new(cap, bands)?,
        small: RealGemmWorkspace::with_pool(cap, cap, pool.clone())?,
        partial,
        sum,
        pending: vec![0.0; bands],
        lens: vec![0.0; cap],
    }));
    Ok(())
}

/// Copy the symmetric matrix defined by the upper triangle of `src` (w x w)
/// into S[:w, :w] and factor it from scratch. Returns whether every pivot
/// passed (the start is then the factored basis).
pub(super) fn factor_start(
    state: &mut LeanState,
    src: &RealMatrix,
    w: usize,
) -> Result<bool, SolverError> {
    let LeanState { scratch, small, .. } = state;
    for j in 0..w {
        for i in 0..=j {
            let v = src[(i, j)];
            scratch.s[(i, j)] = v;
            scratch.s[(j, i)] = v;
        }
    }
    let (ok, ratio) = extend_factor(small, scratch, 0, w)?;
    pivot(ratio);
    if ok {
        count(9, 1);
    }
    Ok(ok)
}

/// Solve context (replicated arguments and closures of `solve_on`).
pub(super) struct Ctx<'a, 'f, C: Collective> {
    pub solver: &'a crate::RealDavidson,
    pub comm: &'a C,
    pub layout: &'a Layout,
    pub apply: &'a mut (dyn FnMut(&RealMatrix, &mut RealMatrix) -> Result<(), SolverError> + 'f),
    pub diag_local: &'a [f64],
    pub policy: BandTolerances,
    pub capacity: usize,
    pub tuning: DistributedTuning,
    pub kept: Option<(&'a mut RealMatrix, &'a mut DistributedKeptStart)>,
    pub gather: bool,
    pub rng: Rng,
    pub panel: usize,
    pub ld: usize,
    pub applications: usize,
    /// Columns of S and T already factored by the start (kept or Gram start).
    pub factored: usize,
    pub use_b5: bool,
}

/// Workspace buffers used by the lean iteration.
pub(super) struct Bufs<'a> {
    pub gemm: &'a mut RealGemmWorkspace,
    pub eig: &'a mut Option<(usize, RealEighWorkspace)>,
    pub b5: &'a mut Option<((usize, usize), Bordered, RealGemmWorkspace)>,
    pub basis: &'a mut RealMatrix,
    pub products: &'a mut RealMatrix,
    pub held: &'a mut Vec<f64>,
    pub reduced: &'a mut RealMatrix,
    pub coefficients: &'a mut RealMatrix,
    pub values: &'a mut Vec<f64>,
    pub vectors: &'a mut RealMatrix,
    pub rotated: &'a mut RealMatrix,
    pub residual: &'a mut RealMatrix,
    pub corrections: &'a mut RealMatrix,
    pub extra: &'a mut RealMatrix,
    pub extra_products: &'a mut RealMatrix,
    pub norms: &'a mut Vec<f64>,
    pub norm_partial: &'a mut Vec<f64>,
    pub converged: &'a mut Vec<bool>,
    pub selected: &'a mut Vec<usize>,
    pub sel_values: &'a mut Vec<f64>,
    pub sel_norms: &'a mut Vec<f64>,
    pub band_residuals: &'a mut Vec<f64>,
    pub ortho: &'a mut OrthoScratch,
    pub lean: &'a mut LeanState,
}

/// Residual columns hv - lambda v (column `slot` from input column
/// `map[slot]`, identity without a map) and their LOCAL squared norms (no
/// reduction: they ride with the next border).
fn residual_local(
    hv: &RealMatrix,
    v: &RealMatrix,
    map: Option<&[usize]>,
    lambdas: &[f64],
    out: &mut RealMatrix,
    partial: &mut [f64],
) -> Result<(), SolverError> {
    let (nl, cols) = (v.rows(), lambdas.len());
    fit(out, nl, cols)?;
    for (j, &lambda) in lambdas.iter().enumerate() {
        let c = map.map_or(j, |m| m[j]);
        let (h, x) = (column_of(hv, c), column_of(v, c));
        let col = &mut out.as_mut_slice()[j * nl..(j + 1) * nl];
        for ((r, &h), &x) in col.iter_mut().zip(h).zip(x) {
            *r = h - lambda * x;
        }
        partial[j] = col.iter().map(|a| a * a).sum::<f64>();
    }
    Ok(())
}

/// Keep the leading `cols` columns of a local-row matrix.
fn truncate(m: &mut RealMatrix, cols: usize) -> Result<(), SolverError> {
    let rows = m.rows();
    m.reshape_within_storage(rows, cols)?;
    Ok(())
}

/// Copy `src` (all columns) into `dst` as its whole content.
fn assign(dst: &mut RealMatrix, src: &RealMatrix) -> Result<(), SolverError> {
    fit(dst, src.rows(), src.cols())?;
    dst.as_mut_slice().copy_from_slice(src.as_slice());
    Ok(())
}

/// One border reduction: S[:, s_from..m] = V^T V[:, s_from..m] and, with
/// `h_from`, H_r[:, h_from..m] = V^T (HV)[:, h_from..m] (local GEMMs into one
/// partial matrix) plus the `pending` local squared residual norms in the
/// last column; one tree sum. S and H_r are written symmetric (the upper
/// triangle defines them). The reduced squared norms go to `norms_out`.
#[allow(clippy::too_many_arguments)]
fn border<C: Collective>(
    comm: &C,
    gemm: &mut RealGemmWorkspace,
    basis: &RealMatrix,
    products: &RealMatrix,
    s_from: usize,
    h_from: Option<usize>,
    pending: &[f64],
    partial: &mut RealMatrix,
    sum: &mut RealMatrix,
    s: &mut RealMatrix,
    held: &mut [f64],
    ld: usize,
    norms_out: &mut Vec<f64>,
) -> Result<(), SolverError> {
    let (nl, m) = (basis.rows(), basis.cols());
    let ks = m - s_from;
    let kh = h_from.map_or(0, |h| m - h);
    let cols = ks + kh + 1;
    fit(partial, m, cols)?;
    partial.as_mut_slice().fill(0.0);
    if nl > 0 && m > 0 {
        let v = &basis.as_slice()[..nl * m];
        if ks > 0 {
            gemm.multiply_slices_into(
                GemmSite::Other,
                v,
                nl,
                m,
                nl,
                RealOp::Transpose,
                &basis.as_slice()[s_from * nl..m * nl],
                ks,
                nl,
                partial.as_mut().subcols_mut(0, ks),
            )?;
        }
        if let Some(h) = h_from.filter(|_| kh > 0) {
            gemm.multiply_slices_into(
                GemmSite::Other,
                v,
                nl,
                m,
                nl,
                RealOp::Transpose,
                &products.as_slice()[h * nl..m * nl],
                kh,
                nl,
                partial.as_mut().subcols_mut(ks, kh),
            )?;
        }
    }
    let last = (ks + kh) * m;
    partial.as_mut_slice()[last..last + pending.len()].copy_from_slice(pending);
    fit(sum, m, cols)?;
    site(if ks > 0 || !pending.is_empty() {
        Site::LeanBorder
    } else {
        Site::ReducedGram
    });
    comm.tree_sum_f64(partial.as_slice(), sum.as_mut_slice())
        .map_err(comm_error)?;
    for j in s_from..m {
        for i in 0..=j {
            let v = sum[(i, j - s_from)];
            s[(i, j)] = v;
            s[(j, i)] = v;
        }
    }
    if let Some(h) = h_from {
        for c in h..m {
            for r in 0..=c {
                let v = sum[(r, ks + c - h)];
                held[r + c * ld] = v;
                held[c + r * ld] = v;
            }
        }
    }
    norms_out.clear();
    norms_out.extend((0..pending.len()).map(|slot| sum[(slot, ks + kh)]));
    count(2, 1);
    Ok(())
}

/// Rebuild of a fresh basis that failed the pivot test (expected never):
/// Gram-Schmidt of the whole basis (serial-order orthonormalization with its
/// rank floor) and H on the result. Returns the applied column count.
#[allow(clippy::too_many_arguments)]
fn rebuild<C: Collective>(
    comm: &C,
    apply: &mut dyn FnMut(&RealMatrix, &mut RealMatrix) -> Result<(), SolverError>,
    tuning: &DistributedTuning,
    panel: usize,
    basis: &mut RealMatrix,
    products: &mut RealMatrix,
    extra: &mut RealMatrix,
    gemm: &mut RealGemmWorkspace,
    ortho: &mut OrthoScratch,
    wasted: usize,
) -> Result<usize, SolverError> {
    count(7, 1);
    count(8, wasted);
    ortho_into(comm, gemm, None, basis, tuning, panel, ortho, extra)?;
    if extra.cols() == 0 {
        return Err(SolverError::Application(
            "distributed lean Davidson: rebuilt basis is empty".into(),
        ));
    }
    assign(basis, extra)?;
    let nl = basis.rows();
    fit(products, nl, basis.cols())?;
    apply(basis, products)?;
    Ok(basis.cols())
}

pub(super) fn run<C: Collective>(
    cx: Ctx<'_, '_, C>,
    b: Bufs<'_>,
) -> Result<RealEigenSolution, SolverError> {
    let Ctx {
        solver,
        comm,
        layout,
        apply,
        diag_local,
        policy,
        capacity,
        tuning,
        mut kept,
        gather,
        mut rng,
        panel,
        ld,
        mut applications,
        factored: start_factored,
        use_b5,
    } = cx;
    let Bufs {
        gemm,
        eig,
        b5,
        basis,
        products,
        held,
        reduced,
        coefficients,
        values,
        vectors,
        rotated,
        residual,
        corrections,
        extra,
        extra_products,
        norms,
        norm_partial,
        converged,
        selected,
        sel_values,
        sel_norms,
        band_residuals,
        ortho,
        lean,
    } = b;
    let LeanState {
        scratch: ls,
        small,
        partial,
        sum,
        pending: pending_local,
        lens,
        ..
    } = lean;
    let (n, nl, bands) = (layout.n, layout.nl(), solver.bands);
    let wide = tuning.wide_gemm;
    count(0, 1);
    // R1-b: previous directions kept at a contraction (0 = the plain restart).
    let keep_r = tuning.restart_keep.directions(bands);
    if keep_r > 0 {
        ls.keep_prepare(bands);
    }
    // Kept direction columns of the current fresh contraction basis [X | P]
    // (until its factor test).
    let mut p_cols = 0usize;
    let threshold = |j: usize| {
        // Ritz-change thresholds: the last mandatory root uses the buffer one.
        if j + 1 < policy.strict_bands {
            solver.tolerance
        } else {
            policy.buffer_tolerance
        }
    };
    let precondition = |d: f64, lambda: f64, r: f64| {
        let delta = d - lambda;
        let denominator = if delta < 0.0 {
            -delta.abs().max(0.1)
        } else {
            delta.max(0.1)
        };
        -r / denominator
    };
    let buffer: Vec<bool> = (0..bands).map(|j| j >= policy.strict_bands).collect();
    let mut history = Vec::with_capacity(solver.max_iterations);
    let (mut restarts, mut contractions) = (0usize, 0usize);
    let mut previous = vec![f64::NAN; bands];
    let mut max_basis = 0usize;
    band_residuals.clear();
    band_residuals.resize(bands, 0.0);
    // Columns of S/T factored, of H_r formed (in `held`), of A bordered.
    let mut factored = start_factored.min(basis.cols());
    let mut h_width = 0usize;
    let mut a_width = 0usize;

    for iteration in 1..=solver.max_iterations {
        call_probe(iteration);
        count(1, 1);
        // 1. H border of the columns applied since the last iteration (and S
        // of a fresh basis: start or rebuild), one tree sum.
        let mut rebuilds_here = 0usize;
        loop {
            let m = basis.cols();
            if factored < m || h_width < m {
                let t0 = dft_numerics::stage_clock::start();
                border(
                    comm,
                    gemm,
                    basis,
                    products,
                    factored,
                    Some(h_width),
                    &[],
                    partial,
                    sum,
                    &mut ls.s,
                    held,
                    ld,
                    norms,
                )?;
                h_width = m;
                stage_add(0, t0);
            }
            if factored < m {
                let t0 = dft_numerics::stage_clock::start();
                let (ok, ratio) = extend_factor(small, ls, factored, m)?;
                pivot(ratio);
                stage_add(1, t0);
                if !ok {
                    rebuilds_here += 1;
                    if rebuilds_here > 1 {
                        return Err(SolverError::Application(
                            "distributed lean Davidson: orthonormalized basis failed the overlap pivot test"
                                .into(),
                        ));
                    }
                    let t0 = dft_numerics::stage_clock::start();
                    applications += rebuild(
                        comm, apply, &tuning, panel, basis, products, extra, gemm, ortho, m,
                    )?;
                    ls.keep_invalidate();
                    (factored, h_width, a_width) = (0, 0, 0);
                    stage_add(4, t0);
                    continue;
                }
                factored = m;
            }
            break;
        }
        let m = basis.cols();
        max_basis = max_basis.max(m);
        // 2. Reduced problem A = T H_r T^T (border only), eigenpairs, C = T^T Y.
        let t0 = dft_numerics::stage_clock::start();
        fit(reduced, m, m)?;
        for c in 0..m {
            reduced.as_mut_slice()[c * m..(c + 1) * m].copy_from_slice(&held[c * ld..c * ld + m]);
        }
        reduced_border(small, reduced, ls, a_width, m)?;
        let prefix = a_width;
        a_width = m;
        fit(coefficients, m, m)?;
        record_solve(m, m - prefix);
        if let Some((_, eig)) = eig.as_mut() {
            solve_reduced(
                eig,
                if use_b5 {
                    b5.as_mut().map(|(_, b, g)| (b, g))
                } else {
                    None
                },
                &ls.a,
                prefix,
                &mut values[..m],
                coefficients,
            )?;
        }
        stage_add(1, t0);
        converged.clear();
        converged.extend(
            (0..bands).map(|j| iteration > 1 && (values[j] - previous[j]).abs() < threshold(j)),
        );
        previous.copy_from_slice(&values[..bands]);
        let finished = converged.iter().all(|&b| b)
            || (iteration == solver.max_iterations
                && converged[..policy.strict_bands].iter().all(|&b| b));
        if finished || iteration == solver.max_iterations {
            // Ritz vectors, H images and residual norms of all bands, once.
            let t0 = dft_numerics::stage_clock::start();
            ritz_coefficients(small, ls, coefficients, m, bands, None)?;
            rotate_ritz(
                gemm,
                basis,
                products,
                &ls.c.as_slice()[..m * bands],
                bands,
                vectors,
                rotated,
                wide,
            )?;
            stage_add(2, t0);
            let t0 = dft_numerics::stage_clock::start();
            residual_into(
                comm,
                rotated,
                vectors,
                &values[..bands],
                true,
                residual,
                norms,
                norm_partial,
            )?;
            stage_add(3, t0);
            let maximum = norms.iter().copied().fold(0.0, f64::max);
            history.push(maximum);
            if !finished {
                return Err(SolverError::NotConverged {
                    iterations: iteration,
                    max_residual: maximum,
                    residuals: norms.clone(),
                });
            }
            let t0 = dft_numerics::stage_clock::start();
            keep_return(kept.as_mut(), vectors, rotated, true, nl, bands)?;
            let full = if gather {
                gather_full(comm, layout, vectors)?
            } else {
                RealMatrix::new(0, bands, Vec::new())?
            };
            stage_add(9, t0);
            return Ok(RealEigenSolution {
                eigenvalues: values[..bands].to_vec(),
                eigenvectors: full,
                residuals: norms.clone(),
                converged: converged.clone(),
                buffer,
                residuals_certified: false,
                iterations: iteration,
                history,
                operator_applications: applications,
                restarts,
                contracted_restarts: contractions,
                certification_attempts: 0,
                max_basis_dimension: max_basis,
            });
        }
        if iteration == 1 && m == n {
            // Complete space: the next iteration confirms the exact pairs.
            history.push(0.0);
            continue;
        }
        // 3. Residuals of the unconverged roots only (norms reduced below).
        selected.clear();
        selected.extend(
            (0..bands)
                .filter(|&j| !converged[j])
                .take(bands.min(capacity - bands)),
        );
        let k = selected.len();
        if k == 0 {
            return Err(SolverError::Application(
                "distributed lean Davidson: no root to correct below a full capacity".into(),
            ));
        }
        sel_values.clear();
        sel_values.extend(selected.iter().map(|&j| values[j]));
        let contracted = m + k > capacity;
        let t0 = dft_numerics::stage_clock::start();
        if contracted {
            ritz_coefficients(small, ls, coefficients, m, bands, None)?;
            // R1-b: kept direction coefficients appended to C (replicated
            // arithmetic on the reduced S; no reduction, no H application).
            let (kept_now, dropped) = if keep_r > 0 {
                ls.keep_directions(m, bands, selected, keep_r.min(capacity - bands - k))?
            } else {
                (0, 0)
            };
            count(13, dropped);
            p_cols = kept_now;
            let cols = bands + kept_now;
            rotate_ritz(
                gemm,
                basis,
                products,
                &ls.c.as_slice()[..m * cols],
                cols,
                vectors,
                rotated,
                wide,
            )?;
        } else {
            ritz_coefficients(small, ls, coefficients, m, k, Some(&selected[..]))?;
            if keep_r > 0 {
                // R1-b history: this iteration's selected Ritz coefficients.
                ls.keep_record(m, selected);
            }
            rotate_ritz(
                gemm,
                basis,
                products,
                &ls.c.as_slice()[..m * k],
                k,
                vectors,
                rotated,
                wide,
            )?;
        }
        stage_add(2, t0);
        let t0 = dft_numerics::stage_clock::start();
        residual_local(
            rotated,
            vectors,
            contracted.then_some(&selected[..]),
            sel_values,
            residual,
            &mut pending_local[..k],
        )?;
        fit(corrections, nl, k)?;
        for (slot, &j) in selected.iter().enumerate() {
            let r = column_of(residual, slot);
            for ((c, &d), &r) in corrections.as_mut_slice()[slot * nl..(slot + 1) * nl]
                .iter_mut()
                .zip(diag_local)
                .zip(r)
            {
                *c = precondition(d, values[j], r);
            }
        }
        stage_add(3, t0);
        // `coefficients` holds the full Ritz decomposition of the basis unless
        // a contraction replaces the basis by X below.
        let mut ritz_base = true;
        if contracted {
            // Contraction restart: the basis becomes X (S-orthonormal Ritz
            // vectors), factored afresh in the border below.
            restarts += 1;
            contractions += 1;
            count(10, 1);
            assign(basis, vectors)?;
            assign(products, rotated)?;
            (factored, h_width, a_width) = (0, 0, 0);
            ritz_base = false;
            if keep_r > 0 {
                // The history of the next iteration: the columns of X.
                ls.keep_identity(bands, selected);
            }
        }
        // 4. Raw corrections appended without projection. ONE tree sum gives
        // the S border (the whole fresh part after a contraction), |E_j|^2 and
        // the selected residual norms; normalization, pivot test and the
        // fallbacks run before H is applied to the accepted block.
        let mut append_base = basis.cols();
        append_in_place(basis, corrections)?;
        count(3, 1);
        let mut stage = 0u8;
        let mut pending = k;
        let mut rebuilt = false;
        loop {
            let m = basis.cols();
            let t0 = dft_numerics::stage_clock::start();
            border(
                comm,
                gemm,
                basis,
                products,
                factored,
                None,
                &pending_local[..pending],
                partial,
                sum,
                &mut ls.s,
                held,
                ld,
                norms,
            )?;
            if pending > 0 {
                sel_norms.clear();
                sel_norms.extend(norms.iter().map(|x| x.sqrt()));
                for (slot, &j) in selected.iter().enumerate() {
                    band_residuals[j] = sel_norms[slot];
                }
                history.push(sel_norms.iter().copied().fold(0.0, f64::max));
                pending = 0;
            }
            if stage == 0 {
                // Unit-norm corrections: scale the new columns of V and S by
                // 1 / |E_j| (|E_j|^2 is the reduced S diagonal).
                lens.clear();
                lens.extend((0..m).map(|j| {
                    let len = ls.s[(j, j)].sqrt();
                    if j >= append_base && len > 0.0 && len.is_finite() {
                        len
                    } else {
                        1.0
                    }
                }));
                for j in append_base..m {
                    if lens[j] != 1.0 {
                        for x in &mut basis.as_mut_slice()[j * nl..(j + 1) * nl] {
                            *x /= lens[j];
                        }
                    }
                    for i in 0..=j {
                        let v = ls.s[(i, j)] / (lens[i] * lens[j]);
                        ls.s[(i, j)] = v;
                        ls.s[(j, i)] = v;
                    }
                }
            }
            stage_add(0, t0);
            let t0 = dft_numerics::stage_clock::start();
            if p_cols > 0 {
                // R1-b: the fresh basis is [X | P]; factor X, then P (pivot
                // rule). Any failure falls back to the plain restart from X
                // alone: P is cut, the raw corrections are appended again and
                // the border is repeated (one extra reduction, counted).
                let (x_ok, _) = extend_factor(small, ls, 0, bands)?;
                let p_ok = x_ok && {
                    let (ok, ratio) = extend_factor(small, ls, bands, append_base)?;
                    pivot(ratio);
                    ok
                };
                let kept_now = std::mem::replace(&mut p_cols, 0);
                if p_ok {
                    factored = append_base;
                    count(11, 1);
                    count(12, kept_now);
                } else {
                    count(14, 1);
                    truncate(basis, bands)?;
                    truncate(products, bands)?;
                    append_base = bands;
                    append_in_place(basis, corrections)?;
                    stage_add(1, t0);
                    continue;
                }
            }
            if factored < append_base {
                // The fresh contraction basis X.
                let (ok, ratio) = extend_factor(small, ls, factored, append_base)?;
                pivot(ratio);
                if !ok {
                    stage_add(1, t0);
                    let t0 = dft_numerics::stage_clock::start();
                    applications += rebuild(
                        comm,
                        apply,
                        &tuning,
                        panel,
                        basis,
                        products,
                        extra,
                        gemm,
                        ortho,
                        append_base,
                    )?;
                    ls.keep_invalidate();
                    (factored, h_width, a_width) = (0, 0, 0);
                    rebuilt = true;
                    stage_add(4, t0);
                    break;
                }
                factored = append_base;
            }
            let (ok, ratio) = extend_factor(small, ls, factored, m)?;
            pivot(ratio);
            stage_add(1, t0);
            if ok {
                factored = m;
                break;
            }
            // The appended block failed the pivot test (H not applied yet).
            let t0 = dft_numerics::stage_clock::start();
            if stage == 2 {
                return Err(SolverError::Application(
                    "distributed lean Davidson: projected corrections failed the overlap pivot test"
                        .into(),
                ));
            }
            let base = append_base;
            let kk = m - base;
            fit(extra, nl, kk)?;
            extra
                .as_mut_slice()
                .copy_from_slice(&basis.as_slice()[base * nl..m * nl]);
            truncate(basis, base)?;
            let mut project = stage == 1;
            if stage == 0 {
                // Stage 1: E' = E - V G with G = T^T T S[:base, base..m].
                count(5, 1);
                let cap = ls.cap;
                ls.z.reshape_within_storage(base, kk)?;
                ls.p.reshape_within_storage(base, kk)?;
                if base > 0 {
                    small.multiply_slices_into(
                        GemmSite::Other,
                        ls.t.as_slice(),
                        base,
                        base,
                        cap,
                        RealOp::Identity,
                        &ls.s.as_slice()[base * cap..],
                        kk,
                        cap,
                        ls.z.as_mut(),
                    )?;
                    small.multiply_slices_into(
                        GemmSite::Other,
                        ls.t.as_slice(),
                        base,
                        base,
                        cap,
                        RealOp::Transpose,
                        ls.z.as_slice(),
                        kk,
                        base,
                        ls.p.as_mut(),
                    )?;
                }
                fit(vectors, nl, kk)?;
                rotate_slices_with(
                    gemm,
                    &basis.as_slice()[..nl * base],
                    nl,
                    base,
                    &ls.p.as_slice()[..base * kk],
                    kk,
                    vectors,
                    wide,
                )?;
                fit(corrections, nl, kk)?;
                for j in 0..kk {
                    let dst = &mut corrections.as_mut_slice()[j * nl..(j + 1) * nl];
                    for ((d, &e), &p) in dst
                        .iter_mut()
                        .zip(column_of(extra, j))
                        .zip(column_of(vectors, j))
                    {
                        *d = e - p;
                    }
                }
                site(Site::OrthoNormalize);
                column_norms_squared_into(comm, corrections, norm_partial, norms)?;
                let mut kept_cols = 0;
                for (j, &n2) in norms.iter().enumerate().take(kk) {
                    let len = n2.sqrt();
                    if !len.is_finite() {
                        return Err(SolverError::Application(
                            "distributed lean Davidson: non-finite re-projected correction".into(),
                        ));
                    }
                    if len > DEPENDENCE_TOLERANCE {
                        if kept_cols != j {
                            corrections
                                .as_mut_slice()
                                .copy_within(j * nl..(j + 1) * nl, kept_cols * nl);
                        }
                        kept_cols += 1;
                    }
                }
                truncate(corrections, kept_cols)?;
                if kept_cols > 0 {
                    ortho_into(comm, gemm, None, corrections, &tuning, panel, ortho, extra)?;
                } else {
                    fit(extra, nl, 0)?;
                }
                if extra.cols() > 0 {
                    stage = 1;
                } else {
                    project = true;
                }
            }
            if project {
                // Stage 2: rotate V to its Ritz basis (unless it is the fresh
                // contraction basis), then the classic projection of the
                // surviving re-projected corrections.
                count(6, 1);
                if ritz_base {
                    ritz_coefficients(small, ls, coefficients, base, base, None)?;
                    rotate_ritz(
                        gemm,
                        basis,
                        products,
                        &ls.c.as_slice()[..base * base],
                        base,
                        vectors,
                        rotated,
                        wide,
                    )?;
                    assign(basis, vectors)?;
                    assign(products, rotated)?;
                    (factored, h_width, a_width) = (0, 0, 0);
                    ritz_base = false;
                    if keep_r > 0 {
                        // R1-b: root j is column j of the full Ritz basis.
                        ls.keep_identity(base, selected);
                    }
                }
                ortho_into(
                    comm,
                    gemm,
                    Some(basis),
                    corrections,
                    &tuning,
                    panel,
                    ortho,
                    extra,
                )?;
                if extra.cols() == 0 {
                    let w = bands.min(capacity - base);
                    let random = RealMatrix::new(n, w, (0..n * w).map(|_| rng.next()).collect())?;
                    ortho_into(
                        comm,
                        gemm,
                        Some(basis),
                        &local_rows(layout, &random)?,
                        &tuning,
                        panel,
                        ortho,
                        extra,
                    )?;
                }
                if extra.cols() == 0 {
                    stage_add(4, t0);
                    return Err(SolverError::NotConverged {
                        iterations: iteration,
                        max_residual: history.last().copied().unwrap_or(f64::NAN),
                        residuals: band_residuals.clone(),
                    });
                }
                stage = 2;
            }
            // After a stage-2 rotation the whole basis is fresh: its S comes
            // with the next border (factored = 0).
            append_base = basis.cols();
            append_in_place(basis, extra)?;
            stage_add(4, t0);
        }
        if !rebuilt {
            // 5. H on the accepted block only.
            let m = basis.cols();
            let kk = m - append_base;
            fit(extra, nl, kk)?;
            extra
                .as_mut_slice()
                .copy_from_slice(&basis.as_slice()[append_base * nl..m * nl]);
            fit(extra_products, nl, kk)?;
            apply(extra, extra_products)?;
            applications += kk;
            append_in_place(products, extra_products)?;
            count(4, kk);
        }
    }
    unreachable!("validated positive iteration limit")
}
