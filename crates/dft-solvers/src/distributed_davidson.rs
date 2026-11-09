//! Row-block distributed real Davidson (ADR 0041, M1c). Interface in and out is
//! replicated (full initial vectors, full diagonal, full eigenvectors on every
//! rank); inside, the basis and the products are plane-wave row blocks, the dense
//! algebra is local GEMM plus the deterministic tree sum, and H is applied to
//! column pairs owned by one rank after gathering the correction block. With one
//! rank the caller should use the threads-only `RealDavidson` (bit-identical by
//! construction); this solver is used for two or more ranks.
use crate::{
    BandTolerances, RealLinearOperator, SolverError,
    distributed::{RowBlock, column_norms_squared, gram, project_out, rotate_local},
    real_davidson::{RealEigenSolution, free_electron_initial_guess},
};
use dft_mpi::{Collective, block_counts, block_displacements};
use dft_numerics::{
    linalg::DEPENDENCE_TOLERANCE,
    pool::WorkerPool,
    real::{RealEighWorkspace, RealGemmWorkspace, RealMatrix},
};

fn comm_error(e: dft_mpi::MpiError) -> SolverError {
    SolverError::Application(e.to_string())
}

use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

const STAGES: [&str; 8] = [
    "reduced_gram",
    "reduced_eigen",
    "rotate",
    "residual",
    "orthonormalize",
    "apply_total",
    "apply_operator",
    "other_solve",
];
static STAGE_NANOS: [AtomicU64; 8] = [const { AtomicU64::new(0) }; 8];

/// Wall seconds of this process in each stage of the distributed solver (cheap
/// per-stage accumulation, no I/O in the loop; read once for the run report).
pub fn stage_seconds() -> Vec<(&'static str, f64)> {
    (0..8)
        .map(|i| {
            (
                STAGES[i],
                STAGE_NANOS[i].load(Ordering::Relaxed) as f64 * 1e-9,
            )
        })
        .collect()
}
fn stage_add(i: usize, start: Instant) {
    STAGE_NANOS[i].fetch_add(start.elapsed().as_nanos() as u64, Ordering::Relaxed);
}

struct Layout {
    n: usize,
    rows: RowBlock,
    counts: Vec<usize>,
    displs: Vec<usize>,
}

/// Rows of a replicated matrix owned by this rank, as a row block.
fn local_rows(layout: &Layout, full: &RealMatrix) -> Result<RealMatrix, SolverError> {
    layout.rows.scatter(full)
}

/// Gather the full n x w matrix on every rank from the row blocks (one allgatherv).
fn gather_full<C: Collective>(
    comm: &C,
    layout: &Layout,
    local: &RealMatrix,
) -> Result<RealMatrix, SolverError> {
    let w = local.cols();
    let counts: Vec<usize> = layout.counts.iter().map(|c| c * w).collect();
    let mut flat = vec![0.0; layout.n * w];
    comm.allgatherv_f64(local.as_slice(), &counts, &mut flat)
        .map_err(comm_error)?;
    let mut data = vec![0.0; layout.n * w];
    let mut at = 0;
    for (r, &rc) in layout.counts.iter().enumerate() {
        for j in 0..w {
            let src = &flat[at + j * rc..at + (j + 1) * rc];
            data[j * layout.n + layout.displs[r]..j * layout.n + layout.displs[r] + rc]
                .copy_from_slice(src);
        }
        at += rc * w;
    }
    Ok(RealMatrix::new(layout.n, w, data)?)
}

/// H applied to the columns of a row-block matrix: gather the block, apply H to
/// the column pairs this rank owns (the pair partition, and therefore every
/// column's result, does not depend on the rank count), gather the products and
/// keep the local rows.
fn apply_distributed<C: Collective>(
    comm: &C,
    layout: &Layout,
    op: &dyn RealLinearOperator,
    x_local: &RealMatrix,
) -> Result<RealMatrix, SolverError> {
    let apply_start = Instant::now();
    let w = x_local.cols();
    let n = layout.n;
    let full = gather_full(comm, layout, x_local)?;
    let pairs = w.div_ceil(2);
    let pair_counts = block_counts(pairs, comm.size());
    let pair_displs = block_displacements(&pair_counts);
    let first = (2 * pair_displs[comm.rank()]).min(w);
    let last = (2 * (pair_displs[comm.rank()] + pair_counts[comm.rank()])).min(w);
    let owned = last - first;
    let mut mine = RealMatrix::zeros(n, owned)?;
    if owned > 0 {
        let x = RealMatrix::new(n, owned, full.as_slice()[first * n..last * n].to_vec())?;
        let op_start = Instant::now();
        op.apply_block(&x, &mut mine)?;
        stage_add(6, op_start);
        if mine.rows() != n
            || mine.cols() != owned
            || mine.as_slice().iter().any(|v| !v.is_finite())
        {
            return Err(SolverError::Application(
                "operator returned an invalid block".into(),
            ));
        }
    }
    let col_counts: Vec<usize> = (0..comm.size())
        .map(|r| {
            let f = (2 * pair_displs[r]).min(w);
            let l = (2 * (pair_displs[r] + pair_counts[r])).min(w);
            (l - f) * n
        })
        .collect();
    let mut products = vec![0.0; n * w];
    comm.allgatherv_f64(mine.as_slice(), &col_counts, &mut products)
        .map_err(comm_error)?;
    let out = local_rows(layout, &RealMatrix::new(n, w, products)?);
    stage_add(5, apply_start);
    out
}

/// Sequential modified Gram-Schmidt inside the panel (two passes per column), the
/// distributed counterpart of the serial solver's in-block orthogonalisation:
/// a column is dropped only when its remaining relative norm is at most the
/// dependence tolerance, so near-dependent corrections are kept exactly as in the
/// serial solver (a Gram-matrix Cholesky would discard them). Each reduction is a
/// deterministic tree sum of a short vector.
fn mgs_panel<C: Collective>(comm: &C, x: &RealMatrix) -> Result<RealMatrix, SolverError> {
    let rows = x.rows();
    let mut kept: Vec<f64> = Vec::new(); // column-major, `count` columns of `rows`
    let mut count = 0usize;
    for j in 0..x.cols() {
        let mut v = column_of(x, j).to_vec();
        for _ in 0..2 {
            if count > 0 {
                let mut local = vec![0.0; count];
                for (c, q) in local.iter_mut().zip(kept.chunks(rows)) {
                    *c = q.iter().zip(&v).map(|(a, b)| a * b).sum();
                }
                let mut sums = vec![0.0; count];
                comm.tree_sum_f64(&local, &mut sums).map_err(comm_error)?;
                for (c, q) in sums.iter().zip(kept.chunks(rows)) {
                    for (vi, qi) in v.iter_mut().zip(q) {
                        *vi -= c * qi;
                    }
                }
            }
        }
        let local = [v.iter().map(|a| a * a).sum::<f64>()];
        let mut total = [0.0];
        comm.tree_sum_f64(&local, &mut total).map_err(comm_error)?;
        let len = total[0].sqrt();
        if !len.is_finite() || len <= DEPENDENCE_TOLERANCE {
            continue;
        }
        kept.extend(v.iter().map(|a| a / len));
        count += 1;
    }
    Ok(RealMatrix::new(rows, count, kept)?)
}

fn column_of(x: &RealMatrix, j: usize) -> &[f64] {
    &x.as_slice()[j * x.rows()..(j + 1) * x.rows()]
}

/// Normalise the columns, project twice against `against` (block Gram-Schmidt),
/// drop columns that vanish, then orthonormalise the rest with two CholeskyQR
/// passes. Returns the kept block (possibly with zero columns).
fn orthonormalize<C: Collective>(
    comm: &C,
    ws: &mut RealGemmWorkspace,
    against: Option<&RealMatrix>,
    block: &RealMatrix,
) -> Result<RealMatrix, SolverError> {
    let rows = block.rows();
    let mut x = block.clone();
    let norms = column_norms_squared(comm, &x)?;
    for (j, nn) in norms.iter().enumerate() {
        let len = nn.sqrt();
        let col = &mut x.as_mut_slice()[j * rows..(j + 1) * rows];
        for v in col.iter_mut() {
            *v = if len == 0.0 { 0.0 } else { *v / len };
        }
    }
    if let Some(q) = against.filter(|q| q.cols() > 0) {
        for _ in 0..2 {
            let mut overlap = RealMatrix::zeros(q.cols(), x.cols())?;
            let mut scratch = RealMatrix::zeros(rows, x.cols())?;
            project_out(comm, ws, q, &mut x, &mut overlap, &mut scratch)?;
        }
    }
    // Drop columns that vanished after the projection, renormalise the others.
    let norms = column_norms_squared(comm, &x)?;
    let survivors: Vec<usize> = (0..x.cols())
        .filter(|&j| norms[j].sqrt() > DEPENDENCE_TOLERANCE)
        .collect();
    if survivors.is_empty() {
        return Ok(RealMatrix::zeros(rows, 0)?);
    }
    let mut y = RealMatrix::zeros(rows, survivors.len())?;
    for (c, &j) in survivors.iter().enumerate() {
        let len = norms[j].sqrt();
        for i in 0..rows {
            y[(i, c)] = x[(i, j)] / len;
        }
    }
    let mut y = mgs_panel(comm, &y)?;
    // Renormalising a nearly dependent column amplifies rounding error along the
    // basis; one more projection plus a panel pass (nothing drops now) restores
    // orthogonality to working precision.
    if let Some(q) = against.filter(|q| q.cols() > 0 && y.cols() > 0) {
        let mut overlap = RealMatrix::zeros(q.cols(), y.cols())?;
        let mut scratch = RealMatrix::zeros(rows, y.cols())?;
        project_out(comm, ws, q, &mut y, &mut overlap, &mut scratch)?;
        y = mgs_panel(comm, &y)?;
    }
    Ok(y)
}

fn append_columns(a: &mut RealMatrix, b: &RealMatrix) -> Result<(), SolverError> {
    let mut data = a.as_slice().to_vec();
    data.extend_from_slice(b.as_slice());
    *a = RealMatrix::new(a.rows(), a.cols() + b.cols(), data)?;
    Ok(())
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        (((z ^ (z >> 31)) >> 11) as f64) * (1. / ((1_u64 << 53) as f64)) - 0.5
    }
}

impl crate::RealDavidson {
    /// Distributed solve over the ranks of `comm` (two or more; with one rank use
    /// `solve_policy`/`solve_ritz_changes`). Every rank must call it with the same
    /// replicated arguments and returns the same replicated solution.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_distributed<C: Collective>(
        &self,
        comm: &C,
        op: &dyn RealLinearOperator,
        initial: Option<&RealMatrix>,
        diagonal: &[f64],
        policy: BandTolerances,
        ritz: bool,
        certify: bool,
        capacity_per_band: usize,
    ) -> Result<RealEigenSolution, SolverError> {
        let n = op.dimension();
        if !(2..=12).contains(&capacity_per_band) {
            return Err(SolverError::InvalidParameter(
                "davidson_capacity_per_band must be 2 to 12",
            ));
        }
        if diagonal.len() != n || self.bands == 0 || self.bands > n {
            return Err(SolverError::InvalidParameter("bands or diagonal length"));
        }
        if !self.tolerance.is_finite()
            || self.tolerance <= 0.0
            || self.max_iterations == 0
            || policy.strict_bands == 0
            || policy.strict_bands > self.bands
            || !policy.buffer_tolerance.is_finite()
            || policy.buffer_tolerance < self.tolerance
        {
            return Err(SolverError::InvalidParameter(
                "invalid tolerance, iteration limit or strict prefix",
            ));
        }
        let bands = self.bands;
        let capacity = self.max_subspace.max(bands * capacity_per_band).min(n);
        let ranks = comm.size();
        let rank = comm.rank();
        let counts = block_counts(n, ranks);
        let displs = block_displacements(&counts);
        let layout = Layout {
            n,
            rows: RowBlock::new(n, rank, ranks),
            counts: counts.clone(),
            displs,
        };
        let nl = layout.rows.local_rows();
        let diag_local = &diagonal[layout.rows.offset()..layout.rows.offset() + nl];
        let pool = match op.worker_pool() {
            Some(p) => p,
            None => WorkerPool::new(op.inner_threads())
                .map_err(|e| SolverError::Application(e.to_string()))?,
        };
        let mut ws = RealGemmWorkspace::with_pool(counts[0].max(1), capacity.max(1), pool)?;
        let mut eig = RealEighWorkspace::new(capacity)?;
        let mut rng = Rng(self.seed);

        // Initial basis.
        let start = match initial {
            Some(x) => x.clone(),
            None => free_electron_initial_guess(diagonal, bands)?,
        };
        let mut basis = orthonormalize(comm, &mut ws, None, &local_rows(&layout, &start)?)?;
        while basis.cols() < bands {
            let w = bands - basis.cols();
            let random = RealMatrix::new(n, w, (0..n * w).map(|_| rng.next()).collect())?;
            let extra =
                orthonormalize(comm, &mut ws, Some(&basis), &local_rows(&layout, &random)?)?;
            if extra.cols() == 0 {
                return Err(SolverError::InvalidParameter(
                    "unable to complete initial basis",
                ));
            }
            append_columns(&mut basis, &extra)?;
        }
        let mut applications = 0usize;
        let mut products = apply_distributed(comm, &layout, op, &basis)?;
        applications += basis.cols();

        let threshold = |j: usize| {
            if j < policy.strict_bands - usize::from(ritz) {
                self.tolerance
            } else {
                policy.buffer_tolerance
            }
        };
        let meets_policy = |r: &[f64]| r.iter().enumerate().all(|(j, &x)| x <= threshold(j));
        let meets_strict = |r: &[f64]| {
            r[..policy.strict_bands]
                .iter()
                .all(|&x| x <= self.tolerance)
        };
        let mut history = Vec::with_capacity(self.max_iterations);
        let (mut restarts, mut contractions, mut certificates) = (0usize, 0usize, 0usize);
        let mut previous = vec![f64::NAN; bands];
        let mut reduced_width = 0usize;
        let mut reduced_prev = RealMatrix::zeros(1, 1)?;
        let mut max_basis = 0usize;
        let buffer: Vec<bool> = (0..bands).map(|j| j >= policy.strict_bands).collect();

        for iteration in 1..=self.max_iterations {
            let m = basis.cols();
            max_basis = max_basis.max(m);
            // Reduced matrix V^T HV (replicated), symmetrised from the upper triangle.
            // Incremental, as in the serial solver: only the columns added since the last
            // iteration are formed (m x new products); the rest is reused.
            let t0 = Instant::now();
            let mut reduced = RealMatrix::zeros(m, m)?;
            for c in 0..reduced_width {
                for r in 0..reduced_width {
                    reduced[(r, c)] = reduced_prev[(r, c)];
                }
            }
            let fresh = m - reduced_width;
            let mut fresh_products = RealMatrix::zeros(nl, fresh)?;
            fresh_products
                .as_mut_slice()
                .copy_from_slice(&products.as_slice()[reduced_width * nl..m * nl]);
            let mut block = RealMatrix::zeros(m, fresh)?;
            gram(comm, &mut ws, &basis, &fresh_products, &mut block)?;
            for c in reduced_width..m {
                for r in 0..=c {
                    let v = block[(r, c - reduced_width)];
                    reduced[(r, c)] = v;
                    reduced[(c, r)] = v;
                }
            }
            reduced_width = m;
            stage_add(0, t0);
            reduced_prev = reduced.clone();
            let mut values = vec![0.0; m];
            let mut coefficients = RealMatrix::zeros(m, m)?;
            let t0 = Instant::now();
            eig.solve_into(&reduced, &mut values, &mut coefficients)?;
            stage_add(1, t0);
            let t0 = Instant::now();
            let lowest = RealMatrix::new(m, bands, coefficients.as_slice()[..m * bands].to_vec())?;
            let mut vectors = RealMatrix::zeros(nl, bands)?;
            let mut rotated = RealMatrix::zeros(nl, bands)?;
            rotate_local(&mut ws, &basis, &lowest, &mut vectors)?;
            rotate_local(&mut ws, &products, &lowest, &mut rotated)?;
            stage_add(2, t0);
            let residual_of = |hv: &RealMatrix,
                               v: &RealMatrix,
                               comm: &C|
             -> Result<(RealMatrix, Vec<f64>), SolverError> {
                let mut r = RealMatrix::zeros(nl, bands)?;
                for j in 0..bands {
                    for i in 0..nl {
                        r[(i, j)] = hv[(i, j)] - values[j] * v[(i, j)];
                    }
                }
                let norms = column_norms_squared(comm, &r)?
                    .into_iter()
                    .map(f64::sqrt)
                    .collect();
                Ok((r, norms))
            };
            let t0 = Instant::now();
            let (mut residual_block, mut residuals) = residual_of(&rotated, &vectors, comm)?;
            stage_add(3, t0);
            let mut certified = false;
            let mut trusted = true;
            if !ritz
                && (certify || iteration == self.max_iterations)
                && (meets_policy(&residuals)
                    || (iteration == self.max_iterations && meets_strict(&residuals)))
            {
                certificates += 1;
                let actual = apply_distributed(comm, &layout, op, &vectors)?;
                applications += bands;
                (residual_block, residuals) = residual_of(&actual, &vectors, comm)?;
                certified = true;
                trusted = meets_policy(&residuals);
            }
            let maximum = residuals.iter().copied().fold(0.0, f64::max);
            history.push(maximum);
            let converged: Vec<bool> = (0..bands)
                .map(|j| {
                    if ritz {
                        iteration > 1 && (values[j] - previous[j]).abs() < threshold(j)
                    } else {
                        residuals[j] <= threshold(j)
                    }
                })
                .collect();
            previous.copy_from_slice(&values[..bands]);
            let finished = converged.iter().all(|&b| b)
                || (iteration == self.max_iterations
                    && converged[..policy.strict_bands].iter().all(|&b| b));
            if finished {
                let full = gather_full(comm, &layout, &vectors)?;
                return Ok(RealEigenSolution {
                    eigenvalues: values[..bands].to_vec(),
                    eigenvectors: full,
                    residuals,
                    converged,
                    buffer,
                    residuals_certified: certified,
                    iterations: iteration,
                    history,
                    operator_applications: applications,
                    restarts,
                    contracted_restarts: contractions,
                    certification_attempts: certificates,
                    max_basis_dimension: max_basis,
                });
            }
            if iteration == self.max_iterations {
                return Err(SolverError::NotConverged {
                    iterations: iteration,
                    max_residual: maximum,
                    residuals,
                });
            }
            if ritz && iteration == 1 && basis.cols() == n {
                continue;
            }
            // Corrections for the unconverged roots (row-local preconditioning).
            let selected: Vec<usize> = (0..bands)
                .filter(|&j| !converged[j])
                .take(bands.min(capacity - bands))
                .collect();
            let mut corrections = RealMatrix::zeros(nl, selected.len())?;
            for (k, &j) in selected.iter().enumerate() {
                for i in 0..nl {
                    let delta = diag_local[i] - values[j];
                    let denominator = if delta < 0.0 {
                        -delta.abs().max(0.1)
                    } else {
                        delta.max(0.1)
                    };
                    corrections[(i, k)] = -residual_block[(i, j)] / denominator;
                }
            }
            if !trusted || basis.cols() + corrections.cols() > capacity {
                restarts += 1;
                reduced_width = 0;
                let mut defect = RealMatrix::zeros(bands, bands)?;
                gram(comm, &mut ws, &vectors, &vectors, &mut defect)?;
                let mut worst = 0.0_f64;
                for j in 0..bands {
                    for i in 0..bands {
                        worst = worst.max((defect[(i, j)] - if i == j { 1.0 } else { 0.0 }).abs());
                    }
                }
                if trusted && worst <= 1e-12 {
                    contractions += 1;
                    basis = vectors.clone();
                    products = rotated.clone();
                } else {
                    basis = orthonormalize(comm, &mut ws, None, &vectors)?;
                    products = apply_distributed(comm, &layout, op, &basis)?;
                    applications += basis.cols();
                }
            }
            let t0 = Instant::now();
            let mut extra = orthonormalize(comm, &mut ws, Some(&basis), &corrections)?;
            stage_add(4, t0);
            if extra.cols() == 0 {
                let w = bands.min(capacity - basis.cols());
                let random = RealMatrix::new(n, w, (0..n * w).map(|_| rng.next()).collect())?;
                extra =
                    orthonormalize(comm, &mut ws, Some(&basis), &local_rows(&layout, &random)?)?;
            }
            if extra.cols() == 0 {
                return Err(SolverError::NotConverged {
                    iterations: iteration,
                    max_residual: maximum,
                    residuals,
                });
            }
            let extra_products = apply_distributed(comm, &layout, op, &extra)?;
            applications += extra.cols();
            append_columns(&mut products, &extra_products)?;
            append_columns(&mut basis, &extra)?;
        }
        unreachable!("validated positive iteration limit")
    }
}
