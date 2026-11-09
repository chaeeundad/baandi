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

/// Pivot floor of the panel Cholesky (relative to unit-norm columns).
const PANEL_PIVOT: f64 = 1e-16;

fn comm_error(e: dft_mpi::MpiError) -> SolverError {
    SolverError::Application(e.to_string())
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
        op.apply_block(&x, &mut mine)?;
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
    local_rows(layout, &RealMatrix::new(n, w, products)?)
}

/// Inverse of an upper-triangular matrix (small, replicated).
fn upper_inverse(r: &RealMatrix) -> RealMatrix {
    let k = r.rows();
    let mut inv = RealMatrix::zeros(k, k).unwrap();
    for j in 0..k {
        inv[(j, j)] = 1.0 / r[(j, j)];
        for i in (0..j).rev() {
            let mut s = 0.0;
            for p in i + 1..=j {
                s += r[(i, p)] * inv[(p, j)];
            }
            inv[(i, j)] = -s / r[(i, i)];
        }
    }
    inv
}

/// One CholeskyQR pass on a block of unit-ish columns: returns the orthonormalised
/// kept columns (columns whose pivot falls below the floor are dropped).
fn cholesky_qr_pass<C: Collective>(
    comm: &C,
    ws: &mut RealGemmWorkspace,
    x: &RealMatrix,
) -> Result<RealMatrix, SolverError> {
    let w = x.cols();
    let mut g = RealMatrix::zeros(w, w)?;
    gram(comm, ws, x, x, &mut g)?;
    let mut kept: Vec<usize> = Vec::new();
    let mut l = RealMatrix::zeros(w, w)?; // l[(a, b)], a, b indexes into `kept`
    for j in 0..w {
        let mut row = vec![0.0; kept.len()];
        for (b, &kb) in kept.iter().enumerate() {
            let mut s = g[(j, kb)];
            for p in 0..b {
                s -= row[p] * l[(b, p)];
            }
            row[b] = s / l[(b, b)];
        }
        let diag = g[(j, j)] - row.iter().map(|v| v * v).sum::<f64>();
        if diag <= PANEL_PIVOT * g[(j, j)].max(PANEL_PIVOT) || !diag.is_finite() {
            continue;
        }
        let a = kept.len();
        for (b, v) in row.iter().enumerate() {
            l[(a, b)] = *v;
        }
        l[(a, a)] = diag.sqrt();
        kept.push(j);
    }
    let k = kept.len();
    let mut r = RealMatrix::zeros(k, k)?;
    for a in 0..k {
        for b in 0..=a {
            r[(b, a)] = l[(a, b)];
        }
    }
    let rinv = upper_inverse(&r);
    let rows = x.rows();
    let mut selected = RealMatrix::zeros(rows, k)?;
    for (c, &j) in kept.iter().enumerate() {
        selected.as_mut_slice()[c * rows..(c + 1) * rows]
            .copy_from_slice(&x.as_slice()[j * rows..(j + 1) * rows]);
    }
    let mut out = RealMatrix::zeros(rows, k)?;
    rotate_local(ws, &selected, &rinv, &mut out)?;
    Ok(out)
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
    let y = cholesky_qr_pass(comm, ws, &y)?;
    if y.cols() == 0 {
        return Ok(y);
    }
    cholesky_qr_pass(comm, ws, &y)
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
        let mut max_basis = 0usize;
        let buffer: Vec<bool> = (0..bands).map(|j| j >= policy.strict_bands).collect();

        for iteration in 1..=self.max_iterations {
            let m = basis.cols();
            max_basis = max_basis.max(m);
            // Reduced matrix V^T HV (replicated), symmetrised from the upper triangle.
            let mut reduced = RealMatrix::zeros(m, m)?;
            gram(comm, &mut ws, &basis, &products, &mut reduced)?;
            for j in 0..m {
                for i in j + 1..m {
                    let v = reduced[(j, i)];
                    reduced[(i, j)] = v;
                }
            }
            let mut values = vec![0.0; m];
            let mut coefficients = RealMatrix::zeros(m, m)?;
            eig.solve_into(&reduced, &mut values, &mut coefficients)?;
            let lowest = RealMatrix::new(m, bands, coefficients.as_slice()[..m * bands].to_vec())?;
            let mut vectors = RealMatrix::zeros(nl, bands)?;
            let mut rotated = RealMatrix::zeros(nl, bands)?;
            rotate_local(&mut ws, &basis, &lowest, &mut vectors)?;
            rotate_local(&mut ws, &products, &lowest, &mut rotated)?;
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
            let (mut residual_block, mut residuals) = residual_of(&rotated, &vectors, comm)?;
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
            let mut extra = orthonormalize(comm, &mut ws, Some(&basis), &corrections)?;
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
