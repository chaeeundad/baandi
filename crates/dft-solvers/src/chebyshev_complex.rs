//! R1 (`docs/roadmap/rpa-crpa-plan.md` 3.2): Chebyshev-filtered subspace
//! iteration with locking for many bands of a fixed complex Hamiltonian on
//! owned plane-wave rows.
//!
//! The block holds `bands` requested columns plus buffer columns that are
//! never certified and are discarded. Each sweep filters the unlocked columns
//! with a scaled Chebyshev polynomial of degree `degree` that damps the
//! interval [largest Ritz value of the block, `upper_bound`] (Zhou, Saad,
//! Tiago, Chelikowsky, J. Comput. Phys. 219, 172 (2006)), orthonormalizes them
//! against the locked columns and among themselves (SVQB, twice), and runs a
//! Rayleigh-Ritz step on the unlocked subspace. Leading contiguous columns
//! whose residual is at most `tolerance` are locked and never touched again.
//!
//! Every inner product is a local GEMM followed by the fixed-association
//! `tree_sum_f64`; results agree bit for bit at a fixed rank count and differ
//! across rank counts only by roundoff. No clocks, I/O or logging. Residuals
//! returned here are the solver's own (from the Rayleigh-Ritz H X); callers
//! certify with a fresh application.
use crate::{SolverError, ref_davidson_complex::DistributedComplexOperator};
use dft_mpi::Collective;
use dft_numerics::linalg::{ComplexMatrix as Matrix, MatrixOp, gemm, hermitian_eigh};
use num_complex::Complex64 as C;

#[derive(Clone, Copy, Debug)]
pub struct Options {
    /// Leading columns that must converge (the requested bands).
    pub bands: usize,
    /// Residual (Ha) for locking a leading column.
    pub tolerance: f64,
    /// Chebyshev polynomial degree per sweep (at least 1).
    pub degree: usize,
    /// Sweeps before reporting failure.
    pub max_iterations: usize,
    /// Upper bound of the operator spectrum (Ha); must not be below it.
    pub upper_bound: f64,
    /// Columns per operator application (the operator's scratch capacity).
    pub chunk: usize,
}

#[derive(Debug)]
pub struct Solution {
    /// The `bands` lowest Ritz values, ascending (Ha).
    pub eigenvalues: Vec<f64>,
    /// Owned rows x `bands`, orthonormal over all ranks.
    pub eigenvectors: Matrix,
    /// Solver residuals ||H x - theta x|| (Ha) at locking time.
    pub residuals: Vec<f64>,
    pub iterations: usize,
    /// Operator columns applied, including the start.
    pub h_columns: u64,
}

fn error(e: impl std::fmt::Display) -> SolverError {
    SolverError::Application(e.to_string())
}

fn reduce<K: Collective>(comm: &K, a: &mut Matrix) -> Result<(), SolverError> {
    let flat: Vec<f64> = a.as_slice().iter().flat_map(|z| [z.re, z.im]).collect();
    let mut sum = vec![0.; flat.len()];
    comm.tree_sum_f64(&flat, &mut sum).map_err(error)?;
    for (z, p) in a.as_mut_slice().iter_mut().zip(sum.as_chunks::<2>().0) {
        *z = C::new(p[0], p[1]);
    }
    Ok(())
}

/// a^H b summed over ranks.
fn gram<K: Collective>(comm: &K, a: &Matrix, b: &Matrix) -> Result<Matrix, SolverError> {
    let mut g = gemm(a, MatrixOp::Adjoint, b, MatrixOp::Identity)?;
    reduce(comm, &mut g)?;
    Ok(g)
}

fn hermitize(a: &mut Matrix) {
    let n = a.rows();
    for j in 0..n {
        a[(j, j)] = C::new(a[(j, j)].re, 0.);
        for i in 0..j {
            let v = (a[(i, j)] + a[(j, i)].conj()) * 0.5;
            a[(i, j)] = v;
            a[(j, i)] = v.conj();
        }
    }
}

fn columns(x: &Matrix, from: usize, to: usize) -> Result<Matrix, SolverError> {
    let mut out = Matrix::zeros(x.rows(), to - from)?;
    for j in from..to {
        out.column_mut(j - from).copy_from_slice(x.column(j));
    }
    Ok(out)
}

fn product(x: &Matrix, c: &Matrix) -> Result<Matrix, SolverError> {
    if x.rows() == 0 {
        return Ok(Matrix::zeros(0, c.cols())?);
    }
    Ok(gemm(x, MatrixOp::Identity, c, MatrixOp::Identity)?)
}

/// H x in column chunks of at most `chunk` columns.
fn apply<K: Collective, O: DistributedComplexOperator>(
    comm: &K,
    op: &mut O,
    x: &Matrix,
    chunk: usize,
    count: &mut u64,
) -> Result<Matrix, SolverError> {
    let mut out = Matrix::zeros(x.rows(), x.cols())?;
    let mut start = 0;
    while start < x.cols() {
        let end = (start + chunk).min(x.cols());
        let block = columns(x, start, end)?;
        let mut hb = Matrix::zeros(x.rows(), end - start)?;
        op.apply_block(comm, &block, &mut hb)?;
        for j in start..end {
            out.column_mut(j).copy_from_slice(hb.column(j - start));
        }
        start = end;
    }
    *count += x.cols() as u64;
    Ok(out)
}

/// SVQB (Stathopoulos and Wu 2002): x <- x D U Lambda^-1/2 with S = x^H x,
/// D its inverse root diagonal and (Lambda, U) the eigenpairs of D S D;
/// eigenvalues below 1e-14 of the largest are clamped (dependent columns
/// become finite noise that the next filter sweep removes).
fn svqb<K: Collective>(comm: &K, x: &mut Matrix) -> Result<(), SolverError> {
    let m = x.cols();
    if m == 0 {
        return Ok(());
    }
    let mut s = gram(comm, x, x)?;
    hermitize(&mut s);
    let d: Vec<f64> = (0..m)
        .map(|i| {
            let v = s[(i, i)].re;
            if v > 0. && v.is_finite() {
                1. / v.sqrt()
            } else {
                1.
            }
        })
        .collect();
    for j in 0..m {
        for i in 0..m {
            s[(i, j)] *= d[i] * d[j];
        }
    }
    let (values, vectors) = hermitian_eigh(&s)?;
    let top = values.last().copied().unwrap_or(1.).max(f64::MIN_POSITIVE);
    let mut t = vectors;
    for (j, &value) in values.iter().enumerate() {
        let scale = 1. / value.max(1e-14 * top).sqrt();
        for i in 0..m {
            t[(i, j)] *= d[i] * scale;
        }
    }
    *x = product(x, &t)?;
    x.validate_finite()?;
    Ok(())
}

/// y <- y - l (l^H y).
fn project_out<K: Collective>(comm: &K, l: &Matrix, y: &mut Matrix) -> Result<(), SolverError> {
    if l.cols() == 0 || y.cols() == 0 {
        return Ok(());
    }
    let c = gram(comm, l, y)?;
    let ly = product(l, &c)?;
    for (a, b) in y.as_mut_slice().iter_mut().zip(ly.as_slice()) {
        *a -= b;
    }
    Ok(())
}

fn residuals<K: Collective>(
    comm: &K,
    x: &Matrix,
    hx: &Matrix,
    theta: &[f64],
) -> Result<Vec<f64>, SolverError> {
    let local: Vec<f64> = (0..x.cols())
        .map(|j| {
            x.column(j)
                .iter()
                .zip(hx.column(j))
                .map(|(a, b)| (b - theta[j] * a).norm_sqr())
                .sum()
        })
        .collect();
    let mut all = vec![0.; local.len()];
    comm.tree_sum_f64(&local, &mut all).map_err(error)?;
    Ok(all.into_iter().map(f64::sqrt).collect())
}

/// Rayleigh-Ritz on orthonormal y: returns ascending Ritz values and rotates
/// y and hy in place.
fn rayleigh_ritz<K: Collective>(
    comm: &K,
    y: &mut Matrix,
    hy: &mut Matrix,
) -> Result<Vec<f64>, SolverError> {
    let mut a = gram(comm, y, hy)?;
    hermitize(&mut a);
    let (values, w) = hermitian_eigh(&a)?;
    *y = product(y, &w)?;
    *hy = product(hy, &w)?;
    Ok(values)
}

/// `start` holds the block (bands plus buffer columns, at most the global
/// basis size); its columns need not be orthonormal but must be linearly
/// independent. All ranks pass the same options.
pub fn solve<K: Collective, O: DistributedComplexOperator>(
    comm: &K,
    op: &mut O,
    start: &Matrix,
    options: Options,
) -> Result<Solution, SolverError> {
    let m = start.cols();
    let Options {
        bands,
        tolerance,
        degree,
        max_iterations,
        upper_bound,
        chunk,
    } = options;
    if bands == 0
        || bands > m
        || m > op.global_len()
        || start.rows() != op.local_rows().len()
        || degree == 0
        || chunk == 0
        || !(tolerance > 0. && tolerance.is_finite())
        || !upper_bound.is_finite()
    {
        return Err(SolverError::InvalidParameter("Chebyshev subspace options"));
    }
    let mut count = 0_u64;
    let mut x = start.clone();
    svqb(comm, &mut x)?;
    svqb(comm, &mut x)?;
    let mut hx = apply(comm, op, &x, chunk, &mut count)?;
    let mut theta = rayleigh_ritz(comm, &mut x, &mut hx)?;
    let mut res = residuals(comm, &x, &hx, &theta)?;
    let mut locked = 0;
    let lock = |locked: &mut usize, res: &[f64]| {
        while *locked < bands && res[*locked] <= tolerance {
            *locked += 1;
        }
    };
    lock(&mut locked, &res);
    let mut iterations = 0;
    while locked < bands {
        if iterations == max_iterations {
            return Err(SolverError::Application(format!(
                "Chebyshev subspace: {locked} of {bands} bands locked after {max_iterations} sweeps; max unlocked residual {:e} Ha",
                res[locked..bands].iter().copied().fold(0., f64::max)
            )));
        }
        iterations += 1;
        let a = theta[m - 1];
        let b = upper_bound;
        if a >= b {
            return Err(SolverError::Application(format!(
                "Chebyshev subspace: Ritz value {a} Ha at or above the spectral bound {b} Ha"
            )));
        }
        // Scale point: the lowest unlocked Ritz value, kept strictly below a.
        let lowest = theta[locked].min(a - 1e-3 * (b - a));
        let e = 0.5 * (b - a);
        let c = 0.5 * (b + a);
        let sigma1 = e / (lowest - c);
        let tau = 2. / sigma1;
        let mut previous = columns(&x, locked, m)?;
        let h = apply(comm, op, &previous, chunk, &mut count)?;
        let mut current = previous.clone();
        for ((y, hv), v) in current
            .as_mut_slice()
            .iter_mut()
            .zip(h.as_slice())
            .zip(previous.as_slice())
        {
            *y = (hv - c * v) * (sigma1 / e);
        }
        let mut sigma = sigma1;
        for _ in 1..degree {
            let sigma2 = 1. / (tau - sigma);
            let h = apply(comm, op, &current, chunk, &mut count)?;
            let mut next = h;
            for ((y, v), p) in next
                .as_mut_slice()
                .iter_mut()
                .zip(current.as_slice())
                .zip(previous.as_slice())
            {
                *y = (*y - c * v) * (2. * sigma2 / e) - (sigma * sigma2) * p;
            }
            previous = std::mem::replace(&mut current, next);
            sigma = sigma2;
        }
        drop(previous);
        let mut y = current;
        let fixed = columns(&x, 0, locked)?;
        project_out(comm, &fixed, &mut y)?;
        svqb(comm, &mut y)?;
        project_out(comm, &fixed, &mut y)?;
        svqb(comm, &mut y)?;
        let mut hy = apply(comm, op, &y, chunk, &mut count)?;
        let values = rayleigh_ritz(comm, &mut y, &mut hy)?;
        let r = residuals(comm, &y, &hy, &values)?;
        for j in locked..m {
            x.column_mut(j).copy_from_slice(y.column(j - locked));
            hx.column_mut(j).copy_from_slice(hy.column(j - locked));
            theta[j] = values[j - locked];
            res[j] = r[j - locked];
        }
        lock(&mut locked, &res);
    }
    Ok(Solution {
        eigenvalues: theta[..bands].to_vec(),
        eigenvectors: columns(&x, 0, bands)?,
        residuals: res[..bands].to_vec(),
        iterations,
        h_columns: count,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use dft_numerics::pool::WorkerPool;
    use std::sync::Arc;

    /// Dense Hermitian test operator on rows distributed round-robin.
    struct Dense {
        h: Matrix,
        rows: Vec<usize>,
    }
    impl DistributedComplexOperator for Dense {
        fn global_len(&self) -> usize {
            self.h.rows()
        }
        fn local_rows(&self) -> &[usize] {
            &self.rows
        }
        fn worker_pool(&self) -> Option<Arc<WorkerPool>> {
            None
        }
        fn apply_block<K: Collective>(
            &mut self,
            comm: &K,
            x: &Matrix,
            out: &mut Matrix,
        ) -> Result<(), SolverError> {
            let n = self.h.rows();
            let mut full = vec![0.; 2 * n * x.cols()];
            for j in 0..x.cols() {
                for (i, &g) in self.rows.iter().enumerate() {
                    full[2 * (j * n + g)] = x[(i, j)].re;
                    full[2 * (j * n + g) + 1] = x[(i, j)].im;
                }
            }
            let mut sum = vec![0.; full.len()];
            comm.tree_sum_f64(&full, &mut sum).map_err(error)?;
            for j in 0..x.cols() {
                for (i, &g) in self.rows.iter().enumerate() {
                    out[(i, j)] = (0..n)
                        .map(|l| {
                            self.h[(g, l)] * C::new(sum[2 * (j * n + l)], sum[2 * (j * n + l) + 1])
                        })
                        .sum();
                }
            }
            Ok(())
        }
    }

    fn matrix(n: usize) -> Matrix {
        let mut h = Matrix::zeros(n, n).unwrap();
        for j in 0..n {
            for i in 0..n {
                h[(i, j)] = if i == j {
                    C::new(0.5 * (i as f64).powf(1.3) - 1.0, 0.)
                } else {
                    let (a, b) = (i.min(j), i.max(j));
                    let z = C::new(
                        0.05 * ((a * 7 + b * 3) as f64).sin(),
                        0.04 * ((a * 5 + b * 11) as f64).cos(),
                    ) / (1. + (b - a) as f64);
                    if i < j { z } else { z.conj() }
                };
            }
        }
        h
    }

    #[test]
    fn chebyshev_matches_dense_and_rank_counts() {
        let n = 60;
        let h = matrix(n);
        let (exact, _) = hermitian_eigh(&h).unwrap();
        let bound = (0..n)
            .map(|i| (0..n).map(|j| h[(i, j)].norm()).sum::<f64>())
            .fold(0., f64::max);
        let bands = 20;
        let block = 26;
        let mut reference: Option<Vec<f64>> = None;
        for ranks in [1, 2, 4] {
            let outputs = dft_mpi::simulate(ranks, |comm| {
                let rows: Vec<usize> = (0..n).filter(|g| g % ranks == comm.rank()).collect();
                let mut start = Matrix::zeros(rows.len(), block).unwrap();
                for j in 0..block {
                    for (i, &g) in rows.iter().enumerate() {
                        start[(i, j)] = C::new(
                            ((g * 31 + j * 17 + 3) as f64).sin(),
                            ((g * 13 + j * 7 + 1) as f64).cos(),
                        );
                    }
                }
                let mut op = Dense { h: h.clone(), rows };
                solve(
                    comm,
                    &mut op,
                    &start,
                    Options {
                        bands,
                        tolerance: 1e-9,
                        degree: 8,
                        max_iterations: 200,
                        upper_bound: bound,
                        chunk: 7,
                    },
                )
                .unwrap()
            });
            for s in outputs {
                let error = s
                    .eigenvalues
                    .iter()
                    .zip(&exact)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0., f64::max);
                assert!(error < 1e-12, "ranks {ranks}: {error:e}");
                assert!(s.residuals.iter().all(|r| *r <= 1e-9));
                match &reference {
                    None => reference = Some(s.eigenvalues.clone()),
                    Some(r) => {
                        let d = r
                            .iter()
                            .zip(&s.eigenvalues)
                            .map(|(a, b)| (a - b).abs())
                            .fold(0., f64::max);
                        assert!(d < 1e-12, "rank-count difference {d:e}");
                    }
                }
            }
        }
    }
}
