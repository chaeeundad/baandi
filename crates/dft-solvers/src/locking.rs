//! Experimental prefix deflation. The legacy solver remains the default.
use super::*;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Locking {
    #[default]
    Off,
    Auto,
}

/// Deterministic solve-local diagnostics, with no clocks or I/O.
#[derive(Debug, Clone, PartialEq)]
pub struct LockingReport {
    pub locked_roots: usize,
    pub unlocks: usize,
    /// (iteration, sorted root, true residual). Reserved before iteration.
    pub failures: Vec<(usize, usize, f64)>,
}

fn buffer(n: usize, capacity: usize) -> Result<ComplexMatrix, LinalgError> {
    let mut result = ComplexMatrix::zeros(n, 0)?;
    result.reserve_columns(capacity)?;
    Ok(result)
}

impl Davidson {
    /// Hard-lock only complete, lowest clusters. Every successful return uses
    /// an all-root Rayleigh–Ritz rotation and a fresh true-residual certificate.
    pub fn solve_locked(
        &self,
        operator: &dyn LinearOperator,
        initial: Option<&ComplexMatrix>,
        preconditioner: Option<&dyn Preconditioner>,
        policy: BandTolerances,
    ) -> Result<EigenSolution, SolverError> {
        let mut trace = if replay::ENABLED {
            replay::Trace::new(self.bands, self.max_iterations)
        } else {
            replay::Trace::default()
        };
        let result = self.solve_locked_inner(operator, initial, preconditioner, policy, &mut trace);
        if replay::ENABLED {
            let dense_values = replay::dense_oracle(operator, self.bands)?;
            replay::store(replay::Record {
                dense_values,
                traces: vec![(
                    self.tolerance,
                    policy.strict_bands,
                    policy.buffer_tolerance,
                    result.is_ok(),
                    trace,
                )],
            });
        }
        result
    }

    fn solve_locked_inner(
        &self,
        operator: &dyn LinearOperator,
        initial: Option<&ComplexMatrix>,
        preconditioner: Option<&dyn Preconditioner>,
        policy: BandTolerances,
        trace: &mut replay::Trace,
    ) -> Result<EigenSolution, SolverError> {
        let n = operator.dimension();
        let m = self.bands;
        validate_bands(n, m)?;
        let capacity = self.max_subspace.min(n);
        if capacity < m
            || (capacity == m && m < n)
            || self.block_size == 0
            || self.max_iterations == 0
            || !self.tolerance.is_finite()
            || self.tolerance <= 0.0
            || policy.strict_bands == 0
            || policy.strict_bands > m
            || !policy.buffer_tolerance.is_finite()
            || policy.buffer_tolerance < self.tolerance
        {
            return Err(SolverError::InvalidParameter(
                "invalid locking solver configuration",
            ));
        }
        let threshold = |j: usize| {
            if j < policy.strict_bands {
                self.tolerance
            } else {
                policy.buffer_tolerance
            }
        };
        // Source-only ablation, never a runtime/default tolerance change.
        let margin_divisor = option_env!("P24A_MARGIN_DIVISOR")
            .map_or(Ok(1.0), str::parse::<f64>)
            .map_err(|_| SolverError::InvalidParameter("invalid margin divisor"))?;
        if !margin_divisor.is_finite() || margin_divisor < 1.0 {
            return Err(SolverError::InvalidParameter("margin divisor must be >= 1"));
        }
        let margin = 1.0 / ((m as f64).sqrt() * margin_divisor);
        let mut gemm = match operator.worker_pool() {
            Some(pool) => GemmWorkspace::with_pool(n, capacity, operator.inner_threads(), pool)?,
            None => GemmWorkspace::new(n, capacity, operator.inner_threads())?,
        };
        let mut projection = ProjectionWorkspace::new(n, capacity)?;
        let mut rng = SplitMix64(self.seed);
        let mut basis = buffer(n, capacity)?;
        if let Some(initial) = initial {
            if initial.rows() != n || initial.cols() > capacity {
                return Err(LinalgError::Dimension.into());
            }
            basis.replace_columns(&projection.orthonormalize(initial, None, &mut gemm)?.vectors)?;
        }
        while basis.cols() < m {
            let random = rng.block(n, m - basis.cols())?;
            let extra = projection
                .orthonormalize(&random, Some(&basis), &mut gemm)?
                .vectors;
            if extra.cols() == 0 {
                return Err(SolverError::InvalidParameter("deficient initial space"));
            }
            basis.append_columns(&extra)?;
        }
        let mut hb = buffer(n, capacity)?;
        hb.resize_columns(basis.cols())?;
        operator.apply_block(&basis, &mut hb)?;
        hb.reserve_columns(capacity)?;
        let mut applications = basis.cols();
        let mut work = work_accounting::Trace::new(
            m,
            policy.strict_bands,
            self.tolerance,
            policy.buffer_tolerance,
        );
        work.apply(basis.cols(), None);
        trace.blocks = 1;
        let mut trace_values = vec![0.0; m];
        let mut locked = buffer(n, m)?;
        let mut hl = buffer(n, m)?;
        let mut x = buffer(n, m)?;
        let mut hx = buffer(n, m)?;
        let mut all = buffer(n, capacity)?;
        let mut hall = buffer(n, capacity)?;
        let mut corrections = buffer(n, m)?;
        let mut preconditioned = buffer(n, m)?;
        let mut h_extra = buffer(n, m)?;
        let mut extra = buffer(n, m)?;
        let mut dropped = Vec::with_capacity(m);
        let mut shifts = Vec::with_capacity(m);
        let mut norms = vec![0.0; m];

        let mut history = Vec::with_capacity(self.max_iterations);
        let mut report = LockingReport {
            locked_roots: 0,
            unlocks: 0,
            failures: Vec::with_capacity(m * self.max_iterations),
        };
        let mut reduced = ComplexMatrix::zeros(capacity, capacity)?;
        let mut cache = ComplexMatrix::zeros(capacity, capacity)?;
        let mut cache_valid = false;
        let mut eigensolver = dft_numerics::linalg::HermitianWorkspace::new(capacity)?;
        let mut coefficients = ComplexMatrix::zeros(capacity, capacity)?;
        let mut values = Vec::with_capacity(capacity);
        let mut final_values = Vec::with_capacity(m);
        let mut final_reduced = ComplexMatrix::zeros(m, m)?;
        let mut rotation = ComplexMatrix::zeros(m, m)?;
        let mut vectors = buffer(n, m)?;
        vectors.resize_columns(m)?;
        let mut residual_block = buffer(n, m)?;
        residual_block.resize_columns(m)?;
        let mut residuals = vec![0.0; m];
        // A failed all-root audit disables subsequent locking in this solve.
        // Retrying the full space avoids repeatedly locking the same bad state.
        let mut allow_lock = true;
        for iteration in 1..=self.max_iterations {
            gemm.accounting_step(iteration);
            let count = locked.cols();
            let active = m - count;
            gemm.reduced_hamiltonian_into(
                &basis,
                &hb,
                cache_valid.then_some(&cache),
                option_env!("P24A_TWO_RECTANGLES").is_none(),
                &mut reduced,
            )?;
            eigensolver.solve_into(&reduced, &mut values, &mut coefficients)?;
            std::mem::swap(&mut reduced, &mut cache);
            cache_valid = true;
            trace_values[count..].copy_from_slice(&values[..active]);
            coefficients.resize_columns(active)?;
            x.resize_columns(active)?;
            hx.resize_columns(active)?;
            gemm.accounting_site("ritz_and_residual");
            gemm.multiply_into(&basis, Identity, &coefficients, Identity, &mut x)?;
            gemm.multiply_into(&hb, Identity, &coefficients, Identity, &mut hx)?;
            corrections.resize_columns(active)?;
            for j in 0..active {
                for i in 0..n {
                    corrections[(i, j)] = hx[(i, j)] - values[j] * x[(i, j)];
                }
                norms[count + j] = norm(corrections.column(j));
            }
            if norms.iter().any(|x| !x.is_finite()) {
                return Err(LinalgError::NonFinite.into());
            }
            let finished = (0..m).all(|j| norms[j] <= threshold(j));
            if finished || iteration == self.max_iterations {
                all.resize_columns(count)?;
                all.append_columns(&x)?;
                hall.replace_columns(&hl)?;
                hall.append_columns(&hx)?;
                // Independent product: no conjugate fill before Hermiticity validation.
                gemm.accounting_site("locking_final_gram");
                gemm.multiply_into(&all, Adjoint, &hall, Identity, &mut final_reduced)?;
                eigensolver.solve_into(&final_reduced, &mut final_values, &mut rotation)?;
                gemm.accounting_site("locking_final_rotation");
                gemm.multiply_into(&all, Identity, &rotation, Identity, &mut vectors)?;
                operator.apply_block(&vectors, &mut residual_block)?;
                for j in 0..m {
                    for i in 0..n {
                        residual_block[(i, j)] -= vectors[(i, j)] * final_values[j];
                    }
                    residuals[j] = norm(residual_block.column(j));
                }
                if residuals.iter().any(|r| !r.is_finite()) {
                    return Err(LinalgError::NonFinite.into());
                }
                applications += m;
                work.apply(m, Some(m - policy.strict_bands));
                trace.blocks += 1;
                if replay::ENABLED {
                    trace.step(
                        basis.cols() + count,
                        applications,
                        &final_values,
                        &residuals,
                    );
                }
                work.residuals(iteration, &residuals);
                let maximum = max_residual(&residuals);
                history.push(maximum);
                if (0..m).all(|j| residuals[j] <= threshold(j)) {
                    work.locking(&report);
                    work.returned();
                    return Ok(EigenSolution {
                        locking: Some(report),
                        spectral_completeness: None,
                        eigenvalues: final_values,
                        eigenvectors: vectors,
                        residuals,
                        converged: vec![true; m],
                        buffer: (0..m).map(|j| j >= policy.strict_bands).collect(),
                        residuals_certified: true,
                        iterations: iteration,
                        history,
                        operator_applications: applications,
                    });
                }
                for (root, &r) in residuals.iter().enumerate() {
                    if r > threshold(root) {
                        report.failures.push((iteration, root, r));
                    }
                }
                if iteration == self.max_iterations {
                    return Err(SolverError::NotConverged {
                        iterations: iteration,
                        max_residual: maximum,
                        residuals,
                    });
                }
                report.unlocks += count;
                allow_lock = false;
                locked.resize_columns(0)?;
                hl.resize_columns(0)?;

                basis.replace_columns(&vectors)?;
                hb.resize_columns(m)?;
                operator.apply_block(&basis, &mut hb)?;
                applications += m;
                work.apply(m, None);
                work.restart();
                trace.blocks += 1;
                trace.restarts += 1;
                cache_valid = false;
                continue;
            }
            if replay::ENABLED {
                trace.step(basis.cols() + count, applications, &trace_values, &norms);
            }
            work.residuals(iteration, &norms);
            history.push(max_residual(&norms));
            // Cluster gap in eigenvalue units: ten requested strict tolerances,
            // with an epsilon-scaled floor. Never split a boundary cluster.
            let newly_locked = if allow_lock {
                lockable_prefix(
                    &values,
                    &norms[count..],
                    count,
                    self.tolerance,
                    policy,
                    margin,
                )
            } else {
                0
            };
            if newly_locked > 0 {
                let old = locked.cols();
                locked.resize_columns(old + newly_locked)?;
                hl.resize_columns(old + newly_locked)?;
                all.resize_columns(old + newly_locked)?;
                for j in 0..newly_locked {
                    locked.column_mut(old + j).copy_from_slice(x.column(j));
                    // Preserve this projection prefix across subsequent steps.
                    all.column_mut(old + j).copy_from_slice(x.column(j));
                    hl.column_mut(old + j).copy_from_slice(hx.column(j));
                }
                report.locked_roots += newly_locked;
            }
            // A lock transition contracts to retained Ritz vectors. Subsequent
            // expansions retain their complete active space until capacity.
            let restart = newly_locked > 0
                || basis.cols() + locked.cols() + self.block_size.min(active) > capacity;
            if restart {
                work.restart();
                trace.restarts += 1;
                basis.resize_columns(active - newly_locked)?;
                hb.resize_columns(active - newly_locked)?;
                for j in newly_locked..active {
                    basis
                        .column_mut(j - newly_locked)
                        .copy_from_slice(x.column(j));
                    hb.column_mut(j - newly_locked)
                        .copy_from_slice(hx.column(j));
                }
                cache_valid = false;
            }
            shifts.clear();
            let mut selected = 0;
            for j in newly_locked..active {
                if norms[count + j] > threshold(count + j)
                    && selected < self.block_size.min(capacity - m)
                {
                    for i in 0..n {
                        corrections[(i, selected)] = hx[(i, j)] - values[j] * x[(i, j)];
                    }
                    shifts.push(values[j]);
                    selected += 1;
                }
            }
            corrections.resize_columns(selected)?;
            preconditioned.resize_columns(selected)?;
            if let Some(p) = preconditioner {
                p.apply(&shifts, &corrections, &mut preconditioned)?;
            } else {
                preconditioned.replace_columns(&corrections)?;
            }
            all.resize_columns(locked.cols())?;
            all.append_columns(&basis)?;
            projection.orthonormalize_into(
                &preconditioned,
                Some(&all),
                &mut gemm,
                &mut extra,
                &mut dropped,
            )?;
            if extra.cols() == 0 {
                let width = self.block_size.min(m).min(capacity - all.cols());
                preconditioned.resize_columns(width)?;
                for j in 0..width {
                    for i in 0..n {
                        preconditioned[(i, j)] = Complex64::new(rng.next(), rng.next());
                    }
                }
                projection.orthonormalize_into(
                    &preconditioned,
                    Some(&all),
                    &mut gemm,
                    &mut extra,
                    &mut dropped,
                )?;
            }
            if extra.cols() == 0 {
                return Err(SolverError::NotConverged {
                    iterations: iteration,
                    max_residual: max_residual(&norms),
                    residuals: norms,
                });
            }
            h_extra.resize_columns(extra.cols())?;
            operator.apply_block(&extra, &mut h_extra)?;
            applications += extra.cols();
            work.apply(extra.cols(), None);
            trace.blocks += 1;
            basis.append_columns(&extra)?;
            hb.append_columns(&h_extra)?;
        }
        unreachable!("nonempty iteration budget validated")
    }
}

// The full reduced spectrum is needed to avoid splitting a cluster at the
// requested-band boundary. Only a contiguous lowest-root prefix may lock.
fn lockable_prefix(
    values: &[f64],
    residuals: &[f64],
    offset: usize,
    tolerance: f64,
    policy: BandTolerances,
    margin: f64,
) -> usize {
    let mut start = 0;
    while start < residuals.len() {
        let mut end = start + 1;
        while end < values.len()
            && values[end] - values[end - 1]
                <= (10.0 * tolerance).max(64.0 * f64::EPSILON * values[end].abs().max(1.0))
        {
            end += 1;
        }
        if end > residuals.len()
            || (start..end).any(|j| {
                residuals[j]
                    > margin
                        * if offset + j < policy.strict_bands {
                            tolerance
                        } else {
                            policy.buffer_tolerance
                        }
            })
        {
            break;
        }
        start = end;
    }
    start
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn clusters_are_atomic_and_never_split_at_requested_boundary() {
        let p = BandTolerances {
            strict_bands: 3,
            buffer_tolerance: 1e-8,
        };
        let v = [-2.0, 0.0, 1e-9, 2.0];
        assert_eq!(lockable_prefix(&v, &[0.0, 0.0, 1e-6], 0, 1e-8, p, 0.5), 1);
        assert_eq!(lockable_prefix(&v, &[0.0, 0.0, 0.0], 0, 1e-8, p, 0.5), 3);
        assert_eq!(lockable_prefix(&v, &[0.0, 0.0], 0, 1e-8, p, 0.5), 1);
        assert_eq!(lockable_prefix(&v, &[1e-6, 0.0, 0.0], 0, 1e-8, p, 0.5), 0);
    }
}
