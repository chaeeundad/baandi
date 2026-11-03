// Frozen pre-O4 implementation, retained only as a regression oracle.
use super::*;
impl Davidson {
    pub(super) fn solve_legacy(
        &self,
        operator: &dyn LinearOperator,
        initial: Option<&ComplexMatrix>,
        preconditioner: Option<&dyn Preconditioner>,
    ) -> Result<EigenSolution, SolverError> {
        let n = operator.dimension();
        validate_bands(n, self.bands)?;
        let capacity = self.max_subspace.min(n);
        if !self.tolerance.is_finite()
            || self.tolerance <= 0.0
            || self.max_iterations == 0
            || self.block_size == 0
            || capacity < self.bands
            || (capacity == self.bands && self.bands < n)
        {
            return Err(SolverError::InvalidParameter(
                "positive tolerance/iterations/block size and subspace larger than bands required (unless bands=n)",
            ));
        }
        let mut rng = SplitMix64(self.seed);
        let mut basis = if let Some(initial) = initial {
            if initial.rows() != n || initial.cols() > capacity {
                return Err(LinalgError::Dimension.into());
            }
            orthonormalize(initial, None)?.vectors
        } else {
            ComplexMatrix::zeros(n, 0)?
        };
        let initial_width = self.bands.max(self.block_size).min(capacity);
        // Fill missing rank with independent seeded columns, including deficient user guesses.
        while basis.cols() < initial_width {
            let random = rng.block(n, initial_width - basis.cols())?;
            let extra = orthonormalize(&random, Some(&basis))?.vectors;
            if extra.cols() == 0 {
                return Err(SolverError::InvalidParameter(
                    "unable to complete initial basis",
                ));
            }
            basis = concatenate(&basis, &extra)?;
        }
        let mut history = Vec::new();
        let mut operator_applications = 0;
        for iteration in 1..=self.max_iterations {
            operator_applications += basis.cols() + self.bands;
            let h_basis = apply(operator, &basis)?;
            let projected = gemm(&basis, Adjoint, &h_basis, Identity)?;
            let (mut values, coefficients) = hermitian_eigh(&projected)?;
            values.truncate(self.bands);
            let vectors = gemm(
                &basis,
                Identity,
                &first_columns(&coefficients, self.bands)?,
                Identity,
            )?;
            // Reapply H: convergence is based on actual operator residuals, not a recurrence.
            let (residual_block, residuals) = residual(operator, &vectors, &values)?;
            let maximum = max_residual(&residuals);
            history.push(maximum);
            if residuals.iter().all(|&r| r <= self.tolerance) {
                return Ok(EigenSolution {
                    eigenvalues: values,
                    eigenvectors: vectors,
                    residuals,
                    iterations: iteration,
                    history,
                    operator_applications,
                });
            }
            if iteration == self.max_iterations {
                return Err(SolverError::NotConverged {
                    iterations: iteration,
                    max_residual: maximum,
                    residuals,
                });
            }
            let selected: Vec<_> = residuals
                .iter()
                .enumerate()
                .filter(|(_, r)| **r > self.tolerance)
                .map(|(j, _)| j)
                .take(self.block_size.min(capacity - self.bands))
                .collect();
            let mut corrections = ComplexMatrix::zeros(n, selected.len())?;
            let mut shifts = Vec::new();
            for (k, &j) in selected.iter().enumerate() {
                corrections
                    .column_mut(k)
                    .copy_from_slice(residual_block.column(j));
                shifts.push(values[j]);
            }
            if let Some(p) = preconditioner {
                let mut output = ComplexMatrix::zeros(n, corrections.cols())?;
                p.apply(&shifts, &corrections, &mut output)?;
                if output.rows() != n || output.cols() != corrections.cols() {
                    return Err(LinalgError::Dimension.into());
                }
                output.validate_finite()?;
                corrections = output;
            }
            if basis.cols() + corrections.cols() > capacity {
                basis = orthonormalize(&vectors, None)?.vectors;
            }
            let mut extra = orthonormalize(&corrections, Some(&basis))?.vectors;
            // A stalled direction can occur with an exact preconditioner or an invariant
            // initial space. Try fresh independent directions, never infer convergence.
            if extra.cols() == 0 {
                let width = self.block_size.min(capacity - basis.cols());
                let random = rng.block(n, width)?;
                extra = orthonormalize(&random, Some(&basis))?.vectors;
            }
            if extra.cols() == 0 {
                return Err(SolverError::NotConverged {
                    iterations: iteration,
                    max_residual: maximum,
                    residuals,
                });
            }
            basis = concatenate(&basis, &extra)?;
        }
        unreachable!("positive iteration limit validated above")
    }
}
