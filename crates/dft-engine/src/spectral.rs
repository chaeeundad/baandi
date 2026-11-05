//! Final-state audit records. I/O and serialization belong to the CLI.
use dft_solvers::{Davidson, DiagonalPreconditioner, EigenSolution, LinearOperator};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SpectralGuard {
    #[default]
    Off,
    Final,
}
impl SpectralGuard {
    pub fn name(self) -> &'static str {
        match self {
            Self::Off => "off",
            Self::Final => "final",
        }
    }
}
#[derive(Clone, Debug)]
pub struct SpectralEvent {
    pub iteration: usize,
    pub spin: Option<usize>,
    pub k_point: usize,
    pub method: &'static str,
    /// True only if the original reported block agreed; a repaired block must
    /// enter a new SCF iteration and fixed-point certificate.
    pub passed: bool,
    pub retries: usize,
    pub operator_applications: usize,
    pub discrepancy_ha: Option<f64>,
    pub error: Option<String>,
}
#[derive(Clone, Debug, Default)]
pub struct SpectralReport {
    pub elapsed: std::time::Duration,
    pub passed: bool,
    pub retries: usize,
    pub operator_applications: usize,
    pub events: Vec<SpectralEvent>,
}
impl SpectralReport {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn audit(
        &mut self,
        iteration: usize,
        spin: Option<usize>,
        k_point: usize,
        solver: &Davidson,
        operator: &dyn LinearOperator,
        preconditioner: &DiagonalPreconditioner,
        original: &EigenSolution,
    ) -> Option<EigenSolution> {
        let (result, work, retries) =
            solver.guard_spectrum_measured(operator, Some(preconditioner), original.clone());
        self.operator_applications += work;
        self.retries += retries;
        let (solution, method, passed, discrepancy, error) = match result {
            Ok(solution) => {
                let meta = solution
                    .spectral_completeness
                    .as_ref()
                    .expect("guard metadata");
                let fields = (meta.method, retries == 0, Some(meta.maximum_discrepancy));
                (Some(solution), fields.0, fields.1, fields.2, None)
            }
            Err(error) => (
                None,
                "independent_spectrum",
                false,
                None,
                Some(error.to_string()),
            ),
        };
        self.events.push(SpectralEvent {
            iteration,
            spin,
            k_point,
            method,
            passed,
            retries,
            operator_applications: work,
            discrepancy_ha: discrepancy,
            error,
        });
        solution
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use dft_physics::data::{Complex64, ComplexMatrix};
    use dft_solvers::SolverError;
    use std::cell::Cell;
    struct Diagonal {
        calls: Cell<usize>,
        fail: bool,
    }
    impl LinearOperator for Diagonal {
        fn dimension(&self) -> usize {
            8
        }
        fn apply(&self, x: &[Complex64], y: &mut [Complex64]) -> Result<(), SolverError> {
            self.calls.set(self.calls.get() + 1);
            if self.fail {
                return Err(SolverError::Application("injected probe failure".into()));
            }
            for ((y, x), e) in y
                .iter_mut()
                .zip(x)
                .zip([-3., -3., -1., 2., 4., 6., 8., 10.])
            {
                *y = *x * e;
            }
            Ok(())
        }
    }
    #[test]
    fn cumulative_retries_and_failed_probes_are_not_lost() {
        let op = Diagonal {
            calls: Cell::new(0),
            fail: false,
        };
        let cfg = Davidson {
            bands: 2,
            tolerance: 1e-10,
            max_iterations: 100,
            max_subspace: 8,
            block_size: 2,
            seed: 2026,
        };
        let mut guess = ComplexMatrix::zeros(8, 2).unwrap();
        guess[(1, 0)] = 1.0.into();
        guess[(2, 1)] = 1.0.into();
        let original = cfg.solve(&op, Some(&guess), None).unwrap();
        // A scalar preconditioner tests retry accounting without the separate
        // exact-diagonal cancellation/rank-loss mechanism.
        let pre = DiagonalPreconditioner::new(vec![0.0; 8], 0.1).unwrap();
        op.calls.set(0);
        let mut report = SpectralReport::default();
        for iteration in [6, 10] {
            let repaired = report.audit(iteration, None, 0, &cfg, &op, &pre, &original);
            assert!(repaired.is_some(), "{:?}", report.events);
            let repaired = repaired.unwrap();
            assert!((repaired.eigenvalues[1] + 3.).abs() < 1e-9);
        }
        assert_eq!(report.retries, 2 * report.events[0].retries);
        assert_eq!(report.operator_applications, op.calls.get());
        assert!(report.events.iter().all(|e| !e.passed && e.retries > 0));
        let fail = Diagonal {
            calls: Cell::new(0),
            fail: true,
        };
        assert!(
            report
                .audit(11, None, 0, &cfg, &fail, &pre, &original)
                .is_none()
        );
        assert_eq!(report.retries, 2 * report.events[0].retries);
        assert!(
            report.events[2]
                .error
                .as_ref()
                .unwrap()
                .contains("injected")
        );
        assert!(!report.passed);
        assert_eq!(
            report.operator_applications,
            report.events.iter().map(|e| e.operator_applications).sum()
        );
    }
}
