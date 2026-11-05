//! Independent-start spectral guard. Residual certification alone cannot count
//! roots. Large operators use a deterministic independent search, not a rigorous
//! inertia certificate; small operators are checked in the complete space.
use super::*;

#[derive(Clone, Debug, PartialEq)]
pub struct SpectralCompleteness {
    pub method: &'static str,
    pub passed: bool,
    pub retries: usize,
    pub extra_operator_applications: usize,
    pub maximum_discrepancy: f64,
}

impl Davidson {
    /// Preserve the supplied solution bit for bit when its spectrum agrees.
    /// A discrepancy triggers an independent wider retry before replacement.
    /// Every failed probe/retry returns an error; it never sets `passed`.
    pub fn guard_spectrum(
        &self,
        operator: &dyn LinearOperator,
        preconditioner: Option<&dyn Preconditioner>,
        original: EigenSolution,
    ) -> Result<EigenSolution, SolverError> {
        self.guard_spectrum_measured(operator, preconditioner, original)
            .0
    }

    /// Count all production H columns, including a probe that returns an error.
    pub fn guard_spectrum_measured(
        &self,
        operator: &dyn LinearOperator,
        preconditioner: Option<&dyn Preconditioner>,
        original: EigenSolution,
    ) -> (Result<EigenSolution, SolverError>, usize, usize) {
        let counted = Counted {
            operator,
            columns: std::cell::Cell::new(0),
        };
        let retries = std::cell::Cell::new(0);
        let native_work = original.operator_applications;
        let mut result = self.guard_inner(&counted, preconditioner, original, &retries);
        if let Ok(solution) = &mut result {
            solution.operator_applications = native_work + counted.columns.get();
            let meta = solution
                .spectral_completeness
                .as_mut()
                .expect("guard metadata");
            meta.extra_operator_applications = counted.columns.get();
            meta.retries = retries.get();
        }
        (result, counted.columns.get(), retries.get())
    }

    fn guard_inner(
        &self,
        operator: &dyn LinearOperator,
        preconditioner: Option<&dyn Preconditioner>,
        mut original: EigenSolution,
        retry_count: &std::cell::Cell<usize>,
    ) -> Result<EigenSolution, SolverError> {
        const CAP: f64 = 1e-9;
        let n = operator.dimension();
        validate_bands(n, self.bands)?;
        if original.eigenvalues.len() != self.bands
            || original.residuals.len() != self.bands
            || original.eigenvectors.rows() != n
            || original.eigenvectors.cols() != self.bands
            || original.eigenvalues.iter().any(|v| !v.is_finite())
            || original.eigenvalues.windows(2).any(|v| v[0] > v[1])
        {
            return Err(SolverError::InvalidParameter(
                "spectral guard requires a finite, sorted, dimensionally consistent solution",
            ));
        }
        original.eigenvectors.validate_finite()?;
        let dense = n <= 128;
        let mut screening_work = 0;
        if !dense && self.bands < n && option_env!("P15C_FULL_GUARD").is_none() {
            let complement = Complement::new(operator, &original);
            let screening = Davidson {
                bands: 1,
                block_size: 4,
                tolerance: 1e-10,
                seed: self.seed ^ 0x94d049bb133111eb,
                ..self.clone()
            };
            let screen = screening.solve_policy_traced::<false>(
                &complement,
                None,
                preconditioner,
                SolveMode {
                    ritz: false,
                    recurrence: true,
                    contract_restart: true,
                    certify: true,
                },
                BandTolerances {
                    strict_bands: 1,
                    buffer_tolerance: 1e-10,
                },
                &mut replay::Trace::default(),
            )?;
            screening_work = screen.operator_applications;
            if screen.eigenvalues[0] >= original.eigenvalues[self.bands - 1] - CAP {
                original.operator_applications += screening_work;
                original.spectral_completeness = Some(SpectralCompleteness {
                    method: "independent_complement_davidson",
                    passed: true,
                    retries: 0,
                    extra_operator_applications: screening_work,
                    maximum_discrepancy: 0.0,
                });
                return Ok(original);
            }
        }
        let mut probe = if dense {
            DenseEigensolver.solve(operator, self.bands)?
        } else {
            self.spectral_probe(operator, preconditioner, 2, self.seed ^ 0x9e3779b97f4a7c15)?
        };
        if probe.residuals.iter().any(|r| !r.is_finite() || *r > 1e-10) {
            return Err(SolverError::Application(
                "spectral completeness: probe residual failed".into(),
            ));
        }
        let discrepancy = difference(&original.eigenvalues, &probe.eigenvalues);
        let mut extra = screening_work + probe.operator_applications;
        let mut retries = 0;
        if discrepancy > CAP {
            // A second independent start is needed before replacing a large
            // operator's state. The full-space small-operator oracle needs none.
            retry_count.set(1);
            if dense {
                // The dense spectrum detects the failure; recover through the
                // same seeded enlarged-subspace policy used by iterative solves.
                let first = self.spectral_probe(
                    operator,
                    preconditioner,
                    2,
                    self.seed ^ 0xd1b54a32d192ed03,
                );
                let repaired = match first {
                    Ok(repaired)
                        if difference(&repaired.eigenvalues, &probe.eigenvalues) <= CAP =>
                    {
                        repaired
                    }
                    _ => {
                        retry_count.set(2);
                        self.spectral_probe(
                            operator,
                            preconditioner,
                            3,
                            self.seed ^ 0x94d049bb133111eb,
                        )?
                    }
                };
                if difference(&repaired.eigenvalues, &probe.eigenvalues) > CAP {
                    return Err(SolverError::Application(
                        "spectral completeness: enlarged solve disagrees with dense spectrum"
                            .into(),
                    ));
                }
                probe = repaired;
            } else {
                let retry = self.spectral_probe(
                    operator,
                    preconditioner,
                    3,
                    self.seed ^ 0xd1b54a32d192ed03,
                )?;
                extra += retry.operator_applications;
                if difference(&retry.eigenvalues, &probe.eigenvalues) > CAP {
                    return Err(SolverError::Application(
                        "spectral completeness: independent searches disagree".into(),
                    ));
                }
            }
            if probe
                .eigenvalues
                .iter()
                .zip(&original.eigenvalues)
                .any(|(p, o)| p - o > CAP)
            {
                return Err(SolverError::Application(
                    "spectral completeness: independent search missed an existing lower root"
                        .into(),
                ));
            }
            retries = retry_count.get();
            let old_work = original.operator_applications;
            let buffer = original.buffer.clone();
            original = probe;
            // The buffer classification is a policy, not a property of the probe.
            original.buffer = buffer;
            original.operator_applications = old_work;
        }
        original.operator_applications += extra;
        original.spectral_completeness = Some(SpectralCompleteness {
            method: if dense {
                "full_space_dense"
            } else {
                "independent_seeded_davidson"
            },
            passed: true,
            retries,
            extra_operator_applications: extra,
            maximum_discrepancy: discrepancy,
        });
        Ok(original)
    }

    fn spectral_probe(
        &self,
        operator: &dyn LinearOperator,
        preconditioner: Option<&dyn Preconditioner>,
        multiplier: usize,
        seed: u64,
    ) -> Result<EigenSolution, SolverError> {
        let solver = Self {
            tolerance: 1e-10,
            seed,
            max_subspace: self.max_subspace.max(multiplier * self.bands + 1),
            ..self.clone()
        };
        let width = (multiplier * self.bands).min(operator.dimension());
        let start = SplitMix64(seed).block(operator.dimension(), width)?;
        solver.solve_policy_traced::<false>(
            operator,
            Some(&start),
            preconditioner,
            SolveMode {
                ritz: false,
                recurrence: true,
                contract_restart: true,
                certify: true,
            },
            BandTolerances {
                strict_bands: self.bands,
                buffer_tolerance: 1e-10,
            },
            &mut replay::Trace::default(),
        )
    }
}

fn difference(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max)
}

/// P H P + shift Q Q*, P=I-Q Q*. Scratch is prepared outside operator calls.
/// One independent lowest root in the complement screens for missed states.
struct Complement<'a> {
    operator: &'a dyn LinearOperator,
    vectors: &'a ComplexMatrix,
    shift: f64,
    scratch: std::cell::RefCell<(Vec<Complex64>, Vec<Complex64>)>,
}
impl<'a> Complement<'a> {
    fn new(operator: &'a dyn LinearOperator, solution: &'a EigenSolution) -> Self {
        let top = solution.eigenvalues[solution.eigenvalues.len() - 1];
        let n = operator.dimension();
        Self {
            operator,
            vectors: &solution.eigenvectors,
            shift: top + top.abs().max(1.0),
            scratch: std::cell::RefCell::new((
                vec![Complex64::default(); n],
                vec![Complex64::default(); n],
            )),
        }
    }
    fn project(&self, x: &mut [Complex64]) {
        for j in 0..self.vectors.cols() {
            let q = self.vectors.column(j);
            let overlap: Complex64 = q.iter().zip(x.iter()).map(|(q, x)| q.conj() * x).sum();
            for (x, q) in x.iter_mut().zip(q) {
                *x -= overlap * q;
            }
        }
    }
}
impl LinearOperator for Complement<'_> {
    fn dimension(&self) -> usize {
        self.operator.dimension()
    }
    fn apply(&self, x: &[Complex64], y: &mut [Complex64]) -> Result<(), SolverError> {
        let mut scratch = self.scratch.borrow_mut();
        let (projected, product) = &mut *scratch;
        projected.copy_from_slice(x);
        self.project(projected);
        self.operator.apply(projected, product)?;
        self.project(product);
        for ((y, h), (x, p)) in y.iter_mut().zip(product).zip(x.iter().zip(projected)) {
            *y = *h + self.shift * (*x - *p);
        }
        Ok(())
    }
}

struct Counted<'a> {
    operator: &'a dyn LinearOperator,
    columns: std::cell::Cell<usize>,
}
impl LinearOperator for Counted<'_> {
    fn dimension(&self) -> usize {
        self.operator.dimension()
    }
    fn worker_pool(&self) -> Option<std::sync::Arc<dft_numerics::pool::WorkerPool>> {
        self.operator.worker_pool()
    }
    fn inner_threads(&self) -> usize {
        self.operator.inner_threads()
    }
    fn apply(&self, x: &[Complex64], y: &mut [Complex64]) -> Result<(), SolverError> {
        self.columns.set(self.columns.get() + 1);
        self.operator.apply(x, y)
    }
    fn apply_block(&self, x: &ComplexMatrix, y: &mut ComplexMatrix) -> Result<(), SolverError> {
        self.columns.set(self.columns.get() + x.cols());
        self.operator.apply_block(x, y)
    }
}
