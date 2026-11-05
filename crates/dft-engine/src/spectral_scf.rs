//! Rebuild the reported ordinary map's Hamiltonian after its fresh fixed-point
//! check. The returned candidate is never compared against the fresh map's H.
use super::*;
impl ScfProblem {
    pub(super) fn audit_reported_spectrum(
        &mut self,
        input: &Array3<f64>,
        points: &mut [KPointSolution],
        iteration: usize,
        report: &mut crate::SpectralReport,
    ) -> Result<(Option<Vec<ComplexMatrix>>, bool), ScfError> {
        let start = Instant::now();
        let h = hartree::evaluate(input, &self.cell, &mut self.fft)?;
        let x = self.options.functional.evaluate_with_floor(
            input,
            self.core_density.as_ref(),
            &mut self.xc_density,
            &self.cell,
            &mut self.fft,
            self.density_floor,
        )?;
        let mut effective = self.ionic_real.clone();
        for ((v, h), x) in effective.iter_mut().zip(&h.potential).zip(&x.potential) {
            *v += h.re + x;
        }
        let bands = self.occupations.len();
        let solver = Davidson {
            bands,
            tolerance: 1e-10,
            max_iterations: self.options.davidson_max_iterations,
            max_subspace: self.options.davidson_subspace.max(6 * bands),
            block_size: bands,
            seed: self.options.seed,
        };
        let first = report.events.len();
        let mut vectors = Vec::with_capacity(points.len());
        let mut failed = false;
        for (k, (state, point)) in self.states.iter_mut().zip(points).enumerate() {
            state.hamiltonian.set_local_potential(effective.clone())?;
            let pre =
                DiagonalPreconditioner::new(state.hamiltonian.kinetic_diagonal().to_vec(), 0.1)?;
            if let Some(checked) = report.audit(
                iteration,
                None,
                k,
                &solver,
                &state.hamiltonian,
                &pre,
                &point.solution,
            ) {
                if checked.spectral_completeness.as_ref().unwrap().retries == 0 {
                    point.solution.spectral_completeness = checked.spectral_completeness.clone();
                }
                vectors.push(checked.eigenvectors);
            } else {
                failed = true;
            }
        }
        let repaired = report.events[first..].iter().any(|e| e.retries > 0);
        report.passed = report.events[first..].iter().all(|e| e.passed);
        report.elapsed += start.elapsed();
        Ok(((!failed && repaired).then_some(vectors), failed))
    }
}
