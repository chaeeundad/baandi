//! Charge-conserving mixing of real-space electron number densities.
use crate::SolverError;
use dft_numerics::linalg::{ComplexMatrix, hermitian_eigh};

/// Input and output are samples on the same uniform grid, in electrons/Bohr^3.
/// R = output - input. `mean` is the prescribed N_e / volume, not a fitted charge.
pub trait DensityMixer {
    fn mix(&mut self, input: &[f64], output: &[f64], mean: f64) -> Result<Vec<f64>, SolverError>;
    fn reset(&mut self);
}
fn validate(input: &[f64], output: &[f64], mean: f64, beta: f64) -> Result<(), SolverError> {
    if input.is_empty()
        || input.len() != output.len()
        || !mean.is_finite()
        || mean <= 0.0
        || !beta.is_finite()
        || beta <= 0.0
        || input.iter().chain(output).any(|v| !v.is_finite())
    {
        return Err(SolverError::InvalidParameter(
            "finite equal nonempty densities, positive mean and beta required",
        ));
    }
    Ok(())
}
fn charge(mut n: Vec<f64>, mean: f64) -> Result<Vec<f64>, SolverError> {
    // Two uniform zero-mode corrections remove summation roundoff without
    // selecting a distinguished spatial point or changing nonzero modes.
    for _ in 0..2 {
        let shift = mean - n.iter().sum::<f64>() / n.len() as f64;
        for v in &mut n {
            *v += shift;
        }
    }
    if n.iter().any(|v| !v.is_finite()) {
        return Err(dft_numerics::linalg::LinalgError::NonFinite.into());
    }
    Ok(n)
}
#[derive(Debug, Clone)]
pub struct LinearMixer {
    pub beta: f64,
}
impl DensityMixer for LinearMixer {
    fn mix(&mut self, input: &[f64], output: &[f64], mean: f64) -> Result<Vec<f64>, SolverError> {
        validate(input, output, mean, self.beta)?;
        charge(
            input
                .iter()
                .zip(output)
                .map(|(a, b)| a + self.beta * (b - a))
                .collect(),
            mean,
        )
    }
    fn reset(&mut self) {}
}
/// Type-II Anderson: x+beta R - (Delta x+beta Delta R) gamma,
/// minimizing ||R-Delta R gamma||. The Gram eigensolve is in dft-numerics.
/// Eigenvalues <= 1e-12 * largest are discarded (spectral regularization).
#[derive(Debug)]
pub struct PulayMixer {
    pub beta: f64,
    capacity: usize,
    metric: Option<Box<dyn ResidualMetric>>,
    metric_vectors: Vec<Vec<f64>>,
    metric_r: Vec<f64>,
    pairs: Vec<(Vec<f64>, Vec<f64>)>,
    mean: Option<f64>,
}
/// A positive semidefinite density metric. Applies its symmetric kernel to
/// a residual; unlike preconditioning, it changes only the least-squares fit.
/// Implementations own reusable transform scratch and have no history.
pub trait ResidualMetric: std::fmt::Debug {
    fn duplicate(&self) -> Box<dyn ResidualMetric>;
    fn apply(&mut self, residual: &[f64], output: &mut [f64]) -> Result<(), SolverError>;
}
impl Clone for PulayMixer {
    fn clone(&self) -> Self {
        Self {
            beta: self.beta,
            capacity: self.capacity,
            pairs: self.pairs.clone(),
            mean: self.mean,
            metric: self.metric.as_ref().map(|m| m.duplicate()),
            metric_vectors: self.metric_vectors.clone(),
            metric_r: self.metric_r.clone(),
        }
    }
}
impl PulayMixer {
    pub fn new(beta: f64, history: usize) -> Result<Self, SolverError> {
        if !beta.is_finite() || beta <= 0.0 || history == 0 || history.checked_add(1).is_none() {
            return Err(SolverError::InvalidParameter(
                "positive beta and history required",
            ));
        }
        Ok(Self {
            beta,
            capacity: history,
            metric: None,
            metric_vectors: Vec::new(),
            metric_r: Vec::new(),
            pairs: Vec::new(),
            mean: None,
        })
    }
    pub fn set_metric(&mut self, metric: Box<dyn ResidualMetric>) {
        self.reset();
        self.metric = Some(metric);
    }
    pub fn history_len(&self) -> usize {
        self.pairs.len()
    }
}
impl Default for PulayMixer {
    fn default() -> Self {
        Self::new(0.3, 6).expect("valid defaults")
    }
}
impl DensityMixer for PulayMixer {
    fn reset(&mut self) {
        self.pairs.clear();
        self.mean = None;
    }
    fn mix(&mut self, input: &[f64], output: &[f64], mean: f64) -> Result<Vec<f64>, SolverError> {
        validate(input, output, mean, self.beta)?;
        if self.mean != Some(mean) || self.pairs.first().is_some_and(|p| p.0.len() != input.len()) {
            self.reset();
        }
        self.mean = Some(mean);
        if self.metric.is_some() && self.metric_r.len() != input.len() {
            // Session initialization validates a zero step before SCF, so this
            // metric scratch is allocated once, not during the SCF hot loop.
            self.metric_vectors = vec![vec![0.0; input.len()]; self.capacity];
            self.metric_r = vec![0.0; input.len()];
        }
        let r: Vec<_> = output.iter().zip(input).map(|(b, a)| b - a).collect();
        self.pairs.push((input.to_vec(), r.clone()));
        if self.pairs.len() > self.capacity + 1 {
            self.pairs.remove(0);
        }
        let mut next: Vec<_> = input
            .iter()
            .zip(&r)
            .map(|(a, r)| a + self.beta * r)
            .collect();
        let m = self.pairs.len() - 1;
        if m > 0 {
            let dx: Vec<Vec<f64>> = self
                .pairs
                .windows(2)
                .map(|p| p[1].0.iter().zip(&p[0].0).map(|(a, b)| a - b).collect())
                .collect();
            let dr: Vec<Vec<f64>> = self
                .pairs
                .windows(2)
                .map(|p| p[1].1.iter().zip(&p[0].1).map(|(a, b)| a - b).collect())
                .collect();
            let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
            let weighted = &mut self.metric_vectors;
            let weighted_r = &mut self.metric_r;
            if let Some(metric) = &mut self.metric {
                for (source, dest) in dr.iter().zip(weighted.iter_mut()) {
                    metric.apply(source, dest)?;
                }
                metric.apply(&r, weighted_r)?;
            }
            let mut gram = ComplexMatrix::zeros(m, m)?;
            for i in 0..m {
                for j in 0..m {
                    gram[(i, j)] = if self.metric.is_some() {
                        (0.5 * (dot(&dr[i], &weighted[j]) + dot(&dr[j], &weighted[i]))).into()
                    } else {
                        dot(&dr[i], &dr[j]).into()
                    };
                }
            }
            let (values, vectors) = hermitian_eigh(&gram)?;
            let largest = values.iter().copied().fold(0.0, f64::max);
            let rhs: Vec<_> = dr
                .iter()
                .map(|d| {
                    dot(
                        d,
                        if self.metric.is_some() {
                            weighted_r
                        } else {
                            &r
                        },
                    )
                })
                .collect();
            let mut gamma = vec![0.0; m];
            for (k, &value) in values.iter().enumerate() {
                if value > 1e-12 * largest {
                    let projection: f64 =
                        (0..m).map(|i| vectors[(i, k)].re * rhs[i]).sum::<f64>() / value;
                    for i in 0..m {
                        gamma[i] += vectors[(i, k)].re * projection;
                    }
                }
            }
            for j in 0..m {
                for (i, n) in next.iter_mut().enumerate() {
                    *n -= gamma[j] * (dx[j][i] + self.beta * dr[j][i]);
                }
            }
        }
        charge(next, mean)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn linear_rate_and_charge() {
        let mut mixer = LinearMixer { beta: 0.4 };
        let mut x = vec![1.4, 0.6];
        for k in 1..=12 {
            let out: Vec<_> = x.iter().map(|v| 1.0 + 0.5 * (v - 1.0)).collect();
            x = mixer.mix(&x, &out, 1.0).unwrap();
            assert!((x[0] - 1.0 - 0.4 * 0.8_f64.powi(k)).abs() < 2e-15);
            assert!((x.iter().sum::<f64>() - 2.0).abs() < 1e-14);
        }
        println!(
            "linear: 12-step error={:.16e}, expected rate=0.8",
            x[0] - 1.0
        );
    }
    #[test]
    fn charge_projection_and_rank_deficient_history() {
        let mut linear = LinearMixer { beta: 0.4 };
        let mut pulay = PulayMixer::default();
        for mixer in [&mut linear as &mut dyn DensityMixer, &mut pulay] {
            let n = mixer.mix(&[0.8, 1.2, 1.0], &[1.1, 1.4, 1.7], 1.0).unwrap();
            let defect = (n.iter().sum::<f64>() - 3.0).abs();
            println!("mixer prescribed charge defect={defect:.9e} electrons (unit voxel volume)");
            assert!(defect < 1e-14);
        }
        pulay.reset();
        let mut x = vec![1.3, 0.7];
        for _ in 0..6 {
            let y: Vec<_> = x.iter().map(|v| 1.0 + 0.5 * (v - 1.0)).collect();
            x = pulay.mix(&x, &y, 1.0).unwrap();
        }
        assert!(x.iter().all(|v| (v - 1.0).abs() < 1e-13));
        pulay.reset();
        let a = pulay.mix(&[1.3, 0.7], &[1.15, 0.85], 1.0).unwrap();
        let b = PulayMixer::default()
            .mix(&[1.3, 0.7], &[1.15, 0.85], 1.0)
            .unwrap();
        assert_eq!(a, b);
        assert!(PulayMixer::new(0.3, usize::MAX).is_err());
    }
    #[test]
    fn pulay_linear_finite_termination_reset_and_singular() {
        let mut p = PulayMixer::new(0.4, 6).unwrap();
        let mut x = vec![1.4, 0.6, 1.2, 0.8];
        // Two independent zero-mean modes: fixed-point dimension two.
        for _ in 0..3 {
            let out = vec![
                1.0 + 0.2 * (x[0] - 1.0),
                1.0 + 0.2 * (x[1] - 1.0),
                1.0 + 0.6 * (x[2] - 1.0),
                1.0 + 0.6 * (x[3] - 1.0),
            ];
            x = p.mix(&x, &out, 1.0).unwrap();
            assert!((x.iter().sum::<f64>() - 4.0).abs() < 1e-14);
        }
        let error = x.iter().map(|v| (v - 1.0).abs()).fold(0.0, f64::max);
        println!("Pulay: dimension=2, steps=3, error={error:.16e}");
        assert!(error < 1e-12);
        p.reset();
        assert_eq!(p.history_len(), 0);
        for _ in 0..10 {
            assert_eq!(p.mix(&[1.0; 4], &[1.0; 4], 1.0).unwrap(), vec![1.0; 4]);
        }
        assert_eq!(p.history_len(), 7);
        p.mix(&[2.0; 2], &[2.0; 2], 2.0).unwrap();
        assert_eq!(p.history_len(), 1);
        assert!(p.mix(&[], &[], 1.0).is_err());
    }
}

/// Versioned, backend-independent Anderson memory. Extraction is explicit and
/// never called by the ordinary mixing path. Pair order is oldest to newest.
#[derive(Clone, Debug)]
pub struct PulayState {
    pub beta: f64,
    pub capacity: usize,
    pub inputs: Vec<Vec<f64>>,
    pub residuals: Vec<Vec<f64>>,
    pub mean: Option<f64>,
}
impl PulayMixer {
    pub fn checkpoint(&self) -> PulayState {
        PulayState {
            beta: self.beta,
            capacity: self.capacity,
            inputs: self.pairs.iter().map(|p| p.0.clone()).collect(),
            residuals: self.pairs.iter().map(|p| p.1.clone()).collect(),
            mean: self.mean,
        }
    }
    pub fn restore(&mut self, state: &PulayState, size: usize) -> Result<(), SolverError> {
        if state.beta.to_bits() != self.beta.to_bits()
            || state.capacity != self.capacity
            || state.inputs.len() != state.residuals.len()
            || state.inputs.len() > self.capacity.saturating_add(1)
            || state.mean.is_some_and(|v| !v.is_finite() || v <= 0.0)
            || (!state.inputs.is_empty() && state.mean.is_none())
            || state
                .inputs
                .iter()
                .chain(&state.residuals)
                .any(|v| v.len() != size || v.iter().any(|x| !x.is_finite()))
        {
            return Err(SolverError::InvalidParameter("invalid Pulay checkpoint"));
        }
        self.pairs = state
            .inputs
            .iter()
            .cloned()
            .zip(state.residuals.iter().cloned())
            .collect();
        self.mean = state.mean;
        Ok(())
    }
}

#[cfg(test)]
mod metric_tests {
    use super::*;
    #[derive(Debug, Clone)]
    struct DiagonalMetric([f64; 3]);
    impl ResidualMetric for DiagonalMetric {
        fn duplicate(&self) -> Box<dyn ResidualMetric> {
            Box::new(self.clone())
        }
        fn apply(&mut self, r: &[f64], out: &mut [f64]) -> Result<(), SolverError> {
            for ((v, r), w) in out.iter_mut().zip(r).zip(self.0) {
                *v = r * w;
            }
            Ok(())
        }
    }
    #[test]
    fn metric_changes_only_the_pulay_least_squares_and_restores_raw_history() {
        let mut p = PulayMixer::new(0.3, 2).unwrap();
        p.set_metric(Box::new(DiagonalMetric([1., 2., 4.])));
        let x0 = [2., 2., 2.];
        let r0 = [0.2, -0.1, -0.1];
        let x1 = [2.03, 1.98, 1.99];
        let r1 = [-0.1, 0.2, -0.1];
        let y0: Vec<_> = x0.iter().zip(r0).map(|(x, r)| x + r).collect();
        let y1: Vec<_> = x1.iter().zip(r1).map(|(x, r)| x + r).collect();
        p.mix(&x0, &y0, 2.).unwrap();
        let state = p.checkpoint();
        let result = p.mix(&x1, &y1, 2.).unwrap();
        let d = [-0.3, 0.3, 0.];
        let gamma = (d[0] * r1[0] + 2. * d[1] * r1[1]) / (d[0] * d[0] + 2. * d[1] * d[1]);
        for i in 0..3 {
            let expected = x1[i] + 0.3 * r1[i] - (x1[i] - x0[i] + 0.3 * d[i]) * gamma;
            assert!((result[i] - expected).abs() < 1e-14);
        }
        p.restore(&state, 3).unwrap();
        assert_eq!(result, p.mix(&x1, &y1, 2.).unwrap());
    }
}
