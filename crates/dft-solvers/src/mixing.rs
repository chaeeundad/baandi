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
    /// Seeded difference pairs (dx, dr), oldest first, from another map (the
    /// cutoff-ladder hand-off). They are stored as differences and form a closed
    /// block: no difference is ever formed between a seed and an absolute pair,
    /// so a cross-map difference never enters the Gram matrix. They precede the
    /// window differences and age out first (seeds + window <= capacity).
    seeds: Vec<(Vec<f64>, Vec<f64>)>,
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
            seeds: self.seeds.clone(),
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
            seeds: Vec::new(),
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
    /// Number of seeded difference pairs still in the window.
    pub fn seeded_len(&self) -> usize {
        self.seeds.len()
    }
    /// Seed an EMPTY history with difference pairs (dx_i, dr_i), oldest first,
    /// taken from a related map (secant recycling across the cutoff-ladder
    /// hand-off). `mean` is the prescribed density mean of the coming `mix`
    /// calls. At most `capacity` pairs are kept (the newest). The seeds enter the
    /// least squares like ordinary differences (same spectral cut) and are
    /// dropped oldest first as new differences arrive.
    pub fn seed_differences(
        &mut self,
        dx: &[Vec<f64>],
        dr: &[Vec<f64>],
        mean: f64,
    ) -> Result<(), SolverError> {
        let size = dx.first().map_or(0, Vec::len);
        if !self.pairs.is_empty()
            || !self.seeds.is_empty()
            || dx.len() != dr.len()
            || !mean.is_finite()
            || mean <= 0.0
            || size == 0
            || dx
                .iter()
                .chain(dr)
                .any(|v| v.len() != size || v.iter().any(|x| !x.is_finite()))
        {
            return Err(SolverError::InvalidParameter(
                "seed_differences: empty history, equal finite pairs and positive mean required",
            ));
        }
        let skip = dx.len().saturating_sub(self.capacity);
        self.seeds = dx[skip..]
            .iter()
            .cloned()
            .zip(dr[skip..].iter().cloned())
            .collect();
        self.mean = Some(mean);
        Ok(())
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
        self.seeds.clear();
        self.mean = None;
    }
    fn mix(&mut self, input: &[f64], output: &[f64], mean: f64) -> Result<Vec<f64>, SolverError> {
        validate(input, output, mean, self.beta)?;
        if self.mean != Some(mean)
            || self.pairs.first().is_some_and(|p| p.0.len() != input.len())
            || self.seeds.first().is_some_and(|p| p.0.len() != input.len())
        {
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
        while !self.seeds.is_empty() && self.seeds.len() + self.pairs.len() - 1 > self.capacity {
            self.seeds.remove(0);
        }
        let mut next: Vec<_> = input
            .iter()
            .zip(&r)
            .map(|(a, r)| a + self.beta * r)
            .collect();
        let m = self.seeds.len() + self.pairs.len() - 1;
        if m > 0 {
            // Seeded differences first (older), then the window differences of the
            // absolute pairs; never a difference between the two blocks.
            let dx: Vec<Vec<f64>> = self
                .seeds
                .iter()
                .map(|s| s.0.clone())
                .chain(
                    self.pairs
                        .windows(2)
                        .map(|p| p[1].0.iter().zip(&p[0].0).map(|(a, b)| a - b).collect()),
                )
                .collect();
            let dr: Vec<Vec<f64>> = self
                .seeds
                .iter()
                .map(|s| s.1.clone())
                .chain(
                    self.pairs
                        .windows(2)
                        .map(|p| p[1].1.iter().zip(&p[0].1).map(|(a, b)| a - b).collect()),
                )
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
    /// Seeded difference pairs (oldest first); empty unless seeded.
    pub seed_dx: Vec<Vec<f64>>,
    pub seed_dr: Vec<Vec<f64>>,
}
impl PulayMixer {
    pub fn checkpoint(&self) -> PulayState {
        PulayState {
            beta: self.beta,
            capacity: self.capacity,
            inputs: self.pairs.iter().map(|p| p.0.clone()).collect(),
            residuals: self.pairs.iter().map(|p| p.1.clone()).collect(),
            mean: self.mean,
            seed_dx: self.seeds.iter().map(|p| p.0.clone()).collect(),
            seed_dr: self.seeds.iter().map(|p| p.1.clone()).collect(),
        }
    }
    pub fn restore(&mut self, state: &PulayState, size: usize) -> Result<(), SolverError> {
        if state.beta.to_bits() != self.beta.to_bits()
            || state.capacity != self.capacity
            || state.inputs.len() != state.residuals.len()
            || state.inputs.len() > self.capacity.saturating_add(1)
            || state.mean.is_some_and(|v| !v.is_finite() || v <= 0.0)
            || ((!state.inputs.is_empty() || !state.seed_dx.is_empty()) && state.mean.is_none())
            || state.seed_dx.len() != state.seed_dr.len()
            || (!state.seed_dx.is_empty()
                && state.seed_dx.len() + state.inputs.len().saturating_sub(1) > self.capacity)
            || state
                .inputs
                .iter()
                .chain(&state.residuals)
                .chain(&state.seed_dx)
                .chain(&state.seed_dr)
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
        self.seeds = state
            .seed_dx
            .iter()
            .cloned()
            .zip(state.seed_dr.iter().cloned())
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

#[cfg(test)]
mod seed_tests {
    use super::*;
    /// Linear fixed-point map with three zero-mean modes (rates 0.2, 0.6, 0.9).
    fn map(x: &[f64]) -> Vec<f64> {
        let rates = [0.2, 0.6, 0.9, 0.9];
        let mut y: Vec<_> = x
            .iter()
            .zip(rates)
            .map(|(v, a)| 1.0 + a * (v - 1.0))
            .collect();
        // Keep the charge: the last sample absorbs the mean change.
        let defect = y.iter().sum::<f64>() - x.iter().sum::<f64>();
        y[3] -= defect;
        y
    }
    fn residual(x: &[f64]) -> Vec<f64> {
        map(x).iter().zip(x).map(|(a, b)| a - b).collect()
    }
    fn error(x: &[f64]) -> f64 {
        x.iter().map(|v| (v - 1.0).abs()).fold(0.0, f64::max)
    }

    #[test]
    fn exact_seeds_of_the_same_map_reach_the_fixed_point_in_one_step() {
        // Differences of the SAME linear map span its three modes, so the first
        // seeded least squares is exact (Anderson on a linear map).
        let probes = [
            [1.3, 0.8, 1.1, 0.8],
            [1.1, 1.2, 0.7, 1.0],
            [0.9, 1.0, 1.25, 0.85],
            [1.05, 0.9, 1.0, 1.05],
        ];
        let dx: Vec<Vec<f64>> = probes
            .windows(2)
            .map(|p| p[1].iter().zip(&p[0]).map(|(a, b)| a - b).collect())
            .collect();
        let dr: Vec<Vec<f64>> = probes
            .windows(2)
            .map(|p| {
                residual(&p[1])
                    .iter()
                    .zip(residual(&p[0]))
                    .map(|(a, b)| a - b)
                    .collect()
            })
            .collect();
        let start = [1.2, 0.85, 0.9, 1.05];
        let mut seeded = PulayMixer::new(0.3, 6).unwrap();
        seeded.seed_differences(&dx, &dr, 1.0).unwrap();
        assert_eq!(seeded.seeded_len(), 3);
        let one = seeded.mix(&start, &map(&start), 1.0).unwrap();
        let mut cold = PulayMixer::new(0.3, 6).unwrap();
        let cold_one = cold.mix(&start, &map(&start), 1.0).unwrap();
        println!(
            "seeded first-step error {:.3e}, cold {:.3e}",
            error(&one),
            error(&cold_one)
        );
        assert!(error(&one) < 1e-12);
        assert!(error(&cold_one) > 1e-2);
        assert!((one.iter().sum::<f64>() - 4.0).abs() < 1e-14);
    }

    #[test]
    fn seeds_never_difference_against_new_pairs_and_age_out() {
        // A seed that would be catastrophic if differenced against the first
        // fine pair (huge absolute offset) has no effect: only its own dx, dr enter.
        let mut p = PulayMixer::new(0.3, 3).unwrap();
        let dx = vec![vec![0.1, -0.1, 0.0, 0.0], vec![0.0, 0.05, -0.05, 0.0]];
        let dr = vec![vec![-0.08, 0.06, 0.02, 0.0], vec![0.0, -0.02, 0.0, 0.02]];
        p.seed_differences(&dx, &dr, 1.0).unwrap();
        let x0 = [1.2, 0.85, 0.9, 1.05];
        let r0 = residual(&x0);
        let out = p.mix(&x0, &map(&x0), 1.0).unwrap();
        // Hand-computed least squares on the two seeds only.
        let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
        let g = [
            [dot(&dr[0], &dr[0]), dot(&dr[0], &dr[1])],
            [dot(&dr[1], &dr[0]), dot(&dr[1], &dr[1])],
        ];
        let b = [dot(&dr[0], &r0), dot(&dr[1], &r0)];
        let det = g[0][0] * g[1][1] - g[0][1] * g[1][0];
        let gamma = [
            (b[0] * g[1][1] - g[0][1] * b[1]) / det,
            (g[0][0] * b[1] - g[1][0] * b[0]) / det,
        ];
        for i in 0..4 {
            let expected = x0[i] + 0.3 * r0[i]
                - gamma[0] * (dx[0][i] + 0.3 * dr[0][i])
                - gamma[1] * (dx[1][i] + 0.3 * dr[1][i]);
            assert!(
                (out[i] - expected).abs() < 1e-13,
                "{i}: {} {expected}",
                out[i]
            );
        }
        // Capacity 3: after the second fine pair (one fine difference) both seeds
        // fit; the third pair (two differences) drops the oldest seed, etc.
        let mut x = out;
        let mut seen = vec![p.seeded_len()];
        for _ in 0..4 {
            let y = map(&x);
            x = p.mix(&x, &y, 1.0).unwrap();
            seen.push(p.seeded_len());
        }
        assert_eq!(seen, vec![2, 2, 1, 0, 0]);
        // Checkpoint/restore carries the seeds exactly.
        let mut q = PulayMixer::new(0.3, 3).unwrap();
        q.seed_differences(&dx, &dr, 1.0).unwrap();
        let state = q.checkpoint();
        assert_eq!(state.seed_dx, dx);
        let mut restored = PulayMixer::new(0.3, 3).unwrap();
        restored.restore(&state, 4).unwrap();
        assert_eq!(
            q.mix(&x0, &map(&x0), 1.0).unwrap(),
            restored.mix(&x0, &map(&x0), 1.0).unwrap()
        );
        // Seeding needs an empty history; reset and a changed mean clear seeds.
        assert!(q.seed_differences(&dx, &dr, 1.0).is_err());
        q.reset();
        assert_eq!(q.seeded_len(), 0);
        q.seed_differences(&dx, &dr, 1.0).unwrap();
        q.mix(&[2.0; 4], &[2.0; 4], 2.0).unwrap();
        assert_eq!(q.seeded_len(), 0);
        // Mismatched or invalid seeds reject.
        let mut r = PulayMixer::new(0.3, 3).unwrap();
        assert!(r.seed_differences(&dx, &dr[..1], 1.0).is_err());
        assert!(r.seed_differences(&dx, &dr, 0.0).is_err());
        // More seeds than capacity keep the newest.
        let many = vec![dx[0].clone(); 5];
        r.seed_differences(&many, &many, 1.0).unwrap();
        assert_eq!(r.seeded_len(), 3);
    }
}
