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
    /// `BAANDI_MIX_INCREMENTAL=on` (ticket MPI-REMAINING): keep the differences
    /// and Gram entries across calls (bit-identical; see `Incremental`).
    incremental: Option<Incremental>,
    /// Diagnostics of the last `mix` call (read-only; never feeds back).
    last: Option<PulayStepDiagnostics>,
    /// SCF-TRAJECTORY Coulomb guard (`BAANDI_MIX_COULOMB_GUARD`): the metric
    /// in which the L2 fit's predicted residual is checked, and whether a
    /// predicted growth refits the step in that metric (false = trace only).
    guard: Option<(Box<dyn ResidualMetric>, bool)>,
}
/// What the last Pulay update did (SCF-TRAJECTORY research output). Reading it
/// never changes a proposal; it costs three grid dot products per call.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PulayStepDiagnostics {
    /// Difference pairs in the least squares (seeds included); 0 = plain step.
    pub differences: usize,
    /// Gram eigen-directions kept by the spectral cut (<= differences).
    pub kept: usize,
    /// Largest / smallest Gram eigenvalue (fit metric); infinite when the
    /// smallest is not positive, 1 without differences.
    pub gram_condition: f64,
    /// Extrapolation coefficients gamma, oldest difference first.
    pub coefficients: Vec<f64>,
    /// L2 norm (plain sum over grid points) of the residual R = out - in.
    pub residual_norm: f64,
    /// Fit-metric norm of the predicted residual R - dR gamma divided by the
    /// fit-metric norm of R (1 without differences; 0 = R fully explained).
    pub predicted_ratio: f64,
    /// L2 norm (same sum) of the proposed step next - input, before the
    /// charge projection; step_norm / residual_norm = beta for a plain step.
    pub step_norm: f64,
    /// Coulomb guard (trace or on): sqrt(P^T C P / R^T C R) with P = R - dR gamma
    /// of the L2 fit and C the guard metric; None without a guard or differences.
    pub coulomb_predicted_ratio: Option<f64>,
    /// The guard refitted this step in the Coulomb metric (ratio > 1, mode on).
    pub coulomb_refit: bool,
}
/// Least-squares summary of one update, for `PulayStepDiagnostics`.
struct FitSummary {
    kept: usize,
    smallest: f64,
    largest: f64,
    gamma: Vec<f64>,
    gamma_rhs: f64,
    rr_fit: f64,
}
/// Cross-call memory of the incremental Pulay update. Every pair and seed gets
/// a serial number when it enters the history; a difference is identified by
/// the serial numbers it is formed from, so a cached difference (and its metric
/// image) is reused only for exactly the same operands, and a Gram entry only
/// for exactly the same two differences. Formulas and summation orders are the
/// ones of the plain path, so results are bit-identical.
#[derive(Debug, Clone, Default)]
struct Incremental {
    next_serial: u64,
    pair_serials: Vec<u64>,
    seed_serials: Vec<u64>,
    /// Cached differences in the order of the current history (seeds first).
    diffs: Vec<CachedDiff>,
    /// Gram entries keyed by the two difference keys.
    gram: Vec<((DiffKey, DiffKey), f64)>,
    /// Buffers of dropped differences, reused for new ones.
    spare: Vec<CachedDiff>,
}
/// (first operand serial, second operand serial); a seed is (serial, u64::MAX).
type DiffKey = (u64, u64);
#[derive(Debug, Clone, Default)]
struct CachedDiff {
    key: DiffKey,
    dx: Vec<f64>,
    dr: Vec<f64>,
    weighted: Vec<f64>,
}
/// Process-wide counters of the Pulay mixer (ticket MPI-REMAINING; relaxed
/// atomics, no allocation): mix calls with history, differences formed (two
/// grid vectors each), Gram dot products, metric applications.
static MIX_CALLS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
static MIX_DIFFS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
static MIX_DOTS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
static MIX_METRIC: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
/// (mix calls with history, differences formed, Gram and right-hand-side dot
/// products, metric applications) since process start.
pub fn mixer_counters() -> [u64; 4] {
    use std::sync::atomic::Ordering::Relaxed;
    [
        MIX_CALLS.load(Relaxed),
        MIX_DIFFS.load(Relaxed),
        MIX_DOTS.load(Relaxed),
        MIX_METRIC.load(Relaxed),
    ]
}
fn bump(c: &std::sync::atomic::AtomicU64, n: usize) {
    c.fetch_add(n as u64, std::sync::atomic::Ordering::Relaxed);
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
            incremental: self.incremental.clone(),
            last: self.last.clone(),
            guard: self.guard.as_ref().map(|(m, act)| (m.duplicate(), *act)),
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
            incremental: incremental_switch().then(Incremental::default),
            last: None,
            guard: None,
        })
    }
    /// SCF-TRAJECTORY Coulomb guard (see `BAANDI_MIX_COULOMB_GUARD`): after the
    /// L2 fit, the predicted residual P = R - dR gamma is measured in `metric`;
    /// with `refit` and P^T C P > R^T C R the step's coefficients are solved
    /// again in that metric (same spectral cut). Only for the unweighted fit
    /// (ignored when a fit metric is set). History is not changed.
    pub fn set_coulomb_guard(&mut self, metric: Box<dyn ResidualMetric>, refit: bool) {
        self.guard = Some((metric, refit));
    }
    /// Diagnostics of the last `mix` call; None before the first call and
    /// after `reset`.
    pub fn last_step(&self) -> Option<&PulayStepDiagnostics> {
        self.last.as_ref()
    }
    /// Force the incremental update on or off (tests; the constructor reads
    /// `BAANDI_MIX_INCREMENTAL`). Clears the history.
    pub fn set_incremental(&mut self, on: bool) {
        self.reset();
        self.incremental = on.then(Incremental::default);
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
        if let Some(inc) = &mut self.incremental {
            inc.forget(self.pairs.len(), self.seeds.len());
        }
        Ok(())
    }
}
thread_local! {
    static INCREMENTAL_OVERRIDE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}
/// Override `BAANDI_MIX_INCREMENTAL` for mixers constructed on the calling
/// thread (simulated ranks in tests); `None` restores the environment.
pub fn set_incremental_override(on: Option<bool>) {
    INCREMENTAL_OVERRIDE.with(|o| o.set(on));
}
fn incremental_switch() -> bool {
    INCREMENTAL_OVERRIDE
        .with(|o| o.get())
        .unwrap_or_else(|| dft_numerics::switches::on(dft_numerics::switches::Id::MIX_INCREMENTAL))
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum::<f64>()
}
impl Incremental {
    fn serial(&mut self) -> u64 {
        self.next_serial += 1;
        self.next_serial
    }
    /// Fresh serial numbers for `pairs` pairs and `seeds` seeds; every cached
    /// difference and Gram entry is dropped (buffers kept for reuse).
    fn forget(&mut self, pairs: usize, seeds: usize) {
        self.pair_serials.clear();
        self.seed_serials.clear();
        for _ in 0..pairs {
            let s = self.serial();
            self.pair_serials.push(s);
        }
        for _ in 0..seeds {
            let s = self.serial();
            self.seed_serials.push(s);
        }
        self.spare.append(&mut self.diffs);
        self.gram.clear();
    }
    /// Renumber when the history changed outside `mix` (seeding, restore).
    fn sync(&mut self, pairs: usize, seeds: usize) {
        if self.pair_serials.len() != pairs || self.seed_serials.len() != seeds {
            self.forget(pairs, seeds);
        }
    }
}
impl PulayMixer {
    /// The plain path's Pulay update with cached differences, metric images and
    /// Gram entries (bit-identical: same operands, formulas and orders).
    fn incremental_update(
        &mut self,
        r: &[f64],
        next: &mut [f64],
        m: usize,
    ) -> Result<FitSummary, SolverError> {
        let Self {
            beta,
            metric,
            metric_r,
            pairs,
            seeds,
            incremental,
            ..
        } = self;
        let inc = incremental.as_mut().expect("incremental state");
        inc.sync(pairs.len(), seeds.len());
        let mut old = std::mem::take(&mut inc.diffs);
        let keys = inc
            .seed_serials
            .iter()
            .map(|&s| (s, u64::MAX))
            .chain(inc.pair_serials.windows(2).map(|w| (w[0], w[1])));
        for (j, key) in keys.enumerate() {
            if let Some(at) = old.iter().position(|d| d.key == key) {
                inc.diffs.push(old.swap_remove(at));
                continue;
            }
            let mut d = inc.spare.pop().unwrap_or_default();
            d.key = key;
            d.dx.clear();
            d.dr.clear();
            if j < seeds.len() {
                d.dx.extend_from_slice(&seeds[j].0);
                d.dr.extend_from_slice(&seeds[j].1);
            } else {
                let w = j - seeds.len();
                let (p0, p1) = (&pairs[w], &pairs[w + 1]);
                d.dx.extend(p1.0.iter().zip(&p0.0).map(|(a, b)| a - b));
                d.dr.extend(p1.1.iter().zip(&p0.1).map(|(a, b)| a - b));
            }
            bump(&MIX_DIFFS, 1);
            if let Some(metric) = metric.as_mut() {
                d.weighted.resize(d.dr.len(), 0.0);
                metric.apply(&d.dr, &mut d.weighted)?;
                bump(&MIX_METRIC, 1);
            }
            inc.diffs.push(d);
        }
        inc.spare.append(&mut old);
        debug_assert_eq!(inc.diffs.len(), m);
        if let Some(metric) = metric.as_mut() {
            metric.apply(r, metric_r)?;
            bump(&MIX_METRIC, 1);
        }
        let weighted_metric = metric.is_some();
        let mut gram = ComplexMatrix::zeros(m, m)?;
        let mut kept = Vec::with_capacity(m * m);
        for i in 0..m {
            for j in 0..m {
                let key = (inc.diffs[i].key, inc.diffs[j].key);
                let value = match inc.gram.iter().find(|(k, _)| *k == key) {
                    Some(&(_, v)) => v,
                    None => {
                        let (a, b) = (&inc.diffs[i], &inc.diffs[j]);
                        bump(&MIX_DOTS, if weighted_metric { 2 } else { 1 });
                        if weighted_metric {
                            0.5 * (dot(&a.dr, &b.weighted) + dot(&b.dr, &a.weighted))
                        } else {
                            dot(&a.dr, &b.dr)
                        }
                    }
                };
                gram[(i, j)] = value.into();
                kept.push((key, value));
            }
        }
        inc.gram = kept;
        let (values, vectors) = hermitian_eigh(&gram)?;
        let largest = values.iter().copied().fold(0.0, f64::max);
        let rhs: Vec<_> = inc
            .diffs
            .iter()
            .map(|d| {
                dot(
                    &d.dr,
                    if weighted_metric {
                        metric_r.as_slice()
                    } else {
                        r
                    },
                )
            })
            .collect();
        bump(&MIX_DOTS, m);
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
        for (d, g) in inc.diffs.iter().zip(&gamma) {
            for (i, n) in next.iter_mut().enumerate() {
                *n -= g * (d.dx[i] + *beta * d.dr[i]);
            }
        }
        let rr_fit = dot(
            r,
            if weighted_metric {
                metric_r.as_slice()
            } else {
                r
            },
        );
        Ok(FitSummary::new(&values, largest, gamma, &rhs, rr_fit))
    }
}
impl FitSummary {
    fn new(values: &[f64], largest: f64, gamma: Vec<f64>, rhs: &[f64], rr_fit: f64) -> Self {
        Self {
            kept: values.iter().filter(|v| **v > 1e-12 * largest).count(),
            smallest: values.iter().copied().fold(f64::INFINITY, f64::min),
            largest,
            gamma_rhs: gamma.iter().zip(rhs).map(|(g, b)| g * b).sum(),
            gamma,
            rr_fit,
        }
    }
}
impl PulayMixer {
    /// Differences (dx, dr) of the current history, seeds first, then the
    /// window differences (the plain path's formulas and order).
    fn differences(&self) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
        let window = |k: usize| {
            self.pairs.windows(2).map(move |p| {
                p[1].0
                    .iter()
                    .zip(&p[0].0)
                    .zip(p[1].1.iter().zip(&p[0].1))
                    .map(|((a, b), (c, d))| if k == 0 { a - b } else { c - d })
                    .collect::<Vec<f64>>()
            })
        };
        let dx = self
            .seeds
            .iter()
            .map(|s| s.0.clone())
            .chain(window(0))
            .collect();
        let dr = self
            .seeds
            .iter()
            .map(|s| s.1.clone())
            .chain(window(1))
            .collect();
        (dx, dr)
    }
    /// SCF-TRAJECTORY Coulomb guard of an unweighted fit: returns the predicted
    /// Coulomb ratio and whether the step was refitted (then `next` and `fit`
    /// are replaced by the Coulomb-metric solution). Runs once per SCF map,
    /// outside every inner loop; allocates m + 3 grid vectors.
    fn coulomb_guard(
        &mut self,
        input: &[f64],
        r: &[f64],
        next: &mut [f64],
        fit: &mut FitSummary,
    ) -> Result<(f64, bool), SolverError> {
        let (dx, dr) = self.differences();
        let beta = self.beta;
        let mut pred = r.to_vec();
        for (d, g) in dr.iter().zip(&fit.gamma) {
            for (p, v) in pred.iter_mut().zip(d) {
                *p -= g * v;
            }
        }
        let (metric, refit) = self.guard.as_mut().expect("guard set");
        let mut cr = vec![0.0; r.len()];
        metric.apply(r, &mut cr)?;
        let mut cp = vec![0.0; r.len()];
        metric.apply(&pred, &mut cp)?;
        let (qr, qp) = (dot(r, &cr), dot(&pred, &cp));
        let ratio = if qr > 0.0 {
            (qp.max(0.0) / qr).sqrt()
        } else {
            0.0
        };
        if !*refit || qp <= qr || qr <= 0.0 {
            return Ok((ratio, false));
        }
        let m = dr.len();
        let mut weighted = vec![vec![0.0; r.len()]; m];
        for (s, w) in dr.iter().zip(weighted.iter_mut()) {
            metric.apply(s, w)?;
        }
        let mut gram = ComplexMatrix::zeros(m, m)?;
        for i in 0..m {
            for j in 0..m {
                gram[(i, j)] =
                    (0.5 * (dot(&dr[i], &weighted[j]) + dot(&dr[j], &weighted[i]))).into();
            }
        }
        let (values, vectors) = hermitian_eigh(&gram)?;
        let largest = values.iter().copied().fold(0.0, f64::max);
        let rhs: Vec<f64> = dr.iter().map(|d| dot(d, &cr)).collect();
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
        for (i, n) in next.iter_mut().enumerate() {
            *n = input[i] + beta * r[i];
        }
        for j in 0..m {
            for (i, n) in next.iter_mut().enumerate() {
                *n -= gamma[j] * (dx[j][i] + beta * dr[j][i]);
            }
        }
        *fit = FitSummary::new(&values, largest, gamma, &rhs, qr);
        Ok((ratio, true))
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
        self.last = None;
        if let Some(inc) = &mut self.incremental {
            inc.forget(0, 0);
        }
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
        if let Some(inc) = &mut self.incremental {
            inc.sync(self.pairs.len() - 1, self.seeds.len());
            let serial = inc.serial();
            inc.pair_serials.push(serial);
        }
        if self.pairs.len() > self.capacity + 1 {
            self.pairs.remove(0);
            if let Some(inc) = &mut self.incremental {
                inc.pair_serials.remove(0);
            }
        }
        while !self.seeds.is_empty() && self.seeds.len() + self.pairs.len() - 1 > self.capacity {
            self.seeds.remove(0);
            if let Some(inc) = &mut self.incremental {
                inc.seed_serials.remove(0);
            }
        }
        let mut next: Vec<_> = input
            .iter()
            .zip(&r)
            .map(|(a, r)| a + self.beta * r)
            .collect();
        let m = self.seeds.len() + self.pairs.len() - 1;
        if m > 0 {
            bump(&MIX_CALLS, 1);
        }
        let mut fit = None;
        if m > 0 && self.incremental.is_some() {
            fit = Some(self.incremental_update(&r, &mut next, m)?);
        } else if m > 0 {
            bump(&MIX_DIFFS, m);
            bump(&MIX_DOTS, m * m + m);
            if self.metric.is_some() {
                bump(&MIX_METRIC, m + 1);
            }
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
            let rr_fit = dot(
                &r,
                if self.metric.is_some() {
                    weighted_r.as_slice()
                } else {
                    r.as_slice()
                },
            );
            fit = Some(FitSummary::new(&values, largest, gamma, &rhs, rr_fit));
        }
        let mut coulomb = (None, false);
        if self.guard.is_some()
            && self.metric.is_none()
            && let Some(f) = fit.as_mut()
        {
            let (ratio, refit) = self.coulomb_guard(input, &r, &mut next, f)?;
            coulomb = (Some(ratio), refit);
        }
        let residual_norm = r.iter().map(|v| v * v).sum::<f64>().sqrt();
        let step_norm = next
            .iter()
            .zip(input)
            .map(|(n, x)| (n - x) * (n - x))
            .sum::<f64>()
            .sqrt();
        self.last = Some(match fit {
            Some(f) => PulayStepDiagnostics {
                differences: m,
                kept: f.kept,
                gram_condition: if f.smallest > 0.0 {
                    f.largest / f.smallest
                } else {
                    f64::INFINITY
                },
                coefficients: f.gamma,
                residual_norm,
                predicted_ratio: if f.rr_fit > 0.0 {
                    ((f.rr_fit - f.gamma_rhs).max(0.0) / f.rr_fit).sqrt()
                } else {
                    0.0
                },
                step_norm,
                coulomb_predicted_ratio: coulomb.0,
                coulomb_refit: coulomb.1,
            },
            None => PulayStepDiagnostics {
                differences: 0,
                kept: 0,
                gram_condition: 1.0,
                coefficients: Vec::new(),
                residual_norm,
                predicted_ratio: 1.0,
                step_norm,
                coulomb_predicted_ratio: None,
                coulomb_refit: false,
            },
        });
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
        if let Some(inc) = &mut self.incremental {
            inc.forget(self.pairs.len(), self.seeds.len());
        }
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

#[cfg(test)]
mod incremental_tests {
    use super::*;
    /// Nonlinear charge-keeping map on a small grid (deterministic).
    fn map(x: &[f64]) -> Vec<f64> {
        let n = x.len();
        let mut y: Vec<_> = (0..n)
            .map(|i| 1.0 + 0.7 * (x[(i + 1) % n] - 1.0) + 0.05 * (x[i] - 1.0).powi(2))
            .collect();
        let defect = (y.iter().sum::<f64>() - x.iter().sum::<f64>()) / n as f64;
        for v in &mut y {
            *v -= defect;
        }
        y
    }
    #[derive(Debug, Clone)]
    struct Weights(Vec<f64>);
    impl ResidualMetric for Weights {
        fn duplicate(&self) -> Box<dyn ResidualMetric> {
            Box::new(self.clone())
        }
        fn apply(&mut self, r: &[f64], out: &mut [f64]) -> Result<(), SolverError> {
            for ((v, r), w) in out.iter_mut().zip(r).zip(&self.0) {
                *v = r * w;
            }
            Ok(())
        }
    }
    fn run(incremental: bool, metric: bool, history: usize, script: &[u8]) -> Vec<Vec<u64>> {
        let n = 17;
        let mut p = PulayMixer::new(0.4, history).unwrap();
        if metric {
            p.set_metric(Box::new(Weights(
                (0..n).map(|i| 1.0 + i as f64 * 0.1).collect(),
            )));
        }
        p.set_incremental(incremental);
        let mut x: Vec<f64> = (0..n)
            .map(|i| 1.0 + 0.3 * ((i * 7 % 5) as f64 - 2.0) / 2.0)
            .collect();
        let mean = x.iter().sum::<f64>() / n as f64;
        let mut out = Vec::new();
        let mut saved = None;
        for &op in script {
            match op {
                b'm' => {
                    x = p.mix(&x, &map(&x), mean).unwrap();
                    out.push(x.iter().map(|v| v.to_bits()).collect());
                }
                b'c' => saved = Some(p.checkpoint()),
                b'r' => p.restore(saved.as_ref().unwrap(), n).unwrap(),
                b'z' => p.reset(),
                b's' => {
                    p.reset();
                    let dx: Vec<Vec<f64>> = (0..3)
                        .map(|k| {
                            (0..n)
                                .map(|i| 1e-3 * ((i + k) % 4) as f64 - 1.5e-3)
                                .collect()
                        })
                        .collect();
                    let dr: Vec<Vec<f64>> = (0..3)
                        .map(|k| {
                            (0..n)
                                .map(|i| 2e-3 * ((i * k + 1) % 3) as f64 - 2e-3)
                                .collect()
                        })
                        .collect();
                    p.seed_differences(&dx, &dr, mean).unwrap();
                }
                _ => unreachable!(),
            }
        }
        out
    }
    #[test]
    fn incremental_update_is_bit_identical_to_the_plain_update() {
        // Window filling and sliding, checkpoint/restore with the same pair
        // count, reset, seeds that age out; l2 and a metric; capacities 2 and 6.
        let scripts: [&[u8]; 4] = [
            b"mmmmmmmmmmmmmm",
            b"mmmmcmmmrmmmm",
            b"mmmzmmmmmmm",
            b"smmmmmmmmsmmm",
        ];
        for script in scripts {
            for metric in [false, true] {
                for history in [2, 6] {
                    let plain = run(false, metric, history, script);
                    let fast = run(true, metric, history, script);
                    assert_eq!(
                        plain,
                        fast,
                        "script {:?} metric {metric} history {history}",
                        std::str::from_utf8(script).unwrap()
                    );
                }
            }
        }
    }
    #[test]
    fn incremental_update_keeps_one_difference_per_history_entry() {
        let mut p = PulayMixer::new(0.4, 6).unwrap();
        p.set_incremental(true);
        let n = 9;
        let mut x: Vec<f64> = (0..n).map(|i| 1.0 + 0.1 * (i as f64 - 4.0) / 4.0).collect();
        let mean = x.iter().sum::<f64>() / n as f64;
        for _ in 0..9 {
            x = p.mix(&x, &map(&x), mean).unwrap();
        }
        // Process-wide counters are shared with concurrent tests; the cache
        // sizes are exact: six differences, a 6 x 6 Gram, at most one spare.
        let inc = p.incremental.as_ref().unwrap();
        assert_eq!(inc.diffs.len(), 6);
        assert_eq!(inc.gram.len(), 36);
        assert!(inc.spare.len() <= 1);
    }

    /// Per-entry weights (a stand-in for the Coulomb metric's 4 pi / G^2).
    #[derive(Debug, Clone)]
    struct Weights4([f64; 4]);
    impl ResidualMetric for Weights4 {
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
    /// Two zero-mean modes: u (entries 0, 1) is screened (output fixed, linear
    /// mixing at 0.5 halves it), v (entries 2, 3) is over-responding (output
    /// deviation -3 v, linear mixing at 0.5 flips its sign: charge sloshing).
    /// With u >> v the one-difference L2 fit gives gamma near -1 and predicts a
    /// tripled v residual: the SCF-TRAJECTORY bounce in two numbers.
    fn sloshing(x: &[f64]) -> Vec<f64> {
        let v = (x[2] - x[3]) / 2.0;
        vec![1.0, 1.0, 1.0 - 3.0 * v, 1.0 + 3.0 * v]
    }
    fn start() -> Vec<f64> {
        vec![1.0 + 1e-2, 1.0 - 1e-2, 1.0 + 1e-4, 1.0 - 1e-4]
    }
    #[test]
    fn diagnostics_describe_the_step_without_changing_it() {
        let mut p = PulayMixer::new(0.5, 6).unwrap();
        assert!(p.last_step().is_none());
        let x0 = start();
        let x1 = p.mix(&x0, &sloshing(&x0), 1.0).unwrap();
        let d = p.last_step().unwrap().clone();
        assert_eq!((d.differences, d.kept, d.coefficients.len()), (0, 0, 0));
        assert_eq!(d.predicted_ratio, 1.0);
        assert!((d.step_norm / d.residual_norm - 0.5).abs() < 1e-14);
        assert!(d.coulomb_predicted_ratio.is_none() && !d.coulomb_refit);
        p.mix(&x1, &sloshing(&x1), 1.0).unwrap();
        let d = p.last_step().unwrap();
        assert_eq!((d.differences, d.kept), (1, 1));
        // Closed form of the one-difference fit: dR = (u/2, 8v) per mode
        // entry, R = (-u/2, 4v): gamma = (-u^2/4 + 32 v^2) / (u^2/4 + 64 v^2).
        let (u, v) = (1e-2f64, 1e-4f64);
        let exact = (-u * u / 4.0 + 32.0 * v * v) / (u * u / 4.0 + 64.0 * v * v);
        assert!(
            (d.coefficients[0] - exact).abs() < 1e-12,
            "gamma {:?}",
            d.coefficients
        );
        assert!(d.predicted_ratio < 1.0 && d.gram_condition == 1.0);
        // The L2 step doubles the plain one (the bounce precursor).
        assert!(d.step_norm / d.residual_norm > 0.9);
        p.reset();
        assert!(p.last_step().is_none());
    }
    #[test]
    fn coulomb_guard_trace_is_inert_and_on_refits_only_a_predicted_growth() {
        let run = |guard: Option<bool>, incremental: bool| {
            let mut p = PulayMixer::new(0.5, 6).unwrap();
            p.set_incremental(incremental);
            if let Some(refit) = guard {
                p.set_coulomb_guard(Box::new(Weights4([1.0, 1.0, 1e6, 1e6])), refit);
            }
            let mut x = start();
            let mut steps = Vec::new();
            for _ in 0..2 {
                x = p.mix(&x, &sloshing(&x), 1.0).unwrap();
                steps.push((x.clone(), p.last_step().unwrap().clone()));
            }
            steps
        };
        for incremental in [false, true] {
            let plain = run(None, incremental);
            let trace = run(Some(false), incremental);
            let on = run(Some(true), incremental);
            for ((a, _), (b, t)) in plain.iter().zip(&trace) {
                let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
                assert_eq!(bits(a), bits(b));
                assert!(!t.coulomb_refit);
            }
            // Step 2: the L2 fit predicts an almost tripled v residual (the
            // predicted v residual is 4v - 8v gamma; the metric weights v by
            // 1e6, so the ratio is |1 - 2 gamma| up to the u part), trace
            // records it, on refits.
            let gamma = trace[1].1.coefficients[0];
            let ratio = trace[1].1.coulomb_predicted_ratio.unwrap();
            assert!((ratio - (1.0 - 2.0 * gamma)).abs() < 1e-3, "ratio {ratio}");
            assert!(ratio > 2.5);
            assert!(on[1].1.coulomb_refit && !on[0].1.coulomb_refit);
            // The refit coefficient removes the predicted v growth: the next
            // input has |v| no larger than the current one, the L2 step had 3x.
            let v = |x: &[f64]| (x[2] - x[3]).abs() / 2.0;
            let x1 = &plain[0].0;
            assert!(v(&on[1].0) <= 1.01 * v(x1), "{} vs {}", v(&on[1].0), v(x1));
            assert!(v(&plain[1].0) > 2.0 * v(x1));
            assert!((on[1].1.coefficients[0] - 0.5).abs() < 0.05);
        }
    }
}
