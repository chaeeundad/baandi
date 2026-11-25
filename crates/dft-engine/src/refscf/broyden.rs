//! Unnormalized difference-history Anderson update (Johnson 1988 / Pulay 1980).
//! F=y-x; dF=F_previous-F, dx=x_previous-x. Slot (k-2)%m receives the newest
//! pair before A and c are formed; k=1 only saves x,F. No clipping or charge
//! renormalization. The symmetric spectral factorization is the existing faer
//! adapter (MKL dsyevd on that backend), not a hand-written elimination.
//! A=Q Lambda Q^T, gamma=Q Lambda^-1 Q^T c. A numerical rank deficiency
//! (lambda_min <= p*eps*lambda_max) is a reported error, never a ridge or
//! pseudoinverse. cond_2=lambda_max/lambda_min is reported on every solve.
use super::memory::{Component as Mem, bytes, observe};
use super::{
    Error, invalid,
    metric::Metric,
    numerical,
    profile::{Preconditioning, Profile},
};
use dft_mpi::Collective;
use dft_numerics::real::{RealEighWorkspace, RealMatrix};
use dft_physics::data::Complex64 as C;
use std::f64::consts::PI;
const MAX_FORMS: usize = 25 * 26 / 2 + 25;
struct Solve {
    a: RealMatrix,
    q: RealMatrix,
    lambda: Vec<f64>,
}
/// Optional scratch for a read-only secant amplification model. Allocated once per
/// opted-in run; no density history is inserted by a prediction.
pub struct PredictionScratch {
    dx: Vec<C>,
    df: Vec<C>,
    forms: [f64; 3 * 25 + 1],
}
impl PredictionScratch {
    pub fn new(len: usize) -> Self {
        let scratch = Self {
            dx: vec![C::default(); len],
            df: vec![C::default(); len],
            forms: [0.; 3 * 25 + 1],
        };
        observe(Mem::Prediction, bytes(&scratch.dx) + bytes(&scratch.df));
        scratch
    }
}
#[derive(Clone, Copy, Debug, Default)]
pub struct Diagnostics {
    pub iteration: usize,
    pub occupancy: usize,
    /// Zero based slot overwritten; None at startup.
    pub inserted_slot: Option<usize>,
    pub gamma_l2: f64,
    pub gamma_max: f64,
    pub condition: Option<f64>,
}
pub struct Broyden {
    profile: Profile,
    iteration: usize,
    previous_x: Vec<C>,
    previous_f: Vec<C>,
    df: Vec<Vec<C>>,
    dx: Vec<Vec<C>>,
    smooth: Vec<bool>,
    screening: Vec<f64>,
    solves: Vec<Solve>,
    factor: RealEighWorkspace,
    forms: [f64; MAX_FORMS],
    gamma: [f64; 25],
    incremental: bool,
    gram: [[f64; 25]; 25],
    gram_metric: Option<u64>,
}
impl Broyden {
    /// ONEPATH-3 charge/magnetization history. TF must not freeze a free
    /// uniform moment: retain plain screening for that spin-only mode.
    pub(super) fn spin(
        profile: Profile,
        g2: &[f64],
        smooth: &[bool],
        volume: f64,
        electrons: f64,
    ) -> Result<Self, Error> {
        let mut mixer = Self::new(profile, &g2.repeat(2), &smooth.repeat(2), volume, electrons)?;
        for (i, &g) in g2.iter().enumerate() {
            if g == 0. {
                mixer.screening[g2.len() + i] = 1.;
            }
        }
        Ok(mixer)
    }

    pub fn new(
        profile: Profile,
        g2: &[f64],
        smooth: &[bool],
        volume: f64,
        electrons: f64,
    ) -> Result<Self, Error> {
        let profile = profile.validate()?;
        if g2.len() != smooth.len()
            || !volume.is_finite()
            || volume <= 0.0
            || !electrons.is_finite()
            || electrons <= 0.0
            || g2.iter().any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(invalid("Broyden geometry/shape"));
        }
        let rs = (3.0 * volume / (4.0 * PI * electrons)).cbrt();
        let a = (12.0 / PI).powf(2.0 / 3.0) / rs; // bohr^-2
        let screening = g2
            .iter()
            .map(|g| match profile.preconditioning {
                Preconditioning::Plain => 1.0,
                Preconditioning::ThomasFermi => g / (g + a),
            })
            .collect();
        let mut solves = Vec::with_capacity(profile.history);
        for p in 1..=profile.history {
            solves.push(Solve {
                a: RealMatrix::zeros(p, p).map_err(numerical)?,
                q: RealMatrix::zeros(p, p).map_err(numerical)?,
                lambda: vec![0.; p],
            });
        }
        let mixer = Self {
            profile,
            iteration: 0,
            previous_x: vec![C::default(); g2.len()],
            previous_f: vec![C::default(); g2.len()],
            df: vec![vec![C::default(); g2.len()]; profile.history],
            dx: vec![vec![C::default(); g2.len()]; profile.history],
            smooth: smooth.to_vec(),
            screening,
            solves,
            factor: RealEighWorkspace::new(profile.history).map_err(numerical)?,
            forms: [0.; MAX_FORMS],
            gamma: [0.; 25],
            incremental: dft_numerics::switches::on(dft_numerics::switches::Id::REFSCF_PERMAP),
            gram: [[0.; 25]; 25],
            gram_metric: None,
        };
        observe(
            Mem::History,
            bytes(&mixer.previous_x)
                + bytes(&mixer.previous_f)
                + mixer.df.iter().chain(&mixer.dx).map(bytes).sum::<usize>(),
        );
        observe(
            Mem::Mixer,
            bytes(&mixer.smooth)
                + bytes(&mixer.screening)
                + mixer
                    .solves
                    .iter()
                    .map(|s| {
                        std::mem::size_of_val(s.a.as_slice())
                            + std::mem::size_of_val(s.q.as_slice())
                            + bytes(&s.lambda)
                    })
                    .sum::<usize>(),
        );
        Ok(mixer)
    }
    pub fn occupancy(&self) -> usize {
        self.iteration.saturating_sub(1).min(self.profile.history)
    }
    pub fn gamma(&self) -> &[f64] {
        &self.gamma[..self.occupancy()]
    }
    /// Estimate the next unmixed residual energy from the largest observed
    /// amplification ||Delta F(x)||_H^2 / ||Delta x||_H^2 in the mixer history.
    /// Include the current secant virtually, with ordinary ring eviction.
    /// This is exact for a scalar linear map, not a nonlinear upper bound.
    /// Mixer history, factor scratch, gamma and iteration are read-only.
    pub fn predict_output_error<K: Collective>(
        &self,
        comm: &K,
        metric: &mut Metric,
        x: &[C],
        f: &[C],
        scratch: &mut PredictionScratch,
    ) -> Result<Option<f64>, Error> {
        if x.len() != self.previous_x.len() || f.len() != x.len() || scratch.dx.len() != x.len() {
            return Err(invalid("prediction vector shape"));
        }
        if self.iteration == 0 {
            return Ok(None);
        }
        let p = self.iteration.min(self.profile.history);
        let slot = (self.iteration - 1) % self.profile.history;
        for j in 0..x.len() {
            scratch.dx[j] = self.previous_x[j] - x[j];
            scratch.df[j] = self.previous_f[j] - f[j];
        }
        let dx = |i: usize| {
            if i == slot {
                &scratch.dx[..]
            } else {
                &self.dx[i][..]
            }
        };
        let df = |i: usize| {
            if i == slot {
                &scratch.df[..]
            } else {
                &self.df[i][..]
            }
        };
        // Observed amplification of the unmixed map F along each stored
        // density secant: ||S + D||_H^2 / ||S||_H^2. Unlike the mixed dr2
        // ratio, this measures F itself, without a fitted multiplier.
        let mut pairs: [(&[C], &[C]); 76] = [(&[], &[]); 76];
        for i in 0..p {
            pairs[3 * i] = (dx(i), dx(i));
            pairs[3 * i + 1] = (dx(i), df(i));
            pairs[3 * i + 2] = (df(i), df(i));
        }
        pairs[3 * p] = (f, f);
        let n = 3 * p + 1;
        if !metric.prediction_batch(comm, &pairs[..n], &mut scratch.forms[..n])? {
            return Ok(None);
        }
        let mut amplification = 0_f64;
        for i in 0..p {
            let [xx, xf, ff] = [
                scratch.forms[3 * i],
                scratch.forms[3 * i + 1],
                scratch.forms[3 * i + 2],
            ];
            let yy = xx + 2. * xf + ff;
            let ratio = yy / xx;
            if xx <= 0. || yy < 0. || !ratio.is_finite() {
                return Ok(None);
            }
            amplification = amplification.max(ratio);
        }
        let predicted = amplification * (scratch.forms[3 * p] / 2.);
        Ok(predicted.is_finite().then_some(predicted))
    }
    /// Call only for an accepted mixing step: neither convergence nor the
    /// first-map redo inserts history. An error terminates this mixer instance.
    pub fn update<K: Collective>(
        &mut self,
        comm: &K,
        metric: &mut Metric,
        x: &[C],
        f: &[C],
        next: &mut [C],
    ) -> Result<Diagnostics, Error> {
        if x.len() != self.previous_x.len() || f.len() != x.len() || next.len() != x.len() {
            return Err(invalid("Broyden vector shape"));
        }
        self.iteration += 1;
        let p = self.occupancy();
        let slot = (self.iteration > 1).then(|| (self.iteration - 2) % self.profile.history);
        if let Some(slot) = slot {
            for j in 0..x.len() {
                self.df[slot][j] = self.previous_f[j] - f[j];
                self.dx[slot][j] = self.previous_x[j] - x[j];
            }
        }
        self.previous_x.copy_from_slice(x);
        self.previous_f.copy_from_slice(f);
        self.gamma.fill(0.0);
        let mut condition = None;
        if p > 0 {
            let mut pairs: [(&[C], &[C]); MAX_FORMS] = [(&[], &[]); MAX_FORMS];
            let mut n = 0;
            let reuse = self.incremental && self.gram_metric == Some(metric.identity());
            for i in 0..p {
                for j in 0..=i {
                    if !reuse || slot == Some(i) || slot == Some(j) {
                        pairs[n] = (&self.df[i], &self.df[j]);
                        n += 1;
                    }
                }
            }
            for i in 0..p {
                pairs[n] = (&self.df[i], f);
                n += 1;
            }
            metric.batch(comm, &pairs[..n], &mut self.forms[..n])?;
            let solve = &mut self.solves[p - 1];
            let mut index = 0;
            for i in 0..p {
                for j in 0..=i {
                    if !reuse || slot == Some(i) || slot == Some(j) {
                        self.gram[i][j] = self.forms[index];
                        index += 1;
                    }
                    solve.a[(i, j)] = self.gram[i][j];
                    solve.a[(j, i)] = self.gram[i][j];
                }
            }
            self.gram_metric = Some(metric.identity());
            self.factor
                .solve_into(&solve.a, &mut solve.lambda, &mut solve.q)
                .map_err(numerical)?;
            let (min, max) = (solve.lambda[0], solve.lambda[p - 1]);
            let cond = if min > 0.0 { max / min } else { f64::INFINITY };
            if min <= p as f64 * f64::EPSILON * max || !cond.is_finite() {
                return Err(Error::Singular {
                    iteration: self.iteration,
                    min,
                    max,
                    condition: cond,
                });
            }
            condition = Some(cond);
            for j in 0..p {
                let projection: f64 = (0..p)
                    .map(|i| solve.q[(i, j)] * self.forms[index + i])
                    .sum();
                for i in 0..p {
                    self.gamma[i] += solve.q[(i, j)] * (projection / solve.lambda[j]);
                }
            }
        }
        for j in 0..x.len() {
            let (mut xx, mut ff) = (x[j], f[j]);
            if self.smooth[j] {
                for i in 0..p {
                    xx -= self.gamma[i] * self.dx[i][j];
                    ff -= self.gamma[i] * self.df[i][j];
                }
                ff *= self.screening[j]; // TF acts on F_tilde only, never dx.
            }
            next[j] = xx + self.profile.beta * ff;
        }
        if next.iter().any(|v| !v.re.is_finite() || !v.im.is_finite()) {
            return Err(numerical("nonfinite Broyden update"));
        }
        Ok(Diagnostics {
            iteration: self.iteration,
            occupancy: p,
            inserted_slot: slot,
            gamma_l2: self.gamma[..p].iter().map(|v| v * v).sum::<f64>().sqrt(),
            gamma_max: self.gamma[..p].iter().map(|v| v.abs()).fold(0., f64::max),
            condition,
        })
    }
}
