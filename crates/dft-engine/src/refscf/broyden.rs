//! Unnormalized difference-history Anderson update (Johnson 1988 / Pulay 1980).
//! F=y-x; dF=F_previous-F, dx=x_previous-x. Slot (k-2)%m receives the newest
//! pair before A and c are formed; k=1 only saves x,F. No clipping or charge
//! renormalization. The symmetric spectral factorization is the existing faer
//! adapter (MKL dsyevd on that backend), not a hand-written elimination.
//! A=Q Lambda Q^T, gamma=Q Lambda^-1 Q^T c. A numerical rank deficiency
//! (lambda_min <= p*eps*lambda_max) is a reported error by default, never a ridge or
//! pseudoinverse. Explicit conditioning opt-in follows proposed ADR 0045.
//! cond_2=lambda_max/lambda_min is reported on every solve.
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
/// Candidate-only metadata, kept separate from the legacy scientific record.
#[derive(Clone, Copy, Debug, Default)]
pub struct ConditioningDiagnostics {
    /// Oldest-first removals, with the pre-removal condition of each solve.
    pub dropped_pairs: usize,
    pub dropped_conditions: [f64; 25],
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
    conditioning: bool,
    conditioning_diagnostics: ConditioningDiagnostics,
    retained: usize,
    /// Scratch: the history-reduced residual F~ of the last update (before
    /// screening), the least-squares remainder of F on the retained secants.
    reduced: Vec<C>,
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
            conditioning: false,
            conditioning_diagnostics: ConditioningDiagnostics::default(),
            retained: 0,
            reduced: vec![C::default(); g2.len()],
        };
        observe(
            Mem::History,
            bytes(&mixer.previous_x)
                + bytes(&mixer.previous_f)
                + bytes(&mixer.reduced)
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
    /// Opt-in before the first update. Ordinary drivers retain rejection.
    pub fn condition_history(&mut self) {
        assert_eq!(self.iteration, 0, "configure history before mixing");
        self.conditioning = true;
    }
    pub fn conditioning_diagnostics(&self) -> ConditioningDiagnostics {
        self.conditioning_diagnostics
    }
    pub fn occupancy(&self) -> usize {
        if self.conditioning {
            self.retained
        } else {
            self.iteration.saturating_sub(1).min(self.profile.history)
        }
    }
    /// Read-only geometry of the most recently assembled history Gram matrix.
    /// Pair squared norms (Ry) and largest absolute pair cosine diagnose scale
    /// collapse and duplicate secants without changing the solve.
    pub fn history_geometry(&self) -> ([f64; 25], f64) {
        let p = self.occupancy();
        let mut norms = [0.; 25];
        let mut cosine = 0_f64;
        for (i, norm) in norms.iter_mut().enumerate().take(p) {
            *norm = self.gram[i][i];
            for j in 0..i {
                let denominator = (self.gram[i][i] * self.gram[j][j]).sqrt();
                if denominator > 0. {
                    cosine = cosine.max(self.gram[i][j].abs() / denominator);
                }
            }
        }
        (norms, cosine)
    }
    pub fn gamma(&self) -> &[f64] {
        &self.gamma[..self.occupancy()]
    }
    /// History-reduced residual F~ = F - sum_i gamma_i dF_i of the last
    /// accepted update, before screening (equal to F without history). Its
    /// metric norm is the part of the residual the secant model leaves
    /// unexplained: GM1 uses it to forecast the next map's residual.
    pub fn reduced_residual(&self) -> &[C] {
        &self.reduced
    }

    /// GV1 (`VERIFY_RESPONSE_PREDICTOR`): estimate the fixed-point residual
    /// E_H[F(y)] at the output y = x + F of the current map, before any
    /// verification map exists. Linearly F(y) = F + J F = K F with K = 1 + J
    /// the output response (dy = K dx). The secant pairs give K exactly on
    /// their span: dy_i = dx_i + dF_i. F is projected onto span(dx_i) in the
    /// Coulomb metric (least squares, F ~ sum c_i dx_i), and the prediction is
    /// ||sum c_i dy_i||^2 (Ry; the caller halves it). The pairs are the
    /// retained history plus the current pair (previous x, F -> x, F) with the
    /// ring eviction of the next update. The part of F outside the span is
    /// not propagated (no assumed response). An estimator, never a bound; it
    /// reads the history only. Returns None without history, for a singular
    /// projection or a nonfinite value. Cold path (stop maps only): it
    /// allocates its scratch.
    pub fn predict_output_residual<K: Collective>(
        &self,
        comm: &K,
        metric: &mut Metric,
        x: &[C],
        f: &[C],
    ) -> Result<Option<f64>, Error> {
        if x.len() != self.previous_x.len() || f.len() != x.len() {
            return Err(invalid("prediction vector shape"));
        }
        if self.iteration == 0 {
            return Ok(None);
        }
        let current_dx: Vec<C> = self.previous_x.iter().zip(x).map(|(a, b)| a - b).collect();
        let current_df: Vec<C> = self.previous_f.iter().zip(f).map(|(a, b)| a - b).collect();
        // Retained pairs plus the current one; the ring (or the conditioned
        // list) evicts its oldest pair when full, as the next update would.
        let p = self.occupancy();
        let evicted = if p == self.profile.history {
            if self.conditioning {
                Some(0)
            } else {
                Some((self.iteration - 1) % self.profile.history)
            }
        } else {
            None
        };
        let mut pairs: Vec<(&[C], &[C])> = (0..p)
            .filter(|i| Some(*i) != evicted)
            .map(|i| (&self.dx[i][..], &self.df[i][..]))
            .collect();
        pairs.push((&current_dx, &current_df));
        let m = pairs.len();
        let mut forms_in: Vec<(&[C], &[C])> = Vec::with_capacity(m * (m + 1) / 2 + m);
        for i in 0..m {
            for j in 0..=i {
                forms_in.push((pairs[i].0, pairs[j].0));
            }
        }
        for pair in &pairs {
            forms_in.push((pair.0, f));
        }
        let mut forms = vec![0.; forms_in.len()];
        if !metric.prediction_batch(comm, &forms_in, &mut forms)? {
            return Ok(None);
        }
        let mut a = RealMatrix::zeros(m, m).map_err(numerical)?;
        let mut index = 0;
        for i in 0..m {
            for j in 0..=i {
                a[(i, j)] = forms[index];
                a[(j, i)] = forms[index];
                index += 1;
            }
        }
        let b = &forms[index..];
        let mut lambda = vec![0.; m];
        let mut q = RealMatrix::zeros(m, m).map_err(numerical)?;
        RealEighWorkspace::new(m)
            .map_err(numerical)?
            .solve_into(&a, &mut lambda, &mut q)
            .map_err(numerical)?;
        let max = lambda[m - 1];
        if !(max.is_finite() && max > 0.) {
            return Ok(None);
        }
        // Least squares on the numerically independent part of the span
        // (directions below the mixer's own rank floor are dropped).
        let mut c = vec![0.; m];
        for (j, &l) in lambda.iter().enumerate() {
            if l <= m as f64 * f64::EPSILON * max {
                continue;
            }
            let projection: f64 = (0..m).map(|i| q[(i, j)] * b[i]).sum();
            for (i, c) in c.iter_mut().enumerate() {
                *c += q[(i, j)] * (projection / l);
            }
        }
        let mut g = vec![C::default(); f.len()];
        for (pair, c) in pairs.iter().zip(&c) {
            for ((g, dx), df) in g.iter_mut().zip(pair.0).zip(pair.1) {
                *g += *c * (dx + df);
            }
        }
        let mut predicted = [0.];
        if !metric.prediction_batch(comm, &[(&g, &g)], &mut predicted)? {
            return Ok(None);
        }
        Ok(Some(predicted[0]))
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
        self.conditioning_diagnostics = ConditioningDiagnostics::default();
        self.iteration += 1;
        let slot = if self.conditioning && self.iteration > 1 {
            if self.retained == self.profile.history {
                self.df.rotate_left(1);
                self.dx.rotate_left(1);
                self.retained -= 1;
            }
            let slot = self.retained;
            self.retained += 1;
            Some(slot)
        } else {
            (self.iteration > 1).then(|| (self.iteration - 2) % self.profile.history)
        };
        let mut p = self.occupancy();
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
        let mut dropped_pairs = 0;
        let mut dropped_conditions = [0.; 25];
        while p > 0 {
            let mut pairs: [(&[C], &[C]); MAX_FORMS] = [(&[], &[]); MAX_FORMS];
            let mut n = 0;
            let reuse = !self.conditioning
                && self.incremental
                && self.gram_metric == Some(metric.identity());
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
            if self.conditioning
                && (max == 0. || min <= (p as f64 * f64::EPSILON).sqrt() * max || !cond.is_finite())
            {
                if !min.is_finite() || !max.is_finite() {
                    return Err(numerical("nonfinite history spectrum"));
                }
                dropped_conditions[dropped_pairs] = cond;
                dropped_pairs += 1;
                // Preserve every event even if a later rebuild/backend fails.
                self.conditioning_diagnostics = ConditioningDiagnostics {
                    dropped_pairs,
                    dropped_conditions,
                };
                self.df.rotate_left(1);
                self.dx.rotate_left(1);
                self.retained -= 1;
                p = self.retained;
                self.gram_metric = None;
                continue;
            }
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
            break;
        }
        for j in 0..x.len() {
            let (mut xx, mut ff) = (x[j], f[j]);
            if self.smooth[j] {
                for i in 0..p {
                    xx -= self.gamma[i] * self.dx[i][j];
                    ff -= self.gamma[i] * self.df[i][j];
                }
                self.reduced[j] = ff;
                ff *= self.screening[j]; // TF acts on F_tilde only, never dx.
            } else {
                self.reduced[j] = ff;
            }
            next[j] = xx + self.profile.beta * ff;
        }
        if next.iter().any(|v| !v.re.is_finite() || !v.im.is_finite()) {
            return Err(numerical("nonfinite Broyden update"));
        }
        self.conditioning_diagnostics = ConditioningDiagnostics {
            dropped_pairs,
            dropped_conditions,
        };
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

impl Broyden {
    /// Scratch eigensolver workspaces are rebuilt; history and the incremental
    /// Gram cache are authoritative and retain their exact bits and slot order.
    pub(super) fn checkpoint(
        &mut self,
        io: &mut super::checkpoint::Exchange<'_>,
        metric: &Metric,
    ) -> Result<(), crate::model::ScfError> {
        use super::checkpoint::invalid;
        io.same(&format!("{:?}", self.profile))?;
        io.same(&self.smooth)?;
        io.same(&self.screening)?;
        io.same(&self.incremental)?;
        io.same(&self.conditioning)?;
        io.value(&mut self.iteration)?;
        io.reciprocal(&mut self.previous_x)?;
        io.reciprocal(&mut self.previous_f)?;
        let lengths: Vec<_> = self.df.iter().chain(&self.dx).map(Vec::len).collect();
        io.value(&mut self.df)?;
        io.value(&mut self.dx)?;
        if self.df.len() != self.profile.history
            || self.dx.len() != self.profile.history
            || !self.df.iter().chain(&self.dx).map(Vec::len).eq(lengths)
            || self
                .df
                .iter()
                .chain(&self.dx)
                .flatten()
                .any(|v| !v.re.is_finite() || !v.im.is_finite())
        {
            return Err(invalid("Broyden history shape/values"));
        }
        io.value(&mut self.gram)?;
        let mut valid = self.gram_metric.is_some();
        io.value(&mut valid)?;
        if io.loading() {
            self.gram_metric = valid.then(|| metric.identity());
        }
        io.value(&mut self.retained)?;
        io.value(&mut self.conditioning_diagnostics)?;
        if self.retained > self.profile.history
            || self.gram.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(invalid("Broyden Gram/retained history"));
        }
        Ok(())
    }
}
