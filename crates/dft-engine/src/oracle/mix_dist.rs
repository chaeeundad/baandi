//! Distributed Pulay mixing of the M2 path (ticket M2-GRID, 2026-11-18;
//! frozen default with numerical eligibility fallback).
//!
//! The replicated `PulayMixer` keeps every input density and residual of its
//! window on every rank and forms its differences, Gram entries and the update
//! over the full grid: grid work that does not shrink with the rank count. This
//! mixer keeps only the rank's block of the grid. The grid is cut into
//! fixed chunks of `CHUNK` points (independent of the rank count); every rank
//! owns a contiguous range of whole chunks. Per `mix` call:
//!
//! * the newest difference pair (dx, dr) is formed on the block;
//! * its Gram row and the right-hand side dr_i . r are partial sums per chunk,
//!   gathered with ONE allgatherv and added in ascending chunk order (a
//!   canonical order: the same bits for every rank count), the earlier Gram
//!   entries are kept;
//! * the small eigenproblem and the coefficients are replicated (the same code
//!   as the replicated mixer: spectral cut 1e-12 of the largest eigenvalue);
//! * the update and the two uniform charge corrections run on the block (one
//!   chunk-sum allgatherv per correction) and ONE allgatherv assembles the new
//!   input density on every rank.
//!
//! Same window rule as `PulayMixer` (at most `history` differences; seeded
//! differences first, dropped oldest first). Rounding class against the
//! replicated mixer (the dot products and charge sums are added in another
//! order); a fixed rank count is bit-reproducible and different rank counts
//! give the same bits. Exact checkpoints reject this mixer (opaque state,
//! `MixerHistory::checkpoint` default): the switch is off by default.
//! Every collective is entered by every rank (replicated control flow).
use dft_mpi::Collective;
use dft_numerics::linalg::{ComplexMatrix, hermitian_eigh};
use dft_numerics::switches::{GlueProfile, RunCounter as Gc};
use dft_solvers::SolverError;
use dft_solvers::kerker::MixerHistory;
use dft_solvers::mixing::{DensityMixer, PulayStepDiagnostics};
use std::collections::VecDeque;
use std::sync::atomic::{AtomicU64, Ordering};

/// Points per canonical chunk (partial sums are formed per chunk).
pub const CHUNK: usize = 4096;

thread_local! {
    static OVERRIDE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}
/// Test override of distributed mixing for sessions initialized on the calling thread.
pub fn set_override(on: Option<bool>) {
    OVERRIDE.with(|o| o.set(on));
}
pub fn requested() -> bool {
    OVERRIDE.with(|o| o.get()).unwrap_or(true)
}

static COUNTS: [AtomicU64; 3] = [const { AtomicU64::new(0) }; 3];
/// (distributed mix calls with a Pulay update, Gram entries formed, mixers
/// built) since process start, summed over the threads of the process.
pub fn counters() -> [u64; 3] {
    std::array::from_fn(|i| COUNTS[i].load(Ordering::Relaxed))
}

struct Diff {
    dx: Vec<f64>,
    dr: Vec<f64>,
}

pub(crate) struct DistributedPulay<K: Collective> {
    last_step: Option<PulayStepDiagnostics>,
    glue_profile: Option<GlueProfile>,
    comm: K,
    beta: f64,
    capacity: usize,
    size: usize,
    lo: usize,
    hi: usize,
    /// Points per canonical chunk: `CHUNK` (row-major chunks) or one slab
    /// plane (plane mode, M2-OWNED).
    chunk: usize,
    /// Plane mode (`BAANDI_M2_LOCAL_DENSITY`): the block is this rank's slab
    /// planes, handed to `mix_block` directly; `mix` (full arrays) is refused.
    #[cfg(test)]
    planes: bool,
    rank: usize,
    /// Chunk counts per rank (allgatherv layout).
    chunk_counts: Vec<usize>,
    mean: Option<f64>,
    /// Latest absolute pair (input, residual) on the block.
    last: Option<(Vec<f64>, Vec<f64>)>,
    /// Absolute pairs in the window (`PulayMixer::history_len`).
    pairs: usize,
    seeds: usize,
    /// Seeds first, then window differences, oldest first.
    diffs: VecDeque<Diff>,
    /// Symmetric Gram entries dr_i . dr_j aligned with `diffs`.
    gram: VecDeque<VecDeque<f64>>,
    spare: Vec<Diff>,
    partial: Vec<f64>,
    gathered: Vec<f64>,
}

/// (guard checks, refits, distributed transforms of the guard) since process
/// start, summed over the threads of the process (simulated ranks add up).
static GUARD_COUNTS: [AtomicU64; 3] = [const { AtomicU64::new(0) }; 3];
pub fn guard_counters() -> [u64; 3] {
    std::array::from_fn(|i| GUARD_COUNTS[i].load(Ordering::Relaxed))
}

/// Replicated spectral least squares of the Pulay step (the cut and the
/// formulas of `PulayMixer`): (eigenvalues, largest, gamma).
fn spectral_solve(
    gram_entries: impl Fn(usize, usize) -> f64,
    rhs: &[f64],
) -> Result<(Vec<f64>, f64, Vec<f64>), SolverError> {
    let m = rhs.len();
    let mut gram = ComplexMatrix::zeros(m, m)?;
    for i in 0..m {
        for j in 0..m {
            gram[(i, j)] = gram_entries(i, j).into();
        }
    }
    let (values, vectors) = hermitian_eigh(&gram)?;
    let largest = values.iter().copied().fold(0.0, f64::max);
    let mut gamma = vec![0.0; m];
    for (k, &value) in values.iter().enumerate() {
        if value > 1e-12 * largest {
            let projection: f64 = (0..m).map(|i| vectors[(i, k)].re * rhs[i]).sum::<f64>() / value;
            for i in 0..m {
                gamma[i] += vectors[(i, k)].re * projection;
            }
        }
    }
    Ok((values, largest, gamma))
}

/// Least-squares summary for `PulayStepDiagnostics`.
struct Fit {
    values: Vec<f64>,
    largest: f64,
    gamma: Vec<f64>,
    gamma_rhs: f64,
    /// R^T R of the plain fit, R^T C R after a Coulomb refit (None: R^T R,
    /// known only after the charge reduction).
    rr_fit: Option<f64>,
}

impl<K: Collective> DistributedPulay<K> {
    pub(crate) fn new(
        comm: K,
        beta: f64,
        history: usize,
        size: usize,
    ) -> Result<Self, SolverError> {
        if !beta.is_finite() || beta <= 0.0 || history == 0 || size == 0 {
            return Err(SolverError::InvalidParameter(
                "positive beta and history required",
            ));
        }
        let (rank, ranks) = (comm.rank(), comm.size());
        let chunks = size.div_ceil(CHUNK);
        let start = |r: usize| r * (chunks / ranks) + r.min(chunks % ranks);
        let chunk_counts = (0..ranks).map(|r| start(r + 1) - start(r)).collect();
        let lo = (start(rank) * CHUNK).min(size);
        let hi = (start(rank + 1) * CHUNK).min(size);
        COUNTS[2].fetch_add(1, Ordering::Relaxed);
        Ok(Self {
            glue_profile: GlueProfile::current(),
            comm,
            beta,
            capacity: history,
            size,
            lo,
            hi,
            chunk: CHUNK,
            #[cfg(test)]
            planes: false,
            rank,
            chunk_counts,
            mean: None,
            last: None,
            pairs: 0,
            seeds: 0,
            diffs: VecDeque::new(),
            gram: VecDeque::new(),
            spare: Vec::new(),
            partial: Vec::new(),
            gathered: Vec::new(),
            last_step: None,
        })
    }

    /// Plane mode (M2-OWNED): the window lives on this rank's slab planes of
    /// `plane_len` points each (`plane_counts` planes per rank, ascending plane
    /// order = rank order); every canonical partial sum is per plane, so the
    /// results do not depend on the rank count.
    #[cfg(test)]
    pub(crate) fn new_planes(
        comm: K,
        beta: f64,
        history: usize,
        size: usize,
        plane_len: usize,
        plane_counts: &[usize],
    ) -> Result<Self, SolverError> {
        let mut m = Self::new(comm, beta, history, size)?;
        if plane_len == 0 || plane_counts.len() != m.comm.size() {
            return Err(SolverError::InvalidParameter("plane layout"));
        }
        m.chunk = plane_len;
        m.planes = true;
        m.chunk_counts = plane_counts.to_vec();
        m.lo = 0;
        m.hi = plane_counts[m.rank] * plane_len;
        Ok(m)
    }

    fn local_chunks(&self) -> usize {
        self.chunk_counts[self.rank]
    }

    /// Length of this rank's block.
    pub(crate) fn block_len(&self) -> usize {
        self.hi - self.lo
    }

    /// Canonical sums of `k` quantities: `fill(c, out)` writes the k partial
    /// sums of local chunk c; the totals come back in chunk order.
    fn chunk_sums(
        &mut self,
        k: usize,
        mut fill: impl FnMut(usize, &mut [f64]),
    ) -> Result<Vec<f64>, SolverError> {
        let glue_span = self
            .glue_profile
            .map(|p| p.span(Gc::GlueDistLocalCalls, Gc::GlueDistLocalNanos));
        let local = self.local_chunks();
        self.partial.clear();
        self.partial.resize(local * k, 0.0);
        for c in 0..local {
            fill(c, &mut self.partial[c * k..(c + 1) * k]);
        }
        let counts: Vec<usize> = self.chunk_counts.iter().map(|c| c * k).collect();
        self.gathered.resize(counts.iter().sum(), 0.0);
        drop(glue_span);
        if let Some(p) = self.glue_profile {
            p.count(Gc::GlueDistSendValues, self.partial.len());
            p.count(Gc::GlueDistRecvValues, self.gathered.len());
        }
        let glue_span = self
            .glue_profile
            .map(|p| p.span(Gc::GlueDistGatherCalls, Gc::GlueDistGatherNanos));
        self.comm
            .allgatherv_f64(&self.partial, &counts, &mut self.gathered)
            .map_err(|e| SolverError::Application(e.to_string()))?;
        drop(glue_span);
        let _glue = self
            .glue_profile
            .map(|p| p.span(Gc::GlueDistReduceCalls, Gc::GlueDistReduceNanos));
        // Rank blocks are ascending chunk ranges: the buffer is in chunk order.
        let mut out = vec![0.0; k];
        for chunk in self.gathered.chunks_exact(k) {
            for (o, v) in out.iter_mut().zip(chunk) {
                *o += v;
            }
        }
        Ok(out)
    }

    fn drop_diff(&mut self, at: usize) {
        if let Some(d) = self.diffs.remove(at) {
            self.spare.push(d);
        }
        self.gram.remove(at);
        for row in &mut self.gram {
            row.remove(at);
        }
    }

    fn push_diff(&mut self, dx: impl Iterator<Item = f64>, dr: impl Iterator<Item = f64>) {
        let mut d = self.spare.pop().unwrap_or(Diff {
            dx: Vec::new(),
            dr: Vec::new(),
        });
        d.dx.clear();
        d.dr.clear();
        d.dx.extend(dx);
        d.dr.extend(dr);
        self.diffs.push_back(d);
        // Placeholder row/column; filled by the next Gram reduction.
        for row in &mut self.gram {
            row.push_back(f64::NAN);
        }
        self.gram
            .push_back(VecDeque::from(vec![f64::NAN; self.diffs.len()]));
    }

    /// Two uniform charge corrections of the block (`mixing::charge`), each
    /// with the canonical global sum. The first reduction also carries the
    /// diagnostics sums |next - x|^2 and |r|^2 (separate accumulators: the
    /// charge sum keeps its bits), returned as (step^2, residual^2).
    fn charge(
        &mut self,
        next: &mut [f64],
        mean: f64,
        x: &[f64],
        r: &[f64],
    ) -> Result<(f64, f64), SolverError> {
        let _glue = self
            .glue_profile
            .map(|p| p.span(Gc::GlueDistChargeCalls, Gc::GlueDistChargeNanos));
        let chunk = self.chunk;
        let mut norms = (0.0, 0.0);
        for pass in 0..2 {
            let block: &[f64] = next;
            let k = if pass == 0 { 3 } else { 1 };
            let sums = self.chunk_sums(k, |c, out| {
                let a = c * chunk;
                let b = ((c + 1) * chunk).min(block.len());
                out[0] = block[a..b].iter().sum::<f64>();
                if pass == 0 {
                    out[1] = block[a..b]
                        .iter()
                        .zip(&x[a..b])
                        .map(|(n, x)| (n - x) * (n - x))
                        .sum::<f64>();
                    out[2] = r[a..b].iter().map(|v| v * v).sum::<f64>();
                }
            })?;
            if pass == 0 {
                norms = (sums[1], sums[2]);
            }
            let shift = mean - sums[0] / self.size as f64;
            for v in next.iter_mut() {
                *v += shift;
            }
        }
        Ok(norms)
    }
}

impl<K: Collective> DensityMixer for DistributedPulay<K> {
    fn reset(&mut self) {
        while !self.diffs.is_empty() {
            self.drop_diff(0);
        }
        self.last = None;
        self.pairs = 0;
        self.seeds = 0;
        self.mean = None;
        self.last_step = None;
    }

    fn mix(&mut self, input: &[f64], output: &[f64], mean: f64) -> Result<Vec<f64>, SolverError> {
        // Replicated validation (the full densities are replicated).
        if input.len() != self.size
            || output.len() != self.size
            || !mean.is_finite()
            || mean <= 0.0
            || input.iter().chain(output).any(|v| !v.is_finite())
        {
            return Err(SolverError::InvalidParameter(
                "finite equal nonempty densities, positive mean and beta required",
            ));
        }
        #[cfg(test)]
        if self.planes {
            return Err(SolverError::InvalidParameter(
                "plane-mode mixer takes owned blocks (mix_block)",
            ));
        }
        let (lo, hi) = (self.lo, self.hi);
        let mut next = self.mix_block_with(&input[lo..hi], &output[lo..hi], mean, None)?;
        if next.iter().any(|v| !v.is_finite()) {
            // Every rank must take the same branch: the flag rides with the gather.
            next.iter_mut().for_each(|v| *v = f64::NAN);
        }
        let counts: Vec<usize> = (0..self.chunk_counts.len())
            .map(|q| {
                let start: usize = self.chunk_counts[..q].iter().sum();
                let a = (start * CHUNK).min(self.size);
                let b = ((start + self.chunk_counts[q]) * CHUNK).min(self.size);
                b - a
            })
            .collect();
        let mut full = vec![0.0; self.size];
        let clock = dft_numerics::stage_clock::start();
        self.comm
            .allgatherv_f64(&next, &counts, &mut full)
            .map_err(|e| SolverError::Application(e.to_string()))?;
        {
            use crate::grid_dist as g;
            g::owned_add(g::O_MIXER_GATHERS, 1);
            g::owned_add(g::O_MIXER_GATHER_VALUES, self.size as u64);
            g::owned_time(g::T_MIXER_GATHER, clock);
        }
        if full.iter().any(|v| !v.is_finite()) {
            return Err(dft_numerics::linalg::LinalgError::NonFinite.into());
        }
        Ok(full)
    }
}

impl<K: Collective> DistributedPulay<K> {
    /// The Pulay update on this rank's block (`input` and `output` are the
    /// block's values, `block_len()` each): the next input density on the
    /// block, without a gather (M2-OWNED plane mode; also the body of `mix`).
    /// Collective (Gram row and charge sums). The inputs are not validated
    /// locally (an early error on one rank would desynchronize the
    /// collectives): a nonfinite value reaches every rank through the charge
    /// sums, and every rank's result is then nonfinite.
    #[cfg(test)]
    pub(crate) fn mix_block(
        &mut self,
        x: &[f64],
        output: &[f64],
        mean: f64,
    ) -> Result<Vec<f64>, SolverError> {
        self.mix_block_with(x, output, mean, None)
    }

    /// `mix_block` with the guard's plane field of the residual when the
    /// caller has it (chunk mode: from the replicated full densities, so no
    /// gather is needed for R).
    fn mix_block_with(
        &mut self,
        x: &[f64],
        output: &[f64],
        mean: f64,
        _r_planes: Option<Vec<f64>>,
    ) -> Result<Vec<f64>, SolverError> {
        let _glue = self
            .glue_profile
            .map(|p| p.span(Gc::GlueDistMixCalls, Gc::GlueDistMixNanos));
        let glue_span = self
            .glue_profile
            .map(|p| p.span(Gc::GlueDistPrepareCalls, Gc::GlueDistPrepareNanos));
        if x.len() != self.block_len()
            || output.len() != x.len()
            || !mean.is_finite()
            || mean <= 0.0
        {
            return Err(SolverError::InvalidParameter(
                "mix_block: block length or mean",
            ));
        }
        if self.mean != Some(mean) {
            self.reset();
        }
        self.mean = Some(mean);
        let (beta, chunk) = (self.beta, self.chunk);
        let r: Vec<f64> = output.iter().zip(x).map(|(b, a)| b - a).collect();
        if let Some((x0, r0)) = self.last.take() {
            self.push_diff(
                x.iter().zip(&x0).map(|(a, b)| a - b),
                r.iter().zip(&r0).map(|(a, b)| a - b),
            );
        }
        self.last = Some((x.to_vec(), r.clone()));
        self.pairs += 1;
        if self.pairs > self.capacity + 1 {
            self.pairs -= 1;
            self.drop_diff(self.seeds);
        }
        while self.seeds > 0 && self.seeds + self.pairs - 1 > self.capacity {
            self.seeds -= 1;
            self.drop_diff(0);
        }
        let mut next: Vec<f64> = x.iter().zip(&r).map(|(a, r)| a + beta * r).collect();
        let m = self.diffs.len();
        drop(glue_span);
        if let Some(p) = self.glue_profile {
            p.count(Gc::GlueDistPoints, x.len());
            p.count(Gc::GlueDistHistory, m);
            p.count(Gc::GlueDistUpdateElements, m * x.len());
        }
        let mut fit: Option<Fit> = None;
        if m > 0 {
            // Missing Gram entries (new rows) and the right-hand side in one
            // canonical reduction.
            let missing: Vec<(usize, usize)> = (0..m)
                .flat_map(|i| (i..m).map(move |j| (i, j)))
                .filter(|&(i, j)| self.gram[i][j].is_nan())
                .collect();
            let k = missing.len() + m;
            if let Some(p) = self.glue_profile {
                p.count(Gc::GlueDistDotElements, k * r.len());
            }
            let diffs = std::mem::take(&mut self.diffs);
            let sums = self.chunk_sums(k, |c, out| {
                let a = c * chunk;
                let b = ((c + 1) * chunk).min(r.len());
                for (o, &(i, j)) in out.iter_mut().zip(&missing) {
                    *o = diffs[i].dr[a..b]
                        .iter()
                        .zip(&diffs[j].dr[a..b])
                        .map(|(p, q)| p * q)
                        .sum::<f64>();
                }
                for (o, d) in out[missing.len()..].iter_mut().zip(&diffs) {
                    *o = d.dr[a..b]
                        .iter()
                        .zip(&r[a..b])
                        .map(|(p, q)| p * q)
                        .sum::<f64>();
                }
            });
            self.diffs = diffs;
            let sums = sums?;
            COUNTS[0].fetch_add(1, Ordering::Relaxed);
            COUNTS[1].fetch_add(missing.len() as u64, Ordering::Relaxed);
            for (&(i, j), &v) in missing.iter().zip(&sums) {
                self.gram[i][j] = v;
                self.gram[j][i] = v;
            }
            let rhs = &sums[missing.len()..];
            let glue_span = self
                .glue_profile
                .map(|p| p.span(Gc::GlueDistSolveCalls, Gc::GlueDistSolveNanos));
            let (values, largest, gamma) = spectral_solve(|i, j| self.gram[i][j], rhs)?;
            drop(glue_span);
            let glue_span = self
                .glue_profile
                .map(|p| p.span(Gc::GlueDistUpdateCalls, Gc::GlueDistUpdateNanos));
            for (d, g) in self.diffs.iter().zip(&gamma) {
                for (i, n) in next.iter_mut().enumerate() {
                    *n -= g * (d.dx[i] + beta * d.dr[i]);
                }
            }
            drop(glue_span);
            fit = Some(Fit {
                gamma_rhs: gamma.iter().zip(rhs).map(|(g, b)| g * b).sum(),
                values,
                largest,
                gamma,
                rr_fit: None,
            });
        }
        let coulomb = (None, false);

        let (step2, rr) = self.charge(&mut next, mean, x, &r)?;
        let residual_norm = rr.sqrt();
        self.last_step = Some(match fit {
            Some(f) => {
                let rr_fit = f.rr_fit.unwrap_or(rr);
                let smallest = f.values.iter().copied().fold(f64::INFINITY, f64::min);
                PulayStepDiagnostics {
                    differences: m,
                    kept: f.values.iter().filter(|v| **v > 1e-12 * f.largest).count(),
                    gram_condition: if smallest > 0.0 {
                        f.largest / smallest
                    } else {
                        f64::INFINITY
                    },
                    coefficients: f.gamma,
                    residual_norm,
                    predicted_ratio: if rr_fit > 0.0 {
                        ((rr_fit - f.gamma_rhs).max(0.0) / rr_fit).sqrt()
                    } else {
                        0.0
                    },
                    step_norm: step2.sqrt(),
                    coulomb_predicted_ratio: coulomb.0,
                    coulomb_refit: coulomb.1,
                }
            }
            None => PulayStepDiagnostics {
                differences: 0,
                kept: 0,
                gram_condition: 1.0,
                coefficients: Vec::new(),
                residual_norm,
                predicted_ratio: 1.0,
                step_norm: step2.sqrt(),
                coulomb_predicted_ratio: None,
                coulomb_refit: false,
            },
        });
        Ok(next)
    }

    /// Test fixture seeding of plane-mode mixer history.
    #[cfg(test)]
    fn seed_block(
        &mut self,
        dx: &[Vec<f64>],
        dr: &[Vec<f64>],
        mean: f64,
    ) -> Result<(), SolverError> {
        if self.pairs != 0
            || self.seeds != 0
            || dx.len() != dr.len()
            || !mean.is_finite()
            || mean <= 0.0
            || dx.iter().chain(dr).any(|v| v.len() != self.block_len())
        {
            return Err(SolverError::InvalidParameter(
                "seed_block: empty history, equal block pairs and positive mean required",
            ));
        }
        let skip = dx.len().saturating_sub(self.capacity);
        for (x, r) in dx[skip..].iter().zip(&dr[skip..]) {
            self.push_diff(x.iter().copied(), r.iter().copied());
            self.seeds += 1;
        }
        self.mean = Some(mean);
        Ok(())
    }
}

impl<K: Collective> MixerHistory for DistributedPulay<K> {
    fn history_len(&self) -> usize {
        self.pairs
    }
    fn pulay_step(&self) -> Option<PulayStepDiagnostics> {
        self.last_step.clone()
    }
    fn seed_differences(
        &mut self,
        dx: &[Vec<f64>],
        dr: &[Vec<f64>],
        mean: f64,
    ) -> Result<(), SolverError> {
        let size = dx.first().map_or(0, Vec::len);
        #[cfg(test)]
        if self.planes {
            return Err(SolverError::InvalidParameter(
                "plane-mode mixer takes owned blocks (seed_block)",
            ));
        }
        if self.pairs != 0
            || self.seeds != 0
            || dx.len() != dr.len()
            || !mean.is_finite()
            || mean <= 0.0
            || size != self.size
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
        let (lo, hi) = (self.lo, self.hi);
        for (x, r) in dx[skip..].iter().zip(&dr[skip..]) {
            self.push_diff(x[lo..hi].iter().copied(), r[lo..hi].iter().copied());
            self.seeds += 1;
        }
        self.mean = Some(mean);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use dft_mpi::{SimRank, simulate};
    use dft_solvers::mixing::PulayMixer;

    #[test]
    fn configured_beta_is_applied_on_distributed_chunks_and_planes() {
        for ranks in 1..=3 {
            for beta in [0.2, 0.5, 0.8] {
                for planes in [false, true] {
                    let runs = simulate(ranks, move |comm| {
                        let size = 2 * CHUNK;
                        let x = vec![1.0; size];
                        let y: Vec<_> = (0..size)
                            .map(|i| if i % 2 == 0 { 1.2 } else { 0.8 })
                            .collect();
                        let mut m = if planes {
                            // Even plane lengths keep each local residual charge zero.
                            let counts: Vec<_> = (0..ranks)
                                .map(|r| CHUNK / ranks + usize::from(r < CHUNK % ranks))
                                .collect();
                            DistributedPulay::new_planes(comm.clone(), beta, 6, size, 2, &counts)
                                .unwrap()
                        } else {
                            DistributedPulay::new(comm.clone(), beta, 6, size).unwrap()
                        };
                        let next = if planes {
                            m.mix_block(&x[..m.block_len()], &y[..m.block_len()], 1.0)
                                .unwrap()
                        } else {
                            m.mix(&x, &y, 1.0).unwrap()
                        };
                        for (i, v) in next.iter().enumerate() {
                            assert!((v - (x[i] + beta * (y[i] - x[i]))).abs() < 2e-13);
                        }
                        let d = m.last_step.as_ref().unwrap();
                        assert!((d.step_norm / d.residual_norm - beta).abs() < 2e-13);
                        assert_eq!(d.differences, 0);
                        d.predicted_ratio
                    });
                    assert!(runs.iter().all(|r| *r == 1.0));
                }
            }
        }
    }

    fn lcg(s: &mut u64) -> f64 {
        *s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (*s >> 11) as f64 / (1u64 << 53) as f64
    }

    /// A contracting affine map with a fixed point of mean 1 (`size` points).
    struct Model {
        diag: Vec<f64>,
        fixed: Vec<f64>,
    }
    impl Model {
        fn new(size: usize, seed: u64) -> Self {
            let mut s = seed;
            let diag = (0..size).map(|_| 0.2 + 0.7 * lcg(&mut s)).collect();
            let mut fixed: Vec<f64> = (0..size).map(|_| 0.5 + lcg(&mut s)).collect();
            let mean = fixed.iter().sum::<f64>() / size as f64;
            fixed.iter_mut().for_each(|v| *v /= mean);
            Self { diag, fixed }
        }
        fn apply(&self, x: &[f64]) -> Vec<f64> {
            x.iter()
                .zip(&self.diag)
                .zip(&self.fixed)
                .map(|((x, d), f)| f + d * (x - f))
                .collect()
        }
    }

    /// Steps of the mixer: (inputs fed, outputs), with a reset at step 5, a
    /// seeded start and a window of `capacity` differences.
    fn trajectory(
        mixer: &mut dyn MixerHistory,
        model: &Model,
        size: usize,
        seeds: &(Vec<Vec<f64>>, Vec<Vec<f64>>),
    ) -> Vec<Vec<f64>> {
        mixer.seed_differences(&seeds.0, &seeds.1, 1.0).unwrap();
        let mut x: Vec<f64> = (0..size)
            .map(|i| 1.0 + 0.1 * ((i % 7) as f64 - 3.0))
            .collect();
        let mut out = Vec::new();
        for step in 0..10 {
            if step == 5 {
                mixer.reset();
            }
            let y = model.apply(&x);
            x = mixer.mix(&x, &y, 1.0).unwrap();
            out.push(x.clone());
        }
        out
    }

    /// Plane mode (M2-OWNED): blocks of whole planes, no gather; the
    /// concatenated blocks stay in the rounding class of the replicated mixer
    /// and have the same bits for every rank count.
    #[test]
    fn plane_mode_matches_the_replicated_mixer_for_every_rank_count() {
        plane_profile_trajectory("off");
    }
    #[test]
    fn glue_profile_plane_mixer_preserves_trajectory_bits() {
        let reference = plane_profile_trajectory("off");
        assert_eq!(reference, plane_profile_trajectory("counts"));
        assert_eq!(reference, plane_profile_trajectory("timed"));
    }
    fn plane_profile_trajectory(mode: &str) -> Vec<Vec<u64>> {
        let (plane, planes) = (517, 9);
        let size = plane * planes;
        let model = Model::new(size, 3);
        let mut s = 11;
        let seeds: (Vec<Vec<f64>>, Vec<Vec<f64>>) = (0..3)
            .map(|_| {
                let dx: Vec<f64> = (0..size).map(|_| 0.01 * (lcg(&mut s) - 0.5)).collect();
                let dr = dx
                    .iter()
                    .zip(&model.diag)
                    .map(|(d, a)| (a - 1.0) * d)
                    .collect();
                (dx, dr)
            })
            .unzip();
        let mut serial = PulayMixer::new(0.5, 3).unwrap();
        let reference = trajectory(&mut serial, &model, size, &seeds);
        let mut first: Option<Vec<Vec<u64>>> = None;
        for ranks in 1..=7 {
            let counts: Vec<usize> = (0..ranks)
                .map(|r| {
                    let st = |r: usize| r * (planes / ranks) + r.min(planes % ranks);
                    st(r + 1) - st(r)
                })
                .collect();
            let cfg = dft_numerics::switches::RunConfig::from_pairs(&[
                ("glue_profile", if mode == "off" { "off" } else { "on" }),
                ("phase_profile", "on"),
                ("stage_timers", if mode == "timed" { "on" } else { "off" }),
            ])
            .unwrap();
            let runs = simulate(ranks, |comm: &SimRank| {
                let _scope = cfg.enter();
                let rank = dft_mpi::Collective::rank(comm);
                let lo = counts[..rank].iter().sum::<usize>() * plane;
                let hi = lo + counts[rank] * plane;
                let mut m =
                    DistributedPulay::new_planes(comm.clone(), 0.5, 3, size, plane, &counts)
                        .unwrap();
                let block = |v: &Vec<f64>| v[lo..hi].to_vec();
                let sx: Vec<Vec<f64>> = seeds.0.iter().map(block).collect();
                let sr: Vec<Vec<f64>> = seeds.1.iter().map(block).collect();
                m.seed_block(&sx, &sr, 1.0).unwrap();
                let mut x: Vec<f64> = (0..size)
                    .map(|i| 1.0 + 0.1 * ((i % 7) as f64 - 3.0))
                    .collect();
                let mut out = Vec::new();
                for step in 0..10 {
                    if step == 5 {
                        m.reset();
                    }
                    let y = model.apply(&x);
                    let next = m.mix_block(&x[lo..hi], &y[lo..hi], 1.0).unwrap();
                    // Assemble the full vector for the model (test only).
                    let mut full = vec![0.0; size];
                    let c: Vec<usize> = counts.iter().map(|c| c * plane).collect();
                    comm.allgatherv_f64(&next, &c, &mut full).unwrap();
                    x = full;
                    out.push(x.clone());
                }
                out
            });
            let c = cfg.counters();
            assert_eq!(
                c[Gc::GlueDistMixCalls as usize],
                if mode == "off" { 0 } else { 10 * ranks as u64 }
            );
            if mode != "timed" {
                assert_eq!(c[Gc::GlueDistMixNanos as usize], 0);
            }
            for run in &runs {
                for (a, b) in run.iter().zip(&reference) {
                    let worst = a
                        .iter()
                        .zip(b)
                        .fold(0.0_f64, |w, (p, q)| w.max((p - q).abs()));
                    assert!(worst < 1e-11, "P={ranks}: {worst:e}");
                }
                assert_eq!(run, &runs[0]);
            }
            let bits: Vec<Vec<u64>> = runs[0]
                .iter()
                .map(|v| v.iter().map(|x| x.to_bits()).collect())
                .collect();
            match &first {
                None => first = Some(bits),
                Some(f) => assert!(&bits == f, "P={ranks}: rank-count dependent bits"),
            }
        }
        first.unwrap()
    }

    #[test]
    fn distributed_pulay_matches_the_replicated_mixer_for_every_rank_count() {
        // Several chunks with a short last one; capacity 3 slides the window.
        let size = 3 * CHUNK + 517;
        let model = Model::new(size, 7);
        let mut s = 9;
        let seeds: (Vec<Vec<f64>>, Vec<Vec<f64>>) = (0..4)
            .map(|_| {
                let dx: Vec<f64> = (0..size).map(|_| 0.01 * (lcg(&mut s) - 0.5)).collect();
                let dr = dx
                    .iter()
                    .zip(&model.diag)
                    .map(|(d, a)| (a - 1.0) * d)
                    .collect();
                (dx, dr)
            })
            .unzip();
        let mut serial = PulayMixer::new(0.5, 3).unwrap();
        let reference = trajectory(&mut serial, &model, size, &seeds);
        let mut first: Option<Vec<Vec<u64>>> = None;
        for ranks in 1..=7 {
            let runs = simulate(ranks, |comm: &SimRank| {
                let mut m = DistributedPulay::new(comm.clone(), 0.5, 3, size).unwrap();
                trajectory(&mut m, &model, size, &seeds)
            });
            for run in &runs {
                // Rounding class: the chunked dot products differ from the
                // serial row order by ~1e-16 relative; the Pulay least squares
                // (near-dependent differences of a converging trajectory) and
                // the fed-back trajectory lift that to at most 1.3e-12 on
                // densities of order 1 here (measured, no growth over 10 steps).
                for (a, b) in run.iter().zip(&reference) {
                    let worst = a
                        .iter()
                        .zip(b)
                        .fold(0.0_f64, |w, (p, q)| w.max((p - q).abs()));
                    assert!(worst < 1e-11, "P={ranks}: {worst:e}");
                }
                assert_eq!(run, &runs[0], "ranks agree");
            }
            let bits: Vec<Vec<u64>> = runs[0]
                .iter()
                .map(|v| v.iter().map(|x| x.to_bits()).collect())
                .collect();
            match &first {
                None => first = Some(bits),
                Some(f) => assert!(&bits == f, "P={ranks}: rank-count dependent bits"),
            }
        }
    }
}
