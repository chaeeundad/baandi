//! Half-sphere Hartree bilinear form: 8 pi Omega sum Re(a* b)/G^2 [Ry].
//! The Ha form is exactly half. G=0 is carried in vectors but has zero weight.
//! Each batch uses ONE allgather of scalar partials, followed by ascending-rank
//! summation. Fixed local G order; across rank counts differences are rounding
//! only (this does not claim bitwise rank-count independence).
use super::memory::{Component as Mem, bytes, observe};
use super::{Error, invalid, numerical};
use dft_mpi::Collective;
use dft_physics::data::Complex64 as C;
use std::f64::consts::PI;

pub struct Metric {
    identity: u64,
    weights: Vec<f64>,
    local: Vec<f64>,
    gathered: Vec<f64>,
    counts: Vec<usize>,
}
impl Metric {
    pub fn new(
        volume: f64,
        g2: &[f64],
        smooth: &[bool],
        ranks: usize,
        max_batch: usize,
    ) -> Result<Self, Error> {
        if !volume.is_finite()
            || volume <= 0.0
            || g2.len() != smooth.len()
            || ranks == 0
            || max_batch == 0
            || g2.iter().any(|g| !g.is_finite() || *g < 0.0)
        {
            return Err(invalid("metric geometry/shape"));
        }
        let weights = g2
            .iter()
            .zip(smooth)
            .map(|(&g, &low)| {
                if low && g > 0.0 {
                    8.0 * PI * volume / g
                } else {
                    0.0
                }
            })
            .collect();
        let metric = Self {
            identity: {
                static NEXT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
                NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
            },
            weights,
            local: vec![0.; max_batch],
            gathered: vec![0.; ranks * max_batch],
            counts: vec![0; ranks],
        };
        observe(
            Mem::Metric,
            bytes(&metric.weights)
                + bytes(&metric.local)
                + bytes(&metric.gathered)
                + bytes(&metric.counts),
        );
        Ok(metric)
    }
    /// Immutable geometry identity for the mixer's exact Gram cache. A new
    /// metric invalidates it even when the vector lengths happen to agree.
    pub(super) fn identity(&self) -> u64 {
        self.identity
    }
    pub fn batch<K: Collective>(
        &mut self,
        comm: &K,
        pairs: &[(&[C], &[C])],
        out: &mut [f64],
    ) -> Result<(), Error> {
        self.prediction_batch(comm, pairs, out)?
            .then_some(())
            .ok_or_else(|| numerical("nonfinite metric"))
    }
    /// Predictor-only form: nonfinite arithmetic is undefined, while shape and
    /// collective errors still propagate. All ranks inspect the same sums.
    pub fn prediction_batch<K: Collective>(
        &mut self,
        comm: &K,
        pairs: &[(&[C], &[C])],
        out: &mut [f64],
    ) -> Result<bool, Error> {
        let n = pairs.len();
        if n > self.local.len()
            || out.len() != n
            || comm.size() != self.counts.len()
            || pairs
                .iter()
                .any(|(a, b)| a.len() != self.weights.len() || b.len() != self.weights.len())
        {
            return Err(invalid("metric batch shape"));
        }
        for (p, (a, b)) in pairs.iter().enumerate() {
            // Neumaier compensation limits partition-dependent rounding without
            // gathering densities or introducing a rank-count-dependent tree.
            let (mut sum, mut correction) = (0.0_f64, 0.0_f64);
            for ((a, b), w) in a.iter().zip(*b).zip(&self.weights) {
                let v = w * (a.re * b.re + a.im * b.im);
                let next = sum + v;
                correction += if sum.abs() >= v.abs() {
                    (sum - next) + v
                } else {
                    (v - next) + sum
                };
                sum = next;
            }
            self.local[p] = sum + correction;
        }
        self.counts.fill(n);
        comm.allgatherv_f64(
            &self.local[..n],
            &self.counts,
            &mut self.gathered[..n * comm.size()],
        )
        .map_err(numerical)?;
        out.fill(0.0);
        for values in self.gathered[..n * comm.size()].chunks(n.max(1)) {
            for (o, v) in out.iter_mut().zip(values) {
                *o += v;
            }
        }
        Ok(out.iter().all(|v| v.is_finite()))
    }
    pub fn norm<K: Collective>(&mut self, comm: &K, f: &[C]) -> Result<f64, Error> {
        let mut out = [0.0];
        self.batch(comm, &[(f, f)], &mut out)?;
        if out[0] < 0.0 {
            return Err(numerical("negative dr2"));
        }
        Ok(out[0])
    }
}

/// Serial collective with the same adapter calls as distributed execution.
#[derive(Clone, Copy, Default)]
pub struct Serial;
impl Collective for Serial {
    fn rank(&self) -> usize {
        0
    }
    fn size(&self) -> usize {
        1
    }
    fn barrier(&self) -> Result<(), dft_mpi::MpiError> {
        Ok(())
    }
    fn allgatherv_f64(
        &self,
        local: &[f64],
        counts: &[usize],
        out: &mut [f64],
    ) -> Result<(), dft_mpi::MpiError> {
        if counts != [local.len()] || out.len() != local.len() {
            return Err(dft_mpi::MpiError::Shape("serial gather"));
        }
        out.copy_from_slice(local);
        Ok(())
    }
    fn tree_sum_f64(&self, local: &[f64], out: &mut [f64]) -> Result<(), dft_mpi::MpiError> {
        self.allgatherv_f64(local, &[local.len()], out)
    }
    fn bcast_f64(&self, _: &mut [f64], root: usize) -> Result<(), dft_mpi::MpiError> {
        if root == 0 {
            Ok(())
        } else {
            Err(dft_mpi::MpiError::Shape("serial root"))
        }
    }
    fn alltoall_f64(
        &self,
        send: &[f64],
        block: usize,
        recv: &mut [f64],
    ) -> Result<(), dft_mpi::MpiError> {
        self.allgatherv_f64(send, &[block], recv)
    }
}
