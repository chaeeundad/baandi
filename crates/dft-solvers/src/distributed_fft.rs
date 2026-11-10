//! G-vector (stick) distributed pruned 3D FFT for packed real pairs (MPI stage M2a).
//!
//! Placement: this module lives in dft-solvers because that crate already depends
//! on both `dft_mpi::Collective` and dft-numerics, and dft-numerics must stay free
//! of the MPI layer. Pairing of two real columns into one complex vector is done
//! upstream (`PairMap`); this plan transforms one complex coefficient vector.
//!
//! Conventions follow `PairFusedPlan`: the first transformed axis is the one with
//! fewer occupied sticks (ties: x), z is last. Sticks are columns along the first
//! axis, keyed by (z, other axis).
//!
//! Partition rule (deterministic in shape, indices, ranks): sticks sorted by
//! length descending (ties: lower (z, other) slot) are assigned one by one to the
//! rank with the smallest accumulated length (ties: lower rank). Each rank owns
//! the sticks assigned to it; its local coefficients are the entries of the global
//! coefficient vector that fall on those sticks, in ascending global position
//! (`local_positions`). Real-space planes are contiguous ranges of the first axis
//! (sizes differ by at most one; ranks beyond n0 own none).
//!
//! Transpose: one `alltoall_f64` per direction with equal blocks of
//! 2*max_sticks*max_planes f64; shorter (rank, plane-range) pairs are zero-padded
//! and the padding is never read. Per-plane work reuses the serial order of
//! operations (clear, fill, batched 1D FFT along the second axis, 1D FFT along
//! z), so each plane's slab is expected to be bit-identical to the serial plan's
//! for any rank count (tests assert bit identity for one rank only). The 1/N
//! factor is applied at the end of the reciprocal direction.
//! No allocation after construction; no unsafe; no logging.
use dft_mpi::{Collective, MpiError};
use dft_numerics::FftError;
use num_complex::Complex64 as C;
use rustfft::{Fft, FftPlanner};
use std::{collections::BTreeSet, sync::Arc};
use thiserror::Error;

#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum DistFftError {
    #[error("FFT: {0}")]
    Fft(#[from] FftError),
    #[error("MPI: {0}")]
    Mpi(#[from] MpiError),
}

struct Partition {
    first: usize,
    n: [usize; 3],
    sticks: Vec<(usize, usize)>,
    planes: Vec<usize>,
    slots: Vec<usize>,
    owner: Vec<usize>,
    local_index: Vec<usize>,
    counts: Vec<usize>,
}

/// Stick partition shared by `new` and `position_owners` (inputs already validated).
fn partition(shape: [usize; 3], indices: &[[usize; 3]], ranks: usize) -> Partition {
    let supports: [_; 2] = std::array::from_fn(|a| {
        indices
            .iter()
            .map(|p| (p[2], p[1 - a]))
            .collect::<BTreeSet<_>>()
    });
    let first = usize::from(supports[1].len() < supports[0].len());
    let n = [shape[first], shape[1 - first], shape[2]];
    let sticks: Vec<(usize, usize)> = supports[first].iter().copied().collect();
    let planes: Vec<usize> = sticks
        .iter()
        .map(|p| p.0)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let mut len = vec![0usize; sticks.len()];
    let slots: Vec<usize> = indices
        .iter()
        .map(|p| {
            let s = sticks.binary_search(&(p[2], p[1 - first])).unwrap();
            len[s] += 1;
            s
        })
        .collect();
    // Balanced greedy assignment by stick length.
    let mut order: Vec<usize> = (0..sticks.len()).collect();
    order.sort_by(|&a, &b| len[b].cmp(&len[a]).then(a.cmp(&b)));
    let mut load = vec![0usize; ranks];
    let mut owner = vec![0usize; sticks.len()];
    for &s in &order {
        let r = (0..ranks).min_by_key(|&r| (load[r], r)).unwrap();
        owner[s] = r;
        load[r] += len[s];
    }
    let mut local_index = vec![0usize; sticks.len()];
    let mut counts = vec![0usize; ranks];
    for s in 0..sticks.len() {
        local_index[s] = counts[owner[s]];
        counts[owner[s]] += 1;
    }
    Partition {
        first,
        n,
        sticks,
        planes,
        slots,
        owner,
        local_index,
        counts,
    }
}

/// Owning rank of every global coefficient position (same partition as `new`).
pub fn position_owners(
    shape: [usize; 3],
    indices: &[[usize; 3]],
    ranks: usize,
) -> Result<Vec<usize>, DistFftError> {
    if ranks == 0 || shape.contains(&0) {
        return Err(FftError::InvalidGrid.into());
    }
    if indices.iter().any(|p| (0..3).any(|d| p[d] >= shape[d])) {
        return Err(FftError::IncompatibleArray.into());
    }
    let part = partition(shape, indices, ranks);
    Ok(part.slots.iter().map(|&s| part.owner[s]).collect())
}

pub struct DistributedPairFft {
    shape: [usize; 3],
    first: usize,
    n: [usize; 3],
    total: usize,
    ranks: usize,
    local_positions: Vec<usize>,
    // local stick slot * n0 + first-axis index, one per local coefficient.
    scatter: Vec<usize>,
    local_sticks: usize,
    max_planes: usize,
    block: usize, // complex elements per rank block
    plane_starts: Vec<usize>,
    rank: usize,
    planes: Vec<usize>,
    // (index in the middle batch, complex offset of the stick in the rank-block layout).
    occupied: Vec<(usize, usize)>,
    sticks: Vec<C>,
    middle: Vec<C>,
    slab: Vec<C>,
    scratch: Vec<C>,
    send: Vec<f64>,
    recv: Vec<f64>,
    inverse: [Arc<dyn Fft<f64>>; 3],
    forward: [Arc<dyn Fft<f64>>; 3],
}

impl DistributedPairFft {
    pub fn new(
        shape: [usize; 3],
        indices: &[[usize; 3]],
        rank: usize,
        ranks: usize,
    ) -> Result<Self, DistFftError> {
        let total = shape
            .iter()
            .try_fold(1usize, |a, &b| a.checked_mul(b))
            .ok_or(FftError::InvalidGrid)?;
        if total == 0 || total > isize::MAX as usize / size_of::<C>() || ranks == 0 || rank >= ranks
        {
            return Err(FftError::InvalidGrid.into());
        }
        if indices.iter().any(|p| (0..3).any(|d| p[d] >= shape[d])) {
            return Err(FftError::IncompatibleArray.into());
        }
        let Partition {
            first,
            n,
            sticks,
            planes,
            slots,
            owner,
            local_index,
            counts,
        } = partition(shape, indices, ranks);
        let max_sticks = counts.iter().copied().max().unwrap_or(0);
        let plane_starts: Vec<usize> = (0..=ranks)
            .map(|r| r * (n[0] / ranks) + r.min(n[0] % ranks))
            .collect();
        let max_planes = plane_starts
            .windows(2)
            .map(|w| w[1] - w[0])
            .max()
            .unwrap_or(0);
        let block = max_sticks * max_planes;
        let mut local_positions = Vec::new();
        let mut scatter = Vec::new();
        for (g, (p, &s)) in indices.iter().zip(&slots).enumerate() {
            if owner[s] == rank {
                local_positions.push(g);
                scatter.push(local_index[s] * n[0] + p[first]);
            }
        }
        let occupied = sticks
            .iter()
            .enumerate()
            .map(|(s, &(z, y))| {
                (
                    planes.binary_search(&z).unwrap() * n[1] + y,
                    owner[s] * block + local_index[s] * max_planes,
                )
            })
            .collect::<Vec<_>>();
        let mut planner = FftPlanner::new();
        let inverse = n.map(|d| planner.plan_fft_inverse(d));
        let forward = n.map(|d| planner.plan_fft_forward(d));
        let scratch = inverse
            .iter()
            .chain(&forward)
            .map(|p| p.get_inplace_scratch_len())
            .max()
            .unwrap();
        let local_sticks = counts[rank];
        Ok(Self {
            shape,
            first,
            n,
            total,
            ranks,
            local_positions,
            scatter,
            local_sticks,
            max_planes,
            block,
            rank,
            plane_starts,
            middle: vec![C::default(); planes.len() * n[1]],
            sticks: vec![C::default(); local_sticks * n[0]],
            slab: vec![C::default(); n[1] * n[2]],
            scratch: vec![C::default(); scratch],
            planes,
            occupied,
            send: vec![0.; 2 * block * ranks],
            recv: vec![0.; 2 * block * ranks],
            inverse,
            forward,
        })
    }

    /// Positions in the global coefficient vector owned by this rank (ascending);
    /// local coefficient vectors follow this order.
    pub fn local_positions(&self) -> &[usize] {
        &self.local_positions
    }
    pub fn local_len(&self) -> usize {
        self.scatter.len()
    }
    pub fn shape(&self) -> [usize; 3] {
        self.shape
    }
    /// Number of owned real-space planes (contiguous along the first FFT axis).
    pub fn plane_count(&self) -> usize {
        self.plane_starts[self.rank + 1] - self.plane_starts[self.rank]
    }
    /// Grid index of owned plane i along grid axis `first_axis()`.
    pub fn plane_x(&self, i: usize) -> usize {
        self.plane_starts[self.rank] + i
    }
    pub fn first_axis(&self) -> usize {
        self.first
    }
    /// Offset into a row-major (shape) grid of the z-line at (plane i, row y); the
    /// line has `shape[2]` entries and matches `plane()[y*nz..]`.
    pub fn grid_start(&self, i: usize, y: usize) -> usize {
        let x = self.plane_x(i);
        if self.first == 0 {
            (x * self.shape[1] + y) * self.n[2]
        } else {
            (y * self.shape[1] + x) * self.n[2]
        }
    }
    /// Extent along the second FFT axis (rows of `plane`).
    pub fn rows(&self) -> usize {
        self.n[1]
    }
    /// Current plane slab, row-major [y][z] (set by `to_real_plane`).
    pub fn plane(&self) -> &[C] {
        &self.slab
    }
    pub fn plane_mut(&mut self) -> &mut [C] {
        &mut self.slab
    }

    fn all_to_all<K: Collective>(&mut self, comm: &K) -> Result<(), DistFftError> {
        if comm.size() != self.ranks || comm.rank() != self.rank {
            return Err(MpiError::Shape("communicator differs from the plan").into());
        }
        if self.block > 0 {
            comm.alltoall_f64(&self.send, 2 * self.block, &mut self.recv)?;
        }
        Ok(())
    }

    /// Step 1 of G->r: 1D transforms on local sticks and the transpose.
    pub fn begin_to_real<K: Collective>(&mut self, comm: &K, c: &[C]) -> Result<(), DistFftError> {
        if c.len() != self.scatter.len() {
            return Err(FftError::IncompatibleArray.into());
        }
        self.sticks.fill(C::default());
        for (&i, &v) in self.scatter.iter().zip(c) {
            self.sticks[i] = v;
        }
        if !self.sticks.is_empty() {
            self.inverse[0].process_with_scratch(&mut self.sticks, &mut self.scratch);
        }
        let (nx, mp) = (self.n[0], self.max_planes);
        for q in 0..self.ranks {
            let (x0, x1) = (self.plane_starts[q], self.plane_starts[q + 1]);
            for sp in 0..self.local_sticks {
                let base = 2 * (q * self.block + sp * mp);
                for (xi, x) in (x0..x1).enumerate() {
                    let v = self.sticks[sp * nx + x];
                    self.send[base + 2 * xi] = v.re;
                    self.send[base + 2 * xi + 1] = v.im;
                }
            }
        }
        self.all_to_all(comm)
    }

    /// Step 2 of G->r for owned plane i: 2D transform; result in `plane()`.
    pub fn to_real_plane(&mut self, i: usize) {
        let [_, ny, nz] = self.n;
        self.middle.fill(C::default());
        for &(m, off) in &self.occupied {
            let k = 2 * (off + i);
            self.middle[m] = C::new(self.recv[k], self.recv[k + 1]);
        }
        if !self.middle.is_empty() {
            self.inverse[1].process_with_scratch(&mut self.middle, &mut self.scratch);
        }
        self.slab.fill(C::default());
        for (k, &z) in self.planes.iter().enumerate() {
            for y in 0..ny {
                self.slab[y * nz + z] = self.middle[k * ny + y];
            }
        }
        self.inverse[2].process_with_scratch(&mut self.slab, &mut self.scratch);
    }

    /// r->G for owned plane i from `plane()`; results are staged for the transpose.
    pub fn to_reciprocal_plane(&mut self, i: usize) {
        let [_, ny, nz] = self.n;
        self.forward[2].process_with_scratch(&mut self.slab, &mut self.scratch);
        for (k, &z) in self.planes.iter().enumerate() {
            for y in 0..ny {
                self.middle[k * ny + y] = self.slab[y * nz + z];
            }
        }
        if !self.middle.is_empty() {
            self.forward[1].process_with_scratch(&mut self.middle, &mut self.scratch);
        }
        for &(m, off) in &self.occupied {
            let k = 2 * (off + i);
            self.send[k] = self.middle[m].re;
            self.send[k + 1] = self.middle[m].im;
        }
    }

    /// Transpose back, 1D transforms on local sticks, 1/N scaling; writes local coefficients.
    pub fn finish_to_reciprocal<K: Collective>(
        &mut self,
        comm: &K,
        out: &mut [C],
    ) -> Result<(), DistFftError> {
        if out.len() != self.scatter.len() {
            return Err(FftError::IncompatibleArray.into());
        }
        self.all_to_all(comm)?;
        let (nx, mp) = (self.n[0], self.max_planes);
        for q in 0..self.ranks {
            let x0 = self.plane_starts[q];
            for sp in 0..self.local_sticks {
                let base = 2 * (q * self.block + sp * mp);
                for xi in 0..self.plane_starts[q + 1] - x0 {
                    self.sticks[sp * nx + x0 + xi] =
                        C::new(self.recv[base + 2 * xi], self.recv[base + 2 * xi + 1]);
                }
            }
        }
        if !self.sticks.is_empty() {
            self.forward[0].process_with_scratch(&mut self.sticks, &mut self.scratch);
        }
        for (o, &i) in out.iter_mut().zip(&self.scatter) {
            *o = self.sticks[i] / self.total as f64;
        }
        Ok(())
    }

    /// Local potential application: G -> r, multiply by `potential` (full row-major
    /// grid, read on owned planes only), r -> G. Equals `PairFusedPlan::apply_local`.
    pub fn apply_local<K: Collective>(
        &mut self,
        comm: &K,
        c: &[C],
        potential: &[f64],
        out: &mut [C],
    ) -> Result<(), DistFftError> {
        if potential.len() != self.total {
            return Err(FftError::IncompatibleArray.into());
        }
        self.begin_to_real(comm, c)?;
        let (ny, nz) = (self.n[1], self.n[2]);
        for i in 0..self.plane_count() {
            self.to_real_plane(i);
            for y in 0..ny {
                let start = self.grid_start(i, y);
                for (a, &v) in self.slab[y * nz..(y + 1) * nz]
                    .iter_mut()
                    .zip(&potential[start..start + nz])
                {
                    *a *= v;
                }
            }
            self.to_reciprocal_plane(i);
        }
        self.finish_to_reciprocal(comm, out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use dft_mpi::simulate;
    use dft_numerics::fused::FusedOperator;
    use dft_numerics::{Fft3d, FftPlan3d, PairFusedPlan};
    use ndarray::Array3;

    fn lcg(s: &mut u64) -> f64 {
        *s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*s >> 11) as f64 / (1u64 << 53) as f64) - 0.5
    }

    fn sphere(shape: [usize; 3], radius2: f64, shift: [f64; 3]) -> Vec<[usize; 3]> {
        let f = |i: usize, n: usize| {
            if i <= n / 2 {
                i as f64
            } else {
                i as f64 - n as f64
            }
        };
        let mut v = Vec::new();
        for i in 0..shape[0] {
            for j in 0..shape[1] {
                for k in 0..shape[2] {
                    let g = [
                        f(i, shape[0]) + shift[0],
                        f(j, shape[1]) + shift[1],
                        f(k, shape[2]) + shift[2],
                    ];
                    if g.iter().map(|x| x * x).sum::<f64>() <= radius2 {
                        v.push([i, j, k]);
                    }
                }
            }
        }
        v
    }

    type Case = ([usize; 3], Vec<[usize; 3]>);
    fn cases() -> Vec<Case> {
        let skew: Vec<[usize; 3]> = sphere([12, 14, 16], 20., [0.5, 0.3, 0.])
            .into_iter()
            .filter(|p| (p[0] + 2 * p[1] + p[2]) % 3 != 0)
            .collect();
        vec![
            ([16, 16, 16], sphere([16, 16, 16], 16., [0.; 3])),
            ([15, 15, 45], sphere([15, 15, 45], 49., [0.; 3])),
            ([12, 14, 16], skew),
        ]
    }

    fn random_c(n: usize, seed: u64) -> Vec<C> {
        let mut s = seed;
        (0..n).map(|_| C::new(lcg(&mut s), lcg(&mut s))).collect()
    }

    fn rel(a: &[C], b: &[C]) -> f64 {
        let num = a
            .iter()
            .zip(b)
            .map(|(x, y)| (x - y).norm())
            .fold(0., f64::max);
        let den = b.iter().map(|x| x.norm()).fold(0., f64::max);
        num / den
    }

    /// Distributed G->r field gathered into a full grid.
    fn dist_field(shape: [usize; 3], idx: &[[usize; 3]], c: &[C], ranks: usize) -> Vec<C> {
        let res = simulate(ranks, |comm| {
            let mut p = DistributedPairFft::new(shape, idx, comm.rank(), ranks).unwrap();
            let local: Vec<C> = p.local_positions().iter().map(|&g| c[g]).collect();
            p.begin_to_real(comm, &local).unwrap();
            let mut lines = Vec::new();
            for i in 0..p.plane_count() {
                p.to_real_plane(i);
                for y in 0..p.rows() {
                    lines.push((
                        p.grid_start(i, y),
                        p.plane()[y * shape[2]..(y + 1) * shape[2]].to_vec(),
                    ));
                }
            }
            lines
        });
        let mut grid = vec![C::default(); shape.iter().product()];
        let mut covered = 0;
        for lines in &res {
            for (s, l) in lines {
                grid[*s..*s + l.len()].copy_from_slice(l);
                covered += l.len();
            }
        }
        assert_eq!(covered, grid.len());
        grid
    }

    #[test]
    fn partition_is_exact_and_balanced() {
        for (shape, idx) in cases() {
            for ranks in 1..=8 {
                let mut seen = vec![0; idx.len()];
                let mut lens = Vec::new();
                for r in 0..ranks {
                    let p = DistributedPairFft::new(shape, &idx, r, ranks).unwrap();
                    lens.push(p.local_len());
                    for &g in p.local_positions() {
                        seen[g] += 1;
                    }
                }
                assert!(seen.iter().all(|&s| s == 1));
                let max = *lens.iter().max().unwrap();
                let min = *lens.iter().min().unwrap();
                assert!(max - min <= shape[0].max(shape[1]), "{lens:?}");
            }
        }
    }

    #[test]
    fn forward_matches_serial_field() {
        for (shape, idx) in cases() {
            let n = shape.iter().product::<usize>() as f64;
            let c = random_c(idx.len(), 7);
            let mut full = Array3::<C>::zeros(shape);
            for (p, v) in idx.iter().zip(&c) {
                full[[p[0], p[1], p[2]]] = *v;
            }
            let mut serial = FftPlan3d::new(shape).unwrap();
            let mut field = Array3::<C>::zeros(shape);
            serial.to_real(&full, &mut field).unwrap();
            let reference: Vec<C> = field.iter().copied().collect();
            let mut worst = 0f64;
            for ranks in 1..=8 {
                let grid = dist_field(shape, &idx, &c, ranks);
                let e = rel(&grid, &reference);
                worst = worst.max(e);
                assert!(e < 8. * n * f64::EPSILON, "{shape:?} ranks {ranks}: {e:e}");
                assert_eq!(
                    grid,
                    dist_field(shape, &idx, &c, ranks),
                    "bit reproducibility at fixed rank count"
                );
            }
            eprintln!(
                "forward {shape:?} worst rel {worst:e} (limit {:e})",
                8. * n * f64::EPSILON
            );
        }
    }

    #[test]
    fn roundtrip_and_local_potential_match_serial() {
        for (shape, idx) in cases() {
            let n = shape.iter().product::<usize>();
            let c = random_c(idx.len(), 11);
            let mut s = 99;
            let pot: Vec<f64> = (0..n).map(|_| 1.0 + lcg(&mut s)).collect();
            let ones = vec![1.0; n];
            let mut serial = PairFusedPlan::new(shape, &idx).unwrap();
            let mut want = vec![C::default(); idx.len()];
            serial.apply_local(&c, &pot, &mut want).unwrap();
            for ranks in 1..=8 {
                let parts = simulate(ranks, |comm| {
                    let mut p = DistributedPairFft::new(shape, &idx, comm.rank(), ranks).unwrap();
                    let local: Vec<C> = p.local_positions().iter().map(|&g| c[g]).collect();
                    let mut rt = vec![C::default(); local.len()];
                    p.apply_local(comm, &local, &ones, &mut rt).unwrap();
                    let mut out = vec![C::default(); local.len()];
                    p.apply_local(comm, &local, &pot, &mut out).unwrap();
                    (p.local_positions().to_vec(), rt, out)
                });
                let mut got = vec![C::default(); idx.len()];
                let mut rt = vec![C::default(); idx.len()];
                for (pos, r, o) in &parts {
                    for (k, &g) in pos.iter().enumerate() {
                        got[g] = o[k];
                        rt[g] = r[k];
                    }
                }
                let e_rt = rel(&rt, &c);
                let e_pot = rel(&got, &want);
                eprintln!(
                    "{shape:?} ranks {ranks}: roundtrip {e_rt:e} potential {e_pot:e} bitwise {}",
                    got == want
                );
                assert!(e_rt < 4. * n as f64 * f64::EPSILON);
                assert!(e_pot < 1e-12);
                if ranks == 1 {
                    assert_eq!(got, want, "one rank is bit-identical to PairFusedPlan");
                }
            }
        }
    }

    #[test]
    fn no_reallocation_after_construction() {
        let (shape, idx) = cases().remove(0);
        simulate(1, |comm| {
            let mut p = DistributedPairFft::new(shape, &idx, 0, 1).unwrap();
            let c = random_c(idx.len(), 3);
            let pot = vec![1.; 4096];
            let mut out = vec![C::default(); idx.len()];
            let ptrs = |p: &DistributedPairFft| {
                [
                    p.send.as_ptr() as usize,
                    p.recv.as_ptr() as usize,
                    p.sticks.as_ptr() as usize,
                    p.middle.as_ptr() as usize,
                ]
            };
            let before = ptrs(&p);
            for _ in 0..3 {
                p.apply_local(comm, &c, &pot, &mut out).unwrap();
            }
            assert_eq!(before, ptrs(&p));
        });
    }
}
