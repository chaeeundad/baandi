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
//! rank with the smallest accumulated length (ties: lower rank); with
//! `BAANDI_M2_STICKS=mirror` (MPI-SCALING-GAP) a stick and its Gamma mirror stick
//! are assigned together as one unit (no pair representative is split, so the
//! mirror exchanges of the real Gamma operator vanish). Each rank owns
//! the sticks assigned to it; its local coefficients are the entries of the global
//! coefficient vector that fall on those sticks, in ascending global position
//! (`local_positions`). Real-space planes are contiguous ranges of the first axis
//! (sizes differ by at most one; ranks beyond n0 own none).
//!
//! Transpose: one `alltoall_f64` per direction with equal blocks of
//! 2*max_sticks*max_planes f64 per transformed vector; shorter (rank, plane-range)
//! pairs are zero-padded and the padding is never read. Batches (M2d): the `_batch`
//! entry points transform np vectors with ONE alltoall per direction (block
//! 2*max_sticks*max_planes*np, destination-major, vector-minor); every vector's
//! arithmetic is exactly that of the single-vector path, so a batch is bit-identical
//! to np single calls. Transpose buffers grow to the widest batch seen and are then
//! reused (no allocation once the widest batch has occurred). Per-plane work reuses the serial order of
//! operations (clear, fill, batched 1D FFT along the second axis, 1D FFT along
//! z), so each plane's slab is expected to be bit-identical to the serial plan's
//! for any rank count (tests assert bit identity for one rank only). The 1/N
//! factor is applied at the end of the reciprocal direction.
//! 1D kernels (MPI entry cost, `BAANDI_M2_FFT`): rustfft by default; `mkl` selects
//! MKL DFTI batch lanes (one descriptor per lane buffer, as `PairFusedPlan<MklLanes>`)
//! on MKL builds, a rounding change; pure builds always use rustfft.
//! No allocation after construction; no unsafe; no logging.
use dft_mpi::{Collective, MpiError};
use dft_numerics::FftError;
use dft_numerics::fused_pair_plan::BackendLanes;
use num_complex::Complex64 as C;
use std::collections::BTreeSet;
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

thread_local! {
    static MIRROR_OVERRIDE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}
/// Force the stick rule of the calling thread (tests with simulated ranks; every
/// rank thread must set the same value); `None` restores `BAANDI_M2_STICKS`.
pub fn set_mirror_sticks_override(mirror: Option<bool>) {
    MIRROR_OVERRIDE.with(|o| o.set(mirror));
}
/// True when the mirror-paired stick partition is in force on this thread.
pub fn mirror_sticks() -> bool {
    MIRROR_OVERRIDE.with(|o| o.get()).unwrap_or_else(|| {
        dft_numerics::switches::choice(dft_numerics::switches::Id::M2_STICKS) == "mirror"
    })
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
    let mut owner = vec![usize::MAX; sticks.len()];
    // `BAANDI_M2_STICKS=mirror` (MPI-SCALING-GAP): a stick (z, o) and its mirror
    // ((-z) mod n2, (-o) mod n1) are assigned to one rank as a unit, so the two
    // positions +G and -G of every Gamma pair representative live on one rank.
    // A stick without a mirror in the set (non-symmetric index sets) or its own
    // mirror is assigned alone, exactly as in the greedy rule.
    let mirror = mirror_sticks();
    let n_other = n[1];
    let n_z = n[2];
    for &s in &order {
        if owner[s] != usize::MAX {
            continue;
        }
        let r = (0..ranks).min_by_key(|&r| (load[r], r)).unwrap();
        owner[s] = r;
        load[r] += len[s];
        if mirror {
            let (z, o) = sticks[s];
            let key = ((n_z - z) % n_z, (n_other - o) % n_other);
            if let Ok(m) = sticks.binary_search(&key)
                && owner[m] == usize::MAX
            {
                owner[m] = r;
                load[r] += len[m];
            }
        }
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
    // (index in the middle batch, owner rank, complex offset of the stick inside
    // one vector's rank block).
    occupied: Vec<(usize, usize, usize)>,
    sticks: Vec<C>,
    middle: Vec<C>,
    slab: Vec<C>,
    send: Vec<f64>,
    recv: Vec<f64>,
    /// 1D kernels: rustfft (default) or MKL DFTI batch lanes (`BAANDI_M2_FFT=mkl`
    /// on MKL builds). Lane buffers: local sticks, middle batch, plane slab.
    lanes: BackendLanes,
    profile: dft_numerics::real::fft_profile::Probe,
    kept_fuse: bool,
}

impl DistributedPairFft {
    /// Plan with the 1D kernels of `BAANDI_M2_FFT` (read once per process).
    pub fn new(
        shape: [usize; 3],
        indices: &[[usize; 3]],
        rank: usize,
        ranks: usize,
    ) -> Result<Self, DistFftError> {
        let vendor = dft_numerics::switches::choice(dft_numerics::switches::Id::M2_FFT) == "mkl";
        Self::with_lanes(shape, indices, rank, ranks, vendor)
    }

    /// Plan with rustfft lanes (`vendor` false) or MKL DFTI batch lanes (`vendor`
    /// true on an MKL build; a pure build has none and uses rustfft, bit-identical).
    pub fn with_lanes(
        shape: [usize; 3],
        indices: &[[usize; 3]],
        rank: usize,
        ranks: usize,
        vendor: bool,
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
                    owner[s],
                    local_index[s] * max_planes,
                )
            })
            .collect::<Vec<_>>();
        let local_sticks = counts[rank];
        let lanes = BackendLanes::plan(n, [local_sticks, planes.len(), n[1]], vendor)?;
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
            planes,
            occupied,
            send: vec![0.; 2 * block * ranks],
            recv: vec![0.; 2 * block * ranks],
            lanes,
            profile: dft_numerics::real::fft_profile::Probe::new(n),
            kept_fuse: dft_numerics::switches::on(dft_numerics::switches::Id::FFT_KEPT_FUSE),
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
    /// True when the MKL DFTI lanes are active (`BAANDI_M2_FFT=mkl` on an MKL build).
    pub fn vendor_lanes(&self) -> bool {
        self.lanes.is_vendor()
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

    fn all_to_all<K: Collective>(&mut self, comm: &K, np: usize) -> Result<(), DistFftError> {
        if comm.size() != self.ranks || comm.rank() != self.rank {
            return Err(MpiError::Shape("communicator differs from the plan").into());
        }
        if self.block > 0 {
            let len = 2 * self.block * self.ranks * np;
            comm.alltoall_f64(
                &self.send[..len],
                2 * self.block * np,
                &mut self.recv[..len],
            )?;
        }
        Ok(())
    }

    /// Grow the transpose buffers to batches of `np` vectors (no-op once grown).
    fn reserve_batch(&mut self, np: usize) {
        let len = 2 * self.block * self.ranks * np;
        if self.send.len() < len {
            self.send.resize(len, 0.);
            self.recv.resize(len, 0.);
        }
    }

    /// Complex offset of vector t's block for rank q in a batch of `np`.
    fn base(&self, q: usize, t: usize, np: usize) -> usize {
        (q * np + t) * self.block
    }

    /// Step 1 of G->r: 1D transforms on local sticks and the transpose.
    pub fn begin_to_real<K: Collective>(&mut self, comm: &K, c: &[C]) -> Result<(), DistFftError> {
        self.begin_to_real_batch(comm, c, 1)
    }

    /// Step 1 of G->r for `np` vectors stored one after another in `c` (each
    /// `local_len` long): ONE transpose for all of them.
    pub fn begin_to_real_batch<K: Collective>(
        &mut self,
        comm: &K,
        c: &[C],
        np: usize,
    ) -> Result<(), DistFftError> {
        let ll = self.scatter.len();
        if c.len() != ll * np {
            return Err(FftError::IncompatibleArray.into());
        }
        self.reserve_batch(np);
        let (nx, mp) = (self.n[0], self.max_planes);
        for t in 0..np {
            let mark = self.profile.start();
            self.sticks.fill(C::default());
            for (&i, &v) in self.scatter.iter().zip(&c[t * ll..(t + 1) * ll]) {
                self.sticks[i] = v;
            }
            self.profile.finish(
                mark,
                0,
                self.scatter.len(),
                16 * self.sticks.len() + 32 * self.scatter.len(),
                None,
            );
            if !self.sticks.is_empty() {
                let mark = self.profile.start();
                self.lanes.run(0, true, &mut self.sticks)?;
                self.profile
                    .finish(mark, 1, self.sticks.len(), 32 * self.sticks.len(), Some(0));
            }
            let mark = self.profile.start();
            for q in 0..self.ranks {
                let (x0, x1) = (self.plane_starts[q], self.plane_starts[q + 1]);
                let qb = self.base(q, t, np);
                for sp in 0..self.local_sticks {
                    let base = 2 * (qb + sp * mp);
                    for (xi, x) in (x0..x1).enumerate() {
                        let v = self.sticks[sp * nx + x];
                        self.send[base + 2 * xi] = v.re;
                        self.send[base + 2 * xi + 1] = v.im;
                    }
                }
            }
            self.profile
                .finish(mark, 2, self.sticks.len(), 32 * self.sticks.len(), None);
        }
        let mark = self.profile.start();
        if self.block > 0 {
            dft_numerics::stage_clock::arrival::site(104);
        }
        let result = self.all_to_all(comm, np);
        if self.block > 0 {
            self.profile.finish(
                mark,
                3,
                self.block * self.ranks * np,
                16 * self.block * self.ranks * np,
                None,
            );
        }
        result
    }

    /// Step 2 of G->r for owned plane i: 2D transform; result in `plane()`.
    pub fn to_real_plane(&mut self, i: usize) -> Result<(), DistFftError> {
        self.to_real_plane_in(0, 1, i)
    }

    /// Step 2 of G->r for vector t of a batch of `np`, owned plane i.
    pub fn to_real_plane_in(&mut self, t: usize, np: usize, i: usize) -> Result<(), DistFftError> {
        let [_, ny, nz] = self.n;
        let mark = self.profile.start();
        self.middle.fill(C::default());
        for &(m, owner, off) in &self.occupied {
            let k = 2 * ((owner * np + t) * self.block + off + i);
            self.middle[m] = C::new(self.recv[k], self.recv[k + 1]);
        }
        self.profile.finish(
            mark,
            4,
            self.occupied.len(),
            16 * self.middle.len() + 32 * self.occupied.len(),
            None,
        );
        if !self.middle.is_empty() {
            let mark = self.profile.start();
            self.lanes.run(1, true, &mut self.middle)?;
            self.profile
                .finish(mark, 5, self.middle.len(), 32 * self.middle.len(), Some(1));
        }
        let mark = self.profile.start();
        self.slab.fill(C::default());
        for (k, &z) in self.planes.iter().enumerate() {
            for y in 0..ny {
                self.slab[y * nz + z] = self.middle[k * ny + y];
            }
        }
        self.profile.finish(
            mark,
            6,
            self.middle.len(),
            16 * self.slab.len() + 32 * self.middle.len(),
            None,
        );
        let mark = self.profile.start();
        self.lanes.run(2, true, &mut self.slab)?;
        self.profile
            .finish(mark, 7, self.slab.len(), 32 * self.slab.len(), Some(2));
        Ok(())
    }

    /// r->G for owned plane i from `plane()`; results are staged for the transpose.
    pub fn to_reciprocal_plane(&mut self, i: usize) -> Result<(), DistFftError> {
        self.to_reciprocal_plane_in(0, 1, i)
    }

    /// r->G of vector t of a batch of `np`, owned plane i (staged for the transpose).
    pub fn to_reciprocal_plane_in(
        &mut self,
        t: usize,
        np: usize,
        i: usize,
    ) -> Result<(), DistFftError> {
        let [_, ny, nz] = self.n;
        let mark = self.profile.start();
        self.lanes.run(2, false, &mut self.slab)?;
        self.profile
            .finish(mark, 9, self.slab.len(), 32 * self.slab.len(), Some(2));
        let mark = self.profile.start();
        for (k, &z) in self.planes.iter().enumerate() {
            for y in 0..ny {
                self.middle[k * ny + y] = self.slab[y * nz + z];
            }
        }
        self.profile
            .finish(mark, 10, self.middle.len(), 32 * self.middle.len(), None);
        if !self.middle.is_empty() {
            let mark = self.profile.start();
            self.lanes.run(1, false, &mut self.middle)?;
            self.profile
                .finish(mark, 11, self.middle.len(), 32 * self.middle.len(), Some(1));
        }
        let mark = self.profile.start();
        for &(m, owner, off) in &self.occupied {
            let k = 2 * ((owner * np + t) * self.block + off + i);
            self.send[k] = self.middle[m].re;
            self.send[k + 1] = self.middle[m].im;
        }
        self.profile.finish(
            mark,
            12,
            self.occupied.len(),
            32 * self.occupied.len(),
            None,
        );
        Ok(())
    }

    /// Transpose back, 1D transforms on local sticks, 1/N scaling; writes local coefficients.
    pub fn finish_to_reciprocal<K: Collective>(
        &mut self,
        comm: &K,
        out: &mut [C],
    ) -> Result<(), DistFftError> {
        self.finish_to_reciprocal_batch(comm, out, 1)
    }

    /// Batched `finish_to_reciprocal`: ONE transpose back for `np` vectors; `out`
    /// receives the np local coefficient vectors one after another.
    pub fn finish_to_reciprocal_batch<K: Collective>(
        &mut self,
        comm: &K,
        out: &mut [C],
        np: usize,
    ) -> Result<(), DistFftError> {
        let ll = self.scatter.len();
        if out.len() != ll * np || self.send.len() < 2 * self.block * self.ranks * np {
            return Err(FftError::IncompatibleArray.into());
        }
        let mark = self.profile.start();
        if self.block > 0 {
            dft_numerics::stage_clock::arrival::site(114);
        }
        self.all_to_all(comm, np)?;
        if self.block > 0 {
            self.profile.finish(
                mark,
                13,
                self.block * self.ranks * np,
                16 * self.block * self.ranks * np,
                None,
            );
        }
        let (nx, mp) = (self.n[0], self.max_planes);
        for t in 0..np {
            let mark = self.profile.start();
            for q in 0..self.ranks {
                let x0 = self.plane_starts[q];
                let qb = self.base(q, t, np);
                for sp in 0..self.local_sticks {
                    let base = 2 * (qb + sp * mp);
                    for xi in 0..self.plane_starts[q + 1] - x0 {
                        self.sticks[sp * nx + x0 + xi] =
                            C::new(self.recv[base + 2 * xi], self.recv[base + 2 * xi + 1]);
                    }
                }
            }
            self.profile
                .finish(mark, 14, self.sticks.len(), 32 * self.sticks.len(), None);
            if !self.sticks.is_empty() {
                let mark = self.profile.start();
                self.lanes.run(0, false, &mut self.sticks)?;
                self.profile
                    .finish(mark, 15, self.sticks.len(), 32 * self.sticks.len(), Some(0));
            }
            let mark = self.profile.start();
            for (o, &i) in out[t * ll..(t + 1) * ll].iter_mut().zip(&self.scatter) {
                *o = self.sticks[i] / self.total as f64;
            }
            self.profile
                .finish(mark, 16, self.scatter.len(), 32 * self.scatter.len(), None);
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
        let _profile_scope = self.profile.apply_scope();
        if potential.len() != self.total {
            return Err(FftError::IncompatibleArray.into());
        }
        self.begin_to_real(comm, c)?;
        self.multiply_planes(potential, 1)?;
        self.finish_to_reciprocal(comm, out)
    }

    /// `apply_local` for `np` vectors stored one after another in `data`
    /// (overwritten by the result) with one transpose per direction; bit-identical
    /// to np calls of `apply_local`.
    pub fn apply_local_batch<K: Collective>(
        &mut self,
        comm: &K,
        data: &mut [C],
        potential: &[f64],
        np: usize,
    ) -> Result<(), DistFftError> {
        let _profile_scope = self.profile.apply_scope();
        if potential.len() != self.total {
            return Err(FftError::IncompatibleArray.into());
        }
        self.begin_to_real_batch(comm, data, np)?;
        self.multiply_planes(potential, np)?;
        self.finish_to_reciprocal_batch(comm, data, np)
    }

    fn multiply_planes(&mut self, potential: &[f64], np: usize) -> Result<(), DistFftError> {
        for t in 0..np {
            for i in 0..self.plane_count() {
                self.to_real_plane_in(t, np, i)?;
                self.multiply_slab(potential, i);
                self.to_reciprocal_plane_in(t, np, i)?;
            }
        }
        Ok(())
    }

    /// Multiply the current slab (owned plane i) by `potential` (full row-major grid).
    fn multiply_slab(&mut self, potential: &[f64], i: usize) {
        let mark = self.profile.start();
        let (ny, nz) = (self.n[1], self.n[2]);
        for y in 0..ny {
            let start = self.grid_start(i, y);
            for (a, &v) in self.slab[y * nz..(y + 1) * nz]
                .iter_mut()
                .zip(&potential[start..start + nz])
            {
                *a *= v;
            }
        }
        self.profile
            .finish(mark, 8, self.slab.len(), 40 * self.slab.len(), None);
    }

    /// MPI-KEPT-START: `apply_local_batch` for `np` vectors of which the first
    /// `kept` have their real-space planes supplied in `planes` (vector t, owned
    /// plane i at `planes[t * plane_points() + i * rows() * shape[2]..]`, the
    /// slab `plane()` that `to_real_plane_in` leaves for that vector), so that
    /// only the r -> G transform runs for them; the remaining `np - kept` vectors
    /// are read from `data[kept * local_len()..]` and transformed both ways. All
    /// np results are written to `data` with ONE transpose back. A kept plane
    /// equal to the slab the inverse transform would produce gives the bits of
    /// `apply_local_batch` (the transforms of one vector do not depend on the
    /// batch). Collectives: one transpose to r when `kept < np`, one back; every
    /// rank must call with the same np and kept.
    pub fn apply_local_batch_kept<K: Collective>(
        &mut self,
        comm: &K,
        data: &mut [C],
        potential: &[f64],
        np: usize,
        kept: usize,
        planes: &[C],
    ) -> Result<(), DistFftError> {
        let _profile_scope = self.profile.apply_scope();
        let ll = self.scatter.len();
        let pp = self.plane_points();
        if potential.len() != self.total
            || data.len() != ll * np
            || kept > np
            || planes.len() < kept * pp
        {
            return Err(FftError::IncompatibleArray.into());
        }
        let fresh = np - kept;
        if fresh > 0 {
            self.begin_to_real_batch(comm, &data[kept * ll..], fresh)?;
        }
        // Keeps the received fresh data (resize preserves the contents).
        self.reserve_batch(np);
        let line = self.n[1] * self.n[2];
        for t in 0..kept {
            for i in 0..self.plane_count() {
                let source = &planes[t * pp + i * line..t * pp + (i + 1) * line];
                if self.kept_fuse {
                    // Same scalar operations as copy then multiply_slab. Read the
                    // immutable kept plane directly, avoiding a slab read/write pass.
                    let mark = self.profile.start();
                    let nz = self.n[2];
                    for y in 0..self.n[1] {
                        let start = self.grid_start(i, y);
                        for ((dst, &src), &v) in self.slab[y * nz..(y + 1) * nz]
                            .iter_mut()
                            .zip(&source[y * nz..(y + 1) * nz])
                            .zip(&potential[start..start + nz])
                        {
                            let mut value = src;
                            value *= v;
                            *dst = value;
                        }
                    }
                    self.profile.finish(mark, 8, line, 40 * line, None);
                } else {
                    let mark = self.profile.start();
                    self.slab.copy_from_slice(source);
                    self.profile.finish(mark, 17, line, 32 * line, None);
                    self.multiply_slab(potential, i);
                }
                self.to_reciprocal_plane_in(t, np, i)?;
            }
        }
        for t in 0..fresh {
            for i in 0..self.plane_count() {
                self.to_real_plane_in(t, fresh, i)?;
                self.multiply_slab(potential, i);
                self.to_reciprocal_plane_in(kept + t, np, i)?;
            }
        }
        self.finish_to_reciprocal_batch(comm, data, np)
    }

    /// Grid points on this rank's planes (`plane_count * rows * shape[2]`): the
    /// length of a plane-local field ordered [plane][row][z].
    pub fn plane_points(&self) -> usize {
        self.plane_count() * self.n[1] * self.n[2]
    }

    /// Assemble a full row-major grid field from every rank's plane-local field
    /// (one allgatherv; the copy into the grid is a pure permutation, so the
    /// result does not depend on the rank count). `scratch` grows to the grid size
    /// on first use.
    pub fn gather_planes<K: Collective>(
        &self,
        comm: &K,
        local: &[f64],
        scratch: &mut Vec<f64>,
        full: &mut [f64],
    ) -> Result<(), DistFftError> {
        if comm.size() != self.ranks || comm.rank() != self.rank {
            return Err(MpiError::Shape("communicator differs from the plan").into());
        }
        if local.len() != self.plane_points() || full.len() != self.total {
            return Err(FftError::IncompatibleArray.into());
        }
        let (ny, nz) = (self.n[1], self.n[2]);
        let counts: Vec<usize> = self
            .plane_starts
            .windows(2)
            .map(|w| (w[1] - w[0]) * ny * nz)
            .collect();
        scratch.resize(self.total, 0.);
        comm.allgatherv_f64(local, &counts, scratch)?;
        // Plane x of the concatenation starts at x*ny*nz for every rank count.
        for x in 0..self.n[0] {
            for y in 0..ny {
                let start = if self.first == 0 {
                    (x * self.shape[1] + y) * nz
                } else {
                    (y * self.shape[1] + x) * nz
                };
                full[start..start + nz].copy_from_slice(&scratch[(x * ny + y) * nz..][..nz]);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use dft_mpi::simulate;
    use dft_numerics::fused::FusedOperator;
    use dft_numerics::{Fft3d, FftPlan3d, PairFusedPlan};
    use ndarray::Array3;

    #[test]
    fn kept_fusion_bits_counters_and_rank_independence() {
        use dft_numerics::real::fft_profile;
        use dft_numerics::switches::RunConfig;
        let mut fixtures = cases();
        // First axis y, ranks with no planes/sticks, and empty reciprocal support.
        fixtures.push(([3, 5, 7], vec![[0, 0, 0], [0, 1, 0], [0, 2, 0], [0, 3, 0]]));
        fixtures.push(([1, 2, 3], vec![]));
        for (shape, idx) in fixtures {
            let potential: Vec<_> = (0..shape.iter().product())
                .map(|i| if i % 5 == 0 { -0.0 } else { (i as f64).cos() })
                .collect();
            for np in [1, 3, 8] {
                let input = random_c(idx.len() * np, 509);
                let mut rank_one = Vec::new();
                for ranks in [1, 2, 3, 8] {
                    let results = simulate(ranks, |comm| {
                        let cfg = RunConfig::from_pairs(&[
                            ("fft_kept_fuse", "off"),
                            ("fft_profile", "on"),
                            ("phase_profile", "on"),
                            ("stage_timers", "off"),
                        ])
                        .unwrap();
                        let _scope = cfg.enter();
                        let mut p =
                            DistributedPairFft::new(shape, &idx, comm.rank(), ranks).unwrap();
                        let local: Vec<_> = (0..np)
                            .flat_map(|t| {
                                p.local_positions()
                                    .iter()
                                    .map(|&g| input[t * idx.len() + g])
                                    .collect::<Vec<_>>()
                            })
                            .collect();
                        let mut plain = local.clone();
                        p.apply_local_batch(comm, &mut plain, &potential, np)
                            .unwrap();
                        p.begin_to_real_batch(comm, &local, np).unwrap();
                        let pp = p.plane_points();
                        let line = p.rows() * shape[2];
                        let mut planes = vec![C::default(); np * pp];
                        for t in 0..np {
                            for i in 0..p.plane_count() {
                                p.to_real_plane_in(t, np, i).unwrap();
                                planes[t * pp + i * line..t * pp + (i + 1) * line]
                                    .copy_from_slice(p.plane());
                            }
                        }
                        let bits = |xs: &[C]| {
                            xs.iter()
                                .map(|v| (v.re.to_bits(), v.im.to_bits()))
                                .collect::<Vec<_>>()
                        };
                        let plane_bits = bits(&planes);
                        for kept in [0, 1, np] {
                            let mut off_counts = fft_profile::Snapshot::default();
                            for enabled in [false, true] {
                                p.kept_fuse = enabled;
                                let pointers = (p.slab.as_ptr(), p.send.as_ptr(), p.recv.as_ptr());
                                let mut out = local.clone();
                                let before = fft_profile::snapshot();
                                p.apply_local_batch_kept(
                                    comm, &mut out, &potential, np, kept, &planes,
                                )
                                .unwrap();
                                let counts = fft_profile::snapshot().since(before);
                                assert_eq!(bits(&out), bits(&plain));
                                assert_eq!(bits(&planes), plane_bits);
                                assert_eq!(
                                    pointers,
                                    (p.slab.as_ptr(), p.send.as_ptr(), p.recv.as_ptr())
                                );
                                if enabled {
                                    assert_eq!(counts.0[17], [0; 5]);
                                    for stage in 0..17 {
                                        assert_eq!(counts.0[stage], off_counts.0[stage]);
                                    }
                                } else {
                                    assert_eq!(counts.0[17][1], (kept * p.plane_count()) as u64);
                                    assert_eq!(counts.0[17][3], (32 * kept * pp) as u64);
                                    off_counts = counts;
                                }
                            }
                        }
                        (p.local_positions.clone(), plain)
                    });
                    let mut global = vec![C::default(); input.len()];
                    for (positions, values) in results {
                        for t in 0..np {
                            for (l, &g) in positions.iter().enumerate() {
                                global[t * idx.len() + g] = values[t * positions.len() + l];
                            }
                        }
                    }
                    let bits: Vec<_> = global
                        .iter()
                        .map(|v| (v.re.to_bits(), v.im.to_bits()))
                        .collect();
                    if ranks == 1 {
                        rank_one = bits;
                    } else {
                        assert_eq!(bits, rank_one);
                    }
                }
            }
        }
    }

    #[test]
    fn lane_fft_instrumentation_preserves_bits_and_counts() {
        use dft_numerics::real::fft_profile;
        use dft_numerics::switches::RunConfig;
        let mut fixtures = cases();
        fixtures.push(([2, 3, 4], vec![]));
        for (shape, idx) in fixtures {
            let pot: Vec<_> = (0..shape.iter().product())
                .map(|i| (i as f64).cos())
                .collect();
            for ranks in [1, 2, 3, 8] {
                simulate(ranks, |comm| {
                    let mut reference = Vec::new();
                    for (arm, timers) in [("off", "off"), ("on", "off"), ("on", "on")] {
                        let cfg = RunConfig::from_pairs(&[
                            ("fft_profile", arm),
                            ("phase_profile", "on"),
                            ("stage_timers", timers),
                        ])
                        .unwrap();
                        let _scope = cfg.enter();
                        let mut p =
                            DistributedPairFft::new(shape, &idx, comm.rank(), ranks).unwrap();
                        let np = 3;
                        let mut data = random_c(p.local_len() * np, 103);
                        let before = fft_profile::snapshot();
                        p.apply_local_batch(comm, &mut data, &pot, np).unwrap();
                        let delta = fft_profile::snapshot().since(before);
                        if arm == "off" {
                            assert_eq!(delta, fft_profile::Snapshot::default());
                            reference = data;
                        } else {
                            for (a, b) in data.iter().zip(&reference) {
                                assert_eq!(
                                    (a.re.to_bits(), a.im.to_bits()),
                                    (b.re.to_bits(), b.im.to_bits())
                                );
                            }
                            assert_eq!(delta.0[0][1], np as u64);
                            assert_eq!(delta.0[3][1], u64::from(p.block > 0));
                            assert_eq!(delta.0[13][1], u64::from(p.block > 0));
                            assert_eq!(delta.0[8][1], (np * p.plane_count()) as u64);
                            assert_eq!(delta.0[8][2], (np * p.plane_points()) as u64);
                            assert_eq!(delta.0[8][3], (40 * np * p.plane_points()) as u64);
                            if timers == "off" {
                                assert!(delta.0.iter().all(|row| row[0] == 0));
                            }
                        }
                    }
                });
            }
        }
    }

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
                p.to_real_plane(i).unwrap();
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
    fn mirror_partition_keeps_gamma_pairs_on_one_rank() {
        // BAANDI_M2_STICKS=mirror: exact partition for every case; on the
        // Gamma-symmetric spheres +G and -G always share an owner, while the greedy
        // rule splits some pairs for two or more ranks.
        for (case, (shape, idx)) in cases().into_iter().enumerate() {
            let mirror_of = |p: &[usize; 3]| {
                [
                    (shape[0] - p[0]) % shape[0],
                    (shape[1] - p[1]) % shape[1],
                    (shape[2] - p[2]) % shape[2],
                ]
            };
            let symmetric = idx.iter().all(|p| idx.contains(&mirror_of(p)));
            for ranks in 1..=8 {
                for mirror in [false, true] {
                    set_mirror_sticks_override(Some(mirror));
                    let owners = position_owners(shape, &idx, ranks).unwrap();
                    let mut seen = vec![0; idx.len()];
                    for r in 0..ranks {
                        let p = DistributedPairFft::new(shape, &idx, r, ranks).unwrap();
                        for &g in p.local_positions() {
                            seen[g] += 1;
                            assert_eq!(owners[g], r);
                        }
                    }
                    set_mirror_sticks_override(None);
                    assert!(seen.iter().all(|&s| s == 1), "case {case} ranks {ranks}");
                    if !symmetric {
                        continue;
                    }
                    let split = idx
                        .iter()
                        .enumerate()
                        .filter(|(g, p)| {
                            let m = idx.iter().position(|q| *q == mirror_of(p)).unwrap();
                            owners[*g] != owners[m]
                        })
                        .count();
                    if mirror {
                        assert_eq!(split, 0, "case {case} ranks {ranks}");
                    } else if ranks >= 2 {
                        assert!(split > 0, "case {case} ranks {ranks}");
                    }
                }
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

    /// `BAANDI_M2_FFT=mkl` lanes: within 1e-14 relative of the rustfft lanes on
    /// every rank count (MKL builds), bit-identical on pure builds (no vendor lanes).
    #[test]
    fn vendor_lanes_match_rust_lanes() {
        for (shape, idx) in cases() {
            let n = shape.iter().product::<usize>();
            let c = random_c(idx.len(), 21);
            let mut s = 17;
            let pot: Vec<f64> = (0..n).map(|_| 1.0 + lcg(&mut s)).collect();
            for ranks in [1usize, 2, 3, 5] {
                let parts = simulate(ranks, |comm| {
                    let mut out = Vec::new();
                    for vendor in [false, true] {
                        let mut p =
                            DistributedPairFft::with_lanes(shape, &idx, comm.rank(), ranks, vendor)
                                .unwrap();
                        assert_eq!(p.vendor_lanes(), vendor && dft_numerics::MKL_BACKEND);
                        let mut data: Vec<C> = p.local_positions().iter().map(|&g| c[g]).collect();
                        let np = 1;
                        p.apply_local_batch(comm, &mut data, &pot, np).unwrap();
                        out.push(data);
                    }
                    out
                });
                for part in &parts {
                    let e = rel(&part[1], &part[0]);
                    assert!(e < 1e-14, "{shape:?} P={ranks}: vendor lanes {e:e}");
                    if !dft_numerics::MKL_BACKEND {
                        assert_eq!(part[0], part[1], "pure build: identical lanes");
                    }
                }
            }
        }
    }

    #[test]
    fn batched_apply_equals_single_calls_bitwise_and_gather_assembles_the_grid() {
        for (shape, idx) in cases() {
            let n = shape.iter().product::<usize>();
            let mut s = 5;
            let pot: Vec<f64> = (0..n).map(|_| 1.0 + lcg(&mut s)).collect();
            let field: Vec<f64> = (0..n).map(|_| lcg(&mut s)).collect();
            let np = 3;
            let vecs: Vec<Vec<C>> = (0..np)
                .map(|t| random_c(idx.len(), 40 + t as u64))
                .collect();
            for ranks in 1..=5 {
                let res = simulate(ranks, |comm| {
                    let mut p = DistributedPairFft::new(shape, &idx, comm.rank(), ranks).unwrap();
                    let ll = p.local_len();
                    let mut single = vec![C::default(); ll * np];
                    for t in 0..np {
                        let local: Vec<C> =
                            p.local_positions().iter().map(|&g| vecs[t][g]).collect();
                        p.apply_local(comm, &local, &pot, &mut single[t * ll..(t + 1) * ll])
                            .unwrap();
                    }
                    let vecs = &vecs;
                    let mut batch: Vec<C> = (0..np)
                        .flat_map(|t| p.local_positions().iter().map(move |&g| vecs[t][g]))
                        .collect();
                    p.apply_local_batch(comm, &mut batch, &pot, np).unwrap();
                    // Plane-local copy of a known field, gathered back.
                    let mut local = Vec::with_capacity(p.plane_points());
                    for i in 0..p.plane_count() {
                        for y in 0..p.rows() {
                            let st = p.grid_start(i, y);
                            local.extend_from_slice(&field[st..st + shape[2]]);
                        }
                    }
                    let mut scratch = Vec::new();
                    let mut full = vec![0.; n];
                    p.gather_planes(comm, &local, &mut scratch, &mut full)
                        .unwrap();
                    (single == batch, full == field)
                });
                assert!(
                    res.iter().all(|r| r.0),
                    "{shape:?} P={ranks}: batch differs"
                );
                assert!(
                    res.iter().all(|r| r.1),
                    "{shape:?} P={ranks}: gather differs"
                );
            }
        }
    }

    /// MPI-KEPT-START: with the planes of the first `kept` vectors taken from
    /// their own inverse transforms, `apply_local_batch_kept` gives the bits of
    /// `apply_local_batch` for every kept count, 1 to 5 ranks.
    #[test]
    fn kept_planes_batch_equals_plain_batch_bitwise() {
        for (shape, idx) in cases() {
            let n = shape.iter().product::<usize>();
            let mut s = 11;
            let pot: Vec<f64> = (0..n).map(|_| lcg(&mut s) - 0.5).collect();
            let np = 3;
            let vecs: Vec<Vec<C>> = (0..np)
                .map(|t| random_c(idx.len(), 70 + t as u64))
                .collect();
            for ranks in 1..=5 {
                let res = simulate(ranks, |comm| {
                    let mut p = DistributedPairFft::new(shape, &idx, comm.rank(), ranks).unwrap();
                    let vecs = &vecs;
                    let input: Vec<C> = (0..np)
                        .flat_map(|t| p.local_positions().iter().map(move |&g| vecs[t][g]))
                        .collect();
                    let mut plain = input.clone();
                    p.apply_local_batch(comm, &mut plain, &pot, np).unwrap();
                    // Planes of every vector from the inverse transform alone.
                    let pp = p.plane_points();
                    let mut planes = vec![C::default(); np * pp];
                    p.begin_to_real_batch(comm, &input, np).unwrap();
                    let line = p.rows() * shape[2];
                    for t in 0..np {
                        for i in 0..p.plane_count() {
                            p.to_real_plane_in(t, np, i).unwrap();
                            planes[t * pp + i * line..t * pp + (i + 1) * line]
                                .copy_from_slice(p.plane());
                        }
                    }
                    (0..=np)
                        .map(|kept| {
                            let mut out = input.clone();
                            p.apply_local_batch_kept(comm, &mut out, &pot, np, kept, &planes)
                                .unwrap();
                            out == plain
                        })
                        .collect::<Vec<_>>()
                });
                for r in &res {
                    assert!(r.iter().all(|&b| b), "{shape:?} P={ranks}: {r:?}");
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
