//! Stick-distributed real Gamma Hamiltonian (MPI stage M2b).
//!
//! Placement: dft-physics, not dft-solvers, because the operator needs the
//! `Hamiltonian`/`NonlocalOperator` snapshot types and `pack_projector`, and
//! dft-physics already depends on dft-solvers (which holds `DistributedPairFft`,
//! `dft_mpi::Collective` comes in through dft-mpi). dft-solvers cannot depend on
//! dft-physics without a cycle. The M1 Davidson can call `apply_block` through a
//! small adapter in place of `apply_distributed` (not wired in this stage).
//!
//! Row partition. The ORC vector has n = 2*nh - 1 real rows: row 0 is G=0, the
//! representative i >= 1 owns rows 2i-1, 2i. The M2a stick partition owns full-sphere
//! complex positions. Representative i is owned by the rank that owns position
//! `plus[i]` (G=0 by the owner of position 0). A rank's local rows are its
//! representatives' ORC rows in ascending ORC row index (`local_rows`), so the
//! permutation "ORC row -> (rank, local row)" is a pure relabelling: every algebra
//! step (Gram sums, rotations, norms) is row-wise or a sum over rows and is
//! independent of it. For one rank the permutation is the identity.
//!
//! Mirror data. A pair (a, b) of real columns is packed into one complex vector
//! with entries at `plus[i]` and `minus[i]`; the two positions can sit on different
//! sticks and ranks. Let L(r -> q) be the representatives with plus on r and minus
//! on q != r (ascending i). Per apply: (1) one alltoall sends the ORC rows of
//! L(r -> q) for all columns from r to q (so q can pack its minus positions);
//! (2) per pair of columns the distributed pair FFT runs (2 alltoalls: G->r, r->G),
//! multiplying by the rank's owned planes of the potential; (3) one alltoall sends
//! h[minus] of L(r -> q) for all pairs from q back to r, which unpacks. Pack and
//! unpack use exactly the formulas of `PairMap::pack/unpack`.
//!
//! Nonlocal. Projector rows are partitioned like the basis rows. P^T X is a local
//! GEMM (count x columns), summed over ranks with one `tree_sum_f64`, the per-atom
//! coupling is replicated, P W is a local GEMM. The kinetic term is local.
//!
//! Collective calls per apply with m columns, P = ceil(m/2) pairs (M2d batched
//! transposes, default): alltoall = 2 + 2 (one FFT transpose per direction for all
//! pairs of the block, the 2 mirror exchanges are skipped when no representative
//! is split, e.g. one rank), tree_sum = 1 (0 without projectors). Unbatched
//! (`set_batched(false)`, the M2b/M2c code): alltoall = 2P + 2. Both are
//! bit-identical (every pair keeps its arithmetic; only data movement differs).
//! See `collectives_per_apply`. All counts are rank independent (no deadlock).
//!
//! M2d additions: `band_densities` (occupation-weighted |psi|^2 of ORC columns on
//! this rank's FFT planes from the same pair packing and one batched G->r
//! transpose), `gather_planes` (plane-local field to the full grid, one
//! allgatherv) and `nonlocal_band_energies` (x^T P D P^T x per column with one
//! tree_sum). Counters and stage timers: `m2_counters`.
//!
//! Numerics. One rank reproduces `RealHamiltonian::apply_block` to rounding; it is
//! bit-identical when the serial and distributed FFT paths agree bitwise (tested,
//! reported by the gate test). Nonlocal GEMMs: the pure-Rust tiled path by default
//! (`BAANDI_MKL_NONLOCAL` is not honoured); `BAANDI_M2_DGEMM=on` takes both projector
//! products through `multiply_slices_wide_into` like the serial MKL path (per-panel
//! dgemm on MKL builds, a rounding change; the same product on pure builds). Inputs are not scanned for NaN/inf: an early error on one rank only
//! would desynchronise the collectives; non-finite values propagate to the output.
//! No allocation after construction (conversion helpers `extract_local` and
//! `scatter_local` allocate and are for the library boundary); no unsafe; no logging.
use crate::hamiltonian::Hamiltonian;
use crate::real_gamma::{RealGammaError, pack_projector};
use dft_basis::{HalfSphere, PlaneWaveBasis};
use dft_mpi::{Collective, MpiError};
use dft_numerics::linalg::LinalgError;
use dft_numerics::{FftError, GemmSite, RealGemmWorkspace, RealMatrix, RealOp, pool::WorkerPool};
use dft_solvers::distributed_fft::{DistFftError, DistributedPairFft, position_owners};
use num_complex::Complex64 as C;
use std::f64::consts::SQRT_2;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::time::Instant;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum DistributedRealError {
    #[error("incompatible basis, rank layout or block dimensions")]
    Shape,
    #[error(transparent)]
    Real(#[from] RealGammaError),
    #[error(transparent)]
    Half(#[from] dft_basis::HalfSphereError),
    #[error(transparent)]
    DistFft(#[from] DistFftError),
    #[error(transparent)]
    Fft(#[from] FftError),
    #[error(transparent)]
    Algebra(#[from] LinalgError),
    #[error(transparent)]
    Mpi(#[from] MpiError),
}

/// Where the packed value at one local FFT position comes from.
#[derive(Clone, Copy)]
enum Src {
    /// G=0: local row slot.
    Zero(usize),
    /// plus position: first of two local row slots.
    Plus(usize),
    MinusLocal(usize),
    /// Item id in the received mirror rows.
    MinusRemote(usize),
}
#[derive(Clone, Copy)]
enum MinusAt {
    None,
    /// Local FFT position.
    Local(usize),
    /// Item id in the received mirror h values.
    Remote(usize),
}
struct Rep {
    slot: usize,
    plus_pos: usize,
    minus: MinusAt,
}

struct Coupling {
    start: usize,
    count: usize,
    values: Vec<f64>,
}

pub struct DistributedRealHamiltonian {
    rank: usize,
    ranks: usize,
    n_global: usize,
    rows: Vec<usize>,
    capacity: usize,
    kinetic: Vec<f64>,
    potential: Vec<f64>,
    fft: DistributedPairFft,
    src: Vec<Src>,
    reps: Vec<Rep>,
    /// (first row slot, item id) of rows sent to another rank for packing.
    out_items: Vec<(usize, usize)>,
    /// (local FFT position of minus, item id) of h values sent back for unpacking.
    in_items: Vec<(usize, usize)>,
    max_items: usize,
    /// One transpose per direction for all pairs of a block (default true).
    batched: bool,
    hbuf: Vec<C>,
    gather_scratch: Vec<f64>,
    xsend: Vec<f64>,
    xrecv: Vec<f64>,
    hsend: Vec<f64>,
    hrecv: Vec<f64>,
    projectors: RealMatrix,
    coupling: Vec<Coupling>,
    count: usize,
    projections: RealMatrix,
    reduced: RealMatrix,
    weights: RealMatrix,
    image: RealMatrix,
    gemm: RealGemmWorkspace,
    /// `BAANDI_M2_DGEMM=on`: projector products through the per-panel dgemm route.
    dgemm: bool,
    source_identity: Arc<()>,
    pool: Arc<WorkerPool>,
}

impl DistributedRealHamiltonian {
    /// Snapshot of the complex Gamma H for rank `rank` of `ranks`.
    pub fn new(
        source: &Hamiltonian,
        basis: &PlaneWaveBasis,
        capacity: usize,
        rank: usize,
        ranks: usize,
        pool: Arc<WorkerPool>,
    ) -> Result<Self, DistributedRealError> {
        if source.basis_momenta() != basis.momenta() || rank >= ranks {
            return Err(DistributedRealError::Shape);
        }
        let grid = source.fft_grid();
        let half = HalfSphere::new(basis, grid)?;
        let (plus, minus) = (half.plus(), half.minus());
        let n = half.full_len();
        let nh = plus.len();
        let dims = grid.dims();
        let owners = position_owners(dims, grid.g_to_fft(), ranks)?;
        let fft = DistributedPairFft::new(dims, grid.g_to_fft(), rank, ranks)?;

        // Local rows: representatives owned by this rank, ascending ORC row.
        let mut slot_of_rep = vec![usize::MAX; nh];
        let mut rows = Vec::new();
        for i in 0..nh {
            if owners[plus[i] as usize] == rank {
                slot_of_rep[i] = rows.len();
                if i == 0 {
                    rows.push(0);
                } else {
                    rows.extend([2 * i - 1, 2 * i]);
                }
            }
        }
        let n_loc = rows.len();
        // Mirror links.
        let mut link = vec![Vec::new(); ranks * ranks];
        for i in 1..nh {
            let (a, b) = (owners[plus[i] as usize], owners[minus[i] as usize]);
            if a != b {
                link[a * ranks + b].push(i);
            }
        }
        let max_items = link.iter().map(Vec::len).max().unwrap_or(0);
        let mut pos_to_local = vec![usize::MAX; n];
        for (k, &g) in fft.local_positions().iter().enumerate() {
            pos_to_local[g] = k;
        }
        // (representative, is_plus) of every full-sphere position.
        let mut rep_of_pos = vec![(0usize, true); n];
        for i in 1..nh {
            rep_of_pos[plus[i] as usize] = (i, true);
            rep_of_pos[minus[i] as usize] = (i, false);
        }
        let item = |a: usize, b: usize, i: usize| {
            a * max_items + link[a * ranks + b].binary_search(&i).unwrap()
        };
        let src = fft
            .local_positions()
            .iter()
            .map(|&g| {
                let (i, is_plus) = rep_of_pos[g];
                if i == 0 {
                    Src::Zero(slot_of_rep[0])
                } else if is_plus {
                    Src::Plus(slot_of_rep[i])
                } else if owners[plus[i] as usize] == rank {
                    Src::MinusLocal(slot_of_rep[i])
                } else {
                    Src::MinusRemote(item(owners[plus[i] as usize], rank, i))
                }
            })
            .collect();
        let mut reps = Vec::new();
        for i in 0..nh {
            if slot_of_rep[i] == usize::MAX {
                continue;
            }
            let minus_at = if i == 0 {
                MinusAt::None
            } else if owners[minus[i] as usize] == rank {
                MinusAt::Local(pos_to_local[minus[i] as usize])
            } else {
                // Block of the returned h values is indexed by the minus owner.
                let b = owners[minus[i] as usize];
                MinusAt::Remote(b * max_items + link[rank * ranks + b].binary_search(&i).unwrap())
            };
            reps.push(Rep {
                slot: slot_of_rep[i],
                plus_pos: pos_to_local[plus[i] as usize],
                minus: minus_at,
            });
        }
        let mut out_items = Vec::new();
        let mut in_items = Vec::new();
        for q in 0..ranks {
            for (j, &i) in link[rank * ranks + q].iter().enumerate() {
                out_items.push((slot_of_rep[i], q * max_items + j));
            }
            for (j, &i) in link[q * ranks + rank].iter().enumerate() {
                in_items.push((pos_to_local[minus[i] as usize], q * max_items + j));
            }
        }

        // Kinetic diagonal in ORC order (as RealHamiltonian), restricted to local rows.
        let mut kinetic_full = vec![0.; n];
        kinetic_full[0] = source.kinetic_diagonal()[0];
        for (i, &p) in plus.iter().enumerate().skip(1) {
            kinetic_full[2 * i - 1] = source.kinetic_diagonal()[p as usize];
            kinetic_full[2 * i] = kinetic_full[2 * i - 1];
        }
        let kinetic = rows.iter().map(|&r| kinetic_full[r]).collect();

        // Projector rows restricted to local rows.
        let count: usize = source
            .nonlocal()
            .atoms()
            .iter()
            .map(|a| a.columns().len())
            .sum();
        let mut projectors = RealMatrix::zeros(n_loc, count)?;
        let mut coupling = Vec::new();
        let mut full = vec![0.; n];
        let mut start = 0;
        for atom in source.nonlocal().atoms() {
            for (j, beta) in atom.columns().iter().enumerate() {
                pack_projector(&half, beta, &mut full)?;
                let col =
                    &mut projectors.as_mut_slice()[(start + j) * n_loc..(start + j + 1) * n_loc];
                for (c, &r) in col.iter_mut().zip(&rows) {
                    *c = full[r];
                }
            }
            coupling.push(Coupling {
                start,
                count: atom.columns().len(),
                values: atom.coupling_hartree().to_vec(),
            });
            start += atom.columns().len();
        }
        let pairs = capacity.div_ceil(2);
        let local_len = fft.local_len();
        Ok(Self {
            rank,
            ranks,
            n_global: n,
            capacity,
            kinetic,
            potential: source.local_potential_values().to_vec(),
            src,
            reps,
            out_items,
            in_items,
            max_items,
            batched: true,
            hbuf: vec![C::default(); pairs * local_len],
            gather_scratch: Vec::new(),
            xsend: vec![0.; ranks * max_items * 2 * capacity],
            xrecv: vec![0.; ranks * max_items * 2 * capacity],
            hsend: vec![0.; ranks * max_items * 2 * pairs],
            hrecv: vec![0.; ranks * max_items * 2 * pairs],
            projections: RealMatrix::zeros(count, capacity)?,
            reduced: RealMatrix::zeros(count, capacity)?,
            weights: RealMatrix::zeros(count, capacity)?,
            image: RealMatrix::zeros(n_loc, capacity)?,
            gemm: RealGemmWorkspace::with_pool(n_loc, count.max(capacity), pool.clone())?,
            dgemm: dft_numerics::switches::on(dft_numerics::switches::Id::M2_DGEMM),
            source_identity: source.operator_identity(),
            pool,
            projectors,
            coupling,
            count,
            rows,
            fft,
        })
    }

    /// Refresh the local potential from the updated complex Hamiltonian of the same
    /// operator (same basis, projectors and kinetic term), as
    /// `RealHamiltonian::update_from_complex` does; false when `source` is another
    /// operator (the caller rebuilds).
    pub fn update_from_complex(&mut self, source: &Hamiltonian, basis: &PlaneWaveBasis) -> bool {
        if !Arc::ptr_eq(&self.source_identity, &source.operator_identity())
            || source.basis_momenta() != basis.momenta()
            || basis.k_fractional() != [0.; 3]
            || self.potential.len() != source.local_potential_values().len()
        {
            return false;
        }
        self.potential
            .copy_from_slice(source.local_potential_values());
        true
    }

    /// ORC row indices owned by this rank, ascending; local row l is ORC row `local_rows()[l]`.
    pub fn local_rows(&self) -> &[usize] {
        &self.rows
    }
    pub fn local_len(&self) -> usize {
        self.rows.len()
    }
    pub fn global_len(&self) -> usize {
        self.n_global
    }

    /// (alltoall calls, tree_sum calls) of one `apply_block` with `columns` columns.
    pub fn collectives_per_apply(&self, columns: usize) -> (usize, usize) {
        if columns == 0 {
            return (0, 0);
        }
        let pairs = columns.div_ceil(2);
        let transposes = if self.batched { 2 } else { 2 * pairs };
        (
            transposes + if self.max_items > 0 { 2 } else { 0 },
            usize::from(self.count > 0),
        )
    }

    /// Choose batched (one transpose per direction per block, default) or per-pair
    /// FFT transposes. Results are bit-identical; only the collective count differs.
    pub fn set_batched(&mut self, batched: bool) {
        self.batched = batched;
    }
    pub fn batched(&self) -> bool {
        self.batched
    }
    /// Test hook: force the dgemm route of the projector products (`BAANDI_M2_DGEMM`).
    pub fn set_dgemm(&mut self, on: bool) {
        self.dgemm = on;
    }

    /// Grid points of this rank's FFT planes (length of a plane-local field).
    pub fn plane_points(&self) -> usize {
        self.fft.plane_points()
    }
    /// FFT grid shape of the operator (row-major).
    pub fn grid_shape(&self) -> [usize; 3] {
        self.fft.shape()
    }
    /// Row-major grid offset of every plane-local point (allocates; setup only).
    pub fn plane_grid_offsets(&self) -> Vec<usize> {
        let nz = self.fft.shape()[2];
        let mut out = Vec::with_capacity(self.fft.plane_points());
        for i in 0..self.fft.plane_count() {
            for y in 0..self.fft.rows() {
                let s = self.fft.grid_start(i, y);
                out.extend(s..s + nz);
            }
        }
        out
    }

    /// Full row-major grid from every rank's plane-local field (one allgatherv).
    pub fn gather_planes<K: Collective>(
        &mut self,
        comm: &K,
        local: &[f64],
        full: &mut [f64],
    ) -> Result<(), DistributedRealError> {
        let start = Instant::now();
        self.fft
            .gather_planes(comm, local, &mut self.gather_scratch, full)?;
        stage_add(Stage::Gather, start);
        Ok(())
    }

    /// Exchange the ORC rows of split representatives so that the minus owner can
    /// pack (one alltoall; skipped when nothing is split).
    fn exchange_mirror_rows<K: Collective>(
        &mut self,
        comm: &K,
        xs: &[f64],
        m: usize,
    ) -> Result<(), DistributedRealError> {
        let mi = self.max_items;
        if mi == 0 {
            return Ok(());
        }
        let n_loc = self.rows.len();
        let ranks = self.ranks;
        let blk = mi * 2 * m;
        self.xsend[..ranks * blk].fill(0.);
        for &(slot, id) in &self.out_items {
            for c in 0..m {
                let off = (id * m + c) * 2;
                self.xsend[off] = xs[c * n_loc + slot];
                self.xsend[off + 1] = xs[c * n_loc + slot + 1];
            }
        }
        comm.alltoall_f64(
            &self.xsend[..ranks * blk],
            blk,
            &mut self.xrecv[..ranks * blk],
        )?;
        Ok(())
    }

    /// Pack column pair t (columns 2t, 2t+1 < m) of the local block into the
    /// complex local FFT vector `out`, exactly as `PairMap::pack`.
    fn pack_pair(
        src: &[Src],
        xrecv: &[f64],
        xs: &[f64],
        n_loc: usize,
        m: usize,
        t: usize,
        out: &mut [C],
    ) {
        let (ca, cb) = (2 * t, (2 * t + 1 < m).then_some(2 * t + 1));
        let get =
            |col: usize, slot: usize| C::new(xs[col * n_loc + slot], xs[col * n_loc + slot + 1]);
        let recv = |col: usize, id: usize| {
            let off = (id * m + col) * 2;
            C::new(xrecv[off], xrecv[off + 1])
        };
        for (phi, &s) in out.iter_mut().zip(src) {
            *phi = match s {
                Src::Zero(slot) => C::new(
                    xs[ca * n_loc + slot],
                    cb.map_or(0., |cb| xs[cb * n_loc + slot]),
                ),
                Src::Plus(slot) | Src::MinusLocal(slot) => {
                    let a = get(ca, slot) / SQRT_2;
                    let b = cb.map_or(C::default(), |cb| get(cb, slot)) / SQRT_2;
                    if matches!(s, Src::Plus(_)) {
                        a + C::new(-b.im, b.re)
                    } else {
                        a.conj() + C::new(b.im, b.re)
                    }
                }
                Src::MinusRemote(id) => {
                    let a = recv(ca, id) / SQRT_2;
                    let b = cb.map_or(C::default(), |cb| recv(cb, id)) / SQRT_2;
                    a.conj() + C::new(b.im, b.re)
                }
            };
        }
    }

    /// Per-column occupation-weighted densities on this rank's FFT planes (MPI
    /// stage M2d): `out[b * plane_points() + p] = weights[b] * psi_b(r_p)^2 / volume`
    /// for the first m = x.cols() columns, where psi_b is the real-space orbital of
    /// ORC column b (the convention of `PairDensityOperator::accumulate_density_pair`).
    /// Plane-local points are ordered [plane][row][z]. Collectives: one mirror
    /// alltoall (when split) and ONE batched G->r transpose. Every rank must call
    /// with the same m and weights.
    pub fn band_densities<K: Collective>(
        &mut self,
        comm: &K,
        x: &RealMatrix,
        weights: &[f64],
        volume: f64,
        out: &mut [f64],
    ) -> Result<(), DistributedRealError> {
        let start = Instant::now();
        let n_loc = self.rows.len();
        let m = x.cols();
        let points = self.fft.plane_points();
        if comm.size() != self.ranks
            || comm.rank() != self.rank
            || x.rows() != n_loc
            || m > self.capacity
            || weights.len() != m
            || out.len() != m * points
            || !(volume.is_finite() && volume > 0.)
        {
            return Err(DistributedRealError::Shape);
        }
        if m == 0 {
            return Ok(());
        }
        let pairs = m.div_ceil(2);
        let ll = self.fft.local_len();
        let xs = x.as_slice();
        self.exchange_mirror_rows(comm, xs, m)?;
        for t in 0..pairs {
            Self::pack_pair(
                &self.src,
                &self.xrecv,
                xs,
                n_loc,
                m,
                t,
                &mut self.hbuf[t * ll..(t + 1) * ll],
            );
        }
        self.fft
            .begin_to_real_batch(comm, &self.hbuf[..pairs * ll], pairs)?;
        let (rows, nz) = (self.fft.rows(), self.fft.shape()[2]);
        let line = rows * nz;
        for t in 0..pairs {
            let (ca, cb) = (2 * t, (2 * t + 1 < m).then_some(2 * t + 1));
            for i in 0..self.fft.plane_count() {
                self.fft.to_real_plane_in(t, pairs, i)?;
                let slab = self.fft.plane();
                let base = i * line;
                let wa = weights[ca];
                for (o, v) in out[ca * points + base..ca * points + base + line]
                    .iter_mut()
                    .zip(slab)
                {
                    *o = wa * v.re * v.re / volume;
                }
                if let Some(cb) = cb {
                    let wb = weights[cb];
                    for (o, v) in out[cb * points + base..cb * points + base + line]
                        .iter_mut()
                        .zip(slab)
                    {
                        *o = wb * v.im * v.im / volume;
                    }
                }
            }
        }
        stage_add(Stage::Density, start);
        COUNTS[Count::Densities as usize].fetch_add(1, Ordering::Relaxed);
        Ok(())
    }

    /// Per-column nonlocal energies x_b^T P D P^T x_b (Ha for normalized columns)
    /// of the first m = x.cols() columns: local GEMM, one tree_sum, replicated
    /// coupling. Every rank receives identical values. No collective without
    /// projectors.
    pub fn nonlocal_band_energies<K: Collective>(
        &mut self,
        comm: &K,
        x: &RealMatrix,
        out: &mut [f64],
    ) -> Result<(), DistributedRealError> {
        let start = Instant::now();
        let n_loc = self.rows.len();
        let m = x.cols();
        if comm.size() != self.ranks
            || comm.rank() != self.rank
            || x.rows() != n_loc
            || m > self.capacity
            || out.len() != m
        {
            return Err(DistributedRealError::Shape);
        }
        out.fill(0.);
        let count = self.count;
        if count == 0 || m == 0 {
            return Ok(());
        }
        self.projections_reduced(comm, x)?;
        for (band, o) in out.iter_mut().enumerate() {
            let mut sum = 0.;
            for atom in &self.coupling {
                for j in 0..atom.count {
                    let mut value = 0.;
                    for k in 0..atom.count {
                        value +=
                            atom.values[j * atom.count + k] * self.reduced[(atom.start + k, band)];
                    }
                    sum += self.reduced[(atom.start + j, band)] * value;
                }
            }
            *o = sum;
        }
        stage_add(Stage::NonlocalEnergy, start);
        COUNTS[Count::NonlocalEnergies as usize].fetch_add(1, Ordering::Relaxed);
        Ok(())
    }

    /// `reduced` = sum over ranks of P_local^T x (count x m), identical on all ranks.
    fn projections_reduced<K: Collective>(
        &mut self,
        comm: &K,
        x: &RealMatrix,
    ) -> Result<(), DistributedRealError> {
        let (count, m) = (self.count, x.cols());
        let n_loc = self.rows.len();
        if n_loc > 0 && self.dgemm {
            self.gemm.multiply_slices_wide_into(
                GemmSite::Other,
                self.projectors.as_slice(),
                n_loc,
                count,
                n_loc,
                RealOp::Transpose,
                &x.as_slice()[..n_loc * m],
                m,
                n_loc,
                self.projections.as_mut().subcols_mut(0, m),
            )?;
        } else if n_loc > 0 {
            self.gemm.multiply_into(
                self.projectors.as_ref(),
                RealOp::Transpose,
                x.as_ref(),
                self.projections.as_mut().subcols_mut(0, m),
            )?;
        } else {
            self.projections.as_mut_slice()[..count * m].fill(0.);
        }
        comm.tree_sum_f64(
            &self.projections.as_slice()[..count * m],
            &mut self.reduced.as_mut_slice()[..count * m],
        )?;
        Ok(())
    }

    /// Library-boundary helper (allocates): local rows of a global ORC block.
    pub fn extract_local(&self, global: &RealMatrix) -> Result<RealMatrix, DistributedRealError> {
        if global.rows() != self.n_global {
            return Err(DistributedRealError::Shape);
        }
        let mut local = RealMatrix::zeros(self.rows.len(), global.cols())?;
        for (g, l) in global
            .as_slice()
            .chunks(self.n_global)
            .zip(local.as_mut_slice().chunks_mut(self.rows.len().max(1)))
        {
            for (v, &r) in l.iter_mut().zip(&self.rows) {
                *v = g[r];
            }
        }
        Ok(local)
    }
    /// Library-boundary helper: write this rank's rows into a global ORC block.
    pub fn scatter_local(
        &self,
        local: &RealMatrix,
        global: &mut RealMatrix,
    ) -> Result<(), DistributedRealError> {
        if global.rows() != self.n_global
            || local.rows() != self.rows.len()
            || local.cols() != global.cols()
        {
            return Err(DistributedRealError::Shape);
        }
        for (g, l) in global
            .as_mut_slice()
            .chunks_mut(self.n_global)
            .zip(local.as_slice().chunks(self.rows.len().max(1)))
        {
            for (&v, &r) in l.iter().zip(&self.rows) {
                g[r] = v;
            }
        }
        Ok(())
    }

    /// Overwrite `y` with H x on this rank's rows (both `local_len` x m, m <= capacity).
    /// Every rank must call with the same m.
    pub fn apply_block<K: Collective>(
        &mut self,
        comm: &K,
        x: &RealMatrix,
        y: &mut RealMatrix,
    ) -> Result<(), DistributedRealError> {
        let n_loc = self.rows.len();
        let m = x.cols();
        if comm.size() != self.ranks
            || comm.rank() != self.rank
            || x.rows() != n_loc
            || y.rows() != n_loc
            || y.cols() != m
            || m > self.capacity
        {
            return Err(DistributedRealError::Shape);
        }
        if m == 0 {
            return Ok(());
        }
        let pairs = m.div_ceil(2);
        let local_len = self.fft.local_len();
        let xs = x.as_slice();
        let mi = self.max_items;
        let ranks = self.ranks;

        // (1) mirror rows to the ranks that own the minus positions.
        let t0 = Instant::now();
        self.exchange_mirror_rows(comm, xs, m)?;
        stage_add(Stage::Mirror, t0);

        // (2) pack every pair, distributed local-potential apply (one transpose per
        // direction for the whole block when batched), stage mirror h.
        let t0 = Instant::now();
        for t in 0..pairs {
            Self::pack_pair(
                &self.src,
                &self.xrecv,
                xs,
                n_loc,
                m,
                t,
                &mut self.hbuf[t * local_len..(t + 1) * local_len],
            );
        }
        stage_add(Stage::PackUnpack, t0);
        let t0 = Instant::now();
        if self.batched {
            self.fft.apply_local_batch(
                comm,
                &mut self.hbuf[..pairs * local_len],
                &self.potential,
                pairs,
            )?;
        } else {
            for t in 0..pairs {
                self.fft.apply_local_batch(
                    comm,
                    &mut self.hbuf[t * local_len..(t + 1) * local_len],
                    &self.potential,
                    1,
                )?;
            }
        }
        stage_add(Stage::Fft, t0);
        let t0 = Instant::now();
        for t in 0..pairs {
            let out = &self.hbuf[t * local_len..(t + 1) * local_len];
            for &(k, id) in &self.in_items {
                let off = (id * pairs + t) * 2;
                self.hsend[off] = out[k].re;
                self.hsend[off + 1] = out[k].im;
            }
        }
        stage_add(Stage::PackUnpack, t0);

        // (3) mirror h back to the representative owners.
        if mi > 0 {
            let t0 = Instant::now();
            let blk = mi * 2 * pairs;
            comm.alltoall_f64(
                &self.hsend[..ranks * blk],
                blk,
                &mut self.hrecv[..ranks * blk],
            )?;
            stage_add(Stage::Mirror, t0);
        }

        // Unpack, kinetic term.
        let t0 = Instant::now();
        let ys = y.as_mut_slice();
        for t in 0..pairs {
            let (ca, cb) = (2 * t, (2 * t + 1 < m).then_some(2 * t + 1));
            let h = &self.hbuf[t * local_len..(t + 1) * local_len];
            for rep in &self.reps {
                let hp = h[rep.plus_pos];
                let hm = match rep.minus {
                    MinusAt::None => {
                        ys[ca * n_loc + rep.slot] = hp.re;
                        if let Some(cb) = cb {
                            ys[cb * n_loc + rep.slot] = hp.im;
                        }
                        continue;
                    }
                    MinusAt::Local(k) => h[k],
                    MinusAt::Remote(id) => {
                        let off = (id * pairs + t) * 2;
                        C::new(self.hrecv[off], self.hrecv[off + 1])
                    }
                };
                let even = (hp + hm.conj()) * 0.5;
                let odd = (hp - hm.conj()) * 0.5;
                ys[ca * n_loc + rep.slot] = even.re * SQRT_2;
                ys[ca * n_loc + rep.slot + 1] = even.im * SQRT_2;
                if let Some(cb) = cb {
                    ys[cb * n_loc + rep.slot] = odd.im * SQRT_2;
                    ys[cb * n_loc + rep.slot + 1] = -odd.re * SQRT_2;
                }
            }
        }
        for (j, (x, y)) in xs.iter().zip(ys.iter_mut()).enumerate() {
            *y += self.kinetic[j % n_loc] * x;
        }
        stage_add(Stage::PackUnpack, t0);

        // Nonlocal: P^T X local + tree_sum, replicated coupling, P W local.
        let t0 = Instant::now();
        let count = self.count;
        if count > 0 {
            self.projections_reduced(comm, x)?;
            for atom in &self.coupling {
                for band in 0..m {
                    for j in 0..atom.count {
                        let mut value = 0.;
                        for k in 0..atom.count {
                            value += atom.values[j * atom.count + k]
                                * self.reduced[(atom.start + k, band)];
                        }
                        self.weights[(atom.start + j, band)] = value;
                    }
                }
            }
            if n_loc > 0 && self.dgemm {
                self.gemm.multiply_slices_wide_into(
                    GemmSite::Other,
                    self.projectors.as_slice(),
                    n_loc,
                    count,
                    n_loc,
                    RealOp::Identity,
                    &self.weights.as_slice()[..count * m],
                    m,
                    count,
                    self.image.as_mut().subcols_mut(0, m),
                )?;
                for (y, v) in y.as_mut_slice().iter_mut().zip(self.image.as_slice()) {
                    *y += v;
                }
            } else if n_loc > 0 {
                self.gemm.multiply_into(
                    self.projectors.as_ref(),
                    RealOp::Identity,
                    self.weights.as_ref().subcols(0, m),
                    self.image.as_mut().subcols_mut(0, m),
                )?;
                for (y, v) in y.as_mut_slice().iter_mut().zip(self.image.as_slice()) {
                    *y += v;
                }
            }
        }
        stage_add(Stage::Nonlocal, t0);
        COUNTS[if self.batched {
            Count::BatchedApplies
        } else {
            Count::UnbatchedApplies
        } as usize]
            .fetch_add(1, Ordering::Relaxed);
        COUNTS[Count::Columns as usize].fetch_add(m, Ordering::Relaxed);
        Ok(())
    }
}

/// Process-wide counters and stage timers of the M2 distributed operator (cheap
/// relaxed atomics; read once at the reporting boundary).
#[derive(Clone, Copy)]
enum Stage {
    Mirror,
    PackUnpack,
    Fft,
    Nonlocal,
    Density,
    NonlocalEnergy,
    Gather,
}
const STAGE_NAMES: [&str; 7] = [
    "apply_mirror_exchange",
    "apply_pack_unpack_kinetic",
    "apply_fft_transposes",
    "apply_nonlocal",
    "band_densities",
    "nonlocal_energy",
    "gather_planes",
];
static STAGE_NANOS: [AtomicU64; 7] = [const { AtomicU64::new(0) }; 7];
fn stage_add(stage: Stage, start: Instant) {
    STAGE_NANOS[stage as usize].fetch_add(start.elapsed().as_nanos() as u64, Ordering::Relaxed);
}
#[derive(Clone, Copy)]
enum Count {
    BatchedApplies,
    UnbatchedApplies,
    Columns,
    Densities,
    NonlocalEnergies,
}
const COUNT_NAMES: [&str; 5] = [
    "batched_applies",
    "unbatched_applies",
    "applied_columns",
    "band_density_builds",
    "nonlocal_energy_evaluations",
];
static COUNTS: [AtomicUsize; 5] = [const { AtomicUsize::new(0) }; 5];

/// (counter name, value) list and (stage name, wall seconds) list.
pub type M2Counters = (Vec<(&'static str, usize)>, Vec<(&'static str, f64)>);

/// (counter name, value) and (stage name, wall seconds) of every M2 distributed
/// operator in this process since start.
pub fn m2_counters() -> M2Counters {
    (
        COUNT_NAMES
            .iter()
            .zip(&COUNTS)
            .map(|(n, c)| (*n, c.load(Ordering::Relaxed)))
            .collect(),
        STAGE_NAMES
            .iter()
            .zip(&STAGE_NANOS)
            .map(|(n, c)| (*n, c.load(Ordering::Relaxed) as f64 * 1e-9))
            .collect(),
    )
}

impl dft_solvers::distributed_davidson::DistributedRealOperator for DistributedRealHamiltonian {
    fn global_len(&self) -> usize {
        self.n_global
    }
    fn local_rows(&self) -> &[usize] {
        &self.rows
    }
    fn worker_pool(&self) -> Option<Arc<WorkerPool>> {
        Some(self.pool.clone())
    }
    fn apply_block<K: Collective>(
        &mut self,
        comm: &K,
        x_local: &RealMatrix,
        y_local: &mut RealMatrix,
    ) -> Result<(), dft_solvers::SolverError> {
        DistributedRealHamiltonian::apply_block(self, comm, x_local, y_local)
            .map_err(|e| dft_solvers::SolverError::Application(e.to_string()))
    }
}
