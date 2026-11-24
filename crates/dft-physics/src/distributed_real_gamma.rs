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
use crate::real_gamma::{RealGammaError, pack_projector_rows};
use dft_basis::{HalfSphere, PlaneWaveBasis};
use dft_mpi::{Collective, MpiError};
use dft_numerics::linalg::LinalgError;
use dft_numerics::{FftError, GemmSite, RealGemmWorkspace, RealMatrix, RealOp, pool::WorkerPool};
use dft_solvers::distributed_fft::{DistFftError, DistributedPairFft, position_owners};
pub use dft_solvers::ref_davidson::{APPLY_H_STAGES, ApplyH, take_apply_h_solves};
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
    /// M2-OWNED (`BAANDI_M2_LOCAL_FIELDS`): Some(epoch) when `potential` holds
    /// only this rank's planes of the owned field of that epoch (every other
    /// entry NaN, never read by the apply); None when it is the full grid copied
    /// from the complex Hamiltonian.
    owned_epoch: Option<u64>,
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
    direct_pack: bool,
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
    single_projection: bool,
    fuse_coupling: bool,
    back_accum: bool,
    reference_profile: bool,
    apply_times: Option<[f64; 9]>,
    lane: lane_profile::Profile,
    source_identity: Arc<()>,
    pool: Arc<WorkerPool>,
}

/// Full-sphere representative indices required by one rank's ORC rows.
/// Sorted for the sparse complex table, independent of the ORC packing order.
pub fn owned_projector_indices(
    basis: &PlaneWaveBasis,
    grid: &dft_basis::FftGridSpec,
    rank: usize,
    ranks: usize,
) -> Result<Vec<usize>, DistributedRealError> {
    if ranks == 0 || rank >= ranks {
        return Err(DistributedRealError::Shape);
    }
    let half = HalfSphere::new(basis, grid)?;
    let owners = position_owners(grid.dims(), grid.g_to_fft(), ranks)?;
    let mut indices: Vec<_> = half
        .plus()
        .iter()
        .map(|&g| g as usize)
        .filter(|&g| owners[g] == rank)
        .collect();
    indices.sort_unstable();
    Ok(indices)
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
        Self::from_parts(
            source.fft_grid(),
            basis,
            source.nonlocal(),
            source.local_potential_values().to_vec(),
            source.operator_identity(),
            capacity,
            rank,
            ranks,
            pool,
        )
    }

    /// Construct without a serial Hamiltonian, serial FFTs or full complex
    /// projector table. Potential entries start as poison; the caller must set
    /// owned planes before applying H. The current FFT apply contract still
    /// requires a full-indexed potential allocation (8*grid_points bytes).
    #[allow(clippy::too_many_arguments)]
    pub fn new_owned(
        basis: &PlaneWaveBasis,
        grid: &dft_basis::FftGridSpec,
        nonlocal: &crate::hamiltonian::NonlocalOperator,
        capacity: usize,
        rank: usize,
        ranks: usize,
        pool: Arc<WorkerPool>,
    ) -> Result<Self, DistributedRealError> {
        if ranks == 0 || rank >= ranks || capacity == 0 || basis.k_fractional() != [0.; 3] {
            return Err(DistributedRealError::Shape);
        }
        Self::from_parts(
            grid,
            basis,
            nonlocal,
            vec![f64::NAN; grid.dims().iter().product()],
            Arc::new(()),
            capacity,
            rank,
            ranks,
            pool,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn from_parts(
        grid: &dft_basis::FftGridSpec,
        basis: &PlaneWaveBasis,
        nl: &crate::hamiltonian::NonlocalOperator,
        potential: Vec<f64>,
        source_identity: Arc<()>,
        capacity: usize,
        rank: usize,
        ranks: usize,
        pool: Arc<WorkerPool>,
    ) -> Result<Self, DistributedRealError> {
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

        // Kinetic diagonal in ORC order (as RealHamiltonian) of the local rows
        // only (ISSUE-20: the same copied values, no full-length vector).
        let kinetic = rows
            .iter()
            .map(|&r| {
                if r == 0 {
                    0.5 * basis.momenta()[0].norm_squared
                } else {
                    0.5 * basis.momenta()[plus[r.div_ceil(2)] as usize].norm_squared
                }
            })
            .collect();

        // Projector rows of the local rows only (ISSUE-20: packed directly,
        // validated over all representatives as `pack_projector`).
        let count = nl.projector_count();
        let mut projectors = RealMatrix::zeros(n_loc, count)?;
        let mut coupling = Vec::new();
        let mut start = 0;
        for atom in nl.atoms() {
            let nc = atom.channel_count();
            for j in 0..nc {
                let col =
                    &mut projectors.as_mut_slice()[(start + j) * n_loc..(start + j + 1) * n_loc];
                if nl.has_owned_rows() {
                    for (out, &r) in col.iter_mut().zip(&rows) {
                        let g = if r == 0 {
                            0
                        } else {
                            half.plus()[r.div_ceil(2)] as usize
                        };
                        let z = nl
                            .owned_value(start + j, g)
                            .ok_or(DistributedRealError::Shape)?;
                        *out = if r == 0 {
                            z.re
                        } else if r % 2 == 1 {
                            SQRT_2 * z.re
                        } else {
                            SQRT_2 * z.im
                        };
                    }
                } else {
                    pack_projector_rows(&half, atom.columns().column(j), &rows, col)?;
                }
            }
            coupling.push(Coupling {
                start,
                count: nc,
                values: atom.coupling_hartree().to_vec(),
            });
            start += nc;
        }
        let pairs = capacity.div_ceil(2);
        let local_len = fft.local_len();
        Ok(Self {
            rank,
            ranks,
            n_global: n,
            capacity,
            kinetic,
            potential,
            owned_epoch: None,
            src,
            reps,
            out_items,
            in_items,
            max_items,
            batched: true,
            direct_pack: dft_numerics::switches::on(dft_numerics::switches::Id::FFT_DIRECT_PACK),
            hbuf: vec![C::default(); pairs * local_len],
            gather_scratch: Vec::new(),
            xsend: vec![0.; ranks * max_items * 2 * capacity],
            xrecv: vec![0.; ranks * max_items * 2 * capacity],
            hsend: vec![0.; ranks * max_items * 2 * pairs],
            hrecv: vec![0.; ranks * max_items * 2 * pairs],
            projections: RealMatrix::zeros(count, capacity)?,
            reduced: RealMatrix::zeros(count, capacity)?,
            weights: RealMatrix::zeros(
                if dft_numerics::switches::on(dft_numerics::switches::Id::NONLOCAL_FUSE_COUPLING) {
                    0
                } else {
                    count
                },
                capacity,
            )?,
            image: RealMatrix::zeros(n_loc, capacity)?,
            gemm: RealGemmWorkspace::with_pool(n_loc, count.max(capacity), pool.clone())?,
            dgemm: dft_numerics::switches::on(dft_numerics::switches::Id::M2_DGEMM),
            single_projection: dft_numerics::switches::on(
                dft_numerics::switches::Id::NONLOCAL_SINGLE_GEMM,
            ),
            fuse_coupling: dft_numerics::switches::on(
                dft_numerics::switches::Id::NONLOCAL_FUSE_COUPLING,
            ),
            back_accum: dft_numerics::switches::on(dft_numerics::switches::Id::NONLOCAL_BACK_ACCUM),
            reference_profile: dft_numerics::switches::on(
                dft_numerics::switches::Id::DAVIDSON_STAGES,
            ),
            apply_times: None,
            lane: lane_profile::Profile::new(coupling.iter().map(|c| c.count)),
            source_identity,
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
        self.owned_epoch = None;
        true
    }

    /// The identity checks of `update_from_complex` without copying the
    /// potential (M2-OWNED: the owned planes are set by `set_owned_potential`).
    pub fn same_operator(&self, source: &Hamiltonian, basis: &PlaneWaveBasis) -> bool {
        Arc::ptr_eq(&self.source_identity, &source.operator_identity())
            && source.basis_momenta() == basis.momenta()
            && basis.k_fractional() == [0.; 3]
            && self.potential.len() == source.local_potential_values().len()
    }

    /// M2-OWNED: set the local potential from this rank's planes (`plane_points()`
    /// values, [plane][row][z], the layout of `band_densities`) of an owned field
    /// tagged `epoch`. The apply reads only these planes, so the result is
    /// bit-identical to a full-grid potential with the same plane values. The
    /// other entries are NaN (set once on entering the owned mode), so a full-grid
    /// reader of `potential()` sees non-finite values instead of a stale field.
    pub fn set_owned_potential(
        &mut self,
        planes: &[f64],
        epoch: u64,
    ) -> Result<(), DistributedRealError> {
        if planes.len() != self.fft.plane_points() {
            return Err(DistributedRealError::Shape);
        }
        if self.owned_epoch.is_none() {
            self.potential.fill(f64::NAN);
        }
        write_planes(&self.fft, planes, &mut self.potential);
        self.owned_epoch = Some(epoch);
        Ok(())
    }

    /// Epoch of the owned field the potential holds (None: full grid).
    pub fn owned_epoch(&self) -> Option<u64> {
        self.owned_epoch
    }

    /// This rank's planes of a full row-major grid field into `out`
    /// (`plane_points()` values, [plane][row][z]).
    pub fn read_planes(&self, full: &[f64], out: &mut Vec<f64>) {
        let (rows, nz) = (self.fft.rows(), self.fft.shape()[2]);
        out.clear();
        for i in 0..self.fft.plane_count() {
            for y in 0..rows {
                let s = self.fft.grid_start(i, y);
                out.extend_from_slice(&full[s..s + nz]);
            }
        }
    }

    /// Write this rank's planes (`plane_points()` values) into a full row-major
    /// grid field (other entries untouched).
    pub fn write_planes(&self, planes: &[f64], full: &mut [f64]) {
        write_planes(&self.fft, planes, full);
    }

    /// First FFT axis (the axis of the owned real-space planes) and the grid
    /// index of this rank's first plane.
    pub fn plane_layout(&self) -> (usize, usize) {
        (self.fft.first_axis(), self.fft.plane_x(0))
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
        let start = dft_numerics::stage_clock::start();
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
        for (phi, &s) in out.iter_mut().zip(src) {
            *phi = Self::packed_value(s, xrecv, xs, n_loc, m, t);
        }
    }

    /// Same ORC arithmetic for both the buffered and direct-stick paths.
    #[inline]
    fn packed_value(s: Src, xrecv: &[f64], xs: &[f64], n_loc: usize, m: usize, t: usize) -> C {
        let (ca, cb) = (2 * t, (2 * t + 1 < m).then_some(2 * t + 1));
        let get =
            |col: usize, slot: usize| C::new(xs[col * n_loc + slot], xs[col * n_loc + slot + 1]);
        let recv = |col: usize, id: usize| {
            let off = (id * m + col) * 2;
            C::new(xrecv[off], xrecv[off + 1])
        };
        match s {
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
        self.band_densities_keeping(comm, x, weights, volume, out, &mut [])
    }

    /// `band_densities` that also copies the real-space slab of each of the
    /// leading `keep.len() / plane_points()` column pairs into `keep` (pair t at
    /// `keep[t * plane_points()..]`, [plane][row][z], re = column 2t, im = column
    /// 2t+1, the unscaled inverse transform) for `apply_local_change`
    /// (MPI-KEPT-START). The density bits do not depend on `keep`.
    pub fn band_densities_keeping<K: Collective>(
        &mut self,
        comm: &K,
        x: &RealMatrix,
        weights: &[f64],
        volume: f64,
        out: &mut [f64],
        keep: &mut [C],
    ) -> Result<(), DistributedRealError> {
        let start = dft_numerics::stage_clock::start();
        if weights.len() != x.cols() || !(volume.is_finite() && volume > 0.) {
            return Err(DistributedRealError::Shape);
        }
        let points = self.fft.plane_points();
        if out.len() != x.cols() * points {
            return Err(DistributedRealError::Shape);
        }
        self.band_planes(comm, x, Some(keep), |c, p, v| {
            out[c * points + p] = weights[c] * v * v / volume
        })?;
        stage_add(Stage::Density, start);
        dft_numerics::work_timers::dist_transforms(
            dft_numerics::work_timers::DistKind::Density,
            start,
            x.cols().div_ceil(2),
            x.cols(),
        );
        bump(Count::Densities, 1);
        Ok(())
    }

    /// The plane-local density sum_b weights[b] psi_b(r)^2 / volume of the first
    /// m = x.cols() ORC columns, accumulated in ascending band order into `out`
    /// (`plane_points()` values, overwritten): bit-identical to summing the
    /// per-band output of `band_densities` in ascending band order, without the
    /// m x plane_points buffer and its second pass (`BAANDI_M2_DENSITY=fused`).
    /// Same collectives as `band_densities`.
    pub fn band_density_sum<K: Collective>(
        &mut self,
        comm: &K,
        x: &RealMatrix,
        weights: &[f64],
        volume: f64,
        out: &mut [f64],
    ) -> Result<(), DistributedRealError> {
        self.band_density_sum_keeping(comm, x, weights, volume, out, &mut [])
    }

    /// `band_density_sum` that also keeps the leading pair slabs in `keep`, as
    /// `band_densities_keeping` (MPI-KEPT-START).
    pub fn band_density_sum_keeping<K: Collective>(
        &mut self,
        comm: &K,
        x: &RealMatrix,
        weights: &[f64],
        volume: f64,
        out: &mut [f64],
        keep: &mut [C],
    ) -> Result<(), DistributedRealError> {
        let start = dft_numerics::stage_clock::start();
        if weights.len() != x.cols()
            || !(volume.is_finite() && volume > 0.)
            || out.len() != self.fft.plane_points()
        {
            return Err(DistributedRealError::Shape);
        }
        out.fill(0.);
        self.band_planes(comm, x, Some(keep), |c, p, v| {
            out[p] += weights[c] * v * v / volume
        })?;
        stage_add(Stage::Density, start);
        dft_numerics::work_timers::dist_transforms(
            dft_numerics::work_timers::DistKind::Density,
            start,
            x.cols().div_ceil(2),
            x.cols(),
        );
        bump(Count::Densities, 1);
        Ok(())
    }

    /// Real-space values psi_b(r) of the first m = x.cols() ORC columns on this
    /// rank's FFT planes, `out[b * plane_points() + p]`, in the normalization of
    /// `band_densities` (density = weight * psi^2 / volume) and of the complex
    /// `to_real` transform of the Hermitian-completed column. Same collectives as
    /// `band_densities` (one mirror alltoall when split, one batched transpose).
    /// Used by the distributed certificate screen (`BAANDI_M2_SCREEN=dist`).
    pub fn band_values<K: Collective>(
        &mut self,
        comm: &K,
        x: &RealMatrix,
        out: &mut [f64],
    ) -> Result<(), DistributedRealError> {
        let start = dft_numerics::stage_clock::start();
        let points = self.fft.plane_points();
        if out.len() != x.cols() * points {
            return Err(DistributedRealError::Shape);
        }
        self.band_planes(comm, x, None, |c, p, v| out[c * points + p] = v)?;
        stage_add(Stage::Values, start);
        bump(Count::Values, 1);
        Ok(())
    }

    /// Shared body of `band_densities`, `band_density_sum` and `band_values`:
    /// transform the ORC columns of x to this rank's planes and call
    /// `sink(column, plane point, psi(r))`, pair by pair in ascending column order
    /// (column 2t before 2t+1) for every point, so a sink that accumulates sums
    /// each point in ascending band order.
    fn band_planes<K: Collective>(
        &mut self,
        comm: &K,
        x: &RealMatrix,
        mut keep: Option<&mut [C]>,
        mut sink: impl FnMut(usize, usize, f64),
    ) -> Result<(), DistributedRealError> {
        let n_loc = self.rows.len();
        let m = x.cols();
        if comm.size() != self.ranks
            || comm.rank() != self.rank
            || x.rows() != n_loc
            || m > self.capacity
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
        bump(Count::InversePairs, pairs);
        let (rows, nz) = (self.fft.rows(), self.fft.shape()[2]);
        let line = rows * nz;
        let pp = self.fft.plane_points();
        for t in 0..pairs {
            let (ca, cb) = (2 * t, (2 * t + 1 < m).then_some(2 * t + 1));
            for i in 0..self.fft.plane_count() {
                self.fft.to_real_plane_in(t, pairs, i)?;
                let slab = self.fft.plane();
                let base = i * line;
                // MPI-KEPT-START: the pair's slab for the next solve's potential change.
                if let Some(k) = keep.as_deref_mut()
                    && (t + 1) * pp <= k.len()
                {
                    k[t * pp + base..t * pp + base + line].copy_from_slice(&slab[..line]);
                }
                for (j, v) in slab.iter().take(line).enumerate() {
                    sink(ca, base + j, v.re);
                }
                if let Some(cb) = cb {
                    for (j, v) in slab.iter().take(line).enumerate() {
                        sink(cb, base + j, v.im);
                    }
                }
            }
        }
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
        let start = dft_numerics::stage_clock::start();
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
        let lane = self.lane.start(
            lane_profile::Stage::Audit,
            self.rows.len(),
            self.count,
            x.cols(),
        );
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
        drop(lane);
        stage_add(Stage::NonlocalEnergy, start);
        bump(Count::NonlocalEnergies, 1);
        Ok(())
    }

    fn apply_clock(&self) -> Option<Instant> {
        self.apply_times.is_some().then(Instant::now)
    }
    fn apply_elapsed(&mut self, slot: usize, start: Option<Instant>) {
        if let (Some(times), Some(start)) = (&mut self.apply_times, start) {
            times[slot] += start.elapsed().as_secs_f64();
        }
    }

    /// `reduced` = sum over ranks of P_local^T x (count x m), identical on all ranks.
    fn projections_reduced<K: Collective>(
        &mut self,
        comm: &K,
        x: &RealMatrix,
    ) -> Result<(), DistributedRealError> {
        let (count, m) = (self.count, x.cols());
        let n_loc = self.rows.len();
        let span = self
            .lane
            .start(lane_profile::Stage::Projection, n_loc, count, m);
        let apply_clock = self.apply_clock();
        let arrival_start = dft_numerics::stage_clock::arrival::clock();
        if n_loc > 0 && self.dgemm {
            // Both entries fit p*b <= max(p, capacity)^2. The single entry
            // removes repeated reduction-panel BLAS calls/packing on MKL;
            // pure builds retain the same deterministic panel arithmetic.
            let multiply = if self.single_projection {
                RealGemmWorkspace::multiply_slices_single_into
            } else {
                RealGemmWorkspace::multiply_slices_wide_into
            };
            multiply(
                &mut self.gemm,
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
        self.apply_elapsed(5, apply_clock);
        dft_numerics::stage_clock::arrival::stage(
            301,
            dft_numerics::stage_clock::elapsed_ns(arrival_start),
            n_loc * m,
            8 * (n_loc * count + n_loc * m + count * m),
        );
        dft_numerics::stage_clock::arrival::shape(n_loc, count, m);
        dft_numerics::stage_clock::arrival::site(301);
        drop(span);
        let _span = self
            .lane
            .start(lane_profile::Stage::Reduce, n_loc, count, m);
        let apply_clock = self.apply_clock();
        comm.tree_sum_f64(
            &self.projections.as_slice()[..count * m],
            &mut self.reduced.as_mut_slice()[..count * m],
        )?;
        self.apply_elapsed(6, apply_clock);
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
        use dft_solvers::ref_davidson::{qe_clock_end, qe_clock_start};
        let qe_h = self.reference_profile.then(qe_clock_start).flatten();
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
        // h_apply parts (BAANDI_DAVIDSON_STAGES=on; off: no clock, no-op).
        let mut parts = dft_numerics::work_timers::DistApply::begin();

        let qe_local = self.reference_profile.then(qe_clock_start).flatten();
        // (1) mirror rows to the ranks that own the minus positions.
        let t0 = dft_numerics::stage_clock::start();
        let apply_clock = self.apply_clock();
        self.exchange_mirror_rows(comm, xs, m)?;
        self.apply_elapsed(1, apply_clock);
        stage_add(Stage::Mirror, t0);

        // (2) Optionally produce each ORC packed coefficient directly at the
        // stick scatter. The returned hbuf and mirror/unpack path are unchanged.
        let t0 = dft_numerics::stage_clock::start();
        let apply_clock = self.apply_clock();
        if !self.direct_pack {
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
        }
        self.apply_elapsed(2, apply_clock);
        stage_add(Stage::PackUnpack, t0);
        let t0 = dft_numerics::stage_clock::start();
        if self.direct_pack {
            let src = &self.src;
            let xrecv = &self.xrecv;
            if self.batched {
                self.fft.apply_local_batch_from(
                    comm,
                    |t, g| Self::packed_value(src[g], xrecv, xs, n_loc, m, t),
                    &mut self.hbuf[..pairs * local_len],
                    &self.potential,
                    pairs,
                )?;
            } else {
                for t in 0..pairs {
                    self.fft.apply_local_batch_from(
                        comm,
                        |_, g| Self::packed_value(src[g], xrecv, xs, n_loc, m, t),
                        &mut self.hbuf[t * local_len..(t + 1) * local_len],
                        &self.potential,
                        1,
                    )?;
                }
            }
        } else if self.batched {
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
        parts.pair_fft(t0, pairs, m);
        bump(Count::InversePairs, pairs);
        bump(Count::ForwardPairs, pairs);
        // (3) mirror h back, unpack; kinetic term.
        self.mirror_back_and_unpack(comm, m, y, Some(xs))?;
        qe_clock_end(1, qe_local);
        parts.lap(dft_numerics::work_timers::DistPart::Local);

        // Nonlocal: P^T X local + tree_sum, replicated coupling, P W local.
        let t0 = dft_numerics::stage_clock::start();
        let count = self.count;
        if count > 0 {
            let qe_projection = self.reference_profile.then(qe_clock_start).flatten();
            self.projections_reduced(comm, x)?;
            qe_clock_end(3, qe_projection);
            let qe_back = self.reference_profile.then(qe_clock_start).flatten();
            parts.lap(dft_numerics::work_timers::DistPart::Projection);
            let span = self
                .lane
                .start(lane_profile::Stage::Coupling, n_loc, count, m);
            let apply_clock = self.apply_clock();
            // The local projection is dead after tree_sum. Fuse the coupling
            // result into that storage, retaining the global Q and the exact
            // atom/band/j/k arithmetic. No copy of W or in-place DQ hazard.
            let weights = if self.fuse_coupling {
                &mut self.projections
            } else {
                &mut self.weights
            };
            for atom in &self.coupling {
                for band in 0..m {
                    for j in 0..atom.count {
                        let mut value = 0.;
                        for k in 0..atom.count {
                            value += atom.values[j * atom.count + k]
                                * self.reduced[(atom.start + k, band)];
                        }
                        weights[(atom.start + j, band)] = value;
                    }
                }
            }
            self.apply_elapsed(7, apply_clock);
            drop(span);
            let apply_clock = self.apply_clock();
            let span = self.lane.start(lane_profile::Stage::Back, n_loc, count, m);
            let weights = if self.fuse_coupling {
                &self.projections
            } else {
                &self.weights
            };
            let accumulated = n_loc > 0
                && self.dgemm
                && self.back_accum
                && self.gemm.try_multiply_slices_single_add(
                    GemmSite::Other,
                    self.projectors.as_slice(),
                    n_loc,
                    count,
                    n_loc,
                    &weights.as_slice()[..count * m],
                    m,
                    count,
                    y,
                )?;
            if accumulated {
                // HX already contains local+kinetic. Beta=1 consumed that value;
                // neither the image buffer nor a separate Add pass is needed.
                drop(span);
            } else if n_loc > 0 && self.dgemm {
                self.gemm.multiply_slices_wide_into(
                    GemmSite::Other,
                    self.projectors.as_slice(),
                    n_loc,
                    count,
                    n_loc,
                    RealOp::Identity,
                    &weights.as_slice()[..count * m],
                    m,
                    count,
                    self.image.as_mut().subcols_mut(0, m),
                )?;
                drop(span);
                let _span = self.lane.start(lane_profile::Stage::Add, n_loc, count, m);
                for (y, v) in y.as_mut_slice().iter_mut().zip(self.image.as_slice()) {
                    *y += v;
                }
            } else if n_loc > 0 {
                self.gemm.multiply_into(
                    self.projectors.as_ref(),
                    RealOp::Identity,
                    weights.as_ref().subcols(0, m),
                    self.image.as_mut().subcols_mut(0, m),
                )?;
                drop(span);
                let _span = self.lane.start(lane_profile::Stage::Add, n_loc, count, m);
                for (y, v) in y.as_mut_slice().iter_mut().zip(self.image.as_slice()) {
                    *y += v;
                }
            }
            self.apply_elapsed(8, apply_clock);
            qe_clock_end(4, qe_back);
        }
        stage_add(Stage::Nonlocal, t0);
        parts.lap(dft_numerics::work_timers::DistPart::Back);
        parts.finish(m);
        COUNTS[if self.batched {
            Count::BatchedApplies
        } else {
            Count::UnbatchedApplies
        } as usize]
            .fetch_add(1, Ordering::Relaxed);
        bump(Count::Columns, m);
        qe_clock_end(0, qe_h);
        Ok(())
    }

    /// Stage the mirror h values of the `ceil(m/2)` pairs in `hbuf`, send them
    /// back to the representative owners (one alltoall when split) and unpack
    /// into `y` (overwritten); with `kinetic_x` the kinetic term of that block is
    /// added (the tail of `apply_block`, same arithmetic and order).
    fn mirror_back_and_unpack<K: Collective>(
        &mut self,
        comm: &K,
        m: usize,
        y: &mut RealMatrix,
        kinetic_x: Option<&[f64]>,
    ) -> Result<(), DistributedRealError> {
        let n_loc = self.rows.len();
        let pairs = m.div_ceil(2);
        let local_len = self.fft.local_len();
        let (mi, ranks) = (self.max_items, self.ranks);
        let t0 = dft_numerics::stage_clock::start();
        let apply_clock = self.apply_clock();
        for t in 0..pairs {
            let out = &self.hbuf[t * local_len..(t + 1) * local_len];
            for &(k, id) in &self.in_items {
                let off = (id * pairs + t) * 2;
                self.hsend[off] = out[k].re;
                self.hsend[off + 1] = out[k].im;
            }
        }
        self.apply_elapsed(3, apply_clock);
        stage_add(Stage::PackUnpack, t0);
        if mi > 0 {
            let t0 = dft_numerics::stage_clock::start();
            let blk = mi * 2 * pairs;
            let apply_clock = self.apply_clock();
            comm.alltoall_f64(
                &self.hsend[..ranks * blk],
                blk,
                &mut self.hrecv[..ranks * blk],
            )?;
            self.apply_elapsed(1, apply_clock);
            stage_add(Stage::Mirror, t0);
        }
        let t0 = dft_numerics::stage_clock::start();
        let apply_clock = self.apply_clock();
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
        self.apply_elapsed(3, apply_clock);
        let apply_clock = self.apply_clock();
        let qe_kinetic = self
            .reference_profile
            .then(dft_solvers::ref_davidson::qe_clock_start)
            .flatten();
        if let Some(xs) = kinetic_x {
            for (j, (x, y)) in xs.iter().zip(y.as_mut_slice().iter_mut()).enumerate() {
                *y += self.kinetic[j % n_loc] * x;
            }
        }
        self.apply_elapsed(4, apply_clock);
        dft_solvers::ref_davidson::qe_clock_end(2, qe_kinetic);
        stage_add(Stage::PackUnpack, t0);
        Ok(())
    }

    /// Local potential of the operator (Ha, full row-major grid).
    pub fn potential(&self) -> &[f64] {
        &self.potential
    }

    /// MPI-KEPT-START: `y` = (local potential `delta`) applied to the columns of
    /// x on this rank's rows (both `local_len` x m, m <= capacity): the local part
    /// of `apply_block` with `delta` in place of the operator's potential, no
    /// kinetic and no nonlocal term. Pairs t < `kept_pairs` take their
    /// real-space slabs from `planes` (`band_densities_keeping` of the same x;
    /// only the r -> G transform runs); the other pairs are packed and
    /// transformed both ways. With planes from the same x the result has the bits
    /// of the plain transform of every pair. Collectives (every rank calls with
    /// the same m and `kept_pairs`): the mirror-row alltoall and one transpose to
    /// r only when some pair is not kept, one transpose back, the mirror-h
    /// alltoall (when split).
    pub fn apply_local_change<K: Collective>(
        &mut self,
        comm: &K,
        x: &RealMatrix,
        planes: &[C],
        kept_pairs: usize,
        delta: &[f64],
        y: &mut RealMatrix,
    ) -> Result<(), DistributedRealError> {
        let n_loc = self.rows.len();
        let m = x.cols();
        let pairs = m.div_ceil(2);
        if comm.size() != self.ranks
            || comm.rank() != self.rank
            || x.rows() != n_loc
            || y.rows() != n_loc
            || y.cols() != m
            || m > self.capacity
            || kept_pairs > m / 2
            || planes.len() < kept_pairs * self.fft.plane_points()
            || delta.len() != self.potential.len()
        {
            return Err(DistributedRealError::Shape);
        }
        if m == 0 {
            return Ok(());
        }
        let local_len = self.fft.local_len();
        let xs = x.as_slice();
        let t0 = dft_numerics::stage_clock::start();
        if kept_pairs < pairs {
            self.exchange_mirror_rows(comm, xs, m)?;
        }
        stage_add(Stage::Mirror, t0);
        let t0 = dft_numerics::stage_clock::start();
        for t in kept_pairs..pairs {
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
        let t0 = dft_numerics::stage_clock::start();
        self.fft.apply_local_batch_kept(
            comm,
            &mut self.hbuf[..pairs * local_len],
            delta,
            pairs,
            kept_pairs,
            planes,
        )?;
        stage_add(Stage::Fft, t0);
        dft_numerics::work_timers::dist_transforms(
            dft_numerics::work_timers::DistKind::Kept,
            t0,
            pairs,
            m,
        );
        bump(Count::InversePairs, pairs - kept_pairs);
        bump(Count::ForwardPairs, pairs);
        bump(Count::KeptPlanePairs, kept_pairs);
        bump(Count::ChangeColumns, m);
        self.mirror_back_and_unpack(comm, m, y, None)
    }
}

/// Copy plane-local values ([plane][row][z]) into a full row-major grid.
fn write_planes(fft: &DistributedPairFft, planes: &[f64], full: &mut [f64]) {
    let (rows, nz) = (fft.rows(), fft.shape()[2]);
    for i in 0..fft.plane_count() {
        for y in 0..rows {
            let s = fft.grid_start(i, y);
            let at = (i * rows + y) * nz;
            full[s..s + nz].copy_from_slice(&planes[at..at + nz]);
        }
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
    Values,
}
const STAGE_NAMES: [&str; 8] = [
    "apply_mirror_exchange",
    "apply_pack_unpack_kinetic",
    "apply_fft_transposes",
    "apply_nonlocal",
    "band_densities",
    "nonlocal_energy",
    "gather_planes",
    "band_values",
];
const STAGE_CALL_NAMES: [&str; 8] = [
    "calls.apply_mirror_exchange",
    "calls.apply_pack_unpack_kinetic",
    "calls.apply_fft_transposes",
    "calls.apply_nonlocal",
    "calls.band_densities",
    "calls.nonlocal_energy",
    "calls.gather_planes",
    "calls.band_values",
];
static STAGE_NANOS: [AtomicU64; 8] = [const { AtomicU64::new(0) }; 8];
static STAGE_CALLS: [AtomicUsize; 8] = [const { AtomicUsize::new(0) }; 8];
/// F21: `start` from `dft_numerics::stage_clock::start` (`None` with
/// `BAANDI_STAGE_TIMERS=off`: the entry is counted, no clock is read).
fn stage_add(stage: Stage, start: Option<Instant>) {
    if dft_numerics::stage_clock::arrival::enabled() {
        dft_numerics::stage_clock::arrival::stage(
            201 + stage as usize,
            dft_numerics::stage_clock::elapsed_ns(start),
            0,
            0,
        );
    }
    STAGE_CALLS[stage as usize].fetch_add(1, Ordering::Relaxed);
    if start.is_some() {
        STAGE_NANOS[stage as usize].fetch_add(
            dft_numerics::stage_clock::elapsed_ns(start),
            Ordering::Relaxed,
        );
    }
}
#[derive(Clone, Copy)]
enum Count {
    BatchedApplies,
    UnbatchedApplies,
    Columns,
    Densities,
    NonlocalEnergies,
    Values,
    InversePairs,
    ForwardPairs,
    KeptPlanePairs,
    ChangeColumns,
}
const COUNT_NAMES: [&str; 10] = [
    "batched_applies",
    "unbatched_applies",
    "applied_columns",
    "band_density_builds",
    "nonlocal_energy_evaluations",
    "band_value_builds",
    // Pair transforms G -> r and r -> G (apply, density/values, potential change).
    "inverse_pair_transforms",
    "forward_pair_transforms",
    // MPI-KEPT-START: pairs of a potential change taken from kept slabs, and
    // columns given (V_new - V_old) X instead of a full apply.
    "kept_plane_pairs",
    "potential_change_columns",
];
static COUNTS: [AtomicUsize; 10] = [const { AtomicUsize::new(0) }; 10];
/// Add `v` to a counter and mirror the work counts (applied columns, potential
/// change columns, pair transforms, kept pairs) per thread for the engine's
/// phase profile (`dft_numerics::stage_clock::local`, off: one relaxed load).
fn bump(c: Count, v: usize) {
    COUNTS[c as usize].fetch_add(v, Ordering::Relaxed);
    use dft_numerics::stage_clock::local::{OP_COUNT, add_n};
    let slot = match c {
        Count::Columns => 0,
        Count::ChangeColumns => 1,
        Count::ForwardPairs => 2,
        Count::InversePairs => 3,
        Count::KeptPlanePairs => 4,
        _ => return,
    };
    add_n(OP_COUNT + slot, 0, v as u64);
}

/// (counter name, value) list and (stage name, wall seconds) list.
pub type M2Counters = (Vec<(&'static str, usize)>, Vec<(&'static str, f64)>);

/// (counter name, value) and (stage name, wall seconds) of every M2 distributed
/// operator in this process since start. The counters end with one
/// `calls.<stage>` entry per stage (entries, both time modes); the seconds are 0
/// with `BAANDI_STAGE_TIMERS=off`.
pub fn m2_counters() -> M2Counters {
    (
        COUNT_NAMES
            .iter()
            .zip(&COUNTS)
            .map(|(n, c)| (*n, c.load(Ordering::Relaxed)))
            .chain(
                STAGE_CALL_NAMES
                    .iter()
                    .zip(&STAGE_CALLS)
                    .map(|(n, c)| (*n, c.load(Ordering::Relaxed))),
            )
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
        use dft_numerics::real::fft_profile;
        use dft_solvers::ref_davidson::{apply_h_active, record_apply_h};
        let profile = self.reference_profile && apply_h_active();
        let before = profile.then(fft_profile::snapshot);
        self.apply_times = profile.then_some([0.; 9]);
        let start = self.apply_clock();
        let result = DistributedRealHamiltonian::apply_block(self, comm, x_local, y_local);
        self.apply_elapsed(0, start);
        if let Some(times) = self.apply_times.take()
            && result.is_ok()
        {
            record_apply_h(
                times,
                fft_profile::snapshot().since(before.unwrap()),
                x_local.cols(),
                self.dgemm,
                self.count > 0 && !self.rows.is_empty(),
            );
        }
        result.map_err(|e| dft_solvers::SolverError::Application(e.to_string()))
    }
}

/// NONLOCAL Phase 0: allocation-free, thread-local diagnostic counters. They are
/// read by the engine at phase boundaries; no numerical decision reads them.
pub mod lane_profile {
    use std::cell::Cell;
    use std::time::Instant;

    pub const STAGES: usize = 6;
    pub const SHAPES: usize = 8;
    #[derive(Clone, Copy)]
    pub enum Stage {
        Projection,
        Reduce,
        Coupling,
        Back,
        Add,
        Audit,
    }

    #[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
    pub struct Snapshot {
        pub nanos: [u64; STAGES],
        pub calls: [u64; STAGES],
        pub flops: [u64; STAGES],
        pub bytes: [u64; STAGES],
        /// Sum of rows, projector columns, bands; calls at widths 1, 2, 3..8,
        /// 9..32, 33+. Only projection calls contribute (not back products).
        pub shapes: [u64; SHAPES],
    }
    thread_local! {
        static DATA: Cell<Snapshot> = const { Cell::new(Snapshot {
            nanos: [0; STAGES], calls: [0; STAGES], flops: [0; STAGES],
            bytes: [0; STAGES], shapes: [0; SHAPES],
        }) };
    }
    pub fn snapshot() -> Snapshot {
        DATA.with(Cell::get)
    }
    impl Snapshot {
        pub fn since(self, earlier: Self) -> Self {
            let mut out = self;
            for i in 0..STAGES {
                out.nanos[i] -= earlier.nanos[i];
                out.calls[i] -= earlier.calls[i];
                out.flops[i] -= earlier.flops[i];
                out.bytes[i] -= earlier.bytes[i];
            }
            for i in 0..SHAPES {
                out.shapes[i] -= earlier.shapes[i];
            }
            out
        }
    }

    /// Cached at operator construction; OFF performs no clock, TLS access,
    /// counter arithmetic, allocation, lock or I/O. Only the enable branch remains.
    #[derive(Clone, Copy)]
    pub(crate) struct Profile {
        enabled: bool,
        coupling_entries: u64,
        audit_nonzero: Option<u64>,
    }
    impl Profile {
        pub(crate) fn new(counts: impl Iterator<Item = usize>) -> Self {
            let enabled = dft_numerics::switches::on(dft_numerics::switches::Id::PHASE_PROFILE);
            Self {
                enabled,
                audit_nonzero: None,
                coupling_entries: if enabled {
                    counts.map(|c| (c as u64) * (c as u64)).sum()
                } else {
                    0
                },
            }
        }
        /// Serial audit skips exactly-zero coupling entries; count them once,
        /// only when instrumentation is enabled, without changing the audit.
        pub(crate) fn with_sparse_audit(mut self, nonzero: impl FnOnce() -> usize) -> Self {
            if self.enabled {
                self.audit_nonzero = Some(nonzero() as u64);
            }
            self
        }
        pub(crate) fn start(
            self,
            stage: Stage,
            rows: usize,
            projectors: usize,
            bands: usize,
        ) -> Span {
            if !self.enabled
                || projectors == 0
                || bands == 0
                || (rows == 0 && matches!(stage, Stage::Projection | Stage::Back | Stage::Add))
            {
                return Span(None);
            }
            let (r, p, b) = (rows as u64, projectors as u64, bands as u64);
            let (flops, bytes) = match stage {
                Stage::Projection | Stage::Back => (2 * r * p * b, 8 * (r * p + r * b + p * b)),
                Stage::Reduce => (0, 8 * p * b), // payload, not network traffic
                Stage::Coupling => (
                    2 * self.coupling_entries * b,
                    8 * (self.coupling_entries + 2 * p * b),
                ),
                Stage::Add => (r * b, 24 * r * b),
                Stage::Audit => (
                    self.audit_nonzero
                        .map_or(2 * self.coupling_entries + 2 * p, |nz| 3 * nz)
                        * b,
                    8 * (self.coupling_entries + p * b + b),
                ),
            };
            Span(Some(Active {
                stage,
                rows: r,
                projectors: p,
                bands: b,
                flops,
                bytes,
                clock: dft_numerics::stage_clock::start(),
            }))
        }
    }
    struct Active {
        stage: Stage,
        rows: u64,
        projectors: u64,
        bands: u64,
        flops: u64,
        bytes: u64,
        clock: Option<Instant>,
    }
    pub(crate) struct Span(Option<Active>);
    impl Drop for Span {
        fn drop(&mut self) {
            let Some(a) = self.0.take() else {
                return;
            };
            let nanos = dft_numerics::stage_clock::elapsed_ns(a.clock);
            DATA.with(|data| {
                let mut s = data.get();
                let i = a.stage as usize;
                s.calls[i] += 1;
                s.nanos[i] += nanos;
                s.flops[i] += a.flops;
                s.bytes[i] += a.bytes;
                if matches!(a.stage, Stage::Projection) {
                    s.shapes[0] += a.rows;
                    s.shapes[1] += a.projectors;
                    s.shapes[2] += a.bands;
                    let bin = match a.bands {
                        1 => 3,
                        2 => 4,
                        3..=8 => 5,
                        9..=32 => 6,
                        _ => 7,
                    };
                    if a.bands > 0 {
                        s.shapes[bin] += 1;
                    }
                }
                data.set(s);
            });
        }
    }
    #[cfg(test)]
    mod tests {
        use super::*;
        #[test]
        fn nonlocal_profile_disabled_is_inert() {
            let before = snapshot();
            let p = Profile {
                enabled: false,
                coupling_entries: u64::MAX,
                audit_nonzero: None,
            };
            drop(p.start(Stage::Projection, usize::MAX, usize::MAX, usize::MAX));
            assert_eq!(snapshot(), before);
        }
        #[test]
        fn nonlocal_profile_empty_products_do_not_report_traffic() {
            let before = snapshot();
            let p = Profile {
                enabled: true,
                coupling_entries: 0,
                audit_nonzero: None,
            };
            drop(p.start(Stage::Projection, 7, 0, 2));
            drop(p.start(Stage::Projection, 7, 3, 0));
            drop(p.start(Stage::Back, 0, 3, 2));
            assert_eq!(snapshot(), before);
        }
        #[test]
        fn nonlocal_profile_count_only_and_delta() {
            let before = snapshot();
            drop(Span(Some(Active {
                stage: Stage::Projection,
                rows: 7,
                projectors: 3,
                bands: 2,
                flops: 84,
                bytes: 328,
                clock: None,
            })));
            let d = snapshot().since(before);
            assert_eq!(d.calls, [1, 0, 0, 0, 0, 0]);
            assert_eq!(d.nanos, [0; STAGES]);
            assert_eq!(d.flops[0], 84);
            assert_eq!(d.bytes[0], 328);
            assert_eq!(d.shapes, [7, 3, 2, 0, 1, 0, 0, 0]);
        }
    }
}

#[cfg(test)]
mod lane_identity_tests {
    use super::*;
    use crate::hamiltonian::NonlocalOperator;
    use dft_geometry::Cell;
    use dft_numerics::switches::RunConfig;
    use dft_units::Hartree;

    #[test]
    fn nonlocal_single_projection_matches_panel_route() {
        let mut worst = 0.0_f64;
        let mut comparisons = 0;
        for lattice in [
            // More than 256 rows on rank 0 at P=1: cross the reduction-panel boundary.
            [[20., 0., 0.], [0., 20., 0.], [0., 0., 20.]],
            [[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]],
            [[2., 0., 0.], [0., 2., 0.], [0., 0., 2.]],
            [[6., 0.4, 0.2], [1., 7., 0.3], [0.5, 0.7, 8.]],
            [[4., 0., 0.], [0., 5., 0.], [0., 0., 14.]],
        ] {
            for ranks in [1, 2, 3, 8] {
                let results = dft_mpi::simulate(ranks, |comm| {
                    let cell = Cell::from_user_rows_bohr(lattice).unwrap();
                    let basis = PlaneWaveBasis::new(&cell, Hartree(1.2)).unwrap();
                    if lattice[0][0] == 20. && ranks == 1 {
                        assert!(basis.vectors().len() > 256);
                    }
                    let grid = dft_basis::FftGridSpec::new(
                        &basis,
                        &basis.density_basis(&cell).unwrap(),
                        None,
                    )
                    .unwrap();
                    let atoms = [0, 1, 3, 5]
                        .into_iter()
                        .map(|nc| {
                            let columns = (0..nc)
                                .map(|j| {
                                    basis
                                        .momenta()
                                        .iter()
                                        .map(|g| {
                                            C::from_polar(
                                                0.1 / (1. + g.norm_squared),
                                                g.cartesian[0] * (j + 1) as f64 * 0.13,
                                            )
                                        })
                                        .collect()
                                })
                                .collect();
                            let d = (0..nc * nc)
                                .map(|i| if i / nc == i % nc { 0.7 } else { -0.03 })
                                .collect();
                            (columns, d)
                        })
                        .collect();
                    let source = Hamiltonian::new(
                        &cell,
                        &basis,
                        ndarray::Array3::zeros(grid.dims()),
                        NonlocalOperator::from_columns(basis.momenta().to_vec(), atoms, 5),
                    )
                    .unwrap();
                    let mut worst = 0.0_f64;
                    for width in [0, 1, 2, 5, 17] {
                        let mut baseline: Option<RealMatrix> = None;
                        let mut unfused: Option<RealMatrix> = None;
                        let mut unfused_energies = Vec::new();
                        for (single, fuse, back) in [
                            ("off", "off", "off"),
                            ("on", "off", "off"),
                            ("on", "on", "off"),
                            ("on", "off", "on"),
                            ("on", "on", "on"),
                        ] {
                            let cfg = RunConfig::from_pairs(&[
                                ("nonlocal_single_gemm", single),
                                ("nonlocal_fuse_coupling", fuse),
                                ("nonlocal_back_accum", back),
                                ("m2_dgemm", "on"),
                                ("mkl_dgemm", "on"),
                            ])
                            .unwrap();
                            let _scope = cfg.enter();
                            let mut op = DistributedRealHamiltonian::new(
                                &source,
                                &basis,
                                17,
                                comm.rank(),
                                ranks,
                                WorkerPool::new(if fuse == "on" { 3 } else { 1 }).unwrap(),
                            )
                            .unwrap();
                            let mut x = RealMatrix::zeros(op.local_len(), width).unwrap();
                            for j in 0..width {
                                for (i, &r) in op.local_rows().iter().enumerate() {
                                    x[(i, j)] = ((r * 17 + j * 13) as f64).sin() / 10.;
                                }
                            }
                            let mut y = RealMatrix::zeros(op.local_len(), width).unwrap();
                            op.image.as_mut_slice().fill(f64::NAN);
                            op.apply_block(comm, &x, &mut y).unwrap();
                            if back == "on" && dft_numerics::MKL_BACKEND {
                                assert!(op.image.as_slice().iter().all(|v| v.is_nan()));
                            }
                            let mut energies = vec![0.; width];
                            op.nonlocal_band_energies(comm, &x, &mut energies).unwrap();
                            if fuse == "on" {
                                assert!(op.weights.as_slice().is_empty());
                                for (&a, &b) in unfused
                                    .as_ref()
                                    .unwrap()
                                    .as_slice()
                                    .iter()
                                    .zip(y.as_slice())
                                {
                                    assert_eq!(a.to_bits(), b.to_bits());
                                }
                                assert_eq!(energies, unfused_energies);
                                // An energy audit overwrites the recycled local
                                // projection. A subsequent H must still rebuild W.
                                let mut again = RealMatrix::zeros(op.local_len(), width).unwrap();
                                op.apply_block(comm, &x, &mut again).unwrap();
                                for (&a, &b) in y.as_slice().iter().zip(again.as_slice()) {
                                    assert_eq!(a.to_bits(), b.to_bits());
                                }
                            } else if single == "on" {
                                unfused = Some(y.clone());
                                unfused_energies = energies;
                            }
                            if let Some(base) = &baseline {
                                for (&a, &b) in base.as_slice().iter().zip(y.as_slice()) {
                                    worst = worst.max((a - b).abs());
                                    assert!((a - b).abs() < 1e-13 * a.abs().max(1.));
                                    if !dft_numerics::MKL_BACKEND {
                                        assert_eq!(a.to_bits(), b.to_bits());
                                    }
                                }
                            } else {
                                baseline = Some(y);
                            }
                        }
                    }
                    worst
                });
                for delta in results {
                    worst = worst.max(delta);
                }
                comparisons += ranks * 5;
            }
        }
        eprintln!(
            "P3-ADDV: {comparisons} rank/block fixtures, projection and coupling arms, max H component delta {worst:e} Ha"
        );
    }

    #[test]
    fn a2a_balance_h_equivalence() {
        type RankResult = (Vec<usize>, usize, Vec<(RealMatrix, (usize, usize))>);
        let mut worst = 0.0f64;
        let mut changed = 0;
        for (lattice, shape) in [
            ([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]], [12, 12, 12]),
            (
                [[6., 0.4, 0.2], [1., 7., 0.3], [0.5, 0.7, 8.]],
                [11, 13, 17],
            ),
            ([[4., 0., 0.], [0., 5., 0.], [0., 0., 12.]], [7, 11, 23]),
        ] {
            for ranks in [1, 2, 3, 8] {
                for mirror in ["greedy", "mirror"] {
                    for channels in [0, 3] {
                        let run = |arm| {
                            dft_mpi::simulate(ranks, |comm| {
                                let config = RunConfig::from_pairs(&[
                                    ("fft_a2a_balance", arm),
                                    ("m2_sticks", mirror),
                                    ("fft_direct_pack", "on"),
                                    ("fft_sparse_clear", "on"),
                                ])
                                .unwrap();
                                let _scope = config.enter();
                                let cell = Cell::from_user_rows_bohr(lattice).unwrap();
                                let basis = PlaneWaveBasis::new(&cell, Hartree(1.2)).unwrap();
                                let atoms = if channels == 0 {
                                    vec![]
                                } else {
                                    [1, channels, 5]
                                        .into_iter()
                                        .map(|nc| {
                                            let columns = (0..nc)
                                                .map(|j| {
                                                    basis
                                                        .momenta()
                                                        .iter()
                                                        .map(|g| {
                                                            C::from_polar(
                                                                0.1 / (1. + g.norm_squared),
                                                                g.cartesian[0]
                                                                    * (j + 1) as f64
                                                                    * 0.13,
                                                            )
                                                        })
                                                        .collect()
                                                })
                                                .collect();
                                            let coupling = (0..nc * nc)
                                                .map(|i| if i / nc == i % nc { 0.7 } else { -0.03 })
                                                .collect();
                                            (columns, coupling)
                                        })
                                        .collect()
                                };
                                let nl = NonlocalOperator::from_columns(
                                    basis.momenta().to_vec(),
                                    atoms,
                                    5,
                                );
                                let potential =
                                    ndarray::Array3::from_shape_fn(shape, |(i, j, k)| {
                                        0.1 * ((i + 2 * j + 3 * k) as f64).sin()
                                    });
                                let source =
                                    Hamiltonian::new(&cell, &basis, potential, nl).unwrap();
                                let mut op = DistributedRealHamiltonian::new(
                                    &source,
                                    &basis,
                                    16,
                                    comm.rank(),
                                    ranks,
                                    WorkerPool::new(1).unwrap(),
                                )
                                .unwrap();
                                let mut results = Vec::new();
                                for batched in [true, false] {
                                    op.set_batched(batched);
                                    for m in [16, 1, 2, 5, 0] {
                                        let x = RealMatrix::new(
                                            op.local_len(),
                                            m,
                                            (0..m)
                                                .flat_map(|t| {
                                                    op.local_rows().iter().map(move |&g| {
                                                        ((g * 17 + t * 13) as f64).sin() * 0.5
                                                    })
                                                })
                                                .collect(),
                                        )
                                        .unwrap();
                                        let mut y = RealMatrix::zeros(op.local_len(), m).unwrap();
                                        op.apply_block(comm, &x, &mut y).unwrap();
                                        results.push((y, op.collectives_per_apply(m)));
                                    }
                                }
                                (op.local_rows().to_vec(), op.global_len(), results)
                            })
                        };
                        let off = run("off");
                        let on = run("on");
                        if off.iter().zip(&on).any(|(a, b)| a.0 != b.0) {
                            changed += 1;
                        }
                        for (j, m) in [16, 1, 2, 5, 0].into_iter().cycle().take(10).enumerate() {
                            let assemble = |parts: &[RankResult]| {
                                let mut out = vec![0.; parts[0].1 * m];
                                for (rows, n, results) in parts {
                                    for t in 0..m {
                                        for (k, &g) in rows.iter().enumerate() {
                                            out[t * n + g] = results[j].0[(k, t)];
                                        }
                                    }
                                }
                                out
                            };
                            let a = assemble(&off);
                            let b = assemble(&on);
                            for (a, b) in a.iter().zip(&b) {
                                worst = worst.max((a - b).abs());
                                assert!(
                                    (a - b).abs() < 1e-12,
                                    "H difference {shape:?} ranks={ranks} {mirror}"
                                );
                                if channels == 0 || ranks == 1 {
                                    assert_eq!(a.to_bits(), b.to_bits());
                                }
                            }
                            if mirror == "mirror" {
                                for (a, b) in off.iter().zip(&on) {
                                    assert_eq!(a.2[j].1, b.2[j].1);
                                }
                            }
                        }
                    }
                }
            }
        }
        assert!(changed > 0, "exercise changed ownership, not just fallback");
        eprintln!(
            "P3-A2A: H max absolute difference {worst:e} Ha; {changed} layouts changed; local+kinetic bit-identical"
        );
    }

    #[test]
    fn fft_pack_switches_preserve_distributed_h_bits() {
        let mut comparisons = 0;
        for lattice in [
            [[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]],
            [[6., 0.4, 0.2], [1., 7., 0.3], [0.5, 0.7, 8.]],
            [[4., 0., 0.], [0., 5., 0.], [0., 0., 12.]],
        ] {
            for ranks in [1, 2, 3, 8] {
                for mirror in ["greedy", "mirror"] {
                    dft_mpi::simulate(ranks, |comm| {
                        let cell = Cell::from_user_rows_bohr(lattice).unwrap();
                        let basis = PlaneWaveBasis::new(&cell, Hartree(1.2)).unwrap();
                        let grid = dft_basis::FftGridSpec::new(
                            &basis,
                            &basis.density_basis(&cell).unwrap(),
                            None,
                        )
                        .unwrap();
                        let potential = ndarray::Array3::from_shape_fn(grid.dims(), |(i, j, k)| {
                            0.1 * ((i + 2 * j + 3 * k) as f64).sin()
                        });
                        let nl =
                            NonlocalOperator::from_columns(basis.momenta().to_vec(), vec![], 0);
                        let source = Hamiltonian::new(&cell, &basis, potential, nl).unwrap();
                        let pool = WorkerPool::new(1).unwrap();
                        let mut reference = Vec::new();
                        for (direct, sparse, rows) in [
                            ("off", "off", "off"),
                            ("on", "off", "off"),
                            ("off", "on", "off"),
                            ("on", "on", "off"),
                            ("off", "off", "on"),
                            ("on", "on", "on"),
                        ] {
                            let config = RunConfig::from_pairs(&[
                                ("fft_direct_pack", direct),
                                ("fft_sparse_clear", sparse),
                                ("fft_fused_rows", rows),
                                ("m2_sticks", mirror),
                            ])
                            .unwrap();
                            let _guard = config.enter();
                            let mut op = DistributedRealHamiltonian::new(
                                &source,
                                &basis,
                                5,
                                comm.rank(),
                                ranks,
                                pool.clone(),
                            )
                            .unwrap();
                            let mut bits = Vec::new();
                            for batched in [true, false] {
                                op.set_batched(batched);
                                for m in [5, 1, 2, 3, 0] {
                                    let n = op.local_len();
                                    let x = RealMatrix::new(
                                        n,
                                        m,
                                        (0..n * m)
                                            .map(|i| {
                                                if i % 7 == 0 {
                                                    -0.0
                                                } else {
                                                    ((i * 17 % 113) as f64 - 56.) / 113.
                                                }
                                            })
                                            .collect(),
                                    )
                                    .unwrap();
                                    let mut y = RealMatrix::zeros(n, m).unwrap();
                                    op.apply_block(comm, &x, &mut y).unwrap();
                                    bits.extend(y.as_slice().iter().map(|v| v.to_bits()));
                                }
                            }
                            if direct == "off" && sparse == "off" && rows == "off" {
                                reference = bits;
                            } else {
                                assert_eq!(bits, reference, "ranks {ranks} mirror {mirror}");
                            }
                        }
                    });
                    comparisons += ranks * 5 * 2 * 5;
                }
            }
        }
        eprintln!(
            "P3-FFT: {comparisons} rank-local H on/off comparisons, max difference 0 Ha (bitwise)"
        );
    }

    /// Both MPI-using tests share one Universe: its Drop finalizes MPI, and MPI
    /// cannot be initialized twice in one process (--features mpi).
    #[test]
    fn mpi_universe_tests() {
        let universe = dft_mpi::Universe::init().unwrap();
        reference_h_apply_profile_is_bit_identical(&universe);
        nonlocal_profile_preserves_operator_and_energy_bits(&universe);
    }

    fn reference_h_apply_profile_is_bit_identical(universe: &dft_mpi::Universe) {
        use dft_solvers::ref_davidson::{Options, solve, take_apply_h_solves};
        let comm = universe.world();
        for bands in [2, 3] {
            let cell =
                Cell::from_user_rows_bohr([[6., 0.4, 0.2], [1., 7., 0.3], [0.5, 0.7, 8.]]).unwrap();
            let basis = PlaneWaveBasis::new(&cell, Hartree(1.2)).unwrap();
            let grid =
                dft_basis::FftGridSpec::new(&basis, &basis.density_basis(&cell).unwrap(), None)
                    .unwrap();
            let columns = (0..3)
                .map(|j| {
                    basis
                        .momenta()
                        .iter()
                        .map(|g| {
                            C::from_polar(
                                0.1 / (1. + g.norm_squared),
                                g.cartesian[0] * (j + 1) as f64 * 0.13,
                            )
                        })
                        .collect()
                })
                .collect();
            let nl = NonlocalOperator::from_columns(
                basis.momenta().to_vec(),
                vec![(
                    columns,
                    vec![0.7, -0.03, 0., -0.03, 0.6, 0.02, 0., 0.02, 0.8],
                )],
                3,
            );
            let potential = ndarray::Array3::from_shape_fn(grid.dims(), |(i, j, k)| {
                0.1 * ((i + 2 * j + 3 * k) as f64).sin()
            });
            let source = Hamiltonian::new(&cell, &basis, potential, nl).unwrap();
            let mut baseline = None;
            for mode in ["off", "on", "fft"] {
                let config = RunConfig::from_pairs(&[
                    ("davidson_stages", if mode == "off" { "off" } else { "on" }),
                    ("phase_profile", if mode == "fft" { "on" } else { "off" }),
                    ("fft_profile", if mode == "fft" { "on" } else { "off" }),
                    ("stage_timers", if mode == "fft" { "on" } else { "off" }),
                ])
                .unwrap();
                let _scope = config.enter();
                take_apply_h_solves();
                let mut op = DistributedRealHamiltonian::new(
                    &source,
                    &basis,
                    bands,
                    0,
                    1,
                    WorkerPool::new(1).unwrap(),
                )
                .unwrap();
                let rows = op.local_len();
                let x = RealMatrix::new(
                    rows,
                    bands,
                    (0..rows * bands)
                        .map(|i| ((i * 17 % 113) as f64 - 56.) / 113.)
                        .collect(),
                )
                .unwrap();
                let diagonal = op.kinetic.clone();
                let mut options = Options::new(bands, 1e-7, bands);
                options.profile = mode != "off";
                options.return_hx = true;
                let qe_scope = (mode != "off").then(dft_solvers::ref_davidson::QeClockScope::begin);
                let fft_before = dft_numerics::real::fft_profile::snapshot();
                let result = solve(&comm, &mut op, &x, &diagonal, options).unwrap();
                if let Some(scope) = &qe_scope {
                    let clocks = scope.snapshot();
                    for slot in [0, 1, 2, 3, 4, 6] {
                        assert_eq!(clocks.calls[slot], result.record.h_calls);
                        assert!(clocks.seconds[slot] > 0.);
                    }
                    assert!(clocks.seconds[0] >= clocks.seconds[1]);
                    assert!(clocks.seconds[1] >= clocks.seconds[2]);
                    assert_eq!(clocks.calls[5], 1);
                } else {
                    assert!(dft_solvers::ref_davidson::qe_clock_start().is_none());
                }
                let bits = result
                    .eigenvalues
                    .iter()
                    .chain(result.eigenvectors.as_slice())
                    .chain(result.hx.as_ref().unwrap().as_slice())
                    .map(|v| v.to_bits())
                    .collect::<Vec<_>>();
                if let Some(previous) = &baseline {
                    assert_eq!(&bits, previous);
                } else {
                    baseline = Some(bits);
                }
                let records = take_apply_h_solves();
                if mode == "off" {
                    assert_eq!(
                        dft_numerics::real::fft_profile::snapshot().since(fft_before),
                        dft_numerics::real::fft_profile::Snapshot::default()
                    );
                    assert!(records.is_empty());
                    assert!(result.record.apply_h.is_none());
                } else {
                    assert_eq!(records.len(), 1);
                    let a = &records[0];
                    assert_eq!(a.columns_per_call.len(), result.record.h_calls);
                    assert_eq!(
                        a.columns_per_call.iter().sum::<usize>(),
                        result.record.h_columns
                    );
                    assert_eq!(2 * a.pairs - a.odd_tails, result.record.h_columns);
                    assert_eq!(a.wide_calls, 2 * result.record.h_calls);
                    assert_eq!(a.direct_calls, 0);
                    if mode != "off" {
                        for i in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16] {
                            assert!(a.fft.0[i][1] > 0, "missing FFT stage {i}");
                            assert!(a.fft.0[i][0] > 0, "missing FFT clock {i}");
                        }
                        assert_eq!(a.fft.0[0][1] as usize, a.pairs);
                        assert_eq!(a.fft.0[16][1] as usize, a.pairs);
                    }
                    assert!(a.seconds[0] > 0. && a.seconds[5] > 0. && a.seconds[8] > 0.);
                    println!(
                        "P3-HAPPLY mode={mode} bands={bands} calls={} columns={} pairs={} odd_tails={} bit-identical",
                        result.record.h_calls, result.record.h_columns, a.pairs, a.odd_tails
                    );
                }
            }
        }
    }

    fn nonlocal_profile_preserves_operator_and_energy_bits(universe: &dft_mpi::Universe) {
        let mut comparisons = 0;
        for lattice in [
            [[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]],
            [[6., 0.4, 0.2], [1., 7., 0.3], [0.5, 0.7, 8.]],
            [[4., 0., 0.], [0., 5., 0.], [0., 0., 14.]],
        ] {
            let cell = Cell::from_user_rows_bohr(lattice).unwrap();
            let basis = PlaneWaveBasis::new(&cell, Hartree(1.2)).unwrap();
            let grid =
                dft_basis::FftGridSpec::new(&basis, &basis.density_basis(&cell).unwrap(), None)
                    .unwrap();
            for channels in [0, 1, 7] {
                let columns = (0..channels)
                    .map(|j| {
                        basis
                            .momenta()
                            .iter()
                            .map(|g| {
                                C::from_polar(
                                    0.1 / (1. + g.norm_squared),
                                    g.cartesian[0] * (j + 1) as f64 * 0.13
                                        + g.cartesian[1] * 0.2
                                        + g.cartesian[2] * 0.3,
                                )
                            })
                            .collect()
                    })
                    .collect();
                let coupling = (0..channels * channels)
                    .map(|i| {
                        if i / channels == i % channels {
                            0.7
                        } else {
                            -0.03
                        }
                    })
                    .collect();
                let nl = NonlocalOperator::from_columns(
                    basis.momenta().to_vec(),
                    vec![(columns, coupling)],
                    channels,
                );
                let source =
                    Hamiltonian::new(&cell, &basis, ndarray::Array3::zeros(grid.dims()), nl)
                        .unwrap();
                let mut op = DistributedRealHamiltonian::new(
                    &source,
                    &basis,
                    17,
                    0,
                    1,
                    WorkerPool::new(1).unwrap(),
                )
                .unwrap();
                let comm = universe.world();
                for m in [0, 1, 2, 3, 17] {
                    let n = op.rows.len();
                    let x = RealMatrix::new(
                        n,
                        m,
                        (0..n * m)
                            .map(|i| ((i * 17 % 113) as f64 - 56.) / 113.)
                            .collect(),
                    )
                    .unwrap();
                    let mut images = Vec::new();
                    let mut energies = Vec::new();
                    for mode in ["off", "on"] {
                        let cfg = RunConfig::from_pairs(&[("phase_profile", mode)]).unwrap();
                        let _guard = cfg.enter();
                        op.lane = lane_profile::Profile::new(op.coupling.iter().map(|c| c.count));
                        let mut image = RealMatrix::zeros(n, m).unwrap();
                        let mut energy = vec![0.; m];
                        op.apply_block(&comm, &x, &mut image).unwrap();
                        op.nonlocal_band_energies(&comm, &x, &mut energy).unwrap();
                        images.push(
                            image
                                .as_slice()
                                .iter()
                                .map(|v| v.to_bits())
                                .collect::<Vec<_>>(),
                        );
                        energies.push(energy.iter().map(|v| v.to_bits()).collect::<Vec<_>>());
                    }
                    assert_eq!(images[0], images[1]);
                    assert_eq!(energies[0], energies[1]);
                    comparisons += 1;
                }
            }
        }
        println!("NONLOCAL_PROFILE {comparisons} operator/energy pairs bit-identical");
    }
}
