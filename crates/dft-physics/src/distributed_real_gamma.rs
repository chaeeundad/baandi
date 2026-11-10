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
//! Collective calls per apply with m columns, P = ceil(m/2) pairs:
//! alltoall = 2P + 2 (the 2 mirror exchanges are skipped when no representative is
//! split, e.g. one rank), tree_sum = 1 (0 without projectors). See
//! `collectives_per_apply`. All counts are rank independent (no deadlock).
//!
//! Numerics. One rank reproduces `RealHamiltonian::apply_block` to rounding; it is
//! bit-identical when the serial and distributed FFT paths agree bitwise (tested,
//! reported by the gate test). Default GEMM path only (`BAANDI_MKL_NONLOCAL` is not
//! honoured). Inputs are not scanned for NaN/inf: an early error on one rank only
//! would desynchronise the collectives; non-finite values propagate to the output.
//! No allocation after construction (conversion helpers `extract_local` and
//! `scatter_local` allocate and are for the library boundary); no unsafe; no logging.
use crate::hamiltonian::Hamiltonian;
use crate::real_gamma::{RealGammaError, pack_projector};
use dft_basis::{HalfSphere, PlaneWaveBasis};
use dft_mpi::{Collective, MpiError};
use dft_numerics::linalg::LinalgError;
use dft_numerics::{FftError, RealGemmWorkspace, RealMatrix, RealOp, pool::WorkerPool};
use dft_solvers::distributed_fft::{DistFftError, DistributedPairFft, position_owners};
use num_complex::Complex64 as C;
use std::f64::consts::SQRT_2;
use std::sync::Arc;
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
    cvec: Vec<C>,
    hbuf: Vec<C>,
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
            cvec: vec![C::default(); local_len],
            hbuf: vec![C::default(); pairs * local_len],
            xsend: vec![0.; ranks * max_items * 2 * capacity],
            xrecv: vec![0.; ranks * max_items * 2 * capacity],
            hsend: vec![0.; ranks * max_items * 2 * pairs],
            hrecv: vec![0.; ranks * max_items * 2 * pairs],
            projections: RealMatrix::zeros(count, capacity)?,
            reduced: RealMatrix::zeros(count, capacity)?,
            weights: RealMatrix::zeros(count, capacity)?,
            image: RealMatrix::zeros(n_loc, capacity)?,
            gemm: RealGemmWorkspace::with_pool(n_loc, count.max(capacity), pool)?,
            projectors,
            coupling,
            count,
            rows,
            fft,
        })
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
        let pairs = columns.div_ceil(2);
        (
            2 * pairs + if self.max_items > 0 { 2 } else { 0 },
            usize::from(self.count > 0),
        )
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
        if mi > 0 {
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
        }

        // (2) per pair: pack, distributed local-potential apply, stage mirror h.
        for t in 0..pairs {
            let (ca, cb) = (2 * t, (2 * t + 1 < m).then_some(2 * t + 1));
            let get = |col: usize, slot: usize| {
                C::new(xs[col * n_loc + slot], xs[col * n_loc + slot + 1])
            };
            let recv = |col: usize, id: usize| {
                let off = (id * m + col) * 2;
                C::new(self.xrecv[off], self.xrecv[off + 1])
            };
            for (phi, &s) in self.cvec.iter_mut().zip(&self.src) {
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
            let out = &mut self.hbuf[t * local_len..(t + 1) * local_len];
            self.fft
                .apply_local(comm, &self.cvec, &self.potential, out)?;
            for &(k, id) in &self.in_items {
                let off = (id * pairs + t) * 2;
                self.hsend[off] = out[k].re;
                self.hsend[off + 1] = out[k].im;
            }
        }

        // (3) mirror h back to the representative owners.
        if mi > 0 {
            let blk = mi * 2 * pairs;
            comm.alltoall_f64(
                &self.hsend[..ranks * blk],
                blk,
                &mut self.hrecv[..ranks * blk],
            )?;
        }

        // Unpack, kinetic term.
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

        // Nonlocal: P^T X local + tree_sum, replicated coupling, P W local.
        let count = self.count;
        if count > 0 {
            if n_loc > 0 {
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
            if n_loc > 0 {
                self.gemm.multiply_into(
                    self.projectors.as_ref(),
                    RealOp::Identity,
                    self.weights.as_ref().subcols(0, m),
                    self.image.as_mut().subcols_mut(0, m),
                )?;
                for (y, v) in ys.iter_mut().zip(self.image.as_slice()) {
                    *y += v;
                }
            }
        }
        Ok(())
    }
}
