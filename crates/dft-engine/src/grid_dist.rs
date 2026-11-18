//! Distributed full-grid stages of the M2 path (ticket M2-GRID, 2026-11-18;
//! `BAANDI_M2_GRID=on|off`, default off, rounding class).
//!
//! Without the switch every rank of the M2 path repeats the same full-grid work
//! per SCF map: the Hartree FFT pair of the input density, the effective
//! potential sum, the Hartree FFT of the output density in the energy audit and
//! the audit's grid sums. None of it shrinks with the rank count. With the switch
//! each rank works on its own real-space planes and reciprocal sticks of a
//! slab-decomposed full-grid FFT (the M2a `DistributedPairFft` with every grid
//! point as an index: sticks along the first FFT axis, contiguous planes of the
//! first axis per rank, two transposes per transform):
//!
//! * potential: forward transform of the (replicated) input density on the
//!   owned planes, 4 pi / G^2 on the owned sticks, inverse transform, LDA XC on
//!   the owned plane points, effective potential ionic + (V_H + V_xc) on the
//!   planes; ONE allgatherv brings the effective potential (all planes) and the
//!   partial sums (E_H per stick, XC energy and integral per plane, error flags)
//!   to every rank. V_H and V_xc stay on the planes for the audit.
//! * energy audit: forward transform of the output density on the planes,
//!   E_H per stick; charge, double counting, the two Bregman linear terms, the
//!   local energy (real-space form w sum n V_ion, Parseval-equal to the
//!   reciprocal form) and the output XC energy per plane; ONE allgatherv.
//! * stop rule (`BAANDI_STOP_METRIC=kept`): the owned stick coefficients of the
//!   input and output densities are kept; the Coulomb estimate is a per-stick
//!   partial sum plus one small allgatherv.
//!
//! Determinism and numeric class. The transforms move data between ranks but
//! never change the arithmetic of a plane (2D lines) or of a stick (1D line), so
//! every transformed value is the same for every rank count. Every reduction is
//! a canonical sum: partial sums per plane (row-major inside the plane) or per
//! stick (ascending first-axis index), gathered and added in ascending plane or
//! stick order. The results therefore do not depend on the rank count (P >= 2),
//! and a fixed rank count is bit-reproducible. Against the replicated path
//! (`FftPlan3d`, row-major full-grid sums) the potential and the energies differ
//! by rounding only: rounding class, hence the switch.
//!
//! Every branch condition is replicated data (switch, functional, communicator,
//! grid) and every rank enters the same collectives. An invalid density on one
//! rank does not leave the collectives early: it travels as a flag in the
//! gathered sums and every rank returns the same error afterwards.
use crate::ScfError;
use dft_mpi::Collective;
use dft_physics::data::Array3;
use dft_physics::data::Complex64 as C;
use dft_physics::{PhysicsError, xc};
use dft_solvers::distributed_fft::DistributedPairFft;
use std::f64::consts::PI;
use std::sync::atomic::{AtomicU64, Ordering};

thread_local! {
    static OVERRIDE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}
/// Override `BAANDI_M2_GRID` for problems used on the calling thread (simulated
/// ranks in tests; every rank thread must set the same value).
pub fn set_override(on: Option<bool>) {
    OVERRIDE.with(|o| o.set(on));
}
/// True when the distributed grid stages are requested on this thread.
pub fn requested() -> bool {
    OVERRIDE
        .with(|o| o.get())
        .unwrap_or_else(|| dft_numerics::switches::on(dft_numerics::switches::Id::M2_GRID))
}

static COUNTS: [AtomicU64; 5] = [const { AtomicU64::new(0) }; 5];
const C_POTENTIAL: usize = 0;
const C_AUDIT: usize = 1;
const C_STOP: usize = 2;
const C_PLANS: usize = 3;
const C_DECLINED: usize = 4;
/// (distributed potential maps, distributed audits, distributed stop-rule
/// estimates, plans built, maps that declined: not eligible) since process
/// start, summed over the threads of the process (simulated ranks add up).
pub fn counters() -> [u64; 5] {
    std::array::from_fn(|i| COUNTS[i].load(Ordering::Relaxed))
}
fn bump(i: usize) {
    COUNTS[i].fetch_add(1, Ordering::Relaxed);
}
pub(crate) fn count_declined() {
    bump(C_DECLINED);
}

/// Pointwise LDA flavour of the distributed XC (the m2 path's eligibility).
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum Lda {
    Pz,
    Pw,
}

/// V_H and V_xc of the latest distributed input density on the owned planes.
pub(crate) struct PlaneFields {
    vh: Vec<f64>,
    vxc: Vec<f64>,
    hartree: f64,
    xc_energy: f64,
}

/// Terms of the energy audit (Ha) from the distributed output density.
pub(crate) struct AuditTerms {
    pub local: f64,
    pub charge: f64,
    pub hartree: f64,
    pub xc: f64,
    pub double_counting: f64,
    pub rh: f64,
    pub rx: f64,
}
impl AuditTerms {
    /// [local, charge, hartree, xc, double counting, rh, rx]
    /// (`energy::SharedAudit::from_distributed`).
    pub(crate) fn values(&self) -> [f64; 7] {
        [
            self.local,
            self.charge,
            self.hartree,
            self.xc,
            self.double_counting,
            self.rh,
            self.rx,
        ]
    }
}

pub(crate) struct GridPlan {
    fft: DistributedPairFft,
    ranks: usize,
    dims: [usize; 3],
    nz: usize,
    rows: usize,
    planes: usize,
    /// Row-major grid offset of every (plane, row) z-line of this rank.
    lines: Vec<usize>,
    /// 4 pi / |G|^2 at the local coefficients (0 at G = 0).
    kernel: Vec<f64>,
    /// Canonical stick slot (0..local_sticks) of every local coefficient.
    stick_of: Vec<u32>,
    local_sticks: usize,
    /// Allgathered layout: per rank [planes * PLANE_VALUES | sticks | flags].
    plane_counts: Vec<usize>,
    stick_counts: Vec<usize>,
    /// Position in the gathered stick block (all ranks concatenated) of the
    /// global stick with canonical index s (ascending (z, other)).
    stick_order: Vec<usize>,
    /// Coulomb-metric weights at the local coefficients (stop rule).
    stop_weights: Vec<f64>,
    lda: Lda,
    volume: f64,
    // Scratch (no allocation per map once sized).
    coef: Vec<C>,
    work: Vec<C>,
    plane_vals: Vec<f64>,
    stick_vals: Vec<f64>,
    send: Vec<f64>,
    recv: Vec<f64>,
    effective_planes: Vec<f64>,
    fields: Option<PlaneFields>,
    /// Owned stick coefficients of the latest input / output densities (the
    /// checked density, `to_reciprocal` normalization) and whether they hold
    /// an unclamped density (stop rule).
    input_hat: Vec<C>,
    output_hat: Vec<C>,
}

const POT_PLANE: usize = 2; // [xc energy, xc integral]
const AUDIT_PLANE: usize = 7; // [charge, dc, linear_h, linear_x, local, xc energy, xc integral]
const FLAGS: usize = 3; // [invalid density, clamped density, xc failure]

impl GridPlan {
    /// Plan for `dims` on this rank of `comm`. `stop_weights_full` is the
    /// Coulomb metric's full-grid weight array (row-major), if any.
    pub(crate) fn new<K: Collective>(
        comm: &K,
        dims: [usize; 3],
        reciprocal_columns: [[f64; 3]; 3],
        volume: f64,
        lda: Lda,
        stop_weights_full: Option<&[f64]>,
    ) -> Result<Self, ScfError> {
        let (rank, ranks) = (comm.rank(), comm.size());
        let total: usize = dims.iter().product();
        let mut indices = Vec::with_capacity(total);
        for i in 0..dims[0] {
            for j in 0..dims[1] {
                for k in 0..dims[2] {
                    indices.push([i, j, k]);
                }
            }
        }
        let err = |e: dft_solvers::distributed_fft::DistFftError| {
            ScfError::Numerical(format!("M2 grid plan: {e}"))
        };
        let fft = DistributedPairFft::new(dims, &indices, rank, ranks).map_err(err)?;
        let owners =
            dft_solvers::distributed_fft::position_owners(dims, &indices, ranks).map_err(err)?;
        drop(indices);
        let first = fft.first_axis();
        let n_other = dims[1 - first];
        let nz = dims[2];
        // Canonical stick id of a position: (z, other) in ascending order.
        let stick_id = |g: usize| {
            let (i, j, k) = (g / (dims[1] * nz), (g / nz) % dims[1], g % nz);
            k * n_other + if first == 0 { j } else { i }
        };
        let sticks_total = nz * n_other;
        let mut stick_owner = vec![usize::MAX; sticks_total];
        for (g, &o) in owners.iter().enumerate() {
            stick_owner[stick_id(g)] = o;
        }
        let mut stick_counts = vec![0usize; ranks];
        let mut slot_in_rank = vec![0usize; sticks_total];
        for (s, &o) in stick_owner.iter().enumerate() {
            slot_in_rank[s] = stick_counts[o];
            stick_counts[o] += 1;
        }
        let stick_offsets: Vec<usize> = std::iter::once(0)
            .chain(stick_counts.iter().scan(0, |a, &c| {
                *a += c;
                Some(*a)
            }))
            .collect();
        let stick_order = (0..sticks_total)
            .map(|s| stick_offsets[stick_owner[s]] + slot_in_rank[s])
            .collect();
        let stick_of = fft
            .local_positions()
            .iter()
            .map(|&g| slot_in_rank[stick_id(g)] as u32)
            .collect();
        let local_sticks = stick_counts[rank];
        // Kernel: the arithmetic of `hartree::squared_frequencies`.
        let b = reciprocal_columns;
        let kernel = fft
            .local_positions()
            .iter()
            .map(|&g| {
                let idx = [g / (dims[1] * nz), (g / nz) % dims[1], g % nz];
                if g == 0 {
                    return 0.0;
                }
                let m: [f64; 3] = std::array::from_fn(|a| {
                    if idx[a] <= (dims[a] - 1) / 2 {
                        idx[a] as f64
                    } else {
                        idx[a] as f64 - dims[a] as f64
                    }
                });
                let g2 = (0..3)
                    .map(|a| (0..3).map(|j| b[a][j] * m[j]).sum::<f64>().powi(2))
                    .sum::<f64>();
                4.0 * PI / g2
            })
            .collect::<Vec<_>>();
        if kernel.iter().any(|k| !k.is_finite()) {
            return Err(ScfError::NonFinite);
        }
        let stop_weights = match stop_weights_full {
            Some(w) if w.len() == total => fft.local_positions().iter().map(|&g| w[g]).collect(),
            _ => Vec::new(),
        };
        // The plan's plane rule: rank r starts at r * (n0 / P) + min(r, n0 % P).
        let n0 = dims[first];
        let plane_start = |r: usize| r * (n0 / ranks) + r.min(n0 % ranks);
        let plane_counts: Vec<usize> = (0..ranks)
            .map(|r| plane_start(r + 1) - plane_start(r))
            .collect();
        let planes = fft.plane_count();
        let rows = fft.rows();
        let mut lines = Vec::with_capacity(planes * rows);
        for i in 0..planes {
            for y in 0..rows {
                lines.push(fft.grid_start(i, y));
            }
        }
        let ll = fft.local_len();
        let points = planes * rows * nz;
        bump(C_PLANS);
        Ok(Self {
            ranks,
            dims,
            nz,
            rows,
            planes,
            lines,
            kernel,
            stick_of,
            local_sticks,
            plane_counts,
            stick_counts,
            stick_order,
            stop_weights,
            lda,
            volume,
            coef: vec![C::default(); ll],
            work: vec![C::default(); ll],
            plane_vals: Vec::new(),
            stick_vals: Vec::new(),
            send: Vec::new(),
            recv: Vec::new(),
            effective_planes: vec![0.0; points],
            fields: None,
            input_hat: Vec::new(),
            output_hat: Vec::new(),
            fft,
        })
    }

    pub(crate) fn dims(&self) -> [usize; 3] {
        self.dims
    }
    fn points(&self) -> usize {
        self.planes * self.rows * self.nz
    }
    fn weight(&self) -> f64 {
        self.volume / self.dims.iter().product::<usize>() as f64
    }

    /// Forward transform of the checked density `n` (replicated full grid) on
    /// the owned planes into `self.coef`; flags[0] invalid, flags[1] clamped.
    fn forward<K: Collective>(
        &mut self,
        comm: &K,
        n: &[f64],
        flags: &mut [f64; FLAGS],
    ) -> Result<(), ScfError> {
        let err = |e: dft_solvers::distributed_fft::DistFftError| {
            ScfError::Numerical(format!("M2 grid transform: {e}"))
        };
        let nz = self.nz;
        for i in 0..self.planes {
            let slab = self.fft.plane_mut();
            for y in 0..self.rows {
                let start = self.lines[i * self.rows + y];
                for (z, &v) in n[start..start + nz].iter().enumerate() {
                    // `checked_density`: invalid below -NEGATIVE_ROUNDOFF or nonfinite.
                    let c = match xc_checked(v) {
                        Some(c) => c,
                        None => {
                            flags[0] = 1.0;
                            0.0
                        }
                    };
                    if c.to_bits() != v.to_bits() {
                        flags[1] = 1.0;
                    }
                    slab[y * nz + z] = C::new(c, 0.0);
                }
            }
            self.fft.to_reciprocal_plane(i).map_err(err)?;
        }
        self.fft
            .finish_to_reciprocal(comm, &mut self.coef)
            .map_err(err)
    }

    /// Per-stick partial sums of `kernel |c|^2` into `self.stick_vals`
    /// (ascending first-axis order inside a stick: local coefficients are in
    /// ascending row-major position).
    fn stick_energy(&mut self, coef: &[C]) {
        self.stick_vals.clear();
        self.stick_vals.resize(self.local_sticks, 0.0);
        for ((z, &k), &s) in coef.iter().zip(&self.kernel).zip(&self.stick_of) {
            if k != 0.0 {
                self.stick_vals[s as usize] += k * z.norm_sqr();
            }
        }
    }

    /// One allgatherv of [prefix | plane_vals | stick_vals | flags] per rank;
    /// `prefix_counts[r]` leading values of rank r (planes of a field) come first.
    fn gather<K: Collective>(
        &mut self,
        comm: &K,
        prefix: Option<&[f64]>,
        per_plane: usize,
        flags: &[f64; FLAGS],
    ) -> Result<(), ScfError> {
        let plane_points = self.rows * self.nz;
        let with_prefix = prefix.is_some();
        self.send.clear();
        if let Some(p) = prefix {
            self.send.extend_from_slice(p);
        }
        self.send.extend_from_slice(&self.plane_vals);
        self.send.extend_from_slice(&self.stick_vals);
        self.send.extend_from_slice(flags);
        let counts: Vec<usize> = (0..self.ranks)
            .map(|r| {
                let planes = self.plane_counts[r];
                (if with_prefix {
                    planes * plane_points
                } else {
                    0
                }) + planes * per_plane
                    + self.stick_counts[r]
                    + FLAGS
            })
            .collect();
        debug_assert_eq!(self.send.len(), counts[comm_rank(comm)]);
        self.recv.resize(counts.iter().sum(), 0.0);
        comm.allgatherv_f64(&self.send, &counts, &mut self.recv)
            .map_err(|e| ScfError::Numerical(format!("M2 grid gather: {e}")))
    }

    /// Canonical sums of the gathered block: plane values (per_plane each, in
    /// ascending plane order), stick values (ascending canonical stick order)
    /// and the flags (max). `prefix` tells whether a field precedes them.
    fn reduce(
        &self,
        with_prefix: bool,
        per_plane: usize,
        plane_sums: &mut [f64],
    ) -> (f64, [f64; FLAGS]) {
        let plane_points = self.rows * self.nz;
        let mut offset = 0;
        let mut stick_base = Vec::with_capacity(self.ranks);
        let mut flags = [0.0; FLAGS];
        plane_sums.iter_mut().for_each(|v| *v = 0.0);
        // Ranks hold ascending contiguous plane ranges: rank order is plane order.
        for r in 0..self.ranks {
            let planes = self.plane_counts[r];
            if with_prefix {
                offset += planes * plane_points;
            }
            for p in 0..planes {
                for (q, s) in plane_sums.iter_mut().enumerate().take(per_plane) {
                    *s += self.recv[offset + p * per_plane + q];
                }
            }
            offset += planes * per_plane;
            stick_base.push(offset);
            offset += self.stick_counts[r];
            for (f, v) in flags.iter_mut().zip(&self.recv[offset..offset + FLAGS]) {
                *f = f64::max(*f, *v);
            }
            offset += FLAGS;
        }
        // Gathered stick position -> (rank, slot): stick_order indexes the
        // concatenation of the stick blocks; map it to the received buffer.
        let mut rank_of = Vec::with_capacity(self.ranks + 1);
        let mut acc = 0;
        for r in 0..self.ranks {
            rank_of.push(acc);
            acc += self.stick_counts[r];
        }
        rank_of.push(acc);
        let mut sticks = 0.0;
        for &pos in &self.stick_order {
            let r = rank_of.partition_point(|&a| a <= pos) - 1;
            sticks += self.recv[stick_base[r] + pos - rank_of[r]];
        }
        (sticks, flags)
    }

    /// LDA XC of the plane points of `n` into `vxc` (if given) and the per-plane
    /// [energy, integral] partial sums (unweighted) appended at `out[2p..]`.
    fn xc_planes(
        &self,
        n: &[f64],
        mut vxc: Option<&mut [f64]>,
        out: &mut [f64],
        stride: usize,
        at: usize,
        flags: &mut [f64; FLAGS],
    ) {
        let (rows, nz) = (self.rows, self.nz);
        let mut block = Vec::new();
        for i in 0..self.planes {
            let (mut e, mut s) = (0.0, 0.0);
            match self.lda {
                Lda::Pz => {
                    for y in 0..rows {
                        let start = self.lines[i * rows + y];
                        for (z, &v) in n[start..start + nz].iter().enumerate() {
                            match xc::lda_pz(v) {
                                Ok(p) => {
                                    let nn = v.max(0.0);
                                    let pot = p.v_xc();
                                    if let Some(vx) = vxc.as_deref_mut() {
                                        vx[(i * rows + y) * nz + z] = pot;
                                    }
                                    e += nn * p.eps_xc();
                                    s += nn * pot;
                                }
                                Err(_) => flags[2] = 1.0,
                            }
                        }
                    }
                }
                Lda::Pw => {
                    block.clear();
                    for y in 0..rows {
                        let start = self.lines[i * rows + y];
                        block.extend_from_slice(&n[start..start + nz]);
                    }
                    let len = block.len();
                    // Volume = point count makes the library weight exactly 1.
                    let arr = Array3::from_shape_vec((1, 1, len), block.clone())
                        .expect("plane block shape");
                    match dft_physics::spin::unpolarized_pw(&arr, len as f64) {
                        Ok(r) => {
                            if let Some(vx) = vxc.as_deref_mut() {
                                vx[i * rows * nz..(i + 1) * rows * nz]
                                    .copy_from_slice(r.potential.as_slice().unwrap());
                            }
                            e = r.energy.0;
                            s = r.potential_integral.0;
                        }
                        Err(_) => flags[2] = 1.0,
                    }
                }
            }
            out[i * stride + at] = e;
            out[i * stride + at + 1] = s;
        }
    }

    /// Distributed Hartree + XC fields and effective potential of the input
    /// density `n` (replicated, row-major). Returns (E_H, E_xc, integral
    /// n v_xc, effective potential on the full grid); V_H and V_xc stay on the
    /// planes for the audit, the owned stick coefficients for the stop rule.
    /// Collective.
    pub(crate) fn potential<K: Collective>(
        &mut self,
        comm: &K,
        n: &[f64],
        ionic: &[f64],
        effective: &mut [f64],
    ) -> Result<(f64, f64, f64), ScfError> {
        let start = dft_numerics::stage_clock::start();
        let (mut vh, mut vxc) = match self.fields.take() {
            Some(f) => (f.vh, f.vxc),
            None => (Vec::new(), Vec::new()),
        };
        self.input_hat.clear();
        let mut flags = [0.0; FLAGS];
        self.forward(comm, n, &mut flags)?;
        let coef = std::mem::take(&mut self.coef);
        self.stick_energy(&coef);
        // V_H(G) = 4 pi / G^2 n(G) on the owned sticks (G = 0 vanishes).
        for ((w, z), &k) in self.work.iter_mut().zip(&coef).zip(&self.kernel) {
            *w = if k != 0.0 { *z * k } else { C::default() };
        }
        self.input_hat.extend_from_slice(&coef);
        self.coef = coef;
        let err = |e: dft_solvers::distributed_fft::DistFftError| {
            ScfError::Numerical(format!("M2 grid transform: {e}"))
        };
        let work = std::mem::take(&mut self.work);
        let begun = self.fft.begin_to_real(comm, &work);
        self.work = work;
        begun.map_err(err)?;
        let points = self.points();
        let (rows, nz) = (self.rows, self.nz);
        vh.resize(points, 0.0);
        vxc.resize(points, 0.0);
        for i in 0..self.planes {
            self.fft.to_real_plane(i).map_err(err)?;
            let slab = self.fft.plane();
            for (o, z) in vh[i * rows * nz..(i + 1) * rows * nz].iter_mut().zip(slab) {
                *o = z.re;
            }
        }
        let mut vals = std::mem::take(&mut self.plane_vals);
        vals.clear();
        vals.resize(self.planes * POT_PLANE, 0.0);
        self.xc_planes(n, Some(&mut vxc), &mut vals, POT_PLANE, 0, &mut flags);
        self.plane_vals = vals;
        // ionic + (V_H + V_xc), the replicated path's expression per point.
        let mut eff = std::mem::take(&mut self.effective_planes);
        for i in 0..self.planes {
            for y in 0..rows {
                let start = self.lines[i * rows + y];
                let at = (i * rows + y) * nz;
                for z in 0..nz {
                    eff[at + z] = ionic[start + z] + (vh[at + z] + vxc[at + z]);
                }
            }
        }
        let gathered = self.gather(comm, Some(&eff), POT_PLANE, &flags);
        self.effective_planes = eff;
        gathered?;
        let mut sums = [0.0; POT_PLANE];
        let (stick_sum, flags) = self.reduce(true, POT_PLANE, &mut sums);
        if flags[0] != 0.0 {
            return Err(PhysicsError::Density(-1.0).into());
        }
        if flags[2] != 0.0 {
            return Err(ScfError::Numerical(
                "M2 grid XC: invalid density on a rank".into(),
            ));
        }
        // Planes of the concatenation in ascending first-axis order.
        self.scatter_planes(effective);
        let weight = self.weight();
        let hartree = 0.5 * self.volume * stick_sum;
        let (xe, xi) = (sums[0] * weight, sums[1] * weight);
        if !hartree.is_finite() || !xe.is_finite() || !xi.is_finite() {
            return Err(ScfError::NonFinite);
        }
        if flags[1] != 0.0 {
            self.input_hat.clear();
        }
        self.fields = Some(PlaneFields {
            vh,
            vxc,
            hartree,
            xc_energy: xe,
        });
        crate::m2::add(crate::m2::Stage::Xc, start);
        bump(C_POTENTIAL);
        Ok((hartree, xe, xi))
    }

    /// Copy the gathered plane field (prefix of every rank block) into the
    /// row-major full grid.
    fn scatter_planes(&self, full: &mut [f64]) {
        let plane_points = self.rows * self.nz;
        let nz = self.nz;
        let first = self.fft.first_axis();
        let mut offset = 0;
        let mut x = 0;
        for r in 0..self.ranks {
            let planes = self.plane_counts[r];
            for p in 0..planes {
                for y in 0..self.rows {
                    let start = if first == 0 {
                        ((x + p) * self.dims[1] + y) * nz
                    } else {
                        (y * self.dims[1] + x + p) * nz
                    };
                    let src = offset + (p * self.rows + y) * nz;
                    full[start..start + nz].copy_from_slice(&self.recv[src..src + nz]);
                }
            }
            x += planes;
            offset += planes * plane_points + planes * POT_PLANE + self.stick_counts[r] + FLAGS;
        }
    }

    /// True when `potential` stored the plane fields of the latest input.
    pub(crate) fn has_fields(&self) -> bool {
        self.fields.is_some()
    }

    /// Audit terms of the output density `out` against the input `input`
    /// whose fields `potential` computed last. Collective.
    pub(crate) fn audit<K: Collective>(
        &mut self,
        comm: &K,
        input: &[f64],
        out: &[f64],
        ionic: &[f64],
    ) -> Result<AuditTerms, ScfError> {
        let start = dft_numerics::stage_clock::start();
        let fields = self
            .fields
            .as_ref()
            .ok_or(ScfError::InvalidInput("M2 grid audit without fields"))?;
        let (e_h_in, e_xc_in) = (fields.hartree, fields.xc_energy);
        self.output_hat.clear();
        let mut flags = [0.0; FLAGS];
        self.forward(comm, out, &mut flags)?;
        let coef = std::mem::take(&mut self.coef);
        self.stick_energy(&coef);
        self.output_hat.extend_from_slice(&coef);
        self.coef = coef;
        let (rows, nz) = (self.rows, self.nz);
        let mut vals = std::mem::take(&mut self.plane_vals);
        vals.clear();
        vals.resize(self.planes * AUDIT_PLANE, 0.0);
        let fields = self.fields.as_ref().expect("checked above");
        for i in 0..self.planes {
            let (mut charge, mut dc, mut lh, mut lx, mut local) = (0.0, 0.0, 0.0, 0.0, 0.0);
            for y in 0..rows {
                let start = self.lines[i * rows + y];
                let at = (i * rows + y) * nz;
                for z in 0..nz {
                    let (o, a) = (out[start + z], input[start + z]);
                    let (h, x) = (fields.vh[at + z], fields.vxc[at + z]);
                    charge += o;
                    dc += o * (h + x);
                    lh += (o - a) * h;
                    lx += (o - a) * x;
                    local += o * ionic[start + z];
                }
            }
            vals[i * AUDIT_PLANE..i * AUDIT_PLANE + 5]
                .copy_from_slice(&[charge, dc, lh, lx, local]);
        }
        self.xc_planes(out, None, &mut vals, AUDIT_PLANE, 5, &mut flags);
        self.plane_vals = vals;
        let gathered = self.gather(comm, None, AUDIT_PLANE, &flags);
        gathered?;
        let mut sums = [0.0; AUDIT_PLANE];
        let (stick_sum, flags) = self.reduce(false, AUDIT_PLANE, &mut sums);
        if flags[0] != 0.0 {
            return Err(PhysicsError::Density(-1.0).into());
        }
        if flags[2] != 0.0 {
            return Err(ScfError::Numerical(
                "M2 grid XC: invalid density on a rank".into(),
            ));
        }
        if flags[1] != 0.0 {
            self.output_hat.clear();
        }
        let w = self.weight();
        let hartree = 0.5 * self.volume * stick_sum;
        let xc = sums[5] * w;
        let linear_h = w * sums[2];
        let linear_x = w * sums[3];
        let terms = AuditTerms {
            local: w * sums[4],
            charge: sums[0] * w,
            hartree,
            xc,
            double_counting: w * sums[1],
            rh: hartree - e_h_in - linear_h,
            rx: xc - e_xc_in - linear_x,
        };
        if [
            terms.local,
            terms.charge,
            terms.hartree,
            terms.xc,
            terms.double_counting,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(ScfError::NonFinite);
        }
        crate::m2::add(crate::m2::Stage::Xc, start);
        bump(C_AUDIT);
        Ok(terms)
    }

    /// True when both kept coefficient sets belong to unclamped densities of
    /// the latest potential / audit (replicated after the gathers).
    pub(crate) fn stop_ready(&self) -> bool {
        !self.stop_weights.is_empty()
            && self.input_hat.len() == self.kernel.len()
            && self.output_hat.len() == self.kernel.len()
    }

    /// Coulomb estimate N sum_G w(G) |out(G) - in(G)|^2 from the kept owned
    /// coefficients, canonical stick order. Collective.
    pub(crate) fn stop_error<K: Collective>(&mut self, comm: &K) -> Result<f64, ScfError> {
        self.stick_vals.clear();
        self.stick_vals.resize(self.local_sticks, 0.0);
        for (((a, b), &w), &s) in self
            .input_hat
            .iter()
            .zip(&self.output_hat)
            .zip(&self.stop_weights)
            .zip(&self.stick_of)
        {
            if w != 0.0 {
                self.stick_vals[s as usize] += w * (b - a).norm_sqr();
            }
        }
        self.plane_vals.clear();
        self.plane_vals.resize(self.planes, 0.0);
        let flags = [0.0; FLAGS];
        self.gather(comm, None, 1, &flags)?;
        let mut unused = [0.0; 1];
        let (sum, _) = self.reduce(false, 1, &mut unused);
        let energy = sum * self.dims.iter().product::<usize>() as f64;
        if !energy.is_finite() {
            return Err(ScfError::NonFinite);
        }
        bump(C_STOP);
        Ok(energy.max(0.0))
    }
}

fn comm_rank<K: Collective>(comm: &K) -> usize {
    comm.rank()
}

/// `xc::checked_density` without the error payload: None for a nonfinite or
/// too negative value, otherwise max(n, 0).
fn xc_checked(n: f64) -> Option<f64> {
    if !n.is_finite() || n < -xc::NEGATIVE_ROUNDOFF {
        None
    } else {
        Some(n.max(0.0))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use dft_geometry::Cell;
    use dft_mpi::{SimRank, simulate};
    use dft_physics::data::{Fft3d, FftPlan3d};
    use dft_physics::hartree;

    fn lcg(s: &mut u64) -> f64 {
        *s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (*s >> 11) as f64 / (1u64 << 53) as f64
    }
    fn field(dims: [usize; 3], seed: u64, base: f64) -> Array3<f64> {
        let mut s = seed;
        Array3::from_shape_fn(dims, |_| base + 0.05 * lcg(&mut s))
    }
    fn rel(a: f64, b: f64) -> f64 {
        (a - b).abs() / b.abs().max(1e-300)
    }

    #[derive(PartialEq, Debug)]
    struct Bits {
        effective: Vec<u64>,
        energies: [u64; 3],
        audit: [u64; 7],
        stop: u64,
    }
    type RankOut = (Vec<f64>, (f64, f64, f64), AuditTerms, f64);
    fn bits(r: &RankOut) -> Bits {
        Bits {
            effective: r.0.iter().map(|v| v.to_bits()).collect(),
            energies: [r.1.0.to_bits(), r.1.1.to_bits(), r.1.2.to_bits()],
            audit: [
                r.2.local,
                r.2.charge,
                r.2.hartree,
                r.2.xc,
                r.2.double_counting,
                r.2.rh,
                r.2.rx,
            ]
            .map(f64::to_bits),
            stop: r.3.to_bits(),
        }
    }

    /// Replicated reference (serial physics routines, row-major full-grid
    /// sums) against the distributed stages on `ranks` simulated ranks; returns
    /// the bits of rank 0 after checking that every rank holds the same bits.
    fn run(cell: &Cell, dims: [usize; 3], ranks: usize, lda: Lda) -> Bits {
        let volume = cell.volume_bohr3();
        let n_in = field(dims, 11, 0.02);
        let n_out = field(dims, 23, 0.02);
        let ionic = field(dims, 37, -0.4);
        let total: usize = dims.iter().product();
        let mut s = 5;
        let weights: Vec<f64> = (0..total)
            .map(|g| if g == 0 { 0.0 } else { lcg(&mut s) })
            .collect();
        let mut fft = FftPlan3d::new(dims).unwrap();
        let h_in = hartree::evaluate(&n_in, cell, &mut fft).unwrap();
        let h_out = hartree::evaluate(&n_out, cell, &mut fft).unwrap();
        let lda_eval = |n: &Array3<f64>| match lda {
            Lda::Pz => xc::evaluate(n, volume).unwrap(),
            Lda::Pw => dft_physics::spin::unpolarized_pw(n, volume).unwrap(),
        };
        let (x_in, x_out) = (lda_eval(&n_in), lda_eval(&n_out));
        let w = volume / total as f64;
        let effective_ref: Vec<f64> = ionic
            .iter()
            .zip(&h_in.potential)
            .zip(&x_in.potential)
            .map(|((i, h), x)| i + (h.re + x))
            .collect();
        let mut sums = [0.0; 5];
        for ((((o, a), h), x), v) in n_out
            .iter()
            .zip(&n_in)
            .zip(&h_in.potential)
            .zip(&x_in.potential)
            .zip(&ionic)
        {
            sums[0] += o;
            sums[1] += o * (h.re + x);
            sums[2] += (o - a) * h.re;
            sums[3] += (o - a) * x;
            sums[4] += o * v;
        }
        let mut hat_in = Array3::zeros(dims);
        let mut hat_out = Array3::zeros(dims);
        fft.to_reciprocal(&n_in.mapv(|v| C::new(v, 0.0)), &mut hat_in)
            .unwrap();
        fft.to_reciprocal(&n_out.mapv(|v| C::new(v, 0.0)), &mut hat_out)
            .unwrap();
        let stop_ref = total as f64
            * weights
                .iter()
                .zip(hat_in.iter().zip(&hat_out))
                .map(|(w, (a, b))| w * (b - a).norm_sqr())
                .sum::<f64>();
        let runs: Vec<RankOut> = simulate(ranks, |comm: &SimRank| {
            let mut plan = GridPlan::new(
                comm,
                dims,
                cell.reciprocal_columns_inverse_bohr(),
                volume,
                lda,
                Some(&weights),
            )
            .unwrap();
            let mut eff = vec![0.0; total];
            let e = plan
                .potential(
                    comm,
                    n_in.as_slice().unwrap(),
                    ionic.as_slice().unwrap(),
                    &mut eff,
                )
                .unwrap();
            let t = plan
                .audit(
                    comm,
                    n_in.as_slice().unwrap(),
                    n_out.as_slice().unwrap(),
                    ionic.as_slice().unwrap(),
                )
                .unwrap();
            assert!(plan.stop_ready());
            let stop = plan.stop_error(comm).unwrap();
            (eff, e, t, stop)
        });
        let scale = effective_ref.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        for (eff, e, t, stop) in &runs {
            let worst = eff
                .iter()
                .zip(&effective_ref)
                .fold(0.0_f64, |m, (a, b)| m.max((a - b).abs()));
            assert!(worst <= 1e-13 * scale, "effective potential {worst:e}");
            assert!(rel(e.0, h_in.energy.0) < 1e-12, "E_H in");
            assert!(rel(e.1, x_in.energy.0) < 1e-12, "E_xc in");
            assert!(rel(e.2, x_in.potential_integral.0) < 1e-12, "int n v_xc");
            assert!(rel(t.hartree, h_out.energy.0) < 1e-12, "E_H out");
            assert!(rel(t.xc, x_out.energy.0) < 1e-12, "E_xc out");
            assert!(rel(t.charge, w * sums[0]) < 1e-13, "charge");
            assert!(rel(t.double_counting, w * sums[1]) < 1e-12, "dc");
            assert!(rel(t.local, w * sums[4]) < 1e-12, "local");
            let rh = h_out.energy.0 - h_in.energy.0 - w * sums[2];
            let rx = x_out.energy.0 - x_in.energy.0 - w * sums[3];
            assert!((t.rh - rh).abs() < 1e-12 * h_out.energy.0.abs(), "rh");
            assert!((t.rx - rx).abs() < 1e-12 * x_out.energy.0.abs(), "rx");
            assert!(rel(*stop, stop_ref) < 1e-12, "stop {stop} {stop_ref}");
        }
        let first = bits(&runs[0]);
        for r in &runs[1..] {
            assert!(bits(r) == first, "ranks hold different bits");
        }
        first
    }

    #[test]
    fn distributed_grid_stages_match_the_replicated_path_for_every_rank_count() {
        let cubic = Cell::from_user_rows_bohr([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]]).unwrap();
        let skew =
            Cell::from_user_rows_bohr([[6.2, 0., 0.], [1.1, 5.7, 0.], [0.4, 0.7, 6.6]]).unwrap();
        for (cell, dims, lda) in [
            (&cubic, [12, 12, 12], Lda::Pz),
            (&skew, [10, 12, 9], Lda::Pz),
            (&skew, [9, 8, 11], Lda::Pw),
        ] {
            let reference = run(cell, dims, 1, lda);
            for ranks in 2..=7 {
                // Rank-count independent bits: canonical plane/stick sums and
                // transforms whose per-line arithmetic ignores the partition.
                assert!(
                    run(cell, dims, ranks, lda) == reference,
                    "{dims:?} P={ranks}"
                );
            }
        }
    }

    #[test]
    fn invalid_density_fails_on_every_rank_after_the_collectives() {
        let cell = Cell::from_user_rows_bohr([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]]).unwrap();
        let dims = [8, 8, 8];
        let mut n = field(dims, 3, 0.02);
        n[[7, 3, 1]] = -1.0; // on the last rank's planes only
        let ionic = field(dims, 4, -0.4);
        let errors = simulate(3, |comm: &SimRank| {
            let mut plan = GridPlan::new(
                comm,
                dims,
                cell.reciprocal_columns_inverse_bohr(),
                cell.volume_bohr3(),
                Lda::Pz,
                None,
            )
            .unwrap();
            let mut eff = vec![0.0; 512];
            plan.potential(
                comm,
                n.as_slice().unwrap(),
                ionic.as_slice().unwrap(),
                &mut eff,
            )
            .is_err()
        });
        assert_eq!(errors, vec![true; 3]);
    }
}
