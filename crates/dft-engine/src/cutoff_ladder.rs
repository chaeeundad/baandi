//! C1 cutoff ladder (design note docs/roadmap/c1-cutoff-ladder.md): the first
//! part of an SCF runs at a reduced cutoff E1 = fraction * E2 on the smaller FFT
//! grid that E1 implies; the state is then prolonged EXACTLY into the full basis
//! and the unchanged strict session finishes at the full cutoff E2.
//!
//! Exactness of the prolongation: plane-wave spheres are nested, so every coarse
//! coefficient c(G) is copied to the same Miller index of the fine basis and the
//! new G vectors get zero (an injection: orthonormality and the Hermitian Gamma
//! symmetry c(-G) = conj c(G) are preserved bit for bit). The density is a
//! trigonometric polynomial on the coarse density sphere; its Fourier
//! coefficients are copied to the fine grid and zero padded (Fourier
//! interpolation, exact for a band-limited field; the G = 0 term and so the
//! charge are unchanged).
//!
//! The result of the calculation is the full-cutoff session's: its energy audit,
//! eigenpairs and fixed-point certificate are all at the full cutoff; the coarse
//! stage supplies only the start state (density and warm orbitals). By default
//! the mixer history starts empty at the full cutoff: absolute (input, residual)
//! pairs do not translate between the two maps, whose fixed points differ.
//!
//! S2 (`secant_pairs > 0`, `BAANDI_LADDER_SECANT`, default off): the DIFFERENCE
//! pairs of the coarse Pulay history do translate to first order (dr_i is the
//! response of the residual map to dx_i, dominated by the low-G dielectric
//! response that the coarse grid holds). The last h of them are prolonged with
//! the same exact Fourier embedding as the density and seeded into the fine
//! mixer as a closed block of differences (never differenced against a fine
//! pair). The first fine step is a plain evaluation at the prolonged density;
//! its mix uses the seeded least squares with the unchanged spectral cut, and
//! the seeds age out first as fine differences arrive.
//!
//! Scope: one exact Gamma point on the real (ORC) path, fixed uniform integer
//! occupations, no symmetry projector, no frozen-Miller strain reference. The
//! strict profile always; the standard profile only with `standard = true`
//! (`BAANDI_CUTOFF_LADDER_STANDARD=on`, S1, default off): its coarse stage runs
//! the standard profile at E1 (uncertified, stopped by the same residual rule)
//! and the full-cutoff standard session is unchanged except for the first-step
//! eigen tolerance (`StandardSeed`). Everything else runs the ordinary SCF with
//! a recorded reason.
//!
//! General rule for the coarse level (no per-system constant): E1 is a stated
//! fraction of the full cutoff (default 0.5, the design's pre-registered arm);
//! the coarse stage stops when its density residual falls to a fixed level
//! (default 1e-4 Bohr^-3/2, pre-registered in the design note) or after a fixed
//! iteration budget. Neither value is tuned on any system.
//!
//! D1 hand-off knobs (work order 2026-11-14, all default to the adopted ladder):
//! `BAANDI_LADDER_COARSE_STOP` (the fixed level), `BAANDI_LADDER_COARSE_RULE`
//! (`fixed`, or `atomic`: stop at the larger of the fixed level and the
//! truncation distance of the atomic start density, `CoarseStopRule`) and
//! `BAANDI_LADDER_FIRST_TOL` (first fine eigen tolerance, `FirstTolerance`).
//! The fraction is `BAANDI_CUTOFF_LADDER=<fraction>` (CLI). Every knob changes
//! only the start state or the first-step tolerance of the full-cutoff session;
//! its stopping rule, certificate and verification are untouched.
use super::*;
use crate::checkpoint::{ElectronicState, Orbitals};
use dft_mpi::Collective;
use dft_numerics::real::RealMatrix;
use dft_solvers::SolverError;
use std::collections::HashMap;

/// Exclusive handoff diagnostics. Enabled only by `BAANDI_PHASE_PROFILE`.
/// Process sums are exported through the existing M2 counters; `handoff_profile`
/// returns this rank's measurements, avoiding sums/maxima of unrelated ranks.
#[derive(Clone, Copy, Debug, Default)]
pub struct HandoffProfile {
    /// Exclusive seconds in `HANDOFF_PART_LABELS` order; sum is handoff wall
    /// on this rank (not the sum of independently selected rank maxima).
    pub seconds: [f64; 10],
    pub collectives: [u64; 5],
    pub materializations: u64,
    /// Delivered orbital gather payload (not MPI wire traffic).
    pub gathered_bytes: u64,
    /// Full complex coefficient storage per rank (coarse plus fine).
    pub full_matrix_bytes: u64,
    /// Upper bound on live coefficient-buffer payload, excluding allocator,
    /// density, solver and immutable operator storage. Not process RSS.
    pub coefficient_peak_bound_bytes: u64,
    /// Collective wall inside the materialization span (subset of seconds[0]).
    pub materialize_collective_seconds: f64,
    /// transfers, bytes sent including transport padding, populated bytes,
    /// ambiguous-check fallbacks, coefficient peak bytes (per handoff).
    pub owned: [u64; 5],
}
impl HandoffProfile {
    pub fn seconds(&self) -> f64 {
        self.seconds.iter().sum()
    }
}
thread_local! {
    static HANDOFF: std::cell::RefCell<HandoffProfile> = Default::default();
    static HANDOFF_CHILD_NS: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
}
/// Cumulative diagnostics on the calling rank; no additional collectives.
pub fn handoff_profile() -> HandoffProfile {
    HANDOFF.with(|h| *h.borrow())
}

#[derive(Clone, Copy)]
pub(crate) enum HandoffPart {
    Materialize,
    OrbitalCopy,
    DensityCopy,
    Encode,
    Setup,
    Import,
    Check,
    ExportDensity,
    ExportCopy,
    Transfer,
}
pub const HANDOFF_PART_LABELS: [&str; 10] = [
    "handoff.materialize",
    "handoff.orbital_copy",
    "handoff.density_copy",
    "handoff.encode",
    "handoff.setup",
    "handoff.import",
    "handoff.check",
    "handoff.export_density",
    "handoff.export_copy",
    "handoff.transfer",
];
use std::sync::atomic::{AtomicU64, Ordering::Relaxed};
static HANDOFF_NS: [AtomicU64; 10] = [const { AtomicU64::new(0) }; 10];
static HANDOFF_COUNTS: [AtomicU64; 9] = [const { AtomicU64::new(0) }; 9];
static HANDOFF_GATHER_NS: AtomicU64 = AtomicU64::new(0);
const HANDOFF_COUNT_NAMES: [&str; 9] = [
    "handoff.allgatherv",
    "handoff.tree_sum",
    "handoff.bcast",
    "handoff.alltoall",
    "handoff.barrier",
    "handoff.materializations",
    "handoff.gathered_bytes",
    "handoff.full_matrix_bytes",
    "handoff.coefficient_peak_bound_bytes",
];
static OWNED_HANDOFF: [AtomicU64; 5] = [const { AtomicU64::new(0) }; 5];
const OWNED_HANDOFF_NAMES: [&str; 5] = [
    "handoff.owned_transfers",
    "handoff.sent_bytes_with_padding",
    "handoff.populated_payload_bytes",
    "handoff.check_fallbacks",
    "handoff.owned_peak_bytes_sum",
];
pub(crate) fn handoff_owned_count(i: usize, n: usize) {
    if crate::profile::enabled() {
        OWNED_HANDOFF[i].fetch_add(n as u64, Relaxed);
        HANDOFF.with(|h| h.borrow_mut().owned[i] += n as u64);
    }
}
pub(crate) fn handoff_owned_peak(bytes: usize) {
    if !crate::profile::enabled() {
        return;
    }
    HANDOFF.with(|h| {
        let mut h = h.borrow_mut();
        let old = h.owned[4];
        h.owned[4] = old.max(bytes as u64);
        OWNED_HANDOFF[4].fetch_add(h.owned[4] - old, Relaxed);
    });
}

pub(crate) fn handoff_stats() -> dft_physics::distributed_real_gamma::M2Counters {
    (
        HANDOFF_COUNT_NAMES
            .iter()
            .zip(&HANDOFF_COUNTS)
            .map(|(&n, v)| (n, v.load(Relaxed) as usize))
            .chain(
                OWNED_HANDOFF_NAMES
                    .iter()
                    .zip(&OWNED_HANDOFF)
                    .map(|(&n, v)| (n, v.load(Relaxed) as usize)),
            )
            .collect(),
        HANDOFF_PART_LABELS
            .iter()
            .zip(&HANDOFF_NS)
            .map(|(&n, v)| (n, v.load(Relaxed) as f64 * 1e-9))
            .chain(std::iter::once((
                "handoff.materialize_collective",
                HANDOFF_GATHER_NS.load(Relaxed) as f64 * 1e-9,
            )))
            .collect(),
    )
}

/// Nested spans subtract child time. Collective deltas are recorded only on
/// the outer span, so import/check collectives are not counted twice.
pub(crate) struct HandoffSpan {
    part: HandoffPart,
    start: Instant,
    child_ns: u64,
    before: dft_numerics::stage_clock::local::Snapshot,
    parent: bool,
}
thread_local! { static HANDOFF_DEPTH: std::cell::Cell<usize> = const { std::cell::Cell::new(0) }; }
impl HandoffSpan {
    pub(crate) fn check() -> Option<Self> {
        if !crate::profile::enabled() || HANDOFF_DEPTH.with(|d| d.get()) == 0 {
            return None;
        }
        Self::enter(HandoffPart::Check)
    }
    pub(crate) fn enter(part: HandoffPart) -> Option<Self> {
        if !crate::profile::enabled() {
            return None;
        }
        let parent = HANDOFF_DEPTH.with(|d| {
            let v = d.get();
            d.set(v + 1);
            v > 0
        });
        Some(Self {
            part,
            start: Instant::now(),
            child_ns: HANDOFF_CHILD_NS.with(|n| n.get()),
            before: dft_numerics::stage_clock::local::snapshot(),
            parent,
        })
    }
}
impl Drop for HandoffSpan {
    fn drop(&mut self) {
        let ns = self.start.elapsed().as_nanos() as u64;
        let exclusive = HANDOFF_CHILD_NS.with(|n| {
            let children = n.get() - self.child_ns;
            n.set(self.child_ns + ns);
            ns.saturating_sub(children)
        });
        let part = self.part as usize;
        HANDOFF_NS[part].fetch_add(exclusive, Relaxed);
        HANDOFF.with(|h| h.borrow_mut().seconds[part] += exclusive as f64 * 1e-9);
        HANDOFF_DEPTH.with(|d| d.set(d.get() - 1));
        let delta = dft_numerics::stage_clock::local::snapshot().since(&self.before);
        if matches!(self.part, HandoffPart::Materialize) {
            let offset = dft_numerics::stage_clock::local::COLLECTIVE;
            let ns: u64 = delta.nanos[offset..offset + 5].iter().sum();
            HANDOFF_GATHER_NS.fetch_add(ns, Relaxed);
            HANDOFF.with(|h| h.borrow_mut().materialize_collective_seconds += ns as f64 * 1e-9);
        }
        if !self.parent {
            HANDOFF.with(|h| {
                let mut h = h.borrow_mut();
                for (i, count) in h.collectives.iter_mut().enumerate() {
                    let n = delta.calls[dft_numerics::stage_clock::local::COLLECTIVE + i];
                    *count += n;
                    HANDOFF_COUNTS[i].fetch_add(n, Relaxed);
                }
            });
        }
    }
}

fn handoff_storage(coarse: usize, fine: usize, bands: usize) {
    if !crate::profile::enabled() {
        return;
    }
    let bytes = (16 * (coarse + fine) * bands) as u64;
    HANDOFF.with(|h| h.borrow_mut().full_matrix_bytes += bytes);
    HANDOFF_COUNTS[7].fetch_add(bytes, Relaxed);
    // Lazy unpack: gathered real block + cached and cloned complex blocks.
    // Prolong/encode: coarse complex + two fine complex-size buffers.
    let peak = (40 * coarse * bands).max(16 * coarse * bands + 32 * fine * bands) as u64;
    HANDOFF.with(|h| {
        let mut h = h.borrow_mut();
        let previous = h.coefficient_peak_bound_bytes;
        h.coefficient_peak_bound_bytes = previous.max(peak);
        // Process sum of rank-local high-water bounds (the simulated harness
        // divides by P). Repeated ladders on a rank do not inflate its peak.
        HANDOFF_COUNTS[8].fetch_add(h.coefficient_peak_bound_bytes - previous, Relaxed);
    });
}
/// Called only at an actual lazy materialization inside a handoff span.
pub(crate) fn handoff_materialized(rows: usize, bands: usize, ranks: usize) {
    if !crate::profile::enabled() || HANDOFF_DEPTH.with(|d| d.get()) == 0 {
        return;
    }
    let bytes = (8 * (ranks + rows + rows * bands)) as u64;
    HANDOFF.with(|h| {
        let mut h = h.borrow_mut();
        h.materializations += 1;
        h.gathered_bytes += bytes;
    });
    HANDOFF_COUNTS[5].fetch_add(1, Relaxed);
    HANDOFF_COUNTS[6].fetch_add(bytes, Relaxed);
}

// The transfer and import policy belong to the ladder. Gamma-real exposes
// only a read-only view and a constructor for its existing owned handle.
impl crate::gamma_real::LazyOrbitals {
    /// Agree a local validation outcome before entering any subsequent
    /// collective. No rank returns early on local shape/nonfinite errors.
    pub(crate) fn handoff_agree(&self, valid: bool) -> Result<(), ScfError> {
        let mut bad = [0.];
        self.handoff_view()
            .comm
            .tree_sum_f64(&[f64::from(!valid)], &mut bad)
            .map_err(|e| ScfError::Numerical(format!("owned handoff agreement: {e}")))?;
        if bad[0] != 0. {
            return Err(ScfError::InvalidInput(
                "owned handoff: agreed validation failure",
            ));
        }
        Ok(())
    }

    /// Validate overlaps in distributed storage using the complex unpack
    /// convention, not raw ORC norms. The threshold remains 1e-8. Close to
    /// that boundary use the existing sequential complex-row oracle so a
    /// changed reduction order cannot change admission.
    pub(crate) fn handoff_orthonormal(&self) -> Result<bool, ScfError> {
        let inner = self.handoff_view();
        self.handoff_agree(
            inner.rows.len() == inner.x.rows() && inner.x.as_slice().iter().all(|v| v.is_finite()),
        )?;
        let b = self.bands();
        let mut local = vec![0.; b * b];
        let s = std::f64::consts::SQRT_2;
        for j in 0..b {
            for i in 0..=j {
                let mut sum = 0.;
                for (l, &row) in inner.rows.iter().enumerate() {
                    let a = inner.x.as_slice()[i * inner.x.rows() + l];
                    let c = inner.x.as_slice()[j * inner.x.rows() + l];
                    sum += if row == 0 {
                        a * c
                    } else {
                        2. * ((a / s) * (c / s))
                    };
                }
                local[j * b + i] = sum;
            }
        }
        crate::cutoff_ladder::handoff_owned_peak(8 * (inner.x.as_slice().len() + 2 * local.len()));
        let mut gram = vec![0.; local.len()];
        inner
            .comm
            .tree_sum_f64(&local, &mut gram)
            .map_err(|e| ScfError::Numerical(format!("owned handoff overlap: {e}")))?;
        let mut ambiguous = false;
        let mut valid = gram.iter().all(|v| v.is_finite());
        for j in 0..b {
            for i in 0..=j {
                let error = (gram[j * b + i] - if i == j { 1. } else { 0. }).abs();
                // Conservative dot-product rounding envelope, scaled by
                // Cauchy-Schwarz norms; no admission tolerance is widened.
                let bound = 64.
                    * f64::EPSILON
                    * self.global_len() as f64
                    * (gram[i * b + i].abs() * gram[j * b + j].abs())
                        .sqrt()
                        .max(1.);
                ambiguous |= (error - 1e-8).abs() <= bound;
                valid &= error <= 1e-8;
            }
        }
        if ambiguous && gram.iter().all(|v| v.is_finite()) {
            crate::cutoff_ladder::handoff_owned_count(3, 1);
            let full = self.materialize()?;
            return crate::scf::session::orthonormal_split(full, Some(inner.comm))
                .map_err(|e| ScfError::Numerical(format!("owned handoff overlap oracle: {e}")));
        }
        Ok(valid)
    }

    /// One precomputed coarse/fine permutation and one equal-block all-to-all.
    /// Only coarse rows are payload; absent fine rows are zeroed locally.
    /// The existing transport requires padding to the largest peer message;
    /// its full size (including padding) is recorded, never hidden as payload.
    pub(crate) fn prolong_owned(
        &self,
        coarse: &PlaneWaveBasis,
        fine: &PlaneWaveBasis,
        grid: &dft_basis::FftGridSpec,
    ) -> Result<Self, ScfError> {
        use std::collections::HashMap;
        let _span =
            crate::cutoff_ladder::HandoffSpan::enter(crate::cutoff_ladder::HandoffPart::Transfer);
        let inner = self.handoff_view();
        let (p, rank, b, nc) = (
            inner.comm.size(),
            inner.comm.rank(),
            self.bands(),
            self.global_len(),
        );
        let err = |e: dft_mpi::MpiError| ScfError::Numerical(format!("owned handoff: {e}"));
        let half = dft_basis::HalfSphere::new(fine, grid);
        self.handoff_agree(
            half.is_ok()
                && coarse.vectors().len() == nc
                && inner.x.rows() == inner.rows.len()
                && b > 0
                && inner.rows.iter().all(|&r| r < nc)
                && inner.x.as_slice().iter().all(|v| v.is_finite()),
        )?;
        let half = half.map_err(|e| ScfError::Numerical(e.to_string()))?;
        // All dimensions and integer layout identity agree before variable
        // messages. Hash words are encoded exactly as two u32-valued f64s.
        let mut hash = 0xcbf29ce484222325u64;
        for word in coarse
            .vectors()
            .iter()
            .chain(fine.vectors())
            .flat_map(|g| g.miller)
            .map(|v| v as u32 as u64)
            .chain(grid.dims().into_iter().map(|v| v as u64))
        {
            hash = (hash ^ word).wrapping_mul(0x100000001b3);
        }
        let header = [
            nc as f64,
            half.full_len() as f64,
            b as f64,
            (hash as u32) as f64,
            (hash >> 32) as f64,
        ];
        let mut headers = vec![0.; p * header.len()];
        inner
            .comm
            .allgatherv_f64(&header, &vec![header.len(); p], &mut headers)
            .map_err(err)?;
        self.handoff_agree(headers.chunks(header.len()).all(|h| h == header))?;
        let mut lengths = vec![0.; p];
        inner
            .comm
            .allgatherv_f64(&[inner.rows.len() as f64], &vec![1; p], &mut lengths)
            .map_err(err)?;
        let counts: Vec<usize> = lengths.iter().map(|&n| n as usize).collect();
        self.handoff_agree(counts.iter().sum::<usize>() == nc)?;
        let mut global_rows = vec![0.; nc];
        let rows: Vec<f64> = inner.rows.iter().map(|&r| r as f64).collect();
        inner
            .comm
            .allgatherv_f64(&rows, &counts, &mut global_rows)
            .map_err(err)?;
        let mut source = vec![usize::MAX; nc];
        let mut offset = 0;
        let mut valid = true;
        for (r, &n) in counts.iter().enumerate() {
            for &g in &global_rows[offset..offset + n] {
                let i = g as usize;
                if i >= nc || source[i] != usize::MAX {
                    valid = false;
                } else {
                    source[i] = r;
                }
            }
            offset += n;
        }
        let owners = dft_solvers::distributed_fft::position_owners(grid.dims(), grid.g_to_fft(), p);
        self.handoff_agree(valid && owners.is_ok())?;
        let owners = owners.map_err(|e| ScfError::Numerical(e.to_string()))?;
        let fine_index: HashMap<_, _> = half
            .plus()
            .iter()
            .enumerate()
            .map(|(i, &g)| (fine.vectors()[g as usize].miller, i))
            .collect();
        let mut mapping = vec![0; nc];
        let mut populated = vec![false; half.plus().len()];
        for (i, (&plus, &minus)) in inner.half.plus().iter().zip(inner.half.minus()).enumerate() {
            let g = coarse.vectors()[plus as usize].miller;
            let Some(&f) = fine_index.get(&g) else {
                valid = false;
                continue;
            };
            valid &= coarse.vectors()[minus as usize].miller
                == fine.vectors()[half.minus()[f] as usize].miller
                && (i == 0) == (f == 0);
            populated[f] = true;
            if i == 0 {
                mapping[0] = 0;
            } else {
                mapping[2 * i - 1] = 2 * f - 1;
                mapping[2 * i] = 2 * f;
            }
        }
        self.handoff_agree(valid)?;
        let nf = half.full_len();
        let target: Vec<usize> = (0..nf)
            .map(|r| owners[half.plus()[r.div_ceil(2)] as usize])
            .collect();
        let fine_rows: Vec<usize> = (0..nf).filter(|&r| target[r] == rank).collect();
        let mut slot = vec![usize::MAX; nf];
        for (i, &r) in fine_rows.iter().enumerate() {
            slot[r] = i;
        }
        // Both sides enumerate edges in global coarse ORC order, independent
        // of source storage order. This is the one pair-of-layouts plan.
        let mut edges = vec![Vec::new(); p * p];
        for (r, &f) in mapping.iter().enumerate() {
            edges[source[r] * p + target[f]].push(r);
        }
        let block = edges.iter().map(Vec::len).max().unwrap_or(0) * b;
        let mut send = vec![0.; p * block];
        let mut recv = vec![0.; p * block];
        let mut coarse_slot = vec![usize::MAX; nc];
        for (l, &r) in inner.rows.iter().enumerate() {
            coarse_slot[r] = l;
        }
        for dest in 0..p {
            for (i, &r) in edges[rank * p + dest].iter().enumerate() {
                for j in 0..b {
                    send[dest * block + i * b + j] =
                        inner.x.as_slice()[j * inner.x.rows() + coarse_slot[r]];
                }
            }
        }
        inner
            .comm
            .alltoall_f64(&send, block, &mut recv)
            .map_err(err)?;
        let mut x = RealMatrix::zeros(fine_rows.len(), b).map_err(SolverError::from)?;
        for src in 0..p {
            for (i, &r) in edges[src * p + rank].iter().enumerate() {
                for j in 0..b {
                    x.as_mut_slice()[j * fine_rows.len() + slot[mapping[r]]] =
                        recv[src * block + i * b + j];
                }
            }
        }
        crate::cutoff_ladder::handoff_owned_count(0, 1);
        crate::cutoff_ladder::handoff_owned_count(1, 8 * send.len());
        crate::cutoff_ladder::handoff_owned_count(2, 8 * inner.x.as_slice().len());
        crate::cutoff_ladder::handoff_owned_peak(
            8 * (inner.x.as_slice().len() + x.as_slice().len() + send.len() + recv.len()),
        );
        Ok(Self::from_owned_handoff(
            x,
            fine_rows,
            half,
            inner.comm.clone(),
            populated
                .iter()
                .enumerate()
                .filter_map(|(i, &yes)| (!yes).then_some(i))
                .collect(),
        ))
    }
}

// Offline experiment only: no production switch or additional SCF work.
#[cfg(test)]
use super::pool_fixture as handoff_fixtures;

/// First-step eigen tolerance of the full-cutoff standard session after a
/// ladder hand-off.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StandardSeed {
    /// The ratchet cap, exactly as a cold standard session (the first-iteration
    /// redo rule then tightens it if the start is already good).
    Cap,
    /// The tolerance the coarse standard session would have used next:
    /// ratchet(last coarse tolerance, last coarse Coulomb error, N_e). The fine
    /// session's own ratchet takes over from its first record.
    Coarse,
}

/// Coarse stop level rule (D1).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CoarseStopRule {
    /// Stop at `coarse_density_residual` (the adopted rule).
    Fixed,
    /// Stop at max(`coarse_density_residual`, d) with d the distance (same norm
    /// as the density residual) between the atomic start density of the full
    /// cutoff and the prolonged atomic start density of the coarse cutoff: the
    /// part of a density the coarse grid cannot hold, so a coarse residual
    /// below d is below the coarse-to-fine model difference. No constant; it can
    /// only stop the coarse stage earlier. Without an atomic start (uniform
    /// start density) d is unavailable and the fixed level applies (recorded).
    Atomic,
}

/// First-step eigen tolerance of the full-cutoff session after a hand-off (D1).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FirstTolerance {
    /// Adopted behaviour: standard follows `standard_seed`, strict starts its
    /// residual schedule at its cold value (1e-2).
    Auto,
    /// Both profiles start cold (standard cap, strict 1e-2).
    Cap,
    /// Both start from the eigen tolerance the last coarse step requested (the
    /// coarse orbitals are converged to it).
    Coarse,
    /// Both start from the tolerance the coarse schedule would use next
    /// (standard ratchet; strict factor x last coarse density residual).
    Next,
}

/// Ladder controls. Off unless a caller runs `run_ladder`.
#[derive(Clone, Debug)]
pub struct CutoffLadderOptions {
    /// Coarse cutoff as a fraction of the full cutoff, in (0, 1).
    pub fraction: f64,
    /// Hand off when the coarse density residual (Bohr^-3/2) is at or below this.
    pub coarse_density_residual: f64,
    /// Budget of coarse iterations; the state reached is handed off either way.
    pub coarse_max_iterations: usize,
    /// S1: the standard profile is eligible (`BAANDI_CUTOFF_LADDER_STANDARD=on`;
    /// default off).
    pub standard: bool,
    /// S1: first-step tolerance rule of the fine standard session
    /// (`BAANDI_LADDER_STANDARD_SEED=cap|coarse`; default coarse).
    pub standard_seed: StandardSeed,
    /// S2: number h of coarse Pulay difference pairs prolonged into the fine
    /// mixer; 0 is off (`BAANDI_LADDER_SECANT=off|on|<h>`, on means 4; default off).
    pub secant_pairs: usize,
    /// D1: coarse stop level rule (`BAANDI_LADDER_COARSE_RULE=fixed|atomic`;
    /// default fixed).
    pub coarse_rule: CoarseStopRule,
    /// D1: first fine eigen tolerance (`BAANDI_LADDER_FIRST_TOL=auto|cap|coarse|next`;
    /// default auto).
    pub first_tolerance: FirstTolerance,
}
impl Default for CutoffLadderOptions {
    /// Pre-registered values (fraction 0.5, residual 1e-4, budget 30) and the
    /// S1/S2 experiment switches read from the environment (an invalid value
    /// means off / the default).
    fn default() -> Self {
        Self::from_run(dft_numerics::switches::RunConfig::current())
    }
}
impl CutoffLadderOptions {
    /// RUN-CONFIG: the pre-registered values with the experiment switches of
    /// `run` (the explicit per-call form of `default`).
    pub fn from_run(run: dft_numerics::switches::RunConfig) -> Self {
        Self {
            standard: run.on(dft_numerics::switches::Id::CUTOFF_LADDER_STANDARD),
            standard_seed: match run.choice(dft_numerics::switches::Id::LADDER_STANDARD_SEED) {
                "cap" => StandardSeed::Cap,
                _ => StandardSeed::Coarse,
            },
            secant_pairs: match run.choice(dft_numerics::switches::Id::LADDER_SECANT) {
                "off" => 0,
                "on" => 4,
                v => v.parse().unwrap_or(0),
            },
            coarse_density_residual: run
                .number(dft_numerics::switches::Id::LADDER_COARSE_STOP)
                .unwrap_or(1e-4),
            coarse_rule: match run.choice(dft_numerics::switches::Id::LADDER_COARSE_RULE) {
                "atomic" => CoarseStopRule::Atomic,
                _ => CoarseStopRule::Fixed,
            },
            first_tolerance: match run.choice(dft_numerics::switches::Id::LADDER_FIRST_TOL) {
                "cap" => FirstTolerance::Cap,
                "coarse" => FirstTolerance::Coarse,
                "next" => FirstTolerance::Next,
                _ => FirstTolerance::Auto,
            },
            ..Self::baseline()
        }
    }
}
impl CutoffLadderOptions {
    /// The pre-registered ladder with both experiment switches off (no
    /// environment read): the strict-only ladder as adopted.
    pub fn baseline() -> Self {
        Self {
            fraction: 0.5,
            coarse_density_residual: 1e-4,
            coarse_max_iterations: 30,
            standard: false,
            standard_seed: StandardSeed::Coarse,
            secant_pairs: 0,
            coarse_rule: CoarseStopRule::Fixed,
            first_tolerance: FirstTolerance::Auto,
        }
    }
    pub fn validate(&self) -> Result<(), ScfError> {
        if !(self.fraction.is_finite() && self.fraction > 0.0 && self.fraction < 1.0)
            || !(self.coarse_density_residual.is_finite() && self.coarse_density_residual > 0.0)
            || self.coarse_max_iterations == 0
            || self.secant_pairs > 64
        {
            return Err(ScfError::InvalidInput(
                "cutoff_ladder: fraction in (0,1), positive residual and budget, secant <= 64 required",
            ));
        }
        Ok(())
    }
}

/// Work and outcome of the coarse stage. H columns are counted in the coarse
/// basis; `coarse_grid_points / fine_grid_points` converts them to full-grid
/// equivalents for an H-cost comparison.
#[derive(Clone, Debug, Default)]
pub struct CutoffLadderReport {
    /// True when the full-cutoff session started from the prolonged state.
    pub used: bool,
    pub fallback_reason: Option<String>,
    pub fraction: f64,
    pub coarse_cutoff_hartree: f64,
    pub coarse_grid: [usize; 3],
    pub fine_grid: [usize; 3],
    pub coarse_plane_waves: usize,
    pub fine_plane_waves: usize,
    pub coarse_iterations: usize,
    pub coarse_final_density_residual: f64,
    pub coarse_davidson_iterations: usize,
    /// Ordinary solve H columns of the coarse history (no certificate ran).
    pub coarse_h_vectors: usize,
    /// H columns of the coarse atomic start (counted in the coarse basis).
    pub coarse_start_h_vectors: usize,
    pub coarse_setup: Duration,
    /// Setup, coarse session and prolongation together (unchanged meaning).
    pub coarse_time: Duration,
    /// P48F diagnostics, parts of `coarse_time`: wall of the coarse session
    /// (initialization and steps) and of the hand-off (orbital and density
    /// prolongation with its checks; S2 secant pairs are not included).
    pub coarse_session: Duration,
    pub coarse_prolong: Duration,
    /// Phase timings summed over the coarse history (certificates included).
    pub coarse_phases: PhaseTimings,
    /// Process clocks accumulated inside the coarse session (atomic start,
    /// step wall, certificate blocks).
    pub coarse_clock: crate::phase_clock::PhaseClock,
    /// The full-cutoff session runs the standard profile (S1).
    pub standard: bool,
    /// First-step eigen tolerance given to the fine standard session (Ha);
    /// None means the cap (cold-start rule).
    pub standard_seed_ha: Option<f64>,
    /// S2: requested and seeded difference pairs, and why none were seeded.
    pub secant_requested: usize,
    pub secant_seeded: usize,
    pub secant_fallback: Option<String>,
    /// D1: rule and level (Bohr^-3/2) that stopped the coarse stage.
    pub coarse_rule: Option<CoarseStopRule>,
    pub coarse_stop_level: f64,
    /// D1: distance between the full-cutoff atomic start density and the
    /// prolonged coarse one (density-residual norm, Bohr^-3/2); None without an
    /// atomic start.
    pub atomic_truncation_distance: Option<f64>,
    /// D1: distance between the converged full-cutoff density and the handed-off
    /// (prolonged coarse) density, same norm: the measured coarse-to-fine
    /// density difference of this run.
    pub handoff_to_final_density_distance: Option<f64>,
    /// D1: first-step tolerance rule and the strict seed given (Ha; None means
    /// the cold schedule start).
    pub first_tolerance: Option<FirstTolerance>,
    pub strict_seed_ha: Option<f64>,
    /// D1: full-cutoff work after the hand-off: SCF steps, ordinary (recorded)
    /// H columns and certificate H columns.
    pub fine_iterations: usize,
    pub fine_h_vectors: usize,
    pub fine_certificate_h_vectors: usize,
    /// SCF-TRAJECTORY: the coarse session's per-map records (the fine
    /// history is the run's own history).
    pub coarse_history: Vec<IterationRecord>,
    /// Effective start of the coarse stage (the first stage of the run): the
    /// density it started from ("atomic" = superposition of the UPF atomic
    /// charges, "uniform") and its atomic start-wavefunction report (None when
    /// `start_wavefunctions` is not atomic). The full-cutoff stage always
    /// starts from the prolonged coarse density and orbitals.
    pub coarse_start_density: Option<&'static str>,
    pub coarse_atomic_start: Option<AtomicStartReport>,
}
impl CutoffLadderReport {
    /// D1: coarse stage cost in fine-step equivalents: full-grid-equivalent
    /// coarse H columns divided by the mean recorded H columns of one fine SCF
    /// step (None before a fine step ran).
    pub fn coarse_fine_step_equivalents(&self) -> Option<f64> {
        (self.fine_iterations > 0 && self.fine_h_vectors > 0).then(|| {
            self.coarse_h_equivalents() / (self.fine_h_vectors as f64 / self.fine_iterations as f64)
        })
    }
    /// Coarse H columns scaled by the FFT grid ratio (H cost per column is
    /// dominated by the grid transforms): full-grid H-column equivalents.
    pub fn coarse_h_equivalents(&self) -> f64 {
        let ratio = self.coarse_grid.iter().product::<usize>() as f64
            / self.fine_grid.iter().product::<usize>().max(1) as f64;
        (self.coarse_h_vectors + self.coarse_start_h_vectors) as f64 * ratio
    }
}

/// Why a problem cannot use the ladder (None means eligible).
pub fn ineligibility(problem: &ScfProblem, ladder: &CutoffLadderOptions) -> Option<String> {
    if ladder.validate().is_err() {
        return Some("invalid ladder options".into());
    }
    if problem.states.len() != 1 {
        return Some("multiple k points".into());
    }
    if let Some(reason) = problem.gamma_real_exclusion {
        return Some(reason.into());
    }
    if problem.states[0].basis.k_fractional() != [0.; 3] {
        return Some("k is not exactly Gamma".into());
    }
    let o = &problem.options;
    if o.gamma_real == crate::GammaReal::Off {
        return Some("gamma_real off".into());
    }
    if o.standard.is_some() && !ladder.standard {
        return Some("standard profile".into());
    }
    if o.locking != dft_solvers::locking::Locking::Off
        || o.spectral_guard != crate::SpectralGuard::Off
        || !matches!(o.eigensolver, dft_solvers::SolverChoice::Davidson)
    {
        return Some("solver options outside the real Gamma path".into());
    }
    if o.strain_reference.is_some() {
        return Some("frozen-Miller strain reference".into());
    }
    if o.compat_mode != CompatMode::None {
        return Some("compatibility mode".into());
    }
    if problem.symmetry.is_some() {
        return Some("density symmetrization enabled".into());
    }
    let occupied = problem.occupations.iter().filter(|&&f| f > 0.).count();
    if occupied == 0
        || problem
            .occupations
            .iter()
            .any(|&f| f != 0. && f != problem.occupations[0])
        || problem.occupations[..occupied].iter().any(|&f| f <= 0.)
    {
        return Some("occupations are not fixed and uniform".into());
    }
    None
}

/// Copy coefficients by Miller index into the fine basis (zero for new G).
/// Fails if a coarse G is missing from the fine basis (bases not nested).
pub fn prolong_orbitals(
    coarse: &PlaneWaveBasis,
    fine: &PlaneWaveBasis,
    orbitals: &ComplexMatrix,
) -> Result<ComplexMatrix, ScfError> {
    let _span = HandoffSpan::enter(HandoffPart::OrbitalCopy);
    if orbitals.rows() != coarse.vectors().len() {
        return Err(ScfError::InvalidInput(
            "cutoff_ladder: orbital rows differ from the coarse basis",
        ));
    }
    let index: HashMap<[i32; 3], usize> = fine
        .vectors()
        .iter()
        .enumerate()
        .map(|(i, g)| (g.miller, i))
        .collect();
    let map = coarse
        .vectors()
        .iter()
        .map(|g| index.get(&g.miller).copied())
        .collect::<Option<Vec<_>>>()
        .ok_or(ScfError::InvalidInput(
            "cutoff_ladder: coarse basis is not a subset of the fine basis",
        ))?;
    let rows = fine.vectors().len();
    let mut out = ComplexMatrix::zeros(rows, orbitals.cols())
        .map_err(|e| ScfError::Numerical(format!("cutoff_ladder: {e}")))?;
    for j in 0..orbitals.cols() {
        let dst = out.column_mut(j);
        for (src, &i) in orbitals.column(j).iter().zip(&map) {
            dst[i] = *src;
        }
    }
    Ok(out)
}

/// Fourier interpolation of a real periodic field from the coarse FFT grid to a
/// finer one (zero padding). The Nyquist plane of an even coarse dimension is
/// dropped: a field band-limited to the coarse density sphere has no weight
/// there because the grid holds that sphere (dims >= 2 max|m| + 1).
pub fn prolong_density(
    coarse: &Array3<f64>,
    fine_dims: [usize; 3],
) -> Result<Array3<f64>, ScfError> {
    let cd = [coarse.dim().0, coarse.dim().1, coarse.dim().2];
    let mut fft = FftPlan3d::new(cd)?;
    let input = coarse.mapv(|v| Complex64::new(v, 0.0));
    let mut reciprocal = Array3::<Complex64>::zeros(cd);
    fft.to_reciprocal(&input, &mut reciprocal)?;
    let mut padded = Array3::<Complex64>::zeros(fine_dims);
    let signed = |i: usize, n: usize| -> Option<i32> {
        if n.is_multiple_of(2) && i == n / 2 {
            None
        } else if i <= (n - 1) / 2 {
            Some(i as i32)
        } else {
            Some(i as i32 - n as i32)
        }
    };
    for ((i, j, k), z) in reciprocal.indexed_iter() {
        let (Some(a), Some(b), Some(c)) = (signed(i, cd[0]), signed(j, cd[1]), signed(k, cd[2]))
        else {
            continue;
        };
        let [p, q, r] = FftGridSpec::map_miller([a, b, c], fine_dims).map_err(|_| {
            ScfError::InvalidInput("cutoff_ladder: coarse frequencies exceed the fine grid")
        })?;
        padded[[p, q, r]] = *z;
    }
    let mut real = Array3::<Complex64>::zeros(fine_dims);
    let mut fine_fft = FftPlan3d::new(fine_dims)?;
    fine_fft.to_real(&padded, &mut real)?;
    Ok(real.mapv(|z| z.re))
}

fn coarse_problem(
    fine: &ScfProblem,
    potentials: &[SpeciesPotential<'_>],
    ladder: &CutoffLadderOptions,
) -> Result<ScfProblem, ScfError> {
    let mut options = fine.options.clone();
    options.cutoff_hartree *= ladder.fraction;
    // A fixed grid belongs to the full cutoff; the coarse grid follows E1.
    options.grid_override = None;
    // ADR 0043: the coarse stage only supplies a start; verification belongs to the fine run.
    options.verification = crate::Verification::Off;
    // RUN-CONFIG: the coarse problem is constructed and run under the fine
    // problem's configuration.
    let mut coarse = {
        let _run = fine.run_config().enter();
        ScfProblem::new(&fine.cell, &fine.atoms, potentials, options)?
    };
    coarse.set_run_config(fine.run_config());
    if let Some(floor) = fine.density_floor {
        coarse.set_density_floor(floor)?;
    }
    Ok(coarse)
}

/// What the coarse stage hands to the full cutoff.
struct CoarseState {
    density: Array3<f64>,
    orbitals: ComplexMatrix,
    owned: Option<crate::gamma_real::LazyOrbitals>,
    basis: PlaneWaveBasis,
    /// Pulay memory of the coarse session (only when S2 asks for it).
    mixer: Option<crate::checkpoint::MixerState>,
    /// Standard profile: the tolerance the coarse ratchet would use next.
    standard_next: Option<f64>,
    /// Strict residual schedule: the tolerance the coarse session would use next.
    strict_next: Option<f64>,
    /// Eigen tolerance requested by the last coarse step.
    last_tolerance: Option<f64>,
}

/// Coarse stage: drive the ordinary session at E1 and return its next input
/// density and the warm orbitals of its last solve.
fn coarse_stage(
    coarse: ScfProblem,
    ladder: &CutoffLadderOptions,
    stop_level: f64,
    report: &mut CutoffLadderReport,
    owned_handoff: bool,
) -> Result<CoarseState, ScfError> {
    // P0.2: Davidson and operator work of this stage is reported as ladder coarse.
    let _run_stage = dft_numerics::work_timers::enter(dft_numerics::work_timers::RunStage::Coarse);
    // Phase profile and snapshot context: the coarse bucket.
    let _profile_stage = crate::profile::enter(crate::profile::Top::Coarse);
    let _context_stage = crate::profile::context_stage("coarse");
    let basis = coarse.states[0].basis.clone();
    let session_start = Instant::now();
    let clock = crate::phase_clock::snapshot();
    let setup_clock = crate::profile::clock();
    let mut session = ScfSession::new(
        coarse,
        CompletionMode::RetainLast,
        crate::events::ExecutionConfig::default(),
    )?;
    session.initialize()?;
    crate::profile::between(crate::profile::BetweenPart::SessionSetup, setup_clock);
    loop {
        let status = session.step()?.status;
        let last = session
            .history()
            .last()
            .map_or(f64::INFINITY, |r| r.density_residual);
        if status != StepStatus::Ready
            || last <= stop_level
            || session.history().len() >= ladder.coarse_max_iterations
        {
            break;
        }
    }
    report.coarse_session = session_start.elapsed();
    report.coarse_clock = crate::phase_clock::snapshot().since(&clock);
    let history = session.history();
    for t in history.iter().flat_map(|r| {
        std::iter::once(&r.timings).chain(r.fixed_point_certificate.as_ref().map(|c| &c.timings))
    }) {
        report.coarse_phases.potential += t.potential;
        report.coarse_phases.diagonalization += t.diagonalization;
        report.coarse_phases.density += t.density;
        report.coarse_phases.energy += t.energy;
        report.coarse_phases.mixing += t.mixing;
    }
    report.coarse_iterations = history.len();
    report.coarse_final_density_residual = history.last().map_or(f64::NAN, |r| r.density_residual);
    report.coarse_davidson_iterations = history.iter().map(|r| r.davidson_iterations).sum();
    report.coarse_h_vectors = history
        .iter()
        .map(|r| {
            r.operator_applications
                + r.fixed_point_certificate
                    .as_ref()
                    .map_or(0, |c| c.operator_applications)
        })
        .sum();
    report.coarse_start_h_vectors = session.problem_ref().atomic_start.initial_h_applications;
    report.coarse_history = history.to_vec();
    let coarse_problem = session.problem_ref();
    report.coarse_start_density = Some(if coarse_problem.start_density.is_some() {
        "atomic"
    } else {
        "uniform"
    });
    report.coarse_atomic_start = (coarse_problem.options.start_wavefunctions
        == StartWavefunctions::Atomic)
        .then(|| coarse_problem.atomic_start.clone());
    let standard_next = history.last().and_then(|r| {
        r.standard_error_ha.map(|e| {
            StandardProfile::ratchet(
                r.requested_eigen_tolerance,
                e,
                session.problem_ref().electrons,
            )
        })
    });
    let strict_next = session.next_strict_tolerance();
    let last_tolerance = history.last().map(|r| r.requested_eigen_tolerance);
    let mixer = if ladder.secant_pairs > 0 {
        session.mixer_state()
    } else {
        None
    };
    let export_span = HandoffSpan::enter(HandoffPart::ExportCopy);
    let owned_state = if owned_handoff {
        session.warm_owned_handoff()?
    } else {
        None
    };
    let owned = owned_state.as_ref().map(|(_, h)| h.clone());
    // MPI-LAZY-ORBITALS: `warm_state` gathers a lazy coarse block (once).
    let (density, orbitals) = if let Some((density, h)) = owned_state {
        (
            density,
            ComplexMatrix::zeros(0, h.bands()).map_err(|e| ScfError::Numerical(e.to_string()))?,
        )
    } else {
        match session.warm_state()? {
            Some((n, Some(mut orbitals))) if orbitals.len() == 1 => {
                (n.clone(), orbitals.pop().expect("one block"))
            }
            _ => {
                // Terminal coarse session (converged or exhausted): use its result.
                let outcome = session.result().ok_or(ScfError::InvalidInput(
                    "cutoff_ladder: coarse state unavailable",
                ))?;
                // LAZY-RESULT: a lazily delivered coarse block is gathered here
                // (collective; every rank takes this branch together).
                let k = &outcome.result.k_points[0];
                let orbitals = match &k.lazy {
                    Some(lazy) => {
                        let _span = HandoffSpan::enter(HandoffPart::Materialize);
                        lazy.materialize()?.clone()
                    }
                    None => k.solution.eigenvectors.clone(),
                };
                (outcome.result.density.clone(), orbitals)
            }
        }
    };
    drop(export_span);
    Ok(CoarseState {
        density,
        orbitals,
        owned,
        basis,
        mixer,
        standard_next,
        strict_next,
        last_tolerance,
    })
}

/// Prolonged difference pairs (dx_i, dr_i), oldest first.
type SecantPairs = (Vec<Vec<f64>>, Vec<Vec<f64>>);

/// S2: the last `h` difference pairs (dx, dr) of the coarse Pulay memory,
/// oldest first, each prolonged to the fine grid by the exact Fourier
/// embedding. Differences are formed inside the coarse history only.
fn prolonged_secant_pairs(
    mixer: Option<&crate::checkpoint::MixerState>,
    coarse_dims: [usize; 3],
    fine_dims: [usize; 3],
    h: usize,
) -> Result<SecantPairs, String> {
    let Some(crate::checkpoint::MixerState::Pulay(state)) = mixer else {
        return Err("coarse mixer is not plain Pulay".into());
    };
    let size = coarse_dims.iter().product::<usize>();
    let pairs = state.inputs.len().min(state.residuals.len());
    if pairs < 2 {
        return Err("coarse history has no difference pair".into());
    }
    let first = (pairs - 1).saturating_sub(h);
    let prolong = |a: &[f64], b: &[f64]| -> Result<Vec<f64>, String> {
        if a.len() != size || b.len() != size {
            return Err("coarse history size differs from the coarse grid".into());
        }
        let d = Array3::from_shape_vec(coarse_dims, a.iter().zip(b).map(|(x, y)| x - y).collect())
            .map_err(|e| e.to_string())?;
        Ok(prolong_density(&d, fine_dims)
            .map_err(|e| e.to_string())?
            .into_iter()
            .collect())
    };
    let mut dx = Vec::with_capacity(pairs - 1 - first);
    let mut dr = Vec::with_capacity(pairs - 1 - first);
    for i in first..pairs - 1 {
        dx.push(prolong(&state.inputs[i + 1], &state.inputs[i])?);
        dr.push(prolong(&state.residuals[i + 1], &state.residuals[i])?);
    }
    Ok((dx, dr))
}

/// Result of a ladder run: the full-cutoff outcome and the coarse-stage report.
pub struct LadderOutcome {
    pub outcome: ScfOutcome,
    pub report: CutoffLadderReport,
}

/// Run the ladder on an eligible full-cutoff problem, or the ordinary SCF with
/// the reason recorded. `potentials` must be the ones `fine` was built from.
/// Only the full-cutoff session's certificate can set convergence.
pub fn run_ladder(
    fine: ScfProblem,
    potentials: &[SpeciesPotential<'_>],
    ladder: &CutoffLadderOptions,
    retain_last: bool,
) -> Result<LadderOutcome, ScfError> {
    let _run = fine.run_config().enter();
    let mut report = CutoffLadderReport {
        standard: fine.options.standard.is_some(),
        secant_requested: ladder.secant_pairs,
        fraction: ladder.fraction,
        coarse_rule: Some(ladder.coarse_rule),
        coarse_stop_level: ladder.coarse_density_residual,
        first_tolerance: Some(ladder.first_tolerance),
        coarse_cutoff_hartree: fine.options.cutoff_hartree * ladder.fraction,
        fine_grid: fine.grid.dims(),
        fine_plane_waves: fine.states[0].basis.vectors().len(),
        ..Default::default()
    };
    if let Some(reason) = ineligibility(&fine, ladder) {
        report.fallback_reason = Some(reason);
        let outcome = direct_min::finish(fine, None, retain_last)?;
        return Ok(LadderOutcome { outcome, report });
    }
    let requested = dft_numerics::switches::on(dft_numerics::switches::Id::M2_OWNED_HANDOFF);
    let eligible = requested
        && fine.options.standard.is_some()
        && fine.owned_declined().is_none()
        && crate::snapshot::dir().is_none()
        && crate::gamma_real::owned_handoff_supported();
    let owned_handoff = if requested {
        if let Some(comm) = crate::parallel::SharedComm::active() {
            let _span = HandoffSpan::enter(HandoffPart::Transfer);
            let mut count = [0.];
            comm.tree_sum_f64(&[f64::from(eligible)], &mut count)
                .map_err(|e| ScfError::Numerical(format!("owned handoff eligibility: {e}")))?;
            if count[0] != 0. && count[0] != comm.size() as f64 {
                return Err(ScfError::InvalidInput(
                    "owned handoff: eligibility differs across ranks",
                ));
            }
            count[0] == comm.size() as f64
        } else {
            false
        }
    } else {
        false
    };
    let setup = Instant::now();
    let start = match coarse_problem(&fine, potentials, ladder) {
        Ok(coarse) => {
            report.coarse_setup = setup.elapsed();
            // BREAKTHROUGH-4: the coarse problem's construction is `init`
            // work (the phase profile's construction bucket).
            crate::profile::add_init(report.coarse_setup.as_secs_f64());
            let handoff_clock = crate::profile::clock();
            report.coarse_grid = coarse.grid.dims();
            report.coarse_plane_waves = coarse.states[0].basis.vectors().len();
            report.atomic_truncation_distance = atomic_truncation_distance(&fine, &coarse);
            crate::profile::between(crate::profile::BetweenPart::LadderTruncation, handoff_clock);
            if ladder.coarse_rule == CoarseStopRule::Atomic
                && let Some(d) = report.atomic_truncation_distance
            {
                report.coarse_stop_level = report.coarse_stop_level.max(d);
            }
            let stop_level = report.coarse_stop_level;
            coarse_stage(coarse, ladder, stop_level, &mut report, owned_handoff).and_then(|c| {
                let prolong = Instant::now();
                let handoff_clock = crate::profile::clock();
                let state = if let Some(h) = &c.owned {
                    h.prolong_owned(
                        &c.basis,
                        &fine.states[0].basis,
                        fine.states[0].hamiltonian.fft_grid(),
                    )
                    .and_then(|owned| {
                        handoff_density(&fine, &c.density).map(|state| (state, Some(owned)))
                    })
                } else {
                    hand_off(&fine, &c.basis, &c.density, &c.orbitals).map(|state| (state, None))
                };
                crate::profile::between(crate::profile::BetweenPart::LadderHandoff, handoff_clock);
                report.coarse_prolong = prolong.elapsed();
                state.map(|state| (state, c))
            })
        }
        Err(e) => Err(e),
    };
    let (electronic, seeds, owned) = match start {
        Ok(((state, owned), coarse)) => {
            report.used = true;
            let standard = fine.options.standard.is_some();
            let seed = match ladder.first_tolerance {
                FirstTolerance::Auto => (standard && ladder.standard_seed == StandardSeed::Coarse)
                    .then_some(coarse.standard_next)
                    .flatten(),
                FirstTolerance::Cap => None,
                FirstTolerance::Coarse => coarse.last_tolerance,
                FirstTolerance::Next if standard => coarse.standard_next,
                FirstTolerance::Next => coarse.strict_next,
            };
            if standard {
                report.standard_seed_ha = seed;
            } else if fine.options.residual_eigen_schedule {
                report.strict_seed_ha = seed;
            }
            let seeds = if ladder.secant_pairs > 0 {
                match prolonged_secant_pairs(
                    coarse.mixer.as_ref(),
                    report.coarse_grid,
                    fine.grid.dims(),
                    ladder.secant_pairs,
                ) {
                    Ok(pairs) => Some(pairs),
                    Err(reason) => {
                        report.secant_fallback = Some(reason);
                        None
                    }
                }
            } else {
                None
            };
            (Some(state), seeds, owned)
        }
        Err(e) => {
            report.fallback_reason = Some(format!("coarse stage: {e}"));
            (None, None, None)
        }
    };
    report.coarse_time = setup.elapsed();
    let Some(electronic) = electronic else {
        let outcome = direct_min::finish(fine, None, retain_last)?;
        return Ok(LadderOutcome { outcome, report });
    };
    let handoff = Array3::from_shape_vec(electronic.grid, electronic.density.clone())
        .map_err(|e| ScfError::Numerical(format!("cutoff_ladder: {e}")))?;
    let (volume, electrons) = (fine.cell.volume_bohr3(), fine.electrons);
    let outcome = finish_from(fine, &electronic, seeds, retain_last, &mut report, owned)?;
    let h = &outcome.result.history;
    report.fine_iterations = h.len();
    report.fine_h_vectors = h.iter().map(|r| r.operator_applications).sum();
    report.fine_certificate_h_vectors = h
        .iter()
        .filter_map(|r| r.fixed_point_certificate.as_ref())
        .map(|c| c.operator_applications)
        .sum();
    report.handoff_to_final_density_distance = Some(density_distance(
        &outcome.result.density,
        &handoff,
        volume,
        electrons,
    ))
    .filter(|d| d.is_finite());
    Ok(LadderOutcome { outcome, report })
}

/// D1: distance (density-residual norm) between the full-cutoff atomic start
/// density and the coarse one prolonged to the full grid; None without an
/// atomic start on both levels.
fn atomic_truncation_distance(fine: &ScfProblem, coarse: &ScfProblem) -> Option<f64> {
    let (f, c) = (
        fine.atomic_superposition_density()?,
        coarse.atomic_superposition_density()?,
    );
    let prolonged = prolong_density(c, fine.grid.dims()).ok()?;
    Some(density_distance(
        f,
        &prolonged,
        fine.cell.volume_bohr3(),
        fine.electrons,
    ))
    .filter(|d| d.is_finite())
}

/// `direct_min::finish` from the prolonged state, plus the optional S1
/// standard first-step tolerance and S2 mixer seeds (each recorded; a rejected
/// seed is a recorded fallback, never an error). Only the session's own
/// convergence rule ends the run.
fn finish_from(
    fine: ScfProblem,
    start: &ElectronicState,
    seeds: Option<SecantPairs>,
    retain_last: bool,
    report: &mut CutoffLadderReport,
    owned: Option<crate::gamma_real::LazyOrbitals>,
) -> Result<ScfOutcome, ScfError> {
    let completion = if retain_last {
        CompletionMode::RetainLast
    } else {
        CompletionMode::Strict
    };
    let setup_clock = crate::profile::clock();
    let handoff_setup = HandoffSpan::enter(HandoffPart::Setup);
    let mut session = ScfSession::new(fine, completion, crate::events::ExecutionConfig::default())?;
    if let Some(owned) = owned {
        session.initialize_owned_ladder(start, owned)?;
    } else {
        session.initialize_ladder_electronic(start)?;
    }
    if let Some(t) = report.standard_seed_ha {
        session.seed_standard_tolerance(t)?;
    }
    if let Some(t) = report.strict_seed_ha {
        session.seed_strict_tolerance(t)?;
    }
    if let Some((dx, dr)) = seeds {
        match session.seed_mixer_differences(&dx, &dr) {
            Ok(()) => report.secant_seeded = dx.len(),
            Err(e) => report.secant_fallback = Some(e.to_string()),
        }
    }
    crate::profile::between(crate::profile::BetweenPart::SessionSetup, setup_clock);
    drop(handoff_setup);
    while session.step()?.status == StepStatus::Ready {}
    if !retain_last && session.result().is_none() {
        return Err(ScfError::NotConverged {
            history: session.history().to_vec(),
        });
    }
    session
        .into_result()
        .ok_or(ScfError::InvalidInput("result unavailable"))
}

/// Prolonged density and orbitals as the start state of the full-cutoff session,
/// with the session's own import checks repeated here so that a rejection is a
/// recorded fallback, not an error.
fn hand_off(
    fine: &ScfProblem,
    coarse_basis: &PlaneWaveBasis,
    density: &Array3<f64>,
    orbitals: &ComplexMatrix,
) -> Result<ElectronicState, ScfError> {
    let fine_basis = &fine.states[0].basis;
    handoff_storage(
        coarse_basis.vectors().len(),
        fine_basis.vectors().len(),
        orbitals.cols(),
    );
    let prolonged = prolong_orbitals(coarse_basis, fine_basis, orbitals)?;
    let state = handoff_density(fine, density)?;
    let _encode = HandoffSpan::enter(HandoffPart::Encode);
    let mut coefficients = Vec::with_capacity(2 * prolonged.rows() * prolonged.cols());
    coefficients.extend(prolonged.as_slice().iter().flat_map(|z| [z.re, z.im]));
    Ok(ElectronicState {
        grid: fine.grid.dims(),
        density: state.density,
        orbitals: vec![Orbitals {
            rows: prolonged.rows(),
            cols: prolonged.cols(),
            coefficients,
        }],
        occupations: fine.occupations.clone(),
    })
}

fn handoff_density(fine: &ScfProblem, density: &Array3<f64>) -> Result<ElectronicState, ScfError> {
    let density_span = HandoffSpan::enter(HandoffPart::DensityCopy);
    let n = prolong_density(density, fine.grid.dims())?;
    if n.iter()
        .any(|v| !v.is_finite() || *v < -xc::NEGATIVE_ROUNDOFF)
    {
        return Err(ScfError::InvalidInput(
            "cutoff_ladder: prolonged density is negative or nonfinite",
        ));
    }
    let charge = n.sum() * fine.cell.volume_bohr3() / n.len() as f64;
    if (charge - fine.electrons).abs() / fine.electrons > 1e-8 {
        return Err(ScfError::InvalidInput(
            "cutoff_ladder: prolonged density charge error exceeds 1e-8",
        ));
    }
    drop(density_span);
    Ok(ElectronicState {
        grid: fine.grid.dims(),
        density: n.iter().copied().collect(),
        orbitals: Vec::new(),
        occupations: fine.occupations.clone(),
    })
}

#[cfg(test)]
mod handoff_experiment {
    use super::*;
    use dft_numerics::{real::RealMatrix, switches::RunConfig, work_timers};
    use dft_physics::real_gamma::RealHamiltonian;
    use std::{collections::BTreeMap, path::Path};

    struct Fixture {
        cell: Cell,
        atoms: Vec<Atom>,
        species: Vec<(String, Pseudopotential)>,
        options: ScfOptions,
    }
    impl Fixture {
        fn load(path: &Path) -> Self {
            let text = std::fs::read_to_string(path).unwrap();
            let mut rows = Vec::new();
            let mut atoms = Vec::new();
            let mut species = Vec::new();
            let mut options = ScfOptions::default();
            for line in text.lines() {
                let p: Vec<_> = line.split_whitespace().collect();
                let number = |i: usize| p[i].parse::<f64>().unwrap();
                match p[0] {
                    "row" => rows.push([number(1), number(2), number(3)]),
                    "atom" => atoms.push(Atom::new(p[1], [number(2), number(3), number(4)])),
                    "cutoff" => options.cutoff_hartree = number(1),
                    "seed" => options.seed = p[1].parse().unwrap(),
                    "bands" => options.bands = Some(p[1].parse().unwrap()),
                    "pbe" => options.functional = Functional::Pbe,
                    "pp" => species.push((
                        p[1].into(),
                        dft_upf::parse_upf(&std::fs::read_to_string(p[2]).unwrap()).unwrap(),
                    )),
                    "synthetic" => species.push((
                        p[1].into(),
                        handoff_fixtures::potential(number(2), number(3), p[4] == "nonlocal", 401),
                    )),
                    _ => panic!("unknown fixture line {line}"),
                }
            }
            Self {
                cell: Cell::from_user_rows_bohr(rows.try_into().unwrap()).unwrap(),
                atoms,
                species,
                options,
            }
        }
        fn refs(&self) -> Vec<SpeciesPotential<'_>> {
            self.species
                .iter()
                .map(|(label, potential)| SpeciesPotential { label, potential })
                .collect()
        }
        fn problem(&self, options: ScfOptions) -> ScfProblem {
            let mut p = ScfProblem::new(&self.cell, &self.atoms, &self.refs(), options).unwrap();
            if p.core_density().is_some() {
                p.set_density_floor(1e-10).unwrap();
            }
            p
        }
    }

    // The high-G projection is onto the complement of the coarse PW space,
    // not merely the complement of the occupied subspace. Full-complex rows
    // include both +/-G, so no extra factor of two is applied.
    fn weighted_residual(
        basis: &PlaneWaveBasis,
        coarse: &PlaneWaveBasis,
        x: &ComplexMatrix,
        hx: &ComplexMatrix,
        occupations: &[f64],
        shift: f64,
    ) -> (f64, f64) {
        let coarse: std::collections::HashSet<_> =
            coarse.vectors().iter().map(|g| g.miller).collect();
        let mut high = 0.;
        let mut low = 0.;
        for (j, &occupation) in occupations.iter().take(x.cols()).enumerate() {
            let projected: Vec<Complex64> = (0..x.cols())
                .map(|k| {
                    x.column(k)
                        .iter()
                        .zip(hx.column(j))
                        .map(|(a, b)| a.conj() * b)
                        .sum()
                })
                .collect();
            for (i, g) in basis.vectors().iter().enumerate() {
                let r = hx.column(j)[i]
                    - projected
                        .iter()
                        .enumerate()
                        .map(|(k, &a)| x.column(k)[i] * a)
                        .sum::<Complex64>();
                let v = occupation * r.norm_sqr() / (0.5 * g.norm_squared + shift);
                if coarse.contains(&g.miller) {
                    low += v;
                } else {
                    high += v;
                }
            }
        }
        (high, low)
    }

    #[test]
    fn high_g_indicator_has_the_known_diagonal_value() {
        let cell = handoff_fixtures::cell();
        let coarse = PlaneWaveBasis::new(&cell, dft_units::Hartree(0.01)).unwrap();
        let fine = PlaneWaveBasis::new(&cell, dft_units::Hartree(2.)).unwrap();
        let mut x = ComplexMatrix::zeros(fine.vectors().len(), 1).unwrap();
        let mut hx = x.clone();
        let zero = fine
            .vectors()
            .iter()
            .position(|g| g.miller == [0; 3])
            .unwrap();
        let high = fine
            .vectors()
            .iter()
            .position(|g| g.miller != [0; 3])
            .unwrap();
        x.column_mut(0)[zero] = Complex64::new(1., 0.);
        hx.column_mut(0)[zero] = Complex64::new(3., 0.);
        hx.column_mut(0)[high] = Complex64::new(0.2, 0.3);
        let (a, b) = weighted_residual(&fine, &coarse, &x, &hx, &[2.], 1.);
        let expected = 2. * 0.13 / (0.5 * fine.vectors()[high].norm_squared + 1.);
        assert!((a - expected).abs() < 1e-15);
        assert_eq!(b, 0.);
        assert_eq!(weighted_residual(&fine, &fine, &x, &hx, &[2.], 1.).0, 0.);
        assert_eq!(
            weighted_residual(&fine, &coarse, &x, &hx, &[0.], 1.),
            (0., 0.)
        );
    }

    // Offline coarse-only schedule experiment. Its factor=1 control must
    // reproduce coarse_stage bit for bit on every fixture before factor=4 is
    // measured. No session policy or production tolerance is modified.
    fn loose_stage(
        mut problem: ScfProblem,
        ladder: &CutoffLadderOptions,
        factor: f64,
        report: &mut CutoffLadderReport,
    ) -> (CoarseState, Vec<f64>, Vec<f64>) {
        use dft_solvers::mixing::DensityMixer;
        let _scope = problem.run_config().enter();
        let _stage = work_timers::enter(work_timers::RunStage::Coarse);
        let _profile = crate::profile::enter(crate::profile::Top::Coarse);
        let _context = crate::profile::context_stage("coarse");
        let clock = Instant::now();
        let Mixing::Pulay { beta, history } = problem.options.mixing else {
            panic!("Pulay required")
        };
        assert!(!problem.options.hartree_mixing);
        assert_eq!(
            problem
                .run_config()
                .choice(dft_numerics::switches::Id::MIX_COULOMB_GUARD),
            "off"
        );
        let mut mixer = PulayMixer::new(beta, history).unwrap();
        let mut n = problem.starting_density();
        let mut warm: Option<Vec<ComplexMatrix>> = None;
        let mut tolerance = StandardProfile::CAP_HA;
        let mut errors = Vec::new();
        let mut tolerances = Vec::new();
        for iteration in 1..=ladder.coarse_max_iterations {
            crate::profile::set_context(crate::profile::Context {
                stage: "coarse",
                kind: "ordinary",
                iteration,
            });
            if let Some(&error) = errors.last() {
                tolerance = (factor
                    * StandardProfile::ratchet(tolerance, error, problem.electrons))
                .min(tolerance)
                .min(StandardProfile::CAP_HA);
            }
            let mut step = problem
                .step_session(
                    &n,
                    warm.as_deref(),
                    &[],
                    tolerance,
                    tolerance <= problem.options.eigen_tolerance,
                )
                .unwrap();
            let mut error = problem.stop_error(&n, &step.density).unwrap().unwrap();
            if iteration == 1 && error < tolerance * problem.electrons.max(1.) {
                tolerance = (factor
                    * StandardProfile::ratchet(tolerance, error, problem.electrons))
                .min(tolerance);
                let guess = vec![step.k_points[0].solution.eigenvectors.clone()];
                crate::profile::set_context(crate::profile::Context {
                    stage: "coarse",
                    kind: "redo",
                    iteration,
                });
                let mut redo = problem
                    .step_session(&n, Some(&guess), &[], tolerance, false)
                    .unwrap();
                for (new, old) in redo.k_points.iter_mut().zip(&step.k_points) {
                    new.solution.iterations += old.solution.iterations;
                    new.solution.operator_applications += old.solution.operator_applications;
                }
                redo.timings.potential += step.timings.potential;
                redo.timings.diagonalization += step.timings.diagonalization;
                redo.timings.density += step.timings.density;
                redo.timings.energy += step.timings.energy;
                step = redo;
                error = problem.stop_error(&n, &step.density).unwrap().unwrap();
            }
            report.coarse_final_density_residual = density_distance(
                &step.density,
                &n,
                problem.cell.volume_bohr3(),
                problem.electrons,
            );
            report.coarse_h_vectors += step
                .k_points
                .iter()
                .map(|k| k.solution.operator_applications)
                .sum::<usize>();
            report.coarse_phases.potential += step.timings.potential;
            report.coarse_phases.diagonalization += step.timings.diagonalization;
            report.coarse_phases.density += step.timings.density;
            report.coarse_phases.energy += step.timings.energy;
            errors.push(error);
            tolerances.push(tolerance);
            let mix_clock = Instant::now();
            let mut mixed = mixer
                .mix(
                    n.as_slice().unwrap(),
                    step.density.as_slice().unwrap(),
                    problem.electrons / problem.cell.volume_bohr3(),
                )
                .unwrap();
            let fraction = n
                .iter()
                .zip(&mixed)
                .filter(|(_, b)| **b < 0.)
                .map(|(a, b)| 0.9 * a / (a - b))
                .fold(1., f64::min);
            if fraction < 1. {
                for (b, a) in mixed.iter_mut().zip(n.iter()) {
                    *b = a + fraction * (*b - a);
                }
                mixer.reset();
            }
            n.as_slice_mut().unwrap().copy_from_slice(&mixed);
            report.coarse_phases.mixing += mix_clock.elapsed();
            warm = Some(
                step.k_points
                    .into_iter()
                    .map(|k| k.solution.eigenvectors)
                    .collect(),
            );
            if report.coarse_final_density_residual <= ladder.coarse_density_residual {
                break;
            }
        }
        report.coarse_iterations = tolerances.len();
        report.coarse_start_h_vectors = problem.atomic_start.initial_h_applications;
        report.coarse_session = clock.elapsed();
        let state = CoarseState {
            basis: problem.states[0].basis.clone(),
            density: n,
            orbitals: warm.unwrap().remove(0),
            owned: None,
            mixer: None,
            standard_next: Some(StandardProfile::ratchet(
                tolerance,
                *errors.last().unwrap(),
                problem.electrons,
            )),
            strict_next: None,
            last_tolerance: Some(tolerance),
        };
        (state, tolerances, errors)
    }

    fn save_f64(path: &Path, values: impl Iterator<Item = f64>) {
        use std::io::Write;
        let mut file = std::io::BufWriter::new(std::fs::File::create(path).unwrap());
        for v in values {
            file.write_all(&v.to_le_bytes()).unwrap();
        }
        file.flush().unwrap();
    }

    /// Run through scripts/ladder_handoff.py. Every job has one geometry and
    /// one process, with a strict reference, candidates, and optional snapshots.
    #[test]
    #[ignore = "offline full-SCF experiment; release only, scripts/ladder_handoff.py"]
    fn measure() {
        let input = std::env::var("HANDOFF_INPUT").expect("HANDOFF_INPUT");
        let fixture = Fixture::load(Path::new(&input));
        let root = std::env::var("HANDOFF_OUTPUT").expect("HANDOFF_OUTPUT");
        let root = Path::new(&root);
        std::fs::create_dir_all(root).unwrap();
        let run = RunConfig::from_pairs(&[
            ("phase_profile", "on"),
            ("stage_timers", "on"),
            ("davidson_stages", "on"),
        ])
        .unwrap();
        let _scope = run.enter();
        work_timers::set_enabled(true);
        crate::snapshot::force_dir(Some(None));
        let mut strict = fixture.options.clone();
        strict.standard = None;
        strict.verification = Verification::Off;
        strict.mixing = ScfOptions::strict_reference().mixing;
        strict.energy_tolerance = 1e-11;
        strict.density_tolerance = 1e-9;
        strict.eigen_tolerance = 1e-10;
        let reference = fixture.problem(strict).run_outcome().unwrap();
        assert!(reference.converged, "strict reference did not converge");
        let e_ref = reference.result.total_hartree();
        let arms = std::env::var("HANDOFF_ARMS").unwrap_or_else(|_| "0.5:0.0001:default".into());
        for (index, arm) in arms.split(',').enumerate() {
            let parts: Vec<_> = arm.split(':').collect();
            let fraction: f64 = parts[0].parse().unwrap();
            let stop: f64 = parts[1].parse().unwrap();
            let mode = parts[2];
            let ladder = CutoffLadderOptions {
                fraction,
                coarse_density_residual: stop,
                standard: true,
                ..CutoffLadderOptions::from_run(run)
            };
            let dir = root.join(format!("arm-{index}"));
            std::fs::create_dir_all(&dir).unwrap();
            let capture = std::env::var("HANDOFF_CAPTURE").as_deref() == Ok("on");
            if capture {
                std::fs::create_dir_all(dir.join("snapshots")).unwrap();
                crate::snapshot::force_dir(Some(Some(dir.join("snapshots"))));
            }
            if mode == "loose4" {
                crate::snapshot::force_dir(Some(None));
                let f = fixture.problem(fixture.options.clone());
                let mut normal = CutoffLadderReport::default();
                let expected = coarse_stage(
                    coarse_problem(&f, &fixture.refs(), &ladder).unwrap(),
                    &ladder,
                    stop,
                    &mut normal,
                    false,
                )
                .unwrap();
                let mut control = CutoffLadderReport::default();
                let (actual, _, _) = loose_stage(
                    coarse_problem(&f, &fixture.refs(), &ladder).unwrap(),
                    &ladder,
                    1.,
                    &mut control,
                );
                assert_eq!(
                    actual.density, expected.density,
                    "coarse factor=1 density control"
                );
                assert_eq!(
                    actual.orbitals.as_slice(),
                    expected.orbitals.as_slice(),
                    "coarse factor=1 orbital control"
                );
                assert_eq!(control.coarse_h_vectors, normal.coarse_h_vectors);
                if capture {
                    crate::snapshot::force_dir(Some(Some(dir.join("snapshots"))));
                }
            }
            let begin = work_timers::snapshot();
            let reuse_begin = crate::gamma_real::reuse_psi_counters();
            let wall = Instant::now();
            let fine = fixture.problem(fixture.options.clone());
            assert!(ineligibility(&fine, &ladder).is_none());
            let fine_grid = fine.grid.dims();
            let mut coarse = coarse_problem(&fine, &fixture.refs(), &ladder).unwrap();
            if mode == "minimal-grid" {
                let mut options = coarse.options.clone();
                let density_basis = coarse.states[0].basis.density_basis(&fixture.cell).unwrap();
                let mut dims = [1; 3];
                for g in density_basis.vectors() {
                    for (i, d) in dims.iter_mut().enumerate() {
                        *d = (*d).max(2 * g.miller[i].unsigned_abs() as usize + 1);
                    }
                }
                options.grid_override = Some(dims);
                coarse = fixture.problem(options);
            } else {
                assert!(mode == "default" || mode == "loose4");
            }
            let coarse_grid = coarse.grid.dims();
            let mut report = CutoffLadderReport {
                standard: true,
                coarse_grid,
                fine_grid,
                ..Default::default()
            };
            let (c, tolerances, errors) = if mode == "loose4" {
                loose_stage(coarse, &ladder, 4., &mut report)
            } else {
                let c = coarse_stage(coarse, &ladder, stop, &mut report, false).unwrap();
                let t = report
                    .coarse_history
                    .iter()
                    .map(|r| r.requested_eigen_tolerance)
                    .collect::<Vec<_>>();
                let e = report
                    .coarse_history
                    .iter()
                    .map(|r| r.standard_error_ha.unwrap())
                    .collect::<Vec<_>>();
                (c, t, e)
            };
            let transfer = Instant::now();
            let state = hand_off(&fine, &c.basis, &c.density, &c.orbitals).unwrap();
            let transfer_s = transfer.elapsed().as_secs_f64();
            report.standard_seed_ha = c.standard_next;
            let outcome = finish_from(fine, &state, None, true, &mut report, None).unwrap();
            let wall_s = wall.elapsed().as_secs_f64();
            let end = work_timers::snapshot();
            let reuse_end = crate::gamma_real::reuse_psi_counters();
            crate::snapshot::force_dir(Some(None));
            assert!(outcome.converged, "candidate {arm} did not converge");
            assert!(
                outcome
                    .result
                    .verification
                    .as_ref()
                    .is_some_and(|v| v.passed())
            );
            if capture {
                assert!(
                    std::fs::read_dir(dir.join("snapshots")).unwrap().any(|e| e
                        .unwrap()
                        .path()
                        .extension()
                        .is_some_and(|s| s == "json")),
                    "snapshot capture produced no manifests"
                );
            }
            save_f64(
                &dir.join("coarse-orbitals.bin"),
                c.orbitals.as_slice().iter().flat_map(|z| [z.re, z.im]),
            );
            save_f64(&dir.join("coarse-density.bin"), c.density.iter().copied());
            // Explicit solved handoff state; the existing snapshots are solver entries.
            save_f64(
                &dir.join("handoff-density.bin"),
                state.density.iter().copied(),
            );
            save_f64(
                &dir.join("handoff-orbitals.bin"),
                state.orbitals[0].coefficients.iter().copied(),
            );
            // A separate fine problem prevents diagnostic H/cache changes from
            // entering the measured SCF or its work counters.
            let probe_clock = Instant::now();
            let mut probe = fixture.problem(fixture.options.clone());
            let probe_build_s = probe_clock.elapsed().as_secs_f64();
            let fields_clock = Instant::now();
            let n = Array3::from_shape_vec(state.grid, state.density.clone()).unwrap();
            let (_, _, potential) = probe.input_fields(&n).unwrap();
            probe.states[0]
                .hamiltonian
                .set_local_potential(potential)
                .unwrap();
            let probe_fields_s = fields_clock.elapsed().as_secs_f64();
            let basis = &probe.states[0].basis;
            let x_all = prolong_orbitals(&c.basis, basis, &c.orbitals).unwrap();
            let occupied = probe.occupations.iter().filter(|&&f| f > 0.).count();
            let mut x = ComplexMatrix::zeros(x_all.rows(), occupied).unwrap();
            for j in 0..occupied {
                x.column_mut(j).copy_from_slice(x_all.column(j));
            }
            let mut real = RealHamiltonian::from_complex(
                &probe.states[0].hamiltonian,
                basis,
                occupied,
                probe.pool.clone(),
            )
            .unwrap();
            let mut xr = RealMatrix::zeros(x.rows(), occupied).unwrap();
            let mut hr = RealMatrix::zeros(x.rows(), occupied).unwrap();
            for j in 0..occupied {
                real.half_sphere()
                    .pack(
                        x.column(j),
                        &mut xr.as_mut_slice()[j * x.rows()..(j + 1) * x.rows()],
                    )
                    .unwrap();
            }
            let apply_clock = Instant::now();
            real.apply_block(&xr, &mut hr).unwrap();
            let probe_apply_s = apply_clock.elapsed().as_secs_f64();
            let mut hx = ComplexMatrix::zeros(x.rows(), occupied).unwrap();
            for j in 0..occupied {
                real.half_sphere()
                    .unpack(
                        &hr.as_slice()[j * x.rows()..(j + 1) * x.rows()],
                        hx.column_mut(j),
                    )
                    .unwrap();
            }
            let indicators: Vec<_> = [0.5, 1., 2.]
                .iter()
                .map(|&shift| {
                    weighted_residual(basis, &c.basis, &x, &hx, &probe.occupations, shift)
                })
                .collect();
            let probe_s = probe_clock.elapsed().as_secs_f64();
            // Independent complex H application validates ORC packing and multiplicity.
            if index == 0 {
                use dft_solvers::LinearOperator;
                let mut independent = ComplexMatrix::zeros(x.rows(), occupied).unwrap();
                probe.states[0]
                    .hamiltonian
                    .apply_block(&x, &mut independent)
                    .unwrap();
                let error = independent
                    .as_slice()
                    .iter()
                    .zip(hx.as_slice())
                    .map(|(a, b)| (*a - *b).norm())
                    .fold(0., f64::max);
                assert!(error < 1e-11, "real/complex probe mismatch {error}");
            }
            let history = &outcome.result.history;
            let verify = outcome.result.verification.as_ref().unwrap();
            let mut row = BTreeMap::<String, String>::new();
            macro_rules! put {
                ($k:expr,$v:expr) => {
                    row.insert($k.into(), format!("{:?}", $v));
                };
            }
            put!("arm", arm);
            put!("capture", capture);
            put!("fraction", fraction);
            put!("stop", stop);
            put!("coarse_maps", report.coarse_iterations);
            put!("fine_maps", history.len());
            put!("coarse_grid", coarse_grid);
            put!("coarse_plane_waves", c.basis.vectors().len());
            put!("fine_plane_waves", x.rows());
            put!("bands", state.orbitals[0].cols);
            put!("occupied", occupied);
            put!("fine_grid", fine_grid);
            put!("coarse_h", report.coarse_h_vectors);
            put!("coarse_start_h", report.coarse_start_h_vectors);
            put!(
                "fine_history_h",
                history
                    .iter()
                    .map(|r| r.operator_applications)
                    .sum::<usize>()
            );
            put!(
                "verification_h",
                verify
                    .attempts
                    .iter()
                    .map(|v| v.operator_applications)
                    .sum::<usize>()
            );
            put!("verification_maps", verify.attempts.len());
            put!("coarse_residual", report.coarse_final_density_residual);
            put!("coarse_tolerances", tolerances);
            put!("coarse_standard_errors", errors);
            put!("loose_driver_identity", mode == "loose4");
            put!(
                "coarse_solve_s",
                report.coarse_phases.diagonalization.as_secs_f64()
            );
            put!(
                "coarse_density_s",
                report.coarse_phases.density.as_secs_f64()
            );
            put!(
                "coarse_potential_s",
                report.coarse_phases.potential.as_secs_f64()
            );
            put!("coarse_energy_s", report.coarse_phases.energy.as_secs_f64());
            put!("coarse_mixing_s", report.coarse_phases.mixing.as_secs_f64());
            put!("coarse_session_s", report.coarse_session.as_secs_f64());
            put!("transfer_s", transfer_s);
            put!("wall_s", wall_s);
            put!("energy_ha", outcome.result.total_hartree());
            put!("reference_ha", e_ref);
            put!("error_ha", (outcome.result.total_hartree() - e_ref).abs());
            put!(
                "indicator_high_ha",
                indicators.iter().map(|v| v.0).collect::<Vec<_>>()
            );
            put!(
                "indicator_low_ha",
                indicators.iter().map(|v| v.1).collect::<Vec<_>>()
            );
            put!("probe_h", occupied);
            put!("probe_calls", 1);
            put!("probe_s", probe_s);
            put!("probe_build_s", probe_build_s);
            put!("probe_fields_s", probe_fields_s);
            put!(
                "delta_apply_calls",
                reuse_end.armed
                    - reuse_begin.armed
                    - (reuse_end.unchanged_potential - reuse_begin.unchanged_potential)
            );
            put!("probe_apply_s", probe_apply_s);
            for (s, label) in ["fine", "coarse", "verification"].iter().enumerate() {
                put!(
                    format!("{label}_apply_calls"),
                    end[s].apply_calls - begin[s].apply_calls
                );
                put!(
                    format!("{label}_apply_h"),
                    end[s].apply_columns - begin[s].apply_columns
                );
                put!(
                    format!("{label}_wave_h"),
                    end[s].pair_columns + end[s].other_columns
                        - begin[s].pair_columns
                        - begin[s].other_columns
                );
                put!(
                    format!("{label}_kept_h"),
                    end[s].kept_columns - begin[s].kept_columns
                );
            }
            // One exact off-path oracle per job, including delivered density.
            if index == 0 && mode == "default" {
                let ordinary = run_ladder(
                    fixture.problem(fixture.options.clone()),
                    &fixture.refs(),
                    &ladder,
                    true,
                )
                .unwrap();
                assert_eq!(
                    ordinary.outcome.result.total_hartree().to_bits(),
                    outcome.result.total_hartree().to_bits()
                );
                assert_eq!(ordinary.outcome.result.density, outcome.result.density);
                assert_eq!(ordinary.report.coarse_h_vectors, report.coarse_h_vectors);
                assert_eq!(ordinary.report.fine_iterations, history.len());
                put!("baseline_identity", true);
            }
            let json = format!(
                "{{{}}}",
                row.iter()
                    .map(|(k, v)| format!("{k:?}:{v}"))
                    .collect::<Vec<_>>()
                    .join(",")
            );
            std::fs::write(dir.join("measurement.json"), format!("{json}\n")).unwrap();
            println!("HANDOFF {json}");
        }
    }
}
