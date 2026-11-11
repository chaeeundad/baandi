//! MPI stage M2d: distributed SCF stages on the M2 path (`BAANDI_MPI_M2=on` with a
//! communicator of two or more ranks; the default and the one-rank paths never
//! enter this module).
//!
//! * Density: occupation-weighted |psi_b|^2 of the distributed ORC eigenvectors on
//!   each rank's FFT planes (stick FFT, one batched transpose), summed over bands
//!   in ascending order on the planes and assembled into the replicated grid with
//!   one allgatherv (a pure permutation: the full density does not depend on how
//!   the planes are split). The per-band plane contributions are kept for the
//!   certificate screen (`weighted_integrals`: local sums plus one tree_sum).
//! * Potential: LDA (PZ, PW) exchange-correlation is pointwise; each rank evaluates
//!   a contiguous block of the row-major grid, the potential is assembled with one
//!   allgatherv (bit-identical values to the replicated evaluation) and the energy
//!   and potential integral are reduced with one deterministic tree_sum. Hartree
//!   (one 3D FFT pair) stays replicated.
//! * Energies: the nonlocal energy of the audit from the distributed projections
//!   (one tree_sum); the output-density XC energy as above.
//! * Mixing stays replicated (0.3 percent of the si32 wall at one rank, ledger
//!   MPI-M1-GATE4; distributing it would add a grid allgatherv per step that costs
//!   about what it saves).
//!
//! Determinism (ADR 0041): every reduction is a fixed-order local sum followed by
//! the binomial tree_sum, so a fixed rank count is bit-reproducible; different
//! rank counts differ by rounding only. Every collective here is entered by all
//! ranks under the same replicated conditions (no data-dependent branch).
use crate::parallel::SharedComm;
use crate::{Functional, Pw92Prefactor, ScfError, ScfOptions};
use dft_mpi::{Collective, block_counts, block_displacements};
use dft_physics::data::{Array3, Hartree};
use dft_physics::{PhysicsError, xc};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::time::Instant;

/// The M2 communicator when the M2 path is requested and two or more ranks run.
pub(crate) fn active() -> Option<SharedComm> {
    SharedComm::active().filter(|_| crate::parallel::m2_enabled())
}

/// `BAANDI_M2_BATCH=off` restores the per-pair FFT transposes of M2b/M2c (A/B).
pub(crate) fn batched_transposes() -> bool {
    dft_numerics::switches::on(dft_numerics::switches::Id::M2_BATCH)
}

/// Per-band plane-local density contributions of the latest M2 density build.
pub(crate) struct M2Density {
    /// bands x plane points, [band][plane][row][z].
    pub(crate) contributions: Vec<f64>,
    pub(crate) bands: usize,
    pub(crate) points: usize,
    /// Row-major grid offset of every plane-local point.
    pub(crate) offsets: Vec<usize>,
    /// True when the latest SCF map built its density here (else the serial
    /// `BandDensity` holds the current contributions).
    pub(crate) valid: bool,
    pub(crate) partial: Vec<f64>,
}

impl M2Density {
    pub(crate) fn new(offsets: Vec<usize>) -> Self {
        Self {
            contributions: Vec::new(),
            bands: 0,
            points: offsets.len(),
            offsets,
            valid: false,
            partial: Vec::new(),
        }
    }

    /// out[b] = sum over grid points of band b's contribution times `field`
    /// (full row-major grid), as `BandDensity::weighted_integrals`; local sums in
    /// plane order plus one tree_sum. Collective.
    pub(crate) fn weighted_integrals<K: Collective>(
        &mut self,
        comm: &K,
        field: &[f64],
        out: &mut [f64],
    ) -> Result<(), ScfError> {
        let start = Instant::now();
        self.partial.clear();
        self.partial.resize(out.len(), 0.);
        for (b, p) in self.partial.iter_mut().enumerate().take(self.bands) {
            let c = &self.contributions[b * self.points..(b + 1) * self.points];
            *p = c
                .iter()
                .zip(&self.offsets)
                .map(|(c, &o)| c * field[o])
                .sum();
        }
        comm.tree_sum_f64(&self.partial, out)
            .map_err(|e| ScfError::Numerical(e.to_string()))?;
        add(Stage::ScreenIntegrals, start);
        count(Count::ScreenIntegrals);
        Ok(())
    }
}

/// Distributed LDA XC of the replicated density `n` when the functional and
/// options allow it (LDA PZ or PW, no core density, no density floor, modified
/// PW92 prefactor); None means: use the replicated evaluation. Collective.
pub(crate) fn xc<K: Collective>(
    comm: &K,
    options: &ScfOptions,
    n: &Array3<f64>,
    core: Option<&Array3<f64>>,
    floor: Option<f64>,
    volume: f64,
) -> Option<Result<xc::XcGrid, PhysicsError>> {
    if core.is_some()
        || floor.is_some()
        || options.pw92_prefactor != Pw92Prefactor::Modified
        || !matches!(options.functional, Functional::LdaPz | Functional::LdaPw)
        || n.as_slice().is_none()
    {
        return None;
    }
    Some(xc_lda(comm, options.functional, n, volume))
}

fn xc_lda<K: Collective>(
    comm: &K,
    functional: Functional,
    n: &Array3<f64>,
    volume: f64,
) -> Result<xc::XcGrid, PhysicsError> {
    let start = Instant::now();
    if !(volume.is_finite() && volume > 0.) || n.is_empty() {
        return Err(PhysicsError::Shape);
    }
    let all = n.as_slice().ok_or(PhysicsError::Shape)?;
    let counts = block_counts(all.len(), comm.size());
    let displs = block_displacements(&counts);
    let (lo, len) = (displs[comm.rank()], counts[comm.rank()]);
    let local = &all[lo..lo + len];
    let mut potential_local = vec![0.; len];
    // [energy sum, integral sum, failure flag] of this rank's points.
    let mut sums = [0.; 3];
    let ok = (|| -> Result<(), PhysicsError> {
        if len == 0 {
            return Ok(());
        }
        match functional {
            Functional::LdaPz => {
                for (v, &n) in potential_local.iter_mut().zip(local) {
                    let p = xc::lda_pz(n)?;
                    let n = n.max(0.0);
                    *v = p.v_xc();
                    sums[0] += n * p.eps_xc();
                    sums[1] += n * *v;
                }
            }
            _ => {
                // Volume = point count makes the library weight exactly 1.
                let block = Array3::from_shape_vec((1, 1, len), local.to_vec())
                    .map_err(|_| PhysicsError::Shape)?;
                let r = dft_physics::spin::unpolarized_pw(&block, len as f64)?;
                potential_local.copy_from_slice(r.potential.as_slice().unwrap());
                sums[0] = r.energy.0;
                sums[1] = r.potential_integral.0;
            }
        }
        Ok(())
    })();
    if ok.is_err() {
        sums = [0., 0., 1.];
    }
    let mut total = [0.; 3];
    comm.tree_sum_f64(&sums, &mut total)
        .map_err(|e| PhysicsError::Xc(e.to_string()))?;
    if total[2] != 0. {
        // Same error on every rank (the local one where it occurred).
        return Err(ok
            .err()
            .unwrap_or_else(|| PhysicsError::Xc("invalid density on another rank".into())));
    }
    let mut potential = Array3::zeros(n.raw_dim());
    comm.allgatherv_f64(&potential_local, &counts, potential.as_slice_mut().unwrap())
        .map_err(|e| PhysicsError::Xc(e.to_string()))?;
    let weight = volume / all.len() as f64;
    let (energy, integral) = (total[0] * weight, total[1] * weight);
    if !energy.is_finite() || !integral.is_finite() || potential.iter().any(|v| !v.is_finite()) {
        return Err(PhysicsError::Unrepresentable);
    }
    add(Stage::Xc, start);
    count(Count::Xc);
    Ok(xc::XcGrid {
        potential,
        energy: Hartree(energy),
        potential_integral: Hartree(integral),
    })
}

#[derive(Clone, Copy)]
pub(crate) enum Stage {
    Density,
    Xc,
    NonlocalEnergy,
    ScreenIntegrals,
}
const STAGE_NAMES: [&str; 4] = ["density", "xc", "nonlocal_energy", "screen_integrals"];
static STAGE_NANOS: [AtomicU64; 4] = [const { AtomicU64::new(0) }; 4];
pub(crate) fn add(stage: Stage, start: Instant) {
    STAGE_NANOS[stage as usize].fetch_add(start.elapsed().as_nanos() as u64, Ordering::Relaxed);
}
#[derive(Clone, Copy)]
pub(crate) enum Count {
    Density,
    Xc,
    NonlocalEnergy,
    ScreenIntegrals,
}
const COUNT_NAMES: [&str; 4] = [
    "density_builds",
    "xc_evaluations",
    "nonlocal_energy_audits",
    "screen_integrals",
];
static COUNTS: [AtomicUsize; 4] = [const { AtomicUsize::new(0) }; 4];
pub(crate) fn count(c: Count) {
    COUNTS[c as usize].fetch_add(1, Ordering::Relaxed);
}

/// Engine-level M2d counters and stage seconds of this process.
pub fn stats() -> dft_physics::distributed_real_gamma::M2Counters {
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
