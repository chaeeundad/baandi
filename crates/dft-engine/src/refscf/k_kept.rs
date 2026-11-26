//! Exact per-k X/HX/potential cache. Verification borrows without publishing.
use super::{kept_hx::Report, numerical};
use crate::{model::ScfError, refstart::KStartState};
use dft_mpi::Collective;
use dft_numerics::linalg::ComplexMatrix;
use dft_physics::{data::Complex64, distributed_complex::DistributedComplexHamiltonian};
use std::sync::Arc;
struct Entry {
    identity: Arc<()>,
    epoch: u64,
    x: ComplexMatrix,
    hx: ComplexMatrix,
    potential: Vec<f64>,
    /// GK17: leading columns of x whose real-space planes are in `slabs[k]`
    /// (stored by the density pass of the map that produced x); 0 = none.
    slab_bands: usize,
}
pub(super) struct KeptHx {
    enabled: bool,
    preserve: bool,
    entries: Vec<Option<Entry>>,
    hits: usize,
    columns: usize,
    /// GK17 `K_KEPT_SLABS`: per-rank slab budget in bytes (None = off).
    slab_budget: Option<usize>,
    /// Per k: planes of the leading bands of the entry's x (operator layout).
    slabs: Vec<Vec<Complex64>>,
    slab_high: usize,
}
impl KeptHx {
    pub fn new(enabled: bool) -> Self {
        use dft_numerics::switches::{Id, number, on};
        let slab_budget = on(Id::K_KEPT_SLABS).then(|| {
            number(Id::K_KEPT_SLAB_MB).map_or(0, |mb| (mb as usize).saturating_mul(1 << 20))
        });
        Self {
            enabled,
            preserve: false,
            entries: Vec::new(),
            hits: 0,
            columns: 0,
            slab_budget,
            slabs: Vec::new(),
            slab_high: 0,
        }
    }
    /// GK17: slab sink for the density pass of k: the buffer and the number of
    /// leading bands to keep. `band_bytes` must be the same on every rank that
    /// shares the k solve (the caller passes 16 x grid points / ranks rounded
    /// up), and `local_k` the k count of this rank, so the band count (which
    /// decides the transforms, hence the collectives, of the next update) is
    /// rank-independent. Nothing is kept while verification borrows the cache.
    pub fn slab_sink(
        &mut self,
        k: usize,
        local_k: usize,
        bands: usize,
        band_bytes: usize,
        plane_points: usize,
    ) -> (&mut [Complex64], usize) {
        let keep = match self.slab_budget {
            Some(budget) if self.enabled && !self.preserve && band_bytes > 0 => {
                bands.min(budget / local_k.max(1) / band_bytes)
            }
            _ => 0,
        };
        if self.slabs.len() <= k {
            self.slabs.resize_with(k + 1, Vec::new);
        }
        let store = &mut self.slabs[k];
        store.resize(keep * plane_points, Complex64::default());
        let held: usize = self.slabs.iter().map(|v| 16 * v.len()).sum();
        self.slab_high = self.slab_high.max(held);
        (&mut self.slabs[k][..], keep)
    }
    /// `retain` for an x whose leading `slab_bands` planes the density pass
    /// just stored through [`Self::slab_sink`] for the same k.
    pub fn retain_with_slabs(
        &mut self,
        k: usize,
        state: &KStartState,
        x: &ComplexMatrix,
        hx: &ComplexMatrix,
        potential: &[f64],
        slab_bands: usize,
    ) {
        self.retain(k, state, x, hx, potential);
        if let Some(Some(e)) = self.entries.get_mut(k) {
            e.slab_bands = slab_bands;
        }
    }
    /// GF1: the density pass of the spin/FD map stored the leading
    /// `slab_bands` planes of k's retained x (retained before the density).
    pub fn mark_slabs(&mut self, k: usize, slab_bands: usize) {
        if let Some(Some(e)) = self.entries.get_mut(k) {
            e.slab_bands = slab_bands;
        }
    }
    pub fn preserve_input(&mut self, p: bool) {
        self.preserve = p;
    }
    pub fn invalidate(&mut self) {
        for e in &mut self.entries {
            *e = None;
        }
    }
    pub fn release_storage(&mut self) {
        self.entries.clear();
        self.slabs.clear();
    }
    pub fn prepare<K: Collective>(
        &mut self,
        k: usize,
        comm: &K,
        op: &mut DistributedComplexHamiltonian,
        x: &ComplexMatrix,
        potential: &[f64],
    ) -> Result<Option<ComplexMatrix>, ScfError> {
        if !self.enabled {
            return Ok(None);
        }
        if self.entries.len() <= k {
            self.entries.resize_with(k + 1, || None);
        }
        let old = self.entries[k].take();
        let hit = old.as_ref().is_some_and(|e| {
            Arc::ptr_eq(&e.identity, &op.identity())
                && e.x.rows() == x.rows()
                && e.x.cols() == x.cols()
                && e.x.as_slice().iter().zip(x.as_slice()).all(|(a, b)| {
                    a.re.to_bits() == b.re.to_bits() && a.im.to_bits() == b.im.to_bits()
                })
                && e.potential.len() == potential.len()
                && op.owned_epoch().is_some_and(|v| v > e.epoch)
        });
        let delta: Vec<_> = match (&old, hit) {
            (Some(e), true) => potential
                .iter()
                .zip(&e.potential)
                .map(|(a, b)| a - b)
                .collect(),
            _ => Vec::new(),
        };
        // The zero-change count (GK10B redo on the same potential) rides in
        // the same sum as the miss count.
        let local_nonzero = f64::from(delta.iter().any(|d| *d != 0.));
        let mut flags = [0.; 2];
        comm.tree_sum_f64(&[f64::from(!hit), local_nonzero], &mut flags)
            .map_err(numerical)?;
        if flags[0] != 0. {
            return Ok(None);
        }
        let e = old.expect("global hit");
        let mut hx = e.hx.clone();
        // A zero change on every rank adds nothing: skip its transforms.
        if flags[1] != 0. {
            // GK17: the planes stored with this exact x skip the inverse
            // transform of the leading bands (same bits; rank-independent
            // band count, see `slab_sink`).
            let kept = e.slab_bands;
            if kept > 0
                && op.kept_planes_exact()
                && self.slabs.get(k).is_some_and(|v| !v.is_empty())
            {
                op.add_local_change_kept(comm, x, &delta, &mut hx, &self.slabs[k], kept)
                    .map_err(numerical)?;
            } else {
                op.add_local_change(comm, x, &delta, &mut hx)
                    .map_err(numerical)?;
            }
        }
        if self.preserve {
            self.entries[k] = Some(e);
        }
        self.hits += 1;
        self.columns += x.cols();
        Ok(Some(hx))
    }
    pub fn retain(
        &mut self,
        k: usize,
        state: &KStartState,
        x: &ComplexMatrix,
        hx: &ComplexMatrix,
        potential: &[f64],
    ) {
        if !self.enabled || self.preserve {
            return;
        }
        if self.entries.len() <= k {
            self.entries.resize_with(k + 1, || None);
        }
        let op = &state.blocks[k].operator;
        self.entries[k] = Some(Entry {
            identity: op.identity(),
            epoch: op.owned_epoch().expect("potential"),
            x: x.clone(),
            hx: hx.clone(),
            potential: potential.to_vec(),
            slab_bands: 0,
        });
    }
    pub fn report<K: Collective>(&self, comm: &K) -> Result<Report, ScfError> {
        Ok(Report {
            cache_hits: self.hits,
            full_h_starting_columns_replaced: self.columns,
            delta_v_columns: self.columns,
            // GK17: this rank's own high-water (other ranks' slots stay 0;
            // no collective is spent on a report field).
            retained_slab_bytes_per_rank: (0..comm.size())
                .map(|r| if r == comm.rank() { self.slab_high } else { 0 })
                .collect(),
        })
    }
}

impl KeptHx {
    pub(super) fn checkpoint(
        &mut self,
        io: &mut super::checkpoint::Exchange<'_>,
        state: &KStartState,
    ) -> Result<(), ScfError> {
        use super::checkpoint::invalid;
        io.same(&self.enabled)?;
        io.value(&mut self.preserve)?;
        if self.preserve {
            return Err(invalid("k HX capture during verification"));
        }
        let mut data = if io.loading() {
            Vec::new()
        } else {
            self.entries
                .iter()
                .map(|e| {
                    e.as_ref()
                        .map(|e| (e.x.clone(), e.hx.clone(), (e.potential.clone(), e.epoch)))
                })
                .collect()
        };
        io.value(&mut data)?;
        if data.len() > state.blocks.len() {
            return Err(invalid("k HX block count"));
        }
        if io.loading() {
            self.entries.clear();
        }
        for (k, entry) in data.into_iter().enumerate() {
            let op = &state.blocks[k].operator;
            let restored = if let Some((x, hx, (potential, epoch))) = entry {
                let bands = state.orbitals.0[k].cols();
                if x.rows() != op.local_rows().len()
                    || x.cols() != bands
                    || hx.rows() != x.rows()
                    || hx.cols() != bands
                    || potential.len() != op.plane_points()
                    || x.as_slice()
                        .iter()
                        .chain(hx.as_slice())
                        .any(|v| !v.re.is_finite() || !v.im.is_finite())
                {
                    return Err(invalid("k HX shape/values"));
                }
                io.finite(&potential)?;
                Some(Entry {
                    identity: op.identity(),
                    x,
                    hx,
                    potential,
                    epoch,
                    // Planes are not checkpointed: the restarted map runs
                    // both transforms (same bits).
                    slab_bands: 0,
                })
            } else {
                None
            };
            if io.loading() {
                self.entries.push(restored);
            }
        }
        io.value(&mut self.hits)?;
        io.value(&mut self.columns)?;
        Ok(())
    }
}
