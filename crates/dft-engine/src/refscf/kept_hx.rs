//! Run-stage-owned starting products. No cache crosses a retry, verification
//! boundary or cutoff handoff. The fixed operator token includes the basis and
//! both row/plane layouts; orbital equality is bitwise, never a checksum.
use super::numerical;
use crate::ScfError;
use dft_mpi::Collective;
use dft_numerics::RealMatrix;
use dft_physics::{data::Complex64, distributed_real_gamma::DistributedRealHamiltonian};
use std::sync::Arc;

/// Maximum retained real-space pair slabs per rank (other arrays excluded).
/// 512 MiB lets the whole start block of the Si64/Si128 MPI8 cells hit; a smaller
/// limit keeps only as many pairs as fit (partial reuse, recorded in kept_hx).
pub const SLAB_LIMIT_BYTES: usize = 512 * 1024 * 1024;

#[derive(Clone, Debug, Default)]
pub struct Report {
    pub cache_hits: usize,
    pub full_h_starting_columns_replaced: usize,
    pub delta_v_columns: usize,
    /// High-water retained slab payload per rank, bytes (not RSS).
    pub retained_slab_bytes_per_rank: Vec<usize>,
}
impl Report {
    pub(super) fn merge(&mut self, other: &Self) {
        self.cache_hits += other.cache_hits;
        self.full_h_starting_columns_replaced += other.full_h_starting_columns_replaced;
        self.delta_v_columns += other.delta_v_columns;
        for (a, b) in self
            .retained_slab_bytes_per_rank
            .iter_mut()
            .zip(&other.retained_slab_bytes_per_rank)
        {
            *a = (*a).max(*b);
        }
    }
}
struct Entry {
    identity: Arc<()>,
    x: RealMatrix,
    hx: RealMatrix,
    potential: Vec<f64>,
    epoch: u64,
    slabs: Vec<Complex64>,
    pairs: usize,
}

/// A single-use cache bound to exact returned Ritz vectors. Disabling it keeps
/// the ordinary density and fresh-H path. Potential arguments are owned planes
/// in the operator's layout and must already have been installed on the operator.
pub struct KeptHx {
    enabled: bool,
    entry: Option<Entry>,
    hits: usize,
    columns: usize,
    slab_bytes: usize,
}
impl KeptHx {
    pub fn new(enabled: bool) -> Self {
        Self {
            enabled,
            entry: None,
            hits: 0,
            columns: 0,
            slab_bytes: 0,
        }
    }
    pub fn invalidate(&mut self) {
        self.entry = None;
    }
    /// Consume the previous entry. One collective agrees eligibility before any
    /// rank enters the delta kernel, including ranks with no rows or planes.
    pub fn prepare<K: Collective>(
        &mut self,
        comm: &K,
        operator: &mut DistributedRealHamiltonian,
        x: &RealMatrix,
        potential: &[f64],
    ) -> Result<Option<RealMatrix>, ScfError> {
        if !self.enabled {
            return Ok(None);
        }
        let old = self.entry.take();
        let hit = old.as_ref().is_some_and(|e| {
            Arc::ptr_eq(&e.identity, &operator.kept_identity())
                && e.x.rows() == x.rows()
                && e.x.cols() == x.cols()
                && e.x
                    .as_slice()
                    .iter()
                    .zip(x.as_slice())
                    .all(|(a, b)| a.to_bits() == b.to_bits())
                && e.potential.len() == potential.len()
                && potential.len() == operator.plane_points()
                && potential.iter().all(|v| v.is_finite())
                && operator.owned_epoch().is_some_and(|epoch| epoch > e.epoch)
        });
        let mut misses = [0.];
        comm.tree_sum_f64(&[f64::from(!hit)], &mut misses)
            .map_err(numerical)?;
        if misses[0] != 0. {
            return Ok(None);
        }
        let mut e = old.expect("agreed cache hit");
        let delta_planes: Vec<_> = potential
            .iter()
            .zip(&e.potential)
            .map(|(a, b)| a - b)
            .collect();
        // The existing kernel consumes full-grid indexing but reads owned planes
        // only. Unowned slots remain NaN to expose accidental remote reads.
        let mut delta = vec![f64::NAN; operator.potential().len()];
        operator.write_planes(&delta_planes, &mut delta);
        let mut change = RealMatrix::zeros(x.rows(), x.cols()).map_err(numerical)?;
        operator
            .apply_local_change(comm, x, &e.slabs, e.pairs, &delta, &mut change)
            .map_err(numerical)?;
        for (h, d) in e.hx.as_mut_slice().iter_mut().zip(change.as_slice()) {
            *h += d;
        }
        self.hits += 1;
        self.columns += x.cols();
        Ok(Some(e.hx))
    }
    /// Build density while retaining bounded leading pair slabs from the exact
    /// returned X. An odd last band (or pairs beyond the cap) is transformed
    /// again by the delta kernel, so all B initial full-H columns are replaced.
    #[allow(clippy::too_many_arguments)]
    pub fn density<K: Collective>(
        &mut self,
        comm: &K,
        operator: &mut DistributedRealHamiltonian,
        x: &RealMatrix,
        hx: &RealMatrix,
        potential: &[f64],
        occupations: &[f64],
        volume: f64,
        out: &mut [f64],
    ) -> Result<(), ScfError> {
        self.invalidate();
        if !self.enabled {
            operator
                .band_density_sum(comm, x, occupations, volume, out)
                .map_err(numerical)?;
            return Ok(());
        }
        let dims = operator.grid_shape();
        let axis = operator.plane_layout().0;
        // A replicated upper bound ensures the same pair count on every rank.
        let max_points =
            dims[axis].div_ceil(comm.size()) * dims.iter().product::<usize>() / dims[axis];
        let pairs = (x.cols() / 2).min(SLAB_LIMIT_BYTES / (16 * max_points.max(1)));
        let mut slabs = vec![Complex64::default(); pairs * operator.plane_points()];
        operator
            .band_density_sum_keeping(comm, x, occupations, volume, out, &mut slabs)
            .map_err(numerical)?;
        self.slab_bytes = self.slab_bytes.max(slabs.len() * 16);
        self.entry = Some(Entry {
            identity: operator.kept_identity(),
            x: x.clone(),
            hx: hx.clone(),
            potential: potential.to_vec(),
            epoch: operator.owned_epoch().expect("owned potential"),
            slabs,
            pairs,
        });
        Ok(())
    }
    pub fn report<K: Collective>(&self, comm: &K) -> Result<Report, ScfError> {
        let mut local = vec![0.; comm.size()];
        local[comm.rank()] = self.slab_bytes as f64;
        let mut all = local.clone();
        comm.tree_sum_f64(&local, &mut all).map_err(numerical)?;
        Ok(Report {
            cache_hits: self.hits,
            full_h_starting_columns_replaced: self.columns,
            delta_v_columns: self.columns,
            retained_slab_bytes_per_rank: all.into_iter().map(|v| v as usize).collect(),
        })
    }
}
