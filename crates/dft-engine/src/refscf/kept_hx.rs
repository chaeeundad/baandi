//! Run-stage-owned starting products. No cache crosses a retry or cutoff
//! handoff. Warm verification borrows the candidate products without publishing
//! its own products or density into the ordinary-map cache. The fixed operator
//! token includes the basis and both row/plane layouts; orbital equality is
//! bitwise, never a checksum.
use super::{
    memory::{Component as Mem, bytes, observe},
    numerical,
};
use crate::model::ScfError;
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
    pairs: usize,
}

/// A cache bound to exact returned Ritz vectors, consumed by ordinary maps
/// and borrowed by warm verification. Disabling it keeps
/// the ordinary density and fresh-H path. Potential arguments are owned planes
/// in the operator's layout and must already have been installed on the operator.
pub struct KeptHx {
    enabled: bool,
    recompute_slabs: bool,
    preserve_input: bool,
    entry: Option<Entry>,
    // Storage belongs to the stage, not the single-use entry. Invalidation
    // drops eligibility, while the next density overwrites the same slabs.
    slabs: Vec<Complex64>,
    hits: usize,
    columns: usize,
    slab_bytes: usize,
}
impl KeptHx {
    pub fn new(enabled: bool) -> Self {
        Self {
            enabled,
            recompute_slabs: dft_numerics::switches::on(
                dft_numerics::switches::Id::REFSCF_KEPT_RECOMPUTE,
            ),
            preserve_input: false,
            entry: None,
            slabs: Vec::new(),
            hits: 0,
            columns: 0,
            slab_bytes: 0,
        }
    }
    /// Verification may read the candidate cache but must not replace it with
    /// the check state: recovery still starts from the candidate Ritz block.
    pub(super) fn preserve_input(&mut self, preserve: bool) {
        self.preserve_input = preserve;
    }
    pub(super) fn release_storage(&mut self) {
        self.entry = None;
        self.slabs = Vec::new();
    }
    pub fn invalidate(&mut self) {
        self.entry = None;
    }
    /// Consume the previous entry unless verification is preserving it. A miss
    /// always invalidates it. One collective agrees eligibility before any
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
        let e = old.expect("agreed cache hit");
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
        observe(
            Mem::KeptDelta,
            bytes(&delta_planes) + bytes(&delta) + std::mem::size_of_val(change.as_slice()),
        );
        operator
            .apply_local_change(comm, x, &self.slabs, e.pairs, &delta, &mut change)
            .map_err(numerical)?;
        let mut hx = if self.preserve_input {
            let hx = e.hx.clone();
            self.entry = Some(e);
            hx
        } else {
            e.hx
        };
        for (h, d) in hx.as_mut_slice().iter_mut().zip(change.as_slice()) {
            *h += d;
        }
        self.hits += 1;
        self.columns += x.cols();
        Ok(Some(hx))
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
        if !self.preserve_input {
            self.invalidate();
        }
        if !self.enabled || self.preserve_input {
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
        // Memory/work tradeoff: retain exact X/HX and the old potential, but
        // regenerate the inverse planes inside apply_local_change. Zero pairs
        // uses the existing f64 fallback and makes density the plain path.
        let pairs = if self.recompute_slabs {
            0
        } else {
            (x.cols() / 2).min(SLAB_LIMIT_BYTES / (16 * max_points.max(1)))
        };
        self.slabs
            .resize(pairs * operator.plane_points(), Complex64::default());
        operator
            .band_density_sum_keeping(comm, x, occupations, volume, out, &mut self.slabs)
            .map_err(numerical)?;
        self.slab_bytes = self.slab_bytes.max(self.slabs.len() * 16);
        observe(Mem::Slabs, bytes(&self.slabs));
        self.entry = Some(Entry {
            identity: operator.kept_identity(),
            x: x.clone(),
            hx: hx.clone(),
            potential: potential.to_vec(),
            epoch: operator.owned_epoch().expect("owned potential"),
            pairs,
        });
        let entry = self.entry.as_ref().expect("installed entry");
        observe(
            Mem::KeptEntry,
            std::mem::size_of_val(entry.x.as_slice())
                + std::mem::size_of_val(entry.hx.as_slice())
                + bytes(&entry.potential),
        );
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

impl KeptHx {
    pub(super) fn checkpoint(
        &mut self,
        io: &mut super::checkpoint::Exchange<'_>,
        operator: &DistributedRealHamiltonian,
        bands: usize,
    ) -> Result<(), ScfError> {
        use super::checkpoint::invalid;
        io.same(&self.enabled)?;
        io.same(&self.recompute_slabs)?;
        io.value(&mut self.preserve_input)?;
        if self.preserve_input {
            return Err(invalid("HX capture during verification"));
        }
        let mut data = if io.loading() {
            None
        } else {
            self.entry.as_ref().map(|e| {
                (
                    e.x.clone(),
                    e.hx.clone(),
                    (e.potential.clone(), e.epoch, e.pairs),
                )
            })
        };
        io.value(&mut data)?;
        if let Some((x, hx, (potential, epoch, pairs))) = data {
            if x.rows() != operator.local_len()
                || x.cols() != bands
                || hx.rows() != x.rows()
                || hx.cols() != bands
                || potential.len() != operator.plane_points()
                || pairs > bands / 2
            {
                return Err(invalid("Gamma HX shape"));
            }
            io.finite(x.as_slice())?;
            io.finite(hx.as_slice())?;
            io.finite(&potential)?;
            if io.loading() {
                self.entry = Some(Entry {
                    identity: operator.kept_identity(),
                    x,
                    hx,
                    potential,
                    epoch,
                    pairs,
                });
            }
        } else if io.loading() {
            self.entry = None;
        }
        io.value(&mut self.slabs)?;
        io.value(&mut self.hits)?;
        io.value(&mut self.columns)?;
        io.value(&mut self.slab_bytes)?;
        if self
            .slabs
            .iter()
            .any(|v| !v.re.is_finite() || !v.im.is_finite())
            || self.entry.as_ref().is_some_and(|e| {
                e.pairs.checked_mul(operator.plane_points()) != Some(self.slabs.len())
            })
        {
            return Err(invalid("Gamma HX slabs"));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slab_storage_survives_consumption_and_invalidation() {
        let pp =
            crate::scf::pool_fixture::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
        let cell =
            dft_geometry::Cell::from_user_rows_bohr([[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]])
                .unwrap();
        let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
        let potentials = [crate::SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let options = crate::ScfOptions {
            cutoff_hartree: 1.,
            bands: Some(3),
            gamma_real: crate::GammaReal::On,
            ..Default::default()
        };
        dft_mpi::simulate(2, |comm| {
            let mut state = crate::refstart::start(
                crate::refstart::StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &options,
                },
                comm,
            )
            .unwrap();
            let x = &state.orbitals;
            let op = &mut state.operator;
            let v = vec![0.03; op.plane_points()];
            op.set_owned_potential(&v, 1).unwrap();
            let mut hx = RealMatrix::zeros(x.rows(), x.cols()).unwrap();
            op.apply_block(comm, x, &mut hx).unwrap();
            let mut density = vec![0.; v.len()];
            let mut kept = KeptHx::new(true);
            kept.density(
                comm,
                op,
                x,
                &hx,
                &v,
                &[2., 0., 0.],
                cell.volume_bohr3(),
                &mut density,
            )
            .unwrap();
            let pointer = kept.slabs.as_ptr();
            let capacity = kept.slabs.capacity();
            assert!(capacity > 0);
            let original: Vec<_> = density.iter().map(|v| v.to_bits()).collect();
            for epoch in 2..=4 {
                op.set_owned_potential(&v, epoch).unwrap();
                if epoch == 3 {
                    kept.invalidate();
                    assert!(kept.prepare(comm, op, x, &v).unwrap().is_none());
                } else {
                    assert!(kept.prepare(comm, op, x, &v).unwrap().is_some());
                }
                // Poison the reused storage: the next map must overwrite every
                // retained value before publishing eligibility again.
                kept.slabs.fill(Complex64::new(f64::NAN, f64::NAN));
                kept.density(
                    comm,
                    op,
                    x,
                    &hx,
                    &v,
                    &[2., 0., 0.],
                    cell.volume_bohr3(),
                    &mut density,
                )
                .unwrap();
                assert_eq!(kept.slabs.as_ptr(), pointer);
                assert_eq!(kept.slabs.capacity(), capacity);
                assert!(
                    kept.slabs
                        .iter()
                        .all(|v| v.re.is_finite() && v.im.is_finite())
                );
                assert_eq!(
                    density.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                    original
                );
            }
        });
    }
}
