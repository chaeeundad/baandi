//! Exact per-k X/HX/potential cache. Verification borrows without publishing.
use super::{kept_hx::Report, numerical};
use crate::{model::ScfError, refstart::KStartState};
use dft_mpi::Collective;
use dft_numerics::linalg::ComplexMatrix;
use dft_physics::distributed_complex::DistributedComplexHamiltonian;
use std::sync::Arc;
struct Entry {
    identity: Arc<()>,
    epoch: u64,
    x: ComplexMatrix,
    hx: ComplexMatrix,
    potential: Vec<f64>,
}
pub(super) struct KeptHx {
    enabled: bool,
    preserve: bool,
    entries: Vec<Option<Entry>>,
    hits: usize,
    columns: usize,
}
impl KeptHx {
    pub fn new(enabled: bool) -> Self {
        Self {
            enabled,
            preserve: false,
            entries: Vec::new(),
            hits: 0,
            columns: 0,
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
        let mut misses = [0.];
        comm.tree_sum_f64(&[f64::from(!hit)], &mut misses)
            .map_err(numerical)?;
        if misses[0] != 0. {
            return Ok(None);
        }
        let e = old.expect("global hit");
        let delta: Vec<_> = potential
            .iter()
            .zip(&e.potential)
            .map(|(a, b)| a - b)
            .collect();
        let mut hx = e.hx.clone();
        op.add_local_change(comm, x, &delta, &mut hx)
            .map_err(numerical)?;
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
        });
    }
    pub fn report<K: Collective>(&self, comm: &K) -> Result<Report, ScfError> {
        Ok(Report {
            cache_hits: self.hits,
            full_h_starting_columns_replaced: self.columns,
            delta_v_columns: self.columns,
            retained_slab_bytes_per_rank: vec![0; comm.size()],
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
