//! P4-MEM: allocation-boundary census of accessible storage. No allocator hook,
//! I/O, JSON, or sampling in elementwise kernels. Each counter is an independent
//! component maximum; adding them does NOT produce a simultaneous peak or RSS.
use super::numerical;
use crate::ScfError;
use dft_mpi::Collective;
use std::{cell::Cell, marker::PhantomData, rc::Rc};

pub const NAMES: [&str; 18] = [
    "start_retained_vectors",
    "start_vectors_released",
    "driver_reciprocal_vectors",
    "driver_plane_vectors",
    "density_plan_vectors",
    "metric_vectors",
    "broyden_history_vectors",
    "broyden_metadata_and_reduced",
    "prediction_vectors",
    "map_potential_vectors",
    "map_density_vectors",
    "field_cache_vectors",
    "verification_reciprocal_vector",
    "kept_slabs",
    "kept_entry",
    "kept_delta_vectors",
    "orbital_gather_buffers",
    "delivered_arrays",
];
#[derive(Clone, Copy)]
#[repr(usize)]
pub(crate) enum Component {
    Start,
    ReleasedStart,
    DriverReciprocal,
    DriverPlanes,
    DensityPlan,
    Metric,
    History,
    Mixer,
    Prediction,
    MapPotential,
    MapDensity,
    FieldCache,
    Verification,
    Slabs,
    KeptEntry,
    KeptDelta,
    Gather,
    Delivered,
}
thread_local! {
    static ACTIVE: Cell<Option<[usize; NAMES.len()]>> = const { Cell::new(None) };
}
pub(crate) struct Scope {
    previous: Option<[usize; NAMES.len()]>,
    solver: dft_solvers::ref_davidson::work::MemoryScope,
    _thread: PhantomData<Rc<()>>,
}
#[derive(Debug)]
pub struct Breakdown {
    /// Rank order; each row follows NAMES. Vec capacity in bytes, excluding the
    /// inline Vec header; matrices count initialized backing elements.
    pub components_per_rank: Vec<Vec<usize>>,
    /// V, HV, X, HX, E, HE, reduced, vectors, returned, cached (bytes).
    pub davidson_per_rank: Vec<Vec<usize>>,
}
pub const DAVIDSON_NAMES: [&str; 10] = [
    "v",
    "hv",
    "x",
    "hx",
    "e",
    "he",
    "reduced_arrays",
    "vectors",
    "returned_backing",
    "cached_rows",
];
impl Scope {
    pub fn begin() -> Self {
        Self {
            previous: ACTIVE.with(|a| a.replace(Some([0; NAMES.len()]))),
            solver: dft_solvers::ref_davidson::work::MemoryScope::begin(),
            _thread: PhantomData,
        }
    }
    pub fn report<K: Collective>(&self, comm: &K) -> Result<Breakdown, ScfError> {
        let engine = ACTIVE.with(|a| a.get().expect("active census"));
        let s = self.solver.snapshot();
        let solver = [
            s.rows[0], s.rows[1], s.rows[2], s.rows[3], s.rows[4], s.rows[5], s.reduced, s.vectors,
            s.returned, s.cached,
        ];
        let local: Vec<_> = engine.into_iter().chain(solver).map(|v| v as f64).collect();
        let width = local.len();
        let mut all = vec![0.; width * comm.size()];
        comm.allgatherv_f64(&local, &vec![width; comm.size()], &mut all)
            .map_err(numerical)?;
        let rows: Vec<Vec<usize>> = all
            .chunks_exact(width)
            .map(|r| r.iter().map(|&v| v as usize).collect())
            .collect();
        Ok(Breakdown {
            components_per_rank: rows.iter().map(|r| r[..NAMES.len()].to_vec()).collect(),
            davidson_per_rank: rows.iter().map(|r| r[NAMES.len()..].to_vec()).collect(),
        })
    }
}
impl Drop for Scope {
    fn drop(&mut self) {
        ACTIVE.with(|a| a.set(self.previous));
    }
}
pub(crate) fn bytes<T>(v: &Vec<T>) -> usize {
    v.capacity() * size_of::<T>()
}
pub(crate) fn observe(component: Component, bytes: usize) {
    ACTIVE.with(|a| {
        if let Some(mut maxima) = a.get() {
            let i = component as usize;
            maxima[i] = maxima[i].max(bytes);
            a.set(Some(maxima));
        }
    });
}
