//! Shared distributed orbital delivery and materialization.
use crate::model::ScfError;
use crate::profile::distributed;
use dft_mpi::Collective;
use dft_numerics::{RealMatrix, linalg::ComplexMatrix};
use dft_physics::data::Complex64;
use dft_solvers::SolverError;
use std::sync::Arc;

thread_local! {
    static DELIVER_ORBITALS: std::cell::Cell<bool> = const { std::cell::Cell::new(true) };
}
/// LAZY-RESULT: tell the engine whether the caller on this thread will read
/// the delivered orbitals of the following SCF runs (`ScfResult::orbitals`,
/// `k_points[..].solution.eigenvectors`: forces, stress, orbital volumetric
/// outputs, checkpoints, warm starts). Default `true` (always gathered). With
/// `false` and `BAANDI_M2_LAZY_RESULT=on` (default) a result whose blocks are
/// lazy on the M2 path keeps them distributed (0 x B placeholders) until
/// `ScfResult::materialize_orbitals` (collective). Every rank must pass the
/// same value (a replicated decision, e.g. from the parsed input).
pub fn set_delivered_orbitals(needed: bool) {
    DELIVER_ORBITALS.with(|d| d.set(needed));
}
/// Whether the delivery point may keep lazy blocks (see `set_delivered_orbitals`).
pub(crate) fn lazy_result_allowed() -> bool {
    !DELIVER_ORBITALS.with(|d| d.get())
        && dft_numerics::switches::on(dft_numerics::switches::Id::M2_LAZY_RESULT)
}

/// MPI-LAZY-ORBITALS: the eigenvectors of one standard M2 solve held as this
/// rank's stick rows (ORC real coordinates, the operator's `local_rows` order),
/// the authoritative SCF state of the M2 path. The complex N x B block is
/// formed only by `materialize` (three allgathervs and the Hermitian completion
/// of the gathered path, so its bits equal the gathered block's) and cached
/// with the handle. Cloning shares the data. Every rank must call
/// `materialize` at the same point of the (replicated) control flow.
#[derive(Clone)]
pub struct LazyOrbitals(pub(crate) Arc<LazyInner>);
pub(crate) struct LazyInner {
    /// n_loc x bands, this rank's rows.
    pub(crate) x: RealMatrix,
    /// Global ORC row index of each local row.
    pub(crate) rows: Arc<[usize]>,
    /// Global ORC length N.
    pub(crate) n: usize,
    pub(crate) half: Arc<dft_basis::HalfSphere>,
    pub(crate) comm: crate::parallel::SharedComm,
    pub(crate) full: std::sync::OnceLock<ComplexMatrix>,
    // Injected fine rows materialize as +0 complex entries, exactly like
    // prolong_orbitals (ordinary Hermitian unpack makes mirror imaginary -0).
    pub(crate) handoff_zero_pairs: Option<Arc<[usize]>>,
}
impl std::fmt::Debug for LazyOrbitals {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("LazyOrbitals")
            .field("n", &self.0.n)
            .field("local_rows", &self.0.x.rows())
            .field("bands", &self.0.x.cols())
            .field("materialized", &self.0.full.get().is_some())
            .finish()
    }
}
/// LAZY-ASSEMBLE-FUSED: the complex N x B block from the gathered rank
/// segments (`flat`: per rank, its `counts[r]` local rows of every column,
/// column-major; `index`: the global ORC coordinate of each gathered row).
/// Every ORC coordinate goes straight to its two complex slots of the column
/// (the representative and its mirror), column by column so one column stays
/// in cache across the rank segments: the arithmetic of `HalfSphere::unpack`
/// (G = 0 real, x / sqrt2 for both parts, the mirror's imaginary part
/// negated), so the block is bit-identical to scattering into a global ORC
/// block first and unpacking it (the former path, `assemble_reference` in the
/// tests), without that N x B intermediate. Same errors: a length that is not
/// the half sphere's, any nonfinite gathered value.
pub(crate) fn assemble_lazy_block(
    half: &dft_basis::HalfSphere,
    counts: &[usize],
    index: &[f64],
    flat: &[f64],
    n: usize,
    bands: usize,
) -> Result<ComplexMatrix, ScfError> {
    let fail = |e: dft_basis::HalfSphereError| ScfError::Numerical(format!("lazy orbitals: {e}"));
    if n != half.full_len() {
        return Err(fail(dft_basis::HalfSphereError::Dimension));
    }
    if flat.iter().any(|x| !x.is_finite()) {
        return Err(fail(dft_basis::HalfSphereError::InvalidCoefficients));
    }
    let (plus, minus) = (half.plus(), half.minus());
    let s = std::f64::consts::SQRT_2;
    let mut full = ComplexMatrix::zeros(n, bands).map_err(SolverError::from)?;
    for j in 0..bands {
        let column = full.column_mut(j);
        let (mut at, mut offset) = (0, 0);
        for &c in counts {
            let rows = &index[offset..offset + c];
            let src = &flat[at + j * c..at + (j + 1) * c];
            for (&g, &v) in rows.iter().zip(src) {
                let g = g as usize;
                if g == 0 {
                    column[0] = Complex64::new(v, 0.);
                    continue;
                }
                let i = g.div_ceil(2);
                let (p, m) = (plus[i] as usize, minus[i] as usize);
                let w = v / s;
                if g % 2 == 1 {
                    column[p].re = w;
                    column[m].re = w;
                } else {
                    column[p].im = w;
                    column[m].im = -w;
                }
            }
            at += c * bands;
            offset += c;
        }
    }
    Ok(full)
}

impl LazyOrbitals {
    /// The caller has agreed shape/partition validity on every rank. Rows
    /// contain raw ORC coordinates; the existing warm-start round trip remains
    /// responsible for the exact complex-unpack/pack arithmetic.
    pub(crate) fn from_owned_handoff(
        x: RealMatrix,
        rows: Vec<usize>,
        half: dft_basis::HalfSphere,
        comm: crate::parallel::SharedComm,
        zero_pairs: Vec<usize>,
    ) -> Self {
        Self(Arc::new(LazyInner {
            n: half.full_len(),
            x,
            rows: rows.into(),
            half: Arc::new(half),
            comm,
            full: std::sync::OnceLock::new(),
            handoff_zero_pairs: Some(zero_pairs.into()),
        }))
    }

    pub fn bands(&self) -> usize {
        self.0.x.cols()
    }
    /// Global plane-wave count (rows of the materialized block).
    pub fn global_len(&self) -> usize {
        self.0.n
    }

    /// The full complex N x B block (gathered once, then cached). Collective:
    /// every rank must call it at the same point.
    pub fn materialize(&self) -> Result<&ComplexMatrix, ScfError> {
        if let Some(m) = self.0.full.get() {
            return Ok(m);
        }
        let start = dft_numerics::stage_clock::start();
        let inner = &*self.0;
        let comm = &inner.comm;
        let size = comm.size();
        let (n, nl, bands) = (inner.n, inner.x.rows(), inner.x.cols());
        let err = |e: dft_mpi::MpiError| ScfError::Numerical(format!("lazy orbitals: {e}"));
        let mut counts_f = vec![0.0; size];
        comm.allgatherv_f64(&[nl as f64], &vec![1; size], &mut counts_f)
            .map_err(err)?;
        let counts: Vec<usize> = counts_f.iter().map(|&c| c as usize).collect();
        if counts.iter().sum::<usize>() != n {
            return Err(ScfError::Numerical(
                "lazy orbitals: row partition does not cover the basis".into(),
            ));
        }
        let mine: Vec<f64> = inner.rows.iter().map(|&r| r as f64).collect();
        let mut index = vec![0.0; n];
        comm.allgatherv_f64(&mine, &counts, &mut index)
            .map_err(err)?;
        let block: Vec<usize> = counts.iter().map(|c| c * bands).collect();
        let mut flat = vec![0.0; n * bands];
        comm.allgatherv_f64(inner.x.as_slice(), &block, &mut flat)
            .map_err(err)?;
        let mut full = assemble_lazy_block(&inner.half, &counts, &index, &flat, n, bands)?;
        if let Some(pairs) = &inner.handoff_zero_pairs {
            for j in 0..bands {
                for &i in pairs.iter() {
                    full.column_mut(j)[inner.half.plus()[i] as usize] = Complex64::new(0., 0.);
                    full.column_mut(j)[inner.half.minus()[i] as usize] = Complex64::new(0., 0.);
                }
            }
        }
        let _ = inner.full.set(full);
        distributed::count(distributed::Count::LazyMaterializations);
        distributed::add(distributed::Stage::LazyMaterialize, start);
        Ok(inner.full.get().expect("set above"))
    }
    /// The full block, owned (gathered if not yet; moved out of the cache when
    /// this is the last handle, cloned otherwise). Collective like `materialize`.
    pub fn into_full(self) -> Result<ComplexMatrix, ScfError> {
        self.materialize()?;
        match Arc::try_unwrap(self.0) {
            Ok(inner) => Ok(inner.full.into_inner().expect("materialized")),
            Err(shared) => Ok(shared.full.get().expect("materialized").clone()),
        }
    }
}
