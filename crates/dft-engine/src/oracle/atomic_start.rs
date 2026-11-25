//! P45A atomic Fourier span and first-density Rayleigh-Ritz initialization.
use crate::atomic_start::{BLOCK_PIVOT_GUARD, cholesky_upper, numerical, upper_inverse};
pub use crate::model::{AtomicStartReport, StartWavefunctions};
use crate::model::{ScfError, ScfOptions};
use dft_basis::HalfSphere;
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_geometry::Cell;
use dft_mpi::Collective;
use dft_numerics::linalg::{MatrixOp, gemm, hermitian_eigh, orthonormalize};
use dft_numerics::switches::{InitProfile, RunCounter as IC};
use dft_physics::{
    data::*,
    hamiltonian::{self, Hamiltonian, Species},
};
use dft_solvers::LinearOperator;
use dft_upf::radial::{projector_transform_with_rule, real_spherical_harmonic};

pub(super) struct AtomicSpecies {
    mesh: dft_upf::radial::RadialMesh,
    orbitals: Vec<dft_upf::AtomicWavefunction>,
    positions_bohr: Vec<[f64; 3]>,
}
pub(super) fn prepare(species: &[Species<'_>]) -> Result<Vec<AtomicSpecies>, ScfError> {
    species
        .iter()
        .map(|s| {
            Ok(AtomicSpecies {
                mesh: dft_upf::radial::RadialMesh::new(
                    s.potential.mesh().r().to_vec(),
                    s.potential.mesh().rab().to_vec(),
                )
                .map_err(numerical)?,
                orbitals: s.potential.atomic_wavefunctions().unwrap().to_vec(),
                positions_bohr: s.positions_bohr.to_vec(),
            })
        })
        .collect()
}

pub(super) struct AtomicSpan {
    q: ComplexMatrix,
    report: AtomicStartReport,
}

/// `pool` (P48F, `BAANDI_PHASE_THREADS=on`): build the columns and run the
/// two-pass Gram-Schmidt on the SCF pool. Every column is computed by the
/// same elementwise arithmetic and `orthonormalize_pooled` reproduces the
/// serial bits, so the span is identical with and without the pool.
///
/// `split` (`BAANDI_M2_ATOMIC=split` on the M2 path): every rank builds only its
/// cyclic share of the columns and the two-pass Gram-Schmidt runs distributed
/// with the serial arithmetic (`orthonormalize_split`); the span is
/// bit-identical to the replicated one on every rank. Collective.
#[allow(clippy::too_many_arguments)]
pub(super) fn span(
    cell: &Cell,
    basis: &PlaneWaveBasis,
    grid: &FftGridSpec,
    species: &[AtomicSpecies],
    options: &ScfOptions,
    bands: usize,
    pool: Option<&dft_numerics::pool::WorkerPool>,
    split: Option<&crate::parallel::SharedComm>,
) -> Result<AtomicSpan, ScfError> {
    let profile = InitProfile::current();
    let _span = profile.map(|p| p.span(IC::InitAtomicSpanCalls, IC::InitAtomicSpanNanos));
    let npw = basis.vectors().len();
    let count = species
        .iter()
        .map(|s| s.positions_bohr.len() * s.orbitals.iter().map(|c| 2 * c.l + 1).sum::<usize>())
        .sum();
    if let Some(p) = profile {
        p.count(IC::InitAtomicRows, npw);
        p.count(IC::InitAtomicColumns, count);
        p.count(IC::InitAtomicInputBytes, 16 * npw * count);
        let [x, y, z] = grid.dims();
        p.count(IC::InitAtomicGridPoints, x * y * z);
        p.count(IC::InitAtomicGridX, x);
        p.count(IC::InitAtomicGridY, y);
        p.count(IC::InitAtomicGridZ, z);
        p.count(IC::InitAtomicBands, bands);
    }
    let mut report = AtomicStartReport {
        atomic_orbitals: count,
        ..Default::default()
    };
    // Radial transforms per orbital (cached by |G|^2), then one descriptor per
    // column in the established column order: species, orbital, position, m.
    use crate::profile::{AtomicPart, atomic_part, clock};
    let part = clock();
    let mut radials = Vec::new();
    let mut columns = Vec::with_capacity(count);
    for sp in species {
        for c in &sp.orbitals {
            let mut cache = std::collections::BTreeMap::new();
            let mut radial = Vec::with_capacity(npw);
            for g in basis.momenta() {
                let key = g.norm_squared.to_bits();
                let value = if let Some(v) = cache.get(&key) {
                    *v
                } else {
                    report.radial_samples += c.radial_times_r.len();
                    let v = projector_transform_with_rule(
                        &sp.mesh,
                        c.l,
                        &c.radial_times_r,
                        g.norm_squared.sqrt(),
                        options.radial_even_rule,
                    )
                    .map_err(numerical)?
                        / cell.volume_bohr3().sqrt();
                    cache.insert(key, v);
                    v
                };
                radial.push(value);
            }
            let il = Complex64::new(0.0, -1.0).powu(c.l as u32);
            for &pos in &sp.positions_bohr {
                for m in -(c.l as isize)..=c.l as isize {
                    columns.push((radials.len(), c.l, m, pos, il));
                }
            }
            radials.push(radial);
        }
    }
    atomic_part(AtomicPart::Radial, part);
    // Rows first..first + column.len() of span column j (the same elementwise
    // arithmetic for any row range).
    let fill_rows = |j: usize, first: usize, column: &mut [Complex64]| -> Result<(), ScfError> {
        let (r, l, m, pos, il) = columns[j];
        let radial = &radials[r];
        let momenta = &basis.momenta()[first..first + column.len()];
        for ((out, g), &rad) in column.iter_mut().zip(momenta).zip(&radial[first..]) {
            let y = if g.norm_squared == 0.0 {
                if l == 0 {
                    1.0 / (4.0 * std::f64::consts::PI).sqrt()
                } else {
                    0.0
                }
            } else {
                real_spherical_harmonic(l, m, g.cartesian).map_err(numerical)?
            };
            let phase = hamiltonian::structure_factor(g.cartesian, &[pos])?;
            *out = il * phase * (rad * y);
        }
        Ok(())
    };
    let fill = |j: usize, column: &mut [Complex64]| fill_rows(j, 0, column);
    let part = clock();
    // BAANDI_M2_ATOMIC_BLOCK=on: row-distributed CholQR2 (rounding class);
    // None (any rank could not use it) continues with the split routine below.
    let distributed_block =
        split.filter(|_| npw > 0 && crate::parallel::tuning::tuning().atomic_block);
    // Local arithmetic reuse, not a new MPI mode. A disabled M2 split path
    // on a multi-rank job must not accidentally enter the local experiment.
    let local_block = npw > 0
        && split.is_none()
        && dft_numerics::switches::on(dft_numerics::switches::Id::ATOMIC_LOCAL_BLOCK)
        && crate::parallel::SharedComm::active().is_none();
    if distributed_block.is_some() || local_block {
        if let Some(p) = profile.filter(|_| local_block) {
            p.count(IC::InitAtomicLocalAttempts, 1);
        }
        let gamma = basis.k_fractional() == [0.; 3];
        if let Some(ortho) = orthonormalize_block(distributed_block, &fill_rows, npw, count, gamma)?
        {
            if let Some(p) = profile.filter(|_| local_block) {
                p.count(IC::InitAtomicLocalAccepted, 1);
            }
            atomic_part(AtomicPart::Orthonormalize, part);
            let part = clock();
            let span = finish_span(basis, grid, bands, ortho, npw, report);
            atomic_part(AtomicPart::Finish, part);
            return span;
        }
        if local_block {
            if let Some(p) = profile {
                p.count(IC::InitAtomicLocalFallbacks, 1);
            }
        } else {
            crate::profile::distributed::count_atomic_block_fallback();
        }
    }
    if let Some(comm) = split.filter(|_| npw > 0) {
        let (rank, ranks) = (comm.rank(), comm.size());
        let mine: Vec<usize> = (rank..count).step_by(ranks).collect();
        let mut own = vec![Complex64::default(); mine.len() * npw];
        // A fill failure is carried as a nonfinite column so that every rank
        // keeps entering the same broadcasts and all of them return the error.
        for (c, &j) in own.chunks_exact_mut(npw).zip(&mine) {
            if fill(j, c).is_err() {
                c[0] = Complex64::new(f64::NAN, 0.0);
            }
        }
        atomic_part(AtomicPart::Fill, part);
        let part = clock();
        let ortho = orthonormalize_split(comm, own, &mine, npw, count)?;
        atomic_part(AtomicPart::Orthonormalize, part);
        let part = clock();
        let span = finish_span(basis, grid, bands, ortho, npw, report);
        atomic_part(AtomicPart::Finish, part);
        return span;
    }
    if let Some(p) = profile {
        p.count(
            IC::InitAtomicMgsMacs,
            2 * npw * count * count.saturating_sub(1),
        );
    }
    let mut raw = ComplexMatrix::zeros(npw, count).map_err(numerical)?;
    let ortho = match pool {
        Some(pool) if npw > 0 => {
            pool.install(|| {
                dft_numerics::pool::for_each_chunk_mut(raw.as_mut_slice(), npw, &|first, column| {
                    fill(first / npw, column)
                })
            })?;
            atomic_part(AtomicPart::Fill, part);
            let part = clock();
            let mgs = profile.map(|p| p.span(IC::InitAtomicMgsCalls, IC::InitAtomicMgsNanos));
            let ortho =
                dft_numerics::linalg::orthonormalize_pooled(&raw, None, pool).map_err(numerical)?;
            drop(mgs);
            atomic_part(AtomicPart::Orthonormalize, part);
            ortho
        }
        _ => {
            for j in 0..count {
                fill(j, raw.column_mut(j))?;
            }
            atomic_part(AtomicPart::Fill, part);
            let part = clock();
            let mgs = profile.map(|p| p.span(IC::InitAtomicMgsCalls, IC::InitAtomicMgsNanos));
            let ortho = orthonormalize(&raw, None).map_err(numerical)?;
            drop(mgs);
            atomic_part(AtomicPart::Orthonormalize, part);
            ortho
        }
    };
    let part = clock();
    let span = finish_span(basis, grid, bands, ortho, npw, report);
    atomic_part(AtomicPart::Finish, part);
    span
}

/// Rank check and the plane-wave completion of the span (replicated).
fn finish_span(
    basis: &PlaneWaveBasis,
    grid: &FftGridSpec,
    bands: usize,
    ortho: dft_numerics::linalg::Orthonormalized,
    npw: usize,
    mut report: AtomicStartReport,
) -> Result<AtomicSpan, ScfError> {
    let count = report.atomic_orbitals;
    if !ortho.dropped.is_empty() {
        return Err(ScfError::InvalidInput(
            "rank deficient atomic starting orbitals",
        ));
    }
    let mut q = ortho.vectors;
    report.mgs_projection_macs = 2 * npw * count * count.saturating_sub(1);
    let half = if basis.k_fractional() == [0.; 3] {
        Some(HalfSphere::new(basis, grid).map_err(numerical)?)
    } else {
        None
    };
    if q.cols() < bands {
        q.reserve_columns(bands).map_err(numerical)?;
        let mut order: Vec<_> = (0..npw).collect();
        order.sort_by(|&a, &b| {
            basis.momenta()[a]
                .norm_squared
                .total_cmp(&basis.momenta()[b].norm_squared)
        });
        for i in order {
            // At Gamma, use cosine and sine columns, retaining kinetic ordering.
            // Both have explicit Hermitian completion, including odd band counts.
            if half.as_ref().is_some_and(|h| i > h.mirror()[i] as usize) {
                continue;
            }
            let components = if half.is_some() && i != 0 { 2 } else { 1 };
            for component in 0..components {
                let mut e = ComplexMatrix::zeros(npw, 1).map_err(numerical)?;
                if let Some(h) = &half
                    && i != 0
                {
                    let v = if component == 0 {
                        Complex64::new(std::f64::consts::FRAC_1_SQRT_2, 0.)
                    } else {
                        Complex64::new(0., std::f64::consts::FRAC_1_SQRT_2)
                    };
                    e[(i, 0)] = v;
                    e[(h.mirror()[i] as usize, 0)] = v.conj();
                } else {
                    e[(i, 0)] = Complex64::new(1., 0.);
                }
                report.mgs_projection_macs += 4 * npw * q.cols();
                let extra = orthonormalize(&e, Some(&q)).map_err(numerical)?;
                if extra.vectors.cols() == 1 {
                    q.append_columns(&extra.vectors).map_err(numerical)?;
                }
                if q.cols() == bands {
                    break;
                }
            }
            if q.cols() == bands {
                break;
            }
        }
    }
    if q.cols() < bands {
        return Err(ScfError::InvalidInput("insufficient atomic starting rank"));
    }
    Ok(AtomicSpan { q, report })
}

/// [`orthonormalize`] (two-pass Gram-Schmidt, no `against`) of an n-column
/// block whose columns are owned cyclically by the ranks (`mine`: this rank's
/// global column indices, ascending; `own`: their values, column-major).
///
/// Bit-identical to the serial routine: every column receives exactly the
/// serial sequence of operations (normalize by its original norm; first-pass
/// projections against the accepted columns in acceptance order; second pass
/// the same; final norm, drop or accept). First passes run right-looking on
/// the owners as each accepted column arrives (one broadcast per column from
/// its owner: a status word and the normalized column); the second pass of a
/// column runs on its owner when its turn comes. Every rank receives every
/// accepted column, so the result is replicated. Collective: all ranks enter
/// the same n broadcasts; a nonfinite column is reported by its owner in the
/// status word and every rank returns the same error.
fn orthonormalize_split(
    comm: &crate::parallel::SharedComm,
    mut own: Vec<Complex64>,
    mine: &[usize],
    m: usize,
    n: usize,
) -> Result<dft_numerics::linalg::Orthonormalized, ScfError> {
    use dft_numerics::linalg::{DEPENDENCE_TOLERANCE, norm, subtract_projection};
    let ranks = comm.size();
    // Per own column: original norm (NaN marks a nonfinite column).
    let originals: Vec<f64> = own
        .chunks_exact_mut(m)
        .map(|v| {
            if v.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
                return f64::NAN;
            }
            let original = norm(v);
            if !original.is_finite() {
                return f64::NAN;
            }
            if original != 0.0 {
                for z in v.iter_mut() {
                    *z /= original;
                }
            }
            original
        })
        .collect();
    let mut data: Vec<Complex64> = Vec::with_capacity(m * n);
    let mut dropped = Vec::new();
    let mut count = 0;
    let mut buffer = vec![0.0; 1 + 2 * m];
    for j in 0..n {
        let owner = j % ranks;
        if owner == comm.rank() {
            let l = j / ranks;
            let v = &mut own[l * m..(l + 1) * m];
            buffer[0] = if originals[l].is_nan() {
                2.0
            } else if originals[l] == 0.0 {
                0.0
            } else {
                for k in 0..count {
                    subtract_projection(v, &data[k * m..(k + 1) * m]);
                }
                let remaining = norm(v);
                if !remaining.is_finite() {
                    2.0
                } else if remaining <= DEPENDENCE_TOLERANCE {
                    0.0
                } else {
                    for (b, z) in buffer[1..].as_chunks_mut::<2>().0.iter_mut().zip(v.iter()) {
                        let z = z / remaining;
                        b[0] = z.re;
                        b[1] = z.im;
                    }
                    1.0
                }
            };
        }
        comm.bcast_f64(&mut buffer, owner)
            .map_err(|e| ScfError::Numerical(e.to_string()))?;
        if buffer[0] == 2.0 {
            return Err(numerical(dft_numerics::linalg::LinalgError::NonFinite));
        }
        if buffer[0] == 0.0 {
            dropped.push(j);
            continue;
        }
        data.extend(
            buffer[1..]
                .as_chunks::<2>()
                .0
                .iter()
                .map(|b| Complex64::new(b[0], b[1])),
        );
        count += 1;
        let q = &data[(count - 1) * m..count * m];
        // Right-looking first pass of this rank's later columns.
        for (l, &g) in mine.iter().enumerate() {
            if g > j && originals[l].is_finite() && originals[l] != 0.0 {
                subtract_projection(&mut own[l * m..(l + 1) * m], q);
            }
        }
    }
    Ok(dft_numerics::linalg::Orthonormalized {
        vectors: ComplexMatrix::new(m, count, data).map_err(numerical)?,
        dropped,
    })
}

/// Local products of the block orthonormalization. At Gamma the span columns
/// are Hermitian-symmetric, so their overlaps and the triangular factors are
/// real: the complex rows x n block is used as the real 2 rows x n block of
/// its interleaved parts (Re A^H A = A_r^T A_r, A T = A_r T for real T), half
/// the work of the complex products.
enum BlockAlgebra {
    Complex,
    Real(dft_numerics::real::RealGemmWorkspace),
}
impl BlockAlgebra {
    fn new(gamma: bool, rows: usize, n: usize) -> Result<Self, ScfError> {
        if !gamma {
            return Ok(Self::Complex);
        }
        let pool = dft_numerics::pool::WorkerPool::new(1).map_err(numerical)?;
        dft_numerics::real::RealGemmWorkspace::with_pool(2 * rows, n, pool)
            .map(Self::Real)
            .map_err(numerical)
    }
    /// Values per Gram matrix in the reduction (real: n^2, complex: 2 n^2).
    fn gram_len(&self, n: usize) -> usize {
        match self {
            Self::Complex => 2 * n * n,
            Self::Real(_) => n * n,
        }
    }
    /// This rank's Gram contribution into `out` (`gram_len` values).
    fn gram(&mut self, a: &ComplexMatrix, out: &mut [f64]) -> Result<(), ScfError> {
        let (rows, n) = (a.rows(), a.cols());
        match self {
            Self::Complex => {
                let g = gemm(a, MatrixOp::Adjoint, a, MatrixOp::Identity).map_err(numerical)?;
                for (b, z) in out.as_chunks_mut::<2>().0.iter_mut().zip(g.as_slice()) {
                    b[0] = z.re;
                    b[1] = z.im;
                }
            }
            Self::Real(ws) => {
                let flat = interleaved(a.as_slice());
                let mut g = dft_numerics::real::RealMatrix::zeros(n, n).map_err(numerical)?;
                ws.multiply_slices_into(
                    dft_numerics::real::GemmSite::Other,
                    &flat,
                    2 * rows,
                    n,
                    2 * rows,
                    dft_numerics::real::RealOp::Transpose,
                    &flat,
                    n,
                    2 * rows,
                    g.as_mut(),
                )
                .map_err(numerical)?;
                out.copy_from_slice(g.as_slice());
            }
        }
        Ok(())
    }
    /// Reduced Gram values back to a complex n x n matrix.
    fn unpack(&self, values: &[f64]) -> Vec<Complex64> {
        match self {
            Self::Complex => values
                .as_chunks::<2>()
                .0
                .iter()
                .map(|b| Complex64::new(b[0], b[1]))
                .collect(),
            Self::Real(_) => values.iter().map(|&v| Complex64::new(v, 0.0)).collect(),
        }
    }
    /// A T for this rank's rows (T n x n; its real part only at Gamma).
    fn update(&mut self, a: &ComplexMatrix, t: &[Complex64]) -> Result<ComplexMatrix, ScfError> {
        let (rows, n) = (a.rows(), a.cols());
        let profile = InitProfile::current();
        let _span = profile.map(|p| {
            p.span(
                IC::InitAtomicBlockUpdateCalls,
                IC::InitAtomicBlockUpdateNanos,
            )
        });
        if let Some(p) = profile {
            p.count(
                IC::InitAtomicBlockUpdateFlops,
                (if matches!(self, Self::Real(_)) { 4 } else { 8 }) * rows * n * n,
            );
        }
        match self {
            Self::Complex => {
                let t = ComplexMatrix::new(n, n, t.to_vec()).map_err(numerical)?;
                gemm(a, MatrixOp::Identity, &t, MatrixOp::Identity).map_err(numerical)
            }
            Self::Real(ws) => {
                let t: Vec<f64> = t.iter().map(|z| z.re).collect();
                let mut q =
                    dft_numerics::real::RealMatrix::zeros(2 * rows, n).map_err(numerical)?;
                ws.multiply_slices_into(
                    dft_numerics::real::GemmSite::Other,
                    &interleaved(a.as_slice()),
                    2 * rows,
                    n,
                    2 * rows,
                    dft_numerics::real::RealOp::Identity,
                    &t,
                    n,
                    n,
                    q.as_mut(),
                )
                .map_err(numerical)?;
                ComplexMatrix::new(
                    rows,
                    n,
                    q.as_slice()
                        .as_chunks::<2>()
                        .0
                        .iter()
                        .map(|b| Complex64::new(b[0], b[1]))
                        .collect(),
                )
                .map_err(numerical)
            }
        }
    }
}

/// The interleaved (re, im) parts of a complex slice (copy; no unsafe view).
fn interleaved(z: &[Complex64]) -> Vec<f64> {
    z.iter().flat_map(|z| [z.re, z.im]).collect()
}

/// Global Gram matrix A^H A of row-distributed columns (`local`: this rank's
/// rows; None for a rank without rows) by one tree_sum; the last summed value
/// counts the ranks that flagged a failure (`failed`, or a local product
/// error). None when any rank failed (the same on every rank).
fn block_gram(
    comm: Option<&crate::parallel::SharedComm>,
    algebra: &mut BlockAlgebra,
    local: Option<&ComplexMatrix>,
    failed: bool,
    n: usize,
) -> Result<Option<Vec<Complex64>>, ScfError> {
    let profile = InitProfile::current();
    let _span = profile.map(|p| p.span(IC::InitAtomicBlockGramCalls, IC::InitAtomicBlockGramNanos));
    if let (Some(p), Some(a)) = (profile, local) {
        p.count(
            IC::InitAtomicBlockGramFlops,
            (if matches!(algebra, BlockAlgebra::Real(_)) {
                4
            } else {
                8
            }) * a.rows()
                * n
                * n,
        );
    }
    let len = algebra.gram_len(n);
    let mut send = vec![0.0; len + 1];
    if failed {
        send[len] = 1.0;
    } else if let Some(a) = local
        && algebra.gram(a, &mut send[..len]).is_err()
    {
        send[..len].fill(0.0);
        send[len] = 1.0;
    }
    let mut out = vec![0.0; send.len()];
    if let Some(comm) = comm {
        comm.tree_sum_f64(&send, &mut out)
            .map_err(|e| ScfError::Numerical(e.to_string()))?;
    } else {
        out.copy_from_slice(&send);
    }
    if out[len] != 0.0 || out.iter().any(|v| !v.is_finite()) {
        return Ok(None);
    }
    Ok(Some(algebra.unpack(&out[..len])))
}

/// Fills rows `first..first + column.len()` of span column `j`.
type FillRows<'a> = dyn Fn(usize, usize, &mut [Complex64]) -> Result<(), ScfError> + 'a;

/// `BAANDI_M2_ATOMIC_BLOCK=on`: orthonormal span by a row-distributed CholQR2.
/// Every rank fills its contiguous block of rows of all `n` columns (the
/// elementwise arithmetic of the replicated fill), then twice: local Gram
/// GEMM, one tree_sum, replicated Cholesky and triangular inverse, local
/// update Q = A R^-1 (the first pass on the column-normalized Gram, as the
/// Gram-Schmidt routine normalizes each column first). One allgatherv returns
/// the full rows. Same span sequence as Gram-Schmidt (column j of Q spans the
/// first j + 1 inputs, positive R diagonal); rounding class. `Ok(None)` (every
/// rank alike) asks for the split Gram-Schmidt: a rank failed to fill, data
/// were nonfinite, or a pivot fell below `BLOCK_PIVOT_GUARD`. Collective.
fn orthonormalize_block(
    comm: Option<&crate::parallel::SharedComm>,
    fill_rows: &FillRows<'_>,
    m: usize,
    n: usize,
    gamma: bool,
) -> Result<Option<dft_numerics::linalg::Orthonormalized>, ScfError> {
    let profile = InitProfile::current();
    let fill_span =
        profile.map(|p| p.span(IC::InitAtomicBlockFillCalls, IC::InitAtomicBlockFillNanos));
    let rank = comm.map_or(0, Collective::rank);
    let counts = dft_mpi::block_counts(m, comm.map_or(1, Collective::size));
    let first: usize = counts[..rank].iter().sum();
    let rows = counts[rank];
    let mut a = vec![Complex64::default(); rows * n];
    let mut filled = true;
    if rows > 0 {
        for (j, column) in a.chunks_exact_mut(rows).enumerate() {
            if fill_rows(j, first, column).is_err() {
                filled = false;
                break;
            }
        }
    }
    let a = ComplexMatrix::new(rows, n, a).map_err(numerical)?;
    drop(fill_span);
    let all_finite = |x: &ComplexMatrix| x.as_slice().iter().all(finite);
    let failed = !filled || !all_finite(&a);
    let mut algebra = BlockAlgebra::new(gamma, rows, n)?;
    let Some(g) = block_gram(comm, &mut algebra, (rows > 0).then_some(&a), failed, n)? else {
        return Ok(None);
    };
    // First pass on the column-normalized Gram: Q1 = A D^-1 R1^-1.
    let factor_span = profile.map(|p| {
        p.span(
            IC::InitAtomicBlockFactorCalls,
            IC::InitAtomicBlockFactorNanos,
        )
    });
    let scale: Vec<f64> = (0..n).map(|j| g[j * n + j].re.sqrt()).collect();
    if scale.iter().any(|s| !s.is_finite() || *s == 0.0) {
        return Ok(None);
    }
    let mut gs = g;
    for j in 0..n {
        for i in 0..n {
            gs[j * n + i] /= scale[i] * scale[j];
        }
    }
    let Some(r1) = cholesky_upper(&gs, n, BLOCK_PIVOT_GUARD) else {
        return Ok(None);
    };
    let mut t1 = upper_inverse(&r1, n);
    for j in 0..n {
        for i in 0..n {
            t1[j * n + i] /= scale[i];
        }
    }
    drop(factor_span);
    let (q1, failed) = if rows > 0 {
        match algebra.update(&a, &t1) {
            Ok(q) => {
                let bad = !all_finite(&q);
                (q, bad)
            }
            Err(_) => (a, true),
        }
    } else {
        (a, false)
    };
    let Some(g2) = block_gram(comm, &mut algebra, (rows > 0).then_some(&q1), failed, n)? else {
        return Ok(None);
    };
    // Second pass: R2 is the identity to rounding for an accepted first pass.
    let factor_span = profile.map(|p| {
        p.span(
            IC::InitAtomicBlockFactorCalls,
            IC::InitAtomicBlockFactorNanos,
        )
    });
    let Some(r2) = cholesky_upper(&g2, n, 0.5) else {
        return Ok(None);
    };
    let t2 = upper_inverse(&r2, n);
    drop(factor_span);
    // Q1 is finite on every rank (checked by the reduction) and R2^-1 is
    // finite (pivots >= 0.5), so this GEMM cannot reject its inputs.
    let q = if rows > 0 {
        algebra.update(&q1, &t2)?
    } else {
        q1
    };
    let vectors = if let Some(comm) = comm {
        let _gather = profile.map(|p| {
            p.span(
                IC::InitAtomicBlockGatherCalls,
                IC::InitAtomicBlockGatherNanos,
            )
        });
        let local: Vec<f64> = q.as_slice().iter().flat_map(|z| [z.re, z.im]).collect();
        let vectors = gather_rows(comm, &counts, &local, n)?;
        crate::profile::distributed::count(crate::profile::distributed::Count::AtomicBlockStarts);
        vectors
    } else {
        q
    };
    Ok(Some(dft_numerics::linalg::Orthonormalized {
        vectors,
        dropped: Vec::new(),
    }))
}

fn finite(z: &Complex64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

/// `BAANDI_M2_ATOMIC=split`: Q^H H Q with the H applications and the Gram
/// columns split over ranks in contiguous column blocks (one allgatherv of the
/// small matrix). Each rank's block product is a GEMM of a column subset, so
/// the Gram matrix is in the rounding class of the replicated product
/// (bit-identical when the GEMM backend's per-entry order does not depend on
/// the column count). Collective.
fn split_gram(
    comm: &crate::parallel::SharedComm,
    q: &ComplexMatrix,
    h: &Hamiltonian,
) -> Result<ComplexMatrix, ScfError> {
    let (rows, n) = (q.rows(), q.cols());
    let counts = dft_mpi::block_counts(n, comm.size());
    let displs = dft_mpi::block_displacements(&counts);
    let (first, width) = (displs[comm.rank()], counts[comm.rank()]);
    let mut local = Vec::with_capacity(2 * n * width);
    if width > 0 {
        let block = ComplexMatrix::new(
            rows,
            width,
            q.as_slice()[first * rows..(first + width) * rows].to_vec(),
        )
        .map_err(numerical)?;
        let mut hq = ComplexMatrix::zeros(rows, width).map_err(numerical)?;
        h.apply_block(&block, &mut hq)?;
        let part = gemm(q, MatrixOp::Adjoint, &hq, MatrixOp::Identity).map_err(numerical)?;
        local.extend(part.as_slice().iter().flat_map(|z| [z.re, z.im]));
    }
    let gathered: Vec<usize> = counts.iter().map(|c| 2 * n * c).collect();
    let mut all = vec![0.0; 2 * n * n];
    comm.allgatherv_f64(&local, &gathered, &mut all)
        .map_err(|e| ScfError::Numerical(e.to_string()))?;
    ComplexMatrix::new(
        n,
        n,
        all.as_chunks::<2>()
            .0
            .iter()
            .map(|b| Complex64::new(b[0], b[1]))
            .collect(),
    )
    .map_err(numerical)
}

/// `BAANDI_M2_ATOMIC_BLOCK=on`: Q C with every rank multiplying its contiguous
/// block of rows (a GEMM of a row subset: rounding class of the replicated
/// product) and one allgatherv of the rows. Collective.
fn rotate_rows(
    comm: &crate::parallel::SharedComm,
    q: &ComplexMatrix,
    c: &ComplexMatrix,
) -> Result<ComplexMatrix, ScfError> {
    let (m, k, n) = (q.rows(), q.cols(), c.cols());
    let counts = dft_mpi::block_counts(m, comm.size());
    let first: usize = counts[..comm.rank()].iter().sum();
    let rows = counts[comm.rank()];
    let mut local = Vec::with_capacity(2 * rows * n);
    if rows > 0 {
        let mut block = Vec::with_capacity(rows * k);
        for j in 0..k {
            block.extend_from_slice(&q.column(j)[first..first + rows]);
        }
        let block = ComplexMatrix::new(rows, k, block).map_err(numerical)?;
        // Q and C are finite (validated by the Gram and the eigensolver), so
        // every rank either succeeds here or fails identically on its inputs.
        let part = gemm(&block, MatrixOp::Identity, c, MatrixOp::Identity).map_err(numerical)?;
        local.extend(part.as_slice().iter().flat_map(|z| [z.re, z.im]));
    }
    gather_rows(comm, &counts, &local, n)
}

/// Full column-major m x n block from every rank's contiguous row block
/// (`local`: this rank's rows x n block as interleaved re, im; `counts`: rows
/// per rank, `dft_mpi::block_counts`). One allgatherv. Collective.
fn gather_rows(
    comm: &crate::parallel::SharedComm,
    counts: &[usize],
    local: &[f64],
    n: usize,
) -> Result<ComplexMatrix, ScfError> {
    let m: usize = counts.iter().sum();
    let gathered: Vec<usize> = counts.iter().map(|r| 2 * r * n).collect();
    let mut all = vec![0.0; 2 * m * n];
    comm.allgatherv_f64(local, &gathered, &mut all)
        .map_err(|e| ScfError::Numerical(e.to_string()))?;
    let mut data = vec![Complex64::default(); m * n];
    let (mut at, mut d) = (0, 0);
    for &r in counts {
        for j in 0..n {
            for (out, b) in data[j * m + d..j * m + d + r]
                .iter_mut()
                .zip(all[at..at + 2 * r].as_chunks::<2>().0)
            {
                *out = Complex64::new(b[0], b[1]);
            }
            at += 2 * r;
        }
        d += r;
    }
    ComplexMatrix::new(m, n, data).map_err(numerical)
}

impl AtomicSpan {
    /// Orthonormal span columns (ATOMIC-REAL applies H to them on the real path).
    pub(super) fn vectors(&self) -> &ComplexMatrix {
        &self.q
    }
}

/// `hq`: H Q already computed (ATOMIC-REAL); `None` applies the complex H.
/// `split`: `BAANDI_M2_ATOMIC=split` on the M2 path (H Q and Gram columns
/// split over ranks, `split_gram`); ignored when `hq` is supplied.
pub(super) fn ritz(
    span: &AtomicSpan,
    h: &Hamiltonian,
    bands: usize,
    gamma: bool,
    hq: Option<ComplexMatrix>,
    split: Option<&crate::parallel::SharedComm>,
) -> Result<(ComplexMatrix, AtomicStartReport), ScfError> {
    use crate::profile::{AtomicPart, atomic_part, clock};
    let part = clock();
    let q = &span.q;
    let hq = hq.filter(|m| m.rows() == q.rows() && m.cols() == q.cols());
    let mut small = match (hq, split) {
        (None, Some(comm)) => split_gram(comm, q, h)?,
        (hq, _) => {
            let hq = match hq {
                Some(hq) => hq,
                None => {
                    let mut hq = ComplexMatrix::zeros(q.rows(), q.cols()).map_err(numerical)?;
                    h.apply_block(q, &mut hq)?;
                    hq
                }
            };
            gemm(q, MatrixOp::Adjoint, &hq, MatrixOp::Identity).map_err(numerical)?
        }
    };
    atomic_part(AtomicPart::Gram, part);
    let part = clock();
    for j in 0..small.cols() {
        for i in 0..=j {
            let v = (small[(i, j)] + small[(j, i)].conj()) * 0.5;
            small[(i, j)] = v;
            small[(j, i)] = v.conj();
        }
    }
    if gamma {
        for v in small.as_mut_slice() {
            if v.im.abs() >= 1e-12 {
                return Err(ScfError::InvalidInput(
                    "atomic Gamma reduced Hamiltonian is not real",
                ));
            }
            v.im = 0.0;
        }
    }
    let (_, rot) = hermitian_eigh(&small).map_err(numerical)?;
    atomic_part(AtomicPart::Eigh, part);
    let part = clock();
    let mut selected = ComplexMatrix::zeros(q.cols(), bands).map_err(numerical)?;
    for j in 0..bands {
        selected.column_mut(j).copy_from_slice(rot.column(j));
    }
    // BAANDI_M2_ATOMIC_BLOCK=on: the rotation by row blocks (one allgatherv).
    let guess = match split.filter(|_| crate::parallel::tuning::tuning().atomic_block) {
        Some(comm) => rotate_rows(comm, q, &selected)?,
        None => gemm(q, MatrixOp::Identity, &selected, MatrixOp::Identity).map_err(numerical)?,
    };
    atomic_part(AtomicPart::Rotate, part);
    let mut report = span.report.clone();
    report.initial_h_applications = q.cols();
    report.gram_rotation_macs = q.rows() * q.cols() * (q.cols() + bands);
    report.evd_orders.push(q.cols());
    Ok((guess, report))
}

#[cfg(test)]
mod tests {
    use super::*;
    use dft_mpi::{SimRank, simulate};
    use std::sync::Arc;

    /// `orthonormalize_split` reproduces `orthonormalize` bit for bit on every
    /// rank (ranks 2 to 5, cyclic ownership), including a zero column and an
    /// exactly dependent column that both routines drop.
    #[test]
    fn split_two_pass_gram_schmidt_is_bit_identical_to_the_serial_routine() {
        let (m, n) = (37, 13);
        let mut s: u64 = 17;
        let mut next = || {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((s >> 11) as f64 / (1u64 << 53) as f64) - 0.5
        };
        let mut block = ComplexMatrix::zeros(m, n).unwrap();
        for j in 0..n {
            for z in block.column_mut(j) {
                *z = Complex64::new(next(), next());
            }
        }
        block.column_mut(4).fill(Complex64::default());
        let (a, b) = (block.column(1).to_vec(), block.column(2).to_vec());
        for ((z, x), y) in block.column_mut(9).iter_mut().zip(&a).zip(&b) {
            *z = x * 0.5 - y * 2.0;
        }
        let want = orthonormalize(&block, None).unwrap();
        assert!(want.dropped.contains(&4) && want.dropped.contains(&9));
        let bits = |v: &ComplexMatrix| {
            v.as_slice()
                .iter()
                .flat_map(|z| [z.re.to_bits(), z.im.to_bits()])
                .collect::<Vec<_>>()
        };
        for ranks in 2..=5 {
            let got = simulate(ranks, |comm: &SimRank| {
                crate::parallel::set_communicator(Some(Arc::new(comm.clone())));
                let shared = crate::parallel::SharedComm::active().unwrap();
                let mine: Vec<usize> = (comm.rank()..n).step_by(ranks).collect();
                let own: Vec<Complex64> = mine
                    .iter()
                    .flat_map(|&j| block.column(j).to_vec())
                    .collect();
                let r = orthonormalize_split(&shared, own, &mine, m, n).unwrap();
                crate::parallel::set_communicator(None);
                r
            });
            for r in &got {
                assert_eq!(r.dropped, want.dropped, "P={ranks}");
                assert_eq!(bits(&r.vectors), bits(&want.vectors), "P={ranks}");
            }
        }
    }

    /// `orthonormalize_block` (BAANDI_M2_ATOMIC_BLOCK): for 2 to 7 ranks
    /// (including ranks without rows) every rank returns the same bits, the
    /// columns are orthonormal and agree with the serial two-pass Gram-Schmidt
    /// to rounding (same span sequence and phase convention); a dependent or
    /// zero column and a failing fill make every rank fall back (None).
    #[test]
    fn block_cholesky_qr_matches_gram_schmidt_and_falls_back_together() {
        let (m, n) = (41, 11);
        let mut s: u64 = 29;
        let mut next = || {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((s >> 11) as f64 / (1u64 << 53) as f64) - 0.5
        };
        let mut block = ComplexMatrix::zeros(m, n).unwrap();
        for j in 0..n {
            for z in block.column_mut(j) {
                *z = Complex64::new(next(), next());
            }
        }
        // Columns of very different length (the block path normalizes first).
        for z in block.column_mut(3) {
            *z *= 1e3;
        }
        // A Gamma-like block: row 0 real, rows 21..41 the conjugates of rows
        // 1..21 (Hermitian-symmetric columns: real overlaps, real branch).
        let mut mirrored = block.clone();
        for j in 0..n {
            let c = mirrored.column_mut(j);
            c[0].im = 0.0;
            for i in 1..=20 {
                c[20 + i] = c[i].conj();
            }
        }
        let bits = |v: &ComplexMatrix| {
            v.as_slice()
                .iter()
                .flat_map(|z| [z.re.to_bits(), z.im.to_bits()])
                .collect::<Vec<_>>()
        };
        let run = |block: &ComplexMatrix, ranks: usize, fail_rank: Option<usize>, gamma: bool| {
            if ranks == 1 {
                let fill = |j: usize, first: usize, column: &mut [Complex64]| {
                    if fail_rank == Some(0) {
                        return Err(ScfError::Numerical("fill".into()));
                    }
                    column.copy_from_slice(&block.column(j)[first..first + column.len()]);
                    Ok(())
                };
                return vec![orthonormalize_block(None, &fill, m, n, gamma).unwrap()];
            }
            simulate(ranks, |comm: &SimRank| {
                crate::parallel::set_communicator(Some(Arc::new(comm.clone())));
                let shared = crate::parallel::SharedComm::active().unwrap();
                let rank = comm.rank();
                let fill = |j: usize, first: usize, column: &mut [Complex64]| {
                    if fail_rank == Some(rank) {
                        return Err(ScfError::Numerical("fill".into()));
                    }
                    column.copy_from_slice(&block.column(j)[first..first + column.len()]);
                    Ok(())
                };
                let r = orthonormalize_block(Some(&shared), &fill, m, n, gamma).unwrap();
                crate::parallel::set_communicator(None);
                r
            })
        };
        for (block, gamma) in [(&block, false), (&mirrored, true)] {
            let want = orthonormalize(block, None).unwrap();
            assert!(want.dropped.is_empty());
            for ranks in 1..=7 {
                let got = run(block, ranks, None, gamma);
                let first = got[0].as_ref().expect("block path accepted");
                for r in &got {
                    let r = r.as_ref().unwrap();
                    assert!(r.dropped.is_empty());
                    assert_eq!(bits(&r.vectors), bits(&first.vectors), "P={ranks}");
                }
                let q = &first.vectors;
                let mut worst_gs = 0.0f64;
                for (x, y) in q.as_slice().iter().zip(want.vectors.as_slice()) {
                    worst_gs = worst_gs.max((x - y).norm());
                }
                assert!(
                    worst_gs < 1e-12,
                    "P={ranks}: |Q_block - Q_gs| = {worst_gs:e}"
                );
                let g = gemm(q, MatrixOp::Adjoint, q, MatrixOp::Identity).unwrap();
                for j in 0..n {
                    for i in 0..n {
                        let e = if i == j { 1.0 } else { 0.0 };
                        assert!((g[(i, j)] - e).norm() < 1e-14, "P={ranks} ({i},{j})");
                    }
                }
            }
        }
        // Exactly dependent column, zero column, failing fill on one rank.
        let mut dependent = block.clone();
        let (a, b) = (block.column(1).to_vec(), block.column(2).to_vec());
        for ((z, x), y) in dependent.column_mut(9).iter_mut().zip(&a).zip(&b) {
            *z = x * 0.5 - y * 2.0;
        }
        let mut zero = block.clone();
        zero.column_mut(4).fill(Complex64::default());
        for ranks in [1, 2, 3, 5] {
            assert!(
                run(&dependent, ranks, None, false)
                    .iter()
                    .all(Option::is_none)
            );
            assert!(run(&zero, ranks, None, false).iter().all(Option::is_none));
            assert!(
                run(&block, ranks, Some(ranks - 1), false)
                    .iter()
                    .all(Option::is_none)
            );
        }
    }
}
