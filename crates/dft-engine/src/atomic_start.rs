//! P45A atomic Fourier span and first-density Rayleigh-Ritz initialization.
use super::*;
use dft_basis::HalfSphere;
use dft_mpi::Collective;
use dft_numerics::linalg::{MatrixOp, gemm, hermitian_eigh, orthonormalize};
use dft_solvers::LinearOperator;
use dft_upf::radial::{projector_transform_with_rule, real_spherical_harmonic};

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum StartWavefunctions {
    FreeElectron,
    /// Default (ADR 0042): atomic-orbital subspace with Rayleigh-Ritz; falls back to
    /// the free-electron start with a recorded reason when it does not apply.
    #[default]
    Atomic,
}
impl StartWavefunctions {
    pub fn name(self) -> &'static str {
        match self {
            Self::FreeElectron => "free_electron",
            Self::Atomic => "atomic",
        }
    }
}

/// Initialization work is separate from Davidson work; MACs are complex MACs.
#[derive(Debug, Clone, Default)]
pub struct AtomicStartReport {
    pub atomic_orbitals: usize,
    pub initial_h_applications: usize,
    pub gram_rotation_macs: usize,
    pub mgs_projection_macs: usize,
    pub radial_samples: usize,
    pub evd_orders: Vec<usize>,
    pub fallback_reason: Option<String>,
}

fn numerical(e: impl std::fmt::Display) -> ScfError {
    ScfError::Numerical(e.to_string())
}

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
    let npw = basis.vectors().len();
    let count = species
        .iter()
        .map(|s| s.positions_bohr.len() * s.orbitals.iter().map(|c| 2 * c.l + 1).sum::<usize>())
        .sum();
    let mut report = AtomicStartReport {
        atomic_orbitals: count,
        ..Default::default()
    };
    // Radial transforms per orbital (cached by |G|^2), then one descriptor per
    // column in the established column order: species, orbital, position, m.
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
    let fill = |j: usize, column: &mut [Complex64]| -> Result<(), ScfError> {
        let (r, l, m, pos, il) = columns[j];
        let radial = &radials[r];
        for (i, g) in basis.momenta().iter().enumerate() {
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
            column[i] = il * phase * (radial[i] * y);
        }
        Ok(())
    };
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
        let ortho = orthonormalize_split(comm, own, &mine, npw, count)?;
        return finish_span(basis, grid, bands, ortho, npw, report);
    }
    let mut raw = ComplexMatrix::zeros(npw, count).map_err(numerical)?;
    let ortho = match pool {
        Some(pool) if npw > 0 => {
            pool.install(|| {
                dft_numerics::pool::for_each_chunk_mut(raw.as_mut_slice(), npw, &|first, column| {
                    fill(first / npw, column)
                })
            })?;
            dft_numerics::linalg::orthonormalize_pooled(&raw, None, pool).map_err(numerical)?
        }
        _ => {
            for j in 0..count {
                fill(j, raw.column_mut(j))?;
            }
            orthonormalize(&raw, None).map_err(numerical)?
        }
    };
    finish_span(basis, grid, bands, ortho, npw, report)
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
    let mut selected = ComplexMatrix::zeros(q.cols(), bands).map_err(numerical)?;
    for j in 0..bands {
        selected.column_mut(j).copy_from_slice(rot.column(j));
    }
    let guess = gemm(q, MatrixOp::Identity, &selected, MatrixOp::Identity).map_err(numerical)?;
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
}
