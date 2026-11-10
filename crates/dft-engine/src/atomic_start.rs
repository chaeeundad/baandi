//! P45A atomic Fourier span and first-density Rayleigh-Ritz initialization.
use super::*;
use dft_basis::HalfSphere;
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

pub(super) fn span(
    cell: &Cell,
    basis: &PlaneWaveBasis,
    grid: &FftGridSpec,
    species: &[AtomicSpecies],
    options: &ScfOptions,
    bands: usize,
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
    let mut raw = ComplexMatrix::zeros(npw, count).map_err(numerical)?;
    let mut col = 0;
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
                    for (i, g) in basis.momenta().iter().enumerate() {
                        let y = if g.norm_squared == 0.0 {
                            if c.l == 0 {
                                1.0 / (4.0 * std::f64::consts::PI).sqrt()
                            } else {
                                0.0
                            }
                        } else {
                            real_spherical_harmonic(c.l, m, g.cartesian).map_err(numerical)?
                        };
                        let phase = hamiltonian::structure_factor(g.cartesian, &[pos])?;
                        raw[(i, col)] = il * phase * (radial[i] * y);
                    }
                    col += 1;
                }
            }
        }
    }
    let ortho = orthonormalize(&raw, None).map_err(numerical)?;
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

impl AtomicSpan {
    /// Orthonormal span columns (ATOMIC-REAL applies H to them on the real path).
    pub(super) fn vectors(&self) -> &ComplexMatrix {
        &self.q
    }
}

/// `hq`: H Q already computed (ATOMIC-REAL); `None` applies the complex H.
pub(super) fn ritz(
    span: &AtomicSpan,
    h: &Hamiltonian,
    bands: usize,
    gamma: bool,
    hq: Option<ComplexMatrix>,
) -> Result<(ComplexMatrix, AtomicStartReport), ScfError> {
    let q = &span.q;
    let hq = match hq.filter(|m| m.rows() == q.rows() && m.cols() == q.cols()) {
        Some(hq) => hq,
        None => {
            let mut hq = ComplexMatrix::zeros(q.rows(), q.cols()).map_err(numerical)?;
            h.apply_block(q, &mut hq)?;
            hq
        }
    };
    let mut small = gemm(q, MatrixOp::Adjoint, &hq, MatrixOp::Identity).map_err(numerical)?;
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
