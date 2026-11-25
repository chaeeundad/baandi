//! Starting density from the superposition of the UPF atomic charges (the
//! reference code's default starting potential), as an opt-in alternative to the
//! uniform density. Unpolarized fixed-occupation problems only.
use super::*;
pub use crate::model::StartDensity;

fn numerical(e: impl std::fmt::Display) -> ScfError {
    ScfError::Numerical(e.to_string())
}

/// rho(G) = (1/Omega) sum_atoms exp(-i G.tau) int 4 pi r^2 rho_at(r) sinc(G r) dr,
/// on the density sphere, transformed to the grid, clipped at zero and rescaled
/// to the electron count. None when a species carries no atomic charge samples.
///
/// INIT-DENSITY-SPLIT: with a registered communicator of two or more ranks
/// every rank evaluates the coefficients of its contiguous block of G vectors
/// (`superposition_rows`, the per-G arithmetic of the full evaluation: the
/// same radial integral per |G| and the same species/atom accumulation order)
/// and one allgatherv assembles them, bit-identical for any rank count. A
/// failure on any rank makes every rank repeat the replicated evaluation
/// (same error). Collective with a communicator.
pub(super) fn atomic_superposition(
    cell: &Cell,
    density_basis: &PlaneWaveBasis,
    dg: &FftGridSpec,
    species: &[Species<'_>],
    fft: &mut FftPlan3d,
    electrons: f64,
) -> Result<Option<Array3<f64>>, ScfError> {
    let volume = cell.volume_bohr3();
    let vectors = density_basis.vectors();
    // Replicated data check (same on every rank, no collective needed).
    for s in species {
        let rho = &s.potential.raw().rhoatom.values;
        if rho.len() != s.potential.mesh().r().len() || rho.iter().all(|v| *v == 0.0) {
            return Ok(None);
        }
    }
    let coefficients = split_superposition(vectors, volume, species)?;
    let reciprocal = dg.scatter(&coefficients).map_err(numerical)?;
    let mut real = Array3::zeros(reciprocal.raw_dim());
    fft.to_real(&reciprocal, &mut real)?;
    let mut density = real.mapv(|z: Complex64| z.re.max(0.0));
    let charge = density.sum() * volume / density.len() as f64;
    if !(charge.is_finite() && charge > 0.0) {
        return Ok(None);
    }
    density.mapv_inplace(|v| v * electrons / charge);
    Ok(Some(density))
}

/// Coefficients of `vectors` (all of them without a communicator, else by
/// contiguous blocks and one allgatherv; see `atomic_superposition`).
fn split_superposition(
    vectors: &[dft_basis::GVector],
    volume: f64,
    species: &[Species<'_>],
) -> Result<Vec<Complex64>, ScfError> {
    use dft_mpi::Collective;
    let Some(comm) = crate::parallel::SharedComm::active() else {
        return superposition_rows(vectors, volume, species);
    };
    let n = vectors.len();
    let counts = dft_mpi::block_counts(n, comm.size());
    let first: usize = counts[..comm.rank()].iter().sum();
    let local = superposition_rows(
        &vectors[first..first + counts[comm.rank()]],
        volume,
        species,
    );
    // Per rank: 2 values per coefficient, then a failure flag.
    let mut send = vec![0.0; 2 * counts[comm.rank()] + 1];
    match &local {
        Ok(values) => {
            for (b, z) in send.as_chunks_mut::<2>().0.iter_mut().zip(values) {
                b[0] = z.re;
                b[1] = z.im;
            }
        }
        Err(_) => *send.last_mut().unwrap() = 1.0,
    }
    let gathered: Vec<usize> = counts.iter().map(|c| 2 * c + 1).collect();
    let mut all = vec![0.0; gathered.iter().sum()];
    comm.allgatherv_f64(&send, &gathered, &mut all)
        .map_err(|e| ScfError::Numerical(e.to_string()))?;
    let mut result = Vec::with_capacity(n);
    let mut at = 0;
    for &c in &counts {
        if all[at + 2 * c] != 0.0 {
            return superposition_rows(vectors, volume, species);
        }
        result.extend(
            all[at..at + 2 * c]
                .as_chunks::<2>()
                .0
                .iter()
                .map(|b| Complex64::new(b[0], b[1])),
        );
        at += 2 * c + 1;
    }
    crate::m2::count(crate::m2::Count::DensitySplits);
    Ok(result)
}

/// Superposition coefficients of `vectors` (any contiguous part of the
/// density sphere): per G, the radial integral of its |G| and the sum over
/// species and atoms in input order. The caller has checked that every
/// species carries atomic charge samples.
fn superposition_rows(
    vectors: &[dft_basis::GVector],
    volume: f64,
    species: &[Species<'_>],
) -> Result<Vec<Complex64>, ScfError> {
    let mut coefficients = vec![Complex64::default(); vectors.len()];
    for s in species {
        let rho = &s.potential.raw().rhoatom.values;
        let mesh = s.potential.mesh();
        let mut cache = std::collections::HashMap::new();
        let mut radial = Vec::with_capacity(vectors.len());
        for g in vectors {
            let key = g.norm_squared.to_bits();
            let value = if let Some(v) = cache.get(&key) {
                *v
            } else {
                let q = g.norm_squared.sqrt();
                let samples: Vec<f64> = mesh
                    .r()
                    .iter()
                    .zip(rho)
                    .map(|(&r, &f)| {
                        if q == 0.0 {
                            f
                        } else {
                            f * dft_upf::radial::spherical_bessel(0, q * r).unwrap_or(0.0)
                        }
                    })
                    .collect();
                let v = mesh.integrate(&samples).map_err(numerical)? / volume;
                cache.insert(key, v);
                v
            };
            radial.push(value);
        }
        for pos in s.positions_bohr {
            for ((c, g), f) in coefficients.iter_mut().zip(vectors).zip(&radial) {
                let phase =
                    -(g.cartesian[0] * pos[0] + g.cartesian[1] * pos[1] + g.cartesian[2] * pos[2]);
                *c += Complex64::new(phase.cos(), phase.sin()) * f;
            }
        }
    }
    Ok(coefficients)
}
