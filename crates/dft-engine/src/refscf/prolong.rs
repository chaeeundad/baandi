//! Exact Miller-index injection and Fourier grid prolongation.
use crate::model::ScfError;
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_numerics::Fft3d;
use dft_physics::data::{Array3, Complex64, ComplexMatrix, FftPlan3d};
use std::collections::HashMap;

/// Copy coefficients by Miller index into the fine basis (zero for new G).
/// Fails if a coarse G is missing from the fine basis (bases not nested).
pub fn prolong_orbitals(
    coarse: &PlaneWaveBasis,
    fine: &PlaneWaveBasis,
    orbitals: &ComplexMatrix,
) -> Result<ComplexMatrix, ScfError> {
    let _span = OrbitalCopySpan::enter();
    if orbitals.rows() != coarse.vectors().len() {
        return Err(ScfError::InvalidInput(
            "cutoff_ladder: orbital rows differ from the coarse basis",
        ));
    }
    let map = prolong_indices(coarse, fine)?;
    let rows = fine.vectors().len();
    let mut out = ComplexMatrix::zeros(rows, orbitals.cols())
        .map_err(|e| ScfError::Numerical(format!("cutoff_ladder: {e}")))?;
    for j in 0..orbitals.cols() {
        let dst = out.column_mut(j);
        for (src, &i) in orbitals.column(j).iter().zip(&map) {
            dst[i] = *src;
        }
    }
    Ok(out)
}

/// Shared exact Miller injection; no coefficient arithmetic.
pub(crate) fn prolong_indices(
    coarse: &PlaneWaveBasis,
    fine: &PlaneWaveBasis,
) -> Result<Vec<usize>, ScfError> {
    let index: HashMap<[i32; 3], usize> = fine
        .vectors()
        .iter()
        .enumerate()
        .map(|(i, g)| (g.miller, i))
        .collect();
    coarse
        .vectors()
        .iter()
        .map(|g| index.get(&g.miller).copied())
        .collect::<Option<Vec<_>>>()
        .ok_or(ScfError::InvalidInput(
            "cutoff_ladder: coarse basis is not a subset of the fine basis",
        ))
}

/// Fourier interpolation of a real periodic field from the coarse FFT grid to a
/// finer one (zero padding). The Nyquist plane of an even coarse dimension is
/// dropped: a field band-limited to the coarse density sphere has no weight
/// there because the grid holds that sphere (dims >= 2 max|m| + 1).
pub fn prolong_density(
    coarse: &Array3<f64>,
    fine_dims: [usize; 3],
) -> Result<Array3<f64>, ScfError> {
    let cd = [coarse.dim().0, coarse.dim().1, coarse.dim().2];
    let mut fft = FftPlan3d::new(cd)?;
    let input = coarse.mapv(|v| Complex64::new(v, 0.0));
    let mut reciprocal = Array3::<Complex64>::zeros(cd);
    fft.to_reciprocal(&input, &mut reciprocal)?;
    let mut padded = Array3::<Complex64>::zeros(fine_dims);
    for ((i, j, k), z) in reciprocal.indexed_iter() {
        let (Some(a), Some(b), Some(c)) = (
            prolong_frequency(i, cd[0]),
            prolong_frequency(j, cd[1]),
            prolong_frequency(k, cd[2]),
        ) else {
            continue;
        };
        let [p, q, r] = FftGridSpec::map_miller([a, b, c], fine_dims).map_err(|_| {
            ScfError::InvalidInput("cutoff_ladder: coarse frequencies exceed the fine grid")
        })?;
        padded[[p, q, r]] = *z;
    }
    let mut real = Array3::<Complex64>::zeros(fine_dims);
    let mut fine_fft = FftPlan3d::new(fine_dims)?;
    fine_fft.to_real(&padded, &mut real)?;
    Ok(real.mapv(|z| z.re))
}

/// Shared Fourier-padding convention; an even-grid Nyquist plane is absent
/// from the representable density sphere. Used by both ladder drivers.
pub(crate) fn prolong_frequency(i: usize, n: usize) -> Option<i32> {
    if n.is_multiple_of(2) && i == n / 2 {
        None
    } else if i <= (n - 1) / 2 {
        Some(i as i32)
    } else {
        Some(i as i32 - n as i32)
    }
}

// Keep the reference ladder's observational copy timing in the results layer.
static ORBITAL_COPY_NS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
struct OrbitalCopySpan(std::time::Instant);
impl OrbitalCopySpan {
    fn enter() -> Option<Self> {
        crate::profile::enabled().then(|| Self(std::time::Instant::now()))
    }
}
impl Drop for OrbitalCopySpan {
    fn drop(&mut self) {
        ORBITAL_COPY_NS.fetch_add(
            self.0.elapsed().as_nanos() as u64,
            std::sync::atomic::Ordering::Relaxed,
        );
    }
}
#[cfg(feature = "strict-oracle")]
pub(crate) fn orbital_copy_seconds() -> f64 {
    ORBITAL_COPY_NS.load(std::sync::atomic::Ordering::Relaxed) as f64 * 1e-9
}
