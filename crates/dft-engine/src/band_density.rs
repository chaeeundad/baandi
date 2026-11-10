//! Reusable intra-k density transforms with an ascending-band reduction.
use crate::scf::ScfError;
use dft_basis::FftGridSpec;
use dft_physics::data::*;
use dft_physics::density::{FusedOperator, FusedPlan, SplitPairDensity};
use std::sync::{Arc, OnceLock};

/// DENSITY-PAIR (opt-in `BAANDI_DENSITY_PAIR=on`, default off): on a real Gamma
/// solve (Hermitian-completed eigenvectors) two occupied bands share one inverse
/// transform of `c_a + i c_b`; each band keeps its own contribution grid (read by
/// the certificate screen), so only the rounding of the per-band values changes.
/// Pairs are fixed as (0,1), (2,3), ... and worker ranges are even-aligned, so the
/// result does not depend on the thread count. Rounding class, not bit-identical.
fn density_pair() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var("BAANDI_DENSITY_PAIR").as_deref() == Ok("on"))
}
/// DENSITY-REDUCE (opt-in `BAANDI_DENSITY_REDUCE=on`, default off): the band sum
/// runs over fixed grid blocks (pool-parallel) with every point still summed in
/// ascending band order, and contributions of zero-occupation bands (exact +0
/// grids) are neither cleared again nor added. Bit-identical to the default.
fn density_reduce() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var("BAANDI_DENSITY_REDUCE").as_deref() == Ok("on"))
}
/// Grid points per reduction block (fixed; independent of the thread budget).
const REDUCE_BLOCK: usize = 4096;

pub(crate) struct BandDensity {
    workers: Vec<Worker>,
    pool: Arc<WorkerPool>,
    density: Array3<f64>,
    threads: usize,
    pair: bool,
    reduce: bool,
}
struct Worker {
    first: usize,
    fft: FusedPlan,
    contributions: Vec<Array3<f64>>,
    /// DENSITY-REDUCE: contribution holds values written by a nonzero weight.
    live: Vec<bool>,
    /// DENSITY-PAIR scratch: combined coefficients (empty when pairing is off).
    phi: Vec<Complex64>,
}
impl BandDensity {
    pub(crate) fn new(
        grid: &FftGridSpec,
        bands: usize,
        threads: usize,
        pool: Arc<WorkerPool>,
    ) -> Result<Self, ScfError> {
        Self::with_switches(grid, bands, threads, pool, density_pair(), density_reduce())
    }
    pub(crate) fn with_switches(
        grid: &FftGridSpec,
        bands: usize,
        threads: usize,
        pool: Arc<WorkerPool>,
        pair: bool,
        reduce: bool,
    ) -> Result<Self, ScfError> {
        // Pair mode: ranges in whole pairs so a pair never spans two workers.
        let units = if pair { bands.div_ceil(2) } else { bands };
        let scale = if pair { 2 } else { 1 };
        let count = threads.max(1).min(units);
        let mut workers = Vec::with_capacity(count);
        for i in 0..count {
            let first = (scale * (units * i / count)).min(bands);
            let end = (scale * (units * (i + 1) / count)).min(bands);
            workers.push(Worker {
                first,
                fft: FusedPlan::new(grid.dims(), grid.g_to_fft())?,
                contributions: (first..end).map(|_| Array3::zeros(grid.dims())).collect(),
                live: vec![true; end - first],
                phi: if pair {
                    vec![Complex64::default(); grid.g_to_fft().len()]
                } else {
                    Vec::new()
                },
            });
        }
        Ok(Self {
            workers,
            pool,
            density: Array3::zeros(grid.dims()),
            threads,
            pair,
            reduce,
        })
    }
    pub(crate) fn compute(
        &mut self,
        coefficients: &[Vec<Complex64>],
        occupations: &[f64],
        grid: &FftGridSpec,
        volume: f64,
    ) -> Result<&Array3<f64>, ScfError> {
        self.compute_with(coefficients, occupations, grid, volume, false)
    }
    /// `hermitian`: every coefficient vector is exactly Hermitian-completed (the
    /// real Gamma path), so each band is a real function and DENSITY-PAIR applies.
    pub(crate) fn compute_with(
        &mut self,
        coefficients: &[Vec<Complex64>],
        occupations: &[f64],
        grid: &FftGridSpec,
        volume: f64,
        hermitian: bool,
    ) -> Result<&Array3<f64>, ScfError> {
        let mode = Mode {
            pair: self.pair && hermitian,
            reduce: self.reduce,
        };
        let workers = &mut self.workers;
        if workers.len() == 1 {
            transform(workers, coefficients, occupations, grid, volume, mode)?;
        } else {
            self.pool
                .install(|| transform(workers, coefficients, occupations, grid, volume, mode))?;
        }
        let density = &mut self.density;
        if self.reduce {
            let live: Vec<&[f64]> = self
                .workers
                .iter()
                .flat_map(|w| w.contributions.iter().zip(&w.live))
                .filter(|(_, live)| **live)
                .map(|(c, _)| c.as_slice().expect("standard layout"))
                .collect();
            let out = density.as_slice_mut().expect("standard layout");
            if self.threads > 1 {
                self.pool.install(|| reduce_blocks(out, 0, &live, true));
            } else {
                reduce_blocks(out, 0, &live, false);
            }
        } else {
            density.fill(0.0);
            for worker in &self.workers {
                for contribution in &worker.contributions {
                    density.zip_mut_with(contribution, |n, c| *n += c);
                }
            }
        }
        if density.iter().any(|n: &f64| !n.is_finite()) {
            return Err(ScfError::NonFinite);
        }
        Ok(density)
    }
    pub(crate) fn density(&self) -> &Array3<f64> {
        &self.density
    }
    /// out[b] = sum over grid points of (band b's occupation-weighted density
    /// from the latest `compute`) times `field`, in ascending band order. The
    /// caller applies the grid weight and divides by the occupation. Read-only.
    pub(crate) fn weighted_integrals(&self, field: &[f64], out: &mut [f64]) {
        let mut band = 0;
        for worker in &self.workers {
            for contribution in &worker.contributions {
                if band < out.len() {
                    out[band] = contribution
                        .as_slice()
                        .map_or(f64::NAN, |c| c.iter().zip(field).map(|(a, b)| a * b).sum());
                }
                band += 1;
            }
        }
    }
}
/// out = sum over `parts` in order, per point, from +0; fixed blocks of
/// REDUCE_BLOCK points split in halves (whole blocks) over the pool.
fn reduce_blocks(out: &mut [f64], offset: usize, parts: &[&[f64]], parallel: bool) {
    let blocks = out.len().div_ceil(REDUCE_BLOCK);
    if parallel && blocks > 1 {
        let split = (blocks / 2) * REDUCE_BLOCK;
        let (left, right) = out.split_at_mut(split);
        WorkerPool::join(
            || reduce_blocks(left, offset, parts, true),
            || reduce_blocks(right, offset + split, parts, true),
        );
        return;
    }
    for (b, chunk) in out.chunks_mut(REDUCE_BLOCK).enumerate() {
        let start = offset + b * REDUCE_BLOCK;
        let len = chunk.len();
        chunk.fill(0.0);
        for part in parts {
            for (n, c) in chunk.iter_mut().zip(&part[start..start + len]) {
                *n += c;
            }
        }
    }
}
#[derive(Clone, Copy)]
struct Mode {
    pair: bool,
    reduce: bool,
}
fn check(band: usize, c: &[Complex64], grid: &FftGridSpec) -> Result<(), ScfError> {
    if c.len() != grid.g_to_fft().len() || c.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
    {
        return Err(ScfError::Numerical(format!(
            "band index {band}: invalid coefficients"
        )));
    }
    Ok(())
}
fn transform(
    workers: &mut [Worker],
    coefficients: &[Vec<Complex64>],
    occupations: &[f64],
    grid: &FftGridSpec,
    volume: f64,
    mode: Mode,
) -> Result<(), ScfError> {
    if workers.len() == 1 {
        let worker = &mut workers[0];
        if !mode.pair && !mode.reduce {
            for (local, contribution) in worker.contributions.iter_mut().enumerate() {
                let band = worker.first + local;
                let c = &coefficients[band];
                check(band, c, grid)?;
                contribution.fill(0.0);
                worker.fft.accumulate_density(
                    c,
                    occupations[band],
                    volume,
                    contribution.as_slice_mut().unwrap(),
                )?;
            }
            return Ok(());
        }
        let count = worker.contributions.len();
        let mut local = 0;
        while local < count {
            let band = worker.first + local;
            let c = &coefficients[band];
            check(band, c, grid)?;
            let w = occupations[band];
            if mode.pair
                && band.is_multiple_of(2)
                && local + 1 < count
                && w != 0.0
                && occupations[band + 1] != 0.0
            {
                let d = &coefficients[band + 1];
                check(band + 1, d, grid)?;
                // phi = c + i d: Re(psi) is band `band`, Im(psi) band `band + 1`.
                for ((p, a), b) in worker.phi.iter_mut().zip(c).zip(d) {
                    *p = Complex64::new(a.re - b.im, a.im + b.re);
                }
                let (left, right) = worker.contributions.split_at_mut(local + 1);
                let (da, db) = (&mut left[local], &mut right[0]);
                da.fill(0.0);
                db.fill(0.0);
                worker.fft.accumulate_density_pair_split(
                    &worker.phi,
                    [w, occupations[band + 1]],
                    volume,
                    da.as_slice_mut().unwrap(),
                    db.as_slice_mut().unwrap(),
                )?;
                worker.live[local] = true;
                worker.live[local + 1] = true;
                local += 2;
                continue;
            }
            let contribution = &mut worker.contributions[local];
            if mode.reduce && w == 0.0 {
                // Exact zero grid; clear only once after a nonzero weight.
                if worker.live[local] {
                    contribution.fill(0.0);
                    worker.live[local] = false;
                }
                local += 1;
                continue;
            }
            contribution.fill(0.0);
            worker
                .fft
                .accumulate_density(c, w, volume, contribution.as_slice_mut().unwrap())?;
            worker.live[local] = true;
            local += 1;
        }
        return Ok(());
    }
    let middle = workers.len() / 2;
    let (left, right) = workers.split_at_mut(middle);
    let first_left = left[0].first;
    let first_right = right[0].first;
    let (left, right) = WorkerPool::join(
        || {
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                transform(left, coefficients, occupations, grid, volume, mode)
            }))
        },
        || {
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                transform(right, coefficients, occupations, grid, volume, mode)
            }))
        },
    );
    left.map_err(|_| {
        ScfError::Numerical(format!("band index {first_left}: density worker panicked"))
    })??;
    right.map_err(|_| {
        ScfError::Numerical(format!("band index {first_right}: density worker panicked"))
    })?
}
#[cfg(test)]
mod tests {
    use super::*;
    use dft_basis::PlaneWaveBasis;
    use dft_geometry::Cell;
    #[test]
    fn ordered_density_and_reused_storage() {
        let cell = Cell::from_user_rows_bohr([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]).unwrap();
        let basis = PlaneWaveBasis::new(&cell, Hartree(2.5)).unwrap();
        let grid = FftGridSpec::new(&basis, &basis.density_basis(&cell).unwrap(), None).unwrap();
        let mut coefficients = Vec::new();
        for band in 0..7 {
            let mut c: Vec<_> = (0..basis.vectors().len())
                .map(|i| {
                    Complex64::new(((i + band * 3) as f64).sin(), ((2 * i + band) as f64).cos())
                })
                .collect();
            let norm = c.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
            for z in &mut c {
                *z /= norm;
            }
            coefficients.push(c);
        }
        let occupations = [2., 2., 2., 0.5, 0., 0., 0.];
        let expected = dft_physics::density::from_orbitals(
            &coefficients,
            &occupations,
            &grid,
            cell.volume_bohr3(),
            &mut FftPlan3d::new(grid.dims()).unwrap(),
        )
        .unwrap();
        let fused = dft_physics::density::from_orbitals_fused(
            &coefficients,
            &occupations,
            &grid,
            cell.volume_bohr3(),
            &mut FusedPlan::new(grid.dims(), grid.g_to_fft()).unwrap(),
        )
        .unwrap();
        assert_eq!(expected, fused);
        for threads in [1, 2, 3, 4, 8, 16] {
            let mut pool =
                BandDensity::new(&grid, 7, threads, WorkerPool::new(threads).unwrap()).unwrap();
            let addresses = |p: &BandDensity| {
                p.workers
                    .iter()
                    .map(|w| {
                        (
                            &w.fft as *const FusedPlan as usize,
                            w.contributions
                                .iter()
                                .map(|c| c.as_ptr() as usize)
                                .collect::<Vec<_>>(),
                        )
                    })
                    .collect::<Vec<_>>()
            };
            let before = addresses(&pool);
            let density_ptr = pool.density.as_ptr();
            for _ in 0..3 {
                let actual = pool
                    .compute(&coefficients, &occupations, &grid, cell.volume_bohr3())
                    .unwrap();
                assert_eq!(
                    actual.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                    expected.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
                );
                assert_eq!(before, addresses(&pool));
                assert_eq!(pool.density.as_ptr(), density_ptr);
            }
            let mut bad = coefficients.clone();
            bad[2][0].re = f64::NAN;
            bad[5][0].re = f64::NAN;
            assert!(
                pool.compute(&bad, &occupations, &grid, cell.volume_bohr3())
                    .unwrap_err()
                    .to_string()
                    .contains("band index 2:")
            );
            assert_eq!(before, addresses(&pool));
            assert_eq!(
                pool.compute(&coefficients, &occupations, &grid, cell.volume_bohr3())
                    .unwrap(),
                &expected
            );
        }
    }
    fn bits(a: &Array3<f64>) -> Vec<u64> {
        a.iter().map(|x| x.to_bits()).collect()
    }
    /// DENSITY-PAIR (rounding class, thread-count independent) and
    /// DENSITY-REDUCE (bit-identical) against the default path.
    #[test]
    fn pair_and_reduce_switches() {
        let cell =
            Cell::from_user_rows_bohr([[4.1, 0.3, 0.], [0., 3.9, 0.2], [0.1, 0., 4.4]]).unwrap();
        let basis = PlaneWaveBasis::new(&cell, Hartree(3.0)).unwrap();
        let grid = FftGridSpec::new(&basis, &basis.density_basis(&cell).unwrap(), None).unwrap();
        let half = dft_basis::HalfSphere::new(&basis, &grid).unwrap();
        let mirror = half.mirror();
        let n = basis.vectors().len();
        let bands = 9;
        // Hermitian-completed vectors: c(-G) = conj(c(G)), c(0) real.
        let coefficients: Vec<Vec<Complex64>> = (0..bands)
            .map(|band| {
                let mut c = vec![Complex64::default(); n];
                for i in 0..n {
                    let m = mirror[i] as usize;
                    if i == m {
                        c[i] = Complex64::new(((i + 5 * band) as f64).sin(), 0.);
                    } else if i < m {
                        c[i] = Complex64::new(
                            ((i + band * 3) as f64 * 0.7).sin(),
                            ((2 * i + band) as f64 * 1.3).cos(),
                        );
                        c[m] = c[i].conj();
                    }
                }
                let norm = c.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
                c.iter().map(|z| z / norm).collect()
            })
            .collect();
        let volume = cell.volume_bohr3();
        let sets: [&[f64]; 3] = [
            &[2., 2., 2., 2., 2., 2., 2., 2., 0.],
            &[2., 0., 2., 2., 0.5, 2., 0., 1., 2.],
            &[0., 2., 2., 2., 2., 0., 0., 0., 0.],
        ];
        for occupations in sets {
            let mut plain = BandDensity::with_switches(
                &grid,
                bands,
                1,
                WorkerPool::new(1).unwrap(),
                false,
                false,
            )
            .unwrap();
            let expected = plain
                .compute_with(&coefficients, occupations, &grid, volume, true)
                .unwrap()
                .clone();
            let mut expected_integrals = vec![0.; bands];
            let field: Vec<f64> = (0..expected.len()).map(|i| 1. + (i % 13) as f64).collect();
            plain.weighted_integrals(&field, &mut expected_integrals);
            let scale = expected.iter().fold(0_f64, |a, b| a.max(b.abs()));
            let mut pair_bits = None;
            for threads in [1, 2, 3, 4, 8] {
                for (pair, reduce) in [(false, true), (true, false), (true, true)] {
                    let mut pool = BandDensity::with_switches(
                        &grid,
                        bands,
                        threads,
                        WorkerPool::new(threads).unwrap(),
                        pair,
                        reduce,
                    )
                    .unwrap();
                    // Repeated calls with changing weights exercise the zero-grid reuse.
                    for _ in 0..2 {
                        pool.compute_with(&coefficients, &[1.; 9], &grid, volume, true)
                            .unwrap();
                        let actual = pool
                            .compute_with(&coefficients, occupations, &grid, volume, true)
                            .unwrap()
                            .clone();
                        if pair {
                            let worst = actual
                                .iter()
                                .zip(&expected)
                                .fold(0_f64, |a, (x, y)| a.max((x - y).abs()));
                            assert!(worst <= 1e-14 * scale, "pair {worst:e} vs {scale:e}");
                            let first = pair_bits.get_or_insert_with(|| bits(&actual));
                            assert_eq!(*first, bits(&actual), "pair threads {threads}");
                            let mut integrals = vec![0.; bands];
                            pool.weighted_integrals(&field, &mut integrals);
                            for (a, b) in integrals.iter().zip(&expected_integrals) {
                                assert!((a - b).abs() <= 1e-13 * b.abs().max(1.));
                            }
                        } else {
                            assert_eq!(bits(&expected), bits(&actual), "reduce {threads}");
                            let mut integrals = vec![0.; bands];
                            pool.weighted_integrals(&field, &mut integrals);
                            assert_eq!(integrals, expected_integrals);
                        }
                        // Non-Hermitian callers never pair: identical to the default.
                        let complex = pool
                            .compute_with(&coefficients, occupations, &grid, volume, false)
                            .unwrap()
                            .clone();
                        assert_eq!(bits(&expected), bits(&complex));
                    }
                }
            }
        }
    }
}
