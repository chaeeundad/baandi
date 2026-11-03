//! Reusable intra-k density transforms with an ascending-band reduction.
use crate::scf::ScfError;
use dft_basis::FftGridSpec;
use dft_physics::data::*;
use dft_physics::density::{FusedOperator, FusedPlan};

pub(crate) struct BandDensity {
    workers: Vec<Worker>,
}
struct Worker {
    first: usize,
    fft: FusedPlan,
    contributions: Vec<Array3<f64>>,
}
impl BandDensity {
    pub(crate) fn new(grid: &FftGridSpec, bands: usize, threads: usize) -> Result<Self, ScfError> {
        let count = threads.max(1).min(bands);
        let mut workers = Vec::with_capacity(count);
        for i in 0..count {
            let first = bands * i / count;
            let end = bands * (i + 1) / count;
            workers.push(Worker {
                first,
                fft: FusedPlan::new(grid.dims(), grid.g_to_fft())?,
                contributions: (first..end).map(|_| Array3::zeros(grid.dims())).collect(),
            });
        }
        Ok(Self { workers })
    }
    pub(crate) fn compute(
        &mut self,
        coefficients: &[Vec<Complex64>],
        occupations: &[f64],
        grid: &FftGridSpec,
        volume: f64,
    ) -> Result<Array3<f64>, ScfError> {
        transform(&mut self.workers, coefficients, occupations, grid, volume)?;
        let mut density = Array3::zeros(grid.dims());
        for worker in &self.workers {
            for contribution in &worker.contributions {
                density.zip_mut_with(contribution, |n, c| *n += c);
            }
        }
        if density.iter().any(|n: &f64| !n.is_finite()) {
            return Err(ScfError::NonFinite);
        }
        Ok(density)
    }
}
fn transform(
    workers: &mut [Worker],
    coefficients: &[Vec<Complex64>],
    occupations: &[f64],
    grid: &FftGridSpec,
    volume: f64,
) -> Result<(), ScfError> {
    if workers.len() == 1 {
        let worker = &mut workers[0];
        for (local, contribution) in worker.contributions.iter_mut().enumerate() {
            let band = worker.first + local;
            let c = &coefficients[band];
            if c.len() != grid.g_to_fft().len()
                || c.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
            {
                return Err(ScfError::Numerical(format!(
                    "band index {band}: invalid coefficients"
                )));
            }
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
    let middle = workers.len() / 2;
    let (left, right) = workers.split_at_mut(middle);
    let first_left = left[0].first;
    let first_right = right[0].first;
    std::thread::scope(|scope| {
        let handle = scope.spawn(|| transform(left, coefficients, occupations, grid, volume));
        let right = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            transform(right, coefficients, occupations, grid, volume)
        }));
        let left = handle.join();
        left.map_err(|_| {
            ScfError::Numerical(format!("band index {first_left}: density worker panicked"))
        })??;
        right.map_err(|_| {
            ScfError::Numerical(format!("band index {first_right}: density worker panicked"))
        })?
    })
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
        for threads in [1, 2, 3, 4, 8] {
            let mut pool = BandDensity::new(&grid, 7, threads).unwrap();
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
            for _ in 0..3 {
                let actual = pool
                    .compute(&coefficients, &occupations, &grid, cell.volume_bohr3())
                    .unwrap();
                assert_eq!(
                    actual.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                    expected.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
                );
                assert_eq!(before, addresses(&pool));
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
                expected
            );
        }
    }
}
