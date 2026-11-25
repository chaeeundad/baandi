//! Analytic tests of the surviving reference prolongation helpers.
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_engine::refscf::prolong as cutoff_ladder;
use dft_geometry::Cell;
use dft_numerics::linalg::ComplexMatrix;
use dft_physics::data::{Array3, Complex64};
use dft_units::Hartree;

#[test]
fn orbital_prolongation_is_an_exact_injection() {
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0.5, 4., 0.], [0.2, 0.4, 6.]]).unwrap();
    let coarse = PlaneWaveBasis::new(&cell, Hartree(2.0)).unwrap();
    let fine = PlaneWaveBasis::new(&cell, Hartree(4.0)).unwrap();
    let n = coarse.vectors().len();
    let mut c = ComplexMatrix::zeros(n, 3).unwrap();
    for j in 0..3 {
        for (i, z) in c.column_mut(j).iter_mut().enumerate() {
            *z = Complex64::new(((i * 7 + j * 3) as f64).sin(), ((i + 5 * j) as f64).cos());
        }
    }
    let p = cutoff_ladder::prolong_orbitals(&coarse, &fine, &c).unwrap();
    assert_eq!(p.rows(), fine.vectors().len());
    // Every coarse coefficient sits at its Miller index; every new G is zero;
    // all inner products are unchanged bit for bit (same terms, zeros added).
    let mut copied = 0;
    for (i, g) in fine.vectors().iter().enumerate() {
        match coarse.vectors().iter().position(|h| h.miller == g.miller) {
            Some(k) => {
                copied += 1;
                for j in 0..3 {
                    assert_eq!(p.column(j)[i], c.column(j)[k]);
                }
            }
            None => (0..3).for_each(|j| assert_eq!(p.column(j)[i], Complex64::default())),
        }
    }
    assert_eq!(copied, n);
    // Bases are not nested the other way round.
    assert!(cutoff_ladder::prolong_orbitals(&fine, &coarse, &p).is_err());
}

#[test]
fn density_prolongation_is_exact_fourier_interpolation() {
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0.5, 4., 0.], [0.2, 0.4, 6.]]).unwrap();
    let wave = PlaneWaveBasis::new(&cell, Hartree(2.0)).unwrap();
    let density = wave.density_basis(&cell).unwrap();
    let coarse = FftGridSpec::new(&density, &density, None).unwrap();
    let fine_wave = PlaneWaveBasis::new(&cell, Hartree(4.0)).unwrap();
    let fine_density = fine_wave.density_basis(&cell).unwrap();
    let fine = FftGridSpec::new(&fine_density, &fine_density, None).unwrap();
    // A real band-limited field on the coarse density sphere, evaluated directly
    // at the sample points of both grids.
    let terms: Vec<_> = density
        .vectors()
        .iter()
        .enumerate()
        .map(|(i, g)| {
            (
                g.miller,
                1.0 / (1.0 + g.norm_squared) * (1.0 + 0.3 * (i as f64).sin()),
            )
        })
        .collect();
    let field = |dims: [usize; 3]| {
        Array3::from_shape_fn(dims, |(a, b, c)| {
            let r = [
                a as f64 / dims[0] as f64,
                b as f64 / dims[1] as f64,
                c as f64 / dims[2] as f64,
            ];
            // Symmetrised so that the field is real: sum over +G and -G.
            terms
                .iter()
                .map(|(m, w)| {
                    let phase = std::f64::consts::TAU
                        * (m[0] as f64 * r[0] + m[1] as f64 * r[1] + m[2] as f64 * r[2]);
                    w * phase.cos()
                })
                .sum::<f64>()
                + 3.0
        })
    };
    let on_coarse = field(coarse.dims());
    let expected = field(fine.dims());
    let prolonged = cutoff_ladder::prolong_density(&on_coarse, fine.dims()).unwrap();
    let scale = expected.iter().fold(0f64, |m, v| m.max(v.abs()));
    let error = prolonged
        .iter()
        .zip(&expected)
        .fold(0f64, |m, (a, b)| m.max((a - b).abs()));
    assert!(error < 1e-12 * scale, "interpolation error {error:e}");
    let mean = |a: &Array3<f64>| a.sum() / a.len() as f64;
    assert!((mean(&prolonged) - mean(&on_coarse)).abs() < 1e-13 * scale);
}
