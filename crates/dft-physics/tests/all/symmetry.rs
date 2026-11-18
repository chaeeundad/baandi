use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_geometry::{Atom, Cell, symmetry::detect};
use dft_numerics::FftPlan3d;
use dft_physics::symmetry::DensitySymmetrizer;
use dft_units::Hartree;
use ndarray::Array3;
use std::f64::consts::TAU;
#[test]
fn translated_density_independent_real_space_oracle() {
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0., 5., 0.], [0., 0., 5.]]).unwrap();
    let ops = detect(&cell, &[Atom::new("X", [0.13, 0.21, 0.37])], 1e-7).unwrap();
    let b = PlaneWaveBasis::new(&cell, Hartree(3.)).unwrap();
    let g = FftGridSpec::new(&b, &b, None).unwrap();
    let dims = g.dims();
    let field = |x: [f64; 3]| 2.0 + 0.1 * (TAU * (x[0] + x[1])).cos() + 0.2 * (TAU * x[2]).sin();
    let mut n = Array3::from_shape_fn(dims, |(i, j, k)| {
        field([
            i as f64 / dims[0] as f64,
            j as f64 / dims[1] as f64,
            k as f64 / dims[2] as f64,
        ])
    });
    // Cubic rotations are orthogonal signed permutations: independently evaluate
    // n(R^T(x-tau)) directly, without the production reciprocal mapping or FFT.
    let oracle = Array3::from_shape_fn(dims, |(i, j, k)| {
        let x = [
            i as f64 / dims[0] as f64,
            j as f64 / dims[1] as f64,
            k as f64 / dims[2] as f64,
        ];
        ops.iter()
            .map(|o| {
                field(std::array::from_fn(|a| {
                    (0..3)
                        .map(|j| o.rotation[j][a] as f64 * (x[j] - o.translation[j]))
                        .sum()
                }))
            })
            .sum::<f64>()
            / ops.len() as f64
    });
    let mut plan = DensitySymmetrizer::new(&b, &g, &ops).unwrap();
    let mut fft = FftPlan3d::new(dims).unwrap();
    plan.apply(&mut n, &mut fft).unwrap();
    let error = n
        .iter()
        .zip(&oracle)
        .map(|(a, b)| (a - b).abs())
        .fold(0., f64::max);
    println!("translated density real-space reconstruction max error={error:e} Bohr^-3");
    assert!(error < 1e-13);
    let previous = n.clone();
    plan.apply(&mut n, &mut fft).unwrap();
    assert!(n.iter().zip(&previous).all(|(a, b)| (a - b).abs() < 1e-13));
    assert!((n.sum() / n.len() as f64 - 2.).abs() < 1e-13);
}
