use super::*;
use dft_numerics::{Fft3d, FftPlan3d};
use dft_units::{Angstrom, Bohr};
use std::collections::BTreeSet;

fn cubic(a: f64) -> Cell {
    Cell::from_user_rows_bohr([[a, 0., 0.], [0., a, 0.], [0., 0., a]]).unwrap()
}
fn silicon() -> Cell {
    let h = Bohr::from(Angstrom(5.43)).0 / 2.;
    Cell::from_user_rows_bohr([[0., h, h], [h, 0., h], [h, h, 0.]]).unwrap()
}
fn indices(b: &PlaneWaveBasis) -> BTreeSet<[i32; 3]> {
    b.vectors().iter().map(|g| g.miller).collect()
}

#[test]
fn exact_shell_and_tolerance() {
    let cell = cubic(TAU); // Reciprocal columns are unit vectors in Bohr^-1.
    let b = PlaneWaveBasis::new(&cell, Hartree(0.5)).unwrap();
    assert_eq!(b.vectors().len(), 7);
    assert_eq!(b.vectors()[0].miller, [0; 3]);
    assert_eq!(b.vectors()[1].miller, [-1, 0, 0]);
    assert_eq!(
        PlaneWaveBasis::new(&cell, Hartree(0.5 * (1. - 0.5e-12)))
            .unwrap()
            .vectors()
            .len(),
        7
    );
    assert_eq!(
        PlaneWaveBasis::new(&cell, Hartree(0.5 * (1. - 2e-12)))
            .unwrap()
            .vectors()
            .len(),
        1
    );
}

#[test]
fn complete_spheres_and_ordering() {
    let cells = [
        cubic(10.),
        silicon(),
        Cell::from_user_rows_bohr([[4., 0., 0.], [-2., 12_f64.sqrt(), 0.], [0., 0., 6.]]).unwrap(),
        Cell::from_user_rows_bohr([[3., 0., 0.], [29., 0.7, 0.], [41., 2.3, 0.2]]).unwrap(),
    ];
    for cell in cells {
        let basis = PlaneWaveBasis::new(&cell, Hartree(2.)).unwrap();
        for b in [&basis, &basis.density_basis(&cell).unwrap()] {
            let set = indices(b);
            assert_eq!(set.len(), b.vectors().len());
            assert_eq!(b.vectors()[0].miller, [0; 3]);
            for g in b.vectors() {
                assert!(0.5 * g.norm_squared <= b.cutoff().0 * (1. + CUTOFF_TOLERANCE));
                assert!(set.contains(&g.miller.map(|m| -m)));
            }
            for pair in b.vectors().windows(2) {
                assert!(
                    pair[0].norm_squared < pair[1].norm_squared
                        || (pair[0].norm_squared == pair[1].norm_squared
                            && pair[0].miller < pair[1].miller)
                );
            }
            // Independent brute force over a box two integers larger on every side.
            let bound = bounding_box(&cell, b.cutoff()).unwrap().map(|x| x + 2);
            let reciprocal = cell.reciprocal_columns_inverse_bohr();
            let mut brute = BTreeSet::new();
            for i in -bound[0]..=bound[0] {
                for j in -bound[1]..=bound[1] {
                    for k in -bound[2]..=bound[2] {
                        let mut norm = 0.;
                        for row in reciprocal {
                            let component = row[0] * f64::from(i)
                                + row[1] * f64::from(j)
                                + row[2] * f64::from(k);
                            norm += component * component;
                        }
                        if norm / 2. <= b.cutoff().0 * (1. + CUTOFF_TOLERANCE) {
                            brute.insert([i, j, k]);
                        }
                    }
                }
            }
            assert_eq!(set, brute);
        }
        let density = indices(&basis.density_basis(&cell).unwrap());
        for g in basis.vectors() {
            for h in basis.vectors() {
                assert!(density.contains(&std::array::from_fn(|i| g.miller[i] - h.miller[i])));
            }
        }
    }
}

#[test]
fn continuum_count() {
    let cell = cubic(TAU);
    let b = PlaneWaveBasis::new(&cell, Hartree(200.)).unwrap(); // Radius 20 in Miller space.
    let estimate = cell.volume_bohr3() * 400_f64.powf(1.5) / (6. * std::f64::consts::PI.powi(2));
    let relative = (b.vectors().len() as f64 / estimate - 1.).abs();
    println!(
        "continuum: count={} estimate={estimate:.12} relative_error={relative:.16e}",
        b.vectors().len()
    );
    // Finite lattice-point count vs sphere volume: surface-shell oscillations
    // vanish relative to volume as radius grows; 5% is a coarse asymptotic check.
    assert!(relative < 0.05);
}

#[test]
fn grids_mapping_and_errors() {
    let cell = silicon();
    let basis = PlaneWaveBasis::new(&cell, Hartree(15.)).unwrap();
    let density = basis.density_basis(&cell).unwrap();
    let grid = FftGridSpec::new(&basis, &density, None).unwrap();
    println!(
        "Si a=5.43 Angstrom Ecut=30 Ry: wavefunctions={} density={} FFT={:?}",
        basis.vectors().len(),
        density.vectors().len(),
        grid.dims()
    );
    for i in 0..3 {
        assert!(grid.dims()[i] >= grid.required_dims()[i]);
        assert!(smooth(grid.dims()[i]));
    }
    let mapped: BTreeSet<_> = density
        .vectors()
        .iter()
        .map(|g| FftGridSpec::map_miller(g.miller, grid.dims()).unwrap())
        .collect();
    assert_eq!(mapped.len(), density.vectors().len());
    for (g, &[i, j, k]) in grid.g_to_fft().iter().enumerate() {
        assert_eq!(
            grid.fft_to_g()[(i * grid.dims()[1] + j) * grid.dims()[2] + k],
            Some(g)
        );
    }
    for axis in 0..3 {
        let mut small = grid.required_dims();
        small[axis] -= 1;
        assert!(matches!(
            FftGridSpec::new(&basis, &density, Some(small)),
            Err(BasisError::GridTooSmall { .. })
        ));
    }
    assert!(FftGridSpec::new(&basis, &density, Some([usize::MAX; 3])).is_err());
    assert!(grid.scatter(&[]).is_err());
    assert!(grid.gather(&Array3::zeros([1; 3])).is_err());
    assert_eq!(
        FftGridSpec::map_miller([-2, -1, 2], [4, 5, 6]).unwrap(),
        [2, 4, 2]
    );
    assert_eq!(
        FftGridSpec::map_miller([1, 2, -3], [4, 5, 6]).unwrap(),
        [1, 2, 3]
    );
    for m in [
        [2, 0, 0],
        [-3, 0, 0],
        [0, 3, 0],
        [0, 0, 3],
        [i32::MIN, 0, 0],
    ] {
        assert!(matches!(
            FftGridSpec::map_miller(m, [4, 5, 6]),
            Err(BasisError::IndexOutOfRange { .. })
        ));
    }
    assert!(FftGridSpec::map_miller([0; 3], [0; 3]).is_err());
    for cutoff in [0., -1., f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(PlaneWaveBasis::new(&cell, Hartree(cutoff)).is_err());
    }
}

#[test]
fn fft_roundtrip_normalization_and_enlargement() {
    let cell = silicon();
    let basis = PlaneWaveBasis::new(&cell, Hartree(3.)).unwrap();
    let density = basis.density_basis(&cell).unwrap();
    let small = FftGridSpec::new(&basis, &density, None).unwrap();
    // Fixed-seed LCG test fixture; no extra random-number dependency.
    let mut seed = 0x1234_5678_9abc_def0_u64;
    let mut random = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        (seed >> 11) as f64 / ((1_u64 << 53) as f64) - 0.5
    };
    let mut coefficients: Vec<_> = basis
        .vectors()
        .iter()
        .map(|_| Complex64::new(random(), random()))
        .collect();
    let norm = coefficients
        .iter()
        .map(|z| z.norm_sqr())
        .sum::<f64>()
        .sqrt();
    for c in &mut coefficients {
        *c /= norm;
    }
    let large = FftGridSpec::new(&basis, &density, Some(small.dims().map(|n| n + 3))).unwrap();
    assert!(large.overridden()); // Non-smooth overrides are accepted too.
    let mut results = Vec::new();
    for grid in [small, large] {
        let reciprocal = grid.scatter(&coefficients).unwrap();
        assert_eq!(grid.gather(&reciprocal).unwrap(), coefficients);
        assert_eq!(
            reciprocal.iter().filter(|z| z.norm() != 0.).count(),
            coefficients.len()
        );
        let mut real = Array3::zeros(grid.dims());
        let mut back = real.clone();
        let mut plan = FftPlan3d::new(grid.dims()).unwrap();
        plan.to_real(&reciprocal, &mut real).unwrap();
        plan.to_reciprocal(&real, &mut back).unwrap();
        let gathered = grid.gather(&back).unwrap();
        let error = gathered
            .iter()
            .zip(&coefficients)
            .map(|(a, b)| (a - b).norm())
            .fold(0., f64::max);
        let psi_norm = real
            .iter()
            .map(|z| (*z / cell.volume_bohr3().sqrt()).norm_sqr())
            .sum::<f64>()
            * cell.volume_bohr3()
            / real.len() as f64;
        println!(
            "FFT {:?}: coefficient_max_abs={error:.16e} normalization_abs={:.16e}",
            grid.dims(),
            (psi_norm - 1.).abs()
        );
        assert!(error < 1e-12); // Dimensionless absolute coefficient error.
        assert!((psi_norm - 1.).abs() < 1e-12); // Dimensionless integrated probability.
        results.push(gathered);
    }
    let difference = results[0]
        .iter()
        .zip(&results[1])
        .map(|(a, b)| (a - b).norm())
        .fold(0., f64::max);
    println!("grid enlargement coefficient_max_abs={difference:.16e}");
    assert!(difference < 1e-12);
}

#[test]
fn shifted_basis_keeps_density_grid_and_handles_skewed_gauge() {
    let cell = Cell::from_user_rows_bohr([
        [TAU, 0.0, 0.0],
        [100.0 * TAU, 0.1 * TAU, 0.0],
        [0.0, 0.0, TAU],
    ])
    .unwrap();
    let cutoff = Hartree(0.001);
    let gamma = PlaneWaveBasis::new(&cell, cutoff).unwrap();
    let shifted = PlaneWaveBasis::at_k(&cell, cutoff, [0.01, 0.0, 0.0]).unwrap();
    let density = shifted.density_basis(&cell).unwrap();
    assert_eq!(
        density.vectors(),
        gamma.density_basis(&cell).unwrap().vectors()
    );
    // q=(0.01,0,0), G=(0,10,0): its periodic part aliases to the sole grid
    // entry, but its density and all represented differences are exactly G=0.
    assert_eq!(shifted.vectors().len(), 1);
    assert_eq!(shifted.vectors()[0].miller, [0, 1, 0]);
    let grid = FftGridSpec::new(&shifted, &density, None).unwrap();
    assert_eq!(
        grid.dims(),
        FftGridSpec::new(&gamma, &density, None).unwrap().dims()
    );
    let c = [Complex64::new(0.6, 0.8)];
    assert_eq!(grid.gather(&grid.scatter(&c).unwrap()).unwrap(), c);
    let error = (shifted.momenta()[0].norm_squared / 2.0 - 0.00005).abs();
    println!(
        "skewed shifted gauge: grid={:?}, kinetic error={error:.16e} Ha",
        grid.dims()
    );
    assert!(error < 1e-12);
    assert!(PlaneWaveBasis::at_k(&cell, cutoff, [f64::NAN, 0.0, 0.0]).is_err());
}

#[test]
fn frozen_millers_preserve_order_under_strain_and_shell_crossing() {
    let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [0., 8., 0.], [0., 0., 8.]]).unwrap();
    let deformed =
        Cell::from_user_rows_bohr([[9., 0.1, 0.], [0., 7., 0.2], [0., 0., 12.]]).unwrap();
    for k in [[0.; 3], [0.1, 0.2, -0.3]] {
        let base = PlaneWaveBasis::at_k(&cell, Hartree(2.), k).unwrap();
        let frozen = base.in_cell(&deformed).unwrap();
        assert_eq!(
            base.vectors().iter().map(|g| g.miller).collect::<Vec<_>>(),
            frozen
                .vectors()
                .iter()
                .map(|g| g.miller)
                .collect::<Vec<_>>()
        );
        let rebuilt = PlaneWaveBasis::at_k(&deformed, Hartree(2.), k).unwrap();
        assert_ne!(rebuilt.vectors().len(), frozen.vectors().len());
        for (g, q) in frozen.vectors().iter().zip(frozen.momenta()) {
            for a in 0..3 {
                assert!((g.cartesian[a] + frozen.k_cartesian()[a] - q.cartesian[a]).abs() < 1e-14);
            }
        }
        let original = base.in_cell(&cell).unwrap();
        assert_eq!(base.vectors(), original.vectors());
    }
}

#[test]
fn frozen_density_support_survives_expansion_shell_crossing() {
    let reference = Cell::from_user_rows_bohr([[5., 0., 0.], [0., 5., 0.], [0., 0., 5.]]).unwrap();
    let expanded = Cell::from_user_rows_bohr([[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]]).unwrap();
    let basis = PlaneWaveBasis::at_k(&reference, Hartree(2.), [0.25; 3]).unwrap();
    let density = basis.density_basis(&reference).unwrap();
    let grid = FftGridSpec::new(&basis, &density, None).unwrap();
    let plain = basis.in_cell(&expanded).unwrap();
    assert!(
        FftGridSpec::new(
            &plain,
            &plain.density_basis(&expanded).unwrap(),
            Some(grid.dims())
        )
        .is_err()
    );
    let frozen = basis
        .in_cell_with_frozen_density(&expanded, &reference)
        .unwrap();
    let moved = frozen.density_basis(&expanded).unwrap();
    assert_eq!(
        density
            .vectors()
            .iter()
            .map(|g| g.miller)
            .collect::<Vec<_>>(),
        moved.vectors().iter().map(|g| g.miller).collect::<Vec<_>>()
    );
    assert!(FftGridSpec::new(&frozen, &moved, Some(grid.dims())).is_ok());
    assert!(
        moved
            .vectors()
            .iter()
            .zip(density.vectors())
            .any(|(a, b)| a.norm_squared != b.norm_squared)
    );
}
