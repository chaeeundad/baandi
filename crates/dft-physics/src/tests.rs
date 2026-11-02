use super::*;
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_geometry::Cell;
use dft_numerics::{Fft3d, FftPlan3d};
use dft_units::{Angstrom, Bohr, Hartree};
use ndarray::Array3;
use std::f64::consts::{PI, TAU};

fn si(scale: f64) -> Cell {
    let a = Bohr::from(Angstrom(5.43)).0 * scale;
    Cell::from_user_rows_bohr([
        [0., a / 2., a / 2.],
        [a / 2., 0., a / 2.],
        [a / 2., a / 2., 0.],
    ])
    .unwrap()
}
fn relative(a: f64, b: f64, tolerance: f64) -> f64 {
    let error = (a - b).abs() / b.abs();
    assert!(
        error < tolerance,
        "{a:.16e} vs {b:.16e}: relative {error:e} >= {tolerance:e}"
    );
    error
}

#[test]
fn density_and_kinetic() {
    let cell = si(1.0);
    let basis = PlaneWaveBasis::new(&cell, Hartree(5.0)).unwrap(); // 10 Ry
    let grid = FftGridSpec::new(&basis, &basis.density_basis(&cell).unwrap(), None).unwrap();
    let mut fft = FftPlan3d::new(grid.dims()).unwrap();
    let mut seed = 0x1234_5678_9abc_def0_u64;
    let mut random = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        (seed >> 11) as f64 / ((1_u64 << 53) as f64) - 0.5
    };
    let coefficients: Vec<Vec<Complex64>> = (0..4)
        .map(|_| {
            let mut c: Vec<_> = basis
                .vectors()
                .iter()
                .map(|_| Complex64::new(random(), random()))
                .collect();
            let norm = c.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
            c.iter_mut().for_each(|z| *z /= norm);
            c
        })
        .collect();
    let occupations = [2., 2., 1., 0.5];
    let n = density::from_orbitals(
        &coefficients,
        &occupations,
        &grid,
        cell.volume_bohr3(),
        &mut fft,
    )
    .unwrap();
    let electrons = n.sum() * cell.volume_bohr3() / n.len() as f64;
    let err = relative(electrons, 5.5, 1e-12);
    let mut wave = vec![Complex64::default(); basis.vectors().len()];
    wave[3] = Complex64::new(1., 0.);
    let t = kinetic::energy(&basis, &[wave.clone()], &[1.]).unwrap().0;
    let single_err = relative(t, 0.5 * basis.vectors()[3].norm_squared, 1e-14);
    wave[3] = Complex64::new(0.6, 0.);
    wave[20] = Complex64::new(0., 0.8);
    let mixture = kinetic::energy(&basis, &[wave], &[2.]).unwrap().0;
    let mix_err = relative(
        mixture,
        0.36 * basis.vectors()[3].norm_squared + 0.64 * basis.vectors()[20].norm_squared,
        1e-14,
    );
    println!(
        "density: basis={} grid={:?} electrons={electrons:.16e} sum_rel={err:.3e}; kinetic single={t:.16e} Ha rel={single_err:.3e}, mixture={mixture:.16e} Ha rel={mix_err:.3e}",
        basis.vectors().len(),
        grid.dims()
    );
}

#[test]
fn hartree_analytic() {
    let cell = si(1.0);
    let dims = [12, 10, 9];
    let mut fft = FftPlan3d::new(dims).unwrap();
    let uniform = Array3::from_elem(dims, 0.3);
    let zero = hartree::evaluate(&uniform, &cell, &mut fft).unwrap();
    let max_zero = zero.potential.iter().map(|z| z.norm()).fold(0., f64::max);
    assert!(zero.energy.0.abs() < 1e-14 && max_zero < 1e-14);
    let amplitude = 0.02;
    let cosine = Array3::from_shape_fn(dims, |(i, j, _)| {
        (TAU * (i as f64 / dims[0] as f64 - j as f64 / dims[1] as f64)).cos()
    });
    let n = cosine.mapv(|c| 0.3 + 2.0 * amplitude * c);
    let result = hartree::evaluate(&n, &cell, &mut fft).unwrap();
    let g2s = hartree::squared_frequencies(&cell, dims).unwrap();
    let g2 = g2s[[1, dims[1] - 1, 0]];
    let expected = 4.0 * PI * cell.volume_bohr3() * amplitude * amplitude / g2;
    let energy_error = relative(result.energy.0, expected, 1e-12);
    let vscale = 8.0 * PI * amplitude / g2;
    let potential_error = result
        .potential
        .iter()
        .zip(&cosine)
        .map(|(v, c)| (*v - Complex64::new(vscale * c, 0.)).norm() / vscale)
        .fold(0., f64::max);
    assert!(potential_error < 1e-12);
    // Exactly compare E_H to 1/2 integral (n - mean(n)) V_H, complex residual included.
    let mean = n.sum() / n.len() as f64;
    let integral: Complex64 = n
        .iter()
        .zip(&result.potential)
        .map(|(n, v)| (n - mean) * v)
        .sum::<Complex64>()
        * (0.5 * cell.volume_bohr3() / n.len() as f64);
    let integral_error = relative(integral.re, result.energy.0, 1e-12);
    assert!(integral.im.abs() / expected < 1e-12);
    let mut laplacian = Array3::zeros(dims);
    let spectral = &result.reciprocal_potential * &g2s;
    fft.to_real(&spectral, &mut laplacian).unwrap();
    let poisson_error = laplacian
        .iter()
        .zip(&n)
        .map(|(v, n)| {
            (*v - Complex64::new(4. * PI * (n - mean), 0.)).norm() / (8. * PI * amplitude)
        })
        .fold(0., f64::max);
    assert!(poisson_error < 1e-12);
    let second = Array3::from_shape_fn(dims, |(_, j, k)| {
        0.1 + 0.01 * (TAU * (j as f64 / 10. + k as f64 / 9.)).cos()
    });
    let h2 = hartree::evaluate(&second, &cell, &mut fft).unwrap();
    let combined = hartree::evaluate(&(&n * 0.7 + &second * 1.3), &cell, &mut fft).unwrap();
    let linear_error = combined
        .potential
        .iter()
        .zip(result.potential.iter().zip(&h2.potential))
        .map(|(v, (v1, v2))| (*v - *v1 * 0.7 - *v2 * 1.3).norm() / vscale)
        .fold(0., f64::max);
    assert!(linear_error < 1e-12);
    println!(
        "Hartree: uniform E={:.3e} Ha maxV={max_zero:.3e} Ha; cosine E={:.16e} Ha energy_rel={energy_error:.3e} potential_rel={potential_error:.3e} integral_rel={integral_error:.3e} Poisson_rel={poisson_error:.3e} linearity_rel={linear_error:.3e}",
        zero.energy.0, result.energy.0
    );
}

#[test]
fn xc_cross_checks() {
    use xcx::{Functional, FunctionalId, Spin, XcInput};
    let exchange = Functional::new(FunctionalId::LdaX, Spin::Unpolarized).unwrap();
    let pw = Functional::new(FunctionalId::LdaCPw, Spin::Unpolarized).unwrap();
    let densities: Vec<f64> = (0..=70)
        .map(|i| 10_f64.powf(-6. + i as f64 / 10.))
        .collect();
    let independent = exchange
        .eval(densities.len(), &XcInput::lda(&densities))
        .unwrap();
    let correlation = pw.eval(densities.len(), &XcInput::lda(&densities)).unwrap();
    let mut exchange_error = 0_f64;
    let mut pw_difference = 0_f64;
    for (i, &n) in densities.iter().enumerate() {
        let p = xc::lda_pz(n).unwrap();
        let closed = -3. / (4. * PI) * (3. * PI * PI * n).cbrt();
        relative(p.eps_x, closed, 1e-12);
        exchange_error = exchange_error.max(relative(p.eps_x, independent.exc[i], 1e-12));
        relative(p.v_x, independent.vrho[i], 1e-12);
        pw_difference = pw_difference.max((p.eps_c - correlation.exc[i]).abs());
    }
    // PW92 is a different correlation fit: a loose 5 mHa/electron sanity bound.
    assert!(pw_difference < 0.005);
    let mut max_derivative = 0_f64;
    let mut rs_values: Vec<f64> = (0..=100)
        .map(|i| 0.1 * (200_f64).powf(i as f64 / 100.))
        .collect();
    rs_values.extend([0.999, 1.001]); // Stencils stay on either side of rounded branch join.
    for rs in rs_values {
        let n = 3. / (4. * PI * rs.powi(3));
        let h = n * 1e-5;
        // Differentiate the actual energy density, never the analytic potential.
        let finite_difference = ((n + h) * xc::lda_pz(n + h).unwrap().eps_xc()
            - (n - h) * xc::lda_pz(n - h).unwrap().eps_xc())
            / (2. * h);
        max_derivative = max_derivative.max(relative(
            xc::lda_pz(n).unwrap().v_xc(),
            finite_difference,
            1e-7,
        ));
    }
    let left = xc::lda_pz(3. / (4. * PI * (1_f64 - 1e-10).powi(3))).unwrap();
    let right = xc::lda_pz(3. / (4. * PI * (1_f64 + 1e-10).powi(3))).unwrap();
    let expected_eps_jump = -0.1423 / (1. + 1.0529 + 0.3334) - (-0.048 - 0.0116);
    let expected_v_jump = -0.1423 / (1. + 1.0529 + 0.3334)
        - 0.1423 * (1.0529 / 2. + 0.3334) / (3. * (1_f64 + 1.0529 + 0.3334).powi(2))
        - (-0.048 - 0.0116 - (0.0311 + 0.0020 - 0.0116) / 3.);
    assert!((right.eps_c - left.eps_c - expected_eps_jump).abs() < 1e-10);
    assert!((right.v_c - left.v_c - expected_v_jump).abs() < 1e-10);
    assert!((right.eps_c - left.eps_c).abs() < 1e-4 && (right.v_c - left.v_c).abs() < 1e-4);
    let n = 0.1;
    let result = xc::evaluate(&Array3::from_elem([4, 5, 6], n), 270.).unwrap();
    relative(
        result.energy.0,
        270. * n * xc::lda_pz(n).unwrap().eps_xc(),
        1e-12,
    );
    relative(
        result.potential_integral.0,
        270. * n * xc::lda_pz(n).unwrap().v_xc(),
        1e-12,
    );
    println!(
        "XC: exchange xcx max_rel={exchange_error:.3e}; energy-density FD max_rel={max_derivative:.3e}; PZ81-PW92 max_abs={pw_difference:.16e} Ha/electron; rs=1 right-left eps_c={:.16e} Ha/electron (parameter limit={expected_eps_jump:.16e}), v_c={:.16e} Ha (parameter limit={expected_v_jump:.16e}); uniform E_xc={:.16e} Ha integral_nv={:.16e} Ha",
        right.eps_c - left.eps_c,
        right.v_c - left.v_c,
        result.energy.0,
        result.potential_integral.0
    );
}

#[test]
fn validation_and_thresholds() {
    for n in [-1e-14, -1e-16, 0., 1e-12, 1e-10] {
        let p = xc::lda_pz(n).unwrap();
        assert_eq!(p.eps_xc(), 0.);
        assert_eq!(p.v_xc(), 0.);
    }
    assert!(xc::lda_pz(1.00001e-10).unwrap().eps_xc() < 0.);
    for n in [
        -1.00001e-14,
        -1.,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
    ] {
        assert!(matches!(xc::lda_pz(n), Err(PhysicsError::Density(_))));
        assert!(xc::evaluate(&Array3::from_elem([1, 1, 1], n), 1.).is_err());
        assert!(
            hartree::evaluate(
                &Array3::from_elem([1, 1, 1], n),
                &si(1.),
                &mut FftPlan3d::new([1; 3]).unwrap()
            )
            .is_err()
        );
    }
    assert!(xc::evaluate(&Array3::zeros([0, 1, 1]), 1.).is_err());
    assert!(xc::evaluate(&Array3::zeros([1; 3]), 0.).is_err());
    let basis = PlaneWaveBasis::new(&si(1.), Hartree(0.5)).unwrap();
    assert!(kinetic::energy(&basis, &[vec![]], &[1.]).is_err());
    assert!(
        kinetic::energy(
            &basis,
            &[vec![Complex64::default(); basis.vectors().len()]],
            &[-1.]
        )
        .is_err()
    );
    let ion = ewald::Ion {
        charge: 4.,
        fractional_position: [0.; 3],
    };
    assert!(matches!(
        ewald::energy(&si(1.), &[ion, ion], Default::default()),
        Err(PhysicsError::CoincidentIons)
    ));
    assert!(
        ewald::energy(
            &si(1.),
            &[ion],
            ewald::EwaldOptions {
                eta: Some(0.),
                ..Default::default()
            }
        )
        .is_err()
    );
}

#[test]
fn ewald_invariances() {
    use ewald::{EwaldOptions, Ion, energy};
    let ions = [
        Ion {
            charge: 4.,
            fractional_position: [0.; 3],
        },
        Ion {
            charge: 4.,
            fractional_position: [0.25; 3],
        },
    ];
    let cell = si(1.);
    let reference = energy(&cell, &ions, Default::default()).unwrap();
    println!(
        "Ewald Si informational: E={:.16e} Ha = {:.16e} Ry; real={:.16e} reciprocal={:.16e} self={:.16e} background={:.16e} Ha; tail_bound={:.3e} Ha Rcut={:.8} Bohr Gcut={:.8} Bohr^-1",
        reference.energy.0,
        2. * reference.energy.0,
        reference.real.0,
        reference.reciprocal.0,
        reference.self_energy.0,
        reference.background.0,
        reference.tail_bound.0,
        reference.real_cutoff_bohr,
        reference.reciprocal_cutoff_inverse_bohr
    );
    for factor in [0.65, 0.85, 1.2, 1.6] {
        let result = energy(
            &cell,
            &ions,
            EwaldOptions {
                eta: Some(reference.eta * factor),
                ..Default::default()
            },
        )
        .unwrap();
        let err = relative(result.energy.0, reference.energy.0, 1e-10);
        println!(
            "Ewald eta={:.10} Bohr^-1 E={:.16e} Ha relative={err:.3e}",
            result.eta, result.energy.0
        );
    }
    let refined = energy(
        &cell,
        &ions,
        EwaldOptions {
            tolerance: 1e-16,
            ..Default::default()
        },
    )
    .unwrap();
    let refinement = (refined.energy.0 - reference.energy.0).abs();
    assert!(refinement < reference.tail_bound.0 + refined.tail_bound.0 + 1e-13);
    let translated = ions.map(|i| Ion {
        fractional_position: std::array::from_fn(|k| {
            i.fractional_position[k] + [-1.31, 0.47, 2.93][k]
        }),
        ..i
    });
    let translation = relative(
        energy(&cell, &translated, Default::default())
            .unwrap()
            .energy
            .0,
        reference.energy.0,
        1e-10,
    );
    for scale in [0.8, 1.3] {
        let result = energy(&si(scale), &ions, Default::default()).unwrap();
        let err = relative(result.energy.0, reference.energy.0 / scale, 1e-10);
        println!(
            "Ewald scale={scale} E={:.16e} Ha relative={err:.3e}",
            result.energy.0
        );
    }
    let mut columns = cell.lattice_columns_bohr();
    for row in &mut columns {
        row[0] *= 2.;
    }
    let supercell = Cell::from_matrix_columns_bohr(columns).unwrap();
    let replicated: Vec<_> = (0..2)
        .flat_map(|replica| {
            ions.map(|i| Ion {
                fractional_position: [
                    (i.fractional_position[0] + replica as f64) / 2.,
                    i.fractional_position[1],
                    i.fractional_position[2],
                ],
                ..i
            })
        })
        .collect();
    let result = energy(&supercell, &replicated, Default::default()).unwrap();
    let super_error = relative(result.energy.0, 2. * reference.energy.0, 1e-9);
    println!(
        "Ewald translation_rel={translation:.3e} refinement_abs={refinement:.3e} Ha supercell={:.16e} Ha relative={super_error:.3e}",
        result.energy.0
    );
}

#[test]
fn hartree_nyquist_and_odd_grids() {
    // Independently construct Cartesian G from the documented signed indices.
    // On a skew cell, mixed Nyquist conjugate partners have different G^2.
    let cell = si(1.);
    for dims in [[6, 5, 7], [7, 5, 7]] {
        let n = Array3::from_shape_fn(dims, |(i, j, _)| {
            0.2 + 0.02 * (TAU * (-3. * i as f64 / dims[0] as f64 + j as f64 / 5.)).cos()
        });
        let mut fft = FftPlan3d::new(dims).unwrap();
        let h = hartree::evaluate(&n, &cell, &mut fft).unwrap();
        let indices = if dims[0] == 6 {
            [[3, 1, 0], [3, 4, 0]]
        } else {
            [[4, 1, 0], [3, 4, 0]]
        };
        let millers = if dims[0] == 6 {
            [[-3., 1., 0.], [-3., -1., 0.]]
        } else {
            [[-3., 1., 0.], [3., -1., 0.]]
        };
        let b = cell.reciprocal_columns_inverse_bohr();
        let mut expected_energy = 0.;
        for (index, m) in indices.into_iter().zip(millers) {
            let norm: f64 = b
                .iter()
                .map(|row| row.iter().zip(m).map(|(b, m)| b * m).sum::<f64>().powi(2))
                .sum();
            let expected = 4. * PI * 0.01 / norm;
            relative(h.reciprocal_potential[index].re, expected, 1e-12);
            assert!(h.reciprocal_potential[index].im.abs() < 1e-14);
            expected_energy += 0.5 * cell.volume_bohr3() * 4. * PI * 0.01_f64.powi(2) / norm;
        }
        relative(h.energy.0, expected_energy, 1e-12);
        let g2 = hartree::squared_frequencies(&cell, dims).unwrap();
        let mut poisson = Array3::zeros(dims);
        fft.to_real(&(&h.reciprocal_potential * &g2), &mut poisson)
            .unwrap();
        let error = poisson
            .iter()
            .zip(&n)
            .map(|(p, n)| (*p - Complex64::new(4. * PI * (n - 0.2), 0.)).norm() / (4. * PI * 0.02))
            .fold(0., f64::max);
        assert!(error < 1e-12);
        println!(
            "Hartree grid={dims:?} Nyquist/signed-mode E={:.16e} Ha Poisson_rel={error:.3e}",
            h.energy.0
        );
    }
}

#[test]
fn ewald_lattice_representation_and_charge_scaling() {
    // The same lattice described by a unimodular skew/left-handed basis.
    // The physical energy must be invariant despite different shell boxes.
    let a = 7.;
    let cubic = Cell::from_user_rows_bohr([[a, 0., 0.], [0., a, 0.], [0., 0., a]]).unwrap();
    let skew = Cell::from_user_rows_bohr([[a, 0., 0.], [2. * a, a, 0.], [0., 0., -a]]).unwrap();
    let ions = [
        ewald::Ion {
            charge: 1.,
            fractional_position: [0.; 3],
        },
        ewald::Ion {
            charge: -1.,
            fractional_position: [0.31, 0.27, 0.43],
        },
    ];
    let remapped = ions.map(|i| ewald::Ion {
        fractional_position: [
            i.fractional_position[0] - 2. * i.fractional_position[1],
            i.fractional_position[1],
            -i.fractional_position[2],
        ],
        ..i
    });
    let original = ewald::energy(&cubic, &ions, Default::default()).unwrap();
    let transformed = ewald::energy(&skew, &remapped, Default::default()).unwrap();
    let error = relative(transformed.energy.0, original.energy.0, 1e-10);
    assert_eq!(original.background.0, 0.);
    let scaled = ions.map(|i| ewald::Ion {
        charge: 2. * i.charge,
        ..i
    });
    relative(
        ewald::energy(&cubic, &scaled, Default::default())
            .unwrap()
            .energy
            .0,
        4. * original.energy.0,
        1e-12,
    );
    println!(
        "Ewald neutral skew/left-handed E={:.16e} Ha representation_rel={error:.3e}",
        original.energy.0
    );
}
