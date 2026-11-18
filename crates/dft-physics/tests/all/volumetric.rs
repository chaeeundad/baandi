use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_geometry::Cell;
use dft_physics::{data::*, volumetric::*};
use std::f64::consts::PI;

#[test]
fn uniform_gas_and_vacuum() {
    let n = Array3::from_elem([3, 4, 5], 0.2_f64);
    let tau = n.mapv(|n| 0.3 * (3.0 * PI * PI).powf(2.0 / 3.0) * n.powf(5.0 / 3.0));
    let e = elf_from_fields(&n, &tau, &Array3::zeros(n.raw_dim())).unwrap();
    assert!(e.iter().all(|&e| e == 0.5));
    let zero = Array3::zeros([1, 1, 1]);
    assert_eq!(
        elf_from_fields(&zero, &zero, &zero).unwrap()[[0, 0, 0]],
        0.0
    );
    assert!(elf_from_fields(&n, &zero, &zero).is_err());
}

#[test]
fn bloch_plane_wave_and_partial_charge() {
    let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [1., 7., 0.], [0., 1., 6.]]).unwrap();
    let basis = PlaneWaveBasis::at_k(&cell, Hartree(1.0), [0.25, -0.125, 0.0]).unwrap();
    let db = basis.density_basis(&cell).unwrap();
    let grid = FftGridSpec::new(&basis, &db, None).unwrap();
    let index = basis
        .vectors()
        .iter()
        .position(|g| g.miller == [1, 0, 0])
        .unwrap();
    let mut c = vec![Complex64::default(); basis.vectors().len()];
    c[index] = Complex64::new(0., 1.);
    let coefficients = vec![c];
    let block = OrbitalBlock {
        basis: &basis,
        coefficients: &coefficients,
        occupations: &[2.],
        weight: 0.3,
    };
    let tau = kinetic_density(&[block], &cell, grid.dims()).unwrap();
    let expected = 0.3 * basis.momenta()[index].norm_squared / cell.volume_bohr3();
    let error = tau.iter().map(|t| (t - expected).abs()).fold(0., f64::max);
    assert!(error < 1e-15);
    let charge = partial_charge(&coefficients, &grid, cell.volume_bohr3()).unwrap();
    let integral = charge.sum() * cell.volume_bohr3() / charge.len() as f64;
    assert!((integral - 1.).abs() < 1e-14);
    println!(
        "plane-wave tau max error={error:e} Ha/Bohr^3; partial integral={integral:.16} electrons"
    );
}

#[test]
fn real_gaussian_orbital_fft_elf() {
    // Real periodic Gaussian-like orbital: product of (1+cos(2pi*x/L))^4.
    // Compact Fourier support gives exact spectral derivatives of psi and n;
    // L=20 Bohr, occupation=2, normalized coefficients. No sampled-tail aliasing.
    let cell = Cell::from_user_rows_bohr([[20., 0., 0.], [0., 20., 0.], [0., 0., 20.]]).unwrap();
    let basis = PlaneWaveBasis::new(&cell, Hartree(3.)).unwrap();
    let db = basis.density_basis(&cell).unwrap();
    let grid = FftGridSpec::new(&basis, &db, None).unwrap();
    let binomial = [1., 8., 28., 56., 70., 56., 28., 8., 1.];
    let mut c: Vec<_> = basis
        .vectors()
        .iter()
        .map(|g| {
            Complex64::new(
                g.miller
                    .iter()
                    .map(|&m| {
                        if m.abs() <= 4 {
                            binomial[(m + 4) as usize]
                        } else {
                            0.
                        }
                    })
                    .product(),
                0.,
            )
        })
        .collect();
    let norm = c.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
    for z in &mut c {
        *z /= norm;
    }
    let coefficients = vec![c];
    let n = dft_physics::density::from_orbitals(
        &coefficients,
        &[2.],
        &grid,
        cell.volume_bohr3(),
        &mut FftPlan3d::new(grid.dims()).unwrap(),
    )
    .unwrap();
    let e = elf_from_orbitals(
        &n,
        &[OrbitalBlock {
            basis: &basis,
            coefficients: &coefficients,
            occupations: &[2.],
            weight: 1.,
        }],
        &cell,
    )
    .unwrap();
    let error = e
        .iter()
        .zip(&n)
        .filter(|(_, n)| **n > ELF_DENSITY_THRESHOLD)
        .map(|(e, _)| (e - 1.).abs())
        .fold(0., f64::max);
    assert!(error < 1e-10, "ELF error {error:e}");
    println!(
        "real spin-paired Gaussian-like orbital: max ELF error={error:e}; charge={:.16} electrons",
        n.sum() * cell.volume_bohr3() / n.len() as f64
    );
}

#[test]
fn potentials_and_rejections() {
    let cell = Cell::from_user_rows_bohr([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]).unwrap();
    let n = Array3::from_elem([5, 5, 5], 0.1);
    let h = hartree_potential(&n, &cell).unwrap();
    assert!(h.iter().all(|v| v.abs() < 1e-14));
    let local = total_local(&n, &h, &n).unwrap();
    assert!(local.iter().all(|v| (v - 0.2).abs() < 1e-14));
    assert!(total_local(&n, &Array3::zeros([1, 1, 1]), &n).is_err());
}
