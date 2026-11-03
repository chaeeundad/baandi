#[allow(dead_code)]
mod common;
use dft_engine::{NscfError, ScfOptions, ScfProblem, SpeciesPotential, gaussian_dos, nscf};
use dft_geometry::{Atom, Cell, kpoints::KPoint};

#[test]
fn fixed_density_and_threads() {
    let cell = Cell::from_user_rows_bohr([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]).unwrap();
    let pp = common::potential(2.0, 1.0, false, 401);
    let species = [SpeciesPotential {
        label: "X",
        potential: &pp,
    }];
    let atoms = [Atom::new("X", [0.13, 0.21, 0.32])];
    let options = ScfOptions {
        cutoff_hartree: 1.5,
        eigen_tolerance: 1e-10,
        density_tolerance: 1e-11,
        energy_tolerance: 1e-12,
        ..Default::default()
    };
    let points = vec![
        KPoint {
            fractional: [0.; 3],
            weight: 0.5,
        },
        KPoint {
            fractional: [0.1, 0., 0.],
            weight: 0.5,
        },
    ];
    let scf = ScfProblem::with_kpoints(
        &cell,
        &atoms,
        &species,
        options.clone(),
        &points,
        dft_engine::OccupationPolicy::Insulator {
            minimum_gap_hartree: 1e-7,
        },
    )
    .unwrap()
    .run()
    .unwrap();
    let density = scf.density.clone();
    let one = nscf(&cell, &atoms, &species, &density, &points, options.clone()).unwrap();
    let four = nscf(
        &cell,
        &atoms,
        &species,
        &density,
        &points,
        ScfOptions {
            threads: 4,
            ..options.clone()
        },
    )
    .unwrap();
    let error = one
        .solutions
        .iter()
        .zip(&scf.k_points)
        .flat_map(|(a, b)| a.solution.eigenvalues.iter().zip(&b.solution.eigenvalues))
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    println!("synthetic fixed-density versus SCF maximum eigenvalue difference {error:.16e} Ha");
    assert!(error < 1e-8);
    assert_eq!(density, scf.density);
    for (a, b) in one.solutions.iter().zip(&four.solutions) {
        assert_eq!(a.solution.eigenvalues, b.solution.eigenvalues);
        assert_eq!(a.solution.residuals, b.solution.residuals);
        assert!(a.solution.residuals.iter().all(|&r| r <= 1e-10));
    }
    assert!(one.converged && four.converged);
    for threads in [0, 2, 9] {
        let other = nscf(
            &cell,
            &atoms,
            &species,
            &density,
            &points,
            ScfOptions {
                threads,
                ..options.clone()
            },
        )
        .unwrap();
        for (a, b) in one.solutions.iter().zip(other.solutions) {
            assert_eq!(a.solution.eigenvalues, b.solution.eigenvalues);
        }
    }
    let failure = nscf(
        &cell,
        &atoms,
        &species,
        &density,
        &points,
        ScfOptions {
            davidson_max_iterations: 1,
            ..options.clone()
        },
    )
    .unwrap_err();
    assert!(matches!(failure, NscfError::KPoint { index: 0, .. }));
    let bad = nscf(
        &cell,
        &atoms,
        &species,
        &density,
        &points,
        ScfOptions {
            bands: Some(1000),
            ..options
        },
    )
    .unwrap_err();
    assert!(matches!(bad, NscfError::KPoint { index: 0, .. }));
}

#[test]
fn gaussian_normalization_peaks_and_weights() {
    let energies: Vec<_> = (-1000..=1000).map(|i| i as f64 * 0.005).collect();
    let spectrum = vec![vec![-1., 1.], vec![-1., 1.]];
    let rows = gaussian_dos(&spectrum, &[0.25, 0.75], &energies, 0.1).unwrap();
    let integral = rows
        .windows(2)
        .map(|p| (p[1][0] - p[0][0]) * (p[0][1] + p[1][1]) / 2.)
        .sum::<f64>();
    assert!((integral - 4.).abs() < 1e-6);
    assert!((rows.last().unwrap()[2] - 4.).abs() < 1e-12);
    for target in [-1., 1.] {
        let peak = rows
            .iter()
            .filter(|r| (r[0] - target).abs() < 0.5)
            .max_by(|a, b| a[1].total_cmp(&b[1]))
            .unwrap();
        assert!((peak[0] - target).abs() < 1e-12);
    }
    assert!(gaussian_dos(&spectrum, &[0.2, 0.7], &energies, 0.1).is_err());
    println!(
        "synthetic Gaussian integral {integral:.16e} states/cell (expected 4; includes spin factor 2)"
    );
}

#[test]
#[ignore = "real Si; DFT_SI_UPF=/path/to/Si.pz-vbc.UPF cargo test --offline --release -p dft-engine --test nscf real_si -- --ignored --nocapture"]
fn real_si() {
    let path =
        std::env::var("DFT_SI_UPF").expect("set DFT_SI_UPF to the read-only Si validation UPF");
    let pp = dft_upf::from_file(path).unwrap();
    let cell = common::cell();
    let atoms = common::atoms();
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: 5.0,
        bands: Some(12),
        eigen_tolerance: 1e-10,
        density_tolerance: 1e-10,
        energy_tolerance: 1e-12,
        davidson_max_iterations: 1200,
        davidson_subspace: 64,
        ..Default::default()
    };
    let scf = ScfProblem::new(&cell, &atoms, &species, options.clone())
        .unwrap()
        .run()
        .unwrap();
    let points = [KPoint {
        fractional: [0.; 3],
        weight: 1.,
    }];
    let one = nscf(
        &cell,
        &atoms,
        &species,
        &scf.density,
        &points,
        options.clone(),
    )
    .unwrap();
    let four = nscf(
        &cell,
        &atoms,
        &species,
        &scf.density,
        &points,
        ScfOptions {
            threads: 4,
            ..options
        },
    )
    .unwrap();
    let delta = one.solutions[0]
        .solution
        .eigenvalues
        .iter()
        .zip(&scf.eigenvalues)
        .map(|(a, b)| (a - b).abs())
        .fold(0., f64::max);
    assert!(delta < 1e-8);
    assert_eq!(
        one.solutions[0].solution.eigenvalues,
        four.solutions[0].solution.eigenvalues
    );
    println!(
        "Si Gamma fixed-density/SCF max difference {delta:.16e} Ha; threads 1/4 exactly equal"
    );
}
