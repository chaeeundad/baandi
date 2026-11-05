//! Discrete functional oracle including clipped points. Synthetic only.
use crate::Functional;
use dft_geometry::Cell;
use dft_physics::data::*;

#[test]
fn floor_chain_rule_and_legacy_identity() {
    let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [1., 9., 0.], [0., 1., 10.]]).unwrap();
    let dims = [8, 8, 8];
    let mut n = Array3::from_elem(dims, 0.03);
    n[[1, 2, 3]] = -0.002;
    let core = Array3::from_elem(dims, 0.001);
    let original = n.clone();
    let mut fft = FftPlan3d::new(dims).unwrap();
    let mut scratch = Array3::zeros(dims);
    let weight = cell.volume_bohr3() / n.len() as f64;
    for f in [Functional::LdaPz, Functional::LdaPw, Functional::Pbe] {
        let out = f
            .evaluate_with_floor(&n, Some(&core), &mut scratch, &cell, &mut fft, Some(1e-10))
            .unwrap();
        assert_eq!(out.potential[[1, 2, 3]], 0.0);
        assert_eq!(n, original);
        for index in [[1, 2, 3], [2, 2, 3]] {
            let expected = weight * out.potential[index];
            let mut errors = Vec::new();
            for h in [1e-5, 5e-6, 2.5e-6] {
                let mut plus = n.clone();
                plus[index] += h;
                let mut minus = n.clone();
                minus[index] -= h;
                let ep = f
                    .evaluate_with_floor(
                        &plus,
                        Some(&core),
                        &mut scratch,
                        &cell,
                        &mut fft,
                        Some(1e-10),
                    )
                    .unwrap()
                    .energy
                    .0;
                let em = f
                    .evaluate_with_floor(
                        &minus,
                        Some(&core),
                        &mut scratch,
                        &cell,
                        &mut fft,
                        Some(1e-10),
                    )
                    .unwrap()
                    .energy
                    .0;
                errors.push(((ep - em) / (2. * h) - expected).abs());
            }
            println!("{f:?} index={index:?} derivative errors Ha Bohr^3: {errors:?}");
            assert!(errors[2] < 1e-7);
            assert!(errors[2] < errors[0] * 0.4 || errors[0] < 1e-8);
        }
        let positive = Array3::from_elem(dims, 0.03);
        let a = f
            .evaluate_with_core(&positive, Some(&core), &mut scratch, &cell, &mut fft)
            .unwrap();
        let b = f
            .evaluate_with_floor(&positive, Some(&core), &mut scratch, &cell, &mut fft, None)
            .unwrap();
        assert_eq!(a.energy.0.to_bits(), b.energy.0.to_bits());
        assert!(
            a.potential
                .iter()
                .zip(&b.potential)
                .all(|(a, b)| a.to_bits() == b.to_bits())
        );
    }
}

#[test]
fn clipped_charge_acceptance_budget() {
    for (electrons, limit) in [(2.0, 1e-3), (16.0, 1.6e-3), (100.0, 1e-2)] {
        assert!(crate::ScfProblem::validate_clipped_charge(limit, electrons).is_ok());
        let error =
            crate::ScfProblem::validate_clipped_charge(limit * 1.001, electrons).unwrap_err();
        assert!(error.contains("xc.density_floor") && error.contains("increase cutoff/grid"));
    }
    assert!(crate::ScfProblem::validate_clipped_charge(f64::NAN, 2.0).is_err());
}
