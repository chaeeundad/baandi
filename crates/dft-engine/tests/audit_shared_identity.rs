//! Bit-for-bit comparison with the retained pre-P28B audit, including weighted
//! fields. Include the module to exercise crate-private preparation without
//! adding a public API solely for tests.
#[allow(dead_code)]
mod common;
#[allow(dead_code)]
#[path = "../src/energy.rs"]
mod energy;
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_engine::{Functional, ScfError};
use dft_geometry::Cell;
use dft_physics::{
    data::*,
    gga,
    hamiltonian::{NonlocalOperator, Species},
    hartree,
    pbe_pw92::Pw92Prefactor,
};
use energy::{AuditInput, EnergyAudit};

fn identical(a: &EnergyAudit, b: &EnergyAudit) {
    // Destructure exhaustively: a new field must extend this oracle.
    let EnergyAudit {
        xc_floor,
        nlcc,
        kinetic,
        local,
        local_zero_mode,
        nonlocal,
        hartree,
        xc,
        ewald,
        band_sum,
        double_counting,
        band_form,
        harris_foulkes,
        functional_difference_bound,
    } = a;
    macro_rules! bit {
        ($x:expr, $y:expr) => {
            assert_eq!(($x).to_bits(), ($y).to_bits(), "{}", stringify!($y));
        };
    }
    bit!(*kinetic, b.kinetic);
    bit!(*local, b.local);
    bit!(*local_zero_mode, b.local_zero_mode);
    bit!(*nonlocal, b.nonlocal);
    bit!(*hartree, b.hartree);
    bit!(*xc, b.xc);
    bit!(*ewald, b.ewald);
    bit!(*band_sum, b.band_sum);
    bit!(*double_counting, b.double_counting);
    bit!(*band_form, b.band_form);
    bit!(*harris_foulkes, b.harris_foulkes);
    bit!(*functional_difference_bound, b.functional_difference_bound);
    bit!(a.total_hartree(), b.total_hartree());
    bit!(a.total_rydberg(), b.total_rydberg());
    match (xc_floor, &b.xc_floor) {
        (None, None) => (),
        (Some(a), Some(b)) => {
            bit!(a.density_floor, b.density_floor);
            bit!(a.input_negative_charge, b.input_negative_charge);
            bit!(a.output_negative_charge, b.output_negative_charge);
        }
        _ => panic!("floor presence"),
    }
    match (nlcc, &b.nlcc) {
        (None, None) => (),
        (Some(a), Some(b)) => {
            bit!(a.core_electrons, b.core_electrons);
            bit!(a.core_potential_integral_ha, b.core_potential_integral_ha);
            bit!(
                a.valence_potential_integral_ha,
                b.valence_potential_integral_ha
            );
        }
        _ => panic!("core presence"),
    }
}

#[test]
fn every_audit_field_matches_legacy_on_eight_fixtures() {
    for fixture in 0..8 {
        let cell = Cell::from_user_rows_bohr(if fixture % 2 == 0 {
            [[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]]
        } else {
            [[6., 0., 0.], [1., 7., 0.], [0.5, 1., 9.]]
        })
        .unwrap();
        let functional = [Functional::LdaPz, Functional::LdaPw, Functional::Pbe][fixture % 3];
        let basis = PlaneWaveBasis::new(&cell, Hartree(1.)).unwrap();
        let grid = FftGridSpec::new(&basis, &basis.density_basis(&cell).unwrap(), None).unwrap();
        let dims = grid.dims();
        let mut fft = FftPlan3d::new(dims).unwrap();
        let mut n = Array3::from_shape_fn(dims, |(i, j, k)| {
            0.02 + 0.01 * ((i + 2 * j + 3 * k) as f64).sin()
        });
        if fixture >= 6 {
            n[[0, 0, 0]] = -5e-15;
            n[[0, 0, 1]] = -0.0;
        }
        let input = n.mapv(|v| v * 0.9);
        let core = (fixture == 3 || fixture == 4).then(|| Array3::from_elem(dims, 0.01));
        let floor = (fixture == 5 || fixture == 7).then_some(0.015);
        let mut scratch = Array3::zeros(dims);
        let h = hartree::evaluate(&input, &cell, &mut fft).unwrap();
        let (eh, _) = hartree::energy_and_density_reciprocal(&input, &cell, &mut fft).unwrap();
        assert_eq!(eh.0.to_bits(), h.energy.0.to_bits());
        let x = functional
            .evaluate_with_floor(&input, core.as_ref(), &mut scratch, &cell, &mut fft, floor)
            .unwrap();
        let output_xc = if fixture == 2 {
            Some(gga::evaluate_pw92(Pw92Prefactor::Original, &n, &cell, &mut fft).unwrap())
        } else if core.is_some() || floor.is_some() {
            Some(
                functional
                    .evaluate_with_floor(&n, core.as_ref(), &mut scratch, &cell, &mut fft, floor)
                    .unwrap(),
            )
        } else {
            None
        };
        let ionic = Array3::from_shape_fn(dims, |(i, j, k)| {
            Complex64::new(0.03 * ((i + j + k) as f64).cos(), 0.02)
        });
        let pp = common::potential(2., 1., true, 101);
        let mut old = Vec::new();
        let mut new = Vec::new();
        let mut shared = None;
        for (ki, k) in [[0.; 3], [0.125, 0.25, 0.], [-0.25, 0., 0.125]]
            .into_iter()
            .enumerate()
        {
            let kb = PlaneWaveBasis::at_k(&cell, Hartree(1.), k).unwrap();
            let nl = NonlocalOperator::new(
                &kb,
                cell.volume_bohr3(),
                &[Species {
                    potential: &pp,
                    positions_bohr: &[[0., 0., 0.], [1., 2., 1.]],
                }],
            )
            .unwrap();
            let mut c = vec![Complex64::new(0., 0.); kb.vectors().len()];
            c[0] = Complex64::new(1., 0.);
            let orbitals = vec![c];
            let occupations = [0.7 + ki as f64 * 0.2];
            let eigenvalues = [-0.3 + ki as f64 * 0.1];
            let a = || AuditInput {
                cell: &cell,
                basis: &kb,
                grid: &grid,
                ionic_g: &ionic,
                nonlocal: &nl,
                orbitals: &orbitals,
                occupations: &occupations,
                eigenvalues: &eigenvalues,
                eigen_residuals: &[1e-9],
                input_density: &input,
                output_density: &n,
                input_hartree: &h,
                input_xc: &x,
                ewald: -0.123,
            };
            if shared.is_none() {
                shared = Some(
                    energy::prepare_shared(
                        &a(),
                        &mut fft,
                        functional,
                        output_xc.as_ref(),
                        core.as_ref(),
                    )
                    .unwrap(),
                );
            }
            let before = energy::evaluate_with_core_xc(
                a(),
                &mut fft,
                functional,
                output_xc.as_ref(),
                core.as_ref(),
            )
            .unwrap();
            let after = energy::evaluate_shared(a(), shared.as_ref().unwrap()).unwrap();
            identical(&before, &after);
            let weight = [0.1, 0.3, 0.6][ki];
            old.push((weight, before));
            new.push((weight, after));
        }
        identical(
            &energy::weighted_audits(&old),
            &energy::weighted_audits(&new),
        );
    }
}

#[test]
fn energy_only_hartree_preserves_even_nyquist_and_density_validation() {
    let cell = Cell::from_user_rows_bohr([[6., 0., 0.], [2., 7., 0.], [1., 2., 9.]]).unwrap();
    for dims in [[4, 5, 6], [5, 6, 7], [8, 8, 8], [3, 3, 3]] {
        let mut fft = FftPlan3d::new(dims).unwrap();
        let mut n = Array3::from_shape_fn(dims, |(i, j, k)| {
            0.03 + 0.02 * ((i + 3 * j + 2 * k) as f64).cos()
        });
        n[[0, 0, 0]] = -5e-15;
        let before = hartree::evaluate(&n, &cell, &mut fft).unwrap();
        let (after, ng) = hartree::energy_and_density_reciprocal(&n, &cell, &mut fft).unwrap();
        assert_eq!(before.energy.0.to_bits(), after.0.to_bits());
        let mut expected = Array3::zeros(dims);
        fft.to_reciprocal(&n.mapv(|v| Complex64::new(v.max(0.0), 0.0)), &mut expected)
            .unwrap();
        for (a, b) in ng.iter().zip(&expected) {
            assert_eq!(a.re.to_bits(), b.re.to_bits());
            assert_eq!(a.im.to_bits(), b.im.to_bits());
        }
        for invalid in [-2e-14, f64::NAN, f64::INFINITY] {
            n[[0, 0, 0]] = invalid;
            assert!(hartree::evaluate(&n, &cell, &mut fft).is_err());
            assert!(hartree::energy_and_density_reciprocal(&n, &cell, &mut fft).is_err());
        }
    }
}
