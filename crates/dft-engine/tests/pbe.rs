#[allow(dead_code)]
mod common;
use dft_engine::*;
use dft_geometry::{Atom, Cell, kpoints::monkhorst_pack};
fn cell() -> Cell {
    Cell::from_user_rows_bohr([[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]]).unwrap()
}
#[test]
fn functional_labels() {
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    for label in ["PBE", " SLA PW PBX PBC ", "SLA PZ"] {
        let pp = common::potential_with_functional(2.0, 1.0, false, 401, label);
        let species = [SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        for allow in [false, true] {
            let options = ScfOptions {
                functional: Functional::Pbe,
                allow_functional_override: allow,
                cutoff_hartree: 1.5,
                ..Default::default()
            };
            assert_eq!(
                ScfProblem::new(&cell(), &atoms, &species, options).is_ok(),
                allow || label != "SLA PZ"
            );
        }
    }
    for label in ["PBE0", "PBESOL", "SLA PW", "PBE NOGX"] {
        assert!(!Functional::Pbe.matches_upf(label));
    }
}
#[test]
#[ignore = "PBE Gamma and k-grid SCF; cargo test --offline --release -p dft-engine --test pbe -- --ignored --nocapture"]
fn synthetic_scf() {
    let a = common::potential_with_functional(1.0, 0.55, true, 401, "PBE");
    let b = common::potential_with_functional(1.0, 1.15, true, 401, "PBE");
    let species = [
        SpeciesPotential {
            label: "A",
            potential: &a,
        },
        SpeciesPotential {
            label: "B",
            potential: &b,
        },
    ];
    let atoms = [
        Atom::new("A", [0.11, 0.21, 0.17]),
        Atom::new("B", [0.39, 0.57, 0.68]),
    ];
    for dims in [[1; 3], [2; 3]] {
        let points =
            monkhorst_pack(dims, if dims == [1; 3] { [0.0; 3] } else { [0.5; 3] }).unwrap();
        let options = ScfOptions {
            functional: Functional::Pbe,
            cutoff_hartree: 2.5,
            bands: Some(3),
            ..Default::default()
        };
        let r = ScfProblem::with_kpoints(
            &cell(),
            &atoms,
            &species,
            options.clone(),
            &points,
            OccupationPolicy::Insulator {
                minimum_gap_hartree: 1e-7,
            },
        )
        .unwrap()
        .run()
        .unwrap();
        let e = &r.energy;
        let last = r.history.last().unwrap();
        let band_bound = 2.0 * last.max_eigen_residual + 1e-12;
        let harris_error = (e.total_hartree() - e.harris_foulkes).abs();
        let band_error = (e.total_hartree() - e.band_form).abs();
        assert!(e.total_hartree().is_finite());
        assert!(band_error <= band_bound);
        assert!(harris_error <= e.functional_difference_bound + 1e-12);
        assert!((r.charge - 2.0).abs() < 1e-10);
        assert!(last.fixed_point_certificate.is_some());
        println!(
            "PBE {dims:?}: iterations={}, E={:.16e} Ha, charge={:.16e} electrons, band error={band_error:.9e} <= {band_bound:.9e} Ha, Harris error={harris_error:.9e} <= {:.9e} Ha, density residual={:.9e}",
            r.history.len(),
            e.total_hartree(),
            r.charge,
            e.functional_difference_bound,
            last.density_residual
        );
    }
}
