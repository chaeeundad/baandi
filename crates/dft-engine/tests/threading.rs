#[allow(dead_code)]
mod common;
use dft_engine::{
    OccupationPolicy, PhaseTimings, ScfOptions, ScfProblem, ScfResult, SpeciesPotential,
};
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use std::time::Duration;

// Debug prints every scalar field in these engine records, with round-trippable
// f64 representations (including signed zero). Arrays/matrices are also checked
// explicitly by bits below so array display truncation cannot hide a mismatch.
fn normalized(mut r: ScfResult) -> ScfResult {
    r.setup_time = Duration::ZERO;
    r.total_time = Duration::ZERO;
    for h in &mut r.history {
        h.timings = PhaseTimings::default();
        if let Some(c) = &mut h.fixed_point_certificate {
            c.timings = PhaseTimings::default();
        }
    }
    r
}
fn assert_identical(a: &ScfResult, b: &ScfResult) {
    assert_eq!(format!("{a:?}"), format!("{b:?}"));
    assert_eq!(a.density.shape(), b.density.shape());
    for (x, y) in a.density.iter().zip(&b.density) {
        assert_eq!(x.to_bits(), y.to_bits());
    }
    for (x, y) in a.k_points.iter().zip(&b.k_points) {
        assert_eq!(x.solution, y.solution);
        for j in 0..x.solution.eigenvectors.cols() {
            for (v, w) in x
                .solution
                .eigenvectors
                .column(j)
                .iter()
                .zip(y.solution.eigenvectors.column(j))
            {
                assert_eq!(v.re.to_bits(), w.re.to_bits());
                assert_eq!(v.im.to_bits(), w.im.to_bits());
            }
        }
    }
}

#[test]
fn scf_thread_counts_are_bit_identical() {
    let pp = common::potential(2.0, 1.0, true, 401);
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let cell = Cell::from_user_rows_bohr([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    for nk in [1, 2, 5] {
        let points: Vec<_> = (0..nk)
            .map(|i| KPoint {
                fractional: [0.025 * i as f64, 0.015 * i as f64, 0.0],
                weight: 1.0 / nk as f64,
            })
            .collect();
        let run = |threads| {
            normalized(
                ScfProblem::with_kpoints(
                    &cell,
                    &atoms,
                    &species,
                    ScfOptions {
                        threads,
                        cutoff_hartree: 2.5,
                        bands: Some(3),
                        ..Default::default()
                    },
                    &points,
                    OccupationPolicy::Insulator {
                        minimum_gap_hartree: 1e-7,
                    },
                )
                .unwrap()
                .run()
                .unwrap(),
            )
        };
        // threads=1 executes the retained pre-threading serial solve loop.
        let reference = run(1);
        assert!(
            reference
                .history
                .iter()
                .any(|h| h.fixed_point_certificate.is_some())
        );
        for threads in [0, 2, 3, 4, 8] {
            assert_identical(&reference, &run(threads));
        }
        println!(
            "{nk} k points: threads 0,1,2,3,4,8 bit-identical, energy {:.17} Ha",
            reference.total_hartree()
        );
    }
}
