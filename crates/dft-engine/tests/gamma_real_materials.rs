//! External-file gates. Never stores a pseudopotential in the repository.
#[allow(dead_code)]
mod common;
use dft_engine::{GammaReal, Mixing, ScfOptions, ScfProblem, SpeciesPotential};
use dft_geometry::{Atom, Cell};

#[test]
#[ignore = "requires DFT_SI_UPF; release only; includes the design SCF-count gate"]
fn real_si_accuracy_determinism_and_iteration_gate() {
    let path = std::env::var("DFT_SI_UPF").expect("set DFT_SI_UPF to the external Si.pz-vbc.UPF");
    let pp = dft_upf::parse_upf(&std::fs::read_to_string(path).unwrap()).unwrap();
    let a = 10.261_212_856_717_933;
    let mut count_failures = Vec::new();
    for case in ["gamma", "si8"] {
        let (cell, positions) = if case == "gamma" {
            (common::cell(), vec![[0.; 3], [0.25; 3]])
        } else {
            (
                Cell::from_user_rows_bohr([[a, 0., 0.], [0., a, 0.], [0., 0., a]]).unwrap(),
                vec![
                    [0., 0., 0.],
                    [0., 0.5, 0.5],
                    [0.5, 0., 0.5],
                    [0.5, 0.5, 0.],
                    [0.25, 0.25, 0.25],
                    [0.25, 0.75, 0.75],
                    [0.75, 0.25, 0.75],
                    [0.75, 0.75, 0.25],
                ],
            )
        };
        let atoms: Vec<_> = positions.into_iter().map(|p| Atom::new("Si", p)).collect();
        let run = |mode, threads| {
            ScfProblem::new(
                &cell,
                &atoms,
                &[SpeciesPotential {
                    label: "Si",
                    potential: &pp,
                }],
                ScfOptions {
                    gamma_real: mode,
                    threads,
                    cutoff_hartree: 15.,
                    max_iterations: 100,
                    mixing: Mixing::Pulay {
                        beta: 0.5,
                        history: 6,
                    },
                    ..Default::default()
                },
            )
            .unwrap()
            .run()
            .unwrap()
        };
        let complex = run(GammaReal::Off, 1);
        let real = run(GammaReal::On, 1);
        let occupied = real.occupations.iter().filter(|&&x| x > 0.).count();
        let de = (real.total_hartree() - complex.total_hartree()).abs();
        let ev = real.eigenvalues[..occupied]
            .iter()
            .zip(&complex.eigenvalues)
            .map(|(x, y)| (x - y).abs())
            .fold(0., f64::max);
        let dn = (real
            .density
            .iter()
            .zip(&complex.density)
            .map(|(x, y)| (x - y).powi(2))
            .sum::<f64>()
            * cell.volume_bohr3()
            / real.density.len() as f64)
            .sqrt();
        let work = |r: &dft_engine::ScfResult| {
            (
                r.history
                    .iter()
                    .map(|h| h.davidson_iterations)
                    .sum::<usize>(),
                r.history
                    .iter()
                    .map(|h| h.operator_applications)
                    .sum::<usize>(),
                r.history
                    .iter()
                    .filter_map(|h| h.fixed_point_certificate.as_ref())
                    .map(|h| h.davidson_iterations)
                    .sum::<usize>(),
                r.history
                    .iter()
                    .filter_map(|h| h.fixed_point_certificate.as_ref())
                    .map(|h| h.operator_applications)
                    .sum::<usize>(),
            )
        };
        // R4 performs at least one bands-wide true-residual application in each
        // solve for which the complex caller would permit recurrence acceptance.
        let extra_lower_bound = real.occupations.len()
            * real
                .history
                .iter()
                .filter(|h| h.requested_eigen_tolerance > 1e-8)
                .count();
        println!(
            "MATERIAL {case} E={:.17e} dE={de:.17e} eigen={ev:.17e} density={dn:.17e} SCF={}/{} work_complex={:?} work_real={:?} extra_certificate_H_lower_bound={extra_lower_bound}",
            real.total_hartree(),
            complex.history.len(),
            real.history.len(),
            work(&complex),
            work(&real)
        );
        // Frozen before real comparisons: design energy, measured complex noise.
        assert!(de < 1e-10 && ev < 3e-8 && dn < 7e-7);
        assert!(real.gamma_real_used && real.gamma_real_fallback_reason.is_none());
        for threads in [2, 4, 8, 12] {
            let other = run(GammaReal::On, threads);
            assert_eq!(
                real.total_hartree().to_bits(),
                other.total_hartree().to_bits()
            );
            assert!(
                real.density
                    .iter()
                    .zip(&other.density)
                    .all(|(x, y)| x.to_bits() == y.to_bits())
            );
            assert_eq!(real.history.len(), other.history.len());
            assert_eq!(work(&real), work(&other));
        }
        // Section 11 explicitly allows only one SCF iteration difference on
        // realistic fixtures. Keep failures visible, without retuning a solver.
        if real.history.len().abs_diff(complex.history.len()) > 1 {
            count_failures.push(format!(
                "{case}: {} -> {}",
                complex.history.len(),
                real.history.len()
            ));
        }
    }
    assert!(
        count_failures.is_empty(),
        "design SCF-count gate failed: {count_failures:?}"
    );
}
