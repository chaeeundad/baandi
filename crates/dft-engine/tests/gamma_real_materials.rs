//! External-file gates. Never stores a pseudopotential in the repository.
#[allow(dead_code)]
mod common;
use dft_engine::{GammaReal, Mixing, ScfOptions, ScfProblem, SpeciesPotential};
use dft_geometry::{Atom, Cell};

#[test]
#[ignore = "requires DFT_SI_UPF; release only; includes the design SCF-count gate"]
fn real_si_accuracy_determinism_and_iteration_gate() {
    material_matrix(false);
}

#[test]
#[ignore = "requires DFT_SI_UPF; release standard-profile matrix"]
fn standard_si_accuracy_and_determinism() {
    material_matrix(true);
}

fn material_matrix(standard: bool) {
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
        let solve = |mode, threads, reference: Option<(f64, f64, f64)>| {
            ScfProblem::new(
                &cell,
                &atoms,
                &[SpeciesPotential {
                    label: "Si",
                    potential: &pp,
                }],
                ScfOptions {
                    gamma_real: mode,
                    standard: (standard && reference.is_none())
                        .then(dft_engine::StandardProfile::default),
                    energy_tolerance: reference
                        .map_or(ScfOptions::default().energy_tolerance, |v| v.0),
                    density_tolerance: reference
                        .map_or(ScfOptions::default().density_tolerance, |v| v.1),
                    eigen_tolerance: reference
                        .map_or(ScfOptions::default().eigen_tolerance, |v| v.2),
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
        let run = |mode, threads| solve(mode, threads, None);
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
            "MATERIAL {case} standard={standard} E={:.17e} dE={de:.17e} eigen={ev:.17e} density={dn:.17e} SCF={}/{} work_complex={:?} work_real={:?} extra_certificate_H_lower_bound={extra_lower_bound}",
            real.total_hartree(),
            complex.history.len(),
            real.history.len(),
            work(&complex),
            work(&real)
        );
        // Frozen before real comparisons: design energy, measured complex noise.
        if standard {
            for r in [&complex, &real] {
                assert!(!r.certified);
                assert!(r.k_points.iter().all(|k| !k.solution.residuals_certified));
                assert!(
                    r.history
                        .iter()
                        .all(|h| h.fixed_point_certificate.is_none())
                );
                assert!(r.history.iter().all(|h| {
                    h.standard_error_ha
                        .is_some_and(|e| e.is_finite() && e >= 0.)
                }));
                assert!(
                    r.history.last().unwrap().standard_error_ha.unwrap()
                        < dft_engine::StandardProfile::default().threshold_ha
                );
                assert!(r.history.iter().skip(1).all(|h| !h.standard_redo));
                for h in r.history.windows(2) {
                    assert_eq!(
                        h[1].requested_eigen_tolerance,
                        dft_engine::StandardProfile::ratchet(
                            h[0].requested_eigen_tolerance,
                            h[0].standard_error_ha.unwrap(),
                            r.occupations.iter().sum()
                        )
                    );
                }
                println!(
                    "STANDARD_FIELDS {case} real={} redo={} error={:.17e} Ha",
                    r.gamma_real_used,
                    r.history[0].standard_redo,
                    r.history.last().unwrap().standard_error_ha.unwrap()
                );
            }
            // Same P17A tight/tighter settings and architect-approved delivered-
            // accuracy gate as the synthetic standard matrix; strict is unchanged.
            let tight = solve(GammaReal::Off, 1, Some((1e-12, 1e-9, 1e-10)));
            let reference = solve(GammaReal::Off, 1, Some((1e-13, 1e-10, 1e-11)));
            for r in [&tight, &reference] {
                assert!(r.certified && r.standard.is_none() && !r.gamma_real_used);
                assert!(r.k_points.iter().all(|k| k.solution.residuals_certified));
            }
            let plateau = (tight.total_hartree() - reference.total_hartree()).abs();
            assert!(plateau <= 1e-12, "reference plateau {case}: {plateau} Ha");
            let errors = |r: &dft_engine::ScfResult| {
                let eigen = r.eigenvalues[..occupied]
                    .iter()
                    .zip(&reference.eigenvalues)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0., f64::max);
                let all_eigen = r
                    .eigenvalues
                    .iter()
                    .zip(&reference.eigenvalues)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0., f64::max);
                let density = (r
                    .density
                    .iter()
                    .zip(&reference.density)
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
                    * cell.volume_bohr3()
                    / r.density.len() as f64)
                    .sqrt();
                (
                    (r.total_hartree() - reference.total_hartree()).abs(),
                    eigen,
                    all_eigen,
                    density,
                )
            };
            let c = errors(&complex);
            let r = errors(&real);
            let relative_cap = (3. * c.0).max(1e-10);
            println!(
                "ACCURACY {case} tight={:.17e} reference={:.17e} plateau={plateau:.17e} complex_E={:.17e} complex_error={:.17e} real_error={:.17e} relative_cap={relative_cap:.17e} complex_eigen={:.17e} real_eigen={:.17e} complex_all_eigen={:.17e} real_all_eigen={:.17e} complex_density={:.17e} real_density={:.17e} reference_SCF={}",
                tight.total_hartree(),
                reference.total_hartree(),
                complex.total_hartree(),
                c.0,
                r.0,
                c.1,
                r.1,
                c.2,
                r.2,
                c.3,
                r.3,
                reference.history.len()
            );
            if r.0 > 5e-9 || r.0 > relative_cap {
                count_failures.push(format!(
                    "{case}: complex error {} Ha, real error {} Ha, relative cap {relative_cap} Ha",
                    c.0, r.0
                ));
            }
        } else {
            assert!(de < 1e-10 && ev < 3e-8 && dn < 7e-7);
        }
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
            assert_eq!(real.k_points[0].solution, other.k_points[0].solution);
            for (a, b) in real.history.iter().zip(&other.history) {
                assert_eq!(a.standard_redo, b.standard_redo);
                assert_eq!(a.standard_error_ha, b.standard_error_ha);
                assert_eq!(a.requested_eigen_tolerance, b.requested_eigen_tolerance);
            }
        }
        // Section 11 explicitly allows only one SCF iteration difference on
        // realistic fixtures. Keep failures visible, without retuning a solver.
        if !standard && real.history.len().abs_diff(complex.history.len()) > 1 {
            count_failures.push(format!(
                "{case}: {} -> {}",
                complex.history.len(),
                real.history.len()
            ));
        }
    }
    assert!(
        count_failures.is_empty(),
        "material acceptance gate failed: {count_failures:?}"
    );
}
