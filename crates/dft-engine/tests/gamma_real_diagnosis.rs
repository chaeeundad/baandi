//! Diagnostic only: no production start, tolerance, or correction order changes.
#[allow(dead_code)]
mod common;
use dft_basis::{FftGridSpec, HalfSphere, PlaneWaveBasis};
use dft_engine::{GammaReal, Mixing, ScfOptions, ScfProblem, SpeciesPotential};
use dft_geometry::Atom;
use dft_numerics::linalg::ComplexMatrix;

#[test]
#[ignore = "external Si UPF; cargo test --release -p dft-engine --test gamma_real_diagnosis -- --ignored --nocapture"]
fn primitive_start_and_history_diagnosis() {
    let pp =
        dft_upf::parse_upf(&std::fs::read_to_string(std::env::var("DFT_SI_UPF").unwrap()).unwrap())
            .unwrap();
    let cell = common::cell();
    let atoms = [Atom::new("Si", [0.; 3]), Atom::new("Si", [0.25; 3])];
    let options = |mode| ScfOptions {
        gamma_real: mode,
        cutoff_hartree: 15.,
        max_iterations: 100,
        mixing: Mixing::Pulay {
            beta: 0.5,
            history: 6,
        },
        ..Default::default()
    };
    let make = |mode| {
        ScfProblem::new(
            &cell,
            &atoms,
            &[SpeciesPotential {
                label: "Si",
                potential: &pp,
            }],
            options(mode),
        )
        .unwrap()
    };
    for mode in [GammaReal::Off, GammaReal::On] {
        let result = make(mode).run().unwrap();
        for h in &result.history {
            println!(
                "HISTORY {mode:?} i={} tol={:.9e} density={:.9e} steps={} H={}",
                h.iteration,
                h.requested_eigen_tolerance,
                h.density_residual,
                h.davidson_iterations,
                h.operator_applications
            );
        }
    }
    let basis = PlaneWaveBasis::new(&cell, dft_units::Hartree(15.)).unwrap();
    let density = PlaneWaveBasis::new(&cell, dft_units::Hartree(60.)).unwrap();
    let grid = FftGridSpec::new(&basis, &density, None).unwrap();
    let half = HalfSphere::new(&basis, &grid).unwrap();
    let kinetic: Vec<_> = basis
        .vectors()
        .iter()
        .map(|g| 0.5 * g.norm_squared)
        .collect();
    let mut diagonal = vec![0.; kinetic.len()];
    for (i, &p) in half.plus().iter().enumerate().skip(1) {
        diagonal[2 * i - 1] = kinetic[p as usize];
        diagonal[2 * i] = kinetic[p as usize];
    }
    let real = dft_solvers::real_davidson::free_electron_initial_guess(&diagonal, 8).unwrap();
    let mut initial = ComplexMatrix::zeros(kinetic.len(), 8).unwrap();
    for j in 0..8 {
        half.unpack(
            &real.as_slice()[j * kinetic.len()..(j + 1) * kinetic.len()],
            initial.column_mut(j),
        )
        .unwrap();
    }
    let mut order: Vec<_> = (0..kinetic.len()).collect();
    order.sort_by(|&i, &j| kinetic[i].total_cmp(&kinetic[j]).then(i.cmp(&j)));
    let edge = kinetic[order[7]];
    let below = kinetic.iter().filter(|&&x| x < edge - 1e-12).count();
    let shell = kinetic
        .iter()
        .filter(|&&x| (x - edge).abs() < 1e-12)
        .count();
    let overlap: f64 = (0..8)
        .map(|j| {
            order[..8]
                .iter()
                .map(|&i| initial.column(j)[i].norm_sqr())
                .sum::<f64>()
        })
        .sum();
    println!(
        "START dimension={} edge_Ha={edge:.17e} below={below} shell={shell} selected_in_shell={} projector_distance_squared={:.17e}",
        kinetic.len(),
        8 - below,
        16. - 2. * overlap
    );
    let mut first_density = None;
    let mut common_start_density = None;
    for (mode, same_start) in [
        (GammaReal::Off, false),
        (GammaReal::On, false),
        (GammaReal::Off, true),
        (GammaReal::On, true),
    ] {
        let mut problem = make(mode);
        let step = problem
            .step(
                &problem.uniform_density(),
                same_start.then_some(&initial),
                1e-2,
            )
            .unwrap();
        if let Some(ref reference) = first_density {
            let difference = step
                .density
                .iter()
                .zip(reference)
                .map(|(x, y): (&f64, &f64)| (x - y).powi(2))
                .sum::<f64>();
            println!(
                "FIRST_DENSITY {mode:?} same_ORC_start={same_start} delta_L2={:.17e}",
                (difference * cell.volume_bohr3() / step.density.len() as f64).sqrt()
            );
        } else {
            first_density = Some(step.density.iter().copied().collect::<Vec<_>>());
        }
        if same_start {
            if let Some(ref reference) = common_start_density {
                let difference = step
                    .density
                    .iter()
                    .zip(reference)
                    .map(|(x, y): (&f64, &f64)| (x - y).powi(2))
                    .sum::<f64>();
                let delta = (difference * cell.volume_bohr3() / step.density.len() as f64).sqrt();
                println!("COMMON_START density_L2={delta:.17e}");
                assert!(delta < 1e-12);
            } else {
                common_start_density = Some(step.density.iter().copied().collect::<Vec<_>>());
            }
        }
        println!(
            "FIRST {mode:?} same_ORC_start={same_start} steps={} H={} eigen={:?} residual={:?} buffer={:?}",
            step.solution.iterations,
            step.solution.operator_applications,
            step.solution.eigenvalues,
            step.solution.residuals,
            step.solution.buffer
        );
    }
}
