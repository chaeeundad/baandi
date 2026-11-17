//! Issue #17: physical occupation boundary and criterion-specific diagnostics.
use dft_engine::first_block::{self, Criterion, Input};
use dft_numerics::{pool::WorkerPool, real::RealMatrix};
use std::sync::Mutex;
static LOCK: Mutex<()> = Mutex::new(());

fn fixture(criterion: Criterion, empty_residual: f64) -> first_block::Record {
    let mut x = RealMatrix::zeros(10, 5).unwrap();
    let mut hx = RealMatrix::zeros(10, 5).unwrap();
    for j in 0..5 {
        x[(j, j)] = 1.;
        hx[(j, j)] = [0., 1., 4., 9., 16.][j];
        hx[(j + 5, j)] = [0.125, 0.25, empty_residual, 0., 0.][j];
    }
    first_block::take_records();
    first_block::evaluate::<dft_mpi::SimRank>(
        None,
        Input {
            x: &x,
            hx: Some(&hx),
            bands: 5,
            occupations: &[2., 2., 0., 0., 0.],
            electrons: 4.,
            criterion,
            tolerance: 1e-7,
            strict_bands: 3,
            buffer_tolerance: 10.,
            capacity: 10,
            pool: WorkerPool::new(1).unwrap(),
        },
        &mut |_, _| panic!("kept HX must be used"),
    )
    .unwrap();
    first_block::take_records().pop().unwrap()
}

#[test]
fn physical_boundary_spectrum_and_criterion() {
    let _lock = LOCK.lock().unwrap();
    let standard = fixture(Criterion::RitzChange, 1.);
    assert_eq!(standard.occupied, 2);
    assert_eq!(standard.electrons, 4.);
    assert_eq!(standard.occupations, vec![2., 2., 0., 0., 0.]);
    assert_eq!(standard.mandatory_roots, 3);
    assert_eq!(standard.empty_buffer_roots, 2);
    assert_eq!(standard.boundary_gap, 3.);
    assert_eq!(standard.singular_values, vec![0.25, 0.125]);
    assert_eq!(
        standard.singular_values,
        fixture(Criterion::RitzChange, 100.).singular_values
    );
    // No preceding RR iteration exists at entry: Ritz changes are unavailable,
    // not inferred from residuals or the previous SCF Hamiltonian.
    assert_eq!(standard.ritz_changes, None);
    assert_eq!(standard.residual_unconverged_all, None);
    assert_eq!(standard.rank_above_residual_tolerance, None);
    assert_eq!(standard.root_thresholds, vec![1e-7, 1e-7, 10., 10., 10.]);
    let strict = fixture(Criterion::ResidualNorm, 1.);
    assert_eq!(strict.root_thresholds, vec![1e-7, 1e-7, 1e-7, 10., 10.]);
    assert_eq!(strict.residual_unconverged_all, Some(3));
    assert_eq!(strict.residual_unconverged_occupied, Some(2));
    assert_eq!(strict.ritz_first_correction_capacity, None);
    assert_eq!(strict.rank_above_residual_tolerance, Some(2));
}

#[allow(dead_code)]
mod common;

#[test]
fn fixed_steps_preserve_production_bits_and_work() {
    use dft_engine::{GammaReal, ScfOptions, ScfProblem, SpeciesPotential};
    use dft_geometry::{Atom, Cell};
    let _lock = LOCK.lock().unwrap();
    let pp = common::potential(4., 1., true, 401);
    let cell = Cell::from_user_rows_bohr([[4., 0., 0.], [0.5, 5., 0.], [0.2, 0.4, 6.]]).unwrap();
    for standard in [false, true] {
        let dir = std::env::temp_dir().join(format!("issue17-{}-{standard}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let run = |on: bool| {
            dft_engine::snapshot::force_dir(Some(on.then(|| dir.clone())));
            first_block::force(Some(on));
            first_block::take_records();
            let mut p = ScfProblem::new(
                &cell,
                &[Atom::new("Xx", [0.13, 0.21, 0.32])],
                &[SpeciesPotential {
                    label: "Xx",
                    potential: &pp,
                }],
                ScfOptions {
                    gamma_real: GammaReal::On,
                    cutoff_hartree: 2.5,
                    bands: Some(5),
                    standard: standard.then(dft_engine::StandardProfile::default),
                    ..ScfOptions::strict_reference()
                },
            )
            .unwrap();
            let mut density = p.uniform_density();
            let mut warm = None;
            let mut steps = Vec::new();
            for _ in 0..3 {
                let step = p.step(&density, warm.as_ref(), 1e-9).unwrap();
                density = step.density.clone();
                warm = Some(step.solution.eigenvectors.clone());
                steps.push(step);
            }
            (steps, first_block::take_records())
        };
        let (plain, absent) = run(false);
        assert!(absent.is_empty());
        let (observed, records) = run(true);
        assert_eq!(records.len(), 2);
        for (a, b) in plain.iter().zip(&observed) {
            assert_eq!(
                a.energy.total_hartree().to_bits(),
                b.energy.total_hartree().to_bits()
            );
            assert_eq!(
                a.density.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                b.density.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
            );
            // Includes residuals, vectors, iteration count and H-application count.
            assert_eq!(a.solution, b.solution);
        }
        dft_engine::snapshot::force_dir(Some(None));
        let manifests: Vec<_> = std::fs::read_dir(&dir)
            .unwrap()
            .map(|e| e.unwrap().path())
            .filter(|p| p.extension().is_some_and(|e| e == "json"))
            .collect();
        assert_eq!(manifests.len(), 3);
        for path in manifests {
            let text = std::fs::read_to_string(path).unwrap();
            for expected in [
                "\"schema\": \"baandi-snapshot-1\"",
                "\"occupation_metadata\": \"physical-v2\"",
                "\"occupied\": 2,",
                "\"electrons\": 4.00000000000000000e0,",
                "\"occupations\": [2.0, 2.0, 0.0, 0.0, 0.0]",
                "\"mandatory_roots\": 3,",
                "\"density_carrying_roots\": 2,",
                "\"empty_buffer_roots\": 2,",
            ] {
                assert!(text.contains(expected), "missing {expected} in {text}");
            }
            assert!(text.contains(&format!(
                "\"convergence_criterion\": \"{}\"",
                if standard {
                    "ritz_change"
                } else {
                    "residual_norm"
                }
            )));
        }
        std::fs::remove_dir_all(&dir).unwrap();
        for r in &records {
            assert_eq!(r.occupied, 2);
            assert_eq!(r.electrons, 4.);
            assert_eq!(r.mandatory_roots, 3);
            assert_eq!(r.singular_values.len(), 2);
            assert_eq!(r.boundary_gap, r.ritz_values[2] - r.ritz_values[1]);
            assert_eq!(
                r.criterion,
                if standard {
                    Criterion::RitzChange
                } else {
                    Criterion::ResidualNorm
                }
            );
            assert_eq!(
                r.diagnostic_operator_applications,
                if r.kept_hx { 0 } else { 5 }
            );
        }
        println!(
            "standard={standard}: 3 fixed steps, energy/density/solver identical; 2 occupied, 4 electrons, 3 mandatory, 2 buffers"
        );
    }
    first_block::force(None);
}
