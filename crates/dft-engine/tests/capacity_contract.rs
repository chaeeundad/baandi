#[allow(dead_code)]
mod common;
#[allow(dead_code)]
#[path = "common/checkpoint.rs"]
mod fixture;
use dft_engine::*;
use dft_geometry::{Atom, Cell, kpoints::KPoint};

fn problem(multiplier: usize) -> Result<ScfProblem, ScfError> {
    let pp = common::potential(2.0, 1.0, true, 401);
    let cell = Cell::from_user_rows_bohr([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]).unwrap();
    ScfProblem::with_kpoints(
        &cell,
        &[Atom::new("Xx", [0.13, 0.21, 0.32])],
        &[SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }],
        ScfOptions {
            davidson_capacity_per_band: multiplier,
            cutoff_hartree: 2.5,
            bands: Some(3),
            ..Default::default()
        },
        &[KPoint {
            fractional: [0.; 3],
            weight: 1.0,
        }],
        OccupationPolicy::Insulator {
            minimum_gap_hartree: 1e-7,
        },
    )
}

#[test]
fn capacity_range_and_default() {
    assert_eq!(ScfOptions::default().davidson_capacity_per_band, 6);
    for value in [0, 1, 13, usize::MAX] {
        let error = problem(value)
            .err()
            .expect("invalid capacity must reject")
            .to_string();
        assert!(error.contains("minimum 2"), "{error}");
        assert!(!bands_fit_subspace(48, value, 20, 2969));
    }
    for value in 2..=12 {
        assert!(problem(value).is_ok());
        assert!(bands_fit_subspace(48, value, 68, 11877));
        assert!(bands_fit_subspace(48, value, 50, 50));
        assert!(!bands_fit_subspace(48, value, 100, 99));
    }
}

#[test]
fn exact_restart_rejects_changed_capacity() {
    use dft_engine::checkpoint::*;
    let mut session = ScfSession::new(
        problem(6).unwrap(),
        CompletionMode::Strict,
        ExecutionConfig::default(),
    )
    .unwrap();
    session.initialize().unwrap();
    let (model, execution) = fixture::identity(0, 1);
    let cp = session
        .checkpoint(
            CheckpointKind::ExactContinuation,
            model.clone(),
            execution.clone(),
        )
        .unwrap();
    assert_eq!(
        ScfSession::restore_exact(problem(2).unwrap(), &cp, &model, &execution).err(),
        Some(Rejection::Options)
    );
    assert!(ScfSession::restore_exact(problem(6).unwrap(), &cp, &model, &execution).is_ok());
}
