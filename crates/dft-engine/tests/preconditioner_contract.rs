#[allow(dead_code)]
mod common;
#[allow(dead_code)]
#[path = "common/checkpoint.rs"]
mod fixture;
use dft_engine::*;
use dft_geometry::{Atom, Cell, kpoints::KPoint};

fn problem(variant: dft_solvers::DavidsonPreconditioner) -> Result<ScfProblem, ScfError> {
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
            preconditioner: variant,
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
fn exact_restart_rejects_changed_preconditioner() {
    use dft_engine::checkpoint::*;
    let mut session = ScfSession::new(
        problem(dft_solvers::DavidsonPreconditioner::Default).unwrap(),
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
        ScfSession::restore_exact(
            problem(dft_solvers::DavidsonPreconditioner::Tpa).unwrap(),
            &cp,
            &model,
            &execution
        )
        .err(),
        Some(Rejection::Options)
    );
    assert!(
        ScfSession::restore_exact(
            problem(dft_solvers::DavidsonPreconditioner::Default).unwrap(),
            &cp,
            &model,
            &execution
        )
        .is_ok()
    );
}
