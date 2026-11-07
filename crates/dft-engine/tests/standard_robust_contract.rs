#[allow(dead_code)]
mod common;
#[allow(dead_code)]
#[path = "common/checkpoint.rs"]
mod fixture;
use dft_engine::{checkpoint::*, *};
use dft_geometry::{Atom, Cell, kpoints::KPoint};

fn problem(standard: bool, density_forcing: bool) -> Result<ScfProblem, ScfError> {
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
            standard: standard.then_some(StandardProfile {
                threshold_ha: 1e-14,
            }),
            standard_density_forcing: density_forcing,
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
fn density_forcing_is_default_off_and_cannot_change_strict_profile() {
    assert!(!ScfOptions::default().standard_density_forcing);
    assert!(matches!(
        problem(false, true),
        Err(ScfError::InvalidInput(_))
    ));
}

#[test]
fn exact_continuation_preserves_density_forcing_and_rejects_switch_mismatch() {
    let mut session = ScfSession::new(
        problem(true, true).unwrap(),
        CompletionMode::Strict,
        ExecutionConfig::default(),
    )
    .unwrap();
    session.initialize().unwrap();
    for _ in 0..6 {
        assert_eq!(session.step().unwrap().status, StepStatus::Ready);
    }
    let (model, execution) = fixture::identity(0, 1);
    let cp = session
        .checkpoint(
            CheckpointKind::ExactContinuation,
            model.clone(),
            execution.clone(),
        )
        .unwrap();
    assert_eq!(
        ScfSession::restore_exact(problem(true, false).unwrap(), &cp, &model, &execution).err(),
        Some(Rejection::Options)
    );
    let restored =
        ScfSession::restore_exact(problem(true, true).unwrap(), &cp, &model, &execution).unwrap();
    let direct = fixture::finish(session);
    let resumed = fixture::finish(restored);
    assert!(!direct.certified && !resumed.certified);
    assert!(direct.standard_density_forcing && resumed.standard_density_forcing);
    fixture::assert_identical(&direct, &resumed);
}
