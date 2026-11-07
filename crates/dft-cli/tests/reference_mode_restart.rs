#[allow(dead_code)]
#[path = "../../dft-engine/tests/common/checkpoint.rs"]
mod checkpoint_fixture;
#[allow(dead_code)]
#[path = "common/mod.rs"]
mod cli_common;
#[allow(dead_code)]
#[path = "../../dft-engine/tests/common/mod.rs"]
mod common;
use dft_engine::{
    CompatMode, Functional, LocalMesh as L, RadialEvenRule as R, ScfOptions, ScfProblem,
    SpeciesPotential,
};
use dft_geometry::{Atom, Cell};

fn small_problem(mode: CompatMode) -> ScfProblem {
    let cell = Cell::from_user_rows_bohr([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]).unwrap();
    let p = common::potential_with_functional(2., 1., false, 401, "PBE");
    ScfProblem::new(
        &cell,
        &[Atom::new("Xx", [0.13, 0.21, 0.32])],
        &[SpeciesPotential {
            label: "Xx",
            potential: &p,
        }],
        ScfOptions {
            functional: Functional::Pbe,
            compat_mode: mode,
            pw92_prefactor: dft_engine::Pw92Prefactor::Original,
            local_mesh: L::QeLike,
            radial_even_rule: R::SimpsonQuadraticEnd,
            cutoff_hartree: 1.5,
            ..Default::default()
        },
    )
    .unwrap()
}

#[test]
fn exact_restart_rejects_both_directions_and_warm_start_is_explicit() {
    use dft_engine::{CompletionMode, ExecutionConfig, ScfSession, checkpoint::*};
    let (mut model, exec) = checkpoint_fixture::identity(0, 1);
    model.cutoff_hartree = 1.5;
    model.xc = "pbe".into();
    for (source, target) in [
        (CompatMode::None, CompatMode::Reference),
        (CompatMode::Reference, CompatMode::None),
    ] {
        let mut s = ScfSession::new(
            small_problem(source),
            CompletionMode::Strict,
            ExecutionConfig::default(),
        )
        .unwrap();
        s.initialize().unwrap();
        s.step().unwrap();
        let exact = s
            .checkpoint(
                CheckpointKind::ExactContinuation,
                model.clone(),
                exec.clone(),
            )
            .unwrap();
        assert_eq!(
            ScfSession::restore_exact(small_problem(target), &exact, &model, &exec)
                .err()
                .unwrap(),
            Rejection::Options
        );
        ScfSession::restore_exact(small_problem(source), &exact, &model, &exec).unwrap();
        let warm = s
            .checkpoint(CheckpointKind::WarmStart, model.clone(), exec.clone())
            .unwrap();
        ScfSession::warm_start(
            small_problem(target),
            &warm,
            &model,
            &exec,
            CompletionMode::Strict,
            ExecutionConfig {
                calc_id: "radial-warm-run".into(),
                ..Default::default()
            },
        )
        .unwrap();
    }
}
