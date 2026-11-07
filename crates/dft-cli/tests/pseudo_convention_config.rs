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
    Functional, LocalMesh as L, RadialEvenRule as R, ScfOptions, ScfProblem, SpeciesPotential,
};
use dft_geometry::{Atom, Cell};
use std::{fs, process::Command};

#[test]
fn config_rejects_unknown_conventions() {
    let dir = std::env::temp_dir().join(format!("p27d-config-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("syn.xml"),
        cli_common::gaussian_with(2., 1., false, 401)
            .replace("SLA PZ", "PBE")
            .replace("element=\"Xx\"", "element=\"He\""),
    )
    .unwrap();
    let base = "[system]\nname='radial'\nlattice_unit='bohr'\nlattice_vectors=[[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0]]\n[[system.atoms]]\nspecies='He'\nfractional=[0.0,0.0,0.0]\n[pseudopotentials]\nHe='syn.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='pbe'\n";
    for (extra, ok) in [
        ("", true),
        (
            "[pseudo]\nlocal_mesh='baandi'\nradial_even_rule='simpson_3_8'",
            true,
        ),
        (
            "[pseudo]\nlocal_mesh='qe_like'\nradial_even_rule='simpson_quadratic_end'",
            true,
        ),
        ("[pseudo]\nlocal_mesh='invalid'", false),
        ("[pseudo]\nradial_even_rule='invalid'", false),
        ("[pseudo]\nunknown='baandi'", false),
    ] {
        fs::write(dir.join("in.toml"), format!("{base}{extra}")).unwrap();
        let r = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .current_dir(&dir)
            .args(["validate", "in.toml"])
            .output()
            .unwrap();
        assert_eq!(
            r.status.success(),
            ok,
            "{extra}: {}",
            String::from_utf8_lossy(&r.stderr)
        );
    }
}

fn small_problem(a: (L, R)) -> ScfProblem {
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
            local_mesh: a.0,
            radial_even_rule: a.1,
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
        ((L::Baandi, R::Simpson3_8), (L::QeLike, R::Simpson3_8)),
        ((L::QeLike, R::Simpson3_8), (L::Baandi, R::Simpson3_8)),
        (
            (L::Baandi, R::Simpson3_8),
            (L::Baandi, R::SimpsonQuadraticEnd),
        ),
        (
            (L::Baandi, R::SimpsonQuadraticEnd),
            (L::Baandi, R::Simpson3_8),
        ),
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
