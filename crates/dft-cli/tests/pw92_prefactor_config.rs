#[allow(dead_code)]
#[path = "../../dft-engine/tests/common/checkpoint.rs"]
mod checkpoint_fixture;
#[allow(dead_code)]
#[path = "common/mod.rs"]
mod cli_common;
#[allow(dead_code)]
#[path = "../../dft-engine/tests/common/mod.rs"]
mod common;
use dft_engine::{Functional, Pw92Prefactor as A, ScfOptions, ScfProblem, SpeciesPotential};
use dft_geometry::{Atom, Cell};
use std::{fs, process::Command};

#[test]
fn config_rejects_unknown_and_unsupported_conventions() {
    let dir = std::env::temp_dir().join(format!("p27b-config-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("syn.xml"),
        cli_common::gaussian_with(2., 1., false, 401)
            .replace("SLA PZ", "PBE")
            .replace("element=\"Xx\"", "element=\"He\""),
    )
    .unwrap();
    let base = "[system]\nname='prefactor'\nlattice_unit='bohr'\nlattice_vectors=[[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0]]\n[[system.atoms]]\nspecies='He'\nfractional=[0.0,0.0,0.0]\n[pseudopotentials]\nHe='syn.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='pbe'\n";
    for (key, functional, spin, ok, message) in [
        ("", "pbe", "", true, ""),
        ("modified", "pbe", "", true, ""),
        ("original", "pbe", "", true, ""),
        ("other", "pbe", "", false, "modified or original"),
        ("original", "lda-pw", "", false, "unpolarized PBE"),
        ("original", "lda-pz", "", false, "unpolarized PBE"),
        (
            "original",
            "pbe",
            "[spin]\nmode='collinear'\ntotal_moment=0.0\n",
            false,
            "unpolarized PBE",
        ),
    ] {
        let mut text = base.replace("functional='pbe'", &format!("functional='{functional}'"));
        if !key.is_empty() {
            text += &format!("pw92_prefactor='{key}'\n");
        }
        text += spin;
        fs::write(dir.join("in.toml"), text).unwrap();
        let r = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .current_dir(&dir)
            .args(["validate", "in.toml"])
            .output()
            .unwrap();
        let msg = String::from_utf8_lossy(&r.stderr);
        assert_eq!(r.status.success(), ok, "key={key} {msg}");
        if !ok {
            assert!(msg.contains(message), "{msg}");
        }
    }
}

fn small_problem(a: A) -> ScfProblem {
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
            pw92_prefactor: a,
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
    for (source, target) in [(A::Original, A::Modified), (A::Modified, A::Original)] {
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
                calc_id: "pw92-warm-run".into(),
                ..Default::default()
            },
        )
        .unwrap();
    }
}

#[test]
#[ignore = "release force and strain SCF oracle: cargo test --release -p dft-cli --test pw92_prefactor_config original_force_stress -- --ignored --nocapture"]
fn original_force_stress() {
    use dft_engine::{
        ElectronicState, ForceState,
        forces::forces,
        stress::{energy_at, frozen_stress},
    };
    let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [0., 8., 0.], [0., 0., 8.]]).unwrap();
    let p = common::potential_with_functional(1., 1.5, true, 801, "PBE");
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &p,
    }];
    let atoms = vec![
        Atom::new("Xx", [0.; 3]),
        Atom::new("Xx", [0.27, 0.23, 0.19]),
    ];
    let options = ScfOptions {
        functional: Functional::Pbe,
        pw92_prefactor: A::Original,
        cutoff_hartree: 1.,
        grid_override: Some([32; 3]),
        density_tolerance: 1e-11,
        energy_tolerance: 1e-13,
        eigen_tolerance: 1e-12,
        max_iterations: 500,
        ..Default::default()
    };
    let run = |atoms: &[Atom]| {
        ScfProblem::new(&cell, atoms, &species, options.clone())
            .unwrap()
            .run()
            .unwrap()
    };
    let result = run(&atoms);
    let state = ForceState {
        cell: &cell,
        atoms: &atoms,
        potentials: &species,
        options: &options,
        electronic: ElectronicState::Fixed(&result),
        density_floor: None,
        symmetry: &[],
    };
    let force = forces(&state).unwrap();
    let stress = frozen_stress(&state).unwrap();
    for axis in 0..3 {
        let mut errors = [0.; 3];
        for (i, h) in [0.01, 0.005, 0.0025].into_iter().enumerate() {
            let mut plus = atoms.clone();
            let mut minus = atoms.clone();
            plus[1].fractional_position[axis] += h / 8.;
            minus[1].fractional_position[axis] -= h / 8.;
            let fd = -(run(&plus).total_hartree() - run(&minus).total_hartree()) / (2. * h);
            errors[i] = (fd - force.values[1][axis]).abs();
        }
        println!("force axis={axis} errors={errors:?} Ha/Bohr");
        assert!(errors[2] < 1e-6);
        assert!(errors[2] < 0.15 * errors[0] || errors.iter().all(|e| *e < 1e-9));
    }
    for a in 0..3 {
        for b in a..3 {
            let mut errors = [0.; 3];
            for (i, h) in [0.002, 0.001, 0.0005].into_iter().enumerate() {
                let strained = |h: f64| {
                    let c = cell.lattice_columns_bohr();
                    let mut d = c;
                    for j in 0..3 {
                        d[a][j] += h * c[b][j];
                        if a != b {
                            d[b][j] += h * c[a][j];
                        }
                    }
                    Cell::from_user_rows_bohr(std::array::from_fn(|i| {
                        std::array::from_fn(|j| d[j][i])
                    }))
                    .unwrap()
                };
                let mut frozen = options.clone();
                frozen.strain_reference = Some(Box::new(cell.clone()));
                let fd = (energy_at(&state, &strained(h), frozen.clone()).unwrap()
                    - energy_at(&state, &strained(-h), frozen).unwrap())
                    / (2. * h * cell.volume_bohr3() * if a == b { 1. } else { 2. });
                errors[i] = (fd - stress.values[a][b]).abs();
            }
            println!("stress {a}{b} errors={errors:?} Ha/Bohr^3");
            assert!(errors[2] < 1e-7);
            assert!(errors[2] < 0.15 * errors[0] || errors.iter().all(|e| *e < 1e-9));
        }
    }
}
