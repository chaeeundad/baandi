#[allow(dead_code)]
mod common;
use dft_engine::{
    CompletionMode, ExecutionConfig, ScfOptions, ScfProblem, ScfSession, SpeciesPotential,
    StepStatus,
};
use dft_geometry::{Atom, Cell};

#[test]
fn rejects_invalid_density_and_changed_model_before_initialization() {
    let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [0., 8., 0.], [0., 0., 8.]]).unwrap();
    let pp = common::potential(1., 1.5, false, 401);
    let refs = [SpeciesPotential {
        label: "X",
        potential: &pp,
    }];
    let atoms = [
        Atom::new("X", [0., 0., 0.]),
        Atom::new("X", [0.27, 0.23, 0.19]),
    ];
    let opts = ScfOptions {
        cutoff_hartree: 1.5,
        bands: Some(4),
        ..Default::default()
    };
    let source = ScfProblem::new(&cell, &atoms, &refs, opts.clone()).unwrap();
    let mut session = ScfSession::new(
        ScfProblem::new(&cell, &atoms, &refs, opts.clone()).unwrap(),
        CompletionMode::Strict,
        ExecutionConfig::default(),
    )
    .unwrap();
    let density = source.uniform_density();
    let mut bad = density.clone();
    bad[[0, 0, 0]] = f64::NAN;
    assert!(session.initialize_density_transfer(&source, &bad).is_err());
    let bad = &density * 1.000001;
    assert!(session.initialize_density_transfer(&source, &bad).is_err());
    let mut changed = opts.clone();
    changed.cutoff_hartree = 2.;
    let changed = ScfProblem::new(&cell, &atoms, &refs, changed).unwrap();
    assert!(
        session
            .initialize_density_transfer(&changed, &density)
            .is_err()
    );
    let other_cell =
        Cell::from_user_rows_bohr([[8.1, 0., 0.], [0., 8., 0.], [0., 0., 8.]]).unwrap();
    let other_pp = common::potential(1., 1.6, false, 401);
    let other_refs = [SpeciesPotential {
        label: "X",
        potential: &other_pp,
    }];
    let renamed = [Atom::new("Y", [0.; 3]), Atom::new("Y", [0.27, 0.23, 0.19])];
    let renamed_refs = [SpeciesPotential {
        label: "Y",
        potential: &pp,
    }];
    let mut floor = ScfProblem::new(&cell, &atoms, &refs, opts.clone()).unwrap();
    floor.set_density_floor(1e-10).unwrap();
    let alternatives = [
        (
            "cell",
            ScfProblem::new(&other_cell, &atoms, &refs, opts.clone()).unwrap(),
        ),
        (
            "PP content",
            ScfProblem::new(&cell, &atoms, &other_refs, opts.clone()).unwrap(),
        ),
        (
            "species",
            ScfProblem::new(&cell, &renamed, &renamed_refs, opts.clone()).unwrap(),
        ),
        (
            "occupation layout",
            ScfProblem::new(
                &cell,
                &atoms,
                &refs,
                ScfOptions {
                    bands: Some(5),
                    ..opts.clone()
                },
            )
            .unwrap(),
        ),
        (
            "functional",
            ScfProblem::new(
                &cell,
                &atoms,
                &refs,
                ScfOptions {
                    functional: dft_engine::Functional::LdaPw,
                    allow_functional_override: true,
                    ..opts.clone()
                },
            )
            .unwrap(),
        ),
        ("density floor", floor),
        (
            "k mesh",
            ScfProblem::with_kpoints(
                &cell,
                &atoms,
                &refs,
                opts.clone(),
                &[dft_geometry::kpoints::KPoint {
                    fractional: [0.1, 0., 0.],
                    weight: 1.,
                }],
                dft_engine::OccupationPolicy::Insulator {
                    minimum_gap_hartree: 1e-7,
                },
            )
            .unwrap(),
        ),
    ];
    for (label, other) in alternatives {
        assert!(
            session
                .initialize_density_transfer(&other, &density)
                .is_err(),
            "accepted changed {label}"
        );
        assert!(
            session.density().is_none(),
            "mutated before rejecting {label}"
        );
    }
    let mut wrong_layout = density.clone();
    wrong_layout.swap_axes(0, 1);
    assert!(
        session
            .initialize_density_transfer(&source, &wrong_layout)
            .is_err()
    );
    let wrong_shape = density.clone();
    // A separately shaped array, preserving total length, must not be accepted.
    let wrong_shape = wrong_shape
        .into_shape_with_order((1, 1, density.len()))
        .unwrap();
    assert!(
        session
            .initialize_density_transfer(&source, &wrong_shape)
            .is_err()
    );
    session
        .initialize_density_transfer(&source, &density)
        .unwrap();
    assert_eq!(session.density().unwrap(), &density);
    assert!(session.history().is_empty());
    assert!(
        session
            .initialize_density_transfer(&source, &density)
            .is_err()
    );
}

#[test]
#[ignore = "three tight fresh/warm SCF comparisons; run in release"]
fn three_displaced_structures_match_fresh() {
    let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [0., 8., 0.], [0., 0., 8.]]).unwrap();
    let pp = common::potential(1., 1.5, false, 801);
    let refs = [SpeciesPotential {
        label: "X",
        potential: &pp,
    }];
    let opts = ScfOptions {
        cutoff_hartree: 1.5,
        bands: Some(4),
        density_tolerance: 1e-11,
        energy_tolerance: 1e-12,
        eigen_tolerance: 1e-11,
        max_iterations: 300,
        ..Default::default()
    };
    let atoms = [
        Atom::new("X", [0., 0., 0.]),
        Atom::new("X", [0.27, 0.23, 0.19]),
    ];
    let source = ScfProblem::new(&cell, &atoms, &refs, opts.clone()).unwrap();
    let initial = ScfProblem::new(&cell, &atoms, &refs, opts.clone())
        .unwrap()
        .run()
        .unwrap();
    for dx in [-0.01, 0.005, 0.02] {
        let atoms = [
            Atom::new("X", [0., 0., 0.]),
            Atom::new("X", [0.27 + dx, 0.23, 0.19]),
        ];
        let fresh = ScfProblem::new(&cell, &atoms, &refs, opts.clone())
            .unwrap()
            .run()
            .unwrap();
        let problem = ScfProblem::new(&cell, &atoms, &refs, opts.clone()).unwrap();
        let mut session =
            ScfSession::new(problem, CompletionMode::Strict, ExecutionConfig::default()).unwrap();
        session
            .initialize_density_transfer(&source, &initial.density)
            .unwrap();
        while session.step().unwrap().status == StepStatus::Ready {}
        let warm = session.into_result().unwrap();
        assert!(warm.converged);
        let delta = (warm.result.total_hartree() - fresh.total_hartree()).abs();
        println!(
            "dx={dx} fractional delta_E={delta:.17e} Ha fresh={:.17e} warm={:.17e}",
            fresh.total_hartree(),
            warm.result.total_hartree()
        );
        assert!(delta <= 1e-10);
    }
}

#[test]
fn affine_cell_transfer_preserves_charge_and_rejects_changed_support() {
    let a = Cell::from_user_rows_bohr([[8., 0., 0.], [0., 8., 0.], [0., 0., 8.]]).unwrap();
    let b = Cell::from_user_rows_bohr([[8.2, 0., 0.], [0., 7.9, 0.], [0., 0., 8.1]]).unwrap();
    let pp = common::potential(1., 1.5, false, 401);
    let refs = [SpeciesPotential {
        label: "X",
        potential: &pp,
    }];
    let atoms = [Atom::new("X", [0.; 3]), Atom::new("X", [0.27, 0.23, 0.19])];
    let mut opts = ScfOptions {
        cutoff_hartree: 1.5,
        bands: Some(4),
        ..Default::default()
    };
    let original = ScfProblem::new(&a, &atoms, &refs, opts.clone()).unwrap();
    opts.grid_override = Some(original.grid_dimensions());
    opts.strain_reference = Some(Box::new(a.clone()));
    let source = ScfProblem::new(&a, &atoms, &refs, opts.clone()).unwrap();
    let density = source.uniform_density();
    let target = ScfProblem::new(&b, &atoms, &refs, opts.clone()).unwrap();
    let mut session =
        ScfSession::new(target, CompletionMode::Strict, ExecutionConfig::default()).unwrap();
    assert!(
        session
            .initialize_density_transfer(&source, &density)
            .is_err()
    );
    assert!(
        session
            .initialize_cell_density_transfer(&original, &density)
            .is_err()
    );
    assert!(
        session
            .initialize_cell_density_transfer(&source, &(&density * 1.00001))
            .is_err()
    );
    session
        .initialize_cell_density_transfer(&source, &density)
        .unwrap();
    let transferred = session.density().unwrap();
    assert!((transferred.sum() * b.volume_bohr3() / transferred.len() as f64 - 2.).abs() < 1e-12);
    let ratio = a.volume_bohr3() / b.volume_bohr3();
    for (actual, source) in transferred.iter().zip(&density) {
        assert_eq!(actual.to_bits(), (source * ratio).to_bits());
    }
    let target = ScfProblem::new(
        &b,
        &atoms,
        &refs,
        ScfOptions {
            cutoff_hartree: 1.4,
            ..opts
        },
    )
    .unwrap();
    let mut session =
        ScfSession::new(target, CompletionMode::Strict, ExecutionConfig::default()).unwrap();
    assert!(
        session
            .initialize_cell_density_transfer(&source, &density)
            .is_err()
    );
}
