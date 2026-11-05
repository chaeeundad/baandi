#[allow(dead_code)]
mod common;
use dft_engine::{
    ElectronicState, ForceState, Functional, ScfOptions, ScfProblem, SpeciesPotential, forces,
};
use dft_geometry::{Atom, Cell};
#[test]
fn retained_model_rejects_independent_input_mutations() {
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0., 5., 0.], [0., 0., 5.]]).unwrap();
    let atoms = vec![
        Atom::new("Xx", [0.; 3]),
        Atom::new("Xx", [0.27, 0.23, 0.19]),
    ];
    let pp = common::potential(1., 1.5, false, 401);
    let refs = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let opts = ScfOptions {
        cutoff_hartree: 1.,
        ..Default::default()
    };
    let mut r = ScfProblem::new(&cell, &atoms, &refs, opts.clone())
        .unwrap()
        .run()
        .unwrap();
    let eval = |cell: &Cell,
                atoms: &[Atom],
                refs: &[SpeciesPotential<'_>],
                opts: &ScfOptions,
                r: &dft_engine::ScfResult| {
        forces(&ForceState {
            cell,
            atoms,
            potentials: refs,
            options: opts,
            electronic: ElectronicState::Fixed(r),
            density_floor: None,
            symmetry: &[],
        })
    };
    eval(&cell, &atoms, &refs, &opts, &r).unwrap();
    let reject = |result: Result<dft_engine::ForceResult, dft_engine::ScfError>| {
        let message = result.unwrap_err().to_string();
        assert!(
            message.contains("observables.forces.model") && message.contains("rerun SCF"),
            "{message}"
        );
    };
    let changed = Cell::from_user_rows_bohr([[5.01, 0., 0.], [0., 5., 0.], [0., 0., 5.]]).unwrap();
    reject(eval(&changed, &atoms, &refs, &opts, &r));
    let mut moved = atoms.clone();
    moved[1].fractional_position[0] += 0.001;
    reject(eval(&cell, &moved, &refs, &opts, &r));
    let mut renamed = atoms.clone();
    renamed[0].species = "other".into();
    reject(eval(&cell, &renamed, &refs, &opts, &r));
    let different = common::potential(1., 1.6, false, 401);
    reject(eval(
        &cell,
        &atoms,
        &[SpeciesPotential {
            label: "Xx",
            potential: &different,
        }],
        &opts,
        &r,
    ));
    // Independently parsed identical content succeeds: no pointer identity shortcut.
    let identical = common::potential(1., 1.5, false, 401);
    eval(
        &cell,
        &atoms,
        &[SpeciesPotential {
            label: "Xx",
            potential: &identical,
        }],
        &opts,
        &r,
    )
    .unwrap();
    reject(eval(
        &cell,
        &atoms,
        &refs,
        &ScfOptions {
            cutoff_hartree: 1.1,
            ..opts.clone()
        },
        &r,
    ));
    reject(eval(
        &cell,
        &atoms,
        &refs,
        &ScfOptions {
            functional: Functional::LdaPw,
            ..opts.clone()
        },
        &r,
    ));
    let old = r.occupations[0];
    r.occupations[0] = 1.9;
    reject(eval(&cell, &atoms, &refs, &opts, &r));
    r.occupations[0] = old;
    r.k_points[0].point.fractional[0] = 0.1;
    reject(eval(&cell, &atoms, &refs, &opts, &r));
}

#[test]
#[ignore = "release-mode FD policy mutation; cargo test --release -p dft-engine --test force_identity -- --include-ignored --nocapture"]
fn finite_temperature_policy_is_bound_to_certificate() {
    use dft_engine::finite_temperature::{FiniteTProblem, OccupationWorkspace};
    use dft_geometry::kpoints::KPoint;
    use dft_physics::occupations::Constraint;
    let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [1., 7., 0.], [0.5, 1., 9.]]).unwrap();
    let atoms = vec![
        Atom::new("Xx", [0.; 3]),
        Atom::new("Xx", [0.27, 0.23, 0.19]),
    ];
    let pp = common::potential_with_functional(1., 1.5, true, 801, "SLA PW");
    let refs = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let opts = ScfOptions {
        cutoff_hartree: 2.,
        functional: Functional::LdaPw,
        grid_override: Some([48; 3]),
        davidson_subspace: 96,
        ..Default::default()
    };
    let points = [KPoint {
        fractional: [0.; 3],
        weight: 1.,
    }];
    let mut r = FiniteTProblem::new(&cell, &atoms, &refs, opts.clone(), &points, 0.05, 48)
        .unwrap()
        .run()
        .unwrap();
    let evaluate = |r: &dft_engine::finite_temperature::FiniteTResult| {
        forces(&ForceState {
            cell: &cell,
            atoms: &atoms,
            potentials: &refs,
            options: &opts,
            electronic: ElectronicState::FiniteTemperature(r),
            density_floor: None,
            symmetry: &[],
        })
    };
    evaluate(&r).unwrap();
    for (theta, constraint, bands) in [
        (0.06, Constraint::Total, 48),
        (0.05, Constraint::FixedMoment(0.), 48),
        (0.05, Constraint::Total, 47),
    ] {
        let original = std::mem::replace(
            &mut r.occupations,
            OccupationWorkspace::new(bands, &points, false, theta, constraint).unwrap(),
        );
        assert!(
            evaluate(&r)
                .unwrap_err()
                .to_string()
                .contains("observables.forces.model")
        );
        r.occupations = original;
    }
    evaluate(&r).unwrap();
}
