//! R4a atomic-orbital subspace projection engine gates (ADR 0059) on a
//! synthetic in-memory pseudopotential with s and p PSWFC channels in a skewed
//! cell (not reference evidence). Checks:
//! - isolated s-band projection: orthonormality <= 1e-12, exact eigenvalue
//!   recovery on every mesh k-point, H(R) Hermiticity, exact interpolation
//! - entangled sp3 projection: sum rule `sum_m p_mk = N_orb`, orthonormality
//! - MPI collective rank independence (`np = 1` bit-identical to serial,
//!   `np = 2` within 1e-12 Ha)
//! - input guards (missing species, spin-polarized density, narrow window)
//!
//! Fast gate: `cargo test -p dft-engine --test response_r4a -- --nocapture`
#[allow(dead_code)]
mod common;
use dft_engine::response::subspace::{
    BandWindow, ProjectionSelector, SubspaceOptions, gamma_mesh, subspace_with,
};
use dft_engine::{refstart::StartInput, *};
use dft_geometry::{Atom, Cell};
use dft_physics::data::Array3;

fn build_system() -> (
    Cell,
    Vec<Atom>,
    dft_upf::Pseudopotential,
    ScfOptions,
    Array3<f64>,
) {
    let cell = Cell::from_user_rows_bohr([[6.0, 0., 0.], [0.7, 5.5, 0.], [0.3, 0.4, 6.5]]).unwrap();
    let p = common::potential_with_atomic_start_dij(
        2.,
        0.8,
        true,
        401,
        "SLA PZ",
        "-1 0.1 0 0.1 -0.5 0 0 0 0.8",
    );
    let atoms = vec![Atom::new("Xx", [0.13, 0.21, 0.37])];
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &p,
    }];
    let options = ScfOptions {
        cutoff_hartree: 4.0,
        bands: Some(8),
        eigen_tolerance: 1e-8,
        functional: Functional::LdaPz,
        ..ScfOptions::strict_reference()
    };
    let points = gamma_mesh([2, 2, 2]).unwrap();
    let problem = ScfProblem::with_kpoints(
        &cell,
        &atoms,
        &potentials,
        options.clone(),
        &points,
        OccupationPolicy::Insulator {
            minimum_gap_hartree: 1e-7,
        },
    )
    .unwrap();
    let density = problem.starting_density();
    (cell, atoms, p, options, density)
}

#[test]
fn isolated_s_projection_recovers_band_and_interpolates_exactly() {
    let (cell, atoms, p, scf_opts, density) = build_system();
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &p,
    }];
    let input = StartInput {
        cell: &cell,
        atoms: &atoms,
        potentials: &potentials,
        options: &scf_opts,
    };
    let mut opts = SubspaceOptions::new(
        [2, 2, 2],
        vec![ProjectionSelector {
            species: "Xx".into(),
            atoms: None,
            l: Some(0),
            channel: None,
            orbitals: vec![],
        }],
    );
    opts.bands = Some(6);
    opts.window = BandWindow::BandRange { first: 0, last: 0 };

    let res = subspace_with(
        input,
        std::slice::from_ref(&density),
        &opts,
        &refscf::metric::Serial,
    )
    .unwrap();
    assert_eq!(res.orbitals.len(), 1);
    assert_eq!(res.orbitals[0].orbital_name, "s");
    assert!(
        res.max_orthonormality_defect < 1e-13,
        "orth defect = {:.3e}",
        res.max_orthonormality_defect
    );
    assert!(
        res.hr_hermiticity_defect < 1e-14,
        "H(R) Hermiticity = {:.3e}",
        res.hr_hermiticity_defect
    );
    // Isolated 1-band manifold: interpolated eigenvalue at every mesh k matches
    // the exact KS band 0 eigenvalue to machine precision, and population = 2.0.
    assert!((res.lowdin_populations[0] - 2.0).abs() < 1e-13);
    for (k, pt) in res.k_points.iter().enumerate() {
        let interp = res
            .real_space
            .interpolate_eigenvalues(pt.fractional)
            .unwrap();
        assert_eq!(interp.len(), 1);
        assert!(
            (interp[0] - res.eigenvalues[k][0]).abs() < 1e-13,
            "k={k}: interp={} ks={}",
            interp[0],
            res.eigenvalues[k][0]
        );
        assert!((res.projections[k].band_weights[0] - 1.0).abs() < 1e-13);
        assert_eq!(res.projections[k].band_weights[1], 0.0);
    }
}

#[test]
fn entangled_sp3_projection_sum_rules_and_rank_independence() {
    let (cell, atoms, p, scf_opts, density) = build_system();
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &p,
    }];
    let make_input = || StartInput {
        cell: &cell,
        atoms: &atoms,
        potentials: &potentials,
        options: &scf_opts,
    };
    let mut opts = SubspaceOptions::new(
        [2, 2, 2],
        vec![ProjectionSelector {
            species: "Xx".into(),
            atoms: None,
            l: None,
            channel: None,
            orbitals: vec!["sp3".into()],
        }],
    );
    opts.bands = Some(8);
    opts.window = BandWindow::BandRange { first: 0, last: 7 };

    let serial = subspace_with(
        make_input(),
        std::slice::from_ref(&density),
        &opts,
        &refscf::metric::Serial,
    )
    .unwrap();
    assert_eq!(serial.orbitals.len(), 4);
    assert!(serial.max_orthonormality_defect < 1e-12);
    assert!(serial.hr_hermiticity_defect < 1e-13);
    for (k, proj) in serial.projections.iter().enumerate() {
        let sum_w: f64 = proj.band_weights.iter().sum();
        assert!(
            (sum_w - 4.0).abs() < 1e-12,
            "k={k}: sum_m p_mk = {sum_w}, expected 4.0"
        );
    }

    // MPI collective test: np = 1 bit-identical; np = 2 agrees within 1e-12 Ha.
    let one = dft_mpi::simulate(1, |comm| {
        subspace_with(make_input(), std::slice::from_ref(&density), &opts, comm).unwrap()
    });
    assert_eq!(
        one[0]
            .lowdin_populations
            .iter()
            .map(|x| x.to_bits())
            .collect::<Vec<_>>(),
        serial
            .lowdin_populations
            .iter()
            .map(|x| x.to_bits())
            .collect::<Vec<_>>()
    );
    let two = dft_mpi::simulate(2, |comm| {
        subspace_with(make_input(), std::slice::from_ref(&density), &opts, comm).unwrap()
    });
    for (a, b) in two[0]
        .lowdin_populations
        .iter()
        .zip(&serial.lowdin_populations)
    {
        assert!((a - b).abs() < 1e-12, "np=2 population diff: {a} vs {b}");
    }
    let e_ser = serial
        .real_space
        .interpolate_eigenvalues([0.25, 0.25, 0.25])
        .unwrap();
    let e_two = two[0]
        .real_space
        .interpolate_eigenvalues([0.25, 0.25, 0.25])
        .unwrap();
    for (a, b) in e_ser.iter().zip(&e_two) {
        assert!(
            (a - b).abs() < 1e-12,
            "np=2 interpolated band diff: {a} vs {b}"
        );
    }
}

#[test]
fn subspace_input_guards_reject() {
    let (cell, atoms, p, scf_opts, density) = build_system();
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &p,
    }];
    let make_input = || StartInput {
        cell: &cell,
        atoms: &atoms,
        potentials: &potentials,
        options: &scf_opts,
    };

    // Spin-polarized density rejected.
    let opts = SubspaceOptions::new(
        [2, 2, 2],
        vec![ProjectionSelector {
            species: "Xx".into(),
            atoms: None,
            l: Some(0),
            channel: None,
            orbitals: vec![],
        }],
    );
    assert!(
        subspace_with(
            make_input(),
            &[density.clone(), density.clone()],
            &opts,
            &refscf::metric::Serial
        )
        .is_err()
    );

    // Window narrower than orbital count rejected.
    let mut narrow = SubspaceOptions::new(
        [2, 2, 2],
        vec![ProjectionSelector {
            species: "Xx".into(),
            atoms: None,
            l: Some(1),
            channel: None,
            orbitals: vec![],
        }],
    );
    narrow.bands = Some(6);
    narrow.window = BandWindow::BandRange { first: 0, last: 1 }; // 2 bands < 3 p orbitals
    assert!(
        subspace_with(
            make_input(),
            std::slice::from_ref(&density),
            &narrow,
            &refscf::metric::Serial
        )
        .is_err()
    );
}
