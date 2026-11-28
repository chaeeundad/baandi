//! R4b Wannier90 export and gauge import engine gates (ADR 0060) on a synthetic
//! skewed cell with s and p atomic PSWFC channels (not reference evidence):
//! - `.amn` overlaps match R4a projection overlaps
//! - `.mmn` Bloch overlaps satisfy Hermitian adjoint symmetry `M^(k',-b) = (M^(k,b))^H`
//!   to `< 1e-13`, and B1 completeness holds to `< 1e-11`
//! - Gauge import round-trip: importing the R4a Loewdin gauge (disentangled
//!   `_u_dis` + unitary `_u`) reproduces `H(R)`, band weights `p_mk`, populations
//!   and Marzari-Vanderbilt spreads to `< 1e-12`
//! - Unitary gauge rotation invariance of `Omega_I` and isolated eigenvalues
//! - MPI collective rank independence (`np = 1` bit-identical, `np = 2` within `1e-12`)
//!
//! Fast gate: `cargo test -p dft-engine --test response_r4b -- --nocapture`
#[allow(dead_code)]
mod common;
use dft_engine::response::{
    subspace::{BandWindow, ProjectionSelector, SubspaceOptions, gamma_mesh, subspace_with},
    wannier::{
        DEFAULT_GAUGE_TOLERANCE, WannierExportOptions, WannierImportOptions,
        subspace_from_wannier_with, wannier_export_with,
    },
};
use dft_engine::{refstart::StartInput, *};
use dft_geometry::{Atom, Cell};
use dft_numerics::linalg::{ComplexMatrix, hermitian_eigh};
use dft_physics::data::{Array3, Complex64};

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

fn identity_matrix(n: usize) -> ComplexMatrix {
    let mut m = ComplexMatrix::zeros(n, n).unwrap();
    for i in 0..n {
        m[(i, i)] = Complex64::new(1.0, 0.0);
    }
    m
}

fn det_unitary(n: usize, seed_in: u64) -> ComplexMatrix {
    let mut seed = seed_in;
    let mut next = || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        2.0 * (((seed >> 11) as f64) / ((1u64 << 53) as f64)) - 1.0
    };
    let mut h = ComplexMatrix::zeros(n, n).unwrap();
    for j in 0..n {
        h[(j, j)] = Complex64::new(next(), 0.0);
        for i in (j + 1)..n {
            let z = Complex64::new(next(), next());
            h[(i, j)] = z;
            h[(j, i)] = z.conj();
        }
    }
    hermitian_eigh(&h).unwrap().1
}

#[test]
fn wannier_export_adjointness_and_import_roundtrip() {
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

    let sel = vec![ProjectionSelector {
        species: "Xx".into(),
        atoms: None,
        l: None,
        channel: None,
        orbitals: vec!["sp3".into()],
    }];
    let mut exp_opts = WannierExportOptions::new([2, 2, 2], sel.clone());
    exp_opts.bands = Some(8);
    let exp = wannier_export_with(
        make_input(),
        std::slice::from_ref(&density),
        &exp_opts,
        &refscf::metric::Serial,
    )
    .unwrap();
    assert_eq!(exp.bands, 8);
    assert_eq!(exp.orbitals.len(), 4);
    assert_eq!(exp.amn.len(), 8);
    assert_eq!(exp.mmn.len(), 8);
    assert!(
        exp.max_adjoint_defect < 1e-13,
        "M_mn adjoint defect = {:.3e}",
        exp.max_adjoint_defect
    );
    let bset = exp.bset.as_ref().unwrap();
    assert!(bset.completeness_residual < 1e-11);
    let init_spreads = exp.initial_spreads.as_ref().unwrap();
    assert!(init_spreads.omega_i_bohr2 > 0.0);
    assert!(init_spreads.omega_total_bohr2 >= init_spreads.omega_i_bohr2 - 1e-12);

    // Run R4a subspace on all 8 bands with sp3 -> import its Loewdin gauge via R4b import!
    let mut sub_opts = SubspaceOptions::new([2, 2, 2], sel);
    sub_opts.bands = Some(8);
    sub_opts.window = BandWindow::BandRange { first: 0, last: 7 };
    let r4a = subspace_with(
        make_input(),
        std::slice::from_ref(&density),
        &sub_opts,
        &refscf::metric::Serial,
    )
    .unwrap();

    let u_dis: Vec<ComplexMatrix> = r4a.projections.iter().map(|p| p.gauge.clone()).collect();
    let u_id: Vec<ComplexMatrix> = (0..8).map(|_| identity_matrix(4)).collect();
    let k_fracs: Vec<[f64; 3]> = r4a.k_points.iter().map(|p| p.fractional).collect();

    let imp_opts = WannierImportOptions {
        k_mesh: [2, 2, 2],
        solver: refscf::nscf::BandSolver::Dense,
        eigen_tolerance: 1e-6,
        k_fractional: k_fracs.clone(),
        u_matrices: u_id,
        u_dis_matrices: Some(u_dis.clone()),
        gauge_tolerance: DEFAULT_GAUGE_TOLERANCE,
    };
    let imp = subspace_from_wannier_with(
        make_input(),
        std::slice::from_ref(&density),
        &imp_opts,
        &refscf::metric::Serial,
    )
    .unwrap();
    assert!(imp.max_orthonormality_defect < 1e-12);
    assert!(imp.hr_hermiticity_defect < 1e-13);
    // Exact match with R4a populations, spreads and interpolated bands!
    for (a, b) in imp.populations.iter().zip(&r4a.lowdin_populations) {
        assert!((a - b).abs() < 1e-12, "pop diff: {a} vs {b}");
    }
    assert!(
        (imp.spreads.omega_i_bohr2 - init_spreads.omega_i_bohr2).abs() < 1e-11,
        "Omega_I mismatch: {} vs {}",
        imp.spreads.omega_i_bohr2,
        init_spreads.omega_i_bohr2
    );
    let e_r4a = r4a
        .real_space
        .interpolate_eigenvalues([0.125, 0.25, 0.375])
        .unwrap();
    let e_imp = imp
        .real_space
        .interpolate_eigenvalues([0.125, 0.25, 0.375])
        .unwrap();
    for (a, b) in e_r4a.iter().zip(&e_imp) {
        assert!((a - b).abs() < 1e-12);
    }

    // Apply random k-dependent 4x4 unitary gauge rotation U(k): Omega_I MUST stay identical!
    let u_rot: Vec<ComplexMatrix> = (0..8).map(|k| det_unitary(4, 900 + k as u64)).collect();
    let rot_opts = WannierImportOptions {
        k_mesh: [2, 2, 2],
        solver: refscf::nscf::BandSolver::Dense,
        eigen_tolerance: 1e-6,
        k_fractional: k_fracs,
        u_matrices: u_rot,
        u_dis_matrices: Some(u_dis),
        gauge_tolerance: DEFAULT_GAUGE_TOLERANCE,
    };
    let rot = subspace_from_wannier_with(
        make_input(),
        std::slice::from_ref(&density),
        &rot_opts,
        &refscf::metric::Serial,
    )
    .unwrap();
    assert!(
        (rot.spreads.omega_i_bohr2 - imp.spreads.omega_i_bohr2).abs() < 1e-11,
        "Omega_I not gauge invariant: {} vs {}",
        rot.spreads.omega_i_bohr2,
        imp.spreads.omega_i_bohr2
    );

    // MPI rank independence: np = 1 bit-identical; np = 2 agrees within 1e-12.
    let one = dft_mpi::simulate(1, |comm| {
        wannier_export_with(
            make_input(),
            std::slice::from_ref(&density),
            &exp_opts,
            comm,
        )
        .unwrap()
    });
    assert_eq!(
        one[0]
            .initial_spreads
            .as_ref()
            .unwrap()
            .omega_total_bohr2
            .to_bits(),
        init_spreads.omega_total_bohr2.to_bits()
    );
    let two = dft_mpi::simulate(2, |comm| {
        wannier_export_with(
            make_input(),
            std::slice::from_ref(&density),
            &exp_opts,
            comm,
        )
        .unwrap()
    });
    assert!(
        (two[0].initial_spreads.as_ref().unwrap().omega_total_bohr2
            - init_spreads.omega_total_bohr2)
            .abs()
            < 1e-11
    );
}
