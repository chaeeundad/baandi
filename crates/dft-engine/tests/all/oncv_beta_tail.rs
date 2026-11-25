//! External SG15 smoke, not a QE reference or cutoff-convergence claim.
//! DFT_ONCV_UPF=/absolute/Si_ONCV_PBE-1.0.upf cargo test --release
//! -p dft-engine --test all oncv_beta_tail -- --ignored --nocapture

use dft_engine::{
    Functional, GammaReal, ScfOptions, ScfProblem, SpeciesPotential,
    refscf::{driver, profile::Profile},
    refstart::StartInput,
};
use dft_geometry::{Atom, Cell};

#[test]
#[ignore = "external ONCV file and two Si2 SCFs; run release command above"]
fn si2_oncv_standard_against_strict_oracle() {
    let path = std::env::var_os("DFT_ONCV_UPF")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| {
            std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("../../../pseudo/Si_ONCV_PBE-1.0.upf")
        });
    let pp = dft_upf::from_file(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    assert_eq!(pp.element().trim(), "Si");
    assert_eq!(pp.z_valence(), 4.0);
    let cell =
        Cell::from_user_rows_bohr([[0., 5.13, 5.13], [5.13, 0., 5.13], [5.13, 5.13, 0.]]).unwrap();
    let atoms = [Atom::new("Si", [0.; 3]), Atom::new("Si", [0.25; 3])];
    // Both drivers receive exactly the same immutable parsed projector prefixes.
    let species = [SpeciesPotential {
        label: "Si",
        potential: &pp,
    }];
    let options = ScfOptions {
        functional: Functional::Pbe,
        gamma_real: GammaReal::On,
        cutoff_hartree: 5.0,
        bands: Some(6),
        max_iterations: 200,
        ..Default::default()
    };
    let strict = ScfOptions {
        functional: options.functional,
        gamma_real: options.gamma_real,
        cutoff_hartree: options.cutoff_hartree,
        bands: options.bands,
        max_iterations: options.max_iterations,
        eigen_tolerance: 1e-11,
        energy_tolerance: 1e-12,
        density_tolerance: 1e-9,
        ..ScfOptions::strict_reference()
    };
    let oracle = ScfProblem::new(&cell, &atoms, &species, strict)
        .unwrap()
        .run()
        .unwrap();
    let run = driver::run(
        StartInput {
            cell: &cell,
            atoms: &atoms,
            potentials: &species,
            options: &options,
        },
        Profile::default(),
    )
    .unwrap();
    let result = &run.outcome.result;
    assert!(
        run.outcome.converged && run.reference_stopped,
        "{}",
        run.exit_reason
    );
    assert!(result.verification.as_ref().unwrap().passed());
    let delta = result.total_hartree() - oracle.total_hartree();
    println!(
        "ONCV Si2 10 Ry: reference={:.17e} Ha strict={:.17e} Ha delta={delta:.17e} Ha",
        result.total_hartree(),
        oracle.total_hartree()
    );
    assert!(delta.abs() < 5e-9, "standard class: {delta:e} Ha");
    assert!((result.charge - 8.0).abs() < 1e-10);
}
