//! ONEPATH-2B regression: the standard path's Coulomb stop-rule metric must use
//! the problem's own density list. In a variable-cell relaxation the Miller set
//! and FFT grid are frozen at the reference cell (`strain_reference`,
//! `grid_override`); when the cell expands, the current cell's 4 x cutoff sphere
//! holds larger Miller indices than the frozen grid can represent. The metric
//! used to be built from the current cell alone and failed with "Miller index
//! ... is outside the canonical frequencies" (dev08, real Si8 PBE, mode cell).
//! Synthetic fixture only.
#[allow(dead_code)]
#[path = "../common/mod.rs"]
mod common;
use dft_engine::{ScfOptions, ScfProblem, SpeciesPotential};
use dft_geometry::{Atom, Cell};

#[test]
fn standard_scf_in_an_expanded_cell_on_the_frozen_reference_grid() {
    let reference = Cell::from_user_rows_bohr([[8., 0., 0.], [0., 8., 0.], [0., 0., 8.]]).unwrap();
    // 15 % expansion: the expanded cell's own density sphere needs more Miller
    // indices than the reference grid holds (checked below, so the fixture
    // cannot silently stop exercising the bug).
    let expanded =
        Cell::from_user_rows_bohr([[9.2, 0., 0.], [0., 9.2, 0.], [0., 0., 9.2]]).unwrap();
    let pp = common::potential(1., 1.5, false, 401);
    let refs = [SpeciesPotential {
        label: "X",
        potential: &pp,
    }];
    let atoms = [Atom::new("X", [0.; 3]), Atom::new("X", [0.27, 0.23, 0.19])];
    let standard = ScfOptions {
        cutoff_hartree: 1.5,
        bands: Some(4),
        ..ScfOptions::strict_reference()
    };
    let grid = ScfProblem::new(&reference, &atoms, &refs, standard.clone())
        .unwrap()
        .grid_dimensions();
    let own = ScfProblem::new(&expanded, &atoms, &refs, standard.clone())
        .unwrap()
        .grid_dimensions();
    assert!(
        (0..3).any(|i| own[i] > grid[i]),
        "fixture must need a larger grid in the expanded cell: own {own:?}, frozen {grid:?}"
    );
    let frozen = |base: ScfOptions| ScfOptions {
        grid_override: Some(grid),
        strain_reference: Some(Box::new(reference.clone())),
        ..base
    };
    let standard_result = ScfProblem::new(&expanded, &atoms, &refs, frozen(standard.clone()))
        .expect("standard problem on the frozen grid")
        .run()
        .expect("standard SCF converges with verification");
    assert!(standard_result.certified);
    let strict = ScfOptions {
        energy_tolerance: 1e-12,
        density_tolerance: 1e-10,
        eigen_tolerance: 1e-10,
        cutoff_hartree: 1.5,
        bands: Some(4),
        ..ScfOptions::strict_reference()
    };
    let strict_result = ScfProblem::new(&expanded, &atoms, &refs, frozen(strict))
        .unwrap()
        .run()
        .unwrap();
    let de = standard_result.total_hartree() - strict_result.total_hartree();
    println!(
        "ONEPATH-2B frozen expanded cell: grid {grid:?} (own {own:?}), dE standard-strict {de:e} Ha"
    );
    // Standard class: the delivered energy is resolved to the threshold.
    assert!(de.abs() <= 1e-10, "{de:e}");
}
