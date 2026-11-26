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

/// CELL-FROZEN-DENSITY regression (dev08, real Si 30 Ry and MgO 50 Ry,
/// variable-cell relaxation from a compressed start): the default reference
/// path (`refscf`) built its density plan from the current cell's 4 x cutoff
/// sphere instead of the frozen reference Miller set, so after the cell
/// expanded every SCF failed with "density sphere reaches FFT Nyquist". The
/// grid here is the tightest canonical grid of the reference density sphere
/// and the cell is expanded isotropically by 4 %; the fixture asserts that the
/// expanded cell's own sphere does not fit, so it keeps exercising the bug.
/// Both the Gamma driver and the complex k driver are checked against the
/// strict oracle on the same frozen protocol.
#[test]
fn reference_path_in_an_expanded_cell_on_a_tight_frozen_grid() {
    use dft_engine::refscf::{api, profile::Profile};
    use dft_geometry::kpoints::KPoint;
    // a * |G|max / (2 pi) = 4.90 in the reference cell and 5.09 after the
    // 4 % expansion: the expanded sphere crosses a Miller shell.
    let reference =
        Cell::from_user_rows_bohr([[8.88, 0., 0.], [0., 8.88, 0.], [0., 0., 8.88]]).unwrap();
    let expanded =
        Cell::from_user_rows_bohr([[9.2352, 0., 0.], [0., 9.2352, 0.], [0., 0., 9.2352]]).unwrap();
    let pp = common::potential(1., 1.5, false, 401);
    let refs = [SpeciesPotential {
        label: "X",
        potential: &pp,
    }];
    let atoms = [Atom::new("X", [0.; 3]), Atom::new("X", [0.27, 0.23, 0.19])];
    let cutoff = 1.5;
    let dense = dft_units::Hartree(4. * cutoff);
    let max_miller = |cell: &Cell| {
        let sphere = dft_basis::PlaneWaveBasis::new(cell, dense).unwrap();
        std::array::from_fn::<usize, 3, _>(|d| {
            sphere
                .vectors()
                .iter()
                .map(|g| g.miller[d].unsigned_abs() as usize)
                .max()
                .unwrap()
        })
    };
    let grid = max_miller(&reference).map(|m| 2 * m + 1);
    let own = max_miller(&expanded);
    assert!(
        (0..3).any(|d| 2 * own[d] >= grid[d]),
        "fixture must put the expanded cell's own sphere beyond the frozen grid: max Miller {own:?}, grid {grid:?}"
    );
    let frozen = |base: ScfOptions| ScfOptions {
        cutoff_hartree: cutoff,
        bands: Some(4),
        grid_override: Some(grid),
        strain_reference: Some(Box::new(reference.clone())),
        ..base
    };
    let strict = ScfOptions {
        energy_tolerance: 1e-12,
        density_tolerance: 1e-10,
        eigen_tolerance: 1e-10,
        ..frozen(ScfOptions::strict_reference())
    };
    let options = frozen(ScfOptions::default());
    let input = || dft_engine::refstart::StartInput {
        cell: &expanded,
        atoms: &atoms,
        potentials: &refs,
        options: &options,
    };
    let profile = Profile {
        threshold_ry: 1e-12,
        ..Profile::default()
    };
    for k in [[0.; 3], [0.25, 0., 0.]] {
        let points = vec![KPoint {
            fractional: k,
            weight: 1.,
        }];
        let run = api::solve(input(), api::Request::new(profile, points.clone()))
            .expect("reference path SCF on the frozen grid")
            .run
            .into_scalar()
            .unwrap();
        assert!(run.outcome.converged);
        let mut oracle = strict.clone();
        if k != [0.; 3] {
            oracle.gamma_real = dft_engine::GammaReal::Off;
        }
        let strict_result = ScfProblem::with_kpoints(
            &expanded,
            &atoms,
            &refs,
            oracle,
            &points,
            dft_engine::OccupationPolicy::Insulator {
                minimum_gap_hartree: 1e-6,
            },
        )
        .unwrap()
        .run()
        .unwrap();
        let de = run.outcome.result.total_hartree() - strict_result.total_hartree();
        println!(
            "CELL-FROZEN-DENSITY k={k:?}: grid {grid:?} (expanded max Miller {own:?}), dE reference-strict {de:e} Ha"
        );
        assert!(de.abs() <= 1e-10, "{de:e}");
    }
}
