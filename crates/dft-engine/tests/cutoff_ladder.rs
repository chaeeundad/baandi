//! C1 cutoff ladder: exact prolongation, full-cutoff result identity within the
//! strict class, fallback, thread-count independence. Fast gates run in debug.
#[allow(dead_code)]
mod common;
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_engine::cutoff_ladder::{self, CutoffLadderOptions};
use dft_engine::{GammaReal, ScfOptions, ScfProblem, SpeciesPotential};
use dft_geometry::{Atom, Cell};
use dft_numerics::linalg::ComplexMatrix;
use dft_physics::data::{Array3, Complex64};
use dft_units::Hartree;

struct Model {
    cell: Cell,
    atoms: Vec<Atom>,
    species: Vec<(String, dft_upf::Pseudopotential)>,
    options: ScfOptions,
}
impl Model {
    fn refs(&self) -> Vec<SpeciesPotential<'_>> {
        self.species
            .iter()
            .map(|(label, p)| SpeciesPotential {
                label,
                potential: p,
            })
            .collect()
    }
    fn problem(&self, options: ScfOptions) -> ScfProblem {
        ScfProblem::new(&self.cell, &self.atoms, &self.refs(), options).unwrap()
    }
}

/// Skewed two-species synthetic fixture (coupled s,s,p projectors on Xx, local Yy).
fn synthetic_skew() -> Model {
    Model {
        cell: Cell::from_user_rows_bohr([[5., 0., 0.], [0.5, 4., 0.], [0.2, 0.4, 6.]]).unwrap(),
        atoms: vec![
            Atom::new("Xx", [0.13, 0.21, 0.32]),
            Atom::new("Yy", [0.6, 0.7, 0.8]),
        ],
        species: vec![
            ("Xx".into(), common::potential(2.0, 1.0, true, 401)),
            ("Yy".into(), common::potential(2.0, 1.3, false, 401)),
        ],
        options: ScfOptions {
            cutoff_hartree: 4.0,
            bands: Some(4),
            energy_tolerance: 1e-10,
            ..Default::default()
        },
    }
}

/// Elongated (slab-like) one-species synthetic fixture.
fn synthetic_slab() -> Model {
    Model {
        cell: Cell::from_user_rows_bohr([[4., 0., 0.], [0., 4., 0.], [0., 0., 10.]]).unwrap(),
        atoms: vec![
            Atom::new("Xx", [0.1, 0.2, 0.0]),
            Atom::new("Yy", [0.6, 0.7, 0.2]),
        ],
        species: vec![
            ("Xx".into(), common::potential(2.0, 1.0, true, 401)),
            ("Yy".into(), common::potential(2.0, 1.3, false, 401)),
        ],
        options: ScfOptions {
            cutoff_hartree: 4.0,
            bands: Some(4),
            energy_tolerance: 1e-10,
            ..Default::default()
        },
    }
}

#[test]
fn orbital_prolongation_is_an_exact_injection() {
    let cell = synthetic_skew().cell;
    let coarse = PlaneWaveBasis::new(&cell, Hartree(2.0)).unwrap();
    let fine = PlaneWaveBasis::new(&cell, Hartree(4.0)).unwrap();
    let n = coarse.vectors().len();
    let mut c = ComplexMatrix::zeros(n, 3).unwrap();
    for j in 0..3 {
        for (i, z) in c.column_mut(j).iter_mut().enumerate() {
            *z = Complex64::new(((i * 7 + j * 3) as f64).sin(), ((i + 5 * j) as f64).cos());
        }
    }
    let p = cutoff_ladder::prolong_orbitals(&coarse, &fine, &c).unwrap();
    assert_eq!(p.rows(), fine.vectors().len());
    // Every coarse coefficient sits at its Miller index; every new G is zero;
    // all inner products are unchanged bit for bit (same terms, zeros added).
    let mut copied = 0;
    for (i, g) in fine.vectors().iter().enumerate() {
        match coarse.vectors().iter().position(|h| h.miller == g.miller) {
            Some(k) => {
                copied += 1;
                for j in 0..3 {
                    assert_eq!(p.column(j)[i], c.column(j)[k]);
                }
            }
            None => (0..3).for_each(|j| assert_eq!(p.column(j)[i], Complex64::default())),
        }
    }
    assert_eq!(copied, n);
    // Bases are not nested the other way round.
    assert!(cutoff_ladder::prolong_orbitals(&fine, &coarse, &p).is_err());
}

#[test]
fn density_prolongation_is_exact_fourier_interpolation() {
    let cell = synthetic_skew().cell;
    let wave = PlaneWaveBasis::new(&cell, Hartree(2.0)).unwrap();
    let density = wave.density_basis(&cell).unwrap();
    let coarse = FftGridSpec::new(&density, &density, None).unwrap();
    let fine_wave = PlaneWaveBasis::new(&cell, Hartree(4.0)).unwrap();
    let fine_density = fine_wave.density_basis(&cell).unwrap();
    let fine = FftGridSpec::new(&fine_density, &fine_density, None).unwrap();
    // A real band-limited field on the coarse density sphere, evaluated directly
    // at the sample points of both grids.
    let terms: Vec<_> = density
        .vectors()
        .iter()
        .enumerate()
        .map(|(i, g)| {
            (
                g.miller,
                1.0 / (1.0 + g.norm_squared) * (1.0 + 0.3 * (i as f64).sin()),
            )
        })
        .collect();
    let field = |dims: [usize; 3]| {
        Array3::from_shape_fn(dims, |(a, b, c)| {
            let r = [
                a as f64 / dims[0] as f64,
                b as f64 / dims[1] as f64,
                c as f64 / dims[2] as f64,
            ];
            // Symmetrised so that the field is real: sum over +G and -G.
            terms
                .iter()
                .map(|(m, w)| {
                    let phase = std::f64::consts::TAU
                        * (m[0] as f64 * r[0] + m[1] as f64 * r[1] + m[2] as f64 * r[2]);
                    w * phase.cos()
                })
                .sum::<f64>()
                + 3.0
        })
    };
    let on_coarse = field(coarse.dims());
    let expected = field(fine.dims());
    let prolonged = cutoff_ladder::prolong_density(&on_coarse, fine.dims()).unwrap();
    let scale = expected.iter().fold(0f64, |m, v| m.max(v.abs()));
    let error = prolonged
        .iter()
        .zip(&expected)
        .fold(0f64, |m, (a, b)| m.max((a - b).abs()));
    assert!(error < 1e-12 * scale, "interpolation error {error:e}");
    let mean = |a: &Array3<f64>| a.sum() / a.len() as f64;
    assert!((mean(&prolonged) - mean(&on_coarse)).abs() < 1e-13 * scale);
}

fn check_ladder(model: &Model) {
    let plain = model.problem(model.options.clone()).run_outcome().unwrap();
    assert!(plain.converged);
    let ladder = CutoffLadderOptions::default();
    let out = cutoff_ladder::run_ladder(
        model.problem(model.options.clone()),
        &model.refs(),
        &ladder,
        true,
    )
    .unwrap();
    let r = &out.report;
    assert!(r.used, "fallback: {:?}", r.fallback_reason);
    assert!(out.outcome.converged);
    assert!(
        out.outcome.result.certified
            || out
                .outcome
                .result
                .history
                .last()
                .unwrap()
                .fixed_point_certificate
                .is_some()
    );
    assert!(r.coarse_iterations > 0 && r.coarse_plane_waves < r.fine_plane_waves);
    assert!(r.coarse_grid.iter().product::<usize>() < r.fine_grid.iter().product::<usize>());
    // The full-cutoff result: same grid, same model, energy within the strict class.
    let a = &plain.result;
    let b = &out.outcome.result;
    assert_eq!(a.density.dim(), b.density.dim());
    let de = (a.total_hartree() - b.total_hartree()).abs();
    assert!(de < 1e-10, "energy difference {de:e} Ha");
    let occupied = a.occupations.iter().filter(|f| **f > 0.).count();
    for (x, y) in a.eigenvalues[..occupied]
        .iter()
        .zip(&b.eigenvalues[..occupied])
    {
        assert!((x - y).abs() < 1e-7, "eigenvalue {x} vs {y}");
    }
    println!(
        "plain {} SCF, {} H; ladder coarse {} SCF ({} H, {:.0} full-grid eq.), fine {} SCF, {} H; dE {de:e} Ha",
        a.history.len(),
        a.history
            .iter()
            .map(|h| h.operator_applications)
            .sum::<usize>(),
        r.coarse_iterations,
        r.coarse_h_vectors,
        r.coarse_h_equivalents(),
        b.history.len(),
        b.history
            .iter()
            .map(|h| h.operator_applications)
            .sum::<usize>(),
    );
}

#[test]
fn ladder_reaches_the_full_cutoff_fixed_point_skewed() {
    check_ladder(&synthetic_skew());
}

#[test]
fn ladder_reaches_the_full_cutoff_fixed_point_slab() {
    check_ladder(&synthetic_slab());
}

#[test]
fn ineligible_problem_falls_back_bit_identically() {
    let model = synthetic_skew();
    let options = ScfOptions {
        gamma_real: GammaReal::Off,
        ..model.options.clone()
    };
    let plain = model.problem(options.clone()).run_outcome().unwrap();
    let out = cutoff_ladder::run_ladder(
        model.problem(options),
        &model.refs(),
        &CutoffLadderOptions::default(),
        true,
    )
    .unwrap();
    assert!(!out.report.used);
    assert_eq!(
        out.report.fallback_reason.as_deref(),
        Some("gamma_real off")
    );
    assert_eq!(
        plain.result.total_hartree().to_bits(),
        out.outcome.result.total_hartree().to_bits()
    );
    assert_eq!(plain.result.density, out.outcome.result.density);
    assert!(
        CutoffLadderOptions {
            fraction: 1.0,
            ..Default::default()
        }
        .validate()
        .is_err()
    );
}

#[test]
fn ladder_result_is_independent_of_thread_count() {
    let model = synthetic_slab();
    let run = |threads| {
        cutoff_ladder::run_ladder(
            model.problem(ScfOptions {
                threads,
                ..model.options.clone()
            }),
            &model.refs(),
            &CutoffLadderOptions::default(),
            true,
        )
        .unwrap()
    };
    let one = run(1);
    let four = run(4);
    assert!(one.report.used && four.report.used);
    assert_eq!(one.report.coarse_iterations, four.report.coarse_iterations);
    assert_eq!(
        one.outcome.result.total_hartree().to_bits(),
        four.outcome.result.total_hartree().to_bits()
    );
    assert_eq!(one.outcome.result.density, four.outcome.result.density);
}
