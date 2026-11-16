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
            ..ScfOptions::strict_reference()
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
            ..ScfOptions::strict_reference()
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

// ---------------------------------------------------------------------------
// S1 (standard-profile ladder) and S2 (secant seeds across the hand-off).

fn standard(options: &ScfOptions) -> ScfOptions {
    ScfOptions {
        standard: Some(dft_engine::StandardProfile::default()),
        ..options.clone()
    }
}
fn strict(options: &ScfOptions) -> ScfOptions {
    options.clone()
}

#[test]
fn standard_ladder_is_opt_in_and_stays_in_the_standard_class() {
    let model = synthetic_skew();
    let strict = model.problem(model.options.clone()).run_outcome().unwrap();
    assert!(strict.converged);
    let reference = strict.result.total_hartree();
    // Off (the adopted default): fallback with the old reason.
    let off = cutoff_ladder::run_ladder(
        model.problem(standard(&model.options)),
        &model.refs(),
        &CutoffLadderOptions::baseline(),
        true,
    )
    .unwrap();
    assert!(!off.report.used);
    assert_eq!(
        off.report.fallback_reason.as_deref(),
        Some("standard profile")
    );
    for seed in [
        cutoff_ladder::StandardSeed::Coarse,
        cutoff_ladder::StandardSeed::Cap,
    ] {
        let on = cutoff_ladder::run_ladder(
            model.problem(standard(&model.options)),
            &model.refs(),
            &CutoffLadderOptions {
                standard: true,
                standard_seed: seed,
                ..CutoffLadderOptions::baseline()
            },
            true,
        )
        .unwrap();
        let r = &on.report;
        assert!(r.used && r.standard, "fallback: {:?}", r.fallback_reason);
        assert_eq!(
            r.standard_seed_ha.is_some(),
            seed == cutoff_ladder::StandardSeed::Coarse
        );
        assert!(on.outcome.converged);
        // Standard stays uncertified; its delivered class is 5e-9 Ha.
        assert!(!on.outcome.result.certified);
        let de = (on.outcome.result.total_hartree() - reference).abs();
        println!(
            "standard ladder {seed:?}: seed {:?} Ha, fine {} SCF, |E - E_strict| {de:e} Ha",
            r.standard_seed_ha,
            on.outcome.result.history.len()
        );
        assert!(de < 5e-9, "standard ladder error {de:e} Ha");
    }
}

#[test]
fn secant_seeds_keep_the_strict_fixed_point() {
    for model in [synthetic_skew(), synthetic_slab()] {
        let plain = model.problem(model.options.clone()).run_outcome().unwrap();
        for h in [2, 6] {
            let out = cutoff_ladder::run_ladder(
                model.problem(model.options.clone()),
                &model.refs(),
                &CutoffLadderOptions {
                    secant_pairs: h,
                    ..CutoffLadderOptions::baseline()
                },
                true,
            )
            .unwrap();
            let r = &out.report;
            assert!(r.used, "fallback: {:?}", r.fallback_reason);
            assert_eq!(r.secant_requested, h);
            assert!(
                r.secant_seeded > 0 && r.secant_seeded <= h,
                "seeded {} ({:?})",
                r.secant_seeded,
                r.secant_fallback
            );
            assert!(out.outcome.converged);
            let de = (plain.result.total_hartree() - out.outcome.result.total_hartree()).abs();
            println!(
                "secant h={h}: seeded {}, fine {} SCF, dE {de:e} Ha",
                r.secant_seeded,
                out.outcome.result.history.len()
            );
            assert!(de < 1e-10, "energy difference {de:e} Ha");
        }
    }
}

#[test]
fn standard_and_secant_ladder_is_independent_of_thread_count() {
    let model = synthetic_slab();
    let ladder = CutoffLadderOptions {
        standard: true,
        secant_pairs: 4,
        ..CutoffLadderOptions::baseline()
    };
    let profiles: [fn(&ScfOptions) -> ScfOptions; 2] = [standard, strict];
    for profile in profiles {
        let run = |threads| {
            cutoff_ladder::run_ladder(
                model.problem(ScfOptions {
                    threads,
                    ..profile(&model.options)
                }),
                &model.refs(),
                &ladder,
                true,
            )
            .unwrap()
        };
        let one = run(1);
        let four = run(4);
        assert!(one.report.used && one.report.secant_seeded > 0);
        assert_eq!(one.report.secant_seeded, four.report.secant_seeded);
        assert_eq!(
            one.outcome.result.total_hartree().to_bits(),
            four.outcome.result.total_hartree().to_bits()
        );
        assert_eq!(one.outcome.result.density, four.outcome.result.density);
    }
}

// ---------------------------------------------------------------------------
// D1 hand-off knobs (coarse stop level and rule, first fine tolerance).

use cutoff_ladder::{CoarseStopRule, FirstTolerance};

/// Cubic synthetic fixture (simple cubic cell, one nonlocal atom; the Gamma
/// fixture of the stepwise-equivalence matrix).
fn synthetic_cubic() -> Model {
    Model {
        cell: Cell::from_user_rows_bohr([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]).unwrap(),
        atoms: vec![Atom::new("Xx", [0.13, 0.21, 0.32])],
        species: vec![("Xx".into(), common::potential(2.0, 1.0, true, 401))],
        options: ScfOptions {
            cutoff_hartree: 4.0,
            bands: Some(3),
            energy_tolerance: 1e-10,
            ..ScfOptions::strict_reference()
        },
    }
}

fn run_with(
    model: &Model,
    options: ScfOptions,
    ladder: &CutoffLadderOptions,
) -> cutoff_ladder::LadderOutcome {
    cutoff_ladder::run_ladder(model.problem(options), &model.refs(), ladder, true).unwrap()
}

fn same_bits(a: &cutoff_ladder::LadderOutcome, b: &cutoff_ladder::LadderOutcome) {
    assert_eq!(
        a.outcome.result.total_hartree().to_bits(),
        b.outcome.result.total_hartree().to_bits()
    );
    assert_eq!(a.outcome.result.density, b.outcome.result.density);
    assert_eq!(
        a.outcome.result.history.len(),
        b.outcome.result.history.len()
    );
    assert_eq!(a.report.coarse_iterations, b.report.coarse_iterations);
    assert_eq!(a.report.coarse_h_vectors, b.report.coarse_h_vectors);
}

#[test]
fn d1_defaults_are_the_adopted_ladder() {
    // The test process leaves the D1 variables unset: Default keeps the
    // pre-registered stop level, the fixed rule and the adopted first step.
    let d = CutoffLadderOptions::default();
    let b = CutoffLadderOptions::baseline();
    assert_eq!(d.fraction, 0.5);
    assert_eq!(d.coarse_density_residual.to_bits(), 1e-4f64.to_bits());
    assert_eq!(b.coarse_density_residual.to_bits(), 1e-4f64.to_bits());
    assert_eq!(d.coarse_rule, CoarseStopRule::Fixed);
    assert_eq!(d.first_tolerance, FirstTolerance::Auto);
    assert_eq!(b.first_tolerance, FirstTolerance::Auto);
    // Auto equals the explicit equivalents bit for bit: strict Auto = Cap
    // (cold schedule), standard Auto with the coarse seed = Next.
    let model = synthetic_skew();
    let auto = run_with(
        &model,
        model.options.clone(),
        &CutoffLadderOptions::baseline(),
    );
    let cap = run_with(
        &model,
        model.options.clone(),
        &CutoffLadderOptions {
            first_tolerance: FirstTolerance::Cap,
            ..CutoffLadderOptions::baseline()
        },
    );
    assert!(auto.report.used && auto.report.strict_seed_ha.is_none());
    same_bits(&auto, &cap);
    let std_ladder = |first_tolerance| CutoffLadderOptions {
        standard: true,
        first_tolerance,
        ..CutoffLadderOptions::baseline()
    };
    let auto = run_with(
        &model,
        standard(&model.options),
        &std_ladder(FirstTolerance::Auto),
    );
    let next = run_with(
        &model,
        standard(&model.options),
        &std_ladder(FirstTolerance::Next),
    );
    assert!(auto.report.standard_seed_ha.is_some());
    assert_eq!(auto.report.standard_seed_ha, next.report.standard_seed_ha);
    same_bits(&auto, &next);
    // The reported fine work matches the fine history.
    let r = &auto.report;
    assert_eq!(r.fine_iterations, auto.outcome.result.history.len());
    assert!(r.fine_h_vectors > 0 && r.coarse_fine_step_equivalents().unwrap() > 0.0);
    assert!(r.handoff_to_final_density_distance.unwrap() > 0.0);
}

/// Each D1 knob changes only its documented quantity, and the full-cutoff
/// session still converges to the strict fixed point (1e-10 Ha) on the cubic,
/// skewed and slab fixtures.
#[test]
fn d1_knobs_change_only_their_quantity_and_the_fine_stage_converges() {
    for model in [synthetic_cubic(), synthetic_skew(), synthetic_slab()] {
        let plain = model.problem(model.options.clone()).run_outcome().unwrap();
        assert!(plain.converged);
        let reference = plain.result.total_hartree();
        let base = run_with(
            &model,
            model.options.clone(),
            &CutoffLadderOptions::baseline(),
        );
        assert!(
            base.report.used,
            "fallback: {:?}",
            base.report.fallback_reason
        );
        // Coarse stop level: only the coarse stage changes.
        for level in [1e-2, 1e-3, 3e-3] {
            let out = run_with(
                &model,
                model.options.clone(),
                &CutoffLadderOptions {
                    coarse_density_residual: level,
                    ..CutoffLadderOptions::baseline()
                },
            );
            let r = &out.report;
            assert!(r.used && out.outcome.converged);
            assert_eq!(r.coarse_stop_level, level);
            assert!(r.coarse_iterations <= base.report.coarse_iterations);
            assert!(
                r.coarse_final_density_residual <= level
                    || r.coarse_iterations == CutoffLadderOptions::baseline().coarse_max_iterations
            );
            assert!(r.strict_seed_ha.is_none());
            let de = (out.outcome.result.total_hartree() - reference).abs();
            println!(
                "stop {level:e}: coarse {} SCF ({:.0} H eq, {:.2} fine-step eq), fine {} SCF, dE {de:e} Ha",
                r.coarse_iterations,
                r.coarse_h_equivalents(),
                r.coarse_fine_step_equivalents().unwrap_or(f64::NAN),
                r.fine_iterations
            );
            assert!(de < 1e-10, "stop {level:e}: energy difference {de:e} Ha");
        }
        // First fine tolerance: the coarse stage is untouched, the first fine
        // step requests the seed (inside the schedule's bounds).
        for first_tolerance in [FirstTolerance::Coarse, FirstTolerance::Next] {
            let out = run_with(
                &model,
                model.options.clone(),
                &CutoffLadderOptions {
                    first_tolerance,
                    ..CutoffLadderOptions::baseline()
                },
            );
            let r = &out.report;
            assert!(r.used && out.outcome.converged);
            assert_eq!(r.coarse_iterations, base.report.coarse_iterations);
            assert_eq!(r.coarse_h_vectors, base.report.coarse_h_vectors);
            let seed = r.strict_seed_ha.expect("strict seed");
            assert!(seed < 1e-2);
            let first = out.outcome.result.history[0].requested_eigen_tolerance;
            assert_eq!(
                first.to_bits(),
                seed.min(1e-2)
                    .max(0.1 * model.options.eigen_tolerance)
                    .to_bits()
            );
            let de = (out.outcome.result.total_hartree() - reference).abs();
            println!(
                "first {first_tolerance:?}: seed {seed:e} Ha, fine {} SCF ({} H) vs {} SCF ({} H), dE {de:e} Ha",
                r.fine_iterations,
                r.fine_h_vectors,
                base.report.fine_iterations,
                base.report.fine_h_vectors
            );
            assert!(
                de < 1e-10,
                "{first_tolerance:?}: energy difference {de:e} Ha"
            );
        }
    }
}

#[test]
fn atomic_coarse_rule_only_raises_the_stop_level() {
    let mut model = synthetic_skew();
    model.species = vec![
        (
            "Xx".into(),
            common::potential_with_atomic_charge(2.0, 1.0, true, 401),
        ),
        (
            "Yy".into(),
            common::potential_with_atomic_charge(2.0, 1.3, false, 401),
        ),
    ];
    // Uniform start: no atomic distance, the fixed level applies.
    let uniform = run_with(
        &model,
        model.options.clone(),
        &CutoffLadderOptions {
            coarse_rule: CoarseStopRule::Atomic,
            ..CutoffLadderOptions::baseline()
        },
    );
    assert!(uniform.report.atomic_truncation_distance.is_none());
    same_bits(
        &uniform,
        &run_with(
            &model,
            model.options.clone(),
            &CutoffLadderOptions::baseline(),
        ),
    );
    let atomic = ScfOptions {
        start_density: dft_engine::StartDensity::Atomic,
        ..model.options.clone()
    };
    let fixed = run_with(&model, atomic.clone(), &CutoffLadderOptions::baseline());
    let rule = run_with(
        &model,
        atomic.clone(),
        &CutoffLadderOptions {
            coarse_rule: CoarseStopRule::Atomic,
            ..CutoffLadderOptions::baseline()
        },
    );
    let d = rule
        .report
        .atomic_truncation_distance
        .expect("atomic distance");
    assert!(d > 0.0 && d.is_finite());
    assert_eq!(fixed.report.atomic_truncation_distance, Some(d));
    assert_eq!(fixed.report.coarse_stop_level, 1e-4);
    assert_eq!(rule.report.coarse_stop_level, d.max(1e-4));
    assert!(rule.report.coarse_iterations <= fixed.report.coarse_iterations);
    assert!(rule.outcome.converged);
    let plain = model.problem(atomic).run_outcome().unwrap();
    let de = (plain.result.total_hartree() - rule.outcome.result.total_hartree()).abs();
    println!(
        "atomic rule: d {d:e}, coarse {} vs {} SCF, measured coarse-to-fine {:e}, dE {de:e} Ha",
        rule.report.coarse_iterations,
        fixed.report.coarse_iterations,
        rule.report
            .handoff_to_final_density_distance
            .unwrap_or(f64::NAN)
    );
    assert!(de < 1e-10, "energy difference {de:e} Ha");
}

#[test]
fn d1_knobs_on_a_k_point_problem_fall_back_bit_identically() {
    // Four k points in a cubic cell (the stepwise-equivalence k fixture).
    let model = Model {
        cell: Cell::from_user_rows_bohr([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]).unwrap(),
        atoms: vec![Atom::new("Xx", [0.13, 0.21, 0.32])],
        species: vec![("Xx".into(), common::potential(2.0, 1.0, true, 401))],
        options: ScfOptions {
            cutoff_hartree: 2.5,
            bands: Some(3),
            ..ScfOptions::strict_reference()
        },
    };
    let points: Vec<_> = (0..4)
        .map(|i| dft_geometry::kpoints::KPoint {
            fractional: [0.025 * i as f64, 0.015 * i as f64, 0.],
            weight: 0.25,
        })
        .collect();
    let problem = || {
        ScfProblem::with_kpoints(
            &model.cell,
            &model.atoms,
            &model.refs(),
            model.options.clone(),
            &points,
            dft_engine::OccupationPolicy::Insulator {
                minimum_gap_hartree: 1e-7,
            },
        )
        .unwrap()
    };
    let plain = problem().run_outcome().unwrap();
    let out = cutoff_ladder::run_ladder(
        problem(),
        &model.refs(),
        &CutoffLadderOptions {
            coarse_density_residual: 1e-2,
            coarse_rule: CoarseStopRule::Atomic,
            first_tolerance: FirstTolerance::Coarse,
            ..CutoffLadderOptions::baseline()
        },
        true,
    )
    .unwrap();
    assert!(!out.report.used);
    assert_eq!(
        out.report.fallback_reason.as_deref(),
        Some("multiple k points")
    );
    assert!(out.report.strict_seed_ha.is_none() && out.report.fine_iterations == 0);
    assert_eq!(
        plain.result.total_hartree().to_bits(),
        out.outcome.result.total_hartree().to_bits()
    );
    assert_eq!(plain.result.density, out.outcome.result.density);
}

#[test]
fn d1_knobs_are_independent_of_thread_count() {
    let model = synthetic_slab();
    let ladder = CutoffLadderOptions {
        coarse_density_residual: 3e-3,
        first_tolerance: FirstTolerance::Coarse,
        ..CutoffLadderOptions::baseline()
    };
    let run = |threads| {
        run_with(
            &model,
            ScfOptions {
                threads,
                ..model.options.clone()
            },
            &ladder,
        )
    };
    let one = run(1);
    let four = run(4);
    assert!(one.report.used && one.report.strict_seed_ha.is_some());
    assert_eq!(
        one.report.strict_seed_ha.map(f64::to_bits),
        four.report.strict_seed_ha.map(f64::to_bits)
    );
    same_bits(&one, &four);
}

// ---------------------------------------------------------------------------
// D1 Mac count scan (release only):
// `cargo test --release -p dft-engine --test cutoff_ladder ladder_d1_scan -- --ignored --nocapture --test-threads=1`
// Env: D1_CASES (default si16,si32), D1_AMP (default f1e-3), D1_SEED (default 601),
// D1_PROFILES (default standard), D1_STOPS (default 1e-4,3e-3,1e-3,1e-2),
// D1_FRACTIONS (default 0.5,0.7), D1_FIRST (default auto,coarse), D1_ATOMIC=on adds the
// atomic rule arm (needs D1_START=atomic; default uniform start density as the CLI),
// LADDER_THREADS (default 4).
#[test]
#[ignore = "release-only D1 scan (minutes on a Mac)"]
fn ladder_d1_scan() {
    let cases = env_list("D1_CASES", "si16,si32");
    let amp = std::env::var("D1_AMP").unwrap_or_else(|_| "f1e-3".into());
    let seed: u64 = std::env::var("D1_SEED")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(601);
    let profiles = env_list("D1_PROFILES", "standard");
    let stops: Vec<f64> = env_list("D1_STOPS", "1e-4,3e-3,1e-3,1e-2")
        .iter()
        .map(|s| s.parse().unwrap())
        .collect();
    let fractions: Vec<f64> = env_list("D1_FRACTIONS", "0.5,0.7")
        .iter()
        .map(|s| s.parse().unwrap())
        .collect();
    let firsts = env_list("D1_FIRST", "auto,coarse");
    let atomic = std::env::var("D1_ATOMIC").as_deref() == Ok("on");
    let threads: usize = std::env::var("LADDER_THREADS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(4);
    let (fractional, amplitude) = match amp.split_at(1) {
        ("b", v) => (false, v.parse::<f64>().unwrap()),
        ("f", v) => (true, v.parse::<f64>().unwrap()),
        _ => panic!("D1_AMP starts with b or f"),
    };
    println!(
        "case | profile | arm | conv | fine_scf | fine_H | fine_cert_H | coarse_scf | coarse_H_eq | coarse_step_eq | total_step_eq | total_H_eq | stop_level | atomic_d | handoff_d | seed_tol_Ha | dE_vs_strict_plain_Ha | wall_s"
    );
    for case in &cases {
        let Some(mut model) = bench_model(case) else {
            println!("{case} | skipped (fixture missing)");
            continue;
        };
        if seed != 0 {
            model.perturb(seed, amplitude, fractional);
        }
        let base = ScfOptions {
            threads,
            start_density: if std::env::var("D1_START").as_deref() == Ok("atomic") {
                dft_engine::StartDensity::Atomic
            } else {
                dft_engine::StartDensity::Uniform
            },
            ..model.options.clone()
        };
        let strict_plain = model.problem(base.clone()).run_outcome().unwrap();
        let e_ref = strict_plain.result.total_hartree();
        for profile in &profiles {
            let options = ScfOptions {
                standard: (profile == "standard").then(dft_engine::StandardProfile::default),
                ..base.clone()
            };
            let started = std::time::Instant::now();
            let plain = model.problem(options.clone()).run_outcome().unwrap();
            let wall = started.elapsed().as_secs_f64();
            let (scf, h) = fine_work(&plain);
            println!(
                "{case} | {profile} | off | {} | {scf} | {h} | - | 0 | 0 | 0 | {scf} | {h} | - | - | - | - | {:.2e} | {wall:.2}",
                plain.converged,
                (plain.result.total_hartree() - e_ref).abs()
            );
            let mut arms = Vec::new();
            for &fraction in &fractions {
                for &stop in &stops {
                    for first in &firsts {
                        arms.push((fraction, stop, first.clone(), CoarseStopRule::Fixed));
                    }
                }
                if atomic {
                    arms.push((fraction, 1e-4, "auto".into(), CoarseStopRule::Atomic));
                }
            }
            for (fraction, stop, first, rule) in arms {
                let ladder = CutoffLadderOptions {
                    fraction,
                    coarse_density_residual: stop,
                    coarse_rule: rule,
                    standard: true,
                    first_tolerance: match first.as_str() {
                        "cap" => FirstTolerance::Cap,
                        "coarse" => FirstTolerance::Coarse,
                        "next" => FirstTolerance::Next,
                        _ => FirstTolerance::Auto,
                    },
                    ..CutoffLadderOptions::baseline()
                };
                let started = std::time::Instant::now();
                let out = cutoff_ladder::run_ladder(
                    model.problem(options.clone()),
                    &model.refs(),
                    &ladder,
                    true,
                )
                .unwrap();
                let wall = started.elapsed().as_secs_f64();
                let r = &out.report;
                let step_eq = r.coarse_fine_step_equivalents().unwrap_or(f64::NAN);
                println!(
                    "{case} | {profile} | E1={fraction} stop={stop:e} first={first} rule={rule:?} | {} | {} | {} | {} | {} | {:.0} | {step_eq:.2} | {:.2} | {:.0} | {:.1e} | {} | {} | {} | {:.2e} | {wall:.2}",
                    out.outcome.converged,
                    r.fine_iterations,
                    r.fine_h_vectors,
                    r.fine_certificate_h_vectors,
                    r.coarse_iterations,
                    r.coarse_h_equivalents(),
                    r.fine_iterations as f64 + step_eq,
                    r.coarse_h_equivalents()
                        + (r.fine_h_vectors + r.fine_certificate_h_vectors) as f64,
                    r.coarse_stop_level,
                    r.atomic_truncation_distance
                        .map_or("-".into(), |d| format!("{d:.1e}")),
                    r.handoff_to_final_density_distance
                        .map_or("-".into(), |d| format!("{d:.1e}")),
                    r.standard_seed_ha
                        .or(r.strict_seed_ha)
                        .map_or("-".into(), |t| format!("{t:.1e}")),
                    (out.outcome.result.total_hartree() - e_ref).abs(),
                );
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Mac count bench for S1/S2 (release only, minutes):
// `cargo test --release -p dft-engine --test cutoff_ladder -- --ignored --nocapture --test-threads=1`
// Env: LADDER_CASES (comma list), LADDER_SEEDS (default 601,602,603; 0 = unperturbed),
// LADDER_AMPS (b1e-7 = +-1e-7 Bohr Cartesian, f1e-3 = +-0.001 fractional; default both),
// LADDER_PROFILES (strict,standard), LADDER_SECANT (h list, default 2,4,6),
// LADDER_THREADS (default 4; results are thread-count independent),
// LADDER_STD_SEED=cap (standard first-step tolerance rule; default coarse),
// DFT_SI_UPF, DFT_PSEUDODOJO_DIR.

struct BenchModel {
    cell: Cell,
    atoms: Vec<Atom>,
    species: Vec<(String, dft_upf::Pseudopotential)>,
    options: ScfOptions,
}
impl BenchModel {
    /// Deterministic xorshift displacement, uniform in [-a, a] per component:
    /// Cartesian Bohr (`fractional = false`) or fractional coordinates.
    fn perturb(&mut self, seed: u64, amplitude: f64, fractional: bool) {
        let mut x = seed ^ 0x2545_f491_4f6c_dd1d;
        let mut uniform = || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            2. * ((x >> 11) as f64 / (1u64 << 53) as f64) - 1.
        };
        for atom in &mut self.atoms {
            let shift: [f64; 3] = std::array::from_fn(|_| amplitude * uniform());
            let df = if fractional {
                shift
            } else {
                self.cell.cartesian_bohr_to_fractional(shift)
            };
            for (f, d) in atom.fractional_position.iter_mut().zip(df) {
                *f += d;
            }
        }
    }
    fn refs(&self) -> Vec<SpeciesPotential<'_>> {
        self.species
            .iter()
            .map(|(label, p)| SpeciesPotential {
                label,
                potential: p,
            })
            .collect()
    }
    /// As `dft-cli run`: an NLCC model gets the XC density floor 1e-10.
    fn problem(&self, options: ScfOptions) -> ScfProblem {
        let mut p = ScfProblem::new(&self.cell, &self.atoms, &self.refs(), options).unwrap();
        if p.core_density().is_some() {
            p.set_density_floor(1e-10).unwrap();
        }
        p
    }
}

fn si_upf() -> Option<dft_upf::Pseudopotential> {
    let path = std::env::var("DFT_SI_UPF")
        .unwrap_or_else(|_| "/Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF".into());
    let text = std::fs::read_to_string(path).ok()?;
    Some(dft_upf::parse_upf(&text).unwrap())
}
fn dojo(element: &str) -> Option<dft_upf::Pseudopotential> {
    let dir = std::env::var("DFT_PSEUDODOJO_DIR")
        .unwrap_or_else(|_| "/Users/mino/Develop/rust-dft-wt/pseudo/pseudodojo/ex".into());
    let text = std::fs::read_to_string(format!("{dir}/{element}.upf")).ok()?;
    Some(dft_upf::parse_upf(&text).unwrap())
}

const A_SI: f64 = 10.261_212_856_717_933;
const DIAMOND: [[f64; 3]; 8] = [
    [0., 0., 0.],
    [0., 0.5, 0.5],
    [0.5, 0., 0.5],
    [0.5, 0.5, 0.],
    [0.25, 0.25, 0.25],
    [0.25, 0.75, 0.75],
    [0.75, 0.25, 0.75],
    [0.75, 0.75, 0.25],
];

fn si_bench(rows: [[f64; 3]; 3], atoms: Vec<[f64; 3]>, ry: f64) -> Option<BenchModel> {
    Some(BenchModel {
        cell: Cell::from_user_rows_bohr(rows).unwrap(),
        atoms: atoms.into_iter().map(|p| Atom::new("Si", p)).collect(),
        species: vec![("Si".into(), si_upf()?)],
        options: ScfOptions {
            cutoff_hartree: ry / 2.,
            max_iterations: 100,
            ..ScfOptions::strict_reference()
        },
    })
}

fn bench_model(id: &str) -> Option<BenchModel> {
    let fcc = |shear: [f64; 2]| {
        let h = A_SI / 2.;
        [[shear[0], h, h], [h, 0., h], [h, h, shear[1]]]
    };
    let synthetic = |m: Model| BenchModel {
        cell: m.cell,
        atoms: m.atoms,
        species: m.species,
        options: ScfOptions {
            max_iterations: 100,
            ..m.options
        },
    };
    match id {
        "si8" => si_bench(
            [[A_SI, 0., 0.], [0., A_SI, 0.], [0., 0., A_SI]],
            DIAMOND.to_vec(),
            30.,
        ),
        "si16" => si_bench(
            [[A_SI, 0., 0.], [0., A_SI, 0.], [0., 0., 2. * A_SI]],
            DIAMOND
                .iter()
                .flat_map(|p| [[p[0], p[1], p[2] / 2.], [p[0], p[1], p[2] / 2. + 0.5]])
                .collect(),
            30.,
        ),
        "si32" => si_bench(
            [[2. * A_SI, 0., 0.], [0., 2. * A_SI, 0.], [0., 0., A_SI]],
            DIAMOND
                .iter()
                .flat_map(|p| {
                    [[0., 0.], [0., 0.5], [0.5, 0.], [0.5, 0.5]]
                        .map(|s| [p[0] / 2. + s[0], p[1] / 2. + s[1], p[2]])
                })
                .collect(),
            30.,
        ),
        "si-skew" => si_bench(fcc([0.7, 0.3]), vec![[0.; 3], [0.25; 3]], 30.),
        "si-low" => si_bench(fcc([0., 0.]), vec![[0.; 3], [0.25; 3]], 8.),
        "mgo8" => {
            // Rocksalt MgO conventional cell, PseudoDojo PBE (O has NLCC); 35 Ha
            // is an algorithmic test cutoff, not a converged value.
            let (mg, o) = (dojo("Mg")?, dojo("O")?);
            let a = 7.956;
            let mut atoms = Vec::new();
            for s in [[0., 0., 0.], [0., 0.5, 0.5], [0.5, 0., 0.5], [0.5, 0.5, 0.]] {
                atoms.push(Atom::new("Mg", s));
                atoms.push(Atom::new("O", [s[0] + 0.5, s[1], s[2]].map(|x| x % 1.)));
            }
            Some(BenchModel {
                cell: Cell::from_user_rows_bohr([[a, 0., 0.], [0., a, 0.], [0., 0., a]]).unwrap(),
                atoms,
                species: vec![("Mg".into(), mg), ("O".into(), o)],
                options: ScfOptions {
                    cutoff_hartree: 35.,
                    functional: dft_engine::Functional::Pbe,
                    max_iterations: 100,
                    ..ScfOptions::strict_reference()
                },
            })
        }
        "syn-skew2" => Some(synthetic(synthetic_skew())),
        "syn-slab" => Some(synthetic(synthetic_slab())),
        "syn-fcc" => Some(BenchModel {
            cell: common::cell(),
            atoms: common::atoms(),
            species: vec![("Xx".into(), common::potential(4.0, 1.0, true, 401))],
            options: ScfOptions {
                cutoff_hartree: 6.,
                max_iterations: 100,
                ..ScfOptions::strict_reference()
            },
        }),
        "syn-lowsym6" => {
            let p = common::potential(2.0, 1.0, true, 401);
            let q = common::potential(2.0, 1.3, false, 401);
            let frac = [
                [0.03, 0.07, 0.11],
                [0.52, 0.18, 0.09],
                [0.21, 0.61, 0.27],
                [0.73, 0.55, 0.41],
                [0.37, 0.29, 0.68],
                [0.81, 0.86, 0.79],
            ];
            Some(BenchModel {
                cell: Cell::from_user_rows_bohr([[9., 0., 0.], [0.9, 8., 0.], [0.4, 0.8, 10.]])
                    .unwrap(),
                atoms: frac
                    .iter()
                    .enumerate()
                    .map(|(i, f)| Atom::new(if i % 2 == 0 { "Xx" } else { "Yy" }, *f))
                    .collect(),
                species: vec![("Xx".into(), p), ("Yy".into(), q)],
                options: ScfOptions {
                    cutoff_hartree: 4.,
                    max_iterations: 100,
                    ..ScfOptions::strict_reference()
                },
            })
        }
        _ => None,
    }
}

/// (SCF steps, H columns recorded + certificate + fine atomic start).
fn fine_work(o: &dft_engine::ScfOutcome) -> (usize, usize) {
    let h = &o.result.history;
    (
        h.len(),
        h.iter()
            .map(|r| {
                r.operator_applications
                    + r.fixed_point_certificate
                        .as_ref()
                        .map_or(0, |c| c.operator_applications)
            })
            .sum::<usize>()
            + o.result.atomic_start.initial_h_applications,
    )
}

fn env_list(key: &str, default: &str) -> Vec<String> {
    std::env::var(key)
        .unwrap_or_else(|_| default.into())
        .split(',')
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect()
}

#[test]
#[ignore = "release-only S1/S2 count bench (minutes on a Mac)"]
fn ladder_standard_secant_bench() {
    let cases = env_list(
        "LADDER_CASES",
        "syn-skew2,syn-slab,syn-fcc,syn-lowsym6,si-low,si-skew,si8,mgo8,si16",
    );
    let seeds: Vec<u64> = env_list("LADDER_SEEDS", "601,602,603")
        .iter()
        .map(|s| s.parse().unwrap())
        .collect();
    let amps = env_list("LADDER_AMPS", "b1e-7,f1e-3");
    let profiles = env_list("LADDER_PROFILES", "strict,standard");
    let secants: Vec<usize> = env_list("LADDER_SECANT", "2,4,6")
        .iter()
        .map(|s| s.parse().unwrap())
        .collect();
    let threads: usize = std::env::var("LADDER_THREADS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(4);
    println!(
        "case | amp | seed | profile | arm | conv | fine_scf | fine_H | coarse_scf | coarse_H_eq | total_H_eq | dE_vs_plain_Ha | dE_vs_strict_plain_Ha | seed_tol_Ha | seeded | wall_s"
    );
    for case in &cases {
        for amp in &amps {
            let (fractional, amplitude) = match amp.split_at(1) {
                ("b", v) => (false, v.parse::<f64>().unwrap()),
                ("f", v) => (true, v.parse::<f64>().unwrap()),
                _ => panic!("LADDER_AMPS entries start with b or f"),
            };
            for &seed in &seeds {
                let Some(mut model) = bench_model(case) else {
                    println!("{case} | skipped (fixture missing)");
                    break;
                };
                if seed != 0 {
                    model.perturb(seed, amplitude, fractional);
                }
                let mut strict_plain = None;
                for profile in &profiles {
                    let options = ScfOptions {
                        threads,
                        standard: (profile == "standard")
                            .then(dft_engine::StandardProfile::default),
                        ..model.options.clone()
                    };
                    let started = std::time::Instant::now();
                    let plain = model.problem(options.clone()).run_outcome();
                    let wall = started.elapsed().as_secs_f64();
                    let Ok(plain) = plain else {
                        println!("{case} | {amp} | {seed} | {profile} | plain | error");
                        continue;
                    };
                    let e_plain = plain.result.total_hartree();
                    if profile == "strict" {
                        strict_plain = Some(e_plain);
                    }
                    let (scf, h) = fine_work(&plain);
                    println!(
                        "{case} | {amp} | {seed} | {profile} | plain | {} | {scf} | {h} | 0 | 0 | {h} | 0 | {:.2e} | - | 0 | {wall:.3}",
                        plain.converged,
                        strict_plain.map_or(f64::NAN, |e| (e_plain - e).abs()),
                    );
                    for h_sec in std::iter::once(0).chain(secants.iter().copied()) {
                        let ladder = CutoffLadderOptions {
                            standard: true,
                            standard_seed: if std::env::var("LADDER_STD_SEED").as_deref()
                                == Ok("cap")
                            {
                                cutoff_ladder::StandardSeed::Cap
                            } else {
                                cutoff_ladder::StandardSeed::Coarse
                            },
                            secant_pairs: h_sec,
                            ..CutoffLadderOptions::baseline()
                        };
                        let started = std::time::Instant::now();
                        let out = cutoff_ladder::run_ladder(
                            model.problem(options.clone()),
                            &model.refs(),
                            &ladder,
                            true,
                        );
                        let wall = started.elapsed().as_secs_f64();
                        let arm = if h_sec == 0 {
                            "ladder".to_string()
                        } else {
                            format!("secant{h_sec}")
                        };
                        let Ok(out) = out else {
                            println!("{case} | {amp} | {seed} | {profile} | {arm} | error");
                            continue;
                        };
                        let r = &out.report;
                        let (scf, h) = fine_work(&out.outcome);
                        let e = out.outcome.result.total_hartree();
                        println!(
                            "{case} | {amp} | {seed} | {profile} | {arm} | {}{} | {scf} | {h} | {} | {:.0} | {:.0} | {:.2e} | {:.2e} | {} | {} | {wall:.3}",
                            out.outcome.converged,
                            if r.used {
                                String::new()
                            } else {
                                format!(
                                    " (fallback: {})",
                                    r.fallback_reason.as_deref().unwrap_or("?")
                                )
                            },
                            r.coarse_iterations,
                            r.coarse_h_equivalents(),
                            h as f64 + r.coarse_h_equivalents(),
                            (e - e_plain).abs(),
                            strict_plain.map_or(f64::NAN, |s| (e - s).abs()),
                            r.standard_seed_ha
                                .map_or("-".into(), |t| format!("{t:.2e}")),
                            r.secant_seeded,
                        );
                    }
                }
            }
        }
    }
}
