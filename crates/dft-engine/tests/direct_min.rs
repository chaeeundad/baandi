//! Track A (direct minimization pre-converger), stages 0 and 1.
//! Fast gates run in debug; the stage-0 count bench is `#[ignore]`:
//! `cargo test --release -p dft-engine --test direct_min -- --ignored --nocapture --test-threads=1`
//! (optional `DFT_SI_UPF`, `DFT_PSEUDODOJO_DIR`, `TRACK_A_CASES=si8,syn-fcc,...`).
#[allow(dead_code)]
mod common;
use dft_engine::direct_min::{
    self, DirectMinOptions, MethodOutcome, ScfMethod, StartPath, StopReason,
};
use dft_engine::{Functional, ScfOptions, ScfOutcome, ScfProblem, SpeciesPotential};
use dft_geometry::{Atom, Cell};
use dft_upf::Pseudopotential;

struct Model {
    cell: Cell,
    atoms: Vec<Atom>,
    species: Vec<(String, Pseudopotential)>,
    options: ScfOptions,
}
impl Model {
    fn perturb(&mut self, seed: u64, bohr: f64) {
        let mut x = seed ^ 0x2545_f491_4f6c_dd1d;
        let mut uniform = || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            2. * ((x >> 11) as f64 / (1u64 << 53) as f64) - 1.
        };
        for atom in &mut self.atoms {
            let shift: [f64; 3] = std::array::from_fn(|_| bohr * uniform());
            let df = self.cell.cartesian_bohr_to_fractional(shift);
            for (f, d) in atom.fractional_position.iter_mut().zip(df) {
                *f += d;
            }
        }
    }
    fn problem(&self, threads: usize) -> ScfProblem {
        let species: Vec<_> = self
            .species
            .iter()
            .map(|(label, p)| SpeciesPotential {
                label,
                potential: p,
            })
            .collect();
        ScfProblem::new(
            &self.cell,
            &self.atoms,
            &species,
            ScfOptions {
                threads,
                ..self.options.clone()
            },
        )
        .unwrap()
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
            cutoff_hartree: 2.5,
            bands: Some(4),
            ..Default::default()
        },
    }
}

fn h_columns(o: &ScfOutcome) -> (usize, usize, usize) {
    let h = &o.result.history;
    (
        h.len(),
        h.iter().map(|r| r.operator_applications).sum(),
        h.iter()
            .filter_map(|r| r.fixed_point_certificate.as_ref())
            .map(|c| c.operator_applications)
            .sum(),
    )
}

#[test]
fn default_method_is_the_ordinary_scf() {
    let m = synthetic_skew();
    let reference = m.problem(1).run().unwrap();
    let MethodOutcome {
        outcome,
        start_path,
        minimizer,
    } = direct_min::run_with_method(
        m.problem(1),
        ScfMethod::default(),
        &DirectMinOptions::default(),
        false,
    )
    .unwrap();
    assert_eq!(start_path, StartPath::Scf);
    assert!(minimizer.is_none() && outcome.converged);
    assert_eq!(
        outcome.result.total_hartree().to_bits(),
        reference.total_hartree().to_bits()
    );
    assert_eq!(outcome.result.density, reference.density);
}

#[test]
fn minimizer_descends_to_the_scf_state_and_hands_off_certified() {
    let m = synthetic_skew();
    let reference = m.problem(1).run().unwrap();
    let tight = DirectMinOptions {
        gradient_tolerance: 1e-7,
        ..Default::default()
    };
    let mut problem = m.problem(1);
    let state = direct_min::minimize(&mut problem, &tight).unwrap();
    let r = &state.report;
    // Near 1e-7 Ha the energy decrease reaches rounding level; a failed line
    // search there is the documented, safe stop (state = last accepted point).
    assert!(
        matches!(
            r.stop,
            StopReason::GradientTolerance | StopReason::LineSearchFailed
        ) && r.projected_gradient < 1e-6,
        "{r:?}"
    );
    assert!(r.energies.windows(2).all(|w| w[1] < w[0]), "{r:?}");
    assert!(r.max_orthonormality_defect < 1e-12);
    // The minimizer's functional is the audit's direct functional (Ha).
    let de = (r.energy_hartree - reference.total_hartree()).abs();
    assert!(de < 1e-10, "minimizer energy differs by {de:e} Ha");
    let occupied = state.ritz_values.len();
    for (a, b) in state
        .ritz_values
        .iter()
        .zip(&reference.eigenvalues[..occupied])
    {
        assert!((a - b).abs() < 1e-6, "{a} {b}");
    }
    let out = direct_min::run_with_method(
        m.problem(1),
        ScfMethod::DirectThenScf,
        &DirectMinOptions::default(),
        false,
    )
    .unwrap();
    assert_eq!(out.start_path, StartPath::DirectMinimization);
    assert!(out.outcome.converged);
    let de = (out.outcome.result.total_hartree() - reference.total_hartree()).abs();
    assert!(de < 1e-9, "certified energies differ by {de:e} Ha");
}

#[test]
fn minimizer_is_thread_count_independent() {
    let m = synthetic_skew();
    let run = |threads| {
        let mut p = m.problem(threads);
        direct_min::minimize(&mut p, &DirectMinOptions::default()).unwrap()
    };
    let a = run(1);
    for threads in [2, 4] {
        let b = run(threads);
        assert_eq!(
            a.report.energy_hartree.to_bits(),
            b.report.energy_hartree.to_bits()
        );
        assert_eq!(a.report.h_columns, b.report.h_columns);
        assert!(
            a.density
                .iter()
                .zip(&b.density)
                .all(|(x, y)| x.to_bits() == y.to_bits())
        );
        assert_eq!(a.orbitals, b.orbitals);
    }
}

#[test]
fn ineligible_problems_fall_back_with_a_reason() {
    let mut m = synthetic_skew();
    m.options.gamma_real = dft_engine::GammaReal::Off;
    let out = direct_min::run_with_method(
        m.problem(1),
        ScfMethod::DirectThenScf,
        &DirectMinOptions::default(),
        false,
    )
    .unwrap();
    assert_eq!(
        out.start_path,
        StartPath::FallbackScf("gamma_real off".into())
    );
    assert!(out.minimizer.is_none() && out.outcome.converged);
    let mut p = m.problem(1);
    assert!(direct_min::minimize(&mut p, &DirectMinOptions::default()).is_err());
}

// ---------------------------------------------------------------------------
// Stage-0 bench (counts first, time second). Release only.

fn si_upf() -> Option<Pseudopotential> {
    let path = std::env::var("DFT_SI_UPF")
        .unwrap_or_else(|_| "/Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF".into());
    let text = std::fs::read_to_string(path).ok()?;
    Some(dft_upf::parse_upf(&text).unwrap())
}
fn dojo(element: &str) -> Option<Pseudopotential> {
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

fn si_model(id: &str) -> Option<Model> {
    let fcc = |shear: [f64; 2]| {
        let h = A_SI / 2.;
        [[shear[0], h, h], [h, 0., h], [h, h, shear[1]]]
    };
    let (rows, atoms, ry): ([[f64; 3]; 3], Vec<[f64; 3]>, f64) = match id {
        "si8" => (
            [[A_SI, 0., 0.], [0., A_SI, 0.], [0., 0., A_SI]],
            DIAMOND.to_vec(),
            30.,
        ),
        "si16" => (
            [[A_SI, 0., 0.], [0., A_SI, 0.], [0., 0., 2. * A_SI]],
            DIAMOND
                .iter()
                .flat_map(|p| [[p[0], p[1], p[2] / 2.], [p[0], p[1], p[2] / 2. + 0.5]])
                .collect(),
            30.,
        ),
        "si32" => (
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
        "si-skew" => (fcc([0.7, 0.3]), vec![[0.; 3], [0.25; 3]], 30.),
        "si-elongated" => (
            [[8., 0., 0.], [0., 8., 0.], [0., 0., 24.]],
            vec![
                [0., 0., 0.],
                [0.25, 0.25, 0.125],
                [0., 0., 0.5],
                [0.25, 0.25, 0.625],
            ],
            8.,
        ),
        "si-slab" => (
            [[8., 0., 0.], [0., 8., 0.], [0., 0., 30.676_713_495_509_24]],
            vec![
                [0., 0., 0.],
                [0.5, 0.5, 0.065_196_032_172_507],
                [0., 0., 0.195_588_096_517_521],
                [0.5, 0.5, 0.260_784_128_690_028],
            ],
            8.,
        ),
        _ => return None,
    };
    Some(Model {
        cell: Cell::from_user_rows_bohr(rows).unwrap(),
        atoms: atoms.into_iter().map(|p| Atom::new("Si", p)).collect(),
        species: vec![("Si".into(), si_upf()?)],
        options: ScfOptions {
            cutoff_hartree: ry / 2.,
            max_iterations: 100,
            ..Default::default()
        },
    })
}

fn other_model(id: &str) -> Option<Model> {
    match id {
        // Non-silicon synthetic multi-projector fixtures (generality rule).
        "syn-skew2" => Some(synthetic_skew()),
        "syn-fcc" => Some(Model {
            cell: common::cell(),
            atoms: common::atoms(),
            species: vec![("Xx".into(), common::potential(4.0, 1.0, true, 401))],
            options: ScfOptions {
                cutoff_hartree: 6.,
                ..Default::default()
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
            Some(Model {
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
                    ..Default::default()
                },
            })
        }
        // Real non-silicon multi-species: rocksalt MgO conventional cell, PseudoDojo
        // PBE (O has NLCC). Cutoff 35 Ha is an algorithmic test, not a converged value.
        "mgo8" => {
            let (mg, o) = (dojo("Mg")?, dojo("O")?);
            let a = 7.956;
            let mut atoms = Vec::new();
            for s in [[0., 0., 0.], [0., 0.5, 0.5], [0.5, 0., 0.5], [0.5, 0.5, 0.]] {
                atoms.push(Atom::new("Mg", s));
                atoms.push(Atom::new("O", [s[0] + 0.5, s[1], s[2]].map(|x| x % 1.)));
            }
            Some(Model {
                cell: Cell::from_user_rows_bohr([[a, 0., 0.], [0., a, 0.], [0., 0., a]]).unwrap(),
                atoms,
                species: vec![("Mg".into(), mg), ("O".into(), o)],
                options: ScfOptions {
                    cutoff_hartree: 35.,
                    functional: Functional::Pbe,
                    max_iterations: 100,
                    ..Default::default()
                },
            })
        }
        _ => None,
    }
}

#[test]
#[ignore = "release-only Track A stage-0 count bench (minutes on a Mac)"]
fn track_a_stage0_bench() {
    let default_cases =
        "syn-skew2,syn-fcc,syn-lowsym6,si-skew,si-elongated,si-slab,si8,si16,mgo8,si32";
    let cases = std::env::var("TRACK_A_CASES").unwrap_or_else(|_| default_cases.into());
    let threads: usize = std::env::var("TRACK_A_THREADS")
        .ok()
        .and_then(|t| t.parse().ok())
        .unwrap_or(1);
    println!(
        "case | occ | min_iter | min_evals | min_cols | H-eq/band | grad_Ha | stop | start | dE_min_Ha | handoff_iter | handoff_cols(hist+cert) | total_A_cols | strict_iter | strict_cols(hist+cert) | strict_atomic_cols | dE_cert_Ha | certified(handoff/strict) | ratio_A/S | wall_min_s | wall_handoff_s | wall_strict_s"
    );
    let bench_options = DirectMinOptions {
        gradient_tolerance: std::env::var("TRACK_A_GTOL")
            .ok()
            .and_then(|t| t.parse().ok())
            .unwrap_or(DirectMinOptions::default().gradient_tolerance),
        ..Default::default()
    };
    println!(
        "gradient_tolerance_ha = {:e}, threads = {threads}",
        bench_options.gradient_tolerance
    );
    // Perturbed starts as in MPI-M1-GATE3: every atom moved by a seeded
    // uniform +-1e-7 Bohr per Cartesian axis (0 = unperturbed).
    let seeds: Vec<u64> = std::env::var("TRACK_A_SEEDS")
        .map(|s| s.split(',').map(|x| x.trim().parse().unwrap()).collect())
        .unwrap_or_else(|_| vec![0]);
    let runs: Vec<(String, u64)> = cases
        .split(',')
        .flat_map(|c| seeds.iter().map(move |&s| (c.to_string(), s)))
        .collect();
    for (base, seed) in &runs {
        let label = if *seed == 0 {
            base.clone()
        } else {
            format!("{base}@{seed}")
        };
        let id = label.as_str();
        let mut model = match si_model(base).or_else(|| other_model(base)) {
            Some(m) => m,
            None => {
                println!("{id} | skipped (fixture unavailable)");
                continue;
            }
        };
        if *seed != 0 {
            model.perturb(*seed, 1e-7);
        }
        let t0 = std::time::Instant::now();
        let strict = model.problem(threads).run_outcome();
        let wall_strict = t0.elapsed().as_secs_f64();
        let (strict_e, s_iter, s_hist, s_cert, s_atomic, s_ok) = match &strict {
            Ok(o) => {
                let (i, h, c) = h_columns(o);
                (
                    o.result.total_hartree(),
                    i,
                    h,
                    c,
                    o.result.atomic_start.initial_h_applications,
                    o.converged,
                )
            }
            Err(e) => {
                println!("{id} | strict reference failed: {e}");
                (f64::NAN, 0, 0, 0, 0, false)
            }
        };
        let t0 = std::time::Instant::now();
        let mut problem = model.problem(threads);
        let state = match direct_min::minimize(&mut problem, &bench_options) {
            Ok(s) => s,
            Err(e) => {
                println!("{id} | minimizer error: {e}");
                continue;
            }
        };
        let wall_min = t0.elapsed().as_secs_f64();
        let electronic = direct_min::hand_off_state(&problem, &state).unwrap();
        let t0 = std::time::Instant::now();
        let handoff = direct_min::finish(problem, Some(&electronic), true);
        let wall_handoff = t0.elapsed().as_secs_f64();
        let r = &state.report;
        let (h_e, h_iter, h_hist, h_cert, h_ok, h_note) = match &handoff {
            Ok(o) => {
                let (i, h, c) = h_columns(o);
                (
                    o.result.total_hartree(),
                    i,
                    h,
                    c,
                    o.converged,
                    String::new(),
                )
            }
            Err(dft_engine::ScfError::Diverged { reason, history }) => (
                f64::NAN,
                history.len(),
                history.iter().map(|r| r.operator_applications).sum(),
                0,
                false,
                format!(
                    " handoff diverged: {reason}; residuals {:?}",
                    history
                        .iter()
                        .map(|r| r.density_residual)
                        .collect::<Vec<_>>()
                ),
            ),
            Err(e) => (f64::NAN, 0, 0, 0, false, format!(" handoff error: {e}")),
        };
        if std::env::var("TRACK_A_TRACE").is_ok() {
            for (label, o) in [
                ("strict", strict.as_ref().ok()),
                ("handoff", handoff.as_ref().ok()),
            ] {
                for rec in o.map_or(&[][..], |o| &o.result.history[..]) {
                    println!(
                        "  {id} {label} it={} tol={:.1e} res={:.2e} dav={} ops={} cert={:?}",
                        rec.iteration,
                        rec.requested_eigen_tolerance,
                        rec.density_residual,
                        rec.davidson_iterations,
                        rec.operator_applications,
                        rec.fixed_point_certificate.as_ref().map(|c| (
                            c.davidson_iterations,
                            c.operator_applications,
                            c.density_residual,
                            c.occupied_eigenvalue_change
                        ))
                    );
                }
            }
        }
        let hash = |d: &dft_physics::data::Array3<f64>| {
            d.iter().fold(0xcbf2_9ce4_8422_2325u64, |h, v| {
                (h ^ v.to_bits()).wrapping_mul(0x0100_0000_01b3)
            })
        };
        println!(
            "  {id} density fnv64: handoff {:016x} strict {:016x} minimizer {:016x}",
            handoff.as_ref().map_or(0, |o| hash(&o.result.density)),
            strict.as_ref().map_or(0, |o| hash(&o.result.density)),
            hash(&state.density)
        );
        let total_a = r.h_columns + h_hist + h_cert;
        let total_s = s_hist + s_cert;
        println!(
            "{id} | {} | {} | {} | {} | {:.1} | {:.2e} | {:?} | {} | {:.3e} (E_min bits {:016x}) | {} | {}+{} | {} | {} | {}+{} | {} | {:.3e} | {}/{} | {:.3} | {:.2} | {:.2} | {:.2}{h_note}",
            r.occupied_bands,
            r.iterations,
            r.energy_evaluations,
            r.h_columns,
            r.h_equivalents_per_band(),
            r.projected_gradient,
            r.stop,
            r.start,
            r.energy_hartree - strict_e,
            r.energy_hartree.to_bits(),
            h_iter,
            h_hist,
            h_cert,
            total_a,
            s_iter,
            s_hist,
            s_cert,
            s_atomic,
            h_e - strict_e,
            h_ok,
            s_ok,
            total_a as f64 / total_s as f64,
            wall_min,
            wall_handoff,
            wall_strict,
        );
    }
}
