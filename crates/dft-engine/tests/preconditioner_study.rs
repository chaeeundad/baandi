//! P47A frozen-rule paired evaluation. See docs/performance/preconditioner-variants.md.
//! Run in release mode, redirecting stdout to durable evidence storage:
//! P47A_SI_UPF=/absolute/Si.pz-vbc.UPF cargo test --release -p dft-engine
//! --test preconditioner_study sweep -- --ignored --nocapture --test-threads=1
//! P47A_CASES optionally selects comma-separated benchmark IDs.
use dft_engine::DavidsonPreconditioner as Kind;
use dft_engine::{
    CompletionMode, ExecutionConfig, Functional, Mixing, OccupationPolicy, ScfOptions, ScfProblem,
    ScfSession, SpeciesPotential, StandardProfile, StepStatus,
};
use dft_geometry::{Atom, Cell, kpoints};
use std::path::PathBuf;

// Only the benchmark suite's documented flat, numeric array fields are read.
// This is not a general TOML parser or a production input adapter.
fn numbers(text: &str, key: &str) -> Option<Vec<f64>> {
    let value = text.lines().find_map(|line| {
        let (name, value) = line.split_once('=')?;
        (name.trim() == key).then_some(value)
    })?;
    Some(
        value
            .split(['[', ']', ',', ' ', '\t'])
            .filter(|s| !s.is_empty())
            .map(|s| s.parse().expect("numeric benchmark field"))
            .collect(),
    )
}

fn noise(seed: u64, len: usize) -> Vec<f64> {
    let mut state = seed;
    let mut values: Vec<_> = (0..len)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 11) as f64 / ((1_u64 << 53) as f64) - 0.5
        })
        .collect();
    let mean = values.iter().sum::<f64>() / len as f64;
    for v in &mut values {
        *v = (*v - mean) * 1e-10;
    }
    values
}

#[test]
fn perturbations_are_small_charge_neutral_and_distinct() {
    let mut previous = Vec::new();
    for seed in 401..=405 {
        let values = noise(seed, 4096);
        assert!(values.iter().all(|v| v.abs() <= 1e-10));
        assert!(values.iter().sum::<f64>().abs() < 1e-20);
        assert_ne!(values, previous);
        previous = values;
    }
}

#[test]
#[ignore = "multi-case SCF study; release only, external UPF and durable stdout required"]
fn sweep() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let pp = dft_upf::parse_upf(
        &std::fs::read_to_string(std::env::var("P47A_SI_UPF").expect("P47A_SI_UPF required"))
            .unwrap(),
    )
    .unwrap();
    let cases = std::env::var("P47A_CASES").unwrap_or_else(|_| {
        "si8-gamma,si16-gamma,si-skew,si-elongated,si-k222,si-k444s,material-si-bulk,material-mgo-bulk".into()
    });
    // Screen one selected system; do not continue past a rejection.
    for case in cases.split(',').take(1) {
        let material = case.starts_with("material-");
        let file = if material {
            format!(
                "tests/reference/materials/{}.toml",
                case.trim_start_matches("material-")
            )
        } else {
            format!("benchmarks/cases/{case}.toml")
        };
        let input = std::env::var("P47A_INPUT")
            .map(PathBuf::from)
            .unwrap_or_else(|_| root.join(file));
        let text = std::fs::read_to_string(input).unwrap();
        let labels: Vec<String> = if material {
            text.lines()
                .filter_map(|l| l.strip_prefix("species = "))
                .map(|s| s.trim_matches('"').to_owned())
                .collect()
        } else {
            vec!["Si".into()]
        };
        let mut potentials = Vec::new();
        if material {
            let dir = PathBuf::from(std::env::var("DFT_NC_DIR").expect("DFT_NC_DIR required"));
            for label in &labels {
                if !potentials.iter().any(|(name, _)| name == label) {
                    potentials.push((
                        label.clone(),
                        dft_upf::parse_upf(
                            &std::fs::read_to_string(dir.join(format!("{label}.upf"))).unwrap(),
                        )
                        .unwrap(),
                    ));
                }
            }
        }
        let species: Vec<_> = if material {
            potentials
                .iter()
                .map(|(label, potential)| SpeciesPotential { label, potential })
                .collect()
        } else {
            vec![SpeciesPotential {
                label: "Si",
                potential: &pp,
            }]
        };
        let lattice = numbers(
            &text,
            if material {
                "lattice_vectors"
            } else {
                "lattice"
            },
        )
        .unwrap();
        assert_eq!(lattice.len(), 9);
        let cell = Cell::from_user_rows_bohr(std::array::from_fn(|i| {
            std::array::from_fn(|j| lattice[3 * i + j])
        }))
        .unwrap();
        let xyz = if material {
            text.lines()
                .filter(|l| l.starts_with("fractional ="))
                .flat_map(|l| numbers(l, "fractional").unwrap())
                .collect()
        } else {
            numbers(&text, "atoms").unwrap()
        };
        assert_eq!(xyz.len() % 3, 0);
        let atoms: Vec<_> = xyz
            .as_chunks::<3>()
            .0
            .iter()
            .enumerate()
            .map(|(i, v)| Atom::new(if material { &labels[i] } else { "Si" }, *v))
            .collect();
        let grid = numbers(&text, "kgrid").unwrap_or(vec![1.; 3]);
        let shift = numbers(&text, "kshift").unwrap_or(vec![0.; 3]);
        let mesh = kpoints::reduce_time_reversal(
            &kpoints::monkhorst_pack(
                std::array::from_fn(|i| grid[i] as usize),
                std::array::from_fn(|i| shift[i]),
            )
            .unwrap(),
        )
        .unwrap();
        let mut opts = ScfOptions {
            cutoff_hartree: numbers(&text, if material { "ecutwfc_ry" } else { "cutoff_ry" })
                .unwrap()[0]
                * 0.5,
            threads: 1,
            ..Default::default()
        };
        if let Some(v) = numbers(&text, "fft_grid") {
            opts.grid_override = Some(std::array::from_fn(|i| v[i] as usize));
        }
        if let Some(v) = numbers(&text, "max_iterations") {
            opts.max_iterations = v[0] as usize;
        }
        if material {
            opts.functional = Functional::Pbe;
            opts.mixing = Mixing::KerkerPulay {
                beta: numbers(&text, "mixing_beta").unwrap()[0],
                history: 6,
                q0_inverse_bohr: numbers(&text, "mixing_q0_inverse_bohr").unwrap()[0],
            };
        }
        if std::env::var("P47A_REAL").as_deref() == Ok("on") {
            opts.gamma_real = dft_engine::GammaReal::On;
        }
        let make = |options: ScfOptions| {
            let mut p = ScfProblem::with_kpoints(
                &cell,
                &atoms,
                &species,
                options,
                &mesh,
                OccupationPolicy::Insulator {
                    minimum_gap_hartree: 1e-6,
                },
            )
            .unwrap();
            if material {
                p.set_density_floor(numbers(&text, "density_floor").unwrap()[0])
                    .unwrap();
            }
            p
        };
        // References precede variant evaluation and retain identical physical inputs.
        let mut tight = opts.clone();
        tight.density_tolerance = 1e-10;
        tight.energy_tolerance = 1e-12;
        tight.eigen_tolerance = 1e-11;
        let oracle = run_one(case, "oracle", 0, make(tight.clone()), make(tight), None);
        assert!(
            oracle.converged && oracle.certified,
            "reference must certify"
        );
        let base_strict = run_one(
            case,
            "default-strict",
            0,
            make(opts.clone()),
            make(opts.clone()),
            None,
        );
        assert!(
            base_strict.converged && base_strict.certified,
            "default strict must certify"
        );
        for variant in [Kind::Tpa, Kind::RitzShift] {
            let mut variant_opts = opts.clone();
            variant_opts.preconditioner = variant;
            let strict = run_one(
                case,
                &format!("{}-strict", variant.name()),
                0,
                make(variant_opts.clone()),
                make(variant_opts.clone()),
                None,
            );
            if !strict.converged || !strict.certified {
                println!("STOP,{case},{},strict-certification", variant.name());
                continue;
            }
            let base_error = (base_strict.energy - oracle.energy).abs();
            let error = (strict.energy - oracle.energy).abs();
            println!(
                "ACCURACY,{case},{},strict,{base_error:.17e},{error:.17e}",
                variant.name()
            );
            if error > (1.5 * base_error).max(1e-12) || strict.scf > base_strict.scf {
                println!(
                    "STOP,{case},{},strict-accuracy-or-SCF,default_scf={},variant_scf={}",
                    variant.name(),
                    base_strict.scf,
                    strict.scf
                );
                continue;
            }
            let mut hs = Vec::new();
            let mut hb = Vec::new();
            let mut rejected = false;
            for seed in 401..=405 {
                let mut base_opts = opts.clone();
                base_opts.standard = Some(StandardProfile::default());
                variant_opts.standard = base_opts.standard;
                let (base, trial) = if seed % 2 == 1 {
                    (
                        run_one(
                            case,
                            "default",
                            seed,
                            make(base_opts.clone()),
                            make(base_opts),
                            Some(seed),
                        ),
                        run_one(
                            case,
                            variant.name(),
                            seed,
                            make(variant_opts.clone()),
                            make(variant_opts.clone()),
                            Some(seed),
                        ),
                    )
                } else {
                    let trial = run_one(
                        case,
                        variant.name(),
                        seed,
                        make(variant_opts.clone()),
                        make(variant_opts.clone()),
                        Some(seed),
                    );
                    (
                        run_one(
                            case,
                            "default",
                            seed,
                            make(base_opts.clone()),
                            make(base_opts),
                            Some(seed),
                        ),
                        trial,
                    )
                };
                let be = (base.energy - oracle.energy).abs();
                let ve = (trial.energy - oracle.energy).abs();
                println!(
                    "ACCURACY,{case},{},seed={seed},{be:.17e},{ve:.17e}",
                    variant.name()
                );
                if (base.converged && !trial.converged)
                    || trial.scf > base.scf
                    || ve > 5e-9
                    || ve > 1.5 * be
                {
                    println!(
                        "STOP,{case},{},seed={seed},convergence-SCF-or-accuracy,default_scf={},variant_scf={}",
                        variant.name(),
                        base.scf,
                        trial.scf
                    );
                    rejected = true;
                    break;
                }
                hb.push(base.h);
                hs.push(trial.h);
            }
            if !rejected {
                hb.sort();
                hs.sort();
                let ratio = hs[2] as f64 / hb[2] as f64;
                println!(
                    "MEDIAN,{case},{},{},{},{ratio:.9}",
                    variant.name(),
                    hb[2],
                    hs[2]
                );
                if ratio > 1.03 {
                    println!("STOP,{case},{},H-regression", variant.name());
                }
            }
        }
    }
}

struct Measurement {
    converged: bool,
    certified: bool,
    energy: f64,
    scf: usize,
    h: usize,
}
fn run_one(
    case: &str,
    label: &str,
    seed: u64,
    problem: ScfProblem,
    source: ScfProblem,
    perturb: Option<u64>,
) -> Measurement {
    let mut density = problem.uniform_density();
    let len = density.len();
    if let Some(seed) = perturb {
        for (v, delta) in density.iter_mut().zip(noise(seed, len)) {
            *v *= 1. + delta;
        }
    }
    // Density transfer uses an independently constructed identical source problem.
    let mut session = ScfSession::new(
        problem,
        CompletionMode::RetainLast,
        ExecutionConfig::default(),
    )
    .unwrap();
    session
        .initialize_density_transfer(&source, &density)
        .unwrap();
    dft_engine::work_accounting::take();
    dft_numerics::linalg::accounting::take();
    let start = std::time::Instant::now();
    let status = loop {
        match session.step() {
            Ok(step) if step.status == StepStatus::Ready => {}
            Ok(step) => break format!("{:?}", step.status),
            Err(e) => break format!("error:{e}"),
        }
    };
    let wall = start.elapsed().as_secs_f64();
    let history = session.history();
    let h: usize = history
        .iter()
        .map(|r| {
            r.operator_applications
                + r.fixed_point_certificate
                    .as_ref()
                    .map_or(0, |c| c.operator_applications)
        })
        .sum();
    let steps: usize = history
        .iter()
        .map(|r| {
            r.davidson_iterations
                + r.fixed_point_certificate
                    .as_ref()
                    .map_or(0, |c| c.davidson_iterations)
        })
        .sum();
    let macs: u64 = dft_numerics::linalg::accounting::take()
        .iter()
        .map(|r| r.flops / 2)
        .sum();
    let energy = session
        .result()
        .map_or(f64::NAN, |r| r.result.energy.total_hartree());
    let certified = session.result().is_some_and(|r| r.result.certified);
    for row in history {
        println!(
            "ROW,{case},{label},{seed},{},{:.17e},{},{},{:.17e}",
            row.iteration,
            row.energy.total_hartree(),
            row.davidson_iterations,
            row.operator_applications,
            row.density_residual
        );
    }
    println!(
        "MEASURE,{case},{label},{seed},{status},{energy:.17e},{},{steps},{h},{macs},{wall:.9},{certified},accounting={}",
        history.len(),
        dft_numerics::linalg::accounting::ENABLED
    );
    Measurement {
        converged: status == "Converged",
        certified,
        energy,
        scf: history.len(),
        h,
    }
}
