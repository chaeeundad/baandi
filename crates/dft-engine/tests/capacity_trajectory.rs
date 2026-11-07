//! P31A measurement harness; no production initialization or tolerance changes.
//! Run in release mode, redirecting stdout to durable evidence storage:
//! P31A_SI_UPF=/absolute/Si.pz-vbc.UPF cargo test --release -p dft-engine
//! --test standard_stall_trajectory sweep -- --ignored --nocapture --test-threads=1
//! P31A_CASES selects comma-separated benchmark IDs (or material-si-20/material-mgo-40).
//! P31A_CAPACITY=2..12, P31A_PROFILE=strict|standard, P31A_REAL=off|on.
//! Runs 0/1 are tight references; 2..6 are five deterministic perturbations.
//! P31A_START_RUN/P31A_END_RUN select a half-open range; material cases need DFT_NC_DIR.
//! The full experiment commands are retained in P31A-evidence/study.py.
use dft_engine::{
    CompletionMode, ExecutionConfig, Functional, Mixing, OccupationPolicy, ScfOptions, ScfProblem,
    ScfSession, SpeciesPotential, StandardProfile, StepStatus,
};
use dft_geometry::{Atom, Cell, kpoints};
use std::{fmt::Write, path::PathBuf};

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
    for seed in 1..=10 {
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
        &std::fs::read_to_string(std::env::var("P31A_SI_UPF").expect("P31A_SI_UPF required"))
            .unwrap(),
    )
    .unwrap();
    let cases = std::env::var("P31A_CASES")
        .unwrap_or_else(|_| "si-gamma,si-k222s,si8-gamma,si16-gamma,si32-gamma".into());
    for case in cases.split(',') {
        let material = case.starts_with("material-");
        let file = if material {
            format!(
                "tests/reference/materials/{}.toml",
                case.trim_start_matches("material-")
            )
        } else {
            format!("benchmarks/cases/{case}.toml")
        };
        let text = std::fs::read_to_string(root.join(file)).unwrap();
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
        // Two strict reference levels, then independent density perturbations.
        // A failed reference is retained as a failure, never an accuracy oracle.
        let first = std::env::var("P31A_START_RUN")
            .unwrap_or("0".into())
            .parse::<u64>()
            .unwrap();
        let last = std::env::var("P31A_END_RUN")
            .unwrap_or("7".into())
            .parse::<u64>()
            .unwrap();
        let multiplier = std::env::var("P31A_CAPACITY")
            .unwrap_or("6".into())
            .parse::<usize>()
            .unwrap();
        let profile = std::env::var("P31A_PROFILE").unwrap_or("standard".into());
        let real = std::env::var("P31A_REAL").is_ok_and(|v| v == "on");
        for run in first..last {
            let started = std::time::Instant::now();
            let mut opts = ScfOptions {
                davidson_capacity_per_band: multiplier,
                gamma_real: if real {
                    dft_engine::GammaReal::On
                } else {
                    dft_engine::GammaReal::Off
                },
                cutoff_hartree: numbers(&text, if material { "ecutwfc_ry" } else { "cutoff_ry" })
                    .unwrap()[0]
                    * 0.5,
                ..Default::default()
            };
            if material {
                opts.functional = Functional::Pbe;
                opts.threads = 1;
                opts.max_iterations = numbers(&text, "max_iterations").unwrap()[0] as usize;
                opts.mixing = Mixing::KerkerPulay {
                    beta: numbers(&text, "mixing_beta").unwrap()[0],
                    history: 6,
                    q0_inverse_bohr: numbers(&text, "mixing_q0_inverse_bohr").unwrap()[0],
                };
            }
            if run < 2 {
                let scale = if run == 0 { 1.0 } else { 0.1 };
                opts.density_tolerance = 1e-9 * scale;
                opts.energy_tolerance = 1e-11 * scale;
                opts.eigen_tolerance = 1e-10 * scale;
            } else if profile == "standard" {
                opts.standard = Some(StandardProfile::default());
                opts.standard_stall_guard = std::env::var("P31A_GUARD").is_ok_and(|v| v == "on");
            }
            let make = || {
                let mut problem = ScfProblem::with_kpoints(
                    &cell,
                    &atoms,
                    &species,
                    opts.clone(),
                    &mesh,
                    OccupationPolicy::Insulator {
                        minimum_gap_hartree: 1e-6,
                    },
                )
                .unwrap();
                if material {
                    problem
                        .set_density_floor(numbers(&text, "density_floor").unwrap()[0])
                        .unwrap();
                }
                problem
            };
            let source = make();
            let mut density = source.uniform_density();
            if run >= 2 {
                let perturbation = noise(run - 1, density.len());
                for (v, delta) in density.iter_mut().zip(perturbation) {
                    *v *= 1.0 + delta;
                }
            }
            let analysis = dft_numerics::pool::analysis::Span::new("scf");
            let mut session = ScfSession::new(
                make(),
                CompletionMode::RetainLast,
                ExecutionConfig::default(),
            )
            .unwrap();
            session
                .initialize_density_transfer(&source, &density)
                .unwrap();
            let mut failed_history = None;
            let status = loop {
                match session.step() {
                    Ok(step) => {
                        if step.status != StepStatus::Ready {
                            break format!("{:?}", step.status);
                        }
                    }
                    Err(e) => {
                        if let dft_engine::ScfError::Diverged { history, .. }
                        | dft_engine::ScfError::NotConverged { history } = &e
                        {
                            failed_history = Some(history.clone());
                        }
                        break format!(
                            "error:{}",
                            e.to_string().split(", history:").next().unwrap()
                        );
                    }
                }
            };
            let wall = started.elapsed().as_secs_f64();
            drop(analysis);
            for record in dft_numerics::linalg::accounting::take() {
                println!(
                    "GEMM,{case},{run},{},{},{},{},{},{},{},{}",
                    record.solve,
                    record.step,
                    record.site,
                    record.m,
                    record.n,
                    record.k,
                    record.calls,
                    record.nanoseconds
                );
            }
            let work = dft_engine::work_accounting::take();
            println!(
                "WORK,{case},{run},{},{},{},{}",
                work.iter().map(|w| w.iterations).sum::<usize>(),
                work.iter().map(|w| w.vector_applications).sum::<usize>(),
                work.iter().map(|w| w.restarts).sum::<usize>(),
                work.iter().filter(|w| !w.returned).count()
            );
            let restarts = work.iter().map(|w| w.restarts).sum::<usize>();
            println!(
                "SUMMARY,{case},{run},{multiplier},{profile},{real},{wall:.9},{restarts},{}",
                dft_engine::work_accounting::ENABLED
            );
            let mut output = String::new();
            for r in failed_history
                .as_deref()
                .unwrap_or_else(|| session.history())
            {
                writeln!(
                    output,
                    "ROW,{case},{run},{},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{},{}",
                    r.iteration,
                    r.energy.total_hartree(),
                    r.density_residual,
                    r.standard_error_ha.unwrap_or(f64::NAN),
                    r.requested_eigen_tolerance,
                    r.max_eigen_residual,
                    r.davidson_iterations,
                    r.operator_applications
                )
                .unwrap();
            }
            // I/O is outside the iteration loop; includes the terminal failed record.
            print!("{output}");
            if let Some(outcome) = session.result() {
                println!(
                    "RESULT,{case},{run},{:.17e},{},{}",
                    outcome.result.energy.total_hartree(),
                    outcome.result.certified,
                    outcome.result.history.len()
                );
            }
            let certificates: Vec<_> = session
                .history()
                .iter()
                .filter_map(|r| r.fixed_point_certificate.as_ref())
                .collect();
            println!(
                "CERT,{case},{run},{},{},{:.17e},{:.17e}",
                certificates
                    .iter()
                    .map(|c| c.davidson_iterations)
                    .sum::<usize>(),
                certificates
                    .iter()
                    .map(|c| c.operator_applications)
                    .sum::<usize>(),
                certificates.last().map_or(f64::NAN, |c| c.density_residual),
                certificates
                    .last()
                    .map_or(f64::NAN, |c| c.occupied_eigenvalue_change)
            );
            println!("END,{case},{run},{status}");
            if let Some(outcome) = session.result() {
                if run >= 2 && profile == "standard" {
                    assert!(!outcome.result.certified);
                }
                println!("GAP,{case},{run},{:?}", outcome.result.gap_hartree);
            }
        }
    }
}
