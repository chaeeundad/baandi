//! P42A boundary distribution and frozen-density spectrum harness; no production initialization or tolerance changes.
//! Run in release mode, redirecting stdout to durable evidence storage:
//! P42A_SI_UPF=/absolute/Si.pz-vbc.UPF cargo test --release -p dft-engine
//! --test band_boundary_scan sweep -- --ignored --nocapture --test-threads=1
//! P42A_CASES optionally selects comma-separated benchmark IDs.
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
        &std::fs::read_to_string(std::env::var("P42A_SI_UPF").expect("P42A_SI_UPF required"))
            .unwrap(),
    )
    .unwrap();
    let cases = std::env::var("P42A_CASES").unwrap_or_else(|_| {
        "si-gamma,si-k222s,si-skew,si-elongated,si-slab,si8-gamma,si16-gamma,si-low".into()
    });
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
        // Two strict reference levels, then ten independent density perturbations.
        // A failed reference is retained as a failure, never an accuracy oracle.
        let first = std::env::var("P42A_START_RUN")
            .unwrap_or("0".into())
            .parse::<u64>()
            .unwrap();
        for run in first..12 {
            let mut opts = ScfOptions {
                cutoff_hartree: numbers(&text, if material { "ecutwfc_ry" } else { "cutoff_ry" })
                    .unwrap()[0]
                    * 0.5,
                ..Default::default()
            };
            opts.threads = 1;
            if let Some(v) = numbers(&text, "fft_grid") {
                opts.grid_override = Some(std::array::from_fn(|i| v[i] as usize));
            }
            if let Some(v) = numbers(&text, "max_iterations") {
                opts.max_iterations = v[0] as usize;
            }
            if material {
                opts.functional = Functional::Pbe;

                opts.max_iterations = numbers(&text, "max_iterations").unwrap()[0] as usize;
                opts.mixing = Mixing::KerkerPulay {
                    beta: numbers(&text, "mixing_beta").unwrap()[0],
                    history: 6,
                    q0_inverse_bohr: numbers(&text, "mixing_q0_inverse_bohr").unwrap()[0],
                };
            }
            // Optional number of extra bands above the occupied ones (default: engine default).
            if let Ok(v) = std::env::var("P42A_EXTRA_BANDS") {
                opts.extra_bands = v.parse().expect("P42A_EXTRA_BANDS must be an integer");
            }
            // Optional real Gamma path (default off: the study was written for the complex
            // path). `P42A_GAMMA_REAL=on` fails closed on ineligible cases, `auto` falls back.
            match std::env::var("P42A_GAMMA_REAL").as_deref() {
                Ok("on") => opts.gamma_real = dft_engine::GammaReal::On,
                Ok("auto") => opts.gamma_real = dft_engine::GammaReal::Auto,
                _ => {}
            }
            if run < 2 {
                let scale = if run == 0 { 1.0 } else { 0.1 };
                opts.density_tolerance = 1e-9 * scale;
                opts.energy_tolerance = 1e-11 * scale;
                opts.eigen_tolerance = 1e-10 * scale;
            } else {
                opts.standard = Some(StandardProfile::default());
                opts.standard_stall_guard = std::env::var("P42A_GUARD").is_ok_and(|v| v == "on");
            }
            let make = |options: ScfOptions| {
                let mut problem = ScfProblem::with_kpoints(
                    &cell,
                    &atoms,
                    &species,
                    options,
                    &mesh,
                    OccupationPolicy::Insulator {
                        minimum_gap_hartree: 1e-6,
                    },
                )?;
                if material {
                    problem
                        .set_density_floor(numbers(&text, "density_floor").unwrap()[0])
                        .unwrap();
                }
                Ok::<_, dft_engine::ScfError>(problem)
            };
            let source = match make(opts.clone()) {
                Ok(p) => p,
                Err(e) => {
                    println!("END,{case},{run},setup-error:{e}");
                    continue;
                }
            };
            let mut density = source.uniform_density();
            if run >= 2 {
                let perturbation = noise(run - 1, density.len());
                for (v, delta) in density.iter_mut().zip(perturbation) {
                    *v *= 1.0 + delta;
                }
            }
            let mut session = ScfSession::new(
                make(opts.clone()).unwrap(),
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
            println!("END,{case},{run},{status}");
            let history = failed_history
                .as_deref()
                .unwrap_or_else(|| session.history());
            let cert_h: usize = history
                .iter()
                .filter_map(|r| r.fixed_point_certificate.as_ref())
                .map(|c| c.operator_applications)
                .sum();
            println!("CERT,{case},{run},{cert_h}");
            if let Some(outcome) = session.result() {
                if run >= 2 {
                    assert!(!outcome.result.certified);
                }
                println!("GAP,{case},{run},{:?}", outcome.result.gap_hartree);
                if status != "Converged" {
                    continue;
                }
                let result = &outcome.result;
                let occupied = result.occupations.iter().filter(|&&v| v > 0.).count();
                println!(
                    "RESULT,{case},{run},{:.17e},{},{occupied},{}",
                    result.total_hartree(),
                    result.certified,
                    result.gamma_real_used
                );
                for (k, point) in result.k_points.iter().enumerate() {
                    println!("EIG|{case}|{run}|{k}|{:?}", point.solution.eigenvalues);
                }
                let out = PathBuf::from(std::env::var("P42A_OUT").expect("P42A_OUT required"));
                std::fs::create_dir_all(&out).unwrap();
                let bytes: Vec<_> = result
                    .density
                    .iter()
                    .flat_map(|v| v.to_le_bytes())
                    .collect();
                std::fs::write(out.join(format!("density-{run}.bin")), bytes).unwrap();
                if run == 1 {
                    // Independent diagnostic at the tight reference density; never feeds SCF.
                    let n = result.eigenvalues.len();
                    let mut probe_opts = opts.clone();
                    probe_opts.extra_bands += 1;
                    probe_opts.eigen_tolerance = 1e-5;
                    match make(probe_opts)
                        .and_then(|mut p| p.step_kpoints(&result.density, None, 1e-5))
                    {
                        Ok(probe) => {
                            for (k, point) in probe.k_points.iter().enumerate() {
                                println!(
                                    "SPECTRUM|{case}|{k}|{occupied}|{n}|{:?}|{:?}|{}|{}",
                                    point.solution.eigenvalues,
                                    point.solution.residuals,
                                    point.solution.residuals_certified,
                                    point.solution.operator_applications
                                );
                            }
                        }
                        Err(e) => println!("PROBE_ERROR|{case}|{e}"),
                    }
                }
            }
        }
    }
}
