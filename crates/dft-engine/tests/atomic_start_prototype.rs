//! P45A measurement harness; no production initialization or tolerance changes.
//! Run in release mode, redirecting stdout to durable evidence storage:
//! P45A_SI_UPF=/absolute/Si.pz-vbc.UPF cargo test --release -p dft-engine
//! --test atomic_start_prototype sweep -- --ignored --nocapture --test-threads=1
//! P45A_CASES optionally selects comma-separated benchmark IDs.
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_engine::{
    CompletionMode, ExecutionConfig, Functional, Mixing, OccupationPolicy, ScfError, ScfOptions,
    ScfProblem, ScfSession, ScfStep, SpeciesPotential, StandardProfile, StepStatus,
    density_distance,
};
use dft_geometry::{Atom, Cell, kpoints};
use dft_physics::data::*;
use dft_solvers::{
    kerker::{Kerker, MixerHistory, PreconditionedPulay},
    mixing::{LinearMixer, PulayMixer},
};
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
        &std::fs::read_to_string(std::env::var("P45A_SI_UPF").expect("P45A_SI_UPF required"))
            .unwrap(),
    )
    .unwrap();
    let cases = std::env::var("P45A_CASES")
        .unwrap_or_else(|_| "si8-gamma,si16-gamma,si-k222,material-si-20".into());
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
        let text = std::fs::read_to_string(if case == "si16-skew" {
            PathBuf::from(std::env::var("P45A_SKEW").unwrap())
        } else {
            root.join(file)
        })
        .unwrap();
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
        let mut reference: Option<Measurement> = None;
        let mut first_outputs = Vec::new();
        for run in 0..12 {
            let mut opts = ScfOptions {
                cutoff_hartree: numbers(&text, if material { "ecutwfc_ry" } else { "cutoff_ry" })
                    .unwrap()[0]
                    * 0.5,
                ..Default::default()
            };
            if material {
                opts.functional = Functional::Pbe;
                opts.threads = numbers(&text, "threads").unwrap()[0] as usize;
                opts.max_iterations = numbers(&text, "max_iterations").unwrap()[0] as usize;
                opts.mixing = Mixing::KerkerPulay {
                    beta: numbers(&text, "mixing_beta").unwrap()[0],
                    history: 6,
                    q0_inverse_bohr: numbers(&text, "mixing_q0_inverse_bohr").unwrap()[0],
                };
            }
            // A path comparison must fail closed rather than silently fall back.
            match std::env::var("P45A_GAMMA_REAL")
                .unwrap_or("off".into())
                .as_str()
            {
                "on" => opts.gamma_real = dft_engine::GammaReal::On,
                "off" => {}
                other => panic!("P45A_GAMMA_REAL must be on/off, got {other}"),
            }
            let profile = if run < 2 { "strict" } else { "standard" };
            if run < 2 {
                let scale = if run == 0 { 1.0 } else { 0.1 };
                opts.density_tolerance = 1e-9 * scale;
                opts.energy_tolerance = 1e-11 * scale;
                opts.eigen_tolerance = 1e-10 * scale;
            } else {
                opts.standard = Some(StandardProfile::default());
            }
            let arm = std::env::var("P45A_ARM").unwrap_or("default".into());
            assert!(["default", "atomic"].contains(&arm.as_str()));
            let amplitude: f64 = std::env::var("P45A_AMPLITUDE")
                .unwrap_or("1e-10".into())
                .parse()
                .unwrap();
            assert!(amplitude.is_finite() && amplitude > 0.0);
            let make_with = |options: ScfOptions| {
                let mut problem = ScfProblem::with_kpoints(
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
                    problem
                        .set_density_floor(numbers(&text, "density_floor").unwrap()[0])
                        .unwrap();
                }
                problem
            };
            let make = || make_with(opts.clone());
            let source = make();
            let mut density = source.uniform_density();
            let seed = if run < 2 { run } else { run + 199 };
            let perturbation = if run < 2 {
                vec![0.0; density.len()]
            } else {
                noise(seed, density.len())
            };
            for (v, delta) in density.iter_mut().zip(perturbation) {
                *v *= 1.0 + delta * (amplitude / 1e-10);
            }
            let mut tight_opts = opts.clone();
            tight_opts.standard = None;
            let (initial, start_h) = if arm == "atomic" {
                let raw: Vec<_> = species
                    .iter()
                    .map(|sp| {
                        let path = if material {
                            PathBuf::from(std::env::var("DFT_NC_DIR").unwrap())
                                .join(format!("{}.upf", sp.label))
                        } else {
                            PathBuf::from(std::env::var("P45A_SI_UPF").unwrap())
                        };
                        parse_chi(
                            &std::fs::read_to_string(path).unwrap(),
                            sp.potential.mesh().r().len(),
                        )
                    })
                    .collect();
                let (v, h) = atomic_initial(
                    &cell,
                    &atoms,
                    &species,
                    &raw,
                    &mesh,
                    &opts,
                    &density,
                    if material {
                        Some(numbers(&text, "density_floor").unwrap()[0])
                    } else {
                        None
                    },
                );
                (Some(v), h)
            } else {
                (None, 0)
            };
            let result = drive(
                make(),
                make_with(tight_opts),
                &cell,
                &opts,
                density.clone(),
                &arm,
                initial,
            );
            // The custom driver must reproduce the actual session, including certificates.
            if arm == "default" {
                let mut session = ScfSession::new(
                    make(),
                    CompletionMode::RetainLast,
                    ExecutionConfig::default(),
                )
                .unwrap();
                session
                    .initialize_density_transfer(&source, &density)
                    .unwrap();
                loop {
                    let step = session.step().unwrap();
                    if step.status != StepStatus::Ready {
                        break;
                    }
                }
                assert_eq!(
                    session.history().len(),
                    result.rows.len(),
                    "native length {case}/{run}"
                );
                for (native, measured) in session.history().iter().zip(&result.rows) {
                    assert_eq!(
                        native.density_residual.to_bits(),
                        measured.residual.to_bits(),
                        "native residual {case}/{run}/{}",
                        native.iteration
                    );
                    assert_eq!(
                        native.requested_eigen_tolerance.to_bits(),
                        measured.tolerance.to_bits()
                    );
                    assert_eq!(native.davidson_iterations, measured.steps);
                    assert_eq!(native.operator_applications, measured.work);
                    assert_eq!(
                        native.energy.total_hartree().to_bits(),
                        measured.energy.to_bits()
                    );
                }
                assert_eq!(session.result().unwrap().converged, result.converged);
            }
            let total_h: usize = result
                .rows
                .iter()
                .map(|r| r.work + r.certificate_work)
                .sum::<usize>()
                + start_h;
            println!(
                "TOTAL,{case},{arm},{run},{seed},{start_h},{total_h},{:.17e}",
                result.rows.last().unwrap().energy
            );
            if run == 0 {
                reference = Some(result.clone());
            }
            if run == 1 {
                println!(
                    "PLATEAU,{case},{arm},{:.17e}",
                    (result.rows.last().unwrap().energy
                        - reference.as_ref().unwrap().rows.last().unwrap().energy)
                        .abs()
                );
                reference = Some(result.clone());
            }
            if run >= 2 {
                let reference = reference.as_ref().unwrap();
                let de = (result.rows.last().unwrap().energy
                    - reference.rows.last().unwrap().energy)
                    .abs();
                let dn = density_distance(
                    &result.density,
                    &reference.density,
                    cell.volume_bohr3(),
                    source.electron_count(),
                );
                let dv = result
                    .eigen
                    .iter()
                    .zip(&reference.eigen)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0, f64::max);
                println!("ACCURACY,{case},{arm},{seed},{de:.17e},{dn:.17e},{dv:.17e}");
                if run < 4 {
                    first_outputs.push(result.first.clone());
                }
                if run == 3 {
                    println!(
                        "FIRST,{case},{arm},{:.17e}",
                        density_distance(
                            &first_outputs[0],
                            &first_outputs[1],
                            cell.volume_bohr3(),
                            source.electron_count()
                        )
                    );
                }
            }
            let path = if opts.gamma_real == dft_engine::GammaReal::On {
                "real"
            } else {
                "complex"
            };
            let mut output = String::new();
            for (i, r) in result.rows.iter().enumerate() {
                writeln!(output, "ROW,{case},{path},{profile},{arm},{amplitude:.1e},{run},{},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{},{},{},{},{},{:.17e},{:.17e}",
                    i+1,r.residual,r.error,r.tolerance,r.max_strict,r.max_strict/r.tolerance,r.steps,r.work,r.oracle_work,r.certificate_work,r.unconverged,r.energy,r.density_noise).unwrap();
                for (k, (steps, work, strict, unconverged)) in r.per_k.iter().enumerate() {
                    writeln!(output,"SOLVE,{case},{path},{profile},{arm},{amplitude:.1e},{run},{},{k},{steps},{work},{strict},{unconverged}",i+1).unwrap();
                }
            }
            print!("{output}");
            println!(
                "END,{case},{path},{profile},{arm},{amplitude:.1e},{run},{},{}",
                result.converged,
                result.rows.len()
            );
        }
    }
}

// Copied from convergence.rs for a test-only observation; no production metric changes.
pub(crate) struct CoulombMetric {
    fft: FftPlan3d,
    source: Array3<Complex64>,
    reciprocal: Array3<Complex64>,
    target: Array3<Complex64>,
    weights: Vec<f64>,
    residual: Vec<f64>,
    mapped: Vec<f64>,
}
impl std::fmt::Debug for CoulombMetric {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CoulombMetric")
            .field("grid", &self.fft.shape())
            .finish()
    }
}
impl CoulombMetric {
    pub(crate) fn new(cell: &Cell, cutoff: f64, dims: [usize; 3]) -> Result<Self, ScfError> {
        let fft = FftPlan3d::new(dims)?;
        let mut weights = vec![0.0; dims.iter().product()];
        let sphere = PlaneWaveBasis::new(cell, Hartree(cutoff))?.density_basis(cell)?;
        let scale = cell.volume_bohr3() / (2.0 * weights.len() as f64);
        for g in sphere.vectors() {
            if g.norm_squared > 0.0 {
                let [i, j, k] = FftGridSpec::map_miller(g.miller, dims)?;
                weights[(i * dims[1] + j) * dims[2] + k] =
                    scale * 4.0 * std::f64::consts::PI / g.norm_squared;
            }
        }
        Ok(Self {
            fft,
            source: Array3::zeros(dims),
            reciprocal: Array3::zeros(dims),
            target: Array3::zeros(dims),
            residual: vec![0.0; weights.len()],
            mapped: vec![0.0; weights.len()],
            weights,
        })
    }
    fn transform(&mut self) -> Result<(), dft_solvers::SolverError> {
        for (z, r) in self.source.iter_mut().zip(&self.residual) {
            *z = (*r).into();
        }
        self.fft
            .to_reciprocal(&self.source, &mut self.reciprocal)
            .map_err(|e| dft_solvers::SolverError::Application(e.to_string()))?;
        for (z, w) in self.reciprocal.iter_mut().zip(&self.weights) {
            *z *= w;
        }
        self.fft
            .to_real(&self.reciprocal, &mut self.target)
            .map_err(|e| dft_solvers::SolverError::Application(e.to_string()))?;
        for (out, z) in self.mapped.iter_mut().zip(&self.target) {
            *out = z.re;
        }
        Ok(())
    }
    pub(crate) fn error(
        &mut self,
        input: &Array3<f64>,
        output: &Array3<f64>,
    ) -> Result<f64, ScfError> {
        for ((r, a), b) in self.residual.iter_mut().zip(input).zip(output) {
            *r = b - a;
        }
        self.transform()?;
        let energy: f64 = self
            .residual
            .iter()
            .zip(&self.mapped)
            .map(|(a, b)| a * b)
            .sum();
        if !energy.is_finite() {
            return Err(ScfError::NonFinite);
        }
        Ok(energy.max(0.0))
    }
}

// Copied default mixer construction from session.rs.
fn make_mixer(
    cell: &Cell,
    _opts: &ScfOptions,
    dims: [usize; 3],
    mixing: &Mixing,
) -> Result<Box<dyn MixerHistory>, ScfError> {
    let mixer: Box<dyn MixerHistory> = match *mixing {
        Mixing::Linear { beta } => Box::new(LinearMixer { beta }),
        Mixing::Pulay { beta, history } => {
            let p = PulayMixer::new(beta, history)?;

            Box::new(p)
        }
        Mixing::KerkerPulay {
            beta,
            history,
            q0_inverse_bohr,
        } => {
            let mut fft = FftPlan3d::new(dims)?;
            let mut source = Array3::<Complex64>::zeros(dims);
            let mut target = source.clone();
            let transform = move |values: &mut [Complex64], inverse: bool| {
                source.as_slice_mut().unwrap().copy_from_slice(values);
                let result = if inverse {
                    fft.to_real(&source, &mut target)
                } else {
                    fft.to_reciprocal(&source, &mut target)
                };
                result.map_err(|e| dft_solvers::SolverError::Application(e.to_string()))?;
                values.copy_from_slice(target.as_slice().unwrap());
                Ok(())
            };
            let filter = Kerker::new(
                dims,
                cell.reciprocal_columns_inverse_bohr(),
                q0_inverse_bohr,
                transform,
            )?;
            let p = PreconditionedPulay::new(beta, history, dims.iter().product(), filter)?;

            Box::new(p)
        }
    };
    Ok(mixer)
}

#[derive(Default, Clone)]
struct Row {
    residual: f64,
    error: f64,
    tolerance: f64,
    max_strict: f64,
    steps: usize,
    work: usize,
    oracle_work: usize,
    certificate_work: usize,
    unconverged: usize,
    energy: f64,
    density_noise: f64,
    per_k: Vec<(usize, usize, usize, usize)>,
}
#[derive(Clone)]
struct Measurement {
    density: Array3<f64>,
    first: Array3<f64>,
    eigen: Vec<f64>,
    rows: Vec<Row>,
    converged: bool,
}
fn warm(step: &ScfStep) -> Vec<ComplexMatrix> {
    step.k_points
        .iter()
        .map(|k| k.solution.eigenvectors.clone())
        .collect()
}
fn work(step: &ScfStep) -> usize {
    step.k_points
        .iter()
        .map(|k| k.solution.operator_applications)
        .sum()
}
// Measurement-only SCF loop. Default histories are checked bit-for-bit against
// ScfSession; tight/density arms keep its mixing and convergence predicates.
fn drive(
    mut problem: ScfProblem,
    mut oracle: ScfProblem,
    cell: &Cell,
    opts: &ScfOptions,
    mut n: Array3<f64>,
    arm: &str,
    initial: Option<Vec<ComplexMatrix>>,
) -> Measurement {
    let electrons = problem.electron_count();
    let occupied = (electrons / 2.0).round() as usize;
    let volume = cell.volume_bohr3();
    let mean = electrons / volume;
    let mut mixer = make_mixer(cell, opts, problem.grid_dimensions(), &opts.mixing).unwrap();
    let uniform = problem.uniform_density();
    mixer
        .mix(
            uniform.as_slice().unwrap(),
            uniform.as_slice().unwrap(),
            mean,
        )
        .unwrap();
    mixer.reset();
    let mut metric =
        CoulombMetric::new(cell, opts.cutoff_hartree, problem.grid_dimensions()).unwrap();
    let mut previous = initial;
    let mut first = n.clone();
    let mut eigen = Vec::new();
    let mut rows: Vec<Row> = Vec::new();
    let mut best = f64::INFINITY;
    let mut growing = 0;
    for iteration in 1..=opts.max_iterations {
        let standard = opts.standard.is_some();
        let mut tolerance = if standard {
            rows.last().map_or(StandardProfile::CAP_HA, |r| {
                StandardProfile::ratchet(r.tolerance, r.error, electrons)
            })
        } else {
            rows.last()
                .map_or(1e-2, |r| 0.1 * r.residual)
                .min(1e-2)
                .max(0.1 * opts.eigen_tolerance)
        };
        let mut step = if arm == "tight" {
            tolerance = 1e-10;
            oracle
                .step_kpoints_with_certification(&n, previous.as_deref(), tolerance, true)
                .unwrap()
        } else {
            problem
                .step_kpoints_with_certification(
                    &n,
                    previous.as_deref(),
                    tolerance,
                    tolerance <= opts.eigen_tolerance,
                )
                .unwrap()
        };
        if iteration == 1 {
            first = step.density.clone();
        }
        eigen = step
            .k_points
            .iter()
            .flat_map(|k| k.solution.eigenvalues.iter().take(occupied).copied())
            .collect();
        let mut error = metric.error(&n, &step.density).unwrap();
        if arm != "tight" && standard && iteration == 1 && error < tolerance * electrons.max(1.0) {
            tolerance = StandardProfile::ratchet(tolerance, error, electrons);
            let mut redo = problem
                .step_kpoints_with_certification(&n, Some(&warm(&step)), tolerance, false)
                .unwrap();
            for (new, old) in redo.k_points.iter_mut().zip(&step.k_points) {
                new.solution.iterations += old.solution.iterations;
                new.solution.operator_applications += old.solution.operator_applications;
            }
            step = redo;
            error = metric.error(&n, &step.density).unwrap();
        }
        let mut row = Row {
            tolerance,
            error,
            ..Default::default()
        };
        if arm == "density" {
            // Identical n => identical input potential; retain loose orbitals for
            // the next ordinary solve, and count the density-only solve separately.
            let tight = oracle
                .step_kpoints_with_certification(&n, Some(&warm(&step)), 1e-10, true)
                .unwrap();
            row.oracle_work = work(&tight);
            row.density_noise = density_distance(&tight.density, &step.density, volume, electrons);
            step.density = tight.density;
            row.error = metric.error(&n, &step.density).unwrap();
        }
        row.residual = density_distance(&step.density, &n, volume, electrons);
        row.energy = step.energy.total_hartree();
        for k in &step.k_points {
            let s = &k.solution;
            let strict = s.buffer.iter().filter(|b| !**b).count();
            let unconverged = s.converged.iter().filter(|c| !**c).count();
            row.per_k
                .push((s.iterations, s.operator_applications, strict, unconverged));
            row.steps += s.iterations;
            row.work += s.operator_applications;
            row.unconverged += unconverged;
            row.max_strict = s
                .residuals
                .iter()
                .zip(&s.buffer)
                .filter(|(_, b)| !**b)
                .map(|(r, _)| *r)
                .fold(row.max_strict, f64::max);
        }
        let occupied_residual = step
            .k_points
            .iter()
            .flat_map(|k| k.solution.residuals.iter().take(occupied))
            .copied()
            .fold(0.0, f64::max);
        let charge_error =
            (step.density.sum() * volume / n.len() as f64 - electrons).abs() / electrons;
        let input_charge_error = (n.sum() * volume / n.len() as f64 - electrons).abs() / electrons;
        let converged = if standard {
            row.error < opts.standard.unwrap().threshold_ha
        } else {
            row.residual <= opts.density_tolerance
                && rows
                    .last()
                    .is_some_and(|r| (row.energy - r.energy).abs() <= opts.energy_tolerance)
                && occupied_residual <= opts.eigen_tolerance
                && tolerance <= opts.eigen_tolerance
                && charge_error <= 1e-8
                && input_charge_error <= 1e-8
        };
        growing = if row.residual > opts.divergence_factor * best.max(opts.density_tolerance) {
            growing + 1
        } else {
            0
        };
        best = best.min(row.residual);
        assert!(
            growing < opts.divergence_patience,
            "diverged at {iteration}"
        );
        let mut certified = converged;
        if converged && !standard {
            let check = problem
                .step_kpoints(
                    &step.density,
                    Some(&warm(&step)),
                    if arm == "tight" {
                        1e-10
                    } else {
                        0.1 * opts.eigen_tolerance
                    },
                )
                .unwrap();
            row.certificate_work = work(&check);
            let dr = density_distance(&check.density, &step.density, volume, electrons);
            let de = check
                .k_points
                .iter()
                .zip(&step.k_points)
                .flat_map(|(a, b)| {
                    a.solution
                        .eigenvalues
                        .iter()
                        .zip(&b.solution.eigenvalues)
                        .take(occupied)
                })
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max);
            certified = dr <= opts.density_tolerance && de <= 1e-7;
        }
        rows.push(row);
        if certified {
            return Measurement {
                density: step.density,
                first,
                eigen,
                rows,
                converged: true,
            };
        }
        let mut mixed = mixer
            .mix(
                n.as_slice().unwrap(),
                step.density.as_slice().unwrap(),
                mean,
            )
            .unwrap();
        let fraction = n
            .iter()
            .zip(&mixed)
            .filter(|(_, b)| **b < 0.0)
            .map(|(a, b)| 0.9 * a / (a - b))
            .fold(1.0, f64::min);
        if fraction < 1.0 {
            for (b, a) in mixed.iter_mut().zip(n.iter()) {
                *b = a + fraction * (*b - a);
            }
            mixer.reset();
        }
        n.as_slice_mut().unwrap().copy_from_slice(&mixed);
        previous = Some(warm(&step));
    }
    Measurement {
        density: n,
        first,
        eigen,
        rows,
        converged: false,
    }
}

#[derive(Debug)]
struct Chi {
    l: usize,
    occupation: f64,
    samples: Vec<f64>,
}
// Deliberately narrow UPF-v2 text scan, confined to PP_PSWFC. No parser changes.
fn parse_chi(raw: &str, mesh_len: usize) -> Vec<Chi> {
    fn attr<'a>(header: &'a str, name: &str) -> &'a str {
        let marker = format!("{name}=\"");
        let (offset, _) = header
            .match_indices(&marker)
            .find(|(i, _)| *i == 0 || header.as_bytes()[i - 1].is_ascii_whitespace())
            .expect("missing CHI attribute");
        header[offset + marker.len()..]
            .split('"')
            .next()
            .unwrap()
            .trim()
    }
    let section = raw
        .split_once("<PP_PSWFC>")
        .expect("PP_PSWFC required")
        .1
        .split_once("</PP_PSWFC>")
        .expect("closed PP_PSWFC")
        .0;
    let mut out = Vec::new();
    let mut rest = section;
    while let Some((_, tail)) = rest.split_once("<PP_CHI.") {
        let (header, body) = tail.split_once('>').expect("CHI header");
        let index = header.split_whitespace().next().unwrap();
        let closing = format!("</PP_CHI.{index}>");
        let (values, next) = body.split_once(&closing).expect("matching CHI close");
        let samples: Vec<f64> = values
            .split_whitespace()
            .map(|s| s.replace(['D', 'd'], "E").parse().unwrap())
            .collect();
        let l = attr(header, "l").parse().unwrap();
        let occupation: f64 = attr(header, "occupation").parse().unwrap();
        assert!(l <= 3 && occupation.is_finite() && occupation >= 0.0);
        assert_eq!(samples.len(), mesh_len, "CHI must cover radial mesh");
        assert!(samples.iter().all(|v| v.is_finite()));
        out.push(Chi {
            l,
            occupation,
            samples,
        });
        rest = next;
    }
    assert!(!out.is_empty(), "no atomic orbitals");
    out
}

// All operations here occur once before the SCF driver. The first potential is
// assembled in the same order and with the same public physics functions as scf.rs.
#[allow(clippy::too_many_arguments)]
fn atomic_initial(
    cell: &Cell,
    atoms: &[Atom],
    potentials: &[SpeciesPotential<'_>],
    chi: &[Vec<Chi>],
    points: &[kpoints::KPoint],
    opts: &ScfOptions,
    density: &Array3<f64>,
    floor: Option<f64>,
) -> (Vec<ComplexMatrix>, usize) {
    use dft_numerics::linalg::{MatrixOp, gemm, hermitian_eigh, orthonormalize};
    use dft_physics::{
        hamiltonian::{self, Hamiltonian, NonlocalOperator, Species},
        hartree,
    };
    use dft_solvers::LinearOperator;
    use dft_upf::radial::{projector_transform_with_rule, real_spherical_harmonic};
    let positions: Vec<Vec<_>> = potentials
        .iter()
        .map(|p| {
            atoms
                .iter()
                .filter(|a| a.species == p.label)
                .map(|a| {
                    cell.fractional_to_cartesian_bohr(
                        a.fractional_position.map(|r| r.rem_euclid(1.0)),
                    )
                })
                .collect()
        })
        .collect();
    let species: Vec<_> = potentials
        .iter()
        .zip(&positions)
        .map(|(p, pos)| Species {
            potential: p.potential,
            positions_bohr: pos,
        })
        .collect();
    let db = PlaneWaveBasis::at_k(cell, Hartree(4.0 * opts.cutoff_hartree), [0.0; 3]).unwrap();
    let dims = std::array::from_fn(|i| density.shape()[i]);
    let dg = FftGridSpec::new(&db, &db, Some(dims)).unwrap();
    let ig = dg
        .scatter(
            &hamiltonian::ionic_local_with_conventions(
                &db,
                cell.volume_bohr3(),
                &species,
                opts.local_mesh,
                opts.radial_even_rule,
            )
            .unwrap(),
        )
        .unwrap();
    let mut fft = FftPlan3d::new(dims).unwrap();
    let core = core_density_with_rule(
        &db,
        &dg,
        cell.volume_bohr3(),
        &species,
        &mut fft,
        opts.radial_even_rule,
    )
    .unwrap();
    let mut ionic = Array3::zeros(dims);
    fft.to_real(&ig, &mut ionic).unwrap();
    let h = hartree::evaluate(density, cell, &mut fft).unwrap();
    let mut scratch = Array3::zeros(dims);
    let x = opts
        .functional
        .evaluate_with_floor(density, core.as_ref(), &mut scratch, cell, &mut fft, floor)
        .unwrap();
    let mut effective = ionic.mapv(|z| z.re);
    for ((v, h), x) in effective.iter_mut().zip(&h.potential).zip(&x.potential) {
        *v += h.re + x;
    }
    let electrons: f64 = species
        .iter()
        .map(|s| s.potential.z_valence() * s.positions_bohr.len() as f64)
        .sum();
    let bands = opts
        .bands
        .unwrap_or(electrons as usize / 2 + opts.extra_bands);
    let mut guesses = Vec::new();
    let mut start_h = 0;
    for point in points {
        let basis =
            PlaneWaveBasis::at_k(cell, Hartree(opts.cutoff_hartree), point.fractional).unwrap();
        let npw = basis.vectors().len();
        let count: usize = species
            .iter()
            .zip(chi)
            .map(|(s, cs)| s.positions_bohr.len() * cs.iter().map(|c| 2 * c.l + 1).sum::<usize>())
            .sum();
        let mut raw = ComplexMatrix::zeros(npw, count).unwrap();
        let mut col = 0;
        let mut radial_samples = 0usize;
        for (sp, cs) in species.iter().zip(chi) {
            for c in cs {
                // Include every supplied orbital, even if occupation is zero.
                assert!(c.occupation >= 0.0);
                let mut cache = std::collections::BTreeMap::new();
                let radial: Vec<_> = basis
                    .momenta()
                    .iter()
                    .map(|g| {
                        *cache.entry(g.norm_squared.to_bits()).or_insert_with(|| {
                            radial_samples += c.samples.len();
                            projector_transform_with_rule(
                                sp.potential.mesh(),
                                c.l,
                                &c.samples,
                                g.norm_squared.sqrt(),
                                opts.radial_even_rule,
                            )
                            .unwrap()
                                / cell.volume_bohr3().sqrt()
                        })
                    })
                    .collect();
                let il = Complex64::new(0.0, -1.0).powu(c.l as u32);
                for &pos in sp.positions_bohr {
                    for m in -(c.l as isize)..=c.l as isize {
                        for (i, g) in basis.momenta().iter().enumerate() {
                            let y = if g.norm_squared == 0.0 {
                                if c.l == 0 {
                                    1.0 / (4.0 * std::f64::consts::PI).sqrt()
                                } else {
                                    0.0
                                }
                            } else {
                                real_spherical_harmonic(c.l, m, g.cartesian).unwrap()
                            };
                            let phase = hamiltonian::structure_factor(g.cartesian, &[pos]).unwrap();
                            raw[(i, col)] = il * phase * (radial[i] * y);
                        }
                        col += 1;
                    }
                }
            }
        }
        let ortho = orthonormalize(&raw, None).unwrap();
        assert!(ortho.dropped.is_empty(), "rank deficient atomic orbitals");
        let mut q = ortho.vectors;
        // No noise is needed when the atomic span covers all requested bands.
        // Otherwise append lowest kinetic-energy waves, rejecting dependent ones.
        if q.cols() < bands {
            let mut order: Vec<_> = (0..npw).collect();
            order.sort_by(|&a, &b| {
                basis.momenta()[a]
                    .norm_squared
                    .total_cmp(&basis.momenta()[b].norm_squared)
            });
            for i in order {
                let mut e = ComplexMatrix::zeros(npw, 1).unwrap();
                e[(i, 0)] = Complex64::new(1.0, 0.0);
                let extra = orthonormalize(&e, Some(&q)).unwrap();
                if extra.vectors.cols() == 1 {
                    q.append_columns(&extra.vectors).unwrap();
                }
                if q.cols() == bands {
                    break;
                }
            }
        }
        assert!(q.cols() >= bands, "insufficient starting rank");
        let ham = Hamiltonian::new(
            cell,
            &basis,
            effective.clone(),
            NonlocalOperator::new_with_rule(
                &basis,
                cell.volume_bohr3(),
                &species,
                opts.radial_even_rule,
            )
            .unwrap(),
        )
        .unwrap();
        let mut hq = ComplexMatrix::zeros(npw, q.cols()).unwrap();
        ham.apply_block(&q, &mut hq).unwrap();
        start_h += q.cols();
        let mut small = gemm(&q, MatrixOp::Adjoint, &hq, MatrixOp::Identity).unwrap();
        for j in 0..small.cols() {
            for i in 0..=j {
                let v = (small[(i, j)] + small[(j, i)].conj()) * 0.5;
                small[(i, j)] = v;
                small[(j, i)] = v.conj();
            }
        }
        // Gamma atomic orbitals represent real functions. Remove only the
        // roundoff imaginary reduced entries before EVD so degenerate Ritz
        // vectors retain conjugate symmetry required by the real bridge.
        if point.fractional == [0.0; 3] {
            for v in small.as_mut_slice() {
                assert!(v.im.abs() < 1e-12, "Gamma reduced Hamiltonian must be real");
                v.im = 0.0;
            }
        }
        let (_, rot) = hermitian_eigh(&small).unwrap();
        let mut selected = ComplexMatrix::zeros(q.cols(), bands).unwrap();
        for j in 0..bands {
            selected.column_mut(j).copy_from_slice(rot.column(j));
        }
        let guess = gemm(&q, MatrixOp::Identity, &selected, MatrixOp::Identity).unwrap();
        println!(
            "START,k={:?},atomic={count},rank={},H={},radial_samples={radial_samples},gram_rotation_complex_mac={},evd_order={}",
            point.fractional,
            q.cols(),
            q.cols(),
            npw * q.cols() * (q.cols() + bands),
            q.cols()
        );
        guesses.push(guess);
    }
    (guesses, start_h)
}

#[test]
fn radial_gaussian_and_chi_contract() {
    use dft_upf::radial::{RadialMesh, projector_transform};
    // Fourier transform of exp(-a*r^2) is (pi/a)^(3/2) exp(-q^2/(4a)).
    let a: f64 = 0.7;
    let dr = 0.005;
    let r: Vec<_> = (0..2001).map(|i| i as f64 * dr).collect();
    let mesh = RadialMesh::new(r.clone(), vec![dr; r.len()]).unwrap();
    let chi: Vec<_> = r.iter().map(|r| r * (-a * r * r).exp()).collect();
    for q in [0.0, 0.5, 2.0] {
        let actual = projector_transform(&mesh, 0, &chi, q).unwrap();
        let expected = (std::f64::consts::PI / a).powf(1.5) * (-q * q / (4.0 * a)).exp();
        assert!((actual - expected).abs() < 1e-10);
    }
    let parsed = parse_chi(
        "<PP_PSWFC><PP_CHI.1 label=\"3P\" l=\"1\" occupation=\" 2.0\">0 1D-2 0</PP_CHI.1></PP_PSWFC>",
        3,
    );
    assert_eq!(parsed[0].l, 1);
    assert_eq!(parsed[0].samples, [0.0, 0.01, 0.0]);
}
