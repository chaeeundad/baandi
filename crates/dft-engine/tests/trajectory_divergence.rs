//! P41A measurement harness; no production initialization or tolerance changes.
//! Run in release mode, redirecting stdout to durable evidence storage:
//! P41A_SI_UPF=/absolute/Si.pz-vbc.UPF cargo test --release -p dft-engine
//! --test trajectory_divergence sweep -- --ignored --nocapture --test-threads=1
//! P41A_CASES optionally selects comma-separated benchmark IDs.
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
        &std::fs::read_to_string(std::env::var("P41A_SI_UPF").expect("P41A_SI_UPF required"))
            .unwrap(),
    )
    .unwrap();
    let cases = std::env::var("P41A_CASES")
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
        for run in 1..=10 {
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
            match std::env::var("P41A_GAMMA_REAL")
                .unwrap_or("off".into())
                .as_str()
            {
                "on" => opts.gamma_real = dft_engine::GammaReal::On,
                "off" => {}
                other => panic!("P41A_GAMMA_REAL must be on/off, got {other}"),
            }
            let profile = std::env::var("P41A_PROFILE").unwrap_or("standard".into());
            assert!(["strict", "standard"].contains(&profile.as_str()));
            if profile == "standard" {
                opts.standard = Some(StandardProfile::default());
            }
            let arm = std::env::var("P41A_ARM").unwrap_or("default".into());
            assert!(["default", "tight", "density"].contains(&arm.as_str()));
            let amplitude: f64 = std::env::var("P41A_AMPLITUDE")
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
            let perturbation = noise(run, density.len());
            for (v, delta) in density.iter_mut().zip(perturbation) {
                *v *= 1.0 + delta * (amplitude / 1e-10);
            }
            let mut tight_opts = opts.clone();
            tight_opts.standard = None;
            let result = drive(
                make(),
                make_with(tight_opts),
                &cell,
                &opts,
                density.clone(),
                &arm,
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

#[derive(Default)]
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
struct Measurement {
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
    let mut previous = None;
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
        rows,
        converged: false,
    }
}

// Stable projection defect, avoiding subtraction of two nearly equal traces.
fn occupied_subspace_defect(a: &ComplexMatrix, b: &ComplexMatrix, occupied: usize) -> f64 {
    let mut squared = 0.0;
    for j in 0..occupied {
        let overlaps: Vec<_> = (0..occupied)
            .map(|k| {
                a.column(k)
                    .iter()
                    .zip(b.column(j))
                    .map(|(u, v)| u.conj() * v)
                    .sum::<Complex64>()
            })
            .collect();
        for i in 0..a.rows() {
            let projected: Complex64 = overlaps
                .iter()
                .enumerate()
                .map(|(k, c)| a.column(k)[i] * c)
                .sum();
            squared += (b.column(j)[i] - projected).norm_sqr();
        }
    }
    squared.sqrt()
}

/// First-map diagnostic, separate from the full ensemble. Same release command
/// as sweep, replacing its test filter with `first_map`.
#[test]
#[ignore = "first-map subspace diagnostic; release only and external UPF required"]
fn first_map() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let pp = dft_upf::parse_upf(
        &std::fs::read_to_string(std::env::var("P41A_SI_UPF").unwrap()).unwrap(),
    )
    .unwrap();
    for case in ["si8-gamma", "si16-gamma"] {
        let text =
            std::fs::read_to_string(root.join(format!("benchmarks/cases/{case}.toml"))).unwrap();
        let lattice = numbers(&text, "lattice").unwrap();
        let cell = Cell::from_user_rows_bohr(std::array::from_fn(|i| {
            std::array::from_fn(|j| lattice[3 * i + j])
        }))
        .unwrap();
        let xyz = numbers(&text, "atoms").unwrap();
        let atoms: Vec<_> = xyz
            .as_chunks::<3>()
            .0
            .iter()
            .map(|v| Atom::new("Si", *v))
            .collect();
        let species = [SpeciesPotential {
            label: "Si",
            potential: &pp,
        }];
        for real in [false, true] {
            for standard in [false, true] {
                let options = ScfOptions {
                    cutoff_hartree: 0.5 * numbers(&text, "cutoff_ry").unwrap()[0],
                    gamma_real: if real {
                        dft_engine::GammaReal::On
                    } else {
                        dft_engine::GammaReal::Off
                    },
                    standard: standard.then(StandardProfile::default),
                    ..Default::default()
                };
                let tolerance = if standard {
                    StandardProfile::CAP_HA
                } else {
                    1e-2
                };
                let mut reference: Option<(ScfStep, ScfStep, Array3<f64>)> = None;
                let mut output = String::new();
                for seed in 1..=10 {
                    let mut loose =
                        ScfProblem::new(&cell, &atoms, &species, options.clone()).unwrap();
                    let mut tight_options = options.clone();
                    tight_options.standard = None;
                    let mut tight =
                        ScfProblem::new(&cell, &atoms, &species, tight_options).unwrap();
                    let electrons = loose.electron_count();
                    let occupied = (electrons / 2.0).round() as usize;
                    let mut n = loose.uniform_density();
                    let perturbation = noise(seed, n.len());
                    for (v, d) in n.iter_mut().zip(perturbation) {
                        *v *= 1.0 + d;
                    }
                    let a = loose
                        .step_kpoints_with_certification(&n, None, tolerance, false)
                        .unwrap();
                    let b = tight
                        .step_kpoints_with_certification(&n, None, 1e-10, true)
                        .unwrap();
                    let (ra, rb, rn) = reference
                        .as_ref()
                        .map_or((&a, &b, &n), |(a, b, n)| (a, b, n));
                    let sa = &a.k_points[0].solution;
                    let sb = &b.k_points[0].solution;
                    let distance = |x, y| density_distance(x, y, cell.volume_bohr3(), electrons);
                    writeln!(output,"FIRST,{case},{real},{standard},{seed},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                        distance(&n,rn),distance(&a.density,&ra.density),distance(&b.density,&rb.density),
                        occupied_subspace_defect(&ra.k_points[0].solution.eigenvectors,&sa.eigenvectors,occupied),
                        occupied_subspace_defect(&rb.k_points[0].solution.eigenvectors,&sb.eigenvectors,occupied),
                        sa.eigenvalues[occupied]-sa.eigenvalues[occupied-1],
                        sb.eigenvalues[occupied]-sb.eigenvalues[occupied-1],
                        sa.eigenvalues.iter().zip(&ra.k_points[0].solution.eigenvalues).take(occupied+1).map(|(a,b)|(a-b).abs()).fold(0.0,f64::max),
                        sa.residuals.iter().zip(&ra.k_points[0].solution.residuals).take(occupied+1).map(|(a,b)|(a-b).abs()).fold(0.0,f64::max)).unwrap();
                    if reference.is_none() {
                        reference = Some((a, b, n));
                    }
                }
                print!("{output}");
            }
        }
    }
}
