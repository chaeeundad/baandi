//! P42A wide fixed-density audit; residual certification is not completeness certification.
//! Run in release mode, redirecting stdout to durable evidence storage:
//! P42A_SI_UPF=/absolute/Si.pz-vbc.UPF P42A_DENSITY=/path/density-1.bin
//! P42A_CASES=si-gamma cargo test --release -p dft-engine
//! --test band_boundary_spectrum -- --ignored --nocapture --test-threads=1
use dft_engine::{Functional, Mixing, OccupationPolicy, ScfOptions, ScfProblem, SpeciesPotential};
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

#[test]
#[ignore = "multi-case SCF study; release only, external UPF and durable stdout required"]
fn frozen_wide_spectrum() {
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
        {
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
            opts.eigen_tolerance = 1e-5;
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
            let input =
                PathBuf::from(std::env::var("P42A_DENSITY").expect("P42A_DENSITY required"));
            let mut chosen = None;
            // The immutable basis can be smaller than occupied+17 in toy fixtures.
            for extra in (1..=17).rev() {
                opts.extra_bands = extra;
                if let Ok(problem) = make(opts.clone()) {
                    chosen = Some((extra, problem));
                    break;
                }
            }
            let (extra, mut problem) = chosen.expect("no supported probe band count");
            let mut density = problem.uniform_density();
            let bytes = std::fs::read(input).unwrap();
            assert_eq!(bytes.len(), density.len() * 8);
            for (value, raw) in density.iter_mut().zip(bytes.as_chunks::<8>().0) {
                *value = f64::from_le_bytes(*raw);
                assert!(value.is_finite());
            }
            match problem.step_kpoints(&density, None, 1e-5) {
                Ok(probe) => {
                    for (k, point) in probe.k_points.iter().enumerate() {
                        println!(
                            "WIDE|{case}|{k}|{extra}|{:?}|{:?}|{}|{}",
                            point.solution.eigenvalues,
                            point.solution.residuals,
                            point.solution.residuals_certified,
                            point.solution.operator_applications
                        );
                    }
                }
                Err(e) => println!("WIDE_ERROR|{case}|{e}"),
            }
        }
    }
}
