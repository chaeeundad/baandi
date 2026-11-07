//! Fixed-density probes for the five ambiguous P38A spectrum classifications.
//! Release only; see docs/performance/band-count-policy.md for input provenance.
//! P38A_SI_UPF, P38A_PD_SI_UPF, P38A_SYN_S2_UPF and P38A_SPECTRUM_OUT are required.
//! cargo test --release -p dft-engine --test band_count_spectrum -- --ignored --nocapture
use dft_engine::{
    Functional, Mixing, OccupationPolicy, ScfOptions, ScfProblem, SpeciesPotential, StandardProfile,
};
use dft_geometry::{Atom, Cell, kpoints};
use std::path::PathBuf;

// Restricted reader for the repository's numeric fixture fields, not a TOML adapter.
fn numbers(text: &str, key: &str) -> Option<Vec<f64>> {
    let value = text.lines().find_map(|line| {
        let (name, value) = line.split_once('=')?;
        (name.trim() == key).then_some(value)
    })?;
    Some(
        value
            .split(['[', ']', ',', ' ', '\t'])
            .filter(|s| !s.is_empty())
            .map(|s| s.parse().expect("numeric fixture field"))
            .collect(),
    )
}

#[test]
#[ignore = "external UPFs and repeated SCF/frozen-spectrum solves; release measurement only"]
fn frozen_spectra_for_ambiguous_rows() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let out =
        PathBuf::from(std::env::var("P38A_SPECTRUM_OUT").expect("P38A_SPECTRUM_OUT required"));
    std::fs::create_dir_all(&out).unwrap();
    for (case, extra, standard, material, upf_env) in [
        ("si16-gamma", 3, true, false, "P38A_SI_UPF"),
        ("si16-gamma", 4, true, false, "P38A_SI_UPF"),
        ("pd-si", 2, true, true, "P38A_PD_SI_UPF"),
        ("syn-s2", 1, false, false, "P38A_SYN_S2_UPF"),
        ("syn-s2", 2, false, false, "P38A_SYN_S2_UPF"),
    ] {
        let path = if material {
            root.join("tests/reference/materials/si-20.toml")
        } else {
            root.join(format!("benchmarks/cases/{case}.toml"))
        };
        let text = std::fs::read_to_string(path).unwrap();
        let lattice = numbers(
            &text,
            if material {
                "lattice_vectors"
            } else {
                "lattice"
            },
        )
        .unwrap();
        let cell = Cell::from_user_rows_bohr(std::array::from_fn(|i| {
            std::array::from_fn(|j| lattice[3 * i + j])
        }))
        .unwrap();
        let xyz = if material {
            text.lines()
                .filter(|line| line.starts_with("fractional ="))
                .flat_map(|line| numbers(line, "fractional").unwrap())
                .collect()
        } else {
            numbers(&text, "atoms").unwrap()
        };
        let atoms: Vec<_> = xyz
            .as_chunks::<3>()
            .0
            .iter()
            .map(|v| Atom::new("Si", *v))
            .collect();
        let pp = dft_upf::parse_upf(
            &std::fs::read_to_string(std::env::var(upf_env).expect("external UPF path required"))
                .unwrap(),
        )
        .unwrap();
        let species = [SpeciesPotential {
            label: "Si",
            potential: &pp,
        }];
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
        let mut options = ScfOptions {
            cutoff_hartree: numbers(&text, if material { "ecutwfc_ry" } else { "cutoff_ry" })
                .unwrap()[0]
                * 0.5,
            extra_bands: extra,
            max_iterations: if material { 120 } else { 100 },
            seed: if material { 2026 } else { 42 },
            mixing: Mixing::Pulay {
                beta: 0.5,
                history: 6,
            },
            standard: standard.then(StandardProfile::default),
            ..Default::default()
        };
        if material {
            options.functional = Functional::Pbe;
            options.mixing = Mixing::KerkerPulay {
                beta: 0.3,
                history: 6,
                q0_inverse_bohr: 1.0,
            };
        }
        let make = |opts| {
            let mut problem = ScfProblem::with_kpoints(
                &cell,
                &atoms,
                &species,
                opts,
                &mesh,
                OccupationPolicy::Insulator {
                    minimum_gap_hartree: 1e-7,
                },
            )
            .unwrap();
            if material {
                problem.set_density_floor(1e-10).unwrap();
            }
            problem
        };
        let original = make(options.clone()).run().unwrap();
        assert_eq!(original.certified, !standard);
        let profile = if standard { "standard" } else { "strict" };
        // Persist the replayed density outside the solver so Python can independently
        // compare its SHA-256 with the original CLI density metadata.
        let bytes: Vec<_> = original
            .density
            .iter()
            .flat_map(|v| v.to_le_bytes())
            .collect();
        std::fs::write(
            out.join(format!("{case}-{profile}-e{extra:02}-density.bin")),
            bytes,
        )
        .unwrap();
        println!(
            "P38A_ORIGIN|{case}|{profile}|{extra}|{:.17e}",
            original.total_hartree()
        );

        // Equal mandatory/buffer targets in this *frozen-density diagnostic*.
        // 1e-5 Ha bounds every residual, below the 1e-4 Ha cluster separation.
        // Neither the original SCF trajectory nor its tolerances are changed.
        const PROBE_TOLERANCE_HA: f64 = 1e-5;
        options.standard = None;
        options.extra_bands = if case == "syn-s2" { 5 } else { 12 };
        options.eigen_tolerance = PROBE_TOLERANCE_HA;
        options.grid_override = Some(std::array::from_fn(|i| original.density.shape()[i]));
        let probe = make(options)
            .step_kpoints(&original.density, None, PROBE_TOLERANCE_HA)
            .unwrap();
        let occupied = original.occupations.iter().filter(|&&v| v > 0.).count();
        for (k, point) in probe.k_points.iter().enumerate() {
            assert!(point.solution.residuals_certified);
            assert!(
                point
                    .solution
                    .residuals
                    .iter()
                    .all(|r| r.is_finite() && *r <= PROBE_TOLERANCE_HA)
            );
            println!(
                "P38A_SPECTRUM|{case}|{profile}|{extra}|{k}|{occupied}|{:?}|{:?}",
                &point.solution.eigenvalues[occupied - 1..],
                &point.solution.residuals[occupied - 1..]
            );
        }
    }
}
