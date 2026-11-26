//! GK14: radial transform tables for the k start (`BAANDI_RADIAL_TABLE_START`)
//! and the k projector factors (`BAANDI_RADIAL_TABLE_NL`), both default on.
//! The tables are not bit-identical to the quadrature; acceptance is the
//! delivered energy within the strict class (1e-10 Ha) of the quadrature arm
//! on the k generality set: cubic, skewed and elongated cells, Gamma in the
//! list, a shifted mesh, collinear fixed spin and Fermi-Dirac. Synthetic data,
//! never QE evidence.
//!
//! Real files (ignored): DFT_NC_DIR=/abs/pseudodojo/ex cargo test --release \
//!   -p dft-engine --test k_radial_table real_gete -- --ignored --nocapture
//! Mac stage timing (ignored): cargo test --release -p dft-engine \
//!   --test k_radial_table timing -- --ignored --nocapture
#[allow(dead_code)]
mod common;
use dft_engine::{
    ScfOptions, SpeciesPotential, Verification,
    collinear::SpinSector,
    refscf::{
        electronic::{self, Policy},
        k_driver,
        profile::Profile,
    },
    refstart::{self, StartInput},
};
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use dft_numerics::switches::{Id, RunConfig};

const STRICT_HA: f64 = 1e-10;

/// The four arms: quadrature reference, start table only, NL table only, both.
const ARMS: [(&str, &[(&str, &str)]); 4] = [
    (
        "quadrature",
        &[("radial_table_start", "off"), ("radial_table_nl", "off")],
    ),
    (
        "start",
        &[("radial_table_start", "on"), ("radial_table_nl", "off")],
    ),
    (
        "nl",
        &[("radial_table_start", "off"), ("radial_table_nl", "on")],
    ),
    ("default", &[]),
];

struct Species2 {
    a: dft_upf::Pseudopotential,
    b: dft_upf::Pseudopotential,
}
fn species() -> Species2 {
    Species2 {
        a: common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ"),
        b: common::potential_with_atomic_start(2., 1.4, true, 201, "SLA PZ"),
    }
}
fn potentials(s: &Species2) -> [SpeciesPotential<'_>; 2] {
    [
        SpeciesPotential {
            label: "A",
            potential: &s.a,
        },
        SpeciesPotential {
            label: "B",
            potential: &s.b,
        },
    ]
}
fn atoms() -> Vec<Atom> {
    vec![
        Atom::new("A", [0.01, 0.02, 0.03]),
        Atom::new("B", [0.52, 0.47, 0.41]),
    ]
}
fn kp(f: [f64; 3], w: f64) -> KPoint {
    KPoint {
        fractional: f,
        weight: w,
    }
}

/// (label, cell rows, k list)
type Case = (&'static str, [[f64; 3]; 3], Vec<KPoint>);
fn cases() -> Vec<Case> {
    let gamma_generic = vec![kp([0.; 3], 0.25), kp([0.13, -0.29, 0.21], 0.75)];
    let shifted = vec![
        kp([0.25, 0.25, 0.25], 0.25),
        kp([0.25, 0.25, -0.25], 0.25),
        kp([0.25, -0.25, 0.25], 0.25),
        kp([-0.25, 0.25, 0.25], 0.25),
    ];
    vec![
        (
            "cubic",
            [[6.5, 0., 0.], [0., 6.5, 0.], [0., 0., 6.5]],
            gamma_generic.clone(),
        ),
        (
            "skewed",
            [[6.6, 0., 0.], [1.1, 6.2, 0.], [0.4, -0.7, 6.9]],
            gamma_generic.clone(),
        ),
        (
            "elongated",
            [[5., 0., 0.], [0., 5., 0.], [0., 0., 11.]],
            gamma_generic,
        ),
        (
            "shifted_mesh",
            [[6.6, 0., 0.], [1.1, 6.2, 0.], [0.4, -0.7, 6.9]],
            shifted,
        ),
    ]
}

fn options() -> ScfOptions {
    ScfOptions {
        cutoff_hartree: 3.,
        bands: Some(6),
        max_iterations: 100,
        verification: Verification::FixedPoint,
        ..Default::default()
    }
}

/// Unpolarized k runs: every arm within the strict class of the quadrature
/// arm, converged; map counts printed.
#[test]
fn unpolarized_k_set_within_strict_class() {
    let s = species();
    let potentials = potentials(&s);
    let atoms = atoms();
    let options = options();
    for (label, rows, points) in cases() {
        let cell = Cell::from_user_rows_bohr(rows).unwrap();
        for ranks in [1, 2] {
            let out = dft_mpi::simulate(ranks, |comm| {
                ARMS.map(|(arm, pairs)| {
                    let config = RunConfig::from_pairs(pairs).unwrap();
                    let _scope = config.enter();
                    let run = k_driver::run_with(
                        StartInput {
                            cell: &cell,
                            atoms: &atoms,
                            potentials: &potentials,
                            options: &options,
                        },
                        Profile::default(),
                        &points,
                        comm,
                    )
                    .unwrap();
                    assert!(config.consulted(Id::RADIAL_TABLE_START), "{arm}");
                    assert!(config.consulted(Id::RADIAL_TABLE_NL), "{arm}");
                    assert!(run.outcome.converged, "{label} {arm}: {}", run.exit_reason);
                    (run.outcome.result.total_hartree(), run.records.len())
                })
            });
            let arms = &out[0];
            let (e0, m0) = arms[0];
            for ((arm, _), (e, maps)) in ARMS.iter().zip(arms) {
                let d = e - e0;
                println!(
                    "GK14 {label} P={ranks} arm={arm} E={e:.15} Ha dE={d:.3e} maps={maps} (quadrature {m0})"
                );
                assert!(d.abs() <= STRICT_HA, "{label} P={ranks} {arm}: {d:e}");
            }
        }
    }
}

/// Collinear fixed spin and Fermi-Dirac (spin-unpolarized) on the skewed cell.
#[test]
fn spin_and_fd_k_within_strict_class() {
    let s = species();
    let potentials = potentials(&s);
    let atoms = atoms();
    let cell =
        Cell::from_user_rows_bohr([[6.6, 0., 0.], [1.1, 6.2, 0.], [0.4, -0.7, 6.9]]).unwrap();
    let points = [kp([0.; 3], 0.25), kp([0.13, -0.29, 0.21], 0.75)];
    let options = ScfOptions {
        functional: dft_engine::Functional::LdaPw,
        allow_functional_override: true,
        gamma_real: dft_engine::GammaReal::Off,
        ..options()
    };
    let policies = [
        (
            "fixed_spin",
            Policy::FixedSpin(SpinSector {
                n_up: 3,
                n_down: 1,
                axis_cartesian: [0., 0., 1.],
                initial_moment_seed: None,
            }),
        ),
        (
            "fermi_dirac",
            Policy::FermiDirac {
                theta_ha: 0.01,
                band_cap: 12,
            },
        ),
    ];
    for (label, policy) in policies {
        let out = dft_mpi::simulate(2, |comm| {
            ARMS.map(|(_, pairs)| {
                let config = RunConfig::from_pairs(pairs).unwrap();
                let _scope = config.enter();
                let r = electronic::run_with_trace(
                    StartInput {
                        cell: &cell,
                        atoms: &atoms,
                        potentials: &potentials,
                        options: &options,
                    },
                    Profile::default(),
                    &points,
                    policy.clone(),
                    comm,
                    None,
                )
                .unwrap();
                assert!(r.converged, "{label}");
                (r.energy.total_hartree(), r.history.len())
            })
        });
        let (e0, _) = out[0][0];
        for ((arm, _), (e, maps)) in ARMS.iter().zip(&out[0]) {
            let d = e - e0;
            println!("GK14 {label} arm={arm} E={e:.15} Ha dE={d:.3e} maps={maps}");
            assert!(d.abs() <= STRICT_HA, "{label} {arm}: {d:e}");
        }
    }
}

/// Real PseudoDojo Ge and Te (PBE) on a skewed 2-atom GeTe-like cell with
/// generic k: delivered energy of every arm against the quadrature arm.
#[test]
#[ignore = "requires DFT_NC_DIR (Ge.upf, Te.upf); release"]
fn real_gete() {
    let dir = std::path::PathBuf::from(std::env::var_os("DFT_NC_DIR").expect("set DFT_NC_DIR"));
    let ge = dft_upf::from_file(dir.join("Ge.upf")).unwrap();
    let te = dft_upf::from_file(dir.join("Te.upf")).unwrap();
    let potentials = [
        SpeciesPotential {
            label: "Ge",
            potential: &ge,
        },
        SpeciesPotential {
            label: "Te",
            potential: &te,
        },
    ];
    let atoms = [
        Atom::new("Ge", [0.02, 0.01, 0.03]),
        Atom::new("Te", [0.48, 0.51, 0.47]),
    ];
    let cell = Cell::from_user_rows_bohr([[8.0, 0., 0.], [3.9, 6.9, 0.], [3.9, 2.3, 6.5]]).unwrap();
    let points = [
        kp([0.; 3], 0.125),
        kp([0.13, -0.29, 0.21], 0.375),
        kp([0.37, 0.11, -0.23], 0.5),
    ];
    let cutoff: f64 = std::env::var("GK14_ECUT").map_or(15., |v| v.parse().unwrap());
    let options = ScfOptions {
        functional: dft_engine::Functional::Pbe,
        cutoff_hartree: cutoff,
        bands: Some((ge.z_valence() + te.z_valence()) as usize / 2 + 4),
        max_iterations: 100,
        verification: Verification::FixedPoint,
        ..Default::default()
    };
    let out = dft_mpi::simulate(2, |comm| {
        ARMS.map(|(_, pairs)| {
            let config = RunConfig::from_pairs(pairs).unwrap();
            let _scope = config.enter();
            let run = k_driver::run_with(
                StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &options,
                },
                Profile::default(),
                &points,
                comm,
            )
            .unwrap();
            assert!(run.outcome.converged, "{}", run.exit_reason);
            let eig: Vec<f64> = run
                .outcome
                .result
                .k_points
                .iter()
                .flat_map(|k| k.solution.eigenvalues.clone())
                .collect();
            (run.outcome.result.total_hartree(), run.records.len(), eig)
        })
    });
    let (e0, _, eig0) = &out[0][0];
    for ((arm, _), (e, maps, eig)) in ARMS.iter().zip(&out[0]) {
        let d = e - e0;
        let de = eig
            .iter()
            .zip(eig0)
            .map(|(a, b)| (a - b).abs())
            .fold(0., f64::max);
        println!(
            "GK14 real_gete ecut={cutoff} arm={arm} E={e:.15} Ha dE={d:.3e} max_d_eig={de:.3e} Ha maps={maps}"
        );
        assert!(d.abs() <= STRICT_HA, "{arm}: {d:e}");
    }
}

/// Observation only: k start stage split (k_radial_tables,
/// k_start_nonlocal_tables, k_start_atomic_columns) with the tables off and on
/// on a GeTe-like synthetic: 1501-point meshes, 3 projectors and 3 orbitals
/// per species, skewed 4-atom cell, ~5000 plane waves per k at 28 Ha, generic
/// k. Environment: GK14_K (default 8), GK14_ECUT (default 28),
/// GK14_RANKS (default 1). With DFT_NC_DIR set, real Ge/Te files are used.
#[test]
#[ignore = "Mac timing observation; run in release"]
fn timing() {
    let nk: usize = std::env::var("GK14_K").map_or(8, |v| v.parse().unwrap());
    let ecut: f64 = std::env::var("GK14_ECUT").map_or(28., |v| v.parse().unwrap());
    let ranks: usize = std::env::var("GK14_RANKS").map_or(1, |v| v.parse().unwrap());
    let (a, b) = match std::env::var_os("DFT_NC_DIR") {
        Some(dir) => {
            let dir = std::path::PathBuf::from(dir);
            (
                dft_upf::from_file(dir.join("Ge.upf")).unwrap(),
                dft_upf::from_file(dir.join("Te.upf")).unwrap(),
            )
        }
        None => (
            common::potential_with_atomic_start(4., 1.8, true, 1501, "SLA PZ"),
            common::potential_with_atomic_start(6., 1.6, true, 1501, "SLA PZ"),
        ),
    };
    let real = std::env::var_os("DFT_NC_DIR").is_some();
    let potentials = [
        SpeciesPotential {
            label: "A",
            potential: &a,
        },
        SpeciesPotential {
            label: "B",
            potential: &b,
        },
    ];
    let atoms = vec![
        Atom::new("A", [0.01, 0.02, 0.03]),
        Atom::new("A", [0.51, 0.47, 0.52]),
        Atom::new("B", [0.26, 0.24, 0.27]),
        Atom::new("B", [0.77, 0.74, 0.76]),
    ];
    let cell =
        Cell::from_user_rows_bohr([[9.6, 0., 0.], [1.3, 9.2, 0.], [0.9, -0.8, 9.9]]).unwrap();
    let options = ScfOptions {
        cutoff_hartree: ecut,
        bands: Some((2. * (a.z_valence() + b.z_valence())) as usize / 2 + 6),
        functional: if real {
            dft_engine::Functional::Pbe
        } else {
            dft_engine::Functional::LdaPz
        },
        ..Default::default()
    };
    let generic = [0.13, -0.29, 0.21];
    let points: Vec<_> = (0..nk)
        .map(|i| {
            kp(
                generic.map(|x| x * (i as f64 + 1.) / nk as f64),
                1. / nk as f64,
            )
        })
        .collect();
    let arms: [(&str, &[(&str, &str)]); 3] = [
        (
            "off",
            &[("radial_table_start", "off"), ("radial_table_nl", "off")],
        ),
        ("on", &[]),
        (
            "off",
            &[("radial_table_start", "off"), ("radial_table_nl", "off")],
        ),
    ];
    for (arm, pairs) in arms {
        let report = dft_mpi::simulate(ranks, |comm| {
            let config = RunConfig::from_pairs(pairs).unwrap();
            let _scope = config.enter();
            let clock = std::time::Instant::now();
            let state = refstart::start_kpoints(
                StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &options,
                },
                &points,
                comm,
            )
            .unwrap();
            (
                clock.elapsed().as_secs_f64(),
                state.blocks[0].basis.vectors().len(),
                state.setup.phases.clone(),
            )
        });
        let (wall, npw, phases) = &report[0];
        eprintln!(
            "GK14_TIMING arm={arm} real={real} ranks={ranks} k={nk} npw0={npw} wall={wall:.3}"
        );
        for name in [
            "k_radial_tables",
            "k_start_nonlocal_tables",
            "k_start_atomic_columns",
            "complex_k_atomic_starts",
        ] {
            eprintln!(
                "GK14_TIMING arm={arm}   {name} {:.4}",
                phases.get(name).copied().unwrap_or(0.)
            );
        }
    }
}
