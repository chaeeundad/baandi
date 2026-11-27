//! GK1 / GK1B / GK14: the k atomic start (memoized radial transforms,
//! rows-only projector tables, shared collinear spin start, radial tables) is
//! the only path since CL1 (2026-11-27). The former A/B tests against the
//! retired arms (`K_START_MEMO`, `K_NL_TABLE_MEMO`, `K_START_SPIN_SHARE`,
//! `RADIAL_TABLE_START`, `RADIAL_TABLE_NL` off: bit-identical for the memo,
//! rows-only tables and shared spin start, strict class for the tables) are
//! recorded in docs/performance/attempts-ledger.md and
//! docs/performance/switch-retirement.md. These tests keep the start shape and
//! its bit reproducibility at 1 and 3 simulated ranks.
//! Mac observation (release, not an assertion):
//! cargo test --release -p dft-engine --test k_start_memo timing -- --ignored --nocapture
#[allow(dead_code)]
mod common;
use dft_engine::{
    ScfOptions, SpeciesPotential,
    refstart::{self, OrbitalStartKind, StartInput},
};
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use dft_numerics::switches::RunConfig;

struct Fixture {
    a: dft_upf::Pseudopotential,
    b: dft_upf::Pseudopotential,
    cell: Cell,
    atoms: Vec<Atom>,
}

fn fixture() -> Fixture {
    Fixture {
        a: common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ"),
        b: common::potential_with_atomic_start(2., 1.4, false, 201, "SLA PZ"),
        // Skewed cell so that |k+G| shells are irregular.
        cell: Cell::from_user_rows_bohr([[7.2, 0., 0.], [1.1, 6.8, 0.], [0.4, -0.7, 7.5]]).unwrap(),
        atoms: vec![
            Atom::new("A", [0.01, 0.02, 0.03]),
            Atom::new("A", [0.51, 0.47, 0.12]),
            Atom::new("A", [0.26, 0.74, 0.61]),
            Atom::new("B", [0.77, 0.21, 0.38]),
            Atom::new("B", [0.33, 0.29, 0.86]),
        ],
    }
}

/// Gamma (degenerate |G| shells), a symmetric line point and a generic point.
fn points() -> Vec<KPoint> {
    vec![
        KPoint {
            fractional: [0.; 3],
            weight: 0.25,
        },
        KPoint {
            fractional: [0.25, 0., 0.],
            weight: 0.25,
        },
        KPoint {
            fractional: [0.13, -0.29, 0.21],
            weight: 0.5,
        },
    ]
}

fn bits(x: &dft_numerics::linalg::ComplexMatrix) -> Vec<u64> {
    x.as_slice()
        .iter()
        .flat_map(|z| [z.re.to_bits(), z.im.to_bits()])
        .collect()
}

type StartBits = (Vec<Vec<u64>>, OrbitalStartKind, f64, Vec<Vec<u64>>);

/// Start orbitals (bits per k block), the start kind, the start wall and the
/// bits of H X of every k operator applied to its start orbitals.
fn run(
    f: &Fixture,
    points: &[KPoint],
    options: &ScfOptions,
    switches: &[(&str, &str)],
    ranks: usize,
) -> Vec<StartBits> {
    let potentials = [
        SpeciesPotential {
            label: "A",
            potential: &f.a,
        },
        SpeciesPotential {
            label: "B",
            potential: &f.b,
        },
    ];
    dft_mpi::simulate(ranks, |comm| {
        // Pools off: the operators below are applied on the world
        // communicator (pooled operators are one-rank; the k-pools default
        // became `auto` after these tests were written).
        let mut pairs = vec![("refscf_k_pools", "off")];
        pairs.extend_from_slice(switches);
        let config = RunConfig::from_pairs(&pairs).unwrap();
        let _scope = config.enter();
        let mut state = refstart::start_kpoints(
            StartInput {
                cell: &f.cell,
                atoms: &f.atoms,
                potentials: &potentials,
                options,
            },
            points,
            comm,
        )
        .unwrap();
        let orbitals: Vec<_> = state.orbitals.0.iter().map(bits).collect();
        let mut hx = Vec::new();
        for (block, x) in state.blocks.iter_mut().zip(&state.orbitals.0) {
            let mut y = dft_numerics::linalg::ComplexMatrix::zeros(x.rows(), x.cols()).unwrap();
            block.operator.apply_block(comm, x, &mut y).unwrap();
            hx.push(bits(&y));
        }
        let seconds = state.setup.phases["complex_k_atomic_starts"];
        (orbitals, state.report.orbitals, seconds, hx)
    })
}

#[test]
fn gk1_k_start_is_atomic_and_bit_reproducible() {
    let f = fixture();
    let options = ScfOptions {
        cutoff_hartree: 2.5,
        bands: Some(8),
        ..Default::default()
    };
    for ranks in [1, 3] {
        let first = run(&f, &points(), &options, &[], ranks);
        let again = run(&f, &points(), &options, &[], ranks);
        for (rank, (a, b)) in first.iter().zip(&again).enumerate() {
            // 5 atoms x (s + p) = 20 atomic columns cover the 8 bands.
            assert_eq!(a.1, OrbitalStartKind::Atomic);
            assert_eq!(a.0.len(), 3);
            assert!(
                a.0.iter()
                    .all(|k| !k.is_empty() && k.iter().any(|&x| x != 0))
            );
            assert_eq!(
                a.0, b.0,
                "ranks {ranks} rank {rank}: start orbital bits differ"
            );
            assert_eq!(a.3.len(), 3);
            assert_eq!(a.3, b.3, "ranks {ranks} rank {rank}: H X bits differ");
        }
    }
}

/// A collinear fixed-spin k run (3 maps, LDA PW, 2 species, 3 k) from the
/// shared spin start is bit-reproducible: history energies, channel densities,
/// eigenvalues and final orbitals.
#[test]
fn gk1b_spin_shared_start_is_bit_reproducible() {
    use dft_engine::{
        collinear::SpinSector,
        refscf::{
            electronic::{self, Policy},
            profile::Profile,
        },
    };
    let f = fixture();
    let potentials = [
        SpeciesPotential {
            label: "A",
            potential: &f.a,
        },
        SpeciesPotential {
            label: "B",
            potential: &f.b,
        },
    ];
    let options = ScfOptions {
        cutoff_hartree: 2.5,
        bands: Some(8),
        max_iterations: 3,
        functional: dft_engine::Functional::LdaPw,
        allow_functional_override: true,
        gamma_real: dft_engine::GammaReal::Off,
        verification: dft_engine::Verification::FixedPoint,
        standard: Some(dft_engine::StandardProfile { threshold_ha: 5e-9 }),
        ..ScfOptions::strict_reference()
    };
    let sector = SpinSector {
        n_up: 6,
        n_down: 4,
        axis_cartesian: [0., 0., 1.],
        initial_moment_seed: None,
    };
    for ranks in [1, 2] {
        let arm = || {
            dft_mpi::simulate(ranks, |comm| {
                let r = electronic::run_with_trace(
                    StartInput {
                        cell: &f.cell,
                        atoms: &f.atoms,
                        potentials: &potentials,
                        options: &options,
                    },
                    Profile::default(),
                    &points(),
                    Policy::FixedSpin(sector.clone()),
                    comm,
                    None,
                )
                .unwrap();
                let mut out: Vec<u64> = r.history.iter().map(|i| i.energy_ha.to_bits()).collect();
                for d in &r.density {
                    out.extend(d.iter().map(|v| v.to_bits()));
                }
                for channel in &r.blocks {
                    for b in channel {
                        out.extend(b.eigenvalues_ha.iter().map(|v| v.to_bits()));
                        out.extend(bits(&b.orbitals));
                    }
                }
                (r.history.len(), out)
            })
        };
        let on = arm();
        let off = arm();
        assert!(on[0].0 > 0);
        assert_eq!(on, off, "ranks {ranks}: spin run bits differ");
    }
}

/// Observation only: start phase split on a larger
/// case (10 atoms x (s + p) = 40 atomic columns, 1201-point meshes, a few k).
/// Environment: GK1_ECUT (Ha, default 15), GK1_K (k count 1..8, default 4).
#[test]
#[ignore = "Mac timing observation; run in release"]
fn timing() {
    let ecut: f64 = std::env::var("GK1_ECUT").map_or(15., |v| v.parse().unwrap());
    let nk: usize = std::env::var("GK1_K").map_or(4, |v| v.parse().unwrap());
    let f = Fixture {
        a: common::potential_with_atomic_start(2., 1.8, true, 1201, "SLA PZ"),
        b: common::potential_with_atomic_start(2., 1.4, true, 1201, "SLA PZ"),
        cell: Cell::from_user_rows_bohr([[12.2, 0., 0.], [1.1, 11.8, 0.], [0.4, -0.7, 12.5]])
            .unwrap(),
        atoms: (0..10)
            .map(|i| {
                let t = i as f64;
                Atom::new(
                    if i < 5 { "A" } else { "B" },
                    [
                        (0.13 * t + 0.01) % 1.,
                        (0.37 * t + 0.05) % 1.,
                        (0.71 * t + 0.11) % 1.,
                    ],
                )
            })
            .collect(),
    };
    let options = ScfOptions {
        cutoff_hartree: ecut,
        bands: Some(24),
        ..Default::default()
    };
    let generic = [0.13, -0.29, 0.21];
    let points: Vec<_> = (0..nk)
        .map(|i| KPoint {
            fractional: generic.map(|x| x * (i as f64 + 1.) / nk as f64),
            weight: 1. / nk as f64,
        })
        .collect();
    let potentials = [
        SpeciesPotential {
            label: "A",
            potential: &f.a,
        },
        SpeciesPotential {
            label: "B",
            potential: &f.b,
        },
    ];
    let arms: [(&str, &[(&str, &str)]); 1] = [("default", &[])];
    let ranks: usize = std::env::var("GK1_RANKS").map_or(4, |v| v.parse().unwrap());
    for (memo, pairs) in arms {
        let report = dft_mpi::simulate(ranks, |comm| {
            let config = RunConfig::from_pairs(pairs).unwrap();
            let _scope = config.enter();
            let clock = std::time::Instant::now();
            let state = refstart::start_kpoints(
                StartInput {
                    cell: &f.cell,
                    atoms: &f.atoms,
                    potentials: &potentials,
                    options: &options,
                },
                &points,
                comm,
            )
            .unwrap();
            let wall = clock.elapsed().as_secs_f64();
            (
                wall,
                state.blocks[0].basis.vectors().len(),
                state.grid.dims(),
                state.setup.phases.clone(),
            )
        });
        let (wall, npw, dims, phases) = &report[0];
        eprintln!(
            "GK1_TIMING arm={memo} ranks={ranks} k={nk} npw0={npw} grid={dims:?} rank0_wall_s={wall:.3}"
        );
        for (name, s) in phases {
            eprintln!("GK1_TIMING arm={memo}   {name} {s:.4}");
        }
    }
}
