//! GK1: the memoized k atomic start (`BAANDI_K_START_MEMO`, default on) gives
//! start orbitals bit-identical to the original per-(row, atom, m) quadrature.
//! The memo lives in `finish_kpoints`, shared by the unpolarized, spin/FD and
//! fixed-density k starts.
//! Mac observation (release, not an assertion):
//! cargo test --release -p dft-engine --test k_start_memo timing -- --ignored --nocapture
#[allow(dead_code)]
mod common;
use dft_engine::{
    ScfOptions, SpeciesPotential,
    refstart::{self, OrbitalStartKind, StartInput},
};
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use dft_numerics::switches::{Id, RunConfig};

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

/// Start orbitals (bits per k block, and the start kind) with the memo switch.
fn run(
    f: &Fixture,
    points: &[KPoint],
    options: &ScfOptions,
    memo: &str,
    ranks: usize,
) -> Vec<(Vec<Vec<u64>>, OrbitalStartKind, f64)> {
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
        let config = RunConfig::from_pairs(&[("k_start_memo", memo)]).unwrap();
        let _scope = config.enter();
        let state = refstart::start_kpoints(
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
        assert!(config.consulted(Id::K_START_MEMO));
        let bits = state
            .orbitals
            .0
            .iter()
            .map(|x| {
                x.as_slice()
                    .iter()
                    .flat_map(|z| [z.re.to_bits(), z.im.to_bits()])
                    .collect()
            })
            .collect();
        let seconds = state.setup.phases["complex_k_atomic_starts"];
        (bits, state.report.orbitals, seconds)
    })
}

#[test]
fn gk1_k_start_memo_is_bit_identical() {
    let f = fixture();
    let options = ScfOptions {
        cutoff_hartree: 2.5,
        bands: Some(8),
        ..Default::default()
    };
    for ranks in [1, 3] {
        let on = run(&f, &points(), &options, "on", ranks);
        let off = run(&f, &points(), &options, "off", ranks);
        for (rank, (a, b)) in on.iter().zip(&off).enumerate() {
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
        }
    }
}

/// Observation only: start phase split with memo on and off on a larger
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
    for memo in ["on", "off", "on"] {
        let config = RunConfig::from_pairs(&[("k_start_memo", memo)]).unwrap();
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
            &dft_engine::refscf::metric::Serial,
        )
        .unwrap();
        let wall = clock.elapsed().as_secs_f64();
        eprintln!(
            "GK1_TIMING memo={memo} k={nk} npw0={} grid={:?} wall_s={wall:.3}",
            state.blocks[0].basis.vectors().len(),
            state.grid.dims()
        );
        for (name, s) in &state.setup.phases {
            eprintln!("GK1_TIMING memo={memo}   {name} {s:.4}");
        }
    }
}
