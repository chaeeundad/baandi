//! GK1: the memoized k atomic start (`BAANDI_K_START_MEMO`, default on) gives
//! start orbitals bit-identical to the original per-(row, atom, m) quadrature.
//! The memo lives in `finish_kpoints`, shared by the unpolarized, spin/FD and
//! fixed-density k starts.
//! GK1B: rows-only projector tables (`BAANDI_K_NL_TABLE_MEMO`) keep the start
//! orbitals and the operators' H X bit-identical; the shared collinear spin
//! start (`BAANDI_K_START_SPIN_SHARE`) keeps a whole spin run bit-identical.
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
        // GK14 radial tables off: these tests compare against the quadrature bits.
        // Pools off: the operators below are applied on the world
        // communicator (pooled operators are one-rank; the k-pools default
        // became `auto` after these tests were written).
        let mut pairs = vec![
            ("radial_table_start", "off"),
            ("radial_table_nl", "off"),
            ("refscf_k_pools", "off"),
        ];
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
        assert!(config.consulted(Id::K_START_MEMO));
        assert!(config.consulted(Id::K_NL_TABLE_MEMO));
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
fn gk1_k_start_memo_is_bit_identical() {
    let f = fixture();
    let options = ScfOptions {
        cutoff_hartree: 2.5,
        bands: Some(8),
        ..Default::default()
    };
    for ranks in [1, 3] {
        let on = run(&f, &points(), &options, &[("k_start_memo", "on")], ranks);
        let off = run(&f, &points(), &options, &[("k_start_memo", "off")], ranks);
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

#[test]
fn gk1b_rows_only_projector_tables_are_bit_identical() {
    let f = fixture();
    let options = ScfOptions {
        cutoff_hartree: 2.5,
        bands: Some(8),
        ..Default::default()
    };
    for ranks in [1, 3] {
        let on = run(&f, &points(), &options, &[("k_nl_table_memo", "on")], ranks);
        let off = run(
            &f,
            &points(),
            &options,
            &[("k_nl_table_memo", "off")],
            ranks,
        );
        for (rank, (a, b)) in on.iter().zip(&off).enumerate() {
            assert_eq!(
                a.0, b.0,
                "ranks {ranks} rank {rank}: start orbital bits differ"
            );
            assert_eq!(a.3.len(), 3);
            assert_eq!(a.3, b.3, "ranks {ranks} rank {rank}: H X bits differ");
        }
    }
}

/// A collinear fixed-spin k run (3 maps, LDA PW, 2 species, 3 k) with the shared
/// spin start equals the per-channel start bit for bit: history energies,
/// channel densities, eigenvalues and final orbitals.
#[test]
fn gk1b_spin_shared_start_is_bit_identical() {
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
        let arm = |share: &str| {
            dft_mpi::simulate(ranks, |comm| {
                let config = RunConfig::from_pairs(&[("k_start_spin_share", share)]).unwrap();
                let _scope = config.enter();
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
                assert!(config.consulted(Id::K_START_SPIN_SHARE));
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
        let on = arm("on");
        let off = arm("off");
        assert!(on[0].0 > 0);
        assert_eq!(on, off, "ranks {ranks}: spin run bits differ");
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
    // GK1B arms: "before" = GK1 state (full-sphere projector tables).
    let arms: [(&str, &[(&str, &str)]); 3] = [
        ("before", &[("k_nl_table_memo", "off")]),
        ("after", &[("k_nl_table_memo", "on")]),
        ("before", &[("k_nl_table_memo", "off")]),
    ];
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
