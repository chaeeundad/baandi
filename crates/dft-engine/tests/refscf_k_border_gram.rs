//! GK10: border-only projected S/H in the complex reference Davidson
//! (`BAANDI_CDAV_BORDER_GRAM`, default on) is a rounding-class change. Whole k
//! SCF runs with the border on and off agree within the strict accuracy class
//! (1e-10 Ha) and the standard class (5e-9 Ha) on cubic, skewed and elongated
//! cells with a Gamma-in-list set and a shifted mesh, on Fermi-Dirac and
//! collinear-spin k runs, at 1 and 3 simulated ranks (8 in the ignored
//! matrix). Map counts and H columns are printed. Synthetic data, never QE
//! evidence.
//! Full matrix: cargo test --release -p dft-engine --test refscf_k_border_gram -- --include-ignored --nocapture
#[allow(dead_code)]
mod common;
use dft_engine::{
    Functional, GammaReal, ScfOptions, SpeciesPotential, StandardProfile, Verification,
    collinear::SpinSector,
    refscf::{
        electronic::{self, Policy},
        k_driver,
        profile::Profile,
    },
    refstart::StartInput,
};
use dft_geometry::{
    Atom, Cell,
    kpoints::{KPoint, monkhorst_pack},
};
use dft_mpi::Collective;
use dft_numerics::switches::{Id, RunConfig};

const STRICT_HA: f64 = 1e-10;

/// "off": the full-recompute solver and the pre-GK10B first map (GK10
/// border Gram, GK10B bordered factor, fused corrections and the kept start
/// and redo HX all off); "on": the defaults.
fn arm(border: &str) -> Vec<(&'static str, &str)> {
    if border == "off" {
        vec![
            ("cdav_border_gram", "off"),
            ("cdav_border_factor", "off"),
            ("cdav_stream", "off"),
            ("k_start_kept_hx", "off"),
            ("k_redo_kept_hx", "off"),
        ]
    } else {
        Vec::new()
    }
}
const STANDARD_HA: f64 = 5e-9;

fn cells() -> [(&'static str, Cell); 3] {
    [
        (
            "cubic",
            Cell::from_user_rows_bohr([[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]]).unwrap(),
        ),
        (
            "skewed",
            Cell::from_user_rows_bohr([[7.2, 0., 0.], [1.1, 6.8, 0.], [0.4, -0.7, 7.5]]).unwrap(),
        ),
        (
            "elongated",
            Cell::from_user_rows_bohr([[5., 0., 0.], [0., 5., 0.], [0., 0., 11.]]).unwrap(),
        ),
    ]
}

/// Gamma in an explicit list with a symmetric and a generic point.
fn gamma_in_list() -> Vec<KPoint> {
    vec![
        KPoint {
            fractional: [0.; 3],
            weight: 0.2,
        },
        KPoint {
            fractional: [0.13, -0.21, 0.07],
            weight: 0.3,
        },
        KPoint {
            fractional: [0.25; 3],
            weight: 0.5,
        },
    ]
}

/// The reference k path always runs the standard stopping rule; its
/// threshold comes from the profile. Strict arm: 1e-12 Ry (5e-13 Ha) so the
/// delivered energies resolve the 1e-10 Ha class; standard arm: the default
/// 1e-8 Ry (5e-9 Ha).
fn options() -> ScfOptions {
    ScfOptions {
        cutoff_hartree: 3.,
        bands: Some(6),
        max_iterations: 100,
        verification: Verification::Off,
        gamma_real: GammaReal::Off,
        standard: Some(StandardProfile {
            threshold_ha: STANDARD_HA,
        }),
        ..ScfOptions::strict_reference()
    }
}

fn profile(strict: bool) -> Profile {
    Profile {
        threshold_ry: if strict { 1e-12 } else { 1e-8 },
        ..Profile::default()
    }
}

struct Outcome {
    energy: f64,
    converged: bool,
    maps: usize,
    h_columns: usize,
    /// Per (channel, k) block: eigenvalues and the count of occupied bands
    /// (converged to the occupied-band threshold; the others only to the
    /// empty-band threshold, so they are reported, not gated).
    blocks: Vec<(Vec<f64>, usize)>,
}

/// (|energy difference|, max occupied eigenvalue difference, max over all
/// bands), Ha.
fn diff(a: &Outcome, b: &Outcome) -> (f64, f64, f64) {
    assert_eq!(a.blocks.len(), b.blocks.len());
    let (mut occupied, mut all) = (0_f64, 0_f64);
    for ((x, nx), (y, ny)) in a.blocks.iter().zip(&b.blocks) {
        assert_eq!((x.len(), nx), (y.len(), ny));
        for (i, (p, q)) in x.iter().zip(y).enumerate() {
            let d = (p - q).abs();
            all = all.max(d);
            if i < *nx {
                occupied = occupied.max(d);
            }
        }
    }
    ((a.energy - b.energy).abs(), occupied, all)
}

fn scalar(
    comm: &dft_mpi::SimRank,
    cell: &Cell,
    points: &[KPoint],
    strict: bool,
    border: &str,
) -> Outcome {
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    let atoms = [
        Atom::new("Xx", [0.13, 0.21, 0.32]),
        Atom::new("Xx", [0.61, 0.55, 0.71]),
    ];
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = options();
    let config = RunConfig::from_pairs(&arm(border)).unwrap();
    let _scope = config.enter();
    let run = k_driver::run_with(
        StartInput {
            cell,
            atoms: &atoms,
            potentials: &potentials,
            options: &options,
        },
        profile(strict),
        points,
        comm,
    )
    .unwrap();
    assert!(config.consulted(Id::CDAV_BORDER_GRAM));
    let r = &run.outcome.result;
    Outcome {
        energy: r.energy.total_hartree(),
        converged: run.outcome.converged,
        maps: run.records.len(),
        h_columns: run.records.iter().map(|x| x.h_columns).sum(),
        // Two atoms of valence 2: two doubly occupied bands.
        blocks: r
            .k_points
            .iter()
            .map(|k| (k.solution.eigenvalues.clone(), 2))
            .collect(),
    }
}

fn electronic_run(comm: &dft_mpi::SimRank, spin: bool, strict: bool, border: &str) -> Outcome {
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0.2, 5., 0.], [0., 0.1, 5.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    let pp = common::potential_with_atomic_start(
        1.,
        1.,
        true,
        161,
        if spin { "SLA PW" } else { "SLA PZ" },
    );
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        // Large enough that the 12-band FD block grows by borders instead
        // of repairing at every step (n > dim / 2 forces repairs).
        cutoff_hartree: 8.,
        functional: if spin {
            Functional::LdaPw
        } else {
            Functional::LdaPz
        },
        bands: Some(if spin { 4 } else { 12 }),
        ..options()
    };
    let policy = if spin {
        Policy::FixedSpin(SpinSector {
            n_up: 1,
            n_down: 0,
            axis_cartesian: [0., 0., 1.],
            initial_moment_seed: Some(0.3),
        })
    } else {
        Policy::FermiDirac {
            theta_ha: 0.01,
            band_cap: 12,
        }
    };
    let points = monkhorst_pack([2; 3], [0.5; 3]).unwrap();
    let config = RunConfig::from_pairs(&arm(border)).unwrap();
    let _scope = config.enter();
    let run = electronic::run_with(
        StartInput {
            cell: &cell,
            atoms: &atoms,
            potentials: &potentials,
            options: &options,
        },
        profile(strict),
        &points,
        policy,
        comm,
    )
    .unwrap();
    Outcome {
        energy: run.thermodynamics.free_energy_ha,
        converged: run.converged,
        maps: run.history.len(),
        h_columns: run
            .history
            .iter()
            .flat_map(|i| &i.calls)
            .map(|(_, _, c)| c.h_columns)
            .sum(),
        // Lowest band of the up channel (spin, one electron) or of the
        // single FD channel; FD tails and the down channel are reported only.
        blocks: run
            .blocks
            .iter()
            .enumerate()
            .flat_map(|(c, channel)| {
                channel
                    .iter()
                    .map(move |b| (b.eigenvalues_ha.clone(), usize::from(c == 0)))
            })
            .collect(),
    }
}

fn check(label: &str, strict: bool, comm: &dft_mpi::SimRank, off: &Outcome, on: &Outcome) {
    let bound = if strict { STRICT_HA } else { STANDARD_HA };
    let (energy, occupied, all) = diff(off, on);
    assert!(off.converged && on.converged, "{label}");
    if comm.rank() == 0 {
        println!(
            "GK10 {label}: E_on={:.15e} Ha |on-off| energy {energy:.2e} occupied eigenvalues {occupied:.2e} all bands {all:.2e} Ha; maps off {} on {}; H columns off {} on {}",
            on.energy, off.maps, on.maps, off.h_columns, on.h_columns
        );
    }
    assert!(energy <= bound, "{label}: energy {energy:e} Ha > {bound:e}");
    // Occupied eigenvalues are first order in the self-consistency error
    // (the energy is second order): GK10B's kept start HX changes the
    // rounding of the first map, and the delivered eigenvalues then differ
    // by up to the same 7e-10 Ha that rank counts already give on the skewed
    // Gamma-in-list fixture. Gated at 1e-8 Ha as a sanity bound only.
    assert!(
        occupied <= 1e-8,
        "{label}: occupied eigenvalues {occupied:e} Ha > 1e-8"
    );
}

fn scalar_matrix(ranks: &[usize], cell_names: &[&str]) {
    for (name, cell) in cells() {
        if !cell_names.contains(&name) {
            continue;
        }
        for (set, points) in [
            ("gamma-in-list", gamma_in_list()),
            ("mp222-shifted", monkhorst_pack([2; 3], [0.5; 3]).unwrap()),
        ] {
            for strict in [true, false] {
                let mut one_rank: Option<(Outcome, Outcome)> = None;
                for &p in ranks {
                    let (off, on) = dft_mpi::simulate(p, |comm| {
                        let off = scalar(comm, &cell, &points, strict, "off");
                        let on = scalar(comm, &cell, &points, strict, "on");
                        let label = format!(
                            "{name} {set} {} P={p}",
                            if strict { "strict" } else { "standard" }
                        );
                        check(&label, strict, comm, &off, &on);
                        (off, on)
                    })
                    .swap_remove(0);
                    // Rank-count independence: delivered energies within the
                    // class; occupied eigenvalues are first order in the
                    // self-consistency error, so their rank spread is
                    // reported for both arms (the full-recompute arm has the
                    // same spread) instead of gated at the energy class.
                    if let Some((first_off, first_on)) = &one_rank {
                        let bound = if strict { STRICT_HA } else { STANDARD_HA };
                        let (e_on, occ_on, _) = diff(first_on, &on);
                        let (e_off, occ_off, _) = diff(first_off, &off);
                        println!(
                            "GK10 {name} {set} {} P={p} vs P=1: energy on {e_on:.2e} off {e_off:.2e}, occupied eigenvalues on {occ_on:.2e} off {occ_off:.2e} Ha",
                            if strict { "strict" } else { "standard" }
                        );
                        assert!(e_on <= bound, "{name} {set} P={p} vs P=1: {e_on:e} Ha");
                    } else {
                        one_rank = Some((off, on));
                    }
                }
            }
        }
    }
}

#[test]
fn border_gram_scf_agrees_on_skewed_k_cell() {
    scalar_matrix(&[1, 3], &["skewed"]);
}

#[test]
#[ignore = "full GK10 matrix; run in release"]
fn border_gram_scf_agrees_on_k_generality_set() {
    scalar_matrix(&[1, 3, 8], &["cubic", "skewed", "elongated"]);
}

#[test]
fn border_gram_fd_and_spin_k_runs_agree() {
    for spin in [false, true] {
        for strict in [true, false] {
            for p in [1, 3] {
                dft_mpi::simulate(p, |comm| {
                    let off = electronic_run(comm, spin, strict, "off");
                    let on = electronic_run(comm, spin, strict, "on");
                    let label = format!(
                        "{} {} P={p}",
                        if spin { "spin" } else { "fd" },
                        if strict { "strict" } else { "standard" }
                    );
                    check(&label, strict, comm, &off, &on);
                });
            }
        }
    }
}

/// Mac observation only (not an assertion beyond the energy class): a k
/// synthetic of B12-like shape (12 atoms, 22 bands, 24 electrons, skewed
/// cell), standard class, one run per arm. Environment: GK10_ECUT (Ha,
/// default 12), GK10_K (generic k count, default 2), GK10_RANKS (default 1).
/// cargo test --release -p dft-engine --test refscf_k_border_gram timing -- --ignored --nocapture
#[test]
#[ignore = "Mac timing observation; run in release"]
fn timing() {
    let env = |k: &str, d: f64| std::env::var(k).map_or(d, |v| v.parse().unwrap());
    let ecut = env("GK10_ECUT", 12.);
    let nk = env("GK10_K", 2.) as usize;
    let ranks = env("GK10_RANKS", 1.) as usize;
    let pp = common::potential_with_atomic_start(2., 1.8, true, 401, "SLA PZ");
    let cell =
        Cell::from_user_rows_bohr([[9.4, 0., 0.], [1.3, 9.1, 0.], [0.6, -0.8, 13.2]]).unwrap();
    let atoms: Vec<_> = (0..12)
        .map(|i| {
            let f = i as f64;
            Atom::new(
                "Xx",
                [
                    (0.07 + 0.31 * f).fract(),
                    (0.11 + 0.53 * f).fract(),
                    (0.05 + 0.083 * f).fract(),
                ],
            )
        })
        .collect();
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let points: Vec<_> = (0..nk)
        .map(|i| KPoint {
            fractional: [0.13 + 0.1 * i as f64, -0.21, 0.07 + 0.05 * i as f64],
            weight: 1. / nk as f64,
        })
        .collect();
    let options = ScfOptions {
        cutoff_hartree: ecut,
        bands: Some(22),
        ..options()
    };
    let mut energies = Vec::new();
    for border in ["off", "on", "off", "on"] {
        let out = dft_mpi::simulate(ranks, |comm| {
            let config = RunConfig::from_pairs(&arm(border)).unwrap();
            let _scope = config.enter();
            let t = std::time::Instant::now();
            let run = k_driver::run_with(
                StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &options,
                },
                profile(false),
                &points,
                comm,
            )
            .unwrap();
            let wall = t.elapsed().as_secs_f64();
            let solve: f64 = run.records.iter().map(|r| r.solve_s).sum();
            let columns: usize = run.records.iter().map(|r| r.h_columns).sum();
            let iterations: usize = run.records.iter().map(|r| r.davidson_iterations).sum();
            (
                run.outcome.result.total_hartree(),
                wall,
                solve,
                run.records.len(),
                columns,
                iterations,
            )
        });
        let (e, wall, solve, maps, columns, iterations) = out[0];
        println!(
            "GK10 timing ecut={ecut} k={nk} P={ranks} border={border}: E={e:.15e} Ha wall {wall:.3} s solve {solve:.3} s maps {maps} H columns {columns} Davidson iterations {iterations}"
        );
        energies.push(e);
    }
    for e in &energies {
        assert!((e - energies[0]).abs() <= STANDARD_HA);
    }
}
