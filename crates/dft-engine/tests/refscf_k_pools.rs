//! GK9 k-point pools (`BAANDI_REFSCF_K_POOLS`, `BAANDI_REFSCF_K_NPOOL`, ADR
//! 0041 option D): stage 1 (one rank per pool) and GK9B stage 2 (pools of
//! several ranks on sub-communicators). Acceptance on simulated ranks
//! (dft_mpi::simulate), synthetic pseudopotentials, never QE evidence:
//! * one rank: pools are never used; on and off are bit-identical;
//! * npool 1 is the distributed layout bit for bit at every rank count;
//! * P = 2, 3 (8 in the ignored matrix) with npool = P (`on`), 2 and 4
//!   (including unequal pool sizes, P = 3 npool 2) and `auto` against P = 1:
//!   delivered energy within strict 1e-10 Ha (strict-threshold runs) and
//!   standard 5e-9 Ha, same k count, map counts equal or within the
//!   distributed layout's own rank-count spread (printed);
//! * a fixed (P, npool) is bit-reproducible (two runs);
//! * generality set: cubic, skewed, elongated and slab-like cells, Gamma in an
//!   explicit list (3 k: fewer k than ranks gives pools of several ranks) and
//!   a shifted mesh, Fermi-Dirac and fixed-spin k runs;
//! * a pooled checkpoint restarts bit-identically at a fixed (P, npool).
//!
//! Full matrix: cargo test --release -p dft-engine --test refscf_k_pools -- --include-ignored --nocapture
#[allow(dead_code)]
mod common;
use dft_engine::{
    Functional, GammaReal, ScfOptions, SpeciesPotential, StandardProfile, Verification,
    collinear::SpinSector,
    refscf::{
        checkpoint::Identity,
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
use dft_mpi::{Collective, MpiError, SimRank};
use dft_numerics::switches::{Id, RunConfig};
use std::sync::atomic::{AtomicUsize, Ordering};

const STRICT_HA: f64 = 1e-10;
const STANDARD_HA: f64 = 5e-9;

fn cells() -> [(&'static str, Cell, [[f64; 3]; 2]); 4] {
    let generic = [[0.13, 0.21, 0.32], [0.61, 0.55, 0.71]];
    [
        (
            "cubic",
            Cell::from_user_rows_bohr([[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]]).unwrap(),
            generic,
        ),
        (
            "skewed",
            Cell::from_user_rows_bohr([[7.2, 0., 0.], [1.1, 6.8, 0.], [0.4, -0.7, 7.5]]).unwrap(),
            generic,
        ),
        (
            "elongated",
            Cell::from_user_rows_bohr([[5., 0., 0.], [0., 5., 0.], [0., 0., 11.]]).unwrap(),
            generic,
        ),
        (
            // Two atoms in a thin layer with vacuum along c.
            "slab",
            Cell::from_user_rows_bohr([[5.5, 0., 0.], [0., 5.5, 0.], [0., 0., 15.]]).unwrap(),
            [[0.1, 0.2, 0.46], [0.6, 0.7, 0.55]],
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

/// Strict arm: 1e-12 Ry (5e-13 Ha) so delivered energies resolve the 1e-10 Ha
/// class; standard arm: the default 1e-8 Ry (5e-9 Ha).
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

#[derive(Clone, Debug)]
struct Outcome {
    energy: f64,
    converged: bool,
    maps: usize,
    h_columns: usize,
    /// Per (channel, k) block: eigenvalues and the count of occupied bands.
    blocks: Vec<(Vec<f64>, usize)>,
    /// Bits of the delivered density, eigenvalues and orbitals.
    bits: Vec<u64>,
    /// The layout the switch took ("pools" / "distributed"), if consulted.
    path: Option<String>,
    /// The delivered blocks were lazy (materialized before `bits`).
    lazy: bool,
}

fn bits_of<'a>(values: impl Iterator<Item = &'a f64>) -> Vec<u64> {
    values.map(|v| v.to_bits()).collect()
}

/// `off`, `auto`, `on`, or `npool=N` (`auto` with `refscf_k_npool = N`).
fn config(pools: &str) -> RunConfig {
    match pools.strip_prefix("npool=") {
        Some(n) => {
            RunConfig::from_pairs(&[("refscf_k_pools", "auto"), ("refscf_k_npool", n)]).unwrap()
        }
        None => RunConfig::from_pairs(&[("refscf_k_pools", pools)]).unwrap(),
    }
}

/// Pool count an arm asks for at `p` ranks and `nk` k points (the engine's
/// rule for `on` and forced counts), and the layout note it must produce.
fn expected(pools: &str, p: usize, nk: usize) -> (usize, String) {
    let npool = match pools.strip_prefix("npool=") {
        Some(n) => n.parse::<usize>().unwrap().clamp(1, p.min(nk)),
        None if pools == "on" => (1..=p.min(nk))
            .filter(|&d| p.is_multiple_of(d))
            .max()
            .unwrap(),
        None => 1,
    };
    let path = if p < 2 || npool < 2 {
        "distributed".to_string()
    } else if npool == p {
        "pools".to_string()
    } else {
        format!("pools:npool={npool}")
    };
    (npool, path)
}

/// Pool arms at `p` ranks: one rank per pool, the distributed layout as one
/// pool, and pools of several ranks (unequal sizes at P = 3).
fn arms(p: usize) -> Vec<String> {
    let mut out = vec!["on".to_string(), "npool=1".into(), "auto".into()];
    for n in [2usize, 4] {
        if n < p {
            out.push(format!("npool={n}"));
        }
    }
    out
}

fn scalar<K: Collective>(
    comm: &K,
    cell: &Cell,
    positions: [[f64; 3]; 2],
    points: &[KPoint],
    strict: bool,
    pools: &str,
) -> Outcome {
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    let atoms = [Atom::new("Xx", positions[0]), Atom::new("Xx", positions[1])];
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = options();
    let config = config(pools);
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
    let r = &run.outcome.result;
    let mut bits = bits_of(r.density.iter());
    for k in &r.k_points {
        bits.extend(bits_of(k.solution.eigenvalues.iter()));
        bits.extend(
            k.solution
                .eigenvectors
                .as_slice()
                .iter()
                .flat_map(|z| [z.re.to_bits(), z.im.to_bits()]),
        );
    }
    bits.push(r.energy.total_hartree().to_bits());
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
        bits,
        path: config.path_note(Id::REFSCF_K_POOLS).map(|n| n.path),
        lazy: false,
    }
}

fn electronic_run<K: Collective>(comm: &K, spin: bool, strict: bool, pools: &str) -> Outcome {
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
    let config = config(pools);
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
    let mut run = run;
    let lazy = !run.orbitals_delivered();
    run.materialize_orbitals().unwrap();
    let mut bits: Vec<u64> = run.density.iter().flat_map(|d| bits_of(d.iter())).collect();
    for b in run.blocks.iter().flatten() {
        bits.extend(bits_of(b.eigenvalues_ha.iter()));
        bits.extend(bits_of(b.occupations.iter()));
        bits.extend(
            b.orbitals
                .as_slice()
                .iter()
                .flat_map(|z| [z.re.to_bits(), z.im.to_bits()]),
        );
    }
    bits.push(run.thermodynamics.free_energy_ha.to_bits());
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
        bits,
        path: config.path_note(Id::REFSCF_K_POOLS).map(|n| n.path),
        lazy,
    }
}

/// (|energy difference|, max occupied eigenvalue difference), Ha.
fn diff(a: &Outcome, b: &Outcome) -> (f64, f64) {
    assert_eq!(a.blocks.len(), b.blocks.len(), "same k count");
    let mut occupied = 0_f64;
    for ((x, nx), (y, ny)) in a.blocks.iter().zip(&b.blocks) {
        assert_eq!((x.len(), nx), (y.len(), ny));
        for (p, q) in x.iter().zip(y).take(*nx) {
            occupied = occupied.max((p - q).abs());
        }
    }
    ((a.energy - b.energy).abs(), occupied)
}

/// Spread of the distributed layout (pools off) against one rank over every
/// rank count from 2 to the largest of a matrix row: (max occupied eigenvalue difference, Ha;
/// max map-count difference). Today's own rank-count sensitivity of a fixture.
#[derive(Clone, Copy, Default)]
struct Spread {
    occupied: f64,
    maps: usize,
}

/// Pools on at `p` ranks against one rank: the accuracy class (the gate),
/// the k count, the layout taken, and the map count and occupied eigenvalues
/// as sanity checks. Map counts must be equal (+-1), unless the distributed
/// layout itself moves the map count of this fixture by more over the tested
/// rank counts (`spread`, the justification, printed). Occupied eigenvalues
/// are first order in the self-consistency error: bound 1e-8 Ha as in the
/// GK10 rank-count checks, or 10 times the distributed layout's own spread.
fn check_class(
    label: &str,
    strict: bool,
    one: &Outcome,
    pooled: &Outcome,
    spread: Option<Spread>,
    expect: Option<&str>,
) {
    let bound = if strict { STRICT_HA } else { STANDARD_HA };
    let (energy, occupied) = diff(one, pooled);
    let spread = spread.unwrap_or_default();
    println!(
        "GK9 {label}: E={:.15e} Ha |pools-P1| energy {energy:.2e} occupied eigenvalues {occupied:.2e} Ha; maps P1 {} pools {}; H columns P1 {} pools {}; distributed-layout spread over the row: occupied {:.2e} Ha, maps {}; layout {:?}",
        pooled.energy,
        one.maps,
        pooled.maps,
        one.h_columns,
        pooled.h_columns,
        spread.occupied,
        spread.maps,
        pooled.path
    );
    assert!(one.converged && pooled.converged, "{label}");
    if let Some(path) = expect {
        assert_eq!(pooled.path.as_deref(), Some(path), "{label}");
    }
    assert!(energy <= bound, "{label}: energy {energy:e} Ha > {bound:e}");
    let sanity = 1e-8_f64.max(10. * spread.occupied);
    assert!(
        occupied <= sanity,
        "{label}: occupied eigenvalues {occupied:e} Ha > {sanity:e}"
    );
    assert!(
        pooled.maps.abs_diff(one.maps) <= spread.maps.max(1),
        "{label}: maps {} vs {} (distributed spread {})",
        pooled.maps,
        one.maps,
        spread.maps
    );
}

fn on_rank0<R: Clone + Send>(p: usize, f: impl Fn(&SimRank) -> R + Sync) -> R {
    dft_mpi::simulate(p, f).swap_remove(0)
}

fn scalar_matrix(ranks: &[usize], cell_names: &[&str]) {
    for (name, cell, positions) in cells() {
        if !cell_names.contains(&name) {
            continue;
        }
        for (set, points) in [
            ("gamma-in-list", gamma_in_list()),
            ("mp222-shifted", monkhorst_pack([2; 3], [0.5; 3]).unwrap()),
        ] {
            for strict in [true, false] {
                let class = if strict { "strict" } else { "standard" };
                let one = on_rank0(1, |c| scalar(c, &cell, positions, &points, strict, "on"));
                // Today's rank-count sensitivity of the fixture: the
                // distributed layout at every rank count 2..=max(ranks).
                let mut spread = Spread::default();
                for p in 2..=ranks.iter().copied().max().unwrap_or(1) {
                    let off = on_rank0(p, |c| scalar(c, &cell, positions, &points, strict, "off"));
                    spread.occupied = spread.occupied.max(diff(&one, &off).1);
                    spread.maps = spread.maps.max(off.maps.abs_diff(one.maps));
                }
                for &p in ranks {
                    let off = on_rank0(p, |c| scalar(c, &cell, positions, &points, strict, "off"));
                    for arm in arms(p) {
                        let pooled =
                            on_rank0(p, |c| scalar(c, &cell, positions, &points, strict, &arm));
                        let (npool, path) = expected(&arm, p, points.len());
                        if npool == 1 && arm != "auto" {
                            assert_eq!(pooled.bits, off.bits, "{name} {set} P={p} {arm}: npool 1");
                        }
                        check_class(
                            &format!(
                                "{name} {set} {class} P={p} nk={} {arm}{}",
                                points.len(),
                                if arm == "auto" {
                                    String::new()
                                } else {
                                    format!(" (npool {npool})")
                                }
                            ),
                            strict,
                            &one,
                            &pooled,
                            Some(spread),
                            (arm != "auto").then_some(path.as_str()),
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn pools_agree_with_one_rank_on_skewed_k_cell() {
    scalar_matrix(&[2, 3], &["skewed"]);
}

#[test]
#[ignore = "full GK9 matrix; run in release"]
fn pools_agree_on_k_generality_set() {
    scalar_matrix(&[2, 3, 8], &["cubic", "skewed", "elongated", "slab"]);
}

#[test]
fn pools_fd_and_spin_k_runs_agree() {
    for spin in [false, true] {
        for strict in [true, false] {
            let one = on_rank0(1, |c| electronic_run(c, spin, strict, "on"));
            for p in [2, 3] {
                for arm in arms(p) {
                    let pooled = on_rank0(p, |c| electronic_run(c, spin, strict, &arm));
                    check_class(
                        &format!(
                            "{} {} P={p} nk=8 {arm}",
                            if spin { "spin" } else { "fd" },
                            if strict { "strict" } else { "standard" }
                        ),
                        strict,
                        &one,
                        &pooled,
                        None,
                        (arm != "auto").then_some(expected(&arm, p, 8).1.as_str()),
                    );
                }
            }
        }
    }
}

#[test]
#[ignore = "FD/spin at eight ranks; run in release"]
fn pools_fd_and_spin_k_runs_agree_at_eight_ranks() {
    for spin in [false, true] {
        for strict in [true, false] {
            let one = on_rank0(1, |c| electronic_run(c, spin, strict, "on"));
            for arm in arms(8) {
                let pooled = on_rank0(8, |c| electronic_run(c, spin, strict, &arm));
                check_class(
                    &format!(
                        "{} {} P=8 nk=8 {arm}",
                        if spin { "spin" } else { "fd" },
                        if strict { "strict" } else { "standard" }
                    ),
                    strict,
                    &one,
                    &pooled,
                    None,
                    (arm != "auto").then_some(expected(&arm, 8, 8).1.as_str()),
                );
            }
        }
    }
}

/// One rank never pools (bit-identical on/off); a fixed rank count is
/// bit-reproducible with pools on (two runs, every rank's delivery equal).
#[test]
fn pools_one_rank_identical_and_fixed_ranks_reproducible() {
    let (_, cell, positions) = cells().into_iter().nth(1).unwrap();
    let points = monkhorst_pack([2; 3], [0.5; 3]).unwrap();
    let off = on_rank0(1, |c| scalar(c, &cell, positions, &points, false, "off"));
    let on = on_rank0(1, |c| scalar(c, &cell, positions, &points, false, "on"));
    assert_eq!(on.path.as_deref(), Some("distributed"));
    assert_eq!(off.bits, on.bits, "one rank: pools on must equal off");
    for p in [2, 3] {
        let runs: Vec<Vec<Outcome>> = (0..2)
            .map(|_| dft_mpi::simulate(p, |c| scalar(c, &cell, positions, &points, false, "on")))
            .collect();
        for run in &runs {
            for rank in run {
                assert_eq!(
                    rank.bits, runs[0][0].bits,
                    "P={p}: ranks deliver the same bits"
                );
            }
        }
        assert_eq!(runs[0][0].path.as_deref(), Some("pools"));
        println!(
            "GK9 reproducible P={p}: two runs bit-identical (E={:.15e} Ha, {} maps)",
            runs[0][0].energy, runs[0][0].maps
        );
    }
    for spin in [false, true] {
        for arm in ["on", "npool=2"] {
            let a = dft_mpi::simulate(3, |c| electronic_run(c, spin, false, arm));
            let b = dft_mpi::simulate(3, |c| electronic_run(c, spin, false, arm));
            for x in a.iter().chain(&b) {
                assert_eq!(
                    x.bits, a[0].bits,
                    "spin={spin} {arm}: fixed P=3 bit-reproducible"
                );
            }
        }
    }
    // Pools of several ranks (2 x 2 and unequal 2 + 1): fixed (P, npool)
    // bit-reproducible, every rank delivers the same bits.
    for (p, arm) in [(4, "npool=2"), (3, "npool=2"), (4, "auto")] {
        let a = dft_mpi::simulate(p, |c| scalar(c, &cell, positions, &points, false, arm));
        let b = dft_mpi::simulate(p, |c| scalar(c, &cell, positions, &points, false, arm));
        for x in a.iter().chain(&b) {
            assert_eq!(x.bits, a[0].bits, "P={p} {arm}: bit-reproducible");
        }
        println!(
            "GK9B reproducible P={p} {arm}: two runs bit-identical, layout {:?} (E={:.15e} Ha, {} maps)",
            a[0].path, a[0].energy, a[0].maps
        );
    }
}

/// GK18: a pooled run whose caller does not read the orbitals delivers lazy k
/// blocks; after `materialize_orbitals` every delivered bit (density,
/// eigenvalues, orbitals, energy, `orbitals`) equals the eager delivery, for
/// one rank per pool and pools of several ranks (2 x 2, unequal 2 + 1), and
/// `M2_LAZY_RESULT=off` keeps the eager delivery.
#[test]
fn pooled_lazy_delivery_materializes_the_eager_bits() {
    use std::sync::Arc;
    let (_, cell, positions) = cells().into_iter().nth(1).unwrap();
    let points = monkhorst_pack([2; 3], [0.5; 3]).unwrap();
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    let atoms = [Atom::new("Xx", positions[0]), Atom::new("Xx", positions[1])];
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = options();
    for (p, npool) in [(2usize, "2"), (3, "3"), (4, "2"), (3, "2")] {
        let ranks = dft_mpi::simulate(p, |comm: &SimRank| {
            let shared: Arc<dyn Collective + Send + Sync> = Arc::new(comm.clone());
            dft_engine::parallel::set_communicator(Some(shared));
            let mut bits = Vec::new();
            for (needed, lazy_switch) in [(true, "on"), (false, "on"), (false, "off")] {
                let config = RunConfig::from_pairs(&[
                    ("refscf_k_pools", "auto"),
                    ("refscf_k_npool", npool),
                    ("m2_lazy_result", lazy_switch),
                ])
                .unwrap();
                let _scope = config.enter();
                dft_engine::refscf::lazy::set_delivered_orbitals(needed);
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
                let mut r = run.outcome.result;
                let lazy = !needed && lazy_switch == "on";
                assert_eq!(r.orbitals_delivered(), !lazy, "P={p} npool={npool}");
                if lazy {
                    assert_eq!(r.orbitals.rows(), 0);
                    assert!(
                        r.k_points
                            .iter()
                            .all(|k| k.solution.eigenvectors.rows() == 0)
                    );
                }
                r.materialize_orbitals().unwrap();
                assert!(r.orbitals_delivered());
                let mut b = bits_of(r.density.iter());
                b.push(r.energy.total_hartree().to_bits());
                for k in &r.k_points {
                    b.extend(bits_of(k.solution.eigenvalues.iter()));
                    b.extend(
                        k.solution
                            .eigenvectors
                            .as_slice()
                            .iter()
                            .flat_map(|z| [z.re.to_bits(), z.im.to_bits()]),
                    );
                }
                b.extend(
                    r.orbitals
                        .as_slice()
                        .iter()
                        .flat_map(|z| [z.re.to_bits(), z.im.to_bits()]),
                );
                bits.push(b);
            }
            dft_engine::refscf::lazy::set_delivered_orbitals(true);
            dft_engine::parallel::set_communicator(None);
            bits
        });
        for rank in &ranks {
            for b in rank {
                assert_eq!(b, &ranks[0][0], "P={p} npool={npool}: lazy = eager bits");
            }
        }
        println!("GK18 lazy pooled delivery P={p} npool={npool}: bits equal eager on every rank");
    }
}

/// GK18B: the spin/FD driver delivers lazy k blocks when the caller does not
/// read the orbitals; after `Run::materialize_orbitals` every delivered bit
/// equals the eager delivery (FD and fixed spin; one rank per pool and pools of
/// several ranks).
#[test]
fn electronic_lazy_delivery_materializes_the_eager_bits() {
    use std::sync::Arc;
    for spin in [false, true] {
        for (p, arm) in [(2usize, "on"), (4, "npool=2")] {
            let ranks = dft_mpi::simulate(p, |comm: &SimRank| {
                let shared: Arc<dyn Collective + Send + Sync> = Arc::new(comm.clone());
                dft_engine::parallel::set_communicator(Some(shared));
                let mut out = Vec::new();
                for needed in [true, false] {
                    dft_engine::refscf::lazy::set_delivered_orbitals(needed);
                    out.push(electronic_run(comm, spin, false, arm));
                }
                dft_engine::refscf::lazy::set_delivered_orbitals(true);
                dft_engine::parallel::set_communicator(None);
                out
            });
            for rank in &ranks {
                assert!(!rank[0].lazy && rank[1].lazy, "spin={spin} P={p} {arm}");
                assert_eq!(rank[0].bits, ranks[0][0].bits);
                assert_eq!(
                    rank[1].bits, ranks[0][0].bits,
                    "spin={spin} P={p} {arm}: lazy = eager"
                );
            }
            println!(
                "GK18B electronic lazy spin={spin} P={p} {arm}: bits equal eager on every rank"
            );
        }
    }
}

/// npool 1 (forced, and `auto` when it picks one pool) is the distributed
/// layout bit for bit; fewer k than ranks gives pools of several ranks.
#[test]
fn npool_one_is_the_distributed_layout_and_few_k_use_sub_pools() {
    let (_, cell, positions) = cells().into_iter().nth(1).unwrap();
    let points = gamma_in_list();
    for p in [2usize, 3, 4] {
        let off = on_rank0(p, |c| scalar(c, &cell, positions, &points, false, "off"));
        let one = on_rank0(p, |c| {
            scalar(c, &cell, positions, &points, false, "npool=1")
        });
        assert_eq!(one.path.as_deref(), Some("distributed"));
        assert_eq!(off.bits, one.bits, "P={p}: npool 1");
    }
    let sub = on_rank0(4, |c| scalar(c, &cell, positions, &points, false, "on"));
    assert_eq!(
        sub.path.as_deref(),
        Some("pools:npool=2"),
        "3 k on 4 ranks: 2 pools of 2"
    );
}

/// Collectives on the world communicator per map, pools off and on (P = 3,
/// 8 k points). With pools every per-k operation uses the serial collective.
#[test]
fn pools_collective_counts_per_map() {
    struct Counting<'a> {
        inner: &'a SimRank,
        calls: [AtomicUsize; 5],
        values: AtomicUsize,
    }
    impl Counting<'_> {
        fn count(&self, kind: usize, values: usize) {
            self.calls[kind].fetch_add(1, Ordering::Relaxed);
            self.values.fetch_add(values, Ordering::Relaxed);
        }
    }
    impl Collective for Counting<'_> {
        fn rank(&self) -> usize {
            self.inner.rank()
        }
        fn size(&self) -> usize {
            self.inner.size()
        }
        fn barrier(&self) -> Result<(), MpiError> {
            self.count(4, 0);
            self.inner.barrier()
        }
        fn allgatherv_f64(&self, a: &[f64], n: &[usize], b: &mut [f64]) -> Result<(), MpiError> {
            self.count(0, b.len());
            self.inner.allgatherv_f64(a, n, b)
        }
        fn tree_sum_f64(&self, a: &[f64], b: &mut [f64]) -> Result<(), MpiError> {
            self.count(1, a.len());
            self.inner.tree_sum_f64(a, b)
        }
        fn bcast_f64(&self, a: &mut [f64], root: usize) -> Result<(), MpiError> {
            self.count(2, a.len());
            self.inner.bcast_f64(a, root)
        }
        fn alltoall_f64(&self, a: &[f64], n: usize, b: &mut [f64]) -> Result<(), MpiError> {
            self.count(3, a.len());
            self.inner.alltoall_f64(a, n, b)
        }
        fn split(&self, color: usize, key: usize) -> Result<dft_mpi::Comm, MpiError> {
            self.inner.split(color, key)
        }
    }
    let (_, cell, positions) = cells().into_iter().nth(1).unwrap();
    let points = monkhorst_pack([2; 3], [0.5; 3]).unwrap();
    // Totals include the start and the delivery; the marginal count per map is
    // the difference between a strict (more maps) and a standard run.
    for pools in ["off", "on"] {
        let mut rows = Vec::new();
        for strict in [false, true] {
            rows.push(on_rank0(3, |c| {
                let counting = Counting {
                    inner: c,
                    calls: Default::default(),
                    values: AtomicUsize::new(0),
                };
                let out = scalar(&counting, &cell, positions, &points, strict, pools);
                (
                    counting.calls.map(|a| a.into_inner()),
                    counting.values.into_inner(),
                    out.maps,
                )
            }));
        }
        for (counts, values, maps) in &rows {
            println!(
                "GK9 collectives P=3 nk=8 pools={pools}: {maps} maps, {} world collectives (allgatherv {} tree_sum {} bcast {} alltoall {} barrier {}), {values} values; start and delivery included",
                counts.iter().sum::<usize>(),
                counts[0],
                counts[1],
                counts[2],
                counts[3],
                counts[4],
            );
        }
        let ((c0, v0, m0), (c1, v1, m1)) = (&rows[0], &rows[1]);
        let maps = (m1 - m0) as f64;
        let per = |i: usize| (c1[i] as f64 - c0[i] as f64) / maps;
        println!(
            "GK9 collectives per map P=3 nk=8 pools={pools}: {:.1} (allgatherv {:.1} tree_sum {:.1} bcast {:.1} alltoall {:.1}), {:.0} values per map",
            (c1.iter().sum::<usize>() as f64 - c0.iter().sum::<usize>() as f64) / maps,
            per(0),
            per(1),
            per(2),
            per(3),
            (*v1 as f64 - *v0 as f64) / maps
        );
    }
}

/// A pooled run checkpoints and restarts bit-identically at a fixed rank
/// count (each rank saves its own whole k points).
#[test]
fn pooled_checkpoint_restarts_bit_identically() {
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0., 5., 0.], [0., 0., 5.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    let pp = common::potential_with_atomic_start(2., 1., true, 161, "SLA PZ");
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: 2.5,
        bands: Some(4),
        gamma_real: GammaReal::Off,
        max_iterations: 30,
        verification: Verification::FixedPoint,
        ..Default::default()
    };
    let points = monkhorst_pack([2; 3], [0.5; 3]).unwrap();
    let identity = |ranks| Identity {
        input_sha256: "1".repeat(64),
        upf_sha256: vec!["2".repeat(64)],
        binary_sha256: "3".repeat(64),
        platform: "synthetic-test".into(),
        threads: 1,
        ranks,
    };
    for (ranks, arm, layout) in [(2usize, "on", "pools"), (4, "npool=2", "pools:npool=2")] {
        dft_mpi::simulate(ranks, |comm| {
            let config = config(arm);
            let _scope = config.enter();
            let build = || {
                k_driver::stepper_with(
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
            };
            let energy = |run: &dft_engine::refscf::driver::Run| {
                let r = &run.outcome.result;
                let mut bits = bits_of(r.density.iter());
                bits.push(r.energy.total_hartree().to_bits());
                bits
            };
            let uninterrupted = build().run_to_end().unwrap();
            assert!(uninterrupted.outcome.converged);
            assert_eq!(
                config
                    .path_note(Id::REFSCF_K_POOLS)
                    .map(|n| n.path)
                    .as_deref(),
                Some(layout)
            );
            let mut s = build();
            s.initialize().unwrap();
            loop {
                let r = s.step().unwrap();
                assert!(!s.is_done());
                if r.iteration >= 2 && r.safe_point {
                    break;
                }
            }
            let cp = s.checkpoint(identity(ranks), comm.rank()).unwrap();
            let captured = s.run_to_end().unwrap();
            assert_eq!(
                energy(&captured),
                energy(&uninterrupted),
                "capture changed trajectory"
            );
            let mut restored = build();
            restored
                .restore(&cp, &identity(ranks), comm.rank())
                .unwrap();
            let resumed = restored.run_to_end().unwrap();
            assert_eq!(energy(&resumed), energy(&uninterrupted), "pooled restart");
            if comm.rank() == 0 {
                println!(
                    "GK9 pooled checkpoint P={ranks} {arm} at map {}: restart bit-identical ({} maps)",
                    cp.iteration,
                    resumed.records.len()
                );
            }
        });
    }
}

/// Mac observation only (not acceptance): a k synthetic of B12/GeTe-like
/// shape (12 atoms, 22 bands, skewed cell, GK9_K generic k points, default
/// 16), standard class, pools off and on at GK9_RANKS simulated ranks
/// (default 8), one run per arm (GK9_ARMS, default `off,on`). GK9_ECUT (Ha,
/// default 8).
/// cargo test --release -p dft-engine --test refscf_k_pools timing -- --ignored --nocapture
#[test]
#[ignore = "Mac timing observation; run in release"]
fn timing() {
    let env = |k: &str, d: f64| std::env::var(k).map_or(d, |v| v.parse().unwrap());
    let ecut = env("GK9_ECUT", 8.);
    let nk = env("GK9_K", 16.) as usize;
    let ranks = env("GK9_RANKS", 8.) as usize;
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
            fractional: [
                0.13 + 0.37 * i as f64 / nk as f64,
                -0.21 + 0.11 * (i % 3) as f64,
                0.07 + 0.05 * (i % 5) as f64,
            ],
            weight: 1. / nk as f64,
        })
        .collect();
    let options = ScfOptions {
        cutoff_hartree: ecut,
        bands: Some(22),
        ..options()
    };
    let mut energies = Vec::new();
    // GK9_ARMS: comma-separated arms (`off`, `on`, `auto`, `npool=N`).
    let arms = std::env::var("GK9_ARMS").unwrap_or_else(|_| "off,on".into());
    for pools in arms.split(',') {
        let out = dft_mpi::simulate(ranks, |comm| {
            let config = config(pools);
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
            (
                run.outcome.result.total_hartree(),
                wall,
                solve,
                run.records.len(),
                columns,
                run.start_s,
            )
        });
        let (e, wall, solve, maps, columns, start) = out[0];
        println!(
            "GK9 timing ecut={ecut} k={nk} P={ranks} pools={pools}: E={e:.15e} Ha wall {wall:.3} s (start {start:.3} s, solve {solve:.3} s) maps {maps} H columns {columns}"
        );
        energies.push(e);
    }
    for e in &energies {
        assert!((e - energies[0]).abs() <= STANDARD_HA);
    }
}
