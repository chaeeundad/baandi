//! GF1 (spin/FD k path): ports of the scalar k driver's kept-HX items to the
//! electronic (spin/FD) driver. Synthetic pseudopotentials, simulated ranks,
//! never QE evidence.
//! * `SPIN_FD_START_KEPT_HX` (GK10B start) and `SPIN_FD_REDO_KEPT_HX` (GK10B
//!   redo): defaults against both off, at P = 1 and pooled P = 2, 3 (8 in
//!   the ignored matrix): delivered free energy within strict 1e-10 Ha
//!   (strict-threshold runs) or standard 5e-9 Ha, same convergence, map count
//!   within one, and map 1 applies H to fewer columns;
//! * `K_KEPT_SLABS` on the electronic map: bit-identical to slabs off;
//! * automatic FD band growth with the ports on and pooled.
//!
//! Full matrix: cargo test --release -p dft-engine --test refscf_fd_kept -- --include-ignored --nocapture
#[allow(dead_code)]
mod common;
use dft_engine::{
    Functional, GammaReal, ScfOptions, SpeciesPotential, StandardProfile, Verification,
    collinear::SpinSector,
    refscf::{
        electronic::{self, BandGrowth, Policy},
        profile::Profile,
    },
    refstart::StartInput,
};
use dft_geometry::{
    Atom, Cell,
    kpoints::{KPoint, monkhorst_pack},
};
use dft_mpi::{Collective, SimRank};
use dft_numerics::switches::RunConfig;

const STRICT_HA: f64 = 1e-10;
const STANDARD_HA: f64 = 5e-9;
const PORTS_OFF: &[(&str, &str)] = &[
    ("spin_fd_start_kept_hx", "off"),
    ("spin_fd_redo_kept_hx", "off"),
];

#[derive(Clone, Copy, Debug, PartialEq)]
enum Case {
    /// One atom, skewed near-cubic cell, FD, MP 2x2x2 shifted (8 k).
    Fd,
    /// Two atoms (two electrons each), elongated cell, FD, Gamma in an
    /// explicit list (3 k); map 1 is redone.
    FdElongated,
    /// One atom, skewed cell, fixed spin sector, MP 2x2x2 shifted.
    Spin,
}

#[derive(Clone, Debug)]
struct Outcome {
    free_energy: f64,
    converged: bool,
    maps: usize,
    map1_h: usize,
    h_columns: usize,
    redo: bool,
    hits: usize,
    growth_events: usize,
    bits: Vec<u64>,
}

fn bits_of<'a>(values: impl Iterator<Item = &'a f64>) -> Vec<u64> {
    values.map(|v| v.to_bits()).collect()
}

fn run<K: Collective>(
    comm: &K,
    case: Case,
    strict: bool,
    pools: &str,
    pairs: &[(&str, &str)],
    growth: bool,
) -> Outcome {
    let spin = case == Case::Spin;
    let (cell, atoms, points) = match case {
        Case::Fd | Case::Spin => (
            Cell::from_user_rows_bohr([[5., 0., 0.], [0.2, 5., 0.], [0., 0.1, 5.]]).unwrap(),
            vec![Atom::new("Xx", [0.13, 0.21, 0.32])],
            monkhorst_pack([2; 3], [0.5; 3]).unwrap(),
        ),
        Case::FdElongated => (
            Cell::from_user_rows_bohr([[4.6, 0., 0.], [0., 4.8, 0.], [0.3, 0., 9.5]]).unwrap(),
            vec![
                Atom::new("Xx", [0.1, 0.2, 0.15]),
                Atom::new("Xx", [0.55, 0.6, 0.62]),
            ],
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
            ],
        ),
    };
    // The two-atom case carries two electrons per atom (an atomic start
    // close enough to the fixed point that map 1 is redone, as on real metals).
    let pp = common::potential_with_atomic_start(
        if case == Case::FdElongated { 2. } else { 1. },
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
        bands: if growth {
            None
        } else {
            Some(if spin { 4 } else { 12 })
        },
        max_iterations: 100,
        verification: Verification::Off,
        gamma_real: GammaReal::Off,
        standard: Some(StandardProfile {
            threshold_ha: STANDARD_HA,
        }),
        ..ScfOptions::strict_reference()
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
    let mut all = vec![("refscf_k_pools", pools)];
    all.extend_from_slice(pairs);
    let config = RunConfig::from_pairs(&all).unwrap();
    let _scope = config.enter();
    let input = StartInput {
        cell: &cell,
        atoms: &atoms,
        potentials: &potentials,
        options: &options,
    };
    let profile = Profile {
        threshold_ry: if strict { 1e-12 } else { 1e-8 },
        ..Profile::default()
    };
    let mut r = if growth {
        electronic::stepper_growing(input, profile, &points, policy, BandGrowth::From(2), comm)
            .run_to_end()
    } else {
        electronic::run_with(input, profile, &points, policy, comm)
    }
    .unwrap();
    r.materialize_orbitals().unwrap();
    let mut bits: Vec<u64> = r.density.iter().flat_map(|d| bits_of(d.iter())).collect();
    for b in r.blocks.iter().flatten() {
        bits.extend(bits_of(b.eigenvalues_ha.iter()));
        bits.extend(bits_of(b.occupations.iter()));
        bits.extend(
            b.orbitals
                .as_slice()
                .iter()
                .flat_map(|z| [z.re.to_bits(), z.im.to_bits()]),
        );
    }
    bits.push(r.thermodynamics.free_energy_ha.to_bits());
    let h = |i: &electronic::Iteration| i.calls.iter().map(|(_, _, c)| c.h_columns).sum::<usize>();
    Outcome {
        free_energy: r.thermodynamics.free_energy_ha,
        converged: r.converged,
        maps: r.history.len(),
        map1_h: h(&r.history[0]),
        h_columns: r.history.iter().map(h).sum(),
        redo: r.history[0].first_iteration_redo,
        hits: r.kept_hx.iter().map(|k| k.cache_hits).sum(),
        growth_events: r.band_growth.as_ref().map_or(0, |g| g.events.len()),
        bits,
    }
}

fn on_rank0<R: Clone + Send>(p: usize, f: impl Fn(&SimRank) -> R + Sync) -> R {
    dft_mpi::simulate(p, f).swap_remove(0)
}

fn pools_arm(p: usize) -> &'static str {
    if p == 1 { "off" } else { "on" }
}

fn ports_matrix(ranks: &[usize]) {
    for case in [Case::Fd, Case::FdElongated, Case::Spin] {
        for strict in [true, false] {
            let bound = if strict { STRICT_HA } else { STANDARD_HA };
            let one = on_rank0(1, |c| run(c, case, strict, "off", &[], false));
            for &p in ranks {
                let pools = pools_arm(p);
                let off = on_rank0(p, |c| run(c, case, strict, pools, PORTS_OFF, false));
                let on = on_rank0(p, |c| run(c, case, strict, pools, &[], false));
                let d_ports = (on.free_energy - off.free_energy).abs();
                let d_one = (on.free_energy - one.free_energy).abs();
                println!(
                    "GF1 {case:?} {} P={p} pools={pools}: F={:.15e} |on-off|={d_ports:.2e} |on-P1|={d_one:.2e} Ha; maps off {} on {}; map-1 H off {} on {} (redo {}); H total off {} on {}; kept hits off {} on {}",
                    if strict { "strict" } else { "standard" },
                    on.free_energy,
                    off.maps,
                    on.maps,
                    off.map1_h,
                    on.map1_h,
                    on.redo,
                    off.h_columns,
                    on.h_columns,
                    off.hits,
                    on.hits,
                );
                assert!(off.converged && on.converged, "{case:?} P={p}");
                assert!(d_ports <= bound, "{case:?} P={p}: ports {d_ports:e}");
                assert!(d_one <= bound, "{case:?} P={p}: vs P=1 {d_one:e}");
                assert!(on.maps.abs_diff(off.maps) <= 1, "{case:?} P={p}: maps");
                assert!(on.map1_h < off.map1_h, "{case:?} P={p}: map-1 H columns");
                assert!(on.hits > off.hits, "{case:?} P={p}: kept hits");
            }
        }
    }
}

#[test]
fn ports_within_class_and_cut_map_one() {
    ports_matrix(&[1, 2, 3]);
}

#[test]
#[ignore = "eight simulated ranks; run in release"]
fn ports_within_class_at_eight_ranks() {
    ports_matrix(&[8]);
}

/// `K_KEPT_SLABS` (default off) on the electronic map: every delivered bit
/// equal to slabs off (whole budget and a partial 1 MiB budget), FD and fixed
/// spin, serial and pooled; a fixed rank count is bit-reproducible.
#[test]
fn slabs_are_bit_identical_and_runs_reproducible() {
    for case in [Case::Fd, Case::Spin] {
        for p in [1, 2] {
            let pools = pools_arm(p);
            let base = dft_mpi::simulate(p, |c| run(c, case, false, pools, &[], false));
            let again = dft_mpi::simulate(p, |c| run(c, case, false, pools, &[], false));
            for (a, b) in base.iter().zip(&again) {
                assert_eq!(a.bits, b.bits, "{case:?} P={p}: reproducible");
            }
            for budget in ["4096", "1"] {
                let slabs = dft_mpi::simulate(p, |c| {
                    run(
                        c,
                        case,
                        false,
                        pools,
                        &[("k_kept_slabs", "on"), ("k_kept_slab_mb", budget)],
                        false,
                    )
                });
                for (rank, (a, b)) in base.iter().zip(&slabs).enumerate() {
                    assert_eq!(a.bits, b.bits, "{case:?} P={p} rank {rank} budget {budget}");
                    assert_eq!(a.h_columns, b.h_columns);
                }
            }
            println!("GF1 slabs {case:?} P={p}: bits equal (budgets 4096, 1 MiB)");
        }
    }
}

/// Automatic FD band growth with the ports on: serial and pooled runs grow,
/// converge and agree with the ports-off arm within the standard class.
#[test]
fn band_growth_with_ports_and_pools() {
    for p in [1, 2] {
        let pools = pools_arm(p);
        let off = on_rank0(p, |c| run(c, Case::Fd, false, pools, PORTS_OFF, true));
        let on = on_rank0(p, |c| run(c, Case::Fd, false, pools, &[], true));
        println!(
            "GF1 growth P={p}: events off {} on {}; |dF|={:.2e} Ha; maps off {} on {}; H off {} on {}",
            off.growth_events,
            on.growth_events,
            (on.free_energy - off.free_energy).abs(),
            off.maps,
            on.maps,
            off.h_columns,
            on.h_columns
        );
        assert!(on.converged && off.converged);
        assert!(on.growth_events > 0 && off.growth_events > 0);
        assert!((on.free_energy - off.free_energy).abs() <= STANDARD_HA);
    }
}

/// Real metals (PseudoDojo Al and Au, `DFT_NC_DIR`; skipped without it): the
/// atomic start is close enough that map 1 is redone, which exercises
/// `SPIN_FD_REDO_KEPT_HX` next to the start port. FD, fcc, generic shifted
/// 3x3x3 mesh (27 k), P = 1 and pooled P = 3: ports on against off within the
/// strict class, map 1 redone in both arms with fewer H columns on.
///
///   DFT_NC_DIR=/abs/pseudodojo/ex cargo test --release -p dft-engine \
///     --test refscf_fd_kept real_metals_redo -- --ignored --nocapture
#[test]
#[ignore = "needs DFT_NC_DIR (real PseudoDojo UPF); run in release"]
fn real_metals_redo() {
    let Some(dir) = std::env::var_os("DFT_NC_DIR") else {
        println!("GF1 real metals: DFT_NC_DIR unset, skipped");
        return;
    };
    let dir = std::path::PathBuf::from(dir);
    for (element, a_bohr, cutoff, bands) in [("Al", 7.65, 10., 8), ("Au", 7.71, 18., 14)] {
        let pp = dft_upf::from_file(dir.join(format!("{element}.upf"))).unwrap();
        let potentials = [SpeciesPotential {
            label: element,
            potential: &pp,
        }];
        let h = a_bohr / 2.;
        let cell = Cell::from_user_rows_bohr([[0., h, h], [h, 0., h], [h, h, 0.]]).unwrap();
        let atoms = [Atom::new(element, [0.; 3])];
        let points = monkhorst_pack([3; 3], [0.5; 3]).unwrap();
        let options = ScfOptions {
            cutoff_hartree: cutoff,
            functional: Functional::Pbe,
            bands: Some(bands),
            max_iterations: 100,
            verification: Verification::Off,
            gamma_real: GammaReal::Off,
            standard: Some(StandardProfile {
                threshold_ha: STANDARD_HA,
            }),
            ..ScfOptions::strict_reference()
        };
        let policy = Policy::FermiDirac {
            theta_ha: 0.0073,
            band_cap: bands,
        };
        let one = |p: usize, pairs: &[(&str, &str)]| {
            on_rank0(p, |comm| {
                let mut all = vec![("refscf_k_pools", pools_arm(p))];
                all.extend_from_slice(pairs);
                let config = RunConfig::from_pairs(&all).unwrap();
                let _scope = config.enter();
                let r = electronic::run_with(
                    StartInput {
                        cell: &cell,
                        atoms: &atoms,
                        potentials: &potentials,
                        options: &options,
                    },
                    Profile {
                        threshold_ry: 1e-12,
                        ..Profile::default()
                    },
                    &points,
                    policy.clone(),
                    comm,
                )
                .unwrap();
                let h = |i: &electronic::Iteration| {
                    i.calls.iter().map(|(_, _, c)| c.h_columns).sum::<usize>()
                };
                (
                    r.thermodynamics.free_energy_ha,
                    r.converged,
                    r.history.len(),
                    r.history[0].first_iteration_redo,
                    h(&r.history[0]),
                    r.history.iter().map(h).sum::<usize>(),
                )
            })
        };
        for p in [1, 3] {
            let off = one(p, PORTS_OFF);
            let on = one(p, &[]);
            let start_only = one(p, &[("spin_fd_redo_kept_hx", "off")]);
            println!(
                "GF1 real {element} P={p}: F={:.12} |on-off|={:.2e} Ha; maps {} / {}; redo {} / {}; map-1 H off {} start-only {} on {}; H total off {} on {}",
                on.0,
                (on.0 - off.0).abs(),
                off.2,
                on.2,
                off.3,
                on.3,
                off.4,
                start_only.4,
                on.4,
                off.5,
                on.5
            );
            assert!(on.1 && off.1 && start_only.1);
            assert!((on.0 - off.0).abs() <= STRICT_HA && (start_only.0 - off.0).abs() <= STRICT_HA);
            assert!(on.3 && off.3, "{element}: map 1 is redone");
            assert!(
                on.4 < start_only.4 && start_only.4 < off.4,
                "{element}: map-1 H"
            );
        }
    }
}
