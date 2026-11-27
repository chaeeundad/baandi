//! GK10B: first-map work of the reference k path. Prints, per fixture, the
//! Davidson iterations and H columns of map 1 split into the coarse first
//! solve (ethr 1e-2 Ry), the first-iteration redo and the in-map retries,
//! against the ordinary maps. Synthetic data, never QE evidence.
//! cargo test --release -p dft-engine --test refscf_k_first_map -- --ignored --nocapture
#[allow(dead_code)]
mod common;
use dft_engine::{
    GammaReal, ScfOptions, SpeciesPotential, StandardProfile, StartDensity, Verification,
    refscf::{k_driver, profile::Profile},
    refstart::StartInput,
};
use dft_geometry::{
    Atom, Cell,
    kpoints::{KPoint, monkhorst_pack},
};
use dft_numerics::switches::RunConfig;

pub struct FirstMap {
    pub energy: f64,
    pub maps: usize,
    /// (iterations, H columns) of map 1: coarse calls, redo calls.
    pub coarse: (usize, usize),
    pub redo: (usize, usize),
    pub coarse_calls: usize,
    pub redo_calls: usize,
    /// Mean (iterations, H columns) of maps 2.. .
    pub ordinary: (f64, f64),
    pub total: (usize, usize),
}

#[allow(clippy::too_many_arguments)]
pub fn run(
    comm: &dft_mpi::SimRank,
    cell: &Cell,
    atoms: &[Atom],
    pp: &dft_upf::Pseudopotential,
    points: &[KPoint],
    bands: usize,
    cutoff: f64,
    threshold_ry: f64,
    start_density: StartDensity,
    pairs: &[(&str, &str)],
) -> FirstMap {
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: cutoff,
        bands: Some(bands),
        max_iterations: 100,
        verification: Verification::Off,
        gamma_real: GammaReal::Off,
        standard: Some(StandardProfile { threshold_ha: 5e-9 }),
        start_density,
        ..ScfOptions::strict_reference()
    };
    let config = RunConfig::from_pairs(pairs).unwrap();
    let _scope = config.enter();
    let run = k_driver::run_with(
        StartInput {
            cell,
            atoms,
            potentials: &potentials,
            options: &options,
        },
        Profile {
            threshold_ry,
            ..Profile::default()
        },
        points,
        comm,
    )
    .unwrap();
    assert!(run.outcome.converged, "{}", run.exit_reason);
    let first = &run.records[0];
    if std::env::var("GK10B_CALLS").is_ok() {
        for (m, r) in run.records.iter().take(3).enumerate() {
            for c in &r.calls {
                println!(
                    "  map {} ethr {:.2e} first {} it {} H {} restarts {} unconverged {} reduced orders {:?}",
                    m + 1,
                    c.ethr_ry,
                    c.first_step,
                    c.iterations,
                    c.h_columns,
                    c.restarts,
                    c.unconverged,
                    c.reduced_orders
                );
            }
            println!("  map {} dr2 {:.3e}", m + 1, r.dr2_ry);
        }
    }
    let requested = first.ethr_requested_ry;
    let (mut coarse, mut redo, mut cc, mut rc) = ((0, 0), (0, 0), 0, 0);
    for c in &first.calls {
        if c.ethr_ry == requested {
            coarse.0 += c.iterations;
            coarse.1 += c.h_columns;
            cc += 1;
        } else {
            redo.0 += c.iterations;
            redo.1 += c.h_columns;
            rc += 1;
        }
    }
    let rest = &run.records[1..];
    let n = rest.len().max(1) as f64;
    FirstMap {
        energy: run.outcome.result.total_hartree(),
        maps: run.records.len(),
        coarse,
        redo,
        coarse_calls: cc,
        redo_calls: rc,
        ordinary: (
            rest.iter().map(|r| r.davidson_iterations).sum::<usize>() as f64 / n,
            rest.iter().map(|r| r.h_columns).sum::<usize>() as f64 / n,
        ),
        total: (
            run.records.iter().map(|r| r.davidson_iterations).sum(),
            run.records.iter().map(|r| r.h_columns).sum(),
        ),
    }
}

pub fn print(label: &str, f: &FirstMap) {
    println!(
        "GK10B {label}: E={:.15e} Ha maps {}; map1 coarse {} calls it/H {}/{}, redo {} calls it/H {}/{}; ordinary mean it/H {:.0}/{:.0}; total it/H {}/{}",
        f.energy,
        f.maps,
        f.coarse_calls,
        f.coarse.0,
        f.coarse.1,
        f.redo_calls,
        f.redo.0,
        f.redo.1,
        f.ordinary.0,
        f.ordinary.1,
        f.total.0,
        f.total.1
    );
}

pub fn cells() -> [(&'static str, Cell); 3] {
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

pub fn gamma_in_list() -> Vec<KPoint> {
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

/// The k fixture matrix of GK10 (2 atoms, 6 bands, 3 Ha) and a 12-atom,
/// 22-band synthetic of B12/GeTe shape, at the standard threshold.
pub fn matrix(pairs: &[(&str, &str)], tag: &str, threshold_ry: f64) -> Vec<(String, FirstMap)> {
    let mut out = Vec::new();
    for (start_name, start) in [
        // The reference start always builds the atomic density; the option
        // does not change it, so one variant is enough.
        ("atomic-start", StartDensity::Uniform),
    ] {
        let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
        let atoms = [
            Atom::new("Xx", [0.13, 0.21, 0.32]),
            Atom::new("Xx", [0.61, 0.55, 0.71]),
        ];
        for (name, cell) in cells() {
            for (set, points) in [
                ("gamma-in-list", gamma_in_list()),
                ("mp222-shifted", monkhorst_pack([2; 3], [0.5; 3]).unwrap()),
            ] {
                let f = dft_mpi::simulate(1, |comm| {
                    run(
                        comm,
                        &cell,
                        &atoms,
                        &pp,
                        &points,
                        6,
                        3.,
                        threshold_ry,
                        start,
                        pairs,
                    )
                })
                .swap_remove(0);
                let label = format!("{tag} {start_name} {name} {set}");
                print(&label, &f);
                out.push((label, f));
            }
        }
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
        let points: Vec<_> = (0..2)
            .map(|i| KPoint {
                fractional: [0.13 + 0.1 * i as f64, -0.21, 0.07 + 0.05 * i as f64],
                weight: 0.5,
            })
            .collect();
        let f = dft_mpi::simulate(1, |comm| {
            run(
                comm,
                &cell,
                &atoms,
                &pp,
                &points,
                22,
                12.,
                threshold_ry,
                start,
                pairs,
            )
        })
        .swap_remove(0);
        let label = format!("{tag} {start_name} b12-shape 12 atoms 22 bands 2 k");
        print(&label, &f);
        out.push((label, f));
    }
    out
}

#[test]
#[ignore = "first-map work report; run in release"]
fn first_map_report() {
    // Arms: default (GK10B kept start and redo HX; the pre-GK10B "before" arm
    // was retired in CL1, its comparison is in the attempts ledger) and empty
    // (default plus the research empty-band rule). Strict arm 1e-12 Ry
    // (class 1e-10 Ha) and standard arm 1e-8 Ry (class 5e-9 Ha, judged
    // against the strict default run of the same fixture).
    let empty = [("k_first_empty_rule", "on")];
    let strict_after = matrix(&[], "strict-default", 1e-12);
    let strict_empty = matrix(&empty, "strict-empty", 1e-12);
    let after = matrix(&[], "default", 1e-8);
    let with_empty = matrix(&empty, "empty", 1e-8);
    for (i, (label, a)) in after.iter().enumerate() {
        let x = &with_empty[i].1;
        let (sa, sx) = (&strict_after[i].1, &strict_empty[i].1);
        let strict_x = (sx.energy - sa.energy).abs();
        let [ea, ex] = [a, x].map(|r| (r.energy - sa.energy).abs());
        let map1 = |r: &FirstMap| (r.coarse.0 + r.redo.0, r.coarse.1 + r.redo.1);
        println!(
            "GK10B delta {label}: strict |empty-default| {strict_x:.2e} Ha; standard vs strict default {ea:.2e} empty {ex:.2e} Ha; map1 it/H default {:?} empty {:?}; maps {}/{}; total H {}/{}; strict maps {}/{}, total H {}/{}",
            map1(a),
            map1(x),
            a.maps,
            x.maps,
            a.total.1,
            x.total.1,
            sa.maps,
            sx.maps,
            sa.total.1,
            sx.total.1
        );
    }
}
