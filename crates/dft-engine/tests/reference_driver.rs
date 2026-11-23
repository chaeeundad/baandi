//! R-SCF-2 full-path evidence. Slow matrix (real external PP, pure Rust):
//! BAANDI_DAVIDSON_STAGES=on DFT_SI_UPF=/absolute/Si.pz-vbc.UPF cargo test --release
//! -p dft-engine --test reference_driver -- --include-ignored --nocapture --test-threads=1
//! The rank-count rule reads per-band Ritz values and thresholds, which are full
//! diagnostics (I26-2): the matrix requires BAANDI_DAVIDSON_STAGES=on.
#[allow(dead_code)]
mod common;
#[path = "common/rank_rule.rs"]
mod rank_rule;
use dft_engine::{
    refscf::{driver, profile::Profile},
    refstart::StartInput,
    *,
};
use dft_geometry::{Atom, Cell};
use dft_mpi::{SimRank, simulate};

fn check(
    name: &str,
    cell: &Cell,
    atoms: &[Atom],
    pp: &dft_upf::Pseudopotential,
    options: ScfOptions,
) -> bool {
    let species = [SpeciesPotential {
        label: "Xx",
        potential: pp,
    }];
    let strict = ScfOptions {
        standard: None,
        verification: Verification::Off,
        eigen_tolerance: 1e-11,
        energy_tolerance: 1e-12,
        density_tolerance: 1e-9,
        mixing: Mixing::Pulay {
            beta: 0.3,
            history: 8,
        },
        ..options.clone()
    };
    let old = ScfProblem::new(cell, atoms, &species, strict)
        .unwrap()
        .run()
        .unwrap_or_else(|e| panic!("strict {name}: {e}"));
    let mut ranks: Vec<usize> = std::env::var("RSCF_RANKS")
        .unwrap_or_else(|_| "1,2,3,8".into())
        .split(',')
        .map(|p| p.parse().unwrap())
        .collect();
    ranks.push(1);
    ranks.sort_unstable();
    ranks.dedup();
    assert!(ranks.iter().all(|&p| p > 0));
    let mut baseline: Option<driver::Run> = None;
    let mut rank_ok = true;
    for p in ranks {
        let run_layout = || {
            simulate(p, |comm: &SimRank| {
                driver::run_with(
                    StartInput {
                        cell,
                        atoms,
                        potentials: &species,
                        options: &options,
                    },
                    Profile::default(),
                    comm,
                )
                .unwrap()
            })
        };
        let runs = run_layout();
        assert!(
            runs.iter().all(|r| r.full_diagnostics),
            "rank-count rule needs full diagnostics: set BAANDI_DAVIDSON_STAGES=on"
        );
        let repeated = run_layout();
        for (run, repeat) in runs.iter().zip(&repeated) {
            assert_eq!(
                rank_rule::fingerprint(run),
                rank_rule::fingerprint(repeat),
                "fixed-layout repeatability {name} P={p}"
            );
            assert!(
                rank_rule::delivered(
                    run,
                    old.total_hartree(),
                    atoms.len() as f64 * pp.z_valence()
                ),
                "delivered state {name} P={p}: {}",
                run.exit_reason
            );
            let base = baseline.as_ref().unwrap_or(&runs[0]);
            let audit = rank_rule::accept(
                base,
                run,
                old.total_hartree(),
                atoms.len() as f64 * pp.z_valence(),
            );
            if let Err(error) = &audit {
                println!("RSCF_RANK_FAILURE case={name} ranks={p}: {error}");
                rank_ok = false;
            }
        }
        let run = &runs[0];
        let result = &run.outcome.result;
        let base = baseline.as_ref().unwrap_or(run);
        let audit = rank_rule::accept(
            base,
            run,
            old.total_hartree(),
            atoms.len() as f64 * pp.z_valence(),
        );
        let first = match &audit {
            Ok(branches) if branches.is_empty() => "none".to_string(),
            Ok(branches) => {
                for b in branches {
                    println!(
                        "RSCF_BRANCH case={name} ranks={p} map={} verification={} call={} correction={} band={} changes_ha={:.17e}/{:.17e} thresholds_ha={:.17e}/{:.17e} margins_ha={:.17e}/{:.17e} rounding_ha={:.17e}",
                        b.map,
                        b.verification,
                        b.call,
                        b.correction,
                        b.band,
                        b.changes[0],
                        b.changes[1],
                        b.thresholds[0],
                        b.thresholds[1],
                        b.margins[0],
                        b.margins[1],
                        b.rounding
                    );
                }
                let b = &branches[0];
                format!(
                    "{}:{}:{}:{}",
                    b.map,
                    if b.verification { "verify" } else { "scf" },
                    b.call,
                    b.correction
                )
            }
            Err(e) => format!("FAIL({e})"),
        };
        let verification = result.verification.as_ref().unwrap();
        let attempt = verification.attempts.last().unwrap();
        println!(
            "RSCF_MATRIX case={name} ranks={p} first={first} steps={} energy_ha={:.17e} strict_energy_ha={:.17e} delta_ha={:.17e} density_max_e_bohr3={:.17e} verify_ha={:.17e} residual_estimate_ha={:.17e} charge_error_e={:.17e} delivered={} repeat_bits=true delivered_wall_s={:.6} verification_wall_s={:.6} attempts={} h_columns={}",
            run.records.len(),
            result.total_hartree(),
            old.total_hartree(),
            result.total_hartree() - old.total_hartree(),
            rank_rule::density_delta(base, run),
            attempt.coulomb_error_ha,
            attempt.residual_bound_ha.unwrap(),
            result.charge - atoms.len() as f64 * pp.z_valence(),
            rank_rule::delivered(
                run,
                old.total_hartree(),
                atoms.len() as f64 * pp.z_valence()
            ),
            result.total_time.as_secs_f64(),
            verification.wall().as_secs_f64(),
            verification.attempts.len(),
            run.records.iter().map(|r| r.h_columns).sum::<usize>()
                + run
                    .verification_calls
                    .iter()
                    .map(|r| r.h_columns)
                    .sum::<usize>()
        );
        let mut rejected = repeated.into_iter().next().unwrap();
        rejected
            .outcome
            .result
            .verification
            .as_mut()
            .unwrap()
            .attempts
            .last_mut()
            .unwrap()
            .passed = false;
        assert_eq!(
            rank_rule::accept(
                base,
                &rejected,
                old.total_hartree(),
                atoms.len() as f64 * pp.z_valence()
            )
            .unwrap_err(),
            "delivered-state checks failed"
        );
        if baseline.is_none() {
            baseline = runs.into_iter().next();
        }
    }
    rank_ok
}
#[test]
#[ignore = "full synthetic SCFs at four rank counts; run release command above"]
fn synthetic_nonlocal_fixed_layout_repeatability() {
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    let cell = Cell::from_user_rows_bohr([[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    assert!(check(
        "synthetic",
        &cell,
        &atoms,
        &pp,
        ScfOptions {
            cutoff_hartree: 2.,
            bands: Some(3),
            start_density: StartDensity::Atomic,
            gamma_real: GammaReal::On,
            max_iterations: 200,
            ..Default::default()
        }
    ));
}
#[test]
#[ignore = "external Si PP, 40 distributed SCFs plus strict oracles; release only"]
fn silicon_generality_matrix() {
    let pp = dft_upf::parse_upf(
        &std::fs::read_to_string(std::env::var("DFT_SI_UPF").expect("DFT_SI_UPF required"))
            .unwrap(),
    )
    .unwrap();
    let base = [
        [0., 0., 0.],
        [0., 0.5, 0.5],
        [0.5, 0., 0.5],
        [0.5, 0.5, 0.],
        [0.25, 0.25, 0.25],
        [0.25, 0.75, 0.75],
        [0.75, 0.25, 0.75],
        [0.75, 0.75, 0.25],
    ];
    let mut rank_ok = true;
    for (name, rx, ry, skew, stretch) in [
        ("si8", 1, 1, false, 1.),
        ("si16", 2, 1, false, 1.),
        ("si32", 2, 2, false, 1.),
        ("skew", 1, 1, true, 1.),
        ("elongated", 3, 1, false, 1.),
    ] {
        if std::env::var("RSCF_MATRIX_CASE")
            .ok()
            .is_some_and(|s| s != name)
        {
            continue;
        }
        let a = 10.261_212_856_717_933;
        let cell = Cell::from_user_rows_bohr([
            [a * rx as f64, 0., 0.],
            [if skew { 0.7 } else { 0. }, a * ry as f64, 0.],
            [0., if skew { 0.3 } else { 0. }, a * stretch],
        ])
        .unwrap();
        let mut atoms = Vec::new();
        for x in 0..rx {
            for y in 0..ry {
                for p in base {
                    let i = atoms.len();
                    atoms.push(Atom::new(
                        "Xx",
                        [
                            (p[0] + x as f64) / rx as f64 + 1e-4 * ((i * 7 + 1) as f64).sin(),
                            (p[1] + y as f64) / ry as f64 + 1e-4 * ((i * 11 + 3) as f64).sin(),
                            p[2] + 1e-4 * ((i * 13 + 5) as f64).sin(),
                        ],
                    ));
                }
            }
        }
        let (cell, atoms) = if skew {
            (
                Cell::from_user_rows_bohr([
                    [0.7, a / 2., a / 2.],
                    [a / 2., 0., a / 2.],
                    [a / 2., a / 2., 0.3],
                ])
                .unwrap(),
                vec![
                    Atom::new(
                        "Xx",
                        [
                            0.9992102689475646,
                            0.9997106070535452,
                            0.0006016594584669142,
                        ],
                    ),
                    Atom::new(
                        "Xx",
                        [0.2499873540117459, 0.24999285393389137, 0.2506104968350948],
                    ),
                ],
            )
        } else {
            (cell, atoms)
        };
        rank_ok &= check(
            name,
            &cell,
            &atoms,
            &pp,
            ScfOptions {
                cutoff_hartree: 15.,
                bands: Some(2 * atoms.len() + 4),
                start_density: StartDensity::Atomic,
                gamma_real: GammaReal::On,
                max_iterations: 200,
                ..Default::default()
            },
        );
    }
    assert!(
        rank_ok,
        "rank acceptance failed; see first divergence and delivered-state evidence above"
    );
}
