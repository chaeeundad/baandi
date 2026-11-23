//! R-SCF-2 full-path evidence. Slow matrix (real external PP, pure Rust):
//! DFT_SI_UPF=/absolute/Si.pz-vbc.UPF cargo test --release -p dft-engine
//! --test reference_driver -- --include-ignored --nocapture --test-threads=1
#[allow(dead_code)]
mod common;
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
    let ranks: Vec<usize> = std::env::var("RSCF_RANKS")
        .unwrap_or_else(|_| "1,2,3,8".into())
        .split(',')
        .map(|p| p.parse().unwrap())
        .collect();
    let mut baseline: Option<driver::Run> = None;
    let mut rank_ok = true;
    for p in ranks {
        let runs = simulate(p, |comm: &SimRank| {
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
        });
        for run in &runs {
            let result = &run.outcome.result;
            let de = result.total_hartree() - old.total_hartree();
            let attempt = result
                .verification
                .as_ref()
                .unwrap()
                .attempts
                .last()
                .unwrap();
            let last = run.records.last().unwrap();
            println!(
                "RSCF case={name} ranks={p} steps={} energy_ha={:.17e} strict_energy_ha={:.17e} delta_ha={de:.17e} dr2_ry={:.17e} verify_ha={:.17e} residual_bound_ha={:.17e} verified={} charge_error_e={:.3e}",
                run.records.len(),
                result.total_hartree(),
                old.total_hartree(),
                last.dr2_ry,
                attempt.coulomb_error_ha,
                attempt.residual_bound_ha.unwrap(),
                attempt.passed,
                result.charge - 2. * (atoms.len() as f64) * pp.z_valence() / 2.
            );
            assert!(run.outcome.converged, "{name} P={p}: {}", run.exit_reason);
            assert!(attempt.passed, "{name} P={p}: {attempt:?}");
            assert!(de.abs() < 5e-9, "standard accuracy {name} P={p}: {de:e} Ha");
            assert!((result.charge - atoms.len() as f64 * pp.z_valence()).abs() < 1e-10);
            assert!(
                run.records
                    .iter()
                    .all(|r| r.energy.total_hartree().is_finite())
            );
            if let Some(base) = &baseline {
                rank_ok &= run.records.len() == base.records.len();
                let mut max_dr: f64 = 0.;
                let mut max_ethr: f64 = 0.;
                for (r, b) in run.records.iter().zip(&base.records) {
                    if std::env::var_os("RSCF_TRACE").is_some() {
                        println!(
                            "TRACE {name} P={p} k={} dr2={:.16e}/{:.16e} ethr={:.16e}/{:.16e} dav={}/{} H={}/{} cond={:?}/{:?} active={:?}/{:?}",
                            r.iteration,
                            r.dr2_ry,
                            b.dr2_ry,
                            r.ethr_ry,
                            b.ethr_ry,
                            r.davidson_iterations,
                            b.davidson_iterations,
                            r.h_columns,
                            b.h_columns,
                            r.mixing.condition,
                            b.mixing.condition,
                            r.calls.iter().map(|c| &c.active_bands).collect::<Vec<_>>(),
                            b.calls.iter().map(|c| &c.active_bands).collect::<Vec<_>>()
                        );
                    }
                    max_dr = max_dr.max((r.dr2_ry - b.dr2_ry).abs());
                    max_ethr = max_ethr.max((r.ethr_ry - b.ethr_ry).abs());
                    if r.davidson_iterations != b.davidson_iterations
                        && std::env::var_os("RSCF_TRACE").is_some()
                    {
                        println!(
                            "RITZ_DIFF {name} P={p} k={} a={:?} b={:?}",
                            r.iteration, r.calls, b.calls
                        );
                    }
                    // Rank-count roundoff amplified by Broyden history, no
                    // physics-scale threshold/energy tolerance change.
                    let matched = (r.dr2_ry - b.dr2_ry).abs() < 1e-11 + 1e-8 * b.dr2_ry.abs()
                        && (r.ethr_ry - b.ethr_ry).abs() < 1e-13 + 1e-8 * b.ethr_ry.abs()
                        && r.first_iteration_redo == b.first_iteration_redo
                        && r.davidson_iterations == b.davidson_iterations
                        && r.h_columns == b.h_columns;
                    rank_ok &= matched;
                    if !matched {
                        println!(
                            "RSCF_RANK_FAILURE case={name} ranks={p} step={} dr2_delta_ry={:.17e} dav={}/{} h={}/{}",
                            r.iteration,
                            r.dr2_ry - b.dr2_ry,
                            r.davidson_iterations,
                            b.davidson_iterations,
                            r.h_columns,
                            b.h_columns
                        );
                    }
                }
                let dn = result
                    .density
                    .iter()
                    .zip(&base.outcome.result.density)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0., f64::max);
                println!(
                    "RSCF_RANK case={name} ranks={p} max_dr2_ry={max_dr:.17e} max_ethr_ry={max_ethr:.17e} density_max_e_bohr3={dn:.17e}"
                );
                rank_ok &= dn < 1e-9;
            }
        }
        if baseline.is_none() {
            baseline = runs.into_iter().next();
        }
    }
    rank_ok
}
#[test]
#[ignore = "full synthetic SCFs at four rank counts; run release command above"]
fn synthetic_nonlocal() {
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
#[ignore = "external Si PP, 20 full distributed SCFs plus strict oracles; release only"]
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
        "rank trajectory audit failed; all differences are printed above; no tolerance relaxed"
    );
}
