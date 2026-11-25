//! Reference v2 state identity, independent uninterrupted trajectory oracle.
//! Full matrix: cargo test --release -p dft-engine --test refscf_checkpoint -- --include-ignored --nocapture
#[allow(dead_code)]
mod common;
use dft_engine::{
    GammaReal, ScfOptions, SpeciesPotential, Verification,
    refscf::{checkpoint::Identity, driver::Run, k_driver, profile::Profile, stepper::Stepper},
    refstart::StartInput,
};
use dft_geometry::{
    Atom, Cell,
    kpoints::{KPoint, monkhorst_pack},
};
use dft_mpi::Collective;

fn identity(ranks: usize) -> Identity {
    Identity {
        input_sha256: "1".repeat(64),
        upf_sha256: vec!["2".repeat(64)],
        binary_sha256: "3".repeat(64),
        platform: "synthetic-test".into(),
        threads: 1,
        ranks,
    }
}
fn bits(run: &Run) -> (u64, Vec<u64>, Vec<u64>) {
    let r = &run.outcome.result;
    (
        r.energy.total_hartree().to_bits(),
        r.density.iter().map(|v| v.to_bits()).collect(),
        r.k_points
            .iter()
            .flat_map(|k| k.solution.eigenvectors.as_slice())
            .flat_map(|v| [v.re.to_bits(), v.im.to_bits()])
            .collect(),
    )
}
fn hash(v: &[u64]) -> String {
    format!(
        "{:016x}",
        v.iter()
            .flat_map(|v| v.to_le_bytes())
            .fold(0xcbf29ce484222325u64, |h, b| (h ^ u64::from(b))
                .wrapping_mul(0x100000001b3))
    )
}
fn boundary(s: &mut Stepper<'_, Run>, iteration: usize) {
    s.initialize().unwrap();
    if iteration == 0 {
        return;
    }
    loop {
        let r = s.step().unwrap();
        assert!(!s.is_done(), "terminal before requested checkpoint");
        if r.iteration >= iteration && r.safe_point {
            break;
        }
    }
}
fn fixture(mode: &str, ranks: usize, budget: usize, points_to_save: &[usize], reuse: bool) {
    let skew = mode == "skew-two-species";
    let cell = Cell::from_user_rows_bohr(if skew {
        [[6., 0., 0.], [0.6, 5.5, 0.], [0.1, 0.3, 6.5]]
    } else {
        [[5., 0., 0.], [0., 5., 0.], [0., 0., 5.]]
    })
    .unwrap();
    let mut atoms = vec![Atom::new("Xx", [0.13, 0.21, 0.32])];
    if skew {
        atoms.push(Atom::new("Yy", [0.63, 0.71, 0.72]));
    }
    let pp = common::potential_with_atomic_start(2., 1., true, 161, "SLA PZ");
    let pp2 = common::potential_with_atomic_start(2., 1.2, true, 161, "SLA PZ");
    let potentials = [
        SpeciesPotential {
            label: "Xx",
            potential: &pp,
        },
        SpeciesPotential {
            label: "Yy",
            potential: &pp2,
        },
    ];
    let options = ScfOptions {
        cutoff_hartree: 2.5,
        bands: Some(if skew { 6 } else { 4 }),
        gamma_real: if mode == "k222s" {
            GammaReal::Off
        } else {
            GammaReal::On
        },
        max_iterations: budget,
        verification: Verification::FixedPoint,
        ..Default::default()
    };
    let points = if mode == "k222s" {
        monkhorst_pack([2; 3], [0.5; 3]).unwrap()
    } else {
        vec![KPoint {
            fractional: [0.; 3],
            weight: 1.,
        }]
    };
    dft_mpi::simulate(ranks, |comm| {
        let _scope = dft_numerics::switches::RunConfig::from_pairs(&[])
            .unwrap()
            .enter();
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
        let uninterrupted = build().run_to_end().unwrap();
        let expected = bits(&uninterrupted);
        if budget > 10 {
            assert!(
                uninterrupted.outcome.converged,
                "{mode}: {}",
                uninterrupted.exit_reason
            );
        }
        for &at in points_to_save {
            let mut s = build();
            boundary(&mut s, at);
            let cp = s.checkpoint(identity(ranks), comm.rank()).unwrap();
            assert_eq!(
                cp.state,
                s.checkpoint(identity(ranks), comm.rank()).unwrap().state,
                "capture changed state"
            );
            let captured = s.run_to_end().unwrap();
            assert_eq!(bits(&captured), expected, "capture changed trajectory");
            if budget == 6 && at == 0 {
                for damage in 0..3 {
                    let mut malformed = cp.clone();
                    match damage {
                        0 => {
                            malformed.state.pop();
                        }
                        1 => malformed.state.extend_from_slice(&[0; 8]),
                        _ => malformed.iteration += 1,
                    }
                    let mut rejected = build();
                    let error = rejected
                        .restore(&malformed, &identity(ranks), comm.rank())
                        .unwrap_err()
                        .to_string();
                    assert!(error.contains("reference checkpoint v2"), "{error}");
                    assert!(rejected.is_done(), "failed import retained live state");
                }
            }
            let mut restored = build();
            restored
                .restore(&cp, &identity(ranks), comm.rank())
                .unwrap();
            assert_eq!(
                cp.state,
                restored
                    .checkpoint(identity(ranks), comm.rank())
                    .unwrap()
                    .state,
                "restore changed saved state"
            );
            let resumed = restored.run_to_end().unwrap();
            assert_eq!(
                bits(&resumed),
                expected,
                "restart at {at}, {mode}, ranks {ranks}"
            );
            assert_eq!(
                resumed.outcome.result.history.len(),
                uninterrupted.outcome.result.history.len()
            );
            if comm.rank() == 0 {
                let row = format!(
                    "{{\"fixture\":\"{mode}\",\"ranks\":{ranks},\"budget\":{budget},\"converged\":{},\"safe_point\":{},\"reuse\":{reuse},\"energy_ha\":{:.17e},\"energy_bits\":\"{:016x}\",\"density_fnv1a\":\"{}\",\"orbitals_fnv1a\":\"{}\",\"identity\":true}}",
                    uninterrupted.outcome.converged,
                    cp.iteration,
                    f64::from_bits(expected.0),
                    expected.0,
                    hash(&expected.1),
                    hash(&expected.2)
                );
                println!("R8B_MEMORY {row}");
                if let Ok(path) = std::env::var("BAANDI_CHECKPOINT_EVIDENCE") {
                    use std::io::Write;
                    writeln!(
                        std::fs::OpenOptions::new()
                            .create(true)
                            .append(true)
                            .open(path)
                            .unwrap(),
                        "{row}"
                    )
                    .unwrap();
                }
            }
        }
    });
}
#[test]
fn short_safe_point_identity() {
    fixture("gamma", 1, 6, &[0, 1, 3], false);
}
#[test]
#[ignore = "full convergence matrix; run in release (module header)"]
fn gamma_k222s_skew_ranks_1_3_8() {
    for mode in ["gamma", "k222s", "skew-two-species"] {
        for ranks in [1, 3, 8] {
            fixture(mode, ranks, 400, &[0, 1, 3], false);
        }
    }
}

#[test]
fn rejects_identity_and_corruption() {
    let cp = dft_engine::refscf::checkpoint::Checkpoint {
        schema_version: 2,
        identity: identity(1),
        rank: 0,
        iteration: 0,
        state: vec![0; 8],
    };
    let mut wrong = identity(1);
    wrong.input_sha256 = "4".repeat(64);
    assert!(
        cp.validate(&wrong, 0)
            .unwrap_err()
            .to_string()
            .contains("identity mismatch")
    );
    assert!(cp.validate(&identity(1), 1).is_err());
    let mut bad = cp.clone();
    bad.schema_version = 3;
    assert!(bad.validate(&identity(1), 0).is_err());
}

#[test]
#[ignore = "additional spin/FD two-map continuation; run in release"]
fn electronic_fixed_spin_and_fd_restart() {
    use dft_engine::refscf::electronic::{self, Policy};
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0.2, 5., 0.], [0., 0.1, 5.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    for spin in [false, true] {
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
            cutoff_hartree: 2.5,
            bands: Some(if spin { 4 } else { 12 }),
            gamma_real: GammaReal::Off,
            functional: if spin {
                dft_engine::Functional::LdaPw
            } else {
                dft_engine::Functional::LdaPz
            },
            max_iterations: 2,
            verification: Verification::Off,
            ..Default::default()
        };
        let points = [KPoint {
            fractional: [0.; 3],
            weight: 1.,
        }];
        dft_mpi::simulate(3, |comm| {
            let _scope = dft_numerics::switches::RunConfig::from_pairs(&[])
                .unwrap()
                .enter();
            let build = || {
                electronic::stepper_with(
                    StartInput {
                        cell: &cell,
                        atoms: &atoms,
                        potentials: &potentials,
                        options: &options,
                    },
                    Profile::default(),
                    &points,
                    if spin {
                        Policy::FixedSpin(dft_engine::collinear::SpinSector {
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
                    },
                    comm,
                )
            };
            let bits = |r: &electronic::Run| {
                (
                    r.energy.total_hartree().to_bits(),
                    r.thermodynamics.free_energy_ha.to_bits(),
                    r.density
                        .iter()
                        .flatten()
                        .map(|v| v.to_bits())
                        .collect::<Vec<_>>(),
                    r.blocks
                        .iter()
                        .flatten()
                        .flat_map(|b| b.orbitals.as_slice())
                        .flat_map(|v| [v.re.to_bits(), v.im.to_bits()])
                        .collect::<Vec<_>>(),
                )
            };
            let expected = build().run_to_end().unwrap();
            assert!(!expected.converged);
            for at in [0, 1] {
                let mut s = build();
                s.initialize().unwrap();
                if at > 0 {
                    loop {
                        let r = s.step().unwrap();
                        assert!(!s.is_done());
                        if r.iteration >= at && r.safe_point {
                            break;
                        }
                    }
                }
                let cp = s.checkpoint(identity(3), comm.rank()).unwrap();
                drop(s);
                let mut s = build();
                s.restore(&cp, &identity(3), comm.rank()).unwrap();
                assert_eq!(bits(&s.run_to_end().unwrap()), bits(&expected));
                if comm.rank() == 0 {
                    println!("R8B_ELECTRONIC spin={spin} iteration={at} identity=true");
                }
            }
        });
    }
}
