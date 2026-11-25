//! R8a identity and boundary contracts. Synthetic data, never QE evidence.
//! Full matrix: cargo test --release -p dft-engine --test refscf_stepper -- --include-ignored --nocapture
#[allow(dead_code)]
mod common;
use dft_engine::{
    Functional, GammaReal, ScfOptions, SpeciesPotential, Verification,
    collinear::SpinSector,
    refscf::{
        driver,
        electronic::{self, Policy},
        k_driver,
        profile::Profile,
        stepper::{StepKind, StepReport, Stepper},
    },
    refstart::StartInput,
};
use dft_geometry::{
    Atom, Cell,
    kpoints::{KPoint, monkhorst_pack},
};
use dft_mpi::Collective;
use dft_numerics::switches::RunConfig;

// FNV-1a over all little-endian IEEE-754 bytes. Also compare the entire bit
// vectors, so a hash collision cannot hide a difference in the stepwise test.
fn hash(bits: &[u64]) -> u64 {
    bits.iter()
        .flat_map(|v| v.to_le_bytes())
        .fold(0xcbf29ce484222325, |h, b| {
            (h ^ u64::from(b)).wrapping_mul(0x100000001b3)
        })
}
fn scalar_bits(run: &driver::Run) -> (u64, Vec<u64>, Vec<u64>) {
    let r = &run.outcome.result;
    (
        r.energy.total_hartree().to_bits(),
        r.density.iter().map(|x| x.to_bits()).collect(),
        r.k_points
            .iter()
            .flat_map(|k| k.solution.eigenvectors.as_slice())
            .flat_map(|x| [x.re.to_bits(), x.im.to_bits()])
            .collect(),
    )
}
fn electronic_bits(run: &electronic::Run) -> (u64, Vec<u64>, Vec<u64>) {
    (
        run.thermodynamics.free_energy_ha.to_bits(),
        run.density.iter().flatten().map(|x| x.to_bits()).collect(),
        run.blocks
            .iter()
            .flatten()
            .flat_map(|b| b.orbitals.as_slice())
            .flat_map(|x| [x.re.to_bits(), x.im.to_bits()])
            .collect(),
    )
}
fn drive<T>(mut s: Stepper<'_, T>) -> (T, Vec<StepReport>) {
    assert!(!s.is_done());
    assert!(s.result().is_none());
    assert!(s.step().is_err());
    s.initialize().unwrap();
    assert!(s.initialize().is_err());
    assert!(!s.is_done());
    let mut reports = Vec::new();
    while !s.is_done() {
        let r = s.step().unwrap();
        assert!(r.iteration > 0);
        assert!(r.dr2_ry.is_finite());
        assert_eq!(r.estimator_ha.to_bits(), (r.dr2_ry / 2.).to_bits());
        assert!(r.ethr_ry > 0.);
        assert!(r.energy_ha.is_none_or(f64::is_finite));
        reports.push(r);
    }
    assert!(s.result().is_some());
    assert!(s.step().is_err());
    assert!(reports.last().unwrap().safe_point);
    (s.finish().unwrap(), reports)
}
fn fixture(mode: &str, ranks: usize, verify: bool, budget: usize, stepped: bool) {
    let electronic_mode = mode == "fd" || mode == "spin";
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0.2, 5., 0.], [0., 0.1, 5.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    let pp = common::potential_with_atomic_start(
        if electronic_mode { 1. } else { 2. },
        1.,
        true,
        161,
        if mode == "spin" { "SLA PW" } else { "SLA PZ" },
    );
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: 2.5,
        functional: if mode == "spin" {
            Functional::LdaPw
        } else {
            Functional::LdaPz
        },
        bands: Some(if mode == "fd" { 12 } else { 4 }),
        gamma_real: if mode == "gamma" {
            GammaReal::On
        } else {
            GammaReal::Off
        },
        max_iterations: budget,
        verification: if verify {
            Verification::FixedPoint
        } else {
            Verification::Off
        },
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
        let _scope = RunConfig::from_pairs(&[
            ("refscf_verify_energy", "off"),
            ("refscf_verify_gate", "off"),
        ])
        .unwrap()
        .enter();
        let input = || StartInput {
            cell: &cell,
            atoms: &atoms,
            potentials: &potentials,
            options: &options,
        };
        let profile = Profile::default();
        let (a, b, reports, expected_h, iterations, attempts, converged) = if electronic_mode {
            let policy = || {
                if mode == "fd" {
                    Policy::FermiDirac {
                        theta_ha: 0.01,
                        band_cap: 12,
                    }
                } else {
                    Policy::FixedSpin(SpinSector {
                        n_up: 1,
                        n_down: 0,
                        axis_cartesian: [0., 0., 1.],
                        initial_moment_seed: Some(0.3),
                    })
                }
            };
            let a = electronic::run_with(input(), profile, &points, policy(), comm).unwrap();
            let (b, reports) = if stepped {
                drive(electronic::stepper_with(
                    input(),
                    profile,
                    &points,
                    policy(),
                    comm,
                ))
            } else {
                (
                    electronic::run_with(input(), profile, &points, policy(), comm).unwrap(),
                    Vec::new(),
                )
            };
            assert_eq!(
                a.energy.total_hartree().to_bits(),
                b.energy.total_hartree().to_bits()
            );
            assert_eq!(
                a.thermodynamics.entropy_term_ha.to_bits(),
                b.thermodynamics.entropy_term_ha.to_bits()
            );
            let h = a
                .history
                .iter()
                .flat_map(|r| &r.calls)
                .chain(a.verification_calls.iter())
                .map(|(_, _, c)| c.h_columns)
                .sum::<usize>();
            let attempts = a.verification.as_ref().map_or(0, |v| v.attempts.len());
            (
                electronic_bits(&a),
                electronic_bits(&b),
                reports,
                h,
                a.history.len(),
                attempts,
                a.converged,
            )
        } else {
            let a = k_driver::run_with(input(), profile, &points, comm).unwrap();
            let (b, reports) = if stepped {
                drive(k_driver::stepper_with(input(), profile, &points, comm))
            } else {
                (
                    k_driver::run_with(input(), profile, &points, comm).unwrap(),
                    Vec::new(),
                )
            };
            let h = a.records.iter().map(|r| r.h_columns).sum::<usize>()
                + a.verification_calls
                    .iter()
                    .map(|c| c.h_columns)
                    .sum::<usize>();
            let attempts = a
                .outcome
                .result
                .verification
                .as_ref()
                .map_or(0, |v| v.attempts.len());
            (
                scalar_bits(&a),
                scalar_bits(&b),
                reports,
                h,
                a.records.len(),
                attempts,
                a.outcome.converged,
            )
        };
        assert_eq!(a, b, "{mode} ranks={ranks}");
        if budget > 1
            && cfg!(all(target_os = "macos", target_arch = "aarch64"))
            && !dft_numerics::MKL_BACKEND
        {
            assert_eq!(
                (a.0, hash(&a.1), hash(&a.2)),
                parent_fingerprint(mode, ranks),
                "pre-refactor identity: {mode} ranks={ranks}"
            );
        }
        if budget > 1 {
            assert!(converged, "{mode} ranks={ranks} did not converge");
        }
        if stepped {
            assert_eq!(
                reports.iter().map(|r| r.h_columns).sum::<usize>(),
                expected_h
            );
            assert_eq!(
                reports
                    .iter()
                    .filter(|r| r.kind == StepKind::Verification)
                    .count(),
                attempts
            );
            assert_eq!(
                reports
                    .iter()
                    .filter(|r| r.kind == StepKind::Ordinary)
                    .count(),
                iterations
            );
            assert!(reports.iter().filter(|r| r.kind == StepKind::Retry).count() <= 1);
            for pair in reports.windows(2) {
                if pair[1].kind == StepKind::Verification {
                    assert!(!pair[0].safe_point);
                    assert_eq!(pair[0].iteration, pair[1].iteration);
                }
            }
        }
        if comm.rank() == 0 {
            println!(
                "R8A {mode} ranks={ranks} verify={verify} iterations={iterations} attempts={attempts} H={expected_h} energy_bits={:016x} density_hash={:016x} orbitals_hash={:016x}",
                a.0,
                hash(&a.1),
                hash(&a.2)
            );
        }
    });
}
#[test]
fn lifecycle_and_exhaustion() {
    fixture("gamma", 1, false, 1, true);
}
#[test]
#[ignore = "full convergence matrix; run in release (module header)"]
fn gamma_k222s_fd_spin_ranks_1_3_8() {
    for mode in ["gamma", "k222s", "fd", "spin"] {
        for ranks in [1, 3, 8] {
            for verify in [false, true] {
                fixture(mode, ranks, verify, 400, true);
            }
        }
    }
}

// Independently captured by running the run-only version of this fixture on
// git archive 81763a54f574b4a63b625303a91667626f89013b, Rust 1.99.1,
// pure Rust macOS/aarch64, release. The first word is E for fixed occupations
// and F for FD; the matrix additionally compares internal E and entropy bits.
fn parent_fingerprint(mode: &str, ranks: usize) -> (u64, u64, u64) {
    match (mode, ranks) {
        ("gamma", 1) => (0xbfed5e1afb07c70a, 0xee7a0c4277b28991, 0x27546b91b340a3da),
        ("gamma", 3) => (0xbfed5e1afb07c708, 0xa33db0b2cdfc408e, 0x2bd87ac97d1ad9d7),
        ("gamma", 8) => (0xbfed5e1afb07c706, 0x2c9ce88efec2ffb2, 0x1163fd048eb28344),
        ("k222s", 1) => (0xbf91303dc3366000, 0x630b05e7e61ad5f4, 0x4643b41b3f7ca877),
        ("k222s", 3) => (0xbf91303dc3366040, 0x1243d4361ed7e5e7, 0xcdc6b2cce6487340),
        ("k222s", 8) => (0xbf91303dc3365fc0, 0xf66a18b577a5dde6, 0xdc366eba092b4b1b),
        ("fd", 1) => (0xbfcc696ea977493c, 0x4aca2b364c55f301, 0x591a6eb3296b6d58),
        ("fd", 3) => (0xbfcc696ea977493f, 0xee5618e54415fd2a, 0xe25859d1b8ca1c06),
        ("fd", 8) => (0xbfcc696ea9774941, 0xa6fd6e005020aa64, 0xb433d50f760bc79e),
        ("spin", 1) => (0xbfcd9844ebb015e8, 0xe1a419a4c0bb6fdf, 0x1ffb73b91074e5c1),
        ("spin", 3) => (0xbfcd9844ebb015f2, 0x741fa18d86289670, 0x89205267cc664428),
        ("spin", 8) => (0xbfcd9844ebb015f0, 0x4dd5b1b829d5201e, 0xfa4f865022efe2ca),
        _ => unreachable!(),
    }
}

#[test]
fn verification_is_a_separate_step() {
    fixture("gamma", 1, true, 400, true);
}

#[test]
fn facade_partial_run_and_invalid_initialization() {
    use dft_engine::refscf::api;
    let cell = Cell::from_user_rows_bohr([[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: 1.,
        bands: Some(2),
        gamma_real: GammaReal::On,
        max_iterations: 100,
        ..Default::default()
    };
    let input = || StartInput {
        cell: &cell,
        atoms: &atoms,
        potentials: &potentials,
        options: &options,
    };
    let request = || {
        api::Request::new(
            Profile::default(),
            vec![KPoint {
                fractional: [0.; 3],
                weight: 1.,
            }],
        )
    };
    let _scope = RunConfig::from_pairs(&[]).unwrap().enter();
    let a = api::solve(input(), request())
        .unwrap()
        .run
        .into_scalar()
        .unwrap();
    let mut s = api::stepper(input(), request());
    s.initialize().unwrap();
    let mut overlapping = api::stepper(input(), request());
    assert!(overlapping.initialize().is_err());
    s.step().unwrap();
    let b = s.run_to_end().unwrap().run.into_scalar().unwrap();
    assert_eq!(scalar_bits(&a), scalar_bits(&b));
    // The rejected initialization did not poison the new session. Dropping a
    // suspended continuation also restores scopes and releases its thread.
    overlapping.initialize().unwrap();
    drop(overlapping);
    let mut bad = request();
    bad.profile.threshold_ry = f64::NAN;
    let mut s = api::stepper(input(), bad);
    assert!(s.initialize().is_err());
    assert!(s.is_done());
    assert!(s.result().is_none());
    assert!(s.step().is_err());
    assert!(s.finish().is_err());
    assert!(api::stepper(input(), request()).finish().is_err());
}

#[test]
#[ignore = "recovery/reuse and ladder; run in release (module header)"]
fn failed_verification_reuse_and_ladder_boundaries() {
    let cell = Cell::from_user_rows_bohr([[11., 0., 0.], [0.2, 11., 0.], [0., 0.1, 11.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    let pp = common::potential_with_atomic_start(4., 1.2, true, 161, "SLA PZ");
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: 2.,
        bands: Some(3),
        gamma_real: GammaReal::On,
        max_iterations: 100,
        ..Default::default()
    };
    let input = || StartInput {
        cell: &cell,
        atoms: &atoms,
        potentials: &potentials,
        options: &options,
    };
    let profile = Profile {
        history: 4,
        ..Default::default()
    };
    for (reuse, ladder) in [("off", "off"), ("on", "off"), ("off", "on")] {
        let _scope = RunConfig::from_pairs(&[
            ("refscf_verify_energy", "off"),
            ("refscf_verify_gate", "off"),
            ("refscf_verify_reuse", reuse),
            ("refscf_ladder", ladder),
        ])
        .unwrap()
        .enter();
        let a = driver::run(input(), profile).unwrap();
        let (b, reports) = drive(driver::stepper(input(), profile));
        assert_eq!(scalar_bits(&a), scalar_bits(&b));
        let attempts = &a.outcome.result.verification.as_ref().unwrap().attempts;
        assert_eq!(
            reports
                .iter()
                .filter(|r| r.kind == StepKind::Verification)
                .count(),
            attempts.len()
        );
        let coarse_h = a.ladder.as_ref().map_or(0, |l| {
            l.coarse_records.iter().map(|r| r.h_columns).sum::<usize>()
        });
        assert_eq!(
            reports.iter().map(|r| r.h_columns).sum::<usize>(),
            coarse_h
                + a.records.iter().map(|r| r.h_columns).sum::<usize>()
                + a.verification_calls
                    .iter()
                    .map(|c| c.h_columns)
                    .sum::<usize>()
        );
        if ladder == "off" {
            assert!(attempts.iter().any(|a| !a.passed));
            assert_eq!(reports.iter().any(|r| r.reused), reuse == "on");
            assert!(
                reports
                    .iter()
                    .filter(|r| r.reused)
                    .all(|r| r.h_columns == 0)
            );
        } else {
            assert!(reports.iter().any(|r| r.stage == "coarse"));
            assert_eq!(reports.last().unwrap().stage, "fine");
        }
    }
}
