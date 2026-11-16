//! ADR 0043: fixed-point verification after a standard stop. Synthetic fixture.
use super::*;
use dft_geometry::kpoints::KPoint;

fn problem(
    potential: &Pseudopotential,
    standard: Option<f64>,
    verification: Verification,
) -> Result<ScfProblem, ScfError> {
    let cell = Cell::from_user_rows_bohr([[4., 0., 0.], [0.3, 4.2, 0.], [0., 0.2, 4.4]]).unwrap();
    ScfProblem::with_kpoints(
        &cell,
        &[Atom::new("Xx", [0.13, 0.21, 0.32])],
        &[SpeciesPotential {
            label: "Xx",
            potential,
        }],
        ScfOptions {
            standard: standard.map(|threshold_ha| StandardProfile { threshold_ha }),
            verification,
            cutoff_hartree: 2.5,
            bands: Some(3),
            ..ScfOptions::strict_reference()
        },
        &[KPoint {
            fractional: [0.; 3],
            weight: 1.0,
        }],
        OccupationPolicy::Insulator {
            minimum_gap_hartree: 1e-7,
        },
    )
}

fn forces_of(
    potential: &Pseudopotential,
    r: &ScfResult,
    options: &ScfOptions,
) -> Result<Vec<[f64; 3]>, ScfError> {
    let cell = Cell::from_user_rows_bohr([[4., 0., 0.], [0.3, 4.2, 0.], [0., 0.2, 4.4]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    let species = [SpeciesPotential {
        label: "Xx",
        potential,
    }];
    crate::forces::forces(&crate::forces::ForceState {
        cell: &cell,
        atoms: &atoms,
        potentials: &species,
        options,
        electronic: crate::forces::ElectronicState::Fixed(r),
        density_floor: None,
        symmetry: &[],
    })
    .map(|f| f.values)
}

#[test]
fn verification_requires_standard_profile() {
    let pp = super::pool_fixture::potential(2.0, 1.0, true, 401);
    assert!(matches!(
        problem(&pp, None, Verification::FixedPoint),
        Err(ScfError::InvalidInput(_))
    ));
    assert_eq!(
        ScfOptions::strict_reference().verification,
        Verification::Off
    );
}

/// V4: verification is an observation. A passing first attempt delivers the
/// converged step's state bit for bit; forces then work, and match strict.
#[test]
fn passed_verification_keeps_state_and_enables_forces() {
    let pp = super::pool_fixture::potential(2.0, 1.0, true, 401);
    let off = problem(&pp, Some(5e-9), Verification::Off)
        .unwrap()
        .run()
        .unwrap();
    let verified_problem = problem(&pp, Some(5e-9), Verification::FixedPoint).unwrap();
    let options = verified_problem.options.clone();
    let on = verified_problem.run().unwrap();
    assert!(off.verification.is_none());
    let report = on.verification.as_ref().expect("verification report");
    assert_eq!(report.attempts.len(), 1, "{:?}", report.attempts);
    let attempt = &report.attempts[0];
    println!("verification attempt {attempt:?}");
    assert!(attempt.passed && report.passed());
    assert!(attempt.coulomb_error_ha < attempt.threshold_ha);
    assert!(attempt.eigen_tolerance_ha <= on.history.last().unwrap().requested_eigen_tolerance);
    assert!(attempt.davidson_iterations > 0 && attempt.operator_applications > 0);
    assert!(
        !on.certified,
        "a verified standard state is not strict-certified"
    );
    assert_eq!(on.history.len(), off.history.len());
    assert_eq!(on.total_hartree().to_bits(), off.total_hartree().to_bits());
    assert!(
        on.density
            .iter()
            .zip(&off.density)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
    let cert = on.history.last().unwrap().fixed_point_certificate.as_ref();
    assert_eq!(
        cert.map(|c| c.operator_applications),
        Some(attempt.operator_applications)
    );

    // Unverified standard state: rejected exactly as before.
    let mut off_options = options.clone();
    off_options.verification = Verification::Off;
    let message = forces_of(&pp, &off, &off_options).unwrap_err().to_string();
    assert!(message.contains("strict convergence"), "{message}");
    // Verified: accepted, and close to the strict forces of the same geometry.
    let standard_forces = forces_of(&pp, &on, &options).unwrap();
    let strict_problem = problem(&pp, None, Verification::Off).unwrap();
    let strict_options = strict_problem.options.clone();
    let strict = strict_problem.run().unwrap();
    let strict_forces = forces_of(&pp, &strict, &strict_options).unwrap();
    let worst = standard_forces
        .iter()
        .flatten()
        .zip(strict_forces.iter().flatten())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    println!("standard {standard_forces:?} strict {strict_forces:?} max |dF| {worst:e} Ha/Bohr");
    assert!(standard_forces.iter().flatten().all(|f| f.is_finite()));
    assert!(worst < 1e-4, "max force difference {worst:e} Ha/Bohr");
}

/// V5: a failed verification continues SCF and verifies again.
#[test]
fn failed_verification_continues_scf_and_verifies_again() {
    let pp = super::pool_fixture::potential(2.0, 1.0, true, 401);
    let first = problem(&pp, Some(5e-9), Verification::FixedPoint)
        .unwrap()
        .run()
        .unwrap();
    let mut forced = problem(&pp, Some(5e-9), Verification::FixedPoint).unwrap();
    forced.force_verification_failures = 1;
    let outcome = forced.run_outcome().unwrap();
    assert!(outcome.converged);
    let r = outcome.result;
    let report = r.verification.as_ref().unwrap();
    assert_eq!(report.attempts.len(), 2, "{:?}", report.attempts);
    assert!(!report.attempts[0].passed && report.attempts[1].passed);
    assert!(report.attempts[1].iteration > report.attempts[0].iteration);
    assert!(r.history.len() > first.history.len());
    assert_eq!(
        r.history
            .iter()
            .filter(|h| h.fixed_point_certificate.is_some())
            .count(),
        2
    );
    assert!(report.wall() >= report.attempts[0].wall);
    // Budget exhaustion with every attempt failing: not converged, as today.
    let mut never = problem(&pp, Some(5e-9), Verification::FixedPoint).unwrap();
    never.force_verification_failures = usize::MAX;
    never.options.max_iterations = first.history.len() + 2;
    let outcome = never.run_outcome().unwrap();
    assert!(!outcome.converged);
    assert!(!outcome.result.verification.as_ref().unwrap().passed());
}

/// VERIFY-REUSE (`BAANDI_VERIFY_REUSE=warm|feed`): a failed verification's solve
/// is reused (orbitals as the next warm start; with feed also its density pair in
/// the mixer). The run still ends only on a passed verification, at the same
/// fixed point as the default (standard class, gate 5e-9 Ha), and off is the
/// default path bit for bit.
#[test]
fn failed_verification_reuse_keeps_the_fixed_point() {
    use super::session::VerifyReuse;
    let pp = super::pool_fixture::potential(2.0, 1.0, true, 401);
    let run = |mode: Option<VerifyReuse>, failures: usize| {
        let mut p = problem(&pp, Some(5e-9), Verification::FixedPoint).unwrap();
        p.force_verification_failures = failures;
        p.force_verify_reuse = mode;
        p.run_outcome().unwrap()
    };
    let default = run(None, 2);
    let off = run(Some(VerifyReuse::Off), 2);
    assert_eq!(
        default.result.energy.total_hartree().to_bits(),
        off.result.energy.total_hartree().to_bits()
    );
    assert!(
        default
            .result
            .density
            .iter()
            .zip(&off.result.density)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
    let reference = default.result.energy.total_hartree();
    for mode in [VerifyReuse::Warm, VerifyReuse::Feed] {
        let o = run(Some(mode), 2);
        assert!(o.converged, "{mode:?}");
        let report = o.result.verification.as_ref().unwrap();
        assert!(report.passed(), "{mode:?}");
        assert_eq!(report.attempts.len(), 3, "{mode:?} {:?}", report.attempts);
        let de = (o.result.energy.total_hartree() - reference).abs();
        println!(
            "{mode:?}: |dE| {de:e} Ha, history {} vs {}",
            o.result.history.len(),
            default.result.history.len()
        );
        assert!(de < 5e-9, "{mode:?} |dE| {de:e} Ha");
    }
    // No failure: the reuse never triggers, bit-identical to the default.
    let plain = run(None, 0);
    let feed = run(Some(VerifyReuse::Feed), 0);
    assert_eq!(
        plain.result.energy.total_hartree().to_bits(),
        feed.result.energy.total_hartree().to_bits()
    );
}

/// ONEPATH-4 (user decision 2026-11-16): the engine default is the standard
/// path with the verification; every other field equals the strict reference,
/// which keeps the former default bit for bit.
#[test]
fn engine_default_is_the_verified_standard_path() {
    let d = ScfOptions::default();
    let s = ScfOptions::strict_reference();
    assert_eq!(d.standard, Some(StandardProfile { threshold_ha: 5e-9 }));
    assert_eq!(d.verification, Verification::FixedPoint);
    assert!(matches!(
        d.mixing,
        Mixing::Pulay {
            beta: 0.5,
            history: 6
        }
    ));
    assert!(!d.is_strict_reference() && s.is_strict_reference());
    assert_eq!((s.standard, s.verification), (None, Verification::Off));
    assert!(matches!(
        s.mixing,
        Mixing::Pulay {
            beta: 0.3,
            history: 6
        }
    ));
    let aligned = ScfOptions {
        standard: None,
        verification: Verification::Off,
        mixing: s.mixing.clone(),
        ..d.clone()
    };
    assert!(crate::checkpoint::options_identical(&aligned, &s));
    // A default problem converges with a passed verification, uncertified.
    let pp = super::pool_fixture::potential(2.0, 1.0, true, 401);
    let cell = Cell::from_user_rows_bohr([[4., 0., 0.], [0.3, 4.2, 0.], [0., 0.2, 4.4]]).unwrap();
    let r = ScfProblem::new(
        &cell,
        &[Atom::new("Xx", [0.13, 0.21, 0.32])],
        &[SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }],
        ScfOptions {
            cutoff_hartree: 2.5,
            bands: Some(3),
            ..ScfOptions::default()
        },
    )
    .unwrap()
    .run()
    .unwrap();
    assert!(r.verification.as_ref().is_some_and(|v| v.passed()));
    assert!(!r.certified && r.standard.is_some());
}

fn identity(
    problem: &ScfProblem,
) -> (
    crate::checkpoint::ModelIdentity,
    crate::checkpoint::ExecutionIdentity,
) {
    (
        crate::checkpoint::ModelIdentity {
            cell_bohr: problem.cell.lattice_columns_bohr(),
            atoms_fractional: vec![[0.13, 0.21, 0.32]],
            species: vec!["Xx".into()],
            pseudopotential_sha256: vec!["0".repeat(64)],
            cutoff_hartree: problem.options.cutoff_hartree,
            kpoints: vec![[0., 0., 0., 1.]],
            electrons: problem.electrons,
            xc: problem.options.functional.name().into(),
            spin: "unpolarized".into(),
            environment: "none".into(),
        },
        crate::checkpoint::ExecutionIdentity {
            build_sha256: "1".repeat(64),
            platform: "synthetic-same-test-build".into(),
            threads: problem.options.threads,
        },
    )
}

fn finish(mut s: ScfSession) -> ScfOutcome {
    while s.snapshot().status == StepStatus::Ready {
        s.step().unwrap();
    }
    s.into_result().unwrap()
}

/// ONEPATH-4: checkpoint v1 carries the verification option and the failed
/// attempts. An exact continuation taken after a failed verification resumes
/// bit for bit (energy, density, history, attempts without wall time) against
/// the live session continued after its checkpoint; tampered attempt records
/// and a different verification option reject.
#[test]
fn exact_checkpoint_after_failed_verification_resumes_bit_for_bit() {
    use crate::checkpoint::{CheckpointKind, Rejection};
    let pp = super::pool_fixture::potential(2.0, 1.0, true, 401);
    let fresh = || problem(&pp, Some(5e-9), Verification::FixedPoint).unwrap();
    let mut forced = fresh();
    forced.force_verification_failures = 1;
    let (model, execution) = identity(&forced);
    let mut live = ScfSession::new(
        forced,
        CompletionMode::RetainLast,
        crate::events::ExecutionConfig::default(),
    )
    .unwrap();
    live.initialize().unwrap();
    while live
        .history()
        .last()
        .is_none_or(|h| h.fixed_point_certificate.is_none())
    {
        assert_eq!(live.step().unwrap().status, StepStatus::Ready);
    }
    let cp = live
        .checkpoint(
            CheckpointKind::ExactContinuation,
            model.clone(),
            execution.clone(),
        )
        .unwrap();
    assert_eq!(cp.verification, Verification::FixedPoint);
    let saved = &cp.continuation.as_ref().unwrap().verification_attempts;
    assert_eq!(saved.len(), 1);
    assert!(!saved[0].passed);
    let restored = ScfSession::restore_exact(fresh(), &cp, &model, &execution).unwrap();
    let a = finish(live);
    let b = finish(restored);
    assert!(a.converged && b.converged);
    assert_eq!(
        a.result.total_hartree().to_bits(),
        b.result.total_hartree().to_bits()
    );
    assert!(
        a.result
            .density
            .iter()
            .zip(&b.result.density)
            .all(|(x, y)| x.to_bits() == y.to_bits())
    );
    assert_eq!(a.result.history.len(), b.result.history.len());
    let strip = |o: &ScfOutcome| {
        o.result
            .verification
            .as_ref()
            .unwrap()
            .attempts
            .iter()
            .map(|a| {
                format!(
                    "{:?}",
                    VerificationAttempt {
                        wall: std::time::Duration::ZERO,
                        ..a.clone()
                    }
                )
            })
            .collect::<Vec<_>>()
    };
    assert_eq!(strip(&a), strip(&b));
    assert_eq!(strip(&a).len(), 2);
    // Tampering: attempts dropped, attempt marked passed, mode mismatch.
    let mut dropped = cp.clone();
    dropped
        .continuation
        .as_mut()
        .unwrap()
        .verification_attempts
        .clear();
    assert_eq!(
        ScfSession::restore_exact(fresh(), &dropped, &model, &execution).err(),
        Some(Rejection::History)
    );
    let mut passed = cp.clone();
    passed.continuation.as_mut().unwrap().verification_attempts[0].passed = true;
    assert_eq!(
        ScfSession::restore_exact(fresh(), &passed, &model, &execution).err(),
        Some(Rejection::History)
    );
    let mut mode = cp.clone();
    mode.verification = Verification::Off;
    assert_eq!(
        ScfSession::restore_exact(fresh(), &mode, &model, &execution).err(),
        Some(Rejection::Options)
    );
    // A target problem with the verification off: options mismatch.
    let off = problem(&pp, Some(5e-9), Verification::Off).unwrap();
    assert_eq!(
        ScfSession::restore_exact(off, &cp, &model, &execution).err(),
        Some(Rejection::Options)
    );
}
