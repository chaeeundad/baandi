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
            ..Default::default()
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
    assert_eq!(ScfOptions::default().verification, Verification::Off);
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
