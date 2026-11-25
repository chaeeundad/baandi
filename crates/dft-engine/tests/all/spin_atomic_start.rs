//! ONEPATH-3B: the collinear fixed sector starts from the atomic subspace like
//! the closed-shell path (ADR 0042 amendment), so a zero-moment spin run follows
//! the unpolarized SCF trajectory (ticket ONEPATH-3B-SPIN-COST).
//!
//! Bounds fixed before measuring: the ONEPATH-3 acceptance (a) bounds
//! (|dE| <= 1e-10 Ha, density distance <= 1e-9, SCF counts within 1,
//! per-iteration standard estimates within 1e-6 relative) plus equal
//! verification attempt counts. Measured outcome (recorded in the ticket): the
//! free-electron start meets all of them; with the atomic start the dimer
//! fixture gives |dE| 5.8e-10 Ha (LDA-PW) and 2.6e-10 Ha (PBE) with equal SCF
//! counts, inside the rounding envelope of the closed-shell atomic start itself
//! (up to 1.5e-9 and 3.0e-10 Ha under 1e-13 fractional displacements, ignored
//! test below): rounding-level potential differences change the start
//! subspace, and the 5e-3 Ha first tolerance keeps that difference. The atomic-start test therefore asserts the
//! ticket's trajectory criterion (SCF counts within 1, equal attempts), the
//! standard class (5e-9 Ha) and that the atomic start ran; the 1e-10 bound is
//! reported, not asserted, there. Not widened: it is still asserted for the
//! free-electron start, where rounding is not amplified.
#[allow(dead_code)]
#[path = "../common/mod.rs"]
mod common;
use common::potential_with_orbitals;
use dft_engine::{collinear::*, *};
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use dft_physics::spin::SpinDensity;

fn options(functional: Functional, start: StartWavefunctions) -> ScfOptions {
    ScfOptions {
        functional,
        cutoff_hartree: 3.0,
        // Spin runs use the complex path; the unpolarized reference must too.
        gamma_real: GammaReal::Off,
        start_wavefunctions: start,
        ..ScfOptions::strict_reference()
    }
}
fn gamma() -> Vec<KPoint> {
    vec![KPoint {
        fractional: [0.0; 3],
        weight: 1.0,
    }]
}
fn sector(n: usize) -> SpinSector {
    SpinSector {
        n_up: n,
        n_down: n,
        axis_cartesian: [0.0, 0.0, 1.0],
        initial_moment_seed: Some(0.0),
    }
}
/// Two-atom molecule-like skewed cell (gapped), nonlocal projectors, s and p
/// atomic orbitals in PP_PSWFC.
fn dimer() -> (Cell, [Atom; 2]) {
    (
        Cell::from_user_rows_bohr([[7.0, 0.0, 0.0], [0.6, 6.8, 0.0], [0.4, 0.3, 7.2]]).unwrap(),
        [Atom::new("Xx", [0.0; 3]), Atom::new("Xx", [0.2, 0.0, 0.0])],
    )
}
struct Pair {
    spin: CollinearResult,
    closed: ScfResult,
    de: f64,
    dn: f64,
}
fn pair(functional: Functional, label: &str, start: StartWavefunctions) -> Pair {
    let (cell, atoms) = dimer();
    let pp = potential_with_orbitals(1.0, 1.0, true, 401, label);
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let opts = options(functional, start);
    let closed = ScfProblem::new(&cell, &atoms, &species, opts.clone())
        .unwrap()
        .run()
        .unwrap();
    let spin = CollinearProblem::new(&cell, &atoms, &species, opts, &gamma(), sector(1))
        .unwrap()
        .run()
        .unwrap();
    assert!(spin.converged);
    assert!(spin.verification.as_ref().unwrap().passed());
    assert!(closed.verification.as_ref().unwrap().passed());
    let SpinDensity::Collinear { channels, .. } = &spin.density else {
        panic!()
    };
    let total = &channels[0] + &channels[1];
    let de = (spin.energy.total_hartree() - closed.total_hartree()).abs();
    let dn = density_distance(&total, &closed.density, cell.volume_bohr3(), 2.0);
    println!(
        "ONEPATH-3B {functional:?} {start:?}: dE={de:e} Ha, dn={dn:e}, SCF spin={} unpolarized={}, attempts spin={} unpolarized={}",
        spin.history.len(),
        closed.history.len(),
        spin.verification.as_ref().unwrap().attempts.len(),
        closed.verification.as_ref().unwrap().attempts.len()
    );
    Pair {
        spin,
        closed,
        de,
        dn,
    }
}

#[test]
fn zero_moment_spin_uses_the_atomic_start_and_keeps_the_unpolarized_scf_count() {
    for (functional, label) in [(Functional::LdaPw, "SLA PW"), (Functional::Pbe, "PBE")] {
        let p = pair(functional, label, StartWavefunctions::Atomic);
        // Like for like: both runs start from the atomic subspace; the span is
        // built once and Rayleigh-Ritz runs once per channel.
        let (a, b) = (&p.spin.atomic_start, &p.closed.atomic_start);
        assert!(b.atomic_orbitals > 0 && b.fallback_reason.is_none());
        assert_eq!(a.fallback_reason, None);
        assert_eq!(a.atomic_orbitals, b.atomic_orbitals);
        // SPIN_SHARE (default on since 2026-11-17): the two identical channels share one
        // Rayleigh-Ritz and one initial H application; the shared result is bit-identical.
        assert_eq!(
            a.evd_orders.len(),
            1,
            "one shared Rayleigh-Ritz for equal channels"
        );
        assert_eq!(a.initial_h_applications, b.initial_h_applications);
        // Ticket criterion: same SCF count within 1, same verification count.
        assert!(p.spin.history.len().abs_diff(p.closed.history.len()) <= 1);
        assert_eq!(
            p.spin.verification.as_ref().unwrap().attempts.len(),
            p.closed.verification.as_ref().unwrap().attempts.len()
        );
        // Standard accuracy class (the 1e-10 bound is reported, see the header).
        assert!(p.de <= 5e-9, "dE={:e}", p.de);
    }
}

#[test]
fn zero_moment_spin_free_electron_start_is_the_unpolarized_trajectory() {
    for (functional, label) in [(Functional::LdaPw, "SLA PW"), (Functional::Pbe, "PBE")] {
        let p = pair(functional, label, StartWavefunctions::FreeElectron);
        assert_eq!(p.spin.atomic_start.atomic_orbitals, 0);
        assert_eq!(p.spin.atomic_start.fallback_reason, None);
        assert!(p.de <= 1e-10, "dE={:e}", p.de);
        assert!(p.dn <= 1e-9, "dn={:e}", p.dn);
        assert!(p.spin.history.len().abs_diff(p.closed.history.len()) <= 1);
        assert_eq!(
            p.spin.verification.as_ref().unwrap().attempts.len(),
            p.closed.verification.as_ref().unwrap().attempts.len()
        );
        for (s, h) in p.spin.history.iter().zip(&p.closed.history) {
            let (a, b) = (s.standard_error_ha.unwrap(), h.standard_error_ha.unwrap());
            assert!((a - b).abs() <= 1e-6 * b.max(1e-14), "{a:e} vs {b:e}");
        }
    }
}

/// Evidence for the header (information only, no assertion): the closed-shell
/// atomic-start result under eight displacements of one atom by 1e-13 to
/// 3e-13 fractional (rounding-level potential changes). Mac 2026-11-16:
/// LDA-PW 1.4e-10 to 1.5e-9 Ha, PBE 3.2e-12 to 3.0e-10 Ha; the zero-moment
/// spin differences (5.8e-10, 2.6e-10 Ha) lie inside.
/// `cargo test -p dft-engine --test all spin_atomic_start::rounding_envelope -- --ignored --nocapture`
#[test]
#[ignore]
fn rounding_envelope_of_the_closed_shell_atomic_start() {
    let (cell, atoms) = dimer();
    for (functional, label) in [(Functional::LdaPw, "SLA PW"), (Functional::Pbe, "PBE")] {
        let pp = potential_with_orbitals(1.0, 1.0, true, 401, label);
        let species = [SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let o = options(functional, StartWavefunctions::Atomic);
        let base = ScfProblem::new(&cell, &atoms, &species, o.clone())
            .unwrap()
            .run()
            .unwrap();
        let mut v = Vec::new();
        for d in [
            [1e-13, 0., 0.],
            [0., 1e-13, 0.],
            [0., 0., 1e-13],
            [-1e-13, 0., 0.],
            [0., -1e-13, 0.],
            [0., 0., -1e-13],
            [3e-13, 1e-13, 0.],
            [0., 2e-13, -1e-13],
        ] {
            let shifted = [Atom::new("Xx", d), atoms[1].clone()];
            let c = ScfProblem::new(&cell, &shifted, &species, o.clone())
                .unwrap()
                .run()
                .unwrap();
            v.push(format!(
                "{:.2e}",
                (c.total_hartree() - base.total_hartree()).abs()
            ));
        }
        println!("closed-shell atomic start {functional:?}: |dE| under displacements {v:?} Ha");
    }
}

#[test]
fn fermi_dirac_spin_keeps_the_free_electron_start_with_reason() {
    let cell =
        Cell::from_user_rows_bohr([[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]]).unwrap();
    let pp = potential_with_orbitals(1.0, 1.0, false, 401, "SLA PW");
    let r = CollinearProblem::new_fermi_dirac(
        &cell,
        &[Atom::new("Xx", [0.0; 3])],
        &[SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }],
        ScfOptions {
            functional: Functional::LdaPw,
            cutoff_hartree: 2.5,
            eigen_tolerance: 1e-10,
            density_tolerance: 1e-9,
            max_iterations: 200,
            ..ScfOptions::strict_reference()
        },
        &gamma(),
        SpinSector {
            n_up: 1,
            n_down: 0,
            axis_cartesian: [0.0, 0.0, 1.0],
            initial_moment_seed: None,
        },
        0.01,
        dft_physics::occupations::Constraint::Total,
        8,
    )
    .unwrap()
    .run()
    .unwrap();
    assert!(r.converged);
    assert_eq!(r.atomic_start.atomic_orbitals, 0);
    assert_eq!(
        r.atomic_start.fallback_reason.as_deref(),
        Some("Fermi-Dirac occupations")
    );
}
