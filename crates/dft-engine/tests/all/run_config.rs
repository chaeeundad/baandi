//! RUN-CONFIG (2026-11-18): two SCF runs with opposite switch settings in one
//! process, concurrently, each with its own `RunConfig` passed per call
//! (`ScfProblem::set_run_config`), no lock and no environment.
//!
//! Checks: each concurrent run is bit-identical (energy, density, eigenvalues,
//! iteration count) to the same configuration run alone; the default
//! configuration's run equals a run under the process configuration (no
//! switch set in this test's environment for the migrated switches, so the
//! default science output is unchanged); the run-local counters are
//! independent (kept-coefficient stop estimates only in the `kept` run,
//! transformed ones only in the `fft` run) and every consulted switch is
//! marked in its own configuration only; the two arms agree in the standard
//! class (5e-9 Ha). The construction is also done under the configuration's
//! thread scope (`RunConfig::enter`), so construction-time switches follow it.
#[allow(dead_code)]
#[path = "../common/mod.rs"]
mod common;
use dft_engine::{GammaReal, ScfOptions, ScfProblem, ScfResult, SpeciesPotential};
use dft_geometry::{Atom, Cell};
use dft_numerics::switches::{Id, RunConfig, RunCounter};

const STANDARD_CLASS_HA: f64 = 5e-9;

fn run(pp: &dft_upf::Pseudopotential, rc: Option<RunConfig>) -> ScfResult {
    let cell = Cell::from_user_rows_bohr([[6.5, 0., 0.], [0.6, 5.8, 0.], [0.3, 0.5, 7.]]).unwrap();
    let atoms = vec![
        Atom::new("Xx", [0.13, 0.21, 0.32]),
        Atom::new("Xx", [0.61, 0.57, 0.74]),
    ];
    let build = || {
        ScfProblem::new(
            &cell,
            &atoms,
            &[SpeciesPotential {
                label: "Xx",
                potential: pp,
            }],
            ScfOptions {
                gamma_real: GammaReal::Auto,
                cutoff_hartree: 3.,
                bands: Some(6),
                threads: 2,
                ..ScfOptions::strict_reference()
            },
        )
        .unwrap()
    };
    let problem = match rc {
        Some(rc) => {
            let mut p = {
                let _scope = rc.enter();
                build()
            };
            // The explicit per-call form (construction already took it).
            p.set_run_config(rc);
            assert_eq!(p.run_config(), rc);
            p
        }
        None => build(),
    };
    problem.run().unwrap()
}

fn bits(r: &ScfResult) -> (u64, Vec<u64>, Vec<u64>, usize) {
    (
        r.energy.total_hartree().to_bits(),
        r.density.iter().map(|x| x.to_bits()).collect(),
        r.k_points[0]
            .solution
            .eigenvalues
            .iter()
            .map(|x| x.to_bits())
            .collect(),
        r.history.len(),
    )
}

#[test]
fn concurrent_runs_with_opposite_settings_are_independent() {
    let pp = common::potential(4., 1., true, 401);
    let kept = || RunConfig::from_pairs(&[("ewald_cache", "on")]).unwrap();
    let fft = || RunConfig::from_pairs(&[("ewald_cache", "off")]).unwrap();
    // References, one at a time.
    let (ka, fa) = (kept(), fft());
    let ref_kept = bits(&run(&pp, Some(ka)));
    let ref_fft = bits(&run(&pp, Some(fa)));
    // Concurrent: two threads, opposite settings, fresh configurations.
    let (kb, fb) = (kept(), fft());
    let (out_kept, out_fft) = std::thread::scope(|s| {
        let a = s.spawn(|| bits(&run(&pp, Some(kb))));
        let b = s.spawn(|| bits(&run(&pp, Some(fb))));
        (a.join().unwrap(), b.join().unwrap())
    });
    assert_eq!(out_kept, ref_kept, "kept run changed under concurrency");
    assert_eq!(out_fft, ref_fft, "fft run changed under concurrency");
    // The default science output: the registry default (kept) equals a run
    // under the process configuration when the environment leaves it unset.
    if std::env::var_os("BAANDI_EWALD_CACHE").is_none() {
        assert_eq!(bits(&run(&pp, None)), ref_kept);
    }
    let de = f64::from_bits(ref_kept.0) - f64::from_bits(ref_fft.0);
    assert!(de.abs() <= STANDARD_CLASS_HA, "arms differ by {de:e} Ha");
    // Independent run-local counters.
    for (rc, label) in [(ka, "kept ref"), (kb, "kept concurrent")] {
        let c = rc.counters();
        assert_eq!(c[RunCounter::Problems as usize], 1, "{label}");
        assert!(c[RunCounter::SessionSteps as usize] > 0, "{label}");
    }
    // Same run, same SCF counts. The Ewald slots may differ: the Ewald cache
    // is per construction thread (the second reference run on the main
    // thread hits it, a fresh thread computes), and hits + computed = 1.
    let scf = |rc: RunConfig| {
        let c = rc.counters();
        let ewald = c[RunCounter::EwaldHits as usize] + c[RunCounter::EwaldComputed as usize];
        assert_eq!(ewald, 1);
        [
            RunCounter::StopKept,
            RunCounter::StopFft,
            RunCounter::Problems,
            RunCounter::SessionSteps,
            RunCounter::Verifications,
        ]
        .map(|k| c[k as usize])
    };
    assert_eq!(scf(kb), scf(ka), "same run, same counts");
    assert_eq!(scf(fb), scf(fa), "same run, same counts");
    // Consulted marks belong to the configuration that was read.
    for rc in [ka, fa, kb, fb] {
        assert!(rc.consulted(Id::EWALD_CACHE));
    }
    let untouched = RunConfig::from_input(&[]).unwrap();
    assert!(!untouched.consulted(Id::EWALD_CACHE));
}
