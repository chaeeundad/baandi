//! BREAKTHROUGH-3 (`BAANDI_SCF_ADVANCE`, design docs/design/bounded-subspace-scf.md):
//! the production standard SCF (fixed-point verification) with intermediate
//! maps run as bounded subspace updates. Own test binary: it reads the
//! process-wide advance counters.
#[allow(dead_code)]
mod common;
use dft_engine::{GammaReal, ScfOptions, ScfProblem, ScfResult, SpeciesPotential};
use dft_geometry::{Atom, Cell};

fn cell(kind: usize) -> Cell {
    Cell::from_user_rows_bohr(if kind == 0 {
        [[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]]
    } else {
        [[6.5, 0., 0.], [0.6, 5.8, 0.], [0.3, 0.5, 7.]]
    })
    .unwrap()
}

/// Two atoms of four valence electrons (four occupied bands, six bands), the
/// production standard path (verification on), serial real Gamma path.
fn problem(kind: usize, pp: &dft_upf::Pseudopotential) -> ScfProblem {
    let atoms = vec![
        Atom::new("Xx", [0.13, 0.21, 0.32]),
        Atom::new("Xx", [0.61, 0.57, 0.74]),
    ];
    ScfProblem::new(
        &cell(kind),
        &atoms,
        &[SpeciesPotential {
            label: "Xx",
            potential: pp,
        }],
        ScfOptions {
            gamma_real: GammaReal::Auto,
            cutoff_hartree: 3.,
            bands: Some(6),
            threads: 1,
            ..ScfOptions::default()
        },
    )
    .unwrap()
}

fn run(kind: usize, pp: &dft_upf::Pseudopotential, advance: Option<usize>) -> ScfResult {
    dft_engine::advance::set_override(Some(advance.map_or((false, 2), |b| (true, b))));
    let r = problem(kind, pp).run();
    dft_engine::advance::set_override(None);
    r.unwrap()
}

fn work(r: &ScfResult) -> (usize, usize, usize, usize) {
    (
        r.history.len(),
        r.history.iter().map(|h| h.davidson_iterations).sum(),
        r.history.iter().map(|h| h.operator_applications).sum(),
        r.verification.as_ref().map_or(0, |v| v.attempts.len()),
    )
}

/// With the switch on (blocks 1 and 2, cubic and skewed cells) the run
/// converges with a passed fixed-point verification, the total energy stays
/// in the standard class (5e-9 Ha, pre-registered) of the switch-off run,
/// intermediate maps carry advance records while the first map does not, the
/// stop step is an exact map or a completed candidate, every accepted
/// advanced map obeys the inner rule (B <= 0.1 E), and the run is
/// reproducible. With the switch off no record exists. On a classic-Davidson
/// process every request declines and the run is the exact path.
#[test]
fn scf_advance_keeps_the_standard_class_and_records_every_advanced_map() {
    use dft_engine::{ADVANCE_INNER_RATIO, CompletionTrigger};
    let pp = common::potential(4., 1., true, 401);
    let lean = dft_numerics::switches::choice(dft_numerics::switches::Id::DAVIDSON) == "lean";
    for kind in 0..2 {
        let off = run(kind, &pp, None);
        assert!(off.verification.as_ref().is_some_and(|v| v.passed()));
        assert!(off.history.iter().all(|h| h.advance.is_none()));
        for blocks in [1usize, 2] {
            let before = dft_engine::advance_counters();
            let on = run(kind, &pp, Some(blocks));
            let again = run(kind, &pp, Some(blocks));
            let after = dft_engine::advance_counters();
            let de = (on.total_hartree() - off.total_hartree()).abs();
            println!(
                "kind {kind} blocks {blocks}: off (SCF, Davidson, H, verification) {:?} on {:?}, dE {de:e} Ha, advanced {} exact {} completions {} declined {}",
                work(&off),
                work(&on),
                after.maps_advanced - before.maps_advanced,
                after.maps_exact - before.maps_exact,
                after.completions() - before.completions(),
                after.declined - before.declined,
            );
            assert!(on.verification.as_ref().is_some_and(|v| v.passed()));
            assert!(de <= 5e-9, "kind {kind} blocks {blocks}: dE {de:e} Ha");
            assert_eq!(
                on.total_hartree().to_bits(),
                again.total_hartree().to_bits(),
                "not reproducible"
            );
            assert_eq!(work(&on), work(&again));
            assert!(on.history[0].advance.is_none(), "the first map is exact");
            // The stop rule fires on exact maps only: the stop step is an
            // ordinary map or an advanced map that was completed (by any
            // trigger; a candidate completion when the advanced map itself
            // was already below the threshold).
            let last = on.history.last().unwrap();
            assert!(
                last.advance.as_ref().is_none_or(|a| a.completion.is_some()),
                "the stop step is an exact map or a completed map"
            );
            assert!(
                on.history
                    .iter()
                    .filter_map(|h| h.advance.as_ref())
                    .filter(|a| a.scf_error_ha < 5e-9)
                    .all(|a| a.completion == Some(CompletionTrigger::Candidate)),
                "an advanced map below the threshold is always completed as a candidate"
            );
            if !lean {
                assert!(on.history.iter().all(|h| h.advance.is_none()));
                assert!(after.declined > before.declined);
                continue;
            }
            let advanced: Vec<_> = on
                .history
                .iter()
                .filter_map(|h| h.advance.as_ref())
                .collect();
            assert!(!advanced.is_empty(), "kind {kind}: no advanced map");
            // The first advanced map uses the starting budget; later ones the
            // adapted budget (one block up after an inner-noise completion,
            // one down after a very small ratio), always inside [1, 8].
            assert_eq!(advanced[0].budget_blocks, blocks);
            for a in &advanced {
                assert!((1..=8).contains(&a.budget_blocks));
                assert!(a.blocks_used >= 1 && a.blocks_used <= a.budget_blocks);
                assert!(a.h_columns >= 1);
                assert!(a.inner_estimate_ha.is_finite() && a.scf_error_ha.is_finite());
                if a.completion.is_none() {
                    assert!(
                        a.inner_estimate_ha <= ADVANCE_INNER_RATIO * a.scf_error_ha,
                        "accepted map violates the inner rule: B {:e} E {:e}",
                        a.inner_estimate_ha,
                        a.scf_error_ha
                    );
                    assert_eq!(a.completion_h_columns, 0);
                } else {
                    assert!(a.completion_h_columns >= 1);
                }
            }
            assert_eq!(
                (after.maps_advanced - before.maps_advanced) as usize,
                2 * advanced.len()
            );
            assert_eq!(after.declined, before.declined);
        }
    }
}
