//! M2-GRID (`BAANDI_M2_GRID=on|off`, default off, rounding class): the
//! distributed full-grid stages of the M2 path (Hartree FFT pair, XC and
//! effective potential on the planes of a slab FFT, the energy audit's
//! output-density transform and grid sums, the stop-rule estimate from kept
//! owned coefficients) on the production standard SCF (fixed-point
//! verification), cubic and skewed cells, simulated ranks.
//!
//! Checks: the standard accuracy class (5e-9 Ha) against the replicated M2 path
//! and the threads-only run, a passed verification, every rank bit-equal, a
//! fixed rank count bit-reproducible, the default (switch unset) bit-identical
//! to the explicit off arm, and the counters prove the distributed stages ran
//! (potential, audit and stop rule in every map). Own test binary with a lock:
//! it reads process-wide counters.
#[allow(dead_code)]
mod common;
use dft_engine::{GammaReal, ScfOptions, ScfProblem, ScfResult, SpeciesPotential, parallel};
use dft_geometry::{Atom, Cell};
use dft_mpi::{Collective, SimRank, simulate};
use std::sync::{Arc, Mutex};

static LOCK: Mutex<()> = Mutex::new(());
const STANDARD_CLASS_HA: f64 = 5e-9;

fn cell(kind: usize) -> Cell {
    Cell::from_user_rows_bohr(if kind == 0 {
        [[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]]
    } else {
        [[6.5, 0., 0.], [0.6, 5.8, 0.], [0.3, 0.5, 7.]]
    })
    .unwrap()
}

/// Two atoms of four valence electrons, six bands, production standard path.
fn run(kind: usize, pp: &dft_upf::Pseudopotential) -> ScfResult {
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
    .run()
    .unwrap()
}

/// `grid`: None leaves the switch at its default.
fn ranks_run(
    ranks: usize,
    kind: usize,
    pp: &dft_upf::Pseudopotential,
    grid: Option<bool>,
) -> Vec<ScfResult> {
    simulate(ranks, |comm: &SimRank| {
        parallel::set_m2_override(Some(true));
        dft_engine::set_stop_metric_override(Some(true));
        dft_engine::grid_dist::set_override(grid);
        let shared: Arc<dyn Collective + Send + Sync> = Arc::new(comm.clone());
        parallel::set_communicator(Some(shared));
        let r = run(kind, pp);
        parallel::set_communicator(None);
        dft_engine::grid_dist::set_override(None);
        dft_engine::set_stop_metric_override(None);
        parallel::set_m2_override(None);
        r
    })
}

fn bits(r: &ScfResult) -> (u64, Vec<u64>, Vec<u64>, usize) {
    (
        r.total_hartree().to_bits(),
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
fn distributed_grid_stages_keep_the_standard_class() {
    let _lock = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let pp = common::potential(4., 1., true, 401);
    for kind in [0, 1] {
        let threads = run(kind, &pp);
        let default = ranks_run(2, kind, &pp, None);
        let off = ranks_run(2, kind, &pp, Some(false));
        assert_eq!(
            bits(&default[0]),
            bits(&off[0]),
            "default path is the off arm"
        );
        let c0 = dft_engine::grid_dist::counters();
        let on = ranks_run(2, kind, &pp, Some(true));
        let c1 = dft_engine::grid_dist::counters();
        let again = ranks_run(2, kind, &pp, Some(true));
        let on3 = ranks_run(3, kind, &pp, Some(true));
        // Every map ran the distributed potential and audit; the stop rule
        // took the kept owned coefficients; nothing declined.
        let maps = on[0].history.len() as u64;
        assert!(
            c1[0] - c0[0] >= 2 * maps,
            "potential maps {:?} -> {:?}",
            c0,
            c1
        );
        assert!(c1[1] - c0[1] >= 2 * maps, "audits");
        assert!(c1[2] > c0[2], "stop-rule estimates");
        assert_eq!(c1[4], c0[4], "no declined map");
        for (p, arm) in [(2, &on), (3, &on3)] {
            for r in arm.iter() {
                assert!(r.verification.as_ref().is_some_and(|v| v.passed()));
                assert_eq!(bits(r), bits(&arm[0]), "ranks agree, P={p}");
            }
            let d_off = (arm[0].total_hartree() - off[0].total_hartree()).abs();
            let d_thr = (arm[0].total_hartree() - threads.total_hartree()).abs();
            println!(
                "cell {kind} P={p}: |dE| vs replicated M2 {d_off:.3e} Ha, vs threads-only {d_thr:.3e} Ha, SCF {} (off {})",
                arm[0].history.len(),
                off[0].history.len()
            );
            assert!(d_off <= STANDARD_CLASS_HA && d_thr <= STANDARD_CLASS_HA);
        }
        assert_eq!(
            bits(&on[0]),
            bits(&again[0]),
            "fixed rank count reproducible"
        );
    }
}
