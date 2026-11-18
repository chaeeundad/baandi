//! M2-OWNED (`BAANDI_M2_LOCAL_FIELDS`, `BAANDI_M2_LOCAL_DENSITY`, default off):
//! the owned effective potential (and density) of the M2 path on the
//! production standard SCF (fixed-point verification), simulated ranks.
//!
//! LOCAL_FIELDS: bit-identical to `BAANDI_M2_GRID=on` alone (the operator reads
//! the same plane values; only the full-grid gather is skipped), on a cubic, a
//! skewed and an elongated cell whose slab plan and distributed operator use
//! crossed plane axes (the owned potential is then redistributed by one
//! alltoall), P = 2 and 3; the counters prove owned maps, epoch checks, fewer
//! potential gathers and no consumer gather. Own test binary with a lock: it
//! reads process-wide counters.
#[allow(dead_code)]
mod common;
use dft_engine::{GammaReal, ScfOptions, ScfProblem, ScfResult, SpeciesPotential, parallel};
use dft_geometry::{Atom, Cell};
use dft_mpi::{Collective, SimRank, simulate};
use std::sync::{Arc, Mutex};

static LOCK: Mutex<()> = Mutex::new(());

fn cell(kind: usize) -> Cell {
    Cell::from_user_rows_bohr(match kind {
        0 => [[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]],
        1 => [[6.5, 0., 0.], [0.6, 5.8, 0.], [0.3, 0.5, 7.]],
        // Elongated along b with a square (x, y) grid: the wavefunction sphere
        // has more sticks along y (operator planes along y), the full grid
        // ties (slab planes along x).
        _ => [[6., 0., 0.], [0., 9., 0.], [0., 0., 6.]],
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
            grid_override: (kind == 2).then_some([18, 18, 12]),
            ..ScfOptions::default()
        },
    )
    .unwrap()
    .run()
    .unwrap()
}

/// (grid, mix, local fields, local density) on `ranks` simulated ranks.
fn ranks_run(
    ranks: usize,
    kind: usize,
    pp: &dft_upf::Pseudopotential,
    arm: (bool, bool, bool, bool),
) -> Vec<ScfResult> {
    simulate(ranks, |comm: &SimRank| {
        parallel::set_m2_override(Some(true));
        dft_engine::set_stop_metric_override(Some(true));
        dft_engine::grid_dist::set_override(Some(arm.0));
        dft_engine::mix_dist::set_override(Some(arm.1));
        dft_engine::grid_dist::set_owned_override(Some((arm.2, arm.3)));
        let shared: Arc<dyn Collective + Send + Sync> = Arc::new(comm.clone());
        parallel::set_communicator(Some(shared));
        let r = run(kind, pp);
        parallel::set_communicator(None);
        dft_engine::grid_dist::set_owned_override(None);
        dft_engine::mix_dist::set_override(None);
        dft_engine::grid_dist::set_override(None);
        dft_engine::set_stop_metric_override(None);
        parallel::set_m2_override(None);
        r
    })
}

type Bits = (u64, Vec<u64>, Vec<u64>, usize, Vec<u64>);
fn bits(r: &ScfResult) -> Bits {
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
        r.history
            .iter()
            .map(|h| h.operator_applications as u64)
            .collect(),
    )
}

fn owned() -> std::collections::BTreeMap<String, usize> {
    dft_engine::grid_dist::owned_stats().0.into_iter().collect()
}
fn delta(
    a: &std::collections::BTreeMap<String, usize>,
    b: &std::collections::BTreeMap<String, usize>,
    key: &str,
) -> usize {
    b[&format!("owned.{key}")] - a[&format!("owned.{key}")]
}

#[test]
fn owned_fields_are_bit_identical_to_the_gathered_grid_path() {
    let _lock = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let pp = common::potential(4., 1., true, 401);
    for kind in [0, 1, 2] {
        for ranks in [2, 3] {
            let c0 = owned();
            let gathered = ranks_run(ranks, kind, &pp, (true, false, false, false));
            let c1 = owned();
            let fields = ranks_run(ranks, kind, &pp, (true, false, true, false));
            let c2 = owned();
            for r in &fields {
                assert!(r.verification.as_ref().is_some_and(|v| v.passed()));
                assert_eq!(bits(r), bits(&fields[0]), "ranks agree");
            }
            assert_eq!(
                bits(&fields[0]),
                bits(&gathered[0]),
                "cell {kind} P={ranks}: owned fields change no bit"
            );
            let maps = delta(&c1, &c2, "fields_owned_maps");
            let gathers_off = delta(&c0, &c1, "potential_gathers");
            let gathers_on = delta(&c1, &c2, "potential_gathers");
            println!(
                "cell {kind} P={ranks}: owned maps {maps}, potential gathers {gathers_off} -> {gathers_on}, epoch checks {}, redistributions {}, consumer gathers {}, declined {}",
                delta(&c1, &c2, "epoch_checks"),
                delta(&c1, &c2, "redistributions"),
                delta(&c1, &c2, "consumer_potential_gathers"),
                delta(&c1, &c2, "fields_declined_maps"),
            );
            assert!(maps > 0, "owned maps ran");
            assert_eq!(
                gathers_off - gathers_on,
                maps,
                "one gather saved per owned map"
            );
            assert_eq!(
                delta(&c1, &c2, "epoch_checks"),
                maps,
                "every owned map checked"
            );
            assert_eq!(delta(&c1, &c2, "consumer_potential_gathers"), 0);
            if kind == 2 {
                assert_eq!(
                    delta(&c1, &c2, "redistributions"),
                    maps,
                    "crossed plane axes: one redistribution per owned map"
                );
            } else {
                assert_eq!(delta(&c1, &c2, "redistributions"), 0);
            }
        }
    }
}
