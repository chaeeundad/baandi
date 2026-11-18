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

/// LOCAL_DENSITY boundary contract: at an active boundary of the owned loop
/// `density()` is None (an exact checkpoint, which reads it, rejects with
/// `Boundary`), `materialize_density()` gathers the input (charge = electrons)
/// and the run continues with the bits of an uninterrupted owned run.
#[test]
fn owned_density_boundary_materializes_on_demand() {
    use dft_engine::{CompletionMode, ScfSession, StepStatus, events::ExecutionConfig};
    let _lock = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let pp = common::potential(4., 1., true, 401);
    let session_run = |interrupt: bool| {
        simulate(2, |comm: &SimRank| {
            parallel::set_m2_override(Some(true));
            dft_engine::set_stop_metric_override(Some(true));
            dft_engine::grid_dist::set_override(Some(true));
            dft_engine::mix_dist::set_override(Some(true));
            dft_engine::grid_dist::set_owned_override(Some((true, true)));
            let shared: Arc<dyn Collective + Send + Sync> = Arc::new(comm.clone());
            parallel::set_communicator(Some(shared));
            let atoms = vec![
                Atom::new("Xx", [0.13, 0.21, 0.32]),
                Atom::new("Xx", [0.61, 0.57, 0.74]),
            ];
            let problem = ScfProblem::new(
                &cell(1),
                &atoms,
                &[SpeciesPotential {
                    label: "Xx",
                    potential: &pp,
                }],
                ScfOptions {
                    gamma_real: GammaReal::Auto,
                    cutoff_hartree: 3.,
                    bands: Some(6),
                    threads: 1,
                    ..ScfOptions::default()
                },
            )
            .unwrap();
            let mut s =
                ScfSession::new(problem, CompletionMode::Strict, ExecutionConfig::default())
                    .unwrap();
            s.initialize().unwrap();
            let mut charge = None;
            let mut steps = 0;
            while s.step().unwrap().status == StepStatus::Ready {
                steps += 1;
                if interrupt && steps == 2 {
                    assert!(s.density().is_none(), "owned boundary: no full density");
                    s.materialize_density().unwrap();
                    let n = s.density().expect("materialized");
                    charge = Some(n.sum());
                }
            }
            parallel::set_communicator(None);
            dft_engine::grid_dist::set_owned_override(None);
            dft_engine::mix_dist::set_override(None);
            dft_engine::grid_dist::set_override(None);
            dft_engine::set_stop_metric_override(None);
            parallel::set_m2_override(None);
            (s.result().map(|o| bits(&o.result)), charge)
        })
    };
    let plain = session_run(false);
    let interrupted = session_run(true);
    assert!(plain[0].0.is_some());
    assert_eq!(
        plain[0].0, interrupted[0].0,
        "materialization changes no bit"
    );
    assert!(interrupted[0].1.is_some_and(|c| c.is_finite() && c > 0.0));
}

/// LOCAL_DENSITY (with LOCAL_FIELDS, M2_GRID, M2_MIX): the output density stays
/// on the slab planes through the audit, the stop rule, the session metrics,
/// the plane-mode Pulay mixer and the next potential. Standard class (5e-9 Ha)
/// against the gathered M2_GRID + M2_MIX path and the threads-only run,
/// verification passed, ranks bit-equal, a fixed rank count reproducible; no
/// per-map density or mixer gather, one materialization per rank (the
/// delivered result).
#[test]
fn owned_density_keeps_the_standard_class_without_density_gathers() {
    let _lock = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let pp = common::potential(4., 1., true, 401);
    for kind in [0, 1, 2] {
        let threads = run(kind, &pp);
        for ranks in [2, 3] {
            let c0 = owned();
            let gathered = ranks_run(ranks, kind, &pp, (true, true, true, false));
            let c1 = owned();
            let local = ranks_run(ranks, kind, &pp, (true, true, true, true));
            let c2 = owned();
            let again = ranks_run(ranks, kind, &pp, (true, true, true, true));
            for r in &local {
                assert!(r.verification.as_ref().is_some_and(|v| v.passed()));
                assert_eq!(bits(r), bits(&local[0]), "ranks agree");
            }
            assert_eq!(bits(&local[0]), bits(&again[0]), "reproducible");
            let d_g = (local[0].total_hartree() - gathered[0].total_hartree()).abs();
            let d_t = (local[0].total_hartree() - threads.total_hartree()).abs();
            let d_n = local[0]
                .density
                .iter()
                .zip(&gathered[0].density)
                .fold(0.0_f64, |m, (a, b)| m.max((a - b).abs()));
            println!(
                "cell {kind} P={ranks}: |dE| vs gathered {d_g:.3e} Ha, vs threads-only {d_t:.3e} Ha, max |dn| {d_n:.3e}, SCF {} (gathered {}); density gathers {} -> {}, mixer gathers {} -> {}, owned density maps {}, materializations {}, owned reductions {}, declined {}",
                local[0].history.len(),
                gathered[0].history.len(),
                delta(&c0, &c1, "density_gathers"),
                delta(&c1, &c2, "density_gathers"),
                delta(&c0, &c1, "mixer_gathers"),
                delta(&c1, &c2, "mixer_gathers"),
                delta(&c1, &c2, "density_owned_maps"),
                delta(&c1, &c2, "density_materializations"),
                delta(&c1, &c2, "owned_reductions"),
                delta(&c1, &c2, "density_declined_maps"),
            );
            assert!(d_g <= 5e-9 && d_t <= 5e-9, "standard class");
            assert!(delta(&c0, &c1, "density_gathers") > 0);
            assert_eq!(delta(&c1, &c2, "density_gathers"), 0, "no density gather");
            // One per rank: the session's one-time mixer validation at
            // initialization (setup); none per map.
            assert_eq!(
                delta(&c1, &c2, "mixer_gathers"),
                ranks,
                "no per-map mixer gather"
            );
            assert_eq!(delta(&c1, &c2, "density_declined_maps"), 0);
            assert!(delta(&c1, &c2, "density_owned_maps") > 0);
            assert_eq!(
                delta(&c1, &c2, "density_materializations"),
                ranks,
                "one materialization per rank: the delivered density"
            );
        }
    }
}
