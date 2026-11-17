//! BREAKTHROUGH-0-INSTRUMENT: the measurement switches (phase profile,
//! first-block residual diagnostic, SCF snapshot capture) leave the SCF bit for
//! bit unchanged, the snapshot of a solve replays its work exactly, and the
//! profile buckets are filled. Production standard SCF (fixed-point
//! verification) on the threads-only path (one rank, REUSE-PSI kept start) and
//! on the M2 stick-row path (two simulated ranks, kept start). Own test binary
//! with a lock (process-wide forced switches and counters).
#[allow(dead_code)]
mod common;
use dft_engine::{GammaReal, ScfOptions, ScfProblem, ScfResult, SpeciesPotential, parallel};
use dft_geometry::{Atom, Cell};
use dft_mpi::{Collective, SimRank, simulate};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

static LOCK: Mutex<()> = Mutex::new(());

fn problem(pp: &dft_upf::Pseudopotential) -> ScfProblem {
    let cell = Cell::from_user_rows_bohr([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]]).unwrap();
    let atoms = vec![
        Atom::new("Xx", [0.13, 0.21, 0.32]),
        Atom::new("Xx", [0.61, 0.57, 0.74]),
    ];
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
            threads: 1,
            ..ScfOptions::default()
        },
    )
    .unwrap()
}

fn run(pp: &dft_upf::Pseudopotential) -> ScfResult {
    let r = problem(pp).run().unwrap();
    dft_engine::profile::finish_rank(0);
    r
}

fn results(ranks: usize, pp: &dft_upf::Pseudopotential) -> Vec<ScfResult> {
    if ranks == 1 {
        return vec![run(pp)];
    }
    simulate(ranks, |comm: &SimRank| {
        parallel::set_m2_override(Some(true));
        dft_engine::gamma_real::set_m2_reuse_override(Some(true));
        let shared: Arc<dyn Collective + Send + Sync> = Arc::new(comm.clone());
        parallel::set_communicator(Some(shared));
        let r = problem(pp).run().unwrap();
        dft_engine::profile::finish_rank(comm.rank());
        parallel::set_communicator(None);
        dft_engine::gamma_real::set_m2_reuse_override(None);
        parallel::set_m2_override(None);
        r
    })
}

fn bits(r: &ScfResult) -> (u64, Vec<u64>) {
    (
        r.total_hartree().to_bits(),
        r.density.iter().map(|x| x.to_bits()).collect(),
    )
}

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("baandi-instrument-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// A number field of a hand-written manifest.
fn field(text: &str, key: &str) -> String {
    let at = text
        .find(&format!("\"{key}\":"))
        .unwrap_or_else(|| panic!("manifest lacks {key}"));
    let rest = &text[at + key.len() + 3..];
    let end = rest.find([',', '}', '\n']).unwrap();
    rest[..end].trim().trim_matches('"').to_string()
}

fn manifests(dir: &Path) -> Vec<PathBuf> {
    let mut v: Vec<PathBuf> = std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|e| e == "json"))
        .collect();
    v.sort();
    v
}

fn replay_manifest(path: &Path, ranks: usize, pp: &dft_upf::Pseudopotential) {
    let text = std::fs::read_to_string(path).unwrap();
    let label = field(&text, "label");
    let tolerance: f64 = field(&text, "tolerance_ha").parse().unwrap();
    let start = field(&text, "start");
    let grid: [usize; 3] = {
        let at = text.find("\"grid\": [").unwrap();
        let rest = &text[at + 9..];
        let v: Vec<usize> = rest[..rest.find(']').unwrap()]
            .split(',')
            .map(|s| s.trim().parse().unwrap())
            .collect();
        [v[0], v[1], v[2]]
    };
    let dir = path.parent().unwrap();
    let density =
        dft_engine::snapshot::read_density_file(&dir.join(format!("{label}.density.bin")), grid)
            .unwrap();
    let captured_iterations: usize = field(&text, "davidson_iterations").parse().unwrap();
    let captured_h: usize = field(&text, "h_columns").parse().unwrap();
    let captured_fnv = field(&text, "eigen_fnv");
    let solve = |comm: Option<&SimRank>| -> (usize, usize, String) {
        let rank = comm.map_or(0, |c| c.rank());
        let mut p = problem(pp);
        let step = if start == "none" {
            p.step_kpoints_lazy(&density, None, &[], tolerance, false)
                .unwrap()
        } else {
            let (x, hx, _, _) =
                dft_engine::snapshot::read_rank_file(&dir.join(format!("{label}.rank{rank}.bin")))
                    .unwrap();
            p.replay_solve(
                &density,
                dft_engine::snapshot::ReplayStart { x, hx },
                tolerance,
            )
            .unwrap()
        };
        let s = &step.k_points[0].solution;
        (
            s.iterations,
            s.operator_applications,
            format!("{:#018x}", dft_engine::snapshot::fnv(&s.eigenvalues)),
        )
    };
    let outs = if ranks == 1 {
        vec![solve(None)]
    } else {
        simulate(ranks, |comm: &SimRank| {
            parallel::set_m2_override(Some(true));
            dft_engine::gamma_real::set_m2_reuse_override(Some(true));
            let shared: Arc<dyn Collective + Send + Sync> = Arc::new(comm.clone());
            parallel::set_communicator(Some(shared));
            let o = solve(Some(comm));
            parallel::set_communicator(None);
            dft_engine::gamma_real::set_m2_reuse_override(None);
            parallel::set_m2_override(None);
            o
        })
    };
    for (iterations, h, fnv) in outs {
        println!(
            "replay {label} ({ranks} rank(s), start {start}): Davidson {iterations} (captured {captured_iterations}), H columns {h} (captured {captured_h}), eigen FNV {fnv} (captured {captured_fnv})"
        );
        assert_eq!(
            iterations, captured_iterations,
            "{label}: Davidson iterations"
        );
        assert_eq!(h, captured_h, "{label}: H columns");
        assert_eq!(fnv, captured_fnv, "{label}: eigenvalues");
    }
}

#[test]
fn instrumentation_is_bit_identical_and_snapshots_replay() {
    let _lock = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let pp = common::potential(4., 1., true, 401);
    for ranks in [1usize, 2] {
        // Plain run: every switch off.
        dft_engine::profile::force(Some(false));
        dft_engine::first_block::force(Some(false));
        dft_engine::snapshot::force_dir(Some(None));
        let plain = results(ranks, &pp);
        let _ = dft_engine::profile::take_registry();
        let _ = dft_engine::first_block::take_records();
        // Instrumented run: profile, first-block diagnostic and capture on.
        let dir = scratch(&format!("p{ranks}"));
        dft_engine::profile::force(Some(true));
        dft_engine::first_block::force(Some(true));
        dft_engine::snapshot::force_dir(Some(Some(dir.clone())));
        let instrumented = results(ranks, &pp);
        dft_engine::profile::force(Some(false));
        dft_engine::first_block::force(Some(false));
        dft_engine::snapshot::force_dir(Some(None));
        for r in 0..ranks {
            assert_eq!(
                bits(&plain[r]),
                bits(&instrumented[r]),
                "rank {r}: results changed"
            );
            assert!(
                instrumented[r]
                    .verification
                    .as_ref()
                    .is_some_and(|v| v.passed())
            );
        }
        // Profile: one bucket set per rank, maps and solve stages filled.
        let profiles = dft_engine::profile::take_registry();
        assert_eq!(profiles.len(), ranks);
        let report = dft_engine::profile::Report::from_ranks(&profiles);
        print!("{}", report.table(None));
        let fine = &report.max[dft_engine::profile::Bucket::Fine as usize];
        let steps = instrumented[0].history.len();
        assert_eq!(
            fine.counts[0] as usize + report.max[3].counts[0] as usize,
            steps
        );
        assert!(fine.seconds[2] > 0.0, "fine solve seconds");
        let verify = report.max[dft_engine::profile::Bucket::VerifyPass as usize].counts[0]
            + report.max[dft_engine::profile::Bucket::VerifyFail as usize].counts[0];
        assert_eq!(
            verify as usize,
            instrumented[0]
                .verification
                .as_ref()
                .unwrap()
                .attempts
                .len()
        );
        if ranks > 1 {
            assert!(
                fine.counts[12] > 0,
                "tree sums counted inside the fine bucket"
            );
            assert!(fine.seconds[9] > 0.0, "solve.apply seconds on the M2 path");
        }
        // First-block records: one per warm solve, ranks within [1, n_o].
        let records = dft_engine::first_block::take_records();
        print!("{}", dft_engine::first_block::table(&records));
        let warm_solves = steps - 1
            + instrumented[0]
                .history
                .iter()
                .filter(|h| h.standard_redo)
                .count()
            + verify as usize;
        assert_eq!(records.len(), warm_solves, "one record per warm solve");
        for r in &records {
            assert_eq!(r.bands, 6);
            assert!(r.rank_1e1 >= 1 && r.rank_1e1 <= r.occupied);
            assert!(r.rank_1e4 >= r.rank_1e1);
            assert_eq!(r.occupied, 4);
            assert_eq!(r.electrons, 8.);
            assert_eq!(r.mandatory_roots, 5);
            assert_eq!(r.empty_buffer_roots, 1);
            assert_eq!(r.criterion, dft_engine::first_block::Criterion::RitzChange);
            assert_eq!(r.residual_unconverged_all, None);
            assert!(
                r.overlap_defect < 1e-8,
                "kept vectors orthonormal: {}",
                r.overlap_defect
            );
        }
        // Snapshots: one manifest per solve, replayed with the same work.
        let files = manifests(&dir);
        assert_eq!(
            files.len(),
            steps + warm_solves - (steps - 1),
            "one manifest per solve"
        );
        let kept: Vec<_> = files
            .iter()
            .filter(|f| {
                std::fs::read_to_string(f)
                    .unwrap()
                    .contains("\"start\": \"kept\"")
            })
            .collect();
        assert!(!kept.is_empty(), "a kept-start snapshot exists");
        // The cold first map, a kept warm map and the verification map.
        replay_manifest(&files[0], ranks, &pp);
        replay_manifest(kept[kept.len() / 2], ranks, &pp);
        replay_manifest(files.last().unwrap(), ranks, &pp);
        let _ = std::fs::remove_dir_all(&dir);
    }
}
