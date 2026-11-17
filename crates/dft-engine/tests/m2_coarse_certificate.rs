//! MPI coarse stage and certificate (docs/performance/mpi-coarse-certificate-analysis.md):
//! the engine M2 switches `BAANDI_M2_SCREEN=dist` (certificate-screen clusters from
//! the distributed band values) and `BAANDI_M2_ATOMIC=split` (atomic start split
//! over ranks), forced per simulated rank with `m2::set_tuning_override`.
//!
//! Contracts checked: every rank agrees bit for bit; a fixed rank count is
//! bit-reproducible; the distributed screen estimates equal the replicated ones
//! within 1e-10 relative (rounding of the plane-local sums); the SCF energy stays
//! in the strict class (1e-10 Ha) of the threads-only run; the split atomic
//! start keeps the Gram-Schmidt arithmetic (bit-identical span) and is in the
//! rounding class of the replicated start through the Gram block GEMM.
#[allow(dead_code)]
#[path = "../../dft-cli/tests/common/mod.rs"]
mod common;
use dft_engine::m2::{EngineTuning, set_tuning_override, take_screen_bounds};
use dft_engine::{GammaReal, ScfOptions, ScfProblem, ScfResult, SpeciesPotential, parallel};
use dft_geometry::{Atom, Cell};
use dft_mpi::{Collective, SimRank, simulate};
use std::sync::Arc;

/// Gaussian norm-conserving fixture with nonlocal projectors and, when
/// `orbitals`, s and p atomic orbitals (PP_PSWFC) for the atomic start.
fn potential(z: f64, orbitals: bool) -> dft_upf::Pseudopotential {
    let mut text = common::gaussian_with(z, 1., true, 401);
    if orbitals {
        let pp = dft_upf::parse_upf(&text).unwrap();
        let channel = |l: u32| {
            pp.mesh()
                .r()
                .iter()
                .map(|r| format!("{}", r.powi(l as i32 + 1) * (-0.7 * r * r).exp()))
                .collect::<Vec<_>>()
                .join(" ")
        };
        text = text.replace("number_of_wfc=\"0\"", "number_of_wfc=\"2\"");
        let chi = format!(
            "<PP_PSWFC><PP_CHI.1 l=\"0\" label=\"1S\" occupation=\"0\">{}</PP_CHI.1><PP_CHI.2 l=\"1\" label=\"2P\" occupation=\"0\">{}</PP_CHI.2></PP_PSWFC>",
            channel(0),
            channel(1)
        );
        text = text.replace("</UPF>", &format!("{chi}</UPF>"));
    }
    dft_upf::parse_upf(&text).unwrap()
}

/// kind 0: one atom in a cubic cell (triply degenerate p-like cluster among the
/// occupied bands); kind 1: two atoms in a skewed cell.
fn run(
    kind: usize,
    orbitals: bool,
    comm: Option<Arc<dyn Collective + Send + Sync>>,
    tuning: EngineTuning,
) -> (ScfResult, Vec<f64>) {
    let m2 = comm.is_some();
    if m2 {
        parallel::set_m2_override(Some(true));
        set_tuning_override(Some(tuning));
    }
    parallel::set_communicator(comm);
    take_screen_bounds();
    let (rows, atoms, z) = if kind == 0 {
        (
            [[5., 0., 0.], [0., 5., 0.], [0., 0., 5.]],
            vec![Atom::new("Xx", [0.13, 0.21, 0.32])],
            8.,
        )
    } else {
        (
            [[6., 0., 0.], [0.5, 5., 0.], [0.2, 0.4, 7.]],
            vec![
                Atom::new("Xx", [0.1, 0.2, 0.3]),
                Atom::new("Xx", [0.55, 0.6, 0.75]),
            ],
            4.,
        )
    };
    let cell = Cell::from_user_rows_bohr(rows).unwrap();
    let pp = potential(z, orbitals);
    let result = ScfProblem::new(
        &cell,
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
            eigen_tolerance: 1e-10,
            ..ScfOptions::strict_reference()
        },
    )
    .unwrap()
    .run()
    .unwrap();
    let bounds = take_screen_bounds();
    parallel::set_communicator(None);
    set_tuning_override(None);
    parallel::set_m2_override(None);
    (result, bounds)
}

fn arm(
    kind: usize,
    orbitals: bool,
    ranks: usize,
    tuning: EngineTuning,
) -> Vec<(ScfResult, Vec<f64>)> {
    simulate(ranks, |comm: &SimRank| {
        run(kind, orbitals, Some(Arc::new(comm.clone())), tuning)
    })
}

fn m2_count(name: &str) -> usize {
    parallel::m2_stats()
        .0
        .into_iter()
        .find(|(n, _)| n == name)
        .map_or(0, |(_, v)| v)
}

fn density_bits(r: &ScfResult) -> Vec<u64> {
    r.density.iter().map(|x| x.to_bits()).collect()
}

/// Every rank equal to rank 0 (energy and density bits) and a second run equal
/// to the first.
fn assert_ranks_and_reproducible(a: &[(ScfResult, Vec<f64>)], b: &[(ScfResult, Vec<f64>)]) {
    for ((r, _), (again, _)) in a.iter().zip(b) {
        assert_eq!(
            r.total_hartree().to_bits(),
            a[0].0.total_hartree().to_bits(),
            "ranks disagree"
        );
        assert_eq!(density_bits(r), density_bits(&a[0].0), "rank densities");
        assert_eq!(
            r.total_hartree().to_bits(),
            again.total_hartree().to_bits(),
            "not bit-reproducible at fixed rank count"
        );
        assert_eq!(density_bits(r), density_bits(again));
    }
}

#[test]
fn distributed_screen_matches_the_replicated_screen() {
    let replicated = EngineTuning::default();
    let dist = EngineTuning {
        screen_dist: true,
        ..replicated
    };
    let before = m2_count("screen_dist_clusters");
    for kind in 0..2 {
        let (serial, _) = run(kind, false, None, replicated);
        for ranks in [2usize, 3] {
            let base = arm(kind, false, ranks, replicated);
            let on = arm(kind, false, ranks, dist);
            let again = arm(kind, false, ranks, dist);
            assert_ranks_and_reproducible(&on, &again);
            for ((r, bounds), (b, base_bounds)) in on.iter().zip(&base) {
                let de = (r.total_hartree() - serial.total_hartree()).abs();
                println!(
                    "screen kind {kind} P={ranks}: dE vs threads-only {de:e} Ha, bounds {bounds:?} vs {base_bounds:?}"
                );
                assert!(de < 1e-10, "kind {kind} P={ranks}: dE {de:e}");
                assert!(!bounds.is_empty(), "screen did not run");
                assert_eq!(bounds.len(), base_bounds.len());
                for (x, y) in bounds.iter().zip(base_bounds) {
                    assert!(
                        (x - y).abs() <= 1e-10 * y.abs().max(1e-12),
                        "screen estimate {x:e} vs replicated {y:e}"
                    );
                }
                // The screen only decides whether the certificate solve runs:
                // equal decisions give the same state bits as the replicated arm.
                assert_eq!(r.history.len(), b.history.len());
                assert_eq!(r.total_hartree().to_bits(), b.total_hartree().to_bits());
            }
        }
    }
    assert!(
        m2_count("screen_dist_clusters") > before,
        "no cluster took the distributed screen"
    );
}

#[test]
fn split_atomic_start_matches_the_replicated_start() {
    let replicated = EngineTuning::default();
    let split = EngineTuning {
        atomic_split: true,
        ..replicated
    };
    let before = m2_count("atomic_split_starts");
    for kind in 0..2 {
        let (serial, _) = run(kind, true, None, replicated);
        assert!(
            serial.atomic_start.initial_h_applications > 0,
            "atomic start not used"
        );
        for ranks in [2usize, 3, 5] {
            let base = arm(kind, true, ranks, replicated);
            let on = arm(kind, true, ranks, split);
            let again = arm(kind, true, ranks, split);
            assert_ranks_and_reproducible(&on, &again);
            for ((r, _), (b, _)) in on.iter().zip(&base) {
                let de = (r.total_hartree() - serial.total_hartree()).abs();
                let db = (r.total_hartree() - b.total_hartree()).abs();
                println!(
                    "atomic kind {kind} P={ranks}: dE vs threads-only {de:e} Ha, vs replicated M2 {db:e} Ha"
                );
                assert!(de < 1e-10, "kind {kind} P={ranks}: dE {de:e}");
                assert_eq!(
                    r.atomic_start.atomic_orbitals,
                    b.atomic_start.atomic_orbitals
                );
                assert_eq!(
                    r.atomic_start.initial_h_applications,
                    b.atomic_start.initial_h_applications
                );
                // Rounding class: the Gram block GEMM of a column subset may
                // round differently from the whole product (observed: one ulp of
                // the total energy at P=5 on the pure build).
                assert!(db < 1e-10, "kind {kind} P={ranks}: dE {db:e}");
            }
        }
    }
    assert!(m2_count("atomic_split_starts") > before);
}

/// INIT-SCALE `BAANDI_M2_ATOMIC_BLOCK=on` (default off): row-distributed CholQR2
/// span and row-split Ritz rotation. Rounding class of the split start: every
/// rank agrees bit for bit, a fixed rank count reproduces its bits, the energy
/// stays in the strict class of the threads-only run and of the split arm, the
/// block path is taken (no fallback on these well-conditioned spans).
#[test]
fn block_atomic_start_stays_in_the_class_of_the_split_start() {
    let split = EngineTuning {
        atomic_split: true,
        ..EngineTuning::default()
    };
    let block = EngineTuning {
        atomic_block: true,
        ..split
    };
    let (starts, fallbacks) = (
        m2_count("atomic_block_starts"),
        m2_count("atomic_block_fallbacks"),
    );
    for kind in 0..2 {
        let (serial, _) = run(kind, true, None, split);
        for ranks in [2usize, 3, 5, 7] {
            let base = arm(kind, true, ranks, split);
            let on = arm(kind, true, ranks, block);
            let again = arm(kind, true, ranks, block);
            assert_ranks_and_reproducible(&on, &again);
            for ((r, _), (b, _)) in on.iter().zip(&base) {
                let de = (r.total_hartree() - serial.total_hartree()).abs();
                let db = (r.total_hartree() - b.total_hartree()).abs();
                println!(
                    "block atomic kind {kind} P={ranks}: dE vs threads-only {de:e} Ha, vs split {db:e} Ha"
                );
                assert!(de < 1e-10, "kind {kind} P={ranks}: dE {de:e}");
                assert!(db < 1e-10, "kind {kind} P={ranks}: dE {db:e}");
                assert_eq!(
                    r.atomic_start.initial_h_applications,
                    b.atomic_start.initial_h_applications
                );
            }
        }
    }
    assert!(m2_count("atomic_block_starts") > starts);
    assert_eq!(m2_count("atomic_block_fallbacks"), fallbacks);
}
