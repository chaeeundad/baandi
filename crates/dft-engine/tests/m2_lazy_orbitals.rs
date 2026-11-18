//! MPI-LAZY-ORBITALS (`BAANDI_M2_LAZY`): the production standard SCF (fixed-point
//! verification) on simulated ranks with the eigenvectors held as distributed
//! stick rows between solves. Own test binary: it reads process-wide M2 and
//! collective counters.
#[allow(dead_code)]
mod common;
use dft_engine::checkpoint::{CheckpointKind, ExecutionIdentity, ModelIdentity};
use dft_engine::{
    CompletionMode, ExecutionConfig, GammaReal, ScfOptions, ScfProblem, ScfResult, ScfSession,
    SpeciesPotential, StepStatus, parallel,
};
use dft_geometry::{Atom, Cell};
use dft_mpi::{Collective, SimRank, simulate};
use std::sync::{Arc, Mutex};

/// Both tests read and bump process-wide counters (lazy materializations,
/// allgatherv bytes): run them one at a time.
static LOCK: Mutex<()> = Mutex::new(());

fn m2_count(name: &str) -> usize {
    parallel::m2_stats()
        .0
        .into_iter()
        .find(|(n, _)| n == name)
        .map_or(0, |(_, v)| v)
}

fn allgatherv_bytes() -> u64 {
    parallel::collective_bytes()
        .into_iter()
        .find(|(n, _)| *n == "allgatherv")
        .map_or(0, |(_, v)| v)
}

fn cell(kind: usize) -> Cell {
    Cell::from_user_rows_bohr(if kind == 0 {
        [[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]]
    } else {
        [[6.5, 0., 0.], [0.6, 5.8, 0.], [0.3, 0.5, 7.]]
    })
    .unwrap()
}

/// Two atoms of four valence electrons (four occupied bands, six bands), the
/// production standard path (verification on), on the M2 path when `comm` has
/// two or more ranks; switches per thread.
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

fn with_switches<R>(
    comm: Option<Arc<dyn Collective + Send + Sync>>,
    lazy: bool,
    reuse: bool,
    f: impl FnOnce() -> R,
) -> R {
    parallel::set_m2_override(Some(true));
    dft_engine::gamma_real::set_m2_lazy_override(Some(lazy));
    dft_engine::gamma_real::set_m2_reuse_override(Some(reuse));
    parallel::set_communicator(comm);
    let r = f();
    parallel::set_communicator(None);
    dft_engine::gamma_real::set_m2_reuse_override(None);
    dft_engine::gamma_real::set_m2_lazy_override(None);
    parallel::set_m2_override(None);
    r
}

fn bits(r: &ScfResult) -> (Vec<u64>, Vec<u64>, Vec<u64>, Vec<u64>) {
    let density = r.density.iter().map(|x| x.to_bits()).collect();
    let orbitals = r
        .orbitals
        .as_slice()
        .iter()
        .flat_map(|z| [z.re.to_bits(), z.im.to_bits()])
        .collect();
    let k_orbitals = r.k_points[0]
        .solution
        .eigenvectors
        .as_slice()
        .iter()
        .flat_map(|z| [z.re.to_bits(), z.im.to_bits()])
        .collect();
    let eigenvalues = r.eigenvalues.iter().map(|x| x.to_bits()).collect();
    (density, orbitals, k_orbitals, eigenvalues)
}

fn work(r: &ScfResult) -> (usize, usize, usize, usize) {
    (
        r.history.len(),
        r.history.iter().map(|h| h.davidson_iterations).sum(),
        r.history.iter().map(|h| h.operator_applications).sum(),
        r.verification.as_ref().map_or(0, |v| v.attempts.len()),
    )
}

/// MPI-LAZY-ORBITALS (docs/design/mpi-lazy-orbitals.md section 4): with the
/// switch on, the production standard SCF on 2 and 3 simulated ranks (cubic and
/// skewed cells, kept start off and on) returns local rows only from every
/// standard solve, starts the next solve from them, gathers each k block once
/// for the result, and delivers bit-identical densities, eigenvalues, orbitals
/// and SCF / Davidson / H / verification counts to the switch-off run; the total
/// energy differs only by the rounding of the distributed kinetic sum (at most
/// 1e-12 Ha here); ranks agree bit for bit and a fixed rank count is
/// bit-reproducible. The allgatherv volume drops.
#[test]
fn m2_lazy_orbitals_match_the_gathered_path() {
    let _lock = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let pp = common::potential(4., 1., true, 401);
    for kind in 0..2 {
        for ranks in [2usize, 3] {
            for reuse in [false, true] {
                let arm = |lazy: bool| {
                    simulate(ranks, |comm: &SimRank| {
                        with_switches(Some(Arc::new(comm.clone())), lazy, reuse, || {
                            problem(kind, &pp).run().unwrap()
                        })
                    })
                };
                let bytes0 = allgatherv_bytes();
                let off = arm(false);
                let bytes_off = allgatherv_bytes() - bytes0;
                let (returns, starts, gathers) = (
                    m2_count("lazy_returns"),
                    m2_count("lazy_starts"),
                    m2_count("lazy_materializations"),
                );
                let bytes0 = allgatherv_bytes();
                let on = arm(true);
                let bytes_on = allgatherv_bytes() - bytes0;
                let returns = m2_count("lazy_returns") - returns;
                let starts = m2_count("lazy_starts") - starts;
                let gathers = m2_count("lazy_materializations") - gathers;
                let again = arm(true);
                let (scf, davidson, h, attempts) = work(&on[0]);
                println!(
                    "lazy kind {kind} P={ranks} reuse {reuse}: (SCF, Davidson, H, verification) ({scf}, {davidson}, {h}, {attempts}), lazy returns {returns} starts {starts} gathers {gathers}, allgatherv bytes off {bytes_off} on {bytes_on}, dE on-off {:e} Ha, kinetic on-off {:e} Ha",
                    on[0].total_hartree() - off[0].total_hartree(),
                    on[0].energy.kinetic - off[0].energy.kinetic,
                );
                for r in 0..ranks {
                    assert_eq!(
                        bits(&on[r]),
                        bits(&off[r]),
                        "rank {r}: state differs from off"
                    );
                    assert_eq!(
                        work(&on[r]),
                        work(&off[r]),
                        "rank {r}: work differs from off"
                    );
                    let de = (on[r].total_hartree() - off[r].total_hartree()).abs();
                    assert!(de <= 1e-12, "rank {r}: dE {de:e} Ha");
                    assert_eq!(
                        on[r].energy.nonlocal.to_bits(),
                        off[r].energy.nonlocal.to_bits()
                    );
                    assert_eq!(bits(&on[r]), bits(&on[0]), "ranks disagree");
                    assert_eq!(
                        on[r].total_hartree().to_bits(),
                        on[0].total_hartree().to_bits()
                    );
                    assert_eq!(bits(&on[r]), bits(&again[r]), "not reproducible");
                    assert_eq!(
                        on[r].total_hartree().to_bits(),
                        again[r].total_hartree().to_bits()
                    );
                }
                // Every standard solve of every rank returned rows only; every
                // warm solve started from them; one gather per rank at the end.
                let solves = scf + attempts;
                assert!(davidson > 0);
                assert!(returns >= ranks * solves, "{returns} < {ranks} x {solves}");
                assert!(starts >= ranks * (solves - 1));
                assert_eq!(gathers, ranks, "one gather per rank for the result");
                // The N x B gathers are gone (in this 6-band fixture the plane
                // and XC gathers of the density dominate the remaining volume).
                assert!(bytes_on < bytes_off, "allgatherv {bytes_on} vs {bytes_off}");
            }
        }
    }
}

/// MPI-LAZY-ORBITALS: a checkpoint taken while the warm orbitals are lazy
/// gathers them (bit-identical to the switch-off capture), and the session
/// continues to the same result afterwards.
#[test]
fn m2_lazy_checkpoint_gathers_the_warm_orbitals() {
    let _lock = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let pp = common::potential(4., 1., true, 401);
    let ranks = 2;
    let capture = |lazy: bool| {
        simulate(ranks, |comm: &SimRank| {
            with_switches(Some(Arc::new(comm.clone())), lazy, false, || {
                let mut problem = problem(1, &pp);
                // M2-OWNED-DECLINE: a boundary checkpoint is declared.
                problem.decline_owned(dft_engine::grid_dist::OwnedDecline::Checkpoint);
                let electrons = problem.electron_count();
                let mut session =
                    ScfSession::new(problem, CompletionMode::Strict, ExecutionConfig::default())
                        .unwrap();
                session.initialize().unwrap();
                for _ in 0..3 {
                    assert_eq!(session.step().unwrap().status, StepStatus::Ready);
                }
                let c = cell(1);
                let model = ModelIdentity {
                    cell_bohr: c.lattice_columns_bohr(),
                    atoms_fractional: vec![[0.13, 0.21, 0.32], [0.61, 0.57, 0.74]],
                    species: vec!["Xx".into(), "Xx".into()],
                    pseudopotential_sha256: vec!["0".repeat(64); 2],
                    cutoff_hartree: 3.,
                    kpoints: vec![[0., 0., 0., 1.]],
                    electrons,
                    xc: "lda-pz".into(),
                    spin: "unpolarized".into(),
                    environment: "none".into(),
                };
                let execution = ExecutionIdentity {
                    build_sha256: "1".repeat(64),
                    platform: "synthetic-same-test-build".into(),
                    threads: 1,
                };
                let cp = session
                    .checkpoint(CheckpointKind::WarmStart, model, execution)
                    .unwrap();
                let orbitals: Vec<u64> = cp
                    .electronic_state
                    .orbitals
                    .iter()
                    .flat_map(|o| o.coefficients.iter().map(|x| x.to_bits()))
                    .collect();
                while session.snapshot().status == StepStatus::Ready {
                    session.step().unwrap();
                }
                let result = session.into_result().unwrap().result;
                (orbitals, bits(&result))
            })
        })
    };
    let gathers = m2_count("lazy_materializations");
    let on = capture(true);
    assert!(m2_count("lazy_materializations") - gathers >= 2 * ranks);
    let off = capture(false);
    for r in 0..ranks {
        assert!(!on[r].0.is_empty());
        assert_eq!(on[r].0, off[r].0, "rank {r}: checkpoint orbitals differ");
        assert_eq!(on[r].1, off[r].1, "rank {r}: continued state differs");
    }
}
