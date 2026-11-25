//! M2-OWNED safe by construction: a run that will need the full density,
//! potential or mixer state at an ACTIVE boundary declines the owned modes
//! (`BAANDI_M2_LOCAL_FIELDS`, `BAANDI_M2_LOCAL_DENSITY`) up front, takes the
//! gathered path with the reason counted (`owned.fields_declined.<reason>`,
//! `owned.density_declined.<reason>`) and never errors at the consumer.
//! Consumers of the DELIVERED result (forces, stress, density output) keep
//! the owned mode: the delivered density and orbitals are gathered at delivery.
//! Simulated ranks; own test binary with a lock (process-wide counters).
#[allow(dead_code)]
mod common;
use dft_engine::checkpoint::{CheckpointKind, ExecutionIdentity, ModelIdentity};
use dft_engine::grid_dist::OwnedDecline;
use dft_engine::{
    CompletionMode, ExecutionConfig, ExecutionMode, Functional, GammaReal, ScfOptions, ScfProblem,
    ScfResult, ScfSession, SpeciesPotential, StepStatus, parallel,
};
use dft_geometry::{Atom, Cell};
use dft_mpi::{Collective, SimRank, simulate};
use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};

static LOCK: Mutex<()> = Mutex::new(());

const POSITIONS: [[f64; 3]; 2] = [[0.13, 0.21, 0.32], [0.61, 0.57, 0.74]];

fn cell() -> Cell {
    Cell::from_user_rows_bohr([[6.5, 0., 0.], [0.6, 5.8, 0.], [0.3, 0.5, 7.]]).unwrap()
}
fn atoms() -> Vec<Atom> {
    POSITIONS.iter().map(|p| Atom::new("Xx", *p)).collect()
}
fn options(functional: Functional) -> ScfOptions {
    ScfOptions {
        gamma_real: GammaReal::Auto,
        cutoff_hartree: 3.,
        bands: Some(6),
        threads: 1,
        functional,
        ..ScfOptions::strict_reference()
    }
}
fn problem(pp: &dft_upf::Pseudopotential, functional: Functional) -> ScfProblem {
    ScfProblem::new(
        &cell(),
        &atoms(),
        &[SpeciesPotential {
            label: "Xx",
            potential: pp,
        }],
        options(functional),
    )
    .unwrap()
}

/// (grid, mix, local fields, local density) on this rank thread around `f`.
fn with_arm<T>(comm: &SimRank, arm: (bool, bool, bool, bool), f: impl FnOnce() -> T) -> T {
    parallel::set_m2_override(Some(true));
    dft_engine::grid_dist::set_override(Some(arm.0));
    dft_engine::mix_dist::set_override(Some(arm.1));
    dft_engine::grid_dist::set_owned_override(Some((arm.2, arm.3)));
    let shared: Arc<dyn Collective + Send + Sync> = Arc::new(comm.clone());
    parallel::set_communicator(Some(shared));
    let out = f();
    parallel::set_communicator(None);
    dft_engine::grid_dist::set_owned_override(None);
    dft_engine::mix_dist::set_override(None);
    dft_engine::grid_dist::set_override(None);
    parallel::set_m2_override(None);
    out
}

const ALL: (bool, bool, bool, bool) = (true, true, true, true);
const GATHERED: (bool, bool, bool, bool) = (true, true, false, false);

type Bits = (u64, Vec<u64>, Vec<u64>, usize);
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
    )
}
fn owned() -> BTreeMap<String, usize> {
    dft_engine::grid_dist::owned_stats().0.into_iter().collect()
}
fn delta(a: &BTreeMap<String, usize>, b: &BTreeMap<String, usize>, key: &str) -> usize {
    b[&format!("owned.{key}")] - a[&format!("owned.{key}")]
}
fn model(electrons: f64) -> (ModelIdentity, ExecutionIdentity) {
    (
        ModelIdentity {
            cell_bohr: cell().lattice_columns_bohr(),
            atoms_fractional: POSITIONS.to_vec(),
            species: vec!["Xx".into(), "Xx".into()],
            pseudopotential_sha256: vec!["0".repeat(64); 2],
            cutoff_hartree: 3.,
            kpoints: vec![[0., 0., 0., 1.]],
            electrons,
            xc: "lda-pz".into(),
            spin: "unpolarized".into(),
            environment: "none".into(),
        },
        ExecutionIdentity {
            build_sha256: "1".repeat(64),
            platform: "synthetic-same-test-build".into(),
            threads: 1,
        },
    )
}
fn drive(s: &mut ScfSession) {
    while s.step().unwrap().status == StepStatus::Ready {}
}

/// Reason `checkpoint`: a declared boundary checkpoint makes both owned modes
/// decline and selects the replicated Pulay mixer, so exact and warm
/// checkpoints of an active boundary succeed (they rejected with `Boundary`
/// or `Mixer` before), an exact restore continues bit for bit, and the
/// result equals the gathered path with the replicated mixer bit for bit.
#[test]
fn checkpoint_request_declines_and_boundary_checkpoints_succeed() {
    let _lock = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let pp = common::potential(4., 1., true, 401);
    let ranks = 2;
    // The same flow (an exact checkpoint drops the kept solver starts, so
    // the reference takes it too) on the gathered path with the replicated
    // mixer and no declaration, and with every owned switch requested.
    let flow = |arm, declare: bool| {
        simulate(ranks, |comm: &SimRank| {
            with_arm(comm, arm, || {
                let mut p = problem(&pp, Functional::LdaPz);
                if declare {
                    p.decline_owned(OwnedDecline::Checkpoint);
                }
                let electrons = p.electron_count();
                let mut s =
                    ScfSession::new(p, CompletionMode::Strict, ExecutionConfig::default()).unwrap();
                s.initialize().unwrap();
                for _ in 0..2 {
                    assert_eq!(s.step().unwrap().status, StepStatus::Ready);
                }
                assert!(s.density().is_some(), "full density at the boundary");
                let (m, e) = model(electrons);
                s.checkpoint(CheckpointKind::WarmStart, m.clone(), e.clone())
                    .expect("warm checkpoint of an active boundary");
                let exact = s
                    .checkpoint(CheckpointKind::ExactContinuation, m.clone(), e.clone())
                    .expect("exact checkpoint of an active boundary");
                drive(&mut s);
                let full = bits(&s.result().unwrap().result);
                // Exact restore on a fresh problem that requests the owned modes:
                // the restore declares the checkpoint itself.
                let mut r =
                    ScfSession::restore_exact(problem(&pp, Functional::LdaPz), &exact, &m, &e)
                        .expect("exact restore");
                drive(&mut r);
                (full, bits(&r.result().unwrap().result))
            })
        })
    };
    let reference = flow((true, false, false, false), false);
    let c0 = owned();
    let runs = flow(ALL, true);
    let c1 = owned();
    for (full, restored) in &runs {
        assert_eq!(full, &runs[0].0, "ranks agree");
        assert_eq!(restored, full, "exact restore continues bit for bit");
        assert_eq!(full, &reference[0].0, "equals the gathered path");
    }
    println!(
        "checkpoint: density_declined.checkpoint {}, fields_declined.checkpoint {}, owned maps {}/{}",
        delta(&c0, &c1, "density_declined.checkpoint"),
        delta(&c0, &c1, "fields_declined.checkpoint"),
        delta(&c0, &c1, "fields_owned_maps"),
        delta(&c0, &c1, "density_owned_maps"),
    );
    // Two sessions per rank (the run and the restore).
    assert_eq!(delta(&c0, &c1, "density_declined.checkpoint"), 2 * ranks);
    assert!(delta(&c0, &c1, "fields_declined.checkpoint") > 0);
    assert_eq!(delta(&c0, &c1, "fields_owned_maps"), 0);
    assert_eq!(delta(&c0, &c1, "density_owned_maps"), 0);
}

/// Reason `runtime_control`: an agent-controlled session declines by itself;
/// the boundary density is readable at every step and the result equals the
/// gathered path bit for bit.
#[test]
fn agent_controlled_session_declines_and_exposes_boundary_density() {
    let _lock = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let pp = common::potential(4., 1., true, 401);
    let ranks = 2;
    let session = |arm, mode| {
        simulate(ranks, |comm: &SimRank| {
            with_arm(comm, arm, || {
                let mut s = ScfSession::new(
                    problem(&pp, Functional::LdaPz),
                    CompletionMode::Strict,
                    ExecutionConfig {
                        mode,
                        ..Default::default()
                    },
                )
                .unwrap();
                s.initialize().unwrap();
                while s.step().unwrap().status == StepStatus::Ready {
                    let n = s.density().expect("boundary density");
                    assert!(n.sum().is_finite());
                }
                bits(&s.result().unwrap().result)
            })
        })
    };
    let gathered = session(GATHERED, ExecutionMode::AgentControlled);
    let c0 = owned();
    let agent = session(ALL, ExecutionMode::AgentControlled);
    let c1 = owned();
    assert_eq!(agent[0], gathered[0], "equals the gathered path");
    assert_eq!(agent[1], agent[0], "ranks agree");
    assert_eq!(delta(&c0, &c1, "density_declined.runtime_control"), ranks);
    assert!(delta(&c0, &c1, "fields_declined.runtime_control") > 0);
    assert_eq!(delta(&c0, &c1, "fields_owned_maps"), 0);
    assert_eq!(delta(&c0, &c1, "density_owned_maps"), 0);
}

/// Reason `functional`: PBE is not on the distributed grid path; with the
/// owned modes requested the run is the gathered run bit for bit.
#[test]
fn pbe_declines_the_owned_mode() {
    let _lock = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let pp = common::potential_with_functional(4., 1., true, 401, "PBE");
    let ranks = 2;
    let arm_run = |arm| {
        simulate(ranks, |comm: &SimRank| {
            with_arm(comm, arm, || {
                bits(&problem(&pp, Functional::Pbe).run().unwrap())
            })
        })
    };
    let g = arm_run(GATHERED);
    let c0 = owned();
    let all = arm_run(ALL);
    let c1 = owned();
    assert_eq!(all[0], g[0], "equals the gathered path");
    assert_eq!(all[1], all[0]);
    assert_eq!(delta(&c0, &c1, "density_declined.functional"), ranks);
    assert!(delta(&c0, &c1, "fields_declined.functional") > 0);
    assert_eq!(delta(&c0, &c1, "fields_owned_maps"), 0);
}

/// Reason `spin`: the collinear driver keeps full densities; with the owned
/// modes requested the run is the gathered run bit for bit.
#[test]
fn spin_declines_the_owned_mode() {
    use dft_engine::collinear::{CollinearProblem, SpinSector};
    let _lock = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let pp = common::potential_with_functional(4., 1., true, 401, "SLA PW");
    let ranks = 2;
    let arm_run = |arm| {
        simulate(ranks, |comm: &SimRank| {
            with_arm(comm, arm, || {
                let r = CollinearProblem::new(
                    &cell(),
                    &atoms(),
                    &[SpeciesPotential {
                        label: "Xx",
                        potential: &pp,
                    }],
                    options(Functional::LdaPw),
                    &[dft_geometry::kpoints::KPoint {
                        fractional: [0.0; 3],
                        weight: 1.0,
                    }],
                    SpinSector {
                        n_up: 5,
                        n_down: 3,
                        axis_cartesian: [0.0, 0.0, 1.0],
                        initial_moment_seed: None,
                    },
                )
                .unwrap()
                .run()
                .unwrap();
                // Debug of f64 is the shortest round-trip form: equal strings
                // are equal bits.
                format!("{:?} {:?} {:?}", r.energy, r.k_points, r.moment_mu_b)
            })
        })
    };
    let g = arm_run(GATHERED);
    let c0 = owned();
    let all = arm_run(ALL);
    let c1 = owned();
    assert_eq!(all[0], g[0], "equals the gathered path");
    assert_eq!(delta(&c0, &c1, "density_declined.spin"), ranks);
    assert_eq!(delta(&c0, &c1, "fields_declined.spin"), ranks);
}
