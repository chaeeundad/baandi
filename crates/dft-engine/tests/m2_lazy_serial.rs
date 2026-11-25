//! ISSUE-20: `BAANDI_M2_LAZY_SERIAL` (per-thread override). An M2 real Gamma
//! solve state owns only the ORC half-sphere descriptor and the distributed
//! operator; the serial real Hamiltonian is built only when a fallback
//! applies it.
//!
//! Production standard SCF (fixed-point verification) on two simulated ranks
//! (M2 stick rows, kept start, lazy orbitals as by default): lazy on and off
//! are bit-identical (energy, density, eigenvalues, history length), the lazy
//! arm builds no serial operator, and with the distributed operator disabled
//! (test hook) the M1 fallback builds it lazily and again matches the eager
//! arm bit for bit. Own test binary with a lock (process-wide counter).
#[allow(dead_code)]
mod common;
use dft_engine::{GammaReal, ScfOptions, ScfProblem, ScfResult, SpeciesPotential, parallel};
use dft_geometry::{Atom, Cell};
use dft_mpi::{Collective, SimRank, simulate};
use std::sync::{Arc, Mutex};

static LOCK: Mutex<()> = Mutex::new(());

fn problem(pp: &dft_upf::Pseudopotential, skew: bool) -> ScfProblem {
    let rows = if skew {
        [[6.2, 0., 0.], [1.1, 5.7, 0.], [0.4, 0.7, 6.6]]
    } else {
        [[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]]
    };
    let cell = Cell::from_user_rows_bohr(rows).unwrap();
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
            ..ScfOptions::strict_reference()
        },
    )
    .unwrap()
}

fn results(
    ranks: usize,
    pp: &dft_upf::Pseudopotential,
    skew: bool,
    _lazy: bool,
    dist_disabled: bool,
) -> Vec<ScfResult> {
    simulate(ranks, |comm: &SimRank| {
        parallel::set_m2_override(Some(true));
        dft_engine::gamma_real::set_m2_dist_disabled_for_tests(dist_disabled);
        let shared: Arc<dyn Collective + Send + Sync> = Arc::new(comm.clone());
        parallel::set_communicator(Some(shared));
        let r = problem(pp, skew).run().unwrap();
        parallel::set_communicator(None);
        dft_engine::gamma_real::set_m2_dist_disabled_for_tests(false);

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
fn lazy_serial_state_is_bit_identical_and_builds_no_serial_operator() {
    let _lock = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let pp = common::potential(4., 1., true, 401);
    for skew in [false, true] {
        let b0 = dft_engine::gamma_real::serial_operator_builds();
        let eager = results(2, &pp, skew, false, false);
        let b1 = dft_engine::gamma_real::serial_operator_builds();
        let lazy = results(2, &pp, skew, true, false);
        let b2 = dft_engine::gamma_real::serial_operator_builds();
        for r in 0..2 {
            assert_eq!(bits(&eager[r]), bits(&lazy[r]), "skew {skew} rank {r}");
            assert!(lazy[r].certified);
        }
        assert_eq!(bits(&lazy[0]), bits(&lazy[1]), "ranks agree");
        assert_eq!(b1, b0, "default path builds no serial operators");
        assert_eq!(b2, b1, "lazy arm builds no serial operator on the M2 path");
    }
}

#[test]
fn lazy_serial_fallback_builds_the_operator_when_needed() {
    let _lock = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let pp = common::potential(4., 1., true, 401);
    let b0 = dft_engine::gamma_real::serial_operator_builds();
    let eager = results(2, &pp, false, false, true);
    let b1 = dft_engine::gamma_real::serial_operator_builds();
    let lazy = results(2, &pp, false, true, true);
    let b2 = dft_engine::gamma_real::serial_operator_builds();
    for r in 0..2 {
        assert_eq!(bits(&eager[r]), bits(&lazy[r]), "M1 fallback rank {r}");
        assert!(lazy[r].certified);
    }
    assert!(b1 > b0);
    assert!(b2 > b1, "the M1 fallback builds the serial operator lazily");
}
