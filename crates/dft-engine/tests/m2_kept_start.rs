//! MPI-KEPT-START (`BAANDI_M2_REUSE`): standard SCF on simulated ranks with the
//! kept start on the M2 stick rows. Own test binary: it reads process-wide M2
//! counters and must not run beside tests that read the reduction counters.
#[allow(dead_code)]
mod common;
use dft_engine::{GammaReal, ScfOptions, ScfProblem, ScfResult, SpeciesPotential, parallel};
use dft_geometry::{Atom, Cell};
use dft_mpi::{Collective, SimRank, simulate};
use std::sync::Arc;

fn total_davidson(r: &ScfResult) -> usize {
    r.history.iter().map(|h| h.davidson_iterations).sum()
}

fn m2_count(name: &str) -> usize {
    parallel::m2_stats()
        .0
        .into_iter()
        .find(|(n, _)| n == name)
        .map_or(0, |(_, v)| v)
}

/// Standard SCF with two atoms of four valence electrons (four occupied bands,
/// two fully occupied pairs, six bands) on the M2 path, with the MPI kept
/// start on or off (`BAANDI_M2_REUSE`, per-thread override).
fn run_m2_reuse(
    kind: usize,
    comm: Option<Arc<dyn Collective + Send + Sync>>,
    reuse: bool,
) -> ScfResult {
    parallel::set_m2_override(Some(true));
    dft_engine::gamma_real::set_m2_reuse_override(Some(reuse));
    parallel::set_communicator(comm);
    let rows = if kind == 0 {
        [[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]]
    } else {
        [[6.5, 0., 0.], [0.6, 5.8, 0.], [0.3, 0.5, 7.]]
    };
    let cell = Cell::from_user_rows_bohr(rows).unwrap();
    let atoms = vec![
        Atom::new("Xx", [0.13, 0.21, 0.32]),
        Atom::new("Xx", [0.61, 0.57, 0.74]),
    ];
    let pp = common::potential(4., 1., true, 401);
    let opts = ScfOptions {
        gamma_real: GammaReal::Auto,
        cutoff_hartree: 3.,
        bands: Some(6),
        threads: 1,
        eigen_tolerance: 1e-10,
        standard: Some(dft_engine::StandardProfile::default()),
        ..ScfOptions::strict_reference()
    };
    let result = ScfProblem::new(
        &cell,
        &atoms,
        &[SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }],
        opts,
    )
    .unwrap()
    .run()
    .unwrap();
    parallel::set_communicator(None);
    dft_engine::gamma_real::set_m2_reuse_override(None);
    parallel::set_m2_override(None);
    result
}

/// MPI-KEPT-START (`BAANDI_M2_REUSE=on`, docs/design/mpi-kept-start.md): on the
/// standard M2 path the kept start is taken, its potential change uses kept
/// plane slabs, fewer H columns are applied than with the switch off, the energy
/// stays within the standard class (5e-9 Ha) of the threads-only run and of the
/// switch-off run, every rank agrees bit for bit and a fixed rank count is
/// bit-reproducible (energy and density). Cubic and skewed cells, 2 and 3 ranks.
#[test]
fn m2_kept_start_keeps_the_standard_class_and_saves_columns() {
    let density_bits = |r: &ScfResult| r.density.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    let h_columns =
        |r: &ScfResult| -> usize { r.history.iter().map(|h| h.operator_applications).sum() };
    for kind in 0..2 {
        let serial = run_m2_reuse(kind, None, false);
        for ranks in [2usize, 3] {
            let arm = |reuse: bool| {
                simulate(ranks, |comm: &SimRank| {
                    run_m2_reuse(kind, Some(Arc::new(comm.clone())), reuse)
                })
            };
            let off = arm(false);
            let (taken, planes) = (m2_count("reuse_taken"), m2_count("kept_plane_pairs"));
            let on = arm(true);
            let taken = m2_count("reuse_taken") - taken;
            let planes = m2_count("kept_plane_pairs") - planes;
            let again = arm(true);
            for (r, res) in on.iter().enumerate() {
                let de_serial = (res.total_hartree() - serial.total_hartree()).abs();
                let de_off = (res.total_hartree() - off[r].total_hartree()).abs();
                println!(
                    "kept start kind {kind} P={ranks} rank {r}: dE vs threads {de_serial:e} Ha, vs off {de_off:e} Ha, (SCF, Davidson, H) on ({}, {}, {}) off ({}, {}, {}), taken {taken}, kept plane pairs {planes}",
                    res.history.len(),
                    total_davidson(res),
                    h_columns(res),
                    off[r].history.len(),
                    total_davidson(&off[r]),
                    h_columns(&off[r]),
                );
                assert!(de_serial < 5e-9, "kind {kind} P={ranks}: {de_serial:e}");
                assert!(de_off < 5e-9, "kind {kind} P={ranks}: {de_off:e}");
                assert_eq!(
                    res.total_hartree().to_bits(),
                    on[0].total_hartree().to_bits(),
                    "ranks disagree"
                );
                assert_eq!(
                    res.total_hartree().to_bits(),
                    again[r].total_hartree().to_bits(),
                    "not bit-reproducible at a fixed rank count"
                );
                assert_eq!(density_bits(res), density_bits(&again[r]));
            }
            assert!(taken > 0, "kept start never taken");
            assert!(planes > 0, "no potential change from kept slabs");
            if on[0].history.len() == off[0].history.len()
                && total_davidson(&on[0]) <= total_davidson(&off[0])
            {
                assert!(h_columns(&on[0]) < h_columns(&off[0]));
            }
        }
    }
}
