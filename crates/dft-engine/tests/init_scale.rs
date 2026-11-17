//! INIT-SCALE (scaling of problem construction with MPI ranks): the ionic local
//! coefficients are evaluated by contiguous G blocks per rank and gathered once
//! (always on, bit-identical). Contract checked here: the real-space ionic
//! potential of a problem constructed on 2 to 7 simulated ranks equals the one
//! constructed without a communicator bit for bit on every rank, for a cubic
//! one-atom cell and a skewed two-atom cell, with the default and the QE-like
//! local-mesh convention; the split path is actually taken.
#[allow(dead_code)]
#[path = "../../dft-cli/tests/common/mod.rs"]
mod common;
use dft_engine::{ScfOptions, ScfProblem, SpeciesPotential, parallel};
use dft_geometry::{Atom, Cell};
use dft_mpi::{SimRank, simulate};
use std::sync::Arc;

fn ionic_bits(kind: usize, mesh: dft_upf::radial::LocalMesh) -> Vec<u64> {
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
    let pp = dft_upf::parse_upf(&common::gaussian_with(z, 1., true, 401)).unwrap();
    let problem = ScfProblem::new(
        &cell,
        &atoms,
        &[SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }],
        ScfOptions {
            cutoff_hartree: 3.,
            bands: Some(6),
            threads: 1,
            local_mesh: mesh,
            ..ScfOptions::default()
        },
    )
    .unwrap();
    problem
        .ionic_potential()
        .iter()
        .map(|v| v.to_bits())
        .collect()
}

fn splits() -> usize {
    parallel::m2_stats()
        .0
        .into_iter()
        .find(|(n, _)| n == "ionic_local_splits")
        .map_or(0, |(_, v)| v)
}

#[test]
fn distributed_ionic_local_is_bit_identical_for_every_rank_count() {
    use dft_upf::radial::LocalMesh;
    let before = splits();
    for kind in 0..2 {
        for mesh in [LocalMesh::Baandi, LocalMesh::QeLike] {
            let serial = ionic_bits(kind, mesh);
            for ranks in 2..=7 {
                let got = simulate(ranks, |comm: &SimRank| {
                    parallel::set_communicator(Some(Arc::new(comm.clone())));
                    let bits = ionic_bits(kind, mesh);
                    parallel::set_communicator(None);
                    bits
                });
                for (rank, bits) in got.iter().enumerate() {
                    assert!(
                        *bits == serial,
                        "kind {kind} {mesh:?} P={ranks} rank {rank}: ionic potential bits differ"
                    );
                }
            }
        }
    }
    assert!(splits() > before, "the split ionic local path was not taken");
}
