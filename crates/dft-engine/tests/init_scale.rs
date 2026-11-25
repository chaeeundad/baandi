//! INIT-SCALE (scaling of problem construction with MPI ranks): the ionic local
//! coefficients (contiguous G blocks per rank, one allgatherv) and the Ewald
//! sum (contiguous candidate blocks per rank, accepted terms gathered and added
//! in the serial order) are split over the ranks, always on, bit-identical.
//! Contract checked here: the real-space ionic potential and every Ewald field
//! of a problem constructed on 2 to 7 simulated ranks equal those constructed
//! without a communicator bit for bit on every rank, for a cubic one-atom cell
//! and a skewed two-atom cell, with the default and the QE-like local-mesh
//! convention (Ewald cache off, so each construction evaluates); both split
//! paths are actually taken. INIT-EWALD-GATE: the Ewald split runs only from
//! `EWALD_SPLIT_MIN_RANKS` ranks on; below it every rank evaluates the serial
//! sum (no split counted), with the same bits on both sides of the gate.
//! INIT-DENSITY-SPLIT: the atomic-superposition start density (opt-in
//! `start_density = Atomic`, a Gaussian atomic charge substituted for the
//! synthetic zero PP_RHOATOM) is split by G blocks too; its grid values are
//! bit-identical on every rank count and the split is counted.
//! INIT-EWALD-FAST / INIT-IONIC-FAST: the same contract with the balanced Ewald
//! scheme and the phase-table structure factor (per-thread overrides set in
//! every rank thread), against the serial construction with the same switches.
#[allow(dead_code)]
#[path = "../../dft-cli/tests/common/mod.rs"]
mod common;
use dft_engine::{ScfOptions, ScfProblem, SpeciesPotential, parallel};
use dft_geometry::{Atom, Cell};
use dft_mpi::{SimRank, simulate};
use std::sync::Arc;

/// The synthetic Gaussian UPF with a normalized Gaussian atomic charge
/// (4 pi r^2 rho(r), rho = z pi^(-3/2) exp(-r^2)) in PP_RHOATOM, on the mesh
/// `common::gaussian_with` uses (r0 = 0.02, end 12 Bohr).
fn with_atomic_charge(z: f64, points: usize) -> String {
    let text = common::gaussian_with(z, 1., true, points);
    let dx = (1.0_f64 + 12.0 / 0.02).ln() / (points - 1) as f64;
    let values: Vec<String> = (0..points)
        .map(|i| {
            let r = 0.02 * (i as f64 * dx).exp_m1();
            let rho = z * std::f64::consts::PI.powf(-1.5) * (-r * r).exp();
            format!("{:e}", 4.0 * std::f64::consts::PI * r * r * rho)
        })
        .collect();
    let open = text.find("<PP_RHOATOM>").unwrap() + "<PP_RHOATOM>".len();
    let close = text.find("</PP_RHOATOM>").unwrap();
    format!("{}{}{}", &text[..open], values.join(" "), &text[close..])
}

type Bits = (Vec<u64>, Vec<u64>, Vec<u64>);

fn construction_bits(kind: usize, mesh: dft_upf::radial::LocalMesh, fast: bool) -> Bits {
    dft_physics::ewald::set_scheme_override(Some(if fast {
        dft_physics::ewald::Scheme::Balanced
    } else {
        dft_physics::ewald::Scheme::Reference
    }));
    dft_physics::hamiltonian::set_ionic_fast_override(Some(fast));
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
    let pp = dft_upf::parse_upf(&with_atomic_charge(z, 401)).unwrap();
    dft_engine::set_ewald_cache_override(Some(false));
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
            start_density: dft_engine::StartDensity::Atomic,
            ..ScfOptions::strict_reference()
        },
    )
    .unwrap();
    dft_engine::set_ewald_cache_override(None);
    dft_physics::ewald::set_scheme_override(None);
    dft_physics::hamiltonian::set_ionic_fast_override(None);
    let e = problem.ewald_terms();
    let ewald = [
        e.energy.0,
        e.real.0,
        e.reciprocal.0,
        e.self_energy.0,
        e.background.0,
        e.eta,
        e.real_cutoff_bohr,
        e.reciprocal_cutoff_inverse_bohr,
        e.tail_bound.0,
    ];
    (
        problem
            .ionic_potential()
            .iter()
            .map(|v| v.to_bits())
            .collect(),
        ewald.iter().map(|v| v.to_bits()).collect(),
        problem
            .atomic_superposition_density()
            .expect("atomic start density")
            .iter()
            .map(|v| v.to_bits())
            .collect(),
    )
}

fn count(name: &str) -> usize {
    parallel::m2_stats()
        .0
        .into_iter()
        .find(|(n, _)| n == name)
        .map_or(0, |(_, v)| v)
}

#[test]
fn distributed_construction_is_bit_identical_for_every_rank_count() {
    use dft_upf::radial::LocalMesh;
    let (ionic_before, ewald_before) = (count("ionic_local_splits"), count("ewald_splits"));
    let density_before = count("density_splits");
    for (kind, fast) in [(0, false), (1, false), (0, true), (1, true)] {
        for mesh in [LocalMesh::Baandi, LocalMesh::QeLike] {
            let serial = construction_bits(kind, mesh, fast);
            for ranks in 2..=7 {
                let ewald_before_p = count("ewald_splits");
                let got = simulate(ranks, |comm: &SimRank| {
                    parallel::set_communicator(Some(Arc::new(comm.clone())));
                    let bits = construction_bits(kind, mesh, fast);
                    parallel::set_communicator(None);
                    bits
                });
                let split = count("ewald_splits") > ewald_before_p;
                assert_eq!(
                    split,
                    dft_engine::ewald_split_applies(ranks),
                    "kind {kind} fast {fast} {mesh:?} P={ranks}: Ewald split taken = {split}"
                );
                assert_eq!(split, ranks >= dft_engine::EWALD_SPLIT_MIN_RANKS);
                for (rank, (ionic, ewald, start)) in got.iter().enumerate() {
                    assert!(
                        *start == serial.2,
                        "kind {kind} {mesh:?} P={ranks} rank {rank}: start density bits differ"
                    );
                    assert!(
                        *ionic == serial.0,
                        "kind {kind} {mesh:?} P={ranks} rank {rank}: ionic potential bits differ"
                    );
                    assert_eq!(
                        *ewald, serial.1,
                        "kind {kind} {mesh:?} P={ranks} rank {rank}: Ewald bits differ"
                    );
                }
            }
        }
    }
    assert!(count("ionic_local_splits") > ionic_before);
    assert!(count("ewald_splits") > ewald_before);
    assert!(count("density_splits") > density_before);
}
