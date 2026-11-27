//! GK17: bit-identity of the k-path FFT and kept-HX changes on simulated
//! ranks with synthetic pseudopotentials (never QE evidence):
//! * `K_KEPT_SLABS` (the density pass keeps psi_b(r) of the leading bands and
//!   the next map's kept-HX update runs only the forward transform for them):
//!   on (whole budget, and a budget that keeps only some bands) against off;
//! * the one-rank direct stick FFT (GK17) is the only path since CL1 (the
//!   former `FFT_ONE_RANK_DIRECT=off` transpose arm was bit-identical; the
//!   plan-level identity test lives in `dft-solvers` distributed_fft);
//!
//! at P = 1 (one-rank operators), P = 3 with pools on (one-rank operators per
//! k) and P = 3 with pools off (distributed operators: the slab band count
//! must be rank-independent or the update's collectives would not match),
//! with verification on (the verification borrows the cache without
//! storing slabs) and off, on a skewed Gamma-in-list and an elongated shifted
//! mesh. Every delivered bit (energy, density, eigenvalues, orbitals) and
//! the map and H-column counts must be equal.
//!
//! cargo test --release -p dft-engine --test refscf_k_kept_slabs -- --include-ignored --nocapture
#[allow(dead_code)]
mod common;
use dft_engine::{
    GammaReal, ScfOptions, SpeciesPotential, StandardProfile, Verification,
    refscf::{k_driver, profile::Profile},
    refstart::StartInput,
};
use dft_geometry::{
    Atom, Cell,
    kpoints::{KPoint, monkhorst_pack},
};
use dft_mpi::Collective;
use dft_numerics::switches::RunConfig;

fn gamma_in_list() -> Vec<KPoint> {
    vec![
        KPoint {
            fractional: [0.; 3],
            weight: 0.2,
        },
        KPoint {
            fractional: [0.13, -0.21, 0.07],
            weight: 0.3,
        },
        KPoint {
            fractional: [0.25; 3],
            weight: 0.5,
        },
    ]
}

#[derive(PartialEq, Debug)]
struct Fingerprint {
    maps: usize,
    h_columns: usize,
    bits: Vec<u64>,
}

fn run<K: Collective>(
    comm: &K,
    cell: &Cell,
    points: &[KPoint],
    verified: bool,
    pairs: &[(&str, &str)],
) -> Fingerprint {
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    let atoms = [
        Atom::new("Xx", [0.13, 0.21, 0.32]),
        Atom::new("Xx", [0.61, 0.55, 0.71]),
    ];
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: 3.,
        bands: Some(6),
        max_iterations: 100,
        verification: if verified {
            Verification::FixedPoint
        } else {
            Verification::Off
        },
        gamma_real: GammaReal::Off,
        standard: Some(StandardProfile { threshold_ha: 5e-9 }),
        ..ScfOptions::strict_reference()
    };
    let config = RunConfig::from_pairs(pairs).unwrap();
    let _scope = config.enter();
    let run = k_driver::run_with(
        StartInput {
            cell,
            atoms: &atoms,
            potentials: &potentials,
            options: &options,
        },
        Profile {
            threshold_ry: 1e-10,
            ..Profile::default()
        },
        points,
        comm,
    )
    .unwrap();
    assert!(run.outcome.converged, "{}", run.exit_reason);
    let r = &run.outcome.result;
    let mut bits: Vec<u64> = r.density.iter().map(|v| v.to_bits()).collect();
    for k in &r.k_points {
        bits.extend(k.solution.eigenvalues.iter().map(|v| v.to_bits()));
        bits.extend(
            k.solution
                .eigenvectors
                .as_slice()
                .iter()
                .flat_map(|z| [z.re.to_bits(), z.im.to_bits()]),
        );
    }
    bits.push(r.energy.total_hartree().to_bits());
    Fingerprint {
        maps: run.records.len(),
        h_columns: run.records.iter().map(|x| x.h_columns).sum(),
        bits,
    }
}

fn fixtures() -> Vec<(&'static str, Cell, Vec<KPoint>)> {
    vec![
        (
            "skewed gamma-in-list",
            Cell::from_user_rows_bohr([[7.2, 0., 0.], [1.1, 6.8, 0.], [0.4, -0.7, 7.5]]).unwrap(),
            gamma_in_list(),
        ),
        (
            "elongated mp222-shifted",
            Cell::from_user_rows_bohr([[5., 0., 0.], [0., 5., 0.], [0., 0., 11.]]).unwrap(),
            monkhorst_pack([2; 3], [0.5; 3]).unwrap(),
        ),
    ]
}

fn check(
    ranks: usize,
    pools: &'static str,
    verified: bool,
    fixtures: &[(&str, Cell, Vec<KPoint>)],
) {
    // Arms: reference (slabs off), slabs with the whole budget, slabs with a
    // budget below one band per k (keeps none) and a mid budget.
    let arms: [&[(&str, &str)]; 4] = [
        &[("k_kept_slabs", "off")],
        &[("k_kept_slabs", "on")],
        &[("k_kept_slabs", "on"), ("k_kept_slab_mb", "0")],
        &[("k_kept_slabs", "on"), ("k_kept_slab_mb", "1")],
    ];
    for (name, cell, points) in fixtures {
        let mut reference: Option<Vec<Fingerprint>> = None;
        for arm in arms {
            let mut pairs = vec![("refscf_k_pools", pools)];
            pairs.extend_from_slice(arm);
            let out = dft_mpi::simulate(ranks, |comm| run(comm, cell, points, verified, &pairs));
            println!(
                "GK17 {name} P={ranks} pools={pools} verified={verified} {arm:?}: maps {} H {}",
                out[0].maps, out[0].h_columns
            );
            match &reference {
                None => reference = Some(out),
                Some(r) => assert_eq!(r, &out, "{name} P={ranks} pools={pools} {arm:?}"),
            }
        }
    }
}

#[test]
fn kept_slabs_and_direct_fft_are_bit_identical_one_rank() {
    let f = fixtures();
    check(1, "off", false, &f[..1]);
}

#[test]
#[ignore = "release matrix: P = 1, 3 pools on/off, verification on/off, two fixtures"]
fn kept_slabs_and_direct_fft_matrix() {
    let f = fixtures();
    for verified in [false, true] {
        check(1, "off", verified, &f);
        check(3, "on", verified, &f);
        check(3, "off", verified, &f);
    }
}
