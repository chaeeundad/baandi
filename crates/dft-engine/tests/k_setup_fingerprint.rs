//! GK16 bit-identity fingerprints (ignored; uses only APIs that exist on main
//! 69107e9e, so the same file builds in a base tree). Prints one line per run:
//! label, ranks, pools, energy bits, map count and a hash of every delivered
//! bit (energy, density, eigenvalues, orbitals). Compare the lines of a base
//! build with those of a candidate build (`diff`). GK16_PAIRS adds switch
//! pairs ("key=value,..."), e.g. `radial_table_fields=off` in the candidate.
//!
//!   cargo test --release -p dft-engine --test k_setup_fingerprint -- --ignored --nocapture
#[allow(dead_code)]
mod common;
use dft_engine::{
    ScfOptions, SpeciesPotential, Verification,
    collinear::SpinSector,
    refscf::{
        electronic::{self, Policy},
        k_driver,
        profile::Profile,
    },
    refstart::StartInput,
};
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use dft_numerics::switches::RunConfig;
use std::hash::{DefaultHasher, Hash, Hasher};

fn kp(f: [f64; 3], w: f64) -> KPoint {
    KPoint {
        fractional: f,
        weight: w,
    }
}

fn pairs(pools: &str) -> Vec<(String, String)> {
    let mut out = vec![("refscf_k_pools".to_string(), pools.to_string())];
    for kv in std::env::var("GK16_PAIRS")
        .unwrap_or_default()
        .split(',')
        .filter(|s| !s.is_empty())
    {
        let (k, v) = kv.split_once('=').expect("key=value");
        out.push((k.into(), v.into()));
    }
    out
}

fn hash_f64<'a>(h: &mut DefaultHasher, values: impl Iterator<Item = &'a f64>) {
    for v in values {
        v.to_bits().hash(h);
    }
}

#[test]
#[ignore = "bit-identity fingerprint; release"]
fn fingerprints() {
    let a = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    let b = common::potential_with_nlcc(true, "SLA PZ");
    let potentials = [
        SpeciesPotential {
            label: "A",
            potential: &a,
        },
        SpeciesPotential {
            label: "B",
            potential: &b,
        },
    ];
    let atoms = [
        Atom::new("A", [0.01, 0.02, 0.03]),
        Atom::new("B", [0.52, 0.47, 0.41]),
        Atom::new("A", [0.27, 0.71, 0.18]),
    ];
    let cells = [
        (
            "skewed",
            [[6.6, 0., 0.], [1.1, 6.2, 0.], [0.4, -0.7, 6.9]],
            vec![
                kp([0.; 3], 0.2),
                kp([0.13, -0.29, 0.21], 0.3),
                kp([0.25, 0.25, -0.25], 0.5),
            ],
        ),
        (
            "elongated_mp",
            [[5., 0., 0.], [0., 5., 0.], [0., 0., 11.]],
            vec![
                kp([0.25, 0.25, 0.25], 0.25),
                kp([0.25, 0.25, -0.25], 0.25),
                kp([0.25, -0.25, 0.25], 0.25),
                kp([-0.25, 0.25, 0.25], 0.25),
            ],
        ),
    ];
    for (label, rows, points) in &cells {
        let cell = Cell::from_user_rows_bohr(*rows).unwrap();
        for (strict, threshold) in [("strict", 1e-12), ("standard", 1e-8)] {
            let options = ScfOptions {
                cutoff_hartree: 3.5,
                bands: Some(7),
                max_iterations: 100,
                verification: Verification::FixedPoint,
                ..Default::default()
            };
            for ranks in [1, 2, 3] {
                for pools in ["off", "on"] {
                    let out = dft_mpi::simulate(ranks, |comm| {
                        let p = pairs(pools);
                        let refs: Vec<(&str, &str)> =
                            p.iter().map(|(k, v)| (k.as_str(), v.as_str())).collect();
                        let config = RunConfig::from_pairs(&refs).unwrap();
                        let _scope = config.enter();
                        let run = k_driver::run_with(
                            StartInput {
                                cell: &cell,
                                atoms: &atoms,
                                potentials: &potentials,
                                options: &options,
                            },
                            Profile {
                                threshold_ry: threshold,
                                ..Profile::default()
                            },
                            points,
                            comm,
                        )
                        .unwrap();
                        let r = &run.outcome.result;
                        let mut h = DefaultHasher::new();
                        hash_f64(&mut h, r.density.iter());
                        for k in &r.k_points {
                            hash_f64(&mut h, k.solution.eigenvalues.iter());
                            for z in k.solution.eigenvectors.as_slice() {
                                z.re.to_bits().hash(&mut h);
                                z.im.to_bits().hash(&mut h);
                            }
                        }
                        let e = r.total_hartree();
                        e.to_bits().hash(&mut h);
                        (e, run.records.len(), h.finish())
                    });
                    for (rank, (e, maps, hash)) in out.iter().enumerate() {
                        println!(
                            "GK16_FP {label} {strict} P={ranks} pools={pools} rank={rank} E={e:.15} bits={:016x} maps={maps} hash={hash:016x}",
                            e.to_bits()
                        );
                    }
                }
            }
        }
    }
    // Spin and Fermi-Dirac (electronic driver), skewed cell, P = 1, 2.
    let a = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PW");
    let b = common::potential_with_nlcc(true, "SLA PW");
    let potentials = [
        SpeciesPotential {
            label: "A",
            potential: &a,
        },
        SpeciesPotential {
            label: "B",
            potential: &b,
        },
    ];
    let atoms = [
        Atom::new("A", [0.01, 0.02, 0.03]),
        Atom::new("B", [0.52, 0.47, 0.41]),
    ];
    let cell =
        Cell::from_user_rows_bohr([[6.6, 0., 0.], [1.1, 6.2, 0.], [0.4, -0.7, 6.9]]).unwrap();
    let points = [kp([0.; 3], 0.25), kp([0.13, -0.29, 0.21], 0.75)];
    let options = ScfOptions {
        cutoff_hartree: 3.,
        bands: Some(6),
        max_iterations: 100,
        verification: Verification::FixedPoint,
        functional: dft_engine::Functional::LdaPw,
        allow_functional_override: true,
        gamma_real: dft_engine::GammaReal::Off,
        ..Default::default()
    };
    let policies = [
        (
            "fixed_spin",
            Policy::FixedSpin(SpinSector {
                n_up: 3,
                n_down: 1,
                axis_cartesian: [0., 0., 1.],
                initial_moment_seed: None,
            }),
        ),
        (
            "fermi_dirac",
            Policy::FermiDirac {
                theta_ha: 0.01,
                band_cap: 12,
            },
        ),
    ];
    for (label, policy) in policies {
        for ranks in [1, 2] {
            for pools in ["off", "on"] {
                let out = dft_mpi::simulate(ranks, |comm| {
                    let p = pairs(pools);
                    let refs: Vec<(&str, &str)> =
                        p.iter().map(|(k, v)| (k.as_str(), v.as_str())).collect();
                    let config = RunConfig::from_pairs(&refs).unwrap();
                    let _scope = config.enter();
                    let r = electronic::run_with_trace(
                        StartInput {
                            cell: &cell,
                            atoms: &atoms,
                            potentials: &potentials,
                            options: &options,
                        },
                        Profile::default(),
                        &points,
                        policy.clone(),
                        comm,
                        None,
                    )
                    .unwrap();
                    let mut h = DefaultHasher::new();
                    for d in &r.density {
                        hash_f64(&mut h, d.iter());
                    }
                    for channel in &r.blocks {
                        for block in channel {
                            format!("{block:?}").hash(&mut h);
                        }
                    }
                    let e = r.energy.total_hartree();
                    e.to_bits().hash(&mut h);
                    (e, r.history.len(), h.finish())
                });
                for (rank, (e, maps, hash)) in out.iter().enumerate() {
                    println!(
                        "GK16_FP {label} P={ranks} pools={pools} rank={rank} E={e:.15} maps={maps} hash={hash:016x}"
                    );
                }
            }
        }
    }
}
