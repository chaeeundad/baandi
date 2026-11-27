//! GK4: deferred per-k energy terms are the only path since CL1 (the eager
//! per-map arm `K_LEAN_ENERGY=off` was retired; its bit identity in
//! everything delivered and in every stop input is recorded in the attempts
//! ledger GK4 and docs/performance/switch-retirement.md). This test keeps the
//! lean path bit-reproducible, converged, with null intermediate slots.
//! Synthetic data, never QE evidence.
#[allow(dead_code)]
mod common;
use dft_engine::{
    ScfOptions, SpeciesPotential, Verification,
    refscf::{driver::Run, k_driver, profile::Profile},
    refstart::StartInput,
};
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use dft_mpi::Collective;

fn bits(v: &[f64]) -> Vec<u64> {
    v.iter().map(|x| x.to_bits()).collect()
}

fn energy_bits(e: &dft_engine::energy::EnergyAudit) -> Vec<u64> {
    bits(&[
        e.kinetic,
        e.local,
        e.local_zero_mode,
        e.nonlocal,
        e.hartree,
        e.xc,
        e.ewald,
        e.band_sum,
        e.double_counting,
        e.band_form,
        e.harris_foulkes,
        e.functional_difference_bound,
        e.total_hartree(),
    ])
}

fn assert_identical(a: &Run, b: &Run, label: &str) {
    let (ar, br) = (&a.outcome.result, &b.outcome.result);
    assert_eq!(a.outcome.converged, b.outcome.converged, "{label}");
    assert_eq!(a.exit_reason, b.exit_reason, "{label}");
    assert_eq!(energy_bits(&ar.energy), energy_bits(&br.energy), "{label}");
    assert!(br.energy.total_hartree().is_finite(), "{label}");
    assert_eq!(
        format!("{:?}", ar.energy.nlcc),
        format!("{:?}", br.energy.nlcc),
        "{label}"
    );
    let density = |r: &dft_engine::ScfResult| bits(r.density.as_slice().unwrap());
    assert_eq!(density(ar), density(br), "{label} density");
    assert_eq!(ar.k_points.len(), br.k_points.len());
    for (x, y) in ar.k_points.iter().zip(&br.k_points) {
        assert_eq!(bits(&x.solution.eigenvalues), bits(&y.solution.eigenvalues));
        assert_eq!(bits(&x.solution.residuals), bits(&y.solution.residuals));
        assert_eq!(x.solution.converged, y.solution.converged);
        assert_eq!(x.solution.eigenvectors, y.solution.eigenvectors);
    }
    assert_eq!(a.records.len(), b.records.len(), "{label} map count");
    assert_eq!(a.first_reference_stop, b.first_reference_stop);
    assert_eq!(a.verification_calls.len(), b.verification_calls.len());
    for (x, y) in a.records.iter().zip(&b.records) {
        assert_eq!(x.dr2_ry.to_bits(), y.dr2_ry.to_bits(), "{label}");
        assert_eq!(x.ethr_ry.to_bits(), y.ethr_ry.to_bits());
        assert_eq!(x.h_columns, y.h_columns);
        assert_eq!(x.davidson_iterations, y.davidson_iterations);
        assert_eq!(x.descf_ha.to_bits(), y.descf_ha.to_bits());
        assert_eq!(
            x.variational_energy_ha.to_bits(),
            y.variational_energy_ha.to_bits()
        );
        assert_eq!(x.energy.band_sum.to_bits(), y.energy.band_sum.to_bits());
        assert_eq!(
            x.energy.double_counting.to_bits(),
            y.energy.double_counting.to_bits()
        );
        assert_eq!(
            x.energy.harris_foulkes.to_bits(),
            y.energy.harris_foulkes.to_bits()
        );
        // The four deferred slots are null (NaN) together; compare the rest.
        if y.energy.kinetic.is_nan() {
            assert!(y.energy.local.is_nan() && y.energy.nonlocal.is_nan());
            assert!(y.energy.local_zero_mode.is_nan());
            assert!(x.energy.kinetic.is_nan());
        } else {
            assert_eq!(energy_bits(&x.energy), energy_bits(&y.energy));
        }
    }
    for (x, y) in ar.history.iter().zip(&br.history) {
        assert_eq!(
            x.max_eigen_residual.to_bits(),
            y.max_eigen_residual.to_bits()
        );
        assert_eq!(x.density_residual.to_bits(), y.density_residual.to_bits());
    }
    assert_eq!(
        format!("{:?}", a.map_diagnostics),
        format!("{:?}", b.map_diagnostics)
    );
    assert_eq!(
        format!("{:?}", ar.verification.as_ref().map(|v| v.attempts.len())),
        format!("{:?}", br.verification.as_ref().map(|v| v.attempts.len()))
    );
    assert_eq!(
        energy_bits(&ar.history.last().unwrap().energy),
        energy_bits(&br.history.last().unwrap().energy)
    );
}

#[test]
fn k_lean_energy_is_bit_reproducible() {
    let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
    let atoms = [
        Atom::new("Xx", [0.13, 0.21, 0.32]),
        Atom::new("Xx", [0.61, 0.55, 0.71]),
    ];
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let points = [
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
    ];
    for rows in [
        [[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]],
        [[5., 0., 0.], [0., 5., 0.], [0., 0., 11.]],
    ] {
        let cell = Cell::from_user_rows_bohr(rows).unwrap();
        for (verification, max_iterations) in [
            (Verification::FixedPoint, 100),
            (Verification::Off, 100),
            (Verification::Off, 2),
        ] {
            for ranks in [1, 3] {
                dft_mpi::simulate(ranks, |comm| {
                    let options = ScfOptions {
                        cutoff_hartree: 1.,
                        bands: Some(4),
                        max_iterations,
                        verification,
                        ..Default::default()
                    };
                    let run = || {
                        k_driver::run_with(
                            StartInput {
                                cell: &cell,
                                atoms: &atoms,
                                potentials: &potentials,
                                options: &options,
                            },
                            Profile::default(),
                            &points,
                            comm,
                        )
                        .unwrap()
                    };
                    let off = run();
                    let on = run();
                    let label = format!("{rows:?} {verification:?} max={max_iterations} P={ranks}");
                    assert_identical(&off, &on, &label);
                    if max_iterations == 100 {
                        assert!(on.outcome.converged, "{label}: {}", on.exit_reason);
                    }
                    // Ordinary intermediate maps carry null slots.
                    if on.records.len() > 1 {
                        assert!(on.records[0].energy.kinetic.is_nan(), "{label}");
                    }
                    if comm.rank() == 0 {
                        println!(
                            "GK4 {label}: E={:.17e} Ha maps={}",
                            on.outcome.result.total_hartree(),
                            on.records.len()
                        );
                    }
                });
            }
        }
    }
}
