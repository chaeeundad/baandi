//! GD2: joint force + stress pass (`DERIV_JOINT`) and the split derivative
//! PBE evaluation (`DERIV_XC_SPLIT`) are bit-identical to the former separate,
//! replicated passes, for every simulated rank count, with and without the
//! negative-density screening. Synthetic PBE + NLCC fixture (not QE evidence).
//! Release: cargo test --release -p dft-engine --test gd2_derivatives -- --nocapture
#[allow(dead_code)]
mod common;
use dft_engine::{
    GammaReal, ScfOptions, SpeciesPotential, StartDensity, StartWavefunctions,
    refscf::{
        api::{Derivatives, Request, solve},
        delivered,
        electronic::Policy,
        k_driver,
        profile::Profile,
    },
    refstart::StartInput,
};
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use dft_numerics::switches::RunConfig;
use std::sync::Arc;

fn options() -> ScfOptions {
    ScfOptions {
        standard: Some(dft_engine::StandardProfile {
            threshold_ha: 5e-13,
        }),
        cutoff_hartree: 3.,
        bands: Some(4),
        gamma_real: GammaReal::Off,
        functional: dft_engine::Functional::Pbe,
        start_density: StartDensity::Atomic,
        start_wavefunctions: StartWavefunctions::Atomic,
        max_iterations: 200,
        ..ScfOptions::default()
    }
}

/// Switch arms: `legacy` = all GD2 switches off (separate replicated passes).
fn arm(legacy: bool, screening: bool) -> RunConfig {
    let screen = if screening { "on" } else { "off" };
    if legacy {
        RunConfig::from_pairs(&[
            ("xc_screen_negative", screen),
            ("deriv_joint", "off"),
            ("deriv_xc_split", "off"),
            ("deriv_shell_tables", "off"),
        ])
    } else {
        RunConfig::from_pairs(&[("xc_screen_negative", screen)])
    }
    .unwrap()
}

type Bits = Vec<u64>;
fn force_bits(f: &dft_engine::forces::ForceResult) -> Bits {
    [&f.values, &f.ewald, &f.local, &f.nonlocal, &f.core]
        .into_iter()
        .flat_map(|v| v.iter().flatten().map(|x| x.to_bits()))
        .collect()
}
fn stress_bits(s: &dft_engine::stress::StressResult) -> Bits {
    s.values
        .iter()
        .flatten()
        .chain(s.terms.iter().flat_map(|(_, t)| t.iter().flatten()))
        .map(|x| x.to_bits())
        .collect()
}

#[test]
fn joint_pass_and_split_xc_match_separate_replicated_passes() {
    let pp = common::potential_with_nlcc(true, "PBE");
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let cell = Cell::from_user_rows_bohr([[7., 0., 0.], [0.4, 6.5, 0.], [0.3, 0.2, 7.5]]).unwrap();
    let atoms = [
        Atom::new("Xx", [0.1, 0.2, 0.3]),
        Atom::new("Xx", [0.37, 0.43, 0.61]),
    ];
    let opts = options();
    let points: Vec<_> = [
        [0., 0., 0.],
        [0.25, 0., 0.],
        [0., 0.25, 0.],
        [0., 0., 0.25],
        [0.25, 0.25, 0.],
    ]
    .into_iter()
    .zip([0.1, 0.2, 0.2, 0.2, 0.3])
    .map(|(fractional, weight)| KPoint { fractional, weight })
    .collect();
    let profile = Profile::from_standard_ha(opts.standard.unwrap().threshold_ha).unwrap();
    let input = || StartInput {
        cell: &cell,
        atoms: &atoms,
        potentials: &species,
        options: &opts,
    };
    for screening in [false, true] {
        let legacy_scope = arm(true, screening);
        let new_scope = arm(false, screening);
        let mut run = {
            let _s = legacy_scope.enter();
            k_driver::run(input(), profile, &points).unwrap()
        };
        let view = delivered::force_state(input(), &mut run.outcome).unwrap();
        let (reference_f, reference_s) = {
            let _s = legacy_scope.enter();
            let f = dft_engine::forces(&view).unwrap();
            let s = dft_engine::stress::frozen_stress(&view).unwrap();
            assert!(
                dft_engine::stress::forces_and_frozen_stress(&view)
                    .unwrap()
                    .is_none()
            );
            (force_bits(&f), stress_bits(&s))
        };
        {
            let _s = new_scope.enter();
            let f = dft_engine::forces(&view).unwrap();
            let s = dft_engine::stress::frozen_stress(&view).unwrap();
            assert_eq!(
                force_bits(&f),
                reference_f,
                "split forces screening={screening}"
            );
            assert_eq!(
                stress_bits(&s),
                reference_s,
                "split stress screening={screening}"
            );
            let (jf, js) = dft_engine::stress::forces_and_frozen_stress(&view)
                .unwrap()
                .unwrap();
            assert_eq!(
                force_bits(&jf),
                reference_f,
                "joint forces screening={screening}"
            );
            assert_eq!(
                stress_bits(&js),
                reference_s,
                "joint stress screening={screening}"
            );
            assert_eq!(jf.band_passes, f.band_passes);
            assert_eq!(js.band_passes, s.band_passes);
            assert_eq!(js.timings, vec![("joint_with_forces", 0.)]);
        }
        for ranks in [2, 3, 7] {
            let results = dft_mpi::simulate(ranks, |comm| {
                let _s = new_scope.enter();
                dft_engine::parallel::set_communicator(Some(Arc::new(comm.clone())));
                let f = dft_engine::forces(&view).unwrap();
                let s = dft_engine::stress::frozen_stress(&view).unwrap();
                let (jf, js) = dft_engine::stress::forces_and_frozen_stress(&view)
                    .unwrap()
                    .unwrap();
                dft_engine::parallel::set_communicator(None);
                (
                    force_bits(&f),
                    stress_bits(&s),
                    force_bits(&jf),
                    stress_bits(&js),
                )
            });
            for (rank, (f, s, jf, js)) in results.into_iter().enumerate() {
                let tag = format!("ranks={ranks} rank={rank} screening={screening}");
                assert_eq!(f, reference_f, "{tag}");
                assert_eq!(s, reference_s, "{tag}");
                assert_eq!(jf, reference_f, "{tag}");
                assert_eq!(js, reference_s, "{tag}");
            }
        }
        println!(
            "GD2 scalar k fixture screening={screening}: legacy, split, joint and 2/3/7 ranks identical"
        );
    }
}

/// The deferred derivative delivery of a Fermi-Dirac reference run (the
/// production path) computes forces and stress in one joint pass; both equal
/// the separate passes of the legacy arm bit for bit.
#[test]
fn fd_delivery_joint_matches_separate() {
    let pp = common::potential_with_nlcc(true, "PBE");
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let cell = Cell::from_user_rows_bohr([[7., 0., 0.], [0.3, 7., 0.], [0., 0.2, 7.]]).unwrap();
    let atoms = [
        Atom::new("Xx", [0.1, 0.2, 0.3]),
        Atom::new("Xx", [0.37, 0.43, 0.49]),
    ];
    let points = vec![
        KPoint {
            fractional: [0.; 3],
            weight: 0.5,
        },
        KPoint {
            fractional: [0.25, 0.25, 0.],
            weight: 0.5,
        },
    ];
    let run = |legacy: bool| {
        let _s = arm(legacy, true).enter();
        let threshold_ha = 5e-11;
        let mut opts = options();
        opts.bands = Some(10);
        opts.standard = Some(dft_engine::StandardProfile { threshold_ha });
        let mut request = Request::new(
            Profile::from_standard_ha(threshold_ha).unwrap(),
            points.clone(),
        );
        request.policy = Some(Policy::FermiDirac {
            theta_ha: 0.004,
            band_cap: 10,
        });
        request.derivatives = Derivatives {
            forces: true,
            stress: true,
            deferred: false,
            ..Default::default()
        };
        let input = StartInput {
            cell: &cell,
            atoms: &atoms,
            potentials: &species,
            options: &opts,
        };
        let d = solve(input, request).unwrap();
        let f = d.derivatives.forces.unwrap();
        let s = d.derivatives.stress.unwrap();
        let r = d.run.into_electronic().unwrap();
        (
            r.thermodynamics.free_energy_ha.to_bits(),
            force_bits(&f),
            stress_bits(&s),
            s.timings.clone(),
        )
    };
    let (e0, f0, s0, t0) = run(true);
    let (e1, f1, s1, t1) = run(false);
    assert_eq!(e0, e1);
    assert_eq!(f0, f1);
    assert_eq!(s0, s1);
    assert_ne!(t0, vec![("joint_with_forces", 0.)]);
    assert_eq!(t1, vec![("joint_with_forces", 0.)]);
}

/// GD2 (`DERIV_SHELL_TABLES`): on a cell with more distinct |G| than table
/// nodes the tabulated shell values are independent of the rank count (bit
/// for bit) and agree with the per-shell quadrature to roundoff; on a cell
/// with fewer shells than nodes the exact values are kept.
#[test]
fn shell_tables_rank_independent_and_close_to_quadrature() {
    use dft_physics::{hamiltonian::Species, radial_shells::ShellKind};
    let pp = common::potential_with_nlcc(true, "PBE");
    let cell = Cell::from_user_rows_bohr([[9., 0., 0.], [0.7, 8.5, 0.], [0.3, 0.9, 10.]]).unwrap();
    let positions = [
        cell.fractional_to_cartesian_bohr([0.1, 0.2, 0.3]),
        cell.fractional_to_cartesian_bohr([0.6, 0.45, 0.8]),
    ];
    let species: Vec<_> = positions
        .iter()
        .map(|p| Species {
            potential: &pp,
            positions_bohr: std::slice::from_ref(p),
        })
        .collect();
    let volume = cell.volume_bohr3();
    let kinds = [
        ShellKind::LocalValue,
        ShellKind::LocalSlope,
        ShellKind::CoreValue,
        ShellKind::CoreSlope,
    ];
    for (cutoff, tabulated) in [(10., true), (2., false)] {
        let db = dft_basis::PlaneWaveBasis::at_k(&cell, dft_units::Hartree(4. * cutoff), [0.; 3])
            .unwrap();
        let qs: Vec<f64> = db
            .vectors()
            .iter()
            .filter(|g| g.miller != [0; 3])
            .map(|g| g.norm_squared.sqrt())
            .collect();
        let exact = {
            let _s = RunConfig::from_pairs(&[("deriv_shell_tables", "off")])
                .unwrap()
                .enter();
            dft_engine::forces::shell_values(&species, &db, volume, true)
                .unwrap()
                .unwrap()
        };
        let serial = dft_engine::forces::shell_values(&species, &db, volume, true)
            .unwrap()
            .unwrap();
        let (mut worst, mut peak, mut identical) = ([0f64; 4], [0f64; 4], true);
        for &q in &qs {
            for (i, kind) in kinds.into_iter().enumerate() {
                let (a, b) = (
                    serial.get(&pp, kind, q).unwrap(),
                    exact.get(&pp, kind, q).unwrap(),
                );
                identical &= a.to_bits() == b.to_bits();
                worst[i] = worst[i].max((a - b).abs());
                peak[i] = peak[i].max(b.abs());
            }
        }
        println!(
            "GD2 shells cutoff {cutoff} Ha: {} values, tables={tabulated}, max |table - quadrature| per kind {worst:?} (peaks {peak:?})",
            serial.len()
        );
        assert_eq!(identical, !tabulated);
        for i in 0..4 {
            assert!(worst[i] <= 1e-13 * peak[i], "{kinds:?} {worst:?}");
        }
        for ranks in [2, 3, 7] {
            let results = dft_mpi::simulate(ranks, |comm| {
                dft_engine::parallel::set_communicator(Some(Arc::new(comm.clone())));
                let v = dft_engine::forces::shell_values(&species, &db, volume, true)
                    .unwrap()
                    .unwrap();
                dft_engine::parallel::set_communicator(None);
                v
            });
            for v in results {
                for &q in &qs {
                    for kind in kinds {
                        assert_eq!(
                            v.get(&pp, kind, q).unwrap().to_bits(),
                            serial.get(&pp, kind, q).unwrap().to_bits()
                        );
                    }
                }
            }
        }
    }
}
