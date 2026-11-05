//! Synthetic NLCC only. Frozen caps: radial 1e-10 relative; directional XC
//! derivative 1e-7 Ha; rigid-translation XC energy 1e-7 Ha (same MF06 energy
//! scale, frozen before this test is run). No material or QE qualification.
use dft_engine::{Functional, ScfOptions, ScfProblem, SpeciesPotential};
use dft_geometry::{Atom, Cell};
use dft_physics::{
    data::*,
    spin::{SpinFunctional, SpinXc, SpinXcResponse},
};
use std::f64::consts::PI;
fn synthetic_text(core: bool, functional: &str) -> String {
    let n = 2001;
    let h = 0.006;
    let a = |f: &dyn Fn(f64) -> f64| {
        (0..n)
            .map(|i| format!("{:.17e}", f(i as f64 * h)))
            .collect::<Vec<_>>()
            .join(" ")
    };
    let nlcc = if core {
        format!(
            "<PP_NLCC>{}</PP_NLCC>",
            a(&|r| (0.5 / PI).powf(1.5) * (-0.5 * r * r).exp())
        )
    } else {
        String::new()
    };
    format!(
        r#"<UPF version="2.0.1"><PP_INFO>Synthetic Gaussian core</PP_INFO>
<PP_HEADER generated="synthetic" author="test" date="2026" comment="synthetic" element="Xx" pseudo_type="NC" relativistic="scalar" is_ultrasoft="F" is_paw="F" is_coulomb="F" has_so="F" has_wfc="F" has_gipaw="F" core_correction="{}" functional="{functional}" z_valence="2" total_psenergy="0" wfc_cutoff="0" rho_cutoff="0" l_max="0" l_max_rho="0" l_local="-1" mesh_size="{n}" number_of_wfc="0" number_of_proj="0"/>
<PP_MESH dx="{h}" mesh="{n}" xmin="0" rmax="12" zmesh="1"><PP_R>{}</PP_R><PP_RAB>{}</PP_RAB></PP_MESH>{nlcc}<PP_LOCAL>{}</PP_LOCAL><PP_NONLOCAL/><PP_RHOATOM>{}</PP_RHOATOM></UPF>"#,
        if core { "T" } else { "F" },
        a(&|r| r),
        a(&|_| h),
        a(&|r| if r == 0.0 {
            -4.0 / PI.sqrt()
        } else {
            -4.0 * erf(r / 2.0) / r
        }),
        a(&|r| 8.0 * PI * r * r * (0.5 / PI).powf(1.5) * (-0.5 * r * r).exp())
    )
}
fn synthetic(core: bool, functional: &str) -> dft_upf::Pseudopotential {
    dft_upf::parse_upf(&synthetic_text(core, functional)).unwrap()
}
fn trend(slopes: &[f64; 3], expected: f64) -> f64 {
    let cap = 1e-7;
    let d0 = (slopes[0] - slopes[1]).abs();
    let d1 = (slopes[1] - slopes[2]).abs();
    assert!(
        d1 <= 0.4 * d0
            || (slopes.iter().all(|v| (v - expected).abs() < cap) && d0.max(d1) < cap / 3.0)
    );
    let error = (slopes[2] - expected).abs();
    assert!(error < cap, "derivative error {error:e} Ha");
    error
}
#[test]
fn core_xc_localized_variations_and_spin_neutrality() {
    let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [1., 7., 0.], [0.5, 1., 9.]]).unwrap();
    let dims = [12, 14, 16];
    let n = Array3::from_shape_fn(dims, |(i, j, k)| {
        0.05 + 0.003 * ((i + 2 * j + 3 * k) as f64).sin()
    });
    let core = Array3::from_shape_fn(dims, |(i, _, _)| {
        0.02 + 0.01 * (2.0 * PI * i as f64 / 12.0).cos()
    });
    let mut scratch = Array3::zeros(dims);
    let mut fft = FftPlan3d::new(dims).unwrap();
    let weight = cell.volume_bohr3() / n.len() as f64;
    let index = [3, 5, 7];
    for f in [Functional::LdaPz, Functional::LdaPw, Functional::Pbe] {
        let x = f
            .evaluate_with_core(&n, Some(&core), &mut scratch, &cell, &mut fft)
            .unwrap();
        let expected = weight * 0.01 * x.potential[index];
        let mut slopes = [0.; 3];
        for (i, h) in [0.01, 0.005, 0.0025].into_iter().enumerate() {
            let mut plus = n.clone();
            plus[index] += h * 0.01;
            let mut minus = n.clone();
            minus[index] -= h * 0.01;
            let ep = f
                .evaluate_with_core(&plus, Some(&core), &mut scratch, &cell, &mut fft)
                .unwrap()
                .energy
                .0;
            let em = f
                .evaluate_with_core(&minus, Some(&core), &mut scratch, &cell, &mut fft)
                .unwrap()
                .energy
                .0;
            slopes[i] = (ep - em) / (2. * h);
        }
        println!(
            "{f:?} localized derivative error={} Ha",
            trend(&slopes, expected)
        );
        let integral = weight * n.iter().zip(&x.potential).map(|(n, v)| n * v).sum::<f64>();
        assert_eq!(integral.to_bits(), x.potential_integral.0.to_bits());
        let old = f.evaluate(&n, &cell, &mut fft).unwrap();
        let off = f
            .evaluate_with_core(&n, None, &mut scratch, &cell, &mut fft)
            .unwrap();
        assert_eq!(old.energy.0.to_bits(), off.energy.0.to_bits());
        assert!(
            old.potential
                .iter()
                .zip(&off.potential)
                .all(|(a, b)| a.to_bits() == b.to_bits())
        );
    }
    for f in [SpinFunctional::LdaPw, SpinFunctional::Pbe] {
        let mut xc = SpinXc::new(f, &cell, dims).unwrap();
        let mut out = SpinXcResponse::new(dims);
        let channels = [n.mapv(|n| 0.7 * n), n.mapv(|n| 0.3 * n)];
        xc.evaluate_with_core(&channels, Some(&core), &mut out)
            .unwrap();
        let expected = weight * 0.01 * out.potentials[0][index];
        let mut slopes = [0.; 3];
        for (i, h) in [0.01, 0.005, 0.0025].into_iter().enumerate() {
            let mut plus = channels.clone();
            plus[0][index] += h * 0.01;
            let mut minus = channels.clone();
            minus[0][index] -= h * 0.01;
            xc.evaluate_with_core(&plus, Some(&core), &mut out).unwrap();
            let ep = out.energy_ha;
            xc.evaluate_with_core(&minus, Some(&core), &mut out)
                .unwrap();
            slopes[i] = (ep - out.energy_ha) / (2. * h);
        }
        println!(
            "spin {f:?} localized derivative error={} Ha",
            trend(&slopes, expected)
        );
        xc.evaluate_with_core(&channels, Some(&core), &mut out)
            .unwrap();
        let energy = out.energy_ha;
        let v = out.potentials.clone();
        xc.evaluate_with_core(
            &[channels[1].clone(), channels[0].clone()],
            Some(&core),
            &mut out,
        )
        .unwrap();
        assert!((energy - out.energy_ha).abs() < 1e-10);
        for (s, potential) in v.iter().enumerate() {
            for (a, b) in potential.iter().zip(&out.potentials[1 - s]) {
                assert!((a - b).abs() < 1e-10);
            }
        }
    }
}
#[test]
fn gaussian_core_translation_and_engine_ledger() {
    let pp = synthetic(true, "SLA PZ");
    for rows in [
        [[8., 0., 0.], [0., 8., 0.], [0., 0., 8.]],
        [[8., 0., 0.], [1., 7., 0.], [0.5, 1., 9.]],
        [[6., 0., 0.], [0., 6., 0.], [0., 0., 18.]],
    ] {
        let cell = Cell::from_user_rows_bohr(rows).unwrap();
        let species = [SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let opts = ScfOptions {
            cutoff_hartree: 4.0,
            bands: Some(2),
            ..Default::default()
        };
        let mut energies = Vec::new();
        for pos in [
            [0., 0., 0.],
            [
                2.0_f64.sqrt() / 10.,
                3.0_f64.sqrt() / 10.,
                5.0_f64.sqrt() / 10.,
            ],
        ] {
            let p =
                ScfProblem::new(&cell, &[Atom::new("Xx", pos)], &species, opts.clone()).unwrap();
            let n = Array3::from_elem(p.grid_dimensions(), 0.05);
            let core = p.core_density().unwrap();
            assert!((core.sum() * cell.volume_bohr3() / core.len() as f64 - 1.).abs() < 1e-10);
            let mut scratch = Array3::zeros(p.grid_dimensions());
            let mut fft = FftPlan3d::new(p.grid_dimensions()).unwrap();
            let mut e = [0.; 5];
            e[..3].copy_from_slice(
                &[Functional::LdaPz, Functional::LdaPw, Functional::Pbe].map(|f| {
                    f.evaluate_with_core(&n, Some(core), &mut scratch, &cell, &mut fft)
                        .unwrap()
                        .energy
                        .0
                }),
            );
            for (i, f) in [SpinFunctional::LdaPw, SpinFunctional::Pbe]
                .into_iter()
                .enumerate()
            {
                let mut xc = SpinXc::new(f, &cell, p.grid_dimensions()).unwrap();
                let mut out = SpinXcResponse::new(p.grid_dimensions());
                xc.evaluate_with_core(
                    &[n.mapv(|n| 0.7 * n), n.mapv(|n| 0.3 * n)],
                    Some(core),
                    &mut out,
                )
                .unwrap();
                e[3 + i] = out.energy_ha;
            }
            energies.push(e);
        }
        for (f, (&a, &b)) in energies[0].iter().zip(&energies[1]).enumerate() {
            let delta = (a - b).abs();
            println!("translation functional={f} delta={delta:e} Ha");
            assert!(delta < 1e-7);
        }
    }
    let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [0., 8., 0.], [0., 0., 8.]]).unwrap();
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let mut p = ScfProblem::new(
        &cell,
        &[Atom::new("Xx", [0.; 3])],
        &species,
        ScfOptions {
            cutoff_hartree: 2.,
            bands: Some(2),
            ..Default::default()
        },
    )
    .unwrap();
    let step = p.step(&p.uniform_density(), None, 1e-8).unwrap();
    assert!(
        (step.density.sum() * cell.volume_bohr3() / step.density.len() as f64 - 2.).abs() < 1e-10
    );
    let core = step.energy.nlcc.as_ref().unwrap();
    assert!((core.core_electrons - 1.).abs() < 1e-10);
    assert!(core.core_potential_integral_ha < 0.);
    assert!((step.energy.total_hartree() - step.energy.band_form).abs() < 1e-7);
}

#[test]
fn five_projectors_coupled_d_against_analytic_plane_wave_energy() {
    use dft_basis::PlaneWaveBasis;
    use dft_physics::hamiltonian::{NonlocalOperator, Species};
    let ls = [0, 0, 1, 1, 2];
    let widths = [1.0_f64, 0.7, 0.9, 0.6, 0.8];
    let mut block = String::from("<PP_NONLOCAL>");
    for i in 0..5 {
        let v = (0..2001)
            .map(|j| {
                let r = j as f64 * 0.006;
                r.powi(ls[i] + 1) * (-widths[i] * r * r).exp()
            })
            .map(|x| format!("{x:.17e}"))
            .collect::<Vec<_>>()
            .join(" ");
        block.push_str(&format!("<PP_BETA.{} index=\"{}\" angular_momentum=\"{}\" cutoff_radius=\"12\" cutoff_radius_index=\"2001\">{v}</PP_BETA.{}>",i+1,i+1,ls[i],i+1));
    }
    let d: Vec<_> = (0..5)
        .flat_map(|i| {
            (0..5).map(move |j| {
                if i == j {
                    2. * (i + 1) as f64
                } else if ls[i] == ls[j] {
                    0.4
                } else {
                    0.
                }
            })
        })
        .collect();
    block.push_str(&format!(
        "<PP_DIJ>{}</PP_DIJ></PP_NONLOCAL>",
        d.iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(" ")
    ));
    let text = synthetic_text(false, "SLA PZ")
        .replace("number_of_proj=\"0\"", "number_of_proj=\"5\"")
        .replace("l_max=\"0\"", "l_max=\"2\"")
        .replace("<PP_NONLOCAL/>", &block);
    let pp = dft_upf::parse_upf(&text).unwrap();
    let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [1., 7., 0.], [0.5, 1., 9.]]).unwrap();
    let basis = PlaneWaveBasis::new(&cell, Hartree(1.)).unwrap();
    let op = NonlocalOperator::new(
        &basis,
        cell.volume_bohr3(),
        &[Species {
            potential: &pp,
            positions_bohr: &[[0.37, 0.52, 0.91]],
        }],
    )
    .unwrap();
    let mut worst = 0.0_f64;
    for state in [0, basis.vectors().len() / 2, basis.vectors().len() - 1] {
        let q = basis.vectors()[state].norm_squared.sqrt();
        let radial: Vec<_> = (0..5)
            .map(|i| {
                4. * PI * PI.sqrt() * q.powi(ls[i]) * (-q * q / (4. * widths[i])).exp()
                    / (2.0_f64.powi(ls[i] + 2) * widths[i].powf(ls[i] as f64 + 1.5))
            })
            .collect();
        let mut expected = 0.;
        for i in 0..5 {
            for j in 0..5 {
                if ls[i] == ls[j] {
                    expected +=
                        radial[i] * (d[5 * i + j] / 2.) * radial[j] * (2 * ls[i] + 1) as f64
                            / (4. * PI * cell.volume_bohr3());
                }
            }
        }
        let mut c = vec![Complex64::default(); basis.vectors().len()];
        c[state] = Complex64::new(1., 0.);
        let mut out = vec![Complex64::default(); c.len()];
        op.apply(&c, &mut out).unwrap();
        let error = (out[state].re - expected).abs() / expected.abs();
        worst = worst.max(error);
        assert!(error < 1e-10);
        assert!(out[state].im.abs() < 1e-12);
    }
    println!("five-projector analytic assembly max relative error={worst:e}");
}

#[test]
fn combined_density_rules_and_zero_spin_reduction() {
    let dims = [6, 6, 6];
    let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [0., 8., 0.], [0., 0., 8.]]).unwrap();
    let n = Array3::from_elem(dims, 0.05);
    let core = Array3::from_elem(dims, 0.02);
    let mut scratch = Array3::zeros(dims);
    let mut fft = FftPlan3d::new(dims).unwrap();
    for (f, s) in [
        (Functional::LdaPw, SpinFunctional::LdaPw),
        (Functional::Pbe, SpinFunctional::Pbe),
    ] {
        let unpolarized = f
            .evaluate_with_core(&n, Some(&core), &mut scratch, &cell, &mut fft)
            .unwrap();
        let mut xc = SpinXc::new(s, &cell, dims).unwrap();
        let mut out = SpinXcResponse::new(dims);
        xc.evaluate_with_core(
            &[n.mapv(|n| n / 2.), n.mapv(|n| n / 2.)],
            Some(&core),
            &mut out,
        )
        .unwrap();
        assert!((unpolarized.energy.0 - out.energy_ha).abs() < 1e-10);
        assert!(
            (unpolarized.potential_integral.0 - out.valence_potential_integral_ha).abs() < 1e-10
        );
        for v in &out.potentials {
            assert!(
                v.iter()
                    .zip(&unpolarized.potential)
                    .all(|(a, b)| (a - b).abs() < 1e-10)
            );
        }
    }
    let negative = Array3::from_elem(dims, -0.1);
    for f in [Functional::LdaPz, Functional::LdaPw, Functional::Pbe] {
        assert!(
            f.evaluate_with_core(&n, Some(&negative), &mut scratch, &cell, &mut fft)
                .is_err()
        );
    }
    // Negative Fourier core samples are not independently clipped: a feasible
    // combined density is accepted and the original scalar core stays intact.
    let small_negative = Array3::from_elem(dims, -0.01);
    for f in [Functional::LdaPz, Functional::LdaPw, Functional::Pbe] {
        assert!(
            f.evaluate_with_core(&n, Some(&small_negative), &mut scratch, &cell, &mut fft)
                .is_ok()
        );
    }
    assert!(small_negative.iter().all(|&x| x == -0.01));
}

/// cargo test --release -p dft-engine --test nlcc synthetic_core_scf -- --ignored --nocapture
#[test]
#[ignore = "release-only converged spin/unpolarized SCF reduction"]
fn synthetic_core_scf() {
    use dft_engine::collinear::{CollinearProblem, SpinSector};
    use dft_geometry::kpoints::KPoint;
    let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [0.4, 7., 0.], [0.2, 0.3, 9.]]).unwrap();
    for f in [Functional::LdaPw, Functional::Pbe] {
        let pp = synthetic(
            true,
            if f == Functional::Pbe {
                "PBE"
            } else {
                "SLA PW"
            },
        );
        let species = [SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let atoms = [Atom::new("Xx", [0.173, 0.227, 0.319])];
        let options = ScfOptions {
            functional: f,
            cutoff_hartree: 2.5,
            bands: Some(3),
            eigen_tolerance: 1e-11,
            density_tolerance: 1e-10,
            energy_tolerance: 1e-11,
            max_iterations: 120,
            ..Default::default()
        };
        let u = ScfProblem::new(&cell, &atoms, &species, options.clone())
            .unwrap()
            .run()
            .unwrap();
        let s = CollinearProblem::new(
            &cell,
            &atoms,
            &species,
            options,
            &[KPoint {
                fractional: [0.; 3],
                weight: 1.,
            }],
            SpinSector {
                n_up: 1,
                n_down: 1,
                axis_cartesian: [0., 0., 1.],
                initial_moment_seed: None,
            },
        )
        .unwrap()
        .run()
        .unwrap();
        assert!(s.converged);
        assert!(s.certificate.as_ref().unwrap().passed);
        let error = (u.energy.total_hartree() - s.energy.total_hartree()).abs();
        assert!(error <= 1e-10);
        let dft_physics::spin::SpinDensity::Collinear { channels, .. } = &s.density else {
            panic!("wrong density layout")
        };
        let distance = (u
            .density
            .iter()
            .zip(&channels[0])
            .zip(&channels[1])
            .map(|((&u, &a), &b)| (u - a - b).powi(2))
            .sum::<f64>()
            * cell.volume_bohr3()
            / u.density.len() as f64)
            .sqrt();
        assert!(distance <= 1e-9);
        assert!(u.energy.nlcc.as_ref().unwrap().core_potential_integral_ha < 0.);
        assert!((pp.mesh().integrate(&pp.raw().rhoatom.values).unwrap() - 2.).abs() <= 2e-10);
        println!(
            "synthetic core SCF {f:?}: unpolarized_energy_ha={:.16e} spin_energy_ha={:.16e} delta_ha={error:e} density_distance={distance:e} iterations={}/{}",
            u.energy.total_hartree(),
            s.energy.total_hartree(),
            u.history.len(),
            s.history.len()
        );
    }
}

/// Frozen before running: total-energy rigid-translation cap 1e-7 Ha, SCF
/// energy/density tightened independently. Two atoms, skew cell, no grid-aligned
/// displacement. This tests the entire ledger, not only point/grid XC.
#[test]
#[ignore = "release-only translated two-atom NLCC SCF oracle"]
fn rigid_translation_total_scf() {
    let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [0.7, 7., 0.], [0.3, 0.2, 9.]]).unwrap();
    let shift = [
        2.0_f64.sqrt() / 10.,
        3.0_f64.sqrt() / 10.,
        5.0_f64.sqrt() / 10.,
    ];
    for f in [Functional::LdaPz, Functional::LdaPw, Functional::Pbe] {
        let pp = synthetic(
            true,
            match f {
                Functional::LdaPz => "SLA PZ",
                Functional::LdaPw => "SLA PW",
                Functional::Pbe => "PBE",
            },
        );
        let species = [SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let options = ScfOptions {
            functional: f,
            cutoff_hartree: 6.,
            bands: Some(4),
            eigen_tolerance: 1e-11,
            density_tolerance: 1e-10,
            energy_tolerance: 1e-11,
            max_iterations: 150,
            ..Default::default()
        };
        let base_atoms = [
            Atom::new("Xx", [0.1, 0.2, 0.3]),
            Atom::new("Xx", [0.37, 0.43, 0.49]),
        ];
        let dims = ScfProblem::new(&cell, &base_atoms, &species, options.clone())
            .unwrap()
            .grid_dimensions();
        for refinement in [1, 2, 3] {
            let mut options = options.clone();
            options.grid_override = Some(dims.map(|n| n * refinement));
            let mut energies = [0.; 2];
            for (i, displacement) in [[0.; 3], shift].into_iter().enumerate() {
                let atoms: [[f64; 3]; 2] = [[0.1, 0.2, 0.3], [0.37, 0.43, 0.49]];
                let atoms: Vec<_> = atoms
                    .into_iter()
                    .map(|a| Atom::new("Xx", std::array::from_fn(|j| a[j] + displacement[j])))
                    .collect();
                let r = ScfProblem::new(&cell, &atoms, &species, options.clone())
                    .unwrap()
                    .run()
                    .unwrap();
                energies[i] = r.energy.total_hartree();
                assert!((r.energy.nlcc.as_ref().unwrap().core_electrons - 2.).abs() < 2e-10);
                assert!((r.energy.total_hartree() - r.energy.band_form).abs() < 1e-7);
            }
            let delta = (energies[0] - energies[1]).abs();
            println!(
                "total SCF translation {f:?}: refinement={refinement} grid={:?} energies_ha={energies:?} delta_ha={delta:e}",
                options.grid_override.unwrap()
            );
            // Coarse PBE quadrature is explicitly retained as failed evidence.
            // Acceptance is on the independently refined grid, same basis/cap.
            if refinement >= 2 {
                assert!(delta <= 1e-7);
            }
        }
    }
}
