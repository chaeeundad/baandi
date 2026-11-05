//! Licensed files stay external. These are low-cutoff isolated-atom smoke runs,
//! NOT cutoff convergence, transferability, or matched external-code references.
use dft_engine::{collinear::*, *};
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use dft_physics::{data::*, spin::SpinDensity};

/// DFT_NC_DIR=/external/path cargo test --release -p dft-engine --test real_nc -- --ignored --nocapture --test-threads=1
#[test]
#[ignore = "external DFT_NC_DIR, release-only atom SCF smoke"]
fn isolated_atoms_pbe_smoke() {
    let dir = std::path::PathBuf::from(std::env::var_os("DFT_NC_DIR").expect("set DFT_NC_DIR"));
    let cell = Cell::from_user_rows_bohr([[10., 0., 0.], [0., 10., 0.], [0., 0., 10.]]).unwrap();
    let cutoff: f64 = std::env::var("DFT_NC_CUTOFF")
        .unwrap_or("20".into())
        .parse()
        .unwrap();
    assert!(
        cutoff == 20. || cutoff == 30.,
        "smoke cutoff must be 20 or 30 Ha"
    );
    for element in ["O", "Ni", "Al", "Mg", "Si"] {
        let start = std::time::Instant::now();
        let pp = dft_upf::from_file(dir.join(format!("{element}.upf"))).unwrap();
        let species = [SpeciesPotential {
            label: element,
            potential: &pp,
        }];
        let atoms = [Atom::new(element, [0.173, 0.227, 0.319])];
        let electrons = pp.z_valence() as usize;
        let options = ScfOptions {
            functional: Functional::Pbe,
            cutoff_hartree: cutoff,
            bands: Some(electrons.div_ceil(2) + 3),
            max_iterations: 60,
            ..Default::default()
        };
        // Explicit fixed sector permits odd-electron Al; no fractional occupation
        // or hidden exception to the legacy even-electron constructor.
        let sector = SpinSector {
            n_up: electrons.div_ceil(2),
            n_down: electrons / 2,
            axis_cartesian: [0., 0., 1.],
            initial_moment_seed: None,
        };
        let run = CollinearProblem::new(
            &cell,
            &atoms,
            &species,
            options,
            &[KPoint {
                fractional: [0.; 3],
                weight: 1.,
            }],
            sector,
        )
        .and_then(|p| p.run());
        match run {
            Ok(result) => {
                let SpinDensity::Collinear { channels, .. } = &result.density else {
                    panic!("wrong spin layout")
                };
                let count = (channels[0].sum() + channels[1].sum()) * cell.volume_bohr3()
                    / channels[0].len() as f64;
                assert!(result.energy.total_hartree().is_finite());
                assert!((count - pp.z_valence()).abs() <= 1e-10 * pp.z_valence());
                assert_eq!(result.energy.nlcc.is_some(), pp.core_density().is_some());
                println!(
                    "SMOKE {element}: cutoff_ha={cutoff} status={} energy_ha={:.16e} electrons={:.16e} iterations={} seconds={:.6}",
                    if result.converged {
                        "converged"
                    } else {
                        "not_converged_iteration_budget"
                    },
                    result.energy.total_hartree(),
                    count,
                    result.history.len(),
                    start.elapsed().as_secs_f64()
                );
            }
            Err(error) => {
                // An explicit solver/occupation/XC failure is smoke evidence, not
                // a passed SCF. No unavailable energy or count is fabricated.
                let reason = error.to_string();
                assert!(
                    matches!(&error, ScfError::NonInsulating { .. })
                        || reason.contains("NLCC combined spin density")
                        || reason
                            .contains("spin density residual exceeded configured growth limit"),
                    "unexpected smoke failure: {reason}"
                );
                let iteration = reason
                    .split("iteration ")
                    .nth(1)
                    .and_then(|s| s.split(':').next())
                    .unwrap_or("unavailable");
                println!(
                    "SMOKE {element}: cutoff_ha={cutoff} status=failed energy_ha=unavailable electrons_target={} iterations={iteration} seconds={:.6} reason={error}",
                    pp.z_valence(),
                    start.elapsed().as_secs_f64()
                );
            }
        }
    }
}

#[test]
#[ignore = "external DFT_NC_DIR; independent D-matrix assembly and radial inspection"]
fn real_nonlocal_and_radial_layout() {
    use dft_basis::PlaneWaveBasis;
    use dft_physics::hamiltonian::{NonlocalOperator, Species};
    let dir = std::path::PathBuf::from(std::env::var_os("DFT_NC_DIR").expect("set DFT_NC_DIR"));
    let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [0.7, 7., 0.], [0.3, 0.2, 9.]]).unwrap();
    let basis = PlaneWaveBasis::new(&cell, Hartree(1.)).unwrap();
    for element in ["O", "Ni", "Al", "Mg", "Si"] {
        let pp = dft_upf::from_file(dir.join(format!("{element}.upf"))).unwrap();
        let volume = cell.volume_bohr3();
        let op = NonlocalOperator::new(
            &basis,
            volume,
            &[Species {
                potential: &pp,
                positions_bohr: &[[0.37, 0.52, 0.91]],
            }],
        )
        .unwrap();
        let np = pp.projectors().len();
        for state in [0, basis.vectors().len() / 2, basis.vectors().len() - 1] {
            let mut c = vec![Complex64::default(); basis.vectors().len()];
            c[state] = Complex64::new(1., 0.);
            let mut out = c.clone();
            out.fill(Complex64::default());
            op.apply(&c, &mut out).unwrap();
            let q = basis.vectors()[state].norm_squared.sqrt();
            let radial: Vec<_> = (0..np)
                .map(|i| pp.projector_transform(i, q).unwrap())
                .collect();
            // Independent addition theorem sum_m |Y_lm|^2=(2l+1)/(4*pi).
            // One plane wave removes atomic phases, isolating full D coupling.
            let mut expected = 0.;
            for i in 0..np {
                for j in 0..np {
                    if pp.projectors()[i].l() == pp.projectors()[j].l() {
                        expected += radial[i]
                            * pp.dij_hartree()[i * np + j]
                            * radial[j]
                            * (2 * pp.projectors()[i].l() + 1) as f64
                            / (4. * std::f64::consts::PI * volume);
                    }
                }
            }
            let error = (out[state].re - expected).abs();
            assert!(error <= 1e-10 * expected.abs().max(1e-12));
            assert!(out[state].im.abs() < 1e-12);
            println!("NONLOCAL {element} q={q:e} energy_ha={expected:e} error_ha={error:e}");
        }
        let r = pp.mesh().r();
        let rab = pp.mesh().rab();
        let rho = &pp.raw().rhoatom.values;
        assert!(r.len() >= 936 && r.len() <= 1876);
        assert!(r.windows(2).all(|w| w[1] > w[0]));
        let integral = pp.mesh().integrate(rho).unwrap();
        // Independent weight-array assembly, including a final 3/8 panel.
        let n = r.len();
        let last = if n % 2 == 0 { n - 4 } else { n - 1 };
        let mut weights = vec![0.; n];
        for i in (0..last).step_by(2) {
            weights[i] += 1. / 3.;
            weights[i + 1] += 4. / 3.;
            weights[i + 2] += 1. / 3.;
        }
        if n % 2 == 0 {
            for (i, w) in [3. / 8., 9. / 8., 9. / 8., 3. / 8.].into_iter().enumerate() {
                weights[last + i] += w;
            }
        }
        let independent: f64 = (0..n).map(|i| weights[i] * rab[i] * rho[i]).sum();
        assert!((integral - independent).abs() <= 1e-10 * integral.abs());
        println!(
            "RADIAL {element} points={n} atomic_electrons={integral:.16e} target={} relative_charge_deficit={:.9e} independent_error={:.9e}; source finite radial table, not renormalized",
            pp.z_valence(),
            (integral - pp.z_valence()) / pp.z_valence(),
            integral - independent
        );
    }
}
