//! Independent strain-energy differences; no SCF or cutoff-convergence claim.
use dft_basis::PlaneWaveBasis;
use dft_geometry::Cell;
use dft_physics::{ewald, kinetic};
use dft_units::Hartree;
use num_complex::Complex64;

fn strained(cell: &Cell, a: usize, b: usize, h: f64) -> Cell {
    let columns = cell.lattice_columns_bohr();
    let mut next = columns;
    for j in 0..3 {
        next[a][j] += h * columns[b][j];
        if a != b {
            next[b][j] += h * columns[a][j];
        }
    }
    Cell::from_user_rows_bohr(std::array::from_fn(|i| std::array::from_fn(|j| next[j][i]))).unwrap()
}

fn check_differences(cell: &Cell, tensor: [[f64; 3]; 3], energy: impl Fn(&Cell) -> f64) {
    let mut maximum = [0.0_f64; 3];
    for (a, row) in tensor.iter().enumerate() {
        for (b, value) in row.iter().enumerate().skip(a) {
            assert!((value - tensor[b][a]).abs() < 1e-10);
            let mut errors = [0.0; 3];
            for (i, h) in [0.002, 0.001, 0.0005].into_iter().enumerate() {
                let fd = (energy(&strained(cell, a, b, h)) - energy(&strained(cell, a, b, -h)))
                    / (2.0 * h * cell.volume_bohr3() * if a == b { 1.0 } else { 2.0 });
                errors[i] = (fd - value).abs();
                maximum[i] = maximum[i].max(errors[i]);
            }
            assert!(errors[2] < 1e-7, "({a},{b}) errors {errors:?} Ha/Bohr^3");
            for pair in errors.windows(2) {
                assert!(
                    pair[1] < 0.3 * pair[0] || pair[1] < 1e-11,
                    "second order or roundoff plateau: {errors:?}"
                );
            }
        }
    }
    println!("strain steps [0.002,0.001,0.0005], max errors {maximum:?} Ha/Bohr^3");
}

#[test]
fn frozen_coefficient_kinetic_strain_metric() {
    for rows in [
        [[8., 0., 0.], [0., 8., 0.], [0., 0., 8.]],
        [[8., 0.3, 0.1], [1.2, 7., 0.6], [0.4, 0.8, 9.]],
        [[7., 0., 0.], [0.2, 8., 0.], [0., 0., 22.]],
    ] {
        let cell = Cell::from_user_rows_bohr(rows).unwrap();
        let basis = PlaneWaveBasis::at_k(&cell, Hartree(0.8), [0.13, -0.21, 0.07]).unwrap();
        let mut c: Vec<_> = (0..basis.vectors().len())
            .map(|i| Complex64::new((i as f64 + 0.7).sin(), (i as f64 * 0.3).cos()))
            .collect();
        let norm = c.iter().map(|v| v.norm_sqr()).sum::<f64>().sqrt();
        for v in &mut c {
            *v /= norm;
        }
        let coefficients = vec![c];
        let occupations = [1.3];
        let tensor =
            kinetic::stress(&basis, &coefficients, &occupations, cell.volume_bohr3()).unwrap();
        let energy = kinetic::energy(&basis, &coefficients, &occupations)
            .unwrap()
            .0;
        let trace: f64 = (0..3).map(|a| tensor[a][a]).sum();
        assert!((trace + 2.0 * energy / cell.volume_bohr3()).abs() < 1e-14);
        assert!(-trace / 3.0 > 0.0);
        // Explicit integer list from the center: do not rebuild a cutoff sphere.
        check_differences(&cell, tensor, |strained| {
            let b = strained.reciprocal_columns_inverse_bohr();
            basis
                .vectors()
                .iter()
                .zip(&coefficients[0])
                .map(|(g, c)| {
                    let q: [f64; 3] = std::array::from_fn(|i| {
                        (0..3)
                            .map(|j| b[i][j] * (f64::from(g.miller[j]) + basis.k_fractional()[j]))
                            .sum()
                    });
                    0.5 * occupations[0] * c.norm_sqr() * q.iter().map(|x| x * x).sum::<f64>()
                })
                .sum()
        });
        assert!(kinetic::stress(&basis, &coefficients, &occupations, 0.0).is_err());
        assert!(kinetic::stress(&basis, &coefficients, &[f64::NAN], cell.volume_bohr3()).is_err());
    }
}

#[test]
#[ignore = "Ewald strain sweep; run cargo test --release -p dft-physics --test stress_terms -- --ignored --nocapture"]
fn ewald_strain_includes_background_and_eta_cancellation() {
    for rows in [
        [[8., 0., 0.], [0., 8., 0.], [0., 0., 8.]],
        [[8., 0.3, 0.1], [1.2, 7., 0.6], [0.4, 0.8, 9.]],
        [[7., 0., 0.], [0.2, 8., 0.], [0., 0., 22.]],
    ] {
        let cell = Cell::from_user_rows_bohr(rows).unwrap();
        let ions = [
            ewald::Ion {
                charge: 1.3,
                fractional_position: [0.1, 0.2, 0.3],
            },
            ewald::Ion {
                charge: 0.7,
                fractional_position: [0.4, 0.65, 0.82],
            },
        ];
        let options = ewald::EwaldOptions::default();
        let tensor = ewald::stress(&cell, &ions, options).unwrap();
        check_differences(&cell, tensor, |c| {
            ewald::energy(c, &ions, options).unwrap().energy.0
        });
        let fixed = ewald::EwaldOptions {
            eta: Some(0.31),
            ..options
        };
        let other = ewald::stress(&cell, &ions, fixed).unwrap();
        for a in 0..3 {
            for b in 0..3 {
                assert!((tensor[a][b] - other[a][b]).abs() < 1e-11);
            }
        }
        check_differences(&cell, other, |c| {
            ewald::energy(c, &ions, fixed).unwrap().energy.0
        });
    }
    let cubic = Cell::from_user_rows_bohr([[8., 0., 0.], [0., 8., 0.], [0., 0., 8.]]).unwrap();
    let ions = [ewald::Ion {
        charge: 1.0,
        fractional_position: [0.; 3],
    }];
    let stress = ewald::stress(&cubic, &ions, ewald::EwaldOptions::default()).unwrap();
    let energy = ewald::energy(&cubic, &ions, ewald::EwaldOptions::default())
        .unwrap()
        .energy
        .0;
    for (a, row) in stress.iter().enumerate() {
        for (b, value) in row.iter().enumerate() {
            let expected = if a == b {
                -energy / (3.0 * cubic.volume_bohr3())
            } else {
                0.0
            };
            assert!((value - expected).abs() < 1e-9);
        }
    }
}

#[allow(dead_code)]
#[path = "../../dft-engine/tests/common/mod.rs"]
mod common;
mod stress_support;

#[test]
#[ignore = "release radial/metric sweep"]
fn isolated_local_core_hartree_nonlocal() {
    use dft_physics::{
        data::*,
        hamiltonian::{NonlocalOperator, Species, ionic_local},
        stress,
    };
    let cell = Cell::from_user_rows_bohr([[8., 0.2, 0.], [1., 7., 0.3], [0.5, 0.7, 9.]]).unwrap();
    let basis = PlaneWaveBasis::at_k(&cell, Hartree(0.8), [0.12, 0.04, -0.1]).unwrap();
    let pp = stress_support::potential(1., 1.5, true, 401);
    let position = cell.fractional_to_cartesian_bohr([0.17, 0.28, 0.39]);
    let species = [Species {
        potential: &pp,
        positions_bohr: std::slice::from_ref(&position),
    }];
    let mut orbitals = ComplexMatrix::zeros(basis.vectors().len(), 2).unwrap();
    for band in 0..2 {
        for g in 0..orbitals.rows() {
            orbitals[(g, band)] =
                Complex64::new((g as f64 * 0.7 + band as f64).sin(), (g as f64 * 0.3).cos());
        }
        let norm = orbitals
            .column(band)
            .iter()
            .map(|c| c.norm_sqr())
            .sum::<f64>()
            .sqrt();
        for g in 0..orbitals.rows() {
            orbitals[(g, band)] /= norm;
        }
    }
    let f = [0.8, 0.3];
    let tensor = stress::nonlocal(&basis, cell.volume_bohr3(), &species, &orbitals, &f).unwrap();
    check_differences(&cell, tensor, |c| {
        let b = basis.in_cell(c).unwrap();
        let position = c.fractional_to_cartesian_bohr([0.17, 0.28, 0.39]);
        let species = [Species {
            potential: &pp,
            positions_bohr: std::slice::from_ref(&position),
        }];
        let op = NonlocalOperator::new(&b, c.volume_bohr3(), &species).unwrap();
        let mut y = vec![Complex64::default(); orbitals.rows()];
        (0..2)
            .map(|band| {
                op.apply(orbitals.column(band), &mut y).unwrap();
                f[band]
                    * orbitals
                        .column(band)
                        .iter()
                        .zip(&y)
                        .map(|(a, b)| (a.conj() * b).re)
                        .sum::<f64>()
            })
            .sum()
    });
    let db = PlaneWaveBasis::new(&cell, Hartree(1.)).unwrap();
    let field: Vec<_> = db
        .vectors()
        .iter()
        .map(|g| Complex64::new(0.02 * (-g.norm_squared).exp(), 0.))
        .collect();
    let local = stress::local_or_core(&db, cell.volume_bohr3(), &species, &field, false).unwrap();
    check_differences(&cell, local, |c| {
        let b = db.in_cell(c).unwrap();
        let position = c.fractional_to_cartesian_bohr([0.17, 0.28, 0.39]);
        let species = [Species {
            potential: &pp,
            positions_bohr: std::slice::from_ref(&position),
        }];
        cell.volume_bohr3()
            * field
                .iter()
                .zip(ionic_local(&b, c.volume_bohr3(), &species).unwrap())
                .map(|(a, b)| (a.conj() * b).re)
                .sum::<f64>()
    });
    let dims = [12; 3];
    let n = Array3::from_shape_fn(dims, |(i, j, k)| {
        0.03 + 0.01 * (std::f64::consts::TAU * (i + 2 * j + k) as f64 / 12.).cos()
    });
    let mut fft = FftPlan3d::new(dims).unwrap();
    let h = stress::hartree(&n, &cell, &mut fft).unwrap();
    check_differences(&cell, h, |c| {
        dft_physics::hartree::evaluate(
            &(&n * (cell.volume_bohr3() / c.volume_bohr3())),
            c,
            &mut FftPlan3d::new(dims).unwrap(),
        )
        .unwrap()
        .energy
        .0
    });
}

use dft_physics::data::erf;
use std::f64::consts::PI;
fn synthetic_text(core: bool, functional: &str) -> String {
    let n = 1001;
    let h = 0.012;
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

#[test]
#[ignore = "release core and XC metric sweeps"]
fn core_and_xc_metric() {
    use dft_physics::{
        data::*,
        hamiltonian::{Species, structure_factor},
        spin::{SpinFunctional, SpinXc, SpinXcResponse},
        stress,
    };
    let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [1., 7., 0.], [0.5, 1., 9.]]).unwrap();
    let pp = synthetic(true, "PBE");
    let db = PlaneWaveBasis::new(&cell, Hartree(2.)).unwrap();
    let position = cell.fractional_to_cartesian_bohr([0.2, 0.3, 0.4]);
    let species = [Species {
        potential: &pp,
        positions_bohr: std::slice::from_ref(&position),
    }];
    let field: Vec<_> = db
        .vectors()
        .iter()
        .map(|g| Complex64::new((-g.norm_squared).exp(), 0.))
        .collect();
    let tensor = stress::local_or_core(&db, cell.volume_bohr3(), &species, &field, true).unwrap();
    check_differences(&cell, tensor, |c| {
        let b = db.in_cell(c).unwrap();
        let p = c.fractional_to_cartesian_bohr([0.2, 0.3, 0.4]);
        cell.volume_bohr3() / c.volume_bohr3()
            * b.vectors()
                .iter()
                .zip(&field)
                .map(|(g, f)| {
                    let radial = dft_upf::radial::spherical_density_transform(
                        pp.mesh(),
                        pp.core_density().unwrap(),
                        g.norm_squared.sqrt(),
                    )
                    .unwrap();
                    (f.conj() * radial * structure_factor(g.cartesian, &[p]).unwrap()).re
                })
                .sum::<f64>()
    });
    let dims = [12; 3];
    let n = Array3::from_shape_fn(dims, |(i, j, k)| {
        0.03 + 0.01 * (std::f64::consts::TAU * (i + j + 2 * k) as f64 / 12.).cos()
    });
    for kind in [0, 1, 2] {
        let eval = |n: &Array3<f64>, c: &Cell| match kind {
            0 => dft_physics::xc::evaluate(n, c.volume_bohr3()).unwrap(),
            1 => dft_physics::spin::unpolarized_pw(n, c.volume_bohr3()).unwrap(),
            _ => dft_physics::gga::evaluate(n, c, &mut FftPlan3d::new(dims).unwrap()).unwrap(),
        };
        let xc = eval(&n, &cell);
        let mut tensor = if kind == 2 {
            dft_physics::gga::stress_metric(&n, &cell, &mut FftPlan3d::new(dims).unwrap()).unwrap()
        } else {
            [[0.; 3]; 3]
        };
        for (a, row) in tensor.iter_mut().enumerate() {
            row[a] += (xc.energy.0 - xc.potential_integral.0) / cell.volume_bohr3();
        }
        check_differences(&cell, tensor, |c| {
            eval(&(&n * (cell.volume_bohr3() / c.volume_bohr3())), c)
                .energy
                .0
        });
    }
    for kind in [SpinFunctional::LdaPw, SpinFunctional::Pbe] {
        let densities = [&n * 0.7, &n * 0.3];
        let mut xc = SpinXc::new(kind, &cell, dims).unwrap();
        let mut out = SpinXcResponse::new(dims);
        xc.evaluate(&densities, &mut out).unwrap();
        let mut tensor = xc.stress_metric().unwrap();
        for (a, row) in tensor.iter_mut().enumerate() {
            row[a] += (out.energy_ha - out.valence_potential_integral_ha) / cell.volume_bohr3();
        }
        check_differences(&cell, tensor, |c| {
            let mut xc = SpinXc::new(kind, c, dims).unwrap();
            let mut out = SpinXcResponse::new(dims);
            let scale = cell.volume_bohr3() / c.volume_bohr3();
            xc.evaluate(&[&densities[0] * scale, &densities[1] * scale], &mut out)
                .unwrap();
            out.energy_ha
        });
    }
}
