//! Independent neutral Gaussian electrostatics. No constants fitted to energies.
mod common;
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_engine::energy::local_energy;
use dft_geometry::Cell;
use dft_physics::{
    data::*,
    ewald::{self, EwaldOptions, Ion},
    hamiltonian::{Species, ionic_local},
    hartree,
};
use std::f64::consts::{PI, TAU};

// Exact C = Gaussian self energies - short-range image/pair deficits.
// For the fixtures, summing to d=12 sqrt(2) rc omits <1e-60 Ha.
// Dual-lattice bounds include all translations inside that sphere.
fn correction(cell: &Cell, ions: &[Ion], rc: f64) -> (f64, f64) {
    let self_energy =
        ions.iter().map(|i| i.charge * i.charge).sum::<f64>() / ((2.0 * PI).sqrt() * rc);
    let radius = 12.0 * 2.0_f64.sqrt() * rc;
    let b = cell.reciprocal_columns_inverse_bohr();
    let bounds: [i32; 3] = std::array::from_fn(|j| {
        ((0..3).map(|i| b[i][j] * b[i][j]).sum::<f64>().sqrt() * radius / TAU + 1.0).ceil() as i32
    });
    let mut pair = 0.0;
    for (i, a) in ions.iter().enumerate() {
        for (j, b) in ions.iter().enumerate() {
            for x in -bounds[0]..=bounds[0] {
                for y in -bounds[1]..=bounds[1] {
                    for z in -bounds[2]..=bounds[2] {
                        if i == j && [x, y, z] == [0; 3] {
                            continue;
                        }
                        let m = [x, y, z];
                        let s = std::array::from_fn(|k| {
                            a.fractional_position[k].rem_euclid(1.0)
                                - b.fractional_position[k].rem_euclid(1.0)
                                + m[k] as f64
                        });
                        let d = cell
                            .fractional_to_cartesian_bohr(s)
                            .iter()
                            .map(|r| r * r)
                            .sum::<f64>()
                            .sqrt();
                        if d <= radius {
                            pair += 0.5 * a.charge * b.charge * erfc(d / (2.0_f64.sqrt() * rc)) / d;
                        }
                    }
                }
            }
        }
    }
    (self_energy - pair, pair)
}
#[test]
#[ignore = "48 independent electrostatic cases and refined reciprocal sphere; run in release"]
fn neutral_gaussian_cross_check() {
    // Also exercise the common synthetic geometry helper in this test binary.
    assert_eq!(common::atoms().len(), 2);
    let cases = [
        (
            common::cell(),
            vec![
                Ion {
                    charge: 4.0,
                    fractional_position: [0.0; 3],
                },
                Ion {
                    charge: 4.0,
                    fractional_position: [0.25; 3],
                },
            ],
        ),
        (
            Cell::from_user_rows_bohr([[7.0, 0.3, 0.1], [1.1, 6.8, 0.4], [0.6, 1.3, 7.2]]).unwrap(),
            vec![
                Ion {
                    charge: 2.0,
                    fractional_position: [0.03, 0.11, 0.08],
                },
                Ion {
                    charge: 3.0,
                    fractional_position: [0.39, 0.48, 0.52],
                },
                Ion {
                    charge: 1.0,
                    fractional_position: [0.71, 0.12, 0.69],
                },
            ],
        ),
    ];
    let mut max_relative = 0.0_f64;
    for (ci, (cell, original)) in cases.iter().enumerate() {
        let volume = cell.volume_bohr3();
        let q = original.iter().map(|i| i.charge).sum::<f64>();
        // Analytic Gaussian reciprocal tails at G>=12 are exponentially below
        // 1e-12 Ha for rc>=0.8 and electron widths>=1.0. A second cutoff checks it.
        for cutoff in [72.0, 90.0] {
            let basis = PlaneWaveBasis::new(cell, Hartree(cutoff)).unwrap();
            let grid = FftGridSpec::new(&basis, &basis, None).unwrap();
            let mut fft = FftPlan3d::new(grid.dims()).unwrap();
            for rc in [0.8, 1.0, 1.2] {
                let potentials: Vec<_> = original
                    .iter()
                    .map(|i| common::potential(i.charge, rc, false, 2001))
                    .collect();
                let mut unshifted = [0.0; 2];
                for translation in [false, true] {
                    let shift = if translation {
                        [0.183, -0.227, 0.319]
                    } else {
                        [0.0; 3]
                    };
                    let ions: Vec<_> = original
                        .iter()
                        .map(|i| Ion {
                            charge: i.charge,
                            fractional_position: std::array::from_fn(|a| {
                                i.fractional_position[a] + shift[a]
                            }),
                        })
                        .collect();
                    let positions: Vec<_> = ions
                        .iter()
                        .map(|i| [cell.fractional_to_cartesian_bohr(i.fractional_position)])
                        .collect();
                    let species: Vec<_> = potentials
                        .iter()
                        .zip(&positions)
                        .map(|(potential, positions_bohr)| Species {
                            potential,
                            positions_bohr,
                        })
                        .collect();
                    let v = ionic_local(&basis, volume, &species).unwrap();
                    let vg = grid.scatter(&v).unwrap();
                    let ew = ewald::energy(
                        cell,
                        &ions,
                        EwaldOptions {
                            tolerance: 1e-15,
                            ..Default::default()
                        },
                    )
                    .unwrap();
                    let (c, pair) = correction(cell, &ions, rc);
                    for (shape, offset) in [[0.07, 0.03, -0.09], [-0.15, 0.17, 0.04]]
                        .into_iter()
                        .enumerate()
                    {
                        let width = 1.0 + shape as f64 * 0.4;
                        let mut density = Vec::new();
                        let mut direct = 0.0;
                        for g in basis.vectors() {
                            let phase = |s: [f64; 3]| {
                                Complex64::from_polar(
                                    1.0,
                                    -TAU * (0..3).map(|j| g.miller[j] as f64 * s[j]).sum::<f64>(),
                                )
                            };
                            let rho: Complex64 = ions
                                .iter()
                                .map(|i| {
                                    i.charge / volume
                                        * (-g.norm_squared * rc * rc / 4.0).exp()
                                        * phase(i.fractional_position)
                                })
                                .sum();
                            let n = if g.miller == [0; 3] {
                                Complex64::new(q / volume, 0.0)
                            } else {
                                ions.iter()
                                    .map(|i| {
                                        i.charge / volume
                                            * (-g.norm_squared * width * width / 4.0).exp()
                                            * phase(std::array::from_fn(|j| {
                                                i.fractional_position[j] + offset[j]
                                            }))
                                    })
                                    .sum()
                            };
                            if g.norm_squared > 0.0 {
                                direct +=
                                    0.5 * volume * 4.0 * PI / g.norm_squared * (rho - n).norm_sqr();
                            }
                            density.push(n);
                        }
                        let mut real = Array3::zeros(grid.dims());
                        fft.to_real(&grid.scatter(&density).unwrap(), &mut real)
                            .unwrap();
                        assert!(real.iter().all(|v| v.re > 0.0 && v.im.abs() < 1e-12));
                        let n = real.mapv(|v| v.re);
                        let h = hartree::evaluate(&n, cell, &mut fft).unwrap().energy.0;
                        let local = local_energy(&n, &vg, volume, &mut fft).unwrap();
                        let engine = h + local + ew.energy.0;
                        let rel = (engine + c - direct).abs() / direct.abs();
                        max_relative = max_relative.max(rel);
                        let zero = q * v[0].re;
                        println!(
                            "Gaussian: cell={ci}, cutoff={cutoff} Ha, rc={rc} Bohr, shape={shape}, translated={translation}, H={h:.12e}, loc={local:.12e}, loc0={zero:.12e}, Ewald={:.12e}, C={c:.12e}, pair={pair:.6e}, direct={direct:.12e} Ha, relative={rel:.6e}",
                            ew.energy.0
                        );
                        assert!(rel < 1e-9, "electrostatic convention mismatch: {rel:e}");
                        assert!((n.sum() * volume / n.len() as f64 - q).abs() / q < 1e-12);
                        if translation {
                            assert!((direct - unshifted[shape]).abs() < 1e-10);
                        } else {
                            unshifted[shape] = direct;
                        }
                    }
                }
            }
        }
    }
    println!("Gaussian cross-check maximum relative error={max_relative:.12e}");
}
