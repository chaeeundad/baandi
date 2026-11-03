//! Independent synthetic references. No real UPF data or QE reference numbers.
#[path = "../../dft-upf/tests/common/mod.rs"]
mod synthetic;
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_geometry::Cell;
use dft_numerics::{
    Fft3d, FftPlan3d,
    linalg::{ComplexMatrix, hermitian_eigh},
};
use dft_physics::hamiltonian::*;
use dft_solvers::{Davidson, DenseEigensolver, DiagonalPreconditioner, LinearOperator};
use dft_units::{Angstrom, Bohr, Hartree};
use dft_upf::{Pseudopotential, parse_upf};
use ndarray::Array3;
use num_complex::Complex64 as C;
use std::{
    collections::BTreeMap,
    f64::consts::{PI, TAU},
};
fn setup() -> (Cell, PlaneWaveBasis, Pseudopotential, Vec<[f64; 3]>) {
    let a = Bohr::from(Angstrom(5.43)).0;
    let cell = Cell::from_user_rows_bohr([
        [0., a / 2., a / 2.],
        [a / 2., 0., a / 2.],
        [a / 2., a / 2., 0.],
    ])
    .unwrap();
    let basis = PlaneWaveBasis::new(&cell, Hartree(3.0)).unwrap();
    let pp = parse_upf(&synthetic::fixture_with_points(801)).unwrap();
    (cell, basis, pp, vec![[0.; 3], [a / 4.; 3]])
}
fn random(n: usize, mut seed: u64) -> Vec<C> {
    let mut next = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        (seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5
    };
    (0..n).map(|_| C::new(next(), next())).collect()
}
fn maxdiff(x: &[C], y: &[C]) -> f64 {
    x.iter()
        .zip(y)
        .map(|(x, y)| (x - y).norm())
        .fold(0., f64::max)
}
fn dot(x: &[C], y: &[C]) -> C {
    x.iter().zip(y).map(|(x, y)| x.conj() * y).sum()
}
fn mv(a: &ComplexMatrix, x: &[C]) -> Vec<C> {
    (0..a.rows())
        .map(|i| (0..a.cols()).map(|j| a[(i, j)] * x[j]).sum())
        .collect()
}
fn apply(op: &dyn LinearOperator, x: &[C]) -> Vec<C> {
    let mut y = vec![C::default(); x.len()];
    op.apply(x, &mut y).unwrap();
    y
}
fn nl(b: &PlaneWaveBasis, cell: &Cell, pp: &Pseudopotential, pos: &[[f64; 3]]) -> NonlocalOperator {
    NonlocalOperator::new(
        b,
        cell.volume_bohr3(),
        &[Species {
            potential: pp,
            positions_bohr: pos,
        }],
    )
    .unwrap()
}
// Independent angular/phase reference specialized to this synthetic s,s,p fixture.
// No access to the production projector tables, coupling expansion or FFT path.
fn reference_columns(
    b: &PlaneWaveBasis,
    cell: &Cell,
    pp: &Pseudopotential,
    pos: &[[f64; 3]],
) -> Vec<Vec<Vec<C>>> {
    pos.iter()
        .map(|r| {
            (0..5)
                .map(|c| {
                    b.vectors()
                        .iter()
                        .map(|g| {
                            let (i, l) = if c < 2 { (c, 0) } else { (2, 1) };
                            let q = g.norm_squared.sqrt();
                            let angular = if l == 0 {
                                1. / (4. * PI).sqrt()
                            } else if q == 0. {
                                0.
                            } else {
                                (3. / (4. * PI)).sqrt() * g.cartesian[[1, 2, 0][c - 2]] / q
                            };
                            let angle: f64 = (0..3).map(|k| g.cartesian[k] * r[k]).sum();
                            let phase = C::new(angle.cos(), -angle.sin());
                            phase
                                * if l == 0 {
                                    C::new(1., 0.)
                                } else {
                                    C::new(0., -1.)
                                }
                                * angular
                                * pp.projector_transform(i, q).unwrap()
                                / cell.volume_bohr3().sqrt()
                        })
                        .collect()
                })
                .collect()
        })
        .collect()
}
fn reference_nl(
    b: &PlaneWaveBasis,
    cell: &Cell,
    pp: &Pseudopotential,
    pos: &[[f64; 3]],
) -> ComplexMatrix {
    let columns = reference_columns(b, cell, pp, pos);
    let n = b.vectors().len();
    let mut a = ComplexMatrix::zeros(n, n).unwrap();
    for row in 0..n {
        for col in 0..n {
            for p in &columns {
                for i in 0..2 {
                    for j in 0..2 {
                        a[(row, col)] += p[i][row] * pp.dij_hartree()[i * 3 + j] * p[j][col].conj();
                    }
                }
                for channel in p.iter().skip(2) {
                    a[(row, col)] += channel[row] * pp.dij_hartree()[8] * channel[col].conj();
                }
            }
        }
    }
    a
}
fn sampled(dims: [usize; 3], coefficients: &BTreeMap<[i32; 3], C>) -> Array3<f64> {
    Array3::from_shape_fn(dims, |(i, j, k)| {
        let f = [
            i as f64 / dims[0] as f64,
            j as f64 / dims[1] as f64,
            k as f64 / dims[2] as f64,
        ];
        coefficients
            .iter()
            .map(|(m, c)| {
                *c * C::from_polar(1., TAU * (0..3).map(|a| m[a] as f64 * f[a]).sum::<f64>())
            })
            .sum::<C>()
            .re
    })
}
fn local_coefficients() -> BTreeMap<[i32; 3], C> {
    BTreeMap::from([
        ([0, 0, 0], C::new(0.17, 0.)),
        ([1, 0, 0], C::new(0.07, 0.03)),
        ([-1, 0, 0], C::new(0.07, -0.03)),
        ([0, 1, -1], C::new(-0.04, 0.02)),
        ([0, -1, 1], C::new(-0.04, -0.02)),
    ])
}
#[test]
fn ionic_phases_and_zero_mode() {
    let (cell, b, pp, pos) = setup();
    let density = b.density_basis(&cell).unwrap();
    let grid = FftGridSpec::new(&density, &density, None).unwrap();
    let delta = [0.31, -0.27, 0.18];
    let shifted: Vec<_> = pos
        .iter()
        .map(|r| std::array::from_fn(|k| r[k] + delta[k]))
        .collect();
    let v = ionic_local(
        &density,
        cell.volume_bohr3(),
        &[Species {
            potential: &pp,
            positions_bohr: &pos,
        }],
    )
    .unwrap();
    let vs = ionic_local(
        &density,
        cell.volume_bohr3(),
        &[Species {
            potential: &pp,
            positions_bohr: &shifted,
        }],
    )
    .unwrap();
    let transform = pp.local_transform(cell.volume_bohr3(), None).unwrap();
    let mut translation: f64 = 0.;
    let mut unit: f64 = 0.;
    let mut direct: f64 = 0.;
    for (i, g) in density.vectors().iter().enumerate() {
        let phase = C::from_polar(1., -(0..3).map(|k| g.cartesian[k] * delta[k]).sum::<f64>());
        translation = translation.max((vs[i] - v[i] * phase).norm());
        unit = unit.max((structure_factor(g.cartesian, &pos[1..]).unwrap().norm() - 1.).abs());
        let radial = if i == 0 {
            transform.alpha_z().unwrap()
        } else {
            transform.at_g(g.norm_squared.sqrt()).unwrap()
        };
        let sf: C = pos
            .iter()
            .map(|r| {
                let t: f64 = (0..3).map(|k| g.cartesian[k] * r[k]).sum();
                C::new(t.cos(), -t.sin())
            })
            .sum();
        direct = direct.max((v[i] - radial * sf).norm());
    }
    assert_eq!(v[0], C::new(2. * transform.alpha_z().unwrap(), 0.));
    let mut real = Array3::zeros(grid.dims());
    FftPlan3d::new(grid.dims())
        .unwrap()
        .to_real(&grid.scatter(&v).unwrap(), &mut real)
        .unwrap();
    let imaginary = real.iter().map(|z| z.im.abs()).fold(0., f64::max)
        / real.iter().map(|z| z.re.abs()).fold(0., f64::max);
    println!(
        "ionic: imaginary relative={imaginary:e}, translation={translation:e} Ha, single-atom modulus={unit:e}, direct={direct:e} Ha, alpha-Z per atom={} Ha",
        transform.alpha_z().unwrap()
    );
    assert!(imaginary < 1e-12 && translation < 1e-12 && unit < 1e-12 && direct < 1e-12);
}
#[test]
fn nonlocal_dense_phases_and_rotation() {
    let (cell, b, pp, pos) = setup();
    let n = b.vectors().len();
    let op = nl(&b, &cell, &pp, &pos);
    let reference = reference_nl(&b, &cell, &pp, &pos);
    let columns = reference_columns(&b, &cell, &pp, &pos);
    let mut column_error: f64 = 0.;
    for (atom, reference) in op.atoms().iter().zip(&columns) {
        for (a, r) in atom.columns().iter().zip(reference) {
            column_error = column_error.max(maxdiff(a, r));
        }
    }
    let delta = [0.31, -0.27, 0.18];
    let mut shifted = pos.clone();
    for (r, d) in shifted[1].iter_mut().zip(delta) {
        *r += d;
    }
    let other = nl(&b, &cell, &pp, &shifted);
    let mut translation: f64 = 0.;
    for (a, r) in other.atoms()[1]
        .columns()
        .iter()
        .zip(op.atoms()[1].columns())
    {
        for (i, g) in b.vectors().iter().enumerate() {
            let phase = C::from_polar(1., -(0..3).map(|k| g.cartesian[k] * delta[k]).sum::<f64>());
            translation = translation.max((a[i] - r[i] * phase).norm());
        }
    }
    let x = random(n, 17);
    let y = random(n, 42);
    let mut ax = vec![C::default(); n];
    let mut ay = ax.clone();
    op.apply(&x, &mut ax).unwrap();
    op.apply(&y, &mut ay).unwrap();
    let dense = maxdiff(&ax, &mv(&reference, &x));
    let herm = (dot(&x, &ay) - dot(&ax, &y)).norm() / dot(&x, &ay).norm().max(dot(&ax, &y).norm());
    let mut rotation: f64 = 0.;
    for (i, g) in b.vectors().iter().enumerate().skip(1) {
        let measured: f64 = op.atoms()[0].columns()[2..]
            .iter()
            .map(|c| c[i].norm_sqr())
            .sum();
        let expected = 3. / (4. * PI * cell.volume_bohr3())
            * pp.projector_transform(2, g.norm_squared.sqrt())
                .unwrap()
                .powi(2);
        rotation = rotation.max((measured - expected).abs() / expected);
    }
    println!(
        "nonlocal: N={n}, columns={column_error:e}, translation={translation:e}, dense={dense:e} Ha, Hermiticity relative={herm:e}, p rotation relative={rotation:e}"
    );
    assert!(
        column_error < 1e-12
            && translation < 1e-12
            && dense < 1e-12
            && herm < 1e-12
            && rotation < 1e-12
    );
}
#[test]
fn convolution_hamiltonian_and_solvers() {
    let (cell, b, pp, pos) = setup();
    let n = b.vectors().len();
    let density = b.density_basis(&cell).unwrap();
    let grid = FftGridSpec::new(&b, &density, None).unwrap();
    let mut coeff = local_coefficients();
    let ionic = ionic_local(
        &density,
        cell.volume_bohr3(),
        &[Species {
            potential: &pp,
            positions_bohr: &pos,
        }],
    )
    .unwrap();
    for (g, v) in density.vectors().iter().zip(ionic) {
        *coeff.entry(g.miller).or_default() += v;
    }
    let potential = sampled(grid.dims(), &coeff);
    let local = Hamiltonian::new(
        &cell,
        &b,
        potential.clone(),
        NonlocalOperator::new(&b, cell.volume_bohr3(), &[]).unwrap(),
    )
    .unwrap();
    let op = Hamiltonian::new(&cell, &b, potential, nl(&b, &cell, &pp, &pos)).unwrap();
    let mut reference = reference_nl(&b, &cell, &pp, &pos);
    let mut local_reference = ComplexMatrix::zeros(n, n).unwrap();
    for (i, g) in b.vectors().iter().enumerate() {
        for (j, h) in b.vectors().iter().enumerate() {
            let m = std::array::from_fn(|k| g.miller[k] - h.miller[k]);
            let v = coeff.get(&m).copied().unwrap_or_default();
            local_reference[(i, j)] = v;
            reference[(i, j)] += v;
            if i == j {
                reference[(i, j)] += 0.5 * g.norm_squared;
                local_reference[(i, j)] += 0.5 * g.norm_squared;
            }
        }
    }
    let x = random(n, 18);
    let y = random(n, 46);
    let ax = apply(&op, &x);
    let ay = apply(&op, &y);
    let convolution = maxdiff(&apply(&local, &x), &mv(&local_reference, &x));
    let mut entries: f64 = 0.;
    for j in 0..n {
        let mut unit = vec![C::default(); n];
        unit[j] = C::new(1., 0.);
        entries = entries.max(maxdiff(&apply(&op, &unit), reference.column(j)));
    }
    let herm = (dot(&x, &ay) - dot(&ax, &y)).norm() / dot(&x, &ay).norm().max(dot(&ax, &y).norm());
    assert_eq!(ax, apply(&op, &x));
    let block = ComplexMatrix::new(n, 2, [x.clone(), y.clone()].concat()).unwrap();
    let mut output = ComplexMatrix::zeros(n, 2).unwrap();
    op.apply_block(&block, &mut output).unwrap();
    let block_error = maxdiff(output.column(0), &ax).max(maxdiff(output.column(1), &ay));
    let constant = Hamiltonian::new(
        &cell,
        &b,
        Array3::from_elem(grid.dims(), 0.23),
        NonlocalOperator::new(&b, cell.volume_bohr3(), &[]).unwrap(),
    )
    .unwrap();
    let expected: Vec<_> = x
        .iter()
        .zip(b.vectors())
        .map(|(x, g)| x * (0.5 * g.norm_squared + 0.23))
        .collect();
    let constant_error = maxdiff(&apply(&constant, &x), &expected);
    let free = DenseEigensolver.solve(&constant, n).unwrap();
    let mut exact: Vec<_> = b
        .vectors()
        .iter()
        .map(|g| 0.5 * g.norm_squared + 0.23)
        .collect();
    exact.sort_by(f64::total_cmp);
    let free_error = free
        .eigenvalues
        .iter()
        .zip(exact)
        .map(|(a, b)| (a - b).abs())
        .fold(0., f64::max);
    let (values, _) = hermitian_eigh(&reference).unwrap();
    let dense = DenseEigensolver.solve(&op, 8).unwrap();
    let dense_error = dense
        .eigenvalues
        .iter()
        .zip(&values)
        .map(|(a, b)| (a - b).abs())
        .fold(0., f64::max);
    let solver = Davidson {
        bands: 8,
        tolerance: 1e-10,
        max_iterations: 300,
        max_subspace: 32,
        block_size: 8,
        seed: 2026,
    };
    let pre = DiagonalPreconditioner::new(op.kinetic_diagonal().to_vec(), 0.1).unwrap();
    let result = solver.solve(&op, None, Some(&pre)).unwrap();
    let eig = result
        .eigenvalues
        .iter()
        .zip(&values)
        .map(|(a, b)| (a - b).abs())
        .fold(0., f64::max);
    let residual = result.residuals.iter().copied().fold(0., f64::max);
    println!(
        "Hamiltonian: grid={:?}, convolution={convolution:e} Ha, entries={entries:e} Ha, Hermiticity relative={herm:e}, block={block_error:e} Ha, constant={constant_error:e} Ha, free eigenvalues={free_error:e} Ha, dense eigenvalues={dense_error:e} Ha",
        grid.dims()
    );
    println!(
        "Davidson: eigenvalue error={eig:e} Ha, residual={residual:e} Ha, iterations={}, eigenvalues={:?} Ha",
        result.iterations, result.eigenvalues
    );
    assert!(
        convolution < 1e-12
            && entries < 1e-11
            && herm < 1e-12
            && block_error < 1e-14
            && constant_error < 1e-12
            && free_error < 1e-12
            && dense_error < 1e-11
            && eig < 1e-9
            && residual < 1e-10
    );
}
#[test]
fn filtering_and_aliasing() {
    let (cell, b, _, _) = setup();
    let density = b.density_basis(&cell).unwrap();
    let grid = FftGridSpec::new(&b, &density, None).unwrap();
    let dims = grid.dims();
    let coeff = local_coefficients();
    let v = sampled(dims, &coeff);
    let filtered = truncate_to_density_sphere(&v, &density).unwrap();
    let retained = v
        .iter()
        .zip(&filtered)
        .map(|(a, b)| (a - b).abs())
        .fold(0., f64::max);
    assert!(retained < 1e-12);
    // A true mode n1+1 is outside the represented range and density sphere.
    // Its samples coincide with mode 1, producing an erroneous low-mode coupling.
    let high = [dims[0] as i32 + 1, 0, 0];
    let high_coeff = BTreeMap::from([(high, C::new(0.5, 0.)), (high.map(|m| -m), C::new(0.5, 0.))]);
    let high_values = sampled(dims, &high_coeff);
    let op = Hamiltonian::new(
        &cell,
        &b,
        high_values.clone(),
        NonlocalOperator::new(&b, cell.volume_bohr3(), &[]).unwrap(),
    )
    .unwrap();
    let mut unit = vec![C::default(); b.vectors().len()];
    unit[0] = C::new(1., 0.);
    let y = apply(&op, &unit);
    let alias = y.iter().map(|z| z.norm()).fold(0., f64::max);
    let filtered_high = truncate_to_density_sphere(&high_values, &density).unwrap();
    let unrecoverable = filtered_high
        .iter()
        .zip(&high_values)
        .map(|(a, b)| (a - b).abs())
        .fold(0., f64::max);
    // A represented mode outside the sphere is removable before multiplication.
    let outside = (0..dims[0])
        .flat_map(|i| {
            (0..dims[1])
                .flat_map(move |j| (0..dims[2]).map(move |k| [i as i32, j as i32, k as i32]))
        })
        .find(|m| {
            m.iter().zip(dims).all(|(&m, n)| m <= ((n - 1) / 2) as i32)
                && !density.vectors().iter().any(|g| g.miller == *m)
        })
        .unwrap();
    let out_values = sampled(
        dims,
        &BTreeMap::from([
            (outside, C::new(0.5, 0.)),
            (outside.map(|m| -m), C::new(0.5, 0.)),
        ]),
    );
    let removed = truncate_to_density_sphere(&out_values, &density)
        .unwrap()
        .iter()
        .map(|x| x.abs())
        .fold(0., f64::max);
    println!(
        "filter: retained error={retained:e} Ha, outside-sphere removal={removed:e} Ha; alias high mode={high:?}, spurious matrix element={alias:e} Ha, filtering already-aliased change={unrecoverable:e} Ha"
    );
    assert!(removed < 1e-12 && (alias - 0.5).abs() < 1e-12 && unrecoverable < 1e-12);
}
#[test]
fn invalid_inputs() {
    let (cell, b, pp, pos) = setup();
    let density = b.density_basis(&cell).unwrap();
    let grid = FftGridSpec::new(&b, &density, None).unwrap();
    assert!(matches!(
        NonlocalOperator::new(&b, 0., &[]),
        Err(HamiltonianError::Volume)
    ));
    assert!(matches!(
        structure_factor([1., 0., 0.], &[[f64::NAN, 0., 0.]]),
        Err(HamiltonianError::NonFinite)
    ));
    let op = nl(&b, &cell, &pp, &pos);
    assert!(matches!(
        op.apply(&[], &mut []),
        Err(HamiltonianError::Shape)
    ));
    let h = Hamiltonian::new(&cell, &b, Array3::zeros(grid.dims()), op).unwrap();
    assert!(h.apply(&[], &mut []).is_err());
    assert!(
        h.apply(
            &vec![C::new(f64::NAN, 0.); h.dimension()],
            &mut vec![C::default(); h.dimension()]
        )
        .is_err()
    );
    assert!(
        Hamiltonian::new(
            &cell,
            &b,
            Array3::from_elem(grid.dims(), f64::NAN),
            nl(&b, &cell, &pp, &pos)
        )
        .is_err()
    );
    let x = ComplexMatrix::zeros(h.dimension() + 1, 1).unwrap();
    let mut y = ComplexMatrix::zeros(h.dimension(), 1).unwrap();
    assert!(h.apply_block(&x, &mut y).is_err());
    // Existing fixture entry point remains valid and its original resolution unchanged.
    assert_eq!(
        parse_upf(&synthetic::fixture()).unwrap().mesh().r().len(),
        601
    );
}
#[test]
#[ignore = "large-cutoff projector norm; run in release mode"]
fn projector_norm_convergence() {
    let (cell, _, pp, _) = setup();
    for cutoff in [8., 32.] {
        let b = PlaneWaveBasis::new(&cell, Hartree(cutoff)).unwrap();
        let op = nl(&b, &cell, &pp, &[[0.; 3]]);
        for (c, col) in op.atoms()[0].columns().iter().enumerate() {
            let p = &pp.projectors()[if c < 2 { c } else { 2 }];
            let radial = pp
                .mesh()
                .integrate(&p.radial_times_r().iter().map(|u| u * u).collect::<Vec<_>>())
                .unwrap();
            let reciprocal: f64 = col.iter().map(|z| z.norm_sqr()).sum();
            println!(
                "projector norm: cutoff={} Ry, N={}, channel={c}, radial={radial:e}, reciprocal={reciprocal:e}, relative difference={:e} (dimensionless; periodic-image overlaps remain)",
                2. * cutoff,
                b.vectors().len(),
                (reciprocal - radial).abs() / radial
            );
            assert!(reciprocal.is_finite() && radial > 0.);
        }
    }
}

// Direct real-space Gaussian overlap, independent of radial transforms, Bessel
// functions, harmonics, and reciprocal normalization in the production path.
// Completing the square gives I_s(R)/I_s(0)=exp(-a R^2/2).
// For a Cartesian p harmonic the shifted Gaussian second moment contributes
// (1-a R_axis^2), where a=1/3 Bohr^-2 and m=-1,0,1 corresponds to y,z,x.
fn gaussian_image_correction(cell: &Cell, channel: usize, extent: i32) -> f64 {
    let a = 1. / (if channel < 2 { channel + 1 } else { 3 }) as f64;
    let mut correction = 0.;
    for i in -extent..=extent {
        for j in -extent..=extent {
            for k in -extent..=extent {
                if [i, j, k] == [0; 3] {
                    continue;
                }
                let r = cell.fractional_to_cartesian_bohr([i as f64, j as f64, k as f64]);
                let r2: f64 = r.iter().map(|x| x * x).sum();
                let polynomial = if channel < 2 {
                    1.
                } else {
                    1. - a * r[[1, 2, 0][channel - 2]].powi(2)
                };
                correction += polynomial * (-0.5 * a * r2).exp();
            }
        }
    }
    correction
}

#[test]
#[ignore = "large-cell projector normalization and image-overlap evidence; run in release mode"]
fn projector_norm_isolation_and_images() {
    let (si, _, pp, _) = setup();
    let radial: Vec<_> = pp
        .projectors()
        .iter()
        .map(|p| {
            pp.mesh()
                .integrate(&p.radial_times_r().iter().map(|u| u * u).collect::<Vec<_>>())
                .unwrap()
        })
        .collect();
    let side = 40.;
    let cubic =
        Cell::from_user_rows_bohr([[side, 0., 0.], [0., side, 0.], [0., 0., side]]).unwrap();
    let support = pp
        .projectors()
        .iter()
        .map(|p| pp.mesh().r()[p.cutoff_index() - 1])
        .fold(0., f64::max);
    // Surface-to-surface distance between support balls, not center distance.
    let margin = side - 2. * support;
    assert!(margin >= 10.);
    println!(
        "isolated support: cubic side={side} Bohr, radial cutoff={support:.15} Bohr, nearest-image support gap={margin:.15} Bohr"
    );
    let mut previous = [f64::INFINITY; 5];
    for cutoff_ry in [8., 16., 32., 64.] {
        let b = PlaneWaveBasis::new(&cubic, Hartree(cutoff_ry / 2.)).unwrap();
        let op = nl(&b, &cubic, &pp, &[[0.; 3]]);
        for (channel, col) in op.atoms()[0].columns().iter().enumerate() {
            let index = if channel < 2 { channel } else { 2 };
            let m = if channel < 2 { 0 } else { channel as isize - 3 };
            let sum: f64 = col.iter().map(|p| p.norm_sqr()).sum();
            let error = (sum / radial[index] - 1.).abs();
            println!(
                "isolated norm: cutoff={cutoff_ry} Ry, N={}, radial_index={index}, m={m}, radial={:.16e}, reciprocal={sum:.16e}, relative_difference={error:.16e}",
                b.vectors().len(),
                radial[index]
            );
            // Independent Cartesian Gaussian Fourier integrals. This checks
            // finite-cutoff deficits separately from normalization mistakes.
            let a = 1. / (index + 1) as f64;
            let analytic: f64 = b
                .vectors()
                .iter()
                .map(|g| {
                    let base = PI.powi(2) / (4. * cubic.volume_bohr3() * a.powi(3))
                        * (-g.norm_squared / (2. * a)).exp();
                    if channel < 2 {
                        base
                    } else {
                        base * 3. * g.cartesian[[1, 2, 0][channel - 2]].powi(2) / (4. * a * a)
                    }
                })
                .sum();
            let analytic_error = (sum - analytic).abs() / radial[index];
            println!(
                "analytic Gaussian: cutoff={cutoff_ry} Ry, channel={channel}, analytic_relative_deficit={:.16e}, numerical_vs_analytic={analytic_error:.16e}",
                1. - analytic / radial[index]
            );
            assert!(analytic_error < 1e-10);
            // The cutoff deficit decreases until radial-quadrature/roundoff
            // dominates. Do not require strict monotonicity below 1e-8.
            assert!(error < previous[channel] || error < 1e-8);
            if cutoff_ry == 64. {
                assert!(
                    error < 1e-8,
                    "isolated normalization defect, channel={channel}: {error:e}"
                );
            }
            previous[channel] = error;
        }
    }
    // At 64 Ry the earlier Si cutoff sweep is saturated. Compare its signed
    // excess/deficit to the independently integrated real-space image sum.
    let b = PlaneWaveBasis::new(&si, Hartree(32.)).unwrap();
    let op = nl(&b, &si, &pp, &[[0.; 3]]);
    for (channel, col) in op.atoms()[0].columns().iter().enumerate() {
        let index = if channel < 2 { channel } else { 2 };
        let sum: f64 = col.iter().map(|p| p.norm_sqr()).sum();
        let observed = sum / radial[index] - 1.;
        let predicted = gaussian_image_correction(&si, channel, 4);
        let refined = gaussian_image_correction(&si, channel, 5);
        assert!((predicted - refined).abs() < 1e-15);
        let relative = (observed - predicted).abs() / predicted.abs();
        println!(
            "Si image correction: channel={channel}, observed={observed:.16e}, predicted={predicted:.16e}, relative_disagreement={relative:.16e}, image_norm_term={:.16e} (dimensionless)",
            predicted * radial[index]
        );
        // s0's image signal (~4e-11) is comparable to absolute quadrature
        // error; the required 1e-3 relative certificate concerns s1 and p.
        if channel > 0 {
            assert!(
                relative < 1e-3,
                "image-overlap prediction failed for channel={channel}"
            );
        } else {
            assert!((observed - predicted).abs() < 1e-8);
        }
    }
}

// Pre-O3 apply oracle: preserve scatter, full FFTs, potential multiplication,
// gather-plus-kinetic and nonlocal accumulation order exactly.
fn full_fft_apply_oracle(
    h: &Hamiltonian,
    grid: &FftGridSpec,
    potential: &Array3<f64>,
    x: &[C],
) -> Vec<C> {
    let mut fft = FftPlan3d::new(grid.dims()).unwrap();
    let mut reciprocal = grid.scatter(x).unwrap();
    let mut real = Array3::zeros(grid.dims());
    fft.to_real(&reciprocal, &mut real).unwrap();
    for (z, v) in real.iter_mut().zip(potential.iter()) {
        *z *= v;
    }
    fft.to_reciprocal(&real, &mut reciprocal).unwrap();
    let mut y: Vec<_> = grid
        .g_to_fft()
        .iter()
        .enumerate()
        .map(|(i, &idx)| reciprocal[idx] + h.kinetic_diagonal()[i] * x[i])
        .collect();
    for atom in h.nonlocal().atoms() {
        let projections: Vec<C> = atom
            .columns()
            .iter()
            .map(|col| col.iter().zip(x).map(|(b, c)| b.conj() * c).sum())
            .collect();
        let nc = projections.len();
        for (i, col) in atom.columns().iter().enumerate() {
            let weight: C = (0..nc)
                .map(|j| atom.coupling_hartree()[i * nc + j] * projections[j])
                .sum();
            for (out, b) in y.iter_mut().zip(col) {
                *out += b * weight;
            }
        }
    }
    y
}

#[test]
fn pruned_apply_bitwise_gamma_and_shifted() {
    let (cell, _, pp, pos) = setup();
    for k in [[0.; 3], [0.25, -0.25, 0.125]] {
        let basis = PlaneWaveBasis::at_k(&cell, Hartree(3.), k).unwrap();
        let grid = FftGridSpec::new(&basis, &basis.density_basis(&cell).unwrap(), None).unwrap();
        let potential = Array3::from_shape_vec(
            grid.dims(),
            random(grid.dims().iter().product(), 12)
                .iter()
                .map(|z| z.re)
                .collect(),
        )
        .unwrap();
        let h = Hamiltonian::new(
            &cell,
            &basis,
            potential.clone(),
            nl(&basis, &cell, &pp, &pos),
        )
        .unwrap();
        for seed in [1, 42, 765] {
            let x = random(basis.vectors().len(), seed);
            let expected = full_fft_apply_oracle(&h, &grid, &potential, &x);
            let actual = apply(&h, &x);
            for (a, b) in actual.iter().zip(expected) {
                assert_eq!(
                    (a.re.to_bits(), a.im.to_bits()),
                    (b.re.to_bits(), b.im.to_bits())
                );
            }
        }
    }
}

#[test]
fn band_workers_bitwise_and_ordered_errors() {
    let (cell, _, pp, pos) = setup();
    for k in [[0.; 3], [0.25, -0.25, 0.125]] {
        let basis = PlaneWaveBasis::at_k(&cell, Hartree(3.), k).unwrap();
        let grid = FftGridSpec::new(&basis, &basis.density_basis(&cell).unwrap(), None).unwrap();
        let potential = Array3::from_shape_vec(
            grid.dims(),
            random(grid.dims().iter().product(), 12)
                .iter()
                .map(|z| z.re)
                .collect(),
        )
        .unwrap();
        let n = basis.vectors().len();
        let block = ComplexMatrix::new(n, 7, random(n * 7, 927)).unwrap();
        let mut expected = ComplexMatrix::zeros(n, 7).unwrap();
        let reference = Hamiltonian::new(
            &cell,
            &basis,
            potential.clone(),
            nl(&basis, &cell, &pp, &pos),
        )
        .unwrap();
        reference.apply_block(&block, &mut expected).unwrap();
        for threads in [1, 2, 3, 4, 8] {
            let h = Hamiltonian::new(
                &cell,
                &basis,
                potential.clone(),
                nl(&basis, &cell, &pp, &pos),
            )
            .unwrap()
            .with_band_threads(threads)
            .unwrap();
            let mut result = ComplexMatrix::zeros(n, 7).unwrap();
            for _ in 0..3 {
                h.apply_block(&block, &mut result).unwrap();
                for (a, b) in result.as_slice().iter().zip(expected.as_slice()) {
                    assert_eq!(
                        (a.re.to_bits(), a.im.to_bits()),
                        (b.re.to_bits(), b.im.to_bits())
                    );
                }
            }
            let mut bad = block.clone();
            bad.column_mut(2)[0].re = f64::NAN;
            bad.column_mut(5)[0].im = f64::INFINITY;
            assert!(
                h.apply_block(&bad, &mut result)
                    .unwrap_err()
                    .to_string()
                    .contains("band index 2:")
            );
            // All scopes have joined and released scratch even on failure.
            h.apply_block(&block, &mut result).unwrap();
            assert_eq!(result, expected);
        }
    }
}

#[test]
#[ignore = "release microbenchmark: cargo test --release -p dft-physics --test hamiltonian conventional_cell_block_timing -- --ignored --nocapture"]
fn conventional_cell_block_timing() {
    // Synthetic potential on the requested conventional-cell geometry. This
    // measures block application only, not a converged real-silicon SCF.
    let a = Bohr::from(Angstrom(5.43)).0;
    let cell = Cell::from_user_rows_bohr([[a, 0., 0.], [0., a, 0.], [0., 0., a]]).unwrap();
    let basis = PlaneWaveBasis::new(&cell, Hartree(15.)).unwrap();
    let grid = FftGridSpec::new(&basis, &basis.density_basis(&cell).unwrap(), None).unwrap();
    let potential =
        Array3::from_shape_fn(grid.dims(), |(i, j, k)| ((i + 2 * j + 3 * k) as f64).sin());
    let block = ComplexMatrix::new(
        basis.vectors().len(),
        20,
        random(basis.vectors().len() * 20, 73),
    )
    .unwrap();
    let mut expected = ComplexMatrix::zeros(block.rows(), block.cols()).unwrap();
    for threads in [1, 2, 4, 8] {
        let h = Hamiltonian::new(
            &cell,
            &basis,
            potential.clone(),
            NonlocalOperator::new(&basis, cell.volume_bohr3(), &[]).unwrap(),
        )
        .unwrap()
        .with_band_threads(threads)
        .unwrap();
        let mut output = expected.clone();
        h.apply_block(&block, &mut output).unwrap();
        if threads == 1 {
            expected = output.clone();
        }
        assert!(
            output
                .as_slice()
                .iter()
                .zip(expected.as_slice())
                .all(|(a, b)| a.re.to_bits() == b.re.to_bits() && a.im.to_bits() == b.im.to_bits())
        );
        let mut samples = Vec::new();
        for _ in 0..11 {
            let start = std::time::Instant::now();
            for _ in 0..10 {
                h.apply_block(&block, &mut output).unwrap();
            }
            samples.push(start.elapsed().as_secs_f64() / 10.);
        }
        samples.sort_by(f64::total_cmp);
        println!(
            "synthetic conventional cell: grid={:?} npw={} bands=20 threads={threads} block_median_s={:.9} bit_identical=true",
            grid.dims(),
            block.rows(),
            samples[5]
        );
    }
}
