//! Component oracles on normalized ORC states; synthetic NC data, no QE claims.
use dft_basis::{FftGridSpec, HalfSphere, PlaneWaveBasis};
use dft_geometry::Cell;
use dft_numerics::{RealMatrix, linalg::ComplexMatrix, pool::WorkerPool};
use dft_physics::{
    RealHamiltonian, RealNonlocal,
    hamiltonian::{Hamiltonian, NonlocalOperator, Species},
};
use dft_solvers::LinearOperator;
use dft_units::Hartree;
use ndarray::Array3;
use num_complex::Complex64;

fn pseudo(ls: &[usize]) -> dft_upf::Pseudopotential {
    let n = 81;
    let a = |f: &dyn Fn(f64) -> f64| {
        (0..n)
            .map(|i| format!("{:.17e}", f(i as f64 * 0.15)))
            .collect::<Vec<_>>()
            .join(" ")
    };
    let mut beta = String::new();
    for (j, &l) in ls.iter().enumerate() {
        let i = j + 1;
        beta.push_str(&format!("<PP_BETA.{i} index=\"{i}\" angular_momentum=\"{l}\" cutoff_radius=\"12\" cutoff_radius_index=\"{n}\">{}</PP_BETA.{i}>",a(&|r| r.powi(l as i32+1)*(-r*r/(1.+j as f64*0.1)).exp())));
    }
    let dij = (0..ls.len())
        .flat_map(|i| {
            (0..ls.len()).map(move |j| {
                if i == j {
                    1. + i as f64 * 0.1
                } else if ls[i] == ls[j] {
                    0.2
                } else {
                    0.
                }
            })
        })
        .map(|v| v.to_string())
        .collect::<Vec<_>>()
        .join(" ");
    dft_upf::parse_upf(&format!(r#"<UPF version="2.0.1"><PP_INFO>Synthetic R3 fixture</PP_INFO>
<PP_HEADER generated="synthetic" author="tests" date="2026" comment="synthetic" element="Xx" pseudo_type="NC" relativistic="scalar" is_ultrasoft="F" is_paw="F" is_coulomb="F" has_so="F" has_wfc="F" has_gipaw="F" core_correction="F" functional="SLA PZ" z_valence="4" total_psenergy="0" wfc_cutoff="0" rho_cutoff="0" l_max="3" l_max_rho="0" l_local="-1" mesh_size="{n}" number_of_wfc="0" number_of_proj="{}"/>
<PP_MESH dx="0.15" mesh="{n}" xmin="0" rmax="12" zmesh="1"><PP_R>{}</PP_R><PP_RAB>{}</PP_RAB></PP_MESH><PP_LOCAL>{}</PP_LOCAL><PP_NONLOCAL>{beta}<PP_DIJ>{dij}</PP_DIJ></PP_NONLOCAL><PP_RHOATOM>{}</PP_RHOATOM></UPF>"#,ls.len(),a(&|r|r),a(&|_|0.15),a(&|_|0.),a(&|_|0.))).unwrap()
}
fn fixture(rows: [[f64; 3]; 3], atoms: usize, cutoff: f64) -> (PlaneWaveBasis, Hamiltonian) {
    fixture_with(rows, atoms, cutoff, &[0, 0, 1, 1, 2, 2, 3, 3], &[0, 1])
}
fn fixture_with(
    rows: [[f64; 3]; 3],
    atoms: usize,
    cutoff: f64,
    first: &[usize],
    second: &[usize],
) -> (PlaneWaveBasis, Hamiltonian) {
    let cell = Cell::from_user_rows_bohr(rows).unwrap();
    let basis = PlaneWaveBasis::new(&cell, Hartree(cutoff)).unwrap();
    let grid = FftGridSpec::new(&basis, &basis.density_basis(&cell).unwrap(), None).unwrap();
    let positions: Vec<_> = (0..atoms)
        .map(|i| {
            cell.fractional_to_cartesian_bohr([
                ((i * 7 + 1) % 23) as f64 / 23.,
                ((i * 11 + 2) % 29) as f64 / 29.,
                ((i * 13 + 3) % 31) as f64 / 31.,
            ])
        })
        .collect();
    let p = pseudo(first);
    let q = pseudo(second);
    let split = atoms.div_ceil(2);
    let nl = NonlocalOperator::new(
        &basis,
        cell.volume_bohr3(),
        &[
            Species {
                potential: &p,
                positions_bohr: &positions[..split],
            },
            Species {
                potential: &q,
                positions_bohr: &positions[split..],
            },
        ],
    )
    .unwrap();
    let v = Array3::from_shape_fn(grid.dims(), |(i, j, k)| {
        0.3 * ((i * 13 + j * 7 + k * 3) as f64).sin()
    });
    let h = Hamiltonian::new(&cell, &basis, v, nl).unwrap();
    (basis, h)
}
fn states(n: usize, m: usize, seed: usize) -> RealMatrix {
    let mut x = RealMatrix::new(
        n,
        m,
        (0..n * m).map(|i| ((i * 17 + seed) as f64).sin()).collect(),
    )
    .unwrap();
    for col in x.as_mut_slice().chunks_mut(n) {
        let norm = col.iter().map(|v| v * v).sum::<f64>().sqrt();
        for v in col {
            *v /= norm;
        }
    }
    x
}
fn complete(half: &HalfSphere, x: &RealMatrix) -> ComplexMatrix {
    let mut y = ComplexMatrix::zeros(x.rows(), x.cols()).unwrap();
    for (a, b) in x
        .as_slice()
        .chunks(x.rows())
        .zip(y.as_mut_slice().chunks_mut(x.rows()))
    {
        half.unpack(a, b).unwrap();
    }
    y
}
fn relative(a: &[Complex64], b: &[Complex64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(a, b)| (a - b).norm_sqr())
        .sum::<f64>()
        .sqrt()
        / b.iter().map(|b| b.norm_sqr()).sum::<f64>().sqrt()
}
#[test]
fn ten_shapes_nonlocal_full_h_hermiticity_and_threads() {
    let mut max_h = 0.0_f64;
    let mut max_nl = 0.0_f64;
    let mut max_sym = 0.0_f64;
    for case in 0..10 {
        let s = case as f64 * 0.13;
        let rows = match case % 4 {
            0 => [[6. + s, 0., 0.], [0., 6. + s, 0.], [0., 0., 6. + s]],
            1 => [[6. + s, 0.4, 0.2], [1., 7. + s, 0.3], [0.5, 0.7, 8. + s]],
            2 => [[4. + s, 0., 0.], [0., 5. + s, 0.], [0., 0., 14. + s]],
            _ => [[7. + s, 0., 0.], [0., 8. + s, 0.], [0.3, 0., 18. + s]],
        };
        let atoms = [1, 2, 4, 8, 12, 16, 20, 3, 6, 10][case];
        let (basis, h) = fixture(rows, atoms, 1.2);
        let n = basis.vectors().len();
        let m = 3 + case % 2;
        let x = states(n, m, 3);
        let z = states(n, m, 37);
        let mut hx = RealMatrix::zeros(n, m).unwrap();
        let mut hz = hx.clone();
        let mut reference_bits = Vec::new();
        for threads in [1, 2, 4, 8, 12] {
            let pool = WorkerPool::new(threads).unwrap();
            let mut real = RealHamiltonian::from_complex(&h, &basis, m + 1, pool.clone()).unwrap();
            real.apply_block(&x, &mut hx).unwrap();
            real.apply_block(&z, &mut hz).unwrap();
            let bits: Vec<_> = hx.as_slice().iter().map(|v| v.to_bits()).collect();
            if threads == 1 {
                reference_bits = bits;
            } else {
                assert_eq!(bits, reference_bits);
            }
            let cx = complete(real.half_sphere(), &x);
            let rhx = complete(real.half_sphere(), &hx);
            let mut chx = ComplexMatrix::zeros(n, m).unwrap();
            h.apply_block(&cx, &mut chx).unwrap();
            let err = relative(rhx.as_slice(), chx.as_slice());
            max_h = max_h.max(err);
            // 1e-13 admits different FFT pairing and ~thousands-term GEMM sums,
            // while retaining the design's fixed component-level accuracy gate.
            assert!(err < 1e-13, "case {case}: H relative={err:e}");
            let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
            let left = dot(x.as_slice(), hz.as_slice());
            let right = dot(z.as_slice(), hx.as_slice());
            let sym = (left - right).abs() / left.abs().max(right.abs());
            max_sym = max_sym.max(sym);
            assert!(sym < 1e-13, "case {case}: Hermiticity={sym:e}");
            let mut nl = RealNonlocal::new(h.nonlocal(), real.half_sphere(), m, pool).unwrap();
            let mut ny = RealMatrix::zeros(n, m).unwrap();
            nl.add_block(&x, &mut ny).unwrap();
            let rny = complete(real.half_sphere(), &ny);
            let mut cny = ComplexMatrix::zeros(n, m).unwrap();
            for (x, y) in cx
                .as_slice()
                .chunks(n)
                .zip(cny.as_mut_slice().chunks_mut(n))
            {
                h.nonlocal().apply(x, y).unwrap();
            }
            let err = relative(rny.as_slice(), cny.as_slice());
            max_nl = max_nl.max(err);
            assert!(err < 1e-13);
        }
        println!("case={case} atoms={atoms} rows={n} bands={m}");
    }
    println!(
        "max H relative={max_h:e}; nonlocal relative={max_nl:e}; Hermiticity relative={max_sym:e}"
    );
}
/// B6 physics half: for X = V C the cached projections (P^T V) C reproduce
/// H X without the first projector product, to rounding (relative 1e-13), on a
/// skewed two-species cell with s, p, d, f projectors; shape errors are typed.
#[test]
fn cached_projections_reproduce_apply() {
    let (basis, h) = fixture([[6., 0.4, 0.2], [1., 7., 0.3], [0.5, 0.7, 8.]], 5, 1.2);
    let n = basis.vectors().len();
    let (m, k) = (6, 4);
    let pool = WorkerPool::new(2).unwrap();
    let mut real = RealHamiltonian::from_complex(&h, &basis, m, pool).unwrap();
    let count = real.projector_count();
    assert!(count > 0);
    let v = states(n, m, 11);
    let mut hv = RealMatrix::zeros(n, m).unwrap();
    real.apply_block(&v, &mut hv).unwrap();
    let pv = real.last_projections(m).unwrap().to_vec();
    assert_eq!(pv.len(), count * m);
    // X = V C with a fixed dense C (m x k, column-major).
    let c: Vec<f64> = (0..m * k).map(|i| ((i * 5 + 1) as f64).cos()).collect();
    let mut x = RealMatrix::zeros(n, k).unwrap();
    for j in 0..k {
        for i in 0..m {
            let cij = c[j * m + i];
            for r in 0..n {
                x.as_mut_slice()[j * n + r] += cij * v.as_slice()[i * n + r];
            }
        }
    }
    let mut pc = vec![0.; count * k];
    for j in 0..k {
        for i in 0..m {
            for p in 0..count {
                pc[j * count + p] += pv[i * count + p] * c[j * m + i];
            }
        }
    }
    let mut reference = RealMatrix::zeros(n, k).unwrap();
    real.apply_block(&x, &mut reference).unwrap();
    let direct = real.last_projections(k).unwrap().to_vec();
    let scale = direct.iter().map(|v| v * v).sum::<f64>().sqrt();
    let diff = direct
        .iter()
        .zip(&pc)
        .map(|(a, b)| (a - b) * (a - b))
        .sum::<f64>()
        .sqrt();
    assert!(diff <= 1e-13 * scale, "projections {diff:e} of {scale:e}");
    let mut cached = RealMatrix::zeros(n, k).unwrap();
    real.apply_block_with_projections(&x, &pc, &mut cached)
        .unwrap();
    let norm = reference
        .as_slice()
        .iter()
        .map(|v| v * v)
        .sum::<f64>()
        .sqrt();
    let err = reference
        .as_slice()
        .iter()
        .zip(cached.as_slice())
        .map(|(a, b)| (a - b) * (a - b))
        .sum::<f64>()
        .sqrt()
        / norm;
    println!("B6 cached-projection H relative={err:e}, projections={count}");
    assert!(err < 1e-13, "relative {err:e}");
    // Supplying the direct projections reproduces apply_block bit for bit.
    real.apply_block_with_projections(&x, &direct, &mut cached)
        .unwrap();
    assert_eq!(cached.as_slice(), reference.as_slice());
    assert!(
        real.apply_block_with_projections(&x, &pc[..count], &mut cached)
            .is_err()
    );
}
/// NL-ENERGY-REAL: real per-column nonlocal expectation values agree with the
/// complex projector sums of the completed vectors to rounding, for any thread
/// budget bit for bit.
#[test]
fn real_nonlocal_expectation_matches_complex() {
    let (basis, h) = fixture([[6., 0.4, 0.2], [1., 7., 0.3], [0.5, 0.7, 8.]], 5, 1.2);
    let n = basis.vectors().len();
    let m = 7;
    let x = states(n, m, 21);
    let mut reference: Option<Vec<f64>> = None;
    for threads in [1, 3] {
        let mut real =
            RealHamiltonian::from_complex(&h, &basis, m, WorkerPool::new(threads).unwrap())
                .unwrap();
        let mut values = vec![0.; m];
        real.nonlocal_expectation_values(&x, &mut values).unwrap();
        let full = complete(real.half_sphere(), &x);
        let columns: Vec<&[Complex64]> = (0..m).map(|j| full.column(j)).collect();
        let mut expected = vec![0.; m];
        h.nonlocal()
            .expectation_values(&columns, &mut expected)
            .unwrap();
        for (a, b) in values.iter().zip(&expected) {
            assert!((a - b).abs() <= 1e-13 * b.abs().max(1e-3), "{a:e} vs {b:e}");
        }
        assert!(values.iter().any(|v| v.abs() > 1e-6));
        match &reference {
            None => reference = Some(values),
            Some(r) => assert_eq!(r, &values),
        }
    }
}
/// APPLY-LEAN (`BAANDI_APPLY_LEAN=on`) is bit-identical to the default apply,
/// for plain and cached-projection applications, and rejects non-finite input
/// with the same error.
#[test]
fn apply_lean_is_bit_identical() {
    let (basis, h) = fixture([[6., 0.4, 0.2], [1., 7., 0.3], [0.5, 0.7, 8.]], 5, 1.2);
    let n = basis.vectors().len();
    for m in [1, 2, 5] {
        let mut real =
            RealHamiltonian::from_complex(&h, &basis, m, WorkerPool::new(2).unwrap()).unwrap();
        let x = states(n, m, 3 + m);
        let mut plain = RealMatrix::zeros(n, m).unwrap();
        let mut lean = RealMatrix::zeros(n, m).unwrap();
        real.set_apply_lean(false);
        real.apply_block(&x, &mut plain).unwrap();
        let projections = real.last_projections(m).unwrap().to_vec();
        let mut plain_cached = RealMatrix::zeros(n, m).unwrap();
        real.apply_block_with_projections(&x, &projections, &mut plain_cached)
            .unwrap();
        real.set_apply_lean(true);
        real.apply_block(&x, &mut lean).unwrap();
        assert_eq!(plain.as_slice(), lean.as_slice());
        assert_eq!(real.last_projections(m).unwrap(), projections.as_slice());
        let mut lean_cached = RealMatrix::zeros(n, m).unwrap();
        real.apply_block_with_projections(&x, &projections, &mut lean_cached)
            .unwrap();
        assert_eq!(plain_cached.as_slice(), lean_cached.as_slice());
        let mut bad = x.clone();
        bad.as_mut_slice()[n / 2] = f64::NAN;
        for on in [false, true] {
            real.set_apply_lean(on);
            let err = real.apply_block(&bad, &mut lean).unwrap_err().to_string();
            assert_eq!(err, dft_physics::RealGammaError::NonFinite.to_string());
        }
    }
}
fn bits(x: &RealMatrix) -> Vec<u64> {
    x.as_slice().iter().map(|v| v.to_bits()).collect()
}
/// HAPPLY-OVERLAP (`BAANDI_HAPPLY_OVERLAP=on`): local part and nonlocal image
/// in one pool join are bit-identical to the serial order, for plain and
/// cached-projection applications, widths 1..8 and pool budgets 1, 2, 4; the
/// stored projections are unchanged and non-finite input is rejected alike.
#[test]
fn apply_overlap_is_bit_identical() {
    for rows in [
        [[6., 0.4, 0.2], [1., 7., 0.3], [0.5, 0.7, 8.]],
        [[4., 0., 0.], [0., 5., 0.], [0., 0., 14.]],
    ] {
        let (basis, h) = fixture(rows, 5, 1.2);
        let n = basis.vectors().len();
        for m in [1, 2, 5, 8] {
            let x = states(n, m, 5 + m);
            let mut reference: Option<(Vec<u64>, Vec<u64>)> = None;
            for threads in [1, 2, 4] {
                let mut real =
                    RealHamiltonian::from_complex(&h, &basis, m, WorkerPool::new(threads).unwrap())
                        .unwrap();
                let mut plain = RealMatrix::zeros(n, m).unwrap();
                let mut over = RealMatrix::zeros(n, m).unwrap();
                real.set_apply_overlap(false);
                real.apply_block(&x, &mut plain).unwrap();
                let projections = real.last_projections(m).unwrap().to_vec();
                let mut plain_cached = RealMatrix::zeros(n, m).unwrap();
                real.apply_block_with_projections(&x, &projections, &mut plain_cached)
                    .unwrap();
                real.set_apply_overlap(true);
                real.apply_block(&x, &mut over).unwrap();
                assert_eq!(bits(&plain), bits(&over), "m={m} T={threads}");
                assert_eq!(real.last_projections(m).unwrap(), projections.as_slice());
                let mut over_cached = RealMatrix::zeros(n, m).unwrap();
                real.apply_block_with_projections(&x, &projections, &mut over_cached)
                    .unwrap();
                assert_eq!(bits(&plain_cached), bits(&over_cached));
                match &reference {
                    None => reference = Some((bits(&over), bits(&over_cached))),
                    Some((a, b)) => {
                        assert_eq!(a, &bits(&over));
                        assert_eq!(b, &bits(&over_cached));
                    }
                }
                let mut bad = x.clone();
                bad.as_mut_slice()[n / 3] = f64::INFINITY;
                let err = real.apply_block(&bad, &mut over).unwrap_err().to_string();
                assert_eq!(err, dft_physics::RealGammaError::NonFinite.to_string());
                let mut bad_p = projections.clone();
                bad_p[0] = f64::NAN;
                assert!(
                    real.apply_block_with_projections(&x, &bad_p, &mut over)
                        .is_err()
                );
            }
        }
    }
}
/// HAPPLY-FOLD (`BAANDI_HAPPLY_FOLD=on`): even blocks are bit-identical to the
/// pair path; the odd tail column agrees with the zero-partner pair and with
/// the complex H to relative 1e-13 (rounding class); output bits do not depend
/// on the pool budget, and fold plus overlap equals fold alone bit for bit.
#[test]
fn apply_fold_matches_pairs() {
    let mut active = 0;
    let mut max_err = 0.0_f64;
    for rows in [
        [[6., 0.4, 0.2], [1., 7., 0.3], [0.5, 0.7, 8.]],
        [[4., 0., 0.], [0., 5., 0.], [0., 0., 14.]],
        [[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]],
        [[7., 0., 0.], [0., 8., 0.], [0.3, 0., 18.]],
    ] {
        let (basis, h) = fixture(rows, 4, 1.2);
        let n = basis.vectors().len();
        for m in [1, 2, 3, 6, 7] {
            let x = states(n, m, 9 + m);
            let mut reference: Option<Vec<u64>> = None;
            for threads in [1, 2, 4] {
                let mut real =
                    RealHamiltonian::from_complex(&h, &basis, m, WorkerPool::new(threads).unwrap())
                        .unwrap();
                real.set_apply_overlap(false);
                let mut pair = RealMatrix::zeros(n, m).unwrap();
                real.set_apply_fold(false).unwrap();
                real.apply_block(&x, &mut pair).unwrap();
                if !real.set_apply_fold(true).unwrap() {
                    assert!(!h.fft_grid().dims()[2].is_multiple_of(2));
                    continue;
                }
                active += 1;
                let mut fold = RealMatrix::zeros(n, m).unwrap();
                real.apply_block(&x, &mut fold).unwrap();
                if m % 2 == 0 {
                    assert_eq!(bits(&pair), bits(&fold));
                } else {
                    // Leading pairs are untouched; only the tail column moves.
                    assert_eq!(bits(&pair)[..(m - 1) * n], bits(&fold)[..(m - 1) * n]);
                    let tail = |a: &RealMatrix| a.as_slice()[(m - 1) * n..].to_vec();
                    let (p, f) = (tail(&pair), tail(&fold));
                    let norm = p.iter().map(|v| v * v).sum::<f64>().sqrt();
                    let err = p
                        .iter()
                        .zip(&f)
                        .map(|(a, b)| (a - b) * (a - b))
                        .sum::<f64>()
                        .sqrt()
                        / norm;
                    max_err = max_err.max(err);
                    assert!(err < 1e-13, "m={m}: relative {err:e}");
                }
                let cx = complete(real.half_sphere(), &x);
                let rhx = complete(real.half_sphere(), &fold);
                let mut chx = ComplexMatrix::zeros(n, m).unwrap();
                h.apply_block(&cx, &mut chx).unwrap();
                assert!(relative(rhx.as_slice(), chx.as_slice()) < 1e-13);
                real.set_apply_overlap(true);
                let mut both = RealMatrix::zeros(n, m).unwrap();
                real.apply_block(&x, &mut both).unwrap();
                assert_eq!(bits(&fold), bits(&both));
                match &reference {
                    None => reference = Some(bits(&fold)),
                    Some(r) => assert_eq!(r, &bits(&fold), "m={m} T={threads}"),
                }
            }
        }
    }
    println!("fold active cases={active}, max tail relative={max_err:e}");
    assert!(active > 0);
}
#[test]
fn empty_nonlocal_and_invalid_blocks() {
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0., 5., 0.], [0., 0., 5.]]).unwrap();
    let b = PlaneWaveBasis::new(&cell, Hartree(1.)).unwrap();
    let g = FftGridSpec::new(&b, &b.density_basis(&cell).unwrap(), None).unwrap();
    let nl = NonlocalOperator::new(&b, cell.volume_bohr3(), &[]).unwrap();
    let h = Hamiltonian::new(&cell, &b, Array3::zeros(g.dims()), nl).unwrap();
    let mut r = RealHamiltonian::from_complex(&h, &b, 3, WorkerPool::new(1).unwrap()).unwrap();
    let x = states(b.vectors().len(), 3, 7);
    let mut y = x.clone();
    let ptr = y.as_slice().as_ptr();
    for _ in 0..3 {
        r.apply_block(&x, &mut y).unwrap();
    }
    assert_eq!(ptr, y.as_slice().as_ptr());
    let cx = complete(r.half_sphere(), &x);
    let cy = complete(r.half_sphere(), &y);
    let mut expected = cy.clone();
    h.apply_block(&cx, &mut expected).unwrap();
    assert!(relative(cy.as_slice(), expected.as_slice()) < 1e-13);
    let mut bad = x.clone();
    bad.as_mut_slice()[0] = f64::NAN;
    assert!(r.apply_block(&bad, &mut y).is_err());
    assert!(
        r.apply_block(&RealMatrix::zeros(x.rows(), 4).unwrap(), &mut y)
            .is_err()
    );
    let mut empty = RealMatrix::zeros(x.rows(), 0).unwrap();
    r.apply_block(&empty.clone(), &mut empty).unwrap();
}

/// Release only: large setup and 30 alternating pairs exceed the debug budget.
#[test]
#[ignore]
fn timing_blocks() {
    use std::{hint::black_box, time::Instant};
    for (atoms, length, cutoff, bands) in [(4, 8., 2., 9), (8, 12., 3.5, 20), (20, 20., 4., 41)] {
        let (b, h) = fixture(
            [[length, 0., 0.], [0.5, length, 0.], [0., 0.3, length]],
            atoms,
            cutoff,
        );
        let mut r =
            RealHamiltonian::from_complex(&h, &b, bands, WorkerPool::new(1).unwrap()).unwrap();
        let x = states(b.vectors().len(), bands, 11);
        let mut y = x.clone();
        let cx = complete(r.half_sphere(), &x);
        let mut cy = cx.clone();
        r.apply_block(&x, &mut y).unwrap();
        h.apply_block(&cx, &mut cy).unwrap();
        let mut rt = Vec::new();
        let mut ct = Vec::new();
        for i in 0..30 {
            for real in if i % 2 == 0 {
                [true, false]
            } else {
                [false, true]
            } {
                let start = Instant::now();
                if real {
                    r.apply_block(black_box(&x), black_box(&mut y)).unwrap();
                    rt.push(start.elapsed().as_secs_f64());
                } else {
                    h.apply_block(black_box(&cx), black_box(&mut cy)).unwrap();
                    ct.push(start.elapsed().as_secs_f64());
                }
            }
        }
        println!(
            "atoms={atoms} rows={} bands={bands} grid={:?} real_s={rt:?} complex_s={ct:?}",
            b.vectors().len(),
            h.fft_grid().dims()
        );
    }
}

/// HAPPLY microbenchmark (Mac observation only; dev08 decides). Widths 1, 2,
/// 16, 68 at T=1 and T=4; variants base / overlap / fold / fold+overlap in
/// rotating order, 30 rounds, median microseconds per apply call. Fixture:
/// cubic 16 Bohr, 32 atoms with 8 projectors each (s,s,p,p), Ec 15.6 Ha.
/// cargo test --release -p dft-physics --test real_gamma_hamiltonian happly_widths -- --ignored --nocapture --test-threads=1
#[test]
#[ignore]
fn happly_widths() {
    use std::{hint::black_box, time::Instant};
    let length = 16.;
    let (b, h) = fixture_with(
        [[length, 0., 0.], [0., length, 0.], [0., 0., length]],
        32,
        15.6,
        &[0, 0, 1, 1],
        &[0, 0, 1, 1],
    );
    let n = b.vectors().len();
    println!("rows={n} grid={:?}", h.fft_grid().dims());
    let names = ["base", "overlap", "fold", "fold+overlap"];
    for threads in [1, 4] {
        let mut r =
            RealHamiltonian::from_complex(&h, &b, 68, WorkerPool::new(threads).unwrap()).unwrap();
        println!("T={threads} projectors={}", r.projector_count());
        for width in [1, 2, 16, 68] {
            let x = states(n, width, 3);
            let mut y = x.clone();
            let mut samples = vec![Vec::new(); 4];
            let rounds = if width >= 16 { 12 } else { 30 };
            for round in 0..rounds {
                for k in 0..4 {
                    let v = (k + round) % 4;
                    r.set_apply_overlap(v % 2 == 1);
                    r.set_apply_fold(v >= 2).unwrap();
                    r.apply_block(&x, &mut y).unwrap();
                    let start = Instant::now();
                    r.apply_block(black_box(&x), black_box(&mut y)).unwrap();
                    samples[v].push(start.elapsed().as_secs_f64() * 1e6);
                }
            }
            let median = |s: &mut Vec<f64>| {
                s.sort_by(f64::total_cmp);
                s[s.len() / 2]
            };
            let base = median(&mut samples[0]);
            let mut line = format!("T={threads} width={width} base={base:.1}us");
            for v in 1..4 {
                let t = median(&mut samples[v]);
                line += &format!(" {}={t:.1}us ({:.3})", names[v], t / base);
            }
            println!("{line}");
        }
    }
}

/// Exercise the actual parallel GEMM branch: >4 fixed panels and >1,048,576 MACs.
/// Run with cargo test --release -p dft-physics --test real_gamma_hamiltonian large_parallel -- --ignored --nocapture.
#[test]
#[ignore]
fn large_parallel() {
    let (b, h) = fixture([[20., 0., 0.], [0.5, 20., 0.], [0., 0.3, 20.]], 20, 4.);
    let n = b.vectors().len();
    let bands = 17;
    assert!(n > 1024);
    let channels: usize = h.nonlocal().atoms().iter().map(|a| a.columns().len()).sum();
    assert!(n * channels * bands > 1_048_576);
    let x = states(n, bands, 19);
    let mut y = x.clone();
    let mut expected = Vec::new();
    for threads in [1, 2, 4, 8, 12] {
        let mut r = RealHamiltonian::from_complex(&h, &b, bands, WorkerPool::new(threads).unwrap())
            .unwrap();
        r.apply_block(&x, &mut y).unwrap();
        let bits: Vec<_> = y.as_slice().iter().map(|v| v.to_bits()).collect();
        if threads == 1 {
            expected = bits;
        } else {
            assert_eq!(bits, expected);
        }
        let cx = complete(r.half_sphere(), &x);
        let ry = complete(r.half_sphere(), &y);
        let mut cy = cx.clone();
        h.apply_block(&cx, &mut cy).unwrap();
        let error = relative(ry.as_slice(), cy.as_slice());
        assert!(error < 1e-13);
        println!(
            "large rows={n} channels={channels} bands={bands} threads={threads} relative={error:e}"
        );
    }
}

#[test]
fn persistent_hamiltonian_refreshes_potential_and_rejects_replacement() {
    let rows = [[5., 0., 0.], [0.3, 6., 0.], [0.1, 0.2, 7.]];
    let (basis, mut h) = fixture(rows, 2, 1.2);
    let pool = WorkerPool::new(2).unwrap();
    let mut cached = RealHamiltonian::from_complex(&h, &basis, 5, pool.clone()).unwrap();
    let x = states(basis.vectors().len(), 3, 71);
    let mut actual = RealMatrix::zeros(x.rows(), x.cols()).unwrap();
    let mut expected = actual.clone();
    cached.apply_block(&x, &mut actual).unwrap();
    let old = actual.clone();
    for offset in [0.125, -0.37, 0.0] {
        let v = Array3::from_shape_fn(h.fft_grid().dims(), |(i, j, k)| {
            offset + 0.2 * ((i + 3 * j + 7 * k) as f64).cos()
        });
        h.set_local_potential(v).unwrap();
        assert!(cached.update_from_complex(&h, &basis));
        let mut fresh = RealHamiltonian::from_complex(&h, &basis, 5, pool.clone()).unwrap();
        cached.apply_block(&x, &mut actual).unwrap();
        fresh.apply_block(&x, &mut expected).unwrap();
        assert_eq!(
            actual
                .as_slice()
                .iter()
                .map(|v| v.to_bits())
                .collect::<Vec<_>>(),
            expected
                .as_slice()
                .iter()
                .map(|v| v.to_bits())
                .collect::<Vec<_>>()
        );
        assert_ne!(actual, old);
    }
    // Even equal reconstructed data have a different identity, preventing stale
    // reuse for new cells/projectors and allocator address recycling.
    for (rows, atoms, cutoff) in [
        (rows, 2, 1.2),
        (rows, 3, 1.2),
        (rows, 2, 1.5),
        ([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]], 2, 1.2),
    ] {
        let (b, replacement) = fixture(rows, atoms, cutoff);
        assert!(!cached.update_from_complex(&replacement, &b));
    }
}
