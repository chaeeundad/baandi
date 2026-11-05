//! Bitwise scalar-versus-block oracle with real projector construction on
//! synthetic NC data: one to three radial projectors for each of s and p.
#[path = "../../dft-upf/tests/common/mod.rs"]
#[allow(dead_code)]
mod synthetic;
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_geometry::Cell;
use dft_numerics::linalg::ComplexMatrix;
use dft_physics::hamiltonian::{Hamiltonian, NonlocalOperator, Species};
use dft_solvers::LinearOperator;
use dft_units::Hartree;
use dft_upf::parse_upf;
use ndarray::Array3;
use num_complex::Complex64;

#[test]
fn physical_projector_variants_and_cell_shapes_keep_bits() {
    for (rows, atoms, multiplicity) in [
        ([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]], 2, 1),
        ([[4., 0., 0.], [0.8, 5., 0.], [0.3, 0.5, 6.]], 8, 2),
        ([[3., 0., 0.], [0., 3., 0.], [0., 0., 12.]], 20, 3),
    ] {
        let n = 81;
        let (r, _) = synthetic::mesh(n, 12.0);
        let mut text = synthetic::fixture_with_points(n);
        let np = 2 * multiplicity;
        let mut projectors = String::from("<PP_NONLOCAL>");
        for i in 0..np {
            let l = i / multiplicity;
            let values: Vec<_> = r
                .iter()
                .map(|r| r.powi(l + 1) * (-r * r / (i + 1) as f64).exp())
                .collect();
            let index = i + 1;
            projectors.push_str(&format!("<PP_BETA.{index} index=\"{index}\" angular_momentum=\"{l}\" cutoff_radius=\"12\" cutoff_radius_index=\"{n}\">{}</PP_BETA.{index}>", synthetic::array(&values)));
        }
        let dij: Vec<_> = (0..np * np)
            .map(|ij| {
                let (i, j) = (ij / np, ij % np);
                if i == j {
                    2.0
                } else if i / multiplicity == j / multiplicity {
                    0.2
                } else {
                    0.0
                }
            })
            .collect();
        projectors.push_str(&format!(
            "<PP_DIJ>{}</PP_DIJ></PP_NONLOCAL>",
            synthetic::array(&dij)
        ));
        let start = text.find("<PP_NONLOCAL>").unwrap();
        let end = text.find("</PP_NONLOCAL>").unwrap() + "</PP_NONLOCAL>".len();
        text.replace_range(start..end, &projectors);
        text = text.replace("number_of_proj=\"3\"", &format!("number_of_proj=\"{np}\""));
        let pp = parse_upf(&text).unwrap();
        let cell = Cell::from_user_rows_bohr(rows).unwrap();
        let basis = PlaneWaveBasis::at_k(&cell, Hartree(1.2), [0.13, 0.07, 0.11]).unwrap();
        let positions: Vec<_> = (0..atoms)
            .map(|a| [a as f64 * 0.17, a as f64 * 0.23, a as f64 * 0.31])
            .collect();
        let grid = FftGridSpec::new(&basis, &basis.density_basis(&cell).unwrap(), None).unwrap();
        let make = || {
            Hamiltonian::new(
                &cell,
                &basis,
                Array3::from_elem(grid.dims(), 0.13),
                NonlocalOperator::new(
                    &basis,
                    cell.volume_bohr3(),
                    &[Species {
                        potential: &pp,
                        positions_bohr: &positions,
                    }],
                )
                .unwrap(),
            )
            .unwrap()
        };
        let scalar = make();
        let n = scalar.dimension();
        let x = ComplexMatrix::new(
            n,
            9,
            (0..n * 9)
                .map(|i| Complex64::new((i as f64).sin(), (i as f64 * 0.7).cos()))
                .collect(),
        )
        .unwrap();
        let mut expected = x.clone();
        for j in 0..9 {
            scalar.apply(x.column(j), expected.column_mut(j)).unwrap();
        }
        for threads in [1, 8] {
            let h = make().with_band_threads(threads).unwrap();
            let mut actual = x.clone();
            h.apply_block(&x, &mut actual).unwrap();
            for (a, b) in actual.as_slice().iter().zip(expected.as_slice()) {
                assert_eq!(
                    (a.re.to_bits(), a.im.to_bits()),
                    (b.re.to_bits(), b.im.to_bits())
                );
            }
        }
    }
}
