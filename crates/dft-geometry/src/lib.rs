//! Lattice columns in Bohr, reciprocal columns in inverse Bohr.
pub mod kpoints;

use nalgebra::{Matrix3, Vector3};
use std::f64::consts::TAU;
use thiserror::Error;

/// Reject |det(A / max_ij |A_ij|)| <= this dimensionless threshold.
pub const MIN_SCALED_DETERMINANT: f64 = 1e-12;

#[derive(Debug, Clone, PartialEq, Error)]
pub enum CellError {
    #[error("lattice contains non-finite entries")]
    NonFinite,
    #[error("lattice is zero or near singular (scaled determinant <= 1e-12)")]
    NearSingular,
    #[error("volume or inverse cannot be represented as finite nonzero f64 values")]
    Unrepresentable,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    pub species: String,
    /// Dimensionless fractional coordinates; not wrapped into [0, 1).
    pub fractional_position: [f64; 3],
}
impl Atom {
    pub fn new(species: impl Into<String>, fractional_position: [f64; 3]) -> Self {
        Self {
            species: species.into(),
            fractional_position,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Cell {
    lattice: Matrix3<f64>,
    inverse: Matrix3<f64>,
    volume: f64,
}
impl Cell {
    /// User/TOML rows each contain one lattice vector, NOT a matrix row of A.
    /// Explicitly places those vectors into the columns of the internal A.
    pub fn from_user_rows_bohr(rows: [[f64; 3]; 3]) -> Result<Self, CellError> {
        let a = Matrix3::from_columns(&rows.map(|v| Vector3::new(v[0], v[1], v[2])));
        Self::from_matrix(a)
    }

    /// Conventional matrix entries indexed [Cartesian component][vector index].
    /// The lattice vectors are the columns of this matrix.
    pub fn from_matrix_columns_bohr(entries: [[f64; 3]; 3]) -> Result<Self, CellError> {
        Self::from_matrix(Matrix3::from_fn(|i, j| entries[i][j]))
    }

    fn from_matrix(lattice: Matrix3<f64>) -> Result<Self, CellError> {
        if !lattice.iter().all(|x| x.is_finite()) {
            return Err(CellError::NonFinite);
        }
        let scale = lattice.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
        if scale == 0.0 {
            return Err(CellError::NearSingular);
        }
        let scaled = lattice / scale;
        let det = scaled.determinant();
        if det.abs() <= MIN_SCALED_DETERMINANT {
            return Err(CellError::NearSingular);
        }
        let inverse = scaled.try_inverse().ok_or(CellError::NearSingular)? / scale;
        let volume = ((det.abs() * scale) * scale) * scale;
        if !volume.is_finite()
            || volume == 0.0
            || !inverse
                .iter()
                .all(|x| x.is_finite() && (TAU * x).is_finite())
        {
            return Err(CellError::Unrepresentable);
        }
        Ok(Self {
            lattice,
            inverse,
            volume,
        })
    }
    /// Matrix entries [Cartesian component][lattice vector], in Bohr.
    pub fn lattice_columns_bohr(&self) -> [[f64; 3]; 3] {
        std::array::from_fn(|i| std::array::from_fn(|j| self.lattice[(i, j)]))
    }
    /// B = 2*pi*A^(-T), entries [Cartesian component][reciprocal vector], Bohr^-1.
    pub fn reciprocal_columns_inverse_bohr(&self) -> [[f64; 3]; 3] {
        std::array::from_fn(|i| std::array::from_fn(|j| TAU * self.inverse[(j, i)]))
    }
    pub fn volume_bohr3(&self) -> f64 {
        self.volume
    }
    pub fn fractional_to_cartesian_bohr(&self, fractional: [f64; 3]) -> [f64; 3] {
        let r = self.lattice * Vector3::from(fractional);
        [r[0], r[1], r[2]]
    }
    pub fn cartesian_bohr_to_fractional(&self, cartesian: [f64; 3]) -> [f64; 3] {
        (self.inverse * Vector3::from(cartesian)).into()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn cell(rows: [[f64; 3]; 3]) -> Cell {
        Cell::from_user_rows_bohr(rows).unwrap()
    }
    #[test]
    fn duality_volumes_and_positions() {
        // Fixture: 5.43 Angstrom converted to Bohr using CODATA 2018.
        // Conversion APIs are independently tested in dft-units.
        let a: f64 = 10.261_212_856_717_933;
        let cases = [
            (
                "cubic",
                [[10., 0., 0.], [0., 10., 0.], [0., 0., 10.]],
                1000.,
            ),
            (
                "fcc",
                [
                    [0., a / 2., a / 2.],
                    [a / 2., 0., a / 2.],
                    [a / 2., a / 2., 0.],
                ],
                a.powi(3) / 4.,
            ),
            (
                "hexagonal",
                [[2., 0., 0.], [-1., 3_f64.sqrt(), 0.], [0., 0., 3.]],
                6. * 3_f64.sqrt(),
            ),
            (
                "skew",
                [[3., 0., 0.], [29., 0.7, 0.], [41., 2.3, 0.2]],
                0.42,
            ),
            ("left", [[0., 10., 0.], [10., 0., 0.], [0., 0., 10.]], 1000.),
        ];
        for (name, rows, volume) in cases {
            let c = cell(rows);
            let b = c.reciprocal_columns_inverse_bohr();
            let mut duality = 0.0_f64;
            for (i, ai) in rows.iter().enumerate() {
                for (j, _) in b[0].iter().enumerate() {
                    let dot: f64 = (0..3).map(|k| ai[k] * b[k][j]).sum();
                    duality = duality.max((dot - if i == j { TAU } else { 0. }).abs());
                }
            }
            let fractional = [-0.31, 0.73, 1.27];
            let back = c.cartesian_bohr_to_fractional(c.fractional_to_cartesian_bohr(fractional));
            let roundtrip = (0..3)
                .map(|i| (back[i] - fractional[i]).abs())
                .fold(0., f64::max);
            let volume_error = (c.volume_bohr3() - volume).abs() / volume;
            println!(
                "{name}: volume={:.16e} Bohr^3 duality_abs={duality:.16e} position_abs={roundtrip:.16e} volume_rel={volume_error:.16e}",
                c.volume_bohr3()
            );
            assert!(duality < 1e-12); // Dimensionless absolute.
            assert!(roundtrip < 1e-12); // Fractional, absolute.
            assert!(volume_error < 1e-12); // Relative.
        }
    }
    #[test]
    fn rows_are_vectors_not_matrix_rows() {
        let c = cell([[2., 1., 0.], [0., 3., 1.], [1., 0., 4.]]);
        assert_eq!(c.fractional_to_cartesian_bohr([1., 0., 0.]), [2., 1., 0.]);
        let matrix = [[2., 0., 1.], [1., 3., 0.], [0., 1., 4.]];
        let from_columns = Cell::from_matrix_columns_bohr(matrix).unwrap();
        assert_eq!(c.lattice_columns_bohr(), matrix);
        assert_eq!(
            c.lattice_columns_bohr(),
            from_columns.lattice_columns_bohr()
        );
        let atom = Atom::new("Si", [0.25; 3]);
        assert_eq!(atom.species, "Si");
        assert_eq!(atom.fractional_position, [0.25; 3]);
    }
    #[test]
    fn invalid_cells() {
        for rows in [
            [[0.; 3]; 3],
            [[1., 0., 0.], [2., 0., 0.], [0., 0., 1.]],
            [[1., 0., 0.], [0., 1., 0.], [0., 0., 1e-13]],
        ] {
            assert_eq!(
                Cell::from_user_rows_bohr(rows).unwrap_err(),
                CellError::NearSingular
            );
        }
        for x in [f64::NAN, f64::INFINITY] {
            assert_eq!(
                Cell::from_user_rows_bohr([[x; 3]; 3]).unwrap_err(),
                CellError::NonFinite
            );
        }
        for x in [1e-200, 1e200] {
            assert_eq!(
                Cell::from_user_rows_bohr([[x, 0., 0.], [0., x, 0.], [0., 0., x]]).unwrap_err(),
                CellError::Unrepresentable
            );
        }
    }
}
