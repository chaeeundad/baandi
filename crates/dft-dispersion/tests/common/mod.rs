//! Shared fixtures of the library tests (feature `dftd3`). The library is
//! selected by `BAANDI_DFTD3_LIB`; tests skip with a message when it is unset.
#![allow(dead_code)]

use std::path::PathBuf;

/// One periodic test structure in Bohr.
pub struct Fixture {
    pub name: &'static str,
    pub numbers: Vec<i32>,
    pub positions: Vec<[f64; 3]>,
    pub lattice_rows: [[f64; 3]; 3],
}

impl Fixture {
    pub fn structure(&self) -> dft_dispersion::Structure<'_> {
        dft_dispersion::Structure {
            numbers: &self.numbers,
            positions_bohr: &self.positions,
            lattice_rows_bohr: self.lattice_rows,
        }
    }
    pub fn volume(&self) -> f64 {
        let [a, b, c] = self.lattice_rows;
        let cross = [
            b[1] * c[2] - b[2] * c[1],
            b[2] * c[0] - b[0] * c[2],
            b[0] * c[1] - b[1] * c[0],
        ];
        (a[0] * cross[0] + a[1] * cross[1] + a[2] * cross[2]).abs()
    }
}

fn cartesian(rows: &[[f64; 3]; 3], f: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|c| (0..3).map(|v| f[v] * rows[v][c]).sum())
}

/// AB graphite (Bernal), a = 4.6563 Bohr, c = 12.6820 Bohr, four carbons.
pub fn graphite_ab() -> Fixture {
    let a = 4.6563;
    let c = 12.682;
    let rows = [
        [a, 0., 0.],
        [-0.5 * a, 0.5 * 3f64.sqrt() * a, 0.],
        [0., 0., c],
    ];
    let frac = [
        [0., 0., 0.25],
        [0., 0., 0.75],
        [1. / 3., 2. / 3., 0.25],
        [2. / 3., 1. / 3., 0.75],
    ];
    Fixture {
        name: "graphite_ab",
        numbers: vec![6; 4],
        positions: frac.iter().map(|f| cartesian(&rows, *f)).collect(),
        lattice_rows: rows,
    }
}

/// Triclinic (skewed) cell with Si, O, B and N at low-symmetry positions.
pub fn skewed() -> Fixture {
    let rows = [[7.1, 0.3, -0.4], [1.9, 6.8, 0.5], [-0.7, 1.1, 8.3]];
    let frac = [
        [0.03, 0.07, 0.11],
        [0.31, 0.22, 0.19],
        [0.57, 0.61, 0.43],
        [0.12, 0.47, 0.71],
        [0.83, 0.29, 0.62],
    ];
    Fixture {
        name: "skewed_triclinic",
        numbers: vec![14, 8, 8, 5, 7],
        positions: frac.iter().map(|f| cartesian(&rows, *f)).collect(),
        lattice_rows: rows,
    }
}

/// Water and methane in a 16 Bohr cubic box (periodic images included, as in
/// every engine run).
pub fn molecules_in_box() -> Fixture {
    let positions = vec![
        // H2O
        [3.000, 3.000, 3.000],
        [4.430, 3.000, 4.110],
        [1.570, 3.000, 4.110],
        // CH4, displaced so that the molecules interact
        [9.000, 8.500, 7.800],
        [10.190, 9.690, 9.000],
        [7.810, 7.310, 9.000],
        [10.190, 7.310, 6.600],
        [7.810, 9.690, 6.600],
    ];
    Fixture {
        name: "h2o_ch4_box",
        numbers: vec![8, 1, 1, 6, 1, 1, 1, 1],
        positions,
        lattice_rows: [[16., 0., 0.], [0., 16., 0.], [0., 0., 16.]],
    }
}

pub fn all() -> Vec<Fixture> {
    vec![graphite_ab(), skewed(), molecules_in_box()]
}

/// The library under test, or `None` (with a message) when not configured.
pub fn library() -> Option<PathBuf> {
    match std::env::var_os(dft_dispersion::LIBRARY_ENV) {
        Some(p) if !p.is_empty() => Some(PathBuf::from(p)),
        _ => {
            eprintln!(
                "SKIP: {} is not set; run scripts/fetch_dftd3.py and point it at libs-dftd3",
                dft_dispersion::LIBRARY_ENV
            );
            None
        }
    }
}

/// The `s-dftd3` program: `BAANDI_DFTD3_BIN`, else `<lib dir>/../bin/s-dftd3`.
pub fn program(library: &std::path::Path) -> Option<PathBuf> {
    if let Some(p) = std::env::var_os("BAANDI_DFTD3_BIN").filter(|p| !p.is_empty()) {
        return Some(PathBuf::from(p));
    }
    let candidate = library.parent()?.parent()?.join("bin").join("s-dftd3");
    if candidate.exists() {
        Some(candidate)
    } else {
        eprintln!("SKIP: s-dftd3 not found next to the library; set BAANDI_DFTD3_BIN");
        None
    }
}
