//! Input parsing and plain-text P1 reports.
pub mod config;

use config::{Config, PreparedInput};
use std::fmt::Write;

pub fn grid_report(config: &Config, input: &PreparedInput) -> String {
    let mut text = String::new();
    writeln!(text, "System: {}", config.system.name).unwrap();
    writeln!(
        text,
        "Cell volume: {:.12} Bohr^3",
        input.cell.volume_bohr3()
    )
    .unwrap();
    writeln!(text, "Reciprocal lattice vectors (Bohr^-1):").unwrap();
    let b = input.cell.reciprocal_columns_inverse_bohr();
    for (i, ((x, y), z)) in b[0].iter().zip(&b[1]).zip(&b[2]).enumerate() {
        writeln!(text, "  b{} = [{:.12}, {:.12}, {:.12}]", i + 1, x, y, z).unwrap();
    }
    writeln!(
        text,
        "Ecut: {:.12} Ry = {:.12} Ha",
        config.basis.ecutwfc_ry,
        input.wavefunctions.cutoff().0
    )
    .unwrap();
    writeln!(
        text,
        "Plane-wave count: {}",
        input.wavefunctions.vectors().len()
    )
    .unwrap();
    writeln!(
        text,
        "Density-sphere count: {}",
        input.density.vectors().len()
    )
    .unwrap();
    let d = input.grid.dims();
    writeln!(
        text,
        "FFT dims: {} x {} x {} ({})",
        d[0],
        d[1],
        d[2],
        if input.grid.overridden() {
            "overridden"
        } else {
            "automatic"
        }
    )
    .unwrap();
    text
}
