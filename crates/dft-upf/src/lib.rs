//! Validated UPF v2 NC data and radial transforms in Hartree atomic units.
mod parser;
pub mod radial;
pub use parser::{Projector, Pseudopotential, UpfError, from_file, parse_upf};
