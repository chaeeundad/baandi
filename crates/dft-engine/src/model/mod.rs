//! Data and contracts shared by the production and strict-oracle SCF paths.
pub mod error;
pub mod input;
pub mod options;
pub mod result;
pub mod runtime;
pub mod spin;
pub mod verify;

pub use error::*;
pub use input::{PreparedSymmetry, SpeciesPotential};
pub use options::*;
pub use result::*;
pub use spin::{SpinSector, ThermodynamicLedger};
pub use verify::*;

pub mod collinear;
pub mod finite_temperature;
