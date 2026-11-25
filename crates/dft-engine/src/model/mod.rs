//! Data and contracts shared by the production and strict-oracle SCF paths.
pub mod error;
pub mod options;
pub mod result;
pub mod spin;
pub mod verify;

pub use error::*;
pub use options::*;
pub use result::*;
pub use spin::{SpinSector, ThermodynamicLedger};
pub use verify::*;
