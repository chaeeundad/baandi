//! Handwritten LP64 oneMKL adapter. No library is linked without feature `mkl`.
//! Linux x86-64 is the supported execution target; other targets can type-check.
#[cfg(feature = "mkl")]
mod adapter;
#[cfg(feature = "mkl")]
pub use adapter::*;
#[cfg(feature = "mkl")]
mod keep;
