//! Handwritten LP64 oneMKL adapter. No library is linked without feature `mkl`.
//! Linux x86-64 is the supported execution target; other targets can type-check.
#[cfg(feature = "mkl")]
mod adapter;
#[cfg(feature = "mkl")]
pub use adapter::*;
#[cfg(feature = "mkl")]
mod keep;

// Shared by the FFI entry and its no-MKL contract tests.
#[cfg(any(feature = "mkl", test))]
fn valid_in_place_len(planned: usize, actual: usize) -> bool {
    planned > 0 && planned == actual
}

#[cfg(test)]
#[path = "../tests/contracts/in_place_contract.rs"]
mod in_place_contract;
