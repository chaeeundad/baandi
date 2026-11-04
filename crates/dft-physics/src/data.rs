//! Types used by the physics APIs, re-exported for orchestration clients.
//! This avoids requiring clients to add backend dependencies merely to construct
//! the arrays, FFT plans, and unit values accepted by the existing physics API.
pub use dft_numerics::{
    Fft3d, FftError, FftPlan3d,
    linalg::{ComplexMatrix, LinalgError},
    pool::WorkerPool,
};
pub use dft_units::Hartree;
pub use ndarray::Array3;
pub use num_complex::Complex64;

/// Gaussian Coulomb screening functions, provided by the adopted libm backend.
/// Arguments are dimensionless. These are also useful for analytic fixtures.
pub use libm::{erf, erfc};
