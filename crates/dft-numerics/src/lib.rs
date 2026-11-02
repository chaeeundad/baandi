//! Numerical backend boundary. P0 validates only the native Accelerate ABI.

#[cfg(all(target_os = "macos", target_arch = "aarch64"))]
#[allow(unsafe_code)]
mod accelerate;

#[cfg(not(all(target_os = "macos", target_arch = "aarch64")))]
compile_error!("P0 Accelerate validation requires aarch64-apple-darwin");
