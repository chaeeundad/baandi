# P0 measured environment

Observed 2026-11-02, without network lookup or installation.

| Command / probe | Actual result |
| --- | --- |
| `uname -m` | arm64 |
| `sw_vers` | macOS 27.0.1, build 26A434 |
| `rustc --version` | rustc 1.99.1 (58cf8f9c3 2026-10-13) |
| `cargo --version` before pin | cargo 1.99.1 (5f94df478 2026-08-27) |
| `rustup show` | stable-aarch64-apple-darwin only; target aarch64-apple-darwin |
| `rustup component list --installed` | cargo, clippy, rust-std, rustc, rustfmt |
| `xcode-select -p` | /Library/Developer/CommandLineTools |
| `xcrun --show-sdk-path` | /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk |
| `xcrun clang --version` | Apple clang 21.0.0 (clang-2100.3.34.2), arm64-apple-darwin27.0.0 |
| SDK Accelerate.framework directory listing | Accelerate.tbd, Frameworks, Headers, Modules, Versions |
| `command -v pw.x`, `command -v fftw-wisdom` | no PATH entry |
| `brew list --versions fftw openblas quantum-espresso` | no output |
| `pkg-config --modversion fftw3` / `openblas` | pkg-config command not found (127) |
| `printenv VECLIB_MAXIMUM_THREADS` / `OMP_NUM_THREADS` | unset (1) |

These discovery probes do not exclude software installed outside PATH/Homebrew.
No QE reference calculation was performed.

## Exact version and restricted sandbox

rust-toolchain.toml pins 1.99.1. Initially `cargo --version` with that pin exited
1: rustup tried syncing and could not create a file under ~/.rustup/tmp
(Operation not permitted). The installed stable directory contains exactly
1.99.1. For this sandbox only, a temporary rustup home aliases that verified
installation under its exact version name; no downloads or global edits:

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
# Verify the source version BEFORE aliasing; do not reuse after a stable update.
/Users/mino/.rustup/toolchains/stable-aarch64-apple-darwin/bin/rustc --version
mkdir -p /private/tmp/rust-dft-p0-rustup/toolchains
# Run once, in a fresh temporary home:
ln -s /Users/mino/.rustup/toolchains/stable-aarch64-apple-darwin \
  /private/tmp/rust-dft-p0-rustup/toolchains/1.99.1-aarch64-apple-darwin
export RUSTUP_HOME=/private/tmp/rust-dft-p0-rustup
export CARGO_NET_OFFLINE=true
export VECLIB_MAXIMUM_THREADS=1
cargo build --workspace
cargo test --workspace
cargo fmt --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test -p dft-numerics -- --nocapture
```

On an unrestricted host, provision the genuine 1.99.1 rustup toolchain and its
rustfmt/clippy components normally. The temporary symlink is not a portable
installation or immutable compiler archive.

## Native ABI contract

Declarations were inspected in SDK vecLib Headers/cblas.h and clapack.h.
The legacy `cblas_ddot` uses C int dimensions/strides by value; `dsyev_` uses
32-bit integer pointers, double buffers, character pointers, and returns C int.
INFO is the LAPACK status; no hidden Fortran character lengths occur in this
SDK C declaration. Rust uses i32/f64/c_char accordingly. This is LP64 (64-bit
pointers, 32-bit BLAS integers), not ILP64. No ILP64 symbol or modern
ACCELERATE_NEW_LAPACK interface is used. SDK headers mark these legacy APIs
deprecated; migration requires a separate ABI test.

The matrix is column-major with lda=2. Upper-triangle selection is tested by
poisoning the unused lower triangle with NaN. Workspace is 3*n-1=5 doubles;
jobz=V checks eigenvectors as well as eigenvalues. Storage is fixed-size,
aligned, initialized where read, non-overlapping and alive across calls.
Unsafe declarations/calls exist only in the private adapter.

P0 uses real ddot/dsyev. Complex calling conventions are **not tested**.
For future zgemm/zheev, verify the SDK's two-double real/imag layout against
an explicit repr(C) type and pointer-based scalar/buffer arguments. Never
assume Rust tuple layout or a complex return ABI. P1 needs a genuinely complex
Hermitian test before production Γ-only complex calculations.

Thread control is process configuration, not a library-side environment edit.
Validation below requests VECLIB_MAXIMUM_THREADS=1 before process startup;
actual worker counts and bitwise reproducibility across OS versions are not
measured. An initial successful smoke run also used the default unset value.

## Numerical evidence

Dimensionless synthetic inputs, not energies or QE reference data:
- ddot([1,2,3],[4,5,6]) = 32.0000000000000000, exact expected 32.
- A=[[2,1],[1,2]]: eigenvalues [1.0,3.0], analytic 2±1; INFO=0.
- max |Av−λv| = 0.0000000000000000e0.
- max |VᵀV−I| = 2.2204460492503131e-16.
- Absolute eigenvalue/residual/orthogonality tolerance: 1e-12, unchanged.

## Validation command exits

All commands used the PATH, temporary RUSTUP_HOME, CARGO_NET_OFFLINE=true
and VECLIB_MAXIMUM_THREADS=1 settings above.

| Command | Exit code |
| --- | --- |
| `cargo build --workspace` | 0 |
| `cargo test --workspace` | 0 (10 unit tests passed; 0 doctests) |
| `cargo fmt --check` | 0 |
| `cargo clippy --workspace --all-targets -- -D warnings` | first 101 (lint group priority), then 0 after setting group priority=-1 |
| `cargo test -p dft-numerics -- --nocapture` | 0 (2 native tests, values above) |

## P1A: Accelerate FFI removed

Per ADR 0004, P1A removes `dft-numerics/src/accelerate.rs` and its native ABI
tests. The P0 ABI discussion and measurements above are historical, not the
current implementation. The foundation now uses safe Rust nalgebra geometry
and reusable ndrustfft handlers over rustfft, with no Accelerate linkage in
workspace source. FFT normalization and tolerances are specified in
`physics-conventions.md`. Offline checks use the verified Rust 1.99.1 temporary
rustup-home procedure above; no global toolchain changes are needed.
