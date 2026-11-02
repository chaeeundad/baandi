# P0 dependencies

External Rust crates: **0**. Inter-workspace dependencies: **0**.
All nine packages use std only; no build dependencies or build scripts.
Native test dependency: Apple Accelerate.framework supplied by macOS/CLT.
Toolchain: Rust 1.99.1, rustfmt, clippy; resolver 3, edition 2024.
No registry fetch is necessary. Set CARGO_NET_OFFLINE=true for validation.

Cargo generates Cargo.lock containing only the nine local packages. Cargo.lock
is outside P0 OWNED_FILES, so it is not included in the commit; a later ticket
should authorize tracking it before external dependencies are added.
FFTW/rustfft are candidates only (unverified); see architecture.md.
Apple SDK terms and eventual dependency licenses need review before distribution.
The project LICENSE is deliberately undecided.
