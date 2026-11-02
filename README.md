# rust-dft

Rust-first plane-wave DFT engine bootstrap. P0 contains nine std-only crates
and Apple Accelerate ABI tests; it does not perform DFT calculations.

Requires Apple Silicon macOS, Xcode Command Line Tools, Rust 1.99.1,
rustfmt and clippy. See [environment](docs/environment.md) for the measured
host and offline toolchain setup.

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
export CARGO_NET_OFFLINE=true
export VECLIB_MAXIMUM_THREADS=1
cargo build --workspace
cargo test --workspace
cargo fmt --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test -p dft-numerics -- --nocapture
```

[Architecture](docs/architecture.md), [dependencies](docs/dependencies.md),
[backend decision](docs/decisions/0001-numerical-backend-isolation.md).
No project license has been selected yet.
