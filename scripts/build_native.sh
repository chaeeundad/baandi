#!/usr/bin/env bash
# Build dft-cli tuned for the machine it is built on (target-cpu=native). The result is faster
# (dev08, 2026-11-05: pure Rust 1.12 to 1.40x, bit-identical energy and density hashes for the
# tested cases) but runs only on CPUs that have the same instruction set. Benchmarks on dev08 use
# this build in addition to the default build and record both. Default builds stay portable.
# Usage: scripts/build_native.sh [extra cargo args, e.g. --features mkl]
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="/opt/homebrew/opt/rustup/bin:$PATH"
TARGET_DIR="${TARGET_DIR:-target/native}"
RUSTFLAGS="${RUSTFLAGS:-} -C target-cpu=native" cargo build --release -p dft-cli --target-dir "$TARGET_DIR" "$@"
echo "built $TARGET_DIR/release/dft-cli (target-cpu=native)"
sha256sum "$TARGET_DIR/release/dft-cli" 2>/dev/null || shasum -a 256 "$TARGET_DIR/release/dft-cli"
