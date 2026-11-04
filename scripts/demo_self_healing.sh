#!/usr/bin/env bash
# Reproduce the frozen three-arm budget-recovery demo with no network or LLM.
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="/opt/homebrew/opt/rustup/bin:$PATH"
export CARGO_NET_OFFLINE=true
export DFT_HEAL_OUTPUT="${DFT_HEAL_OUTPUT:-$PWD/tests/reference/runtime/generated}"
mkdir -p "$DFT_HEAL_OUTPUT"
DFT_HEAL_OUTPUT="$(cd "$DFT_HEAL_OUTPUT" && pwd)"
{
  rustc -Vv
  cargo -V
  uname -sm
  git rev-parse HEAD
  git status --short
} > "$DFT_HEAL_OUTPUT/build.txt"
cargo test --offline --release -p dft-cli --test self_healing three_arm_demo -- --ignored --exact --nocapture --test-threads=1
