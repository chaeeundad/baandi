#!/usr/bin/env bash
# Single entry point for the acceptance gates. Usage: scripts/check.sh [crate]
# Prints one PASS/FAIL line per gate and exits nonzero on any failure.
set -uo pipefail
export PATH="/opt/homebrew/opt/rustup/bin:$PATH"
pkg=(); [[ -n "${1:-}" ]] && pkg=(-p "$1") || pkg=(--workspace)
fail=0
run() { local name="$1"; shift; if "$@" >/tmp/check.$$.log 2>&1; then echo "PASS $name"; else echo "FAIL $name"; tail -20 /tmp/check.$$.log; fail=1; fi; }
run build   cargo build "${pkg[@]}"
run test    cargo test "${pkg[@]}"
run fmt     cargo fmt --check
run clippy  cargo clippy "${pkg[@]}" --all-targets -- -D warnings
rm -f /tmp/check.$$.log
exit $fail
