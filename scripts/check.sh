#!/usr/bin/env bash
# Single entry point for the acceptance gates. Usage: scripts/check.sh [crate]
# Prints one PASS/FAIL line per gate (build, test, fmt, clippy) and exits nonzero on any failure.
#
# Speed (docs/performance/slow-tests.md): tests are compiled once, then the test
# binaries run concurrently (scripts/run_test_binaries.py, CHECK_TEST_JOBS workers,
# default 4) while clippy runs; doctests run with `cargo test --doc`. Coverage is the
# same as `cargo test --no-fail-fast`. CHECK_TEST_JOBS=1 runs the binaries one at a time.
set -uo pipefail
export PATH="/opt/homebrew/opt/rustup/bin:$PATH"
here="$(cd "$(dirname "$0")" && pwd)"
pkg=(); [[ -n "${1:-}" ]] && pkg=(-p "$1") || pkg=(--workspace)
jobs="${CHECK_TEST_JOBS:-4}"
tmp="$(mktemp -d "${TMPDIR:-/tmp}/check.XXXXXX")"; trap 'rm -rf "$tmp"' EXIT
fail=0
report() { local name="$1" rc="$2" log="$3"; if [[ "$rc" == 0 ]]; then echo "PASS $name"; else echo "FAIL $name"; tail -40 "$log"; fail=1; fi; }

cargo build "${pkg[@]}" >"$tmp/build.log" 2>&1
report build $? "$tmp/build.log"

target_dir="$(cargo metadata --format-version 1 --no-deps 2>/dev/null \
  | python3 -I -c 'import json,sys; print(json.load(sys.stdin)["target_directory"])' 2>/dev/null)"

# test: compile once, then run the binaries; clippy runs concurrently with them in
# its own target subdirectory, so it never holds the build lock that tests spawning
# `cargo run` (dft-upf parser) need.
test_rc=0; clippy_pid=""
if cargo test "${pkg[@]}" --no-run --message-format=json >"$tmp/test.json" 2>"$tmp/test.log"; then
  clippy_dir="${target_dir:+$target_dir/check-clippy}"
  CARGO_TARGET_DIR="${clippy_dir:-${CARGO_TARGET_DIR:-target}}" \
    cargo clippy "${pkg[@]}" --all-targets -- -D warnings >"$tmp/clippy.log" 2>&1 &
  clippy_pid=$!
  times_args=(); [[ -n "$target_dir" ]] && times_args=(--times "$target_dir/check-test-times.json")
  # Under `cargo test` the rustup proxy passes RUSTUP_TOOLCHAIN to the test binaries; tests
  # that spawn `cargo run` (dft-upf parser) need it, or registry crates are compiled by the
  # machine's default rustc instead of the pinned one (mixed-toolchain errors, rebuilds).
  # Such a `cargo run -p dft-cli` also resolves a different dependency graph than the
  # workspace build; sharing one target dir made each invalidate the other (faer and five
  # workspace crates rebuilt on every gate, about 40 s), so they get their own target dir.
  toolchain="${RUSTUP_TOOLCHAIN:-$(rustup show active-toolchain 2>/dev/null | awk 'NR==1{print $1}')}"
  nested_dir="${target_dir:+$target_dir/check-nested}"
  env ${toolchain:+RUSTUP_TOOLCHAIN="$toolchain"} ${nested_dir:+CARGO_TARGET_DIR="$nested_dir"} python3 -I "$here/run_test_binaries.py" --jobs "$jobs" ${times_args[@]+"${times_args[@]}"} <"$tmp/test.json" >"$tmp/bins.log" 2>&1 || test_rc=1
  cargo test "${pkg[@]}" --doc >"$tmp/doc.log" 2>&1 || { test_rc=1; tail -40 "$tmp/doc.log" >>"$tmp/test.log"; }
  # On failure the report shows the failing binaries' tails and the runner summary last.
  cat "$tmp/bins.log" >>"$tmp/test.log"
else
  test_rc=1
fi
report test $test_rc "$tmp/test.log"

cargo fmt --check >"$tmp/fmt.log" 2>&1
report fmt $? "$tmp/fmt.log"

if [[ -n "$clippy_pid" ]]; then wait "$clippy_pid"; clippy_rc=$?
else cargo clippy "${pkg[@]}" --all-targets -- -D warnings >"$tmp/clippy.log" 2>&1; clippy_rc=$?; fi
report clippy $clippy_rc "$tmp/clippy.log"
exit $fail
