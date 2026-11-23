#!/usr/bin/env bash
# Single entry point for the acceptance gates. Prints one PASS/FAIL line per gate
# (build, test, fmt, clippy) and exits nonzero on any failure.
#
#   scripts/check.sh fast [crate...]  iteration tier: the crates changed against the
#                                     merge base (scripts/changed_crates.py) and their
#                                     direct dependents, or the crates named
#   scripts/check.sh full             pre-push gate: the whole workspace (same as no argument)
#   scripts/check.sh <crate>          one crate, every gate (older form, kept)
#
# Both tiers run the same gates, the same tests (tests ignored for runtime are skipped
# in both; they run in release, docs/performance/slow-tests.md) and the same lints;
# `fast` only narrows the crate set. Protocol: docs/performance/fast-loop-protocol.md
# ("Two lanes").
#
# Speed (docs/performance/slow-tests.md): tests are compiled once, then the test
# binaries run concurrently (scripts/run_test_binaries.py, CHECK_TEST_JOBS workers,
# default 4) while clippy runs; doctests run with `cargo test --doc`. Coverage is the
# same as `cargo test --no-fail-fast`. CHECK_TEST_JOBS=1 runs the binaries one at a time.
#
# Environment:
#   CHECK_BASE=<rev>   base of the `fast` diff (default: nearer merge base with main / origin/main)
#   CHECK_NICE=1       run the test phase under `nice -n 10` (default 1; 0 = no nice), so a
#                      gate does not starve other sessions on a shared machine
#   CHECK_TEST_JOBS=4  concurrent test binaries
#   TMPDIR_990=<dir>   scratch for logs (default /Volumes/990-2TB/Develop/rust-dft-scratch/tmp
#                      when that volume exists; user rule 2026-11-19: no scratch on the
#                      internal disk)
#   BAANDI_TEST_TMP=<dir>  writable absolute test scratch root OUTSIDE any checkout;
#                      defaults to TMPDIR_990 (or its external-volume default).
#                      A private gate directory under it is exported as
#                      BAANDI_TEST_TMP and TMPDIR (also for the source-local test).
#                      For direct cargo tests, set both variables to a private root.
#                      A sandbox may use $PWD/target/tmp, outside Cargo packages,
#                      but cache success tests will still reject that checkout path.
set -uo pipefail
export PATH="/opt/homebrew/opt/rustup/bin:$PATH"
here="$(cd "$(dirname "$0")" && pwd)"
scratch="${BAANDI_TEST_TMP:-${TMPDIR_990:-/Volumes/990-2TB/Develop/rust-dft-scratch/tmp}}"
# Fail rather than silently writing test scratch to the internal disk.
mkdir -p "$scratch" || exit 1
export BAANDI_TEST_TMP="$scratch" TMPDIR="$scratch"
# sccache (user-approved 2026-11-19): shared compile cache on the 990-2TB volume, safe across worktrees
# because it keys on the full compiler input. CHECK_SCCACHE=0 disables. sccache cannot cache incremental
# crates, so incremental is switched off for gate builds when it is used.
if [[ "${CHECK_SCCACHE:-1}" != "0" ]] && command -v sccache >/dev/null 2>&1; then
  export SCCACHE_DIR="${SCCACHE_DIR:-/Volumes/990-2TB/Develop/rust-dft-cache/sccache}"
  mkdir -p "$SCCACHE_DIR" 2>/dev/null && export RUSTC_WRAPPER=sccache CARGO_INCREMENTAL=0
fi

tier=full; pkg=(--workspace); fmt_pkg=()
case "${1:-}" in
  ""|full) tier=full ;;
  fast)
    tier=fast; shift
    if [[ $# -gt 0 ]]; then crates=("$@")
    else crates=(); while IFS= read -r c; do [[ -n "$c" ]] && crates+=("$c"); done \
      < <(python3 -I "$here/changed_crates.py" --explain 2>"${TMPDIR:-/tmp}/check-select.$$"); fi
    sed 's/^/  /' "${TMPDIR:-/tmp}/check-select.$$" 2>/dev/null; rm -f "${TMPDIR:-/tmp}/check-select.$$"
    if [[ ${#crates[@]} -eq 1 && "${crates[0]}" == WORKSPACE ]]; then pkg=(--workspace)
    elif [[ ${#crates[@]} -eq 0 ]]; then pkg=()
    else pkg=(); for c in "${crates[@]}"; do pkg+=(-p "$c"); done; fmt_pkg=("${pkg[@]}"); fi ;;
  *) tier=crate; pkg=(-p "$1"); fmt_pkg=(-p "$1") ;;
esac
jobs="${CHECK_TEST_JOBS:-4}"
nice_cmd=(); [[ "${CHECK_NICE:-1}" != 0 ]] && nice_cmd=(nice -n 10)

# One-line load report: a gate wall measured under load is an upper bound, not a regression.
python3 -I -c '
import os, subprocess, sys
l = os.getloadavg()
p = subprocess.run(["ps", "-A", "-o", "comm="], stdout=subprocess.PIPE, text=True).stdout.split()
n = sum(1 for c in p if os.path.basename(c) in ("rustc", "cargo", "clippy-driver", "rustdoc"))
what = sys.argv[2] or "no crate affected"
print(f"LOAD {l[0]:.1f} {l[1]:.1f} {l[2]:.1f} on {os.cpu_count()} cpus; {n} cargo/rustc processes running; "
      f"tier {sys.argv[1]}: {what}")
' "$tier" "${pkg[*]:-}"

start=$(date +%s)
tmp="$(mktemp -d "${TMPDIR:-/tmp}/check.XXXXXX")" || exit 1
trap 'rm -rf "$tmp"' EXIT
# Isolate cache fixtures per gate, including the source-local fetch test that
# uses std::env::temp_dir(). The existing trap cleans up this private directory.
export BAANDI_TEST_TMP="$tmp" TMPDIR="$tmp"
fail=0
report() { local name="$1" rc="$2" log="$3"; if [[ "$rc" == 0 ]]; then echo "PASS $name"; else echo "FAIL $name"; tail -40 "$log"; fail=1; fi; }

if [[ ${#pkg[@]} -eq 0 ]]; then
  # fast tier, no crate affected (docs, scripts, tickets): only the format gate has work.
  echo "PASS build (no crate affected)"; echo "PASS test (no crate affected)"
  cargo fmt --check >"$tmp/fmt.log" 2>&1; report fmt $? "$tmp/fmt.log"
  echo "PASS clippy (no crate affected)"
  echo "WALL $(( $(date +%s) - start )) s"
  exit $fail
fi

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
  # The dft-upf `parser` test runs the CLI: when this gate built dft-cli it gets that
  # binary (DFT_CLI_BIN) and no nested build happens; otherwise it falls back to
  # `cargo run -p dft-cli` in check-nested.
  toolchain="${RUSTUP_TOOLCHAIN:-$(rustup show active-toolchain 2>/dev/null | awk 'NR==1{print $1}')}"
  nested_dir="${target_dir:+$target_dir/check-nested}"
  cli_bin=""
  if [[ -n "$target_dir" && -x "$target_dir/debug/dft-cli" ]] && [[ " ${pkg[*]} " == *" --workspace "* || " ${pkg[*]} " == *" dft-cli "* ]]; then
    cli_bin="$target_dir/debug/dft-cli"
  fi
  env ${toolchain:+RUSTUP_TOOLCHAIN="$toolchain"} ${nested_dir:+CARGO_TARGET_DIR="$nested_dir"} ${cli_bin:+DFT_CLI_BIN="$cli_bin"} \
    ${nice_cmd[@]+"${nice_cmd[@]}"} python3 -I "$here/run_test_binaries.py" --jobs "$jobs" \
    ${times_args[@]+"${times_args[@]}"} <"$tmp/test.json" >"$tmp/bins.log" 2>&1 || test_rc=1
  ${nice_cmd[@]+"${nice_cmd[@]}"} cargo test "${pkg[@]}" --doc >"$tmp/doc.log" 2>&1 \
    || { test_rc=1; tail -40 "$tmp/doc.log" >>"$tmp/test.log"; }
  # On failure the report shows the failing binaries' tails and the runner summary last.
  cat "$tmp/bins.log" >>"$tmp/test.log"
else
  test_rc=1
fi
report test $test_rc "$tmp/test.log"
tail -1 "$tmp/bins.log" 2>/dev/null | sed 's/^/  /'

cargo fmt --check ${fmt_pkg[@]+"${fmt_pkg[@]}"} >"$tmp/fmt.log" 2>&1
report fmt $? "$tmp/fmt.log"

if [[ -n "$clippy_pid" ]]; then wait "$clippy_pid"; clippy_rc=$?
else cargo clippy "${pkg[@]}" --all-targets -- -D warnings >"$tmp/clippy.log" 2>&1; clippy_rc=$?; fi
report clippy $clippy_rc "$tmp/clippy.log"
echo "WALL $(( $(date +%s) - start )) s ($tier)"
exit $fail
