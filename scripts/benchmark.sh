#!/usr/bin/env bash
# macOS measurement helper, not an acceptance gate. Outputs override TOML directory.
set -uo pipefail
if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "usage: $0 input.toml [output-directory (default: results)]" >&2
  exit 1
fi
export PATH="/opt/homebrew/opt/rustup/bin:$PATH"
input=$1
output=${2:-results}
mkdir -p "$output" || exit 5
report=$(mktemp "$output/benchmark.XXXXXX") || exit 5
printf 'Measurement helper (not a gate). Includes cargo build/startup overhead.\nInput: %s\nOutput: %s\n' "$input" "$output" > "$report"
/usr/bin/time -l cargo run --offline --release -p dft-cli -- run "$input" --output-dir "$output" 2>> "$report"
code=$?
printf '\nExit code: %s\nmacOS maximum resident set size is in bytes; real time is in seconds.\n' "$code" >> "$report"
printf 'Measurement report: %s\n' "$report"
exit "$code"
