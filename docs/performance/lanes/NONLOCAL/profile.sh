#!/usr/bin/env bash
# dev08 only; one CASE/NP per job; run through the master in its reserved window.
set -euo pipefail
case_name=${1:?si16, si32 or si64}
ranks=${2:?1 or 8}
out=${3:?fresh output directory}
upf=${4:?external Si.pz-vbc.UPF}
case "$case_name:$ranks" in si16:8|si32:8|si64:8|si32:1|si64:1) ;; *) exit 2;; esac
[[ ! -e "$out" ]] || { echo 'Output directory must be fresh' >&2; exit 2; }
mkdir -p "$out"
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
binary=target/release/dft-cli
python3 -I scripts/fast_screen_dav.py --upf "$upf" --case "$case_name" \
  --seed 601 --amp 0.001 --threads 1 --emit-input "$out/input.toml"
sha256sum "$binary" "$upf" "$out/input.toml" > "$out/sha256.txt"
git rev-parse HEAD > "$out/commit.txt"
git diff > "$out/source.patch"
rustc --version > "$out/rustc.txt"
# off:0 is the discarded warm-up. Retained orders are crossed across repeats.
for item in off:0 off:1 counts:1 timed:1 timed:2 counts:2 off:2 counts:3 off:3 timed:3; do
  mode=${item%:*}
  rep=${item#*:}
  case "$mode" in
    off) profile=off; clocks=off;;
    counts) profile=on; clocks=off;;
    timed) profile=on; clocks=on;;
  esac
  prefix="$out/$mode-$rep"
  uptime > "$prefix.load-before.txt"
  launcher=(taskset -c 0)
  if [[ "$ranks" = 8 ]]; then launcher=(mpirun -np 8 --bind-to core --map-by core --report-bindings); fi
  env BAANDI_PHASE_PROFILE="$profile" BAANDI_STAGE_TIMERS="$clocks" \
    "${launcher[@]}" "$binary" run "$out/input.toml" --output-dir "$prefix" \
    > "$prefix.stdout" 2> "$prefix.stderr"
  uptime > "$prefix.load-after.txt"
done
python3 -I scripts/lane_nonlocal_parse.py --machine dev08 \
  "$out"/timed-{1,2,3}/results.json > "$out/nonlocal.csv"
python3 -I scripts/lane_nonlocal_parse.py --machine dev08 --format json \
  "$out"/counts-{1,2,3}/results.json > "$out/nonlocal-counts.json"
