#!/usr/bin/env bash
set -euo pipefail
case_name=${1:?case si16, si32, or si64 required}
ranks=${2:?rank count 1 or 8 required}
case "$case_name:$ranks" in si16:8|si32:8|si64:8|si32:1|si64:1) ;; *) exit 2 ;; esac
: "${BIN:?set absolute binary path}" "${UPF:?set approved UPF path}" "${OUT:?set absolute output root}"
job="$OUT/$case_name-np$ranks"
mkdir -p "$job"
python3 -I scripts/fast_screen_dav.py --upf "$UPF" --case "$case_name" \
  --seed 601 --amp 0.001 --threads 1 --emit-input "$job/input.toml"
sha256sum "$BIN" "$UPF" "$job/input.toml" > "$job/hashes.txt"
uptime > "$job/load-before.txt"
for pair in 1 2 3; do
  order="off on"
  if [[ "$pair" = 2 ]]; then order="on off"; fi
  for arm in $order; do
    dest="$job/$arm-$pair"
    launcher=(taskset -c 0)
    if [[ "$ranks" = 8 ]]; then
      launcher=(taskset -c 0-7 mpirun -np 8 --bind-to core --map-by core)
    fi
    set +e
    env MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 BAANDI_DAVIDSON_PROFILE="$arm" \
      BAANDI_PHASE_PROFILE=on BAANDI_STAGE_TIMERS=on BAANDI_DAVIDSON_STAGES=on \
      "${launcher[@]}" "$BIN" run "$job/input.toml" --output-dir "$dest" \
      > "$job/$arm-$pair.stdout" 2> "$job/$arm-$pair.stderr"
    rc=$?
    set -e
    echo "$rc" > "$job/$arm-$pair.exit-code"
    if [[ "$rc" != 0 ]]; then exit "$rc"; fi
  done
  python3 -I docs/performance/lanes/DAVIDSON/verify_identity.py \
    "$job/off-$pair" "$job/on-$pair" > "$job/identity-$pair.json"
done
uptime > "$job/load-after.txt"
python3 -I scripts/lane_davidson_parse.py "$job/"on-*/results.json \
  --json-out "$job/tables.json" > "$job/tables.md"
