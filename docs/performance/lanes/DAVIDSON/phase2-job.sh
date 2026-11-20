#!/usr/bin/env bash
# dev08 only. Each invocation is one case/model/stage and one <=600s job.
set -euo pipefail
cell=${1:?si16|si32|si64}
np=${2:?1|8}
mode=${3:?solve|screen|full|profile}
case "$cell:$np" in si16:8|si32:8|si64:8|si32:1|si64:1) ;; *) exit 2 ;; esac
case "$mode" in solve|screen|full|profile) ;; *) exit 2 ;; esac
: "${BIN:?absolute dft-cli path}" "${BENCH:?absolute davidson_bench path}"
: "${UPF:?absolute approved UPF path}" "${OUT:?absolute output root}"
job="$OUT/$cell-np$np-$mode"
mkdir -p "$job"
python3 -I scripts/fast_screen_dav.py --upf "$UPF" --case "$cell" \
  --seed 601 --amp 0.001 --threads 1 --emit-input "$job/input.toml"
if [[ "$mode" = screen ]]; then
  sed -i 's/max_iterations = 100/max_iterations = 10/' "$job/input.toml"
fi
sha256sum "$BIN" "$BENCH" "$UPF" "$job/input.toml" > "$job/hashes.txt"
uptime > "$job/load-before.txt"
launcher=(taskset -c 0)
if [[ "$np" = 8 ]]; then launcher=(taskset -c 0-7 mpirun -np 8 --bind-to core --map-by core); fi
# AA first, then crossed AB, in the same placement/window. Profile is separate.
phases='aa ab'
pairs=3
profile=off
if [[ "$mode" = profile ]]; then phases=ab; pairs=1; profile=on; fi
for phase in $phases; do
  for ((pair=1; pair<=pairs; pair++)); do
    order='a b'
    if (( pair % 2 == 0 )); then order='b a'; fi
    for arm in $order; do
      switch=off
      if [[ "$phase:$arm" = ab:b ]]; then switch=on; fi
      stem="$job/$phase-$pair-$arm"
      cmd=("$BIN" run "$job/input.toml" --output-dir "$stem")
      if [[ "$mode" = solve ]]; then
        cmd=("$BENCH" --input "$job/input.toml" --repeats 1 --kinds warm)
      fi
      set +e
      /usr/bin/time -f '%e' -o "$stem.external-s" env \
        OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
        BAANDI_DAVIDSON_UPPER_BORDER="$switch" BAANDI_DAVIDSON_PROFILE="$profile" \
        BAANDI_PHASE_PROFILE="$profile" BAANDI_STAGE_TIMERS="$profile" \
        BAANDI_DAVIDSON_STAGES="$profile" BAANDI_DGEMM_TIMERS=off \
        "${launcher[@]}" "${cmd[@]}" > "$stem.stdout" 2> "$stem.stderr"
      rc=$?
      set -e
      echo "$rc" > "$stem.exit-code"
      if [[ "$rc" != 0 && !( "$mode" = screen && "$rc" = 3 ) ]]; then exit "$rc"; fi
    done
  done
done
uptime > "$job/load-after.txt"
