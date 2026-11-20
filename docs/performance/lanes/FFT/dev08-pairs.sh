#!/usr/bin/env bash
# dev08 only; syntax-checked on Mac. One case/model per <=600 s queue job.
set -euo pipefail
CASE=${1:?si16, si32 or si64}
MODE=${2:?solve, steps or full}
KIND=${3:?counts, profile or time}
case "$CASE" in si16|si32|si64) ;; *) exit 2;; esac
case "$MODE" in solve|steps|full) ;; *) exit 2;; esac
case "$KIND" in counts|profile|time) ;; *) exit 2;; esac
: "${DFT_SI_UPF:?same external UPF as baseline}"
: "${FFT_WORK:?fresh absolute output directory}"
: "${FFT_TARGET:?absolute target directory containing MPI+MKL binaries}"
mkdir -p "$FFT_WORK"
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export BAANDI_FFT_PROFILE=off BAANDI_PHASE_PROFILE=off
export BAANDI_STAGE_TIMERS=off BAANDI_DAVIDSON_STAGES=off
if [[ "$KIND" != time ]]; then
    export BAANDI_FFT_PROFILE=on BAANDI_PHASE_PROFILE=on
fi
if [[ "$KIND" == profile ]]; then export BAANDI_STAGE_TIMERS=on; fi
python3 -I scripts/fast_screen_dav.py --case "$CASE" --upf "$DFT_SI_UPF" \
    --threads 1 --seed 601 --amp 0.001 --emit-input "$FFT_WORK/input.toml"
if [[ "$MODE" == steps ]]; then
    python3 - "$FFT_WORK/input.toml" <<'PY'
import pathlib, re, sys
p = pathlib.Path(sys.argv[1])
s, n = re.subn(r'(?m)^max_iterations\s*=.*$', 'max_iterations = 10', p.read_text())
assert n == 1, 'expected one SCF iteration limit'
p.write_text(s)
PY
fi
sha256sum "$FFT_TARGET/release/dft-cli" "$FFT_TARGET/release/examples/davidson_bench" \
    "$DFT_SI_UPF" "$FFT_WORK/input.toml" > "$FFT_WORK/hashes.txt"
uptime > "$FFT_WORK/load-before.txt"
lscpu -e=CPU,CORE,SOCKET,NODE > "$FFT_WORK/topology.txt"
# Open MPI 5: report and check actual binding in every log; abort if not CPUs 0..7.
MPI=(mpirun -np 8 --map-by ppr:8:socket --bind-to core --report-bindings)
run_arm() {
    local tag=$1 setting=$2 rc=0
    export BAANDI_FFT_KEPT_FUSE="$setting"
    if [[ "$MODE" == solve ]]; then
        "${MPI[@]}" "$FFT_TARGET/release/examples/davidson_bench" \
            --input "$FFT_WORK/input.toml" --kinds warm --repeats 5 \
            > "$FFT_WORK/$tag.log" 2>&1 || rc=$?
    else
        /usr/bin/time -f '%e' -o "$FFT_WORK/$tag.wall" \
            "${MPI[@]}" "$FFT_TARGET/release/dft-cli" run "$FFT_WORK/input.toml" \
            --output-dir "$FFT_WORK/$tag" > "$FFT_WORK/$tag.log" 2>&1 || rc=$?
    fi
    printf '%s\n' "$rc" > "$FFT_WORK/$tag.exit"
    if [[ "$MODE" == steps && "$rc" == 3 ]]; then
        test -s "$FFT_WORK/$tag/results.json" # expected fixed-step non-convergence
    elif [[ "$rc" != 0 ]]; then
        return "$rc"
    fi
}
if [[ "$KIND" == counts ]]; then
    run_arm counter-A off
    run_arm counter-B on
else
    # Same-window A/A floor FIRST, then crossed order-alternated A/B.
    comparisons="AA AB"
    if [[ "$KIND" == profile ]]; then comparisons=AB; fi
    for comparison in $comparisons; do
        for pair in 0 1 2; do
            order='A B'
            if (( pair % 2 )); then order='B A'; fi
            for arm in $order; do
                setting=off
                if [[ "$comparison" == AB && "$arm" == B ]]; then setting=on; fi
                run_arm "$comparison-$pair-$arm" "$setting"
            done
        done
    done
fi
uptime > "$FFT_WORK/load-after.txt"
