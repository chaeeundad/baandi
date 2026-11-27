#!/bin/sh
# Run a benchmark matrix sequentially (one relaxation at a time, identical threads).
#   run_matrix.sh PY EXE PSEUDO_ROOT THREADS WORKROOT OUTDIR MAX_EVALS CASE_JSON "SEEDS" "METHOD:START ..."
# Example (dev08, one core):
#   sh run_matrix.sh /path/venv/bin/python target/release/baandi /path/pseudo 1 /scratch/rp results/dev08-1core 200 \
#      cases/si8_AlP-pbe-k222.json "1 2 3" "baandi-bfgs:warm baandi-bfgs:cold baandi-fire:cold ase-bfgs:cold \
#      ase-lbfgs:cold ase-fire:cold ase-precon-lbfgs:cold a3-lbfgs:cold"
# LAUNCHER (environment, optional) is passed as --launcher, identical for every arm.
# SCF_THRESHOLD (environment, optional; default force-tied = 4.73e-11 Ha) is passed as
# --scf-threshold-ha to every arm; 'none' reproduces the Phase 1 setting (case text unchanged).
# Phase 2 native arm (ADR 0058): baandi-precon:warm (and :cold).
# Existing result files are skipped, so an interrupted matrix can be resumed.
set -u
PY=$1; EXE=$2; PSEUDO=$3; THREADS=$4; WORK=$5; OUT=$6; MAXE=$7; CASE=$8; SEEDS=$9; shift 9; ARMS=$1
HERE=$(cd "$(dirname "$0")" && pwd)
NAME=$(basename "$CASE" .json)
mkdir -p "$OUT"
for seed in $SEEDS; do
  for arm in $ARMS; do
    m=${arm%%:*}; st=${arm##*:}
    f="$OUT/$NAME-s$seed-$m-$st.json"
    [ -f "$f" ] && continue
    rm -rf "$WORK/$NAME-s$seed-$m-$st"
    "$PY" -I "$HERE/run_bench.py" --case "$CASE" --seed "$seed" --method "$m" --start "$st" --exe "$EXE" \
      --pseudo-root "$PSEUDO" --threads "$THREADS" --max-evals "$MAXE" \
      --launcher "${LAUNCHER:-}" --scf-threshold-ha "${SCF_THRESHOLD:-force-tied}" --workdir "$WORK/$NAME-s$seed-$m-$st" --out "$f" 2>>"$OUT/stderr.log"
  done
done
