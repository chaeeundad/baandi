# G-OWNED-ROWS: default-off owned input bridge A/B

Status: DEV08_AB_REQUEST after Mac correctness gates in validation.md. Master
executes on dev08; this task has made no remote changes or timing claim.
Candidate: **M2_OWNED_INPUT**, off baseline, on candidate. Do not toggle the
already-adopted I20 lazy serial/owned projectors. No QE run requested.

## Hypothesis / non-overlapping mechanism gate

Remove the intermediate full ORC input matrix and redundant local extraction,
keeping full complex input, full validation and the solver's local start storage.
At a successful eligible entry with G global rows, L_r local rows and B bands:
A bridge allocation requests = 8(G+L_r)B, B = 8L_r B bytes. Saving = 8GB.
At fine Si64 25,207,776 B/rank/entry; fine Si128 99,530,080 B/rank/entry.
No >=0.15 s forecast is defensible: setup/entry spans include other work and
most warm starts are already lazy. Memory qualifies the experiment.
If every start is already local, calls=0 makes this A/B **inapplicable**; do not
force materialization in production arms just to manufacture a win.

Mandatory diagnostics for **every rank**, JSON `run.counters` keys:

- `owned_lane.input.calls`: same positive number in A/B.
- `owned_lane.input.direct_calls`: A=0, B=calls.
- `owned_lane.input.baseline_bytes`: same in A/B; sum 8(G_i+L_ri)B_i.
- `owned_lane.input.allocated_bytes`: A=baseline; B=sum 8L_ri B_i.
- `owned_lane.input.saved_bytes`: A=0; B=sum 8G_i B_i > 0.
- B allocated+saved=baseline exactly. Compare per-rank then sum; do not compare
  rank-max allocated with independently maximized baseline/saved.
- Existing init.problem, state_grid, hamiltonian, band_density, owned_rows_build
  construction calls must remain equal. Same SCF/Davidson/H/verification work,
  energy bits, density bytes, eigenvalues and full normalized scientific history.

Counters are allocation **requests for these matrices**, not a malloc hook,
DRAM traffic, total allocated bytes or RSS. Record separate per-rank Linux peak
RSS. A reduction of virtual calloc payload need not equal peak RSS reduction.

## Build and prepare (outside timed window)

Use the master-managed MPI/MKL toolchain, record build exits and commit/diff:

```bash
MKLROOT=/opt/intel/oneapi/mkl/2026.1 CARGO_TARGET_DIR=target/owned-mpi \
  cargo build --release -p dft-cli --features mpi,mkl
export OWNED_BIN="$PWD/target/owned-mpi/release/dft-cli"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# UPF: master's existing absolute Si.pz-vbc.UPF path; never downloaded here.
sha256sum "$OWNED_BIN" "$UPF"
git rev-parse HEAD
lscpu -e=CPU,CORE,SOCKET,ONLINE
```

Run **one case per job** in order si32, si64, si16, then si128. Each of these
is an exact separate invocation of the following recipe, with CASE set once:

```bash
CASE=si32  # next independent jobs: si64, si16, si128
export OWNED_JOB="/data/sig-dft/lanes/OWNED/$CASE-np8"
mkdir -p "$OWNED_JOB"
python3 -I docs/performance/lanes/OWNED/prepare_input.py \
  --case "$CASE" --upf "$UPF" --output "$OWNED_JOB/input.toml"
"$OWNED_BIN" validate "$OWNED_JOB/input.toml" --json > "$OWNED_JOB/validate.json"
"$OWNED_BIN" grid "$OWNED_JOB/input.toml" > "$OWNED_JOB/grid.txt"
sha256sum "$OWNED_BIN" "$UPF" "$OWNED_JOB/input.toml" > "$OWNED_JOB/identities.txt"
```

Si128 uses **benchmarks/cases/si128-gamma.toml**, 512 electrons, 260 bands,
47,851 fine PW, 72x72x144 FFT at 30 Ry. Confirm with actual preflight, input and
output, not only this expected table. Record coarse stage separately. One
physical core per rank, eight ranks, CPUs 0-7 socket 0, same binding/input/UPF/
backend/caches/output filesystem and verification in both arms. Clear unrelated
BAANDI experimental environment variables before the window and preserve the
full effective switch manifest. Warm binary-hash cache equally for both arms.

## External wall: A/A first, then crossed pairs

Exclusive window, load below 0.30, no builds or other jobs. Normal job cap 600 s;
Si128 full verification may use its fixture's 7200 s final-job cap after master
memory qualification. Never silently drop a timeout/OOM. Launch and exit,
verification and output stay inside wall. These functions use bash and GNU time
on dev08; Mac does not execute them.

```bash
run_owned() {
  local name="$1" setting="$2"
  env BAANDI_M2_OWNED_INPUT="$setting" BAANDI_INIT_PROFILE=off \
      BAANDI_PHASE_PROFILE=off BAANDI_STAGE_TIMERS=off \
    /usr/bin/time -f '%e' -o "$OWNED_JOB/$name.wall" \
    mpirun -np 8 --cpu-set 0-7 --map-by core --bind-to core --report-bindings \
      "$OWNED_BIN" run "$OWNED_JOB/input.toml" --output-dir "$OWNED_JOB/$name" \
      > "$OWNED_JOB/$name.stdout" 2> "$OWNED_JOB/$name.stderr"
}
run_owned warm-a off
run_owned warm-b on
# First pair discarded by pre-registration, six retained A/A pairs.
for i in 0 1 2 3 4 5 6; do
  run_owned "aa-$i-a" off || exit 1
  run_owned "aa-$i-b" off || exit 1
done
```

Master checks A/A **before** starting B: mean drift must be unresolved
abs(mean)/SE <3 (zero SE is not statistical evidence); otherwise stop and
repeat a whole quiet window. Floor = max absolute paired A/A difference on
retained pairs. Then 12 retained crossed pairs after discarding pair 0:

```bash
for i in $(seq 0 12); do
  if (( i % 2 == 0 )); then
    run_owned "ab-$i-a" off || exit 1
    run_owned "ab-$i-b" on || exit 1
  else
    run_owned "ab-$i-b" on || exit 1
    run_owned "ab-$i-a" off || exit 1
  fi
done
```

Savings d_i=A_i-B_i. Report all raw times, means/medians, sample SD,
SE=SD/sqrt(12), abs(mean)/SE, A/A floor and elapsed-clock resolution. Require
positive mean, lower B median, abs(mean)/SE>=3 and effect>=3*floor. If GNU time's
0.01 s wall resolution dominates, use the master's high-resolution external
launcher for **both entire A/A and A/B series**, never combine clocks. Unresolved
stays off. Any resolved general gain qualifies; no minimum percentage. Memory
improvement alone is documented separately and is not speed adoption evidence.

## Separate count/RSS diagnostic (never acceptance wall)

Three crossed pairs per case. Turn INIT_PROFILE on; PHASE_PROFILE and
STAGE_TIMERS on for attribution, or clocks off for count-only corroboration.
Use GNU time **inside** mpirun for each rank, not launcher RSS:

```bash
export OWNED_LABEL=diag-a BAANDI_M2_OWNED_INPUT=off
# Repeat with diag-b/on; alternate order for subsequent diagnostic pairs.
mkdir -p "$OWNED_JOB/$OWNED_LABEL"
env BAANDI_INIT_PROFILE=on BAANDI_PHASE_PROFILE=on BAANDI_STAGE_TIMERS=on \
  mpirun -np 8 --cpu-set 0-7 --map-by core --bind-to core --report-bindings \
  bash -c 'exec /usr/bin/time -v -o "$OWNED_JOB/$OWNED_LABEL/rss-rank${OMPI_COMM_WORLD_RANK}.txt" \
    "$OWNED_BIN" run "$OWNED_JOB/input.toml" --output-dir "$OWNED_JOB/$OWNED_LABEL"' \
  > "$OWNED_JOB/$OWNED_LABEL.stdout" 2> "$OWNED_JOB/$OWNED_LABEL.stderr"
```

Return eight results.json and density.bin files per run (root and rank1..7),
RSS in KiB per rank and maximum, minor/major faults, profiles, call/byte counters,
source/build/binary/input/UPF hashes and binding. Same-rank A/B normalized
histories must match; no cross-rank-count bit identity assertion (MPI reductions
have a different fixed association). Retain simulation rank-one negative control.
Before default adoption extend to skew/elongated/slab and synthetic/non-Si PPs,
strict controls, M1 fallback, and snapshot/owned-handoff switch combinations.

Withdraw on any correctness/fallback/counter mismatch, new allocation outside
the stated sites, or resolved loss. If counters move but RSS/time do not, record
that result and keep OFF; do not credit historical I20 savings. Remaining idle
BandDensity storage is a separate owner ticket, not a reason to expand this A/B.
