# MIXER: cell-filter A/B on dev08 (master relay only)

Candidate C: `BAANDI_MIX_CELL_KERKER=on`; A: off. Default remains off.
This is a trajectory-changing Stage 1/2 request, **to convergence**, with the
unchanged standard stop and fixed-point verification included. Mac counts
motivate the test; Mac durations provide no speed evidence. No QE runs or
cross-parallel-model ratios are requested. The first version intentionally
uses the existing replicated Kerker mixer even when M2 is active; measure its
density gathers and FFT cost explicitly. No claimed distributed-filter speedup.

## Inputs and build identity

Use the exact archived si128 np8 input from `/data/sig-dft/lanes/convergence-history.md`
and its source run, then si216 np8, si64 np8 and si32 np8 as controls. Do not
substitute the Mac displaced inputs for those geometries. One input/cell/seed
per job, one rank per physical core 0-7, one thread per rank. Keep pseudopotential,
cutoff, FFT grid, bands (si128: 260), start, ladder, threshold and binding identical.
Retain verification ON, `STOP_PREDICTOR=off`, standard threshold 5e-9 Ha.

The master records baseline commit `6a060639934b60859c1a70f6bd8a918b9aac315f`,
its current main and the reviewed candidate commit/diff, build commands and
exit codes, Rust/MPI/MKL versions and binary/input/UPF hashes. Source changes
are deliberately uncommitted by the implementer. Build the reviewed candidate
once with the approved MPI/MKL configuration, e.g.:

```sh
MKLROOT=/opt/intel/oneapi/mkl/2026.1 cargo build --release -p dft-cli --features mpi,mkl
sha256sum target/release/dft-cli
```

Before timing, compare candidate-off with current-main on each archived input:
energy/density bits and scientific histories/work counters must match. Obtain
or reuse a **certified tight strict reference for the same physical input**;
record its binary/input/UPF hashes and energy. Tight input uses energy 1e-12 Ha,
density 1e-10 Bohr^-3/2, eigen residual 1e-11 Ha with
`BAANDI_REFERENCE_STRICT=on`. Never use the last energy of a failed reference.
Check geometry/cutoff/grid/PP/bands identity manually against both inputs; the
runner checks certification, not cross-file physical identity.

## A/A first, then crossed pairs

Fill these task variables with the archived paths; no regeneration of the
geometry, UPF or production settings is authorized by this script:

```sh
MIXER_BIN=/absolute/path/to/reviewed/target/release/dft-cli
MIXER_INPUT=/absolute/path/to/archived/si128.toml
MIXER_UPF=/absolute/path/to/the/same/Si.pz-vbc.UPF
MIXER_STRICT=/absolute/path/to/certified/si128-strict/results.json
MIXER_RUNS=/data/sig-dft/lanes/MIXER
```

Each command below is a separate exclusive queue job. The runner refuses an
existing output directory, records exact commands/switches, and bounds the
**whole job**, not just each child. It runs MPI8 with `taskset -c 0-7`,
`--bind-to core --map-by core`, OMP/MKL/OpenBLAS threads 1. Confirm this matches
the archived affinity and record the rank-to-CPU map. No concurrent benchmark,
build, sample/perf or clock-instrumented alternative in the timing window.

```sh
python3 -I docs/performance/lanes/MIXER/dev08_pairs.py --binary "$MIXER_BIN" \
  --input "$MIXER_INPUT" --upf "$MIXER_UPF" --strict-result "$MIXER_STRICT" \
  --sequence AA --job-name mixer-si128-aa-01 --budget 600 --out "$MIXER_RUNS/si128-aa-01"
python3 -I docs/performance/lanes/MIXER/dev08_pairs.py --binary "$MIXER_BIN" \
  --input "$MIXER_INPUT" --upf "$MIXER_UPF" --strict-result "$MIXER_STRICT" \
  --sequence AB --job-name mixer-si128-ab-01 --budget 600 --out "$MIXER_RUNS/si128-ab-01"
python3 -I docs/performance/lanes/MIXER/dev08_pairs.py --binary "$MIXER_BIN" \
  --input "$MIXER_INPUT" --upf "$MIXER_UPF" --strict-result "$MIXER_STRICT" \
  --sequence BA --job-name mixer-si128-ba-01 --budget 600 --out "$MIXER_RUNS/si128-ba-01"
```

Repeat A/A in three separate windows before interpreting A/B; complete at least
three AB and three BA jobs per cell. Alternate AB/BA windows. Apply the exact
commands separately to si216, si64 and si32 by changing the archived paths,
job names and output directories. Never combine cells into a job. If the
measured A/A durations cannot fit two runs plus output in 600 s, declare that
before launching A/B and use a `final-` job **with a 7200 s total cap**, e.g.:

```sh
python3 -I docs/performance/lanes/MIXER/dev08_pairs.py --binary "$MIXER_BIN" \
  --input "$MIXER_INPUT" --upf "$MIXER_UPF" --strict-result "$MIXER_STRICT" \
  --sequence AB --job-name final-mixer-si216-ab-01 --budget 7200 \
  --out "$MIXER_RUNS/final-si216-ab-01"
```

No capped/failed/uncertified pair counts as a completed timing pair; retain all
its logs. No retries chosen from favorable timings. Extend to the registered
geometry seeds after the first screen and before default consideration. A real
MPI8 job is mandatory: simulated rank threads on the Mac are correctness only.

## Pre-registered counters and verdict

For every coarse and fine map return full `convergence.history` and ladder
`coarse_history`, `standard_error_ha`, normalized `density_residual`,
`delta_energy_ha`, `requested_eigen_tolerance_ha`, `mixing_fraction`,
`mixer_reset`, `mixing_resets`, `mixing_history_used`, `mixing_raw_recovery`, and
`pulay.{predicted_ratio,coulomb_predicted_ratio,coulomb_refit,coefficients,kept}`.
A **divergent stretch** is a maximal consecutive increase in standard error:
record inclusive map endpoints, its length in intervals and final/start ratio;
also report maximum error/best-previous error and the L2 ratio separately.

C's existing Kerker wrapper returns null Pulay diagnostics. Report that as
**unavailable under C**, never zero or no divergence. A's full predicted trace
is mandatory; native filtered-trace forwarding requires the wrapper ownership
extension noted in the design. Resets in C remain feasibility/raw-recovery
resets; C is not the proposed session-level recovery controller.

Primary mechanism: fewer total **fine maps to passed verification**, especially
a shorter early divergent stretch on si128. No preregistered numerical promise
is made from the different Mac geometry. The si32/si64/si216 controls must have
no new convergence/certificate failures; no-regression maps/H and full wall
are reported separately. Record H vectors by fine/coarse/certificate/start,
Davidson iterations, verification attempts (failed and passed), charge error,
final energy, maximum RSS, mixing time, owned-density declines/gathers and all
effective switches. Certificate H is already present in the fine history;
do not add verification aliases twice. Separate per-map cost and map count.

For paired full-process walls use `d_i = log(t_B/t_A)`, aligned by AB or BA,
report mean, sample SD, `SE = SD/sqrt(n)` and exp(mean). Require a favorable
mean with `abs(mean)/SE >= 3`, beyond the observed A/A floor, and consistent
AB/BA signs. No minimum percentage gain. If unresolved, report unresolved and
preserve the OFF default. Also report medians and raw pairs; do not select only
fast seeds. Standard delivered energy must be within 5e-9 Ha of the strict
reference, verification must pass; any strict candidate check uses 1e-10 Ha.

Withdraw/revise if saved maps are outweighed by replicated-filter/gather cost,
controls regress beyond noise, convergence/certificate/energy fails, charge
conservation fails, or arithmetic depends on thread count. Rank-count arithmetic
changes must be separately reported, not represented as bit-identical. A result
on si128 alone is not default approval. G/E designs are not executable arms in
this patch; no nonexistent switch may be silently timed as baseline.
