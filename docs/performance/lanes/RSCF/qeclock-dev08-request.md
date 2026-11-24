# P3-QECLOCK: one diagnostic pair per cell, MPI8

Execute through the master session and `baandi-dev08` after review (current
performance brief; the old `dft5090` relay is retired). No remote execution,
QE installation, source reading, or server changes by this implementer. This
request is **diagnosis only**, under the Diagnostics level rule in
[comparison-rules.md](../../comparison-rules.md). No headline timing, median,
speed ratio, adoption decision, or additional A/B campaign is requested.

## Four executions total

| Job (7200 s limit) | Cell | Program | Ranks | Threads/rank |
|---|---|---|---:|---:|
| final-qeclock-baandi-si64-np8 | si64 | reviewed BAANDI, stages on | 8 | 1 |
| final-qeclock-qe-si64-np8 | si64 | QE 7.6, matched Davidson | 8 | 1 |
| final-qeclock-baandi-si128-np8 | si128 | reviewed BAANDI, stages on | 8 | 1 |
| final-qeclock-qe-si128-np8 | si128 | QE 7.6, matched Davidson | 8 | 1 |

Exactly ONE profiled BAANDI run and ONE QE run for each cell. No warmups or
repetitions. Keep failure output; report failures rather than silently replacing
samples. One case per scheduler job; use existing approved launch wrappers and
MPI/MKL installation. Compare MPI with MPI only.

## Identity before launch

Resolve the current matched inputs through the dev08 manifests (do not invent
paths or silently reuse partial mixing snapshots). Record complete inputs,
actual paths and SHA-256, cell/coordinates including displacement seed, cutoff
in Ha/Ry, explicit FFT dimensions, identical UPF bytes/hash, bands, Gamma,
fixed occupations, functional, no symmetry, atomic density/wavefunction start,
threshold (`5e-9 Ha` / `conv_thr=1e-8 Ry`), plain modified Broyden beta 0.7 and
history 8. Starting-wavefunction random perturbation/completion differences
must be reported: a common starting *policy* does not prove identical arrays.
If a condition cannot be matched, return the mismatch; do not label the pair
matched. Use the same eight physical cores/socket, rank order, core affinity,
MPI launcher configuration and one thread per rank for both programs.

BAANDI uses the reference driver, final cutoff only, verification on, and
`BAANDI_DAVIDSON_STAGES=on`. No other diagnostic switch is required. Use the
reviewed release `mpi,mkl` build, record source commit/dirty diff and executable
SHA-256 plus Rust/MPI/MKL versions. Keep the numerical switches identical to
current approved defaults; save the resolved run manifest.

QE uses the matched Davidson input (`diagonalization='david'`,
`diago_david_ndim=2`), Gamma real path and `verbosity='high'`; record full QE
version/executable hash and its printed parallel layout. Use one pool, one task
group, no band groups and serial reduced diagonalization (`-nk 1 -ntg 1 -nb 1
-ndiag 1`) so the printed routine is `regterg`, not `pregterg`. Record all flags.
Do not change or instrument QE source. The standard rank-zero clock table is
the input to this diagnosis.

Launch template (replace placeholders with reviewed manifest values, not new
cases; scheduler submits each invocation independently):

```sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
BAANDI_DAVIDSON_STAGES=on taskset -c "$cores" \
  mpirun -np 8 --map-by core --bind-to core \
  "$baandi_binary" run "$baandi_input" --output-dir "$baandi_output" \
  > "$baandi_stdout" 2> "$baandi_stderr"
# In the separate matching QE job, with exactly the same core allocation:
taskset -c "$cores" mpirun -np 8 --map-by core --bind-to core \
  "$qe_binary" -nk 1 -ntg 1 -nb 1 -ndiag 1 -in "$qe_input" \
  > "$qe_stdout" 2> "$qe_stderr"
```

Preserve actual commands, scheduler IDs, exit codes, start/end times, binding
output, complete stdout/stderr, results.json on every rank, density hashes and
verification pass/fail. No numerical reference is asserted from a clock alone.

## Parse and return

After confirming both complete runs and the identity checklist, for each cell:

```sh
python3 -I scripts/qe_clock_compare.py "$baandi_output/results.json" "$qe_stdout" > "$cell-clocks.md"
python3 -I scripts/qe_clock_compare.py "$baandi_output/results.json" "$qe_stdout" --json > "$cell-clocks.json"
```

Return the raw files and both tables. Values are **rank 0 wall**, not rank
maxima or sums; existing BAANDI max/min remain separate imbalance evidence.
CPU times appear in QE but BAANDI does not collect them. Per-call time is total
wall/calls, with missing or zero-call entries unavailable; preserve QE printed
time resolution. Check H/projection/back calls, startup inclusion and separate
verification calls against BAANDI records. An unqualified QE `calbec` is not
assumed equivalent to H-scoped `h_psi:calbec`.

Use the [exact/inexact map](../../../design/reference-scf.md#qe-matched-clock-table).
Unavailable `init_run`, `electrons`, `c_bands`, `fftw`, `rdiaghg`, Davidson
remainder, `sum_band`, `v_of_rho` and `mix_rho` entries are expected, not zeros:
BAANDI alternatives have distinct names. Do not subtract global `h_psi` or
`rdiaghg` from `regterg` (startup contributes to those global children). Do not
compare a standard reduced solve to QE's generalized solve. Resolving those
holes needs separately owned map/FFT/initialization boundary instrumentation;
this request does not authorize changing numerical work or QE code.

Mac validation is in the task report. No dev08 measurements have been performed
by P3-QECLOCK; no new README speed result is justified. The shared dashboard,
session registry and performance ledger remain with the orchestrator because
they are outside this ticket's owned files.
