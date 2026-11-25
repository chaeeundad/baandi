# R8a stepper: dev08 identity and wall screen

Prepared for the master session / `baandi-dev08`; no remote command was run by
the implementer. Question: does splitting the reference loops into synchronous
map steps preserve the complete trajectory and delivered result, and add any
visible wall cost? This is a control-flow refactor, not a numerical optimization.

## Exactly four full-convergence runs

Per the R8A-STEPPER ticket, run **ONE run per arm**, with one thread per MPI
rank, si128 at MPI8. Do not widen the matrix or silently add repetitions.

| Job | Source | Verification | Runs |
| --- | --- | --- | ---: |
| final-r8a-main-si128-np8-off | main baseline | off | 1 |
| final-r8a-stepper-si128-np8-off | reviewed R8a | off | 1 |
| final-r8a-main-si128-np8-on | main baseline | fixed_point | 1 |
| final-r8a-stepper-si128-np8-on | reviewed R8a | fixed_point | 1 |

Use the unchanged si128 fixture already present on dev08, including its exact
geometry/displacement seed, cutoff, explicit FFT grid, band count, UPF bytes,
atomic start, threshold, beta/history, and all switches. Only `scf.verify`
differs between the two verification pairs. `scf.driver = "reference"`,
`scf.threads = 1`, and the usual full convergence budget apply. Record the actual
input paths and SHA-256; do not invent a remote fixture path or substitute a
nominal geometry. Pure/MKL/backend/compiler flags must match across each pair.

The local starting commit is
`81763a54f574b4a63b625303a91667626f89013b`. Record the exact reviewed candidate
commit or patch SHA-256 (the implementer leaves the work uncommitted). If main
has advanced, retain the starting-commit identity comparison separately from
any comparison to that newer main; unrelated main changes cannot be attributed
to R8a. Build each source in a separate checkout/target using the established
MPI+MKL setup, with the same toolchain and flags:

```sh
cargo build --release -p dft-cli --features mpi,mkl
```

Record command, exit code, source/patch hash, binary SHA-256, rustc version,
MPI/MKL versions, full switch manifest and input/UPF SHA-256. No installs or
server configuration changes are requested. Use the existing launcher and an
exclusive socket-0 window, physical cores 0–7, with actual bindings recorded:

```sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# input, output and binary are absolute paths recorded in the job manifest.
taskset -c 0-7 mpirun -np 8 --map-by core --bind-to core --report-bindings \
  "$binary" run "$input" --output-dir "$output"
```

Each job is a `final-` full-convergence job with a **7200 s** limit, queued by
the master. Record idle/load state and actual binding; serialize with builds,
gates and other measurement jobs. Use the existing external delivered-wall
wrapper (process startup through output and MPI finalization). Keep verification
and all failed verification/recovery work inside wall time. The final ordinary
candidate remains the delivered state, just as on main.

## Return evidence and interpretation

Return raw `results.json`, SCF records, stdout/stderr and process exit status for
every arm, including failures. Return one table row per arm with:

- Delivered external wall in seconds; verification wall separately in seconds.
- Ordinary iteration count, actual SCF maps (first-map redo and verification
  included, cached maps identified separately), verification attempts, and H
  columns including verification/recovery. Preserve per-map counts for diagnosis.
- Convergence/verification outcome, delivered energy in Ha and its exact IEEE-754
  bits, and density SHA-256 over ndarray iteration order with little-endian f64
  bytes. Hash the same raw density representation in both arms, not rounded text.
- Source, binary, input and UPF hashes, switches, rank count and binding.

Within each verify pair require identical energy bits, density hash, map counts,
H columns, retry decisions and verification outcomes. Any difference blocks R8a
acceptance; locate the first divergent record without changing tolerances.
Report candidate/main wall for that same verify setting only. One sample per
arm is a screen, not a median or a speedup claim. Compare any observed increase
to an already measured matching A/A floor if available; if that floor is absent
or the result is inconclusive, report the limitation to the master. Do not claim
"no regression beyond noise" from an unmeasured noise floor. This request makes
no QE comparison and no comparison between parallel models.

The local rank 1/3/8 stepwise tests use simulated collectives. They cannot replace
this real MPI8 identity and wall check. Checkpoint, runtime-policy and transport
validation belongs to R8b–R8d.
