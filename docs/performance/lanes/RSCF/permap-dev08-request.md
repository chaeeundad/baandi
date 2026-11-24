# P3-PERMAP: length/batch FFT dispatch and per-map reductions

Status: experimental, both new switches default off. No dev08 result is claimed.
Base checkout: `796258987d380bdacf57b82c1871b338c11d014f`.
Motivation (ticket-relayed, not reproduced here): P3s `e5ec8ae`, si128 MPI8,
verification off, 15 maps: 1D FFT kernels about 4.9 s; non-solve/non-density
map work about 0.66 s (energy_correction 0.25, mix 0.10,
density_transform 0.07, energy 0.20, potential 0.04). Those are aggregate
phase seconds across the maps, not seconds per individual map.

## Mechanisms and accounting

`BAANDI_FFT_KERNEL_BY_LENGTH=on` overrides `M2_FFT` only at the existing
`BackendLanes` kernel dispatcher. Plan layout, ownership, sharing, transposes,
normalization and batch construction are untouched. Pure builds keep rustfft.
For each actual (run ID, length, batch count, direction), an MKL build warms
both kernels, measures seven alternating-order samples, and selects the lower
median seconds per batch (ties use rustfft). A sample sums timed transforms
until 0.2 ms or 64 calls. Input restoration and allocation are outside each
timed call, and all calibration, scratch allocation and locking are confined
to construction. Input is deterministic and bounded; no SCF data, material,
cell label, cutoff or fitted length threshold enters the decision. Zero
batches are skipped. A process-wide registry reuses the first decision for
all later plans of that run with the same key, including worker-created plans.
Both kernel plans remain resident; execution only branches to the selected one.

This changes the premise of the rejected global M2_FFT/P3e experiment: short
and long lengths need not choose the same backend. It does not repeat lane
batching, FFT_DIRECT_PACK/SPARSE_CLEAR or FFT_FUSED_ROWS. Calibration is charged
to setup and delivered wall. `FFT_FUSED_ROWS` is conservatively bypassed for
calibrated dispatch so its independent rustfft row kernels cannot silently
ignore a choice. It is explicitly off in every requested arm.

`results.json.fft_kernel_choices` records process-local unique keys, chosen
kernel, both medians in seconds/batch, and calibration wall seconds. Off/pure
choices have null measurements and zero calibration time. Each MPI process
writes its own choices (`rankN/results.json` for nonroot ranks); collect all
eight files. No JSON, timing or registry lookup occurs in FFT execution.
These are distributed pair-lane choices, not full-box field-FFT choices.

Choices are stable within a run, but independent calibrations may choose
different kernels across processes or runs near a timing tie. Rounding-class
transform/SCF gates therefore apply; fresh-run bit reproducibility and
thread-count independence are not established for this opt-in experiment.
Do not promote a default based on the Mac pure fallback tests or one screen.

`BAANDI_REFSCF_PERMAP=on` makes two changes:

* Broyden retains unchanged Gram entries. Only entries touching the overwritten
  history slot and the RHS are recomputed, preserving operand orientation,
  compensated local sum and ascending-rank reduction order. At history p,
  forms fall from p(p+1)/2+p to 2p (44 to 16 at p=8). A newly supplied Metric
  invalidates the entire Gram cache, even with equal vector lengths. Ring
  eviction, startup, failure and prediction rules are unchanged.
* The descf integral uses exactly the existing LEAN_ENERGY cached mixed-density
  HXC fields. Its loop also accumulates the input/output charges and density
  difference norm. One four-scalar tree sum replaces separate one- and
  three-scalar tree sums. On a continuing successful map this saves three
  real-grid passes and one collective; each scalar's summation order is kept.
  Failed mixing and delivered maps keep the existing diagnostic path.

The mixed fields are still necessary for descf and the next Hamiltonian;
their two FFTs and XC evaluation are not eliminated. The output-density
forward transform and inverse transform of the mixer proposal are also
necessary and unchanged. `density_transform` includes the output forward
transform/residual subtraction; `mix` includes Broyden and the inverse
proposal transform; `energy_correction` includes mixed fields, descf and,
when fused, the diagnostic integral pass and combined reduction. Compare
both the phase split and their total: attribution alone is not a speedup.
The existing exact-density field cache remains responsible for next-map reuse.
No stopping rule, threshold, energy definition or verification policy changes.

## Exactly sixteen profiled SCF executions

Explicit ticket protocol: **ONE profiled run per arm per case**, no warmup
SCFs, repeated runs, A/A campaign, medians or new QE jobs. Kernel calibration
inside an opted-in run is part of that run, not an extra SCF execution.

| Case | Model | D | K | P | KP |
| --- | --- | --- | --- | --- | --- |
| si64 | MPI8, one thread/rank | off/off | on/off | off/on | on/on |
| si128 | MPI8, one thread/rank | off/off | on/off | off/on | on/on |
| elongated | MPI8, one thread/rank | off/off | on/off | off/on | on/on |
| rskew | MPI8, one thread/rank | off/off | on/off | off/on | on/on |

Columns give FFT_KERNEL_BY_LENGTH / REFSCF_PERMAP. All arms: verify off,
REFSCF_LEAN_ENERGY=on, M2_FFT=rust, FFT_FUSED_ROWS=off; other current adopted
switches remain identical. Use the original P3s/P3q physical input manifests,
including the established `/data/sig-dft/vc/` geometry/UPF provenance. Resolve
actual paths through the master; names are not sufficient to reconstruct a
fixture. Freeze cell, coordinates, cutoff, FFT grid, bands, XC, UPF bytes,
start and convergence threshold. Only switches and output destination differ.

Build one reviewed binary for all arms with the campaign's existing native
flags (record exact flags, compiler, MPI/MKL versions and exit codes):

```sh
CARGO_TARGET_DIR=target/p3-permap cargo build --release -p dft-cli --features mpi,mkl
cargo test --release -p dft-numerics --features mkl --lib kernel_by_length -- --nocapture
cargo test --release -p dft-solvers --features dft-numerics/mkl --lib kernel_by_length -- --nocapture
cargo test --release -p dft-engine --features dft-numerics/mkl --test refscf permap -- --nocapture
```

Run by `baandi-dev08` through the master, never from this implementer session.
Record source commit plus patch SHA-256, binary SHA-256 and each resolved
input/UPF/start hash. One case/arm per fast job, maximum 600 s, retain all
failed/timed-out output without automatically repeating it. Use socket-0
physical cores 0-7, preserve actual binding output and machine load; no
concurrent builds/benchmarks. Resolved command skeleton for D:

```sh
BAANDI_FFT_KERNEL_BY_LENGTH=off BAANDI_REFSCF_PERMAP=off \
BAANDI_REFSCF_LEAN_ENERGY=on BAANDI_M2_FFT=rust BAANDI_FFT_FUSED_ROWS=off \
BAANDI_DAVIDSON_STAGES=on BAANDI_PHASE_PROFILE=on BAANDI_FFT_PROFILE=on \
BAANDI_STAGE_TIMERS=on OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
mpirun -np 8 --map-by core --bind-to core --report-bindings \
  <reviewed-binary> run <frozen-reference-driver-input-with-verify-off>
```

Check the effective manifest (environment overrides input performance tables). Change only the two arm switches and destination for K/P/KP.

## Required output

* Each rank's full results JSON, stdout/stderr, completion status and exit code;
  resolved manifests and provenance above; density.bin hashes and numerical
  differences; delivered energy in Ha and difference from the existing strict
  reference (only if its provenance and physical identity are valid).
* FFT kernel choices and calibration seconds per key per rank; sum unique
  calibration costs within each process. Report selected MKL versus rustfft
  keys with length, batch and direction. Report H-scoped first/middle/last
  kernel seconds, calls and points, both directions, using FFT_PROFILE;
  retain packing, transpose, multiplication and nonlocal controls separately.
  Do not mistake calibration medians for production-kernel stage times.
* Map count, redo count, Davidson iterations/H calls/H columns, delivery status;
  `energy_correction`, `mix`, `density_transform`, `metric`, `energy`, `potential`
  seconds per map and aggregate, setup/calibration/output and external wall.
  Preserve same-rank rows and the chosen critical rank; never sum independent
  rank maxima. Present full wall including calibration beside steady map work.
* P vs D should retain bit-identical energies, densities, descf, dr2 and counts.
  K/KP are rounding class: strict 1e-10 Ha / standard 5e-9 Ha, unchanged physical
  definitions, report every trajectory/count difference. These verify-off
  screens cannot establish independent fixed-point verification or promote
  the adaptive selector's cross-run determinism. No new QE ratio is requested.

## Local validation record

Mac commands and real outcomes are recorded below after execution. No speed
conclusion, MKL runtime validation, new QE reference or adoption is inferred.
Ledger/README/dashboard files are outside this ticket's OWNED_FILES; master
reconciliation should link this request after review and dev08 evidence.

Mac pure-build results (2026-11-25, this worktree; logs in ignored `target/`):

| Command | Exit | Outcome |
| --- | --- | --- |
| `cargo check --workspace --tests` | 0 | Workspace and tests type-check |
| `cargo test -p dft-solvers` | 0 | 186 passed, 25 ignored |
| `cargo test -p dft-engine --test refscf` | 0 | 35 passed, 5 ignored |
| `cargo test -p dft-numerics` | 0 | 91 passed, 22 ignored |
| `cargo fmt --check` | 0 | Clean |
| `cargo clippy -p dft-solvers -p dft-numerics -p dft-engine --all-targets -- -D warnings` | 0 | Clean |
| `cargo check -p dft-solvers -p dft-numerics -p dft-engine --tests --features dft-numerics/mkl` | 0 | Type-check only; MKL unavailable |
| Same clippy command with `--features dft-numerics/mkl` | 0 | Initial exit 101: calibration sample range-loop lint; iterator rewrite fixed it |
| `CHECK_SCCACHE=0 BAANDI_TEST_TMP="$PWD/target/tmp" scripts/check.sh dft-numerics` | 0 | build/test/fmt/clippy PASS; sandbox blocked optional `ps` load inventory |

`kernel_by_length_preserves_transforms`: five grid/support families x ranks
1/3/8; unity-potential round trip and nonconstant-potential application,
relative tolerance 2e-14 (dimensionless), pure on/off bit equality and repeated
execution bit equality. `kernel_by_length_records_and_reuses_choices` checks
empty batches, per-run isolation and reused decisions; the MKL-gated assertion
checks selection against measured medians (not executed on Mac).

`permap_preserves_delivery_and_counts`: 12 synthetic two-rank A/B cases over
cubic/skewed/elongated cells, verify off, Coulomb verification, energy
verification and max-iteration delivery. Delivered energy difference **0 Ha**,
density/descf/dr2 bit equality and equal solver/verification counts.
`permap_gram_cache_ring_and_metric_invalidation`: 12 updates, history 4,
complex residuals, excluded high-G/zero modes, ring overwrite and a new
Metric after seven updates: **128 -> 82** computed/reduced scalar forms,
identical mixer outputs/coefficients/condition estimates. At full history,
14 -> 8 forms per update; metric replacement correctly restores all 14.
These are work-count/correctness observations, not wall-time measurements.
The registry count is 165 -> 167; switch tests were searched and no default
change is proposed. Ignored long-running tests were not executed.
