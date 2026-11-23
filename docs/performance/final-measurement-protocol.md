# Final performance measurement protocol

## Current state — 2026-11-20..24 (issue #26)

The new reference SCF path (issue #25) is the active work; the threads lane is
dropped. Defaults are unchanged by this reconciliation. dev08 is operated by
`baandi-dev08`; the `dft5090` relay is retired. See the
[current operations and protocol map](fast-loop-protocol.md#current-state--2026-11-2024-issues-25-and-26).
Use fixed-H replay for kernels, common-input replay for density-update
discrepancies, and full convergence for tail, verification or trajectory changes.
Timing uses [the same median for both engines](comparison-rules.md#statistics-and-delivered-wall);
fast jobs are limited to 600 s, authorized `final-` jobs to 7200 s.

## Current campaign instructions

Use the [protocol map](fast-loop-protocol.md#protocol-map-live) to choose the
experiment: fixed-H replay for kernels, common-input replay for density-update
discrepancies, full convergence for tail/verification/trajectory changes.
Only the full-convergence campaign supplies final delivered-result tables;
replay or ten-step numbers must stay labelled diagnostics.

- Pre-register the source commit, switches, build features and binary hashes,
  cases/geometries, matched inputs and accuracy references. Keep defaults
  unchanged until the multi-system promotion gates pass. Issue #25 being active
  does not mean its reference path is promoted.
- Run through `baandi-dev08` under the
  [current operations and job classes](fast-loop-protocol.md#job-classes-and-promotion-live):
  socket 0, np8 cores 0-7 or np16 cores 0-15; one core as appropriate.
  Fast screens stay within 600 s. Full runs or longer diagnostics use authorized
  `final-` jobs up to 7200 s, one case per job. The socket-1 gate is serialized
  through `final-a0-gate-<hash>` from `/scr/mino/gate`.
- Compare only matching parallel models and identical physical/numerical
  conditions under [comparison rules](comparison-rules.md). No threads-lane
  campaign or BAANDI-threads/QE-MPI column is authorized by the old P10 text.
- Use [one statistic and timing boundary](comparison-rules.md#statistics-and-delivered-wall):
  at least three samples, **median for both engines**, spread and raw records.
  Include all verification, failed attempts and recovery in delivered wall;
  show verification separately. Retain failures/timeouts and noisy batches.
- Require full convergence, independent fixed-point verification, strict
  **1e-10 Ha** / standard **5e-9 Ha** delivered-accuracy gates and applicable
  determinism checks. Report energy differences to independent tight references
  in Ha, density comparisons, SCF/Davidson/H work and verification attempts.
  The 1e-12 Ha unique-solution cross-path check is a separate regression check.
- Default promotion requires the full generality set, not Si alone: varied cell
  shapes, sizes, k meshes, band counts and projector structures. Report measured
  limits for missing/failed cases; never relax tolerances to fill a table.
- Preserve commands, versions, input/UPF SHA-256, source/binary hashes, queue
  IDs, actual binding, load, raw results and relayed provenance. No QE source.
  Generate README evidence/figures through the established generator only after
  verified evidence is available; this documentation change adds no results.

## Historical P10 campaign

The original request below is preserved for provenance. Its mixed parallel
models, unequal statistics, old prerequisites and timing/accuracy definitions
are superseded by the current instructions above. Do not rerun it verbatim.

<details>
<summary>Historical — original post-P10 final-measurement instructions; superseded, not live instructions</summary>

# Final performance measurement protocol (after the P10 tickets)

Purpose: the README performance section is rewritten from one coherent measurement campaign on the reference node with Intel MKL available, once the P10 work (MKL backend, fused FFT path, and the solver/mixing modules that were merged) is in `main`. Until then the README keeps the pure-Rust numbers and says so.

## Preconditions
- P10B (MKL backend) and P10C (fused FFT call structure) merged, plus whatever of P10A/P10D was merged; the commit under test is recorded.
- Binaries built on dev08 (2 x Xeon Gold 6230, AVX2 and AVX-512 present, Rocky Linux 9.7) with `MKLROOT=/opt/intel/oneapi/mkl/2026.1`: (a) default pure-Rust build, (b) `--features mkl` build. Binary sha256 recorded.
- Node idle (load below 0.3 before each batch), uptime before and after, nothing else running, cores pinned (`taskset -c 0-(T-1)` on one socket), MKL sequential layer, `MKL_NUM_THREADS=1`.

## Cases and layouts
- Cases: Gamma, 2x2x2 shifted (4 k), 2x2x2 (8 k), 4x4x4 shifted (32 k), the 8-atom silicon cell at Gamma (default settings, after P7PC), and the cases of the general suite (ticket P10G) flagged as "default set", so that the numbers are not silicon-specific (AGENTS.md generality rule).
- Layouts: 1, 4, 8, 16 threads for BAANDI; QE 7.6 with the same number of MPI processes (one thread each, `--map-by core --bind-to core`).
- Accuracy classes: QE at `conv_thr` 1e-8 (energy already equal to the 1e-12 run to the printed digits) and 1e-12; BAANDI at its defaults and with its energy difference to the previous state reported.
- Repeats: 3 runs per cell, median (BAANDI) and mean (QE), spread reported; any batch with a spread above 5 percent is repeated on an idle node or discarded.

## Reported quantities
- Wall time per case and layout for: BAANDI pure Rust, BAANDI with MKL, QE 7.6 (both thresholds); speedup MKL over pure Rust; ratio BAANDI over QE for the same cores.
- Numerical gates for every cell: total energy difference between pure Rust and MKL builds (gate 1e-12 Ha), QE agreement scripts unchanged and passing, density hashes identical across thread counts within each backend.
- Kernel level (`crates/dft-numerics/examples/lib_bench.rs`, `tools/lib_bench`) for the same build, to relate kernel and SCF speedups.

## Outputs
- `tests/reference/perf/dev08_final.json` (relayed attestation, binary hashes, uptimes), `docs/validation.md` R7 update, README performance block and figures regenerated by `scripts/make_readme_figures.py` (extended with a backend dimension), `docs/feature-status.json` speed line updated, dashboard test entries.
- The README states the build used (feature set, MKL version), the node, and that QE also used MKL.

</details>
