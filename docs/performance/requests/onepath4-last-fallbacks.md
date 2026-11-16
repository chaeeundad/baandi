# dev08 request: ONEPATH-4 last strict fallbacks and engine default

Status: pre-registered (2026-11-16), not yet run. Ticket `tickets/ONEPATH-4-LAST-FALLBACKS.md`,
ADR 0043 amendment ONEPATH-4. BAANDI-only: no QE run and no QE ratio is requested, so the
parallel-model rule (`comparison-rules.md`) is satisfied trivially; every table is one model
at one thread count. Accuracy classes unchanged (strict 1e-10 Ha, standard 5e-9 Ha).

## What is being checked
The last strict fallbacks are gone: `scf.locking = "auto"` is an input error on the default
path (strict reference only), the observe/agent runtime modes and the non-`run` commands
(bands, dos, checkpoint, restart, explain) run the standard path with the fixed-point
verification, checkpoint v1 carries the verification option and attempts, and the engine
(`ScfOptions::default()`) and Python (`Options()`) defaults are the verified standard path.
Mac evidence (synthetic, ledger ONEPATH-4) meets every pre-registered bound; dev08 confirms
on real pseudopotentials and the MKL build that (a) the new paths deliver the standard class,
(b) exactness claims hold, (c) the reference switch is bit-identical to the former default.

## Build
Branch HEAD of ONEPATH-4 after merge with main. MKL build (`--features mkl`,
`MKLROOT=/opt/intel/oneapi/mkl/2026.1`, `MKL_NUM_THREADS=1`). Binary sha256, build command
and exit code recorded. Same binary for every arm. For C5 also the Python wheel
(`maturin build --release --offline --locked` from `python/baandi`, same tree).

## Cases
si8-gamma and si16-gamma (`benchmarks/cases/`), si-skew (generality set), each with the
PseudoDojo Si NC file of the benchmark suite (SHA-256 recorded; never committed). T=1
(taskset core 0) unless stated; node idle (load below 0.30), uptime before and after, arms
interleaved per case.

## Checks (pre-registered bounds, never widened)
- C1 locking auto. si8-gamma with `[scf] locking = "auto"`: (a) default environment: exit 1,
  no `results.json`, stderr names `scf.locking` and `BAANDI_REFERENCE_STRICT=on`; (b) with
  `BAANDI_REFERENCE_STRICT=on`: exit 0, `single_path.status = "reference_strict"`,
  `scf.locking = "auto"`, energy bit-identical to the former default build's locking run
  (main before ONEPATH-4, same input) or, if that build is not available, |dE| against the
  strict reference without locking <= 1e-10 Ha.
- C2 observe/agent. si8-gamma and si-skew with `[scf] cutoff_ladder = "off"`: classic,
  `--observe --events events.jsonl` and `--agent` (socket in the output directory) in the
  default environment. Pass: all exit 0, `single_path.status = "applied"`, `verified = true`;
  energy bits, density SHA-256 and verification attempt count identical across the three
  modes; the last `certificate_evaluated` event has `passed = true`. Against REF (the same
  modes under `BAANDI_REFERENCE_STRICT=on`): |dE| <= 5e-9 Ha. Report wall per mode.
- C3 checkpoint/restart exactness with verification. si8-gamma and si16-gamma, default
  environment: `checkpoint input exact.cp --iteration K` for K = 2 and K = (iterations of the
  uninterrupted run minus 1), `restart input exact.cp out.cp`, and
  `checkpoint input base.cp --kind result`. With `BAANDI_REUSE_PSI=off`: `dft-cli status` of
  `out.cp` and `base.cp` and the decoded result (energy bits, eigenvalues, history, verification
  attempts without `wall`) identical. With REUSE_PSI on (default): |dE| <= 1e-10 Ha and
  `verified = true` on both. Also one case where a verification fails before the checkpoint
  (si-skew at `standard_threshold_ha = 5e-11`; report whether any attempt failed; if none
  fails naturally, record that and skip).
- C4 bands/dos. si8 with the benchmark band path: `bands` and `dos` in the default environment
  and under `BAANDI_REFERENCE_STRICT=on`. Pass: exit 0; `bands.json.scf_path.verified = true`;
  max |eigenvalue difference| over all k and bands between the two <= 1e-6 Ha (first order in
  the density error at 5e-9 Ha; report the value) and the gap difference <= 1e-6 Ha.
- C5 Python session default. si8-gamma Gamma through `baandi.Session` with default `Options()`
  and with `reference_strict = True`, `mixing_beta = 0.3`: (a) default `converged = True`,
  energy bit-identical to the CLI default run with `cutoff_ladder = "off"`; (b) reference
  energy bit-identical to the CLI under `BAANDI_REFERENCE_STRICT=on` with
  `cutoff_ladder = "off"`; (c) |E_default - E_reference| <= 5e-9 Ha; (d) `stress()` raises
  ValueError at the default threshold and returns a tensor at `standard_threshold_ha = 5e-11`.
- C6 strict-reference identity. si8-gamma and si16-gamma plain `run` under
  `BAANDI_REFERENCE_STRICT=on` with this build and with main before ONEPATH-4 (same input,
  same machine): energy bits, density SHA-256, SCF iterations and certificate record
  identical (the reference path must be untouched). Same for `checkpoint --kind result`.

## Report per run
exit code, wall (s), `timings.verification_s`, SCF iterations, verification attempts and the
last attempt's `coulomb_error_ha`, `residual_bound_ha`; total energy (17 digits); density
SHA-256; `single_path` and `execution` blocks; for C3 the `dft-cli status` JSON of every
checkpoint.

## Not requested
No speed claim. Wall times are observations for the ledger; the verification cost is already
part of every default timed run (ADR 0043 V6).
