# P3-VERWARM: one D/W screen on Si64 and Si128 MPI8

Status: default-off implementation; dev08 measurements pending. For the
master session and baandi-dev08, with no direct remote action by Codex.
Use the reviewed candidate commit containing P3-VERWARM for **both** arms;
record its full commit and binary SHA-256. The motivating ac6aee9c figures
are historical single runs, never the new D arm.

## Exactly four full-convergence executions

| Cell | Order | Parallel model | Verification | Executions |
| --- | --- | --- | --- | ---: |
| Archived work-order vc Si64, input SHA prefix af97b69b | D then W | MPI8, 1 thread/rank | fixed_point | 2 |
| Archived ac6aee9c Si128 | W then D | MPI8, 1 thread/rank | fixed_point | 2 |

D: `BAANDI_REFSCF_VERIFY_WARM=off`; W: `on`. Both explicitly set
`BAANDI_REFSCF_KEPT_HX=on`, `BAANDI_REFSCF_VERIFY_ENERGY=off`,
`BAANDI_REFSCF_VERIFY_REUSE=off`, `BAANDI_REFSCF_VERIFY_GATE=on`,
`BAANDI_REFSCF_VERIFY_GATE_FIRST=off`, `BAANDI_REFSCF_RMM=off` and
`BAANDI_REFSCF_LADDER=off`. Input `[performance]` overrides must agree with
the environment; inspect the effective manifest, not just the launcher.
Keep all other switches identical to the archived work-order configuration.

Recover exact input, full input and UPF SHA-256 and launcher from
`/data/sig-dft/wo/RESULT.md` and the archived Si128 run. A hash prefix is not
sufficient identity. Identical coordinates, cutoff, FFT grid, pseudopotential,
bands, start, threshold, beta/history, backend and diagnostics are required.
Use the ordinary pinned release `mpi,mkl` build (no new target-cpu flags),
OMP/MKL/engine threads 1, socket-0 physical cores 0–7 with archived binding.
Run sequentially without concurrent builds. No additional repetitions,
warmup SCFs, verify-off arm, QE execution or seed sweep is requested.
Si64 retains the 120 s timeout; Si128 uses the 600 s fast-job cap. Archive
failures/timeouts without altering parameters or rerunning.

Invocation template (fill only from the archived manifests):

```sh
BAANDI_REFSCF_VERIFY_WARM=off BAANDI_REFSCF_KEPT_HX=on \
  OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /usr/bin/time -p mpirun -np 8 <archived-binding> <candidate> run <D-input> --output-dir <D-output>
BAANDI_REFSCF_VERIFY_WARM=on BAANDI_REFSCF_KEPT_HX=on \
  OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /usr/bin/time -p mpirun -np 8 <archived-binding> <candidate> run <W-input> --output-dir <W-output>
<candidate> refscf-check --stop-to-verified <D-output>/results.json
<candidate> refscf-check --stop-to-verified <W-output>/results.json
```

The remaining shared settings above must be in each input or environment.
Store source/diff, exact build command, binary hash, compiler/MPI/backend
versions, all input/PP hashes, resolved manifest, binding/load, stdout/stderr,
results and exit codes. External wall includes setup, every failed check and
recovery, output serialization and MPI exit. Preserve rank-local phases;
never add componentwise rank maxima to infer a wall.

## Required return table

For each cell and arm report:

- verified external wall (s), total verification seconds (s), stop-to-verified
  interval (s), ordinary recovery seconds (s), and reader closure error (s);
- ordinary maps, first stop, each attempt's map/pass/warm flag and full-H
  columns (including retries), total verification H columns, gated skips;
- ordinary full-H columns and total kept-HX hits/replaced/delta-V columns;
- delivered energy (Ha), **signed E - E_strict (Ha)**, density hash,
  verification Coulomb errors (Ha), eigen tolerances (Ha) and final verdict;
- retained slab bytes/rank and peak RSS/rank with units, if available from the
  existing work-order instrumentation. No new full-profile run is authorized.

Si64's archived strict value is -253.3781677212637 Ha **only for its exact
work-order input**. Recover the matching strict result and provenance for
Si128; missing provenance must be reported as unavailable, never inferred or
filled by an extra run. Both arms must pass the unchanged fixed-point check
and standard 5e-9 Ha accuracy class against their matching strict oracle.
Compare D/W delivered state within rounding (local regression uses 1e-12 Ha),
all decisions, requested tolerances, attempts and work counts. Explain any
trajectory change instead of hiding it in totals.

For B=260 with two eligible attempts, expected verification starting-H saving
is 520 columns plus 260 for the first ordinary recovery map; confirm actual
counts. A warm flag records eligibility used, not a pass decision. New-density
and delta-V transforms remain charged. A single-run W/D ratio is a screen
observation, not a median, a noise-resolved improvement or default promotion.
This request has no QE comparison column and no cross-parallel-model ratio.

## Task-local attempt record and Mac evidence

Mechanism: preserve the candidate cache across a separately solved verification
map, then reuse it at the mixed recovery potential. This differs from
I25-VERIFY-COST V1 (reusing the failed map as an ordinary iterate), which remains
off, and from P3r RMM (changing the solver/trajectory), which remains off.
KEPTHX2's adopted fused density retention is reused without a second slab cache.
The independence and si128-sized work accounting are in the
[design section](../../../design/reference-scf.md#warm-verification-maps).

Commands/results are recorded below after execution. The shared attempts
ledger, optimization log, session/dashboard JSON and README generator are
outside OWNED_FILES; this local record is provided for master reconciliation.
No QE source was read, no UPF was committed, and no remote run was performed.

### Mac numerical results

`cargo test -p dft-engine --test refscf verify_warm_ -- --nocapture`
passed all five tests. Each row below was checked at **1, 3 and 8 simulated
ranks**, including identical attempt indices/verdicts, first estimator stop,
ordinary-map count, threshold and the unchanged next-ethr formula. Units of
energy are Ha; density differences are Bohr^-3.

| Synthetic fixture | Ordinary maps | Attempts (map: verdict) | Verification full-H columns D -> W | Max abs delivered energy difference |
| --- | ---: | --- | ---: | ---: |
| cubic, local | 5 | 5: pass | 6 -> 3 | 0 |
| skew, local | 5 | 5: pass | 6 -> 3 | 0 |
| elongated, local | 5 | 5: pass | 6 -> 3 | 0 |
| skew, nonlocal (s/s/p projectors) | 24 | 24: pass | 13 -> 10 | 0 |
| failed-check recovery, nonlocal | 28 | 23: fail, 28: pass | 30 -> 24 | 8.881784197001252e-16 |

These counts are the same at all three rank counts. Every W attempt records
an actual warm hit; every D attempt records cold. The first recovery after
failure is warm in W and cold in D. Maximum delivered density difference
across the matrix is 9.506e-16 Bohr^-3, solely in the recovery fixture; the
other four fixtures deliver bit-identical densities and energies. The tests
also check rank-count energy agreement to 1e-12 Ha and stop-to-verified
partition closure. No tolerance was relaxed and no wall-speed claim is made.

The initial workspace check found a new-test type inference error, then an
out-of-scope rank-rule test constructing `SolveRecord`. Both were corrected:
the explicit `Option<f64>` fixes inference; recording warm state on the map
rather than changing `SolveRecord` preserves its existing constructors. No
out-of-scope test was edited. The subsequent workspace check exits 0.

A live CLI round trip additionally ran the synthetic failed-check fixture
with D/W at both lightweight and full diagnostic levels. Per-attempt full-H
columns were **18 -> 15** at failed map 23 and **12 -> 9** at passing map 28.
The recorded warm flags were false/false in D and true/true in W, and every
`h_columns` equalled the existing `h_vectors`. All four reader outputs were
VERIFIED with **closure_error_s = 0**. This directly tests serialized warm
flags rather than relying only on engine records. Reproducer and raw output
are local scratch: `python3 target/tmp/verwarm-live.py`,
`target/tmp/verwarm-live.log`, `target/tmp/verwarm-live/*-reader.json`.
The reader unit test also confirms legacy files without the optional flag
produce the same verdict and clock fields as false/true flagged files.
The live fixture's first ordinary recovery (map 24) also reduces **18 -> 15
full-H columns** in both diagnostic modes. Total replaced starting columns
increase **78 -> 87**, exactly two 3-column verification hits plus one
3-column recovery hit. This measures the recovery saving separately from
the verification saving.

### Validation commands and environment

All commands use the pinned toolchain with `/opt/homebrew/opt/rustup/bin`
on PATH. Scratch variables `BAANDI_TEST_TMP`, `TMPDIR` and `TMPDIR_990` point
to this worktree's `target/tmp`; `CHECK_SCCACHE=0` avoids shared-cache writes.
The task's five requested command groups passed:

| Command | Exit | Result |
| --- | ---: | --- |
| `cargo check --workspace --tests` | 0 | Workspace and all test constructors type-check |
| `cargo test -p dft-engine --test refscf` | 0 | 53 passed, 5 existing ignored |
| `cargo test -p dft-cli stop_to_verified` | 0 | 10 passed across binaries, including warm-flag compatibility and live CLI fixtures |
| `cargo test -p dft-numerics switches` | 0 | 14 passed across binaries; 174 registered switches |
| `cargo fmt --check && cargo clippy -p dft-engine -p dft-cli -p dft-numerics --all-targets -- -D warnings` | 0 | Format and warnings-as-errors checks pass |
| `python3 target/tmp/verwarm-live.py` | 0 | Four D/W/diagnostics CLI round trips, actual failed attempts and reader closure |

The additional `scripts/check.sh fast dft-engine dft-cli dft-numerics` test
phase records **585 passed, 3 failed, 152 ignored across 58 binaries**. The
three failures are unchanged CLI cache-path fixtures rejecting checkout-local
scratch (`pseudopotential cache must be outside repositories`):
`pseudo::cache_resolution_and_no_network_side_effects`,
`pseudo::cached_fetch_checks_hash_and_provenance`, and
`pseudo::tests::fetch_boundary_verifies_before_publication_and_never_overwrites`.
No cache policy or out-of-scope source was changed to bypass that constraint.
The optional gate load probe cannot invoke `ps` in the sandbox; it continues
to the real gates. This is not a claim of a passing combined gate.
Raw command exits/logs are `target/tmp/verwarm-validation.json` and
`target/tmp/verwarm-*.log`; the preserved gate-test detail is
`target/tmp/verwarm-gate-bins.log`. Changes remain uncommitted for review.

Combined-gate completion: **exit 1**, build/fmt/clippy PASS, test FAIL only
for the three cache-path fixtures above (216 s gate wall, not a performance
measurement). The final workspace type-check also exits 0. All task acceptance
checks pass locally; combined-gate success requires rerunning the existing
cache tests with permitted scratch outside any checkout. No dev08 result or
default adoption is claimed.
