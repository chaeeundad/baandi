# P3-KEPTHX2: one D/H screen on Si64 and Si128, MPI8

For the master session and `baandi-dev08`. No remote run was performed by
this task. Candidate base is `b69571dbf512a323c17b03732f168c33e2683ae7`
plus the reviewed P3-KEPTHX2 diff. `BAANDI_REFSCF_KEPT_HX` remains default off.

## Premise and change

The ledger's **Work-order screen** at `bb416c73` (raw evidence on dev08:
`/data/sig-dft/wo/RESULT.md`) reports Si64 MPI8 D/H wall 6.41/6.60 s,
density 0.371/0.850 s, fine full-H columns 6094/4642, 11 cache hits,
1452 replaced starting columns, and peak RSS 206/259 MB per rank. Solve saved
0.15 s but retaining slabs erased the gain. The earlier adopted M2_REUSE
mechanism establishes the potential-change identity, not the performance of
this reference driver. This follow-up removes the per-map slab allocation
and separate copy scan; it does not change the eigensolver or default.

[Density accounting and layout](../../../design/reference-scf.md#p3-kepthx2-density-pass-accounting-and-storage-lifetime)
identify equal FFT work before and after: ceil(B/2) inverse pairs, no forward
pairs per density call. The retained branch now stores directly while
accumulating density. The stage keeps its slab allocation across maps,
including cache invalidations; snapshot X/HX/potential costs remain.

## Exactly four runs

Use the same candidate binary for D (`BAANDI_REFSCF_KEPT_HX=off`) and H (`on`).
One execution per arm per cell; no repeat, warmup SCF, seed sweep or extra
verification-off arm. Do not use the historical wall as the new D measurement.

| Cell | Order | Parallel model | Verification | Runs |
| --- | --- | --- | --- | ---: |
| Archived work-order vc Si64 (input SHA starts `af97b69b`) | D then H | MPI8, one thread/rank | on | 2 |
| Archived work-order Si128 | H then D | MPI8, one thread/rank | on | 2 |

Recover the **exact** commands, full input hashes and manifests from the
archived work-order run; a hash prefix is not an input identity. Same cell,
coordinates, PP bytes, cutoff, FFT grid, bands, start, beta/history, threshold,
backend, diagnostics and verification judgment in both arms. Keep ladder and
RMM off, as in D/H. Keep the work-order socket-0 physical cores 0–7 and rank
binding. OMP/MKL threads and engine threads are all one. Run sequentially
with no concurrent build. Si64 retains its 120 s timeout; Si128 uses the
600 s fast-job cap. Report failures/timeouts without rerunning or changing
accuracy. No QE execution or mixed-model comparison is requested.

Invocation pattern (substitute the archived launcher/binding and input paths;
retain any additional archived settings rather than guessing them):

```sh
BAANDI_REFSCF_KEPT_HX=off OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /usr/bin/time -p mpirun -np 8 <archived-binding> <candidate> run <D-input> --output-dir <D-output>
BAANDI_REFSCF_KEPT_HX=on OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /usr/bin/time -p mpirun -np 8 <archived-binding> <candidate> run <H-input> --output-dir <H-output>
```

Use separate output directories. Archive source commit and diff, binary SHA,
build command/features, toolchain/MPI/backend versions, full input and PP
SHA-256, environment and resolved switches, binding, load, stdout/stderr and
results. Charge setup, failed verification/recovery, output and MPI exit to
external wall. Keep the existing work-order diagnostics identical in D/H;
additional heavy profiling is outside these four runs and needs a separate
request. Rank-wise phase values must be preserved, not only sums of maxima.

## Return table and decision

For each cell report D/H external wall (s), density time (s), solve time (s),
verification/recovery time (s), peak RSS per rank (MB with unit definition),
retained slab payload per rank (bytes), fine map/solve counts, full-H columns,
cache hits, replaced columns, delta-V columns, and density inverse/forward
pair counts if available from the existing instrumentation. Compare matched
maps; full-run FFT totals also contain the intentionally saved warm-start
inverse transforms. Include complete energy (Ha), density hash, verification
verdict and attempts, exit status, and exact run provenance.

The screen's target is **H density time <= D density time**, with the expected
full-H saving and no wall regression. Report the raw H/D wall ratio as one
screen observation, not a median or a performance adoption. If the density
regression persists, separate first allocation from subsequent maps using
existing per-map records; do not attribute it to extra FFTs without counts.

Both arms must converge and pass unchanged independent verification. Use the
matching archived strict oracle (Si64 -253.3781677212637 Ha, only for the exact
archived input); report |E-E_strict| <= 5e-9 Ha. Missing Si128 oracle provenance
is a qualification limit, not permission to fabricate one or add runs.
D/H need not be bit-identical because kept HX is already rounding class;
**new H versus old H** must be bit-identical at fixed configuration. Compare
archived old-H hashes only when build/backend/input conditions match; otherwise
label that dev08 comparison unavailable. Mac pre/post bit checks are separate
internal evidence. Defaults remain off until the master reviews measured
benefit and the wider correctness/generality gates.

## Local evidence

Mac commands and results are recorded below after execution. This request
contains the task-local attempt record because the shared attempts ledger,
optimization log, dashboard/session JSON and README generator are outside
this ticket's OWNED_FILES. The master session can reconcile those artifacts.

Pre/post check on this Mac: the new delivery test was compiled once against
both production files from base HEAD and once against the candidate files,
with the same new test source and toolchain. Both runs exited 0 and produced
fingerprint `eda186a682b49798`, energy **-0.839575255872316961 Ha**, 28 maps,
27 hits and 81 replaced starting columns. The fingerprint covers every
returned density/eigenvalue/occupation/orbital bit, total energy, charge and
the named counts (not timings); each version also repeats the calculation.
This is a synthetic skew-cell/nonlocal two-rank internal oracle, not Si64,
MKL, real MPI, a portable fixed hash, or QE evidence. The pre/post harness and
raw logs remain in `target/tmp/kepthx2-prepost.py` and
`target/tmp/kepthx2-{before,after}.log` in this worktree.

The counter/bit test covers 1/2/5 simulated ranks, 3/4 bands and 0/1/all complete
retained pairs, including a zero-occupation band. Every density call has
**2 inverse pair transforms and 0 forward pair transforms** on each rank,
plain and keeping; every compared density/slab/delta-HX bit agrees. These
counts are isolated from concurrent tests in a child process. Existing
`fast_order_kept_hx_delta_and_invalidation` still reports relative product
error **4.39708409238569948e-17**, fresh/hit/invalidated starting full-H
columns **3/0/3**, and retained slab payload **4096 bytes/rank**.

| Local command | Exit | Result |
| --- | ---: | --- |
| `cargo test -p dft-engine --test refscf fast_order_ -- --nocapture` | 0 | 4 passed |
| `cargo check --workspace --tests` | 0 | All workspace tests type-check |
| `cargo test -p dft-engine --test refscf` | 0 | 48 passed, 5 ignored |
| `cargo test -p dft-physics` | 0 | 97 passed, 18 ignored across binaries |
| `cargo test -p dft-engine --lib refscf::kept_hx::tests` | 0 | 1 passed; allocation pointer/capacity retained across 3 subsequent maps |
| `cargo fmt --check && cargo clippy -p dft-engine -p dft-physics --all-targets -- -D warnings` | 0 | Formatting/lints pass |
| `python3 target/tmp/kepthx2-prepost.py` | 0 | Base/candidate verified delivery fingerprints identical |

The check harness uses `BAANDI_TEST_TMP="$PWD/target/tmp"`,
`TMPDIR="$PWD/target/tmp"`, `TMPDIR_990="$PWD/target/tmp"` and
`CHECK_SCCACHE=0` for workspace-local scratch and no shared-cache writes.
Raw logs and command exit codes: `target/tmp/kepthx2-check-*.log` and
`target/tmp/kepthx2-checks.json` (local scratch, not tracked evidence).

The additional full-workspace `scripts/check.sh` test phase reports **924
passed, 3 failed, 196 ignored across 95 binaries**. All three failures are
unchanged CLI pseudopotential-cache fixtures rejecting this task's mandated
worktree-local scratch (`pseudopotential cache must be outside repositories`):
`pseudo::cache_resolution_and_no_network_side_effects`,
`pseudo::cached_fetch_checks_hash_and_provenance`, and
`pseudo::tests::fetch_boundary_verifies_before_publication_and_never_overwrites`.
No cache rule or unrelated CLI file was changed to bypass that restriction.
The full gate is therefore not a passing workspace claim. Raw failure detail
is saved in `target/tmp/kepthx2-workspace-tests.log`. The script's optional
load probe also cannot invoke `ps` in this sandbox; it continues to the gates.

Full gate completion: `scripts/check.sh` **exit 1**; build, fmt and workspace
clippy PASS, test FAIL only for the three cache-path fixtures above. All five
ticket TEST_COMMANDS exit 0. No dev08 wall/RSS result or default promotion is
claimed; changes remain uncommitted for review.
