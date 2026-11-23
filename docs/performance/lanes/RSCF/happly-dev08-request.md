# P3-HAPPLY: reference H application profile and route audit

Status: implemented; all ticket Mac checks passed after the approved FFT activation change. No speed claim and no
new default. The master owns dev08 execution, full gates, ledger/dashboard/README
integration and the commit. No QE source was read. The source base is
`8aa5c1eec9827ac0981f0512e4f16a1a7713a315`; changes are uncommitted.

## Approved FFT activation

The user approved the activation-condition-only change to
`crates/dft-numerics/src/real.rs`. `DAVIDSON_STAGES=on` now enables the existing
FFT probe and clocks independently of FFT_PROFILE, PHASE_PROFILE and
STAGE_TIMERS. With DAVIDSON_STAGES off, the previous FFT/PHASE/STAGE switch
behavior is unchanged. When all profiling is off, probe sites retain the
cached disabled branch: no clock, TLS, counter update or allocation.
No `distributed_fft.rs` or FFT numerical operation changed.

## Measurement contract

The reference `davidson_stages.apply_h` object contains rank-local totals and
one record per completed Davidson solve, in execution order, including SCF,
retries and verification. The existing outer `over_ranks` merge takes maxima
and provides the corresponding tree of minima under `min_over_ranks`. Do not
sum different ranks' maxima to reconstruct wall time. Preserve rank-local
JSON as well as merged output for critical-rank closure.

`stage_s.total` is inclusive. Exclusive children are mirror exchange/packing,
ORC pair pack/unpack, kinetic, projector projection, projector reduction,
D coupling, and back-product plus image addition. FFT children separate local
stick (1D) work, local plane (two 1D axes) work, transpose pack/unpack,
alltoall, grid scatter/gather/scaling and local potential multiplication.
The FFT stages are the existing `fft_profile::NAMES`, not a second numerical
implementation. Their individual clocks do not include surrounding loop or
timer overhead; report the nonnegative unattributed remainder on each rank.
Mirror includes staging of outgoing ORC rows. Back-and-add includes the
existing GEMM finite-output check. No new barrier or collective is added.

Counts: `columns_per_call` is the exact sequence, `complex_pairs_per_direction`
is sum ceil(columns/2), `single_band_tails` is sum columns mod 2, and
`two_band_pairs` is their difference. Each pair uses one inverse and one forward
complex 3D FFT. Always check columns = 2*pairs - tails and calls/columns against
the parent reference records. QE's h_psi calls count blocks, not these columns.

Storage is reserved before each solve (at most max_iterations + 1 H calls),
updated in memory, and copied to the report queue only after the solve.
Diagnostics off: no added clock, FFT snapshot, counter update, or allocation
in H application; cached branches remain. This is not a measured zero-cycle
claim. Early solve errors clear the active scope; incomplete solves do not
masquerade as completed solve records. CLI drains stale records before a run.

## GEMM audit: no new switch justified

Both H nonlocal products already use `multiply_slices_wide_into` under the
adopted, default-on `BAANDI_M2_DGEMM`. The alternative `multiply_into` is
`pure_direct`. Repeating `REF_DAV_WIDE_ROTATE` here would be a redundant arm.

Let r = local ORC rows, p = projectors, b = current block columns and
s = max(p, operator capacity). The workspace is constructed with long = r,
short = s. Projection has output p x b and reduction r; p*b <= s*s, so the
whole-output fallback cannot occur. Back-application has output r x b and
reduction p; its wide Identity entry already removes the whole-output bound.
Remaining fallback candidates are `pure_panel_bound` when min(256,r)*b > s*s,
or a disabled MKL dgemm switch; a pure build necessarily reports `pure_build`.
These are not repaired by calling the same wide entry again. Zero-local-row
ranks skip GEMM. Shape/layout/long/short bounds otherwise follow construction.

`nonlocal_wide_entry_calls` and `nonlocal_direct_pure_calls` identify requests
inside these reference solves. Confirm effective backend with the existing
`gemm_dispatch` rows and pure-build `rows_as_if_mkl`; `other` includes other
callers, so do not attribute all of its process FLOPs to H. Use caller-resolved
profiles or the existing source-only GEMM-accounting build for that attribution.
The local fixture is correctness/work evidence, not a large-share profile.
No new switch or switch-count change is made before dev08 identifies a real
fallback with a material share. If it does, return the exact site/shape/route
and time to the master for a separately scoped candidate.

The old-path rejected FFT batching/threading, blocked nonlocal traversal,
kept-plane fusion and arrival-ownership ideas are not retried. Read ledger
A11/A28/A34/A43 and the LANE-FFT/LANE-NONLOCAL closing notes before proposing
follow-up changes. The present work is measurement of the reference path.

## dev08 jobs and arms

Use the exact Si64 and Si128 MPI8 inputs/UPFs from `/data/sig-dft/vc/` and their
setup/table/raw-input manifests. Resolve absolute paths on dev08; record full
source commits/diff hashes, compiler/MPI/MKL versions, build commands, binary,
input, UPF and captured-start SHA-256. Never infer grid, bands or coordinates
from the cell name. Keep thresholds, mixing, start, output, adopted verification
settings and independent verification identical. Verification/recovery cost
stays inside delivered wall. Build baseline, current main and reviewed branch
separately; report main drift separately from the candidate effect.

```sh
CARGO_TARGET_DIR=target/happly cargo build --release -p dft-cli --features mpi,mkl
```

Master runs `scripts/check.sh full` and the native rank matrix on dev08, not
Mac. One case/rank count per fast job <=600 s including output; split Si128
campaigns into bounded jobs. If a full converged pair cannot fit, request a
named authorized `final-` job <=7200 s; do not truncate SCF or verification.
Use an exclusive socket-0 window, eight physical cores 0-7, one thread per
rank, OMP_NUM_THREADS=1, MKL_NUM_THREADS=1 and established binding flags.
Record actual CPU mapping, `--report-bindings`, load, raw stdout/stderr and
exit status. No builds overlap timing.

1. One profiled full Si64 and Si128 MPI8 run per available arm. Current A is
   the unchanged default (`M2_DGEMM=on`); D is unset and is an identity control,
   not another optimization. Set `BAANDI_DAVIDSON_STAGES=on`. Retain ordinary,
   retry and verification solve splits and combined totals. A current-main
   profile supplies the pre-instrumentation baseline. At least 30 samples in
   a claimed hot function; lower counts are directional only.
2. Confirm parent/child work closure, positive FFT stage clocks on nonempty
   ranks and per-rank extrema. Compare diagnostic-on/off energy, density,
   orbitals, H work, iteration decisions and verification bit for bit.
3. Establish diagnostic-off A/A noise with disposable warmup and initially
   three crossed pairs; extend to 20 pairs for a promotion campaign. Compare
   pre-change main versus candidate A to check instrumentation-off overhead.
4. For EACH future independently implemented switch, define A=off, B=on,
   default=off; run one profiled full run per arm/cell first, then fixed-H
   cold/warm/verification replay. Only survivors proceed to crossed A-B/B-A
   full converged runs: initially three pairs, 30 for a close/surviving effect,
   interleaved A/A controls. Current new-switch list is empty: do not fabricate
   B arms or claim an A/B optimization result for identical configurations.

Launch template (substitute frozen paths and established binding wrapper):

```sh
BAANDI_DAVIDSON_STAGES=on OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /usr/bin/time -p mpirun -np 8 --report-bindings \
  /absolute/reviewed/dft-cli run /absolute/frozen/input.toml \
  --output-dir /absolute/unique/profile-output
```

Acceptance timing sets diagnostics off. Report external wall medians for both
engines, paired effects/uncertainty and A/A floor, including every failure.
No minimum gain threshold; unresolved noise is inconclusive. No default
promotion on silicon alone: include skew/elongated/slab and different projector
structures/band counts, native rank 1/2/3/8 and applicable thread determinism.
Strict accuracy class 1e-10 Ha and standard 5e-9 Ha require independent strict
oracles and passing delivered-density verification, not just mutual agreement.

## Table requested next to QE

These ticket-supplied numbers are motivation, not a new QE reference. Recover
original QE 7.6 commands, version, raw output, input and UPF SHA-256 before
publishing comparisons. Use QE MPI8 only, same physical and numerical conditions,
cores/binding and fastest valid QE configuration. No mixed parallel model.

| Quantity | BAANDI Si128 MPI8 supplied | QE MPI8 supplied | New profile |
| --- | ---: | ---: | --- |
| H apply / h_psi | 16.5–17.0 s, 14203 columns | 11.14 s, 82 calls | pending |
| local work / vloc_psi | not split | 7.98 s | pending |
| local FFT + transposes / fft | not split | 6.70 s | pending |
| projection + rank reduction / calbec | not split | 1.41 s | pending |
| coupling + back + add / add_vuspsi | not split | 1.52 s | pending |

QE child timers overlap (`fft` is inside local work); do not sum them all.
Report FFT compute and communication separately, inverse/forward directions,
pack/unpack, potential and kinetic, rank min/max, solve counts and column widths.
Compare timer scopes explicitly; these internal lines are not a delivered-wall
speed ratio. Si64 needs its own matched QE output, not scaled Si128 numbers.

## Local validation

Executed on Mac, pure backend, 2026-11-24. No release/full gate, native
MPI/MKL execution, remote installation or git commit was performed.

| Command | Exit | Evidence |
| --- | ---: | --- |
| `cargo check --workspace --tests` | 0 | Workspace/test targets compile. |
| `cargo test -p dft-physics` | 0 | 95 passed, 18 ignored, 0 failed across targets. |
| `cargo test -p dft-engine --test refscf` | 0 | 19 passed, 1 ignored, 0 failed. |
| `cargo test -p dft-numerics switches` | 0 | 14 passed, 0 failed; no registry changes. |
| `cargo test -p dft-physics --lib reference_h_apply_profile_is_bit_identical -- --nocapture` | 0 | Four comparisons to diagnostics-off preserve all eigenvalue/vector/HX bits. |
| `cargo test -p dft-engine --test refscf setup_phase_reuse_preserves_verified_science_at_both_diagnostic_levels -- --nocapture` | 0 | Four fully verified SCFs preserve scientific bits. |
| `cargo test -p dft-cli --lib davidson_stages_rank_merge_takes_max_and_min` | 0 | Nested apply_h totals and solve arrays retain correct extrema at 3 simulated ranks. |
| `git diff --check` | 0 | No whitespace errors. |

`cargo fmt --check && cargo clippy -p dft-physics -p dft-solvers -p dft-cli
-p dft-numerics --all-targets -- -D warnings` passed with exit 0. Initial workspace
check exited 101 because CLI has no direct solver dependency; fixed by a
physics-layer re-export without manifest changes. Initial clippy exited 101
for a collapsible conditional; corrected without changing arithmetic.

The new fixed-H fixture is a skew cell with three coupled projectors and a
nonconstant potential. It checks two and three bands against diagnostics-off,
using stage-only and stage+existing-FFT-diagnostics arms. It intentionally
compares the returned state at the default 20-correction budget; it is not a
converged-eigenpair or SCF accuracy claim. Both enabled arms have identical
calls/columns/pair/tail counts:

| Bands | Calls | Columns | Complex pairs per direction | Single tails | Eigenvalue delta (Ha) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 2 | 21 | 31 | 21 | 11 | 0 (all output bits identical) |
| 3 | 21 | 49 | 28 | 7 | 0 (all output bits identical) |

Both enabled arms require positive counts and clocks in all 17 applicable
FFT stages; kept_copy is intentionally absent. The stage-only arm sets
FFT_PROFILE, PHASE_PROFILE and STAGE_TIMERS off, proving DAVIDSON_STAGES
alone activates the split. The fully disabled arm also asserts that the
entire FFT snapshot is unchanged (zero added counts, work and nanoseconds). All nonlocal
products in these solves request the wide entry (42 calls per solve), with
zero direct-pure requests; this pure-backend test cannot prove native MKL use.

The independent complete small reference-SCF fixture passes its verification
in all four combinations of setup reuse and diagnostics on/off: **6 maps**,
energy **-1.10503701485504457 Ha**, delta **0 Ha**, bit-identical science
(energy, density, orbitals and work/verification fields covered by the existing
test). It uses synthetic pseudopotential data, not a QE comparison.
Scratch logs: `target/p3-happly-{check,physics,refscf,switches,identity,scf-identity,ranks,clippy}.log`.

Draft ledger entry for the master: P3-HAPPLY adds per-reference-solve H clocks,
columns/pairs and existing FFT snapshot deltas, with no numerical/default
change; two fixed-H widths and four verified SCFs preserve bits. Both nonlocal
products already use wide dispatch by default; no redundant optimization arm
added. Single-switch FFT activation is implemented in the approved real.rs scope. Native MPI/MKL performance and cross-system promotion remain pending.
Dashboard/session registry, README figures and central ledger are outside this
ticket's ownership and remain with the master; the checkout registry has no
P3-HAPPLY run entry. No recurring reporting timer was installed.

### Approved activation follow-up

All ticket Mac commands were rerun after the activation change and the
disabled-probe snapshot assertion; every exit code was 0. Clippy additionally
covered dft-numerics because real.rs is now touched. The targeted fixed-H
and fully verified SCF identity tests were rerun with unchanged counts and
energy above. Exact final command/exit records are in
`target/p3-happly-final-commands.json`. No ownership blocker remains; dev08
measurement is the next master-owned step. Changes remain uncommitted.
