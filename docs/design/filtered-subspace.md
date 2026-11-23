# G-FILTER-REPLAY: resident-block filter, stage 1

Decision: **STOP production integration of this prototype.** The best fixed-degree
11-map estimate at Si64 MPI8 saves **-9.8693 s optimistically** (a loss), with a
conditional pessimistic estimate of **-35.5381 s**. These are cost-model outputs,
not measured MPI timings. All five degree policies have negative optimistic
savings. Correctness also prevents adoption: the captured early-map stopping
profiles do not imply the absolute frozen-band-energy gates, for either solver.
No tolerance, SCF path, switch default, or production caller changed.

## Delivery boundary and provenance

The ticket owns `crates/dft-solvers/examples/davidson_bench.rs`, which does not
exist. The actual harness is `crates/dft-cli/examples/davidson_bench.rs`. Ownership
correction was requested; no answer has arrived. Consequently the new mode is a
reviewable [bench.patch](../performance/lanes/FILTER/bench.patch), not an edit to
that unowned file. [build_replay.py](../performance/lanes/FILTER/build_replay.py)
applies it only to scratch and builds the prepared example against the existing
CLI dependencies. It was compiled, linted, and exercised. The solver module and
its declaration are installed; neither has a production caller. The overall
ticket remains blocked on this integration boundary, despite a scientific STOP
verdict. Do not describe the checked-in CLI as supporting `--filter-replay` yet.

`gh issue view 24` exited 1: API connection unavailable. The web fallback also
failed. The issue body was **not read**; the issue's scope/exit rule here comes
from the user-supplied ticket. Read locally: architecture-redesign brief
2026-11-23, Si64 gap attribution 2026-11-20, DAVIDSON closing report, lane charter,
fast-loop snapshot contracts, performance brief/backlog, ADRs 0003/0004 and the
Davidson ledger rows. No QE source or remote machine was accessed.

The new premise relative to the older Chebyshev warning is the measured MPI8
algebra budget (2.953 s): a fixed-width resident block might remove growing
subspace Gram/eigen/rotation work. It does not retain restart directions, truncate
SCF maps as Advance did, slice an EVD as REL-SLICES did, or minimize the nonlinear
energy with PR-CG. Those rejections remain closed. Existing wide local kernels
and Gamma pairing are retained; the rejected UPPER_BORDER/ROTATE_COLUMNS GEMM
shape experiments are not repeated.

## Numerical contract

Input is the captured warm ORC block, current fixed Hamiltonian, optional kept
HX, and the exact captured tolerance/strict prefix/buffer tolerance/Ritz-change
flag. Input row ordering is unchanged. The module uses `DistributedRealOperator`;
the serial replay adapter calls the existing `RealHamiltonian::apply_block`, and
the distributed adapter calls `DistributedRealHamiltonian::apply_block`.

1. Check the start Gram defect against the existing 1e-12 rule. Reuse HX only
   for an unchanged orthonormal start and validated operator; otherwise QR and H.
   A setup Rayleigh–Ritz (RR), including both rotations, is charged separately
   from subsequent sweeps.
2. Estimate the upper edge by 12 shifted power steps on two deterministic
   globally indexed random columns, preserving paired application. This costs
   **24 H columns** per solve, plus reductions/normalizations. Lower shift is
   the smallest warm Ritz value minus the maximum warm residual. Upper padding
   is a dimensionless 20%; this is an **estimate, not a rigorous bound**. The
   independent strict reference remains necessary. Cutoff is the largest warm
   Ritz value plus the greater of its residual and the mean spacing across the
   warm Ritz spectrum. None of these choices depends on atom type/count.
3. With c=(upper+cutoff)/2, e=(upper-cutoff)/2 and alpha=(lower-c)/e,
   compute the scaled three-term recurrence. sigma_1=1/alpha and
   sigma_(j+1)=1/(2 alpha-sigma_j). Two reusable recurrence buffers hold the
   previous/current block; HX is the application scratch. The first polynomial
   action uses the resident HX. There are no allocations in the recurrence.
4. QR uses the existing two-pass panel Gram-Schmidt routine, after column
   normalization, with its shared 1e-12 dependence threshold. This is CGS2 QR,
   not CholeskyQR2. It avoids reopening the rejected correction-block CholQR
   experiment, but its many reductions and streaming passes are charged.
5. Apply H to Q, form Q^T H Q, solve one order-B symmetric problem through
   `RealEighWorkspace`, rotate Q and HQ through the existing wide kernels, and
   measure residuals. A degree-d sweep therefore costs d full-width H blocks
   including the post-QR application. One extra true-residual H application is
   charged at a successful strict return. All matrices/journals are preallocated
   before the sweep loop; backend-owned EVD packing can still allocate.
6. Standard thresholds exactly use the captured strict prefix minus the first
   empty mandatory root, which uses the buffer tolerance. All roots must pass;
   strict returns additionally recompute true residuals. Successful returns
   require a separated occupied edge (Ritz gap exceeds the two edge residuals).
7. Nonfinite output, rank loss, invalid interval, unresolvable degree/edge,
   three consecutive nonimproving all-band residual maxima, or 24 sweeps request
   fallback. The harness **executes** Davidson again from the original captured
   start/profile; it does not substitute a previously computed result. Its work
   and every attempted filter sweep remain in the charged totals. This conservative
   restart avoids pretending that unusable filter vectors are a free warm start.

Explicit degrees 4, 8, 12, 16 are experiment arms, not recommended defaults.
`auto` derives an even degree from acosh(2 max(error/tolerance,1)) divided by
acosh((c-occupied_edge)/e), with a safety budget of 64. Using the captured Ritz
threshold as an attenuation target is conservative: auto declines 13/15 sampled
maps before a sweep and falls back. No material-dependent retuning was done.
The sweep log exposes every chosen degree and decline. An estimated bound does
not establish spectral completeness, even when residuals are small.

## Capture and replay evidence

All five default SCFs converged and produced snapshots after creating the
snapshot directory explicitly. The first capture attempt omitted that mkdir;
the CLI exited 0 but logged write failures. It was discarded and repeated.
[Capture index](../performance/lanes/FILTER/capture-index.json) stores manifest
hashes, source/binary/UPF hashes, selected maps and exact group frequencies.
Binary snapshot files and UPFs are not committed; snapshots are retained under
`target/scratch/<case>-snapshots/` in this worktree.

Fine ordinary maps were split chronologically into three groups; first early,
central middle and last late maps were replayed. Frequencies are actual group
sizes, not three equally weighted solves. The 11-map estimate multiplies each
sample by 11 times its group size divided by the actual fine-map count. This is
explicitly a **three-snapshot approximation**, not a full candidate SCF or an
exact sum of all map costs. Coarse, redo and verification maps are not replaced.

| Case | Actual fine maps | Selected map numbers | Actual early/mid/late frequencies |
|---|---:|---|---|
| Si8 | 11 | 1 / 7 / 11 | 4 / 4 / 3 |
| Si16 | 10 | 1 / 6 / 10 | 4 / 3 / 3 |
| Si32 | 10 | 1 / 6 / 10 | 4 / 3 / 3 |
| Si64 | 13 | 1 / 8 / 13 | 5 / 4 / 4 |
| Si skew | 8 | 1 / 5 / 8 | 3 / 3 / 2 |

Every serial engine replay and reconstructed serial Davidson reproduced the
captured H count and eigenvalue FNV exactly. Reference is a separate fixed-H
Davidson residual solve at **1e-11 Ha for every band**, not QE. Principal angles
are computed from singular values of `(I-X_ref X_ref^T) X` through its small
Gram eigenproblem; no bandwise eigenvector comparison is used. Density comes
from unpacked occupied ORC vectors with occupation 2 and the existing FFT reader
of the basis, in electrons/Bohr^3. Its L2 error is sqrt(integral delta-rho^2 dV),
in electrons/Bohr^(3/2). The angle test includes an arbitrary rotation of a
fully degenerate occupied block and a known 0.6-radian leakage angle.

The reported energy is **2 times the sum of occupied fixed-H eigenvalues**.
It is not an SCF total energy and is not a fixed-point certificate. We apply
1e-10/5e-9 Ha gates to that reported number and retain failures. Across the five
cells and five policies: early standard gate 5/25, middle 20/25, late 25/25;
strict band-sum gate early 2/25, middle 4/25, late 22/25. Thus there is no claim
that the requested absolute accuracy acceptance criteria are fully met. For
example, Si64 early Davidson itself has band-sum error 6.38157e-5 Ha; degree 4
has 4.85340e-5 Ha. The late Si64 errors are 6.21e-12, 9.22e-13, 1.24e-13 and
8.90e-15 Ha for degrees 4/8/12/16. Any additional terminal refinement needed to
meet an absolute gate would add cost; none is silently supplied for free.

Simulated ranks 2 and 4 exercised the existing M2 Hamiltonian, collectives and
fallback. All rank-local records are returned in `per_rank`; reductions include
count and doubles payload, with H-internal collectives in the wrapper totals.
The engine exact replay needed `BAANDI_M2_LOCAL_FIELDS=off`: capture disables
owned local fields, while default replay otherwise hashes the owned/sentinel
layout and rejects it. The standalone reconstructed MPI operator uses serial
Hartree rounding, so its potential bits differ; maximum HX difference was
3.89e-16/2.60e-16 Ha. It discards stale HX and is explicitly marked
`rebuilt_baseline_exact=false`, while the original engine replay remains exact.
These are rounding-class algorithm comparisons, **not exact MPI filter replay**.
An exact production MPI operator hook would require separate engine ownership.

[replay-summary.csv](../performance/lanes/FILTER/replay-summary.csv) is the plain
scientific table. Individual `*-early/mid/late.json` contain complete H-width
histograms, QR/Gram/rotation FLOPs, RR orders, bounds, sweep/fallback journals,
principal angles and density errors. Fallback Davidson exposes count and sum of
order cubed through existing counters; its full order histogram is unavailable.
Mac wall fields are labeled observations only, including simulated-rank walls.
Optional `--filter-profile` observes disjoint setup/bounds/filter/QR/terminal/RR
regions via a replay callback; clocks stay in the bench. H time by block width
is inclusive within those phases and fallback, never an extra additive cost.
The profile-on/off smoke has identical science, counts, residuals and sweep
journals; `profile-smoke.json` stores its result. Ordinary calls have no observer.

### Measured H columns per snapshot (filter plus any executed fallback)

| Case | Phase | Davidson | d4 | d8 | d12 | d16 | auto |
|---|---|---:|---:|---:|---:|---:|---:|
| si8 | early | 114 | 1324 | 1004 | 1004 | 1004 | 158 |
| si8 | mid | 20 | 104 | 184 | 264 | 344 | 44 |
| si8 | late | 42 | 264 | 824 | 1224 | 1304 | 66 |
| si16 | early | 194 | 2076 | 1788 | 1788 | 1788 | 254 |
| si16 | mid | 56 | 312 | 600 | 888 | 1176 | 80 |
| si16 | late | 93 | 1176 | 1176 | 2184 | 2328 | 117 |
| si32 | early | 403 | 4716 | 3900 | 3356 | 3356 | 495 |
| si32 | mid | 112 | 840 | 1112 | 1656 | 2200 | 136 |
| si32 | late | 98 | 568 | 1112 | 1656 | 2200 | 122 |
| si64 | early | 593 | 5436 | 4380 | 4908 | 6492 | 12828 |
| si64 | mid | 289 | 1608 | 3192 | 4776 | 6360 | 313 |
| si64 | late | 235 | 1608 | 2136 | 4776 | 6360 | 259 |
| si-skew | early | 45 | 416 | 352 | 416 | 416 | 77 |
| si-skew | mid | 8 | 56 | 88 | 120 | 152 | 472 |
| si-skew | late | 9 | 56 | 152 | 216 | 280 | 33 |

## Break-even and the Si64 MPI8 exit rule

Relayed dev08 prices, from the requested brief/attribution and DAVIDSON closing
report: exclusive Gram 1.002 + eigen 0.988 + rotation 0.963 = **2.953 s**, across
**all stages**, including **183** reduced solves; H 2.782 s / 6025 columns =
**0.000461743 s/column**. Fine stage has **3504 columns**. DGEMM envelope
64–108 GFLOP/s and STREAM 12.97 GB/s/core are relayed, not measured here.
The 0.573 s kept-update stage is not additionally credited: the filter still
needs the current-H warm state. No unexplained outside-solver gap becomes free
setup savings.

Net saving = avoided exclusive Gram/eigen/rotation/setup seconds minus
incremental H (including bounds, post-QR/strict terminal products and fallback),
QR, new RR/rotations, scalar recurrence/residual/bound work, state copies,
workspace initialization and solver reductions. Initial RR and all failed work
are included. Reference solves are oracle-only and excluded from both arms.
Coarse and verification work remain unchanged. No comparison crosses parallel
models; Mac **counts**, not timings or speed ratios, parameterize this MPI8 model.

The optimistic case deliberately over-credits all 2.953 s to the fine-stage
replacement, even though coarse/verification contribute to it. It prices all
candidate H columns at 0.893333 of the measured average (the entire 26% nonlocal
share gets the 64-to-108 GFLOP/s improvement; the 74% FFT/local share stays).
QR is allowed the 108 GFLOP/s arithmetic floor with ideal cache reuse. Other
local algebra is divided by eight ranks at 108 GFLOP/s. Explicit copies/zeroed
workspace use 12.97 GB/s/core. New EVD cost is 0.988 times candidate sum(m^3)
over `183*(3B)^3`, using the largest permitted MPI order as the baseline divisor.
This is an optimistic **order-envelope estimate**, not a claim that syevd runs
at a DGEMM rate. Solver reductions cost an assumed 1 microsecond each.

The conditional pessimistic scenario uses 1.2 times average H-column cost,
64 GFLOP/s, half the single-core STREAM rate, streaming QR traffic
`40*N*B*(B-1)+48*N*B` bytes per QR (global), 10 microseconds per solver sum,
and EVD baseline divisor `183*B^3`. Local work is divided by eight. Reduction
latencies/cache effects and bandwidth sharing must be priced on dev08; this
scenario is a sensitivity calculation, not a guaranteed worst-case bound.
Fallback maps lose their proportional dense-algebra credit and retain all H
cost. Opaque backend allocation costs and unperformed accuracy-repair solves
are nonnegative unpriced terms; giving them zero cost only favors the filter.
Neither can rescue a negative optimistic result.

| Policy | 11-map charged H columns | New RR calls | Optimistic net saving (s) | Conditional pessimistic (s) |
|---|---:|---:|---:|---:|
| 4 | 33883.38 | 73.62 | -9.8693 | -35.5381 |
| 8 | 36564.00 | 44.85 | -10.8026 | -27.3088 |
| 12 | 53094.46 | 44.00 | -17.6246 | -36.1870 |
| 16 | 70518.46 | 44.00 | -24.8206 | -45.8464 |
| auto | 56208.31 | 19.46 | -20.5793 | -31.4433 |

The least negative fixed policy is d4: the optimistic H increment alone is
12.3702 s; QR 0.1210, Gram/rotate 0.2232, recurrence/scalar 0.0059, EVD 0.0147,
copies/initialization 0.0540, solver sums 0.0332 s. Total new/incremental work
12.8223 s minus the very generous 2.953 s credit gives **-9.8693 s**.
Even replacing every one of the 11 maps by the cheapest observed Si64 snapshot
(1608 columns) costs more than the entire algebra credit on H alone. Therefore
oversampling the expensive first fine map is not the sole reason for rejection.

The sampled serial Davidson 11-map estimate is 4282.38 columns, whereas the
actual measured MPI8 fine budget is 3504. They are reported separately; the
model uses the latter for its baseline. A measured MPI8 capture may alter the
counts, but it must overturn an order-of-magnitude H penalty, not resolve a
small timing ambiguity. Stop integration under the issue's optimistic-exit
rule. Retain only the replay experiment and evidence; no new default proposal.

Recompute all arithmetic with
`python3 -I docs/performance/lanes/FILTER/summarize.py`. The full breakdown is
[summary.json](../performance/lanes/FILTER/summary.json).

## Checks, unresolved items and next boundary

- Solver unit tests: 4 passed, including ranks 1/2/5/25 (empty local rows),
  degrees 4/8/12/16, degenerate roots, strict/Ritz profiles, rank-loss fallback,
  and scale/shift-invariant degree selection. No tolerance was relaxed.
- Prepared bench angle oracle: 1 passed; prepared-source clippy `-D warnings`
  passed. The actual CLI file is untouched pending the ownership correction.
- Full check: exit 1; build/fmt/clippy passed, **799 passed, 3 failed, 181 ignored**.
  Three CLI pseudo-cache tests reject the required in-repository TMPDIR, as the
  ticket anticipated. The load probe also cannot run `ps` in this sandbox.
- Five default SCFs before/after: energies and density SHA-256 identical; all
  verified. [default-identity.json](../performance/lanes/FILTER/default-identity.json).
- Full MPI8 fixed-operator identity, actual MPI8 filter prices, full fallback EVD
  order histograms, rigorous spectral enclosure, and total-SCF/certificate
  acceptance are **not established**. This failed cost screen does not justify
  production integration to obtain them.
- Dashboard/session registry and README generators are outside this ticket's
  owned files. No registry entry for G-FILTER-REPLAY or delegate report exists in
  this checkout. The master must own those updates. No commit was made.

The [dev08 request](../performance/lanes/FILTER/dev08-request.md) is a bounded
counter/price falsification protocol, not a request to integrate or benchmark
full candidate SCFs. Read issue #24 in a network-enabled master session before
any follow-up beyond this ticket.
