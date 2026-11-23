# P3-DAV-DENSE: reference Davidson dense work

Status: candidate implementation, default off, dev08 measurement pending.
No performance acceptance is inferred from Mac tests. The master session owns
remote execution through `baandi-dev08`, full gates, rank matrix, promotion,
ledger/dashboard/README integration and the commit. No QE source was used.
[Cost accounting](../../../design/reference-scf.md#inner-solver-cost-accounting).

## Premise and mechanism

Ticket-relayed main `2aee083` MPI8 profile, median of three, maximum rank,
SCF plus verification solves: Si128 restart 8.500 s, border 5.495 s,
corrections 5.953 s, return 3.874 s, factor 0.384 s, reduced eigen 2.647 s,
H apply 16.850 s; 20 solver calls, 83 corrections, 14203 H columns.
Si64 restart/border/corrections/return: 0.696/0.707/0.906/0.378 s.
These are supplied motivation, not newly reproduced numbers. Preserve the
original raw profiles/hashes before using them in a final comparison.

The current reference solver already selects active bands, borders only new
columns, rotates V/HV at restart and does only the necessary final rotations.
It nonetheless multiplies by identity C in the correction immediately after
restart. The candidate omits those two GEMMs and coefficient packing, using
the selected retained V/HV columns directly. All preconditioning, norm/drop
rules, repaired-rank logic, active decisions and fresh H applications remain.
Per affected correction: 4*r*n*u fewer FLOPs per rank and
16*(r*n + n*u + r*u) + 16*n*u fewer logical operand/copy bytes.
There is no material-dependent knob or threshold change.

This is not A52 CholeskyQR2, hard locking, smaller capacity, REL-SLICES,
Advance, or the rejected upper-border tiling (resolved Si64 slowdown). The
reference solver's previously rejected exact diagonal-H/unit-S restart reset
is not repeated: measured S/H refresh and factor remain. The private
reference workspace is outside the original OWNED_FILES, so no buffer-swap
or reduction-fusion change is silently made there or to the old workspace.

## Arms, build and frozen inputs

One independent switch is implemented. Consequently the all-candidates arm
is B1 itself; do not invent or time duplicate configurations as independent
observations.

| Arm | BAANDI_REF_DAV_RESTART_CORRECTIONS | Meaning |
| --- | --- | --- |
| D | unset | candidate default, currently off |
| A | off | explicit baseline arithmetic |
| B1 / ALL | on | omit the post-restart identity rotations |

Input key is `[performance] ref_dav_restart_corrections='on'|'off'`.
Environment/input precedence and MPI agreement hash use the existing registry.
Resolve switches on every rank before entering the solve. Confirm consulted
state in the result manifest. Keep all other effective switches equal,
especially adopted `REFSCF_VERIFY_GATE=on`, `REFSCF_VERIFY_REUSE=off`, and
`REFSCF_VERIFY_GATE_FIRST=off`. Do not use the historical both-off gate arm
from the original verify-cost proposal as current default.

Use the **exact Si64 and Si128 MPI8 input/UPF files from the verify-cost A/B**:
remote provenance `/data/sig-dft/vc/{si64,si128}/setup.md`, `TABLE.md`,
`table.json` and their raw-input manifests. Do not regenerate structures or
infer FFT grids/bands from atom counts. The implementer cannot read those
remote manifests: the operator must resolve absolute paths and record SHA-256
for input, UPF, start/snapshot and any reference result before launching.
Retain identical coordinates, cell, cutoff, FFT grid, band count, atomic start,
reference driver, thresholds, mixing, verification and output options.
Independent fixed-point verification remains on and charged to every full run.

Build baseline `2aee083`, current main and reviewed candidate separately using
the same pinned compiler and `mpi,mkl` release features. Record full commits,
source diff/tree hash (candidate is uncommitted here), actual build commands,
MPI/MKL/compiler versions and each binary SHA-256. Example candidate build:

```sh
CARGO_TARGET_DIR=target/dav-dense cargo build --release -p dft-cli --features mpi,mkl
```

Use the established dev08 build environment; no installs or remote changes are
requested. Run the full gate on socket 1 through the master, serialized with
measurement jobs. Run `scripts/check.sh full` and the existing reference rank
matrix there, never on Mac. Compare pre-change main against candidate A and D
for bit-identical science and work; instrumentation clocks/manifest additions
are the only exclusions. If current main moved, report the baseline-to-main
movement separately from candidate effect.

## Sequence and timing protocol

1. Correctness replay first: capture one warm reference solve that actually
   restarts on each cell, plus a cold and a late/verification solve. Use the
   established reference fixed-H replay harness and identical captured inputs.
   Verify exact A replay before A/B. Check incoming active sets, every Ritz
   operand/threshold, H calls/columns, repairs/drops/restarts and actual GEMM
   savings. A replay without a restart is a negative control, not a speed test
   of this change. Keep final-iteration-budget exits as explicit nonconvergence.
2. In one exclusive socket-0 window, MPI8 on the same eight physical cores
   0-7, one thread per rank, library threads 1. Capture actual CPU mapping and
   `--report-bindings`, host load and versions. No builds overlap measurements.
3. Establish A/A floor before A/B: disposable warmup, initially three crossed
   pairs (A0-A1, A1-A0) per case. Accumulate 20 A/A pairs if progressing to
   full acceptance, interleaving controls in subsequent windows.
4. B1/ALL independently against A: initially three crossed A-B1 / B1-A pairs,
   then 30 crossed pairs/case for a close or surviving effect, using the same
   frozen full-SCF inputs. Do not select the best sample or fastest subset.
5. Separate profiled full run for **each arm A and B1/ALL**, with
   `BAANDI_DAVIDSON_STAGES=on` and all other settings identical. Timing samples
   use `BAANDI_DAVIDSON_STAGES=off`. Profile ordinary and verification calls
   separately and combined. D identity only needs a correctness control.

One case and rank count per job, fast class <=600 s including output. Si128
starts with one pair/job; split the campaign into bounded jobs. If a complete
convergence/verification pair cannot fit, request a named authorized `final-`
job (<=7200 s) through the master; never truncate the tail to fit. Failed or
censored samples remain failures and are retained.

Launch template (substitute frozen absolute paths and the existing binding
wrapper, never change the scientific input between arms):

```sh
BAANDI_REF_DAV_RESTART_CORRECTIONS=off BAANDI_DAVIDSON_STAGES=off \
  OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /usr/bin/time -p mpirun -np 8 --report-bindings \
  /absolute/reviewed/dft-cli run /absolute/frozen/input.toml \
  --output-dir /absolute/unique/A-output
```

Actual mapping/binding flags must come from verify-cost setup and be recorded.
Replace only the candidate switch for B1 and only stage profiling for profile
runs. Time the whole launcher process including output, failed verification
attempts and recovery, and record exit codes. Report both engine medians,
median ratio, paired log-effect with uncertainty, A/A floor and raw samples.
No minimum percentage gain; unresolved noise is inconclusive. This request
compares BAANDI MPI8 arms only. Any later QE table must use QE MPI8 under
identical physical/input/threshold/core conditions, documented 7.6 provenance
and the same delivered-wall statistic, never mixed parallel models.

## Evidence and acceptance

Required table per cell and arm: external wall (s); SCF and verification map
counts; solver calls, Davidson iterations, H columns, restarts, repairs,
dropped corrections; active-width sequence; actual GEMM FLOPs and collective
counts/payload; eight stage seconds (apply_h, border, factor, reduced_eigen,
corrections, repair, restart, return). Retain per-rank times, rank maxima and
critical-rank closure separately; summing unrelated maxima is not wall.
First-order expectation is reduced corrections FLOPs/time, unchanged H work,
unchanged restart/border/return operations; reductions stay identical.

Preserve results JSON, density.bin, full profiled Ritz traces, scientific
fingerprint, commands, stdout/stderr and raw paths/SHA-256. Require exact
fixed-layout decision/H-column counts. For a differing predicate, retain
first map/phase/call/correction/band, both previous/current Ritz values,
changes and actual thresholds, signed margins in Ha. Apply the existing
rank-rule envelope `32*EPSILON*max(1 Ha, abs(four Ritz operands))` to every
straddling band; no missing earlier trace or unrelated difference is excused.
No new rounding exception is asserted by the Mac fixtures.

Every full run must converge, independently verify its delivered density,
and pass the unchanged strict 1e-10 Ha / standard 5e-9 Ha accuracy classes
against converged strict oracles with recorded provenance. A/B agreement
alone is not an oracle. Verify charge and finite outputs under the existing
rank-rule checks. Confirm candidate A and D preserve old behavior, and the
legacy distributed solver is unaffected even with the new switch on.
Run the established rank matrix (1/2/3/8), different shapes/band counts and
closed-shell synthetic projector fixtures before promotion. A Si64/Si128 win
alone cannot enable this switch by default. No standalone cost threshold
rejects a reproducible small gain with passing correctness.

## Local validation and ledger draft

Local commands/results will be appended after execution. Only the ticket's
Mac compilation, targeted tests, fmt and clippy are run; no full check or
release matrix. Tests cover nondegenerate/clustered spectra, cyclic row
partitions at 1/2/3/8 simulated ranks, root/replicated EVD, one-step exhaustion
and converged restarts. Simulation does not establish native MPI/MKL timing
or backend bit identity.

Draft for the orchestrator-owned attempts ledger (fill dev08 cells only from
actual measurements):

| ID | Mechanism | Status | Evidence | Lesson / limits |
| --- | --- | --- | --- | --- |
| P3-DAV-DENSE | Reference post-restart corrections directly from retained V/HV; identity C rotations and selected-coefficient packing omitted | Implemented, default off; dev08 pending | 4*r*n*u FLOPs removed per affected correction; unchanged measured Gram refresh; local bit/work tests below | Already-active-only corrections need no locking. Exact H/S reset remains rejected. Native MPI/MKL, multi-system full verification and delivered wall still required; no Mac speed claim. Derived from davidson-iteration-structure-gamma.md and engine inspection. |

### Executed Mac checks (2026-11-24)

Base HEAD: `07389881f37ea1d42811884ec33306c16a91d1d7`; the reference solver
and its workspace are unchanged between ticket-profile `2aee083` and this
HEAD (`git diff 2aee083..HEAD --` for those two paths is empty). Candidate
changes remain uncommitted. Original ownership was sufficient for the
implemented mechanism; no expanded workspace ownership is required.

| Command | Exit | Result |
| --- | ---: | --- |
| `cargo check --workspace --tests` | 0 | All workspace/test targets compile. |
| `cargo test -p dft-solvers` | 0 | 161 passed, 25 ignored, 0 failed across targets. |
| `cargo test -p dft-engine --test refscf` | 0 | 16 passed, 1 ignored (external/full-run fixture). |
| `cargo test -p dft-numerics switches` | 0 | 14 passed, 0 failed; registry count 137. |
| `cargo fmt --check && cargo clippy -p dft-solvers -p dft-numerics --all-targets -- -D warnings` | 0 | Formatting and all-target warnings-as-errors pass. |
| `git diff --check` | 0 | No whitespace errors. |
| `cargo test -p dft-solvers --lib restart_corrections_preserve -- --nocapture` | 0 (final) | Counter/bit-equivalence fixture passes. Two initial development invocations exited 101: a missing helper closing brace, then use of an environment name instead of a lowercase RunConfig input key. Both corrected before the full checks; neither was a numerical mismatch. |
| `BAANDI_REF_DAV_RESTART_CORRECTIONS=on target/debug/deps/ref_davidson-4f96658c9cc927aa --nocapture` | 0 | Existing compiled debug oracle target: all 5 tests pass with candidate enabled, including dense/legacy and complete degenerate occupied-edge checks. |

Scratch logs: `target/tmp/dav-dense-{solvers,refscf,switches,clippy,enabled-oracles}.log`.
The oracle executable suffix is build-specific; the reproducible Cargo
invocation is `BAANDI_REF_DAV_RESTART_CORRECTIONS=on cargo test -p dft-solvers
--test ref_davidson -- --nocapture`. No release execution was used.

The new unit test covers **48 off/on layout-option pairs**: three spectra and
band/dimension combinations, four rank counts, root/replicated EVD and one-step
or converged budgets. It verifies eigenvalue/vector/HX bits, Ritz operands,
active selections, iteration counts, per-step H columns, repairs, drops,
restarts, reduced orders, collectives and exact eliminated FLOP/operand-byte
formulas. The 7-row case forces repair/drops and includes an empty rank at
np8. The one-step controls never take the new branch. No divergence or new
rounding exception occurred. The test takes about 0.05 s, not a speed claim.

| Synthetic dimension / bands | Completed correction iterations (both arms) | H columns (both arms) | Restarts (both arms) | Eliminated nominal FLOPs on rank 0 at ranks 1 / 2 / 3 / 8 | Max off/on eigenvalue delta (Ha) |
| --- | ---: | ---: | ---: | --- | ---: |
| 12 / 3, nonclustered | 4 | 13 | 3 | 1008 / 504 / 336 / 168 | 0, bit-identical |
| 17 / 5, clustered | 4 | 21 | 3 | 3740 / 1980 / 1320 / 660 | 0, bit-identical |
| 7 / 4, capacity/repair | 2 | 7 | 1 | 448 / 256 / 192 / 64 | 0, bit-identical |

These are cyclic partitions, so rank-0 rows are ceil(dimension/ranks); these
FLOPs are not global sums. Both EVD placements have the same table. The
separate enabled independent dense oracles (32 rows, 8 bands, ranks 1/2/3/8,
both placements) have maximum eigenvalue error **7.327e-14 Ha** and maximum
occupied-projector error **2.392e-9** (dimensionless), with the existing
unchanged tolerances. They test absolute eigenpairs, not SCF fixed-point
accuracy. Native MPI/MKL equivalence and speed remain pending.

No remote action, full `check.sh`, release matrix, UPF installation, commit,
or out-of-scope README/dashboard/session/ledger edit was performed. Those
integration artifacts remain with the master session. Draft ledger evidence
can now cite 48 bit-identical pairs and 4rnu saved FLOPs per affected
correction; delivered-wall cells remain pending.
