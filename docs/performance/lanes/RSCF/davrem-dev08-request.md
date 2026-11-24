# P3-DAVREM: reference Davidson remainder

Status: default-off candidates; Mac correctness/counter checks only. Master
reviews/commits and requests remote execution through `baandi-dev08`. No QE
source was read. Derived from the own-words Davidson GEMM/iteration notes;
implementation against this engine's workspace. See
[accounting and inventory](../../../design/reference-scf.md#p3-davrem-full-diagnostic-remainder-attribution-and-copy-removal).

## Question and arms

Ticket-relayed si128 MPI8 QE-matched clock: remainder 9.25 s / 18 solves
versus QE 5.29 s / 16 solves (approximately 0.51 versus 0.33 s/solve).
These are one-run diagnosis numbers, not reproduced acceptance results.
The GEMMs already use adopted WIDE_ROTATE/SINGLE_GEMM. New attribution
separates scalar preconditioning, normalizing, appending, H validation,
small-matrix transforms/coefficients, borders, factor, restart/return and glue.

| Arm | REF_DAV_SWAP_BUFFERS | REF_DAV_LEAN_RESET | Mechanism |
| --- | --- | --- | --- |
| A | off | off | unchanged arithmetic and memory operations |
| B | on | off | exchange rotation/basis buffers after restart and repair |
| C | off | on | omit five overwritten small-matrix resets per restart |
| D | on | on | both independent removals |

Switch names above have the `BAANDI_` prefix in the environment, or lower-case
names under `[performance]`. Both default off. There is no new tuning knob.
B removes 16*r*n destination-copy bytes per restart (32*r*n logical traffic),
and 16*r*m destination bytes per repair. C removes 40*n*n written bytes per
restart. GEMMs, H applications, collectives and every band decision must match A.
This does not repeat the rejected GEMM tiling, changed rotation shapes, or
idealized unmeasured restart Gram. All measured refreshes remain.

## Development protocol: exactly one profiled run per arm and cell

| Case | MPI model | A | B | C | D |
| --- | --- | --- | --- | --- | --- |
| si64 | 8 ranks, 1 thread/rank | 1 | 1 | 1 | 1 |
| si128 | 8 ranks, 1 thread/rank | 1 | 1 | 1 | 1 |
| elongated | 8 ranks, 1 thread/rank | 1 | 1 | 1 | 1 |
| rskew | 8 ranks, 1 thread/rank | 1 | 1 | 1 | 1 |

Sixteen runs total. No warmup runs, A/A, repetitions, medians, crossed pairs
or automatic replacement of failed samples in this development request.
One case/arm per fast job, maximum 600 s. Timeout/failure is preserved with
its exit code. Precise consolidated-default measurement is separate work.

Use the frozen converging reference inputs/UPF/start manifests of the current
P3 campaign, under `/data/sig-dft/vc/` and `/data/sig-dft/p3g/` as appropriate.
Resolve actual paths and SHA-256 on dev08: cell labels are not specifications.
Preserve cell, coordinates, cutoff, FFT grid, UPF, bands, atomic start,
threshold, verification policy and output policy. Independent verification
and all recovery stay inside the delivered wall. Keep all effective switches
identical except the two candidate settings; explicitly keep adopted
REF_DAV_WIDE_ROTATE, REF_DAV_SINGLE_GEMM and REF_DAV_RESTART_CORRECTIONS on.
FFT_FUSED_ROWS and REFSCF_VERIFY_ENERGY remain off. Record the full manifest.

Build one reviewed candidate binary using the established native flags and
MPI/MKL environment, no installs:

```sh
CARGO_TARGET_DIR=target/p3-davrem cargo build --release -p dft-cli --features mpi,mkl
```

Record compiler/MPI/MKL versions, build exit, full source commit and patch
hash, binary SHA-256, frozen input/UPF/start hashes, environment and exact
commands. Master handles required full dev08 gates before the screen; Mac
checks here are the ticket's limited commands only. No builds or other
measurements overlap. Use the established socket-0 physical cores 0-7 and
capture actual bindings; `--map-by core` alone is not binding proof.
Invocation skeleton (master resolves frozen input and separate output paths):

```sh
BAANDI_REF_DAV_SWAP_BUFFERS=off BAANDI_REF_DAV_LEAN_RESET=off \
BAANDI_REF_DAV_WIDE_ROTATE=on BAANDI_REF_DAV_SINGLE_GEMM=on \
BAANDI_REF_DAV_RESTART_CORRECTIONS=on \
BAANDI_DAVIDSON_STAGES=on OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
mpirun -np 8 --map-by core --bind-to core --report-bindings \
  <reviewed-binary> run <frozen-reference-driver-input-with-unique-output-path>
```

B/C/D change only their two values and output destination. Capture the
standard caller-resolved CPU stack profile during this same execution,
especially hypot, matrix copies/zeros, GEMM validation/packing, shared factor
and collective waits. State profile tool/rate and sample counts; do not make
quantitative claims from fewer than 30 relevant samples. If the established
profiler cannot attach, retain counters and report that limitation without
adding another execution.

## Required report, including per-stage throughput

Save raw stdout/stderr, exit codes, external wall, results.json, density and
orbital/eigenvalue fingerprints, and the existing strict-reference provenance.
A/B/C/D must retain maps, retries/verification outcomes, active-band sequences,
Ritz values/threshold operands, repairs/drops/restarts, H calls/columns and
GEMM shapes. Expect bit identity; any difference is a correctness failure
requiring investigation, never explained away as a faster trajectory. Keep
strict 1e-10 Ha / standard 5e-9 Ha and independent verification gates intact.

For **each case/arm/stage**, populate this table from the single run:

| case | arm | stage | rank0 s | max-rank s | nominal GFLOP rank0 | copied B rank0 | logical streamed+operand B rank0 | GEMM calls/shapes | collective calls/payload B | rank0 GFLOP/s | rank0 logical GB/s | max per-rank GFLOP/s | max per-rank logical GB/s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| pending | A/B/C/D | all 19 stages | pending | pending | pending | pending | pending | see histogram | pending | pending | pending | pending | pending |

Extractor (results file path as the first argument), no repeated run needed:

```python
import json, sys
w = json.load(open(sys.argv[1]))['convergence']['work']['davidson_stages']['work']
print('stage,rank0_s,max_s,nominal_GFLOP,copied_B,logical_B,gemms,collectives,payload_B,GFLOPs,GBs,max_GFLOPs,max_GBs')
for stage, a in w['rank0']['stages'].items():
    b = w['max'][stage]
    coll = a['collectives'].values()
    print(','.join(map(str, [stage, a['seconds'], b['seconds'],
        a['nominal_flops']/1e9, a['copied_bytes'],
        a['streamed_bytes']+a['gemm_operand_bytes'], a['gemm_calls'],
        sum(c['calls'] for c in coll), sum(c['payload_bytes'] for c in coll),
        a['nominal_gflop_s'], a['logical_gb_s'],
        b['nominal_gflop_s'], b['logical_gb_s']])))
print('rank0 solves:', w['rank0']['solves'])
print('rank0 GEMM histogram:', json.dumps(w['rank0']['gemm_shapes']))
```

Rates are rank-local before taking maxima; never divide a maximum work count
by a different rank's maximum time. Never sum maximum stages into wall.
Collective seconds are subsets of stages; logical input payload is not wire
traffic. Copy destination bytes and read+write traffic are different columns.
GEMM shapes count adapter calls, not BLAS microtiles. H/EVD internal arithmetic,
backend packing/scans/allocations and diagnostic storage are excluded from
nominal counts. Hypot/sqrt are separate operations, not assigned fake FLOPs.
Full factor panels are nominal; rejected panels are flagged and overestimate
partial scalar work. Do not interpret logical GB/s as measured DRAM bandwidth.

Use existing QE-matched clock totals for the original remainder comparison;
new exclusive stages help decompose it but have different eigen/H boundaries.
Full `work` includes verification; do not compare it to a QE unverified interval
without matching solve scope. No new QE run is requested. Any QE table uses
MPI8 against MPI8 only, matched conditions and genuine saved provenance.
The screen can locate the cost and justify a later switch decision, but does
not establish a precise speedup or prove the full 1.5x gap removed.

## Mac evidence and handoff ownership

`cargo test -p dft-solvers davrem -- --nocapture` passed the two new solver-fixture tests; the full suite also runs the new collective-forwarding test.
Synthetic fixed-H probes use dimensions/bands/capacities 12/3/6, 17/5/15,
7/4/7; cyclic 1/3/8-rank partitions, Root and Replicated EVD, each independent
switch and their combination. All eigenvalue/vector/HX bits, flags, Ritz
histories, H counts, repair/drop decisions and logical GEMMs match baseline.
Diagnostic off publishes no work; full diagnostics preserve numerical results.
The small fixture forces repair/drop and includes a zero-row rank at P=8.
These are synthetic correctness results, not material or QE timing evidence.

| fixture (dim/bands) | ranks | H columns | corrections | restarts | B restart copy destinations saved, rank0 B | C reset stores saved, rank0 B |
| --- | --- | --- | --- | --- | --- | --- |
| 12/3 | 1 / 3 / 8 | 13 | 4 | 3 | 1728 / 576 / 288 | 1080 |
| 17/5 | 1 / 3 / 8 | 20 | 3 | 1 | 1360 / 480 / 240 | 1000 |
| 7/4 | 1 / 3 / 8 | 7 | 2 | 1 | 448 / 192 / 64 | 640 |

The one-correction analytic probe checks 36 hypot calls, 468 nominal scalar
FLOPs and 1176 logical bytes in precondition, three sqrt calls and one 24-byte
norm reduction; error exits clear the diagnostic scope. Collective stage sums
and GEMM FLOPs reconcile with the existing independent solve record.
Raw local command logs are in `target/tmp/davrem-*.log` (uncommitted scratch).
No Mac timing is used to accept a performance change.

Canonical ledger/optimization-log, README generation, session registry and
dashboard publication are outside OWNED_FILES; master integrates this request's
attempt record and verified evidence after review. No recurring timer or
external dashboard upload was installed by this implementer. Changes remain
uncommitted as requested. Both switches remain default off until dev08 review.

### Verified Mac commands

| Command | Exit | Result |
| --- | --- | --- |
| `cargo check --workspace --tests` | 0 | Final tree checked; the first attempt exited 101 because CLI has no production dft-solvers dependency. Fixed using the already exported ApplyH diagnostic bridge; no manifest edit. |
| `cargo test -p dft-solvers` | 0 | 171 passed, 25 ignored across 18 test/doc binaries; includes three new solver/work tests. |
| `cargo test -p dft-engine --test refscf` | 0 | 21 passed, 1 ignored. |
| `BAANDI_REF_DAV_SWAP_BUFFERS=on BAANDI_REF_DAV_LEAN_RESET=on cargo test -p dft-engine --test refscf` | 0 | 21 passed, 1 ignored with both candidates on. |
| `cargo test -p dft-numerics switches` | 0 | 14 passed; registry count 147. |
| `cargo test -p dft-cli --lib davidson_work_json_tests` | 0 | Unequal per-rank shape lists, zero-row-like empty shapes, rank0/max values, null rates, nested collective fields, exact integers above 2^53 and serial reporting. |
| `cargo fmt --check` | 0 | Passed. |
| `cargo clippy -p dft-solvers -p dft-numerics -p dft-cli --all-targets -- -D warnings` | 0 | Passed. |

The additional collective-forwarding test exercises all five collective kinds
at three simulated ranks: exactly one call each, payload bytes
`[0, 8, 16, 8, 24]` for barrier/sum/broadcast/gather/alltoall. Results are
checked independently of the counters, and stage totals retain the exact calls.
No ignored release/material campaign or `scripts/check.sh full` was run on Mac,
consistent with the ticket's restricted Mac check list. No local performance
acceptance is inferred from test duration. Existing stored README/reference
artifacts were not modified.
