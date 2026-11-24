# P3-NL: nonlocal projection and lazy final orbital delivery

Implementation and requested Mac gates complete; numerical fa055919
shape/GFLOP/s accounting is blocked on the absent raw profile. Both switches
default **off**. No dev08 execution or
performance acceptance is claimed. Derived from the matrix-product/coupling
structure in `docs/research/qe-study/qe-hpsi-and-fft.md` and the engine's own
adapter audit; no QE source was accessed. Source base: `d5937bd33a209e8a1c978d17d2303bb56d85413e`.
Changes are left uncommitted for review.

## What the existing path actually does

Let r be this rank's owned ORC row count, p the total projector count across
all atoms/species, b the width of this H call, C the operator capacity, and
s=max(p,C). These are runtime shapes, not constants inferred from atom count.
P is already column-major r by p with contiguous owned rows, built once at
operator construction. There is no projector conversion/repacking in the H
loop at the Rust layer. BLAS may pack its operands internally.

| Stage | Logical GEMM (m,n,k), column-major operands | Default route on eligible MKL | Work per H call |
| --- | --- | --- | --- |
| Projection P^T X | (p,b,r); P lda=r, X ldb=r, output p by b | `other/mkl`, wide entry, ceil(r/256) reduction-panel BLAS calls | 2 r p b real FLOPs; combine partial p by b outputs |
| Rank reduction | p by b | one deterministic `tree_sum_f64` | 8 p b bytes of payload; no GEMM FLOPs |
| Coupling W=D Q | one row-major c_a by c_a D block per atom | scalar loop, atom/band/j/k order | 2 b sum_a(c_a^2) nominal FLOPs, including zero entries |
| Back P W | (r,b,p); P lda=r, W ldb=p, image r by b | `other/mkl_wide` when r*b>s*s, otherwise `other/mkl`; ceil(r/256) output-row BLAS calls | 2 r p b real FLOPs |
| Image addition | Y += image | contiguous scalar pass | r b additions; 24 r b nominal bytes |

Both projector products already cover **all** projectors in one logical
adapter call. M2_DGEMM off uses `pure_direct`; pure builds use `pure_build`.
The back product can hit `pure_panel_bound` if min(256,r)*b>s*s. Confirm actual
routes in `gemm_dispatch`, including `rows_as_if_mkl` for Mac pure builds.
Zero-owned-row ranks skip the products but participate in the rank reduction.

The existing P3-HAPPLY split is `projection`, `projection_reduce`, `coupling`,
and `back_and_add` (use the actual emitted stage names); calbec encloses the
first two and add_vuspsi encloses the last two. Back-and-add includes the
finite-output check and image addition. Do not call its rate pure BLAS speed.
For each rank, accumulate F=sum_calls(2*r*p*b). The achieved projection rate
is F/(1e9*projection_s); back-and-add effective rate is
F/(1e9*back_and_add_s). A separate backend-only rate needs the attributed
DGEMM clocks. Report total FLOPs/total seconds, not the mean of per-call rates.
Coupling rate uses sum_calls(2*b*sum(c_a^2))/coupling_s/1e9. Report reduction
seconds/call and payload bytes, not GFLOP/s. calbec per-call time is
(projection_s+projection_reduce_s)/H_calls; add_vuspsi is
(coupling_s+back_and_add_s)/H_calls. Keep rank 0 and each rank's extrema separate;
never divide rank-0 work by a componentwise max clock from another rank.

### Available numbers and missing evidence

The ticket supplies one **relayed** fa055919 Si128 MPI8 rank-0 run, verify off:

| Inclusive stage | BAANDI s/H call | QE rmm-davidson s/H call | BAANDI approximate total (83 calls) |
| --- | ---: | ---: | ---: |
| calbec | 0.0208 | 0.0172 | 1.7264 s |
| add_vuspsi | 0.0245 | 0.0185 | 2.0335 s |
| Output | 0.83 s total | outside QE electrons boundary | not a per-H quantity |

Totals above multiply rounded ticket rates, not recovered raw clocks. The
0.6 s orbital interval and 0.3 s nested collectives overlap: do not add them.
The supplied data do **not** include r,p,the 83 b values, exclusive split,
route rows or the raw profile. No P3k/fa055919 raw profile or P3k ledger row is
present in this checkout. Thus numerical achieved GFLOP/s and numerical GEMM
shapes for this run are **unavailable**, not zero or estimates. The older
LANE-NONLOCAL bad02ee Si64 relay (34.15 projection / 27.90 back GFLOP/s) is a
different path/run and must not fill this table. Recover the fa055919 records
or fill the requested table from the new runs below.

## Candidate N: BAANDI_NONLOCAL_SINGLE_GEMM=on

Only the local P^T X product changes: use the existing
`multiply_slices_single_into` instead of the 256-term reduction-panel entry.
The entire p by b output fits the existing s by s scratch by construction.
There is one sequential MKL call per projection, retaining the same P and X
storage. No scratch allocation, logging, environment reads, or I/O is added
inside H. Resolve the switch once at operator construction; M2_DGEMM off
still overrides it. On pure builds both entries use the same panel arithmetic.

Removed work: ceil(r/256)-1 backend dispatches per projection, repeated
per-panel operand packing as implemented by BLAS, and the Rust summation of
partial projection tiles. Backend packing savings are a hypothesis, not a
measured byte count. This is a **rounding change on MKL**. D coupling, the
single deterministic p*b rank reduction, P W and image addition are unchanged.
The independent nonlocal energy audit uses the same projection helper and
therefore takes the selected entry too; include its cost/work separately.

Back cannot use this single API without increasing scratch: its r*b output
usually exceeds s*s. We do not inflate every GEMM tile to accommodate it.
There is no adjacent same-dependency collective that can absorb P^T X's
reduction: the back product needs the completed global coefficients. No
asynchronous collective or reduction-order change is introduced.

P3-HAPPLY's solver-side two-entry counter predates this switch. CLI reports
one `nonlocal_single_entry_calls` plus one wide entry per eligible H when N
is on; off retains its old JSON fields. `gemm_dispatch` is authoritative for
effective backend routes. Its `other` site also includes non-H consumers;
do not attribute its process total wholly to H.

This does not repeat A43/A58b traversal, column grouping or PACKED_REDUCE.
LANE-NONLOCAL's closing assessment predates the reference path and the
successful P3 single-GEMM Gram change; the new premise is the supplied
reference-path per-call gap and an existing bounded single-call adapter.
It is a screen candidate, not an asserted speedup.

## Candidate O: BAANDI_REFSCF_LAZY_OUTPUT=on

`driver::run` can now retain final owned ORC rows in the existing LazyOrbitals
handle when the CLI says orbitals are unnecessary and M2_LAZY_RESULT is on.
The real block moves into the handle; the two public complex arrays are
0 by bands placeholders. `ScfResult::materialize_orbitals` remains collective
and restores both arrays on demand. Requested orbital output remains eager.
The exact active SharedComm is passed into the internal driver; the generic
borrowed `run_with` API stays eager, so no unrelated registered communicator
can be attached accidentally. Serial runs also stay eager.

This removes one full N*b tree sum, the two full real delivery buffers, complex
unpacking and both public complex orbital copies from ordinary MPI CLI output.
The owned block is retained until the result is dropped. Density/charge,
eigenvalues, occupations, energies, SCF/verification and solver work are not
changed. The conservative pre-existing peak envelope is retained.

Density already has exactly one final gather and no second CLI gather.
`density.bin`, density metadata and optional cube/XSF writers need the full
C-order grid. The established CLI writes complete results and density under
each rank's directory, not just root. Writing disjoint owned segments or
root-only density would change that file contract and requires writer/MPI
adapter scope not owned here. Therefore the existing one density gather and
all rank files are deliberately preserved. No unsupported owned-layout writer
or silent rank-file removal is claimed.

Scientific results.json content and density bytes are preserved. Full JSON
bytes cannot match across runs: timings, hashes, run IDs, switch provenance
and collective telemetry change. The new output path consequently stays
behind its own default-off switch. No schema keys are removed. O never turns
verification off and its wall cost remains inside delivered wall.

## dev08 request: ONE profiled run per arm per cell

Master/baandi-dev08 only; implementer must not log in or install MPI/MKL.
Use frozen si64, si128, elongated and rskew MPI8 inputs from the current
campaign, exact UPFs/grids/bands/starts and socket-0 eight-core binding. Record
input/UPF SHA-256, source commit plus patch hash, binary SHA-256, compiler,
MPI/MKL versions, actual rank CPU bindings, load, commands and exit codes.
Never infer geometry or dimensions from a case name. One thread per rank.

| Arm | NONLOCAL_SINGLE_GEMM | REFSCF_LAZY_OUTPUT | Purpose |
| --- | --- | --- | --- |
| A | off | off | Same reviewed binary control |
| N | on | off | Projection dispatch only |
| O | off | on | Delivery only |
| NO | on | on | Stack/interactions |

Exactly **16 profiled runs**: four arms by four cells, no warmups/repeats or
precision campaign in this request. Keep other switches identical, including
M2_DGEMM=on, MKL_DGEMM=on and M2_LAZY_RESULT=on. Use the campaign's verify-off
inputs to reproduce the motivating boundary; label them **unverified** and
compare with its independently produced strict oracle. If the master instead
uses verification on, hold it on for every arm and split ordinary/verification
work; do not splice those rates into the fa055919 verify-off row. A later
promotion campaign must validate independent fixed-point verification and
strict 1e-10 Ha / standard 5e-9 Ha classes; no tolerance change.

Build separately from current main (baseline fa055919 drift is not the candidate
effect):

```sh
CARGO_TARGET_DIR=target/p3-nl cargo build --release -p dft-cli --features mpi,mkl
```

Launch template (resolve frozen paths/binding with the master):

```sh
BAANDI_NONLOCAL_SINGLE_GEMM=off BAANDI_REFSCF_LAZY_OUTPUT=off \
BAANDI_M2_DGEMM=on BAANDI_MKL_DGEMM=on BAANDI_M2_LAZY_RESULT=on \
BAANDI_DAVIDSON_STAGES=on BAANDI_STAGE_TIMERS=on BAANDI_DGEMM_TIMERS=1 \
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/usr/bin/time -p mpirun -np 8 --map-by core --bind-to core --report-bindings \
/absolute/reviewed/dft-cli run /absolute/frozen/input.toml \
--output-dir /absolute/unique/p3-nl/case/arm
```

One case/arm per fast job, 600 s maximum; no truncation of converged runs.
A longer run needs the master's named authorized final job (7200 s maximum).
No builds during profiling. Store complete stdout/stderr and rank directories.

Return a table per case/arm/rank with: r,p,C; exact columns_per_call; H calls
and columns; GEMM shapes/routes/logical FLOPs; projection_s, projection_reduce_s,
coupling_s, back_and_add_s; achieved projection and effective back GFLOP/s;
calbec/add_vuspsi s/call; density_gather, orbital_gather_unpack and nested
collective seconds; engine output_s, CLI output_s and external delivered wall.
Attribute native GEMM timers to H versus startup/audit; include unmatched
records explicitly. Preserve raw per-rank clocks alongside merged max/min.
Check sum(b), pair/tail closure, call counts and route transitions A -> N.
Report O's omitted final tree-sum payload (8*N*bands bytes per rank), density
hash equality and scientific JSON equality against A, and eigenvalues/H work.
Do not sum overlapping output/collective timers. All cells get their own rows.

QE comparisons are MPI8 against MPI8 only, identical conditions and fastest
valid QE configuration. The ticket's rounded QE timers are motivation only;
retrieve real command, version, input and UPF hashes before publishing a QE
reference. QE internal electrons excludes output; compare external delivered
wall separately. These one-run profiles screen mechanisms, not uncertainty or
a precision speed claim. Native rank/thread reproducibility and multi-projector
MKL correctness remain required before promotion.

## Local validation and handoff

Mac pure backend, simulated collectives only. Final command records are in
`target/p3-nl-commands.json`, logs in `target/p3-nl-*-verified.log` (switches:
`target/p3-nl-switches-final.log`). No benchmark timing is used on Mac.

| Command | Exit | Result |
| --- | ---: | --- |
| `cargo check --workspace --tests` | 0 | Workspace and test targets compile. |
| `cargo test -p dft-physics` | 0 | 97 passed, 18 ignored, 0 failed. |
| `cargo test -p dft-engine --test refscf` | 0 | 24 passed, 1 ignored, 0 failed. |
| `cargo test -p dft-engine --lib reference_lazy_delivery -- --nocapture` | 0 | 2/3 simulated ranks; eager, lazy and explicit-orbital-request arms converge and pass independent verification. Materialized orbitals, energy and density bits agree. |
| `cargo test -p dft-cli reference` | 0 | 10 passed, 9 ignored, 0 failed; two-rank output fixture has equal scientific JSON/history/work and byte-identical density.bin/density.json. |
| `cargo test -p dft-numerics switches` | 0 | 14 passed; registry count 156. |
| `cargo fmt --check && cargo clippy -p dft-physics -p dft-engine -p dft-cli --all-targets -- -D warnings` | 0 | Format and warning gates pass. |
| `git diff --check` | 0 | No whitespace errors. |

The nonlocal switch fixture covers cubic, skewed and elongated cells, empty
and mixed 1/3/5-channel atom blocks, 0/1/2/5/17 bands and 1/2/3/8 simulated
ranks. A larger cubic case explicitly crosses the 256-row panel boundary.
There are **280 rank/block comparisons**, maximum H component difference
**0 Ha** (bit-identical pure backend). Native MKL rounding is
not exercised here. The CLI fixture intentionally stops at two maps (tests
exhaustion delivery); the separate engine fixture tests fully verified
convergence and on-demand materialization. Neither is external-code evidence.
Full results.json bytes are not asserted identical; its scientific fields,
complete clock-stripped iteration history, work counts and density files are.

Earlier attempts: workspace check exited 101 for missing test/import names
(fixed); CLI test exited 101 because synthetic Xx lacked an explicit element
(fixture changed to He); clippy exited 101 for duplicate test-fixture module
imports (reused the existing fixture modules). Every affected command was
rerun successfully, with no tolerance changes. A temporary unexpected-cfg
warning in the physics test was fixed by using `dft_numerics::MKL_BACKEND`.

 The task explicitly limits checks to the
listed Mac commands; full `scripts/check.sh`/native MPI/MKL gates remain with
the master. Central ledger/dashboard/session registry and README-generated
results are outside OWNED_FILES. This document is the scoped draft ledger and
handoff record; no fabricated remote numbers are added to those artifacts.
The checkout registry has no P3-NL entry or .delegate/P3-NL/report.json; no
recurring reporting timer is active in this implementer session.
