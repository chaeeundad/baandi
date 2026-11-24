# P4-DENSKEEP: remove slab storage, retain starting HX

For the master session and `baandi-dev08`. Base:
`1a3fc4347ce65e0933a6f03c189b8621ae0cebfa` plus the reviewed uncommitted
P4-DENSKEEP diff. No remote run, default adoption, or QE source access occurred.
The new `BAANDI_REFSCF_KEPT_RECOMPUTE` switch defaults **off**.

## Mechanism and decision to test

The supplied fa4cfad Si128 MPI8 snapshot reports density 2.327 s / 16 calls,
versus QE sum_band 1.78 s. These are task motivation, not a new matched timing
comparison. Ledger KEPTHX2 records the earlier solve saving of about 2.1 s
and increased density time/RSS. KEPTHX2 already removed per-map allocation and
fused copying with accumulation. Repeating those changes cannot remove the
remaining 16*K*P-byte store (K kept complete pairs, P local plane points).

The candidate keeps exact f64 X/HX/old potential and all eligibility and
verification rules, but retains **zero** real-space pairs. Existing
`apply_local_change` reconstructs inverse planes in f64. Density makes no
slab writes or allocation; every eligible hit still replaces B full-H columns
and avoids the kinetic/nonlocal parts. Each hit adds K inverse pair transforms
relative to storing K pairs. For even B with full retention, it also restores
the mirror-row exchange/inverse transpose previously avoided altogether.
Forward transforms are unchanged. This deliberately trades FFT work for
memory and density traffic; **the old 2.1 s solve saving is not promised**.
No guessed machine/material crossover is built into the default.

## Exactly four remote executions

Use one reviewed candidate binary, built once with the campaign's unchanged
`--release -p dft-cli --features mpi,mkl` flags/toolchain/backend. Keep both
arms' build flags identical; do not introduce target-cpu tuning.

| Cell | First arm | Second arm | MPI ranks | Threads/rank | Verify |
| --- | --- | --- | ---: | ---: | --- |
| archived precise-campaign Si128 | H | R | 8 | 1 | off |
| archived precise-campaign Si216 | R | H | 8 | 1 | off |

- H: `BAANDI_REFSCF_KEPT_HX=on BAANDI_REFSCF_KEPT_RECOMPUTE=off`.
- R: `BAANDI_REFSCF_KEPT_HX=on BAANDI_REFSCF_KEPT_RECOMPUTE=on`.

**One run per arm per cell**, sequential, no warmup SCFs, repeats, added
strict runs, or verification-on arms. Use the archived inputs/manifests under
`/data/sig-dft/precise-2026-11-26/`, identified by the Si128/Si216 MPI8 verify-off
rows of `table.json`/`TABLE.md`. Recover the full input, UPF and strict-oracle
hashes/commands from that archive; this checkout does not contain them.
Keep cell/coordinates, cutoff, FFT grid, bands, starting orbitals/density,
threshold, beta/history, backend and all other resolved switches identical.
Explicitly keep ladder and RMM off. Use socket-0 physical cores 0–7, the same
archived rank binding, engine/OMP/MKL threads all one, no concurrent builds.
Each job has the 600 s fast-job limit. Preserve failures/timeouts; do not rerun
or adjust accuracy. No QE execution or QE speed ratio is requested.

Keep **lightweight diagnostics** in both arms (`BAANDI_DAVIDSON_STAGES=off`,
`[profile] davidson_stages=false`): the reference driver retains coarse
`Record::density_s` and `solve_s` in this mode. Do not add heavy profiling or
Ritz traces to these four runs. `Record::density_s` includes field redistribution;
retain that boundary identically in both arms. Record rank-local values and
per-map records rather than adding maxima from unrelated ranks. Capture
external wall including setup/output/MPI exit; retain engine wall separately.
Capture peak RSS **per rank** using the campaign's rank wrapper; the maximum
RSS of the mpirun launcher alone is not the rank RSS. Record units explicitly.

Invocation pattern, using the exact archived launcher/binding/options and
fresh output directories, with `[scf] verify="off"` in both input copies:

```sh
BAANDI_REFSCF_KEPT_HX=on BAANDI_REFSCF_KEPT_RECOMPUTE=off \
BAANDI_DAVIDSON_STAGES=off OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /usr/bin/time -p mpirun -np 8 <archived-binding-and-rank-RSS-wrapper> \
  <candidate> run <H-input> --output-dir <H-output>
BAANDI_REFSCF_KEPT_HX=on BAANDI_REFSCF_KEPT_RECOMPUTE=on \
BAANDI_DAVIDSON_STAGES=off OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /usr/bin/time -p mpirun -np 8 <archived-binding-and-rank-RSS-wrapper> \
  <candidate> run <R-input> --output-dir <R-output>
```

Archive commit plus patch/source hashes, exact build command, binary SHA-256,
MPI/backend versions, input and UPF SHA-256, environment/resolved-switch
manifest, actual binding, load, exit status, stdout/stderr/results and the
existing matching strict reference with provenance. Missing reference/input
identity is a qualification limit, never permission to substitute numbers.

Return this table (one observation per cell/arm, **not a median**):

| Cell/arm | Density s | Solve s | External wall s | Peak RSS/rank bytes | Retained slabs/rank bytes | E Ha | E - E_strict Ha | E_R - E_H Ha |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Si128 H | | | | | | | | |
| Si128 R | | | | | 0 expected | | | |
| Si216 H | | | | | | | | |
| Si216 R | | | | | 0 expected | | | |

Also return density hashes, convergence/exit reason, map counts, Davidson
iterations, full-H columns, kept-HX hits/replaced columns, per-hit and density
inverse/forward counts where existing diagnostics expose them. Separate the
first density map from subsequent ones to distinguish first touch. Expected:
identical density FFT work, zero slab payload in R, same full-H replacements,
extra inverse transforms only in the local-change path. Map/count changes
must be reported, not hidden by a per-call average.

Precision gate: `abs(E_R-E_H) <= 1e-10 Ha`, with density/count identity expected
on the pure Mac fixtures; MKL identity remains unmeasured. Report both signed
absolute strict-reference errors too. The archived verify-off campaign used
the standard stopping threshold and can already exceed 1e-10 Ha versus its
strict oracle: relative equality does not certify absolute strict accuracy.
Do not relabel a standard result strict or change its threshold to pass this
screen. Default promotion requires the master to evaluate density + solve +
wall, observed RSS, delivered accuracy and the broader validation gates.
A memory win alone does not establish that recomputation is faster; a single
run cannot resolve a small gain against the recorded roughly 1 s Si128 noise.

## Mac measurements and task-local attempt record

macOS 27.0.1 arm64, 16,384-byte pages; pure Rust FFT, simulated collectives,
not real MPI. The diagnostic tests are ignored in ordinary checks and run in
release. This task's narrow OWNED_FILES excludes the shared ledger/log,
README generator/evidence, dashboard JSON and delegate registry. This section
is the attempt record for master-session reconciliation. The local registry
has no P4-DENSKEEP entry or `.delegate/P4-DENSKEEP/report.json`; no registry
status or dashboard upload is fabricated.

### Common-input density, real operator

`cargo test --release -p dft-engine --test refscf denskeep_density_probe -- --ignored --nocapture`
(or run the compiled test executable with those test arguments). Two simulated
ranks, skew one-atom synthetic nonlocal fixture, 16 bands, grid 64x64x64,
16 density calls on exactly the same X. One new cache per arm, H then R.
Counters: **128 inverse pairs, 0 forward pairs per rank in each arm**;
all output density bits match. The 16,777,216-byte/rank slab is initialized
once and written per map in H; R has zero bytes. Logical eliminated stores
are 268,435,456 bytes/rank over 16 calls, excluding the initial zeroing.

| Rank/arm | First density s | Median of next 15 calls s | Stored bytes/map |
| --- | ---: | ---: | ---: |
| 0 H | 0.003661667 | 0.002714458 | 16,777,216 |
| 0 R | 0.001965667 | 0.002528208 | 0 |
| 1 H | 0.003670792 | 0.002696708 | 16,777,216 |
| 1 R | 0.002561917 | 0.002473667 | 0 |

These sequential Mac observations include FFT/collective/snapshot overhead;
they do not establish a dev08 speedup or isolate hardware cache behavior.

### Large memory replay: first touch and store traversal

`denskeep_store_probe` replays the production plane traversal with synthetic
f64 values, no FFT: 128 pairs, 4096 complex points per plane, contiguous
16-byte stores. Requested 185/500 MiB sizes round up to whole planes,
**192/504 MiB** actual payloads. The plain arm uses the no-store two-pass
accumulation; fused uses current store + two accumulations per point; split
uses a separate contiguous copy followed by two accumulation passes.
Density bits agree in all arms. No production split-path change is proposed.

Build once, then run **sequentially with no overlapping builds**:

```sh
cargo test --release -p dft-engine --test refscf denskeep_store_probe --no-run
DENSKEEP_PROBE_ARM=fused DENSKEEP_PROBE_MIB=185 DENSKEEP_PROBE_MAPS=16 \
  <printed-refscf-test-executable> --exact denskeep_store_probe --ignored --nocapture
```

Sweep arms plain/fused/split, sizes 185/500, maps 1/16 (12 processes). Collect
process faults and RSS by launching each executable from a fresh Python
process with `subprocess.run(..., env=...)`, then
`resource.getrusage(resource.RUSAGE_CHILDREN)`; record `ru_minflt`, `ru_majflt`,
`ru_maxrss` (bytes on this Mac). `/usr/bin/time -l true` failed with exit 1:
`sysctl kern.clockrate: Operation not permitted`; Python collection succeeded.
These are **whole-process minor faults**, not in-kernel hardware page-fault
samples; controls include the same harness and source/density allocations.

| MiB | Arm | Maps | Init s | First scan s | Median scan s | Minor faults | Peak RSS bytes |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 192 | plain | 1 | 0.000000042 | 0.010120291 | 0.010120291 | 474 | 3,833,856 |
| 192 | fused | 1 | 0.010657083 | 0.012825166 | 0.012825166 | 12,762 | 205,193,216 |
| 192 | split | 1 | 0.010743584 | 0.018222458 | 0.018222458 | 12,761 | 205,193,216 |
| 192 | plain | 16 | 0 | 0.012463458 | 0.011721459 | 474 | 3,850,240 |
| 192 | fused | 16 | 0.010258583 | 0.011323125 | 0.011192750 | 12,760 | 205,193,216 |
| 192 | split | 16 | 0.010093958 | 0.012699958 | 0.011294291 | 12,763 | 205,193,216 |
| 504 | plain | 1 | 0.000000041 | 0.031320667 | 0.031320667 | 551 | 5,111,808 |
| 504 | fused | 1 | 0.024716125 | 0.038139291 | 0.038139291 | 32,807 | 533,610,496 |
| 504 | split | 1 | 0.023759000 | 0.044185750 | 0.044185750 | 32,807 | 533,610,496 |
| 504 | plain | 16 | 0.000000042 | 0.037772417 | 0.028495250 | 554 | 5,144,576 |
| 504 | fused | 16 | 0.024853667 | 0.040813583 | 0.030103750 | 32,811 | 533,643,264 |
| 504 | split | 16 | 0.024793084 | 0.044482750 | 0.035100625 | 32,809 | 533,626,880 |

All processes exit 0; major faults are zero. The fused-minus-plain faults for
one map are **12,288 / 32,256**, exactly payload/16,384. Sixteen-map processes
have virtually the same counts: pages are reused, not faulted anew per map.
Logical stored bytes/map are **201,326,592 / 528,482,304**; 16-map fused medians
give **17.987 / 17.555 GB/s** payload divided by full scan time, **not measured
DRAM bandwidth** (source reads, density updates and write allocation excluded).
The density working sets are 786,432 / 2,064,384 bytes; streaming slabs can
compete for cache, but cache misses/evictions were not measured. The 192 MiB
steady observations do not even show a stable fused penalty over plain.
The 504 MiB observations show extra store/traversal cost; they cannot assign
it uniquely to bandwidth versus cache pollution. No Mac timing is adoption
evidence. Larger remote first-map/steady-map profiles are still needed.

### Alternatives considered

- Local-only storage is already implemented: `apply_local_batch_kept` reads
  every retained owned-plane value. Removing a subregion would require
  reconstruction; buffer-band occupations do not make those values unused.
- Pre-faulting relocates the one-off initialization cost and retains the RSS;
  existing reusable storage plus the fault counts rule out repeated slab
  first touch as the steady regression. No redundant pre-fault switch added.
- Split copy: measured above as a control for traversal/cache interaction;
  no evidence to overturn KEPTHX2's fused implementation, so not shipped.
- f32: could halve payload but changes HX precision and needs an FFT consumer
  for compressed planes (the existing consumer is in unowned dft-solvers).
  Full decompression would recreate transient f64 slab memory. Not implemented
  or claimed accurate; the implemented option is entirely f64.
- Non-temporal stores: require an architecture-specific adapter/measurement;
  unsafe intrinsics are outside this ticket's permitted physics/engine layer.
  No portability or speed claim is made.
- Recompute: retained as the explicit default-off memory option. It eliminates
  slab stores and payload; remote measurement decides whether saved memory
  traffic compensates for added inverse transforms/communication.

### Delivered accuracy and work

`cargo test -p dft-engine --test refscf fast_order_ -- --nocapture --test-threads=1`
passes 9 tests. New product tests cover 1/2/9 simulated ranks, B=3/4,
unoccupied tail/buffer bands and ranks without planes. R HX equals H HX bit
for bit. Maximum relative Frobenius error versus direct H is
**1.20867983040270716e-16**, below the unchanged 1e-13 bound. Density inverse
pairs stay ceil(B/2); per-hit inverse pairs change 1->2 (B=3), 0->2 (B=4).
Forward pairs stay 2; full-H replaced columns stay B; R slabs are zero.

Four full SCFs, two simulated ranks, threshold 1e-12 Ry, verification on:

| Fixture | Delivered E (Ha), both arms | Maps | Ordinary Davidson iterations | Ordinary full-H columns | Hits including verification | Replaced columns | H slabs/rank bytes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| cubic/local | -1.10301826312976980 | 7 | 11 | 30 | 7 | 21 | 4096, 4096 |
| skew/nonlocal | -0.839575255872316961 | 28 | 60 | 125 | 28 | 84 | 4096, 4096 |
| elongated/local | -1.17570725935517029 | 7 | 10 | 31 | 7 | 21 | 2160, 1440 |
| slab/local | -1.07379706580614176 | 16 | 33 | 74 | 16 | 48 | 3600, 2880 |

All energy differences are **0 Ha**, all density bits match, ordinary and
verification Davidson/H-column counts match, independent verification passes,
and R slabs are zero on every rank. These are internal synthetic comparisons
against f64 H, not absolute external strict-oracle or material validation.
The unchanged existing delivery fingerprint is `43cc717a519f09dd` both before
and after the production edit (28 maps, 28 hits, 84 replaced columns).

An initial 8x8x4 Bohr slab fixture failed before density with `rank deficient
atomic starting orbitals`: at cutoff 1 Ha its thin direction excludes p_z.
Two test invocations exited 101 (the first attempted fixture edit did not
match the formatted source, so the second still used 8x8x4). The final
10x10x5 Bohr fixture resolves that component and passes. No tolerance,
production start rule or accuracy gate was changed.

### Verification commands

All commands run locally with the pinned toolchain/offline vendored crates.

| Command | Exit | Result |
| --- | ---: | --- |
| `cargo test -p dft-engine --test refscf fast_order_ -- --nocapture` before production edits | 0 | 4 passed; existing fingerprint recorded |
| same fast-order test, initial invalid slab fixture (two invocations) | 101 | pre-density start failure documented above |
| `cargo test -p dft-engine --test refscf fast_order_ -- --nocapture --test-threads=1` after fixture correction | 0 | 9 passed |
| `cargo check --workspace --tests` | 0 | workspace/test type check |
| `cargo test -p dft-engine --test refscf` | 0 | 58 passed, 7 ignored |
| `cargo test -p dft-physics` | 0 | 97 passed, 18 ignored across targets |
| `cargo fmt --check` | 0 | clean |
| `cargo clippy -p dft-engine -p dft-physics -p dft-numerics --all-targets -- -D warnings` | 0 | clean |
| release memory probe, 12 direct test-binary executions | 0 each | counts/times/RSS above |
| release real-operator density probe | 0 | two-rank counts and bits above |
| `python3 -I scripts/tests/test_readme_figures.py` | 0 | 9 passed, generated evidence remains current |

No README evidence/figure change is made: no new external comparison or
promoted performance result was produced, and the shared files are outside
this task's ownership. Temporary raw logs are under `/private/tmp/p4-*.log`;
the tables above preserve the evidence needed for review in the owned request.

The real-operator density diagnostic executable SHA-256 was
`e7478c7abc6546326ca604ce63084be44ed2d5e529806281494a154f1d8b8d76`
(the final slab-fixture-only correction was made after that release build).
The earlier memory-only replay was compiled before the production patch; its
binary hash was not captured before rebuilding. Its test body and exact
parameters remain in the appended test for reproduction; these local timings
are diagnostics, not an acceptance benchmark.

Final standard gate:
`BAANDI_TEST_TMP=/private/tmp/p4-denskeep-check CHECK_SCCACHE=0 scripts/check.sh fast dft-engine dft-physics dft-numerics dft-cli`
exited **0**, with **PASS build/test/fmt/clippy**, 63 binaries,
**690 passed, 0 failed, 172 ignored**, 262 s wall (cold worktree).
The sandbox denied the auxiliary `ps` load report and `nice` priority change;
these warnings did not prevent test execution or alter the successful gates.
`git diff --check` also exited 0. All changes remain uncommitted and inside
the seven owned files listed by `git status --short`.
