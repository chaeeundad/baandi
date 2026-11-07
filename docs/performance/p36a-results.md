# P36A: persistent real Gamma solve state and paired Ritz rotation

**BLOCKED:** implementation and identity checks pass, but the required zero-allocation
Davidson loop still fails inside the existing EVD backend. No performance acceptance.

Baseline: `4bce6acbe563fb7aec04dd3db876765b05b0110f`, branch
`task/P36A-REAL-PERSISTENT-STATE`. Candidate is uncommitted. Raw binaries,
inputs, snapshots, result directories, commands and timings are under
`P36A-evidence/`, outside the single `target/`. No artifacts JSON, thread-budget
selection, symmetry mapping, defaults, complex kernels or QE source were edited.

## Motivation and contract

This implements B21 and only the V/HV part of B22. The supplied
[thread-width study](thread-width-analysis.md), sections 3 and 5, measured
per-solve real Hamiltonian reconstruction at 4.9% of Si8 T=8 wall and 4.4% of
Si16 wall (Mac observations). The [dev08 relay](../../tests/reference/perf/dev08_d0_profile_relayed.md)
reports 6.98% memset and 5.65% memcpy in the post-P32A T=1 profile. Those are
relayed sample shares, not new timings or a predicted saving. A42/A45 support
retaining the dispatch guard and pairing independent rotations; A47/THR-3 rule
out another spin-policy change. A50/THR-4 concern P36B, which is not changed here.

A `KState` now lazily owns its real operator and Davidson scratch. The real
operator retains all pair FFT plans, local workspace, real projectors, nonlocal
image and nonlocal GEMM storage. Davidson retains basis/H images, projection,
GEMM/EVD scratch, corrections, residual blocks, reduced matrix and selection
storage. Each solve resets the reduced prefix width, previous Ritz values and
selection state, reads the current tolerance/policy and constructs fresh random
state from the requested seed. Output arrays/history remain independently owned
by the returned solution and are allocated before the Davidson step loop.
Packing/unpacking at the complex boundary still occurs per solve.

A Hamiltonian has a retained `Arc<()>` identity token. Cell, basis, grid, kinetic
and nonlocal data are immutable through its public API; constructing a replacement
creates a distinct token, even if it has identical numerical content. Keeping the
token alive prevents allocator address recycling from masquerading as identity.
`RealHamiltonian::update_from_complex` verifies the token, exact Gamma basis and
potential shape, then copies only local samples. It returns false without mutation
for a replacement. All kinetic/projector changes therefore rebuild the snapshot.
This avoids a projector scan/hash or projector copy on every iteration.

The engine also keys scratch on band count, effective capacity, requested thread
budget, pool identity, Davidson iteration/subspace/capacity settings, seed,
standard-profile threshold and real-path option. Ineligible dispatch drops the
cache. Per-iteration tolerance and certification policy are supplied afresh;
ratcheting them must not rebuild the same storage. Neither old eigenvectors nor
convergence decisions are implicit warm starts: only the caller's initial matrix
is used. Solver workspace shape/executor checks also protect direct library use.

Paired rotation calls the existing `rotate` tree once per operand under one
`install` and one outer join. Panel boundaries, 4096-MAC bounded backend call
shapes, reduction order and the existing parallel guard are unchanged. Si8's
2969 rows produce 12 panels per operand, hence 24 independent pieces. Other
basis sizes retain `2*ceil(rows/256)` pieces; no fixed 24-way repartition changes
summation. Small products remain serial.

Additive plumbing outside the explicitly named implementation files:
`crates/dft-physics/src/hamiltonian.rs` adds the immutable identity token and its
read-only accessor. SCF changes only add/cache/pass state and include a unit-test
module; no outer/inner budget lines were changed. No new unsafe was introduced.

## Validation protocol

The snapshot harness extends the existing P33A harness before production edits:
four synthetic R5 cells (cubic, skew, elongated/local-only, slab), Si primitive,
Si8 and Si16; strict and standard; T=1,2,4,8,12. The six historical R5 fixtures
already include primitive Si and Si8, so this is seven distinct systems, 70 files.
Snapshots contain final energy/density/eigenvalue bits, per-iteration Davidson/H counts
and complete solver results, excluding wall clocks. The retained environment
variable is `P33A_SNAPSHOTS`; this name preserves compatibility with old commands.
Si16 uses the documented 1x1x2 conventional-cell orientation.

Tests compare reuse with forced reconstruction on successive changing-density
SCF maps, compare freshly reconstructed and reused fixed operators after option
and band-count changes, and check Hamiltonian potential refresh and rejection of
replacement cells/bases/projectors. An engine unit test checks actual retained
mapping storage and invalidation for option, bands, pool, capacity setting and
operator replacement. The paired-GEMM oracle compares every output bit with two
original separate multiplications at panel boundaries and at Si8 dimensions,
including NaN-poisoned output, zero sizes and all five budgets.

The P33A zero-allocation test remains unchanged. Its existing EVD dependency
limitation is documented in [P33A](real-gamma-r9-workspace-results.md); no backend
replacement or tolerance relaxation is part of this ticket.

## Mac measurement limits

Each initial timing batch has one warmup and five measured runs, strict/standard Si8/Si16,
T=1,2,4,8, recording SCF/Davidson/H counts, energy and density hash. The
workspace-only binary was saved before adding paired rotation. No task-owned
build overlaps a timing batch. This is a shared Mac: high load and large outliers were observed, and load
averages were recorded during baseline and per condition in
later batches. The initial batches are sequential, not quiet-node paired acceptance runs;
small differences and the individual rotation benefit cannot be established
from these batches. A further alternating three-binary experiment is reported below.
All samples, including regressions/outliers, are retained.

## dev08 request

After review, compare starting HEAD, workspace-only and final candidate on one
quiet socket, strict/standard Si8 and Si16, T=1,2,4,8,12, pure and MKL, with fixed
accuracy and thread pinning. Alternate at least five runs per condition (30 pairs
for a 1–3% claim), preserve binary/input/UPF hashes and scientific outputs, and
report SCF/Davidson/H counts, density identity, medians and spread. Profile real
Hamiltonian construction, memset/memcpy and paired rotation with at least 30
samples per reported site. No remote command or new QE reference run occurred
here; dev08 decides performance acceptance. The backend allocation gate needs
its own dependency work if it still fails.

## Gate results

**Status: BLOCKED on the unchanged zero-allocation requirement.** Implementation
and numerical identity are reviewable, but this is not completion or performance
acceptance.

- Baseline, workspace-only and final snapshots: **70/70** byte-identical for
  each comparison. Energy difference **0 Ha**; density, all eigenvalues, solver
  results and SCF/Davidson/H counts unchanged.
- All **14 fixture/profile groups** also match bytewise across T=1/2/4/8/12.
  The **16 Si8/Si16 measurement conditions** retain complete scientific SCF
  histories across baseline/workspace/final after recursively excluding timing
  objects (including certificate timings). Per-iteration density/orbital arrays
  are not separately saved by the CLI; the snapshot array comparison is final-state.
- Complex/option-off benchmark: **25/25** comparisons pass `--require-identical`,
  energy difference **0 Ha**, identical density hashes and identical SCF,
  Davidson/H and certificate work counts. Each benchmark command exits 0.
- Complete SCF reuse versus forced rebuild: **8 combinations** (two synthetic
  cells, strict/standard, T=1/4), all identical. A test-only field clears caches
  at every map, including certificate maps; no production flag exists. Separate
  changing-density maps exercise mid-run solver setting changes.
- Real Davidson fixed-operator suite: **9 passed**, including release-only cases.
  Physics suite: **3 passed, 2 ignored**. Paired rotation: **1 passed** across all tested
  shapes/budgets. Final engine unit suite: **24 passed, 3 ignored**; the full
  reuse/rebuild test was also run explicitly and passed.
- `scripts/check.sh`: **exit 0**, build/test/fmt/clippy all PASS.
- README generator and nine freshness/specification tests: **exit 0**, no generated
  README/figure changes. No new accepted reference claim was added.
- The first new engine unit-test compile failed (exit 101) because of a test-only
  import and omitted cache argument; both were corrected and the full engine
  unit suite rerun successfully. The failed log is retained.

### Allocation requirement: still fails after reuse

The existing P33A assertion is unmodified and returns **exit 101**. For example,
T=1 counts are **91 allocations** for 2969 rows/20 bands and **1535** for
11877 rows/68 bands. A new test first completes a solve, retains its workspace,
then counts from the next initial-H return through completion:

| Shape | T | Allocations after persistent initial H | Davidson steps |
|---|---:|---:|---:|
| 2969 x 20 | 1 | 70 | 21 |
| 2969 x 20 | 8 | 72 | 21 |
| 11877 x 68 | 1 | 1499 | 21 |
| 11877 x 68 | 8 | 1502 | 21 |

Its exact zero assertion also returns **exit 101**. Independent localization
reproduces **0/32/88/121 allocations** for EVD sizes 20/40/88/160 with prebuilt
scratch, while direct rotation and Gram GEMMs each count **0**. This is the
existing faer EVD/general-GEMM packing boundary described by P33A; counts also
include runtime/backend initialization and are not all assigned to one source
line. Persisting engine/solver ownership cannot supply that opaque library
storage. Fixing it requires a separately reviewed backend/dependency change
preserving arithmetic, outside this ticket's owned files and library policy.
No allocator replacement, unsafe, relaxed assertion or changed EVD algorithm
was introduced. Full-SCF allocation freedom is not claimed either.

## Reproduction

Prefix cargo with `/opt/homebrew/opt/rustup/bin` in PATH. External UPF is
`../pseudo/Si.pz-vbc.UPF`; SHA-256 is recorded in `P36A-evidence/environment.json`.

```sh
DFT_SI_UPF=/external/Si.pz-vbc.UPF P33A_SNAPSHOTS=/absolute/output cargo test --release -p dft-engine --test gamma_workspace_snapshot -- --ignored --nocapture
cargo test --release -p dft-engine --lib real_state_full_scf -- --ignored --nocapture
cargo test --release -p dft-solvers --test real_workspace no_allocations_after_initial_h -- --ignored --nocapture --test-threads=1
cargo test --release -p dft-solvers --test real_workspace no_allocations_after_persistent_initial_h -- --ignored --nocapture --test-threads=1
cargo test --release -p dft-solvers --test real_workspace localize_allocations -- --ignored --nocapture --test-threads=1
python3 -I P36A-evidence/measure.py before
python3 -I P36A-evidence/measure.py workspace
python3 -I P36A-evidence/measure.py final
scripts/check.sh
```

`commands.json` and `followup-commands.json` retain command/exit-code history;
`verify.py` and `verify_followup.py` reproduce the validation sequences. The
baseline binary was copied before production edits, and the workspace-only
binary before paired rotation. Raw evidence must be archived separately before
worktree cleanup; no git commit was created, and artifacts JSON was left intact.
Peak retained RSS, native MKL behavior and dev08 performance were not measured.

## Final Mac timing observations

The repeated experiment alternates all three saved binaries, rotating their order
each round, with one warmup and five measured samples per binary per condition.
No task-owned build overlaps either timing experiment. These five triples reduce
time-order bias but do not establish a 1–3% effect; this is a shared host, not dev08.
All 16 conditions keep identical energy/density/work counts in every run.
Full raw samples and per-sample load are in `paired-timings.json`.

| Case / profile / T | Before (s) | Workspace only (s) | Final (s) | SCF (all) | Davidson (all) | H vectors (all) | Same trajectory |
|---|---:|---:|---:|---:|---:|---:|---|
| si8-strict-T1 | 1.194820 | 1.166913 | 1.160916 | 15 | 253 | 2001 | Yes |
| si8-strict-T2 | 0.997001 | 0.955405 | 0.949825 | 15 | 253 | 2001 | Yes |
| si8-strict-T4 | 0.784887 | 0.719752 | 0.727408 | 15 | 253 | 2001 | Yes |
| si8-strict-T8 | 0.752295 | 0.711260 | 0.699368 | 15 | 253 | 2001 | Yes |
| si8-standard-T1 | 0.507100 | 0.472371 | 0.469970 | 17 | 89 | 877 | Yes |
| si8-standard-T2 | 0.390327 | 0.353438 | 0.354803 | 17 | 89 | 877 | Yes |
| si8-standard-T4 | 0.348365 | 0.311618 | 0.312287 | 17 | 89 | 877 | Yes |
| si8-standard-T8 | 0.337523 | 0.300506 | 0.295843 | 17 | 89 | 877 | Yes |
| si16-strict-T1 | 5.391738 | 5.353987 | 5.293381 | 14 | 264 | 3254 | Yes |
| si16-strict-T2 | 4.026223 | 3.895578 | 3.905112 | 14 | 264 | 3254 | Yes |
| si16-strict-T4 | 2.889066 | 2.777441 | 2.750283 | 14 | 264 | 3254 | Yes |
| si16-strict-T8 | 2.939219 | 2.806166 | 2.779107 | 14 | 264 | 3254 | Yes |
| si16-standard-T1 | 1.967928 | 1.866178 | 1.847443 | 16 | 88 | 1521 | Yes |
| si16-standard-T2 | 1.454688 | 1.353469 | 1.324543 | 16 | 88 | 1521 | Yes |
| si16-standard-T4 | 1.547732 | 1.435421 | 1.402975 | 16 | 88 | 1521 | Yes |
| si16-standard-T8 | 1.255008 | 1.134854 | 1.135642 | 16 | 88 | 1521 | Yes |

Here “same trajectory” means unchanged compared scientific histories/work and final
arrays for the real path, not equality with the complex representation. Work columns
follow P33A: sums over SCF history; certificate work remains separately recorded in
raw results and included in the independent snapshot comparisons.

### Retained initial sequential batches

These original five-sample batches are retained, not silently discarded when the
shared-host timing became unstable. Parentheses give max-minus-min divided by the
median, in percent. Large spreads and shifting load preclude interpreting their
before/after differences as an isolated speedup.

| Case / profile / T | Before s (spread %) | Workspace s (spread %) | Final s (spread %) |
|---|---:|---:|---:|
| si8-strict-T1 | 1.238183 (9.2%) | 1.159328 (2.8%) | 1.625548 (13.9%) |
| si8-strict-T2 | 1.217795 (36.4%) | 0.969096 (1.7%) | 1.852066 (69.5%) |
| si8-strict-T4 | 0.896640 (19.8%) | 0.807497 (2.5%) | 2.726584 (26.1%) |
| si8-strict-T8 | 0.915583 (48.4%) | 0.794968 (172.6%) | 0.778197 (120.7%) |
| si8-standard-T1 | 0.557126 (15.3%) | 0.468443 (0.9%) | 0.501478 (1.5%) |
| si8-standard-T2 | 0.391447 (3.0%) | 0.352237 (0.8%) | 0.406688 (2.1%) |
| si8-standard-T4 | 0.349183 (2.0%) | 0.313985 (1.7%) | 0.342275 (1.2%) |
| si8-standard-T8 | 0.340028 (2.6%) | 0.302823 (1.7%) | 0.326305 (3.3%) |
| si16-strict-T1 | 5.079188 (16.2%) | 5.444892 (112.8%) | 4.964811 (30.0%) |
| si16-strict-T2 | 3.509863 (0.7%) | 3.517435 (43.6%) | 3.391341 (2.0%) |
| si16-strict-T4 | 2.884367 (0.2%) | 2.815411 (52.4%) | 2.742652 (0.8%) |
| si16-strict-T8 | 2.697119 (84.6%) | 2.894511 (101.7%) | 2.616321 (279.0%) |
| si16-standard-T1 | 1.972351 (0.7%) | 1.970106 (4.8%) | 2.031792 (3.1%) |
| si16-standard-T2 | 1.445106 (0.7%) | 1.566131 (17.3%) | 1.519788 (60.7%) |
| si16-standard-T4 | 1.239082 (2.2%) | 1.295939 (17.9%) | 3.165797 (83.1%) |
| si16-standard-T8 | 1.158329 (11.5%) | 1.207770 (18.9%) | 1.224323 (116.5%) |

### Paired-batch spread and ablation limits

| Case / profile / T | Before spread % | Workspace spread % | Final spread % | Final / before | Final / workspace |
|---|---:|---:|---:|---:|---:|
| si8-strict-T1 | 1.0 | 2.2 | 0.7 | 0.9716 | 0.9949 |
| si8-strict-T2 | 1.2 | 4.3 | 1.6 | 0.9527 | 0.9942 |
| si8-strict-T4 | 9.1 | 11.4 | 8.3 | 0.9268 | 1.0106 |
| si8-strict-T8 | 1.7 | 5.5 | 2.9 | 0.9296 | 0.9833 |
| si8-standard-T1 | 1.0 | 1.5 | 1.5 | 0.9268 | 0.9949 |
| si8-standard-T2 | 0.5 | 1.9 | 1.3 | 0.9090 | 1.0039 |
| si8-standard-T4 | 0.7 | 1.5 | 2.1 | 0.8964 | 1.0021 |
| si8-standard-T8 | 3.1 | 2.8 | 3.7 | 0.8765 | 0.9845 |
| si16-strict-T1 | 3.8 | 13.1 | 3.0 | 0.9818 | 0.9887 |
| si16-strict-T2 | 14.4 | 40.1 | 13.4 | 0.9699 | 1.0024 |
| si16-strict-T4 | 1.0 | 9.9 | 9.8 | 0.9520 | 0.9902 |
| si16-strict-T8 | 12.2 | 15.0 | 11.6 | 0.9455 | 0.9904 |
| si16-standard-T1 | 3.5 | 8.8 | 42.9 | 0.9388 | 0.9900 |
| si16-standard-T2 | 28.1 | 32.7 | 75.7 | 0.9105 | 0.9786 |
| si16-standard-T4 | 9.1 | 18.7 | 10.6 | 0.9065 | 0.9774 |
| si16-standard-T8 | 9.4 | 5.4 | 6.8 | 0.9049 | 1.0007 |

Ratios below 1 denote lower observed time. The individual rotation effect is mixed
and often smaller than spread; no acceptance claim is made for it or the combined
change. The implementation remains reviewable on this branch pending the explicit
allocation gate and dev08 ablation/acceptance. No default has changed.
