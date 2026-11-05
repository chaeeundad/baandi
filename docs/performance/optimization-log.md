# Optimization log (newest first)

## P14A — work accounting and rejected schedules (2026-11-05)

Derived from study note 1 H2/H6/H7; existing notes and outputs only, no QE source.
[Study note 3](../research/qe-study/davidson-work-accounting.md) and
[complete measured tables](../../tests/reference/perf/p14a/work-tables.md) contain
commands, binary/input hashes, the relayed QE calibration and all rejections.
No numerical policy or default is adopted. Optional source-only counters add
no numerical changes; the temporary tolerance patch is archived and removed.

Default local si8: 15 SCF iterations, 340 ordinary Davidson steps, 410 H block
calls and 2416 H vector applications including certification, 14 restarts.
The relay's 83 QE steps / 98 H blocks and 290 BAANDI steps / 2070 H vectors
are different units/builds and are kept separate. For factors 1/3/10/30/100,
local si8 SCF counts are 15/16/16/21/57 and total H vectors are
2416/2072/1794/1722/2280. No loosened factor preserves SCF count on any of the
four tested systems. Conditional H6/H2 implementation is therefore deferred;
this is not proof against all adaptive schedules.

One extra band gives si8/k444s/skewed/elongated H counts 2444/8310/3260/164;
skewed does not converge in 100 iterations. Two extras give 3249/9068/295/149;
si8 needs 17 iterations and skewed changes energy by +29.3839772182878 Ha.
Reject both as general defaults. Detailed per-iteration root and phase counts
are stored; buffer-origin fractions exclude the first empty band and explicitly
bound unattributed rank-loss/fallback directions. Timings are shared-host
observations, not accepted speedups. No new sampling profile was collected.

Final check.sh and README gates pass (0). Default Gamma/k222/k444s match the
pre-change run with exactly equal energies and density bytes; same-input
robustness comparison passes (0). Historical robustness comparison exits 5
because the stored synthetic Xx inputs/UPFs differ from current Si-labelled
fixtures; the old baseline is not changed. All four accounting cases are
bit-identical in density and energy at pure-Rust T=1/4. No MKL claim is made.

## Status index — 2026-11-05

See [current P10A–P10M status](README.md#reconciled-status--2026-11-05) and
[claim-to-evidence map](../validation/evidence-map.md). Entries below are historical
measurements and handoffs; “uncommitted”, “blocked” and “pending” describe their
recorded stage. P10I–M retained changes are now merged, but rejected QR and unmet
worker-utilization goals remain rejected/unmet. Current dev08 results are the
fdce016 Part 4 relay, not a local reproduction or blanket acceptance of all goals.

Every entry: profile finding, change, numbers before and after (same machine, `scripts/bench_scf.py`), equivalence evidence, rejected ideas. Rules: [README.md](README.md).

## 2026-11-03: profile of the starting point (macOS Apple Silicon, release build with symbols, 8-k-point case, `sample` for 8 s, 6445 samples)
Top of stack by category: FFT in `ndrustfft` lane transforms 48.0 percent; allocation (`malloc`, `free`, `memset`, `bzero`) 15.3 percent plus 15.0 percent attributed to `mach_absolute_time` (called inside the macOS `free` path, so it is allocation overhead too); `ndarray` iterator copies 5.9 percent; Hamiltonian apply and nonlocal projector accumulation 3.4 percent; `hypot` 2.1 percent (complex norm in a hot loop); orthonormalization 1.8 percent; GEMM 1.3 percent. About 30 percent of the time is allocation overhead, and the stacks show it comes from `ndrustfft::ifft_lane` and the per-call temporaries (`Fft::process` allocates its scratch on each call). Diagonalization is 91 to 97 percent of the engine time (R7).
Baseline (commit `6049217`, release, 1 thread): see baseline.json: Gamma 1.30 s, 4 k points 6.55 s, 8 k points 12.1 s (median of 3).
Planned steps (each gated by `bench_scf.py --require-identical` unless it changes the algorithm): O1 allocation-free FFT plan (preallocated scratch, `process_with_scratch`); O2 Hamiltonian apply without per-call temporaries and without `hypot`; O3 pruned FFT that skips lanes that are exactly zero for the wavefunction sphere (about 40 percent fewer lane transforms, bit-identical); O4 Davidson work: fewer Hamiltonian applications (tolerance schedule, subspace, restarts) only with the energy and certificate gates and measured equivalence; then re-profile.

## 2026-11-03: where the 10 to 100 times against QE comes from (code and output reading, before O3/O4)
Counts, not timings. QE 7.6, 8 k points, per SCF step: average Davidson iterations 3.9, 2.0, 3.0, 2.8, 6.5, 4.6, 5.0, 4.5, 1.5 (about 34 per k point in total; first threshold ethr 1e-2). rust-dft, same case: Davidson iterations per SCF step (summed over 8 k points) 320, 98, 95, 76, 80, 62, 67, 55, 32, 16, 19 (about 115 per k point; the first step alone 40 per k point; first eigen tolerance 1e-5, random initial guess). In addition `Davidson::solve` applies the operator to the whole basis in every iteration (`apply(operator, &basis)`) and again to the Ritz vectors for the residuals, so with capacity 48 and block 8 it performs about 40 applications per iteration where 16 suffice. Hypothesis ranking (to be confirmed by measuring after each step): (1) redundant applications in Davidson (bit-identical fix, ticket P7P3 stage 1), (2) residual from stored products (stage 2), (3) SCF tolerance schedule and initial guess (about 3 times fewer iterations if QE-like), (4) cost per application (O1 done, O3 pruned FFT in progress).

## 2026-11-03: re-profile after O1 (merged as 78de2eb; same method, 5020 samples)
FFT lane transforms and butterflies 75.4 percent; Hamiltonian apply and nonlocal 4.8 percent; orthonormalization 2.6 percent; malloc, free, memset and copies 2.5 percent (was about 30); hypot 2.5 percent; GEMM 2.1 percent. Allocation overhead is gone; the FFT now dominates, so the next step is O3 (pruned FFT, ticket P7P2) and, after it, a batching of strided lanes if the profile still shows cache misses; O2 (hypot, temporaries) is small now (about 5 percent together).

## 2026-11-03: O1 / P7P1-FFT-ALLOC — persistent FFT scratch

Starting evidence: the pre-change profile above (6445 samples, 8 s) attributes
48.0% to FFT lanes and about 30% to allocation, including `Fft::process` scratch
and ndrustfft lane temporaries. This task uses that supplied profile and the
stored baseline at `6049217`; it does not claim a fresh pre-change profile.
Implementation worktree starts at `55c2565c68d60ec2225cd7058720095513e7b6b3`,
macOS arm64 (Apple Silicon), release build, no threading changes.

Change: plan both directions once with rustfft; reserve the maximum in-place
scratch requirement across all six plans and one maximum-axis-length lane.
Copy input to output once, gather/transform/scatter strided axes, and transform
contiguous lanes directly. Axis order remains 0, 1, 2. Inverse multiplication
by `1 / axis_length` after each axis and final multiplication by grid size
reproduce the former ndrustfft operation sequence. Forward division by grid
size is unchanged. The two full-grid work arrays are no longer needed.
Public API, manifests, lockfile and numerical tolerances are unchanged.
ndrustfft remains only as the test oracle; its manifest disposition is deferred
to Claude (see dependency inventory).

Allocation verification uses the permitted documented argument plus debug
assertions, not a counting global allocator (which would require unsafe):
all Vec allocation and planning occurs in `new`; `transform` uses only borrowed
slices, fixed-length copies, indexing, arithmetic and rustfft's
`process_with_scratch`. Scratch length is asserted against each plan's
`get_inplace_scratch_len`; lane length is asserted against the current axis.
No resizing, planning, logging or I/O occurs in a transform. The regression test
also checks scratch/lane pointers, lengths and capacities across repeated calls;
these checks alone are not claimed to count allocations inside dependencies.

Validation: the historical ndrustfft implementation is retained in the test
module. Seeded random complex inputs on 1x1x1, 1x3x5, 3x1x5, 3x5x1, 7x9x11,
13x5x17, 15x20x24, 24^3, 25^3, 27^3 and 36^3 produce bit-identical real and
imaginary components in both directions, twice per direction (zero ulp).
All existing tests are unchanged. `cargo test -p dft-numerics -- --nocapture`
exited 0: 11 tests passed in 0.55 s; no heavy cases need ignoring. Maximum
relative round-trip error was 1.7644611924521969e-15; maximum relative Parseval
error was 3.5528748979513846e-15 (both dimensionless).

Rejected/deferred: collapsing inverse normalization into one scale would change
rounding; retained the original sequence. A counting allocator was not added
because this ticket requires no unsafe. Pruning, threading and Hamiltonian or
solver changes remain outside O1. README figures and artifact JSON are outside
this ticket's explicit owned files and were not modified.

Verification commands (all exit 0):

```sh
export PATH="/opt/homebrew/opt/rustup/bin:$PATH"
cargo test -p dft-numerics -- --nocapture
cargo test --release -p dft-numerics legacy_bitwise_equivalence_and_workspace_reuse -- --nocapture
scripts/check.sh dft-numerics
scripts/check.sh
cargo build --release -p dft-cli
python3 -I scripts/bench_scf.py --binary /Users/mino/Develop/rust-dft-wt/P7P1-FFT-ALLOC/target/release/dft-cli --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF --cases gamma,k222s --repeats 3 --baseline docs/performance/baseline.json --require-identical --record /tmp/P7P1-FFT-ALLOC-bench.json
```

Both check scripts passed build, test, fmt and Clippy with warnings denied.
Release oracle test passed in 0.01 s. Benchmark ran after local build/check jobs
finished, on macOS 27.0.1 arm64; the provided baseline identifies this machine
as an Apple Silicon Mac mini. Hardware-model sysctl was sandbox-denied.
Other jobs may run on this machine: Claude re-measures on a quiet machine at
acceptance. Ratios below are observations against the stored same-machine
baseline, not controlled speedup claims. Baseline was not modified.

| Case (30 Ry) | Before median wall (s) | After median wall (s) | Before engine / diag (s) | After engine / diag (s) | Iterations before / after | Energy difference (Ha) | Density | Observed ratio |
| --- | ---: | ---: | --- | --- | --- | ---: | --- | ---: |
| Gamma | 1.304109792 | 0.960475208 | 1.278177542 / 1.207626833 | 0.925481750 / 0.864972624 | 11 / 11 | 0 | identical | 1.36x |
| Shifted 2x2x2 (4 k) | 6.548742542 | 4.622876292 | 6.759375458 / 6.588469207 | 4.722267875 / 4.586680250 | 11 / 11 | 0 | identical | 1.42x |

Engine and diagonalization times are from the final repeat (script behavior),
not medians. After wall samples (s): Gamma `[1.268224375, 0.960475208,
0.929638417]`; shifted `[4.609614333, 4.622876292, 4.726622250]`.
Final energies (Ha): Gamma `-7.295530468210927`, shifted
`-7.9194234124690706`. Density SHA-256 values, identical to baseline:

- Gamma: `c4fde645365f4e7ed609486d76233335b3894540239b0443880351eea82257cd`
- Shifted: `7b7e9c5aaf457f4d6ca5de52c8f2bc9aca531e27c4f0cf7f334cad08e2515c72`

The UPF was read only from the supplied external path, never copied or sent.
These are rust-dft baseline comparisons, not new QE references. Changes are
left uncommitted for Claude's review and acceptance.

## 2026-11-03: P7P3 / O4a and O4b Davidson product reuse

Pre-change evidence (recorded before implementation): the supplied post-O1
profile above attributes 75.4% to FFTs and 4.8% to Hamiltonian/nonlocal work.
The supplied application-count analysis identifies repeated whole-basis and
Ritz applications. Starting worktree commit: `11cd2f1280cf402abc450bc86b31892ae389ce9b`.
Use this existing profile, plus a fresh local pre-change benchmark below;
no new sampling profile is claimed. Scope excludes FFT/Hamiltonian changes.

Implementation: stage 1 preserves products of unchanged columns. A fixed-capacity
cache and correction-output blocks of widths 1 through min(block_size,
capacity-bands) are allocated once per solve. Cache extraction replaces the
previous per-iteration application-result allocation; append application itself
adds no allocations. Restart products are computed from the reorthonormalized
basis, never rotated from old products. Existing allocation-owning GEMM, EVD,
orthonormalization and residual APIs remain; this is not an allocation-free
solver. Extra persistent storage is n*capacity plus n*b*(b+1)/2 complex entries.
Custom operator implementations remain responsible for their own allocations.

Stage 2 uses `(H V) C` for residuals, followed by a true block application when
all recurrence norms pass. Failed certification replaces that history entry
with the true maximum and continues with true residual corrections. The
additive `solve_with_recurrence(..., false)` selects stage 1, while `true`
selects stage 2. No Davidson fields, solver parameters or tolerances changed.
The result's additive `operator_applications` counts vector columns, not block
calls. No timers, logging or I/O remain in production solver code.

Timing evidence: local release builds, sequential execution, same Apple Silicon
host, three repeats per case. The pre-change local run and stage 1 both passed
`--require-identical` against the immutable baseline. Stored data, wall samples,
energies, hashes, UPF hash and detailed comparisons are in
`crates/dft-solvers/benches/p7p3-evidence.json`. No UPF was copied. These are
rust-dft comparisons, not QE reference measurements. Other jobs may share the
machine; these observations are not quiet-machine speedup claims. The initial
baseline run overlapped initial test compilation; acceptance timing must be
remeasured by Claude on a quiet machine.

| Stage 1 case | Local before wall median (s) | Stage 1 wall median (s) | Before engine / diag (s) | Stage 1 engine / diag (s) | SCF iterations | Energy delta vs baseline (Ha) | Density |
| --- | ---: | ---: | --- | --- | --- | ---: | --- |
| Gamma | 0.938335083 | 0.490146708 | 0.926858292 / 0.867829084 | 0.483907792 / 0.425435664 | 11 / 11 | 0 | identical |
| k222s | 4.828421834 | 2.224247959 | 5.052977416 / 4.915385669 | 2.217291375 / 2.086770875 | 11 / 11 | 0 | identical |
| k222 | 10.036017375 | 4.273151833 | 10.030164334 / 9.774817584 | 4.269255958 / 4.045469958 | 11 / 11 | 0 | identical |

| Stage 2 case | Stage 1 wall median (s) | Stage 2 wall median (s) | Stage 2 engine / diag (s) | SCF iterations | Absolute energy delta vs stage 1 / baseline (Ha) | Density |
| --- | ---: | ---: | --- | --- | ---: | --- |
| Gamma | 0.490146708 | 0.338361042 | 0.333877416 / 0.275249874 | 11 / 11 | 7.993605777301127e-15 | differs |
| k222s | 2.224247959 | 1.426848916 | 1.433234168 / 1.301882875 | 11 / 11 | 8.881784197001252e-16 | differs |
| k222 | 4.273151833 | 2.733060250 | 2.712224958 / 2.488674375 | 11 / 11 | 8.881784197001252e-16 | differs |

Engine/diagonalization times are the last repeat, per benchmark script behavior.
Application counts include final fixed-point certification solves:

| Case | Solves | Old applications | Stage 1 | Stage 2 | Old / stage 1 | Davidson iterations stage 1 / stage 2 |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| Gamma | 12 | 3281 | 1497 | 833 | 2.19172 | 95 / 95 |
| k222s | 48 | 17037 | 7351 | 3863 | 2.31764 | 484 / 484 |
| k222 | 96 | 33490 | 14309 | 7437 | 2.34049 | 955 / 955 |

SCF count collection used separate, disposable measurement builds: at each
successful solve, emit the result counter; for stage 1 also sum
`basis.cols() + bands` over iterations to reconstruct the old count. This sum
is exact for the old implementation on the identical trajectory; it is not a
separately instrumented pre-change SCF run. The test-only Counted wrapper
independently verifies both old and new actual counts. Measurement logging was
removed before final builds and is absent from all timing binaries. Detail
runs confirmed the same energies and density hashes as the timing runs.

| Stage 2 vs stage 1 | Max eigenvalue delta (Ha), all k points | Final max-residual delta (Ha) | Max pointwise density delta (electrons/Bohr^3) | Certificate density residual (Bohr^-3/2) | Certificate occupied eigenvalue change (Ha) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Gamma | 5.051514762044462e-15 | 2.4567275272909227e-14 | 1.27675647831893e-15 | 5.584022603033996e-10 | 6.126193441424732e-9 |
| k222s | 2.3869795029440866e-15 | 1.612687882389246e-16 | 3.608224830031759e-16 | 2.422697054270338e-10 | 7.370652144178536e-10 |
| k222 | 6.217248937900877e-15 | 1.7932287464015435e-15 | 2.914335439641036e-16 | 2.767961299901316e-10 | 1.1329614468813531e-9 |

All three detail runs report converged and retain their fixed-point certificates.
Synthetic oracle tests compare every floating-point bit of eigenvalues, vectors,
residuals and histories, plus iteration counts, to the frozen old implementation:

| Synthetic case (dimensionless) | Old / stage 1 / stage 2 applications | Stage 1 / stage 2 iterations | Stage 2 max eigenvalue delta vs stage 1 |
| --- | --- | --- | ---: |
| N=36 random dense, decaying diagonal | 216 / 72 / 40 | 9 / 9 | 6.21725e-15 |
| N=18 restart, deficient guess, capacity 9 | 424 / 284 / 161 | 42 / 42 | 1.77636e-15 |
| N=48 sparse kinetic, preconditioned | 280 / 128 / 72 | 15 / 15 | 1.21969e-13 |
| N=8 threefold degeneracy | 26 / 17 / 11 | 3 / 3 | 1.11022e-15 |

The stale-product mock first returns an identity product, then diag(2,3).
A zero recurrence residual is rejected by a true residual of 1. With a one-step
budget it returns NotConverged; with sufficient budget it converges at step 3,
history `[1, 1, 2.220446049250313e-16]`, 7 applications. All returned stage 2
synthetic eigenpairs are independently compared to dense solves and reapplied.

Rejected ideas and limitations: deriving restart products by rotation changes
rounding and was not used. A first oracle attempt using the existing test GEMM
operator failed bitwise comparison: GEMM changes its kernel with block width.
The bitwise oracle therefore uses genuinely column-independent dense products,
as required by the ticket's premise; custom block implementations with
width-dependent rounding cannot claim bitwise equivalence. Existing GEMM-based
tests remain unchanged and test numerical accuracy. No full-workspace
allocation rewrite, tolerance schedule, parameter, initial-guess or FFT change
was attempted. README, artifact JSON, tickets and baseline are explicitly outside
ownership and were not changed. Changes remain uncommitted.

Acceptance: stage 2 is enabled by default. All three SCF energy differences are
below 1e-12 Ha, certificates remain present, and all workspace convergence gates
pass. Stage 1 remains explicitly selectable for bitwise comparisons.

Verification commands and real exit codes:

- `cargo build --release -p dft-cli`: 0 for pre-change, stage 1 and stage 2.
- `cargo test -p dft-solvers product_tests -- --nocapture`: initially 101 for
  the width-dependent GEMM oracle described above; corrected column-independent
  oracle included in the following passing run.
- `cargo test -p dft-solvers -- --nocapture`: 0, 10 passed / 2 ignored,
  0.36 s test execution (stage 1 default at this point).
- `cargo test --release -p dft-solvers product_tests -- --nocapture`: 0,
  both oracle/certification tests pass (0.01 s, both stages explicitly selected).
- `scripts/check.sh dft-solvers`: 0, build/test/fmt/clippy all PASS (stage 2).
- `scripts/check.sh`: 0, build/test/fmt/clippy all PASS (stage 2).
- `cargo test --release -p dft-solvers -- --ignored --nocapture`: 0, both
  release numerical tests pass in 3.78 s (stage 2). N=200 max eigenvalue error
  6.039613253960852e-14, residual 9.632521390673855e-11, 94 iterations.
  N=300 preconditioned residual 1.4413223917690354e-11 / 52 iterations;
  unpreconditioned residual 9.904873174714745e-10 / 2269 iterations.
  Synthetic quantities are dimensionless.
- `python3 /tmp/p7p3-details.py stage1 /tmp/p7p3-count-stage1-cli` and
  corresponding `stage2 /tmp/p7p3-count-stage2-cli`: both 0; one independent
  run per case retaining results/density for the detailed comparison above.
  These are temporary measurement helpers, not repository entry points.
- `git diff --check`: 0.

Exact SCF benchmark invocation (run separately with pre-change, stage 1 and
stage 2 default; all three exited 0):

```sh
python3 -I scripts/bench_scf.py \
  --binary /Users/mino/Develop/rust-dft-wt/P7P3-DAVIDSON-HV/target/release/dft-cli \
  --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF \
  --cases gamma,k222s,k222 --repeats 3 \
  --baseline docs/performance/baseline.json \
  --require-identical --record /tmp/p7p3-stage1.json
```

For pre-change use `--record /tmp/p7p3-before.json`; for stage 2 remove
`--require-identical` and use `--record /tmp/p7p3-stage2.json`. To reproduce
stage 1 through the CLI, temporarily select `false` in the delegating `solve`
method and rebuild; the retained default selects `true`. No CLI or engine
source changes are necessary. The final release CLI is rebuilt without any
measurement instrumentation. No new dependency or unsafe code was introduced.

Final review reused the already-extracted Ritz coefficients for the residual
GEMM, eliminating a duplicate coefficient-block allocation. After this cleanup,
`scripts/check.sh`, the release product tests and the release CLI rebuild all
again exited 0. A final three-repeat benchmark with the same command above,
`--baseline /tmp/p7p3-stage2.json --require-identical --record
/tmp/p7p3-final.json`, exited 0: all three energies and density hashes are exactly
identical to the measured stage 2 run. Final medians were 0.345, 1.487 and
2.736 s; the short release product-test build overlapped its first Gamma repeat.
These final verification samples are also stored in the evidence JSON.

## 2026-11-03: O3 / P7P2-PRUNED-FFT — wavefunction lane pruning

Pre-change evidence: supplied O1 profile (5020 samples) attributes 75.4% to
FFT lanes/butterflies. Starting commit: `6fcbfe8a26d24fc399209bd213d9951fdd04a759`.
This task uses that measured profile, not a fresh sampling claim. A local
pre-change release benchmark is collected before modifying the implementation.

Implementation: `SpherePattern` precomputes C-layout lane starts for both
axis-0/1/2 passes and a support mask. Hamiltonian construction creates it from
`g_to_fft`; only the two wavefunction transforms in `apply_with` change.
Inverse skips absent (j,k) lanes and then absent k planes; forward skips
unrequested i planes and then unrequested (i,j) lanes. Dependencies cannot
cross those coordinates in subsequent passes. Forward entries outside support
are explicitly zero, not equal to full-transform coefficients. Inverse input
outside support must be canonical (+0,+0), matching scatter/fill; debug
assertions are fused into the required copy, with no separate scan. Empty
support tests verify that skipping zero lanes matches full FFT output bits.
No hot-loop allocation, logging, timers, unsafe, new dependencies, axis-order
changes or normalization changes. Full-grid methods and the Fft3d trait are
unchanged; density (including its wavefunction transform), Hartree, XC and ionic
paths are untouched under the ticket's explicit exclusions.

Bitwise tests: seeded random complex coefficients on radius sweeps
(empty, 0, 1, 3.5, 6, full grid), shapes 24^3, 25^3, 27^3, 15x20x24, 7x9x11,
and 1^3, twice per direction with scratch reuse. All inverse entries and all
requested forward entries match to_bits (including zero signs); forward
non-support entries are canonical zero. The Gamma and shifted-k synthetic
Hamiltonian tests retain the pre-O3 full FFT apply sequence as an oracle,
including nonlocal accumulation, and compare three seeded random vectors each.
Debug FFT tests took 0.89 s, Hamiltonian test 0.06 s; no new heavy test requires
ignoring. Release FFT tests took 0.04 s. Test-only counters increment at actual
pruned rustfft calls: on 24^3 radius 6, inverse and forward each execute
1001 of the full transform's 1728 lanes (576 per axis), ratio 0.579282407,
a 42.0717593% reduction. Counters are absent in production builds.

Rejected/deferred ideas: changing normalization or axis order would change
rounding and was rejected. Optional adjacent-lane batching was not implemented
or benchmarked; its cache benefit is unestablished and remains a separate
experiment. No batching performance claim is made. README generation and
artifact/session updates are excluded by this ticket's owned-file restrictions.

Benchmark environment: macOS 27.0.1 arm64, release, same worktree and host
as the local pre-change run, no threading changes. Other jobs may share the
machine; the final pre-change k222 repeat overlapped local test compilation.
The after run started after build/check jobs finished. Claude re-measures on
a quiet machine at acceptance. These are observed ratios, not controlled
speedup claims. Engine/diagonalization times are the final repeat, not medians.

| Case (30 Ry) | Stored baseline wall (s) | Local before wall median (s) | After wall median (s) | Before engine / diag (s) | After engine / diag (s) | Iterations before / after | dE vs baseline (Ha) | Density | Local before / after |
| --- | ---: | ---: | ---: | --- | --- | --- | ---: | --- | ---: |
| gamma | 1.304109792 | 0.907344667 | 0.643498083 | 0.902268791 / 0.843729832 | 0.638620542 / 0.577091207 | 11 / 11 | 0 | identical | 1.410x |
| k222s | 6.548742542 | 4.498554500 | 3.006559750 | 4.487862834 / 4.357106377 | 3.000548250 / 2.866770332 | 11 / 11 | 0 | identical | 1.496x |
| k222 | 12.080713459 | 8.870671458 | 5.798963250 | 12.010770708 / 11.752801791 | 5.793543875 / 5.564493375 | 11 / 11 | 0 | identical | 1.530x |

Stored raw measurements (wall samples in seconds, energies in Ha):

- gamma: before `[1.248156333, 0.907344667, 0.9067223329999998]`; after `[0.958626375, 0.6270296660000001, 0.6434980829999999]`;
  energy `-7.295530468210927`; density SHA-256 `c4fde645365f4e7ed609486d76233335b3894540239b0443880351eea82257cd`.
- k222s: before `[4.499430542, 4.498554500000001, 4.492647792]`; after `[3.082523917, 2.998330083999999, 3.006559749999999]`;
  energy `-7.9194234124690706`; density SHA-256 `7b7e9c5aaf457f4d6ca5de52c8f2bc9aca531e27c4f0cf7f334cad08e2515c72`.
- k222: before `[8.870671458, 8.793498875000001, 12.015688333]`; after `[5.841755833000001, 5.795453083000002, 5.79896325]`;
  energy `-7.8312300571210764`; density SHA-256 `a1ca9c0697400090be911e275e88ebe660c0416eb0eb6438d18ce68603910334`.

Verification commands (each exited 0):

```sh
export PATH="/opt/homebrew/opt/rustup/bin:$PATH"
cargo build --release -p dft-cli # before and after implementation
cargo test -p dft-numerics pruned -- --nocapture
cargo test --release -p dft-numerics pruned -- --nocapture
cargo test -p dft-physics --test hamiltonian pruned_apply -- --nocapture
scripts/check.sh # PASS build, test, fmt, clippy
git diff --check
python3 -I scripts/bench_scf.py --binary /Users/mino/Develop/rust-dft-wt/P7P2-PRUNED-FFT/target/release/dft-cli --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF --cases gamma,k222s,k222 --repeats 3 --baseline docs/performance/baseline.json --require-identical --record /tmp/P7P2-before.json
python3 -I scripts/bench_scf.py --binary /Users/mino/Develop/rust-dft-wt/P7P2-PRUNED-FFT/target/release/dft-cli --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF --cases gamma,k222s,k222 --repeats 3 --baseline docs/performance/baseline.json --require-identical --record /tmp/P7P2-after.json
```

The UPF was read only at the supplied external path, never copied or committed.
Comparisons are against rust-dft, not new QE references. Baseline, manifests,
lockfile and other ticket-owned files were not modified. Changes are left
uncommitted for Claude review.

## 2026-11-03: P7P9 FFT lane batching (H8): measured and NOT merged
Four adjacent strided lanes gathered into one staging block (W = 4; W = 8 and 16 measured and worse or not better), per-lane arithmetic unchanged, bit-identical to the current transform on all test shapes. Transform alone (Mac, shared): 24^3 1.02 times, 25^3 1.03, 27^3 1.01, 36^3 1.02, 15x20x24 1.04 (inverse; the forward transform 1.02 to 1.03); tiny grids 0.85 to 0.96. End to end on the 8-k-point case, alternating A/B of five runs each on the same machine: minimum 1.023 s against 1.019 s, median 1.069 s against 1.051 s (1.00 to 1.02 times). The gain is within the run-to-run noise and the added code (a staging path plus an oracle) is not worth it; the branch is kept in the worktree history of the task only. Reconsider if the FFT becomes dominant again or on a CPU where the strided gather is much more expensive (the dev08 x86 test was not run because the effect on the Mac is already negligible).

## 2026-11-03: P7P4 / O4c SCF schedule and initial subspace

Pre-change evidence: the supplied post-O1 sampling profile attributes 75.4% to
FFT work; the application/iteration analysis above identifies excessive early
Davidson work after product reuse. This experiment uses that existing profile;
no fresh sampling profile is claimed. Counters are collected in memory and
serialized after SCF, including separate certificate and per-k counters.
Local pre-change measurements are collected before changing the algorithm.

Decisions: retain both measured algorithm changes. The residual schedule uses
c=0.1 with the existing normalized density L2 measure in atomic units, floor
0.1*eigen_tolerance and cap 1e-2 Ha. The loose first request avoids solving an
unconverged potential precisely; 1e-2 is motivated by the supplied QE starting
threshold, but our residual-norm criterion is not asserted to be QE's energy
threshold. c=0.1 tightens the request to eigen_tolerance when the density residual
reaches the default 1e-7; it is an empirical scheduler, not a unit conversion.
No sweep or optimality claim is made for c. The existing requested-tolerance gate
and the actual occupied residual gate both still apply. The independent solve
still requests 0.1*eigen_tolerance and checks the same density/eigenvalue bounds.
The new explicit Boolean mode defaults to true; false restores the original
initial_eigen_tolerance and eigen_tolerance_decay semantics. No existing numeric
control or convergence gate changed.

Cold subspaces use exactly bands lowest-kinetic coordinate vectors, sorted by
kinetic energy (equivalently |k+G|^2), then original basis index. These columns
already have full rank, so zero random completion columns are normally needed;
existing seeded Davidson rank completion remains available. Existing
Rayleigh–Ritz handles this initial subspace. No atomic/UPF guess is used. Warm
starts are unchanged. The additive free_electron_initial_guess switch restores
the previous random cold start when false.

Capacity 48 and block size bands are retained. Davidson already selects only
unconverged bands for correction. Capacity 32 increases total applications;
64 produces exactly the same iterations, applications, energies and density
hashes as 48 on these cases, with larger allocated capacity. Block 4 reduces
applications slightly for some k cases but increases iterations, worsens Gamma,
and does not improve measured wall time consistently. These are rejected ideas,
not default changes. No solver implementation, dependency or unsafe changes.

All counts below include certificate solves unless labeled per-step. Independent
ablation runs use the same release build configuration and supplied external
read-only UPF. Temporary detail runs retain results at /tmp/p7p4/<stage>/<case>;
the tables here preserve their key numerical evidence. Before = old schedule
and random start, schedule = residual schedule only, guess = kinetic guess only,
combined = both, capacity32/64 and block4 = both with that one parameter changed.
The parameter experiments temporarily changed engine defaults and were restored.

| Stage | Case | SCF steps | Total Davidson iterations | Total vector applications | E (Ha) |
| --- | --- | ---: | ---: | ---: | ---: |
| before | gamma | 11 | 95 | 833 | -7.295530468210919 |
| before | k222s | 11 | 484 | 3863 | -7.919423412469071 |
| before | k222 | 11 | 955 | 7437 | -7.831230057121077 |
| before | k444s | 11 | 4207 | 31963 | -7.925976126648079 |
| schedule | gamma | 11 | 71 | 592 | -7.29553046821087 |
| schedule | k222s | 10 | 318 | 2604 | -7.919423412469079 |
| schedule | k222 | 10 | 577 | 4924 | -7.83123005712114 |
| schedule | k444s | 10 | 2547 | 21078 | -7.925976126648085 |
| guess | gamma | 11 | 78 | 690 | -7.295530468210897 |
| guess | k222s | 11 | 419 | 3168 | -7.919423412469081 |
| guess | k222 | 11 | 787 | 6145 | -7.831230057121135 |
| combined | gamma | 10 | 50 | 422 | -7.295530468210769 |
| combined | k222s | 10 | 235 | 1883 | -7.919423412469065 |
| combined | k222 | 10 | 453 | 3661 | -7.83123005712114 |
| combined | k444s | 10 | 1921 | 15015 | -7.925976126648079 |
| capacity32 | gamma | 10 | 53 | 462 | -7.295530468210799 |
| capacity32 | k222s | 10 | 225 | 2032 | -7.919423412469083 |
| capacity32 | k222 | 10 | 481 | 4086 | -7.831230057121113 |
| capacity64 | gamma | 10 | 50 | 422 | -7.295530468210769 |
| capacity64 | k222s | 10 | 235 | 1883 | -7.919423412469065 |
| capacity64 | k222 | 10 | 453 | 3661 | -7.83123005712114 |
| block4 | gamma | 12 | 98 | 510 | -7.295530468210892 |
| block4 | k222s | 10 | 375 | 1858 | -7.919423412469069 |
| block4 | k222 | 10 | 715 | 3653 | -7.831230057121113 |

Timing environment: macOS Apple Silicon, same worktree, release, no threading
changes, shared machine. Pre-change three-repeat benchmark overlapped the
separate count run; final three-repeat timing ran after local build/check work
finished. These are observations, not controlled quiet-machine speedup claims.
Stage detail wall times include cold-start effects and are not used for ratios.
The k444s before/combined single-run walls were 8.911642208 / 4.239271917 s.

| Case | Before median wall (s) | After median wall (s) | After samples (s) | dE vs post-O4 (Ha) | dE vs original (Ha) |
| --- | ---: | ---: | --- | ---: | ---: |
| gamma | 0.266296208 | 0.156326875 | [0.6753115000000001, 0.15632687499999998, 0.1558603329999999] | 1.501021529e-13 | 1.580957587e-13 |
| k222s | 1.103937459 | 0.568264000 | [0.5678403750000001, 0.5682639999999999, 0.5683251249999999] | 6.217248938e-15 | 5.329070518e-15 |
| k222 | 2.089890833 | 1.062432875 | [1.0624328750000003, 1.0621882919999996, 1.064506625] | -6.306066780e-14 | -6.394884622e-14 |

All three strict 1e-12 Ha gates pass. Neither stored baseline contains k444s;
its independently captured pre-change E=-7.925976126648079 Ha is identical to
the combined result (delta 0 Ha). This is a local pre-change comparison, not a
claim to a missing stored baseline. Every run converged below the 100-step
budget. Density hashes change as expected for an algorithm change. Gamma's
1.50e-13 Ha shift exceeds the roughly 1e-14–3e-14 Ha tight-tolerance noise in
validation.md R2 but is below 1e-12; its larger final density residual is
consistent with a trajectory-dependent, second-order SCF energy difference.
No threshold was relaxed and no failing energy comparison was hidden.

Per-step sequences below sum over k points and exclude certificate work.
Certificates (including any failed attempts) are listed separately. The final
counts above include all such attempts.

### gamma count and certificate evidence

- before Davidson iterations: `[29, 9, 8, 8, 7, 6, 6, 6, 5, 5, 2]`
- before vector applications: `[241, 83, 80, 73, 72, 56, 56, 46, 36, 33, 24]`
- before step 11 certificate: 4 iterations, 33 applications, density residual 5.5840226030339960e-10 Bohr^(-3/2), occupied eigenvalue change 6.1261934414247321e-09 Ha.
- combined Davidson iterations: `[5, 4, 5, 4, 4, 4, 4, 4, 4, 5]`
- combined vector applications: `[47, 35, 41, 34, 33, 40, 36, 33, 34, 40]`
- combined step 10 certificate: 7 iterations, 49 applications, density residual 6.9320805331347533e-08 Bohr^(-3/2), occupied eigenvalue change 3.2212220946181702e-08 Ha.
- Final requested tolerance 8.3773221955763188e-09 Ha, occupied residual 3.9046074818507481e-09 Ha, density residual 4.1634599737292189e-08 Bohr^(-3/2), dE -4.8849813083506888e-14 Ha, relative charge error 3.5527136788005009e-15.

First-step effect per k point (existing CLI k order):

- before: iterations `[29]`; applications `[241]`.
- schedule: iterations `[22]`; applications `[180]`.
- guess: iterations `[13]`; applications `[102]`.
- combined: iterations `[5]`; applications `[47]`.

### k222s count and certificate evidence

- before Davidson iterations: `[142, 54, 50, 42, 52, 31, 25, 24, 21, 9, 11]`
- before vector applications: `[1217, 395, 383, 329, 347, 271, 211, 189, 185, 97, 91]`
- before step 11 certificate: 23 iterations, 148 applications, density residual 2.4226970542703381e-10 Bohr^(-3/2), occupied eigenvalue change 7.3706521441785355e-10 Ha.
- combined Davidson iterations: `[25, 20, 20, 20, 20, 20, 20, 24, 20, 23]`
- combined vector applications: `[206, 172, 175, 151, 175, 160, 163, 180, 171, 177]`
- combined step 10 certificate: 23 iterations, 153 applications, density residual 7.9514072114933366e-10 Bohr^(-3/2), occupied eigenvalue change 1.4192803182710989e-09 Ha.
- Final requested tolerance 1.9944572382027465e-09 Ha, occupied residual 1.1071737237635544e-09 Ha, density residual 6.8382774302015153e-09 Bohr^(-3/2), dE 1.5987211554602254e-14 Ha, relative charge error 2.2204460492503131e-15.

First-step effect per k point (existing CLI k order):

- before: iterations `[32, 35, 37, 38]`; applications `[283, 312, 312, 310]`.
- schedule: iterations `[24, 27, 26, 30]`; applications `[218, 239, 236, 236]`.
- guess: iterations `[20, 19, 19, 19]`; applications `[129, 131, 131, 131]`.
- combined: iterations `[7, 6, 6, 6]`; applications `[53, 51, 51, 51]`.

### k222 count and certificate evidence

- before Davidson iterations: `[320, 98, 95, 76, 80, 62, 67, 55, 32, 16, 19]`
- before vector applications: `[2288, 739, 719, 642, 648, 532, 551, 425, 291, 155, 186]`
- before step 11 certificate: 35 iterations, 261 applications, density residual 2.7679612999013162e-10 Bohr^(-3/2), occupied eigenvalue change 1.1329614468813531e-09 Ha.
- combined Davidson iterations: `[57, 39, 40, 39, 40, 35, 40, 45, 40, 47]`
- combined vector applications: `[406, 331, 345, 310, 341, 297, 333, 345, 326, 364]`
- combined step 10 certificate: 31 iterations, 263 applications, density residual 3.6636979692294761e-10 Bohr^(-3/2), occupied eigenvalue change 3.1084336216125052e-09 Ha.
- Final requested tolerance 1.7768811374767566e-09 Ha, occupied residual 1.7142994349745229e-09 Ha, density residual 5.2367282855205105e-09 Bohr^(-3/2), dE -4.4408920985006262e-15 Ha, relative charge error 1.1102230246251565e-16.

First-step effect per k point (existing CLI k order):

- before: iterations `[29, 41, 42, 41, 42, 40, 42, 43]`; applications `[241, 288, 289, 305, 287, 294, 296, 288]`.
- schedule: iterations `[22, 23, 23, 22, 23, 22, 23, 25]`; applications `[180, 212, 214, 209, 209, 211, 212, 212]`.
- guess: iterations `[13, 16, 16, 26, 16, 26, 26, 16]`; applications `[102, 117, 117, 143, 117, 143, 143, 117]`.
- combined: iterations `[5, 7, 7, 8, 7, 8, 8, 7]`; applications `[47, 50, 50, 53, 50, 53, 53, 50]`.

### k444s count and certificate evidence

- before Davidson iterations: `[1264, 507, 458, 390, 413, 240, 256, 167, 183, 92, 76]`
- before vector applications: `[10247, 3345, 3106, 2641, 2776, 2163, 2172, 1335, 1457, 811, 727]`
- before step 11 certificate: 161 iterations, 1183 applications, density residual 2.9030864194408356e-10 Bohr^(-3/2), occupied eigenvalue change 4.1669112604836300e-10 Ha.
- combined Davidson iterations: `[225, 178, 173, 172, 170, 170, 167, 172, 168, 180]`
- combined vector applications: `[1666, 1385, 1415, 1259, 1400, 1288, 1356, 1337, 1379, 1372]`
- combined step 10 certificate: 146 iterations, 1158 applications, density residual 7.6732956439607595e-10 Bohr^(-3/2), occupied eigenvalue change 3.1090208046924417e-09 Ha.
- Final requested tolerance 4.8209214642495654e-09 Ha, occupied residual 4.8066611153274375e-09 Ha, density residual 5.5647399322741725e-09 Bohr^(-3/2), dE -8.1712414612411521e-14 Ha, relative charge error 1.5543122344752192e-15.

First-step effect per k point (existing CLI k order):

- before: iterations `[31, 39, 39, 33, 34, 40, 41, 38, 34, 40, 47, 38, 35, 37, 38, 37, 41, 34, 41, 37, 38, 40, 46, 46, 41, 45, 49, 42, 38, 43, 40, 42]`; applications `[277, 316, 310, 295, 294, 314, 339, 331, 296, 336, 332, 332, 311, 324, 334, 316, 318, 296, 333, 327, 312, 291, 331, 336, 340, 330, 330, 339, 323, 334, 334, 316]`.
- schedule: iterations `[24, 26, 26, 24, 25, 26, 27, 27, 25, 27, 30, 26, 26, 27, 27, 26, 26, 25, 27, 27, 26, 23, 28, 27, 27, 27, 29, 27, 27, 28, 26, 26]`; applications `[216, 240, 233, 227, 229, 236, 259, 257, 229, 248, 240, 248, 233, 247, 260, 241, 242, 233, 248, 258, 234, 216, 243, 239, 259, 239, 238, 259, 248, 243, 245, 235]`.
- combined: iterations `[5, 7, 8, 7, 7, 7, 7, 7, 8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8]`; applications `[47, 54, 52, 54, 54, 53, 52, 51, 52, 52, 51, 51, 54, 51, 51, 54, 54, 53, 52, 51, 53, 53, 52, 51, 52, 52, 52, 52, 51, 51, 52, 52]`.

Validation: synthetic small schedule energy difference 8.8817841970012523e-16 Ha;
full nonlocal synthetic comparison 9.2370555648813024e-14 Ha (both <=1e-12).
Tests check the final requested tolerance and unchanged convergence gates,
certificate presence, and deterministic kinetic ordering including a tie.
The unchanged divergence/iteration-budget test passes. The existing seeded
positivity-safeguard audit explicitly selects both legacy controls, because its
seed-specific assertion tests a particular trajectory. Its assertion is retained;
a new full synthetic test covers the new defaults and the same energy audit.
An exploratory capacity-64 release test initially failed that seed-specific
assertion (6 passed / 1 failed); the final release SCF suite passes 8/8. One
initial build failed for a missing LinalgError conversion in the new helper;
it was corrected by using the existing SolverError conversion.

QE evidence: compare_kpoints_qe.py with unchanged defaults passes for new Gamma,
k222 and k222s results against tests/reference/si30_auto/qe.out and
kpoints/{k222_c,k222s_c}/qe.out. Energy offsets are respectively
1.2808946152631506e-6, 1.2814394296789544e-6, 1.2812654675009583e-6 Ha/atom;
max eigenvalue differences 1.2836681463664567e-6, 1.9904196993442533e-6,
1.7723341061606671e-6 Ha; normalized weight differences zero. Against the
stored rust-dft k222/k222s reference results, maximum eigenvalue changes are
3.3778383978777526e-9 / 5.212579257118932e-9 Ha, and total-energy changes
-6.394884621840902e-14 / 5.329070518200751e-15 Ha: all below 1e-8 Ha.
These reuse real stored QE 7.6 evidence; no new QE run is claimed.

Reproduction and verification (real exit codes):

- `cargo build --release -p dft-cli`: 0 for all completed measurement stages
  and the final build; the initial helper-conversion error exited 101.
- `python3 /tmp/p7p4-details.py <stage> [gamma,k222s,k222]`: 0 for before,
  schedule, guess, combined, capacity32, capacity64, block4 and final. The helper
  imports `scripts/bench_scf.py`'s TEMPLATE and CASES, writes temporary inputs,
  runs the release CLI once per case and retains results.json/density.bin.
  Before, schedule, combined and final also include k444s. No UPF was copied.
- `scripts/check.sh`: 0, build/test/fmt/clippy all PASS.
- `cargo test --release -p dft-engine --test scf -- --include-ignored --nocapture`:
  final 0, 8 passed; exploratory capacity64 version 101 as explained above.
- `cargo test --release --workspace -- --ignored`: 0, all 14 ignored tests pass,
  including existing solver, engine, physics and CLI release cases.
- `python3 -I scripts/compare_kpoints_qe.py --ours
  /tmp/p7p4/final/k222/results.json --qe-out
  tests/reference/kpoints/k222_c/qe.out`: 0. Substituting k222s/k222s_c: 0;
  gamma/tests/reference/si30_auto/qe.out: 0. No overrides.
- `git diff --check`: 0. Final detail runs reproduce combined energies and
  counts exactly for all four cases; capacity64 densities are byte-identical
  to capacity48 for all three tested cases.

Before and final benchmark commands each exited 0 (substitute `before` for
`final` in `--record` for the pre-change run):

```sh
python3 -I scripts/bench_scf.py \
  --binary "$PWD/target/release/dft-cli" \
  --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF \
  --cases gamma,k222s,k222 --repeats 3 \
  --baseline docs/performance/baseline-after-o4.json \
  --record /tmp/p7p4-final.json
```

Changes remain uncommitted for review. README, artifacts (including the session
registry), tickets, manifests, baselines and all other crates remain untouched
under the ticket's explicit exclusions. Final k444s count verification overlapped
release tests; use the earlier isolated combined detail wall quoted above.

## 2026-11-03: P7P5 O2 investigation

Pre-change evidence: supplied post-O4c 32-k-point profile, 2330 samples,
hypot 12.2%, orthonormalization 8.3%, GEMM 8.7%. No fresh profile claimed.
Residual products already call the shared norm, so require no solver edit.

Candidate: four independent norm accumulators, sqrt only when the squared sum
is finite and exceeds `MIN_POSITIVE * max(1, 2*n)`. Each subnormal product loses
at most MIN_POSITIVE*EPSILON, so this bounds aggregate underflow loss relative
to the sum by EPSILON. Otherwise the original hypot fold is used. Rounding
changes at roughly machine precision (1e-16 order for short vectors; accumulation
error depends on length); results are not bit-identical. Projection dots use
four independent complex accumulators; subtraction and both MGS passes remain.
No solver edit is needed: residual_from_products already imports this norm.
No new dependencies, unsafe, production allocations or public API changes.

Separate measurements (Apple Silicon shared host, release, 3 repeats; seconds):

| case | local before | norm only | norm + projection |
| --- | ---: | ---: | ---: |
| gamma | 0.157 | 0.180 | 0.162 |
| k222s | 0.572 | 0.531 | 0.511 |
| k222 | 1.066 | 1.736 | 0.955 |

Raw benchmark records and final detail results are retained in
`crates/dft-numerics/benches/p7p5-evidence.json`. Shared-host noise prevents a
controlled SCF speedup claim; Claude must remeasure on quiet dev08. An ignored
release kernel test (2048 complex entries, 20000 calls) measured old/new norm
0.689414/0.014653 s and old/new projection 0.046215/0.042230 s; this run overlapped
Clippy and is indicative only. No fresh sampling profile was taken.

Oracle tests retain the complete previous MGS and hypot implementations.
Fourier columns have analytic condition numbers 1e6 and 1e10 before floating
point construction; degenerate, random, rank-zero against an existing block,
1e-300/1e300 and mixed-scale fixtures pass. Dropped sets match the oracle;
maximum reported Fourier Gram defect is 6.662e-16 (limit 2e-15).
NaN/infinity rejection remains in the original finite-input validation.
Tiny vectors still use hypot before initial scaling; normalized arithmetic can
change rounding, so tiny output columns are not claimed bit-identical.

Combined SCF energy deltas against O4c: gamma 2.665e-15, k222s 1.776e-15,
k222 4.441e-15 Ha (absolute); all retain 10 iterations. The supplied O4c JSON
has no k444s entry: that specific baseline gate cannot be evaluated. k444s
energy is -7.925976126648081 Ha, differing from stored reference by 8.882e-15 Ha.
QE script PASS for k222/k222s/k444s with unchanged tolerances: energy offsets
1.28143942746e-6 / 1.28126546661e-6 / 1.28167595959e-6 Ha/atom; maximum eigenvalue
offsets 1.99041970013e-6 / 1.77233410771e-6 / 2.10013302574e-6 Ha. Changes from
stored rust-dft eigenvalues are at most 3.378e-9 / 5.213e-9 / 3.656e-9 Ha,
respectively, all below 1e-8. These reuse stored QE evidence, not new QE runs.

Rejected/deferred: block GEMM requires copies/allocations through the present
adapter, so was not implemented; skipping the second pass has no established
bound here and was not attempted. The accumulator variants were retained as
review candidates based on kernel savings and the larger SCF cases, pending
quiet-host acceptance. Gamma has no measured wall-time improvement.

Commands: release builds and norm-only/combined benchmark scripts exited 0;
`python3 /tmp/p7p5-details.py` ran all three requested QE comparisons, exit 0.
`cargo test --release -p dft-numerics -- --include-ignored --nocapture` exited 0
(16 tests). Initial `cargo fmt --all --check` exited 1 before formatting.
Initial `scripts/check.sh` exited 1: build/test/fmt passed, Clippy required
`as_chunks` instead of `chunks_exact`; corrected without changing arithmetic.
A diagnostic baseline lookup exited 1 because k444s is absent, as noted above.
Explicit ownership excludes README, session artifacts and baselines; untouched.
Changes remain uncommitted.

Final `scripts/check.sh`: exit 0, all four gates PASS. Final release rebuild,
three-repeat benchmark and all three QE comparisons: exit 0. Final wall medians
(gamma/k222s/k222) 0.147/0.541/1.129 s; energies and QE offsets exactly reproduce
the preceding combined run. This variability reinforces that SCF performance
acceptance is pending quiet-dev08 measurement. Final evidence is stored alongside
the intermediate stages. The missing k444s O4c baseline remains an open gate.

## P7P6: deterministic threads over independent k-point solves (2026-11-03)

Pre-change release build with debug symbols, commit 672d03f, shared Mac:
`bench_scf.py --cases gamma,k222s,k222,k444s --repeats 1 --baseline
docs/performance/baseline-after-o4c.json --require-identical` passed (exit 0).
Wall / diagonalization seconds: Gamma 0.498 / 0.097, k222s 0.567 / 0.439,
k222 1.093 / 0.866, k444s 4.314 / 3.556. The first three match the stored
baseline energy and density exactly; that baseline has no k444s entry.
The phase measurements locate 77–82% of multi-k wall time in diagonalization,
consistent with the earlier sampling profiles in this log. A fresh `sample`
attempt was denied by the sandbox (cannot examine the running process).
`sysctl hw.perflevel0.logicalcpu` was also denied (Operation not permitted).
Raw pre-change measurements are in ephemeral `/tmp/p7p6-before.json`.

Implementation: balanced contiguous disjoint mutable k-state groups, scoped
standard-library threads, ordered result collection after all joins. Only the
Davidson solve and its per-k preconditioner/initial-guess preparation run in
workers. The default serial loop retains the previous numerical operations.
Panic conversion adds no floating-point operations. No Send/Sync changes were
needed: Hamiltonian scratch remains private and non-Sync. All density, energy,
weight, mixing and certificate reductions remain sequential. The final
certificate calls the same threaded step. No unsafe or dependencies added.

Rejected: parallel density/energy reductions (change summation order), nested
FFT/faer parallelism (oversubscription and a different scope), defaulting to
all cores (changes defaults), global mutex around Hamiltonians (serializes work).

Final measurements use the release binary with debug symbols, with no concurrent
build jobs launched by this task. The machine is shared, not quiet; these are
indicative observed wall-time ratios, not accepted quiet-machine speedups.
Claude must remeasure on dev08. An earlier sweep overlapping debug compilation
was discarded for timing purposes. Gamma always uses one worker; its apparent
speedup below measures machine noise, not a parallel benefit. Core count remains
unavailable because the requested sysctl is sandbox-denied.

| Case | Threads 1 wall (s) | Threads 2 wall (s) | Threads 4 wall (s) | 1/2 ratio | 1/4 ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| gamma | 0.226836 | 0.157412 | 0.156562 | 1.441 | 1.449 |
| k222s | 0.713213 | 0.357061 | 0.265382 | 1.997 | 2.687 |
| k222 | 1.099054 | 0.674109 | 0.494169 | 1.630 | 2.224 |
| k444s | 5.425408 | 2.717299 | 1.882429 | 1.997 | 2.882 |

Gamma/k222s/k222: median of three wall samples per requested count. k444s:
one final sample per count. Diagonalization is elapsed region wall time including
spawn/join, not CPU time summed over workers. Energy differences are exactly
0 Ha and density hashes identical across all requested counts. All runs took
10 SCF iterations. Threads=1 also matches the unmodified starting binary;
Gamma/k222s/k222 match baseline-after-o4c.json, while k444s uses the pre-change
measurement because the immutable baseline contains no k444s entry.

Stored numerical evidence (seconds, Ha; raw benchmark files also in /tmp):

```json
{
  "gamma": {
    "energy_ha": -7.295530468210769,
    "density_sha256": "7a3cb36c709815636c659a2f88222c5de79a48a7d4e07e4213bd0e8697f99166",
    "threads": {
      "1": {
        "wall_samples_s": [
          0.28477916600000003,
          0.219016,
          0.22683616699999998
        ],
        "last_diagonalization_s": 0.14100746
      },
      "2": {
        "wall_samples_s": [
          0.172999167,
          0.15741175000000002,
          0.15667891699999997
        ],
        "last_diagonalization_s": 0.094446
      },
      "4": {
        "wall_samples_s": [
          0.16599124999999998,
          0.15656245800000002,
          0.15492520900000006
        ],
        "last_diagonalization_s": 0.094390377
      }
    }
  },
  "k222s": {
    "energy_ha": -7.919423412469065,
    "density_sha256": "ff47b9b53eed8c1fbf0b125186c652adff7f75770bd27c45ad164fa86776d61c",
    "threads": {
      "1": {
        "wall_samples_s": [
          0.723669875,
          0.695584166,
          0.7132131670000001
        ],
        "last_diagonalization_s": 0.553315667
      },
      "2": {
        "wall_samples_s": [
          0.357060917,
          0.3562234999999999,
          0.36276733400000016
        ],
        "last_diagonalization_s": 0.231115543
      },
      "4": {
        "wall_samples_s": [
          0.27199808299999995,
          0.265065958,
          0.26538183299999996
        ],
        "last_diagonalization_s": 0.135304458
      }
    }
  },
  "k222": {
    "energy_ha": -7.83123005712114,
    "density_sha256": "5d950d17264dab08a5e073667413a35d8e9bfacf5458b150d1a12800910f4c15",
    "threads": {
      "1": {
        "wall_samples_s": [
          1.1913966670000002,
          1.0990541250000003,
          1.0962359160000004
        ],
        "last_diagonalization_s": 0.86888554
      },
      "2": {
        "wall_samples_s": [
          0.6582013329999998,
          0.6876955000000002,
          0.6741085410000003
        ],
        "last_diagonalization_s": 0.451637457
      },
      "4": {
        "wall_samples_s": [
          0.4941693330000001,
          0.5002586250000001,
          0.48912695799999995
        ],
        "last_diagonalization_s": 0.266567333
      }
    }
  },
  "k444s": {
    "energy_ha": -7.925976126648079,
    "density_sha256": "d401e8861bbf4968c78c3e613874c89af212aaa90b3f74246dc264b62523738f",
    "threads": {
      "1": {
        "wall_samples_s": [
          5.425408417
        ],
        "last_diagonalization_s": 4.534472999
      },
      "2": {
        "wall_samples_s": [
          2.71729925
        ],
        "last_diagonalization_s": 1.9511803749999999
      },
      "4": {
        "wall_samples_s": [
          1.882429459
        ],
        "last_diagonalization_s": 1.071298748
      }
    }
  }
}
```

Reproduction (UPF is read-only; never copied):

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
CARGO_PROFILE_RELEASE_DEBUG=true cargo build --release -p dft-cli
python3 -I scripts/bench_scf.py --binary target/release/dft-cli --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF --cases gamma,k222s,k222 --threads 1 --repeats 3 --baseline docs/performance/baseline-after-o4c.json --require-identical --record /tmp/p7p6-final-t1.json
# Repeat with --threads 2 and 4, comparing against /tmp/p7p6-final-t1.json.
# k444s: --cases k444s --repeats 1, first compare threads 1 with the pre-change
# /tmp/p7p6-before.json, then compare threads 2 and 4 with threads 1.
```

All six final benchmark invocations printed OK; the benchmark sweep exited 0.
The synthetic engine sweep compares every numerical result field (timings
zeroed), untruncated density bits and every complex orbital component, across
requests 0,1,2,3,4,8 and 1,2,5 k points, including converged certificates.
Release energies: -0.62379772639925035, -0.59844111627273688,
-0.45953915387519850 Ha respectively. Serial reference uses the retained
pre-threading solve loop. Error injection checks failures at indices 2 and 3,
delays index 2, and verifies that index 2 wins and every solve completes;
a panic at index 2 becomes an indexed error. The normal debug sweep took
1.23 s; no new slow ignored test is required.

Scope exceptions: README, generated figures, artifacts and session registry
were intentionally untouched as explicitly required by this ticket. The ADR
0001 append uses its actual filename, 0001-numerical-backend-isolation.md.

Final gates: `scripts/check.sh` exit 0 (build, tests, fmt, clippy all PASS).
The first gate exited 1 because the new CLI test used unweighted explicit
multi-k rows, which the existing input contract correctly rejects; adding
explicit 0.5 weights fixed the test fixture. No production tolerance changed.
CLI end-to-end coverage verifies absent/0/1/2/8 requests, Gamma capping,
two-point capping, additive metadata and identical JSON numerical fields and
density bytes. `git diff --check` exit 0. Changes remain uncommitted.


## P7PA: deterministic per-k sequential phases (2026-11-03)

Step 0, before implementation: commit `66fbaaab55d5b3424994c484afdc8868d597e24a`,
shared macOS Apple Silicon host, 12 logical CPUs reported by Python. Release
build with `CARGO_PROFILE_RELEASE_DEBUG=true`. These are indicative timings,
not quiet-node speedup claims. `sample <pid> 1` could not examine the process
under the sandbox; profiling by stack sampling is unavailable. Existing phase
boundary timers identify setup, density and energy as the remaining large costs.
Certificate phases are included in phase totals; the certificate column below
is their sum, not an additional cost. Mixing is a shared-grid operation.

Commands: `python3 -I /tmp/p7pa-evidence/measure.py before T gamma,k222s,k222,k444s 3`
for T=1,8 (both exit 0). This temporary wrapper invokes `bench_scf.py.main`
with the requested binary, UPF, `--threads`, repeats and baseline; it preserves
results.json before benchmark temporary directories disappear. All four cases
have zero energy difference and identical density hashes against
`baseline-after-o2.json`, with 10 SCF iterations. Evidence is under
`/tmp/p7pa-evidence/before/` (temporary local evidence).

| case | threads | engine s | setup s | diagonalization s | density s | energy s | potential s | mixing s | certificate s |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| gamma | 1 | 0.137509 | 0.032735 | 0.080739 | 0.010742 | 0.007213 | 0.004707 | 0.001105 | 0.012708 |
| k222s | 1 | 0.516830 | 0.045352 | 0.390881 | 0.043455 | 0.029407 | 0.004946 | 0.002490 | 0.039487 |
| k222 | 1 | 0.975471 | 0.061845 | 0.755646 | 0.089469 | 0.060210 | 0.005360 | 0.002619 | 0.060935 |
| k444s | 1 | 3.911493 | 0.165713 | 3.117550 | 0.370894 | 0.247211 | 0.007006 | 0.002603 | 0.277586 |
| gamma | 8 | 0.139885 | 0.032369 | 0.082943 | 0.010935 | 0.007392 | 0.004867 | 0.001105 | 0.015489 |
| k222s | 8 | 0.255730 | 0.046011 | 0.126373 | 0.044319 | 0.030774 | 0.005240 | 0.002676 | 0.017129 |
| k222 | 8 | 0.376302 | 0.061207 | 0.157528 | 0.090035 | 0.059665 | 0.005099 | 0.002442 | 0.024326 |
| k444s | 8 | 1.394210 | 0.172382 | 0.611530 | 0.366458 | 0.234029 | 0.006838 | 0.002486 | 0.101098 |

Source inspection: density already forms per-k partials and folds them in k
order; this task need not change its arithmetic association. Certificate solves
already call the threaded iteration step. Preserve the Davidson call sites.

Implementation: scoped workers construct disjoint k states (basis, Hamiltonian,
projector tables), update per-k local potentials, extract orbitals and build
per-k densities, and evaluate per-k energy audits. Certificate work inherits
all these changes through its existing `step_kpoints` call. No Davidson call
site, solver option, physics API or backend implementation changed. No new
Send/Sync bound on a physics type, dependency or unsafe code was necessary.

One-thread/default execution spawns no worker. All counts use the same
ascending-k density and weighted energy folds. Workers return independent
partials; thread completion order never enters either sum. This preserves the
pre-change floating-point association: contrary to the ticket's initial
assumption, the merged code already used per-k partial density arrays.

Memory bound: W = min(max(threads, 1), nk). At most W partial density arrays,
each Ngrid f64 values (8 * W * Ngrid bytes), plus the single output array.
Each chunk is joined and folded before the next starts. W FFT plans are built
once and reused across density/audit phases and SCF/certificate steps; their
scratch and the per-worker audit/FFT temporaries add O(W * Ngrid) storage.
Existing per-k Hamiltonians and orbital copies remain O(nk). No new allocation,
logging or timer was introduced inside a per-band numerical loop; existing
physics routines retain their allocation behavior. Worker orchestration and
partial allocation occur at phase/k boundaries as required by the ticket.

Final release measurements below use three repeats for gamma,k222s,k222 and
one for k444s at T=1,4,8. Phase/engine values are from the last repeat (the
benchmark's convention), while wall values are medians. Builds/tests were not
running concurrently with this final sweep. Sixteen threads were omitted:
this host reports only 12 logical CPUs. Claude must measure the quiet dev08
node separately. Evidence: `/tmp/p7pa-evidence/after/t{1,2,3,4,8}/`.

| case | threads | wall median s | engine s | setup s | diagonalization s | density s | energy s | potential s | mixing s | certificate s |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| gamma | 1 | 0.140649 | 0.136564 | 0.032210 | 0.080262 | 0.010843 | 0.007216 | 0.004671 | 0.001097 | 0.012619 |
| k222s | 1 | 0.496518 | 0.492212 | 0.044398 | 0.368701 | 0.042912 | 0.028660 | 0.004841 | 0.002420 | 0.036831 |
| k222 | 1 | 0.926564 | 0.920209 | 0.059671 | 0.707968 | 0.087097 | 0.057617 | 0.005119 | 0.002439 | 0.060348 |
| k444s | 1 | 3.813364 | 3.807423 | 0.157285 | 3.010818 | 0.388579 | 0.240652 | 0.007016 | 0.002561 | 0.269538 |
| gamma | 4 | 0.141299 | 0.137272 | 0.032274 | 0.080614 | 0.010978 | 0.007326 | 0.004705 | 0.001103 | 0.012756 |
| k222s | 4 | 0.190203 | 0.185620 | 0.033581 | 0.118531 | 0.015471 | 0.010158 | 0.005143 | 0.002450 | 0.012588 |
| k222 | 4 | 0.319258 | 0.314982 | 0.038663 | 0.217621 | 0.030571 | 0.020210 | 0.005159 | 0.002453 | 0.020114 |
| k444s | 4 | 1.157091 | 1.151026 | 0.066341 | 0.870511 | 0.124302 | 0.081245 | 0.005725 | 0.002469 | 0.088866 |
| gamma | 8 | 0.142938 | 0.136663 | 0.032322 | 0.080186 | 0.010857 | 0.007252 | 0.004676 | 0.001103 | 0.012626 |
| k222s | 8 | 0.190885 | 0.186234 | 0.033583 | 0.119407 | 0.015198 | 0.010127 | 0.005169 | 0.002462 | 0.012592 |
| k222 | 8 | 0.237410 | 0.232455 | 0.035668 | 0.152525 | 0.022133 | 0.013794 | 0.005468 | 0.002515 | 0.013772 |
| k444s | 8 | 0.821099 | 0.814546 | 0.056650 | 0.605931 | 0.088522 | 0.054464 | 0.005977 | 0.002504 | 0.058581 |

At T=8, before/after time outside diagonalization (engine minus diagonalization;
includes setup, mixing and orchestration) is:

| case | before s | after s | reduction |
|---|---:|---:|---:|
| gamma | 0.056942 | 0.056477 | 0.82% |
| k222s | 0.129357 | 0.066827 | 48.34% |
| k222 | 0.218774 | 0.079930 | 63.46% |
| k444s | 0.782680 | 0.208615 | 73.35% |

Validation on the final binary:

- `CARGO_PROFILE_RELEASE_DEBUG=true cargo build --release -p dft-cli`: exit 0.
- `scripts/check.sh dft-engine`: exit 0 (initial implementation).
- `scripts/check.sh`: exit 0 on the final code; build, tests, fmt and Clippy
  with warnings denied all PASS.
- `cargo test -p dft-engine --test threading -- --nocapture`: exit 0, 1.06 s
  debug runtime. Synthetic nk=1,2,5 problems agree bit-for-bit for T=0,1,2,3,4,8
  in densities, energies, complete numerical histories and certificate data.
  Energies: -0.62379772639925035, -0.59844111627273666,
  -0.45953915387519950 Ha respectively. Timings are intentionally excluded.
- Added an order-sensitive fold test using [1e16, -1e16, 1, 2, 4]: every
  T=1,2,3,4,8 returns exactly 7.0 in global k order, including odd chunks and
  nk < threads. It also checks global k error indices in later chunks.
- Final benchmark wrapper invocations: `measure.py after T gamma,k222s,k222 3`
  and `measure.py after T k444s 1`, T=1,4,8; plus
  `measure.py after T gamma,k222s,k222,k444s 1`, T=2,3. All exit 0.
  Each invokes `bench_scf.py.main` with `--threads T --baseline
  docs/performance/baseline-after-o2.json`, worktree release binary and the
  read-only external Si UPF. No UPF was copied. Temporary wrapper only adds
  retention of results.json and benchmark JSON.
- `python3 -I /tmp/p7pa-evidence/verify.py`: exit 0. Across all four benchmark
  cases and T=1,2,3,4,8, maximum baseline energy difference is **0 Ha**,
  density hashes match each other AND the baseline, all iteration counts stay
  at 10, and complete numerical histories and energy records match T=1.

| case | energy Ha | density SHA-256 (all five thread counts) |
|---|---:|---|
| gamma | -7.2955304682107718 | `968fdd04bce7c2bfe30e227a9ea35e95a8971b203e20f18a412e552eff4454eb` |
| k222s | -7.919423412469067 | `40b373200b0be1668f781a809e81fb08fc50ded555948bf040742ea7b873874f` |
| k222 | -7.8312300571211448 | `402c5135f02d806c2420c5ef79260fbab09b90a8e103a581817a431b4f04112d` |
| k444s | -7.9259761266480808 | `4d6a4b558a8fd82f8bdd854e16bacfcf3a3e4255bc4f09550408e99d95bfeaed` |

QE comparisons use existing real-run evidence only, without reading QE source
or making a new QE reference claim. Commands all exit 0 and print PASS:
`python3 -I scripts/compare_kpoints_qe.py --ours
/tmp/p7pa-evidence/after/t8/CASE-results.json --qe-out
 tests/reference/kpoints/CASE_c/qe.out` for k222,k222s;
for k444s use `--qe-summary tests/reference/kpoints/k444s_c/qe_summary.json`.
Default tolerances were unchanged.

| case | energy difference Ha/atom | maximum eigenvalue difference Ha | largest change from stored comparison (Ha or Ha/atom) |
|---|---:|---:|---:|
| k222 | 1.2814394274585084e-6 | 1.990419700127612e-6 | 2.9275148293630678e-9 |
| k222s | 1.2812654666127798e-6 | 1.772334107711065e-6 | 5.197538273549009e-10 |
| k444s | 1.2816759595857263e-6 | 2.1001330257428846e-6 | 1.0373617408751258e-10 |

All comparison changes are below 1e-8 Ha, and all k weights match exactly.

Rejected/deferred ideas: completion-order sums and thread-group reductions
would break bit identity; retaining all nk density arrays would violate the
memory bound; nested FFT/solver threading is out of scope. Hoisting common
full-grid terms out of the per-k audit could remove additional redundant work,
but is deferred to keep this change to independent execution of the existing
validated audit. Shared-grid mixing has no independent k dimension and remains
sequential (about milliseconds here); no new mixing algorithm was introduced.
A thread pool could reduce spawn overhead but is outside the requested scoped
thread design. No QE study note was used.

Limitations: shared-host timings are indicative and peak RSS was not measured;
the memory bound above is structural, not an RSS claim. Stack sampling was
blocked by the sandbox. Raw local evidence under /tmp is temporary; the tables
and numerical checks above are the durable record. README, artifacts, manifests,
baselines, tickets and other tickets' source files were deliberately left
untouched under this ticket's explicit ownership. Changes are uncommitted.

## P7P7: Davidson hypotheses (2026-11-03)

Starting commit 1decf9c. Pre-change profile supplied with the ticket: FFT about
54%, orthonormalization about 8%, hypot about 12%, GEMM about 9%; approximately
53 Davidson iterations per k point across SCF and 8.1 vector applications per
iteration. This is existing measured evidence, not a fresh sampling claim.
Only the prose/mathematics study note was read; no QE source was accessed.
Stages below are derived from note sections 1/5 (restart and capacity), 2/5
(preconditioner), and 1/5 (Ritz convergence). Shared-Mac timings are indicative;
quiet dev08 remeasurement remains the reviewer's responsibility.

### Measurements and decisions

Raw medians, energies, density hashes, every SCF step's per-k-point Davidson
iteration/application counts, certificate counts, final bands, and QE checks
are stored in `crates/dft-solvers/benches/p7p7-evidence.json`. Counts in the
following table sum SCF steps and k points, excluding the separately recorded
certificate. Every stage starts from kept predecessors; S3/S4 use capacity 48.
Times are median seconds, three repeats, ordered Gamma / k222s / k222.

| Stage | Wall medians (s) | Davidson iterations | Vector applications | Max abs dE vs O2 (Ha) | Decision |
| --- | --- | --- | --- | ---: | --- |
| Before | 0.1404 / 0.4997 / 0.9269 | 43 / 212 / 422 | 373 / 1730 / 3398 | 0 | baseline |
| S1 contraction | 0.1464 / 0.5076 / 1.0986 | 43 / 212 / 422 | 373 / 1730 / 3398 | 0 | retain for restart savings; default cases neutral |
| S2 capacity 16 (2x) | 0.2430 / 0.5944 / 1.1100 | 62 / 340 / 667 | 452 / 2018 / 4076 | 1.19904e-13 | reject |
| S2 capacity 24 (3x) | 0.1540 / 0.4926 / 0.9399 | 52 / 260 / 522 | 428 / 1791 / 3615 | 1.30562e-13 | reject |
| S2 capacity 32 (4x) | 0.1400 / 0.4710 / 1.1011 | 46 / 202 / 447 | 381 / 1687 / 3408 | 3.46390e-14 | reject |
| S3 full smooth diagonal | 0.1545 / 0.5328 / 0.9619 | 52 / 243 / 438 | 427 / 1799 / 3461 | 1.16351e-13 | reject: Gamma 11 SCF steps |
| S4 Ritz selection + strict refinement | 0.1492 / 0.5218 / 0.9823 | 45 / 221 / 441 | 373 / 1731 / 3360 | 2.13163e-14 | reject: more iterations, slower |
| Final (S1 only) | 0.1509 / 0.5030 / 0.9298 | 43 / 212 / 422 | 373 / 1730 / 3398 | 0 | verified, no default-case speedup |

All stages except S3 Gamma have 10 SCF steps. No stage approaches the 5e-13 Ha
stop threshold. Smaller subspaces do not consistently improve wall time and
convergence together. The public `davidson_subspace` field remains unchanged.
Final k444s: 3.71665 s (one benchmark sample), 10 SCF steps,
-7.925976126648081 Ha; energy difference 0 Ha and density hash identical to O2.
The other final energies are -7.295530468210772 / -7.919423412469067 /
-7.831230057121145 Ha, also bit-identical energies and density hashes to O2.
Shared-machine variability, especially S1/S2 k222, precludes a speedup claim.

S1 is **derived from note section 1 (items 3/5), hypothesis H3 in section 5**.
The guarded restart reuses H products, saving one H application per retained
Ritz column at an accepted restart. Random n=18 and degenerate n=12 Hermitian
oracle fixtures measure net applications 165 -> 114 and 92 -> 64 (51/28 saved),
with unchanged iterations 47/21 respectively. Returned Gram defects are
3.33067e-15 and 4.21885e-15. The 1e-12 guard is over 200 times these measured
rounding defects; it is an accuracy ceiling, not a convergence tolerance.
Tests explicitly exercise the Gram failure and untrusted-product fallback.
An initial solver-test run exited 101: the existing deliberately stale operator
fixture exposed the need to invalidate products after failed certification.
That safeguard was added; the unchanged fixture now passes. The old recurrence
path remains an oracle through `solve_with_restart(..., true, false)`, and the
pre-O4 legacy oracle and its bitwise tests remain. S1 has zero application
savings on the default SCF cases because those solves fit without restarting.

S3 is **derived from note section 2, H5 in section 5**, rejected and removed.
The trial engine wiring assembled kinetic + grid mean of the effective local
potential + diag(B D B^H), including off-diagonal projector couplings using
existing physics accessors. No physics-crate edit was needed. A separate
preconditioner implemented the existing trait. The tested smoothing scale was
one Hartree (the note's expression evaluated directly in internal energy units),
which measurement rejected; this value is not adopted by the engine.
For x=d-lambda and t=x-1, D=1+(t+sqrt(1+t^2))/2 is strictly greater than 1
in exact arithmetic, hence has no finite zero/singularity. D/x tends to 1
as x tends to positive infinity; D tends to 1 from above as x tends to negative
infinity. For negative t, rationalization uses t+sqrt(1+t^2)=1/(sqrt(1+t^2)-t)
to avoid cancellation. Tests covered monotonicity, agreement with the direct
formula, finite extremes through +/-f64::MAX, asymptotes, output and invalid
input. Fourteen solver tests passed; no production preconditioner change remains.

S4 is **derived from note section 1 (item 4), H1 in section 5**, rejected and
removed. The trial considered root j stationary when the absolute change of
its sorted Ritz value between successive iterations was <= requested tolerance
(in Ha). The first iteration had no converged roots. This is an energy-error
estimator, not a rigorous error bound: for an isolated eigenpair near convergence,
eigenvalue error scales roughly as residual squared divided by the spectral gap;
small Ritz changes alone can also indicate stagnation. Once all roots were
stationary the trial resumed residual-based correction selection. Every return
still required a fresh true-residual certification with every band <= the
original tolerance. The fixed-point certificate retained residual selection
throughout; no empty-band exception or SCF threshold change was introduced.
The deliberately clustered slow-band test passed with actual residual
9.98312e-10 and eigenvalue error 1.13087e-12 (dimensionless), 90 iterations.
Fourteen solver tests passed. The extra refinement erased the prospective
iteration savings; k222 saved 38 applications but added 19 iterations. Rejected
rather than weakening certification. No convergence-option change remains.

S5 (optional H4) was not attempted: a generalized reduced solve would be a
separate high-risk extension; the earlier trials supplied no demonstrated
need to introduce it. No generalized solver/conditioning option is retained.
H2 empty-band policy and the SCF threshold schedule were explicitly out of scope.

### Final numerical/reference gates and reproduction

QE comparison script PASS for k222, k222s, k444s using the stored real QE runs
with unchanged tolerances, respectively:
energy offsets 1.28143942746e-6 / 1.28126546661e-6 / 1.28167595959e-6 Ha/atom;
maximum all-band eigenvalue offsets 1.99041970013e-6 / 1.77233410771e-6 /
2.10013302574e-6 Ha. Compared directly with stored rust-dft reference results,
maximum all-band changes are 3.37783961912e-9 / 5.21257520480e-9 /
3.65524610579e-9 Ha; energy changes <=6.83898e-14 Ha. Both are below 1e-8 Ha.
This reuses existing QE evidence, not a new QE run. Gap and fixed-point
certificate code are untouched. The benchmark success alone does not assert
iteration acceptance: S3's Gamma regression was checked separately and rejected.

```sh
cargo build --release -p dft-cli
python3 -I scripts/bench_scf.py --binary target/release/dft-cli \
  --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF \
  --cases gamma,k222s,k222 --repeats 3 \
  --baseline docs/performance/baseline-after-o2.json --record /tmp/p7p7-final.json
# Final large case: --cases k444s --repeats 1 --record /tmp/p7p7-final-k444s.json
python3 /tmp/p7p7-details.py final gamma,k222s,k222,k444s
python3 -I scripts/compare_kpoints_qe.py --ours /tmp/p7p7/final/k222/out/results.json \
  --qe-out tests/reference/kpoints/k222_c/qe.out
# k222s uses k222s_c/qe.out; k444s uses --qe-summary k444s_c/qe_summary.json.
scripts/check.sh
cargo test --release --workspace -- --ignored
cargo test --release -p dft-solvers -- --include-ignored --nocapture
```

The ephemeral detail helper uses bench_scf.py's unchanged TEMPLATE/CASES, runs
`dft-cli run <input.toml> --output-dir <out>` once per case, and retains the
results in /tmp/p7p7/<stage>/<case>/out. It does not copy the UPF. Permanent
JSON evidence retains the required histories and bands if /tmp is cleared.
Every completed release build, stage benchmark, detail helper and final QE
comparison exited 0. S1's first test failure (101) is described above; after
its correction S1/S3/S4 targeted release suites exited 0. The final workspace
build/test/fmt/clippy gates all print PASS. Final test counts are recorded below.

Acceptance limitation: the requested reduction of both applications and
Davidson iterations on the default benchmark is **not achieved**. Only restart
work is reduced on fixtures that actually restart. Rejected changes were removed;
no claim of a complete performance win or quiet-host acceptance is made.
README, artifacts, tickets, baselines, manifests, the parallel FFT work and all
other crates remain untouched under the explicit ownership exclusions.
Changes remain uncommitted for review.

Final release validation: all 15 ignored workspace tests passed; the solver
suite including ignored tests passed 15/15. Direct contraction product errors
on random/degenerate fixtures are 3.66384e-16 / 9.30536e-16 (dimensionless),
with measured restart Gram defects 8.88178e-16 / 4.44089e-16. Final complete
workspace check prints PASS build, test, fmt and clippy. The shell sequence
running that check, ignored workspace tests and complete solver tests exited 0;
logs show no failed test. `git diff --check` exited 0.
An initial evidence-assembly helper raised KeyError because Gamma has no
`kpoints` object; optional Gamma handling fixed it and the corrected helper
exited 0. This did not affect any numerical run.

## P7PB — intra-k scoped band workers (2026-11-03)

Starting commit: `20cb4168c4fea3c2886e23c0848e5d39a15d30a8`.
Before code changes, built release with debug symbols in
`/private/tmp/p7pb-before-target` (exit 0). Fresh k444s phase measurement:
diagonalization 3.014673292 s, density 0.369364208 s, engine 3.788528583 s.
`sample <owned-child-pid> 3 -file /private/tmp/p7pb-profile/sample.txt`
failed with exit 255 (sandbox cannot examine process); no fresh sampling profile
is claimed. The existing O1 profile above identifies FFTs as the dominant work;
this fresh phase measurement supports distributing independent band transforms.
Before medians (three runs, shared Apple Silicon Mac): Gamma 0.141 s,
k222 0.930 s, k222s 0.500 s, k444s 3.727 s. Energies exactly match
baseline-after-o2.json. The requested unmodified si8 input already fails in the
BEFORE binary: Davidson 400 iterations, residual 1.90890249000622e-4 Ha.
The combined benchmark exits 3 before writing its record. This is a pre-existing
acceptance blocker; tolerances and defaults will not be changed to conceal it.

### Implementation and invariants

The total budget is split deterministically into outer=min(max(N,1),nk) and
inner=max(1,max(N,1)/outer); unused remainders never become extra workers.
Immutable Hamiltonian data is shared inside each scoped call; each worker owns
FFT/grid/projector scratch and an output vector, allocated during construction.
The existing matrix interface cannot split mutable column ranges, so a mutex
protects only the completed-column copy. Recursion uses the calling thread as
one leaf, with at most inner active numerical workers; no join-handle Vec is
allocated. Scoped OS-thread bookkeeping still allocates; no numerical scratch
allocation occurs per block call. Threads=1 retains the original arithmetic.
Davidson already calls apply_block, so the solver file is unchanged.

Density workers reuse FFT/grid storage and per-band real contribution arrays,
then fold ascending bands and k points. Memory is O(inner*grid+bands*grid) per
k state, separate from existing orbital copies. Inner=1 keeps the prior density
function. Empty bands are not skipped. Errors are selected in ascending band
or k order after joins. No tolerance, algorithm, manifest, lockfile or production
FFT change; no unsafe, OpenMP, rayon calls or new dependency. This change is
based on the repository's own implementation and measurements, not QE source.

### Final paired measurements

Shared 12-logical-CPU Apple Silicon Mac, release binaries; no task-local build
or test ran concurrently with this matrix. Other users' jobs are uncontrolled,
so these are indicative measurements, not quiet-node acceptance. Median of
three full CLI runs per row. Both binaries use the starting checkout's same
compiler/dependencies; BEFORE was built before edits with a separate target
folder and `CARGO_PROFILE_RELEASE_DEBUG=true`; AFTER used the same profile.
UPF was read in place, SHA-256:
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
Records, commands with exit codes and logs are in `/private/tmp/p7pb-evidence`;
the tables below retain the key evidence if temporary files are cleared.
B and A denote BEFORE and AFTER; phase entries are medians too, in seconds.
Speedup B(N)/A(N) isolates the change at a fixed request; B(1)/A(N) also includes
the existing k-point parallelism.

| Case | N | B wall s | A wall s | B(N)/A(N) | B(1)/A(N) | B diag s | A diag s | B density s | A density s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| gamma | 1 | 0.139655 | 0.140629 | 0.993 | 0.993 | 0.079739 | 0.080354 | 0.010740 | 0.010858 |
| gamma | 2 | 0.140236 | 0.107186 | 1.308 | 1.303 | 0.079949 | 0.051194 | 0.010760 | 0.006105 |
| gamma | 4 | 0.140692 | 0.096700 | 1.455 | 1.444 | 0.080346 | 0.041663 | 0.010852 | 0.004531 |
| gamma | 8 | 0.140021 | 0.094070 | 1.488 | 1.485 | 0.079639 | 0.038982 | 0.010711 | 0.003924 |
| k222 | 1 | 0.924151 | 0.927180 | 0.997 | 0.997 | 0.707949 | 0.709786 | 0.086472 | 0.086747 |
| k222 | 2 | 0.498396 | 0.506250 | 0.984 | 1.825 | 0.363297 | 0.369889 | 0.046846 | 0.046652 |
| k222 | 4 | 0.324154 | 0.324703 | 0.998 | 2.846 | 0.220191 | 0.219915 | 0.030940 | 0.031256 |
| k222 | 8 | 0.240648 | 0.242731 | 0.991 | 3.807 | 0.155027 | 0.156542 | 0.022385 | 0.022172 |
| k222s | 1 | 0.498994 | 0.498479 | 1.001 | 1.001 | 0.370931 | 0.369851 | 0.043186 | 0.043422 |
| k222s | 2 | 0.275705 | 0.280145 | 0.984 | 1.781 | 0.188461 | 0.192006 | 0.022912 | 0.023000 |
| k222s | 4 | 0.190494 | 0.193798 | 0.983 | 2.575 | 0.118175 | 0.120059 | 0.015489 | 0.015453 |
| k222s | 8 | 0.192605 | 0.159950 | 1.204 | 3.120 | 0.119218 | 0.089168 | 0.015576 | 0.011966 |
| k444s | 1 | 3.701550 | 3.697095 | 1.001 | 1.001 | 2.954853 | 2.943547 | 0.346786 | 0.349675 |
| k444s | 2 | 1.918251 | 1.969330 | 0.974 | 1.880 | 1.500144 | 1.543661 | 0.187539 | 0.188197 |
| k444s | 4 | 1.181079 | 1.195208 | 0.988 | 3.097 | 0.890606 | 0.901763 | 0.125197 | 0.126499 |
| k444s | 8 | 0.826387 | 0.836033 | 0.988 | 4.428 | 0.615397 | 0.617486 | 0.087001 | 0.088015 |

All 96 successful CLI runs (2 binaries x 4 budgets x 4 cases x 3 repeats)
have identical densities within each case. All four energies are exactly equal
to baseline-after-o2.json (delta 0 Ha), with 10 SCF iterations in every run.

| Case | Total energy Ha | Density SHA-256 |
| --- | ---: | --- |
| gamma | -7.2955304682107718 | `968fdd04bce7c2bfe30e227a9ea35e95a8971b203e20f18a412e552eff4454eb` |
| k222 | -7.8312300571211448 | `402c5135f02d806c2420c5ef79260fbab09b90a8e103a581817a431b4f04112d` |
| k222s | -7.919423412469067 | `40b373200b0be1668f781a809e81fb08fc50ded555948bf040742ea7b873874f` |
| k444s | -7.9259761266480808 | `4d6a4b558a8fd82f8bdd854e16bacfcf3a3e4255bc4f09550408e99d95bfeaed` |

### si8 acceptance blocker and supplementary measurement

The added input has the conventional cubic lattice a=5.43 Angstrom and eight
diamond sites; no band override (20 bands), cutoff 30 Ry. The actual automatically
selected grid is 36^3, not the approximately 48^3 estimate in the ticket.
BEFORE and every AFTER budget fail with CLI exit 4 at Davidson iteration 400,
with the same maximum residual `1.90890249000622e-4` Ha. Three runs per row:

| Binary | N | Median time to failure s | CLI exit |
| --- | ---: | ---: | ---: |
| before | 1 | 5.150228 | 4 |
| after | 1 | 5.150535 | 4 |
| after | 2 | 4.094765 | 4 |
| after | 4 | 3.704369 | 4 |
| after | 8 | 3.576224 | 4 |

These are times to the same failure, **not** converged SCF speedups. There is
no si8 density/hash, total energy or phase report; its requested gates cannot
be certified. Changing solver defaults, subspace logic or tolerance to obtain
a passing result would exceed this ticket. Full acceptance remains blocked
on a separate investigation of the starting solver on this input.

For diagnostic evidence only, the ignored `conventional_cell_block_timing`
test uses the same cell/cutoff/grid (2969 plane waves, 20 random bands), a
synthetic local potential and no nonlocal projectors. Its output is bit-identical
across N=1,2,4,8. Median of 11 batches of ten block applications, seconds:

| N | Block median s | Speedup vs N=1 |
| ---: | ---: | ---: |
| 1 | 0.007561308 | 1.000 |
| 2 | 0.003976538 | 1.901 |
| 4 | 0.002514742 | 3.007 |
| 8 | 0.001662129 | 4.549 |

This is not an si8 SCF or a QE comparison. The full real-UPF si8 failure above
remains the acceptance result.

### Optional FFT-lane decision and rejected ideas

Rejected the packed-lane scoped-thread prototype; production FFTs stay serial
at every size (no enabling threshold). The ignored `fft_lanes` test gathers
independent lanes into reusable packed buffers, executes identical rustfft
arithmetic in scoped threads, then scatters them in each pass. It compares
full inverse transforms bitwise against the production plan before timing.
Final measurements below ran without another task-local build/test; an earlier
measurement overlapped check.sh and is excluded. Median of 11 batches of ten
transforms; serial and lane times are seconds.

| Grid | Lane threads | Serial s | Packed lanes s | Speedup |
| --- | ---: | ---: | ---: | ---: |
| 25^3 | 2 | 0.000111746 | 0.000144288 | 0.774 |
| 25^3 | 4 | 0.000111746 | 0.000167379 | 0.668 |
| 25^3 | 8 | 0.000111746 | 0.000255563 | 0.437 |
| 36^3 | 2 | 0.000302187 | 0.000271279 | 1.114 |
| 36^3 | 4 | 0.000302187 | 0.000316650 | 0.954 |
| 36^3 | 8 | 0.000302187 | 0.000368675 | 0.820 |
| 48^3 | 2 | 0.000808467 | 0.000829717 | 0.974 |
| 48^3 | 4 | 0.000808467 | 0.000860388 | 0.940 |
| 48^3 | 8 | 0.000808467 | 0.000853592 | 0.947 |

25^3 and 48^3 show no benefit; the actual si8 36^3 grid improves only 1.114x
with two lane threads and regresses with four/eight. This is less than the
1.901x two-band-worker block throughput improvement. No useful threshold or
additional benefit beyond band workers is established. Nested lanes would
also need to consume the same budget, not multiply it; no integration is kept.
Only the ignored reproducible experiment is retained, not a selectable backend.
This rejects this particular packing/spawning strategy, not all possible FFT
parallel implementations. Parallel small GEMM/orthonormalization was not added:
the measured diagonalization gain already comes from independent vectors.
A persistent custom executor was not introduced; scope lifetime and ownership
remain explicit. No unsafe mutable-matrix slicing was used to avoid the copy
mutex. Skipping empty density bands or re-associating sums was rejected because
bit identity, including the original zero contributions, is the contract.

### Reference and verification evidence

Synthetic complete SCFs compare all recorded fields after removing timing,
plus explicit bit checks of density, energies and eigenvectors; 1,2,3,5 k points,
3 bands, thread requests 0,1,2,3,4,8. Random seven-column Hamiltonian blocks at
Gamma and shifted k compare bitwise at 1,2,3,4,8 workers, including reuse after
multiple failures. Unit tests verify pool pointer stability and lowest failing
band 2 before band 5. Existing tests cover lowest failing k and panic joining.
CLI tests verify both additive budget fields and outer*inner<=max(N,1).
The 5 generated benchmark TOMLs validate (CLI exit 0); old default case inputs
remain byte-identical. The benchmark now verifies every repeat's energy/hash,
retains the old fields, and adds phase/engine medians.

Fresh AFTER N=8 outputs compared with the stored real QE outputs at unchanged
tolerances (not new QE runs):

| Case | Energy difference Ha/atom | Max eigenvalue difference Ha | Exit |
| --- | ---: | ---: | ---: |
| k222 | 1.2814394274585084e-06 | 1.9904197001276118e-06 | 0 |
| k222s | 1.2812654666127798e-06 | 1.772334107711065e-06 | 0 |
| k444s | 1.2816759595857263e-06 | 2.1001330257428846e-06 | 0 |

`check.sh` passed build, test, fmt and clippy (exit 0). Targeted physics/engine
suites passed (exit 0). The first targeted run failed to compile a test-only
`dft_units` import (101); using the existing physics data re-export fixed it
without a dependency. An initial FFT test launch preceded creation of its test
directory and exited 101; the corrected ignored test passed. A Python TOML
probe exited 1 because system Python lacks tomllib; validation through the real
CLI then passed for all five inputs. These were harness issues, not numerical
mismatches. Every paired benchmark command exited 0; the initial benchmark
including si8 exited 3, correctly reporting the CLI's exit 4. Profiling remains
limited by sample's exit 255 permission failure, as recorded above.

Reproduction (run from this worktree; preserve the starting binary before edits):

```sh
CARGO_TARGET_DIR=/private/tmp/p7pb-before-target CARGO_PROFILE_RELEASE_DEBUG=true cargo build --release -p dft-cli --offline
CARGO_PROFILE_RELEASE_DEBUG=true cargo build --release -p dft-cli --offline
python3 -I scripts/bench_scf.py --binary /private/tmp/p7pb-before-target/release/dft-cli \
  --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF \
  --cases gamma,k222,k222s,k444s --repeats 3 --threads 1 \
  --baseline docs/performance/baseline-after-o2.json --require-identical \
  --record /private/tmp/p7pb-evidence/before-1.json
# Repeat for each binary with --threads 1,2,4,8 separately, using before-1.json
# as the baseline and a separate record for each binary/budget.
python3 -I scripts/bench_scf.py --binary target/release/dft-cli \
  --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF --cases si8 --threads 8
cargo test --release --offline -p dft-physics --test hamiltonian conventional_cell_block_timing -- --ignored --nocapture
cargo test --release --offline -p dft-numerics --test fft_lanes -- --ignored --nocapture
python3 -I scripts/compare_kpoints_qe.py --ours /private/tmp/p7pb-evidence/k222/out/results.json --qe-out tests/reference/kpoints/k222_c/qe.out
python3 -I scripts/compare_kpoints_qe.py --ours /private/tmp/p7pb-evidence/k222s/out/results.json --qe-out tests/reference/kpoints/k222s_c/qe.out
python3 -I scripts/compare_kpoints_qe.py --ours /private/tmp/p7pb-evidence/k444s/out/results.json --qe-summary tests/reference/kpoints/k444s_c/qe_summary.json
scripts/check.sh
```

The temporary detail helper generates retained inputs using `case_input` and
runs `dft-cli run INPUT --output-dir OUT`. It records three direct si8 attempts
per budget so failures are measured without disguising them as passing results.
No remote actions, UPF copies, or writes to README, artifacts, tickets, baselines,
manifests or lockfile were made under the explicit ticket exclusions. Changes
remain uncommitted for review. Quiet-node remeasurement remains with Claude.

Final release ignored-workspace suite: 17 tests passed, zero failures,
`cargo test --release --offline --workspace -- --ignored` exit 0.
`git diff --check` exit 0. Full acceptance status remains blocked by the requested real-UPF si8 SCF
convergence/density gates. Quiet-host timing remains pending with Claude.

## P7PC: buffer-band diagnosis (2026-11-03)

Starting commit 18a55472d7c74518d9bef0b6a038764466d32a85. Derived from
`docs/research/qe-study/davidson-and-scf-thresholds.md`, H2; only the note was
read, never QE source. The policy thresholds below are local experimental
choices, not copied reference-code parameters. Existing P7P7 profile identifies
FFT 54%, orthonormalization 8%, hypot 12%, GEMM 9%; fresh sampling is attempted
separately. No timing claim uses instrumentation.

Unmodified release reproduces the supplied si8 input exactly: exit 4, iteration
400, maximum residual 4.238824601895987e-6 Ha. Temporary in-memory diagnostic
arrays (removed before implementation) record Ritz values and per-band residuals;
printing occurs after the solver loop. At failing SCF solve the tolerance is
2.285479797266155e-6 Ha. Iterations 1/10/50/100/200/400 have maximum residuals
8.214317027772144e-5 / 3.497194079315726e-5 / 9.263035888829632e-5 /
1.64303119215965e-5 / 1.33699801659869e-5 / 4.238824601895987e-6 Ha;
maximum bands (one-based) are 20/15/17/20/20/20. At iteration 400 the first
17 bands have maximum 9.060650089045546e-7 Ha, below the requested tolerance.
Thus strict bands have converged at the failure; buffer work blocks success.

Re-solving the SAME Hamiltonian with 24 bands and otherwise identical settings
converges. Bands 17--22 are 0.2507886854604968, 0.250791994322361,
0.250793791719534, 0.2507952913561576, 0.250796940751843,
0.250800679419172 Ha. Band 20 cuts a near-degenerate six-member cluster;
its lower/upper gaps are 1.4996366236e-6 / 1.6493956854e-6 Ha. Band 23 is
0.32277088203316945 Ha. The failed 20-band Ritz value is
0.2507968838486795 Ha, displaced within that cluster. This supports the
multiplet-cut hypothesis, without claiming a general theorem about stagnation.
The unconverged projected upper Ritz value (0.3214663826599302 Ha) was NOT a
reliable estimate of the missing neighbor; the expanded solve is necessary.

QE comparison inputs explicitly specify nbnd=8. The benchmark CLI template does
NOT explicitly specify the band count: it depends on the default extra_bands=4.
The default must therefore remain unchanged for this ticket.

Fresh `sample PID 1 -file /private/tmp/p7pc/profile.txt` exited 255: the sandbox
could not examine the process. The sampled unmodified CLI separately exited 4.
No privilege escalation was attempted. Profiling is therefore limited to the
previous measured profile and the new solver application/iteration counters.

Rejected trials (all preserve strict residual tolerances):
- `max(100*tolerance,1e-5)` with the original 48-dimensional space fails later
  on si8 (exit 4, max residual 1.215557559625791e-7 Ha).
- Correcting buffers down to the strict tolerance despite their loose stopping
  rule also fails (exit 4, max 9.786814290123058e-6 Ha); band 17 residual
  5.864361458095459e-6 exceeds the request 2.2939918872732104e-6 Ha.
- Floor-only `max(tolerance,1e-5)` at capacity 48 also fails: with strict
  correction selection max 1.548021051664704e-4 Ha; with per-band selection max
  1.3251144170067985e-5 Ha, strict band 17 at 4.411086909297155e-6 Ha versus
  request 6.37614022741499e-7 Ha. Loosening only the buffers cannot certify this.
- Capacity 60 plus the 100-times policy converges: -31.32180360103162 Ha,
  17 SCF steps, 701 Davidson iterations, 3195 main vector applications.
- Capacity 120 plus the 100-times policy converges: -31.32180360103174 Ha,
  17 SCF steps, 311 Davidson iterations, 2307 main vector applications.
- Capacity 120 plus floor-only policy converges: -31.321803601031377 Ha,
  15 SCF steps, 283 Davidson iterations, 2340 main vector applications.

Keep the floor-only target and six columns per band (minimum configured capacity
48). This is a measured extension needed beyond the initial tolerance-only
hypothesis, not a change to the requested number of bands. Existing eight-band
cases retain their capacity. The larger factor saved more applications in the
k cases but increased Gamma SCF steps 10->11 and si8 steps 15->17, so it was
rejected in favor of the more conservative threshold. All these numerical
choices were measured on rust-dft. No QE source was opened.

### Final policy, measurements and evidence

Raw diagnosis histories, complete final band residuals/certificates, benchmark
repeats/counts/hashes and comparison metrics are stored in
`crates/dft-solvers/benches/p7pc-evidence.json`. Release on the shared Apple
Silicon Mac, one thread, starting commit stated above, three repeats. Retained
final timing series did not overlap this task's builds/tests/other benchmarks;
other users' workloads are uncontrolled. Timings are indicative, not quiet-host
speedup claims. No remote server was modified; Claude remeasures on dev08.

Main-solve counts sum all SCF iterations and k points. Certificate applications
are separate; add both columns for total solver applications. SCF iteration
count differs from Davidson iterations.

| Case | Main applications before -> after | Certificate applications before -> after | Davidson iterations before -> after | SCF steps before -> after | Wall median before -> after (s) | Final energy (Ha) | abs dE vs O2 (Ha) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| gamma | 373 -> 281 | 49 -> 38 | 43 -> 41 | 10 -> 10 | 0.140166 -> 0.123313 | -7.2955304682107691 | 2.664535e-15 |
| k222s | 1730 -> 1266 | 153 -> 113 | 212 -> 204 | 10 -> 10 | 0.500443 -> 0.409511 | -7.919423412469083 | 1.598721e-14 |
| k222 | 3398 -> 2543 | 263 -> 208 | 422 -> 408 | 10 -> 10 | 0.935331 -> 0.787666 | -7.8312300571211368 | 7.993606e-15 |
| si8 | unavailable (failure) -> 2060 | unavailable -> 217 | failing solve 400 -> 283 total | unavailable -> 15 | failure 6.020170 -> 4.630403 | -31.321803601031377 | no O2 entry |

Every before-si8 repeat exits 4 with the reproduced residual. Its time is time
to failure, not converged SCF time; the old error exposes neither cumulative
application/SCF counters nor an energy, and these are recorded as null rather
than fabricated. The final si8 energy differs from the supplied
-31.3218036010316 Ha by 2.23821e-13 Ha (gate 1e-9 Ha). Seeds 1 and 7 at
threads 1 and 4 all produce -31.32180360103137673 Ha and 15 SCF steps;
the density files are bit-identical between thread counts for each seed.
The regression uses the full 30 Ry case (ignored in debug with its release
command), not an easier surrogate.

The floor-policy-only stage has main application counts 353/1554/3119/2340
and the final stage 281/1266/2543/2060 (Gamma/k222s/k222/si8). Skipping the
intermediate true-residual block therefore saves 72/288/576/280 applications,
with bit-identical energies and density hashes relative to that policy stage.
Only scheduled tolerances above the final eigen_tolerance use recurrence-only
acceptance. Final-accuracy maps, public fixed-point maps, fixed-point
certificates and all iteration-limit acceptance keep true residual checks.
The standalone solver reports whether residuals are certified.

The final si8 median was 4.630 s versus 3.724 s in the earlier skip trial;
a subsequent retained detail run was 3.724 s with identical energy/counts.
This observed shared-host variability is why no quiet-host speedup is claimed.

### Extra-band experiment (final policy, three repeats)

| Case | Extra bands | Main applications | Certificate applications | Davidson iterations | SCF steps | Median wall (s) | Energy (Ha) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| gamma | 1 | 236 | 65 | 59 | 11 | 0.119001 | -7.2955304682109068 |
| gamma | 2 | 245 | 37 | 54 | 10 | 0.118190 | -7.2955304682108348 |
| gamma | 4 | 281 | 38 | 41 | 10 | 0.123313 | -7.2955304682107691 |
| k222s | 1 | 905 | 96 | 214 | 10 | 0.300858 | -7.9194234124690688 |
| k222s | 2 | 1034 | 104 | 222 | 10 | 0.342279 | -7.9194234124690759 |
| k222s | 4 | 1266 | 113 | 204 | 10 | 0.409511 | -7.919423412469083 |
| k222 | 1 | 1766 | 212 | 463 | 10 | 0.577861 | -7.8312300571211413 |
| k222 | 2 | 1964 | 197 | 401 | 10 | 0.609176 | -7.8312300571211466 |
| k222 | 4 | 2543 | 208 | 408 | 10 | 0.787666 | -7.8312300571211368 |
| si8 | 1 | 2198 | 291 | 782 | 15 | 5.597247 | -31.321803601031441 |
| si8 | 2 | 2567 | 543 | 641 | 17 | 6.461023 | -31.321803601031618 |
| si8 | 4 | 2060 | 217 | 283 | 15 | 4.630403 | -31.321803601031377 |

Keep extra_bands=4. One/two extras reduce small-cell costs, but are slower on
si8 and one extra adds a Gamma SCF step. The QE inputs explicitly use nbnd=8;
the benchmark CLI template relies on the default instead of setting nbands.
`--extra-toml 'extra_bands = N'` appends an explicit override without changing
other default benchmark inputs. The existing si8 generator now preserves the
supplied input's atom order; the earlier generator's different atom order had
a different floating-point SCF trajectory. All one/two-extra small-cell energy
differences stay below 1.36e-13 Ha against O2; every si8 variant passes 1e-9 Ha.

### Reference and accuracy gates

No new QE runs: these are comparisons against the existing stored QE 7.6 runs,
using the unchanged comparison script/tolerances. One-based bands 1--5 are
strict and 6--8 are buffers; the ticket's four highest bands include strict
band 5, so both the three-buffer and four-empty maxima are reported.

| Case | Energy difference vs QE (Ha/atom) | Strict eigenvalue difference vs QE (Ha) | Buffer difference vs QE (Ha) | Four-empty difference vs QE (Ha) | Strict change vs stored rust-dft (Ha) | Script exit |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| k222 | 1.281439431455e-06 | 1.301039467504e-06 | 1.989492584469e-06 | 1.989492584469e-06 | 4.207432124925e-09 | 0 |
| k222s | 1.281265458619e-06 | 1.772201724705e-06 | 1.551901905068e-06 | 1.551901905068e-06 | 5.348791964899e-09 | 0 |
| k444s | 1.281675954257e-06 | 2.100463423490e-06 | 2.060307613584e-06 | 2.100463423490e-06 | 4.195325031819e-09 | 0 |

All buffer comparisons are below the stated 1e-5 Ha accuracy scale. This is an
observed comparison, not a theorem that residual norm bounds an indexed
cluster eigenvalue. The unchanged script also passes its original 1e-3 Ha
all-band eigenvalue gate. All strict changes are below 1e-8 Ha. k444s total
energy is -7.9259761266480915 Ha, 1.06581e-14 Ha from O2. No tested default or
extra-band result exceeds the 5e-13 Ha small-cell stop threshold.

The synthetic diagonal Hermitian fixture cuts an exactly degenerate multiplet
and has nearby levels. All-strict fails at 400 iterations (last maximum
1.4170027203180743e-11, requested 1e-12, dimensionless); buffer policy returns
in one iteration with strict residuals zero and buffer residual
2.7807830091299275e-8. Tests also require true residual recomputation, per-band
flags, iteration-limit acceptance with a failed buffer, rejection of a failed
strict band, invalid policies, and four saved applications for recurrence-only
intermediate acceptance. The small engine fixture verifies the strict/buffer
certificate partition and true final residuals. CLI tests verify additive
fields and preserve all existing keys.

### Verification and reproduction

`check.sh` passed build/test/fmt/clippy. The release ignored workspace suite
passed all 18 tests, including the four real-UPF si8 runs. Initial failures
were recorded and fixed: debug and ignored CLI exact-key assertions needed the
additive band_accuracy key (check exit 1, ignored suite exit 101); the first
external regression fixture omitted mandatory system.name (test exit 101).
None was a numerical tolerance failure. No tolerance was loosened to pass a
test. No dependencies, unsafe, physics/FFT/Hamiltonian code, baseline, ticket,
README or artifacts changes were made. An unrelated .DS_Store change appeared
during the run and was not edited or reverted by this task. Task changes remain
uncommitted under the worktree restriction.

```sh
cargo build --release --offline -p dft-cli
python3 -I scripts/bench_scf.py --binary target/release/dft-cli \
  --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF \
  --cases gamma,k222s,k222,si8 --repeats 3 \
  --baseline docs/performance/baseline-after-o2.json --record /private/tmp/p7pc/final.json
# Repeat with --extra-toml 'extra_bands = 1' or 2 and a different --record path.
DFT_SI_UPF=/Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF \
  cargo test --release --offline --workspace -- --ignored --nocapture
cargo test --release --offline -p dft-solvers degenerate_buffer_policy_and_iteration_limit -- --nocapture
scripts/check.sh
python3 -I scripts/compare_kpoints_qe.py \
  --ours /private/tmp/p7pc/final/k222/out/results.json \
  --qe-out tests/reference/kpoints/k222_c/qe.out
# k222s: corresponding k222s paths; k444s: --qe-summary tests/reference/kpoints/k444s_c/qe_summary.json.
```

For retained comparison outputs, generate TOML using `case_input` from
`scripts/bench_scf.py`, then run `dft-cli run INPUT --output-dir OUT`; this is
what the temporary detail helper did. The evidence JSON preserves those final
bands and certificates when temporary paths are cleared. The UPF was only read,
never copied. Session/dashboard/README updates are excluded by this ticket's
explicit ownership and left to the orchestrator.

## 2026-11-03: P10C pre-change investigation (not an accepted optimization)

Starting commit `0142ebf`; derived from the prose note
`docs/research/qe-study/fft-call-structure.md`, without reading QE source.
The current inverse and forward both run x/y/z, so the proposed symmetric
z-stick/xy-plane fusion changes inverse arithmetic order. Investigate it as
an explicit alternative, and retain a same-order implementation for the strict
bit-identity gate. Historical post-O4c profile (ADR 0019) attributes 54 percent
to FFTs. A fresh `sample <owned-child-pid> 1 -file
/tmp/p10c-evidence/sample.txt` on the symbol-bearing release pre-change binary
failed with exit 255 (sandbox cannot examine process); its k444s SCF exited 0.
No fresh sampling evidence is claimed. Pre-change benchmark and implementation
measurements are recorded in `batched-fft-study.md`; no speedup is accepted here.

### P10C implementation and local evidence

Implemented `FusedOperator`, default `OrderedFusedPlan`, and explicit
`PlaneFusedPlan` experiment. The default uses compact x sticks, active-z
inverse intermediates, and per-y xz planes that fuse inverse z / potential /
forward x. Requested x columns feed compact forward y/z passes. This preserves
both historical x/y/z orders without a full real wavefunction grid or
transform copy. The initial full-grid ordered prototype was superseded.
Density uses the inverse prefix; zero occupations skip FFTs. One plan is owned
by each existing worker; no internal parallelism, unsafe or dependency added.

Full oracle sweep: default local and density component bits are identical;
symmetric maximum relative errors are 6.12098052203517064e-16 and
1.15185970029673752e-15. Both sequential and 16-thread SCF benchmarks against
baseline-after-o2 pass `--require-identical`; complete Hamiltonian tests cover
Gamma/shifted bases, six grids and 0/3/5 synthetic projector channels. Workspace
build/test/fmt/clippy and the ignored release tests pass. Exact commands and
raw arrays are in `batched-fft-study.md`.

Local indicative medians (microseconds/call, before -> ordered): 25^3 local
145.975 -> 138.511, density 129.367 -> 78.674; 36^3 local 399.018 -> 373.790,
density 343.839 -> 196.438; 15x20x24 local 38.737 -> 30.097, density
37.934 -> 20.521. Complete Hamiltonian medians: Si Gamma 165.853 -> 146.687,
Si shifted 156.880 -> 146.576, Si8 505.183 -> 456.658. The Si8 default SCF
already fails before changes (max residual 1.90890249000622e-4 Ha after 400
Davidson iterations); no convergence setting was modified.

Five alternating single-thread SCF pairs: Gamma median 0.159570 -> 0.148165 s,
k222s 0.526594 -> 0.485518 s, k222 0.991584 -> 0.904057 s; synthetic skew
0.005434 -> 0.005471 s and slab 0.005057 -> 0.004637 s. Energies, all eigenvalues,
convergence histories (excluding timing) and density hashes match exactly.
These one-thread measurements use the generic density fallback with zero-band
skipping: `scf.rs` conditionally constructs fused density workers only for
`inner > 1`. Authorization for the minimal always-construct wiring is pending;
this file is outside the listed ownership. The 16-thread tests exercise both
fused operators together. Do not claim default single-thread fused density is
connected until that constructor change is authorized and verified.

Rejected variants: symmetric axis reversal passes the energy bound but changes
density hashes; it stays explicit. Actual stage-grouped batches 1/2/4/8/16/32
show no consistent gain across 25^3, 36^3 and 15x20x24; batching remains test-only.
Plane pruning was compared only at <=60% active columns and helped those sparse
supports; no tuned threshold or batching default was adopted. Shared-host
variability is material, a fresh profile was denied by the sandbox, and quiet
pinned dev08 measurements remain necessary before performance acceptance.

### P10C final integration on main 9d6fef1

The orchestrator approved the `scf.rs` constructor-only exception. Fused density
plans are now constructed for an inner budget of one as well as multiple
workers; P7PC tolerances, RMM-DIIS, mixing and thread policy are unchanged.
Main `9d6fef14ce86dce821bd359c6f58fef0aecad8ab` was integrated at source level
with backups and a three-way merge because Git could not write ORIG_HEAD in
its external worktree metadata (merge exit 128). The original branch metadata
is unchanged; the main-relative review patch is `/tmp/p10c-resume/task.patch`.
ADR 0026 supersedes this task's uncommitted 0021 numbering to avoid main's ADR.

The `FusedPlan` alias preserves P10B: optimized ordered Rust by default and a
reusable `GridFusedPlan<FftPlan3d>` compatibility adapter for the selected MKL
3D backend. It does not silently run Rust FFTs in an MKL build or claim an MKL
fusion speedup. Generic-adapter oracle and pointer-reuse tests pass. The optional
MKL-feature check encounters unchanged main's missing `Accum`/`matmul` imports
in `linalg.rs` (exit 101); no native execution or successful feature check is
claimed and that out-of-scope file is unchanged.

Final workspace `scripts/check.sh`: build/test/fmt/clippy PASS, exit 0. Release
fused integration tests: 5 passed including the full sweep; fused unit tests:
3 passed. Against `baseline-after-o4c.json`, final energy errors are 0,
1.7763568394002505e-14 and 3.552713678800501e-15 Ha (Gamma/k222s/k222), all within
1e-12 Ha. Historical density hashes already differ in unmodified main (strict
comparison exit 5). Against the freshly built main snapshot, both 1-thread
and 16-thread-budget runs pass `--require-identical` (exit 0), including Si8.

Three alternating main/final single-thread pairs also compare every eigenvalue
and convergence-history float, apart from timings, and density SHA-256 values:
all match for Gamma, k222s, k222, Si8, synthetic skew and synthetic slab. P7PC
now permits Si8 convergence at the default settings (15 SCF steps). Median
wall seconds, main -> final: Gamma 0.125360 -> 0.113772; k222s 0.404885 ->
0.359317; k222 0.757255 -> 0.670668; Si8 3.725412 -> 3.634331; synthetic skew
0.005038 -> 0.004968; synthetic slab 0.004410 -> 0.004790. Density phase medians:
Gamma 0.010644 -> 0.003619 s; k222s 0.042536 -> 0.014617 s; k222 0.084887 ->
0.029029 s; Si8 0.104090 -> 0.052569 s. Both synthetic engine/density phases
also decrease, but the slab process wall time increases; no universal wall
speedup is asserted. Shared-host samples are indicative and quiet-node profiling
remains with the orchestrator. Raw samples, hashes, exact commands, input
adaptation for the new preflight rules and remaining limitations are in the
latest section of `batched-fft-study.md`. The prior ownership blocker is resolved.

### P11A: incremental SCF, no numerical or measured performance regression

Base: main/task HEAD `9495d9fbaaf0e0bf6bf141ea61675e1d803eeeec`,
2026-11-03, arm64 macOS 27.0.1, rustc 1.99.1, default Rust backend.
No algorithm, numerical default, tolerance, reduction order or solver/mixer
implementation changed. The original loop body now advances a persistent
`ScfSession` once. Classic constructs no events/queue/timestamps; observation
publishes bounded scalar records at completed boundaries. Strict versus retained
budget termination preserves the original final-mixing difference. No speedup
or new scientific accuracy claim is made.

Before editing, built `cargo build --release -p dft-cli` (exit 0), saved the
binary as `/tmp/p11-before`, and recorded Gamma/k222s with
`python3 -I scripts/bench_scf.py --binary /tmp/p11-before --upf
/Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF --cases gamma,k222s
--repeats 5 --record /tmp/p11-baseline.json` (exit 0). A subsequent attempt to
profile the preserved baseline k444s process with `sample <pid> 1 -file
/tmp/p11-profile.txt` was denied by the sandbox (exit 255); the SCF itself exited
0. Thus no new stack-profile evidence is claimed. This is state extraction,
not a profile-derived kernel optimization.

Independent equivalence: temporarily compiled the original HEAD `scf.rs` and
`lib.rs` against the new test's identical synthetic fixture generator, captured
scientific fingerprints, then restored the final source. The temporary oracle
runner `/tmp/p11-oracle.py` and its release test exited 0. The final 18/18
fingerprints match the old loop. They cover three fixtures, three mixers and
threads 1/4. Both thread counts give the following same FNV-1a fingerprints of
normalized Debug scalars/history/certificates plus explicit complete density and
orbital bits (durations excluded):

| Fixture | Linear | Pulay | KerkerPulay |
| --- | --- | --- | --- |
| Cubic Gamma, nonlocal ion | b98923c3dfc5550f | d5725b937b95e719 | b1fba17f7710cadb |
| Cubic four k points | eaff564979a8d804 | 202f52bf37064600 | 1202a2f23a08e261 |
| Skew two species, local + nonlocal | b6c6edacfbc8454c | de548bc4df3e32ef | 5db516b4353c01aa |

The permanent release matrix additionally compares each of these 18 cases in
Classic, Observed and AgentControlled (54 comparisons) to the conventional
API, checking scientific records and each array element. Example Pulay energies:
-0.62379772639925168, -0.51412521505579167, -1.52491411485654815 Ha,
respectively; 11, 18, 13 steps. No inter-machine hash portability is claimed.

`cargo test --release -p dft-engine --test stepwise_equivalence --
--include-ignored --nocapture --test-threads=1` exited 0: seven tests passed,
including the matrix, retained/strict exhaustion, divergence, invalid mixer
initialization, slow-consumer loss accounting and critical-capacity exhaustion
before numerical work. A failed step and its events use one committed version;
terminal/repeated calls cannot advance it. The core has no durable journal:
critical events remain in bounded memory and resource exhaustion is explicit.
AgentControlled without intervention is measured; intervention/checkpoint work
is outside P11A.

Final `bench_scf.py --cases gamma,k222s --repeats 5 --baseline
/tmp/p11-baseline.json --require-identical --record /tmp/p11-after.json` exited
0: energy differences **0 Ha**, identical density SHA-256, ten steps each.
Gamma: -7.295530468210769 Ha,
`74963159f8eb058d58220f3bbde70a806a3622ca3c568b7d7d1ef068207610b0`.
Four k points: -7.919423412469083 Ha,
`c1de1d09cb8dc2c53500b10ce90e7ee8dc1bf692e1a87e2f2fec59cc64315617`.
These are internal before/after comparisons, not new QE reference runs.
Both old and new binaries also ran against `docs/performance/baseline.json`:
both strict comparisons exit 5 because historical density hashes already
changed on main; energy differences remain 1.58e-13 and 1.24e-14 Ha, within
1e-12 Ha. The historical baseline was not overwritten.

Timing: the first five-run measurement overlapped workspace compilation
(k222s wall median 0.490 s) and is not used for a performance conclusion.
After compilation completed, `/tmp/p11-paired-bench.py` called the repository
benchmark's `run_case` with alternating old/new order, one warmup pair and eight
measured pairs per case (exit 0; all energies and hashes identical):

| Case | Old/new wall median (s) | Change | Old/new engine median (s) |
| --- | --- | --- | --- |
| Gamma | 0.112248480 / 0.111934146 | -0.280% | 0.107682229 / 0.107233709 |
| Four k points | 0.359220604 / 0.359966521 | +0.208% | 0.354042292 / 0.354897500 |

Wall ranges overlap: Gamma old [0.111012125, 0.113199125] s, new
[0.111467375, 0.113269167] s; four-k old [0.357676541, 0.360752250] s,
new [0.357610583, 0.362651917] s. Engine median changes are -0.417%/+0.242%.
No measurable regression is resolved by these samples.

For a short synthetic Gamma workload, the release `mode_timings` test rotates
mode order, uses ten warmup rounds and 101 measured rounds, and includes
read-only snapshot assertions and explicit queue consumption. Median
Classic/Observed/AgentControlled times: 0.000375250 / 0.000375709 /
0.000375500 s; observed overhead versus Classic +0.122%/+0.067%. Ranges:
[0.000358000,0.000401916], [0.000358000,0.000390125],
[0.000357583,0.000420500] s. These agree within shared-host noise; they do not
establish universal sub-percent overhead or replace quiet-node P11J measurements.

Evidence logs and raw samples are temporary files `/tmp/p11-oracle.log`,
`/tmp/p11-tests-final.log`, `/tmp/p11-paired-bench.json`. The fixture/matrix and
interleaved mode timing command are retained in `stepwise_equivalence.rs`.
No new result dataset was added outside ticket ownership. Running
`python3 -I scripts/make_readme_figures.py` exited 0 and left generated README
and figures unchanged. No new README accuracy claim is warranted by this
behavior-preserving extraction. No default or optimization was selected on Si
alone, no kernel was modified, and no rejected numerical alternative was tried.

Verification gates: `scripts/check.sh dft-engine` and `scripts/check.sh` both
exit 0 (build/test/fmt/clippy PASS). Final engine check after adding snapshot
metadata and failure/queue tests also exits 0. README generator consistency:
`python3 -I scripts/tests/test_readme_figures.py` exits 0, three tests passed
(existing ResourceWarnings only). `git diff --check` exits 0. Main and task HEAD
remain identical at `9495d9f`; no rebase is necessary and all changes are left
uncommitted as required. Session/dashboard JSON and dataset publication are
outside this ticket's OWNED_FILES and remain with the orchestrator.

## P11J — matched runtime benchmark, local evidence only (2026-11-04)

No optimization or default change is adopted. Added an offline standalone driver
and paired Python harness, independently building pre-P11A `9495d9f` and candidate
`43868d9` with Rust 1.99.1/release on the same arm64 Mac. Three synthetic fixtures,
one-step exhausted and complete SCF, three modes, three warm-ups and 20 pairs:
360 measured pairs and 240 separate event/control/checkpoint/saturation probes.
All arrays and normalized scientific results match the pre-extraction reference.

The final collection's 1-minute load was 11.64 -> 14.50.
Across 18 mode/case cells, internal total timing yields {'not_met': 7, 'inconclusive': 10, 'met': 1};
external process wall timing yields {'inconclusive': 14, 'not_met': 3, 'met': 1}.
A local cell pass does not establish the target across cases or quiet-host
acceptance. No mode passes all required cases. Exact capture medians range
2.25–3.12 microseconds; synchronous persistence medians
6.69–7.53 milliseconds. Routine queue saturation and exact
checkpoint restoration preserve final numerical identity.

A first pilot unnecessarily initialized probe statistics inside normal timing.
That fixed adapter cost was removed before the final full collection. Complete
pilot data is retained, as is an intermediate collection explicitly invalidated
for editing its driver during build. Final source hashes are stable and verified;
no individual sample is excluded and physics/run counts/thresholds are unchanged.

Rejected inference: historical Si times cannot be denominators for synthetic
inputs or different hardware. Historical baseline JSON is embedded for provenance;
independent builds supply the matched denominator. Subtracting RSS does not prove
zero Classic allocations, and an entire step does not isolate enqueue/boundary
cost. These remain separate tracing requirements. No kernel was changed, no
profile-driven optimization is claimed, and no QE source or study note was used.
Detailed metrics and gaps: [validation results](../validation/runtime-equivalence.md#p11j-local-matched-measurements-2026-11-04)
and `tests/reference/runtime/benchmark.json`. Claude must remeasure on quiet dev08
via dft5090; full P11J acceptance remains blocked.

Historical Si follow-up: strict O4c density comparison exits 5 despite energy
errors 0 and 3.55e-15 Ha. Independently rebuilt pre-P11A and current CLIs have
identical Si Gamma/k222 energies and density hashes (exit 0), so the historical
density difference predates extraction. No gate or baseline was relaxed.
Evidence: `tests/reference/runtime/si-*.json` and logs. Workspace gates pass;
release checkpoint/runtime/stepwise suites: 31 passed, 0 failed.


## P10I — baseline audit; profiling blocked (2026-11-05 KST)

Status: BLOCKED; no numerical code, algorithm, default, or tolerance changed.
Base commit `ec37408cb08248cc2816068dee45d122f5509482`, pure Rust, Rust 1.99.1,
macOS arm64 (12 logical CPUs). Release build includes debug symbols. Binary
SHA-256: `8b1eff3f9b7ed69e1422c3c9e2ab587b464aade83dce7efe83bf7c05b3c31da7`.

The ticket's premise needs qualification: `Hamiltonian::apply_block` already
uses intra-k band workers, and `BandDensity` already transforms bands in
parallel and combines their contributions in ascending band order. SCF already
splits the total budget as outer=min(threads,nk), inner=max(threads/outer,1).
Davidson calls the block operator for its correction blocks and certificates.
Thus adding these mechanisms again would not address an established bottleneck.
The current solver also has sequential orthogonalization and linear algebra,
but source inspection alone does not establish their runtime share.

Attempted `sample <owned-child-pid> 1 -file /tmp/p10i/sample-tN.txt`
for N=1 and 8, each while running the unmodified si8 CLI. Both profiler calls
exited **255**: "sample cannot examine process ...". Both CLI runs exited 0.
The sandbox cannot supply the required stack profiles; escalation is unavailable.
Existing CLI phase timers combine FFT, projectors, GEMM, orthogonalization and
EVD in diagonalization, so they cannot justify selecting one for optimization.
No speculative kernel change or speedup acceptance was made.

### Local baseline evidence (not dev08 performance acceptance)

Three runs per case and budget. Builds/checks ran concurrently for part of this
collection, and this is a shared Mac; timings are diagnostic only and MUST NOT
be used for scaling or regression acceptance. All five cases have identical
energy, iteration counts and density hashes across budgets. Against historical
`baseline.json`, available energy differences are <=1.58e-13 Ha, but historical
density hashes differ (already documented for prior changes). The freshly built
T=1 binary is the matched identity baseline; historical baselines are unchanged.

| Case | T=1 wall median s | T=8 wall median s | T=1 / T=8 diagonalization median s | Energy Ha | Iterations |
| --- | ---: | ---: | ---: | ---: | ---: |
| gamma | 0.215141250 | 0.094937334 | 0.108906584 / 0.033509873 | -7.2955304682107691 | 10 |
| k222s | 0.699425958 | 0.151113750 | 0.494480331 / 0.078403168 | -7.919423412469083 | 10 |
| k222 | 1.138275709 | 0.196635792 | 0.845015124 / 0.120031917 | -7.8312300571211368 | 10 |
| k444s | 3.063932458 | 0.692012625 | 2.466737748 / 0.507720124 | -7.9259761266480915 | 10 |
| si8 | 4.539551917 | 2.918878584 | 4.205260290 / 2.658038043 | -31.321803601031377 | 15 |

Density SHA-256 (same for T=1 and T=8):

- gamma: `74963159f8eb058d58220f3bbde70a806a3622ca3c568b7d7d1ef068207610b0`
- k222s: `c1de1d09cb8dc2c53500b10ce90e7ee8dc1bf692e1a87e2f2fec59cc64315617`
- k222: `a38c265d381a9b918f667fbb14d43e210d3513387c9d3c41b1533930cc2515b3`
- k444s: `c081ea55e2265c416240e48f9bf275bd8e53d1a161648853a13b08f802d612c2`
- si8: `484f4f279ea968d7b0466f85130c7e33552cc37595becf09d29e306359c26698`

si8 uses 2969 plane waves, 20 bands and a 36^3 grid; energy is
-31.321803601031377 Ha. This is an engine result, not a new QE reference.
The >=3x dev08 target remains unverified, as do MKL and quiet many-k timing gates.

### Generality baseline probes

Synthetic Gamma fixtures use `scripts/bench_suite.py` and its existing synthetic
UPF generator, valence=2, channels=[0,1,2], shape="soft", cutoff=3 Ry,
seed=1, default tolerances, Pulay beta=0.5, max_iterations=100. For atom index
`i=0..n-1`, fractional positions are
`[(i*0.61803398875)%1, (i*0.41421356237)%1, (i*0.73205080757)%1]`.
Cubic 4-atom lattice rows: [[6,0,0],[0,6,0],[0,0,6]] Bohr.
Skew 20-atom lattice rows: [[8,0,0],[1,7,0],[0.5,1,6]] Bohr.
These deliberately small nonphysical fixtures test identity only; they do not
establish large-system performance or general default selection.

Both converge in 7 iterations at T=1 and T=8. Independent direct `density.bin`
byte comparison and exact energy equality assertions exit 0:

| Fixture | Energy Ha | Density SHA-256 |
| --- | ---: | --- |
| syn-20 | 1677.666321506571 | `bfafad504200d2be9ad5ddfe3450fe583257184388e268000252e5c3a10f62ac` |
| syn-4 | 20.108668260299748 | `cb637326288b6b4600227fb6dd16508a23c2e27be6b77d5c1bea826a4027c90b` |

The first 4-atom probe used a 4-Bohr cubic cell and failed input validation:
8 requested bands exceeded 7 plane waves. The invalid probe is retained under
`/tmp/p10i/general-t1`; enlarging the fixture to 6 Bohr makes the basis valid.
No solver setting/tolerance was relaxed. The valid T=1 suite exits 0.
The T=8 suite with `--baseline ... --require-identical` exits 5 because that
harness intentionally rejects different thread-count input fingerprints;
all individual calculations pass. Direct byte/energy comparisons above supply
the cross-thread check, without altering harness policy or baseline metadata.

### Commands and handoff

- `CARGO_PROFILE_RELEASE_DEBUG=true cargo build --release -p dft-cli --offline`: exit 0.
- `scripts/check.sh`: exit 0; build, test, fmt, clippy all PASS.
- `python3 -I scripts/bench_scf.py --binary target/release/dft-cli --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF --cases gamma,k222s,k222,k444s,si8 --threads 1 --repeats 3 --baseline docs/performance/baseline.json --record /tmp/p10i/baseline-t1.json`: exit 0.
- Same benchmark with `--threads 8 --baseline /tmp/p10i/baseline-t1.json --require-identical --record /tmp/p10i/baseline-t8.json`: exit 0.
- `python3 -I scripts/bench_suite.py --binary target/release/dft-cli --case-dir /tmp/p10i/cases --include-long --threads 1 --output /tmp/p10i/general-valid-t1.json --work-dir /tmp/p10i/general-valid-t1`: exit 0.
- Same suite with `--threads 8 --baseline /tmp/p10i/general-valid-t1.json --require-identical --output /tmp/p10i/general-valid-t8.json --work-dir /tmp/p10i/general-valid-t8`: exit 5, fingerprint mismatch described above.

Temporary evidence is under `/tmp/p10i` (inputs, reports, full CLI results and
logs); numeric summaries and fixture reproduction parameters are preserved here.
Claude/remote dft5090 session must collect si8 T=1/T=8 `perf` stacks on dev08,
separating H-apply FFT, projector operations, orthogonalization, GEMM/EVD and
density, then return the profiles for implementation. Do not change defaults
or relax across-thread density identity. Recheck both backends, many-k cases
and the 4/8/20-atom cubic/skew coverage after implementation.

Session registry and delegate report are absent for this task in this checkout.
Dashboard/session JSON publication belongs to the orchestrator and lies outside
this ticket's OWNED_FILES. No recurring orchestrator timer is active in this run.
No implementation was verified, so no README validation result was added.
Changes remain uncommitted as required.


### P10I resumed: architect-supplied profile (2026-11-05)

The architect supplied a Mac `sample` profile of release si8 Gamma at T=1
(~1550 samples); no new profiling is requested or performed. Self samples:
GEMM neon microkernel 553 (~36%), subtract_projection 264 (~17%), FFT ~270
(~18%), nonlocal add 93, fused local 60, GEMM wrapper 42, EVD ~26,
copies ~25. Reported diagonalization is 3.18/2.54 s of total 3.44/2.77 s
at T=1/8. This authorizes deterministic tiled GEMM and blocked projection work;
small EVD remains sequential. Reduction tiles must not depend on thread budget.
Changing projection rounding is authorized, with energy error <=1e-12 Ha and
unchanged across-thread density identity. No tolerance or default change.


### P10I implementation and final local validation (2026-11-05 KST)

Implemented against the supplied profile, without re-profiling or reading QE
source. Default pure Rust; baseline commit remains `ec37408`. All changes are
uncommitted. No SCF option, tolerance, band count, or physics default changed.

- `GemmWorkspace`: fixed 256-plane-wave tiles, borrowed faer submatrix views,
  or reusable packed MKL operands through `gemm_into`. Rotation tiles own
  disjoint result rows; Gram tiles own partial sums combined in ascending order.
  Tile boundaries, backend calls and combination order are independent of
  threads, including the one-thread path. EVD remains sequential.
- The existing Hamiltonian intra-k budget is exposed through an additive trait
  method (default 1), forwarded by the RMM counting adapter. Davidson's projected,
  rotation, residual-product and contracted-restart GEMMs and RMM's GEMMs use it.
  Outer k parallelism is unchanged; H applications and GEMMs do not overlap.
- Block projection against the established basis uses two GEMM projection
  passes, with internal two-pass MGS and normalization BETWEEN passes. The small
  incoming block and RMM's internal MGS remain sequential for stability.
- Numerical tile/projection buffers are allocated at solver-workspace creation
  and reused. `multiply_into` reuses caller output. Legacy Davidson's owning
  result allocations remain; this is not a wholly allocation-free Davidson
  rewrite. There are no new per-panel/per-band numerical allocations, logging,
  unsafe code, dependencies, or global backend thread settings. Scoped thread
  bookkeeping is transient, matching existing Hamiltonian workers.
- Panels with fewer than 65,536 complex multiply-adds per panel run sequentially.
  This is a scheduling-only overhead guard; it preserves exactly the same
  backend calls and reduction order. All five SCF cases at all four budgets
  match the pre-guard candidate's density hashes and energies exactly.

Rejected prototypes and failures (none hidden by changing tolerances):

1. Packing Rust operands into contiguous panels incurred ~5-8% many-k overhead.
   Borrowed views removed those copies; MKL still packs because its existing
   safe adapter requires contiguous operands. Small-panel worker launch overhead
   was measured separately: k222s T=8 median 0.129123458 s before the guard.
2. Projecting twice before normalizing the incoming block amplified near-rank
   roundoff. The existing RMM benchmark's 64-band quadratic-spectrum Davidson
   case failed at 800 iterations (max residual 5.934679923644799e-9 versus the
   unchanged 1e-9 threshold); an independently built HEAD archive passed.
   Interleaving normalization and reprojection fixed this. The full 8/16/32/64/128
   band benchmark then passed, and a direct nearly-dependent-column regression
   checks both cross-basis orthogonality and Gram defect <1e-12.
3. The frozen legacy product oracle is still untouched. Its former bitwise
   comparisons against the new solver are no longer valid after the explicitly
   authorized rounding change. The test now checks both solvers against the
   independent dense/residual/orthogonality oracle at existing tolerances and
   requires old/new eigenvalue agreement <1e-12. Product-count checks remain.
4. The optional release `mixing_scf` test fails before numerics because its old
   `Xx` synthetic species label is rejected by current input preflight. The
   unchanged baseline CLI reproduces that same failure. The fixture is outside
   ticket ownership and was not modified; it is not a passing gate.

### Final Mac timings

Same arm64 Mac (12 logical CPUs), Rust 1.99.1, release with debug symbols,
three CLI repetitions per case/budget. No build/tests ran concurrently with
these timing collections. Shared-host measurements are local evidence, not
quiet dev08 acceptance. Baseline is the preserved pre-change CLI; final is the
scheduling-guard candidate. Raw baseline/earlier candidate/final runs are retained
in `/tmp/p10i-impl`; compact raw wall samples are preserved below. The final
collection follows the paired baseline/borrowed-view collection, so small
percent differences must not be treated as universal speedups.

| Threads | Before wall s | After wall s | Before diagonalization s | After diagonalization s | After engine total s |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 3.537210500 | 3.754060500 | 3.255350291 | 3.472181499 | 3.744632666 |
| 2 | 3.100470667 | 2.622319292 | 2.832585957 | 2.348060039 | 2.604787584 |
| 4 | 2.920329542 | 2.025889708 | 2.663807334 | 1.768112831 | 2.015014209 |
| 8 | 2.837652167 | 1.811374875 | 2.583551210 | 1.560777418 | 1.801405124 |

si8: 8 atoms, Gamma, 20 bands, 2969 plane waves, 36^3 grid, 15 SCF iterations.
T=8 wall improves by 1.567x relative to the old T=8 run; final T=1/T=8 scaling
is 2.072x (old scaling 1.247x). **T=1 is 6.13% slower**, reflecting stronger
reorthogonalization and tiled arithmetic. The >=3x target is NOT established
on this Mac and remains for the remote dev08 session, including MKL execution.

Other si8 phase medians (seconds; do not sum medians to infer a measured total):

| Version | Threads | Density | Potential | Mixing | Total-energy audit |
| --- | ---: | ---: | ---: | ---: | ---: |
| before | 1 | 0.050949126 | 0.018816501 | 0.013368831 | 0.058602332 |
| before | 2 | 0.034198749 | 0.019253460 | 0.013719164 | 0.058786124 |
| before | 4 | 0.024026833 | 0.019063164 | 0.013776000 | 0.058746168 |
| before | 8 | 0.017388667 | 0.019132752 | 0.013764208 | 0.058855293 |
| after | 1 | 0.051450167 | 0.018968917 | 0.013451667 | 0.058550086 |
| after | 2 | 0.034102416 | 0.019132543 | 0.013884625 | 0.058702917 |
| after | 4 | 0.023874627 | 0.019084253 | 0.013930376 | 0.058681957 |
| after | 8 | 0.017071125 | 0.019105834 | 0.013835585 | 0.058728000 |

Many-k wall medians (seconds). No median regression is observed in this final
collection; differences of a few percent need quiet-node confirmation.

| Case | Threads | Before | After | Before / after |
| --- | ---: | ---: | ---: | ---: |
| k222s | 1 | 0.358561542 | 0.348447083 | 1.029027 |
| k222s | 2 | 0.206233917 | 0.203835625 | 1.011766 |
| k222s | 4 | 0.147871334 | 0.145943125 | 1.013212 |
| k222s | 8 | 0.126543375 | 0.124297375 | 1.018070 |
| k222 | 1 | 0.665486750 | 0.651468250 | 1.021518 |
| k222 | 2 | 0.368135709 | 0.362344000 | 1.015984 |
| k222 | 4 | 0.249241834 | 0.248104791 | 1.004583 |
| k222 | 8 | 0.178711500 | 0.178436250 | 1.001543 |
| k444s | 1 | 2.605961000 | 2.552044667 | 1.021127 |
| k444s | 2 | 1.373905791 | 1.344825583 | 1.021624 |
| k444s | 4 | 0.849184583 | 0.838687792 | 1.012516 |
| k444s | 8 | 0.603230959 | 0.588216167 | 1.025526 |

### Equivalence and generality

All five Si cases have exactly equal energies, density hashes and iteration
counts at T=1,2,4,8 within the new pure-Rust backend. The changed projection
rounding changes historical density hashes, as authorized. Energy differences
from the matched old build remain <=1e-12 Ha:

| Case | Final energy Ha | Absolute old/new delta Ha | Density SHA-256, all budgets |
| --- | ---: | ---: | --- |
| gamma | -7.2955304682107647 | 4.440892099e-15 | `1854d1ea91d96933ec85e908ec62b4466ddb044abe16b060c1cbe376daea0d3f` |
| k222s | -7.919423412469083 | 0.000000000e+00 | `2b310246ec91cf7ea8b845e03fcee0fe6c7769e5556c4af16ab5c8e4916e3ee5` |
| k222 | -7.8312300571211368 | 0.000000000e+00 | `f203190f8c973771213a5074eb769c3ceecf264945aa805ea2da6edad5988840` |
| k444s | -7.9259761266480915 | 0.000000000e+00 | `ad4d3eabe1b8823365908f6647c4470831a51416e8c738cd3f2aaabf7549b220` |
| si8 | -31.321803601031256 | 1.207922651e-13 | `9152c7a9c64db67f623379a343a7b8c84b532ccc08736f8870292ad287e3e4de` |

The stored `docs/performance/baseline.json` Gamma/k222 comparison also passes
(absolute energy differences 1.63e-13 / 6.04e-14 Ha). No new QE reference run is
claimed. Synthetic fixtures use the exact 4-atom cubic and 20-atom skew cell,
positions, cutoff, channels and generator parameters recorded in the earlier
P10I baseline audit above. Both converge in 7 iterations. These are deliberately
nonphysical small-basis identity fixtures, not realistic 20-atom timing claims.
The large/partial-panel GEMM tests and synthetic 8..128-band eigensolver matrix
supply separate large-band coverage. Direct density-file byte comparisons across
all four budgets exit 0 (the suite's normal fingerprint comparison includes the
thread count and is intentionally not used as a cross-budget comparator).

| Fixture | Final energy Ha | Absolute old/new delta Ha | Density SHA-256 |
| --- | ---: | ---: | --- |
| syn-20 | 1677.6663215065712 | 2.273736754e-13 | `4b88a4db6022ecc5375b78ae8a64b8213de14859d2330d23645f00f74966c5bc` |
| syn-4 | 20.108668260299716 | 3.197442311e-14 | `9971b193278611fad89aa2319a248979733825442421cee4a78ab8cd9053ac6c` |

### Reproduction and gates

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
CARGO_PROFILE_RELEASE_DEBUG=true cargo build --release -p dft-cli --offline
# Preserve the old executable BEFORE editing; never overwrite the reference.
# For t in 1 2 4 8, run each version with:
python3 -I scripts/bench_scf.py --binary target/release/dft-cli \
  --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF \
  --cases gamma,k222s,k222,k444s,si8 --threads 8 --repeats 3 \
  --baseline /tmp/p10i-impl/final-t8.json --require-identical \
  --record /tmp/p10i-impl/scheduled-t8.json
# final-tN is the pre-scheduling-guard candidate; matched-before-tN is old code.
# Cross-version energy comparison uses matched-before without --require-identical.
scripts/check.sh
cargo test --release -p dft-solvers --lib --test rmm -- --include-ignored
cargo check -p dft-cli --features mkl --all-targets --offline
python3 -I scripts/make_readme_figures.py
python3 -I scripts/tests/test_readme_figures.py
git diff --check
```

Final workspace check: exit 0, build/test/fmt/clippy PASS. Numerical unit tests
cover all adjoint combinations, scalar multiplication oracle, empty/partial
panels, nonfinite/dimension errors, stable scratch pointers, budgets 0/1/2/3/4/8,
rank filtering, extreme scales, and near-dependent reprojection. Every final
benchmark command (five cases x four budgets x three repeats) exits 0. Each
synthetic suite command for T=1,2,4,8 exits 0; direct old/new energy and
cross-thread byte assertions also exit 0. Historical baseline command exits 0.
README generation and its three consistency tests exit 0, with no generated diff.
No new reference JSON/figure type was added outside ticket ownership. Final release solver command exits 0 (16 library tests and 9 RMM tests,
including the full ignored benchmark). Final CLI MKL all-targets type check
exits 0; MKLROOT is unset, so no MKL executable was linked or run.

Optional `cargo test --release -p dft-solvers -- --include-ignored` exits 101
because of the pre-existing mixing input fixture described above; repeating
`--test mixing_scf -- --ignored` with
`DFT_CLI_BINARY=/tmp/p10i-impl/before-cli` also exits 101. This is not hidden by
excluding it from the record. No thresholds were relaxed.

Limits: real MKL linking/running and within-MKL density identity require dev08;
Mac feature checks cannot establish them. Quiet-node >=3x scaling and many-k
regression acceptance remain remote gates. Scheduler cutoff portability and
legacy Davidson result allocations remain review considerations. No task registry
entry/delegate report exists in this checkout; session/dashboard publication and
README new-result ingestion remain orchestrator-owned. No commit was made.

### Compact retained timing evidence

Each JSON line gives case, budget, baseline/final raw wall samples in seconds;
phase medians and energies are tabulated above. All repetitions within each
benchmark command had identical energy and density hash (enforced by the runner).

```jsonl
{"case": "gamma", "threads": 1, "before_wall_s": [0.112580291, 0.11098308300000001, 0.11055337500000001], "after_wall_s": [0.5364219160000001, 0.10966995800000001, 0.10905683300000002]}
{"case": "k222s", "threads": 1, "before_wall_s": [0.35578316699999996, 0.3585615420000001, 0.364549875], "after_wall_s": [0.34844708300000005, 0.34839395799999995, 0.3487854170000002]}
{"case": "k222", "threads": 1, "before_wall_s": [0.6663425000000003, 0.6645673330000004, 0.6654867499999999], "after_wall_s": [0.6514682500000002, 0.654110583, 0.650691042]}
{"case": "k444s", "threads": 1, "before_wall_s": [2.614540375, 2.605873334, 2.6059610000000006], "after_wall_s": [2.552627375, 2.550915916, 2.5520446670000005]}
{"case": "si8", "threads": 1, "before_wall_s": [3.529003125000001, 3.5418400420000005, 3.5372105000000005], "after_wall_s": [3.7540604999999996, 3.7519456669999993, 3.7804411249999994]}
{"case": "gamma", "threads": 2, "before_wall_s": [0.092941584, 0.092758333, 0.09242975], "after_wall_s": [0.09228108299999999, 0.09170329199999999, 0.09170712499999997]}
{"case": "k222s", "threads": 2, "before_wall_s": [0.206233917, 0.250535375, 0.20613662499999985], "after_wall_s": [0.20303008300000003, 0.20383562499999996, 0.20403445900000006]}
{"case": "k222", "threads": 2, "before_wall_s": [0.36813570900000014, 0.36772291700000004, 0.3682974160000001], "after_wall_s": [0.36402791700000003, 0.362344, 0.3619150000000002]}
{"case": "k444s", "threads": 2, "before_wall_s": [1.3719460000000003, 1.3739057910000003, 1.4282253330000003], "after_wall_s": [1.344825583, 1.3447189169999998, 1.3509109170000002]}
{"case": "si8", "threads": 2, "before_wall_s": [3.1004706669999997, 3.097067417, 3.1027445], "after_wall_s": [2.6142488330000004, 2.6718101250000004, 2.6223192920000002]}
{"case": "gamma", "threads": 4, "before_wall_s": [0.085582042, 0.08480141600000002, 0.08487020799999997], "after_wall_s": [0.083779291, 0.083415083, 0.083755042]}
{"case": "k222s", "threads": 4, "before_wall_s": [0.14842975000000003, 0.14787133400000008, 0.14782833299999998], "after_wall_s": [0.14587529100000002, 0.14594312499999995, 0.146066417]}
{"case": "k222", "threads": 4, "before_wall_s": [0.249241834, 0.250314417, 0.24885725000000014], "after_wall_s": [0.24810479099999994, 0.2494019999999999, 0.24698100000000012]}
{"case": "k444s", "threads": 4, "before_wall_s": [0.8508494579999999, 0.8446654160000002, 0.8491845829999995], "after_wall_s": [0.863383625, 0.838687792, 0.8341287500000005]}
{"case": "si8", "threads": 4, "before_wall_s": [2.9160355410000003, 2.9315147909999997, 2.920329542000001], "after_wall_s": [2.0258897079999993, 2.0115028329999998, 2.028634167]}
{"case": "gamma", "threads": 8, "before_wall_s": [0.08316754199999998, 0.08288324999999999, 0.083088], "after_wall_s": [0.08201462500000001, 0.081721291, 0.08108395900000001]}
{"case": "k222s", "threads": 8, "before_wall_s": [0.126543375, 0.12806995799999998, 0.126352291], "after_wall_s": [0.12455462499999997, 0.12429737499999999, 0.12332058299999993]}
{"case": "k222", "threads": 8, "before_wall_s": [0.17809175, 0.1787114999999999, 0.18030250000000003], "after_wall_s": [0.17939087499999995, 0.17843624999999996, 0.1772950419999999]}
{"case": "k444s", "threads": 8, "before_wall_s": [0.5951217499999999, 0.6032309590000002, 0.6094902909999997], "after_wall_s": [0.596923667, 0.5830745000000002, 0.5882161670000001]}
{"case": "si8", "threads": 8, "before_wall_s": [2.837652167, 2.822592041, 2.839838499999999], "after_wall_s": [1.813366125, 1.8113748750000003, 1.8000496249999998]}
```

Binary SHA-256: final `31532965a6ac416075527b22c0769b90787513c9d929d9f5d228d870d01d3e06`; baseline `8b1eff3f9b7ed69e1422c3c9e2ab587b464aade83dce7efe83bf7c05b3c31da7`.

Final additional feature check: `cargo check -p dft-numerics --features mkl --all-targets --offline` exits 0, including the MKL-configured numerical test code.

## P10J: persistent workers and MKL logical products (2026-11-05)

Pre-change evidence supplied in the ticket: dev08 P10I `494bece`, si8 pure
T=1/4/8/16 wall 6.89/6.20/7.61/10.15 s; 10/22 T=8 samples stop in
`tiled::run -> std::thread::scope -> pthread_clockjoin`. T=16 has 12/29
samples with only one live thread. Diagonalization alone increases beyond T=4.
Source inspection confirms recursive OS thread creation for every dispatched
GEMM and Hamiltonian band block. MKL additionally packs both operands for every
256-element panel and calls sequential BLAS separately for each panel, then
copies/sums outputs. These are measured scheduling evidence and a source-based
MKL regression hypothesis, respectively; Mac cannot establish MKL timings.
The implementation will retain pure-Rust panel arithmetic and use a single
sequential MKL call per logical GEMM (a fixed decomposition at all budgets).
The latter changes historical MKL rounding and requires remote energy validation.

Additional architect-supplied dev08 MKL T=1 profile: tile glue grows from
0.325 to 3.874 s, memcpy/memset from 0.433 to 1.107 s, while MKL kernel
samples fall from 1.137 to 0.664 s. Projection improves from 1.354 to 0.443 s;
blocked projection is therefore retained. Nonlocal samples grow from 0.325
to 0.885 s. `git diff fe0492b 494bece -- crates/dft-physics/src/hamiltonian.rs`
shows only the `inner_threads` accessor was added: NonlocalOperator arithmetic
is identical in those commits. Local si8 application counts rise 2060 -> 2212
(7.38%), insufficient to explain a 2.7x attribution increase. Symbol attribution,
inlining/cache effects, and remote per-call timing remain distinct from a
verified source-level root cause. A cache-local output traversal is being
measured against an independent copy of the old traversal; every output keeps
its original atom/channel arithmetic order and scratch is reserved once.

Implementation: an SCF-owned Rayon pool is shared by k solves, GEMMs,
projection and Hamiltonian band tasks, including nested calls. No global pool
configuration and no unsafe lifetime erasure. The T=1 executor creates zero
workers and does not enter Rayon. A startup counter records actual workers;
it is written only on worker startup, never during submission. Existing ordered
error/panic handling is retained and both branches finish before propagation.
The ticket explicitly offers Rayon; after the architect's continuation direction,
the existing vendored dependency is promoted with the minimum numerics manifest
and lockfile changes. No workspace manifest or new package is introduced.

Pure rotations write fixed 256-row panels directly through disjoint strided
output views. Gram products keep the original 256-plane-wave partial matrices
and ascending accumulation; removing those partials would change arithmetic.
All input operands remain borrowed views. MKL instead uses one sequential
logical product at every budget, without panel scratch, operand packing,
assembly, or per-panel finite scans/anchor reads. This trades intra-GEMM MKL
parallelism for efficient fixed-backend arithmetic; MKL still parallelizes
independent k solves and Hamiltonian band batches. No tolerances change.

### P10J local verification and acceptance limits

Mac macOS 27.0.1 arm64, 12 logical CPUs, Rust 1.99.1. Shared host: these are
local measurements, not quiet-dev08 acceptance. The final benchmark driver ran
commands sequentially, without overlapping its own builds/tests and measurements.
Baseline HEAD `c943279` has no crate-source diff from P10I `494bece`. No QE
source was read. No new QE reference is claimed. Three wall samples per entry;
the same five-case protocol and UPF were used before/after.

Binary SHA-256 before: `1a479abe5c71defa99f52237fe027fe3c8d7299106672dd91108d276a3cf7283`;
final: `42bdc235ccbcad279ab4c45d14110a97495f6d85a4ecb0fbd84f9c2af857ebf7`.
UPF SHA-256: `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.

| Case | T | Before wall s | After wall s | Before/after |
| --- | ---: | ---: | ---: | ---: |
| si8 | 1 | 3.773557417 | 3.744093500 | 1.007869 |
| si8 | 2 | 2.574230791 | 2.489679959 | 1.033961 |
| si8 | 4 | 2.035547625 | 1.975378750 | 1.030459 |
| si8 | 8 | 1.773360875 | 1.857173167 | 0.954871 |
| si8 | 12 | 1.739875000 | 1.737495458 | 1.001370 |
| gamma | 1 | 0.110252667 | 0.110652833 | 0.996384 |
| gamma | 2 | 0.091013959 | 0.091617791 | 0.993409 |
| gamma | 4 | 0.083401333 | 0.084285042 | 0.989515 |
| gamma | 8 | 0.079923750 | 0.082358625 | 0.970436 |
| gamma | 12 | 0.080952083 | 0.083011000 | 0.975197 |
| k222s | 1 | 0.361803458 | 0.352917833 | 1.025178 |
| k222s | 2 | 0.203144417 | 0.203048500 | 1.000472 |
| k222s | 4 | 0.146064875 | 0.149054333 | 0.979944 |
| k222s | 8 | 0.125120958 | 0.126042500 | 0.992689 |
| k222s | 12 | 0.123187625 | 0.123380916 | 0.998433 |
| k222 | 1 | 0.654386750 | 0.658657750 | 0.993516 |
| k222 | 2 | 0.362499959 | 0.366410792 | 0.989327 |
| k222 | 4 | 0.246409333 | 0.246600125 | 0.999226 |
| k222 | 8 | 0.193606334 | 0.178503375 | 1.084609 |
| k222 | 12 | 0.181971541 | 0.180408500 | 1.008664 |
| k444s | 1 | 2.557629000 | 2.577893875 | 0.992139 |
| k444s | 2 | 1.347024666 | 1.371418375 | 0.982213 |
| k444s | 4 | 0.834103291 | 0.830045292 | 1.004889 |
| k444s | 8 | 0.596984125 | 0.590573375 | 1.010855 |
| k444s | 12 | 0.502824334 | 0.519489916 | 0.967919 |

Final si8 T=8 speedup over final T=1 is **2.016x**, below 3x.
The T=8 before/after regression and smaller timing regressions above are not
hidden or accepted as wins. T=1 differences on the five cases range from about
-2.5% to +0.8%; strict no-regression acceptance is not established on this shared
host. k444s T=8 improves, but not every budget improves. The dev08 3x target,
MKL old/new >=0.98 on every protocol combination, and return of si8 MKL T=1
to ~6.3 s all remain **unverified remote gates**.

Intra-k density is deliberately serial in this review version to eliminate the
remaining per-batch scoped thread creation without editing outside the ticket.
It remains parallel across k points on the SCF pool. In si8 T=8, density time
changes from 0.016815 to 0.051618 s (about 0.035 s added); diagonalization changes
from 1.522996 to 1.571365 s. This fallback is not a claim of density acceleration.
The architect can authorize the small `band_density.rs` pool-handle change to
restore intra-k density parallelism. That file has not been edited.

### Numerical equivalence

All five cases, all measured T=1/2/4/8/12: absolute energy delta from P10I is
**0 Ha**, iteration counts and density SHA-256 are identical. Each of the ten
final before/after `bench_scf.py` invocations passed `--require-identical`.
The historical `docs/performance/baseline.json` gamma/k222 comparison also
passes (1.63e-13 / 6.04e-14 Ha); historical density differs as it already did in P10I.

| Case | Final energy Ha | Iterations | Density SHA-256 |
| --- | ---: | ---: | --- |
| gamma | -7.2955304682107647 | 10 | `1854d1ea91d96933ec85e908ec62b4466ddb044abe16b060c1cbe376daea0d3f` |
| k222s | -7.919423412469083 | 10 | `2b310246ec91cf7ea8b845e03fcee0fe6c7769e5556c4af16ab5c8e4916e3ee5` |
| k222 | -7.8312300571211368 | 10 | `f203190f8c973771213a5074eb769c3ceecf264945aa805ea2da6edad5988840` |
| k444s | -7.9259761266480915 | 10 | `ad4d3eabe1b8823365908f6647c4470831a51416e8c738cd3f2aaabf7549b220` |
| si8 | -31.321803601031256 | 15 | `9152c7a9c64db67f623379a343a7b8c84b532ccc08736f8870292ad287e3e4de` |

Six nonphysical synthetic Gamma fixtures also pass direct energy and density.bin
byte comparisons against the old executable at T=1/2/4/8 (24 comparisons).
Cubic cell rows are [[8,0,0],[0,8,0],[0,0,8]] Bohr; skew rows are
[[8,0,0],[1,7,0],[0.5,1,6]] Bohr. For N=4,8,20, atom i=0..N-1 is at
[(i*0.61803398875)%1,(i*0.41421356237)%1,(i*0.73205080757)%1].
Use bench_suite.py synthetic valence=2, channels=[0,1,2], radius=1, shape="soft",
cutoff_ry=3, tier="long", timeout_s=120, deterministic=true, seed=1. The suite
uses a Si label only for input identity validation: these are not silicon data
or realistic large-system performance evidence. Exact inputs and results are
in `/tmp/p10j/general-{before,final}-t*/`; no UPF is committed.

| Fixture | Energy Ha | Iterations | Density SHA-256 |
| --- | ---: | ---: | --- |
| syn-20-cubic | 1208.5691875011141 | 7 | `9548037d81e930094440aa4822869d7ee52ed964179c76167ad180a49eb27c00` |
| syn-20-skew | 1677.6663215065712 | 7 | `4b88a4db6022ecc5375b78ae8a64b8213de14859d2330d23645f00f74966c5bc` |
| syn-4-cubic | -0.31779893646649704 | 9 | `183bdb1387c1ca3b359d1c36d5b2b4bac29fbda171bbc6023955d52754036bc1` |
| syn-4-skew | 2.9423497853369449 | 9 | `a79d14e927db5a366b1351cc7b4d87b37f06c52df6920ecbf453521ce2321d79` |
| syn-8-cubic | 46.806203281047821 | 7 | `8f383670ad613439dcd37f02e3d4c42f425e4dbac010b57da78cd2482a20c162` |
| syn-8-skew | 77.082482692123506 | 7 | `26142b6bcfecafa3a5b7934e4f62f3bb634dcaba635770dfd30b7114de372e3d` |

### Granularity experiment

Gram products of rows x width matrices; five batches of 50 repetitions per
mode, median microseconds per call, same fixed decomposition and ordered sum.
The pool mode forces dispatch to expose the crossover; it is not the production
threshold path. Inputs cover multiple band counts and long/short dimensions.

| Rows | Width | Complex multiply-adds | Inline us | Pool T=8 us |
| ---: | ---: | ---: | ---: | ---: |
| 257 | 8 | 16448 | 14.561 | 20.222 |
| 257 | 16 | 65792 | 31.728 | 32.276 |
| 257 | 32 | 263168 | 88.862 | 68.238 |
| 257 | 64 | 1052672 | 172.687 | 175.646 |
| 1025 | 8 | 65600 | 14.903 | 20.942 |
| 1025 | 16 | 262400 | 52.940 | 31.163 |
| 1025 | 32 | 1049600 | 180.514 | 72.098 |
| 1025 | 64 | 4198400 | 684.451 | 229.558 |
| 4097 | 8 | 262208 | 59.841 | 41.023 |
| 4097 | 16 | 1048832 | 211.115 | 84.810 |
| 4097 | 32 | 4195328 | 719.628 | 219.471 |
| 4097 | 64 | 16781312 | 2726.339 | 609.138 |

Selected conservative dispatch guard: at least 262,144 complex multiply-adds
(2,097,152 real flops) AND at least four panels. The guard depends on shape,
never T. Small 257-row cases can contain one full panel and a one-row tail;
keeping them inline avoids dispatch with essentially no parallel work. This is
a Mac measurement, not a claim of a portable optimum.

### Nonlocal traversal experiment

Five batches of 50 repetitions, median us; same compiler and independent old
traversal in the test. All output bits match, including nonzero initial output.
Weights use reusable scratch. Only traversal changes: atom and channel sums
for each output coefficient stay in the old order.

| Atoms | Rows | Channels/atom | Old us | Candidate us |
| ---: | ---: | ---: | ---: | ---: |
| 4 | 257 | 1 | 4.985 | 5.045 |
| 8 | 4097 | 3 | 197.409 | 117.521 |
| 20 | 16385 | 9 | 2883.383 | 2925.173 |

The 8-atom synthetic shape improves substantially, while the other shapes are
roughly neutral/slightly slower. This is a review candidate, not evidence that
the dev08 2.7x profile attribution has been recovered. Remote measurement and
generality review are still required; no universal nonlocal speedup is claimed.

### Commands and worker lifetime

All final commands in `/tmp/p10j/commands.json` exit 0, including:

```sh
cargo build --release -p dft-cli --offline
cargo test --release -p dft-numerics --lib granularity_table -- --ignored --nocapture
cargo test --release -p dft-physics --lib nonlocal_traversal -- --include-ignored --nocapture
scripts/check.sh
cargo check -p dft-cli --features mkl --all-targets --offline
cargo clippy -p dft-numerics --features mkl --all-targets --offline -- -D warnings
cargo test -p dft-engine --lib scf_reuses_workers -- --nocapture
```

Benchmark invocation template (T=1,2,4,8,12; binaries before-cli and final):

```sh
python3 -I scripts/bench_scf.py --binary target/release/dft-cli \
  --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF \
  --cases gamma,k222s,k222,k444s,si8 --repeats 3 --threads 8 \
  --baseline /tmp/p10j/before-t1.json --require-identical --record /tmp/p10j/final-t8.json
python3 -I scripts/bench_suite.py --binary target/release/dft-cli \
  --case-dir /tmp/p10j/cases --include-long --threads 8 \
  --output /tmp/p10j/general-final-t8.json --work-dir /tmp/p10j/general-final-t8
```

The complete-SCF startup-counter test runs 11 iterations: T=1 creates 0 workers,
T=4 creates 4 workers; counts stay unchanged through completion. The executor
test also reuses four workers over 100 nested regions and recovers after a
worker panic. The SCF production path no longer calls std::thread::scope; the
legacy density helper still contains that code for its existing tests, but SCF
constructs it with one worker. Pure matrix tests include budgets 0/1/2/3/4/8/16,
partial panels, all adjoints, scalar oracle, stable scratch, and rank tests.

Initial MKL Clippy found a needless return (101); corrected and final feature
Clippy passes. An intermediate workspace check found an unused test import (1);
removed and final workspace build/test/fmt/clippy all PASS (0). An initial
benchmark with a relative UPF path failed (3); all retained measurements use
the absolute path. A recursive closure-reference type caused a release build
failure (101); passing a shared reference without adding reference levels
resolved it. These failures are not reported as successful gates.

README figure regeneration and its three specification tests exit 0, with no
generated diff. New reference ingestion/dashboard JSON publication remain
outside this ticket ownership; the authoritative external delegate registry
records this task, while the checkout registry lacks it. CLI threading metadata
still names scoped threads (an out-of-scope hard-coded label); actual execution
is as described here. MKL was type-checked only, never linked/run on this Mac.
No remote server changes, no numerical tolerance/default changes, no commit.

### Retained raw wall samples (seconds)

```jsonl
{"case": "gamma", "threads": 1, "before_wall_s": [0.12100762500000001, 0.11009787500000001, 0.11025266699999997], "after_wall_s": [0.46546025, 0.11054066699999998, 0.11065283299999995]}
{"case": "k222s", "threads": 1, "before_wall_s": [0.36334587500000004, 0.3537758750000001, 0.361803458], "after_wall_s": [0.35171025, 0.352917833, 0.35405395900000003]}
{"case": "k222", "threads": 1, "before_wall_s": [0.6553941250000002, 0.65438675, 0.6533040420000003], "after_wall_s": [0.6604214160000001, 0.658069459, 0.6586577500000002]}
{"case": "k444s", "threads": 1, "before_wall_s": [2.55627625, 2.5576289999999995, 2.5654136670000014], "after_wall_s": [2.5967527079999995, 2.5762810420000006, 2.577893874999999]}
{"case": "si8", "threads": 1, "before_wall_s": [3.773557417000001, 3.7810254590000003, 3.7537350829999987], "after_wall_s": [3.7336026669999995, 3.7440935000000017, 3.7525112500000013]}
{"case": "gamma", "threads": 2, "before_wall_s": [0.091535125, 0.091004417, 0.09101395899999998], "after_wall_s": [0.09197733400000002, 0.09082537499999999, 0.09161779099999998]}
{"case": "k222s", "threads": 2, "before_wall_s": [0.20314441700000002, 0.2029570420000001, 0.20339029200000003], "after_wall_s": [0.20304850000000008, 0.20401195800000005, 0.20259024999999997]}
{"case": "k222", "threads": 2, "before_wall_s": [0.362499959, 0.36239287500000006, 0.3627051250000002], "after_wall_s": [0.36525675000000013, 0.36641079199999993, 0.4352149169999999]}
{"case": "k444s", "threads": 2, "before_wall_s": [1.3466667500000002, 1.3470246659999998, 1.3510660420000002], "after_wall_s": [1.374684916, 1.3691084589999996, 1.3714183750000002]}
{"case": "si8", "threads": 2, "before_wall_s": [2.5501285830000002, 2.608148625, 2.5742307909999997], "after_wall_s": [2.488680208, 2.4914503749999994, 2.489679959]}
{"case": "gamma", "threads": 4, "before_wall_s": [0.0840775, 0.0833045, 0.08340133300000002], "after_wall_s": [0.08475608300000001, 0.084151958, 0.084285042]}
{"case": "k222s", "threads": 4, "before_wall_s": [0.14641983400000003, 0.14606487499999998, 0.145532333], "after_wall_s": [0.146772042, 0.14905433299999993, 0.14939816699999997]}
{"case": "k222", "threads": 4, "before_wall_s": [0.24881366700000007, 0.24460108299999983, 0.24640933300000012], "after_wall_s": [0.24660012500000006, 0.24276020900000006, 0.2511605830000001]}
{"case": "k444s", "threads": 4, "before_wall_s": [0.8341032910000001, 0.8357893750000001, 0.8300177079999997], "after_wall_s": [0.8374634159999998, 0.8274792500000001, 0.8300452920000003]}
{"case": "si8", "threads": 4, "before_wall_s": [2.0441743750000003, 2.0179336250000004, 2.0355476249999995], "after_wall_s": [1.9746299159999996, 1.9773984169999999, 1.9753787500000008]}
{"case": "gamma", "threads": 8, "before_wall_s": [0.07973770799999999, 0.07992375000000002, 0.08147412499999998], "after_wall_s": [0.08294599999999999, 0.08235862499999999, 0.08230633299999998]}
{"case": "k222s", "threads": 8, "before_wall_s": [0.12518666700000003, 0.12512095800000006, 0.12400866700000002], "after_wall_s": [0.1260425, 0.12572470800000002, 0.12660233399999998]}
{"case": "k222", "threads": 8, "before_wall_s": [0.1794037919999999, 0.19360633399999994, 0.204145625], "after_wall_s": [0.17936050000000003, 0.17646687500000002, 0.17850337500000002]}
{"case": "k444s", "threads": 8, "before_wall_s": [0.5972275, 0.593997041, 0.5969841250000001], "after_wall_s": [0.5889834170000001, 0.590573375, 0.5942791669999998]}
{"case": "si8", "threads": 8, "before_wall_s": [1.7653575829999997, 1.7733608749999998, 1.8090862910000007], "after_wall_s": [1.8588701250000002, 1.857173167, 1.8423446660000007]}
{"case": "gamma", "threads": 12, "before_wall_s": [0.08074712499999999, 0.08095208300000001, 0.08251716599999998], "after_wall_s": [0.083281416, 0.082665917, 0.08301099999999997]}
{"case": "k222s", "threads": 12, "before_wall_s": [0.11830174999999998, 0.12318762500000002, 0.132602542], "after_wall_s": [0.12322725000000001, 0.12338091600000001, 0.12354158300000007]}
{"case": "k222", "threads": 12, "before_wall_s": [0.204564458, 0.1819715409999999, 0.17983700000000002], "after_wall_s": [0.17990887499999997, 0.18227279200000002, 0.18040849999999997]}
{"case": "k444s", "threads": 12, "before_wall_s": [0.502824334, 0.5152213749999999, 0.4999750829999998], "after_wall_s": [0.5194899159999999, 0.532255792, 0.5170357910000001]}
{"case": "si8", "threads": 12, "before_wall_s": [1.7380063750000003, 1.7398750000000005, 1.7491301669999997], "after_wall_s": [1.7374954579999997, 1.7331418750000003, 1.737531166000001]}
```

### P10J resumed after architect decisions (2026-11-05)

The preceding local verification is retained as the rejected/intermediate
version, not the final acceptance result. The architect authorized
`band_density.rs` integration and accepted pinned Rayon 1.12.0. Intra-k density
is now parallel on the SCF's existing pool; each band retains its own independent
contribution and the reduction order is unchanged. The density result is now a
reused plan-owned array, avoiding a new partial array per call. SCF folds those
completed arrays in ascending k order. No new pool is constructed by density.

The Hamiltonian output mutex and copy were also removed: recursive tasks receive
disjoint mutable column-major output slices and write their bands directly.
No thread spawn, std thread scope, or Mutex remains in numerics/solvers/physics/
engine source. Scratch/pool allocation occurs at construction. `GemmWorkspace::multiply_into`,
Hamiltonian apply_block and BandDensity::compute allocate no task or numerical
buffer on successful repeated calls in our code; errors may format messages.
Legacy solver result ownership (`multiply`, orthonormalized result matrices,
SCF orbital/result vectors) still allocates as before and is not misrepresented
as an allocation-free complete solver.

Source audit of vendored rayon-core 1.13.0 confirms `join` uses StackJob and
SpinLatch, while ThreadPool::install's same-registry path invokes its closure
inline. An external installation uses a reusable thread-local LockLatch. K-solve
regions are now installed once at their outer boundary, including Gamma, so
repeated GEMM/band tasks stay resident instead of repeatedly entering from the
external caller. No task boxing or global pool configuration is used. Dependency
queue growth/backend-internal allocation is not independently instrumented; our
reused-buffer tests and source audit establish the stated adapter boundary.

An initial seven-pair alternating si8 T=8 comparison after density restoration
still showed 1.800989 s before vs 1.822750 s after (1.21% regression).
Density was restored to 0.017271 s vs 0.017021 s before; remaining diagonalization
was 1.571341 vs 1.550524 s. The macOS sample attempts on both binaries exited
255 (cannot examine process; sandbox/process attachment permission), so no
new profile is claimed. Numerical CLI runs used for sampling both exited 0.

The next controlled scheduling change replaced worker-sized groups of fixed
GEMM tiles with one leaf job per fixed tile. The old scheduling exposed only
T chunks and could not steal the remaining tile from a busy chunk; the new tree
exposes all independent tiles to the existing workers. Arithmetic and summation
are unchanged, and the tree itself no longer depends on T. A five-run pilot
si8 T=8 fell to 1.724437708 s (see exact final measurements below); the final matrix
and alternating comparison, rather than this pilot, determine acceptance.

### Final P10J review version: resident phases and individual tile jobs

The final energy-audit loop is installed once for the entire bounded chunk
sequence, rather than once per chunk. This is the same single-boundary strategy
used for diagonalization and density. A targeted k444s T=2 experiment reduced
audit time from 0.150048 to 0.116903 s and wall time from 1.372 to 1.334 s.
The earlier 1.4% k444s T=2 regression was therefore investigated, not accepted.
An alternative small-column Nonlocal streaming branch did not reduce that audit
time and was removed. The final production kernel retains the measured fixed
256-row output traversal, original atom/channel arithmetic, and reusable weights.

Final production binary SHA-256: `231358500a3e030d4c1cbb8aef2c15d3d6aab4708330ae5b883220aa03b5f56c`.
Baseline remains `1a479abe5c71defa99f52237fe027fe3c8d7299106672dd91108d276a3cf7283`
(P10I source). Host/compiler/UPF and synthetic fixture recipes are unchanged
from the preceding record. Tables below supersede the intermediate timing tables.
Commands were run sequentially without overlapping our own builds and timings;
the Mac is still a shared host, not a quiet-node acceptance environment.

| Case | T | Before median s | Final median s | Before/final |
| --- | ---: | ---: | ---: | ---: |
| si8 | 1 | 3.783201750 | 3.743456375 | 1.010617 |
| si8 | 2 | 2.556416625 | 2.471017500 | 1.034560 |
| si8 | 4 | 2.032091875 | 1.942979208 | 1.045864 |
| si8 | 8 | 1.825680708 | 1.770901167 | 1.030933 |
| si8 | 12 | 1.747488000 | 1.712823875 | 1.020238 |
| gamma | 1 | 0.110966583 | 0.111276208 | 0.997218 |
| gamma | 2 | 0.091486458 | 0.091299750 | 1.002045 |
| gamma | 4 | 0.084447667 | 0.083732167 | 1.008545 |
| gamma | 8 | 0.081330291 | 0.081284250 | 1.000566 |
| gamma | 12 | 0.081662250 | 0.080821584 | 1.010402 |
| k222s | 1 | 0.355575250 | 0.359573042 | 0.988882 |
| k222s | 2 | 0.209448708 | 0.203619625 | 1.028627 |
| k222s | 4 | 0.150472458 | 0.149502375 | 1.006489 |
| k222s | 8 | 0.124470500 | 0.125272125 | 0.993601 |
| k222s | 12 | 0.118810084 | 0.121759208 | 0.975779 |
| k222 | 1 | 0.656463166 | 0.653020375 | 1.005272 |
| k222 | 2 | 0.363295334 | 0.363720209 | 0.998832 |
| k222 | 4 | 0.247602458 | 0.248784209 | 0.995250 |
| k222 | 8 | 0.179291542 | 0.187391292 | 0.956776 |
| k222 | 12 | 0.180120083 | 0.178490250 | 1.009131 |
| k444s | 1 | 2.571535542 | 2.560646750 | 1.004252 |
| k444s | 2 | 1.354912541 | 1.339387792 | 1.011591 |
| k444s | 4 | 0.837692583 | 0.823189125 | 1.017619 |
| k444s | 8 | 0.613585416 | 0.581817250 | 1.054602 |
| k444s | 12 | 0.526887958 | 0.504945333 | 1.043455 |

Independent alternating-order seven-pair si8 T=8 check: before 1.830482375 s,
final 1.761006167 s (medians); 3.796% less wall time. Final is faster
in 7/7 pairs. The original Mac 4.7% regression is no longer present.
Final T=1/T=8 scaling from the three-sample table is 2.114x; this is
reported as a Mac result, not the dev08 >=3x gate.

### Dispatch latency and fixed threshold

Nine batches of 1000 iterations, median us/region, after 100 warmups. Ingress
is empty install; external region is install plus a binary join tree with T
empty leaves; resident region runs that same tree already inside the pool.
Empty leaves may execute on the submitting worker; this is dispatch overhead,
not a forced all-worker wakeup/barrier latency or a cold-pool benchmark.

| T | Empty ingress us | External empty region us | Resident empty region us |
| ---: | ---: | ---: | ---: |
| 1 | 0.002 | 0.004 | 0.003 |
| 2 | 6.111 | 4.365 | 0.017 |
| 4 | 4.213 | 4.309 | 0.035 |
| 8 | 6.687 | 7.172 | 0.266 |
| 12 | 8.452 | 9.753 | 0.699 |
| 16 | 8.601 | 10.155 | 1.632 |

Work threshold remains **262,144 complex multiply-adds = 2,097,152 real flops**
and at least four fixed panels, independent of T. Below it, identical panels
run inline. Above it, each fixed panel is independently stealable. The final
size sweep below uses the same operation and reduction in both modes (five
batches x 50 repetitions, median us); forced pool mode includes external
ingress and does not apply the production scheduling guard.

| Rows | Width | Complex multiply-adds | Inline us | Forced pool T=8 us |
| ---: | ---: | ---: | ---: | ---: |
| 257 | 8 | 16448 | 3.328 | 11.422 |
| 257 | 16 | 65792 | 13.408 | 23.483 |
| 257 | 32 | 263168 | 45.921 | 53.000 |
| 257 | 64 | 1052672 | 174.184 | 180.572 |
| 1025 | 8 | 65600 | 14.906 | 21.476 |
| 1025 | 16 | 262400 | 53.406 | 31.644 |
| 1025 | 32 | 1049600 | 182.070 | 72.266 |
| 1025 | 64 | 4198400 | 692.085 | 243.684 |
| 4097 | 8 | 262208 | 61.002 | 36.858 |
| 4097 | 16 | 1048832 | 214.539 | 72.207 |
| 4097 | 32 | 4195328 | 730.856 | 186.652 |
| 4097 | 64 | 16781312 | 2755.042 | 575.781 |

### Nonlocal evidence and limits

Independent old traversal versus the final cache-blocked traversal; median us
over five batches of 50, exact output-bit assertions. Production source is the
same for both numerical backends. No native MKL execution is involved here.

| Atoms | Rows | Channels | Old us | Final us |
| ---: | ---: | ---: | ---: | ---: |
| 4 | 257 | 1 | 5.926 | 5.047 |
| 2 | 769 | 4 | 21.317 | 21.565 |
| 4 | 1025 | 4 | 39.358 | 35.554 |
| 4 | 4097 | 4 | 109.922 | 81.107 |
| 8 | 2969 | 4 | 93.032 | 94.478 |
| 8 | 4097 | 3 | 96.534 | 97.305 |
| 20 | 16385 | 9 | 2892.298 | 2924.296 |

These timings establish behavior over different projector counts and working
sets, not a universal speedup. Small changes in the standalone benchmark case
list changed compiler-generated benchmark timings while the CLI executable was
byte-identical (verified SHA-256 after reverting the streaming experiment).
Consequently the dev08 Nonlocal 2.7x sample attribution is not declared recovered
from these microbenchmarks. The final candidate improves locality without
changing sums; remote per-call/SCF measurement must decide its platform benefit.

### Final correctness and checks

For gamma/k222s/k222/k444s/si8 at T=1/2/4/8/12/16, all final bench_scf.py
commands exit 0 with --require-identical: energy delta **0 Ha**, same density
SHA-256 and iteration counts as P10I. The six 4/8/20-atom cubic/skew fixtures
also have exact energy and density.bin bytes at T=1/2/4/8 (24 comparisons).
The earlier tables give the identical energies/hashes and fixture recipe.
No tolerance, occupation, band policy, or numerical default changed.

The complete-SCF lifetime test still reports T=1: 0 workers, T=4: 4 workers,
11 iterations each, with no additional starts before completion. Density
checks cover T=1/2/3/4/8/16, stable output and scratch addresses, indexed errors
and reuse after failure. Hamiltonian output now uses safe disjoint slices;
there is no output mutex/copy. No thread::scope, spawn call, or Mutex remains
in the numerics/solvers/physics/engine source. Local joins are stack/spin work
stealing; one external installation/completion occurs per logical phase.

Final verification driver: `python3 /tmp/p10j-final/verify.py` (exit 0). Its
`commands.json` records exact arguments and all exit codes, including workspace
build/test/fmt/clippy, MKL all-targets type check and Clippy, microbenchmarks,
SCF protocol matrix, synthetic suite, and startup-counter test. The MKL code was
not linked or run on this Mac. Final README generation/specification checks
are recorded separately below. Raw measured wall samples are retained here
so the evidence survives removal of temporary logs.

Outstanding remote acceptance: dev08 pure si8 >=3x at T=8; MKL T=1 old/new
>=0.98 for every protocol case and si8 near 6.3 s; within-MKL thread identity
and cross-backend energy; Nonlocal platform benefit. These are unverified, not
claimed passes. The Mac thread-creation mechanism and T=8 regression are fixed.
CLI backend/thread-count metadata still has the pre-existing scoped-thread
label outside ticket ownership; actual engine execution and measured requested
budgets are described above. Dashboard publication/new README result ingestion
remain orchestrator-owned. No commit or remote server change was made.

### Final retained raw wall samples (seconds)

```jsonl
{"case": "gamma", "threads": 1, "before_wall_s": [0.14096162499999998, 0.11096658299999998, 0.11079783400000004], "final_wall_s": [0.28231775, 0.11094362500000005, 0.11127620799999999]}
{"case": "k222s", "threads": 1, "before_wall_s": [0.35557524999999995, 0.35697424999999994, 0.3540107909999999], "final_wall_s": [0.35957304199999995, 0.36537991700000005, 0.350715417]}
{"case": "k222", "threads": 1, "before_wall_s": [0.6564631660000002, 0.6591534169999997, 0.6560864999999998], "final_wall_s": [0.6536765, 0.6524307500000002, 0.6530203749999997]}
{"case": "k444s", "threads": 1, "before_wall_s": [2.5715355420000003, 2.5684963329999997, 2.5733644160000004], "final_wall_s": [2.56064675, 2.5645413329999993, 2.550131792]}
{"case": "si8", "threads": 1, "before_wall_s": [3.7845609999999983, 3.773529667, 3.78320175], "final_wall_s": [3.747937875, 3.743456375000001, 3.742245624999999]}
{"case": "gamma", "threads": 2, "before_wall_s": [0.091118041, 0.09148645799999999, 0.09779854100000002], "final_wall_s": [0.091864667, 0.09115737500000001, 0.09129974999999999]}
{"case": "k222s", "threads": 2, "before_wall_s": [0.24953941699999999, 0.20944870800000004, 0.20374308400000007], "final_wall_s": [0.20307754200000006, 0.20361962499999997, 0.204605417]}
{"case": "k222", "threads": 2, "before_wall_s": [0.363295334, 0.3630570000000002, 0.36481641600000003], "final_wall_s": [0.365686583, 0.363720209, 0.36317850000000007]}
{"case": "k444s", "threads": 2, "before_wall_s": [1.3569943330000003, 1.354912541, 1.3540071669999998], "final_wall_s": [1.4454587080000003, 1.3393877919999997, 1.338081291]}
{"case": "si8", "threads": 2, "before_wall_s": [2.560337874999999, 2.551449917000001, 2.556416624999999], "final_wall_s": [2.5312980420000004, 2.4709661250000003, 2.4710174999999985]}
{"case": "gamma", "threads": 4, "before_wall_s": [0.084202, 0.08461395800000002, 0.08444766700000003], "final_wall_s": [0.084082, 0.083732167, 0.08310945900000002]}
{"case": "k222s", "threads": 4, "before_wall_s": [0.15047245800000003, 0.14688524999999997, 0.150534041], "final_wall_s": [0.146827667, 0.14950237500000002, 0.15311079100000002]}
{"case": "k222", "threads": 4, "before_wall_s": [0.24845675, 0.24760245800000003, 0.245195834], "final_wall_s": [0.24738654199999988, 0.24878420899999987, 0.250806667]}
{"case": "k444s", "threads": 4, "before_wall_s": [0.8377554999999999, 0.8358927500000002, 0.8376925829999999], "final_wall_s": [0.8281713329999998, 0.8215156669999999, 0.8231891249999999]}
{"case": "si8", "threads": 4, "before_wall_s": [2.0320918750000008, 2.036958041000001, 2.018762584000001], "final_wall_s": [1.9434200000000001, 1.9429792080000006, 1.9389556250000002]}
{"case": "gamma", "threads": 8, "before_wall_s": [0.081330291, 0.08090641699999998, 0.08165158399999997], "final_wall_s": [0.08128425, 0.080797167, 0.082610708]}
{"case": "k222s", "threads": 8, "before_wall_s": [0.12447049999999998, 0.12451329200000005, 0.12425345899999996], "final_wall_s": [0.12363383300000003, 0.12527212499999996, 0.1256085]}
{"case": "k222", "threads": 8, "before_wall_s": [0.17997462500000005, 0.17875400000000008, 0.17929154200000008], "final_wall_s": [0.17701608299999994, 0.18739129200000004, 0.20016316599999984]}
{"case": "k444s", "threads": 8, "before_wall_s": [0.5938195, 0.613585416, 0.627432292], "final_wall_s": [0.5818172500000001, 0.589766375, 0.5761040829999997]}
{"case": "si8", "threads": 8, "before_wall_s": [1.8256807080000002, 1.7933995829999994, 1.8287805830000003], "final_wall_s": [1.7709011670000003, 1.752785125, 1.7726685420000008]}
{"case": "gamma", "threads": 12, "before_wall_s": [0.08247112499999998, 0.08166224999999999, 0.08108516699999999], "final_wall_s": [0.080949125, 0.080821584, 0.08040820900000001]}
{"case": "k222s", "threads": 12, "before_wall_s": [0.11881008399999998, 0.11913895900000004, 0.1181693330000001], "final_wall_s": [0.120459708, 0.12175920799999995, 0.12304704199999994]}
{"case": "k222", "threads": 12, "before_wall_s": [0.18035829199999998, 0.18012008300000004, 0.17861574999999985], "final_wall_s": [0.17997275000000001, 0.17849024999999985, 0.17739883300000003]}
{"case": "k444s", "threads": 12, "before_wall_s": [0.5049285000000001, 0.5346220830000001, 0.5268879580000001], "final_wall_s": [0.49927966700000015, 0.504945333, 0.5074201249999999]}
{"case": "si8", "threads": 12, "before_wall_s": [1.7496089159999997, 1.7382919590000006, 1.7474879999999997], "final_wall_s": [1.7128238750000002, 1.7292924579999998, 1.6958403329999996]}
{"case": "si8", "threads": 8, "paired_check": 0, "before": 1.832585249999994, "after": 1.7544765000000098}
{"case": "si8", "threads": 8, "paired_check": 1, "after": 1.7610061669999766, "before": 1.8402260840000224}
{"case": "si8", "threads": 8, "paired_check": 2, "before": 1.845558249999982, "after": 1.7751604579999878}
{"case": "si8", "threads": 8, "paired_check": 3, "after": 1.7394630840000218, "before": 1.8171394170000212}
{"case": "si8", "threads": 8, "paired_check": 4, "before": 1.7937623340000073, "after": 1.7641034170000012}
{"case": "si8", "threads": 8, "paired_check": 5, "after": 1.7448447499999702, "before": 1.8304823749999741}
{"case": "si8", "threads": 8, "paired_check": 6, "before": 1.806771000000026, "after": 1.7612048749999758}
```

Final short-case T=1 cross-check (nine alternating-order pairs, same executable):
| Case | Before median s | Final median s | Change |
| --- | ---: | ---: | ---: |
| gamma | 0.106985708 | 0.107269042 | +0.265% |
| k222s | 0.345119250 | 0.343244750 | -0.543% |

Gamma differs by less than 0.3 ms; these sub-percent short-run changes are
reported rather than claimed as a universal zero-cost guarantee. k222s reverses
the small apparent regression of its three-sample table under alternating order.
Both retain exact P10I energy/density hashes in every run. The requested si8
T=8 paired comparison remains faster in all seven pairs; k444s is faster at
all five measured budgets in the final matrix. Quiet dev08 measurements remain
the performance acceptance source.

Final README/figure regeneration and its three specification tests exit 0;
README.md and docs/img have no generated diff. Final git diff --check exits 0.
All code changes remain uncommitted on task/P10J-PERSISTENT-POOL.

```jsonl
{"case": "gamma", "threads": 1, "pair": 0, "version": "before", "wall_s": 0.10918866599999999}
{"case": "gamma", "threads": 1, "pair": 0, "version": "after", "wall_s": 0.10997512500000003}
{"case": "gamma", "threads": 1, "pair": 1, "version": "after", "wall_s": 0.11191879200000004}
{"case": "gamma", "threads": 1, "pair": 1, "version": "before", "wall_s": 0.109778875}
{"case": "gamma", "threads": 1, "pair": 2, "version": "before", "wall_s": 0.10700658299999993}
{"case": "gamma", "threads": 1, "pair": 2, "version": "after", "wall_s": 0.10801120799999997}
{"case": "gamma", "threads": 1, "pair": 3, "version": "after", "wall_s": 0.10874258299999995}
{"case": "gamma", "threads": 1, "pair": 3, "version": "before", "wall_s": 0.10698570800000007}
{"case": "gamma", "threads": 1, "pair": 4, "version": "before", "wall_s": 0.10656266700000006}
{"case": "gamma", "threads": 1, "pair": 4, "version": "after", "wall_s": 0.10726904199999998}
{"case": "gamma", "threads": 1, "pair": 5, "version": "after", "wall_s": 0.10698300000000005}
{"case": "gamma", "threads": 1, "pair": 5, "version": "before", "wall_s": 0.10677875000000014}
{"case": "gamma", "threads": 1, "pair": 6, "version": "before", "wall_s": 0.10635745800000018}
{"case": "gamma", "threads": 1, "pair": 6, "version": "after", "wall_s": 0.10685845800000005}
{"case": "gamma", "threads": 1, "pair": 7, "version": "after", "wall_s": 0.10630841599999985}
{"case": "gamma", "threads": 1, "pair": 7, "version": "before", "wall_s": 0.10704599999999997}
{"case": "gamma", "threads": 1, "pair": 8, "version": "before", "wall_s": 0.106514}
{"case": "gamma", "threads": 1, "pair": 8, "version": "after", "wall_s": 0.10653591600000012}
{"case": "k222s", "threads": 1, "pair": 0, "version": "before", "wall_s": 0.3420332500000003}
{"case": "k222s", "threads": 1, "pair": 0, "version": "after", "wall_s": 0.34176075000000017}
{"case": "k222s", "threads": 1, "pair": 1, "version": "after", "wall_s": 0.341950416}
{"case": "k222s", "threads": 1, "pair": 1, "version": "before", "wall_s": 0.3432762079999998}
{"case": "k222s", "threads": 1, "pair": 2, "version": "before", "wall_s": 0.34178700000000006}
{"case": "k222s", "threads": 1, "pair": 2, "version": "after", "wall_s": 0.3399159169999999}
{"case": "k222s", "threads": 1, "pair": 3, "version": "after", "wall_s": 0.3414910830000002}
{"case": "k222s", "threads": 1, "pair": 3, "version": "before", "wall_s": 0.3448559170000003}
{"case": "k222s", "threads": 1, "pair": 4, "version": "before", "wall_s": 0.34517016599999994}
{"case": "k222s", "threads": 1, "pair": 4, "version": "after", "wall_s": 0.3432447500000002}
{"case": "k222s", "threads": 1, "pair": 5, "version": "after", "wall_s": 0.3476579590000002}
{"case": "k222s", "threads": 1, "pair": 5, "version": "before", "wall_s": 0.34901445900000017}
{"case": "k222s", "threads": 1, "pair": 6, "version": "before", "wall_s": 0.34693754200000004}
{"case": "k222s", "threads": 1, "pair": 6, "version": "after", "wall_s": 0.35830987500000067}
{"case": "k222s", "threads": 1, "pair": 7, "version": "after", "wall_s": 0.3470639999999996}
{"case": "k222s", "threads": 1, "pair": 7, "version": "before", "wall_s": 0.346961125}
{"case": "k222s", "threads": 1, "pair": 8, "version": "before", "wall_s": 0.34511924999999977}
{"case": "k222s", "threads": 1, "pair": 8, "version": "after", "wall_s": 0.3488626669999997}
```

## P10K investigation (2026-11-05)

Starting source: `0cdcf41` (P10J implementation merged at `7224a0f`).
The ticket supplies the dev08 sample: pure si8 T=8 main-thread latch wait
32%, multiply_into 18.8%, subtract_projection 13.8%; MKL gemm_into self
1.299 s versus old gemm 0.323 s, and NonlocalOperator::add 0.487 versus
0.323 s. These are supplied remote measurements, not new local measurements.

Source inspection before changes rejects repeated **per-tile** validation:
pure inputs are checked once before scheduling and output once after reduction;
MKL already performs one borrowed, sequential BLAS product. Matrix::new checks
only dimensions. The three finite scans still short-circuit per complex value,
which can inhibit SIMD; an equivalent Boolean reduction is worth measuring.
There is no MKL tile assembly to remove. Projection and legacy Davidson still
copy workspace/result vectors. The fixed-panel recursion passes its full budget
to both children (intentionally individually stealable panels in P10J); a bounded
leaf experiment must be compared, not presumed faster. P10J already installs
whole solve/audit phases, so nested install is not external pool ingress.

The Nonlocal diff against fe0492b replaces channel-major streaming with fixed
256-row blocks and precomputed weights; its per-element arithmetic order is
unchanged. Existing independent legacy traversal tests cover this comparison.
No MKL runtime is available on this Apple Silicon host; attributing the extra
SCF iteration to blocked projection requires a remote controlled ablation.

New dev08 evidence supplied during this run: 8,044 main-thread FUTEX_WAKE calls,
1,686 blocking waits, 15 SCF / 290 Davidson iterations / 2,070 operator
applications. Of sampled latch waits: Hamiltonian 34.8%, multiply 25.0%,
ordered_chunks 15.2%, construct 13.0%, projection 10.9%. The supplied 184 us
per wake and 879 us per blocking wait describe attributed wait time, not the
wall duration of every logical region. Futex wakes and pool installations are
separate observables and are counted separately in the local experiment.

**Correction to the preceding inspection and historical P10J claim:**
`ScfProblem::step` has a separate serial-k branch selected by
`min(threads, states.len()) == 1`. Gamma takes that branch even at T>1, bypassing
`ordered_solves` and its resident pool installation. Therefore its Hamiltonian
and GEMM phases really do enter the pool externally. Select that branch only
for a total budget <=1; a threaded Gamma solve now uses the existing ordered
path with exclusive mutable KState ownership. No Sync requirement is added to
the operator, no unsafe reference transfer, and no thread budget is enlarged.
This directly addresses the supplied call-site evidence. Construction already
uses one outer region for all k states; the master waits for actual setup work,
not one dispatch per atom/shell.

### Local hypotheses and rejected prototypes

Host: Apple Silicon arm64, rustc 1.99.1 (58cf8f9c3 2026-10-13), release,
shared Mac. `sample <own-child-pid> 1 1 -file ...` returned 255 (sandbox
cannot examine the process); no local sampling profile is claimed. Initial
`scripts/check.sh` and `cargo check --features mkl --all-targets` both exited 0.

- Finite Boolean reduction microbenchmark (median of nine batches of 100,
  microseconds, old early-exit / reduction): 1024 elements 0.328/0.321;
  16384 5.393/5.362; 131072 42.963/42.902; 524288 176.470/171.810.
  This is at most 2.6% here and does not explain the remote 0.976 s gap.
  More useful: projection validates immutable q at entry; each of two passes
  formerly rescanned q twice. Private finite-input entry points reuse that
  validation, scan each newly normalized block once, and retain every product's
  output/shape checks. Coefficients checked on output need no immediate rescan.
  Public GEMM entry points still reject nonfinite inputs. No unsafe code added.
- Rotation products for basis and H*basis share a fixed-panel task tree. Each
  panel computes the same separate products, with no fused arithmetic or changed
  reduction. MKL and T=1 use the separate-call path. No GEMM packing introduced.
- Limiting recursive leaves to min(T, panels) was rejected: combined prototype
  si8 T=8 wall median 1.815 s versus baseline 1.727 s (roughly 5% slower).
  Independent fixed panels remain stealable. Bounding pieces to workers worsens
  tails; it does not guarantee that Rayon's other workers remain asleep.
- Existing Nonlocal old-streaming/blocked benchmark (us; independent bit checks):
  (atoms,rows,channels)=(4,257,1): 3.626/3.681; (2,769,4): 21.289/17.229;
  (4,1025,4): 36.220/36.519; (4,4097,4): 100.875/76.404;
  (8,2969,4): 93.159/94.031; (8,4097,3): 96.598/97.767;
  (20,16385,9): 2895.517/2902.473. No general Mac regression justifies
  restoring streaming globally. Leave production Nonlocal unchanged pending
  dev08 measurement, where the supplied 0.164 s loss remains unexplained.
- No spin policy change: supplied sched_yield counts and tail waits do not
  justify burning additional cores. Pool ingress reduction is addressed first.

### Mac production comparison

Three runs per case/budget, median wall seconds. Uninstrumented executables;
no builds or other task benchmarks overlapped these runs. Shared-host timing
noise remains. Before is HEAD/P10J; after includes validation reuse, paired
rotations and the single-k resident branch. No arithmetic/default changes.

| Case | T | Before s | After s | Before/after |
| --- | ---: | ---: | ---: | ---: |
| si8 | 1 | 3.774887584 | 3.666024542 | 1.029695 |
| si8 | 2 | 2.494354667 | 2.445633250 | 1.019922 |
| si8 | 4 | 1.935665875 | 1.844976292 | 1.049155 |
| si8 | 8 | 1.726885541 | 1.647487750 | 1.048193 |
| si8 | 12 | 1.706870500 | 1.590182250 | 1.073380 |
| gamma | 1 | 0.110782209 | 0.109270167 | 1.013838 |
| gamma | 2 | 0.090828583 | 0.090389500 | 1.004858 |
| gamma | 4 | 0.083097250 | 0.083152708 | 0.999333 |
| gamma | 8 | 0.080732625 | 0.080689250 | 1.000538 |
| gamma | 12 | 0.086510667 | 0.080824791 | 1.070348 |
| k222s | 1 | 0.351281167 | 0.348871250 | 1.006908 |
| k222s | 2 | 0.201997375 | 0.201806000 | 1.000948 |
| k222s | 4 | 0.150384750 | 0.146936125 | 1.023470 |
| k222s | 8 | 0.124678584 | 0.124057000 | 1.005010 |
| k222s | 12 | 0.121810583 | 0.119047459 | 1.023210 |
| k222 | 1 | 0.652517417 | 0.650152500 | 1.003637 |
| k222 | 2 | 0.364066875 | 0.361900584 | 1.005986 |
| k222 | 4 | 0.246716750 | 0.241078083 | 1.023389 |
| k222 | 8 | 0.178062375 | 0.177101166 | 1.005427 |
| k222 | 12 | 0.176486625 | 0.173612583 | 1.016554 |
| k444s | 1 | 2.551488417 | 2.535976833 | 1.006117 |
| k444s | 2 | 1.337836500 | 1.352080875 | 0.989465 |
| k444s | 4 | 0.825829583 | 0.813585125 | 1.015050 |
| k444s | 8 | 0.626953042 | 0.569873500 | 1.100162 |
| k444s | 12 | 0.513669833 | 0.495419625 | 1.036838 |

All 25 after case/budget comparisons to the independently executed before T=1
record pass `bench_scf.py --require-identical`: energy difference 0 Ha, identical
density SHA-256 and iteration counts. Pure si8 remains 15 SCF iterations.
The after T=1/T=8 wall speedup is **2.225x**, below 3x on this Mac.

| Case | Energy Ha | Iterations | Density SHA-256 |
| --- | ---: | ---: | --- |
| gamma | -7.2955304682107647 | 10 | `1854d1ea91d96933ec85e908ec62b4466ddb044abe16b060c1cbe376daea0d3f` |
| k222s | -7.919423412469083 | 10 | `2b310246ec91cf7ea8b845e03fcee0fe6c7769e5556c4af16ab5c8e4916e3ee5` |
| k222 | -7.8312300571211368 | 10 | `f203190f8c973771213a5074eb769c3ceecf264945aa805ea2da6edad5988840` |
| k444s | -7.9259761266480915 | 10 | `ad4d3eabe1b8823365908f6647c4470831a51416e8c738cd3f2aaabf7549b220` |
| si8 | -31.321803601031256 | 15 | `9152c7a9c64db67f623379a343a7b8c84b532ccc08736f8870292ad287e3e4de` |

The initial k444s T=2 table has a 1.1% slowdown. Seven alternating-order
pairs resolve this as noise here: before median 1.355740417 s, after 1.346649333 s
(1.006751x). One after outlier was 1.413047459 s; it was retained.

### Controlled change ablations (Mac)

Three repetitions, wall median s, executed sequentially. These small
differences are not additive causal attributions or predictions for dev08.

| Variant | T=1 s | T=8 s |
| --- | ---: | ---: |
| before | 3.771342416 | 1.708182542 |
| validation | 3.724640583 | 1.681641916 |
| resident | 3.765390209 | 1.765678625 |
| both | 3.697109291 | 1.677927750 |
| final | 3.704783375 | 1.669834500 |

validation = finite scan/reuse only; resident = only the SCF branch fix;
both = those two; final adds paired rotations. Resident alone does not show
a Mac wall-time win in this batch despite removing external ingress.
Validation alone saves about 0.047 s at T=1 / 0.027 s at T=8; pairing on top
of both saves about 0.008 s at T=8, within shared-host noise. No hypothesis
is claimed to explain a measured number of remote seconds without rerunning
dev08. All ablations retained the same energies and density hashes.

### Generality

The six existing P10J synthetic recipes above (4/8/20 atoms, cubic/skew,
channels [0,1,2]) were rerun using the before binary at T=1 and after at
T=1/2/4/8. All 24 comparisons have exact energy, iteration count and density
SHA-256 equality. These small-basis, nonphysical fixtures establish identity,
not broad realistic-system scaling. Large/partial-panel and adjoint-shape
numerics tests additionally cover different widths and thread budgets.
The suite comparison at T=2 initially returned 5 because its input fingerprint
includes threads; all six runs themselves passed. Cross-thread equality was
then checked explicitly; T=4/8 suite commands exit 0 without that incompatible
input comparison. No tolerance or fingerprint rule was relaxed.

### Test-counter region measurements

Temporary instrumented release builds of HEAD and the final source counted
`WorkerPool::install` at entry/exit and `WorkerPool::join`. Each pool owned four
AtomicU64 values: external count/ns and resident count/ns (Relaxed); residence
was `rayon::current_thread_index().is_some()` before installation. One global
Relaxed counter counted join calls. A temporary RAII probe around
`Davidson::solve_policy` captured differences for each solve. Counters printed
only at probe/pool destruction. Timers and atomics affect timings; these are
diagnostic builds, not the production table. All instrumentation was removed.
The exact temporary patch generator and logs are `/tmp/p10k_instrument.py`,
`/tmp/p10k-counts.log`, `/tmp/p10k-solver-{before,after}.log`; production source
contains neither P10K_METRICS nor P10K_SOLVER.

The current Mac case has **340 main + 68 certificate = 408 Davidson iterations**,
16 solve calls and 15 SCF iterations, unchanged before/after. This differs from
the supplied dev08 290-iteration workload; do not mix their denominators.

| Solver-only quantity | Before | After |
| --- | ---: | ---: |
| Installed GEMM/band regions | 2463.000000 | 2055.000000 |
| Regions / Davidson iteration | 6.036765 | 5.036765 |
| Regions / SCF iteration | 164.200000 | 137.000000 |
| Mean installed region duration, us | 294.187353 | 327.391934 |
| External ingress inside solve | 2463.000000 | 0.000000 |
| Join calls | 25313.000000 | 20825.000000 |

The fused rotation saves exactly one region per Davidson iteration (408 total).
The outer SCF adds 16 resident whole-solve boundaries instead of 2,463
external entries from inside those solves. Whole-process counts (separate
instrumented execution) are 2,512 external + 16 resident installs before,
65 external + 2,071 resident after; joins 25,425 -> 20,937. The 65 final
external regions include setup, potential updates, whole solves, density and
energy audits. Setup has one installation covering all k states; Nonlocal
construction and worker FFT-plan allocation run inside it, with no pool
installation per atom or shell. The supplied 0.19 s setup latch attribution
therefore includes real setup computation, not 0.19 s empty dispatch.

Whole-process external mean duration: 323.505 us before versus 24,788.548 us
after; the increase is expected because a region now includes a complete solve.
The inline T=1 path has zero join calls and zero parallel regions (49 no-pool
install calls simply execute closures). These counts are logical regions, not
futex syscalls. Rayons internal stealing/wake policy may issue multiple wakes
per tree. No syscall count or worker-wakeup count is inferred from them.

Dispatch microbenchmark, median of nine batches of 1000 after 100 warmups,
unchanged WorkerPool implementation; empty tasks need not wake every worker:

| T | Empty ingress us | External tree us | Resident tree us |
| ---: | ---: | ---: | ---: |
| 1 | 0.000 | 0.001 | 0.001 |
| 2 | 3.779 | 3.683 | 0.011 |
| 4 | 4.160 | 4.193 | 0.030 |
| 8 | 6.827 | 7.291 | 0.188 |
| 12 | 8.529 | 9.985 | 0.674 |
| 16 | 9.119 | 10.709 | 1.146 |

An extra master spin was not tested or adopted. Rayon offers no per-install
worker-parking cap here; the number of useful fixed pieces limits work, but
this is not a guarantee that other pool workers never wake. The rejected
bounded-leaf prototype and lack of a universal wakeup cap remain explicit.

### MKL iteration count and outstanding remote checks

No native MKL is linked/executed on Apple Silicon. The supplied dev08 result
**16 SCF iterations versus the historical 15** is not resolved or hidden.
These changes do not change MKL's BLAS decomposition or projection arithmetic;
validation reuse only removes read-only scans, and paired rotation falls back
to the original two calls for MKL. Therefore this patch is not expected to
restore 15 iterations by itself. The P10I blocked projection and P10J logical
MKL Gram product both changed rounding relative to fe0492b; the observation
alone cannot identify which caused a convergence-boundary crossing.

The remote session should retain the exact same UPF SHA, input, BLAS thread
settings and binary provenance and compare P10J versus this patch first:
energies within 1e-12 Ha, within-backend density identity across T, full residual
history and fixed-point certificate, pure T=8 >=3x, MKL T=1 <=6.5 s. To attribute
16 vs 15, use separate temporary ablations of established-basis projection
(legacy two-pass MGS versus blocked projection) and Gram decomposition; do not
change both at once, tolerances, defaults or mixing. Compare the first differing
residual/history point and stopping conditions; do not accept a historical
iteration count as a correctness criterion. Any ablation that changes arithmetic
must still meet the energy/certificate gates. Do not merge those ablations as
part of this patch. The 0.976 s GEMM and 0.164 s Nonlocal remote gaps cannot be
quantitatively declared recovered from Mac timing or static reasoning.

### Reproduction and gate record

Production benchmark command (T=1,2,4,8,12; before/final executables saved under
`/tmp/p10k-{before,final}`), absolute UPF path required because inputs are created
in a temporary directory:

```sh
python3 -I scripts/bench_scf.py --binary /tmp/p10k-final \
  --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF \
  --cases gamma,k222s,k222,k444s,si8 --repeats 3 --threads 8 \
  --baseline /tmp/p10k-before-1.json --require-identical \
  --record /tmp/p10k-final-8.json
```

All five final commands exited 0. The early trial using a relative UPF path
exited 3; it was corrected to the absolute path, not counted as passing.
`scripts/check.sh`: build/test/fmt/clippy PASS, exit 0 after final code changes.
`cargo check --features mkl --all-targets`: exit 0 (type checking only).
`python3 -I scripts/make_readme_figures.py`: exit 0, generated README/figures
unchanged because existing reference evidence is unchanged.
`python3 -I scripts/tests/test_readme_figures.py`: exit 0, 3 tests passed;
pre-existing Python unclosed-file ResourceWarnings remain.

The task's owned-file list excludes dashboard state, reference datasets and
figure-generator changes. No new scientific/QE reference numbers were invented
or added to those files. The authoritative parent run registry (outside this
worktree) records P10K#2 running; the checked-in worktree registry predates this
run. Dashboard upload/registry completion remains with the orchestrator. No
remote host was changed; no commit was created.

Executable/input provenance SHA-256:

- `/tmp/p10k-before`: `2b5bdc7b5f28bb421aa832eb3b29962c02c27a96d3e3c13f9c32028f5f67ae51`
- `/tmp/p10k-final`: `fe04f720a3fd7f46c1f26f7192d09d5b12e872158519fedc54ce7e3fcf221c47`
- `/Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF`: `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`

Raw production wall-time samples (s), retained for audit:

```json
{"case": "gamma", "threads": 1, "before_wall_s": [0.48323425, 0.11078220900000002, 0.11020425], "after_wall_s": [0.49267112499999993, 0.10927016699999992, 0.10913333300000005]}
{"case": "k222s", "threads": 1, "before_wall_s": [0.3515085, 0.351281167, 0.35116141700000014], "after_wall_s": [0.3488712500000001, 0.3489778749999999, 0.3488233329999999]}
{"case": "k222", "threads": 1, "before_wall_s": [0.6525699580000002, 0.6519921660000003, 0.6525174169999999], "after_wall_s": [0.6509671669999999, 0.6483848329999997, 0.6501524999999999]}
{"case": "k444s", "threads": 1, "before_wall_s": [2.551763541, 2.5501540840000008, 2.551488417], "after_wall_s": [2.535976833, 2.537482292, 2.5303685410000014]}
{"case": "si8", "threads": 1, "before_wall_s": [3.774887584, 3.753132375, 3.7824475], "after_wall_s": [3.671732125, 3.6660245420000006, 3.650778457999998]}
{"case": "gamma", "threads": 2, "before_wall_s": [0.091345041, 0.09082858300000002, 0.09080170800000001], "after_wall_s": [0.09080120800000002, 0.09038949999999998, 0.09038354099999998]}
{"case": "k222s", "threads": 2, "before_wall_s": [0.20247670799999995, 0.201689958, 0.20199737500000003], "after_wall_s": [0.20124779200000004, 0.20323629199999993, 0.20180599999999993]}
{"case": "k222", "threads": 2, "before_wall_s": [0.3684926249999999, 0.361383958, 0.364066875], "after_wall_s": [0.3604152920000001, 0.36193633300000005, 0.36190058400000025]}
{"case": "k444s", "threads": 2, "before_wall_s": [1.3504679579999999, 1.3378364999999999, 1.3369852909999995], "after_wall_s": [1.5056060839999996, 1.352080875, 1.3323814160000005]}
{"case": "si8", "threads": 2, "before_wall_s": [2.552047375000001, 2.473412167000001, 2.4943546669999996], "after_wall_s": [2.47867675, 2.417811624999999, 2.44563325]}
{"case": "gamma", "threads": 4, "before_wall_s": [0.08368479199999998, 0.08309725000000001, 0.082921875], "after_wall_s": [0.08346283400000001, 0.08315270800000002, 0.082957417]}
{"case": "k222s", "threads": 4, "before_wall_s": [0.15123837499999998, 0.15038474999999996, 0.14942408299999999], "after_wall_s": [0.145409458, 0.15383629100000001, 0.14693612499999997]}
{"case": "k222", "threads": 4, "before_wall_s": [0.253353959, 0.24671675000000004, 0.245226041], "after_wall_s": [0.2410780829999999, 0.23919241599999985, 0.24155266600000003]}
{"case": "k444s", "threads": 4, "before_wall_s": [0.8258295829999998, 0.8303068329999999, 0.823436708], "after_wall_s": [0.8217623329999999, 0.8107909580000001, 0.8135851249999999]}
{"case": "si8", "threads": 4, "before_wall_s": [1.933264958, 1.9366081669999993, 1.9356658750000006], "after_wall_s": [1.8864376249999997, 1.8448487919999996, 1.8449762919999992]}
{"case": "gamma", "threads": 8, "before_wall_s": [0.08072358300000002, 0.080746542, 0.08073262499999997], "after_wall_s": [0.08068924999999999, 0.080690292, 0.07958854199999998]}
{"case": "k222s", "threads": 8, "before_wall_s": [0.12476654200000004, 0.12333475000000005, 0.12467858399999998], "after_wall_s": [0.12326958300000002, 0.12405700000000003, 0.13051558399999996]}
{"case": "k222", "threads": 8, "before_wall_s": [0.17821545800000005, 0.17718687500000008, 0.17806237499999988], "after_wall_s": [0.19541254200000002, 0.17464899999999994, 0.17710116599999992]}
{"case": "k444s", "threads": 8, "before_wall_s": [0.6397413750000001, 0.626953042, 0.590805875], "after_wall_s": [0.5789400409999998, 0.5636456669999999, 0.5698734999999999]}
{"case": "si8", "threads": 8, "before_wall_s": [1.7408905419999998, 1.7268855409999997, 1.7217335829999998], "after_wall_s": [1.6583395410000001, 1.6474877500000007, 1.6364545000000001]}
{"case": "gamma", "threads": 12, "before_wall_s": [0.08656325000000001, 0.08651066699999999, 0.08498775000000003], "after_wall_s": [0.08082479099999999, 0.08026370900000002, 0.08179245799999998]}
{"case": "k222s", "threads": 12, "before_wall_s": [0.12209637500000003, 0.12181058300000003, 0.12180170800000001], "after_wall_s": [0.11904745899999997, 0.11885908299999998, 0.11933204100000006]}
{"case": "k222", "threads": 12, "before_wall_s": [0.176486625, 0.17621529199999997, 0.1797105000000001], "after_wall_s": [0.17261666600000003, 0.17361258299999993, 0.17600962500000006]}
{"case": "k444s", "threads": 12, "before_wall_s": [0.513669833, 0.5118042500000002, 0.5208842090000001], "after_wall_s": [0.49366829199999995, 0.5021435840000001, 0.4954196249999998]}
{"case": "si8", "threads": 12, "before_wall_s": [1.7043947499999996, 1.7088302079999993, 1.7068705], "after_wall_s": [1.579857916, 1.599937208, 1.5901822499999998]}
```

Final additional checks: `cargo test --release -p dft-numerics --lib` exit 0
(28 passed, 5 ignored); `cargo test --release -p dft-engine --test
stepwise_equivalence -- --include-ignored --nocapture --test-threads=1` exit 0
(7 passed, including the execution-mode replay matrix). The final uninstrumented
`target/release/dft-cli` is byte-identical to the measured `/tmp/p10k-final`.
A fresh historical baseline.json gamma/k222 comparison exited 0: absolute
energy differences 1.63e-13 / 6.04e-14 Ha. Historical density differs, as it
already did in P10J; this is separate from exact P10J density preservation.
`git diff --check` exits 0. Final source changes are restricted to the ticket's
owned files. All changes remain uncommitted for architect review.

## P10L initial measurement (2026-11-05; before optimization)

Architect-provided dev08 evidence: si8 pure T=8, 12 merged eu-stack runs,
293 samples, code equal to main. Outside-pool main-thread time 6.5%
(0.312/4.805 s), worker utilization 19.5%, busy-worker distribution 0/1/2/3/4
= 9.2/39.6/38.2/11.9/1.0%, no samples with >=5 busy workers. Per-run
max/mean 1.76--3.20. Supplied ideal estimate 1.27 s (5.04x), target utilization
33% for 2x and 55% for 3x. These are supplied remote observations, not runs
performed by this implementer. Pool residency includes serial Davidson work;
it must not be equated with the parallelizable fraction in Amdahl's law.

Source-only `worker-analysis` cfg was enabled without changing numerical
operations. Mac measurements below are single instrumented runs, not quiet-node
acceptance or uninstrumented performance claims. Build: `RUSTFLAGS='--cfg
feature="worker-analysis"' cargo build --release -p dft-cli`. Raw traces, exact
inputs, results and density files: `/tmp/p10l-before/{si8,skew4}-{T}/`;
summary `/tmp/p10l-before/summary.json`. The skew4 fixture doubles the first
primitive silicon lattice vector, with fractional sites (0,0,0), (.125,.25,.25),
(.5,0,0), (.625,.25,.25), same 30 Ry cutoff and UPF as si8.

| Case | T | SCF s | Logical kernel regions s | Measured leaf busy / (T * SCF) |
|---|---:|---:|---:|---:|
| si8 | 1 | 3.491367417 | 2.795968472 | .792 |
| si8 | 2 | 2.2463 | 1.5640 | .637 |
| si8 | 4 | 1.7026 | 1.0313 | .470 |
| si8 | 8 | 1.511036166 | .834697366 | .304 |
| si8 | 12 | 1.4375 | .7595 | .236 |
| skew4 | 1 | .301098625 | .220706908 | .728 |
| skew4 | 2 | .2006 | .1223 | .560 |
| skew4 | 4 | .1641 | .0888 | .389 |
| skew4 | 8 | .147450542 | .073735711 | .247 |
| skew4 | 12 | .1408 | .0671 | .193 |

At si8 T=8: 1192 Gram regions each contain 12 reduction panels, 798 single
rotations and 408 paired rotations each contain 12 panels; 410 Hamiltonian
regions have 1--8 static groups. Aggregate per-phase max/mean worker times:
Gram 1.21, rotation 1.41, paired rotation 1.03, Hamiltonian 1.29. These aggregate
ratios can hide severe individual-region imbalance and are NOT an acceptance
claim of per-region imbalance <1.5. Projection consumes .3051 s (inclusive),
EVD .1417 s, concatenate .0360 s, first_columns .0221 s, residual norms .0399 s,
preconditioner .0137 s. Pool installs cover 1.4939 s but logical kernel regions
cover only .8347 s: residence in the pool does not imply useful parallel work.
Outside logical kernels at T=1 is .6954 s (19.92%); treating that as fixed
serial work gives ideal T=8 time 1.0449 s and speedup 3.34x. This is only a
Mac decomposition: not evidence that dev08 cannot reach 3x.

Inspection confirms an internal two-pass MGS chain, even though projection
against the established basis already uses blocked GEMM. First-pass projections
onto each newly finalized vector can be distributed over future columns while
preserving each column's original arithmetic order; its second pass and rank
normalization remain dependent. Trial this before changing numerical methods.

### P10L experiments and retained candidate

The source retains only arithmetic-preserving scheduling and copy changes:

- First-pass internal MGS work is right-looking for tall blocks (>=2048 rows),
  with fixed one-column pieces. Each column retains the original scalar dot
  grouping and ascending projection order. For >=16384 future complex elements,
  the next column's dependent second pass runs alongside first-pass work on
  later columns. The next normalized basis vector is appended only after both
  branches finish. Small blocks keep the original two-pass loop; changing their
  traversal provided no reliable benefit in the small-cell measurements.
- Gamma/single-inner-pool Hamiltonian calls with more bands than workers use
  one dynamically stealable piece per band, with preallocated scratch slots.
  Scratch acquisition precedes measured piece time. Outer-k calls retain static
  disjoint groups; the first modulo-scratch prototype regressed k222s by causing
  collisions when the inner budget was smaller than the whole pool.
- Concatenation reserves the combined capacity before copying either operand.
  The profile measured this copy phase (.036 s at si8 T=8); this removes an
  avoidable growth/reallocation without changing the output or reduction order.

No thread, tolerance, band, mixing, iteration-limit or numerical-method default
changed. No new dependency or unsafe code was introduced. The dimensional guards
are scheduling choices; small and large/partial dimensions, deficient columns,
zero columns, duplicate columns, and 1e-200/1e200 scales are tested against the
original scalar MGS bit patterns. Generality evidence also includes a twenty-atom
skew/slab-like fixture with local-only and multi-projector potentials below.

Rejected experiments (raw inputs/results/traces retained under `/tmp/p10l-*`):

1. Four-column Gram subtiles and 64-row rotations (`fine`, `band`): retained
   256-term Gram reductions, but faer's smaller calls changed rounding. All five
   benchmark energies stayed within 1e-12 Ha and densities were identical across
   T. However si8 main/certificate Davidson counts changed from 340+68=408 to
   293+69=362 (SCF count stayed 15). The apparent 1.659 -> 1.518 s T=8 improvement
   therefore mixes scheduling with a different workload. T=1 regressed from
   3.679 to 3.820 s. This is not a justified replacement for the P10K arithmetic.
2. Eight-column subtiles (`wide`): si8 retained 408 Davidson iterations but
   instrumented T=8 improvement was only about 2%, with T=1/small-cell losses.
3. Finer decomposition only above 2048 rows (`guard`): small-cell arithmetic
   became identical again, but si8 still changed iteration count. Rejected for
   the same attribution/generality reason; the final code has **no** finer GEMM
   arithmetic and keeps P10K's 256-element panels.
4. Lower GEMM scheduling cutoff 262144 -> 65536 complex products (`cutoff`):
   exact P10K results, but gains/losses varied with T and shape; cutoff unchanged.
5. Unconditional dynamic Hamiltonian pieces under outer-k parallelism: rejected
   in favor of static groups for that case, as described above.

The retained pipeline preserves all 408 Davidson iterations, 15 SCF iterations,
and every benchmark energy/density bit. The exact three-arm rotating-order
measurement (`/tmp/p10l-paired.json`, five repeats at T=8, three at T=1) separates
this result from the discarded smaller-matrix rounding effects:

| Case | T | P10K wall median s | First-pass/band/copy changes s | With dependent-pass overlap s |
|---|---:|---:|---:|---:|
| si8 | 8 | 1.663557583 | 1.622848667 | 1.579725125 |
| k444s | 8 | .577412750 | .577696000 | .571007958 |
| si8 | 1 | 3.687254375 | 3.684710584 | 3.667978458 |

This is a 5.0% si8 T=8 wall reduction (1.053x relative improvement), not a 2x
improvement over P10K. T=1/T=8 speedup for the retained candidate is 2.322x on
this Mac; the supplied dev08 baseline is different and remains to be measured.
The paired k444s result shows no regression; shared-host timings remain indicative.

### Final uninstrumented equivalence matrix

`/tmp/p10l-base-{T}.json` and `/tmp/p10l-final-{T}.json`, median of three runs;
commands use `scripts/bench_scf.py`, all five cases, the absolute UPF path,
`--threads T --repeats 3 --baseline /tmp/p10l-base-T.json --require-identical`.
All five final commands exited 0. Every case has exactly one energy bit pattern
and one density hash across T=1,2,4,8,12, identical to P10K. No tolerance exception
was needed for the retained code.

| T | si8 P10K s | si8 final s | k444s P10K s | k444s final s |
|---:|---:|---:|---:|---:|
| 1 | 3.679080500 | 3.637859209 | 2.555586083 | 2.523534042 |
| 2 | 2.468224792 | 2.369165000 | 1.341805625 | 1.327141208 |
| 4 | 1.845511500 | 1.777932167 | 0.817114000 | 0.812098792 |
| 8 | 1.659135666 | 1.577597542 | 0.585230458 | 0.576186375 |
| 12 | 1.593232750 | 1.524704083 | 0.502583084 | 0.495796084 |

| Case | Final energy Ha | SCF iterations | Density SHA-256 |
|---|---:|---:|---|
| gamma | -7.2955304682107647 | 10 | `1854d1ea91d96933ec85e908ec62b4466ddb044abe16b060c1cbe376daea0d3f` |
| k222s | -7.919423412469083 | 10 | `2b310246ec91cf7ea8b845e03fcee0fe6c7769e5556c4af16ab5c8e4916e3ee5` |
| k222 | -7.8312300571211368 | 10 | `f203190f8c973771213a5074eb769c3ceecf264945aa805ea2da6edad5988840` |
| k444s | -7.9259761266480915 | 10 | `ad4d3eabe1b8823365908f6647c4470831a51416e8c738cd3f2aaabf7549b220` |
| si8 | -31.321803601031256 | 15 | `9152c7a9c64db67f623379a343a7b8c84b532ccc08736f8870292ad287e3e4de` |

### Nonlocal and memset investigation; MKL limits

The existing independent `nonlocal_traversal` release microbenchmark compares
fe0492b's channel-major traversal against the current cache-blocked traversal.
It exited 0 and preserved bits on every fixture. Measured Mac old/current us:

| Atoms | PW rows | Channels | Historical us | Current us |
|---:|---:|---:|---:|---:|
| 4 | 257 | 1 | 4.993 | 3.670 |
| 2 | 769 | 4 | 21.348 | 18.500 |
| 4 | 1025 | 4 | 39.243 | 36.188 |
| 4 | 4097 | 4 | 103.207 | 76.707 |
| 8 | 2969 | 4 | 91.582 | 89.639 |
| 8 | 4097 | 3 | 92.493 | 93.428 |
| 20 | 16385 | 9 | 2834.007 | 2829.397 |

Thus a universal restoration is not supported by this machine's multi-shape
measurements. No Nonlocal arithmetic/traversal change is retained. This does
not explain away the supplied dev08 MKL regression (0.98 versus 0.323 s): the
native host needs the same isolated comparison and matching compiler/binary
provenance. A per-block BLAS rewrite would change scalar association and needs
its own measured benefit and numerical validation; it was not silently adopted.

A definite per-call full-grid clear exists in `GridFusedPlan::inverse`:
`reciprocal.fill(0)` precedes scatter. It is **not redundant in the present
layout**: `apply_local` reuses that array for the dense forward FFT output.
Removing the clear would leak coefficients outside the input sphere into the
next inverse. Eliminating it safely would require a separate persistent sparse
input and dense forward-output buffer (an extra 16*Ngrid bytes per worker), or
a measured backend-specific gather interface. Matrix result zero-initialization
and grow/reshape fills are other possible memset callers. The supplied aggregate
0.382 s cannot identify one of them without caller stacks. No unsafe uninitialized
buffer or speculative zero-fill removal is retained.

Source provenance matters: in this worktree `FusedPlan` selects
`OrderedFusedPlan` for pure and `GridFusedPlan<FftPlan3d>` for MKL. The latter is
also a generic test adapter for Rust. Remote frames named GridFusedPlan should
therefore be checked against the actual selected backend/build before attributing
them to this pure production path.

MKL's 16 versus historical 15 SCF iterations remains **unresolved**, not asserted
to be a rounding-only effect. The rejected pure four-column experiment shows
that backend call shape alone can change Davidson stopping trajectories while
energy remains within 1e-12 Ha, but it is not proof of the MKL cause. Required
remote ablations remain separate restoration of historical established-basis
projection and historical MKL Gram grouping, with full residual/certificate
histories and the first differing stopping decision. Final retained pure code
is bit-identical to P10K; it does not purport to recover MKL's historical path.
No native MKL or dev08 executable was run by this implementer.

### Twenty-atom generality fixture

Six retained iterations, not converged/QE reference values. Skew lattice rows
(16,0,0), (1.5,12,0), (.5,.7,10) Bohr; 4x5 fractional xy sites, alternating
z=.35/.65; 20 Gaussian synthetic ions, valence 2, rc=1, 101 radial points,
cutoff 10 Ha, 24 bands. Local-only and s/s/p projector fixtures use the existing
test generator. Full density and energy bits agree across all five budgets.
Single-run wall times include construction; before/after isolates the dependent
MGS overlap from the otherwise identical retained scheduling/copy candidate.

| Nonlocal | T | Without overlap s | With overlap s | Energy Ha |
|---|---:|---:|---:|---:|
| false | 1 | 0.805417500 | 0.810090250 | -23.44535292535889681 |
| false | 2 | 0.494505792 | 0.490638708 | -23.44535292535889681 |
| false | 4 | 0.370452875 | 0.365748292 | -23.44535292535889681 |
| false | 8 | 0.331065542 | 0.314533167 | -23.44535292535889681 |
| false | 12 | 0.333021917 | 0.308875667 | -23.44535292535889681 |
| true | 1 | 3.377506458 | 3.426797042 | 6.23449291357981039 |
| true | 2 | 2.083691959 | 2.067152333 | 6.23449291357981039 |
| true | 4 | 1.480985542 | 1.460948584 | 6.23449291357981039 |
| true | 8 | 1.268658958 | 1.197445334 | 6.23449291357981039 |
| true | 12 | 1.219533333 | 1.167459084 | 6.23449291357981039 |

T=8 wall reductions are 5.0% (local-only) and 5.6% (multi-projector);
these bounded-work measurements support the scheduling change beyond silicon,
but do not establish convergence, a universal default speedup, or dev08 targets.

### Final counters and the unmet utilization goal

Single analysis runs, `/tmp/p10l-final-analysis-3/summary.json`, with complete
per-region rows in each `regions.csv`. SCF scope excludes initial construction.
Inline logical kernels are distinguished from scheduled parallel regions.
Region busy time includes timed serial work *inside* a parallel region; nested
serial events carry a flag so the per-worker SCF total counts that time once.
The baseline did not time internal MGS pieces/serial subphases separately;
aggregate busy fractions before/after are not a like-for-like CPU-usage
comparison. Acceptance speedups come from uninstrumented equal-work runs above.

| Case | T | SCF s | Scheduled parallel regions s | Between logical regions s |
|---|---:|---:|---:|---:|
| si8 | 1 | 3.424605542 | 0.000000000 | 0.605520507 |
| si8 | 2 | 2.230545667 | 1.557531814 | 0.569535877 |
| si8 | 4 | 1.641624208 | 0.983106487 | 0.556465733 |
| si8 | 8 | 1.411440416 | 0.752220099 | 0.555594001 |
| si8 | 12 | 1.395442709 | 0.726998972 | 0.563300710 |
| skew4 | 1 | 0.294820666 | 0.000000000 | 0.078859974 |
| skew4 | 2 | 0.197661208 | 0.112450328 | 0.076868960 |
| skew4 | 4 | 0.159526125 | 0.077134459 | 0.074173957 |
| skew4 | 8 | 0.141163125 | 0.059721708 | 0.073172412 |
| skew4 | 12 | 0.152169417 | 0.065534954 | 0.077949078 |

si8 T=8 measured busy fractions by physical pool worker (includes explicit
serial phases but excludes unclassified glue): 0.685466, 0.334529, 0.350771, 0.324671, 0.315340, 0.305952, 0.297281, 0.280141.
Pool residency covers 1.377504749 s, but scheduled
parallel kernels cover 0.752220099 s of 1.411440416 s.
Serial Davidson work inside the resident pool is real; the supplied 6.5%
outside-pool sample fraction does not by itself measure Amdahl's serial fraction.

| si8 T=8 phase | Inclusive wall s | Intersecting scheduled parallel s |
|---|---:|---:|
| projection | 0.260345453 | 0.133919039 |
| mgs_second_pass | 0.073817513 | 0.036319305 |
| mgs_first_pass | 0.062500011 | 0.050539807 |
| evd | 0.138514276 | 0.000000000 |
| hamiltonian | 0.310303127 | 0.248703837 |
| gram | 0.288288724 | 0.270342651 |
| gram_assembly | 0.012495899 | 0.012427814 |
| rotation_pair | 0.139588301 | 0.139588301 |

Projection spends 51.4% of its wall interval inside scheduled parallel regions.
49.2% of the dependent MGS phase overlaps such regions; this is elapsed
interval overlap, not an assertion that the dependent dot products are parallel.

Qualifying scheduled regions with more than eight pieces; denominator is all
eight workers, including idle ones. Report maxima as well as medians:

| Phase | Regions | Median max/mean | Worst max/mean | Fraction below 1.5 |
|---|---:|---:|---:|---:|
| mgs_first_pass | 1352 | 2.242879 | 5.548690 | 0.179734 |
| hamiltonian | 85 | 1.202447 | 1.410508 | 1.000000 |
| gram | 888 | 1.287531 | 2.698167 | 0.870495 |
| rotation_pair | 408 | 1.329977 | 1.541656 | 0.995098 |
| rotation | 494 | 1.338344 | 2.946802 | 0.759109 |

Hamiltonian now has up to 20 pieces instead of at most eight static groups;
all 85 qualifying large regions meet 1.5. Gram/rotation panels remain 12 for
si8 and six for skew4. Internal MGS has at most 19 future-column pieces and
a genuine dependency chain. Thus the requested 32--64 pieces and max/mean
below 1.5 in every phase are **not achieved**. Single/few-band correction calls
remain intrinsically short; inflating their work merely to populate workers
would not be an optimization. Refining GEMM reduced imbalance but changed
backend rounding and workload; it was rejected, not concealed.

For final si8 T=8, the 1.301732290 s Davidson wall interval decomposes as:

- 0.549512191 s outside scheduled parallel kernel regions (including inline
  kernels, copies, validation and EVD).
- 0.459437420 s sum of per-region measured worker-time means.
- 0.200245089 s sum of per-region max-minus-mean measured worker time.
- 0.092537590 s region wall time not accounted for by the busiest measured
  worker (scheduling, waits and unclassified glue; not all pure scheduler cost).

These four terms sum to Davidson wall time. The max-minus-mean term includes
serial dependency work inside a parallel region, not just unequal independent
pieces. Counter coverage and waiting prevent interpreting it as a CPU-sample
histogram. About 42.2% of *T=8* Davidson wall time is outside scheduled regions;
this is not a 42.2% T=1 immutable serial fraction.

The initial T=1 partition has .6954 s outside logical kernels and 2.7960 s
in them. If the former remains fixed, perfect eight-way kernel scaling gives
1.0449 s (3.34x SCF speedup); actual initial T=8 is 1.5110 s. Kernel time
.8347 s exceeds the ideal .3495 s by .4852 s, with .0191 s less nonkernel time.
This excess combines imbalance, scheduling, reduction/assembly and cache effects.
The outside-kernel remainder includes uninstrumented non-Davidson work and is
a conditional model, not a proof that it can never be parallelized. Including
setup reduces whole-process ideal scaling further. Neither this Mac model nor
the supplied dev08 ideal 5.04x establishes that 3x is structurally impossible
on dev08. The retained patch supplies a measured ~5% same-work improvement,
not evidence that dev08 >=2x/3x or MKL <=6.5 s targets have been reached.

### Gate record, provenance, and handoff

- `scripts/check.sh`: exit 0; build, test, fmt and clippy all PASS. Final log:
  `/tmp/p10l-check-complete.log`.
- `cargo check --features mkl --all-targets`: exit 0; type checking only,
  `/tmp/p10l-mkl-complete.log`.
- Analysis-feature context/attribution test: exit 0, one selected test passed.
  It checks stolen-work attribution, nested serial charging without increasing
  piece count, and restoration after panic. Log: `/tmp/p10l-analysis-nested-test.log`.
- Analysis-feature clippy for numerics/solvers/physics/engine, all targets and
  `-D warnings`: exit 0, `/tmp/p10l-analysis-nested-clippy.log`.
- Original scalar-MGS bit-pattern oracle and synthetic20 release test: exit 0.
  Full ordinary tests also pass through `scripts/check.sh`; synthetic20 is ignored
  in debug and its release command is documented in code and the performance README.
- Nonlocal multi-shape release benchmark: exit 0, `/tmp/p10l-nonlocal.log`.
- README regeneration: exit 0, generated files unchanged. README specification
  tests: exit 0, three passed; existing unclosed-file ResourceWarnings remain.
- `git diff --check`: exit 0. No commit created; branch remains
  `task/P10L-WORKER-UTILIZATION`, based on `2e12ab95acd5b4f8245733aaf1635f35b1b08e40`.

Early checks caught no-op guard `drop_non_drop`, two analysis-only clippy issues,
and exporter division by zero for clock-resolution-zero empty regions. These
were fixed and rechecked; undefined zero-duration busy fractions are now empty
CSV fields. The final counter adds Gram assembly and nested dependent-pass work
to region busy time, and explicitly avoids double-counting nested serial events
in SCF totals. Earlier `/tmp/p10l-final-analysis*` runs are superseded by `-3`.

The k222s T=2 three-run table contained a timing outlier. Seven rotating-order
paired checks resolve it: P10K .201413875 s, final .201784000 s (0.18% difference),
`/tmp/p10l-k222s-paired.json`; no numerical difference. No further scheduling
threshold was tuned to this noise.

Binary/UPF SHA-256:

- P10K `/tmp/p10k-final`:
  `fe04f720a3fd7f46c1f26f7192d09d5b12e872158519fedc54ce7e3fcf221c47`.
- Retained pipeline binary used for paired and three-repeat tables:
  `a3f02fe0973496e5749dec72aa2d0a6a3168a61bdfdfd8069810c02ea349b637`.
- Final rebuilt production binary after adding the disabled-by-default Gram
  assembly diagnostic scope:
  `123b01cc810e53a59d911497b5bb271d999e66987a654e95689868a597fd77a3`.
- Si UPF, external and not committed:
  `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.

Production binaries do not contain the `P10L,` trace marker. The final added
scope compiles to a zero-sized no-op; it does not change numerical work. Timing
claims above identify the measured binary rather than conflating binary hashes.

Remote acceptance remains blocked: no dev08 access/native MKL runtime in this
implementer environment, no >=2x dev08 result, no MKL <=6.5 s result, no resolved
16-versus-15 iteration cause, and no all-phase imbalance <1.5 result. The next
remote measurement should preserve full histories and caller-resolved memset/
Nonlocal profiles, not infer a cause from aggregate sample labels. The retained
code is a locally verified ~5% improvement and reproducible diagnostic apparatus,
not completion of all P10L targets.

Dashboard generation was prepared in `/tmp/p10l-dashboard` using a copy of
`scripts/gen_agents_json.py` and the authoritative parent registry (73 rows).
The four dashboard JSON datasets are present there, but private Claude Artifact
upload/dataset updates and authoritative session completion require the
orchestrator. Those files and new README reference-data types are outside this
ticket's owned files. No new scientific reference or fabricated dashboard
completion was written. The regular scheduled-report timer was not active;
progress updates were sent during this run. All changes remain uncommitted.

Final rebuilt production binary recheck: all five `bench_scf.py` commands
(T=1,2,4,8,12, one repeat each, `--require-identical`) exited 0;
`/tmp/p10l-current-{T}.json`. All energies and density hashes remain identical
to P10K and across budgets. This recheck follows the final diagnostic-only
scope addition; the three-repeat and paired timing tables retain their stated
binary provenance.

## P10M investigation (2026-11-05; copy changes verified, QR acceptance blocked)

Base: `3e1c350a00aa5d67d3ddd07ad1a2bf2d36f80939`. The supplied dev08
profile attributes 233 idle-worker samples and 59 busy samples to scalar
`subtract_projection` called by internal orthonormalization. Source inspection
confirms a column dependency: a new normalized column must be finalized before
its first-pass updates to later columns, and the next column's second pass
visits all accepted columns sequentially. Each scalar projection reduces four
lanes on one worker, then updates the vector on that worker. This is not a
projector-block loop. Existing external-basis projection already uses two GEMMs.

P10L stored Mac evidence reports 0.260345453 s inclusive projection wall,
0.133919039 s overlapping scheduled parallel regions (51.4%), and first-pass
max/mean worker busy imbalance median 2.242879, worst 5.548690. These are prior
measurements, not newly measured P10M results. First-pass pieces peak at 19;
2969 plane waves produce only 12 fixed 256-row GEMM panels. The planned trial
uses CholeskyQR2 for internal full-rank blocks and fixed 64-plane-wave projection
panels (47 pieces for 2969 rows), retaining an ordered partial sum and the
historical rank-revealing path when Gram conditioning is unsuitable. Rank
thresholds and solver policies must remain unchanged. Acceptance is pending
fresh SCF energy, iteration-count and timing gates.

Fresh baseline `scripts/check.sh` passed all four gates. A release-with-symbols
`sample PID 1 -file /tmp/p10m-base.sample` attempt exited 255: macOS refused
process inspection in the sandbox; no sampled profile was produced. Use the
supplied profile plus the compile-time P10L counter build below, without
requesting elevated access. Fresh five-case T=1 baseline is stored at
`/tmp/p10m-base-1.json`; raw reference binary is `/tmp/p10m-base`.

### Current disposition and new dev08 evidence

The architect supplied six dev08 MKL T=1 runs with approximately 230 samples
per run, comparing main `9fcb356` and `fe0492b`: subtract_projection self time
0.276 vs 1.674 s, zgemm kernels 1.516 vs 1.005 s, memcpy 0.662 vs 0.167 s,
and wall 6.589 vs 6.334 s. These are supplied remote measurements, not Mac
measurements or a new QE reference. `9fcb356` is already an ancestor of this
worktree's base `3e1c350`; no Hamiltonian edits were made here. Main subsequently
advanced through `3715d13` and `cfe298b` (dashboard and log correction). Git
metadata is read-only to the implementer: the architect must merge current
main before accepting the branch; the numerical revert is already present.

The default retained change is storage-only: one reserved live-prefix basis,
borrowed live H basis, in-place coefficient truncation, new-column-only append,
and direct contracted-restart output into an alternate H cache. The operator
API still requires owned contiguous blocks, so a single reusable narrow product
buffer remains for appended columns. No full-subspace concatenate or cached
first_columns copy remains in production Davidson. The legacy concatenate is
compiled only for the independent test oracle. All physical defaults and the
P10L orthogonalization schedule are retained. No new dependency or unsafe code.

The QR candidate uses reusable Gram/inverse/rollback buffers and faer's Cholesky
factor and triangular solver. The `sqrt(epsilon)*max(diag(Gram))` squared-pivot
check only rejects an unsuitable Cholesky trial; it does not change the MGS rank
threshold. CholeskyQR2 is interleaved with the existing two established-basis
projection passes. It preserves column order and was checked against scalar
MGS, including zero/duplicate/extreme-scale/nearly-dependent inputs. Pure Rust
projection pieces are fixed at 64 rows independent of T; 2969 rows yield 47.
MKL retains one unpadded sequential call per logical GEMM. CholeskyQR2 with
256-row panels was also measured to isolate QR from panel refinement.

The QR experiment is disabled through a compile-time constant unless
`RUSTFLAGS='--cfg feature="block-qr"'` is explicitly supplied. Default builds
reserve no QR-specific numerical buffers and retain 256-row GEMM storage.
It is not accepted: the six-step 20-atom skew synthetic iteration count rises
233 -> 329 (64-row panels), or 233 -> 311 (256-row panels). Neither variant is
an acceptable general default. No tolerance, SCF policy, eigensolver policy,
or system-specific predicate was adjusted to conceal the failed gate.

### Numerical gates and iteration counts

Default storage changes: five actual inputs at T=1,2,4,8,12 have exactly the
baseline energy, density SHA-256, main Davidson count and certificate count.
Baseline binary: `/tmp/p10m-base`; final storage binary: `/tmp/p10m-copy-final`.
The baseline is current-main numerical state (including `9fcb356`). Evidence:
`/tmp/p10m-base-1.json`, `/tmp/p10m-copy-{1,2,4,8,12}.json` and their logs.

| Case | Energy Ha (baseline = storage change) | Main Davidson | Certificate Davidson | SCF steps |
|---|---:|---:|---:|---:|
| gamma | -7.2955304682107647 | 41 | 8 | 10 |
| k222s | -7.919423412469083 | 204 | 21 | 10 |
| k222 | -7.8312300571211368 | 408 | 24 | 10 |
| k444s | -7.9259761266480915 | 1686 | 141 | 10 |
| si8 | -31.321803601031256 | 340 | 68 | 15 |

The QR64 trial's five-case maximum energy difference was 2.415845301584341e-13
Ha (si8). Four small cases retained their counts; si8 changed from main/cert
340/68 to 308/66. The final experimental switch preserves the baseline small
case path below 2048 rows because the all-size QR trial was slower there in
exploratory measurements. No firm timing ratio is inferred from those three-run
samples. Cross-T density identity held for QR in its diagnostic Si/skew runs.

The reproducible synthetic driver generates nonphysical fixtures with valence
2, one s projector (cubic) or two s plus one p projector (skew), 101 radial
points, 20 Ry, n+4 bands, and six SCF steps at T=1,2,4,8. These are bounded
nonconverged outcomes (CLI exit 3), not converged energy certificates. All 24
storage-change outcomes are bit-identical to baseline; reference and final
summaries are `/tmp/p10m-synthetic-base/summary.json` and
`/tmp/p10m-synthetic-copy/summary.json`. The checked-in driver reproduced them
at `/tmp/p10m-synthetic-reproduce/summary.json` with assertions. Certificates
are absent in these bounded runs.

| Synthetic case | Baseline/main Davidson | Storage change | QR64 | QR256 |
|---|---:|---:|---:|---:|
| syn4-cubic | 30 | 30 | 30 | 30 |
| syn4-skew | 91 | 91 | 91 | 91 |
| syn8-cubic | 34 | 34 | 34 | 34 |
| syn8-skew | 126 | 126 | 126 | 126 |
| syn20-cubic | 40 | 40 | 40 | 40 |
| syn20-skew | 233 | 233 | 329 | 311 |

QR64 synthetic20 cubic/skew six-step energy changes were -3.3707786713e-6 /
+2.0423614069e-4 Ha; those unfinished SCFs are not covered by the five converged
Si energy claim. Both variants still had identical outputs across their own
thread budgets. This distinguishes deterministic rounding from unchanged
nonlinear SCF/Davidson behavior.

### Exact work and explicit copy accounting

The iteration-local counter uses Arc context propagated with pool tasks;
it does not subtract process-global snapshots across concurrently running k
solves. The unit test checks inherited nested work and exact summation.
Production builds inline the counter calls to no-ops.

| si8 T=1 or 8 | Iterations | Logical 2mnk | Backend-call 2mnk | Explicit copied bytes | Mean copied bytes/iteration |
|---|---:|---:|---:|---:|---:|
| copy_before | 408 | 21834221954 | 21834221954 | 3083884288 | 7558539.922 |
| copy_after | 408 | 21834221954 | 21834221954 | 666653344 | 1633954.275 |
| qr_cost | 374 | 21107612646 | 21107612646 | 1194068784 | 3192697.283 |

The executed/logical GEMM ratio is exactly 1 for every counted iteration, not
an estimate from timing samples. The retained default has identical GEMM work
and 2,417,230,944 fewer explicit copy bytes over 408 iterations. Count evidence:
`/tmp/p10m-copy-{before,after}-counters/*/trace.csv`; complete per-iteration
records remain in P10M_COST lines. The QR-cost trial includes the copy changes
and has a different iteration count; its lower total GEMM count is not a
same-work kernel comparison and does not waive the synthetic failure.

No panel is padded at the adapter boundary, no plane-wave range is repeated,
and MKL has one direct full-shape call per logical product. Pure Rust reduces
in fixed ascending panel order, including the short tail. Extra Gram assembly
adds, scalar projections, factorization and opaque backend padding are excluded
from the 2mnk count. Copy counters include the instrumented EVD/orthogonalization and
Davidson sites; unmodified Hamiltonian/FFT copies are outside their scope.
GEMM counters cover GemmWorkspace products. The counters exclude zero filling
and backend-internal traffic, and count
bytes copied once (not read+write bytes). See README for exact definitions.
These structural checks satisfy the measurable MKL T=1 no-more-GEMM-work /
no-more-explicit-copy constraint for the default change, not a remote wall-time
guarantee. Shape-dependent MKL efficiency remains a remote profiling question.

### Projection-only worker counters (exploratory diagnostic runs)

These instrumented runs are not thirty independent SCF timing samples; do not
interpret their wall ratios as accepted speedups. Individual-region imbalance
summaries below each contain over 1400 regions, and still describe these runs
only. Idle-or-unmeasured mass = T*projection wall - instrumented worker busy
seconds, excluding nested serial records already charged to a region. It is
not a CPU-sample idle histogram. Concurrent builds affected some later counter
runs; only exact operation/byte counts from them are used for acceptance.

| Variant | Projection wall s | Idle-or-unmeasured worker s | Max pieces | Qualifying regions | Median max/mean | Worst max/mean |
|---|---:|---:|---:|---:|---:|---:|
| baseline | 0.267038538 | 1.539600755 | 19 | 2312 | 1.673728 | 5.509513 |
| qr256 | 0.214251005 | 1.154498498 | 12 | 1452 | 1.305537 | 3.270874 |
| qr64 | 0.227786746 | 1.192473612 | 47 | 1564 | 1.225035 | 2.605039 |
| copy_after | 0.282241040 | 1.633152970 | 19 | 2312 | 1.679886 | 8.000000 |

Per-worker busy seconds for the projection phase at T=8 (worker indices 0..7):

- baseline: 0.210527978, 0.056795983, 0.056569124, 0.055686320, 0.055976937, 0.053715539, 0.055420143, 0.052015525.
- qr256: 0.137545460, 0.060708678, 0.062058267, 0.061325070, 0.061337260, 0.061081930, 0.059573135, 0.055879742.
- qr64: 0.157963224, 0.067609208, 0.067961451, 0.067348098, 0.067835294, 0.066219513, 0.067589752, 0.067293816.

Fresh P10L-style baseline projection covered 0.137422054 s of scheduled
parallel work in 0.267038538 s inclusive wall. QR64 covered 0.124120132 s in
0.227786746 s. The latter has more finely divided work but also different
Davidson work and extra Gram/rotation boundaries. Merely obtaining 47 pieces
does not establish all-worker utilization or the dev08 2x target.

The pre-copy QR64 prototype also had the following **exploratory three-run**
Mac medians. They are preserved for experiment provenance, not firm speedup
ratios or acceptance evidence. Binaries: `/tmp/p10m-base`, `/tmp/p10m-qr64`;
raw samples: `/tmp/p10m-paired.json`. QR changed si8 Davidson work (408 -> 374).

| T | Baseline si8 s (3 runs) | QR64 si8 s (3 runs) |
|---|---:|---:|
| 1 | 3.681320250 | 3.464171167 |
| 2 | 2.369899125 | 2.195221750 |
| 4 | 1.761237042 | 1.638604750 |
| 8 | 1.574253167 | 1.446550917 |
| 12 | 1.530826542 | 1.416553500 |

The all-size QR trial's k444s medians increased in this exploratory sweep;
small blocks were restored to the existing MGS path before further experiments.
That guard is a data-size boundary already present in P10L, not a silicon label
or thread-count-dependent numerical choice. The separate 20-atom failure still
prevents enabling the tall-block candidate by default.

### Thirty-sample Mac timing: first storage version (superseded below)

Uninstrumented alternating before/after pairs; one warmup per binary/budget
was discarded, followed by 30 samples of each. No concurrent build or benchmark
was launched by this session. This is the shared Apple Silicon Mac, not quiet
dev08 acceptance. Raw data: `/tmp/p10m-times-30.json`; driver is checked in at
`crates/dft-numerics/src/pool/orthonormalize_bench.py`. Every one of the 720
measured SCFs matched baseline energy, density hash and both Davidson counts.
Ratios below use 30 samples per condition, unlike the earlier QR explorations.

| si8 T | Before median s | After median s | Before / after |
|---|---:|---:|---:|
| 1 | 3.673281958 | 3.697323437 | 0.993498 |
| 2 | 2.339866562 | 2.303031667 | 1.015994 |
| 3 | 1.974924750 | 1.932406208 | 1.022003 |
| 4 | 1.795078771 | 1.750488480 | 1.025473 |
| 5 | 1.710852270 | 1.672836020 | 1.022726 |
| 6 | 1.601592292 | 1.564837229 | 1.023488 |
| 7 | 1.619370397 | 1.587577063 | 1.020026 |
| 8 | 1.624617397 | 1.592291104 | 1.020302 |
| 9 | 1.610427000 | 1.577656375 | 1.020772 |
| 10 | 1.582111188 | 1.544229146 | 1.024531 |
| 11 | 1.584237833 | 1.545799604 | 1.024866 |
| 12 | 1.582336479 | 1.541288542 | 1.026632 |

The retained change is a small multithread wall improvement on this Mac, not
the requested orthogonalization redesign. T=1 median is slightly slower despite
the explicit copy reduction; no T=1 wall-time win is claimed. The structural
MKL gate is equal backend-call work and fewer explicit copy bytes, not an
inference from Mac wall time. The dev08 utilization and >=2x target remain
unverified. Thirty repeated samples improve this local comparison's evidence;
they do not eliminate shared-machine load, thermal or cross-budget effects.

### Follow-up: k444s regression and eliminating redundant restart work

The first storage version's 30-pair k444s gate did **not** pass at T=1:
median 2.548991146 -> 2.569225875 s; 28 of 30 paired samples were slower.
Raw data: `/tmp/p10m-k444-times-30.json`. T=8 medians were
0.597748395 -> 0.598572187 s and T=12 0.509781667 -> 0.507254771 s.
Thus the preceding si8 table is an intermediate version, not acceptance of
that version. The initial T=1 regression was not dismissed as noise.

Two further storage changes address unnecessary work. Empty live-prefix H
buffers now reserve capacity without zeroing unused columns, and replacement
copies extend directly into reserved storage without a preceding zero fill.
A contracted restart formerly recomputed the exact H-basis/coefficient product
already produced by the paired rotation. When the basis can overflow, that
product is saved in the reserved restart buffer before residual formation and
reused if the existing product-trust and orthogonality guards permit it. The
same backend result, shapes and floating-point order are preserved; no H action,
residual, tolerance, restart criterion or solver algorithm is changed. The
fallback still recomputes into reserved output. This is reuse of an existing
engine result; no new reference-code source was consulted.

The final candidate binary is `/tmp/p10m-final2`. All five real cases at
T=1,2,4,8,12 and all 24 bounded synthetic cases again match baseline energy,
density and Davidson/certificate counts exactly. Logs:
`/tmp/p10m-final2-gates.log`, `/tmp/p10m-synthetic-final2.log`.
Fresh final counters: `/tmp/p10m-final2-counters/*/trace.csv`.

Over 408 si8 iterations, backend-call and logical 2mnk both fall from
21834221954 to 21636961594;
14 duplicate restart GEMMs are removed. Explicit copy bytes are
736959264 total (1806272.706 bytes/iteration),
versus 3083884288 baseline (7558539.922 bytes/iteration).
The intermediate copy-only version copied fewer bytes but performed more GEMM
work; the final candidate trades a narrow saved-product copy for removal of
identical recomputation. Each individual iteration still has executed/logical
GEMM count exactly one.

Final-candidate k444s remeasurement, 30 alternating pairs at each budget,
after discarded warmups and with no concurrent local build/benchmark:

| k444s T | Before median s | Final median s | Before / final |
|---|---:|---:|---:|
| 1 | 2.536909562 | 2.532314750 | 1.001814 |
| 8 | 0.586903520 | 0.584165750 | 1.004687 |
| 12 | 0.506513229 | 0.501538563 | 1.009919 |

All three medians are no slower in this remeasurement; the previous T=1
regression is resolved in these observations. This is a small effect on a
shared Mac, not a large or cross-machine speedup claim. Raw paired samples,
phases and equality checks: `/tmp/p10m-final2-k444-30.json`.

A per-call executed/logical ratio of one does not detect duplicate *requests*
from a caller. The removed 14 restart requests account for 197,260,360 units
of 2mnk. The iteration trace comparison, together with unchanged density and
iteration counts, detects that solver-level redundancy separately from tile
padding/recomputation within one request. No claim is made about the theoretical
minimum total work or MKL's internal instruction count.

### Final verification commands and scope

All commands below were actually executed. Logs are local diagnostic evidence;
source, reproducer commands and summarized numbers are retained in this entry.

| Command | Exit | Evidence |
|---|---:|---|
| `scripts/check.sh` | 0 | `/tmp/p10m-final2-check.log`: build/test/fmt/clippy PASS |
| `cargo check --features mkl --all-targets` | 0 | `/tmp/p10m-final2-mkl.log`; compile only, no MKL runtime on Mac |
| `cargo test --release -p dft-numerics block_qr_si8_shapes -- --ignored --nocapture` | 0 | `/tmp/p10m-final-shapes.log`: one selected test passed; 2969 rows, widths 20/40/60, T=1/2/4/8/12 |
| `RUSTFLAGS='--cfg feature="worker-analysis"' cargo test -p dft-numerics iteration_costs_follow_workers -- --nocapture` | 0 | `/tmp/p10m-cost-test.log`: exact nested worker cost aggregation test passed |
| `RUSTFLAGS='--cfg feature="worker-analysis"' cargo clippy -p dft-numerics -p dft-solvers --all-targets -- -D warnings` | 0 | `/tmp/p10m-analysis-clippy.log` |
| `python3 -I /tmp/p10m_final2_gates.py` | 0 | `/tmp/p10m-final2-gates.log`: bench_scf.py, five cases, five budgets, baseline identity and iteration assertions |
| `python3 -I crates/dft-numerics/src/pool/orthonormalize_synthetic.py /tmp/p10m-final2 /tmp/p10m-synthetic-final2 /tmp/p10m-synthetic-base/summary.json` | 0 | 24 matched bounded outcomes; individual CLI exits are 3 (nonconverged), not falsely reported as converged |
| `python3 -I scripts/make_readme_figures.py` | 0 | Regenerated validation block/figures; no content changes |
| `python3 -I scripts/tests/test_readme_figures.py` | 0 | Three tests passed; existing ResourceWarnings about open files were emitted |

The paired timing commands are reproduced in the README and their JSON keeps
all 30 samples, phase timings, energies, hashes and iteration counts per binary
and budget. The debug guards for matrix prefix capacity and the production
solver's original numerical/certification tests are included in check.sh.
Analysis counters do not exist in ordinary builds. Explicit copy-site totals
exclude unmodified Hamiltonian/FFT sites and opaque library copies; they are
not interposed libc memcpy measurements.

The full P10M goal remains **BLOCKED**: CholeskyQR2 failed the bounded 20-atom
skew iteration regression, so the 47-piece orthogonalization path is opt-in
only and the default dependency chain remains. The default storage/reuse change
preserves numerical behavior and reduces measured work/copies, but it does not
establish dev08 >=2x/3x scaling or remote MKL wall time. Synthetic runs are six
steps, not convergence certificates. The architect must merge current main and
review/commit the uncommitted owned-file changes; git metadata is read-only in
this worktree. Session/dashboard files and reference datasets were not edited
outside ticket ownership. No UPF, QE source, manifest, tolerance, physics default
or nonlocal-Hamiltonian change is included.

### Final retained candidate: complete si8 T=1..12 timing

Thirty alternating pairs per budget, discarded warmups, uninstrumented
`/tmp/p10m-base` versus `/tmp/p10m-final2`. No local build or other benchmark ran
concurrently. All 720 measured runs passed energy/density/Davidson-count
identity assertions. Raw data: `/tmp/p10m-final2-si8-30.json` (driver exit 0).
This table supersedes the intermediate storage-only timing table above.

| T | Baseline median s | Final median s | Baseline / final |
|---|---:|---:|---:|
| 1 | 3.667245479 | 3.603509583 | 1.017687 |
| 2 | 2.338728145 | 2.317020833 | 1.009369 |
| 3 | 1.983455000 | 1.938644042 | 1.023115 |
| 4 | 1.795927646 | 1.762801312 | 1.018792 |
| 5 | 1.711750876 | 1.670288354 | 1.024824 |
| 6 | 1.598668833 | 1.565457459 | 1.021215 |
| 7 | 1.615405833 | 1.581073000 | 1.021715 |
| 8 | 1.626799208 | 1.585099250 | 1.026307 |
| 9 | 1.608345000 | 1.569704000 | 1.024617 |
| 10 | 1.580744979 | 1.542475771 | 1.024810 |
| 11 | 1.584643791 | 1.542218229 | 1.027509 |
| 12 | 1.580154188 | 1.539370937 | 1.026493 |

The final candidate has lower observed medians at all twelve si8 budgets and
at k444s T=1,8,12 in the separate 30-pair sweep (180 more measured runs, driver
exit 0). These are local observations; they do not establish dev08 scaling or
MKL runtime. The failed QR prototype remains disabled. The full ticket status
is BLOCKED, even though the retained default's numerical, iteration, operation
count, copy-count and sampled Mac no-regression checks pass.

Reproduction commands for the final paired sweeps:

```sh
python3 -I crates/dft-numerics/src/pool/orthonormalize_bench.py   /tmp/p10m-base /tmp/p10m-final2 /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF   /tmp/p10m-final2-si8-30.json --cases si8 --repeats 30
python3 -I crates/dft-numerics/src/pool/orthonormalize_bench.py   /tmp/p10m-base /tmp/p10m-final2 /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF   /tmp/p10m-final2-k444-30.json --cases k444s --threads 1,8,12 --repeats 30
```

SHA-256 provenance (binaries built with `CARGO_PROFILE_RELEASE_DEBUG=true`):

- `/tmp/p10m-base`: `0a53062534cfd52090387e7206514985cb24924dc80709ee7c8c7d8d8c093d4e`.
- `/tmp/p10m-final2`: `292185da8edb6a938bbdc1d1c6eba4680ee887a1b5eb259f6ce98f395e085abc`.
- `/Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF`: `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
### Correction, 2026-11-04: NonlocalOperator::add comparison value
Commit 9fcb356 ("restore original NonlocalOperator::add loop, blocked traversal measured 2.9x slower on dev08") rests on a wrong reference value. The remote measurement session corrected it: fe0492b `NonlocalOperator::add` is 0.474 s (7.5 percent, 17 of 227 samples, 6 runs) on dev08 MKL T=1 si8, not 0.323 s (6 of 117 samples, a count too small for the stated precision). The P10I regression is about 2.0x, not 2.9x. The revert of 9fcb356 is bit-identical and reduced the function from 0.929 s to 0.800 s (239 samples) but the wall time did not move beyond noise (mkl T1 6.612 to 6.589 s, pure T8 4.567 to 4.540 s). The remaining MKL T=1 gap to fe0492b (6.589 versus 6.334 s) is spread: MKL zgemm kernel +0.512 s, memcpy +0.494 s, NonlocalOperator::add +0.325 s, tile multiply_into +0.221 s, apply_local +0.219 s, concatenate +0.192 s, offset by subtract_projection -1.398 s and old gemm -0.307 s. Rule applied from now on: never quote a profile-derived ratio as a firm number when the counts are below about 30 samples.

## P13C: correctness extension for NLCC, no optimization claim

NLCC-off retains the original XC evaluator and summation order. New combined
XC-density storage and spin scratch are prepared at construction; scalar core
arrays are immutable. Radial shell caching is setup-only. No core allocation,
logging or I/O occurs in orbital/Hamiltonian loops. Existing allocating XC
library/output boundaries remain; this is not a wholly allocation-free SCF.
NLCC evaluates output XC once for the common density before k audits and records
core integrals as audit-only values. No numerical defaults or solver policies
were changed, and no silicon-only performance default was adopted.

The original HEAD CLI was archived/built for paired legacy checks. Benchmark
commands and results are in tests/reference/datasets/p13c-status.json. A macOS
`sample PID 1 1 -file /private/tmp/p13c-before-profile.txt` attempt on that control
binary failed with exit 255 (process inspection unavailable in the sandbox);
the underlying k444s run exited 0. No profile finding or speedup is claimed.
The older historical baseline's density bits already differ from the control
HEAD, despite energies within 1e-12 Ha; acceptance uses the recorded same-HEAD
pair as well as explicitly reporting the historical-gate failure.
