# Optimization log (newest first)

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
