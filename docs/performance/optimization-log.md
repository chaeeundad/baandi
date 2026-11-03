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
