# P33A: reusable real Davidson storage and fixed column work

**Status: BLOCKED.** Bit-identity and workspace checks pass; the strict
zero-allocation gate fails inside the unchanged EVD backend. This is a reviewable
partial implementation with a measured dependency blocker, not acceptance.

Baseline main: `73f59e347e5b6eb74d8cfdbbe8ce87d55f7132c2` (this worktree's
starting HEAD), not the earlier `72574f9` profiled on dev08. Candidate is
uncommitted on `task/P33A-REAL-WORKSPACE-PARALLEL`. Raw commands, binaries,
inputs, snapshots, stdout and timings are in `P33A-evidence/`, outside target.
No QE source was read, no remote command was executed, no artifacts were edited.
The task's single-target instruction overrides the brief's separate-target rule:
reference binaries were copied before rebuilding, and their hashes retained.

## Motivation and scope

The supplied **relayed** dev08 profile of `72574f9` places 48.3% of working
samples in algebra, 6.98% in memset and 5.65% in memcpy, with poor T=8 worker
utilization. This motivates storage lifetime and independent-column work, not
another arithmetic-changing GEMM shape experiment (P32A/A49) or CholeskyQR2
(A52). Capacity, defaults, tolerances, complex kernels and pair FFT scheduling
are unchanged. The protected capacity expression was compared verbatim.

The controlled synthetic operator is a symmetric tridiagonal with diagonal
`1 + 0.01*i` Ha and adjacent entries `0.03` Ha. Dimensions are (2969,20) and
(11877,68); five solves per thread budget 1/2/4/8/12, strict residual threshold
1e-10 Ha. Each solve takes 21 iterations and respectively 60/156 H-vector
applications. This cheaply repeatable operator does not reproduce the full
Si8 subspace-width/rank distribution and is not a material performance oracle.
Timers are exclusively in test wrappers; nested durations are subtracted, and
no production solver/kernel contains a timer, logger or file write.

`tests/real_workspace.rs` freezes the baseline solver control flow and compares
its eigenvalues, vectors, residuals, history and work counts with production.
It shares numerical adapters; independent saved SCF bytes and the P32A-shaped
kernel oracle are the before/after arithmetic gates. A generated current
wrapper (`P33A-evidence/current-probe.rs`, generator `make_probe.py`) provides
current routine attribution without instrumenting production. Its temporary
integration-test copy is removed after measurement.

## What `matmul_1_4_dyn` means

The proposed interpretation as necessarily a one-column or dimension-one
product is incorrect. In vendored `nano-gemm-codegen/src/lib.rs`, x86 f64 AVX
has N=4 lanes; `nano-gemm/src/lib.rs::Plan::from_masked_impl` computes
`mr = MR_DIV_N * N`. The symbol's first number counts SIMD row vectors, not
scalar rows. A **4-by-4 output with dynamic reduction length** selects the
1-by-4-vector kernel. P32A's ordinary rotation microcalls are exactly 4x4xK
(K up to 256), explaining this symbol without a one-column caller.

Both Ritz rotations and projection reconstruction already pass whole panels to
`multiply_into`; `bounded_matmul` subdivides them to stay within the 4096-MAC
allocation-free backend path. Reduced/projection Gram products can legitimately
have one new/surviving column late in convergence, and final row/column tiles
can be partial. There is no serial per-band GEMM caller to combine into a panel.
Changing the bounded shape would reopen P32A's rounding/packing tradeoff and is
not part of this bit-identical task.

## Storage and arithmetic contracts

| Site | Before | Candidate |
|---|---|---|
| Basis V / HV append | allocate vector, copy old and new columns | initialize capacity once, expose active shape, copy only new columns |
| Projection result prefix | `to_vec` allocation/copy | swap initialized output storage; active prefix is a borrowed view |
| Ritz V / HV | allocate and zero N x bands each step | reuse capacity storage, existing GEMM Replace writes output |
| Residual | clone HV plus allocate norms | overwrite reserved residual and norms |
| Correction / selected roots | zero new block / collect vector | overwrite active columns / clear and extend reserved index vector |
| Reduced matrix | allocate/zero m x m and copy prefix | keep prefix, expand compact stride backwards in place, overwrite new rectangle |
| EVD | allocate vectors, diagonal and backend scratch each call | reserve all before loop; same contiguous input/output and sequential faer call |
| Random fallback / flags / history | per-step blocks/flags | reserved before loop; random draw order unchanged |
| Gram assembly | clear output then serial ordered sums | fixed four-column pieces, overwrite with `0 + first`, then same ordered sums |

`RealMatrix` keeps initialized backing storage and an active contiguous shape;
views, indexing, equality and cloning see active entries only. Reshaping never
allocates or initializes. EVD scratch covers the maximum requirement over **all
positive active sizes**, not an assumed monotonic requirement. Size zero needs
no backend scratch. Tests poison reused EVD/GEMM outputs with NaNs to catch
incomplete overwrites. Reduced stride repacking retains the exact old layout
rather than introducing a possibly different strided backend path.

Column work has four-column leaves, identical at every budget. Residual norms
keep their old four-accumulator order; corrections and projection normalization
keep their scalar expressions. Gram partials remain the same 256-row pieces,
and each output entry adds them in ascending order, including its initial +0.
MGS order, norm fallback, rank decisions, sign conventions and EVD are unchanged.

### Allocation gate: FAILED, backend boundary isolated

A test-only counting allocator starts at the return of the initial H block,
which now occurs **after all solver workspace construction**, and stays active
through every step, restart, certificate and successful return. It counts
alloc, alloc_zeroed and realloc across all threads, and delegates unchanged to
System. Its only unsafe code is this documented test allocator adapter, as
explicitly permitted by the ticket's request for a counting-allocator test;
production contains no new unsafe or instrumentation.

The zero-allocation assertion fails. First T=1 runs counted 91 allocations for
2969/20 and 1499 for 11877/68, despite all solver-owned buffers being reserved.
Warm-start 97/8 cases also cover strict restarts/certification and standard
Ritz-change returns. All budget results, including the failing assertion, are
retained in `final-allocation.log` (exit 101). The ignored test is not weakened:
its zero assertion remains the unmet acceptance criterion.

A separate diagnostic isolates `RealEighWorkspace::solve_into`, with all of its
inputs/outputs/scratch already constructed: sizes 20/40/88/160 count
0/32/88/121 allocations. Both direct real GEMM operations count **0** at
2969x40 with 20 output columns. Source inspection identifies EVD's internal
Householder/tridiagonal/divide-and-conquer matrix products, which can reach
faer's general GEMM path. The supplied EVD MemStack does not supply GEMM packing
or bookkeeping. In vendored `gemm-common/src/gemm.rs`, `did_pack_lhs` is an
allocated vector in general GEMM even with the thread-local packing slab warm.
Pool/runtime initialization can contribute additional thread-dependent counts;
the full-solve counts are not all attributed to a single backend line.

Thus a source-only no-allocation argument would have been incorrect. Solver
storage is bounded (history reserves max_iterations; selected indices reserve
bands; matrix shapes are bounded by initialized capacity; result values/flags
are moved), but **the entire loop is not allocation-free**. The final result is
BLOCKED on gate 5, not a completed optimization or performance acceptance.

A library change exposing caller-owned packing/bookkeeping, preserving the
current kernels' arithmetic, needs a separate dependency/backend ticket and
orchestrator review. Global allocator replacement, global packing thresholds,
copying an eigensolver, or changing EVD/GEMM algorithms to hide the allocation
would exceed this task's contracts. None was done. The reusable EVD boundary
is explicitly documented as potentially allocating internally.

## Validation and reproduction

Baseline outputs were saved before production edits. The additive engine test
`gamma_workspace_snapshot` records energy/density/all eigenvalue bits, per-SCF
Davidson/H counts and complete per-k solver results for the four R5 synthetic
cells and Si primitive/Si8, at T=1,2,4,8,12. It does not relax the separate
historical real-versus-complex SCF-count gates. Snapshot comparison is bytewise.

Commands (prefix cargo commands with `/opt/homebrew/opt/rustup/bin` in PATH):

```
cargo test --release -p dft-solvers --test real_workspace routine_split -- --ignored --nocapture --test-threads=1
cargo test --release -p dft-solvers --test real_davidson -- --include-ignored --nocapture --test-threads=1
DFT_SI_UPF=/external/Si.pz-vbc.UPF P33A_SNAPSHOTS=/absolute/output cargo test --release -p dft-engine --test gamma_workspace_snapshot -- --ignored --nocapture
python3 -I P33A-evidence/measure.py before
python3 -I P33A-evidence/measure.py final
scripts/check.sh
```

The measurement driver invokes bench_scf for gamma/k222s/k222/k444s/si8 at all
five budgets, with `--require-identical` on the final comparisons, then keeps
one warmup and five measured real-on runs per Si8/Si16 profile/budget. All raw
results and densities are retained. Si16 uses the benchmark suite's 1x1x2
orientation, not a substituted dev08 orientation. Timings are sequential shared
Mac observations, without task-owned builds overlapping; five samples cannot
resolve a 1–3% improvement. No dev08 speed or utilization claim is made.

“Same trajectory” below means optimization before/after **with gamma_real on**.
It does not mean gamma_real on equals the complex off path: those trajectories
already differ in main. The off path is checked separately against its own
same-HEAD baseline.

## dev08 request

After review, compare the saved starting HEAD and candidate on the same quiet
node, strict/standard Si8 and Si16, T=1/2/4/8/12, equal accuracy and unchanged
work counts. Preserve binary/input/UPF hashes and raw scientific outputs.
Profile working and waiting samples separately, resolve the 4x4 dynamic kernel
callers, and measure allocator/copy shares and MGS/EVD/Gram assembly shares.
At least 30 samples per named site; do not infer scaling from this Mac. Compare
against QE's best Gamma configuration only with its documented accuracy/input
provenance. No default or algorithm change is requested.

## Gate outcomes and limits

- Gate 1: 30/30 saved binary snapshots match byte for byte (four R5 synthetic
  fixtures plus Si primitive/Si8, five budgets). Energy difference is 0 Ha;
  density/eigenvalues/SCF/Davidson/H work and solver state match.
- Gate 2: all 25 final complex/off bench comparisons pass, each bench command
  exit 0, energy difference 0 Ha and identical density. Work counts are also
  compared from the recorded JSON.
- Gate 3: `scripts/check.sh` exit 0: build, test, fmt, clippy all PASS.
- Gate 4: eight unchanged RealDavidson fixed-operator tests pass, including the
  ignored release cases; added real algebra/poisoned-buffer tests pass.
- Gate 5: **FAIL**, counting-allocator test exit 101; the exact zero requirement
  remains asserted. No numerical tolerance or allocation threshold was relaxed.
- README generator and freshness test exit 0; generated content is unchanged.
  No new accepted reference result was invented to update the README.

The only additive plumbing outside the owned implementation/test files is
`crates/dft-engine/tests/gamma_workspace_snapshot.rs`, an ignored output-snapshot
harness reusing the R5 fixtures. It changes no production engine code.
Raw evidence is untracked under P33A-evidence and should be archived before
worktree cleanup. Capacity-sized backing storage trades peak reserved memory
for stable ownership; peak RSS was not measured. The total serial fraction on
SCF and dev08 allocation/utilization changes remain unmeasured.

<!-- P33A measured tables -->
## Controlled routine attribution

Exclusive accumulated seconds over five solves. `run`/`run_into` is projection
outside its separately timed GEMMs/MGS. “Other” includes setup, corrections,
convergence logic, validation, timers and allocations not assigned a wrapper.
These are region-wall measurements, not CPU utilization or dev08 sample shares.

| Shape / T | Routine | Before calls | Before (s) | After calls | After (s) | After share |
|---|---|---:|---:|---:|---:|---:|
| 2969x20 / 1 | gram | 305 | 0.026738 | 305 | 0.027397 | 22.25% |
| 2969x20 / 1 | rotate | 410 | 0.065684 | 410 | 0.067021 | 54.43% |
| 2969x20 / 1 | eigh | 105 | 0.005270 | 105 | 0.005060 | 4.11% |
| 2969x20 / 1 | mgs | 200 | 0.007472 | 200 | 0.007450 | 6.05% |
| 2969x20 / 1 | append | 200 | 0.002761 | 200 | 0.000173 | 0.14% |
| 2969x20 / 1 | residual | 110 | 0.003891 | 110 | 0.003429 | 2.78% |
| 2969x20 / 1 | multiply | 210 | 0.000720 | 210 | 0.000011 | 0.01% |
| 2969x20 / 1 | reduced_matrix | 105 | 0.000086 | 105 | 0.000027 | 0.02% |
| 2969x20 / 1 | run | 100 | 0.003235 | 100 | 0.003289 | 2.67% |
| 2969x20 / 8 | gram | 305 | 0.015562 | 305 | 0.017520 | 22.37% |
| 2969x20 / 8 | rotate | 410 | 0.032441 | 410 | 0.035358 | 45.14% |
| 2969x20 / 8 | eigh | 105 | 0.003200 | 105 | 0.003468 | 4.43% |
| 2969x20 / 8 | mgs | 200 | 0.004938 | 200 | 0.005770 | 7.37% |
| 2969x20 / 8 | append | 200 | 0.001389 | 200 | 0.000044 | 0.06% |
| 2969x20 / 8 | residual | 110 | 0.002726 | 110 | 0.002949 | 3.76% |
| 2969x20 / 8 | multiply | 210 | 0.000594 | 210 | 0.000022 | 0.03% |
| 2969x20 / 8 | reduced_matrix | 105 | 0.000029 | 105 | 0.000018 | 0.02% |
| 2969x20 / 8 | run | 100 | 0.002109 | 100 | 0.008514 | 10.87% |
| 11877x68 / 1 | gram | 305 | 0.206004 | 305 | 0.201279 | 14.76% |
| 11877x68 / 1 | rotate | 410 | 1.037902 | 410 | 1.024073 | 75.08% |
| 11877x68 / 1 | eigh | 105 | 0.025162 | 105 | 0.024708 | 1.81% |
| 11877x68 / 1 | mgs | 200 | 0.021616 | 200 | 0.022321 | 1.64% |
| 11877x68 / 1 | append | 200 | 0.032755 | 200 | 0.000459 | 0.03% |
| 11877x68 / 1 | residual | 110 | 0.036726 | 110 | 0.032851 | 2.41% |
| 11877x68 / 1 | multiply | 210 | 0.005028 | 210 | 0.000012 | 0.00% |
| 11877x68 / 1 | reduced_matrix | 105 | 0.000273 | 105 | 0.000106 | 0.01% |
| 11877x68 / 1 | run | 100 | 0.008250 | 100 | 0.008777 | 0.64% |
| 11877x68 / 8 | gram | 305 | 0.131413 | 305 | 0.119923 | 20.74% |
| 11877x68 / 8 | rotate | 410 | 0.358598 | 410 | 0.343124 | 59.34% |
| 11877x68 / 8 | eigh | 105 | 0.025901 | 105 | 0.024820 | 4.29% |
| 11877x68 / 8 | mgs | 200 | 0.021744 | 200 | 0.022944 | 3.97% |
| 11877x68 / 8 | append | 200 | 0.018118 | 200 | 0.000916 | 0.16% |
| 11877x68 / 8 | residual | 110 | 0.036853 | 110 | 0.012800 | 2.21% |
| 11877x68 / 8 | multiply | 210 | 0.005162 | 210 | 0.000009 | 0.00% |
| 11877x68 / 8 | reduced_matrix | 105 | 0.000285 | 105 | 0.000124 | 0.02% |
| 11877x68 / 8 | run | 100 | 0.008391 | 100 | 0.010125 | 1.75% |

### What remains serial

| Shape / T | EVD | MGS | EVD + MGS measured serial share |
|---|---:|---:|---:|
| 2969x20 / 1 | 4.11% | 6.05% | 10.16% |
| 2969x20 / 8 | 4.43% | 7.37% | 11.79% |
| 11877x68 / 1 | 1.81% | 1.64% | 3.45% |
| 11877x68 / 8 | 4.29% | 3.97% | 8.26% |

EVD and dependent two-pass MGS remain sequential. Restart orthogonality checks,
finite validation, compact reduced-stride repacking, small flag/index work and
ordered scalar norms in H validation also remain serial. This synthetic cold
operator has no restart, so it cannot assign a measured share to restart work.
Gram/rotate timings mix pooled kernels with serial input validation and region
overhead. Consequently EVD+MGS is a measured serial **lower bound**, not a full
Amdahl fraction or a claim that everything else is parallel. Warm restart paths
are covered separately by fixed-operator and allocation tests.

## Final Mac observations

| Case / profile / T | Before median (s) | Workspace-only (s) | Final median (s) | SCF before/after | Davidson before/after | H vectors before/after | Same trajectory |
|---|---:|---:|---:|---:|---:|---:|---|
| si8-strict-T1 | 1.184494 | 1.184805 | 1.182609 | 15/15 | 253/253 | 2001/2001 | True |
| si8-strict-T2 | 0.874180 | 0.862015 | 0.868063 | 15/15 | 253/253 | 2001/2001 | True |
| si8-strict-T4 | 0.769906 | 0.756071 | 0.767473 | 15/15 | 253/253 | 2001/2001 | True |
| si8-strict-T8 | 0.762807 | 0.746050 | 0.751955 | 15/15 | 253/253 | 2001/2001 | True |
| si8-standard-T1 | 0.508547 | 0.514513 | 0.506753 | 17/17 | 89/89 | 877/877 | True |
| si8-standard-T2 | 0.392897 | 0.392050 | 0.390304 | 17/17 | 89/89 | 877/877 | True |
| si8-standard-T4 | 0.348995 | 0.349637 | 0.349471 | 17/17 | 89/89 | 877/877 | True |
| si8-standard-T8 | 0.336735 | 0.338961 | 0.343691 | 17/17 | 89/89 | 877/877 | True |
| si16-strict-T1 | 5.160023 | 5.077582 | 5.067058 | 14/14 | 264/264 | 3254/3254 | True |
| si16-strict-T2 | 3.593589 | 3.518632 | 3.506389 | 14/14 | 264/264 | 3254/3254 | True |
| si16-strict-T4 | 2.981895 | 2.912192 | 2.865164 | 14/14 | 264/264 | 3254/3254 | True |
| si16-strict-T8 | 2.710891 | 2.650844 | 2.614851 | 14/14 | 264/264 | 3254/3254 | True |
| si16-standard-T1 | 1.983218 | 1.980625 | 1.966025 | 16/16 | 88/88 | 1521/1521 | True |
| si16-standard-T2 | 1.454854 | 1.450029 | 1.452491 | 16/16 | 88/88 | 1521/1521 | True |
| si16-standard-T4 | 1.248036 | 1.245086 | 1.231696 | 16/16 | 88/88 | 1521/1521 | True |
| si16-standard-T8 | 1.151732 | 1.156270 | 1.142311 | 16/16 | 88/88 | 1521/1521 | True |

All 16 rows retain energy bits, density hashes and work counts. Full eigenvalue
and solver-state identity is independently covered by the 30 binary snapshots.
The workspace-only ablation and final column scheduling effects are small and
mixed on this Mac; these batches do not establish a general speedup. No change
is recommended as a performance acceptance until dev08 measurement and the
unmet backend allocation gate are resolved.
