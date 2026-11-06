# P26D / R4: real ORC Davidson library

Base: `995f53479d27af07b7315336db4f77718ec3f030`, branch
`task/P26D-REAL-GAMMA-R4`. Changes intentionally uncommitted for architect review.
Derived from the approved `real-gamma-design.md`, sections 4, 8 and 11, and its
own-words study-note inputs. No QE source read, copied or used as an oracle.

## Delivered contract

- `RealLinearOperator` lives in dft-solvers, with scalar/block application and an
  optional shared persistent pool. R5 supplies the RealHamiltonian adapter.
- `RealDavidson` uses bands-wide correction blocks and capacity
  `min(dimension, max(max_subspace, 6*bands))`. There is no independently tunable
  block size and no retuned capacity. Warm starts are orthonormalized; missing
  rank is completed with seeded real directions. Cold starts are coordinate
  units sorted by kinetic energy, with ORC half-list order breaking ties.
- Cached HV, one old/new rectangle plus new square for the incremental symmetric
  reduced matrix, real symmetric faer eigen solve, and real GEMM Ritz rotations.
  The old reduced prefix is retained bit for bit. Only new basis columns need H.
- Corrections use `-r/(d-lambda)`, with sign-preserving 0.1 Ha floor (zero positive).
  Projection uses two blocked passes with MGS/normalization between them, and
  two MGS passes within each incoming block. No CholeskyQR2 or scoped OS threads.
- Restart contracts V and HV when the rotated vectors remain orthonormal;
  otherwise it orthonormalizes and reapplies H. A failed certificate invalidates
  stored products and uses the actual residual for the next correction.
- Strict prefix and buffer semantics match complex Davidson, including the
  iteration-limit exception for unconverged buffers. Every successful strict
  return is truly certified. Standard Ritz-change mode is implemented entirely
  here: the first empty mandatory root uses the buffer threshold, consecutive
  changes use strict `<`, and results are explicitly uncertified. No R6 solver
  deferral is needed; engine policy/dispatch still belongs to later tickets.
- H sees normalized columns, never tiny raw residuals alongside unit vectors.
  Debug assertions enforce the design's pair norm ratio of at most 1e3. The large
  test operator also independently checks each supplied norm against unity.

The existing real numerical adapter is reused unchanged. A complex MAC costs
8 real flops, so the former 262,144-MAC work guard represents 2,097,152 flops.
A real MAC costs 2 flops, giving 1,048,576 real MACs for the same nominal work.
This is a derived dispatch heuristic, not a measured optimal threshold or an
assertion of equal real/complex throughput. The adapter requires four 256-row
tiles before parallel dispatch and combines partial sums in increasing tile
order. Tests cross this threshold and confirm no new workers are created during
solves. Projection scratch and GEMM tile storage are reserved before iteration;
solver-owned matrix results and the existing eigen wrapper still allocate at
orchestration boundaries, as does the complex solver. This is not a claim of an
allocation-free complete eigensolve. No logging, file I/O or allocation is added
inside scalar arithmetic loops or GEMM/projection kernels.

## Tolerances (unchanged)

| Quantity | Value | Reason |
|---|---|---|
| Kinetic denominator floor | 0.1 Ha | Existing complex preconditioner, prevents singular correction |
| Dependent normalized column | <=1e-12 | Existing shared dimensionless rank threshold |
| Restart contraction orthogonality | <=1e-12 | Existing complex cache validity rule |
| Fixed-H test residual | <=1e-11 Ha | Tight iterative solve, below the energy accuracy required of these bounded operators |
| Dense/known-spectrum eigenvalue error | <1e-12 Ha | R4/design fixed-operator acceptance gate |
| Test orthogonality defect | <1e-12 | Same dimensionless contraction gate |
| Buffer comparison | strict 1e-9 Ha, buffer 1e-6 Ha | Exercises distinct per-band contracts without changing solver defaults |
| Incremental matrix comparison | <1e-14 Ha | Unit-scale algebra rounding check; old prefix and symmetry additionally bit-exact |
| Pair norm ratio | <=1e3 | Design section 6 protects paired FFT application |

The solver has no built-in convergence tolerance default; the caller supplies
one. Tests do not loosen any solver threshold on failure.

## Fixed-operator observations

Apple Silicon arm64, macOS 27.0.1 (26A434), ordinary pinned release build.
Single-run times below are **Mac observations**, not speed acceptance, not a
statistical comparison, and not dev08 evidence. Commands run from this worktree:

```sh
cargo test -p dft-solvers --test real_davidson -- --nocapture
cargo test --release -p dft-solvers --test real_davidson -- --include-ignored --nocapture --test-threads=1
```

The 49-dimensional dense fixtures are (0) deterministic random Householder
rotations of a known nondegenerate spectrum, (1) diagonally dominant kinetic
plus potential, and (2) random rotations of triply degenerate shells, with the
4-band cut inside a shell. DenseEigensolver on the complexified real operator
is the dense oracle; fixtures 0 and 2 also have analytically prescribed spectra.
Both iterative solvers start from the exact same real coordinate units.

| Fixture | Bands | Steps real/complex | H vectors real/complex | Real restarts / contractions | Max eigenvalue error (Ha) | True residual max (Ha) | Orthogonality | Release seconds real/complex (O) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Random symmetric | 3 | 173 / 177 | 384 / 388 | 25 / 25 | 2.5091040e-14 | 9.1118786e-12 | 3.1086245e-14 | 0.006798 / 0.007222 |
| Kinetic + potential | 3 | 17 / 17 | 48 / 48 | 2 / 2 | 2.6645353e-15 | 6.8196185e-12 | 2.2204460e-15 | 0.000785 / 0.000840 |
| Degenerate shells | 4 | 130 / 130 | 468 / 480 | 22 / 22 | 2.4424907e-15 | 7.7885113e-12 | 1.8651747e-14 | 0.008453 / 0.009157 |

The random case's four-step difference cannot be attributed to the cold start:
the starts are identical. Real vs complex EVD, real reduction association and
one-rectangle symmetric filling introduce different rounding paths; individual
roots cross residual thresholds on different iterations and change correction
widths/restart trajectories. The degenerate case similarly has equal iteration
counts but different H-vector totals; a degenerate eigenspace has no unique
Ritz basis. These differences are consistent with those arithmetic changes,
not isolated proof that one particular kernel caused them. All final eigenvalues
agree with complex Davidson within 1e-12 Ha. No parameter was changed to force
matching counts. The three-case step totals are 320 real vs 324 complex; this is
fixed-operator evidence only, not an SCF-wide iteration claim.

Additional behavior comparisons:

| Test | Real / complex outcome |
|---|---|
| Strict prefix + loose buffers | 8 / 8 steps, 31 / 31 H vectors, equal convergence flags, certified |
| Ritz-change | 6 / 6 steps, 17 / 17 H vectors, equal convergence flags, both uncertified |
| Deliberately stale first H block | 2 / 2 steps, 10 / 10 H vectors; real makes 2 certificate attempts and a noncontracting restart |
| Limit with exact strict roots, unconverged buffer | Both return certified with `[true,true,false]`; making buffer mandatory fails |
| Full-space diagonal | Strict 1 step / 10 H vectors; standard 2 steps, uncertified |

The stale-block probe is intentionally not an exact fixed operator on its first
call; it tests the defensive certificate path, not a physical approximation.
Additional tests check finite/dimension errors, sign at the denominator floor,
rank-deficient/tiny input columns, incremental matrix preservation and cold-start
ordering. Tiny columns of scale 1e-200 exercise the stable norm fallback.

## Determinism

A 65-dimensional fixture checks eigenvectors, values, residuals and histories
bit for bit at T=1,2,4,8, with identical step and H counts. The ignored release
test uses a sparse symmetric 1025-dimensional kinetic/potential operator with
17 bands, capacity 102, local and distant symmetric couplings. It asserts that
`1025 * max_basis_dimension * 17 >= 1048576`, so it exercises actual parallel
real GEMM dispatch, not just the serial guard. It checks worker-creation counts
before and after solving and the same full-state bit identity.

| T | Steps | H vectors | Release seconds (O) |
|---:|---:|---:|---:|
| 1 | 13 | 91 | 0.018275 |
| 2 | 13 | 91 | 0.013473 |
| 4 | 13 | 91 | 0.011758 |
| 8 | 13 | 91 | 0.010321 |

No timing acceptance or threshold retuning follows from these observations.

## Isolation and verification

The baseline source was exported with `git archive HEAD` to
`/tmp/p26d-baseline`, built with the same release profile and vendored crates,
and its CLI copied to `/tmp/p26d-baseline-cli`. An initial concurrent baseline
build failed with exit 101 because the module declaration preceded creation of
the new file; the isolated HEAD export/build succeeded with exit 0. It is the
successful isolated binary, not the failed attempt, used for the baseline.

Baseline and final release CLI SHA-256 are identical:
`15e8e8aee9d3807a1e25022df0b1f8d5629050ef72517953d253670634e871ae`.
External Si UPF SHA-256:
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
No UPF is committed. The engine has no reference to the new solver.

The baseline command is:

```sh
python3 -I scripts/bench_scf.py --binary /tmp/p26d-baseline-cli \
  --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF \
  --cases gamma,k222s,k222,k444s,si8 --threads 1 --repeats 1 \
  --record /tmp/p26d-baseline-t1.json
```

Final comparisons use `target/release/dft-cli`, the same cases and UPF,
`--baseline /tmp/p26d-baseline-t1.json --require-identical`, at each T=1,2,4,8,12.
The benchmark's printed speedups are ignored: these are isolation checks and
some overlap workspace validation. Only the isolated solver times above are
reported as observations.

Verified results (all exit codes from completed processes):

| Command / check | Exit | Result |
|---|---:|---|
| `scripts/check.sh dft-solvers` | 0 | build, test, fmt, clippy PASS |
| `scripts/check.sh` | 0 | workspace build, test, fmt, clippy PASS |
| Release integration test command above | 0 | 7 passed, including ignored parallel-dispatch test |
| `python3 -I scripts/tests/test_readme_figures.py` | 0 | Existing generated evidence remains current |
| `git diff --check` | 0 | No whitespace errors |
| Baseline benchmark command above | 0 | All five cases recorded |
| Final benchmark, T=1,2,4,8,12 (five invocations) | 0 each | 25/25 energies and densities identical |

A separate Python check of recorded JSON verifies energy **bits**, density
SHA-256 and SCF iteration counts are equal in all 25 comparisons (exit 0).
Reference values are this engine's same-HEAD baseline, **not QE results**:

| Case | Energy (Ha), unchanged | SCF iterations, unchanged |
|---|---:|---:|
| gamma | -7.29553046821076 | 10 |
| k222s | -7.9194234124690865 | 10 |
| k222 | -7.831230057121135 | 10 |
| k444s | -7.9259761266480915 | 10 |
| si8 | -31.32180360103127 | 15 |

Local raw records: `/tmp/p26d-baseline-t1.json`,
`/tmp/p26d-final-t{1,2,4,8,12}.json`, `/tmp/p26d-bench-exits.json`,
`/tmp/p26d-final-release-tests.log`, `/tmp/p26d-final-solvers.log`, and
`/tmp/p26d-check-all.log`. These temporary records are not durable repository
artifacts; the numerical tables and reproducible tests above are the owned-file
record. A source comparison stripping only the permitted additions proves all
pre-existing lib.rs text is byte-identical to HEAD (exit 0).

## Scope limits and handoff

No engine, dispatcher, options, CLI, locking, spectral guard or RMM changes.
No dft-physics dependency, unsafe code or new crate. The lib.rs diff contains
only the new module/export and RealLinearOperator trait. Existing complex
items and shared numerical code are unchanged.

A true residual certificate proves eigenpair accuracy, not that all lowest
roots were found (the existing P15C caveat). Dense/known-spectrum tests check
that additional property only for their fixed fixtures. No real-path SCF,
material, MKL or dev08 acceptance is claimed. R5 supplies adapter/engine
integration; R6 still owns standard-profile engine wiring; R7 owns broader
performance evidence.

Ledger handoff (owned-file restriction prevents editing the shared ledger,
backlog, generated README evidence or dashboard): R4 carries over A12 restart
contraction and A55 incremental reduction in an isolated real path. A33/A35
scoped threads and A52 CholeskyQR2 remain rejected. Fixed spectra pass, all
three contraction fixtures exercise restarts, and real/complex counts need not
match exactly even with the same start. Keep the complex default unchanged;
no rejected optimization is revived and no default adoption is requested.
