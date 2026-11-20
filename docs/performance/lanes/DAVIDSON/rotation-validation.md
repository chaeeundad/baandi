# Follow-up validation: serial ROTATE_COLUMNS

Base: 2bde16a (UPPER_BORDER candidate committed by the orchestrator). The
rejected UPPER_BORDER implementation and phase2-job.sh are retained unchanged, default off. New work
is left uncommitted. No Mac timing is used as performance evidence.

## Implementation and invariants

`real_davidson::rotation_pair` replaces only the non-fused serial pair inside
`ritz_and_residuals` when DAVIDSON_ROTATE_COLUMNS=on. It calls the existing
`RealGemmWorkspace::multiply_column_groups` twice (32 fixed output columns,
whole PW rows), retaining X and H X and the unchanged residual calculation.
The OFF branch calls the original rotate_pair_slices with identical inputs.
Fused paths retain priority. The distributed workspace does not call this
helper, so its rotation algorithm and collectives are unchanged.

The column-group adapter writes output columns directly on MKL. Its groups
and inner dimension do not depend on worker count; on pure builds it uses the
existing library/pure rotation adapter. No new FFI, unsafe block, dependency,
workspace buffer, residual shortcut, tolerance or physical threshold. Changes
can round differently on MKL; Mac bit identity is not a claim about MKL.

`rotation_pair` is an additive labelled record under DAVIDSON_PROFILE; no new
clock reads when profile is off. Shape-derived MKL subcall/copy counts are
explicit models. Existing GEMM dispatch counters move FLOP from rotate to
other (adapter label); total FLOP remains unchanged. This must not be reported
as eliminated floating-point work.

## Commands and real exits

All build/test scratch used `TMPDIR=$PWD/target/scratch` on the 990 volume.
PATH includes `/opt/homebrew/opt/rustup/bin`.

| Command | Exit | Result |
|---|---:|---|
| cargo build --release -p dft-cli | 0 | final candidate built |
| cargo test -p dft-solvers --lib column_rotation | 0 | 1 new test |
| BAANDI_DAVIDSON_ROTATE_COLUMNS=on cargo test --release -p dft-solvers --test lean_davidson --test distributed_lean_algebra | 0 | 13 tests |
| CHECK_SCCACHE=0 TMPDIR_990="$PWD/target/scratch" scripts/check.sh fast dft-solvers dft-numerics dft-engine dft-cli | 1 | 616 passed, 3 known failures, 160 ignored; build/fmt/clippy PASS |
| python3 -I scripts/tests/test_lane_davidson_parse.py | 0 | 2 tests |
| python3 -I docs/performance/lanes/DAVIDSON/test_tools.py | 0 | 3 synthetic schema/statistics tests |
| bash -n docs/performance/lanes/DAVIDSON/rotation-job.sh | 0 | syntax only, not a dev08 run |
| git diff --check | 0 | clean whitespace |

Gate failures are the same three pseudo cache tests that reject repository
TMPDIR: cache_resolution_and_no_network_side_effects,
cached_fetch_checks_hash_and_provenance, and
fetch_boundary_verifies_before_publication_and_never_overwrites. No skips or
looser tolerances added. The optional ps load report is sandbox-denied; gates
continue. Full gate is the orchestrator's pre-push step, not claimed here.
Raw logs: `target/scratch/phase2-rotation/{build,unit,solvers-on,check}.log`.

The new unit test covers n=1/257/513, m=37, k=65 (both row-panel and column-group
tails), both products against the original kernel, independent scalar dots,
workers 1/3 bit identity, and non-finite rejection. Existing solver tests cover
cubic, skew, clustered/degenerate and synthetic multi-projector operators,
independent dense eigensolutions, fallbacks, workers and simulated ranks 1..4.
The MPI candidate is inactive; these tests ensure serial-reference comparisons
and existing distributed results retain their contracts.

## Scientific evidence

`rotation-correctness.json` records hashes, energies, statuses and paths.
Inputs: fast_screen_dav.py --emit-input, seed601, fractional displacement0.001,
T=1, same approved external Si.pz-vbc.UPF (not committed). Standard entry:

```
env BAANDI_DAVIDSON_ROTATE_COLUMNS=off|on BAANDI_DAVIDSON_PROFILE=off \
  BAANDI_PHASE_PROFILE=off BAANDI_STAGE_TIMERS=off BAANDI_DAVIDSON_STAGES=off \
  BINARY run INPUT --output-dir OUTPUT
```

UPPER_BORDER remains at its default off throughout. Strict adds
BAANDI_REFERENCE_STRICT=on. Saved pre-edit binary, new OFF and new ON all match
scientific bits and density bytes on **Si8, Si32, Si64**. Skew OFF/ON also
converges and matches. Si8 strict converges with identical science, history
including fixed-point certificates, and density bytes. abs(dE)=0 Ha in all
pairs. The Phase 0 standard-only verifier lacks strict's `verified` field,
so its strict invocation returned 1; explicit strict-schema comparison of
energy/eigenvalues/occupations/history/certificates/density passed. This is
not a numerical mismatch.

Slab exits 3 in BOTH arms after 100 steps, with identical history, final
energy and density. It is not a converged generality pass. No tolerance was
relaxed. Si32 ON with workers 1/3 also has identical science and density.

Additional Si32 profile OFF/ON files validate only the extended schema and
route/counter arithmetic; their clocks were not analyzed. Both match the
non-profiled science. Parser and rotation-counter exit 0. There are 87 pairs
in either arm with unchanged 23,845,356,496 FLOP. The MKL call model changes
6138 -> 382 and extra-copy model 1,178,558,720 B ->0. These are **Mac-shape
schema checks**, not actual MKL calls, memory-controller measurements or dev08
speed evidence. A real MKL run must also show the registered actual dispatch
route/FLOP movement before the model is applicable.

## Remaining scope

Diagnosis uses relayed dev08 counts/dominant shapes and rates, with corrected
pair/EVD arithmetic and separate serial/MPI time bases. Exact accumulated
roofline gaps need joint per-shape counts/times; mode x count is not used as
measured total work. The implementation/request is no longer blocked on that
refinement. Real MKL/MPI A/B and adoption remain the orchestrator's work. The
new request is the active alternative; UPPER_BORDER is rejected and must not be rerun.


## Rejection-response validation reuse

The numerical ROTATE_COLUMNS implementation has not changed since the tests
above. No new tile-narrowing variant was implemented in response to the
UPPER_BORDER loss. Reuse the recorded binary hash/scientific evidence instead
of repeating numerical tests on unchanged code. This response changes the
ledger/diagnosis/request and extends the Python parser with exact per-map
EVD work accounting. Parser tests now total3, including multiple maps, full
EVD and B5 attempts, plus the original malformed-input cases. The lane helper
tests remain3. No Mac wall comparisons were made.
