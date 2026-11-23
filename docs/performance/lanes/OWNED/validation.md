# G-OWNED-ROWS validation (Mac correctness only)

Base f526588, changes uncommitted. No Mac timings support a performance claim.
No real MPI, MKL or dev08 execution; no QE comparison was performed. Evidence:
[mac-evidence.json](mac-evidence.json), source unit tests in
`crates/dft-engine/src/gamma_real.rs::owned_input_tests`.

## Executed commands

Environment for gates: `TMPDIR_990="$PWD/target/scratch" CHECK_SCCACHE=0`.
Cargo commands use `/opt/homebrew/opt/rustup/bin` in PATH and `--offline`.

| Command | Exit | Result |
|---|---:|---|
| `gh issue view 20`, `gh issue view 24` | 1 each | Network unavailable; both issues subsequently read through connected GitHub API. |
| `cargo test --offline -p dft-engine --lib owned_input_tests -- --nocapture` | 0 | 2 fast tests pass, 2 long tests ignored. |
| `DFT_SI_UPF="$PWD/../pseudo/Si.pz-vbc.UPF" cargo test --release --offline -p dft-engine --lib owned_input_tests -- --include-ignored --nocapture --test-threads=1` | 0 | 3 tests passed at execution: full real-Si matrix, synthetic/fallback/strict matrix, packing oracle. Descriptor invalidation test subsequently added and passed in debug/full gates. |
| `scripts/check.sh full` (final) | 1 | Build/fmt/clippy pass; 797 tests pass, 3 pre-existing CLI scratch-policy failures, 183 ignored, 87 binaries. |
| `scripts/check.sh dft-engine` (final) | 0 | All four gates pass; 265 tests pass, 0 fail, 77 ignored, 33 binaries. |
| `python3 -I scripts/make_readme_figures.py` | 0 | Regenerated outputs; no tracked README/figure changes. |
| `python3 -I scripts/tests/test_readme_figures.py` | 1 | 8/9 pass. Unmodified `test_structural_scope_and_sources` expects obsolete text `standard profile is unsupported`; generator no longer emits it. Freshness assertions pass. |
| `python3 -I docs/performance/lanes/OWNED/prepare_input.py --case si128 --upf ../pseudo/Si.pz-vbc.UPF --output target/scratch/owned-si128.toml` | 0 | Uses canonical benchmark renderer; no hand-written Si128 geometry. |
| `target/debug/dft-cli validate target/scratch/owned-si128.toml --json` | 0 | Si128 input preflight. |
| `target/debug/dft-cli grid target/scratch/owned-si128.toml` | 0 | Si128 basis/grid preflight; no SCF. |
| `git diff --check` | 0 | Whitespace clean. |

Full-gate failure names:
`pseudo::tests::fetch_boundary_verifies_before_publication_and_never_overwrites`,
`pseudo::cache_resolution_and_no_network_side_effects`,
`pseudo::cached_fetch_checks_hash_and_provenance`.
All reject scratch caches inside repositories, exactly the known lane-charter
limitation with worktree-local scratch. No test, tolerance or cache policy was
changed to suppress these failures. The ancillary `ps` load probe is sandbox
blocked; all four gates still ran. Thus **the full gate is not claimed green**.
No edits outside OWNED_FILES were made to global registry/dashboard, README
sources or these failing tests. This worktree has no G-OWNED-ROWS registry entry
or `.delegate/G-OWNED-ROWS/report.json`; master owns that orchestration update.

## Identity coverage

Real external Si UPF SHA-256:
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
Si8/16/32 conventional supercells and the existing si-skew primitive fixture,
30 Ry, standard path with fixed-point verification, capacity-per-band 3,
ranks 1/2/3/8 simulated. Six synthetic two-rank off/on pairs cover local-only
and nonlocal pseudopotentials, materialized standard warm starts, forced M1
fallback, and strict-reference controls. **22 A/B pairs, 68 per-rank records**.

For each same-rank-count A/B pair, compare energy and density bits, complete
IterationRecord Debug serialization after clearing only PhaseTimings (including
certificate timing), and the complete delivered EigenSolution. Float Debug
roundtrips finite binary64 values; no scientific history field is omitted.
All standard runs pass verification. Maximum energy difference **0 Ha**;
density/history/eigen-solution comparisons exact. Differences across rank counts
are permitted by ADR 0041 and are not used as an A/B identity claim.

Rank 1 deliberately uses the serial path (`SharedComm::active` filters size 1):
zero eligible calls/bytes, identical output. Direct single-owner packing is
separately covered by the packing oracle. Real Si multi-rank runs have one
eligible atomic initial block each; later warm starts are already lazy. Synthetic
PPs lack atomic orbitals, so the harness disables lazy returns for those cases
only to exercise materialized warm inputs. Fallback/strict cases have zero
eligible bridge counts. No production default changed for these controls.

For every successful eligible counter record the test independently reconstructs
the basis and stick-owner map from cell/cutoff/ranks, counts L_r, reads delivered
band count B, and checks saved=8GB*calls and allocated=8L_r B*calls.
Counter versus model equality is exact, not merely A allocated minus B allocated.

| Fixture | MPI ranks with active switch | Saved per rank per eligible entry (B) |
|---|---|---:|
| Si8 | 2/3/8 | 475,040 |
| Si16 | 2/3/8 | 1,703,520 |
| Si32 | 2/3/8 | 6,461,088 |
| si-skew | 2/3/8 | 43,456 |

The scalar pack oracle checks cubic/skew/slab shapes, two cutoffs, actual stick
row sets for 1/2/3/8 ranks, empty owners, signed zero/subnormals, zero-mode
threshold equality/failure, off-owner NaN/Inf, conjugacy failure, overflow and
out-of-range row errors, against the existing full HalfSphere::pack. Validation
must complete before writing output. Descriptor tests change cell, positions,
charge/PP, cutoff and FFT layout; potential-only updates reuse the descriptor.
New rank layouts partition the current descriptor exactly; no cross-session cache
or mid-session communicator migration was introduced. Existing lazy-serial and
persistent-state tests also remain in the gate, including lazy M1 fallback.

## Intermediate failures retained honestly

- Initial compile: RealMatrix has no `column_mut`; replaced by its explicit
  column-major mutable slice. Final compile/gates use the corrected code.
- First full-matrix harness incorrectly demanded active counts at rank 1 and
  on synthetic inputs that never materialized; corrected path expectations and
  explicit synthetic materialized-start coverage. No numeric tolerance changed.
- A newly invented sheared conventional Si8 fixture diverged in **baseline A**
  under the standard defaults. Replaced by the repository's existing si-skew
  fixture, rather than relaxing divergence/convergence tolerances. Its failure
  remains in `target/scratch/owned-matrix.log`; it is not a qualified fixture.
- Initial full gate also reported duplicate test-fixture module import. The
  existing scf fixture module is private and out of scope; a narrow test-only
  duplicate_mod allowance documents reuse without changing unrelated visibility.
- Final full gate still has the three scratch-policy failures above. README
  suite still has its unrelated expected-phrase mismatch. No green full gate.

No runtime memory allocator instrumentation or Mac peak RSS was collected.
The counters cover requested matrix payloads only. Si64/Si128 byte savings are
static predictions from the validated formula, pending dev08 RSS and wall.
