TASK_ID: G-PHASE-VALIDATE
GOAL: |
  plan-2026-11-19-followup.md section 4.2: the owned projector constructor still scans, for every atom, every plane-wave momentum with a three-term finite-dot check (crates/dft-physics/src/hamiltonian.rs ~L2356; engine call scf.rs ~L1568; phase materialization ~L342). Si64 at fraction 0.5: 64 x (23871 + 8409) = 2,065,920 pair checks per rank (source-derived; confirm by instrumenting).
TASK: (1) Instrument call count and exclusive time of this validation per rank (si64 on the Mac one core and 2-rank MPI if available; report seconds and fraction of init). (2) Implement the exact fast path behind a default-off switch: scan components of G and atom positions once for finiteness and max |.|; if both maxima <= 2^500, each product <= 2^1000 and a three-term sum cannot overflow, so every finite-dot test passes; otherwise run the original loop with its exact error behavior. Keep factor/norm checks, the phase computation, lazy-table guarantees and collective ordering. Tests: NaN, infinity, huge imbalanced inputs, signed zero, multiple species, large finite components with cancellation; acceptance = identical stored coefficients (bit) and identical error behavior. Count pair checks, fast accepts, fallbacks. Only gate this on numerical range, never on physics. (3) Report the measured saving (isolated timer) and whether it is resolvable (>= 3x timer noise); if not resolvable say so. (4) docs/design/phase-validate.md and a ledger row. Defaults unchanged until dev08 A/B; update switch-count tests (switches.rs tests, dft-cli run.rs).
OWNED_FILES:
  - crates/dft-physics/src/hamiltonian.rs
  - crates/dft-engine/src/scf.rs   # call site only
  - crates/dft-numerics/src/switches.rs
  - crates/dft-numerics/tests/all/switches.rs
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/design/phase-validate.md
  - docs/performance/attempts-ledger.md
READ_ONLY_CONTEXT: AGENTS.md, docs/performance/plan-2026-11-19-followup.md
TEST_COMMANDS:
  - scripts/check.sh fast dft-physics dft-engine dft-numerics dft-cli   # 3 CLI cache tests reject in-repo TMPDIR
ACCEPTANCE_CRITERIA: bit-identical coefficients; edge-case tests; measured instrumentation numbers; default off.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
