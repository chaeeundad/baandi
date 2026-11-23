TASK_ID: I25-STOP-TO-VERIFIED
DEPENDS_ON: I26-1-RANKRULE and I26-2-DIAGLEVEL merged (same files; start only after both are on main)
GOAL: |
  Issue #25 comment (2026-11-23, quoted because the sandbox has no network): "For the planned Si128 pure-MPI8 run and its Si64 control, report the first map that meets the stopping estimator, the map that passes independent verification, and wall seconds between those events. Break that interval into (1) all verification attempts, including the passing one, (2) ordinary SCF work after a failed attempt, and (3) any unassigned time. Use non-overlapping intervals on the same wall clock and report the attempt count. If the first verification passes, recovery work is zero. If none passes, mark the run unverified." Also: keep the passing criterion unchanged (fixed-point density check); keep verification on; report the diagnostic on/off delta (I26-2) separately from this interval because they can overlap.
TASK: |
  1. In the reference driver record, add per-map wall timestamps on one clock (seconds since driver start): map start, verification start/end (when an attempt runs), map end. Lightweight: present at every diagnostic level. Version the record if the reader contract requires it.
  2. Add a reader (in refscf-check, or a small `refscf-check --stop-to-verified` mode) that prints: first stop map, passing map (or UNVERIFIED), attempt count, interval wall, and the three non-overlapping parts (verification attempts, ordinary SCF work after a failed attempt, unassigned = interval minus both), with a check that the parts sum to the interval.
  3. Tests on Mac: si8 displaced (the case that fails first verification at map 10 and passes at map 12) and one case whose first attempt passes (recovery must be 0). Defaults and scientific output bit-identical.
  4. Write docs/performance/lanes/RSCF/stop-to-verified-dev08-request.md: si128 and si64 pure MPI8, same binary, verification on, diagnostics at production level, all output in external wall; job names per gate-on-dev08.md, one case per fast job (<=600 s).
OWNED_FILES:
  - crates/dft-engine/src/refscf/   # record timestamps only
  - crates/dft-cli/src/refscf_check.rs
  - crates/dft-cli/tests/   # new test file only
  - crates/dft-engine/tests/refscf.rs
  - docs/design/reference-scf.md   # a "Stop-to-verified interval" section only
  - docs/performance/lanes/RSCF/stop-to-verified-dev08-request.md
TEST_COMMANDS:   # Mac rule (user, 2026-11-24): only these; never scripts/check.sh full or release matrices on the Mac
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-cli <your new test name>
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-cli --all-targets -- -D warnings
  - (orchestrator, on dev08) scripts/check.sh full and the release si8 displaced run
ACCEPTANCE_CRITERIA: timestamps recorded at all diagnostic levels, parts sum to the interval, si8 displaced shows nonzero recovery and a first-pass case shows zero, defaults bit-identical, dev08 request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
