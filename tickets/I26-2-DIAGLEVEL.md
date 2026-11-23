TASK_ID: I26-2-DIAGLEVEL
GOAL: |
  Read `gh issue view 26` is not possible in the Codex sandbox (no network), so the issue text relevant to you is quoted below. Issue #26 (user-filed 2026-11-24, reviewed on main 2f91471a): "Align SCF acceptance tests, profiling and benchmark instructions with the current design". Gates for every item: retain independent fixed-point verification and charge failed attempts and recovery to delivered wall; keep the strict 1e-10 Ha and standard 5e-9 Ha classes; multi-system evidence before any default promotion; the independent-source policy (no QE source). Defaults unchanged. Scratch inside the worktree. Do not edit scripts/delegate.sh. Leave changes uncommitted; return docs/schemas/report.schema.json.
  Item: "Make the diagnostic detail level selectable on the reference driver. The map forces options.profile = true (crates/dft-engine/src/refscf/map.rs about lines 144-157). Each inner solve allocates per-band Ritz histories (crates/dft-solvers/src/ref_davidson.rs about 145-164), the driver clones them, and the CLI encodes the iterations in two places (crates/dft-cli/src/run.rs about 1693-1705 history JSON and 1775-1786 top-level JSON). Make inner clocks and full Ritz traces honor the run-owned diagnostic setting. Preserve required production output, lightweight counters, the verification record and completion marker. Preserve or version the trace-reader contract. Demonstrate identical scientific output and verification on the same input, complete diagnostic-on records, and report trace bytes, output bytes and paired external wall. Include all required output in the final timing."
TASK: implement it (default for production: lightweight; full traces only when the run asks, through the existing run config/profile mechanism), version the trace-reader contract if its shape changes (update refscf-check/refscf-map readers and tests), prove identical scientific output and verification on/off on si8/si16/si32 (Mac), report trace bytes and output bytes on/off, and write docs/performance/lanes/RSCF/diaglevel-dev08-request.md for the paired external-wall A/B on dev08 (si64 and si128 np8, A/A floor first, crossed pairs, |mean|/SE >= 3, all output included in the timing).
OWNED_FILES:
  - crates/dft-engine/src/refscf/
  - crates/dft-solvers/src/ref_davidson.rs
  - crates/dft-cli/src/run.rs   # reference-driver output section only
  - crates/dft-cli/src/refscf_check.rs
  - crates/dft-cli/src/refscf_map.rs
  - crates/dft-cli/tests/
  - crates/dft-engine/tests/refscf.rs
  - docs/design/reference-scf.md   # a "Diagnostics level" section only
  - docs/performance/lanes/RSCF/diaglevel-dev08-request.md
TEST_COMMANDS:
  - scripts/check.sh full
ACCEPTANCE_CRITERIA: selectable level, identical science on/off, byte numbers reported, dev08 A/B request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
