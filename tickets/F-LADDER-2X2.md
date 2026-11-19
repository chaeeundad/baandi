TASK_ID: F-LADDER-2X2
GOAL: |
  Causal test of the si32 handoff cliff (docs/performance/plan-2026-11-19-followup.md section 3; read it). Source facts: the first standard fine solve seed is min(previous requested tol, 0.1 * last_coarse_std_error_Ha / n_electrons) clamped (crates/dft-engine/src/cutoff_ladder.rs ~L1423); si32 128 electrons, std_err 6.8e-5 / 3.15e-3 Ha give 5.31e-8 / 2.46e-6 (matches the dev08 records). The Mac study (docs/design/ladder-fraction.md) explained the si32 loss by both the earlier coarse stop and the looser seed; a clean 2x2 intervention is still missing.
TASK: (1) Mac release, one core, si32 (and si64 if affordable, else si16 as held-out) seed 601: capture the two natural handoff states at fractions 0.50 and 0.40 (input density + warm orbitals). Continue each twice with the tight and the loose scalar observed in the two natural runs (use the existing test injection point: production seed setter session.rs ~L1934; the offline driver cutoff_ladder.rs ~L2035 assigns c.standard_next unconditionally, so add an explicit seed override to the ignored offline driver only; never as a production default). Same fine basis, Hamiltonian construction, state import, solver, mixer and certificate. Diagonal arms (0.50+tight, 0.40+loose) must reproduce the unsplit natural controls (report bit/within-1e-12 agreement). Keep the 0.45 natural trajectory as cross-check. (2) Per map report effective requested tolerance, first redo, H and block calls, mixer resets/refits, Coulomb residual, charge, rejected attempts, full verified seconds; fine-only time and the full price of constructing each starting state. Add geometries (3 seeds) for the interaction question. (3) Decision table from the plan: tightening the 0.40 state recovers maps -> handoff-aware tolerance candidate (information available during the run only); 0.40 stays expensive at both tolerances -> representation/mixer history, report common-low-G density error and occupied principal angles; interaction -> report all four cells, keep default. If a tolerance rule is supported, implement it behind a default-off switch and evaluate on >=7 systems x 3 seeds with total H and |dE| to the strict reference. (4) Also fix two harness faults from the review: scripts/ladder_handoff.py alias `si64` falls through to benchmarks/cases/si64.toml (fixture is si64-gamma.toml; add alias and a test); the offline minimal-grid arm starts the timer before the original coarse problem and builds a second coarse problem (asymmetric setup): fix or document the timing. (5) docs/design/ladder-2x2.md with all four cells and the dev08 command lines (crossed pairs, order-alternated, A/A floor). Defaults unchanged, bit-identical. Scratch in worktree (TMPDIR_990="$PWD/target/scratch").
OWNED_FILES:
  - crates/dft-engine/src/cutoff_ladder.rs
  - crates/dft-engine/src/session.rs   # seed override accessor only if needed, default off
  - crates/dft-numerics/src/switches.rs
  - crates/dft-numerics/tests/all/switches.rs
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - scripts/ladder_handoff.py
  - scripts/ladder_2x2_*.py
  - docs/design/ladder-2x2.md
  - docs/performance/attempts-ledger.md
READ_ONLY_CONTEXT: AGENTS.md, docs/performance/plan-2026-11-19-followup.md, docs/design/ladder-fraction.md, tickets/F-LADDER-FRACTION-msg1.md
TEST_COMMANDS:
  - scripts/check.sh fast dft-engine dft-numerics dft-cli   # 3 CLI cache tests reject in-repo TMPDIR, pass on dev08
ACCEPTANCE_CRITERIA: four cells reported with controls reproduced; decision-table branch stated with data; defaults bit-identical; |dE| to strict reference reported (strict 1e-10 / standard 5e-9).
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
