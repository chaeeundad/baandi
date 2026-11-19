TASK_ID: F-LADDER-HANDOFF
GOAL: |
  dev08 (master, si64 MPI np=8, 6c31696 defaults): the coarse stage of the cutoff ladder costs 1.2007 s = 31 % of the time outside ordinary fine maps (3.93 s of a total budget 4.16 s = 0.90 x the other code); coarse solve 1.0406 s, density 0.0982, energy 0.0495, potential 0.0216. The ladder is adopted and stays (ledger LADDER-REEVAL-DEV08: 3.480 s with vs 4.083 s without at Si32 np8), but its handoff rule is a fixed coarse residual threshold. GitHub issue #23 (text: tickets/issues/ISSUE-23.md) asks whether the CURRENT orbitals' unresolved fine-basis error can predict a better handoff. TASK (measurement first, a policy only if the evidence supports it): (1) offline experiment on Mac small cells (si8, si16, si32, si-skew, mgo8 and synthetic-PP fixtures: scripts/fast_screen.py, dft-cli on benchmark inputs; several seeds): capture coarse states at several candidate stopping points (use the existing snapshot capture BAANDI_SNAPSHOT_DIR and the new run-config), prolong the orbitals to the SAME final cutoff and grid, evaluate the residual component outside the coarse orbital basis (kinetic-inverse-weighted indicator: sum over G of |r(G)|^2/(G^2/2 + shift)) and calibrate against an accurately completed fine solve (the energy error or the number of fine maps/H columns needed to reach the standard class after handoff); report for each candidate handoff point the total cost coarse + handoff + all fine + verification in H columns and operator applications (not just coarse iterations) so that the optimum of the current fixed rule can be compared with an indicator-based rule; (2) cost the probe itself (one fine-grid H block at a sparse trigger) and say whether an indicator can beat the fixed rule on total cost; (3) ONLY if the indicator-based rule wins on at least 4 different systems (cubic, skewed, elongated/slab-like, non-Si PP) implement it as BAANDI_HANDOFF_INDICATOR=off|trace|on (default off; 'trace' records the indicator without changing behavior; add the switch to crates/dft-numerics/src/switches.rs and update the switch-count tests crates/dft-numerics/tests/all/switches.rs and crates/dft-cli/tests/all/run.rs); final cutoff, grid, PP, final tolerance and verification stay unchanged; (4) also record how the coarse stage splits its time on the Mac cases (solve/density/potential/energy, maps, tolerance schedule) and whether a cheaper coarse configuration (lower coarse cutoff fraction, coarser FFT grid, looser inner tolerance) reduces the total without more fine work: report any such candidate with counts, but do not change defaults without evidence on >= 4 systems. Success may be a validated indicator and a CLOSED experiment (negative result): state it plainly. Use TMPDIR_990 under /Volumes/990-2TB/Develop/rust-dft-scratch (never /tmp).
OWNED_FILES:
  - crates/dft-engine/src/cutoff_ladder.rs
  - crates/dft-engine/src/session.rs   # handoff indicator trace only
  - crates/dft-engine/src/snapshot.rs
  - crates/dft-numerics/src/switches.rs
  - crates/dft-numerics/tests/all/switches.rs
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/design/ladder-handoff.md   # new
  - docs/performance/attempts-ledger.md
READ_ONLY_CONTEXT: AGENTS.md, tickets/issues/ISSUE-23.md, ledger rows LADDER-*, BT3-COARSE, D1-LADDER-*, SCF-TRAJECTORY, docs/performance/plan-2026-11-18-digest.md
IMPLEMENTATION_STEPS:
  1. Run the offline experiment (script under scripts/, results table in docs/design/ladder-handoff.md).
  2. Decide go/no-go with the pre-registered rule above; implement only on go.
TEST_COMMANDS:
  - scripts/check.sh fast dft-engine dft-numerics dft-cli
ACCEPTANCE_CRITERIA: table of total cost per candidate handoff on >= 5 systems; explicit go/no-go; default behavior bit-identical; any switch default off.
NON_GOALS: removing the ladder, changing final-stage tolerances or verification.
GIT: leave changes uncommitted; Claude commits after review.
RETURN_FORMAT: docs/schemas/report.schema.json
