TASK_ID: G4G-OATOM-STALL
GOAL: |
  dev08 gen28b: O atom (O.pbe-hgh, PBE, fixed collinear sector up 4 / down 2, isolated in a box, ADR 0046 screening on) on the default reference driver does NOT converge in 200 maps: F stays within ~1e-7 Ha of -15.6967565, the stop estimate dr2 stalls at ~1.9e-8 Ha from map ~50 (threshold 5e-9), the Broyden step fraction collapses to 1e-16..1e-12 with history condition numbers 1e5..5e10, M = 2.0 throughout. Screening is tiny (max 1.4e-10 e). The same atom with PseudoDojo O.upf (PBE+NLCC) converged in 9 maps (gen19). Real molecules/atoms in vacuum must converge; find the cause and fix it with a general rule.
READ FIRST: docs/design/reference-scf.md (Broyden mixer, stop rule, Spin and Fermi-Dirac, vacuum densities); docs/decisions/0045 (history conditioning, today applied only on the collinear FD entry) and 0046 (screening); crates/dft-engine/src/refscf/{broyden.rs,metric.rs,electronic.rs,electronic_map.rs}; dev08 report /data/sig-dft/gen28/RESULT.md is not accessible to you: the numbers above are the evidence.
TASK: |
  1. Reproduce with a synthetic analogue on the Mac: an isolated open-shell atom with a GGA (PBE) functional and no core charge, HGH-like analytic projectors (long Gaussian tails), fixed sector, in a vacuum box, until the stall pattern appears (dr2 floor, fraction collapse, ill-conditioned history). Record per map dr2, step fraction, history condition number, screened points, and whether the screening mask changes between maps (mask flicker makes the map non-smooth, which can stall quasi-Newton mixing).
  2. Candidate causes to test: (a) screening-mask flicker near the 1e-10 threshold between maps; (b) history ill-conditioning (ADR 0045 rule not applied on the fixed-sector entry); (c) GGA gradient terms in near-vacuum regions amplifying noise; (d) something else.
  3. Fix with a general rule (no atom-specific constants): e.g. apply the ADR 0045 history conditioning to every reference entry (not only collinear FD) if (b); make the screening mask smooth or hysteretic if (a); document the rule (amend ADR 0045/0046). Acceptance: the synthetic analogue converges and verifies; all baselines and the G4/G4B/G4C/G4F matrices bit-identical where the change is not active (state which cases change and why).
OWNED_FILES:
  - crates/dft-engine/src/refscf/
  - crates/dft-physics/src/   # screening hook only
  - crates/dft-engine/tests/refscf.rs   # new tests appended
  - docs/decisions/0045-broyden-history-conditioning.md
  - docs/decisions/0046-xc-small-and-negative-density.md
  - docs/design/reference-scf.md   # relevant sections only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - scripts/check.sh (CHECK_SCCACHE=0, BAANDI_TEST_TMP outside the checkout)
  - python3 -I scripts/bench_scf.py --require-identical (both baselines)
  - cargo fmt --check && cargo clippy --workspace --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: cause identified with per-map evidence; general fix; synthetic analogue converges; unaffected results bit-identical.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
