TASK_ID: R1B-ROUTING-AND-SCREEN
GOAL: |
  Two follow-ups from dev08 gen13/gen14 after the reference driver became the default (R2):
  A. Routing: an input that explicitly requests driver = "reference" and also contains a legacy-only key currently falls back to legacy silently (or, before the L9 fix 5d787b70, failed with a misleading "unknown mixer"). An explicit request must never be silently overridden: fail with an error naming the legacy-only key(s) and the reason. Fallback stays allowed only when the driver came from the default or the environment. Also `dft-cli explain` must print the driver selection before it runs validation (today it dies in preflight without printing a selection).
  B. XC screening default: README QE case PBE Gamma (Si2, PBE override on Si.pz-vbc, ecut 30, 25^3, 8 bands) fails on the default reference path because a Broyden map transiently produces a density of -5.6e-6 Bohr^-3 at one grid point (the atom site); with XC_SCREEN_NEGATIVE=on it converges in 7 iterations, VERIFIED, 2 points screened in total (<= 1e-7 e), E - E_legacy = -6.5e-11 Ha, E - E_QE unchanged from the README value. Adopt the ADR 0046 rule as the default: flip XC_SCREEN_NEGATIVE to "on", mark ADR 0046 accepted with this evidence, and keep the rejection behaviour reachable with the switch off.
READ FIRST: crates/dft-cli/src/driver_select.rs, config.rs (parse-time mixing adapter around preliminary()), inspect.rs (explain); docs/decisions/0046-xc-small-and-negative-density.md; crates/dft-physics/src/xc.rs (switch read); tests that assert rejection of negative densities (crates/dft-engine/tests/refscf.rs around lines 6179 and 6693 and any others: grep "nonnegative", "XC failed", "negative total density").
TASK: |
  A. In driver_select: when selection.source == Input and a LEGACY_ONLY predicate matches, return an error path (e.g. Selection with an error reason that config validation turns into InvalidInput naming the key(s)); add routing fixtures for explicit reference + each legacy-only key; explain prints "driver: <driver> (source, reason)" first and then validation errors.
  B. Flip XC_SCREEN_NEGATIVE default to on (switch table and default tests); pin every test that asserts rejection of negative densities to xc_screen_negative=off; ADR 0046 status accepted with the dev08 gen14 evidence above (cite /data/sig-dft/gen14/RESULT.md); docs/design/reference-scf.md vacuum paragraph updated. Positive-density results must stay bit-identical (reference_identity, bench_scf --require-identical against baseline-reference.json and baseline-oracle.json).
OWNED_FILES:
  - crates/dft-cli/src/driver_select.rs
  - crates/dft-cli/src/config.rs   # routing error surfacing only
  - crates/dft-cli/src/inspect.rs
  - crates/dft-cli/tests/driver_routing.rs
  - crates/dft-cli/tests/fixtures/routing/
  - crates/dft-numerics/src/switches.rs   # XC_SCREEN_NEGATIVE entry only
  - crates/dft-numerics/tests/all/switches.rs
  - crates/dft-engine/tests/   # pins only
  - crates/dft-physics/tests/   # pins only
  - crates/dft-cli/tests/   # pins only
  - docs/decisions/0046-xc-small-and-negative-density.md
  - docs/design/reference-scf.md   # vacuum densities paragraph only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - scripts/check.sh (BAANDI_TEST_TMP outside the checkout)
  - python3 -I scripts/bench_scf.py --require-identical (release; reference baseline default, and with BAANDI_REFERENCE_STRICT=on against baseline-oracle.json)
  - cargo fmt --check && cargo clippy --workspace --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: explicit reference requests never silently fall back; explain shows the selection; screening on by default with ADR accepted; baselines bit-identical.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
