TASK_ID: G7B-UPF-FULL-SUPPORT
GOAL: |
  dev08 gen26 showed that G7's rule (truncate at cutoff_radius_index, accept only print-rounding tails) still rejects the real files: SG15 Si_ONCV_PBE-1.0 PP_BETA.3 has 6 nonzero points after index 360 (a smooth continuation through a sign change, 4.0e-6 relative 3.2e-6) and O.pbe-hgh PP_BETA.1 has 125 analytic Gaussian tail points after index 773 (2.8e-13 and decaying). Architect decision: do not truncate and do not use a tolerance. Integrate each projector over its full stored support, i.e. up to its last nonzero sample (or the declared index if that is larger), so the file's data is used exactly as given; the declared index is only recorded. Reject only non-finite samples. This is the rule that needs no threshold and cannot lose projector norm.
READ FIRST: crates/dft-upf/src/parser.rs (G7 tail logic), docs/upf-support.md (G7 section), crates/dft-upf/tests/all/beta_tail.rs, crates/dft-engine/tests/all/oncv_beta_tail.rs; how the projector radial support length is used by the nonlocal setup (crates/dft-physics).
TASK: |
  1. Replace the G7 truncation/tolerance with the full-support rule; record per projector the declared index, the effective support index and the largest tail beyond the declared index (inspect-upf --json shows them).
  2. Existing fixtures and files without tails: bit-identical (bench_scf both baselines, reference_identity).
  3. Tests: synthetic ONCV-like (smooth continuation beyond the declared index) and HGH-like (Gaussian tail) files accepted; projector integrals equal a direct full-mesh quadrature to 1e-14; non-finite rejected. Update docs/upf-support.md (supersede the G7 tolerance text, keep the history line).
OWNED_FILES:
  - crates/dft-upf/
  - crates/dft-physics/src/   # only if the support length must be taken from the new field
  - crates/dft-cli/src/   # inspect-upf JSON fields only
  - docs/upf-support.md
  - crates/dft-engine/tests/all/oncv_beta_tail.rs
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-upf
  - scripts/check.sh (CHECK_SCCACHE=0, BAANDI_TEST_TMP outside the checkout)
  - python3 -I scripts/bench_scf.py --require-identical (both baselines)
  - cargo fmt --check && cargo clippy --workspace --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: ONCV/HGH-style tails accepted by the full-support rule; recorded support metadata; existing results bit-identical.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
