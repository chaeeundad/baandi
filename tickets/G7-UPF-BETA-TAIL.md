TASK_ID: G7-UPF-BETA-TAIL
GOAL: |
  dev08 gen25: the UPF reader rejects SG15 ONCV (Si_ONCV_PBE-1.0/1.1) and HGH (O.pbe-hgh) files in preflight with "invalid field PP_BETA: non-finite samples or nonzero samples beyond cutoff". These files store beta projectors on the full radial mesh; samples beyond the projector cutoff (cutoff_radius_index / the PP_BETA cutoff_radius attribute) are tiny nonzero tails (rounding) rather than physical support. ONCV (SG15, PseudoDojo older sets) and HGH are widely used NC families, so the generality rule requires reading them.
READ FIRST: crates/dft-upf/src/ (PP_BETA parsing and the validation that produces this message), docs/upf-support.md, how the nonlocal projector Fourier transforms use the radial support (crates/dft-physics nonlocal setup), the UPF v2 format attributes cutoff_radius_index and cutoff_radius (public format documentation only; never copy QE source).
TASK: |
  1. Define the rule from the file's own metadata: when PP_BETA declares cutoff_radius_index (or cutoff_radius), integrate the projector only up to that index (truncate the tail), and accept tails beyond it if their magnitude is below a stated relative tolerance of the projector's maximum (derive the tolerance from the file's printed precision, not from a material); still reject non-finite samples and real support beyond the cutoff (tails above the tolerance), with a message naming the projector, the index and the largest tail value. Document in docs/upf-support.md.
  2. Both drivers use the same parsed projector (single reader).
  3. Tests: a synthetic UPF with full-mesh beta and 1e-12-level tails (accepted, truncated, projector integrals unchanged within 1e-12 vs the clean file), one with a real tail (rejected with the precise message), and the existing UPF fixtures bit-identical; if Si_ONCV_PBE-1.0.upf exists under ../pseudo (DFT_ONCV_UPF env), an ignored test that parses it and runs a tiny Si2 reference SCF against the strict oracle (standard class).
OWNED_FILES:
  - crates/dft-upf/
  - crates/dft-physics/src/   # only if the projector radial integration must honour the cutoff index
  - docs/upf-support.md
  - crates/dft-engine/tests/   # new ignored real-file test only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-upf
  - scripts/check.sh (CHECK_SCCACHE=0, BAANDI_TEST_TMP outside the checkout)
  - python3 -I scripts/bench_scf.py --require-identical (both baselines)
  - cargo fmt --check && cargo clippy --workspace --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: ONCV/HGH-style full-mesh projectors accepted by a documented metadata-based rule; real tails rejected precisely; existing fixtures and baselines bit-identical.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
