TASK_ID: P3-ENERGY
GOAL: |
  Per-map work outside the eigensolver and density costs 1.86 s on dev08 si128 MPI8 (one run at aa6a062, verify off, 15 maps): energy 1.08, energy_correction 0.25, potential 0.29, mix 0.10, density_transform 0.07, output 0.07. QE spends ~0.27 s on everything in electrons() outside c_bands and sum_band. The reference map (crates/dft-engine/src/refscf/map.rs) evaluates a full EnergyAudit every map: kinetic energy over all rows and bands, nonlocal band energies (a second projection of all occupied bands), local and double-counting integrals, and fields(state, output density) (Hartree + XC of the output density: FFTs); the driver then computes fields(next input) for the mixed-energy correction, and the next map may compute fields of its input again.
READ FIRST: docs/performance/AGENT-BRIEF.md, attempts-ledger.md (P3 rows), docs/design/reference-scf.md (map, energy terms, "Verification", stopping rule), docs/research/qe-study/qe-scf-cycle.md (what QE evaluates per iteration: band energy, double counting from the input potential, Harris-Foulkes/variational correction with the new potential). Never read QE source.
TASK: |
  1. List every per-map field/energy computation with its cost (FFTs, projections, passes) and which result is actually consumed: by the stopping estimator, the ethr rule, the mixing, the mixed-energy correction (descf), the verification, the trace/records, or the final delivered result.
  2. Behind a switch REFSCF_LEAN_ENERGY (default off): compute per map only what is consumed (QE-style: band sum, double counting with the input fields already known, the variational energy using fields of the mixed density that the correction step computes anyway, reused as the next map's input fields instead of recomputing); compute the full audit (kinetic, nonlocal, local, output-density fields, functional difference bound) only for the delivered map and for verification maps. Records/trace keep their schema; fields not computed in lean mode are null with a documented meaning.
  3. Remove duplicate field evaluations (same density, same fields computed twice) even when the switch is off if bit-identical (prove it with a test), otherwise behind the switch.
  4. Mac checks only: cargo check --workspace --tests; cargo test -p dft-engine --test refscf; cargo test -p dft-cli reference; cargo test -p dft-numerics switches; fmt; clippy -D warnings. Delivered energies and densities identical within rounding to the default; counts unchanged on fixtures; stop-to-verified and rank-rule tests still pass.
  5. docs/performance/lanes/RSCF/energy-dev08-request.md: ONE profiled run per arm, si64/si128/elongated/rskew MPI8 (verify on and off for si128).
OWNED_FILES:
  - crates/dft-engine/src/refscf/
  - crates/dft-engine/src/energy.rs   # only if a lean constructor is needed
  - crates/dft-engine/tests/refscf.rs
  - crates/dft-engine/tests/common/rank_rule.rs   # only if record fields become optional
  - crates/dft-cli/src/run.rs   # reference-driver record serialization only
  - crates/dft-numerics/src/switches.rs   # append at the end
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/design/reference-scf.md   # a "Per-map energy cost" section only
  - docs/performance/lanes/RSCF/energy-dev08-request.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-cli reference
  - cargo test -p dft-numerics switches
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-cli -p dft-numerics --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: consumer table written; lean per-map energy behind a switch with full audit at delivery/verification; duplicate field evaluations removed; delivered results within rounding; dev08 one-run request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
