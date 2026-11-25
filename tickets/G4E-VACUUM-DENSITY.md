TASK_ID: G4E-VACUUM-DENSITY
GOAL: |
  Isolated systems with vacuum and NLCC fail in both drivers (dev08 gen9, PseudoDojo O.upf PBE+NLCC, O atom in a 12 bohr box, ecut 84 Ry, grid 72^3, collinear): the per-channel XC input (valence channel + half core) is slightly negative in the vacuum (about -2e-7 Bohr^-3, Fourier truncation ringing of the core and the atomic start). Legacy fails at iteration 1 with a precise message; the reference start fails with the generic "reference start: XC failed". Molecules, slabs and surfaces with NLCC pseudopotentials (most PseudoDojo files) will hit this, so a principled, documented policy is needed.
READ FIRST: docs/decisions on density positivity / signed-density / NLCC (grep docs/decisions for "signed", "negative", "NLCC", "clipping"); docs/design/reference-scf.md (NLCC, PBE, Spin and Fermi-Dirac: "vacuum screening" already exists for distributed spin PBE in G4C); crates/dft-physics/src/xc.rs, spin.rs, gga.rs, gga/distributed.rs (current LDA/GGA policies for small and negative densities, unpolarized and polarized); crates/dft-engine/src/functional.rs evaluate_with_core; refstart.rs start XC. Never read QE source.
TASK: |
  1. Write docs/decisions/0046-xc-small-and-negative-density.md (status proposed): inventory the current rules per evaluator (unpolarized LDA-PZ/PW, PBE; polarized LDA-PW, PBE; with/without core) and where they differ; propose ONE general rule: below an absolute density threshold (derived from the evaluator's own numerical limits, stated with reasoning, not tuned to a material), the XC energy density and potential are zero for that point/channel and the point is recorded as screened; the density itself is never clipped (Hartree, kinetic and the electron count use the real density); report the number of screened points, the integrated |n| over screened points and the most negative value per map, and fail only if the screened charge exceeds a stated fraction of the total (a real valence infeasibility, not ringing). Explain the energy/potential consistency (screening must be applied identically to energy and potential so the functional derivative stays consistent) and the effect on verification.
  2. Implement behind a registry switch XC_SCREEN_NEGATIVE (default off) in both drivers (legacy XC library adapters and the reference distributed evaluators), including the polarized path with core split.
  3. Error messages: the reference start and map XC failures must name the stage (start/map/verification), the channel, the minimum value and its grid point, like legacy.
  4. Tests (Mac): synthetic isolated atom with core charge in a box at low cutoff that reproduces negative total density: off -> both drivers fail with the precise message; on -> both converge, reference vs legacy within the standard class and 1e-10 Ha tight, screened statistics recorded and identical across ranks 1/3/8; positive-density systems bit-identical with the switch on and off.
OWNED_FILES:
  - docs/decisions/0046-xc-small-and-negative-density.md   # new
  - crates/dft-physics/src/xc.rs
  - crates/dft-physics/src/spin.rs
  - crates/dft-physics/src/gga.rs
  - crates/dft-physics/src/gga/distributed.rs
  - crates/dft-engine/src/functional.rs
  - crates/dft-engine/src/refstart.rs   # XC error reporting and screening hook only
  - crates/dft-engine/src/refscf/map.rs   # XC screening records only
  - crates/dft-engine/src/refscf/electronic_map.rs   # XC screening records only
  - crates/dft-numerics/src/switches.rs   # append at the end
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - crates/dft-cli/src/run.rs   # screening record JSON only
  - crates/dft-engine/tests/refscf.rs   # new tests appended
  - crates/dft-physics/tests/   # new test file
  - docs/design/reference-scf.md   # NLCC section, a "Vacuum densities" paragraph only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-physics
  - cargo test -p dft-engine --test refscf
  - cargo fmt --check && cargo clippy -p dft-physics -p dft-engine -p dft-cli -p dft-numerics --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: ADR proposed with the inventory and rule; switch default off; isolated NLCC fixture converges in both drivers with agreement when on; positive-density paths bit-identical; precise error messages.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
