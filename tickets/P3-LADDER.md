TASK_ID: P3-LADDER
GOAL: |
  Bring the C1 cutoff ladder to the reference SCF driver (switch REFSCF_LADDER, default off). The reference driver's trajectory now matches QE map by map (si128 15 maps to QE's stop, si64 11), and its per-map cost is within 5-10% of QE, so the remaining lever to beat QE is doing fewer full-cutoff maps. On the legacy path the ladder was adopted (ledger C1-LADDER-DEV08): coarse stage at E1 = fraction * E2 on the smaller grid, exact prolongation of density and orbitals into the full basis, then the full-cutoff SCF; si32 strict T=1 at realistic displacements -28% to -36% (full-cutoff SCF 23-29 -> 11-17), zero correctness violations in 160 runs; the reference driver bypasses it today (issue 25 patch 2).
READ FIRST: docs/roadmap/c1-cutoff-ladder.md; crates/dft-engine/src/cutoff_ladder.rs (prolongation, exactness argument, secant seeding S2); ledger rows C1-DESIGN, C1-LADDER, C1-LADDER-DEV08, ALGEBRA-STREAM-B1-MAC; docs/design/reference-scf.md (start, ethr rule, stopping rule, verification); docs/performance/AGENT-BRIEF.md. Never read QE source.
TASK: |
  1. With REFSCF_LADDER on: run the reference driver at the coarse cutoff (same driver, same rules, a coarse stopping threshold chosen by a documented general rule, e.g. a fixed ratio above the fine threshold or a residual target tied to the cutoff ratio, never tuned to Si), prolong the density and orbitals exactly (reuse cutoff_ladder's injection/Fourier padding; Gamma real layout, owned rows), start the fine reference SCF from them with an empty Broyden history (optionally the S2 secant seeding behind a second switch), and keep the fine ethr rule, stopping estimator and verification exactly as now. The delivered result, energy audit and verification are the fine stage's.
  2. Records: mark coarse and fine maps separately (stage field), charge coarse time to the run's wall, and keep refscf-check --stop-to-verified working on the fine stage. Choose the coarse fraction from the legacy design (documented); do not tune it on Si.
  3. Tests (Mac only: cargo check, touched crates' tests): fixtures cubic/skew/elongated/synthetic nonlocal at simulated ranks 1/3/8: converge and verify, delivered energy within the standard class of the strict reference and within 1e-10 Ha of the non-ladder run where both are tight, fine-map counts per arm reported, prolongation exactness test (charge, orthonormality, Gamma symmetry bitwise).
  4. docs/performance/lanes/RSCF/ladder-dev08-request.md: ONE run per arm (D vs L) on si64/si128/elongated/rskew MPI8, verify on and off, plus si32 and si64 with a realistic displacement (fractional +-0.001 and +-0.01 as in C1-LADDER-DEV08); report wall, coarse/fine maps, coarse stage seconds, E - E_strict.
OWNED_FILES:
  - crates/dft-engine/src/refscf/
  - crates/dft-engine/src/refstart.rs   # only to accept a prolonged start
  - crates/dft-engine/src/cutoff_ladder.rs   # reuse of the prolongation helpers only; legacy behaviour unchanged
  - crates/dft-engine/tests/refscf.rs   # new tests only
  - crates/dft-cli/src/run.rs   # reference ladder wiring and records only
  - crates/dft-cli/src/refscf_check.rs   # stage field only
  - crates/dft-numerics/src/switches.rs   # append at the end
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/design/reference-scf.md   # a "Cutoff ladder" section only
  - docs/performance/lanes/RSCF/ladder-dev08-request.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-engine cutoff_ladder
  - cargo test -p dft-cli reference
  - cargo test -p dft-numerics switches
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-cli -p dft-numerics --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: REFSCF_LADDER default off; fine stage delivers and verifies with unchanged rules; exact prolongation tests; counts per arm; dev08 one-run request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
