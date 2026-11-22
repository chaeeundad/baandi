TASK_ID: G-FIXTURE-SI128
GOAL: |
  The benchmark fixtures stop at 64 atoms (si8/si16/si32/si64). dev08 measurements show the BAANDI/QE ratio worsening with size (si32 np8 1.72 -> si64 np8 2.0-2.43), and the long-term target systems are much larger, so the scaling trend must be measured beyond 64 atoms. Add larger diamond-silicon Gamma fixtures that follow exactly the existing convention of benchmarks/cases/si64-gamma.toml (same potential "si", cutoff_ry 30, same lattice constant 10.261212856717932 Bohr per conventional cell edge i.e. si64 = 2x2x2 conventional cells of 8 atoms, same fields, deterministic, phase scf).
TASK: (1) create benchmarks/cases/si128-gamma.toml (2x2x4 conventional cells: 128 atoms; tier long; timeout_s generous, e.g. 7200) and benchmarks/cases/si216-gamma.toml (3x3x3: 216 atoms) and, if memory allows on dev08, note si512-gamma (4x4x4) as documented-only (do not create a runnable case unless the 64 GB scale is verified; state the estimated memory and plane-wave count); generate atoms by replicating the existing 8-atom conventional cell (verify the coordinates against si64-gamma.toml, no overlapping or missing atoms, minimum interatomic distance 4.6 Bohr = 2.43 A (diamond bond), count equals the target). (2) Compute and write down for each case: lattice, number of atoms, electrons, bands, plane waves at the fine cutoff, fine FFT grid (use the existing basis code through a small ignored test or a script that calls the CLI with a dry-run/preflight mode if one exists), estimated memory of the wavefunctions and projector tables per rank for np 1, 8, 16, and a rough operation-count scaling relative to si64. Mac check only that the files parse (the CLI reads them in preflight/dry-run mode) and that preflight succeeds; do not run SCF at 128+ atoms on the Mac. (3) Update docs/performance/benchmark-suite.md and tests that enumerate cases (if any test counts the fixtures, update the count) and write docs/performance/lanes/fixtures-si128.md with the measurement plan for dev08: np8 and np16 cases, QE comparison arms (same cell, same cutoff, 72-style FFT grid derived by the engine: the same grid in the QE input), A/A floor first, one case per job within the time cap (final- jobs 7200 s), what to return. Do not commit UPF files. Scratch inside the worktree.
OWNED_FILES:
  - benchmarks/cases/si128-gamma.toml
  - benchmarks/cases/si216-gamma.toml
  - docs/performance/benchmark-suite.md
  - docs/performance/lanes/fixtures-si128.md
  - docs/performance/attempts-ledger.md
READ_ONLY_CONTEXT: AGENTS.md, benchmarks/cases/si64-gamma.toml, benchmarks/cases/si32-gamma.toml, docs/performance/benchmark-suite.md
TEST_COMMANDS:
  - scripts/check.sh full   # fast on the Mac now (about 40 s warm); three CLI cache tests may reject the in-repo TMPDIR
ACCEPTANCE_CRITERIA: fixtures parse and preflight on the Mac; atom counts and distances verified; memory and size table written; dev08 plan written; no change to existing fixtures or code.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
