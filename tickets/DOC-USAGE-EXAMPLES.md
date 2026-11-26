TASK_ID: DOC-USAGE-EXAMPLES
GOAL: |
  User request (2026-11-26): add usage guidance to README and a set of runnable examples. Users should be able to go from "cargo build" to each supported workflow in minutes, with inputs that pass `dft-cli validate` on the current single-path engine (ADR 0047).
READ FIRST: README.md (current structure; sections "Development Status", "What BAANDI can do today"; the status and validation blocks are generated: do not hand-edit them); examples/ (README.md, si_lda*.toml); docs/getting-started.md; docs/cli.md; docs/user/*.md; docs/decisions/0047-reference-driver-single-path.md; crates/baandi-py (Python API) and crates/baandi-mcp (MCP tools); scripts/gen_synthetic_upf.py; docs/upf-support.md (which pseudopotential files are supported; never commit UPF files).
TASK: |
  1. examples/: one directory per workflow, each with input.toml, a short README.md (what it computes, the command, expected runtime on a laptop, what to look at in the output, the UPF file it needs and where to obtain it: PseudoDojo/QE pseudo library/SG15 by name and URL, never committed), and nothing generated. Workflows:
     - 01-si-gamma-lda (total energy, eigenvalues)
     - 02-si-kmesh (Monkhorst-Pack, symmetry on)
     - 03-si-pbe (PBE; PseudoDojo NLCC note)
     - 04-bands-dos (bands path and DOS, output formats)
     - 05-al-metal-fd (Fermi-Dirac occupations)
     - 06-fe-spin (collinear spin with fixed moment; note G4D: free-moment PBE uses the strict oracle)
     - 07-forces-relax (displaced Si8 forces/stress and ion relaxation)
     - 08-cell-relax (variable cell)
     - 09-checkpoint-restart (checkpoint at a map, restart; output is an export container)
     - 10-observe (run --observe event stream; agent mode pointer)
     - 11-python (baandi-py script: run, forces, stepwise)
     - 12-mpi (mpirun with --features mpi; same input; rank-count notes)
     - 13-strict-oracle (BAANDI_REFERENCE_STRICT=on for validation)
     Keep the existing si_lda*.toml files (move into 01/02 if equivalent, keeping old paths working via the examples/README index or leave them in place).
  2. A test (crates/dft-cli/tests/examples_validate.rs) that runs `dft-cli validate` on every examples/*/input.toml with a synthetic UPF substituted (scripts/gen_synthetic_upf.py or the existing test fixture helper) so the examples cannot rot; mark heavy ones as validate-only, never run SCF in debug tests.
  3. README.md: add a "Usage" section after "Development Status" (above the generated status block): build (cargo build --release -p dft-cli; features mpi/mkl), the smallest run command and what files it writes, pseudopotentials (where to get them; PSEUDO_DIR), a table linking each examples/ directory with one line, Python and MCP pointers, and where to read more (docs/getting-started.md, docs/cli.md, docs/user/). Keep it short (about 40-60 lines). Do not edit the generated blocks; run make_readme_figures and its test.
  4. docs/getting-started.md and examples/README.md: point to the new examples; examples/README.md is the index.
OWNED_FILES:
  - examples/
  - README.md   # new Usage section only; generated blocks untouched
  - docs/getting-started.md
  - crates/dft-cli/tests/examples_validate.rs   # new
  - scripts/gen_synthetic_upf.py   # only if a helper is needed
TEST_COMMANDS:
  - cargo test -p dft-cli --test examples_validate
  - python3 -I scripts/make_readme_figures.py && python3 -I scripts/tests/test_readme_figures.py
  - cargo fmt --check && cargo clippy -p dft-cli --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: 13 example directories that validate; README Usage section; generated blocks fresh.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
