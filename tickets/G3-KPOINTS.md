TASK_ID: G3-KPOINTS
GOAL: |
  Generalization stage G3 (docs/roadmap/reference-generalization.md): k-point meshes and explicit k lists in the reference SCF driver. Today the driver requires real Gamma (driver.rs: "reference SCF requires real Gamma, standard class, no symmetry..."). The legacy path supports Gamma, Monkhorst-Pack meshes (shifted or not) and k lists with time-reversal reduction, validated against QE on Si.
READ FIRST: docs/roadmap/reference-generalization.md; docs/design/reference-scf.md (whole: map, Broyden, ethr rule, verification, kept HX, stop rule); crates/dft-solvers/src/ref_davidson.rs (real Gamma) and the legacy complex distributed Davidson (crates/dft-solvers/src/distributed_davidson/, complex path) for reusable pieces; the legacy k-point SCF (weights, time-reversal pairs, density sum over k, Fermi level not needed here: fixed occupations); docs/research/qe-study/davidson-iteration-structure-gamma.md (the complex version differs only in Hermitian algebra); docs/performance/AGENT-BRIEF.md (hot-path rules). Never read QE source.
TASK: |
  Stage it; each stage compiles, passes tests and is documented:
  1. A complex reference Davidson (same algorithm, rules and switches as ref_davidson, Hermitian products, complex preconditioner), either generic over the scalar or a sibling module sharing the workspace code; Gamma real stays the special case and bit-identical.
  2. Driver over a k list: per-k orbital blocks and kept state, per-k Davidson with the same ethr rule (QE applies one ethr to all k), density = weighted sum over k (time-reversal weights as the legacy path), dr2/Broyden/verification unchanged (they act on the density); records per k.
  3. Distribution: k points processed sequentially on all ranks (G distributed as now) in this stage; k-parallel pools are a later optimization (note it).
  4. Unsupported combinations stay rejected with clear messages (spin, Fermi-Dirac: G4; symmetry reduction: G6).
  5. Tests (Mac): Si2 and a skewed cell with 2x2x2 shifted and unshifted meshes and a k list: reference vs legacy delivered energy within the standard class and 1e-10 Ha tight, eigenvalues per k, ranks 1/3/8; Gamma-only results bit-identical to main.
OWNED_FILES:
  - crates/dft-solvers/src/ref_davidson.rs
  - crates/dft-solvers/src/ref_davidson/
  - crates/dft-solvers/src/ref_davidson_complex.rs   # new, if a sibling module is chosen
  - crates/dft-solvers/src/lib.rs   # module declaration only
  - crates/dft-engine/src/refscf/
  - crates/dft-engine/src/refstart.rs   # k-aware start only
  - crates/dft-physics/src/   # only complex distributed operator pieces the legacy path lacks
  - crates/dft-engine/tests/refscf.rs   # new tests only
  - crates/dft-cli/src/run.rs   # reference-driver k wiring and records only
  - docs/design/reference-scf.md   # a "k points" section only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-solvers
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-cli reference
  - cargo fmt --check && cargo clippy -p dft-solvers -p dft-engine -p dft-physics -p dft-cli --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: k meshes and lists run in the reference driver with legacy agreement; Gamma path bit-identical; unsupported combos rejected; design documented.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
