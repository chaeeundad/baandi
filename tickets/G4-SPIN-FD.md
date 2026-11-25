TASK_ID: G4-SPIN-FD
GOAL: |
  Generalization stage G4 (docs/roadmap/reference-generalization.md): Fermi-Dirac occupations and collinear spin (fixed sector and FD, with or without a fixed moment) in the reference SCF driver, Gamma and k meshes (G3 is merged: refscf/k_driver.rs, k_map.rs, ref_davidson_complex.rs). Today the driver rejects spin and FD.
READ FIRST: docs/design/onepath3-spin-metals.md (the approved estimator for two spin densities, FD occupations on the standard path, acceptance bounds); docs/design/reference-scf.md (whole, incl. the new "K points" section); crates/dft-engine/src/finite_temperature.rs, collinear.rs and crates/dft-physics/src/occupations.rs (legacy oracles: Fermi level solve, entropy term, Mermin free energy, per-channel number solves); docs/research/qe-study/ notes on mixing if any. Never read QE source.
TASK: |
  Stage it; each stage compiles, passes tests and is documented:
  1. FD occupations (spin-unpolarized): per map after all k solves, solve the Fermi level from the band energies (legacy occupations module), weight the density by f_ik w_k, report the Mermin free energy and entropy term as the legacy path does; the band count must leave empty bands above E_F as in legacy (same default rule). The Davidson convergence set covers the bands with non-negligible occupation exactly as the legacy FD path defines it. Stop estimator and verification unchanged in meaning (onepath3 section on FD).
  2. Collinear spin: two density channels; per channel the k loop and kept HX state; V_s = V_H[n] + V_xc,s[n_up, n_down] + V_loc (LDA-PZ/PW spin from the legacy evaluators; spin PBE only if the distributed GGA already supports two channels, otherwise reject it clearly); Broyden on the stacked (n_up, n_down) residual (or (n, m), whichever the onepath3 design fixes) with the Coulomb metric of the design; the onepath3 spin estimator for the stop; verification on both channels. Fixed sector (integer occupation per channel) and FD with free or fixed moment.
  3. Records: per channel per k eigenvalues and occupations, E_F (per channel when the moment is fixed), magnetization total and absolute, entropy term.
  4. Rejections kept: noncollinear, other smearing kinds not in legacy.
  5. Tests (Mac): FD on a small metal-like synthetic fixture and Si2 with a gap (FD must reproduce the integer result within the class), Gamma and a 2x2x2 mesh; spin: fixed sector on a synthetic open-shell atom, FD free moment and fixed moment; reference vs legacy delivered (free) energy within the standard class and 1e-10 Ha tight, magnetization and E_F agreement, ranks 1/3/8; existing closed-shell Gamma results bit-identical to main.
OWNED_FILES:
  - crates/dft-engine/src/refscf/
  - crates/dft-engine/src/refstart.rs   # spin/FD-aware start only
  - crates/dft-solvers/src/ref_davidson_complex.rs   # occupation-aware convergence set only
  - crates/dft-solvers/src/ref_davidson.rs   # occupation-aware convergence set only, real Gamma unchanged when occupations are integer
  - crates/dft-physics/src/distributed_complex.rs   # spin channel plumbing only
  - crates/dft-physics/src/distributed_real_gamma.rs   # spin channel plumbing only
  - crates/dft-engine/tests/refscf.rs   # new tests appended at the end only
  - crates/dft-cli/src/run.rs   # reference-driver spin/FD wiring and records only
  - docs/design/reference-scf.md   # a "Spin and Fermi-Dirac" section only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-solvers
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-cli reference
  - cargo fmt --check && cargo clippy -p dft-solvers -p dft-engine -p dft-physics -p dft-cli --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: FD and collinear spin run in the reference driver with legacy agreement; closed-shell paths bit-identical; unsupported combinations rejected; design documented.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
