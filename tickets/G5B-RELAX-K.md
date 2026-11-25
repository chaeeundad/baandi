TASK_ID: G5B-RELAX-K
GOAL: |
  Close the remaining gaps found while generalizing the reference driver (docs/roadmap/reference-generalization.md "Open gaps"; G5 and G6 reports):
  (a) relaxation with k meshes: warm start of each ionic step from the previous step's per-k orbitals and density in the reference k driver (G5 did Gamma only);
  (b) CLI scf.verify = "off" for the reference driver, with forces/stress labelled "unverified" (engine side exists, config.rs wiring missing);
  (c) stress with symmetry: forward the symmetry to stress.rs so the optional Pulay cutoff correction works on symmetric runs (G6 returned the frozen-basis tensor only);
  (d) `dft-cli refscf-check --stop-to-verified` must read the k/FD/spin results schema (dev08 gen6: "missing field map_start_s");
  (e) the binary hash cache test race (dev08 gen7, crates/dft-cli/src/run.rs binary_hash_cache): a same-length rewrite within one coarse ctime tick reports a stale "hit". Fix the cache the racy-git way: when the file's mtime or ctime is within a safety window of the time the entry was recorded (or of now), treat the entry as unverified and rehash; make the test deterministic.
READ FIRST: docs/design/reference-scf.md (K points, Forces/stress/relaxation, Symmetry, Spin and Fermi-Dirac); crates/dft-cli/src/relax.rs; crates/dft-engine/src/refscf/k_driver.rs and delivered.rs; crates/dft-engine/src/stress.rs; crates/dft-cli/src/config.rs; crates/dft-cli/src/refscf_check.rs; run.rs binary_hash_cache. Never read QE source.
TASK: |
  Implement (a)-(e). Tests (Mac): 3-step Si2 relaxation on a 2x2x2 mesh reference vs legacy (tightened thresholds as in G5: dE <= 5e-9 Ha, dF <= 1e-6 Ha/Bohr, dr <= 1e-6 Bohr), with warm transfer recorded per step; verify=off CLI run writes "unverified" labels and the same delivered state as verify on minus the verification map; symmetric Si2 stress with Pulay correction vs no-symmetry full mesh (1e-8 Ha/Bohr^3 tight); refscf-check on k/FD/spin reference outputs; cache test deterministic (run it 50 times in a loop in the test or with a forced ctime collision).
OWNED_FILES:
  - crates/dft-cli/src/relax.rs
  - crates/dft-cli/src/config.rs   # scf.verify for the reference driver only
  - crates/dft-cli/src/refscf_check.rs
  - crates/dft-cli/src/run.rs   # reference-driver relax/verify wiring and binary_hash_cache only
  - crates/dft-engine/src/refscf/k_driver.rs   # warm start seed only
  - crates/dft-engine/src/refscf/delivered.rs
  - crates/dft-engine/src/refscf/electronic.rs   # warm start seed only
  - crates/dft-engine/src/stress.rs   # symmetry forwarding only
  - crates/dft-engine/tests/refscf_forces.rs   # new tests appended
  - crates/dft-cli/tests/   # new tests only
  - docs/design/reference-scf.md   # "Forces, stress and relaxation" section only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf_forces
  - cargo test -p dft-cli relax
  - cargo test -p dft-cli binary_hash_cache
  - cargo test -p dft-cli stop_to_verified
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-cli --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: (a)-(e) done with the stated tests passing; Gamma and existing paths unchanged.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
