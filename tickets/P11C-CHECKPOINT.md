TASK_ID: P11C-CHECKPOINT
GOAL: Versioned checkpoint/restart with three explicit kinds: exact continuation (bit-identical continuation for the same build and threads), warm start (compatible initial guess, different settings allowed), result export. Directive item A4; the payload is shared with ESLE state (ADR 0016).
OWNED_FILES: crates/dft-engine/src/checkpoint.rs (new), crates/dft-cli/src/checkpoint_io.rs (new; file format lives in the CLI layer), crates/dft-cli/src/main.rs (subcommands only), tests under both crates
DEPENDENCIES: P11A, P11B merged
IMPLEMENTATION_STEPS: serialize mixer history, density, previous eigenvectors, iteration counters, options, policy log; compatibility checks (cell, atoms, species, pseudopotential identity hash, cutoff, k mesh, electron count, XC, spin, schema version) with named rejection reasons; large arrays in a binary section (little-endian f64, documented layout) with JSON header; checksum.
TEST_COMMANDS: scripts/check.sh; interrupt at iteration k, restore, finish: bit-identical to uninterrupted (threads 1 and 4); warm start from a different mixing converges to the same energy within 1e-10 Ha; incompatible metadata rejected for each field.
ACCEPTANCE_CRITERIA: as above; corrupted checkpoint rejected; file format documented.
NON_GOALS: remote storage.
GIT: branch task/P11C-CHECKPOINT, uncommitted
RETURN_FORMAT: docs/schemas/report.schema.json
