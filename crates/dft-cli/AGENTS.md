# dft-cli

**Purpose:** CLI, TOML input, output files. No numerical logic.

**Invariants to preserve and test:** run/validate/inspect-upf subcommands; nonzero exit with cause on failure; outputs carry reproducibility metadata.

**Before editing:** read ../../AGENTS.md, ../../ownership.toml (allowed dependencies) and the ADRs in ../../docs/decisions/.
**Verify:** `scripts/check.sh dft-cli`, then `scripts/check.sh` for the whole workspace.
**Rules:** keep public API changes out of unrelated tickets; every behavior change adds or updates a test with stated tolerance and units; prefer general-purpose crates (ADR 0003); record design decisions in an ADR.

## Collinear input/results (P13B / MF03)

`SpinConfig` and shared preflight validate explicit integer sectors; `spin_run`
serializes the additive version 2 result profile, proposal-schema extension,
full spin/k/band axes, fixed-sector certificate and conventional density cubes.
Unpolarized results remain version 1. Never route a collinear result through a
scalar legacy orbital view or unpolarized runtime/checkpoint/postprocessing.
Use `tests/spin.rs` and its independent standard-library schema/cube checker.
See `docs/user/collinear-spin.md` for supported combinations and format losses.
