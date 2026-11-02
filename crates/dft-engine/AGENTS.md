# dft-engine

**Purpose:** Orchestrates SCF, occupations, total energy and its audit (ion-ion/Ewald, double counting).

**Invariants to preserve and test:** Every energy term has a unit test and appears in results; no hidden constants.

**Before editing:** read ../../AGENTS.md, ../../ownership.toml (allowed dependencies) and the ADRs in ../../docs/decisions/.
**Verify:** `scripts/check.sh dft-engine`, then `scripts/check.sh` for the whole workspace.
**Rules:** keep public API changes out of unrelated tickets; every behavior change adds or updates a test with stated tolerance and units; prefer general-purpose crates (ADR 0003); record design decisions in an ADR.
