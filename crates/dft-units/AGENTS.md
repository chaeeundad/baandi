# dft-units

**Purpose:** Hartree atomic units and constants. The only place for unit conversions (1 Ha = 2 Ry); newtypes over bare f64.

**Invariants to preserve and test:** Ry<->Ha and Bohr<->Angstrom round trips; no mixed-unit APIs.

**Before editing:** read ../../AGENTS.md, ../../ownership.toml (allowed dependencies) and the ADRs in ../../docs/decisions/.
**Verify:** `scripts/check.sh dft-units`, then `scripts/check.sh` for the whole workspace.
**Rules:** keep public API changes out of unrelated tickets; every behavior change adds or updates a test with stated tolerance and units; prefer general-purpose crates (ADR 0003); record design decisions in an ADR.
