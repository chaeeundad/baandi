# dft-physics

**Purpose:** Kinetic, local ionic, nonlocal KB, Hartree, LDA XC, density; matrix-free H*psi.

**Invariants to preserve and test:** Hartree G=0 convention; XC v = d(n*eps)/dn finite-difference test; H Hermitian; density integrates to N_valence.

**Before editing:** read ../../AGENTS.md, ../../ownership.toml (allowed dependencies) and the ADRs in ../../docs/decisions/.
**Verify:** `scripts/check.sh dft-physics`, then `scripts/check.sh` for the whole workspace.
**Rules:** keep public API changes out of unrelated tickets; every behavior change adds or updates a test with stated tolerance and units; prefer general-purpose crates (ADR 0003); record design decisions in an ADR.
