# dft-solvers

**Purpose:** Eigensolvers (dense baseline, then Davidson/LOBPCG), mixers (linear, then Pulay), convergence checks.

**Invariants to preserve and test:** Residual and orthonormality reported; convergence needs density, energy and eigen residual, never energy alone.

**Before editing:** read ../../AGENTS.md, ../../ownership.toml (allowed dependencies) and the ADRs in ../../docs/decisions/.
**Verify:** `scripts/check.sh dft-solvers`, then `scripts/check.sh` for the whole workspace.
**Rules:** keep public API changes out of unrelated tickets; every behavior change adds or updates a test with stated tolerance and units; prefer general-purpose crates (ADR 0003); record design decisions in an ADR.
