# dft-upf

**Purpose:** UPF v2 norm-conserving reading and validation layer, radial-to-reciprocal transforms. Never copy QE (GPL) code.

**Invariants to preserve and test:** Unsupported files fail with explicit errors; no silent field skipping; Hermitian D_ij; finite, monotonic radial mesh.

**Before editing:** read ../../AGENTS.md, ../../ownership.toml (allowed dependencies) and the ADRs in ../../docs/decisions/.
**Verify:** `scripts/check.sh dft-upf`, then `scripts/check.sh` for the whole workspace.
**Rules:** keep public API changes out of unrelated tickets; every behavior change adds or updates a test with stated tolerance and units; prefer general-purpose crates (ADR 0003); record design decisions in an ADR.
