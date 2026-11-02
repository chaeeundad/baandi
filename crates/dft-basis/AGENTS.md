# dft-basis

**Purpose:** Plane-wave basis: G-vectors with 0.5*|G|^2 <= Ecut, FFT grid choice, G-index <-> FFT-index mapping.

**Invariants to preserve and test:** Negative-frequency wrap; cutoff boundary; round trip; Parseval; aliasing policy documented.

**Before editing:** read ../../AGENTS.md, ../../ownership.toml (allowed dependencies) and the ADRs in ../../docs/decisions/.
**Verify:** `scripts/check.sh dft-basis`, then `scripts/check.sh` for the whole workspace.
**Rules:** keep public API changes out of unrelated tickets; every behavior change adds or updates a test with stated tolerance and units; prefer general-purpose crates (ADR 0003); record design decisions in an ADR.
