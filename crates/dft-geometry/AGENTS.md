# dft-geometry

**Purpose:** Cell (lattice columns a1..a3), atoms, fractional/Cartesian conversion, reciprocal lattice B = 2*pi*inverse(A)^T.

**Invariants to preserve and test:** a_i . b_j = 2*pi*delta_ij; volume = |det A|; fractional<->Cartesian round trip.

**Before editing:** read ../../AGENTS.md, ../../ownership.toml (allowed dependencies) and the ADRs in ../../docs/decisions/.
**Verify:** `scripts/check.sh dft-geometry`, then `scripts/check.sh` for the whole workspace.
**Rules:** keep public API changes out of unrelated tickets; every behavior change adds or updates a test with stated tolerance and units; prefer general-purpose crates (ADR 0003); record design decisions in an ADR.
