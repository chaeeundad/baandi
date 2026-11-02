# dft-numerics

**Purpose:** Adapters to general-purpose numerical libraries (linear algebra, FFT). Thin traits so backends stay swappable; the only crate allowed to contain unsafe.

**Invariants to preserve and test:** Column-major complex layout stated; FFT normalization stated; plans and scratch reused; unsafe preconditions documented and tested.

**Before editing:** read ../../AGENTS.md, ../../ownership.toml (allowed dependencies) and the ADRs in ../../docs/decisions/.
**Verify:** `scripts/check.sh dft-numerics`, then `scripts/check.sh` for the whole workspace.
**Rules:** keep public API changes out of unrelated tickets; every behavior change adds or updates a test with stated tolerance and units; prefer general-purpose crates (ADR 0003); record design decisions in an ADR.
