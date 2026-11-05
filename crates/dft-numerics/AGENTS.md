# dft-numerics

## Current ownership (2026-11-05)

Owns FFT adapters (`fused.rs`, optional `mkl_backend.rs`), linear algebra and
`tiled.rs` workspaces, and the persistent `pool.rs` worker executor. Pure-Rust
reductions use fixed pieces independent of worker count; MKL remains optional,
with sequential library calls inside the engine budget. Analysis scripts/counters
under `pool/` are diagnostics, not production telemetry. P10M retained copy and
duplicate-product removal; its CholeskyQR2 trial was rejected. See the performance
README and optimization log before modifying these paths.


**Purpose:** Adapters to general-purpose numerical libraries (linear algebra, FFT). Thin traits so backends stay swappable; unsafe is confined to adapters here and in dft-mkl under ownership.toml.

**Invariants to preserve and test:** Column-major complex layout stated; FFT normalization stated; plans and scratch reused; unsafe preconditions documented and tested.

**Before editing:** read ../../AGENTS.md, ../../ownership.toml (allowed dependencies) and the ADRs in ../../docs/decisions/.
**Verify:** `scripts/check.sh dft-numerics`, then `scripts/check.sh` for the whole workspace.
**Rules:** keep public API changes out of unrelated tickets; every behavior change adds or updates a test with stated tolerance and units; prefer general-purpose crates (ADR 0003); record design decisions in an ADR.
