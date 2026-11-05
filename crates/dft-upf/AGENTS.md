# dft-upf

**Purpose:** UPF v2 norm-conserving reading and validation layer, radial-to-reciprocal transforms. Never copy QE (GPL) code.

**Invariants to preserve and test:** Unsupported files fail with explicit errors; no silent field skipping; Hermitian D_ij; finite, monotonic radial mesh.

**Before editing:** read ../../AGENTS.md, ../../ownership.toml (allowed dependencies) and the ADRs in ../../docs/decisions/.
**Verify:** `scripts/check.sh dft-upf`, then `scripts/check.sh` for the whole workspace.
**Rules:** keep public API changes out of unrelated tickets; every behavior change adds or updates a test with stated tolerance and units; prefer general-purpose crates (ADR 0003); record design decisions in an ADR.

P13C: generic numbered projector/wavefunction tags are normalized only after
original suffix/index/order validation; every array and D entry is retained.
Missing informational attributes have reader-only defaults; PP_R/PP_RAB alone
control radial numerics. `raw()` is the normalized compatibility model;
`header_facts` preserves original metadata and must never authorize SCF.
PP_NLCC is a density in electrons/Bohr^3, transformed by
`radial::spherical_density_transform` with the shared Simpson/3/8 quadrature.
Synthetic parser/Gaussian tests run by default. `tests/real_nc.rs` is ignored,
reads DFT_NC_DIR and records the exact third-party error plus adapter support.
See docs/cli.md for full conventions, reproduction and qualification limits.
