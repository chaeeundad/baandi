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

## P26A real Gamma foundation

`real.rs` owns column-major real matrices and borrowed faer views, reusable
256-row tiled AB / A-transpose-B products, and sequential symmetric EVD.
Gram partials are summed in ascending panel order at every pool budget; real
GEMM dispatch uses 1,048,576 multiply-adds (the complex guard's real-flop budget).
This heuristic is not a speed measurement. Construction owns scratch allocation;
`multiply_into` does not allocate numerical storage. EVD allocates and fixes each
vector's sign by its first largest-magnitude component; degenerate eigenspaces
have no unique basis. Real kernels use faer even under the optional MKL feature.
This is library-only groundwork from real-gamma-design sections 3-5/11, not an
adopted engine optimization or a repeat of rejected complex-kernel experiments.
Within each 256-row panel, fixed 16-by-16-by-16 backend microtiles keep faer
0.24 on its allocation-free nano-GEMM path; large generic GEMM may allocate
packing buffers. Reduction microtiles run in ascending order. This is an R1
storage-contract choice, not a revival of rejected complex-path panel tuning;
no engine arithmetic or performance default changes. Recheck backend allocation
behavior when upgrading faer.

R1 local validation (2026-11-06): `cargo test -p dft-numerics real::tests --
--nocapture` covers 0/1/255/256/257/1025 rows, including dispatch above the
parallel threshold, with pool budgets 1/2/4/8 and nested resident calls.
All output bits agree across budgets. Maximum absolute naive-GEMM difference:
1.1102230246251565e-15; complex-on-real tiled-GEMM comparisons pass 1e-14.
EVD sizes 0/1/2/9/24 have maximum entrywise residual 1.0658141036401503e-14
and orthogonality defect 1.9984014443252818e-15 (both gates 1e-13). These are
synthetic dimensionless matrix tests, not material energies or speed evidence.

SCF isolation was checked against a release build of parent
`c99eb2683e9d92846394bf1702b658492d5ac964`, using `scripts/bench_scf.py`
with `--cases gamma,k222s,k222,k444s,si8 --repeats 1 --require-identical`.
All five energy differences were 0 Ha and density SHA-256 values matched;
iterations were 10/10/10/10/15. UPF SHA-256:
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
Parent binary SHA-256:
`cdda20ce5dcd5ea2e348dd72de88591bd3cb334f4588922fc3e555a9d8e6557d`;
final R1 binary SHA-256:
`ac818152ad0e04b4963afc02a11a899801fe0729f2cc147db6023f796d61705b`.
This is Mac correctness evidence only; no QE run or performance acceptance.
