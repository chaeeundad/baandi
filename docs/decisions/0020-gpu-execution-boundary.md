# ADR 0020: Optional GPU execution behind resident numerical resources

- Date: 2026-11-03
- Task: P10E-GPU-ARCH
- Status: proposed, design only; no implementation or hardware authorization
- Number: 0020 was the next free number in this worktree at `0142ebf`.
  Concurrent P10 branches may propose the same number; reviewer must reconcile
  numbering and links before merge rather than overwrite another ADR.

## Context

The current CPU implementation keeps orbitals in host `ComplexMatrix` values
and grids in `Array3`. Replacing only FFT/GEMM libraries cannot preserve device
residency across Davidson, density and energy auditing. Stored profiles identify
FFT and then block algebra as useful targets; small Si runs do not establish a
GPU crossover. Apple development must remain supported without proprietary SDKs.

ADR 0001 isolates numerical backends; ADR 0004 prefers maintained libraries;
ADR 0019 approves optional compile-time MKL selection. ADRs 0016 and 0017 keep
state/coupling semantics independent of execution and persistence outside kernels.
This proposal extends those boundaries without superseding their CPU contracts.

## Proposed decision

1. Keep pure Rust CPU default, APIs and numerical results unchanged. Accelerator
   support requires explicit optional compile-time features and a separate
   implementation ticket; no mandatory CUDA, Metal, Vulkan, HIP or SYCL dependency.
2. Put memory/plan/queue/library ownership behind `dft-numerics` and approved
   private FFI adapters. Physics and SCF see neutral blocks/grids and semantic
   operations, never vendor device types. No new crate now.
3. Prepare persistent C/V/HV, residual/projector buffers and bounded workspaces.
   Explicit transfers, completion/error tracking, resource lifetimes and potential
   generations are part of the contract. No implicit per-operation host copies,
   paging fallback, silent precision conversion or runtime backend plugin loader.
4. Share P10C's local-potential and density semantics with CPU/MKL; add resident
   views only through a later additive API proposal after its final signature is
   known. Existing host-slice interfaces remain available.
5. Preserve current FFT normalization and physical conventions. FP64 is the
   production target; cross-backend comparisons use numerical tolerances, including
   the existing <=1e-12 Ha benchmark-energy gate. CPU bit identity remains intact.
   GPU deterministic execution requires separately measured fixed reduction and
   library/stream configurations; no cross-device bit identity is promised.
6. Prefer an optional CUDA FP64 feasibility experiment on an approved existing
   Linux machine because documented FFT/BLAS/EVD coverage is available. This is
   a study priority, not library adoption. Apple CPU stays supported; Metal mixed
   precision is a separately approved research track. HIP, Vulkan/wgpu and SYCL
   remain candidates, with the limitations recorded in the ecosystem document.
7. Physical state/coupling identity excludes device handles and keeps units,
   environment, energy ownership and convergence semantics. Snapshot arrays only
   after synchronization; record execution backend/precision in provenance. Warm
   continuation is not an exact-replay promise.

## Alternatives and reasons

- GPU FFT behind unchanged host-array calls only: useful first oracle, but repeated
  transfers limit performance; insufficient architecture for full SCF.
- Expose CUDA/Metal buffers through physics/controller interfaces: rejected,
  creates backend coupling and breaks state/coupling neutrality.
- Mandatory portable wgpu backend: rejected as a production choice now; portable
  FP64 and the required numerical library stack are not established.
- Mandatory Metal parity: deferred; native double precision is not established,
  and no full accuracy-qualified numerical path is available in this review.
- Write custom FFT/GEMM/EVD: rejected absent an approved library gap study;
  increases maintenance and numerical validation burden.
- Multi-GPU/MPI first: deferred until resident single-device evidence exists.
- Runtime-loaded optional plugins: not adopted; follow ADR 0019 compile-time
  policy and explicit user-provided SDK linkage.

## Consequences and validation

Detailed [architecture and calculated capacity/transfer model](../architecture/gpu-backend.md),
[dated ecosystem/source review](../research/rust-gpu-ecosystem.md), and
[roadmap with numerical/hardware/timing gates](../roadmap/gpu.md) accompany this ADR.
The model estimates 0.906 MiB to 3,539.333 MiB per orbital block over an illustrative
8–500 atom family; this is not measured performance or universal material scaling.
No hardware crossover, GPU precision or GPU speedup has been verified.

This document-only task runs the unchanged CPU `scripts/check.sh`; all four gates
passed, exit 0. Future promotion requires diverse P10G fixtures, CPU numerical
oracles, existing QE agreements, true residual and fixed-point certification,
repeatability, and transfer-inclusive release timing. A failed gate preserves CPU
and does not change a tolerance. The optional hardware prototype is skipped:
local Metal enumeration is not a ready FP64 numerical stack.

Pending user decisions for later work: machine identity/access, whether Metal is
a research or production target, and proprietary SDK distribution/build policy.
No decision is needed to review these documents. No code, manifests, dependencies,
shared artifacts or accepted ADRs changed; no superseded ADR. Leave uncommitted
for reviewer reconciliation and acceptance.
