# Rust GPU ecosystem feasibility assessment

Date/access date: 2026-11-03. Task: P10E-GPU-ARCH. Documentation and package
metadata review only; no third-party implementation source was inspected or
copied, no package was fetched into the workspace, and no GPU library was run.
Version observations below are the versions displayed by upstream documentation
on this date, not proposed pins. Mutable `latest` pages must be rechecked before
an implementation ADR. Publication year is stated where known; `n.d.` means the
page did not establish a publication year (2026 is the access year).

## 1. Decision matrix

The engine needs complex FP64, 3D FFT, GEMM, Hermitian reduced solves and explicit
memory lifetime control. A safe dispatch wrapper alone is not a numerical stack.

| Candidate | Rust integration and numerical coverage | Apple Silicon / Linux | Feasibility and maintenance evidence |
| --- | --- | --- | --- |
| CUDA driver + cudarc | Context, stream, device slices/views; cuBLAS, cuSOLVER and cuFFT modules; safe/result/raw layers. cuFFT module documents complex-double 3D/batched transforms. | NVIDIA Linux target; no Apple GPU execution. Documentation compiling for a host does not mean hardware support. | Best first FP64 experiment if a machine is approved. cudarc 0.19.10 dated 2026-09-24 and several 2026 releases show active publication [S1]. Exact workspace controls and solver routine coverage need a pinned-version compile/test, not inference from module names. |
| CUDA driver + cust | Driver context/streams/events/device memory, raw interop; numerical libraries need additional bindings/adapters. | NVIDIA Linux; no Apple GPU. | 0.3.2 is the published version observed [S2]. No current release cadence or complete cuFFT/cuBLAS/cuSOLVER wrapper suite verified. More adapter work than cudarc; do not infer abandonment solely from version. |
| CUDA runtime bindings | `cuda-runtime-sys` exposes low-level runtime entry points; version 0.3.0-alpha.1, zero documented items shown [S3]. | NVIDIA runtime/driver required on supported Linux; no Apple GPU. | Raw bindings do not own asynchronous lifetimes. Current-toolkit ABI coverage and maintenance not verified. Prefer one driver-based context owner over mixing runtime and driver ownership accidentally. |
| Apple Metal | `objc2-metal` framework bindings; MPS/MPSGraph for selected operations; custom shaders otherwise. | Native Apple platform; no Linux Metal runtime. | `metal` 0.33.0 explicitly deprecated in favor of objc2; `objc2-metal` 0.3.2 dated 2025-10-04 [S7,S8]. FP64 GPU path is not established; Metal precision limits block a drop-in backend. |
| wgpu / Vulkan | wgpu resource/queue validation and compute; Vulkan via ash gives lower-level control. Neither supplies this whole scientific numerical stack. | wgpu has Metal on Apple and Vulkan on Linux. Native FP64 feature is Vulkan/SPIR-V only in the reviewed docs. | wgpu 30.0.1 dated 2026-08-22 is active; ash observed 0.38.0+1.3.281 [S11,S12]. Portability of dispatch does not imply portability of FP64 FFT or eigensolvers. |
| SYCL / oneMath | C++ heterogeneous model, buffers/USM and events; oneMath lists BLAS, DFT and LAPACK including Hermitian EVD. Rust would need a small C ABI bridge and exception conversion. | Linux vendor toolchains; AdaptiveCpp documents an Apple Metal backend, with limitations. | Active documentation, but complete maintained Rust coverage not verified. Adds C++ compiler/runtime/backend matrix. Keep as an alternative for a concrete multi-vendor user [S13]. |
| ROCm / HIP | HIP runtime, hipFFT to rocFFT or cuFFT; rocBLAS/hipBLAS and rocSOLVER/hipSOLVER family for linear algebra. Rust needs pinned bindings plus ownership wrappers. | Supported AMD Linux hardware; NVIDIA route depends on selected HIP libraries/toolchain; not an Apple Metal route. | AMD publishes current library docs (hipFFT 1.0.26, rocSOLVER 3.37.0 observed). Complete maintained safe Rust stack was not verified. Viable second FP64 candidate if AMD hardware is the available target [S14]. |

## 2. CUDA: adopt libraries, own the boundary

The driver API explicitly owns contexts/modules/streams; the runtime API is a
higher-level CUDA interface commonly used by C/C++ and numerical libraries.
Choose a single context policy, verify library handle interoperation and keep
raw handles private. cudarc provides documented safe wrappers for driver buffers
and stream submission; raw calls and custom kernel launches still require
adapter-side safety proofs [S1]. A Rust lifetime on a buffer is insufficient
unless last-use completion is tracked. Guard host staging, stream ownership and
plan destruction through completion tokens.

cuFFT supplies complex-double transforms, reusable plans and workspace controls;
normalize explicitly [S4]. cuBLAS supplies complex GEMM with column-major layout
and documents multi-stream reproducibility restrictions [S5]. cuSOLVER supplies
Hermitian dense solves and workspace/status APIs [S6]; distinguish enqueue success
from the numerical `info` result. Test ABI layout of complex scalars, leading
dimensions, checked integer conversions and workspace sizes. Do not send each
small Gram matrix through a fresh handle or plan. Small reduced EVD on CPU is a
valid explicit capability until measurements justify device EVD.

Reject managed-memory demand paging as a substitute for residency analysis.
Unified virtual addressing does not imply uniform access cost. Preallocate device
and pinned staging pools; use copies/events explicitly. Consumer NVIDIA GPUs may
have weak FP64 throughput despite strong FP32/AI throughput: query the exact
hardware rather than select by product family or tensor peak numbers.

Build implication: cudarc defaults to runtime dynamic loading [S1], conflicting
with the project's compile-time selection policy. A future experiment should
explicitly enable dynamic linking and a fixed CUDA version feature; no auto-SDK
download, default loader, hidden fallback or toolkit detection in CPU builds.
No custom Rust GPU compiler is required to call cuFFT/cuBLAS; custom physics
glue kernels still need a separately reviewed compilation strategy (build-time
SDK compilation or explicit NVRTC, whose cold cost and packaging must be counted).

## 3. Apple: separate CPU acceleration from GPU precision research

Accelerate/vDSP offers CPU single- and double-precision FFT APIs [S9]. It is a
possible CPU backend, not evidence of GPU FFT support. Three-dimensional use
would require validated axis/batch handling, layout and normalization; Accelerate
BLAS/LAPACK can also support a separately ticketed CPU study. Keep the pure Rust
path available without adding any requirement.

Metal command buffers/queues and shared buffers fit persistent ownership. However,
shared physical memory still needs synchronization and correct access lifetimes.
`objc2-metal` explicitly warns about unsafe shader calls, bounds, synchronization
and resource retention [S8]. Prefer that maintained direction over deprecated
`metal-rs`; do not call all framework methods safe because Rust bindings exist.

MPSGraph does expose FFT operations, confirmed by Apple API documentation and
WWDC24 [S10]. It is incorrect to claim Apple has no GPU FFT option. What was not
verified is a complete Rust-accessible complex-FP64 3D FFT/GEMM/EVD pipeline that
meets this engine's accuracy gates. The official MSL PDF could not be fetched by
the web tool (size limit); the implementation project's own AdaptiveCpp Metal
backend documentation explicitly states lack of native double precision [S13].
Treat native FP64 Metal support as unavailable for this proposal, and recheck
an approved device's published capabilities before implementation. The lack of
a usable `metal` compiler in this environment is a separate tooling observation,
not a claim that every Metal operation requires offline shader compilation.

Do not substitute f32 for f64 to satisfy portability. Mixed precision FFTs plus
FP64 CPU refinement, software double emulation and MPSGraph experiments require
separate numerical and timing evidence on skewed/vacuum/hard-projector systems.
They cannot inherit the 1e-12 Ha energy gate by assertion. Near-term recommendation:
CPU on Apple; optional precision experiment only if the user makes Metal a target.

## 4. Portable alternatives and missing pieces

wgpu is useful for resource management and safe compute dispatch across native
APIs [S11]. Its documented `SHADER_F64` capability is native Vulkan/SPIR-V only;
portable WGSL/Metal execution cannot be assumed to provide it. Buffer sizes,
workgroup limits and mapping rules must be queried. Implementing a new FFT,
GEMM and stable eigensolver from scratch would violate the general-purpose-library
policy and create a large verification burden. A portable FFT package could be
surveyed in a separate ticket; none is selected or license-cleared here.

ash exposes Vulkan with much more explicit unsafe synchronization/resource
management [S12]. Apple Vulkan translation still inherits the underlying Metal
precision constraints; it does not manufacture native double support. Neither
API makes GPU math bit-identical across vendors.

SYCL is a C++ standard, not itself a Rust library or an assurance of identical
backends. oneMath lists the relevant numerical domains, but coverage varies by
backend/version. Its Apple route must not be inferred from a SYCL logo:
AdaptiveCpp's documented Metal support has FP64 and event-cost limitations [S13].
A Rust C ABI shim would need opaque context/buffer ownership, asynchronous errors
and an exception barrier; do not expose C++ templates in engine contracts.

HIP/ROCm makes an AMD route plausible without changing the architecture. hipFFT
marshals to rocFFT or cuFFT and recommends device-visible resident data; rocSOLVER
provides LAPACK routines [S14]. Need exact hardware support, ROCm/driver matrix,
complex ABI tests and supported solver routines before selection. No assurance
of current safe Rust ownership wrappers or Apple execution is made.

## 5. License and CI implications

These are dependency-screening observations, not a completed redistribution
review. Verify licenses and transitive components at the chosen pinned versions
before updating `docs/dependencies.md` through a future ticket.

| Component | License evidence / risk | Build and CI consequence |
| --- | --- | --- |
| cudarc, cust, cuda-runtime-sys | MIT OR Apache-2.0 package metadata [S1–S3]; wrapper license does not cover NVIDIA libraries. | Optional Rust adapters only; pinned versions, no default SDK dependency. |
| CUDA driver/toolkit, cuFFT/cuBLAS/cuSOLVER | Proprietary NVIDIA SDK license and CUDA supplement [S18]. | User-provided SDK/runtime; no binaries/headers redistribution by default. Linux GPU runner needed for real tests. |
| objc2-metal | Zlib OR Apache-2.0 OR MIT package metadata [S8]. Deprecated metal wrapper must not be a fresh default. | Target-gated Apple adapter. SDK/frameworks remain Apple-supplied; check minimum OS availability. |
| Apple Metal/MPS/Accelerate frameworks | Proprietary platform frameworks; Rust binding licenses do not cover SDK redistribution. Exact redistribution terms not reviewed. | Use installed SDK/frameworks; separate Apple hardware job, no framework linking on Linux. |
| wgpu, ash | MIT OR Apache-2.0 upstream package/project metadata [S11,S12]. | Bigger optional dependency graph and shader validation/compiler versions; platform matrix still needs hardware. |
| SYCL / oneMath / implementation | SYCL is an open royalty-free standard [S13]; that is not an implementation license. Compiler, oneMath, plugins and vendor libraries need separate pinned license review; not verified completely here. | C++ toolchain and bridge add ABI/compiler CI. Do not adopt on portability claims alone. |
| HIP/ROCm libraries | hipFFT explicitly MIT; rocSOLVER's page has BSD-style two-clause AMD terms plus separate LAPACK/MAGMA notices. Full ROCm/component and selected Rust-wrapper license audit remains open [S14]. | User-provided supported ROCm stack; no mandatory download/build. NVIDIA backend retains NVIDIA terms. |
| FFTW / copied external implementations | FFTW is GPL under existing ADR 0019 policy; excluded. No copyleft source adoption proposed. | No QE/other project source copying; use published APIs and independent tests. |

Default CI must build/test on Apple Silicon and Linux with all accelerator
features off, no GPU SDK installed and offline vendored dependencies. Optional
compile jobs check supported targets and feature conflicts, but compilation is
not a GPU correctness test. Hardware jobs must fail or explicitly report skipped
coverage when the requested accelerator is missing; never silently run CPU and
label it GPU PASS. Add synchronization, delayed error, OOM, aliasing, invalid
context and device-loss tests. Long tests are ignored in debug and have release
commands. Recheck MSRV against the pinned repository toolchain when choosing a
binding; no compatibility was proved by this document.

## 6. Source register

All accessed 2026-11-03. Only API/reference docs, release/package metadata and
prose documentation are used; implementation source links are not evidence.

- **S1 — cudarc crate/API and release metadata** (2026, 0.19.10, release 2026-09-24): https://docs.rs/cudarc/0.19.10/cudarc/ ; https://docs.rs/crate/cudarc/0.19.10 ; **cufft module** (2026): https://docs.rs/cudarc/0.19.10/cudarc/cufft/index.html . Driver ownership, linking modes, library coverage. License advertised MIT/Apache-2.0; pin audit still required.
- **S2 — cust API and package metadata** (n.d., 0.3.2): https://docs.rs/cust/0.3.2/cust/ ; https://docs.rs/crate/cust/0.3.2 . Driver wrapper and CUDA development-library requirement. Current maintenance cadence unverified.
- **S3 — cuda-runtime-sys API/package metadata** (n.d., 0.3.0-alpha.1): https://docs.rs/cuda-runtime-sys/0.3.0-alpha.1/cuda_runtime_sys/ . Raw runtime bindings and MIT/Apache-2.0 metadata; current toolkit compatibility unverified.
- **S4 — NVIDIA cuFFT API Reference** (2026 observed CUDA 13.4 docs): https://docs.nvidia.com/cuda/cufft/index.html . Transform types, plans, normalization and workspaces.
- **S5 — NVIDIA cuBLAS API Reference, Data Layout and Results Reproducibility** (2026 observed CUDA 13.4 docs): https://docs.nvidia.com/cuda/cublas/index.html#results-reproducibility . Column-major and stream/workspace reproducibility conditions.
- **S6 — NVIDIA cuSOLVER API Reference** (2026 observed CUDA 13.4 docs): https://docs.nvidia.com/cuda/cusolver/index.html . Dense eigensolve, workspace and numerical status contracts.
- **S7 — metal crate deprecation notice** (n.d., 0.33.0): https://docs.rs/metal/0.33.0/metal/ . Recommends objc2/objc2-metal for new development.
- **S8 — objc2-metal, Safety considerations and package metadata** (2025 release 0.3.2, 2025-10-04): https://docs.rs/objc2-metal/0.3.2/objc2_metal/ ; https://docs.rs/crate/objc2-metal/0.3.2 . License metadata: https://crates.io/crates/objc2-metal/0.3.2 . Shader/buffer/command lifetime risks.
- **S9 — Apple vDSP.FFT** (n.d.): https://developer.apple.com/documentation/accelerate/vdsp/fft . Single/double CPU FFT API; not a GPU implementation promise.
- **S10 — Apple MPSGraph API; Accelerate machine learning with Metal, WWDC24** (API n.d.; video 2024): https://developer.apple.com/documentation/metalperformanceshadersgraph/mpsgraph ; https://developer.apple.com/videos/play/wwdc2024/10218/ . FFT existence and GPU graph usage, not complex-FP64 validation.
- **S11 — wgpu project, package releases and Features** (2026, 30.0.1 release 2026-08-22): https://wgpu.rs/ ; https://docs.rs/crate/wgpu/30.0.1 ; https://docs.rs/wgpu/30.0.1/wgpu/struct.Features.html#associatedconstant.SHADER_F64 ; https://github.com/gfx-rs/wgpu . Native Vulkan-only FP64 feature, platform coverage and project license metadata.
- **S12 — ash API and package description** (n.d., 0.38.0+1.3.281): https://docs.rs/ash/0.38.0+1.3.281/ash/ ; https://docs.rs/crate/ash/0.38.0+1.3.281 . Vulkan loader/FFI boundary; maintenance cadence not independently audited.
- **S13 — Khronos SYCL overview; UXL oneMath documentation; AdaptiveCpp Metal backend guide** (living docs, publication year unverified): https://www.khronos.org/sycl/ ; https://uxlfoundation.github.io/oneMath/ ; https://github.com/AdaptiveCpp/AdaptiveCpp/blob/develop/doc/install-metal.md . C++ execution model, numerical domains, Apple backend's documented double-precision/event limitations.
- **S14 — AMD HIP, hipFFT and rocSOLVER documentation and licenses** (2026 docs observed): https://rocm.docs.amd.com/projects/HIP/en/latest/what_is_hip.html ; https://rocm.docs.amd.com/projects/hipFFT/en/latest/ ; https://rocm.docs.amd.com/projects/hipFFT/en/latest/license.html ; https://rocm.docs.amd.com/projects/rocSOLVER/en/latest/ ; https://rocm.docs.amd.com/projects/rocSOLVER/en/latest/license.html . Runtime/library architecture; hipFFT MIT verified, full stack licensing unverified.
- **S15 — NVIDIA H100 specifications; Hopper Tuning Guide** (n.d., current docs observed 2026): https://www.nvidia.com/en-gb/data-center/h100/ ; https://docs.nvidia.com/cuda/hopper-tuning-guide/ . H100 SXM NVLink 900 GB/s bidirectional; PCIe Gen5 128 GB/s aggregate. Topology matters.
- **S16 — Apple MacBook Pro (16-inch, 2024), Tech Specs** (2024): https://support.apple.com/en-gb/121554 . M4 Pro 273 GB/s shared-memory bandwidth; illustrative specification, not this host's measurement.
- **S17 — Apple Metal Shading Language Specification** (publication year not verified): https://developer.apple.com/metal/Metal-Shading-Language-Specification.pdf . Fetch failed (web tool content-size limit); no direct quote or inspected-PDF claim. Precision assessment uses S13 and S11; this remains a verification limit.
- **S18 — NVIDIA Software License Agreement and CUDA Supplement** (n.d., current page accessed 2026): https://docs.nvidia.com/cuda/eula/index.html . Proprietary SDK conditions; no redistribution approval inferred.

## 7. Local feasibility observation

Read-only probe on 2026-11-03: `uname -m` reports arm64;
`system_profiler SPDisplaysDataType` reports an integrated Apple GPU with Metal
supported (host-provided label: Apple M6, 12 GPU cores). No performance or
capability conclusion beyond that enumeration is inferred from the label.
`command -v nvidia-smi` found nothing; no NVIDIA device is established.
`xcrun --find metal` failed with “unable to find utility”. The offline vendor
inventory contains none of cudarc/cust/metal/objc2-metal/wgpu. Hardware exists,
but a ready-to-run FP64 numerical GPU stack was not found. Metal runtime graph
experiments might be possible without that compiler; they were not proven usable
and would not satisfy the FP64 path. No install, remote access or prototype is
claimed. The ticket's optional prototype is skipped; the design is complete
without hardware results.
