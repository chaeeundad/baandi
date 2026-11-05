# GPU execution architecture (proposal)

Date: 2026-11-03. Task: P10E-GPU-ARCH. Repository inspected: `0142ebf`.
Design only; no backend, dependency, default or scientific result changes.
See [ecosystem and sources](../research/rust-gpu-ecosystem.md),
[implementation gates](../roadmap/gpu.md) and [ADR proposal 0020](../decisions/0020-gpu-execution-boundary.md).

## 1. Recommendation and evidence boundary

Keep the pure Rust CPU implementation as the default and reference. Design for
persistent backend-owned blocks, then test FP64 CUDA on an approved Linux GPU
before selecting a production accelerator. Apple Silicon remains a first-class
CPU development platform; a Metal experiment is a distinct precision research
track, not an interchangeable FP64 implementation. Neither CUDA nor Metal is
mandatory. No universal crossover atom count or GPU speedup has been measured.
Small Si cases remain on CPU pending transfer-inclusive evidence.

All repository facts below refer to this checkout, not concurrent task branches.
P10B (MKL), P10C (fused FFT), P10A (RMM-DIIS), P10D (mixing) and P10G (general
benchmarks) have separate tickets/worktrees. Their pending changes are not
assumed implemented. ADR 0019 is an accepted design; the inspected numerical
implementation is still rustfft/faer. The session registry snapshot contains no
P10 entries; this is not evidence that the other worktrees are idle.

## 2. Verified CPU inventory

Paths are relative to the repository root; named symbols are searchable anchors.

| Data / operation | Current storage and data flow | GPU suitability / host responsibility |
| --- | --- | --- |
| Orbitals | `crates/dft-numerics/src/linalg.rs::ComplexMatrix`: owned column-major `Vec<Complex64>`, Npw rows by bands. Public columns are host slices. `scf.rs::KPointSolution` owns per-k `EigenSolution`; `run_impl` retains warm starts. | Keep each active k block resident; export only on explicit result/snapshot requests. Current host-slice API cannot provide residency by swapping a library. |
| Davidson basis and products | `crates/dft-solvers/src/lib.rs::solve_with_restart`: basis V, fixed-capacity cached HV, correction output blocks for each possible width. Computes V†HV, reduced eigensolve, VC and (HV)C, residuals, preconditioning, orthonormalization, append/restart. Fresh H application certifies recurrence residuals. | GEMM, Gram matrices, vector updates, norms and orthonormalization suit large blocks. Keep V/HV/residuals resident. Host chooses active roots, restart and convergence using small summaries. |
| Dense numerics | `linalg.rs::gemm`, `hermitian_eigh` convert owned matrices to/from faer; calls use `Par::Seq`. `orthonormalize` performs two projection passes and dependent-column filtering. | Do not promise current Davidson is allocation-free. New workspace APIs must preallocate and preserve rank/conditioning checks. Small reduced EVD can stay on CPU. |
| Local Hamiltonian | `crates/dft-physics/src/hamiltonian.rs::HamiltonianData::apply_with`: scatter sphere coefficients; pruned inverse FFT; multiply by real potential; pruned forward FFT; gather plus kinetic diagonal; nonlocal add. `Hamiltonian` owns immutable data and a `RefCell` pool of per-worker FFT/grid/projection/output scratch. | High-value resident pipeline. Physics defines ordering, units and coefficients; execution adapter owns device operations. No host round trip between its stages. |
| Nonlocal operator | Same file, `NonlocalProjectors::new_bloch`: precomputed per-atom complex columns, real coupling matrices and k-dependent phases; `add` computes overlaps and reconstruction in loops. | Upload once per immutable basis/geometry. Backend may compute B†C and B(D(B†C)) with GEMM, respecting distinct atom/channel blocks. No change to projector physics. |
| FFTs | `crates/dft-numerics/src/lib.rs::Fft3d`, `FftPlan3d`, `SpherePattern`: C-order Array3, last index contiguous; reusable rustfft plans, scratch and lane buffer; pruned patterns. | 3D or batched/stick FFT behind backend. Plans keyed by grid, precision, batch, layout, pattern, context. Scratch is never shared by overlapping executions. |
| Density | `crates/dft-physics/src/density.rs::from_orbitals` scatters each orbital, full inverse FFT, adds f times modulus squared divided by volume. `scf.rs::step_kpoints` copies columns into vectors. `band_density.rs::BandDensity` reuses worker grids and per-band real contributions for threaded path. | Resident inverse FFT and density accumulation; bounded band tiles rather than storing every band's real-space density. Preserve fixed band/k reduction order. |
| Potentials | `scf.rs::ScfProblem` owns ionic reciprocal and real arrays; `step_kpoints` evaluates Hartree/XC from host Array3 density, forms effective potential, clones into each Hamiltonian. | Initially CPU Hartree/XC and mixing; transfer one updated real potential and one output density per SCF step. Later optional resident field operations, unchanged G=0/XC contracts. |
| Audit and results | `crates/dft-engine/src/energy.rs::evaluate_with_functional` receives orbital vectors, input/output density and energy terms; `scf.rs` aggregates weighted audits, runs fixed-point certificate, builds host `ScfResult`. | Resident kinetic/nonlocal expectation reductions must replace audit's orbital downloads before claiming full GPU SCF. Host owns energy ledger, acceptance, occupations/gap checks, history and CLI-facing output. |
| Threads | `scf.rs::ordered_solves`, `ordered_chunks`; outer=min(max(threads,1),nk), inner=max(1,max(threads,1)/outer). Scoped workers own mutable k states; Hamiltonian non-Sync. Density folds ascending band then k; errors selected in order after joining. | Backend queues are not one stream per OS thread. Preserve CPU budget behavior; accelerator owns a separately bounded stream/workspace pool. |

### Normalization is an invariant, not a vendor default

Current `Fft3d` states `to_real = sum_G c_G exp(+iG.r)` (unscaled),
`to_reciprocal = sum_r f_r exp(-iG.r)/M`. With normalized sphere coefficients,
physical psi includes `1/sqrt(Omega)`; density adds `f |sum c exp|²/Omega`.
The earlier ADR 0004 spike measured a library convention (inverse /M), not this
current public contract. Match current signs, axis order, k-specific sphere
indices and Hartree G=0 semantics. cuFFT is unnormalized in both directions;
the adapter supplies the forward 1/M factor [S4]. Do not assume Gamma implies
real coefficients: this engine uses complex Gamma orbitals.

### What the existing measurements establish

[Optimization log](../performance/optimization-log.md): the initial 6,445-sample
profile attributed 48.0% to FFT and about 30% to allocation. Post-O1 FFT was
75.4%; those are historical profiles, not today's universal fractions. The
P7P7 supplied profile reports FFT about 54%, orthonormalization 8%, hypot 12%,
GEMM 9%; ADR 0019 also cites FFT/GEMM as the large remaining kernels. This
supports FFT/Hamiltonian first, not a GPU reduced eigensolver first. Amdahl's
limit for accelerating only a 54% region is 1/(1-0.54)=2.17x even with zero
kernel time and zero transfers; it is a model, not a forecast.

P7PB's shared-Mac record has Gamma wall 0.140629 s at one thread and 0.094070 s
at eight, k444s 3.697095 s at one thread. Its synthetic 8-atom-shaped fixture
has 2,969 plane waves, 20 bands and a 36³ grid: H block times 7.561308 ms
(one worker) and 1.662129 ms (eight). These are not GPU measurements. Its real
si8 SCF fails at 400 Davidson iterations, residual 1.90890249000622e-4 Ha;
time-to-failure is not a converged speedup. `baseline-after-o2.json` and the
other baseline files describe distinct historical revisions; do not mix their
timings into a new before/after comparison.

Rejected CPU lane packing and band batching in P7P9/P7PB do not justify repeating
those CPU experiments. A GPU full-batch FFT is a different hypothesis (vendor
FFT throughput versus launch/transfer costs) and must be measured independently.

## 3. Backend boundary and lifecycle

Extend the ADR 0001 boundary in `dft-numerics`, with private optional FFI adapters
following ADR 0019. No crate is added by this proposal. A future adapter crate
requires its own ownership/license approval. Physics, solvers and the SCF
controller see only backend-neutral resource types; no CUDA pointers, Metal
objects, Vulkan handles or library enums cross the boundary.

Use several capability traits under one execution session, not a giant
calculator trait. The following is a contract inventory, not a Rust ABI:

| Capability | Required contract |
| --- | --- |
| Session / capabilities | Explicit construction for a compiled backend and device ordinal; FP64 complex support, FFT shapes, memory limits, reproducibility mode. Reject unsupported combinations before allocating. CPU remains default. No silent downcast or mid-run CPU fallback. |
| Storage | Opaque owned RealGrid / ComplexBlock / IndexMap resources with shape, strides, logical length, precision, context id and version. Host storage is a distinct owned value; no `Deref` or implicit conversion to host slices. Checked read-only / exclusive mutable subviews; prohibit aliasing outputs with inputs unless an operation explicitly supports it. |
| Transfers | Explicit upload/download into preallocated staging, returning completion tokens. Pinned host storage is adapter-owned and bounded. Host buffers cannot be mutated or released while DMA is outstanding. Mapping unified memory still requires completion and exclusive access. |
| FFT and fused local operator | Prepared plans own scratch; shape and sphere-pattern identity checked. Execute unscaled inverse and normalized forward, local-potential apply, weighted density accumulation; fixed buffer layouts and explicit overwrite/add semantics. |
| Dense algebra | Complex FP64 GEMM with identity/adjoint, Gram, column selection/rotation, dot/norm, axpy, preconditioner and rank-aware orthonormalization; workspace query/reserve before the loop. Small Hermitian EVD reports convergence/status and supports an explicit host staging path. |
| Scheduling | Queue submission, event dependency, wait and finish; resources remain alive until their last use completes. Separate plan/workspace/BLAS handle per concurrently active queue. No global device synchronization for every kernel. |
| Failure | Typed Unsupported, Shape, ContextMismatch, OutOfMemory, NonFinite, LibraryFailure, DeviceLost, NotConverged. Attach operation/k/band identifiers outside loops. Completion may fail after enqueue succeeds. Poison affected outputs and invalidate caches; never publish partial data as converged. |

Public methods can synchronously return a neutral handle whose completion is
tracked privately; advanced asynchronous entry points live at the numerics
boundary. Trait-associated device types remain private implementation details.
The CPU implementation can complete immediately. A host-only mock with delayed
completion is sufficient to test ownership/error sequencing before hardware.

Lifecycle: validate metadata -> query/reserve total workspace -> upload immutable
basis/projectors and initial C -> update potential version -> enqueue solve ->
read small residual/eigenvalue summaries -> accumulate density/audit -> wait at
SCF decision boundary -> repeat -> explicit host result snapshot -> drain and
release session. Potential updates cannot race with applications; HV caches
carry potential/basis generations and are invalidated at each changed operator.
OOM is a planning failure: reduce explicitly configured k/band concurrency
before beginning a solve, or return required/available bytes. Never page large
Davidson blocks implicitly and call the result a resident algorithm.

### Relationship to P10C and MKL

P10C defines semantic operations: apply local V to sphere coefficients and
accumulate weighted squared modulus into density, with a grid/pattern plan and
owned scratch. P10C is now merged (`crates/dft-numerics/src/fused.rs`, status 2026-11-05).
Use its current CPU contracts without duplicating its implementation. A host-slice fused interface is useful for CPU/MKL and testing,
but cannot alone keep GPU C/HV resident. A later additive resident-view overload
or companion capability must use the same normalization tests. The GPU may use
full cuFFT grids where CPU uses pruned sticks; zero outside the specified output
sphere and exact input-support assumptions remain part of the contract.

ADR 0019 selects CPU libraries at compile time. Follow that policy: optional,
explicit accelerator feature builds, user-supplied SDKs, no runtime plugin
loader. Device ordinal/stream count within that compiled backend is runtime
configuration, not library switching. cudarc's default dynamic loading must be
disabled if chosen; select its documented dynamic-linking mode [S1]. CPU/MKL
features continue their existing semantics. Feature combinations need explicit
compile-time diagnostics and CI; `--all-features` must not become a mandatory
GPU build on the Mac. Do not repurpose CPU `threads` as a GPU launch dimension.

### State and coupling semantics

ADR 0016 and [ElectronicState](electronic-state.md) define physical model,
basis/order, convergence and provenance, not an address-space format. Archive
arrays are synchronized host snapshots; device buffers and queues are never
serialized as state. Backend/library/device/precision/reduction mode belong to
execution provenance, not the physical model fingerprint. Equal physical models
do not imply byte-identical state payloads or exact trajectory replay. Compare
densities and occupied subspaces (phase/degenerate-rotation invariant), not raw
orbital equality. Mixed precision needs a distinct accuracy qualification.

[Scientific coupling](scientific-coupling.md) and ADR 0017 keep field role, units,
origin, environment identity and energy ownership independent of execution.
External field upload is explicit at a request boundary. GPU errors do not
become zero observables; unavailable results and failed inner convergence remain
typed. Core execution adds no serialization, networking or file formats.

## 4. Precision and reproducibility

Default CPU arithmetic and thread-count bit identity are preserved. Across
CPU/MKL/GPU, reduction association and FFT implementation differ: numerical
equivalence is required, bit identity is not promised (ADR 0019). GPU promotion
still requires the existing benchmark energy delta <=1e-12 Ha and unchanged
SCF/QE gates. Failure is a failed gate, not permission to loosen it.

For a selected GPU backend, require repeatability on a fixed device architecture,
library/compiler version, launch geometry and deterministic reduction plan;
record all of them. Ordered k/band tiles and fixed reduction trees must be
independent of host worker count. Avoid floating atomic density accumulation;
use bounded partial arrays and an ordered fold. Disable reduced-precision math
and automatic mixed modes for the FP64 path. cuBLAS documents reproducibility
conditions and multi-stream workspace caveats [S5]; separate stream workspaces
are necessary but not proof of whole-SCF determinism. If stream-count independence
cannot be certified, restrict the validated mode to one fixed stream configuration
and mark broader support experimental. Do not advertise the CPU's stronger
thread contract for an untested accelerator. Metal mixed precision, software
FP64 and refinement are separate research, not an automatic fallback.

## 5. Transfer and capacity model

These are calculated planning examples, not implemented support for 500 atoms.
One k point, complex FP64 (16 bytes), real FP64 (8 bytes). Atom count alone is
insufficient: Npw scales approximately as volume times Ecut^(3/2); occupied bands
as valence electrons/2; vacuum and harder projectors raise memory sharply.
The family below extrapolates the measured synthetic si8 shape: Npw=ceil(2969 A/8),
b=ceil(2.5 A), L=ceil(36 (A/8)^(1/3)), M=L³. The 500-atom grid is illustrative,
not the engine's automatically selected smooth grid. No material generality is
claimed from this sizing family.

W=16 Npw b (one orbital block), D=8 M (one density OR potential).
Use proposed Davidson capacity m=3b solely for budgeting, not a solver default.
A planning subtotal includes V/HV (6W), four b-wide blocks (4W), two complex FFT
grids for each of q=8 band slots (32Mq), four real grids (4D), P=4A projector
columns (16 Npw P), and two reduced matrices (32m²). Vendor FFT/EVD workspace,
allocator padding, staging, CPU mirrors, mixer histories and runtime allocations
are additional and must be queried. Thus the subtotal is not peak memory.

| Atoms A | Npw | b | Grid | W MiB | D MiB | Subtotal GiB | C upload + HC download at 25 GB/s, ms | One density down + potential up, ms |
| ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 8 | 2,969 | 20 | 36³ | 0.906 | 0.356 | 0.023 | 0.076 | 0.030 |
| 64 | 23,752 | 160 | 72³ | 57.988 | 2.848 | 0.764 | 4.864 | 0.239 |
| 216 | 80,163 | 540 | 108³ | 660.523 | 9.611 | 7.899 | 55.409 | 0.806 |
| 500 | 185,563 | 1,250 | 143³ | 3,539.333 | 22.310 | 41.297 | 296.901 | 1.871 |

MiB=2^20 bytes, GiB=2^30 bytes, GB/s=10^9 bytes/s. Times are 2W/B or 2D/B
with no overlap/latency. The 25 GB/s sustained one-way link is an assumption,
not a measured machine. Multiple k points multiply per-k blocks/projectors;
shared potential/density need not be duplicated. At 500 atoms four simultaneous
k states can exceed 160 GiB before library workspaces: use bounded active-k
scheduling. Doubling valence/bands doubles W and quadruples reduced matrices;
doubling cutoff roughly multiplies Npw and M by 2.83 at fixed volume. Slabs
increase M with vacuum without increasing bands. P10G must cover these cases.

| Link / memory example | Published or assumed bandwidth | Correct interpretation |
| --- | --- | --- |
| PCIe 4 x16 | ~31.5 GB/s theoretical each direction; assume 25 GB/s useful | Derived from 16 GT/s, 128/130 encoding and 16 lanes; actual pinned transfers must be measured. |
| PCIe 5 x16 | ~63 GB/s theoretical each direction; assume 50 GB/s useful | Published H100 figure 128 GB/s is rounded bidirectional aggregate [S15]; do not use it for one-way copies. |
| H100 SXM NVLink | 900 GB/s aggregate bidirectional [S15] | GPU-peer topology capability, not a normal CPU-to-GPU PCIe rate. Assuming 450 GB/s for one direction only gives an optimistic bound. |
| Apple M4 Pro unified memory | 273 GB/s published memory bandwidth [S16] | Shared CPU/GPU memory traffic, not a copy link and not this host's measured bandwidth. No PCIe DMA for shared allocations; still synchronization, cache and bandwidth costs. |

For the 500-atom W round trip, idealized 50/450 GB/s gives 148.451/16.495 ms.
The latter only applies to an appropriate peer path, not this workstation.
Reading and writing W once consumes 2W of shared-memory traffic: at the M4 Pro
peak that is at least 27.189 ms; sharing an allocation avoids a separate copy,
not the algorithm's memory traffic. None of these figures predicts FP64 Metal
performance.

### Crossover is a measured surface, not an atom-count constant

For one region, GPU wins only if
`Tcpu > Tgpu_compute + bytes/B + Nlaunch*latency + Tsync + Tplan/R`.
R is the reuse count. Include initial/final uploads and all SCF certificates and
audits in whole-run timing. For a resident solve the large W exchange is amortized
across iterations; only small reduced matrices and convergence summaries cross
inside Davidson. At m=60 the reduced complex matrix is 57,600 bytes; at m=3750
it is 225,000,000 bytes, so host EVD is no longer automatically cheap.

A transparent sensitivity example: assume GPU compute is 4x faster for the region,
12 launches at 10 microseconds each, 20 microseconds synchronization, and no
planning cost after warmup. Break-even CPU-region time is
`(transfer + 0.140 ms)/(1-1/4)`. Using the W round trips in the table yields
0.288, 6.672, 74.066 and 396.055 ms for 8/64/216/500 atoms. These assumed kernel
speed/latency values are unverified; the small region often cannot achieve 4x.
This is a testable model, not a size threshold. The si8-shaped H block could
benefit if these assumptions hold, so asserting that *every* 8-atom workload
loses would be unsupported. The tiny primitive Si benchmarks are CPU-only
production candidates because their total wall time is already ~0.1 s and no
transfer-inclusive win exists. Use 8/64/216/500 as measurement bins, not dispatch
thresholds. Only repeated end-to-end measurements can establish a crossover.

Keep C, V, HV, corrections, kinetic, projectors, sphere maps, FFT scratch, current
Veff and density partials resident. Downloading all orbitals each audit defeats
this design. Initially transfer density and Veff once per SCF step; later move
Hartree/XC/mixing only when measured. Final exported orbitals are a deliberate,
accounted transfer required by the legacy result API.

## 6. Review and verification scope

Acceptance mapping: A is section 2; B is sections 3–4; C is section 5; D is the
linked ecosystem assessment; E/F are the linked roadmap. All four deliverables
are proposals. CPU checks are recorded in the task report/roadmap, not claimed
as GPU tests. No source files, manifests, README evidence, artifact JSON or
other tickets are edited: the explicit OWNED_FILES restriction takes precedence
over the general dashboard/README update workflow. Dashboard publication and
session-registry reconciliation belong to the orchestrator after review.

Source identifiers S1–S18 refer to the dated bibliography in the ecosystem
assessment. No source implementation from QE or another project was used.
