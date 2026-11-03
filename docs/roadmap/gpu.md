# GPU roadmap and acceptance gates (proposal)

Date: 2026-11-03. Task: P10E-GPU-ARCH. This is a sequence of future tickets,
not authorization to acquire machines, install SDKs or implement GPU code.
Architecture: [gpu-backend.md](../architecture/gpu-backend.md).
Candidate evidence: [rust-gpu-ecosystem.md](../research/rust-gpu-ecosystem.md).
Decision proposal: [ADR 0020](../decisions/0020-gpu-execution-boundary.md).

## 1. Build first, defer deliberately

First prepare neutral resource/workspace contracts and CPU/delayed-completion
mock tests in a separately owned ticket, reconciled with P10B/P10C. Then, if an
approved FP64-capable Linux NVIDIA machine is available, validate single-GPU
cuFFT and a resident local-potential pipeline. This yields the earliest useful
answer about precision, transfers and library integration. Do not start with
MPI, multi-GPU FFT decomposition, custom GPU eigensolvers or mixed precision.
Small Si workloads keep CPU as the default; no automatic atom-count dispatch
threshold is justified. CPU MKL gains from P10B are part of the future comparison,
not a reason to benchmark GPU against an obsolete CPU implementation.

Stage order separates correctness from adoption: a correct isolated FFT can
pass its scientific gate while failing its performance gate. Failure means keep
CPU and record the result. No stage changes tolerances, solver defaults or the
physics model to produce a speedup. Metal research does not gate CPU or CUDA work.

## 2. Shared numerical and measurement contract

All stages retain `scripts/check.sh` (default CPU, actual exit code), relevant
ignored release tests, existing fixed-point certification and unchanged QE
comparison scripts. GPU integration must also run default builds on a host with
no GPU SDK. External-library bit equality is not the cross-backend oracle.
Existing CPU overhead-only changes retain exact density hashes. Repeated GPU
runs compare bits within the supported fixed deterministic configuration and
across host thread budgets; unverified configurations are not production modes.

Before a hardware stage, fix a reproducible fixture manifest: geometry in Bohr,
cutoff in Ha, basis/grid/band/k sizes, pseudopotential SHA-256, precision,
reduction mode, driver/runtime/library/toolchain versions, git revision, thread
and stream counts. Store measured residuals, orthogonality, electron count,
eigenvalues and energy, not just PASS. Preserve existing QE provenance; no new
QE result is claimed without an actual authorized run.

Proposed additional FP64 adapter gates (dimensionless unless specified) are
initial acceptance targets, not measured GPU errors:

- FFT round trip and CPU transform comparison: relative L2 <=1e-12; Parseval
  relative error <=1e-12. Define relative L2 as norm(error)/max(norm(reference),
  1) for dimensionless fixtures; also report absolute maxima. Exact zero and
  signed single-mode tests check sign, scale and sphere support independently.
- H application: for normalized input, max-column norm(Hgpu C-Hcpu C) <=1e-10 Ha
  on explicitly scaled fixtures; Hermiticity defect
  abs(x†Hy-(Hx)†y)/max(1 Ha,abs(x†Hy),abs((Hx)†y)) <=1e-12.
- GEMM/Gram: relative Frobenius error <=1e-12 against CPU; orthonormal columns
  max abs(Q†Q-I) <=1e-10. Include dependent/nearly dependent vectors and rank
  outcomes, not just random well-conditioned matrices.
- Eigensolver: compare eigenvalues <=1e-10 Ha on bounded small fixtures; all
  fresh true residuals <= the requested solver tolerance (e.g. 1e-10 Ha),
  orthogonality <=1e-10. Degenerate subspaces compared by singular values or
  projector distance, not vector signs. Include dense CPU oracle checks for
  missing low roots and warm-start traps.
- Density: relative integrated electron-count error <=1e-12. For same cell,
  report scaled L2 `sqrt(Omega) * density_distance` and maximum density error in
  electron/Bohr³; proposed scaled L2 <=1e-8 for fully converged matching cases.
  Keep existing fixture-specific convergence tests if stricter.
- Full SCF: total-energy difference <=1e-12 Ha on existing benchmark cases,
  unchanged real-QE comparison tolerances, and every existing eigen-residual,
  density, charge, energy-ledger and fixed-point gate. Adapter unit gates alone
  do not imply this stronger full-run gate. Large/new fixtures need documented
  measured CPU noise floors and independently reviewed accuracy budgets before
  promotion; do not scale away an energy mismatch merely by atom count.

P10G supplies diversity: cubic/skewed/elongated/slab-like cells, 8–500 atom sizing
bins where practical, Gamma and multiple shifted k points, few/many bands,
small/large insulating gaps, synthetic projector structures including no
projectors, multiple channels and non-diagonal couplings, varied valence/cutoff.
Keep real Si as one regression, not the adoption benchmark. Current si8 real-UPF
nonconvergence must be resolved in its owning solver ticket before its SCF timing
can be an oracle. No metal/smearing case is invented under insulating occupations.

Timing protocol: release builds, same quiet host and physical input, fastest
validated CPU baseline with explicit thread budget, at least five repeats after
one warmup. Report median, min/max and variability for cold setup, warm kernels,
transfers/synchronization, peak host/device memory and full wall SCF including
certificate, audit and final export. Report transfer bytes by category. Device
events measure kernels only; a host monotonic clock brackets enqueue through
completion. Initial plan/JIT setup appears in cold end-to-end numbers. Synthetic
kernel timing cannot certify full SCF speedup. A proposed adoption gate is >=1.2x
median transfer-inclusive improvement with non-overlapping measured timing ranges
on at least two distinct supported families; otherwise retain an experimental
option. This is an engineering gate, not a predicted gain or changed default.

## 3. Stages

| Stage | Work and exit evidence | Hardware / oracle / timing |
| --- | --- | --- |
| G0: contracts and readiness | Additive neutral buffers, plans, capabilities, failure/completion state machine; CPU/mock implementation; reconcile merged P10C signature and P10B selection. Document budgets and alias rules; no new SDK required. | Apple/Linux CPU. Test delayed errors, stale potential/HV generations, OOM, unsupported precision, disjoint views, cancellation and explicit snapshot semantics. Default CPU result/hash gates unchanged. |
| G1: single-GPU FFT | Optional FP64 adapter; 3D and batched complex transforms, normalization/scatter/gather; queried persistent workspace. Test tiny, mixed-radix, prime and non-cubic grids, reuse and invalid dimensions. Compare full vs pruned semantic outputs. | Approved Linux NVIDIA GPU with installed compatible driver/cuFFT and FP64; CPU FFT plus small direct-DFT oracle, shared FFT gates. Measure upload+FFT+download and resident repeated FFT separately, including plan cost. No full-SCF claim. |
| G2: GPU Hamiltonian | Resident local V pipeline through P10C semantics, kinetic and nonlocal block projections, density accumulation; ordered reductions and potential-generation invalidation. Explicit CPU solver bridge may transfer correction blocks at this stage. | Same GPU, approved cuBLAS if used; CPU H and density references, H/Hermiticity/charge gates, multiple projector fixtures and shifted k. Transfer-inclusive H-block and complete SCF bridge timing. Retain CPU bridge if useful only as a test oracle; do not label it resident Davidson. |
| G3: GPU linear algebra | Resident GEMM/Gram, residual updates, diagonal preconditioning, orthonormalization and rotations; small EVD can stay CPU with explicit reduced-matrix staging. Bounded reusable output views replace width-indexed allocation patterns. | FP64 GPU plus cuBLAS; CPU faer/MKL and analytic/rank-deficient cases. GEMM/Gram/orthogonality gates. Time complete subspace iteration including small-matrix transfers; measure size-dependent EVD policy, never assume GPU wins for m=16–48. |
| G4: GPU eigensolver | Backend-neutral Davidson state with C/V/HV/corrections resident; true residual certificate and restart fallback preserved. Test optional device EVD through cuSOLVER only when reduced size warrants it. RMM-DIIS remains an independent switchable module after P10A, not a second mandatory GPU solver. | Same hardware with validated solver library; small dense CPU reference, clustered/degenerate spectra, failed recurrence certification and warm-start missing-root tests. Full solver gate; transfer-inclusive solve, not just GEMM timing. Host receives scalar summaries and small reduced outputs only. |
| G5: full GPU SCF | Resident density/audit contractions; initially host Hartree/XC/mixer with one density down and Veff up per step. Move these grid operations later if profiling warrants it. Host controller retains occupation/gap/ledger/certificate decisions. Explicit final legacy ScfResult snapshot. | Same GPU, P10G suite plus existing real-Si/QE evidence. Full SCF gates and adoption timing protocol; measured memory must fit without hidden paging. Report partial host residency honestly. Separate optional Metal precision study cannot borrow this stage's FP64 certification. |
| G6: multiple GPUs on one node | First distribute independent k points with ordered density/energy folds. Add band decomposition or distributed FFT only when one-device capacity/throughput proves insufficient. Budget duplicated projectors and reduced collectives. | At least two approved GPUs; actual PCIe/NVLink topology and peer support measured. Compare with one-GPU and CPU references at unchanged gates; deterministic order independent of partition count. Report peer bytes, imbalance and full wall timing. No NVLink assumption on PCIe-only hosts. |
| G7: MPI across nodes | Separate ownership/ADR and dependency approvals. Process/rank placement, GPU-aware collectives or staged copies, failure policy and distributed checkpoints; physical state contracts unchanged. | At least two approved nodes, fabric, MPI and compatible drivers. Same scientific gates across rank counts; full wall strong/weak scaling and transfer-inclusive collectives. No claim of MPI portability from a single-node test. |

No fixed delivery dates are assigned to unavailable hardware stages. Each future
implementation ticket must estimate duration once hardware, fixture scope and
SDK policy are known; this design does not authorize remote changes.

## 4. Decisions reserved for the user

No answer is needed to finish this document-only ticket. These are prerequisites
for future implementation authorization:

| Decision | Options | Recommendation and reason |
| --- | --- | --- |
| Hardware | Existing approved Linux NVIDIA FP64 machine; another vendor's available accelerator; defer hardware work | Use an existing NVIDIA FP64-capable machine with its already installed stack if available, with explicit identity/access approval. Otherwise stop at G0. Confirm actual memory, FP64 throughput and interconnect, not only GPU brand. No purchase/install requested here. |
| Apple GPU target | CPU-supported Apple only; optional Metal accuracy experiment; mandatory production Metal GPU | Keep CPU-supported Apple and allow a separately approved Metal precision experiment if useful. Native FP64 parity is not established, so mandatory production Metal is premature. |
| Proprietary SDK/build policy | User-supplied optional dynamically linked SDK; redistributed SDK/runtime bundle; avoid proprietary GPU libraries | Follow ADR 0019's user-supplied optional SDK policy for the first CUDA experiment, compile-time selected and no bundled runtime. If proprietary libraries are disallowed, reassess HIP/Vulkan hardware and library readiness before G1. |

A later reviewer must also approve pinned dependency licenses/MSRV, FFI ownership,
feature interactions and each experimental precision mode. None is approved by
this ADR proposal alone.

## 5. Risks and stop conditions

- Precision: Metal/native portable FP64 missing or cross-backend 1e-12 Ha gate
  fails -> keep CPU, investigate normalization/reductions and conditioning;
  do not relax tolerances.
- Capacity: 500-atom model subtotal is already 41.297 GiB per active k before
  library workspaces -> query peak, bound concurrency, reject OOM early.
- Small jobs: launch/plan/transfer cost or modern CPU performance erases gains
  -> retain CPU; publish measured crossover surface rather than atom cutoff.
- Determinism: concurrent reductions/library algorithm changes -> fixed execution
  configuration, independent repeated tests; label unsupported combinations.
- API leakage: a host-slice compatibility shim downloads every vector -> stop
  residency claim until solver and audit use neutral resident resources.
- Maintenance/build: binding API churn, unsupported hardware/SDK, incomplete
  solver coverage -> pin and isolate adapter; never force SDK into default CI.
- Scope: multi-GPU/MPI and mixed precision increase both numerical and operational
  risk -> separate gates/authorization after a validated single-device path.

## 6. P10E verification record

On the inspected Apple arm64 checkout, `scripts/check.sh` exited **0** with
**PASS build, PASS test, PASS fmt, PASS clippy**. These are existing CPU tests,
not GPU correctness or timing evidence. Sizing equations were evaluated using
Python standard-library arithmetic; they are models with explicit assumptions.
Local probes found Metal-capable hardware but no established usable FP64 stack;
prototype skipped (details in ecosystem section 7). No remote connection,
installation, new QE run, dependency or code modification occurred.

Only the three requested documents and ADR proposal are owned. README numerical
results/figures are unchanged because there is no verified implementation or
new scientific evidence. Session/dashboard artifact generation and upload would
write outside this ticket's ownership; left to the orchestrator rather than
silently modifying shared state. Changes remain uncommitted as required.
