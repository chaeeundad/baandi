# External coupling efficiency: design and measurement plan

Status: **Designed**, P11I, 2026-11-03. No benchmark results or implementation claims. This is a repository-derived transport and scheduling proposal, not a survey of external package implementations. Normative contracts: [workflow API](../architecture/workflow-api.md), [numerical coupling API](../architecture/coupling-api-lowlevel.md), [active learning](../architecture/active-learning-contract.md).

## Separate call frequency and data volume

| Path | Payload/cadence | Designed transport | Current status |
| --- | --- | --- | --- |
| Scientific workflow | Requests, status, model/provenance, artifact IDs per calculation/batch | Typed library/Python calls; MCP only for agent control | Designed; P11G/P11F pending |
| Runtime control | Small scalar snapshots and boundary policy proposals | Bounded nonblocking event delivery outside numerics | Designed; P11R/P11 common addendum |
| Numerical operator | Complex vector/block per external eigensolver apply | In-process borrowed buffers, preallocated outputs and scratch | Internal LinearOperator Implemented; external lease Designed |
| Coupled outer iteration | Density, source/field arrays per completed inner solve | Immutable shared-memory leases locally, otherwise bounded binary IPC | Designed contract; transport implementation Deferred |
| Archive/visualization | Explicit on-demand export | Binary payload/manifest plus standard cube/XSF when supported | Existing density.bin Implemented; proposed richer exchange Deferred |

Remote/HTTP/LLM control is never an operator or SCF inner-loop callback. Large arrays never ride in events, JSON, MCP tool text or repeated status snapshots. An external environment solve is performed outside the completed electronic solve, with a new immutable request. Shared memory is a possible local IPC optimization, not proof that a remote solver has zero-copy access.

## Cost model and ownership

For N grid nodes, a real scalar f64 payload is **8N bytes**, a Cartesian real vector is **24N bytes**, and a complex block of p plane waves and b columns is **16pb bytes** (derived storage formulas, not measurements). One scalar field and one density transfer move at least **16N bytes** per outer exchange, excluding descriptors and any retained copies. Peak memory additionally includes simultaneous producer/consumer buffers, solver state, scratch, queued frames and retained snapshots. Report these separately; do not label payload size as peak memory.

Prefer a borrow only when dtype, shape, contiguity, alignment, native endianness, coordinate frame and lifetime already match. Otherwise make the conversion explicit and count its allocation and bytes. Never obtain zero-copy by lending mutable live SCF state to another process. Freeze the owner, synchronize completion, then publish a read-only bounded lease; release before reuse. Hash once at the immutable boundary, not per Hamiltonian application. A descriptor is not evidence of immutability: enforce process/shared-memory access and crash invalidation.

Batch apply can reduce call overhead for external eigensolvers, but block size and scratch must be measured across small/large band counts. No default is selected here. Session mutation cannot occur while an operator lease is live. Whole-session copying to avoid lifetime validation is not the default; an explicit retained snapshot is a budgeted alternative.

Use bounded bytes and lease counts in addition to bounded message counts. Routine events may coalesce/drop with a dropped-count record as specified by P11 common; critical certificate, state, policy and error records must be retained without blocking numerics. Implementation must reconcile queue capacity and durable critical storage outside the fast path before claiming this guarantee. Slow consumers cannot accumulate unbounded density snapshots. Backpressure on bulk exchange is allowed at an outer boundary, never hidden inside FFT/apply. Timeout/crash produces a typed failed exchange and preserves the last immutable completed result.

## Benchmark plan before implementation acceptance

1. Freeze hardware, OS, compiler, source/dirty patch/lockfile, numerical backend, thread topology, seed, model/grid/k/band identity and workload. Reuse P10G's diverse benchmark plan when available: cubic, skew, elongated and slab-like periodic cells, atom/band counts, Gamma/multiple k, gaps and projector structures. Synthetic fixtures are valid for transport/lifetime behavior; silicon-only timings do not set defaults.
2. Measure the existing conventional path, an adapter with no observation/intervention, read-only observation and bounded agent control with no intervention. Compare same-thread numerical energy, density, history and certificate bit-for-bit for overhead-only work. Use `scripts/bench_scf.py` and applicable stored baseline per repository performance rules; do not run or claim these implementation benchmarks for this design ticket.
3. Measure in-process borrow, explicit copy, local shared memory and binary IPC separately for scalar density, vector field and orbital blocks. Include validation/hash/conversion cost, startup, actual bytes copied, allocations, peak resident memory, lease count, throughput and end-to-end outer-iteration time. Repeated release runs report count, median, range and tail latency with queue saturation and slow-consumer/crash cases.
4. Measure short SCF cases where fixed overhead matters as well as larger arrays. Separate setup from steady state, but include full calculation time. Report callback/operator-call cost and total application counts so algorithm changes are not misreported as transport improvements.
5. Test disabled observers, saturated bounded queues, expired control proposals and crashed array producers; verify unchanged numerics and explicit missing-result status. Distinguish safe-point control latency from a full inner/outer solver latency.
6. Record experiments, rejected alternatives and accepted changes in `docs/performance/optimization-log.md` under a future implementation ticket. Runtime targets Classic <0.5%, Observed <1%, AgentControlled without intervention <1% are targets from P11 common, **not measurements or guarantees**. No new coupling throughput or speedup target is asserted.

## Compatibility and deferred decisions

Use existing P11R state/driver ownership, not a new coupling runtime, scheduler or serialization path. The P11R draft reviewed here predates the later three-mode/stale-request addendum; [workflow API](../architecture/workflow-api.md) records precedence. P11R's external-solver transport deferral is refined to semantic buffer rules only; concrete OS shared-memory mechanism, IPC framing implementation, ABI, package dependencies, device memory and remote bulk transport remain **Deferred**. Provider-specific performance claims require direct experiments and versioned evidence.

No coupling algorithm, QM/MM implementation, validated electrostatic reconstruction, force/stress calculation, restart guarantee or efficiency improvement follows from these documents. The [mandatory coupling gates](../architecture/coupling-api-lowlevel.md#mandatory-gate-before-any-qmmm-or-conservative-coupling-claim) must pass before scientific capability claims; transport benchmarks cannot substitute for them.
