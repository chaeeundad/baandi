# ADR 0021: Split scientific workflows from numerical coupling

- Date: 2026-11-03
- Status: proposed; **Designed**, no implementation authorization or completion claim.
- Number check: local checkout ends at 0019; sibling P11R worktree already contains `0020-agent-native-runtime.md`. Reserve 0021 here to avoid that known collision; recheck numbering when integrating concurrent branches.
- Context: P11I requires concrete interfaces extending [ADR 0017](0017-coupling-commitments.md) without changing numerical code. [ADR 0016](0016-esle-layering.md) governs electronic-state identity. P11R's draft runtime/state/control/integration/solver documents were read from its sibling worktree, are not merged here, and require integration review before implementation.

## Decisions

1. **Designed:** separate a high-level [workflow API](../architecture/workflow-api.md) for structures, qualified results, capabilities, identity/provenance, datasets, Python/ASE/MLIP and agents from a [low-level API](../architecture/coupling-api-lowlevel.md) for typed fields and apply-only numerical operators. They share model/environment identity, status, availability and artifact types; neither defines an alternative SCF implementation.
2. **Designed:** reuse P11R's single runtime and checkpoint envelope. ESLE remains its payload. The later P11 common three-mode, bounded-event and stale-control addendum takes precedence over the reviewed P11R draft. P11G owns Python/ASE; P11F owns the separate single MCP server consumed by P9A. No new core dependency or crate is introduced here. Future binding/server crates require their ticketed ownership and dependency reviews.
3. **Designed:** distinguish phi (Ha/e), electron energy U (Ha), number density, signed charge and electric field; U=-phi. Specify cell columns, origin, FFT/basis normalization, boundary/reference, energy ownership, buffer lifetimes and version negotiation. Existing Hamiltonian setters or Hartree solves are not external coupling capabilities.
4. **Designed:** arrays use leased in-process memory or binary exchange with validated descriptors, immutable identities and bounded resource use. Durable arrays reuse P11R ArrayRef, not another checkpoint format. No JSON arrays or LLM/HTTP/Python callbacks in numerical loops. Operator leases freeze the potential and basis, never expose a mutable matrix or imply nonlinear response.
5. **Designed:** the [active-learning contract](../architecture/active-learning-contract.md) carries every stage from simulation to controlled deployment. Optional uncertainties require a calibrated flag and evidence when true. State labels remain distinct at equal geometry; unresolved states are never merged by energy. Model changes, dataset partitions, alignment and deployment require explicit provenance; no silent model mixing or missing-force zeros.
6. **Designed:** inner certification and outer coupled convergence are separate. Complete ledger/functional and validated derivatives are prerequisites for coupled energy/force/stress claims. The low-level document lists energy, sign, grid, boundary, force, stress, coupling, lifetime and transport tests that must precede QM/MM claims. No test result is asserted by specifying a test.
7. **Implemented** labels in the gap table apply only to inspected existing internal building blocks. All proposed public APIs remain **Designed**; unsupported physics and concrete implementations named below are **Deferred**. Runtime capability flags describe only tested implementation combinations, never this design roadmap.

## Deferred

Concrete ABI, shared-memory/IPC implementation, Python/MCP packages, array import, external operator leases, workflow executors, dataset/training/deployment services, exact restart guarantees and canonical digest fixtures. Forces/stress, QM/MM partitions/link atoms, charged/background/open/slab models, uniform electric fields, reconstructed electrostatics, density/geometry/time-dependent and nonlocal embedding remain separately gated. A callback sketch does not authorize any of these physical models. No target performance, calibration or scientific interoperability is established.

## Alternatives rejected and consequences

Reject one oversized calculator trait mixing control, physics and transport; arbitrary executable callbacks in SCF; JSON field arrays; exporting V_eff as electrostatics; zero-filled unsupported derivatives; deduplicating labels by geometry/energy; one mutable MLIP alias standing in for provenance. These hide physical or ownership assumptions and make verification impossible.

Keep all existing APIs/defaults unchanged. Implementers must resolve P11R merge differences, canonical identity encoding and concrete transport/dependency choices under later tickets, then run the specified numerical and resource gates. The [efficiency note](../research/external-coupling-efficiency.md) supplies a measurement plan, not an optimization result. This ADR extends 0017's semantic commitments while retaining its implementation deferrals; accepted earlier ADR text is not edited.
