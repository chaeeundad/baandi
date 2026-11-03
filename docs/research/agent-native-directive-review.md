# Review of the Agent-Native Scientific Computing Runtime directive (2026-11-03)
Source: user's master directive (five layers A numerical core, B typed state, C runtime control, D integration, E agent tools; phases 1 observable, 2 controllable, 3 agent tools, 4 interoperability, 5 provenance/workflow; first demo: Self-Healing SCF).

## Verdict
Adopted as the product direction. The layering (core free of LLM/MCP/HTTP) matches ADR 0017. The directive is sound; the points below adjust it to the live repository.

## Repository grounding
- `ScfProblem::run_impl` (crates/dft-engine/src/scf.rs) is a monolithic loop. It already has `step_kpoints` (one eigen solve + density), `IterationRecord` history, a fixed-point certificate, divergence detection and `run_outcome`. A1 therefore means extracting loop state into a session, not writing SCF anew; bit-identity with `run` is the gate.
- Mixers are already objects (`DensityMixer`); mixing/tolerance policy changes are feasible at iteration boundaries. Mixer history semantics on change (reset vs keep) must be specified.
- ESLE phase B already plans checkpoint/state; A4 must share one format (ESLE state = payload). Serialization lives in dft-cli (ADR 0016/0017), so dft-engine stays serde-free.
- P9A (user interface design) planned an agent/MCP layer; there must be exactly one MCP server (P11F); P9A consumes it.

## Conflicts and sequencing
1. scf.rs is being edited now by P10A (RMM-DIIS), P10D (mixing), P7PC (buffer bands), P10C. P11A (the refactor) must wait for their merge; otherwise a four-way conflict in the hot file. Design tickets run now.
2. Self-healing needs real recovery levers: P10D (Kerker/preconditioned mixing) and P7PC (buffer bands) give them, so the demo gains from finishing those first.
3. Determinism: policy changes at safe points are deterministic only if the change log is replayed at the same iteration and thread count; this is a stated acceptance test.
4. The directive's "agent-supervised" arm must be reproducible without an LLM: the demo replays a recorded decision file; live LLM runs are optional extras.
5. Intervention must never touch convergence: only the certificate sets `converged`; a policy may not relax the certificate threshold silently (recorded as `unsafe_relaxed` if ever allowed).
6. Python/ASE: forces and stress do not exist in BAANDI yet; the calculator must raise NotImplemented, never return zeros.
7. Dependencies (MCP SDK, PyO3) need license/maintenance/Apple Silicon review and ADRs; new crates sit outside the numerical core.

## Ticket mapping
P11R design docs + ADR (now); P11A A1-A3 (after P10A/D/C, P7PC); P11B B1-B3; P11C A4; P11D B4; P11E demo; P11F C1-C2; P11G D1; D2/D3/E1/E2 are design sections in P11R only (implementation later, no stubs).
Dependency order: P11R -> P11A -> {P11B, P11C} -> P11D -> P11E; P11F after B/D; P11G after A.
