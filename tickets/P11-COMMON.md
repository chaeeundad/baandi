# Common rules for the P11 tickets: Agent-Native Scientific Computing Runtime (read first)
Source: the user's master directive "BAANDI - Agent-Native Scientific Computing Runtime" (summary and review: docs/research/agent-native-directive-review.md). Read AGENTS.md, tickets/P10-COMMON.md (generality, correctness, workflow; all of it applies), docs/architecture/{electronic-state,scientific-coupling}.md, docs/decisions/0016 and 0017, docs/architecture/user-interface.md when present.

Layers: A numerical core (no LLM/MCP/HTTP/Python dependency, ever); B typed scientific state (units, metadata, schema versions); C runtime control (stepwise execution, state machine, pause/resume, checkpoint, events, safe-point policy changes); D integration (Python/ASE, MLIP, external solvers); E agent tools (MCP). Dependencies point downward only: E -> D -> C -> B -> A.

Hard rules for every P11 ticket:
1. Equivalence first: the stepwise path and `ScfProblem::run`/`run_outcome` must produce bit-identical energy, density, history and certificate for the same input and thread count (check against docs/performance/baseline*.json). The conventional path stays; it becomes a thin driver over the stepwise one only after bit-identity is proven.
2. No LLM, network, MCP, serde_json or file I/O inside numerical hot loops or in dft-engine numerics. Events/snapshots are built from already computed scalars and small vectors; no large arrays in events (density/wavefunctions only by explicit request, in a binary or array format).
3. A policy change only takes effect at a safe point (between SCF iterations), is validated against a schema and a physical-invariant list before being applied, is transactional (rejected = state unchanged) and is recorded (who, when, why, before/after). Restricted changes (XC, pseudopotential, species, electron count, structure, k mesh, cutoff) are rejected while a calculation is running; they need a new calculation.
4. Events and checkpoints carry calc id, monotonically increasing sequence number, wall-clock timestamp, schema version. Determinism: a run with a recorded policy-change log replays to the same result for the same thread count.
5. Diagnostics distinguish evidence (measured numbers) from hypothesis (a suggested cause); never claim a cause as fact.
6. Convergence is only ever declared by the existing fixed-point certificate; nothing in the runtime layer may set `converged`.
7. Real UPF files are never sent anywhere; demos use synthetic fixtures or the internal Si validation file read-only on this machine.
8. Public workflow: branch task/<TASK>, Claude opens PR; no QE source; all docs in English.
9. Sequencing constraint: crates/dft-engine/src/scf.rs is being edited by P10A, P10D, P7PC and P10C. P11A starts only after those merge; design/doc tickets (P11R, P11H) start now.
