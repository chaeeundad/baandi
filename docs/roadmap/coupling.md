# Coupling roadmap (gated; design only now)

Status: plan from [../research/coupling-directive.md](../research/coupling-directive.md), reconciled with [../research/coupling-codex-review.md](../research/coupling-codex-review.md). Gates advance on evidence for their physics, never because a transport adapter works. Every gate needs the user's approval before implementation starts.
| Gate | Content | Entry | Exit evidence |
|---|---|---|---|
| 0 Design (now) | semantics, contracts, test plan, ADR 0017 | directive received | docs reviewed; no code; k-point work untouched |
| 1 Baseline | baseline SCF with validated k-points (done for Si on three grids) and the ESLE state infrastructure (phase B) as needed | R6 recorded (done) | all existing gates stay green; real-Si Gamma density bit-identical |
| 2 Fixed-field validation | first scope: neutral fixed-N, exact-grid, real periodic electron-energy field; ledger with direct and band forms; reference policy; metadata rejection | gate 1 and user approval of scope and policies | constant-shift, grid and origin, sign, weak-field response, double-counting, serialization and rejection tests with measured tolerances; zero-field path unchanged |
| 3 Observables and derivatives | typed observables, forces and then stress as separate physics tickets (ordinary terms first, then coupling terms); energy-only datasets before force data | gate 2 | finite-difference validation of a complete stationary functional; no MLIP force data before this |
| 4 External drivers | request, result and lifecycle types in an adapter; one-way fixed-environment driver | gate 3 | idempotence, cancellation, partial-result and separate inner and outer status tests |
| 5 Coupled-state exploration | two-way drivers, density-dependent embedding, ESLE under an environment | gate 4 and ESLE phase C | documented coupled functional with derivative consistency; ranking only by complete free energy |
## Decisions needed from the user before gate 2
1. First scope: electron-only field, or also a supplied electrostatic ion contribution.
2. Potential reference (default recommended: keep the supplied mean and record it), source regularization, charged-cell policy (initial rejection recommended).
3. Energy ownership and output names (electronic partial, interaction, source and environment self energy, coupled total or free energy); whether potential-only models are rejected or diagnostic-only.
4. Joint ESLE fingerprint extension and the warm-transfer versus replay promise; branch labeling governance.
5. Observable order (fixed-field energy and density first; forces and stress as separate tickets).
6. Optional adapter dependencies only for a concrete consumer, each with license, maintenance and Apple Silicon checks.
7. Calibration fixtures and response acceptance envelopes after the current k milestone.
