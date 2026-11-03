TASK_ID: P11I-COUPLING-API-SPEC
GOAL: Design-only interface specifications for external coupling, extending ADR 0017 and docs/architecture/scientific-coupling.md; no implementation claims.
OWNED_FILES: docs/architecture/workflow-api.md, docs/architecture/coupling-api-lowlevel.md, docs/architecture/active-learning-contract.md, docs/research/external-coupling-efficiency.md, docs/decisions/<next free>-coupling-api-split.md (check next number)
IMPLEMENTATION_STEPS:
 1. High-level scientific workflow API (ASE, Python, MLIP, active learning, research agents): structures, energies, forces and stress when supported (explicit capability flags, never zeros), calculation status, model identity, provenance, validation metadata, dataset references; Python binding plan (P11G) and the separate MCP adapter (P11F).
 2. Low-level coupling API (QM/MM, DFT-CES, embedding, external numerical solvers): typed external potentials (electrostatic potential vs electron potential energy, sign), density/electrostatic field exchange, Hamiltonian operator access (apply-only, ownership and lifetime), boundary conditions, energy/force coupling terms, callbacks only where mathematically valid, memory-efficient exchange (shared memory/zero-copy in process, binary IPC otherwise; no JSON for arrays, no LLM/HTTP in tight loops). Specify units, basis, coordinate frame, FFT normalization, ownership, lifetime, versioning.
 3. Active learning contract: MLIP simulation -> candidate screening -> DFT labeling -> dataset validation -> model training -> model evaluation -> controlled deployment. No hard-coded MLIP (SevenNet etc. only as examples); uncertainty estimators optional with a `calibrated` flag and calibration reference; state labels preserved (ESLE) so distinct electronic states at the same geometry are not merged; no silent model mixing; dataset provenance (code version, input hash, pseudopotential identity hash, convergence certificate).
 4. Test plan required before any QM/MM claim: energy/force/boundary/coupling consistency tests listed explicitly; statement of what is NOT claimed.
 5. Mark each interface Implemented / Designed / Deferred.
ACCEPTANCE_CRITERIA: documents concrete (types, fields, units), consistent with P11R docs and ADR 0016/0017, a gap table versus current code.
NON_GOALS: code.
GIT: branch task/P11I-COUPLING-API-SPEC, uncommitted
RETURN_FORMAT: docs/schemas/report.schema.json
