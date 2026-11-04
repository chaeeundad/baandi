TASK_ID: P11L-CASE-STUDIES
GOAL: Detailed, reproducible case studies for every external control and coupling path once it works: MCP agent control, Python/ASE workflow, active-learning loop with an MLIP, QM/MM-style or external-solver coupling (low-level API), Self-Healing SCF, agent-assisted preflight. Each is a document a new user can follow and a reviewer can verify.
OWNED_FILES: docs/case-studies/*.md (new), docs/case-studies/data/ (small inputs, logs, provenance JSON, figures; no UPF files), docs/INDEX.md, scripts/ only for a case-study runner script
DEPENDENCIES: starts per case after its feature is merged: P11E (self-healing), P11F (MCP), P11G (Python/ASE), P11K (preflight), P11I implementation of the coupling API when it exists
STRUCTURE of each case study: 1 problem and why it matters; 2 setup (system, parameters, versions, pseudopotential SHA-256, exact commands); 3 baseline run (Classic) and the same run through the external path; 4 what the agent/external tool did, step by step, with the event log, control requests, accepted/rejected decisions and audit excerpt; 5 results with real numbers and matched measurements (energy difference, iterations, wall time, overhead versus Classic); 6 failure cases and stale or rejected requests shown honestly; 7 limits (what was not tested, system size, single machine); 8 how to reproduce in one command.
RULES: numbers only from runs executed for the document (command, build, date recorded); no claims beyond measured systems (generality rule: use at least two different systems per case where feasible); the agent is shown as optional and supervised; Korean summary for the user, English for the documents.
TESTS: a script re-runs each case on small systems and checks the recorded numbers within stated tolerances.
ACCEPTANCE_CRITERIA: every external path has one case study with a verified one-command reproduction; INDEX links them; README links the list.
NON_GOALS: new features, marketing claims, performance claims not measured on the recorded machine.
GIT: branch task/P11L-CASE-STUDIES, uncommitted
RETURN_FORMAT: docs/schemas/report.schema.json
