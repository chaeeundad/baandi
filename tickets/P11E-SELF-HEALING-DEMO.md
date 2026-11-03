TASK_ID: P11E-SELF-HEALING-DEMO
GOAL: End-to-end Self-Healing SCF demo: real calculation -> events -> diagnostics -> controller (deterministic rules, then a scripted/recorded agent decision source) -> validated policy change at a safe point -> continue -> verified result -> provenance report; three arms compared (baseline, deterministic recovery, agent-supervised).
OWNED_FILES: crates/dft-cli/src/heal.rs (new), crates/dft-cli/tests/self_healing.rs, scripts/demo_self_healing.sh, docs/demos/self-healing-scf.md, docs/validation/intervention-safety.md, docs/validation/runtime-equivalence.md (fill results), tests/reference/runtime/ (generated evidence)
DEPENDENCIES: P11A, P11B, P11D merged; P11C optional
ACCEPTANCE_CRITERIA: reproducible command; baseline fails or is slower, recovery converges to the same total energy within 1e-10 Ha of a fully converged reference; every intervention recorded; rejected unsafe change demonstrated; no network or LLM needed to reproduce (agent arm replays a recorded decision file); results summarized in README only through the generated results block (extend scripts/make_readme_figures.py).
GIT: branch task/P11E-SELF-HEALING-DEMO, uncommitted
RETURN_FORMAT: docs/schemas/report.schema.json
