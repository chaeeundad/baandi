TASK_ID: P11J-RUNTIME-BENCH
GOAL: Matched benchmark and regression report of the three modes: total wall time, per-iteration overhead, event generation/delivery cost, memory overhead, checkpoint overhead, control-request latency, reproducibility (bit identity), queue saturation behavior; short SCF workloads included. Compare against the pre-refactor path from docs/performance/baseline*.json.
OWNED_FILES: scripts/bench_runtime.py, docs/validation/runtime-equivalence.md (results), docs/performance/optimization-log.md (append), tests/reference/runtime/
DEPENDENCIES: P11A, P11B merged
ACCEPTANCE_CRITERIA: repeated runs with variability; same inputs/hardware per mode; targets (0.5%/1%/1%) reported as met or not met with numbers; Claude re-measures on quiet dev08 via the dft5090 session.
GIT: branch task/P11J-RUNTIME-BENCH, uncommitted
RETURN_FORMAT: docs/schemas/report.schema.json
