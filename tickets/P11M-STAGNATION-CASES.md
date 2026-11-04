TASK_ID: P11M-STAGNATION-CASES
GOAL: Find and document real SCF failure modes (stagnation, oscillation, slow charge-sloshing-like drift, divergence) on several different systems so that the Self-Healing SCF demo and the case studies rest on genuine cases, not a budget cut. Output: a reproducible set of input fixtures plus a measured table.
OWNED_FILES: tests/reference/stagnation/ (new: inputs, synthetic UPF generator calls, results JSON, scripts), scripts/find_stagnation.py (new), docs/validation/scf-stagnation-cases.md (new), docs/INDEX.md
DEPENDENCIES: P11D (diagnostics), P11E (demo), P10D (mixers), P7PC merged
METHOD: sweep cell shape (cubic, skewed, elongated, slab-like with vacuum), atom count, k mesh and Gamma, band count (small/large), gap size (small to large, including near-metallic synthetic models), mixing beta (0.1 to 1.0), mixing kind (Linear, Pulay, Kerker-preconditioned), eigensolver (Davidson, Hybrid RMM). Use synthetic UPF fixtures from scripts/gen_synthetic_upf.py (never commit real UPF). Record per case: converged or not, iterations, final residual, diagnostics classification from the P11D module with evidence, and whether a P11B-valid intervention (beta reduction, history reset, mixer switch) recovers it with the same final energy within 1e-10 Ha.
RULES: report honestly if no case fails with defaults (that is a result); do not tune solver defaults in this ticket; no claim of generality beyond the swept systems; keep runs small enough for the local machine (each case under about 2 minutes in release).
ACCEPTANCE_CRITERIA: at least three distinct failing or badly converging cases found (or a documented statement of what was swept without finding them), each reproducible with one command, each with before/after intervention numbers; table of everything swept including non-failing cases.
NON_GOALS: changing solvers, new physics, performance claims.
GIT: branch task/P11M-STAGNATION-CASES, uncommitted
RETURN_FORMAT: docs/schemas/report.schema.json
