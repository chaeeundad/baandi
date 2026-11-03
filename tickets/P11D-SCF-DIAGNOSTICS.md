TASK_ID: P11D-SCF-DIAGNOSTICS
GOAL: Deterministic SCF diagnostics over the event stream: stagnation, oscillation, divergence, non-finite values, charge sloshing hint; each finding has evidence (numbers) and an optional hypothesis (labelled) and a list of candidate policy changes. Directive item B4.
OWNED_FILES: crates/dft-engine/src/diagnostics.rs (new), lib.rs, tests/diagnostics.rs
DEPENDENCIES: P11A merged (uses events); P11B for candidate policy types
IMPLEMENTATION_STEPS: pure functions of the iteration history; exact numeric definitions documented (window sizes, ratios) and justified on a set of different systems (generality rule), no tuning on Si only; synthetic residual histories as unit tests plus real runs (including the 8-atom stagnation recorded in P7PC, large-beta divergence).
ACCEPTANCE_CRITERIA: each rule triggers on its synthetic pattern and not on healthy runs of >= 4 different systems (zero false positives on the generality set); deterministic.
GIT: branch task/P11D-SCF-DIAGNOSTICS, uncommitted
RETURN_FORMAT: docs/schemas/report.schema.json
