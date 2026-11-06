TASK_ID: P25A-BENCH-LARGE
GOAL: Extend the benchmark set so the 8 to 16 core comparison against QE 7.6 is measurable (QE times on si8 and k444s are below its 0.3 s timer resolution). Add larger cases and a second-phase workload group (forces, stress, relaxation step). Inputs only; no QE run here (QE runs go through dft5090 on dev08).
OWNED_FILES: benchmarks/cases/*.toml (new files only), benchmarks/README.md, docs/performance/benchmark-suite.md, scripts/bench_suite.py (only to support a new `phase` field, backward compatible), scripts/tests/ (tests for the new field)
DEPENDENCIES: none
IMPLEMENTATION_STEPS: (1) add si32-gamma (2x2x2 conventional cells, 64 atoms is too large: use 32 atoms, 1 core measured time target 20 to 120 s on dev08), si16-k222, mgo-40 variants with at least 0.3 s QE-expected time per SCF iteration; follow the existing case schema and tier rules (long tier if over five minutes on the Mac). (2) add optional field `phase` (`scf` default, `forces`, `stress`, `relax1`) so the suite can time force, stress and one relaxation step separately; reports list the phase. (3) document in benchmark-suite.md the rule that QE timings under 0.3 s are invalid and which new cases are above that line (state estimates as estimates until dev08 measures them). (4) no potential files, no numbers recorded as QE reference.
TEST_COMMANDS: scripts/check.sh; python3 -I -m unittest discover scripts/tests
ACCEPTANCE_CRITERIA: new cases parse with the existing validator; existing cases and reports unchanged; phase field tested; docs state which cases await QE reference runs.
NON_GOALS: running QE, changing solver code, new potentials.
GIT: branch task/P25A-BENCH-LARGE, uncommitted
RETURN_FORMAT: docs/schemas/report.schema.json
