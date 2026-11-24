TASK_ID: P3-QECLOCK
GOAL: |
  A stage-by-stage speed comparison with QE 7.6 needs both codes timed at the same boundaries. QE always prints its clock table (routine, CPU and wall seconds, calls) from rank 0. Make the reference driver print and record a table with the same boundaries, statistic and call counts, so one profiled run of each code can be compared line by line (diagnosis only; headline timings stay profiling-off).
READ FIRST: docs/research/qe-study/qe-stage-by-stage.md, qe-scf-cycle.md, qe-hpsi-and-fft.md, davidson-iteration-structure-gamma.md (own-words notes; never read QE source), docs/performance/comparison-rules.md ("Diagnostics level").
TASK: |
  1. Define the common boundary map in docs/design/reference-scf.md section "QE-matched clock table": only lines whose boundaries can be matched exactly in meaning: setup (init_run), SCF total (electrons), eigensolver (c_bands), H application (h_psi) with local potential incl. FFTs (vloc_psi), FFT (fftw), nonlocal projections (calbec) and back-application (add_vuspsi), reduced eigensolve (rdiaghg), Davidson remainder (regterg minus h_psi minus rdiaghg), density (sum_band), potential (v_of_rho), mixing (mix_rho), plus BAANDI-only lines (verification, output). For each line state what is inside on both sides; where a QE line has no exact BAANDI boundary, say so instead of approximating.
  2. Collect these clocks in the reference driver when BAANDI_DAVIDSON_STAGES=on, with call counts, on rank 0 (and keep max/min as now). Print a QE-style table at the end of the run to stdout and write it to results.json (reference_scf.qe_clock: name, wall_s, calls, rank). Zero cost when off.
  3. Mac checks only (user rule 2026-11-24/25): cargo check --workspace --tests; cargo test -p dft-engine --test refscf; cargo test -p dft-solvers; cargo test -p dft-physics; fmt; clippy -D warnings. Default results bit-identical with the switch off and on.
  4. Write docs/performance/lanes/RSCF/qeclock-dev08-request.md: ONE profiled run each of BAANDI (stages on) and QE (its matched Davidson input, same cell, cutoff, grid, UPF, bands, start, threshold, ranks, binding) on si64 and si128 MPI8, and a parser that puts both tables side by side (per-call time and calls per line).
OWNED_FILES:
  - crates/dft-engine/src/refscf/record.rs   # clock-table fields only
  - crates/dft-engine/src/refscf/driver.rs   # clock collection only; I25-VERIFY-ENERGY edits verification code in this file in parallel, keep your edits local to timing
  - crates/dft-solvers/src/ref_davidson.rs   # counters only
  - crates/dft-solvers/src/ref_davidson/accounting.rs
  - crates/dft-physics/src/distributed_real_gamma.rs   # counters only
  - crates/dft-cli/src/run.rs   # qe_clock table output only
  - docs/design/reference-scf.md   # "QE-matched clock table" section only
  - docs/performance/lanes/RSCF/qeclock-dev08-request.md
  - scripts/qe_clock_compare.py   # new parser
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-solvers
  - cargo test -p dft-physics
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-solvers -p dft-physics -p dft-cli --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: boundary map documented with exact/inexact marked; QE-style table printed and serialized on rank 0 with calls; zero cost and bit-identical when off; parser and dev08 request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
